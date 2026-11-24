"""
📂 analyse.py – wertet automatisch alle FASTA-Dateien in einem Ordner aus.

Für JEDE Sequenz:      Länge, Basenzusammensetzung, GC-Gehalt, ORFs (mögliche Gene) + Proteine
Für JEDES Paar:        Alignment (Needleman-Wunsch), Identität, Lücken; bei gleicher Länge Mutationsliste
Für ALLE zusammen:     Distanzmatrix und Stammbaum (UPGMA + Neighbor-Joining), ab 3 Sequenzen

Große Sequenzen: Statistik und ORF-Suche gehen auch für ganze Bakteriengenome (Millionen Basen).
Alignments brauchen Zeit ~ Länge1 × Länge2 und werden deshalb nur bis MAX_ALIGN_LAENGE gemacht.
Für den Stammbaum wird dann automatisch eine alignment-freie k-mer-Distanz benutzt.

Alle Ergebnisse landen als Dateien im Ausgabe-Ordner.

Benutzung (in Colab oder lokal):
    from analyse import analysiere_ordner
    analysiere_ordner("daten")
"""

import os
import time

from werkzeuge import (lese_ordner, gc_gehalt, zusammensetzung, finde_orfs, needleman_wunsch,
                       alignment_statistik, alignment_als_text, mutationen, reverse_complement,
                       distanz_aus_alignment, kmer_menge, kmer_distanz, upgma, neighbor_joining,
                       newick, baum_als_text, zeichne_baum)

STANDARD = {
    "ORF_MIN_AMINOSAEUREN": 100,          # kürzere ORFs ignorieren
    "STARTCODONS": ("ATG",),              # für Bakterien realistischer: ("ATG", "GTG", "TTG")
    "MAX_ALIGN_LAENGE": 5000,             # längere Sequenzen werden nicht aligniert
    "MAX_ALIGN_ZELLEN_GESAMT": 150_000_000,  # Summe aller Länge1×Länge2; ca. 1–2 Minuten in Colab
    "MATCH": 1, "MISMATCH": -1, "GAP": -2,
    "KMER_K": None,                       # None = automatisch nach Sequenzlänge
    "AUSSENGRUPPE": None,                 # Liste von Namen für den NJ-Baum, z. B. ["Cypripedium_calceolus"]
    "ZEIGE_GRAFIK": True,
}


def _k_automatisch(min_laenge):
    if min_laenge >= 100_000:
        return 21          # Standard von Mash für Genome
    if min_laenge >= 5_000:
        return 15
    return 11


def _ausgabe(bericht, text=""):
    print(text)
    bericht.append(text)


def analysiere_ordner(ordner, ausgabe_ordner="ergebnisse", **einstellungen):
    E = dict(STANDARD)
    for schluessel, wert in einstellungen.items():
        if schluessel not in E:
            raise ValueError(f"Unbekannte Einstellung: {schluessel}. Möglich: {', '.join(E)}")
        E[schluessel] = wert

    os.makedirs(ausgabe_ordner, exist_ok=True)
    bericht = []
    startzeit = time.time()

    # ------------------------------------------------------------ 1. Einlesen
    sequenzen, herkunft = lese_ordner(ordner)
    sequenzen = {n: s for n, s in sequenzen.items() if len(s) > 0}
    namen = list(sequenzen)
    if not namen:
        print(f"⚠️ Keine FASTA-Sequenzen in '{ordner}' gefunden. Erlaubte Endungen: .fasta .fa .fna .fas .ffn .txt")
        return None

    _ausgabe(bericht, "=" * 78)
    _ausgabe(bericht, f"1. ÜBERSICHT: {len(namen)} Sequenzen aus '{ordner}'")
    _ausgabe(bericht, "=" * 78)
    _ausgabe(bericht, f"{'Name':<32} {'Datei':<22} {'Länge':>11} {'GC %':>6} {'N':>7}")
    for n in namen:
        s = sequenzen[n]
        _ausgabe(bericht, f"{n[:32]:<32} {herkunft[n][:22]:<22} {len(s):>11,} {gc_gehalt(s):>6.1f} {s.count('N'):>7,}")
    with open(os.path.join(ausgabe_ordner, "uebersicht.csv"), "w") as f:
        f.write("name;datei;laenge;gc_prozent;A;C;G;T;N\n")
        for n in namen:
            z = zusammensetzung(sequenzen[n])
            f.write(f"{n};{herkunft[n]};{len(sequenzen[n])};{gc_gehalt(sequenzen[n]):.2f};"
                    f"{z['A']};{z['C']};{z['G']};{z['T']};{z['N']}\n")

    # ------------------------------------------------------------ 2. ORFs
    _ausgabe(bericht)
    _ausgabe(bericht, "=" * 78)
    _ausgabe(bericht, f"2. ORF-SUCHE (6 Leseraster, Start {'/'.join(E['STARTCODONS'])}, "
                      f"mind. {E['ORF_MIN_AMINOSAEUREN']} Aminosäuren)")
    _ausgabe(bericht, "=" * 78)
    with open(os.path.join(ausgabe_ordner, "orfs.csv"), "w") as csv, \
            open(os.path.join(ausgabe_ordner, "orf_proteine.fasta"), "w") as faa:
        csv.write("sequenz;nr;strang;rahmen;start;ende;laenge_bp;laenge_as\n")
        for n in namen:
            orfs = finde_orfs(sequenzen[n], E["STARTCODONS"], min_as=E["ORF_MIN_AMINOSAEUREN"])
            for nr, o in enumerate(orfs, 1):
                csv.write(f"{n};{nr};{o['strang']};{o['rahmen']};{o['start']};{o['ende']};"
                          f"{len(o['dna'])};{len(o['protein'])}\n")
                faa.write(f">{n}_ORF{nr} {o['strang']}{o['rahmen']} {o['start']}-{o['ende']} {len(o['protein'])}aa\n")
                for i in range(0, len(o["protein"]), 70):
                    faa.write(o["protein"][i:i + 70] + "\n")
            laengste = max((len(o["protein"]) for o in orfs), default=0)
            abdeckung = sum(len(o["dna"]) for o in orfs) / len(sequenzen[n]) * 100
            _ausgabe(bericht, f"{n[:32]:<32} {len(orfs):>6} ORFs   längster: {laengste:>6} AS   "
                              f"Anteil der Sequenz in ORFs: {min(abdeckung, 100):.0f} %")
    _ausgabe(bericht, "→ Alle ORFs: orfs.csv, Proteinsequenzen: orf_proteine.fasta")

    # ------------------------------------------------------------ 3. Paarweise Alignments
    _ausgabe(bericht)
    _ausgabe(bericht, "=" * 78)
    _ausgabe(bericht, "3. PAARWEISER VERGLEICH (globales Alignment, Needleman-Wunsch)")
    _ausgabe(bericht, "=" * 78)
    paare = [(namen[i], namen[j]) for i in range(len(namen)) for j in range(i + 1, len(namen))]
    alignierbar = [p for p in paare
                   if len(sequenzen[p[0]]) <= E["MAX_ALIGN_LAENGE"] and len(sequenzen[p[1]]) <= E["MAX_ALIGN_LAENGE"]]
    zellen = sum(len(sequenzen[a]) * len(sequenzen[b]) for a, b in alignierbar)
    distanzen_aln = {}

    if not paare:
        _ausgabe(bericht, "Nur eine Sequenz, kein Vergleich möglich.")
    elif zellen > E["MAX_ALIGN_ZELLEN_GESAMT"]:
        _ausgabe(bericht, f"⚠️ Übersprungen: {len(alignierbar)} Alignments wären zusammen {zellen:,} Rechenschritte "
                          f"(Limit {E['MAX_ALIGN_ZELLEN_GESAMT']:,}). Weniger oder kürzere Sequenzen verwenden, "
                          f"oder MAX_ALIGN_ZELLEN_GESAMT erhöhen.")
        alignierbar = []
    if len(alignierbar) < len(paare) and paare and alignierbar:
        _ausgabe(bericht, f"Hinweis: {len(paare) - len(alignierbar)} Paare übersprungen "
                          f"(Sequenz länger als {E['MAX_ALIGN_LAENGE']:,} bp).")

    if alignierbar:
        k_test = 8
        mengen = {n: kmer_menge(sequenzen[n], k_test) for n in namen if len(sequenzen[n]) <= E["MAX_ALIGN_LAENGE"]}
        _ausgabe(bericht, f"{'Sequenz 1':<26} {'Sequenz 2':<26} {'Strang':<7} {'Identität':>9} {'Lücken':>7} {'Score':>7}")
        mutationszeilen = []
        with open(os.path.join(ausgabe_ordner, "alignments.txt"), "w") as f:
            for nr, (a, b) in enumerate(alignierbar, 1):
                s1, s2 = sequenzen[a], sequenzen[b]
                # Liegt Sequenz 2 auf dem Gegenstrang? Dann passt ihr Reverse Complement besser.
                rc = reverse_complement(s2)
                strang = "+"
                if len(mengen[a] & kmer_menge(rc, k_test)) > len(mengen[a] & mengen[b]):
                    s2, strang = rc, "−(rc)"
                erg = needleman_wunsch(s1, s2, E["MATCH"], E["MISMATCH"], E["GAP"])
                st = alignment_statistik(erg["a1"], erg["a2"])
                distanzen_aln[(a, b)] = distanz_aus_alignment(erg["a1"], erg["a2"])
                _ausgabe(bericht, f"{a[:26]:<26} {b[:26]:<26} {strang:<7} {st['identitaet']:>8.1f}% "
                                  f"{st['luecken']:>7} {erg['score']:>7}")
                f.write(f"### {a}  vs  {b}  (Strang {strang})\n")
                f.write(f"Score {erg['score']}, Identität {st['identitaet']:.1f} %, Mismatches {st['mismatches']}, "
                        f"Lücken {st['luecken']}\n\n")
                f.write(alignment_als_text(erg) + "\n\n")
                if len(s1) == len(s2) and st["luecken"] == 0:
                    for m in mutationen(s1, s2):
                        mutationszeilen.append(f"{a};{b};{m['position']};{m['von']};{m['zu']};{m['typ']};"
                                               f"{m['codon']};{m['aminosaeure']};{m['wirkung']}\n")
        if mutationszeilen:
            with open(os.path.join(ausgabe_ordner, "mutationen.csv"), "w") as f:
                f.write("sequenz1;sequenz2;position;von;zu;typ;codon;aminosaeure;wirkung\n")
                f.writelines(mutationszeilen)
            _ausgabe(bericht, f"→ {len(mutationszeilen)} Punktmutationen zwischen gleich langen Sequenzen: mutationen.csv")
        _ausgabe(bericht, "→ Vollständige Alignments: alignments.txt")
        _ausgabe(bericht, "   Achtung: Ein globales Alignment findet IMMER ein Ergebnis. Bei DNA sind Identitäten unter")
        _ausgabe(bericht, "   etwa 50 % mit vielen Lücken auch zwischen zufälligen Sequenzen üblich, also kein Hinweis auf Verwandtschaft.")
        _ausgabe(bericht, "   Strang −(rc): Sequenz 2 passte als Reverse Complement besser und wurde so aligniert.")

    # ------------------------------------------------------------ 4. Stammbaum
    _ausgabe(bericht)
    _ausgabe(bericht, "=" * 78)
    _ausgabe(bericht, "4. STAMMBAUM")
    _ausgabe(bericht, "=" * 78)
    baeume = None
    if len(namen) < 3:
        _ausgabe(bericht, "Für einen Stammbaum braucht man mindestens 3 Sequenzen.")
    else:
        D = {n: {n: 0.0} for n in namen}
        if len(distanzen_aln) == len(paare):
            methode = "Jukes-Cantor-Distanz aus den Alignments"
            for (a, b), d in distanzen_aln.items():
                D[a][b] = D[b][a] = d
        else:
            k = E["KMER_K"] or _k_automatisch(min(len(s) for s in sequenzen.values()))
            methode = f"k-mer-Distanz (alignment-frei, k = {k}), weil nicht alle Paare aligniert werden konnten"
            mengen = {n: kmer_menge(sequenzen[n], k) for n in namen}
            for a, b in paare:
                rc_menge = kmer_menge(reverse_complement(sequenzen[b]), k)
                d = min(kmer_distanz(mengen[a], mengen[b], k), kmer_distanz(mengen[a], rc_menge, k))
                D[a][b] = D[b][a] = d
        endlich = [D[a][b] for a, b in paare if D[a][b] != float("inf")]
        if len(endlich) < len(paare):
            ersatz = 2 * max(endlich) if endlich else 1.0
            _ausgabe(bericht, f"⚠️ Manche Paare sind so verschieden, dass keine Distanz schätzbar ist. "
                              f"Sie bekommen den Ersatzwert {ersatz:.3f}. Der Baum ist dort nicht aussagekräftig.")
            for a, b in paare:
                if D[a][b] == float("inf"):
                    D[a][b] = D[b][a] = ersatz
        _ausgabe(bericht, f"Methode: {methode}")
        laengen = [len(sequenzen[n]) for n in namen]
        if max(laengen) > 3 * min(laengen):
            _ausgabe(bericht, "⚠️ Die Sequenzen sind sehr unterschiedlich lang. Ein Stammbaum ist nur sinnvoll, wenn alle")
            _ausgabe(bericht, "   Sequenzen dieselbe Region (z. B. dasselbe Gen) aus verschiedenen Organismen sind.")

        with open(os.path.join(ausgabe_ordner, "distanzmatrix.csv"), "w") as f:
            f.write(";" + ";".join(namen) + "\n")
            for a in namen:
                f.write(a + ";" + ";".join(f"{D[a][b]:.5f}" for b in namen) + "\n")

        aussen = E["AUSSENGRUPPE"]
        if aussen and not set(aussen) <= set(namen):
            _ausgabe(bericht, f"⚠️ Außengruppe {aussen} nicht in den Daten gefunden, NJ wird am ersten Namen bewurzelt.")
            aussen = None
        b_upgma = upgma(namen, D)
        b_nj = neighbor_joining(namen, D, aussen)
        baeume = {"upgma": b_upgma, "nj": b_nj}

        _ausgabe(bericht, "\nUPGMA:")
        _ausgabe(bericht, baum_als_text(b_upgma))
        _ausgabe(bericht, "\nNeighbor-Joining" + (f" (Außengruppe: {', '.join(aussen)})" if aussen else
                                                  f" (bewurzelt an {namen[0]}, keine Außengruppe angegeben)") + ":")
        _ausgabe(bericht, baum_als_text(b_nj))
        with open(os.path.join(ausgabe_ordner, "baum_upgma.nwk"), "w") as f:
            f.write(newick(b_upgma) + ";\n")
        with open(os.path.join(ausgabe_ordner, "baum_nj.nwk"), "w") as f:
            f.write(newick(b_nj) + ";\n")
        _ausgabe(bericht, "→ Distanzmatrix: distanzmatrix.csv, Bäume im Newick-Format: baum_upgma.nwk, baum_nj.nwk")

        try:
            import matplotlib.pyplot as plt
            hoehe = max(4, 0.35 * len(namen))
            fig, achsen = plt.subplots(1, 2, figsize=(15, hoehe))
            zeichne_baum(b_upgma, "UPGMA", achsen[0])
            zeichne_baum(b_nj, "Neighbor-Joining", achsen[1])
            plt.tight_layout()
            fig.savefig(os.path.join(ausgabe_ordner, "stammbaum.png"), dpi=150)
            _ausgabe(bericht, "→ Grafik: stammbaum.png")
            if E["ZEIGE_GRAFIK"]:
                plt.show()
            else:
                plt.close(fig)
        except ImportError:
            _ausgabe(bericht, "(matplotlib fehlt, keine Grafik)")

    _ausgabe(bericht)
    _ausgabe(bericht, f"Fertig in {time.time() - startzeit:.1f} s. Alle Ergebnisse liegen in '{ausgabe_ordner}/'.")
    with open(os.path.join(ausgabe_ordner, "bericht.txt"), "w") as f:
        f.write("\n".join(bericht) + "\n")
    return {"sequenzen": sequenzen, "baeume": baeume}


if __name__ == "__main__":
    import sys
    analysiere_ordner(sys.argv[1] if len(sys.argv) > 1 else "daten",
                      sys.argv[2] if len(sys.argv) > 2 else "ergebnisse")
