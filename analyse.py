"""
📂 analyse.py – wertet automatisch alle FASTA-Dateien in einem Ordner aus.

Für JEDE Sequenz:      Länge, Basenzusammensetzung, GC-Gehalt, ORFs (mögliche Gene) + Proteine
Für JEDES Paar:        Alignment (Needleman-Wunsch), Identität, Lücken
Gegen eine REFERENZ:   Variantenanalyse: Mutationen (auch Insertionen/Deletionen), Transkription + Translation,
                       Wirkung auf das Protein (stumm, Missense, Nonsense, Frameshift ...)
Für ALLE zusammen:     Distanzmatrix und Stammbaum (UPGMA + Neighbor-Joining), ab 3 Sequenzen

Große Sequenzen: Statistik und ORF-Suche gehen auch für ganze Bakteriengenome (Millionen Basen).
Alignments brauchen Zeit ~ Länge1 × Länge2 und werden deshalb nur bis MAX_ALIGN_LAENGE gemacht.
Für den Stammbaum wird dann automatisch eine alignment-freie k-mer-Distanz benutzt.

Benutzung:
    from analyse import analysiere_ordner, analysiere_alle
    analysiere_ordner("daten")                 # ein Ordner → Textbericht + Dateien + bericht.html
    analysiere_alle("beispiele")               # jeder Unterordner einzeln → EIN gemeinsamer bericht.html

Ein Ordner kann zusätzlich enthalten:
    einstellungen.json   eigene Einstellungen für diesen Ordner, z. B. {"REFERENZ": "...", "AUSSENGRUPPE": [...]}
    accessions.txt       NCBI-Nummern (eine pro Zeile); fehlende Sequenzen werden automatisch heruntergeladen
    LIESMICH.md          kurze Beschreibung, erscheint im HTML-Bericht
"""

import json
import os
import time

from werkzeuge import FASTA_ENDUNGEN
from werkzeuge import (lese_ordner, gc_gehalt, zusammensetzung, finde_orfs, ausrichten,
                       alignment_statistik, alignment_als_text, reverse_complement,
                       distanz_aus_alignment, kmer_menge, kmer_distanz, upgma, neighbor_joining,
                       newick, baum_als_text, transkribiere, uebersetze, finde_cds,
                       finde_varianten, wirkung_der_varianten, variante_auf_referenz_position,
                       vergleiche_proteine, lade_von_ncbi, normalisiere_indels)

STANDARD = {
    "ORF_MIN_AMINOSAEUREN": 100,          # kürzere ORFs ignorieren
    "STARTCODONS": ("ATG",),              # für Bakterien realistischer: ("ATG", "GTG", "TTG")
    "MAX_ALIGN_LAENGE": 20_000,           # längere Sequenzen werden nicht aligniert (z. B. mtDNA mit 16.569 bp geht noch)
    "MAX_ALIGN_ZELLEN_GESAMT": 1_000_000_000,  # Summe aller Länge1×Länge2; ca. 40 Sekunden in Colab
    # Alignment-Punkte wie bei NCBI BLASTN; eine Lücke der Länge L kostet LUECKE_OEFFNEN + L × LUECKE_VERLAENGERN
    "MATCH": 2, "MISMATCH": -3, "LUECKE_OEFFNEN": 5, "LUECKE_VERLAENGERN": 2,
    "KMER_K": None,                       # None = automatisch nach Sequenzlänge
    "AUSSENGRUPPE": None,                 # Liste von Namen für den NJ-Baum
    "REFERENZ": None,                     # Referenz für die Variantenanalyse (None = erste Sequenz)
    "CDS": None,                          # (start, ende) des Gens in der Referenz, ab 1; None = automatisch
    "MIN_IDENTITAET_VARIANTEN": 80,       # nur Sequenzen, die der Referenz mindestens so ähnlich sind (%)
    "ZEIGE_GRAFIK": True,                 # Stammbaum direkt anzeigen
    "TEXT_AUSGABE": True,                 # Textbericht ausgeben
}


def _k_automatisch(min_laenge):
    if min_laenge >= 100_000:
        return 21          # Standard von Mash für Genome
    if min_laenge >= 5_000:
        return 15
    return 11


def _einstellungen(einstellungen, ordner=None):
    E = dict(STANDARD)
    einstellungen = {k: v for k, v in einstellungen.items() if v is not None}   # None = "nicht gesetzt"
    if ordner and os.path.exists(os.path.join(ordner, "einstellungen.json")):
        with open(os.path.join(ordner, "einstellungen.json")) as f:
            einstellungen = {**json.load(f), **einstellungen}          # Werte aus dem Aufruf haben Vorrang
    for schluessel, wert in einstellungen.items():
        if schluessel not in E:
            raise ValueError(f"Unbekannte Einstellung: {schluessel}. Möglich: {', '.join(E)}")
        E[schluessel] = tuple(wert) if schluessel in ("STARTCODONS", "CDS") and wert else wert
    return E


def _aligniere(s1, s2, E):
    """Globales Alignment mit affiner Lückenstrafe (Gotoh).
    Unterschiedlich lange mRNA-Enden werden so als EINE Lücke am Ende erklärt; die Identität
    wird danach nur über den überlappenden Bereich berechnet (siehe alignment_statistik)."""
    return ausrichten(s1, s2, "global", E["MATCH"], E["MISMATCH"], E["LUECKE_OEFFNEN"],
                      E["LUECKE_VERLAENGERN"])


class _Bericht:
    def __init__(self, drucken):
        self.zeilen = []
        self.drucken = drucken

    def __call__(self, text=""):
        if self.drucken:
            print(text)
        self.zeilen.append(text)

    def titel(self, text):
        self("")
        self("=" * 78)
        self(text)
        self("=" * 78)


# ======================================================================
# Hauptfunktion für einen Ordner
# ======================================================================

def analysiere_ordner(ordner, ausgabe_ordner="ergebnisse", html=True, **einstellungen):
    E = _einstellungen(einstellungen, ordner)
    os.makedirs(ausgabe_ordner, exist_ok=True)
    B = _Bericht(E["TEXT_AUSGABE"])
    startzeit = time.time()
    R = {"titel": os.path.basename(os.path.normpath(ordner)), "ordner": ordner, "ausgabe_ordner": ausgabe_ordner,
         "einstellungen": E, "beschreibung": "", "warnungen": []}
    liesmich = os.path.join(ordner, "LIESMICH.md")
    if os.path.exists(liesmich):
        with open(liesmich) as f:
            R["beschreibung"] = f.read().strip()

    # ------------------------------------------------------------ 0. Fehlende Daten von NCBI holen
    acc_datei = os.path.join(ordner, "accessions.txt")
    if os.path.exists(acc_datei):
        with open(acc_datei) as f:
            accessions = [z.split("#")[0].strip() for z in f if z.split("#")[0].strip()]
        fehlend = [a for a in accessions if not os.path.exists(os.path.join(ordner, a + ".fasta"))]
        if fehlend:
            B(f"Lade {len(fehlend)} Sequenz(en) von NCBI: {', '.join(fehlend)}")
            for a, fehler in lade_von_ncbi(fehlend, ordner):
                if fehler:
                    R["warnungen"].append(f"{a} konnte nicht geladen werden: {fehler}")
                    B(f"⚠️ {a} konnte nicht geladen werden: {fehler}")

    # ------------------------------------------------------------ 1. Einlesen
    sequenzen, herkunft = lese_ordner(ordner)
    sequenzen = {n: s for n, s in sequenzen.items() if len(s) > 0}
    namen = list(sequenzen)
    R["namen"] = namen
    if not namen:
        B(f"⚠️ Keine FASTA-Sequenzen in '{ordner}' gefunden. Erlaubte Endungen: .fasta .fa .fna .fas .ffn .txt")
        R["sequenzen"] = []
        return R

    B.titel(f"1. ÜBERSICHT: {len(namen)} Sequenz(en) aus '{ordner}'")
    B(f"{'Name':<48} {'Datei':<22} {'Länge':>11} {'GC %':>6} {'N':>7}")
    R["sequenzen"] = []
    for n in namen:
        s = sequenzen[n]
        z = zusammensetzung(s)
        R["sequenzen"].append({"name": n, "datei": herkunft[n], "laenge": len(s), "gc": gc_gehalt(s), **z})
        B(f"{n[:48]:<48} {herkunft[n][:22]:<22} {len(s):>11,} {gc_gehalt(s):>6.1f} {s.count('N'):>7,}")
    with open(os.path.join(ausgabe_ordner, "uebersicht.csv"), "w") as f:
        f.write("name;datei;laenge;gc_prozent;A;C;G;T;N\n")
        for x in R["sequenzen"]:
            f.write(f"{x['name']};{x['datei']};{x['laenge']};{x['gc']:.2f};{x['A']};{x['C']};{x['G']};{x['T']};{x['N']}\n")

    # ------------------------------------------------------------ 2. ORFs
    B.titel(f"2. ORF-SUCHE (6 Leseraster, Start {'/'.join(E['STARTCODONS'])}, mind. {E['ORF_MIN_AMINOSAEUREN']} Aminosäuren)")
    with open(os.path.join(ausgabe_ordner, "orfs.csv"), "w") as csv, \
            open(os.path.join(ausgabe_ordner, "orf_proteine.fasta"), "w") as faa:
        csv.write("sequenz;nr;strang;rahmen;start;ende;laenge_bp;laenge_as\n")
        for x in R["sequenzen"]:
            n = x["name"]
            orfs = finde_orfs(sequenzen[n], E["STARTCODONS"], min_as=E["ORF_MIN_AMINOSAEUREN"])
            x["orfs"] = [{k: o[k] for k in ("strang", "rahmen", "start", "ende")} | {"laenge_as": len(o["protein"])}
                         for o in orfs]
            for nr, o in enumerate(orfs, 1):
                csv.write(f"{n};{nr};{o['strang']};{o['rahmen']};{o['start']};{o['ende']};"
                          f"{len(o['dna'])};{len(o['protein'])}\n")
                faa.write(f">{n}_ORF{nr} {o['strang']}{o['rahmen']} {o['start']}-{o['ende']} {len(o['protein'])}aa\n")
                for i in range(0, len(o["protein"]), 70):
                    faa.write(o["protein"][i:i + 70] + "\n")
            x["laengster_orf"] = max((len(o["protein"]) for o in orfs), default=0)
            x["orf_anteil"] = min(100.0, sum(len(o["dna"]) for o in orfs) / len(sequenzen[n]) * 100)
            B(f"{n[:48]:<48} {len(orfs):>6} ORFs   längster: {x['laengster_orf']:>6} AS   "
              f"Anteil der Sequenz in ORFs: {x['orf_anteil']:.0f} %")
    B("→ Alle ORFs: orfs.csv, Proteinsequenzen: orf_proteine.fasta")

    # ------------------------------------------------------------ 3. Paarweise Alignments
    B.titel("3. PAARWEISER VERGLEICH (globales Alignment mit affiner Lückenstrafe)")
    paare = [(namen[i], namen[j]) for i in range(len(namen)) for j in range(i + 1, len(namen))]
    alignierbar = [p for p in paare
                   if len(sequenzen[p[0]]) <= E["MAX_ALIGN_LAENGE"] and len(sequenzen[p[1]]) <= E["MAX_ALIGN_LAENGE"]]
    zellen = sum(len(sequenzen[a]) * len(sequenzen[b]) for a, b in alignierbar)
    R["paare"] = []
    R["alignment_hinweis"] = ""
    identitaeten = {}
    distanzen_aln = {}

    if not paare:
        R["alignment_hinweis"] = "Nur eine Sequenz, kein Vergleich möglich."
    elif zellen > E["MAX_ALIGN_ZELLEN_GESAMT"]:
        R["alignment_hinweis"] = (f"Übersprungen: Die Alignments wären zusammen {zellen:,} Rechenschritte "
                                  f"(Limit {E['MAX_ALIGN_ZELLEN_GESAMT']:,}).")
        alignierbar = []
    elif len(alignierbar) < len(paare):
        R["alignment_hinweis"] = (f"{len(paare) - len(alignierbar)} Paar(e) übersprungen, weil eine Sequenz länger "
                                  f"als {E['MAX_ALIGN_LAENGE']:,} bp ist.")
    if R["alignment_hinweis"]:
        B(R["alignment_hinweis"])

    if alignierbar:
        k_test = 8
        mengen = {n: kmer_menge(sequenzen[n], k_test) for n in namen if len(sequenzen[n]) <= E["MAX_ALIGN_LAENGE"]}
        B(f"{'Sequenz 1':<40} {'Sequenz 2':<40} {'Strang':<7} {'Identität':>9} {'Lücken':>7}")
        with open(os.path.join(ausgabe_ordner, "alignments.txt"), "w") as f:
            for a, b in alignierbar:
                s1, s2 = sequenzen[a], sequenzen[b]
                rc = reverse_complement(s2)
                strang = "+"
                if len(mengen[a] & kmer_menge(rc, k_test)) > len(mengen[a] & mengen[b]):
                    s2, strang = rc, "−"
                erg = _aligniere(s1, s2, E)
                st = alignment_statistik(erg["a1"], erg["a2"])
                d = distanz_aus_alignment(erg["a1"], erg["a2"])
                distanzen_aln[(a, b)] = d
                identitaeten[(a, b)] = identitaeten[(b, a)] = st["identitaet"]
                R["paare"].append({"a": a, "b": b, "strang": strang, "identitaet": st["identitaet"],
                                   "luecken": st["luecken"], "mismatches": st["mismatches"], "score": erg["score"],
                                   "distanz": d})
                B(f"{a[:40]:<40} {b[:40]:<40} {strang:<7} {st['identitaet']:>8.1f}% {st['luecken']:>7}")
                f.write(f"### {a}  vs  {b}  (Strang {strang})\n")
                f.write(f"Score {erg['score']}, Identität {st['identitaet']:.1f} %, Mismatches {st['mismatches']}, "
                        f"Lücken {st['luecken']}\n\n")
                f.write(alignment_als_text(erg) + "\n\n")
        B("→ Vollständige Alignments: alignments.txt")
        B("   Strang −: Sequenz 2 passte als Reverse Complement (Gegenstrang) besser und wurde so aligniert.")
        B("   Achtung: Ein globales Alignment findet IMMER ein Ergebnis. Bei DNA sind Identitäten unter etwa")
        B("   50 % mit vielen Lücken auch zwischen zufälligen Sequenzen üblich, also kein Hinweis auf Verwandtschaft.")

    # ------------------------------------------------------------ 4. Variantenanalyse
    B.titel("4. VARIANTENANALYSE: Mutationen gegenüber einer Referenz und ihre Wirkung")
    R["varianten"] = _variantenanalyse(sequenzen, namen, identitaeten, E, ausgabe_ordner, B)

    # ------------------------------------------------------------ 5. Stammbaum
    B.titel("5. STAMMBAUM")
    R["baum"] = _stammbaum(sequenzen, namen, paare, distanzen_aln, E, ausgabe_ordner, B)

    R["dauer"] = time.time() - startzeit
    B("")
    B(f"Fertig in {R['dauer']:.1f} s. Alle Ergebnisse liegen in '{ausgabe_ordner}/'.")
    R["text"] = "\n".join(B.zeilen)
    with open(os.path.join(ausgabe_ordner, "bericht.txt"), "w") as f:
        f.write(R["text"] + "\n")
    if html:
        from bericht import erstelle_html_bericht
        pfad = erstelle_html_bericht([R], os.path.join(ausgabe_ordner, "bericht.html"))
        if E["TEXT_AUSGABE"]:
            print(f"→ Bericht mit Tabellen und Diagrammen: {pfad}")
    return R


# ======================================================================
# Variantenanalyse
# ======================================================================

def _variantenanalyse(sequenzen, namen, identitaeten, E, ausgabe_ordner, B):
    V = {"referenz": None, "cds": None, "ergebnisse": [], "hinweis": "", "uebersprungen": []}
    if len(namen) < 2:
        V["hinweis"] = "Braucht mindestens 2 Sequenzen (Referenz + Variante)."
        B(V["hinweis"])
        return V
    ref_name = E["REFERENZ"] or namen[0]
    if ref_name not in sequenzen:
        V["hinweis"] = f"Referenz '{ref_name}' nicht gefunden. Vorhandene Namen: {', '.join(namen)}"
        B("⚠️ " + V["hinweis"])
        return V
    referenz = sequenzen[ref_name]
    kandidaten = [n for n in namen if n != ref_name and identitaeten.get((ref_name, n)) is not None
                  and identitaeten[(ref_name, n)] >= E["MIN_IDENTITAET_VARIANTEN"]]
    V["uebersprungen"] = [n for n in namen if n != ref_name and n not in kandidaten]
    V["referenz"] = ref_name
    if not kandidaten:
        V["hinweis"] = (f"Keine Sequenz ist der Referenz zu mindestens {E['MIN_IDENTITAET_VARIANTEN']} % ähnlich "
                        f"(oder die Sequenzen sind zu lang zum Alignieren). Die Variantenanalyse ist für Varianten "
                        f"DESSELBEN Gens gedacht, z. B. gesund gegen Mutation oder dasselbe Gen nah verwandter Arten.")
        B("Übersprungen: " + V["hinweis"])
        return V

    cds = tuple(E["CDS"]) if E["CDS"] else finde_cds(referenz)
    B(f"Referenz: {ref_name} ({len(referenz):,} bp)")
    if cds and not E["CDS"] and (cds[1] - cds[0] + 1) < 0.5 * len(referenz):
        V["hinweis"] = (f"Kein eindeutiges Gen erkannt: Der längste ORF ({cds[0]}–{cds[1]}) deckt weniger als die "
                        f"Hälfte der Referenz ab (nicht codierende DNA oder mehrere Gene). Es werden nur die "
                        f"Positionen gezeigt. Falls es doch ein Gen ist: CDS = (start, ende) einstellen.")
        B(V["hinweis"])
        cds = None
    elif cds is None:
        V["hinweis"] = "Kein Gen (ORF ab 30 Aminosäuren) in der Referenz gefunden. Es werden nur die Positionen gezeigt."
        B(V["hinweis"])

    protein = mrna = cds_dna = ""
    if cds:
        cds_dna = referenz[cds[0] - 1:cds[1]]
        mrna = transkribiere(cds_dna)
        protein = uebersetze(cds_dna)
        V.update({"cds": cds, "cds_quelle": "eingestellt" if E["CDS"] else "automatisch (längster ORF)",
                  "dna": cds_dna, "mrna": mrna, "protein": protein})
        B(f"Gen (CDS): Position {cds[0]}–{cds[1]} ({V['cds_quelle']}), {len(protein)} Aminosäuren")
        B("")
        B("Vom Gen zum Protein (Referenz):")
        B(f"  DNA  (codierender Strang)  5'-{cds_dna[:45]}{'...' if len(cds_dna) > 45 else ''}-3'")
        B(f"  mRNA (Transkription T→U)   5'-{mrna[:45]}{'...' if len(mrna) > 45 else ''}-3'")
        codons = " ".join(mrna[i:i + 3] for i in range(0, min(len(mrna), 30), 3))
        aminos = " ".join(f"{a:^3}" for a in protein[:10])
        B(f"  Codons                        {codons}{' ...' if len(mrna) > 30 else ''}")
        B(f"  Protein (Translation)         {aminos}{' ...' if len(protein) > 10 else ''}")
    if V["uebersprungen"]:
        B(f"\nNicht analysiert (unter {E['MIN_IDENTITAET_VARIANTEN']} % Identität oder zu lang): "
          f"{', '.join(V['uebersprungen'])}")

    with open(os.path.join(ausgabe_ordner, "varianten.csv"), "w") as csv, \
            open(os.path.join(ausgabe_ordner, "proteine_varianten.fasta"), "w") as faa, \
            open(os.path.join(ausgabe_ordner, "mrna_varianten.fasta"), "w") as fmrna:
        csv.write("referenz;variante;art;position_referenz;ref;alt;c_notation;p_notation;wirkung\n")
        if cds:
            faa.write(f">{ref_name} (Referenz)\n{protein}\n")
            fmrna.write(f">{ref_name} (Referenz)\n{mrna}\n")
        for name in kandidaten:
            variante = sequenzen[name]
            rc = reverse_complement(variante)
            ref_k = kmer_menge(referenz, 8)
            if len(ref_k & kmer_menge(rc, 8)) > len(ref_k & kmer_menge(variante, 8)):
                variante = rc
            erg = _aligniere(referenz, variante, E)
            varianten = normalisiere_indels(referenz, finde_varianten(erg["a1"], erg["a2"]))
            if cds:
                wirkung_der_varianten(referenz, varianten, cds)
            eintrag = {"name": name, "varianten": varianten, "identitaet": identitaeten[(ref_name, name)]}
            B(f"\n▶ {name}: {len(varianten)} Unterschied(e) zur Referenz")
            for nr_v, v in enumerate(varianten):
                v["position"] = (f"{v['start']}" if v["start"] == v["ende"] or v["art"] == "Insertion"
                                 else f"{v['start']}-{v['ende']}")
                if v["art"] == "SNV":
                    v["aenderung"] = f"{v['ref']}→{v['alt']}"
                elif v["art"] == "Deletion":
                    v["aenderung"] = f"−{v['ref'][:12]}{'…' if len(v['ref']) > 12 else ''} ({len(v['ref'])} bp)"
                else:
                    v["aenderung"] = f"+{v['alt'][:12]}{'…' if len(v['alt']) > 12 else ''} ({len(v['alt'])} bp)"
                if nr_v < 25:
                    B(f"   {v['art']:<9} Pos. {v['position']:<11} {v['aenderung']:<22} {v.get('c', ''):<18} "
                      f"{v.get('p', ''):<14} {v.get('wirkung', '')}")
                elif nr_v == 25:
                    B(f"   ... und {len(varianten) - 25} weitere (alle in varianten.csv)")
                csv.write(f"{ref_name};{name};{v['art']};{v['position']};{v['ref']};{v['alt']};{v.get('c', '')};"
                          f"{v.get('p', '')};{v.get('wirkung', '')}\n")
            if cds:
                vstart = variante_auf_referenz_position(erg["a1"], erg["a2"], cds[0])
                if vstart is None:
                    eintrag["protein_vergleich"] = "Startcodon fehlt in der Variante, kein Protein bestimmbar."
                    eintrag["protein"] = ""
                else:
                    v_protein = uebersetze(variante[vstart - 1:])
                    v_mrna = transkribiere(variante[vstart - 1:vstart - 1 + 3 * (len(v_protein) + 1)])
                    eintrag["protein"] = v_protein
                    eintrag["protein_vergleich"] = vergleiche_proteine(protein, v_protein)
                    faa.write(f">{name}\n{v_protein}\n")
                    fmrna.write(f">{name}\n{v_mrna}\n")
                B(f"   Protein: {eintrag['protein_vergleich']}")
            V["ergebnisse"].append(eintrag)
    B("\n→ Alle Varianten: varianten.csv, Proteine: proteine_varianten.fasta, mRNAs: mrna_varianten.fasta")
    B("   c. = Position im Gen (c.1 = A des Startcodons), p. = Aminosäure-Änderung (HGVS-Schreibweise)")
    return V


# ======================================================================
# Stammbaum
# ======================================================================

def _stammbaum(sequenzen, namen, paare, distanzen_aln, E, ausgabe_ordner, B):
    if len(namen) < 3:
        B("Für einen Stammbaum braucht man mindestens 3 Sequenzen.")
        return None
    T = {"warnungen": []}
    D = {n: {n: 0.0} for n in namen}
    if len(distanzen_aln) == len(paare):
        T["methode"] = "Jukes-Cantor-Distanz aus den Alignments"
        for (a, b), d in distanzen_aln.items():
            D[a][b] = D[b][a] = d
    else:
        k = E["KMER_K"] or _k_automatisch(min(len(s) for s in sequenzen.values()))
        T["methode"] = f"k-mer-Distanz (alignment-frei, k = {k}), weil nicht alle Paare aligniert werden konnten"
        mengen = {n: kmer_menge(sequenzen[n], k) for n in namen}
        for a, b in paare:
            rc_menge = kmer_menge(reverse_complement(sequenzen[b]), k)
            D[a][b] = D[b][a] = min(kmer_distanz(mengen[a], mengen[b], k), kmer_distanz(mengen[a], rc_menge, k))
    endlich = [D[a][b] for a, b in paare if D[a][b] != float("inf")]
    if len(endlich) < len(paare):
        ersatz = 2 * max(endlich) if endlich else 1.0
        T["warnungen"].append(f"Manche Paare sind so verschieden, dass keine Distanz schätzbar ist. Sie bekommen den "
                              f"Ersatzwert {ersatz:.3f}; der Baum ist dort nicht aussagekräftig.")
        for a, b in paare:
            if D[a][b] == float("inf"):
                D[a][b] = D[b][a] = ersatz
    laengen = [len(sequenzen[n]) for n in namen]
    if max(laengen) > 3 * min(laengen):
        T["warnungen"].append("Die Sequenzen sind sehr unterschiedlich lang. Ein Stammbaum ist nur sinnvoll, wenn alle "
                              "Sequenzen dieselbe Region (z. B. dasselbe Gen) aus verschiedenen Organismen sind.")
    B(f"Methode: {T['methode']}")
    for w in T["warnungen"]:
        B("⚠️ " + w)
    T["distanzen"] = D

    with open(os.path.join(ausgabe_ordner, "distanzmatrix.csv"), "w") as f:
        f.write(";" + ";".join(namen) + "\n")
        for a in namen:
            f.write(a + ";" + ";".join(f"{D[a][b]:.5f}" for b in namen) + "\n")

    aussen = E["AUSSENGRUPPE"]
    if aussen and not set(aussen) <= set(namen):
        T["warnungen"].append(f"Außengruppe {aussen} nicht in den Daten gefunden; NJ wird am ersten Namen bewurzelt.")
        B("⚠️ " + T["warnungen"][-1])
        aussen = None
    T["aussengruppe"] = aussen
    T["upgma"] = upgma(namen, D)
    T["nj"] = neighbor_joining(namen, D, aussen)
    B("\nUPGMA:")
    B(baum_als_text(T["upgma"]))
    B("\nNeighbor-Joining" + (f" (Außengruppe: {', '.join(aussen)})" if aussen else
                              f" (bewurzelt an {namen[0]}, keine Außengruppe angegeben)") + ":")
    B(baum_als_text(T["nj"]))
    for art in ("upgma", "nj"):
        with open(os.path.join(ausgabe_ordner, f"baum_{art}.nwk"), "w") as f:
            f.write(newick(T[art]) + ";\n")
    B("→ Distanzmatrix: distanzmatrix.csv, Bäume im Newick-Format: baum_upgma.nwk, baum_nj.nwk")

    try:
        from bericht import baum_abbildung
        fig = baum_abbildung(T)
        fig.savefig(os.path.join(ausgabe_ordner, "stammbaum.png"), dpi=150, bbox_inches="tight")
        B("→ Grafik: stammbaum.png")
        if E["ZEIGE_GRAFIK"] and E["TEXT_AUSGABE"]:
            try:
                from IPython.display import display
                display(fig)
            except ImportError:
                pass
    except ImportError:
        B("(matplotlib fehlt, keine Grafik)")
    return T


# ======================================================================
# Mehrere Ordner → ein gemeinsamer Bericht
# ======================================================================

def projekte_im_ordner(eltern_ordner, ausgabe_ordner=None):
    """Liste der Projektordner: jeder Unterordner ist ein Projekt; FASTA-Dateien direkt im
    Ordner bilden zusätzlich ein eigenes Projekt."""
    unterordner = sorted(d for d in os.listdir(eltern_ordner)
                         if os.path.isdir(os.path.join(eltern_ordner, d)) and not d.startswith((".", "_"))
                         and (ausgabe_ordner is None or
                              os.path.abspath(os.path.join(eltern_ordner, d)) != os.path.abspath(ausgabe_ordner)))
    eigene_dateien = [f for f in os.listdir(eltern_ordner)
                      if f.lower().endswith(FASTA_ENDUNGEN) and not f.startswith(".")
                      and f.lower() not in ("liesmich.txt", "readme.txt")]
    if eigene_dateien or not unterordner:
        unterordner = ["."] + unterordner
    return [os.path.normpath(os.path.join(eltern_ordner, d)) for d in unterordner]


def analysiere_projekte(projekt_ordner, ausgabe_ordner="ergebnisse_alle", **einstellungen):
    """Analysiert mehrere Projektordner einzeln und fasst alles in EINEM HTML-Bericht zusammen.
    Gibt (Liste der Ergebnisse, Pfad zum Bericht) zurück."""
    from bericht import erstelle_html_bericht
    ergebnisse = []
    for pfad in projekt_ordner:
        ziel = os.path.join(ausgabe_ordner, os.path.basename(os.path.abspath(pfad)))
        print(f"▶ Analysiere {pfad} ...", end=" ", flush=True)
        start = time.time()
        R = analysiere_ordner(pfad, ziel, html=False, **{"ZEIGE_GRAFIK": False, "TEXT_AUSGABE": False, **einstellungen})
        for w in R.get("warnungen", []):
            print(f"\n   ⚠️ {w}", end="")
        if not R.get("sequenzen"):
            print("\n   keine Sequenzen, übersprungen.")
            if R.get("warnungen"):
                ergebnisse.append(R)
            continue
        print(f"{len(R['sequenzen'])} Sequenz(en), {time.time() - start:.1f} s")
        ergebnisse.append(R)
    pfad = erstelle_html_bericht(ergebnisse, os.path.join(ausgabe_ordner, "bericht.html"))
    print(f"\n✅ Gemeinsamer Bericht: {pfad}")
    print(f"   Einzelergebnisse (CSV, FASTA, Bäume) je Projekt im Ordner '{ausgabe_ordner}/'.")
    return ergebnisse, pfad


def analysiere_alle(eltern_ordner, ausgabe_ordner="ergebnisse_alle", **einstellungen):
    """Analysiert jeden Unterordner von eltern_ordner einzeln und fasst alles in EINEM HTML-Bericht zusammen."""
    return analysiere_projekte(projekte_im_ordner(eltern_ordner, ausgabe_ordner), ausgabe_ordner, **einstellungen)


if __name__ == "__main__":
    import sys
    analysiere_ordner(sys.argv[1] if len(sys.argv) > 1 else "daten",
                      sys.argv[2] if len(sys.argv) > 2 else "ergebnisse")
