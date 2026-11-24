# 🧬 ORF-Finder: mögliche Gene in einer DNA-Sequenz finden
#
# Ein ORF (Open Reading Frame, offenes Leseraster) ist ein Abschnitt, der mit einem
# Startcodon beginnt und ohne Unterbrechung bis zu einem Stoppcodon läuft.
# Ein Gen kann auf jedem der beiden Stränge liegen und in jedem von 3 Leserastern
# beginnen → 6 Leseraster werden durchsucht.
#
# So benutzt du es: Ganzen Code in eine Colab-Zelle kopieren und ausführen.
# Standardmäßig wird das Plasmid pPCP1 geladen und die gefundenen ORFs werden mit
# den 10 bei NCBI annotierten Genen verglichen.

import os
import urllib.request

# ============================================================
# EINSTELLUNGEN
# ============================================================
MIN_AMINOSAEUREN = 60            # kürzere ORFs werden ignoriert (zufällige ORFs sind meist kurz)
STARTCODONS = ["ATG"]            # für Bakterien realistischer: ["ATG", "GTG", "TTG"]
STOPPCODONS = ["TAA", "TAG", "TGA"]

# Eigene Sequenz statt pPCP1? Dann hier eintragen (sonst leer lassen):
EIGENE_SEQUENZ = ""
# ============================================================

BASIS_URL = "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/"
DATEIEN = {"pPCP1_plasmid.fasta": "NC_005816.fna", "pPCP1_gene.fasta": "NC_005816.ffn"}

basen = "TCAG"
aminosaeuren = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_tabelle = {}
nummer = 0
for b1 in basen:
    for b2 in basen:
        for b3 in basen:
            codon_tabelle[b1 + b2 + b3] = aminosaeuren[nummer]   # hier direkt mit DNA-Codons (T statt U)
            nummer = nummer + 1


def lese_fasta(dateiname):
    sequenzen = {}
    name = None
    with open(dateiname) as f:
        for zeile in f:
            zeile = zeile.strip()
            if zeile == "":
                continue
            if zeile.startswith(">"):
                name = zeile[1:]
                sequenzen[name] = ""
            else:
                sequenzen[name] = sequenzen[name] + zeile.upper()
    return sequenzen


def reverse_complement(dna):
    komplement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    ergebnis = ""
    for base in dna:
        ergebnis = komplement.get(base, "N") + ergebnis
    return ergebnis


def uebersetze(dna):
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        protein = protein + codon_tabelle.get(dna[i:i + 3], "X")
    return protein


def finde_orfs_in_strang(seq, startcodons, stoppcodons, min_as):
    """Sucht ORFs in den 3 Leserastern eines Strangs.
    Gibt Liste mit (rahmen, start, ende) zurück; start/ende als Python-Index (ab 0, ende exklusiv).
    Pro Stoppcodon wird nur der längste ORF genommen (erstes Startcodon)."""
    orfs = []
    for rahmen in range(3):
        start = None
        for i in range(rahmen, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if start is None and codon in startcodons:
                start = i
            elif start is not None and codon in stoppcodons:
                ende = i + 3                         # Stoppcodon gehört zum ORF dazu
                laenge_as = (ende - start) // 3 - 1  # ohne Stopp
                if laenge_as >= min_as:
                    orfs.append((rahmen, start, ende))
                start = None
    return orfs


def finde_orfs(dna, startcodons=STARTCODONS, stoppcodons=STOPPCODONS, min_as=MIN_AMINOSAEUREN):
    """Durchsucht alle 6 Leseraster. Positionen werden biologisch angegeben (ab 1, Ende eingeschlossen)."""
    L = len(dna)
    ergebnisse = []

    for rahmen, s, e in finde_orfs_in_strang(dna, startcodons, stoppcodons, min_as):
        ergebnisse.append({"strang": "+", "rahmen": rahmen + 1,
                           "start": s + 1, "ende": e,
                           "dna": dna[s:e]})

    gegenstrang = reverse_complement(dna)
    for rahmen, s, e in finde_orfs_in_strang(gegenstrang, startcodons, stoppcodons, min_as):
        # Position im Gegenstrang zurückrechnen: Index k im Gegenstrang = Position L - k im Originalstrang (ab 1)
        ergebnisse.append({"strang": "-", "rahmen": rahmen + 1,
                           "start": L - s, "ende": L - e + 1,
                           "dna": gegenstrang[s:e]})

    for orf in ergebnisse:
        orf["protein"] = uebersetze(orf["dna"]).rstrip("*")
    ergebnisse.sort(key=lambda orf: min(orf["start"], orf["ende"]))
    return ergebnisse


def zeige_orfs(orfs):
    print(f"{'Nr':>3}  {'Strang':<6} {'Rahmen':<6} {'Position':<14} {'bp':>5} {'AS':>5}  Protein (Anfang)")
    for nr, orf in enumerate(orfs, start=1):
        if orf["strang"] == "+":
            position = f"{orf['start']}-{orf['ende']}"
        else:
            position = f"c{orf['start']}-{orf['ende']}"     # c = complement, wie bei NCBI
        print(f"{nr:>3}  {orf['strang']:<6} {orf['rahmen']:<6} {position:<14} {len(orf['dna']):>5} "
              f"{len(orf['protein']):>5}  {orf['protein'][:30]}...")


def lade_pPCP1():
    for lokal, online in DATEIEN.items():
        if not os.path.exists(lokal):
            urllib.request.urlretrieve(BASIS_URL + online, lokal)
    plasmid = list(lese_fasta("pPCP1_plasmid.fasta").values())[0]
    gene = lese_fasta("pPCP1_gene.fasta")
    return plasmid, gene


def vergleiche_mit_annotation(orfs, gene):
    """Vergleicht gefundene ORFs mit den NCBI-Genen. Gleiches Ende (Stoppcodon) = gleiches Gen."""
    print("\nVergleich mit den 10 bei NCBI annotierten Genen:")
    print(f"{'NCBI-Position':<14} {'Start':<5} {'Ergebnis':<44} Gen")
    gefunden = 0
    for kopf in gene:
        position = kopf.split(":")[1].split(" ")[0]          # z. B. "87-1109" oder "c5888-4815"
        beschreibung = kopf.split(" ", 1)[1].split(" [")[0]
        strang = "-" if position.startswith("c") else "+"
        a, b = position.lstrip("c").split("-")
        start, ende = int(a), int(b)
        startcodon = gene[kopf][:3]

        treffer = None
        for orf in orfs:
            if orf["strang"] == strang and orf["ende"] == ende:
                treffer = orf
        if treffer is None:
            ergebnis = "❌ nicht gefunden"
        elif treffer["start"] == start:
            ergebnis = "✅ exakt gefunden"
            gefunden = gefunden + 1
        else:
            differenz = abs(treffer["start"] - start)
            richtung = "länger" if len(treffer["dna"]) > len(gene[kopf]) else "kürzer"
            ergebnis = f"🟡 gefunden, Start {differenz} bp anders ({richtung})"
            gefunden = gefunden + 1
        print(f"{position:<14} {startcodon:<5} {ergebnis:<44} {beschreibung}")
    print(f"\n{gefunden} von {len(gene)} Genen wiedergefunden. ORFs insgesamt: {len(orfs)}")


if EIGENE_SEQUENZ.strip():
    dna = "".join(EIGENE_SEQUENZ.split()).upper().replace("U", "T")
    print(f"Eigene Sequenz: {len(dna)} bp\n")
    orfs = finde_orfs(dna)
    zeige_orfs(orfs)
else:
    plasmid, gene = lade_pPCP1()
    print(f"Plasmid pPCP1 (NCBI NC_005816.1): {len(plasmid)} bp")
    print(f"Startcodons: {STARTCODONS}, Mindestlänge: {MIN_AMINOSAEUREN} Aminosäuren\n")
    orfs = finde_orfs(plasmid)
    zeige_orfs(orfs)
    vergleiche_mit_annotation(orfs, gene)
