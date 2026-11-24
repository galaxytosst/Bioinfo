"""
🧰 werkzeuge.py – alle Bioinformatik-Funktionen aus diesem Projekt an einem Ort.

Dieses Modul führt beim Laden nichts aus, es stellt nur Funktionen bereit.
Benutzt wird es von analyse.py (Ordner-Analyse) und kann auch direkt importiert werden:

    from werkzeuge import lese_fasta, gc_gehalt, finde_orfs, needleman_wunsch
"""

import math
import os

# ======================================================================
# Grundlagen (Teil 0)
# ======================================================================

FASTA_ENDUNGEN = (".fasta", ".fa", ".fna", ".fas", ".ffn", ".faa", ".txt")

_basen = "TCAG"
_aminosaeuren = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABELLE = {}
_nr = 0
for _b1 in _basen:
    for _b2 in _basen:
        for _b3 in _basen:
            CODON_TABELLE[_b1 + _b2 + _b3] = _aminosaeuren[_nr]   # DNA-Codons (T statt U)
            _nr += 1

KOMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def saeubern(sequenz):
    """Großbuchstaben, keine Leerzeichen, U → T, alle anderen Zeichen (IUPAC-Codes wie R, Y, ...) → N."""
    sequenz = "".join(sequenz.split()).upper().replace("U", "T")
    return "".join(b if b in "ACGT" else "N" for b in sequenz)


def kurzer_name(kopfzeile):
    """Erstes Wort der Kopfzeile. NCBI-Kopfzeilen wie 'gi|45478711|ref|NC_005816.1|' werden zu 'NC_005816.1'."""
    wort = kopfzeile.split()[0] if kopfzeile.split() else "ohne_namen"
    if "|" not in wort:
        return wort
    teile = [t for t in wort.split("|") if t]
    datenbanken = {"gi", "ref", "gb", "emb", "dbj", "sp", "tr", "lcl", "pdb"}
    rest = []
    for nr, t in enumerate(teile):
        if t in datenbanken:
            continue
        if nr > 0 and teile[nr - 1] == "gi":       # GI-Nummer überspringen
            continue
        rest.append(t)
    return "".join(rest) if rest else wort


import re as _re
_ACCESSION = _re.compile(r"^[A-Z]{1,2}_?[0-9]{5,9}(\.[0-9]+)?(:\S+)?$")


_ABKUERZUNGEN = {
    "Severe acute respiratory syndrome coronavirus 2": "SARS-CoV-2",
    "Severe acute respiratory syndrome-related coronavirus": "SARS-CoV",
    "Severe acute respiratory syndrome coronavirus": "SARS-CoV",
    "Middle East respiratory syndrome coronavirus": "MERS-CoV",
    "Human immunodeficiency virus 1": "HIV-1",
}


def beschriftung(kopfzeile, max_zeichen=45):
    """Lesbarer Name aus der Kopfzeile.
    NCBI-Kopfzeile '>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    → 'Severe acute respiratory syndrome… (NC_045512.2)'. Andere Kopfzeilen → erstes Wort."""
    name = kurzer_name(kopfzeile)
    teile = kopfzeile.split(maxsplit=1)
    if not _ACCESSION.match(name) or len(teile) < 2:
        return name
    text = _re.sub(r"\[.*?\]", "", teile[1])                     # [Organismus]-Angaben am Ende weg
    text = text.split(",")[0]                                       # ", complete genome" usw. weg
    text = _re.split(r" (isolate|strain|str\.|clone|chromosome|genomic sequence|complete) ", " " + text + " ")[0]
    text = text.strip() or teile[1].strip()
    for lang, kurz in _ABKUERZUNGEN.items():
        text = text.replace(lang, kurz)
    if len(text) > max_zeichen:
        text = text[:max_zeichen - 1].rstrip() + "…"
    return f"{text} ({name})"


def lese_fasta(dateiname):
    """Liest eine FASTA-Datei. Gibt Dictionary {Name: Sequenz} zurück (Name = erstes Wort der Kopfzeile)."""
    sequenzen = {}
    name = None
    with open(dateiname) as f:
        for zeile in f:
            zeile = zeile.strip()
            if zeile == "":
                continue
            if zeile.startswith(">"):
                name = beschriftung(zeile[1:]) if zeile[1:].strip() else f"ohne_namen_{len(sequenzen) + 1}"
                sequenzen[name] = []
            elif name is not None:
                sequenzen[name].append(zeile)        # Liste statt String: viel schneller bei großen Dateien
    return {n: saeubern("".join(teile)) for n, teile in sequenzen.items()}


def lese_ordner(ordner):
    """Liest alle FASTA-Dateien eines Ordners. Gibt {Name: Sequenz} und {Name: Dateiname} zurück."""
    alle = {}
    herkunft = {}
    for datei in sorted(os.listdir(ordner)):
        if not datei.lower().endswith(FASTA_ENDUNGEN):
            continue
        for name, seq in lese_fasta(os.path.join(ordner, datei)).items():
            eindeutig = name
            zaehler = 2
            while eindeutig in alle:                   # gleiche Namen in verschiedenen Dateien
                eindeutig = f"{name}_{zaehler}"
                zaehler += 1
            alle[eindeutig] = seq
            herkunft[eindeutig] = datei
    return alle, herkunft


def transkribiere(dna):
    return dna.replace("T", "U")


def reverse_complement(dna):
    return "".join(KOMPLEMENT.get(b, "N") for b in reversed(dna))


def gc_gehalt(dna):
    """GC-Anteil in Prozent, bezogen auf alle eindeutigen Basen (N wird nicht mitgezählt)."""
    bekannt = len(dna) - dna.count("N")
    if bekannt == 0:
        return 0.0
    return (dna.count("G") + dna.count("C")) / bekannt * 100


def uebersetze(dna, bis_stopp=True):
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aminosaeure = CODON_TABELLE.get(dna[i:i + 3], "X")
        if aminosaeure == "*" and bis_stopp:
            break
        protein.append(aminosaeure)
    return "".join(protein)


def zusammensetzung(dna):
    return {b: dna.count(b) for b in "ACGTN"}


# ======================================================================
# Alignment (Teil 1), speichersparend
# ======================================================================
# Statt einer Tabelle voller Zahlen wird nur die aktuelle und die vorige Zeile gespeichert.
# Für den Traceback merkt sich jede Zelle nur, woher sie kam (1 Byte):
#   0 = Ende (nur lokal), 1 = diagonal, 2 = oben (Lücke in s2), 3 = links (Lücke in s1)
# Speicher: ca. Länge1 × Länge2 Bytes, z. B. 5.000 × 5.000 → 25 MB.

def _traceback(s1, s2, richtung, i, j, lokal):
    a1, a2 = [], []
    while i > 0 or j > 0:
        r = richtung[i][j]
        if r == 1:
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif r == 2:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        elif r == 3:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
        else:
            break
    return "".join(reversed(a1)), "".join(reversed(a2)), i, j


def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    """Globales Alignment. Gibt dict mit score, a1, a2, start1, start2 zurück."""
    n, m = len(s1), len(s2)
    richtung = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        richtung[0][j] = 3
    vorher = [j * gap for j in range(m + 1)]
    for i in range(1, n + 1):
        aktuell = [i * gap] + [0] * m
        richtung[i][0] = 2
        r_zeile = richtung[i]
        a = s1[i - 1]
        for j in range(1, m + 1):
            diag = vorher[j - 1] + (match if a == s2[j - 1] else mismatch)
            oben = vorher[j] + gap
            links = aktuell[j - 1] + gap
            if diag >= oben and diag >= links:
                aktuell[j] = diag; r_zeile[j] = 1
            elif oben >= links:
                aktuell[j] = oben; r_zeile[j] = 2
            else:
                aktuell[j] = links; r_zeile[j] = 3
        vorher = aktuell
    a1, a2, _, _ = _traceback(s1, s2, richtung, n, m, False)
    return {"score": vorher[m], "a1": a1, "a2": a2, "start1": 1, "start2": 1}


def smith_waterman(s1, s2, match=1, mismatch=-1, gap=-2):
    """Lokales Alignment. Gibt dict mit score, a1, a2, start1, start2 (ab 1) zurück."""
    n, m = len(s1), len(s2)
    richtung = [bytearray(m + 1) for _ in range(n + 1)]
    vorher = [0] * (m + 1)
    bester, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        aktuell = [0] * (m + 1)
        r_zeile = richtung[i]
        a = s1[i - 1]
        for j in range(1, m + 1):
            diag = vorher[j - 1] + (match if a == s2[j - 1] else mismatch)
            oben = vorher[j] + gap
            links = aktuell[j - 1] + gap
            wert = max(0, diag, oben, links)
            if wert == 0:
                continue                           # aktuell[j] = 0, Richtung 0 = Ende
            aktuell[j] = wert
            if wert == diag:
                r_zeile[j] = 1
            elif wert == oben:
                r_zeile[j] = 2
            else:
                r_zeile[j] = 3
            if wert > bester:
                bester, bi, bj = wert, i, j
        vorher = aktuell
    a1, a2, i, j = _traceback(s1, s2, richtung, bi, bj, True)
    return {"score": bester, "a1": a1, "a2": a2, "start1": i + 1, "start2": j + 1}


def alignment_statistik(a1, a2):
    matches = mismatches = luecken = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            luecken += 1
        elif x == y:
            matches += 1
        else:
            mismatches += 1
    laenge = len(a1)
    return {"laenge": laenge, "matches": matches, "mismatches": mismatches, "luecken": luecken,
            "identitaet": matches / laenge * 100 if laenge else 0.0}


def alignment_als_text(ergebnis, breite=60):
    a1, a2 = ergebnis["a1"], ergebnis["a2"]
    mitte = "".join("|" if x == y else (" " if "-" in (x, y) else ".") for x, y in zip(a1, a2))
    zeilen = []
    p1, p2 = ergebnis["start1"], ergebnis["start2"]
    for k in range(0, len(a1), breite):
        t1, t2 = a1[k:k + breite], a2[k:k + breite]
        zeilen += [f"Seq1 {p1:>8}  {t1}", f"               {mitte[k:k + breite]}", f"Seq2 {p2:>8}  {t2}", ""]
        p1 += len(t1) - t1.count("-")
        p2 += len(t2) - t2.count("-")
    return "\n".join(zeilen)


# ======================================================================
# Strang-Vergleich / Mutationen (aus straenge_vergleichen.py)
# ======================================================================

def mutationstyp(b1, b2):
    if "N" in (b1, b2):
        return "unklar"
    return "Transition" if (b1 in "AG") == (b2 in "AG") else "Transversion"


def mutationen(s1, s2):
    """Für zwei gleich lange Sequenzen: Liste aller Unterschiede mit Wirkung im Leseraster 1."""
    ergebnis = []
    for i in range(min(len(s1), len(s2))):
        if s1[i] == s2[i]:
            continue
        c = (i // 3) * 3
        k1, k2 = s1[c:c + 3], s2[c:c + 3]
        if len(k1) < 3:
            wirkung, as1, as2 = "kein vollständiges Codon", "", ""
        else:
            as1, as2 = CODON_TABELLE.get(k1, "X"), CODON_TABELLE.get(k2, "X")
            if as1 == as2:
                wirkung = "stumm"
            elif as2 == "*":
                wirkung = "Nonsense"
            elif as1 == "*":
                wirkung = "Stopp verloren"
            else:
                wirkung = "Missense"
        ergebnis.append({"position": i + 1, "von": s1[i], "zu": s2[i], "typ": mutationstyp(s1[i], s2[i]),
                         "codon": f"{k1}->{k2}", "aminosaeure": f"{as1}->{as2}", "wirkung": wirkung})
    return ergebnis


# ======================================================================
# ORF-Finder (Teil 2)
# ======================================================================

def _orfs_in_strang(seq, startcodons, stoppcodons, min_as):
    orfs = []
    for rahmen in range(3):
        start = None
        for i in range(rahmen, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if start is None:
                if codon in startcodons:
                    start = i
            elif codon in stoppcodons:
                if (i + 3 - start) // 3 - 1 >= min_as:
                    orfs.append((rahmen, start, i + 3))
                start = None
    return orfs


def finde_orfs(dna, startcodons=("ATG",), stoppcodons=("TAA", "TAG", "TGA"), min_as=100):
    """Alle ORFs in 6 Leserastern. Positionen biologisch (ab 1, Ende eingeschlossen).
    Auf dem Gegenstrang ist start > ende (wie bei NCBI 'complement')."""
    startcodons, stoppcodons = set(startcodons), set(stoppcodons)
    L = len(dna)
    orfs = []
    for rahmen, s, e in _orfs_in_strang(dna, startcodons, stoppcodons, min_as):
        orfs.append({"strang": "+", "rahmen": rahmen + 1, "start": s + 1, "ende": e, "dna": dna[s:e]})
    gegen = reverse_complement(dna)
    for rahmen, s, e in _orfs_in_strang(gegen, startcodons, stoppcodons, min_as):
        orfs.append({"strang": "-", "rahmen": rahmen + 1, "start": L - s, "ende": L - e + 1, "dna": gegen[s:e]})
    for orf in orfs:
        orf["protein"] = uebersetze(orf["dna"])
    orfs.sort(key=lambda o: min(o["start"], o["ende"]))
    return orfs


# ======================================================================
# Distanzen und Stammbaum (Teil 3)
# ======================================================================

def distanz_aus_alignment(a1, a2, korrektur="jukes-cantor"):
    vergleichbar = verschieden = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-" or x == "N" or y == "N":
            continue
        vergleichbar += 1
        if x != y:
            verschieden += 1
    if vergleichbar == 0:
        return float("inf")
    p = verschieden / vergleichbar
    if korrektur == "p":
        return p
    if p >= 0.75:
        return float("inf")
    return -0.75 * math.log(1 - 4 / 3 * p)


def kmer_menge(dna, k):
    return {dna[i:i + k] for i in range(len(dna) - k + 1) if "N" not in dna[i:i + k]}


def kmer_distanz(menge1, menge2, k):
    """Alignment-freie Distanz (Mash-Formel, Ondov et al. 2016): aus dem Anteil gemeinsamer k-mere
    (Jaccard-Index J) wird die Mutationsrate geschätzt: d = -1/k · ln(2J / (1 + J)).
    Schnell auch für lange Sequenzen, aber nur eine Näherung und unzuverlässig bei großen Distanzen (> ~0,25)."""
    gemeinsam = len(menge1 & menge2)
    vereinigung = len(menge1 | menge2)
    if vereinigung == 0 or gemeinsam == 0:
        return float("inf")
    j = gemeinsam / vereinigung
    return -1 / k * math.log(2 * j / (1 + j))


def _blatt(name):
    return {"name": name, "kinder": []}


def upgma(namen, D):
    gruppen = {n: {"knoten": _blatt(n), "groesse": 1, "hoehe": 0.0} for n in namen}
    D = {a: dict(D[a]) for a in D}
    zaehler = 0
    while len(gruppen) > 1:
        liste = list(gruppen)
        a, b = min(((x, y) for i, x in enumerate(liste) for y in liste[i + 1:]), key=lambda p: D[p[0]][p[1]])
        hoehe = D[a][b] / 2
        A, B = gruppen.pop(a), gruppen.pop(b)
        zaehler += 1
        neu = f"#U{zaehler}"
        knoten = {"name": None, "kinder": [(A["knoten"], hoehe - A["hoehe"]), (B["knoten"], hoehe - B["hoehe"])]}
        D[neu] = {neu: 0.0}
        for k in gruppen:
            d = (D[a][k] * A["groesse"] + D[b][k] * B["groesse"]) / (A["groesse"] + B["groesse"])
            D[neu][k] = D[k][neu] = d
        gruppen[neu] = {"knoten": knoten, "groesse": A["groesse"] + B["groesse"], "hoehe": hoehe}
    return list(gruppen.values())[0]["knoten"]


def neighbor_joining(namen, D, aussengruppe=None):
    """NJ-Baum. Mit Außengruppe (Liste von Namen) wird dort bewurzelt, sonst am ersten Namen."""
    D = {a: dict(D[a]) for a in D}
    aktiv = list(namen)
    nachbarn = {n: [] for n in namen}
    zaehler = 0
    if len(aktiv) == 2:
        a, b = aktiv
        return {"name": None, "kinder": [(_blatt(a), D[a][b] / 2), (_blatt(b), D[a][b] / 2)]}
    while len(aktiv) > 2:
        n = len(aktiv)
        r = {i: sum(D[i][k] for k in aktiv) for i in aktiv}
        i, j = min(((x, y) for a, x in enumerate(aktiv) for y in aktiv[a + 1:]),
                   key=lambda p: (n - 2) * D[p[0]][p[1]] - r[p[0]] - r[p[1]])
        li = D[i][j] / 2 + (r[i] - r[j]) / (2 * (n - 2))
        lj = D[i][j] - li
        zaehler += 1
        u = f"#N{zaehler}"
        nachbarn[u] = [(i, li), (j, lj)]
        nachbarn[i].append((u, li))
        nachbarn[j].append((u, lj))
        D[u] = {u: 0.0}
        for k in aktiv:
            if k not in (i, j):
                D[u][k] = D[k][u] = (D[i][k] + D[j][k] - D[i][j]) / 2
        aktiv = [k for k in aktiv if k not in (i, j)] + [u]
    a, b = aktiv
    nachbarn[a].append((b, D[a][b]))
    nachbarn[b].append((a, D[a][b]))

    def teilbaum(knoten, herkunft):
        return {"name": knoten if knoten in namen else None,
                "kinder": [(teilbaum(nb, knoten), l) for nb, l in nachbarn[knoten] if nb != herkunft]}

    def blaetter_hinter(knoten, herkunft):
        if knoten in namen:
            return {knoten}
        menge = set()
        for nb, _ in nachbarn[knoten]:
            if nb != herkunft:
                menge |= blaetter_hinter(nb, knoten)
        return menge

    ziel = set(aussengruppe) if aussengruppe else {namen[0]}
    for u in nachbarn:
        for v, laenge in nachbarn[u]:
            if blaetter_hinter(v, u) == ziel:
                return {"name": None, "kinder": [(teilbaum(v, u), laenge / 2), (teilbaum(u, v), laenge / 2)]}
    erstes = sorted(ziel)[0]
    eltern, laenge = nachbarn[erstes][0]
    return {"name": None, "kinder": [(_blatt(erstes), laenge / 2), (teilbaum(eltern, erstes), laenge / 2)]}


def newick_name(name):
    """Newick verbietet Leerzeichen und Zeichen wie , ( ) : ;"""
    return _re.sub(r"[\s,;:()\[\]']+", "_", name).strip("_")


def newick(knoten, laenge=None):
    text = "(" + ",".join(newick(k, l) for k, l in knoten["kinder"]) + ")" if knoten["kinder"] else newick_name(knoten["name"])
    return text + (f":{laenge:.5f}" if laenge is not None else "")


def baum_als_text(knoten):
    zeilen = []

    def schreibe(k, praefix, laenge, letzter, wurzel):
        name = k["name"] if k["name"] else "┐"
        zeilen.append(name if wurzel else praefix + ("└── " if letzter else "├── ") + name + f"  ({laenge:.4f})")
        neu = praefix + ("" if wurzel else ("    " if letzter else "│   "))
        for nr, (kind, l) in enumerate(k["kinder"]):
            schreibe(kind, neu, l, nr == len(k["kinder"]) - 1, False)

    schreibe(knoten, "", None, True, True)
    return "\n".join(zeilen)


def blaetter(knoten):
    if not knoten["kinder"]:
        return [knoten["name"]]
    return [b for kind, _ in knoten["kinder"] for b in blaetter(kind)]


def zeichne_baum(knoten, titel, ax):
    y_pos = {n: i for i, n in enumerate(blaetter(knoten))}

    def zeichne(k, x):
        if not k["kinder"]:
            ax.text(x, y_pos[k["name"]], "  " + (k["name"] if "(" in k["name"] else k["name"].replace("_", " ")), va="center", fontsize=9)
            return y_pos[k["name"]]
        ys = []
        for kind, l in k["kinder"]:
            y = zeichne(kind, x + l)
            ax.plot([x, x + l], [y, y], color="black", lw=1.2)
            ys.append(y)
        ax.plot([x, x], [min(ys), max(ys)], color="black", lw=1.2)
        return sum(ys) / len(ys)

    zeichne(knoten, 0.0)
    ax.set_title(titel)
    ax.set_yticks([])
    ax.set_xlabel("Distanz (Substitutionen pro Position)")
    for seite in ("left", "right", "top"):
        ax.spines[seite].set_visible(False)
    ax.set_xlim(0, ax.get_xlim()[1] * 1.6)
