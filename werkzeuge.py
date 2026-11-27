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
    "Middle East respiratory syndrome-related coronavirus": "MERS-CoV",
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
    text = _re.sub(r"^(PREDICTED|UNVERIFIED|TPA):\s*", "", text)
    text = text.split(",")[0]                                       # ", complete genome" usw. weg
    text = _re.split(r" (isolate|strain|str\.|clone|chromosome|genomic sequence|complete) ", " " + text + " ")[0]
    text = text.strip() or teile[1].strip()
    for lang, kurz in _ABKUERZUNGEN.items():
        text = text.replace(lang, kurz)
    # "Homo sapiens hemoglobin subunit beta (HBB)" → "Homo sapiens HBB"
    ohne_praefix = _re.sub(r"^(PREDICTED|UNVERIFIED|TPA):\s*", "", teile[1])
    gen = _re.match(r"^([A-Z][a-z]+ [a-z]+)\b[^()]*\(([A-Za-z0-9-]+)\),? *(mRNA|gene|transcript)?", ohne_praefix)
    if gen:
        text = f"{gen.group(1)} {gen.group(2)}"
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


def lade_von_ncbi(accessions, ordner):
    """Lädt Nukleotidsequenzen per NCBI-Accession als FASTA in den Ordner (eine Datei pro Nummer).
    Gibt Liste von (accession, fehlertext_oder_None) zurück. Braucht Internet (z. B. in Colab)."""
    import time
    import urllib.request
    os.makedirs(ordner, exist_ok=True)
    ergebnis = []
    for acc in accessions:
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
               f"&id={acc}&rettype=fasta&retmode=text")
        try:
            with urllib.request.urlopen(url, timeout=60) as antwort:
                text = antwort.read().decode()
            if not text.startswith(">"):
                raise ValueError("Antwort ist keine FASTA-Datei")
            with open(os.path.join(ordner, acc + ".fasta"), "w") as f:
                f.write(text)
            ergebnis.append((acc, None))
        except Exception as fehler:
            ergebnis.append((acc, str(fehler)))
        time.sleep(0.4)                       # NCBI erlaubt ohne Anmeldung ca. 3 Anfragen pro Sekunde
    return ergebnis


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


# ----------------------------------------------------------------------
# Schnelles Alignment mit affiner Lückenstrafe (Gotoh-Algorithmus, mit numpy)
# ----------------------------------------------------------------------
# Warum? Bei der einfachen (linearen) Strafe kostet jede Lückenposition gleich viel. Dann ist es
# für den Algorithmus egal, ob er EINE Lücke von 3 Basen macht oder drei verstreute Lücken von je
# 1 Base. Verstreute 1-Basen-Lücken sehen in einem Gen aber aus wie Leseraster-Verschiebungen
# (Frameshifts), die es in Wirklichkeit gar nicht gibt.
# Affine Strafe: Eine Lücke der Länge L kostet  OEFFNEN + L × VERLAENGERN.
# Das Öffnen ist teuer, das Verlängern billig, also entstehen wenige, zusammenhängende Lücken,
# so wie echte Insertionen/Deletionen meist aussehen.
# Standardwerte = Voreinstellung von NCBI BLASTN: Match +2, Mismatch −3, Öffnen 5, Verlängern 2.
#
# freie_enden=True: Lücken ganz am Anfang oder Ende kosten nichts ("semiglobal"). Das passt nur
# für SEHR ähnliche Sequenzen, von denen eine länger ist. Vorsicht: Bei entfernt verwandten
# Sequenzen gewinnt dann oft eine winzige Überlappung von wenigen Basen. Die Ordner-Analyse
# benutzt deshalb das normale globale Alignment (freie_enden=False).

def ausrichten(s1, s2, modus="global", match=2, mismatch=-3, luecke_oeffnen=5, luecke_verlaengern=2,
               freie_enden=False, matrix=None):
    """Alignment mit affiner Lückenstrafe. modus = "global" oder "lokal".
    matrix: optional eine Substitutionsmatrix (z. B. BLOSUM62 für Proteine) statt match/mismatch.
    Gibt dict mit score, a1, a2, start1, start2 (ab 1) zurück, wie needleman_wunsch."""
    import numpy as np
    tabelle = None
    if matrix is not None:                          # Punkte für jedes Buchstabenpaar nachschlagen
        tabelle = np.full((128, 128), min(min(z.values()) for z in matrix.values()), dtype=float)
        for x, zeile in matrix.items():
            for y, wert in zeile.items():
                tabelle[ord(x), ord(y)] = wert
    lokal = modus == "lokal"
    n, m = len(s1), len(s2)
    o, e = float(luecke_oeffnen), float(luecke_verlaengern)
    NEG = -1e18
    b2 = np.frombuffer(s2.encode(), dtype=np.uint8)
    k_e = np.arange(m + 1) * e                     # für die Lücken innerhalb einer Zeile

    # Zeiger pro Zelle (1 Byte): Bit 0-1 woher H kommt (0 Ende, 1 diagonal, 2 von oben, 3 von links),
    # Bit 2: Lücke von oben wird verlängert, Bit 3: Lücke von links wird verlängert
    zeiger = np.zeros((n + 1, m + 1), dtype=np.uint8)
    if lokal or freie_enden:
        H = np.zeros(m + 1)
    else:
        H = -(o + k_e); H[0] = 0.0
        zeiger[0, 2:] = 3 | 8
        if m:
            zeiger[0, 1] = 3
    F = np.full(m + 1, NEG)
    bester, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        # Lücke von oben (Base aus s1 gegen "-"): öffnen oder verlängern
        F_oeffnen = H - o - e
        F_verl = F - e
        F_ext = F_verl > F_oeffnen
        F = np.where(F_ext, F_verl, F_oeffnen)
        # diagonal (Match / Mismatch)
        diag = np.full(m + 1, NEG)
        if tabelle is None:
            diag[1:] = H[:-1] + np.where(b2 == ord(s1[i - 1]), match, mismatch)
        else:
            diag[1:] = H[:-1] + tabelle[ord(s1[i - 1]), b2]
        Ht = np.maximum(diag, F)
        quelle = np.where(diag >= F, 1, 2).astype(np.uint8)
        if lokal:
            null = Ht <= 0
            Ht = np.where(null, 0.0, Ht); quelle[null] = 0
        if lokal or freie_enden:
            Ht[0] = 0.0; quelle[0] = 0
        else:
            Ht[0] = -(o + i * e); quelle[0] = 2
        # Lücke von links (innerhalb der Zeile): E[j] = max über k<j von Ht[k] − o − (j−k)·e
        Eg = np.full(m + 1, NEG)
        Eg[1:] = np.maximum.accumulate(Ht[:-1] + k_e[:-1]) - o - k_e[1:]
        E_ext = np.zeros(m + 1, dtype=bool)
        E_ext[2:] = Eg[1:-1] - e > Ht[1:-1] - o - e
        links = Eg > Ht
        H = np.where(links, Eg, Ht)
        quelle[links] = 3
        zeiger[i] = quelle | (F_ext.astype(np.uint8) << 2) | (E_ext.astype(np.uint8) << 3)
        if lokal:
            j = int(np.argmax(H))
            if H[j] > bester:
                bester, bi, bj = float(H[j]), i, j
        elif freie_enden and m > 0 and (i == 1 or H[m] > bester):
            bester, bi, bj = float(H[m]), i, m      # bestes Ende in der letzten Spalte
    if lokal:
        i, j, score = bi, bj, bester
    elif freie_enden and n > 0 and m > 0:
        j_best = int(np.argmax(H))                  # bestes Ende in der letzten Zeile
        if H[j_best] >= bester:
            i, j, score = n, j_best, float(H[j_best])
        else:
            i, j, score = bi, bj, bester
    elif freie_enden:
        i, j, score = n, m, 0.0
    else:
        i, j, score = n, m, float(H[m])

    a1, a2 = [], []
    if freie_enden and not lokal:                  # überstehendes Ende ohne Strafe
        for k in range(n, i, -1):
            a1.append(s1[k - 1]); a2.append("-")
        for k in range(m, j, -1):
            a1.append("-"); a2.append(s2[k - 1])
    zustand = 0                                    # 0 = H, 2 = Lücke von oben, 3 = Lücke von links
    while i > 0 or j > 0:
        z = int(zeiger[i, j])
        if zustand == 0:
            q = z & 3
            if q == 0:                             # Anfang erreicht (lokal oder freie Enden)
                break
            if q == 1:
                a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
                continue
            zustand = q
        if zustand == 2:
            a1.append(s1[i - 1]); a2.append("-")
            zustand = 2 if (z >> 2) & 1 else 0
            i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1])
            zustand = 3 if (z >> 3) & 1 else 0
            j -= 1
    if freie_enden and not lokal:                  # Anfang ohne Strafe
        while i > 0:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        while j > 0:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
    a1 = "".join(reversed(a1)); a2 = "".join(reversed(a2))
    if score == int(score):
        score = int(score)
    return {"score": score, "a1": a1, "a2": a2, "start1": i + 1 if lokal else 1, "start2": j + 1 if lokal else 1}


def alignment_statistik(a1, a2):
    """Identität, Mismatches und Lücken im überlappenden Bereich.
    Lücken ganz am Anfang/Ende zählen nicht: Dort ist eine Sequenz nur länger (z. B. längere mRNA-Enden)."""
    spalten = [k for k in range(len(a1)) if a1[k] != "-" and a2[k] != "-"]
    if spalten:
        a1, a2 = a1[spalten[0]:spalten[-1] + 1], a2[spalten[0]:spalten[-1] + 1]
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


# ======================================================================
# Variantenanalyse (Teil 4): Mutationen gegenüber einer Referenz und ihre Wirkung
# ======================================================================

DREIBUCHSTABEN = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
                  "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
                  "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}


def finde_cds(dna, min_as=30):
    """Codierende Sequenz (CDS) automatisch finden: längster ORF auf dem + Strang (ATG bis Stopp).
    Gibt (start, ende) biologisch (ab 1, inkl. Stoppcodon) zurück oder None."""
    orfs = [o for o in finde_orfs(dna, ("ATG",), min_as=min_as) if o["strang"] == "+"]
    if not orfs:
        return None
    bester = max(orfs, key=lambda o: len(o["dna"]))
    return bester["start"], bester["ende"]


def finde_varianten(a1, a2):
    """Liest aus einem Alignment (a1 = Referenz, a2 = Variante) alle Unterschiede.
    Lücken am Anfang/Ende zählen nicht (die Sequenzen sind dort nur unterschiedlich lang).
    Positionen beziehen sich auf die Referenz (ab 1)."""
    spalten = [k for k in range(len(a1)) if a1[k] != "-" and a2[k] != "-"]
    if not spalten:
        return []
    erste, letzte = spalten[0], spalten[-1]
    varianten = []
    rpos = 0
    k = 0
    while k < len(a1):
        r, v = a1[k], a2[k]
        if r != "-":
            rpos += 1
        if k < erste or k > letzte:
            k += 1
            continue
        if r != "-" and v != "-":
            if r != v and "N" not in (r, v):
                varianten.append({"art": "SNV", "start": rpos, "ende": rpos, "ref": r, "alt": v})
            k += 1
        elif v == "-":                                   # Deletion: Basen fehlen in der Variante
            start, geloescht = rpos, r
            while k + 1 < len(a1) and a2[k + 1] == "-" and a1[k + 1] != "-":
                k += 1
                rpos += 1
                geloescht += a1[k]
            varianten.append({"art": "Deletion", "start": start, "ende": rpos, "ref": geloescht, "alt": ""})
            k += 1
        else:                                            # Insertion: zusätzliche Basen in der Variante
            eingefuegt = v
            while k + 1 < len(a1) and a1[k + 1] == "-":
                k += 1
                eingefuegt += a2[k]
            varianten.append({"art": "Insertion", "start": rpos, "ende": rpos + 1, "ref": "", "alt": eingefuegt})
            k += 1
    return varianten


def normalisiere_indels(referenz, varianten):
    """HGVS-Regel: Liegt eine Insertion/Deletion in einer Wiederholung (z. B. TTCTTT), ist ihre Position
    mehrdeutig. Man gibt dann die am weitesten 3' (rechts) liegende Position an."""
    for v in varianten:
        if v["art"] == "Deletion":
            a, b = v["start"], v["ende"]                     # ab 1, inklusive
            while b < len(referenz) and referenz[b] == referenz[a - 1]:
                a, b = a + 1, b + 1
            v["start"], v["ende"], v["ref"] = a, b, referenz[a - 1:b]
        elif v["art"] == "Insertion":
            nach, ins = v["start"], v["alt"]                 # Insertion nach Position 'nach'
            while nach < len(referenz) and referenz[nach] == ins[0]:
                ins = ins[1:] + ins[0]
                nach += 1
            v["start"], v["ende"], v["alt"] = nach, nach + 1, ins
    return varianten


def _c_position(p, cds):
    """Position in HGVS-c.-Schreibweise: c.1 = A des Startcodons, c.-5 = 5 Basen davor, c.*3 = 3 Basen nach dem Stopp."""
    start, ende = cds
    if p < start:
        return f"-{start - p}"
    if p > ende:
        return f"*{p - ende}"
    return str(p - start + 1)


def wirkung_der_varianten(referenz, varianten, cds):
    """Ergänzt jede Variante um c.-Notation, p.-Notation und Wirkung (stumm, Missense, Nonsense, Frameshift ...)."""
    start, ende = cds
    snv_im_codon = {}
    for v in varianten:
        if v["art"] == "SNV" and start <= v["start"] <= ende:
            snv_im_codon.setdefault((v["start"] - start) // 3, []).append(v)

    for v in varianten:
        if v["art"] == "SNV":
            v["c"] = f"c.{_c_position(v['start'], cds)}{v['ref']}>{v['alt']}"
            if not (start <= v["start"] <= ende):
                v["p"], v["wirkung"] = "", "außerhalb des Gens (nicht codierend)"
                continue
            nr = (v["start"] - start) // 3
            codon_start = start - 1 + nr * 3
            ref_codon = referenz[codon_start:codon_start + 3]
            alt_codon = list(ref_codon)
            for andere in snv_im_codon[nr]:              # alle SNVs im selben Codon gemeinsam anwenden
                alt_codon[andere["start"] - 1 - codon_start] = andere["alt"]
            alt_codon = "".join(alt_codon)
            as_ref, as_alt = CODON_TABELLE.get(ref_codon, "X"), CODON_TABELLE.get(alt_codon, "X")
            v["codon"] = f"{ref_codon}>{alt_codon}"
            name = f"{DREIBUCHSTABEN[as_ref]}{nr + 1}"
            if as_ref == as_alt:
                v["p"], v["wirkung"] = f"p.{name}=", "stumm (gleiche Aminosäure)"
            elif as_alt == "*":
                v["p"], v["wirkung"] = f"p.{name}*", "Nonsense (vorzeitiges Stoppcodon)"
            elif as_ref == "*":
                v["p"], v["wirkung"] = f"p.{name}{DREIBUCHSTABEN[as_alt]}ext", "Stoppverlust (Protein wird verlängert)"
            elif nr == 0:
                v["p"], v["wirkung"] = "p.Met1?", "Startcodon verändert (Protein entsteht evtl. nicht)"
            else:
                v["p"], v["wirkung"] = f"p.{name}{DREIBUCHSTABEN[as_alt]}", "Missense (andere Aminosäure)"
        else:
            laenge = len(v["ref"]) + len(v["alt"])
            if v["art"] == "Deletion":
                a, b = _c_position(v["start"], cds), _c_position(v["ende"], cds)
                v["c"] = f"c.{a}del" if v["start"] == v["ende"] else f"c.{a}_{b}del"
                im_gen = v["start"] <= ende and v["ende"] >= start
            else:
                v["c"] = f"c.{_c_position(v['start'], cds)}_{_c_position(v['start'] + 1, cds)}ins{v['alt']}"
                im_gen = start <= v["start"] < ende
            v["p"] = ""
            if not im_gen:
                v["wirkung"] = "außerhalb des Gens (nicht codierend)"
            elif laenge % 3 == 0:
                art = "fehlen" if v["art"] == "Deletion" else "kommen hinzu"
                v["wirkung"] = f"In-frame ({laenge // 3} Aminosäure(n) {art}, Leseraster bleibt)"
            else:
                v["wirkung"] = "Frameshift (Leserasterverschiebung, Protein ab hier verändert)"

    # Zwei Lücken kurz hintereinander, die sich gegenseitig ausgleichen (z. B. −2 bp und +2 bp):
    # Das Leseraster ist nur auf wenigen Basen verschoben. Bei verwandten Arten ist das meist
    # ein Alignment-Artefakt (statt einiger Mismatches hat das Alignment zwei Lücken gesetzt).
    indels = sorted((v for v in varianten if v["art"] != "SNV" and v["wirkung"].startswith("Frameshift")),
                    key=lambda v: v["start"])
    k = 0
    while k < len(indels):
        verschiebung, gruppe = 0, []
        for w in indels[k:]:
            if gruppe and w["start"] - gruppe[0]["start"] > 30:
                break
            verschiebung += len(w["alt"]) - len(w["ref"])
            gruppe.append(w)
            if verschiebung % 3 == 0:
                break
        if len(gruppe) > 1 and verschiebung % 3 == 0:
            abstand = gruppe[-1]["start"] - gruppe[0]["start"]
            for w in gruppe:
                w["wirkung"] = (f"In-frame insgesamt (Frameshift nach {abstand} bp wieder ausgeglichen; "
                                f"oft ein Alignment-Artefakt statt echter Lücken)")
            k += len(gruppe)
        else:
            k += 1
    return varianten


def variante_auf_referenz_position(a1, a2, ref_position):
    """Welche Position (ab 1) in der Variante entspricht der Referenzposition? None, wenn dort eine Lücke ist."""
    rpos = vpos = 0
    for r, v in zip(a1, a2):
        if r != "-":
            rpos += 1
        if v != "-":
            vpos += 1
        if rpos == ref_position and r != "-":
            return vpos if v != "-" else None
    return None


def vergleiche_proteine(p_ref, p_var):
    if p_ref == p_var:
        return "identisch"
    if len(p_ref) == len(p_var):
        unterschiede = [i for i in range(len(p_ref)) if p_ref[i] != p_var[i]]
        return f"{len(unterschiede)} Aminosäure(n) verschieden, gleiche Länge ({len(p_ref)} AS)"
    i = 0
    while i < min(len(p_ref), len(p_var)) and p_ref[i] == p_var[i]:
        i += 1
    return f"Länge {len(p_ref)} → {len(p_var)} AS, erste Abweichung an Aminosäure {i + 1}"


# ======================================================================
# Proteine vergleichen
# ======================================================================
# BLOSUM62: Punkte für jedes Aminosäurepaar, abgeleitet aus echten, verwandten Proteinen
# (Henikoff & Henikoff 1992). Positiv = der Austausch kommt in der Natur oft vor (ähnliche
# Aminosäuren, z. B. I↔V +3), negativ = selten (z. B. W↔G −2). Standard bei BLASTP.
# B = D oder N, Z = E oder Q, X = unbekannt, * = Stopp.
_BLOSUM62_TEXT = """
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  A   4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
  R  -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
  N  -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
  D  -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
  C   0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
  Q  -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
  E  -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
  G   0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
  H  -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
  I  -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
  L  -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
  K  -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
  M  -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
  F  -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
  P  -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
  S   1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
  T   0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
  W  -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
  Y  -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
  V   0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
  B  -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
  Z  -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
  X   0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
  *  -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _lese_matrix(text):
    zeilen = [z.split() for z in text.strip().split("\n")]
    kopf = zeilen[0]
    return {z[0]: {b: int(w) for b, w in zip(kopf, z[1:])} for z in zeilen[1:]}


BLOSUM62 = _lese_matrix(_BLOSUM62_TEXT)

# Vereinfachte Lehrbuch-Einteilung der Aminosäuren nach ihrer Seitenkette
EIGENSCHAFT = {}
for _gruppe, _buchstaben in [("unpolar", "GAVLIMFWP"), ("polar", "STCYNQ"), ("basisch (+)", "KRH"),
                             ("sauer (−)", "DE")]:
    for _b in _buchstaben:
        EIGENSCHAFT[_b] = _gruppe


def ist_protein(seq):
    """True, wenn die Sequenz eher Aminosäuren als Basen enthält."""
    if not seq:
        return False
    basen = sum(seq.count(b) for b in "ACGTUN")
    return basen < 0.9 * len(seq)


def protein_aus_dna(dna, min_as=30):
    """Protein des Gens in einer DNA/mRNA: längster ORF ab ATG auf dem + Strang.
    Gibt (protein, (start, ende)) oder (None, None) zurück."""
    cds = finde_cds(dna, min_as=min_as)
    if not cds:
        return None, None
    return uebersetze(dna[cds[0] - 1:cds[1]]), cds


def ausrichten_protein(p1, p2, modus="global", luecke_oeffnen=11, luecke_verlaengern=1):
    """Protein-Alignment mit BLOSUM62 und affinen Lücken (Standardwerte von NCBI BLASTP)."""
    p1, p2 = p1.upper().rstrip("*"), p2.upper().rstrip("*")
    return ausrichten(p1, p2, modus, luecke_oeffnen=luecke_oeffnen, luecke_verlaengern=luecke_verlaengern,
                      matrix=BLOSUM62)


def protein_statistik(a1, a2):
    """Identität (gleiche Aminosäure) und Ähnlichkeit (BLOSUM62 > 0, also ähnliche Aminosäure)
    im überlappenden Bereich."""
    spalten = [k for k in range(len(a1)) if a1[k] != "-" and a2[k] != "-"]
    if not spalten:
        return {"laenge": 0, "gleich": 0, "aehnlich": 0, "luecken": 0, "identitaet": 0.0, "aehnlichkeit": 0.0}
    a1, a2 = a1[spalten[0]:spalten[-1] + 1], a2[spalten[0]:spalten[-1] + 1]
    gleich = aehnlich = luecken = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            luecken += 1
        elif x == y:
            gleich += 1; aehnlich += 1
        elif BLOSUM62.get(x, {}).get(y, -4) > 0:
            aehnlich += 1
    n = len(a1)
    return {"laenge": n, "gleich": gleich, "aehnlich": aehnlich, "luecken": luecken,
            "identitaet": gleich / n * 100, "aehnlichkeit": aehnlich / n * 100}


def protein_unterschiede(a1, a2):
    """Alle Unterschiede im Alignment, Positionen in Protein 1 (ab 1), mit Eigenschaftsänderung."""
    spalten = [k for k in range(len(a1)) if a1[k] != "-" and a2[k] != "-"]
    if not spalten:
        return []
    erste, letzte = spalten[0], spalten[-1]
    liste, pos = [], 0
    for k, (x, y) in enumerate(zip(a1, a2)):
        if x != "-":
            pos += 1
        if k < erste or k > letzte or x == y:
            continue
        if x == "-":
            liste.append({"position": pos, "art": "Insertion", "ref": "", "alt": y, "punkte": None,
                          "bewertung": "zusätzliche Aminosäure"})
        elif y == "-":
            liste.append({"position": pos, "art": "Deletion", "ref": x, "alt": "", "punkte": None,
                          "bewertung": "Aminosäure fehlt"})
        else:
            punkte = BLOSUM62.get(x, {}).get(y, -4)
            e1, e2 = EIGENSCHAFT.get(x, "?"), EIGENSCHAFT.get(y, "?")
            if e1 == e2:
                text = f"ähnlich: beide {e1}"
            else:
                text = f"{e1} → {e2}"
            if punkte > 0:
                text += ", häufiger Austausch (konservativ)"
            elif punkte < 0:
                text += ", seltener Austausch"
            liste.append({"position": pos, "art": "Austausch", "ref": x, "alt": y, "punkte": punkte,
                          "bewertung": text})
    # zusammenhängende Lücken zusammenfassen
    zusammen = []
    for d in liste:
        letzte_d = zusammen[-1] if zusammen else None
        if (letzte_d and d["art"] == letzte_d["art"] and d["art"] != "Austausch" and
                d["position"] - letzte_d["ende"] <= (1 if d["art"] == "Deletion" else 0)):
            letzte_d["ref"] += d["ref"]; letzte_d["alt"] += d["alt"]; letzte_d["ende"] = d["position"]
        else:
            zusammen.append(dict(d, ende=d["position"]))
    return zusammen


def protein_alignment_text(erg, breite=60):
    """Protein-Alignment als Text; Mittelzeile: | gleich, + ähnlich (BLOSUM62 > 0), Leerzeichen sonst."""
    a1, a2 = erg["a1"], erg["a2"]
    mitte = "".join("|" if x == y else ("+" if x != "-" and y != "-" and BLOSUM62.get(x, {}).get(y, -4) > 0
                                        else " ") for x, y in zip(a1, a2))
    zeilen, p1, p2 = [], erg["start1"], erg["start2"]
    for k in range(0, len(a1), breite):
        t1, t2 = a1[k:k + breite], a2[k:k + breite]
        zeilen += [f"Prot1 {p1:>6}  {t1}", f"              {mitte[k:k + breite]}", f"Prot2 {p2:>6}  {t2}", ""]
        p1 += len(t1) - t1.count("-")
        p2 += len(t2) - t2.count("-")
    return "\n".join(zeilen)
