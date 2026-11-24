# 🧬 Paarweises Sequenz-Alignment: Needleman-Wunsch (global) und Smith-Waterman (lokal)
#
# So benutzt du es: Den ganzen Code in eine Colab-Zelle kopieren, unten bei
# "MEINE EINGABE" zwei Sequenzen eintragen und die Zelle ausführen.
#
# Bewertung (Scoring), lineare Lückenstrafe:
#   gleiche Zeichen (Match)        → +MATCH
#   verschiedene Zeichen (Mismatch)→ +MISMATCH (negativ)
#   Lücke (Gap, Insertion/Deletion)→ +GAP      (negativ, pro Lückenposition)
#
# Funktioniert für DNA, RNA und Proteine (Proteine aber nur mit einfacher
# Match/Mismatch-Bewertung, nicht mit BLOSUM-Matrix).


# ============================================================
# MEINE EINGABE
# ============================================================
sequenz_1 = "ATGGCCATTGTAATGGGCCGCTGA"
sequenz_2 = "ATGGCCATTGTATGGGCCGCTGA"    # gleiche Sequenz, aber ein A fehlt (Deletion)

MATCH = 1
MISMATCH = -1
GAP = -2

MODUS = "global"      # "global" (Needleman-Wunsch) oder "lokal" (Smith-Waterman)
# ============================================================


def saeubern(sequenz):
    return "".join(sequenz.split()).upper()


def bewerte(a, b, match, mismatch):
    if a == b:
        return match
    return mismatch


def leere_matrix(zeilen, spalten):
    # Eine Tabelle (Liste von Listen) voller Nullen
    matrix = []
    for i in range(zeilen):
        matrix.append([0] * spalten)
    return matrix


def needleman_wunsch(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP):
    """Globales Alignment: beide Sequenzen werden vollständig von Anfang bis Ende aligniert."""
    n = len(s1)
    m = len(s2)
    F = leere_matrix(n + 1, m + 1)   # F[i][j] = bester Score für s1[:i] gegen s2[:j]

    # 1. Erste Spalte und erste Zeile: nur Lücken
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap

    # 2. Matrix füllen: jede Zelle kommt von links, oben oder schräg oben-links
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = F[i - 1][j - 1] + bewerte(s1[i - 1], s2[j - 1], match, mismatch)
            oben = F[i - 1][j] + gap        # Lücke in s2
            links = F[i][j - 1] + gap       # Lücke in s1
            F[i][j] = max(diagonal, oben, links)

    # 3. Traceback: von unten rechts zurück nach oben links
    a1 = ""
    a2 = ""
    i = n
    j = m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + bewerte(s1[i - 1], s2[j - 1], match, mismatch):
            a1 = s1[i - 1] + a1
            a2 = s2[j - 1] + a2
            i = i - 1
            j = j - 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1 = s1[i - 1] + a1
            a2 = "-" + a2
            i = i - 1
        else:
            a1 = "-" + a1
            a2 = s2[j - 1] + a2
            j = j - 1

    return {"score": F[n][m], "a1": a1, "a2": a2, "matrix": F, "start1": 1, "start2": 1}


def smith_waterman(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP):
    """Lokales Alignment: findet den am besten passenden Abschnitt beider Sequenzen."""
    n = len(s1)
    m = len(s2)
    H = leere_matrix(n + 1, m + 1)   # erste Zeile/Spalte bleiben 0

    bester_score = 0
    beste_i = 0
    beste_j = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = H[i - 1][j - 1] + bewerte(s1[i - 1], s2[j - 1], match, mismatch)
            oben = H[i - 1][j] + gap
            links = H[i][j - 1] + gap
            H[i][j] = max(0, diagonal, oben, links)   # Unterschied zu NW: nie unter 0
            if H[i][j] > bester_score:
                bester_score = H[i][j]
                beste_i = i
                beste_j = j

    # Traceback: vom besten Feld zurück, bis eine 0 erreicht ist
    a1 = ""
    a2 = ""
    i = beste_i
    j = beste_j
    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i - 1][j - 1] + bewerte(s1[i - 1], s2[j - 1], match, mismatch):
            a1 = s1[i - 1] + a1
            a2 = s2[j - 1] + a2
            i = i - 1
            j = j - 1
        elif H[i][j] == H[i - 1][j] + gap:
            a1 = s1[i - 1] + a1
            a2 = "-" + a2
            i = i - 1
        else:
            a1 = "-" + a1
            a2 = s2[j - 1] + a2
            j = j - 1

    # i, j zeigen jetzt auf die Stelle direkt vor dem Alignment → +1 für biologische Zählung ab 1
    return {"score": bester_score, "a1": a1, "a2": a2, "matrix": H, "start1": i + 1, "start2": j + 1}


def zeige_matrix(s1, s2, matrix):
    """Druckt die Score-Matrix (nur sinnvoll für kurze Sequenzen)."""
    print("       " + "".join(f"{b:>4}" for b in "-" + s2))
    for i in range(len(matrix)):
        zeilen_base = "-" if i == 0 else s1[i - 1]
        print(f"{zeilen_base:>4}   " + "".join(f"{wert:>4}" for wert in matrix[i]))


def zeige_alignment(ergebnis, breite=60):
    a1 = ergebnis["a1"]
    a2 = ergebnis["a2"]
    mitte = ""
    for x, y in zip(a1, a2):
        if x == y:
            mitte = mitte + "|"      # Match
        elif x == "-" or y == "-":
            mitte = mitte + " "      # Lücke
        else:
            mitte = mitte + "."      # Mismatch

    pos1 = ergebnis["start1"]
    pos2 = ergebnis["start2"]
    for k in range(0, len(a1), breite):
        teil1 = a1[k:k + breite]
        teil2 = a2[k:k + breite]
        print(f"Seq1 {pos1:>6}  {teil1}")
        print(f"            {mitte[k:k + breite]}")
        print(f"Seq2 {pos2:>6}  {teil2}")
        print()
        pos1 = pos1 + len(teil1) - teil1.count("-")
        pos2 = pos2 + len(teil2) - teil2.count("-")

    laenge = len(a1)
    matches = mitte.count("|")
    mismatches = mitte.count(".")
    luecken = mitte.count(" ")
    print(f"Score:       {ergebnis['score']}")
    if laenge > 0:
        print(f"Länge:       {laenge}")
        print(f"Identität:   {matches}/{laenge} ({matches / laenge * 100:.1f} %)")
        print(f"Mismatches:  {mismatches}")
        print(f"Lücken:      {luecken}")


def aligniere(s1, s2, modus="global", match=MATCH, mismatch=MISMATCH, gap=GAP):
    s1 = saeubern(s1)
    s2 = saeubern(s2)
    if len(s1) == 0 or len(s2) == 0:
        print("⚠️ Eine der Sequenzen ist leer.")
        return None
    if len(s1) * len(s2) > 25_000_000:
        print("⚠️ Achtung: sehr lange Sequenzen, das kann in reinem Python mehrere Minuten dauern.")

    if modus == "global":
        print("Globales Alignment (Needleman-Wunsch)")
        ergebnis = needleman_wunsch(s1, s2, match, mismatch, gap)
    elif modus == "lokal":
        print("Lokales Alignment (Smith-Waterman)")
        ergebnis = smith_waterman(s1, s2, match, mismatch, gap)
    else:
        print('⚠️ MODUS muss "global" oder "lokal" sein.')
        return None
    print(f"Bewertung: Match {match:+}, Mismatch {mismatch:+}, Lücke {gap:+}\n")

    if len(s1) <= 15 and len(s2) <= 15:
        print("Score-Matrix:")
        zeige_matrix(s1, s2, ergebnis["matrix"])
        print()
    zeige_alignment(ergebnis)
    return ergebnis


ergebnis = aligniere(sequenz_1, sequenz_2, MODUS)
