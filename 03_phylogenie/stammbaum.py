# 🌳 Stammbaum aus DNA-Sequenzen: Alignment → Distanzen → Baum (UPGMA und Neighbor-Joining)
#
# Ablauf:
#   1. Jedes Sequenzpaar wird mit Needleman-Wunsch aligniert (dein Code aus Teil 1).
#   2. Aus jedem Alignment wird eine evolutionäre Distanz berechnet (p-Distanz, Jukes-Cantor-Korrektur).
#   3. Aus der Distanzmatrix wird ein Baum gebaut: UPGMA und Neighbor-Joining.
#
# Beispieldaten: 5.8S-rRNA + ITS1/ITS2 von 12 Frauenschuh-Orchideen (Unterfamilie Cypripedioideae),
# 4 Gattungen. Quelle: GenBank (Accession-Nummern stehen in der FASTA-Datei), ausgewählt aus dem
# Beispieldatensatz "ls_orchid.fasta" des Biopython-Projekts.
#
# So benutzt du es: Ganzen Code in eine Colab-Zelle kopieren und ausführen (dauert ca. 15–30 s).

import math
import os
import urllib.request

# ============================================================
# EINSTELLUNGEN
# ============================================================
FASTA_DATEI = "orchideen_ITS.fasta"
FASTA_URL = "https://raw.githubusercontent.com/galaxytosst/Bioinfo/main/03_phylogenie/orchideen_ITS.fasta"
AUSSENGRUPPE = ["Cypripedium_calceolus", "Cypripedium_acaule", "Cypripedium_reginae"]   # zum Bewurzeln des NJ-Baums
KORREKTUR = "jukes-cantor"               # "jukes-cantor" oder "p" (unkorrigierte p-Distanz)

MATCH, MISMATCH, GAP = 1, -1, -2
# ============================================================


# ---------- Sequenzen einlesen ----------
def lese_fasta(dateiname):
    sequenzen = {}
    name = None
    with open(dateiname) as f:
        for zeile in f:
            zeile = zeile.strip()
            if zeile == "":
                continue
            if zeile.startswith(">"):
                name = zeile[1:].split()[0]      # nur das erste Wort als Name
                sequenzen[name] = ""
            else:
                sequenzen[name] = sequenzen[name] + zeile.upper()
    return sequenzen


# ---------- 1. Alignment (Needleman-Wunsch aus Teil 1, kompakt) ----------
def needleman_wunsch(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP):
    n, m = len(s1), len(s2)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            F[i][j] = max(diag, F[i - 1][j] + gap, F[i][j - 1] + gap)
    a1, a2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
            a1, a2, i, j = s1[i - 1] + a1, s2[j - 1] + a2, i - 1, j - 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1, a2, i = s1[i - 1] + a1, "-" + a2, i - 1
        else:
            a1, a2, j = "-" + a1, s2[j - 1] + a2, j - 1
    return a1, a2


# ---------- 2. Distanz ----------
def distanz(a1, a2, korrektur=KORREKTUR):
    """p-Distanz = Anteil verschiedener Basen an allen Positionen ohne Lücke.
    Jukes-Cantor korrigiert für Mehrfachmutationen an derselben Stelle: d = -3/4 · ln(1 - 4/3 · p)."""
    vergleichbar = 0
    verschieden = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            continue                      # Lücken werden hier nicht gezählt
        vergleichbar = vergleichbar + 1
        if x != y:
            verschieden = verschieden + 1
    p = verschieden / vergleichbar
    if korrektur == "p":
        return p
    if p >= 0.75:
        return float("inf")              # Sequenzen sind "gesättigt", Distanz nicht schätzbar
    return -0.75 * math.log(1 - 4 / 3 * p)


def distanzmatrix(sequenzen):
    namen = list(sequenzen)
    D = {a: {a: 0.0} for a in namen}
    paare = len(namen) * (len(namen) - 1) // 2
    nr = 0
    for i in range(len(namen)):
        for j in range(i + 1, len(namen)):
            a, b = namen[i], namen[j]
            a1, a2 = needleman_wunsch(sequenzen[a], sequenzen[b])
            D[a][b] = D[b][a] = distanz(a1, a2)
            nr = nr + 1
            if nr % 11 == 0 or nr == paare:
                print(f"Aligniere Paar {nr}/{paare} ...")
    print("fertig.\n")
    return namen, D


def kurz(name):
    # "Cypripedium_calceolus" → "C.calceolus", "Mexipedium_xerophyticum" → "M.xerophyticum"
    teile = name.split("_")
    if len(teile) == 2:
        gattung = teile[0]
        kuerzel = {"Paphiopedilum": "Paph", "Phragmipedium": "Phrag", "Cypripedium": "Cyp", "Mexipedium": "Mex"}
        return kuerzel.get(gattung, gattung[0]) + "." + teile[1]
    return name


def zeige_distanzmatrix(namen, D):
    breite = max(len(kurz(n)) for n in namen)
    print("Distanzmatrix (" + KORREKTUR + "):")
    print(" " * (breite + 2) + "".join(f"{i + 1:>6}" for i in range(len(namen))))
    for i, a in enumerate(namen):
        print(f"{i + 1:>2} {kurz(a):<{breite}}" + "".join(f"{D[a][b]:>6.3f}" for b in namen))
    print()


# ---------- 3a. UPGMA ----------
# Ein Knoten im Baum ist ein Dictionary: {"name": ..., "kinder": [(kindknoten, astlaenge), ...]}
def blatt(name):
    return {"name": name, "kinder": []}


def upgma(namen, D):
    """Fasst immer die zwei ähnlichsten Gruppen zusammen. Annahme: konstante Mutationsrate (molekulare Uhr)."""
    gruppen = {n: {"knoten": blatt(n), "groesse": 1, "hoehe": 0.0} for n in namen}
    D = {a: dict(D[a]) for a in D}                       # Kopie, damit das Original erhalten bleibt
    zaehler = 0
    while len(gruppen) > 1:
        liste = list(gruppen)
        a, b = min(((x, y) for i, x in enumerate(liste) for y in liste[i + 1:]), key=lambda p: D[p[0]][p[1]])
        hoehe = D[a][b] / 2
        A, B = gruppen.pop(a), gruppen.pop(b)
        zaehler = zaehler + 1
        neu = f"U{zaehler}"
        knoten = {"name": None, "kinder": [(A["knoten"], hoehe - A["hoehe"]), (B["knoten"], hoehe - B["hoehe"])]}
        D[neu] = {neu: 0.0}
        for k in gruppen:                                # Durchschnitt, gewichtet nach Gruppengröße
            d = (D[a][k] * A["groesse"] + D[b][k] * B["groesse"]) / (A["groesse"] + B["groesse"])
            D[neu][k] = D[k][neu] = d
        gruppen[neu] = {"knoten": knoten, "groesse": A["groesse"] + B["groesse"], "hoehe": hoehe}
    return list(gruppen.values())[0]["knoten"]


# ---------- 3b. Neighbor-Joining ----------
def neighbor_joining(namen, D, aussengruppe):
    """Braucht keine molekulare Uhr. Ergibt einen unbewurzelten Baum, der hier an der Außengruppe bewurzelt wird."""
    D = {a: dict(D[a]) for a in D}
    aktiv = list(namen)
    nachbarn = {n: [] for n in namen}                    # ungerichteter Graph: Knoten → [(Nachbar, Länge)]
    zaehler = 0
    while len(aktiv) > 2:
        n = len(aktiv)
        r = {i: sum(D[i][k] for k in aktiv) for i in aktiv}
        i, j = min(((x, y) for a, x in enumerate(aktiv) for y in aktiv[a + 1:]),
                   key=lambda p: (n - 2) * D[p[0]][p[1]] - r[p[0]] - r[p[1]])
        li = D[i][j] / 2 + (r[i] - r[j]) / (2 * (n - 2))
        lj = D[i][j] - li
        zaehler = zaehler + 1
        u = f"N{zaehler}"
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
        kinder = [(teilbaum(nb, knoten), l) for nb, l in nachbarn[knoten] if nb != herkunft]
        return {"name": knoten if knoten in namen else None, "kinder": kinder}

    def blaetter_hinter(knoten, herkunft):
        # alle Blätter, die man von 'knoten' aus erreicht, ohne zurück zu 'herkunft' zu gehen
        if knoten in namen:
            return {knoten}
        menge = set()
        for nb, _ in nachbarn[knoten]:
            if nb != herkunft:
                menge = menge | blaetter_hinter(nb, knoten)
        return menge

    # Wurzel in die Mitte des Asts setzen, der die Außengruppe vom Rest trennt
    ziel = set(aussengruppe)
    for u in nachbarn:
        for v, laenge in nachbarn[u]:
            if blaetter_hinter(v, u) == ziel:
                return {"name": None, "kinder": [(teilbaum(v, u), laenge / 2), (teilbaum(u, v), laenge / 2)]}
    print("⚠️ Die Außengruppe bildet im Baum keine eigene Gruppe. Bewurzelt wird an ihrem ersten Mitglied.")
    erstes = aussengruppe[0]
    elternteil, laenge = nachbarn[erstes][0]
    return {"name": None, "kinder": [(blatt(erstes), laenge / 2), (teilbaum(elternteil, erstes), laenge / 2)]}


# ---------- Ausgabe ----------
def newick(knoten, laenge=None):
    if knoten["kinder"]:
        text = "(" + ",".join(newick(k, l) for k, l in knoten["kinder"]) + ")"
    else:
        text = knoten["name"]
    if laenge is not None:
        text = text + f":{laenge:.4f}"
    return text


def zeige_baum_text(knoten, praefix="", laenge=None, letzter=True, wurzel=True):
    beschriftung = kurz(knoten["name"]) if knoten["name"] else "┐"
    if wurzel:
        print(beschriftung)
    else:
        print(praefix + ("└── " if letzter else "├── ") + beschriftung + f"  ({laenge:.3f})")
    neues_praefix = praefix + ("" if wurzel else ("    " if letzter else "│   "))
    for nr, (kind, l) in enumerate(knoten["kinder"]):
        zeige_baum_text(kind, neues_praefix, l, nr == len(knoten["kinder"]) - 1, False)


def blaetter(knoten):
    if not knoten["kinder"]:
        return [knoten["name"]]
    ergebnis = []
    for kind, _ in knoten["kinder"]:
        ergebnis = ergebnis + blaetter(kind)
    return ergebnis


def zeichne_baum(knoten, titel, ax):
    """Zeichnet einen rechtwinkligen Baum mit matplotlib."""
    farben = {"Cypripedium": "tab:green", "Paphiopedilum": "tab:purple", "Phragmipedium": "tab:orange", "Mexipedium": "tab:red"}
    y_position = {name: i for i, name in enumerate(blaetter(knoten))}

    def zeichne(k, x):
        if not k["kinder"]:
            y = y_position[k["name"]]
            farbe = farben.get(k["name"].split("_")[0], "black")
            ax.text(x + 0.002, y, " " + k["name"].replace("_", " "), va="center", fontsize=9, color=farbe, style="italic")
            return y
        ys = []
        for kind, l in k["kinder"]:
            y_kind = zeichne(kind, x + l)
            ax.plot([x, x + l], [y_kind, y_kind], color="black", lw=1.2)
            ys.append(y_kind)
        ax.plot([x, x], [min(ys), max(ys)], color="black", lw=1.2)
        return sum(ys) / len(ys)

    zeichne(knoten, 0.0)
    ax.set_title(titel)
    ax.set_yticks([])
    ax.set_xlabel("Distanz (Substitutionen pro Position)")
    for seite in ["left", "right", "top"]:
        ax.spines[seite].set_visible(False)
    x_max = ax.get_xlim()[1]
    ax.set_xlim(0, x_max * 1.6)                         # Platz für die Namen


# ---------- Hauptprogramm ----------
if not os.path.exists(FASTA_DATEI):
    urllib.request.urlretrieve(FASTA_URL, FASTA_DATEI)
sequenzen = lese_fasta(FASTA_DATEI)
print(f"{len(sequenzen)} Sequenzen geladen:")
for name, seq in sequenzen.items():
    print(f"  {name:<32} {len(seq)} bp")
print()

namen, D = distanzmatrix(sequenzen)
zeige_distanzmatrix(namen, D)

baum_upgma = upgma(namen, D)
baum_nj = neighbor_joining(namen, D, AUSSENGRUPPE)

print("=== UPGMA ===")
zeige_baum_text(baum_upgma)
print("\nNewick:", newick(baum_upgma) + ";")

print("\n=== Neighbor-Joining (bewurzelt mit Außengruppe: " + ", ".join(kurz(n) for n in AUSSENGRUPPE) + ") ===")
zeige_baum_text(baum_nj)
print("\nNewick:", newick(baum_nj) + ";")

try:
    import matplotlib.pyplot as plt
    fig, achsen = plt.subplots(1, 2, figsize=(15, 6))
    zeichne_baum(baum_upgma, "UPGMA", achsen[0])
    zeichne_baum(baum_nj, "Neighbor-Joining", achsen[1])
    plt.tight_layout()
    plt.show()
except ImportError:
    print("\n(matplotlib nicht installiert, daher keine Grafik.)")
