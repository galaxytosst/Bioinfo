# 🧬 Zwei DNA-Stränge vergleichen
# Trag unten zwei Sequenzen ein (DNA oder RNA, je 5'→3') und führe die Zelle aus.
#
# Der Code beantwortet zwei Fragen:
#  1. Passen die beiden Stränge als Doppelstrang zusammen (sind sie komplementär)?
#  2. Wo unterscheiden sie sich (z. B. Mutationen) und was bewirkt das im Protein?

strang_1 = "ATGGCCATTGTAATGGGCCGCTGA"
strang_2 = "ATGGCTATTGTAATGAGCCGCTGA"

basen = "UCAG"
aminosaeuren = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_tabelle = {}
nummer = 0
for b1 in basen:
    for b2 in basen:
        for b3 in basen:
            codon_tabelle[b1 + b2 + b3] = aminosaeuren[nummer]
            nummer = nummer + 1


def saeubern(sequenz):
    # Leerzeichen entfernen, Großbuchstaben, RNA wird als DNA behandelt (U → T)
    return "".join(sequenz.split()).upper().replace("U", "T")


def reverse_complement(dna):
    komplement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    ergebnis = ""
    for base in dna:
        ergebnis = komplement[base] + ergebnis
    return ergebnis


def uebersetze(dna):
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        aminosaeure = codon_tabelle.get(dna[i:i+3].replace("T", "U"), "X")
        protein = protein + aminosaeure      # Stopp wird hier als * mit angezeigt
    return protein


def mutationstyp(base_1, base_2):
    # Transition: Purin ↔ Purin (A↔G) oder Pyrimidin ↔ Pyrimidin (C↔T)
    # Transversion: Purin ↔ Pyrimidin
    purine = "AG"
    if (base_1 in purine) == (base_2 in purine):
        return "Transition"
    return "Transversion"


def vergleiche(s1, s2):
    s1 = saeubern(s1)
    s2 = saeubern(s2)
    print(f"Strang 1: {len(s1)} Basen")
    print(f"Strang 2: {len(s2)} Basen")

    # --- Frage 1: Doppelstrang? ---
    print("\n=== 1. Sind die Stränge komplementär (bilden sie einen Doppelstrang)? ===")
    if s2 == reverse_complement(s1):
        print("✅ Ja. Strang 2 ist das Reverse Complement von Strang 1.")
        print("   5'-" + s1 + "-3'")
        print("   3'-" + s2[::-1] + "-5'")
        return
    print("Nein, die beiden sind kein perfekter Doppelstrang. Ich vergleiche sie deshalb als zwei Versionen derselben Sequenz.")

    # --- Frage 2: Unterschiede ---
    print("\n=== 2. Unterschiede Position für Position ===")
    if len(s1) != len(s2):
        print("⚠️ Die Stränge sind unterschiedlich lang. Ein Vergleich Base für Base ist dann nicht sinnvoll,")
        print("   weil eine Insertion oder Deletion alles danach verschiebt. Dafür braucht man ein Alignment")
        print("   (Needleman-Wunsch), das ist unser nächstes Thema.")
        return

    markierung = ""
    unterschiede = []
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            markierung = markierung + "|"
        else:
            markierung = markierung + " "
            unterschiede.append(i)

    print("Strang 1: " + s1)
    print("          " + markierung)
    print("Strang 2: " + s2)

    identitaet = (len(s1) - len(unterschiede)) / len(s1) * 100
    print(f"\nUnterschiede (Hamming-Distanz): {len(unterschiede)}")
    print(f"Identität: {identitaet:.1f} %")

    if not unterschiede:
        print("Die beiden Sequenzen sind identisch.")
        return

    print("\nPosition  Änderung  Typ           Codon       Aminosäure   Wirkung")
    for i in unterschiede:
        codon_start = (i // 3) * 3                 # Beginn des Codons, in dem die Mutation liegt
        codon_1 = s1[codon_start:codon_start + 3]
        codon_2 = s2[codon_start:codon_start + 3]
        if len(codon_1) < 3:
            print(f"{i + 1:>8}  {s1[i]}→{s2[i]}       {mutationstyp(s1[i], s2[i]):<12}  (kein vollständiges Codon)")
            continue
        as_1 = codon_tabelle.get(codon_1.replace("T", "U"), "X")
        as_2 = codon_tabelle.get(codon_2.replace("T", "U"), "X")
        if as_1 == as_2:
            wirkung = "stumm (gleiche Aminosäure)"
        elif as_2 == "*":
            wirkung = "Nonsense (neues Stoppcodon!)"
        elif as_1 == "*":
            wirkung = "Stoppcodon verloren"
        else:
            wirkung = "Missense (andere Aminosäure)"
        print(f"{i + 1:>8}  {s1[i]}→{s2[i]}       {mutationstyp(s1[i], s2[i]):<12}  {codon_1}→{codon_2}     {as_1}→{as_2}          {wirkung}")

    print("\nProtein 1: " + uebersetze(s1))
    print("Protein 2: " + uebersetze(s2))
    print("(Positionen zählen ab 1, wie in der Biologie. Gelesen wird ab der ersten Base im Leseraster 1.)")


vergleiche(strang_1, strang_2)
