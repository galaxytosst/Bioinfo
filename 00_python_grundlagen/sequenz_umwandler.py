# 🧬 Sequenz-Umwandler: DNA ↔ RNA und Protein, alles in einem Schritt
# Trag unten bei "eingabe" deine Sequenz ein (DNA oder RNA) und führe die Zelle aus.

eingabe = "ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"

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
    # Leerzeichen und Zeilenumbrüche entfernen, alles in Großbuchstaben
    return "".join(sequenz.split()).upper()


def dna_zu_rna(dna):
    return dna.replace("T", "U")


def rna_zu_dna(rna):
    return rna.replace("U", "T")


def reverse_complement(dna):
    komplement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    ergebnis = ""
    for base in dna:
        ergebnis = komplement[base] + ergebnis
    return ergebnis


def uebersetze(rna):
    protein = ""
    for i in range(0, len(rna) - 2, 3):
        aminosaeure = codon_tabelle.get(rna[i:i+3], "X")   # X = unbekanntes Codon (z. B. mit N)
        if aminosaeure == "*":
            break
        protein = protein + aminosaeure
    return protein


def umwandeln(sequenz):
    sequenz = saeubern(sequenz)
    if "U" in sequenz and "T" in sequenz:
        print("⚠️ Die Sequenz enthält T und U. Das ist weder DNA noch RNA, bitte prüfen.")
        return
    erlaubt = set("ACGTUN")
    falsch = set(sequenz) - erlaubt
    if falsch:
        print("⚠️ Unerlaubte Zeichen:", "".join(sorted(falsch)))
        return

    if "U" in sequenz:
        print("Erkannt: RNA")
        rna = sequenz
        dna = rna_zu_dna(rna)
    else:
        print("Erkannt: DNA (als codierender Strang gelesen)")
        dna = sequenz
        rna = dna_zu_rna(dna)

    protein = uebersetze(rna)
    gc = (dna.count("G") + dna.count("C")) / len(dna) * 100

    print(f"Länge:               {len(dna)} Basen")
    print(f"GC-Gehalt:           {gc:.1f} %")
    print(f"DNA (codierend):     {dna}")
    print(f"RNA (mRNA):          {rna}")
    print(f"Reverse Complement:  {reverse_complement(dna)}")
    print(f"Protein:             {protein}  ({len(protein)} Aminosäuren)")
    if len(dna) % 3 != 0:
        print(f"Hinweis: Die Länge ist nicht durch 3 teilbar, die letzten {len(dna) % 3} Basen werden nicht übersetzt.")
    if not rna.startswith("AUG"):
        print("Hinweis: Die Sequenz beginnt nicht mit dem Startcodon AUG. Übersetzt wurde trotzdem ab der ersten Base.")


umwandeln(eingabe)
