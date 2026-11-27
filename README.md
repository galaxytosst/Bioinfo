# Bioinformatik – Lernprojekt

Schritt für Schritt Bioinformatik mit Python lernen, an echten biologischen Beispielen.

| Ordner | Inhalt |
|---|---|
| `00_python_grundlagen/` | Teil 0: Python-Grundlagen (Basen zählen, Transkription, GC-Gehalt, Reverse Complement, Translation, FASTA) + Sequenz-Umwandler + Strang-Vergleich |
| `01_alignment/` | Teil 1: Paarweises Alignment, Needleman-Wunsch (global) und Smith-Waterman (lokal) |
| `02_orf_finder/` | Teil 2: ORF-Finder für alle 6 Leseraster, getestet an den 10 Genen von pPCP1 |
| `03_phylogenie/` | Teil 3: Stammbaum aus Distanzen (UPGMA, Neighbor-Joining) für 12 Frauenschuh-Orchideen |
| `beispiele/` | 72 Beispielprojekte (Liste mit Beschreibung: `beispiele/UEBERSICHT.md`), je ein Unterordner mit `LIESMICH.md` und `einstellungen.json`. Themen: Prionen, Erbkrankheiten, Variation beim Menschen, ein Gen in vielen Arten, Genfamilien, Viren, Bakterien, Mitochondrien-DNA, Pflanzen. Sequenzen stehen als FASTA im Ordner oder als NCBI-Nummern in `accessions.txt` (werden beim ersten Start geladen). `varianten.txt` baut bekannte Mutationen (HGVS, z. B. `c.20A>T`) in eine echte Referenz ein |
| `daten/` | Platz für eigene Sequenzen (ein Unterordner pro Projekt) |
| `werkzeuge.py` | Alle Funktionen als Modul (FASTA, Transkription, Translation, ORFs, Alignment inkl. affiner Lückenstrafe, Protein-Alignment mit BLOSUM62, Distanzen, Bäume, Varianten) |
| `analyse.py` | Automatische Auswertung: `analysiere_ordner()` für ein Projekt, `analysiere_alle()` für viele Projekte; mit `SCHRITTE` wählbar: alles, ORFs, DNA-Vergleich, Mutationen, Protein-Vergleich, 3D-Struktur, Stammbaum (3D: Struktur aus der PDB, automatisch gesucht oder mit `STRUKTUR` vorgegeben) |
| `bericht.py` | Gemeinsamer HTML-Bericht mit Tabellen und Diagrammen |
| `App.ipynb` | Bioinformatik-App: alle Werkzeuge als Colab-Formulare (Auswahllisten, Textfelder), Code ausgeblendet. Box 1 = Analyse mit Auswahl der Analyseart (inkl. 3D-Struktur), Box 6 = Protein-Vergleich und Eigenschaften, Box 7 = 3D-Struktur |
| `struktur.py` | 3D-Strukturen: laden (PDB, AlphaFold), Bindestellen, Oberfläche je Aminosäure (Shrake-Rupley), UniProt-Stellen |
| `oberflaeche.py` | Klick-Oberfläche für Colab: Analyse, Sequenz umwandeln, ORFs, Vergleichen (DNA und Protein), Mutations-Labor, NCBI-Download |

**Bedienung ohne Code (Formulare):** [App.ipynb in Colab öffnen](https://colab.research.google.com/github/galaxytosst/Bioinfo/blob/main/App.ipynb)

**Alle Werkzeuge mit Code und Erklärungen:** [Werkzeuge.ipynb in Colab öffnen](https://colab.research.google.com/github/galaxytosst/Bioinfo/blob/main/Werkzeuge.ipynb)

**Andere Notebooks öffnen:** colab.research.google.com → *Datei → Notebook öffnen* → Reiter **GitHub** → dieses Repository auswählen.

Beispieldaten: *Yersinia pestis*-Plasmid pPCP1 (NCBI RefSeq NC_005816.1).
