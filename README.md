# Bioinformatik – Lernprojekt

Schritt für Schritt Bioinformatik mit Python lernen, an echten biologischen Beispielen.

| Ordner | Inhalt |
|---|---|
| `00_python_grundlagen/` | Teil 0: Python-Grundlagen (Basen zählen, Transkription, GC-Gehalt, Reverse Complement, Translation, FASTA) + Sequenz-Umwandler + Strang-Vergleich |
| `01_alignment/` | Teil 1: Paarweises Alignment, Needleman-Wunsch (global) und Smith-Waterman (lokal) |
| `02_orf_finder/` | Teil 2: ORF-Finder für alle 6 Leseraster, getestet an den 10 Genen von pPCP1 |
| `03_phylogenie/` | Teil 3: Stammbaum aus Distanzen (UPGMA, Neighbor-Joining) für 12 Frauenschuh-Orchideen |
| `beispiele/` | Sechs Beispielprojekte, je ein Unterordner mit FASTA-Dateien, `LIESMICH.md` und `einstellungen.json`: Sichelzellanämie, Hämoglobin von 7 Arten, Coronaviren, Frauenschuh-Orchideen, Pest-Plasmid, Chloroplasten-Genom |
| `daten/` | Platz für eigene Sequenzen (ein Unterordner pro Projekt) |
| `werkzeuge.py` | Alle Funktionen als Modul (FASTA, Transkription, Translation, ORFs, Alignment inkl. affiner Lückenstrafe, Distanzen, Bäume, Varianten) |
| `analyse.py` | Automatische Auswertung: `analysiere_ordner()` für ein Projekt, `analysiere_alle()` für viele Projekte |
| `bericht.py` | Gemeinsamer HTML-Bericht mit Tabellen und Diagrammen |
| `App.ipynb` | Bioinformatik-App: alle Werkzeuge als Colab-Formulare (Auswahllisten, Textfelder), Code ausgeblendet |
| `oberflaeche.py` | Klick-Oberfläche für Colab: Analyse, Sequenz umwandeln, ORFs, Vergleichen, Mutations-Labor, NCBI-Download |

**Bedienung ohne Code (Formulare):** [App.ipynb in Colab öffnen](https://colab.research.google.com/github/galaxytosst/Bioinfo/blob/main/App.ipynb)

**Alle Werkzeuge mit Code und Erklärungen:** [Werkzeuge.ipynb in Colab öffnen](https://colab.research.google.com/github/galaxytosst/Bioinfo/blob/main/Werkzeuge.ipynb)

**Andere Notebooks öffnen:** colab.research.google.com → *Datei → Notebook öffnen* → Reiter **GitHub** → dieses Repository auswählen.

Beispieldaten: *Yersinia pestis*-Plasmid pPCP1 (NCBI RefSeq NC_005816.1).
