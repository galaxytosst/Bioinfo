"""
🖥️ oberflaeche.py – Klick-Oberfläche für alle Werkzeuge (Google Colab oder Jupyter).

Start in einer Notebook-Zelle:

    import oberflaeche
    oberflaeche.starte()

Reiter:
  📊 Analyse          ganze Projektordner auswerten → ein Bericht mit Tabellen und Diagrammen
  🧬 Eine Sequenz      umwandeln (DNA, mRNA, Gegenstrang, Protein) und ORFs suchen
  ⚖️ Vergleichen       zwei Sequenzen alignieren, Unterschiede und ihre Wirkung
  🩸 Mutations-Labor   selbst eine Mutation in ein Gen einbauen und die Folgen sehen
  🌐 NCBI              Sequenzen per Accession-Nummer herunterladen
"""

import base64
import html
import io
import os
import re
import shutil
import traceback

import ipywidgets as W
from IPython.display import HTML, display

import analyse
import werkzeuge as wz

HIER = os.path.dirname(os.path.abspath(__file__))
BEISPIELE = os.path.join(HIER, "beispiele")
EIGENE = os.path.abspath("daten")                     # in Colab: /content/daten
DRIVE_STANDARD = "/content/drive/MyDrive/Bioinformatik/daten"

KURZES_TESTGEN = "ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
MAX_ZELLEN_VERGLEICH = 400_000_000                    # Länge1 × Länge2, ca. 15 s

STARTCODON_OPTIONEN = [("nur ATG", ("ATG",)), ("ATG, GTG, TTG (Bakterien)", ("ATG", "GTG", "TTG"))]


# ======================================================================
# Kleine Helfer für die Anzeige
# ======================================================================

def _e(text):
    return html.escape(str(text))


def _kasten(inhalt):
    """Weißer Kasten, damit Ergebnisse auch im dunklen Colab-Design lesbar bleiben."""
    return (f'<div style="background:#fff;color:#1d1d1f;border:1px solid #e3e3e3;border-radius:8px;'
            f'padding:12px 16px;margin:6px 0;font-family:system-ui,-apple-system,Segoe UI,sans-serif;'
            f'font-size:14px;line-height:1.45;">{inhalt}</div>')


def _hinweis(text, art="info"):
    farben = {"info": ("#eef5fb", "#0072B2"), "ok": ("#eaf7f1", "#009E73"),
              "warnung": ("#fff6e6", "#E69F00"), "fehler": ("#fdecea", "#D55E00")}
    hg, rand = farben[art]
    return (f'<div style="background:{hg};color:#1d1d1f;border-left:4px solid {rand};padding:8px 12px;'
            f'margin:6px 0;border-radius:4px;font-family:system-ui,sans-serif;font-size:14px;">{text}</div>')


def _pre(text):
    return (f'<pre style="background:#f6f7f9;color:#1d1d1f;padding:10px;border-radius:6px;overflow-x:auto;'
            f'font-size:12.5px;line-height:1.35;white-space:pre;">{_e(text)}</pre>')


def _umbrechen(seq, breite=60, max_zeichen=3000):
    teil = seq[:max_zeichen]
    zeilen = "\n".join(teil[i:i + breite] for i in range(0, len(teil), breite))
    if len(seq) > max_zeichen:
        zeilen += f"\n… ({len(seq) - max_zeichen:,} weitere Zeichen, vollständig im Download)"
    return zeilen


def _tabelle(kopf, zeilen, rechts=()):
    th = "".join(f'<th style="text-align:{"right" if i in rechts else "left"};padding:5px 10px;'
                 f'border-bottom:2px solid #ddd;background:#f6f7f9;">{_e(k)}</th>' for i, k in enumerate(kopf))
    tr = ""
    for z in zeilen:
        tr += "<tr>" + "".join(
            f'<td style="text-align:{"right" if i in rechts else "left"};padding:4px 10px;'
            f'border-bottom:1px solid #eee;">{w}</td>' for i, w in enumerate(z)) + "</tr>"
    return f'<table style="border-collapse:collapse;font-size:13px;margin:6px 0;"><tr>{th}</tr>{tr}</table>'


def _bild(fig):
    puffer = io.BytesIO()
    fig.savefig(puffer, format="png", dpi=110, bbox_inches="tight", facecolor="white")
    daten = base64.b64encode(puffer.getvalue()).decode()
    return f'<img style="max-width:100%;" src="data:image/png;base64,{daten}">'


def _download_link(pfad, text, mime="application/octet-stream"):
    """Download-Knopf mit der Datei direkt im Link (funktioniert ohne Zusatzpakete)."""
    groesse = os.path.getsize(pfad)
    if groesse > 40 * 1024 * 1024:
        return (f"{_e(text)}: zu groß für einen Download-Link ({groesse / 1e6:.0f} MB), "
                f"liegt im Dateibereich links unter <code>{_e(pfad)}</code>")
    with open(pfad, "rb") as f:
        daten = base64.b64encode(f.read()).decode()
    return (f'<a download="{_e(os.path.basename(pfad))}" href="data:{mime};base64,{daten}" '
            f'style="display:inline-block;background:#0072B2;color:#fff;padding:6px 12px;border-radius:5px;'
            f'text-decoration:none;margin:4px 8px 4px 0;font-family:system-ui,sans-serif;font-size:13px;">'
            f'⬇ {_e(text)} ({groesse / 1e6:.1f} MB)</a>')


def _text_download(text, dateiname, beschriftung):
    daten = base64.b64encode(text.encode()).decode()
    return (f'<a download="{_e(dateiname)}" href="data:text/plain;base64,{daten}" '
            f'style="display:inline-block;background:#0072B2;color:#fff;padding:5px 10px;border-radius:5px;'
            f'text-decoration:none;margin:4px 8px 4px 0;font-size:13px;">⬇ {_e(beschriftung)}</a>')


def _bericht_einbetten(pfad, hoehe="900px"):
    with open(pfad, encoding="utf-8") as f:
        inhalt = f.read()
    return (f'<iframe srcdoc="{html.escape(inhalt, quote=True)}" '
            f'style="width:100%;height:{hoehe};border:1px solid #ccc;border-radius:6px;background:#fff;"></iframe>')


def _zeige(ausgabe, inhalt):
    with ausgabe:
        display(HTML(inhalt))


def _fehler_anzeigen(ausgabe, fehler):
    with ausgabe:
        display(HTML(_hinweis(f"<b>Fehler:</b> {_e(fehler)}", "fehler")))
        if not isinstance(fehler, ValueError):        # unerwartete Fehler: Details für die Fehlersuche
            print(traceback.format_exc())


def _lese_eingabe(text):
    """Liest eine Sequenz aus einem Textfeld: roh, FASTA oder GenBank-Stil mit Zahlen.
    Gibt (name, dna, war_rna, hinweise) zurück; U wird zu T."""
    text = text.strip()
    if not text:
        raise ValueError("Das Feld ist leer. Füge eine Sequenz ein oder wähle ein Beispiel.")
    name, hinweise = "Eingabe", []
    if text.startswith(">") and "\n" not in text:
        # aus einem einzeiligen Feld eingefügt: Kopfzeile und Sequenz stehen in einer Zeile
        treffer = re.match(r"^>(.*?)\s((?:[ACGTUNacgtun0-9]+\s*)+)$", text)
        if treffer:
            text = ">" + treffer.group(1) + "\n" + treffer.group(2)
    if text.startswith(">"):
        eintraege = [t for t in text.split(">") if t.strip()]
        kopf, _, rumpf = eintraege[0].partition("\n")
        kopf = kopf.strip()
        if re.match(r"^[A-Z]{1,2}_?\d+(\.\d+)?\s", kopf + " "):   # NCBI-Kopfzeile, z. B. "NM_000518.5 Homo ..."
            name = wz.beschriftung(kopf)
        elif kopf:
            name = kopf[:70]
        text = rumpf
        if len(eintraege) > 1:
            hinweise.append(f"{len(eintraege)} FASTA-Einträge gefunden, verwendet wird der erste.")
    seq = re.sub(r"[^A-Za-z]", "", text).upper()
    if not seq:
        raise ValueError("Keine Basen gefunden.")
    war_rna = "U" in seq and "T" not in seq
    seq = seq.replace("U", "T")
    fremd = sorted(set(seq) - set("ACGTN"))
    if fremd:
        if len(re.sub(r"[ACGTN]", "", seq)) > 0.1 * len(seq):
            raise ValueError(f"Das sieht nicht nach DNA/RNA aus (Zeichen: {' '.join(fremd[:12])}). "
                             f"Proteinsequenzen werden hier nicht unterstützt.")
        seq = re.sub(r"[^ACGT]", "N", seq)
        hinweise.append(f"Mehrdeutige Zeichen ({' '.join(fremd)}) wurden als N gelesen.")
    return name, seq, war_rna, hinweise


def _erste_sequenz(datei, name_enthaelt=None):
    seqs = wz.lese_fasta(datei)
    for name, seq in seqs.items():
        if name_enthaelt is None or name_enthaelt in name:
            return name, seq
    raise ValueError(f"{name_enthaelt} nicht in {datei} gefunden")


def _als_fasta(name, seq):
    return f">{name}\n{_umbrechen(seq, 70, len(seq))}\n"


def _hochgeladene_dateien(uploader):
    """Dateien aus dem Upload-Feld, für ipywidgets 7 (dict) und 8 (tuple)."""
    wert = uploader.value
    if isinstance(wert, dict):
        return [(name, bytes(info["content"])) for name, info in wert.items()]
    return [(d["name"], bytes(d["content"])) for d in wert]


def _projektname(text):
    name = re.sub(r"[^A-Za-z0-9_\-äöüÄÖÜß]+", "_", text.strip()).strip("_")
    return name or "mein_projekt"


def _in_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def _ausfuehren(ausgabe, funktion):
    """Führt eine Werkzeug-Funktion aus und zeigt ihr HTML (oder den Fehler) im Ausgabefeld."""
    ausgabe.clear_output()
    try:
        _zeige(ausgabe, funktion())
    except Exception as fehler:
        _fehler_anzeigen(ausgabe, fehler)


# ======================================================================
# Die Werkzeuge als einfache Funktionen (geben HTML zurück).
# Die Reiter der Oberfläche und die Formular-Zellen im Notebook rufen genau diese auf.
# ======================================================================

def sequenz_html(text):
    """Sequenz umwandeln: Länge, GC, Gen, mRNA, Gegenstrang, Protein. Gibt HTML zurück."""
    name, dna, war_rna, hinweise = _lese_eingabe(text)
    z = wz.zusammensetzung(dna)
    mrna = wz.transkribiere(dna)
    gegen = wz.reverse_complement(dna)
    protein_r1 = wz.uebersetze(dna, bis_stopp=False)
    cds = wz.finde_cds(dna, min_as=10)
    teile = [f"<h3 style='margin:4px 0'>{_e(name)}</h3>"]
    teile += [_hinweis(_e(h), "warnung") for h in hinweise]
    if war_rna:
        teile.append(_hinweis("RNA erkannt (U statt T). Für die Werkzeuge wird sie als DNA (T) gelesen.",
                              "info"))
    teile.append(_tabelle(["Länge", "GC-Gehalt", "A", "C", "G", "T", "N"],
                          [[f"{len(dna):,} Basen", f"{wz.gc_gehalt(dna):.1f} %"] +
                           [f"{z.get(b, 0):,}" for b in "ACGTN"]], rechts=range(7)))
    if cds:
        s, e = cds
        gen = dna[s - 1:e]
        protein = wz.uebersetze(gen)
        m = wz.transkribiere(gen)
        n = min(20, len(gen) // 3)
        kette = ("DNA (codierend)  5'-" + gen[:n * 3] + ("…" if len(gen) > n * 3 else "") + "\n" +
                 "mRNA (T → U)     5'-" + m[:n * 3] + ("…" if len(m) > n * 3 else "") + "\n" +
                 "Codons              " + "".join(m[i:i + 3] + " " for i in range(0, n * 3, 3)) + "\n" +
                 "Protein             " + "".join(f" {a}  " for a in protein[:n]))
        teile.append(f"<b>Gen gefunden</b> (längster ORF ab ATG auf dem + Strang): Position {s}–{e}, "
                     f"{len(protein)} Aminosäuren")
        teile.append(_pre(kette))
        teile.append(f"<b>Protein ({len(protein)} AS)</b>" + _pre(_umbrechen(protein)))
    else:
        teile.append(_hinweis("Kein ORF ab 10 Aminosäuren auf dem + Strang gefunden. Unten steht die "
                              "Übersetzung ab Position 1.", "info"))
    teile.append(f"<b>mRNA</b> (Transkription, T → U)" + _pre(_umbrechen(mrna)))
    teile.append(f"<b>Gegenstrang</b> (Reverse Complement, 5'→3')" + _pre(_umbrechen(gegen)))
    teile.append(f"<b>Übersetzung ab Position 1</b> (Leseraster +1, * = Stopp)" +
                 _pre(_umbrechen(protein_r1)))
    downloads = (_text_download(_als_fasta(name + " mRNA", mrna), "mrna.fasta", "mRNA (FASTA)") +
                 _text_download(_als_fasta(name + " Gegenstrang", gegen), "gegenstrang.fasta",
                                "Gegenstrang (FASTA)"))
    if cds:
        downloads += _text_download(_als_fasta(name + " Protein", wz.uebersetze(dna[cds[0] - 1:cds[1]])),
                                    "protein.fasta", "Protein (FASTA)")
    teile.append(f"<div>{downloads}</div>")
    return _kasten("".join(teile))


def orfs_html(text, min_as=100, startcodons=("ATG",)):
    """ORFs in allen 6 Leserastern, mit Karte und Tabelle. Gibt HTML zurück."""
    name, dna, _r, hinweise = _lese_eingabe(text)
    orfs = wz.finde_orfs(dna, startcodons, min_as=min_as)
    teile = [f"<h3 style='margin:4px 0'>ORFs in {_e(name)}</h3>"]
    teile += [_hinweis(_e(h), "warnung") for h in hinweise]
    if not orfs:
        teile.append(_hinweis(f"Kein ORF ab {min_as} Aminosäuren gefunden. "
                              f"Schieb den Regler nach links.", "info"))
        return _kasten("".join(teile))
    from bericht import orf_karte_abbildung
    fig = orf_karte_abbildung({"sequenzen": [{"name": name, "laenge": len(dna), "orfs": orfs}]})
    teile.append(f"{len(orfs)} ORF(s) ab {min_as} Aminosäuren, Startcodons "
                 f"{', '.join(startcodons)}. Nicht jeder ORF ist ein echtes Gen.")
    teile.append(_bild(fig))
    zeilen = []
    for nr, o in enumerate(orfs[:200], 1):
        pos = f"{o['start']}–{o['ende']}" if o["strang"] == "+" else f"c{o['start']}–{o['ende']}"
        zeilen.append([nr, o["strang"], o["rahmen"], pos, f"{len(o['dna']):,}", len(o["protein"]),
                       f"<code>{_e(o['protein'][:40])}{'…' if len(o['protein']) > 40 else ''}</code>"])
    teile.append(_tabelle(["Nr", "Strang", "Rahmen", "Position", "bp", "AS", "Protein (Anfang)"], zeilen,
                          rechts=(0, 4, 5)))
    if len(orfs) > 200:
        teile.append(f"… und {len(orfs) - 200} weitere (alle im Download).")
    fasta = "".join(_als_fasta(f"ORF{nr} {o['strang']}{o['rahmen']} {o['start']}-{o['ende']} "
                               f"{len(o['protein'])}AS", o["protein"]) for nr, o in enumerate(orfs, 1))
    teile.append(f'<div>{_text_download(fasta, "orf_proteine.fasta", "Alle ORF-Proteine (FASTA)")}</div>')
    return _kasten("".join(teile))


def vergleich_html(text1, text2, modus="global", gegenstrang=True):
    """Zwei Sequenzen alignieren, Unterschiede und Wirkung. Gibt HTML zurück."""
    n1, a, _r1, h1 = _lese_eingabe(text1)
    n2, b, _r2, h2 = _lese_eingabe(text2)
    if len(a) * len(b) > MAX_ZELLEN_VERGLEICH:
        raise ValueError(f"Zu lang für ein Alignment ({len(a):,} × {len(b):,} Basen). "
                         f"Nimm kürzere Abschnitte oder die Analyse im Reiter 📊.")
    teile = [_hinweis(_e(h), "warnung") for h in h1 + h2]
    if gegenstrang:
        k1 = wz.kmer_menge(a, 8)
        gegen = wz.reverse_complement(b)
        if len(k1 & wz.kmer_menge(gegen, 8)) > len(k1 & wz.kmer_menge(b, 8)):
            b = gegen
            teile.append(_hinweis("Sequenz 2 passt als <b>Gegenstrang</b> besser und wurde umgedreht "
                                  "(Reverse Complement).", "info"))
    erg = wz.ausrichten(a, b, modus)
    st = wz.alignment_statistik(erg["a1"], erg["a2"])
    teile.append(f"<h3 style='margin:4px 0'>{_e(n1)} ↔ {_e(n2)}</h3>")
    teile.append(_tabelle(
        ["Modus", "Länge 1", "Länge 2", "Score", "Identität", "Mismatches", "Lücken"],
        [[modus, f"{len(a):,}", f"{len(b):,}", erg["score"], f"{st['identitaet']:.1f} %",
          st["mismatches"], st["luecken"]]], rechts=range(1, 7)))
    teile.append("<div><small>Identität im überlappenden Bereich. Bewertung wie NCBI BLASTN: Match +2, "
                 "Mismatch −3, Lücke 5 + 2 pro Base.</small></div>")
    if modus == "global" and min(len(a), len(b)) < 0.5 * max(len(a), len(b)):
        teile.append(_hinweis("Die Sequenzen sind sehr unterschiedlich lang. Ein globales Alignment "
                              "erzwingt große Lücken. Für ein kurzes Stück in einer langen Sequenz ist "
                              "<b>Lokal</b> besser.", "warnung"))
    if modus == "lokal":
        teile.append(f"Bester gemeinsamer Abschnitt: Sequenz 1 ab Position {erg['start1']}, "
                     f"Sequenz 2 ab Position {erg['start2']}.")
    if modus == "global":
        teile.append(_varianten_html(a, erg))
    text = wz.alignment_als_text(erg)
    zeilen = text.split("\n")
    if len(zeilen) > 160:
        text = "\n".join(zeilen[:160]) + f"\n… gekürzt ({len(zeilen) // 4} Blöcke insgesamt, " \
                                          f"vollständig im Download)"
    teile.append("<b>Alignment</b> (| gleich, . verschieden, Leerzeichen = Lücke)" + _pre(text))
    teile.append(f'<div>{_text_download(wz.alignment_als_text(erg), "alignment.txt", "Alignment (Text)")}</div>')
    return _kasten("".join(teile))


def mutation_html(art="snv", position=20, neue_base="T", anzahl=1, einfuegen="A", sequenz_text=None):
    """Mutation in ein Gen einbauen (Standard: HBB des Menschen). art = "snv", "del" oder "ins".
    position = c.-Position (c.1 = A des ATG). Gibt HTML zurück."""
    if sequenz_text is None:
        name, seq = _erste_sequenz(os.path.join(BEISPIELE, "01_sichelzellanaemie", "1_HbA_normal_NM_000518.5.fasta"))
    else:
        name, seq, _r, _h = _lese_eingabe(sequenz_text)
    cds = wz.finde_cds(seq, min_as=10)
    if not cds:
        raise ValueError("In dieser Sequenz wurde kein Gen (ORF ab ATG, mind. 10 AS) gefunden.")
    s, e = cds
    laenge = e - s + 1
    c = position
    if c > laenge:
        raise ValueError(f"Das Gen ist nur {laenge} Basen lang (c.1 bis c.{laenge}).")
    p = s - 1 + c                                 # Position in der ganzen Sequenz (ab 1)
    art = art
    if art == "snv":
        ref, alt = seq[p - 1], neue_base
        if ref == alt:
            raise ValueError(f"An c.{c} steht schon {ref}. Wähle eine andere Base.")
        v = {"art": "SNV", "start": p, "ende": p, "ref": ref, "alt": alt}
        mutiert = seq[:p - 1] + alt + seq[p:]
    elif art == "del":
        n = anzahl
        v = {"art": "Deletion", "start": p, "ende": p + n - 1, "ref": seq[p - 1:p - 1 + n], "alt": ""}
        mutiert = seq[:p - 1] + seq[p - 1 + n:]
    else:
        ins = re.sub(r"[^ACGT]", "", einfuegen.upper().replace("U", "T"))
        if not ins:
            raise ValueError("Bitte Basen zum Einfügen angeben (A, C, G, T).")
        v = {"art": "Insertion", "start": p, "ende": p + 1, "ref": "", "alt": ins}
        mutiert = seq[:p] + ins + seq[p:]
    original = dict(v)
    v = wz.normalisiere_indels(seq, [v])[0]
    wz.wirkung_der_varianten(seq, [v], cds)
    p_ref = wz.uebersetze(seq[s - 1:e])
    p_mut = wz.uebersetze(mutiert[s - 1:])
    teile = [f"<h3 style='margin:4px 0'>{_e(name)}: <code>{_e(v['c'])}</code> "
             f"{('<code>' + _e(v['p']) + '</code>') if v.get('p') else ''}</h3>"]
    kategorie = v["wirkung"].split(" (")[0]
    farbe = "fehler" if kategorie in ("Nonsense", "Frameshift", "Startcodon verändert") else \
        "ok" if kategorie == "stumm" else "warnung"
    teile.append(_hinweis(f"<b>Wirkung:</b> {_e(v['wirkung'])}", farbe))
    if original["start"] != v["start"]:
        teile.append(_hinweis("Hinweis: Nach der HGVS-Regel wird eine Lücke in einer Wiederholung "
                              "möglichst weit Richtung 3'-Ende benannt, deshalb weicht die Position "
                              "von deiner Eingabe ab. Das Ergebnis ist dieselbe Sequenz.", "info"))
    if art == "snv":
        k = (c - 1) // 3
        codon_alt = seq[s - 1 + 3 * k:s - 1 + 3 * k + 3]
        codon_neu = mutiert[s - 1 + 3 * k:s - 1 + 3 * k + 3]
        teile.append(_tabelle(["Codon Nr.", "DNA vorher", "DNA nachher", "mRNA vorher", "mRNA nachher",
                               "Aminosäure vorher", "nachher"],
                              [[k + 1, codon_alt, codon_neu, wz.transkribiere(codon_alt),
                                wz.transkribiere(codon_neu),
                                wz.DREIBUCHSTABEN.get(wz.CODON_TABELLE.get(codon_alt, "X"), "?"),
                                wz.DREIBUCHSTABEN.get(wz.CODON_TABELLE.get(codon_neu, "X"), "?")]]))
    teile.append(f"<b>Protein:</b> {_e(wz.vergleiche_proteine(p_ref, p_mut))}")
    teile.append(_protein_vergleich(p_ref, p_mut))
    teile.append("<div>" + _text_download(_als_fasta(f"{name} {v['c']}", mutiert), "mutiert.fasta",
                                          "mutierte Sequenz (FASTA)") + "</div>")
    return _kasten("".join(teile))


# ----------------------------------------------------------------------
# Zwei Proteine vergleichen (BLOSUM62)
# ----------------------------------------------------------------------

MAX_ZELLEN_PROTEIN = 25_000_000                       # Länge1 × Länge2 in Aminosäuren
_AS_BUCHSTABEN = set("ACDEFGHIKLMNPQRSTVWYBZX")


def _lese_protein(text, min_as=30):
    """Liest ein Protein aus einem Textfeld. DNA/mRNA wird erkannt und das Gen (längster ORF ab ATG)
    übersetzt. Gibt (name, protein, quelle, hinweise) zurück."""
    text = text.strip()
    if not text:
        raise ValueError("Das Feld ist leer. Füge eine Sequenz ein oder wähle ein Beispiel.")
    name, hinweise = "Eingabe", []
    if text.startswith(">") and "\n" not in text:
        # einzeiliges Feld: Kopfzeile und Sequenz in einer Zeile. Von hinten alle Wörter aus
        # Großbuchstaben oder Zahlen nehmen, die gehören zur Sequenz.
        woerter = text.split()
        k = len(woerter)
        while k > 1 and re.fullmatch(r"[A-Z*\-]+|\d+", woerter[k - 1]):
            k -= 1
        text = " ".join(woerter[:k]) + "\n" + " ".join(woerter[k:])
    if text.startswith(">"):
        eintraege = [t for t in text.split(">") if t.strip()]
        kopf, _, rumpf = eintraege[0].partition("\n")
        kopf = kopf.strip()
        if re.match(r"^[A-Z]{1,2}_?\d+(\.\d+)?\s", kopf + " "):
            name = wz.beschriftung(kopf)
        elif kopf:
            name = kopf[:70]
        text = rumpf
        if len(eintraege) > 1:
            hinweise.append(f"{len(eintraege)} FASTA-Einträge gefunden, verwendet wird der erste.")
    seq = re.sub(r"[^A-Za-z*]", "", text).upper()
    if not seq:
        raise ValueError("Keine Sequenz gefunden.")
    if not wz.ist_protein(seq):
        dna = seq.replace("U", "T")
        protein, cds = wz.protein_aus_dna(dna, min_as=min_as)
        if not protein:
            raise ValueError(f"„{name}“ ist DNA/RNA, aber darin steht kein Gen (ORF ab ATG mit mindestens "
                             f"{min_as} Aminosäuren auf dem angegebenen Strang). Füge die mRNA/CDS ein "
                             f"oder direkt die Proteinsequenz.")
        quelle = f"übersetzt aus DNA, Gen an Position {cds[0]}–{cds[1]}"
        return name, protein.rstrip("*"), quelle, hinweise
    protein = seq.rstrip("*")
    fremd = sorted(set(protein) - _AS_BUCHSTABEN - {"*"})
    if fremd:
        if len([c for c in protein if c in fremd]) > 0.05 * len(protein):
            raise ValueError(f"Unbekannte Zeichen in „{name}“: {' '.join(fremd)}. "
                             f"Erlaubt sind die 20 Aminosäure-Buchstaben (und B, Z, X).")
        hinweise.append(f"Ungewöhnliche Zeichen ({' '.join(fremd)}) in „{name}“ wurden als X gelesen.")
        protein = "".join("X" if c in fremd else c for c in protein)
    if "*" in protein:
        hinweise.append(f"„{name}“ enthält ein Stoppzeichen (*) mitten in der Sequenz.")
    return name, protein, "Proteinsequenz (eingegeben)", hinweise


def _blosum_marke(punkte):
    """BLOSUM62-Wert farbig: grün = häufiger Austausch, rot = seltener, grau = 0."""
    farbe = "#009E73" if punkte > 0 else ("#D55E00" if punkte < 0 else "#6b6b6b")
    return f'<span style="color:{farbe};font-weight:600">{punkte:+d}</span>' if punkte else \
        f'<span style="color:{farbe};font-weight:600">0</span>'


def _kennzahlen_html(liste):
    """Tabelle der Protein-Kennzahlen und Diagramm der Zusammensetzung für [(name, protein), ...]."""
    import bericht
    werte = [(name, wz.protein_kennzahlen(p)) for name, p in liste]
    kopf, zeilen = bericht.kennzahlen_tabelle(werte)
    teile = ["<h3 style='margin:12px 0 4px'>Eigenschaften (berechnet aus der Sequenz)</h3>"]
    teile.append(_tabelle(kopf, [[_e(z) for z in zeile] for zeile in zeilen]))
    if len(werte) == 2:
        (n1, k1), (n2, k2) = werte
        d = k2["ladung_ph74"] - k1["ladung_ph74"]
        if abs(d) >= 0.5:
            teile.append(_hinweis(
                f"Die Ladung bei pH 7,4 unterscheidet sich um {d:+.1f}. Solche Proteine wandern in der "
                f"Gel-Elektrophorese unterschiedlich weit. Genau so wurde 1949 entdeckt, dass sich das "
                f"Sichelzell-Hämoglobin (HbS) vom normalen (HbA) unterscheidet.", "info"))
        if k1["stabil"] != k2["stabil"]:
            teile.append(_hinweis("Der Instabilitätsindex stuft die beiden Proteine verschieden ein. Das ist nur "
                                  "ein grober Hinweis, keine Messung.", "warnung"))
    if any(k["extinktion"] == 0 for _n, k in werte):
        teile.append(_hinweis("Ein Protein ohne Trp und Tyr absorbiert bei 280 nm kaum. Seine Konzentration lässt "
                              "sich so nicht messen.", "warnung"))
    teile.append(_bild(bericht.zusammensetzung_abbildung(werte)))
    teile.append("<div><small>Wie ExPASy ProtParam (pK-Werte nach Bjellqvist, Hydropathie nach Kyte-Doolittle, "
                 "Instabilitätsindex nach Guruprasad 1990), gegen Biopython geprüft. Gilt für die nackte Kette, so wie "
                 "sie eingegeben wurde: ohne Häm, Zucker oder Phosphatgruppen und mit Start-Methionin, das in der "
                 "Zelle oft abgespalten wird. Echte Proteine bestehen oft aus mehreren Ketten (Hämoglobin: 2 α + 2 β)."
                 "</small></div>")
    return "".join(teile)


def protein_eigenschaften_html(text):
    """Kennzahlen eines einzelnen Proteins (Eingabe als Protein oder DNA/mRNA). Gibt HTML zurück."""
    name, protein, quelle, hinweise = _lese_protein(text)
    teile = [_hinweis(_e(h), "warnung") for h in hinweise]
    teile.append(f"<h3 style='margin:4px 0'>{_e(name)}</h3><div>{_e(quelle)} · {len(protein)} Aminosäuren</div>")
    teile.append(_kennzahlen_html([(name, protein)]))
    teile.append("<b>Protein</b>" + _pre(_umbrechen(protein, 60, len(protein))))
    teile.append("<div>" + _text_download(_als_fasta(name, protein), "protein.fasta", "Protein (FASTA)") + "</div>")
    return _kasten("".join(teile))


def protein_vergleich_html(text1, text2, modus="global"):
    """Zwei Proteine mit BLOSUM62 alignieren (wie BLASTP): Identität, Ähnlichkeit, jede Änderung mit
    chemischer Bewertung. Eingabe als Protein oder als DNA/mRNA (dann wird das Gen übersetzt)."""
    from bericht import DREI
    n1, p1, q1, h1 = _lese_protein(text1)
    n2, p2, q2, h2 = _lese_protein(text2)
    if n1 == n2:
        n1, n2 = n1 + " (1)", n2 + " (2)"
    if len(p1) * len(p2) > MAX_ZELLEN_PROTEIN:
        raise ValueError(f"Zu lang ({len(p1):,} × {len(p2):,} Aminosäuren). Nimm kürzere Proteine.")
    erg = wz.ausrichten_protein(p1, p2, modus)
    st = wz.protein_statistik(erg["a1"], erg["a2"])
    teile = [_hinweis(_e(h), "warnung") for h in h1 + h2]
    teile.append(f"<h3 style='margin:4px 0'>{_e(n1)} ↔ {_e(n2)}</h3>")
    teile.append(_tabelle(["", "Name", "Länge (AS)", "Herkunft"],
                          [["1 (Referenz)", _e(n1), len(p1), _e(q1)], ["2", _e(n2), len(p2), _e(q2)]],
                          rechts=(2,)))
    teile.append(_tabelle(
        ["Modus", "Score", "Identität", "Ähnlichkeit", "Lücken (Spalten)"],
        [[modus, int(erg["score"]), f"{st['identitaet']:.1f} %", f"{st['aehnlichkeit']:.1f} %", st["luecken"]]],
        rechts=(1, 2, 3, 4)))
    teile.append("<div><small>Bewertung wie NCBI BLASTP: BLOSUM62, Lücke öffnen 11, verlängern 1. "
                 "<b>Identität</b> = gleiche Aminosäure. <b>Ähnlichkeit</b> = gleiche oder chemisch ähnliche "
                 "Aminosäure (BLOSUM62 &gt; 0). Werte beziehen sich auf den überlappenden Bereich.</small></div>")
    if st["identitaet"] < 25 and st["laenge"] > 50:
        teile.append(_hinweis("Unter etwa 25 % Identität ist unsicher, ob die Proteine überhaupt verwandt sind "
                              "(„Twilight Zone“). Das Alignment kann dann zufällig sein.", "warnung"))
    if modus == "global" and min(len(p1), len(p2)) < 0.5 * max(len(p1), len(p2)):
        teile.append(_hinweis("Die Proteine sind sehr unterschiedlich lang. Für eine Domäne in einem großen "
                              "Protein ist <b>lokal</b> besser.", "warnung"))
    if modus == "lokal":
        teile.append(f"Bester gemeinsamer Abschnitt: Protein 1 ab Aminosäure {erg['start1']}, "
                     f"Protein 2 ab Aminosäure {erg['start2']}.")

    teile.append(_kennzahlen_html([(n1, p1), (n2, p2)]))
    teile.append("<h3 style='margin:12px 0 4px'>Unterschiede</h3>")
    unterschiede = wz.protein_unterschiede(erg["a1"], erg["a2"])
    versatz = erg.get("start1", 1) - 1                 # lokal: Positionen in Protein 1 zählen ab dessen Anfang
    if not unterschiede:
        teile.append(_hinweis("Im verglichenen Bereich sind die Proteine identisch.", "ok"))
    else:
        zaehl = {}
        for d in unterschiede:
            schluessel = d["art"] if d["art"] != "Austausch" else (
                "Austausch, ähnlich (BLOSUM62 > 0)" if d["punkte"] > 0 else "Austausch, nicht ähnlich")
            zaehl[schluessel] = zaehl.get(schluessel, 0) + 1
        teile.append("<p><b>" + str(len(unterschiede)) + " Unterschied(e):</b> " +
                     " · ".join(f"{_e(k)}: {n}" for k, n in zaehl.items()) + "</p>")
        zeilen = []
        for d in unterschiede[:200]:
            pos, ende = d["position"] + versatz, d["ende"] + versatz
            if d["art"] == "Austausch":
                text = f"{DREI.get(d['ref'], d['ref'])}{pos}{DREI.get(d['alt'], d['alt'])}"
                punkte = _blosum_marke(d["punkte"])
            elif d["art"] == "Insertion":
                text, punkte = f"+{d['alt']} (nach {pos})", ""
            else:
                text = f"−{d['ref']} ({pos if pos == ende else f'{pos}–{ende}'})"
                punkte = ""
            zeilen.append([d["art"], f'<code>{_e(text)}</code>', punkte, _e(d["bewertung"])])
        teile.append(_tabelle(["Art", "Änderung (Position in Protein 1)", "BLOSUM62", "Bewertung"], zeilen,
                              rechts=(2,)))
        if len(unterschiede) > 200:
            teile.append(f"<small>… und {len(unterschiede) - 200} weitere (im Alignment unten).</small>")
        teile.append("<div><small>Gruppen (vereinfacht): unpolar G A V L I M F W P · polar S T C Y N Q · "
                     "basisch (+) K R H · sauer (−) D E. Positiver BLOSUM62-Wert = Austausch kommt in echten "
                     "Proteinfamilien oft vor (meist harmlos), negativ = selten (eher störend).</small></div>")
    text = wz.protein_alignment_text(erg)
    zeilen_txt = text.split("\n")
    if len(zeilen_txt) > 120:
        text = "\n".join(zeilen_txt[:120]) + "\n… gekürzt, vollständig im Download"
    teile.append("<b>Alignment</b> (| gleich, + ähnlich, Leerzeichen = verschieden oder Lücke)" + _pre(text))
    teile.append("<div>" + _text_download(wz.protein_alignment_text(erg), "protein_alignment.txt",
                                          "Alignment (Text)") + " " +
                 _text_download(_als_fasta(n1, p1) + _als_fasta(n2, p2), "proteine.fasta",
                                "beide Proteine (FASTA)") + "</div>")
    return _kasten("".join(teile))


# ----------------------------------------------------------------------
# 3D-Struktur: Bindestelle, Lage von Mutationen, 3D-Ansicht
# ----------------------------------------------------------------------

DREIDMOL_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"
STRUKTUR_ORDNER = os.path.abspath("strukturen")


def _js_auswahl(S, indizes):
    """3Dmol-Auswahl {chain, resi} für eine Liste von Rest-Indizes, nach Kette gruppiert."""
    nach_kette = {}
    for j in indizes:
        r = S["reste"][j]
        nach_kette.setdefault(r["kette"], []).append(r["nr"])
    return [{"chain": k if k != "_" else "", "resi": sorted(set(v))} for k, v in nach_kette.items()]


def _ansicht_3d(S, pdb_text, kette=None, bindestelle=(), positionen=(), merkmale=(), hoehe=520):
    """HTML mit 3Dmol.js: Cartoon, Liganden als Stäbchen, Bindestelle hellblau, Positionen rot."""
    import json as _json
    import uuid
    zeilen = [z for z in pdb_text.splitlines()
              if z.startswith(("ATOM", "HETATM", "TER", "HELIX", "SHEET", "MODEL", "ENDMDL"))
              and z[17:20].strip() not in ("HOH", "WAT", "DOD")]
    modell = []
    for z in zeilen:                         # nur das erste Modell
        modell.append(z)
        if z.startswith("ENDMDL"):
            break
    kennung = "v" + uuid.uuid4().hex[:10]
    labels = []
    for j, text in positionen:
        r = S["reste"][j]
        labels.append({"text": text, "sel": {"chain": r["kette"], "resi": r["nr"], "atom": "CA"}})
    befehle = []
    if S["alphafold"]:
        befehle.append("v.setStyle({}, {cartoon: {colorfunc: function(a) { return a.b > 90 ? '#0053D6' : "
                       "a.b > 70 ? '#65CBF3' : a.b > 50 ? '#FFDB13' : '#FF7D45'; }}});")
    else:
        befehle.append("v.setStyle({}, {cartoon: {color: '#d0d0d0'}});")
        if kette:
            befehle.append(f"v.setStyle({{chain: {_json.dumps(kette)}}}, {{cartoon: {{color: 'spectrum'}}}});")
    ligand_namen = sorted({S["reste"][j]["name"] for j in S["liganden"]})
    if ligand_namen:
        befehle.append(f"v.addStyle({{resn: {_json.dumps(ligand_namen)}}}, "
                       "{stick: {colorscheme: 'greenCarbon', radius: 0.25}});")
    for sel in _js_auswahl(S, bindestelle):
        befehle.append(f"v.addStyle({_json.dumps(sel)}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.18}}}});")
    for sel in _js_auswahl(S, [j for j, _t in merkmale]):
        befehle.append(f"v.addStyle({_json.dumps(sel)}, {{stick: {{colorscheme: 'orangeCarbon', radius: 0.22}}}});")
    for sel in _js_auswahl(S, [j for j, _t in positionen]):
        befehle.append(f"v.addStyle({_json.dumps(sel)}, {{sphere: {{color: '#D55E00', scale: 0.75}}}});")
    befehle.append("L.forEach(function(l) { v.addLabel(l.text, {fontSize: 13, fontColor: 'white', "
                   "backgroundColor: '#D55E00', backgroundOpacity: 0.85, inFront: true}, l.sel); });")
    zoom = _json.dumps({"chain": kette}) if kette else "{}"
    skript = f"""
<div id="{kennung}" style="width:100%;height:{hoehe}px;position:relative;border:1px solid #ddd;border-radius:6px;
background:#fff;"><div style="padding:12px;color:#666;font-family:sans-serif;">3D-Ansicht wird geladen …</div></div>
<script>
(function() {{
  var pdb = {_json.dumps(chr(10).join(modell))};
  var L = {_json.dumps(labels)};
  function zeichne() {{
    var el = document.getElementById("{kennung}");
    el.innerHTML = "";
    var v = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
    v.addModel(pdb, "pdb");
    {chr(10).join("    " + b for b in befehle)}
    v.zoomTo({zoom});
    v.render();
  }}
  if (window.$3Dmol) {{ zeichne(); return; }}
  var s = document.createElement("script");
  s.src = "{DREIDMOL_URL}";
  s.onload = zeichne;
  s.onerror = function() {{
    document.getElementById("{kennung}").innerHTML = "<div style='padding:12px'>Die 3D-Ansicht konnte nicht " +
      "geladen werden (keine Verbindung zu cdn.jsdelivr.net). Die Tabellen oben gelten trotzdem.</div>";
  }};
  document.head.appendChild(s);
}})();
</script>"""
    return skript


def _zugaenglichkeit_abbildung(S, kette, bindestelle=(), positionen=()):
    """Balken: relative Zugänglichkeit jedes Restes der Kette; blau = berührt Ligand, rot = gewählte Position."""
    from matplotlib.figure import Figure
    reste = S["ketten"][kette]["reste"]
    werte = [S["reste"][j].get("relativ") or 0 for j in reste]
    nummern = [S["reste"][j]["nr"] for j in reste]
    farben = ["#0072B2" if j in bindestelle else "#bdbdbd" for j in reste]
    hoehe = 3.6 if S["alphafold"] else 2.8
    fig = Figure(figsize=(11, hoehe))
    if S["alphafold"]:
        ax, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    else:
        ax = fig.add_subplot()
    ax.bar(nummern, werte, width=1.0, color=farben)
    ax.axhline(20, color="#888", lw=0.8, ls="--")
    ax.text(nummern[0], 21, "unter 20 % = vergraben", fontsize=7, color="#666")
    for j, text in positionen:
        r = S["reste"][j]
        if r["kette"] == kette:
            ax.bar([r["nr"]], [max(r.get("relativ") or 0, 3)], width=1.6, color="#D55E00")
            ax.annotate(text, (r["nr"], max(r.get("relativ") or 0, 3)), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=8, color="#D55E00")
    ax.set_ylabel("frei zugänglich (%)", fontsize=8)
    ax.set_ylim(0, 110)
    ax.set_title(f"Kette {kette}: Oberfläche je Aminosäure (blau = berührt einen Liganden, rot = gewählte Position)",
                 fontsize=9, loc="left")
    for a in [ax] + ([ax2] if S["alphafold"] else []):
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=8)
    if S["alphafold"]:
        plddt = [float(S["b"][S["reste"][j]["atome"]].mean()) for j in reste]
        ax2.plot(nummern, plddt, color="#0053D6", lw=1.2)
        ax2.axhline(70, color="#888", lw=0.8, ls="--")
        ax2.set_ylim(0, 100)
        ax2.set_ylabel("pLDDT", fontsize=8)
        ax2.set_xlabel("Nummer in der Struktur", fontsize=8)
    else:
        ax.set_xlabel("Nummer in der Struktur", fontsize=8)
    fig.tight_layout()
    return fig


def struktur_html(struktur, protein_text=None, positionen="", uniprot="", grenze=4.0, ordner=None):
    """3D-Struktur laden und auswerten. struktur = PDB-ID, UniProt-Nummer (AlphaFold) oder Dateipfad.
    protein_text = Protein/DNA, auf das sich die Positionen beziehen (None = Nummern der Struktur).
    positionen = z. B. "7" oder "E7V, 64". uniprot = UniProt-Nummer für bekannte aktive Stellen/Bindestellen."""
    import numpy as np
    import struktur as st
    text, quelle = st.lade_struktur(struktur, ordner or STRUKTUR_ORDNER)
    S = st.lese_pdb(text)
    st.zugaenglichkeit(S)
    teile, warnungen = [], []

    # --- Kopf
    kennung = st.kennung_aus_text(struktur) if not os.path.isfile(struktur) else os.path.basename(struktur)
    teile.append(f"<h3 style='margin:4px 0'>{_e(kennung)}: {_e(S['titel'] or 'ohne Titel')}</h3>")
    info = [_e(quelle)]
    if S["alphafold"]:
        info.append("Vorhersage (AlphaFold), nicht gemessen")
    elif S["methode"]:
        info.append(_e(S["methode"].title()))
    if S["aufloesung"]:
        info.append(f"Auflösung {S['aufloesung']:.2f} Å")
    teile.append("<div>" + " · ".join(info) + "</div>")
    if S["aufloesung"]:
        teile.append("<div><small>Auflösung: je kleiner, desto genauer. Unter etwa 2,5 Å sind einzelne Seitenketten "
                     "gut erkennbar, über 3,5 Å nur noch grob.</small></div>")
    zeilen = []
    gruppen = {}
    for k, inf in S["ketten"].items():
        gruppen.setdefault(inf["sequenz"], []).append(k)
    for seq, ks in gruppen.items():
        zeilen.append([", ".join(ks), len(seq), f'<span style="font-family:monospace;font-size:12px">'
                       f'{_e(seq[:40])}{"…" if len(seq) > 40 else ""}</span>'])
    teile.append(_tabelle(["Kette(n)", "Aminosäuren", "Sequenz (Anfang)"], zeilen, rechts=(1,)))
    liganden = {}
    for j in S["liganden"]:
        liganden.setdefault(S["reste"][j]["name"], []).append(j)
    if liganden:
        teile.append("<div><b>Liganden</b> (gebundene Moleküle, ohne Wasser): " + ", ".join(
            f"{_e(st.ligand_name(n))} × {len(v)}" for n, v in liganden.items()) + "</div>")
    elif S["alphafold"]:
        warnungen.append("AlphaFold-Modelle enthalten keine Liganden (kein Häm, kein Hemmstoff). Für die Bindestelle "
                         "eine gemessene Struktur aus der PDB nehmen, z. B. 4HHB für Hämoglobin.")
    else:
        warnungen.append("Diese Struktur enthält keine Liganden. Eine Bindestelle lässt sich hier nicht ablesen.")
    if S["alphafold"]:
        mittel = float(np.mean([S["b"][r["atome"]].mean() for r in S["reste"] if r["typ"] == "aminosaeure"]))
        teile.append(f"<div>AlphaFold-Sicherheit (pLDDT) im Mittel: <b>{mittel:.0f}</b> von 100 · Farben in 3D: "
                     "<span style='color:#0053D6'>■</span> über 90 sehr sicher, <span style='color:#65CBF3'>■</span> "
                     "70–90 sicher, <span style='color:#d4b000'>■</span> 50–70 unsicher, "
                     "<span style='color:#FF7D45'>■</span> unter 50 sehr unsicher (oft ungeordnete Bereiche)</div>")

    # --- Protein zuordnen
    zuordnung, kette, protein, pname = None, None, None, None
    uni = None
    if uniprot.strip():
        try:
            uni = st.uniprot_merkmale(uniprot)
        except Exception as fehler:          # Netz, unbekannte Nummer oder geändertes Format
            warnungen.append(f"UniProt {uniprot.strip()} konnte nicht gelesen werden ({type(fehler).__name__}). "
                             f"Der Rest der Auswertung ist davon nicht betroffen.")
    if protein_text and protein_text.strip():
        pname, protein, _q, hinweise = _lese_protein(protein_text)
        warnungen += hinweise
    if protein:
        z = st.ordne_protein_zu(S, protein)
        if z["identitaet"] < 30 or z["abgedeckt"] < 20:
            warnungen.append(f"„{pname}“ passt zu keiner Kette dieser Struktur (beste: Kette {z['kette']}, "
                             f"{z['identitaet']:.0f} % gleich). Die Positionen werden deshalb als Nummern der "
                             f"Struktur gelesen.")
            protein = None
        else:
            zuordnung, kette = z["zuordnung"], z["kette"]
            gleich = f" Dieselbe Sequenz haben auch Kette {', '.join(z['gleiche_ketten'])}." if z["gleiche_ketten"] else ""
            beispiel = next(iter(sorted(zuordnung.items())), None)
            versatz = ""
            if beispiel:
                pos, j = beispiel
                if S["reste"][j]["nr"] != pos:
                    versatz = (f" Achtung, die Nummern unterscheiden sich: Position {pos} im Protein ist "
                               f"{_e(st.rest_text(S['reste'][j]))} in der Struktur.")
            teile.append(_hinweis(f"„{_e(pname)}“ ({len(protein)} AS) entspricht <b>Kette {_e(kette)}</b>: "
                                  f"{z['identitaet']:.0f} % gleich über {z['abgedeckt']} Positionen.{_e(gleich)}"
                                  f"{versatz}", "info"))
    if kette is None:
        kette = next(iter(S["ketten"]), None)
    z_uni, uni_varianten = None, None
    if uni and uni["sequenz"]:
        z_uni = st.ordne_protein_zu(S, uni["sequenz"])
        if z_uni["identitaet"] < 30:
            warnungen.append(f"Die UniProt-Sequenz {uniprot.strip().upper()} passt nicht zu dieser Struktur "
                             f"({z_uni['identitaet']:.0f} % gleich). Prüfe die Nummer.")
            z_uni = None
        else:
            uni_varianten = {}
            for v in uni["varianten"]:
                j = z_uni["zuordnung"].get(v["position"])
                if j is not None:
                    uni_varianten.setdefault(j, []).append(v)

    # --- gewählte Positionen
    gewaehlt, ergebnisse = [], []
    if positionen.strip():
        for pos, alt, neu in st.lies_positionen(positionen):
            if zuordnung is not None:
                j = zuordnung.get(pos)
                if j is None:
                    warnungen.append(f"Position {pos} des Proteins ist in der Struktur nicht enthalten "
                                     f"(nicht aufgelöst oder außerhalb).")
                    continue
            else:
                treffer = [j for j in S["ketten"].get(kette, {"reste": []})["reste"] if S["reste"][j]["nr"] == pos]
                if not treffer:
                    warnungen.append(f"Kette {kette} hat keinen Rest Nummer {pos}.")
                    continue
                j = treffer[0]
            r = S["reste"][j]
            if alt and alt != r["buchstabe"]:
                warnungen.append(f"Position {pos}: Du hast {alt} angegeben, in der Struktur steht dort "
                                 f"{st.rest_text(r)}. Prüfe die Nummerierung.")
            b = st.bewerte_position(S, j, neu, grenze)
            b["eingabe"] = f"{alt or ''}{pos}{neu or ''}"
            b["index"] = j
            ergebnisse.append(b)
            gewaehlt.append((j, st.rest_text(r, False) + (f"→{neu}" if neu else "")))
    for w in warnungen:
        teile.append(_hinweis(_e(w), "warnung"))
    if ergebnisse:
        teile.append("<h3 style='margin:12px 0 4px'>Wo liegen die gewählten Positionen?</h3>")
        zeilen = []
        for b in ergebnisse:
            r = b["rest"]
            lig = (f"{_e(st.ligand_name(S['reste'][b['ligand'][0]]['name']))}: {b['ligand'][1]:.1f} Å"
                   if b["ligand"] else "–")
            zeile = [_e(b["eingabe"]), _e(st.rest_text(r)),
                     f"{b['relativ']:.0f} %" if b["relativ"] is not None else "–", _e(b["lage"]), lig,
                     _e(", ".join(b["andere_ketten"]) or "–"),
                     "<br>".join(_e(h) for h in b["hinweise"]) or "–"]
            if uni_varianten is not None:
                bekannt = uni_varianten.get(b["index"], [])
                zeile.append("<br>".join(f"<b>{_e(v['ref'])}→{_e(v['alt'])}</b>: {_e(v['text'])}"
                                         for v in bekannt[:8]) or "keine")
            zeilen.append(zeile)
        teile.append(_tabelle(["Eingabe", "In der Struktur", "Frei zugänglich", "Lage", "Nächster Ligand",
                               "Kontakt zu Kette", "Einschätzung"]
                              + (["Bekannte Varianten (UniProt)"] if uni_varianten is not None else []),
                              zeilen, rechts=(2,)))
        teile.append("<div><small><b>Frei zugänglich</b>: wie viel der Oberfläche dieser Aminosäure Wasser "
                     "berühren kann (Shrake-Rupley, mit allen Ketten und Liganden der Datei). Unter 20 % gilt als "
                     "vergraben, über 50 % als Oberfläche (übliche Faustregeln). Die Einschätzung ist ein Hinweis aus "
                     "Lage und Chemie, keine Vorhersage der Stabilität.</small></div>")

    # --- Bindestellen
    bindestelle = set()
    zeige = []
    for name, js in liganden.items():
        if name in st.ZUSAETZE:
            continue
        for j in js:
            k = st.kontakte(S, j, grenze)
            if any(S["reste"][a]["kette"] == kette for a, _d in k):
                zeige.append((j, k))
                break
        else:
            if js:
                zeige.append((js[0], st.kontakte(S, js[0], grenze)))
    if zeige:
        teile.append("<h3 style='margin:12px 0 4px'>Bindestellen: welche Aminosäuren berühren den Liganden?</h3>")
        teile.append("<div><small>Gezählt wird jede Aminosäure mit einem Atom höchstens "
                     f"{grenze:.1f} Å vom Liganden entfernt (Wasserstoffbrücken sind etwa 2,6–3,3 Å lang). "
                     "Bei Enzymen liegt hier das <b>aktive Zentrum</b>, bei Hämoglobin die Tasche, die das Häm "
                     "und damit den Sauerstoff hält.</small></div>")
        rueck = {j: p for p, j in (zuordnung or {}).items()}
        for j, k in zeige:
            L = S["reste"][j]
            bindestelle.update(a for a, _d in k)
            zeilen = []
            for a, d in k:
                r = S["reste"][a]
                zeile = [_e(st.rest_text(r)), f"{d:.1f} Å", _e(wz.EIGENSCHAFT.get(r["buchstabe"], "?")),
                         f"{r['relativ']:.0f} %" if r.get("relativ") is not None else "–"]
                if zuordnung:
                    zeile.insert(2, rueck.get(a, "–"))
                zeilen.append(zeile)
            anzahl = len(liganden[L["name"]])
            offen = " open" if len(zeige) <= 2 else ""
            teile.append(f"<details{offen}><summary><b>{_e(st.ligand_name(L['name']))}</b> an "
                         f"{_e(L['kette'])}:{_e(L['kurz'])}: {len(k)} Aminosäuren"
                         + (f" (dieser Ligand kommt {anzahl}× vor, gezeigt ist der an Kette {_e(kette)})"
                            if anzahl > 1 else "") + "</summary>"
                         + _tabelle(["Aminosäure", "Abstand"] + (["Position im Protein"] if zuordnung else [])
                                    + ["Gruppe", "Frei zugänglich"], zeilen,
                                    rechts=(1, 2, 4) if zuordnung else (1, 3)) + "</details>")

    # --- UniProt
    merkmale_3d = []
    if uni and uni["merkmale"]:
        zeilen = []
        for m in uni["merkmale"]:
            orte = []
            if z_uni:
                for p in ([m["start"], m["ende"]] if m["art"] == "Disulfidbrücke" else
                          range(m["start"], m["ende"] + 1)):
                    j = z_uni["zuordnung"].get(p)
                    if j is not None:
                        orte.append(st.rest_text(S["reste"][j]))
                        merkmale_3d.append((j, m["art"]))
            pos = f"{m['start']}" if m["start"] == m["ende"] else f"{m['start']}–{m['ende']}"
            zeilen.append([_e(m["art"]), pos, _e(m["text"] or "–"), _e(", ".join(orte[:6]) or "nicht in der Struktur")])
        teile.append(f"<h3 style='margin:12px 0 4px'>Bekannt aus UniProt {_e(uniprot.strip().upper())}"
                     f"{': ' + _e(uni['name']) if uni['name'] else ''}</h3>")
        teile.append(_tabelle(["Art", "Position (UniProt)", "Beschreibung", "In der Struktur"], zeilen))
        teile.append("<div><small>Diese Stellen wurden in Experimenten oder durch Vergleich mit verwandten Proteinen "
                     "bestimmt (UniProt-Annotation). In 3D orange.</small></div>")
    elif uni is not None:
        teile.append(_hinweis("UniProt nennt für dieses Protein keine aktiven Stellen oder Bindestellen.", "info"))

    # --- 3D und Diagramm
    teile.append("<h3 style='margin:12px 0 4px'>3D-Ansicht</h3><div><small>Maus: ziehen = drehen, Mausrad = zoomen, "
                 "Rechtsklick ziehen = verschieben. "
                 + ("Farbe = AlphaFold-Sicherheit. " if S["alphafold"] else
                    f"Kette {_e(kette)} in Regenbogenfarben (blau = Anfang, rot = Ende), andere Ketten grau. ")
                 + "Grüne Stäbchen = Liganden, hellblau = Bindestelle, orange = UniProt-Stellen, "
                 "rote Kugel = gewählte Position.</small></div>")
    teile.append(_ansicht_3d(S, text, None if S["alphafold"] else kette, sorted(bindestelle), gewaehlt, merkmale_3d))
    if kette:
        teile.append(_bild(_zugaenglichkeit_abbildung(S, kette, bindestelle, gewaehlt)))

    # --- Downloads
    csv = ["kette;nummer;aminosaeure;flaeche_A2;frei_prozent;lage;beruehrt_ligand"]
    for j, r in enumerate(S["reste"]):
        if r["typ"] == "aminosaeure":
            rel = f"{r['relativ']:.1f}" if r.get("relativ") is not None else ""
            csv.append(f"{r['kette']};{r['kurz']};{r['name']};{r['flaeche']:.1f};{rel};{st.lage(r.get('relativ'))};"
                       f"{'ja' if j in bindestelle else ''}")
    teile.append("<div>" + _text_download("\n".join(csv) + "\n", f"{kennung}_oberflaeche.csv",
                                          "Oberfläche je Aminosäure (CSV)") +
                 _text_download(text, f"{kennung}.pdb", "Struktur (PDB-Datei)") + "</div>")
    return _kasten("".join(teile))


BEISPIEL_SEQUENZEN = {
    "kurzes Testgen": ("kurz", None),
    "HBB-mRNA des Menschen": ("01_sichelzellanaemie/1_HbA_normal_NM_000518.5.fasta", None),
    "HbS (Sichelzelle)": ("01_sichelzellanaemie/2_varianten.fasta", "HbS"),
    "β-Thalassämie Codon 41/42": ("01_sichelzellanaemie/2_varianten.fasta", "Codon41"),
    "HBB Schimpanse": ("02_haemoglobin_arten/haemoglobin_beta_7_arten.fasta", "Pan troglodytes"),
    "HBB Maus": ("02_haemoglobin_arten/haemoglobin_beta_7_arten.fasta", "Mus musculus"),
    "HBB Rind": ("02_haemoglobin_arten/haemoglobin_beta_7_arten.fasta", "Bos taurus"),
    "Pest-Plasmid pPCP1": ("05_pest_plasmid/pPCP1_plasmid.fasta", None),
}


def beispiel_fasta(name):
    """FASTA-Text eines eingebauten Beispiels (Namen siehe BEISPIEL_SEQUENZEN)."""
    datei, teil = BEISPIEL_SEQUENZEN[name]
    if datei == "kurz":
        return ">kurzes Testgen\n" + KURZES_TESTGEN
    return _als_fasta(*_erste_sequenz(os.path.join(BEISPIELE, datei), teil))


def bericht_html(pfad, einbetten=True):
    """Download-Knöpfe (Bericht als HTML, Ordner als ZIP) und der eingebettete Bericht."""
    ordner = os.path.dirname(pfad)
    zip_pfad = shutil.make_archive(ordner, "zip", ordner)
    teile = (_download_link(pfad, "Bericht (HTML)", "text/html") +
             _download_link(zip_pfad, "Alle Ergebnisse (ZIP)", "application/zip"))
    return teile + (_bericht_einbetten(pfad) if einbetten else "")


# ======================================================================
# Reiter 1: Analyse ganzer Projektordner
# ======================================================================

class _AnalyseReiter:
    def __init__(self, app):
        self.app = app
        breit = W.Layout(width="100%")
        self.quelle = W.ToggleButtons(
            options=[("Beispiele", "beispiele"), ("Eigene Daten (Colab)", "eigene"), ("Google Drive", "drive")],
            value="beispiele", style={"button_width": "170px"})
        self.drive_pfad = W.Text(value=DRIVE_STANDARD, description="Drive-Ordner:",
                                 layout=W.Layout(width="65%"), style={"description_width": "initial"})
        self.drive_knopf = W.Button(description="Drive verbinden", icon="link", layout=W.Layout(width="180px"))
        self.drive_box = W.HBox([self.drive_pfad, self.drive_knopf])

        self.projekte = W.SelectMultiple(options=[], rows=7, layout=breit)
        self.aktualisieren = W.Button(description="Liste aktualisieren", icon="refresh",
                                      layout=W.Layout(width="180px"))
        self.alle_knopf = W.Button(description="Alle auswählen", icon="check-square", layout=W.Layout(width="180px"))
        self.projekt_info = W.HTML()

        # Hochladen
        self.upload_name = W.Text(value="mein_projekt", description="Projektname:",
                                  style={"description_width": "initial"})
        self.uploader = W.FileUpload(accept=".fasta,.fa,.fna,.fas,.ffn,.txt,.json,.md", multiple=True,
                                     description="Dateien wählen")
        self.upload_knopf = W.Button(description="In Projekt speichern", icon="save", layout=W.Layout(width="200px"))
        self.upload_box = W.VBox([
            W.HTML("<b>Eigene Dateien hochladen</b> (FASTA; optional einstellungen.json, accessions.txt, "
                   "LIESMICH.md). Große Dateien besser über Google Drive."),
            W.HBox([self.upload_name, self.uploader, self.upload_knopf])])

        # Einstellungen
        stil = {"description_width": "160px"}
        mittel = W.Layout(width="460px")
        self.orf_min = W.Dropdown(description="ORF-Mindestlänge:", style=stil, layout=mittel, options=[
            ("wie im Projekt (meist 100 AS)", None), ("30 Aminosäuren", 30), ("60 Aminosäuren", 60),
            ("100 Aminosäuren", 100), ("150 Aminosäuren", 150), ("300 Aminosäuren", 300)])
        self.startcodons = W.Dropdown(description="Startcodons:", style=stil, layout=mittel,
                                      options=[("wie im Projekt", None)] + STARTCODON_OPTIONEN)
        self.min_ident = W.Dropdown(description="Varianten ab Identität:", style=stil, layout=mittel, options=[
            ("wie im Projekt (meist 80 %)", None), ("50 %", 50), ("70 %", 70), ("80 %", 80), ("90 %", 90)])
        self.referenz = W.Dropdown(description="Referenz:", style=stil, options=[("wie im Projekt", None)],
                                   layout=W.Layout(width="90%"), disabled=True)
        self.aussengruppe = W.SelectMultiple(description="Außengruppe:", style=stil, options=[], rows=4,
                                             layout=W.Layout(width="90%"), disabled=True)
        self.einzel_hinweis = W.HTML("<i>Referenz und Außengruppe kannst du wählen, wenn genau ein Projekt "
                                     "markiert ist.</i>")
        einstellungen = W.VBox([self.orf_min, self.startcodons, self.min_ident, self.einzel_hinweis,
                                self.referenz, self.aussengruppe,
                                W.HTML("<small>Außengruppe: mit Strg/Cmd mehrere wählen. Leer = wie im "
                                       "Projekt.</small>")])
        self.einstellungen = W.Accordion(children=[einstellungen])
        self.einstellungen.set_title(0, "⚙️ Einstellungen (optional)")
        self.einstellungen.selected_index = None

        self.start = W.Button(description="Analyse starten", button_style="primary", icon="play",
                              layout=W.Layout(width="200px", height="38px"))
        self.anzeigen = W.Checkbox(value=True, description="Bericht hier anzeigen", indent=False,
                                   layout=W.Layout(width="auto", margin="8px 0 0 20px"))
        self.protokoll = W.Output(layout=W.Layout(max_height="220px", overflow="auto",
                                                  border="1px solid #ddd", padding="4px"))
        self.ergebnis = W.Output()

        self.widget = W.VBox([
            W.HTML("<b>1. Woher kommen die Daten?</b>"), self.quelle, self.drive_box,
            W.HTML("<b>2. Welche Projekte?</b> (nichts markiert = alle; mehrere mit Strg/Cmd)"),
            self.projekte, W.HBox([self.alle_knopf, self.aktualisieren]), self.projekt_info,
            self.upload_box, self.einstellungen,
            W.HTML("<b>3. Los geht's</b>"), W.HBox([self.start, self.anzeigen]),
            self.protokoll, self.ergebnis])

        self.quelle.observe(lambda _: self.liste_aktualisieren(), names="value")
        self.projekte.observe(lambda _: self.auswahl_geaendert(), names="value")
        self.aktualisieren.on_click(lambda _: self.liste_aktualisieren())
        self.alle_knopf.on_click(lambda _: setattr(self.projekte, "value", tuple(self.projekte.options)))
        self.drive_knopf.on_click(self.drive_verbinden)
        self.upload_knopf.on_click(self.hochladen)
        self.start.on_click(self.starten)
        self.liste_aktualisieren()

    # --------------------------------------------------------------
    def basis_ordner(self):
        q = self.quelle.value
        if q == "beispiele":
            return BEISPIELE
        if q == "drive":
            return self.drive_pfad.value.rstrip("/")
        return EIGENE

    def ausgabe_ordner(self):
        q = self.quelle.value
        if q == "beispiele":
            return os.path.abspath("ergebnisse_beispiele")
        if q == "drive":
            return self.drive_pfad.value.rstrip("/") + "_ergebnisse"
        return os.path.abspath("ergebnisse")

    def liste_aktualisieren(self, markieren=None):
        q = self.quelle.value
        self.drive_box.layout.display = "flex" if q == "drive" else "none"
        self.upload_box.layout.display = "none" if q == "beispiele" else "flex"
        basis = self.basis_ordner()
        if q == "drive" and not os.path.isdir("/content/drive/MyDrive"):
            self.projekte.options = []
            self.projekt_info.value = _hinweis("Google Drive ist noch nicht verbunden. Klick auf "
                                               "<b>Drive verbinden</b>.", "warnung")
            return
        os.makedirs(basis, exist_ok=True)
        projekte = analyse.projekte_im_ordner(basis)
        if projekte == [os.path.normpath(basis)] and not os.listdir(basis):
            projekte = []
        self.pfade = {os.path.basename(p) if p != os.path.normpath(basis) else "(Dateien direkt im Ordner)": p
                      for p in projekte}
        self.projekte.options = list(self.pfade)
        if markieren in self.pfade:
            self.projekte.value = (markieren,)
        if not projekte:
            self.projekt_info.value = _hinweis(
                f"Der Ordner <code>{_e(basis)}</code> ist noch leer. Lade unten Dateien hoch "
                f"(oder hol sie im Reiter 🌐 NCBI).", "info")
        else:
            self.auswahl_geaendert()

    def auswahl_geaendert(self):
        gewaehlt = self.projekte.value
        einzeln = len(gewaehlt) == 1
        self.referenz.disabled = self.aussengruppe.disabled = not einzeln
        self.einzel_hinweis.layout.display = "none" if einzeln else "block"
        if not einzeln:
            self.referenz.options = [("wie im Projekt", None)]
            self.aussengruppe.options = []
            self.projekt_info.value = ""
            return
        pfad = self.pfade[gewaehlt[0]]
        try:
            namen = list(wz.lese_ordner(pfad)[0])
        except Exception:
            namen = []
        self.referenz.options = [("wie im Projekt (sonst erste Sequenz)", None)] + [(n, n) for n in namen]
        self.aussengruppe.options = namen
        info = f"<b>{_e(gewaehlt[0])}</b>: {len(namen)} Sequenz(en)"
        liesmich = os.path.join(pfad, "LIESMICH.md")
        if os.path.exists(liesmich):
            with open(liesmich, encoding="utf-8") as f:
                absaetze = [a.strip() for a in f.read().split("\n\n") if a.strip() and not a.startswith("#")]
            if absaetze:
                info += f"<br>{_e(absaetze[0][:500])}"
        if os.path.exists(os.path.join(pfad, "accessions.txt")) and not namen:
            info += "<br>Die Sequenzen werden beim Start automatisch von NCBI geladen."
        self.projekt_info.value = _kasten(info)

    def drive_verbinden(self, _):
        self.protokoll.clear_output()
        with self.protokoll:
            if not _in_colab():
                print("Google Drive gibt es nur in Google Colab.")
                return
            try:
                from google.colab import drive
                drive.mount("/content/drive")
            except Exception as fehler:
                print("Verbinden hat nicht geklappt:", fehler)
                print("Alternative: Ganz oben im Notebook DRIVE_VERBINDEN = True setzen und die Zelle neu starten.")
                return
        self.liste_aktualisieren()

    def hochladen(self, _):
        self.protokoll.clear_output()
        dateien = _hochgeladene_dateien(self.uploader)
        with self.protokoll:
            if not dateien:
                print("Erst auf 'Dateien wählen' klicken und Dateien aussuchen.")
                return
            name = _projektname(self.upload_name.value)
            ziel = os.path.join(self.basis_ordner(), name)
            os.makedirs(ziel, exist_ok=True)
            for dateiname, inhalt in dateien:
                with open(os.path.join(ziel, os.path.basename(dateiname)), "wb") as f:
                    f.write(inhalt)
                print(f"gespeichert: {name}/{os.path.basename(dateiname)} ({len(inhalt):,} Bytes)")
        self.liste_aktualisieren(markieren=name)

    def starten(self, _):
        self.protokoll.clear_output()
        self.ergebnis.clear_output()
        if not self.projekte.options:
            _zeige(self.ergebnis, _hinweis("Keine Projekte gefunden.", "warnung"))
            return
        gewaehlt = self.projekte.value or tuple(self.projekte.options)
        einstellungen = {"ORF_MIN_AMINOSAEUREN": self.orf_min.value, "STARTCODONS": self.startcodons.value,
                         "MIN_IDENTITAET_VARIANTEN": self.min_ident.value}
        if len(gewaehlt) == 1:
            einstellungen["REFERENZ"] = self.referenz.value
            einstellungen["AUSSENGRUPPE"] = list(self.aussengruppe.value) or None
        self.start.disabled = True
        self.start.description = "läuft …"
        try:
            with self.protokoll:
                ausgabe = self.ausgabe_ordner()
                if os.path.isdir(ausgabe):
                    shutil.rmtree(ausgabe)            # alte Ergebnisse dieses Ordners ersetzen
                ergebnisse, pfad = analyse.analysiere_projekte([self.pfade[g] for g in gewaehlt], ausgabe,
                                                               **einstellungen)
                zip_pfad = shutil.make_archive(ausgabe, "zip", ausgabe)
            links = (_download_link(pfad, "Bericht (HTML)", "text/html") +
                     _download_link(zip_pfad, "Alle Ergebnisse (ZIP)", "application/zip"))
            n_seq = sum(len(R.get("sequenzen", [])) for R in ergebnisse)
            text = _hinweis(f"✅ Fertig: {len(ergebnisse)} Projekt(e), {n_seq} Sequenzen. "
                            f"Ergebnisse liegen in <code>{_e(ausgabe)}</code>.", "ok") + links
            if self.anzeigen.value:
                text += _bericht_einbetten(pfad)
            _zeige(self.ergebnis, text)
        except Exception as fehler:
            _fehler_anzeigen(self.ergebnis, fehler)
        finally:
            self.start.disabled = False
            self.start.description = "Analyse starten"


# ======================================================================
# Reiter 2: eine Sequenz umwandeln und ORFs suchen
# ======================================================================

class _SequenzReiter:
    def __init__(self, app):
        self.app = app
        self.eingabe = W.Textarea(placeholder="DNA oder RNA hier einfügen (roh, FASTA oder mit Zahlen wie bei NCBI)",
                                  layout=W.Layout(width="100%", height="130px"))
        self.beispiel = W.Dropdown(layout=W.Layout(width="400px"), options=[
                                   ("Beispiel laden …", None), ("kurzes Testgen (51 Basen)", "kurz"),
                                            ("HBB-mRNA des Menschen (628 Basen)", "hbb"),
                                            ("Pest-Plasmid pPCP1 (9.609 Basen)", "ppcp1")])
        self.umwandeln = W.Button(description="Umwandeln & übersetzen", button_style="primary", icon="exchange",
                                  layout=W.Layout(width="230px"))
        self.orf_knopf = W.Button(description="ORFs suchen", button_style="info", icon="search",
                                  layout=W.Layout(width="160px"))
        self.orf_min = W.IntSlider(value=100, min=10, max=400, step=10, description="ORF ab (AS):",
                                   style={"description_width": "initial"})
        self.startcodons = W.Dropdown(options=STARTCODON_OPTIONEN, description="Startcodons:",
                                      style={"description_width": "initial"})
        self.ausgabe = W.Output()
        self.widget = W.VBox([
            W.HTML("Füge eine Sequenz ein oder lade ein Beispiel."), self.beispiel, self.eingabe,
            W.HBox([self.umwandeln, self.orf_knopf]), W.HBox([self.orf_min, self.startcodons]),
            self.ausgabe])
        self.beispiel.observe(self.beispiel_laden, names="value")
        self.umwandeln.on_click(self.zeige_umwandlung)
        self.orf_knopf.on_click(self.zeige_orfs)

    def beispiel_laden(self, aenderung):
        wahl = aenderung["new"]
        if wahl == "kurz":
            self.eingabe.value = ">kurzes Testgen\n" + KURZES_TESTGEN
            self.orf_min.value = 10
        elif wahl == "hbb":
            name, seq = _erste_sequenz(os.path.join(BEISPIELE, "01_sichelzellanaemie",
                                                    "1_HbA_normal_NM_000518.5.fasta"))
            self.eingabe.value = _als_fasta(name, seq)
            self.orf_min.value = 100
        elif wahl == "ppcp1":
            name, seq = _erste_sequenz(os.path.join(BEISPIELE, "05_pest_plasmid", "pPCP1_plasmid.fasta"))
            self.eingabe.value = _als_fasta(name, seq)
            self.orf_min.value = 60
            self.startcodons.value = STARTCODON_OPTIONEN[1][1]

    def zeige_umwandlung(self, _):
        _ausfuehren(self.ausgabe, lambda: sequenz_html(self.eingabe.value))

    def zeige_orfs(self, _):
        _ausfuehren(self.ausgabe, lambda: orfs_html(self.eingabe.value, self.orf_min.value, self.startcodons.value))


# ======================================================================
# Reiter 3: zwei Sequenzen vergleichen
# ======================================================================

class _VergleichReiter:
    def __init__(self, app):
        self.app = app
        feld = dict(layout=W.Layout(width="100%", height="110px"))
        self.s1 = W.Textarea(placeholder="Sequenz 1 = Referenz (z. B. das normale Gen)", **feld)
        self.s2 = W.Textarea(placeholder="Sequenz 2 = Variante (z. B. mit Mutation oder aus einer anderen Art)",
                             **feld)
        self.beispiel = W.Dropdown(layout=W.Layout(width="400px"), options=[
            ("Beispiel laden …", None), ("Sichelzelle: HbA gegen HbS", "hbs"),
            ("β-Thalassämie: HbA gegen Codon 41/42 (Deletion)", "thal"),
            ("Mensch gegen Schimpanse (HBB)", "pan"), ("Mensch gegen Maus (HBB)", "mus")])
        self.modus = W.ToggleButtons(options=[("Global (ganze Sequenzen)", "global"),
                                              ("Lokal (bester gemeinsamer Abschnitt)", "lokal")],
                                     style={"button_width": "260px"})
        self.gegenstrang = W.Checkbox(value=True, description="Gegenstrang automatisch prüfen", indent=False,
                                      layout=W.Layout(width="auto", margin="8px 0 0 20px"))
        self.knopf = W.Button(description="Vergleichen", button_style="primary", icon="balance-scale",
                              layout=W.Layout(width="180px"))
        self.ausgabe = W.Output()
        self.widget = W.VBox([self.beispiel, W.HTML("<b>Sequenz 1 (Referenz)</b>"), self.s1,
                              W.HTML("<b>Sequenz 2</b>"), self.s2, self.modus,
                              W.HBox([self.knopf, self.gegenstrang]), self.ausgabe])
        self.beispiel.observe(self.beispiel_laden, names="value")
        self.knopf.on_click(self.vergleichen)

    def beispiel_laden(self, aenderung):
        sichel = os.path.join(BEISPIELE, "01_sichelzellanaemie")
        arten = os.path.join(BEISPIELE, "02_haemoglobin_arten", "haemoglobin_beta_7_arten.fasta")
        hba = _erste_sequenz(os.path.join(sichel, "1_HbA_normal_NM_000518.5.fasta"))
        wahl = aenderung["new"]
        if wahl is None:
            return
        self.s1.value = _als_fasta(*hba)
        if wahl == "hbs":
            self.s2.value = _als_fasta(*_erste_sequenz(os.path.join(sichel, "2_varianten.fasta"), "HbS"))
        elif wahl == "thal":
            self.s2.value = _als_fasta(*_erste_sequenz(os.path.join(sichel, "2_varianten.fasta"), "Codon41"))
        elif wahl == "pan":
            self.s2.value = _als_fasta(*_erste_sequenz(arten, "Pan troglodytes"))
        elif wahl == "mus":
            self.s2.value = _als_fasta(*_erste_sequenz(arten, "Mus musculus"))
        self.modus.value = "global"

    def vergleichen(self, _):
        _ausfuehren(self.ausgabe, lambda: vergleich_html(self.s1.value, self.s2.value, self.modus.value,
                                                          self.gegenstrang.value))



def _varianten_html(referenz, erg):
    varianten = wz.normalisiere_indels(referenz, wz.finde_varianten(erg["a1"], erg["a2"]))
    if not varianten:
        return _hinweis("Keine Unterschiede im überlappenden Bereich.", "ok")
    cds = wz.finde_cds(referenz, min_as=10)
    if cds and (cds[1] - cds[0] + 1) < 0.5 * len(referenz) and len(referenz) > 300:
        cds = None                                  # längster ORF ist wohl kein echtes Gen
    if cds:
        wz.wirkung_der_varianten(referenz, varianten, cds)
        kopf = f"Gen in Sequenz 1: Position {cds[0]}–{cds[1]}. "
    else:
        kopf = "Kein Gen in Sequenz 1 erkannt, daher nur Positionen. "
    zeilen = []
    for v in varianten[:150]:
        if v["art"] == "SNV":
            aend = f"{v['ref']}→{v['alt']}"
        elif v["art"] == "Deletion":
            aend = f"−{v['ref'][:15]}{'…' if len(v['ref']) > 15 else ''} ({len(v['ref'])} bp)"
        else:
            aend = f"+{v['alt'][:15]}{'…' if len(v['alt']) > 15 else ''} ({len(v['alt'])} bp)"
        pos = str(v["start"]) if v["start"] == v["ende"] or v["art"] == "Insertion" else f"{v['start']}–{v['ende']}"
        code = lambda t: f"<code>{_e(t)}</code>" if t else ""
        zeilen.append([v["art"], pos, aend, code(v.get("c", "")), code(v.get("p", "")),
                       _e(v.get("wirkung", ""))])
    mehr = f"<br>… und {len(varianten) - 150} weitere." if len(varianten) > 150 else ""
    protein = ""
    if cds:
        p_ref = wz.uebersetze(referenz[cds[0] - 1:cds[1]])
        p_var = wz.uebersetze(_variante_rekonstruieren(erg, cds))
        protein = f"<br><b>Protein:</b> {_e(wz.vergleiche_proteine(p_ref, p_var))}"
        protein += _protein_vergleich(p_ref, p_var)
    zaehler = {}
    for v in varianten:
        k = v.get("wirkung", v["art"]).split(" (")[0]
        zaehler[k] = zaehler.get(k, 0) + 1
    zusammen = " · ".join(f"{k}: {n}" for k, n in zaehler.items())
    tabelle = _tabelle(["Art", "Position", "Änderung", "c. (Gen)", "p. (Protein)", "Wirkung"], zeilen) + mehr
    offen = " open" if len(varianten) <= 15 else ""
    return (f"<p style='margin:8px 0 2px'><b>{len(varianten)} Unterschied(e)</b> ({_e(zusammen)}). "
            f"{kopf}{protein}</p>"
            f"<details{offen}><summary style='cursor:pointer;color:#0072B2'>Tabelle aller Unterschiede"
            f"</summary>{tabelle}</details>")


def _variante_rekonstruieren(erg, cds):
    """Teil der Variante ab dem Startcodon der Referenz bis zu ihrem Ende (für die Translation)."""
    a1, a2 = erg["a1"], erg["a2"]
    rpos = 0
    for k, x in enumerate(a1):
        if x != "-":
            rpos += 1
            if rpos == cds[0]:
                return a2[k:].replace("-", "")
    return ""


# ======================================================================
# Reiter 4: Mutations-Labor
# ======================================================================

BEKANNTE_MUTATIONEN = [
    ("eigene Mutation", None),
    ("HbS (Sichelzellanämie): c.20A>T", ("snv", 20, "T", 1, "")),
    ("HbC: c.19G>A", ("snv", 19, "A", 1, "")),
    ("β-Thalassämie Codon 39: c.118C>T", ("snv", 118, "T", 1, "")),
    ("β-Thalassämie Codon 41/42: c.126_129del (4 bp)", ("del", 126, "", 4, "")),
    ("Startcodon zerstört: c.2T>C", ("snv", 2, "C", 1, "")),
    ("3 Basen löschen (In-frame): c.22_24del", ("del", 22, "", 3, "")),
    ("1 Base einfügen (Frameshift): nach c.30", ("ins", 30, "", 1, "A")),
]


class _MutationsReiter:
    def __init__(self, app):
        self.app = app
        stil = {"description_width": "120px"}
        self.referenz = W.Dropdown(description="Gen:", style=stil, layout=W.Layout(width="60%"), options=[
            ("HBB des Menschen (β-Globin)", "hbb"), ("Sequenz aus dem Reiter 🧬 Eine Sequenz", "feld")])
        self.bekannt = W.Dropdown(description="Vorlage:", style=stil, layout=W.Layout(width="60%"),
                                  options=BEKANNTE_MUTATIONEN)
        self.art = W.ToggleButtons(options=[("Base austauschen", "snv"), ("Basen löschen", "del"),
                                            ("Basen einfügen", "ins")], style={"button_width": "150px"})
        self.position = W.BoundedIntText(value=20, min=1, max=100000, description="Position c.:", style=stil)
        self.neue_base = W.Dropdown(options=list("ACGT"), value="T", description="neue Base:", style=stil)
        self.anzahl = W.BoundedIntText(value=1, min=1, max=300, description="Anzahl Basen:", style=stil)
        self.einfuegen = W.Text(value="A", description="einfügen:", style=stil)
        self.knopf = W.Button(description="Mutation einbauen", button_style="primary", icon="flask",
                              layout=W.Layout(width="200px"))
        self.ausgabe = W.Output()
        self.widget = W.VBox([
            W.HTML("Baue eine Mutation in ein Gen ein und sieh, was mit Codon und Protein passiert. "
                   "<b>c.1</b> ist das A des Startcodons ATG."),
            self.referenz, self.bekannt, self.art, self.position, self.neue_base, self.anzahl, self.einfuegen,
            self.knopf, self.ausgabe])
        self.bekannt.observe(self.vorlage, names="value")
        self.art.observe(lambda _: self.felder(), names="value")
        self.knopf.on_click(self.einbauen)
        self.felder()

    def felder(self):
        a = self.art.value
        self.neue_base.layout.display = "flex" if a == "snv" else "none"
        self.anzahl.layout.display = "flex" if a == "del" else "none"
        self.einfuegen.layout.display = "flex" if a == "ins" else "none"
        self.position.description = "nach Position c.:" if a == "ins" else "Position c.:"

    def vorlage(self, aenderung):
        if aenderung["new"] is None:
            return
        art, pos, base, anzahl, ins = aenderung["new"]
        self.referenz.value = "hbb"
        self.art.value = art
        self.position.value = pos
        if base:
            self.neue_base.value = base
        self.anzahl.value = anzahl
        if ins:
            self.einfuegen.value = ins

    def _gen(self):
        if self.referenz.value == "hbb":
            return _erste_sequenz(os.path.join(BEISPIELE, "01_sichelzellanaemie", "1_HbA_normal_NM_000518.5.fasta"))
        name, seq, _r, _h = _lese_eingabe(self.app.sequenz.eingabe.value)
        return name, seq

    def einbauen(self, _):
        text = None if self.referenz.value == "hbb" else self.app.sequenz.eingabe.value
        _ausfuehren(self.ausgabe, lambda: mutation_html(self.art.value, self.position.value, self.neue_base.value,
                                                         self.anzahl.value, self.einfuegen.value, text))


def _protein_vergleich(p_ref, p_mut, breite=60):
    """Beide Proteine untereinander. Gleicher Anfang und gleiches Ende werden erkannt, nur der
    veränderte Mittelteil wird markiert (so verschiebt eine In-frame-Deletion nicht alles)."""
    p = 0
    while p < min(len(p_ref), len(p_mut)) and p_ref[p] == p_mut[p]:
        p += 1
    s = 0
    while s < min(len(p_ref), len(p_mut)) - p and p_ref[-1 - s] == p_mut[-1 - s]:
        s += 1
    mitte_r, mitte_m = p_ref[p:len(p_ref) - s], p_mut[p:len(p_mut) - s]
    L = max(len(mitte_r), len(mitte_m))
    oben = p_ref[:p] + mitte_r.ljust(L, "-") + p_ref[len(p_ref) - s:]
    unten = p_mut[:p] + mitte_m.ljust(L, "-") + p_mut[len(p_mut) - s:]
    zeilen = []
    for start in range(0, len(oben), breite):
        a, b = oben[start:start + breite], unten[start:start + breite]
        o = u = ""
        for x, y in zip(a, b):
            if x == y:
                o += _e(x); u += _e(y)
            elif y == "-":
                o += _e(x); u += '<span style="color:#D55E00;font-weight:bold;">-</span>'
            elif x == "-":
                o += _e(x); u += f'<span style="background:#0072B2;color:#fff;">{_e(y)}</span>'
            else:
                o += _e(x); u += f'<span style="background:#D55E00;color:#fff;">{_e(y)}</span>'
        zeilen.append(f"{start + 1:>4} vorher   {o}\n{start + 1:>4} nachher  {u}\n")
    return ('<pre style="background:#f6f7f9;color:#1d1d1f;padding:10px;border-radius:6px;overflow-x:auto;'
            'font-size:12.5px;">' + "\n".join(zeilen) +
            '</pre><small>rot = andere Aminosäure, <b>-</b> = fehlt, blau = zusätzlich. Die Nummern zählen '
            'die Spalten dieser Darstellung.</small>')


# ======================================================================
# Reiter 5: NCBI
# ======================================================================

class _NcbiReiter:
    def __init__(self, app):
        self.app = app
        self.nummern = W.Textarea(placeholder="Eine Accession-Nummer pro Zeile, z. B.\nNM_000518.5\nNC_045512.2",
                                  layout=W.Layout(width="50%", height="120px"))
        self.projekt = W.Text(value="ncbi_download", description="Projektname:",
                              style={"description_width": "initial"})
        self.knopf = W.Button(description="Herunterladen", button_style="primary", icon="download",
                              layout=W.Layout(width="180px"))
        self.ziel = W.HTML()
        self.ausgabe = W.Output()
        self.widget = W.VBox([
            W.HTML("Sequenzen direkt von NCBI laden. Accession-Nummern findest du auf "
                   "<a href='https://www.ncbi.nlm.nih.gov/nuccore' target='_blank'>ncbi.nlm.nih.gov/nuccore</a>."),
            self.nummern, self.projekt, self.ziel, self.knopf, self.ausgabe])
        self.knopf.on_click(self.laden)

    def ziel_ordner(self):
        analyse_reiter = self.app.analyse
        if analyse_reiter.quelle.value == "drive" and os.path.isdir("/content/drive/MyDrive"):
            return analyse_reiter.basis_ordner(), "drive"
        return EIGENE, "eigene"

    def aktualisieren(self):
        basis, _q = self.ziel_ordner()
        self.ziel.value = f"Ziel: <code>{_e(os.path.join(basis, _projektname(self.projekt.value)))}</code>"

    def laden(self, _):
        self.ausgabe.clear_output()
        nummern = [n.strip() for n in re.split(r"[\s,;]+", self.nummern.value) if n.strip()]
        if not nummern:
            _zeige(self.ausgabe, _hinweis("Bitte mindestens eine Accession-Nummer eintragen.", "warnung"))
            return
        basis, quelle = self.ziel_ordner()
        name = _projektname(self.projekt.value)
        ziel = os.path.join(basis, name)
        self.knopf.disabled = True
        try:
            with self.ausgabe:
                for acc, fehler in wz.lade_von_ncbi(nummern, ziel):
                    print(f"{acc}: {'✅ geladen' if fehler is None else '⚠️ ' + fehler}")
            _zeige(self.ausgabe, _hinweis(f"Gespeichert in <code>{_e(ziel)}</code>. Im Reiter 📊 Analyse ist "
                                          f"das Projekt <b>{_e(name)}</b> schon ausgewählt.", "ok"))
            self.app.analyse.quelle.value = quelle
            self.app.analyse.liste_aktualisieren(markieren=name)
        finally:
            self.knopf.disabled = False


# ======================================================================
# Zusammenbauen
# ======================================================================

class Oberflaeche:
    def __init__(self):
        self.analyse = _AnalyseReiter(self)
        self.sequenz = _SequenzReiter(self)
        self.vergleich = _VergleichReiter(self)
        self.mutation = _MutationsReiter(self)
        self.ncbi = _NcbiReiter(self)
        self.reiter = W.Tab(children=[self.analyse.widget, self.sequenz.widget, self.vergleich.widget,
                                      self.mutation.widget, self.ncbi.widget])
        for i, titel in enumerate(["📊 Analyse", "🧬 Eine Sequenz", "⚖️ Vergleichen", "🩸 Mutations-Labor",
                                   "🌐 NCBI"]):
            self.reiter.set_title(i, titel)
        self.reiter.observe(lambda _: self.ncbi.aktualisieren(), names="selected_index")
        self.ncbi.projekt.observe(lambda _: self.ncbi.aktualisieren(), names="value")
        kopf = W.HTML("<h2 style='margin:4px 0'>🧰 Bioinformatik-Werkzeuge</h2>"
                      "<div style='margin-bottom:6px'>Reiter oben wählen, Felder ausfüllen, Knopf drücken. "
                      "Ergebnisse erscheinen jeweils darunter.</div>")
        self.widget = W.VBox([kopf, self.reiter])


def starte():
    """Zeigt die Oberfläche in der aktuellen Notebook-Zelle an."""
    print("Falls unter dieser Zeile keine Reiter erscheinen, nimm das Notebook App.ipynb (Formulare, "
          "funktioniert ohne Zusatz-Technik).")
    app = Oberflaeche()
    display(app.widget)
    return app
