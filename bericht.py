"""
📊 bericht.py – macht aus den Ergebnissen von analyse.py einen HTML-Bericht mit Tabellen und Diagrammen.

Der Bericht ist EINE Datei (bericht.html), die Bilder sind eingebettet. Man kann sie im Browser öffnen,
herunterladen oder weitergeben. In Colab wird sie direkt unter der Zelle angezeigt.
"""

import base64
import html
import io
import os
from datetime import datetime

from matplotlib.figure import Figure

# Farbenblind-freundliche Palette (Okabe & Ito)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#8C6D1F", "#000000"]

WIRKUNG_FARBEN = {
    "stumm": "#009E73",
    "Missense": "#E69F00",
    "Nonsense": "#D55E00",
    "Frameshift": "#CC79A7",
    "In-frame": "#0072B2",
    "Stoppverlust": "#8C6D1F",
    "Startcodon": "#7A3E9D",
    "außerhalb": "#9A9A9A",
    "ohne Gen": "#C4C4C4",
}


def wirkung_kategorie(v):
    w = v.get("wirkung", "")
    if not w:
        return "ohne Gen"
    for k in WIRKUNG_FARBEN:
        if w.startswith(k):
            return k
    return "ohne Gen"


def _e(text):
    return html.escape(str(text))


def _bild(fig, alt):
    puffer = io.BytesIO()
    fig.savefig(puffer, format="png", dpi=130, bbox_inches="tight", facecolor="white")
    daten = base64.b64encode(puffer.getvalue()).decode()
    return f'<img class="abb" alt="{_e(alt)}" src="data:image/png;base64,{daten}">'


def _achsen_stil(ax):
    for seite in ("top", "right"):
        ax.spines[seite].set_visible(False)
    ax.tick_params(labelsize=9)


def _gattung(name):
    return name.split("_")[0].split(" ")[0]


def _gattungsfarben(namen):
    gattungen = []
    for n in namen:
        g = _gattung(n)
        if g not in gattungen:
            gattungen.append(g)
    if len(gattungen) < 2 or len(gattungen) == len(namen) or len(gattungen) > len(PALETTE):
        return {n: "#1a1a1a" for n in namen}
    return {n: PALETTE[gattungen.index(_gattung(n))] for n in namen}


# ======================================================================
# Diagramme
# ======================================================================

def laenge_gc_abbildung(R):
    seqs = R["sequenzen"]
    namen = [s["name"] for s in seqs]
    hoehe = max(2.2, 0.32 * len(seqs) + 1)
    fig = Figure(figsize=(11, hoehe))
    ax1, ax2 = fig.subplots(1, 2, sharey=True)
    y = list(range(len(seqs)))[::-1]
    ax1.barh(y, [s["laenge"] for s in seqs], color=PALETTE[0], height=0.6)
    ax1.set_yticks(y)
    ax1.set_yticklabels([n[:45] for n in namen], fontsize=9)
    ax1.set_title("Länge (Basen)", fontsize=10)
    if max(s["laenge"] for s in seqs) > 50 * min(s["laenge"] for s in seqs):
        ax1.set_xscale("log")
        ax1.set_title("Länge (Basen, logarithmisch)", fontsize=10)
    ax2.barh(y, [s["gc"] for s in seqs], color=PALETTE[2], height=0.6)
    ax2.axvline(50, color="#888", lw=0.8, ls="--")
    ax2.set_xlim(0, 100)
    ax2.set_title("GC-Gehalt (%)", fontsize=10)
    for yy, s in zip(y, seqs):
        ax2.text(s["gc"] + 1, yy, f"{s['gc']:.1f}", va="center", fontsize=8)
    for ax in (ax1, ax2):
        _achsen_stil(ax)
    fig.tight_layout()
    return fig


def orf_karte_abbildung(R, max_sequenzen=8):
    seqs = [s for s in R["sequenzen"] if s.get("orfs")][:max_sequenzen]
    if not seqs:
        return None
    fig = Figure(figsize=(11, 1.25 * len(seqs) + 0.6))
    achsen = fig.subplots(len(seqs), 1, squeeze=False)[:, 0]
    raster = ["+1", "+2", "+3", "−1", "−2", "−3"]
    for ax, s in zip(achsen, seqs):
        for o in s["orfs"]:
            zeile = raster.index(("+" if o["strang"] == "+" else "−") + str(o["rahmen"]))
            links, rechts = min(o["start"], o["ende"]), max(o["start"], o["ende"])
            farbe = PALETTE[0] if o["strang"] == "+" else PALETTE[4]
            ax.barh(5 - zeile, rechts - links + 1, left=links, height=0.7, color=farbe)
        ax.set_yticks(range(6))
        ax.set_yticklabels(raster[::-1], fontsize=7)
        ax.set_xlim(0, s["laenge"])
        ax.set_title(f"{s['name']}  ({len(s['orfs'])} ORFs, {s['laenge']:,} bp)", fontsize=9, loc="left")
        ax.axhline(2.5, color="#bbb", lw=0.6)
        _achsen_stil(ax)
    achsen[-1].set_xlabel("Position (bp)  ·  blau = + Strang, rot = Gegenstrang", fontsize=9)
    fig.tight_layout()
    return fig


def identitaet_heatmap(R, titel="Paarweise Identität aus dem Alignment (%)", wert="identitaet"):
    namen = R["namen"]
    paare = R.get("paare") or []
    if len(namen) < 2 or not paare:
        return None
    index = {n: i for i, n in enumerate(namen)}
    import numpy as np
    M = np.full((len(namen), len(namen)), np.nan)
    for i in range(len(namen)):
        M[i, i] = 100
    for p in paare:
        M[index[p["a"]], index[p["b"]]] = M[index[p["b"]], index[p["a"]]] = p[wert]
    n = len(namen)
    groesse = min(9.0, 2.2 + 0.45 * n)
    fig = Figure(figsize=(groesse + 4.2, groesse))
    ax = fig.subplots()
    vmin = max(0, np.nanmin(M) - 2)
    bild = ax.imshow(M, cmap="viridis", vmin=vmin, vmax=100)
    # Spalten nur mit Nummer beschriften, Zeilen mit Nummer und Name: so bleibt alles lesbar
    ax.set_xticks(range(n))
    ax.set_xticklabels([str(i + 1) for i in range(n)], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{i + 1}  {nm[:45]}" for i, nm in enumerate(namen)], fontsize=8)
    if n <= 14:
        for i in range(n):
            for j in range(n):
                v = M[i, j]
                if np.isnan(v):
                    continue
                if v >= 99.995:
                    text = "100"
                elif n <= 9:
                    text = f"{v:.1f}"
                else:
                    text = f"{min(v, 99.4):.0f}"      # nie auf "100" runden, wenn es nicht 100 ist
                ax.text(j, i, text, ha="center", va="center", fontsize=7 if n <= 9 else 6.5,
                        color="white" if v < (vmin + 100) / 2 else "black")
    fig.colorbar(bild, ax=ax, shrink=0.7, label="Prozent")
    ax.set_title(titel, fontsize=10)
    fig.tight_layout()
    return fig


def protein_unterschiede_abbildung(P):
    """Wo im Protein liegen die Unterschiede zur Referenz? Ein Strich pro Unterschied."""
    vergleiche = [v for v in P.get("vergleiche", []) if v["unterschiede"]][:15]
    if not vergleiche:
        return None
    ref_laenge = next((len(e["protein"]) for e in P["proteine"] if e["name"] == P["referenz"]), None)
    fig = Figure(figsize=(11, 0.45 * len(vergleiche) + 1.4))
    ax = fig.subplots()
    arten = [("Austausch, ähnlich (BLOSUM62 > 0)", "#009E73"), ("Austausch, unähnlich (≤ 0)", "#D55E00"),
             ("Lücke (Insertion/Deletion)", "#7A7A7A")]
    gezeigt = set()
    for zeile, v in enumerate(vergleiche):
        y = len(vergleiche) - 1 - zeile
        for d in v["unterschiede"]:
            if d["art"] == "Austausch":
                label, farbe = arten[0] if d["punkte"] > 0 else arten[1]
            else:
                label, farbe = arten[2]
            ax.plot([d["position"], d["position"]], [y - 0.32, y + 0.32], color=farbe, lw=2.2,
                    label=None if label in gezeigt else label)
            gezeigt.add(label)
    ax.set_yticks(range(len(vergleiche)))
    ax.set_yticklabels([v["name"][:45] for v in vergleiche][::-1], fontsize=9)
    if ref_laenge:
        ax.set_xlim(0, ref_laenge + 1)
    ax.set_xlabel(f"Position im Referenzprotein ({P['referenz'][:50]})", fontsize=9)
    ax.set_ylim(-0.7, len(vergleiche) - 0.3)
    ax.legend(fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    _achsen_stil(ax)
    fig.tight_layout()
    return fig


def varianten_abbildung(V):
    erg = V.get("ergebnisse") or []
    if not erg:
        return None
    kategorien = [k for k in WIRKUNG_FARBEN if any(wirkung_kategorie(v) == k for e in erg for v in e["varianten"])]
    if not kategorien:
        return None
    fig = Figure(figsize=(11, max(2, 0.4 * len(erg) + 1.2)))
    ax = fig.subplots()
    y = list(range(len(erg)))[::-1]
    links = [0] * len(erg)
    for k in kategorien:
        werte = [sum(1 for v in e["varianten"] if wirkung_kategorie(v) == k) for e in erg]
        ax.barh(y, werte, left=links, color=WIRKUNG_FARBEN[k], label=k, height=0.6)
        links = [a + b for a, b in zip(links, werte)]
    ax.set_yticks(y)
    ax.set_yticklabels([e["name"][:45] for e in erg], fontsize=9)
    ax.set_xlabel("Anzahl Unterschiede zur Referenz", fontsize=9)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(fontsize=8, ncol=min(5, len(kategorien)), loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    _achsen_stil(ax)
    fig.tight_layout()
    return fig


def _blaetter(k):
    return [k["name"]] if not k["kinder"] else [b for kind, _ in k["kinder"] for b in _blaetter(kind)]


def _zeichne_baum(ax, knoten, titel, farben):
    y_pos = {n: i for i, n in enumerate(_blaetter(knoten))}
    x_max = [0.0]

    def zeichne(k, x):
        if not k["kinder"]:
            ax.text(x, y_pos[k["name"]], "  " + k["name"].replace("_", " ") if "(" not in k["name"] else "  " + k["name"],
                    va="center", fontsize=9, color=farben.get(k["name"], "black"))
            x_max[0] = max(x_max[0], x)
            return y_pos[k["name"]]
        ys = []
        for kind, laenge in k["kinder"]:
            l = max(laenge, 0.0)                      # negative NJ-Astlängen als 0 zeichnen
            y = zeichne(kind, x + l)
            ax.plot([x, x + l], [y, y], color="#333", lw=1.3)
            ys.append(y)
        ax.plot([x, x], [min(ys), max(ys)], color="#333", lw=1.3)
        return sum(ys) / len(ys)

    zeichne(knoten, 0.0)
    ax.set_title(titel, fontsize=11)
    ax.set_yticks([])
    ax.set_xlabel("Distanz (Substitutionen pro Position)", fontsize=9)
    for seite in ("left", "right", "top"):
        ax.spines[seite].set_visible(False)
    ax.set_xlim(0, (x_max[0] or 1) * 1.9)
    ax.set_ylim(-0.7, len(y_pos) - 0.3)
    ax.invert_yaxis()


def baum_abbildung(T):
    namen = _blaetter(T["upgma"])
    farben = _gattungsfarben(namen)
    fig = Figure(figsize=(13, max(3, 0.42 * len(namen) + 1)))
    ax1, ax2 = fig.subplots(1, 2)
    _zeichne_baum(ax1, T["upgma"], "UPGMA", farben)
    _zeichne_baum(ax2, T["nj"], "Neighbor-Joining", farben)
    fig.tight_layout()
    return fig


# ======================================================================
# HTML
# ======================================================================

CSS = """
:root { --text:#1d1d1f; --leise:#5f6368; --linie:#e3e3e8; --flaeche:#f6f7f9; --akzent:#0072B2; }
* { box-sizing: border-box; }
body { margin:0; background:#ffffff; color:var(--text);
       font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
.seite { max-width: 1120px; margin: 0 auto; padding: 28px 18px 60px; }
h1 { font-size: 28px; margin: 0 0 4px; }
h2 { font-size: 22px; margin: 44px 0 6px; padding-top: 14px; border-top: 3px solid var(--akzent); }
h3 { font-size: 16px; margin: 26px 0 8px; }
.untertitel { color: var(--leise); margin: 0 0 18px; }
.box { background: var(--flaeche); border: 1px solid var(--linie); border-radius: 10px; padding: 12px 16px; margin: 10px 0; }
.warnung { background: #fff6e5; border: 1px solid #f0c36d; border-radius: 10px; padding: 10px 14px; margin: 10px 0; }
.kacheln { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; margin: 12px 0; }
.kachel { background: var(--flaeche); border: 1px solid var(--linie); border-radius: 10px; padding: 10px 14px; }
.kachel .zahl { font-size: 22px; font-weight: 650; }
.kachel .text { color: var(--leise); font-size: 13px; }
.tabelle-huelle { overflow-x: auto; margin: 8px 0 14px; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--linie); vertical-align: top; }
th { background: var(--flaeche); font-weight: 600; white-space: nowrap; }
td.zahl, th.zahl { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
tr:hover td { background: #fafbfc; }
.abb { max-width: 100%; height: auto; display: block; margin: 8px 0 4px; }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 13px; }
.kette { overflow-x: auto; white-space: pre; background: var(--flaeche); border: 1px solid var(--linie);
         border-radius: 8px; padding: 10px 12px; }
.marke { display: inline-block; padding: 1px 8px; border-radius: 999px; color: #fff; font-size: 12px; white-space: nowrap; }
details { margin: 8px 0; } summary { cursor: pointer; color: var(--akzent); }
nav ol { margin: 6px 0 0 18px; padding: 0; } nav a { color: var(--akzent); text-decoration: none; }
.leise { color: var(--leise); font-size: 13px; }
"""


def _tabelle(kopf, zeilen, zahlspalten=()):
    t = ['<div class="tabelle-huelle"><table><thead><tr>']
    for i, k in enumerate(kopf):
        t.append(f'<th class="{"zahl" if i in zahlspalten else ""}">{_e(k)}</th>')
    t.append("</tr></thead><tbody>")
    for z in zeilen:
        t.append("<tr>" + "".join(
            f'<td class="{"zahl" if i in zahlspalten else ""}">{w}</td>' for i, w in enumerate(z)) + "</tr>")
    t.append("</tbody></table></div>")
    return "".join(t)


def _marke(kategorie, text=None):
    return f'<span class="marke" style="background:{WIRKUNG_FARBEN[kategorie]}">{_e(text or kategorie)}</span>'


def _projekt_html(R, nr):
    h = [f'<h2 id="p{nr}">{nr}. {_e(R["titel"])}</h2>']
    if not R.get("sequenzen"):
        if R.get("beschreibung"):
            h.append(f"<p>{_e(R['beschreibung'].split(chr(10) + chr(10))[0])}</p>")
        for w in R.get("warnungen", []):
            h.append(f'<div class="warnung">⚠️ {_e(w)}</div>')
        h.append('<div class="box">Keine Sequenzen vorhanden, daher keine Analyse.</div>')
        return "\n".join(h)
    if R.get("beschreibung"):
        for absatz in R["beschreibung"].split("\n\n"):
            if absatz.startswith("#"):
                continue
            h.append(f"<p>{_e(absatz)}</p>")
    for w in R.get("warnungen", []):
        h.append(f'<div class="warnung">⚠️ {_e(w)}</div>')
    seqs = R["sequenzen"]
    S = set(R.get("schritte") or ("orfs", "alignment", "varianten", "protein", "stammbaum"))
    if len(S) < 5:
        from analyse import SCHRITT_NAMEN
        h.append(f'<p class="leise">Gewählte Analyse: Übersicht + '
                 f'{", ".join(SCHRITT_NAMEN[s] for s in ("orfs", "alignment", "varianten", "protein", "stammbaum") if s in S)}</p>')
    anzahl_orfs = sum(len(s.get("orfs", [])) for s in seqs)
    anzahl_var = sum(len(e["varianten"]) for e in (R.get("varianten") or {}).get("ergebnisse", []))
    h.append('<div class="kacheln">')
    P = R.get("protein") or {}
    kacheln = [(len(seqs), "Sequenzen"), (f"{sum(s['laenge'] for s in seqs):,}", "Basen insgesamt")]
    if "orfs" in S:
        kacheln.append((anzahl_orfs, "ORFs gefunden"))
    if "alignment" in S:
        kacheln.append((len(R.get("paare") or []), "Alignments"))
    if "varianten" in S:
        kacheln.append((anzahl_var, "Varianten"))
    if "protein" in S:
        kacheln.append((sum(1 for e in P.get("proteine", []) if e["protein"]), "Proteine"))
    if "stammbaum" in S:
        kacheln.append(("ja" if R.get("baum") else "–", "Stammbaum"))
    for zahl, text in kacheln:
        h.append(f'<div class="kachel"><div class="zahl">{_e(zahl)}</div><div class="text">{text}</div></div>')
    h.append("</div>")

    # Sequenzen
    h.append("<h3>Sequenzen, GC-Gehalt und ORFs</h3>")
    h.append(_tabelle(["Name", "Datei", "Länge (bp)", "GC %", "N", "ORFs", "Längster ORF (AS)", "Anteil in ORFs"],
                      [[_e(s["name"]), _e(s["datei"]), f"{s['laenge']:,}", f"{s['gc']:.1f}", f"{s['N']:,}",
                        len(s.get("orfs", [])), s.get("laengster_orf", 0), f"{s.get('orf_anteil', 0):.0f} %"]
                       for s in seqs], zahlspalten=(2, 3, 4, 5, 6, 7)))
    h.append(_bild(laenge_gc_abbildung(R), "Länge und GC-Gehalt"))
    E = R["einstellungen"]
    karte = orf_karte_abbildung(R) if "orfs" in S else None
    if karte:
        h.append(f'<p class="leise">ORF-Karte: jedes Kästchen ist ein offenes Leseraster ab '
                 f'{E["ORF_MIN_AMINOSAEUREN"]} Aminosäuren (Startcodon {"/".join(E["STARTCODONS"])}), '
                 f'aufgeteilt auf die 6 Leseraster. Nicht jeder ORF ist ein echtes Gen.</p>')
        h.append(_bild(karte, "ORF-Karte"))

    # Alignments
    if "alignment" in S:
        h.append(_alignment_html(R))
    if "varianten" in S:
        h.append(_varianten_html(R, E))
    if "protein" in S:
        h.append(_protein_html(P))
    if "stammbaum" in S:
        h.append(_stammbaum_html(R))
    h.append(f'<p class="leise">Dateien zu diesem Projekt: {_e(R["ausgabe_ordner"])}/ · '
             f'Rechenzeit {R.get("dauer", 0):.1f} s</p>')
    return "\n".join(h)


def _alignment_html(R):
    h = ["<h3>Paarweiser Vergleich (Alignment)</h3>"]
    if R.get("alignment_hinweis"):
        h.append(f'<div class="box">{_e(R["alignment_hinweis"])}</div>')
    hm = identitaet_heatmap(R)
    if hm:
        h.append(_bild(hm, "Identitäts-Heatmap"))
        zeilen = [[_e(p["a"]), _e(p["b"]), "Gegenstrang" if p["strang"] == "−" else "gleich",
                   f"{p['identitaet']:.1f}", p["mismatches"], p["luecken"]] for p in R["paare"]]
        h.append(f"<details><summary>Tabelle aller {len(zeilen)} Paare</summary>"
                 + _tabelle(["Sequenz 1", "Sequenz 2", "Strang", "Identität %", "Mismatches", "Lücken"], zeilen,
                            zahlspalten=(3, 4, 5)) + "</details>")
        h.append('<p class="leise">Unter etwa 50 % Identität mit vielen Lücken sind sich DNA-Sequenzen nicht ähnlicher '
                 'als Zufallssequenzen.</p>')

    return "\n".join(h)


def _varianten_html(R, E):
    h = []
    V = R.get("varianten") or {}
    h.append("<h3>Variantenanalyse: Mutationen und ihre Wirkung</h3>")
    if V.get("hinweis"):
        h.append(f'<div class="box">{_e(V["hinweis"])}</div>')
    if V.get("referenz") and V.get("ergebnisse"):
        h.append(f'<p>Referenz: <b>{_e(V["referenz"])}</b>')
        if V.get("cds"):
            h.append(f' · Gen (CDS) an Position {V["cds"][0]}–{V["cds"][1]} ({_e(V["cds_quelle"])}), '
                     f'{len(V["protein"])} Aminosäuren</p>')
            mrna, protein = V["mrna"], V["protein"]
            n_codons = min(15, len(protein))
            codons = " ".join(mrna[i * 3:i * 3 + 3] for i in range(n_codons))
            aminos = " ".join(f"{a:^3}" for a in protein[:n_codons])
            h.append('<div class="kette mono">'
                     f'DNA (codierend)   5\'-{_e(V["dna"][:n_codons * 3])}…\n'
                     f'mRNA (T → U)      5\'-{_e(mrna[:n_codons * 3])}…\n'
                     f'Codons               {_e(codons)} …\n'
                     f'Protein              {_e(aminos)} …</div>')
        else:
            h.append("</p>")
        abb = varianten_abbildung(V)
        if abb:
            h.append(_bild(abb, "Varianten nach Wirkung"))
        uebersicht = []
        for e in V["ergebnisse"]:
            kat = {}
            for v in e["varianten"]:
                kat[wirkung_kategorie(v)] = kat.get(wirkung_kategorie(v), 0) + 1
            marken = " ".join(_marke(k, f"{k}: {n}") for k, n in kat.items()) or "keine Unterschiede"
            uebersicht.append([_e(e["name"]), f"{e['identitaet']:.1f}", len(e["varianten"]), marken,
                               _e(e.get("protein_vergleich", "–"))])
        h.append(_tabelle(["Variante", "Identität %", "Unterschiede", "Wirkung", "Protein im Vergleich zur Referenz"],
                          uebersicht, zahlspalten=(1, 2)))
        for e in V["ergebnisse"]:
            if not e["varianten"]:
                continue
            zeilen = [[v["art"], v["position"], _e(v["aenderung"]), f'<span class="mono">{_e(v.get("c", ""))}</span>',
                       f'<span class="mono">{_e(v.get("p", ""))}</span>',
                       _marke(wirkung_kategorie(v), v.get("wirkung") or "ohne Gen")] for v in e["varianten"][:60]]
            rest = (f'<p class="leise">… und {len(e["varianten"]) - 60} weitere in varianten.csv</p>'
                    if len(e["varianten"]) > 60 else "")
            offen = " open" if len(e["varianten"]) <= 6 else ""
            h.append(f"<details{offen}><summary>Alle Unterschiede: {_e(e['name'])} ({len(e['varianten'])})</summary>"
                     + _tabelle(["Art", "Position", "Änderung", "c. (Gen)", "p. (Protein)", "Wirkung"], zeilen,
                                zahlspalten=(1,)) + rest + "</details>")
        if V.get("uebersprungen"):
            h.append(f'<p class="leise">Nicht in der Variantenanalyse (unter {E["MIN_IDENTITAET_VARIANTEN"]} % '
                     f'Identität zur Referenz oder zu lang): {_e(", ".join(V["uebersprungen"]))}</p>')

    return "\n".join(h)


def _stammbaum_html(R):
    h = []
    T = R.get("baum")
    h.append("<h3>Stammbaum</h3>")
    if not T:
        h.append('<div class="box">Für einen Stammbaum braucht man mindestens 3 Sequenzen.</div>')
    else:
        h.append(f"<p>Methode: {_e(T['methode'])}"
                 + (f" · Neighbor-Joining bewurzelt mit Außengruppe: {_e(', '.join(T['aussengruppe']))}"
                    if T.get("aussengruppe") else "") + "</p>")
        for w in T["warnungen"]:
            h.append(f'<div class="warnung">⚠️ {_e(w)}</div>')
        h.append(_bild(baum_abbildung(T), "Stammbaum"))
        from werkzeuge import newick
        h.append("<details><summary>Newick-Format (z. B. für iTOL)</summary>"
                 f'<div class="kette mono">{_e(newick(T["upgma"]))};\n\n{_e(newick(T["nj"]))};</div></details>')
    return "\n".join(h)


def _protein_html(P):
    if not P:
        return ""
    h = ["<h3>Protein-Vergleich</h3>"]
    h.append('<p class="leise">Aus jeder Sequenz wird das Gen (längster ORF ab ATG) in ein Protein übersetzt. '
             'Die Proteine werden mit der BLOSUM62-Matrix aligniert, wie bei BLASTP. <b>Identität</b> = gleiche '
             'Aminosäure, <b>Ähnlichkeit</b> = gleiche oder chemisch ähnliche Aminosäure (BLOSUM62-Wert über 0).</p>')
    vergleich = {v["name"]: v for v in P.get("vergleiche", [])}
    zeilen = []
    for e in P.get("proteine", []):
        v = vergleich.get(e["name"])
        if e["name"] == P.get("referenz"):
            werte = ["Referenz", "", ""]
        elif v:
            werte = [f"{v['identitaet']:.1f}", f"{v['aehnlichkeit']:.1f}", len(v["unterschiede"])]
        else:
            werte = ["–", "–", "–"]
        zeilen.append([_e(e["name"]), _e(e["quelle"]), len(e["protein"]) if e["protein"] else "–"] + werte)
    h.append(_tabelle(["Name", "Protein aus", "Länge (AS)", "Identität %", "Ähnlichkeit %", "Unterschiede"],
                      zeilen, zahlspalten=(2, 3, 4, 5)))
    if P.get("hinweis"):
        h.append(f'<div class="box">{_e(P["hinweis"])}</div>')
    abb = protein_unterschiede_abbildung(P)
    if abb:
        h.append(_bild(abb, "Unterschiede entlang des Proteins"))
    namen = [e["name"] for e in P.get("proteine", []) if e["protein"]]
    if len(namen) >= 2 and P.get("paare"):
        hm = identitaet_heatmap({"namen": namen, "paare": P["paare"]}, "Protein: paarweise Identität (%)")
        h.append(_bild(hm, "Protein-Identität"))
    for v in P.get("vergleiche", []):
        if not v["unterschiede"]:
            continue
        zeilen = []
        for d in v["unterschiede"][:150]:
            if d["art"] == "Austausch":
                text = f"{DREI.get(d['ref'], d['ref'])}{d['position']}{DREI.get(d['alt'], d['alt'])}"
                punkte = d["punkte"]
                farbe = "#009E73" if punkte > 0 else "#D55E00"
                punkte = f'<span style="color:{farbe};font-weight:600">{punkte:+d}</span>'
            else:
                if d["art"] == "Insertion":
                    text = f"+{d['alt']} (nach Pos. {d['position']})"
                else:
                    pos = d["position"] if d["position"] == d["ende"] else f"{d['position']}–{d['ende']}"
                    text = f"−{d['ref']} (Pos. {pos})"
                punkte = ""
            zeilen.append([d["art"], f'<span class="mono">{_e(text)}</span>', punkte, _e(d["bewertung"])])
        offen = " open" if len(v["unterschiede"]) <= 10 else ""
        h.append(f"<details{offen}><summary>Unterschiede: {_e(v['name'])} ({len(v['unterschiede'])})</summary>"
                 + _tabelle(["Art", "Änderung", "BLOSUM62", "Bewertung"], zeilen) + "</details>")
    return "\n".join(h)


DREI = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
        "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
        "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}


def erstelle_html_bericht(ergebnisse, pfad):
    ergebnisse = [R for R in ergebnisse if R and (R.get("sequenzen") or R.get("warnungen"))]
    teile = ['<!doctype html><html lang="de"><head><meta charset="utf-8">',
             '<meta name="viewport" content="width=device-width, initial-scale=1">',
             f"<title>Bioinformatik-Bericht</title><style>{CSS}</style></head><body><div class='seite'>",
             "<h1>🧬 Bioinformatik-Bericht</h1>",
             f'<p class="untertitel">Erstellt am {datetime.now().strftime("%d.%m.%Y um %H:%M")} · '
             f'{len(ergebnisse)} Projekt(e)</p>']
    if len(ergebnisse) > 1:
        zeilen = []
        for i, R in enumerate(ergebnisse, 1):
            V = R.get("varianten")
            P = R.get("protein")
            if not R.get("sequenzen"):
                zeilen.append([f'<a href="#p{i}">{i}. {_e(R["titel"])}</a>', 0, "–", "–", "–", "–", "–"])
                continue
            zeilen.append([f'<a href="#p{i}">{i}. {_e(R["titel"])}</a>', len(R["sequenzen"]),
                           f"{sum(s['laenge'] for s in R['sequenzen']):,}",
                           sum(len(s.get("orfs", [])) for s in R["sequenzen"]),
                           sum(len(e["varianten"]) for e in V.get("ergebnisse", [])) if V else "–",
                           sum(1 for e in P.get("proteine", []) if e["protein"]) if P else "–",
                           "ja" if R.get("baum") else "–"])
        teile.append("<h3>Übersicht aller Projekte</h3>")
        teile.append(_tabelle(["Projekt", "Sequenzen", "Basen", "ORFs", "Varianten", "Proteine", "Stammbaum"],
                              zeilen, zahlspalten=(1, 2, 3, 4, 5)))
    for i, R in enumerate(ergebnisse, 1):
        teile.append(_projekt_html(R, i))
    teile.append("""
<h2>Kurz erklärt</h2>
<div class="box">
<p><b>ORF</b> (offenes Leseraster): Abschnitt von einem Startcodon bis zum nächsten Stoppcodon im selben Raster. Mögliche Gene, aber nicht jeder ORF ist ein Gen.</p>
<p><b>Identität</b>: Anteil gleicher Basen im Alignment. <b>Lücken</b> stehen für Insertionen oder Deletionen.</p>
<p><b>c.</b>-Position: Position im Gen, c.1 ist das A des Startcodons ATG. <b>p.</b>: Änderung im Protein, z. B. p.Glu7Val = Glutamat an Position 7 wird zu Valin (HGVS-Schreibweise).</p>
<p><b>stumm</b>: gleiche Aminosäure · <b>Missense</b>: andere Aminosäure · <b>Nonsense</b>: vorzeitiges Stoppcodon · <b>Frameshift</b>: Insertion/Deletion, die das Leseraster verschiebt · <b>In-frame</b>: Insertion/Deletion von ganzen Codons.</p>
<p><b>Protein-Identität / Ähnlichkeit</b>: Identität zählt nur gleiche Aminosäuren. Ähnlichkeit zählt auch Austausche zwischen chemisch ähnlichen Aminosäuren (BLOSUM62 &gt; 0, z. B. Ile↔Val). <b>BLOSUM62</b>: Punktetabelle aus echten Proteinfamilien; positiv = häufiger Austausch, negativ = seltener.</p>
<p><b>UPGMA</b> nimmt eine gleichmäßige Mutationsrate an (molekulare Uhr), <b>Neighbor-Joining</b> nicht. Ein Stammbaum ist eine Schätzung.</p>
</div>
</div></body></html>""")
    os.makedirs(os.path.dirname(os.path.abspath(pfad)), exist_ok=True)
    with open(pfad, "w", encoding="utf-8") as f:
        f.write("\n".join(teile))
    return pfad
