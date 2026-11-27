"""3D-Strukturen von Proteinen: laden, lesen und auswerten.

- lade_struktur("4HHB")            Struktur aus der Proteindatenbank PDB (rcsb.org)
- lade_struktur("P68871")          AlphaFold-Modell zu einer UniProt-Nummer
- lade_struktur("meine.pdb")       eigene Datei
- lese_pdb(text)                   PDB-Text in Atome, Reste, Ketten und Liganden zerlegen
- zugaenglichkeit(S)               Wie viel Oberfläche jeder Aminosäure frei liegt (Shrake-Rupley)
- kontakte(S, rest, abstand)       Welche Aminosäuren einen Liganden (z. B. Häm) berühren
- ordne_protein_zu(S, protein)     Positionen eines Proteins den Resten der Struktur zuordnen
- uniprot_merkmale("P68871")       bekannte aktive Stellen und Bindestellen aus UniProt

Gesichert: Abstände und Kontakte sind Geometrie aus der gemessenen (oder vorhergesagten) Struktur.
Faustregeln: "vergraben" (unter 20 % frei) und Kontakt bis 4 Å sind übliche Grenzwerte, keine Naturgesetze.
"""

import gzip
import json
import math
import os
import re
import urllib.error
import urllib.request

import numpy as np

from werkzeuge import ausrichten_protein, BLOSUM62

DREI_ZU_EINS = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O"}
# häufige chemisch veränderte Aminosäuren (stehen in PDB-Dateien als HETATM, gehören aber zur Kette)
VERAENDERT = {"MSE": "M", "CSO": "C", "CME": "C", "OCS": "C", "CSD": "C", "SEP": "S", "TPO": "T", "PTR": "Y",
              "HYP": "P", "MLY": "K", "KCX": "K", "LLP": "K", "M3L": "K", "ALY": "K", "PCA": "E", "CGU": "E",
              "MLE": "L", "NLE": "L", "SAC": "S", "CSS": "C", "SCH": "C", "TYS": "Y", "HIC": "H", "NEP": "H"}
WASSER = {"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"}
# kleine Moleküle, die meist aus der Kristallisation stammen (Puffer, Frostschutz, Salze)
ZUSAETZE = {"SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "ACT", "FMT", "DMS", "MPD", "BME",
            "TRS", "EPE", "MES", "CIT", "IMD", "CL", "NA", "K", "BR", "IOD", "NO3", "SCN", "ACY", "IPA", "EOH",
            "MOH", "BU3", "PE4", "2PE", "12P", "15P", "MRD", "NH4", "AZI", "UNX"}
LIGAND_NAMEN = {"HEM": "Häm", "HEC": "Häm c", "HEA": "Häm a", "FAD": "FAD", "FMN": "FMN", "NAD": "NAD⁺",
                "NAP": "NADP⁺", "NDP": "NADPH", "ATP": "ATP", "ADP": "ADP", "ANP": "ATP-Analogon",
                "GTP": "GTP", "GDP": "GDP", "SAM": "S-Adenosylmethionin", "SAH": "S-Adenosylhomocystein",
                "COA": "Coenzym A", "PLP": "Pyridoxalphosphat", "TPP": "Thiaminpyrophosphat", "BTN": "Biotin",
                "CLA": "Chlorophyll a", "CHL": "Chlorophyll b", "RET": "Retinal", "OXY": "Sauerstoff (O₂)",
                "CMO": "Kohlenmonoxid", "ZN": "Zink-Ion", "MG": "Magnesium-Ion", "CA": "Calcium-Ion",
                "FE": "Eisen-Ion", "FE2": "Eisen(II)-Ion", "CU": "Kupfer-Ion", "MN": "Mangan-Ion",
                "SF4": "Eisen-Schwefel-Cluster", "FES": "Eisen-Schwefel-Cluster", "STR": "Progesteron",
                "XK2": "Hemmstoff XK263", "GLC": "Glucose", "BGC": "Glucose", "NAG": "N-Acetylglucosamin",
                "MAN": "Mannose"}

# Atomradien in Å (wie Biopython Bio.PDB.SASA), unbekannte Elemente 2.0
RADIEN = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "F": 1.47, "NA": 2.27, "MG": 1.73, "P": 1.8, "S": 1.8,
          "CL": 1.75, "K": 2.75, "CA": 2.31, "NI": 1.63, "CU": 1.4, "ZN": 1.39, "SE": 1.9, "BR": 1.85,
          "CD": 1.58, "I": 1.98, "HG": 1.55, "FE": 2.0, "MN": 2.0, "CO": 2.0}
# größtmögliche Oberfläche je Aminosäure in Å² (Tien et al. 2013, theoretisch), für die relative Zugänglichkeit
MAX_FLAECHE = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0, "E": 223.0, "G": 104.0,
               "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0,
               "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}
VERGRABEN_BIS = 20.0          # % relative Zugänglichkeit
OBERFLAECHE_AB = 50.0

_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_UNIPROT = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


# ----------------------------------------------------------------------
# Laden
# ----------------------------------------------------------------------

def _hole(url, zeitlimit=60):
    anfrage = urllib.request.Request(url, headers={"User-Agent": "bioinfo-lernprojekt"})
    try:
        with urllib.request.urlopen(anfrage, timeout=zeitlimit) as antwort:
            daten = antwort.read()
    except urllib.error.HTTPError:
        raise
    except (urllib.error.URLError, TimeoutError, OSError) as fehler:
        raise ConnectionError(f"Keine Verbindung zu {url.split('/')[2]} ({fehler}). "
                              f"Ist Internet da? In Colab sollte es gehen.") from None
    if url.endswith(".gz"):
        daten = gzip.decompress(daten)
    return daten.decode("utf-8", errors="replace")


def kennung_aus_text(text):
    """Holt aus einem Auswahl-Text wie "4HHB – Hämoglobin" oder "AF-P68871 …" die Kennung."""
    erstes = text.strip().split()[0] if text.strip() else ""
    return erstes.rstrip(",;:–-")


def lade_struktur(eingabe, ordner="strukturen"):
    """Lädt eine Struktur. eingabe = PDB-ID (4HHB), UniProt-Nummer für AlphaFold (P68871 oder AF-P68871)
    oder Pfad zu einer .pdb/.pdb.gz-Datei. Gibt (pdb_text, quelle) zurück; Downloads werden im ordner gespeichert."""
    eingabe = eingabe.strip()
    if os.path.isfile(eingabe):
        oeffnen = gzip.open if eingabe.endswith(".gz") else open
        with oeffnen(eingabe, "rt", errors="replace") as f:
            return f.read(), f"Datei {os.path.basename(eingabe)}"
    kennung = kennung_aus_text(eingabe)
    os.makedirs(ordner, exist_ok=True)
    if kennung.upper().startswith("AF-"):
        kennung = kennung[3:].split("-")[0]
    if _PDB_ID.match(kennung):
        kennung = kennung.upper()
        pfad = os.path.join(ordner, f"{kennung}.pdb")
        if not os.path.exists(pfad):
            try:
                text = _hole(f"https://files.rcsb.org/download/{kennung}.pdb")
            except urllib.error.HTTPError as fehler:
                if fehler.code == 404:
                    raise ValueError(f"{kennung}: keine Datei im PDB-Format. Entweder gibt es die ID nicht, oder die "
                                     f"Struktur ist so groß, dass es sie nur als mmCIF gibt (wird hier nicht "
                                     f"unterstützt).") from None
                raise
            with open(pfad, "w") as f:
                f.write(text)
        with open(pfad) as f:
            return f.read(), f"PDB {kennung} (rcsb.org)"
    if _UNIPROT.match(kennung.upper()):
        kennung = kennung.upper()
        pfad = os.path.join(ordner, f"AF-{kennung}.pdb")
        if not os.path.exists(pfad):
            try:
                info = json.loads(_hole(f"https://alphafold.ebi.ac.uk/api/prediction/{kennung}"))
            except urllib.error.HTTPError as fehler:
                if fehler.code == 404:
                    raise ValueError(f"Für {kennung} gibt es kein AlphaFold-Modell (oder die UniProt-Nummer "
                                     f"stimmt nicht).") from None
                raise
            eintrag = info[0] if isinstance(info, list) else info
            with open(pfad, "w") as f:
                f.write(_hole(eintrag["pdbUrl"]))
        with open(pfad) as f:
            return f.read(), f"AlphaFold-Modell zu UniProt {kennung} (alphafold.ebi.ac.uk)"
    raise ValueError(f"„{eingabe}“ ist weder eine PDB-ID (4 Zeichen, z. B. 4HHB), noch eine UniProt-Nummer "
                     f"(z. B. P68871), noch eine vorhandene Datei.")


# ----------------------------------------------------------------------
# Lesen
# ----------------------------------------------------------------------

def _element(zeile, atomname):
    el = zeile[76:78].strip().upper() if len(zeile) >= 78 else ""
    if el:
        return el
    name = re.sub(r"[^A-Za-z]", "", atomname).upper()
    return name[:2] if name[:2] in RADIEN and len(name) > 1 and atomname[0] != " " else name[:1]


def lese_pdb(text):
    """Zerlegt eine PDB-Datei (nur das erste Modell). Gibt ein Dictionary S zurück:
    atome (Listen und numpy-Koordinaten), reste (je Aminosäure/Ligand), ketten (Sequenzen), liganden, Kopfdaten."""
    S = {"titel": "", "methode": "", "aufloesung": None, "alphafold": False}
    titel, modres = [], {}
    namen, resn, kette, resi, icode, xyz, bfak, elem, het = [], [], [], [], [], [], [], [], []
    gesehen = {}
    modelle = 0
    for zeile in text.splitlines():
        art = zeile[:6]
        if art == "TITLE ":
            titel.append(zeile[10:80].strip())
        elif art == "EXPDTA":
            S["methode"] = zeile[10:79].strip()
        elif art.startswith("REMARK") and zeile[6:10] == "   2" and "RESOLUTION." in zeile:
            treffer = re.search(r"(\d+\.\d+)", zeile)
            if treffer:
                S["aufloesung"] = float(treffer.group(1))
        elif art == "MODRES":
            modres[zeile[12:15].strip()] = zeile[24:27].strip()
        elif art == "MODEL ":
            modelle += 1
            if modelle > 1:
                break
        elif art == "ENDMDL":
            break
        elif art in ("ATOM  ", "HETATM"):
            name = zeile[12:16]
            schluessel = (zeile[21], zeile[22:27], name)
            besetzung = float(zeile[54:60]) if zeile[54:60].strip() else 1.0
            if schluessel in gesehen:          # alternative Positionen: die mit der höchsten Besetzung nehmen
                i, alt = gesehen[schluessel]
                if besetzung > alt:
                    xyz[i] = (float(zeile[30:38]), float(zeile[38:46]), float(zeile[46:54]))
                    bfak[i] = float(zeile[60:66]) if zeile[60:66].strip() else 0.0
                    gesehen[schluessel] = (i, besetzung)
                continue
            gesehen[schluessel] = (len(namen), besetzung)
            namen.append(name.strip())
            resn.append(zeile[17:20].strip())
            kette.append(zeile[21].strip() or "_")
            resi.append(int(zeile[22:26]))
            icode.append(zeile[26].strip())
            xyz.append((float(zeile[30:38]), float(zeile[38:46]), float(zeile[46:54])))
            bfak.append(float(zeile[60:66]) if zeile[60:66].strip() else 0.0)
            elem.append(_element(zeile, name))
            het.append(art == "HETATM")
    S["titel"] = re.sub(r"\s+", " ", " ".join(titel)).strip().capitalize()
    S["alphafold"] = "ALPHAFOLD" in text[:3000].upper()
    if not namen:
        raise ValueError("In der Datei stehen keine Atome (ATOM/HETATM-Zeilen).")
    S["atom_name"], S["atom_element"] = namen, elem
    S["xyz"] = np.array(xyz, dtype=float)
    S["b"] = np.array(bfak, dtype=float)
    # Reste bilden
    reste, index = [], {}
    for i in range(len(namen)):
        k = (kette[i], resi[i], icode[i], resn[i])
        if k not in index:
            index[k] = len(reste)
            n = resn[i]
            buchstabe = DREI_ZU_EINS.get(n) or VERAENDERT.get(n) or DREI_ZU_EINS.get(modres.get(n, ""))
            if n in WASSER:
                typ = "wasser"
            elif buchstabe and (not het[i] or n in VERAENDERT or n in modres):
                typ = "aminosaeure"
            else:
                typ = "ligand"
            reste.append({"kette": kette[i], "nr": resi[i], "icode": icode[i], "name": n, "typ": typ,
                          "buchstabe": buchstabe or "X", "atome": []})
        reste[index[k]]["atome"].append(i)
    for r in reste:
        r["atome"] = np.array(r["atome"])
        r["kurz"] = f"{r['nr']}{r['icode']}"
    S["reste"] = reste
    ketten = {}
    for j, r in enumerate(reste):
        if r["typ"] == "aminosaeure":
            ketten.setdefault(r["kette"], []).append(j)
    S["ketten"] = {k: {"reste": v, "sequenz": "".join(reste[j]["buchstabe"] for j in v)}
                   for k, v in ketten.items() if len(v) >= 5}
    S["liganden"] = [j for j, r in enumerate(reste) if r["typ"] == "ligand"]
    return S


def rest_text(r, mit_kette=True):
    """z. B. "B:His92" für Rest 92 (Histidin) in Kette B."""
    from werkzeuge import DREIBUCHSTABEN
    name = DREIBUCHSTABEN.get(r["buchstabe"], r["name"].title()) if r["typ"] == "aminosaeure" else r["name"]
    return (f"{r['kette']}:" if mit_kette else "") + f"{name}{r['kurz']}"


def ligand_name(code):
    if code in LIGAND_NAMEN:
        return f"{LIGAND_NAMEN[code]} ({code})"
    if code in ZUSAETZE:
        return f"{code} (vermutlich aus der Kristallisation)"
    return code


# ----------------------------------------------------------------------
# Geometrie
# ----------------------------------------------------------------------

def _schwer(S, atome):
    """Nur Schweratome (ohne Wasserstoff)."""
    return np.array([a for a in atome if S["atom_element"][a] not in ("H", "D")], dtype=int)


def _kugelpunkte(n):
    """n gleichmäßig verteilte Punkte auf der Einheitskugel (goldene Spirale, wie Biopython)."""
    k = np.arange(n)
    z = 1 - 1.0 / n - k * 2.0 / n
    r = np.sqrt(1 - z * z)
    laenge = k * math.pi * (3 - 5 ** 0.5)
    return np.column_stack([np.cos(laenge) * r, np.sin(laenge) * r, z])


def zugaenglichkeit(S, n_punkte=100, sonde=1.4):
    """Lösungsmittel-zugängliche Oberfläche (Shrake & Rupley 1973) für jeden Rest, in Å² und in % des
    Maximums. Wasser und Wasserstoffe werden weggelassen, alle Ketten und Liganden zählen mit
    (also so, wie die Struktur in der Datei vorliegt). Schreibt r["flaeche"] und r["relativ"] in die Reste."""
    from scipy.spatial import cKDTree
    atome = np.array([a for r in S["reste"] if r["typ"] != "wasser" for a in r["atome"]
                      if S["atom_element"][a] not in ("H", "D")], dtype=int)
    xyz = S["xyz"][atome]
    radien = np.array([RADIEN.get(S["atom_element"][a], 2.0) for a in atome]) + sonde
    kugel = _kugelpunkte(n_punkte)
    baum = cKDTree(xyz)
    flaeche = np.zeros(len(atome))
    nachbarn_alle = baum.query_ball_point(xyz, 2 * radien.max())
    for i, nachbarn in enumerate(nachbarn_alle):
        nachbarn = np.array([j for j in nachbarn if j != i], dtype=int)
        if len(nachbarn):
            d = np.sqrt(((xyz[nachbarn] - xyz[i]) ** 2).sum(1))
            nachbarn = nachbarn[d < radien[i] + radien[nachbarn]]      # nur überlappende Kugeln
        punkte = xyz[i] + kugel * radien[i]
        if len(nachbarn):
            # Punkt ist verdeckt, wenn er in der Kugel eines Nachbarn liegt
            d2 = ((punkte[:, None, :] - xyz[nachbarn][None, :, :]) ** 2).sum(-1)
            frei = (d2 > radien[nachbarn][None, :] ** 2).all(axis=1)
        else:
            frei = np.ones(n_punkte, dtype=bool)
        flaeche[i] = frei.sum() / n_punkte * 4 * math.pi * radien[i] ** 2
    pro_atom = dict(zip(atome.tolist(), flaeche))
    for r in S["reste"]:
        r["flaeche"] = float(sum(pro_atom.get(int(a), 0.0) for a in r["atome"]))
        maximum = MAX_FLAECHE.get(r["buchstabe"]) if r["typ"] == "aminosaeure" else None
        r["relativ"] = min(100.0, 100 * r["flaeche"] / maximum) if maximum else None
    return S


def lage(relativ):
    """Einordnung nach relativer Zugänglichkeit (Faustregel)."""
    if relativ is None:
        return "–"
    if relativ < VERGRABEN_BIS:
        return "vergraben (im Inneren)"
    if relativ < OBERFLAECHE_AB:
        return "teilweise frei"
    return "an der Oberfläche"


def abstand(S, rest_a, rest_b):
    """Kleinster Abstand (Å) zwischen den Schweratomen zweier Reste."""
    a = S["xyz"][_schwer(S, rest_a["atome"])]
    b = S["xyz"][_schwer(S, rest_b["atome"])]
    if not len(a) or not len(b):
        return float("inf")
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def kontakte(S, rest_index, grenze=4.0, nur_aminosaeuren=True):
    """Alle Reste mit einem Schweratom höchstens grenze Å von rest_index entfernt.
    Gibt [(index, abstand)] nach Abstand sortiert zurück."""
    ziel = S["reste"][rest_index]
    mitte = S["xyz"][ziel["atome"]].mean(axis=0)
    ausdehnung = np.sqrt(((S["xyz"][ziel["atome"]] - mitte) ** 2).sum(1)).max()
    liste = []
    for j, r in enumerate(S["reste"]):
        if j == rest_index or r["typ"] == "wasser" or (nur_aminosaeuren and r["typ"] != "aminosaeure"):
            continue
        # grobe Vorauswahl über den Abstand der Mittelpunkte
        if np.sqrt(((S["xyz"][r["atome"]].mean(axis=0) - mitte) ** 2).sum()) > ausdehnung + grenze + 12:
            continue
        d = abstand(S, ziel, r)
        if d <= grenze:
            liste.append((j, d))
    return sorted(liste, key=lambda t: t[1])


# ----------------------------------------------------------------------
# Protein-Positionen den Resten der Struktur zuordnen
# ----------------------------------------------------------------------

def ordne_protein_zu(S, protein):
    """Richtet das Protein an jeder Kette aus (BLOSUM62, lokal) und nimmt die beste Kette.
    Gibt {"kette", "zuordnung": {position_im_protein: rest_index}, "identitaet", "gleiche_ketten"} zurück."""
    protein = protein.upper().rstrip("*")
    beste = None
    for k, info in S["ketten"].items():
        erg = ausrichten_protein(protein, info["sequenz"], "lokal")
        if beste is None or erg["score"] > beste[1]["score"]:
            beste = (k, erg)
    if beste is None:
        raise ValueError("Die Struktur enthält keine Proteinkette.")
    k, erg = beste
    zuordnung, gleich, verglichen = {}, 0, 0
    p1, p2 = erg["start1"], erg["start2"]
    for x, y in zip(erg["a1"], erg["a2"]):
        if x != "-" and y != "-":
            zuordnung[p1] = S["ketten"][k]["reste"][p2 - 1]
            verglichen += 1
            gleich += x == y
        if x != "-":
            p1 += 1
        if y != "-":
            p2 += 1
    seq = S["ketten"][k]["sequenz"]
    return {"kette": k, "zuordnung": zuordnung, "identitaet": 100 * gleich / max(1, verglichen),
            "abgedeckt": verglichen, "gleiche_ketten": [a for a, i in S["ketten"].items()
                                                         if a != k and i["sequenz"] == seq]}


def lies_positionen(text):
    """"7", "E7V", "Glu7Val", "p.Glu7Val, 64, H93" → [(7, "E", "V"), (64, None, None), (93, "H", None)]."""
    from werkzeuge import DREIBUCHSTABEN
    eins = {v.upper(): k for k, v in DREIBUCHSTABEN.items()}
    liste = []
    for teil in re.split(r"[,;\s]+", text.strip()):
        teil = teil.strip()
        if teil.lower().startswith("p."):
            teil = teil[2:]
        if not teil:
            continue
        treffer = re.fullmatch(r"([A-Za-z]{0,3})(\d+)([A-Za-z]{0,3}|=)", teil)
        if not treffer:
            raise ValueError(f"„{teil}“ verstehe ich nicht. Beispiele: 7, E7V, Glu7Val")

        def buchstabe(t):
            if not t or t == "=":
                return None
            return eins.get(t.upper()) if len(t) == 3 else t.upper()
        liste.append((int(treffer.group(2)), buchstabe(treffer.group(1)), buchstabe(treffer.group(3))))
    return liste


# ----------------------------------------------------------------------
# UniProt
# ----------------------------------------------------------------------

UNIPROT_ARTEN = {"Active site": "Aktive Stelle", "Binding site": "Bindestelle", "Metal binding": "Metall-Bindung",
                 "Site": "Wichtige Stelle", "Disulfide bond": "Disulfidbrücke"}


def uniprot_merkmale(nummer):
    """Liest Sequenz und bekannte Stellen (aktive Stelle, Bindestellen …) aus UniProt.
    Gibt {"name", "sequenz", "merkmale": [{art, start, ende, text}]} zurück."""
    daten = json.loads(_hole(f"https://rest.uniprot.org/uniprotkb/{nummer.strip().upper()}.json"))
    name = ""
    try:
        name = daten["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        pass
    merkmale = []
    for m in daten.get("features", []):
        art = UNIPROT_ARTEN.get(m.get("type"))
        if not art or (m.get("description") or "").startswith("(Microbial infection)"):
            continue
        try:
            start = int(m["location"]["start"]["value"])
            ende = int(m["location"]["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        teile = [(m.get("ligand") or {}).get("name", "")]
        if (m.get("ligandPart") or {}).get("name"):
            teile[0] += f" ({m['ligandPart']['name']})"
        teile.append(m.get("description") or "")
        merkmale.append({"art": art, "start": start, "ende": ende, "text": " – ".join(t for t in teile if t)})
    varianten = []
    for m in daten.get("features", []):
        if m.get("type") != "Natural variant":
            continue
        try:
            start = int(m["location"]["start"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        alt = m.get("alternativeSequence") or {}
        varianten.append({"position": start, "ref": alt.get("originalSequence", ""),
                          "alt": "/".join(alt.get("alternativeSequences") or []) or "fehlt",
                          "text": re.sub(r";\s*dbSNP:\S+", "", m.get("description") or "")})
    return {"name": name, "sequenz": (daten.get("sequence") or {}).get("value", ""), "merkmale": merkmale,
            "varianten": varianten}


# ----------------------------------------------------------------------
# Bewertung einer Position
# ----------------------------------------------------------------------

def bewerte_position(S, j, neu=None, ligand_grenze=4.0):
    """Lage eines Restes: Zugänglichkeit, nächster Ligand, Kontakte zu anderen Ketten, Hinweise."""
    from werkzeuge import EIGENSCHAFT
    r = S["reste"][j]
    ergebnis = {"rest": r, "relativ": r.get("relativ"), "lage": lage(r.get("relativ"))}
    naechster = None
    for L in S["liganden"]:
        if S["reste"][L]["name"] in ZUSAETZE:
            continue
        d = abstand(S, r, S["reste"][L])
        if naechster is None or d < naechster[1]:
            naechster = (L, d)
    ergebnis["ligand"] = naechster
    andere = sorted({S["reste"][k]["kette"] for k, _d in kontakte(S, j, 4.0) if S["reste"][k]["kette"] != r["kette"]})
    ergebnis["andere_ketten"] = andere
    hinweise = []
    alt = r["buchstabe"]
    if neu and neu != alt:
        e_alt, e_neu = EIGENSCHAFT.get(alt, "?"), EIGENSCHAFT.get(neu, "?")
        punkte = BLOSUM62.get(alt, {}).get(neu, -4)
        if r.get("relativ") is not None and r["relativ"] < VERGRABEN_BIS:
            if e_alt == "unpolar" and e_neu != "unpolar":
                hinweise.append("Im Inneren wird eine wasserabweisende durch eine polare oder geladene Aminosäure "
                                "ersetzt. Das stört die Faltung meist stark (wahrscheinlich destabilisierend).")
            elif neu == "P" or alt == "G":
                hinweise.append("Im Inneren: Prolin oder der Verlust von Glycin verändert die Form der Kette, "
                                "das stört die Faltung oft.")
            elif punkte < 0:
                hinweise.append("Seltener Austausch im Inneren des Proteins: eher destabilisierend.")
            else:
                hinweise.append("Ähnlicher Austausch im Inneren: oft gut verträglich, kann aber die Packung stören.")
        elif r.get("relativ") is not None and r["relativ"] >= OBERFLAECHE_AB:
            text = "An der Oberfläche ist die Faltung meist wenig betroffen."
            if e_alt != "unpolar" and e_neu == "unpolar":
                text += (" Eine neue wasserabweisende Stelle an der Oberfläche kann aber Proteine aneinander "
                         "kleben lassen. Genau das passiert bei HbS (Glu → Val): Die Moleküle bilden Fasern.")
            hinweise.append(text)
    if naechster and naechster[1] <= ligand_grenze:
        hinweise.append(f"Berührt direkt {ligand_name(S['reste'][naechster[0]]['name'])}: "
                        f"eine Änderung hier kann die Bindung verändern.")
    if andere:
        hinweise.append(f"Liegt an der Kontaktfläche zu Kette {', '.join(andere)}: kann das Zusammenspiel der "
                        f"Ketten beeinflussen.")
    if S["alphafold"]:
        plddt = float(S["b"][r["atome"]].mean())
        ergebnis["plddt"] = plddt
        if plddt < 70:
            hinweise.append(f"AlphaFold ist sich hier unsicher (pLDDT {plddt:.0f}). Die Lage ist nur grob.")
    ergebnis["hinweise"] = hinweise
    return ergebnis
