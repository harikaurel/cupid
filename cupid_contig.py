#!/usr/bin/env python3


from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import Counter
import sys
import csv
import json
import colorsys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


csv.field_size_limit(sys.maxsize)


MIN_OVERLAP_MOTIFS_DEFAULT = 1
ATOL_TIE_DEFAULT = 0


CARBAPENEMASE_REGEX_DEFAULT = r"KPC|OXA-48|OXA-244|OXA-181|OXA-232|NDM|VIM|IMP|GES"

TAXID_RE = re.compile(r"\(taxid\s+(\d+)\)")


MAIN_RANK_ORDER = {
    "R": 0,
    "D": 1,
    "K": 2,
    "P": 3,
    "C": 4,
    "O": 5,
    "F": 6,
    "G": 7,
    "S": 8,
}


def find_one(folder: Path, patterns: list[str], label: str) -> Path:
    hits: list[Path] = []
    for pat in patterns:
        hits.extend(folder.glob(pat))
    hits = [h for h in hits if h.is_file()]
    if not hits:
        raise SystemExit(f"[STOP] Could not find {label} in {folder}")
    if len(hits) > 1:
        print(f"[WARN] Multiple {label} files found, using {hits[0].name}")
    return hits[0]


def parse_taxid_from_kraken_field(field: str) -> str:
    if not isinstance(field, str):
        return ""
    s = field.strip()
    m = TAXID_RE.search(s)
    if m:
        return m.group(1)
    if s.isdigit():
        return s
    return ""


def clean_taxon(label: str) -> str:
    if not isinstance(label, str):
        return "Unclassified"
    s = label.strip()
    if not s or s.lower() in ("nan", "unclassified"):
        return "Unclassified"
    s = re.sub(r"\([^)]*\)", "", s)
    s = re.sub(r"\s+", " ", s).strip()

    if not s or s.lower() == "unclassified":
        return "Unclassified"
    return s


def canonical_species_label(taxon: str) -> str:
    s = clean_taxon(taxon)
    if s == "Unclassified":
        return s
    s = re.sub(r"^\s*(s__|g__|p__|k__|c__|o__|f__)\s*", "", s).strip()
    parts = s.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return s


def full_taxon_label(taxon: str) -> str:
    s = clean_taxon(taxon)
    if s == "Unclassified":
        return s
    s = re.sub(r"^\s*(s__|g__|p__|k__|c__|o__|f__)\s*", "", s).strip()
    return s


def rank_depth(rank_code: str):
    if not rank_code:
        return None
    m = re.match(r"^([A-Z])(\d*)$", rank_code)
    if not m:
        return None
    main, num = m.group(1), m.group(2)
    if main not in MAIN_RANK_ORDER:
        return None
    return (MAIN_RANK_ORDER[main], int(num) if num else 0)


def load_kraken_out(path: Path) -> dict[str, dict]:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] < 3:
        raise ValueError(f"Kraken OUT has <3 columns: {path}")

    out = {}
    contigs = df.iloc[:, 1].str.strip()
    raws = df.iloc[:, 2].fillna("Unclassified").str.strip()

    for contig, raw in zip(contigs, raws):
        raw = raw if raw else "Unclassified"
        taxid = parse_taxid_from_kraken_field(raw)
        out[contig] = {
            "raw": raw,
            "taxid": taxid,
            "species_can": canonical_species_label(raw),
            "full_label": full_taxon_label(raw),
        }
    return out


def load_kraken_report(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    taxid_to_rank: dict[str, str] = {}
    taxid_to_genus: dict[str, str] = {}
    stack: list[tuple[tuple, str, str, str]] = []

    with path.open("r", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                parts = re.split(r"\s+", line.strip(), maxsplit=5)
                if len(parts) < 6:
                    continue
            rank_code = parts[3].strip()
            taxid = parts[4].strip()
            name = parts[5].strip()

            if not taxid:
                continue

            taxid_to_rank[taxid] = rank_code

            depth = rank_depth(rank_code)
            if depth is None:

                taxid_to_genus[taxid] = ""
                continue


            while stack and stack[-1][0] >= depth:
                stack.pop()


            genus_name = ""
            for d, rc, tid, nm in reversed(stack):
                if rc == "G":
                    genus_name = nm
                    break
            taxid_to_genus[taxid] = genus_name

            stack.append((depth, rank_code, taxid, name))

    return taxid_to_rank, taxid_to_genus


def read_table_skip_hash(path: Path) -> pd.DataFrame:
    skip = 0
    with path.open("r", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str)


def read_amr_map(amr_file: Path) -> dict[str, set[str]]:
    df = read_table_skip_hash(amr_file)
    contig_col = next(c for c in df.columns if "contig" in c.lower())
    gene_col = next(c for c in df.columns if ("element" in c.lower()) or ("gene" in c.lower()))

    out: dict[str, set[str]] = {}
    for ctg, gene in zip(df[contig_col].astype(str).str.strip(), df[gene_col].astype(str).str.strip()):
        if ctg and gene and ctg.lower() != "nan":
            out.setdefault(ctg, set()).add(gene)
    return out


def extract_contig_types(contig_report: Path) -> pd.DataFrame:
    rep = pd.read_csv(contig_report, sep="\t", dtype=str)
    rep["contig_id"] = rep["contig_id"].astype(str).str.strip()
    rep["molecule_type"] = rep["molecule_type"].astype(str).str.lower().str.strip()
    return rep[["contig_id", "molecule_type"]].drop_duplicates()


def load_motif_scores_long(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "median" not in df.columns:

        df["median"] = pd.to_numeric(df["methylation_value"], errors="coerce")
        if df["median"].isna().all():
            raise ValueError(f"'median' column not found in {path}. Columns: {list(df.columns)}")

    df["median"] = pd.to_numeric(df["median"], errors="coerce")
    df["contig"] = df["contig"].astype(str).str.strip()
    df["motif"] = df["motif"].astype(str).str.strip()
    df["mod_type"] = df["mod_type"].astype(str).str.strip()
    df["mod_position"] = df["mod_position"].astype(str).str.strip()
    df["motif"] = df["motif"] + "_" + df["mod_type"] + "_" + df["mod_position"]

    df = df.groupby(["contig", "motif"], as_index=False)["median"].median()
    df = df.dropna(subset=["median"])
    return df


def build_contig_arrays(df_long: pd.DataFrame):
    contig2motifs: dict[str, np.ndarray] = {}
    contig2vals: dict[str, np.ndarray] = {}

    for contig, sub in df_long.groupby("contig", sort=False):
        m = sub["motif"].to_numpy(dtype=object)
        v = sub["median"].to_numpy(dtype=float)
        order = np.argsort(m)
        contig2motifs[contig] = m[order]
        contig2vals[contig] = v[order]

    return contig2motifs, contig2vals


def rmsd_nm_shared(contigA: str, contigB: str, contig2motifs, contig2vals) -> tuple[float, int]:
    ma = contig2motifs.get(contigA)
    mb = contig2motifs.get(contigB)
    if ma is None or mb is None:
        return np.nan, 0

    shared, ia, ib = np.intersect1d(ma, mb, assume_unique=False, return_indices=True)
    nm = int(shared.size)
    if nm == 0:
        return np.nan, 0

    va = contig2vals[contigA][ia]
    vb = contig2vals[contigB][ib]
    rmsd = float(np.sqrt(np.mean((va - vb) ** 2)))
    return rmsd, nm


def score_candidates(
    query_id: str,
    candidate_ids: list[str],
    contig2motifs,
    contig2vals,
    kraken_out: dict[str, dict],
    min_overlap: int,
):
    rows = []
    for c in candidate_ids:
        if c == query_id:
            continue

        rmsd, nm = rmsd_nm_shared(query_id, c, contig2motifs, contig2vals)
        if nm < min_overlap or not np.isfinite(rmsd):
            continue

        rmss = max(0.0, 1.0 - rmsd)
        final_score = rmss * nm

        kinfo = kraken_out.get(c, {})
        raw = kinfo.get("raw", "Unclassified")
        species_can = canonical_species_label(raw)
        full_label = kinfo.get("full_label", full_taxon_label(raw))
        taxid = kinfo.get("taxid", "")

        rows.append((c, rmsd, nm, rmss, final_score, species_can, full_label, taxid))

    df = pd.DataFrame(
        rows,
        columns=["candidate", "rmsd", "shared_motifs", "rmss", "final_score",
                 "species", "full_label", "taxid"],
    )
    if df.empty:
        return df

    df = df.sort_values(["final_score", "shared_motifs", "rmsd"], ascending=[False, False, True]).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


def tie_species_summary_one_line(df_at_max: pd.DataFrame) -> str:
    counts = df_at_max["species"].astype(str).value_counts()
    return "; ".join([f"{sp} ({int(n)})" for sp, n in counts.items()])


def pick_winner_species_from_top_ties(top_df: pd.DataFrame) -> str:
    counts = Counter(top_df["species"].astype(str).tolist())
    return counts.most_common(1)[0][0] if counts else "Unclassified"


def is_valid_candidate(c: str, kraken_out: dict) -> bool:
    raw = kraken_out.get(c, {}).get("raw", "")
    if not raw:
        return False
    return clean_taxon(raw) != "Unclassified"




GENUS_BASE_COLORS = {
    "Klebsiella":      "#e6194b",
    "Escherichia":     "#4363d8",
    "Citrobacter":     "#3cb44b",
    "Enterobacter":    "#f58231",
    "Morganella":      "#17becf",
    "Corynebacterium": "#911eb4",
    "Enterococcus":    "#c724b1",
    "Streptomyces":    "#bcbd22",


    "Pseudomonas":     "#16a085",
    "Acinetobacter":   "#8e6f3e",
    "Psychrobacter":   "#566573",
    "Aeromonas":       "#884ea0",
    "Pectobacterium":  "#d4a017",
    "Staphylococcus":  "#795548",
}


FIXED_SPECIES_COLORS = {

    "Klebsiella pneumoniae":               "#e6194b",
    "Klebsiella oxytoca":                  "#ecb6c8",
    "Klebsiella michiganensis":            "#8b048b",
    "Klebsiella variicola":                "#f3b9c2",
    "Klebsiella aerogenes":                "#52000d",
    "Klebsiella quasipneumoniae":          "#eb48d0",
    "Klebsiella":                          "#a32941",
    "Klebsiella sp.":                      "#f60aab",

    "Escherichia coli":                    "#4363d8",
    "Escherichia":                         "#b6c7ec",
    "Escherichia sp.":                     "#2cd4e0",
    "Escherichia albertii":              "#742cd8",

    "Citrobacter freundii":                "#3cb44b",
    "Citrobacter arsenatis":               "#b6ecb8",
    "Citrobacter braakii":                 "#D0D38E",
    "Citrobacter europaeus":               "#deed56",
    "Citrobacter koseri":                  "#005212",
    "Citrobacter portucalensis":           "#047904",
    "Citrobacter youngae":                 "#29a395",
    "Citrobacter":                         "#45bc8b",
    "Citrobacter sp.":                     "#8aeaca",

    "Enterobacter hormaechei":             "#f58231",
    "Enterobacter cancerogenus":           "#ecc8b6",
    "Enterobacter cloacae":                "#8b4404",
    "Enterobacter ludwigii":               "#deb9a5",
    "Enterobacter sp.":                    "#522900",
    "Enterobacter":                        "#c1855a",

    "Morganella morganii":                 "#17cf1d",

    "Aeromonas hydrophila":                "#881965",
    "Aeromonas allosaccharophila":         "#d7b6ec",
    "Aeromonas caviae":                    "#6c048b",
    "Aeromonas media":                     "#bd63ee",
    "Aeromonas veronii":                   "#0E06B2",
    "Aeromonas sp.":                       "#eadef2",
    "Aeromonas":                           "#e008ce",

    "Pseudomonas aeruginosa":              "#16a085",
    "Pseudomonas putida":                  "#b6ecdd",

    "Acinetobacter baumannii":             "#8e6f3e",
    "Acinetobacter johnsonii":             "#ecd2b6",

    "Psychrobacter sp.":                   "#566573",

    "Pectobacterium parvum":               "#d4a017",

    "Staphylococcus saprophyticus":        "#795548",

    "Corynebacterium tuberculostearicum":  "#911eb4",
    "Corynebacterium jeikeium":            "#dbb6ec",
    "Corynebacterium":                     "#431CB9",

    "Enterococcus faecium":                "#c724b1",
    "Enterococcus avium":                  "#b6dfec",

    "Streptomyces paludis":                "#bcbd22",

    "cellular organisms":                  "#6ea9da",
    "Enterobacteriaceae":                  "#ffee00",
    "root":                                "#9e9e9e",
    "Gammaproteobacteria":                 "#D44C0D",
    "Caudoviricetes sp.":                 "#7acdc1",
}


def _hex_to_rgb01(hx: str):
    hx = hx.lstrip("#")
    return tuple(int(hx[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _rgb01_to_hex(r, g, b) -> str:
    return "#%02x%02x%02x" % (
        int(round(max(0.0, min(1.0, r)) * 255)),
        int(round(max(0.0, min(1.0, g)) * 255)),
        int(round(max(0.0, min(1.0, b)) * 255)),
    )


def _genus_of(species: str) -> str:
    s = (species or "").strip()
    return s.split()[0] if s else ""


_GENUS_SHADE_STEPS = [
    (0.000, None, None),
    (-0.015, 0.82, 0.58),
    (+0.010, 0.28, 0.95),
    (-0.010, 0.66, 0.80),
    (+0.015, 0.16, 1.00),
    (-0.015, 0.91, 0.45),
    (+0.008, 0.40, 0.60),
    (-0.008, 0.55, 0.92),
    (+0.012, 0.73, 0.70),
    (-0.006, 0.34, 0.78),
]


def _genus_base_hls(genus: str):
    anchor = GENUS_BASE_COLORS.get(genus)
    if anchor is None:
        import hashlib
        hh = (int(hashlib.md5(genus.encode("utf-8")).hexdigest(), 16) % 360) / 360.0
        return (hh, 0.50, 0.70)
    return colorsys.rgb_to_hls(*_hex_to_rgb01(anchor))


def _genus_shade(genus: str, idx: int) -> str:
    h0, l0, s0 = _genus_base_hls(genus)
    dh, L, S = _GENUS_SHADE_STEPS[idx % len(_GENUS_SHADE_STEPS)]
    if L is None:
        anchor = GENUS_BASE_COLORS.get(genus)
        if anchor is not None:
            return anchor
        L, S = l0, s0
    extra = idx // len(_GENUS_SHADE_STEPS)
    h = (h0 + dh) % 1.0
    L = max(0.0, min(1.0, L * (0.88 ** extra)))
    S = max(0.0, min(1.0, S))
    return _rgb01_to_hex(*colorsys.hls_to_rgb(h, L, S))

PCOA_GRAY = "#cfcfcf"


_ABUNDANCE_POS = {
    "upper left":    (0.02, 0.98, "left",   "top"),
    "upper right":   (0.98, 0.98, "right",  "top"),
    "upper center":  (0.50, 0.98, "center", "top"),
    "center left":   (0.02, 0.50, "left",   "center"),
    "center":        (0.50, 0.50, "center", "center"),
    "center right":  (0.98, 0.50, "right",  "center"),
    "lower left":    (0.02, 0.02, "left",   "bottom"),
    "lower right":   (0.98, 0.02, "right",  "bottom"),
    "lower center":  (0.50, 0.02, "center", "bottom"),
}


def _abundance_xy(loc: str):
    return _ABUNDANCE_POS.get((loc or "lower left").strip().lower(),
                              _ABUNDANCE_POS["lower left"])


_INSET_SIZE = 0.345
_INSET_MARGIN = 0.02
_INSET_FAR = 1.0 - _INSET_MARGIN - _INSET_SIZE
_INSET_MID = (1.0 - _INSET_SIZE) / 2.0
_INSET_POS = {
    "upper left":   [_INSET_MARGIN, _INSET_FAR, _INSET_SIZE, _INSET_SIZE],
    "upper right":  [_INSET_FAR,    _INSET_FAR, _INSET_SIZE, _INSET_SIZE],
    "upper center": [_INSET_MID,    _INSET_FAR, _INSET_SIZE, _INSET_SIZE],
    "lower left":   [_INSET_MARGIN, _INSET_MARGIN, _INSET_SIZE, _INSET_SIZE],
    "lower right":  [_INSET_FAR,    _INSET_MARGIN, _INSET_SIZE, _INSET_SIZE],
    "lower center": [_INSET_MID,    _INSET_MARGIN, _INSET_SIZE, _INSET_SIZE],
}


def _inset_bbox(loc: str):
    return _INSET_POS.get((loc or "upper right").strip().lower(),
                          _INSET_POS["upper right"])


_SPECIES_COLOR_REGISTRY: dict[str, str] = dict(FIXED_SPECIES_COLORS)


def load_species_color_registry(path) -> None:
    global _SPECIES_COLOR_REGISTRY
    reg = dict(FIXED_SPECIES_COLORS)
    if path and Path(path).is_file():
        try:
            with open(path) as f:
                saved = json.load(f)
            if isinstance(saved, dict):
                reg = {str(k): str(v) for k, v in saved.items()}
                reg.update(FIXED_SPECIES_COLORS)
                print(f"      -> loaded {len(reg)} species colours from {path}")
        except Exception as e:
            print(f"[color-map] could not read {path} ({e}); starting fresh")
    _SPECIES_COLOR_REGISTRY = reg


def save_species_color_registry(path) -> None:
    if not path:
        return
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(_SPECIES_COLOR_REGISTRY, f, indent=2, sort_keys=True)
        print(f"      -> saved {len(_SPECIES_COLOR_REGISTRY)} species colours to {path}")
    except Exception as e:
        print(f"[color-map] could not write {path} ({e})")


def color_for_species(species: str) -> str:
    if species in ("Unclassified", "other", "", None):
        return PCOA_GRAY
    if species in _SPECIES_COLOR_REGISTRY:
        return _SPECIES_COLOR_REGISTRY[species]
    genus = _genus_of(species)
    if not genus:
        return PCOA_GRAY


    same_genus = [sp for sp in _SPECIES_COLOR_REGISTRY if _genus_of(sp) == genus]
    used_in_genus = {_SPECIES_COLOR_REGISTRY[sp] for sp in same_genus}
    idx = len(same_genus)
    color = _genus_shade(genus, idx)

    guard = 0
    while color in used_in_genus and guard < 3 * len(_GENUS_SHADE_STEPS):
        idx += 1
        guard += 1
        color = _genus_shade(genus, idx)
    _SPECIES_COLOR_REGISTRY[species] = color
    return color


def _taxon_label_for_contig(contig, kraken_out):
    kinfo = kraken_out.get(contig, {})
    full = kinfo.get("full_label", "") or full_taxon_label(kinfo.get("raw", ""))
    if not full or full == "Unclassified":
        return "Unclassified"
    parts = full.split()


    if len(parts) >= 2 and parts[1][:1].islower():
        return f"{parts[0]} {parts[1]}"
    return full


def _backbone_finalscore_matrix(ids, contig2motifs, contig2vals, min_overlap):
    n = len(ids)
    FS = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            rmsd, nm = rmsd_nm_shared(ids[i], ids[j], contig2motifs, contig2vals)
            if nm < min_overlap or not np.isfinite(rmsd):
                continue
            FS[i, j] = FS[j, i] = max(0.0, 1.0 - rmsd) * nm
    fs_max = np.nanmax(FS) if np.any(np.isfinite(FS)) else 1.0
    np.fill_diagonal(FS, fs_max)
    return FS, fs_max


def _finalscore_to_distance(FS, fs_max):
    denom = fs_max if (fs_max and fs_max > 0) else 1.0
    D = 1.0 - (FS / denom)
    D[~np.isfinite(D)] = 1.0
    D = np.clip(D, 0.0, 1.0)
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    return D


def _project_supplementary(d_to_backbone, backbone_coords, backbone_D, eigvals):
    D2 = backbone_D ** 2
    row_means = D2.mean(axis=1)
    d2 = np.asarray(d_to_backbone, dtype=float) ** 2
    safe_lambda = np.where(eigvals > 0, eigvals, np.nan)
    U = backbone_coords / np.sqrt(safe_lambda)
    contrib = (row_means - d2)[:, None] * U
    f = np.nansum(contrib, axis=0) / (2.0 * np.sqrt(safe_lambda))
    return f


def _canonicalize_axis_signs(coords):
    out = np.asarray(coords, dtype=float).copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        if col.size == 0:
            continue
        k = int(np.argmax(np.abs(col)))
        if col[k] < 0:
            out[:, j] = -col
    return out


def _ordinate_chroms(chrom_list, contig2motifs, contig2vals, min_overlap):
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa as skbio_pcoa
    n = len(chrom_list)
    if n < 3:
        return None, None, None, None
    FS, fs_max = _backbone_finalscore_matrix(
        chrom_list, contig2motifs, contig2vals, min_overlap)
    D = _finalscore_to_distance(FS, fs_max)
    denom = fs_max if (fs_max and fs_max > 0) else 1.0
    try:
        dm = DistanceMatrix(D, ids=list(chrom_list))
        ndim = min(10, n - 1)
        try:
            ordr = skbio_pcoa(dm, number_of_dimensions=ndim)
        except TypeError:
            ordr = skbio_pcoa(dm)
        coords = ordr.samples.to_numpy()
        eig = np.asarray(ordr.eigvals)
        keep = min(ndim, coords.shape[1])
        coords = _canonicalize_axis_signs(coords[:, :keep])
        return coords, eig[:keep], denom, D
    except Exception as e:
        print(f"[PCoA] candidate ordination failed ({e}); inset skipped.")
        return None, None, None, None


def _pcoa_topN_for_query(
    q_id: str,
    mol: str,
    q_motif_vector_ids: np.ndarray,
    q_motif_vector_vals: np.ndarray,
    all_chroms: list,
    contig2motifs: dict,
    contig2vals: dict,
    fs_to_all: np.ndarray,
    fs_max: float,
    bb_taxa_all: list,
    line_n: int,
    top_n: int,
    color_of_global: dict,
    out_png: Path,
    out_pdf: Path,
    legend_loc: str = "lower right",
    legend_fontsize: float = 12.0,
    full_coords2: np.ndarray | None = None,
    full_query_xy: tuple | None = None,
    topN_species_set_full: set | None = None,
    full_taxa: list | None = None,
    draw_inset: bool = True,
    inset_redbox: bool = True,
    inset_loc: str = "upper right",
    abundance_loc: str = "lower left",
):
    try:
        from skbio import DistanceMatrix
        from skbio.stats.ordination import pcoa as skbio_pcoa
    except ImportError:
        return None, set(), []

    n_all = len(all_chroms)
    if n_all < 3:
        return None, set(), []


    fs_arr = np.asarray(fs_to_all)
    order_by_score = np.argsort(-fs_arr, kind="stable")
    sel = [int(i) for i in order_by_score[:top_n] if fs_arr[i] > 0]
    if len(sel) < 3:
        return None, set(), []
    sub_chroms = [all_chroms[i] for i in sel]
    sub_taxa = [bb_taxa_all[i] for i in sel]
    sub_fs = fs_arr[sel]


    n_sub = len(sub_chroms)
    FS_sub = np.zeros((n_sub, n_sub), dtype=float)
    for i in range(n_sub):
        for j in range(i + 1, n_sub):
            rmsd, nm = rmsd_nm_shared(sub_chroms[i], sub_chroms[j],
                                       contig2motifs, contig2vals)
            if nm < 1 or not np.isfinite(rmsd):
                continue
            FS_sub[i, j] = FS_sub[j, i] = max(0.0, 1.0 - rmsd) * nm
    fs_max_sub = float(np.nanmax(FS_sub)) if np.any(np.isfinite(FS_sub)) else 1.0
    np.fill_diagonal(FS_sub, fs_max_sub)
    denom = fs_max_sub if fs_max_sub > 0 else 1.0
    D_sub = 1.0 - (FS_sub / denom)
    D_sub[~np.isfinite(D_sub)] = 1.0
    D_sub = np.clip(D_sub, 0.0, 1.0)
    D_sub = (D_sub + D_sub.T) / 2.0
    np.fill_diagonal(D_sub, 0.0)


    dm = DistanceMatrix(D_sub, ids=sub_chroms)
    ndim = min(10, n_sub - 1)
    try:
        ord_res = skbio_pcoa(dm, number_of_dimensions=ndim)
    except TypeError:
        ord_res = skbio_pcoa(dm)
    samples = ord_res.samples
    coords_full = samples.to_numpy()
    eigvals = np.asarray(ord_res.eigvals)
    ve = np.asarray(ord_res.proportion_explained) * 100
    keep_dims = min(ndim, coords_full.shape[1])
    coords_full = _canonicalize_axis_signs(coords_full[:, :keep_dims])
    eigvals = eigvals[:keep_dims]
    ve = ve[:keep_dims]
    coords2 = coords_full[:, :2]


    d_to_sub = []
    for c in sub_chroms:
        rmsd, nm = rmsd_nm_shared(q_id, c, contig2motifs, contig2vals)
        if nm < 1 or not np.isfinite(rmsd):
            d_to_sub.append(1.0)
        else:
            s = max(0.0, 1.0 - rmsd) * nm
            d_to_sub.append(float(np.clip(1.0 - s / denom, 0.0, 1.0)))
    f = _project_supplementary(np.array(d_to_sub), coords_full, D_sub, eigvals)
    px, py = float(f[0]), float(f[1])


    top_score_species = None
    for k in np.argsort(-sub_fs, kind="stable"):
        t = sub_taxa[int(k)]
        if t not in ("Unclassified", "other"):
            top_score_species = t
            break


    N_LEGEND_TAXA = 5
    best_score_by_sp: dict[str, float] = {}
    for k in range(n_sub):
        t = sub_taxa[k]
        if t in ("Unclassified", "other"):
            continue
        best_score_by_sp[t] = max(best_score_by_sp.get(t, -np.inf), float(sub_fs[k]))
    ranked_species = sorted(best_score_by_sp, key=lambda s: -best_score_by_sp[s])
    if top_score_species in ranked_species:
        ranked_species = [top_score_species] + [s for s in ranked_species
                                                if s != top_score_species]
    top_species = ranked_species[:N_LEGEND_TAXA]
    top_species_set = set(top_species)


    topN_species_set = top_species_set
    point_color = [color_for_species(sub_taxa[k]) if sub_taxa[k] in top_species_set
                   else PCOA_GRAY for k in range(n_sub)]
    has_unclassified = any(point_color[k] == PCOA_GRAY for k in range(n_sub))


    ring_set = {int(np.argmax(sub_fs))} if len(sub_fs) else set()


    LABEL_FS = 22
    TICK_FS = 16
    LEGEND_FS = legend_fontsize

    fig, ax = plt.subplots(figsize=(16, 11))

    for k in range(n_sub):
        if k in ring_set:
            continue
        ax.scatter(coords2[k, 0], coords2[k, 1],
                   c=[point_color[k]],
                   s=90, edgecolor="white", linewidth=0.5,
                   alpha=1.0, zorder=2)

    for k in sorted(ring_set):
        ax.scatter(coords2[k, 0], coords2[k, 1],
                   c=[point_color[k]],
                   s=90, edgecolor="black", linewidth=1.8,
                   alpha=1.0, zorder=5)
    ax.scatter(px, py, marker="*", s=480, facecolor="none",
               edgecolor="black", linewidth=2.2, zorder=6)
    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)
    ax.set_xlabel(f"PCo1 ({ve[0]:.1f}%)", fontsize=LABEL_FS)
    ax.set_ylabel(f"PCo2 ({ve[1]:.1f}%)", fontsize=LABEL_FS)
    ax.tick_params(axis="both", labelsize=TICK_FS)
    ax.set_box_aspect(1)


    abund_src = full_taxa if full_taxa is not None else bb_taxa_all
    pop = [t for t in abund_src if t not in ("Unclassified", "other")]
    abundance_text = None
    if top_score_species is not None and pop:
        n_chr = len(pop)
        n_top_bg = sum(1 for t in pop if t == top_score_species)
        pct_bg = 100.0 * n_top_bg / n_chr if n_chr else 0.0
        abundance_text = f"top taxon abundance: {pct_bg:.1f}%\n(of {n_chr} contigs)"


    inset_drawn = False
    inset_bbox = None
    if draw_inset and full_coords2 is not None and full_query_xy is not None:
        fc2 = np.asarray(full_coords2)


        itaxa = full_taxa if full_taxa is not None else bb_taxa_all
        col_set = topN_species_set_full if topN_species_set_full is not None else topN_species_set
        qxf, qyf = float(full_query_xy[0]), float(full_query_xy[1])
        inset_bbox = _inset_bbox(inset_loc)
        axin = ax.inset_axes(inset_bbox)
        inset_drawn = True
        axin.set_facecolor("white")
        fcol = [color_for_species(itaxa[k]) if itaxa[k] in col_set
                else PCOA_GRAY for k in range(len(itaxa))]
        g_idx = [k for k in range(len(fcol)) if fcol[k] == PCOA_GRAY]
        c_idx = [k for k in range(len(fcol)) if fcol[k] != PCOA_GRAY]
        if g_idx:
            axin.scatter(fc2[g_idx, 0], fc2[g_idx, 1], c=PCOA_GRAY, s=6,
                         linewidths=0, alpha=0.75, zorder=1)
        if c_idx:
            axin.scatter(fc2[c_idx, 0], fc2[c_idx, 1],
                         c=[fcol[k] for k in c_idx], s=9,
                         linewidths=0, alpha=0.95, zorder=2)


        if inset_redbox:
            d2q = np.hypot(fc2[:, 0] - qxf, fc2[:, 1] - qyf)
            n_box = min(top_n, len(fc2))
            box_idx = [int(i) for i in np.argsort(d2q, kind="stable")[:n_box]]
        else:
            box_idx = []
        if box_idx:
            xs = np.append(fc2[box_idx, 0], qxf)
            ys = np.append(fc2[box_idx, 1], qyf)
            bx0, bx1 = float(np.min(xs)), float(np.max(xs))
            by0, by1 = float(np.min(ys)), float(np.max(ys))
            bpx = 0.03 * (bx1 - bx0) if bx1 > bx0 else 0.01
            bpy = 0.03 * (by1 - by0) if by1 > by0 else 0.01
            axin.add_patch(Rectangle((bx0 - bpx, by0 - bpy),
                                     (bx1 - bx0) + 2 * bpx, (by1 - by0) + 2 * bpy,
                                     fill=False, edgecolor="red", linewidth=1.6,
                                     zorder=5))
        axin.scatter([qxf], [qyf], marker="*", s=80, facecolor="none",
                     edgecolor="black", linewidth=1.2, zorder=6)
        axin.set_xticks([])
        axin.set_yticks([])
        for sp in axin.spines.values():
            sp.set_edgecolor("0.5")
            sp.set_linewidth(0.8)


    if abundance_text is not None:
        abx = dict(boxstyle="round,pad=0.4", facecolor="white",
                   alpha=0.6, edgecolor="0.6", linewidth=0.8)
        if inset_drawn and inset_bbox is not None:
            ix0, iy0, iw, ih = inset_bbox
            icx = ix0 + iw / 2.0
            if iy0 > 0.15:
                ty, tva = iy0 - 0.02, "top"
            else:
                ty, tva = iy0 + ih + 0.02, "bottom"
            ax.text(icx, ty, abundance_text, transform=ax.transAxes,
                    ha="center", va=tva, fontsize=LEGEND_FS, color="#666666",
                    zorder=8, bbox=abx)
        else:
            ax_x, ax_y, aha, ava = _abundance_xy(abundance_loc)
            ax.text(ax_x, ax_y, abundance_text, transform=ax.transAxes,
                    ha=aha, va=ava, fontsize=LEGEND_FS, color="#666666",
                    zorder=8, bbox=abx)


    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor=color_for_species(t), markeredgecolor="none",
                   markersize=11, alpha=1.0, label=t)
        for t in top_species
    ]
    n_taxa_entries = len(handles)
    if has_unclassified:
        handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor=PCOA_GRAY, markeredgecolor="none",
                                  markersize=11, alpha=1.0, label="other"))
    handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                              markerfacecolor="white", markeredgecolor="black",
                              markersize=11, label="top hit"))
    handles.append(plt.Line2D([0], [0], marker="*", linestyle="",
                              markerfacecolor="black", markeredgecolor="white",
                              markersize=18, label=f"query {mol}"))

    legend_loc_l = (legend_loc or "").strip().lower()
    show_legend = legend_loc_l not in ("none", "off", "hidden", "")
    right_margin = 0.97
    if show_legend:
        if legend_loc_l == "outside":
            leg = ax.legend(handles=handles, loc="center left",
                            bbox_to_anchor=(1.02, 0.5),
                            frameon=False, fontsize=LEGEND_FS)
            right_margin = 0.70
        else:
            leg = ax.legend(handles=handles, loc=legend_loc_l,
                            frameon=True, framealpha=0.6, fontsize=LEGEND_FS)
            leg.get_frame().set_facecolor("white")
            leg.get_frame().set_edgecolor("0.6")
            leg.get_frame().set_linewidth(0.8)
        for txt in leg.get_texts()[:n_taxa_entries]:
            txt.set_fontstyle("italic")

    fig.subplots_adjust(left=0.07, right=right_margin, top=0.97, bottom=0.10)
    out_png.parent.mkdir(parents=True, exist_ok=True)


    fig.canvas.draw()
    _bbox = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.1)
    fig.savefig(out_png, dpi=300, bbox_inches=_bbox)
    fig.savefig(out_pdf, bbox_inches=_bbox)
    plt.close(fig)


    return top_score_species, topN_species_set, sel


_NOBS_CANDIDATE_COLS = (
    "n_motif_obs", "N_motif_obs", "n_motif", "N_motif",
    "n_mod", "N_mod", "motif_obs", "n_obs", "N_obs", "count",
)


def build_nobs_lookup(meth_path: Path) -> dict:
    try:
        df = pd.read_csv(meth_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"[nobs] could not read {meth_path} for counts: {e}")
        return {}
    df.columns = [c.strip() for c in df.columns]
    ncol = next((c for c in _NOBS_CANDIDATE_COLS if c in df.columns), None)
    if ncol is None:
        low = {c.lower(): c for c in df.columns}
        ncol = next((low[k] for k in low
                     if ("obs" in k and ("motif" in k or k[:1] == "n"))
                     or ("motif" in k and "count" in k)), None)
    if ncol is None:
        print("[nobs] no observation-count column found; cells will show 'n=?'.")
        return {}
    for c in ("contig", "motif", "mod_type", "mod_position"):
        if c not in df.columns:
            print(f"[nobs] missing '{c}'; cannot build count lookup.")
            return {}
    contig = df["contig"].astype(str).str.strip()
    key = (df["motif"].astype(str).str.strip() + "_"
           + df["mod_type"].astype(str).str.strip() + "_"
           + df["mod_position"].astype(str).str.strip())
    n = pd.to_numeric(df[ncol], errors="coerce")
    g = pd.DataFrame({"contig": contig, "key": key, "n": n}).dropna(subset=["n"])
    g = g.groupby(["contig", "key"], as_index=False)["n"].sum()
    print(f"[nobs] using '{ncol}' as the observation count.")
    return {(c, k): float(v) for c, k, v in zip(g["contig"], g["key"], g["n"])}


def _italic_taxon(s: str) -> str:
    esc = str(s).replace("\\", " ").replace("$", "")
    for ch in ("_", "^", "{", "}"):
        esc = esc.replace(ch, " ")
    esc = esc.replace(" ", r"\ ")
    return "$\\it{" + esc + "}$"


def select_top_distinct_taxa(cand_df: pd.DataFrame, n_taxa: int,
                             kraken_rank: dict) -> list:
    picked: list = []
    seen: set = set()
    for _, row in cand_df.iterrows():
        sp = str(row["species"]).strip()
        if sp in ("", "Unclassified", "other"):
            continue

        taxid = str(row.get("taxid", "")).strip()
        rank_code = kraken_rank.get(taxid, "") if taxid else ""
        if not rank_code.startswith("S"):
            continue
        if sp in seen:
            continue
        seen.add(sp)
        picked.append({
            "contig": str(row["candidate"]),
            "species": sp,
            "full_label": str(row.get("full_label", sp)),
            "taxid": str(row.get("taxid", "")),
            "rmsd": float(row["rmsd"]),
            "shared": int(row["shared_motifs"]),
            "final_score": float(row["final_score"]),
        })
        if len(picked) >= n_taxa:
            break
    return picked


def plot_top_taxa_heatmap(
    query_id: str,
    query_mol: str,
    picked: list,
    contig2motifs,
    contig2vals,
    out_png: Path,
    out_pdf: Path,
    nobs: dict | None = None,
    all_motifs: bool = False,
    max_motifs: int = 0,
    order: str = "methylation",
    orient: str = "motifs-rows",
    title: str | None = None,
    fontsize: float = 11.0,
):
    from matplotlib.colors import LinearSegmentedColormap

    nobs = nobs or {}

    def _meth_dict(cid):
        ms = contig2motifs.get(cid)
        vs = contig2vals.get(cid)
        if ms is None or vs is None:
            return {}
        return dict(zip(ms.tolist(), vs.tolist()))

    q_meth = _meth_dict(query_id)
    chrom_meth = {r["contig"]: _meth_dict(r["contig"]) for r in picked}

    chrom_union = set().union(*[set(d) for d in chrom_meth.values()]) if chrom_meth else set()


    keys = (set(q_meth) | chrom_union) if all_motifs else set(q_meth)
    if not keys:
        print(f"[heatmap] {query_id}: no motifs to plot "
              f"({'union empty' if all_motifs else 'plasmid has no motifs'}); skipping.")
        return
    keys = list(keys)

    if order == "name":
        keys.sort()
    elif order == "nobs":
        def _row_nobs(k):
            tot = nobs.get((query_id, k), 0.0) or 0.0
            for r in picked:
                tot += nobs.get((r["contig"], k), 0.0) or 0.0
            return -tot
        keys.sort(key=_row_nobs)
    else:
        def _row_mean(k):
            vals = [q_meth.get(k)] + [chrom_meth[r["contig"]].get(k) for r in picked]
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            return -(np.mean(vals) if vals else -1.0)
        keys.sort(key=_row_mean)
    if max_motifs and len(keys) > max_motifs:


        top_meth = chrom_meth.get(picked[0]["contig"], {}) if picked else {}
        have = [k for k in keys if k in top_meth]
        if len(have) >= max_motifs:
            selected = set(have[:max_motifs])
        else:
            missing = [k for k in keys if k not in top_meth]
            selected = set(have) | set(missing[:max_motifs - len(have)])
        keys = [k for k in keys if k in selected]

    col_contigs = [query_id] + [r["contig"] for r in picked]
    col_meth = [q_meth] + [chrom_meth[r["contig"]] for r in picked]
    ncols, nrows = len(col_contigs), len(keys)

    M_val = np.full((nrows, ncols), np.nan)
    N_val = np.full((nrows, ncols), np.nan)
    for i, k in enumerate(keys):
        for j, (cid, md) in enumerate(zip(col_contigs, col_meth)):
            if k in md:
                M_val[i, j] = md[k]
            nn = nobs.get((cid, k))
            if nn is not None:
                N_val[i, j] = nn

    q_tag = "plasmid query" if query_mol == "plasmid" else f"{query_mol} query"
    col_labels = [f"{query_id}\n({q_tag})"]
    for idx, r in enumerate(picked):


        is_top = (idx == 0)
        sp = r["species"]
        lab = f"{r['contig']}"
        if is_top:
            if sp and sp != "Unclassified":
                lab += f"\n{sp}"
            elif sp == "Unclassified":
                lab += "\nUnclassified"
        else:
            if sp and sp != "Unclassified":
                lab += "\n" + _italic_taxon(sp)
            elif sp == "Unclassified":
                lab += "\nUnclassified"
        col_labels.append(lab)

    fs = fontsize
    motif_labels = list(keys)
    cmap = LinearSegmentedColormap.from_list(
        "meth", ["#f7fbff", "#3182bd", "#08306b"])
    cmap.set_bad("#e5e7eb")

    if orient == "contigs-rows":
        disp, dispN = M_val.T, N_val.T
        fig_w = max(8.0, 0.95 * nrows + 3.0)
        fig_h = max(3.6, 0.62 * ncols + 1.6)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(np.ma.masked_invalid(disp), cmap=cmap, vmin=0.0, vmax=1.0,
                       aspect="auto")
        ax.set_yticks(range(ncols))
        ax.set_yticklabels(col_labels, fontsize=fs * 0.9)
        if ncols >= 2:
            _tl = ax.get_yticklabels()[1]
            _tl.set_fontweight("bold"); _tl.set_fontstyle("italic")
        ax.set_xticks(range(nrows))
        ax.set_xticklabels(motif_labels, fontsize=fs * 0.78, family="monospace",
                           rotation=45, ha="right")
        ax.axhline(0.5, color="0.25", lw=1.2)
        gr, gc = ncols, nrows
    else:
        disp, dispN = M_val, N_val
        fig_h = max(4.2, 0.42 * nrows + 2.0)
        fig_w = max(7.0, 1.7 * ncols + 1.5)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(np.ma.masked_invalid(disp), cmap=cmap, vmin=0.0, vmax=1.0,
                       aspect="auto")
        ax.set_xticks(range(ncols))
        ax.set_xticklabels(col_labels, fontsize=fs * 0.7)
        if ncols >= 2:
            _tl = ax.get_xticklabels()[1]
            _tl.set_fontweight("bold"); _tl.set_fontstyle("italic")
        ax.set_yticks(range(nrows))
        ax.set_yticklabels(motif_labels, fontsize=fs, family="monospace")
        ax.axvline(0.5, color="0.25", lw=1.2)
        gr, gc = nrows, ncols

    ax.tick_params(top=False, bottom=False, left=False, right=False)

    for i in range(gr):
        for j in range(gc):
            v = disp[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center",
                        fontsize=fs * 0.8, color="#6b7280")
                continue
            tcol = "white" if v > 0.55 else "#1f2937"
            ax.text(j, i - 0.16, f"{v:.2f}", ha="center", va="center",
                    fontsize=fs * 0.95, color=tcol)
            nn = dispN[i, j]
            ntxt = f"n={int(nn)}" if np.isfinite(nn) else "n=?"
            ax.text(j, i + 0.24, ntxt, ha="center", va="center",
                    fontsize=fs * 0.72, color=tcol, alpha=0.85)

    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("methylation value", fontsize=fs)
    cbar.ax.tick_params(labelsize=fs * 0.85)
    if title:
        ax.set_title(title, fontsize=fs * 1.1, pad=10)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    species_list = ", ".join((r["species"] or "?") for r in picked)
    print(f"[heatmap] {query_id} ({query_mol}): {nrows} motifs x {ncols} cols "
          f"({'union' if all_motifs else 'all-plasmid-motifs'}); "
          f"hosts: {species_list} -> {out_png}")


def run_top_taxa_heatmaps(
    contig_types: pd.DataFrame,
    amr_map: dict,
    kraken_out: dict,
    kraken_rank: dict,
    contig2motifs,
    contig2vals,
    outdir: Path,
    score_path: Path,
    min_overlap_motifs: int = 1,
    n_taxa: int = 5,
    all_motifs: bool = False,
    max_motifs: int = 0,
    order: str = "methylation",
    orient: str = "motifs-rows",
    fontsize: float = 11.0,
    only_contig: set | None = None,
    prefix: str = "",
):
    if only_contig is None:
        print("      -> no --only-contig: skipping top-taxa heatmaps "
              "(figures are only written for requested contigs).")
        return

    pfx = (f"{prefix}." if prefix else "")
    base = outdir / "top_taxa_heatmap" / "per_contig"

    def _safe_name(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))


    chrom_ids = [
        c for c in contig_types.query("molecule_type == 'chromosome'")["contig_id"].astype(str)
        if c in contig2motifs
    ]
    if not chrom_ids:
        print("[heatmap] no chromosome candidates with motif vectors; skipping.")
        return

    mol_of = dict(zip(contig_types["contig_id"].astype(str),
                      contig_types["molecule_type"].astype(str)))


    queries = [ctg for ctg in contig_types["contig_id"].astype(str)
               if ctg in contig2motifs and amr_map.get(ctg)]
    queries = [q for q in queries if q in only_contig]
    missing = sorted(only_contig - set(queries))
    if missing:
        print(f"[heatmap] NOTE: {len(missing)} requested contig(s) are not AMR "
              f"queries with motifs and will be skipped: {missing}")
    if not queries:
        print("[heatmap] nothing to do after --only-contig filter.")
        return

    nobs = build_nobs_lookup(score_path)

    for q in queries:
        mol = mol_of.get(q, "unknown")
        sub = "plasmids" if mol == "plasmid" else "chromosomes"
        cands = [c for c in chrom_ids if c != q]
        cand_df = score_candidates(
            query_id=q, candidate_ids=cands,
            contig2motifs=contig2motifs, contig2vals=contig2vals,
            kraken_out=kraken_out, min_overlap=int(min_overlap_motifs),
        )
        if cand_df.empty:
            print(f"[heatmap] {q}: no chromosome shared >= {min_overlap_motifs} "
                  f"motif(s); skipping.")
            continue
        picked = select_top_distinct_taxa(cand_df, int(n_taxa), kraken_rank)
        if not picked:
            print(f"[heatmap] {q}: no species-level taxa among candidates "
                  f"(genus-only / unclassified); skipping.")
            continue
        stem = f"{pfx}{_safe_name(q)}.{mol}.top{n_taxa}_taxa_heatmap"
        plot_top_taxa_heatmap(
            query_id=q, query_mol=mol, picked=picked,
            contig2motifs=contig2motifs, contig2vals=contig2vals,
            out_png=base / sub / f"{stem}.png",
            out_pdf=base / sub / f"{stem}.pdf",
            nobs=nobs, all_motifs=all_motifs, max_motifs=max_motifs,
            order=order, orient=orient, fontsize=fontsize,
        )


def run_pcoa(
    contig_types: pd.DataFrame,
    amr_map: dict[str, set[str]],
    kraken_out: dict[str, dict],
    kraken_rank: dict[str, str],
    kraken_genus: dict[str, str],
    contig2motifs,
    contig2vals,
    outdir: Path,
    min_overlap_motifs: int,
    carbapenemase_regex: str,
    max_backbone: int = 0,
    color_k: int = 12,
    line_n: int = 10,
    weighted_alpha: float = 0.5,
    weighted_beta: float = 0.5,
    topn_pcoa: int = 50,
    only_contig: set | None = None,
    legend_loc: str = "lower right",
    legend_fontsize: float = 12.0,
    color_map_path=None,
    abundance_loc: str = "lower left",
    inset_loc: str = "upper right",
    inset_topn: int = 0,
    prefix: str = "",
):


    if only_contig is None:
        print("      -> no --only-contig: skipping PCoA stage entirely "
              "(no figures and no neighbours TSV are needed).")
        return

    try:
        from skbio import DistanceMatrix
        from skbio.stats.ordination import pcoa
    except ImportError:
        print("[PCoA] scikit-bio not installed; skipping PCoA stage. "
              "Install with: pip install scikit-bio")
        return


    pcoa_dir = outdir / "pcoa"


    chroms = [
        c for c in contig_types.query("molecule_type == 'chromosome'")["contig_id"].astype(str)
        if c in contig2motifs and is_valid_candidate(c, kraken_out)
    ]
    if max_backbone and len(chroms) > max_backbone:
        chroms = chroms[:max_backbone]
    if len(chroms) < 3:
        print(f"[PCoA] Only {len(chroms)} backbone chromosomes (<3); skipping PCoA.")
        return
    print(f"      -> PCoA backbone: {len(chroms)} chromosomes")


    carb_re = re.compile(carbapenemase_regex, re.IGNORECASE)
    mol_of = dict(zip(contig_types["contig_id"].astype(str),
                      contig_types["molecule_type"].astype(str)))
    queries = []
    for ctg in contig_types["contig_id"].astype(str):
        if ctg not in contig2motifs:
            continue
        genes = sorted(amr_map.get(ctg, set()))
        if not genes:
            continue
        carb_hits = sorted(g for g in genes if carb_re.search(g))
        queries.append((ctg, mol_of.get(ctg, "unknown"), genes, carb_hits))
    print(f"      -> PCoA queries: {len(queries)} AMR contigs "
          f"({sum(1 for _,m,_,_ in queries if m=='plasmid')} plasmid, "
          f"{sum(1 for _,m,_,_ in queries if m=='chromosome')} chromosome; "
          f"{sum(1 for _,_,_,ch in queries if ch)} carry a carbapenemase)")


    if only_contig is not None:
        before_n = len(queries)
        queries = [q for q in queries if q[0] in only_contig]
        missing = sorted(only_contig - {q[0] for q in queries})
        print(f"      -> --only-contig filter: keeping {len(queries)} / {before_n} queries")
        if missing:
            print(f"      -> NOTE: {len(missing)} requested contig(s) not in AMR query set "
                  f"and will be skipped: {missing}")
        if not queries:
            print("      -> nothing to do after --only-contig filter; "
                  "no figures will be written.")
            return


    topN_dir = pcoa_dir.parent / f"pcoa_top{topn_pcoa}" / "per_contig"
    (topN_dir / "plasmids").mkdir(parents=True, exist_ok=True)
    (topN_dir / "chromosomes").mkdir(parents=True, exist_ok=True)


    taxa = [_taxon_label_for_contig(c, kraken_out) for c in chroms]
    LINE_N = line_n

    pfx = (f"{prefix}." if prefix else "")


    load_species_color_registry(color_map_path)


    def _safe_name(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))

    for q, mol, all_genes, carb_hits in queries:
        sub = "plasmids" if mol == "plasmid" else "chromosomes"


        topN_stem = f"{pfx}{_safe_name(q)}.{mol}.pcoa_top{topn_pcoa}"


        fs_to_all_q = []
        for c in chroms:
            if c == q:
                fs_to_all_q.append(0.0)
                continue
            rmsd_q, nm_q = rmsd_nm_shared(q, c, contig2motifs, contig2vals)
            if nm_q < min_overlap_motifs or not np.isfinite(rmsd_q):
                fs_to_all_q.append(0.0)
            else:
                fs_to_all_q.append(max(0.0, 1.0 - rmsd_q) * nm_q)


        cand_idx    = [i for i, fs in enumerate(fs_to_all_q) if fs > 0.0]


        if inset_topn and inset_topn > 0 and len(cand_idx) > inset_topn:
            cand_idx = sorted(cand_idx, key=lambda i: fs_to_all_q[i],
                              reverse=True)[:inset_topn]
        cand_chroms = [chroms[i] for i in cand_idx]
        cand_taxa   = [taxa[i] for i in cand_idx]

        cand_coords, cand_eig, cand_denom, cand_D = _ordinate_chroms(
            cand_chroms, contig2motifs, contig2vals, min_overlap_motifs)


        qxy_full = None
        if cand_coords is not None:
            d_to_cand = []
            for c in cand_chroms:
                rmsd_b, nm_b = rmsd_nm_shared(q, c, contig2motifs, contig2vals)
                if nm_b < min_overlap_motifs or not np.isfinite(rmsd_b):
                    d_to_cand.append(1.0)
                else:
                    s_b = max(0.0, 1.0 - rmsd_b) * nm_b
                    d_to_cand.append(float(np.clip(1.0 - s_b / cand_denom, 0.0, 1.0)))
            fproj = _project_supplementary(np.array(d_to_cand), cand_coords,
                                           cand_D, cand_eig)
            qxy_full = (float(fproj[0]), float(fproj[1]))
        full_coords2_arg = cand_coords[:, :2] if cand_coords is not None else None


        top_score_species, topN_species_set, topN_indices = _pcoa_topN_for_query(
            q_id=q, mol=mol,
            q_motif_vector_ids=contig2motifs.get(q),
            q_motif_vector_vals=contig2vals.get(q),
            all_chroms=chroms,
            contig2motifs=contig2motifs,
            contig2vals=contig2vals,
            fs_to_all=np.asarray(fs_to_all_q),
            fs_max=1.0,
            bb_taxa_all=taxa,
            line_n=LINE_N,
            top_n=topn_pcoa,
            color_of_global={},
            out_png=topN_dir / sub / f"{topN_stem}.png",
            out_pdf=topN_dir / sub / f"{topN_stem}.pdf",
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            full_coords2=full_coords2_arg,
            full_query_xy=qxy_full,
            full_taxa=cand_taxa,
            draw_inset=True,
            inset_redbox=True,
            inset_loc=inset_loc,
            abundance_loc=abundance_loc,
        )

    print(f"      -> wrote {len(queries)} per-contig top-{topn_pcoa} figure(s) under "
          f"pcoa_top{topn_pcoa}/per_contig/{{plasmids,chromosomes}}/")


    save_species_color_registry(color_map_path)


def run_assignment(
    contig_types: pd.DataFrame,
    amr_map: dict[str, set[str]],
    kraken_out: dict[str, dict],
    contig2motifs,
    contig2vals,
    outdir: Path,
    min_overlap_motifs: int,
    atol_tie: float,
    prefix: str = "",
):
    outdir.mkdir(parents=True, exist_ok=True)


    chroms_all = [
        c for c in contig_types.query("molecule_type == 'chromosome'")["contig_id"].astype(str)
        if c in contig2motifs and is_valid_candidate(c, kraken_out)
    ]
    chroms_vec = [c for c in chroms_all if c in contig2motifs]


    plasmids_all = [
        c for c in contig_types.query("molecule_type == 'plasmid'")["contig_id"].astype(str)
        if c in contig2motifs and amr_map.get(c)
    ]
    chromosomes_all = [
        c for c in contig_types.query("molecule_type == 'chromosome'")["contig_id"].astype(str)
        if c in contig2motifs and amr_map.get(c)
    ]

    print(f"[ASSIGN] candidate chromosome pool: {len(chroms_vec)} contigs with motif vectors")
    print(f"[ASSIGN] AMR queries: {len(plasmids_all)} plasmids, {len(chromosomes_all)} chromosomes total")

    def base_row(query_id: str, mol: str) -> dict:
        genes = sorted(amr_map.get(query_id, set()))
        return {
            "query_contig": query_id,
            "molecule_type": mol,
            "amr_genes": ";".join(genes),
            "n_amr_genes": len(genes),
        }

    plasmid_rows = []


    for p in plasmids_all:
        row = base_row(p, "plasmid")

        cand_df = score_candidates(
            query_id=p,
            candidate_ids=chroms_vec,
            contig2motifs=contig2motifs,
            contig2vals=contig2vals,
            kraken_out=kraken_out,
            min_overlap=min_overlap_motifs,
        )

        if cand_df.empty:
            row.update({
                "best_host_contig": "",
                "best_host_species": "Unclassified",
                "best_host_full_label": "Unclassified",
                "best_final_score": np.nan,
                "best_rmsd": np.nan,
                "best_shared_motifs": 0,
                "n_candidates": 0,
                "winner_species_top_ties": "Unclassified",
                "tie_species_summary": "",
            })
            plasmid_rows.append(row)
            continue

        best = cand_df.iloc[0]
        max_score = float(best["final_score"])
        at_max = cand_df[np.isclose(cand_df["final_score"], max_score, atol=atol_tie)]
        winner_species = pick_winner_species_from_top_ties(at_max)

        row.update({
            "best_host_contig": str(best["candidate"]),
            "best_host_species": str(best["species"]),
            "best_host_full_label": str(best["full_label"]),
            "best_final_score": max_score,
            "best_rmsd": float(best["rmsd"]),
            "best_shared_motifs": int(best["shared_motifs"]),
            "n_candidates": int(len(cand_df)),
            "winner_species_top_ties": winner_species,
            "tie_species_summary": tie_species_summary_one_line(at_max),
        })
        plasmid_rows.append(row)


    plasmid_df = pd.DataFrame(plasmid_rows)

    pfx = (f"{prefix}." if prefix else "")
    plasmid_df.to_csv(outdir / f"{pfx}AMR_plasmid_assignment.tsv", sep="\t", index=False)

    print("[DONE] Outputs written to:", str(outdir))
    print(f"  - {pfx}AMR_plasmid_assignment.tsv")


def main():
    ap = argparse.ArgumentParser(
        description="Assign AMR-carrying plasmids to host chromosomes by "
                    "Nanomotif methylation similarity, and render PCoA / "
                    "heatmap figures.")
    ap.add_argument("--nanomotif-dir", required=True,
                    help="Nanomotif output dir (contains motifs-scored-read-methylation.tsv)")
    ap.add_argument("--mobsuite-dir", required=True,
                    help="MobSuite output dir (contains contig_report.txt)")
    ap.add_argument("--amr-dir", required=True,
                    help="AMRFinderPlus output dir")
    ap.add_argument("--kraken-dir", required=True,
                    help="Kraken2 output dir (contains the *.out and *.report/*.kreport files)")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for the assignment tables and figures")

    ap.add_argument("--min-overlap-motifs", type=int, default=MIN_OVERLAP_MOTIFS_DEFAULT,
                    help=f"Minimum shared motifs to score a pair (default {MIN_OVERLAP_MOTIFS_DEFAULT})")
    ap.add_argument("--atol-tie", type=float, default=ATOL_TIE_DEFAULT,
                    help=f"Absolute tolerance for treating final_scores as tied (default {ATOL_TIE_DEFAULT})")
    ap.add_argument("--carbapenemase-regex", default=CARBAPENEMASE_REGEX_DEFAULT,
                    help="Regex flagging carbapenemase genes (annotation only)")
    ap.add_argument("--prefix", default="",
                    help="Optional filename prefix (e.g. sample id) for all outputs")


    ap.add_argument("--pcoa-max-backbone", type=int, default=0,
                    help="Cap the chromosome backbone size (0 = no cap)")
    ap.add_argument("--topn-pcoa", type=int, default=100,
                    help="Number of top candidates for the focused per-contig PCoA (default 100)")
    ap.add_argument("--only-contig", default=None,
                    help="Comma-separated contig ids: only render figures for these "
                         "(assignment tables are always computed for all)")
    ap.add_argument("--legend-loc", default="lower right",
                    help="Legend location for PCoA figures ('outside'/'none' allowed)")
    ap.add_argument("--legend-fontsize", type=float, default=12.0,
                    help="Legend font size for PCoA figures")
    ap.add_argument("--abundance-loc", default="lower left",
                    help="Location of the relative-abundance annotation box")
    ap.add_argument("--inset-loc", default="upper right",
                    help="Corner for the full-ordination inset")
    ap.add_argument("--inset-topn", type=int, default=0,
                    help="Cap the inset's candidate set to the top-N by score (0 = all)")
    ap.add_argument("--color-map", default=None,
                    help="Path to a persistent species->colour JSON (loaded and updated)")


    ap.add_argument("--top-taxa-n", type=int, default=5,
                    help="Number of distinct species-level taxa in the heatmap")
    ap.add_argument("--top-taxa-all-motifs", action="store_true",
                    help="Use the union of query+host motifs (default: query motifs only)")
    ap.add_argument("--top-taxa-max-motifs", type=int, default=0,
                    help="Cap the number of motif rows in the heatmap (0 = no cap)")
    ap.add_argument("--top-taxa-order", default="methylation",
                    choices=["methylation", "name", "nobs"],
                    help="Row ordering for the heatmap")
    ap.add_argument("--top-taxa-orient", default="motifs-rows",
                    choices=["motifs-rows", "contigs-rows"],
                    help="Heatmap orientation")
    ap.add_argument("--top-taxa-fontsize", type=float, default=11.0,
                    help="Base font size for the heatmap")

    args = ap.parse_args()

    nanomotif_dir = Path(args.nanomotif_dir)
    mobsuite_dir = Path(args.mobsuite_dir)
    amr_dir = Path(args.amr_dir)
    kraken_dir = Path(args.kraken_dir)
    outdir = Path(args.outdir)

    only_contig = None
    if args.only_contig:
        only_contig = {c.strip() for c in args.only_contig.split(",") if c.strip()}


    print("[1/8] Locating input files ...")
    motif_path = find_one(
        nanomotif_dir,
        ["*motifs-scored-read-methylation.tsv", "motifs-scored-read-methylation.tsv",
         "*motifs-scored*.tsv"],
        "Nanomotif motif-scored methylation TSV")
    contig_report = find_one(
        mobsuite_dir, ["contig_report.txt", "*contig_report.txt"],
        "MobSuite contig_report.txt")
    amr_file = find_one(
        amr_dir, ["*.tsv", "*amrfinder*.txt", "*.amrfinder.tsv"],
        "AMRFinderPlus results table")
    kraken_out_path = find_one(
        kraken_dir, ["*.out", "*.kraken", "*kraken*.out"],
        "Kraken2 per-contig OUT file")
    kraken_report_path = find_one(
        kraken_dir, ["*.report", "*.kreport", "*report*.txt"],
        "Kraken2 report file")


    print("[2/8] Loading Kraken2 classifications ...")
    kraken_out = load_kraken_out(kraken_out_path)
    kraken_rank, kraken_genus = load_kraken_report(kraken_report_path)

    print("[3/8] Loading AMRFinderPlus gene map ...")
    amr_map = read_amr_map(amr_file)

    print("[4/8] Loading MobSuite contig types ...")
    contig_types = extract_contig_types(contig_report)

    print("[5/8] Loading Nanomotif methylation and building contig arrays ...")
    df_long = load_motif_scores_long(motif_path)
    contig2motifs, contig2vals = build_contig_arrays(df_long)


    print("[6/8] Assigning AMR plasmids to host chromosomes ...")
    run_assignment(
        contig_types=contig_types,
        amr_map=amr_map,
        kraken_out=kraken_out,
        contig2motifs=contig2motifs,
        contig2vals=contig2vals,
        outdir=outdir,
        min_overlap_motifs=args.min_overlap_motifs,
        atol_tie=args.atol_tie,
        prefix=args.prefix,
    )


    print("[7/8] Rendering PCoA figures ...")
    run_pcoa(
        contig_types=contig_types,
        amr_map=amr_map,
        kraken_out=kraken_out,
        kraken_rank=kraken_rank,
        kraken_genus=kraken_genus,
        contig2motifs=contig2motifs,
        contig2vals=contig2vals,
        outdir=outdir,
        min_overlap_motifs=args.min_overlap_motifs,
        carbapenemase_regex=args.carbapenemase_regex,
        max_backbone=args.pcoa_max_backbone,
        topn_pcoa=args.topn_pcoa,
        only_contig=only_contig,
        legend_loc=args.legend_loc,
        legend_fontsize=args.legend_fontsize,
        color_map_path=args.color_map,
        abundance_loc=args.abundance_loc,
        inset_loc=args.inset_loc,
        inset_topn=args.inset_topn,
        prefix=args.prefix,
    )

    print("[8/8] Rendering top-taxa heatmaps ...")
    run_top_taxa_heatmaps(
        contig_types=contig_types,
        amr_map=amr_map,
        kraken_out=kraken_out,
        kraken_rank=kraken_rank,
        contig2motifs=contig2motifs,
        contig2vals=contig2vals,
        outdir=outdir,
        score_path=motif_path,
        min_overlap_motifs=args.min_overlap_motifs,
        n_taxa=args.top_taxa_n,
        all_motifs=args.top_taxa_all_motifs,
        max_motifs=args.top_taxa_max_motifs,
        order=args.top_taxa_order,
        orient=args.top_taxa_orient,
        fontsize=args.top_taxa_fontsize,
        only_contig=only_contig,
        prefix=args.prefix,
    )

    print("[ALL DONE]")


if __name__ == "__main__":
    main()