#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import re
import sys
import colorsys
from pathlib import Path
from multiprocessing import Pool
from collections import Counter

import numpy as np
import pandas as pd
from numba import njit

csv.field_size_limit(sys.maxsize)


MIN_OVERLAP_MOTIFS_DEFAULT = 1
ATOL_TIE_DEFAULT = 0


CARBAPENEMASE_REGEX_DEFAULT = r"KPC|OXA-48|OXA-244|OXA-181|OXA-232|NDM|VIM|IMP|GES"

TAXID_RE = re.compile(r"\(taxid\s+(\d+)\)")

OUTPUT_COLUMNS = [
    "read_id",
    "molecule_type",
    "amr_genes",
    "has_nanomotif_vector",
    "kraken_raw",
    "kraken_rank_code",
    "mean_top_score_species_counts",
]


_read2motifs = None
_read2mean = None
_kraken_out = None
_kraken_rank = None
_amr_map = None
_chroms_vec = None
_chrom_motif_sets = None
_min_overlap = None
_atol_tie = None


def _worker_init(r2motifs, r2mean, k_out, k_rank, a_map, c_vec, c_sets, min_ov, atol):
    global _read2motifs, _read2mean, _kraken_out, _kraken_rank
    global _amr_map, _chroms_vec, _chrom_motif_sets, _min_overlap, _atol_tie
    _read2motifs    = r2motifs
    _read2mean      = r2mean
    _kraken_out     = k_out
    _kraken_rank    = k_rank
    _amr_map        = a_map
    _chroms_vec     = c_vec
    _chrom_motif_sets = c_sets
    _min_overlap    = min_ov
    _atol_tie       = atol

    _warmup_numba()


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


def normalize_read_id(s: str) -> str:
    return re.split(r"[\s|]", s.strip())[0]


def _detect_sep(path) -> str:
    try:
        with open(path, "r", errors="replace") as f:
            first = f.readline()
    except Exception:
        return "\t"
    n_comma, n_tab = first.count(","), first.count("\t")
    if n_tab >= n_comma and n_tab > 0:
        return "\t"
    if n_comma > 0:
        return ","
    return "\t"


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


def load_kraken_out(path: Path) -> dict[str, dict]:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] < 3:
        raise ValueError(f"Kraken OUT has <3 columns: {path}")

    out = {}
    contigs = df.iloc[:, 1].astype(str).apply(normalize_read_id)
    raws = df.iloc[:, 2].fillna("Unclassified").str.strip()

    for contig, raw in zip(contigs, raws):
        raw = raw if raw else "Unclassified"
        taxid = parse_taxid_from_kraken_field(raw)
        out[contig] = {
            "raw": raw,
            "taxid": taxid,
            "species_can": canonical_species_label(raw),
        }

    print(f"[INFO] Kraken OUT: {len(out)} reads loaded")
    return out


def load_kraken_report_rank_map(path: Path) -> dict[str, str]:
    rep: dict[str, str] = {}
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
            if taxid:
                rep[taxid] = rank_code
    return rep


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

    if df.empty:
        print(f"[WARN] AMR file is empty: {amr_file}")
        return {}

    lower_to_orig = {c.lower(): c for c in df.columns}

    if "read_id" in df.columns:
        id_col = "read_id"
    elif "contig_id" in df.columns:
        id_col = "contig_id"
    elif "contig id" in lower_to_orig:
        id_col = lower_to_orig["contig id"]
    else:
        read_candidates = [c for c in df.columns if "read" in c.lower()]
        contig_candidates = [c for c in df.columns if "contig" in c.lower()]
        candidates = read_candidates + contig_candidates
        if not candidates:
            raise ValueError(
                f"Could not find a read-id/contig-id column in AMR file: {amr_file}\n"
                f"Columns: {list(df.columns)}"
            )
        id_col = candidates[0]

    gene_candidates = [c for c in df.columns if ("element" in c.lower()) or ("gene" in c.lower())]
    if not gene_candidates:
        raise ValueError(
            f"Could not find a gene/element column in AMR file: {amr_file}\n"
            f"Columns: {list(df.columns)}"
        )

    if "Element symbol" in df.columns:
        gene_col = "Element symbol"
    else:
        gene_col = gene_candidates[0]

    print(f"[INFO] AMR ID column: '{id_col}'")
    print(f"[INFO] AMR gene column: '{gene_col}'")

    out: dict[str, set[str]] = {}
    for rid, gene in zip(df[id_col].astype(str), df[gene_col].astype(str)):
        rid = normalize_read_id(rid)
        gene = gene.strip()
        if rid and gene and rid.lower() != "nan" and gene.lower() != "nan":
            out.setdefault(rid, set()).add(gene)

    return out


def _carbapenemase_label(symbol, name, carb_re):
    symbol = (symbol or "").strip()
    name = (name or "").strip()

    ms = carb_re.search(symbol)
    if ms:
        return symbol if re.match(r"(?i)^bla", symbol) else "bla" + ms.group(0)

    mn = carb_re.search(name)
    if mn:
        return "bla" + mn.group(0)
    return None


def read_carbapenemase_map(amr_file: Path, carbapenemase_regex: str) -> dict[str, set[str]]:
    df = read_table_skip_hash(amr_file)
    if df.empty:
        print(f"[WARN] AMR file is empty: {amr_file}")
        return {}

    lower_to_orig = {c.lower(): c for c in df.columns}


    if "read_id" in df.columns:
        id_col = "read_id"
    elif "contig_id" in df.columns:
        id_col = "contig_id"
    elif "contig id" in lower_to_orig:
        id_col = lower_to_orig["contig id"]
    else:
        cand = ([c for c in df.columns if "read" in c.lower()]
                + [c for c in df.columns if "contig" in c.lower()])
        if not cand:
            raise ValueError(f"No read/contig id column in AMR file: {amr_file}")
        id_col = cand[0]


    if "Element symbol" in df.columns:
        symbol_col = "Element symbol"
    elif "Gene symbol" in df.columns:
        symbol_col = "Gene symbol"
    else:
        gc = [c for c in df.columns
              if ("element" in c.lower() and "symbol" in c.lower())
              or ("gene" in c.lower() and "symbol" in c.lower())
              or c.lower() in ("element", "gene")]
        symbol_col = gc[0] if gc else None


    if "Element name" in df.columns:
        name_col = "Element name"
    elif "Sequence name" in df.columns:
        name_col = "Sequence name"
    else:
        nc = [c for c in df.columns if "name" in c.lower() and "id" not in c.lower()]
        name_col = nc[0] if nc else None

    if symbol_col is None and name_col is None:
        raise ValueError(
            f"No Element symbol / Element name column in AMR file: {amr_file}\n"
            f"Columns: {list(df.columns)}"
        )

    print(f"[INFO] carbapenemase: symbol column '{symbol_col}', name column '{name_col}'")

    carb_re = re.compile(carbapenemase_regex, re.IGNORECASE)
    ids = df[id_col].astype(str)
    syms = df[symbol_col].astype(str) if symbol_col else pd.Series([""] * len(df))
    names = df[name_col].astype(str) if name_col else pd.Series([""] * len(df))

    out: dict[str, set[str]] = {}
    for rid, symbol, name in zip(ids, syms, names):
        rid = normalize_read_id(rid)
        if not rid or rid.lower() == "nan":
            continue
        label = _carbapenemase_label(symbol, name, carb_re)
        if label:
            out.setdefault(rid, set()).add(label)
    return out


def extract_read_types(read_report: Path) -> pd.DataFrame:
    rep = pd.read_csv(read_report, sep="\t", dtype=str)

    if "molecule_type" not in rep.columns:
        raise ValueError(
            f"'molecule_type' column not found in MOBsuite report: {read_report}\n"
            f"Columns: {list(rep.columns)}"
        )

    rep = rep[rep["molecule_type"].str.lower().isin(["chromosome", "plasmid"])]

    if "contig_id" not in rep.columns:
        raise ValueError(
            f"'contig_id' column not found in MOBsuite report: {read_report}\n"
            f"Columns: {list(rep.columns)}"
        )

    rep = rep.copy()
    rep["contig_id"] = rep["contig_id"].astype(str).apply(normalize_read_id)
    rep["molecule_type"] = rep["molecule_type"].str.lower().str.strip()


    if "filtering_reason" in rep.columns:
        fr = rep["filtering_reason"].fillna("").str.strip()
        n_flagged = int((~fr.isin(["-", ""])).sum())
        if n_flagged:
            print(f"[INFO] MOBsuite: kept {n_flagged} contigs that carry a filtering_reason "
                  f"(molecule_type is the sole inclusion criterion)")

    result = rep[["contig_id", "molecule_type"]].drop_duplicates().rename(
        columns={"contig_id": "read_id"}
    )

    n_plasmid = (result["molecule_type"] == "plasmid").sum()
    n_chrom = (result["molecule_type"] == "chromosome").sum()
    print(f"[INFO] MOBsuite: {n_chrom} chromosomes, {n_plasmid} plasmids "
          f"(molecule_type only; filtering_reason ignored)")

    return result


def load_read_scores_long(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=_detect_sep(path), dtype=str)

    required = ["read_id", "motif_mod_position", "mean_prob"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {missing}. "
            f"Columns: {list(df.columns)}"
        )

    df["read_id"] = df["read_id"].astype(str).apply(normalize_read_id)
    df["motif"] = df["motif_mod_position"].astype(str).str.strip()
    df["mean_prob"] = pd.to_numeric(df["mean_prob"], errors="coerce")

    df = df.dropna(subset=["mean_prob"])
    df = df.groupby(["read_id", "motif"], as_index=False)["mean_prob"].median()
    return df


def build_read_arrays(df_long: pd.DataFrame):
    read2motifs: dict[str, np.ndarray] = {}
    read2mean: dict[str, np.ndarray] = {}

    for rid, sub in df_long.groupby("read_id", sort=False):
        m = sub["motif"].to_numpy(dtype=object)
        meanv = sub["mean_prob"].to_numpy(dtype=float)

        order = np.argsort(m)
        read2motifs[rid] = m[order]
        read2mean[rid] = meanv[order]

    return read2motifs, read2mean


def encode_motifs(read2motifs_str, read2mean):
    print("      -> Building global motif vocabulary...")
    all_motifs = sorted(set(m for arr in read2motifs_str.values() for m in arr))
    motif_to_int = {m: np.int32(i) for i, m in enumerate(all_motifs)}
    print(f"      -> {len(all_motifs)} unique motifs in vocabulary")

    read2motifs_int = {}
    for rid, arr in read2motifs_str.items():
        read2motifs_int[rid] = np.array([motif_to_int[m] for m in arr], dtype=np.int32)

    return read2motifs_int, read2mean


@njit
def _rmsd_sorted_int(ma, mb, va, vb):
    i = 0
    j = 0
    sum_sq = 0.0
    nm = 0

    while i < len(ma) and j < len(mb):
        if ma[i] == mb[j]:
            diff = va[i] - vb[j]
            sum_sq += diff * diff
            nm += 1
            i += 1
            j += 1
        elif ma[i] < mb[j]:
            i += 1
        else:
            j += 1

    if nm == 0:
        return np.nan, 0

    return (sum_sq / nm) ** 0.5, nm


def _warmup_numba():
    dummy = np.array([0, 1, 2], dtype=np.int32)
    dummy_f = np.array([0.5, 0.6, 0.7], dtype=np.float64)
    _rmsd_sorted_int(dummy, dummy, dummy_f, dummy_f)


def is_valid_candidate(c, kraken_out, kraken_rank):
    raw = kraken_out.get(c, {}).get("raw", "")
    if not raw:
        return False
    if clean_taxon(raw) == "Unclassified":
        return False
    taxid = kraken_out.get(c, {}).get("taxid", "")
    rank_code = kraken_rank.get(taxid, "")
    if rank_code in ("R", "R1"):
        return False
    return True


def get_motif_filtered_candidates(query_id, chroms_vec, read2motifs, chrom_motif_sets):
    qm = read2motifs.get(query_id)
    if qm is None:
        return chroms_vec
    qm_set = set(qm)
    return [c for c in chroms_vec if qm_set & chrom_motif_sets[c]]


def score_candidates(query_id, candidate_ids, read2motifs, read2mean,
                     kraken_out, min_overlap):
    ma = read2motifs.get(query_id)
    if ma is None:
        return pd.DataFrame(columns=[
            "candidate", "rmsd_mean", "shared_motifs",
            "rmss_mean", "mean_final_score", "species",
        ])

    va_mean = read2mean[query_id]
    rows = []

    for c in candidate_ids:
        if c == query_id:
            continue
        mb = read2motifs.get(c)
        if mb is None:
            continue

        rmsd_mean, nm = _rmsd_sorted_int(ma, mb, va_mean, read2mean[c])

        if nm < min_overlap or not np.isfinite(rmsd_mean):
            continue

        rmss_mean = max(0.0, 1.0 - rmsd_mean)
        raw = kraken_out.get(c, {}).get("raw", "Unclassified")

        rows.append((
            c, rmsd_mean, nm,
            rmss_mean, rmss_mean * nm,
            canonical_species_label(raw),
        ))

    return pd.DataFrame(rows, columns=[
        "candidate", "rmsd_mean", "shared_motifs",
        "rmss_mean", "mean_final_score", "species",
    ])


def tie_species_summary_one_line(df_at_max):
    counts = df_at_max["species"].astype(str).value_counts()
    return "; ".join([f"{sp} ({int(n)})" for sp, n in counts.items()])


def compute_top_summary(cand_df, score_col, atol_tie):
    if cand_df.empty:
        return "NA"
    max_score = float(cand_df[score_col].max())
    top_df = cand_df[cand_df[score_col] >= (max_score - atol_tie)].copy()
    return tie_species_summary_one_line(top_df)


def pick_top_species(cand_df, score_col, atol_tie):
    if cand_df.empty:
        return "no host candidate"
    max_score = float(cand_df[score_col].max())
    top_df = cand_df[cand_df[score_col] >= (max_score - atol_tie)]
    counts = Counter(top_df["species"].astype(str).tolist())
    return counts.most_common(1)[0][0] if counts else "Unclassified"


def rmsd_nm_shared_read(a_id, b_id, read2motifs, read2mean):
    ma = read2motifs.get(a_id)
    mb = read2motifs.get(b_id)
    if ma is None or mb is None:
        return np.nan, 0
    va = read2mean[a_id]
    vb = read2mean[b_id]
    rmsd, nm = _rmsd_sorted_int(ma, mb, va, vb)
    return rmsd, int(nm)


def _taxon_label_for_read(read_id, kraken_out):
    kinfo = kraken_out.get(read_id, {})
    full = full_taxon_label(kinfo.get("raw", ""))
    if not full or full == "Unclassified":
        return "Unclassified"
    parts = full.split()
    if len(parts) >= 2 and parts[1][:1].islower():
        return f"{parts[0]} {parts[1]}"
    return full


def _project_supplementary(d_to_backbone, backbone_coords, backbone_D, eigvals):
    D2 = backbone_D ** 2
    row_means = D2.mean(axis=1)
    d2 = np.asarray(d_to_backbone, dtype=float) ** 2
    safe_lambda = np.where(eigvals > 0, eigvals, np.nan)
    U = backbone_coords / np.sqrt(safe_lambda)
    contrib = (row_means - d2)[:, None] * U
    f = np.nansum(contrib, axis=0) / (2.0 * np.sqrt(safe_lambda))
    return f


def _pcoa_topN_for_query(
    q_id, mol, all_chroms,
    read2motifs, read2mean,
    fs_to_all, bb_taxa_all,
    top_n, out_png, out_pdf,
    legend_loc="lower right",
    legend_fontsize=12.0,
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from skbio import DistanceMatrix
        from skbio.stats.ordination import pcoa as skbio_pcoa
    except ImportError:
        return False

    n_all = len(all_chroms)
    if n_all < 3:
        return False


    fs_arr = np.asarray(fs_to_all)
    order_by_score = np.argsort(-fs_arr, kind="stable")
    sel = [int(i) for i in order_by_score[:top_n] if fs_arr[i] > 0]
    if len(sel) < 3:
        return False
    sub_chroms = [all_chroms[i] for i in sel]
    sub_taxa = [bb_taxa_all[i] for i in sel]
    sub_fs = fs_arr[np.asarray(sel, dtype=int)]


    n_sub = len(sub_chroms)
    FS_sub = np.zeros((n_sub, n_sub), dtype=float)
    for i in range(n_sub):
        for j in range(i + 1, n_sub):
            rmsd, nm = rmsd_nm_shared_read(sub_chroms[i], sub_chroms[j],
                                           read2motifs, read2mean)
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
    coords_full = ord_res.samples.to_numpy()
    eigvals = np.asarray(ord_res.eigvals)
    ve = np.asarray(ord_res.proportion_explained) * 100
    keep_dims = min(ndim, coords_full.shape[1])
    coords_full = coords_full[:, :keep_dims]
    eigvals = eigvals[:keep_dims]
    ve = ve[:keep_dims]
    coords2 = coords_full[:, :2]


    d_to_sub = []
    for c in sub_chroms:
        rmsd, nm = rmsd_nm_shared_read(q_id, c, read2motifs, read2mean)
        if nm < 1 or not np.isfinite(rmsd):
            d_to_sub.append(1.0)
        else:
            s = max(0.0, 1.0 - rmsd) * nm
            d_to_sub.append(float(np.clip(1.0 - s / denom, 0.0, 1.0)))
    f = _project_supplementary(np.array(d_to_sub), coords_full, D_sub, eigvals)
    px, py = float(f[0]), float(f[1])


    N_DISTINCT_TAXA = 5
    top_taxa_set = []
    for k in np.argsort(-sub_fs, kind="stable"):
        t = sub_taxa[int(k)]
        if t in ("Unclassified", "other"):
            continue
        if t not in top_taxa_set:
            top_taxa_set.append(t)
        if len(top_taxa_set) >= N_DISTINCT_TAXA:
            break

    q_color = {t: color_for_species(t) for t in top_taxa_set}
    GRAY = "#cfcfcf"
    q_label = [sub_taxa[k] if sub_taxa[k] in q_color else "other"
               for k in range(n_sub)]
    q_color["other"] = GRAY


    ring_set = {int(np.argmax(sub_fs))} if len(sub_fs) else set()
    ring_idx = sorted(ring_set)


    LABEL_FS = 22
    TICK_FS = 16
    LEGEND_FS = legend_fontsize

    fig, ax = plt.subplots(figsize=(16, 11))
    for k in range(n_sub):
        if k in ring_set:
            continue
        ax.scatter(coords2[k, 0], coords2[k, 1],
                   c=[q_color[q_label[k]]],
                   s=90, edgecolor="white", linewidth=0.5,
                   alpha=1.0, zorder=2)

    for k in ring_idx:
        ax.scatter(coords2[k, 0], coords2[k, 1],
                   c=[q_color[q_label[k]]],
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


    legend_taxa = [t for t in top_taxa_set if t in q_color]
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor=q_color[t], markeredgecolor="none",
                   markersize=11, alpha=1.0, label=t)
        for t in legend_taxa
    ]
    n_taxa_entries = len(handles)
    if "other" in q_label:
        handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor=GRAY, markeredgecolor="none",
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
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def run_pcoa(
    read_types, amr_map, kraken_out, kraken_rank,
    read2motifs, read2mean,
    outdir, min_overlap_motifs,
    only_reads=None,
    topn_pcoa=100,
    legend_loc="lower right",
    legend_fontsize=12.0,
):
    if not only_reads:
        print("      -> no --read-ids: skipping PCoA stage entirely "
              "(per-read figures are only produced for requested reads).")
        return

    try:
        from skbio import DistanceMatrix
        from skbio.stats.ordination import pcoa
    except ImportError:
        print("[PCoA] scikit-bio not installed; skipping PCoA stage. "
              "Install with: pip install scikit-bio")
        return


    chroms = [
        r for r in read_types.query("molecule_type == 'chromosome'")["read_id"].astype(str)
        if r in read2motifs and is_valid_candidate(r, kraken_out, kraken_rank)
    ]
    if len(chroms) < 3:
        print(f"[PCoA] Only {len(chroms)} backbone chromosomes (<3); skipping PCoA.")
        return
    print(f"      -> PCoA backbone: {len(chroms)} classified chromosome reads")

    mol_of = dict(zip(read_types["read_id"].astype(str),
                      read_types["molecule_type"].astype(str)))


    queries = []
    missing = []
    for rid in sorted(only_reads):
        if rid not in amr_map:
            missing.append((rid, "no AMR gene"))
            continue
        if rid not in read2motifs:
            missing.append((rid, "no methylation vector"))
            continue
        if rid not in mol_of:
            missing.append((rid, "not chromosome/plasmid in MOBsuite"))
            continue
        queries.append((rid, mol_of[rid]))

    print(f"      -> PCoA queries: {len(queries)} requested read(s) "
          f"({sum(1 for _, m in queries if m == 'plasmid')} plasmid, "
          f"{sum(1 for _, m in queries if m == 'chromosome')} chromosome)")
    if missing:
        print(f"      -> NOTE: {len(missing)} requested read(s) skipped for PCoA:")
        for rid, why in missing:
            print(f"           {rid}: {why}")
    if not queries:
        print("      -> nothing to plot.")
        return

    topN_dir = outdir / f"pcoa_top{topn_pcoa}" / "per_read"
    (topN_dir / "plasmids").mkdir(parents=True, exist_ok=True)
    (topN_dir / "chromosomes").mkdir(parents=True, exist_ok=True)

    taxa = [_taxon_label_for_read(c, kraken_out) for c in chroms]

    def _safe_name(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))

    n_written = 0
    for q, mol in queries:
        sub = "plasmids" if mol == "plasmid" else "chromosomes"
        stem = f"{_safe_name(q)}.{mol}.pcoa_top{topn_pcoa}"

        fs_to_all_q = []
        for c in chroms:
            if c == q:
                fs_to_all_q.append(0.0)
                continue
            rmsd_q, nm_q = rmsd_nm_shared_read(q, c, read2motifs, read2mean)
            if nm_q < min_overlap_motifs or not np.isfinite(rmsd_q):
                fs_to_all_q.append(0.0)
            else:
                fs_to_all_q.append(max(0.0, 1.0 - rmsd_q) * nm_q)
        ok = _pcoa_topN_for_query(
            q_id=q, mol=mol,
            all_chroms=chroms,
            read2motifs=read2motifs,
            read2mean=read2mean,
            fs_to_all=np.asarray(fs_to_all_q),
            bb_taxa_all=taxa,
            top_n=topn_pcoa,
            out_png=topN_dir / sub / f"{stem}.png",
            out_pdf=topN_dir / sub / f"{stem}.pdf",
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
        )
        if ok:
            n_written += 1
        else:
            print(f"      -> {q}: <3 scoring candidates; no figure written.")

    print(f"      -> wrote {n_written} per-read top-{topn_pcoa} figure(s) under "
          f"{topN_dir.relative_to(outdir)}/{{plasmids,chromosomes}}/")


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
    "Escherichia albertii":                "#742cd8",

    "Citrobacter freundii":                "#3cb44b",
    "Citrobacter arsenatis":               "#b6ecb8",
    "Citrobacter braakii":                 "#D0D38E",
    "Citrobacter europaeus":               "#deed56",
    "Citrobacter koseri":                  "#005212",
    "Citrobacter portucalensis":           "#047904",
    "Citrobacter youngae":                 "#29a395",
    "Citrobacter":                         "#0e684d",
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
    "Caudoviricetes sp.":                  "#7acdc1",
}

PCOA_GRAY = "#cfcfcf"


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


_SPECIES_COLOR_REGISTRY: dict = dict(FIXED_SPECIES_COLORS)


def load_species_color_registry(path) -> None:
    import json
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
    import json
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


_NOBS_CANDIDATE_COLS = (
    "motif_count", "motif_counts",
    "n_motif_obs", "N_motif_obs", "n_motif", "N_motif",
    "n_mod", "N_mod", "motif_obs", "n_obs", "N_obs", "count",
    "N_motif_observation", "n_motif_observation",
)


def build_nobs_lookup(meth_path) -> dict:
    try:
        df = pd.read_csv(meth_path, sep=_detect_sep(meth_path), dtype=str)
    except Exception as e:
        print(f"[nobs] could not read {meth_path} for counts: {e}")
        return {}
    df.columns = [c.strip() for c in df.columns]
    if "read_id" not in df.columns or "motif_mod_position" not in df.columns:
        print("[nobs] missing read_id/motif_mod_position; cannot build count lookup.")
        return {}
    ncol = next((c for c in _NOBS_CANDIDATE_COLS if c in df.columns), None)
    if ncol is None:
        low = {c.lower(): c for c in df.columns}
        ncol = next((low[k] for k in low
                     if ("obs" in k and ("motif" in k or k[:1] == "n"))
                     or ("motif" in k and "count" in k)), None)
    if ncol is None:
        print("[nobs] no observation-count column found; cells will show 'n=?'.")
        return {}
    rid = df["read_id"].astype(str).apply(normalize_read_id)
    key = df["motif_mod_position"].astype(str).str.strip()
    n = pd.to_numeric(df[ncol], errors="coerce")
    g = pd.DataFrame({"read_id": rid, "key": key, "n": n}).dropna(subset=["n"])
    g = g.groupby(["read_id", "key"], as_index=False)["n"].sum()
    print(f"[nobs] using '{ncol}' as the observation count.")
    return {(r, k): float(v) for r, k, v in zip(g["read_id"], g["key"], g["n"])}


def _italic_taxon(s: str) -> str:
    esc = str(s).replace("\\", " ").replace("$", "")
    for ch in ("_", "^", "{", "}"):
        esc = esc.replace(ch, " ")
    esc = esc.replace(" ", r"\ ")
    return "$\\it{" + esc + "}$"


def _rmsd_nm_mean_read(a_id, b_id, read2motifs, read2mean):
    ma = read2motifs.get(a_id)
    mb = read2motifs.get(b_id)
    if ma is None or mb is None:
        return np.nan, 0
    va = read2mean[a_id]
    vb = read2mean[b_id]
    rmsd, nm = _rmsd_sorted_int(ma, mb, va, vb)
    return rmsd, int(nm)


def _score_candidates_for_heatmap(
    query_id, candidate_ids,
    read2motifs, read2mean,
    kraken_out, min_overlap,
):
    rows = []
    for c in candidate_ids:
        if c == query_id:
            continue
        rmsd, nm = _rmsd_nm_mean_read(query_id, c, read2motifs, read2mean)
        if nm < min_overlap or not np.isfinite(rmsd):
            continue
        rmss = max(0.0, 1.0 - rmsd)
        final_score = rmss * nm
        kinfo = kraken_out.get(c, {})
        raw = kinfo.get("raw", "Unclassified")
        species_can = canonical_species_label(raw)
        full_label = full_taxon_label(raw)
        taxid = kinfo.get("taxid", "")
        rows.append((c, rmsd, nm, rmss, final_score, species_can, full_label, taxid))

    df = pd.DataFrame(
        rows,
        columns=["candidate", "rmsd", "shared_motifs", "rmss", "final_score",
                 "species", "full_label", "taxid"],
    )
    if df.empty:
        return df
    df = df.sort_values(
        ["final_score", "shared_motifs", "rmsd"], ascending=[False, False, True]
    ).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


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
    id2motifs,
    id2vals,
    out_png,
    out_pdf,
    nobs: dict | None = None,
    all_motifs: bool = False,
    max_motifs: int = 0,
    order: str = "methylation",
    orient: str = "motifs-rows",
    title: str | None = None,
    fontsize: float = 11.0,
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    nobs = nobs or {}

    def _meth_dict(cid):
        ms = id2motifs.get(cid)
        vs = id2vals.get(cid)
        if ms is None or vs is None:
            return {}
        return dict(zip(ms.tolist(), vs.tolist()))

    q_meth = _meth_dict(query_id)
    chrom_meth = {r["contig"]: _meth_dict(r["contig"]) for r in picked}

    chrom_union = set().union(*[set(d) for d in chrom_meth.values()]) if chrom_meth else set()
    keys = (set(q_meth) | chrom_union) if all_motifs else set(q_meth)
    if not keys:
        print(f"[heatmap] {query_id}: no motifs to plot "
              f"({'union empty' if all_motifs else 'query has no motifs'}); skipping.")
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

    col_labels = [f"({q_tag})"]
    for idx, r in enumerate(picked):
        is_top = (idx == 0)
        sp = r["species"]
        if not sp or sp == "Unclassified":
            col_labels.append("Unclassified")
        elif is_top:
            col_labels.append(sp)
        else:
            col_labels.append(_italic_taxon(sp))

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
        _tls = ax.get_yticklabels()
        if ncols >= 2:
            _tls[1].set_fontweight("bold"); _tls[1].set_fontstyle("italic")
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
        _tls = ax.get_xticklabels()
        if ncols >= 2:
            _tls[1].set_fontweight("bold"); _tls[1].set_fontstyle("italic")
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
    cbar.set_label("mean methylation value", fontsize=fs)
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
          f"({'union' if all_motifs else 'all-query-motifs'}); "
          f"hosts: {species_list} -> {out_png}")


def run_top_taxa_heatmaps(
    read_types, amr_map, kraken_out, kraken_rank,
    read2motifs, read2motifs_str, read2mean,
    outdir, score_path,
    min_overlap_motifs: int = 1,
    n_taxa: int = 5,
    all_motifs: bool = False,
    max_motifs: int = 0,
    order: str = "methylation",
    orient: str = "motifs-rows",
    fontsize: float = 11.0,
    only_reads: set | None = None,
    color_map_path=None,
):
    if not only_reads:
        print("      -> no --read-ids: skipping top-taxa heatmaps "
              "(figures are only written for requested reads).")
        return

    base = outdir / "top_taxa_heatmap" / "per_read"

    def _safe_name(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))

    chrom_ids = [
        c for c in read_types.query("molecule_type == 'chromosome'")["read_id"].astype(str)
        if c in read2motifs and is_valid_candidate(c, kraken_out, kraken_rank)
    ]
    if not chrom_ids:
        print("[heatmap] no chromosome candidates with motif vectors; skipping.")
        return

    mol_of = dict(zip(read_types["read_id"].astype(str),
                      read_types["molecule_type"].astype(str)))

    queries = [r for r in read_types["read_id"].astype(str)
               if r in read2motifs and amr_map.get(r)]
    queries = [q for q in queries if q in only_reads]
    missing = sorted(set(only_reads) - set(queries))
    if missing:
        print(f"[heatmap] NOTE: {len(missing)} requested read(s) are not AMR "
              f"queries with motifs and will be skipped: {missing}")
    if not queries:
        print("[heatmap] nothing to do after --read-ids filter.")
        return

    load_species_color_registry(color_map_path)
    nobs = build_nobs_lookup(score_path)

    for q in queries:
        mol = mol_of.get(q, "unknown")
        sub = "plasmids" if mol == "plasmid" else "chromosomes"
        cands = [c for c in chrom_ids if c != q]
        cand_df = _score_candidates_for_heatmap(
            query_id=q, candidate_ids=cands,
            read2motifs=read2motifs, read2mean=read2mean,
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
        stem = f"{_safe_name(q)}.{mol}.top{n_taxa}_taxa_heatmap"
        plot_top_taxa_heatmap(
            query_id=q, query_mol=mol, picked=picked,
            id2motifs=read2motifs_str, id2vals=read2mean,
            out_png=base / sub / f"{stem}.png",
            out_pdf=base / sub / f"{stem}.pdf",
            nobs=nobs, all_motifs=all_motifs, max_motifs=max_motifs,
            order=order, orient=orient, fontsize=fontsize,
        )

    save_species_color_registry(color_map_path)


def _score_plasmid_worker(p):
    kout = _kraken_out.get(p, {})
    taxid = kout.get("taxid", "")
    rank_code = _kraken_rank.get(taxid, "")
    raw = kout.get("raw", "Unclassified")

    row = {
        "read_id": p,
        "molecule_type": "plasmid",
        "amr_genes": "; ".join(sorted(_amr_map.get(p, set()))),
        "has_nanomotif_vector": (p in _read2motifs),
        "kraken_raw": raw,
        "kraken_rank_code": rank_code if rank_code else "NA",
        "mean_top_score_species_counts": "NA",
    }

    if p not in _read2motifs:
        return row

    filtered = get_motif_filtered_candidates(p, _chroms_vec, _read2motifs, _chrom_motif_sets)
    cand_df = score_candidates(
        p, filtered, _read2motifs, _read2mean, _kraken_out, _min_overlap
    )

    if not cand_df.empty:
        row["mean_top_score_species_counts"] = compute_top_summary(
            cand_df, "mean_final_score", _atol_tie
        )

    return row


def run_assignment(
    read_types, amr_map, kraken_out, kraken_rank,
    read2motifs, read2mean,
    outdir, min_overlap_motifs, atol_tie, n_workers,
    only_reads=None,
):
    outdir.mkdir(parents=True, exist_ok=True)

    valid_mob_reads = set(read_types["read_id"].tolist())
    plasmids_all    = read_types.query("molecule_type == 'plasmid'")["read_id"].astype(str).tolist()
    chromosomes_all = read_types.query("molecule_type == 'chromosome'")["read_id"].astype(str).tolist()


    if only_reads:
        def _keep(r):
            return r in only_reads
        not_amr   = sorted(r for r in only_reads if r not in amr_map)
        not_inmob = sorted(r for r in only_reads if r not in valid_mob_reads)
        print(f"      -> --read-ids: restricting associations to {len(only_reads)} requested read(s)")
        if not_amr:
            print(f"      -> NOTE: {len(not_amr)} requested read(s) have no AMR gene "
                  f"and will not appear in the summaries: {not_amr}")
        if not_inmob:
            print(f"      -> NOTE: {len(not_inmob)} requested read(s) absent from MOBsuite "
                  f"output and will be skipped: {not_inmob}")
    else:
        def _keep(r):
            return True

    chroms_vec = [
        r for r in chromosomes_all
        if r in read2motifs and is_valid_candidate(r, kraken_out, kraken_rank)
    ]

    print(f"      -> {len(chromosomes_all)} chromosomes in MOBsuite output")
    print(f"      -> {len(plasmids_all)} plasmids in MOBsuite output")
    print(f"      -> {len(chroms_vec)} classified chromosomes in candidate pool (excl. R/R1)")

    if len(chromosomes_all) > 0:
        frac = len(chroms_vec) / len(chromosomes_all)
        if frac < 0.1:
            print(f"[WARN] Only {frac:.1%} of chromosomes in candidate pool")

    amr_not_in_mob = [r for r in amr_map if r not in valid_mob_reads]
    if amr_not_in_mob:
        print(f"[WARN] {len(amr_not_in_mob)} AMR+ reads NOT in MOBsuite output — will be skipped")

    print("      -> Precomputing chromosome motif sets...")
    chrom_motif_sets = {c: set(read2motifs[c]) for c in chroms_vec}
    print(f"      -> Done ({len(chrom_motif_sets)} sets)")

    pool_kwargs = dict(
        processes=n_workers,
        initializer=_worker_init,
        initargs=(
            read2motifs, read2mean,
            kraken_out, kraken_rank, amr_map,
            chroms_vec, chrom_motif_sets,
            min_overlap_motifs, atol_tie,
        ),
    )


    amr_plasmids = [p for p in plasmids_all if p in amr_map and _keep(p)]
    n_plasmid_no_amr = len([p for p in plasmids_all if _keep(p)]) - len(amr_plasmids)
    n_plasmid_written = 0

    plasmid_out_path = outdir / "AMR_plasmid_read_assignment_summary.tsv"
    print(f"\n[PLASMIDS] {len(amr_plasmids)} AMR+ plasmids → {plasmid_out_path}")

    with plasmid_out_path.open("w", newline="") as pf:
        writer = csv.DictWriter(pf, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        pf.flush()

        with Pool(**pool_kwargs) as pool:
            for row in pool.imap_unordered(_score_plasmid_worker, amr_plasmids, chunksize=20):
                writer.writerow(row)
                pf.flush()
                n_plasmid_written += 1
                if n_plasmid_written % 100 == 0:
                    print(f"      -> {n_plasmid_written}/{len(amr_plasmids)} plasmids scored...")

    print(f"      -> {n_plasmid_no_amr} plasmids skipped (no AMR genes)")
    print(f"      -> {n_plasmid_written} AMR+ plasmids written")

    print("\n[DONE] Outputs written to:", str(outdir))
    print(f"  - AMR_plasmid_read_assignment_summary.tsv                  ({n_plasmid_written} rows)")


def run_carbapenemase_summary(
    read_types, carbapenemase_map, kraken_out, kraken_rank,
    read2motifs, read2mean,
    outdir, min_overlap_motifs, atol_tie,
    dataset,
):
    chromosomes_all = read_types.query("molecule_type == 'chromosome'")["read_id"].astype(str).tolist()
    chroms_vec = [
        r for r in chromosomes_all
        if r in read2motifs and is_valid_candidate(r, kraken_out, kraken_rank)
    ]
    chrom_motif_sets = {c: set(read2motifs[c]) for c in chroms_vec}

    # Only PLASMID-borne carbapenemase reads are summarised (plasmid -> host);
    # carbapenemase reads on chromosomes would be chromosome-vs-chromosome and
    # are excluded.
    plasmid_reads = set(read_types.query("molecule_type == 'plasmid'")["read_id"].astype(str))

    gene_hits: dict[str, list[str]] = {}

    for rid in sorted(carbapenemase_map):
        if rid not in plasmid_reads:
            continue
        carb_genes = sorted(carbapenemase_map[rid])
        if not carb_genes:
            continue

        if rid not in read2motifs:
            top_sp = "no methylation vector"
        else:
            cand_ids = [c for c in chroms_vec if c != rid]
            filtered = get_motif_filtered_candidates(rid, cand_ids, read2motifs, chrom_motif_sets)
            cand_df = score_candidates(rid, filtered, read2motifs, read2mean,
                                       kraken_out, min_overlap_motifs)
            top_sp = pick_top_species(cand_df, "mean_final_score", atol_tie)

        for g in carb_genes:
            gene_hits.setdefault(g, []).append(top_sp)

    rows = []
    for gene in sorted(gene_hits):
        species = gene_hits[gene]
        counts = Counter(species)

        parts = [f"{sp} ({n})" for sp, n in
                 sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
        rows.append({
            "dataset": dataset,
            "carbapenemase": gene,
            "n_reads": len(species),
            "top_read_similarity": ", ".join(parts),
        })

    out_df = pd.DataFrame(
        rows, columns=["dataset", "carbapenemase", "n_reads", "top_read_similarity"]
    )
    out_path = outdir / "carbapenemase_read_host_summary.tsv"


    cols = ["dataset", "carbapenemase", "n_reads", "top_read_similarity"]
    with out_path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    total_hits = sum(len(v) for v in gene_hits.values())
    print(f"\n[CARBAPENEMASE SUMMARY] {len(rows)} carbapenemase gene(s) over "
          f"{total_hits} read-hit(s) → {out_path}")
    return out_df


def main():
    ap = argparse.ArgumentParser(
        description=(
            "AMR read association summary using Nanomotif-style score = (1 - RMSD) * NM\n"
            "on the MEAN methylation probability. Uses numba JIT + multiprocessing\n"
            "Pool for fast scoring.\n\n"
            "Filtering rules:\n"
            "  - A read is included if MOBsuite labelled it chromosome or plasmid\n"
            "    (molecule_type is the SOLE inclusion gate; filtering_reason is ignored,\n"
            "    so repetitive/IS-element reads such as ISKpn7-blaKPC are kept)\n"
            "  - Reads absent from MOBsuite contig_report.txt are excluded entirely\n"
            "  - Plasmids without AMR genes are skipped\n"
            "  - Chromosomes without AMR genes are used only as candidates\n"
            "  - Chromosomes at Kraken rank R/R1 are excluded from candidate pool\n"
            "  - Per-query motif pre-filtering applied before scoring\n"
            "  - Results written to disk immediately after each association"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--read-motif-dir", type=Path, required=True)
    ap.add_argument("--read-types",     type=Path, required=True)
    ap.add_argument("--amr-dir",        type=Path, required=True)
    ap.add_argument("--kraken-dir",     type=Path, required=True)
    ap.add_argument("--outdir",         type=Path, required=True)
    ap.add_argument(
        "--min-overlap-motifs", type=int, default=MIN_OVERLAP_MOTIFS_DEFAULT,
        help=f"Minimum shared motifs required (default: {MIN_OVERLAP_MOTIFS_DEFAULT})"
    )
    ap.add_argument(
        "--atol-tie", type=float, default=ATOL_TIE_DEFAULT,
        help=f"Tolerance for tie-breaking at top score (default: {ATOL_TIE_DEFAULT})"
    )
    ap.add_argument(
        "--workers", type=int, default=8,
        help="Number of parallel worker processes (default: 8)"
    )
    ap.add_argument(
        "--dataset", type=str, default=None,
        help="Dataset / sample label written in the 'dataset' column of the "
             "carbapenemase summary (default: the --outdir folder name)."
    )
    ap.add_argument(
        "--carbapenemase-regex", type=str, default=CARBAPENEMASE_REGEX_DEFAULT,
        help="Case-insensitive regex matching carbapenemase gene symbols for the "
             f"per-gene read summary (default: {CARBAPENEMASE_REGEX_DEFAULT})"
    )


    ap.add_argument(
        "--read-ids", type=str, default=None,
        help="Comma-separated read IDs to focus on (e.g. 'read123,read456'). "
             "When given, associations are computed and written ONLY for these "
             "reads, and a per-read top-N PCoA figure is produced for each. "
             "This is the fast path: the full sample is still used as the "
             "candidate/backbone pool, but only the listed reads are scored."
    )
    ap.add_argument(
        "--read-ids-file", type=Path, default=None,
        help="Path to a file with one read ID per line (blank lines and lines "
             "starting with '#' are ignored). Merged with --read-ids if both given."
    )
    ap.add_argument(
        "--associate-all", action="store_true",
        help="By default, when --read-ids is given the association TSVs are "
             "restricted to those reads (the fast path). With this flag, the "
             "association is computed for ALL AMR-carrying reads regardless of "
             "--read-ids (slower). Either way, PCoA + heatmap figures are still "
             "produced ONLY for the --read-ids reads."
    )
    ap.add_argument(
        "--topn-pcoa", type=int, default=100,
        help="Number of top candidate chromosomes (by methylation similarity) to "
             "include in each per-read PCoA. Output goes to "
             "pcoa_top{N}/per_read/... (default: 100)"
    )
    ap.add_argument(
        "--no-pcoa", action="store_true",
        help="Skip the PCoA figures even when --read-ids is given "
             "(only write the restricted association summaries)."
    )
    ap.add_argument(
        "--legend-loc", type=str, default="lower right",
        choices=["best", "upper right", "upper left", "lower left",
                 "lower right", "right", "center left", "center right",
                 "lower center", "upper center", "center",
                 "outside", "none", "off"],
        help="Where to place the per-read PCoA legend (default: 'lower right'). "
             "'outside' places it to the right of the plot; 'none'/'off' hides it."
    )
    ap.add_argument(
        "--legend-fontsize", type=float, default=12.0,
        help="Font size of the per-read PCoA legend entries (default: 12)."
    )


    ap.add_argument(
        "--top-taxa-heatmap", action="store_true",
        help="Also write a motif x read MEAN-methylation heatmap per requested "
             "read: the query next to one chromosome from each of the top-N "
             "DISTINCT SPECIES-LEVEL taxa (best final_score per species; Kraken2 "
             "rank 'S'), with species-coloured host labels. Genus-only / family / "
             "unclassified candidates are skipped. Only written for reads named "
             "by --read-ids. Output: <outdir>/top_taxa_heatmap/per_read/..."
    )
    ap.add_argument(
        "--heatmap-n-taxa", type=int, default=5,
        help="Number of distinct species-level taxa / chromosome columns in the "
             "top-taxa heatmap (default 5)."
    )
    ap.add_argument(
        "--heatmap-min-overlap", type=int, default=None,
        help="Minimum shared motifs for a chromosome to qualify as a heatmap "
             "candidate (default: same as --min-overlap-motifs)."
    )
    ap.add_argument(
        "--heatmap-all-motifs", action="store_true",
        help="Heatmap rows = motifs in the query OR any chosen chromosome (adds "
             "host-only motifs, which show n/a in the query column). Default: "
             "every motif the query read carries."
    )
    ap.add_argument(
        "--heatmap-max-motifs", type=int, default=0,
        help="Cap the number of heatmap rows (0 = all)."
    )
    ap.add_argument(
        "--heatmap-order", choices=["methylation", "nobs", "name"],
        default="methylation",
        help="Heatmap row ordering, top to bottom (default methylation = by mean "
             "methylation)."
    )
    ap.add_argument(
        "--heatmap-orient", choices=["motifs-rows", "contigs-rows"],
        default="motifs-rows",
        help="Heatmap layout (default motifs-rows: motifs down the side, reads "
             "as columns)."
    )
    ap.add_argument(
        "--heatmap-fontsize", type=float, default=11.0,
        help="Base font size for the top-taxa heatmap (default 11)."
    )
    ap.add_argument(
        "--heatmap-color-map", type=Path, default=None,
        help="OPTIONAL. Species colours are already baked into the script "
             "(FIXED_SPECIES_COLORS), so by default no external file is used. "
             "Pass a JSON path here only to PERSIST colours of NEW (non-pinned) "
             "species across separate runs/samples -- point ALL runs at the same "
             "json (e.g. the one your contig PCoA uses) to keep them identical."
    )

    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)


    only_reads = set()
    if args.read_ids:
        only_reads.update(normalize_read_id(x) for x in args.read_ids.split(",") if x.strip())
    if args.read_ids_file:
        if not args.read_ids_file.is_file():
            raise SystemExit(f"[STOP] --read-ids-file not found: {args.read_ids_file}")
        with args.read_ids_file.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                only_reads.add(normalize_read_id(line))
    only_reads = only_reads or None
    if only_reads:
        if args.associate_all:
            print(f"[CONFIG] --associate-all: association computed for ALL AMR reads; "
                  f"PCoA + heatmap figures restricted to {len(only_reads)} requested read(s).")
        else:
            print(f"[CONFIG] association restricted to {len(only_reads)} requested read(s) "
                  f"(default; use --associate-all for every AMR read). "
                  f"PCoA + heatmap figures for the same read(s).")
    else:
        print("[CONFIG] no --read-ids: association for ALL AMR reads; "
              "no PCoA/heatmap figures will be produced.")

    score            = find_one(args.read_motif_dir, ["*.tsv", "*.csv"], "read motif summary (TSV/CSV)")
    amr_file         = find_one(args.amr_dir, ["*.tsv", "*.txt"], "AMR read output")
    kraken_out_file  = find_one(args.kraken_dir, ["*.out", "*kraken*.out", "*kraken*.txt"], "Kraken OUT")
    kraken_rep_file  = find_one(args.kraken_dir, ["*.report", "*kraken*.report"], "Kraken REPORT")

    print("[1/7] Loading read-level motif table...")
    df_long = load_read_scores_long(score)
    print(f"      -> {len(df_long)} rows loaded")

    print("[2/7] Building read arrays...")
    read2motifs_str, read2mean = build_read_arrays(df_long)
    print(f"      -> {len(read2motifs_str)} reads with motif vectors")

    print("[3/7] Encoding motifs to integers for numba...")
    read2motifs, read2mean = encode_motifs(read2motifs_str, read2mean)

    print("[4/7] Warming up numba JIT in main process...")
    _warmup_numba()
    print("      -> Numba JIT ready")

    print("[5/7] Loading MOBsuite read types...")
    read_types = extract_read_types(args.read_types)
    print(f"      -> {len(read_types)} reads after MOBsuite filtering")

    print("[6/7] Loading AMR map + Kraken2...")
    amr_map     = read_amr_map(amr_file)
    carbapenemase_map = read_carbapenemase_map(amr_file, args.carbapenemase_regex)
    kraken_out  = load_kraken_out(kraken_out_file)
    kraken_rank = load_kraken_report_rank_map(kraken_rep_file)
    print(f"      -> {len(amr_map)} AMR+ reads")
    print(f"      -> {len(carbapenemase_map)} reads carrying a carbapenemase")
    print(f"      -> {len(kraken_out)} reads in Kraken2 OUT")
    print(f"      -> {len(kraken_rank)} taxids in Kraken2 REPORT")

    print(f"[7/7] Running assignment with {args.workers} workers...")
    run_assignment(
        read_types=read_types,
        amr_map=amr_map,
        kraken_out=kraken_out,
        kraken_rank=kraken_rank,
        read2motifs=read2motifs,
        read2mean=read2mean,
        outdir=args.outdir,
        min_overlap_motifs=args.min_overlap_motifs,
        atol_tie=args.atol_tie,
        n_workers=args.workers,


        only_reads=(None if args.associate_all else only_reads),
    )


    print("\n[SUMMARY] Building per-carbapenemase read-host summary...")
    run_carbapenemase_summary(
        read_types=read_types,
        carbapenemase_map=carbapenemase_map,
        kraken_out=kraken_out,
        kraken_rank=kraken_rank,
        read2motifs=read2motifs,
        read2mean=read2mean,
        outdir=args.outdir,
        min_overlap_motifs=args.min_overlap_motifs,
        atol_tie=args.atol_tie,
        dataset=(args.dataset if args.dataset else args.outdir.name),
    )


    if only_reads and not args.no_pcoa:
        print("\n[PCoA] Generating per-read top-N PCoA figures...")
        run_pcoa(
            read_types=read_types,
            amr_map=amr_map,
            kraken_out=kraken_out,
            kraken_rank=kraken_rank,
            read2motifs=read2motifs,
            read2mean=read2mean,
            outdir=args.outdir,
            min_overlap_motifs=args.min_overlap_motifs,
            only_reads=only_reads,
            topn_pcoa=args.topn_pcoa,
            legend_loc=args.legend_loc,
            legend_fontsize=args.legend_fontsize,
        )
    elif args.no_pcoa and only_reads:
        print("\n[PCoA] --no-pcoa set; skipping figures (associations only).")


    if only_reads and args.top_taxa_heatmap:
        print("\n[HEATMAP] Generating top-taxa methylation heatmaps...")
        hm_overlap = (int(args.heatmap_min_overlap)
                      if args.heatmap_min_overlap is not None
                      else int(args.min_overlap_motifs))
        hm_color_map = args.heatmap_color_map
        run_top_taxa_heatmaps(
            read_types=read_types,
            amr_map=amr_map,
            kraken_out=kraken_out,
            kraken_rank=kraken_rank,
            read2motifs=read2motifs,
            read2motifs_str=read2motifs_str,
            read2mean=read2mean,
            outdir=args.outdir,
            score_path=score,
            min_overlap_motifs=hm_overlap,
            n_taxa=int(args.heatmap_n_taxa),
            all_motifs=bool(args.heatmap_all_motifs),
            max_motifs=int(args.heatmap_max_motifs),
            order=args.heatmap_order,
            orient=args.heatmap_orient,
            fontsize=float(args.heatmap_fontsize),
            only_reads=only_reads,
            color_map_path=hm_color_map,
        )
    elif args.top_taxa_heatmap and not only_reads:
        print("\n[HEATMAP] --top-taxa-heatmap set but no --read-ids; "
              "no heatmaps written (figures are per-requested-read only).")


if __name__ == "__main__":
    main()