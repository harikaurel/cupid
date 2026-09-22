#!/usr/bin/env python3
"""Aggregate per-site epimetheus_main.tsv into the per-read motif summary
expected by cupid_read.py (--read-motif-dir input).

Output columns: read_id, motif_mod_position, mean_prob, motif_count
  mean_prob   = mean(quality / 255) over all occurrences of the motif in the read
  motif_count = number of occurrences (picked up by cupid_read.py for the
                'n=' annotations in the top-taxa heatmaps)
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", type=Path, required=True,
                    help="epimetheus_main.tsv (per-site)")
    ap.add_argument("--out", type=Path, required=True,
                    help="per-read summary TSV for cupid_read.py")
    ap.add_argument("--min-sites", type=int, default=1,
                    help="drop read/motif pairs with fewer occurrences (default 1)")
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t",
                     usecols=["read_id", "motif_seq", "mod_type", "mod_pos", "quality"],
                     dtype={"read_id": str, "motif_seq": str, "mod_type": str})

    df["quality"] = pd.to_numeric(df["quality"], errors="coerce")
    df = df.dropna(subset=["quality"])

    df["motif_mod_position"] = (df["motif_seq"].str.strip() + "_"
                                + df["mod_type"].str.strip() + "-"
                                + df["mod_pos"].astype(str).str.strip())
    df["prob"] = df["quality"] / 255.0

    g = (df.groupby(["read_id", "motif_mod_position"], sort=False)
           .agg(mean_prob=("prob", "mean"), motif_count=("prob", "size"))
           .reset_index())

    if args.min_sites > 1:
        g = g[g["motif_count"] >= args.min_sites]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] {len(g)} read/motif rows from {df['read_id'].nunique()} reads -> {args.out}")


if __name__ == "__main__":
    main()