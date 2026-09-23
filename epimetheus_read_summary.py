#!/usr/bin/env python3
"""Aggregate per-site epimetheus_main.tsv into the per-read motif summary
expected by cupid_read.py (--read-motif-dir input).

Uses polars (multithreaded, streaming). Install with: pip install polars

Output columns: read_id, motif_mod_position, mean_prob, motif_count
  motif_mod_position = motif_modtype_modposition (e.g. GATC_a_1)
  mean_prob          = mean(quality / 255) over occurrences in the read
  motif_count        = number of occurrences (used for 'n=' heatmap labels)
"""

import argparse
from pathlib import Path

import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", type=Path, required=True,
                    help="epimetheus_main.tsv (per-site)")
    ap.add_argument("--out", type=Path, required=True,
                    help="per-read summary TSV for cupid_read.py")
    ap.add_argument("--min-sites", type=int, default=1,
                    help="drop read/motif pairs with fewer occurrences (default 1)")
    args = ap.parse_args()

    lf = (
        pl.scan_csv(
            args.inp, separator="\t",
            schema_overrides={"read_id": pl.Utf8, "motif_seq": pl.Utf8,
                              "mod_type": pl.Utf8, "mod_pos": pl.Utf8,
                              "quality": pl.Float64},
        )
        .select(["read_id", "motif_seq", "mod_type", "mod_pos", "quality"])
        .drop_nulls("quality")
        .with_columns(
            pl.concat_str(
                [pl.col("motif_seq").str.strip_chars(),
                 pl.col("mod_type").str.strip_chars(),
                 pl.col("mod_pos").str.strip_chars()],
                separator="_",
            ).alias("motif_mod_position"),
            (pl.col("quality") / 255.0).alias("prob"),
        )
        .group_by(["read_id", "motif_mod_position"])
        .agg(
            pl.col("prob").mean().alias("mean_prob"),
            pl.len().alias("motif_count"),
        )
        .filter(pl.col("motif_count") >= args.min_sites)
        .select(["read_id", "motif_mod_position", "mean_prob", "motif_count"])
        .sort("read_id")
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    lf.sink_csv(args.out, separator="\t")

    n = pl.scan_csv(args.out, separator="\t").select(pl.len()).collect().item()
    print(f"[OK] {n:,} read/motif rows -> {args.out}")


if __name__ == "__main__":
    main()