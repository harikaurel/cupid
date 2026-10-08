#!/usr/bin/env python3
"""
CUPID contig fragmentation, two modes.

Given a contig of interest (the query), this script fragments MobSuite
chromosomes so that each fragment is comparable to the query, then writes
the fragmented assembly and remapped annotation tables.

Modes (--mode):
  motif   (default) Each fragment carries the same total number of motif
          occurrences as the query (n_target = sum of the query's motif
          counts, from n_motif_obs if present, otherwise a sequence scan).
          Only chromosomes sharing >= 1 motif with the query are cut; they
          are cut on all of their own motifs every n_target occurrences.
          Chromosomes whose total counts are below n_target are left whole.
  length  Each fragment has the same length as the query (n_target = query
          length in bp). All chromosomes are cut, no motif logic. A trailing
          piece shorter than --min-frac * query length is merged into the
          previous fragment. Chromosomes shorter than the query are left whole.

Plasmids, unclassified contigs and the query itself are never fragmented.

Outputs in --outdir (prefix = query contig name):
  <prefix>_fragmented.fasta
  <prefix>_fragments.tsv                fragment, parent, start, end, length,
                                        n_motifs, n_target, mode
  <prefix>_query_summary.tsv            (if --motifs given)
  <prefix>_contig_report.frag.txt
  <prefix>_kraken2.frag.output          (if --kraken2 given)
  <prefix>_amrfinder.frag.tsv           (if --amrfinder given)

In length mode, n_motifs is NA and n_target is the query length (bp).
Use a separate --outdir per mode; filenames are the same in both.

Usage:
  python fragmentation.py --mode motif \
      --motifs epimetheus_table.tsv --assembly RS7.fasta \
      --mobsuite contig_report.txt --kraken2 RS7.kraken2.out \
      --amrfinder RS7.amrfinder.tsv --contig ctg2602 --outdir frag_motif/ctg2602

  python fragmentation.py --mode length --min-frac 0.5 \
      --assembly RS7.fasta --mobsuite contig_report.txt \
      --kraken2 RS7.kraken2.out --amrfinder RS7.amrfinder.tsv \
      --contig ctg2602 --outdir frag_length/ctg2602
"""

import argparse
import os
import re
import sys
from collections import defaultdict

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

COMP = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


def revcomp(seq: str) -> str:
    return "".join(COMP[b] for b in reversed(seq))


def motif_regex(motif: str) -> re.Pattern:
    return re.compile("(?=" + "".join(IUPAC[b] for b in motif) + ")")


def motif_positions(seq: str, motif: str) -> list:
    pos = [m.start() for m in motif_regex(motif).finditer(seq)]
    rc = revcomp(motif)
    if rc != motif:
        pos += [m.start() for m in motif_regex(rc).finditer(seq)]
    return pos


def read_fasta(path: str) -> dict:
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_motif_table(path: str, min_obs: int = 1):
    """Returns ({contig: set of motif tuples}, {contig: {motif tuple: n_obs}}).

    Works with nanomotif motifs.tsv and with epimetheus methylation-pattern
    output. If an n_motif_obs column is present (epimetheus), a motif counts
    as present on a contig only when n_motif_obs >= min_obs, and the obs
    dict is populated; otherwise the obs dict is empty.
    """
    per_contig = defaultdict(set)
    per_contig_obs = defaultdict(dict)
    n_dropped = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("contig", "motif", "mod_type", "mod_position"):
            if col not in idx:
                sys.exit(f"missing column '{col}' in {path}")
        obs_i = idx.get("n_motif_obs")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[idx["motif"]].upper(), f[idx["mod_type"]],
                   f[idx["mod_position"]])
            if obs_i is not None:
                try:
                    obs = int(float(f[obs_i]))
                except ValueError:
                    obs = 0
                if obs < min_obs:
                    n_dropped += 1
                    continue
                per_contig_obs[f[idx["contig"]]][key] = obs
            per_contig[f[idx["contig"]]].add(key)
    if obs_i is not None:
        print(f"[INFO] epimetheus-style table: presence = n_motif_obs >= "
              f"{min_obs}, {n_dropped} rows dropped", file=sys.stderr)
    return per_contig, per_contig_obs


def norm_id(cid: str) -> str:
    return cid.split("|")[-1].split()[0]


def read_mobsuite_types(path: str) -> dict:
    types = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("contig_id", "molecule_type"):
            if col not in idx:
                sys.exit(f"missing column '{col}' in {path}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            types[norm_id(f[idx["contig_id"]])] = f[idx["molecule_type"]].lower()
    return types


def fragment_by_motifs(seq, hit_positions, n_target):
    """-> [(start, end, n_motifs)], cut every n_target occurrences."""
    if not hit_positions:
        return [(0, len(seq), 0)]
    hits = sorted(hit_positions)
    frags = []
    start, count = 0, 0
    for i, p in enumerate(hits):
        count += 1
        if count == n_target:
            end = hits[i + 1] if i + 1 < len(hits) else len(seq)
            frags.append((start, end, count))
            start, count = end, 0
    if start < len(seq):
        frags.append((start, len(seq), count))
    return frags


def fragment_by_length(seq_len, win, min_frac):
    """-> [(start, end)], windows of `win` bp; a trailing piece shorter than
    min_frac * win is merged into the previous window."""
    if seq_len <= win:
        return [(0, seq_len)]
    frags = [(s, min(s + win, seq_len)) for s in range(0, seq_len, win)]
    if len(frags) > 1 and (frags[-1][1] - frags[-1][0]) < min_frac * win:
        last = frags.pop()
        prev = frags.pop()
        frags.append((prev[0], last[1]))
    return frags


def remap_table(infile, outfile, frag_map, column=None, col_index=None,
                no_header=False):
    n_in, n_out, n_exp = 0, 0, 0
    with open(infile) as fh, open(outfile, "w") as out:
        if not no_header:
            header = fh.readline()
            out.write(header)
            if column is not None:
                cols = header.rstrip("\n").split("\t")
                if column not in cols:
                    sys.exit(f"column '{column}' not in {infile}")
                ci = cols.index(column)
            else:
                ci = col_index
        else:
            ci = col_index
        for line in fh:
            n_in += 1
            f = line.rstrip("\n").split("\t")
            cid = norm_id(f[ci])
            frags = frag_map.get(cid)
            if not frags:
                out.write(line)
                n_out += 1
                continue
            n_exp += 1
            for frag in frags:
                f2 = list(f)
                f2[ci] = frag
                out.write("\t".join(f2) + "\n")
                n_out += 1
    print(f"  {os.path.basename(outfile)}: {n_in} rows in, {n_out} out, "
          f"{n_exp} parent rows expanded")


def main():
    ap = argparse.ArgumentParser(
        description="Fragment chromosomes to match a query contig, by motif "
                    "occurrences (--mode motif) or by length (--mode length).")
    ap.add_argument("--mode", choices=["motif", "length"], default="motif",
                    help="fragmentation mode [motif]")
    ap.add_argument("--motifs", default=None,
                    help="nanomotif or epimetheus motif table "
                         "(required for --mode motif; optional for --mode "
                         "length, then only used for the query summary)")
    ap.add_argument("--assembly", required=True, help="assembly FASTA")
    ap.add_argument("--mobsuite", required=True, help="MobSuite contig_report.txt")
    ap.add_argument("--kraken2", default=None, help="Kraken2 output (headerless)")
    ap.add_argument("--amrfinder", default=None, help="AMRFinder TSV")
    ap.add_argument("--amrfinder-column", default="Contig id",
                    help="contig id column in AMRFinder table [Contig id]")
    ap.add_argument("--contig", required=True, help="contig of interest")
    ap.add_argument("--min-obs", type=int, default=1,
                    help="for epimetheus tables: motif is present on a contig "
                         "only if n_motif_obs >= this [1]")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="length mode: a trailing piece shorter than this "
                         "fraction of the query length is merged into the "
                         "previous fragment [0.5]")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if args.mode == "motif" and not args.motifs:
        sys.exit("--mode motif requires --motifs")
    if not 0.0 <= args.min_frac <= 1.0:
        sys.exit("--min-frac must be between 0 and 1")

    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.contig)

    seqs = read_fasta(args.assembly)
    if args.contig not in seqs:
        sys.exit(f"{args.contig} not found in assembly")
    mob_types = read_mobsuite_types(args.mobsuite)
    query_seq = seqs[args.contig]
    query_len = len(query_seq)

    # --- motif counts for the query (needed in motif mode, summary otherwise) ---
    per_contig_motifs, per_contig_obs = {}, {}
    query_motifs, query_counts, use_obs = [], {}, False
    if args.motifs:
        per_contig_motifs, per_contig_obs = read_motif_table(args.motifs,
                                                             args.min_obs)
        if args.contig not in per_contig_motifs:
            if args.mode == "motif":
                sys.exit(f"{args.contig} not found in motif table")
            print(f"warning: {args.contig} not in motif table; "
                  f"no query summary", file=sys.stderr)
        else:
            query_motifs = sorted(per_contig_motifs[args.contig])
            use_obs = bool(per_contig_obs)
            if use_obs:
                query_counts = {m: per_contig_obs[args.contig].get(m, 0)
                                for m in query_motifs}
                count_src = "n_motif_obs"
            else:
                query_counts = {m: len(motif_positions(query_seq, m[0]))
                                for m in query_motifs}
                count_src = "sequence scan"
            for m, n in query_counts.items():
                if n == 0:
                    print(f"warning: {m[0]}_{m[1]}_{m[2]} has 0 {count_src} "
                          f"counts in {args.contig}", file=sys.stderr)

            qs_path = f"{prefix}_query_summary.tsv"
            with open(qs_path, "w") as qs:
                qs.write("query\tlength\tmotif\tn_occurrences\n")
                for m in query_motifs:
                    qs.write(f"{args.contig}\t{query_len}\t"
                             f"{m[0]}_{m[1]}_{m[2]}\t{query_counts[m]}\n")
                qs.write(f"{args.contig}\t{query_len}\tTOTAL\t"
                         f"{sum(query_counts.values())}\n")
            print(f"query {args.contig}: {query_len} bp, "
                  f"{len(query_motifs)} motifs, "
                  f"{sum(query_counts.values())} occurrences ({count_src}) "
                  f"-> {qs_path}")

    if args.mode == "motif":
        n_target = sum(query_counts.values())
        if n_target == 0:
            sys.exit(f"query {args.contig} has no counted motif sites")
        print(f"[motif mode] n_target = {n_target} occurrences")
    else:
        n_target = query_len
        if n_target == 0:
            sys.exit(f"query {args.contig} has length 0")
        print(f"[length mode] window = {n_target} bp, "
              f"min trailing piece = {int(args.min_frac * n_target)} bp")

    # --- fragmentation ---
    fa_path = f"{prefix}_fragmented.fasta"
    tab_path = f"{prefix}_fragments.tsv"
    frag_map = {}  # parent -> [fragment names], only for truly fragmented parents
    n_fragmented = n_no_overlap = n_too_small = 0
    mode = args.mode

    with open(fa_path, "w") as fa, open(tab_path, "w") as tab:
        tab.write("fragment\tparent\tstart\tend\tlength\tn_motifs\t"
                  "n_target\tmode\n")

        def write_seq(name, seq):
            fa.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")

        def write_row(name, parent, s, e, nm):
            tab.write(f"{name}\t{parent}\t{s}\t{e}\t{e - s}\t{nm}\t"
                      f"{n_target}\t{mode}\n")

        for contig, seq in seqs.items():
            if contig == args.contig or mob_types.get(contig) != "chromosome":
                write_seq(contig, seq)
                continue

            if mode == "length":
                frags = fragment_by_length(len(seq), n_target, args.min_frac)
                if len(frags) == 1:
                    # shorter than (or about) one query length -> left whole
                    n_too_small += 1
                    write_seq(contig, seq)
                    write_row(contig, contig, 0, len(seq), "NA")
                    continue
                n_fragmented += 1
                names = []
                for i, (s, e) in enumerate(frags, 1):
                    name = f"{contig}_{i}"
                    names.append(name)
                    write_seq(name, seq[s:e])
                    write_row(name, contig, s, e, "NA")
                frag_map[contig] = names
                continue

            # motif mode
            own = per_contig_motifs.get(contig, set())
            shared = own & set(query_motifs)
            # rule 1: no overlapping motif -> not fragmented
            if not shared:
                n_no_overlap += 1
                write_seq(contig, seq)
                continue
            # rule 2: chromosome's own total counts must reach n_target
            if use_obs:
                chrom_total = sum(per_contig_obs[contig].values())
            else:
                chrom_total = sum(len(motif_positions(seq, m[0])) for m in own)
            if chrom_total < n_target:
                n_too_small += 1
                write_seq(contig, seq)
                write_row(contig, contig, 0, len(seq), chrom_total)
                continue
            # rule 3: cut on all of the chromosome's own motifs,
            # every n_target occurrences
            positions = []
            for m in own:
                positions += motif_positions(seq, m[0])
            frags = fragment_by_motifs(seq, positions, n_target)
            if len(frags) == 1:
                write_seq(contig, seq)
                s, e, c = frags[0]
                write_row(contig, contig, s, e, c)
                continue
            n_fragmented += 1
            names = []
            for i, (s, e, c) in enumerate(frags, 1):
                name = f"{contig}_{i}"
                names.append(name)
                write_seq(name, seq[s:e])
                write_row(name, contig, s, e, c)
            frag_map[contig] = names

    if mode == "motif":
        print(f"{n_fragmented} chromosomes fragmented, "
              f"{n_no_overlap} skipped (no shared motif), "
              f"{n_too_small} skipped (total counts < n_target) -> {fa_path}")
    else:
        print(f"{n_fragmented} chromosomes fragmented, "
              f"{n_too_small} left whole (not longer than query) -> {fa_path}")

    # --- remapping ---
    remap_table(args.mobsuite, f"{prefix}_contig_report.frag.txt",
                frag_map, column="contig_id")
    if args.kraken2:
        remap_table(args.kraken2, f"{prefix}_kraken2.frag.output",
                    frag_map, col_index=1, no_header=True)
    if args.amrfinder:
        remap_table(args.amrfinder, f"{prefix}_amrfinder.frag.tsv",
                    frag_map, column=args.amrfinder_column)


if __name__ == "__main__":
    main()