#!/usr/bin/env python3
"""
CUPID contig fragmentation, all in one.

Given a contig of interest, this script:
  1. Fragments MobSuite chromosomes so each fragment carries the same total
     number of shared-motif occurrences as the query contig (n_target =
     occurrences, in the query sequence, of motifs present on BOTH the
     query's and that chromosome's nanomotif lists). Chromosomes sharing no
     motif with the query are left whole. Plasmids, unclassified contigs and
     the query are never fragmented.
  2. Writes the full fragmented assembly (every input sequence exactly once,
     fragments named parent_1, parent_2, ...).
  3. Remaps the MobSuite contig report, Kraken2 output and AMRFinder table to
     fragment names by duplicating each fragmented parent's rows per fragment
     (mapping from the fragment table, no suffix parsing).

Outputs in --outdir (prefix = query contig name):
  <prefix>_fragmented.fasta
  <prefix>_fragments.tsv                fragment, parent, start, end, length,
                                        n_motifs, n_target
  <prefix>_contig_report.frag.txt
  <prefix>_kraken2.frag.output          (if --kraken2 given)
  <prefix>_amrfinder.frag.tsv           (if --amrfinder given)

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
    ap = argparse.ArgumentParser()
    ap.add_argument("--motifs", required=True, help="nanomotif contig motifs TSV")
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
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.contig)

    per_contig_motifs, per_contig_obs = read_motif_table(args.motifs, args.min_obs)
    if args.contig not in per_contig_motifs:
        sys.exit(f"{args.contig} not found in motif table")

    seqs = read_fasta(args.assembly)
    if args.contig not in seqs:
        sys.exit(f"{args.contig} not found in assembly")

    mob_types = read_mobsuite_types(args.mobsuite)

    query_seq = seqs[args.contig]
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
            print(f"warning: {m[0]}_{m[1]}_{m[2]} has 0 {count_src} counts "
                  f"in {args.contig}", file=sys.stderr)
    n_target = sum(query_counts.values())
    if n_target == 0:
        sys.exit(f"query {args.contig} has no counted motif sites")

    # query summary: length and motif counts
    qs_path = f"{prefix}_query_summary.tsv"
    with open(qs_path, "w") as qs:
        qs.write("query\tlength\tmotif\tn_occurrences\n")
        for m in query_motifs:
            qs.write(f"{args.contig}\t{len(query_seq)}\t"
                     f"{m[0]}_{m[1]}_{m[2]}\t{query_counts[m]}\n")
        qs.write(f"{args.contig}\t{len(query_seq)}\tTOTAL\t"
                 f"{sum(query_counts.values())}\n")
    print(f"query {args.contig}: {len(query_seq)} bp, "
          f"{len(query_motifs)} motifs, n_target = {n_target} "
          f"({count_src}) -> {qs_path}")

    # --- fragmentation ---
    fa_path = f"{prefix}_fragmented.fasta"
    tab_path = f"{prefix}_fragments.tsv"
    frag_map = {}  # parent -> [fragment names], only for truly fragmented parents
    n_fragmented = n_no_overlap = n_low_obs = 0

    with open(fa_path, "w") as fa, open(tab_path, "w") as tab:
        tab.write("fragment\tparent\tstart\tend\tlength\tn_motifs\tn_target\n")

        def write_seq(name, seq):
            fa.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")

        for contig, seq in seqs.items():
            if contig == args.contig or mob_types.get(contig) != "chromosome":
                write_seq(contig, seq)
                continue
            own = per_contig_motifs.get(contig, set())
            shared = own & set(query_motifs)
            # rule 1: no overlapping motif -> not fragmented (and not a
            # candidate in the assignment anyway)
            if not shared:
                n_no_overlap += 1
                write_seq(contig, seq)
                continue
            # rule 2: chromosome's own total counts must reach n_target,
            # otherwise it cannot fill a single fragment -> left whole
            # (still a candidate in the assignment script)
            if use_obs:
                chrom_total = sum(per_contig_obs[contig].values())
            else:
                chrom_total = sum(len(motif_positions(seq, m[0])) for m in own)
            if chrom_total < n_target:
                n_low_obs += 1
                write_seq(contig, seq)
                tab.write(f"{contig}\t{contig}\t0\t{len(seq)}\t{len(seq)}\t"
                          f"{chrom_total}\t{n_target}\n")
                continue
            # rule 3: cut on ALL of the chromosome's own motifs, every
            # n_target occurrences (n_target = query total, fixed)
            positions = []
            for m in own:
                positions += motif_positions(seq, m[0])
            frags = fragment_by_motifs(seq, positions, n_target)
            if len(frags) == 1:
                write_seq(contig, seq)
                s, e, c = frags[0]
                tab.write(f"{contig}\t{contig}\t{s}\t{e}\t{e - s}\t{c}\t{n_target}\n")
                continue
            n_fragmented += 1
            names = []
            for i, (s, e, c) in enumerate(frags, 1):
                name = f"{contig}_{i}"
                names.append(name)
                write_seq(name, seq[s:e])
                tab.write(f"{name}\t{contig}\t{s}\t{e}\t{e - s}\t{c}\t{n_target}\n")
            frag_map[contig] = names

    print(f"{n_fragmented} chromosomes fragmented, "
          f"{n_no_overlap} skipped (no shared motif), "
          f"{n_low_obs} skipped (total counts < n_target) -> {fa_path}")

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