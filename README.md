# Nanopore AMR-Carrying Plasmid–Host Association Pipeline

A **modular, end-to-end pipeline** that takes Oxford Nanopore sequencing data all the way from **raw electrical signals** to **AMR-carrying plasmid → host associations**.

Host association is inferred from **methylation-pattern similarity** (Nanomotif): each AMR-carrying plasmid is linked to the chromosome contig whose methylation profile it most closely matches. **Taxonomic classification** (Kraken2) is used **strictly for annotation** — never for distance calculation or inference — and **AMR genes** are detected with **AMRFinderPlus**.

Each step is a **standalone script**, but the steps are designed to run **in the order below**.

---

## Contents

- [Requirements](#requirements)
- [Input data](#input-data)
- [Running the pipeline](#running-the-pipeline)
  - [Phase 1 · Basecalling & preprocessing](#phase-1--basecalling--preprocessing-steps-15)
  - [Phase 2 · Assembly & polishing](#phase-2--assembly--polishing-steps-69)
  - [Phase 3 · Methylation & annotation](#phase-3--methylation--annotation-steps-1015)
  - [Phase 4 · Host association](#phase-4--host-association-step-16)
- [Output tables](#output-tables)
- [Read-level host association](#read-level-host-association)
- [Mock-community evaluation](#mock-community-evaluation)
- [Citation](#citation)

---

## Requirements

### External tools

| Tool | Role in the pipeline |
| --- | --- |
| `dorado` | Basecalling, read alignment, assembly polishing |
| `samtools` | BAM manipulation |
| `chopper` | Read quality & length filtering |
| `nanoMDBG`| Metagenome assembly |
| `modkit` | Methylation pileup |
| `nanomotif` | Methylation motif discovery & host association |
| `epimetheus` | Per-read methylation scoring (read-level route) |
| `mob_suite` | Plasmid / chromosome classification |
| `amrfinder` | AMR gene detection |
| `kraken2` | Taxonomic annotation |

### Python

- Python ≥ 3.9
- Packages: `pandas`, `numpy`, `pysam`
- Host-association scripts additionally use `scikit-bio` (PCoA) and `matplotlib` (figures); the read-level script also uses `numba` (JIT-compiled scoring)

---

## Input data

You typically start with:

- Raw Nanopore `.pod5` files
- The sequencing kit name (required for Dorado basecalling)
- A Kraken2 database
- *(Optional)* a reference FASTA for polishing

---

## Running the pipeline

Each step below lists its **script**, its **input → output**, and the command to run it.

### Phase 1 · Basecalling & preprocessing (steps 1–5)

**1 · Basecalling** *(GPU recommended)*
`basecall_dorado_sup.sh` — POD5 directory → basecalled BAM
```bash
bash basecall_dorado_sup.sh
```

**2 · Demultiplexing**
`demux.sh` — basecalled BAM → per-barcode BAMs
```bash
bash demux.sh
```

**3 · BAM → FASTQ**
`bam2fastq.sh` — BAM → FASTQ
```bash
bash bam2fastq.sh
```

**4 · Read quality & length filtering**
`chopper_filter.sh` — FASTQ → filtered FASTQ
```bash
bash chopper_filter.sh
```

**5 · Filter BAM by filtered read IDs**
`filter_bam_by_fastq.py` — filtered FASTQ + original BAM → filtered BAM
Ensures the BAM and FASTQ contain **exactly the same reads**.
```bash
python filter_bam_by_fastq.py
```

### Phase 2 · Assembly & polishing (steps 6–9)

**6 · Assembly** *(choose one assembler)*

*Option A — nanoMDBG*
`nanomdbg_assembly.sh` — filtered FASTQ → assembled contigs FASTA
```bash
bash nanomdbg_assembly.sh
```

**7 · Align reads to the unpolished assembly**
`dorado_align.sh` — filtered BAM + assembly FASTA → aligned BAM
```bash
bash dorado_align.sh
```

**8 · Polish the assembly**
`dorado_polish.sh` — aligned BAM + assembly FASTA → polished assembly FASTA
```bash
bash dorado_polish.sh
```

**9 · Align reads to the polished assembly**
`dorado_align.sh` — filtered BAM + polished assembly FASTA → aligned BAM
```bash
bash dorado_align.sh
```

### Phase 3 · Methylation & annotation (steps 10–15)

**10 · Modification pileup**
`modkit_pileup.sh` — aligned BAM + polished assembly → modification pileup BED
```bash
bash modkit_pileup.sh
```

**11 · Contigs as bins**
`make_contig_bins.py` — assembly FASTA → contig-bin TSV
Used to force **contig-level motif discovery**.
```bash
python make_contig_bins.py
```

**12 · Motif discovery (Nanomotif)**
`nanomotif.sh` — assembly FASTA + modification pileup BED + contig-bin TSV → `bin-motifs.tsv`
Discovers the methylation motifs. It is **contig-level** because each contig is fed as its own bin. `bin-motifs.tsv` (columns `bin`, `motif`, `mod_position`, `mod_type`, …) is the source of the motif list used by both host-association routes.
```bash
bash nanomotif.sh
```

**12b · Methylation scoring (epimetheus)**
`epimetheus4contig.sh` — sorted + bgzipped modkit pileup + polished assembly + motif list → `motifs-scored-read-methylation.tsv`
Builds the motif list from `bin-motifs.tsv` (joining `motif`, `mod_type`, `mod_position` into `motif_modtype_modposition`) and scores per-contig methylation with `epimetheus methylation-pattern contig`. This table is the `--nanomotif-dir` input to `similarity_scores.py` (step 16).
```bash
bash epimetheus4contig.sh
```

**13 · Plasmid / chromosome identification**
`mobsuite.sh` — assembly FASTA → `contig_report.txt`
```bash
bash mobsuite.sh
```

**14 · AMR detection**
`amrfinder.sh` — assembly FASTA → AMRFinderPlus results
```bash
bash amrfinder.sh
```

**15 · Taxonomic classification**
`kraken2.sh` — assembly FASTA → Kraken2 contig classification
```bash
bash kraken2.sh
```

### Phase 4 · Host association (step 16)

**16 · AMR host association** *(final inference step)*
`similarity_scores.py` — combines the Nanomotif, MobSuite, AMRFinderPlus and Kraken2 outputs into a single host-association table.

Each AMR-carrying plasmid is scored against every classified chromosome contig by methylation-vector similarity and assigned to the best-scoring host:

```
final_score = max(0, 1 − RMSD) × n_shared_motifs
```

where RMSD is the root-mean-square deviation between the two contigs' Nanomotif methylation values over their shared motifs. A higher `final_score` means a better host match.

**Principles**
- Uses only **AMR-carrying** plasmids as queries
- Plasmids are linked only to **chromosome** contigs
- Similarity is computed on **Nanomotif methylation vectors**
- Kraken2 taxonomy is used only for **annotation** — never to place or score contigs

```bash
python similarity_scores.py \
  --nanomotif-dir /path/to/nanomotif_output \
  --mobsuite-dir  /path/to/mobsuite_output \
  --amr-dir       /path/to/amrfinder_output \
  --kraken-dir    /path/to/kraken_output \
  --outdir        /path/to/final_results
```

The assignment table is written for **every** AMR-carrying plasmid on every run. To also generate per-contig figures for one or more specific contigs, add `--only-contig` (and optionally `--prefix` to tag outputs with a sample ID):

```bash
python similarity_scores.py \
  ... \
  --only-contig ctg123 \
  --prefix S1
```

For each requested contig this produces:
- a **top-N PCoA** plot (default `--topn-pcoa 100`) with a full-candidate inset, and
- a **top-5 species heatmap** — the plasmid next to one chromosome from each of the top 5 distinct host species.

---

## Output tables

### `AMR_plasmid_assignment.tsv`

One row per **AMR-carrying plasmid** — the single table produced on every run (prefixed with `--prefix` if supplied).

| Column | Description |
| --- | --- |
| `query_contig` | Plasmid contig ID |
| `molecule_type` | Always `plasmid` |
| `amr_genes` | Semicolon-separated AMR genes carried by the plasmid |
| `n_amr_genes` | Number of AMR genes |
| `best_host_contig` | Highest-scoring chromosome contig |
| `best_host_species` | Kraken2 species of the best host |
| `best_host_full_label` | Full Kraken2 label of the best host |
| `best_final_score` | `max(0, 1 − RMSD) × n_shared_motifs` for the best host |
| `best_rmsd` | RMSD to the best host over shared motifs |
| `best_shared_motifs` | Number of motifs shared with the best host |
| `n_candidates` | Chromosome candidates scored |
| `winner_species_top_ties` | Consensus species among score-tied top hits |
| `tie_species_summary` | Per-species counts among the tied top hits |

### Optional per-contig figures (with `--only-contig`)

- `pcoa_top{N}/per_contig/{plasmids,chromosomes}/` — top-N PCoA with a full-candidate inset (PNG + PDF)
- `top_taxa_heatmap/per_contig/{plasmids,chromosomes}/` — top-5 species methylation heatmap (PNG + PDF)

A persistent `species_colors.json` is also written/updated so taxa keep consistent colours across runs (override the path with `--color-map`).

---

## Read-level host association

`read_similarity_scores.py` — a **per-read** alternative to step 16. Instead of assembled contigs, it works directly on **individual reads**: each AMR-carrying read is scored against the pool of classified chromosome reads by methylation similarity, using the same score

```
final_score = max(0, 1 − RMSD) × n_shared_motifs
```

computed on the **mean** methylation probability per motif. A numba-compiled inner loop and a multiprocessing pool let it scale to very large read sets. Use this flow when MobSuite / AMRFinderPlus / Kraken2 were run **per read** (quasi-metagenomic reads classified individually) rather than on the assembly.

### Generating the per-read methylation table

The read-level route needs a **per-read** methylation table (the `--read-motif-dir` input). This is produced by **epimetheus** from the MM/ML base-modification tags in the filtered BAM, scored against the motifs discovered by **Nanomotif**.

**1 · Motif list** — Nanomotif's motif-discovery step (step 12) writes `bin-motifs.tsv` with `motif`, `mod_type` and `mod_position` columns. epimetheus wants these joined as `motif_modtype_modposition`, one per line:

```
AAGNNNNNGTNG_a_1
AGCACC_a_3
ANNANC_a_3
ATNNC_a_0
CACNNNNNTANG_a_1
CCAG_m_1
CCCGGG_m_0
```

(`a` = 6mA, `m` = 5mC; the trailing number is the modified position within the motif.)

**2 · epimetheus** — `epimetheus4read.sh` builds that list from `bin-motifs.tsv` (skip this if you already have `motifs.txt`), extracts the MM/ML reads from `filtered.bam`, and runs `epimetheus methylation-pattern read-fastq`:

```bash
bash epimetheus4read.sh
```

→ `epimetheus/epimetheus_main.tsv` — the per-read methylation table, pointed to by `--read-motif-dir`.

### Inputs

| Flag | Contents |
| --- | --- |
| `--read-motif-dir` | epimetheus per-read methylation table (`read_id`, `motif_mod_position`, `mean_prob`), e.g. `epimetheus/epimetheus_main.tsv`; `.tsv` or `.csv` |
| `--read-types` | MobSuite `contig_report.txt` giving each read a `molecule_type` (chromosome / plasmid) |
| `--amr-dir` | AMRFinderPlus results (reads the `Element symbol` and `Element name` columns) |
| `--kraken-dir` | Kraken2 per-read `*.out` and `*.report` |
| `--outdir` | Output directory |

```bash
python read_similarity_scores.py \
  --read-motif-dir /path/to/read_motif \
  --read-types     /path/to/mobsuite/contig_report.txt \
  --amr-dir        /path/to/amrfinder \
  --kraken-dir     /path/to/kraken \
  --outdir         /path/to/read_results \
  --dataset        M3 \
  --workers        8
```

### Per-read association table

**`AMR_plasmid_read_assignment_summary.tsv`** — one row per AMR-carrying **plasmid** read (plasmid → host chromosome). Chromosome-vs-chromosome associations are not produced.

Columns:

| Column | Description |
| --- | --- |
| `read_id` | Read identifier |
| `molecule_type` | Always `plasmid` |
| `amr_genes` | Semicolon-separated AMR genes on the read |
| `has_nanomotif_vector` | Whether the read has a methylation vector |
| `kraken_raw` | Kraken2 label |
| `kraken_rank_code` | Kraken2 rank code |
| `mean_top_score_species_counts` | Species tally of the read's top-scoring host reads |

### Carbapenemase summary

**`carbapenemase_read_host_summary.tsv`** — one row per carbapenemase gene, aggregating the top host of every **plasmid** read that carries it (chromosome-borne carbapenemase reads are excluded, since those would be chromosome-vs-chromosome):

| Column | Description |
| --- | --- |
| `dataset` | Sample label (`--dataset`; defaults to the outdir name) |
| `carbapenemase` | Carbapenemase gene, e.g. `blaOXA-48` |
| `n_reads` | Number of reads carrying that gene |
| `top_read_similarity` | Tally of each read's top host species, e.g. `Citrobacter freundii (1), Citrobacter koseri (1), Obesumbacterium proteus (1)` |

Carbapenemases are detected from **both** the AMRFinderPlus `Element symbol` and `Element name` columns, so an OXA carbapenemase reported at family level (symbol `blaOXA`, allele `OXA-48` only in the name) is captured as `blaOXA-48`, while non-carbapenemase oxacillinases such as `blaOXA-1` are excluded. The recognised alleles are controlled by `--carbapenemase-regex`.

### Optional per-read figures (with `--read-ids`)

By default only the tables are written. Pass `--read-ids read1,read2` (or `--read-ids-file`) to also produce, for each requested read:
- a **top-N PCoA** plot (`--topn-pcoa`, default 100), and
- a **top-5 species heatmap** (`--top-taxa-heatmap`) — the read next to one chromosome from each of the top 5 distinct host species.

When `--read-ids` is given the association tables are restricted to those reads (the fast path); add `--associate-all` to score every AMR read while still limiting figures to the requested reads.

---

## Mock-community evaluation

A separate workflow that benchmarks Nanomotif-based host association against **known isolate composition**. The mock community is built by pooling isolate datasets; read-level isolate labels are mapped to contigs by alignment to produce **contig-level ground truth**.

### Step 1 · Read-level ground truth

Extract `read_id → isolate_label` from the **unaligned Dorado basecalling BAMs** (one BAM per isolate/species barcode; the filename stem is the isolate label).

→ `read_isolate.tsv`

### Step 2 · Read → contig mapping

Align reads to the **mock assembly** and keep **primary alignments only** (secondary, supplementary and unmapped reads are discarded).

→ `read_to_ctg.tsv`

### Step 3 · Contig-level ground truth

Aggregate read labels per contig by joining `read_isolate.tsv` with `read_to_ctg.tsv`.

```bash
python build_contig_gt.py \
  --read-isolate  /path/to/read_isolate.tsv \
  --read-ctg      /path/to/read_to_contig_map.tsv \
  --mobsuite-root /path/to/mobsuite \
  --kraken-root   /path/to/kraken \
  --out           /path/to/contig_eval_table.tsv
```

→ `contig_ground_truth.tsv`

### Step 4 · AMR host association

Run `similarity_scores.py` on the mock assembly to generate Nanomotif-based host predictions.

### Step 5 · Evaluation (Nanomotif vs ground truth)

Extract the top taxon from `nanomotif_taxonomic_association_profile` and compare it against the ground truth. Comparison collapses to species level, e.g.:

```
Escherichia coli O157                    → Escherichia coli
Klebsiella pneumoniae subsp. pneumoniae  → Klebsiella pneumoniae
```

→ `plasmid_nanotax_vs_gt.tsv`

### Input files and formats

**Dorado basecalling BAMs (unaligned)** — filename stem = isolate label:
```
ecoli_37_38.bam
smarcescens_35_36.bam
```

**Reads FASTQ:** `reads.fastq.gz`
**Mock assembly FASTA:** `mock_assembly.fasta`

**Nanomotif plasmid host association TSV** — required columns `ctg_id`, `nanomotif_taxonomic_association_profile`, e.g.:
```
Citrobacter(18); Citrobacter freundii(14); Bacteria(1)
```

**Contig-level ground truth TSV** — required columns: `ctg_id`, plus one of `gt_isolate_top` / `gt_species` / `isolate_label`.

### Output file definitions

**`read_isolate.tsv`**

| Column | Description |
| --- | --- |
| `read_id` | Read identifier |
| `isolate_label` | Isolate inferred from BAM filename |

**`read_to_ctg.tsv`** *(primary alignments only)*

| Column | Description |
| --- | --- |
| `read_id` | Read identifier |
| `ctg_id` | Contig ID |

**`contig_ground_truth.tsv`**

| Column | Description |
| --- | --- |
| `ctg_id` | Contig identifier |
| `gt_isolate_top` | Dominant isolate |
| `gt_isolate_top_pct` | Fraction supporting the dominant isolate |
| `gt_isolate_dist` | Full isolate distribution |
| `gt_isolate_n` | Number of reads |

**`plasmid_nanotax_vs_gt.tsv`**

| Column | Description |
| --- | --- |
| `ctg_id` | Plasmid contig ID |
| `gt_species` | Ground-truth species |
| `nanomotif_taxon` | Nanomotif-predicted taxon |
| `correct` | Boolean correctness flag |

---

## Citation

If you use this pipeline, please cite:

> *(citation to be added)*