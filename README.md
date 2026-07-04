# Nanopore AMR Host Association Pipeline

This repository provides a **modular, end-to-end pipeline** for processing Nanopore sequencing data from **raw electrical signals** to **AMR-carrying plasmid and chromosome host association**.

Host association is inferred using **methylation pattern similarity (Nanomotif)**, while **taxonomic classification (Kraken2)** is used strictly for **annotation**, not for distance or inference. AMR genes are detected using **AMRFinder**.

Each step is implemented as a **standalone script**, but the steps are intended to be run **in the order described below**.

---

## Requirements

### External tools

You should have the following tools available in your environment:

* `dorado`
* `samtools`
* `chopper`
* `metaMDBG` **or** `flye` (metaFlye)
* `modkit`
* `nanomotif`
* `mob_suite`
* `amrfinder`
* `kraken2`

### Python

* Python ≥ 3.9
* Required packages:

  * `pandas`
  * `numpy`
  * `pysam`

---

## Input data

You typically start with:

* Raw Nanopore `.pod5` files
* Sequencing kit name (for Dorado basecalling)
* Kraken2 database
* (Optional) reference FASTA for polishing

---

### 1. Basecalling (GPU recommended)

**Script:** `basecall_dorado_sup.sh`
**Input:** `.pod5` directory
**Output:** basecalled BAM

```bash
bash basecall_dorado_sup.sh
```

---

### 2. Demultiplexing

**Script:** `demux.sh`
**Input:** basecalled BAM
**Output:** per-barcode BAMs

```bash
bash demux.sh
```

---

### 3. BAM → FASTQ

**Script:** `bam2fastq.sh`
**Input:** BAM
**Output:** FASTQ

```bash
bash bam2fastq.sh
```

---

### 4. Read quality & length filtering

**Script:** `chopper_filter.sh`
**Input:** FASTQ
**Output:** filtered FASTQ

```bash
bash chopper_filter.sh
```

---

### 5. Filter BAM using filtered FASTQ read IDs

Ensures that BAM and FASTQ contain **exactly the same reads**.

**Script:** `filter_bam_by_fastq.py`
**Input:** filtered FASTQ + original BAM
**Output:** filtered BAM

```bash
python filter_bam_by_fastq.py
```

---

### 6. Assembly

Choose **one** assembler.

#### Option A — metaMDBG

**Script:** `nanomdbg_assembly.sh`

```bash
bash nanomdbg_assembly.sh
```

#### Option B — metaFlye

**Script:** `metaflye_assembly.sh`

```bash
bash metaflye_assembly.sh
```

**Output:** assembled contigs FASTA

---

### 7. Align reads to unpolished assembly

**Script:** `dorado_align.sh`
**Input:** filtered BAM + assembly FASTA
**Output:** aligned BAM

```bash
bash dorado_align.sh
```

---

### 8. Assembly polishing

**Script:** `dorado_polish.sh`
**Input:** aligned BAM + assembly FASTA
**Output:** polished assembly FASTA

```bash
bash dorado_polish.sh
```

---

### 9. Align reads to polished assembly

**Script:** `dorado_align.sh`
**Input:** filtered BAM + polished assembly FASTA
**Output:** aligned BAM

```bash
bash dorado_align.sh
```

---

### 10. Modification pileup

**Script:** `modkit_pileup.sh`
**Input:** aligned BAM + polished assembly
**Output:** modification pileup BED

```bash
bash modkit_pileup.sh
```

---

### 11. Contigs as bins

Used forcing **contig-level motif discovery**.

**Script:** `make_contig_bins.py`
**Input:** assembly FASTA
**Output:** contig-bin TSV

```bash
python make_contig_bins.py
```

---

### 12. Nanomotif analysis

**Script:** `nanomotif.sh`
**Input:**

* assembly FASTA
* modification pileup BED
* contig-bin TSV

**Output includes:**

* `motifs-scored-read-methylation.tsv`

```bash
bash nanomotif.sh
```

---

### 13. Plasmid / chromosome identification

**Script:** `mobsuite.sh`
**Input:** assembly FASTA
**Output:** `contig_report.txt`

```bash
bash mobsuite.sh
```

---

### 14. AMR detection

**Script:** `amrfinder.sh`
**Input:** assembly FASTA
**Output:** AMRFinder results

```bash
bash amrfinder.sh
```

---

### 15. Taxonomic classification

**Script:** `kraken2.sh`
**Input:** assembly FASTA
**Output:** Kraken2 contig classification

```bash
bash kraken2.sh
```

---

### 16. AMR host association (final inference step)

**Script:** `amr_host_association.py`

**Key principles:**

* Uses only AMR-carrying plasmid
* Plasmids are mapped only to chromosome contigs
* RMSD is computed on Nanomotif methylation vectors
* Kraken taxonomy is used only for annotation

**Inputs:**

* Nanomotif output directory
* MobSuite output directory
* AMRFinder output directory
* Kraken2 output directory

**Outputs:**

* AMR plasmid host association table

```
python amr_host_association.py \
  --nanomotif-dir /path/to/nanomotif_output \
  --mobsuite-dir  /path/to/mobsuite_output \
  --amr-dir       /path/to/amrfinder_output \
  --kraken-dir    /path/to/kraken_output \
  --outdir        /path/to/final_results
```

---

## Final output tables

### AMR plasmid host association

* One row per **AMR-carrying plasmid**
* Host inferred via **plasmid → chromosome Nanomotif RMSD**

#### More information: docs/amr_host_association_outputs.md

---

# Mock community pipeline: contig-level ground truth and evaluation

This section describes a **mock-community evaluation workflow** used to benchmark Nanomotif-based host association against **known isolate composition**.

The mock community is generated by pooling isolate datasets. Read-level isolate labels are mapped to contigs via alignment, producing contig-level ground truth.

---

## Pipeline overview

### Step 1 — Read-level ground truth

Extract `read_id → isolate_label` from **unaligned Dorado basecalling BAMs**.

* Each BAM corresponds to one isolate (or species barcode)
* BAM filename stem is treated as the isolate label

**Output**

* `read_isolate.tsv`

---

### Step 2 — Read → contig mapping

Align reads to the **mock assembly** and extract **primary alignments only**.

* Secondary, supplementary, and unmapped reads are discarded

**Output**

* `read_to_ctg.tsv`

---

### Step 3 — Contig-level ground truth generation

Generate contig-level host ground truth by aggregating read labels.

* Join:

  * `read_isolate.tsv`
  * `read_to_ctg.tsv`

python /home/haicu/harika.uerel/github_nanomotif/build_contig_gt.py \
 --read-isolate /path/to/read_isolate.tsv \
 --read-ctg /path/to/read_to_contig_map.tsv \
 --mobsuite-root /path/to/mobsuite \
 --kraken-root /path/to/kraken \
 --out /path/to/contig_eval_table.tsv

**Output**

* `contig_ground_truth.tsv`

---

### Step 4 — AMR host association (Nanomotif inference)

Run `amr_host_association.py` on the mock assembly to generate Nanomotif-based host predictions.

---

### Step 5 — Evaluation (Nanomotif vs ground truth)

Evaluate Nanomotif host association against contig-level ground truth.

* Extract top taxon from:

  * `nanomotif_taxonomic_association_profile`
* Compare to GT

Examples:

```
Escherichia coli O157 → Escherichia coli
Klebsiella pneumoniae subsp. pneumoniae → Klebsiella pneumoniae
```

**Output**

* `plasmid_nanotax_vs_gt.tsv`

---

## Input files and formats

### Dorado basecalling BAMs (unaligned)

```
ecoli_37_38.bam
smarcescens_35_36.bam
```

Assumption: BAM filename stem = isolate label.

---

### Reads FASTQ

```
reads.fastq.gz
```

---

### Mock assembly FASTA

```
mock_assembly.fasta
```

---

### Nanomotif plasmid host association TSV

Required columns:

* `ctg_id`
* `nanomotif_taxonomic_association_profile`

Example:

```
Citrobacter(18); Citrobacter freundii(14); Bacteria(1)
```

---

### Contig-level ground truth TSV

Required:

* `ctg_id`
* one of:

  * `gt_isolate_top`
  * `gt_species`
  * `isolate_label`

---

## Output file definitions

### `read_isolate.tsv`

| Column          | Description                        |
| --------------- | ---------------------------------- |
| `read_id`       | Read identifier                    |
| `isolate_label` | Isolate inferred from BAM filename |

---

### `read_to_ctg.tsv`

| Column    | Description     |
| --------- | --------------- |
| `read_id` | Read identifier |
| `ctg_id`  | Contig ID       |

Primary alignments only.

---

### `contig_ground_truth.tsv`

| Column               | Description                          |
| -------------------- | ------------------------------------ |
| `ctg_id`             | Contig identifier                    |
| `gt_isolate_top`     | Dominant isolate                     |
| `gt_isolate_top_pct` | Fraction supporting dominant isolate |
| `gt_isolate_dist`    | Full isolate distribution            |
| `gt_isolate_n`       | Number of reads                      |

---

### `plasmid_nanotax_vs_gt.tsv`

| Column            | Description                          |
| ----------------- | ------------------------------------ |
| `ctg_id`          | Plasmid contig ID                    |
| `gt_species`      | Ground-truth species                 |
| `nanomotif_taxon` | Nanomotif-predicted taxon            |
| `correct`         | Boolean correctness flag             |

