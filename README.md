# adVNTR — Modernized Fork

This is a modernized fork of [adVNTR](https://github.com/mehrdadbakhtiari/adVNTR), a tool for genotyping Variable Number Tandem Repeats (VNTRs) from short-read (Illumina) and long-read (PacBio, Nanopore) sequencing data using Hidden Markov Models.

**What this fork fixes:** The original adVNTR 1.5.0 fails to run on modern Python environments due to incompatible dependency APIs. This fork updates those API calls without changing any core algorithm or output format.

---

## Changes from upstream

| Component | Issue | Fix |
|-----------|-------|-----|
| `pomegranate/hmm.pyx` | `networkx.topological_sort()` no longer accepts `nbunch=` argument (removed in networkx ≥ 2.6) | Removed deprecated `nbunch` keyword |
| `advntr/vntr_graph.py` | `nx.connected_component_subgraphs()` removed in networkx 2.4 | Replaced with `(G.subgraph(c) for c in nx.connected_components(G))` |
| `advntr/vntr_graph.py` | `from networkx import graphviz_layout` removed in networkx ≥ 2.x | Replaced with `nx.drawing.nx_agraph.graphviz_layout` with fallback to spring layout |
| `advntr/vntr_finder.py`, `advntr/reference_vntr.py`, `advntr/pattern_clustering.py`, `advntr/pairwise_aln_generator.py` | `Bio.pairwise2` deprecated and removed in Biopython ≥ 1.80 | Migrated all calls to `Bio.Align.PairwiseAligner` |
| `requirements.txt` / `setup.py` | Pinned `networkx==1.11`, `biopython==1.76`, included `enum34` (built-in since Python 3.4) | Updated to `networkx>=2.6`, `biopython>=1.76`, removed `enum34` |

---

## Installation

**Requirements:** Python 3.8+, conda recommended.

```bash
# Clone this fork
git clone https://github.com/VenithaB/adVNTR.git
cd adVNTR

# Create and activate environment
conda create -n advntr_env python=3.11
conda activate advntr_env

# Install dependencies
pip install -r requirements.txt
```

### Running without installing

Since this fork is not published to PyPI, run it directly by setting `PYTHONPATH`:

```bash
PYTHONPATH=/path/to/adVNTR conda run -n advntr_env python -m advntr genotype [options]
```

This takes precedence over any installed `advntr` package in the environment.

---

## VNTR database

A pre-built database of VNTR models is required. It is not included in this repository.

| Database | Loci | Reference | Format |
|----------|------|-----------|--------|
| `hg38_selected_VNTRs_Illumina.db` | 10,264 | GRCh38 | SQLite `.db` |
| `vntr_data_recommended_loci.zip` | 6,719 (Illumina) / 8,960 (PacBio) | hg19 | Download from upstream |
| `vntr_data_recommended_loci_hg38.zip` | — | GRCh38 | Download from upstream |

Pass the database path with `-m/--models`.

---

## Input files

| File | Flag | Format | Notes |
|------|------|--------|-------|
| Alignment file | `-a/--alignment_file` | SAM, BAM, or CRAM | Must be coordinate-sorted and indexed (`.bai` / `.crai`). **Do not use `-f/--fasta` for alignment files** — that flag is for raw FASTA reads only. |
| Reference genome | `-r/--reference_filename` | FASTA (`.fa` / `.fasta`) | **Required for CRAM input.** Optional for BAM/SAM. |
| VNTR models database | `-m/--models` | SQLite `.db` | Defaults to `vntr_data/hg19_selected_VNTRs_Illumina.db` if not specified. |

---

## Usage

### Basic Illumina genotyping (BAM)

```bash
advntr genotype \
    --alignment_file sample.bam \
    --models /path/to/hg38_selected_VNTRs_Illumina.db \
    --working_directory ./working_dir/ \
    --outfile results.vcf \
    --outfmt vcf
```

### CRAM input (reference required)

```bash
advntr genotype \
    --alignment_file sample.cram \
    --reference_filename /path/to/GRCh38.fa \
    --models /path/to/hg38_selected_VNTRs_Illumina.db \
    --working_directory ./working_dir/ \
    --outfile results.vcf \
    --outfmt vcf
```

### Genotype a single locus

```bash
advntr genotype \
    --alignment_file sample.bam \
    --models /path/to/hg38_selected_VNTRs_Illumina.db \
    --working_directory ./working_dir/ \
    --vntr_id 201
```

### PacBio long reads

```bash
advntr genotype \
    --alignment_file sample_pacbio.bam \
    --models /path/to/hg38_selected_VNTRs_Illumina.db \
    --working_directory ./working_dir/ \
    --pacbio
```

---

## CLI reference — `advntr genotype`

### Input/output options

| Flag | Description |
|------|-------------|
| `-a/--alignment_file <file>` | Alignment file in SAM/BAM/CRAM format |
| `-r/--reference_filename <file>` | Reference FASTA — required for CRAM, overrides path in CRAM header |
| `-f/--fasta <file>` | Raw reads in FASTA format (not for alignment files) |
| `-o/--outfile <file>` | Output file path. Defaults to stdout if not specified. |
| `-of/--outfmt <format>` | Output format: `text`, `bed`, or `vcf` (default: `text`) |
| `--disable_logging` | Suppress log file output except critical errors |

### Algorithm options

| Flag | Description |
|------|-------------|
| `-p/--pacbio` | Input contains PacBio reads |
| `-n/--nanopore` | Input contains Nanopore MinION reads |
| `--accuracy_filter` | Genotype using spanning reads only (higher precision, more `None` calls) |
| `-e/--expansion` | Detect long expansions from PCR-free data (requires `-c/--coverage`) |
| `-c/--coverage <float>` | Average sequencing coverage for expansion detection |
| `--haploid` | Organism is haploid |
| `-fs/--frameshift` | Search for frameshifts instead of copy number. Supported VNTR IDs: 25561, 519759 |
| `-naive/--naive` | Use naive approach for PacBio reads |
| `-u/--update` | Iteratively update the VNTR model |

### Other options

| Flag | Description |
|------|-------------|
| `--working_directory <path>` | Directory for temporary files (log, intermediate FASTA, filtered reads) |
| `-m/--models <file>` | VNTR models database (default: `vntr_data/hg19_selected_VNTRs_Illumina.db`) |
| `-t/--threads <int>` | Number of threads (default: 1). **Note: only affects PacBio mode — Illumina genotyping is single-threaded regardless of this value.** |
| `-vid/--vntr_id <text>` | Comma-separated list of VNTR IDs to genotype. Runs all loci in the database if not specified. |

---

## Output formats

### `text` (default)
Tab-separated lines, one per locus:
```
VID     Chromosome  NumberOfSupportingReads  MeanCoverage  Genotype
201     chr1        30                       28.5          2/2
```

### `vcf`
Standard VCF 4.1 with a header block. Each record encodes the VNTR genotype in the `RU` INFO field (repeat unit count).

### `bed`
BED-format coordinates with genotype appended as extra columns.

---

## Working directory contents

After a run, the working directory contains:

| File | Description |
|------|-------------|
| `log_<sample>.log` | Detailed per-locus log. See `docs/README_log_format.md` for a full description of every message type. |
| `<sample>.unmapped.fasta` | Unmapped reads extracted from the alignment file |
| `keywords_<sample>.unmapped.fasta.txt` | Keyword index used for read pre-filtering |
| `filtering_out_<sample>.unmapped.fasta.txt` | Filtered read IDs per locus |

---

## Performance notes

- **Illumina, all loci (~10k):** ~4 hours on a single core for a 30× CRAM
- **Per-locus breakdown:** ~15 min unmapped read extraction, ~90 sec keyword filtering, then ~1 sec/locus HMM genotyping
- **Threads (`-t`)** only speed up PacBio mode; increasing thread count does not reduce Illumina runtime
- Memory usage is modest (~4–8 GB for a 30× genome)

---

## Running on a SLURM cluster

```bash
#!/bin/bash
#SBATCH --job-name=advntr_sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/path/to/output/advntr_%j.out
#SBATCH --error=/path/to/output/advntr_%j.err

ADVNTR_SRC="/path/to/adVNTR"

source /path/to/conda/etc/profile.d/conda.sh
conda activate advntr_env

PYTHONPATH="$ADVNTR_SRC" python -m advntr genotype \
    --alignment_file sample.cram \
    --reference_filename /path/to/GRCh38.fa \
    --models /path/to/hg38_selected_VNTRs_Illumina.db \
    --working_directory ./working_dir/ \
    --outfile results.vcf \
    --outfmt vcf
```

> Use `--ntasks=1 --cpus-per-task=<N>` (not `--ntasks=<N>`) for threaded jobs.
> The `--output` and `--error` log directories must exist **before** `sbatch` is called.

---

## Citation

- **Original adVNTR:**
  Bakhtiari et al. [Targeted genotyping of variable number tandem repeats with adVNTR](https://genome.cshlp.org/content/28/11/1709). *Genome Research* 28(11), 2018.

- **adVNTR-NN (v1.4.0):**
  Bakhtiari et al. [Variable Number Tandem Repeats mediate the expression of proximal genes](https://doi.org/10.1038/s41467-021-22206-z). *Nature Communications* 12, 2075, 2021.

- **code-adVNTR:**
  Park et al. [Detecting tandem repeat variants in coding regions using code-adVNTR](https://doi.org/10.1016/j.isci.2022.104785). *iScience* 25(8), 104785, 2022.
