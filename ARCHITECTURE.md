# adVNTR Architecture

This document explains how adVNTR works internally — what the main files are, how data flows through the tool, and what each layer does. It is intended as a living reference; update it as your understanding deepens.

---

## What adVNTR does

adVNTR genotypes **Variable Number Tandem Repeats (VNTRs)** from sequencing data. Given a BAM file of aligned reads and a database of known VNTR loci, it answers: *how many copies of each repeat unit does this individual carry at each VNTR?*

For diploid organisms the output is a pair of allele copy numbers (e.g., `3/5`), one per haplotype. For haploid organisms it outputs a single count.

---

## Entry point

**`advntr/__main__.py`** — `main()`

The CLI is built with `argparse`. Four subcommands are defined:

| Command | Purpose |
|---------|---------|
| `genotype` | Core command — find repeat counts and mutations in VNTRs |
| `viewmodel` | Inspect VNTR entries in the model database |
| `addmodel` | Add a custom VNTR locus to the database |
| `delmodel` | Remove a VNTR from the database |

`genotype` is the primary command. The others are database management utilities.

---

## CLI → execution flow

```
advntr/__main__.py          # parse args, dispatch
    └── advntr/advntr_commands.py   # genotype()
            ├── load VNTR models from SQLite (models.py)
            └── advntr/genome_analyzer.py   # GenomeAnalyzer
                    └── advntr/vntr_finder.py   # VNTRFinder (one per VNTR)
                            ├── read recruitment
                            ├── HMM construction  ← advntr/hmm_utils.py
                            │       └── pomegranate/hmm.pyx  (vendored HMM library)
                            ├── Viterbi decoding
                            └── statistical genotyping
```

---

## Layer-by-layer description

### 1. `advntr/advntr_commands.py` — Command dispatch

`genotype()` is the top-level function for the main workflow. It:
- Configures logging (written to `working_dir/log_<input>.log`)
- Sets the error rate for the read type (`--pacbio`/`--nanopore` → 0.30; Illumina → 0.05)
- Loads the VNTR model database (SQLite file pointed to by `settings.TRAINED_MODELS_DB`)
- Constructs a `GenomeAnalyzer` and calls the appropriate analysis method

The branch chosen depends on input type:

| Condition | Method called |
|-----------|--------------|
| Illumina BAM/SAM/CRAM | `find_repeat_counts_from_alignment_file()` |
| PacBio BAM | `find_repeat_counts_from_pacbio_alignment_file()` |
| Illumina, frameshift mode | `find_frameshift_from_alignment_file()` |

---

### 2. `advntr/genome_analyzer.py` — `GenomeAnalyzer`

Orchestrates the run across all target VNTRs. On construction it builds one `VNTRFinder` per target VNTR. Then for each VNTR it:
1. Calls the appropriate `VNTRFinder` method
2. Catches errors and logs them (this is where "Skipping genotyping for this VNTR" messages come from)
3. Formats and prints output in the requested format (text / BED / VCF)

Output formats:
- **text**: human-readable, one line per VNTR
- **BED**: `CHROM START END VNTR_ID GENE MOTIF REF_COPY R1 R2`
- **VCF**: standard VCF 4.2 with INFO fields `END`, `VID`, `RU`, `RC` and FORMAT fields `GT`, `DP`, `SR`, `FR`, `ML`

---

### 3. `advntr/vntr_finder.py` — `VNTRFinder` (core engine)

This is where the algorithmic work happens, one VNTR at a time.

#### Step 1: Read recruitment
Scans the BAM file for reads that overlap the VNTR locus using `pysam`. Also extracts unmapped reads whose k-mers match the VNTR sequence. Each candidate read is scored by a trained Keras classifier (`deep_recruitment.py`). Reads below the minimum score threshold are discarded.

Reads are classified as:
- **Spanning reads** — cover the entire VNTR region (most informative)
- **Flanking reads** — cover one flank and part of the VNTR (used when spanning reads are scarce)

#### Step 2: HMM construction
Calls `hmm_utils.py` to build a profile HMM for this VNTR. The composite HMM has three chained segments:

```
[left-flank suffix HMM] → [repeat unit profile HMM × N] → [right-flank prefix HMM]
```

Each repeat unit profile HMM has one **M** (match), **I** (insert), and **D** (delete) state per base in the repeat pattern, with transition probabilities set from the observed error rate. Multiple repeat unit HMMs are concatenated and the transition matrix is modified to allow variable copy numbers.

The HMM is built via `pomegranate` (see below) and finalized with `model.bake()`.

#### Step 3: Viterbi decoding
Each recruited read is aligned to the composite HMM via the Viterbi algorithm. The optimal state path through the HMM tells you how many repeat-unit states the read traversed — i.e., the read's implied copy number.

#### Step 4: Statistical genotyping
The distribution of per-read copy numbers across all recruited reads is modelled statistically using a binomial model (`scipy.stats.binom`). For diploid organisms, this yields the two most likely allele copy numbers (e.g., `3/5`). For haploid organisms, a single copy number is returned.

---

### 4. `advntr/hmm_utils.py` — HMM construction utilities

Builds pomegranate `HiddenMarkovModel` objects. Key functions:

| Function | Purpose |
|----------|---------|
| `get_prefix_matcher_hmm(pattern)` | HMM for right flanking region |
| `get_suffix_matcher_hmm(pattern)` | HMM for left flanking region |
| `get_constant_number_of_repeats_matcher_hmm(patterns, copies)` | HMM for a fixed number of repeat units |
| `get_variable_number_of_repeats_matcher_hmm(patterns, copies)` | Modifies the transition matrix to allow variable copy numbers |
| `get_read_matcher_model(left, right, patterns, copies)` | Assembles the full composite HMM |
| `build_reference_repeat_finder_hmm(patterns, copies)` | HMM for scanning the reference genome |

Emission probabilities for each state are derived from `profile_hmm.py`, which runs a multiple-sequence alignment (via MUSCLE) across observed repeat instances to compute position-specific substitution rates.

---

### 5. `pomegranate/` — Vendored HMM library

adVNTR ships its own copy of **pomegranate 0.6.1** (a Cython-compiled HMM library), rather than depending on the PyPI package. This is because the modern pomegranate 1.x has a completely different API.

The vendored copy is compiled at import time via `pyximport` (see `pomegranate/__init__.py`). Key files:

| File | Role |
|------|------|
| `hmm.pyx` | `HiddenMarkovModel` class — graph construction, `bake()`, Viterbi, forward/backward |
| `distributions.pyx` | `DiscreteDistribution` and other emission distributions |
| `base.pyx` | `Model`, `State`, `GraphModel` base classes |
| `utils.pyx` | Math utilities (log-sum-exp, etc.) |

The graph topology of the HMM is stored as a `networkx.DiGraph`. The `bake()` method finalizes the graph into a flat array representation for efficient computation — this is where `topological_sort()` is called.

---

### 6. Supporting modules

| File | Role |
|------|------|
| `advntr/models.py` | SQLite I/O for the VNTR model database |
| `advntr/reference_vntr.py` | `ReferenceVNTR` data class — chromosome, position, pattern, flanking regions |
| `advntr/profile_hmm.py` | Builds per-position emission profiles from multiple sequence alignments |
| `advntr/deep_recruitment.py` | Keras-based classifier for read recruitment scoring |
| `advntr/sam_utils.py` | `pysam` wrappers for BAM/SAM access |
| `advntr/coverage_bias.py` | Corrects for coverage non-uniformity near VNTR loci |
| `advntr/vntr_graph.py` | Graph-based analysis of VNTR overlap (uses networkx) |
| `advntr/settings.py` | Global constants (error rates, file paths, thread count) |

---

## Inputs and outputs

### Inputs

| Input | CLI flag | Format | Notes |
|-------|----------|--------|-------|
| Aligned reads | `-a` | BAM / SAM / CRAM | Primary input for Illumina/PacBio/Nanopore |
| VNTR model database | `-m` | SQLite `.db` file | Defaults to bundled database |
| Reference genome | `-r` | FASTA | Only required for CRAM input |
| Target VNTR IDs | `-vid` | Comma-separated integers | Defaults to all VNTRs in the database |
| Average coverage | `-c` | Float | Required only with `--expansion` flag |

### Outputs

| Format | CLI flag | Description |
|--------|----------|-------------|
| Text (default) | `-of text` | VNTR ID → allele copy numbers |
| BED | `-of bed` | Genomic coordinates + copy numbers |
| VCF | `-of vcf` | VCF 4.2 with GT/DP/SR/FR/ML fields |
| Log file | automatic | Written to `<working_dir>/log_<input>.log` |

### Key output fields (VCF FORMAT)

| Tag | Meaning |
|-----|---------|
| `GT` | Genotype (e.g., `3/5` = 3 copies on one allele, 5 on the other) |
| `DP` | Total read depth at the locus |
| `SR` | Number of spanning reads used |
| `FR` | Number of flanking reads used |
| `ML` | Maximum likelihood of the called genotype |

---

## Data flow diagram

```
BAM/SAM/CRAM file
        │
        ▼
  Read recruitment ──────────────── Keras classifier (deep_recruitment.py)
        │                           scores candidate reads
        ▼
  HMM construction ──────────────── profile_hmm.py (per-position emission probs)
        │                           hmm_utils.py (assemble composite HMM)
        │                           pomegranate (HiddenMarkovModel.bake())
        ▼
  Viterbi decoding ──────────────── per read: optimal state path → copy count
        │
        ▼
  Statistical genotyping ─────────── scipy.stats.binom → diploid allele calls
        │
        ▼
  Output (text / BED / VCF)
```

---

## Where errors come from

| Error message | Root cause | File |
|---------------|-----------|------|
| `topological_sort() got an unexpected keyword argument 'nbunch'` | networkx ≥ 2.6 removed `nbunch` parameter | `pomegranate/hmm.pyx:874` |
| `connected_component_subgraphs` AttributeError | networkx ≥ 2.4 removed this function | `advntr/vntr_graph.py:27` |
| `from keras import ...` ImportError | TensorFlow 2.x moved Keras under `tf.keras` | `advntr/vntr_finder.py:8`, `deep_recruitment.py` |
| `Skipping genotyping for this VNTR` | Any exception in `VNTRFinder` is caught and logged; genotyping continues for other VNTRs | `advntr/genome_analyzer.py` |
