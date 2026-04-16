# adVNTR Modernization Fork

## Project Overview
This is a modernized fork of [adVNTR](https://github.com/mehrdadbakhtiari/adVNTR), a tool for identifying and genotyping Variable Number Tandem Repeats (VNTRs) from short-read sequencing data using Hidden Markov Models (HMMs).

## Refactoring Goals
- Fix dependency incompatibilities to allow adVNTR to run without errors
- **Do NOT** replace the core algorithm or HMM-based genotyping approach
- Update API usage for packages that have changed their interfaces (networkx, pomegranate, etc.)
- Ensure compatibility with modern Python (3.8+) and current package versions

## Known Issues to Fix

### 1. networkx / pomegranate incompatibility
**Error:** `topological_sort() got an unexpected keyword argument 'nbunch'`
- Root cause: adVNTR vendors its own copy of pomegranate (see `pomegranate/` below); that vendored code calls `networkx.topological_sort(graph, nbunch=...)` which was removed in networkx ≥ 2.6
- Fix target: `pomegranate/hmm.pyx` — update the `topological_sort` call to the modern networkx API

### 2. Other potential incompatibilities to audit
- Check all `import` statements against current package APIs
- Review `setup.py` / `requirements.txt` for version constraints
- Check any other deprecated function calls in networkx, scipy, numpy, etc.

## Architecture Notes
- HMM-based genotyping lives primarily in files that call pomegranate
- The `topological_sort` error occurs during `find_repeat_count_from_alignment_file` for Illumina data
- VNTR ID 78547 was the specific failing case reported
- adVNTR ships a **vendored, Cython-compiled fork of pomegranate** in `pomegranate/` rather than depending on the PyPI package — fixes must go there, not into an installed package

## Development Constraints
- Preserve the original algorithmic approach (HMM-based genotyping, profile HMMs, etc.)
- Only update API calls / imports / dependencies — not the logic
- Do not add new features or change output formats

## Conventions
- Tests must pass before committing
- Use snake_case for functions and variables
- All public functions need docstrings and type annotations
- Document every change made, and especially every bug fixed or package updated

## Commit Workflow
- **Before refactoring**: commit the current state so there is a clean baseline to diff against
- **After making changes**: show the user a `git diff` (or equivalent) and seek approval before committing
- **Once approved**: commit with a descriptive message documenting what was changed and why
- **Prompt the user to commit** at these moments even if not actively refactoring:
  1. After tests pass
  2. After `make check` passes
  3. Before switching tasks
  4. At the end of a session

## Validation
- After any code change, run `make check` to validate. All checks must pass before considering work complete.

## Validation Feedback Loops
- Use black, flake8, mypy and pytest, in the same order while refactoring the code.
- Format first, then lint, then type-check, then test — in that order.
- Clean up the code first, then change structure.
- Aim for high test coverage.
- flake8 is run with `--max-line-length 100` to accommodate existing long strings in CLI help text and to avoid conflicts with black's 88-char default.

## Project Structure

```
adVNTR/
├── .gitignore
├── .travis.yml
├── LICENSE
├── Makefile
├── README.md
├── build_advntr_filtering.sh
├── requirements-linux.txt
├── requirements.txt
├── setup.py
│
├── advntr/                         # Main Python package
│   ├── __init__.py
│   ├── __main__.py                 # Entry point
│   ├── acgt_filter.py              # Nucleotide filtering
│   ├── advntr_commands.py          # CLI command dispatch
│   ├── coverage_bias.py            # Coverage bias correction
│   ├── deep_recruitment.py         # Deep learning read recruitment
│   ├── distance.py                 # Sequence distance utilities
│   ├── genome_analyzer.py          # Genome-level analysis
│   ├── hierarchical_clustering.py  # Clustering for VNTR patterns
│   ├── hmm_utils.py                # HMM construction/utilities
│   ├── models.py                   # Data models
│   ├── pacbio_haplotyper.py        # PacBio long-read haplotyping
│   ├── pairwise_aln_generator.py   # Pairwise alignment generation
│   ├── pattern_clustering.py       # VNTR pattern clustering
│   ├── plot.py                     # Visualization
│   ├── profile_hmm.py              # Profile HMM construction
│   ├── profiler.py                 # Performance profiling
│   ├── reference_editor.py         # Reference genome editing
│   ├── reference_vntr.py           # Reference VNTR data structures
│   ├── sam_utils.py                # SAM/BAM file utilities
│   ├── settings.py                 # Global settings
│   ├── utils.py                    # General utilities
│   ├── vntr_annotation.py          # VNTR annotation
│   ├── vntr_finder.py              # Core VNTR finding logic
│   └── vntr_graph.py               # Graph-based VNTR analysis
│
├── pomegranate/                    # Vendored Cython fork of pomegranate (HMM library)
│   ├── __init__.py
│   ├── base.pxd / base.pyx         # Base classes
│   ├── distributions.pxd / .pyx    # Probability distributions
│   ├── hmm.pyx                     # HMM implementation — contains networkx calls to fix
│   └── utils.pxd / utils.pyx       # Utility functions
│
├── filtering/                      # C++ filtering component
│   ├── README.md
│   └── main.cc
│
├── docs/                           # Sphinx documentation
│   ├── conf.py
│   ├── index.rst
│   ├── installation.rst
│   ├── quickstart.rst
│   ├── tutorial.rst
│   ├── faq.rst
│   └── publication.rst
│
└── tests/                          # Test suite
    ├── __init__.py
    ├── test_frameshift_identification.py
    ├── test_genotyping.py
    ├── test_hmm_utils.py
    └── data/
        └── hmm_utils.json
```

## Changelog
<!-- Document every bug fix and package update here as refactoring progresses -->

### 2026-04-16

#### Bug fixes
- **`pomegranate/hmm.pyx` line 874**: Removed deprecated `nbunch` keyword from
  `networkx.topological_sort()` call (removed in networkx ≥ 2.6). New call:
  `list(networkx.topological_sort(silent_subgraph))`. The `nbunch` argument was
  redundant because `silent_subgraph` is already the subgraph of silent states.
- **`advntr/vntr_graph.py` line 23**: Replaced removed
  `nx.connected_component_subgraphs(G)` (removed in networkx 2.4) with the
  modern equivalent `(G.subgraph(c) for c in nx.connected_components(G))`.
- **`advntr/vntr_graph.py` lines 14–17**: Replaced removed
  `from networkx import graphviz_layout` import (moved in networkx ≥ 2.x) with
  `nx.drawing.nx_agraph.graphviz_layout` inside a try/except that falls back to
  spring layout when PyGraphviz is unavailable.

#### Dependency updates
- **`requirements.txt`**: `networkx==1.11` → `networkx>=2.6`; `biopython==1.76`
  → `biopython>=1.76`; removed `enum34` (built-in since Python 3.4).
- **`setup.py`**: `networkx==1.11` → `networkx>=2.6`;
  `biopython` → `biopython>=1.76` in `install_requires`.

#### Test data
- Added `tests/data/HG00096.test.bam` + `.bai` index — subsampled from
  `HG00096.final.cram` (GRCh38, chr1:560000–700000, ~27k reads, 1.8 MB).
- Added `tests/data/HG00096.test.cram` + `.crai` index — same region in CRAM
  format (725 KB). Contains VNTR ID 201 (chr1:627933) at full coverage.

#### New tests (`tests/test_advntr_integration.py`)
- **`TestCliValidation`** (8 tests): CLI argument validation — missing alignment
  file, missing working directory, threads ≤ 0, `--expansion` without
  `--coverage`, unsupported file extension, `--frameshift` with invalid VNTR ID.
- **`TestAdVNTRCommandHelpers`** (5 tests): Unit tests for
  `valid_vntr_for_frameshift` and `get_default_vntrs`.
- **`TestGenotypeBAM`** (7 tests): End-to-end genotyping on the subsampled BAM
  using the hg38 VNTR database — text/VCF/BED output formats, `--haploid` flag,
  `--outfile`, and missing-file error reporting.
- **`TestGenotypeCRAM`** (2 tests): CRAM genotyping with explicit reference FASTA;
  BAM/CRAM output agreement. Skipped when reference FASTA is unavailable.
