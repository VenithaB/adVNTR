# adVNTR Log File Format

This document describes every message type written to the adVNTR log file during
`advntr genotype` on Illumina short-read data (BAM/CRAM input).

The log is written to `<working_directory>/log_<sample>.log`.

---

## Log levels

adVNTR uses standard Python logging levels. The log contains `INFO` and `DEBUG`
messages only — no `WARNING` or `ERROR` messages are expected during a clean run.

---

## Stage 1 — Startup

| Level | Message | Meaning |
|-------|---------|---------|
| INFO | `adVNTR 1.5.0` | Tool version |
| INFO | `Reference VNTR DB: <path>` | Path to the `.db` model file |
| INFO | `Running adVNTR for <N> VNTRs` | Total number of loci to genotype |

---

## Stage 2 — Unmapped read extraction

| Level | Message | Meaning |
|-------|---------|---------|
| INFO | `extract_unmapped_reads_to_fasta_file executed in <t>s` | Time to extract all unmapped reads from the CRAM/BAM to a FASTA file. This is the slowest upfront step (~15 min for a 30x CRAM). |

---

## Stage 3 — Keyword filtering

| Level | Message | Meaning |
|-------|---------|---------|
| INFO | `get_filtered_read_ids executed in <t>s` | Time to run keyword-based pre-screening across all loci |
| INFO | `get_vntr_filtered_reads_map executed in <t>s` | Total time to build the per-VNTR read candidate map |
| INFO | `filtering selected <N> reads for <VNTR_ID>` | Number of unmapped reads passing keyword filter for this locus. Most loci return 0 (normal — most VNTRs are covered by mapped reads). Capped at 2000. |

---

## Stage 4 — HMM genotyping (per locus)

Each VNTR is processed serially. For Illumina data the `-t/--threads` flag has no
effect — parallelism is only implemented for PacBio reads.

### 4a. Locus header

| Level | Message | Meaning |
|-------|---------|---------|
| DEBUG | `finding repeat count from alignment file for <VNTR_ID>` | Start of processing for this locus |
| DEBUG | `accuracy filter is False` | Using both spanning and flanking reads for genotyping. Set `--accuracy_filter` to restrict to spanning reads only. |
| DEBUG | `left and right min flanking size: 5 5` | Minimum bp of flanking sequence required on each side for a read to be counted |
| INFO | `Using read length <N>` | Read length detected from the alignment file header |

### 4b. HMM construction

| Level | Message | Meaning |
|-------|---------|---------|
| INFO | `get_suffix_matcher_hmm executed in <t>s` | Built the HMM for the left flanking region |
| INFO | `build_profile_hmm_pseudocounts_for_alignment executed in <t>s` | Computed pseudocounts for the repeat unit profile HMM |
| INFO | `build_profile_hmm_for_repeats executed in <t>s` | Built the repeat unit profile HMM |
| INFO | `get_constant_number_of_repeats_matcher_hmm executed in <t>s` | Built the fixed-copy-number HMM component |
| INFO | `get_variable_number_of_repeats_matcher_hmm executed in <t>s` | Built the variable-copy-number HMM |
| INFO | `get_prefix_matcher_hmm executed in <t>s` | Built the HMM for the right flanking region |
| INFO | `get_read_matcher_model executed in <t>s` | Assembled the combined read matcher |
| INFO | `build_vntr_matcher_hmm executed in <t>s` | Final combined HMM (suffix + repeat + prefix) is ready |
| INFO | `select_illumina_reads executed in <t>s` | Ran candidate reads through the HMM to select those that align |

### 4c. Per-read evaluation

For every read evaluated against the HMM:

| Level | Message | Meaning |
|-------|---------|---------|
| DEBUG | `logp of read: <value>` | Log-probability of this read under the HMM model |
| DEBUG | `repeats: <N>` | Number of repeat units this read spans |
| DEBUG | `left flanking size: <N>` / `right flanking size: <N>` | Flanking coverage on each side (bp) |
| DEBUG | `<sequence>` | The actual flanking sequence string for this read |
| DEBUG | `get_flanking_regions_matching_rate: accuracy filter T <t> F <f>` | Alignment quality scores for flanking regions under strict (T) and permissive (F) thresholds |

Each read is then classified:

| Level | Message | Meaning |
|-------|---------|---------|
| DEBUG | `spanning read <readID> sourced from MAPPED/UNMAPPED visited states: [...]` | Read fully crosses the repeat region. The state list is the Viterbi path through the HMM (e.g. `M1_0`, `unit_end_0`, `unit_start_1`, ...). |
| DEBUG | `flanking read <readID> sourced from MAPPED/UNMAPPED visited states: [...]` | Read partially overlaps the repeat — begins or ends within it. Also includes full Viterbi path. |
| DEBUG | `Rejected Aligned Read` | Mapped read that failed HMM scoring |
| DEBUG | `Rejecting duplicated read` | PCR duplicate, skipped |
| DEBUG | `Rejecting read for poor mapping quality` | MAPQ below threshold |
| DEBUG | `Rejecting read for having so many low quality base pairs` | Too many low base-quality positions |
| DEBUG | `Rejecting read for having long run of low quality base pairs` | Consecutive run of low-quality bases |

### 4d. Genotype result

| Level | Message | Meaning |
|-------|---------|---------|
| DEBUG | `vntr base pairs in unmapped reads: <N>` | Total VNTR-overlapping bp contributed by unmapped reads |
| DEBUG | `vntr base pairs in mapped reads: <N>` | Total VNTR-overlapping bp contributed by mapped reads |
| DEBUG | `number of recruited reads: <N>` | Total reads used for genotyping at this locus |
| INFO | `covered repeats: [N, N, ...]` | Repeat unit counts observed from spanning reads |
| INFO | `flanking repeats: [N, N, ...]` | Repeat unit counts inferred from flanking reads |
| INFO | `Maximum probability for genotyping: <p>` | Confidence of the genotype call. Most well-covered loci produce values >0.999. |
| INFO | `RU count lower bounds: <A>/<B>` | Final diploid genotype — repeat units on each haplotype |
| INFO | `Average coverage is not set` | No `--coverage` flag provided; long-expansion detection disabled |
| INFO | `find_repeat_count_from_alignment_file executed in <t>s` | Total wall time for this locus |

---

## Notes

- Read name prefixes (e.g. `A00132`, `SNA00132`) reflect different sequencing runs
  merged into the same CRAM — both are used normally.
- The `visited states` lists in spanning/flanking read entries show the Viterbi
  decoding path: `M<pos>_<unit>` = match state at position in repeat unit,
  `I<pos>_<unit>` = insertion, `unit_end/unit_start` = repeat unit boundaries,
  `prefix_*` / `suffix_*` = flanking HMM states.
- `RU count lower bounds` is a lower bound because flanking reads can only confirm
  that *at least* N repeat units are present. The true count may be higher for
  very long expansions without `--expansion` mode.
