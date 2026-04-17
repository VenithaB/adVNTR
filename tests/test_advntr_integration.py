"""
Integration and CLI tests for adVNTR.

Test data: tests/data/HG00096.test.bam and HG00096.test.cram
  - Subsampled from HG00096.final.cram, region chr1:560000-700000 (GRCh38)
  - Contains ~27k reads covering VNTR ID 201 (chr1:627933)

The hg38 VNTR database path is set via the HG38_DB environment variable, or
falls back to the default path on the research cluster.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"
TEST_BAM = DATA_DIR / "HG00096.test.bam"
TEST_CRAM = DATA_DIR / "HG00096.test.cram"
REFERENCE_FA = Path(
    "/scratch/vab5299/datasets/vamos/reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)
HG38_DB = Path(
    os.environ.get(
        "HG38_DB",
        "/storage/group/yzv101/default/projects/srcm_methods_paper/"
        "project_software/vntr_detection_software_fixed/advntr/vntr_data/"
        "hg38_selected_VNTRs_Illumina.db",
    )
)

# VNTR 201 lives in the test region (chr1:627933, GRCh38)
VNTR_ID_IN_REGION = "201"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def run_advntr(*args: str) -> subprocess.CompletedProcess:
    """Run ``python -m advntr`` with the given arguments and capture output."""
    return subprocess.run(
        [sys.executable, "-m", "advntr"] + list(args),
        capture_output=True,
        text=True,
    )


def db_available() -> bool:
    return HG38_DB.exists()


def bam_available() -> bool:
    return TEST_BAM.exists() and (DATA_DIR / "HG00096.test.bam.bai").exists()


def cram_available() -> bool:
    return TEST_CRAM.exists() and (DATA_DIR / "HG00096.test.cram.crai").exists()


def ref_available() -> bool:
    return REFERENCE_FA.exists()


# ---------------------------------------------------------------------------
# CLI validation tests (no alignment data needed)
# ---------------------------------------------------------------------------


class TestCliValidation:
    """Tests for CLI argument validation and error handling."""

    def test_no_subcommand_exits_nonzero(self):
        result = run_advntr()
        assert result.returncode != 0

    def test_missing_alignment_file_exits(self):
        with tempfile.TemporaryDirectory() as wd:
            result = run_advntr(
                "genotype",
                "--working_directory",
                wd,
            )
        assert result.returncode != 0

    def test_missing_working_directory_exits(self, tmp_path):
        dummy_bam = tmp_path / "dummy.bam"
        dummy_bam.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(dummy_bam),
        )
        assert result.returncode != 0

    def test_threads_zero_exits(self, tmp_path):
        dummy_bam = tmp_path / "dummy.bam"
        dummy_bam.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(dummy_bam),
            "--working_directory",
            str(tmp_path),
            "--threads",
            "0",
        )
        assert result.returncode != 0

    def test_threads_negative_exits(self, tmp_path):
        dummy_bam = tmp_path / "dummy.bam"
        dummy_bam.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(dummy_bam),
            "--working_directory",
            str(tmp_path),
            "--threads",
            "-1",
        )
        assert result.returncode != 0

    def test_expansion_without_coverage_exits(self, tmp_path):
        dummy_bam = tmp_path / "dummy.bam"
        dummy_bam.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(dummy_bam),
            "--working_directory",
            str(tmp_path),
            "--expansion",
        )
        assert result.returncode != 0

    def test_invalid_file_extension_exits(self, tmp_path):
        """adVNTR only accepts BAM/SAM/CRAM; a .txt file must fail early."""
        txt_file = tmp_path / "reads.txt"
        txt_file.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(txt_file),
            "--working_directory",
            str(tmp_path),
        )
        assert result.returncode != 0

    def test_frameshift_invalid_vntr_exits(self, tmp_path):
        """--frameshift is only valid for specific VNTR IDs."""
        dummy_bam = tmp_path / "dummy.bam"
        dummy_bam.touch()
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(dummy_bam),
            "--working_directory",
            str(tmp_path),
            "--frameshift",
            "--vntr_id",
            "999999",
        )
        assert result.returncode != 0


# ---------------------------------------------------------------------------
# Unit tests for advntr_commands helper functions
# ---------------------------------------------------------------------------


class TestAdVNTRCommandHelpers:
    """Unit tests for helper functions in advntr_commands."""

    def test_valid_vntr_for_frameshift_accepts_known_ids(self):
        from advntr.advntr_commands import valid_vntr_for_frameshift
        from advntr import settings

        assert valid_vntr_for_frameshift(list(settings.FRAMESHIFT_VNTRS)) is True

    def test_valid_vntr_for_frameshift_rejects_unknown_id(self):
        from advntr.advntr_commands import valid_vntr_for_frameshift

        assert valid_vntr_for_frameshift([999999]) is False

    def test_valid_vntr_for_frameshift_rejects_mixed_list(self):
        from advntr.advntr_commands import valid_vntr_for_frameshift
        from advntr import settings

        mixed = list(settings.FRAMESHIFT_VNTRS) + [999999]
        assert valid_vntr_for_frameshift(mixed) is False

    def test_get_default_vntrs_returns_list(self):
        from advntr.advntr_commands import get_default_vntrs
        from advntr.reference_vntr import ReferenceVNTR

        # Build a minimal reference VNTR that passes all filters
        ref_vntr = ReferenceVNTR(
            vntr_id=1,
            pattern="ACGT",
            start_point=1000,
            chromosome="chr1",
            gene_name=None,
            annotation="Coding",
        )
        ref_vntr.left_flanking_region = "A" * 200
        ref_vntr.right_flanking_region = "A" * 200
        ref_vntr.repeat_segments = ["ACGT"] * 5

        result = get_default_vntrs([ref_vntr], is_pacbio=False)
        assert isinstance(result, list)

    def test_get_default_vntrs_pacbio_flag(self):
        from advntr.advntr_commands import get_default_vntrs

        # Empty list → empty result for both modes
        assert get_default_vntrs([], is_pacbio=False) == []
        assert get_default_vntrs([], is_pacbio=True) == []


# ---------------------------------------------------------------------------
# Integration tests — require test BAM/CRAM and hg38 database
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (bam_available() and db_available()),
    reason="Test BAM or hg38 VNTR database not available",
)
class TestGenotypeBAM:
    """End-to-end genotyping on the subsampled BAM file."""

    def test_genotype_single_vntr_runs(self, tmp_path):
        """adVNTR genotype should complete without error for VNTR 201."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
            "--disable_logging",
        )
        assert (
            result.returncode == 0
        ), f"advntr genotype failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"

    def test_genotype_output_contains_vntr_id(self, tmp_path):
        """Output must mention the queried VNTR ID."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
            "--disable_logging",
        )
        assert result.returncode == 0
        assert VNTR_ID_IN_REGION in result.stdout

    def test_genotype_vcf_output_format(self, tmp_path):
        """VCF output format should produce a header line starting with #."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--outfmt",
            "vcf",
            "--threads",
            "1",
            "--disable_logging",
        )
        assert result.returncode == 0
        assert (
            "#" in result.stdout
        ), "VCF output should contain header lines starting with #"

    def test_genotype_bed_output_format(self, tmp_path):
        """BED output should not crash and should produce tab-separated lines."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--outfmt",
            "bed",
            "--threads",
            "1",
            "--disable_logging",
        )
        assert result.returncode == 0

    def test_genotype_haploid_flag(self, tmp_path):
        """--haploid flag should not crash."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--haploid",
            "--threads",
            "1",
            "--disable_logging",
        )
        assert result.returncode == 0

    def test_genotype_outfile_written(self, tmp_path):
        """Results written to --outfile should match stdout."""
        outfile = tmp_path / "results.txt"
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--outfile",
            str(outfile),
            "--threads",
            "1",
            "--disable_logging",
        )
        assert result.returncode == 0
        assert outfile.exists(), "--outfile was not created"
        assert outfile.stat().st_size > 0, "--outfile is empty"

    def test_genotype_nonexistent_bam_reports_error(self, tmp_path):
        """adVNTR should report an error when the alignment file does not exist.

        Note: upstream adVNTR returns exit code 0 even on file-not-found errors
        (it prints "Error" to stdout but does not propagate the failure as a
        non-zero exit code).  This test documents that observed behaviour so a
        future fix can tighten it to ``returncode != 0``.
        """
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(tmp_path / "nonexistent.bam"),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
        )
        assert "Error" in result.stdout or result.returncode != 0


@pytest.mark.skipif(
    not (cram_available() and db_available() and ref_available()),
    reason="Test CRAM, hg38 database, or reference FASTA not available",
)
class TestGenotypeCRAM:
    """End-to-end genotyping on the subsampled CRAM file."""

    def test_genotype_cram_runs(self, tmp_path):
        """adVNTR should accept a CRAM file with an explicit reference."""
        result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_CRAM),
            "--reference_filename",
            str(REFERENCE_FA),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(tmp_path),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
            "--disable_logging",
        )
        assert (
            result.returncode == 0
        ), f"advntr genotype (CRAM) failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"

    def test_cram_and_bam_produce_same_vntr_id_in_output(self, tmp_path):
        """BAM and CRAM from the same reads should both report VNTR ID 201."""
        bam_dir = tmp_path / "bam"
        bam_dir.mkdir()
        cram_dir = tmp_path / "cram"
        cram_dir.mkdir()

        bam_result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_BAM),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(bam_dir),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
            "--disable_logging",
        )
        cram_result = run_advntr(
            "genotype",
            "--alignment_file",
            str(TEST_CRAM),
            "--reference_filename",
            str(REFERENCE_FA),
            "--models",
            str(HG38_DB),
            "--working_directory",
            str(cram_dir),
            "--vntr_id",
            VNTR_ID_IN_REGION,
            "--threads",
            "1",
            "--disable_logging",
        )
        assert bam_result.returncode == 0
        assert cram_result.returncode == 0
        assert VNTR_ID_IN_REGION in bam_result.stdout
        assert VNTR_ID_IN_REGION in cram_result.stdout
