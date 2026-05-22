"""Regression test for the doubled-reference circular-collapse logic.

When minimap2 aligns against a 2× reference, some reads land in copy 2.
Their alignment_start is ≥ reference_length, but they represent biologically
real coverage at the equivalent position in copy 1. The pileup module must
collapse copy-2 positions back into [0, reference_length) so the
linearization breakpoint has no dead zone.

This test bypasses minimap2 entirely: we construct a BAM whose header
declares a doubled-length contig, write some reads that align ONLY in copy 2,
and verify that compute_pileup returns counts at the equivalent positions
in [0, reference_length).
"""
from pathlib import Path

import pytest


def test_pileup_collapses_doubled_reference(tmp_path: Path):
    pysam = pytest.importorskip("pysam")

    ref_seq = "ACGTACGTAC"  # 10 bp original reference
    ref_len = len(ref_seq)
    doubled_len = 2 * ref_len  # what minimap2 saw

    bam_path = tmp_path / "circ.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "plasmid", "LN": doubled_len}],  # doubled-length contig
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # 3 reads aligned entirely in copy 2 (positions 10..20 in doubled coords).
        for i in range(3):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = f"copy2_only_{i}"
            a.query_sequence = ref_seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = ref_len  # entirely in copy 2
            a.mapping_quality = 60
            a.cigarstring = f"{ref_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * ref_len)
            bam.write(a)
        # 2 reads aligned in copy 1 to position 0.
        for i in range(2):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = f"copy1_{i}"
            a.query_sequence = ref_seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigarstring = f"{ref_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * ref_len)
            bam.write(a)

    pysam.sort("-o", str(bam_path) + ".sorted", str(bam_path))
    Path(str(bam_path) + ".sorted").replace(bam_path)
    pysam.index(str(bam_path))

    from fastq2ab1.pileup import compute_pileup

    pile = compute_pileup(
        bam_path=str(bam_path),
        contig="plasmid",
        reference_length=ref_len,
    )

    # All 5 reads should contribute to every position in [0, ref_len).
    # If copy-2 reads were dropped, depth would be 2 everywhere.
    assert pile.total_depth == [5] * ref_len

    # And the per-base counts should match the reference at every position.
    for i, base in enumerate(ref_seq):
        bucket = {"A": pile.A, "T": pile.T, "G": pile.G, "C": pile.C}[base]
        assert bucket[i] == 5, (i, base, bucket[i])
