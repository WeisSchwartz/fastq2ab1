"""End-to-end integration of pileup → consensus → AB1 writer.

We don't go through `pipeline.run()` here (which would require minimap2). We
build a synthetic BAM in pure Python via pysam, run the internal stack, and
verify the resulting AB1 has the correct trace shape, sequence, and quality.
"""

from pathlib import Path

import pytest


def _make_reference(path: Path) -> str:
    ref = "ACGTACGTACGTACGTACGT"  # 20 bp
    with open(path, "w") as fh:
        fh.write(">plasmid\n")
        fh.write(ref + "\n")
    return ref


def _make_bam(bam_path: Path, ref_name: str, ref_seq: str, n_reads: int = 20) -> None:
    pysam = pytest.importorskip("pysam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": ref_name, "LN": len(ref_seq)}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(n_reads):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = f"r{i}"
            a.query_sequence = ref_seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigarstring = f"{len(ref_seq)}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * len(ref_seq))
            bam.write(a)
    pysam.sort("-o", str(bam_path) + ".sorted", str(bam_path))
    Path(str(bam_path) + ".sorted").replace(bam_path)
    pysam.index(str(bam_path))


def test_internal_stack_pileup_to_ab1(tmp_path: Path):
    """Run the internal stack: synthetic BAM → pileup → consensus → AB1 writer
    → BioPython read-back. Verifies trace shape, consensus, and quality."""
    pytest.importorskip("pysam")
    pytest.importorskip("Bio")
    from Bio import SeqIO

    from fastq2ab1.ab1_writer import write_ab1
    from fastq2ab1.consensus import (
        build_consensus,
        peak_locations,
        scaled_channel_traces,
    )
    from fastq2ab1.pileup import compute_pileup

    ref_path = tmp_path / "ref.fasta"
    ref_seq = _make_reference(ref_path)
    bam_path = tmp_path / "reads.bam"
    _make_bam(bam_path, "plasmid", ref_seq, n_reads=25)

    pile = compute_pileup(
        bam_path=str(bam_path),
        contig="plasmid",
        reference_length=len(ref_seq),
    )
    assert pile.n_reads_aligned == 25

    cons = build_consensus(pile, min_coverage=5)
    assert cons.sequence == ref_seq
    assert cons.global_max == 25
    # depth=25, divisor=20 → q=1 at every position
    assert all(q == 1 for q in cons.quality)

    G, A, T, C = scaled_channel_traces(pile, cons.scale)
    ploc = peak_locations(len(ref_seq))

    out = tmp_path / "out.ab1"
    write_ab1(
        output_path=str(out),
        sequence=cons.sequence, quality=cons.quality,
        trace_G=G, trace_A=A, trace_T=T, trace_C=C,
        peak_locations=ploc,
        sample_name="unit_test",
    )

    rec = SeqIO.read(str(out), "abi")
    raw = rec.annotations["abif_raw"]
    base_to_data = {
        "G": list(raw["DATA9"]),
        "A": list(raw["DATA10"]),
        "T": list(raw["DATA11"]),
        "C": list(raw["DATA12"]),
    }
    for i, base in enumerate(ref_seq):
        # Plasmidsaurus-style sharp peaks: [0, 0, v, v]
        assert base_to_data[base][4 * i] == 0
        assert base_to_data[base][4 * i + 1] == 0
        assert base_to_data[base][4 * i + 2] == 1000
        assert base_to_data[base][4 * i + 3] == 1000
        for other in "ACGT":
            if other != base:
                for k in range(4):
                    assert base_to_data[other][4 * i + k] == 0

    assert str(rec.seq) == ref_seq
