"""AB1 round-trip: write our binary, read it back with BioPython, assert match."""

from pathlib import Path

import pytest

from fastq2ab1.ab1_writer import write_ab1


def _make_inputs(n=20, global_max_count=100):
    # synthetic per-base counts with a clear winner each position
    seq = "ACGT" * (n // 4)
    seq = seq[:n]
    quality = [min(50, (i * 50) // (n - 1)) for i in range(n)]
    # 4x oversampled flat traces, peak at 4*i+3
    def expand(vals):
        out = []
        for v in vals:
            out.extend([v, v, v, v])
        return out

    g_counts = [global_max_count if b == "G" else 0 for b in seq]
    a_counts = [global_max_count if b == "A" else 0 for b in seq]
    t_counts = [global_max_count if b == "T" else 0 for b in seq]
    c_counts = [global_max_count if b == "C" else 0 for b in seq]

    # Apply scale so global max → 1000.
    scale = 1000.0 / global_max_count
    g = [int(round(v * scale)) for v in g_counts]
    a = [int(round(v * scale)) for v in a_counts]
    t = [int(round(v * scale)) for v in t_counts]
    c = [int(round(v * scale)) for v in c_counts]

    ploc = [4 * i + 3 for i in range(n)]
    return seq, quality, expand(g), expand(a), expand(t), expand(c), ploc


def test_roundtrip_with_biopython(tmp_path: Path):
    Bio = pytest.importorskip("Bio")
    from Bio import SeqIO

    seq, qual, tG, tA, tT, tC, ploc = _make_inputs(n=24, global_max_count=200)
    out = tmp_path / "test.ab1"
    write_ab1(
        output_path=str(out),
        sequence=seq,
        quality=qual,
        trace_G=tG, trace_A=tA, trace_T=tT, trace_C=tC,
        peak_locations=ploc,
        sample_name="sample_test",
        comment="fastq2ab1 unit test",
    )

    rec = SeqIO.read(str(out), "abi")
    raw = rec.annotations["abif_raw"]

    def _to_str(v):
        return v.decode() if isinstance(v, (bytes, bytearray)) else v

    # Sequence and quality
    assert str(rec.seq) == seq
    assert _to_str(raw["PBAS1"]) == seq
    assert _to_str(raw["PBAS2"]) == seq
    assert list(raw["PCON1"]) == qual
    assert list(raw["PCON2"]) == qual

    # Traces
    assert list(raw["DATA9"]) == tG
    assert list(raw["DATA10"]) == tA
    assert list(raw["DATA11"]) == tT
    assert list(raw["DATA12"]) == tC

    # Peak locations and base order
    assert list(raw["PLOC1"]) == ploc
    assert _to_str(raw["FWO_1"]) == "GATC"


def test_roundtrip_large_reference_ploc_overflow(tmp_path: Path):
    """Reference > ~8 kb has PLOC values that overflow int16. Writer must
    auto-promote PLOC1 to int32 so the file remains valid and BioPython can
    read the values back exactly."""
    pytest.importorskip("Bio")
    from Bio import SeqIO

    n = 12000  # peak[-1] = 4*11999 + 3 = 47999, well over int16 max 32767
    seq = ("ACGT" * (n // 4 + 1))[:n]
    quality = [40] * n
    flat = [100] * (4 * n)  # 4x oversampled, same shape on all channels
    ploc = [4 * i + 3 for i in range(n)]

    out = tmp_path / "big.ab1"
    write_ab1(
        output_path=str(out),
        sequence=seq, quality=quality,
        trace_G=flat, trace_A=flat, trace_T=flat, trace_C=flat,
        peak_locations=ploc,
        sample_name="big_ref", comment="overflow regression",
    )

    rec = SeqIO.read(str(out), "abi")
    raw = rec.annotations["abif_raw"]
    assert list(raw["PLOC1"]) == ploc
    assert max(raw["PLOC1"]) > 32767
