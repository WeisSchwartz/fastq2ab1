from fastq2ab1.consensus import (
    build_consensus,
    compute_quality,
    peak_locations,
    scaled_channel_traces,
)
from fastq2ab1.pileup import PileupResult


def make_pileup(A, T, G, C):
    n = len(A)
    return PileupResult(
        length=n,
        A=list(A), T=list(T), G=list(G), C=list(C),
        total_depth=[a + t + g + c for a, t, g, c in zip(A, T, G, C)],
        n_reads_total=1, n_reads_aligned=1,
    )


def test_scaling_global_max_equals_1000():
    p = make_pileup(A=[100, 50], T=[0, 50], G=[0, 0], C=[0, 0])
    cons = build_consensus(p, min_coverage=1)
    assert cons.global_max == 100
    assert abs(cons.scale - 10.0) < 1e-9
    G, A, T, C = scaled_channel_traces(p, cons.scale)
    # Plasmidsaurus-style sharp peaks: [0, 0, v, v] per base.
    assert A[0:4] == [0, 0, 1000, 1000]   # base 0: A=100 → 1000 peak
    assert A[4:8] == [0, 0, 500, 500]     # base 1: A=50 → 500 peak
    assert T[4:8] == [0, 0, 500, 500]
    assert G[0:4] == [0, 0, 0, 0]


def test_peak_locations():
    assert peak_locations(4) == [3, 7, 11, 15]


def test_quality_fixed_divisor_20():
    # Plasmidsaurus: q = min(depth // 20, 50). Verified against real AB1.
    assert compute_quality([0, 19, 20, 21, 400, 999, 1000, 5000]) == [
        0, 0, 1, 1, 20, 49, 50, 50,
    ]


def test_quality_independent_of_max_depth():
    """Critical: quality does NOT normalize by max_depth — proves the
    fixed-divisor formula. A sample with max depth 100 and a sample with
    max depth 5000 must give the same quality at the same absolute depth."""
    p1 = make_pileup(A=[40, 100], T=[0, 0], G=[0, 0], C=[0, 0])
    p2 = make_pileup(A=[40, 5000], T=[0, 0], G=[0, 0], C=[0, 0])
    c1 = build_consensus(p1, min_coverage=1)
    c2 = build_consensus(p2, min_coverage=1)
    assert c1.quality[0] == c2.quality[0] == 2  # depth=40 → 40//20 = 2


def test_consensus_calls_argmax_and_n_for_low_coverage():
    p = make_pileup(A=[5, 0, 1], T=[1, 0, 1], G=[1, 0, 1], C=[0, 0, 1])
    cons = build_consensus(p, min_coverage=3, ambiguous="N")
    # position 0: depth 7, argmax A
    # position 1: depth 0, low cov → N
    # position 2: depth 4 but tied A/T/G/C → N
    assert cons.sequence == "ANN"


def test_consensus_iupac_tie():
    p = make_pileup(A=[5, 5], T=[0, 5], G=[0, 0], C=[0, 0])
    cons = build_consensus(p, min_coverage=1, ambiguous="iupac")
    assert cons.sequence[1] == "W"  # A/T → W
