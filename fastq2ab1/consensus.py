"""Consensus and quality computation from per-position counts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from .pileup import PileupResult

# Empirically reverse-engineered from a real Plasmidsaurus AB1:
# PCON values follow `q = min(depth // QUALITY_DIVISOR, 50)` with a fixed
# divisor of 20 (not `max_depth / 50` as the original spec speculated).
# At depth=20 → q=1, depth=1000 → q=50 (saturated). Verified against
# G399MF_1_LP_CAG_1.ab1: depth/PCON ratio was 20.000 ± rounding across all
# 7,647 non-saturated positions.
QUALITY_DIVISOR = 20
QUALITY_MAX = 50

# IUPAC ambiguity codes for ties.
_IUPAC = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}


@dataclass
class ConsensusResult:
    sequence: str
    quality: List[int]   # PCON values, 0–50
    global_max: int      # max count across all four channels & positions
    global_max_depth: int  # max total_depth (informational only)
    scale: float         # 1000 / global_max — applied to channel counts


def _call_base(
    a: int, t: int, g: int, c: int, depth: int, min_coverage: int,
    ambiguous: str
) -> str:
    if depth < min_coverage or depth == 0:
        return "N"
    counts = {"A": a, "T": t, "G": g, "C": c}
    top = max(counts.values())
    winners = [b for b, v in counts.items() if v == top]
    if len(winners) == 1:
        return winners[0]
    if ambiguous == "iupac":
        return _IUPAC[frozenset(winners)]
    return "N"


def build_consensus(
    pileup: PileupResult,
    min_coverage: int = 10,
    ambiguous: str = "N",
) -> ConsensusResult:
    if ambiguous not in ("N", "iupac"):
        raise ValueError(f"ambiguous must be 'N' or 'iupac', got {ambiguous!r}")

    n = pileup.length
    seq_chars: List[str] = []
    for i in range(n):
        seq_chars.append(
            _call_base(
                pileup.A[i],
                pileup.T[i],
                pileup.G[i],
                pileup.C[i],
                pileup.total_depth[i],
                min_coverage,
                ambiguous,
            )
        )
    sequence = "".join(seq_chars)

    global_max = 0
    for arr in (pileup.A, pileup.T, pileup.G, pileup.C):
        for v in arr:
            if v > global_max:
                global_max = v

    if global_max == 0:
        raise ValueError(
            "No aligned reads contribute to any reference position "
            "(global_max == 0). Cannot generate AB1."
        )

    scale = 1000.0 / global_max
    global_max_depth = max(pileup.total_depth) if pileup.total_depth else 0
    if global_max_depth == 0:
        raise ValueError("global_max_depth == 0; no coverage anywhere.")

    quality = compute_quality(pileup.total_depth)

    return ConsensusResult(
        sequence=sequence,
        quality=quality,
        global_max=global_max,
        global_max_depth=global_max_depth,
        scale=scale,
    )


def compute_quality(total_depth: Sequence[int]) -> List[int]:
    """Plasmidsaurus-style PCON quality from coverage depth.

    `q = min(depth // 20, 50)` — saturates at depth ≥ 1000. A fixed divisor
    (not `max_depth / 50`) matches what real Plasmidsaurus AB1 files do.
    """
    out: List[int] = []
    for d in total_depth:
        if d <= 0:
            out.append(0)
            continue
        q = int(d) // QUALITY_DIVISOR
        if q > QUALITY_MAX:
            q = QUALITY_MAX
        out.append(q)
    return out


def scaled_channel_traces(
    pileup: PileupResult, scale: float
) -> Tuple[List[int], List[int], List[int], List[int]]:
    """Return (G, A, T, C) scaled int16-ready trace arrays, 4× oversampled.

    Order matches DATA9/10/11/12 in the AB1 spec: G, A, T, C.

    **Trace shape (empirically matched to Plasmidsaurus):** each base
    occupies 4 samples but only the last two carry the value — pattern is
    `[0, 0, v, v]`. The leading two zeros produce the sharp-peak look in
    SnapGene. (The original spec described `[v, v, v, v]`, which produces
    flat-topped plateaus and looks visually inferior.)
    """

    def expand(arr):
        out = []
        for v in arr:
            s = int(round(v * scale))
            if s < 0:
                s = 0
            if s > 32767:
                s = 32767
            out.extend([0, 0, s, s])
        return out

    return expand(pileup.G), expand(pileup.A), expand(pileup.T), expand(pileup.C)


def peak_locations(n_bases: int) -> List[int]:
    return [4 * i + 3 for i in range(n_bases)]
