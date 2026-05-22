"""BAM → per-position A/T/G/C count arrays via pysam.

Substitution-level counts only: insertions, deletions, and N base calls are
excluded. For a doubled reference, only the first `reference_length` positions
of the BAM are returned, but reads whose alignment START is in the second copy
are filtered upstream (they would otherwise contribute to no positions in the
returned window anyway).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List

import pysam

log = logging.getLogger("fastq2ab1.pileup")


@dataclass
class PileupResult:
    length: int
    A: List[int]
    T: List[int]
    G: List[int]
    C: List[int]
    total_depth: List[int]
    n_reads_total: int
    n_reads_aligned: int

    def alignment_rate(self) -> float:
        if self.n_reads_total == 0:
            return 0.0
        return self.n_reads_aligned / self.n_reads_total


def compute_pileup(
    bam_path: str,
    contig: str,
    reference_length: int,
    min_mapping_quality: int = 0,
) -> PileupResult:
    """Walk the BAM and tally A/T/G/C at each reference position.

    Handles the doubled-reference trick used for circular plasmids: we walk
    the pileup across the **entire** mapped region (which may be up to 2× the
    original reference length) and collapse positions modulo `reference_length`.
    That way a read that aligned to the second copy contributes its bases to
    the equivalent positions in the first copy, eliminating the dead zone at
    the linearization breakpoint.

    For BAMs whose header contig length already equals `reference_length`
    (i.e. not from a doubled reference), the collapse is a no-op since
    `position % reference_length == position`.
    """
    A = [0] * reference_length
    T = [0] * reference_length
    G = [0] * reference_length
    C = [0] * reference_length

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            stats = bam.get_index_statistics()
            n_aligned = sum(s.mapped for s in stats)
            n_unmapped = sum(s.unmapped for s in stats)
            n_total = n_aligned + n_unmapped
        except Exception:
            n_total = 0
            n_aligned = 0

        # Determine the contig's actual mapped length in this BAM. May be
        # 2*reference_length (doubled, from our aligner) or reference_length
        # (user-supplied BAM aligned to the original linear ref).
        contig_length = None
        for ref in bam.header.to_dict().get("SQ", []):
            if ref["SN"] == contig:
                contig_length = int(ref["LN"])
                break
        if contig_length is None:
            contig_length = reference_length
        # Don't scan past 2× — anything more would alias multiple copies.
        stop = min(contig_length, 2 * reference_length)

        for col in bam.pileup(
            contig=contig,
            start=0,
            stop=stop,
            truncate=True,
            min_mapping_quality=min_mapping_quality,
            ignore_overlaps=False,
            min_base_quality=0,
        ):
            raw_pos = col.reference_pos
            pos = raw_pos % reference_length  # collapse doubled coord
            for read in col.pileups:
                aln = read.alignment
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                if read.is_del or read.is_refskip:
                    continue
                qpos = read.query_position
                if qpos is None:
                    continue
                base = aln.query_sequence[qpos].upper()
                if base == "A":
                    A[pos] += 1
                elif base == "T":
                    T[pos] += 1
                elif base == "G":
                    G[pos] += 1
                elif base == "C":
                    C[pos] += 1
                # N and others are dropped

    total = [A[i] + T[i] + G[i] + C[i] for i in range(reference_length)]
    return PileupResult(
        length=reference_length,
        A=A,
        T=T,
        G=G,
        C=C,
        total_depth=total,
        n_reads_total=n_total,
        n_reads_aligned=n_aligned,
    )


def low_coverage_positions(
    pileup: PileupResult, min_coverage: int
) -> List[int]:
    return [i for i, d in enumerate(pileup.total_depth) if d < min_coverage]
