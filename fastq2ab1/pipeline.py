"""End-to-end pipeline orchestration.

`pipeline.run(...)` is callable as a library (no CLI required), and is the
function that a future web wrapper will call.
"""

from __future__ import annotations

import datetime as _dt
import logging
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import __version__
from .ab1_writer import write_ab1
from .aligner import align_fastq, load_single_reference
from .consensus import (
    build_consensus,
    peak_locations,
    scaled_channel_traces,
)
from .pileup import compute_pileup, low_coverage_positions

log = logging.getLogger("fastq2ab1.pipeline")


@dataclass
class RunStats:
    output_path: Path
    reference_length: int
    n_reads_total: int
    n_reads_aligned: int
    alignment_rate: float
    mean_coverage: float
    n_low_coverage: int
    global_max: int
    global_max_depth: int


def run(
    reference: str | os.PathLike,
    output: str | os.PathLike,
    fastq: str | os.PathLike,
    sample_name: Optional[str] = None,
    min_coverage: int = 10,
    threads: int = 4,
    ambiguous_bases: str = "N",
    write_consensus: str | os.PathLike | None = None,
    strict: bool = False,
    keep_temp: bool = False,
) -> RunStats:
    """End-to-end FASTQ → AB1 conversion.

    The pipeline always aligns FASTQ reads internally with minimap2 against a
    2×-doubled reference so circular plasmids have no coverage taper at the
    linearization breakpoint. Pre-aligned BAM input is intentionally not
    supported: external BAMs aligned to the original linear reference produce
    a bell-curve coverage profile that cannot be recovered post-hoc (the
    soft-clipped bases are simply absent from the BAM).
    """
    reference = Path(reference)
    output = Path(output)
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")
    fastq_path = Path(fastq)
    if not fastq_path.is_file():
        raise FileNotFoundError(f"FASTQ not found: {fastq_path}")

    ref_record = load_single_reference(reference)
    ref_length = len(ref_record.seq)
    log.info("Reference: %s (%d bp)", ref_record.id, ref_length)

    tmp_path = Path(tempfile.mkdtemp(prefix="fastq2ab1_"))
    try:
        bam_path = tmp_path / "aligned.bam"
        log.info("Aligning FASTQ → BAM (threads=%d)…", threads)
        align_fastq(
            fastq=fastq_path,
            reference_fasta=reference,
            output_bam=bam_path,
            threads=threads,
            workdir=tmp_path,
        )
        contig = ref_record.id

        log.info("Computing pileup on contig '%s'…", contig)
        pile = compute_pileup(
            bam_path=str(bam_path),
            contig=contig,
            reference_length=ref_length,
        )

        if pile.n_reads_aligned == 0:
            raise RuntimeError("No aligned reads — aborting.")

        aln_rate = pile.alignment_rate()
        if pile.n_reads_total > 0 and aln_rate < 0.5:
            msg = f"Alignment rate is {aln_rate:.1%} (< 50%)."
            if strict:
                raise RuntimeError(msg + " --strict set; aborting.")
            log.warning(msg)

        cons = build_consensus(
            pile, min_coverage=min_coverage, ambiguous=ambiguous_bases
        )

        low = low_coverage_positions(pile, min_coverage)
        if low:
            log.warning(
                "%d/%d positions below min_coverage=%d (first few: %s)",
                len(low), ref_length, min_coverage, low[:10],
            )

        traces_G, traces_A, traces_T, traces_C = scaled_channel_traces(
            pile, cons.scale
        )
        ploc = peak_locations(ref_length)

        if sample_name is None:
            sample_name = fastq_path.stem

        mean_cov = (
            sum(pile.total_depth) / ref_length if ref_length else 0.0
        )
        comment = (
            f"fastq2ab1 v{__version__}; "
            f"reference={reference.name}; "
            f"length={ref_length}; "
            f"mean_coverage={mean_cov:.1f}; "
            f"aligned_reads={pile.n_reads_aligned}; "
            f"date={_dt.date.today().isoformat()}"
        )

        log.info(
            "Writing AB1 (global_max=%d, max_depth=%d, mean_cov=%.1f)",
            cons.global_max, cons.global_max_depth, mean_cov,
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        write_ab1(
            output_path=str(output),
            sequence=cons.sequence,
            quality=cons.quality,
            trace_G=traces_G,
            trace_A=traces_A,
            trace_T=traces_T,
            trace_C=traces_C,
            peak_locations=ploc,
            sample_name=sample_name,
            comment=comment,
        )

        # --- optional consensus FASTA ---
        if write_consensus is not None:
            out_fa = Path(write_consensus)
            rec = SeqRecord(
                Seq(cons.sequence),
                id=ref_record.id + "_consensus",
                description=f"fastq2ab1 consensus (mean_cov={mean_cov:.1f})",
            )
            with open(out_fa, "w") as fh:
                SeqIO.write([rec], fh, "fasta")
            log.info("Wrote consensus FASTA: %s", out_fa)

        return RunStats(
            output_path=output,
            reference_length=ref_length,
            n_reads_total=pile.n_reads_total,
            n_reads_aligned=pile.n_reads_aligned,
            alignment_rate=aln_rate,
            mean_coverage=mean_cov,
            n_low_coverage=len(low),
            global_max=cons.global_max,
            global_max_depth=cons.global_max_depth,
        )
    finally:
        if keep_temp:
            log.info("Keeping temp dir for debugging: %s", tmp_path)
        else:
            import shutil
            shutil.rmtree(tmp_path, ignore_errors=True)
