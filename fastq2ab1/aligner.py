"""minimap2 wrapper with circular-reference handling.

For a single-sequence circular plasmid reference, we duplicate the reference
(2× length) before alignment so reads spanning the linearization breakpoint
can align end-to-end. The pileup stage then collapses coordinates modulo the
original reference length, so reads landing in either copy contribute to the
correct position in copy 1 — eliminating the coverage taper at the breakpoint.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import require_binary

log = logging.getLogger("fastq2ab1.aligner")


def load_single_reference(reference_fasta: Path) -> SeqRecord:
    """Parse a single-record FASTA and return its SeqRecord.

    Shared between the aligner (needs the sequence to write the doubled
    reference) and the pipeline (needs the length/id for downstream tags).
    """
    records = list(SeqIO.parse(str(reference_fasta), "fasta"))
    if len(records) == 0:
        raise ValueError(f"Reference FASTA is empty: {reference_fasta}")
    if len(records) > 1:
        raise ValueError(
            "Reference must be a single sequence. Use --contig to specify."
        )
    return records[0]


def _write_doubled_reference(record: SeqRecord, out_path: Path) -> None:
    """Write a 2×-concatenated FASTA to `out_path` for circular alignment."""
    doubled = SeqRecord(
        Seq(str(record.seq) + str(record.seq)),
        id=record.id,
        description="doubled for circular alignment",
    )
    with open(out_path, "w") as fh:
        SeqIO.write([doubled], fh, "fasta")


def align_fastq(
    fastq: Path,
    reference_fasta: Path,
    output_bam: Path,
    threads: int = 4,
    workdir: Path | None = None,
) -> None:
    """Align FASTQ → sorted, indexed BAM against a 2×-doubled reference.

    Side-effect: writes `output_bam` and `output_bam.bai`. The doubled
    reference is written to `workdir/ref_doubled.fasta`.
    """
    minimap2 = require_binary("minimap2")
    samtools = require_binary("samtools")

    workdir = workdir or output_bam.parent
    workdir.mkdir(parents=True, exist_ok=True)

    record = load_single_reference(reference_fasta)
    doubled_path = workdir / "ref_doubled.fasta"
    _write_doubled_reference(record, doubled_path)

    log.info("Indexing doubled reference and aligning with minimap2 (map-ont).")

    mm2 = subprocess.Popen(
        [minimap2, "-ax", "map-ont", "-t", str(threads),
         str(doubled_path), str(fastq)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    sort = subprocess.Popen(
        [samtools, "sort", "-@", str(threads), "-o", str(output_bam)],
        stdin=mm2.stdout,
        stderr=subprocess.PIPE,
    )
    assert mm2.stdout is not None
    mm2.stdout.close()
    _, sort_err = sort.communicate()
    _, mm2_err = mm2.communicate()

    if mm2.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed (exit {mm2.returncode}):\n"
            f"{mm2_err.decode(errors='replace')}"
        )
    if sort.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed (exit {sort.returncode}):\n"
            f"{sort_err.decode(errors='replace')}"
        )

    idx = subprocess.run([samtools, "index", str(output_bam)], capture_output=True)
    if idx.returncode != 0:
        raise RuntimeError(
            f"samtools index failed (exit {idx.returncode}):\n"
            f"{idx.stderr.decode(errors='replace')}"
        )
