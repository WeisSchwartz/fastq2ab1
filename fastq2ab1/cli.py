"""argparse front-end. All logic lives in pipeline.run()."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .pipeline import run
from .utils import configure_logging


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fastq2ab1",
        description="Convert nanopore sequencing data to synthetic AB1 chromatogram.",
    )
    p.add_argument("-f", "--fastq", required=True,
                   help="Unaligned FASTQ file (gzipped or plain).")
    p.add_argument("-r", "--reference", required=True,
                   help="Reference plasmid FASTA (single sequence).")
    p.add_argument("-o", "--output", required=True, help="Output .ab1 file path.")

    p.add_argument("--sample-name", default=None,
                   help="Sample name embedded in AB1 SMPL1 tag (default: input filename stem).")
    p.add_argument("--min-coverage", type=int, default=10,
                   help="Minimum read depth per position for base calling (default: 10).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for minimap2 (default: 4).")
    p.add_argument("--ambiguous-bases", choices=("N", "iupac"), default="N",
                   help="Tie handling: 'N' or 'iupac' (default: N).")
    p.add_argument("--write-consensus", default=None, metavar="FASTA",
                   help="Also write consensus sequence to this FASTA path.")
    p.add_argument("--strict", action="store_true",
                   help="Abort if alignment rate < 50%% (default: warn only).")
    p.add_argument("--keep-temp", action="store_true",
                   help="Retain temporary files for debugging.")
    p.add_argument("--verbose", action="store_true",
                   help="Verbose logging to stderr.")
    p.add_argument("--version", action="version", version=f"fastq2ab1 {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    log = configure_logging(verbose=args.verbose)

    try:
        stats = run(
            reference=args.reference,
            output=args.output,
            fastq=args.fastq,
            sample_name=args.sample_name,
            min_coverage=args.min_coverage,
            threads=args.threads,
            ambiguous_bases=args.ambiguous_bases,
            write_consensus=args.write_consensus,
            strict=args.strict,
            keep_temp=args.keep_temp,
        )
    except FileNotFoundError as e:
        log.error("File not found: %s", e)
        return 2
    except ValueError as e:
        log.error("%s", e)
        return 3
    except RuntimeError as e:
        log.error("%s", e)
        return 4

    log.info(
        "Done. Wrote %s — %d bp, %d aligned reads, mean cov %.1fx, %d low-cov positions.",
        stats.output_path,
        stats.reference_length,
        stats.n_reads_aligned,
        stats.mean_coverage,
        stats.n_low_coverage,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
