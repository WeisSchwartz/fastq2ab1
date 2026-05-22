# fastq2ab1

[![test](https://github.com/WeisSchwartz/fastq2ab1/actions/workflows/test.yml/badge.svg)](https://github.com/WeisSchwartz/fastq2ab1/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Convert nanopore FASTQ reads into a synthetic AB1 chromatogram (`.ab1`)
that opens in SnapGene, ApE, Benchling, and other standard
sequence-analysis tools.

The output AB1 is a **synthetic** chromatogram. Signal intensity represents
*read frequency*, not fluorescence. A mixed signal at a position means
heterogeneity across reads (e.g. mixed plasmid population, indel hotspot) —
not instrument noise.

Inspired by Plasmidsaurus Technical Documentation (https://plasmidsaurus.com).
Code written with the help of Claude Code (https://claude.ai).

---

## What it does

For each reference position, `fastq2ab1`:

1. Aligns reads against the reference with minimap2 (`map-ont` preset).
2. Counts the number of reads supporting each of A / T / G / C
   (substitution-level only; insertions, deletions, and N base calls are
   excluded).
3. Encodes those counts as a four-channel chromatogram trace, scaled so the
   single highest count across all positions and channels maps to 1000.
4. Calls a consensus base at each position (argmax; ties → `N` or IUPAC).
5. Writes an ABIF/AB1 binary with the minimal Plasmidsaurus-compatible tag
   set (`DATA9/10/11/12`, `PBAS1/2`, `PCON1/2`, `PLOC1`, `FWO_1`, plus
   optional `SMPL1` and `CMNT1`).

For circular plasmid references the reference is duplicated (2× length)
before alignment so reads spanning the linearization breakpoint align
end-to-end; the pileup then collapses coordinates modulo the original
reference length, so coverage stays flat across the breakpoint.

---

## Prerequisites

### Python
- Python ≥ 3.9
- [`pysam`](https://github.com/pysam-developers/pysam) (pileup over the
  intermediate BAM that minimap2 produces)
- [`biopython`](https://biopython.org/) (FASTA parsing, AB1 round-trip
  validation)

### External binaries
- [`minimap2`](https://github.com/lh3/minimap2) — long-read aligner
- [`samtools`](https://www.htslib.org/) — BAM sort & index

Both must be on your `PATH`.

### Recommended: conda/mamba environment
```bash
mamba create -n fastq2ab1 -c bioconda -c conda-forge \
    python=3.11 pysam biopython minimap2 samtools pytest
mamba activate fastq2ab1
```

---

## Installation

From the repo root:
```bash
pip install -e .
```
This installs the `fastq2ab1` console script and exposes
`from fastq2ab1 import run` as a library entry point.

---

## How to use

### Command line

```text
usage: fastq2ab1 [-h] -r REFERENCE -f FASTQ -o OUTPUT
               [--sample-name NAME] [--min-coverage INT] [--threads INT]
               [--ambiguous-bases {N,iupac}] [--write-consensus FASTA]
               [--strict] [--keep-temp] [--verbose] [--version]
```

**Required**

| Flag | Description |
|---|---|
| `-f, --fastq` | Unaligned FASTQ (gzipped or plain). |
| `-r, --reference` | Reference plasmid FASTA. Must be a single sequence. |
| `-o, --output` | Output `.ab1` path. |

**Options**

| Flag | Default | Description |
|---|---|---|
| `--sample-name` | input stem | Sample label written to AB1 `SMPL1`. |
| `--min-coverage` | `10` | Positions below this depth are called `N` and logged. |
| `--threads` | `4` | Threads for minimap2. |
| `--ambiguous-bases` | `N` | Tie handling. `iupac` uses IUPAC ambiguity codes (R, Y, W, …). |
| `--write-consensus` | off | Also write consensus to this FASTA path. |
| `--strict` | off | Abort if alignment rate < 50%. |
| `--keep-temp` | off | Keep temp files (doubled-reference FASTA, aligned BAM) for debugging. |
| `--verbose` | off | DEBUG logging to stderr. |

> **Why FASTQ-only?** Earlier versions also accepted pre-aligned BAMs, but
> external BAMs aligned to the original linear plasmid reference produce a
> bell-curve coverage profile — reads spanning the linearization breakpoint
> get soft-clipped at alignment time, and those bases cannot be recovered
> downstream. The FASTQ pipeline doubles the reference internally before
> alignment, then collapses positions modulo `reference_length`, which keeps
> coverage flat across the breakpoint.

### Examples

```bash
# Minimal
fastq2ab1 -r plasmid.fasta -f reads.fastq.gz -o result.ab1

# With a sample label
fastq2ab1 -r plasmid.fasta -f reads.fastq.gz -o result.ab1 \
        --sample-name "Clone_A3"

# Also dump the consensus and bump the coverage threshold
fastq2ab1 -r plasmid.fasta -f reads.fastq -o result.ab1 \
        --min-coverage 20 --write-consensus consensus.fasta
```

### Exit codes

| Code | Meaning |
|---|---|
| `0` | Success |
| `2` | File not found |
| `3` | Bad input (multi-record FASTA, no aligned reads, …) |
| `4` | Runtime failure (alignment failure, `--strict` alignment-rate abort, …) |

### Library usage

`pipeline.run()` is fully usable without the CLI — useful for web wrappers
or notebook workflows:

```python
from fastq2ab1 import run

stats = run(
    reference="plasmid.fasta",
    fastq="reads.fastq.gz",
    output="result.ab1",
    sample_name="Clone_A3",
    min_coverage=15,
    write_consensus="consensus.fasta",
)
print(stats.mean_coverage, stats.n_low_coverage)
```

`run()` returns a `RunStats` dataclass with reference length, read counts,
alignment rate, mean coverage, number of low-coverage positions, and the
global max values used for trace and quality scaling.

---

## Encoding details

These match the encoding observed in real Plasmidsaurus-produced AB1 files
and are what makes the output open cleanly in SnapGene:

- **Trace scaling.** `scale = 1000 / global_max`, where `global_max` is the
  single highest count across all four channels and all positions. Each
  channel value at position *i* is `round(count[i] * scale)`.
- **4× oversampling, sharp-peak shape.** Each base contributes 4 samples
  with values `[0, 0, v, v]` — the leading two zeros produce the sharp-peak
  look in SnapGene. (The original spec described flat `[v, v, v, v]`
  plateaus, but real Plasmidsaurus AB1 files use the sharper shape.) Total
  trace length = `4 × reference_length`. Peak index for base *i* is
  `4*i + 3` (stored in `PLOC1`).
- **Channel → base.** `DATA9 = G`, `DATA10 = A`, `DATA11 = T`, `DATA12 = C`;
  `FWO_1 = b"GATC"`.
- **Quality (PCON1/2).** *Not* Phred. Fixed-divisor mapping of coverage
  depth to 0–50: `quality[i] = min(depth[i] // 20, 50)`. Saturates at
  depth ≥ 1000. This formula was reverse-engineered from a real
  Plasmidsaurus AB1 — the `depth / PCON` ratio was 20.000 ± rounding across
  all 7,647 non-saturated positions.
- **Consensus.** argmax(A, T, G, C) per position. Ties default to `N`, or
  IUPAC ambiguity codes with `--ambiguous-bases iupac`. Positions below
  `--min-coverage` are called `N`.

> **Spec note.** The AB1 spec says PCON tags use ABIF element type `1`
> (byte). In practice, real Plasmidsaurus AB1 files declare PCON1/2 as
> type `2` (char) — BioPython's `AbiIO` reader requires this and will
> crash on type-1 quality arrays. `fastq2ab1` writes type `2`; the on-disk
> byte sequence is identical either way.
>
> For references > ~8 kb, the spec's `PLOC1` type 4 (int16) overflows
> (peak index 4*i+3 exceeds 32767). `fastq2ab1` auto-promotes `PLOC1` to
> type 5 (int32) in that case. BioPython reads both types transparently.

---

## Project layout

```
fastq2ab1/
  cli.py          # argparse only; calls pipeline.run()
  pipeline.py     # orchestrates the full workflow; importable as library
  aligner.py      # minimap2 subprocess wrapper + circularity handling
  pileup.py       # BAM → per-position count arrays via pysam
  consensus.py    # count arrays → consensus + quality + scaled traces
  ab1_writer.py   # from-scratch ABIF binary writer (big-endian struct)
  utils.py        # logging, binary checks
tests/
  test_consensus.py            # scaling, PLOC, quality formula, tie handling
  test_ab1_roundtrip.py        # write AB1 → read with BioPython → match (+ int16 overflow)
  test_circularity.py          # doubled-reference collapse for circular plasmids
  test_pileup_integration.py   # synthetic BAM → pileup → consensus → AB1
```

Run the suite:
```bash
pytest tests/
```

---

## Out of scope (v1)

- FAST5 / POD5 raw signal input (use Dorado upstream).
- Batch / multi-sample processing.
- GUI.
- Indel representation in traces.
- VCF output.
- Methylation signal from BAM `MM`/`ML` tags.

See [`fastq2ab1_spec.md`](fastq2ab1_spec.md) for the full specification.
