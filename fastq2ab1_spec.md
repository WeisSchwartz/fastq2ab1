# Spec Sheet: `fastq2ab1` — Nanopore BAM to Synthetic AB1 Converter

> Empirically verified against a real Plasmidsaurus AB1 file.

---

## 1. Overview

`fastq2ab1` is a Python command-line tool that converts nanopore sequencing data into synthetic AB1 chromatogram files. It computes per-position base counts from a read pileup and encodes these as trace signal intensities in the AB1 format, enabling visualization in standard sequence analysis software (SnapGene, ApE, Benchling, etc.).

The encoding method is confirmed by reverse-engineering a real Plasmidsaurus AB1 output:
- Trace values = raw pileup counts, globally scaled so the max depth position = 1000
- Quality values = linear function of coverage depth (not Phred)
- 4 trace points per base (flat/repeated), peak at index `4*i + 3`

**Important interpretive note for users:** The output AB1 is a *synthetic* chromatogram. Signal intensity represents *read frequency*, not fluorescence. A mixed signal means heterogeneity across reads (e.g., mixed plasmid population, indels), not instrument noise.

---

## 2. Inputs and Outputs

### Inputs

| Required files |
|----------------|
| Unaligned FASTQ (gzipped or plain) + reference FASTA |

> Pre-aligned BAM input was considered but removed. External BAMs aligned
> against the original linear reference produce a bell-curve coverage
> profile on circular plasmids (reads spanning the linearization breakpoint
> get soft-clipped at alignment time, and those bases are unrecoverable
> post-hoc). The FASTQ pipeline avoids this by aligning to a doubled
> reference internally.

### Output

| File | Description |
|------|-------------|
| `.ab1` | Synthetic chromatogram — confirmed compatible with SnapGene |
| `.fasta` (optional) | Consensus sequence derived from pileup majority base |
| Alignment stats | Coverage depth, alignment rate, low-coverage positions (stdout or log) |

---

## 3. Processing Pipeline

### Step 1: Alignment

- Use **minimap2** (subprocess) with `map-ont` preset
- Reference: user-supplied FASTA (single-sequence plasmid)
- **Circularity handling:** Before alignment, concatenate the reference with itself (2× length). Align reads against this doubled reference. The pileup stage collapses positions modulo the original reference length, so reads landing in either copy contribute to the equivalent position in copy 1 — eliminating the coverage taper at the linearization breakpoint.
- Output: sorted, indexed temp BAM
- Command pattern:
  ```
  minimap2 -ax map-ont -t <threads> ref_doubled.fasta reads.fastq | samtools sort -o aligned.bam
  samtools index aligned.bam
  ```

### Step 2: Pileup

- Use **pysam** to compute per-position base counts over the reference
- For each position, count A, T, G, C across all reads
- Exclude N calls, insertions, and deletions from counts (count only substitution-level base calls)
- Collapse the doubled-reference coordinate space modulo the original reference length: `count[i] += count_raw[j]` where `i = j % reference_length`
- Record `total_depth[i] = A[i] + T[i] + G[i] + C[i]` at each position
- Flag positions below minimum coverage threshold (default: 10×; configurable)

### Step 3: Consensus Sequence Generation

- Consensus base at each position = argmax(A, T, G, C)
- Ties: default to `N`; `--ambiguous-bases iupac` flag enables IUPAC ambiguity codes
- Low-coverage positions (below threshold): call as `N`
- The consensus becomes the PBAS2 sequence field in the AB1

### Step 4: AB1 Encoding

#### Trace scaling (confirmed from Plasmidsaurus data)

1. Find `global_max` = the single highest count value across all four channels and all positions
2. Compute `scale = 1000.0 / global_max`
3. For each channel, scaled value at position i = `round(raw_count[i] * scale)`
4. Each base position occupies exactly **4 trace points** with a `[0, 0, v, v]` shape (peak in the second half of the window — empirically matched to Plasmidsaurus AB1 files):
   - Trace indices for base i: `[4*i, 4*i+1, 4*i+2, 4*i+3]`
   - Values: `[0, 0, scaled_count, scaled_count]`
   - This produces the sharp-peak look in SnapGene; the previously documented `[v, v, v, v]` shape produced visually inferior flat-topped plateaus.
5. Total trace array length = `4 * reference_length`

#### PLOC1 (peak locations, confirmed)

- Peak index for base i = `4*i + 3`
- Array: `[3, 7, 11, 15, ..., 4*(N-1)+3]`
- Element type: `4` (short / int16). **Implementation note:** for references > ~8 kb, peak indices overflow int16 (max 32767). The writer auto-promotes to element type `5` (long / int32) in that case; BioPython reads both types transparently.

#### PCON2 (quality values, confirmed — NOT Phred)

- Quality at position i = `min(floor(total_depth[i] / 20), 50)`
- **Fixed divisor of 20** — empirically reverse-engineered from a real Plasmidsaurus AB1, where the `depth / PCON` ratio was 20.000 ± rounding across all 7,647 non-saturated positions. (An earlier guess of `depth / (max_depth/50)` was incorrect; that formula only coincides with the real one when `max_depth ≈ 1000`.)
- Saturates at depth ≥ 1000 (since `1000 // 20 = 50`)
- Zero-coverage positions: quality = 0

#### Channel-to-base mapping (confirmed)

| Tag | Base |
|-----|------|
| DATA9 | G |
| DATA10 | A |
| DATA11 | T |
| DATA12 | C |
| FWO_1 | `b"GATC"` (fixed) |

---

## 4. AB1 File Format — Required Tags

Confirmed minimal tag set from Plasmidsaurus (only these tags are required):

| Tag | Type | Description | Value |
|-----|------|-------------|-------|
| `DATA9` | short array | G channel trace | scaled G counts, 4x oversampled |
| `DATA10` | short array | A channel trace | scaled A counts, 4x oversampled |
| `DATA11` | short array | T channel trace | scaled T counts, 4x oversampled |
| `DATA12` | short array | C channel trace | scaled C counts, 4x oversampled |
| `PBAS1` | string | Called bases (copy of PBAS2) | consensus sequence |
| `PBAS2` | string | Called bases | consensus sequence |
| `PCON1` | char array | Quality (copy of PCON2) | same as PCON2 |
| `PCON2` | char array | Quality values | coverage-scaled, 0–50 |
| `PLOC1` | short array | Peak locations | `[3, 7, 11, ..., 4*(N-1)+3]` |
| `FWO_1` | string | Base order | `b"GATC"` |

Optional tags (add for usability, not required for SnapGene compatibility):
- `SMPL1`: sample name string
- `CMNT1`: free-text comment — record tool version, minimap2 version, mean coverage, reference filename, date

**Note:** Both PBAS1/PBAS2 and PCON1/PCON2 must be present (duplicate fields); SnapGene reads PBAS2/PCON2, older software reads PBAS1/PCON1.

---

## 5. AB1 Binary Format

AB1 is a Tagged Data File (TDF). All multi-byte integers are **big-endian**.

### File structure

```
[128-byte file header]
[Directory entries — one per tag]
[Data blocks]
```

### File header (128 bytes)

| Offset | Size | Value |
|--------|------|-------|
| 0 | 4 | Magic: `b"ABIF"` |
| 4 | 2 | Version: `101` (0x0065) |
| 6 | 4 | Tag name: `b"tdir"` |
| 10 | 4 | Tag number: `1` |
| 14 | 2 | Element type: `1023` (directory type) |
| 16 | 2 | Element size: `28` |
| 18 | 4 | Number of elements (= number of tags) |
| 22 | 4 | Data size (= num_tags * 28) |
| 26 | 4 | Data offset (byte offset to directory start) |
| 30 | 98 | Reserved (zeros) |

### Directory entry (28 bytes each)

| Offset | Size | Field |
|--------|------|-------|
| 0 | 4 | Tag name (4 ASCII chars, e.g. `b"DATA"`) |
| 4 | 4 | Tag number (e.g. `9` for DATA9) |
| 8 | 2 | Element type (see table below) |
| 10 | 2 | Element size in bytes |
| 12 | 4 | Number of elements |
| 16 | 4 | Data size in bytes (= element_size * num_elements) |
| 20 | 4 | Data offset (if data_size > 4, byte offset from file start; else data stored inline here) |
| 24 | 4 | Reserved (0) |

### Element types

| Code | Type | Size |
|------|------|------|
| 1 | byte | 1 |
| 2 | char (string) | 1 |
| 4 | short (int16) | 2 |
| 5 | long (int32) | 4 |
| 10 | date | 4 |
| 11 | time | 4 |
| 18 | pString (Pascal string) | variable |
| 19 | cString (C string) | variable |
| 1023 | directory | 28 |

### Tag-specific encoding

- `DATA9/10/11/12`: element type `4` (int16), big-endian array of scaled counts
- `PBAS1/2`: element type `2` (char), raw ASCII sequence bytes
- `PCON1/2`: element type `2` (char), raw quality byte array. The ABIF spec calls this type "byte" (1) for quality, but real Plasmidsaurus AB1 files declare it as `char` (2), and BioPython's `AbiIO` reader **requires** type 2 — it calls `.decode()` on PCON2 and crashes on type 1. The on-disk bytes are identical regardless.
- `PLOC1`: element type `4` (int16), peak location array
- `FWO_1`: element type `2` (char), `b"GATC"`
- `SMPL1`: element type `18` (pString) — length-prefixed: 1 byte length + ASCII chars
- `CMNT1`: element type `18` (pString)

---

## 6. Edge Cases and Error Handling

| Situation | Behavior |
|-----------|----------|
| Position below min coverage | Trace = 0 all channels, quality = 0, base = N; log warning with position index |
| Zero aligned reads | Abort with informative error, non-zero exit code |
| Alignment rate < 50% | Warn; continue unless `--strict` flag set |
| Reference FASTA has multiple records | Abort: "Reference must be a single sequence. Use --contig to specify." |
| FASTQ mean read length < 200bp | Warn: tool designed for nanopore long reads |
| Ties in consensus base call | Default: N; `--ambiguous-bases iupac` for IUPAC codes |
| global_max = 0 (no aligned reads) | Abort before division by zero |

---

## 7. Implementation Notes

### Language and dependencies

| Dependency | Purpose |
|------------|---------|
| Python ≥ 3.9 | Language |
| `pysam` | BAM parsing and pileup |
| `biopython` | FASTA parsing, IUPAC, AB1 read-back validation |
| `minimap2` (system binary) | Alignment |
| `samtools` (system binary) | BAM sort/index |
| `struct` (stdlib) | AB1 binary writer |
| `argparse` (stdlib) | CLI |
| `logging` (stdlib) | Structured logging |

### AB1 writer

- Implement from scratch using Python `struct` — BioPython has no AB1 writer
- Write a round-trip validation test: write AB1 → read back with `Bio.SeqIO.read(..., "abi")` → assert all tag values match inputs
- All integers big-endian (`>` format in struct)
- For data blocks > 4 bytes: write data after directory, store byte offset in directory entry
- For data ≤ 4 bytes: store inline in the offset field

### Modularity (required for future web deployment)

Keep business logic fully separated from CLI:

```
fastq2ab1/
  cli.py          # argparse only; calls pipeline.run()
  pipeline.py     # orchestrates full workflow; importable as library
  aligner.py      # minimap2 subprocess wrapper + circularity handling
  pileup.py       # BAM → per-position count arrays via pysam
  consensus.py    # count arrays → consensus sequence + quality values
  ab1_writer.py   # count arrays + consensus → AB1 binary
  utils.py        # logging, temp file management
```

`pipeline.run(reference, output, fastq, **kwargs)` must be callable without CLI invocation.

### Temporary files

- Use Python `tempfile.TemporaryDirectory()` as context manager
- Cleaned up on exit and on failure
- `--keep-temp` flag disables cleanup for debugging

---

## 8. CLI Interface

```
usage: fastq2ab1 [-h] -r REFERENCE -f FASTQ -o OUTPUT
               [--sample-name NAME] [--min-coverage INT]
               [--threads INT] [--ambiguous-bases {N,iupac}]
               [--write-consensus FASTA] [--strict] [--keep-temp] [--verbose]

Convert nanopore sequencing data to synthetic AB1 chromatogram.

Required:
  -f, --fastq           Unaligned FASTQ file (gzipped or plain)
  -r, --reference       Reference plasmid FASTA (single sequence)
  -o, --output          Output .ab1 file path

Options:
  --sample-name         Sample name embedded in AB1 SMPL1 tag (default: input filename stem)
  --min-coverage        Minimum read depth per position for base calling (default: 10)
  --threads             Threads for minimap2 (default: 4)
  --ambiguous-bases     Tie handling: 'N' or 'iupac' (default: N)
  --write-consensus     Also write consensus sequence to this FASTA path
  --strict              Abort if alignment rate < 50% (default: warn only)
  --keep-temp           Retain temporary files for debugging
  --verbose             Verbose logging to stderr
```

### Example usage

```bash
# Minimal
fastq2ab1 -r plasmid.fasta -f reads.fastq.gz -o result.ab1

# With a sample label
fastq2ab1 -r plasmid.fasta -f reads.fastq.gz -o result.ab1 --sample-name "Clone_A3"

# With consensus output and custom coverage threshold
fastq2ab1 -r plasmid.fasta -f reads.fastq -o result.ab1 \
        --min-coverage 20 --write-consensus consensus.fasta
```

---

## 9. Testing Requirements

| Test | Description |
|------|-------------|
| Unit: pileup | Synthetic BAM with known reads → assert expected count arrays |
| Unit: scaling | Known counts → assert scaled values, global max = 1000, trace shape `[0,0,v,v]` per base |
| Unit: quality | Known depth values → assert quality = `min(depth // 20, 50)`, independent of max depth |
| Unit: PLOC | N-base consensus → PLOC = [3, 7, 11, ..., 4*(N-1)+3] |
| Unit: AB1 writer round-trip | Write AB1, read back with BioPython, assert all tag values match |
| Unit: AB1 writer overflow | Reference > 8 kb → PLOC1 auto-promoted to int32, round-trips correctly |
| Unit: circularity | Reads aligned to copy 2 of a doubled reference collapse into copy-1 positions |
| Integration: pileup → consensus → AB1 | Synthetic BAM → full internal stack → AB1 readable by BioPython with expected tag values |
| Edge: zero-coverage positions | Correct zero encoding, warning logged |
| Edge: multi-record FASTA | Non-zero exit, informative error message |
| Regression: Plasmidsaurus | Given the provided sample FASTQ + FASTA, output AB1 should match the Plasmidsaurus-provided AB1 in sequence (100%), PLOC1, FWO_1; quality and DATA trace shape match by formula |

---

## 10. Out of Scope (v1)

- FAST5 / POD5 raw signal input (requires Dorado basecaller; document as external prerequisite)
- Batch / multi-sample processing
- GUI (modular design enables future wrapping)
- Indel representation in traces
- VCF output
- Methylation signal from BAM MM/ML tags

---

## 11. Future: Web Deployment

The `pipeline.run()` function as a library enables a thin web wrapper:
- FastAPI or Flask endpoint accepts file uploads (BAM or FASTQ + FASTA reference)
- Calls `pipeline.run(...)` with user-supplied parameters
- Returns AB1 as a file download
- Stateless: no shared state between requests
- Use per-request temp directories, cleaned up after response completes
