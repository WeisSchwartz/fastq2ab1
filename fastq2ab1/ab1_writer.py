"""Minimal AB1 (ABIF) writer.

Emits the tag set required by SnapGene/Plasmidsaurus-compatible chromatograms:
DATA9/10/11/12, PBAS1/2, PCON1/2, PLOC1, FWO_1, and optional SMPL1, CMNT1.

All multi-byte integers are big-endian. The file layout is:

    [128-byte header]
    [data blocks, in tag-write order]
    [directory entries (28 bytes each)]

The directory is written after the data so that data-block offsets are known
before we serialize directory entries.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

# ABIF element type codes
ET_BYTE = 1
ET_CHAR = 2
ET_SHORT = 4
ET_LONG = 5
ET_PSTRING = 18


@dataclass
class _DirEntry:
    tag_name: bytes   # 4 ASCII
    tag_number: int   # int32
    element_type: int # int16
    element_size: int # int16
    num_elements: int # int32
    data_size: int    # int32
    data_offset: int  # int32 — file offset, or inline value if data_size <= 4
    reserved: int = 0


@dataclass
class AB1Builder:
    entries: List[Tuple[_DirEntry, bytes]] = field(default_factory=list)

    # ---- low-level append ----
    def _add(
        self,
        tag_name: bytes,
        tag_number: int,
        element_type: int,
        element_size: int,
        num_elements: int,
        data: bytes,
    ) -> None:
        assert len(tag_name) == 4
        data_size = element_size * num_elements
        assert data_size == len(data), (
            f"data_size mismatch for {tag_name!r}{tag_number}: "
            f"declared {data_size}, got {len(data)}"
        )
        entry = _DirEntry(
            tag_name=tag_name,
            tag_number=tag_number,
            element_type=element_type,
            element_size=element_size,
            num_elements=num_elements,
            data_size=data_size,
            data_offset=0,  # filled in later
        )
        self.entries.append((entry, data))

    # ---- typed tag helpers ----
    def add_short_array(
        self, tag_name: bytes, tag_number: int, values: Sequence[int]
    ) -> None:
        data = struct.pack(f">{len(values)}h", *values)
        self._add(tag_name, tag_number, ET_SHORT, 2, len(values), data)

    def add_int_array_auto(
        self, tag_name: bytes, tag_number: int, values: Sequence[int]
    ) -> None:
        """Write an integer array using the smallest ABIF type that fits.

        Used for tags like PLOC1 whose spec type is `short` (int16) but whose
        values overflow int16 on long references (peak index = 4*i+3, so any
        reference > ~8192 bp triggers overflow). When values don't fit in
        int16 we fall back to element type 5 (long / int32), which is part
        of the ABIF type table and is parsed transparently by BioPython.
        """
        if not values:
            self._add(tag_name, tag_number, ET_SHORT, 2, 0, b"")
            return
        vmax = max(values)
        vmin = min(values)
        if -32768 <= vmin and vmax <= 32767:
            data = struct.pack(f">{len(values)}h", *values)
            self._add(tag_name, tag_number, ET_SHORT, 2, len(values), data)
        else:
            data = struct.pack(f">{len(values)}i", *values)
            self._add(tag_name, tag_number, ET_LONG, 4, len(values), data)

    def add_byte_array(
        self, tag_name: bytes, tag_number: int, values: Sequence[int]
    ) -> None:
        data = bytes(int(v) & 0xFF for v in values)
        self._add(tag_name, tag_number, ET_BYTE, 1, len(values), data)

    def add_quality_array(
        self, tag_name: bytes, tag_number: int, values: Sequence[int]
    ) -> None:
        """PCON tags. On disk these are raw 8-bit values, but the AB1 directory
        entry must declare element type 2 ('char') — BioPython's AbiIO reader
        calls `.decode()` on PCON2 and will crash if the type is 1 (byte). The
        bytes themselves are identical to a byte array; only the type code
        differs."""
        data = bytes(int(v) & 0xFF for v in values)
        self._add(tag_name, tag_number, ET_CHAR, 1, len(values), data)

    def add_char_string(
        self, tag_name: bytes, tag_number: int, value: str | bytes
    ) -> None:
        if isinstance(value, str):
            value = value.encode("ascii", errors="replace")
        self._add(tag_name, tag_number, ET_CHAR, 1, len(value), value)

    def add_pstring(
        self, tag_name: bytes, tag_number: int, value: str | bytes
    ) -> None:
        if isinstance(value, str):
            value = value.encode("ascii", errors="replace")
        if len(value) > 255:
            value = value[:255]
        data = bytes([len(value)]) + value
        # ABIF spec: pString element_size is 1; num_elements is total byte
        # length INCLUDING the length prefix. (This is the standard
        # convention BioPython expects when round-tripping.)
        self._add(tag_name, tag_number, ET_PSTRING, 1, len(data), data)

    # ---- serialization ----
    def build(self) -> bytes:
        num_tags = len(self.entries)
        header_size = 128
        # Data blocks come immediately after the 128-byte header.
        # Inline (data_size <= 4) blocks contribute nothing to the data
        # region; only larger blocks consume space.
        cursor = header_size
        data_regions: List[bytes] = []
        for entry, data in self.entries:
            if entry.data_size <= 4:
                # Inline: store left-justified in a 4-byte field.
                padded = data + b"\x00" * (4 - len(data))
                # data_offset holds the inline 4 bytes interpreted big-endian.
                entry.data_offset = struct.unpack(">i", padded)[0]
            else:
                entry.data_offset = cursor
                data_regions.append(data)
                cursor += len(data)
        # Directory follows the last data block.
        dir_offset = cursor
        dir_size = num_tags * 28

        out = bytearray()

        # ---- 128-byte file header ----
        header = bytearray(128)
        header[0:4] = b"ABIF"
        struct.pack_into(">h", header, 4, 101)        # version
        header[6:10] = b"tdir"                          # tag name
        struct.pack_into(">i", header, 10, 1)           # tag number
        struct.pack_into(">h", header, 14, 1023)        # element type (directory)
        struct.pack_into(">h", header, 16, 28)          # element size
        struct.pack_into(">i", header, 18, num_tags)    # num elements
        struct.pack_into(">i", header, 22, dir_size)    # data size
        struct.pack_into(">i", header, 26, dir_offset)  # data offset
        # bytes 30..128 stay zero
        out.extend(header)

        # ---- data blocks ----
        for region in data_regions:
            out.extend(region)

        # ---- directory entries ----
        for entry, _ in self.entries:
            dirent = struct.pack(
                ">4siHHiiii",
                entry.tag_name,
                entry.tag_number,
                entry.element_type,
                entry.element_size,
                entry.num_elements,
                entry.data_size,
                entry.data_offset,
                entry.reserved,
            )
            out.extend(dirent)

        return bytes(out)


def write_ab1(
    output_path: str,
    sequence: str,
    quality: Sequence[int],
    trace_G: Sequence[int],
    trace_A: Sequence[int],
    trace_T: Sequence[int],
    trace_C: Sequence[int],
    peak_locations: Sequence[int],
    sample_name: Optional[str] = None,
    comment: Optional[str] = None,
) -> None:
    """Serialize the AB1 file with the required Plasmidsaurus-compatible tag set."""
    if not (len(trace_G) == len(trace_A) == len(trace_T) == len(trace_C)):
        raise ValueError("Trace channels must have equal length.")
    if len(sequence) != len(quality):
        raise ValueError("Sequence and quality must be the same length.")
    if len(sequence) != len(peak_locations):
        raise ValueError("Sequence and peak_locations must be the same length.")

    b = AB1Builder()
    # Order chosen so trace data blocks come first (largest), then strings.
    b.add_short_array(b"DATA", 9, trace_G)
    b.add_short_array(b"DATA", 10, trace_A)
    b.add_short_array(b"DATA", 11, trace_T)
    b.add_short_array(b"DATA", 12, trace_C)
    # PLOC1: short array per the spec, but auto-promote to int32 if peak
    # indices overflow int16 (happens for references > ~8 kb, since each
    # base consumes 4 trace samples and the peak sits at 4*i+3).
    b.add_int_array_auto(b"PLOC", 1, peak_locations)
    b.add_char_string(b"PBAS", 1, sequence)
    b.add_char_string(b"PBAS", 2, sequence)
    b.add_quality_array(b"PCON", 1, quality)
    b.add_quality_array(b"PCON", 2, quality)
    b.add_char_string(b"FWO_", 1, b"GATC")
    if sample_name:
        b.add_pstring(b"SMPL", 1, sample_name)
    if comment:
        b.add_pstring(b"CMNT", 1, comment)

    blob = b.build()
    with open(output_path, "wb") as fh:
        fh.write(blob)
