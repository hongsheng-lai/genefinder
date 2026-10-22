#!/usr/bin/env python3
"""Build a positional offset index for a BGZF-compressed PAF, so that only the records overlapping a gene
need to be fetched (HTTP range requests on the remote PAF) instead of downloading the whole file.

Usage: curl -s PAF_URL | build_paf_index.py > offsets.bed      (or: build_paf_index.py local.paf.gz)
then:  LC_ALL=C sort -k1,1 -k2,2n offsets.bed | bgzip > paf_offsets.bed.gz && tabix -p bed paf_offsets.bed.gz

One line per PAF record (BED-like, 0-based, sorted by the target = CHM13 interval):
  target  tstart  tend  cstart  cend  ustart  ulen
  cstart  compressed offset of the BGZF block where the record starts
  cend    compressed offset just after the BGZF block where the record ends -> fetch bytes [cstart, cend)
  ustart  offset of the record inside the uncompressed data of the block at cstart
  ulen    record length in bytes (without the newline)
The input is read once, sequentially; nothing but the index is written.
"""
import struct
import sys
import zlib

PREFIX_MAX = 512  # bytes of each line kept to parse the first 9 PAF columns


def blocks(fh):
    """Yield (compressed_offset, uncompressed_bytes, next_compressed_offset) for each BGZF block."""
    coff = 0
    while True:
        head = fh.read(18)
        if not head:
            return
        if len(head) < 18 or head[:4] != b"\x1f\x8b\x08\x04" or head[12:14] != b"BC":
            sys.exit(f"not a BGZF file (bad block header at compressed offset {coff})")
        bsize = struct.unpack("<H", head[16:18])[0] + 1
        rest = fh.read(bsize - 18)
        data = zlib.decompress(rest[:-8], -15)
        yield coff, data, coff + bsize
        coff += bsize


def main():
    src = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != "-" else None
    fh = open(src, "rb") if src else sys.stdin.buffer
    out = sys.stdout
    start = None          # (cstart, ustart) of the current line
    prefix = bytearray()
    length = 0
    n = 0
    for coff, data, nxt in blocks(fh):
        if not data:
            continue      # EOF marker block
        i = 0
        while True:
            if start is None:
                start = (coff, i)
            j = data.find(b"\n", i)
            seg_end = len(data) if j < 0 else j
            if len(prefix) < PREFIX_MAX:
                prefix += data[i:min(seg_end, i + PREFIX_MAX - len(prefix))]
            length += seg_end - i
            if j < 0:
                break
            f = prefix.split(b"\t", 9)
            if len(f) >= 9:
                out.write(f"{f[5].decode()}\t{int(f[7])}\t{int(f[8])}\t{start[0]}\t{nxt}\t{start[1]}\t{length}\n")
                n += 1
            start, prefix, length = None, bytearray(), 0
            i = j + 1
            if i >= len(data):
                break
    sys.stderr.write(f"indexed {n} PAF records\n")


if __name__ == "__main__":
    main()
