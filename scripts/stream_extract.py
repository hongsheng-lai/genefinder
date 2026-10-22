#!/usr/bin/env python3
"""Extract regions from a FASTA stream in one sequential pass (the assembly is never stored).

Usage: curl URL | bgzip -dc | stream_extract.py --regions regions.txt [--rc] > out.fa
regions.txt: one samtools-style region per line, contig:start-end (1-based, inclusive),
optionally followed by a tab and "-" to reverse-complement that region.
Output records are named exactly like the region, with "/rc" appended for reverse-complemented ones;
--rc reverse-complements every region and does not append "/rc" (like `samtools faidx -i --mark-strand no`). Only contigs that carry a region are kept
in memory, one at a time.
"""
import argparse
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", required=True)
    ap.add_argument("--rc", action="store_true")
    a = ap.parse_args()
    want = defaultdict(list)
    for line in open(a.regions):
        f = line.split()
        if f:
            reg = f[0]
            c, se = reg.rsplit(":", 1)
            s, e = se.split("-")
            want[c].append((int(s), int(e), reg, len(f) > 1 and f[1] == "-"))
    out = sys.stdout
    found = set()

    def emit(name, chunks):
        seq = "".join(chunks)
        for s, e, reg, rc in want[name]:
            sub = seq[s - 1:e]
            if a.rc or rc:
                sub = sub.translate(COMP)[::-1]
            out.write(f">{reg}{'/rc' if rc and not a.rc else ''}\n")
            for i in range(0, len(sub), 60):
                out.write(sub[i:i + 60] + "\n")
        found.add(name)

    name, chunks, keep = None, [], False
    for line in sys.stdin:
        if line.startswith(">"):
            if keep:
                emit(name, chunks)
            name = line[1:].split()[0]
            keep = name in want
            chunks = []
        elif keep:
            chunks.append(line.rstrip())
    if keep:
        emit(name, chunks)
    lost = set(want) - found
    if lost:
        sys.exit(f"contigs not found in stream: {','.join(sorted(lost)[:5])}")


if __name__ == "__main__":
    main()
