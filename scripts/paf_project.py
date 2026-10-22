#!/usr/bin/env python3
"""Project CHM13 (PAF target) BED intervals onto haplotype contigs (PAF query) by walking cg:Z CIGARs.

Fallback for impg, and the engine of the "stream" mode
(curl | bgzip -dc | awk-prefilter | paf_project.py).

Usage: paf_project.py --bed regions.bed [PAF|-]  > projected.paf

Output is PAF (12 columns + tags), one line per (alignment record x BED interval) overlap:
  query = haplotype contig, projected interval
  target = CHM13 interval actually covered by that record (clipped to the BED interval)
  col 10/11 = matches / alignment block length within the clipped piece
  tags: cg:Z (clipped CIGAR), gi:f (gap-compressed identity of clipped piece), ch:Z (chain id), rn:Z (BED name)
This is the same layout as `impg query -o paf`, so downstream code handles both.
"""
import argparse
import re
import sys
from collections import defaultdict
from bisect import bisect_right

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def load_bed(path):
    by_chrom = defaultdict(list)
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.rstrip("\n").split("\t")
        by_chrom[f[0]].append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}"))
    idx = {}
    for c, iv in by_chrom.items():
        iv.sort()
        # max end prefix for overlap search
        starts = [s for s, _, _ in iv]
        idx[c] = (iv, starts, max(e - s for s, e, _ in iv))
    return idx


def overlaps(idx, chrom, s, e):
    if chrom not in idx:
        return []
    iv, starts, maxlen = idx[chrom]
    lo = bisect_right(starts, s - maxlen - 1)
    hi = bisect_right(starts, e - 1)
    return [x for x in iv[lo:hi] if x[1] > s and x[0] < e]


class _Piece:
    """Accumulates the clipped projection of one target interval [cs, ce) during a CIGAR walk."""
    __slots__ = ("cs", "ce", "q0", "q1", "t0", "t1", "ops", "m", "x", "ins", "dele", "gaps")

    def __init__(self, cs, ce):
        self.cs, self.ce = cs, ce
        self.q0 = self.q1 = self.t0 = self.t1 = None
        self.ops = []
        self.m = self.x = self.ins = self.dele = self.gaps = 0

    def step(self, op, n, t, off):
        cs, ce = self.cs, self.ce
        if op in "M=X":
            a, b = max(t, cs), min(t + n, ce)
            if a < b:
                if self.q0 is None:
                    self.q0 = off + (a - t); self.t0 = a
                self.q1 = off + (b - t); self.t1 = b
                self.ops.append((b - a, op))
                if op == "X":
                    self.x += b - a
                else:
                    self.m += b - a
        elif op in "DN":
            a, b = max(t, cs), min(t + n, ce)
            if a < b and self.q0 is not None:  # deletions interior to the clipped piece
                self.ops.append((b - a, "D")); self.dele += b - a; self.gaps += 1
        elif op == "I":
            if cs < t < ce and self.q0 is not None:
                self.ops.append((n, "I")); self.ins += n; self.gaps += 1

    def result(self, qs, qe, strand):
        if self.q0 is None:
            return None
        ops = self.ops
        # trailing deletions/insertions after the last aligned base do not belong to the piece
        while ops and ops[-1][1] in "DI":
            n, op = ops.pop()
            if op == "D":
                self.dele -= n
            else:
                self.ins -= n
            self.gaps -= 1
        if strand == "+":
            q_lo, q_hi = qs + self.q0, qs + self.q1
        else:
            q_lo, q_hi = qe - self.q1, qe - self.q0
        merged = []
        for n, op in ops:
            if merged and merged[-1][1] == op:
                merged[-1][0] += n
            else:
                merged.append([n, op])
        return (q_lo, q_hi, self.t0, self.t1, "".join(f"{n}{op}" for n, op in merged),
                self.m, self.x, self.ins, self.dele, self.gaps)


def project_many(cigar, tstart, qs, qe, strand, intervals):
    """One CIGAR walk for many target intervals [(cs, ce), ...].
    Returns a list (same order) of (q_lo, q_hi, t_lo, t_hi, clipped_cigar, matches, mism, ins, dele, n_gapopen)
    or None where the interval has no aligned base."""
    order = sorted(range(len(intervals)), key=lambda i: intervals[i][0])
    pieces = [_Piece(cs, ce) for cs, ce in intervals]
    t = tstart
    off = 0  # query offset from the alignment's start in alignment direction
    j = 0
    active = []
    for mo in CIG.finditer(cigar):
        n = int(mo.group(1)); op = mo.group(2)
        span = n if op in "M=XDN" else 0
        while j < len(order) and intervals[order[j]][0] <= t + span:
            active.append(pieces[order[j]]); j += 1
        if active:
            active = [p for p in active if p.ce > t]
        if not active and j == len(order):
            break
        for p in active:
            p.step(op, n, t, off)
        if op in "M=X":
            t += n; off += n
        elif op in "DN":
            t += n
        elif op == "I":
            off += n
        # S/H/P ignored (not used in these PAFs)
    return [p.result(qs, qe, strand) for p in pieces]


def project(cigar, tstart, qs, qe, strand, cs, ce):
    return project_many(cigar, tstart, qs, qe, strand, [(cs, ce)])[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("paf", nargs="?", default="-")
    a = ap.parse_args()
    idx = load_bed(a.bed)
    fh = sys.stdin if a.paf == "-" else open(a.paf)
    out = sys.stdout
    for line in fh:
        f = line.rstrip("\n").split("\t", 12)  # don't split the (large) tag/CIGAR part unless needed
        if len(f) < 12:
            continue
        tname, ts, te = f[5], int(f[7]), int(f[8])
        hits = overlaps(idx, tname, ts, te)
        if not hits:
            continue
        tags = {t[:2]: t[5:] for t in (f[12].split("\t") if len(f) > 12 else [])}
        cg = tags.get("cg")
        if cg is None:
            sys.stderr.write(f"no cg:Z on record {f[0]}:{f[2]}-{f[3]}, skipped\n")
            continue
        qs, qe, strand = int(f[2]), int(f[3]), f[4]
        results = project_many(cg, ts, qs, qe, strand, [(max(bs, ts), min(be, te)) for bs, be, _ in hits])
        for (bs, be, name), r in zip(hits, results):
            if r is None:
                continue
            q_lo, q_hi, t_lo, t_hi, ccig, m, x, ins, dele, gaps = r
            alen = m + x + ins + dele
            gi = m / (m + x + gaps) if (m + x + gaps) else 0.0
            out.write(f"{f[0]}\t{f[1]}\t{q_lo}\t{q_hi}\t{strand}\t{tname}\t{f[6]}\t{t_lo}\t{t_hi}"
                      f"\t{m}\t{alen}\t{f[11]}\tgi:f:{gi:.6f}\tch:Z:{tags.get('ch', 'NA')}\trn:Z:{name}\tcg:Z:{ccig}\n")


if __name__ == "__main__":
    main()
