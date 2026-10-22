#!/usr/bin/env python3
"""Project CHM13 BED intervals onto haplotypes by fetching ONLY the overlapping PAF records from the
remote (BGZF) PAF with HTTP range requests, located through the offset index from build_paf_index.py.

Usage: remote_project.py --bed query.bed --index paf_offsets.bed.gz --paf PAF_URL_OR_PATH [-j 16] > projected.paf
Output is the same PAF layout as paf_project.py / `impg query -o paf` (rn:Z = BED name).
For one gene this moves a few MB instead of the whole 5.6 GB alignment file.
"""
import argparse
import os
import subprocess
import sys
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_project import load_bed, project_line  # noqa: E402

MERGE_GAP = 32768  # merge byte ranges closer than this into one request


def fetch_range(src, a, b, retries=4):
    """Bytes [a, b) of a URL (HTTP Range) or local file."""
    if not src.startswith(("http://", "https://")):
        with open(src, "rb") as fh:
            fh.seek(a)
            return fh.read(b - a)
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(src, headers={"Range": f"bytes={a}-{b - 1}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) == b - a:
                return data
            raise IOError(f"short read {len(data)} != {b - a}")
        except Exception as e:  # noqa: BLE001 - retry any transport error
            if attempt == retries:
                raise
            sys.stderr.write(f"retry {attempt} for bytes {a}-{b}: {e}\n")
            time.sleep(2 * attempt)


def inflate_blocks(raw, base):
    """Decompress concatenated BGZF blocks; return (uncompressed bytes, {compressed_offset: position})."""
    pos, out, where = 0, [], {}
    total = 0
    while pos < len(raw):
        bsize = int.from_bytes(raw[pos + 16:pos + 18], "little") + 1
        where[base + pos] = total
        data = zlib.decompress(raw[pos + 18:pos + bsize - 8], -15)
        out.append(data)
        total += len(data)
        pos += bsize
    return b"".join(out), where


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--index", required=True, help="paf_offsets.bed.gz (+ .tbi), local path or URL")
    ap.add_argument("--paf", required=True, help="the BGZF PAF the index was built from (URL or local path)")
    ap.add_argument("-j", "--jobs", type=int, default=16)
    a = ap.parse_args()

    # 1. records overlapping any BED interval (tabix handles many regions in one call)
    p = subprocess.run(["tabix", "-R", a.bed, a.index], capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"tabix failed: {p.stderr.strip()}")
    recs = sorted({tuple(int(x) for x in l.split("\t")[3:7]) for l in p.stdout.splitlines() if l})

    # 2. merge into byte ranges and fetch them in parallel
    ranges = []
    for cs, ce, us, ul in recs:
        if ranges and cs <= ranges[-1][1] + MERGE_GAP:
            ranges[-1][1] = max(ranges[-1][1], ce)
            ranges[-1][2].append((cs, us, ul))
        else:
            ranges.append([cs, ce, [(cs, us, ul)]])
    t0 = time.time()

    def get(rg):
        cs, ce, members = rg
        buf, where = inflate_blocks(fetch_range(a.paf, cs, ce), cs)
        return [buf[where[c] + u: where[c] + u + n].decode() for c, u, n in members]

    idx = load_bed(a.bed)
    out = sys.stdout
    n_out = 0
    with ThreadPoolExecutor(a.jobs) as ex:
        for lines in ex.map(get, ranges):
            for line in lines:
                for o in project_line(line, idx):
                    out.write(o)
                    n_out += 1
    nbytes = sum(ce - cs for cs, ce, _ in ranges)
    sys.stderr.write(f"remote projection: {len(recs)} PAF records in {len(ranges)} range requests, "
                     f"{nbytes / 1e6:.1f} MB fetched in {time.time() - t0:.1f}s -> {n_out} projected pieces\n")


if __name__ == "__main__":
    main()
