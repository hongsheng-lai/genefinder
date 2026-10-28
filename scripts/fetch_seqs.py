#!/usr/bin/env python3
"""Fetch projected regions from the remote HPRC assemblies and write per-gene FASTAs.

Usage: fetch_seqs.py --regions regions.tsv --workdir work/batch/fetch --outdir genes/
                     [-j 16] [--retries 3] [--method auto|faidx|stream] [--stream-threshold 3000] [--keep-cache]

One job per (assembly, orientation) gathers every region of the batch for that assembly:
  faidx  : `samtools faidx URL -r regions.txt [-i]` (bgzip + .fai + .gzi, HTTP range requests).
           Costs ~0.2 MB of network and ~0.4 s latency per region, plus ~0.75 MB for the .gzi per job.
  stream : `curl URL | bgzip -dc | stream_extract.py` - one sequential pass over the assembly
           (~0.9 GB compressed), extracting all regions on the fly; the assembly is never stored.
  auto   : stream when an assembly job has more than --stream-threshold regions, else faidx.
           A local assembly (a path rather than a URL, from --fasta-dir) is always read with faidx.
A gene region inside the padded region of the same hit is sliced from it rather than fetched.
Resumable: a job whose <workdir>/<hap>.<orient>.fa exists is skipped (outputs written atomically).
After all jobs succeed and the per-gene FASTAs are written, the per-assembly cache is removed
(unless --keep-cache).
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
CONTROLS = ("CHM13#0", "GRCh38#0")


def parse_region(reg):
    c, se = reg.rsplit(":", 1)
    s, e = se.split("-")
    return c, int(s), int(e)


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def fetch_one(hap, url, orient, regions, workdir, retries, method):
    tag = hap.replace("#", "_")
    out = os.path.join(workdir, f"{tag}.{orient}.fa")
    if os.path.exists(out):
        return hap, orient, "cached", 0, 0.0, method
    rfile = os.path.join(workdir, f"{tag}.{orient}.regions.txt")
    with open(rfile, "w") as fh:
        fh.write("\n".join(regions) + "\n")
    if method == "faidx":
        cmd = ["samtools", "faidx", url, "-r", rfile, "-o", out + ".tmp"]
        if orient == "minus":
            cmd[2:2] = ["-i", "--mark-strand", "no"]
        shell = False
    else:  # orient == "both": regions carry a "-" column when they need reverse-complementing
        reader = f"curl -sS --retry 5 '{url}'" if "://" in url else f"cat '{url}'"
        cmd = (f"set -o pipefail; {reader} | bgzip -dc | "
               f"python3 '{HERE}/stream_extract.py' --regions '{rfile}' > '{out}.tmp'")
        shell = True
    err = [""]
    for attempt in range(1, retries + 1):
        t0 = time.time()
        # cwd=workdir so any index files htslib might cache stay out of the way
        p = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, shell=shell,
                           executable="/bin/bash" if shell else None)
        if p.returncode == 0:
            os.replace(out + ".tmp", out)
            return hap, orient, "ok", attempt, time.time() - t0, method
        err = p.stderr.strip().splitlines()[-1:] or [f"exit {p.returncode}"]
        time.sleep(5 * attempt)
    return hap, orient, f"FAILED: {err[0]}", retries, 0.0, method


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("-j", "--jobs", type=int, default=16)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--method", choices=("auto", "faidx", "stream"), default="auto")
    ap.add_argument("--stream-threshold", type=int, default=3000,
                    help="auto: stream an assembly when its job has more regions than this (default 3000)")
    ap.add_argument("--keep-cache", action="store_true")
    a = ap.parse_args()
    a.workdir = os.path.abspath(a.workdir)
    os.makedirs(a.workdir, exist_ok=True)

    rows = list(csv.DictReader(open(a.regions), delimiter="\t"))
    # A gene region contained in a padded region of the same hit is not fetched separately:
    # it is sliced out of the padded sequence (halves the number of range requests).
    padded = defaultdict(list)  # (hap, rc, locus, contig) -> [(start, end, region)]
    for r in rows:
        if r["kind"] == "padded":
            c, s, e = parse_region(r["region"])
            padded[(r["haplotype"], r["rc"], r["locus"], c)].append((s, e, r["region"]))
    for r in rows:
        r["src"], r["off"] = r["region"], None
        if r["kind"] == "gene":
            c, s, e = parse_region(r["region"])
            for ps, pe, preg in padded[(r["haplotype"], r["rc"], r["locus"], c)]:
                if ps <= s and e <= pe:
                    # slice offsets inside the fetched (possibly reverse-complemented) padded sequence
                    r["src"] = preg
                    r["off"] = (pe - e, pe - s + 1) if r["rc"] == "1" else (s - ps, e - ps + 1)
                    break
        r["orient"] = "minus" if r["rc"] == "1" else "plus"
    # method per assembly; a streamed assembly is read once for both orientations ("both" job)
    per_hap = defaultdict(set)
    for r in rows:
        per_hap[r["haplotype"]].add((r["src"], r["orient"]))
    hap_url = {r["haplotype"]: r["url"] for r in rows}
    hap_method = {h: (a.method if a.method != "auto" else
                      ("stream" if len(v) > a.stream_threshold and "://" in hap_url[h] else "faidx"))
                  for h, v in per_hap.items()}
    jobs = defaultdict(set)  # (hap, url, orient) -> regions to fetch
    for r in rows:
        if hap_method[r["haplotype"]] == "stream":
            r["job"] = "both"
            r["key"] = r["src"] + ("/rc" if r["orient"] == "minus" else "")
            jobs[(r["haplotype"], r["url"], "both")].add(r["src"] + ("\t-" if r["orient"] == "minus" else ""))
        else:
            r["job"] = r["orient"]
            r["key"] = r["src"]
            jobs[(r["haplotype"], r["url"], r["orient"])].add(r["src"])
    methods = {k: hap_method[k[0]] for k in jobs}
    sys.stderr.write(f"{len(rows)} output regions -> {sum(len(v) for v in jobs.values())} fetched regions in "
                     f"{len(jobs)} assembly jobs ({sum(m == 'stream' for m in methods.values())} streamed)\n")

    log = open(os.path.join(a.workdir, "fetch.log"), "a")
    t_fetch = time.time()
    failed = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(fetch_one, h, u, o, sorted(rs), a.workdir, a.retries, methods[(h, u, o)])
                for (h, u, o), rs in jobs.items()]
        for i, f in enumerate(as_completed(futs), 1):
            hap, orient, st, att, dt, meth = f.result()
            log.write(f"{hap}\t{orient}\t{st}\tattempts={att}\t{dt:.1f}s\t{meth}\n"); log.flush()
            if st.startswith("FAILED"):
                failed.append((hap, orient, st))
            if i % 50 == 0:
                sys.stderr.write(f"  fetched {i}/{len(futs)}\n")
    log.close()
    t_write = time.time()

    # Write per-gene FASTAs one assembly job at a time (bounded memory): truncate, then append
    # in haplotype order with the CHM13 and GRCh38 controls first.
    loci = sorted({r["locus"] for r in rows})
    for locus in loci:
        os.makedirs(os.path.join(a.outdir, locus), exist_ok=True)
        for kind in ("gene", "padded"):
            open(os.path.join(a.outdir, locus, f"{locus}.{kind}.fa"), "w").close()
    by_hap = defaultdict(list)  # rows per haplotype, in regions.tsv order (same output for faidx/stream)
    for r in rows:
        by_hap[r["haplotype"]].append(r)
    order = sorted(by_hap, key=lambda h: (CONTROLS.index(h) if h in CONTROLS else len(CONTROLS), h))
    missing = 0
    buf, buf_bytes = defaultdict(list), 0

    def flush():
        for (locus, kind), recs in buf.items():
            with open(os.path.join(a.outdir, locus, f"{locus}.{kind}.fa"), "a") as fh:
                fh.writelines(recs)
        buf.clear()

    for h in order:
        seqs = {}
        for job in {r["job"] for r in by_hap[h]}:
            path = os.path.join(a.workdir, f"{h.replace('#', '_')}.{job}.fa")
            if os.path.exists(path):
                seqs[job] = dict(read_fasta(path))
        for r in by_hap[h]:
            s = seqs.get(r["job"], {}).get(r["key"])
            if s is None:
                missing += 1
                continue
            if r["off"]:
                s = s[r["off"][0]:r["off"][1]]
            rec = f">{r['header']}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60))
            buf[(r["locus"], r["kind"])].append(rec)
            buf_bytes += len(rec)
        if buf_bytes > 512_000_000:
            flush(); buf_bytes = 0
    flush()
    sys.stderr.write(f"{len(jobs)} assembly jobs, {len(failed)} failed, {missing} regions without sequence; "
                     f"fetch phase {t_write - t_fetch:.0f}s, write phase {time.time() - t_write:.0f}s\n")
    for x in failed:
        sys.stderr.write("  FAILED " + "\t".join(x) + "\n")
    if not failed and not missing and not a.keep_cache:
        for p in glob.glob(os.path.join(a.workdir, "*.fa")) + glob.glob(os.path.join(a.workdir, "*.regions.txt")):
            os.remove(p)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
