#!/usr/bin/env python3
"""Collect GNU `time -v` logs (time.<step>.txt) in a directory into a TSV.

Usage: resources.py LOGDIR [LOGDIR ...]   (prints TSV to stdout)
"""
import glob
import os
import re
import sys


def parse(path):
    d = {}
    for line in open(path):
        line = line.strip()
        if ":" not in line:
            continue
        k, v = line.rsplit(":", 1) if line.startswith("Elapsed") else line.split(":", 1)
        d[k.strip()] = v.strip()
    el = [l for l in open(path) if "Elapsed" in l]
    if el:
        t = el[0].split("):", 1)[1].strip().split(":")
        d["wall_s"] = sum(float(x) * 60 ** i for i, x in enumerate(reversed(t)))
    return d


if __name__ == "__main__":
    print("logdir\tstep\twall_s\tuser_s\tsys_s\tcpu_pct\tmax_rss_mb\tfs_out_mb\twork_dir_mb\texit")
    for logdir in sys.argv[1:]:
        files = sorted(glob.glob(os.path.join(logdir, "time.*.txt")), key=os.path.getmtime)
        for p in files:
            d = parse(p)
            step = re.sub(r"^time\.|\.txt$", "", os.path.basename(p))
            wd = d.get("Work dir bytes")
            print("\t".join([
                os.path.basename(logdir.rstrip("/")), step, f"{d.get('wall_s', 0):.1f}",
                d.get("User time (seconds)", "NA"), d.get("System time (seconds)", "NA"),
                d.get("Percent of CPU this job got", "NA"),
                f"{int(d.get('Maximum resident set size (kbytes)', 0)) / 1024:.1f}",
                f"{int(d.get('File system outputs', 0)) * 512 / 1e6:.1f}",
                f"{int(wd) / 1e6:.1f}" if wd else "NA",
                d.get("Exit status", "NA")]))
