#!/usr/bin/env python3
"""Turn projected PAF pieces (impg or paf_project.py) into per-gene, per-haplotype hits.

Input PAF lines must carry rn:Z: (paf_project.py) or an:Z: (impg) = <LOCUS>|gene or <LOCUS>|padded.

Usage:
  summarize_projection.py --paf proj.paf --genes batch.genes.tsv --index assemblies.csv \
      --outdir genes/ --regions-out work/batch/regions.tsv [--flank 5000]

For each locus x haplotype, gene pieces on the same contig+strand are chained into hits
(a new hit starts when the next piece is > --merge-dist away on the contig or re-covers
>50% of CHM13 bases already covered by the current hit, i.e. a tandem copy).
Per-hit: gene_cov_frac (union of CHM13 gene bases covered / gene length), n_blocks,
identity (matches / (matches + mismatches + gap opens), over pieces).

Hit classes:   primary (best cov, then identity) | copy (cov >= --copy-cov) | fragment
Haplotype status flags (comma-separated, never dropped):
  ok | missing | partial (primary cov < --partial) | split_contigs (pieces on >1 contig needed)
  | split_blocks (primary hit made of >1 alignment block) | multi_copy (>1 hit with cov >= --copy-cov)
  | no_fasta (haplotype in PAF but absent from the assembly index)
Extracted sequences: primary + copies; when the primary is partial, fragments are also extracted
(so a gene split across contigs is fully represented).
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

# s3://human-pangenomics/... paths in the assembly index are read over anonymous HTTPS (config/sources.env)
S3 = os.environ.get("HPRC_BUCKET", "https://human-pangenomics.s3.amazonaws.com").rstrip("/") + "/"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def s3_to_https(u):
    return S3 + u[len("s3://human-pangenomics/"):] if u.startswith("s3://human-pangenomics/") else u


def load_index(path):
    hap = {}
    for r in csv.DictReader(open(path)):
        key = f"{r['sample_id']}#{r['haplotype']}"
        hap[key] = dict(url=s3_to_https(r["assembly"]), name=r["assembly_name"])
    return hap


def union_len(ivs):
    tot = 0; cur_s = cur_e = None
    for s, e in sorted(ivs):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def parse_piece(f):
    tags = {t[:2]: t[5:] for t in f[12:]}
    m = x = gaps = 0
    cg = tags.get("cg")
    if cg:
        for n, op in CIG.findall(cg):
            n = int(n)
            if op in "=M":
                m += n
            elif op == "X":
                x += n
            elif op in "ID":
                gaps += 1
    else:  # no CIGAR: fall back to PAF matches / block length
        m = int(f[9]); x = int(f[10]) - m
    return dict(qname=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                ts=int(f[7]), te=int(f[8]), m=m, x=x, gaps=gaps, rn=tags.get("rn") or tags.get("an", ""))


def chain_hits(pieces, merge_dist):
    """pieces: gene pieces of one haplotype. Returns list of hits."""
    groups = defaultdict(list)
    for p in pieces:
        groups[(p["qname"], p["strand"])].append(p)
    hits = []
    for (qname, strand), ps in groups.items():
        ps.sort(key=lambda p: p["qs"])
        cur = None
        for p in ps:
            if cur is not None:
                covered = union_len(cur["tiv"])
                newcov = union_len(cur["tiv"] + [(p["ts"], p["te"])]) - covered
                redundant = (p["te"] - p["ts"]) > 0 and newcov < 0.5 * (p["te"] - p["ts"])
                if p["qs"] - cur["qe"] > merge_dist or redundant:
                    hits.append(cur); cur = None
            if cur is None:
                cur = dict(qname=qname, qlen=p["qlen"], strand=strand, qs=p["qs"], qe=p["qe"],
                           tiv=[], m=0, x=0, gaps=0, n=0)
            cur["qs"] = min(cur["qs"], p["qs"]); cur["qe"] = max(cur["qe"], p["qe"])
            cur["tiv"].append((p["ts"], p["te"]))
            cur["m"] += p["m"]; cur["x"] += p["x"]; cur["gaps"] += p["gaps"]; cur["n"] += 1
        if cur is not None:
            hits.append(cur)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--regions-out", required=True)
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--merge-dist", type=int, default=20000,
                    help="max contig gap between pieces chained into one hit (default 20 kb)")
    ap.add_argument("--copy-cov", type=float, default=0.5)
    ap.add_argument("--partial", type=float, default=0.95)
    a = ap.parse_args()

    index = load_index(a.index)
    genes = {r["locus"]: r for r in csv.DictReader(open(a.genes), delimiter="\t")}
    for g in genes.values():
        for k in ("start", "end", "pad_start", "pad_end"):
            g[k] = int(g[k])

    # pieces[locus][kind][haplotype] -> list
    pieces = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    paf_haps = set()
    with open(a.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            p = parse_piece(f)
            locus, _, kind = p["rn"].rpartition("|")
            if locus not in genes:
                continue
            sample, hp = p["qname"].split("#")[:2]
            h = f"{sample}#{hp}"
            if h.startswith("CHM13#"):
                continue  # self
            paf_haps.add(h)
            pieces[locus][kind][h].append(p)

    # haplotypes attempted = everything in the assembly index except CHM13 (handled as control) + PAF-only ones
    all_haps = sorted((set(index) - {"CHM13#0"}) | paf_haps)
    reg_rows = []
    os.makedirs(a.outdir, exist_ok=True)
    for locus, g in genes.items():
        glen = g["end"] - g["start"]
        gstrand = g["strand"]
        gdir = os.path.join(a.outdir, locus)
        os.makedirs(gdir, exist_ok=True)
        rows = []
        # CHM13 control: the annotation interval itself
        chm = dict(qname=g["chrom"], qs=g["start"], qe=g["end"], strand="+", n=1, cov=1.0, ident=1.0,
                   ps=g["pad_start"], pe=g["pad_end"], cls="primary", extract=True)
        hap_hits = {"CHM13#0": ([chm], "control")}
        for h in all_haps:
            hits = chain_hits(pieces[locus]["gene"].get(h, []), a.merge_dist)
            for t in hits:
                t["cov"] = union_len(t["tiv"]) / glen
                denom = t["m"] + t["x"] + t["gaps"]
                t["ident"] = t["m"] / denom if denom else 0.0
            hits.sort(key=lambda t: (-t["cov"], -t["ident"]))
            flags = []
            if not hits:
                flags.append("missing")
            else:
                prim = hits[0]
                big = [t for t in hits if t["cov"] >= a.copy_cov]
                for i, t in enumerate(hits):
                    t["cls"] = "primary" if i == 0 else ("copy" if t["cov"] >= a.copy_cov else "fragment")
                    t["extract"] = t["cls"] != "fragment" or prim["cov"] < a.partial
                if prim["cov"] < a.partial:
                    flags.append("partial")
                    tot_cov = union_len([iv for t in hits for iv in t["tiv"]]) / glen
                    if len({t["qname"] for t in hits}) > 1 and tot_cov > prim["cov"] + 0.01:
                        flags.append("split_contigs")
                if prim["n"] > 1:
                    flags.append("split_blocks")
                if len(big) > 1:
                    flags.append("multi_copy")
                # padded interval: gene hit extended by projected flank pieces on the same contig/strand
                pads = pieces[locus]["padded"].get(h, [])
                for t in hits:
                    w = 2 * a.flank + 1000
                    near = [p for p in pads if p["qname"] == t["qname"] and p["strand"] == t["strand"]
                            and p["qe"] > t["qs"] - w and p["qs"] < t["qe"] + w]
                    t["ps"] = min([t["qs"]] + [max(p["qs"], t["qs"] - w) for p in near])
                    t["pe"] = max([t["qe"]] + [min(p["qe"], t["qe"] + w) for p in near])
            if h not in index:
                flags.append("no_fasta")
            hap_hits[h] = (hits, ",".join(flags) if flags else "ok")

        for h, (hits, status) in hap_hits.items():
            sample, hp = h.split("#")
            url = index.get(h, {}).get("url", "NA")
            if not hits:
                rows.append([locus, sample, hp, "NA", "NA", "NA", "NA", "0", "0", "NA", "NA", "NA", "NA",
                             "none", status, "no", "NA"])
                continue
            for rank, t in enumerate(hits, 1):
                rc = (t["strand"] == "-") != (gstrand == "-")
                sense = "-" if rc else "+"
                hdr = f"{locus}|{t['qname']}:{t['qs'] + 1}-{t['qe']}({sense})"
                phdr = f"{locus}|{t['qname']}:{t['ps'] + 1}-{t['pe']}({sense})"
                if rank > 1:
                    hdr += f" {t['cls']}{rank}"; phdr += f" {t['cls']}{rank}"
                do = t["extract"] and url != "NA"
                rows.append([locus, sample, hp, t["qname"], str(t["qs"]), str(t["qe"]), sense,
                             f"{t['cov']:.4f}", str(t["n"]), f"{t['ident']:.5f}", str(t["ps"]), str(t["pe"]),
                             str(t["qe"] - t["qs"]), t["cls"], status, "yes" if do else "no", hdr])
                if do:
                    for kind, s, e, hh in (("gene", t["qs"], t["qe"], hdr), ("padded", t["ps"], t["pe"], phdr)):
                        reg_rows.append([h, url, f"{t['qname']}:{s + 1}-{e}", "1" if rc else "0", locus, kind, hh])

        with open(os.path.join(gdir, f"{locus}.coords.tsv"), "w") as out:
            out.write("gene\tsample\thap\tcontig\tstart\tend\tstrand\tgene_cov_frac\tn_blocks\tidentity"
                      "\tpadded_start\tpadded_end\tlength\thit_class\thap_status\textracted\tfasta_header\n")
            for r in rows:
                out.write("\t".join(r) + "\n")

    with open(a.regions_out, "w") as out:
        out.write("haplotype\turl\tregion\trc\tlocus\tkind\theader\n")
        for r in reg_rows:
            out.write("\t".join(r) + "\n")
    sys.stderr.write(f"{len(genes)} loci, {len(all_haps)} haplotypes, {len(reg_rows)} regions to fetch\n")


if __name__ == "__main__":
    main()
