#!/usr/bin/env python3
"""Write genes/<LOCUS>/<LOCUS>.summary.txt from coords.tsv + gene.fa, and append a row to a batch table.

Usage: gene_summary.py --outdir genes/ --genes batch.genes.tsv [--table batch_summary.tsv]
"""
import argparse
import csv
import os
import statistics
from collections import defaultdict

FLAGS = ("missing", "partial", "split_contigs", "split_blocks", "multi_copy", "no_fasta")
CONTROLS = ("CHM13#0", "GRCh38#0")


def fasta_lengths(path):
    lens = {}
    if not os.path.exists(path):
        return lens
    name = None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip(); lens[name] = 0
        else:
            lens[name] += len(line.strip())
    return lens


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(q * (len(xs) - 1))))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--table")
    a = ap.parse_args()
    table_rows = []
    for g in csv.DictReader(open(a.genes), delimiter="\t"):
        locus = g["locus"]
        d = os.path.join(a.outdir, locus)
        rows = list(csv.DictReader(open(os.path.join(d, f"{locus}.coords.tsv")), delimiter="\t"))
        lens = fasta_lengths(os.path.join(d, f"{locus}.gene.fa"))
        plens = fasta_lengths(os.path.join(d, f"{locus}.padded.fa"))
        status = {}
        extracted = defaultdict(int)
        prim_len = {}
        for r in rows:
            h = f"{r['sample']}#{r['hap']}"
            status[h] = r["hap_status"]
            if r["extracted"] == "yes" and r["fasta_header"] in lens:
                extracted[h] += 1
                if r["hit_class"] == "primary":
                    prim_len[h] = lens[r["fasta_header"]]
        hprc = [h for h in status if h not in CONTROLS]
        counts = {f: sum(1 for h in hprc if f in status[h].split(",")) for f in FLAGS}
        ok = sum(1 for h in hprc if status[h] == "ok")
        n_ext = sum(1 for h in hprc if extracted[h])
        L = [prim_len[h] for h in hprc if h in prim_len]
        glen = int(g["end"]) - int(g["start"])
        out = []
        out.append(f"gene\t{locus}\t{g['gene_name']}\t{g['gene_id']}\t{g['source_gene']}\t{g['biotype']}")
        out.append(f"chm13_locus\t{g['chrom']}:{int(g['start']) + 1}-{g['end']}({g['strand']})\tlength={glen}"
                   f"\tpadded={g['pad_start']}-{g['pad_end']}")
        out.append(f"haplotypes_attempted\t{len(hprc)}\t(HPRC R2 haplotypes; controls CHM13 and GRCh38 reported separately)")
        out.append(f"haplotypes_extracted\t{n_ext}")
        out.append(f"haplotypes_clean_ok\t{ok}\t(single contig, single block, >=95% covered, single copy)")
        for f in FLAGS:
            out.append(f"{f}\t{counts[f]}")
        out.append(f"gene_sequences_written\t{len(lens)}\tpadded_sequences_written\t{len(plens)}")
        if L:
            out.append("primary_length_distribution\tn={}\tmin={}\tp5={}\tmedian={}\tmean={:.0f}\tp95={}\tmax={}".format(
                len(L), min(L), pct(L, .05), int(statistics.median(L)), statistics.mean(L), pct(L, .95), max(L)))
            near = sum(1 for x in L if abs(x - glen) <= 0.05 * glen)
            out.append(f"primary_within_5pct_of_CHM13\t{near}/{len(L)}")
            # coarse histogram, 2% bins relative to CHM13 length
            bins = defaultdict(int)
            for x in L:
                bins[round((x / glen - 1) * 50) * 2] += 1
            out.append("length_hist_pct_vs_CHM13\t" + " ".join(f"{k:+d}%:{v}" for k, v in sorted(bins.items())))
        for c in CONTROLS:
            if c in status:
                out.append(f"control_{c}\t{status[c]}\tlength={prim_len.get(c, 'NA')}")
        flagged = sorted(h for h in hprc if status[h] != "ok")
        out.append(f"flagged_haplotypes\t{len(flagged)}")
        for h in flagged:
            hr = [r for r in rows if f"{r['sample']}#{r['hap']}" == h]
            desc = "; ".join(f"{r['hit_class']} {r['contig']}:{r['start']}-{r['end']}({r['strand']}) "
                             f"cov={r['gene_cov_frac']} blocks={r['n_blocks']} id={r['identity']}"
                             for r in hr if r["contig"] != "NA")
            out.append(f"  {h}\t{status[h]}\t{desc or '-'}")
        open(os.path.join(d, f"{locus}.summary.txt"), "w").write("\n".join(out) + "\n")
        table_rows.append([locus, g["gene_name"], str(glen), str(len(hprc)), str(n_ext), str(ok)]
                          + [str(counts[f]) for f in FLAGS]
                          + [str(int(statistics.median(L))) if L else "NA"])
    if a.table:
        new = not os.path.exists(a.table)
        with open(a.table, "a") as fh:
            if new:
                fh.write("locus\tgene_name\tchm13_len\tattempted\textracted\tok\t" + "\t".join(FLAGS)
                         + "\tmedian_len\n")
            for r in table_rows:
                fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
