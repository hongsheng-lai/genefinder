#!/usr/bin/env python3
"""Look up genes in a CHM13 v2.0 GFF3 annotation and write CHM13 BED files.

Usage:
  gene2bed.py [GENE ...] [--gene-list FILE] [--annotation gencode|refseq|PATH]
              [--flank 5000] [--biotype protein_coding] --out PREFIX

"ALL" (as a gene name) selects every feature of type "gene" (optionally
filtered by --biotype). Genes are matched against gene_name, Name, gene,
gene_id, ID, source_gene (with or without .version) and source_gene_common_name.

Writes:
  PREFIX.gene.bed    BED6, gene body (chrom = CHM13#0#chrN, 0-based start)
  PREFIX.padded.bed  BED6, gene body +/- flank (clipped to chromosome)
  PREFIX.genes.tsv   one row per selected locus with metadata
  PREFIX.lookup.txt  report of not-found and multi-hit queries
Exit code 0 even if some genes are missing (they are reported).
"""
import argparse
import gzip
import os
import re
import sys
import urllib.parse
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA = os.environ.get("GENEFINDER_DATA", os.path.join(os.path.dirname(HERE), "data"))
ANNOTATIONS = {  # file names as in config/sources.env
    "gencode": os.environ.get("GENCODE_GFF3", "chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz"),
    "refseq": os.environ.get("REFSEQ_GFF3", "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"),
}
ANN_URL = os.environ.get("ANNOTATION_BASE_URL",
                         "https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/annotation")
KEYS = ("gene_name", "Name", "gene", "gene_id", "ID", "source_gene", "source_gene_common_name")


def parse_attrs(s):
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = urllib.parse.unquote(v)
    return d


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def load_chrom_sizes(data_dir):
    sizes = {}
    fai = os.path.join(data_dir, "reference", "chm13.fa.gz.fai")
    if os.path.exists(fai):
        for line in open(fai):
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="*")
    ap.add_argument("--gene-list")
    ap.add_argument("--annotation", default="gencode")
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--biotype", help="keep only genes with this biotype/gene_biotype (e.g. protein_coding)")
    ap.add_argument("--chrom-prefix", default="CHM13#0#")
    ap.add_argument("--data-dir", default=DEFAULT_DATA, help="directory populated by run.sh --setup")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    queries = list(a.genes)
    if a.gene_list:
        queries += [l.split()[0] for l in open(a.gene_list) if l.strip() and not l.startswith("#")]
    if not queries:
        sys.exit("no genes given")
    want_all = any(q.upper() == "ALL" for q in queries)
    queries = [q for q in dict.fromkeys(queries) if q.upper() != "ALL"]
    qset = set(queries)

    ann = ANNOTATIONS.get(a.annotation, a.annotation)
    if not os.path.exists(ann):
        ann = os.path.join(a.data_dir, "annotation", ann)
    if not os.path.exists(ann):
        sys.exit(f"annotation not found: {ann} (run ./run.sh --setup, or download from {ANN_URL}/)")

    hits = defaultdict(list)  # query -> list of records
    n_genes = 0
    with gzip.open(ann, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            n_genes += 1
            at = parse_attrs(f[8])
            biotype = at.get("biotype") or at.get("gene_biotype") or "NA"
            if a.biotype and biotype != a.biotype:
                continue
            rec = dict(chrom=f[0], start=int(f[3]) - 1, end=int(f[4]), strand=f[6],
                       gene_name=at.get("gene_name") or at.get("Name") or at.get("gene") or at.get("ID"),
                       gene_id=at.get("gene_id") or at.get("ID") or "NA",
                       source_gene=at.get("source_gene") or at.get("db_xref") or "NA",
                       biotype=biotype)
            if want_all:
                hits[rec["gene_name"]].append(rec)
                continue
            vals = set()
            for k in KEYS:
                v = at.get(k)
                if v:
                    vals.add(v)
                    if k == "source_gene":
                        vals.add(v.split(".")[0])
            for v in vals & qset:
                hits[v].append(rec)

    sizes = load_chrom_sizes(a.data_dir)
    missing = [q for q in queries if q not in hits]
    multi = {q: r for q, r in hits.items() if len(r) > 1}

    with open(a.out + ".gene.bed", "w") as gb, open(a.out + ".padded.bed", "w") as pb, \
            open(a.out + ".genes.tsv", "w") as gt:
        gt.write("locus\tquery\tgene_name\tgene_id\tsource_gene\tbiotype\tchrom\tstart\tend\tstrand"
                 "\tpad_start\tpad_end\tn_copies_in_annotation\n")
        order = sorted(hits) if want_all else [q for q in queries if q in hits]
        for q in order:
            recs = sorted(hits[q], key=lambda r: (r["chrom"], r["start"]))
            for i, r in enumerate(recs):
                # unique, filesystem-safe locus name; extra copies get a suffix
                locus = safe(q) if i == 0 else f"{safe(q)}_copy{i + 1}"
                chrom = a.chrom_prefix + r["chrom"]
                ps = max(0, r["start"] - a.flank)
                pe = r["end"] + a.flank
                if chrom in sizes:
                    pe = min(pe, sizes[chrom])
                gb.write(f"{chrom}\t{r['start']}\t{r['end']}\t{locus}\t0\t{r['strand']}\n")
                pb.write(f"{chrom}\t{ps}\t{pe}\t{locus}\t0\t{r['strand']}\n")
                gt.write(f"{locus}\t{q}\t{r['gene_name']}\t{r['gene_id']}\t{r['source_gene']}\t{r['biotype']}"
                         f"\t{chrom}\t{r['start']}\t{r['end']}\t{r['strand']}\t{ps}\t{pe}\t{len(recs)}\n")

    with open(a.out + ".lookup.txt", "w") as lk:
        lk.write(f"annotation\t{ann}\ngene_features_scanned\t{n_genes}\nflank\t{a.flank}\n")
        lk.write(f"queries\t{'ALL' if want_all else len(queries)}\nfound\t{len(hits)}\n")
        lk.write(f"not_found\t{len(missing)}\t{','.join(missing)}\n")
        lk.write(f"multi_hit\t{len(multi)}\n")
        for q, recs in sorted(multi.items()):
            lk.write("  " + q + "\t" + "; ".join(f"{r['chrom']}:{r['start']}-{r['end']}({r['strand']}) {r['gene_id']}"
                                                  for r in recs) + "\n")
    sys.stderr.write(open(a.out + ".lookup.txt").read() if len(multi) < 50 else
                     f"found {len(hits)}, not found {len(missing)}, multi-hit {len(multi)}\n")


if __name__ == "__main__":
    main()
