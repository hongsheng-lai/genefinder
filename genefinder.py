#!/usr/bin/env python3
"""genefinder: extract one CHM13 gene from every HPRC Release 2 haplotype.

  genefinder.py setup online     one-time: CHM13 gene table, assembly index, PAF offset index (~3-5 min)
  genefinder.py setup offline    one-time: index of gene lines from the per-haplotype Liftoff GFF3s
                                 in the local mirror (~2 min, Salz cluster)
  genefinder.py get GENE         extract GENE -> OUTDIR/GENE.fa + GENE.tsv

Two ways to find the gene in each haplotype:
  online   project the CHM13 gene interval through the HPRC Release 2 alignments (PAF, cg:Z CIGARs).
           Only the PAF records overlapping the gene are fetched (HTTP range requests); sequences are
           read from the bgzipped assemblies on S3 with samtools faidx.
  offline  look the gene up in the per-haplotype Liftoff annotations of the local mirror
           (MIRROR/annotation.liftoff) and read sequences from MIRROR/fasta. No network.
"""
import argparse
import csv
import glob
import gzip
import os
import re
import shutil
import struct
import subprocess
import sys
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BUCKET = "https://human-pangenomics.s3.amazonaws.com"
PAF_URL = os.environ.get("PAF_URL", f"{BUCKET}/pangenomes/freeze/release2/impg/pafs/hprc465vschm13.aln.paf.gz")
ASSEMBLY_INDEX_URL = os.environ.get("ASSEMBLY_INDEX_URL", (
    "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
    "2a7e3a5b1ca10454762699176b1a0a1027ccbf21/data_tables/assemblies_release2_v1.0.index.csv"))
# RefSeq genes on CHM13 v2.0 (Liftoff): the same gene namespace as the mirror's per-haplotype Liftoff files.
CHM13_GFF = os.environ.get("CHM13_GFF", f"{BUCKET}/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
MIRROR = os.environ.get("HPRC_MIRROR", "/ccb/salz8/HPRC/assemblies")

OK_COV, COPY_COV, CHAIN_GAP = 0.95, 0.5, 20000


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def s3_to_https(u):
    return u.replace("s3://human-pangenomics/", BUCKET + "/")


def download(url, dest):
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, 1 << 20)
    os.replace(tmp, dest)


def attrs(col9):
    return dict(kv.split("=", 1) for kv in col9.rstrip(";").split(";") if "=" in kv)


# ---------------------------------------------------------------- setup

def setup_common(data):
    """CHM13 gene table (name, chrom, start0, end, strand), shared by both modes."""
    out = os.path.join(data, "chm13_genes.tsv")
    if os.path.exists(out):
        return
    src = CHM13_GFF
    if "://" in src:
        src = os.path.join(data, os.path.basename(CHM13_GFF))
        log(f"download {CHM13_GFF}")
        download(CHM13_GFF, src)
    with gzip.open(src, "rt") as fh, open(out + ".tmp", "w") as w:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = attrs(f[8].strip())
            name = a.get("gene_name") or a.get("gene") or a.get("Name")
            if name:
                w.write(f"{name}\t{f[0]}\t{int(f[3]) - 1}\t{f[4]}\t{f[6]}\n")
    os.replace(out + ".tmp", out)
    if src != CHM13_GFF:
        os.remove(src)


def bgzf_blocks(fh):
    """Yield (compressed_offset, data, next_offset) for each BGZF block of a stream."""
    coff = 0
    while True:
        head = fh.read(18)
        if not head:
            return
        if len(head) < 18 or head[:4] != b"\x1f\x8b\x08\x04" or head[12:14] != b"BC":
            sys.exit(f"not BGZF at offset {coff}")
        bsize = struct.unpack("<H", head[16:18])[0] + 1
        yield coff, zlib.decompress(fh.read(bsize - 18)[:-8], -15), coff + bsize
        coff += bsize


def build_paf_index(fh, out):
    """One BED-like line per PAF record: target tstart tend cstart cend ustart ulen
    (BGZF block offsets of the record, so it can be range-requested later)."""
    start, prefix, length, n = None, bytearray(), 0, 0
    for coff, data, nxt in bgzf_blocks(fh):
        i = 0
        while data:
            if start is None:
                start = (coff, i)
            j = data.find(b"\n", i)
            seg = len(data) if j < 0 else j
            if len(prefix) < 512:
                prefix += data[i:min(seg, i + 512 - len(prefix))]
            length += seg - i
            if j < 0:
                break
            f = prefix.split(b"\t", 9)
            out.write(f"{f[5].decode()}\t{int(f[7])}\t{int(f[8])}\t{start[0]}\t{nxt}\t{start[1]}\t{length}\n")
            n += 1
            start, prefix, length = None, bytearray(), 0
            i = j + 1
            if i >= len(data):
                break
    return n


def setup_online(data):
    d = os.path.join(data, "online")
    os.makedirs(d, exist_ok=True)
    setup_common(data)
    idx = os.path.join(d, "assemblies.tsv")
    if not os.path.exists(idx):
        with urllib.request.urlopen(ASSEMBLY_INDEX_URL, timeout=120) as r:
            rows = list(csv.DictReader(r.read().decode().splitlines()))
        with open(idx, "w") as w:
            for row in rows:
                w.write(f"{row['sample_id']}#{row['haplotype']}\t{s3_to_https(row['assembly'])}\n")
    off = os.path.join(d, "paf_offsets.bed.gz")
    if not os.path.exists(off + ".tbi"):
        log(f"stream {PAF_URL} once to build its offset index (5.6 GB read, ~120 MB kept)")
        t0 = time.time()
        sort = subprocess.Popen(f"LC_ALL=C sort -k1,1 -k2,2n -S 2G -T '{d}' | bgzip -@4 > '{off}.tmp'",
                                shell=True, stdin=subprocess.PIPE, text=True)
        with urllib.request.urlopen(PAF_URL, timeout=300) as r:
            n = build_paf_index(r, sort.stdin)
        sort.stdin.close()
        if sort.wait():
            sys.exit("sort/bgzip failed")
        os.replace(off + ".tmp", off)
        subprocess.run(["tabix", "-f", "-p", "bed", off], check=True)
        log(f"indexed {n} PAF records in {time.time() - t0:.0f}s")


# Liftoff gene lines -> name, 1, seqid, start(1-based), end, strand, coverage, sequence_ID, extra_copy_number
AWK_GENES = r'''BEGIN{FS=OFS="\t"} $3=="gene"{
  n=""; c="NA"; s="NA"; k="0"; m=split($9,a,";")
  for(i=1;i<=m;i++){p=index(a[i],"="); key=substr(a[i],1,p-1); v=substr(a[i],p+1)
    if(key=="gene_name")n=v; else if(key=="coverage")c=v; else if(key=="sequence_ID")s=v
    else if(key=="extra_copy_number")k=v}
  if(n!="")print n,1,$1,$4,$5,$7,c,s,k}'''


def first_seqname(fai):
    with open(fai) as fh:
        line = fh.readline()
    return line.split("\t")[0] if line else None


def setup_offline(data, jobs):
    d = os.path.join(data, "offline")
    os.makedirs(d, exist_ok=True)
    setup_common(data)
    # 1. assemblies in the mirror, keyed by PanSN sample#hap (file names mix hap1/2 and mat/pat)
    tab = os.path.join(d, "assemblies.tsv")
    with open(tab, "w") as w:
        for fa in sorted(glob.glob(os.path.join(MIRROR, "fasta", "*.fa.gz"))):
            fai = fa + ".fai"
            name = first_seqname(fai) if os.path.getsize(fai) else None
            if name is None:  # unusable .fai in the (read-only) mirror: rebuild a private copy
                fai = os.path.join(d, os.path.basename(fa) + ".fai")
                if not os.path.exists(fai) or not os.path.getsize(fai):
                    log(f"rebuilding empty index {fa}.fai -> {fai}")
                    subprocess.run(["samtools", "faidx", fa, "--fai-idx", fai, "--gzi-idx", fai[:-4] + ".gzi"],
                                   check=True)
                name = first_seqname(fai)
            w.write(f"{'#'.join(name.split('#')[:2])}\t{fa}\t{fai}\n")
    # 2. gene lines of every per-haplotype Liftoff GFF3 -> one name-sorted, tabix-indexed table
    out = os.path.join(d, "liftoff_genes.tsv.gz")
    if os.path.exists(out + ".tbi"):
        return
    gffs = sorted(glob.glob(os.path.join(MIRROR, "annotation.liftoff", "*.gff3.gz")))
    if not gffs:
        sys.exit(f"no Liftoff GFF3 in {MIRROR}/annotation.liftoff")
    parts = os.path.join(d, "parts")
    os.makedirs(parts, exist_ok=True)
    log(f"scan {len(gffs)} Liftoff GFF3 with {jobs} jobs")

    def scan(g):
        p = os.path.join(parts, os.path.basename(g) + ".tsv")
        subprocess.run(f"gzip -dc '{g}' | awk '{AWK_GENES}' > '{p}'", shell=True, check=True)
        return p

    with ThreadPoolExecutor(jobs) as ex:
        files = list(ex.map(scan, gffs))
    with open(os.path.join(d, "annotated.txt"), "w") as w:  # HG00097_hap1.hprc... -> HG00097#1
        for g in gffs:
            m = re.match(r"(.+)_hap(\d)\.", os.path.basename(g))
            w.write(f"{m.group(1)}#{m.group(2)}\n")
    subprocess.run(f"cat '{parts}'/*.tsv | LC_ALL=C sort -k1,1 -S 4G --parallel=8 -T '{d}' "
                   f"| bgzip -@8 > '{out}.tmp'", shell=True, check=True)
    os.replace(out + ".tmp", out)
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", out], check=True)
    for p in files:
        os.remove(p)
    os.rmdir(parts)


# ---------------------------------------------------------------- online: PAF projection

CIG = re.compile(r"(\d+)([MIDN=X])")


def project(cigar, tstart, qs, qe, strand, cs, ce):
    """Project target interval [cs, ce) through one alignment. Returns
    (q_lo, q_hi, t_lo, t_hi, matches, mismatches, gap_opens) or None if no aligned base falls inside."""
    t, off = tstart, 0
    q0 = q1 = t0 = t1 = None
    m = x = gaps = pend = 0  # pend: gaps not yet followed by an aligned base
    for mo in CIG.finditer(cigar):
        n, op = int(mo.group(1)), mo.group(2)
        if t >= ce:
            break
        if op in "M=X":
            a, b = max(t, cs), min(t + n, ce)
            if a < b:
                if q0 is None:
                    q0, t0 = off + a - t, a
                q1, t1 = off + b - t, b
                gaps += pend
                pend = 0
                if op == "X":
                    x += b - a
                else:
                    m += b - a
            t += n
            off += n
        elif op in "DN":
            if q0 is not None and max(t, cs) < min(t + n, ce):
                pend += 1
            t += n
        elif op == "I":
            if q0 is not None and cs < t < ce:
                pend += 1
            off += n
    if q0 is None:
        return None
    if strand == "+":
        return qs + q0, qs + q1, t0, t1, m, x, gaps
    return qe - q1, qe - q0, t0, t1, m, x, gaps


def fetch_range(url, a, b, retries=4):
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b - 1}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) == b - a:
                return data
            raise IOError(f"short read {len(data)} != {b - a}")
        except Exception as e:  # noqa: BLE001
            if attempt == retries:
                raise
            log(f"retry {attempt} bytes {a}-{b}: {e}")
            time.sleep(2 * attempt)


def paf_records(url, index, chrom, start, end, jobs):
    """PAF lines overlapping chrom:[start, end), fetched by byte range."""
    p = subprocess.run(["tabix", index, f"{chrom}:{start + 1}-{end}"], capture_output=True, text=True, check=True)
    recs = sorted({tuple(map(int, l.split("\t")[3:7])) for l in p.stdout.splitlines() if l})
    ranges = []  # merge nearby records into one request
    for cs, ce, us, ul in recs:
        if ranges and cs <= ranges[-1][1] + 32768:
            ranges[-1][1] = max(ranges[-1][1], ce)
            ranges[-1][2].append((cs, us, ul))
        else:
            ranges.append([cs, ce, [(cs, us, ul)]])

    def get(rg):
        cs, ce, members = rg
        raw, pos, buf, where = fetch_range(url, cs, ce), 0, [], {}
        total = 0
        while pos < len(raw):
            bsize = int.from_bytes(raw[pos + 16:pos + 18], "little") + 1
            where[cs + pos] = total
            blk = zlib.decompress(raw[pos + 18:pos + bsize - 8], -15)
            buf.append(blk)
            total += len(blk)
            pos += bsize
        buf = b"".join(buf)
        return [buf[where[c] + u: where[c] + u + n].decode() for c, u, n in members]

    with ThreadPoolExecutor(jobs) as ex:
        lines = [l for chunk in ex.map(get, ranges) for l in chunk]
    log(f"  {len(recs)} PAF records, {sum(r[1] - r[0] for r in ranges) / 1e6:.1f} MB in {len(ranges)} requests")
    return lines


def chain(pieces, gene_len):
    """Group projected pieces of one haplotype into hits: same contig and strand, < CHAIN_GAP apart on the
    contig, and not re-covering CHM13 bases already in the hit (that would be a tandem copy)."""
    hits = []
    for key in sorted({(p[0], p[3]) for p in pieces}):
        cur = None
        for p in sorted((p for p in pieces if (p[0], p[3]) == key), key=lambda p: p[1]):
            contig, q_lo, q_hi, strand, t_lo, t_hi, m, x, g = p
            if cur and q_lo - cur["end"] <= CHAIN_GAP and all(t_hi <= a or t_lo >= b for a, b in cur["t"]):
                cur["end"] = max(cur["end"], q_hi)
                cur["t"].append((t_lo, t_hi))
                cur["mxg"] = [u + v for u, v in zip(cur["mxg"], (m, x, g))]
            else:
                cur = dict(contig=contig, start=q_lo, end=q_hi, aln_strand=strand, t=[(t_lo, t_hi)], mxg=[m, x, g])
                hits.append(cur)
    out = []
    for h in hits:
        m, x, g = h["mxg"]
        out.append(dict(contig=h["contig"], start=h["start"], end=h["end"], aln_strand=h["aln_strand"],
                        cov=sum(b - a for a, b in h["t"]) / gene_len, ident=m / (m + x + g) if m + x + g else 0.0))
    return out


def find_online(gene, data, jobs):
    chrom, gs, ge, gstrand = gene["chrom"], gene["start"], gene["end"], gene["strand"]
    lines = paf_records(PAF_URL, os.path.join(data, "online", "paf_offsets.bed.gz"), f"CHM13#0#{chrom}", gs, ge, jobs)
    pieces = {}
    for line in lines:
        f = line.split("\t")
        ts, te = int(f[7]), int(f[8])
        cg = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), None)
        if cg is None or te <= gs or ts >= ge:
            continue
        r = project(cg.strip(), ts, int(f[2]), int(f[3]), f[4], max(gs, ts), min(ge, te))
        if r:
            hap = "#".join(f[0].split("#")[:2])
            pieces.setdefault(hap, []).append((f[0], r[0], r[1], f[4]) + r[2:])
    hits = {}
    for hap, ps in pieces.items():
        hits[hap] = chain(ps, ge - gs)
        for h in hits[hap]:  # gene orientation on the contig
            h["strand"] = gstrand if h.pop("aln_strand") == "+" else flip(gstrand)
    return hits


# ---------------------------------------------------------------- offline: Liftoff lookup

def find_offline(name, data):
    p = subprocess.run(["tabix", os.path.join(data, "offline", "liftoff_genes.tsv.gz"), "{" + name + "}"],
                       capture_output=True, text=True, check=True)
    hits = {}
    for l in p.stdout.splitlines():
        _, _, seqid, s, e, strand, cov, ident, copy = l.split("\t")
        hap = "#".join(seqid.split("#")[:2])
        hits.setdefault(hap, []).append(dict(
            contig=seqid, start=int(s) - 1, end=int(e), strand=strand,
            cov=float(cov) if cov != "NA" else 0.0, ident=float(ident) if ident != "NA" else 0.0, copy=int(copy)))
    return hits


# ---------------------------------------------------------------- shared: status, fetch, output

def flip(s):
    return "-" if s == "+" else "+"


def status(hits):
    """ok = one hit covering >= 95% of the CHM13 gene; otherwise comma-separated flags."""
    if not hits:
        return "missing"
    major = sum(h["cov"] >= COPY_COV for h in hits)
    flags = []
    if max(h["cov"] for h in hits) < OK_COV:
        flags.append("partial")
    if major > 1:
        flags.append("multi_copy")
    if len(hits) > max(major, 1):
        flags.append("fragments")
    return ",".join(flags) or "ok"


COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def faidx(fasta, fai, regions, workdir):
    """[(header, seq)] for regions 'contig:start1-end1', in order."""
    # a private .fai (rebuilt by setup offline) comes with a private .gzi next to it
    cmd = ["samtools", "faidx", fasta] + (["--fai-idx", fai, "--gzi-idx", fai[:-4] + ".gzi"] if fai else []) + regions
    for attempt in range(1, 4):
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
        if p.returncode == 0:
            break
        if attempt == 3:
            raise RuntimeError(f"samtools faidx {fasta}: {p.stderr.strip()[-300:]}")
        time.sleep(3 * attempt)
    out = []
    for rec in p.stdout.split(">")[1:]:
        head, _, seq = rec.partition("\n")
        out.append(seq.replace("\n", ""))
    return out


def get(args):
    data = args.data
    mode = args.mode
    if mode == "auto":
        mode = "offline" if os.path.exists(os.path.join(data, "offline", "liftoff_genes.tsv.gz.tbi")) else "online"
    need = os.path.join(data, mode, "assemblies.tsv")
    if not os.path.exists(need):
        sys.exit(f"run `genefinder.py setup {mode}` first ({need} missing)")
    t0 = time.time()
    loci = [l.split("\t") for l in open(os.path.join(data, "chm13_genes.tsv")) if l.split("\t", 1)[0] == args.gene]
    if not loci:
        sys.exit(f"{args.gene}: not in the CHM13 RefSeq annotation")
    c, s, e, st = loci[0][1:5]
    gene = dict(chrom=c, start=int(s), end=int(e), strand=st.strip())
    if len(loci) > 1:
        log(f"{args.gene}: {len(loci)} CHM13 loci, using {c}:{int(s) + 1}-{e}")
    log(f"{args.gene} CHM13 {c}:{int(s) + 1}-{e}({gene['strand']}), mode {mode}")

    asm = dict(tuple(l.rstrip("\n").split("\t", 1)) for l in open(need))
    if mode == "online":
        hits, annotated = find_online(gene, data, args.jobs), set(asm)
    else:
        hits = find_offline(args.gene, data)
        annotated = set(open(os.path.join(data, "offline", "annotated.txt")).read().split())
    annotated.add("CHM13#0")
    hits["CHM13#0"] = [dict(contig=f"CHM13#0#{c}", start=gene["start"], end=gene["end"], strand=gene["strand"],
                            cov=1.0, ident=1.0)]
    t_find = time.time() - t0

    # one samtools call per assembly, all its hits at once
    os.makedirs(args.outdir, exist_ok=True)
    work = os.path.join(args.outdir, f".{args.gene}.work")
    os.makedirs(work, exist_ok=True)
    flank = args.flank

    def fetch(hap):
        src = asm[hap].split("\t")
        fasta, fai = src[0], (src[1] if len(src) > 1 and src[1] != src[0] + ".fai" else None)
        regs = [f"{h['contig']}:{max(1, h['start'] + 1 - flank)}-{h['end'] + flank}" for h in hits[hap]]
        return hap, faidx(fasta, fai, regs, work)

    todo = [h for h in hits if h in asm and hits[h]]
    seqs = {}
    with ThreadPoolExecutor(args.jobs) as ex:
        for hap, ss in ex.map(fetch, todo):
            seqs[hap] = ss
    shutil.rmtree(work, ignore_errors=True)
    t_fetch = time.time() - t0 - t_find

    order = ["CHM13#0"] + sorted(h for h in set(asm) | set(hits) if h != "CHM13#0")
    fa = open(os.path.join(args.outdir, f"{args.gene}.fa"), "w")
    tsv = open(os.path.join(args.outdir, f"{args.gene}.tsv"), "w")
    tsv.write("gene\thap\tstatus\tcontig\tstart\tend\tstrand\tcov\tidentity\tlength\n")
    counts = {}
    for hap in order:
        hs = hits.get(hap)
        if hap not in asm:
            stat = "no_assembly"
        elif hap not in annotated:
            stat = "no_annotation"
        else:
            stat = "control" if hap == "CHM13#0" else status(hs or [])
        counts[stat] = counts.get(stat, 0) + 1
        if not hs or hap not in seqs:
            tsv.write(f"{args.gene}\t{hap}\t{stat}\t\t\t\t\t\t\t\n")
            continue
        for i, (h, seq) in enumerate(sorted(zip(hs, seqs[hap]), key=lambda z: -z[0]["cov"])):
            s1 = max(1, h["start"] + 1 - flank)
            e1 = s1 + len(seq) - 1
            if h["strand"] == "-":
                seq = seq.translate(COMP)[::-1]
            tag = "" if i == 0 else f" copy{i + 1}"
            fa.write(f">{args.gene}|{h['contig']}:{s1}-{e1}({h['strand']}){tag}\n")
            fa.write("\n".join(seq[k:k + 60] for k in range(0, len(seq), 60)) + "\n")
            tsv.write(f"{args.gene}\t{hap}\t{stat}\t{h['contig']}\t{h['start']}\t{h['end']}\t{h['strand']}"
                      f"\t{h['cov']:.3f}\t{h['ident']:.4f}\t{h['end'] - h['start']}\n")
    fa.close()
    tsv.close()
    log(f"{args.gene}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])))
    log(f"time: find {t_find:.1f}s, fetch {t_fetch:.1f}s ({len(todo)} assemblies, -j {args.jobs}), "
        f"total {time.time() - t0:.1f}s -> {args.outdir}/{args.gene}.fa")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("GENEFINDER_DATA", os.path.join(HERE, "data")))
    ap.add_argument("-j", "--jobs", type=int, default=32)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("setup", help="one-time pre-build")
    s.add_argument("which", choices=["online", "offline"])
    g = sub.add_parser("get", help="extract one gene")
    g.add_argument("gene", help="gene_name in the CHM13 RefSeq annotation, e.g. MALAT1")
    g.add_argument("--mode", choices=["auto", "online", "offline"], default="auto",
                   help="auto = offline if `setup offline` was run, else online")
    g.add_argument("--flank", type=int, default=0, help="bp added on each side (default 0)")
    g.add_argument("-o", "--outdir", default="genefinder_out")
    a = ap.parse_args()
    a.data = os.path.abspath(a.data)  # samtools runs from a scratch dir
    os.makedirs(a.data, exist_ok=True)
    if a.cmd == "setup":
        t0 = time.time()
        setup_online(a.data) if a.which == "online" else setup_offline(a.data, a.jobs)
        log(f"setup {a.which}: {time.time() - t0:.0f}s, data in {a.data}")
    else:
        get(a)


if __name__ == "__main__":
    main()
