# genefinder

Extract the sequence of any CHM13 v2.0-annotated gene from **every HPRC Release 2 haplotype assembly** (464 haplotypes, plus CHM13 and GRCh38 as controls), **without downloading whole assemblies**.

```
gene name(s) ──► CHM13 BED ──► project through HPRC R2 all-vs-CHM13 alignments (impg) ──► interval per haplotype
                                                                                                │
   per-gene FASTA + coords + QC flags ◄── remote samtools faidx / single-pass stream of each assembly
```

It works in four steps:
1. **Look up the gene.** Genes are found in the CHM13 v2.0 GENCODE (or RefSeq) annotation using the annotation's own coordinates, not GRCh38.
2. **Project.** The gene interval is projected onto every haplotype through the Release 2 alignments (`hprc465vschm13.aln.paf.gz`).
3. **Fetch.** The matching interval is read from each haplotype's bgzipped, indexed assembly on the public HPRC bucket, in the gene's sense orientation.
4. **Flag.** Each haplotype is reported as missing, partial, split across contigs or alignment blocks, or multi-copy. Nothing is dropped silently.

All inputs are public HTTPS. **[docs/DATA_SOURCES.md](docs/DATA_SOURCES.md)** lists every file, URL, size and checksum, and they are all set in [`config/sources.env`](config/sources.env).

## Requirements

- Linux with `bash`, `curl` and Python ≥ 3.9 (standard library only).
- [samtools/htslib](https://www.htslib.org/) ≥ 1.16 **built with libcurl**, [impg](https://github.com/pangenome/impg) ≥ 0.5.0, and [minimap2](https://github.com/lh3/minimap2) for validation. `./run.sh --setup` installs these automatically into `tools/env` with micromamba if they are not on your `PATH`. Alternatively, `mamba env create -f environment.yml`.
- About 7 GB of disk for the one-time data (alignment file plus index), and outbound HTTPS to `human-pangenomics.s3.amazonaws.com`.
- Optional: GNU time (`/usr/bin/time`) for per-step resource logs.

## Quick start

```bash
git clone https://github.com/hongsheng-lai/genefinder.git
cd genefinder
scripts/check_sources.sh            # optional: verify the remote data is reachable / unchanged
./run.sh --setup                    # one time, ~20 min: tools, annotations, assembly index, PAF (5.6 GB), impg index
./run.sh MALAT1                     # results in ./genefinder_out/genes/MALAT1/
scripts/validate.sh MALAT1          # minimap2 a sample of the sequences against CHM13 (expect + strand, >99% id)
```

More examples:

```bash
./run.sh MALAT1 NEAT1 GAPDH --outdir results -j 32       # several genes, 32 parallel jobs
./run.sh --gene-list genes.txt --annotation refseq        # a list (first column = gene name or ID)
./run.sh --biotype protein_coding ALL --fetch stream -j 64   # the whole annotation (see "Scaling")
./run.sh --mode stream MALAT1                             # no local PAF at all: stream it from S3 (slow, ~12 min)
```

- **Lookup:** genes are matched by `gene_name`, `Name`, `gene`, `gene_id`, `ID`, or Ensembl ID (with or without version).
- **Reporting:** genes that aren't found, or that match several annotated loci, are reported in `work/<batch>/batch.lookup.txt`. Extra loci become `<GENE>_copy2` and so on.
- **Rerunning:** reruns resume. Assemblies already fetched in a batch are skipped, and failures are retried.

### Options (`./run.sh --help`)

| Option | Default | Meaning |
|---|---|---|
| `--outdir DIR` | `./genefinder_out` | results |
| `--data-dir DIR` | `<repo>/data` (`$GENEFINDER_DATA`) | one-time downloads and impg index |
| `--annotation gencode\|refseq\|PATH` | `gencode` | CHM13 v2.0 GFF3 (GENCODE v35 CAT/Liftoff, or RefSeq Liftoff) |
| `--flank N` | 5000 | flank on each side for `*.padded.fa` |
| `--biotype B` | – | keep only genes with this biotype (e.g. `protein_coding`) |
| `--mode impg\|local\|stream` | `impg` | projection engine (see below) |
| `--fetch auto\|faidx\|stream` | `auto` | sequence access (see below) |
| `-j N` | 16 | parallel impg processes and assembly jobs; fetching is latency-bound, so 32–64 is faster |
| `--batch NAME` | first gene / list name | names `work/` and `logs/` subdirectories |

## Outputs

Each gene gets a folder, `OUTDIR/genes/<GENE>/`:

| File | Content |
|---|---|
| `<GENE>.gene.fa` | gene body from each haplotype, **in the gene's sense orientation** (reverse-complemented where needed). Header: `>GENE\|sample#hap#contig:start-end(strand)`, 1-based inclusive, with strand = the gene's orientation on that contig. CHM13 comes first, then GRCh38, then HPRC. Extra copies or fragments get a suffix (` copy2`, ` fragment3`). |
| `<GENE>.padded.fa` | the same hits ± flank |
| `<GENE>.coords.tsv` | one row per hit and per missing haplotype; columns listed below |
| `<GENE>.summary.txt` | counts (attempted, extracted, ok, missing, partial, split, multi-copy), length distribution, controls, and every flagged haplotype with its hits |

Columns of `coords.tsv`:
- `gene`, `sample`, `hap`, `contig`
- `start`, `end`: 0-based, half-open
- `strand`
- `gene_cov_frac`: CHM13 gene bases covered
- `n_blocks`
- `identity`: matches / (matches + mismatches + gap opens)
- `padded_start`, `padded_end`, `length`
- `hit_class`: primary, copy or fragment
- `hap_status`, `extracted`, `fasta_header`

`OUTDIR/work/<batch>/` holds the intermediates: BEDs, `projection.paf`, `regions.tsv`, the fetch log and `batch_summary.tsv`. `OUTDIR/logs/<batch>/` holds the GNU `time -v` log of every step and `resources.tsv`.

### Haplotype status flags

Flags are comma-separated.

| Flag | Meaning |
|---|---|
| `ok` | one hit, one contig, one alignment block, ≥ 95% of the gene covered |
| `missing` | no alignment covers the gene. Expected for chrX genes in male paternal haplotypes, and for chrY-PAR genes (masked in the CHM13 target). |
| `partial` | best hit covers < 95% of the gene; all fragments are then extracted |
| `split_contigs` | the gene is only covered by combining hits on more than one contig |
| `split_blocks` | the primary hit is made of more than one alignment piece, usually an indel or SV inside the gene. The sequence is still complete. |
| `multi_copy` | more than one hit covers ≥ 50% of the gene; every copy is extracted |
| `control` | CHM13 itself |

## How it works

| Step | Script | Notes |
|---|---|---|
| setup | `scripts/setup.sh` | tools, environment report, downloads, `impg index` |
| genes → BED6 | `scripts/gene2bed.py` | `CHM13#0#chrN`, 0-based; gene and padded BEDs; not-found / multi-hit report |
| projection | `scripts/impg_project.sh` | `impg query --no-merge -o paf`, BED split across `-j` processes |
| | `scripts/paf_project.py`, `stream_project.sh` | pure-Python fallback: one CIGAR walk per alignment record for all overlapping genes (`--mode local` / `stream`) |
| hits and flags | `scripts/summarize_projection.py` | chains pieces per haplotype, contig and strand. A new hit starts if pieces are > 20 kb apart or re-cover CHM13 bases (a tandem copy). |
| fetch | `scripts/fetch_seqs.py`, `stream_extract.py` | one job per assembly per batch; the gene is sliced from its padded interval; retries; resume |
| summaries | `scripts/gene_summary.py` | `summary.txt` and `batch_summary.tsv` |
| validation | `scripts/validate.sh GENE [N] [OUTDIR]` | minimap2 of sampled sequences (plus GRCh38 and reverse-strand hits) against the CHM13 gene |
| utilities | `scripts/check_sources.sh`, `scripts/s3ls.py`, `scripts/resources.py` | check URLs, list the bucket, collect time logs |

### Choosing a projection mode

| `--mode` | needs | MALAT1 | best for |
|---|---|---|---|
| `impg` (default) | local PAF + impg index (6.2 GB, built once in ~17 min) | 0.2 s | anything beyond a couple of genes |
| `local` | local PAF only | ~20 s | when impg isn't available |
| `stream` | nothing stored | ~12 min, 5.6 GB downloaded each run | a single gene on a machine without disk to spare |

All three produce the same intervals. impg can differ by 1 bp at a gene end when an insertion sits exactly on the boundary.

### Choosing a fetch method

- **`faidx`:** HTTP range requests to each remote assembly. It's latency-bound at about 0.4 s and 0.2 MB per region, so cost grows with the number of genes. Best for up to a few hundred genes.
- **`stream`:** reads each assembly once (~0.9 GB, ~20–40 s) and extracts every region of the batch in that pass. The cost is fixed at ~420 GB of transfer per batch, however many genes it holds. Best for thousands of genes.
- **`auto` (default):** streams an assembly when the batch needs more than 3,000 regions from it.

## Performance

Reference run: 128 CPUs, NFS storage, `-j 16` unless noted.

| | wall | peak RAM | network | disk |
|---|---|---|---|---|
| One-time setup | ~18 min (14.4 min of it `impg index`) | 5.5 GB | 5.8 GB | 7.6 GB |
| One gene (MALAT1, 466 assemblies) | ~95 s (~35 s with `-j 64`) | < 50 MB | ~0.6 GB | 13 MB |
| 200 protein-coding genes (`-j 64`, faidx) | ~15 min | ~0.75 GB¹ | ~20 GB | ~16 GB (~4 GB bgzipped) |

¹ Measured for the bounded-memory fetch step on a 64-assembly × 200-gene batch.

For MALAT1, 464/464 haplotypes were extracted. The median length was 8,829 bp (equal to CHM13), and all 466 sequences align to CHM13 on the + strand at ≥ 99% identity. The one outlier, HG00558#1 at 9,190 bp, carries a polymorphic Alu insertion and is flagged `split_blocks`.

### Scaling to all genes

| Gene set | Sequence per haplotype (gene + padded) | Uncompressed FASTA, all haplotypes |
|---|---|---|
| All ~64k GENCODE genes | ~4.6 Gb | ~2.2 TB |
| Protein-coding (~20k) | ~2.9 Gb | ~1.4 TB |

- **Projection** with `-j 64` takes minutes.
- **Use `--fetch stream`** with large batches, so that each ~420 GB pass over the assemblies serves many genes.
- **Expect output volume to be the bottleneck.** The final write and summary steps are currently single-threaded (~60 MB/s), and the hits/flags step holds ~1.4 MB per gene in memory, so keep batches to ≤ 5–10k genes.

## Caveats

- **Boundaries.** Gene boundaries come from the alignments, so assembly errors and alignment breaks show up as `partial`, `split_*` or `multi_copy`. Check `summary.txt` before using sequences from difficult loci (LPA, AMY1, SMN, NBPF, …).
- **Identity metrics.** `identity` in `coords.tsv` is gap-compressed. minimap2 identity in `validate.sh` counts every indel base, so carriers of structural variants score lower.
- **chrY.** The CHM13 alignment target has the chrY PARs hard-masked, so PAR genes are projected only via their chrX copy.

## Tests

```bash
python3 -m unittest discover -s tests -v   # offline unit tests: CIGAR projection, stream extraction, hit chaining
```

## Acknowledgements and citation

This tool only re-uses public data and tools. Please cite:
- the **Human Pangenome Reference Consortium** (Release 2 assemblies and alignments)
- the **T2T Consortium** (CHM13 v2.0 and its CAT/Liftoff annotations)
- **impg**, **samtools/htslib** and **minimap2**

Use of HPRC data is subject to the [HPRC data use policy](https://human-pangenomics.s3.amazonaws.com/submissions/HPRC_DATA_USE_POLICY_v1.0).
