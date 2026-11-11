# genefinder

Extract one CHM13 gene from every HPRC Release 2 haplotype (464 haplotypes, plus CHM13 as control).
One script, `genefinder.py`, using Python's standard library plus samtools and tabix.

```bash
conda env create -f environment.yml && conda activate genefinder

# one-time pre-build: pick the one you need
./genefinder.py setup online      # anywhere with HTTPS (~3-5 min, keeps 120 MB)
./genefinder.py setup offline     # Salz cluster, reads /ccb/salz8/HPRC/assemblies (~2-3 min, keeps 400 MB)

./genefinder.py get MALAT1        # -> genefinder_out/MALAT1.fa + MALAT1.tsv
```

`get` options: `--mode auto|online|offline`, default auto: offline if `setup offline` was run, else online. Also `--flank N` (bp on each side, default 0), `-o DIR`, and `-j N` before `get` (parallel jobs, default 32). Data lives in `./data`; change it with `--data DIR` or `GENEFINDER_DATA`.

## Online vs offline

The two modes locate the gene in each haplotype differently.

| | online | offline (Salz) |
|---|---|---|
| gene → haplotype coordinates | the CHM13 gene interval is **projected through the Release 2 whole-genome alignments** (`hprc465vschm13.aln.paf.gz`, `cg:Z` CIGARs) | **looked up in each haplotype's own Liftoff annotation** (`annotation.liftoff/*.gff3.gz` in the mirror) |
| sequences | `samtools faidx` on the bgzipped assemblies in S3 (HTTP range requests) | `samtools faidx` on `/ccb/salz8/HPRC/assemblies/fasta` |
| pre-build | stream the 5.6 GB PAF once and keep a 120 MB offset index of its records | scan the 462 Liftoff GFF3s (24 GB) once and keep a name-indexed table of their gene lines |
| per gene | about 466 PAF records (~13 MB) fetched by byte range | one tabix lookup |
| coverage | 464 haplotypes, plus GRCh38 (it is in the PAF) | 462 haplotypes: HG002 has no Liftoff annotation, so it is reported as `no_annotation`, as is GRCh38 |
| `cov` / `identity` | CHM13 gene bases covered by the alignment / gap-compressed identity | Liftoff `coverage` / `sequence_ID` |

Both modes look genes up by `gene_name` in the same RefSeq namespace. Online uses the CHM13 RefSeq Liftoff v5.2 annotation, and the per-haplotype files were lifted from CHM13 RefSeq. The CHM13 control comes from that CHM13 annotation in both modes.

The boundaries could differ in principle: projection reproduces the CHM13 gene ends, while Liftoff places the lifted gene model's ends. On the genes benchmarked below they are identical. `tests/compare.py` reports any difference.

## Output

`OUTDIR/GENE.fa`: one record per hit, **in the gene's sense orientation** (reverse-complemented where needed), CHM13 first.
The header is `>GENE|sample#hap#contig:start-end(strand)`: 1-based inclusive coordinates, and the gene's strand on that contig. Extra copies get ` copy2`, ` copy3`, …

`OUTDIR/GENE.tsv`: one row per hit and per haplotype without one.
- `gene`, `hap`, `status`, `contig`
- `start`, `end`: 0-based half-open
- `strand`, `cov`, `identity`, `length`

| status | meaning |
|---|---|
| `ok` | one hit covering ≥ 95% of the gene |
| `partial` | best hit covers < 95% |
| `multi_copy` | ≥ 2 hits covering ≥ 50% each (all extracted) |
| `fragments` | extra hits covering < 50% (all extracted) |
| `missing` | gene not found (e.g. chrX genes in male paternal haplotypes) |
| `no_annotation` | offline only: no Liftoff file for this haplotype (HG002, GRCh38) |
| `control` | CHM13 |

## Benchmark

```bash
tests/benchmark.sh MALAT1 32      # fresh pre-build of both modes, then the gene online and offline, twice each
python3 -m unittest discover -s tests
```

Reference run, Salz cluster, MALAT1 (8.8 kb, chr11), `-j 32`, starting from an empty data dir:

| step | online | offline |
|---|---|---|
| pre-build (one time) | 165 s (5.6 GB PAF streamed from S3; 160–280 s over three runs, depending on S3 throughput) | 133–135 s (scan of 462 Liftoff GFF3s, 24 GB) |
| data kept | 119 MB | 402 MB |
| query MALAT1 (8.8 kb) | 47.7 s (repeat: 45.3 s) | 0.6 s (repeat: 0.7 s) |

Online time is mostly network latency: 6.7 s to locate the gene (466 PAF records, 12.8 MB in range requests) and ~40 s to fetch 466 sequences from S3 with samtools.

How the two methods compare (`tests/compare.py`; haplotypes found by both modes):

| gene | strand | haplotypes in both | identical coordinates and sequence | only online | extra online hits |
|---|---|---|---|---|---|
| MALAT1 | + | 462 | 462 | GRCh38, HG002 ×2 | none |
| TP53 | − | 462 | 462 | GRCh38, HG002 ×2 | HG01074#1: a second, partial (20%) copy 100 kb away, `fragments` |

The two modes agree base for base wherever both find the gene. Online also covers GRCh38 and HG002, which have no Liftoff annotation. It can also flag extra partial copies that Liftoff does not place.

## Data sources

All are public; nothing needs credentials.
- Alignments: `https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/impg/pafs/hprc465vschm13.aln.paf.gz` (5.6 GB, BGZF; override with `PAF_URL`).
- Assembly index: `hprc_intermediate_assembly` commit 2a7e3a5, `data_tables/assemblies_release2_v1.0.index.csv` (override with `ASSEMBLY_INDEX_URL`).
- CHM13 genes: `T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` (56 MB; override with `CHM13_GFF`, a URL or a local path).
- Offline mirror: `/ccb/salz8/HPRC/assemblies` (override with `HPRC_MIRROR`), which provides `fasta/` and `annotation.liftoff/`. The mirror's `NA18940_hap2…fa.gz.fai` is empty, so `setup offline` rebuilds a private copy in `data/offline/`.
