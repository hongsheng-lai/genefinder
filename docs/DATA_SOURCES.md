# Data sources

genefinder only reads **public** data over anonymous HTTPS. No AWS account, credentials, or `aws` CLI are needed.

- **Where the URLs live:** every URL is defined once, in [`config/sources.env`](../config/sources.env). Any of them can be overridden with an environment variable of the same name.
- **Checking them:** run `scripts/check_sources.sh` to confirm they are reachable and unchanged. It compares size, ETag and Last-Modified, and spot-checks assemblies for `.fai` and `.gzi`.
- **Browsing the bucket:** `scripts/s3ls.py <prefix> [--delimiter /]` lists `s3://human-pangenomics` over HTTPS. The S3 web view is at <https://human-pangenomics.s3.amazonaws.com/index.html>.

Sizes and ETags below were recorded on 2026-10-21.

## Downloaded once by `./run.sh --setup`

### Release 2 all-haplotypes vs CHM13 alignments (`PAF_URL`)

- **URL:** https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/impg/pafs/hprc465vschm13.aln.paf.gz
- **Size:** 5,590,221,127 bytes, ETag `fbacddb3eb2fcb399a9202e8eab2108a-667`
- **Contents:** 7.15 M records (17.4 GB uncompressed) covering 465 haplotypes: 464 HPRC plus GRCh38.
- **Naming:** query = haplotype contig in PanSN naming (`HG00438#1#CM...`); target = `CHM13#0#chrN`.
- **Tags:** `cg:Z` CIGAR (`=`/`X`/`I`/`D`), `gi`/`bi`/`md` identity, `ch:Z` chain id.
- **Format:** already BGZF (the gzip header has the `BC` field), so `impg index` reads it directly.
- **Used for:** projecting CHM13 gene coordinates onto each haplotype.
- **Without a local copy:** `--mode stream` streams it per run instead.

### Release 2 assembly index (`ASSEMBLY_INDEX_URL`)

- **URL:** https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/2a7e3a5b1ca10454762699176b1a0a1027ccbf21/data_tables/assemblies_release2_v1.0.index.csv
- **Version:** pinned to commit `2a7e3a5` of [human-pangenomics/hprc_intermediate_assembly](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_release2_v1.0.index.csv); md5 `90c9e66a58580484facd59d3ed28c52c`.
- **Contents:** 466 rows (232 samples × 2 haplotypes, plus CHM13 and GRCh38) with S3 paths of each `*.fa.gz`, `.fai` and `.gzi`.
- **Used for:** mapping each PanSN `sample#hap` to its FASTA URL.

### CHM13 v2.0 annotation, GENCODE (`GENCODE_GFF3`)

- **URL:** https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz (plus `.tbi`)
- **Size:** 158,318,946 bytes
- **Contents:** GENCODE v35 mapped to CHM13 with CAT/Liftoff. This is the default (`--annotation gencode`).

### CHM13 v2.0 annotation, RefSeq (`REFSEQ_GFF3`)

- **URL:** https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz
- **Size:** 56,368,414 bytes
- **Contents:** RefSeq genes lifted to CHM13 (`--annotation refseq`).

### CHM13 analysis set `.fai` (`CHM13_FASTA_URL`)

- **URL:** https://human-pangenomics.s3.amazonaws.com/working/HPRC_PLUS/CHM13/assemblies/analysis_set/pansn/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz.fai
- **Size:** 988 bytes
- **Used for:** chromosome sizes, so padded intervals are clipped at chromosome ends.

**Disk after setup:** about 6.5 GB, made up of the PAF (5.6 GB), the impg index (0.66 GB) and the annotations (0.2 GB).

## Read remotely on every run (never downloaded whole)

The 466 assemblies are listed in the index above, under `https://human-pangenomics.s3.amazonaws.com/working/HPRC/<SAMPLE>/assemblies/release2/<SAMPLE>_<mat|pat|hap1|hap2>_hprc_r2_v1.0.1.fa.gz`, plus the two references under `working/HPRC_PLUS/`.

- **Format:** each is bgzipped, about 0.9 GB, with `.fai` (~4 KB) and `.gzi` (~0.75 MB). Contig names are already PanSN.
- **`--fetch faidx`:** `samtools faidx <URL> -r regions.txt` issues HTTP range requests. That costs ~0.2 MB and ~0.4 s per region, plus the `.gzi` once per assembly per batch.
- **`--fetch stream`:** `curl <URL> | bgzip -dc` reads each assembly once, extracting all regions on the fly. That costs ~0.9 GB per assembly per batch, and the assembly is never stored.

## Related files not used by default

All are in the same Release 2 folder (`https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/`):

| File | Notes |
|---|---|
| `impg/pafs/hprc465vsgrch38.aln.paf.gz` | the same haplotypes aligned to GRCh38 (6.3 GB) |
| `impg/tpas/*.tpa` | compact per-haplotype alignments |
| `minigraph-cactus/hprc-v2.0-mc-chm13.gbz` | the Minigraph-Cactus graph (5.7 GB) |

## Citing the data

When you use outputs of this pipeline, please cite the data providers:
- the Human Pangenome Reference Consortium (Release 2 assemblies and alignments)
- the T2T Consortium (CHM13 v2.0 and its annotations)
- the tools: impg, samtools/htslib, minimap2

Follow the [HPRC data use policy](https://human-pangenomics.s3.amazonaws.com/submissions/HPRC_DATA_USE_POLICY_v1.0).
