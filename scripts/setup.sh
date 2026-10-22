#!/usr/bin/env bash
# One-time setup: check tools, environment report, small downloads, PAF offset index.
# Idempotent: existing files are kept. Each download/index step is timed into LOG_DIR/time.<step>.txt
# when GNU time is available.
# Usage: setup.sh REPO_DIR DATA_DIR LOG_DIR [FULL=0|1]
#   default : GENCODE GFF3 + assembly index + CHM13 .fai (~0.16 GB) and the PAF offset index (~0.12 GB),
#             built by streaming the 5.6 GB PAF once (nothing but the index is stored)
#   FULL=1  : also store the PAF (5.6 GB) and build the impg index (0.66 GB) for --mode impg/local
# URLs come from config/sources.env (override via environment variables).
set -euo pipefail
ROOT=$1; DATA=$2; LOG=$3; FULL=${4:-0}
source "$ROOT/config/sources.env"
source "$ROOT/scripts/compat.sh"
mkdir -p "$LOG" "$DATA"/{paf,annotation,assembly_index,reference}
T() {  # time a step if GNU time is available
  local s=$1; shift
  if [[ -n $GNU_TIME ]]; then "$GNU_TIME" -v -o "$LOG/time.$s.txt" "$@"; else "$@"; fi
}

# --- tools come from the user's environment (conda env create -f environment.yml && conda activate genefinder)
missing=()
for t in samtools bgzip tabix python3 curl; do command -v "$t" >/dev/null || missing+=("$t"); done
if ((FULL)) && ! command -v impg >/dev/null; then missing+=(impg); fi
if ((${#missing[@]})); then
  echo "ERROR: missing tools: ${missing[*]}"
  echo "  conda env create -f $ROOT/environment.yml && conda activate genefinder"
  echo "  (or: conda install -c conda-forge -c bioconda samtools htslib minimap2 impg)"
  exit 1
fi
samtools --version | grep -q 'libcurl=yes' || { echo "ERROR: samtools/htslib lacks libcurl (needed to read remote assemblies)"; exit 1; }

echo "== environment ($(date '+%Y-%m-%dT%H:%M:%S%z'))"
echo "host: $(hostname)  os: $(os_name)  kernel: $(uname -r)  arch: $(uname -m)"
echo "cpus: $(ncpu)"
echo "ram_total_gb: $(ram_gb)"
df -h "$DATA" | awk 'NR==2{print "disk ("$NF"): size "$2", used "$3", avail "$4}'
echo "samtools: $(samtools --version | head -1) / $(samtools --version | sed -n 2p) / libcurl=yes"
echo "bgzip: $(bgzip --version | head -1); tabix: $(tabix --version | head -1)"
echo "impg: $(impg --version 2>/dev/null || echo 'not installed (optional: only for --mode impg)')"
echo "minimap2: $(minimap2 --version 2>/dev/null || echo 'not installed (optional: only for scripts/validate.sh)')"
echo "python3: $(python3 --version 2>&1)"
echo "GNU time: ${GNU_TIME:-not found (per-step resource logs disabled; on macOS: brew install gnu-time)}"
echo "curl: $(curl --version | head -1)"

# --- small downloads (URLs: config/sources.env)
[[ -s $DATA/annotation/$GENCODE_GFF3 ]] || T annotation_gencode \
  curl -sS --fail --retry 5 -o "$DATA/annotation/$GENCODE_GFF3" "$ANNOTATION_BASE_URL/$GENCODE_GFF3"
IDX=$DATA/assembly_index/assemblies_release2_v1.0.index.csv
[[ -s $IDX ]] || T assembly_index curl -sS --fail --retry 5 -o "$IDX" "$ASSEMBLY_INDEX_URL"
[[ -s $DATA/reference/chm13.fa.gz.fai ]] || curl -sS --fail -o "$DATA/reference/chm13.fa.gz.fai" "$CHM13_FASTA_URL.fai"

# --- PAF offset index: where every alignment record sits in the remote PAF, so a gene only needs the
#     few MB of records that overlap it (scripts/remote_project.py). Built by reading the PAF once.
PAF=$DATA/paf/hprc465vschm13.aln.paf.gz
OFF=$DATA/paf/hprc465vschm13.offsets.bed.gz
if ((FULL)); then
  # The published .paf.gz is already BGZF (gzip header carries the BC extra field), so impg reads it directly.
  [[ -s $PAF ]] || T paf_download curl -sS --fail --retry 5 -C - -o "$PAF" "$PAF_URL"
  [[ -s $PAF.impg ]] || T impg_index impg index -a "$PAF" -t "$(ncpu)"
fi
if [[ ! -s $OFF.tbi ]]; then
  if [[ -s $PAF ]]; then
    echo "== building PAF offset index from the local PAF"
    reader="cat '$PAF'"
  else
    echo "== building PAF offset index: streaming the 5.6 GB PAF once over HTTPS (only the ~0.12 GB index is kept)"
    reader="curl -sS --fail --retry 5 '$PAF_URL'"
  fi
  T paf_offset_index bash -c "set -o pipefail; $reader | python3 '$ROOT/scripts/build_paf_index.py' - \
    | LC_ALL=C sort -k1,1 -k2,2n -S 1G | bgzip > '$OFF.tmp'"
  mv "$OFF.tmp" "$OFF"
  tabix -p bed "$OFF"
fi
echo "== setup complete; data in $DATA ($(du -sh "$DATA" | cut -f1))"
ls -l "$DATA"/*/
