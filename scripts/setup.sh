#!/usr/bin/env bash
# One-time setup: tools, environment report, downloads (annotations, assembly index, PAF), impg index.
# Idempotent: existing files are kept. Each download/index step is timed into LOG_DIR/time.<step>.txt.
# Usage: setup.sh REPO_DIR DATA_DIR LOG_DIR
# URLs come from config/sources.env (override via environment variables).
set -euo pipefail
ROOT=$1; DATA=$2; LOG=$3
source "$ROOT/config/sources.env"
mkdir -p "$LOG" "$DATA"/{paf,annotation,assembly_index,reference}
T() {  # time a step if GNU time is available
  local s=$1; shift
  if [[ -x /usr/bin/time ]]; then /usr/bin/time -v -o "$LOG/time.$s.txt" "$@"; else "$@"; fi
}

# --- tools: use what is on PATH, else install a private micromamba env into REPO_DIR/tools/env
[[ -x $ROOT/tools/env/bin/impg ]] && export PATH=$ROOT/tools/env/bin:$PATH
missing=()
for t in samtools bgzip impg minimap2 python3 curl; do command -v "$t" >/dev/null || missing+=("$t"); done
if ((${#missing[@]})); then
  echo "== installing tools into $ROOT/tools/env (missing: ${missing[*]})"
  mkdir -p "$ROOT/tools"
  [[ -x $ROOT/tools/bin/micromamba ]] || (cd "$ROOT/tools" && curl -sL "$MICROMAMBA_URL" | tar -xj bin/micromamba)
  MAMBA_ROOT_PREFIX=$ROOT/tools/mamba T conda_env "$ROOT/tools/bin/micromamba" create -y -p "$ROOT/tools/env" \
    -f "$ROOT/environment.yml"
  export PATH=$ROOT/tools/env/bin:$PATH
fi
samtools --version | grep -q 'libcurl=yes' || { echo "ERROR: samtools/htslib lacks libcurl (needed for remote faidx)"; exit 1; }

echo "== environment ($(date -Is))"
echo "host: $(hostname)  os: $(. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-$(uname -s)}")  kernel: $(uname -r)"
echo "cpus: $(nproc)"
free -g 2>/dev/null | awk 'NR==2{print "ram_total_gb: "$2"  ram_available_gb: "$7}'
df -h "$DATA" | awk 'NR==2{print "disk ("$6"): size "$2", used "$3", avail "$4}'
echo "samtools: $(samtools --version | head -1) / $(samtools --version | sed -n 2p) / libcurl=yes"
echo "bgzip: $(bgzip --version | head -1)"
echo "impg: $(impg --version)"
echo "minimap2: $(minimap2 --version)"
echo "python3: $(python3 --version 2>&1)"
echo "GNU time: $( [[ -x /usr/bin/time ]] && /usr/bin/time --version 2>&1 | head -1 || echo 'not found (resource logs disabled)')"
echo "curl: $(curl --version | head -1)"

# --- downloads (URLs: config/sources.env)
for f in "$GENCODE_GFF3" "$GENCODE_GFF3.tbi" "$REFSEQ_GFF3"; do
  [[ -s $DATA/annotation/$f ]] || T "annotation_${f%%.*}" curl -sS --fail --retry 5 -o "$DATA/annotation/$f" "$ANNOTATION_BASE_URL/$f"
done
IDX=$DATA/assembly_index/assemblies_release2_v1.0.index.csv
[[ -s $IDX ]] || T assembly_index curl -sS --fail --retry 5 -o "$IDX" "$ASSEMBLY_INDEX_URL"
[[ -s $DATA/reference/chm13.fa.gz.fai ]] || curl -sS --fail -o "$DATA/reference/chm13.fa.gz.fai" "$CHM13_FASTA_URL.fai"
PAF=$DATA/paf/hprc465vschm13.aln.paf.gz
if [[ ${SKIP_PAF:-0} == 1 ]]; then
  echo "== SKIP_PAF=1: not downloading the PAF (use run.sh --mode stream)"
else
  # The published .paf.gz is already BGZF (gzip header carries the BC extra field), so impg reads it directly.
  [[ -s $PAF ]] || T paf_download curl -sS --fail --retry 5 -C - -o "$PAF" "$PAF_URL"
  [[ -s $PAF.impg ]] || T impg_index impg index -a "$PAF" -t "$(nproc)"
fi
echo "== setup complete; data in $DATA"
ls -l "$DATA"/*/
