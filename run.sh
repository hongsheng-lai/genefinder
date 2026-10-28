#!/usr/bin/env bash
# genefinder: extract CHM13-annotated gene sequences from every HPRC Release 2 haplotype assembly.
#
#   ./run.sh --setup [--data-dir DIR]            # one time (~0.3 GB): annotation, assembly index, PAF offset index
#   ./run.sh --setup --full                      # optional: also store the PAF + impg index (6.3 GB) for huge batches
#   ./run.sh [options] GENE [GENE ...]
#   ./run.sh [options] --gene-list genes.txt
#   ./run.sh [options] ALL                       # every gene in the annotation (use --biotype to restrict)
#
# Options:
#   --outdir DIR                       results directory (default ./genefinder_out)
#   --data-dir DIR                     one-time data from --setup (default <repo>/data, or $GENEFINDER_DATA)
#   --annotation gencode|refseq|PATH   CHM13 v2.0 GFF3 (default gencode = GENCODEv35 CAT/Liftoff)
#   --flank N                          flank for *.padded.fa (default 5000)
#   --biotype B                        only genes with this biotype (e.g. protein_coding)
#   --mode auto|remote|impg|local|stream  projection engine. auto (default): impg if the PAF is stored
#                                      (setup --full) and impg is installed, else local if the PAF is stored,
#                                      else remote. remote: fetch only the alignment records overlapping the
#                                      genes from the S3 PAF (HTTP ranges, via the offset index); impg / local:
#                                      use the stored PAF; stream: read the whole PAF from S3 each run
#   --fasta-dir DIR|none               local copies of the assemblies (bgzip + .fai + .gzi, index file names);
#                                      default $LOCAL_FASTA_DIR (config/sources.env) when it exists, else none (S3)
#   --offline                          no network: require a stored PAF (impg/local) and every assembly in
#                                      --fasta-dir (e.g. on the Salz cluster after setup --full)
#   --fetch auto|faidx|stream          sequence fetch: faidx = remote samtools faidx (range requests, best
#                                      for few genes); stream = one curl|bgzip pass per assembly (best for
#                                      thousands of genes); auto = stream if >3000 regions/assembly
#   --batch NAME                       work/log subdirectory name (default: first gene or list basename)
#   -j N                               parallel jobs for projection and assembly fetch (default 16)
#
# Outputs: OUTDIR/genes/<GENE>/{<GENE>.gene.fa,<GENE>.padded.fa,<GENE>.coords.tsv,<GENE>.summary.txt}
#          OUTDIR/work/<batch>/ intermediates; OUTDIR/logs/<batch>/ GNU time -v per step + resources.tsv
# Data sources (all public HTTPS): config/sources.env, docs/DATA_SOURCES.md
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
S=$ROOT/scripts
source "$ROOT/config/sources.env"
source "$S/compat.sh"

ANN=gencode; FLANK=5000; BIOTYPE=; MODE=auto; FASTA_DIR=$LOCAL_FASTA_DIR; OFFLINE=0; BATCH=; JOBS=16; LIST=; SETUP=0; FULL=0; GENES=(); FETCH=auto
OUT=$PWD/genefinder_out; DATA=${GENEFINDER_DATA:-$ROOT/data}
while (($#)); do
  case $1 in
    --outdir) OUT=$2; shift 2;;
    --data-dir) DATA=$2; shift 2;;
    --annotation) ANN=$2; shift 2;;
    --flank) FLANK=$2; shift 2;;
    --biotype) BIOTYPE=$2; shift 2;;
    --mode) MODE=$2; shift 2;;
    --batch) BATCH=$2; shift 2;;
    -j) JOBS=$2; shift 2;;
    --fetch) FETCH=$2; shift 2;;
    --gene-list) LIST=$2; shift 2;;
    --fasta-dir) FASTA_DIR=$2; shift 2;;
    --offline) OFFLINE=1; shift;;
    --setup) SETUP=1; shift;;
    --full) FULL=1; shift;;
    -h|--help) sed -n '2,33p' "$0"; exit 0;;
    -*) echo "unknown option $1 (see --help)"; exit 1;;
    *) GENES+=("$1"); shift;;
  esac
done
mkdir -p "$DATA"; DATA=$(cd "$DATA" && pwd)
PAF=$DATA/paf/hprc465vschm13.aln.paf.gz
OFF=$DATA/paf/hprc465vschm13.offsets.bed.gz
IDX=$DATA/assembly_index/assemblies_release2_v1.0.index.csv

if ((SETUP)); then
  mkdir -p "$DATA/logs"
  bash "$S/setup.sh" "$ROOT" "$DATA" "$DATA/logs" "$FULL" | tee "$DATA/logs/environment.txt"
  exit 0
fi

[[ ${#GENES[@]} -gt 0 || -n $LIST ]] || { sed -n '2,33p' "$0"; exit 1; }
[[ -s $IDX ]] || { echo "[run] data not set up in $DATA: run ./run.sh --setup first"; exit 1; }
for t in samtools bgzip tabix python3 curl; do command -v $t >/dev/null || { echo "[run] $t not on PATH (run --setup)"; exit 1; }; done
if [[ $FASTA_DIR == none || ! -d $FASTA_DIR ]]; then FASTA_DIR=; fi
if [[ $MODE == auto ]]; then
  if [[ -f $PAF.impg ]] && command -v impg >/dev/null; then MODE=impg
  elif [[ -s $PAF ]]; then MODE=local
  else MODE=remote; fi
fi
if ((OFFLINE)); then
  [[ $MODE == impg || $MODE == local ]] || { echo "[run] --offline needs a stored PAF (./run.sh --setup --full) and --mode impg|local"; exit 1; }
  [[ -n $FASTA_DIR ]] || { echo "[run] --offline needs --fasta-dir (local assemblies); $LOCAL_FASTA_DIR not found"; exit 1; }
  [[ $ANN != refseq || -s $DATA/annotation/$REFSEQ_GFF3 ]] || { echo "[run] --offline: RefSeq annotation not downloaded"; exit 1; }
fi
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
BATCH=${BATCH:-$( [[ -n $LIST ]] && basename "${LIST%.*}" || echo "${GENES[0]}" )}
W=$OUT/work/$BATCH; LOG=$OUT/logs/$BATCH; mkdir -p "$W" "$LOG" "$OUT/genes"
echo "[run] batch=$BATCH mode=$MODE fetch=$FETCH annotation=$ANN flank=$FLANK jobs=$JOBS outdir=$OUT"
echo "[run] assemblies: ${FASTA_DIR:-S3 ($HPRC_BUCKET)}"

timed() {  # timed STEP cmd...  -> LOG/time.STEP.txt (GNU time -v, when available)
  local step=$1; shift
  if [[ -n $GNU_TIME ]]; then
    "$GNU_TIME" -v -o "$LOG/time.$step.txt" "$@"
    echo "Work dir bytes: $(dir_bytes "$W")" >> "$LOG/time.$step.txt"
  else
    "$@"
  fi
}

# RefSeq annotation is fetched on first use (setup only downloads GENCODE)
if [[ $ANN == refseq && ! -s $DATA/annotation/$REFSEQ_GFF3 ]]; then
  echo "[run] downloading $REFSEQ_GFF3 (56 MB)"
  curl -sS --fail --retry 5 -o "$DATA/annotation/$REFSEQ_GFF3.tmp" "$ANNOTATION_BASE_URL/$REFSEQ_GFF3"
  mv "$DATA/annotation/$REFSEQ_GFF3.tmp" "$DATA/annotation/$REFSEQ_GFF3"
fi

# 1. genes -> CHM13 BED
args=(--annotation "$ANN" --flank "$FLANK" --data-dir "$DATA" --out "$W/batch")
[[ -n $LIST ]] && args+=(--gene-list "$LIST")
[[ -n $BIOTYPE ]] && args+=(--biotype "$BIOTYPE")
# ${arr[@]+...}: an empty array is "unbound" under set -u in bash 3.2 (macOS)
timed gene2bed python3 "$S/gene2bed.py" ${GENES[@]+"${GENES[@]}"} "${args[@]}"
[[ -s $W/batch.gene.bed ]] || { echo "[run] no genes found"; exit 1; }
awk -v OFS='\t' '{$4=$4"|gene"; print}' "$W/batch.gene.bed" > "$W/query.bed"
awk -v OFS='\t' '{$4=$4"|padded"; print}' "$W/batch.padded.bed" >> "$W/query.bed"

# 2. projection CHM13 -> haplotypes
case $MODE in
  remote)
    [[ -s $OFF.tbi ]] || { echo "[run] PAF offset index missing in $DATA: run ./run.sh --setup"; exit 1; }
    timed project python3 "$S/remote_project.py" --bed "$W/query.bed" --index "$OFF" --paf "$PAF_URL" -j "$JOBS" \
      > "$W/projection.paf";;
  impg)
    [[ -f $PAF.impg ]] || { echo "[run] impg index missing in $DATA: run ./run.sh --setup --full (or use --mode remote)"; exit 1; }
    timed project bash "$S/impg_project.sh" "$W/query.bed" "$W/projection.paf" "$PAF" "$JOBS";;
  local)
    [[ -s $PAF ]] || { echo "[run] local PAF missing in $DATA: run ./run.sh --setup --full"; exit 1; }
    timed project bash "$S/stream_project.sh" "$W/query.bed" "$W/projection.paf" "$PAF";;
  stream) timed project bash "$S/stream_project.sh" "$W/query.bed" "$W/projection.paf" "$PAF_URL";;
  *) echo "unknown mode $MODE"; exit 1;;
esac

# 3. hits, flags, fetch regions
timed summarize python3 "$S/summarize_projection.py" --paf "$W/projection.paf" --genes "$W/batch.genes.tsv" \
  --index "$IDX" --outdir "$OUT/genes" --regions-out "$W/regions.tsv" --flank "$FLANK" ${FASTA_DIR:+--fasta-dir "$FASTA_DIR"}
if ((OFFLINE)) && awk -F'\t' 'NR > 1 && $2 ~ /:\/\// {f = 1} END {exit !f}' "$W/regions.tsv"; then
  echo "[run] --offline: some assemblies are not in $FASTA_DIR:"
  awk -F'\t' 'NR > 1 && $2 ~ /:\/\// {print $1}' "$W/regions.tsv" | sort -u | head; exit 1
fi

# 4. sequence fetch (one job per assembly for the whole batch, parallel, retried, resumable)
timed fetch python3 "$S/fetch_seqs.py" --regions "$W/regions.tsv" --workdir "$W/fetch" --outdir "$OUT/genes" \
  -j "$JOBS" --method "$FETCH" \
  || echo "[run] WARNING: some assembly fetches failed (see $W/fetch/fetch.log); rerun the same command to retry"

# 5. per-gene summaries
rm -f "$W/batch_summary.tsv"
timed gene_summary python3 "$S/gene_summary.py" --outdir "$OUT/genes" --genes "$W/batch.genes.tsv" \
  --table "$W/batch_summary.tsv"
if [[ -n $GNU_TIME ]]; then python3 "$S/resources.py" "$LOG" > "$LOG/resources.tsv"; fi
column -t -s$'\t' "$W/batch_summary.tsv"
echo "[run] done: per-gene outputs in $OUT/genes/, step resources in $LOG/"
