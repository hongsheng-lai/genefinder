#!/usr/bin/env bash
# Stream the HPRC R2 vs CHM13 PAF from S3 (or a local copy) and project BED regions
# without storing the PAF:  curl | bgzip -dc | awk (target chrom + overlap) | paf_project.py
# Usage: stream_project.sh regions.bed out.paf [PAF_URL_OR_PATH]
set -euo pipefail
BED=$1; OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/../config/sources.env"
SRC=${3:-$PAF_URL}
if [[ $SRC == http* ]]; then reader=(curl -sS --retry 5 "$SRC"); else reader=(cat "$SRC"); fi
if (( $(wc -l < "$BED") <= 50 )); then
  # few regions: a cheap awk prefilter (linear in #regions) keeps python off most lines
  filter=(awk -F'\t' -v OFS='\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next}
      { for(i=1;i<=n;i++) if($6==c[i] && $8<e[i] && $9>s[i]) {print; next} }' "$BED" -)
else
  filter=(cat)  # many regions: paf_project.py does a bisect overlap test per line
fi
"${reader[@]}" | bgzip -dc -@ 4 | "${filter[@]}" \
  | python3 "$HERE/paf_project.py" --bed "$BED" - > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
