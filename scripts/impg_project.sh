#!/usr/bin/env bash
# Project a BED through the impg index with J concurrent impg processes.
# (impg 0.5.0 handles BED rows serially, ~0.1 s each, whatever -t is set to.)
# Usage: impg_project.sh query.bed out.paf PAF J
set -euo pipefail
BED=$1; OUT=$2; PAF=$3; J=${4:-16}
TMP=$OUT.parts; rm -rf "$TMP"; mkdir -p "$TMP"
n=$(wc -l < "$BED")
(( J > n )) && J=$n
# round-robin split balances long and short genes across chunks
awk -v J="$J" -v D="$TMP" '{print > (D "/chunk." (NR % J) ".bed")}' "$BED"
ls "$TMP"/chunk.*.bed | xargs -P "$J" -I{} sh -c \
  'impg query -a "$1" -i "$1.impg" -b "$2" --no-merge --min-transitive-len 1 -o paf -t 1 -v 0 > "$2.paf"' _ "$PAF" {}
cat "$TMP"/chunk.*.bed.paf > "$OUT.tmp"
# every chunk must have produced output (impg exits non-zero on error, but double-check)
for b in "$TMP"/chunk.*.bed; do [[ -f $b.paf ]] || { echo "missing output for $b" >&2; exit 1; }; done
mv "$OUT.tmp" "$OUT"; rm -rf "$TMP"
