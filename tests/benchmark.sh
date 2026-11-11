#!/usr/bin/env bash
# Benchmark: pre-build times and one-gene query times, online vs offline, plus how the two results differ.
# Usage: tests/benchmark.sh [GENE=MALAT1] [JOBS=32] [WORKDIR=bench]
# Starts from an empty data dir in WORKDIR (so pre-build times are real) and writes WORKDIR/report.md.
set -euo pipefail
GENE=${1:-MALAT1}; J=${2:-32}; W=${3:-bench}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
GF="python3 $ROOT/genefinder.py --data $W/data -j $J"
rm -rf "$W"; mkdir -p "$W"
T="$W/times.tsv"; printf 'step\tseconds\n' > "$T"

step() {  # step NAME CMD...
    local name=$1; shift
    local t0; t0=$(date +%s.%N)
    "$@" 2> "$W/$name.log" || { cat "$W/$name.log"; exit 1; }
    printf '%s\t%.1f\n' "$name" "$(echo "$(date +%s.%N) - $t0" | bc)" >> "$T"
    tail -n 3 "$W/$name.log"
}

step setup_online   $GF setup online
step setup_offline  $GF setup offline
step get_online     $GF get "$GENE" --mode online  -o "$W/online"
step get_online_2   $GF get "$GENE" --mode online  -o "$W/online"
step get_offline    $GF get "$GENE" --mode offline -o "$W/offline"
step get_offline_2  $GF get "$GENE" --mode offline -o "$W/offline"

{
    echo "## genefinder benchmark: $GENE, -j $J, $(hostname -s), $(date -I)"
    echo
    echo "| step | seconds |"; echo "|---|---|"
    tail -n +2 "$T" | awk -F'\t' '{print "| "$1" | "$2" |"}'
    echo
    echo "Data kept after setup: online $(du -sh "$W/data/online" | cut -f1), offline $(du -sh "$W/data/offline" | cut -f1)," \
         "shared $(du -sh "$W/data/chm13_genes.tsv" | cut -f1)"
    echo "\`_2\` = the same query again (warm page cache / connection reuse)."
    echo
    python3 "$ROOT/tests/compare.py" "$W/online" "$W/offline" "$GENE"
} > "$W/report.md"
cat "$W/report.md"
