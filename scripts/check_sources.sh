#!/usr/bin/env bash
# Check that every remote input of genefinder is reachable, and show its size / ETag / Last-Modified,
# so you can compare against the values recorded in config/sources.env and docs/DATA_SOURCES.md.
# Also checks, for a sample of assemblies from the Release 2 index, that FASTA + .fai + .gzi exist.
# Usage: scripts/check_sources.sh [N_ASSEMBLIES_TO_SAMPLE=5]   (no downloads beyond the small index CSV)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/../config/sources.env"
source "$HERE/compat.sh"
N=${1:-5}

head_info() {  # url -> status, size, etag, last-modified
  curl -sSI --retry 3 "$1" | tr -d '\r' | awk -v u="$1" '
    /^HTTP/ {st=$2} tolower($1)=="content-length:" {sz=$2} tolower($1)=="etag:" {et=$2}
    tolower($1)=="last-modified:" {sub(/^[^:]*: /,""); lm=$0}
    END {printf "%-4s %15s  %-40s %-30s %s\n", st, sz, et, lm, u}'
}

printf "%-4s %15s  %-40s %-30s %s\n" HTTP bytes ETag Last-Modified URL
for u in "$PAF_URL" "$ASSEMBLY_INDEX_URL" "$ANNOTATION_BASE_URL/$GENCODE_GFF3" "$ANNOTATION_BASE_URL/$GENCODE_GFF3.tbi" \
         "$ANNOTATION_BASE_URL/$REFSEQ_GFF3" "$CHM13_FASTA_URL" "$CHM13_FASTA_URL.fai" "$CHM13_FASTA_URL.gzi"; do
  head_info "$u"
done

echo
echo "== assembly index: $ASSEMBLY_INDEX_URL"
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
curl -sS --fail --retry 3 -o "$tmp" "$ASSEMBLY_INDEX_URL"
echo "md5 $(md5_of "$tmp")  (expected 90c9e66a58580484facd59d3ed28c52c for the pinned commit)"
python3 - "$tmp" "$N" "$HPRC_BUCKET" <<'EOF' | while read -r u; do for ext in "" .fai .gzi; do head_info "$u$ext"; done; done
import csv, random, sys
rows = list(csv.DictReader(open(sys.argv[1])))
print(f"{len(rows)} assemblies listed; sampling {sys.argv[2]}", file=sys.stderr)
bucket = sys.argv[3].rstrip("/") + "/"
for r in random.sample(rows, min(int(sys.argv[2]), len(rows))):
    print(r["assembly"].replace("s3://human-pangenomics/", bucket))
EOF
