#!/usr/bin/env bash
# Validate extracted gene sequences: align a sample of them to the CHM13 gene sequence with minimap2
# and report identity and strand (expect >99% identity, '+' strand).
# Usage: validate.sh LOCUS [N_SAMPLE=10] [OUTDIR=./genefinder_out]
set -euo pipefail
LOCUS=$1; N=${2:-10}; OUT=${3:-$PWD/genefinder_out}
REPO=$(cd "$(dirname "$0")/.." && pwd)
[[ -x $REPO/tools/env/bin/minimap2 ]] && export PATH=$REPO/tools/env/bin:$PATH
D=$OUT/genes/$LOCUS
V=$D/validation; mkdir -p "$V"
# CHM13 control is the first record of gene.fa
awk '/^>/{n++} n==1' "$D/$LOCUS.gene.fa" > "$V/chm13.fa"
# random sample of HPRC primary sequences (+ GRCh38 and every reverse-strand one if few)
python3 - "$D/$LOCUS.gene.fa" "$V/sample.fa" "$N" <<'EOF'
import sys, random
src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs, name = [], None
for l in open(src):
    if l.startswith(">"):
        recs.append([l[1:].strip(), []])
    else:
        recs[-1][1].append(l.strip())
recs = [r for r in recs if "CHM13#" not in r[0] and " " not in r[0]]  # primary hits only
random.seed(1)
pick = [r for r in recs if "GRCh38#" in r[0]] + [r for r in recs if r[0].endswith("(-)")][:3]
rest = [r for r in recs if r not in pick]
pick += random.sample(rest, min(n, len(rest)))
with open(dst, "w") as fh:
    for h, s in pick:
        fh.write(f">{h}\n{''.join(s)}\n")
EOF
minimap2 -c -x asm20 --secondary=no "$V/chm13.fa" "$V/sample.fa" 2>/dev/null > "$V/sample_vs_chm13.paf"
{
  printf "query\tqlen\tq_aligned_frac\tstrand\tt_aligned_frac\tidentity\n"
  awk -F'\t' -v OFS='\t' '{printf "%s\t%d\t%.4f\t%s\t%.4f\t%.5f\n",$1,$2,($4-$3)/$2,$5,($9-$8)/$7,$10/$11}' "$V/sample_vs_chm13.paf"
} > "$V/validation.tsv"
column -t -s$'\t' "$V/validation.tsv"
awk -F'\t' 'NR>1{n++; if($4=="+") p++; if($6>0.99) h++} END{printf "summary: %d aligned, %d on + strand, %d with identity > 0.99\n", n, p, h}' "$V/validation.tsv" | tee "$V/validation_summary.txt"
