#!/usr/bin/env bash

# Benchmark raw chr20 calls against GIAB HG001 V4.2.1 with RTG vcfeval
# Run with gx-bio active

set -euo pipefail

REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set
SDF=${REF}.sdf
T=data/truth
TRUTH=$T/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
BED=$T/HG001_GRCh38_1_22_v4.2.1_benchmark.chr20.bed
CALLS=results/calls/NA12878.chr20.raw.vcf.gz
OUT=results/bench/raw

# One-time conversion of the reference to RTG's SDF format
[ -d $SDF ] || rtg format -o $SDF ${REF}.fna > logs/06_rtg_format.log 2>&1

# Confident regions on chr20, and how much of chr20 they cover
awk '$1=="chr20"' $T/HG001_GRCh38_1_22_v4.2.1_benchmark.bed > $BED
awk '{s+=$3-$2} END {printf "Confident regions: %.1f Mb (%.1f%% of chr20)\n", s/1e6, 100*s/64444167}' $BED

# vcfeval cannot overwrite an existing output folder
mkdir -p results/bench
rm -rf $OUT
rtg vcfeval -b $TRUTH -c $CALLS -t $SDF -o $OUT \
    --evaluation-regions $BED --vcf-score-field QUAL -T 8 \
    > logs/06_vcfeval.log 2>&1

echo "== vcfeval summary (all variants) =="
cat $OUT/summary.txt

echo "== By type =="
for type in SNP non-SNP; do
    if [ $type = SNP ]; then sel="-v snps"; else sel="-V snps"; fi
    tpq=$(bcftools view -H $sel $OUT/tp.vcf.gz | wc -l)
    tpb=$(bcftools view -H $sel $OUT/tp-baseline.vcf.gz | wc -l)
    fp=$(bcftools view -H $sel $OUT/fp.vcf.gz | wc -l)
    fn=$(bcftools view -H $sel $OUT/fn.vcf.gz | wc -l)
    awk -v t="$type" -v tpq="$tpq" -v tpb="$tpb" -v fp="$fp" -v fn="$fn" 'BEGIN {
        p=tpq/(tpq+fp); r=tpb/(tpb+fn)
        printf "%-8s TP=%d FP=%d FN=%d precision=%.4f recall=%.4f F1=%.4f\n", t, tpq, fp, fn, p, r, 2*p*r/(p+r) }'
done