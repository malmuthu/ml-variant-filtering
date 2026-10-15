#!/usr/bin/env bash

# Measure base-quality calibration with BaseRecalibrator (BQSR is not applied)
# Run with gx-bio active

set -euo pipefail

REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
BAM=results/align/NA12878.chr20.bam
K=data/known
TABLE=results/align/NA12878.chr20.recal.table

gatk --java-options "-Xmx8g" BaseRecalibrator \
    -R $REF -I $BAM -L chr20 \
    --known-sites $K/1000G_phase1.snps.high_confidence.hg38.chr20.vcf.gz \
    --known-sites $K/Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf.gz \
    --known-sites $K/Homo_sapiens_assembly38.known_indels.chr20.vcf.gz \
    -O $TABLE 2> logs/04d_bqsr_measure.log

echo "== Reported vs. empirical quality, base substitutions (M) =="
awk '/^#:GATKTable:RecalTable1/ {f=1; next} f && NF==0 {exit} f && ($1=="ReadGroup" || $3=="M")' $TABLE \
    | column -t
    