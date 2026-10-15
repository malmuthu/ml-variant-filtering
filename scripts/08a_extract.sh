#!/usr/bin/env bash

# Extract labeled SNP calls (with annotations) from the raw benchmark for the ML table.
# Label 1 = TP; 0 = FP 
# Run with gx-bio active

set -euo pipefail

B=results/bench/raw
OUT=results/ml
mkdir -p $OUT

FMT='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/BaseQRankSum\t%INFO/DP\t[%GT]\t[%GQ]\t[%AD]\n'

{
    printf "CHROM\tPOS\tREF\tALT\tQUAL\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tBaseQRankSum\tDP\tGT\tGQ\tAD\tlabel\n"
    bcftools view -v snps $B/tp.vcf.gz | bcftools query -f "$FMT" | awk '{print $0"\t1"}'
    bcftools view -v snps $B/fp.vcf.gz | bcftools query -f "$FMT" | awk '{print $0"\t0"}'
} > $OUT/snp_calls.tsv

bcftools query -f '%POS\n' results/calls/NA12878.chr20.raw.vcf.gz > $OUT/all_call_positions.txt

echo "Rows by label (1 = TP, 0 = FP):"
cut -f17 $OUT/snp_calls.tsv | tail -n +2 | sort | uniq -c
echo "All raw call positions: $(wc -l < $OUT/all_call_positions.txt)"