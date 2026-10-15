#!/usr/bin/env bash

# Apply GATK's recommended hard filters and benchmark the calls that pass
# SNP thresholds apply to SNPs; indel thresholds apply to all non-SNP records.
# Run with gx-bio active

set -euo pipefail

REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set
RAW=results/calls/NA12878.chr20.raw.vcf.gz
HF=results/calls/NA12878.chr20.hardfilter.vcf.gz
TRUTH=data/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
BED=data/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.chr20.bed
OUT=results/bench/hardfilter

gatk VariantFiltration -R ${REF}.fna -V $RAW -O $HF \
    --filter-name SNP_QD2 --filter-expression "vc.isSNP() && QD < 2.0" \
    --filter-name SNP_QUAL30 --filter-expression "vc.isSNP() && QUAL < 30.0" \
    --filter-name SNP_SOR3 --filter-expression "vc.isSNP() && SOR > 3.0" \
    --filter-name SNP_FS60 --filter-expression "vc.isSNP() && FS > 60.0" \
    --filter-name SNP_MQ40 --filter-expression "vc.isSNP() && MQ < 40.0" \
    --filter-name SNP_MQRS --filter-expression "vc.isSNP() && MQRankSum < -12.5" \
    --filter-name SNP_RPRS --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" \
    --filter-name IND_QD2 --filter-expression "!vc.isSNP() && QD < 2.0" \
    --filter-name IND_QUAL30 --filter-expression "!vc.isSNP() && QUAL < 30.0" \
    --filter-name IND_FS200 --filter-expression "!vc.isSNP() && FS > 200.0" \
    --filter-name IND_RPRS --filter-expression "!vc.isSNP() && ReadPosRankSum < -20.0" \
    2> logs/07_variantfiltration.log

echo "== Calls per FILTER value (all of chr20) =="
bcftools query -f '%FILTER\n' $HF | sort | uniq -c | sort -k1,1nr | sed -n '1,12p'

rm -rf $OUT
rtg vcfeval -b $TRUTH -c $HF -t ${REF}.sdf -o $OUT \
    --evaluation-regions $BED --vcf-score-field QUAL -T 8 \
    > logs/07_vcfeval.log 2>&1

metrics() {
    local dir=$1
    for type in SNP non-SNP; do
    if [ $type = SNP ]; then sel="-v snps"; else sel="-V snps"; fi
    tpq=$(bcftools view -H $sel $dir/tp.vcf.gz | wc -l)
    tpb=$(bcftools view -H $sel $dir/tp-baseline.vcf.gz | wc -l)
    fp=$(bcftools view -H $sel $dir/fp.vcf.gz | wc -l)
    fn=$(bcftools view -H $sel $dir/fn.vcf.gz | wc -l)
    awk -v t="$type" -v tpq="$tpq" -v tpb="$tpb" -v fp="$fp" -v fn="$fn" 'BEGIN {
        p=tpq/(tpq+fp); r=tpb/(tpb+fn)
        printf "%-8s TP=%d FP=%d FN=%d precision=%.4f recall=%.4f F1=%.4f\n", t, tpq, fp, fn, p, r, 2*p*r/(p+r) }'
    done
}

echo "== Raw calls =="
metrics results/bench/raw
echo "== Hard-filtered (PASS only) =="
metrics $OUT