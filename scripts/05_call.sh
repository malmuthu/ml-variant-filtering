#!/usr/bin/env bash

# Call germline variants on chr20 with GATK HaplotypeCaller (single-sample VCF mode).
# Run with gx-bio active

set -euo pipefail

REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
BAM=results/align/NA12878.chr20.bam
OUT=results/calls
VCF=$OUT/NA12878.chr20.raw.vcf.gz
mkdir -p $OUT

gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R $REF -I $BAM -L chr20 \
    -O $VCF \
    2> logs/05_haplotypecaller.log

# Summary statistics
bcftools stats -s - $VCF > $OUT/NA12878.chr20.raw.stats.txt
echo "== Counts =="
grep "^SN" $OUT/NA12878.chr20.raw.stats.txt | cut -f3-4
echo "== Ti/Tv =="
grep "^TSTV" $OUT/NA12878.chr20.raw.stats.txt | cut -f3-5
echo "== Genotypes (hom-ref, hom-alt, het) =="
grep "^PSC" $OUT/NA12878.chr20.raw.stats.txt | cut -f4-6