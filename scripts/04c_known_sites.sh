#!/usr/bin/env bash

# Dtream chr20 of the BQSR known-sites resources from the GATK resource bundle.
# Run with gx-bio active.

set -euo pipefail

BUNDLE=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
mkdir -p data/known && cd data/known

for f in 1000G_phase1.snps.high_confidence.hg38 \
        Mills_and_1000G_gold_standard.indels.hg38 \
        Homo_sapiens_assembly38.known_indels; do 
    bcftools view -r chr20 -Oz -o ${f}.chr20.vcf.gz "${BUNDLE}/${f}.vcf.gz"
    bcftools index -t ${f}.chr20.vcf.gz
    rm -f ${f}.vcf.gz.tbi
    echo "${f}: $(bcftools view -H ${f}.chr20.vcf.gz | wc -l) sites"
done
