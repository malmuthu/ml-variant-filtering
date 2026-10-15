#!/usr/bin/env bash

# Alignment QC:summary statistics, insert size, and coverage
# Run with gx-bio

set -euo pipefail

BAM=results/align/NA12878.chr20.bam
OUT=results/align/NA12878.chr20

samtools stats -@ 4 $BAM > ${OUT}.stats.txt
grep "^SN" ${OUT}.stats.txt | cut -f2-3 \
    | grep -E "raw total|reads mapped:|reads MQ0|error rate|average length|insert size|properly paired|different chromosomes"

samtools coverage -r chr20 $BAM | tee ${OUT}.coverage.txt