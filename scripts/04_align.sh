#!/usr/bin/env bash

# Align NA12878 chr20 reads with BWA-MEM, mark duplicates, then sort and index.
# Run with gx-bio active

set -euo pipefail

REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
SRC=data/raw/NA12878.chr20.source.bam
IN=data/raw/NA12878.chr20
OUT=results/align
PREFIX=$OUT/NA12878.chr20
THREADS=10
mkdir -p $OUT

# All source read groups share one library
SM=$(samtools view -H $SRC | grep '^@RG' | tr '\t' '\n' | grep '^SM:' | sort -u | sed 's/^SM://')
LB=$(samtools view -H $SRC | grep '^@RG' | tr '\t' '\n' | grep '^LB:' | sort -u | sed 's/^LB://')
if [ "$(echo "$LB" | wc -l)" -ne 1 ]; then echo "Expected one library, found: $LB"; exit 1; fi
RG="@RG\tID:${SM}.chr20\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA\tPU:merged"
echo "Read group: $RG"

# Align step
# -K 100000000: fixed batch size, so output does not depend on thread count
# -Y: soft-clip supplementary alignments
# Output stays grouped by read name, which MarkDuplicates expects.
bwa mem -t $THREADS -K 100000000 -Y -R "$RG" $REF \
    ${IN}_R1.fastq.gz ${IN}_R2.fastq.gz 2> logs/04_bwa.log \
    | samtools view -b -o ${PREFIX}.unsorted.bam -

# Mark duplicates step
gatk --java-options "-Xmx8g" MarkDuplicates \
    -I ${PREFIX}.unsorted.bam \
    -O ${PREFIX}.markdup.bam \
    -M ${PREFIX}.markdup_metrics.txt \
    --ASSUME_SORT_ORDER queryname \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    2> logs/04_markdup.log

# Sort and index step
samtools sort -@ 4 -m 2G -o ${PREFIX}.bam ${PREFIX}.markdup.bam
samtools index ${PREFIX}.bam

# Clean up and summarize
rm ${PREFIX}.unsorted.bam ${PREFIX}.markdup.bam
samtools flagstat -@ 4 ${PREFIX}.bam > ${PREFIX}.flagstat.txt
cat ${PREFIX}.flagstat.txt