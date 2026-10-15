#!/usr/bin/env bash

# Stream NA12878 chr20 from the 10000 Genomes 30x CRAM and convert to paired FASTQ.
set -euo pipefail

CRAM=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
REF=data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set
SRC=data/raw/NA12878.chr20.source.bam
OUT=data/raw/NA12878.chr20
mkdir -p data/raw

# Confirm chr20 in the CRAM matches the local reference
samtools view -H "$CRAM" | grep -w "SN:chr20"
grep -w "SN:chr20" ${REF}.dict

# Stream chr20
samtools view -@ 4 -b -T ${REF}.fna -o $SRC "$CRAM" chr20
samtools quickcheck $SRC 
samtools flagstat -@ 4 $SRC

# Convert back to raw paired FASTQ
samtools collate -@ 4 -u -O -T data/raw/collate_tmp $SRC \
    | samtools fastq -@ 4 -n -O -F 0x900 \
        -1 ${OUT}_R1.fastq.gz -2 ${OUT}_R2.fastq.gz \
        -s ${OUT}_singletons.fastq.gz -0 ${OUT}_other.fastq.gz -

# Read counts (R1 must equal R2)
for f in ${OUT}_R1 ${OUT}_R2 ${OUT}_singletons ${OUT}_other; do 
    echo "$(basename $f): $(( $(gunzip -c ${f}.fastq.gz | wc -l) / 4))"
done

