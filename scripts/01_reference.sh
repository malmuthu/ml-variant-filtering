#!/usr/bin/env bash

# Verify and prepare the GRCh38 no-alt analysis set.
# Expects the four files listed in the README to have been downloaded manually or using script into data/ref/

set -euo pipefail

REF=GCA_000001405.15_GRCh38_no_alt_analysis_set
cd data/ref

# Verify checksums (`md5 -q` on macOS only)
for f in ${REF}.fna.fai ${REF}.fna.gz ${REF}.fna.bwa_index.tar.gz; do
    expected=$(grep "/${f}$" md5checksums.txt | awk '{print $1}')
    actual=$(md5 -q "$f")
    [ "$expected" = "$actual" ] && echo "OK   $f" || { echo "FAIL $f"; exit 1; }
done

tar -xzf ${REF}.fna.bwa_index.tar.gz && rm ${REF}.fna.bwa_index.tar.gz
gunzip ${REF}.fna.gz
samtools dict ${REF}.fna -o ${REF}.dict