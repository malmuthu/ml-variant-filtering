#!/usr/bin/env bash

set -euo pipefail

mkdir -p results/qc 
fastqc -t 4 -o results/qc data/raw/NA12878.chr20_R1.fastq.gz data/raw/NA12878.chr20_R2.fastq.gz
(cd results/qc && multiqc . -f)