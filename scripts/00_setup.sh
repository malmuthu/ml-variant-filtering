#!/usr/bin/env/ bash

#Build conda environments and record tool versions

set -euo pipefail

conda env create -f envs/bio.yml
conda env create -f envs/annot.yml
conda env create -f envs/ml.yml

{ echo "## gx-bio"
    conda run -n gx-bio samtools --version | sed -n '1p'
    conda run -n gx-bio bcftools --version | sed -n '1p'
    conda run -n gx-bio bwa 2>&1 | grep Version || true
    conda run -n gx-bio fastqc --version 
    conda run -n gx-bio multiqc --version
    conda run -n gx-bio gatk --version 2>&1 | grep -E "GATK\)|HTSJDK|Picard"
    conda run -n gx-bio rtg version | sed -n '1p'
    conda run -n gx-bio bedtools --version 
    conda run -n gx-bio java -version 2>&1 | sed -n '1p'
    echo "## gx-annot"
    conda run -n gx-annot snpEff -version 2>&1 | sed -n '1p'
    conda run -n gx-annot java -version 2>&1 | sed -n '1p'
    echo "## gx-ml"
    conda run -n gx-ml python -c "
import pandas, numpy, sklearn, xgboost, shap, matplotlib, seaborn, pysam
for m in [pandas, numpy, sklearn, xgboost, shap, matplotlib, seaborn, pysam]:
    print(m.__name__, m.__version__)"
} > versions.txt