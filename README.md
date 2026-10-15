# An ML Variant-Quality Classifier Benchmarked Against GIAB

## Learning to filter variant calls (NA12878, chr20)

A germline variant-calling pipeline built from raw reads, benchmarked against the Genome in a Bottle (GIAB) truth set, followed by a classifier that separates true from false SNP calls and is compared against GATK's recommended hard filters.

## Key results

- **Pipeline accuracy (raw calls, GIAB confident regions):** SNP F1 0.9935, non-SNP F1 0.9930.
- **Learned filter vs. hard filters (held-out regions):** at the same recall, the hard filters kept 24 false-positive SNPs; the learned models kept 10.
- **Simpler was enough:** logistic regression performed on par with XGBoost.
- **What drives predictions:** the five most important features were mapping quality (MQ), call confidence (QUAL), MQRankSum, distance to the nearest other call, and ReadPosRankSum.
- **Beyond hard filters:** about 29% of the model's signal comes from features the hard filters don't use: properties of the call itself (~15%) and reference context (~14%).

## Approach

| Step | Tools | Key decision |
|---|---|---|
| Reads | samtools | Streamed chr20 from the 1000 Genomes CRAM; converted back to FASTQ to realign from scratch |
| QC | FastQC, MultiQC | No trimming: clean reads, and BWA-MEM soft-clips poor read ends |
| Alignment | BWA-MEM, GATK MarkDuplicates | GRCh38 no-alt reference, which avoids mapping-quality loss from alternate contigs |
| Recalibration | GATK BaseRecalibrator | Measured first: qualities were already calibrated, so BQSR was not applied |
| Calling | GATK HaplotypeCaller | Single-sample mode |
| Benchmarking | RTG vcfeval | Haplotype-aware matching against GIAB, inside confident regions |
| Baseline | GATK VariantFiltration | GATK's recommended hard filters |
| Model | scikit-learn, XGBoost, SHAP | Train/test split by 1 Mb genomic blocks to prevent spatial leakage |

Indels were excluded because there were too few false positives to model.

## Data

All open access.

- **Reads:** 1000 Genomes 30x high-coverage, NA12878 
- **Reference:** GRCh38 no-alt analysis set 
- **Truth set:** GIAB HG001 v4.2.1 benchmark, GRCh38
- **Known sites:** GATK resource bundle hg38

## Reproducing

Conda environments with pinned versions are in `envs/`; scripts in `scripts/` run in order from the project root. The reference files are downloaded manually from NCBI before running `01_reference.sh`. Rerunning the pipeline from scratch reproduced byte-identical intermediate data.

## Limitations

- One sample and one chromosome; results may not generalize.
- GIAB confident regions exclude the hardest regions, so accuracy is optimistic.
- Hard filters used GATK's default thresholds; tuning them on this data would likely narrow the gap.

## Planned extensions

- Error analysis by genomic context (GIAB stratifications)
- Functional annotation of errors (SnpEff)
- Robustness to lower coverage (15x, 10x)


## References
- Byrska-Bishop, Marta, et al. "High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios." Cell 185.18 (2022): 3426-3440.
- Regier, Allison A., et al. "Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects." Nature communications 9.1 (2018): 4038.
- Fairley, Susan, et al. "The International Genome Sample Resource (IGSR) collection of open human genomic variation resources." Nucleic acids research 48.D1 (2020): D941-D947.
- Van der Auwera, Geraldine A., and Brian D. O'Connor. Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, 2020.
- Wagner, Justin, et al. "Benchmarking challenging small variants with linked and long reads." Cell genomics 2.5 (2022).
- Sayers, Eric W et al. “Database resources of the National Center for Biotechnology Information in 2025.” Nucleic acids research vol. 53,D1 (2025): D20-D29. doi:10.1093/nar/gkae979
- Krusche, Peter, et al. "Best practices for benchmarking germline small-variant calls in human genomes." Nature biotechnology 37.5 (2019): 555-560.

