#!/usr/bin/env python3

"""
Build the SNP feature table

Run with gx-ml active
"""
import numpy as np
import pandas as pd
import pysam

REF = "data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CALLS = "results/ml/snp_calls.tsv"
ALLPOS = "results/ml/all_call_positions.txt"
OUT = "results/ml/snp_features.csv.gz"

BLOCK_SIZE = 1_000_000
GC_FLANK = 50
HP_FLANK = 10

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

FEATURES = [
    "QUAL", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum",
    "BaseQRankSum", "DP", "GQ", "alt_frac", "is_het", "is_transition",
    "is_multiallelic", "gc_100bp", "homopolymer_max", "nearest_call_dist",
]

def alt_fraction(ad):
    """Fraction of reads supporting any alternate allele, from AD"""
    if not isinstance(ad, str):
        return np.nan
    counts = [int(x) for x in ad.split(",") if x != "."]
    total = sum(counts)
    return (total - counts[0]) / total if total > 0 else np.nan

def gc_fraction(seq, i):
    """GC fraction of A/C/G/T bases within +/-GC_FLANK of 0-based position i"""
    w = seq[max(0, i - GC_FLANK): i + GC_FLANK + 1]
    acgt = sum(c in "ACGT" for c in w)
    return (w.count("G") + w.count("C")) / acgt if acgt else np.nan

def longest_homopolymer(seq, i):
    w = seq[max(0, i - HP_FLANK): i + HP_FLANK + 1]
    best = run = 1
    for a, b in zip(w, w[1:]):
        run = run + 1 if (a == b and a != "N") else 1
        best = max(best, run)
    return best

def main():
    df = pd.read_csv(CALLS, sep="\t", na_values=".", dtype={"GT": str, "AD": str, "ALT": str, "REF": str})
    print(f"Calls read: {len(df):,}")
    print(f"Duplicate positions: {df.duplicated(['CHROM', 'POS']).sum()}")

    gt = df["GT"].str.replace("|", "/", regex=False).str.split("/")
    df["is_het"] = gt.map(lambda a: int(len(set(a)) > 1))
    df["alt_frac"] = df["AD"].map(alt_fraction)
    first_alt = df["ALT"].str.split(",").str[0]
    df["is_transition"] = [int((r, a) in TRANSITIONS) for r, a in zip(df["REF"], first_alt)]
    df["is_multiallelic"] = df["ALT"].str.contains(",").astype(int)

    # Reference sequence context
    seq = pysam.FastaFile(REF).fetch("chr20").upper()
    idx = (df["POS"] - 1).to_numpy()
    df["gc_100bp"] = [gc_fraction(seq, i) for i in idx]
    df["homopolymer_max"] = [longest_homopolymer(seq, i) for i in idx]

    # Distance to the nearest other raw call
    allpos = np.sort(np.loadtxt(ALLPOS, dtype=np.int64))
    pos = df["POS"].to_numpy()
    li = np.searchsorted(allpos, pos, side="left") - 1
    ri = np.searchsorted(allpos, pos, side="right")
    left = np.where(li >= 0, pos - allpos[np.clip(li, 0, None)], np.inf)
    right = np.where(ri < len(allpos), allpos[np.clip(ri, None, len(allpos) - 1)] - pos, np.inf)
    df["nearest_call_dist"] = np.minimum(left, right)

    df["block"] = (df["POS"] - 1) // BLOCK_SIZE

    df.to_csv(OUT, index=False)
    print(f"Saved {OUT}")

    # Summary
    print("\nRows by label (1 = TP, 0 = FP):")
    print(df["label"].value_counts().to_string())
    print(f"\nBlocks: {df['block'].nunique()}, of which "
          f"{df.loc[df['label'] == 0, 'block'].nunique()} contain at least one FP")
    print("\nMissing values (% of rows):")
    print((df[FEATURES].isna().mean() * 100).round(1).to_string())
    print("\nMedian by label (0 = FP, 1 = TP):")
    print(df.groupby("label")[FEATURES].median().T.round(3).to_string())

if __name__ == "__main__":
    main()