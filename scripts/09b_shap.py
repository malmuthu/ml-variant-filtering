#!/usr/bin/env python3

"""
Explain the XGBoost SNP filter on the held-out test blocks with SHAP

Run with gx-ml active
"""

import importlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shap
import xgboost as xgb

sys.path.insert(0, "scripts")
train = importlib.import_module("09_train")

DATA = "results/ml/snp_features.csv.gz"
MODEL = "results/ml/xgb_model.json"
METRICS = "results/ml/metrics_test.json"
FIGDIR = "docs/figures"

GROUPS = {
    "caller annotations": ["QUAL", "QD", "FS", "SOR", "MQ", "MQRankSum",
                           "ReadPosRankSum", "BaseQRankSum", "DP", "GQ"],
    "derived from call": ["alt_frac", "is_het", "is_transition", "is_multiallelic"],
    "reference context": ["gc_100bp", "homopolymer_max", "nearest_call_dist"],
}

def main():
    df = pd.read_csv(DATA)
    with open(METRICS) as f:
        test_blocks = json.load(f)["test_blocks"]
    test = df[df["block"].isin(test_blocks)]
    Xte = train.prepare(test)
    yte = test["label"].to_numpy()

    booster = xgb.Booster()
    booster.load_model(MODEL)
    contribs = booster.predict(xgb.DMatrix(Xte, missing=np.nan), pred_contribs=True)
    sv = contribs[:, :-1]

    # Global importance: mean |SHAP| per feature
    imp = pd.Series(np.abs(sv).mean(axis=0), index=Xte.columns).sort_values(ascending=False)
    print("Mean |SHAP| per feature (log-odds), test set:")
    print(imp.round(3).to_string())

    # Share of total importance by feature group
    print("\nShare of total importance by group:")
    for g, cols in GROUPS.items():
        print(f" {g:20s} {100 * imp[cols].sum() / imp.sum():5.1f}%")

    fp_mean = pd.Series(sv[yte == 0].mean(axis=0), index=Xte.columns).sort_values()
    print("\nMean SHAP among false positives:")
    print(fp_mean.head(6).round(3).to_string())

    # Beeswarm plot
    shap.summary_plot(sv, Xte, max_display=len(Xte.columns), show=False)
    plt.title("SHAP values, XGBoost SNP filter (test blocks)")
    plt.savefig(f"{FIGDIR}/shap_beeswarm.png", dpi=150, bbox_inches="tight")
    print(f"\nSaved {FIGDIR}/shap_beeswarm.png")

if __name__ == "__main__":
    main()