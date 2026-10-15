#!/usr/bin/env python3

"""
Train and evaluate SNP call filters on chr20: QUAL alone, GATK hard filters, 
logistic regression and XGBoost. 30% of 1 Mb blocks are held out as the test set.

Run with gx-ml active
"""

import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

DATA = "results/ml/snp_features.csv.gz"
OUTDIR = "results/ml"
FIGDIR = "docs/figures"
SEED = 42
TEST_FRACTION = 0.3

FEATURES = [
    "QUAL", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum",
    "BaseQRankSum", "DP", "GQ", "alt_frac", "is_het", "is_transition",
    "is_multiallelic", "gc_100bp", "homopolymer_max", "nearest_call_dist",
]

# Data preparation
def prepare(df):
    """Log-scale skewed features"""
    X = df[FEATURES].copy()
    X["QUAL"] = np.log10(X["QUAL"] + 1)
    X["DP"] = np.log10(X["DP"] + 1)
    X["nearest_call_dist"] = np.log10(X["nearest_call_dist"].clip(1, 1e6))
    return X

def hard_filter_fail(df):
    return ((df["QD"] < 2.0) | (df["QUAL"] < 30.0) | (df["SOR"] > 3.0) |
            (df["FS"] > 60.0) | (df["MQ"] < 40.0) |
            (df["MQRankSum"] < -12.5) | (df["ReadPosRankSum"] < -8.0))

def make_models():
    return {
        "logistic": make_pipeline(
            SimpleImputer(strategy="median", add_indicator=True),
            StandardScaler(),
            LogisticRegression(max_iter=2000)),
        "xgboost": XGBClassifier(
            n_estimators=400, max_depth=4, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8,
            eval_metric="aucpr", random_state=SEED, n_jobs=8),
    }

# Metrics
def fp_ap(y, score_true):
    """Average precision for detecting FPs"""
    return average_precision_score(1 - y, -np.asarray(score_true))

def curve(y, score_true):
    """Keep calls from highest score down."""
    ys = y[np.argsort(-np.asarray(score_true), kind="stable")]
    tp, fp = np.cumsum(ys), np.cumsum(1 - ys)
    return tp / ys.sum(), tp / (tp + fp), fp

def at_recall(y, score, r):
    rec, prec, fpk = curve(y, score)
    i = np.argmax(rec >= r)
    return prec[i], int(fpk[i])

def at_precision(y, score, p):
    rec, prec, _ = curve(y, score)
    ok = prec >= p
    return float(rec[ok].max()) if ok.any() else 0.0

def main():
    os.makedirs(FIGDIR, exist_ok=True)
    df = pd.read_csv(DATA)
    y = df["label"].to_numpy()
    X = prepare(df)
    hf_fail = hard_filter_fail(df).to_numpy()

    print(f"Hard filters (Python), all labeled SNPs: remove {(hf_fail & (y == 0)).sum()} FP, "
          f"lose {(hf_fail & (y == 1)).sum()} TP")

    rng = np.random.default_rng(SEED)
    blocks = np.sort(df["block"].unique())
    test_blocks = rng.choice(blocks, size=int(round(TEST_FRACTION * len(blocks))), replace=False)
    is_test = df["block"].isin(test_blocks).to_numpy()
    Xtr, ytr, gtr = X[~is_test], y[~is_test], df["block"][~is_test]
    Xte, yte = X[is_test], y[is_test]
    print(f"Train: {len(ytr):,} calls ({(ytr == 0).sum()} FP) | "
          f"Test: {len(yte):,} calls ({(yte == 0).sum()} FP) in {len(test_blocks)} blocks")

    # Grouped cross-validation
    cv = {"QUAL": [], "logistic": [], "xgboost": []}
    for tr, va in GroupKFold(n_splits=5).split(Xtr, ytr, gtr):
        cv["QUAL"].append(fp_ap(ytr[va], Xtr["QUAL"].to_numpy()[va]))
        for name, m in make_models().items():
            m.fit(Xtr.iloc[tr], ytr[tr])
            cv[name].append(fp_ap(ytr[va], m.predict_proba(Xtr.iloc[va])[:, 1]))
    print("\nCross-validation, FP-detection average precision:")
    for name, vals in cv.items():
        print(f" {name:9s} {np.mean(vals):.3f} +/- {np.std(vals):.3f}")

    # Final models
    scores = {"QUAL": Xte["QUAL"].to_numpy()}
    models = make_models()
    for name, m in models.items():
        m.fit(Xtr, ytr)
        scores[name] = m.predict_proba(Xte)[:, 1]
    models["xgboost"].save_model(f"{OUTDIR}/xgb_model.json")

    n_tp = (yte == 1).sum()
    raw_prec = n_tp / len(yte)
    kept = ~hf_fail[is_test]
    hf_rec = (kept & (yte == 1)).sum() / n_tp
    hf_prec = (kept & (yte == 1)).sum() / kept.sum()
    hf_fp = int((kept & (yte == 0)).sum())

    print("\nTest set:")
    print(f" No filter precision {raw_prec:.4f} FP kept {(yte == 0).sum()}")
    print(f" Hard filters precision {hf_prec:.4f} FP kept {hf_fp} filter recall {hf_rec:.4f}")
    results = {"test_blocks": sorted(int(b) for b in test_blocks), "cv": cv,
               "hard_filter": {"recall": hf_rec, "precision": hf_prec, "fp_kept": hf_fp}}
    for name, s in scores.items():
        p_at, fp_at = at_recall(yte, s, hf_rec)
        r_at = at_precision(yte, s, hf_prec)
        ap = fp_ap(yte, s)
        print(f" {name:9s} FP-AP {ap:.3f} | at hard-filter recall: precision {p_at:.4f}, "
              f"FP kept {fp_at} | at hard-filter precision: recall {r_at:.4f}")
        results[name] = {"fp_ap": ap, "precision_at_hf_recall": p_at,
                         "fp_kept_at_hf_recall": fp_at, "recall_at_hf_precision": r_at}

    with open(f"{OUTDIR}/metrics_test.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    df.loc[is_test, ["POS", "block", "label"]].assign(
        hf_fail=hf_fail[is_test], **{f"score_{k}": v for k, v in scores.items()}
    ).to_csv(f"{OUTDIR}/test_predictions.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for name, s in scores.items():
        rec, prec, _ = curve(yte, s)
        ax.plot(rec, prec, label=name)
    ax.scatter([hf_rec], [hf_prec], color="black", zorder=3, label="GATK hard filters")
    ax.axhline(raw_prec, ls=":", color="gray", label="no filter")
    ax.set_xlim(max(0.0, hf_rec - 0.03), 1.0005)
    ax.set_ylim(raw_prec - 0.002, 1.0005)
    ax.set_xlabel("Filter recall (fraction of true SNP calls kept)")
    ax.set_ylabel("Precision")
    ax.set_title("SNP filters on held-out chr20 blocks")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/pr_test.png", dpi=150)
    print(f"\nSaved {FIGDIR}/pr_test.png")

if __name__ == "__main__":
    main()