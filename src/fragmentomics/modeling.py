"""Assemble per-sample feature JSONs into a model-ready matrix and (optionally)
train a simple classifier.

This is the analysis layer that turns the fragmentomic features - including the
NuPEM coupling score - into a discriminative signal (e.g. cancer vs healthy). It
is deliberately small and dependency-light; swap in your own estimator/CV scheme
for real work.

Usage::

    python -m fragmentomics.modeling \\
        --features results/features/*.features.json \\
        --labels labels.csv \\
        --out model_report.json
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


def feature_vector(record: dict) -> dict[str, float]:
    """Flatten one sample's feature JSON into a named scalar feature dict."""
    f: dict[str, float] = {}
    L = record["lengths"]
    f["len_mean"] = L["mean"]
    f["len_median"] = L["median"]
    f["short_long_ratio"] = L["short_long_ratio"]
    f["frac_sub_nucleosomal"] = L["frac_sub_nucleosomal"]
    f["mds"] = record["motifs"]["mds"]
    f["nupem_coupling"] = record["nupem"]["coupling_score"]
    f["n_nucleosomes"] = record["n_nucleosomes"]
    return f


def build_matrix(paths: list[str]):
    rows, ids = [], []
    for p in paths:
        rec = json.loads(Path(p).read_text())
        ids.append(rec["sample_id"])
        rows.append(feature_vector(rec))
    names = sorted({k for r in rows for k in r})
    X = np.array([[r.get(n, np.nan) for n in names] for r in rows], dtype=float)
    return ids, names, X


def main() -> None:
    ap = argparse.ArgumentParser(description="Build feature matrix / train classifier.")
    ap.add_argument("--features", nargs="+", required=True, help="feature JSON files or globs")
    ap.add_argument("--labels", help="CSV: sample_id,label (enables training)")
    ap.add_argument("--out", default="model_report.json")
    args = ap.parse_args()

    paths: list[str] = []
    for pat in args.features:
        paths.extend(sorted(glob.glob(pat)) or [pat])

    ids, names, X = build_matrix(paths)
    report: dict = {"n_samples": len(ids), "feature_names": names, "sample_ids": ids}

    if args.labels:
        try:
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.model_selection import cross_val_score
        except ImportError:
            raise SystemExit("Install scikit-learn to enable training: pip install scikit-learn") from None

        label_map = {}
        for line in Path(args.labels).read_text().splitlines():
            if not line.strip() or line.startswith("sample_id"):
                continue
            sid, lab = line.split(",")[:2]
            label_map[sid.strip()] = lab.strip()
        y = np.array([label_map[i] for i in ids])

        # median-impute then 5-fold CV ROC-AUC
        col_median = np.nanmedian(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(col_median, inds[1])

        clf = RandomForestClassifier(n_estimators=400, random_state=0)
        scoring = "roc_auc" if len(set(y)) == 2 else "accuracy"
        scores = cross_val_score(clf, X, y, cv=5, scoring=scoring)
        clf.fit(X, y)
        report["cv_metric"] = scoring
        report["cv_scores"] = scores.tolist()
        report["cv_mean"] = float(scores.mean())
        report["feature_importances"] = dict(
            sorted(zip(names, clf.feature_importances_.tolist(), strict=True), key=lambda kv: kv[1], reverse=True)
        )

    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
