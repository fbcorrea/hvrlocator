#!/usr/bin/env python3
import argparse, json, os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.model_selection import StratifiedShuffleSplit
import joblib
from datetime import datetime

"""
Train + evaluate a Random Forest to detect primer presence.

Inputs:
  --no-primers  TSV with features for no-primer runs
  --with-primers TSV with features for primer runs

Outputs (under --outdir, default: eval/):
  - metrics/classification_report.csv|json   <-- contains accuracy & recall
  - metrics/confusion_matrix.csv
  - predictions/test_predictions.csv
  - model/rf_model.pkl  (trained on full data or just train split if --no-refit)
  - model/feature_cols.json
"""

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def main():
    ap = argparse.ArgumentParser(description="Train + evaluate RF model for primer detection")
    ap.add_argument("--no-primers", required=True, help="TSV for no-primer samples")
    ap.add_argument("--with-primers", required=True, help="TSV for primer samples")
    ap.add_argument("--id-col", default="Sample_ID", help="Identifier column name")
    ap.add_argument("--outdir", default="eval", help="Output root directory")
    ap.add_argument("--model-out", default=None, help="Path to save model (overrides default under outdir)")
    ap.add_argument("--test-size", type=float, default=0.2, help="Test fraction (default 0.2)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--n-estimators", type=int, default=100, help="RF trees (default 100)")
    ap.add_argument("--no-refit", action="store_true",
                    help="Do NOT refit on full data after evaluation; keep train-split model")
    args = ap.parse_args()

    np.random.seed(args.seed)

    # Load
    df_no   = pd.read_csv(args.no_primers, sep="\t")
    df_with = pd.read_csv(args.with_primers, sep="\t")
    df_no["label"] = 0
    df_with["label"] = 1
    df = pd.concat([df_no, df_with], ignore_index=True)

    # Features = all numeric except id + label
    drop_cols = {args.id_col, "label"}
    feature_cols = [c for c in df.columns if c not in drop_cols]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values
    y = df["label"].values
    ids = df[args.id_col] if args.id_col in df.columns else pd.Series(np.arange(len(df)), name="ID")

    # Split (stratified 80/20)
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=args.test_size, random_state=args.seed)
    (train_idx, test_idx) = next(splitter.split(X, y))
    Xtr, Xte = X[train_idx], X[test_idx]
    ytr, yte = y[train_idx], y[test_idx]
    id_tr, id_te = ids.iloc[train_idx].astype(str).values, ids.iloc[test_idx].astype(str).values

    # Train
    rf = RandomForestClassifier(n_estimators=args.n_estimators, random_state=args.seed, n_jobs=-1)
    rf.fit(Xtr, ytr)

    # Evaluate
    ypred = rf.predict(Xte)
    acc = accuracy_score(yte, ypred)
    rep = classification_report(yte, ypred, target_names=["no-primer","primer"], output_dict=True, digits=4)
    cm  = confusion_matrix(yte, ypred, labels=[0,1])

    # Prepare output dirs
    metrics_dir = os.path.join(args.outdir, "metrics"); ensure_dir(metrics_dir)
    preds_dir   = os.path.join(args.outdir, "predictions"); ensure_dir(preds_dir)
    model_dir   = os.path.join(args.outdir, "model"); ensure_dir(model_dir)

    # Save predictions
    pred_df = pd.DataFrame({
        "id": id_te,
        "true_label": yte,
        "pred_label": ypred,
        "proba_primer": rf.predict_proba(Xte)[:,1]
    })
    pred_df.to_csv(os.path.join(preds_dir, "test_predictions.csv"), index=False)

    # Save metrics (JSON + CSV) and confusion matrix
    with open(os.path.join(metrics_dir, "classification_report.json"), "w") as f:
        json.dump({
            "generated_at": datetime.utcnow().isoformat()+"Z",
            "accuracy": acc,
            "report": rep
        }, f, indent=2)

    # Flatten to CSV
    rep_df = pd.DataFrame(rep).T
    rep_df.loc["accuracy","precision"] = ""
    rep_df.loc["accuracy","recall"] = acc
    rep_df.to_csv(os.path.join(metrics_dir, "classification_report.csv"))

    cm_df = pd.DataFrame(cm, index=["true:no-primer","true:primer"],
                            columns=["pred:no-primer","pred:primer"])
    cm_df.to_csv(os.path.join(metrics_dir, "confusion_matrix.csv"))

    # Optionally refit on full data for final model
    model_path = args.model_out or os.path.join(model_dir, "rf_model.pkl")
    if args.no_refit:
        joblib.dump({"model": rf, "feature_cols": feature_cols, "seed": args.seed}, model_path)
    else:
        rf_full = RandomForestClassifier(n_estimators=args.n_estimators, random_state=args.seed, n_jobs=-1)
        rf_full.fit(X, y)
        joblib.dump({"model": rf_full, "feature_cols": feature_cols, "seed": args.seed}, model_path)

    # Save feature columns
    with open(os.path.join(model_dir, "feature_cols.json"), "w") as f:
        json.dump({"feature_cols": feature_cols}, f, indent=2)

    print(f"[OK] Accuracy: {acc:.6f}")
    print(f"[OK] Metrics written to: {metrics_dir}")
    print(f"[OK] Predictions written to: {os.path.join(preds_dir, 'test_predictions.csv')}")
    print(f"[OK] Model saved to: {model_path}")

if __name__ == "__main__":
    main()
