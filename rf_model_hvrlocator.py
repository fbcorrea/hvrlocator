#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import joblib
import argparse

"""
This script trains a Random Forest classifier to distinguish
sequences with primers vs. without, using QC metrics, and
saves the model as rf_model.pkl.
Usage:
  python rf_model_train.py \
      --no-primers summary_matrix_no_primers.tsv \
      --with-primers summary_matrix_withprimers.tsv \
      --output rf_model.pkl
"""

def main():
    parser = argparse.ArgumentParser(description="Train RF model for primer detection")
    parser.add_argument('--no-primers', required=True,
                        help='Path to summary matrix TSV for no-primer samples')
    parser.add_argument('--with-primers', required=True,
                        help='Path to summary matrix TSV for with-primer samples')
    parser.add_argument('--output', default='rf_model.pkl',
                        help='Output path for the trained model')
    args = parser.parse_args()

    # Load data
    df_no = pd.read_csv(args.no_primers, sep='\t')
    df_with = pd.read_csv(args.with_primers, sep='\t')

    # Label and combine
    df_no['label'] = 0
    df_with['label'] = 1
    df = pd.concat([df_no, df_with], ignore_index=True)

    # Feature selection: all numeric columns except sample ID and label
    feature_cols = [c for c in df.columns if c not in ['Sample_ID', 'label']]
    X = df[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0).values
    y = df['label'].values

    # Train Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X, y)

    # Save the model
    joblib.dump(rf, args.output)
    print(f"Model saved to: {args.output}")

if __name__ == '__main__':
    main()
