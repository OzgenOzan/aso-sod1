"""
Phase 5 -- Train/Test Splitting
================================
70/30 cluster-aware stratified split.

Outputs:
  - train.csv
  - test.csv
  - split_report.md
"""

import pandas as pd
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    DATA_DIR, REPORT_DIR, RANDOM_SEED, TEST_SIZE,
    INHIBITION_BINS, INHIBITION_BIN_LABELS
)


def cluster_aware_stratified_split(df, test_size=TEST_SIZE, seed=RANDOM_SEED):
    """
    Split by sequence_cluster_id, not by individual row.
    Approximate stratification by inhibition_percent bins.
    """
    rng = np.random.RandomState(seed)

    # Assign inhibition bin to each cluster based on median inhibition
    cluster_stats = df.groupby("sequence_cluster_id").agg(
        median_inh=("inhibition_percent", "median"),
        n_rows=("inhibition_percent", "count"),
    ).reset_index()

    cluster_stats["inh_bin"] = pd.cut(
        cluster_stats["median_inh"],
        bins=INHIBITION_BINS,
        labels=INHIBITION_BIN_LABELS,
        include_lowest=True,
    )

    # Stratified split by bin
    test_clusters = []
    train_clusters = []

    for bin_label in INHIBITION_BIN_LABELS:
        bin_clusters = cluster_stats[cluster_stats["inh_bin"] == bin_label]["sequence_cluster_id"].values.copy()
        rng.shuffle(bin_clusters)

        # Calculate how many clusters to put in test
        n_test = max(1, int(len(bin_clusters) * test_size))
        test_clusters.extend(bin_clusters[:n_test])
        train_clusters.extend(bin_clusters[n_test:])

    test_clusters = set(test_clusters)
    train_clusters = set(train_clusters)

    # Verify no overlap
    assert len(test_clusters & train_clusters) == 0, "Cluster leakage detected!"

    df_train = df[df["sequence_cluster_id"].isin(train_clusters)].copy()
    df_test = df[df["sequence_cluster_id"].isin(test_clusters)].copy()

    return df_train, df_test


def report_split_balance(df_train, df_test, df_full):
    """Generate balance report for the split."""
    report = {}

    # Inhibition distribution
    for name, subset in [("train", df_train), ("test", df_test), ("full", df_full)]:
        report[f"{name}_n"] = len(subset)
        report[f"{name}_inh_mean"] = subset["inhibition_percent"].mean()
        report[f"{name}_inh_std"] = subset["inhibition_percent"].std()
        report[f"{name}_inh_median"] = subset["inhibition_percent"].median()

    # Per-category balance
    cat_cols_to_check = []
    for col in df_full.columns:
        if col.startswith("cell_") or col.startswith("transfection_") or col.startswith("primer_"):
            cat_cols_to_check.append(col)

    return report


def write_split_report(df_train, df_test, df_full):
    """Write split_report.md."""
    path = os.path.join(REPORT_DIR, "split_report.md")
    lines = []
    lines.append("# Phase 5 -- Train/Test Split Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Split Configuration\n\n")
    lines.append(f"- **Split ratio:** {1-TEST_SIZE:.0%} train / {TEST_SIZE:.0%} test\n")
    lines.append(f"- **Split unit:** sequence_cluster_id (cluster-aware)\n")
    lines.append(f"- **Random seed:** {RANDOM_SEED}\n")
    lines.append(f"- **Stratification:** By median inhibition_percent bin\n\n")

    lines.append("## Dataset Sizes\n\n")
    lines.append("| Set | Rows | Clusters | % of Total |\n|---|---|---|---|\n")
    n_train_clusters = df_train["sequence_cluster_id"].nunique()
    n_test_clusters = df_test["sequence_cluster_id"].nunique()
    n_total = len(df_full)
    lines.append(f"| Train | {len(df_train)} | {n_train_clusters} | {len(df_train)/n_total*100:.1f}% |\n")
    lines.append(f"| Test | {len(df_test)} | {n_test_clusters} | {len(df_test)/n_total*100:.1f}% |\n")
    lines.append(f"| Total | {n_total} | {n_train_clusters + n_test_clusters} | 100% |\n\n")

    # Inhibition distribution
    lines.append("## Inhibition Distribution\n\n")
    lines.append("| Statistic | Train | Test | Full |\n|---|---|---|---|\n")
    for stat in ["mean", "std", "median", "min", "max"]:
        train_val = getattr(df_train["inhibition_percent"], stat)()
        test_val = getattr(df_test["inhibition_percent"], stat)()
        full_val = getattr(df_full["inhibition_percent"], stat)()
        lines.append(f"| {stat} | {train_val:.2f} | {test_val:.2f} | {full_val:.2f} |\n")

    # Inhibition bin distribution
    lines.append("\n## Inhibition Bin Distribution\n\n")
    lines.append("| Bin | Train | Test | Train % | Test % |\n|---|---|---|---|---|\n")
    train_bins = pd.cut(df_train["inhibition_percent"], bins=INHIBITION_BINS, labels=INHIBITION_BIN_LABELS, include_lowest=True).value_counts().sort_index()
    test_bins = pd.cut(df_test["inhibition_percent"], bins=INHIBITION_BINS, labels=INHIBITION_BIN_LABELS, include_lowest=True).value_counts().sort_index()
    for b in INHIBITION_BIN_LABELS:
        tr = train_bins.get(b, 0)
        te = test_bins.get(b, 0)
        lines.append(f"| {b} | {tr} | {te} | {tr/len(df_train)*100:.1f}% | {te/len(df_test)*100:.1f}% |\n")

    # Sequence length balance
    lines.append("\n## Sequence Length Distribution\n\n")
    lines.append("| Length | Train | Test |\n|---|---|---|\n")
    for sl in sorted(df_full["seq_length"].unique()):
        tr = (df_train["seq_length"] == sl).sum()
        te = (df_test["seq_length"] == sl).sum()
        lines.append(f"| {sl} | {tr} | {te} |\n")

    # Leakage verification
    lines.append("\n## Leakage Control\n\n")
    train_clusters = set(df_train["sequence_cluster_id"].unique())
    test_clusters = set(df_test["sequence_cluster_id"].unique())
    overlap = train_clusters & test_clusters
    lines.append(f"- **Cluster overlap between train and test:** {len(overlap)}\n")
    lines.append(f"- **Leakage status:** {'[X] LEAKAGE DETECTED' if overlap else '[OK] No leakage'}\n")

    # Sequence overlap check
    train_seqs = set(df_train["sequence"].unique())
    test_seqs = set(df_test["sequence"].unique())
    seq_overlap = train_seqs & test_seqs
    lines.append(f"- **Unique sequence overlap:** {len(seq_overlap)}\n")
    if seq_overlap:
        lines.append(f"  - Note: These sequences share cluster IDs and should not cross partitions.\n")
        lines.append(f"  - This may occur if the same sequence was assigned different clusters.\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 5] Split report saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 5 -- TRAIN/TEST SPLITTING")
    print("=" * 70)

    # Load feature matrix
    fm_path = os.path.join(DATA_DIR, "feature_matrix.csv")
    df = pd.read_csv(fm_path)
    print(f"[Phase 5] Loaded feature matrix: {df.shape}")

    # Split
    df_train, df_test = cluster_aware_stratified_split(df)

    # Save
    train_path = os.path.join(DATA_DIR, "train.csv")
    test_path = os.path.join(DATA_DIR, "test.csv")
    df_train.to_csv(train_path, index=False)
    df_test.to_csv(test_path, index=False)
    print(f"[Phase 5] Train: {len(df_train)} rows -> {train_path}")
    print(f"[Phase 5] Test: {len(df_test)} rows -> {test_path}")

    # Report
    write_split_report(df_train, df_test, df)

    print(f"\n[Phase 5] COMPLETE [OK] -- Train: {len(df_train)}, Test: {len(df_test)}\n")
    return df_train, df_test


if __name__ == "__main__":
    main()
