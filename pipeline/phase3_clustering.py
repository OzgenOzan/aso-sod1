"""
Phase 3 -- Redundancy Removal / High-Identity Clustering
=========================================================
Clusters ASO sequences at >90% nucleotide identity using a fallback
pairwise comparison approach (external tools like CD-HIT-EST not available).

Outputs:
  - aso_sequences.fasta
  - sequence_clusters.tsv
  - redundancy_report.md
"""

import pandas as pd
import numpy as np
import os
import sys
from collections import defaultdict
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    DATA_DIR, REPORT_DIR, CLUSTER_IDENTITY_THRESHOLD, RANDOM_SEED
)


def load_cleaned_data():
    """Load cleaned dataset."""
    path = os.path.join(DATA_DIR, "cleaned_dataset.csv")
    df = pd.read_csv(path)
    print(f"[Phase 3] Loaded {len(df)} cleaned records")
    return df


def export_fasta(df):
    """Export unique ASO sequences to FASTA format."""
    unique_seqs = df.drop_duplicates(subset="sequence")[["sequence", "aso_group_id"]].copy()
    unique_seqs["fasta_id"] = [f"ASO_{i:04d}" for i in range(len(unique_seqs))]

    fasta_path = os.path.join(DATA_DIR, "aso_sequences.fasta")
    with open(fasta_path, "w") as f:
        for _, row in unique_seqs.iterrows():
            f.write(f">{row['fasta_id']}|{row['sequence']}\n")
            f.write(f"{row['sequence']}\n")

    print(f"[Phase 3] FASTA exported: {len(unique_seqs)} unique sequences -> {fasta_path}")
    return unique_seqs


def compute_pairwise_identity(seq1, seq2):
    """
    Compute nucleotide identity between two sequences.
    For equal-length sequences: identity = matches / length.
    For unequal-length sequences: use sliding window alignment.
    """
    if len(seq1) == len(seq2):
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / len(seq1)
    else:
        # Sliding window for unequal lengths
        shorter, longer = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)
        best_identity = 0.0
        for start in range(len(longer) - len(shorter) + 1):
            subseq = longer[start:start + len(shorter)]
            matches = sum(a == b for a, b in zip(shorter, subseq))
            identity = matches / len(shorter)
            best_identity = max(best_identity, identity)
        return best_identity


def cluster_sequences(unique_seqs, threshold=CLUSTER_IDENTITY_THRESHOLD):
    """
    Cluster sequences using pairwise identity and connected components.
    Sequences with identity > threshold are placed in the same cluster.
    """
    sequences = unique_seqs["sequence"].tolist()
    n = len(sequences)
    print(f"[Phase 3] Computing pairwise identity for {n} sequences...")

    # Group by length for efficient comparison
    length_groups = defaultdict(list)
    for i, seq in enumerate(sequences):
        length_groups[len(seq)].append(i)

    # Build adjacency list (connected components)
    adj = defaultdict(set)

    # Compare within same-length groups (most common case for ASOs)
    total_comparisons = 0
    high_identity_pairs = 0

    for length, indices in length_groups.items():
        for i, j in combinations(indices, 2):
            total_comparisons += 1
            identity = compute_pairwise_identity(sequences[i], sequences[j])
            if identity > threshold:
                adj[i].add(j)
                adj[j].add(i)
                high_identity_pairs += 1

    # Also compare across different lengths (for sequences differing by 1-3 nt)
    lengths = sorted(length_groups.keys())
    for li in range(len(lengths)):
        for lj in range(li + 1, len(lengths)):
            if abs(lengths[li] - lengths[lj]) > 3:
                continue  # Skip very different lengths
            for i in length_groups[lengths[li]]:
                for j in length_groups[lengths[lj]]:
                    total_comparisons += 1
                    identity = compute_pairwise_identity(sequences[i], sequences[j])
                    if identity > threshold:
                        adj[i].add(j)
                        adj[j].add(i)
                        high_identity_pairs += 1

    print(f"[Phase 3] Total comparisons: {total_comparisons}, High-identity pairs: {high_identity_pairs}")

    # Find connected components using BFS
    visited = set()
    clusters = {}
    cluster_id = 0

    for i in range(n):
        if i in visited:
            continue
        # BFS
        queue = [i]
        visited.add(i)
        component = []
        while queue:
            node = queue.pop(0)
            component.append(node)
            for neighbor in adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        for node in component:
            clusters[node] = cluster_id
        cluster_id += 1

    unique_seqs = unique_seqs.copy()
    unique_seqs["sequence_cluster_id"] = [clusters[i] for i in range(n)]

    n_clusters = len(set(clusters.values()))
    print(f"[Phase 3] Formed {n_clusters} clusters from {n} sequences")

    return unique_seqs, {
        "total_comparisons": total_comparisons,
        "high_identity_pairs": high_identity_pairs,
        "n_clusters": n_clusters,
        "n_sequences": n,
    }


def assign_clusters_to_dataset(df, unique_seqs):
    """Map sequence_cluster_id back to the full dataset."""
    seq_to_cluster = dict(zip(unique_seqs["sequence"], unique_seqs["sequence_cluster_id"]))
    df["sequence_cluster_id"] = df["sequence"].map(seq_to_cluster)

    # Verify no missing
    n_missing = df["sequence_cluster_id"].isna().sum()
    if n_missing > 0:
        print(f"[Phase 3] WARNING: {n_missing} rows have no cluster assignment")

    return df


def write_cluster_outputs(unique_seqs, df, stats):
    """Write sequence_clusters.tsv and redundancy_report.md."""
    # TSV
    tsv_path = os.path.join(DATA_DIR, "sequence_clusters.tsv")
    out = unique_seqs[["fasta_id", "sequence", "sequence_cluster_id"]].copy()
    out.to_csv(tsv_path, sep="\t", index=False)
    print(f"[Phase 3] Cluster assignments saved -> {tsv_path}")

    # Report
    report_path = os.path.join(REPORT_DIR, "redundancy_report.md")
    lines = []
    lines.append("# Phase 3 -- Redundancy Removal / High-Identity Clustering Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Method\n\n")
    lines.append("**Fallback pairwise comparison** was used because CD-HIT-EST, UCLUST, and VSEARCH\n")
    lines.append("are not available in the current environment.\n\n")
    lines.append("### Algorithm\n")
    lines.append(f"1. Identity threshold: **>{CLUSTER_IDENTITY_THRESHOLD * 100:.0f}%**\n")
    lines.append("2. For equal-length sequences: identity = matches / length\n")
    lines.append("3. For unequal-length sequences: sliding window alignment, best identity\n")
    lines.append("4. Sequences with identity > threshold are connected.\n")
    lines.append("5. Connected components (BFS) define clusters.\n\n")

    lines.append("## Results\n\n")
    lines.append(f"- **Total unique sequences:** {stats['n_sequences']}\n")
    lines.append(f"- **Total pairwise comparisons:** {stats['total_comparisons']}\n")
    lines.append(f"- **High-identity pairs (>{CLUSTER_IDENTITY_THRESHOLD*100:.0f}%):** {stats['high_identity_pairs']}\n")
    lines.append(f"- **Number of clusters:** {stats['n_clusters']}\n\n")

    # Cluster size distribution
    cluster_sizes = unique_seqs["sequence_cluster_id"].value_counts()
    lines.append("## Cluster Size Distribution\n\n")
    lines.append("| Cluster Size | Count |\n|---|---|\n")
    for size, count in sorted(cluster_sizes.value_counts().items()):
        lines.append(f"| {size} | {count} |\n")

    lines.append(f"\n- **Singletons (size=1):** {(cluster_sizes == 1).sum()}\n")
    lines.append(f"- **Largest cluster size:** {cluster_sizes.max()}\n")
    lines.append(f"- **Mean cluster size:** {cluster_sizes.mean():.2f}\n")

    lines.append("\n## Impact on Dataset\n\n")
    lines.append(f"- **Dataset rows with cluster assignment:** {len(df)}\n")
    lines.append(f"- **Unique clusters in dataset:** {df['sequence_cluster_id'].nunique()}\n")

    lines.append("\n## Note on Short ASOs\n\n")
    lines.append("- For a 20-mer: 1 mismatch = 95% identity, 2 mismatches = 90% identity.\n")
    lines.append("- For a 17-mer: 1 mismatch = 94.1% identity, 2 mismatches = 88.2% identity.\n")
    lines.append("- The 90% threshold therefore groups sequences differing by ~2 nucleotides for 20-mers.\n")

    with open(report_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 3] Redundancy report saved -> {report_path}")


def main():
    print("=" * 70)
    print("PHASE 3 -- REDUNDANCY REMOVAL / HIGH-IDENTITY CLUSTERING")
    print("=" * 70)

    df = load_cleaned_data()

    # Export FASTA
    unique_seqs = export_fasta(df)

    # Cluster sequences
    unique_seqs, stats = cluster_sequences(unique_seqs)

    # Assign clusters back to dataset
    df = assign_clusters_to_dataset(df, unique_seqs)

    # Save updated dataset
    df.to_csv(os.path.join(DATA_DIR, "cleaned_dataset_clustered.csv"), index=False)

    # Write outputs
    write_cluster_outputs(unique_seqs, df, stats)

    print(f"\n[Phase 3] COMPLETE [OK] -- {stats['n_clusters']} clusters from {stats['n_sequences']} sequences\n")
    return df, unique_seqs


if __name__ == "__main__":
    main()
