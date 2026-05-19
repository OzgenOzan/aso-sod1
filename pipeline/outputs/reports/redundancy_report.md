# Phase 3 -- Redundancy Removal / High-Identity Clustering Report

**Generated:** 2026-05-13T20:17:26.452229

## Method

**Fallback pairwise comparison** was used because CD-HIT-EST, UCLUST, and VSEARCH
are not available in the current environment.

### Algorithm
1. Identity threshold: **>90%**
2. For equal-length sequences: identity = matches / length
3. For unequal-length sequences: sliding window alignment, best identity
4. Sequences with identity > threshold are connected.
5. Connected components (BFS) define clusters.

## Results

- **Total unique sequences:** 688
- **Total pairwise comparisons:** 230828
- **High-identity pairs (>90%):** 415
- **Number of clusters:** 396

## Cluster Size Distribution

| Cluster Size | Count |
|---|---|
| 1 | 341 |
| 2 | 6 |
| 3 | 18 |
| 4 | 2 |
| 5 | 16 |
| 6 | 3 |
| 7 | 1 |
| 8 | 2 |
| 14 | 1 |
| 17 | 2 |
| 20 | 1 |
| 21 | 1 |
| 22 | 1 |
| 41 | 1 |

- **Singletons (size=1):** 341
- **Largest cluster size:** 41
- **Mean cluster size:** 1.74

## Impact on Dataset

- **Dataset rows with cluster assignment:** 2155
- **Unique clusters in dataset:** 396

## Note on Short ASOs

- For a 20-mer: 1 mismatch = 95% identity, 2 mismatches = 90% identity.
- For a 17-mer: 1 mismatch = 94.1% identity, 2 mismatches = 88.2% identity.
- The 90% threshold therefore groups sequences differing by ~2 nucleotides for 20-mers.
