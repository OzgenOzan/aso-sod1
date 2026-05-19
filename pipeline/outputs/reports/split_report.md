# Phase 5 -- Train/Test Split Report

**Generated:** 2026-05-13T20:17:31.863484

## Split Configuration

- **Split ratio:** 70% train / 30% test
- **Split unit:** sequence_cluster_id (cluster-aware)
- **Random seed:** 42
- **Stratification:** By median inhibition_percent bin

## Dataset Sizes

| Set | Rows | Clusters | % of Total |
|---|---|---|---|
| Train | 1533 | 280 | 71.1% |
| Test | 622 | 116 | 28.9% |
| Total | 2155 | 396 | 100% |

## Inhibition Distribution

| Statistic | Train | Test | Full |
|---|---|---|---|
| mean | 43.85 | 37.27 | 41.95 |
| std | 29.15 | 27.74 | 28.89 |
| median | 40.00 | 33.00 | 38.00 |
| min | 0.00 | 0.00 | 0.00 |
| max | 98.00 | 97.00 | 98.00 |

## Inhibition Bin Distribution

| Bin | Train | Test | Train % | Test % |
|---|---|---|---|---|
| 0-20 | 428 | 206 | 27.9% | 33.1% |
| 20-40 | 339 | 169 | 22.1% | 27.2% |
| 40-60 | 279 | 109 | 18.2% | 17.5% |
| 60-80 | 240 | 73 | 15.7% | 11.7% |
| 80-100 | 247 | 65 | 16.1% | 10.5% |

## Sequence Length Distribution

| Length | Train | Test |
|---|---|---|
| 16 | 180 | 52 |
| 17 | 934 | 320 |
| 18 | 116 | 52 |
| 20 | 303 | 198 |

## Leakage Control

- **Cluster overlap between train and test:** 0
- **Leakage status:** [OK] No leakage
- **Unique sequence overlap:** 0
