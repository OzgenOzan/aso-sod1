# Phase 8 -- Model Validation Report

**Generated:** 2026-05-13T20:17:55.123877

## 1. Held-Out Test Set Performance

| Metric | Value |
|---|---|
| AUPRC | 0.8391 |
| AUROC | 0.9461 |
| Calibration_intercept | 3.0131 |
| Calibration_slope | 1.0153 |
| Confusion_FN | 58 |
| Confusion_FP | 6 |
| Confusion_TN | 509 |
| Confusion_TP | 49 |
| F1 | 0.6049 |
| MAE | 12.4168 |
| MedianAE | 10.1673 |
| Pearson_r | 0.8301 |
| R2 | 0.6726 |
| RMSE | 15.8569 |
| Sensitivity | 0.4579 |
| Spearman_r | 0.8025 |
| Specificity | 0.9883 |
| Within_10 | 48.3923 |
| Within_15 | 68.0064 |
| Within_5 | 26.5273 |

## 2. Internal Cross-Validation

| Metric | Mean +/- Std |
|---|---|
| MAE | 16.6285 +/- 6.6688 |
| RMSE | 21.0304 +/- 7.3087 |
| R2 | 0.3691 +/- 0.3600 |
| Pearson_r | 0.7426 +/- 0.0603 |
| Spearman_r | 0.7310 +/- 0.0572 |

## 4. Interpretation Notes

- Feature importance and SHAP plots are available in the figures/ directory.
- Applicability domain warnings are included with predictions.
- Y-scrambling test recommended for final validation (computationally expensive).
