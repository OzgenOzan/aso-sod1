# Phase 7 -- Model Training Report

**Generated:** 2026-05-13T20:17:48.719823

## Training Configuration

- **Random seed:** 42
- **Target variable:** inhibition_percent (bounded regression, 0-100)
- **Predictions clipped to [0, 100]**
- **Cross-validation:** Cluster-aware GroupKFold

## Model Comparison (Test Set)

| Model | MAE | RMSE | R2 | Pearson_r | Spearman_r | Within_5 | Within_10 | Within_15 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mean Predictor | 24.5613 | 28.4828 | -0.0563 | nan | nan | 9.8071 | 21.2219 | 31.0289 |
| Linear Regression | 16.0957 | 20.6722 | 0.4436 | 0.7004 | 0.6939 | 24.5981 | 40.6752 | 54.0193 |
| Ridge | 15.9506 | 20.3009 | 0.4634 | 0.7080 | 0.7000 | 23.4727 | 42.2830 | 54.3408 |
| Elastic Net | 17.8780 | 22.2916 | 0.3530 | 0.6205 | 0.6131 | 19.6141 | 33.1190 | 49.1961 |
| Random Forest | 17.1999 | 20.9417 | 0.4290 | 0.7030 | 0.7052 | 15.9164 | 33.6013 | 47.7492 |
| Gradient Boosting | 12.6092 | 16.6020 | 0.6411 | 0.8201 | 0.7954 | 28.7781 | 50.3215 | 67.6849 |
| XGBoost | 12.3547 | 16.0018 | 0.6666 | 0.8250 | 0.7921 | 26.5273 | 50.8039 | 68.9711 |
| LightGBM | 12.4168 | 15.8569 | 0.6726 | 0.8301 | 0.8025 | 26.5273 | 48.3923 | 68.0064 |
| SVR | 19.3111 | 24.4433 | 0.2221 | 0.5352 | 0.5874 | 15.7556 | 31.9936 | 49.0354 |
| MLP (PyTorch) | nan | nan | nan | nan | nan | nan | nan | nan |

## Best Model: **LightGBM**

- R^2 = 0.6726
- MAE = 12.4168
- RMSE = 15.8569

## Models Trained

1. **Mean Predictor** -- Baseline (predicts training mean)
2. **Linear Regression** -- OLS on scaled features
3. **Ridge** -- L2-regularized linear regression (α=10)
4. **Elastic Net** -- L1+L2 regularization (α=1.0, l1_ratio=0.5)
5. **Random Forest** -- 500 trees, max_depth=15
6. **Gradient Boosting** -- 500 trees, max_depth=5, lr=0.05
7. **XGBoost** -- 500 trees, max_depth=6, lr=0.05
8. **LightGBM** -- 500 trees, max_depth=6, lr=0.05
9. **SVR** -- RBF kernel, C=10, ε=2.0
10. **MLP (PyTorch)** -- 3-layer (256->128->64->1), ReLU, BatchNorm, Dropout
