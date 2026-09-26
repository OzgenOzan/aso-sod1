"""
Phase 7 -- Modeling
====================
Train baseline and advanced models for inhibition_percent prediction.

Remediation (ASO1-1 / ASO1-4):
  - Model selection is now performed with cluster-aware GroupKFold
    cross-validation on the TRAINING split only (grouped by
    sequence_cluster_id persisted by phase5 in train.csv).
  - The held-out test set is touched exactly ONCE, for the final
    evaluation of the already-selected model ("final held-out evaluation").
  - Previously the best model was chosen by test-set R2 (test leakage),
    and cross_validate() was defined but never called while the training
    report claimed cluster-aware CV.

Models:
  - Mean predictor (baseline)
  - Linear regression
  - Ridge regression
  - Elastic net
  - Random forest
  - Gradient boosting (sklearn)
  - XGBoost
  - LightGBM
  - SVR
  - MLP (PyTorch)

Outputs:
  - model_comparison_table.csv   (cluster-aware CV results on training split)
  - best_model.pkl
  - preprocessing_pipeline.pkl
  - training_report.md
"""

import pandas as pd
import numpy as np
import os
import sys
import pickle
import json
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    DATA_DIR, MODEL_DIR, REPORT_DIR, FIGURE_DIR,
    RANDOM_SEED, HIGH_EFFICACY_THRESHOLD
)

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression, Ridge, ElasticNet
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.metrics import (
    mean_absolute_error, mean_squared_error, r2_score,
    median_absolute_error
)
from scipy.stats import pearsonr, spearmanr
import xgboost as xgb
import lightgbm as lgb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def get_feature_columns(df):
    """Get numeric feature columns (exclude identifiers and target)."""
    exclude = {"sequence", "aso_group_id", "sequence_cluster_id", "inhibition_percent"}
    feature_cols = [c for c in df.columns if c not in exclude and df[c].dtype in [np.float64, np.int64, np.float32, np.int32, float, int]]
    return feature_cols


def prepare_data(df_train, df_test):
    """Prepare feature matrices and targets."""
    feature_cols = get_feature_columns(df_train)

    X_train = df_train[feature_cols].values.astype(np.float32)
    y_train = df_train["inhibition_percent"].values.astype(np.float32)
    X_test = df_test[feature_cols].values.astype(np.float32)
    y_test = df_test["inhibition_percent"].values.astype(np.float32)

    # Handle inf/nan
    X_train = np.nan_to_num(X_train, nan=0, posinf=0, neginf=0)
    X_test = np.nan_to_num(X_test, nan=0, posinf=0, neginf=0)

    # Scale features (scaler fit on TRAINING data only; test is only transformed)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Groups for cluster-aware CV (persisted by phase5 in train.csv)
    groups = df_train["sequence_cluster_id"].values if "sequence_cluster_id" in df_train.columns else None

    return X_train, X_train_scaled, y_train, X_test, X_test_scaled, y_test, scaler, feature_cols, groups


def clip_predictions(y_pred, low=0, high=100):
    """Clip predictions to valid range."""
    return np.clip(y_pred, low, high)


def compute_metrics(y_true, y_pred):
    """Compute comprehensive regression metrics."""
    y_pred = clip_predictions(y_pred)
    metrics = {}
    metrics["MAE"] = mean_absolute_error(y_true, y_pred)
    metrics["RMSE"] = np.sqrt(mean_squared_error(y_true, y_pred))
    metrics["R2"] = r2_score(y_true, y_pred)
    metrics["MedianAE"] = median_absolute_error(y_true, y_pred)

    if len(y_true) > 2:
        metrics["Pearson_r"], metrics["Pearson_p"] = pearsonr(y_true, y_pred)
        metrics["Spearman_r"], metrics["Spearman_p"] = spearmanr(y_true, y_pred)
    else:
        metrics["Pearson_r"] = metrics["Spearman_r"] = np.nan

    # Calibration
    try:
        coeffs = np.polyfit(y_pred, y_true, 1)
        metrics["Calibration_slope"] = coeffs[0]
        metrics["Calibration_intercept"] = coeffs[1]
    except Exception:
        metrics["Calibration_slope"] = np.nan
        metrics["Calibration_intercept"] = np.nan

    # Accuracy bins
    residuals = np.abs(y_true - y_pred)
    metrics["Within_5"] = (residuals <= 5).mean() * 100
    metrics["Within_10"] = (residuals <= 10).mean() * 100
    metrics["Within_15"] = (residuals <= 15).mean() * 100

    return metrics


class MLPWrapper:
    """Minimal fit/predict wrapper around the PyTorch MLP so it can be
    used with the same CV / final-training code path as the sklearn models.

    Note: for CV model selection the MLP is trained with fewer epochs
    (cv_epochs) to keep the per-fold cost reasonable; the final model is
    trained with the full epoch budget (epochs).
    """

    def __init__(self, input_dim, epochs=100, seed=RANDOM_SEED):
        self.input_dim = input_dim
        self.epochs = epochs
        self.seed = seed

    def _build(self):
        import torch
        import torch.nn as nn

        class MLPRegressor(nn.Module):
            def __init__(self, input_dim):
                super().__init__()
                self.net = nn.Sequential(
                    nn.Linear(input_dim, 256),
                    nn.ReLU(),
                    nn.BatchNorm1d(256),
                    nn.Dropout(0.3),
                    nn.Linear(256, 128),
                    nn.ReLU(),
                    nn.BatchNorm1d(128),
                    nn.Dropout(0.2),
                    nn.Linear(128, 64),
                    nn.ReLU(),
                    nn.Linear(64, 1),
                )

            def forward(self, x):
                return self.net(x).squeeze(-1)

        torch.manual_seed(self.seed)
        return MLPRegressor(self.input_dim)

    def fit(self, X, y):
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        self.model_ = self._build()
        optimizer = torch.optim.Adam(self.model_.parameters(), lr=0.001, weight_decay=1e-4)
        criterion = nn.MSELoss()

        X_tensor = torch.FloatTensor(np.asarray(X, dtype=np.float32))
        y_tensor = torch.FloatTensor(np.asarray(y, dtype=np.float32))
        dataset = TensorDataset(X_tensor, y_tensor)
        loader = DataLoader(dataset, batch_size=64, shuffle=True)

        self.model_.train()
        for epoch in range(self.epochs):
            for batch_X, batch_y in loader:
                optimizer.zero_grad()
                pred = self.model_(batch_X)
                loss = criterion(pred, batch_y)
                loss.backward()
                optimizer.step()
        return self

    def predict(self, X):
        import torch
        self.model_.eval()
        with torch.no_grad():
            return self.model_(torch.FloatTensor(np.asarray(X, dtype=np.float32))).numpy()


def get_model_candidates(n_features, cv_epochs=30):
    """Candidate model registry. Each entry: display name, internal key,
    factory (returns a fresh model), and whether it consumes
    StandardScaler-scaled features."""
    return [
        {"name": "Mean Predictor", "key": "mean_predictor", "scaled": True,
         "factory": lambda: DummyRegressor(strategy="mean")},
        {"name": "Linear Regression", "key": "linear_regression", "scaled": True,
         "factory": lambda: LinearRegression()},
        {"name": "Ridge", "key": "ridge", "scaled": True,
         "factory": lambda: Ridge(alpha=10.0, random_state=RANDOM_SEED)},
        {"name": "Elastic Net", "key": "elastic_net", "scaled": True,
         "factory": lambda: ElasticNet(alpha=1.0, l1_ratio=0.5, max_iter=10000, random_state=RANDOM_SEED)},
        {"name": "Random Forest", "key": "random_forest", "scaled": False,
         "factory": lambda: RandomForestRegressor(
             n_estimators=500, max_depth=15, min_samples_leaf=5,
             max_features="sqrt", random_state=RANDOM_SEED, n_jobs=-1)},
        {"name": "Gradient Boosting", "key": "gradient_boosting", "scaled": False,
         "factory": lambda: GradientBoostingRegressor(
             n_estimators=500, max_depth=5, learning_rate=0.05,
             subsample=0.8, min_samples_leaf=5, random_state=RANDOM_SEED)},
        {"name": "XGBoost", "key": "xgboost", "scaled": False,
         "factory": lambda: xgb.XGBRegressor(
             n_estimators=500, max_depth=6, learning_rate=0.05,
             subsample=0.8, colsample_bytree=0.8,
             reg_alpha=0.1, reg_lambda=1.0,
             random_state=RANDOM_SEED, n_jobs=-1, verbosity=0)},
        {"name": "LightGBM", "key": "lightgbm", "scaled": False,
         "factory": lambda: lgb.LGBMRegressor(
             n_estimators=500, max_depth=6, learning_rate=0.05,
             subsample=0.8, colsample_bytree=0.8,
             reg_alpha=0.1, reg_lambda=1.0,
             random_state=RANDOM_SEED, n_jobs=-1, verbose=-1)},
        {"name": "SVR", "key": "svr", "scaled": True,
         "factory": lambda: SVR(kernel="rbf", C=10.0, epsilon=2.0)},
        {"name": "MLP (PyTorch)", "key": "mlp", "scaled": True,
         "factory": lambda: MLPWrapper(n_features, epochs=cv_epochs)},
    ]


def cross_validate(model_fn, X, y, groups, n_splits=5, scale=False):
    """Cluster-aware cross-validation on the given (training) data.

    If scale=True, a fresh StandardScaler is fit inside each fold on the
    fold-training rows only (avoids scaling leakage across folds).
    """
    if groups is None:
        from sklearn.model_selection import KFold
        cv = KFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
        splits = list(cv.split(X, y))
    else:
        unique_groups = np.unique(groups)
        n_splits = min(n_splits, len(unique_groups))
        cv = GroupKFold(n_splits=n_splits)
        splits = list(cv.split(X, y, groups))

    fold_metrics = []
    for fold_idx, (train_idx, val_idx) in enumerate(splits):
        model = model_fn()
        X_tr, y_tr = X[train_idx], y[train_idx]
        X_val, y_val = X[val_idx], y[val_idx]

        if scale:
            fold_scaler = StandardScaler()
            X_tr = fold_scaler.fit_transform(X_tr)
            X_val = fold_scaler.transform(X_val)

        model.fit(X_tr, y_tr)
        y_pred = clip_predictions(model.predict(X_val))
        metrics = compute_metrics(y_val, y_pred)
        metrics["fold"] = fold_idx
        fold_metrics.append(metrics)

    return pd.DataFrame(fold_metrics)


def select_model_by_cv(X_train, y_train, groups):
    """Select the best candidate by MEAN cluster-aware CV R2 on the
    TRAINING split only. The test set is not used anywhere here.

    Returns (best_candidate_dict, cv_summary_df, cv_folds_df).
    """
    candidates = get_model_candidates(n_features=X_train.shape[1])
    summary_rows = []
    all_folds = []

    for cand in candidates:
        print(f"  CV (cluster-aware GroupKFold, train split): {cand['name']}...")
        try:
            folds_df = cross_validate(
                cand["factory"], X_train, y_train, groups,
                n_splits=5, scale=cand["scaled"],
            )
            folds_df["Model"] = cand["name"]
            all_folds.append(folds_df)
            summary_rows.append({
                "Model": cand["name"],
                "MAE": folds_df["MAE"].mean(),
                "RMSE": folds_df["RMSE"].mean(),
                "R2": folds_df["R2"].mean(),          # mean CV R2 (selection metric)
                "R2_std": folds_df["R2"].std(),
                "MedianAE": folds_df["MedianAE"].mean(),
                "Pearson_r": folds_df["Pearson_r"].mean(),
                "Spearman_r": folds_df["Spearman_r"].mean(),
                "Within_5": folds_df["Within_5"].mean(),
                "Within_10": folds_df["Within_10"].mean(),
                "Within_15": folds_df["Within_15"].mean(),
                "n_folds": len(folds_df),
                "Notes": "",
            })
        except Exception as e:
            print(f"    CV for {cand['name']} FAILED: {e}")
            summary_rows.append({
                "Model": f"{cand['name']} [FAILED]",
                "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "R2_std": np.nan,
                "Notes": f"CV failed: {e}",
            })

    cv_summary = pd.DataFrame(summary_rows)
    cv_folds = pd.concat(all_folds, ignore_index=True) if all_folds else pd.DataFrame()

    # Selection: best MEAN CV R2 on the training split, excluding the
    # mean-predictor baseline. No test-set metric is consulted.
    eligible = [c for c, row in zip(candidates, summary_rows)
                if c["key"] != "mean_predictor" and not np.isnan(row["R2"])]
    if not eligible:
        raise RuntimeError("No candidate model completed cluster-aware CV; cannot select a model.")

    best = max(
        eligible,
        key=lambda c: next(r["R2"] for r in summary_rows if r["Model"] == c["name"]),
    )
    best_r2 = next(r["R2"] for r in summary_rows if r["Model"] == best["name"])
    print(f"[Phase 7] Selected by mean CV R^2 (train split): {best['name']} (CV R^2={best_r2:.4f})")
    return best, cv_summary, cv_folds


def train_final_model(candidate, X_train, y_train):
    """Train the CV-selected candidate once on the FULL training split."""
    print(f"  Training final model on full training split: {candidate['name']}...")
    if candidate["key"] == "mlp":
        model = MLPWrapper(X_train.shape[1], epochs=100)  # full epoch budget
    else:
        model = candidate["factory"]()
    model.fit(X_train, y_train)
    return model


def evaluate_on_test(model, candidate, X_test, y_test):
    """Final held-out evaluation: the ONLY use of the test set in phase 7."""
    y_pred = clip_predictions(model.predict(X_test))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = candidate["name"]
    return metrics


def plot_model_comparison(results_df):
    """Plot model comparison charts (mean cluster-aware CV metrics)."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Sort by R2
    results_sorted = results_df.sort_values("R2", ascending=True)

    # R2
    colors = sns.color_palette("viridis", len(results_sorted))
    axes[0].barh(results_sorted["Model"], results_sorted["R2"], color=colors)
    axes[0].set_xlabel("Mean CV R^2 (GroupKFold, train split)")
    axes[0].set_title("Model Comparison -- Mean CV R^2")
    axes[0].axvline(x=0, color="red", linestyle="--", alpha=0.5)

    # MAE
    results_sorted2 = results_df.sort_values("MAE", ascending=False)
    axes[1].barh(results_sorted2["Model"], results_sorted2["MAE"], color=colors)
    axes[1].set_xlabel("Mean CV MAE")
    axes[1].set_title("Model Comparison -- Mean CV MAE")

    # RMSE
    results_sorted3 = results_df.sort_values("RMSE", ascending=False)
    axes[2].barh(results_sorted3["Model"], results_sorted3["RMSE"], color=colors)
    axes[2].set_xlabel("Mean CV RMSE")
    axes[2].set_title("Model Comparison -- Mean CV RMSE")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "model_comparison.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[Phase 7] Model comparison plot saved")


def save_best_model(candidate, best_model, scaler, feature_cols):
    """Save the CV-selected model, scaler, and metadata."""
    best_name = candidate["name"]
    best_key = candidate["key"]
    needs_scaling = candidate["scaled"]

    # Save model
    if best_key == "mlp":
        import torch
        torch.save(best_model.model_.state_dict(), os.path.join(MODEL_DIR, "best_model_mlp.pth"))
        # Also save architecture info
        with open(os.path.join(MODEL_DIR, "mlp_config.json"), "w") as f:
            json.dump({"input_dim": len(feature_cols)}, f)
    else:
        with open(os.path.join(MODEL_DIR, "best_model.pkl"), "wb") as f:
            pickle.dump(best_model, f)

    # Save scaler + metadata (keys kept compatible with phase8/phase9/web tool)
    with open(os.path.join(MODEL_DIR, "preprocessing_pipeline.pkl"), "wb") as f:
        pickle.dump({
            "scaler": scaler,
            "feature_cols": feature_cols,
            "best_model_name": best_name,
            "best_model_key": best_key,
            "needs_scaling": needs_scaling,
        }, f)

    print(f"[Phase 7] Best model (by mean CV R^2): {best_name}")
    return best_name


def write_training_report(cv_summary, best_name, test_metrics):
    """Write training_report.md."""
    path = os.path.join(REPORT_DIR, "training_report.md")
    lines = []
    lines.append("# Phase 7 -- Model Training Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Training Configuration\n\n")
    lines.append(f"- **Random seed:** {RANDOM_SEED}\n")
    lines.append(f"- **Target variable:** inhibition_percent (bounded regression, 0-100)\n")
    lines.append(f"- **Predictions clipped to [0, 100]**\n")
    lines.append("- **Model selection:** cluster-aware GroupKFold cross-validation on the "
                 "TRAINING split only (grouped by sequence_cluster_id); best model chosen "
                 "by mean CV R^2.\n")
    lines.append("- **Test set usage:** evaluated exactly ONCE, after selection, as the "
                 "final held-out evaluation.\n\n")

    lines.append("## Model Comparison (Cluster-Aware GroupKFold CV on Training Split)\n\n")
    cols_to_show = ["Model", "MAE", "RMSE", "R2", "R2_std", "Pearson_r", "Spearman_r", "Within_5", "Within_10", "Within_15"]
    available_cols = [c for c in cols_to_show if c in cv_summary.columns]
    lines.append("| " + " | ".join(available_cols) + " |\n")
    lines.append("| " + " | ".join(["---"] * len(available_cols)) + " |\n")
    for _, row in cv_summary.iterrows():
        vals = []
        for c in available_cols:
            v = row.get(c, np.nan)
            if isinstance(v, float):
                vals.append(f"{v:.4f}" if abs(v) < 1000 else f"{v:.1f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |\n")

    if best_name:
        lines.append(f"\n## Best Model (by mean CV R^2): **{best_name}**\n\n")

    if test_metrics is not None:
        lines.append("\n## Final Held-Out Evaluation (Test Set, evaluated once)\n\n")
        lines.append("These metrics are the unbiased estimate of generalization performance; "
                     "they were NOT used for model selection.\n\n")
        lines.append("| Metric | Value |\n|---|---|\n")
        for k in ["MAE", "RMSE", "R2", "MedianAE", "Pearson_r", "Spearman_r",
                  "Within_5", "Within_10", "Within_15"]:
            if k in test_metrics:
                lines.append(f"| {k} | {test_metrics[k]:.4f} |\n")

    lines.append("\n## Models Trained\n\n")
    lines.append("1. **Mean Predictor** -- Baseline (predicts training mean)\n")
    lines.append("2. **Linear Regression** -- OLS on scaled features\n")
    lines.append("3. **Ridge** -- L2-regularized linear regression (alpha=10)\n")
    lines.append("4. **Elastic Net** -- L1+L2 regularization (alpha=1.0, l1_ratio=0.5)\n")
    lines.append("5. **Random Forest** -- 500 trees, max_depth=15\n")
    lines.append("6. **Gradient Boosting** -- 500 trees, max_depth=5, lr=0.05\n")
    lines.append("7. **XGBoost** -- 500 trees, max_depth=6, lr=0.05\n")
    lines.append("8. **LightGBM** -- 500 trees, max_depth=6, lr=0.05\n")
    lines.append("9. **SVR** -- RBF kernel, C=10, epsilon=2.0\n")
    lines.append("10. **MLP (PyTorch)** -- 3-layer (256->128->64->1), ReLU, BatchNorm, Dropout\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 7] Training report saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 7 -- MODEL TRAINING")
    print("=" * 70)

    # Load data (train.csv/test.csv written by phase5; both carry sequence_cluster_id)
    df_train = pd.read_csv(os.path.join(DATA_DIR, "train.csv"))
    df_test = pd.read_csv(os.path.join(DATA_DIR, "test.csv"))
    print(f"[Phase 7] Train: {len(df_train)}, Test: {len(df_test)}")

    # Prepare
    X_train, X_train_scaled, y_train, X_test, X_test_scaled, y_test, scaler, feature_cols, groups = prepare_data(df_train, df_test)
    print(f"[Phase 7] Features: {len(feature_cols)}")

    # -- Model selection: cluster-aware GroupKFold CV on TRAINING split only --
    # NOTE (ASO1-1 fix): the test set is NOT used for selection. Raw X_train is
    # passed; candidates that need scaling get a per-fold scaler (scale=True).
    best_candidate, cv_summary, cv_folds = select_model_by_cv(X_train, y_train, groups)

    # Save comparison table (CV results on training split)
    cv_summary.to_csv(os.path.join(DATA_DIR, "model_comparison_table.csv"), index=False)
    if len(cv_folds) > 0:
        cv_folds.to_csv(os.path.join(DATA_DIR, "model_comparison_cv_folds.csv"), index=False)

    # Plot comparison
    plot_model_comparison(cv_summary)

    # -- Train the selected model once on the full training split ------------
    X_final = X_train_scaled if best_candidate["scaled"] else X_train
    best_model = train_final_model(best_candidate, X_final, y_train)

    # -- Final held-out evaluation: the ONLY use of the test set -------------
    X_test_final = X_test_scaled if best_candidate["scaled"] else X_test
    test_metrics = evaluate_on_test(best_model, best_candidate, X_test_final, y_test)
    print(f"[Phase 7] Final held-out evaluation: R^2={test_metrics['R2']:.4f}, MAE={test_metrics['MAE']:.2f}")

    # Save best model + preprocessing pipeline
    best_name = save_best_model(best_candidate, best_model, scaler, feature_cols)

    # Write report
    write_training_report(cv_summary, best_name, test_metrics)

    print(f"\n[Phase 7] COMPLETE [OK]\n")
    return cv_summary, best_model, scaler, feature_cols


if __name__ == "__main__":
    main()
