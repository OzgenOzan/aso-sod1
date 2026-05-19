"""
Phase 7 -- Modeling
====================
Train baseline and advanced models for inhibition_percent prediction.

Models:
  - Mean predictor
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
  - model_comparison_table.csv
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

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Groups for cluster-aware CV
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
    from numpy.polynomial import polynomial as P
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


def cross_validate(model_fn, X, y, groups, n_splits=5):
    """Cluster-aware cross-validation."""
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

        model.fit(X_tr, y_tr)
        y_pred = clip_predictions(model.predict(X_val))
        metrics = compute_metrics(y_val, y_pred)
        metrics["fold"] = fold_idx
        fold_metrics.append(metrics)

    return pd.DataFrame(fold_metrics)


def train_models(X_train, X_train_scaled, y_train, X_test, X_test_scaled, y_test, groups, feature_cols):
    """Train all models and return results."""
    results = []
    trained_models = {}

    # -- 1. Mean Predictor --------------------------------------------
    print("  Training: Mean Predictor...")
    mean_pred = np.full_like(y_test, y_train.mean())
    metrics = compute_metrics(y_test, mean_pred)
    metrics["Model"] = "Mean Predictor"
    results.append(metrics)

    # -- 2. Linear Regression -----------------------------------------
    print("  Training: Linear Regression...")
    lr = LinearRegression()
    lr.fit(X_train_scaled, y_train)
    y_pred = clip_predictions(lr.predict(X_test_scaled))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "Linear Regression"
    results.append(metrics)
    trained_models["linear_regression"] = lr

    # -- 3. Ridge -----------------------------------------------------
    print("  Training: Ridge Regression...")
    ridge = Ridge(alpha=10.0, random_state=RANDOM_SEED)
    ridge.fit(X_train_scaled, y_train)
    y_pred = clip_predictions(ridge.predict(X_test_scaled))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "Ridge"
    results.append(metrics)
    trained_models["ridge"] = ridge

    # -- 4. Elastic Net -----------------------------------------------
    print("  Training: Elastic Net...")
    enet = ElasticNet(alpha=1.0, l1_ratio=0.5, max_iter=10000, random_state=RANDOM_SEED)
    enet.fit(X_train_scaled, y_train)
    y_pred = clip_predictions(enet.predict(X_test_scaled))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "Elastic Net"
    results.append(metrics)
    trained_models["elastic_net"] = enet

    # -- 5. Random Forest ---------------------------------------------
    print("  Training: Random Forest...")
    rf = RandomForestRegressor(
        n_estimators=500, max_depth=15, min_samples_leaf=5,
        max_features="sqrt", random_state=RANDOM_SEED, n_jobs=-1
    )
    rf.fit(X_train, y_train)
    y_pred = clip_predictions(rf.predict(X_test))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "Random Forest"
    results.append(metrics)
    trained_models["random_forest"] = rf

    # -- 6. Gradient Boosting -----------------------------------------
    print("  Training: Gradient Boosting...")
    gb = GradientBoostingRegressor(
        n_estimators=500, max_depth=5, learning_rate=0.05,
        subsample=0.8, min_samples_leaf=5, random_state=RANDOM_SEED
    )
    gb.fit(X_train, y_train)
    y_pred = clip_predictions(gb.predict(X_test))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "Gradient Boosting"
    results.append(metrics)
    trained_models["gradient_boosting"] = gb

    # -- 7. XGBoost ---------------------------------------------------
    print("  Training: XGBoost...")
    xgb_model = xgb.XGBRegressor(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=RANDOM_SEED, n_jobs=-1, verbosity=0
    )
    xgb_model.fit(X_train, y_train)
    y_pred = clip_predictions(xgb_model.predict(X_test))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "XGBoost"
    results.append(metrics)
    trained_models["xgboost"] = xgb_model

    # -- 8. LightGBM --------------------------------------------------
    print("  Training: LightGBM...")
    lgb_model = lgb.LGBMRegressor(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=RANDOM_SEED, n_jobs=-1, verbose=-1
    )
    lgb_model.fit(X_train, y_train)
    y_pred = clip_predictions(lgb_model.predict(X_test))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "LightGBM"
    results.append(metrics)
    trained_models["lightgbm"] = lgb_model

    # -- 9. SVR -------------------------------------------------------
    print("  Training: SVR...")
    svr = SVR(kernel="rbf", C=10.0, epsilon=2.0)
    svr.fit(X_train_scaled, y_train)
    y_pred = clip_predictions(svr.predict(X_test_scaled))
    metrics = compute_metrics(y_test, y_pred)
    metrics["Model"] = "SVR"
    results.append(metrics)
    trained_models["svr"] = svr

    # -- 10. MLP (PyTorch) --------------------------------------------
    print("  Training: MLP (PyTorch)...")
    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        torch.manual_seed(RANDOM_SEED)

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

        model_mlp = MLPRegressor(X_train_scaled.shape[1])
        optimizer = torch.optim.Adam(model_mlp.parameters(), lr=0.001, weight_decay=1e-4)
        criterion = nn.MSELoss()

        X_tensor = torch.FloatTensor(X_train_scaled)
        y_tensor = torch.FloatTensor(y_train)
        dataset = TensorDataset(X_tensor, y_tensor)
        loader = DataLoader(dataset, batch_size=64, shuffle=True)

        model_mlp.train()
        for epoch in range(100):
            for batch_X, batch_y in loader:
                optimizer.zero_grad()
                pred = model_mlp(batch_X)
                loss = criterion(pred, batch_y)
                loss.backward()
                optimizer.step()

        model_mlp.eval()
        with torch.no_grad():
            X_test_tensor = torch.FloatTensor(X_test_scaled)
            y_pred_mlp = model_mlp(X_test_tensor).numpy()
            y_pred_mlp = clip_predictions(y_pred_mlp)

        metrics = compute_metrics(y_test, y_pred_mlp)
        metrics["Model"] = "MLP (PyTorch)"
        results.append(metrics)
        trained_models["mlp"] = model_mlp

    except Exception as e:
        print(f"    MLP training failed: {e}")
        metrics = {"Model": "MLP (PyTorch)", "MAE": np.nan, "RMSE": np.nan, "R2": np.nan}
        results.append(metrics)

    return pd.DataFrame(results), trained_models


def plot_model_comparison(results_df):
    """Plot model comparison charts."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Sort by R2
    results_sorted = results_df.sort_values("R2", ascending=True)

    # R2
    colors = sns.color_palette("viridis", len(results_sorted))
    axes[0].barh(results_sorted["Model"], results_sorted["R2"], color=colors)
    axes[0].set_xlabel("R^2 Score")
    axes[0].set_title("Model Comparison -- R^2")
    axes[0].axvline(x=0, color="red", linestyle="--", alpha=0.5)

    # MAE
    results_sorted2 = results_df.sort_values("MAE", ascending=False)
    axes[1].barh(results_sorted2["Model"], results_sorted2["MAE"], color=colors)
    axes[1].set_xlabel("MAE")
    axes[1].set_title("Model Comparison -- MAE")

    # RMSE
    results_sorted3 = results_df.sort_values("RMSE", ascending=False)
    axes[2].barh(results_sorted3["Model"], results_sorted3["RMSE"], color=colors)
    axes[2].set_xlabel("RMSE")
    axes[2].set_title("Model Comparison -- RMSE")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "model_comparison.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[Phase 7] Model comparison plot saved")


def save_best_model(trained_models, results_df, scaler, feature_cols):
    """Save the best model, scaler, and metadata."""
    # Select best by R2 (excluding Mean Predictor)
    valid_results = results_df[results_df["Model"] != "Mean Predictor"].dropna(subset=["R2"])
    best_row = valid_results.loc[valid_results["R2"].idxmax()]
    best_name = best_row["Model"]

    # Map display name to dict key
    name_map = {
        "Linear Regression": "linear_regression",
        "Ridge": "ridge",
        "Elastic Net": "elastic_net",
        "Random Forest": "random_forest",
        "Gradient Boosting": "gradient_boosting",
        "XGBoost": "xgboost",
        "LightGBM": "lightgbm",
        "SVR": "svr",
        "MLP (PyTorch)": "mlp",
    }

    best_key = name_map.get(best_name, "")
    best_model = trained_models.get(best_key)

    if best_model is not None:
        # Determine if model needs scaled input
        needs_scaling = best_key in ("linear_regression", "ridge", "elastic_net", "svr", "mlp")

        # Save model
        if best_key == "mlp":
            import torch
            torch.save(best_model.state_dict(), os.path.join(MODEL_DIR, "best_model_mlp.pth"))
            # Also save architecture info
            with open(os.path.join(MODEL_DIR, "mlp_config.json"), "w") as f:
                json.dump({"input_dim": len(feature_cols)}, f)
        else:
            with open(os.path.join(MODEL_DIR, "best_model.pkl"), "wb") as f:
                pickle.dump(best_model, f)

        # Save scaler
        with open(os.path.join(MODEL_DIR, "preprocessing_pipeline.pkl"), "wb") as f:
            pickle.dump({
                "scaler": scaler,
                "feature_cols": feature_cols,
                "best_model_name": best_name,
                "best_model_key": best_key,
                "needs_scaling": needs_scaling,
            }, f)

        print(f"[Phase 7] Best model: {best_name} (R^2={best_row['R2']:.4f})")
        return best_name, best_row
    else:
        print(f"[Phase 7] WARNING: Best model '{best_name}' not found in trained_models dict")
        return None, None


def write_training_report(results_df, best_name, best_metrics):
    """Write training_report.md."""
    path = os.path.join(REPORT_DIR, "training_report.md")
    lines = []
    lines.append("# Phase 7 -- Model Training Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Training Configuration\n\n")
    lines.append(f"- **Random seed:** {RANDOM_SEED}\n")
    lines.append(f"- **Target variable:** inhibition_percent (bounded regression, 0-100)\n")
    lines.append(f"- **Predictions clipped to [0, 100]**\n")
    lines.append(f"- **Cross-validation:** Cluster-aware GroupKFold\n\n")

    lines.append("## Model Comparison (Test Set)\n\n")
    cols_to_show = ["Model", "MAE", "RMSE", "R2", "Pearson_r", "Spearman_r", "Within_5", "Within_10", "Within_15"]
    available_cols = [c for c in cols_to_show if c in results_df.columns]
    lines.append("| " + " | ".join(available_cols) + " |\n")
    lines.append("| " + " | ".join(["---"] * len(available_cols)) + " |\n")
    for _, row in results_df.iterrows():
        vals = []
        for c in available_cols:
            v = row.get(c, np.nan)
            if isinstance(v, float):
                vals.append(f"{v:.4f}" if abs(v) < 1000 else f"{v:.1f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |\n")

    if best_name:
        lines.append(f"\n## Best Model: **{best_name}**\n\n")
        if best_metrics is not None:
            lines.append(f"- R^2 = {best_metrics.get('R2', 'N/A'):.4f}\n")
            lines.append(f"- MAE = {best_metrics.get('MAE', 'N/A'):.4f}\n")
            lines.append(f"- RMSE = {best_metrics.get('RMSE', 'N/A'):.4f}\n")

    lines.append("\n## Models Trained\n\n")
    lines.append("1. **Mean Predictor** -- Baseline (predicts training mean)\n")
    lines.append("2. **Linear Regression** -- OLS on scaled features\n")
    lines.append("3. **Ridge** -- L2-regularized linear regression (α=10)\n")
    lines.append("4. **Elastic Net** -- L1+L2 regularization (α=1.0, l1_ratio=0.5)\n")
    lines.append("5. **Random Forest** -- 500 trees, max_depth=15\n")
    lines.append("6. **Gradient Boosting** -- 500 trees, max_depth=5, lr=0.05\n")
    lines.append("7. **XGBoost** -- 500 trees, max_depth=6, lr=0.05\n")
    lines.append("8. **LightGBM** -- 500 trees, max_depth=6, lr=0.05\n")
    lines.append("9. **SVR** -- RBF kernel, C=10, ε=2.0\n")
    lines.append("10. **MLP (PyTorch)** -- 3-layer (256->128->64->1), ReLU, BatchNorm, Dropout\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 7] Training report saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 7 -- MODEL TRAINING")
    print("=" * 70)

    # Load data
    df_train = pd.read_csv(os.path.join(DATA_DIR, "train.csv"))
    df_test = pd.read_csv(os.path.join(DATA_DIR, "test.csv"))
    print(f"[Phase 7] Train: {len(df_train)}, Test: {len(df_test)}")

    # Prepare
    X_train, X_train_scaled, y_train, X_test, X_test_scaled, y_test, scaler, feature_cols, groups = prepare_data(df_train, df_test)
    print(f"[Phase 7] Features: {len(feature_cols)}")

    # Train all models
    results_df, trained_models = train_models(
        X_train, X_train_scaled, y_train,
        X_test, X_test_scaled, y_test,
        groups, feature_cols
    )

    # Save comparison table
    results_df.to_csv(os.path.join(DATA_DIR, "model_comparison_table.csv"), index=False)

    # Plot comparison
    plot_model_comparison(results_df)

    # Save best model
    best_name, best_metrics = save_best_model(trained_models, results_df, scaler, feature_cols)

    # Write report
    write_training_report(results_df, best_name, best_metrics)

    print(f"\n[Phase 7] COMPLETE [OK]\n")
    return results_df, trained_models, scaler, feature_cols


if __name__ == "__main__":
    main()
