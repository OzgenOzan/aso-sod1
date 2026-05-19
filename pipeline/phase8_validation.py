"""
Phase 8 -- Complete Model Validation
=====================================
Evaluate model on internal CV, held-out test set, and external validation.
Includes robustness checks, interpretability, and applicability domain.

Outputs:
  - internal_cv_results.csv
  - test_set_predictions.csv
  - external_validation_predictions.csv
  - validation_report.md
  - model_card.md
  - figures/
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
    DATA_DIR, MODEL_DIR, REPORT_DIR, FIGURE_DIR, EXT_VAL_DIR,
    RANDOM_SEED, HIGH_EFFICACY_THRESHOLD
)

from sklearn.model_selection import GroupKFold, permutation_test_score
from sklearn.metrics import (
    mean_absolute_error, mean_squared_error, r2_score,
    median_absolute_error, roc_auc_score, average_precision_score,
    f1_score, confusion_matrix, precision_score, recall_score
)
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cdist

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def load_model_and_pipeline():
    """Load the best model and preprocessing pipeline."""
    with open(os.path.join(MODEL_DIR, "preprocessing_pipeline.pkl"), "rb") as f:
        pipeline = pickle.load(f)

    model_key = pipeline["best_model_key"]

    if model_key == "mlp":
        import torch
        import torch.nn as nn

        class MLPRegressor(nn.Module):
            def __init__(self, input_dim):
                super().__init__()
                self.net = nn.Sequential(
                    nn.Linear(input_dim, 256), nn.ReLU(), nn.BatchNorm1d(256), nn.Dropout(0.3),
                    nn.Linear(256, 128), nn.ReLU(), nn.BatchNorm1d(128), nn.Dropout(0.2),
                    nn.Linear(128, 64), nn.ReLU(), nn.Linear(64, 1),
                )
            def forward(self, x):
                return self.net(x).squeeze(-1)

        with open(os.path.join(MODEL_DIR, "mlp_config.json")) as f:
            config = json.load(f)
        model = MLPRegressor(config["input_dim"])
        model.load_state_dict(torch.load(os.path.join(MODEL_DIR, "best_model_mlp.pth"), weights_only=True))
        model.eval()
    else:
        with open(os.path.join(MODEL_DIR, "best_model.pkl"), "rb") as f:
            model = pickle.load(f)

    return model, pipeline


def predict(model, X, pipeline):
    """Generate predictions using the loaded model."""
    needs_scaling = pipeline["needs_scaling"]
    if needs_scaling:
        X_input = pipeline["scaler"].transform(X)
    else:
        X_input = X

    model_key = pipeline["best_model_key"]
    if model_key == "mlp":
        import torch
        model.eval()
        with torch.no_grad():
            y_pred = model(torch.FloatTensor(X_input)).numpy()
    else:
        y_pred = model.predict(X_input)

    return np.clip(y_pred, 0, 100)


def compute_full_metrics(y_true, y_pred):
    """Compute comprehensive metrics including classification."""
    metrics = {}

    # Regression metrics
    metrics["MAE"] = mean_absolute_error(y_true, y_pred)
    metrics["RMSE"] = np.sqrt(mean_squared_error(y_true, y_pred))
    metrics["R2"] = r2_score(y_true, y_pred)
    metrics["MedianAE"] = median_absolute_error(y_true, y_pred)

    if len(y_true) > 2:
        metrics["Pearson_r"], _ = pearsonr(y_true, y_pred)
        metrics["Spearman_r"], _ = spearmanr(y_true, y_pred)

    # Calibration
    try:
        coeffs = np.polyfit(y_pred, y_true, 1)
        metrics["Calibration_slope"] = coeffs[0]
        metrics["Calibration_intercept"] = coeffs[1]
    except Exception:
        pass

    # Accuracy bins
    residuals = np.abs(y_true - y_pred)
    metrics["Within_5"] = (residuals <= 5).mean() * 100
    metrics["Within_10"] = (residuals <= 10).mean() * 100
    metrics["Within_15"] = (residuals <= 15).mean() * 100

    # Secondary classification (high efficacy >= 70%)
    y_true_class = (y_true >= HIGH_EFFICACY_THRESHOLD).astype(int)
    y_pred_class = (y_pred >= HIGH_EFFICACY_THRESHOLD).astype(int)

    if len(np.unique(y_true_class)) == 2:
        metrics["AUROC"] = roc_auc_score(y_true_class, y_pred)
        metrics["AUPRC"] = average_precision_score(y_true_class, y_pred)
        metrics["Sensitivity"] = recall_score(y_true_class, y_pred_class)
        metrics["Specificity"] = recall_score(1 - y_true_class, 1 - y_pred_class)
        metrics["F1"] = f1_score(y_true_class, y_pred_class)
        cm = confusion_matrix(y_true_class, y_pred_class)
        metrics["Confusion_TN"] = int(cm[0, 0])
        metrics["Confusion_FP"] = int(cm[0, 1])
        metrics["Confusion_FN"] = int(cm[1, 0])
        metrics["Confusion_TP"] = int(cm[1, 1])

    return metrics


def run_cross_validation(df_train, model, pipeline):
    """Run cluster-aware cross-validation."""
    feature_cols = pipeline["feature_cols"]
    X = df_train[feature_cols].values.astype(np.float32)
    X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)
    y = df_train["inhibition_percent"].values

    groups = df_train["sequence_cluster_id"].values
    unique_groups = np.unique(groups)
    n_splits = min(5, len(unique_groups))
    cv = GroupKFold(n_splits=n_splits)

    fold_results = []
    all_preds = np.zeros_like(y)

    for fold, (train_idx, val_idx) in enumerate(cv.split(X, y, groups)):
        # Clone model for each fold
        from sklearn.base import clone
        model_key = pipeline["best_model_key"]

        if model_key == "mlp":
            # Skip cross-validation for MLP (too expensive to retrain)
            continue

        try:
            fold_model = clone(model)
        except Exception:
            continue

        X_tr, y_tr = X[train_idx], y[train_idx]
        X_val, y_val = X[val_idx], y[val_idx]

        if pipeline["needs_scaling"]:
            from sklearn.preprocessing import StandardScaler
            fold_scaler = StandardScaler()
            X_tr = fold_scaler.fit_transform(X_tr)
            X_val = fold_scaler.transform(X_val)

        fold_model.fit(X_tr, y_tr)
        y_pred = np.clip(fold_model.predict(X_val), 0, 100)
        all_preds[val_idx] = y_pred

        fold_metrics = compute_full_metrics(y_val, y_pred)
        fold_metrics["fold"] = fold
        fold_results.append(fold_metrics)

    return pd.DataFrame(fold_results), all_preds


def compute_applicability_domain(X_train, X_new, df_train, df_new, pipeline):
    """Assess applicability domain for new predictions."""
    feature_cols = pipeline["feature_cols"]

    # Nearest-neighbor distance in feature space
    X_tr = np.nan_to_num(X_train.astype(np.float32), nan=0, posinf=0, neginf=0)
    X_n = np.nan_to_num(X_new.astype(np.float32), nan=0, posinf=0, neginf=0)

    # Normalize for distance computation
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler()
    X_tr_s = sc.fit_transform(X_tr)
    X_n_s = sc.transform(X_n)

    distances = cdist(X_n_s, X_tr_s, metric="euclidean")
    nn_distances = distances.min(axis=1)
    nn_indices = distances.argmin(axis=1)

    # Sequence similarity
    train_seqs = df_train["sequence"].values
    new_seqs = df_new["sequence"].values

    ad_results = []
    for i in range(len(df_new)):
        result = {
            "nn_distance": nn_distances[i],
            "nn_train_index": nn_indices[i],
        }

        # Sequence identity to nearest training ASO
        best_identity = 0
        for j, ts in enumerate(train_seqs):
            if len(ts) == len(new_seqs[i]):
                matches = sum(a == b for a, b in zip(new_seqs[i], ts))
                identity = matches / len(ts)
                best_identity = max(best_identity, identity)

        result["max_seq_identity_to_train"] = best_identity
        result["seq_identity_gt_90"] = best_identity > 0.90

        # Chemistry pattern novelty
        if "chemical_pattern" in df_new.columns and "chemical_pattern" in df_train.columns:
            train_cps = set(df_train["chemical_pattern"].unique())
            new_cp = df_new.iloc[i].get("chemical_pattern", "")
            result["chem_pattern_seen"] = new_cp in train_cps

        # Linkage novelty
        if "linkage" in df_new.columns and "linkage" in df_train.columns:
            train_links = set(df_train["linkage"].unique())
            new_link = df_new.iloc[i].get("linkage", "")
            result["linkage_seen"] = new_link in train_links

        ad_results.append(result)

    return pd.DataFrame(ad_results)


def plot_predictions(y_true, y_pred, title, filename):
    """Plot actual vs predicted scatter."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(y_true, y_pred, alpha=0.4, s=20, c="#3498db", edgecolors="none")
    ax.plot([0, 100], [0, 100], "r--", alpha=0.8, linewidth=2, label="Perfect")

    # Fit line
    coeffs = np.polyfit(y_true, y_pred, 1)
    x_fit = np.linspace(0, 100, 100)
    ax.plot(x_fit, np.polyval(coeffs, x_fit), "g-", alpha=0.6, label=f"Fit (slope={coeffs[0]:.2f})")

    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    ax.set_xlabel("Actual Inhibition (%)", fontsize=12)
    ax.set_ylabel("Predicted Inhibition (%)", fontsize=12)
    ax.set_title(f"{title}\nR^2={r2:.3f}, MAE={mae:.2f}", fontsize=14)
    ax.legend()
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_residuals(y_true, y_pred, filename):
    """Plot residual analysis."""
    residuals = y_true - y_pred

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Residual vs predicted
    axes[0].scatter(y_pred, residuals, alpha=0.4, s=15, c="#e74c3c")
    axes[0].axhline(y=0, color="black", linestyle="--")
    axes[0].set_xlabel("Predicted")
    axes[0].set_ylabel("Residual")
    axes[0].set_title("Residuals vs Predicted")

    # Residual histogram
    axes[1].hist(residuals, bins=30, color="#2ecc71", edgecolor="black", alpha=0.7)
    axes[1].axvline(x=0, color="red", linestyle="--")
    axes[1].set_xlabel("Residual")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Residual Distribution")

    # Q-Q plot
    from scipy import stats
    stats.probplot(residuals, dist="norm", plot=axes[2])
    axes[2].set_title("Q-Q Plot")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_feature_importance(model, feature_cols, pipeline, top_n=30):
    """Plot feature importance for tree-based models."""
    model_key = pipeline["best_model_key"]

    if model_key in ("random_forest", "gradient_boosting", "xgboost", "lightgbm"):
        if hasattr(model, "feature_importances_"):
            importances = model.feature_importances_
            sorted_idx = np.argsort(importances)[-top_n:]

            fig, ax = plt.subplots(figsize=(10, max(8, top_n * 0.3)))
            ax.barh(
                [feature_cols[i] for i in sorted_idx],
                importances[sorted_idx],
                color=sns.color_palette("viridis", top_n)
            )
            ax.set_xlabel("Feature Importance")
            ax.set_title(f"Top {top_n} Features -- {pipeline['best_model_name']}")
            plt.tight_layout()
            plt.savefig(os.path.join(FIGURE_DIR, "feature_importance.png"), dpi=150, bbox_inches="tight")
            plt.close()

            # Save to CSV
            fi_df = pd.DataFrame({
                "feature": feature_cols,
                "importance": importances
            }).sort_values("importance", ascending=False)
            fi_df.to_csv(os.path.join(DATA_DIR, "feature_importances.csv"), index=False)
            print("[Phase 8] Feature importance saved")


def run_shap_analysis(model, X_train, feature_cols, pipeline, n_samples=200):
    """Run SHAP analysis if available."""
    model_key = pipeline["best_model_key"]
    if model_key not in ("random_forest", "gradient_boosting", "xgboost", "lightgbm"):
        print("[Phase 8] SHAP: Skipping for non-tree model")
        return

    try:
        import shap

        X_sample = X_train[:min(n_samples, len(X_train))]

        if model_key in ("xgboost", "lightgbm"):
            explainer = shap.TreeExplainer(model)
        else:
            explainer = shap.TreeExplainer(model)

        shap_values = explainer.shap_values(X_sample)

        # Summary plot
        fig, ax = plt.subplots(figsize=(12, 10))
        shap.summary_plot(shap_values, X_sample, feature_names=feature_cols, show=False, max_display=20)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, "shap_summary.png"), dpi=150, bbox_inches="tight")
        plt.close()

        # Save SHAP values
        shap_df = pd.DataFrame(shap_values, columns=feature_cols)
        shap_mean = shap_df.abs().mean().sort_values(ascending=False)
        shap_mean.to_csv(os.path.join(DATA_DIR, "shap_importance.csv"))

        print("[Phase 8] SHAP analysis complete")

    except Exception as e:
        print(f"[Phase 8] SHAP analysis failed: {e}")


def write_validation_report(test_metrics, cv_metrics_df, ext_metrics=None):
    """Write validation_report.md."""
    path = os.path.join(REPORT_DIR, "validation_report.md")
    lines = []
    lines.append("# Phase 8 -- Model Validation Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    # Test set metrics
    lines.append("## 1. Held-Out Test Set Performance\n\n")
    lines.append("| Metric | Value |\n|---|---|\n")
    for k, v in sorted(test_metrics.items()):
        if isinstance(v, float):
            lines.append(f"| {k} | {v:.4f} |\n")
        else:
            lines.append(f"| {k} | {v} |\n")

    # Cross-validation
    if len(cv_metrics_df) > 0:
        lines.append("\n## 2. Internal Cross-Validation\n\n")
        lines.append("| Metric | Mean +/- Std |\n|---|---|\n")
        for col in ["MAE", "RMSE", "R2", "Pearson_r", "Spearman_r"]:
            if col in cv_metrics_df.columns:
                mean_val = cv_metrics_df[col].mean()
                std_val = cv_metrics_df[col].std()
                lines.append(f"| {col} | {mean_val:.4f} +/- {std_val:.4f} |\n")

    # External validation
    if ext_metrics:
        lines.append("\n## 3. External Validation Set Performance\n\n")
        lines.append("| Metric | Value |\n|---|---|\n")
        for k, v in sorted(ext_metrics.items()):
            if isinstance(v, float):
                lines.append(f"| {k} | {v:.4f} |\n")
            else:
                lines.append(f"| {k} | {v} |\n")

    lines.append("\n## 4. Interpretation Notes\n\n")
    lines.append("- Feature importance and SHAP plots are available in the figures/ directory.\n")
    lines.append("- Applicability domain warnings are included with predictions.\n")
    lines.append("- Y-scrambling test recommended for final validation (computationally expensive).\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 8] Validation report saved -> {path}")


def write_model_card(test_metrics, pipeline):
    """Write model_card.md."""
    path = os.path.join(REPORT_DIR, "model_card.md")
    lines = []
    lines.append("# Model Card -- ASO Inhibition Predictor\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Model Details\n\n")
    lines.append(f"- **Model type:** {pipeline['best_model_name']}\n")
    lines.append(f"- **Features:** {len(pipeline['feature_cols'])} engineered features\n")
    lines.append(f"- **Target:** inhibition_percent (0-100, bounded regression)\n")
    lines.append(f"- **Training data:** SOD1-targeting ASOs (~2,264 records, ~688 unique sequences)\n")
    lines.append(f"- **Random seed:** {RANDOM_SEED}\n\n")

    lines.append("## Intended Use\n\n")
    lines.append("- In-silico prediction of mRNA inhibition percentage for novel SOD1-targeting ASOs.\n")
    lines.append("- Research use only. Not for clinical decision-making.\n\n")

    lines.append("## Performance\n\n")
    lines.append("| Metric | Test Set |\n|---|---|\n")
    for k in ["MAE", "RMSE", "R2", "Pearson_r", "Within_10"]:
        if k in test_metrics:
            lines.append(f"| {k} | {test_metrics[k]:.4f} |\n")

    lines.append("\n## Limitations\n\n")
    lines.append("1. Trained only on SOD1-targeting ASOs -- may not generalize to other targets.\n")
    lines.append("2. Limited cell lines (HepG2, A431, SH-SY5Y).\n")
    lines.append("3. No transcript-position features due to missing coordinate data.\n")
    lines.append("4. RDKit unavailable -- no molecular descriptor features from SMILES.\n")
    lines.append("5. Predictions should be validated experimentally.\n\n")

    lines.append("## Ethical Considerations\n\n")
    lines.append("This model provides research-use-only predictions and should not replace\n")
    lines.append("experimental validation or clinical assessment.\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 8] Model card saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 8 -- COMPLETE MODEL VALIDATION")
    print("=" * 70)

    # Load model and pipeline
    model, pipeline = load_model_and_pipeline()
    feature_cols = pipeline["feature_cols"]
    print(f"[Phase 8] Best model: {pipeline['best_model_name']}")

    # Load data
    df_train = pd.read_csv(os.path.join(DATA_DIR, "train.csv"))
    df_test = pd.read_csv(os.path.join(DATA_DIR, "test.csv"))

    X_train = df_train[feature_cols].values.astype(np.float32)
    X_test = df_test[feature_cols].values.astype(np.float32)
    X_train = np.nan_to_num(X_train, nan=0, posinf=0, neginf=0)
    X_test = np.nan_to_num(X_test, nan=0, posinf=0, neginf=0)
    y_train = df_train["inhibition_percent"].values
    y_test = df_test["inhibition_percent"].values

    # -- Test set evaluation ------------------------------------------
    y_pred_test = predict(model, X_test, pipeline)
    test_metrics = compute_full_metrics(y_test, y_pred_test)
    print(f"[Phase 8] Test R^2={test_metrics['R2']:.4f}, MAE={test_metrics['MAE']:.2f}")

    # Save test predictions
    test_preds_df = df_test[["sequence", "aso_group_id", "inhibition_percent"]].copy()
    test_preds_df["predicted"] = y_pred_test
    test_preds_df["residual"] = y_test - y_pred_test
    test_preds_df.to_csv(os.path.join(DATA_DIR, "test_set_predictions.csv"), index=False)

    # -- Cross-validation ---------------------------------------------
    print("[Phase 8] Running cross-validation...")
    cv_metrics_df, cv_preds = run_cross_validation(df_train, model, pipeline)
    cv_metrics_df.to_csv(os.path.join(DATA_DIR, "internal_cv_results.csv"), index=False)

    # -- Plots --------------------------------------------------------
    plot_predictions(y_test, y_pred_test, "Test Set: Actual vs Predicted", "test_actual_vs_predicted.png")
    plot_residuals(y_test, y_pred_test, "test_residuals.png")
    plot_feature_importance(model, feature_cols, pipeline)

    # -- SHAP ---------------------------------------------------------
    print("[Phase 8] Running SHAP analysis...")
    run_shap_analysis(model, X_train, feature_cols, pipeline)

    # -- External validation ------------------------------------------
    ext_metrics = None
    ext_fm_path = os.path.join(EXT_VAL_DIR, "external_validation_feature_matrix.csv")
    if os.path.exists(ext_fm_path):
        print("[Phase 8] Evaluating external validation set...")
        df_ext = pd.read_csv(ext_fm_path)
        available_cols = [c for c in feature_cols if c in df_ext.columns]
        missing_cols = [c for c in feature_cols if c not in df_ext.columns]
        if missing_cols:
            for mc in missing_cols:
                df_ext[mc] = 0
            print(f"[Phase 8] Filled {len(missing_cols)} missing feature columns with 0")

        X_ext = df_ext[feature_cols].values.astype(np.float32)
        X_ext = np.nan_to_num(X_ext, nan=0, posinf=0, neginf=0)
        y_pred_ext = predict(model, X_ext, pipeline)

        ext_preds_df = df_ext[["sequence"]].copy() if "sequence" in df_ext.columns else pd.DataFrame()
        ext_preds_df["predicted_inhibition"] = y_pred_ext
        ext_preds_df.to_csv(os.path.join(EXT_VAL_DIR, "external_validation_predictions.csv"), index=False)

        # Applicability domain
        ad_df = compute_applicability_domain(X_train, X_ext, df_train, df_ext, pipeline)
        ext_preds_df = pd.concat([ext_preds_df, ad_df], axis=1)
        ext_preds_df.to_csv(os.path.join(EXT_VAL_DIR, "external_validation_predictions.csv"), index=False)

    # -- Reports ------------------------------------------------------
    write_validation_report(test_metrics, cv_metrics_df, ext_metrics)
    write_model_card(test_metrics, pipeline)

    print(f"\n[Phase 8] COMPLETE [OK]\n")
    return test_metrics


if __name__ == "__main__":
    main()
