import os
import sys
import pickle
import numpy as np
import pandas as pd

base_dir = r"c:\Users\oozgen\Desktop\Udemy\Gemini\Codebase\20260513_ASO_4"
models_dir = os.path.join(base_dir, "pipeline", "outputs", "models")

# Load preprocessing pipeline
with open(os.path.join(models_dir, "preprocessing_pipeline.pkl"), "rb") as f:
    prep = pickle.load(f)

feature_cols = prep["feature_cols"]
best_model_name = prep["best_model_name"]
best_model_key = prep["best_model_key"]

print(f"Best model: {best_model_name}")

# Load best model
if best_model_key == "mlp":
    print("Best model is MLP, extracting importances using permutation importance or weights is complex.")
    # For now, let's just fall back to Random Forest or XGBoost if available, 
    # but the instructions say "fetch the most important 100 features", so maybe the model is tree-based.
else:
    with open(os.path.join(models_dir, "best_model.pkl"), "rb") as f:
        model = pickle.load(f)

importances = None
if hasattr(model, "feature_importances_"):
    importances = model.feature_importances_
elif hasattr(model, "coef_"):
    importances = np.abs(model.coef_)
    if importances.ndim > 1:
        importances = importances.flatten()

if importances is not None:
    df_imp = pd.DataFrame({
        "feature": feature_cols,
        "importance": importances
    })
    df_imp = df_imp.sort_values(by="importance", ascending=False).head(100)
    df_imp.to_csv(os.path.join(base_dir, "top_100_features.csv"), index=False)
    print("Saved top 100 features to top_100_features.csv")
    for idx, row in df_imp.head(10).iterrows():
        print(f"{row['feature']}: {row['importance']:.4f}")
else:
    print("Model does not have feature_importances_ or coef_")
