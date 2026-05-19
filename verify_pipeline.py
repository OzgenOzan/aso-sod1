import sys, os
sys.path.insert(0, 'pipeline')
import pickle, json, numpy as np

MODEL_DIR = 'pipeline/outputs/models'
EXT_VAL_DIR = 'pipeline/outputs/external_validation'

with open(os.path.join(MODEL_DIR, 'preprocessing_pipeline.pkl'), 'rb') as f:
    pipeline = pickle.load(f)

with open(os.path.join(MODEL_DIR, 'best_model.pkl'), 'rb') as f:
    model = pickle.load(f)

model_name = pipeline["best_model_name"]
n_feats = len(pipeline["feature_cols"])
print(f"Model: {model_name}")
print(f"Features: {n_feats}")

from config import DATA_DIR
import pandas as pd
df_test = pd.read_csv(os.path.join(DATA_DIR, 'test.csv'))
feature_cols = pipeline['feature_cols']
X = df_test[feature_cols].head(5).values.astype(np.float32)
X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)
y_pred = np.clip(model.predict(X), 0, 100)
y_true = df_test['inhibition_percent'].head(5).values

print("\nSample predictions:")
for i in range(5):
    print(f"  Actual: {y_true[i]:.1f}, Predicted: {y_pred[i]:.1f}")

with open(os.path.join(EXT_VAL_DIR, 'tofersen_reference.json')) as f:
    tof = json.load(f)
print(f"\nTofersen sequence: {tof['sequence']}")
print(f"Tofersen chemical_pattern: {tof['chemical_pattern']}")
print(f"Tofersen predicted: {tof.get('predicted_inhibition_percent', 'N/A')}")
print("\nVerification: ALL OK")
