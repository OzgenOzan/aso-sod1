# ASO-SOD1 — ASO Potency Prediction Pipeline for SOD1

Machine-learning pipeline for predicting the in-vitro mRNA inhibition (potency) of
antisense oligonucleotides (ASOs) targeting the human **SOD1** transcript, including a
Streamlit web tool that benchmarks new predictions against a tofersen reference.

## Installation

Requires **Python 3.12**.

```bash
pip install -r requirements.txt
```

## How to run

The pipeline is organized into sequential phases. Run them in order:

```bash
cd pipeline
python run_pipeline.py        # runs all phases end-to-end
```

Or run individual phases in this order:

| Phase | Script | Purpose |
|---|---|---|
| 1 | `phase1_audit.py` | Raw data audit of `dataset/sod-1_data.xlsx` |
| 2 | `phase2_cleaning.py` | Cleaning, deduplication, replicate aggregation |
| 3 | `phase3_clustering.py` | Sequence clustering for leakage-safe splitting |
| 4 | `phase4_features.py` | Feature engineering (sequence, chemistry, conditions) |
| 5 | `phase5_split.py` | Cluster-aware train/test split |
| 6 | `phase6_external.py` | External reference data preparation (tofersen) |
| 7 | `phase7_modeling.py` | Model training & comparison (RF/GBM/XGBoost/LightGBM/MLP) |
| 8 | `phase8_validation.py` | Validation, cluster-aware CV, SHAP, model card |
| 9 | `phase9_tofersen.py` | Tofersen benchmark prediction under standardized conditions |

Sanity-check the committed artifacts with:

```bash
python verify_pipeline.py
```

Launch the web tool:

```bash
streamlit run pipeline/outputs/web_tool/app.py
```

## Repository layout

```
dataset/                  Input data (sod-1_data.xlsx, tofersen.json)
pipeline/                 Pipeline source code
  config.py               Shared paths/constants
  model_integrity.py      SHA-256-verified pickle loading helpers
  phase1_audit.py ... phase9_tofersen.py
  run_pipeline.py         End-to-end runner
  outputs/
    data/                 Processed datasets, metrics tables
    models/               Trained model + preprocessing pipeline (pickle)
    figures/              Generated figures
    reports/              Per-phase reports
    external_validation/  Tofersen reference outputs
    web_tool/             Streamlit app (app.py)
documentation/            Final pipeline report and model card
articles/                 Reference articles (PDF)
tests/                    Static remediation regression tests (no heavy deps)
verify_pipeline.py        Smoke check of committed model artifacts
```

See `documentation/model_card.md` for the model card (intended use, training data,
performance, limitations) and `documentation/final_pipeline_report.md` for the full
phase-by-phase report.

## Notes

- Model artifacts are distributed as pickle files (`pipeline/outputs/models/*.pkl`);
  only load them from a trusted checkout of this repository.
- **Model artifact integrity:** pickle deserialization can execute arbitrary code.
  All pipeline phases, `verify_pipeline.py`, and the web tool load artifacts through
  `pipeline/model_integrity.load_pickle_verified()`, which verifies the file's
  SHA-256 hash before deserializing when the expected hash is configured via
  environment variables:

  ```bash
  export ASO_BEST_MODEL_SHA256=<sha256 of pipeline/outputs/models/best_model.pkl>
  export ASO_PREPROCESSING_SHA256=<sha256 of pipeline/outputs/models/preprocessing_pipeline.pkl>
  ```

  Compute the hashes with e.g. `sha256sum pipeline/outputs/models/*.pkl` after a
  trusted checkout or retraining. If a variable is unset, a loud warning is emitted
  and the artifact is loaded unverified (backward compatible). If a hash is set and
  does not match, loading is refused.
- The tofersen benchmark is an in-silico model self-consistency reference, **not**
  experimental validation.
- No license file is currently provided; all rights reserved by the author.
