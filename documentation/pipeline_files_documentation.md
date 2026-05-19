# Pipeline Files Documentation

> Complete reference for every file in the `pipeline/` directory.

---

## Overview

The `pipeline/` folder contains 11 Python files and an `outputs/` directory. The scripts execute sequentially (Phase 1→9), each reading the previous phase's output.

```
pipeline/
├── config.py              # Configuration & constants
├── run_pipeline.py        # Master runner
├── phase1_audit.py        # Data audit
├── phase2_cleaning.py     # Data cleaning
├── phase3_clustering.py   # Sequence clustering
├── phase4_features.py     # Feature extraction
├── phase5_split.py        # Train/test split
├── phase6_external.py     # External validation
├── phase7_modeling.py      # Model training
├── phase8_validation.py   # Model evaluation
├── phase9_tofersen.py     # Tofersen benchmark
└── outputs/               # All generated files
```

---

## 1. `config.py` — Central Configuration

**Purpose:** Single source of truth for all paths, constants, and settings. Every other script imports from here.

### Constants

| Variable | Value | Purpose |
|---|---|---|
| `RANDOM_SEED` | `42` | Reproducibility across all random operations |
| `TEST_SIZE` | `0.30` | 30% of clusters reserved for testing |
| `CLUSTER_IDENTITY_THRESHOLD` | `0.90` | Sequences >90% similar are clustered together |
| `HIGH_EFFICACY_THRESHOLD` | `70.0` | ASOs above 70% inhibition = "high efficacy" |
| `TOFERSEN_COMPARABLE_MARGIN` | `10.0` | ±10% margin for tofersen comparison |
| `INHIBITION_BINS` | `[0, 20, 40, 60, 80, 100]` | Bins for stratified splitting |

### Paths

- `BASE_DIR` — Project root (parent of `pipeline/`)
- `DATASET_DIR` → `dataset/` (raw Excel + tofersen JSON)
- `OUTPUT_DIR` → `pipeline/outputs/`
- Sub-directories: `DATA_DIR`, `REPORT_DIR`, `FIGURE_DIR`, `MODEL_DIR`, `EXT_VAL_DIR`, `WEB_TOOL_DIR`

### Column Mapping (`COLUMN_MAP`)

Maps raw Excel headers to snake_case Python names:

```python
"ISIS" → "no"
"ASO_volume(nM)" → "aso_concentration_nm"
"Inhibition(%)" → "inhibition_percent"
"Chemical_Pattern" → "chemical_pattern"
# ... 17 columns total
```

### Chemistry Codes

```python
CHEMISTRY_CODES = {
    "M": "2'-MOE (2'-O-methoxyethyl) modified nucleotide",
    "C": "cEt (constrained ethyl) modified nucleotide",
    "d": "2'-deoxy (unmodified DNA) nucleotide",
}
```

### Auto-creates output directories on import via `os.makedirs(d, exist_ok=True)`.

---

## 2. `run_pipeline.py` — Master Pipeline Runner

**Purpose:** Executes all 9 phases sequentially with timing.

**How it works:**
1. Imports each phase's `main()` function
2. Calls them in order: `phase1()` → `phase2()` → ... → `phase9()`
3. Prints total elapsed time

**Usage:** `python run_pipeline.py` (from inside `pipeline/` directory)

**Key detail:** Uses `sys.path.insert(0, ...)` to ensure imports resolve correctly regardless of working directory.

---

## 3. `phase1_audit.py` — Data Loading & Initial Audit

**Purpose:** Load raw Excel data and perform read-only quality checks.

### Functions

| Function | What it does |
|---|---|
| `load_raw_data()` | Reads `sod-1_data.xlsx` into a pandas DataFrame |
| `audit_dataset(df)` | Runs all quality checks, returns dict of findings |
| `write_audit_report(audit)` | Writes `raw_data_audit_report.md` |
| `write_column_dictionary()` | Writes `column_dictionary.md` with column descriptions |
| `main()` | Orchestrates above functions, prints summary |

### Audit Checks Performed

1. **Basic counts** — rows, columns, unique sequences, unique chemical patterns
2. **Missingness** — null count per column
3. **Inhibition statistics** — mean, std, min, max, quartiles
4. **Sequence length distribution** — min/max/distribution table
5. **Target gene validation** — confirms all records are SOD1
6. **Inhibition range** — checks all values are within [0, 100]
7. **Sequence length verification** — compares `seq_length` column vs `len(Sequence)`; finds 12 mismatches
8. **Duplicate analysis** — exact duplicates (20), duplicated sequences, duplicated seq+chemistry combos
9. **Conflicting inhibition** — groups with different inhibition under identical conditions
10. **Chemical pattern catalog** — lists all patterns with gapmer architecture (wing-gap-wing)

### Outputs
- `reports/raw_data_audit_report.md`
- `reports/column_dictionary.md`

### Returns: `(df, audit)` — raw DataFrame and audit dictionary

---

## 4. `phase2_cleaning.py` — Data Cleaning & Validation

**Purpose:** Normalize, validate, deduplicate, and aggregate the dataset.

### Functions

| Function | What it does |
|---|---|
| `load_and_rename()` | Loads Excel, renames columns via `COLUMN_MAP` |
| `normalize_sequences(df)` | Uppercases sequences, converts U→T, validates ACGT alphabet |
| `validate_seq_length(df)` | Compares stated vs computed length, corrects 12 mismatches |
| `normalize_categoricals(df)` | Standardizes cell_line, transfection, primer, modification, linkage |
| `validate_numerics(df)` | Coerces to numeric, checks ranges (density, concentration, etc.) |
| `validate_chemical_pattern(df)` | Checks pattern length matches sequence, valid codes only (M/C/d) |
| `parse_location(df)` | Parses `Location` column into numeric position indices |
| `parse_linkage_location(df)` | Parses `Linkage_Location` into structured linkage positions |
| `validate_smiles(df)` | Checks non-empty (RDKit unavailable for chemical validation) |
| `handle_duplicates(df)` | Removes 20 exact duplicates, aggregates biological replicates by mean |
| `write_cleaning_report(...)` | Writes `cleaning_report.md` |
| `main()` | Orchestrates all steps, saves outputs |

### Key Decisions

- **Transfection merge:** `"uptake"` and `"free uptake"` → `"free_uptake"` (considered equivalent)
- **Sequence length:** Trusts computed length over stated length (12 corrections)
- **Replicate aggregation:** Groups by 8 condition columns, takes mean inhibition
- **Unhashable columns** (`parsed_mod_positions`, `parsed_linkage_positions`) are dropped before deduplication, then excluded from CSV output

### Outputs
- `data/cleaned_dataset.csv` — 2,155 rows
- `data/excluded_or_flagged_records.csv`
- `reports/cleaning_report.md`

### Returns: `df_clean` (cleaned DataFrame)

---

## 5. `phase3_clustering.py` — Sequence Redundancy Clustering

**Purpose:** Group similar sequences to prevent data leakage in train/test splitting.

### Functions

| Function | What it does |
|---|---|
| `load_cleaned_data()` | Reads `cleaned_dataset.csv` |
| `export_fasta(df)` | Writes unique sequences to `aso_sequences.fasta` |
| `compute_pairwise_identity(s1, s2)` | Calculates nucleotide identity between two sequences |
| `cluster_sequences(unique_seqs)` | Groups sequences >90% identity using BFS connected components |
| `assign_clusters_to_dataset(df, unique_seqs)` | Maps cluster IDs back to full dataset |
| `write_cluster_outputs(...)` | Writes TSV and report |
| `main()` | Orchestrates above |

### Algorithm Detail

1. **Group sequences by length** for efficient comparison
2. **Same-length pairs:** identity = exact matches / length
3. **Cross-length pairs** (difference ≤3 nt): sliding window best-match identity
4. **Build adjacency graph:** edges connect pairs with identity >90%
5. **BFS connected components:** each component = one cluster

### Results
- 688 unique sequences → 396 clusters
- 230,828 pairwise comparisons, 415 high-identity pairs

### Outputs
- `data/cleaned_dataset_clustered.csv` — adds `sequence_cluster_id` column
- `data/aso_sequences.fasta`
- `data/sequence_clusters.tsv`
- `reports/redundancy_report.md`

---

## 6. `phase4_features.py` — Feature Extraction

**Purpose:** Transform raw ASO data into 335 numerical features for machine learning.

### Feature Extraction Functions

| Function | Group | # Features | Description |
|---|---|---|---|
| `extract_sequence_features(df)` | A | 282 | Base counts/fractions, di/trinucleotide counts, homopolymers, terminal motifs |
| `extract_thermo_features(df)` | B | 6 | Tm estimates, self-complementarity, palindrome score, GC clamp |
| `extract_chemistry_features(df)` | C | 23 | Chemistry code counts, gapmer architecture, wing lengths, transitions |
| `extract_modification_features(df)` | D | 6 | Binary flags: contains_MOE, contains_cEt, contains_5mc, etc. |
| `extract_location_features(df)` | E | 5 | Number of location groups, position indices, min/max positions |
| `extract_linkage_features(df)` | F | 10 | PS/PO counts and fractions, linkage transitions |
| `extract_smiles_features(df)` | G | 3 | SMILES string length, stereo count, ring count (limited) |
| `extract_experimental_features(df)` | H | 10 | Log-concentration, log-density, one-hot cell_line/transfection/primer |

### Assembly Function: `build_feature_matrix(df, skip_constant_removal=False)`

1. Calls all 8 extraction functions
2. Concatenates into single DataFrame
3. Adds identifier columns (`sequence`, `aso_group_id`, `sequence_cluster_id`)
4. Adds target (`inhibition_percent`)
5. Removes constant columns (zero variance) — unless `skip_constant_removal=True` (used for single-sample prediction)
6. Fills NaN values with 0

**The `skip_constant_removal` parameter** is critical: when predicting a single ASO, many one-hot features will be 0, but they shouldn't be removed because the model expects them.

### Outputs
- `data/feature_matrix.csv` — 2,155 × 339
- `reports/feature_dictionary.md`
- `reports/feature_extraction_report.md`

---

## 7. `phase5_split.py` — Train/Test Splitting

**Purpose:** Split data 70/30 by cluster with stratification by inhibition.

### Core Function: `cluster_aware_stratified_split(df)`

1. Compute median inhibition per cluster
2. Bin clusters into 5 ranges (0–20, 20–40, ..., 80–100)
3. Within each bin, shuffle clusters (using `RANDOM_SEED`)
4. Assign first 30% of clusters to test, rest to train
5. **Assert zero cluster overlap** — prevents data leakage

### Verification: `write_split_report()`

Reports include:
- Dataset sizes (train: 1,533 rows / 277 clusters, test: 622 / 119)
- Inhibition distribution balance (mean, std, median per set)
- Inhibition bin distribution per set
- Sequence length distribution per set
- Leakage control: confirms 0 cluster overlap and checks sequence overlap

### Outputs
- `data/train.csv`, `data/test.csv`
- `reports/split_report.md`

---

## 8. `phase6_external.py` — External Validation Curation

**Purpose:** Parse tofersen from FDA GSRS JSON for external benchmarking.

### Key Functions

| Function | What it does |
|---|---|
| `parse_tofersen_json()` | Extracts sequence, chemical pattern, and linkage from FDA JSON |
| `curate_external_asos()` | Builds external validation set (currently only tofersen) |
| `transform_external_to_feature_schema()` | Fills missing experimental conditions with training medians/modes |
| `write_bibliography()` | Writes `bibliography.bib` with academic citations |
| `write_curation_report()` | Documents curation decisions and limitations |

### Tofersen Parsing Detail

- **Sequence:** Extracted from `nucleicAcid.subunits[0].sequence`, converted U→T
- **Chemical pattern:** Built from `nucleicAcid.sugars` — MOE positions → `M`, dR positions → `d`. Also checks `structuralModifications` for additional MOE positions
- **Linkage:** Identified from `nucleicAcid.linkages` — phosphorothioate + phosphodiester mix
- **Result:** 5-10-5 MOE gapmer (`MMMMMddddddddddMMMMM`)

### Standard Conditions (imputed from training data)
HepG2, electroporation, 3000 nM, 16h, 20000 cells/well, RTS3898

### Outputs
- `external_validation/external_validation_raw.csv`
- `external_validation/external_validation_clean.csv`
- `external_validation/external_validation_feature_matrix.csv`
- `external_validation/tofersen_reference.json`
- `external_validation/bibliography.bib`
- `external_validation/external_validation_curation_report.md`

---

## 9. `phase7_modeling.py` — Model Training

**Purpose:** Train 9 regression models and select the best one.

### Key Functions

| Function | What it does |
|---|---|
| `get_feature_columns(df)` | Returns numeric columns excluding identifiers and target |
| `prepare_data(df_train, df_test)` | Extracts X/y arrays, fits StandardScaler, gets cluster groups |
| `clip_predictions(y_pred)` | Clips to [0, 100] range |
| `compute_metrics(y_true, y_pred)` | Calculates MAE, RMSE, R², Pearson r, Spearman r, accuracy bins |
| `train_models(...)` | Trains all 9 models, returns results DataFrame + model dict |
| `plot_model_comparison(results_df)` | Creates 3-panel bar chart (R², MAE, RMSE) |
| `save_best_model(...)` | Saves best model as `.pkl` + preprocessing pipeline |
| `write_training_report(...)` | Documents all results |

### Model Details

- **Models needing scaled input:** Linear Regression, Ridge, Elastic Net, SVR, MLP
- **Models using raw input:** Random Forest, Gradient Boosting, XGBoost, LightGBM
- **MLP:** 3-layer PyTorch network (256→128→64→1), ReLU, BatchNorm, Dropout. Wrapped in try/except (fails gracefully if PyTorch unavailable)

### Best Model Selection

Selects by highest R² (excluding Mean Predictor). Saves:
1. `models/best_model.pkl` — the trained LightGBM model
2. `models/preprocessing_pipeline.pkl` — dict containing scaler, feature column list, model metadata, and `needs_scaling` flag

### Outputs
- `data/model_comparison_table.csv`
- `models/best_model.pkl`
- `models/preprocessing_pipeline.pkl`
- `figures/model_comparison.png`
- `reports/training_report.md`

---

## 10. `phase8_validation.py` — Model Evaluation

**Purpose:** Comprehensive evaluation including CV, plots, SHAP, and applicability domain.

### Key Functions

| Function | What it does |
|---|---|
| `load_model_and_pipeline()` | Loads best model + preprocessing pipeline from pkl files |
| `predict(model, X, pipeline)` | Generates predictions (handles scaling and MLP vs sklearn) |
| `compute_full_metrics(y_true, y_pred)` | Extended metrics including AUROC, sensitivity, specificity for high-efficacy classification |
| `run_cross_validation(...)` | 5-fold cluster-aware GroupKFold CV with per-fold model cloning |
| `compute_applicability_domain(...)` | Nearest-neighbor distance, sequence identity, chemistry/linkage novelty checks |
| `plot_predictions(...)` | Actual vs Predicted scatter with fit line |
| `plot_residuals(...)` | 3-panel: residuals vs predicted, histogram, Q-Q plot |
| `plot_feature_importance(...)` | Top-30 feature importance bar chart (tree models only) |
| `run_shap_analysis(...)` | SHAP TreeExplainer beeswarm plot + importance CSV |
| `write_validation_report(...)` | Full validation metrics report |
| `write_model_card(...)` | Standardized model documentation card |

### Applicability Domain Assessment

For each new prediction, checks:
- **Nearest-neighbor distance** in feature space (Euclidean after StandardScaler)
- **Maximum sequence identity** to any training ASO
- **Chemistry pattern novelty** — has this pattern been seen in training?
- **Linkage novelty** — has this linkage type been seen?

### Outputs
- `data/internal_cv_results.csv`
- `data/test_set_predictions.csv`
- `data/feature_importances.csv`
- `data/shap_importance.csv`
- `external_validation/external_validation_predictions.csv`
- `figures/test_actual_vs_predicted.png`
- `figures/test_residuals.png`
- `figures/feature_importance.png`
- `figures/shap_summary.png`
- `reports/validation_report.md`
- `reports/model_card.md`

---

## 11. `phase9_tofersen.py` — Tofersen Benchmarking

**Purpose:** Predict tofersen's inhibition and establish the comparison framework.

### Key Functions

| Function | What it does |
|---|---|
| `predict_tofersen()` | Builds tofersen features inline, aligns to training schema, predicts |
| `categorize_vs_tofersen(predicted, tofersen_pred)` | Returns "higher", "comparable", or "lower" |
| `save_tofersen_reference(tofersen_pred)` | Updates JSON with prediction, copies to web_tool dir |
| `write_tofersen_report(...)` | Documents prediction, conditions, and caveats |

### `predict_tofersen()` Detail

1. Loads preprocessing pipeline and best model from pkl files
2. Parses tofersen JSON (calls `phase6_external.parse_tofersen_json()`)
3. Gets standard conditions from training data medians/modes
4. Builds a single-row DataFrame with all required columns
5. Calls `phase4_features.build_feature_matrix()` with `skip_constant_removal=True`
6. Aligns to training feature columns (fills missing with 0)
7. Applies scaling if needed, runs prediction, clips to [0, 100]

**Result:** Tofersen predicted inhibition = **57.54%**

### Outputs
- `external_validation/tofersen_reference.json` (updated with prediction)
- `web_tool/tofersen_reference.json` (copy for web app)
- `reports/tofersen_prediction_report.md`

---

## Output Directory Structure

```
outputs/
├── data/
│   ├── cleaned_dataset.csv              # 2,155 cleaned records
│   ├── cleaned_dataset_clustered.csv    # + cluster IDs
│   ├── aso_sequences.fasta              # Unique sequences in FASTA
│   ├── sequence_clusters.tsv            # Cluster assignments
│   ├── feature_matrix.csv               # 2,155 × 339 features
│   ├── train.csv                        # 1,533 training records
│   ├── test.csv                         # 622 test records
│   ├── model_comparison_table.csv       # 9-model results
│   ├── test_set_predictions.csv         # Predictions + residuals
│   ├── internal_cv_results.csv          # 5-fold CV metrics
│   ├── feature_importances.csv          # LightGBM importances
│   ├── shap_importance.csv              # SHAP mean |values|
│   └── excluded_or_flagged_records.csv  # Flagged records
├── models/
│   ├── best_model.pkl                   # Trained LightGBM
│   └── preprocessing_pipeline.pkl       # Scaler + feature list + metadata
├── reports/
│   ├── raw_data_audit_report.md         # Phase 1
│   ├── column_dictionary.md             # Phase 1
│   ├── cleaning_report.md              # Phase 2
│   ├── redundancy_report.md            # Phase 3
│   ├── feature_dictionary.md           # Phase 4
│   ├── feature_extraction_report.md    # Phase 4
│   ├── split_report.md                 # Phase 5
│   ├── training_report.md              # Phase 7
│   ├── validation_report.md            # Phase 8
│   ├── model_card.md                   # Phase 8
│   └── tofersen_prediction_report.md   # Phase 9
├── figures/
│   ├── model_comparison.png            # R²/MAE/RMSE bar charts
│   ├── test_actual_vs_predicted.png    # Scatter plot
│   ├── test_residuals.png              # Residual analysis (3-panel)
│   ├── feature_importance.png          # Top-30 features
│   └── shap_summary.png               # SHAP beeswarm
├── external_validation/
│   ├── external_validation_raw.csv
│   ├── external_validation_clean.csv
│   ├── external_validation_feature_matrix.csv
│   ├── external_validation_predictions.csv
│   ├── tofersen_reference.json
│   ├── bibliography.bib
│   └── external_validation_curation_report.md
└── web_tool/
    ├── app.py                          # Streamlit web application
    ├── requirements.txt                # Python dependencies
    ├── tofersen_reference.json         # Copy for deployment
    ├── example_input.csv
    └── example_output.csv
```
