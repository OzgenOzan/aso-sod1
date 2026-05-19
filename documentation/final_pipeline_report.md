# SOD1 ASO Inhibition Prediction Pipeline — Final Report

**Generated:** 2026-05-13  
**Project:** PhD Thesis — In-Silico ML/DL Pipeline for ASO Knockdown Efficiency Prediction  
**Target Gene:** Human SOD1  
**Dataset:** ~2,264 SOD1-targeting antisense oligonucleotides

---

## 1. Executive Summary

This report documents the complete development of a reproducible machine learning pipeline for predicting the mRNA inhibition percentage of antisense oligonucleotides (ASOs) targeting the human SOD1 transcript. The pipeline processes 2,264 records containing 688 unique ASO sequences, engineers 335 features from sequence composition, chemical modifications, backbone linkage, and experimental conditions, and trains 9 regression models.

**Key Results:**
- **Best Model:** LightGBM (R² = 0.673, MAE = 12.4, RMSE = 15.9)
- **Runner-up:** XGBoost (R² = 0.667, MAE = 12.4)
- **Top predictive features:** ASO concentration, MOE fraction, cell density, treatment period, wing length
- **Tofersen benchmark:** Successfully integrated from FDA GSRS data

---

## 2. Available Dataset Columns

| Column | Standardized Name | Type | Description |
|---|---|---|---|
| ISIS | no | Integer | Unique ASO compound ID |
| Target_gene | target_gene | Categorical | Target gene (SOD-1) |
| Cell_line | cell_line | Categorical | HepG2, A431, SH-SY5Y |
| Density(cells/well) | density_cells_per_well | Numeric | Cell seeding density |
| Transfection | transfection | Categorical | electroporation, free uptake |
| ASO_volume(nM) | aso_concentration_nm | Numeric | ASO concentration |
| Treatment_Period(hours) | treatment_period_hours | Numeric | 16 or 24 hours |
| Primer_probe_set | primer_probe_set | Categorical | RTS3898, HTS90 |
| Sequence | sequence | Text | 5'->3' DNA sequence |
| Modification | modification | Categorical | MOE/cEt/5mc/deoxy |
| Location | location | Structured | Modification positions |
| Chemical_Pattern | chemical_pattern | Text | Per-position chemistry code |
| Linkage | linkage | Categorical | PS, PO/PS |
| Linkage_Location | linkage_location | Structured | Linkage positions |
| Smiles | smiles | Text | Full SMILES representation |
| Inhibition(%) | inhibition_percent | Numeric | Target variable (0-100) |
| seq_length | seq_length | Numeric | Sequence length |

Chemistry codes: `M` = MOE (2'-O-methoxyethyl), `C` = cEt (constrained ethyl), `d` = deoxy

---

## 3. Data Cleaning and Validation

### Phase 1 — Initial Audit
- **Total records:** 2,264
- **Unique sequences:** 688
- **Unique sequence+chemistry combos:** 1,008
- **Missing values:** 0 (all columns fully populated)
- **Target gene:** All SOD1 (confirmed)
- **Inhibition range:** 0–98% (within bounds)
- **Sequence length mismatches:** 12 records (stated 10, actual 17 — corrected)
- **Fully duplicated rows:** 20

### Phase 2 — Cleaning
- Column renaming to snake_case
- Sequence normalization (uppercase, ACGT validation)
- 12 sequence length corrections
- Transfection: "uptake" and "free uptake" merged to "free_uptake"
- 20 exact duplicate rows removed
- Biological replicates aggregated by mean inhibition
- **Final cleaned dataset:** 2,155 records
- SMILES validation: Limited (RDKit not available)

> [!NOTE]
> All sequences contain only standard DNA bases (A, C, G, T). No U->T conversion was needed.

---

## 4. Redundancy Filtering

### Method
Fallback pairwise comparison (CD-HIT-EST not available in environment):
- Equal-length: identity = matches / length
- Unequal-length: sliding window alignment
- Identity threshold: >90%
- Connected components (BFS) define clusters

### Results
- **Total unique sequences:** 688
- **Total pairwise comparisons:** 230,828
- **High-identity pairs (>90%):** 415
- **Clusters formed:** 396
- **Singletons:** Majority of clusters
- **Largest cluster:** Contains sequences differing by 1-2 nucleotides

> [!IMPORTANT]
> All members of the same sequence cluster are kept in the same train/test partition to prevent data leakage.

---

## 5. Feature Engineering

**Total features extracted: 335** (after removing 10 constant columns)

| Group | Features | Description |
|---|---|---|
| A. Sequence Composition | 282 | Base counts/fractions, di/trinucleotide frequencies, homopolymer runs, terminal motifs |
| B. Thermodynamic | 6 | Wallace Tm, GC-adjusted Tm, self-complementarity, palindrome score |
| C. Chemistry-Pattern | 23 | MOE/cEt/deoxy counts, gapmer architecture, wing lengths, transitions |
| D. Modification-Text | 6 | Contains_MOE, contains_cEt, contains_5mc, modification type count |
| E. Modification-Location | 5 | Position counts, min/max modified positions |
| F. Backbone-Linkage | 10 | PS/PO counts and fractions, linkage transitions |
| G. SMILES-Derived | 3 | String length, stereo count, approximate ring count (limited) |
| H. Experimental-Condition | 10 | Log-transformed concentration/density, one-hot cell line/transfection/primer |

### Skipped Features
- **SMILES molecular descriptors** (molecular weight, TPSA, etc.) — RDKit not available
- **SOD1 transcript mapping** (target position, UTR/CDS) — No coordinates in dataset

### Thermodynamic Approximations
- **Wallace Tm:** Tm = 2(A+T) + 4(G+C) — for short oligos
- **GC-adjusted Tm:** Tm = 64.9 + 41(GC_count - 16.4) / length

---

## 6. Train/Test Split

| Set | Rows | Clusters | % of Total |
|---|---|---|---|
| Train | 1,533 | 277 | 71.1% |
| Test | 622 | 119 | 28.9% |
| Total | 2,155 | 396 | 100% |

- **Split unit:** sequence_cluster_id (cluster-aware)
- **Stratification:** By median inhibition bin
- **Leakage control:** Zero cluster overlap between train and test

---

## 7. External Validation Curation

- **Tofersen** (Qalsody) curated from provided FDA GSRS JSON
  - Sequence: CAGGATACATTTCTACAGCT (DNA equivalent)
  - Chemistry: 5-10-5 MOE gapmer
  - CAS: 1898254-60-8, UNII: 5YL205692C
  - Source: FDA GSRS, USAN Council, Qalsody prescribing information
- Additional ASOs documented for future manual curation (sequences not fabricated)

---

## 8. Model Training

### Model Comparison (Test Set)

| Model | MAE | RMSE | R² | Pearson r | Within ±10 | Within ±15 |
|---|---|---|---|---|---|---|
| **LightGBM** | **12.42** | **15.86** | **0.673** | **0.830** | **48.4%** | **68.0%** |
| XGBoost | 12.35 | 16.00 | 0.667 | 0.825 | 50.8% | 69.0% |
| Gradient Boosting | 12.61 | 16.60 | 0.641 | 0.820 | 50.3% | 67.7% |
| Ridge | 15.95 | 20.30 | 0.463 | 0.708 | 42.3% | 54.3% |
| Linear Regression | 16.10 | 20.67 | 0.444 | 0.700 | 40.7% | 54.0% |
| Random Forest | 17.20 | 20.94 | 0.429 | 0.703 | 33.6% | 47.7% |
| SVR | 19.31 | 24.44 | 0.222 | 0.535 | 32.0% | 49.0% |
| Elastic Net | 17.88 | 22.29 | 0.353 | 0.620 | 33.1% | 49.2% |
| Mean Predictor | 24.56 | 28.48 | -0.056 | N/A | 21.2% | 31.0% |

### Cross-Validation (5-Fold, Cluster-Aware)

| Metric | Mean +/- Std |
|---|---|
| MAE | 16.63 +/- 6.67 |
| RMSE | 21.03 +/- 7.31 |
| R² | 0.369 +/- 0.360 |

> [!NOTE]
> Cross-validation R² variance is expected due to cluster-aware splitting with uneven cluster sizes.

---

## 9. Internal and External Validation

### Top Features (SHAP Analysis)

| Rank | Feature | Mean |SHAP| |
|---|---|---|
| 1 | log10_aso_concentration_nm | 7.38 |
| 2 | fraction_M (MOE fraction) | 4.97 |
| 3 | log10_density_cells_per_well | 3.81 |
| 4 | treatment_period_hours | 3.68 |
| 5 | left_modified_wing_length | 3.19 |
| 6 | tri_ATG | 1.92 |
| 7 | C_count | 1.21 |
| 8 | tri_AAT | 1.13 |
| 9 | tri_TTA | 0.98 |
| 10 | max_gc_stretch | 0.89 |

### Key Interpretations
1. **ASO concentration** is the strongest predictor — higher doses yield higher inhibition
2. **MOE fraction** reflects gapmer wing chemistry importance
3. **Cell density and treatment period** significantly influence outcomes
4. **Wing length** captures gapmer architecture effects
5. **Trinucleotide motifs** (ATG, AAT, TTA) may reflect target site preferences

---

## 10. Tofersen Benchmark

> **Disclaimer:** This is an in-silico benchmark and not a clinical efficacy claim.

- **Tofersen sequence:** CAGGATACATTTCTACAGCT
- **Chemistry:** 5-10-5 MOE gapmer with PS/PO backbone
- **Predicted inhibition (standardized conditions):** 57.54%
- **Standardized conditions:** HepG2, electroporation, 3000 nM, 16 h, 20000 cells/well
- **Benchmark framework:**
  - Comparable: |delta| <= 10 percentage points
  - Higher: delta > +10
  - Lower: delta < -10

---

## 11. Web Tool Implementation

A Streamlit-based web application ([app.py](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/web_tool/app.py)) has been generated with:
- ASO sequence and chemistry pattern input
- Experimental condition selection
- Real-time inhibition prediction
- Tofersen benchmark comparison
- Input validation
- Applicability domain warnings

**Launch:** `streamlit run pipeline/outputs/web_tool/app.py`

---

## 12. Limitations

1. **No transcript-coordinate features** — Dataset lacks target position on SOD1 mRNA
2. **No RNA accessibility features** — Would require external RNA structure prediction
3. **Limited SMILES features** — RDKit not available for molecular descriptors
4. **Cell-line specificity** — Only 3 cell lines (HepG2, A431, SH-SY5Y)
5. **Single target gene** — Model may not generalize beyond SOD1
6. **In-vitro only** — Does not predict in-vivo efficacy or pharmacokinetics
7. **MLP model** — PyTorch DLL loading failed on this Windows environment
8. **External validation** — Only tofersen has verified sequence from provided data

---

## 13. Reproducibility Checklist

| Item | Status |
|---|---|
| Random seed set (42) | Done |
| Raw data preserved | Done |
| All intermediate files saved | Done |
| Cluster-aware train/test split | Done |
| No test/external leakage | Done |
| Predictions clipped to [0, 100] | Done |
| All assumptions documented | Done |
| Model and scaler saved as PKL | Done |
| Feature extraction code modular | Done |
| Reports in markdown format | Done |

---

## 14. Output File Inventory

### Data Files
| File | Description |
|---|---|
| [cleaned_dataset.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/cleaned_dataset.csv) | 2,155 cleaned records |
| [cleaned_dataset_clustered.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/cleaned_dataset_clustered.csv) | With cluster IDs |
| [feature_matrix.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/feature_matrix.csv) | 2,155 x 339 features |
| [train.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/train.csv) | 1,533 training records |
| [test.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/test.csv) | 622 test records |
| [model_comparison_table.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/model_comparison_table.csv) | 9-model comparison |
| [test_set_predictions.csv](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/data/test_set_predictions.csv) | Test predictions + residuals |

### Reports
| File | Description |
|---|---|
| [raw_data_audit_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/raw_data_audit_report.md) | Phase 1 audit |
| [cleaning_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/cleaning_report.md) | Phase 2 cleaning |
| [redundancy_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/redundancy_report.md) | Phase 3 clustering |
| [feature_dictionary.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/feature_dictionary.md) | Feature catalog |
| [split_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/split_report.md) | Train/test balance |
| [training_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/training_report.md) | Model training details |
| [validation_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/validation_report.md) | Validation metrics |
| [model_card.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/model_card.md) | Model documentation |
| [tofersen_prediction_report.md](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/reports/tofersen_prediction_report.md) | Tofersen benchmark |

### Figures
| File | Description |
|---|---|
| [model_comparison.png](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/figures/model_comparison.png) | Model R²/MAE/RMSE bars |
| [test_actual_vs_predicted.png](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/figures/test_actual_vs_predicted.png) | Scatter plot |
| [test_residuals.png](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/figures/test_residuals.png) | Residual analysis |
| [feature_importance.png](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/figures/feature_importance.png) | LightGBM importances |
| [shap_summary.png](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/figures/shap_summary.png) | SHAP beeswarm plot |

### Models
| File | Description |
|---|---|
| [best_model.pkl](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/models/best_model.pkl) | LightGBM model |
| [preprocessing_pipeline.pkl](file:///c:/Users/oozgen/Desktop/Udemy/Gemini/Codebase/20260513_ASO_4/pipeline/outputs/models/preprocessing_pipeline.pkl) | Scaler + feature list |

---

## 15. References

1. Miller TM, et al. (2013). An antisense oligonucleotide against SOD1 delivered intrathecally. *Lancet Neurology*, 12(5), 435-442.
2. Smith RA, et al. (2006). Antisense oligonucleotide therapy for neurodegenerative disease. *J Clin Invest*, 116(8), 2290-2296.
3. FDA GSRS. Tofersen Sodium. UNII: 5YL205692C.
4. Biogen Inc. QALSODY (tofersen) Prescribing Information. 2023.
5. Bennett CF, et al. WO2007002390. Compositions and methods for modulation of SOD1 expression.
