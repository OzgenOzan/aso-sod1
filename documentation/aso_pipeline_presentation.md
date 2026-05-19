# SOD1 ASO Inhibition Prediction Pipeline
## A Complete Step-by-Step Guide

> **Audience:** Anyone — no prior knowledge of genetics, cell biology, or programming is assumed.
> By the end of this guide you will understand every component and be able to recreate the entire codebase from scratch.

---

## Part 1 — Background: The Science

### 1.1 What is ALS?

**Amyotrophic Lateral Sclerosis (ALS)**, also known as Lou Gehrig's disease, is a progressive neurodegenerative disease that attacks motor neurons — the nerve cells in the brain and spinal cord that control voluntary muscle movement. Patients gradually lose the ability to move, speak, eat, and eventually breathe.

### 1.2 The SOD1 Gene

About 2% of ALS cases are caused by mutations in the **SOD1** gene (Superoxide Dismutase 1). This gene normally produces an enzyme that protects cells from damage. When mutated, it produces a toxic protein that kills motor neurons.

**The therapeutic strategy:** If we can stop the SOD1 gene from producing this toxic protein, we can slow down ALS progression.

### 1.3 What is an ASO?

An **Antisense Oligonucleotide (ASO)** is a short, synthetic strand of modified DNA (typically 15–25 nucleotides long) designed to bind to a specific messenger RNA (mRNA) and prevent it from being translated into protein.

```
How ASOs Work (simplified):

DNA  ──transcription──►  mRNA  ──translation──►  Protein (toxic SOD1)
                           ▲
                           │
                      ASO binds here
                      and triggers
                      degradation
                           │
                           ▼
                      mRNA destroyed
                      = less toxic protein
```

### 1.4 What is a Gapmer?

A **gapmer** is a specific ASO design with three regions:

```
 5'─[WING]──[GAP]──[WING]─3'
     MOE    deoxy    MOE
   modified  DNA   modified
   (stable) (active) (stable)
```

- **Wings** (flanking regions): Chemically modified nucleotides (e.g., MOE = 2'-O-methoxyethyl) that protect the ASO from being broken down by enzymes.
- **Gap** (central region): Unmodified DNA nucleotides that recruit an enzyme called **RNase H** to cut the target mRNA.

A common design is **5-10-5**: 5 modified nucleotides on each wing, 10 deoxy nucleotides in the gap.

### 1.5 What is Tofersen (Qalsody)?

**Tofersen** is the first FDA-approved ASO drug for SOD1-ALS (approved April 2023, brand name Qalsody). It is a **5-10-5 MOE gapmer** with the sequence `CAGGATACATTTCTACAGCT`. In this project, tofersen serves as a **benchmark** — we compare our predicted ASOs against it.

### 1.6 What is "Inhibition Percentage"?

When scientists test an ASO in the lab, they measure how much it reduces ("knocks down") the target mRNA. This is expressed as **inhibition percentage (0–100%)**:
- **0%** = No effect (mRNA levels unchanged)
- **100%** = Complete knockdown (all mRNA destroyed)
- **≥70%** = Generally considered "high efficacy"

---

## Part 2 — Project Goal

### What Does This Pipeline Do?

This project builds a **machine learning model** that predicts how effective a new ASO will be at knocking down SOD1 mRNA, **without** needing to synthesize and test it in the lab.

```
Input:  ASO sequence + chemistry + experimental conditions
  ↓
Model:  Machine learning (LightGBM)
  ↓
Output: Predicted inhibition percentage (0–100%)
        + comparison against tofersen
```

### Why Is This Useful?

Synthesizing and testing ASOs is expensive and time-consuming. A predictive model lets researchers:
1. Screen thousands of candidate ASOs computationally
2. Prioritize the most promising ones for lab testing
3. Understand which features make an ASO effective

---

## Part 3 — Project Structure

```
20260513_ASO_4/
├── dataset/
│   ├── sod-1_data.xlsx          # Raw experimental data (2,264 ASO records)
│   └── tofersen.json            # FDA GSRS data for tofersen
├── explore_data.py              # Quick data exploration script
├── verify_pipeline.py           # Pipeline verification script
├── final_pipeline_report.md     # Summary report
└── pipeline/
    ├── config.py                # Central configuration (paths, constants)
    ├── run_pipeline.py          # Master runner (executes all phases)
    ├── phase1_audit.py          # Data loading & initial audit
    ├── phase2_cleaning.py       # Data cleaning & validation
    ├── phase3_clustering.py     # Sequence redundancy clustering
    ├── phase4_features.py       # Feature extraction (335 features)
    ├── phase5_split.py          # Train/test splitting
    ├── phase6_external.py       # External validation curation
    ├── phase7_modeling.py       # Model training (9 models)
    ├── phase8_validation.py     # Model evaluation & SHAP analysis
    ├── phase9_tofersen.py       # Tofersen benchmarking
    └── outputs/
        ├── data/                # CSV outputs (cleaned data, features, splits)
        ├── models/              # Saved model files (.pkl)
        ├── reports/             # Markdown reports for each phase
        ├── figures/             # Plots (model comparison, SHAP, etc.)
        ├── external_validation/ # Tofersen reference data
        └── web_tool/            # Streamlit web application
            ├── app.py
            ├── requirements.txt
            └── tofersen_reference.json
```

---

## Part 4 — Environment Setup

### 4.1 Prerequisites

You need **Python 3.10+** installed. Then install the required packages:

```bash
pip install pandas numpy scikit-learn xgboost lightgbm shap matplotlib seaborn streamlit openpyxl scipy
```

### 4.2 Dataset

Place the raw dataset file `sod-1_data.xlsx` in the `dataset/` folder. This Excel file contains 2,264 rows × 17 columns of experimentally measured ASO data. Also place the `tofersen.json` file (from the FDA GSRS database) in the same folder.

### 4.3 Running the Full Pipeline

```bash
cd pipeline
python run_pipeline.py
```

This sequentially executes Phases 1–9. Total runtime is approximately 2–5 minutes.

---

## Part 5 — The Configuration File (`config.py`)

This file is the **single source of truth** for all settings. Every other script imports from here.

**Key settings:**

| Setting | Value | Purpose |
|---|---|---|
| `RANDOM_SEED` | 42 | Ensures reproducibility — same results every run |
| `TEST_SIZE` | 0.30 | 30% of data reserved for testing |
| `CLUSTER_IDENTITY_THRESHOLD` | 0.90 | Sequences >90% similar are grouped together |
| `INHIBITION_BINS` | [0,20,40,60,80,100] | Bins for stratified splitting |
| `TOFERSEN_COMPARABLE_MARGIN` | 10.0 | ±10% margin for tofersen comparison |

**Column mapping** renames raw Excel column names (e.g., `ASO_volume(nM)`) to clean Python-friendly names (e.g., `aso_concentration_nm`).

**Chemistry codes:** `M` = MOE, `C` = cEt, `d` = deoxy.

---

## Part 6 — Phase 1: Data Audit (`phase1_audit.py`)

### What It Does
Loads the raw Excel file and performs a comprehensive quality check **without modifying** the data.

### Key Checks
1. **Row/column counts** — 2,264 rows × 17 columns
2. **Missing values** — 0 (all columns fully populated)
3. **Target gene validation** — Confirms all records target SOD1
4. **Inhibition range** — All values within 0–98% (valid)
5. **Sequence length verification** — Finds 12 records where stated length (10) differs from actual length (17)
6. **Duplicate detection** — Finds 20 fully duplicated rows
7. **Chemical pattern catalog** — Lists all gapmer architectures observed

### How the Code Works

```python
def load_raw_data():
    df = pd.read_excel(RAW_EXCEL)  # Load Excel file into a DataFrame
    return df

def audit_dataset(df):
    audit = {}
    audit["n_rows"] = len(df)                          # Count rows
    audit["n_unique_sequences"] = df["Sequence"].nunique()  # Count unique ASOs
    # ... more checks ...
    return audit
```

### Output Files
- `reports/raw_data_audit_report.md` — Full audit findings
- `reports/column_dictionary.md` — Description of every column

---

## Part 7 — Phase 2: Data Cleaning (`phase2_cleaning.py`)

### What It Does
Normalizes, validates, and deduplicates the dataset.

### Processing Steps

| Step | Action | Example |
|---|---|---|
| 1 | Rename columns to snake_case | `ASO_volume(nM)` → `aso_concentration_nm` |
| 2 | Uppercase all sequences | `acgt` → `ACGT` |
| 3 | Fix sequence length mismatches | 12 records corrected from 10 → 17 |
| 4 | Normalize transfection labels | `"uptake"` and `"free uptake"` → `"free_uptake"` |
| 5 | Validate numeric ranges | Concentration, density, inhibition |
| 6 | Validate chemical patterns | Check length matches sequence, valid codes only |
| 7 | Parse Location column | Extract modification position indices |
| 8 | Parse Linkage_Location | Extract linkage position indices |
| 9 | SMILES validation | Check non-empty (RDKit not available) |
| 10 | Remove 20 exact duplicates | Drop identical rows |
| 11 | Aggregate biological replicates | Same conditions → mean inhibition |

### Key Code Pattern — Biological Replicate Aggregation

When the same ASO is tested multiple times under identical conditions, results are averaged:

```python
condition_cols = ["sequence", "chemical_pattern", "linkage", "cell_line",
                  "transfection", "aso_concentration_nm", ...]

grouped = df.groupby(condition_cols)
for name, group in grouped:
    rec = group.iloc[0].to_dict()
    if len(group) > 1:
        rec["inhibition_percent"] = group["inhibition_percent"].mean()
    agg_records.append(rec)
```

### Output
- `data/cleaned_dataset.csv` — **2,155 cleaned records** (down from 2,264)

---

## Part 8 — Phase 3: Sequence Clustering (`phase3_clustering.py`)

### Why Cluster Sequences?

If two ASO sequences are nearly identical (e.g., differ by only 1 nucleotide), they likely target the same region and have similar effects. If one goes into training and the other into testing, the model would "cheat" by memorizing similar sequences. This is called **data leakage**.

### How It Works

1. **Extract unique sequences** — 688 unique ASOs
2. **Compare all pairs** — Calculate nucleotide identity

```python
def compute_pairwise_identity(seq1, seq2):
    # For equal-length sequences
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / len(seq1)
    # Example: ACGTACGT vs ACGTACGA = 7/8 = 87.5% identity
```

3. **Group by >90% identity** — Using connected-components (BFS graph algorithm)
4. **Result:** 396 clusters from 688 sequences

### What "90% Identity" Means for ASOs
- For a 20-nucleotide ASO: 2 mismatches = 90% identity
- For a 17-nucleotide ASO: 2 mismatches = 88.2% (below threshold, separate clusters)

### Output
- `data/cleaned_dataset_clustered.csv` — Dataset with `sequence_cluster_id` column
- `data/sequence_clusters.tsv` — Cluster assignments

---

## Part 9 — Phase 4: Feature Extraction (`phase4_features.py`)

### What Are Features?

Features are **numerical measurements** extracted from each ASO that the model uses to learn patterns. Think of them as "descriptors" — they translate complex biological information into numbers a computer can process.

### Feature Groups (335 Total)

#### A. Sequence Composition (282 features)
```
Base counts:     A_count=5, C_count=6, G_count=5, T_count=4
Base fractions:  A_fraction=0.25, GC_fraction=0.55
Dinucleotides:   di_CG=2 (CpG count), di_GG=1 ...
Trinucleotides:  tri_ATG=1, tri_CCA=0 ...
Terminal bases:  five_prime_base_C=1 (starts with C)
Homopolymers:    max_homopolymer_run=3 (e.g., "AAA")
```

#### B. Thermodynamic (6 features)
```
Tm_wallace:           Melting temperature estimate = 2×(A+T) + 4×(G+C)
self_complementarity: Can the ASO fold back on itself?
palindrome_score:     Does it contain palindromic motifs?
terminal_gc_clamp:    Does it have G/C at both ends? (stabilizing)
```

#### C. Chemistry Pattern (23 features)
```
count_M=5, count_d=10         # How many MOE vs deoxy positions
fraction_M=0.25               # Proportion of MOE
left_modified_wing_length=5   # Left wing size
central_deoxy_gap_length=10   # Gap size
gapmer_5-10-5=1               # One-hot: is it a 5-10-5 design?
symmetry_flag=1               # Are wings equal length?
```

#### D–H. Other Features
- **D. Modification text** (6): Binary flags for MOE, cEt, 5-methylcytosine
- **E. Location** (5): Modification position indices
- **F. Linkage** (10): PS/PO backbone chemistry counts and fractions
- **G. SMILES** (3): Basic molecular string features
- **H. Experimental conditions** (10): Log-concentration, cell line, transfection (one-hot encoded)

### Key Code — Feature Extraction Pattern

```python
def extract_sequence_features(df):
    feats = pd.DataFrame(index=df.index)
    seqs = df["sequence"].values

    feats["seq_length"] = [len(s) for s in seqs]
    for base in "ACGT":
        feats[f"{base}_count"] = [s.count(base) for s in seqs]
        feats[f"{base}_fraction"] = feats[f"{base}_count"] / feats["seq_length"]

    # GC content — one of the most important biological features
    feats["GC_fraction"] = (feats["G_count"] + feats["C_count"]) / feats["seq_length"]
    return feats
```

### Output
- `data/feature_matrix.csv` — 2,155 rows × 339 columns

---

## Part 10 — Phase 5: Train/Test Split (`phase5_split.py`)

### Why Not Just Random Split?

A random split might put similar sequences in both train and test sets (leakage). Instead, we split by **cluster** — all sequences in the same cluster go to the same set.

### How It Works

1. Calculate median inhibition per cluster
2. Bin clusters into 5 inhibition ranges (0–20, 20–40, ..., 80–100)
3. Within each bin, randomly assign 70% of clusters to train, 30% to test
4. **Verify zero cluster overlap** between sets

```python
def cluster_aware_stratified_split(df, test_size=0.30, seed=42):
    for bin_label in INHIBITION_BIN_LABELS:
        bin_clusters = ...  # clusters in this inhibition range
        n_test = max(1, int(len(bin_clusters) * test_size))
        test_clusters.extend(bin_clusters[:n_test])
        train_clusters.extend(bin_clusters[n_test:])

    assert len(test_clusters & train_clusters) == 0  # No leakage!
```

### Result

| Set | Rows | Clusters | % |
|---|---|---|---|
| Train | 1,533 | 277 | 71.1% |
| Test | 622 | 119 | 28.9% |

---

## Part 11 — Phase 6: External Validation (`phase6_external.py`)

### What It Does

Curates **tofersen** (the FDA-approved ASO drug) from the provided FDA GSRS JSON file as an external benchmark.

### Tofersen Extraction

The code parses the FDA's structured JSON to extract:
- **Sequence:** Built from the `nucleicAcid.subunits` field, converting RNA (U) to DNA (T)
- **Chemical pattern:** Built from `nucleicAcid.sugars` — maps MOE and dR positions to `M` and `d`
- **Linkage:** From `nucleicAcid.linkages` — identifies PS and PO positions

```python
def parse_tofersen_json():
    subunits = data["nucleicAcid"]["subunits"]
    sequence = subunits[0]["sequence"].replace("U", "T")

    chem_pattern = ['d'] * len(sequence)  # Default: deoxy
    for sugar in sugars:
        if sugar["sugar"] == "MOE":
            for site in sugar["sites"]:
                chem_pattern[site["residueIndex"] - 1] = 'M'

    return tofersen_record
```

### Missing Experimental Conditions

Tofersen's JSON has no lab conditions (cell line, concentration, etc.). The pipeline imputes these using the **most common values** from the training dataset (e.g., HepG2 cell line, 3000 nM, 16 hours).

---

## Part 12 — Phase 7: Model Training (`phase7_modeling.py`)

### Models Trained (9 total)

| # | Model | Type | Key Parameters |
|---|---|---|---|
| 1 | Mean Predictor | Baseline | Predicts the training mean for everything |
| 2 | Linear Regression | Linear | OLS on scaled features |
| 3 | Ridge | Linear | L2 regularization, α=10 |
| 4 | Elastic Net | Linear | L1+L2, α=1.0, l1_ratio=0.5 |
| 5 | Random Forest | Ensemble | 500 trees, max_depth=15 |
| 6 | Gradient Boosting | Ensemble | 500 trees, max_depth=5, lr=0.05 |
| 7 | XGBoost | Ensemble | 500 trees, max_depth=6, lr=0.05 |
| 8 | LightGBM | Ensemble | 500 trees, max_depth=6, lr=0.05 |
| 9 | SVR | Kernel | RBF kernel, C=10, ε=2.0 |

### Feature Scaling

Some models (linear, SVR) require features to be on the same scale. A `StandardScaler` is applied:

```python
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)   # Learn mean/std from train
X_test_scaled = scaler.transform(X_test)          # Apply same transform to test
```

> [!IMPORTANT]
> Tree-based models (Random Forest, XGBoost, LightGBM) use **unscaled** data because they split on thresholds and don't care about scale.

### Results

| Model | R² | MAE | RMSE |
|---|---|---|---|
| **LightGBM** | **0.673** | **12.42** | **15.86** |
| XGBoost | 0.667 | 12.35 | 16.00 |
| Gradient Boosting | 0.641 | 12.61 | 16.60 |
| Ridge | 0.463 | 15.95 | 20.30 |
| Linear Regression | 0.444 | 16.10 | 20.67 |
| Random Forest | 0.429 | 17.20 | 20.94 |
| SVR | 0.222 | 19.31 | 24.44 |
| Elastic Net | 0.353 | 17.88 | 22.29 |
| Mean Predictor | -0.056 | 24.56 | 28.48 |

**Winner: LightGBM** with R²=0.673, meaning it explains ~67% of the variance in inhibition.

### What Gets Saved

```python
# The best model itself
pickle.dump(best_model, open("models/best_model.pkl", "wb"))

# The preprocessing pipeline (scaler + feature list + metadata)
pickle.dump({
    "scaler": scaler,
    "feature_cols": feature_cols,        # Which features, in order
    "best_model_name": "LightGBM",
    "best_model_key": "lightgbm",
    "needs_scaling": False,              # LightGBM doesn't need scaling
}, open("models/preprocessing_pipeline.pkl", "wb"))
```

---

## Part 13 — Phase 8: Validation (`phase8_validation.py`)

### What It Does

Rigorously evaluates the best model through:

1. **Test set evaluation** — Predictions on held-out data, scatter plots, residual analysis
2. **Cross-validation** — 5-fold cluster-aware CV to check stability
3. **Feature importance** — Which features matter most (built-in + SHAP)
4. **Applicability domain** — How "similar" is a new ASO to the training data?
5. **External validation** — Predict tofersen and compare

### Top 5 Most Important Features (SHAP)

| Rank | Feature | Meaning |
|---|---|---|
| 1 | `log10_aso_concentration_nm` | Higher drug dose → higher inhibition |
| 2 | `fraction_M` (MOE fraction) | Wing chemistry matters |
| 3 | `log10_density_cells_per_well` | Lab conditions affect results |
| 4 | `treatment_period_hours` | Longer treatment → more knockdown |
| 5 | `left_modified_wing_length` | Gapmer architecture effects |

### Output Figures
- Actual vs. Predicted scatter plot
- Residual distribution
- Feature importance bar chart
- SHAP beeswarm plot

---

## Part 14 — Phase 9: Tofersen Benchmark (`phase9_tofersen.py`)

### What It Does

Uses the trained model to predict tofersen's inhibition under standardized conditions, then sets up the comparison framework.

### Standardized Conditions

| Parameter | Value |
|---|---|
| Cell line | HepG2 |
| Transfection | electroporation |
| Concentration | 3000 nM |
| Treatment | 16 hours |
| Cell density | 20,000 cells/well |
| Primer/probe | RTS3898 |

### Result
- **Tofersen predicted inhibition:** 57.54%

### Comparison Framework

Any new ASO is categorized relative to tofersen:
- **Higher:** delta > +10 percentage points
- **Comparable:** |delta| ≤ 10
- **Lower:** delta < -10

---

## Part 15 — The Web Application (`app.py`)

### What It Does

A **Streamlit** web interface that lets users input ASO parameters and get instant predictions.

### How to Run

```bash
streamlit run pipeline/outputs/web_tool/app.py
```

### Architecture

```
User Input (sidebar)          Processing              Output
─────────────────          ──────────────          ──────────
Sequence                →  Feature extraction  →  Predicted inhibition %
Chemical pattern        →  (inline, 335 features) Tofersen comparison
Modification            →  Model prediction     →  Delta category
Linkage                 →  (loaded from .pkl)      Detailed results table
Cell line, conc, etc.
```

### Key Design Decision — Inline Feature Extraction

The web app contains its own `extract_all_features()` function (not imported from phase4) so it can be deployed independently without the pipeline directory.

### Path Resolution

The app loads model files from sibling directories:

```python
APP_DIR = os.path.dirname(os.path.abspath(__file__))    # .../web_tool/
MODEL_DIR = os.path.join(APP_DIR, "..", "models")        # .../outputs/models/
EXT_VAL_DIR = os.path.join(APP_DIR, "..", "external_validation")
```

---

## Part 16 — How to Recreate From Scratch

### Step 1: Set Up the Project

```bash
mkdir aso_pipeline && cd aso_pipeline
mkdir dataset pipeline pipeline/outputs
mkdir pipeline/outputs/{data,models,reports,figures,external_validation,web_tool}
```

### Step 2: Install Dependencies

```bash
pip install pandas numpy scikit-learn xgboost lightgbm shap matplotlib seaborn streamlit openpyxl scipy
```

### Step 3: Create `pipeline/config.py`

Define all paths, constants, column mappings, and chemistry codes as shown in Part 5. The key constants are:
- `RANDOM_SEED = 42`
- `TEST_SIZE = 0.30`
- `CLUSTER_IDENTITY_THRESHOLD = 0.90`
- `COLUMN_MAP` dictionary mapping raw → standardized column names

### Step 4: Create Phase Scripts (1–9)

Create each `phaseN_*.py` file following the logic described in Parts 6–14. Each script:
1. Imports from `config.py`
2. Has a `main()` function
3. Reads input from the previous phase's output
4. Writes its own outputs (data + reports)

### Step 5: Create `pipeline/run_pipeline.py`

```python
from phase1_audit import main as phase1
from phase2_cleaning import main as phase2
# ... import all phases ...

def run_pipeline():
    phase1()
    phase2()
    # ... call all phases sequentially ...

if __name__ == "__main__":
    run_pipeline()
```

### Step 6: Create the Web App

Create `pipeline/outputs/web_tool/app.py` with:
1. Inline feature extraction function (mirrors phase4 logic)
2. Model/pipeline loading from `../models/`
3. Streamlit UI with sidebar inputs and result cards
4. Input validation and tofersen benchmarking

### Step 7: Run Everything

```bash
cd pipeline
python run_pipeline.py          # Train the model
streamlit run outputs/web_tool/app.py  # Launch the web app
```

---

## Part 17 — Key Concepts Glossary

| Term | Meaning |
|---|---|
| **ASO** | Short synthetic DNA that silences a gene |
| **Gapmer** | ASO with modified wings and a deoxy gap |
| **MOE** | 2'-O-methoxyethyl — a common wing modification |
| **cEt** | Constrained ethyl — another wing modification |
| **PS** | Phosphorothioate — backbone modification for stability |
| **PO** | Phosphodiester — natural backbone linkage |
| **RNase H** | Enzyme that cuts RNA in a DNA-RNA hybrid |
| **Inhibition %** | How much an ASO reduces target mRNA (0–100) |
| **R²** | How much variance the model explains (0–1, higher = better) |
| **MAE** | Mean Absolute Error — average prediction error |
| **RMSE** | Root Mean Squared Error — penalizes large errors more |
| **SHAP** | Method to explain which features drive each prediction |
| **Data leakage** | When test data information leaks into training |
| **Feature engineering** | Creating numerical inputs from raw data |
| **One-hot encoding** | Converting categories to binary columns |
| **Cluster-aware split** | Keeping similar sequences in the same partition |
| **Cross-validation** | Testing model stability across multiple data splits |
| **Applicability domain** | How "similar" a new input is to training data |

---

## Part 18 — Limitations and Future Work

### Current Limitations
1. **No transcript coordinates** — We don't know where on the SOD1 mRNA each ASO binds
2. **No RNA structure features** — Secondary structure accessibility is not modeled
3. **Limited SMILES features** — RDKit not available for molecular descriptors
4. **Only 3 cell lines** — HepG2, A431, SH-SY5Y
5. **Single target gene** — Trained only on SOD1
6. **In-vitro only** — Does not predict in-vivo efficacy

### Possible Improvements
- Add RNA accessibility features (requires external tools)
- Install RDKit for full molecular descriptor extraction
- Map ASO sequences to SOD1 transcript coordinates
- Expand to multi-target datasets
- Implement deep learning models (transformer-based)
- Add dose-response curve modeling

---

> [!NOTE]
> **Reproducibility:** Every run with `RANDOM_SEED=42` produces identical results. All intermediate files are saved, and every assumption is documented in phase-specific reports.
