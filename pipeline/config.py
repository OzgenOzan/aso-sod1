"""
Central configuration for the ASO inhibition prediction pipeline.
All paths, constants, and reproducibility settings are defined here.
"""

import os

# -- Reproducibility --------------------------------------------------
RANDOM_SEED = 42
TEST_SIZE = 0.30
CLUSTER_IDENTITY_THRESHOLD = 0.90

# -- Paths ------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASET_DIR = os.path.join(BASE_DIR, "dataset")
OUTPUT_DIR = os.path.join(BASE_DIR, "pipeline", "outputs")

RAW_EXCEL = os.path.join(DATASET_DIR, "sod-1_data.xlsx")
TOFERSEN_JSON = os.path.join(DATASET_DIR, "tofersen.json")

DATA_DIR = os.path.join(OUTPUT_DIR, "data")
REPORT_DIR = os.path.join(OUTPUT_DIR, "reports")
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures")
MODEL_DIR = os.path.join(OUTPUT_DIR, "models")
EXT_VAL_DIR = os.path.join(OUTPUT_DIR, "external_validation")
WEB_TOOL_DIR = os.path.join(OUTPUT_DIR, "web_tool")

# -- Column name mapping (raw -> standardized) ------------------------
COLUMN_MAP = {
    "ISIS": "no",                              # Note: dataset uses ISIS, not "No"
    "Target_gene": "target_gene",
    "Cell_line": "cell_line",
    "Density(cells/well)": "density_cells_per_well",
    "Transfection": "transfection",
    "ASO_volume(nM)": "aso_concentration_nm",
    "Treatment_Period(hours)": "treatment_period_hours",
    "Primer_probe_set": "primer_probe_set",
    "Sequence": "sequence",
    "Modification": "modification",
    "Location": "location",
    "Chemical_Pattern": "chemical_pattern",
    "Linkage": "linkage",
    "Linkage_Location": "linkage_location",
    "Smiles": "smiles",
    "Inhibition(%)": "inhibition_percent",
    "seq_length": "seq_length",
}

# -- Chemistry code meanings ------------------------------------------
CHEMISTRY_CODES = {
    "M": "2'-MOE (2'-O-methoxyethyl) modified nucleotide",
    "C": "cEt (constrained ethyl / S-cEt BNA) modified nucleotide",
    "d": "2'-deoxy (unmodified DNA) nucleotide",
}

# -- Inhibition bins for stratification -------------------------------
INHIBITION_BINS = [0, 20, 40, 60, 80, 100]
INHIBITION_BIN_LABELS = ["0-20", "20-40", "40-60", "60-80", "80-100"]

# -- High-efficacy threshold ------------------------------------------
HIGH_EFFICACY_THRESHOLD = 70.0

# -- Tofersen benchmark -----------------------------------------------
TOFERSEN_COMPARABLE_MARGIN = 10.0  # +/-10 inhibition percentage points

# Ensure output directories exist
for d in [DATA_DIR, REPORT_DIR, FIGURE_DIR, MODEL_DIR, EXT_VAL_DIR, WEB_TOOL_DIR]:
    os.makedirs(d, exist_ok=True)
