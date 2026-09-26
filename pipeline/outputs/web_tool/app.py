"""
ASO Inhibition Predictor -- Streamlit Web Application
======================================================
Predicts inhibition percentage for novel SOD1-targeting ASOs
and compares them against an in-silico tofersen model reference
(model prediction under imputed conditions; not experimental validation).

Usage: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import pickle
import json
import os
import sys
import hashlib
import warnings

# -- Page config ------------------------------------------------------
st.set_page_config(
    page_title="SOD1 ASO Inhibition Predictor",
    page_icon="\U0001f9ec",
    layout="wide",
    initial_sidebar_state="expanded",
)

# -- Paths ------------------------------------------------------------
APP_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.join(APP_DIR, "..", "..", "..")  # project root / pipeline
MODEL_DIR = os.path.join(APP_DIR, "..", "models")
DATA_DIR = os.path.join(APP_DIR, "..", "data")
EXT_VAL_DIR = os.path.join(APP_DIR, "..", "external_validation")

# Also check if models are in current dir (deployed scenario)
if not os.path.exists(MODEL_DIR):
    MODEL_DIR = os.path.join(APP_DIR, "models")
    DATA_DIR = os.path.join(APP_DIR, "data")
    EXT_VAL_DIR = os.path.join(APP_DIR, "external_validation")

sys.path.insert(0, PIPELINE_DIR)

# ASO1-2: integrity-verified pickle loading. In a full repo checkout the
# shared helper lives at pipeline/model_integrity.py; in a standalone
# deployment of this directory it may be absent, so provide an identical
# local fallback.
PIPELINE_SRC_DIR = os.path.abspath(os.path.join(APP_DIR, "..", ".."))
if os.path.isdir(PIPELINE_SRC_DIR):
    sys.path.insert(0, PIPELINE_SRC_DIR)

try:
    from model_integrity import load_pickle_verified
except ImportError:
    def load_pickle_verified(path, expected_sha256=None):
        """Fallback copy of pipeline.model_integrity.load_pickle_verified
        (keep in sync with that module)."""
        if expected_sha256:
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                    h.update(chunk)
            actual = h.hexdigest()
            if actual.lower() != expected_sha256.strip().lower():
                raise ValueError(
                    f"SHA-256 mismatch for {path}: expected {expected_sha256}, "
                    f"got {actual}. Refusing to deserialize a possibly tampered artifact."
                )
        else:
            warnings.warn(
                f"Loading {path} WITHOUT integrity verification "
                "(no expected SHA-256 configured). Set the corresponding "
                "ASO_*_SHA256 environment variable to enable verification.",
                stacklevel=2,
            )
        with open(path, "rb") as fh:
            return pickle.load(fh)


# =======================================================================
#  Feature Extraction (inline for deployment)
#  NOTE: this mirrors pipeline/phase4_features.py. The Location and
#  Linkage sections below were fixed (ASO1-3) to match the training-time
#  semantics of extract_location_features() and extract_linkage_features().
# =======================================================================

CHEMISTRY_CODES = {"M": "MOE", "C": "cEt", "d": "deoxy"}


def extract_all_features(seq, chemical_pattern, modification, linkage,
                          linkage_location, smiles, cell_line, transfection,
                          aso_concentration_nm, treatment_period_hours,
                          density_cells_per_well, primer_probe_set,
                          feature_cols):
    """Extract features for a single ASO and return aligned to feature_cols."""
    feats = {}
    seq = seq.upper().strip()
    n = len(seq)

    # -- Sequence features --------------------------------------------
    feats["seq_length"] = n
    for base in "ACGT":
        feats[f"{base}_count"] = seq.count(base)
        feats[f"{base}_fraction"] = seq.count(base) / n if n > 0 else 0

    feats["GC_fraction"] = (seq.count('G') + seq.count('C')) / n if n > 0 else 0
    feats["AT_fraction"] = (seq.count('A') + seq.count('T')) / n if n > 0 else 0
    feats["purine_fraction"] = (seq.count('A') + seq.count('G')) / n if n > 0 else 0
    feats["pyrimidine_fraction"] = (seq.count('C') + seq.count('T')) / n if n > 0 else 0

    from itertools import product as iprod
    for di in [''.join(p) for p in iprod("ACGT", repeat=2)]:
        feats[f"di_{di}"] = seq.count(di)
        feats[f"di_{di}_freq"] = seq.count(di) / max(n - 1, 1)

    for tri in [''.join(p) for p in iprod("ACGT", repeat=3)]:
        feats[f"tri_{tri}"] = seq.count(tri)

    feats["CpG_count"] = seq.count("CG")
    feats["GG_count"] = seq.count("GG")
    feats["CC_count"] = seq.count("CC")

    # Homopolymer
    max_run = 1
    current = 1
    for i in range(1, n):
        if seq[i] == seq[i-1]:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 1
    feats["max_homopolymer_run"] = max_run

    gc_stretch = 0
    current = 0
    for c in seq:
        if c in "GC":
            current += 1
            gc_stretch = max(gc_stretch, current)
        else:
            current = 0
    feats["max_gc_stretch"] = gc_stretch

    # Terminal bases
    for base in "ACGT":
        feats[f"five_prime_base_{base}"] = 1 if seq[0] == base else 0
        feats[f"three_prime_base_{base}"] = 1 if seq[-1] == base else 0

    for di in [''.join(p) for p in iprod("ACGT", repeat=2)]:
        feats[f"first_2nt_{di}"] = 1 if seq[:2] == di else 0
        feats[f"last_2nt_{di}"] = 1 if seq[-2:] == di else 0

    for tri in [''.join(p) for p in iprod("ACGT", repeat=3)]:
        feats[f"first_3nt_{tri}"] = 1 if seq[:3] == tri else 0
        feats[f"last_3nt_{tri}"] = 1 if seq[-3:] == tri else 0

    # -- Thermodynamic ------------------------------------------------
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    feats["Tm_wallace"] = 2 * at + 4 * gc
    feats["Tm_gc_adjusted"] = 64.9 + 41.0 * (gc - 16.4) / n if n > 0 else 0

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rc = ''.join(complement.get(b, 'N') for b in reversed(seq))
    max_self_comp = 0
    for start in range(n):
        for end in range(start + 3, n + 1):
            sub = seq[start:end]
            if sub in rc:
                max_self_comp = max(max_self_comp, len(sub))
    feats["self_complementarity_score"] = max_self_comp

    rc_match = sum(1 for a, b in zip(seq, rc) if a == b)
    feats["rc_internal_match_score"] = rc_match / n if n > 0 else 0

    palindrome_score = 0
    for plen in range(4, min(n + 1, 9)):
        for start in range(n - plen + 1):
            sub = seq[start:start + plen]
            sub_rc = ''.join(complement.get(b, 'N') for b in reversed(sub))
            if sub == sub_rc:
                palindrome_score += 1
    feats["palindrome_score"] = palindrome_score
    feats["terminal_gc_clamp"] = int(seq[0] in "GC" and seq[-1] in "GC")

    # -- Chemistry pattern --------------------------------------------
    cp = str(chemical_pattern).strip()
    feats["chemical_pattern_length"] = len(cp)
    for code in CHEMISTRY_CODES:
        feats[f"count_{code}"] = cp.count(code)
        feats[f"fraction_{code}"] = cp.count(code) / len(cp) if len(cp) > 0 else 0

    left_wing = 0
    for c in cp:
        if c != 'd':
            left_wing += 1
        else:
            break
    right_wing = 0
    for c in reversed(cp):
        if c != 'd':
            right_wing += 1
        else:
            break
    gap = cp.count('d')
    feats["left_modified_wing_length"] = left_wing
    feats["right_modified_wing_length"] = right_wing
    feats["central_deoxy_gap_length"] = gap
    feats["symmetry_flag"] = int(left_wing == right_wing)
    feats["n_chemistry_transitions"] = sum(1 for j in range(1, len(cp)) if cp[j] != cp[j-1])

    arch = f"{left_wing}-{gap}-{right_wing}"

    # -- Modification -------------------------------------------------
    mod_lower = str(modification).lower()
    feats["contains_MOE"] = int("moe" in mod_lower)
    feats["contains_cEt"] = int("cet" in mod_lower)
    feats["contains_5mc"] = int("5-methylcytosine" in mod_lower or "5-methylc" in mod_lower)
    feats["contains_deoxy"] = int("deoxy" in mod_lower)
    feats["contains_LNA"] = int("lna" in mod_lower)
    types = [t.strip() for t in mod_lower.split("/") if t.strip()]
    feats["n_modification_types"] = len(types)

    # -- Location (ASO1-3 fix) -----------------------------------------
    # Training semantics (phase4_features.extract_location_features) parse the
    # MODIFICATION Location field. The app has no separate Location input, so
    # the modification location is derived from the chemical pattern exactly as
    # phase9_tofersen.predict_tofersen() does for tofersen: indices of all
    # non-deoxy ('d') positions.
    mod_positions = [str(i) for i, c in enumerate(cp) if c != 'd']
    location_str = "?".join(mod_positions) + "/C/else" if mod_positions else "else"
    loc = location_str.strip()
    try:
        parts = loc.split("/")
        feats["n_location_groups"] = len(parts)
        all_pos = []
        for p in parts:
            all_pos.extend([int(x) for x in p.split("?") if x.strip().isdigit()])
        feats["n_position_indices"] = len(all_pos)
        if all_pos:
            feats["min_modified_position"] = min(all_pos)
            feats["max_modified_position"] = max(all_pos)
            # matches phase4: validity flag = all positions within sequence
            feats["mod_positions_valid"] = int(all(0 <= p < n for p in all_pos))
        else:
            # phase4 stores NaN here, later filled with 0 in build_feature_matrix
            feats["min_modified_position"] = 0
            feats["max_modified_position"] = 0
            feats["mod_positions_valid"] = 0
    except Exception:
        feats["n_location_groups"] = 0
        feats["n_position_indices"] = 0
        feats["min_modified_position"] = 0
        feats["max_modified_position"] = 0
        feats["mod_positions_valid"] = 0

    # -- Linkage (ASO1-3 fix) ------------------------------------------
    # Faithful single-row replica of phase4_features.extract_linkage_features:
    # PO/PS counts are computed from the linkage_location pattern, where listed
    # positions are PO (phosphodiester) and the remaining internucleotide
    # linkages are PS (phosphorothioate). "else"/"" means all-PS.
    link_lower = str(linkage).lower()
    feats["contains_PS"] = int("phosphorothioate" in link_lower)
    feats["contains_PO"] = int("phosphodiester" in link_lower)
    link_types = [t.strip() for t in link_lower.split("/") if t.strip()]
    feats["linkage_type_count"] = len(link_types)

    n_linkages = n - 1  # N-1 internucleotide linkages (as in phase4)
    link_loc = str(linkage_location).strip()
    try:
        if link_loc.lower() in ("else", "nan", ""):
            # All linkages are the same type (pure PS)
            feats["predicted_PS_count"] = n_linkages
            feats["predicted_PO_count"] = 0
        else:
            # Listed positions are PO; the rest are PS
            parts = link_loc.split("/")
            po_positions = []
            for p in parts:
                if p.strip().lower() == "else":
                    continue
                po_positions.extend([int(x) for x in p.split("?") if x.strip().isdigit()])
            n_po = len(po_positions)
            n_ps = n_linkages - n_po
            feats["predicted_PO_count"] = n_po
            feats["predicted_PS_count"] = max(0, n_ps)
    except Exception:
        feats["predicted_PS_count"] = n_linkages
        feats["predicted_PO_count"] = 0

    feats["predicted_PS_fraction"] = feats["predicted_PS_count"] / n_linkages if n_linkages > 0 else 0
    feats["predicted_PO_fraction"] = feats["predicted_PO_count"] / n_linkages if n_linkages > 0 else 0
    feats["linkage_positions_valid"] = 1

    # Linkage position indices / transition count (as in phase4)
    try:
        if link_loc.lower() not in ("else", "nan", ""):
            parts = link_loc.split("/")
            positions = []
            for p in parts:
                if p.strip().lower() == "else":
                    continue
                positions.extend([int(x) for x in p.split("?") if x.strip().isdigit()])
            feats["n_linkage_position_indices"] = len(positions)
            feats["linkage_transition_count"] = min(len(positions) * 2, n_linkages)
        else:
            feats["n_linkage_position_indices"] = 0
            feats["linkage_transition_count"] = 0
    except Exception:
        feats["n_linkage_position_indices"] = 0
        feats["linkage_transition_count"] = 0

    # -- SMILES -------------------------------------------------------
    smiles_str = str(smiles) if smiles else ""
    feats["smiles_length"] = len(smiles_str)
    feats["smiles_has_stereo"] = int("@" in smiles_str)
    feats["smiles_ring_count_approx"] = sum(1 for c in smiles_str if c.isdigit()) // 2

    # -- Experimental -------------------------------------------------
    feats["log10_aso_concentration_nm"] = np.log10(max(aso_concentration_nm, 1))
    feats["treatment_period_hours"] = treatment_period_hours
    feats["log10_density_cells_per_well"] = np.log10(max(density_cells_per_well, 1))

    # One-hot for cell_line, transfection, primer_probe_set
    known_cells = ["A431", "HepG2", "SH-SY5Y"]
    for cl in known_cells:
        feats[f"cell_{cl}"] = 1 if cell_line == cl else 0

    known_trans = ["electroporation", "free_uptake"]
    for tf in known_trans:
        feats[f"transfection_{tf}"] = 1 if transfection == tf else 0

    known_primers = ["HTS90", "RTS3898"]
    for pp in known_primers:
        feats[f"primer_{pp}"] = 1 if primer_probe_set == pp else 0

    # Align to expected feature columns
    result = {}
    for col in feature_cols:
        if col in feats:
            result[col] = feats[col]
        elif col.startswith("gapmer_"):
            expected_arch = col.replace("gapmer_", "")
            result[col] = 1 if arch == expected_arch else 0
        else:
            result[col] = 0

    return result


# =======================================================================
#  App UI
# =======================================================================

def load_resources():
    """Load model, pipeline, and tofersen reference (integrity-verified)."""
    try:
        # ASO1-2: verify SHA-256 before deserializing when env vars are set
        pipeline = load_pickle_verified(
            os.path.join(MODEL_DIR, "preprocessing_pipeline.pkl"),
            os.environ.get("ASO_PREPROCESSING_SHA256"),
        )

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
            model = load_pickle_verified(
                os.path.join(MODEL_DIR, "best_model.pkl"),
                os.environ.get("ASO_BEST_MODEL_SHA256"),
            )

        # Tofersen reference
        tof_path = os.path.join(EXT_VAL_DIR, "tofersen_reference.json")
        if not os.path.exists(tof_path):
            tof_path = os.path.join(APP_DIR, "tofersen_reference.json")
        with open(tof_path) as f:
            tofersen_ref = json.load(f)

        return model, pipeline, tofersen_ref
    except Exception as e:
        st.error(f"Failed to load model resources: {e}")
        return None, None, None


def main():
    # -- Custom CSS ---------------------------------------------------
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    .stApp {
        font-family: 'Inter', sans-serif;
    }
    .main-header {
        background: linear-gradient(135deg, #0f0c29 0%, #302b63 50%, #24243e 100%);
        padding: 2rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        color: white;
        text-align: center;
    }
    .main-header h1 {
        font-size: 2.5rem;
        font-weight: 700;
        margin: 0;
    }
    .main-header p {
        font-size: 1.1rem;
        opacity: 0.85;
        margin-top: 0.5rem;
    }
    .result-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 12px;
        color: white;
        text-align: center;
        margin: 0.5rem 0;
    }
    .result-card h2 {
        margin: 0;
        font-size: 2.5rem;
        font-weight: 700;
    }
    .result-card p {
        margin: 0.3rem 0 0;
        opacity: 0.9;
    }
    .benchmark-card {
        padding: 1.5rem;
        border-radius: 12px;
        text-align: center;
        margin: 0.5rem 0;
        color: white;
    }
    .higher { background: linear-gradient(135deg, #11998e, #38ef7d); }
    .comparable { background: linear-gradient(135deg, #f093fb, #f5576c); }
    .lower { background: linear-gradient(135deg, #ff6a00, #ee0979); }
    .warning-box {
        background: #fff3cd;
        border-left: 4px solid #ffc107;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
        color: #856404;
    }
    .disclaimer {
        background: #f8f9fa;
        border-left: 4px solid #6c757d;
        padding: 1rem;
        border-radius: 8px;
        margin: 1rem 0;
        font-size: 0.85rem;
        color: #495057;
    }
    </style>
    """, unsafe_allow_html=True)

    # -- Header -------------------------------------------------------
    st.markdown("""
    <div class="main-header">
        <h1>\U0001f9ec SOD1 ASO Inhibition Predictor</h1>
        <p>In-silico prediction of antisense oligonucleotide knockdown efficiency</p>
    </div>
    """, unsafe_allow_html=True)

    # Load resources
    model, pipeline, tofersen_ref = load_resources()
    if model is None:
        st.stop()

    feature_cols = pipeline["feature_cols"]
    # ASO1-6 fix (companion to phase9): no fabricated 50.0 default -- if the
    # reference has no prediction, the benchmark is shown as unavailable.
    tofersen_pred = tofersen_ref.get("predicted_inhibition_percent")
    if tofersen_pred is not None:
        try:
            tofersen_pred = float(tofersen_pred)
        except (TypeError, ValueError):
            tofersen_pred = None

    # -- Sidebar: Input -----------------------------------------------
    with st.sidebar:
        st.header("\U0001f52c ASO Input Parameters")

        sequence = st.text_input("ASO Sequence (5'->3')", value="CCGTCGCCCTTCAGCACGCA",
                                  help="DNA nucleotides only: A, C, G, T")

        chemical_pattern = st.text_input("Chemical Pattern",
                                          value="MMMMMddddddddddMMMMM",
                                          help="M=MOE, C=cEt, d=deoxy. Must match sequence length.")

        modification = st.selectbox("Modification",
                                     ["MOE/5-methylcytosines/deoxy", "MOE/cEt/5-methylcytosines/deoxy"])

        linkage = st.selectbox("Linkage",
                                ["phosphorothioate", "phosphodiester/phosphorothioate"])

        linkage_location = st.text_input("Linkage Location", value="else",
                                          help="Position indices for PO linkages, or 'else' for all PS")

        smiles = st.text_area("SMILES (optional)", value="", height=68)

        st.divider()
        st.subheader("\U0001f9ea Experimental Conditions")

        cell_line = st.selectbox("Cell Line", ["HepG2", "A431", "SH-SY5Y"])
        transfection = st.selectbox("Transfection Method", ["electroporation", "free_uptake"])
        aso_conc = st.number_input("ASO Concentration (nM)", min_value=1.0, max_value=50000.0, value=3000.0)
        treatment_hours = st.selectbox("Treatment Period (hours)", [16, 24], index=0)
        density = st.number_input("Cell Density (cells/well)", min_value=1000, max_value=100000, value=20000)
        primer_probe = st.selectbox("Primer/Probe Set", ["RTS3898", "HTS90"])

        predict_btn = st.button("\U0001f680 Predict Inhibition", type="primary", use_container_width=True)

    # -- Validation ---------------------------------------------------
    if predict_btn:
        # Input validation
        errors = []
        seq_clean = sequence.upper().strip()
        invalid_chars = set(seq_clean) - set("ACGT")
        if invalid_chars:
            errors.append(f"Invalid characters in sequence: {invalid_chars}")
        if len(seq_clean) < 10 or len(seq_clean) > 30:
            errors.append(f"Sequence length ({len(seq_clean)}) should be 10-30 nt")
        if len(chemical_pattern) != len(seq_clean):
            errors.append(f"Chemical pattern length ({len(chemical_pattern)}) != sequence length ({len(seq_clean)})")
        invalid_cp = set(chemical_pattern) - set("MCd")
        if invalid_cp:
            errors.append(f"Invalid chemistry codes: {invalid_cp}")

        if errors:
            for e in errors:
                st.error(e)
        else:
            # -- Feature extraction & prediction ----------------------
            feat_dict = extract_all_features(
                seq=seq_clean, chemical_pattern=chemical_pattern,
                modification=modification, linkage=linkage,
                linkage_location=linkage_location, smiles=smiles,
                cell_line=cell_line, transfection=transfection,
                aso_concentration_nm=aso_conc, treatment_period_hours=treatment_hours,
                density_cells_per_well=density, primer_probe_set=primer_probe,
                feature_cols=feature_cols,
            )

            X_new = np.array([[feat_dict.get(c, 0) for c in feature_cols]], dtype=np.float32)
            X_new = np.nan_to_num(X_new, nan=0, posinf=0, neginf=0)

            # Predict
            if pipeline["needs_scaling"]:
                X_input = pipeline["scaler"].transform(X_new)
            else:
                X_input = X_new

            model_key = pipeline["best_model_key"]
            if model_key == "mlp":
                import torch
                model.eval()
                with torch.no_grad():
                    y_pred = model(torch.FloatTensor(X_input)).numpy()[0]
            else:
                y_pred = model.predict(X_input)[0]

            predicted = float(np.clip(y_pred, 0, 100))

            # -- Results ----------------------------------------------
            st.markdown("---")
            col1, col2, col3 = st.columns(3)

            with col1:
                st.markdown(f"""
                <div class="result-card">
                    <h2>{predicted:.1f}%</h2>
                    <p>Predicted Inhibition</p>
                </div>
                """, unsafe_allow_html=True)

            if tofersen_pred is not None:
                with col2:
                    st.markdown(f"""
                    <div class="result-card" style="background: linear-gradient(135deg, #2193b0, #6dd5ed);">
                        <h2>{tofersen_pred:.1f}%</h2>
                        <p>Tofersen Reference (in-silico model prediction)</p>
                    </div>
                    """, unsafe_allow_html=True)

                delta = predicted - tofersen_pred
                if delta > 10:
                    cat_class = "higher"
                    cat_text = "Higher than tofersen in-silico reference"
                    cat_emoji = "\U0001f7e2"
                elif delta < -10:
                    cat_class = "lower"
                    cat_text = "Lower than tofersen in-silico reference"
                    cat_emoji = "\U0001f534"
                else:
                    cat_class = "comparable"
                    cat_text = "Comparable to tofersen in-silico reference"
                    cat_emoji = "\U0001f7e1"

                with col3:
                    st.markdown(f"""
                    <div class="benchmark-card {cat_class}">
                        <h2>{cat_emoji} {delta:+.1f}%</h2>
                        <p>{cat_text}</p>
                    </div>
                    """, unsafe_allow_html=True)
            else:
                # No fabricated benchmark: reference prediction unavailable
                delta = None
                cat_text = "No tofersen benchmark (reference prediction unavailable)"
                with col2:
                    st.markdown("""
                    <div class="result-card" style="background: linear-gradient(135deg, #6c757d, #adb5bd);">
                        <h2>N/A</h2>
                        <p>Tofersen reference unavailable</p>
                    </div>
                    """, unsafe_allow_html=True)
                with col3:
                    st.markdown("""
                    <div class="benchmark-card comparable">
                        <h2>\u2014</h2>
                        <p>No tofersen benchmark available</p>
                    </div>
                    """, unsafe_allow_html=True)

            # -- Detailed Results -------------------------------------
            st.markdown("### \U0001f4ca Detailed Results")
            results_data = {
                "Parameter": [
                    "Sequence", "Length", "Chemical Pattern",
                    "Gapmer Architecture", "GC Content",
                    "Predicted Inhibition", "Tofersen in-silico reference",
                    "Delta vs tofersen reference", "Category",
                    "Model", "Cell Line", "Concentration (nM)"
                ],
                "Value": [
                    seq_clean, len(seq_clean), chemical_pattern,
                    f"{feat_dict.get('left_modified_wing_length', '?')}-{feat_dict.get('central_deoxy_gap_length', '?')}-{feat_dict.get('right_modified_wing_length', '?')}",
                    f"{feat_dict.get('GC_fraction', 0):.1%}",
                    f"{predicted:.2f}%",
                    f"{tofersen_pred:.2f}%" if tofersen_pred is not None else "unavailable",
                    f"{delta:+.2f}%" if delta is not None else "N/A",
                    cat_text,
                    pipeline["best_model_name"], cell_line, f"{aso_conc}"
                ]
            }
            st.table(pd.DataFrame(results_data))

            # -- Disclaimer -------------------------------------------
            st.markdown("""
            <div class="disclaimer">
                <strong>\u26a0\ufe0f Disclaimer:</strong> This tool provides research-use-only in-silico predictions.
                It does not provide clinical, therapeutic, or regulatory advice.
                Experimental validation is required. The tofersen value shown is this model's own
                prediction under imputed experimental conditions (no measured ground truth):
                it is a model self-consistency reference, NOT experimental validation,
                and not a clinical efficacy claim.
            </div>
            """, unsafe_allow_html=True)

    else:
        # -- Welcome content ------------------------------------------
        st.markdown("### \U0001f448 Enter ASO parameters in the sidebar and click **Predict Inhibition**")

        st.markdown("#### About This Tool")
        st.markdown("""
        This web application predicts the mRNA inhibition percentage of novel antisense
        oligonucleotides (ASOs) targeting the human **SOD1** transcript. The model was trained
        on ~2,264 SOD1-targeting ASOs with various chemical modifications, backbone linkages,
        and experimental conditions.

        **Features:**
        - \U0001f9ec Sequence composition analysis
        - \u2697\ufe0f Chemistry pattern recognition (MOE/cEt/deoxy gapmers)
        - \U0001f517 Backbone linkage analysis
        - \U0001f4ca Experimental condition adjustment
        - \U0001f3c6 Tofersen in-silico reference comparison (model prediction, not experimental)
        """)


if __name__ == "__main__":
    main()
