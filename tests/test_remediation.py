"""
Static remediation regression tests (ASO1).

These tests are intentionally STATIC (source-text assertions). They do NOT
import torch/sklearn/streamlit and do NOT deserialize any pickle artifacts,
so they can run in minimal CI without the heavy dependency stack.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


def test_remediation_phase7_no_test_set_model_selection():
    """ASO1-1/ASO1-4: phase7 must select by cluster-aware CV on the training
    split, not by test-set R2."""
    src = read("pipeline/phase7_modeling.py")

    # The old leakage pattern: selecting the best row by test-set R2
    assert 'valid_results["R2"].idxmax()' not in src
    assert "idxmax" not in src, "phase7 must not rank models with idxmax on result tables"

    # CV-based selection entry point must exist and be called in main()
    assert "def select_model_by_cv(" in src
    assert "select_model_by_cv(" in src.split("def main(")[1]

    # The old dead cross_validate must now actually be used for selection
    assert "cross_validate(" in src.split("def select_model_by_cv(")[1]

    # The single test-set evaluation must be clearly labeled
    assert "final held-out evaluation" in src.lower()


def test_remediation_phase9_no_fabricated_tofersen_fallback():
    """ASO1-6 (phase9): the fabricated 50.0 fallback must be gone and failure
    must raise instead of writing tofersen_reference.json."""
    src = read("pipeline/phase9_tofersen.py")

    assert "tofersen_pred = 50.0" not in src
    assert "Fallback" not in src or "refusing to write fallback" in src
    assert "refusing to write fallback benchmark" in src
    assert "raise RuntimeError" in src


def test_remediation_app_feature_extraction_matches_training():
    """ASO1-3: web tool must parse the modification Location (derived from the
    chemical pattern) for position features, and compute PO/PS counts from the
    linkage_location pattern like phase4_features.extract_linkage_features."""
    src = read("pipeline/outputs/web_tool/app.py")

    # Old bug: Location features were parsed from linkage_location
    assert not re.search(r"(?<![\w])loc = str\(linkage_location\)", src)

    # New behavior: modification location derived from the chemical pattern
    assert "location_str" in src
    assert "mod_positions" in src

    # Old bug: PO/PS hardcoded
    assert 'feats["predicted_PS_count"] = n_linkages\n    feats["predicted_PO_count"] = 0\n    feats["predicted_PS_fraction"] = 1.0' not in src
    # New behavior: PO count computed from parsed PO positions
    assert "po_positions" in src
    assert 'feats["predicted_PO_count"] = n_po' in src


def test_remediation_torch_pin_cve_fix():
    """ASO1-7: torch must be pinned at >= 2.6.0 (CVE-2025-32434 fix)."""
    src = read("requirements.txt")
    m = re.search(r"^torch==(\d+)\.(\d+)\.(\d+)", src, flags=re.MULTILINE)
    assert m, "torch pin not found in requirements.txt"
    major, minor, patch = (int(x) for x in m.groups())
    assert (major, minor, patch) >= (2, 6, 0), f"torch pin {m.group(0)} is below 2.6.0"


def test_remediation_pickle_loads_routed_through_integrity_helper():
    """ASO1-2: pickle.load call sites must go through model_integrity."""
    mi = read("pipeline/model_integrity.py")
    assert "def sha256_of(" in mi
    assert "def load_pickle_verified(" in mi

    for rel in [
        "pipeline/phase8_validation.py",
        "pipeline/phase9_tofersen.py",
        "verify_pipeline.py",
        "pipeline/outputs/web_tool/app.py",
    ]:
        src = read(rel)
        assert "load_pickle_verified" in src, f"{rel} does not use load_pickle_verified"
        # No bare pickle.load( call remaining at module level
        assert not re.search(r"(?<!\w)pickle\.load\(", src.replace("return pickle.load(fh)", "").replace("return pickle.load(f)", "")), \
            f"{rel} still contains a bare pickle.load call"
