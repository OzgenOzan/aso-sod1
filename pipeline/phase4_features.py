"""
Phase 4 -- Feature Extraction
==============================
Constructs a comprehensive feature matrix from available dataset columns.

Feature groups:
  A. Sequence composition
  B. Approximate thermodynamic
  C. Chemistry-pattern
  D. Modification-text
  E. Modification-location
  F. Backbone-linkage
  G. SMILES-derived (skipped -- RDKit unavailable)
  H. Experimental-condition
  I. Optional SOD1 transcript mapping

Outputs:
  - feature_matrix.csv
  - feature_dictionary.md
  - feature_extraction_report.md
"""

import pandas as pd
import numpy as np
import os
import sys
import re
from itertools import product
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DATA_DIR, REPORT_DIR, CHEMISTRY_CODES


# =======================================================================
#  A. Sequence Composition Features
# =======================================================================

def extract_sequence_features(df):
    """Extract all sequence-based features."""
    feats = pd.DataFrame(index=df.index)
    seqs = df["sequence"].values

    # Basic counts
    feats["seq_length"] = [len(s) for s in seqs]
    for base in "ACGT":
        feats[f"{base}_count"] = [s.count(base) for s in seqs]
        feats[f"{base}_fraction"] = feats[f"{base}_count"] / feats["seq_length"]

    feats["GC_fraction"] = (feats["G_count"] + feats["C_count"]) / feats["seq_length"]
    feats["AT_fraction"] = (feats["A_count"] + feats["T_count"]) / feats["seq_length"]
    feats["purine_fraction"] = (feats["A_count"] + feats["G_count"]) / feats["seq_length"]
    feats["pyrimidine_fraction"] = (feats["C_count"] + feats["T_count"]) / feats["seq_length"]

    # Dinucleotide frequencies
    dinucs = [''.join(p) for p in product("ACGT", repeat=2)]
    for di in dinucs:
        feats[f"di_{di}"] = [s.count(di) for s in seqs]
        feats[f"di_{di}_freq"] = feats[f"di_{di}"] / (feats["seq_length"] - 1).clip(lower=1)

    # Trinucleotide frequencies
    trinucs = [''.join(p) for p in product("ACGT", repeat=3)]
    for tri in trinucs:
        feats[f"tri_{tri}"] = [s.count(tri) for s in seqs]

    # Specific counts
    feats["CpG_count"] = [s.count("CG") for s in seqs]
    feats["GG_count"] = [s.count("GG") for s in seqs]
    feats["CC_count"] = [s.count("CC") for s in seqs]

    # Homopolymer runs
    def max_homopolymer(seq):
        max_run = 1
        current = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i - 1]:
                current += 1
                max_run = max(max_run, current)
            else:
                current = 1
        return max_run

    def max_gc_stretch(seq):
        max_run = 0
        current = 0
        for c in seq:
            if c in "GC":
                current += 1
                max_run = max(max_run, current)
            else:
                current = 0
        return max_run

    feats["max_homopolymer_run"] = [max_homopolymer(s) for s in seqs]
    feats["max_gc_stretch"] = [max_gc_stretch(s) for s in seqs]

    # Terminal bases
    feats["five_prime_base_A"] = [1 if s[0] == 'A' else 0 for s in seqs]
    feats["five_prime_base_C"] = [1 if s[0] == 'C' else 0 for s in seqs]
    feats["five_prime_base_G"] = [1 if s[0] == 'G' else 0 for s in seqs]
    feats["five_prime_base_T"] = [1 if s[0] == 'T' else 0 for s in seqs]
    feats["three_prime_base_A"] = [1 if s[-1] == 'A' else 0 for s in seqs]
    feats["three_prime_base_C"] = [1 if s[-1] == 'C' else 0 for s in seqs]
    feats["three_prime_base_G"] = [1 if s[-1] == 'G' else 0 for s in seqs]
    feats["three_prime_base_T"] = [1 if s[-1] == 'T' else 0 for s in seqs]

    # First/last 2nt and 3nt (one-hot)
    all_di = [''.join(p) for p in product("ACGT", repeat=2)]
    for di in all_di:
        feats[f"first_2nt_{di}"] = [1 if s[:2] == di else 0 for s in seqs]
        feats[f"last_2nt_{di}"] = [1 if s[-2:] == di else 0 for s in seqs]

    all_tri = [''.join(p) for p in product("ACGT", repeat=3)]
    for tri in all_tri:
        feats[f"first_3nt_{tri}"] = [1 if s[:3] == tri else 0 for s in seqs]
        feats[f"last_3nt_{tri}"] = [1 if s[-3:] == tri else 0 for s in seqs]

    print(f"[Phase 4A] Sequence features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  B. Approximate Thermodynamic Features
# =======================================================================

def extract_thermo_features(df):
    """Approximate Tm, self-complementarity, palindrome scores."""
    feats = pd.DataFrame(index=df.index)
    seqs = df["sequence"].values

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    for i, seq in enumerate(seqs):
        n = len(seq)
        gc = seq.count('G') + seq.count('C')
        at = seq.count('A') + seq.count('T')

        # Wallace rule Tm (short oligos): Tm = 2*(A+T) + 4*(G+C)
        feats.at[df.index[i], "Tm_wallace"] = 2 * at + 4 * gc

        # GC-adjusted Tm (Marmur-Doty for longer oligos):
        # Tm = 64.9 + 41*(gc-16.4)/n
        feats.at[df.index[i], "Tm_gc_adjusted"] = 64.9 + 41.0 * (gc - 16.4) / n if n > 0 else 0

        # Reverse complement
        rc = ''.join(complement.get(b, 'N') for b in reversed(seq))

        # Self-complementarity: max contiguous match of seq with its reverse complement
        max_self_comp = 0
        for start in range(n):
            for end in range(start + 3, n + 1):
                subseq = seq[start:end]
                if subseq in rc:
                    max_self_comp = max(max_self_comp, len(subseq))
        feats.at[df.index[i], "self_complementarity_score"] = max_self_comp

        # Reverse-complement internal matching score
        # Count matching positions when seq is aligned with its RC
        rc_match = sum(1 for a, b in zip(seq, rc) if a == b)
        feats.at[df.index[i], "rc_internal_match_score"] = rc_match / n if n > 0 else 0

        # Palindromic motif check (is any 4+ mer in the sequence a palindrome?)
        palindrome_score = 0
        for plen in range(4, min(n + 1, 9)):
            for start in range(n - plen + 1):
                sub = seq[start:start + plen]
                sub_rc = ''.join(complement.get(b, 'N') for b in reversed(sub))
                if sub == sub_rc:
                    palindrome_score += 1
        feats.at[df.index[i], "palindrome_score"] = palindrome_score

        # Terminal GC clamp (GC at both 5' and 3' ends)
        gc_5prime = seq[0] in "GC"
        gc_3prime = seq[-1] in "GC"
        feats.at[df.index[i], "terminal_gc_clamp"] = int(gc_5prime and gc_3prime)

    print(f"[Phase 4B] Thermodynamic features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  C. Chemistry-Pattern Features
# =======================================================================

def extract_chemistry_features(df):
    """Parse chemical_pattern and extract gapmer architecture features."""
    feats = pd.DataFrame(index=df.index)
    cps = df["chemical_pattern"].values

    for i, cp in enumerate(cps):
        idx = df.index[i]
        cp = str(cp).strip()

        feats.at[idx, "chemical_pattern_length"] = len(cp)

        # Count codes
        for code in CHEMISTRY_CODES:
            feats.at[idx, f"count_{code}"] = cp.count(code)
            feats.at[idx, f"fraction_{code}"] = cp.count(code) / len(cp) if len(cp) > 0 else 0

        # Parse gapmer architecture: left_wing (modified) - gap (deoxy) - right_wing (modified)
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

        gap_length = cp.count('d')

        feats.at[idx, "left_modified_wing_length"] = left_wing
        feats.at[idx, "right_modified_wing_length"] = right_wing
        feats.at[idx, "central_deoxy_gap_length"] = gap_length
        feats.at[idx, "gapmer_architecture"] = f"{left_wing}-{gap_length}-{right_wing}"
        feats.at[idx, "symmetry_flag"] = int(left_wing == right_wing)

        # Chemistry transitions (number of changes between adjacent positions)
        transitions = sum(1 for j in range(1, len(cp)) if cp[j] != cp[j - 1])
        feats.at[idx, "n_chemistry_transitions"] = transitions

    # One-hot encode gapmer_architecture
    archs = feats["gapmer_architecture"].unique()
    for arch in archs:
        feats[f"gapmer_{arch}"] = (feats["gapmer_architecture"] == arch).astype(int)
    feats = feats.drop(columns=["gapmer_architecture"])

    print(f"[Phase 4C] Chemistry features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  D. Modification-Text Features
# =======================================================================

def extract_modification_features(df):
    """Extract binary flags from the Modification text column."""
    feats = pd.DataFrame(index=df.index)
    mods = df["modification"].values

    for i, mod_str in enumerate(mods):
        idx = df.index[i]
        mod_lower = str(mod_str).lower()

        feats.at[idx, "contains_MOE"] = int("moe" in mod_lower)
        feats.at[idx, "contains_cEt"] = int("cet" in mod_lower)
        feats.at[idx, "contains_5mc"] = int("5-methylcytosine" in mod_lower or "5-methylc" in mod_lower)
        feats.at[idx, "contains_deoxy"] = int("deoxy" in mod_lower)
        feats.at[idx, "contains_LNA"] = int("lna" in mod_lower)

        # Count modification types
        types = [t.strip() for t in mod_lower.split("/") if t.strip()]
        feats.at[idx, "n_modification_types"] = len(types)

    print(f"[Phase 4D] Modification features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  E. Modification-Location Features
# =======================================================================

def extract_location_features(df):
    """Extract features from the Location column."""
    feats = pd.DataFrame(index=df.index)
    locs = df["location"].values
    seq_lens = df["seq_length"].values

    for i, loc_str in enumerate(locs):
        idx = df.index[i]
        loc = str(loc_str).strip()
        sl = int(seq_lens[i])

        try:
            parts = loc.split("/")
            n_groups = len(parts)
            feats.at[idx, "n_location_groups"] = n_groups

            # Extract all numeric positions
            all_positions = []
            for p in parts:
                positions = [int(x) for x in p.split("?") if x.strip().isdigit()]
                all_positions.extend(positions)

            feats.at[idx, "n_position_indices"] = len(all_positions)

            if all_positions:
                feats.at[idx, "min_modified_position"] = min(all_positions)
                feats.at[idx, "max_modified_position"] = max(all_positions)
                feats.at[idx, "mod_positions_valid"] = int(all(0 <= p < sl for p in all_positions))
            else:
                feats.at[idx, "min_modified_position"] = np.nan
                feats.at[idx, "max_modified_position"] = np.nan
                feats.at[idx, "mod_positions_valid"] = 0

        except Exception:
            feats.at[idx, "n_location_groups"] = 0
            feats.at[idx, "n_position_indices"] = 0
            feats.at[idx, "min_modified_position"] = np.nan
            feats.at[idx, "max_modified_position"] = np.nan
            feats.at[idx, "mod_positions_valid"] = 0

    print(f"[Phase 4E] Location features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  F. Backbone-Linkage Features
# =======================================================================

def extract_linkage_features(df):
    """Extract features from Linkage and Linkage_Location columns."""
    feats = pd.DataFrame(index=df.index)

    for i in range(len(df)):
        idx = df.index[i]
        linkage = str(df.iloc[i]["linkage"]).lower()
        link_loc = str(df.iloc[i]["linkage_location"]).strip()
        sl = int(df.iloc[i]["seq_length"])
        n_linkages = sl - 1  # N-1 internucleotide linkages

        feats.at[idx, "contains_PS"] = int("phosphorothioate" in linkage)
        feats.at[idx, "contains_PO"] = int("phosphodiester" in linkage)

        # Count unique linkage types
        types = [t.strip() for t in linkage.split("/") if t.strip()]
        feats.at[idx, "linkage_type_count"] = len(types)

        # Parse linkage positions
        try:
            if link_loc.lower() in ("else", "nan", ""):
                # All linkages are the same type (pure PS)
                if "phosphorothioate" in linkage and "phosphodiester" not in linkage:
                    feats.at[idx, "predicted_PS_count"] = n_linkages
                    feats.at[idx, "predicted_PO_count"] = 0
                else:
                    feats.at[idx, "predicted_PS_count"] = n_linkages
                    feats.at[idx, "predicted_PO_count"] = 0
            else:
                # Parse positions for PO (listed positions), rest are PS
                parts = link_loc.split("/")
                po_positions = []
                for p in parts:
                    if p.strip().lower() == "else":
                        continue
                    positions = [int(x) for x in p.split("?") if x.strip().isdigit()]
                    po_positions.extend(positions)

                n_po = len(po_positions)
                n_ps = n_linkages - n_po
                feats.at[idx, "predicted_PO_count"] = n_po
                feats.at[idx, "predicted_PS_count"] = max(0, n_ps)

        except Exception:
            feats.at[idx, "predicted_PS_count"] = n_linkages
            feats.at[idx, "predicted_PO_count"] = 0

        feats.at[idx, "predicted_PS_fraction"] = feats.at[idx, "predicted_PS_count"] / n_linkages if n_linkages > 0 else 0
        feats.at[idx, "predicted_PO_fraction"] = feats.at[idx, "predicted_PO_count"] / n_linkages if n_linkages > 0 else 0

        # Linkage validity
        feats.at[idx, "linkage_positions_valid"] = 1

        # Linkage transition count
        try:
            if link_loc.lower() not in ("else", "nan", ""):
                parts = link_loc.split("/")
                positions = []
                for p in parts:
                    if p.strip().lower() == "else":
                        continue
                    positions.extend([int(x) for x in p.split("?") if x.strip().isdigit()])
                feats.at[idx, "n_linkage_position_indices"] = len(positions)
                # Count transitions in linkage type across the backbone
                feats.at[idx, "linkage_transition_count"] = min(len(positions) * 2, n_linkages)
            else:
                feats.at[idx, "n_linkage_position_indices"] = 0
                feats.at[idx, "linkage_transition_count"] = 0
        except Exception:
            feats.at[idx, "n_linkage_position_indices"] = 0
            feats.at[idx, "linkage_transition_count"] = 0

    print(f"[Phase 4F] Linkage features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  G. SMILES-derived features (SKIPPED -- RDKit unavailable)
# =======================================================================

def extract_smiles_features(df):
    """SMILES features -- skipped because RDKit is not available."""
    feats = pd.DataFrame(index=df.index)
    feats["smiles_length"] = df["smiles"].astype(str).str.len()
    feats["smiles_has_stereo"] = df["smiles"].astype(str).str.contains(r"@@|@", regex=True).astype(int)
    feats["smiles_ring_count_approx"] = df["smiles"].astype(str).apply(
        lambda s: sum(1 for c in s if c.isdigit()) // 2
    )
    print(f"[Phase 4G] SMILES features: {len(feats.columns)} columns (limited -- RDKit unavailable)")
    return feats


# =======================================================================
#  H. Experimental-Condition Features
# =======================================================================

def extract_experimental_features(df):
    """One-hot encode experimental conditions and log-transform numerics."""
    feats = pd.DataFrame(index=df.index)

    # Log-transform numeric conditions
    feats["log10_aso_concentration_nm"] = np.log10(df["aso_concentration_nm"].clip(lower=1))
    feats["treatment_period_hours"] = df["treatment_period_hours"]
    feats["log10_density_cells_per_well"] = np.log10(df["density_cells_per_well"].clip(lower=1))

    # One-hot: cell_line
    for cl in sorted(df["cell_line"].unique()):
        feats[f"cell_{cl}"] = (df["cell_line"] == cl).astype(int)

    # One-hot: transfection
    for tf in sorted(df["transfection"].unique()):
        feats[f"transfection_{tf}"] = (df["transfection"] == tf).astype(int)

    # One-hot: primer_probe_set
    for pp in sorted(df["primer_probe_set"].unique()):
        feats[f"primer_{pp}"] = (df["primer_probe_set"] == pp).astype(int)

    print(f"[Phase 4H] Experimental features: {len(feats.columns)} columns")
    return feats


# =======================================================================
#  Main assembly
# =======================================================================

def build_feature_matrix(df, skip_constant_removal=False):
    """Assemble all feature groups into a single feature matrix.
    
    Args:
        df: DataFrame with required columns
        skip_constant_removal: If True, keep all columns even if constant.
            Set True when extracting features for single-sample prediction.
    """
    print("\n[Phase 4] Building feature matrix...")

    feat_groups = {}
    feat_groups["A_sequence"] = extract_sequence_features(df)
    feat_groups["B_thermo"] = extract_thermo_features(df)
    feat_groups["C_chemistry"] = extract_chemistry_features(df)
    feat_groups["D_modification"] = extract_modification_features(df)
    feat_groups["E_location"] = extract_location_features(df)
    feat_groups["F_linkage"] = extract_linkage_features(df)
    feat_groups["G_smiles"] = extract_smiles_features(df)
    feat_groups["H_experimental"] = extract_experimental_features(df)

    # Concatenate
    feature_matrix = pd.concat(list(feat_groups.values()), axis=1)

    # Add identifiers and target
    feature_matrix.insert(0, "sequence", df["sequence"].values)
    feature_matrix.insert(1, "aso_group_id", df["aso_group_id"].values)
    if "sequence_cluster_id" in df.columns:
        feature_matrix.insert(2, "sequence_cluster_id", df["sequence_cluster_id"].values)
    feature_matrix["inhibition_percent"] = df["inhibition_percent"].values

    # Remove any constant columns (zero variance) -- skip for single-sample prediction
    if not skip_constant_removal:
        numeric_cols = feature_matrix.select_dtypes(include=[np.number]).columns
        constant_cols = [c for c in numeric_cols if feature_matrix[c].nunique() <= 1 and c != "inhibition_percent"]
        if constant_cols:
            print(f"[Phase 4] Removing {len(constant_cols)} constant-value columns")
            feature_matrix = feature_matrix.drop(columns=constant_cols)

    # Handle NaN
    n_nan = feature_matrix.select_dtypes(include=[np.number]).isna().sum().sum()
    if n_nan > 0:
        print(f"[Phase 4] Filling {n_nan} NaN values with 0")
        feature_matrix = feature_matrix.fillna(0)

    return feature_matrix, feat_groups


def write_feature_reports(feature_matrix, feat_groups):
    """Write feature_dictionary.md and feature_extraction_report.md."""

    # Feature dictionary
    dict_path = os.path.join(REPORT_DIR, "feature_dictionary.md")
    lines = ["# Feature Dictionary\n\n"]
    lines.append(f"**Total features:** {len(feature_matrix.columns) - 3}  ")
    lines.append(f"(excluding sequence, aso_group_id, sequence_cluster_id, and inhibition_percent)\n\n")

    group_names = {
        "A_sequence": "A. Sequence Composition",
        "B_thermo": "B. Approximate Thermodynamic",
        "C_chemistry": "C. Chemistry-Pattern",
        "D_modification": "D. Modification-Text",
        "E_location": "E. Modification-Location",
        "F_linkage": "F. Backbone-Linkage",
        "G_smiles": "G. SMILES-Derived (Limited)",
        "H_experimental": "H. Experimental-Condition",
    }

    for key, name in group_names.items():
        cols = list(feat_groups[key].columns)
        lines.append(f"## {name} ({len(cols)} features)\n\n")
        for c in cols[:50]:  # Limit for readability
            lines.append(f"- `{c}`\n")
        if len(cols) > 50:
            lines.append(f"- ... and {len(cols) - 50} more\n")
        lines.append("\n")

    with open(dict_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    # Extraction report
    report_path = os.path.join(REPORT_DIR, "feature_extraction_report.md")
    lines = ["# Phase 4 -- Feature Extraction Report\n\n"]
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Summary\n\n")
    lines.append(f"- **Total samples:** {len(feature_matrix)}\n")
    lines.append(f"- **Total features:** {len(feature_matrix.columns)}\n")
    lines.append(f"- **Feature groups:** {len(feat_groups)}\n\n")

    lines.append("| Group | Features |\n|---|---|\n")
    for key, name in group_names.items():
        lines.append(f"| {name} | {len(feat_groups[key].columns)} |\n")

    lines.append("\n## Skipped Features\n\n")
    lines.append("### SMILES Molecular Descriptors\n")
    lines.append("- **Reason:** RDKit is not installed in the current environment.\n")
    lines.append("- **Impact:** Molecular weight, TPSA, rotatable bonds, and ring count from SMILES are unavailable.\n")
    lines.append("- **Mitigation:** Basic SMILES string-level features (length, stereo count, approximate ring count) were extracted.\n\n")

    lines.append("### SOD1 Transcript Mapping\n")
    lines.append("- **Reason:** Dataset does not contain transcript coordinates. External mapping was not performed.\n")
    lines.append("- **Impact:** Target position, UTR/CDS/intron annotation, and RNA accessibility features are unavailable.\n")
    lines.append("- **Mitigation:** None -- these features are documented as absent.\n\n")

    lines.append("## Thermodynamic Approximation Formulas\n\n")
    lines.append("- **Wallace Tm:** Tm = 2x(A+T) + 4x(G+C) -- suitable for short oligos <20nt\n")
    lines.append("- **GC-adjusted Tm (Marmur-Doty):** Tm = 64.9 + 41x(GC_count - 16.4) / length\n")
    lines.append("- Note: These are rough approximations. Nearest-neighbor Tm would be more accurate but requires specialized libraries.\n")

    with open(report_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    print(f"[Phase 4] Reports saved -> {dict_path}, {report_path}")


def main():
    print("=" * 70)
    print("PHASE 4 -- FEATURE EXTRACTION")
    print("=" * 70)

    # Load clustered dataset
    path = os.path.join(DATA_DIR, "cleaned_dataset_clustered.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(os.path.join(DATA_DIR, "cleaned_dataset.csv"))
        print("[Phase 4] WARNING: Using non-clustered dataset")

    print(f"[Phase 4] Loaded {len(df)} records")

    feature_matrix, feat_groups = build_feature_matrix(df)

    # Save
    fm_path = os.path.join(DATA_DIR, "feature_matrix.csv")
    feature_matrix.to_csv(fm_path, index=False)
    print(f"[Phase 4] Feature matrix saved -> {fm_path} ({feature_matrix.shape})")

    write_feature_reports(feature_matrix, feat_groups)

    print(f"\n[Phase 4] COMPLETE [OK] -- {feature_matrix.shape[0]} samples x {feature_matrix.shape[1]} features\n")
    return feature_matrix


if __name__ == "__main__":
    main()
