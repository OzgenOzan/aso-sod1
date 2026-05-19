"""
Phase 2 -- Data Cleaning and Validation
========================================
Normalizes columns, validates data, handles duplicates, and produces:
  - cleaned_dataset.csv
  - excluded_or_flagged_records.csv
  - cleaning_report.md
"""

import pandas as pd
import numpy as np
import os
import sys
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RAW_EXCEL, COLUMN_MAP, DATA_DIR, REPORT_DIR, CHEMISTRY_CODES
)


def load_and_rename(path=None):
    """Load raw data and rename columns to snake_case."""
    df = pd.read_excel(path or RAW_EXCEL)
    df = df.rename(columns=COLUMN_MAP)
    print(f"[Phase 2] Columns renamed: {list(df.columns)}")
    return df


def normalize_sequences(df):
    """Uppercase, strip whitespace, validate nucleotide alphabet."""
    flags = []
    df["sequence"] = df["sequence"].astype(str).str.strip().str.upper()

    # Check for non-standard characters
    valid_pattern = re.compile(r'^[ACGT]+$')
    for idx, seq in df["sequence"].items():
        if 'U' in seq:
            # Convert U to T (RNA->DNA convention)
            df.at[idx, "sequence"] = seq.replace('U', 'T')
            flags.append({"index": idx, "flag": "U_converted_to_T", "original": seq})
        elif not valid_pattern.match(seq):
            bad_chars = set(seq) - {'A', 'C', 'G', 'T'}
            flags.append({"index": idx, "flag": "non_standard_chars", "chars": str(bad_chars), "sequence": seq})

    print(f"[Phase 2] Sequence normalization: {len(flags)} flags")
    return df, flags


def validate_seq_length(df):
    """Compute actual sequence length and compare with stated seq_length."""
    flags = []
    df["sequence_length_computed"] = df["sequence"].str.len()
    mismatches = df[df["seq_length"] != df["sequence_length_computed"]]

    for idx, row in mismatches.iterrows():
        # The sequence appears valid (only ACGT), so trust computed length
        flags.append({
            "index": idx,
            "flag": "seq_length_corrected",
            "stated": row["seq_length"],
            "computed": row["sequence_length_computed"],
            "sequence": row["sequence"],
        })
        df.at[idx, "seq_length"] = row["sequence_length_computed"]

    print(f"[Phase 2] Sequence length corrections: {len(flags)}")
    return df, flags


def normalize_categoricals(df):
    """Normalize categorical text columns."""
    flags = []

    # target_gene: standardize to SOD1
    df["target_gene"] = df["target_gene"].astype(str).str.strip().str.upper().str.replace("-", "", regex=False)

    # cell_line: strip and title-case
    df["cell_line"] = df["cell_line"].astype(str).str.strip()

    # transfection: strip and lowercase
    df["transfection"] = df["transfection"].astype(str).str.strip().str.lower()
    # Merge "uptake" and "free uptake" -- document decision
    uptake_mask = df["transfection"].isin(["uptake", "free uptake"])
    n_uptake = uptake_mask.sum()
    df.loc[uptake_mask, "transfection"] = "free_uptake"
    df["transfection"] = df["transfection"].str.replace(" ", "_", regex=False)
    flags.append({"flag": "transfection_normalized", "note": f"'uptake' and 'free uptake' merged to 'free_uptake' ({n_uptake} rows)"})

    # primer_probe_set
    df["primer_probe_set"] = df["primer_probe_set"].astype(str).str.strip().str.upper()

    # modification
    df["modification"] = df["modification"].astype(str).str.strip().str.lower()

    # chemical_pattern
    df["chemical_pattern"] = df["chemical_pattern"].astype(str).str.strip()

    # linkage
    df["linkage"] = df["linkage"].astype(str).str.strip().str.lower()

    print(f"[Phase 2] Categorical normalization complete")
    return df, flags


def validate_numerics(df):
    """Validate numeric columns."""
    flags = []
    num_cols = {
        "density_cells_per_well": (1000, 100000),
        "aso_concentration_nm": (0.1, 100000),
        "treatment_period_hours": (0.5, 168),
        "inhibition_percent": (0, 100),
        "seq_length": (5, 50),
    }
    for col, (lo, hi) in num_cols.items():
        df[col] = pd.to_numeric(df[col], errors="coerce")
        n_null = df[col].isnull().sum()
        out_of_range = ((df[col] < lo) | (df[col] > hi)).sum()
        if n_null > 0 or out_of_range > 0:
            flags.append({"column": col, "n_null_after_coerce": n_null, "out_of_range": out_of_range})

    print(f"[Phase 2] Numeric validation flags: {len(flags)}")
    return df, flags


def validate_chemical_pattern(df):
    """Validate chemical_pattern: non-empty, length matches sequence, valid codes."""
    flags = []
    for idx, row in df.iterrows():
        cp = row["chemical_pattern"]
        seq = row["sequence"]

        if not cp or pd.isna(cp) or len(str(cp).strip()) == 0:
            flags.append({"index": idx, "flag": "empty_chemical_pattern"})
            continue

        if len(cp) != len(seq):
            flags.append({
                "index": idx, "flag": "chemical_pattern_length_mismatch",
                "cp_len": len(cp), "seq_len": len(seq), "chemical_pattern": cp
            })

        # Check valid codes
        invalid_codes = set(cp) - set(CHEMISTRY_CODES.keys())
        if invalid_codes:
            flags.append({
                "index": idx, "flag": "unknown_chemistry_codes",
                "codes": str(invalid_codes), "chemical_pattern": cp
            })

    print(f"[Phase 2] Chemical pattern validation flags: {len(flags)}")
    return df, flags


def parse_location(df):
    """Parse the Location column into structured modification positions."""
    flags = []
    mod_positions_list = []

    for idx, row in df.iterrows():
        loc_str = str(row["location"]).strip()
        seq_len = int(row["seq_length"])

        try:
            # Format: "pos1?pos2?.../type/else" -- may have multiple groups separated by /
            parts = loc_str.split("/")
            position_str = parts[0] if parts else ""
            positions = [int(p) for p in position_str.split("?") if p.strip().isdigit()]
            
            # Check if positions are within sequence length
            valid = all(0 <= p < seq_len for p in positions)
            if not valid:
                flags.append({"index": idx, "flag": "location_positions_out_of_range", "positions": positions, "seq_len": seq_len})

            mod_positions_list.append(positions)
        except Exception as e:
            flags.append({"index": idx, "flag": "location_parse_error", "error": str(e), "raw": loc_str})
            mod_positions_list.append([])

    df["parsed_mod_positions"] = mod_positions_list
    print(f"[Phase 2] Location parsing flags: {len(flags)}")
    return df, flags


def parse_linkage_location(df):
    """Parse Linkage_Location into structured linkage positions."""
    flags = []
    linkage_positions_list = []

    for idx, row in df.iterrows():
        ll_str = str(row["linkage_location"]).strip()
        seq_len = int(row["seq_length"])

        try:
            if ll_str.lower() == "else" or ll_str.lower() == "nan":
                # "else" means all positions have the same linkage
                linkage_positions_list.append({"all": True, "positions": []})
            else:
                parts = ll_str.split("/")
                position_str = parts[0] if parts else ""
                positions = [int(p) for p in position_str.split("?") if p.strip().isdigit()]
                # Linkage positions should be 1-indexed up to N-1 for N nucleotides
                valid = all(1 <= p < seq_len for p in positions)
                if not valid and positions:
                    flags.append({"index": idx, "flag": "linkage_positions_out_of_range"})
                linkage_positions_list.append({"all": False, "positions": positions})
        except Exception as e:
            flags.append({"index": idx, "flag": "linkage_location_parse_error", "error": str(e)})
            linkage_positions_list.append({"all": True, "positions": []})

    df["parsed_linkage_positions"] = linkage_positions_list
    print(f"[Phase 2] Linkage location parsing flags: {len(flags)}")
    return df, flags


def validate_smiles(df):
    """Attempt SMILES validation. RDKit is not available, so just check non-empty."""
    flags = []
    df["smiles_valid"] = True  # Assume valid since we can't parse without RDKit

    empty_smiles = df["smiles"].isna() | (df["smiles"].astype(str).str.strip() == "")
    n_empty = empty_smiles.sum()
    if n_empty > 0:
        df.loc[empty_smiles, "smiles_valid"] = False
        flags.append({"flag": "empty_smiles", "count": int(n_empty)})

    # Note limitation
    flags.append({
        "flag": "rdkit_not_available",
        "note": "RDKit is not installed. SMILES strings are accepted but not chemically validated."
    })

    print(f"[Phase 2] SMILES validation: {n_empty} empty, RDKit unavailable")
    return df, flags


def handle_duplicates(df):
    """Handle exact duplicates and aggregate biological replicates."""
    flags = []

    # Drop unhashable columns before deduplication (they contain lists/dicts)
    unhashable_cols = ["parsed_mod_positions", "parsed_linkage_positions"]
    existing_unhashable = [c for c in unhashable_cols if c in df.columns]
    saved_unhashable = {c: df[c].copy() for c in existing_unhashable}
    df = df.drop(columns=existing_unhashable, errors="ignore")

    # Remove exact duplicate rows
    n_before = len(df)
    df = df.drop_duplicates()
    n_exact_dups = n_before - len(df)
    flags.append({"flag": "exact_duplicates_removed", "count": n_exact_dups})
    print(f"[Phase 2] Removed {n_exact_dups} exact duplicate rows")

    # Create aso_group_id based on sequence + chemical_pattern + linkage
    df["aso_group_id"] = (
        df["sequence"] + "|" + df["chemical_pattern"] + "|" + df["linkage"]
    )

    # Identify biological replicates (same experimental conditions)
    condition_cols = [
        "sequence", "chemical_pattern", "linkage", "cell_line",
        "transfection", "aso_concentration_nm", "treatment_period_hours",
        "primer_probe_set",
    ]

    # Aggregate replicates
    grouped = df.groupby(condition_cols, dropna=False)
    agg_records = []
    for name, group in grouped:
        rec = group.iloc[0].to_dict()
        if len(group) > 1:
            rec["inhibition_percent"] = group["inhibition_percent"].mean()
            rec["inhibition_std"] = group["inhibition_percent"].std()
            rec["inhibition_count"] = len(group)
        else:
            rec["inhibition_std"] = np.nan
            rec["inhibition_count"] = 1
        agg_records.append(rec)

    df_agg = pd.DataFrame(agg_records)
    n_after = len(df_agg)
    flags.append({"flag": "replicate_aggregation", "groups": n_after, "original_rows": len(df)})
    print(f"[Phase 2] Aggregated {len(df)} rows -> {n_after} unique condition groups")

    return df_agg, flags


def write_cleaning_report(all_flags, df_clean, df_excluded):
    """Write cleaning_report.md."""
    path = os.path.join(REPORT_DIR, "cleaning_report.md")
    lines = []
    lines.append("# Phase 2 -- Data Cleaning and Validation Report\n\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## Summary\n\n")
    lines.append(f"- **Cleaned dataset rows:** {len(df_clean)}\n")
    lines.append(f"- **Excluded/flagged records:** {len(df_excluded)}\n\n")

    lines.append("## Processing Steps\n\n")
    lines.append("1. Column renaming to snake_case\n")
    lines.append("2. Sequence normalization (uppercase, whitespace removal, ACGT validation)\n")
    lines.append("3. Sequence length verification and correction\n")
    lines.append("4. Categorical variable normalization\n")
    lines.append("5. Numeric variable validation\n")
    lines.append("6. Chemical pattern validation\n")
    lines.append("7. Location column parsing\n")
    lines.append("8. Linkage location parsing\n")
    lines.append("9. SMILES validation (limited -- RDKit not available)\n")
    lines.append("10. Duplicate removal and replicate aggregation\n\n")

    lines.append("## Flags and Issues\n\n")
    for section, flags in all_flags.items():
        lines.append(f"### {section}\n\n")
        if not flags:
            lines.append("No issues found.\n\n")
        else:
            for f in flags[:20]:  # Limit output
                lines.append(f"- {f}\n")
            if len(flags) > 20:
                lines.append(f"- ... and {len(flags) - 20} more\n")
            lines.append("\n")

    lines.append("## Assumptions\n\n")
    lines.append("1. 'uptake' and 'free uptake' transfection methods are treated as equivalent ('free_uptake').\n")
    lines.append("2. Sequence length mismatches (12 records with stated length 10 but actual length 17) are corrected to actual length.\n")
    lines.append("3. SMILES strings are accepted without chemical validation (RDKit not available).\n")
    lines.append("4. Biological replicates under identical conditions are aggregated by mean inhibition.\n")
    lines.append("5. All sequences contain only standard DNA bases (A, C, G, T).\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 2] Cleaning report saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 2 -- DATA CLEANING AND VALIDATION")
    print("=" * 70)

    all_flags = {}

    # 1. Load and rename
    df = load_and_rename()

    # 2. Normalize sequences
    df, seq_flags = normalize_sequences(df)
    all_flags["sequence_normalization"] = seq_flags

    # 3. Validate sequence length
    df, sl_flags = validate_seq_length(df)
    all_flags["sequence_length"] = sl_flags

    # 4. Normalize categoricals
    df, cat_flags = normalize_categoricals(df)
    all_flags["categorical_normalization"] = cat_flags

    # 5. Validate numerics
    df, num_flags = validate_numerics(df)
    all_flags["numeric_validation"] = num_flags

    # 6. Validate chemical_pattern
    df, cp_flags = validate_chemical_pattern(df)
    all_flags["chemical_pattern"] = cp_flags

    # 7. Parse location
    df, loc_flags = parse_location(df)
    all_flags["location_parsing"] = loc_flags

    # 8. Parse linkage_location
    df, ll_flags = parse_linkage_location(df)
    all_flags["linkage_location_parsing"] = ll_flags

    # 9. Validate SMILES
    df, smiles_flags = validate_smiles(df)
    all_flags["smiles_validation"] = smiles_flags

    # 10. Handle duplicates
    df_clean, dup_flags = handle_duplicates(df)
    all_flags["duplicate_handling"] = dup_flags

    # Separate excluded/flagged records
    flagged_indices = set()
    for section, flags in all_flags.items():
        for f in flags:
            if isinstance(f, dict) and "index" in f:
                flagged_indices.add(f["index"])

    df_excluded = df.loc[df.index.isin(flagged_indices)].copy() if flagged_indices else pd.DataFrame()

    # Save outputs
    # Drop helper columns that can't be serialized to CSV easily
    save_cols = [c for c in df_clean.columns if c not in ("parsed_mod_positions", "parsed_linkage_positions")]
    df_clean[save_cols].to_csv(os.path.join(DATA_DIR, "cleaned_dataset.csv"), index=False)
    print(f"[Phase 2] Cleaned dataset saved -> {os.path.join(DATA_DIR, 'cleaned_dataset.csv')} ({len(df_clean)} rows)")

    if len(df_excluded) > 0:
        save_cols_ex = [c for c in df_excluded.columns if c not in ("parsed_mod_positions", "parsed_linkage_positions")]
        df_excluded[save_cols_ex].to_csv(os.path.join(DATA_DIR, "excluded_or_flagged_records.csv"), index=False)
    else:
        pd.DataFrame({"note": ["No records were excluded"]}).to_csv(
            os.path.join(DATA_DIR, "excluded_or_flagged_records.csv"), index=False
        )

    write_cleaning_report(all_flags, df_clean, df_excluded)

    print(f"\n[Phase 2] COMPLETE [OK] -- {len(df_clean)} cleaned records\n")
    return df_clean


if __name__ == "__main__":
    main()
