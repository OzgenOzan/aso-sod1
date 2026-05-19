"""
Phase 1 -- Data Loading and Initial Audit
=========================================
Loads the raw dataset, performs an initial audit, and produces:
  - raw_data_audit_report.md
  - column_dictionary.md
"""

import pandas as pd
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RAW_EXCEL, COLUMN_MAP, REPORT_DIR, DATA_DIR,
    CHEMISTRY_CODES, INHIBITION_BINS, INHIBITION_BIN_LABELS
)


def load_raw_data():
    """Load the raw Excel dataset."""
    df = pd.read_excel(RAW_EXCEL)
    print(f"[Phase 1] Loaded {len(df)} rows x {len(df.columns)} columns from {RAW_EXCEL}")
    return df


def audit_dataset(df):
    """Run comprehensive audit checks and return audit dict."""
    audit = {}

    # Basic counts
    audit["n_rows"] = len(df)
    audit["n_columns"] = len(df.columns)
    audit["columns"] = list(df.columns)

    # Unique counts
    audit["n_unique_sequences"] = df["Sequence"].nunique()
    audit["n_unique_seq_chem"] = df.groupby(["Sequence", "Chemical_Pattern"]).ngroups
    audit["n_unique_cell_lines"] = df["Cell_line"].nunique()
    audit["cell_lines"] = sorted(df["Cell_line"].unique().tolist())
    audit["n_unique_transfection"] = df["Transfection"].nunique()
    audit["transfection_methods"] = sorted(df["Transfection"].unique().tolist())
    audit["n_unique_chemical_patterns"] = df["Chemical_Pattern"].nunique()
    audit["chemical_patterns"] = sorted(df["Chemical_Pattern"].unique().tolist())
    audit["n_unique_linkage"] = df["Linkage"].nunique()
    audit["linkage_types"] = sorted(df["Linkage"].unique().tolist())

    # Missingness
    audit["missingness"] = df.isnull().sum().to_dict()

    # Inhibition summary statistics
    inh = df["Inhibition(%)"]
    audit["inhibition_stats"] = {
        "count": int(inh.count()),
        "mean": round(float(inh.mean()), 2),
        "std": round(float(inh.std()), 2),
        "min": float(inh.min()),
        "25%": float(inh.quantile(0.25)),
        "50%": float(inh.quantile(0.50)),
        "75%": float(inh.quantile(0.75)),
        "max": float(inh.max()),
    }

    # Sequence length distribution
    sl = df["seq_length"]
    audit["seq_length_stats"] = {
        "min": int(sl.min()),
        "max": int(sl.max()),
        "mean": round(float(sl.mean()), 2),
        "distribution": sl.value_counts().sort_index().to_dict(),
    }

    # Target gene check
    target_genes = df["Target_gene"].unique().tolist()
    audit["target_genes"] = target_genes
    audit["all_sod1"] = all(g.upper().replace("-", "") == "SOD1" for g in target_genes)

    # Inhibition range check
    audit["inhibition_all_numeric"] = pd.api.types.is_numeric_dtype(df["Inhibition(%)"])
    audit["inhibition_out_of_range"] = int(
        ((df["Inhibition(%)"] < 0) | (df["Inhibition(%)"] > 100)).sum()
    )

    # Sequence length validation
    df_temp = df.copy()
    df_temp["computed_len"] = df_temp["Sequence"].str.len()
    mismatches = df_temp[df_temp["seq_length"] != df_temp["computed_len"]]
    audit["seq_length_mismatches"] = len(mismatches)
    if len(mismatches) > 0:
        audit["seq_length_mismatch_details"] = mismatches[
            ["ISIS", "Sequence", "seq_length", "computed_len"]
        ].to_dict("records")

    # Duplicates
    audit["n_fully_duplicated_rows"] = int(df.duplicated().sum())
    audit["n_duplicated_sequences"] = int(df["Sequence"].duplicated().sum())
    audit["n_dup_seq_chem"] = int(
        df.duplicated(subset=["Sequence", "Chemical_Pattern"]).sum()
    )
    audit["n_dup_seq_chem_link"] = int(
        df.duplicated(subset=["Sequence", "Chemical_Pattern", "Linkage"]).sum()
    )

    # Conflicting inhibition under identical conditions
    condition_cols = [
        "Sequence", "Chemical_Pattern", "Linkage", "Cell_line",
        "Transfection", "ASO_volume(nM)", "Treatment_Period(hours)",
        "Primer_probe_set",
    ]
    grouped = df.groupby(condition_cols)["Inhibition(%)"]
    conflict_groups = grouped.filter(lambda x: x.nunique() > 1 if len(x) > 1 else False)
    audit["n_conflicting_inhibition_groups"] = int(
        grouped.apply(lambda x: x.nunique() > 1 if len(x) > 1 else False).sum()
    )
    audit["n_conflicting_inhibition_rows"] = len(conflict_groups)

    return audit


def write_audit_report(audit):
    """Write raw_data_audit_report.md."""
    path = os.path.join(REPORT_DIR, "raw_data_audit_report.md")
    lines = []
    lines.append("# Phase 1 -- Raw Data Audit Report\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().isoformat()}\n\n")

    lines.append("## 1. Dataset Overview\n")
    lines.append(f"- **Rows:** {audit['n_rows']}\n")
    lines.append(f"- **Columns:** {audit['n_columns']}\n")
    lines.append(f"- **Unique ASO sequences:** {audit['n_unique_sequences']}\n")
    lines.append(f"- **Unique sequence + chemical pattern combinations:** {audit['n_unique_seq_chem']}\n")
    lines.append(f"- **Unique cell lines:** {audit['n_unique_cell_lines']} -- {audit['cell_lines']}\n")
    lines.append(f"- **Unique transfection methods:** {audit['n_unique_transfection']} -- {audit['transfection_methods']}\n")
    lines.append(f"- **Unique chemical patterns:** {audit['n_unique_chemical_patterns']}\n")
    lines.append(f"- **Unique linkage patterns:** {audit['n_unique_linkage']} -- {audit['linkage_types']}\n")

    lines.append("\n## 2. Missingness per Column\n")
    lines.append("| Column | Missing |\n|---|---|\n")
    for col, val in audit["missingness"].items():
        lines.append(f"| {col} | {val} |\n")

    lines.append("\n## 3. Inhibition (%) Summary Statistics\n")
    lines.append("| Statistic | Value |\n|---|---|\n")
    for k, v in audit["inhibition_stats"].items():
        lines.append(f"| {k} | {v} |\n")

    lines.append("\n## 4. Sequence Length Distribution\n")
    lines.append(f"- **Range:** {audit['seq_length_stats']['min']} - {audit['seq_length_stats']['max']}\n")
    lines.append(f"- **Mean:** {audit['seq_length_stats']['mean']}\n")
    lines.append("| Length | Count |\n|---|---|\n")
    for length, count in audit["seq_length_stats"]["distribution"].items():
        lines.append(f"| {length} | {count} |\n")

    lines.append("\n## 5. Target Gene Validation\n")
    lines.append(f"- **Target genes found:** {audit['target_genes']}\n")
    lines.append(f"- **All target SOD1:** {'[OK] Yes' if audit['all_sod1'] else '[X] No -- requires investigation'}\n")

    lines.append("\n## 6. Inhibition Range Validation\n")
    lines.append(f"- **Numeric type:** {'[OK] Yes' if audit['inhibition_all_numeric'] else '[X] No'}\n")
    lines.append(f"- **Values outside [0, 100]:** {audit['inhibition_out_of_range']}\n")

    lines.append("\n## 7. Sequence Length Verification\n")
    lines.append(f"- **Mismatches (seq_length != len(Sequence)):** {audit['seq_length_mismatches']}\n")
    if audit["seq_length_mismatches"] > 0:
        lines.append("\n| ISIS | Sequence | Stated Length | Computed Length |\n|---|---|---|---|\n")
        for d in audit.get("seq_length_mismatch_details", []):
            lines.append(f"| {d['ISIS']} | {d['Sequence']} | {d['seq_length']} | {d['computed_len']} |\n")

    lines.append("\n## 8. Duplicate Analysis\n")
    lines.append(f"- **Fully duplicated rows:** {audit['n_fully_duplicated_rows']}\n")
    lines.append(f"- **Rows with duplicated sequences (same Sequence, different conditions):** {audit['n_duplicated_sequences']}\n")
    lines.append(f"- **Rows with duplicated Sequence + Chemical_Pattern:** {audit['n_dup_seq_chem']}\n")
    lines.append(f"- **Rows with duplicated Sequence + Chemical_Pattern + Linkage:** {audit['n_dup_seq_chem_link']}\n")

    lines.append("\n## 9. Conflicting Inhibition Values\n")
    lines.append(f"- **Groups with conflicting inhibition under identical conditions:** {audit['n_conflicting_inhibition_groups']}\n")
    lines.append(f"- **Rows involved:** {audit['n_conflicting_inhibition_rows']}\n")

    lines.append("\n## 10. Chemical Patterns Observed\n")
    lines.append("| Pattern | Length | Wing-Gap-Wing |\n|---|---|---|\n")
    for cp in audit["chemical_patterns"]:
        # Parse gapmer architecture
        left = 0
        for c in cp:
            if c != 'd':
                left += 1
            else:
                break
        gap = cp.count('d')
        right = len(cp) - left - gap
        lines.append(f"| `{cp}` | {len(cp)} | {left}-{gap}-{right} |\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 1] Audit report saved -> {path}")


def write_column_dictionary():
    """Write column_dictionary.md."""
    path = os.path.join(REPORT_DIR, "column_dictionary.md")
    entries = [
        ("ISIS (-> no)", "Integer", "Unique ISIS identifier for the ASO compound"),
        ("Target_gene (-> target_gene)", "Categorical", "Target gene symbol (SOD-1 in this dataset)"),
        ("Cell_line (-> cell_line)", "Categorical", "Cell line used in the assay (HepG2, A431, SH-SY5Y)"),
        ("Density(cells/well) (-> density_cells_per_well)", "Numeric", "Cell seeding density per well"),
        ("Transfection (-> transfection)", "Categorical", "Transfection or delivery method (electroporation, uptake, free uptake)"),
        ("ASO_volume(nM) (-> aso_concentration_nm)", "Numeric", "ASO concentration in nanomolar"),
        ("Treatment_Period(hours) (-> treatment_period_hours)", "Numeric", "Duration of ASO treatment in hours"),
        ("Primer_probe_set (-> primer_probe_set)", "Categorical", "qPCR primer/probe set used for mRNA quantification"),
        ("Sequence (-> sequence)", "Text/Sequence", "ASO nucleotide sequence (5'->3', DNA alphabet A/C/G/T)"),
        ("Modification (-> modification)", "Categorical", "High-level chemical modification description"),
        ("Location (-> location)", "Text/Structured", "Position indices for modifications, format: pos1?pos2?.../type/else"),
        ("Chemical_Pattern (-> chemical_pattern)", "Text/Categorical", "Per-position chemistry code string: M=MOE, C=cEt, d=deoxy"),
        ("Linkage (-> linkage)", "Categorical", "Backbone linkage type(s): phosphorothioate and/or phosphodiester"),
        ("Linkage_Location (-> linkage_location)", "Text/Structured", "Position indices for linkage types"),
        ("Smiles (-> smiles)", "Text/Chemical", "SMILES representation of the full ASO molecule"),
        ("Inhibition(%) (-> inhibition_percent)", "Numeric (0-100)", "Target variable: mRNA inhibition percentage"),
        ("seq_length (-> seq_length)", "Numeric", "Stated sequence length in nucleotides"),
    ]

    lines = []
    lines.append("# Column Dictionary\n\n")
    lines.append("| Original Column (-> Standardized) | Type | Description |\n")
    lines.append("|---|---|---|\n")
    for name, dtype, desc in entries:
        lines.append(f"| {name} | {dtype} | {desc} |\n")

    lines.append("\n## Chemistry Code Definitions\n\n")
    lines.append("| Code | Meaning |\n|---|---|\n")
    for code, meaning in CHEMISTRY_CODES.items():
        lines.append(f"| `{code}` | {meaning} |\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[Phase 1] Column dictionary saved -> {path}")


def main():
    print("=" * 70)
    print("PHASE 1 -- DATA LOADING AND INITIAL AUDIT")
    print("=" * 70)

    df = load_raw_data()
    audit = audit_dataset(df)
    write_audit_report(audit)
    write_column_dictionary()

    print(f"\n[Phase 1] Summary:")
    print(f"  Rows: {audit['n_rows']}")
    print(f"  Unique sequences: {audit['n_unique_sequences']}")
    print(f"  Unique seq+chem combos: {audit['n_unique_seq_chem']}")
    print(f"  Missing values: {sum(audit['missingness'].values())} total")
    print(f"  Seq length mismatches: {audit['seq_length_mismatches']}")
    print(f"  Fully duplicated rows: {audit['n_fully_duplicated_rows']}")
    print(f"  Inhibition range: {audit['inhibition_stats']['min']}-{audit['inhibition_stats']['max']}")
    print(f"  All SOD1: {audit['all_sod1']}")
    print("[Phase 1] COMPLETE [OK]\n")

    return df, audit


if __name__ == "__main__":
    main()
