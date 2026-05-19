# Phase 9 -- Tofersen Benchmarking Report

**Generated:** 2026-05-13T20:17:55.530717

> ⚠️ **Disclaimer:** This is an in-silico benchmark and not a clinical efficacy claim.

## Tofersen Reference

- **Drug name:** Tofersen (Qalsody™)
- **Manufacturer:** Biogen
- **FDA approval:** April 2023 (accelerated approval)
- **Target:** SOD1 mRNA
- **Mechanism:** RNase H-mediated mRNA degradation
- **Chemistry:** 5-10-5 MOE gapmer
- **Sequence:** CAGGATACATTTCTACAGCT (DNA equivalent)
- **CAS:** 1898254-60-8
- **UNII:** 5YL205692C
- **Source:** FDA GSRS database, USAN Council, Qalsody prescribing information

## Standardized Prediction Conditions

The following conditions were derived from the training dataset to enable
fair comparison:

| Parameter | Value |
|---|---|
| cell_line | HepG2 |
| transfection | electroporation |
| aso_concentration_nm | 3000.0 |
| treatment_period_hours | 16.0 |
| density_cells_per_well | 20000 |
| primer_probe_set | RTS3898 |

## Prediction Result

- **Predicted inhibition:** 57.54%

## Benchmark Usage

For any novel ASO prediction, the output will include:

1. **Predicted inhibition (%)** for the novel ASO
2. **Tofersen predicted inhibition:** 57.54% (under standardized conditions)
3. **Delta vs tofersen:** predicted_ASO - predicted_tofersen
4. **Category:**
   - **Lower than tofersen:** delta < -10.0
   - **Comparable to tofersen:** |delta| <= 10.0
   - **Higher than tofersen:** delta > +10.0

## Important Caveats

1. Tofersen's clinical efficacy was demonstrated through intrathecal delivery to ALS patients.
2. Our model predicts in-vitro mRNA inhibition in cell lines (HepG2, A431, SH-SY5Y).
3. In-vitro potency does not directly translate to in-vivo efficacy.
4. The comparison is meaningful only within the model's domain of applicability.
5. This benchmark does not account for pharmacokinetic, safety, or delivery differences.
