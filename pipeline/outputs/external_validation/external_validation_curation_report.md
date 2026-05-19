# Phase 6 -- External Validation Curation Report

**Generated:** 2026-05-13T20:17:32.317395

## Summary

- **External ASOs curated:** 1
- **Tofersen included:** Yes

## Curation Rules

1. Only ASOs with verified sequences from published sources are included.
2. No inhibition values are fabricated.
3. Missing experimental conditions are imputed using training dataset medians/modes.
4. Each entry includes source citation.

## Standardized Conditions (for missing values)

| Parameter | Value |
|---|---|
| cell_line | HepG2 |
| transfection | electroporation |
| aso_concentration_nm | 3000.0 |
| treatment_period_hours | 16.0 |
| density_cells_per_well | 20000 |
| primer_probe_set | RTS3898 |

## Additional ASOs for Future Curation

The following ASOs should be added when verified sequences become available:

- ISIS 333611 -- Early SOD1 ASO from Ionis, Phase 1 (Miller et al. 2013 Lancet Neurol). Sequence available in patent WO2007002390.
- ISIS 666853 -- Predecessor to tofersen with similar gapmer design. Referenced in Ionis pipeline.
- Nusinersen (Spinraza) -- FDA-approved splice-switching ASO (non-SOD1, different mechanism, 2'-MOE uniform).
- Mipomersen (Kynamro) -- FDA-approved RNase H gapmer (targets ApoB, 5-10-5 MOE).
- Inotersen (Tegsedi) -- FDA-approved RNase H gapmer (targets TTR, 5-10-5 MOE).
- Volanesorsen -- Approved RNase H gapmer (targets ApoC-III).
- Pelacarsen -- Clinical-stage Lp(a)-targeting gapmer (cEt wings).

## Limitations

1. Only tofersen has a fully verified sequence from the provided FDA GSRS JSON.
2. Other SOD1 ASOs from literature require manual sequence verification.
3. External ASO inhibition values may not be directly comparable to the internal dataset assay.
4. Clinical efficacy endpoints differ from in-vitro mRNA knockdown measured in our dataset.
