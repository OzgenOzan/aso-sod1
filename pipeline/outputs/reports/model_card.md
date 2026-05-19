# Model Card -- ASO Inhibition Predictor

**Generated:** 2026-05-13T20:17:55.125914

## Model Details

- **Model type:** LightGBM
- **Features:** 335 engineered features
- **Target:** inhibition_percent (0-100, bounded regression)
- **Training data:** SOD1-targeting ASOs (~2,264 records, ~688 unique sequences)
- **Random seed:** 42

## Intended Use

- In-silico prediction of mRNA inhibition percentage for novel SOD1-targeting ASOs.
- Research use only. Not for clinical decision-making.

## Performance

| Metric | Test Set |
|---|---|
| MAE | 12.4168 |
| RMSE | 15.8569 |
| R2 | 0.6726 |
| Pearson_r | 0.8301 |
| Within_10 | 48.3923 |

## Limitations

1. Trained only on SOD1-targeting ASOs -- may not generalize to other targets.
2. Limited cell lines (HepG2, A431, SH-SY5Y).
3. No transcript-position features due to missing coordinate data.
4. RDKit unavailable -- no molecular descriptor features from SMILES.
5. Predictions should be validated experimentally.

## Ethical Considerations

This model provides research-use-only predictions and should not replace
experimental validation or clinical assessment.
