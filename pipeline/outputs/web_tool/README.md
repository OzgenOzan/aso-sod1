# SOD1 ASO Inhibition Predictor — Web Tool

## Overview
This Streamlit web application predicts the mRNA inhibition percentage of novel
antisense oligonucleotides (ASOs) targeting the human SOD1 transcript and
benchmarks predictions against tofersen (Qalsody).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

## Input Parameters
- **ASO Sequence** — DNA nucleotides (A, C, G, T), 5' to 3'
- **Chemical Pattern** — Per-position chemistry (M=MOE, C=cEt, d=deoxy)
- **Modification** — Chemical modification description
- **Linkage** — Backbone linkage type
- **Cell Line** — HepG2, A431, or SH-SY5Y
- **Transfection** — Electroporation or free uptake
- **ASO Concentration** — In nanomolar (nM)
- **Treatment Period** — In hours

## Model
- **Algorithm:** LightGBM regressor
- **Features:** 335 engineered features from sequence, chemistry, linkage, and conditions
- **Training data:** ~2,155 SOD1-targeting ASO records
- **Test performance:** R^2 = 0.67, MAE = 12.4

## Disclaimer
This tool provides research-use-only in-silico predictions. It does not provide
clinical, therapeutic, or regulatory advice. Experimental validation is required.
