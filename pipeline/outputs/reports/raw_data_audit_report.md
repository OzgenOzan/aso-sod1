# Phase 1 -- Raw Data Audit Report
**Generated:** 2026-05-13T20:17:22.855596

## 1. Dataset Overview
- **Rows:** 2264
- **Columns:** 17
- **Unique ASO sequences:** 688
- **Unique sequence + chemical pattern combinations:** 1008
- **Unique cell lines:** 3 -- ['A431', 'HepG2', 'SH-SY5Y']
- **Unique transfection methods:** 3 -- ['electroporation', 'free uptake', 'uptake']
- **Unique chemical patterns:** 20
- **Unique linkage patterns:** 2 -- ['phosphodiester/phosphorothioate', 'phosphorothioate']

## 2. Missingness per Column
| Column | Missing |
|---|---|
| ISIS | 0 |
| Target_gene | 0 |
| Cell_line | 0 |
| Density(cells/well) | 0 |
| Transfection | 0 |
| ASO_volume(nM) | 0 |
| Treatment_Period(hours) | 0 |
| Primer_probe_set | 0 |
| Sequence | 0 |
| Modification | 0 |
| Location | 0 |
| Chemical_Pattern | 0 |
| Linkage | 0 |
| Linkage_Location | 0 |
| Smiles | 0 |
| Inhibition(%) | 0 |
| seq_length | 0 |

## 3. Inhibition (%) Summary Statistics
| Statistic | Value |
|---|---|
| count | 2264 |
| mean | 42.05 |
| std | 28.96 |
| min | 0.0 |
| 25% | 17.0 |
| 50% | 38.0 |
| 75% | 67.0 |
| max | 98.0 |

## 4. Sequence Length Distribution
- **Range:** 10 - 20
- **Mean:** 17.72
| Length | Count |
|---|---|
| 10 | 12 |
| 16 | 232 |
| 17 | 1259 |
| 18 | 169 |
| 20 | 592 |

## 5. Target Gene Validation
- **Target genes found:** ['SOD-1']
- **All target SOD1:** [OK] Yes

## 6. Inhibition Range Validation
- **Numeric type:** [OK] Yes
- **Values outside [0, 100]:** 0

## 7. Sequence Length Verification
- **Mismatches (seq_length != len(Sequence)):** 12

| ISIS | Sequence | Stated Length | Computed Length |
|---|---|---|---|
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |
| 666867 | AGTGTTTAATGTTTATC | 10 | 17 |

## 8. Duplicate Analysis
- **Fully duplicated rows:** 20
- **Rows with duplicated sequences (same Sequence, different conditions):** 1576
- **Rows with duplicated Sequence + Chemical_Pattern:** 1256
- **Rows with duplicated Sequence + Chemical_Pattern + Linkage:** 1241

## 9. Conflicting Inhibition Values
- **Groups with conflicting inhibition under identical conditions:** 43
- **Rows involved:** 149

## 10. Chemical Patterns Observed
| Pattern | Length | Wing-Gap-Wing |
|---|---|---|
| `CMCMddddddddMCMC` | 16 | 4-8-4 |
| `MCMCddddddddCMCMM` | 17 | 4-8-5 |
| `MCMCddddddddMMMMM` | 17 | 4-8-5 |
| `MCddddddddMCMCMM` | 16 | 2-8-6 |
| `MMCCddddddddCCMMM` | 17 | 4-8-5 |
| `MMCCddddddddMMMMM` | 17 | 4-8-5 |
| `MMCCdddddddddCCMM` | 17 | 4-9-4 |
| `MMMCCdddddddCCMMM` | 17 | 5-7-5 |
| `MMMCddddddddCMMMM` | 17 | 4-8-5 |
| `MMMMCddddddddCMMM` | 17 | 5-8-4 |
| `MMMMMMddddddddMMMMMM` | 20 | 6-8-6 |
| `MMMMMMdddddddddMMMMM` | 20 | 6-9-5 |
| `MMMMMddddddddCCMM` | 17 | 5-8-4 |
| `MMMMMddddddddMMMMM` | 18 | 5-8-5 |
| `MMMMMddddddddMMMMMMM` | 20 | 5-8-7 |
| `MMMMMddddddddddMMMMM` | 20 | 5-10-5 |
| `MMMMddddddddCCMMM` | 17 | 4-8-5 |
| `MMMMddddddddCMCMM` | 17 | 4-8-5 |
| `MMMMddddddddMMMMM` | 17 | 4-8-5 |
| `MMMMdddddddddCCMM` | 17 | 4-9-4 |
