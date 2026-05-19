# Phase 4 -- Feature Extraction Report

**Generated:** 2026-05-13T20:17:31.513533

## Summary

- **Total samples:** 2155
- **Total features:** 339
- **Feature groups:** 8

| Group | Features |
|---|---|
| A. Sequence Composition | 282 |
| B. Approximate Thermodynamic | 6 |
| C. Chemistry-Pattern | 23 |
| D. Modification-Text | 6 |
| E. Modification-Location | 5 |
| F. Backbone-Linkage | 10 |
| G. SMILES-Derived (Limited) | 3 |
| H. Experimental-Condition | 10 |

## Skipped Features

### SMILES Molecular Descriptors
- **Reason:** RDKit is not installed in the current environment.
- **Impact:** Molecular weight, TPSA, rotatable bonds, and ring count from SMILES are unavailable.
- **Mitigation:** Basic SMILES string-level features (length, stereo count, approximate ring count) were extracted.

### SOD1 Transcript Mapping
- **Reason:** Dataset does not contain transcript coordinates. External mapping was not performed.
- **Impact:** Target position, UTR/CDS/intron annotation, and RNA accessibility features are unavailable.
- **Mitigation:** None -- these features are documented as absent.

## Thermodynamic Approximation Formulas

- **Wallace Tm:** Tm = 2x(A+T) + 4x(G+C) -- suitable for short oligos <20nt
- **GC-adjusted Tm (Marmur-Doty):** Tm = 64.9 + 41x(GC_count - 16.4) / length
- Note: These are rough approximations. Nearest-neighbor Tm would be more accurate but requires specialized libraries.
