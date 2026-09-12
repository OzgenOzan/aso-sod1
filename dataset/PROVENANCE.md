# Dataset Provenance

> **Status: TEMPLATE — owner must complete the TODO fields before relying on or
> redistributing these data.** Added by the 2026-08 repository audit (finding F5).

## `sod-1_data.xlsx` (training data, 2,264 records)

| Field | Value |
|---|---|
| Primary source | **TODO** — e.g. Ionis Pharmaceuticals internal assay compendium? Published dataset? |
| Citation / accession | **TODO** |
| Download / export date | **TODO** |
| License / usage terms | **TODO** — the `ISIS#######` compound identifiers follow Ionis Pharmaceuticals naming; confirm redistribution rights for this public repository |
| Preprocessing / filtering chain | **TODO** — document how 2,264 records were derived from the source |
| Contact for questions | **TODO** |

## `tofersen.json` (reference compound)

| Field | Value |
|---|---|
| Structure source | **TODO** — e.g. published tofersen (BIIB067 / Qalsody) sequence/chemistry reference |
| Measured inhibition value | **None exists** — `pipeline/outputs/external_validation/tofersen_reference.json`
  contains a *model self-prediction* under imputed conditions (HepG2, electroporation,
  3000 nM, 16 h — training-set medians/modes), not an experimental measurement. |

## Why this matters

- The training data's origin, license, and redistribution rights are currently undocumented
  in a public repository (audit finding F5 — reproducibility and legal/attribution risk).
- Do not cite or reuse `sod-1_data.xlsx` until the TODO fields are filled in.
