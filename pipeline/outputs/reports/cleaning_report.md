# Phase 2 -- Data Cleaning and Validation Report

**Generated:** 2026-05-13T20:17:24.685744

## Summary

- **Cleaned dataset rows:** 2155
- **Excluded/flagged records:** 12

## Processing Steps

1. Column renaming to snake_case
2. Sequence normalization (uppercase, whitespace removal, ACGT validation)
3. Sequence length verification and correction
4. Categorical variable normalization
5. Numeric variable validation
6. Chemical pattern validation
7. Location column parsing
8. Linkage location parsing
9. SMILES validation (limited -- RDKit not available)
10. Duplicate removal and replicate aggregation

## Flags and Issues

### sequence_normalization

No issues found.

### sequence_length

- {'index': 526, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1001, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1022, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1613, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1914, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1936, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1957, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1978, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 1999, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 2209, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 2230, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}
- {'index': 2251, 'flag': 'seq_length_corrected', 'stated': 10, 'computed': 17, 'sequence': 'AGTGTTTAATGTTTATC'}

### categorical_normalization

- {'flag': 'transfection_normalized', 'note': "'uptake' and 'free uptake' merged to 'free_uptake' (555 rows)"}

### numeric_validation

No issues found.

### chemical_pattern

No issues found.

### location_parsing

No issues found.

### linkage_location_parsing

No issues found.

### smiles_validation

- {'flag': 'rdkit_not_available', 'note': 'RDKit is not installed. SMILES strings are accepted but not chemically validated.'}

### duplicate_handling

- {'flag': 'exact_duplicates_removed', 'count': 20}
- {'flag': 'replicate_aggregation', 'groups': 2155, 'original_rows': 2244}

## Assumptions

1. 'uptake' and 'free uptake' transfection methods are treated as equivalent ('free_uptake').
2. Sequence length mismatches (12 records with stated length 10 but actual length 17) are corrected to actual length.
3. SMILES strings are accepted without chemical validation (RDKit not available).
4. Biological replicates under identical conditions are aggregated by mean inhibition.
5. All sequences contain only standard DNA bases (A, C, G, T).
