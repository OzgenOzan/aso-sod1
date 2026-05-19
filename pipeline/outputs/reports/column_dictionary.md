# Column Dictionary

| Original Column (-> Standardized) | Type | Description |
|---|---|---|
| ISIS (-> no) | Integer | Unique ISIS identifier for the ASO compound |
| Target_gene (-> target_gene) | Categorical | Target gene symbol (SOD-1 in this dataset) |
| Cell_line (-> cell_line) | Categorical | Cell line used in the assay (HepG2, A431, SH-SY5Y) |
| Density(cells/well) (-> density_cells_per_well) | Numeric | Cell seeding density per well |
| Transfection (-> transfection) | Categorical | Transfection or delivery method (electroporation, uptake, free uptake) |
| ASO_volume(nM) (-> aso_concentration_nm) | Numeric | ASO concentration in nanomolar |
| Treatment_Period(hours) (-> treatment_period_hours) | Numeric | Duration of ASO treatment in hours |
| Primer_probe_set (-> primer_probe_set) | Categorical | qPCR primer/probe set used for mRNA quantification |
| Sequence (-> sequence) | Text/Sequence | ASO nucleotide sequence (5'->3', DNA alphabet A/C/G/T) |
| Modification (-> modification) | Categorical | High-level chemical modification description |
| Location (-> location) | Text/Structured | Position indices for modifications, format: pos1?pos2?.../type/else |
| Chemical_Pattern (-> chemical_pattern) | Text/Categorical | Per-position chemistry code string: M=MOE, C=cEt, d=deoxy |
| Linkage (-> linkage) | Categorical | Backbone linkage type(s): phosphorothioate and/or phosphodiester |
| Linkage_Location (-> linkage_location) | Text/Structured | Position indices for linkage types |
| Smiles (-> smiles) | Text/Chemical | SMILES representation of the full ASO molecule |
| Inhibition(%) (-> inhibition_percent) | Numeric (0-100) | Target variable: mRNA inhibition percentage |
| seq_length (-> seq_length) | Numeric | Stated sequence length in nucleotides |

## Chemistry Code Definitions

| Code | Meaning |
|---|---|
| `M` | 2'-MOE (2'-O-methoxyethyl) modified nucleotide |
| `C` | cEt (constrained ethyl / S-cEt BNA) modified nucleotide |
| `d` | 2'-deoxy (unmodified DNA) nucleotide |
