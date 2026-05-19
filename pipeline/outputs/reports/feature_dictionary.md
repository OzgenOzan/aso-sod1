# Feature Dictionary

**Total features:** 336  (excluding sequence, aso_group_id, sequence_cluster_id, and inhibition_percent)

## A. Sequence Composition (282 features)

- `seq_length`
- `A_count`
- `A_fraction`
- `C_count`
- `C_fraction`
- `G_count`
- `G_fraction`
- `T_count`
- `T_fraction`
- `GC_fraction`
- `AT_fraction`
- `purine_fraction`
- `pyrimidine_fraction`
- `di_AA`
- `di_AA_freq`
- `di_AC`
- `di_AC_freq`
- `di_AG`
- `di_AG_freq`
- `di_AT`
- `di_AT_freq`
- `di_CA`
- `di_CA_freq`
- `di_CC`
- `di_CC_freq`
- `di_CG`
- `di_CG_freq`
- `di_CT`
- `di_CT_freq`
- `di_GA`
- `di_GA_freq`
- `di_GC`
- `di_GC_freq`
- `di_GG`
- `di_GG_freq`
- `di_GT`
- `di_GT_freq`
- `di_TA`
- `di_TA_freq`
- `di_TC`
- `di_TC_freq`
- `di_TG`
- `di_TG_freq`
- `di_TT`
- `di_TT_freq`
- `tri_AAA`
- `tri_AAC`
- `tri_AAG`
- `tri_AAT`
- `tri_ACA`
- ... and 232 more

## B. Approximate Thermodynamic (6 features)

- `Tm_wallace`
- `Tm_gc_adjusted`
- `self_complementarity_score`
- `rc_internal_match_score`
- `palindrome_score`
- `terminal_gc_clamp`

## C. Chemistry-Pattern (23 features)

- `chemical_pattern_length`
- `count_M`
- `fraction_M`
- `count_C`
- `fraction_C`
- `count_d`
- `fraction_d`
- `left_modified_wing_length`
- `right_modified_wing_length`
- `central_deoxy_gap_length`
- `symmetry_flag`
- `n_chemistry_transitions`
- `gapmer_5-7-5`
- `gapmer_5-10-5`
- `gapmer_4-8-5`
- `gapmer_4-9-4`
- `gapmer_2-8-6`
- `gapmer_4-8-4`
- `gapmer_5-8-5`
- `gapmer_5-8-7`
- `gapmer_5-8-4`
- `gapmer_6-8-6`
- `gapmer_6-9-5`

## D. Modification-Text (6 features)

- `contains_MOE`
- `contains_cEt`
- `contains_5mc`
- `contains_deoxy`
- `contains_LNA`
- `n_modification_types`

## E. Modification-Location (5 features)

- `n_location_groups`
- `n_position_indices`
- `min_modified_position`
- `max_modified_position`
- `mod_positions_valid`

## F. Backbone-Linkage (10 features)

- `contains_PS`
- `contains_PO`
- `linkage_type_count`
- `predicted_PS_count`
- `predicted_PO_count`
- `predicted_PS_fraction`
- `predicted_PO_fraction`
- `linkage_positions_valid`
- `n_linkage_position_indices`
- `linkage_transition_count`

## G. SMILES-Derived (Limited) (3 features)

- `smiles_length`
- `smiles_has_stereo`
- `smiles_ring_count_approx`

## H. Experimental-Condition (10 features)

- `log10_aso_concentration_nm`
- `treatment_period_hours`
- `log10_density_cells_per_well`
- `cell_A431`
- `cell_HepG2`
- `cell_SH-SY5Y`
- `transfection_electroporation`
- `transfection_free_uptake`
- `primer_HTS90`
- `primer_RTS3898`

