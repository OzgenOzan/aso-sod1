# ASO Inhibition Pipeline: Complete Feature Documentation

This document explains **all** features extracted by the pipeline (`phase4_features.py`). Features are logically grouped into categories. For each feature (or feature group), we describe what it represents, how it is calculated, and the specific python logic used for extraction.

---

## Group A: Sequence Composition

These features describe the nucleotide makeup, specific structural motifs, and binding anchors of the ASO sequence.

### Basic Length and Base Fractions
* **`seq_length`**
  * **Significance:** Overall length of the ASO. Affects binding affinity and potential for secondary structures.
  * **Code:** `len(seq)`
* **`{base}_count` & `{base}_fraction` (for A, C, G, T)**
  * **Significance:** The absolute count and relative frequency of each nucleotide. Different balances affect backbone rigidity and immunostimulatory properties.
  * **Code:** `seq.count(base)` and `count / seq_length`
* **`GC_fraction` & `AT_fraction`**
  * **Significance:** Total fraction of strong (GC, 3 hydrogen bonds) vs weak (AT, 2 hydrogen bonds) bases. Higher GC increases binding affinity but can increase off-target toxicity.
  * **Code:** `(G_count + C_count) / seq_length`
* **`purine_fraction` & `pyrimidine_fraction`**
  * **Significance:** Fraction of two-ring (Purines: A, G) vs one-ring (Pyrimidines: C, T) bases. High purine content can alter backbone flexibility.
  * **Code:** `(A_count + G_count) / seq_length`

### Di- and Tri-Nucleotide Motifs
* **`di_{motif}` & `di_{motif}_freq`** (16 combinations, e.g., `di_CG`)
  * **Significance:** Frequencies of 2-mer sliding windows. Specifically, `CpG` (CG) dinucleotides can trigger innate immune responses (TLR9).
  * **Code:** `seq.count(di)` and `count / (seq_length - 1)`
* **`tri_{motif}`** (64 combinations, e.g., `tri_TCA`)
  * **Significance:** Absolute counts of 3-mer motifs. Specific motifs might be recognition sites for nucleases or RNA-binding proteins.
  * **Code:** `seq.count(tri)`
* **Specific Motif Counts (`CpG_count`, `GG_count`, `CC_count`)**
  * **Significance:** `CpG` is immunostimulatory. Poly-G (`GG`) runs can form G-quadruplex structures which aggregate and are highly toxic.
  * **Code:** `seq.count("CG")`, `seq.count("GG")`, `seq.count("CC")`

### Contiguous Structural Blocks
* **`max_homopolymer_run`**
  * **Significance:** The longest contiguous stretch of a single repeated base (e.g., "AAAA"). Long runs cause non-specific protein binding, manufacturing difficulties, and off-target hybridization.
  * **Code:** Iterates through the sequence incrementing a counter when `seq[i] == seq[i-1]`.
* **`max_gc_stretch`**
  * **Significance:** The longest contiguous run of G or C bases. Represents the most rigid, strongest binding "anchor" segment of the ASO.
  * **Code:** Iterates through the sequence incrementing a counter whenever `c in "GC"`.

### Terminal Anchors
* **`five_prime_base_{X}` & `three_prime_base_{X}`** (One-hot, 8 features)
  * **Significance:** The exact base at the 5' and 3' ends. The ends are critical for exonuclease resistance and dictate the initial nucleation of RNA-binding.
  * **Code:** `1 if seq[0] == 'X' else 0` and `1 if seq[-1] == 'X' else 0`
* **`first_2nt_{XX}`, `last_2nt_{XX}`, `first_3nt_{XXX}`, `last_3nt_{XXX}`**
  * **Significance:** Extended sequence context at the extreme ends.
  * **Code:** `1 if seq[:2] == 'XX' else 0`, etc.

---

## Group B: Approximate Thermodynamic

Features approximating the ASO's thermal stability and propensity to form stable, unintended secondary structures (hairpins or self-dimers).

* **`Tm_wallace`**
  * **Significance:** Wallace rule melting temperature ($T_m$). Best for short oligos (<20nt). Higher $T_m$ means stronger binding.
  * **Code:** `2 * (A+T) + 4 * (G+C)`
* **`Tm_gc_adjusted`**
  * **Significance:** Marmur-Doty formula for $T_m$. Better suited for longer sequences and accounts for length dynamically.
  * **Code:** `64.9 + 41.0 * (gc - 16.4) / n`
* **`self_complementarity_score`**
  * **Significance:** Propensity to form a self-dimer or hairpin loop. High scores indicate the ASO might bind itself instead of the target mRNA, reducing efficacy.
  * **Code:** Maximum contiguous match of any substring (length $\ge3$) against the sequence's reverse complement.
* **`rc_internal_match_score`**
  * **Significance:** The overall fraction of the sequence that aligns with its own reverse complement when laid end-to-end.
  * **Code:** `sum(1 for a, b in zip(seq, rc) if a == b) / len(seq)`
* **`palindrome_score`**
  * **Significance:** The number of internal 4- to 8-mer motifs that are exact palindromes, indicating micro-hairpin potential.
  * **Code:** Sliding window (4 to 8 nt); increments score if the window equals its reverse complement.
* **`terminal_gc_clamp`**
  * **Significance:** Binary indicator of whether both ends are anchored by strong G/C bonds, preventing end-fraying.
  * **Code:** `int(seq[0] in "GC" and seq[-1] in "GC")`

---

## Group C: Chemistry-Pattern (Gapmer Architecture)

Features derived from the `chemical_pattern` string (e.g., `MMMMddddddddddMMMM`), outlining the gapmer design.

* **`chemical_pattern_length`**
  * **Significance:** Total length of the chemistry string (should match `seq_length`).
  * **Code:** `len(cp)`
* **`count_M`, `count_C`, `count_d` & `fraction_M`, `fraction_C`, `fraction_d`**
  * **Significance:** Counts and fractions of specific modifications. `d` = Unmodified DNA (Gap), `M` = 2'-MOE (Wing), `C` = cEt (Wing). The balance dictates RNase H recruitment (via DNA gap) vs binding affinity/stability (via wings).
  * **Code:** `cp.count('M')`, etc.
* **`left_modified_wing_length` & `right_modified_wing_length`**
  * **Significance:** The size of the 5' and 3' protective wings. Longer wings increase nuclease resistance but reduce the size of the catalytic gap.
  * **Code:** Counts contiguous non-'d' characters from the start and end of the string.
* **`central_deoxy_gap_length`**
  * **Significance:** The size of the unmodified DNA gap. Must typically be $\ge 8$ nt to efficiently recruit RNase H.
  * **Code:** `cp.count('d')`
* **`symmetry_flag`**
  * **Significance:** 1 if the left and right wings are equal in length, 0 otherwise.
  * **Code:** `int(left_wing == right_wing)`
* **`n_chemistry_transitions`**
  * **Significance:** Number of times the modification changes. For a standard gapmer, this is exactly 2. Higher numbers indicate alternating architectures (e.g., gap-mixmers).
  * **Code:** `sum(1 for j in range(1, len(cp)) if cp[j] != cp[j - 1])`
* **`gapmer_{arch}`** (One-hot encoded)
  * **Significance:** Specific architecture classes (e.g., `gapmer_5-10-5`).
  * **Code:** Dummy variables generated from `f"{left}-{gap}-{right}"`

---

## Group D: Modification-Text & Group E: Modification-Location

Features parsing explicit textual descriptions of chemical modifications.

* **`contains_MOE`, `contains_cEt`, `contains_5mc`, `contains_deoxy`, `contains_LNA`**
  * **Significance:** Binary flags indicating the presence of specific chemical families. E.g., 5-methylcytosine (`5mc`) reduces immunostimulation.
  * **Code:** Substring checks on lowercased `modification` text (e.g., `"moe" in text`).
* **`n_modification_types`**
  * **Significance:** Complexity of the ASO design (how many different chemistries are combined).
  * **Code:** Splits modification text by `/` and counts elements.
* **`n_location_groups` & `n_position_indices`**
  * **Significance:** Number of distinct clusters of modifications vs the absolute number of modified bases.
  * **Code:** Splits the `location` text by `/` and parses integer digits.
* **`min_modified_position` & `max_modified_position`**
  * **Significance:** Identifies the spatial bounds of the modifications (typically spanning the entire sequence).
  * **Code:** `min(all_positions)` and `max(all_positions)`
* **`mod_positions_valid`**
  * **Significance:** Data quality flag ensuring modifications do not fall out of the sequence bounds.
  * **Code:** `1` if `0 <= pos < seq_length` for all parsed positions.

---

## Group F: Backbone-Linkage

Features determining the internucleotide linkage chemistry (Phosphorothioate vs Phosphodiester).

* **`contains_PS` & `contains_PO`**
  * **Significance:** PS (Phosphorothioate) replaces an oxygen with sulfur, drastically improving nuclease resistance but increasing toxicity and reducing binding affinity. PO is natural DNA.
  * **Code:** Checks for `"phosphorothioate"` or `"phosphodiester"` in the `linkage` text.
* **`linkage_type_count`**
  * **Significance:** Usually 1 (fully PS) or 2 (mixed PS/PO backbone).
  * **Code:** Counts segments split by `/`.
* **`predicted_PS_count`, `predicted_PO_count` & their fractions**
  * **Significance:** The absolute number of linkages of each type. A fully PS backbone has `seq_length - 1` PS linkages. Mixed backbones intersperse PO to reduce toxicity.
  * **Code:** Parses numeric indices from `linkage_location` to count PO positions; PS count is `(seq_length - 1) - PO_count`.
* **`linkage_transition_count`**
  * **Significance:** The number of times the backbone switches between PS and PO.
  * **Code:** Approximated by multiplying the number of PO linkage clusters by 2 (entering and exiting the PO region).

---

## Group G: SMILES-Derived Proxies

*Because RDKit is not assumed to be installed, these are basic string-level proxy features.*

* **`smiles_length`**
  * **Significance:** Strong proxy for total molecular weight and atom count.
  * **Code:** `len(smiles)`
* **`smiles_has_stereo`**
  * **Significance:** 1 if the SMILES defines stereochemistry (chirality). Certain chiral linkages affect RNase H activity differently.
  * **Code:** Regex search for `@@` or `@` in the string.
* **`smiles_ring_count_approx`**
  * **Significance:** Proxy for the number of cyclic structures (ribose/base rings).
  * **Code:** Counts numbers in the SMILES string divided by 2.

---

## Group H: Experimental Conditions

Features detailing the biological assay setup.

* **`log10_aso_concentration_nm`**
  * **Significance:** The dosage of ASO. Log-transformed due to massive scale differences (e.g., 1nM vs 10,000nM). Crucial for dose-response predictions.
  * **Code:** `np.log10(concentration)`
* **`treatment_period_hours`**
  * **Significance:** Time cells were incubated with the ASO. Efficacy changes dynamically over time as ASO is taken up and mRNA is degraded.
  * **Code:** Direct value from `treatment_period_hours`.
* **`log10_density_cells_per_well`**
  * **Significance:** Cell confluence. High density can alter ASO uptake efficiency and cell proliferation rates.
  * **Code:** `np.log10(density)`
* **`cell_{cell_line}`** (One-hot)
  * **Significance:** Different cell lines (e.g., SH-SY5Y neuroblastoma vs HepG2 hepatoma) have different internal concentrations of RNase H and different surface receptors for ASO uptake.
  * **Code:** Dummy columns generated from `cell_line`.
* **`transfection_{tf}`** (One-hot)
  * **Significance:** Method of delivery (e.g., Gymnosis/free-uptake vs Lipofectamine). Transfection agents artificially boost ASO uptake, inflating apparent efficacy.
  * **Code:** Dummy columns generated from `transfection`.
* **`primer_{pp}`** (One-hot)
  * **Significance:** The specific qPCR primer probe set used to measure target mRNA. Different primers might target different mRNA isoforms or splice variants.
  * **Code:** Dummy columns generated from `primer_probe_set`.
