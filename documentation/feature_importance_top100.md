# ASO Inhibition Pipeline: Top 100 Features Documentation

This document explains the top 100 most important features identified by the LightGBM model in the ASO Knockdown Prediction Pipeline. Features are grouped logically by their biological, structural, or experimental category. For each feature, we detail what it represents, how it is calculated, and how the underlying codebase performs the extraction.

## 1. SMILES-Derived Proxies (Molecular Descriptors)
*SMILES (Simplified Molecular-Input Line-Entry System) strings represent the full 3D chemical structure in 1D text. Since RDKit is not used in the pipeline, these are approximate proxy features derived purely from string parsing.*

* **`smiles_length`** (Rank 1)
  * **What it represents:** A proxy for the overall molecular weight and size of the ASO.
  * **How to calculate:** The total character count of the SMILES string.
  * **How the code works:** `df["smiles"].astype(str).str.len()`

* **`smiles_ring_count_approx`** (Rank 3)
  * **What it represents:** An approximate count of cyclical ring structures in the molecule (such as ribose rings and nucleobase rings).
  * **How to calculate:** Number of digits present in the SMILES string divided by 2. (In SMILES, digits are used in pairs to denote the opening and closing of a ring).
  * **How the code works:** Uses a lambda function to count numeric characters and applies integer division by 2: `sum(1 for c in s if c.isdigit()) // 2`.

## 2. Experimental Conditions
*These features encapsulate the assay conditions under which the ASO inhibition was measured.*

* **`log10_aso_concentration_nm`** (Rank 2)
  * **What it represents:** The concentration of the ASO dosed in the experiment. Log-transformed due to vast scaling differences across experiments.
  * **How to calculate:** `log10` of the concentration in nanomolar (nM), clipped to a minimum of 1.
  * **How the code works:** `np.log10(df["aso_concentration_nm"].clip(lower=1))`

* **`log10_density_cells_per_well`** (Rank 6)
  * **What it represents:** The number of cells plated per well, reflecting cell confluence.
  * **How to calculate:** `log10` of the cell density.
  * **How the code works:** `np.log10(df["density_cells_per_well"].clip(lower=1))`

* **`treatment_period_hours`** (Rank 7)
  * **What it represents:** The duration cells were exposed to the ASO before mRNA knockdown was measured.
  * **How to calculate:** Taken directly from the dataset in hours.
  * **How the code works:** Directly copied from `df["treatment_period_hours"]`.

* **`primer_HTS90`** (Rank 29)
  * **What it represents:** Indicates if the HTS90 primer/probe set was used for qPCR quantification.
  * **How to calculate:** One-hot encoded binary indicator (1 if true, 0 if false).
  * **How the code works:** Boolean equality check against `"primer_probe_set"`, cast to an integer.

* **`cell_SH-SY5Y`** (Rank 52) & **`cell_HepG2`** (Rank 57)
  * **What it represents:** Indicates if the assay was performed in the SH-SY5Y neuroblastoma cell line or the HepG2 hepatoma cell line.
  * **How to calculate:** One-hot encoded binary indicator.
  * **How the code works:** Boolean equality check against `"cell_line"`, cast to an integer.

## 3. Thermodynamics & Structure
*Approximations of the ASO's binding affinity to the target mRNA and its propensity to form secondary structures or self-dimerize.*

* **`rc_internal_match_score`** (Rank 4)
  * **What it represents:** A measure of how much the sequence aligns with its own reverse complement. High scores suggest the ASO might fold onto itself or bind to other ASO molecules instead of the target mRNA.
  * **How to calculate:** The fraction of positions that perfectly match when the sequence is aligned start-to-end with its reverse complement.
  * **How the code works:** Computes `sum(1 for a, b in zip(seq, rc) if a == b) / n`, where `rc` is the reverse complement.

* **`self_complementarity_score`** (Rank 27)
  * **What it represents:** The maximum contiguous hairpin or self-dimer potential.
  * **How to calculate:** The length of the longest contiguous substring in the sequence that also exists in the reverse complement (minimum length 3).
  * **How the code works:** Uses nested loops to check all substrings of length $\ge 3$ against the reverse complement string, keeping track of the max length.

* **`Tm_wallace`** (Rank 28)
  * **What it represents:** An approximate melting temperature ($T_m$) using the Wallace rule, best for short oligos.
  * **How to calculate:** $T_m = 2 \times (A + T) + 4 \times (G + C)$
  * **How the code works:** Simply multiplies the pre-calculated counts of AT by 2 and GC by 4.

* **`Tm_gc_adjusted`** (Rank 89)
  * **What it represents:** Melting temperature approximation using the Marmur-Doty formula, better suited for longer sequences.
  * **How to calculate:** $T_m = 64.9 + 41.0 \times (GC\_count - 16.4) / length$
  * **How the code works:** Straightforward mathematical application of the formula using the ASO length and GC count.

* **`palindrome_score`** (Rank 41)
  * **What it represents:** The propensity to form palindromic secondary structures.
  * **How to calculate:** The number of 4- to 8-mer motifs in the sequence that are exact palindromes (i.e., sequence equals its reverse complement).
  * **How the code works:** Uses a sliding window of sizes 4 through 8, extracts the sub-sequence, generates its reverse complement, and increments the score if they match.

* **`terminal_gc_clamp`** (Rank 81)
  * **What it represents:** Indicates if the sequence has strong binding anchors (G or C) at both the extreme 5' and 3' ends.
  * **How to calculate:** 1 if both the first and last nucleotides are G or C, otherwise 0.
  * **How the code works:** Checks `seq[0] in "GC"` and `seq[-1] in "GC"`, returning the boolean AND.

## 4. Sequence Composition
*Features describing the base constitution and sub-motifs within the nucleotide sequence.*

### Base Fractions and Counts
* **What they represent:** The raw abundance of specific nucleotides. For instance, high GC content generally improves binding affinity but may increase off-target effects. High purine content can alter backbone rigidity.
* **How they are calculated:** 
  * Counts: Number of occurrences of the specific base(s).
  * Fractions: The count divided by the total sequence length.
* **Code extraction:** Uses Python's native `string.count(base)` and division by `seq_length`.
* **Top Features in this sub-group:**
  * **`purine_fraction`** (Rank 8) (A + G)
  * **`GC_fraction`** (Rank 10) (G + C)
  * **`G_fraction`** (Rank 13)
  * **`A_count`** (Rank 14), **`T_count`** (Rank 15), **`C_fraction`** (Rank 16), **`A_fraction`** (Rank 17), **`C_count`** (Rank 20), **`T_fraction`** (Rank 21), **`G_count`** (Rank 26), **`count_C`** (Rank 86), **`pyrimidine_fraction`** (Rank 69), **`fraction_C`** (Rank 100).

### Structural Composition
* **`max_homopolymer_run`** (Rank 9)
  * **What it represents:** The longest stretch of a single repeated nucleotide (e.g., "AAAA"). Long runs can cause off-target binding or manufacturing issues.
  * **How the code works:** Iterates through the sequence, keeping a running counter that increments when the current base matches the previous one, and resets otherwise.
* **`max_gc_stretch`** (Rank 19)
  * **What it represents:** The longest contiguous run of G or C bases, creating a highly rigid, strongly-binding segment.
  * **How the code works:** Iterates through the sequence, incrementing a counter whenever `c in "GC"`.

### Di- and Tri-nucleotide Frequencies
* **What they represent:** The prevalence of specific 2-mer and 3-mer motifs. Certain motifs (like CpG) are immunostimulatory, while others are recognized by specific nucleases or RNA binding proteins.
* **How they are calculated:**
  * Counts: Sliding window count of the motif.
  * Frequencies (for dinucleotides): Count divided by $(length - 1)$.
* **Code extraction:** Uses Python's `itertools.product` to generate all combinations and `string.count(motif)`.
* **Top Features in this sub-group:**
  * **Dinucleotide Frequencies:** `di_AC_freq` (Rank 22), `di_GT_freq` (Rank 23), `di_TC_freq`, `di_AG_freq`, `di_CT_freq`, `di_TT_freq`, `di_AT_freq`, `di_CA_freq`, `di_TA_freq`, `di_TG_freq`, `di_GA_freq`, `di_AA_freq`, `di_GC_freq`, `di_GG_freq`, `di_CC_freq`.
  * **Dinucleotide Counts:** `di_AG`, `di_TT`, `di_AC`, `di_CA`, `di_TG`, `di_GT`, `di_AT`, `di_TA`, `di_GG`, `di_GA`, `di_GC`, `di_CC`, `di_CT`.
  * **Trinucleotide Counts:** `tri_TCA`, `tri_ATT`, `tri_TGC`, `tri_ATG`, `tri_TTT`, `tri_AGG`, `tri_ATC`, `tri_ACA`, `tri_GCT`, `tri_GGG`, `tri_CTG`, `tri_AAT`, `tri_AGA`, `tri_AAG`, `tri_CCA`, `tri_AGT`, `tri_TAT`, `tri_GAC`, `tri_TTA`, `tri_CCC`, `tri_GGA`, `tri_CAC`, `tri_CAG`.

### Positional Anchors
* **What they represent:** The exact sequence context at the extreme 5' and 3' ends. The ends dictate nuclease resistance and initial binding kinetics to the target.
* **How they are calculated:** Binary flags checking if the start or end of the string matches a specific motif.
* **Code extraction:** Checking `seq[0]`, `seq[-1]`, `seq[:2]`, `seq[:3]`, etc.
* **Top Features in this sub-group:**
  * **Terminal Base:** `five_prime_base_C` (Rank 36), `five_prime_base_T` (Rank 37), `three_prime_base_A` (Rank 61).
  * **Terminal Motifs:** `first_3nt_GCT`, `first_2nt_GC`, `first_3nt_GTT`, `first_3nt_TTC`, `first_2nt_AT`, `last_3nt_ATT`, `first_2nt_TA`, `first_2nt_CC`.

## 5. Chemistry & Modifications (Gapmer Architecture)
*ASOs are heavily chemically modified. These features describe the "pattern" of modifications, typically following a gapmer design (modified ends with an unmodified DNA center).*

* **`count_d`** (Rank 11) & **`fraction_d`** (Rank 38)
  * **What it represents:** The number (and fraction) of unmodified deoxyribonucleotides (DNA). This usually constitutes the central "gap" that recruits RNase H.
  * **How the code works:** Counts occurrences of 'd' in the `chemical_pattern` string.
* **`count_M`** (Rank 12) & **`fraction_M`** (Rank 18)
  * **What it represents:** The amount of 2'-MOE modifications, which increase binding affinity and nuclease resistance.
  * **How the code works:** Counts occurrences of 'M' in the pattern string.
* **`n_chemistry_transitions`** (Rank 25)
  * **What it represents:** The number of times the chemical modification changes from one nucleotide to the next. In a standard gapmer, this is exactly 2 (Wing $\rightarrow$ Gap, and Gap $\rightarrow$ Wing).
  * **How to calculate:** Summing instances where `cp[j] != cp[j - 1]`.
* **`left_modified_wing_length`** (Rank 74)
  * **What it represents:** The size of the 5' modified "wing" of the gapmer.
  * **How the code works:** Iterates forward from the start of the `chemical_pattern`, counting until a 'd' (deoxy) is encountered.
* **`symmetry_flag`** (Rank 92)
  * **What it represents:** Indicates if the left wing and right wing are exactly the same length (e.g., a 5-10-5 gapmer is symmetric, 3-10-7 is asymmetric).
  * **How the code works:** `int(left_wing == right_wing)`

## 6. Backbone & Linkage
*Modifications to the internucleotide linkages.*

* **`contains_PO`** (Rank 76)
  * **What it represents:** Indicates if there is at least one phosphodiester (PO - the natural, unmodified) linkage. Most ASOs use phosphorothioate (PS) linkages for resistance, but strategic PO linkages can reduce toxicity.
  * **How the code works:** Checks if `"phosphodiester"` is a substring in the `linkage` text.
* **`predicted_PO_count`** (Rank 87)
  * **What it represents:** The estimated count of PO linkages.
  * **How the code works:** Parses the `linkage_location` string. It splits by `/` and counts all numeric integers, ignoring `"else"`.
