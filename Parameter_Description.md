# Parameter Descriptions

Once you have installed the libraries, you can execute the following scripts on the provided Coffee dataset and you should be able to see the following views:

### 1. ungrid

<img src="images\documentation\1_ungrid.png" alt="ungrid" style="zoom: 85%;" />

The arguments for the `ungrid` command are:

1. **study**: the .sqlite file containing the study to be analyzed.
2. **PPM**: m/z tolerance in parts per million.
3. **RT_TOL**: retention time tolerance in minutes.
4. **MIN_SIGNAL**: minimum signal filter applied to all MS1 peaks.
5. **MIN_RANGE**: minimum ratio required between max and min signal within a feature. This is a form of minimum dynamic range requirement which eliminates "flat" features which do not resemble a valid peak. 

<div style="page-break-after: always;"></div>

### 2. dewey

<img src="images\documentation\2_dewey.png" alt="dewey" style="zoom:85%;" />

The arguments for the `dewey` command are:

1. **study**: the .sqlite file containing the study to be analyzed.
2. **Metlin**: a True/False variable indicating whether the Metlin library will be used.
3. **NIST**: a True/False variable indicating whether the NIST library will be used.
4. **LipidBLAST**: a True/False variable indicating whether the LipidBlast library will be used.
5. **decoy**: a True/False variable indicating whether to search the associated decoy libraries for every selected library.
6. **HITS**: number of top matches to keep per scan (these are ordered by NIST score).

<div style="page-break-after: always;"></div>

### 3. refine

<img src="images\documentation\3_refine.png" alt="refine" style="zoom:80%;" />

The arguments for the `refine` command are:

1. **unrefined**: the .features file containing the identification s to be refined.
2. **filter_by**: which dewey generated score should be used to refine the identification list. The options are: score, mcrl_score, dot, rev_dot, prob and percentile (the latter is used for ungrid-like, i.e. non-ID based features, where the strength of the feature is characterized by intensity percentile rather than a match score).
3. **equivalence_by**: defines equivalence for features being refined. The options are: inchikey, inchik, name and mzrt. InChIK refers to the first 14 characters of the InChIKey (the structural core) and mzrt is the option to use for quant-based features.
3. **FDR**: a True/False variable indicating whether or not to estimate a False Discovery Rate using the decoy hits present in the .feature file.
4. **MAX_FDR**: an integer number representing the maximum acceptable (estimated) False Discovery Rate. By default this is set to 100 meaning that all refined identifications will be conserved.
5. **IF_NO_FDR_WHAT_MIN_SCORE**: a numeric threshold used to filter out identifications that do no have a target score (specified by the `filter_by` option). This option is typically used when the FDR strategy is not available.
6. **RT_TOL**: retention time tolerance in minutes: two identically named entries which are within this range will be treated as redundant identifications of the same metabolite.
7. **MatchPolarity**: a True/False variable indicating whether or not feature identifications must match the polarity of the scan they supposedly identifying.
8. **Keep_Decoy_Hits**: a True/False variable indicating whether or not decoy hits should remain in the refined output list.

<div style="page-break-after: always;"></div>

### 4. skeleton

<img src="images\documentation\4_skeleton.png" alt="skeleton" style="zoom:95%;" />

The arguments for the `skeleton` command are:

1. **study**: the .sqlite file containing the study to be analyzed.
2. **features**: the .feature file containing the list of features to be quantified.
3. **MZ_TOLERANCE**: the m/z tolerance (in parts per million) of the final feature quantification.
4. **RT_TOLERANCE_**: the retention time tolerance (in minutes) of the final feature quantification.
5. **RT_WINDOW**: the window of retention time (in minutes) around the identified feature, that the algorithm is willing to scan in pursuit of the "anchor" identification (the most intense identification around which the final quantified feature will be defined).

<div style="page-break-after: always;"></div>

### 5. neutral_loss

<img src="images\documentation\5_neutral_loss.png" alt="neutral_loss" style="zoom:95%;" />

The arguments for the `neutral_loss` command are:

1. **study**: the .sqlite file containing the study to be analyzed.
2. **neutral_loss**: the target neutral loss.
3. **abs_mz_**: a True/False variable indicating whether or not the difference between the fragment and precursor is considered in absolute terms or whether the fragment must be strictly speaking lighter than the precursor by the specified neutral loss.
4. **PPM**: the m/z tolerance (in parts per million) of the neutral loss assignment.
5. **RT_TOL**: the retention time tolerance (in minutes) within which neutral loss occurrences are considered equivalent across rawfiles (i.e. across samples).

<div style="page-break-after: always;"></div>

### 6. precursor_ion

<img src="images\documentation\6_precursor_ion.png" alt="precursor" style="zoom:100%;" />

The arguments for the `precursor_ion` command are:

1. **study**: the .sqlite file containing the study to be analyzed.
2. **fragment**: the target fragment.
3. **abs_mz_**: a True/False variable indicating whether or not the difference between the fragment and precursor is considered in absolute terms or whether the fragment must be strictly speaking lighter than the precursor by the specified neutral loss.
4. **PPM**: the m/z tolerance (in parts per million) of the fragment identification.
5. **RT_TOL**: the retention time tolerance (in minutes) within which precursors yielding the target fragment are considered equivalent across rawfiles (i.e. across samples).

<div style="page-break-after: always;"></div>

### 7. timothee

<img src="images\documentation\7_timothee.png" alt="timothee" style="zoom:100%;" />

The arguments for the `timothee` command are:

1. **study**: the .sqlite file containing the study from which a spectral library will be extracted.
