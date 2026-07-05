# metaConvert 1.1.0
- Added new effect size measures: risk difference, hazard ratio (user-input only), Cronbach's alpha, intraclass correlation coefficient, and single-group proportion
- Added a quality-flag system to summary() ('flags' argument) that detects inconsistent input data and implausible effect sizes
- Added a guidance system to summary() ('guidance' argument) that indicates which columns are missing when an effect size cannot be estimated
- Added new converters: es_from_spearman_rho(), es_from_cronbach_alpha(), es_from_icc(), es_from_prop_single_group(), es_from_rd_se(), es_from_rd_ci(), es_from_rd_pval()
- Added psychometric utilities: es_disattenuate(), compute_sem(), compute_sdc(), reliability_change_score()
- Added pool_arms() to pool the arms of multi-arm studies
- Added single-group pre/post converters and the within-group measures dw, gw and mdw
- Standard errors and confidence intervals are now provided for the NNT
- The Yates correction for chi-square inputs can now be set per row
- Improved the handling of paired/pre-post designs

# metaConvert 1.0.3
- Updated the citation
- Improved the checkings for the 95% CI asymmetry

# metaConvert 1.0.2
- Fixed a bug for the online app
- Updated the vignette
- Corrected typos in the documentation

# metaConvert 1.0.1
- Added an 'auto' argument in the hierarchy of convert_df()
- Improved the way Chi-sq and Phi are converted to other measures
- Improved the compareDF function to handle different rows ordering
- Improved the aggregate_df() function to handle different time-points
- Corrected typos in the documentation

# metaConvert 1.0.0
- First version released on CRAN
