#' Automatically compute effect sizes from a well formatted dataset
#'
#' @param x a well formatted dataset
#' @param measure the effect size measure that will be estimated from the information stored in the dataset. See details.
#' @param main_es a logical variable indicating whether a main effect size should be selected when overlapping data are present. See details.
#' @param split_adjusted a logical value indicating whether crude and adjusted effect sizes should be presented separately. See details.
#' @param format_adjusted presentation format of the adjusted effect sizes. See details.
#' @param hierarchy a character string indicating the hierarchy in the information to be prioritized for the effect size calculations. See details.
#' @param selection_auto a character string giving details on the best "auto" hierarchy to use (only useful when \code{es_selected="auto"} and \code{measure= "d", "g" or "md"}). See details.
#' @param es_selected the method used to select the main effect size when several information allows to estimate an effect size for the same association/comparison. Must be either "auto" (the default; the effect size computed from the information specified highest in the built-in hierarchy described below will be selected), "hierarchy" (same, using the hierarchy passed to the \code{hierarchy} argument), "minimum" (the smallest effect size will be selected) or "maximum" (the largest effect size will be selected). See details.
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from a 2x2
#'   contingency table. **Only \code{"tetrachoric"} is available, by design rather than
#'   as a temporary restriction.**
#'
#'   The tetrachoric correlation treats both binary variables as *dichotomised continua*
#'   and estimates the correlation of the two latent variables. It is the only member of
#'   this family that is invariant to the margins, which is what makes it poolable: the
#'   obvious alternative, the phi coefficient, has an attainable range bounded by the
#'   margins (with balanced exposure, max \eqn{|\phi|} is 1.00 at a 50% event rate but
#'   0.33 at 10% and 0.10 at 1%). Two studies of the same association but different event
#'   rates therefore report different phi values, so pooling phi manufactures
#'   heterogeneity that is pure margin artefact: at a fixed latent correlation of 0.40,
#'   \eqn{I^2} for pooled phi rises from 14% to 88% as the primary studies get *larger*,
#'   because the between-study variance is pinned by the margins while the within-study
#'   variance falls as \eqn{1/n}. Jacobs & Viechtbauer (2017) make the same point.
#'
#'   No argument is offered to select the alternative, because the choice cannot be made
#'   correctly: a 2x2 table with fixed *n* has three free parameters and the
#'   dichotomised-bivariate-normal family also has three, so the latent-normal model is
#'   saturated and admits no goodness-of-fit test. Nothing in the data can tell a user
#'   whether their variables are dichotomised or genuinely dichotomous.
#'
#'   **If your variables are genuinely dichotomous** (allocation, genotype, sex) rather
#'   than thresholded severity, a correlation is the wrong effect size for them: use a
#'   binary measure instead (\code{measure = "logor"}, \code{"rr"} or \code{"rd"}), all
#'   of which metaConvert estimates from the same table.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio.
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio.
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient.
#' @param pre_post_to_smd formula used to obtain a SMD from pre/post means and SD of two independent groups.
#' @param pool_sd a logical value indicating whether the standardizing SD should be pooled across the two
#'   groups (default \code{FALSE}). The two options target the same estimand when the arms' true SDs are
#'   equal (as randomization implies at baseline) and differ otherwise; the literature does not agree on
#'   which to prefer, so this is a deliberate choice and not a technical detail.
#'   \itemize{
#'     \item \code{FALSE} (default): each arm's change is standardized by that arm's own SD and the two
#'       within-group values are subtracted, their variances adding because the arms are independent. This
#'       is Morris's (2008) \eqn{d_{ppc1}}, from Becker (1988). It makes no assumption that the arms' true
#'       SDs are equal, and Viechtbauer (see the metafor-project Morris 2008 page) describes it as the more
#'       broadly applicable of the two.
#'     \item \code{TRUE}: the difference in mean change is divided by a single SD pooled across arms. This
#'       is Morris's (2008) \eqn{d_{ppc2}} (his eq. 8-9), which he recommends: it is more efficient, but it
#'       assumes the two arms' true standardizing SDs are equal.
#'   }
#' @param r_pre_post pre-post correlation across the two groups (use this argument only if the precise correlation in each group is unknown). Note that when this correlation is defaulted rather than reported, the pre/post standard errors and confidence intervals are unreliable: at a true correlation of 0.3 with the default 0.8, the reported pre/post variance is about 68% too small (95% CI coverage ~0.75). Supply \code{r_pre_post_exp}/\code{r_pre_post_nexp} when available, or run a sensitivity analysis over plausible values.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient
#'   correlation. Two things differ between the options, and they are independent.
#'   First, which correlation is estimated: "viechtbauer" (the default) estimates the
#'   biserial correlation and "lipsey_cooper" the point-biserial. Second, which
#'   transform lands in the \code{z} column when \code{measure = "z"}: "lipsey_cooper"
#'   reports Fisher's z, \code{atanh(r)}, whereas "viechtbauer" reports a
#'   variance-stabilising transform, a different function of the correlation. At a
#'   point-biserial \eqn{\rho = 0.75} it returns 1.0925 where \code{atanh()} of its own
#'   r is 1.7468.
#'   Every other family (\code{pearson_r}, \code{fisher_z}, the OR routes, 2x2) reports
#'   Fisher's z. So under the default a \code{measure = "z"} review holding both SMD
#'   studies and correlation studies pools two transforms; flag E8 reports it, and
#'   \code{smd_to_cor = "lipsey_cooper"} puts the whole pool on Fisher's z.
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param smd_denom standardizer for the standardized mean difference. "pooled" (default) uses the pooled endpoint SD (Cohen's d / Hedges' g); "glass" (alias "control") uses the control (non-experimental) endpoint SD (Glass's delta); "glass_robust" (alias "control_robust") is Glass's delta with a heteroscedasticity-consistent sampling variance. Only the endpoint means family (es_from_means_sd/se/ci) honours this argument: rows whose effect size comes from any other method (t/F, cohen_d/hedges_g, eta-squared, point-biserial r, medians/ranges, plots, ANCOVA, raw mean differences) always use the pooled-SD standardizer, and a message lists the scoping when a non-pooled value is requested ("glass_robust" additionally has a single variance form, so smd_var is ignored for it).
#' @param cor_to_smd formula used to convert a correlation coefficient value into a SMD.
#' @param prop_to_es method used to compute the effect size from the proportion. Must be either "raw", "logit" or "freeman_tukey" (see \code{\link{es_from_prop_single_group}}).
#' @param alpha_to_es method used to compute the effect size from Cronbach's alpha. One of "bonett" (default, the variance-stabilising log(1 - alpha) transform), "raw", or "hakstian_whalen" (the cube-root transform Rodriguez & Maeda (2006) recommend and most published reliability-generalization syntheses use; stored in metafor's increasing orientation, matching measure = "AHW"). See \code{\link{es_from_cronbach_alpha}}.
#' @param omega_to_es method used to compute the effect size from McDonald's omega. One of "bonett" (default), "raw" or "hakstian_whalen" (see \code{\link{es_from_omega}}). Unlike alpha, omega has no published sampling variance in (n, k), so by default the standard error is taken from the reported omega_se or confidence interval; see \code{omega_se_source} to borrow alpha's closed form.
#' @param alpha_se_source where the standard error of Cronbach's alpha comes from when the study reported neither an SE nor a confidence interval. "closed_form" (default, the pre-existing behaviour) uses the (n, k) formula belonging to the selected \code{alpha_to_es} scale; "reported" refuses it and leaves \code{alpha_se = NA}, so the row stays visible but is left out of the pool. The closed form inherits Bonett's multivariate-normality assumption. Measured on ordinary Likert items it is well calibrated (se_ratio 0.95-1.06, 95% coverage .94-.96), but under a severe floor effect it is ~22% too small (coverage .89) and on dichotomous items ~37% too small (coverage .78, and NOT improving with n). Flag V44 identifies the affected rows from \code{scale_mean} / \code{n_response_categories} / \code{scale_min}.
#' @param omega_se_source where the standard error of McDonald's omega comes from when the study reported neither an SE nor a confidence interval. "reported" (default) leaves \code{omega_se = NA}, so the row stays visible but unpooled. "closed_form" borrows alpha's (n, k) formula on the selected \code{omega_to_es} scale - which is what every published omega reliability generalisation does, and what metafor's \code{measure = "ABT"} computes, so it is required to reproduce that literature. It is opt-in because switching it on makes previously unpooled rows enter the pool. Simulation support (\code{simulations/studies/11_reliability_se.R}) is for one narrow claim - the closed form is no more wrong for omega than for alpha - and NOT for the broader claim that it is safe for reliability generalisation, where estimator mixture moves omega_h by far more than any variance question.
#' @param icc_to_es method used to compute the effect size from an ICC. Must be either "bonett" or "raw" (see \code{\link{es_from_icc}}).
#' @param icc_agreement_se what to do about the closed-form standard error of an absolute-agreement ICC (ICC(2,1)) - which is the default \code{icc_type}, and therefore also what an absent one resolves to. One of:
#' \itemize{
#'  \item \bold{"compute"} (default): emit it. This is the pre-existing behaviour and is accurate when the raters are exchangeable (negligible rater variance).
#'  \item \bold{"drop"}: return \code{NA} for it, so the row is left out of the pool by \code{summary()} instead of entering it over-weighted. A standard error the study itself reported (\code{icc_se}, or \code{icc_ci_lo}/\code{icc_ci_up}) is kept either way, since the option targets the approximation rather than the row.
#' }
#' The closed-form standard error is a one-way approximation assuming negligible between-rater variance. With real rater variance its measured 95% CI coverage is 0.82 at n = 20, 0.67 at n = 50, 0.32 at n = 200 and 0.14 at n = 1000: it degrades as studies get bigger, because ICC(2,1) inherits the between-rater mean square (k - 1 df) while the reported standard error shrinks like \eqn{1/\sqrt{n}}, so such a row can carry up to 50x too much weight. Agreement rows are flagged by V31 whichever option is used.
#' The default is \code{"compute"} in 2.1.0 for backward compatibility and is expected to become \code{"drop"} in a future release; set it explicitly if you care which you get.
#' Consistency ICCs (ICC(3,1)) are unaffected either way: for them the same formula is exact at leading order.
#' @param yates_chisq a logical value indicating whether the Chi square has been performed using Yates' correction for continuity. Can also be given as a column of the dataset when studies differ (rows left NA use this argument).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Use '"sd"' when the increase is expressed in standard deviations of the independent variable, and '"raw_scale"' when it is in the raw units of that variable. '"value"' and '"raw_data"' are accepted synonyms of '"raw_scale"'. Defaults to '"raw_scale"'. Read only by \code{cor_to_smd = "mathur"}; any other value is rejected rather than silently treated as raw units.
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds. Asymmetric CIs are only flagged, not modified (see \code{\link{summary.metaConvert}}).
#' @param verbose a logical variable indicating whether text outputs and messages should be generated. We recommend turning this option to FALSE only after having carefully read all the generated messages.
#' @param correct_inputs a logical value indicating whether invalid input values (negative SDs, inverted CI bounds, etc.) should be set to NA before the calculations (default TRUE). If FALSE, invalid values are only flagged in the 'flags' column of summary() and kept in the data.
#' @param flag_options a named list of thresholds used by the quality flags (see \code{\link{summary.metaConvert}} for the list of options and defaults).
#'
#' @details
#' This function automatically computes or converts between 21 effect sizes
#' measures from any relevant type of input data stored in the
#' dataset you pass to this function.
#'
#' ## Input validation
#' Before the calculations, the input data are checked (negative sample
#' sizes/SD/SE/p-values, inverted CI bounds, values outside their CI,
#' asymmetric CI). By default invalid values are set to NA; issues appear
#' in the 'flags' column of summary().
#'
#' ## Effect size measures
#' Possible effect size measures are:
#' 1. Cohen's d ("d")
#' 2. Hedges' g ("g")
#' 3. mean difference ("md")
#' 4. within-group Cohen's d ("dw")
#' 5. within-group Hedges' g ("gw")
#' 6. within-group mean difference ("mdw")
#' 7. (log) odds ratio ("or" and "logor")
#' 8. (log) risk ratio ("rr" and "logrr")
#' 9. (log) incidence rate ratio ("irr" and "logirr")
#' 10. (log) hazard ratio ("hr" and "loghr")
#' 11. correlation coefficient ("r")
#' 12. transformed r-to-z correlation coefficient ("z")
#' 13. partial correlation coefficient ("rp")
#' 14. Fisher's z of partial correlation ("zp")
#' 15. log variability ratio ("logvr")
#' 16. log coefficient of variation ("logcvr")
#' 17. number needed to treat ("nnt")
#' 18. risk difference ("rd")
#' 19. proportion ("prop")
#' 20. Cronbach's alpha ("alpha")
#'
#' 21. McDonald's omega ("omega")
#' 21. intraclass correlation coefficient ("icc")
#'
#' The hazard ratio is estimated from a user-reported effect size only
#' (\code{user_es_crude}/\code{user_es_adj}), since its computation requires
#' individual time-to-event data that a wide-format dataset cannot store.
#'
#' ## Computation of a main effect size
#' If you enter multiple types of input data
#' (e.g., means/sd of two groups and a student t-test value)
#' for the same comparison i.e., for the same row of the dataset,
#' the \code{convert_df()} function can have two behaviours.
#' If you set:
#' - \code{main_es = FALSE} the function will estimate all possible effect sizes from all
#' types of input data (which implies that if a comparison has **several types of input data**,
#' it will result in **multiple rows** in the dataframe returned by the function)
#' - \code{main_es = TRUE} the function will select one effect size per comparison
#' (which implies that if a comparison has **several types of input data**,
#' it will result in a **unique row** in the dataframe returned by the function)
#'
#' ## Selection of input data for the computation of the main effect size
#' If you choose to estimate one main effect size (i.e., by setting \code{main_es = TRUE}),
#' you have several options to select this main effect size.
#' If you set:
#' - \code{es_selected = "auto"}: the main effect size will be **automatically** selected, by prioritizing
#' specific types of input data over other (see next section "Hierarchy").
#' - \code{es_selected = "hierarchy"}: the main effect size will be selected, by prioritizing
#' specific types of input data over other (see next section "Hierarchy").
#' - \code{es_selected = "minimum"}: the main effect size will be selected, by selecting
#' the lowest effect size available.
#' - \code{es_selected = "maximum"}: the main effect size will be selected, by selecting
#' the highest effect size available.
#'
#' ## Hierarchy
#' More than 70 different combinations of input data can be used to estimate an effect size.
#' You can retrieve the effect size measures estimated by each combination of input data
#' in the \code{\link{see_input_data}()} function and online \code{https://metaconvert.org/input.html}.
#'
#' You have two options to use a hierarchy in the types of input data.
#' * an automatic way (\code{es_selected = "auto"})
#' * an manual way (\code{es_selected = "hierarchy"})
#'
#' ### Automatic
#' If you select an automatic hierarchy, here are the types of input data that will be prioritized.
#'
#' #### Crude SMD or MD (\code{measure=c("d", "g", "md")} and \code{selection_auto="crude"})
#' 1. User's input effect size value
#' 2. SMD value
#' 3. Means at post-test
#' 4. ANOVA/Student's t-test/point biserial correlation statistics
#' 5. Linear regression estimates
#' 6. Mean difference values
#' 7. Quartiles/median/maximum values
#' 8. Post-test means extracted from a plot
#' 9. Pre-test+post-test means or mean change
#' 10. Paired ANOVA/t-test statistics
#' 11. Odds ratio value
#' 12. Contingency table
#' 13. Correlation coefficients
#' 14. Phi/chi-square value
#'
#' #### Paired SMD or MD (\code{measure=c("d", "g", "md")} and \code{selection_auto="paired"})
#' 1. User's input effect size value
#' 2. Paired SMD value
#' 3. Pre-test+post-test means or mean change
#' 4. Paired ANOVA/t-test statistics
#' 5. Means at post-test
#' 6. ANOVA/Student's t-test/point biserial correlation
#' 7. Linear regression estimates
#' 8. Mean difference values
#' 9. Quartiles/median/maximum values
#' 10. Odds ratio value
#' 11. Contingency table
#' 12. Correlation coefficients
#' 13. Phi/chi-square value
#'
#' #### Adjusted SMD or MD (\code{measure=c("d", "g", "md")} and \code{selection_auto="adjusted"})
#' 1. User's input adjusted effect size value
#' 2. Adjusted SMD value
#' 3. Estimated marginal means from ANCOVA, with a standard deviation
#' 4. Adjusted mean difference from ANCOVA, with the residual standard deviation
#' 5. Estimated marginal means from ANCOVA, with a standard error or a confidence interval
#' 6. F- or t-test value (or partial eta-squared) from ANCOVA
#' 7. Adjusted mean difference from ANCOVA, with a standard error, a confidence interval or a p-value
#' 8. Estimated marginal means from ANCOVA extracted from a plot
#'
#' Items 3 to 7 are ordered by how much the covariate imbalance attenuates the estimate
#' (Lai & Kelley 2012). A route handed a standardizer directly (3, 4) is unaffected; one
#' that must recover it from a reported standard error, confidence interval or test
#' statistic sets the leverage term to zero and shrinks both the effect size and its
#' standard error. Item 5 carries half that term, items 6 and 7 the whole of it.
#'
#' #### Odds Ratio (\code{measure=c("or")})
#' 1.	User's input effect size value
#' 2.	Odds ratio value
#' 3.	Contingency table
#' 4.	Risk ratio values
#' 5.	Phi/chi-square value
#' 6.	Correlation coefficients
#' 7.	(Then hierarchy as for "d" or "g" option crude)
#'
#' #### Risk Ratio (\code{measure=c("rr")})
#' 1.	User's input effect size value
#' 2.	Risk ratio values
#' 3.	Contingency table
#' 4.	Odds ratio values
#' 5.	Phi/chi-square value
#'
#' #### Incidence rate ratio (\code{measure=c("irr")})
#' 1.	User's input effect size value
#' 2.	Number of cases and time of disease free observation time
#'
#' #### Correlation (\code{measure=c("r", "z")})
#' 1. User's input effect size value
#' 2. Correlation coefficients
#' 3. Contingency table
#' 4. Odds ratio value
#' 5. Phi/chi-square value
#' 6. SMD value
#' 7. Means at post-test
#' 8. ANOVA/Student's t-test/point biserial correlation
#' 9. Linear regression estimates
#' 10. Mean difference values
#' 11. Quartiles/median/maximum values
#' 12. Post-test means extracted from a plot
#' 13. Pre-test+post-test means or mean change
#' 14. Paired ANOVA/t-test
#'
#' #### Variability ratios (\code{measure=c("vr", "cvr")})
#' 1.	User's input effect size value
#' 2. means/variability indices at post-test
#' 3. means/variability indices at post-test extracted from a plot
#'
#' #### Number needed to treat (\code{measure=c("nnt")})
#' 1.	User's input effect size value
#' 2.	Contingency table
#' 3.	Odds ratio values
#' 4.	Risk ratio values
#' 5.	Incidence rate ratio (person-time NNT, requires baseline_rate)
#' 6.	Phi/chi-square value
#'
#' NNT values should not be pooled directly (the NNT confidence interval is
#' disjoint when the risk difference crosses zero). You should pool RD/OR/RR
#' values and convert the pooled estimate to NNT (Deeks, 2002; Cochrane
#' Handbook, Chapter 15).
#'
#' #### Risk difference (\code{measure=c("rd")})
#' 1.	User's input effect size value
#' 2.	Contingency table
#' 3.	Odds ratio values
#' 4.	Risk ratio values
#' 5.	Incidence rate ratio (requires baseline_rate)
#' 6.	Phi/chi-square value
#'
#' ### Manual
#'
#' If you select a manual hierarchy, you can specify the order in which you want to
#' use each type of input data.
#' You can prioritize some types of input data by placing them at the begining of the
#' hierarchy argument, and you must separate all input data with a ">" separator.
#' For example, if you set:
#' - \code{hierarchy = "means_sd > means_se > student_t"}, the convert_df function will prioritize
#' the means + SD, then the means + SE, then the Student's t-test to estimate the main effect
#' size.
#' - \code{hierarchy = "2x2 > or_se > phi"}, the convert_df function will prioritize
#' the contingency table, then the odds ratio value + SE, then the phi coefficient to estimate
#' the main effect size.
#'
#' Importantly, if none of the types of input data indicated in the \code{hierarchy} argument
#' can be used to estimate the target effect size measure,
#' the \code{convert_df()} function will automatically try to use other types of input
#' data to estimate an effect size.
#'
#' ## Adjusted effect sizes
#' Some datasets will be composed of crude (i.e., non-adjusted) types of input data
#' (such as standard means + SD, Student's t-test, etc.) and adjusted types of input data
#' (such as means + SE from an ANCOVA model, a t-test from an ANCOVA, etc.).
#'
#' In these situations, you can decide to:
#' - treat crude and adjusted input data the same way \code{split_adjusted = FALSE}
#' - split calculations for crude and adjusted types of input data \code{split_adjusted = TRUE}
#'
#' If you want to split the calculations, you can decide to present the final dataset:
#' - in a long format (i.e., crude and adjusted effect sizes presented in separate rows \code{format_adjusted = "long"})
#' - in a wide format (i.e., crude and adjusted effect sizes presented in separate columns \code{format_adjusted = "wide"})
#'
#' @return
#' The \code{convert_df()} function returns a list of
#' more than 70 dataframes
#' (one for each function automatically applied to the dataset).
#' These dataframes systematically contain the columns described in
#' \code{\link{metaConvert-package}}.
#' The list of dataframes can be easily converted to a single,
#' calculations-ready dataframe
#' using the summary function (see \code{\link{summary.metaConvert}}).
#'
#' @md
#'
#' @export convert_df
#'
#' @examples
#' res <- convert_df(df.haza,
#'   measure = "g",
#'   split_adjusted = TRUE,
#'   es_selected = "minimum",
#'   format_adjusted = "long"
#' )
#' summary(res)
convert_df <- function(x, measure = c("d", "g", "md", "dw", "gw", "mdw",
                                      "logor", "logrr", "logirr", "loghr",
                                      "nnt", "rd", "r", "z", "rp", "zp",
                                      "logvr", "logcvr", "prop", "alpha", "omega", "icc"),
                       main_es = TRUE,
                       es_selected = c("auto", "hierarchy", "minimum", "maximum"),
                       selection_auto = c("crude", "paired", "adjusted"),
                       split_adjusted = TRUE,
                       format_adjusted = c("wide", "long"),
                       verbose = TRUE,
                       max_asymmetry = 50,
                       hierarchy = "means_sd > means_se > means_ci",
                       table_2x2_to_cor = "tetrachoric",
                       rr_to_or = "metaumbrella",
                       or_to_rr = "metaumbrella_cases",
                       or_to_cor = "bonett",
                       smd_to_cor = "viechtbauer",
                       smd_var = "borenstein",
                       smd_denom = "pooled",
                       pre_post_to_smd = "bonett",
                       r_pre_post = 0.8,
                       cor_to_smd = "viechtbauer",
                       unit_type = "raw_scale",
                       yates_chisq = FALSE,
                       pool_sd = FALSE,
                       prop_to_es = "raw",
                       alpha_to_es = "bonett",
                       omega_to_es = "bonett",
                       alpha_se_source = "closed_form",
                       omega_se_source = "reported",
                       icc_to_es = "bonett",
                       icc_agreement_se = "compute",
                       correct_inputs = TRUE,
                       flag_options = list()) {
  # @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table.
  # table_2x2_to_cor = "tetrachoric",
  # x = dat
  # measure = "nnt";
  # main_es = TRUE;
  # es_selected = "hierarchy";
  # split_adjusted = TRUE;
  # format_adjusted = "wide";
  # verbose = TRUE;
  # hierarchy = "means_sd > means_se > means_ci";
  # rr_to_or = "metaumbrella";
  # or_to_rr = "dipietrantonj";
  # or_to_cor = "bonett";
  # table_2x2_to_cor = "lipsey";
  # smd_to_cor = "viechtbauer";
  # chisq_to_cor = "tetrachoric";
  # phi_to_cor = "tetrachoric";
  # pre_post_to_smd = "bonett";
  # cor_to_smd = "viechtbauer";
  # unit_type = "raw_scale";
  # yates_chisq = FALSE
  # max_asymmetry = 5
  # selection_auto = c("crude", "paired", "adjusted")
  measure = measure[1]
  es_selected = es_selected[1]
  selection_auto = selection_auto[1]
  format_adjusted = format_adjusted[1]

  x <- .check_data(x, split_adjusted = split_adjusted,
                   main_es = main_es, format = format_adjusted)
  # rows with defaulted r_pre_post
  .r_defaulted_exp <- is.na(x[, "r_pre_post_exp"])
  .r_defaulted_nexp <- is.na(x[, "r_pre_post_nexp"])
  x[.r_defaulted_exp, "r_pre_post_exp"] <- r_pre_post
  x[.r_defaulted_nexp, "r_pre_post_nexp"] <- r_pre_post

  # An unreported r that is imputed without comment is a substantive assumption. How
  # it bites depends on the standardizer: for change-SD methods (morris_dz,
  # cooper/morris_drm) the assumed r scales the point estimate, while for baseline-SD
  # (bonett) and average-SD (morris_dav) the point estimate is r-free and only the SE
  # and CI move. Warn once, and only for rows that actually feed an r-consuming route
  # (pre/post, mean-change or paired), since datasets with no such data never reach
  # the formulas.
  #
  # Which rows carry r-consuming data IN THAT ARM, over .r_consuming_columns().
  # Each mask must be ANDed with its own arm's .r_defaulted_*, never
  # ORed blindly: on a single-group row -- and on any two-group row whose pre/post,
  # mean-change or paired data sits in one arm only -- the absent arm's r_pre_post is
  # structurally NA, so the plain OR reports an imputation on a row that supplied its
  # correlation and whose arithmetic used the supplied value. Defined once and called at
  # both sites, so the verbose note and flag V6 cannot drift apart.
  .r_consuming_arms <- function(x) {
    cols <- intersect(.r_consuming_columns(), colnames(x))
    arm <- function(pattern) {
      cc <- grep(pattern, cols, value = TRUE)
      if (length(cc) == 0) return(rep(FALSE, nrow(x)))
      rowSums(!is.na(x[, cc, drop = FALSE])) > 0
    }
    list(exp = arm("_exp$"), nexp = arm("_nexp$"))
  }

  if (verbose) {
    # Single source of truth for "does this row feed an r-consuming route?", shared
    # with flag V6 below. Evaluated here against the input as the user supplied it;
    # V6 re-evaluates it after validation, which may have NA'd an invalid value.
    .arms <- .r_consuming_arms(x)
    .imputed_rows <- (.r_defaulted_exp & .arms$exp) | (.r_defaulted_nexp & .arms$nexp)
    .n_imputed <- sum(.imputed_rows)
    if (.n_imputed > 0) {
      # Change-SD standardizers make the point estimate itself depend on r; the
      # mean-change and paired-t/F routes are always standardized by the change SD
      # (coerced to cooper), so their point estimates depend on r regardless of the
      # requested method. Resolved per row: pre_post_to_smd is one of the method
      # arguments convert_df reads from a column of the dataset, and that merge (a few
      # lines below) has not run yet, so the scalar alone would name a standardizer the
      # affected rows never received.
      .pp_row <- if ("pre_post_to_smd" %in% colnames(x)) {
        .m <- as.character(x[, "pre_post_to_smd"])
        .m[is.na(.m)] <- as.character(pre_post_to_smd)[1]
        .m
      } else {
        rep(as.character(pre_post_to_smd)[1], nrow(x))
      }
      .pp_imputed <- sort(unique(.pp_row[.imputed_rows]))
      .pp_label <- paste0("'", .pp_imputed, "'", collapse = "/")
      .r_scales_point <- all(.pp_imputed %in% c("morris_dz", "morris_drm", "cooper"))
      .impact <- if (.r_scales_point) {
        paste0("Under the ", .pp_label, " standardizer the assumed\n",
               "  correlation scales the POINT estimate, not only the SE.")
      } else {
        paste0("Under the ", .pp_label, " standardizer the point estimate is\n",
               "  r-free and only the SE/CI depend on the assumed correlation; note that\n",
               "  mean-change and paired-t/F rows are standardized by the change SD, whose\n",
               "  POINT estimate does depend on r.")
      }
      message(
        "Note: r_pre_post was not reported for ", .n_imputed, " row(s) with pre/post, ",
        "mean-change or paired data;\n  the default r_pre_post = ", r_pre_post,
        " was assumed. ", .impact,
        "\n  Supply r_pre_post_exp/r_pre_post_nexp when available, or run a ",
        "sensitivity analysis\n  over plausible values."
      )
    }
  }

  # unit_type is read in exactly one place, as `unit_type == "sd"`, so every other
  # value selects raw units, a typo included. The argument is rejected outright; a
  # column only warns, because one bad cell must never abort a whole run, which is the
  # convention the dipietrantonj and proportion routes already follow.
  # See R/internal_guards.R.
  .validate_unit_type(unit_type)

  r_pre_post = rep(r_pre_post, nrow(x))
  for (i in c("rr_to_or",
              "or_to_rr",
              "or_to_cor",
              # "table_2x2_to_cor",
              "smd_to_cor",
              "smd_var",
              "smd_denom",
              "pre_post_to_smd",
              "cor_to_smd",
              "unit_type")) {
    if (i %in% colnames(x)) {
      if (length(sapply(mget(i), function(x) x)) > 1) {
        stop("The length of the argument '", i, "' should be 1, or should be specified in a column of the dataset.")
      }
      x[, i][is.na(x[, i])] <- sapply(mget(i), function(x) x)
    } else {
      if (length(sapply(mget(i), function(x) x)) > 1) {
        stop("The length of the argument '", i, "' should be 1, or should be specified in a column of the dataset.")
      }
      x[, i] <- sapply(mget(i), function(x) x)
    }
  }

  # Per-row values, after the column has been merged with the argument default.
  # Unrecognised cells are neutralised to NA, not just reported: they are handed on
  # to the exported es_from_*() routes, whose own argument check would stop.
  x[, "unit_type"] <- .validate_unit_type(x[, "unit_type"], arg_name = "unit_type",
                                          column = TRUE)

  # P10: smd_denom is honoured only by the endpoint means family; rows reaching
  # the SMD via any other route keep the pooled-SD standardizer. Say so once,
  # or a "Glass" analysis silently mixes estimands across methods.
  .denom_used <- .normalize_smd_denom(x[, "smd_denom"])
  if (any(!is.na(.denom_used) & .denom_used != "pooled")) {
    message(
      "Note: smd_denom is honoured by the endpoint means family only ",
      "(es_from_means_sd, es_from_means_se, es_from_means_ci).\n  Rows whose effect size ",
      "comes from any other method (t/F statistics, cohen_d/hedges_g, eta-squared,\n  ",
      "point-biserial r, medians/ranges/quartiles, plots, ANCOVA, raw mean differences) ",
      "use the pooled-SD\n  standardizer. Check 'info_used' in summary() to see which ",
      "standardizer each row received."
    )
  }

  # input validation
  # Drop the entries it rejected, so the default really is what gets used.
  .bad_fo <- .validate_flag_options(flag_options, context = "convert_df")   # roadmap 2.5
  if (length(.bad_fo)) flag_options <- flag_options[setdiff(names(flag_options), .bad_fo)]
  enable_info <- if (!is.null(flag_options$enable_informational)) {
    flag_options$enable_informational
  } else {
    .default_flag_options()$enable_informational
  }
  sd_ratio_max_val <- if (!is.null(flag_options$sd_ratio_max)) {
    flag_options$sd_ratio_max
  } else {
    .default_flag_options()$sd_ratio_max
  }
  sd_ratio_bl_ep_min_val <- if (!is.null(flag_options$sd_ratio_bl_ep_min)) {
    flag_options$sd_ratio_bl_ep_min
  } else {
    .default_flag_options()$sd_ratio_bl_ep_min
  }
  baseline_imb_max_val <- if (!is.null(flag_options$baseline_imbalance_max)) {
    flag_options$baseline_imbalance_max
  } else {
    .default_flag_options()$baseline_imbalance_max
  }
  enable_cross_row_val <- if (!is.null(flag_options$enable_cross_row)) {
    flag_options$enable_cross_row
  } else {
    .default_flag_options()$enable_cross_row
  }
  templated_min_match_val <- if (!is.null(flag_options$templated_min_match)) {
    flag_options$templated_min_match
  } else {
    .default_flag_options()$templated_min_match
  }
  paired_as_indep_tol_val <- if (!is.null(flag_options$paired_as_indep_tol)) {
    flag_options$paired_as_indep_tol
  } else {
    .default_flag_options()$paired_as_indep_tol
  }
  margin_range_min_val <- if (!is.null(flag_options$margin_range_min)) {
    flag_options$margin_range_min
  } else {
    .default_flag_options()$margin_range_min
  }
  floor_pos_severe_val <- if (!is.null(flag_options$floor_position_severe)) {
    flag_options$floor_position_severe
  } else {
    .default_flag_options()$floor_position_severe
  }
  floor_pos_mild_val <- if (!is.null(flag_options$floor_position_mild)) {
    flag_options$floor_position_mild
  } else {
    .default_flag_options()$floor_position_mild
  }
  range_sd_lo_mult_val <- if (!is.null(flag_options$range_sd_lo_mult)) {
    flag_options$range_sd_lo_mult
  } else {
    .default_flag_options()$range_sd_lo_mult
  }
  range_sd_hi_mult_val <- if (!is.null(flag_options$range_sd_hi_mult)) {
    flag_options$range_sd_hi_mult
  } else {
    .default_flag_options()$range_sd_hi_mult
  }
  # Normalise the reliability estimand and provenance columns before validation and
  # before they are stored, so the flags, the route and the data.frame the user
  # meta-regresses on all refer to the same levels. Without this, summary() returned
  # the raw strings ("psych", "bifactor CFA") while V39 reported the normalised ones,
  # and rma(mods = ~ omega_estimator) aborted on singleton levels.
  # .normalise_omega_type() returns NA for NA (unlike .normalise_icc_type() below), so a
  # blank cell survives this step and V38's own rescue line can still read it as the
  # documented default rather than as "unspecified".
  if ("omega_type" %in% colnames(x))
    x$omega_type <- .normalise_omega_type(x$omega_type, warn = FALSE)
  if ("omega_estimator" %in% colnames(x))
    x$omega_estimator <- .normalise_omega_estimator(x$omega_estimator, warn = FALSE)
  # icc_type too, for the same reason and one more: V31 keys on this column, and
  # without normalisation here it would compare the raw cell against the literal
  # "agreement", so 8 of the 10 spellings that es_from_icc() resolves to agreement
  # would never fire the flag. NA is left as NA. .normalise_icc_type() would map it to
  # "agreement", but V31 and es_from_icc() already treat an absent value as the
  # default, and rewriting NA into the column would misreport what the user supplied.
  if ("icc_type" %in% colnames(x)) {
    keep_na <- is.na(x$icc_type)
    x$icc_type <- .normalise_icc_type(x$icc_type, warn = FALSE)
    x$icc_type[keep_na] <- NA
  }

  validation <- .validate_input_data(x, max_asymmetry = max_asymmetry, verbose = verbose,
                                      enable_informational = enable_info,
                                      sd_ratio_max = sd_ratio_max_val,
                                      sd_ratio_bl_ep_min = sd_ratio_bl_ep_min_val,
                                      baseline_imbalance_max = baseline_imb_max_val,
                                      enable_cross_row = enable_cross_row_val,
                                      templated_min_match = templated_min_match_val,
                                      measure = measure,
                                      alpha_to_es = alpha_to_es,
                                      icc_to_es = icc_to_es,
                                      omega_to_es = omega_to_es,
                                      paired_as_indep_tol = paired_as_indep_tol_val,
                                      floor_position_severe = floor_pos_severe_val,
                                      floor_position_mild = floor_pos_mild_val,
                                      margin_range_min = margin_range_min_val,
                                      range_sd_lo_mult = range_sd_lo_mult_val,
                                      range_sd_hi_mult = range_sd_hi_mult_val,
                                      # The Tier-1 cross-row checks (V35-V43) assert a
                                      # fact about "the rest of the pool", so they must
                                      # see the same strata the Tier-2 ones do; and the
                                      # SE-source switches decide which reliability rows
                                      # will have an SE at all, which V44 reports on.
                                      flag_group = flag_options$flag_group,
                                      alpha_se_source = alpha_se_source,
                                      omega_se_source = omega_se_source,
                                      correct_inputs = correct_inputs)
  x <- validation$data

  # V11's bounded-column check may have set an out-of-range r_pre_post to NA, AFTER the
  # default fill near the top of this function. Such a row is operationally a defaulted
  # row and must be treated as one, in both directions:
  #   * it is re-filled here, because otherwise it reaches the es_from_* route with
  #     r = NA and picks up that route's own formal default (0.8), silently overriding
  #     the r_pre_post the caller passed to convert_df(). Measured on one row:
  #     g = 0.0727 under the route default against 0.2887 for the requested r = 0.3.
  #   * it is folded into .r_defaulted_*, because those masks were computed before the
  #     invalidation, so the row would otherwise carry r_defaulted = FALSE and get
  #     neither the V6 note below nor the "(r-sensitive: ...)" annotation summary()
  #     attaches from attr(res, "r_defaulted").
  # Inert when correct_inputs = FALSE, which preserves the value instead of NA-ing it.
  # NOTE: r_pre_post was expanded to one value per row above, so it is INDEXED here,
  # not recycled.
  .r_voided_exp  <- is.na(x[, "r_pre_post_exp"])  & !.r_defaulted_exp
  .r_voided_nexp <- is.na(x[, "r_pre_post_nexp"]) & !.r_defaulted_nexp
  if (any(.r_voided_exp))  x[.r_voided_exp,  "r_pre_post_exp"]  <- r_pre_post[.r_voided_exp]
  if (any(.r_voided_nexp)) x[.r_voided_nexp, "r_pre_post_nexp"] <- r_pre_post[.r_voided_nexp]
  .r_defaulted_exp  <- .r_defaulted_exp  | .r_voided_exp
  .r_defaulted_nexp <- .r_defaulted_nexp | .r_voided_nexp

  # V6: flag rows where r_pre_post used the default value and pre-post data is
  # present. Reuses .r_consuming_arms(), defined once above over
  # .r_consuming_columns(). Keeping one list matters here because intersect() drops
  # any name that is not a column, so a list naming columns that do not exist would
  # leave paired-t/F rows -- the ones where the assumed r scales the point estimate --
  # without a V6 flag, without the r_defaulted attribute, and without the
  # (r-sensitive: ...) annotations in summary().
  .v6_arms <- .r_consuming_arms(x)
  .r_defaulted <- (.r_defaulted_exp & .v6_arms$exp) | (.r_defaulted_nexp & .v6_arms$nexp)
  for (i in which(.r_defaulted)) {
    # The quoted column name is load-bearing, not decoration: .v_flag_matches_scope()
    # routes a Tier-1 message on the FIRST quoted token it finds, and returns TRUE for
    # both scopes when it finds none -- so without it this crude-scope note is copied
    # into flags_adjusted, where the covariate correlation is cov_outcome_r and
    # r_pre_post plays no part.
    # Name the arm(s) actually defaulted, not always the exposed one: a row that
    # supplied r_pre_post_exp and omitted r_pre_post_nexp was told 'r_pre_post_exp'
    # was not provided, which is false about the user's own data. The crude-scope
    # column stays FIRST either way, so .v_flag_matches_scope() still routes to
    # flags_crude.
    .def_cols <- c(if (.r_defaulted_exp[i]) "r_pre_post_exp",
                   if (.r_defaulted_nexp[i]) "r_pre_post_nexp")
    if (!length(.def_cols)) .def_cols <- "r_pre_post_exp"
    msg <- sprintf("[INFO] Default r_pre_post = %s used (%s not provided by user)",
                   r_pre_post[1],
                   paste0("'", .def_cols, "'", collapse = " and "))
    if (nzchar(validation$issues[i])) {
      validation$issues[i] <- paste(validation$issues[i], msg, sep = "; ")
    } else {
      validation$issues[i] <- msg
    }
  }

  if (verbose) message("Calculations in progress, it may take up to 30 sec...")
  measure = tolower(measure)
  if (!measure %in% c("d", "g", "md", "dw", "gw", "mdw", "r", "z", "rp", "zp", "or", "rr", "irr", "hr", "logor", "logrr", "logirr", "loghr", "logvr", "logcvr", "nnt", "rd", "prop", "alpha", "omega", "icc")) {
    stop(paste0("'", measure, "' not in tolerated measures. Possible inputs are: 'md', 'd', 'g', 'dw', 'gw', 'mdw', 'or', 'rr', 'irr', 'hr', 'logor', 'logrr', 'logirr', 'loghr', 'r', 'z', 'rp', 'zp', 'logvr', 'logcvr', 'nnt', 'rd', 'prop', 'alpha', 'omega', 'icc'"))
  } else if (!split_adjusted %in% c(TRUE, FALSE)) {
    stop(paste0("'", split_adjusted, "' not in tolerated values for the 'split_adjusted' argument. Should be a logical value (TRUE/FALSE)"))
  } else if (!format_adjusted %in% c("wide", "long")) {
    stop(paste0("'", format_adjusted, "' not in tolerated values for the 'format_adjusted' argument. Should be either 'long' or 'wide'."))
  } else if (!es_selected %in% c("auto", "minimum", "maximum", "hierarchy")) {
    stop(paste0("'", es_selected, "' not in tolerated values for the 'es_selected' argument. Possible inputs are: 'auto', 'hierarchy', 'minimum', 'maximum'"))
  } else if (!main_es %in% c(TRUE, FALSE)) {
    stop(paste0("'", main_es, "' not in tolerated values for the 'main_es' argument. Should be a logical value (TRUE/FALSE)"))
  }

  # 'es_selected' picks one effect size per comparison, whereas 'main_es = FALSE'
  # returns every estimation route. Combining them without this guard overwrites each
  # route row with the same min/max value, leaving the row count inflated k-fold.
  if (!main_es && es_selected %in% c("minimum", "maximum")) {
    stop(paste0("'es_selected = \"", es_selected, "\"' selects a single effect size per comparison ",
                "and is incompatible with 'main_es = FALSE', which returns every estimation route. ",
                "Use es_selected = 'auto' or 'hierarchy' with main_es = FALSE."))
  }

  # table_2x2_to_cor is commented out of the per-row method-column loop and of all
  # three es_from_2x2*() call sites, so without a check here
  # convert_df(x, measure = "r", table_2x2_to_cor = "banana") would run without error
  # and return byte-identical output, leaving the guard to fire only on a direct
  # es_from_2x2() call. Validate it on this path too, so a typo is reported rather
  # than absorbed.
  if (!all(table_2x2_to_cor %in% c("tetrachoric"))) {
    stop(paste0(
      "'",
      paste(unique(table_2x2_to_cor[!table_2x2_to_cor %in% c("tetrachoric")]),
            collapse = "', '"),
      "' not in tolerated values for the 'table_2x2_to_cor' argument. ",
      "The only possible input is 'tetrachoric' -- see ?convert_df for why this is ",
      "by design rather than a temporary restriction."
    ), call. = FALSE)
  }

  if (measure %in% c("or", "rr", "irr", "hr")) {
    exp <- TRUE
    if (measure == "or") {
      measure <- "logor"
    } else if (measure == "rr") {
      measure <- "logrr"
    } else if (measure == "irr") {
      measure <- "logirr"
    } else if (measure == "hr") {
      measure <- "loghr"
    }
  } else {
    exp <- FALSE
  }

  # SMD ----------------------------------------------
  es_cohen_d <- with(x, es_from_cohen_d(
    cohen_d = cohen_d, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_d = reverse_d)
  )

  es_cohen_d_adj <- with(x, es_from_cohen_d_adj(
    cohen_d_adj = cohen_d_adj, n_cov_ancova = n_cov_ancova, cov_outcome_r = cov_outcome_r,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_d = reverse_d
  ))
  es_hedges_g <- with(x, es_from_hedges_g(
    hedges_g = hedges_g, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_g = reverse_g
  ))


  # OR ----------------------------------------------
  # es_from_or() is the OR-alone route: it has no reported uncertainty to work from,
  # so it imputes var(logOR) by enumerating every 2x2 compatible with the case/control
  # margins and averaging (.se_from_or, R/internal_multiple_formulas.R:399). That
  # average runs about 1.4x wide, and its inflation grows with study size.
  #
  # On a row that also reports an SE or a CI, this route is therefore a dominated
  # duplicate of es_from_or_se()/es_from_or_ci(): identical point estimate, worse
  # variance. Left in, it inflates n_estimations, becomes eligible for selection under
  # es_selected = "minimum"/"maximum", and enters the Category-E cross-method
  # dispersion and CI-overlap checks as a spurious second method. Suppress it there.
  #
  # p-values are excluded from the suppression, even though or_list_L2 now ranks
  # es_odds_ratio_pval above es_odds_ratio: the hierarchy already prefers the exact
  # p-value SE, so the reconstruction is no longer selected on those rows, and leaving
  # it running keeps it available as a cross-method comparator (and as a fallback when
  # the p-value is unusable -- see the or_pval <= 0 guard below).
  #
  # The suppression must also never remove a route the user explicitly asked for. With
  # es_selected = "hierarchy" and "or" named in `hierarchy`, the marginal
  # reconstruction is the requested estimator rather than a redundant duplicate, and
  # dropping it would leave the row with no estimate at all. Token-matched on the
  # split hierarchy so that "or_se" and "or_ci" do not count as a request for "or".
  .or_route_requested <- "or" %in% trimws(strsplit(hierarchy, ">", fixed = TRUE)[[1]])
  .or_uncertainty_reported <- if (.or_route_requested) {
    rep(FALSE, nrow(x))
  } else {
    with(x,
      !is.na(logor_se) |
        (!is.na(or_ci_lo) & !is.na(or_ci_up)) |
        (!is.na(logor_ci_lo) & !is.na(logor_ci_up))
    )
  }
  es_odds_ratio <- with(x, es_from_or(
    or = ifelse(.or_uncertainty_reported, NA, or),
    logor = ifelse(.or_uncertainty_reported, NA, logor),
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample, small_margin_prop = small_margin_prop,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or = reverse_or
  ))
  es_odds_ratio_se <- with(x, es_from_or_se(
    or = or, logor = logor, logor_se = logor_se, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
    n_cases = n_cases, n_controls = n_controls,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor,
    reverse_or = reverse_or
  ))
  es_odds_ratio_ci <- with(x, es_from_or_ci(
    or = or, or_ci_lo, or_ci_up, logor = logor, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    logor_ci_lo = logor_ci_lo, logor_ci_up = logor_ci_up,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or = reverse_or,
    max_asymmetry = max_asymmetry
  ))
  # A p-value of 0 -- what "p < 0.001" becomes when it is transcribed as a number --
  # sends qnorm(p/2) to -Inf and the recovered SE to exactly 0, i.e. an infinite
  # inverse-variance weight. Its twin p = 1 is degenerate in the same way and in the
  # same route: qnorm(1/2) is 0, so the recovered SE is +Inf and the CI is (-Inf, Inf).
  # Both are neutralised here rather than in .positive_columns(), which sees 0 and 1 as
  # in range for a probability. Harmless while this route ranked below the marginal
  # reconstruction; it is now the selected one, so a degenerate p would take the row
  # while a perfectly good reconstruction SE sat unused beside it.
  .or_pval_usable <- with(x, ifelse(!is.na(or_pval) & (or_pval <= 0 | or_pval >= 1),
                                    NA_real_, or_pval))
  es_odds_ratio_pval <- with(x, es_from_or_pval(
    or = or, logor = logor, or_pval = .or_pval_usable, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or_pval = reverse_or_pval
  ))
  es_logreg_t <- with(x, es_from_logreg_t(
    or = or, logor = logor, rr = rr, logrr = logrr,
    logreg_t = logreg_t, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor,
    rr_to_or = rr_to_or,
    reverse_logreg_t = reverse_logreg_t
  ))
  # R ----------------------------------------------
  es_pearson_r <- with(x, es_from_pearson_r(
    pearson_r = pearson_r, n_sample = n_sample,
    unit_type = unit_type, sd_iv = sd_iv,
    unit_increase_iv = unit_increase_iv,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_pearson_r = reverse_pearson_r
  ))
  es_fisher_z <- with(x, es_from_fisher_z(
    fisher_z = fisher_z, n_sample = n_sample,
    unit_type = unit_type, sd_iv = sd_iv,
    unit_increase_iv = unit_increase_iv,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_fisher_z = reverse_fisher_z
  ))
  es_spearman_r <- with(x, es_from_spearman_rho(
    spearman_r = spearman_r, n_sample = n_sample,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv,
    unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_spearman_r = reverse_spearman_r
  ))

  # MEANS ----------------------------------------------
  es_means_sd_raw <- with(x, es_from_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp,
    smd_var = smd_var, smd_denom = smd_denom, smd_to_cor = smd_to_cor, reverse_means = reverse_means
  )
  )
  es_means_se_raw <- with(
    x,
    es_from_means_se(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_exp = mean_exp, mean_se_exp = mean_se_exp,
      mean_nexp = mean_nexp, mean_se_nexp = mean_se_nexp,
      smd_var = smd_var, smd_denom = smd_denom, smd_to_cor = smd_to_cor, reverse_means = reverse_means
    )
  )
  es_means_ci_raw <- with(
    x,
    es_from_means_ci(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_exp = mean_exp, mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
      mean_nexp = mean_nexp, mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
      smd_var = smd_var, smd_denom = smd_denom, smd_to_cor = smd_to_cor, reverse_means = reverse_means,
      max_asymmetry = max_asymmetry
    )
  )
  es_means_sd_pooled <- with(
    x,
    es_from_means_sd_pooled(
      n_exp = n_exp, n_nexp = n_nexp, mean_exp = mean_exp, mean_nexp = mean_nexp,
      mean_sd_pooled = mean_sd_pooled, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_means = reverse_means
    )
  )
  # plot
  es_plot_means_raw <- with(x, es_from_plot_means(
    n_exp = n_exp, n_nexp = n_nexp,
    plot_mean_exp = plot_mean_exp, plot_mean_nexp = plot_mean_nexp,
    plot_mean_sd_lo_exp = plot_mean_sd_lo_exp, plot_mean_sd_lo_nexp = plot_mean_sd_lo_nexp,
    plot_mean_sd_up_exp = plot_mean_sd_up_exp, plot_mean_sd_up_nexp = plot_mean_sd_up_nexp,
    plot_mean_se_lo_exp = plot_mean_se_lo_exp, plot_mean_se_lo_nexp = plot_mean_se_lo_nexp,
    plot_mean_se_up_exp = plot_mean_se_up_exp, plot_mean_se_up_nexp = plot_mean_se_up_nexp,
    plot_mean_ci_lo_exp = plot_mean_ci_lo_exp, plot_mean_ci_lo_nexp = plot_mean_ci_lo_nexp,
    plot_mean_ci_up_exp = plot_mean_ci_up_exp, plot_mean_ci_up_nexp = plot_mean_ci_up_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor,
    reverse_plot_means = reverse_plot_means
  ))
  # PRE POST MEANS ----------------------------------------------
  es_means_sd_pre_post <- with(
    x,
    es_from_means_sd_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_sd_exp = mean_pre_sd_exp, mean_sd_exp = mean_sd_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_sd_nexp = mean_pre_sd_nexp, mean_sd_nexp = mean_sd_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd,
      pool_sd = pool_sd
    )
  )
  es_means_se_pre_post <- with(
    x,
    es_from_means_se_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_se_exp = mean_pre_se_exp, mean_se_exp = mean_se_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_se_nexp = mean_pre_se_nexp, mean_se_nexp = mean_se_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd,
      pool_sd = pool_sd
    )
  )
  es_means_ci_pre_post <- with(
    x,
    es_from_means_ci_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_ci_lo_exp = mean_pre_ci_lo_exp, mean_pre_ci_up_exp = mean_pre_ci_up_exp,
      mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_ci_lo_nexp = mean_pre_ci_lo_nexp, mean_pre_ci_up_nexp = mean_pre_ci_up_nexp,
      mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd,
      pool_sd = pool_sd,
      max_asymmetry = max_asymmetry
    )
  )

  # bonett/morris_dav need separate pre/post SDs, use cooper for mean change and paired
  # t/f. Resolved PER ROW: pre_post_to_smd is one of the nine method arguments merged
  # with a per-row column of the same name above, and all 13 receiving routes accept one
  # method per row (es_from_paired_t -> .paired_t_to_smd, es_from_mean_change_* -> the
  # kernel), so restricting from the length-1 formal argument discarded the column for
  # the whole mean-change and paired-t/F family: d_z rows silently got d_rm, which at
  # r = 0.6 is 11% off in the point estimate and 13% in the SE.
  pre_post_to_smd_restricted <- ifelse(
    x[, "pre_post_to_smd"] %in% c("bonett", "morris_dav"),
    "cooper",
    x[, "pre_post_to_smd"]
  )

  # The coercion above changes the estimand for the affected rows without any other
  # signal: they are standardized by the change SD (cooper/d_rm) while means_sd_pre_post rows
  # keep the requested baseline-SD (bonett) or average-SD (morris_dav)
  # standardizer. Tell the user whenever data that actually routes through a
  # coerced method is present.
  # Tested per row now that both sides are length-nrow vectors: the scalar comparison
  # `!identical(vector, scalar)` would be TRUE on every verbose run, including under
  # morris_dz, where no coercion happens at all.
  .coerced_row <- !is.na(x[, "pre_post_to_smd"]) &
    pre_post_to_smd_restricted != x[, "pre_post_to_smd"]
  if (verbose && any(.coerced_row)) {
    .coerced_cols <- c(
      "mean_change_exp", "mean_change_nexp", "mean_change_sd_exp",
      "mean_change_sd_nexp", "mean_change_se_exp", "mean_change_se_nexp",
      "mean_change_pval_exp", "mean_change_pval_nexp",
      "paired_t_exp", "paired_t_nexp", "paired_f_exp", "paired_f_nexp",
      "paired_t_pval_exp", "paired_t_pval_nexp",
      "paired_f_pval_exp", "paired_f_pval_nexp"
    )
    .coerced_cols <- intersect(.coerced_cols, colnames(x))
    .has_coerced_data <- length(.coerced_cols) > 0 &&
      any(!is.na(x[.coerced_row, .coerced_cols, drop = FALSE]))
    if (.has_coerced_data) {
      # Name the values actually coerced, not the scalar argument: a sheet may request
      # a different standardizer on every row.
      .coerced_lbl <- paste0("'", sort(unique(x[.coerced_row, "pre_post_to_smd"])), "'",
                             collapse = "/")
      message(
        "Note: pre_post_to_smd = ", .coerced_lbl, " requires separate pre/post SDs, ",
        "which mean-change and paired t/F data do not carry.\n",
        "  Those rows use 'cooper' (morris_drm: the change SD, rescaled by sqrt(2(1 - r)) ",
        "onto the\n  raw-score metric) instead. Both standardizers are on the raw-score ",
        "scale, so the rows\n  remain combinable, but 'morris_drm' additionally assumes ",
        "equal pre/post SDs and depends\n  on r_pre_post. Supply pre/post SDs to use ",
        .coerced_lbl, " throughout."
      )
    }
  }

  es_means_change_sd <- with(x, es_from_mean_change_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_sd_exp = mean_change_sd_exp,
    mean_change_nexp = mean_change_nexp, mean_change_sd_nexp = mean_change_sd_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change,
    pre_post_to_smd = pre_post_to_smd_restricted,
    pool_sd = pool_sd
  ))
  es_mean_change_se <- with(x, es_from_mean_change_se(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_se_exp = mean_change_se_exp,
    mean_change_nexp = mean_change_nexp, mean_change_se_nexp = mean_change_se_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change,
    pre_post_to_smd = pre_post_to_smd_restricted,
    pool_sd = pool_sd
  ))
  es_mean_change_ci <- with(x, es_from_mean_change_ci(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp,
    mean_change_ci_lo_exp = mean_change_ci_lo_exp, mean_change_ci_up_exp = mean_change_ci_up_exp,
    mean_change_nexp = mean_change_nexp,
    mean_change_ci_lo_nexp = mean_change_ci_lo_nexp, mean_change_ci_up_nexp = mean_change_ci_up_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change,
    pre_post_to_smd = pre_post_to_smd_restricted,
    pool_sd = pool_sd,
    max_asymmetry = max_asymmetry
  ))
  es_mean_change_pval <- with(x, es_from_mean_change_pval(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_pval_exp = mean_change_pval_exp,
    mean_change_nexp = mean_change_nexp, mean_change_pval_nexp = mean_change_pval_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change,
    pre_post_to_smd = pre_post_to_smd_restricted,
    pool_sd = pool_sd
  ))

  # SINGLE-GROUP WITHIN-GROUP EFFECTS (experimental group only) --------------
  es_means_pp_sg <- with(x, es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_sd_exp = mean_pre_sd_exp, mean_sd_exp = mean_sd_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd, smd_to_cor = smd_to_cor,
    reverse_means_pre_post = reverse_means_pre_post
  ))

  es_means_se_pp_sg <- with(x, es_from_means_se_pre_post_single_group(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_se_exp = mean_pre_se_exp, mean_se_exp = mean_se_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd, smd_to_cor = smd_to_cor,
    reverse_means_pre_post = reverse_means_pre_post
  ))

  es_means_ci_pp_sg <- with(x, es_from_means_ci_pre_post_single_group(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_ci_lo_exp = mean_pre_ci_lo_exp, mean_pre_ci_up_exp = mean_pre_ci_up_exp,
    mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd, smd_to_cor = smd_to_cor,
    max_asymmetry = max_asymmetry,
    reverse_means_pre_post = reverse_means_pre_post
  ))

  es_mean_change_sg <- with(x, es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change_exp, mean_change_sd_exp = mean_change_sd_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd_restricted,
    reverse_mean_change = reverse_mean_change
  ))

  es_mean_change_se_sg <- with(x, es_from_mean_change_se_single_group(
    mean_change_exp = mean_change_exp, mean_change_se_exp = mean_change_se_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd_restricted,
    reverse_mean_change = reverse_mean_change
  ))

  es_mean_change_ci_sg <- with(x, es_from_mean_change_ci_single_group(
    mean_change_exp = mean_change_exp,
    mean_change_ci_lo_exp = mean_change_ci_lo_exp, mean_change_ci_up_exp = mean_change_ci_up_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor, max_asymmetry = max_asymmetry,
    pre_post_to_smd = pre_post_to_smd_restricted,
    reverse_mean_change = reverse_mean_change
  ))

  es_mean_change_pval_sg <- with(x, es_from_mean_change_pval_single_group(
    mean_change_exp = mean_change_exp, mean_change_pval_exp = mean_change_pval_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd_restricted,
    reverse_mean_change = reverse_mean_change
  ))

  es_paired_t_sg <- with(x, es_from_paired_t_single_group(
    paired_t_exp = paired_t_exp,
    n_exp = n_exp, r_pre_post_exp = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd_restricted, smd_to_cor = smd_to_cor,
    reverse_paired_t = reverse_paired_t
  ))

  # SINGLE-GROUP PROPORTIONS --------------
  es_prop_sg <- with(x, es_from_prop_single_group(
    prop = prop,
    n_sample = n_sample,
    prop_to_es = prop_to_es,
    reverse_prop = reverse_prop
  ))

  es_prop_counts_sg <- with(x, es_from_prop_single_group_counts(
    n_cases = n_cases,
    n_sample = n_sample,
    prop_to_es = prop_to_es,
    reverse_prop = reverse_prop
  ))

  # CRONBACH'S ALPHA -----------------------------------------
  es_omega_sg <- with(x, es_from_omega(
    omega = omega, omega_se = omega_se,
    omega_ci_lo = omega_ci_lo, omega_ci_up = omega_ci_up,
    n_sample = n_sample, n_items = n_items,
    omega_type = omega_type, omega_estimator = omega_estimator,
    omega_to_es = omega_to_es, omega_se_source = omega_se_source
  ))
  es_alpha_sg <- with(x, es_from_cronbach_alpha(
    cronbach_alpha = cronbach_alpha, n_sample = n_sample,
    n_items = n_items,
    cronbach_alpha_se = cronbach_alpha_se,
    cronbach_alpha_ci_lo = cronbach_alpha_ci_lo,
    cronbach_alpha_ci_up = cronbach_alpha_ci_up,
    alpha_to_es = alpha_to_es, alpha_se_source = alpha_se_source
  ))

  # ICC -------------------------------------------------------
  es_icc_sg <- with(x, es_from_icc(
    icc = icc, n_sample = n_sample,
    n_measurements = n_measurements, icc_type = icc_type,
    icc_se = icc_se, icc_ci_lo = icc_ci_lo, icc_ci_up = icc_ci_up,
    icc_to_es = icc_to_es, agreement_se = icc_agreement_se
  ))

  es_paired_t <- with(x, es_from_paired_t(
    paired_t_exp = paired_t_exp, paired_t_nexp = paired_t_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_t = reverse_paired_t,
    pre_post_to_smd = pre_post_to_smd_restricted
  ))

  es_paired_t_pval <- with(x, es_from_paired_t_pval(
    paired_t_pval_exp = paired_t_pval_exp,
    paired_t_pval_nexp = paired_t_pval_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_t_pval = reverse_paired_t_pval,
    reverse_paired_t_pval_exp = reverse_paired_t_pval_exp,
    reverse_paired_t_pval_nexp = reverse_paired_t_pval_nexp,
    pre_post_to_smd = pre_post_to_smd_restricted
  ))

  es_paired_f <- with(x, es_from_paired_f(paired_f_exp, paired_f_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_f = reverse_paired_f,
    reverse_paired_f_exp = reverse_paired_f_exp,
    reverse_paired_f_nexp = reverse_paired_f_nexp,
    pre_post_to_smd = pre_post_to_smd_restricted
  ))

  es_paired_f_pval <- with(x, es_from_paired_f_pval(
    paired_f_pval_exp = paired_f_pval_exp,
    paired_f_pval_nexp = paired_f_pval_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_f_pval = reverse_paired_f_pval,
    reverse_paired_f_pval_exp = reverse_paired_f_pval_exp,
    reverse_paired_f_pval_nexp = reverse_paired_f_pval_nexp,
    pre_post_to_smd = pre_post_to_smd_restricted
  ))
  # ANOVA, Student t-test  ----------------------------------------------
  es_t_student <- with(x, es_from_student_t(
    student_t = student_t, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_student_t = reverse_student_t
  ))
  es_t_student_pval <- with(x, es_from_student_t_pval(
    student_t_pval = student_t_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_student_t_pval = reverse_student_t_pval
  ))
  es_anova_f <- with(x, es_from_anova_f(
    anova_f = anova_f, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_anova_f = reverse_anova_f
  ))
  es_anova_f_pval <- with(x, es_from_anova_pval(
    anova_f_pval = anova_f_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_anova_f_pval = reverse_anova_f_pval
  ))
  es_etasq <- with(x, es_from_etasq(
    etasq = etasq, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_etasq = reverse_etasq
  ))
  # etasq_adj gets its own reverse column, like the other adjusted routes
  # (reverse_ancova_means / reverse_ancova_md). It is independent of reverse_etasq and,
  # like every other reverse_* column, defaults to FALSE when absent or NA.
  es_etasq_adj <- with(x, es_from_etasq_adj(
    etasq_adj = etasq_adj, n_exp = n_exp, n_nexp = n_nexp, n_cov_ancova = n_cov_ancova,
    cov_outcome_r = cov_outcome_r, smd_var = smd_var, smd_to_cor = smd_to_cor,
    reverse_etasq = reverse_etasq_adj
  ))

  # MD   ----------------------------------------------
  es_md_sd <- with(x, es_from_md_sd(
    md = md, md_sd = md_sd, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_md = reverse_md
  ))
  es_md_se <- with(x, es_from_md_se(
    md = md, md_se = md_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_md = reverse_md
  ))
  es_md_ci <- with(x, es_from_md_ci(
    md = md, md_ci_lo = md_ci_lo, md_ci_up = md_ci_up, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_md = reverse_md,
    max_asymmetry = max_asymmetry
  ))
  es_md_pval <- with(x, es_from_md_pval(
    md = md, md_pval = md_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_md = reverse_md
  ))

  # RANGE/QUARTILES  ----------------------------------------------
  es_med_quarts <- with(x, es_from_med_quarts(
    q1_exp = q1_exp, med_exp = med_exp, q3_exp = q3_exp, n_exp = n_exp,
    q1_nexp = q1_nexp, med_nexp = med_nexp, q3_nexp = q3_nexp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))
  es_med_min_max <- with(x, es_from_med_min_max(
    min_exp = min_exp, med_exp = med_exp, max_exp = max_exp, n_exp = n_exp,
    min_nexp = min_nexp, med_nexp = med_nexp, max_nexp = max_nexp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))
  es_med_min_max_quarts <- with(x, es_from_med_min_max_quarts(
    min_exp = min_exp, q1_exp = q1_exp, med_exp = med_exp, q3_exp = q3_exp, max_exp = max_exp, n_exp = n_exp,
    min_nexp = min_nexp, q1_nexp = q1_nexp, med_nexp = med_nexp, q3_nexp = q3_nexp, max_nexp = max_nexp,
    n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))

  # ANCOVA MEANS ------------------------------------------------------------------
  es_ancova_means_sd <- with(x, es_from_ancova_means_sd(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_exp = ancova_mean_sd_exp, ancova_mean_sd_nexp = ancova_mean_sd_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_se <- with(x, es_from_ancova_means_se(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_se_exp = ancova_mean_se_exp, ancova_mean_se_nexp = ancova_mean_se_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_ci <- with(x, es_from_ancova_means_ci(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_ci_lo_exp = ancova_mean_ci_lo_exp, ancova_mean_ci_up_exp = ancova_mean_ci_up_exp,
    ancova_mean_ci_lo_nexp = ancova_mean_ci_lo_nexp, ancova_mean_ci_up_nexp = ancova_mean_ci_up_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_sd_pooled_adj <- with(x, es_from_ancova_means_sd_pooled_adj(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_pooled = ancova_mean_sd_pooled,
    n_exp = n_exp, n_nexp = n_nexp, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_sd_pooled <- with(x, es_from_ancova_means_sd_pooled_crude(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    mean_sd_pooled = mean_sd_pooled, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  # plot
  es_plot_ancova_means <- with(x, es_from_plot_ancova_means(
    n_exp = n_exp, n_nexp = n_nexp,
    plot_ancova_mean_exp = plot_ancova_mean_exp, plot_ancova_mean_nexp = plot_ancova_mean_nexp,
    plot_ancova_mean_sd_lo_exp = plot_ancova_mean_sd_lo_exp, plot_ancova_mean_sd_lo_nexp = plot_ancova_mean_sd_lo_nexp,
    plot_ancova_mean_sd_up_exp = plot_ancova_mean_sd_up_exp, plot_ancova_mean_sd_up_nexp = plot_ancova_mean_sd_up_nexp,
    plot_ancova_mean_se_lo_exp = plot_ancova_mean_se_lo_exp, plot_ancova_mean_se_lo_nexp = plot_ancova_mean_se_lo_nexp,
    plot_ancova_mean_se_up_exp = plot_ancova_mean_se_up_exp, plot_ancova_mean_se_up_nexp = plot_ancova_mean_se_up_nexp,
    plot_ancova_mean_ci_lo_exp = plot_ancova_mean_ci_lo_exp, plot_ancova_mean_ci_lo_nexp = plot_ancova_mean_ci_lo_nexp,
    plot_ancova_mean_ci_up_exp = plot_ancova_mean_ci_up_exp, plot_ancova_mean_ci_up_nexp = plot_ancova_mean_ci_up_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_plot_ancova_means = reverse_plot_ancova_means
  ))
  # ANCOVA MD ------------------------------------------------------------------
  es_ancova_md_sd <- with(x, es_from_ancova_md_sd(ancova_md = ancova_md, ancova_md_sd = ancova_md_sd,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_var = smd_var, smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_se <- with(x, es_from_ancova_md_se(ancova_md = ancova_md, ancova_md_se = ancova_md_se,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_var = smd_var, smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_ci <- with(x, es_from_ancova_md_ci(ancova_md = ancova_md,
         ancova_md_ci_lo = ancova_md_ci_lo, ancova_md_ci_up = ancova_md_ci_up,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_var = smd_var, smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_pval <- with(x, es_from_ancova_md_pval(
        ancova_md = ancova_md,
        ancova_md_pval = ancova_md_pval,
        cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
        n_exp = n_exp, n_nexp = n_nexp,
        smd_var = smd_var, smd_to_cor = smd_to_cor,
        reverse_ancova_md = reverse_ancova_md))
  # ANCOVA F, T, P ------------------------------------------------------------------
  es_ancova_t <- with(x, es_from_ancova_t(
    ancova_t = ancova_t, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_t = reverse_ancova_t
  ))
  es_ancova_f <- with(x, es_from_ancova_f(
    ancova_f = ancova_f, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_f = reverse_ancova_f
  ))
  es_ancova_t_pval <- with(x, es_from_ancova_t_pval(
    ancova_t_pval = ancova_t_pval, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_t_pval = reverse_ancova_t_pval
  ))
  es_ancova_f_pval <- with(x, es_from_ancova_f_pval(
    ancova_f_pval = ancova_f_pval, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_ancova_f_pval = reverse_ancova_f_pval
  ))

  # CHI-SQ + PHI --------------------------------------------------------
  # yates_chisq column overrides the argument
  yates_chisq_vec <- x$yates_chisq
  yates_chisq_vec[is.na(yates_chisq_vec)] <- as.logical(yates_chisq)[1]

  es_chisq <- with(x, es_from_chisq(
    chisq = chisq, n_sample = n_sample,
    n_cases = n_cases, n_exp = n_exp,
    yates_chisq = yates_chisq_vec,
    reverse_chisq = reverse_chisq
  ))

  es_chisq_pval <- with(x, es_from_chisq_pval(
    chisq_pval = chisq_pval,
    n_cases = n_cases, n_exp = n_exp,
    yates_chisq = yates_chisq_vec,
    n_sample = n_sample, reverse_chisq_pval = reverse_chisq_pval
  ))

  es_phi <- with(x, es_from_phi(
    phi = phi, n_sample = n_sample,
    n_cases = n_cases, n_exp = n_exp,
    reverse_phi = reverse_phi
  ))

  # COR-PB   --------------------------------------------------------
  es_r_point_bis <- with(x, es_from_pt_bis_r(
    pt_bis_r = pt_bis_r, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor,
    reverse_pt_bis_r = reverse_pt_bis_r
  ))

  es_r_point_bis_pval <- with(x, es_from_pt_bis_r_pval(
    pt_bis_r_pval = pt_bis_r_pval,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor,
    reverse_pt_bis_r_pval = reverse_pt_bis_r_pval
  ))

  # 2x2, PROP --------------------------------------------------------
  es_2x2 <- with(x, es_from_2x2(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp,
    # table_2x2_to_cor = table_2x2_to_cor,
    reverse_2x2 = reverse_2x2
  ))
  es_2x2_sum <- with(x, es_from_2x2_sum(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp, n_exp = n_exp, n_nexp = n_nexp,
    # table_2x2_to_cor = table_2x2_to_cor,
    reverse_2x2 = reverse_2x2
  ))
  es_prop <- with(x, es_from_2x2_prop(
    prop_cases_exp = prop_cases_exp, prop_cases_nexp = prop_cases_nexp, n_exp = n_exp, n_nexp = n_nexp,
    # table_2x2_to_cor = table_2x2_to_cor,
    reverse_prop = reverse_prop
  ))

  # RR --------------------------------------------------------
  es_rr_se <- with(x, es_from_rr_se(
    rr = rr, logrr = logrr, logrr_se = logrr_se,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, reverse_rr = reverse_rr
  ))
  es_rr_ci <- with(x, es_from_rr_ci(
    rr = rr, logrr = logrr,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    rr_ci_lo = rr_ci_lo, logrr_ci_lo = logrr_ci_lo, rr_ci_up = rr_ci_up, logrr_ci_up = logrr_ci_up, n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, reverse_rr = reverse_rr,
    max_asymmetry = max_asymmetry
  ))
  es_rr_pval <- with(x, es_from_rr_pval(
    rr = rr, logrr = logrr, rr_pval = rr_pval,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, reverse_rr = reverse_rr_pval
  ))

  # RD (risk difference) --------------------------------------------------------
  es_rd_se <- with(x, es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    reverse_rd = reverse_rd
  ))
  es_rd_ci <- with(x, es_from_rd_ci(
    rd = rd,
    rd_ci_lo = rd_ci_lo, rd_ci_up = rd_ci_up,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    reverse_rd = reverse_rd
  ))
  es_rd_pval <- with(x, es_from_rd_pval(
    rd = rd, rd_pval = rd_pval,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    reverse_rd_pval = reverse_rd_pval
  ))

  # regression
  es_std_beta <- with(x, es_from_beta_std(
    beta_std = beta_std, sd_dv = sd_dv, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_beta_std = reverse_beta_std
  ))
  es_unstd_beta <- with(x, es_from_beta_unstd(
    beta_unstd = beta_unstd, sd_dv = sd_dv, n_exp = n_exp, n_nexp = n_nexp,
    smd_var = smd_var, smd_to_cor = smd_to_cor, reverse_beta_unstd = reverse_beta_unstd
  ))

  # regression t-statistic (partial correlation)
  es_linreg_t <- with(x, es_from_linreg_t(
    linreg_t = linreg_t, n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_t = reverse_linreg_t
  ))

  # regression coefficient (partial correlation)
  es_linreg_b_se <- with(x, es_from_linreg_b_se(
    linreg_b = linreg_b, linreg_b_se = linreg_b_se,
    n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_b = reverse_linreg_b
  ))
  es_linreg_b_ci <- with(x, es_from_linreg_b_ci(
    linreg_b = linreg_b, linreg_b_ci_lo = linreg_b_ci_lo, linreg_b_ci_up = linreg_b_ci_up,
    n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_b = reverse_linreg_b
  ))
  es_linreg_b_pval <- with(x, es_from_linreg_b_pval(
    linreg_b = linreg_b, linreg_b_pval = linreg_b_pval,
    n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_b_pval = reverse_linreg_b_pval
  ))

  # user input
  es_user_crude <- es_from_user_crude(
    user_es_original_measure_crude = x$user_es_original_measure_crude,
    user_es_crude = x$user_es_crude,
    user_se_crude = x$user_se_crude,
    user_ci_lo_crude = x$user_ci_lo_crude,
    user_ci_up_crude = x$user_ci_up_crude,
    user_es_target_measure_crude = measure,
    n_exp = x$n_exp, n_nexp = x$n_nexp, n_sample = x$n_sample,
    n_cases = x$n_cases, n_controls = x$n_controls,
    baseline_risk = x$baseline_risk,
    small_margin_prop = x$small_margin_prop,
    # The merged per-row columns, not the scalar arguments: these five are read from
    # a column of the dataset when one is present (the loop above), and every native
    # route picks the column up through with(x, ...). These two calls sit outside
    # with(), so passing the scalars gave every user-input row the dataset default
    # with no message -- silently substituting, for instance, the biserial correlation
    # for the point-biserial one the row asked for. .dispatch_user_conversion() takes
    # a vector or a scalar and subsets it by the row indices each branch keeps
    # (.user_method_by_row()), so the column arrives aligned.
    or_to_rr = x$or_to_rr, or_to_cor = x$or_to_cor,
    smd_to_cor = x$smd_to_cor, cor_to_smd = x$cor_to_smd,
    rr_to_or = x$rr_to_or,
    alpha_to_es = alpha_to_es, icc_to_es = icc_to_es, prop_to_es = prop_to_es,
    omega_to_es = omega_to_es
  )
  es_user_adj <- es_from_user_adj(
    user_es_original_measure_adj = x$user_es_original_measure_adj,
    user_es_adj = x$user_es_adj,
    user_se_adj = x$user_se_adj,
    user_ci_lo_adj = x$user_ci_lo_adj,
    user_ci_up_adj = x$user_ci_up_adj,
    user_es_target_measure_adj = measure,
    n_exp = x$n_exp, n_nexp = x$n_nexp, n_sample = x$n_sample,
    n_cases = x$n_cases, n_controls = x$n_controls,
    baseline_risk = x$baseline_risk,
    small_margin_prop = x$small_margin_prop,
    # The merged per-row columns, as in es_from_user_crude() above.
    or_to_rr = x$or_to_rr, or_to_cor = x$or_to_cor,
    smd_to_cor = x$smd_to_cor, cor_to_smd = x$cor_to_smd,
    rr_to_or = x$rr_to_or,
    alpha_to_es = alpha_to_es, icc_to_es = icc_to_es, prop_to_es = prop_to_es,
    omega_to_es = omega_to_es
  )
  # survival
  es_cases_time <- with(x, es_from_cases_time(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    time_exp = time_exp, time_nexp = time_nexp,
    baseline_rate = baseline_rate,
    reverse_irr = reverse_irr
  ))
  # variability
  var_means_sd <- with(x, es_variab_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = mean_sd_exp, mean_sd_nexp = mean_sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  var_means_se <- with(x, es_variab_from_means_se(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_se_exp = mean_se_exp, mean_se_nexp = mean_se_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  var_means_ci <- with(x, es_variab_from_means_ci(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
    mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  # ========================================= STEP 2. SYNTHESIS OF CALCULATIONS =========================================== #

  smd_list_L1 = list(es_cohen_d = es_cohen_d,
                es_hedges_g = es_hedges_g)
  # Reported precision (a logOR SE, a reported CI, or a reported p-value) outranks
  # es_odds_ratio, whose SE is *simulated* from the cell marginals (es_from_or) and runs
  # ~1.5x wide. A p-value recovers the SE exactly -- |logOR / qnorm(p/2)| -- so it belongs
  # with the SE and the CI, not below the reconstruction: on the (40,160,20,180) table the
  # reconstruction gives 0.4390 against the exact 0.2946, i.e. 2.22x too little weight, and
  # the inflation grows with study size. This is what the @note on es_from_or() already
  # tells the user to do. It also aligns the OR order with the RR list
  # (es_rr_se, es_rr_ci, es_rr_pval first).
  or_list_L2 = list(es_odds_ratio_se = es_odds_ratio_se,
                 es_odds_ratio_ci = es_odds_ratio_ci,
                 es_odds_ratio_pval = es_odds_ratio_pval,
                 es_odds_ratio = es_odds_ratio,
                 es_logreg_t = es_logreg_t)
  rr_list_L3 = list(es_rr_se = es_rr_se,
                     es_rr_ci = es_rr_ci,
                     es_rr_pval = es_rr_pval)
  rd_list_L30 = list(es_rd_se = es_rd_se,
                     es_rd_ci = es_rd_ci,
                     es_rd_pval = es_rd_pval)
  cor_list_L4 = list(es_pearson_r = es_pearson_r,
                  es_fisher_z = es_fisher_z,
                  es_spearman_r = es_spearman_r)
  irr_list_L5 = list(es_cases_time = es_cases_time)

  var_list_L6 = list(var_means_sd = var_means_sd,
                     var_means_se = var_means_se,
                     var_means_ci = var_means_ci)

  contingency_list_L7 = list(es_2x2 = es_2x2,
                         es_2x2_sum = es_2x2_sum,
                         es_prop = es_prop)

  phi_chisq_list_L8 = list(es_chisq = es_chisq,
                           es_phi = es_phi,
                           es_chisq_pval = es_chisq_pval)

  means_post_list_L9 = list(es_means_sd_raw = es_means_sd_raw,
                         es_means_sd_pooled = es_means_sd_pooled,
                         es_means_se_raw = es_means_se_raw,
                         es_means_ci_raw = es_means_ci_raw)
  md_post_list_L10 = list(es_md_sd = es_md_sd,
                      es_md_se = es_md_se,
                      es_md_ci = es_md_ci,
                      es_md_pval = es_md_pval)
  anova_list_L11 = list(es_t_student = es_t_student,
                    es_anova_f = es_anova_f,
                    es_r_point_bis = es_r_point_bis,
                    es_t_student_pval = es_t_student_pval,
                    es_anova_f_pval = es_anova_f_pval,
                    es_r_point_bis_pval = es_r_point_bis_pval,
                    es_etasq = es_etasq)
  medians_list_L12 = list(es_med_quarts = es_med_quarts,
                      es_med_min_max = es_med_min_max,
                      es_med_min_max_quarts = es_med_min_max_quarts)
  regression_list_L13 = list(es_std_beta = es_std_beta,
                           es_unstd_beta = es_unstd_beta)
  partial_cor_list_L27 = list(es_linreg_t = es_linreg_t,
                               es_linreg_b_se = es_linreg_b_se,
                               es_linreg_b_ci = es_linreg_b_ci,
                               es_linreg_b_pval = es_linreg_b_pval)

  md_paired_list_L14 = list(es_means_change_sd = es_means_change_sd,
                   es_mean_change_se = es_mean_change_se,
                   es_mean_change_ci = es_mean_change_ci,
                   es_mean_change_pval = es_mean_change_pval)
  means_paired_list_L15 = list(es_means_sd_pre_post = es_means_sd_pre_post,
                           es_means_se_pre_post = es_means_se_pre_post,
                           es_means_ci_pre_post = es_means_ci_pre_post)
  anova_paired_list_L16 = list(es_paired_t = es_paired_t,
                           es_paired_t_pval = es_paired_t_pval,
                           es_paired_f = es_paired_f,
                           es_paired_f_pval = es_paired_f_pval)

  es_cohen_d_adj_L17 = list(es_cohen_d_adj = es_cohen_d_adj)
  ancova_adjusted_list_L18 = list(es_ancova_t = es_ancova_t,
                                  es_ancova_f = es_ancova_f,
                                  es_ancova_t_pval = es_ancova_t_pval,
                                  es_ancova_f_pval = es_ancova_f_pval,
                                  es_etasq_adj = es_etasq_adj)
  # Split by attenuation tier, as md_adjusted_list_L20a/b is, and for the same reason.
  # L19a's three routes are handed a standardizer directly -- per-arm residual SDs, a
  # marginal pooled SD, or an adjusted pooled SD -- so their point estimate is unaffected
  # by covariate imbalance. L19b's two must instead recover the standardizer from a
  # reported SE or CI, which sets the Lai & Kelley leverage term D to zero and attenuates
  # d and its SE together (measured at delta_x = 1, n = 60/60, rho = 0.5: d = 1.0960 for
  # ancova_means_se against 1.1635 for ancova_md_sd, the closed form
  # 1/sqrt(1 + n_j*D/4) = 0.941923).
  #
  # The split has to be carried into the assembly below, not merely into this list: with
  # the whole of L19 ranked ahead of md_adjusted_list_L20a, a row reporting adjusted
  # means with their SEs alongside an adjusted MD with the residual SD had the attenuated
  # estimate selected while the unattenuated one sat in the same row. L19b carries only
  # half the D contribution, so it belongs between md_sd and the full-attenuation
  # statistic routes of L18.
  means_adjusted_list_L19a = list(es_ancova_means_sd = es_ancova_means_sd,
                             es_ancova_means_sd_pooled = es_ancova_means_sd_pooled,
                             es_ancova_means_sd_pooled_adj = es_ancova_means_sd_pooled_adj)
  means_adjusted_list_L19b = list(es_ancova_means_se = es_ancova_means_se,
                             es_ancova_means_ci = es_ancova_means_ci)

  # The adjusted MD routes split into two attenuation tiers under covariate imbalance
  # (Lai & Kelley 2012; see the @note blocks of the es_from_ancova_* functions).
  # es_ancova_md_sd receives the adjusted MD together with the residual SD, so its
  # point estimate is unaffected by imbalance. The other three must invert a reported
  # SE / CI / p-value, which was built with the leverage term D = (xbar1-xbar2)^2/SS_x
  # that the wide format cannot carry, so they attenuate d (and its SE) by
  # 1/sqrt(1 + D/(1/n_exp + 1/n_nexp)), the same tier as the ANCOVA test statistics.
  # They are therefore ranked separately: L20a above ancova_adjusted_list_L18, L20b below.
  md_adjusted_list_L20a = list(es_ancova_md_sd = es_ancova_md_sd)
  md_adjusted_list_L20b = list(es_ancova_md_se = es_ancova_md_se,
                               es_ancova_md_ci = es_ancova_md_ci,
                               es_ancova_md_pval = es_ancova_md_pval)


  plot_raw_L21 = list(es_plot_means_raw)
  plot_adjusted_L22 = list(es_plot_ancova_means)
  user_raw_L23 = list(es_user_crude)
  es_user_adj_L24 = list(es_user_adj)

  within_group_list_L25 = list(es_means_pp_sg = es_means_pp_sg,
                               es_means_se_pp_sg = es_means_se_pp_sg,
                               es_means_ci_pp_sg = es_means_ci_pp_sg,
                               es_mean_change_sg = es_mean_change_sg,
                               es_mean_change_se_sg = es_mean_change_se_sg,
                               es_mean_change_ci_sg = es_mean_change_ci_sg,
                               es_mean_change_pval_sg = es_mean_change_pval_sg,
                               es_paired_t_sg = es_paired_t_sg)

  prop_list_L26 = list(es_prop_sg = es_prop_sg,
                       es_prop_counts_sg = es_prop_counts_sg)

  alpha_list_L28 = list(es_alpha_sg = es_alpha_sg)
  omega_list_L31 = list(es_omega_sg = es_omega_sg)

  icc_list_L29 = list(es_icc_sg = es_icc_sg)

  USER_crude = user_raw_L23
  USER_adjusted = es_user_adj_L24

  SMD_post = c(smd_list_L1, means_post_list_L9,
               anova_list_L11, regression_list_L13, md_post_list_L10,
               medians_list_L12, plot_raw_L21)

  SMD_paired = c(means_paired_list_L15, anova_paired_list_L16, md_paired_list_L14)

  SMD_adjusted = c(es_cohen_d_adj_L17, means_adjusted_list_L19a, md_adjusted_list_L20a,
                   means_adjusted_list_L19b, ancova_adjusted_list_L18,
                   md_adjusted_list_L20b, plot_adjusted_L22)

  OR = c(or_list_L2)

  CONT = contingency_list_L7

  RR = rr_list_L3

  RD_stand = rd_list_L30

  PHI = phi_chisq_list_L8

  COR = cor_list_L4

  IRR = irr_list_L5

  VAR = var_list_L6


  if (es_selected == "auto") {
    if (measure %in% c("d", "g", "md")) {
      # OR and the raw 2x2 table (CONT) are the sole gateway into the SMD
      # family via the Cox transform d = log(OR) * sqrt(3)/pi. RD_stand is
      # not used here: RD -> SMD is not identified without an assumed
      # baseline_risk, and its SE is anti-conservative because baseline_risk
      # is treated as a fixed constant. RD keeps its legitimate ratio
      # conversions below.

      if (selection_auto == "crude") {
        res = c(USER_crude, SMD_post, SMD_paired, OR,
                CONT, COR, PHI, SMD_adjusted,
                USER_adjusted,
                #not used
                RR, RD_stand, IRR, VAR)

      } else if (selection_auto == "adjusted") {
        res = c(USER_adjusted, SMD_adjusted,
                USER_crude, SMD_post,
                SMD_paired, OR, CONT, COR,
                PHI,
                #not used
                RR, RD_stand, IRR, VAR)
      } else if (selection_auto == "paired") {
        res = c(USER_crude, SMD_paired, SMD_post, SMD_adjusted,
                USER_adjusted, OR, CONT, COR,
                PHI,
                #not used
                RR, RD_stand, IRR, VAR)
      } else {
        stop("selection_auto should be one of 'crude', 'adjusted' or 'paired'")
      }

    } else if (measure %in% c("logor", "or")) {
      res = c(USER_crude, OR, CONT, RR, RD_stand, PHI, COR, SMD_post, USER_adjusted,
              #not used
              SMD_paired, SMD_adjusted, IRR, VAR)

    } else if (measure %in% c("logrr", "rr")) {
      res = c(USER_crude, RR, CONT, OR, RD_stand, PHI, USER_adjusted,
              #not used
              SMD_post, SMD_paired,
              COR, SMD_adjusted,
              IRR, VAR)

    } else if (measure %in% c("logirr", "irr")) {
      res = c(USER_crude, IRR, USER_adjusted,
              #not used
              SMD_post, SMD_paired, OR,
              CONT, COR, PHI, SMD_adjusted,
              RR, RD_stand, VAR)

    } else if (measure %in% c("loghr", "hr")) {
      # HR: user input only
      res = c(USER_crude, USER_adjusted)

    } else if (measure %in% c("nnt")) {
      res = c(USER_crude, RD_stand, CONT, OR, RR, IRR, PHI, USER_adjusted,
              #not used
              SMD_post, SMD_paired, COR, SMD_adjusted,
              VAR)

    } else if (measure %in% c("rd")) {
      res = c(USER_crude, RD_stand, CONT, OR, RR, IRR, PHI, USER_adjusted,
              #not used
              SMD_post, SMD_paired, COR, SMD_adjusted,
              VAR)

    } else if (measure %in% c("r", "z")) {
      # As for the SMD family, OR and the raw 2x2 (CONT) are the sole gateway
      # into the correlation family (r/z are derived downstream of d_se).
      # RD_stand stays #not used for the same identification/SE reason.
      res = c(USER_crude, COR, CONT, OR, PHI, SMD_post,
              SMD_paired, USER_adjusted,
              #not used
              SMD_adjusted,
              RR, RD_stand, IRR, VAR)

    } else if (measure %in% c("logvr", "logcvr")) {
      res = c(USER_crude, VAR, SMD_post, SMD_paired, USER_adjusted,
              #not used
              OR, CONT, RD_stand, COR, PHI, SMD_adjusted,
              RR, IRR)

    } else if (measure %in% c("dw", "gw", "mdw")) {
      res = c(USER_crude, within_group_list_L25)

    } else if (measure == "prop") {
      res = c(USER_crude, prop_list_L26)

    } else if (measure == "alpha") {
      res = c(USER_crude, alpha_list_L28)
    } else if (measure == "omega") {
      res = c(USER_crude, omega_list_L31)

    } else if (measure == "icc") {
      res = c(USER_crude, icc_list_L29)

    } else if (measure %in% c("rp", "zp")) {
      res = c(USER_crude, partial_cor_list_L27)

    }
  } else {
    if (measure %in% c("dw", "gw", "mdw")) {
      res = c(USER_crude, within_group_list_L25)
    } else if (measure == "prop") {
      res = c(USER_crude, prop_list_L26)
    } else if (measure == "alpha") {
      res = c(USER_crude, alpha_list_L28)
    } else if (measure == "omega") {
      res = c(USER_crude, omega_list_L31)
    } else if (measure == "icc") {
      res = c(USER_crude, icc_list_L29)
    } else if (measure %in% c("rp", "zp")) {
      res = c(USER_crude, partial_cor_list_L27)
    } else if (measure %in% c("loghr", "hr")) {
      res = c(USER_crude, USER_adjusted)
    } else if (measure %in% c("logor", "or")) {
      # Non-auto modes mirror the auto-mode ordering for every measure, so the
      # fallback selection (e.g. hierarchy mode with a default hierarchy) agrees
      # with es_selected = "auto": here, preferring the directly reported OR
      # over a reconstructed 2x2.
      res = c(USER_crude, OR, CONT, RR, RD_stand, PHI, COR, SMD_post, USER_adjusted,
              #not used
              SMD_paired, SMD_adjusted, IRR, VAR)
    } else if (measure %in% c("logrr", "rr")) {
      # Prefer the directly reported RR over a reconstructed 2x2/OR (mirrors auto).
      res = c(USER_crude, RR, CONT, OR, RD_stand, PHI, USER_adjusted,
              #not used
              SMD_post, SMD_paired,
              COR, SMD_adjusted,
              IRR, VAR)
    } else if (measure %in% c("logirr", "irr")) {
      res = c(USER_crude, IRR, USER_adjusted,
              #not used
              SMD_post, SMD_paired, OR,
              CONT, COR, PHI, SMD_adjusted,
              RR, RD_stand, VAR)
    } else if (measure %in% c("nnt", "rd")) {
      # Ratio-family targets: RD_stand is a legitimate source (RD -> OR/RR/NNT).
      res = c(USER_crude, RD_stand, CONT, OR, RR, IRR, PHI, USER_adjusted,
              #not used
              SMD_post, SMD_paired, COR, SMD_adjusted,
              VAR)
    } else if (measure %in% c("r", "z")) {
      # OR and the raw 2x2 are the sole gateway into the correlation family;
      # RD_stand stays #not used (RD -> r/z not identified without an assumed
      # baseline_risk; anti-conservative SE). Mirrors auto.
      res = c(USER_crude, COR, CONT, OR, PHI, SMD_post,
              SMD_paired, USER_adjusted,
              #not used
              SMD_adjusted,
              RR, RD_stand, IRR, VAR)
    } else if (measure %in% c("logvr", "logcvr")) {
      res = c(USER_crude, VAR, SMD_post, SMD_paired, USER_adjusted,
              #not used
              OR, CONT, RD_stand, COR, PHI, SMD_adjusted,
              RR, IRR)
    } else {
      # d/g/md targets. selection_auto (crude/adjusted/paired) does not apply
      # outside es_selected = "auto", so the crude ordering is used. OR and the
      # raw 2x2 are the sole gateway into the SMD family; RD_stand is demoted
      # (RD -> SMD is not identified without an assumed baseline_risk and yields
      # anti-conservative SEs).
      res = c(USER_crude, SMD_post, SMD_paired, OR,
              CONT, COR, PHI, SMD_adjusted,
              USER_adjusted,
              #not used
              RR, RD_stand, IRR, VAR)
    }
  }

  df_sq_list = list(
    es_odds_ratio_pval,  es_rr_pval, es_chisq, es_chisq_pval,
    es_md_pval, es_anova_f,  es_anova_f_pval, es_r_point_bis_pval,
    es_etasq, es_mean_change_pval, es_paired_t_pval, es_paired_f,
    es_paired_f_pval, es_ancova_f, es_ancova_f_pval, es_etasq_adj, es_ancova_md_pval
  )
  warning_reverse = FALSE
  for (i in seq_along(df_sq_list)) {
    df <- df_sq_list[[i]]
    ci_lo_cols <- grep("ci_lo", names(df), value = TRUE)

    if (length(ci_lo_cols) > 0) {
      for (col in ci_lo_cols) {
        if (any(!is.na(df[[col]]))) {
          warning_reverse = TRUE
          break
        }
      }
    }
  }

  if (warning_reverse & verbose) {
    warning("When you enter input data that cannot be negative (F-test, eta-squared, p-value, or chi-square values), do not forget to properly set up the direction of the generated effect size using corresponding reverse_* argument!")
  }
  expected_length <- if (measure %in% c("dw", "gw", "mdw")) {
    9
  } else if (measure == "prop") {
    3
  } else if (measure %in% c("alpha", "omega", "icc", "loghr", "hr")) {
    2
  } else if (measure %in% c("rp", "zp")) {
    5
  } else {
    76
  }
  if (length(res) != expected_length) {
    stop ("failure to estimate all effect sizes (not expected, you trigerred a bug, please contact cgosling@parisnanterre.fr). Expected ", expected_length, ", got ", length(res))
  }
  class(res) <- "metaConvert"
  attr(res, "raw_data") <- x
  attr(res, "exp") <- exp
  attr(res, "measure") <- measure
  attr(res, "split_adjusted") <- split_adjusted
  attr(res, "es_selected") <- es_selected
  attr(res, "format_adjusted") <- format_adjusted
  attr(res, "hierarchy") <- hierarchy
  attr(res, "main_es") <- main_es
  flag_opts <- .default_flag_options()
  flag_opts[names(flag_options)] <- flag_options
  attr(res, "flag_options") <- flag_opts
  attr(res, "input_validation") <- validation$issues
  # rows where r_pre_post was defaulted (used by the summary flags)
  attr(res, "r_defaulted") <- .r_defaulted
  attr(res, "alpha_to_es") <- alpha_to_es
  attr(res, "omega_to_es") <- omega_to_es
  attr(res, "icc_to_es") <- icc_to_es
  attr(res, "icc_agreement_se") <- icc_agreement_se
  # Where a reliability SE is allowed to come from. Stored beside icc_agreement_se so
  # summary() can see it: which columns a row still needs (the `guidance` output)
  # depends on it -- under omega_se_source = "closed_form" an SE-less omega is missing
  # n_sample, under the default "reported" it is missing omega_se, and naming the wrong
  # one sends the user after a statistic that route cannot use.
  attr(res, "alpha_se_source") <- alpha_se_source
  attr(res, "omega_se_source") <- omega_se_source
  attr(res, "prop_to_es") <- prop_to_es
  # standardizer used by the pre/post routes (drives the E6/E7 mixing flags)
  attr(res, "pre_post_to_smd") <- pre_post_to_smd
  # Which transform each route put in the `z` column: Fisher's atanh(r) for the
  # correlation and binary families, a variance-stabilising transform for the SMD
  # family. Derived from the frames themselves rather than from a list of route
  # names, and consumed by the E8 mixing flag. See .z_transform_by_route().
  attr(res, "z_transform") <- .z_transform_by_route(res)
  attr(res, "pool_sd") <- pool_sd
  # per-row endpoint-SMD standardizer ("pooled"/"control"/"control_robust"),
  # used by the summary A6 CI-width check (Glass CIs are on qt(.975, n_nexp-1))
  attr(res, "smd_denom_used") <- .normalize_smd_denom(x[, "smd_denom"])
  # Conversion-formula choices, recorded so es_formulas() can re-run the
  # conversion under each alternative formula without the user restating the call.
  attr(res, "conversion_args") <- list(
    table_2x2_to_cor = table_2x2_to_cor, rr_to_or = rr_to_or, or_to_rr = or_to_rr,
    or_to_cor = or_to_cor, smd_to_cor = smd_to_cor, smd_var = smd_var,
    smd_denom = smd_denom, pre_post_to_smd = pre_post_to_smd,
    r_pre_post = r_pre_post[1], cor_to_smd = cor_to_smd, unit_type = unit_type,
    yates_chisq = yates_chisq, pool_sd = pool_sd, prop_to_es = prop_to_es,
    alpha_to_es = alpha_to_es, omega_to_es = omega_to_es, icc_to_es = icc_to_es,
    # Where a reliability SE comes from is part of the analysis es_formulas() has to
    # reproduce: without these three, a re-run silently reverts to the defaults and
    # prints standard errors the user deliberately refused (alpha_se_source =
    # "reported", omega_se_source = "reported", icc_agreement_se = "drop" all exist to
    # leave a row visible but unpooled).
    alpha_se_source = alpha_se_source, omega_se_source = omega_se_source,
    icc_agreement_se = icc_agreement_se,
    max_asymmetry = max_asymmetry, correct_inputs = correct_inputs,
    selection_auto = selection_auto
  )
  return(res)
}
# x_save2 = x; list_df = df_es; ordering = ordering_crude; digits = digits;
# suffix = "_crude"; measure = measure
#
# x = dat; exp = TRUE; es_selected = "hierarchy"
# split_adjusted = TRUE; measure = "d_to_or, smd_to
# hierarchy = "means_sd"#hierarchy#"user_input_crude"
# digits = 3
# x$user_es_measure_crude = "ROR"
# x$user_es_measure_adj = "MOR"
# x$split_adjustedusted = x$se_adjusted = 2
# load_all(); View(convert_df(df.haza));


# x = dat
# measure = "d"
# main_es = TRUE
# split_adjusted = TRUE
# format_adjusted = "wide"
# verbose = TRUE
# es_selected = "hierarchy"
# phi_to_cor = "tetrachoric";
# chisq_to_cor = "tetrachoric";
# rr_to_or = "metaumbrella"
# or_to_rr = "metaumbrella_cases"
# or_to_cor = "bonett"
# cor_to_smd = "cooper"
# table_2x2_to_cor = "lipsey"
# yates_chisq = FALSE
# smd_to_cor = "viechtbauer"
# pre_post_to_smd = "morris"
# unit_type = "raw_scale"
# hierarchy = "rr_se"
# start_time <- Sys.time()
# end_time <- Sys.time()
# end_time - start_time
