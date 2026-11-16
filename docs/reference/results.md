# Results reference

This page describes what each command stores, how unavailable effects and uncertainty are reported, and how results are saved. The help file of each command lists its stored results in Stata's format. [Working with results](../user-results.md) shows how results move between commands.

## Results of cdid with gvar()

`cdid` posts its results in `e()`. The display reports the available cells; every requested cell, including unavailable ones, is kept with the results and written by `cdid_export, type(cells)`.

| Scalars | Content |
| --- | --- |
| `e(N)` | number of units in the study population |
| `e(N_full)` | number of units in the study population, the scale of `e(V_asymptotic)` |
| `e(N_cells)` | number of group-time cells of the fit |
| `e(N_targets)` | number of rows of the results table, including unavailable and reference cells |
| `e(N_posted)` | number of estimates in `e(b)` |
| `e(clustered)` | 1 if cluster() was specified |
| `e(N_clust)` | number of clusters, with cluster() |
| `e(point_only)` | 1 if only point estimates were computed (noinffunc) |
| `e(level)` | confidence level |
| `e(alpha)` | 1 - `e(level)`/100 |
| `e(B_requested)` | bootstrap draws requested |
| `e(bstrap_requested)` | value of bstrap() |
| `e(cband_requested)` | value of cband() |
| `e(critical_value)` | critical value of the reported intervals or bands |
| `e(normal_critical_value)` | normal critical value at `e(level)` |
| `e(wald_statistic)` | pretrend Wald statistic |
| `e(wald_df)` | its degrees of freedom |
| `e(wald_p)` | its p-value |
| `e(wald_N_clust)` | clusters used by the Wald test, with cluster() |
| `e(anticipation)` | periods of anticipation |

| Macros | Content |
| --- | --- |
| `e(cmd)` | cdid |
| `e(cmdline)` | command as typed |
| `e(title)` | title of the output |
| `e(mode)` | chained |
| `e(depvar)` | name of the outcome |
| `e(ivar)` | unit identifier |
| `e(time)` | date variable |
| `e(gvar)` | cohort variable |
| `e(covariates)` | covariates, as typed |
| `e(clustvar)` | cluster variable, with cluster() |
| `e(method)` | treatment adjustment |
| `e(control_group)` | comparison group |
| `e(base_period)` | base period |
| `e(pt_assumption)` | parallel trends assumption |
| `e(gmm_weighting)` | GMM weighting |
| `e(gmm_covariance)` | covariance regularization |
| `e(outcome_pooling)` | outcome-regression pooling |
| `e(imputation_pooling)` | imputation pooling, when active |
| `e(selection)` | assumption about missing outcomes, as specified |
| `e(selection_canonical)` | the same, in canonical form |
| `e(selection_est_method)` | correction for missing outcomes |
| `e(selection_model)` | observation model |
| `e(selection_pooling)` | observation-model pooling |
| `e(selection_universe)` | units that fit the observation model |
| `e(selection_justification)` | text of seljustification(), when specified |
| `e(cross_pair)` | crosspair, when specified |
| `e(preset)` | cohort_time_fe, when specified |
| `e(missing_reference)` | exclude or error (not with the preset) |
| `e(reference_status)` | restricted, auxiliary_unavailable, or unrestricted (not with the preset) |
| `e(reference_restricted)` | cohorts whose comparisons were restricted, when any |
| `e(reference_variables)` | covariates whose reference values were checked, when comparisons were restricted |
| `e(inactive_options)` | options that have no effect under this specification, specified or not |
| `e(callback)` | callback(), with method(custom) |
| `e(payload)` | payload(), with method(custom) |
| `e(wtype)` | weight type, with weights |
| `e(wexp)` | weight expression, with weights |
| `e(nuisance)` | availability of the reports kept, with keepnuisance |
| `e(vce)` | analytic, bootstrap, or none |
| `e(vcetype)` | Bootstrap, with the multiplier bootstrap |
| `e(band)` | pointwise, simultaneous, or none |
| `e(inference_rule)` | how bstrap() and cband() were applied, in words |
| `e(rng)` | source of the bootstrap multipliers, in words |
| `e(result_kind)` | cells |
| `e(result_status)` | status of the results |
| `e(inference_status)` | status of the inference |
| `e(wald_status)` | status of the pretrend Wald test |
| `e(wald_reason)` | why the Wald test was not computed, when it was not |
| `e(warning)` | warnings, when any |
| `e(posting_mode)` | matrix or state |
| `e(posting_reason)` | description of what `e(b)` holds |
| `e(coefficient_map)` | description of `e(coef_map)` |
| `e(covariance_scale)` | description of `e(V)` and `e(V_asymptotic)` |
| `e(se_scale)` | description of `e(se)` |
| `e(sample_scope)` | description of `e(sample)` |
| `e(state_token)` | reference to the stored results used by estat and the postestimation commands |
| `e(state_file)` | file written by saving(), when specified |
| `e(estat_cmd)` | cdid_estat |
| `e(properties)` | b V, or b with noinffunc |

| Matrices | Content |
| --- | --- |
| `e(b)` | available estimates |
| `e(V)` | analytical covariance of `e(b)` |
| `e(V_asymptotic)` | `e(N_full)` times `e(V)` |
| `e(se)` | reported standard errors |
| `e(analytic_se)` | analytical standard errors |
| `e(ci)` | reported intervals or bands (2 x k) |
| `e(ci_pointwise)` | pointwise intervals (2 x k) |
| `e(coef_map)` | cell of each estimate in `e(b)` |
| `e(reference_eligibility)` | eligible and excluded units by target cohort, when comparisons were restricted |
| `e(reference_contexts)` | excluded units by reference date, when comparisons were restricted |
| `e(reference_aux_unavailable)` | unavailable pretrend comparisons, when any |

| Functions | Content |
| --- | --- |
| `e(sample)` | marks the input rows of the units in the study population, including rows with a missing outcome; units excluded from the comparisons of a cohort, for example under missingreference(exclude), remain marked |

With noinffunc, `e(V)` and `e(V_asymptotic)` are not posted and `e(point_only)` is 1. With posting(matrix) and at least one available cell, `e(se)`, `e(analytic_se)`, `e(ci)`, and `e(ci_pointwise)` have missing entries and `e(properties)` is b; with posting(state), no matrix is posted and `e(properties)` is empty.

`e(sample)` marks the observations used: in chained mode it includes rows whose outcome is missing, because they determine when each unit is observed. `cdid_load` does not set it.

### Ordinary mode

In ordinary mode, `cdid` stores the following; the chained-mode settings macros are not posted.

| Scalars | Content |
| --- | --- |
| `e(N)` | number of units used |
| `e(N_full)` | the same, the scale of `e(V_asymptotic)` |
| `e(N_cells)` | number of group-time cells |
| `e(N_targets)` | number of rows of the results table |
| `e(N_posted)` | number of estimates in `e(b)` |
| `e(clustered)` | 1 when the reported inference uses clusters |
| `e(N_clust)` | number of clusters, when clustered |
| `e(point_only)` | 1 if only point estimates were computed |
| `e(level)` | confidence level |
| `e(alpha)` | 1 - `e(level)`/100 |
| `e(B_requested)` | bootstrap draws requested |
| `e(bstrap)` | value of bstrap() |
| `e(cband)` | value of cband() |
| `e(critical_value)` | critical value of the reported intervals or bands |
| `e(normal_critical_value)` | normal critical value at `e(level)` |
| `e(wald_statistic)` | pretrend Wald statistic |
| `e(wald_df)` | its degrees of freedom |
| `e(wald_p)` | its p-value |
| `e(wald_p_unrounded)` | its p-value before rounding |

| Macros | Content |
| --- | --- |
| `e(cmd)` | cdid |
| `e(cmdline)` | command as typed |
| `e(title)` | title of the output |
| `e(mode)` | ordinary |
| `e(depvar)` | name of the outcome |
| `e(ivar)` | unit identifier |
| `e(time)` | date variable |
| `e(gvar)` | cohort variable |
| `e(covariates)` | covariates, as typed |
| `e(clustvar)` | cluster variable used for inference, when clustered |
| `e(method)` | treatment adjustment |
| `e(inactive_options)` | accepted options that have no effect in ordinary mode, specified or not |
| `e(callback)` | callback(), with method(custom) |
| `e(payload)` | payload(), with method(custom) |
| `e(vce)` | analytic, bootstrap, or none |
| `e(vcetype)` | Bootstrap, with the multiplier bootstrap |
| `e(band)` | pointwise, simultaneous, or none |
| `e(inference_rule)` | how bstrap() and cband() were applied, in words |
| `e(result_kind)` | ordinary_cells (ordinary_simple, ... after estat) |
| `e(wald_status)` | status of the pretrend Wald test |
| `e(wald_reason)` | cells used by the Wald test, or why it was not computed |
| `e(conditional_test_status)` | status of the conditional pretest in this result |
| `e(warning)` | warnings, when any |
| `e(posting_mode)` | matrix or state |
| `e(covariance_scale)` | description of `e(V)`, `e(V_asymptotic)`, and `e(V_did)` |
| `e(se_scale)` | description of `e(se)` |
| `e(sample_scope)` | description of `e(sample)` |
| `e(state_token)` | reference to the stored results used by estat |
| `e(estat_cmd)` | cdid_estat |
| `e(properties)` | b V |

| Matrices | Content |
| --- | --- |
| `e(b)` | available estimates |
| `e(V)` | covariance of `e(b)` |
| `e(V_asymptotic)` | `e(N_full)` times `e(V)` |
| `e(V_did)` | covariance R's att_gt() reports as V_analytical, divided by `e(N_full)` |
| `e(V_did_asymptotic)` | covariance R's att_gt() reports as V_analytical |
| `e(se)` | reported standard errors |
| `e(analytic_se)` | analytical standard errors |
| `e(ci)` | reported intervals or bands (2 x k) |
| `e(ci_pointwise)` | pointwise intervals (2 x k) |
| `e(coef_map)` | cell of each estimate in `e(b)` |

| Functions | Content |
| --- | --- |
| `e(sample)` | marks the observations used |

## Unavailable effects

An unavailable cell has missing estimates and inference and a reason, written by `cdid_export, type(cells)` in the variable `reason` and listed by `cdid_diagnostics, detail`. With `posting(matrix)`, only available estimates enter `e(b)`; `e(coef_map)` gives the cell of each one. An aggregate that needs an unavailable cell is unavailable unless `estat` is given `narm`, in which case it describes the restricted target and records the excluded cells.

### Covariates missing at a reference date

Under `missingreference(exclude)`, chained cohort fits post:

| Result | Content |
| --- | --- |
| `e(missing_reference)` | `exclude` or `error`; not posted with the preset. |
| `e(reference_status)` | `restricted` (some cohort comparisons were restricted), `auxiliary_unavailable` (only pretrend comparisons became unavailable) or `unrestricted`. |
| `e(reference_eligibility)` | Restricted fits: one row per target cohort (rows `g3`, `g4`, ...) with columns `cohort`, `reference_date`, `n_treated`, `n_treated_eligible`, `n_treated_excluded`, `n_never_treated_excluded`, `n_other_cohort_excluded` and `restricted`. The excluded counts are taken across all supplied units at the cohort's reference date, including units outside the cohort's comparison groups. |
| `e(reference_contexts)` | Restricted fits: one row per restricted reference date, columns `reference_date` and `n_excluded`. |
| `e(reference_aux_unavailable)` | One row per unavailable pretrend comparison, columns `cohort`, `s` and `t`. |
| `e(reference_variables)` | The covariates whose reference values were checked. |
| `e(reference_restricted)` | The restricted cohorts. |

The record is saved by `saving()` and restored by `cdid_load`, and `estat` aggregations use the restricted cohort masses.

## Aggregated results

`estat simple|group|dynamic|calendar` replaces `e()` with the aggregation. The four types average the post-treatment cells overall, by cohort, by time since adoption and by calendar date; `estat dynamic` also reports the pretreatment event times, and each type reports an overall estimate. Cohort shares are computed over the study population.

| Scalars | Content |
| --- | --- |
| `e(N_cells)` | number of group-time cells of the fit |
| `e(N_targets)` | number of aggregation bins |
| `e(N_posted)` | number of estimates in `e(b)` |
| `e(bstrap)` | bootstrap flag used |
| `e(cband)` | band flag used |
| `e(overall_att)` | overall estimate |
| `e(overall_se)` | its reported standard error |
| `e(overall_analytic_se)` | its analytical standard error |
| `e(overall_ci_low)` | lower limit of its pointwise interval |
| `e(overall_ci_high)` | upper limit of its pointwise interval |
| `e(overall_n_draws)` | bootstrap draws for the overall estimate |
| `e(overall_draws_used)` | 1 if bootstrap draws were used for it |
| `e(n_draws)` | bootstrap draws for the bins |
| `e(draws_used)` | 1 if bootstrap draws were used for the bins |
| `e(min_e)` | value of mine() |
| `e(max_e)` | value of maxe() |
| `e(balance_e)` | value of balancee() |
| `e(na_rm)` | 1 with narm |

| Macros | Content |
| --- | --- |
| `e(title)` | title of the fit, followed by the aggregation type |
| `e(estat_cmdline)` | the estat command as typed |
| `e(agg_type)` | simple, group, dynamic, or calendar |
| `e(subcmd)` | the estat subcommand |
| `e(result_kind)` | aggregation |
| `e(overall_status)` | status of the overall estimate |
| `e(overall_reason)` | why it is unavailable, when it is |
| `e(overall_inference_status)` | status of its inference |
| `e(overall_inference_kind)` | type of its standard error |
| `e(band_family)` | type of band |
| `e(target_definition)` | definition of the aggregated targets |
| `e(weight_reduction)` | how unit weights were reduced |
| `e(weight_source)` | source of the weights |
| `e(selector_warning)` | note about mine(), maxe(), or balancee(), when any |

In chained mode, estat also posts the results of cdid that describe the data, the inference and the posting (`e(N)`, `e(N_full)`, `e(clustered)`, `e(N_clust)`, `e(level)`, `e(alpha)`, `e(critical_value)`, e(wald_*), `e(vce)`, `e(band)`, `e(rng)`, `e(b)`, `e(V)`, `e(se)`, `e(ci)`, `e(coef_map)`, ...) and the settings macros of the fit, except `e(anticipation)`, `e(inactive_options)`, `e(missing_reference)`, `e(reference_status)`, and `e(nuisance)`; `e(cmdline)` and `e(sample)` are those of the fit. In ordinary mode, estat posts the results listed for ordinary fits in help cdid, except `e(inactive_options)`, `e(V_did)`, and `e(V_did_asymptotic)`, and adds `e(estat_cmdline)`; `e(N_targets)` and `e(N_posted)` count the bins and the overall estimate, which is the last element of `e(b)`. The overall-effect scalars, `e(agg_type)`, and `e(subcmd)` are posted in chained mode only.

## Uncertainty

| Result | Content |
| --- | --- |
| `e(V)` | Analytical covariance of `e(b)`, from the unit influence functions. |
| `e(V_asymptotic)` | `e(N_full)` times `e(V)`. |
| `e(analytic_se)` | Square roots of the diagonal of `e(V)`. |
| `e(se)` | Reported standard errors: analytical with `bstrap(0)`, from the multiplier bootstrap (interquartile range of the draws) with `bstrap(1)`. |
| `e(ci)` | Reported intervals: simultaneous bands with `cband(1)`, pointwise intervals otherwise. |
| `e(ci_pointwise)` | Pointwise intervals. |
| `e(critical_value)` | Critical value of the reported intervals or bands. |

Because the reported standard errors can come from the bootstrap, `sqrt(e(V)[k,k])` need not equal `e(se)[1,k]`. The influence functions include the estimation of every fitted model (treatment, outcome, observation and imputation models, cohort shares). With `cluster()`, scores are summed within clusters. The multiplier bootstrap draws one Rademacher sign per cluster from Stata's random-number stream.

### Point estimates only

With `noinffunc`, `e(V)` and `e(V_asymptotic)` are not posted and `e(point_only)` is 1. With `posting(matrix)` and at least one available cell, `e(b)` is posted, `e(se)`, `e(analytic_se)`, `e(ci)` and `e(ci_pointwise)` are posted with missing entries, and `e(properties)` is `b`; with `posting(state)`, no matrix is posted and `e(properties)` is empty. `estat` is refused.

## Pretrend tests

The joint Wald test of the pretreatment cells is stored in `e(wald_statistic)`, `e(wald_df)` and `e(wald_p)`, with `e(wald_status)` and, when it is not computed, `e(wald_reason)`. It uses the analytical covariance; failure to reject does not establish parallel trends.

`cdid_pretest` returns the conditional test in `r()`:

| Scalars | Content |
| --- | --- |
| `r(CvM)` | Cramér-von Mises statistic |
| `r(CvMcval)` | bootstrap critical value at level() |
| `r(CvMpval)` | bootstrap p-value (share of draws strictly above `r(CvM)`) |
| `r(N)` | number of units |
| `r(M)` | number of evaluation points (one per unit in chained mode) |
| `r(K)` | number of pretreatment cells tested |
| `r(C)` | number of clusters |
| `r(B)` | number of bootstrap draws |
| `r(alpha)` | 1 - level/100 |
| `r(N_input)` | observations supplied |
| `r(N_requested)` | observations requested |
| `r(N_rows)` | rows used |
| `r(N_distinct)` | distinct units |
| `r(multipliers_used)` | 1 if bootstrap multipliers were drawn or supplied |
| `r(bstrap_requested)` | value of bstrap() |
| `r(cband_requested)` | value of cband() |
| `r(n_callbacks)` | number of callback calls, with method(custom) |
| `r(allowunbalanced_supplied)` | 1 if allowunbalanced() was specified |
| `r(rcs_dr_corrected_if)` | 1 when the ordinary doubly robust test with repeated cross-sections uses the corrected influence function |
| `r(formula_defaulted)` | 1 when the ordinary test used the intercept-only model by default |
| `r(id_defaulted)` | 1 when repeated cross-sections used observation identifiers by default |

| Macros | Content |
| --- | --- |
| `r(status)` | available, unavailable, or not_applicable |
| `r(reason)` | reason, when the test is not available |
| `r(mode)` | chained or ordinary |
| `r(method)` | treatment adjustment |
| `r(formula)` | covariates |
| `r(X_columns)` | columns of the treatment design |
| `r(Z_columns)` | columns of the evaluation grid |
| `r(outcome_pooling)` | outcome-regression pooling, when active |
| `r(clustervars)` | cluster variable, with cluster() |
| `r(callback)` | callback(), with method(custom) |
| `r(payload)` | payload(), with method(custom) |
| `r(allowunbalanced_requested)` | value of allowunbalanced(), when specified |
| `r(inactive_options)` | accepted options that have no effect in ordinary mode, specified or not |
| `r(posting)` | matrix or state |
| `r(state_token)` | name of the store() object, with store() |
| `r(rng)` | source of the multipliers, in words |
| `r(description)` | construction of the test, in words |
| `r(grid_rule)` | evaluation grid, in words |
| `r(influence_function)` | influence function used, in words |

| Matrices | Content |
| --- | --- |
| `r(CvMb)` | B x 1 bootstrap draws of the statistic |
| `r(family)` | tested pretreatment cells |
| `r(grid)` | evaluation grid: first_grid_row and multiplicity |

When the test is unavailable or not applicable, the statistics are missing and `r(CvMb)` is not returned. cdid_pretest_load returns the same results, except `r(callback)`, `r(payload)`, `r(n_callbacks)`, and `r(inactive_options)`, and adds `r(state_file)`.

## Results of cdid with dvar()

`cdid` with `dvar()` stores the following in `e()`; `estat` does not apply.

| Scalars | Content |
| --- | --- |
| `e(N)` | number of units |
| `e(T)` | number of dates |
| `e(N_clust)` | number of clusters, with cluster() |
| `e(att)` | switcher-weighted average effect |
| `e(se_att)` | its reported standard error |
| `e(se_att_analytic)` | its analytical standard error |
| `e(se_att_bootstrap)` | its bootstrap standard error |
| `e(n_switchers)` | number of switchers in the requested cells |
| `e(N_missing)` | number of missing outcomes |
| `e(critical_value)` | critical value of the reported intervals or band |
| `e(level)` | confidence level |
| `e(alpha)` | 1 - `e(level)`/100 |
| `e(biters)` | bootstrap draws requested (biters()) |
| `e(draws)` | 1 if bootstrap draws were used |

| Macros | Content |
| --- | --- |
| `e(cmd)` | cdid |
| `e(cmdline)` | command as typed |
| `e(title)` | title of the output |
| `e(design)` | static |
| `e(mode)` | chained |
| `e(depvar)` | outcome |
| `e(ivar)` | unit identifier |
| `e(time)` | date variable |
| `e(dvar)` | treatment status variable |
| `e(covariates)` | covariates |
| `e(method)` | treatment adjustment |
| `e(switchers)` | both, entry, or exit |
| `e(switchtimes)` | dates of switchtimes(), when specified |
| `e(outcome_pooling)` | outcome-regression pooling |
| `e(gmm_weighting)` | GMM weighting, with missing outcomes or selection() |
| `e(selection)` | assumption about missing outcomes |
| `e(selection_est_method)` | correction for missing outcomes, when it applies |
| `e(selection_model)` | observation model, when it applies |
| `e(selection_pooling)` | observation-model pooling, when it applies |
| `e(selection_universe)` | units that fit the observation model, when it applies |
| `e(selection_rho)` | density-ratio variable, when specified |
| `e(selection_justification)` | text of seljustification(), when specified |
| `e(selection_covariates)` | observation-model covariates, when they apply |
| `e(imputation_pooling)` | imputation pooling, when it applies |
| `e(imputation_scale)` | imputation scaling, when it applies |
| `e(vce)` | analytic, bootstrap, or none |
| `e(vcetype)` | Bootstrap, with the multiplier bootstrap |
| `e(band)` | pointwise or simultaneous |
| `e(band_status)` | ok, not_requested, unavailable_target_cell, or unavailable_degenerate_family |
| `e(clustvar)` | cluster variable, with cluster() |
| `e(nuisance)` | availability of the reports kept, with keepnuisance |
| `e(state_token)` | reference to the kept reports, with keepnuisance |
| `e(estat_cmd)` | cdid_estat |
| `e(properties)` | b V, or b with noinffunc |

| Matrices | Content |
| --- | --- |
| `e(b)` | estimable cell effects, named direction:tdate (for example `_b[entry:t2003]`) |
| `e(V)` | analytical covariance of `e(b)`, unless noinffunc |
| `e(cells)` | one row per requested cell; see the column definitions below |

| Functions | Content |
| --- | --- |
| `e(sample)` | marks the observations used |

`e(b)` and `e(V)` are not posted when no requested cell is available. The columns of `e(cells)` are: time, the date t; tprev, the preceding study date; d, the treatment status before the switch (0 for entries, 1 for exits); sign, the sign applied to the switcher-stayer contrast (1 for entries, -1 for exits); available, 1 if the cell is estimated; switchers and stayers, the numbers of units; att, the cell effect; theta, the switcher-stayer contrast before the sign is applied; se, lower and upper, the reported standard error and interval limits; analytic_se and bootstrap_se; and status, 0 when the cell is available, 1 with no switchers, 2 with no stayers, 3 when the observed dates of the switchers do not connect the preceding date and t, and 4 when those of the stayers do not.

## Chain and model reports

`cdid_chain` returns:

| Scalars | Content |
| --- | --- |
| `r(n_cells)` | number of cells reported |
| `r(n_available)` | number of available cells |
| `r(n_chained)` | number of available cells that require chaining |
| `r(n_moments)` | number of comparisons used by GMM |
| `r(overid_df)` | number of overidentifying restrictions |

| Macros | Content |
| --- | --- |
| `r(kind)` | crosspair or samepair |

`cdid_nuisance` returns the following; the block results refer to the selected block:

| Scalars | Content |
| --- | --- |
| `r(N)` | number of units |
| `r(n_blocks)` | number of blocks |
| `r(n_units)` | units in the selected block |
| `r(n_fitting_rows)` | fitting rows of the selected block |
| `r(n_comparisons)` | comparisons that use the selected block |

| Macros | Content |
| --- | --- |
| `r(component)` | component displayed |
| `r(status)` | available, not_applicable, not_captured, or unavailable |
| `r(reason)` | why a component is not available |
| `r(status_component)` | status of each component, for example `r(status_selection)` |
| `r(method)` | treatment adjustment |
| `r(selmethod)` | correction for missing outcomes |
| `r(block)` | identifier of the selected block |
| `r(label)` | label of the selected block |
| `r(equation)` | equation of the selected block |
| `r(terms)` | labels of the coefficients of the selected block |
| `r(kept_slopes)` | slopes kept, after the preset |
| `r(aliased_slopes)` | slopes dropped as exact aliases, after the preset |
| `r(roster_terms)` | labels of `r(roster)` |

| Matrices | Content |
| --- | --- |
| `r(blocks)` | one row per block: parameters, fitting rows, units or weights, comparisons, and related counts |
| `r(table)` | coefficients of the selected block, with standard errors, z statistics, p-values, and 95% confidence limits |
| `r(roster)` | roster probabilities, after the preset with seluniverse(A) or seluniverse(B) |

## Diagnostics

| Scalars | Content |
| --- | --- |
| `r(N_units)` | number of units |
| `r(N_treated)` | treated units |
| `r(N_never_treated)` | never-treated units |
| `r(N_unknown_cohort)` | units with unknown cohort |
| `r(N_with_rows)` | units with data rows |
| `r(N_frame_only)` | units supplied only through frame() |
| `r(N_with_outcome)` | units with at least one observed outcome |
| `r(N_with_two_outcomes)` | units with at least two observed outcomes |
| `r(N_cluster_ids)` | distinct cluster identifiers retained from the fit |
| `r(N_clust)` | clusters used for cell inference, missing when unclustered |
| `r(weight_sum)` | sum of unit weights |
| `r(weight_min)` | smallest unit weight |
| `r(weight_max)` | largest unit weight |
| `r(N_cells)` | number of cells |
| `r(N_post)` | post-treatment cells |
| `r(N_pre)` | pretreatment cells |
| `r(N_reference)` | reference cells |
| `r(N_available)` | available cells |
| `r(N_unavailable)` | unavailable cells |
| `r(N_with_se)` | cells with a standard error |
| `r(B_requested)` | number of multiplier draws requested |
| `r(wald_statistic)` | pretrend Wald statistic |
| `r(wald_df)` | its degrees of freedom |
| `r(wald_p)` | its p-value |
| `r(N_aggregates)` | aggregation bins, after estat |
| `r(N_aggregates_available)` | available bins, after estat |

| Macros | Content |
| --- | --- |
| `r(mode)` | chained or ordinary |
| `r(result_kind)` | cells or aggregation |
| `r(agg_type)` | aggregation type, after estat |
| `r(method)` | treatment adjustment |
| `r(control_group)` | comparison group |
| `r(base_period)` | base period |
| `r(pt_assumption)` | parallel trends assumption |
| `r(gmm_weighting)` | GMM weighting |
| `r(selection)` | assumption about missing outcomes |
| `r(wald_status)` | status of the pretrend Wald test |
| `r(inference_status)` | status of the inference |
| `r(band)` | intervals for the cells: pointwise, simultaneous, or none |
| `r(vce)` | standard errors for the cells: analytic, bootstrap, or none |
| `r(clustvar)` | cluster variable retained from the ordinary fit |
| `r(rng)` | source of multiplier draws, in chained mode |
| `r(inactive_options)` | inactive options in the chained specification |
| `r(wald_reason)` | cells used by the test or why it was unavailable |
| `r(overall_status)` | availability of the overall summary, after estat |
| `r(overall_reason)` | reason for an unavailable overall summary, after estat |
| `r(cells_matrix)` | note posted instead of `r(cells)` when the table exceeds c(max_matdim); use cdid_export |
| `r(aggregate_matrix)` | the same for `r(aggregate)`, after estat |

| Matrices | Content |
| --- | --- |
| `r(cells)` | cells: cohort, date, role, availability, estimate, standard error, and bounds |
| `r(aggregate)` | aggregation bins, after estat |

cdid_diagnostics returns every displayed quantity in r(); the lists above are the main results. A missing numeric quantity is accompanied by `r(name_status)`, which explains why. The matrices are returned only when they fit within c(max_matdim).

## Saving and exporting

`saving()` writes the complete results of `cdid` with `gvar()`, of `estat` after such a fit, or of `cdid_pretest` (it is not available with `dvar()`); `cdid_load` and `cdid_pretest_load` restore them, with `e(state_file)` or `r(state_file)` naming the file. `cdid_export` writes the cells, the aggregation, the units, the influence functions or a summary as Stata datasets; see `help cdid_export` for their variables.

---

**Reading path, step 19 of 22** · Previous: [Options reference](options.md) · Next: [Diagnose errors and unavailable results](../help/troubleshooting.md)

[Guide home](../README.md) · [Package overview](../../README.md)
