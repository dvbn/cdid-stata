# Moving from R cdid

The original Stata port followed [R `cdid` 1.3.0.0](https://github.com/dvbn/cdid/tree/v1.3.0.0). The chained cohort weighting defaults follow [R `cdid` 1.4.0.0](https://github.com/dvbn/cdid/commit/5f554b848038c9613208164032794ae19e7e5d04). The [differences below](#differences-to-keep-in-mind) describe the behavior of each interface. This page maps R functions and arguments to Stata commands and options. "Corresponds" means the same estimator with the same defaults. Estimates and standard errors were compared with R on the configurations of the Stata test suite and agree there up to numerical precision; agreement is not established beyond them, and multiplier-bootstrap draws come from Stata's random-number generator and differ from R's.

## Functions

| R `cdid` | Stata `cdid` |
| --- | --- |
| `att_gt(gname = ...)` | `cdid depvar [indepvars], ivar() time() gvar()` |
| `att_gt(dname = ...)` | `cdid depvar [indepvars], ivar() time() dvar()`; see `help cdid_dvar` |
| `summary(fit)`, `print(fit)` | `cdid` (replay) |
| `aggte(fit, type = "simple"\|"group"\|"dynamic"\|"calendar")` | `estat simple\|group\|dynamic\|calendar` |
| `ggdid()` | `cdid_plot` (`twoway` options in place of the `ggplot` arguments) |
| `tidy()`, `glance()` | `cdid_export, type(cells)` or `type(aggregate)`, `type(summary)`, and `e()` |
| `conditional_did_pretest()` | `cdid_pretest` (results in `r()`) |
| `chain_report(fit, cohort, time, include_pre, tol)` (`group` and `t` are aliases of `cohort` and `time`); its print method with `show_gmm_weights = TRUE` | `cdid_chain, cohort() time() pre tol()`; `weights` |
| `as.data.frame(report, what = "summary"\|"witness"\|"gmm_weights")` | `cdid_chain, saving() table(summary\|witness\|weights)` |
| `nuisance_models(fit, component, block)` | `cdid_nuisance component, block()` |
| `saveRDS()`, `readRDS()` | `saving()`, `cdid_load`, `cdid_pretest_load`, for `gvar()` fits, their aggregations and pretests; `dvar()` results cannot be saved |
| `fit$aggregate` (switching model) | `e(att)`, `e(se_att)` after `dvar()` |

## Arguments of att_gt() with gname

| R argument | Stata option |
| --- | --- |
| `yname` | `depvar` |
| `tname` | `time()` |
| `idname` | `ivar()` |
| `gname` | `gvar()` |
| `xformla` | `indepvars` or `covars()`; `~ -1 + ...` is `noconstant` |
| `data` | the data in memory |
| `panel = FALSE` | `rcs` |
| `allow_unbalanced_panel` | `allowunbalanced` |
| `control_group` | `controlgroup()` |
| `anticipation` | `anticipation()` |
| `weightsname` | `[iweight = ...]` or `wvar()` |
| `fix_weights` | `fixweights()` |
| `alp` | `level()`, with `alp = 1 - level/100` |
| `na.rm` (in `aggte()`) | `narm` of `estat` |
| `bstrap`, `cband`, `biters` | `bstrap()`, `cband()`, `biters()` (0 or 1 for the first two) |
| `clustervars` | `cluster()` (one variable besides the unit identifier) |
| `est_method` | `method()`; an R function corresponds to `method(custom) callback()`, a Mata function, and the arguments passed through `...` to `payload()` |
| `base_period` | `baseperiod()` |
| `compute_inffunc = FALSE` | `noinffunc` |
| `mode` | `mode()` |
| `pt_assumption`, `pt_start` | `ptassumption()`, `ptstart()` |
| `gmm_weighting`, `gmm_covariance` | `gmmweighting()`, `gmmcovariance()` |
| `cross_pair = TRUE` | `crosspair` |
| `outcome_pooling` | `outcomepooling()` |
| `frame` | `frame()`, with `frameivar()`, `framegvar()`, `framecovars()`, `framewvar()`, `frameselcovars()`, `framecluster()` naming its columns |
| `selection` | `selection()` |
| `selection_est_method` | `selmethod()` |
| `selection_model` | `selmodel()` |
| `selection_pooling` | `selpooling()` |
| `selection_xformla` | `selcovars()`; `~ 1` is `selconstant`, `~ -1 + ...` is `selnoconstant`, `~ 0` is `selzero` |
| `selection_universe` | `seluniverse()` |
| `selection_rho` | `selrho()` |
| `selection_justification` | `seljustification()` |
| `imputation_pooling` | `imppooling()` |
| `missing_reference` | `missingreference()` |
| `control$preset` | `preset()` |
| `control$keep_nuisance` | `keepnuisance` |
| `control$numeric_workspace_bytes` | `workspace()` |
| `control$chunk_rows` | `chunk()` |
| `control$counted_work_limit` | `workcap()` |
| `control$condition_trigger` | `condition()` |
| `control$ridge_relative` | `ridge()` |

Values are spelled as in R (`"nevertreated"` is `nevertreated`), and a Stata option that is omitted takes the R default.

## Arguments of att_gt() with dname

| R argument | Stata option |
| --- | --- |
| `dname` | `dvar()` |
| `switchers`, `switch_times` | `switchers()`, `switchtimes()` |
| `outcome_pooling` | `outcomepooling()`, including blocks of cells |
| `imputation_scale` | `impscale()` |
| `selection_pooling`, `imputation_pooling` | `selpooling()`, `imppooling()` (`separate`, `within_group`, `pooled`) |
| other arguments | as with `gname`, where accepted; see `help cdid_dvar` |

## Arguments of aggte() and conditional_did_pretest()

| R argument | Stata option |
| --- | --- |
| `type` | the `estat` subcommand |
| `min_e`, `max_e`, `balance_e` | `mine()`, `maxe()`, `balancee()` |
| `na.rm` | `narm` |
| `bstrap`, `cband`, `biters`, `alp`, `clustervars` | `bstrap()`, `cband()`, `biters()`, `level()`, `cluster()` |
| `pl`, `cores` (pretest) | `pl()`, `cores()`; only `pl(0) cores(1)` is accepted, and the Stata test runs serially |

## R arguments without a Stata option

None of these arguments changes the estimator.

| R function | Argument | Reason |
| --- | --- | --- |
| `att_gt()` | `faster_mode` | Selects the computational path of `did` in ordinary mode; Stata's ordinary mode is computed by `csdid`. |
| `att_gt()`, `conditional_did_pretest()` | `print_details` | Progress printing. |
| `att_gt()`, `conditional_did_pretest()` | `pl`, `cores` | Parallel bootstrap. `cdid` has neither option; `cdid_pretest` accepts only `pl(0) cores(1)`. |
| `att_gt()`, `conditional_did_pretest()` | `control$workers` | Number of parallel workers in chained mode. |
| print method of `chain_report()` | `digits`, `max_cells` | Printing. |
| `ggdid()` | `ylim`, `xlab`, `ylab`, `title`, `xgap`, `ncol`, `legend`, `ref_line`, `theme_choice`, ... | Graphics; `cdid_plot` accepts `twoway` options instead. |
| `tidy()`, `glance()` | | Summaries; `cdid_export` writes the tables. |

## Stata options without an R argument

| Command | Option | Purpose |
| --- | --- | --- |
| `cdid` | `targets()`, `targetcohorts()` | Restrict the target cohorts; R estimates all of them. |
| `cdid`, `estat`, `cdid_pretest`, `cdid_load`, `cdid_pretest_load` | `posting(matrix\|state)`, `store()` | How results are posted in Stata, given its matrix size limit. |
| `cdid`, `estat`, `cdid_pretest` | `saving()` | Save the results to a file; R uses `saveRDS()`. |
| `cdid_pretest` | `signs()` | Supply the bootstrap multipliers. |
| `estat`, `cdid_load`, `cdid_pretest_load` | `workspace()`, `chunk()` (`estat` only) | Memory limits; R reuses the settings of the fit. |
| `cdid_chain` | `saving()`, `table()` | Save a chain table; R uses `as.data.frame()`. |
| `cdid_plot` | `group()`, `cells`, `saving()`, `export()`, `plotdata()`, `replace` | Graph control. |
| `cdid_export` | `type()`, `replace` | Tables as Stata datasets. |
| `cdid_diagnostics` | `detail` | Lists unavailable cells; in R, inspect the result object. |

## Chained cohort weighting

Built-in IPW, regression and DR default to `within` in both languages. For the earlier default, use `gmm_weighting = "optimal"` in R or `gmmweighting(optimal)` in Stata. Switching fits, custom estimators and the IPW and DR branches of the cohort fixed-effects preset default to `optimal`. The regression branch of the preset uses within-unit least squares. See [GMM weighting](methods/gmm.md) for details.

## Differences to keep in mind

- **Ordinary mode.** R `cdid` delegates ordinary mode to `did` (version 2.5.1 in the comparisons). Stata `cdid` uses `csdid` 2.0.2, installed separately. The ordinary fit lists all accepted inactive options in `e(inactive_options)`, including `selection()`, the `sel*()` options and `keepnuisance`. Ordinary `cdid` refuses `covars()`, `noconstant`, `panel()`, `frame()` and its associated options. Ordinary `cdid_pretest` uses `covars()`, `noconstant` and `panel()`, and treats `frame()` and its associated options as inactive. Its `r(inactive_options)` lists all accepted inactive options. Each command displays a note naming those supplied.
- **Clustering of aggregations.** `estat` inherits the clustering of the fit. After an ordinary fit, `cluster(none)` or the `ivar()` variable requests unit-level inference. A request for an unavailable cluster variable warns and uses unit-level inference for both analytical and bootstrap standard errors, as in R's ordinary `aggte()`. After a chained fit, a request to change clustering stops with an error.
- **Printed output.** `cdid_pretest` displays the test, as R's print method; its results are in `r()`.
- **Results.** Stata results are posted in `e()` and `r()` rather than in a fitted object; [Results](reference/results.md) lists them, and `cdid_export` writes the tables that `tidy()` returns in R.
- **Factor covariates.** Factor-variable notation (`i.`, `ib#.`, `#`, `##`) corresponds to R's `factor()` treatment contrasts with the same base. A factor interacted with a continuous covariate without the main effect, `i.cf#c.x` or `ibn.cf#c.x`, keeps one slope per level, as R's `factor(cf):x`; with the main effect of the continuous covariate, `c.x i.cf#c.x` drops the base level, as R's `x + factor(cf):x`.
- **Comparisons with other software.** [Estimate familiar DiD models with cdid](from-other-packages.md) compares with Stata commands where the R guide compares with R packages, and reports a numerical match only where both sides were run in Stata.
