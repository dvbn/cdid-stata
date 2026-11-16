# Diagnose errors and unavailable results

`cdid` checks the syntax and the data first, then estimates the models used for adjustment, determines which effects are identified, and finally computes the inference. A refused request stops with a message that names the condition that failed. An effect that the data cannot identify does not stop the command: it is reported as unavailable with a reason. A successful estimation does not show that its identifying assumptions hold.

The examples on this page use the data of the [hands-on do-file](../hands-on.md): run its Sections 0 and 1 first.

## Read the error

| Return code | Usual meaning |
| --- | --- |
| 198 | Invalid syntax or an option that does not apply, for example an option of `gvar()` estimation used with `dvar()`. |
| 459 | Data or a combination of options that the estimator cannot use: duplicate unit-date rows, weights that vary within units without `fixweights()`, an incomplete specification, a rank-deficient design, or a covariate missing at a reference date under `missingreference(error)`. |
| 498 | Results that do not fit the command, for example `cdid_chain` when the current results are not a chained fit, or a file of the wrong type. |
| 908 | A limit set by `workspace()`, `chunk()` or `workcap()`, or by Stata's matrix size, would be exceeded. |

Two incomplete specifications, as in Section 11 of the hands-on do-file:

```stata
capture noisily cdid y x_bin, ivar(id) time(period) gvar(first_treat) ptassumption(all)
display _rc
capture noisily cdid y x_bin, ivar(id) time(period) gvar(first_treat) selection(none) selmethod(dr)
display _rc
```

```text
ptassumption(all) requires ptstart(#), which must be one of the periods in the data; ptstart() is not allowed with ptassumption(post).
```

```text
selmethod(dr) requires selection(conditional_mtar).
```

| Problem | What to check |
| --- | --- |
| Invalid combination of options or invalid data | The mode, the missing-outcome options, the `frame()` options, `ptstart()` with `ptassumption(all)`, and `fixweights()` when weights vary over time. The [options reference](../reference/options.md) lists the estimators that accept each option. |
| No results, or results of another command | Rerun `cdid`, restore stored results with `estimates restore`, or load saved results with `cdid_load`. Datasets written by `cdid_export` cannot be loaded as results. |
| Existing output file | `cdid_export`, `cdid_plot` and `cdid_chain` accept `replace`; the `saving()` option of `cdid`, `estat` and `cdid_pretest` never overwrites a file, so choose a new name. The `.dta` rule applies to `cdid_export` and to `cdid_plot, plotdata()`: a name with another extension is refused. `cdid_plot, saving()` writes a `.gph` file and `export()` the graph formats of `graph export`. Names created by `tempfile` are accepted by all of them. |
| Memory or size limit | `workspace()`, `chunk()` and `workcap()` concern the main numerical arrays, not the total memory or time used by Stata, and do not relax any check of support or convergence. `posting(state)` avoids Stata's matrix size limit. |

## Inspect a completed fit

`cdid_diagnostics` summarizes the units, the cells by role and availability, the GMM weighting, the settings and the inference; `cdid_diagnostics, detail` also lists each unavailable cell with its reason. `cdid_chain` shows which comparisons identify each available cell, and `cdid_nuisance` the fitted models after `keepnuisance`.

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) bstrap(0) cband(0)
cdid_diagnostics, detail
cdid_chain
```

## An effect is unavailable

The results keep every requested cell; unavailable cells have missing estimates and inference and a stated reason. With `posting(matrix)`, only the available estimates are posted in `e(b)`, and `e(coef_map)` gives the cell of each coefficient. `cdid_export, type(cells)` saves all cells with their reasons. Common reasons are that no comparison unit is observed at the dates involved, that a cohort g has no study date d with d + anticipation < g, and hence no reference date, that the observed changes of a cohort do not connect its reference date to the target date, or, under `missingreference(exclude)`, that no eligible unit remains. Regularizing the GMM covariance cannot create a missing comparison. Remedies, when their assumptions are acceptable, are not-yet-treated controls, `crosspair` for cohorts observed over different windows, or `ptassumption(all)` or `ptassumption(ctrl)`; see [Cross-pair chaining](../methods/cross-pair.md).

`narm` in `estat` changes how an aggregation treats unavailable cells; it does not identify the missing effects. Check which cells were included and excluded in the dataset written by `cdid_export, type(aggregate)` before interpreting the aggregation.

## A logit model has no finite solution

An observation logit with a free cell intercept has no finite estimate in a cell where every unit, or no unit, is observed; a treatment logit has none when a covariate value perfectly predicts cohort membership within a comparison. Pooling the slopes keeps the cell intercepts and does not solve the first problem. Pure imputation (`selmethod(imputation)`) fits no observation logit, so it can be used when every outcome of a cell is observed and imputation is justified; it still needs enough observed changes and support for prediction. `cdid` does not drop observations, change the pooling or switch estimators on its own. That another correction runs on the same data does not fix the original specification.

## An outcome regression or imputation design is refused

A regression design that is rank deficient on its fitting rows is refused with error 459 and a message naming the regression. When the redundant columns are exactly collinear on every row used to fit and to predict (for example a covariate and a copy of it), the redundant column is dropped, as in R. Remove the redundant covariate, coarsen a factor variable, or use a coarser `imppooling()`.

## Missing covariates or weights versus missing outcomes

A missing outcome is an expected input. A missing covariate or weight is not: time-invariant covariates must be complete for the units used, and covariates or weights that vary over time are needed at the dates where the comparisons read them. Under the default `missingreference(exclude)`, a covariate missing at a cohort's reference date excludes the unit from that cohort's comparisons with a warning; under `missingreference(error)`, the command stops and lists the missing values. See [Prepare your panel and unit table](prepare-data.md#reference-dates-for-covariates-that-change-within-a-unit).

## Pretrend tests are unavailable

The Wald test computed by `cdid` requires every pretreatment cell and a valid joint covariance. Centered sums over $C$ clusters have rank at most $C - 1$, so too many restrictions, a singular covariance or a missing cell make the test unavailable; `e(wald_status)` and `e(wald_reason)` say why. Under `missingreference(exclude)`, an unavailable pretrend comparison also makes the joint test unavailable. A subset of cells is never tested instead.

The conditional test of `cdid_pretest` needs every pretreatment cell of the request and the covariates of every unit of the study population at the relevant dates; it reports `r(status)` and `r(reason)`. It stops with error 459 when a required covariate is missing.

## Standard errors look unusual

`e(se)` holds the reported standard errors: with `bstrap(1)`, the default, they come from the multiplier bootstrap and differ from the analytical `e(analytic_se)`, and they change with the random-number stream. With few clusters, bootstrap and analytical standard errors can be unstable. Simultaneous bands are generally wider than pointwise intervals, but the bootstrap critical value of a band need not exceed the normal critical value of a pointwise interval. See [Inference](../examples/inference.md).

## Results differ from ordinary mode or from csdid

Chained and ordinary estimates agree only in the settings described in [Moving from did and csdid](../from-did.md#when-chained-and-ordinary-estimates-should-agree). Ordinary mode drops the units that are not observed at every date unless `allowunbalanced` is specified, and options specific to chained mode have no effect there.

## Commands not available after some results

```stata
quietly cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) noinffunc
capture noisily estat simple
display _rc
```

```text
estat aggregations need influence functions; refit cdid without noinffunc.
```

`estat` needs the influence functions that `noinffunc` omits. After `estat`, `cdid_chain` reports the comparisons of the fit computed in the same session, but `cdid_nuisance` needs the fit: restore it with `estimates restore`. After `dvar()`, `cdid_plot`, `cdid_chain` and `cdid_diagnostics` are not available, and `estat` stops with a message that points to `e(att)` and `e(se_att)`.

---

**Reading path, step 20 of 22** · Previous: [Results reference](../reference/results.md) · Next: [Cross-pair chaining, control parallel trends and within-unit weighting](../methods/cross-pair.md)

[Guide home](../README.md) · [Package overview](../../README.md)
