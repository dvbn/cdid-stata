# Hands-on guide to cdid in Stata

This guide walks through a complete analysis with `cdid`, using the do-file `examples/cdid_hands_on.do`. Run the do-file from top to bottom, or run Section 0 and Section 1, which creates the data, and then run later sections from the Do-file Editor. Later sections use the simulated data already in memory.

- **Stata package repository**: <https://github.com/dvbn/cdid-stata>
- **R package**: <https://github.com/dvbn/cdid>

---

## 0. Settings

Before running `examples/cdid_hands_on.do`, start a new Stata session or save the data in memory. The settings are global macros; values you set before running the do-file are kept. Section 0 gives each unset `CDID_DEMO_*` global its default value and defines a small program, `cdid_demo_init`, so that individual sections can be run from the Do-file Editor with your settings:

```stata
global CDID_DEMO_CLEAR_OK      0          // Set to 1 to allow Section 1 to clear memory
global CDID_DEMO_N             800        // Number of units (600-1000 recommended)
global CDID_DEMO_T             6          // Number of calendar periods (5 or 6 recommended)
global CDID_DEMO_SEED          20261103   // Random-number seed for the simulated data
global CDID_DEMO_BITERS        199        // Bootstrap draws (package default is 1000)
global CDID_DEMO_MISS_LAW      "conditional_mtar" // "mcar", "conditional_mtar", or "rotating_window"
global CDID_DEMO_MISS_STRENGTH 0.85       // Strength of covariate-driven observation log-odds
global CDID_DEMO_TREND_HET     0.35       // Covariate slope in untreated outcome changes
global CDID_DEMO_EFF_BASE      1.20       // Baseline post-treatment causal effect
global CDID_DEMO_EFF_SLOPE     0.30       // Effect growth per event-time period after treatment
global CDID_DEMO_EFF_HET       0.25       // Effect heterogeneity by binary covariate x_bin
global CDID_DEMO_RUN_ADVANCED  0          // Set to 1 to run optional advanced comparisons
global CDID_DEMO_OUTDIR        ""         // Folder for graphs, datasets and saved results ("" = temporary files only)
global CDID_DEMO_REPLACE       0          // Set to 1 to allow cdid_plot/cdid_export to overwrite files
```

### Conventions
1. **Bootstrap draw count**: The package default in both `cdid` and `cdid_pretest` is `biters(1000)`. The do-file sets `CDID_DEMO_BITERS = 199` so that the examples run quickly.
2. **Memory and files**: Section 1 clears the data in memory only if `CDID_DEMO_CLEAR_OK` is 1. Set `CDID_DEMO_OUTDIR` to keep graphs, exported datasets and saved results in a folder. Without it, the do-file exports datasets and saves results to temporary files and skips graph saving. Existing graphs and datasets are replaced only if `CDID_DEMO_REPLACE` is 1.

---

## 1. Data requirements and the simulated panel

`cdid` expects a long panel dataset where each row represents a unique `(ivar, time)` pair:

- **Unit identifier (`ivar(varname)`)**: Numeric or string identifier uniquely identifying each unit in each period.
- **Calendar time (`time(varname)`)**: Numeric calendar date (for example, `1, 2, ..., 6` or `2018, 2019, ..., 2023`).
- **Treatment cohort (`gvar(varname)`)**: Time-constant numeric first-treatment date for each unit. Never-treated control units must be coded as `0`; treated units carry the numeric calendar period in which treatment begins (for example, `3` or `5`).
- **Missing outcomes in long data**: Every unit can be represented in the long dataset across all dates even when the outcome `y` is missing (`y == .`). Rows with missing `y` keep their unit identifier, calendar date, cohort `gvar()`, covariates, cluster identifier, and sampling weights so `cdid` can use full unit-level covariate and observation-history information.
- **Separate unit-level population table (`frame()`)**: When some units in the target population have zero observed outcomes (or do not appear in the outcome panel at all), you can supply a one-row-per-unit Stata frame via `frame(name)` along with `frameivar()`, `framegvar()`, `framecovars()`, `frameselcovars()`, `framewvar()`, and `framecluster()`.

### The simulated data
In Section 1 of `examples/cdid_hands_on.do`, we generate $N = 800$ units over $T = 6$ periods with staggered cohorts `first_treat` $\in \{0, 3, 5\}$, unit fixed effects $\alpha_i$, binary covariate `x_bin`, continuous covariate `x_cont`, and selection predictor `z_sel`. Untreated potential outcomes are generated from their **changes**:

$$\Delta Y_{it}(0) = Y_{it}(0) - Y_{i,t-1}(0) = \delta_t + \theta_1 \text{x\_bin}_i + \theta_2 \text{x\_cont}_i + (\varepsilon_{it} - \varepsilon_{i,t-1})$$

Because chained DiD operates on within-unit outcome differences $Y_{it} - Y_{is}$, the unit fixed effects $\alpha_i$ in outcome levels cancel. Identification does not require observed and unobserved units to have equal outcome levels.

---

## 2. Units observed at both dates, the study population, and missingness checks

### Units observed at both dates and the study population
To estimate the group-time average treatment effect $ATT(g,t)$ through a chain of transitions $(s, t)$, `cdid` forms elementary within-unit differences $\Delta Y_{i,s,t} = Y_{it} - Y_{is}$ among the units in cohort $g$ (or the control group) whose outcomes are observed at **both** dates $s$ and $t$ ($R_{is} R_{it} = 1$).

- **Mean-change representativeness**: Under the stated conditioning (unconditional in baseline `selection(none)`, or conditional on the specified covariates in `selection(conditional_mtar)`), **units observed at both endpoints represent the cohort with respect to the mean outcome change** between $s$ and $t$.
- **Study population**: although each comparison $(s,t)$ is estimated from the units observed at both $s$ and $t$, cohort shares and aggregation weights are computed over the study population of $N$ units (all units of the unit table), not over the units observed at both dates.

### Practical Missingness Checks Before Estimation
Before fitting `cdid`, always inspect observation patterns in the long panel (Section 2):

```stata
* 2a. Count observed outcome dates per unit across all units in the long dataset
tempvar n_obs_unit tag_unit
bysort id: egen long `n_obs_unit' = count(y)
bysort id: generate byte `tag_unit' = (_n == 1)
tabulate `n_obs_unit' first_treat if `tag_unit', missing

* 2b. Inspect outcome observation rates by calendar period and cohort
tabulate period first_treat, summarize(obs_y) means

* 2c. Compare covariate composition across all units vs repeat-observer units (>=2 dates)
summarize x_bin x_cont z_sel if `tag_unit'
summarize x_bin x_cont z_sel if `tag_unit' & `n_obs_unit' >= 2
```

---

## 3. Chained estimation: control groups, base periods, parallel trends and treatment adjustment

Section 3 demonstrates the core model specification options of `cdid`:

```stata
* 3a. Doubly robust (DR) chained estimation with never-treated controls and varying base period
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) controlgroup(nevertreated) baseperiod(varying) ///
    biters($CDID_DEMO_BITERS)

* Redisplay the results without refitting
cdid

* 3b. Not-yet-treated control group and universal base period
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) controlgroup(notyettreated) baseperiod(universal) ///
    biters($CDID_DEMO_BITERS)

* 3c. Parallel trends over all periods from `ptstart(1)`
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) controlgroup(nevertreated) ptassumption(all) ptstart(1) ///
    biters($CDID_DEMO_BITERS)
```

### Key Choices
- **`controlgroup(nevertreated|notyettreated)`**: `nevertreated` uses units with `first_treat == 0`. `notyettreated` expands the comparison group for transition $(s,t)$ to include cohorts not yet treated by $\max(s,t) + \text{anticipation}$.
- **`baseperiod(varying|universal)`**: Controls the reporting baseline for pre-treatment cells (`varying` uses adjacent pre-treatment transitions; `universal` anchors to the reference date $r(g)$, the latest study date before $g - \text{anticipation}$, which is $g - 1 - \text{anticipation}$ when the dates are consecutive integers).
- **`ptassumption(post|all)`**: `baseperiod()` and `ptassumption()` are independent settings. `ptassumption(post)` (default) enforces parallel trends only from $r(g)$ onward, leaving pre-treatment differences unrestricted for pretrend testing. `ptassumption(all)` imposes parallel trends from the common calendar date given in `ptstart(#)`, adding comparisons before the reference period to the GMM system.
- **`method(dr|reg|ipw)` and `outcomepooling(separate|pooled)`**: Selects how treatment nuisance models adjust for covariates $X$:
  - `method(ipw)` reweights control outcome changes using propensity scores (`outcomepooling()` does not apply).
  - `method(reg)` fits linear regressions for untreated **outcome changes** $\Delta Y_{i,s,t}$, fitted separately per comparison block by default (`outcomepooling(separate)`) or with shared duration-scaled covariate slopes across blocks (`outcomepooling(pooled)`). It models outcome changes, not outcome levels.
  - `method(dr)` (default) combines both for doubly robust treatment adjustment.

---

### Compare GMM weighting on the same panel

Built-in chained cohort fits default to `gmmweighting(within)`. The working covariance uses the dates and unit weights of the comparisons. Run the same fit with `gmmweighting(optimal)` to reproduce the earlier default and compare the ATT and analytical standard errors. Both fits compute uncertainty from the full study-unit influence functions.

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) bstrap(0) cband(0)
display e(gmm_weighting)

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    gmmweighting(optimal) bstrap(0) cband(0)
display e(gmm_weighting)
```

## 4. Observation-Model Selection (`selection(conditional_mtar)`), Logit Specifications, and Poolings

When outcome missingness depends on covariates that also drive untreated outcome changes, `selection(none)` can be biased: the units observed at both dates then over-represent strata with steeper or flatter trends, so their mean change no longer represents that of the cohort. `selection(conditional_mtar)` adjusts for this with a model of observation.

### Direct Pair Logit vs. Independent Period Logit
Under `selmethod(ipw)` (the default), `cdid` estimates pair observation probabilities $q_{g,s,t}(Z) = \Pr(R_{is}=1, R_{it}=1 \mid G_i=g, Z_i)$ using one of two logit families:

1. **`selmodel(pair_logit)` (default)**: Directly models the joint indicator $R_{is} R_{it}$ for each pair $(s,t)$.
2. **`selmodel(independent_period_logit)`**: Models marginal date observation probabilities $\pi_{g,t}(Z) = \Pr(R_{it}=1 \mid G_i=g, Z_i)$ and forms $q_{g,s,t}(Z) = \pi_{g,s}(Z)\,\pi_{g,t}(Z)$.

**When both agree vs. when `independent_period_logit` fails**:
- If observation indicators are independent across dates conditional on a binary predictor and each probability model is saturated in that predictor, the two specifications can both be correct. The simulated data also make observation depend on the continuous `z_sel` and include units observed once or never. Section 4a omits `z_sel`, so it compares two specifications without implying that either is correct. A product of date logits need not have a linear pair-logit form.
- **Do not use `selmodel(independent_period_logit)`** when units are observed in rotating two-date windows (`CDID_DEMO_MISS_LAW = "rotating_window"`), staggered entry/exit panels, or absorbing dropout. In those designs, observing a unit at date $s$ strongly predicts whether it can be observed at date $t$, violating conditional independence between dates even when marginal date probabilities are well-behaved. Use `selmodel(pair_logit)` instead.

### Selection Poolings (`selpooling`) and Separate Covariates (`covars` vs. `selcovars`)
- **`selpooling(separate|within_cohort|across_cohorts)`**:
  - `separate` (default) fits unrestricted parameters for each cohort and pair (or date).
  - `within_cohort` and `across_cohorts` pool slope parameters across transitions or across cohorts. More pooled $q$ models impose additional parametric restrictions to stabilize estimation in sparse cells; they do **not** automatically produce better estimates when cell-level support is already sufficient.
- **Separate outcome and selection covariates**: Use `covars()` for covariates $X$ entering the treatment/outcome-change model and `selcovars()` for covariates $Z$ entering the observation model (along with `framecovars()` and `frameselcovars()` when using `frame()`):

```stata
cdid y, ivar(id) time(period) gvar(first_treat) ///
    covars(i.x_bin c.x_cont) ///
    selection(conditional_mtar) selmodel(pair_logit) selpooling(separate) ///
    selcovars(i.x_bin c.x_cont c.z_sel) ///
    frame(cdid_pop) frameivar(id) framegvar(first_treat) ///
    framecovars(x_bin x_cont) frameselcovars(x_bin x_cont z_sel) ///
    biters($CDID_DEMO_BITERS)
```

---

## 5. Pure Outcome-Change Imputation vs. Augmented Correction (`selmethod`)

In addition to inverse pair-observation probability weighting (`selmethod(ipw)`), `cdid` supports two regression-based observation corrections under `selection(conditional_mtar)` (Section 5):

1. **Pure outcome-change imputation (`selmethod(imputation)`)**:
   Imputes missing within-cohort outcome changes $\Delta Y_{i,s,t}$ over the cohort population using linear outcome-change regressions controlled by `imppooling(separate|within_cohort|across_cohorts)`. Because no pair-observation logit is fitted under pure imputation, `selmodel()`, `selpooling()`, and `selcovars()` do not apply and are listed in `e(inactive_options)`.
2. **Augmented observation correction (`selmethod(dr)`)**:
   Requires `selection(conditional_mtar)` with `seluniverse(F)` (and no `selrho()`) and works with `method(ipw|reg|dr)`; `method(custom)` is not supported. It combines inverse pair-observation probability weighting (`selmodel()`, `selpooling()`, `selcovars()`) with linear outcome-change imputation (`imppooling()`) so cohort mean-change moments remain consistent if either the pair-observation model or the outcome-change imputation model is correctly specified (and doubly robust in both layers when paired with `method(dr)`).

```stata
* 5a. Pure outcome-change imputation (`selmethod(imputation)`)
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) selection(conditional_mtar) ///
    selmethod(imputation) imppooling(separate) ///
    biters($CDID_DEMO_BITERS)

* 5b. Doubly robust observation + treatment adjustment (`selmethod(dr)` with `method(dr)`)
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) selection(conditional_mtar) ///
    selmethod(dr) selmodel(pair_logit) selpooling(separate) ///
    selcovars(i.x_bin c.x_cont c.z_sel) imppooling(separate) ///
    biters($CDID_DEMO_BITERS)
```

---

## 6. Clustering and sampling weights

- **Clustering (`cluster(varname)` / `framecluster(name)`)**: Sums the unit influence functions within clusters to compute the analytical covariance `e(V)`, multiplier bootstrap standard errors `e(se)`, simultaneous confidence bands `e(ci)`, and Wald tests. **Clustering affects inference only**: it changes neither the GMM weighting matrix nor the point estimates.
- **Sampling weights (`wvar(varname)` / `[iweight=...]` and `fixweights()`)**: When panel weights vary across calendar dates, `fixweights(first_period)` uses each unit's weight at the first study date, where it must be supplied, while `fixweights(base_period)` uses the weight at each target cohort's reference date:

```stata
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    wvar(w_tv) fixweights(first_period) cluster(cluster_id) ///
    biters($CDID_DEMO_BITERS)
```

---

## 7. Postestimation Aggregations (`estat`) and Pointwise vs. Simultaneous Bands

After `cdid` (unless `noinffunc` was specified), `estat` computes four aggregations from the stored results, without using the data in memory or refitting any model (Section 7):

- `estat simple`: Overall weighted average of post-treatment $ATT(g,t)$ cells.
- `estat group`: Cohort-specific average treatment effects across post-treatment periods, plus an overall cohort-weighted ATT.
- `estat dynamic`: Event-study aggregation by relative event time $e = t - g$.
- `estat calendar`: Calendar-period aggregation across cohorts treated by date $t$.

### Inference settings, event windows and confidence bands
By default, `estat` uses the `bstrap`, `cband`, `biters`, `level` and `cluster` settings of the `cdid` estimation. You can customize event-study windows and compare simultaneous confidence bands (`cband(1)`) against pointwise intervals (`cband(0)`):

```stata
estat simple
estat group
estat calendar
estat dynamic

* 7b. Dynamic event-study selectors (`mine`, `maxe`, `balancee`, `narm`)
estat dynamic, mine(-2) maxe(2) narm
estat dynamic, balancee(1) narm

* 7c. Comparing simultaneous confidence bands (`cband(1)`) vs pointwise intervals (`cband(0)`)
estat dynamic, cband(1)
estat dynamic, cband(0)
```

`e(se)` holds the reported standard errors (from the multiplier bootstrap when `bstrap(1)`), whereas `e(V)` and `e(analytic_se)` hold the analytical covariance and standard errors.

---

## 8. Pretrend tests: the Wald test and `cdid_pretest`

`cdid` provides two pretrend tests (Section 8):

1. **Joint Wald test (`e(wald_*)`)**:
   Computed by `cdid` from all requested pretreatment cells other than reference cells, and not computed if any of them or their joint covariance is unavailable, using the analytical covariance `e(V)`. Stored in `e(wald_statistic)`, `e(wald_df)`, `e(wald_p)`, and `e(wald_status)`.
2. **Conditional Cramér-von Mises test (`cdid_pretest`)**:
   Tests conditional parallel trends over the covariate distribution with a Cramér-von Mises (CvM) statistic and a multiplier bootstrap. `cdid_pretest` returns its results in `r()`, uses `method(ipw)` by default (unlike `cdid`, whose default is `method(dr)`), and leaves the current `e()` results unchanged:

```stata
* 8a. The pretrend Wald test stored in `e()`
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) biters($CDID_DEMO_BITERS)
display as text "Wald test status: " as result "`e(wald_status)'" ///
    as text "  chi2(" as result e(wald_df) as text ") = " as result %9.4f e(wald_statistic) ///
    as text "  p = " as result %6.4f e(wald_p)

* 8b. Conditional pretrend CvM test (`cdid_pretest`)
*     The default treatment adjustment is `ipw`; `e()` from `cdid` is unchanged.
cdid_pretest y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    biters($CDID_DEMO_BITERS)
return list
```

---

## 9. Diagnostics, graphs and datasets

`cdid_diagnostics` summarizes the study population, the comparisons available, the model and the inference behind the current results; with `detail`, it also lists the unavailable cells. Quantities that were not stored with the results are reported as not retained, never as zero. The test of `cdid_pretest` is reported by that command, not by `cdid_diagnostics`.

`cdid_plot` graphs the current results. After `estat dynamic`, it graphs the event study; `cdid_plot, cells` graphs the group-time estimates instead. Graph options such as `title()` can be added. `saving()` writes a `.gph` file, `export()` writes the graph in another format with `graph export`, and `plotdata()` saves the plotted values as a dataset. Existing files are replaced only with `replace`.

`cdid_export` saves the stored results as a Stata dataset, with `type(cells|aggregate|units|influence|summary)`; `cells` is the default. It does not refit or draw random numbers. After `estat`, `type(cells)` still saves every group-time cell, and `type(aggregate)` saves the aggregation. Reported and analytical standard errors are saved in separate variables. After an aggregation in ordinary mode, the pointwise intervals of the group-time cells can be missing; the reported intervals are always saved.

```stata
* 9a. Diagnostics
cdid_diagnostics
cdid_diagnostics, detail

* 9b. Export cell and summary tables to temporary files and inspect them in frames
* cdid_export accepts temporary file names as they are; Stata erases the files
* when the do-file ends. Files are kept on disk only in 9c and Section 10.
tempfile demo_cells_dta demo_summary_dta
cdid_export using "`demo_cells_dta'", type(cells)
capture frame drop demo_cells_tbl
frame create demo_cells_tbl
frame demo_cells_tbl: use "`demo_cells_dta'"
frame demo_cells_tbl: list in 1/10, abbrev(16)

cdid_export using "`demo_summary_dta'", type(summary)
capture frame drop demo_summary_tbl
frame create demo_summary_tbl
frame demo_summary_tbl: use "`demo_summary_dta'"
frame demo_summary_tbl: list, abbrev(20)

* 9c. Graphs and datasets kept on disk (requires CDID_DEMO_OUTDIR)
if `"$CDID_DEMO_OUTDIR"' != "" {
    local repopt = cond($CDID_DEMO_REPLACE == 1, "replace", "")
    cdid_plot, title("Group-Time ATT Estimates") ///
        saving("$CDID_DEMO_OUTDIR/cdid_group_time.gph") `repopt'
    estat dynamic
    cdid_plot, title("Dynamic Event-Study ATT") ///
        saving("$CDID_DEMO_OUTDIR/cdid_dynamic.gph") `repopt'
    cdid_export using "$CDID_DEMO_OUTDIR/cdid_dynamic.dta", type(aggregate) `repopt'
}
```

Section 9 writes temporary datasets and opens them in frames. Stata erases the files at the end of the do-file. Set `CDID_DEMO_OUTDIR` to an existing folder to keep the `.gph` and `.dta` files. For CSV or other formats, use Stata's export commands on the saved dataset.

---

## 10. Saving and loading results (`saving()`, `cdid_load`, `cdid_pretest_load`)

`cdid`, `estat` and `cdid_pretest` accept `saving(filename)`, which writes the complete results to a file. `saving()` never overwrites an existing file.

You can restore saved results in a later Stata session, even when memory is empty or holds an unrelated dataset, using `cdid_load` and `cdid_pretest_load` (Section 10):

```stata
cdid_demo_init

if `"$CDID_DEMO_OUTDIR"' != "" {
    local cells_file "$CDID_DEMO_OUTDIR/demo_cells.cdidstate"
    local dyn_file   "$CDID_DEMO_OUTDIR/demo_dynamic.cdidstate"
    local pt_file    "$CDID_DEMO_OUTDIR/demo_pretest.cdidpt"
}
else {
    tempfile cells_file dyn_file pt_file
}

* 10a. Save the group-time results, the event study, and the pretrend test
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    cluster(cluster_id) biters($CDID_DEMO_BITERS) saving("`cells_file'")

estat dynamic, saving("`dyn_file'")

cdid_pretest y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    biters($CDID_DEMO_BITERS) saving("`pt_file'")

* 10b. Replace the data with a small unrelated dataset and load the results
preserve
clear
set obs 5
generate byte unrelated_x = _n

* Load the group-time results and use `estat` and `cdid_diagnostics`
cdid_load using "`cells_file'"
count if e(sample)
assert r(N) == 0
display as text "e(sample_scope): " as result "`e(sample_scope)'"
display as text "e(cmdline): " as result `"`e(cmdline)'"'
cdid_diagnostics
estat group

* Load the event study
cdid_load using "`dyn_file'"

* Load with `posting(state)` (complete results, without e(b) and e(V))
cdid_load using "`cells_file'", posting(state)
display as text "Coefficients posted with posting(state): " as result e(N_posted)

* Load the pretrend test results (e() is unchanged)
cdid_pretest_load using "`pt_file'"
return list
restore


* ==============================================================================
```

### What loading does
- **No re-estimation**: `cdid_load` restores in `e()` exactly the saved estimates, standard errors `e(se)`, covariance `e(V)`, intervals and bands, and group-time and unit tables; `cdid_pretest_load` restores the test results in `r()`. They draw no random numbers; `c(rngstate)` and `c(sortrngstate)` are unchanged.
- **Command line.** `e(cmdline)` holds the original `cdid` command of the saved fit or aggregation.
- **Size limits**: `posting(matrix)`, the default, posts the results as matrices and fails with `r(908)` when they exceed `c(max_matdim)`; `posting(state)` loads the complete results without posting `e(b)` and `e(V)` (`cdid_pretest_load` then requires `store(name)`).
- **Errors**: damaged files and files of the wrong type (for example, a file saved by `cdid` given to `cdid_pretest_load`) cause error `r(498)`; the data and the current results are then unchanged.

---

## 11. Common errors

Section 11 shows errors caused by incomplete specifications or insufficient data. Each command stops with an error message and leaves the data and the current results unchanged:

```stata
* 11a. `ptassumption(all)` without `ptstart(#)` (rc = 459)
capture noisily cdid y x_bin, ivar(id) time(period) gvar(first_treat) ptassumption(all)
assert _rc == 459

* 11b. `selmethod(dr)` with `selection(none)` (rc = 459)
capture noisily cdid y x_bin, ivar(id) time(period) gvar(first_treat) ///
    selection(none) selmethod(dr)
assert _rc == 459

* 11c. Insufficient support: restricting to a single calendar period leaves no (s, t) pairs
capture noisily cdid y x_bin if period == 1, ivar(id) time(period) gvar(first_treat)
assert _rc != 0

* 11d. Wrong file type: a file saved by `cdid` given to `cdid_pretest_load` (rc = 498)
tempfile tmp_cell_state
cdid y x_bin, ivar(id) time(period) gvar(first_treat) noinffunc saving("`tmp_cell_state'")
capture noisily cdid_pretest_load using "`tmp_cell_state'"
assert _rc == 498
```

---

## References

- Bellégo, C., Benatia, D., & Dortet-Bernadet, V. (2025). *The Chained Difference-in-Differences*. Journal of Econometrics, 248, 105783. <https://doi.org/10.1016/j.jeconom.2024.105783>
- Callaway, B., & Sant'Anna, P. H. C. (2021). *Difference-in-Differences with multiple time periods*. Journal of Econometrics, 225(2), 200-230. <https://doi.org/10.1016/j.jeconom.2020.12.001>

[Guide home](README.md)
