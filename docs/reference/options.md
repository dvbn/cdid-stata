# Options reference

This page lists every option of every command with its default, its accepted values and the estimators that accept it. The help files (`help cdid`, `help cdid_dvar`, `help cdid postestimation` and the help file of each command) give the same information in Stata's format, and [Moving from R cdid](../from-r-cdid.md) gives the R argument of each option.

In the column **Estimators**, *chained* is chained cohort estimation with `gvar()` (the default), *cross-pair* is chained estimation with `crosspair`, *preset* is `preset(cohort_time_fe)`, and *ordinary* is `mode(ordinary)` or `rcs`. Options of `dvar()` estimation have their own section.

## Signatures

```text
cdid depvar [indepvars] [if] [in] [iweight], time(varname) gvar(varname) [ivar(varname) options]
cdid depvar [indepvars] [if] [in], ivar(varname) time(varname) dvar(varname) [dvar_options]
cdid
estat {simple|group|dynamic|calendar} [, estat_options]
cdid_pretest depvar [indepvars] [if] [in] [iweight], time(varname) gvar(varname) [ivar(varname) options]
cdid_chain [, cohort(numlist) time(numlist) pre weights tol(#) saving(filename[, replace]) table(name)]
cdid_nuisance [selection|outcome|treatment|imputation] [, block(id)]
cdid_diagnostics [, detail]
cdid_plot [, group(numlist) cells saving(filename[, replace]) export(filename[, replace options]) plotdata(filename[, replace]) replace twoway_options]
cdid_export using filename [, type(cells|aggregate|units|influence|summary) replace]
cdid_load using filename [, posting(matrix|state) workspace(#)]
cdid_pretest_load using filename [, store(name) posting(matrix|state) workspace(#)]
```

Minimal abbreviations: `nocons(tant)`, `control(group)`, `anticip(ation)`, `base(period)`, `l(evel)`, `cross(pair)`, `missingref(erence)`, `selnocons(tant)`, `selcons(tant)`.

## Data, timing and comparison group

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `depvar` | required | all | Numeric outcome; missing values are unobserved outcomes. |
| `time(varname)` | required | all | Numeric study date. |
| `gvar(varname)` | required | all | First treatment date, constant within unit, nonnegative and nonmissing; 0 for never-treated units. |
| `ivar(varname)` | required in chained mode | all | Unit identifier, numeric or string; numeric in ordinary mode; may be omitted with `rcs`. |
| `controlgroup(nevertreated\|notyettreated)` | `nevertreated` | all | Never-treated units, or never-treated units together with units not yet treated at the dates compared. |
| `anticipation(#)` | `0` | all | Periods of anticipation; the reference date of cohort $g$ is the latest study date $d$ with $d + \# < g$. The preset requires 0. |
| `targetcohorts(numlist)` | all treated cohorts | chained, cross-pair, preset | Target cohorts; not with `targets()`. |
| `targets(matname)` | all treated cohorts | chained, cross-pair, preset | Target cohorts as a one-column matrix; not with `targetcohorts()`. |

By default, every treated cohort whose treatment date is on or before the last study date is a target; the effects of a cohort g with no study date d such that d + anticipation < g are reported as unavailable.

## Treatment adjustment

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `indepvars`, `covars(varlist)` | none | all; `covars()` not in ordinary mode | Covariates of the treatment, outcome and imputation models; one or the other. Factor variables allowed. |
| `noconstant` | off | chained, cross-pair, preset | Omit the intercept from the treatment and outcome designs. |
| `method(dr\|ipw\|reg\|custom)` | `dr` | all; preset: `reg`, `ipw`, `dr` | Doubly robust, inverse probability weighting, outcome regression, or a Mata function. |
| `outcomepooling(separate\|pooled)` | `separate` | chained, cross-pair; preset: `separate` only | One regression of the controls' changes per comparison, or slopes shared across comparisons with local intercepts; applies to `reg` and `dr`. |
| `callback(name)`, `payload(name)` | none | chained, ordinary | Mata function and object for `method(custom)`; see [Custom estimators](../custom-estimators.md). |

## Covariates that vary over time

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `missingreference(exclude\|error)` | `exclude` | chained, cross-pair; ignored in ordinary mode; `exclude` refused with the preset | When a covariate that varies over time is missing at the date where a cohort's comparisons read it: exclude the unit from that cohort's comparisons with a warning, or stop with error 459 and list the missing values. Identical results when nothing is missing. Weights are not covered: with `fixweights()`, every unit used needs a positive weight at its anchor date. |
| `fixweights(first_period\|base_period)` | none | chained, cross-pair; ordinary also `varying` | Weight at the first study date or at each target cohort's reference date when weights vary within units. Without it, weights must be constant within units. |

Covariates that vary over time are read at the reference date of each target cohort for its post-treatment comparisons and at the earlier date of each pretrend comparison; with `crosspair`, every comparison of cohort $g$ reads them at $g$'s reference date. See [Prepare your panel](../help/prepare-data.md#reference-dates-for-covariates-that-change-within-a-unit).

## Handling missing outcomes

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `selection(none\|conditional_mtar)` | `none` | chained, cross-pair, preset | Assumption about missing outcomes; `corollary1` is a synonym of `conditional_mtar`. |
| `selmethod(ipw\|imputation\|dr)` | `ipw` | chained, cross-pair, preset | Observation weighting, change imputation or the augmented correction. `method(custom)` supports `ipw` only. |
| `selmodel(pair_logit\|independent_period_logit)` | `pair_logit` | chained, cross-pair, preset | Logit for observing both dates of a pair, or product of logits for each date. |
| `selpooling(separate\|within_cohort\|across_cohorts)` | `separate` | chained, cross-pair, preset | Slopes of the observation logits. |
| `imppooling(separate\|within_cohort\|across_cohorts)` | `separate` | chained, cross-pair, preset | Slopes of the imputation regressions. |
| `selcovars(varlist)` | the outcome covariates | chained, cross-pair, preset | Observation-model covariates; factor variables allowed. |
| `selnoconstant` | off | chained, cross-pair, preset | No intercept in the observation logits. |
| `selconstant` | off | chained, cross-pair, preset | Intercept-only observation logits; not with `selcovars()`, `frameselcovars()`, `selnoconstant` or `selzero`. |
| `selzero` | off | chained, cross-pair, preset | No observation design. |
| `seluniverse(F\|B\|A)` | `F` | chained, cross-pair, preset | Units that fit the observation model: study population, units with at least one, or at least two, observed outcomes. `B` and `A` require `selmethod(ipw)`, `selmodel(pair_logit)` and `seljustification()`. |
| `selrho(varname)` | none | chained, cross-pair, preset | Known density ratio with `B` or `A`. |
| `seljustification(string)` | none | chained, cross-pair, preset | Stated justification; required with `B` or `A`. |

See [Handling missing outcomes](../methods/missing-outcomes.md) and [Sampling information](../methods/sampling-information.md).

## Parallel trends and reporting

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `ptassumption(post\|ctrl\|all)` | `post` | chained, cross-pair, preset | PT from each reference date onward; `ctrl` adds common trends among eligible control cohorts from `ptstart()`, which is optional and defaults to the first study date; `all` adds common trends of target and controls from `ptstart()`. |
| `ptstart(#)` | none | chained, cross-pair, preset | A study date; required with `all`, optional with `ctrl`, not allowed with `post`. |
| `baseperiod(varying\|universal)` | `varying` | all; preset: `universal` | Reporting of pretreatment cells: relative to the preceding date, or to the reference date. |

## Chaining and designs

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `crosspair` | off | chained | Separate within-unit changes of target and control cohorts linked by GMM; not in ordinary mode, with `method(custom)` or with the preset. See [Cross-pair chaining](../methods/cross-pair.md). |
| `preset(cohort_time_fe)` | none | chained | Linear cohort-by-date fixed effects; requires chained mode without `crosspair`, unit weights, `anticipation(0)`, `baseperiod(universal)`, `ptassumption(post)` with never-treated controls or `ptassumption(all) ptstart(first date)` with not-yet-treated controls, and positive consecutive integer dates. See [Cohort-by-date fixed effects](../methods/cohort-time-fe.md). |
| `keepnuisance` | off | chained, cross-pair, preset | Keep reports of the fitted models for `cdid_nuisance`. |

## GMM weighting and covariance regularization

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `gmmweighting(optimal\|two_step\|within\|identity)` | `optimal` | chained, cross-pair, preset | `two_step` is a synonym of `optimal`. `within` is not available with `method(custom)`. With the preset, `method(reg)` allows only `within`, and `method(ipw)` or `method(dr)` allow `optimal`, `two_step` or `identity`, not `within`. |
| `gmmcovariance(tikhonov)` | `tikhonov` | chained, cross-pair, preset | Only value. |
| `condition(#)` | `1e-10` | chained, cross-pair, preset | Relative eigenvalue threshold in (0, 1) that triggers regularization. |
| `ridge(#)` | `1e-4` | chained, cross-pair, preset | Relative Tikhonov multiplier, > 0. |

See [GMM weighting](../methods/gmm.md).

## Study population and weights

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `frame(name)` | none | chained, cross-pair, preset | Frame with one row per study unit. |
| `frameivar(name)` | the `ivar()` name | chained, cross-pair, preset | Unit identifier in the frame. |
| `framegvar(name)`, `framecovars(varlist)`, `frameselcovars(varlist)`, `framewvar(name)`, `framecluster(name)` | not read | chained, cross-pair, preset | Frame columns; read only when named; require `frame()`. |
| `iweight`, `wvar(varname)` | none | chained, cross-pair, ordinary | Positive sampling weights; one or the other; refused with the preset. |

## Inference

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `bstrap(#)` | `1` | all | 1: multiplier bootstrap; 0: analytical standard errors. |
| `cband(#)` | `1` | all | 1: simultaneous bands; 0: pointwise intervals. |
| `biters(#)` | `1000` | all | Bootstrap draws, at least 1. |
| `level(#)` | `95` | all | Confidence level. |
| `cluster(varname)` | none | all | One cluster variable, constant within units. |
| `noinffunc` | off | all | Point estimates only; `estat` is then refused. |

## Mode

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `mode(chained\|ordinary)` | `chained` with panel data | all | Estimator family. Ordinary mode uses `csdid` 2.0.2. |
| `rcs` | off | ordinary | Repeated cross sections; implies ordinary mode; refused with `mode(chained)`. |
| `panel(true\|1)` | | chained, cross-pair, preset | Accepted; redundant. Refused in ordinary mode. |
| `allowunbalanced`, `allowunbalanced(true\|false)` | | all | In ordinary mode, keep units not observed at every date (repeated cross-section estimators). In chained mode only `true` is accepted, and it changes nothing. |

Ordinary estimation takes covariates from the varlist and refuses `covars()`, `noconstant`, `panel()`, `frame()` and its associated options. `e(inactive_options)` lists every accepted option that has no effect in ordinary mode: `ptassumption()`, `ptstart()`, `targets()`, `targetcohorts()`, `gmmweighting()`, `gmmcovariance()`, `workcap()`, `condition()`, `ridge()`, `selection()`, the `sel*()` options, `outcomepooling()`, `imppooling()`, `missingreference()` and `keepnuisance`. The note names those supplied in the command. `ivar()` is required for panel data and may be omitted with `rcs`.

## Computation and reporting

| Option | Default | Estimators | Description |
| --- | --- | --- | --- |
| `workspace(#)` | `1073741824` | all | Memory limit for the main numerical arrays, in bytes. |
| `chunk(#)` | automatic; `1024` in ordinary mode | chained, cross-pair, preset, ordinary | Rows per processing block. |
| `workcap(#)` | unlimited | chained, cross-pair, preset | Limit on counted operations. |
| `posting(matrix\|state)` | `matrix` | chained, cross-pair, preset, ordinary | Post `e(b)` and `e(V)`, or keep the results without posting them. |
| `saving(filename)` | none | chained, cross-pair, preset, ordinary | Save the complete results for `cdid_load`; never overwrites. |

## Options with dvar()

| Option | Default | Description |
| --- | --- | --- |
| `dvar(varname)` | required | Treatment status, 0 or 1, no missing values. |
| `ivar()`, `time()` | required | Unit identifier and date. |
| `indepvars`, `covars(varlist)` | none | Covariates, read at $t-1$; factor variables allowed. |
| `method(dr\|ipw\|reg)` | `dr` | Treatment adjustment. |
| `switchers(both\|entry\|exit)` | `both` | Directions of the requested cells. |
| `switchtimes(numlist)` | dates with switchers | Dates of the requested cells. |
| `outcomepooling(separate\|pooled\|`*blocks*`)` | `separate` | One stayer regression per cell, common slopes, or blocks such as `entry:2003 entry:2004 \| exit:2003 exit:2004`. |
| `gmmweighting(optimal\|two_step\|identity)` | `optimal` | With chained comparisons under missing outcomes. |
| `gmmcovariance(tikhonov)`, `condition(#)`, `ridge(#)` | `tikhonov`, `1e-10`, `1e-4` | As for `gvar()`. |
| `selection(none\|conditional_mtar)` | `none` | Assumption about missing outcomes. |
| `selmethod(ipw\|imputation\|dr)` | `ipw` | `imputation` and `dr` require `conditional_mtar`. |
| `selmodel(pair_logit\|independent_period_logit)` | `pair_logit` | Observation model. |
| `selpooling()`, `imppooling()` | `separate` | `separate`, `within_group` or `pooled`. |
| `impscale(per_period\|per_change)` | `per_period` | Scaling of the imputation covariates by the length of a pair. |
| `selcovars(varlist)`, `selconstant` | outcome covariates | Not both. |
| `seluniverse(F\|B\|A)` | `F` | `B` and `A` require `conditional_mtar`, `ipw`, `pair_logit` and `seljustification()`. |
| `selrho(varname)`, `seljustification(string)` | none | With `B` or `A`. |
| `cluster(varname)` | none | At least two clusters. |
| `noinffunc`, `bstrap()`, `cband()`, `biters()`, `level()`, `workspace()` | as for `gvar()` | |
| `keepnuisance` | off | Keep reports for `cdid_nuisance`. |
| `mode(chained)` | | Accepted; the only mode. |
| `missingreference(error)` | | Accepted; `exclude` is refused. |

Weights, `gvar()` and every other option, including `saving()`, `crosspair`, `controlgroup()`, `ptassumption()`, `anticipation()`, `method(custom)` and `gmmweighting(within)`, are refused with error 198.

## estat options

| Option | Default | Description |
| --- | --- | --- |
| `mine(#)`, `maxe(#)` | none | `estat dynamic`: bounds of the event time t - g, in date units. `maxe()` also limits `estat simple` and `estat group` to post-treatment cells whose horizon, counted in positions of the sorted dates, is at most #. |
| `balancee(#)` | none | `estat dynamic`: keep the cohorts observed for at least # periods after treatment. |
| `narm` | off | Aggregate the available cells when some are unavailable. |
| `bstrap(#)`, `cband(#)`, `biters(#)`, `level(#)` | those of the fit | Inference of the aggregation. |
| `cluster(varname)` | that of the fit | After an ordinary fit, `cluster(none)` or the `ivar()` variable requests unit-level inference. A request for an unavailable cluster variable warns and uses unit-level inference. After a chained fit, a request to change clustering stops with an error. |
| `workspace(#)` | `1073741824` | Memory limit in bytes. |
| `chunk(#)` | `1024` | Rows per processing block. |
| `posting(matrix\|state)` | `matrix` | As for `cdid`. |
| `saving(filename)` | none | Save the aggregation for `cdid_load`. |

`mine()` and `balancee()` have no effect on `estat simple` and `estat group`; `estat calendar` ignores all three, with a note. `estat` is not available after `dvar()` estimation or after `noinffunc`.

## cdid_pretest options

`cdid_pretest` accepts the model, mode, missing-outcome, study-population and weight options of `cdid` (`ivar()`, `covars()`, `noconstant`, `mode()`, `rcs`, `controlgroup()`, `anticipation()`, `baseperiod()`, `outcomepooling()`, `callback()`, `payload()`, `wvar()`, `fixweights()`, `frame()`, `frameivar()`, `framegvar()`, `framecovars()`, `framewvar()`, `framecluster()`, `frameselcovars()`, the `sel*()` options, `imppooling()`, `cluster()`, `workspace()`, `chunk()`, `workcap()`), with these differences:

| Option | Default | Description |
| --- | --- | --- |
| `method(ipw\|reg\|dr\|custom)` | `ipw` | Treatment adjustment. |
| `panel(0\|1)` | | 0 requests repeated cross sections. |
| `allowunbalanced(0\|1)` | | Chained mode refuses 0; ordinary mode keeps units observed at every date unless 1. |
| `bstrap(0\|1)`, `cband(0\|1)` | | Accepted; the test always uses the bootstrap. |
| `biters(#)`, `level(#)` | `1000`, `95` | Bootstrap draws and level of the critical value. |
| `pl(0\|1)`, `cores(#)` | `0`, `1` | Only `pl(0) cores(1)` is accepted; the test runs serially. |
| `signs(matname)` | none | Bootstrap multipliers supplied as a matrix of -1 and +1. |
| `posting(matrix\|state)`, `store(name)` | `matrix` | `state` requires `store()`. |
| `saving(filename)` | none | Save for `cdid_pretest_load`. |

The ordinary pretest accepts `anticipation()`, `baseperiod()`, `fixweights()`, `frame()` and its associated options, the missing-outcome options, `imppooling()`, `outcomepooling()`, `workspace()`, `chunk()` and `workcap()` as inactive. `r(inactive_options)` lists them whether specified or omitted. The note names those supplied in the command.

## cdid_chain options

| Option | Default | Description |
| --- | --- | --- |
| `cohort(numlist)`, `time(numlist)` | all | Cells reported. |
| `pre` | off | Include pretreatment and reference cells. |
| `weights` | off | List the GMM weight of each comparison. |
| `tol(#)` | `1e-10` | Weights below # in absolute value are treated as zero. |
| `saving(filename[, replace])` | none | Save a table as a dataset. |
| `table(summary\|witness\|weights)` | `summary` | Table to save; requires `saving()`. |

`cdid_chain` is refused after `dvar()` and in ordinary mode (error 198) and after the preset (error 498).

## Other commands

| Command and option | Default | Description |
| --- | --- | --- |
| `cdid_nuisance` *component*, `block(id)` | all blocks | Component and block to display. |
| `cdid_diagnostics, detail` | off | List each unavailable cell and bin. |
| `cdid_plot, group(numlist)` | all | Cohorts to plot. |
| `cdid_plot, cells` | off | Plot the cells after `estat`. |
| `cdid_plot, saving()`, `export()`, `plotdata()`, `replace` | none | Save the graph, export it with `graph export`, or save the plotted values. |
| `cdid_export, type()` | `cells` | `cells`, `aggregate`, `units`, `influence` or `summary`; after `dvar()`, `cells` or `aggregate`. |
| `cdid_export, replace` | off | Overwrite the file. |
| `cdid_load, posting()`, `workspace()` | `matrix`, `1073741824` | As for `cdid`. |
| `cdid_pretest_load, store()`, `posting()`, `workspace()` | none, `matrix`, `1073741824` | `posting(state)` requires `store()`. |

## Refusals

Unless noted, refusals stop with error 198:

- `crosspair` requires chained mode; `rcs` with `mode(chained)` is refused.
- `iweight` together with `wvar()`, and `covars()` together with a covariate varlist, are refused.
- `method(custom)` is refused with `selmethod(imputation)` or `selmethod(dr)` (error 459), with `crosspair` and with `gmmweighting(within)`.
- Ordinary `cdid` refuses `frame()` and its associated options, `noconstant`, `covars()` and `panel()`.
- The preset requires the settings listed under [Chaining and designs](#chaining-and-designs): unmet option settings are refused with error 198, and dates that are not consecutive positive integers with error 459. It also refuses `missingreference(exclude)`.

The following are also refused: `targets()` with `targetcohorts()`; `selconstant` with `selcovars()`, `frameselcovars()`, `selnoconstant` or `selzero`; `ptstart()` with `ptassumption(post)`; and `ptassumption(all)` without `ptstart()` (error 459).

---

**Reading path, step 18 of 22** · Previous: [Clustered inference, confidence bands, and pretrend tests](../examples/inference.md) · Next: [Results reference](results.md)

[Guide home](../README.md) · [Package overview](../../README.md)
