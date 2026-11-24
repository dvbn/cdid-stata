# Estimate familiar DiD models with cdid

Several widely used difference-in-differences estimators correspond to particular settings of `cdid`. This page states, for each, the settings, the sample and the estimand under which the point estimates coincide, and how the correspondence was checked for this guide. A numerical match is reported only where both sides were run in Stata on the same simulated data; the other correspondences are stated as specifications. Choose the comparison groups and treatment history first, then the treatment adjustment, then whether to correct for selective observation.

## Matched specifications

For a numerical match, use the same sample, outcome and covariates, the same untreated reference, and the same weights on the reported effects. Match the variance convention separately. Here `g` is the first treatment date of an absorbing treatment, and `d` is the current status of a binary treatment that can switch on and off.

| To reproduce | Settings in `cdid` | Where the match holds | Checked for this guide |
| --- | --- | --- | --- |
| Callaway and Sant'Anna, Stata `csdid` | `mode(ordinary)`, or a chained fit in the [matched cases](from-did.md#when-chained-and-ordinary-estimates-should-agree) | Same two-date samples, comparison group, base period and adjustment. | Run in Stata with `csdid` 2.0.2: [below](#callaway-and-santanna-with-csdid). |
| Wooldridge's linear cohort-by-date regression (Stata `jwdid`, `xthdidregress twfe`; R `etwfe`) | `preset(cohort_time_fe) method(reg)` | Unit and date fixed effects with one effect per treated cohort and date; never-treated or not-yet-treated controls; time-invariant covariates. | Run in Stata against the same regression fitted with `areg`; `jwdid` was not run. [Below](#the-linear-cohort-by-date-regression). |
| Sun and Abraham, Stata `eventstudyinteract` | Same preset with `controlgroup(nevertreated)` | The saturated cohort-by-event-time regression with never-treated controls and the same reference period. | Same regression as the row above; `eventstudyinteract` was not run. |
| Borusyak, Jaravel and Spiess, Stata `did_imputation`; Gardner, Stata `did2s` | Same preset with `controlgroup(notyettreated) ptassumption(all)` and `ptstart()` at the first date | Complete, unweighted panel; untreated outcomes modeled with unit and date fixed effects. | Run in Stata against the imputation computed with `regress`; `did_imputation` and `did2s` were not run. [Below](#two-stage-and-imputation-event-studies). |
| Two-date clean-control comparison (local-projection DiD) | `method(reg)` on the two dates, or the preset | Complete panel, one adoption cohort, never-treated controls, no covariates. | Run in Stata against `regress` of the long difference. [Below](#the-two-date-comparison). |
| Harmon's stepwise DiD | One two-date `cdid` fit per consecutive step, summed | Each step uses the cohorts still untreated at its later date. | Run in Stata against the step differences of means. [Below](#stepwise-did). |
| Balla-Elliott and Norwich's Cohort-Chained DiD | `controlgroup(notyettreated) baseperiod(universal) crosspair ptassumption(ctrl) gmmweighting(within)` | Stacked target-cohort regressions on the same sample and weights, without covariates. | Tested in Stata against `areg` and `regress`: [below](#cohort-chained-did). |
| de Chaisemartin et al., static binary switchers, Stata `did_multiplegt_stat` | `dvar(d) method(reg)` | Complete panel, no covariates, switcher-weighted average of entries and exits. | Not run in Stata for this guide; an R result is [cited below](#binary-static-switcher-summaries). |

## Callaway and Sant'Anna with csdid

The examples on this page use a complete panel of 72 units over six dates: a third never treated, a third first treated at date 3 and a third at date 5. Each unit has a fixed effect, the untreated path is $\sin(t)$, and the effect grows with time since adoption. The built-in chained comparisons specify `gmmweighting(optimal)` to match the recorded output.

```stata
clear
set seed 14603
set obs 72
generate id = _n
generate g = cond(mod(id - 1, 3) == 0, 0, cond(mod(id - 1, 3) == 1, 3, 5))
generate a = rnormal()
expand 6
bysort id: generate t = _n
generate y = a + sin(t) + rnormal() + cond(g > 0 & t >= g, (t - g + 1) / 2, 0)
generate byte treated = g > 0 & t >= g
```

Ordinary mode computes its estimates with `csdid` 2.0.2. Run directly, `csdid` with not-yet-treated controls gives the same post-treatment cells; the pretreatment cells differ only because `csdid` 2.0.2 uses a universal base period by default and `cdid` a varying one. A chained fit with IPW, no covariates and never-treated controls also gives the same post-treatment cells on this complete panel, because the cohort means then telescope.

```stata
csdid y, ivar(id) time(t) gvar(g) notyet method(ipw)
matrix b_csdid = e(b)
cdid y, ivar(id) time(t) gvar(g) mode(ordinary) controlgroup(notyettreated) method(ipw) bstrap(0) cband(0)
matrix b_ordinary = e(b)
csdid y, ivar(id) time(t) gvar(g) nevertreated method(ipw)
matrix b_never = e(b)
cdid y, ivar(id) time(t) gvar(g) method(ipw) gmmweighting(optimal) bstrap(0) cband(0)
matrix b_chained = e(b)
mata: post = (2, 3, 4, 5, 9, 10)
mata: max(abs(st_matrix("b_csdid")[post] - st_matrix("b_ordinary")[post]))
mata: max(abs(st_matrix("b_never")[post] - st_matrix("b_chained")[post]))
```

The two lines of output are the largest absolute differences over the six post-treatment cells, for ordinary mode against `csdid` and for the chained fit against `csdid` with never-treated controls:

```text
  4.44089e-16
```

```text
  1.77636e-15
```

Columns 2 to 5, 9 and 10 of `e(b)` are the cells (3,3) to (3,6), (5,5) and (5,6) in both commands. Standard errors were not compared. In an incomplete panel, `csdid`'s `bal(pair)` keeps the units with both outcomes in each two-date comparison, whereas chained estimation can also combine other supported pairs into an effect, so the estimates generally differ.

## The linear cohort-by-date regression

The preset with never-treated controls fits the regression of the outcome on unit fixed effects, date fixed effects and one indicator for each treated cohort and date other than the date before adoption. The same regression can be fitted with `areg`:

```stata
cdid y, ivar(id) time(t) gvar(g) preset(cohort_time_fe) method(reg) baseperiod(universal) controlgroup(nevertreated) bstrap(0) cband(0)
matrix b_preset = e(b)
foreach c in 3 5 {
    forvalues s = 1/6 {
        if `s' != `c' - 1 generate byte D`c'_`s' = g == `c' & t == `s'
    }
}
areg y i.t D3_* D5_*, absorb(id)
mata: max(abs(st_matrix("b_preset") - st_matrix("e(b)")[7..16]))
```

```text
  1.09912e-14
```

The ten coefficients, including the four pretreatment deviations, match up to the printed difference. This is the specification of Wooldridge's extended two-way fixed-effects regression (Stata `jwdid`, `xthdidregress twfe`) and, with never-treated controls, the saturated cohort-by-event-time regression of Sun and Abraham (`eventstudyinteract`); those commands were not run for this guide. The R documentation of `cdid` reports matches of the R preset with R `etwfe` and Stata `jwdid`. The preset's standard errors use the unit score sandwich without a degrees-of-freedom adjustment, so they differ from the `areg` standard errors shown; they were not compared.

The match covers unit fixed effects, never-treated or not-yet-treated controls and time-invariant covariates, which let the untreated calendar path and the cohort-by-date effects vary with the covariates; each effect is then evaluated at the covariate mean of its cohort in the study population. The preset requires a never-treated cohort, unit weights, no anticipation and positive consecutive integer dates. When observation dates vary within a cohort, compare with a regression that has unit fixed effects. Specifications outside the preset, such as time-varying regressors, sampling weights or a common effect per event time (`jwdid`'s `hettype(event)`), have no counterpart in it. See [Cohort-by-date fixed effects](methods/cohort-time-fe.md).

## How event-study averages are weighted

An event-study point at horizon 1 averages effects measured one period after adoption. If the 2020 cohort has 60 units and an estimated effect of 1 in 2021, and the 2022 cohort has 40 units and an effect of 3 in 2023, the cohort-size-weighted point is $0.6 \times 1 + 0.4 \times 3 = 1.8$. If only 20 units of the later cohort have an observed outcome in 2023, an average weighted by observed outcomes is $(60 \times 1 + 20 \times 3) / 80 = 1.5$, with the same two cohort effects.

To compare event studies, first compare the cohort-by-date effects, then average the same cohorts at each horizon with the same weights. `estat dynamic` weights each cohort by its number of units in the study population, and averages over available cells only with `narm`. Regression postestimation commands often average over observed treated rows. On a complete panel the two weightings coincide when the samples and supported cohorts are the same.

## Two-stage and imputation event studies

On a complete panel, fitting unit and date effects on the untreated observations and averaging the residuals of the treated observations by cohort and date gives the cohort-by-date effects of the imputation estimator of Borusyak, Jaravel and Spiess and of Gardner's two-stage estimator with one indicator per treated cell. The preset with not-yet-treated controls and PT-all from the first date fits the same untreated model:

```stata
cdid y, ivar(id) time(t) gvar(g) preset(cohort_time_fe) method(reg) baseperiod(universal) controlgroup(notyettreated) ptassumption(all) ptstart(1) bstrap(0) cband(0)
matrix b_preset = e(b)
estat dynamic
matrix b_dynamic = e(b)
quietly regress y i.id i.t if !treated
predict double y0
generate double tau = y - y0 if treated
preserve
collapse (mean) tau if treated, by(g t)
mkmat tau, matrix(b_imputation)
mata: max(abs(st_matrix("b_preset") - st_matrix("b_imputation")'))
restore
preserve
generate e = t - g
collapse (mean) tau if treated, by(e)
mkmat tau, matrix(b_horizon)
mata: max(abs(st_matrix("b_dynamic") - st_matrix("b_horizon")'))
restore
```

```text
  4.44089e-16
```

```text
  3.33067e-16
```

The first difference compares the six post-treatment cohort-by-date effects, the second the four event-study points; on this complete panel with equal cohort sizes, averaging the treated rows by horizon gives the same cohort weights as `estat dynamic`. `did_imputation` and `did2s` themselves were not run for this guide. The R documentation of `cdid` reports the corresponding match with R `didimputation` and `did2s`, with and without a time-invariant covariate interacted with the dates. With incomplete observations, joint fixed-effects regression and untreated-only imputation can use different information to fit the untreated path, and the weights and inference can differ.

## The two-date comparison

Take a complete panel with one adoption cohort and never-treated controls, keep a baseline and one follow-up date, and code them 1 and 2 with adoption at 2. With no covariates, the effect equals the coefficient on the treated indicator in a regression of the long outcome change, the clean-control comparison used by local-projection DiD at that horizon:

```stata
preserve
keep if inlist(g, 0, 3) & inlist(t, 2, 5)
replace t = cond(t == 2, 1, 2)
replace g = cond(g == 3, 2, 0)
cdid y, ivar(id) time(t) gvar(g) method(reg) baseperiod(universal) gmmweighting(optimal) bstrap(0) cband(0)
matrix b = e(b)
scalar att_cdid = b[1, colsof(b)]
bysort id (t): generate double dy = y[2] - y[1]
regress dy i.g if t == 2
display "difference: " %10.3e abs(att_cdid - _b[2.g])
restore
```

```text
difference:  1.998e-15
```

With covariates, `method(reg)`, `method(ipw)` and `method(dr)` give the corresponding adjusted comparisons. With more dates, the standard chained estimator combines the supported changes under the [chosen parallel-trends assumption](methods/parallel-trends.md), and `cdid_chain` shows which pairs contribute to each effect.

## Stepwise DiD

For cohort 3 and the effect at date 6, Harmon's stepwise estimator adds the comparisons for dates 2 to 3, 3 to 4, 4 to 5 and 5 to 6. Cohort 5 is a control for the first two steps; the later steps use the never-treated units. The loop fits each step with `cdid`, takes the post-treatment cell, the last column of `e(b)`, and, as a check, computes the same difference in mean changes directly:

```stata
scalar total_cdid = 0
scalar total_means = 0
forvalues s = 3/6 {
    preserve
    quietly keep if inlist(t, `s' - 1, `s') & (g == 3 | g == 0 | g > `s')
    quietly replace g = cond(g == 3, 2, 0)
    quietly replace t = t - `s' + 2
    quietly cdid y, ivar(id) time(t) gvar(g) method(reg) baseperiod(universal) gmmweighting(optimal) bstrap(0) cband(0)
    matrix b = e(b)
    scalar total_cdid = total_cdid + b[1, colsof(b)]
    bysort id (t): generate double dy = y[2] - y[1]
    quietly summarize dy if t == 2 & g == 2
    scalar m1 = r(mean)
    quietly summarize dy if t == 2 & g == 0
    scalar total_means = total_means + m1 - r(mean)
    restore
}
display "stepwise sum: " total_cdid "   difference: " %10.3e abs(total_cdid - total_means)
```

```text
stepwise sum: 2.3419796   difference:  8.882e-16
```

This reproduces the consecutive-link point estimate. On an incomplete panel, use the target units eligible for the requested horizon and the same observed controls in each step. Uncertainty for the sum must account for the covariance across steps, because units contribute to several comparisons. A single chained fit can also combine longer pairs and additional control pools, and its estimate can then differ from the stepwise sum.

## Cohort-Chained DiD

Balla-Elliott and Norwich apply the chaining argument across later-treated control cohorts whose pretreatment windows are short: one untreated cohort supplies the calendar change from dates 1 to 3 and another the change from 3 to 5, so a stacked regression recovers the change from 1 to 5 through the date effects even if no untreated cohort spans it. In `cdid`, `crosspair` lets the target cohort and the controls supply changes over different pairs, `ptassumption(ctrl)` imposes the shared untreated path among eligible controls, and `gmmweighting(within)` gives the regression weighting.

The test `tests/stata/test_crosspair_stacked.do` checks this correspondence in Stata. With not-yet-treated controls, `baseperiod(universal)`, no covariates and no correction for missing outcomes, the `cdid` estimates equal, to within $10^{-11}$ (largest observed difference $2.5 \times 10^{-13}$), the indicator coefficients of stacked target-cohort regressions fitted with `areg` (unit fixed effects) on a short-window panel: with complete cohort windows, where cohort fixed effects (`regress`) also match; with sampling weights; with monotone attrition within cohorts, where only unit fixed effects match; with dates missing from the calendar; and with one period of anticipation. The same cells are unavailable in both. See [Cross-pair chaining](methods/cross-pair.md) and the [worked example with short pretreatment windows](examples/short-pretreatment-windows.md).

## Binary static switcher summaries

For a complete panel with current binary status `d` and no covariates, `cdid y, ivar(id) time(time) dvar(d) method(reg)` reports a switcher-weighted average of entry and exit effects, expressed in the treated-minus-untreated direction. This corresponds to the WAS estimator of `did_multiplegt_stat`. That command is not installed in the environment used for this guide, so the match was not run in Stata here. The R documentation of `cdid` reports that, on a complete panel of 384 units and eight dates, R `cdid`, R `DIDmultiplegtSTAT` and Stata `did_multiplegt_stat` give the same average, 1.211393; standard errors depend on the variance normalization. With missing outcomes, `cdid` also uses changes over other supported dates in the static model. See [When treatment can stop and restart](methods/switching.md).

## How software uses incomplete panels

The chained DiD paper formalizes a missing-trends assumption: the units observed at both dates of a pair must be representative of their cohort with respect to the mean change over that pair, which this guide calls missing trends completely at random. It is weaker than outcomes missing completely at random: response can depend on a persistent difference in levels that cancels in changes, but missingness that depends on characteristics related to outcome growth selects the changes.

The usual options of the commands above change which observations enter the calculation. A repeated-cross-section DiD (`csdid, bal(none)`) compares outcome levels from separate samples at each date and needs those samples to represent the population. A pair-balanced DiD (`csdid, bal(pair)`) uses changes of units observed at both dates and needs mean-change representativeness. A linear regression on the available rows (`areg`, `reghdfe`, `jwdid`) links them through its outcome model and relies on the corresponding restrictions on the observed rows. None of these choices fits a model of outcome observation.

`cdid` combines within-unit chaining with an explicit correction for covariate-dependent observation. Under conditional missing trends at random, the observed changes must represent the cohort's mean changes at the same covariate values; the corrections average over the study units, including those with one or both outcomes missing for a pair, and require overlap or prediction support in addition to parallel trends.

| How to handle observation | Option in `cdid` |
| --- | --- |
| Missing trends completely at random | `selection(none)` |
| Weight changes by their probability of observation | `selection(conditional_mtar) selmethod(ipw)` |
| Predict the missing changes | `selection(conditional_mtar) selmethod(imputation)` |
| Combine both | `selection(conditional_mtar) selmethod(dr)` |

These choices combine with `method(reg)`, `method(ipw)` and `method(dr)`, in the adoption and switching models and with the cohort-by-date preset. The imputation of missing changes predicts actual changes, which can include treatment effects; this differs from the untreated counterfactual imputed by `did_imputation`. See [Handling missing outcomes](methods/missing-outcomes.md).

## Beyond these specifications

Dynamic switching with carryover, as in `did_multiplegt_dyn`, and nonlinear cohort-by-date models, such as the Poisson option of `jwdid`, are not available in `cdid`; the `dvar()` model targets a static effect. Synthetic control and synthetic DiD (`sdid`), and methods built on an estimated event study, such as the sensitivity analysis of Rambachan and Roth, start from different designs or targets.

## References

- Balla-Elliott, D. and Norwich, I. (2025). Cohort-Chained DiD: Long-Run Effects with Limited Pre-Treatment Data. <https://pdfs.dballaelliott.com/ccdid.pdf>
- Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained Difference-in-Differences. *Journal of Econometrics*, 248, 105783. <https://doi.org/10.1016/j.jeconom.2024.105783>
- Borusyak, K., Jaravel, X. and Spiess, J. (2024). Revisiting Event-Study Designs: Robust and Efficient Estimation. *Review of Economic Studies*. <https://doi.org/10.1093/restud/rdae007>
- Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with Multiple Time Periods. *Journal of Econometrics*, 225(2), 200-230. <https://doi.org/10.1016/j.jeconom.2020.12.001>
- de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F., Sow, D. and Vazquez-Bare, G. (2026). Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period. arXiv:2201.06898v7. <https://arxiv.org/abs/2201.06898v7>
- Dube, A., Girardi, D., Jordà, Ò. and Taylor, A. M. (2025). A Local Projections Approach to Difference-in-Differences. *Journal of Applied Econometrics*. <https://doi.org/10.1002/jae.70000>
- Gardner, J. (2022). Two-Stage Differences in Differences. arXiv:2207.05943. <https://arxiv.org/abs/2207.05943>
- Harmon, N. A. (2024). Difference-in-Differences and Efficient Estimation of Treatment Effects. <https://web.econ.ku.dk/nharmon/docs/harmon2022difference.pdf>
- Sun, L. and Abraham, S. (2021). Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects. *Journal of Econometrics*, 225(2), 175-199. <https://doi.org/10.1016/j.jeconom.2020.09.006>
- Wooldridge, J. M. (2025). Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-Differences Estimators. *Empirical Economics*. <https://doi.org/10.1007/s00181-025-02807-z>

---

[Guide home](README.md) · [Package overview](../README.md)
