# Cohort-by-date fixed effects

Use `preset(cohort_time_fe)` for an absorbing treatment and a linear fixed-effects specification of untreated outcomes with unit effects and common date effects. Each treated cohort can have a different effect at each date. The fit uses the changes across each unit's observed dates, including changes that span gaps. The preset corresponds to `control = list(preset = "cohort_time_fe")` of [R `cdid` 1.3.0.0](https://github.com/dvbn/cdid).

## Fit the model

Supply long data with one row per unit and date, an outcome and the cohort variable, coded with the adoption date and zero for never-treated units. Include at least one never-treated cohort, use positive consecutive integer study dates, and place every adoption after the first study date. Rows with missing outcomes can stay in the data, and `frame()` can supply units with no observed outcome.

The example uses the panel of [Get started](../get-started.md): 600 units observed from 2001 to 2008, with one outcome in five missing.

```stata
clear
set seed 20261116
set obs 600
generate id = _n
generate G = cond(runiform() < .4, 0, cond(runiform() < .5, 2004, 2006))
generate x = rnormal()
expand 8
bysort id: generate year = 2000 + _n
generate y = x + .2*(year - 2000) + (G > 0 & year >= G)*(1 + .25*(year - G)) + rnormal()
replace y = . if runiform() < .2
```

```stata
cdid y, ivar(id) time(year) gvar(G) preset(cohort_time_fe) method(reg) baseperiod(universal) controlgroup(nevertreated) bstrap(0) cband(0)
```

```text
Chained difference-in-differences            Number of units       =       600
Method: linear regression with fixed effects Number of cells       =        14
Cohort-by-date fixed effects (within)
Control group: never treated; base period: universal
------------------------------------------------------------------------------
   Group     Time  Role             ATT  Std. err.     [95% conf. interval]
------------------------------------------------------------------------------
    2004     2001  pre        -.1003211   .1614878    -.4168314    .2161892
    2004     2002  pre        .04890602   .1603049    -.2652857    .3630978
    2004     2004  post        .6800819   .1596719     .3671307    .9930331
    2004     2005  post        1.275351   .1669978     .9480411     1.60266
    2004     2006  post        1.345052    .159823     1.031805    1.658299
    2004     2007  post        1.652869   .1643671     1.330715    1.975022
    2004     2008  post        1.981949   .1596912      1.66896    2.294938
    2006     2001  pre          .003069   .1588206    -.3082137    .3143517
    2006     2002  pre         .1196553   .1610641    -.1960245    .4353351
    2006     2003  pre       -.05484205   .1554356    -.3594902    .2498061
    2006     2004  pre        -.3042621   .1624833    -.6227235   .01419928
    2006     2006  post        .7819229   .1453564     .4970295    1.066816
    2006     2007  post        1.095659   .1693169     .7638035    1.427514
    2006     2008  post        1.730351   .1625521     1.411755    2.048948
------------------------------------------------------------------------------
Std. err.: analytic; pointwise 95% intervals.
Pretrend Wald test (pre-treatment cells): chi2(6) = 8.83, Prob > chi2 = 0.1835
```

With never-treated controls, the date immediately before adoption is the reference, and the other pretreatment cells are estimated deviations; the footer tests them jointly when their covariance allows it. With `method(reg)`, this is the linear regression of the outcome on unit fixed effects, date fixed effects and a separate effect for each treated cohort and date, fitted on the observed rows by least squares within units.

For a common untreated path across both never-treated and not-yet-treated cohorts, use `controlgroup(notyettreated) ptassumption(all)` with `ptstart()` at the first study date. That model imposes the common path on all untreated cells from the first date and reports post-treatment effects only; its pretreatment restrictions are part of the fitted model.

```stata
cdid y, ivar(id) time(year) gvar(G) preset(cohort_time_fe) method(reg) baseperiod(universal) controlgroup(notyettreated) ptassumption(all) ptstart(2001) bstrap(0) cband(0)
```

The preset requires `baseperiod(universal)`, `anticipation(0)` and no sampling weights. With `method(reg)`, `gmmweighting()` may be omitted or set to `within`. `cluster()` and the multiplier-bootstrap options apply as usual, and `estat` aggregates the cells. With `selection(none)`, the default, the observed changes must satisfy the model's within-unit mean restriction; see [Handling missing outcomes](missing-outcomes.md) when observation may be selective.

## Adjust for covariates

Covariates in the varlist let the untreated calendar path and each cohort-by-date effect vary with time-invariant covariates: the design has common date-specific slopes and separate cohort-by-date deviations, with `outcomepooling(separate)`. Supply the covariates for every unit of the study population, including units observed once or never. The reported effect averages the fitted contrast over the cohort in the study population, so a unit without an observed change contributes to that covariate average and to the influence functions. The fit reports rank deficiency or numerical singularity of the covariate design, and cells that the observed changes do not connect to the calendar path remain unavailable.

```stata
cdid y x, ivar(id) time(year) gvar(G) preset(cohort_time_fe) method(reg) baseperiod(universal) bstrap(0) cband(0)
estat dynamic
```

## Treatment weighting and doubly robust adjustment

`method(ipw)` reweights clean comparison units to each treated cohort's covariate distribution. `method(dr)` combines a shared untreated change model with reweighted residuals. Both estimate the average untreated change separately for each cohort over its whole study population. Propensity models use the target and donor cohorts and include an intercept. With never-treated controls only the never-treated cohort is a donor; with not-yet-treated controls, each other cohort is a donor for changes that end before its adoption. The DR outcome model fits a shared linear calendar path to clean changes. The default `gmmweighting(optimal)` uses the covariance of the unit moment scores, and `gmmweighting(identity)` gives equal weight to each retained moment.

These fits require overlap between each treated cohort and the clean donors. With `selection(none)`, a sufficient condition is that the observations behind each cohort and date pair preserve its covariate distribution and conditional mean change. Given those conditions and parallel trends, DR is consistent when either the treatment propensity or the untreated change model is correct.

```stata
cdid y x, ivar(id) time(year) gvar(G) preset(cohort_time_fe) method(dr) baseperiod(universal) bstrap(0) cband(0)
```

## Adjust for selective observation

`selection(conditional_mtar)` applies when the observed changes represent the cohort's changes conditional on the covariates of the observation or imputation model. Choose the treatment method with `method()` and the observation correction with `selmethod()`:

| Observation correction | Response for each unit and supported pair of dates |
| --- | --- |
| `selmethod(ipw)` | The observed change divided by its probability of being observed, or zero when either outcome is missing. |
| `selmethod(imputation)` | The fitted conditional mean change for that cohort and pair. |
| `selmethod(dr)` | The fitted conditional mean change plus the observed residual divided by its observation probability. A unit with a missing change contributes the fitted mean. |

`selmodel(independent_period_logit)`, which assumes that observation at the two dates of a pair is independent given the covariates, and both imputation corrections use the whole study population (`seluniverse(F)`); `seluniverse(B)` and `seluniverse(A)` require `selmodel(pair_logit)` and `selmethod(ipw)`, and the probability that a unit of the cohort belongs to the fitting units is then estimated by its share in the cohort. Fitted observation probabilities are bounded below by $10^{-6}$. Covariates must be constant within units and nonmissing for every unit; `missingreference(exclude)` is not available with the preset. The standard errors include the estimation of the observation, imputation, outcome and propensity models.

```stata
cdid y x, ivar(id) time(year) gvar(G) preset(cohort_time_fe) method(dr) selection(conditional_mtar) selmethod(dr) baseperiod(universal) bstrap(0) cband(0) keepnuisance
cdid_nuisance
```

`cdid_nuisance` lists the calendar regression, the treatment logits and, with a correction for missing outcomes, the observation logits and imputation regressions; see `help cdid_nuisance`. `cdid_chain` is not available after the preset, whose estimates do not come from chained comparisons.

## Comparison with regression commands

The documentation of the R package `cdid` reports that the regression fit of its preset matches R `etwfe` and Stata `jwdid` on the same inputs. The Stata preset was compared with R `cdid`, not directly with `jwdid`. Standard errors use the unit score sandwich without a degrees-of-freedom adjustment, so they differ from regression standard errors that apply one. On incomplete panels, check that each command uses the same observations, and compare aggregations only after aligning their cell weights: `estat` weights cohorts by their size in the study population.

## How the changes are combined

For the regression fit, eliminating a unit intercept gives the same least-squares objective as forming consecutive differences along that unit's observed dates and weighting them by the inverse of their joint covariance under uncorrelated level errors with equal variance. A unit observed at dates 1, 3 and 4 contributes changes from 1 to 3 and from 3 to 4, which share the outcome at date 3; the fit accounts for that dependence. Inference uses the empirical unit or cluster scores. Treatment IPW and DR instead combine observed changes and cohort-standardized untreated changes by GMM.

---

[Guide home](../README.md) · [Package overview](../../README.md)
