# Prepare your panel and unit table

Before running `cdid`, check that the data in memory and the optional unit table in `frame()` have the structure the estimators expect: identifiers, dates, cohorts, covariates and sampling weights. To check the patterns of missing outcomes and the coverage of pairs of dates, see [Investigate missingness and data structure](investigate-missingness.md).

## If treatment can stop and restart

This page describes staggered adoption with `gvar()`. For a binary treatment that units enter and leave, supply `dvar()`, the treatment status at every date, with one row for every unit at every date, including dates with missing outcomes; see [When treatment can stop and restart](../methods/switching.md#what-your-input-needs).

## Core requirements

`cdid` expects long data with one observation per unit and date:

| Element | Option | Requirement |
| --- | --- | --- |
| Outcome | `depvar` | Numeric. A missing value (`.`, `.a` to `.z`) is an unobserved outcome, like an omitted row. |
| Unit identifier | `ivar()` | Numeric or string in chained mode, numeric in ordinary mode. Nonmissing, and unique with the date (`isid id year`). |
| Date | `time()` | Numeric and finite. Dates need not be consecutive (for example 2014, 2016, 2018). Durations used by pooled outcome regressions count positions in the sorted list of study dates, not differences between date values. |
| Cohort | `gvar()` | First treatment date, constant within unit, nonnegative and nonmissing. Code units never treated as `0`, not as missing. |

By default, every treated cohort whose treatment date is on or before the last study date is a target. A cohort treated at or before the first study date has no pretreatment date, and its effects are reported as unavailable.

### Never impute zero for missing outcomes

Keep a missing outcome missing. Replacing it with zero, a carried-forward value or an interpolation creates changes that were never observed, and no option of `cdid` can undo that. Rows with a missing outcome may be kept or dropped: both representations define the same observed outcomes. Keep them when they carry covariates or weights needed at that date (see [Reference dates](#reference-dates-for-covariates-that-change-within-a-unit)).

## The unit table (`frame()`) and unit observation roles

### Study population versus the data in memory

The study population is the set of units whose effects the analysis describes. Without `frame()`, the units present in the data form it, so a unit whose rows were all dropped disappears from it. `frame(name)` supplies a Stata frame with one row per study unit, including units with one or no observed outcome:

| Option | Default | Requirement |
| --- | --- | --- |
| `frame(name)` | none | Existing frame with exactly one row per study unit. |
| `frameivar(name)` | the `ivar()` name | Unique, nonmissing identifier; every unit of the estimation sample must appear in the frame. |
| `framegvar(name)` | not read | Cohort in the frame; must agree with `gvar()` for units in the data. |
| `framecovars(varlist)` | not read | Time-invariant covariates in the frame. |
| `frameselcovars(varlist)` | not read | Time-invariant observation-model covariates in the frame. |
| `framewvar(name)` | not read | Positive unit weights in the frame. |
| `framecluster(name)` | not read | Cluster identifiers in the frame. |

Only `frameivar()` has a default. The other columns are read from the frame only when you name them, even if the frame has variables with the same names as in the data.

### Roles of units with zero, one or repeated observations

`cdid` distinguishes three nested sets of units: the study population $F$; the units with at least one observed outcome, $B$; and the units with at least two, $A$, with $A \subseteq B \subseteq F$. `cdid_diagnostics` reports their sizes. A unit observed once or never supplies no change, but it belongs to its cohort: it enters the cohort shares used by `estat`, the observation models fitted with `seluniverse(F)`, the averages of imputed changes and the influence functions ([Handling missing outcomes](../methods/missing-outcomes.md)).

## Panel mode and `allowunbalanced`

Chained mode, the default with panel data, uses every unit with gaps; `allowunbalanced` changes nothing there. Ordinary mode keeps only units observed at every date unless `allowunbalanced` is specified, in which case it uses repeated cross-section estimators; `rcs` declares repeated cross sections. See [Moving from did and csdid](../from-did.md#incomplete-panels-in-ordinary-and-chained-mode).

## Static versus changing covariates and sampling weights

### Static covariates and weights

Covariates and weights that are constant within units need no date. They must be nonmissing for every unit used in estimation; with `frame()`, supply them for every unit of the frame. Factor-variable notation (`i.region`, `c.x##c.z`) is accepted in the covariates and in `selcovars()`.

### Reference dates for covariates that change within a unit

A cohort's reference date is the latest study date before its treatment date, less the `anticipation()` periods: the latest date $d$ with $d + \text{anticipation} < g$. Covariates that vary over time are read at the reference date of each target cohort for its post-treatment comparisons, and at the earlier date of each pretrend comparison. With `crosspair`, every comparison of cohort $g$, including those of its control cohorts, reads covariates at $g$'s reference date. Covariates are never carried forward, substituted, imputed or read from another date. A known covariate on a row whose outcome is missing is valid input.

When a required covariate is missing at the date where a cohort's comparisons read it, `missingreference()` decides:

- `missingreference(exclude)`, the default, warns and excludes the affected units from that cohort's comparisons only, in treated and comparison roles. The unit still contributes to other cohorts whose reference covariates it has. It keeps its row in the influence functions, and the number of units is unchanged. If treated units are excluded, the cohort's effects describe its remaining treated units, and `estat` weights the cohort by those units. A pretrend comparison whose own earlier-date covariate is missing is reported as unavailable; the main comparisons are unchanged, and the joint pretrend Wald test is then not computed.
- `missingreference(error)` stops with error 459 and lists up to five missing values with the unit, the date, the variable and the comparison that needs it.

`missingreference()` covers covariates only. With `fixweights()`, every unit used needs a positive weight at its anchor date, whatever `missingreference()` says; a missing weight there stops `cdid` with error 459 and a message naming the unit, the date and the variable.

When no required value is missing, both settings give the same results. Time-invariant covariates keep their own rule: they must be complete for the units used. `cdid_pretest` always stops when a required covariate is missing.

The example below has 300 units over years 1 to 4, cohorts first treated in years 3 and 4, and never-treated units. Cohort 3 reads covariates in year 2, cohort 4 in year 3.

```stata
clear
set seed 20261116
set obs 300
generate id = _n
generate G = cond(mod(id, 3) == 0, 3, cond(mod(id, 3) == 1, 4, 0))
expand 4
bysort id: generate year = _n
generate x = rnormal() + 0.2 * year
generate y = 0.5 * x + 0.3 * year + (G > 0 & year >= G) + rnormal()
* Unit 3 (cohort 3) has no year-2 outcome, but its year-2 covariate is known:
* the row stays, with y missing and x supplied.
replace y = . if id == 3 & year == 2
* Unit 6 (cohort 3) has an unknown year-2 covariate.
replace x = . if id == 6 & year == 2
```

```stata
cdid y x, ivar(id) time(year) gvar(G) bstrap(0) cband(0)
matrix list e(reference_eligibility)
```

```text
Warning: Excluded units with missing reference covariates from the affected cohort comparisons.
         Counts across the supplied units at each reference date, including units outside the
         comparison groups: cohort 3 (date 2): 1 treated, 0 never-treated, 0 from other cohorts.
         See e(reference_eligibility) for details. If you know a missing covariate, supply it on
         its unit-date row, with a missing value for a missing outcome.

Chained difference-in-differences            Number of units       =       300
Method: doubly robust (dr)                   Number of cells       =         6
GMM weighting: optimal
Control group: never treated; base period: varying
------------------------------------------------------------------------------
   Group     Time  Role             ATT  Std. err.     [95% conf. interval]
------------------------------------------------------------------------------
       3        2  pre        -.1787712   .2121536    -.5945847    .2370422
       3        3  post        1.171452   .2155867       .74891    1.593994
       3        4  post          1.1487   .2301844     .6975466    1.599853
       4        2  pre        .06346298   .2005939    -.3296938    .4566197
       4        3  pre       -.09102793   .2350894    -.5517947    .3697388
       4        4  post        1.090736   .2152742     .6688067    1.512666
------------------------------------------------------------------------------
Std. err.: analytic; pointwise 95% intervals.
Pretrend Wald test (pre-treatment cells): chi2(3) = 1.79, Prob > chi2 = 0.6163
Restricted comparisons for cohort(s) 3. See e(reference_eligibility).
Effects for cohort 3 cover 99 of 100 treated units.
```

```text
e(reference_eligibility)[2,8]
          cohort  reference_~e     n_treated  n_treated_~e  n_treated_~d  n_never_tr~d
g3             3             2           100            99             1             0
g4             4             3           100           100             0             0

    n_other_co~d    restricted
g3             0             1
g4             0             0
```

Unit 3 is used: its year-2 covariate is known, and its changes between the years in which its outcome is observed enter the comparisons. Unit 6 is excluded from cohort 3's comparisons. In `e(reference_eligibility)`, `n_treated_excluded` counts the cohort-3 units that lack a required covariate in year 2; the effects of cohort 3 describe the `n_treated_eligible` remaining units. `n_never_treated_excluded` and `n_other_cohort_excluded` count units that cannot serve as comparisons for cohort 3; the counts include units outside the comparison groups. Cohort 4 is unrestricted. [Results](../reference/results.md#covariates-missing-at-a-reference-date) lists the other stored results.

The strict rule stops instead:

```stata
capture noisily cdid y x, ivar(id) time(year) gvar(G) bstrap(0) cband(0) missingreference(error)
```

```text
Missing required reference covariate values (time-varying covariates are never carried forward, substituted or deleted): id 6, date 2, variable x, support (1->2, cohort 3) of post link 1 (target 3). Supply each missing covariate on its unit-date row; the outcome may be missing on that row.
```

### Add a known reference covariate row

If you know a missing covariate, supply it on its unit-date row, with a missing outcome if the outcome is missing. If the row exists, fill the covariate. If the row is absent, add it, with every other required field (cohort, weights, clusters). Continuing the example:

```stata
* Suppose unit 6's year-2 covariate becomes known from another source.
local x_year2 = 0.42
replace x = `x_year2' if id == 6 & year == 2
* Unit 9 (cohort 3) has no year-2 row at all, and its year-2 covariate is known.
drop if id == 9 & year == 2
local x9_year2 = 0.15
set obs `=_N + 1'
replace id = 9 in L
replace year = 2 in L
replace G = 3 in L
replace y = . in L
replace x = `x9_year2' in L
cdid y x, ivar(id) time(year) gvar(G) bstrap(0) cband(0)
display "`e(reference_status)'"
```

The fit now has no restricted cohort (`e(reference_status)` is `unrestricted`), and no warning is displayed. For a variable such as age, a recorded date of birth can supply the reference value if it determines age at the date the analysis uses. If a value is unknown, the default excludes the unit from the comparisons of the cohort that needs it.

### Changing sampling weights and `fixweights()`

Supply positive sampling weights with `[iweight = w]` or `wvar(w)`, not both, and with `framewvar()` when `frame()` is used. Chained estimation uses one weight per unit. If weights vary within units, declare the rule with `fixweights(first_period)`, which uses each unit's weight at the first study date, or `fixweights(base_period)`, which uses its weight at each target cohort's reference date. The weight must be observed at that date; `cdid` does not use the weight of another date. Without `fixweights()`, weights that vary within units are refused. Weights are not allowed with `preset(cohort_time_fe)` or `dvar()`.

## Outcome covariates versus observation covariates

The varlist (or `covars()`) holds the covariates of the treatment and outcome models and of change imputation. `selcovars()` gives the observation logits of `selection(conditional_mtar)` their own covariates; by default they are the outcome covariates. `selconstant` requests intercept-only observation logits, `selnoconstant` removes the intercept, and `selzero` uses no observation design. See [Choose an observation model](../examples/observation-models.md).

## Separation, collinearity and covariance

- **Separation.** A treatment or observation logit with no finite solution, for example because a covariate value perfectly predicts treatment within a comparison, stops the fit or makes the affected cells unavailable. Coarsen the covariate, pool the slopes (`selpooling()`), or remove it.
- **Collinearity.** Exactly redundant covariates in an outcome or imputation regression are dropped when they are exactly collinear on every row used to fit and to predict, as in R; other rank deficiencies are refused with error 459 and a message naming the regression.
- **Covariance regularization.** `condition()` and `ridge()` regularize a nearly singular moment covariance in the GMM weighting. This concerns weighting, not identification: it does not make an unidentified effect available ([GMM weighting](../methods/gmm.md)).

## Resource controls

`workspace()` limits the memory used, `chunk()` the number of rows processed at once and `workcap()` the counted work. A fit that would exceed a limit is refused rather than truncated. `posting(state)` keeps results too large for Stata matrices; see [Working with results](../user-results.md).

## Checklist

On the data of the example above:

```stata
assert !missing(id, year, G)
isid id year
bysort id (year): assert G == G[1]
count if missing(y)
count if missing(x)
```

`count if missing(x)` counts the observations with a missing covariate, which may be needed at a reference date. With a unit table, also check it inside the frame, for example `frame units: isid id`.

---

**Reading path, step 6 of 22** · Previous: [How to compare precision](../examples/precision.md) · Next: [Treatment adjustment](../methods/treatment.md)

[Guide home](../README.md) · [Package overview](../../README.md)
