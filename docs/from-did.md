# Moving from did and csdid

`cdid` keeps the workflow of Callaway and Sant'Anna (2021), implemented in R by `did` and in Stata by `csdid`, and extends it to incomplete panels with the chained difference-in-differences estimator of Bellégo, Benatia and Dortet-Bernadet (2025). Users of R `cdid` will find the argument-by-argument mapping in [Moving from R cdid](from-r-cdid.md).

## What stays the same

Cohorts $g$, dates $t$ and group-time effects $\text{ATT}(g,t)$ have the same meaning. The treatment adjustments (`dr`, `reg`, `ipw`), the control groups (`nevertreated`, `notyettreated`), the base periods (`varying`, `universal`), anticipation, multiplier-bootstrap bands and the four aggregations work as in `did`. The defaults are those of `did`: doubly robust adjustment, never-treated controls, no anticipation, a varying base period and 1000 bootstrap draws with simultaneous bands.

| Concept | Stata `csdid` 2.0.2 | Stata `cdid` |
| --- | --- | --- |
| Group-time estimation | `csdid y x, ivar(id) time(t) gvar(g)` | `cdid y x, ivar(id) time(t) gvar(g)` |
| Estimator family | Ordinary DiD | `mode(chained)` (default) or `mode(ordinary)` |
| Repeated cross sections | omit `ivar()` | `rcs` (implies ordinary mode) |
| Treatment adjustment | `method(dr)` (default), `method(reg)`, `method(ipw)` | `method(dr)` (default), `method(reg)`, `method(ipw)` |
| Comparison group | not-yet-treated (default); `nevertreated` | never-treated (default); `controlgroup(notyettreated)` |
| Base period | universal (default); `base_period(varying)` | varying (default); `baseperiod(universal)` |
| Aggregations | `estat simple\|group\|event\|calendar` | `estat simple\|group\|dynamic\|calendar` |
| Plot | `csdid_plot` | `cdid_plot` |
| Binary treatment with entry and exit | not available | `dvar()`; see `help cdid_dvar` |

## Two estimation modes

- **Chained mode**, the default with panel data, uses the outcome changes of each unit between the dates at which it is observed and combines them by GMM.
- **Ordinary mode**, requested with `mode(ordinary)` or `rcs`, computes the estimators of Callaway and Sant'Anna (2021) with `csdid` 2.0.2, installed separately ([Installation](installation.md)). `cdid` prepares the sample, sets the options of `csdid`, and stores the results in the same form as in chained mode, so `estat`, `cdid_plot`, `cdid_export`, `cdid_load` and `cdid_diagnostics` work in both modes.

Ordinary estimation with panel data requires a numeric `ivar()`. With `rcs`, it may be omitted. Ordinary estimation takes covariates from the varlist and refuses `covars()`, `noconstant`, `panel()`, `frame()` and its associated options. Chained-mode options such as `selection()`, the `sel*()` options and the GMM options have no effect on an ordinary fit. `e(inactive_options)` lists all accepted inactive options, and the note names those supplied in the command.

## Incomplete panels in ordinary and chained mode

By default, ordinary mode keeps only the units observed at every date, as `did` does with a balanced panel, and a note reports how many units were dropped. With `allowunbalanced`, ordinary mode keeps all units and estimates each effect with the repeated cross-section estimators, which compare outcome levels; standard errors are clustered by `ivar()`. Those comparisons are valid when the observed units at each date represent the population at that date.

Chained mode uses every unit observed at least twice and keeps units observed once or never in the study population. It requires the units observed at both dates of a comparison to represent their cohort's mean change, possibly after a correction for missing outcomes ([Handling missing outcomes](methods/missing-outcomes.md)). Fixed unit components cancel from these changes. In [Get started](get-started.md#compare-with-ordinary-did), ordinary mode keeps 98 of 600 units of a panel with one outcome in five missing, while chained mode uses all of them.

## When chained and ordinary estimates should agree

In a balanced panel, chained and ordinary point estimates agree up to rounding with built-in IPW, no covariates, no missing-outcome correction, never-treated controls, PT-post, weights that are constant over time, and the same anticipation and reporting settings: the cohort means then telescope. This equality does not extend in general to OR and DR, standard errors, bands or aggregations.

## Why balanced results can differ

Chained OR and DR regress the controls' outcome changes on covariates, one regression per comparison by default (`outcomepooling(separate)`) or with slopes shared across comparisons (`outcomepooling(pooled)`). Even with separate regressions, chained estimation combines multi-step and direct comparisons by GMM, so identical variable lists do not imply identical fitted systems.

With not-yet-treated controls, chained estimation can also use the same pair of dates again at a later horizon with a different control group. For ATT(3,3) with reference date 2 and cohorts 0, 3 and 5, ordinary estimation uses cohort 3 against `{0,5}` on the pair (2,3). Chained estimation also uses that pair against `{0}` at a later horizon, after cohort 5 is treated. Both comparisons require not-yet-treated parallel trends and overlap at their dates; a claim only about the combined `{0,5}` pool at date 3 does not justify the `{0}` reuse. The additional comparisons can change finite-sample estimates; they do not by themselves make the estimates more precise.

Which to use: with a complete panel and no reason to expect selective observation, both are valid under their assumptions, and agreement in the settings listed above is a useful check. With missing outcomes, chained mode is the safer choice unless observation can be taken to be unrelated to persistent differences between units, provided its own assumptions about the changes hold.

## What chained mode adds

- `frame()` supplies the study population, so that units observed once or never still count in the cohort shares, covariate distributions and influence functions ([Prepare your panel](help/prepare-data.md)).
- `selection(conditional_mtar)` separates the assumption about missing outcomes from the treatment adjustment, with observation weighting, change imputation or the augmented correction (`selmethod()`).
- `ptassumption(all)` and `ptassumption(ctrl)` add parallel trends restrictions before the reference period, independently of how `baseperiod()` reports pretreatment cells ([Parallel trends](methods/parallel-trends.md)).
- `crosspair`, `gmmweighting(within)` and `cdid_chain` link cohorts observed over different windows and show which comparisons identify each effect ([Cross-pair chaining](methods/cross-pair.md)).
- `preset(cohort_time_fe)` fits the linear cohort-by-date fixed-effects specification ([Cohort-by-date fixed effects](methods/cohort-time-fe.md)).
- `cdid_pretest` runs the conditional pretrend test of Callaway and Sant'Anna (2021), with IPW by default; it leaves the current `e()` results unchanged.
- `missingreference()` sets what happens when a covariate that varies over time is missing at a cohort's reference date.

## Ordinary-mode details

- The doubly robust estimator for repeated cross sections, also used for unbalanced panels in ordinary mode, uses the influence function of `DRDID` 1.4.0, which corrected the sign of a pre-period term of earlier versions. Point estimates do not depend on this correction; standard errors do.
- The ordinary conditional pretrend test follows `did` 2.5.1, including its defaults when the model has no covariates or the data have no unit identifier; see `help cdid_pretest`.

---

**Reading path, step 2 of 22** · Previous: [Get started: estimate effects in an incomplete panel](get-started.md) · Next: [Estimate treatment effects from incomplete panels](examples/incomplete-panels.md)

[Guide home](README.md) · [Package overview](../README.md)
