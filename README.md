# cdid: chained difference-in-differences for Stata <img src="docs/assets/logo.png" align="right" width="120" alt="cdid logo" />

`cdid` brings widely used difference-in-differences estimators into a common framework for complete and incomplete panels. It reproduces documented specifications of established methods and extends them with chained comparisons, observation weighting, change imputation and doubly robust corrections for selective missingness, by leveraging the chained difference-in-differences (chained DiD).

[User guide](docs/README.md) · [Release notes](NEWS.md) · [R package](https://github.com/dvbn/cdid) · [Chained DiD paper (arXiv)](https://arxiv.org/abs/2301.01085) · [Author website](https://www.davidbenatia.com/research/)

This guide describes `cdid` **1.0.0.1** for Stata. The first Stata release was based on [R `cdid` 1.3.0.0](https://github.com/dvbn/cdid/tree/v1.3.0.0). See [Correspondence with R](#correspondence-with-r) and the [release notes](NEWS.md).

## Installation

`cdid` requires Stata 19.5 (StataNow) or later.

### Install from SSC

Version **1.0.0.0** is available directly from SSC. Run the following in Stata.

```stata
ssc install cdid
mata: mata mlib index
```

### Install from GitHub

Version **1.0.0.1** is available from the GitHub repository. Run the following in Stata.

```stata
net install cdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main")
mata: mata mlib index
```

Chained mode, the default, needs nothing else. Ordinary mode (`mode(ordinary)` or `rcs`) uses `csdid` 2.0.2 by Callaway, Rios-Avila and Sant'Anna, distributed with `cdid` and installed separately so that an existing `csdid` is never replaced without your consent:

```stata
net install csdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main/csdid")
```

[Installation](docs/installation.md) explains how to keep another `csdid` version, install into a project folder, and copy the example do-files with `net get cdid`.

## Get started

The example simulates a staggered panel: 600 units over 2001 to 2008, cohorts first treated in 2004 or 2006, others never treated, and one outcome in five missing. `G` holds each unit's first treatment date, with 0 for units never treated; it is the cohort variable `gvar()`.

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

cdid y x, ivar(id) time(year) gvar(G)
estat dynamic
cdid_plot
```

`cdid` reports one row per cohort and date: post-treatment cells and pretreatment cells, with standard errors and 95% simultaneous bands from the multiplier bootstrap, followed by a joint Wald test of the pretreatment cells. `estat dynamic` averages the cells by time since adoption without refitting, and `cdid_plot` draws the event study. The defaults are those of R's `att_gt()`: doubly robust treatment adjustment (`method(dr)`), never-treated controls, varying base periods, no anticipation, and the multiplier bootstrap with 1,000 draws and simultaneous bands. `csdid` defaults to not-yet-treated controls instead. Typing `cdid` alone replays the current results.

[Get started](docs/get-started.md) shows the output of these commands and how to read it, the four aggregations, export and saved results, and a comparison with ordinary DiD. The [hands-on do-file](examples/cdid_hands_on.do) runs a longer session, and the [user guide](docs/README.md) orders the remaining pages as a reading path.

## What the package does

A treatment cohort is the set of units first treated in the same period. For cohort g and period t, ATT(g,t) is the average effect of the treatment in period t on the units of cohort g. [Callaway and Sant'Anna (2021, *Journal of Econometrics*)](https://doi.org/10.1016/j.jeconom.2020.12.001) define these effects for staggered adoption, estimate each one by comparing the cohort with units never treated or not yet treated, and aggregate them into summaries such as event studies. Their estimators, implemented by `did` in R and `csdid` in Stata, adjust for covariates by outcome regression (OR), inverse probability weighting (IPW) or a doubly robust (DR) combination of the two.

`cdid` extends this framework to incomplete panels with the chained DiD of [Bellégo, Benatia and Dortet-Bernadet (2025, *Journal of Econometrics*)](https://doi.org/10.1016/j.jeconom.2024.105783). Ordinary estimation handles an unbalanced panel either by keeping only the units observed at every date or, with `allowunbalanced`, through repeated cross-section (RCS) estimators that compare outcome levels. Chained DiD uses the changes of each unit between the dates at which it is observed and combines all available short comparisons into the effects by the generalized method of moments (GMM). The targets remain the ATT(g,t) of the cohorts of the original study population, including their units with missing outcomes, and the assumptions concern the missing changes rather than the missing levels. The one exception concerns covariates: under the default `missingreference(exclude)`, a cohort whose treated units lack a required covariate at the reference date is estimated for its eligible treated units (see [Prepare your data](docs/help/prepare-data.md)).

### Why changes within units

Write the outcome of unit i at date t as $Y_{it} = \alpha_i + u_{it}$, where $\alpha_i$ is a fixed individual effect and $u_{it}$ holds everything that varies over time: trends, the treatment effect and noise. The change between two observed dates, $Y_{it} - Y_{is} = u_{it} - u_{is}$, contains no $\alpha_i$. Suppose treated units with persistently low outcomes are more likely to leave the panel while controls leave at random. The treated units still observed later have a higher mean $\alpha_i$, so comparing their outcome levels with earlier ones adds this change in composition to the treatment effect. If attrition depends only on $\alpha_i$ and $\alpha_i$ is unrelated to $u_{it} - u_{is}$, the treated units that remain still have the mean change of their cohort.

This is the main advantage of chained DiD for incomplete panels: observation may depend on the unobservable individual effect, as long as that effect is unrelated to the changes, and the observed levels need not stay representative of the population. Removing the fixed effect does not by itself make missing changes ignorable. If the persistent component is related to trends or to heterogeneous treatment effects, observation that depends on it also selects the changes.

### The assumptions

Within a cohort, the units observed at both dates of a comparison must be representative of that cohort with respect to the mean outcome change over that pair, directly or conditional on covariates with an appropriate adjustment. This equates population expectations of changes. It does not require every unit to have the same change or realized sample means to agree. Parallel trends, no anticipation and overlap must also hold. [Handling missing outcomes](docs/methods/missing-outcomes.md), [Parallel trends](docs/methods/parallel-trends.md) and [Treatment adjustment](docs/methods/treatment.md) state these assumptions.

Neither approach is better in every design. When the sample at each date represents the population, comparisons of levels can be valid too. [Identify effects when the observed units change composition](docs/examples/identification.md) works through a design in which changes identify the effects and levels do not. Differencing removes the persistent unit component, while the transitory noise of each short change accumulates along a chain, so which approach is more precise depends on the design, as Bellégo, Benatia and Dortet-Bernadet (2025, Section 2.1.3 and Proposition 1) analyze. [How to compare precision](docs/examples/precision.md) explains how to compare them on matched targets.

### Two modes

- **Chained mode** is the default for panel data. It estimates the usable short comparisons between pairs of dates at which units are observed and combines them by GMM. An effect that these comparisons cannot recover, given cohort and control support and the assumptions, is reported as [unavailable with a reason](docs/help/troubleshooting.md#an-effect-is-unavailable), never as zero.
- **Ordinary mode**, requested with `mode(ordinary)` and used with `rcs`, estimates the effects of Callaway and Sant'Anna (2021) with `csdid` 2.0.2. [Moving from did and csdid](docs/from-did.md) explains the differences between the modes.

## Data and model decisions

The decisions below come in the order you meet them. Each links to the guide page that states the assumption behind the choice.

### 1. Prepare the panel and the unit table

Supply long data with one row per unit and date: a unit identifier `ivar()`, a numeric date `time()`, a numeric outcome with missing values kept as missing, and the cohort variable `gvar()`. Never replace a missing outcome with zero or a carried-forward value. The units in the data form the study population. To include study units with no observed outcome, pass a frame with one row per unit through `frame()`. Units observed once or never supply no change but stay in the cohort shares, the observation models and the standard errors.

Covariates that change over time are read at each cohort's reference date from the dated rows. A unit can contribute its change from years 3 to 4 even if its year-2 outcome is missing. If a comparison uses a covariate from year 2 and you know its value, supply it on the year-2 row with a missing outcome. If the covariate is unknown, `cdid` warns and excludes the unit from that cohort's comparisons under the default `missingreference(exclude)`. The unit can still contribute to other cohorts for which its required covariates are known. If treated units are excluded, the effect describes the remaining treated units. Specify `missingreference(error)` to stop and list the missing inputs. Sampling weights that change within a unit need an explicit `fixweights()` rule.

[Prepare your panel and unit table](docs/help/prepare-data.md) · [Add a known reference covariate row](docs/help/prepare-data.md#add-a-known-reference-covariate-row) · [Investigate missingness and data structure](docs/help/investigate-missingness.md)

### 2. Choose the comparison group and the parallel trends assumption

`controlgroup()` selects never-treated or not-yet-treated controls, as in `csdid`. `ptassumption(post)`, the default, assumes parallel trends from each cohort's reference period onward. `ptassumption(all)` adds the same assumption from a chosen start date `ptstart()` up to each reference period, so that earlier dates can also serve as baselines. `ptassumption(ctrl)` keeps PT-post and adds common untreated trends among the clean control cohorts only. Valid added restrictions can make more effects available or improve precision, but neither gain is guaranteed, and false restrictions bias the estimates. `baseperiod()` changes only how cells before treatment are reported, not the assumption.

When a cohort and its eligible controls are observed over different windows and share no date pair, `crosspair` links their separate within-unit changes through calendar dates where the trend assumption applies. `cdid_chain` shows which comparisons identify each effect and the weights that combine them. The default `gmmweighting(optimal)` estimates the moment covariance from the data; `gmmweighting(within)` weights comparisons by their within-unit covariance under uncorrelated level shocks.

[Parallel trends](docs/methods/parallel-trends.md) · [Add PT-all restrictions](docs/examples/pt-all-comparison.md) · [Cross-pair chaining, PT-ctrl and within-unit weighting](docs/methods/cross-pair.md) · [Short pretreatment windows](docs/examples/short-pretreatment-windows.md) · [GMM weighting and covariance regularization](docs/methods/gmm.md)

### 3. Adjust for treatment

`method()` makes each cohort comparable with its controls by IPW (`ipw`), OR (`reg`) or DR (`dr`, the default). In chained mode, OR and DR regress the controls' outcome changes on covariates. With the default `outcomepooling(separate)`, each distinct comparison has its own coefficients; `outcomepooling(pooled)` shares slopes per calendar step across comparisons with local intercepts. `method(custom)` accepts a treatment estimator written in Mata.

[Treatment adjustment](docs/methods/treatment.md) · [Choose IPW, OR or DR in an incomplete panel](docs/examples/treatment-choice.md) · [Custom estimators](docs/custom-estimators.md)

### 4. Adjust for missing outcomes, separately

`selection()` states the assumption about missing outcomes. The default `selection(none)` fits no observation model and requires the observed units of each pair to be representative of their cohort's mean change. `selection(conditional_mtar)` allows observation to depend on covariates and adjusts for it. `selmethod()` then chooses observation weighting (`ipw`), which weights each observed change by the inverse of its fitted probability of being observed at both dates, change imputation (`imputation`), which predicts the changes with regressions fitted to the observed changes, or the augmented correction (`dr`), which combines the two. All three work with treatment IPW, OR and DR, and `selcovars()` can give the observation model its own covariates.

Treatment DR with the augmented correction is doubly robust in both layers: for the treatment, either the propensity model or the outcome regression may be misspecified, and for the missing outcomes, either the observation model or the imputation model, provided the other model of each layer and the remaining assumptions are correct. Augmenting the missing-outcome layer does not give treatment IPW or OR the treatment-DR guarantee.

[Handling missing outcomes](docs/methods/missing-outcomes.md) · [Choose an observation model](docs/examples/observation-models.md) · [Impute actual changes](docs/examples/change-imputation.md) · [Sampling information for the observation model](docs/methods/sampling-information.md) · [Robustness](docs/examples/robustness.md)

### 5. Choose the treatment design

The examples above use an absorbing treatment: units adopt once and stay treated, and effects may evolve with time since adoption. Choose the switching model for entry and exit effects, or the linear fixed-effects specification for an adoption panel.

| Design | What to supply | What you estimate | Aggregate |
| --- | --- | --- | --- |
| Adoption cohorts (default) | `gvar()`, the first treatment date, with 0 for units never treated | ATT(g,t) for each cohort and date, so effects can evolve with time since adoption | `estat`, including event studies |
| Binary treatment that can stop and restart | `dvar()`, the treatment status at every study date, including dates with missing outcomes | The immediate effect at each entry or exit, under a static model with no carryover of past treatment | `e(att)`, a switcher-weighted average |
| Linear cohort-by-date fixed effects | `gvar()` and `preset(cohort_time_fe)` | Cohort-by-date effects in a linear specification with unit fixed effects and a common untreated calendar path | `estat` |

The `dvar()` model assumes that past treatment has no carryover beyond current status. Identification also depends on parallel trends and support for each entry or exit comparison. Entrants are compared with units that stay untreated and exiters with units that stay treated, and both estimates are reported on the same treated-minus-untreated scale. Treatment and missing-outcome adjustments are chosen as above. Supplying both `gvar()` and `dvar()` is an error, and absorbing data supplied through `dvar()` still request the static model.

With observation weighting or the augmented correction, the specified observation logits need finite estimates. For example, a separate logit may be separated if everyone with a given covariate value has an observed outcome or change. See [Troubleshooting](docs/help/troubleshooting.md).

[When treatment can stop and restart](docs/methods/switching.md) · [Cohort-by-date fixed effects](docs/methods/cohort-time-fe.md) · [Options reference](docs/reference/options.md)

## Diagnostics and inference

- **Availability.** The table reports the available cells. `cdid_export, type(cells)` writes every requested cell with its availability and the reason an unavailable cell has no estimate. An aggregate that needs an unavailable cell is unavailable unless `estat` is given `narm`. See [Diagnose errors and unavailable results](docs/help/troubleshooting.md).
- **Identification paths.** `cdid_chain` lists the comparisons behind each effect and the weights used to combine them. See [Cross-pair chaining](docs/methods/cross-pair.md).
- **Fitted models.** For a chained fit, specify `keepnuisance`, then type `cdid_nuisance` to see the fitted equations, coefficient tables and standard errors of the treatment, outcome, observation and imputation models. `cdid_diagnostics` summarizes the sample, the settings and the inference.
- **Pretrend tests.** The fit reports a joint Wald test of the cells before treatment when its covariance allows it. `cdid_pretest` runs the conditional test. Failure to reject does not establish parallel trends. See [Inference](docs/examples/inference.md).
- **Uncertainty.** `bstrap()`, `cband()`, `biters()` and `level()` follow `did` and `csdid`: multiplier bootstrap standard errors and simultaneous bands by default, analytical standard errors with `bstrap(0)`. `cluster()` accepts one cluster variable at or above the unit. Standard errors use unit-level influence functions that include the estimation of the fitted models. `noinffunc` returns point estimates only. See [Results](docs/reference/results.md).

## User guide

The [guide home](docs/README.md) orders these pages as a reading path with stages.

| Topic | Pages |
| --- | --- |
| Start with your panel | [Get started](docs/get-started.md), [prepare your data](docs/help/prepare-data.md), [inspect missingness](docs/help/investigate-missingness.md), [hands-on guide](docs/hands-on.md) |
| Work with gaps | [Incomplete panels](docs/examples/incomplete-panels.md), [changing sample composition](docs/examples/identification.md), [short pretreatment windows](docs/examples/short-pretreatment-windows.md) |
| Choose comparisons | [Treatment adjustment](docs/methods/treatment.md), [compare IPW, OR and DR](docs/examples/treatment-choice.md), [parallel trends](docs/methods/parallel-trends.md), [PT-all example](docs/examples/pt-all-comparison.md), [cross-pair chaining](docs/methods/cross-pair.md) |
| Adjust for missing outcomes | [Overview](docs/methods/missing-outcomes.md), [observation models](docs/examples/observation-models.md), [change imputation](docs/examples/change-imputation.md), [sampling information](docs/methods/sampling-information.md), [robustness](docs/examples/robustness.md) |
| Choose the treatment design | [When treatment can stop and restart](docs/methods/switching.md), [cohort-by-date fixed effects](docs/methods/cohort-time-fe.md) |
| Read results | [Options](docs/reference/options.md), [results](docs/reference/results.md), [working with results](docs/user-results.md), [GMM weighting](docs/methods/gmm.md), [inference](docs/examples/inference.md), [compare precision](docs/examples/precision.md), [troubleshooting](docs/help/troubleshooting.md) |
| Compare software | [Familiar DiD models for incomplete panels](docs/comparisons.md), [estimate familiar DiD models with cdid](docs/from-other-packages.md), [moving from did and csdid](docs/from-did.md), [moving from R cdid](docs/from-r-cdid.md), [custom estimators](docs/custom-estimators.md) |

Each command has a help file: `help cdid`, `help cdid postestimation`, `help cdid_dvar`, `help cdid_pretest`, `help cdid_chain`, `help cdid_nuisance`, `help cdid_diagnostics`, `help cdid_plot`, `help cdid_export`, `help cdid_load`, `help cdid_pretest_load` and `help cdid_custom`.

## Correspondence with R

`cdid` for Stata follows the estimators and defaults of [R `cdid` 1.3.0.0](https://github.com/dvbn/cdid/tree/v1.3.0.0). The [R-to-Stata guide](docs/from-r-cdid.md) lists the supported option mappings and interface differences. `cdid` corresponds to `att_gt()`, `estat` to `aggte()`, `cdid_plot` to `ggdid()`, `cdid_pretest` to `conditional_did_pretest()`, `cdid_chain` to `chain_report()` and `cdid_nuisance` to `nuisance_models()`. Estimates and standard errors were compared with those of R on the configurations of the Stata test suite and agree there up to numerical precision; agreement is not established beyond them, and multiplier-bootstrap draws come from Stata's random-number generator and differ from R's. Ordinary mode uses `csdid` 2.0.2 where R `cdid` delegates to `did`. [Moving from R cdid](docs/from-r-cdid.md) maps the R arguments to Stata options and lists those without a counterpart in either direction.

## References

- Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained Difference-in-Differences. *Journal of Econometrics*, 248, article 105783. <https://doi.org/10.1016/j.jeconom.2024.105783>. arXiv version: <https://arxiv.org/abs/2301.01085>
- Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with Multiple Time Periods. *Journal of Econometrics*, 225(2), pages 200 to 230. <https://doi.org/10.1016/j.jeconom.2020.12.001>
- de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F., Sow, D. and Vazquez-Bare, G. (2026). Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period. Working paper, arXiv:2201.06898v7. <https://arxiv.org/abs/2201.06898v7>

Chained IPW (`method(ipw)`) and observation weighting for missing outcomes follow Bellégo, Benatia and Dortet-Bernadet (2025). The package extends these methods with chained OR and DR, change imputation, the augmented correction, PT-all and PT-ctrl, cross-pair chaining, within-unit GMM weighting and the linear cohort-by-date preset. The `dvar()` model adapts the static binary specialization of de Chaisemartin, D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare to complete and incomplete panels. The methods pages of the guide state the assumptions each extension needs.

## Questions and bug reports

Use [Issues](https://github.com/dvbn/cdid-stata/issues) to ask about using the package or choosing a method and to report reproducible bugs. Include the package version, the Stata version, your operating system and a self-contained example.

## Support the project

If cdid is useful to your work, you can [buy me a coffee](https://buymeacoffee.com/dvbn).

## Contributing

See [Contributing](CONTRIBUTING.md) for reporting bugs, suggesting changes and working on a scoped issue.
