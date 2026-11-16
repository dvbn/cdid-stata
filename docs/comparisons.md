# Familiar difference-in-differences models, with incomplete panels

A survey follows households for three years. Some households are interviewed in years 1 and 2, others in years 2 and 3, and few in both years 1 and 3. The two shorter changes still say something about the longer interval. Chained difference-in-differences uses those changes to estimate treatment effects when the observation patterns and the parallel trends assumptions support the comparisons. [Worked example](examples/identification.md) · [Chained DiD paper](https://doi.org/10.1016/j.jeconom.2024.105783)

`cdid` relates several familiar estimators to one Stata command. Its ordinary mode computes the Callaway and Sant'Anna estimator with `csdid`. Its cohort-by-date preset fits the linear regression with unit and date fixed effects and one effect per treated cohort and date, the specification of Wooldridge's extended two-way fixed effects (`jwdid`, `xthdidregress twfe`); with never-treated controls it is also the saturated regression of Sun and Abraham (`eventstudyinteract`), and with not-yet-treated controls and a common untreated path it gives, on a complete panel, the cohort-by-date effects of the imputation and two-stage estimators (`did_imputation`, `did2s`). Consecutive two-date fits give Harmon's stepwise DiD, cross-pair chaining gives the stacked regressions of Cohort-Chained DiD, and the `dvar()` model gives static binary switcher effects (`did_multiplegt_stat`). [Matched specifications and the checks run for this guide](from-other-packages.md)

## What missing trends at random means

Taking changes within households removes persistent differences in their outcome levels. The chained DiD paper asks whether the households observed at both ends of a pair are representative of their cohort with respect to the **mean change** over that interval. This guide calls that baseline missing trends completely at random. It is weaker than outcomes missing completely at random: households with persistently low income can respond less often while their income changes remain representative. If response depends on a characteristic that also predicts income growth, the observed changes can be selective. [Chained DiD paper, Assumptions 2 and 6](https://doi.org/10.1016/j.jeconom.2024.105783)

Repeated cross-section DiD compares levels in the samples observed at different dates; pair-balanced DiD compares changes for the units observed at both dates of each comparison. Each needs its own representativeness assumption. Chaining can use additional observed changes, and the observation correction addresses selection into those changes. [How the commands use incomplete panels](from-other-packages.md#how-software-uses-incomplete-panels)

When representativeness holds **conditional on measured characteristics**, `cdid` can weight changes by their probability of observation, predict the missing changes, or combine both. This choice is separate from the treatment adjustment by outcome regression, inverse probability weighting or doubly robust estimation; for example, the cohort-by-date preset can be combined with doubly robust treatment adjustment and the augmented correction for missing outcomes. [Handling missing outcomes](methods/missing-outcomes.md) · [Cohort-by-date example](methods/cohort-time-fe.md#adjust-for-selective-observation)

The `dvar()` model is for settings where past treatment has no lingering effect. Dynamic treatment histories and nonlinear cohort-by-date models are outside the current package. To try the supported specifications, start with the [matched commands](from-other-packages.md), and use [the data guide](help/investigate-missingness.md) to see which changes your panel supplies.

See the [references for the matched methods](from-other-packages.md#references).

---

[Guide home](README.md) · [Package overview](../README.md)
