# How to compare precision

A smaller standard error does not by itself establish a better estimator. Before comparing chained and ordinary mode, control groups, parallel trends choices or corrections for missing outcomes, match the target population and parameter, then state the assumptions that justify each fit. Balance alone decides neither validity nor a precision ranking. [Moving from did](../from-did.md) explains how to choose between the two estimation modes.

This page explains the covariance of chained estimates and the diagnostics that help compare precision. It reports no simulation results. The examples elsewhere in this guide show single estimations; they do not compare precision over repeated samples.

## More comparisons do not mean independent information

Chained estimation combines two-date comparisons. If a vector of link estimates has covariance matrix $\Sigma$ and its coefficients $a$ are fixed, the variance of the linear combination is $a'\Sigma a$. The off-diagonal covariances matter as much as the individual link variances. This is an exact identity for fixed coefficients. Feasible GMM estimates its weights from the same data: inserting fitted, data-dependent coefficients does not give an exact unconditional finite-sample variance identity. Under the appropriate GMM regularity conditions, a first-order calculation uses the limiting combination and complete influence functions. The package's inference follows its complete influence functions, rather than treating fitted weights as an exact finite-sample variance formula.

For the same unit, consecutive changes telescope:

$$
(Y_{i,4}-Y_{i,3})+(Y_{i,3}-Y_{i,2})=Y_{i,4}-Y_{i,2}.
$$

The shared intermediate outcome cancels. If time-specific measurement errors are mutually independent, their contribution to the covariance of these adjacent changes is minus the variance of the shared error. The total covariance need not be negative: serial dependence, heterogeneous trends and other shared components also matter. With different observed units or weights at each link, estimated link averages need not telescope.

In a balanced panel, chained estimates agree with ordinary DiD up to rounding under the following settings: built-in IPW, no covariates, no correction for missing outcomes, never treated controls, weights fixed across comparisons, PT-post, and matching anticipation and reporting settings. Splitting the same units' long changes into consecutive steps adds comparison rows without independent information. This equality of point estimates does not extend to every standard error, band or aggregation. See [the precise comparison](../from-did.md#when-chained-and-ordinary-estimates-should-agree).

## What GMM weighting and standard errors measure

Overlapping comparison paths can make the moment covariance singular even when an ATT is supported. Lack of support is a different issue: if the available comparisons and trend restrictions cannot recover a cell, weighting cannot make it available. Inspect support and reasons with `cdid_diagnostics, detail` or the `type(cells)` export.

Built-in chained cohort fits default to `gmmweighting(within)`, which uses a working covariance of within-unit changes and accounts for shared dates. With `gmmweighting(optimal)`, feasible two-step GMM uses an estimated moment covariance. `cdid` first removes comparisons that are exactly redundant. With `gmmcovariance(tikhonov)`, it then adds a ridge proportional to the trace when the remaining covariance needs regularization; a covariance that is exactly zero has a separate least-squares rule. Regularization affects the weighting only: it adds no identifying information and does not guarantee a gain in precision. `cdid_diagnostics` reports how the covariance was handled; [GMM weighting](../methods/gmm.md) gives the rules.

The inverse weighting matrix is not the reported estimator covariance. With influence functions computed:

- The influence functions have one row per unit of the study population, including units without an observed change.
- `e(V)` stores analytical estimator covariance, while `e(V_asymptotic) = e(N_full) * e(V)`.
- `e(analytic_se)` is the square root of the diagonal of `e(V)`. `e(se)` holds the reported standard errors, analytical or from the multiplier bootstrap, so it need not equal `e(analytic_se)`.

The influence functions account for the estimation of the adjustment models and, in aggregations, of the cohort shares. Clustering sums their rows within clusters for inference; it does not change the unit-based GMM weights or point estimates. Check the inference settings with `cdid_diagnostics`, and compare pointwise intervals with pointwise intervals, and simultaneous bands with bands of the same aggregation. [Inference](inference.md) and the [results reference](../reference/results.md#uncertainty) explain the endpoints and `estat` outputs.

## Incomplete panels: match the population and sampling conditions

Complete-case analysis, repeated cross sections and chained changes each have settings where they are valid. Their standard errors answer a common precision question only when their targets and assumptions align.

| Approach | What to check before comparing precision |
| --- | --- |
| Complete-case balancing | It uses units observed at every date. Their ATT can differ from the full population ATT when observation relates to heterogeneous effects. The targets can coincide under suitable representativeness conditions; a rotating panel may contain no complete cases. |
| Outcome-level RCS comparisons | Each date's sample must satisfy the level-sampling conditions. Changing composition in an unobserved fixed component can invalidate that comparison even when observed within-unit changes remain representative. Fixed differences alone do not invalidate representative RCS sampling. |
| Chained changes | The observed changes, possibly after correction, must represent the study population under the assumptions chosen. Overlapping pairs can identify effects without any complete unit; when few units link the dates, estimates can still be imprecise. |

Ordinary mode with `allowunbalanced` uses the repeated cross-section estimators while clustering by unit. [Moving from did](../from-did.md) describes the options, and [Handling missing outcomes](../methods/missing-outcomes.md) states the chained conditions.

## Restrictions can change information and robustness

| Choice | What can change, holding the population target fixed |
| --- | --- |
| Never treated versus not yet treated controls | Under not-yet-treated parallel trends and overlap, chained estimation can use the same pair of dates again with other control groups at later horizons. Combining these comparisons by GMM can change estimates and uncertainty even with complete outcomes. More valid comparisons do not guarantee a gain in every finite sample. |
| PT-post versus PT-all | Credible extra untreated-trend restrictions from the chosen `ptstart()` can add identification or precision. False restrictions can bias estimates; neither bias nor a gain follows from the option label alone. PT-all may add no restriction in a particular calendar. `baseperiod()` separately controls how pretreatment cells are reported. |
| Treatment OR or DR (`outcomepooling()`) | Default `outcomepooling(separate)` fits comparison-block-specific counterfactual change slopes $\beta_b$, whereas `outcomepooling(pooled)` imposes one shared slope vector $\beta$ across comparison blocks. Credible pooling can share information; unjustified pooling can invalidate the outcome regression. DR remains consistent if the propensity model is correct, under the other required conditions. |
| Observation or imputation pooling | `selpooling()` shares coefficients of observation-probability models; `imppooling()` shares regressions of actual cohort changes, including treatment effects. Neither argument controls `outcomepooling()` in treatment OR/DR. Their distinct restrictions require justification. |

Small pair probabilities can amplify observation weights; rare links and weak covariate support can worsen uncertainty. Neither a large unit table nor a short chain alone determines precision. [Robustness](robustness.md) separates the model-correctness alternatives, and [Parallel trends](../methods/parallel-trends.md) explains which restrictions are active.

## Read a fit before making a ranking

Use `e(N_full)`, the dataset written by `cdid_export, type(cells)` and `cdid_diagnostics` to compare the population, target cells, support, restrictions and inference actually used. An unavailable effect has no precision ranking. An estimated standard error differs from the standard deviation of estimates over repeated samples. Bias, RMSE and coverage can be assessed only by simulation, with known population values, every replication counted (including those that fail), and Monte Carlo standard errors.

`workspace()` limits the memory used for the main numerical arrays, and `workcap()` optionally limits the number of operations; neither limits the total memory or computation time of Stata.

---

**Reading path, step 5 of 22** · Previous: [Identify effects when the observed units change composition](identification.md) · Next: [Prepare your panel and unit table](../help/prepare-data.md)

[Guide home](../README.md) · [Package overview](../../README.md)
