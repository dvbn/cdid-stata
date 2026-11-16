# When treatment can stop and restart

With `dvar(varname)` in place of `gvar()`, `cdid` estimates the effects of a binary treatment that units may enter, leave, and enter again. `varname` records the treatment status $D_{it} \in \{0, 1\}$ of unit $i$ at date $t$. This page defines the estimated cells and their average, the adjustment methods, the handling of missing outcomes and the inference. The syntax is in `help cdid_dvar`; the design corresponds to `att_gt(dname = ...)` of the R package `cdid` 1.3.0.0.

## First treatment date or current treatment status?

Use `gvar()` when treatment is absorbing and you want effects that may evolve with time since adoption, aggregated by `estat`. Use `dvar()` when units can leave and re-enter treatment and the effect of interest is the immediate effect of current status, under the assumption that past treatment has no lingering effect. Supplying both is an error, and absorbing data supplied through `dvar()` still request the static model and its targets.

## What your input needs

`dvar()` requires one row for every unit at every date, a treatment status, 0 or 1, in every row, including rows whose outcome is missing, and nonmissing covariates at $t-1$ for the units at risk of switching at $t$. Outcomes may be missing. Covariates may vary over time. Weights, `saving()` and the options for comparison groups, anticipation, base periods, parallel-trends restrictions and cross-pair chaining do not apply; `missingreference(exclude)` is refused.

## Assumptions and cells

The design assumes that outcomes depend on the current treatment status only, $Y_{it} = Y_{it}(D_{it})$ (no carryover), and that, conditional on the covariates at the earlier date, units that change status and units that keep it would have followed the same untreated or treated outcome path.

For a date $t$ after the first and a prior status $d \in \{0, 1\}$, cell $c = (t, d)$ is defined on the risk set $R_c = 1\{D_{i,t-1} = d\}$:

- switchers $G_c = R_c\, 1\{D_{it} = 1 - d\}$ (entries for $d = 0$, exits for $d = 1$);
- stayers $C_c = R_c\, 1\{D_{it} = d\}$.

With $\Delta Y_{ic} = Y_{it} - Y_{i,t-1}$ and covariates $X_{ic} = X_{i,t-1}$, the cell contrast is

$$
\theta_c = E[\Delta Y_{c} \mid G_c = 1] - E\big[\, E[\Delta Y_{c} \mid C_c = 1, X_c] \mid G_c = 1 \big] = E[Y_t(1-d) - Y_t(d) \mid G_c = 1],
$$

and the reported effect is $\tau_c = \sigma_c \theta_c$ with $\sigma_c = 1 - 2d$. Every reported effect is therefore the effect on the switchers of being treated rather than untreated at date $t$: for entries, $E[Y_t(1) - Y_t(0) \mid G_c = 1]$; for exits, the same contrast for the units that leave treatment. `e(cells)` reports both $\theta_c$ (`theta`) and $\tau_c$ (`att`).

## Adjustment methods

Let $\pi_c(x)$ be the probability of switching given $R_c = 1$ and $X_c = x$, fitted by a logit on the switchers and stayers of the cell, and $r_c(x) = \pi_c(x) / (1 - \pi_c(x))$ its odds. Let $m_c(x) = \alpha_c + x'\beta$ be a least-squares regression of the stayers' outcome changes on the covariates.

| Method | Estimator of $\theta_c$ |
| --- | --- |
| `method(ipw)` | $\bar{\Delta}_{G_c} - \sum_{C_c} r_c(X_i)\, \Delta Y_{ic} \big/ \sum_{C_c} r_c(X_i)$ |
| `method(reg)` | mean over the switchers of $\Delta Y_{ic} - m_c(X_{ic})$ |
| `method(dr)` (default) | mean over the switchers of $\Delta Y_{ic} - m_c(X_{ic})$, minus the odds-weighted mean over the stayers of the same residual |

With `outcomepooling(separate)` (default), each cell has its own regression. With `outcomepooling(pooled)`, the regression has an intercept $\alpha_c$ for each requested cell and slopes $\beta$ common to all requested cells, fitted on their stayers. A partition such as `outcomepooling(entry:2003 entry:2004 | exit:2003 exit:2004)` fits common slopes $\beta_b$ within each block $b$ of cells, again with an intercept for each cell; every requested cell must belong to one block. A block may include cells that are not requested (for example the exit cells with `switchers(entry)`), whose stayers then contribute only to the fit of $\beta_b$; this is available in complete panels only. Covariates may vary over time; each cell uses their values at $t - 1$.

## Missing outcomes

When outcomes are missing, the change of each role from $t-1$ to $t$ is chained through the dates at which its units are observed. For cell $c$ and role $a \in \{G, C\}$, every pair of dates $s < u$ at which at least one unit of the role has observed outcomes is an edge of the role graph, with moment

$$
\hat\psi^a_{c,(s,u)} = \frac{\sum_i A_{ic} O_{is} O_{iu} h_{ic} (Y_{iu} - Y_{is})}{\sum_i A_{ic} O_{is} O_{iu} h_{ic}},
$$

where $A_{ic}$ is the role indicator, $O_{is}$ indicates an observed outcome, and $h_{ic} = r_c(X_{ic})$ for the stayers under `method(ipw)` and `method(dr)` ($h_{ic} = 1$ otherwise). The levels of every role at its observed dates, relative to $t-1$, are estimated jointly for all requested cells by GMM, with `gmmweighting(optimal)` (default; Tikhonov regularization controlled by `condition()` and `ridge()`) or `gmmweighting(identity)`; each cell contrast is the switchers' level at $t$ minus the stayers' level at $t$.

With `method(reg)` or `method(dr)`, the stayer regression is reconstructed from the graph: the stayers' mean change $m_c$ and its covariance with the covariates combine the observed edges along the least-squares path from $t-1$ to $t$, the slopes solve the pooled normal equations with the covariate moments of all stayers, and $\alpha_c = m_c - \bar z_c'\beta$. Every edge that spans the switch ($s < t \le u$) is then centered by the predicted change $\alpha_c + \bar z'\beta$, evaluated at the switchers' covariate mean for the switchers and at the odds-weighted stayers' covariate mean for the stayers; under `method(reg)` only the switchers' links enter the GMM.

The design assumes that the units of a role observed at a pair of dates have the same expected outcome change as the whole role. A cell is unavailable (`status` 3 or 4 in `e(cells)`) when the observed dates of its switchers or of its stayers do not connect $t-1$ and $t$.

### Observation models

With `selection(conditional_mtar)`, the probability that outcomes are observed may depend on the covariates at $t-1$ within each role. Let $S_{i,(s,u)} = O_{is} O_{iu}$, $q_{ic,(s,u)}$ its fitted probability, and $j_{ic,(s,u)}$ the predicted change from a regression of the observed changes on the covariates. The edge moments become

| `selmethod()` | Weight on unit $i$ | Change |
| --- | --- | --- |
| `ipw` (default) | $A_{ic} h_{ic} S_{i,(s,u)} / \max(q_{ic,(s,u)}, 10^{-6})$ | $Y_{iu} - Y_{is}$ |
| `imputation` | $A_{ic} h_{ic}$ | $j_{ic,(s,u)}$ |
| `dr` | $A_{ic} h_{ic}$ | $j_{ic,(s,u)} + S_{i,(s,u)} (Y_{iu} - Y_{is} - j_{ic,(s,u)}) / \max(q_{ic,(s,u)}, 10^{-6})$ |

`selmodel(pair_logit)` (default) fits $q$ by a logit for each pair; `selmodel(independent_period_logit)` fits a logit for each date and multiplies the two probabilities; it assumes that, given the covariates, observation at the two dates is independent. Every logit and imputation regression has an intercept for each role and pair (or date) and slopes shared within a block: `separate` (each role and pair), `within_group` (each role of a transition), or `pooled` (all), set by `selpooling()` and `imppooling()`. With `impscale(per_period)` (default), the covariates of a pair are multiplied by its number of periods. The observation logits use `selcovars()` (default: the outcome covariates) or an intercept only with `selconstant`. `seluniverse(B)` or `seluniverse(A)` fits the pair logits on the units with at least one or at least two observed outcomes, with an optional known ratio `selrho()` constant within units; it requires `selmethod(ipw)`, `selmodel(pair_logit)` and `seljustification()`. The influence functions of all edge moments include the estimation of these models.

### Fitted models

With `keepnuisance`, `cdid_nuisance` reports the switcher logits of the transitions (`treatment`), the stayer outcome-change regressions (`outcome`), the observation logits (`selection`) and the imputation regressions (`imputation`), with standard errors from the same original-unit influence functions as the effects.

## Requested cells and the average

`switchers(both|entry|exit)` selects the directions. By default, the requested cells are those with at least one switcher in the requested directions. `switchtimes(numlist)` requests the cells at the listed dates; a listed cell without switchers or stayers is reported as unavailable.

The average weights the effects of the requested cells by their numbers of switchers $N_c$:

$$
\delta = \sum_{c} \omega_c \tau_c, \qquad \omega_c = \frac{N_c}{\sum_{c'} N_{c'}}.
$$

It is reported in `e(att)` only when every requested cell is estimable; it is never renormalized over a subset of cells.

## Inference

The influence function of each $\hat\tau_c$ accounts for the estimated logit and regression; that of $\hat\delta$ also accounts for the estimated switcher shares. Analytical standard errors sum the influence functions within clusters (`cluster()`; units otherwise), center the cluster sums, and divide by the number of units. With `bstrap(1)` (default), the multiplier bootstrap draws one Rademacher sign per cluster and reports the interquartile-range standard error; with `cband(1)` (default), the bounds of the cells are simultaneous, from the bootstrap distribution of the maximum $t$-statistic over the estimable cells. The average uses the same draws, has a pointwise interval, and does not enter the band. `e(V)` holds the analytical covariance of the cells in `e(b)`.

## Example

```stata
clear
set seed 2027
set obs 300
generate long id = _n
generate double x = rnormal()
expand 6
bysort id: generate year = 2000 + _n
generate byte treated = runiform() < invlogit(-1 + .5*x)
generate double y = 0.3*(year - 2000) + 0.5*x + treated + rnormal()

* entries and exits at every date, doubly robust
cdid y x, ivar(id) time(year) dvar(treated)

* entries only, outcome regression with common slopes, analytical inference
cdid y x, ivar(id) time(year) dvar(treated) switchers(entry) method(reg) outcomepooling(pooled) bstrap(0) cband(0)

* cells at two dates and a comparison of entry and exit effects
cdid y x, ivar(id) time(year) dvar(treated) switchtimes(2003 2004) bstrap(0) cband(0)
lincom _b[entry:t2003] - _b[exit:t2003]
cdid_export using cells.dta, type(cells) replace

* outcomes missing at some dates: the changes are chained through the observed dates
replace y = . if runiform() < .15
cdid y x, ivar(id) time(year) dvar(treated) bstrap(0) cband(0)

* augmented correction using x, and the fitted models
cdid y x, ivar(id) time(year) dvar(treated) selection(conditional_mtar) selmethod(dr) selpooling(pooled) imppooling(pooled) keepnuisance bstrap(0) cband(0)
cdid_nuisance
cdid_nuisance selection
```

## After estimation

`cdid_export` saves `type(cells)` and `type(aggregate)`, and `cdid_nuisance` reports the fitted models when `keepnuisance` is specified. `estat`, `cdid_plot`, `cdid_chain` and `cdid_diagnostics` are not available. The stored results are listed in the [results reference](../reference/results.md#results-of-cdid-with-dvar).

## References

- de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F., Sow, D. and Vazquez-Bare, G. (2026). Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period. Working paper, arXiv:2201.06898v7. <https://arxiv.org/abs/2201.06898v7>

The `dvar()` model adapts their static binary specialization to complete and incomplete panels.

---

[Guide home](../README.md) · [Package overview](../../README.md)
