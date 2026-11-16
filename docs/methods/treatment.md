# Treatment adjustment

When baseline covariates shift both which units enter treatment and how untreated outcomes evolve over time, `cdid` adjusts each two-date comparison for covariates through the `method()` option (and `outcomepooling()` when an outcome regression is fitted):

- `method(dr)` (default): doubly robust adjustment combining treatment propensity weighting and control outcome-change regression,
- `method(reg)`: outcome regression (OR) adjustment using control outcome-change regression,
- `method(ipw)`: inverse probability weighting (IPW) using risk-set treatment propensity logits.

This page defines the treatment adjustment models in `mode(chained)`, explains how `cdid` parameterizes covariates, comparison-block outcome pooling (`outcomepooling(separate|pooled)`), and calendar duration steps $L(e)$, and shows how to specify and compare the methods in Stata.

## 1. Specifying treatment covariates in Stata

The Stata snippets on this page use the simulated panel and the global `CDID_DEMO_BITERS` created by Sections 0 and 1 of the [hands-on do-file](../hands-on.md); run those sections first.

Supply baseline covariates either as `indepvars` immediately after `depvar` or inside `covars(varlist)` (do not pass non-empty lists to both simultaneously), and pass `framecovars(varlist)` when using a population frame `frame(name)`:

```stata
* Equivalent ways to specify covariates x_bin and x_cont:
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) method(dr) controlgroup(nevertreated) baseperiod(varying) biters($CDID_DEMO_BITERS)
cdid y, ivar(id) time(period) gvar(first_treat) covars(x_bin x_cont) method(dr) bstrap(0) cband(0)
```

- **Intercept handling:** By default, both the treatment propensity logits and the control outcome-change regression include local intercepts. Specify `noconstant` to omit the intercepts; the support and rank checks still apply.
- **Factor variables:** Standard Stata factor-variable notation (`i.group`, `c.x_bin#c.x_cont`) is supported.
- **Separate observation-model covariates:** `indepvars` / `covars()` define the design $z(X)$ used by the treatment propensity model $p$, the counterfactual control outcome regression $m$, and, when it is used, the actual-change imputer $j$. Observation-probability logits $q$ under `selection(conditional_mtar)` inherit `covars()` by default or use a separate list via `selcovars(varlist)` ([Handling missing outcomes](missing-outcomes.md#separate-covariate-lists-covars-vs-selcovars)).

## 2. The three treatment adjustment methods

For each comparison $l$ pairing a target cohort $g$, an ordered calendar pair $e = (s, t)$ with $s < t$, and a control cohort set $K$, let $\tilde{\Delta}_{i,e}$ and $\tilde{w}_{i,e}$ denote the outcome change and unit weight produced by the correction for missing outcomes ([Handling missing outcomes](missing-outcomes.md)). Treatment adjustment forms the comparison estimate $\hat{\mu}_l$ from $(\tilde{\Delta}_{i,e}, \tilde{w}_{i,e})$ as follows:

### Inverse probability weighting (`method(ipw)`)

For each distinct risk set $(g, K, e)$, `cdid` fits a weighted binary logit $p_{g,K,e}(X)$ for membership in target cohort $g$ versus control pool $K$, and reweights control changes by the estimated odds $\frac{\hat{p}_{g,K,e}(X_i)}{1 - \hat{p}_{g,K,e}(X_i)}$ using Hajek (self-normalized) weights.

### Control outcome-change regression (`method(reg)`) and `outcomepooling()`

Rather than refitting the same outcome regression for every comparison that shares a control risk set, `mode(chained)` groups the comparisons into blocks $b = (e, K)$: the distinct pairs of a date pair $e = (s, t)$ and a control group $K$. How `cdid` parameterizes the covariate slopes across those blocks is controlled by `outcomepooling(separate|pooled)` (stored in `e(outcome_pooling)`):

- **`outcomepooling(separate)` (default):** Fits a separate weighted least squares regression on each comparison block $b = (e, K)$:

$$m_b(X_i) = \alpha_b + L(e)\, z(X_i)'\beta_b$$

- **`outcomepooling(pooled)`:** Fits one joint weighted least squares regression across all comparison blocks $b = (e, K)$ simultaneously with a shared per-step slope vector $\beta$:

$$m_b(X_i) = \alpha_b + L(e)\, z(X_i)'\beta$$

where:

- $\alpha_b$ is a local intercept specific to block $b = (e, K)$ (present unless `noconstant` is specified), allowing arbitrary baseline untreated time trends across every pair $(s, t)$ and control set $K$,
- $z(X_i)$ is the non-constant covariate vector from `covars()` / `indepvars`,
- $\beta_b$ (`separate`) is a block-specific slope vector for $(e, K)$, whereas $\beta$ (`pooled`) is a single slope vector shared across all blocks $b$,
- $L(e) = i_t - i_s \ge 1$ is the **number of calendar steps** between the index $i_s$ of date $s$ and the index $i_t$ of date $t$ in the sorted stored calendar.

> [!IMPORTANT]
> **Calendar step count $L(e)$ vs raw date subtraction $t - s$:** `cdid` always defines $L(e) = i_t - i_s$ as the number of transitions in the stored calendar, never the arithmetic difference $t - s$ of the numeric labels. For example, if your panel is biennial with stored dates 2014, 2016, 2018 and 2020, the pair $(2014, 2016)$ has $L(e) = 2 - 1 = 1$ step and $(2014, 2018)$ has $L(e) = 3 - 1 = 2$ steps.

After fitting the block regression(s), `method(reg)` subtracts the predicted counterfactual change $\hat{m}_b(X_i)$ from the treated cohort's change on block $b$. Under `method(ipw)` or `method(custom)`, `outcomepooling()` does not apply (`e(outcome_pooling)` is empty and `outcomepooling()` is listed first in `e(inactive_options)`).

### Doubly robust adjustment (`method(dr)`, default)

`method(dr)` combines both models: it subtracts the control outcome-change regression prediction $\hat{m}_b(X_i)$ (fitted under `outcomepooling(separate)` by default, or `outcomepooling(pooled)`) from both treated and control changes and applies inverse propensity odds weights from $\hat{p}_{g,K,e}(X_i)$ to the control residuals. Under `selection(none)` (or when combined with a valid missing-outcome adjustment), `method(dr)` is consistent for each comparison link if **either** the treatment propensity logit $p$ **or** the control outcome-change regression $m$ is correctly specified ([Compare model robustness](../examples/robustness.md)).

## 3. Change imputation uses a different regression

Do not confuse the counterfactual control regression $m_b(X)$ of `method(reg|dr)` with the actual-change imputation model $j_c(X)$ of `selmethod(imputation|dr)` ([Handling missing outcomes](missing-outcomes.md#2-change-imputation-selmethodimputation)):

| Feature | Counterfactual control regression $m_b(X)$ | Actual-change imputer $j_c(X)$ |
| --- | --- | --- |
| Used with | `method(reg)` or `method(dr)` | `selection(conditional_mtar)` with `selmethod(imputation)` or `selmethod(dr)` |
| Fitted on | Control cohorts $G_i \in K$ on comparison blocks $b = (e, K)$ | Observed pairs ($S_e = 1$) within each individual cohort $k$ (both treated and control) |
| Predicts | Counterfactual **untreated** change $\Delta_e Y(0)$ | **Actual** cohort change $\Delta_e Y$ (including treatment effects for treated cohorts) |
| Slope sharing | Controlled by `outcomepooling(separate\|pooled)` (default `separate` per block $b$; `pooled` shares $\beta$ across blocks) | Controlled by `imppooling(separate\|within_cohort\|across_cohorts)` (default `separate`) |

## 4. Propensity safeguard and collinearity handling

- **Propensity score bound (`1 - 1e-6`):** As in R `cdid`, fitted treatment propensities $\hat{p}_{g,K,e}(X_i)$ are capped above at $1 - 10^{-6} = 0.999999$ so that inverse odds $\hat{p}/(1 - \hat{p})$ remain finite without dropping observations. The number of bounded propensities is shown when it was stored with the results; otherwise `cdid_diagnostics` reports it as not retained.
- **Collinearity and rank checks:** In the control outcome-change regression $m$ (whether fitted per block under `outcomepooling(separate)` or jointly under `outcomepooling(pooled)`), exactly collinear covariates are dropped only when the counterfactual predictions do not depend on them, on the support of the block's control and target units; if a required prediction is not identified or a treatment logit has no interior estimate, estimation stops with a support, rank or convergence error (separation: error 459; a boundary estimate or nonconvergence: error 430); `cdid` does not change the estimator ([Diagnose errors and unavailable results](../help/troubleshooting.md)).

## 5. Comparing the built-in methods in Stata

Sections 3a and 3d of the hands-on do-file estimate the same specification with DR (with the default `outcomepooling(separate)` and with `outcomepooling(pooled)`), IPW and OR. The do-file uses few bootstrap draws; the results neither rank the methods by precision nor show that any model is correct:

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) biters($CDID_DEMO_BITERS)

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) outcomepooling(pooled) biters($CDID_DEMO_BITERS)

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(ipw) biters($CDID_DEMO_BITERS)

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(reg) biters($CDID_DEMO_BITERS)
```

See [Treatment choice](../examples/treatment-choice.md) for the algebra of the risk-set logits and the qualifications on shared slopes.

---

**Reading path, step 7 of 22** · Previous: [Prepare your panel and unit table](../help/prepare-data.md) · Next: [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md)

[Guide home](../README.md) · [Package overview](../../README.md)
