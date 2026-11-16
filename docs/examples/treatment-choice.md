# Choose IPW, OR or DR in an incomplete panel

When baseline covariates $X_i$ influence both which cohort a unit joins and how untreated outcomes evolve over time, unadjusted difference-in-differences comparisons are confounded by covariate differences across cohorts. In `mode(chained)`, `cdid` provides three treatment-adjustment methods: inverse probability weighting (`method(ipw)`), regression of the controls' outcome changes (`method(reg)`), and doubly robust estimation (`method(dr)`), with `outcomepooling(separate|pooled)` governing slope sharing across comparison blocks for `method(reg|dr)` ([Treatment adjustment](../methods/treatment.md)).

The model of cohort membership below is a theoretical example; it does not describe how the simulated data of the hands-on do-file assign cohorts.

## Question

Holding the incomplete panel, the study-population frame `frame(cdid_pop)`, the covariates `x_bin x_cont`, `ptassumption(post)`, `baseperiod(varying)`, and `selection(none)` fixed, how do `method(ipw)`, `method(reg)`, and `method(dr)` construct their nuisance fits, and how do their cell and aggregate estimates compare?

## Model algebra for the two nuisance models

### 1. Control outcome-change regression $m_b(X)$ (`method(reg)` and `method(dr)`)

Suppose untreated potential outcomes follow

$$Y_{i,t}(0) = \alpha_i + \delta_t + t\, X_i'\beta_0 + \varepsilon_{i,t}, \qquad E[\varepsilon_{i,t} \mid G_i, X_i, S_i] = 0$$

Across any ordered calendar pair $e = (s, t)$ on a consecutive calendar $t \in \{1, \dots, T\}$ where $L(e) = i_t - i_s = t - s$, the conditional mean untreated change is:

$$E[Y_{i,t}(0) - Y_{i,s}(0) \mid X_i, G_i \in K] = (\delta_t - \delta_s) + L(e)\, X_i'\beta_0$$

Under default `outcomepooling(separate)`, `cdid` fits a separate regression $m_b(X_i) = \alpha_b + L(e)\, z(X_i)'\beta_b$ on each comparison block $b = (e, K)$, which accommodates both a constant per-step slope $\beta_b = \beta_0$ and block-varying covariate trend slopes $\beta_b$. Under `outcomepooling(pooled)`, `cdid` imposes the constant-per-step restriction across all blocks simultaneously:

$$m_b(X_i) = \alpha_b + L(e)\, z(X_i)'\beta$$

with local block intercept $\alpha_b = \delta_t - \delta_s$ for each comparison block $b = (e, K)$ and shared slope vector $\beta = \beta_0$.

### 2. Risk-set treatment propensity logit $p_{g,K}(X)$ (`method(ipw)` and `method(dr)`)

Suppose cohort assignment across $G_i \in \{0, 3, 5\}$ follows a multinomial logit with linear predictors $\eta_g(X_i) = a_g + X_i'\gamma_g$ relative to $G_i = 0$ ($\eta_0(X_i) = 0$):

- **Under `controlgroup(nevertreated)`:** Every comparison for target cohort $g \in \{3, 5\}$ uses the single control cohort $K = \{0\}$. Conditional on belonging to the risk set $G_i \in \{g, 0\}$, the log-odds of belonging to cohort $g$ versus control $K = \{0\}$ are:
  $$\log\frac{P(G_i = g \mid X_i,\, G_i \in \{g, 0\})}{P(G_i = 0 \mid X_i,\, G_i \in \{g, 0\})} = a_g + X_i'\gamma_g$$
  Because this conditional log-odds is exactly linear in $X_i$, the binary logit $p_{g,\{0\}}(X_i)$ fitted by `cdid` is exact for every cohort $g$.
- **Under `controlgroup(notyettreated)`:** For cohort $g = 3$ on pairs $e = (s, t)$ with $\max(s, t) \le 4$ (when `anticipation(0)`), the control pool is $K = \{0, 5\}$. Conditional on $G_i \in \{3, 0, 5\}$, the log-odds of $G_i = 3$ versus $G_i \in \{0, 5\}$ are:
  $$\log\frac{P(G_i = 3 \mid X_i)}{P(G_i = 0 \mid X_i) + P(G_i = 5 \mid X_i)} = a_3 + X_i'\gamma_3 - \log\!\left(1 + \exp(a_5 + X_i'\gamma_5)\right)$$
  This is strictly linear in $X_i$ when $\gamma_5 = 0$ (or when $X_i$ is binary/saturated), and otherwise is a smooth log-sum-exp function of $X_i$. Under `method(dr)`, even if the pooled risk-set logit is only an approximation, consistency is preserved whenever the control outcome regression $m_b(X_i)$ is correctly specified.

## Stata comparison and diagnostics

Run the three estimations of [Treatment adjustment](../methods/treatment.md#5-comparing-the-built-in-methods-in-stata), which are Sections 3a and 3d of the hands-on do-file. After each, `cdid_diagnostics, detail` reports the number of treatment models and pooled regressions fitted. These numbers depend on the comparisons available and on the options.

The multinomial example above gives sufficient conditions for correct models, not finite-sample guarantees. A binary pair logit can be saturated in a binary covariate, while a combined not-yet-treated pool generally has nonlinear log-sum-exp odds. DR is consistent when either the propensity or the outcome model is correct only under valid trend, observation, overlap and regularity conditions.

---

**Reading path, step 8 of 22** · Previous: [Treatment adjustment](../methods/treatment.md) · Next: [Handling missing outcomes](../methods/missing-outcomes.md)

[Guide home](../README.md) · [Package overview](../../README.md)
