# Handling missing outcomes

In an incomplete panel, only a subset of units is observed at both endpoints of any ordered calendar pair $e = (s, t)$. If the units observed on pair $e$ differ systematically in their outcome **changes** $\Delta_e Y = Y_t - Y_s$ from the rest of their cohort in the study population $F$, unadjusted cohort pair averages no longer represent the full study population.

`cdid` separates **the identifying assumption about missing outcomes** (`selection()`) from **the estimator and the models used to adjust for missing outcomes** (`selmethod()`, `selmodel()`, `selpooling()`, `imppooling()`).

## The study population

Every chained `cdid` fit targets the average treatment effect over the **study population $F$** (`frame(name)`, or all units in the estimation sample when `frame()` is omitted), except that under `missingreference(exclude)` a cohort whose treated units lack a required covariate at the reference date is estimated for its eligible treated units ([Prepare your data](../help/prepare-data.md)):

- **$F$ (`e(N_full)`):** All study units, including singletons observed at only one date ($B \setminus A$) and unit-table-only units observed at zero dates ($F \setminus B$).
- **$B$:** Units observed at $\ge 1$ date.
- **$A$ :** Units observed at $\ge 2$ dates.

Even though a unit in $F \setminus A$ (observed once or never) supplies no raw change $\Delta_e Y_i$, it remains a member of its cohort in $F$: its baseline covariates in `frame()` enter observation-logit fitting on `seluniverse(F)`, actual-change imputation predictions on $F$, cohort weights in `estat`, cluster sums, and the complete $N_F$-row influence functions.

## Without an observation model (`selection(none)`, default)

With `selection(none)`, `cdid` fits no observation-probability model $q$ and no change-imputation model $j$. For each cohort $k$ and date pair $e = (s, t)$, observed changes on $S_e = 1$ are scaled by the cohort's empirical pair response rate on $F$ and passed directly to treatment adjustment (`method(dr|reg|ipw)`).

This requires **unconditional mean-change representativeness within each cohort**:

$$E[Y_t - Y_s \mid G = k, S_{s,t} = 1] = E[Y_t - Y_s \mid G = k]$$

As shown in [Identify effects when the observed units change composition](../examples/identification.md), this condition concerns outcome **changes** $Y_t - Y_s$, not outcome levels $Y_t$ or unit fixed effects $\alpha_i$. However, if a covariate $W_i$ shifts both the probability of observing pair $e$ ($S_{s,t} = 1$) and the outcome change $Y_t - Y_s$ within cohort $k$, `selection(none)` is generally biased for the population $F$.

## Conditional mean trends missing at random (`selection(conditional_mtar)`)

Setting `selection(conditional_mtar)` invokes the **Conditional Mean Trends Missing at Random** assumption (Corollary 1 of Bellégo, Benatia and Dortet-Bernadet 2025): conditional on covariates $W_i$ and cohort $G_i = k$, the mean outcome change over pair $e = (s, t)$ is the same among units observed at both endpoints ($S_{s,t} = 1$) and in the full cohort $G_i = k$ on $F$:

$$E[Y_t - Y_s \mid W, G = k, S_{s,t} = 1] = E[Y_t - Y_s \mid W, G = k]$$

This condition must hold for **every** cohort $k$ (both treated and control cohorts, including actual post-treatment changes of treated cohorts) and **every** date pair $e = (s, t)$ used in the comparisons.

> [!NOTE]
> **Distinguishing the assumption from the estimator:** `selection(conditional_mtar)` states the conditional mean-independence assumption above. How `cdid` estimates the cohort population changes under that assumption is chosen separately by `selmethod(ipw|imputation|dr)`.

## Three estimators under `selection(conditional_mtar)`

### 1. Observation weighting (`selmethod(ipw)`, default)

`selmethod(ipw)` estimates the probability $q_c(W_i) = P(S_e = 1 \mid G_i = k, W_i)$ of observing both endpoints of pair $e = (s, t)$ for each needed cell $c = (k, e)$, and reweights observed changes ($S_{i,e} = 1$) by $1 / \max(\hat{q}_c(W_i), 10^{-6})$ (combined with the cohort's empirical pair rate and any transport factor `selrho()`).

You choose the parametric observation model $q_c(W)$ with two options:

- **`selmodel(pair_logit|independent_period_logit)`:**
  - `pair_logit` (default): directly fits a binary logit for the pair indicator $S_{i,e} = R_{i,s} R_{i,t} \in \{0, 1\}$.
  - `independent_period_logit`: fits binary logits for single-date observation $R_{i,s}$ and $R_{i,t}$ and multiplies the two fitted date probabilities: $\hat{q}_{k,(s,t)}(W_i) = \hat{r}_{k,s}(W_i)\,\hat{r}_{k,t}(W_i)$. This requires conditional independence of the two endpoint observation indicators given cohort and the scientifically sufficient covariates. It requires `seluniverse(F)`.
- **`selpooling(separate|within_cohort|across_cohorts)`:**
  All three pooling choices retain a **separate intercept for every cell** (whenever an intercept is included) and impose different restrictions on the covariate slope vector of $q$:
  - `separate` (default): separate slope vector for every cell $(k, e)$ (under `pair_logit`) or every cohort-date $(k, d)$ (under `independent_period_logit`).
  - `within_cohort`: one shared slope vector per cohort $k$ across its date pairs or dates.
  - `across_cohorts`: one global slope vector shared across all cohorts and pairs/dates.

Under `selmethod(ipw)`, `imppooling()` does not apply (`e(inactive_options)` lists it, together with `outcomepooling()` under `method(ipw|custom)`), and no imputation model $j$ is fitted.

### 2. Change imputation (`selmethod(imputation)`)

Instead of fitting observation probabilities $q_c$, `selmethod(imputation)` fits a weighted least squares linear model $j_c(X_i)$ for the **actual outcome change** $\Delta_{i,e} = Y_{i,t} - Y_{i,s}$ on units of cohort $k$ observed at both dates ($S_{i,e} = 1$):

$$j_{k,e}(X_i) = \alpha_{k,e} + L(e)\, z(X_i)'\beta_r$$

where $L(e) = i_t - i_s$ is the step count in the stored calendar, $\alpha_{k,e}$ is a local cell intercept (unless `noconstant` is specified), and $\beta_r$ is the slope vector of imputation pooling group $r$ controlled by **`imppooling(separate|within_cohort|across_cohorts)`**:

- `imppooling(separate)` (default): separate slope vector $\beta_{k,e}$ for every cohort and date pair cell $(k, e)$.
- `imppooling(within_cohort)`: one shared slope vector $\beta_k$ per cohort across all of its date pairs.
- `imppooling(across_cohorts)`: one global slope vector $\beta$ shared across all cohorts and date pairs.

For each cell $(k, e)$, `cdid` evaluates $\hat{j}_{k,e}(X_i)$ for **every unit of cohort $k$ in the full study population $F$** (including units whose change was observed, singletons, and unobserved unit-table units) and averages those predictions over $F$ before applying treatment adjustment (`method(dr|reg|ipw)`).

- Under `selmethod(imputation)`, `selmodel()`, `selpooling()`, and `selcovars()`/`frameselcovars()` do **not** apply (`e(inactive_options)` lists them, together with `outcomepooling()` under `method(ipw)`), and no observation logit is fitted.
- Pure imputation requires `seluniverse(F)` and does not allow transport weights (`selrho()`).

### 3. Augmented doubly robust correction (`selmethod(dr)`)

`selmethod(dr)` combines the actual-change imputer $j_c(X_i)$ and the observation-probability model $q_c(W_i)$ on `seluniverse(F)` (without `selrho()`) and works with all three built-in treatment methods: `method(ipw)`, `method(reg)`, and `method(dr)` (`method(custom)` is not supported).

For each needed cell $c = (k, e)$, every unit $i$ of cohort $k$ in $F$ receives the reconstructed change:

$$Z_{i,c} = \hat{j}_c(X_i) + S_{i,e}\,\frac{\Delta_{i,e} - \hat{j}_c(X_i)}{\max(\hat{q}_c(W_i),\, 10^{-6})}$$

The residual term is not renormalized. Built-in treatment adjustment then runs on the reconstructed changes $Z_{i,c}$:

- **`method(ipw)`:** fits $q_c, j_c$ and risk-set treatment logits $p$ (no $m$; `outcomepooling()` does not apply). Consistent if ($q_c$ or $j_c$) and $p$ are correctly specified.
- **`method(reg)`:** fits $q_c, j_c$ and counterfactual control regressions $m_b$ (with `outcomepooling(separate|pooled)`, no $p$). Consistent if ($q_c$ or $j_c$) and $m$ are correctly specified.
- **`method(dr)`:** fits $q_c, j_c, p$, and $m_b$ (doubly robust in both layers). Consistent if ($q_c$ or $j_c$) is correctly specified on the common conditioning set $W$ and ($p$ or $m$) is correctly specified ([Compare model robustness](../examples/robustness.md)).

## Separate covariate lists (`covars()` vs `selcovars()`)

`cdid` lets you specify separate covariate lists for the treatment/imputation design $X$ and the observation-probability design $Q$:

| Option | Models using this design | Default / intercept controls |
| --- | --- | --- |
| `indepvars` or `covars(varlist)`, `framecovars(varlist)` | Treatment logit $p$, counterfactual control regression $m$, actual-change imputer $j$ | Includes an intercept unless `noconstant` is specified. |
| `selcovars(varlist)`, `frameselcovars(varlist)` | Observation-probability logit $q$ (`selmethod(ipw)` or `selmethod(dr)`) | Inherits `covars()` / `framecovars()` when omitted. Intercept controls: `selnoconstant` (omit intercept), `selconstant` (intercept-only `~ 1`), `selzero` (empty design `~ 0`, fixing pair probability to $0.5$ under `pair_logit` or $0.25$ under `independent_period_logit`). |

## Probability safeguards

Two fixed numerical safeguards, the same as in R `cdid`, prevent division by zero in inverse probability weights without dropping observations:

1. **Treatment propensity cap (`1 - 1e-6`):** Fitted risk-set propensities $\hat{p}_{g,K,e}(X_i)$ are capped above at $0.999999$ .
2. **Observation probability floor (`1e-6`):** Fitted pair observation probabilities $\hat{q}_c(W_i)$ entering inverse observation weights (`selmethod(ipw)` and `selmethod(dr)`) are floored below at $10^{-6}$. With independent-period logits, the floor applies to the product for the pair, not to each date. Positive empirical observation rates used as denominators have the same floor. Scores and information matrices use the unbounded probabilities, and the derivatives of the weights account for the bounds where they bind. The number of bounded probabilities is displayed when it was stored.

The robustness statements about the original target require that the bounds do not bind in the population (or other conditions that preserve the target), and that crossings of the bounds be negligible for the usual influence-function argument. Binding bounds can change the target and the robustness properties. The bounds provide no identification and no guarantee of nominal coverage.

## With `preset(cohort_time_fe)`

The cohort-by-date preset applies the same three estimators to every pair of dates $e = (s, t)$ observed for at least one unit of cohort $k$. Each unit $i$ of the cohort in $F$, including units observed at fewer than two dates, receives a corrected change

$$Z_i(k; s, t) = \begin{cases} S_{i,e}\, \Delta_e Y_i / \hat q_{k,e}(W_i) & \text{selmethod(ipw)} \\ \hat j_{k,e}(X_i) & \text{selmethod(imputation)} \\ \hat j_{k,e}(X_i) + S_{i,e}\,\{\Delta_e Y_i - \hat j_{k,e}(X_i)\} / \hat q_{k,e}(W_i) & \text{selmethod(dr)} \end{cases}$$

where $\hat q_{k,e}$ is the pair logit (or the product of the two date logits) and $\hat j_{k,e}(X) = a_{k,e} + (t - s)\, X'b$ is the imputation regression, pooled as set by `selpooling()` and `imppooling()`. With `method(reg)`, the cohort-by-date and date levels are fitted by least squares of $Z_i(k; s, t)$ on the difference of the levels of $t$ and $s$, with the covariates entering as in the regression without missing outcomes; with `method(ipw)` or `method(dr)`, the corrected changes replace the observed changes in the treatment moments, each averaged over all units of its cohort. With `seluniverse(B)` or `seluniverse(A)` (only with `selmodel(pair_logit)` and `selmethod(ipw)`), the denominator is $\hat\pi_k\, \rho_i \max\{\hat q_{k,e}(W_i), 10^{-6}\}$, where $\hat\pi_k$ is the share of cohort $k$ in the fitting roster. The standard errors include the estimation of $\hat q$, $\hat\pi$ and $\hat j$. `keepnuisance` keeps all fitted models for `cdid_nuisance`.

## Summary of the 9 supported built-in combinations under `selection(conditional_mtar)`

| `selmethod()` | `method(ipw)` | `method(reg)` | `method(dr)` | Models fitted | Allowed `seluniverse()` |
| --- | --- | --- | --- | --- | --- |
| `selmethod(ipw)` | Supported | Supported | Supported | $q$ + ($p$, $m$, or $p,m$) | `F` (default); `B` or `A` only for `pair_logit`, with `seljustification()` |
| `selmethod(imputation)` | Supported | Supported | Supported | $j$ + ($p$, $m$, or $p,m$) | `F` only (`selrho()` not allowed) |
| `selmethod(dr)` | Supported | Supported | Supported | $q, j$ + ($p$, $m$, or $p,m$) | `F` only (`selrho()` not allowed) |

Other assumptions about missing outcomes, such as sequential missingness with permanent dropout, are not implemented in `cdid`.

---

**Reading path, step 9 of 22** · Previous: [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md) · Next: [Sampling information for the observation model](sampling-information.md)

[Guide home](../README.md) · [Package overview](../../README.md)
