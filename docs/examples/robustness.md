# Compare model robustness across treatment and missing outcome adjustments

When baseline covariates influence which units enter treatment, how untreated outcomes evolve over time, and which outcomes are recorded across dates, chained difference-in-differences (`cdid`, Bellégo, Benatia and Dortet-Bernadet 2025; Callaway and Sant'Anna 2021) separates estimation into two modeling layers: **treatment adjustment** (`method(ipw|reg|dr)`) and **handling missing outcomes** (`selection(none|conditional_mtar)` and `selmethod(ipw|imputation|dr)`).

This article defines the four nuisance models across the two layers, summarizes the large-sample consistency conditions for each supported combination, and walks through five model-correctness scenarios. For full mathematical definitions, see [Treatment adjustment](../methods/treatment.md) and [Handling missing outcomes](../methods/missing-outcomes.md).

## 1. The four labeled nuisance models across two layers

To recover average treatment effects for the study population $F$ in `frame()`, `cdid` uses up to four labeled nuisance models across two layers:

| Symbol | Model name | Layer and Stata options | What it models and where it is fitted |
| --- | --- | --- | --- |
| $p$ | Treatment propensity model | Treatment adjustment (`method(ipw)` or `method(dr)`, using `covars()` / `indepvars`) | For each distinct risk set, a logit model of membership in target cohort $g$ versus `controlgroup()`, weighted according to the correction for missing outcomes. |
| $m$ | Counterfactual untreated change regression | Treatment adjustment (`method(reg)` or `method(dr)`, with `outcomepooling(separate\|pooled)`, using `covars()` / `indepvars`) | Weighted least squares regression fitted on control units per comparison block $b = (e, K)$ by default (`outcomepooling(separate)`: $m_b(X) = \alpha_b + L(e)\, z(X)'\beta_b$) or jointly across blocks (`outcomepooling(pooled)`: $m_b(X) = \alpha_b + L(e)\, z(X)'\beta$). |
| $q$ | Observation pair probability model | Handling missing outcomes (`selection(conditional_mtar)` with `selmethod(ipw)` or `selmethod(dr)`, using `selcovars()`) | For each needed cohort and date pair cell $c = (k, e)$, a logit model $q_c(W)$ for observing both dates of pair $e$, either directly (`selmodel(pair_logit)`) or as the product of two date logits (`selmodel(independent_period_logit)`). |
| $j$ | Actual cohort change imputer | Handling missing outcomes (`selection(conditional_mtar)` with `selmethod(imputation)` or `selmethod(dr)`, using `covars()` / `indepvars`) | For each needed cohort and date pair cell $c = (k, e)$, a weighted least squares linear model $j_c(X) = \alpha_{k,e} + L(e)\, z(X)'\beta_r$ fitted to observed changes of cohort $k$, **treatment effects included**, and evaluated across every unit of cohort $k$ in $F$. |

## 2. Correctness conditions across both layers

Every consistency result below is a large-sample property that requires:

- **Scientifically sufficient conditioning:** Supplying variables in `covars()` or `selcovars()` does not prove sufficiency; the covariates must justify conditional parallel trends and conditional mean trends missing at random ($E[Y_t - Y_s \mid W, G = k, S_{s,t} = 1] = E[Y_t - Y_s \mid W, G = k]$) in the population.
- **Overlap and support in the population, with the probability bounds not binding:** Treatment propensities stay strictly below one, observation probabilities stay strictly positive where $q$ is fitted, imputation cells have positive training mass and prediction support where $j$ is fitted, the fitted models converge to regular interior limits, and all comparisons combined by GMM satisfy valid moment conditions for the same parameter.

| Correction for missing outcomes | `method(ipw)` | `method(reg)` | `method(dr)` |
| --- | --- | --- | --- |
| `selection(none)` with representative pair mean changes (no $q$ or $j$) | $p$ correct | $m$ correct across all blocks | $p$ correct **or** $m$ correct |
| `selection(conditional_mtar) selmethod(ipw)` (fits $q$) | $q$ correct in every cell **and** $p$ correct | $q$ correct in every cell **and** $m$ correct | $q$ correct in every cell **and** ($p$ correct **or** $m$ correct) |
| `selection(conditional_mtar) selmethod(imputation)` (fits $j$) | $j$ correct in every cell/group **and** $p$ correct | $j$ correct in every cell/group **and** $m$ correct | $j$ correct in every cell/group **and** ($p$ correct **or** $m$ correct) |
| `selection(conditional_mtar) selmethod(dr)` (fits $q$ and $j$) | For every needed cell $c$, ($q_c$ correct **or** $j_c$ correct), **and** $p$ correct | For every needed cell $c$, ($q_c$ correct **or** $j_c$ correct), **and** $m$ correct | For every needed cell $c$, ($q_c$ correct **or** $j_c$ correct), **and** ($p$ correct **or** $m$ correct) |

### Four qualifications

1. **Treatment DR supplies no general correction for a wrong imputer:** Under pure imputation (`selmethod(imputation)`), the sufficient condition requires a correct $j$ within each needed cell/pooling group. Choosing `method(dr)` provides the $p$-or-$m$ alternative *after* cohort changes are reconstructed; it does not generally repair biased reconstructed changes from a misspecified $j$. (These conditions are sufficient, not necessary: special projection or aggregation cancellations can still occur.)
2. **Coupling through shared fits:** When `outcomepooling(pooled)` shares one slope vector $\beta$ across all comparison blocks in $m$, and when `imppooling()`, `selpooling()`, or `selmodel(independent_period_logit)` shares parameters across cells or dates, misspecification in one cell or block can shift a shared slope everywhere it is used (whereas default `outcomepooling(separate)` isolates $\beta_b$ to block $b$).
3. **Common conditioning set $W$ and covariates only in `selcovars()`:** Under the augmented correction (`selmethod(dr)`), `covars()` builds $j$, $m$, and $p$ from covariates $X$, whereas `selcovars()` builds $q$ from $W = (X, Z)$. Variables $Z$ listed only in `selcovars()` never enter $j$, $m$, or $p$. If $Z$ shifts the conditional mean of actual changes given $W$ but is omitted from `covars()`, $j_c(X)$ cannot generally represent that $W$-conditional mean, so the condition that the imputer is correct fails; a correct $q_c(W)$ together with valid treatment adjustment on $W$ provides the alternative condition. Adding extra variables to `selcovars()` is therefore not automatically harmless.
4. **Every model used must be estimable:** the logits must pass their rank checks and have interior estimates, and the regressions must pass their rank, conditioning and prediction-support checks, even when the other model of a doubly robust pair is correct; `cdid` does not fall back to another estimator. A completely observed or completely unobserved cell with a free intercept in $q$ stops estimation with `r(498)`. Pure imputation (`selmethod(imputation)`) runs on fully observed cells because it fits no $q$, whereas `selmethod(dr)` stops whenever any model it uses ($q$, $j$, and $p$ or $m$ as required by `method()`) cannot be estimated.

## 3. Five misspecification scenarios

These are theoretical examples of sufficient conditions; they do not describe the simulated data of the hands-on do-file. Consider five specifications on a study population with treatment covariates `x_bin x_cont` and an additional observation-shifting covariate `z_sel`:

- `fit_none_dr`: `selection(none) method(dr)`
- `fit_obs_or`: `selection(conditional_mtar) selmethod(ipw) method(reg) selcovars(x_bin x_cont z_sel)`
- `fit_obs_dr`: `selection(conditional_mtar) selmethod(ipw) method(dr) selcovars(x_bin x_cont z_sel)`
- `fit_imp_dr`: `selection(conditional_mtar) selmethod(imputation) imppooling(separate) method(dr)`
- `fit_aug_dr`: `selection(conditional_mtar) selmethod(dr) imppooling(separate) method(dr) selcovars(x_bin x_cont z_sel)`

| Scenario | Correct vs misspecified models | `fit_none_dr` | `fit_obs_or` | `fit_obs_dr` | `fit_imp_dr` | `fit_aug_dr` |
| --- | --- | --- | --- | --- | --- | --- |
| **1. Baseline design:** `z_sel` shifts observation only; linear changes in `x_bin x_cont`; logit $p(x)$ and pair-observation logits $q(x, z)$ (the default `selmodel(pair_logit)`) | $p, m, q, j$ all correct | Generally biased | Consistent | Consistent | Consistent | Consistent |
| **2. Nonlinear untreated trend:** untreated changes depend on $x_{\mathrm{cont}}^2$, omitted from `covars(x_bin x_cont)` | $p, q$ correct; linear $m, j$ misspecified | No guarantee | No guarantee | Consistent | No guarantee | Consistent |
| **3. Wrong observation model $q$:** observation depends on $z_{\mathrm{sel}}^2$, omitted from `selcovars(x_bin x_cont z_sel)` | $j$ correct; $p, m$ valid on $W$; fitted $q$ misspecified | No guarantee | No guarantee | No guarantee | Consistent | Consistent if the fitted $q$ has a regular interior limit |
| **4. Mispooled imputer $j$:** set `imppooling(within_cohort)` or `across_cohorts` when treatment effects vary with $x$ | $q$ correct and $m, p$ valid; imposed $j$ pooling misspecified | No guarantee | Consistent | Consistent | No guarantee | Consistent |
| **5. `z_sel` also shifts outcomes:** `z_sel` shifts outcome trends but is omitted from `covars(x_bin x_cont)` while included in `selcovars(x_bin x_cont z_sel)` | $q$ and $p$ correct on $W$; $j$ and $m$ omit outcome-shifting `z_sel` | No guarantee | No guarantee from the conditions on $W$ | Consistent through $q$ and $p$ | No guarantee | Consistent through $q$ and $p$ |

“Consistent” means that a sufficient large-sample condition holds under all stated assumptions and regular limits of the fitted models; “No guarantee” means those sufficient conditions no longer establish consistency (it does not prove bias for a specific cell or exclude special cancellation).

## 4. Specifications in Stata

Sections 3, 4 and 5 of the hands-on do-file estimate the models without correction, with observation IPW, with pure imputation and with the augmented correction. See [Treatment adjustment](../methods/treatment.md), [Observation models](observation-models.md) and [Change imputation](change-imputation.md) for the commands. The scenarios above are reasoning about sufficient conditions, not results of a simulation. The qualifications above on shared fits and on models that cannot be estimated apply.

---

**Reading path, step 16 of 22** · Previous: [GMM weighting and covariance regularization](../methods/gmm.md) · Next: [Clustered inference, confidence bands, and pretrend tests](inference.md)

[Guide home](../README.md) · [Package overview](../../README.md)
