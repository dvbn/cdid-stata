# Choose an observation model for missing outcomes

When outcome observation depends on covariates that also shift outcome changes or treatment effects, `selection(conditional_mtar)` with `selmethod(ipw)` (or `selmethod(dr)`) estimates the probability $q_c(W_i)$ that a unit in cohort $k$ is observed at both endpoints of date pair $e = (s, t)$ ([Handling missing outcomes](../methods/missing-outcomes.md)).

This page compares the two observation models and the three pooling restrictions on their slopes. The examples below are Section 4 of the hands-on do-file.

## 1. Two observation model families (`selmodel()`)

For each needed cohort $k$ and ordered date pair $e = (s, t)$, let $R_{i,t} \in \{0, 1\}$ indicate whether outcome $Y_{i,t}$ is observed at date $t$, so pair observation is the product $S_{i,e} = R_{i,s} R_{i,t} \in \{0, 1\}$.

### Direct pair logit (`selmodel(pair_logit)`, default)

`selmodel(pair_logit)` directly models the pair indicator $S_{i,e} = R_{i,s} R_{i,t}$ using a logistic link:

$$q_{k,e}(W_i) = P(S_{i,e} = 1 \mid G_i = k, W_i) = \Lambda\!\left(\alpha_{k,e} + z_Q(W_i)'\gamma_r\right)$$

where $\Lambda(u) = (1 + e^{-u})^{-1}$, $\alpha_{k,e}$ is a cell-specific intercept for cohort $k$ and pair $e$ (unless `selnoconstant` or `selzero` is specified), and $\gamma_r$ is the slope vector for pooling group $r$.

- **Advantage:** Does not require $R_{i,s}$ and $R_{i,t}$ to be conditionally independent given $(G_i = k, W_i)$; it models their joint pair indicator directly.
- **When it is exact vs approximate:** If single-date observation $R_{i,t}$ follows independent date logits $\Lambda(\theta_{k,t} + W_i'\delta)$, their product $\Lambda(\theta_{k,s} + W_i'\delta)\,\Lambda(\theta_{k,t} + W_i'\delta)$ is not itself a single logistic function when $W_i$ is continuous (though with a single binary covariate and `selpooling(separate)`, a saturated two-parameter pair logit still represents every two-point probability exactly).

### Product of date logits (`selmodel(independent_period_logit)`)

`selmodel(independent_period_logit)` fits a binary logit for single-date observation $R_{i,d}$ at each needed cohort $k$ and date $d \in \{s, t\}$:

$$r_{k,d}(W_i) = P(R_{i,d} = 1 \mid G_i = k, W_i) = \Lambda\!\left(\alpha_{k,d} + z_Q(W_i)'\gamma_r\right)$$

and constructs the pair observation probability as the product of the two endpoint probabilities:

$$q_{k,(s,t)}(W_i) = r_{k,s}(W_i) \, r_{k,t}(W_i)$$

- **Advantage:** Directly matches wave-by-wave independent response designs and fits one logit per needed cohort-date $(k, d)$ rather than one per cohort-pair used.
- **Requirement:** Requires that $R_{i,s}$ and $R_{i,t}$ are conditionally independent given $(G_i = k, W_i)$ on the fitting universe and requires `seluniverse(F)`. A mean-outcome restriction alone does not justify replacing the joint response probability by a product.

## 2. Three slope-pooling choices (`selpooling()`)

Across both `pair_logit` and `independent_period_logit`, every cell (each cohort-pair $(k, e)$ or cohort-date $(k, d)$) keeps its **own local intercept** whenever an intercept is included. The option `selpooling()` governs how the **covariate slope vector** $\gamma_r$ is shared across cells:

| `selpooling()` | Slope restriction on $q$ | Fit structure | When appropriate |
| --- | --- | --- | --- |
| `separate` (default) | Separate slope vector $\gamma_{k,e}$ per cohort-pair (or $\gamma_{k,d}$ per cohort-date) | One fit per cell used | Covariate effects on observation can vary across cohorts and across calendar dates/pairs. |
| `within_cohort` | One shared slope vector $\gamma_k$ per cohort across all of its pairs or dates (with separate intercepts per pair/date) | One joint fit per cohort | Covariate effects on observation are stable over time within each cohort, but may differ between cohorts. |
| `across_cohorts` | One global slope vector $\gamma$ shared across all cohorts and all pairs or dates (with separate intercepts per cell) | One joint fit across all cells used | Covariates shift log-odds of observation identically across all cohorts and periods, while baseline response rates shift via cell intercepts. |

> [!WARNING]
> **Slope pooling does not remove cell intercepts or cure complete separation:** Because `within_cohort` and `across_cohorts` keep a free intercept $\alpha_c$ for every cell $c$, a cell in which 100% or 0% of units are observed still has no finite maximum likelihood estimate of its intercept, and estimation stops with `r(498)`. `cdid` does not switch to another pooling ([Diagnose errors and unavailable results](../help/troubleshooting.md)).

## 3. Examples in Stata

Section 4a compares the two models with a binary covariate. Section 4b uses separate covariates for the outcome and observation models and the unit table of the study population:

```stata
cdid y i.x_bin, ivar(id) time(period) gvar(first_treat) ///
    selection(conditional_mtar) selmodel(pair_logit) selpooling(separate) ///
    biters($CDID_DEMO_BITERS)

cdid y i.x_bin, ivar(id) time(period) gvar(first_treat) ///
    selection(conditional_mtar) selmodel(independent_period_logit) selpooling(separate) ///
    biters($CDID_DEMO_BITERS)

* 4b. Separate outcome covariates `covars()` and selection covariates `selcovars()`
*     with the unit table `frame(cdid_pop)` of the study population
cdid y, ivar(id) time(period) gvar(first_treat) ///
    covars(i.x_bin c.x_cont) ///
    selection(conditional_mtar) selmodel(pair_logit) selpooling(separate) ///
    selcovars(i.x_bin c.x_cont c.z_sel) ///
    frame(cdid_pop) frameivar(id) framegvar(first_treat) ///
    framecovars(x_bin x_cont) frameselcovars(x_bin x_cont z_sel) ///
    biters($CDID_DEMO_BITERS)
```

The first two models omit `z_sel` on purpose; the simulated data also include units observed once or never. That these models can be estimated says nothing about whether they are correct or give the same estimates. Section 4c pools the pair logits within and across cohorts. Pooling keeps each cell intercept, so a perfectly separated cell still has no finite estimate.

---

**Reading path, step 11 of 22** · Previous: [Sampling information for the observation model](../methods/sampling-information.md) · Next: [Impute actual changes for missing outcomes](change-imputation.md)

[Guide home](../README.md) · [Package overview](../../README.md)
