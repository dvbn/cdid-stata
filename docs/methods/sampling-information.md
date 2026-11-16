# Sampling information for the observation model

Observation weighting under `selection(conditional_mtar) selmethod(ipw)` fits an observation model on a set of units chosen with `seluniverse()`. This page explains the choice of that set and the use of known transport weights. Changing the set changes the probability being modeled, not the target population, and not only the number of observations in the likelihood.

## 1. The three fitting universes

The nested sets are $A \subseteq B \subseteq F$: $F$ is the study population, $B$ has at least one observed outcome and $A$ has at least two. A unit observed on a particular pair belongs to all three. Singletons and units never observed remain in the target population even when the observation logit is fitted on a subset.

Only `selmethod(ipw) selmodel(pair_logit)` offers `seluniverse(B)` and `seluniverse(A)`. Independent-period logits, pure imputation and the augmented correction require `F`; imputation and the augmented correction also do not allow `selrho()`.

## 2. What a restricted fitting set assumes

Let $W$ be a set of covariates sufficient for the assumption on missing outcomes and let $U$ be $B$ or $A$. Define

$$q_U(e,W,G)=P(S_e=1\mid W,G,U),\qquad \pi_U(G)=P(U\mid G).$$

If no `selrho()` column is supplied, the required `seljustification()` states the assumption that $U$ and the cohort have the same covariate distribution in the weighted population. With known transport, the column is the fixed, known, normalized density ratio

$$\rho_U(W,G)=\frac{f_P(W\mid G,U)}{f_P(W\mid G)},$$

with the **fitting-set density in the numerator**. Since $S_e=1$ implies membership in $U$,

$$P(S_e=1\mid W,G)=\pi_U(G)\,\rho_U(W,G)\,q_U(e,W,G).$$

Within a cohort, the observation adjustment is proportional to the inverse of this probability. Up to cohort-rate normalizations, its unit factor is

$$w_i S_{i,e}\big/\{\rho_U(W_i,G_i)\max(\hat q_U(e,W_i,G_i),10^{-6})\}.$$

The ratio divides the weights. `cdid` neither estimates nor renormalizes it. It must be positive and finite on the fitting set; estimated or unnormalized ratios are not supported. A valid ratio also requires the density of the fitting set to be positive wherever the target population needs it. A logit that fits well establishes none of these assumptions.

For example, let a cohort have population shares $(0.5,0.5)$ at binary $X=(0,1)$, fitting-set mass $P(A\mid G)=0.4$, fitting-set shares $(0.25,0.75)$, ratios $(0.5,1.5)$ and $q_A=0.5$ in both strata. The population pair probabilities are $(0.1,0.3)$. Their inverse factors $(10,10/3)$ restore equal population shares. Assuming a ratio of one would not.

## 3. Supported settings

| Setting | Observation model | Requirement |
| --- | --- | --- |
| `seluniverse(F)` | `pair_logit` or `independent_period_logit` | Omit `selrho()`; justification is optional. |
| `seluniverse(B)` or `seluniverse(A)` | `selmethod(ipw) selmodel(pair_logit)` only | Nonempty `seljustification()`; without `selrho()`, the covariate distribution of the set is assumed representative. |
| Restricted universe with `selrho(varname)` | `selmethod(ipw) selmodel(pair_logit)` only | Known positive normalized ratio as defined above, and `seljustification()`. |

With covariates that are constant within units, the observation covariates and the ratio must be nonmissing on the fitting set. The other models used keep their own requirements on the covariates of the whole population. The identity above concerns covariates that are constant within units; it does not extend automatically to covariates that change over time.

## 4. Example

With the data of the hands-on do-file in memory, a pair logit fitted on $B$ is specified as:

```stata
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    selection(conditional_mtar) selmethod(ipw) ///
    selmodel(pair_logit) seluniverse(B) ///
    seljustification("Substantively justified representative fitting set B") ///
    frame(cdid_pop) framecovars(x_bin x_cont) bstrap(0) cband(0)
```

Use `A` instead of `B` only with a justification that applies to $A$. An analysis with known transport also supplies the known ratio through `selrho()`; the simulated data contain no such ratio, and one should not be estimated or invented to run this example. A fit on a restricted set can still fail the checks of support, rank or convergence. Fitting on fewer observations implies neither faster computation nor smaller standard errors.

---

**Reading path, step 10 of 22** · Previous: [Handling missing outcomes](missing-outcomes.md) · Next: [Choose an observation model for missing outcomes](../examples/observation-models.md)

[Guide home](../README.md) · [Package overview](../../README.md)
