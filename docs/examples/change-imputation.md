# Impute actual changes for missing outcomes

Under `selection(conditional_mtar)`, observation weighting (`selmethod(ipw)`) reweights the changes of units observed at both endpoints of a date pair. Change imputation (`selmethod(imputation)`) instead predicts the actual outcome change of every unit in the cohort from its covariates, and the augmented doubly robust correction (`selmethod(dr)`) combines prediction with inverse observation-probability weighting ([Handling missing outcomes](../methods/missing-outcomes.md)).

This page distinguishes the imputation of actual changes from the regression of counterfactual changes and illustrates the pooling choices. Sections 5a and 5b of the hands-on do-file apply pure imputation and the augmented correction.

## 1. Two outcome regressions: $j$ versus $m$

Change imputation fits a linear model $j$ of **actual cohort changes**. It is distinct from the counterfactual control regression $m$ of `method(reg)` and `method(dr)`:

| Feature | Actual-change model $j_c(X)$ | Counterfactual control regression $m_b(X)$ |
| --- | --- | --- |
| Layer & option | Handling missing outcomes: `selmethod(imputation)` or `selmethod(dr)` | Treatment adjustment: `method(reg)` or `method(dr)` |
| Predicts | The **actual** change $\Delta_e Y = Y_t - Y_s$ of cohort $k$ over pair $e$, **treatment effects included** | The **untreated** change $\Delta_e Y(0)$ |
| Fitted to | Units of cohort $k$ observed at both endpoints of pair $e$ ($S_{i,e} = 1$) | Reconstructed control changes on deduplicated comparison blocks $b = (e, K)$ |
| Slopes | Shared according to `imppooling(separate\|within_cohort\|across_cohorts)` | Controlled by `outcomepooling(separate\|pooled)` (default `separate` per block $b$; `pooled` shares one $\beta$ across blocks) |

Both $j$ and $m$ read the treatment covariate list (`indepvars` or `covars(varlist)`, and `framecovars(varlist)`). By contrast, `selcovars(varlist)` and `frameselcovars(varlist)` define only the observation-probability logit $q$.

## 2. How change imputation works (`selmethod(imputation)`)

For each needed cohort $k$ and ordered date pair $e = (s, t)$, `cdid` fits $j_{k,e}(X_i)$ by weighted least squares to the units of cohort $k$ observed at both $s$ and $t$ ($S_{i,e} = 1$):

$$j_{k,e}(X_i) = \alpha_{k,e} + L(e)\, z(X_i)'\beta_r$$

where $L(e) = i_t - i_s$ counts calendar steps in the stored calendar, $\alpha_{k,e}$ is a local cell intercept (unless `noconstant` is specified), and $\beta_r$ is the slope vector of imputation pooling group $r$.

The fitted model $\hat{j}_{k,e}(X_i)$ is then evaluated for **every unit of cohort $k$ in the study population $F$** (`frame()`), including units whose change was observed, singletons observed once ($B \setminus A$), and unit-table-only units observed never ($F \setminus B$). The average of these predictions over $F$ replaces the observed cohort mean:

- You do not fill in `y` in your dataset; `y` remains `.` where missing.
- A singleton unit observed at only one date supplies no training pair ($S_{i,e} = 0$), and the method does not condition on its single observed outcome level to impute the missing endpoint; it evaluates the predicted change $\hat{j}_{k,e}(X_i)$ from $X_i$.
- Under pure imputation (`selmethod(imputation)`), `selmodel()`, `selpooling()`, and `selcovars()` do not apply (`e(inactive_options)` lists them, together with `outcomepooling()` under `method(ipw)`), and no observation logit $q$ is fitted.

## 3. Three kinds of slope sharing (`imppooling()`)

Because $j_{k,e}(X_i)$ predicts **actual** changes, including treatment effects for treated cohorts when $t \ge g$, it must represent how treatment effects vary with $X_i$:

- Suppose untreated outcome changes have slope $\beta_0\, L(e)$ on $X_i$, while treatment adds an extra effect $\delta_X\, X_i$ at post-treatment dates $t \ge g$.
- For any pair $e = (s, t)$ of treated cohort $g$ that **crosses treatment onset** ($s < g \le t$), the actual change $\Delta_e Y_i$ has total slope $\beta_0\, L(e) + \delta_X$ on $X_i$, so the duration-scaled slope $\beta_r$ in $L(e)\, X_i'\beta_r$ equals $\beta_0 + \delta_X / L(e)$.
- For pre-treatment pairs ($s < t < g$), post-onset pairs where the effect on $X_i$ is constant after onset ($g \le s < t$), and control cohort pairs ($k = 0$), the duration-scaled slope is $\beta_0$.

This makes the three `imppooling()` options distinct restrictions:

| `imppooling()` | Slopes of $j$ | Implication when treatment effects vary with $X_i$ ($\delta_X \ne 0$) |
| --- | --- | --- |
| `separate` (default) | Separate slope vector $\beta_{k,e}$ for each cohort and pair $(k, e)$ | Correctly accommodates pair-specific slope shifts $\beta_0 + \delta_X / L(e)$ on treatment-crossing pairs. |
| `within_cohort` | One shared slope vector $\beta_k$ per cohort across all of its pairs | Forces pre-treatment pairs, treatment-crossing pairs, and post-treatment pairs of cohort $g$ to share one duration-scaled slope $\beta_g$, which is misspecified when $\delta_X \ne 0$. |
| `across_cohorts` | One global slope vector $\beta$ across all cohorts and pairs | Forces treated and control cohorts to share one duration-scaled slope $\beta$, which is misspecified when $\delta_X \ne 0$. |

Each cohort and pair keeps its local intercept $\alpha_{k,e}$ whenever `noconstant` is not specified.

## 4. The augmented doubly robust correction (`selmethod(dr)`)

When `selmethod(dr)` is specified alongside built-in `method(ipw|reg|dr)` on `seluniverse(F)` (with `selrho()` omitted), `cdid` fits **both** $j_c(X_i)$ and $q_c(W_i)$ and reconstructs for every unit $i$ of cohort $k$ in $F$:

$$Z_{i,c} = \hat{j}_c(X_i) + S_{i,e}\,\frac{\Delta_{i,e} - \hat{j}_c(X_i)}{\max(\hat{q}_c(W_i),\, 10^{-6})}$$

and then runs the chosen built-in treatment adjustment (`method(ipw)`, `method(reg)`, or `method(dr)`) on $Z_{i,c}$.

| Correction | Missing-outcome layer requirement | Treatment layer requirement (`method(ipw\|reg\|dr)`) |
| --- | --- | --- |
| `selmethod(ipw)` | $q_c$ correct in every needed cell | `ipw`: $p$ correct; `reg`: $m$ correct; `dr`: $p$ **or** $m$ correct |
| `selmethod(imputation)` | $j_c$ correct in every needed cell/pooling group | `ipw`: $p$ correct; `reg`: $m$ correct; `dr`: $p$ **or** $m$ correct |
| `selmethod(dr)` | For each needed cell $c$, $q_c$ correct **or** $j_c$ correct | `ipw`: $p$ correct; `reg`: $m$ correct; `dr`: $p$ **or** $m$ correct |

Note that under pure imputation (`selmethod(imputation)`), `method(dr)` does **not** repair a misspecified imputer $j_c$: if $j_c$ is misspecified, the reconstructed cohort changes fed into treatment adjustment are already biased unless repaired by the observation-weighted residual term in `selmethod(dr)` ([Compare model robustness](../examples/robustness.md)).

## 5. Examples in Stata

Section 5 of the hands-on do-file runs pure imputation and the augmented correction, with separate pooling:

```stata
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) selection(conditional_mtar) ///
    selmethod(imputation) imppooling(separate) ///
    biters($CDID_DEMO_BITERS)
display as text "Inactive options under selmethod(imputation): " as result "`e(inactive_options)'"

* 5b. Doubly robust observation + treatment adjustment (`selmethod(dr)` with `method(dr)`)
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) selection(conditional_mtar) ///
    selmethod(dr) selmodel(pair_logit) selpooling(separate) ///
    selcovars(i.x_bin c.x_cont c.z_sel) imppooling(separate) ///
    biters($CDID_DEMO_BITERS)
```

To compare other imputation restrictions, change `imppooling(separate)` to `within_cohort` or `across_cohorts` while holding the remaining specification and target fixed. The restrictions can be false when covariates modify treatment effects. A small SE or a successful fit does not establish their validity.

---

**Reading path, step 12 of 22** · Previous: [Choose an observation model for missing outcomes](observation-models.md) · Next: [Parallel trends](../methods/parallel-trends.md)

[Guide home](../README.md) · [Package overview](../../README.md)
