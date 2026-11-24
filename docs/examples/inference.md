# Clustered inference, confidence bands, and pretrend tests

After estimating group-time average treatment effects on an incomplete panel, applied analyses typically need standard errors that account for within-cluster dependence, simultaneous confidence bands across group-time cells or event-study horizons, and diagnostic pretrend tests. `cdid` propagates complete $N_F$-row unit influence functions through every output without refitting the underlying models.

This page explains the standard errors, intervals and bands reported by `cdid`, clustering, and the two pretrend tests. Sections 6 through 8 of the hands-on do-file apply them.

## 1. Pointwise intervals versus simultaneous bands across output families

In Stata `cdid`, `bstrap(#)` (default `1`) and `cband(#)` (default `1`) control standard errors and confidence bands across group-time cells and `estat` aggregations:

| `bstrap()` | `cband()` | Reported inference for `cdid` and `estat group\|dynamic\|calendar` | Group-time cells of `cdid` | Scalar `estat simple` and `Overall ATT` summaries | Separate bins of `estat group`, `dynamic`, `calendar` |
| --- | --- | --- | --- | --- | --- |
| `0` | `0` | `analytical` | Analytical sandwich SEs, pointwise normal CIs | Analytical sandwich SE, pointwise normal CI | Analytical sandwich SEs, pointwise normal CIs |
| `1` | `0` | `multiplier_bootstrap_se` | Multiplier IQR SEs, pointwise normal CIs | Multiplier IQR SE, pointwise normal CI | Multiplier IQR SEs, pointwise normal CIs |
| `1` | `1` | `multiplier_bootstrap_cband` | Multiplier IQR SEs, simultaneous band (`e(critical_value)`) | Multiplier IQR SE (`multiplier_bootstrap_se`), pointwise normal CI | Multiplier IQR SEs, simultaneous band (`e(critical_value)`) |
| `0` | `1` | `analytical` (on `cdid`) / `analytical_cband` (on multi-bin `estat`) | Analytical sandwich SEs, pointwise normal CIs (`0` draws) | Analytical sandwich SE, pointwise normal CI (`0` draws) | Analytical sandwich SEs with a multiplier simultaneous critical value `e(critical_value)` |

Three distinctions matter:

1. `e(se)` holds the reported standard errors, analytical or from the multiplier bootstrap; `e(V)` is the analytical covariance. `e(ci)` holds the reported intervals or bands, and `e(ci_pointwise)` the pointwise intervals when they are available. In ordinary mode, some pointwise intervals can be missing; missing values are not recomputed on export.
2. Bands over the bins of an aggregation differ from the interval of the overall estimate, which is always a pointwise normal interval. `e(overall_inference_kind)` states the inference used for the overall estimate.
3. Each cell fit and aggregation forms its own simultaneous family; separate commands do not establish one global band.

## 2. How the clustered multiplier bootstrap and influence functions work

- **No resampling or refitting:** `cdid` computes an $N_F \times K$ influence-function matrix with one row for every unit in the study population $F$ (`N_F = e(N_full)`), including singletons ($B \setminus A$) and unit-table-only units ($F \setminus B$). These influence functions account for the estimation of every model used ($p, m, q, j$), of the observation rates, and (in `estat`) of the cohort shares.
- **Clustered row sums and Rademacher draws:** When `cluster(varname)` (and `framecluster(name)`) is specified, influence-function rows are summed within each of the $C = \text{e(N\_clusters)}$ clusters on $F$. Each of the `biters(#)` multiplier bootstrap draws multiplies the $C$ cluster sums (or the $N_F$ unit rows when unclustered) by independent Rademacher signs ($\pm 1$ with probability $0.5$) from Stata's random-number stream (`set seed #` reproduces them exactly).
- **Interquartile range (IQR) scale and symmetric bands:** Multiplier bootstrap standard errors use the bootstrap draw interquartile range divided by the standard normal interquartile range ($\Phi^{-1}(0.75) - \Phi^{-1}(0.25)$), and confidence bands are symmetric: $\hat{\theta}_j \pm c \times \text{SE}_j$.
- **Clustering changes inference:** `gmmweighting(optimal)` bases the weights on the covariance of study-unit influence functions on $F$. `gmmweighting(within)` uses the working covariance of within-unit changes. Tikhonov regularization, when needed, acts on the covariance chosen for weighting. Adding or removing `cluster()` leaves `e(b)` identical and can change standard errors, confidence intervals, and Wald statistics.
- **Point estimates only (`noinffunc`):** `cdid` uses the same weighting and returns the same point estimates. Standard errors, intervals and tests are unavailable (`e(wald_status) == "not_computed"`). `estat` cannot aggregate these results.

## 3. Unconditional Wald pretrend test vs standalone conditional `cdid_pretest`

`cdid` provides two distinct pretrend tests that test different null hypotheses:

| Feature | Unconditional joint Wald test (`cdid`) | Conditional Cramér-von Mises test (`cdid_pretest`) |
| --- | --- | --- |
| Command and results | Computed by `cdid` (results in `e()`) | Separate command `cdid_pretest` (results in `r()`; `e()` is unchanged) |
| Default `method()` | `method(dr)` (default of `cdid`) | `method(ipw)` (default of `cdid_pretest`; pass `method(dr)` or `method(reg)` to override) |
| Parallel trends & GMM options | Uses the `ptassumption()` and `gmmweighting()` of the `cdid` estimation | Evaluates direct pre-treatment pairs; has **no** `ptassumption()`, `ptstart()`, `gmmweighting()`, or `gmmcovariance()` options |
| Null hypothesis | Jointly tests that all pretreatment cells other than reference cells (`role == "pre"`), as defined by `baseperiod()`, equal zero **unconditionally** | Tests that, for every pretreatment cell other than reference cells, as defined by `baseperiod()`, the **conditional mean comparison** given $X$ at the earlier date equals zero over the distribution of the units of $F$ |
| Stored results | `e(wald_statistic)`, `e(wald_df)`, `e(wald_p)`, and `e(wald_status)` (`"computed"`, `"unavailable"`, `"not_applicable"`, `"not_computed"`) | `r(CvM)`, `r(CvMb)` (a B-by-1 draw vector), `r(CvMcval)`, `r(CvMpval)`, and `r(status)` (`"available"`, `"unavailable"`, `"not_applicable"`), `r(reason)` |
| Availability & rank requirements | Requires **every** such pretreatment cell to be available, their joint covariance to be well conditioned, and $\text{df} \le C - 1$ (or $N_F - 1$ without clustering). Never tests a subset. | Requires **every** such pretreatment cell to be identified and `covars()` to be nonmissing for **every** unit of $F$. |

> [!NOTE]
> **Cluster rank bound for the joint Wald test:** Centered influence-function sums across $C$ clusters have covariance rank at most $C - 1$. If the number of pretreatment restrictions exceeds $C - 1$, the joint Wald test is reported as `e(wald_status) = "unavailable"` even if floating-point rounding yields positive Cholesky pivots. Point estimates and cell/aggregate confidence intervals remain available.

## 4. Conditional process interpretation

The conditional test evaluates the covariates of the treatment adjustment at the earlier date of each cell, for every unit of the study population. Covariates used only in the observation model are not used for this evaluation. The residuals follow `method(ipw|reg|dr)` and the correction for missing outcomes; they are not residuals of a regression of outcome levels. If covariate values are missing, or if a pretreatment cell is not identified, the chained test is unavailable. Neither test replaces the joint null hypothesis with a subset.

`cdid_pretest` always uses multiplier draws for its CvM distribution, including when its validated `bstrap()`/`cband()` flags are zero. Its default is IPW. The CvM statistic is unstandardized; KS, covariance and Wald slots are absent. A pretest p-value does not establish parallel trends, and selecting a specification after testing requires its own inference argument.

## 5. Examples

Section 8 executes the fitted-model Wald inspection and the separate conditional test:

```stata
* 8a. The pretrend Wald test stored in `e()`
cdid y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) biters($CDID_DEMO_BITERS)
display as text "Wald test status: " as result "`e(wald_status)'" ///
    as text "  chi2(" as result e(wald_df) as text ") = " as result %9.4f e(wald_statistic) ///
    as text "  p = " as result %6.4f e(wald_p)

* 8b. Conditional pretrend CvM test (`cdid_pretest`)
*     The default treatment adjustment is `ipw`; `e()` from `cdid` is unchanged.
cdid_pretest y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    biters($CDID_DEMO_BITERS)
return list

* 8c. Conditional pretrend test under `selection(conditional_mtar)`
cdid_pretest y i.x_bin c.x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) selection(conditional_mtar) selmodel(pair_logit) selpooling(separate) ///
    selcovars(i.x_bin c.x_cont c.z_sel) biters($CDID_DEMO_BITERS)
```

---

**Reading path, step 17 of 22** · Previous: [Compare model robustness across treatment and missing outcome adjustments](robustness.md) · Next: [Options reference](../reference/options.md)

[Guide home](../README.md) · [Package overview](../../README.md)
