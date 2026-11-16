# Parallel trends

In chained difference-in-differences (`cdid`, Bellégo, Benatia and Dortet-Bernadet 2025; Callaway and Sant'Anna 2021), identifying treatment effects requires assuming that, in the absence of treatment, treated and control cohorts would have followed parallel average outcome trajectories conditional on covariates $X_i$.

In Stata, four options govern how parallel trends and pre-treatment comparisons are defined:

- `ptassumption(post|ctrl|all)` (default `post`): selects the parallel-trends restrictions imposed in the GMM system; `ctrl` adds common trends among the eligible control cohorts, and `all` common trends of the target and control cohorts, from `ptstart()`,
- `ptstart(#)`: specifies the common calendar start date, required with `ptassumption(all)`, optional with `ptassumption(ctrl)` (where it defaults to the first study date), and not allowed with `ptassumption(post)`,
- `anticipation(#)` (default `0`): sets the number of pre-treatment calendar periods in which units may anticipate treatment, fixing each cohort's reference period $r(g)$ at the latest study date $d$ with $d < g - \text{anticipation}$ ($g - 1 - \text{anticipation}$ when the dates are consecutive integers),
- `baseperiod(varying|universal)` (default `varying`): independently controls how unrestricted pre-treatment comparisons (`role == "pre"`) are parameterized and reported.

## 1. Post-treatment parallel trends (`ptassumption(post)`, default)

Under `ptassumption(post)` (`PT-post`), conditional parallel trends are assumed **from each cohort's reference period $r(g)$ onward**. Specifically, for every post-reference date pair $(s, t)$ with $r(g) \le s < t$ and control pool $K$:

$$E[Y_t(0) - Y_s(0) \mid X, G = g] = E[Y_t(0) - Y_s(0) \mid X, G \in K]$$

What about date pairs touching periods **before** $r(g)$?

- Under `ptassumption(post)`, pre-reference untreated differentials $\delta_{g, d \to d+1}$ for $d < r(g)$ are **unrestricted** free parameters in the GMM system.
- Under the other assumptions on the comparisons, on observation and on the adjustment models, earlier differentials can enter the model without being restricted to zero. This is a statement about identification in the population, not about exact invariance in finite samples.

### Comparisons through earlier dates under `ptassumption(post)`

Even though pre-reference differentials are unrestricted under `ptassumption(post)`, two-date comparisons starting before $r(g)$ are still used whenever their unrestricted differentials before the reference period cancel through an intermediate date:

- Suppose cohort $g = 3$ has reference period $r(3) = 2$, and units are observed on pairs $(1, 2)$ and $(1, 3)$.
- In the GMM design matrix $A$, the comparison on $(1, 2)$ loads on $\delta_{3, 1 \to 2}$, while the comparison on $(1, 3)$ loads on $\delta_{3, 1 \to 2} + \text{ATT}(3, 3)$.
- Taking the difference between the two moments cancels $\delta_{3, 1 \to 2}$:
  $$E[Y_3 - Y_1] - E[Y_2 - Y_1] = E[Y_3 - Y_2]$$
  and identifies $\text{ATT}(3, 3)$ without assuming parallel trends between $t = 1$ and $t = 2$.

## 2. Adding pre-reference restrictions (`ptassumption(all) ptstart(#)`)

The Stata snippets on this page use the simulated panel and the global `CDID_DEMO_BITERS` created by Sections 0 and 1 of the [hands-on do-file](../hands-on.md); run those sections first.

When domain knowledge justifies parallel trends starting at an earlier calendar date $t_{\text{start}}$ common to all cohorts, specify `ptassumption(all)` together with `ptstart(#)`:

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) controlgroup(nevertreated) ptassumption(all) ptstart(1) ///
    biters($CDID_DEMO_BITERS)
```

What `ptassumption(all) ptstart(#)` changes, and what it keeps unchanged:

1. **Keeps every PT-post condition and each cohort's reference period $r(g)$:** Post-treatment effects $\text{ATT}(g, t)$ ($t \ge g$) are still defined relative to each cohort's reference period $r(g)$ ($t = 2$ for $g = 3$, $t = 4$ for $g = 5$).
2. **Constrains the differentials from `ptstart(#)` up to $r(g)$ to zero:** Every differential on an untreated step between `ptstart(#)` and $r(g)$ is constrained to zero in the GMM design matrix $A$. As a result, a comparison from $t = 1$ to $t = 3$ for cohort $g = 3$ directly loads on $\text{ATT}(3, 3)$ without needing to subtract a separate $(1, 2)$ comparison.
3. **Remains a chained within-unit GMM estimator:** `ptassumption(all)` adds zero restrictions on the differentials to the chained GMM system, which uses the same within-unit comparisons. It does **not** switch to a cross-sectional estimator, and it does **not** guarantee smaller standard errors in every cell or in every finite sample draw ([How to compare precision](../examples/precision.md), [Add PT-all restrictions](../examples/pt-all-comparison.md)).
4. **Validation:**
   - `ptstart(#)` is **required** with `ptassumption(all)` and must be one of the periods of the data.
   - `ptstart(#)` is **not allowed** with `ptassumption(post)`.

### Choosing `ptstart(#)`

Always choose `ptassumption(all)` and `ptstart(#)` substantively from institutional or domain knowledge **before** inspecting point estimates, standard errors, or pretrend test $p$-values. Data-mining `ptstart(#)` after viewing pretrend tests invalidates subsequent inference.

## 3. Base-period parameterization (`baseperiod(varying|universal)`)

Separately from `ptassumption(post|all)`, the option `baseperiod(varying|universal)` controls how unrestricted pre-treatment comparisons (`role == "pre"`, where $t < g - \text{anticipation}$) are defined and reported in `e(b)`:

- **`baseperiod(varying)` (default):** Reports consecutive pre-treatment step comparisons (comparing each pre-treatment period $t$ to the immediately preceding period in the stored calendar).
- **`baseperiod(universal)`:** Reports every non-reference pre-treatment period $t$ relative to the cohort's fixed reference period $r(g)$.

Crucially, **pre-treatment cells (`role == "pre"`) are never forced to zero in `e(b)` by `ptassumption(all)`**: under both `ptassumption(post)` and `ptassumption(all)`, the reported pre-treatment cells remain direct unrestricted pre-period diagnostics.

## 4. How the row-space condition determines cell availability

The availability of each post-treatment effect is checked against the row space of the design. Pretreatment cells, which are direct comparisons, have their own support rules. Unsupported cells stay in the results with missing estimates and inference and a stated reason; `e(b)` holds only the available estimates, mapped to cells by `e(coef_map)`. Use `cdid_diagnostics, detail` or `cdid_export, type(cells)` to see every cell. Regularizing the covariance cannot create support.

## 5. Pretrend tests under each choice

- **Unconditional joint Wald pretrend test (`e(wald_statistic)`, `e(wald_df)`, `e(wald_p)`):** Jointly tests that all unrestricted pretreatment cells (`role == "pre"`) equal zero. Because `role == "pre"` cells remain unrestricted under `ptassumption(all)`, the Wald test evaluates the exact same unrestricted pre-period null under both `ptassumption(post)` and `ptassumption(all)`.
- **Standalone conditional pretrend test (`cdid_pretest`):** Tests that conditional pre-treatment mean comparisons equal zero across the covariate distribution on $F$ (`r(CvM)`, `r(CvMcval)`, `r(CvMpval)`). `cdid_pretest` evaluates direct pre-treatment pairs and therefore has no `ptassumption()` or `ptstart()` options ([Inference](../examples/inference.md)).

## 6. Example

Section 3c of the hands-on do-file estimates the model under PT-all from period 1:

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) controlgroup(nevertreated) ptassumption(all) ptstart(1) ///
    biters($CDID_DEMO_BITERS)
```

`ptstart()` names one of the periods of the data. PT-all leaves the pretreatment cells, which serve as diagnostics, unrestricted. See [Add PT-all restrictions](../examples/pt-all-comparison.md) for an example in which PT-all is false.

---

**Reading path, step 13 of 22** · Previous: [Impute actual changes for missing outcomes](../examples/change-imputation.md) · Next: [Add PT-all restrictions](../examples/pt-all-comparison.md)

[Guide home](../README.md) · [Package overview](../../README.md)
