# Cross-pair chaining, control parallel trends and within-unit weighting

If a target cohort and its controls are observed over different intervals, they may have no pair of dates in common, and the comparisons that `cdid` uses by default cannot identify the effect. Cross-pair chaining uses the changes of each cohort separately and links them through calendar dates where the trend assumptions apply. Specify `crosspair` to allow these comparisons, and use `cdid_chain` to see how they identify an effect.

## What is estimated

For cohort $k$ and a pair of dates $s < t$ at which some units of cohort $k$ are observed, the link moment is the weighted mean outcome change of those units,

$$
\hat m_{k,(s,t)} = \frac{\sum_i 1\{G_i = k\} O_{is} O_{it}\, \omega_i\, (Y_{it} - Y_{is})}{\sum_i 1\{G_i = k\} O_{is} O_{it}\, \omega_i},
$$

where $O_{it}$ indicates an observed outcome and $\omega_i$ the sampling weight. For each target cohort $g$, the links of $g$ itself are changes of its level path $\gamma_g(t)$, and the links of each clean control cohort are changes of the untreated path $\lambda_g(t)$ that the target would have followed. Both paths are anchored at the reference date $b_g$, the last date before $g - \delta$ (with $\delta$ the anticipation), where $\gamma_g(b_g) = \lambda_g(b_g) = 0$. The parameters are estimated jointly by GMM from all links, and

$$
\text{ATT}(g,t) = \gamma_g(t) - \lambda_g(t), \qquad t \geq g.
$$

Clean controls are the never-treated units and, with `controlgroup(notyettreated)`, the cohorts $k$ not yet treated at the later date $t$ of the pair, allowing for anticipation ($t + \delta < k$). Pretreatment cells compare the target path with a control path on the same dates and are reported as pretrend contrasts. When several links of one cohort share the same set of observed units, the moments of a cycle of dates are linearly dependent; the estimator keeps a spanning forest of each such set and carries the dropped links in the GMM metric.

## When it helps

### Contiguous short windows need only PT-post

In staggered designs where each cohort is observed for a few dates before adoption and the clean untreated dates of each control cohort form a contiguous interval, every step between a cohort's reference date and a later date is usually covered by some cohort that is still untreated at both dates. Same-pair chained estimation with `controlgroup(notyettreated)` and `ptassumption(post)` then already identifies the post-treatment effects, including long horizons that a two-period comparison cannot reach. `crosspair` is not needed for these designs; on complete panels both approaches use the same information, and `crosspair` mainly changes how the comparisons are weighted.

### Non-overlapping observation pairs

Suppose the target cohort 3 is observed at dates {1, 2, 4}, control cohort 4 at {2, 3} and control cohort 5 at {3, 4}. No control unit has an observed change over (2, 4), and no target unit is observed at date 3, so same-pair comparisons cannot identify ATT(3, 4). With `crosspair`, the effect is identified by

$$
\mathrm{ATT}(3, 4) = \mathbb{E}[Y_4 - Y_2 \mid G = 3] - \mathbb{E}[Y_3 - Y_2 \mid G = 4] - \mathbb{E}[Y_4 - Y_3 \mid G = 5],
$$

under parallel trends between the target and the two control cohorts over the dates involved. `cdid_chain` reports such a combination as the witness path `2->4 [g=3] - 2->3 [ctrl{4}] - 3->4 [ctrl{5}]`.

## Parallel trends choices

For a target cohort $g$ with reference date $b_g$:

- **`ptassumption(post)`** (default) assumes conditional parallel trends between cohort $g$ and the clean control cohorts from $b_g$ onward. Differences before $b_g$ are left free for each control pool.
- **`ptassumption(ctrl)`** keeps PT-post and assumes that the clean control cohorts share identical untreated trends from `ptstart()` onward (by default the first study date). The target cohort's own trend before $b_g$ remains unrestricted, so its pretreatment cells remain checks of its pretrends against the pooled controls.
- **`ptassumption(all)`** requires `ptstart()` and restricts the control cohorts and the target cohort to common untreated trends from `ptstart()` onward.

Across all target cohorts, not-yet-treated PT-post for the earliest treated cohort already implies equal untreated trends among the clean cohorts from that cohort's reference date onward, so the additional restriction of PT-ctrl concerns earlier dates. If a control cohort has an idiosyncratic trend shift before $b_g$, PT-ctrl transmits that bias into post-treatment estimates, while PT-post leaves the shift unrestricted.

## Options

| Option | Effect with `crosspair` |
| --- | --- |
| `method(ipw\|reg\|dr)` | Covariates adjust each control link toward the target cohort: a treatment logit of the target cohort against the control units of the link (`ipw`), a regression of the control changes on the covariates averaged over the target cohort (`reg`), or both (`dr`). Without covariates the three coincide. `method(custom)` is not available. |
| `outcomepooling(separate\|pooled)` | Under `method(reg\|dr)`, a separate regression for each control cohort and pair of dates, or slopes common to all of them with separate intercepts. |
| Covariates and weights that vary over time | Every comparison of target cohort $g$, including its control links, reads covariates at $g$'s reference date. Weights are read at the first study date (`fixweights(first_period)`) or at the reference date (`fixweights(base_period)`). |
| `selection(conditional_mtar)` | `selmethod(ipw)` weights each observed change by the inverse of its fitted observation probability, with every `selmodel()`, `selpooling()` and `seluniverse()`; `selmethod(imputation)` and `selmethod(dr)` predict the changes of every unit of the cohort and require `seluniverse(F)`. The estimation of the observation and imputation models enters the influence functions. |
| `gmmweighting(optimal\|within\|identity)` | Inverse of the estimated covariance of the link moments, inverse of the analytical covariance of within-unit changes under uncorrelated level shocks, or the identity matrix. |
| `noconstant`, factor-variable covariates, `frame()`, `baseperiod()`, `anticipation()`, `controlgroup()`, `cluster()` | As without `crosspair`. |

`crosspair` requires chained mode and is not available with `method(custom)`, `preset()` or `dvar()`. Outcome-regression designs that are rank deficient on a control cohort without an exact alias that also holds on the target cohort are refused, because the counterfactual mean would not be identified without extrapolation.

## The stacked cohort fixed-effects benchmark

Balla-Elliott and Norwich proposed linking successive control cohorts for cohorts with short pretreatment windows; their stacked cohort fixed-effects regression reproduces the calendar trend by linking those windows ([Cohort-Chained DiD](https://pdfs.dballaelliott.com/ccdid.pdf); see also [Estimate familiar DiD models with cdid](../from-other-packages.md#cohort-chained-did) and [Short pretreatment windows](../examples/short-pretreatment-windows.md)). With `controlgroup(notyettreated)`, `baseperiod(universal)`, `crosspair`, `ptassumption(ctrl)` without `ptstart()`, `gmmweighting(within)`, no covariates, no correction for missing outcomes, the same unit weights and the same target and control sample for each slice, the point estimates equal the unit fixed-effects regression coefficients when observation windows are nested within each cohort, and the cohort fixed-effects coefficients when all units of a cohort share one window. The Stata test `tests/stata/test_crosspair_stacked.do` compares these estimates with stacked regressions fitted with `areg` (unit fixed effects) and `regress` (cohort fixed effects) on short-window panels without never-treated units, with sampling weights, monotone attrition within cohorts (where only the unit fixed-effects coefficients coincide), dates missing from the calendar and one period of anticipation; the largest difference is $2.5 \times 10^{-13}$. Inference in `cdid` uses its unit influence functions, not the regression standard errors.

## Inspect the identification chains

`cdid_chain` shows which effects need chaining, a combination of comparisons that identifies each effect, and the GMM weights of all comparisons. For every available cell, the contributions (weight times comparison estimate) sum to the reported estimate up to rounding.

The example builds a panel in which each cohort is observed only during a short window before its treatment date, with no never-treated units. The last cohort has no not-yet-treated controls, so its cells are unavailable, and `narm` restricts the event-study aggregation to the available cells.

```stata
clear
set seed 20261111
set obs 240
generate long id = _n
generate int G = 2005 + mod(_n, 6)
generate double alpha = rnormal(0.7 * (G - 2005), 1.1)
expand 10
bysort id: generate int period = 2000 + _n
keep if period >= max(2001, G - 4)
generate double y = alpha + 0.35 * (period - 2001) + cond(period >= G, 2 + 0.5 * (period - G), 0) + rnormal(0, 0.45)

cdid y, ivar(id) time(period) gvar(G) controlgroup(notyettreated) baseperiod(universal) crosspair ptassumption(ctrl) gmmweighting(within) bstrap(0) cband(0)
cdid_chain, cohort(2006) time(2007)
estat dynamic, narm
```

The chain of cohort 2006 at 2007:

```text
Identification chains (cross-pair, ptassumption(ctrl), gmmweighting(within))
Comparisons used by GMM: 84    Overidentifying restrictions: 20
Cells: 1 (1 estimable, 0 chained without a direct comparison, 0 unavailable)

   Group     Time  Role      Chain               ATT   Std. err.
---------------------------------------------------------------
    2006     2007  post      direct       2.45420282  .118657047
           witness: 2005->2006 [g=2006] + 2006->2007 [g=2006] - 2005->2006 [ctrl{2007}] - 2006->2007 [ctrl{2008}]
```

The event study over the available cells:

```text
Chained difference-in-differences            Number of units       =       240
Method: doubly robust (dr)                   Number of event times =        15
Cross-pair comparisons; GMM weighting: within
Control group: not yet treated; base period: universal
Aggregation: dynamic
Event-time selection: narm
------------------------------------------------------------------------------
Event time                  ATT  Std. err.     [95% conf. interval]
------------------------------------------------------------------------------
Overall                2.854236   .0753267     2.706598    3.001874
------------------------------------------------------------------------------
                  -9          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2010 at date 2001
      (pre-treatment)
                  -8          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2009 at date 2001
      (pre-treatment); cohort 2010 at date 2002 (pre-treatment)
                  -7          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2008 at date 2001
      (pre-treatment); cohort 2009 at date 2002 (pre-treatment); cohort 2010 at date 2003
      (pre-treatment)
                  -6          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2007 at date 2001
      (pre-treatment); cohort 2008 at date 2002 (pre-treatment); cohort 2009 at date 2003
      (pre-treatment); cohort 2010 at date 2004 (pre-treatment)
                  -5          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2006 at date 2001
      (pre-treatment); cohort 2007 at date 2002 (pre-treatment); cohort 2008 at date 2003
      (pre-treatment); cohort 2009 at date 2004 (pre-treatment); cohort 2010 at date 2005
      (pre-treatment)
                  -4          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2005 at date 2001
      (pre-treatment); cohort 2006 at date 2002 (pre-treatment); cohort 2007 at date 2003
      (pre-treatment); cohort 2008 at date 2004 (pre-treatment); cohort 2009 at date 2005
      (pre-treatment); cohort 2010 at date 2006 (pre-treatment)
                  -3  .04145277  .06193165   -.07993103    .1628366
                  -2 -.02531473  .05681671    -.1366734   .08604398
                  -1          0          0            0           0
                   0   1.931259  .05383278     1.825749     2.03677
                   1   2.443667  .07053394     2.305423    2.581911
                   2   2.883516  .08547695     2.715984    3.051048
                   3   3.356337   .1003161     3.159722    3.552953
                   4   3.656401   .1434481     3.375247    3.937554
                   5          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2005 at date 2010
      (post-treatment)
------------------------------------------------------------------------------
Std. err.: analytic; pointwise 95% intervals.
The overall summary uses a pointwise interval.
Pretrend Wald test: not computed; some pre-treatment cells are unavailable, and the test is not
      computed on a subset of them.
```

## References

- Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The chained difference-in-differences. *Journal of Econometrics*, 248, 105783.
- Balla-Elliott, D. and Norwich, A. Cohort-Chained DiD. Working paper. <https://pdfs.dballaelliott.com/ccdid.pdf>

---

**Reading path, step 21 of 22** · Previous: [Diagnose errors and unavailable results](../help/troubleshooting.md) · Next: [Short pretreatment windows and cross-pair chaining](../examples/short-pretreatment-windows.md)

[Guide home](../README.md) · [Package overview](../../README.md)
