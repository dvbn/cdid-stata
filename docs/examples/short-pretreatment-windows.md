# Short pretreatment windows and cross-pair chaining

Panels with staggered adoption often observe each unit only over a limited window around treatment: administrative registers that follow workers over fixed windows around an event, rotating household surveys, or longitudinal surveys that switch from annual to biennial waves. Without never-treated units, two-period comparisons with not-yet-treated controls run out of control cohorts at long horizons: by the time an early cohort $g$ reaches $g + 4$, the cohorts observed at its reference date $g - 1$ have already been treated, and the cohorts still untreated at $g + 4$ were not yet observed at $g - 1$. This page shows how chained estimation recovers these effects, how to inspect the control chains with `cdid_chain`, and how to choose between `ptassumption(post)` and `ptassumption(ctrl)` and between `gmmweighting(optimal)` and `gmmweighting(within)`.

## A staggered short-window panel without never-treated units

The panel covers 2001 to 2010 and six adoption cohorts, 2005 to 2010, of 100 units each, with no never-treated units. Each cohort $g$ is observed from $\max(2001, g - 4)$ to 2010. The unit effect is related to the cohort, and the effect is $1 + (t - g)$ from the adoption date onward, so that ATT(2005, 2009) is 5.

```stata
clear
set seed 20261111
set obs 600
generate long id = _n
generate int G = 2005 + floor((id - 1) / 100)
generate double alpha = rnormal(0.5 * (G - 2005), 1.2)
generate int first = max(2001, G - 4)
expand 2010 - first + 1
bysort id: generate int period = first + _n - 1
generate double y = alpha + 0.4 * (period - 2001) + cond(period >= G, 1 + (period - G), 0) + rnormal()
```

For cohort 2005, the reference date is 2004. Only cohorts 2009 and 2010 are untreated in 2008, and only cohort 2010 in 2009, but cohort 2009 enters the panel in 2005 and cohort 2010 in 2006: neither is observed in 2004. No single comparison group spans 2004 to 2008 or 2004 to 2009.

## Within-unit weighting and optimal weighting

In this design all units of a cohort share the same dates, and neither fit adjusts for covariates or for missing outcomes. The first fit uses cross-pair links, control parallel trends and within-unit weighting; its point estimates equal the coefficients of stacked target-cohort regressions on the same clean-control samples ([Cohort-Chained DiD](../from-other-packages.md#cohort-chained-did)). The second uses the weaker PT-post with the default optimal weighting.

```stata
cdid y, ivar(id) time(period) gvar(G) controlgroup(notyettreated) baseperiod(universal) crosspair ptassumption(ctrl) gmmweighting(within) bstrap(0) cband(0)
estimates store within
estat dynamic, narm
cdid y, ivar(id) time(period) gvar(G) controlgroup(notyettreated) crosspair ptassumption(post) bstrap(0) cband(0)
cdid y, ivar(id) time(period) gvar(G) controlgroup(notyettreated) bstrap(0) cband(0)
```

The within-weighted fit:

```text
Chained difference-in-differences            Number of units       =       600
Method: doubly robust (dr)                   Number of cells       =        60
Cross-pair comparisons; GMM weighting: within
Control group: not yet treated; base period: universal
------------------------------------------------------------------------------
   Group     Time  Role             ATT  Std. err.     [95% conf. interval]
------------------------------------------------------------------------------
    2005     2001  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2005     2002  pre        -.3205567   .1873138     -.687685   .04657163
    2005     2003  pre        .04395606    .162133    -.2738187    .3617308
    2005     2004  reference          0          0            0           0
    2005     2005  post        .8395668   .1596788     .5266021    1.152532
    2005     2006  post        2.123606   .1619789     1.806133    2.441078
    2005     2007  post        2.913114   .1614665     2.596646    3.229583
    2005     2008  post        4.134776   .1826639     3.776762    4.492791
    2005     2009  post        5.208695   .1974261     4.821746    5.595643
    2005     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
    2006     2001  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2006     2002  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2006     2003  pre        -.2362628   .1981004    -.6245324    .1520069
    2006     2004  pre        -.1639655   .1674376    -.4921373    .1642062
    2006     2005  reference          0          0            0           0
    2006     2006  post        1.235331   .1558061     .9299563    1.540705
    2006     2007  post         2.00499   .1643578     1.682855    2.327126
    2006     2008  post        3.121677   .1921316     2.745106    3.498248
    2006     2009  post        4.158905   .1891561     3.788166    4.529644
    2006     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
    2007     2001  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2007     2002  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2007     2003  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2007     2004  pre         .2055044    .183633    -.1544095    .5654184
    2007     2005  pre         .2483647   .1643524   -.07376004    .5704895
    2007     2006  reference          0          0            0           0
    2007     2007  post        1.404148    .153659     1.102982    1.705314
    2007     2008  post        2.273054   .1871435     1.906259    2.639849
    2007     2009  post        3.357806   .1979437     2.969844    3.745769
    2007     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
    2008     2001  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2008     2002  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2008     2003  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2008     2004  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2008     2005  pre       -.04116385   .1894952    -.4125676    .3302399
    2008     2006  pre        -.1185174   .1746093    -.4607454    .2237106
    2008     2007  reference          0          0            0           0
    2008     2008  post        1.300169   .1789554      .949423    1.650915
    2008     2009  post        2.391802   .1939501     2.011667    2.771937
    2008     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
    2009     2001  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2009     2002  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2009     2003  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2009     2004  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2009     2005  pre                .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this pre-treatment cell to the
      reference date
    2009     2006  pre       -.08856543   .2006478    -.4818279     .304697
    2009     2007  pre        .01184403   .2250112    -.4291699    .4528579
    2009     2008  reference          0          0            0           0
    2009     2009  post        .9958624   .2013366       .60125    1.390475
    2009     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
    2010     2001  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2001 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2002  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2002 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2003  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2003 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2004  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2004 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2005  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2005 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2006  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2006 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2007  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2007 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2008  pre                .          .            .           .
    Unavailable: no control cohort is untreated at date 2009, so the comparison of dates 2008 and
      2009 has no controls; cohorts in the data: 2005, 2006, 2007, 2008, 2009, 2010
    2010     2009  reference          0          0            0           0
    2010     2010  post               .          .            .           .
    Unavailable: the cross-pair comparisons do not connect this cell to the reference date
------------------------------------------------------------------------------
Std. err.: analytic; pointwise 95% intervals.
Reference cells are zero by construction.
Pretrend Wald test: not computed; some pre-treatment cells are unavailable, and the test is not
      computed on a subset of them.
```

Its event study, over the available cells:

```text
Chained difference-in-differences            Number of units       =       600
Method: doubly robust (dr)                   Number of event times =        15
Cross-pair comparisons; GMM weighting: within
Control group: not yet treated; base period: universal
Aggregation: dynamic
Event-time selection: narm
------------------------------------------------------------------------------
Event time                  ATT  Std. err.     [95% conf. interval]
------------------------------------------------------------------------------
Overall                3.167956   .1054286      2.96132    3.374592
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
                  -3 -.09620886  .08625618    -.2652679   .07285015
                  -2  .00433638  .08023198    -.1529154    .1615882
                  -1          0          0            0           0
                   0   1.155015  .07762391     1.002875    1.307156
                   1   2.198363  .09820502     2.005885    2.390841
                   2   3.130866   .1198263     2.896011    3.365721
                   3   4.146841   .1430075     3.866551     4.42713
                   4   5.208695   .1974261     4.821746    5.595643
                   5          .          .            .           .
    Unavailable: every cell of this aggregate is unavailable: cohort 2005 at date 2010
      (post-treatment)
------------------------------------------------------------------------------
Std. err.: analytic; pointwise 95% intervals.
The overall summary uses a pointwise interval.
Pretrend Wald test: not computed; some pre-treatment cells are unavailable, and the test is not
      computed on a subset of them.
```

The other two fits, not shown, report the same set of available cells. All three fits estimate the 15 post-treatment effects of cohorts 2005 to 2009, including ATT(2005, 2008) and ATT(2005, 2009). On contiguous windows like these, same-pair chaining already links the reference date to later dates through overlapping control cohorts, so `crosspair` is not needed for availability. The effects at 2010 are unavailable because every unit is treated by 2010, which leaves no untreated comparison.

## Inspect the control chain

`cdid_chain` shows how the effects of cohort 2005 are identified, including ATT(2005, 2009), which links 2004 to 2009 although no control cohort is observed at both dates:

```stata
estimates restore within
cdid_chain, cohort(2005)
```

```text
Identification chains (cross-pair, ptassumption(ctrl), gmmweighting(within))
Comparisons used by GMM: 84    Overidentifying restrictions: 20
Cells: 6 (5 estimable, 2 chained without a direct comparison, 1 unavailable)

   Group     Time  Role      Chain               ATT   Std. err.
---------------------------------------------------------------
    2005     2005  post      direct       .839566816  .159678826
           witness: 2004->2005 [g=2005] - 2004->2005 [ctrl{2006}]
    2005     2006  post      direct       2.12360559  .161978876
           witness: 2004->2005 [g=2005] + 2005->2006 [g=2005] - 2004->2005 [ctrl{2006}] - 2005->2006 [ctrl{2007}]
    2005     2007  post      direct        2.9131141  .161466526
           witness: 2004->2005 [g=2005] + 2005->2006 [g=2005] + 2006->2007 [g=2005] - 2004->2005 [ctrl{2006}] - 2005->2006 [ctrl{2007}] - 2006->2007 [ctrl{2008}]
    2005     2008  post      chained      4.13477612  .182663873
           witness: 2004->2005 [g=2005] + 2005->2006 [g=2005] + 2006->2007 [g=2005] + 2007->2008 [g=2005] - 2004->2005 [ctrl{2006}] - 2005->2006 [ctrl{2007}] - 2006->2007 [ctrl{2008}] - 2007->2008 [ctrl{2009}]
    2005     2009  post      chained      5.20869451  .197426141
           witness: 2004->2005 [g=2005] + 2005->2006 [g=2005] + 2006->2007 [g=2005] + 2007->2008 [g=2005] + 2008->2009 [g=2005] - 2004->2005 [ctrl{2006}] - 2005->2006 [ctrl{2007}] - 2006->2007 [ctrl{2008}] - 2007->2008 [ctrl{2009}] - 2008->2009 [ctrl{2010}]
    2005     2010  post      unavailable           .           .
           the cross-pair comparisons do not connect this cell to the reference date
```

ATT(2005, 2005) is identified directly by control cohorts observed in both 2004 and 2005; ATT(2005, 2008) and ATT(2005, 2009) require chaining, and the witness path shows one combination of control changes that links the dates. `cdid_chain, weights` lists the GMM weight of every comparison.

## Practical recommendations

1. **Keep irregular waves.** `cdid` chains over the stored study dates, so a switch from annual to biennial waves does not require dropping dates.
2. **Keep unit identifiers.** When attrition is related to persistent outcome levels, cohort fixed-effects regressions mix changes in composition with treatment effects, while within-unit differencing removes the persistent levels. The representativeness of the changes and parallel trends are still required.
3. **Choose the assumption and the weighting.** `gmmweighting(within) ptassumption(ctrl)` reproduces the stacked regressions. `ptassumption(post)` with the default `gmmweighting(optimal)` avoids the PT-ctrl restriction before the reference date, whose violation by a control cohort would bias the post-treatment estimates, and estimates the covariance of the comparisons from the data.
4. **Adjust when needed.** Covariates, sampling weights that vary over time (`fixweights()`), and covariate-dependent attrition (`selection(conditional_mtar)` with any `selmethod()`) are supported with `crosspair`; see [Cross-pair chaining](../methods/cross-pair.md).

---

**Reading path, step 22 of 22** · Previous: [Cross-pair chaining, control parallel trends and within-unit weighting](../methods/cross-pair.md)

[Guide home](../README.md) · [Package overview](../../README.md)
