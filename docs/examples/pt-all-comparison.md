# Add PT-all restrictions

`ptassumption(post)` leaves the differential trends before the reference period unrestricted. `ptassumption(all) ptstart(#)` keeps PT-post and each cohort's reference period and adds the restriction that the differential trends are zero on every step from the common calendar date `ptstart()` onward. `baseperiod()` separately controls how pretreatment cells are reported.

This page compares the two assumptions in a design where PT-all is false. When parallel trends hold in all periods, both sets of restrictions are valid. PT-all can add information or identification, but no finite-sample precision improvement is guaranteed.

## A pre-reference pulse

Take cohorts first treated at 3 and 5, no anticipation, and reference dates 2 and 4. Add a common $+0.5$ untreated pulse to the treated cohorts at date 1 only. All outcomes at dates 2 onward, all post-treatment causal effects, and PT-post remain unchanged. The differential on step $1\to2$ changes by $-0.5$, so PT-all from date 1 imposes a false restriction. Under PT-post, the earlier differential trends remain free parameters and cancel in the valid comparisons through earlier dates.

These statements concern population moments. They do not imply that all finite-sample post estimates are exactly unchanged, that a Wald test must reject, or that every PT-all estimate moves in one direction. Shared adjustment models, the estimated covariance of the comparisons and sampling variation all matter.

## Estimation with and without the pulse

Run Sections 0 and 1 of the hands-on do-file first. Sections 3a and 3c estimate the model under PT-post and PT-all; the following two calls add the pulse:

```stata
generate double y_pulse = y + 0.5 * (first_treat > 0 & period == 1) if !missing(y)

cdid y_pulse x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) ptassumption(post) bstrap(0) cband(0)

cdid y_pulse x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) ptassumption(all) ptstart(1) bstrap(0) cband(0)
```

The simulated observation process is not accounted for by these specifications, which use `selection(none)`. The pulse leaves the observed outcomes in place; the example shows how sensitive the estimates are to the added restrictions and gives no evidence of consistency or coverage for this sample. Choose the period over which parallel trends are assumed from substantive knowledge, before looking at the results.

---

**Reading path, step 14 of 22** · Previous: [Parallel trends](../methods/parallel-trends.md) · Next: [GMM weighting and covariance regularization](../methods/gmm.md)

[Guide home](../README.md) · [Package overview](../../README.md)
