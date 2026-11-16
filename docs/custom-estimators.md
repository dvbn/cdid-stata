# Custom estimators

`cdid` and `cdid_pretest` can use your own estimator of each comparison instead of the built-in IPW, OR and DR estimators. The estimator is a Mata function, the *callback*, named in the options

```text
method(custom) callback(my_att) [payload(my_settings)]
```

where `my_settings` is an optional Mata object passed to the callback by pointer.

## 1. Before compiling a callback

The structures `cdid_callback_request` and `cdid_callback_result` and the function `cdid_callback_result_make()` are defined in the library `lcdid_native.mlib`, which is installed with `cdid`. In a new Stata session, index the library before compiling your callback:

```stata
mata: mata mlib index
```

## 2. Writing a callback

The callback receives a `struct cdid_callback_request scalar`, describing one comparison, and a `pointer() scalar` to the payload (`NULL` when `payload()` is not specified). It returns the result of `cdid_callback_result_make(att, psi, gradient)`, where `att` is the estimate, `psi` is the `m x 1` influence function of the `m` observations passed, and `gradient` is the `m x 1` derivative of `att` with respect to the log weights, or `J(0, 1, .)` when it is not supplied. The callback signals that it cannot estimate the comparison by raising an error, for example `_error(459)`. The error stops chained estimation; in ordinary mode, the cell is reported as unavailable ("the callback stopped with an error"), and in the ordinary-mode `cdid_pretest`, the test is reported as unavailable.

The two callbacks below are used in the [custom-estimator do-file](../examples/cdid_custom_estimators.do). `my_att` estimates the weighted difference in mean outcome changes for panel comparisons; `my_rcs_att` estimates the difference in differences of weighted means for comparisons of repeated cross sections. Both count their calls in the payload.

```stata
version 19.5
mata:
struct cdid_callback_result scalar my_att(
    struct cdid_callback_request scalar r, pointer() scalar payload)
{
    real scalar m, mt, mc
    real colvector dy, wt, wc, psi
    if (payload != NULL) *payload = *payload + 1
    if (r.panel != 1) _error(498)
    m = rows(r.D)
    dy = r.y1 - r.y0
    wt = r.i_weights :* r.D
    wc = r.i_weights :* (1 :- r.D)
    if (sum(wt) <= 0 | sum(wc) <= 0) _error(459)
    mt = sum(wt :* dy) / sum(wt)
    mc = sum(wc :* dy) / sum(wc)
    psi = m :* (wt :* (dy :- mt) / sum(wt) - wc :* (dy :- mc) / sum(wc))
    return(cdid_callback_result_make(mt - mc, psi, psi / m))
}

struct cdid_callback_result scalar my_rcs_att(
    struct cdid_callback_request scalar r, pointer() scalar payload)
{
    real scalar m, d, p, sign, mu, att
    real colvector weight, psi
    if (payload != NULL) *payload = *payload + 1
    if (r.panel != 0) _error(498)
    m = rows(r.D); att = 0; psi = J(m, 1, 0)
    for (d = 0; d <= 1; d++) for (p = 0; p <= 1; p++) {
        weight = r.i_weights :* (r.D :== d) :* (r.post :== p)
        if (sum(weight) <= 0) _error(459)
        mu = sum(weight :* r.y) / sum(weight)
        sign = (d == p ? 1 : -1)
        att = att + sign * mu
        psi = psi + sign * m :* weight :* (r.y :- mu) / sum(weight)
    }
    return(cdid_callback_result_make(att, psi, J(0, 1, .)))
}
end
```

Use the callbacks on a simulated panel, after running the block above:

```stata
clear
set seed 20261116
set obs 600
generate id = _n
generate G = cond(runiform() < .4, 0, cond(runiform() < .5, 2004, 2006))
generate x = rnormal()
expand 8
bysort id: generate year = 2000 + _n
generate y = x + .2*(year - 2000) + (G > 0 & year >= G)*(1 + .25*(year - G)) + rnormal()
replace y = . if runiform() < .2
mata: my_calls = 0
cdid y, ivar(id) time(year) gvar(G) method(custom) callback(my_att) payload(my_calls) bstrap(0) cband(0)
mata: my_calls
estat simple
```

The callbacks ignore covariates, so the varlist has none. The [custom-estimator do-file](../examples/cdid_custom_estimators.do) also runs them in ordinary mode and with `cdid_pretest`.

## 3. Fields of the request

| Field | Type | Description |
| --- | --- | --- |
| `mode` | `string scalar` | `"chained"` or `"ordinary"` |
| `route` | `string scalar` | Label of the type of comparison, for information: `"baseline_static"`, `"baseline_dated"` (covariates or weights that vary over time), `"conditional_mtar_static"` or `"conditional_mtar_dated"` in chained mode; `"panel"`, `"unbalanced_panel"` or `"RCS"` in ordinary mode. Branch on `mode` and `panel` instead |
| `panel` | `real scalar` | `1` when outcomes are passed as pairs (`y0`, `y1`); `0` when they are stacked (`y`, `post`) |
| `inffunc` | `real scalar` | `1` when `psi` is required; `0` when only the estimate is required (ordinary mode with `noinffunc`) |
| `gradient_required` | `real scalar` | `1` when the gradient is required (chained observation weighting); `0` otherwise |
| `g`, `s`, `t` | `real scalar` | Cohort, base period, and comparison period |
| `D` | `real colvector` | `m x 1` indicator: `1` for units of cohort `g`, `0` for comparison units |
| `y0`, `y1` | `real colvector` | `m x 1` outcomes in the base and comparison periods (`panel == 1`) |
| `y`, `post` | `real colvector` | `m x 1` outcome and indicator of the comparison period (`panel == 0`) |
| `i_weights` | `real colvector` | `m x 1` observation weights |
| `covariates` | `real matrix` | `m x p` covariate matrix |
| `covariate_names` | `string rowvector` | `1 x p` names of the covariates |
| `rows_F`, `rows_A` | `real colvector` | Positions of the observations passed among the units of the study population and of the estimation sample; `rows_A` is empty in ordinary mode |
| `call_index` | `real scalar` | Number of this call, starting at 1 |

## 4. Where callbacks are used

1. **Chained estimation without a missing-outcome correction** (`selection(none)`). `att` and `psi` must be finite. The gradient is optional; when it is omitted, `cdid` uses `psi / m`. When it is supplied, it must be a finite `m x 1` vector.
2. **Chained estimation with observation weighting** (`selection(conditional_mtar) selmethod(ipw)`), with any observation model, pooling and fitting universe. The gradient is required: the total derivative of `att` with respect to the log weight of each observation passed, in the units of `att` (not divided by `m` or `N`). `selmethod(imputation)` and `selmethod(dr)` do not support custom estimators.
3. **Ordinary estimation** (`mode(ordinary)` or `rcs`). The callback estimates each group-time effect from the two periods compared, in balanced panels (`panel == 1`), unbalanced panels and repeated cross sections (`panel == 0`). With `noinffunc`, `inffunc` is 0. The gradient is ignored. The influence function is used as returned; in particular, the adjustment that the built-in doubly robust estimator applies for repeated cross sections is not applied.
4. **Conditional pretrend test, ordinary mode.** `cdid_pretest` calls the callback for each evaluation point and pretreatment cell, with the outcome replaced by the indicator-transformed outcome used by the test.
5. **Conditional pretrend test, chained mode.** The test is not defined for custom estimators. `cdid_pretest` returns `r(status) = "not_applicable"` with an explanation in `r(reason)` and does not call the callback.
6. **Covariates missing at a reference date.** When `missingreference(exclude)` restricts a cohort, the callback receives only the units eligible for that cohort's comparisons, with weights computed on those units.
7. **Not supported.** `method(custom)` is not available with `crosspair`, `gmmweighting(within)`, `preset()` or `dvar()`. `cdid_chain` reports the comparisons and GMM weights of custom fits.
8. **After estimation.** `cdid` stores the estimates, the influence functions, and the names of the callback and payload in `e(callback)` and `e(payload)`; `cdid_pretest` returns the names in `r(callback)` and `r(payload)`, with the number of calls in `r(n_callbacks)`, and `cdid_pretest_load` does not restore them. `estat`, `cdid_plot`, `cdid_export`, `cdid_load` and `cdid_pretest_load` use the stored results and never call the callback again; it need not exist after estimation.

See `help cdid_custom` for the same information in Stata.

[Guide home](README.md)
