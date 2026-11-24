* Run in a fresh Stata session with cdid installed in a project ado directory.
* From outside the source checkout: do /path/to/validation/weighting-default.do
version 19.5
clear all
set more off
which cdid, all
findfile lcdid_native.mlib
quietly mata: mata mlib index

* The staggered panel in docs/from-other-packages.md, with sampling weights
* and fixed missing-outcome patterns for the clustered comparison.
set seed 14603
set obs 144
generate long id = _n
generate byte g = cond(mod(id - 1, 3) == 0, 0, cond(mod(id - 1, 3) == 1, 3, 5))
generate double a = rnormal()
generate double x = rnormal()
generate double w = exp(.2 * x)
generate long cl = ceil(id / 6)
expand 6
bysort id: generate byte t = _n
generate double y = a + sin(t) + .15*x*t + rnormal() + cond(g > 0 & t >= g, (t - g + 1) / 2, 0)

local panel_opts ivar(id) time(t) gvar(g) bstrap(0) cband(0)
tempname b V se map
tempfile omitted_if within_if

foreach panel in complete incomplete {
    local model_opts controlgroup(nevertreated)
    if "`panel'" == "incomplete" {
        quietly replace y = . if mod(id + 2*t, 7) == 0
        local model_opts controlgroup(notyettreated) wvar(w) cluster(cl) ///
            selection(conditional_mtar) selmethod(dr)
    }
    foreach method in ipw reg dr {
        local opts `panel_opts' `model_opts' method(`method')
        quietly cdid y x, `opts'
        assert "`e(gmm_weighting)'" == "within"
        matrix `b' = e(b)
        matrix `V' = e(V)
        matrix `se' = e(se)
        matrix `map' = e(coef_map)
        if "`panel'" == "incomplete" & "`method'" == "dr" {
            quietly cdid_export using "`omitted_if'", type(influence)
        }
        quietly cdid y x, `opts' gmmweighting(within)
        assert mreldif(`b', e(b)) == 0
        assert mreldif(`V', e(V)) == 0
        assert mreldif(`se', e(se)) == 0
        assert mreldif(`map', e(coef_map)) == 0
        if "`panel'" == "incomplete" & "`method'" == "dr" {
            assert e(N_clust) == 24
            quietly cdid_export using "`within_if'", type(influence)
            preserve
            quietly use "`within_if'", clear
            cf _all using "`omitted_if'"
            restore
        }
    }
}

* The two spellings of feasible two-step weighting agree.
quietly cdid y x, `panel_opts' gmmweighting(optimal)
assert "`e(gmm_weighting)'" == "optimal"
matrix `b' = e(b)
matrix `V' = e(V)
quietly cdid y x, `panel_opts' gmmweighting(two_step)
assert "`e(gmm_weighting)'" == "two_step"
assert mreldif(`b', e(b)) == 0
assert mreldif(`V', e(V)) == 0
quietly cdid y x, `panel_opts' gmmweighting(identity)
assert "`e(gmm_weighting)'" == "identity"

* Weighted difference in mean changes, as in cdid_custom_estimators.do.
mata:
struct cdid_callback_result scalar weighting_test_att(
    struct cdid_callback_request scalar r, pointer() scalar payload)
{
    real scalar treated_mean, control_mean
    real colvector dy, wt, wc, gradient
    dy = r.y1 - r.y0
    wt = r.i_weights :* r.D
    wc = r.i_weights :* (1 :- r.D)
    treated_mean = sum(wt :* dy) / sum(wt)
    control_mean = sum(wc :* dy) / sum(wc)
    gradient = wt :* (dy :- treated_mean) / sum(wt) -
        wc :* (dy :- control_mean) / sum(wc)
    return(cdid_callback_result_make(treated_mean - control_mean, rows(r.D) :* gradient, gradient))
}
end

local custom `panel_opts' method(custom) callback(weighting_test_att)
quietly cdid y, `custom'
assert "`e(gmm_weighting)'" == "optimal"
capture noisily cdid y, `custom' gmmweighting(within)
assert _rc == 198

local preset `panel_opts' preset(cohort_time_fe) baseperiod(universal)
foreach method in ipw dr reg {
    local expected = cond("`method'" == "reg", "within", "optimal")
    local refused = cond("`method'" == "reg", "optimal", "within")
    quietly cdid y x, `preset' method(`method')
    assert "`e(gmm_weighting)'" == "`expected'"
    capture noisily cdid y x, `preset' method(`method') gmmweighting(`refused')
    assert _rc == 198
}

* The four two-date switching histories in static-observation-separation.do.
quietly keep if t <= 2
generate byte role = mod(id, 4)
generate byte d = cond(t == 1, inlist(role, 2, 3), inlist(role, 1, 2))
local switching ivar(id) time(t) dvar(d) method(reg) bstrap(0) cband(0)
quietly cdid y x, `switching'
assert "`e(gmm_weighting)'" == "optimal"
matrix `b' = e(b)
matrix `V' = e(V)
quietly cdid y x, `switching' gmmweighting(optimal)
assert mreldif(`b', e(b)) == 0
assert mreldif(`V', e(V)) == 0
capture noisily cdid y x, `switching' gmmweighting(within)
assert _rc == 198

* Four working-covariance moments can be full rank with three study units.
mata:
real matrix weighting_test_if(pointer() scalar context, real scalar lo, real scalar hi)
{
    return((*context)[lo..hi, .])
}

void weighting_test_rank_bound()
{
    struct cdid_gmm_result scalar within, empirical
    struct cdid_gmm_unit_gram scalar source
    struct cdid_xp_fit scalar fit
    struct cdid_xp_design scalar design
    real colvector r, A, expected_if
    real rowvector loadings
    real matrix Vw, phi
    real scalar cap

    r = (1 \ 2 \ 6 \ 10)
    A = J(4, 1, 1)
    Vw = diag((1, 2, 4, 8))
    phi = (1, 0, 1, 1 \ 0, 1, -1, 1 \ -1, -1, 0, -2)
    loadings = (8, 4, 2, 1) / 15
    expected_if = phi * loadings'
    cap = 10000
    within = cdid_gmm_fit_within_norm(r, A, Vw, I(4), 4, 1, 1e-10, 1e-4,
        3, 3, 2, &weighting_test_if(), &phi, cap)
    if (within.lambda != 0 | within.handler_invoked | within.covariance_rank != 4 |
        within.sample_rank_upper != 2)
        _cdid_gmm_fail("Full-rank within weighting used the empirical rank bound.", 9)
    if (abs(within.targets - 38 / 15) > 1e-12 |
        max(abs(within.target_loadings - loadings)) > 1e-12 |
        !all(within.retained_if :== phi) |
        abs(within.covariance - quadcross(expected_if, expected_if) / 3) > 1e-12)
        _cdid_gmm_fail("Within estimate or original-unit influence functions disagree with weighted least squares.", 9)

    empirical = cdid_gmm_fit_unit_norm(r, A, I(4), 4, 1, "optimal", 1e-10, 1e-4,
        3, 3, 2, &weighting_test_if(), &phi, cap)
    if (empirical.lambda <= 0 | empirical.lambda >= . | !empirical.handler_invoked)
        _cdid_gmm_fail("Singular empirical weighting covariance did not receive its ridge.", 9)

    source = _cdid_gmm_unit_gram_norm(3, 4, 3, 2, &weighting_test_if(), &phi, cap)
    source.covariance = Vw
    fit.mhat = r
    fit.V = Vw
    fit.Phi = phi
    fit.att = fit.avail = J(1, 1, .)
    fit.Lam = J(4, 1, .)
    fit.IF = J(3, 1, .)
    design.A = A
    _cdid_xp_gmm(fit, design, 1, 1, I(4), 4, "within", 1e-10, 1e-4, source, 1)
    if (fit.diag[1] != 0 | abs(fit.att - 38 / 15) > 1e-12 |
        max(abs(fit.IF - expected_if)) > 1e-12)
        _cdid_gmm_fail("Cross-pair within weighting disagrees with weighted least squares.", 9)
}
weighting_test_rank_bound()
end

display "CDID_WEIGHTING_DEFAULT_PASS"
