*! Display and tables of cdid_chain.
version 19.5
mata:

real scalar _cdid_chain_rrank(real scalar role)
{
    return(role == 1 ? 1 : (role == 3 ? 2 : 3))
}

string scalar _cdid_chain_tag(struct cdid_chain_info scalar I, real scalar k)
{
    if (I.kind == "crosspair") {
        if (I.links[k, 4] == 1) return("g=" + _cdid_chain_lbl(I.links[k, 1]))
        return("ctrl" + I.pool[k])
    }
    return("K=" + I.pool[k])
}

/* Witness rows of reported cell r, ordered by role, s and t: (link, weight). */
real matrix _cdid_chain_wrows(struct cdid_chain_info scalar I, real scalar r)
{
    real colvector h, rr
    real matrix W
    h = selectindex(I.wit[., 1] :== r)
    if (rows(h) == 0) return(J(0, 2, .))
    W = I.wit[h, 2..3]
    rr = J(rows(W), 1, .)
    for (h = 1; h <= rows(W); h++) rr[h] = _cdid_chain_rrank(I.links[W[h, 1], 4])
    return(W[order((rr, I.links[W[., 1], 2..3]), (1, 2, 3)), .])
}

/* GMM weights of reported cell r above tol, ordered by other cohort, role, s and t. */
real matrix _cdid_chain_grows(struct cdid_chain_info scalar I, real scalar r, real scalar tol)
{
    real colvector lam, act, rr, og, h
    real scalar k
    if (I.cells[r, 8] < 1) {
        h = selectindex((I.wit[., 1] :== r) :& (I.wit[., 2] :> I.n_moments))
        return(rows(h) ? I.wit[h, 2..3] : J(0, 2, .))
    }
    lam = I.Lam[., I.cells[r, 8]]
    act = selectindex(abs(lam) :> tol)
    if (rows(act) == 0) return(J(0, 2, .))
    rr = J(rows(act), 1, .)
    for (k = 1; k <= rows(act); k++) rr[k] = _cdid_chain_rrank(I.links[act[k], 4])
    og = (I.links[act, 1] :!= I.cells[r, 1])
    act = act[order((og, rr, I.links[act, 2..3], act), (1, 2, 3, 4, 5))]
    return((act, lam[act]))
}

string scalar _cdid_chain_path(struct cdid_chain_info scalar I, real matrix W)
{
    string scalar out, sg
    real scalar h, k
    out = ""
    for (h = 1; h <= rows(W); h++) {
        k = W[h, 1]
        sg = (W[h, 2] >= 0 ? "+" : "-")
        if (h > 1) out = out + " "
        if (!(h == 1 & sg == "+")) out = out + sg + " "
        out = out + _cdid_chain_lbl(I.links[k, 2]) + "->" + _cdid_chain_lbl(I.links[k, 3]) + " [" + _cdid_chain_tag(I, k) + "]"
    }
    return(out)
}

real colvector _cdid_chain_select(struct cdid_chain_info scalar I, string scalar cohorts,
    string scalar times, real scalar pre)
{
    real colvector keep, hit
    real rowvector cv, tv
    real scalar i
    keep = J(rows(I.cells), 1, 1)
    if (rows(I.cells) == 0) return(J(0, 1, .))
    if (!pre & strtrim(times) == "") keep = keep :& (I.cells[., 5] :== 3)
    if (strtrim(cohorts) != "") {
        cv = strtoreal(tokens(cohorts))
        hit = J(rows(I.cells), 1, 0)
        for (i = 1; i <= cols(cv); i++) hit = hit :| (I.cells[., 1] :== cv[i])
        keep = keep :& hit
    }
    if (strtrim(times) != "") {
        tv = strtoreal(tokens(times))
        hit = J(rows(I.cells), 1, 0)
        for (i = 1; i <= cols(tv); i++) hit = hit :| (I.cells[., 2] :== tv[i])
        keep = keep :& hit
    }
    return(selectindex(keep))
}

struct cdid_chain_info scalar _cdid_chain_current()
{
    struct cdid_chain_info scalar I
    if (!_cdid_chain_get(st_global("e(state_token)"), I)) {
        errprintf("chain information is not available for the current estimation results\n")
        errprintf("cdid_chain needs the results of cdid in chained mode, the estat results computed from them in\n")
        errprintf("this session, or results loaded with cdid_load from a file saved by cdid; refit with cdid\n")
        exit(498)
    }
    return(I)
}

real scalar _cdid_chain_overid(struct cdid_chain_info scalar I)
{
    return(I.n_params < . ? I.n_moments - I.n_params : .)
}

void cdid_chain_show(string scalar cohorts, string scalar times, real scalar pre, real scalar tol,
    real scalar weights)
{
    struct cdid_chain_info scalar I
    struct cdid_ur_view scalar V
    struct cdid_ur_cells scalar C
    real colvector sel
    real matrix W, Gw
    real scalar i, r, h, k, nav, nch
    string scalar role, chain
    I = _cdid_chain_current()
    V = cdid_ur_open()
    C = cdid_ur_cells_view(V)
    if (rows(C.estimate) != rows(I.cells)) _cdid_rs_fail("Chain information does not align with the reported cells.", 498)
    sel = _cdid_chain_select(I, cohorts, times, pre)
    nav = sum(I.cells[sel, 6] :== 1)
    nch = sum((I.cells[sel, 6] :== 1) :& (I.cells[sel, 7] :== 1))
    printf("\n{txt}Identification chains (%s, ptassumption(%s), gmmweighting(%s))\n",
        (I.kind == "crosspair" ? "cross-pair" : "same-pair"), I.pt, I.weighting)
    printf("{txt}Comparisons used by GMM: {res}%g{txt}", I.n_moments)
    if (_cdid_chain_overid(I) < .) printf("{txt}    Overidentifying restrictions: {res}%g{txt}", _cdid_chain_overid(I))
    printf("\n{txt}Cells: {res}%g{txt} ({res}%g{txt} estimable, {res}%g{txt} chained without a direct comparison, {res}%g{txt} unavailable)\n\n",
        rows(sel), nav, nch, rows(sel) - nav)
    if (rows(sel) == 0) return
    printf("{txt}%8s %8s  %-9s %-11s %11s %11s\n", "Group", "Time", "Role", "Chain", "ATT", "Std. err.")
    printf("{hline 63}\n")
    for (i = 1; i <= rows(sel); i++) {
        r = sel[i]
        role = (I.cells[r, 5] == 1 ? "reference" : (I.cells[r, 5] == 2 ? "pre" : "post"))
        chain = (I.cells[r, 6] != 1 ? "unavailable" : (I.cells[r, 7] ? "chained" : "direct"))
        printf("{txt}%8s %8s  %-9s %-11s {res}%11.0g %11.0g\n", _cdid_chain_lbl(I.cells[r, 1]),
            _cdid_chain_lbl(I.cells[r, 2]), role, chain, C.estimate[r], C.se[r])
        W = _cdid_chain_wrows(I, r)
        if (rows(W)) printf("{txt}%10s witness: %s\n", "", _cdid_chain_path(I, W))
        if (I.cells[r, 6] != 1 & C.reason[r] != "") printf("{p 11 13 2}{txt}%s{p_end}\n", _cdid_ur_smcl(cdid_text_unavailable(C.reason[r])))
        if (weights & I.cells[r, 6] == 1) {
            Gw = _cdid_chain_grows(I, r, tol)
            if (rows(Gw)) {
                printf("{txt}%10s %6s %6s %6s  %-11s %10s %10s %10s\n", "", "Cohort", "From", "To", "Comparison", "Weight", "Estimate", "Contrib.")
                for (h = 1; h <= rows(Gw); h++) {
                    k = Gw[h, 1]
                    printf("{txt}%10s %6s %6s %6s  %-11s {res}%10.0g %10.0g %10.0g\n", "",
                        _cdid_chain_lbl(I.links[k, 1]), _cdid_chain_lbl(I.links[k, 2]), _cdid_chain_lbl(I.links[k, 3]),
                        _cdid_chain_tag(I, k), Gw[h, 2], I.mhat[k], Gw[h, 2] * I.mhat[k])
                }
            }
        }
    }
    st_numscalar("r(n_cells)", rows(sel))
    st_numscalar("r(n_available)", nav)
    st_numscalar("r(n_chained)", nch)
    st_numscalar("r(n_moments)", I.n_moments)
    if (_cdid_chain_overid(I) < .) st_numscalar("r(overid_df)", _cdid_chain_overid(I))
    st_global("r(kind)", I.kind)
}

/* Fill the (empty) current frame with one chain table: summary, witness or weights. */
void cdid_chain_table(string scalar which, string scalar cohorts, string scalar times, real scalar pre,
    real scalar tol)
{
    struct cdid_chain_info scalar I
    struct cdid_ur_view scalar V
    struct cdid_ur_cells scalar C
    real colvector sel
    real matrix W, out
    string matrix txt
    real scalar i, r, h, k, n
    I = _cdid_chain_current()
    V = cdid_ur_open()
    C = cdid_ur_cells_view(V)
    sel = _cdid_chain_select(I, cohorts, times, pre)
    if (which == "summary") {
        n = rows(sel)
        out = J(n, 11, .)
        txt = J(n, 3, "")
        for (i = 1; i <= n; i++) {
            r = sel[i]
            W = _cdid_chain_wrows(I, r)
            out[i, .] = (I.cells[r, 1], I.cells[r, 2], I.cells[r, 2] - I.cells[r, 1], I.cells[r, 3], I.cells[r, 4],
                I.cells[r, 6] == 1, (I.cells[r, 6] == 1) & I.cells[r, 7], rows(W), rows(_cdid_chain_grows(I, r, tol)),
                C.estimate[r], C.se[r])
            txt[i, .] = ((I.cells[r, 5] == 1 ? "reference" : (I.cells[r, 5] == 2 ? "pre" : "post")),
                _cdid_chain_path(I, W), cdid_text_unavailable(C.reason[r]))
        }
        _cdid_chain_store(("g", "time", "event_time", "ref", "base", "available", "requires_chaining",
            "n_witness", "n_active", "estimate", "se"), out, ("role", "chain_path", "reason"), txt)
        return
    }
    out = J(0, 8, .)
    txt = J(0, 2, "")
    for (i = 1; i <= rows(sel); i++) {
        r = sel[i]
        W = (which == "witness" ? _cdid_chain_wrows(I, r) : _cdid_chain_grows(I, r, tol))
        for (h = 1; h <= rows(W); h++) {
            k = W[h, 1]
            out = out \ (I.cells[r, 1], I.cells[r, 2], (k <= I.n_moments ? k : .), I.links[k, 1..3], W[h, 2], I.mhat[k])
            txt = txt \ ((I.links[k, 4] == 1 ? "target" : (I.links[k, 4] == 2 ? "control" : (I.links[k, 4] == 3 ? "did_pair" : "auxiliary_did_pair"))), I.pool[k])
        }
    }
    out = out, out[., 7] :* out[., 8]
    _cdid_chain_store(("g", "time", "link", "link_g", "s", "t", (which == "witness" ? "weight" : "gmm_weight"),
        "link_estimate", "contribution"), out, ("link_role", "control_pool"), txt)
}

void _cdid_chain_store(string rowvector nn, real matrix X, string rowvector sn, string matrix S)
{
    real scalar j, n
    n = rows(X)
    st_addobs(n)
    for (j = 1; j <= cols(nn); j++) {
        (void) st_addvar("double", nn[j])
        if (n) st_store(., nn[j], X[., j])
    }
    for (j = 1; j <= cols(sn); j++) {
        (void) st_addvar("strL", sn[j])
        if (n) st_sstore(., sn[j], S[., j])
    }
}

end
