*! Within-unit GMM weighting for same-pair chained comparisons.
version 19.5
mata:

struct _cdid_wsp_vec {
    real colvector v
}

/* Analysis-unit data for the within covariance: observed-outcome grid on the catalogue
   calendar, cohort indices, weights, cohort weight mass, and the covariate design. */
struct _cdid_wsp_data {
    real scalar nA, sumw
    real matrix O, Z
    real colvector cidx, w, mass
}

struct _cdid_wsp_data scalar _cdid_wsp_prepare(struct cdid_cat scalar C,
    struct cdid_population scalar P, real matrix X, real colvector w_owner)
{
    struct _cdid_wsp_data scalar d
    real colvector aF, k, rr
    real scalar T, t, j, same
    aF = P.A_to_F
    d.nA = rows(aF)
    T = rows(C.calendar)
    same = (rows(P.calendar) == T & rows(P.obs_grid) == P.N & cols(P.obs_grid) == T)
    if (same) same = all(P.calendar :== C.calendar)
    d.O = J(d.nA, T, 0)
    for (t = 1; t <= T; t++) {
        rr = (same ? P.obs_grid[aF, t] : cdid_population_rows(P, aF, J(d.nA, 1, C.calendar[t])))
        k = selectindex(rr :> 0)
        if (rows(k)) d.O[k, t] = (P.y[rr[k]] :< .)
    }
    d.w = w_owner[aF]
    d.cidx = _cdid_cohort_index_vec(C.cohorts, C.nt_index, P.cohort[aF], J(d.nA, 1, 1))
    d.mass = J(rows(C.cohorts), 1, 0)
    for (j = 1; j <= rows(C.cohorts); j++) d.mass[j] = sum(select(d.w, d.cidx :== j))
    d.sumw = sum(d.w)
    d.Z = (cols(X) ? X[aF, .] : J(d.nA, 1, 1))
    return(d)
}

real rowvector _cdid_wsp_pool(struct cdid_cat scalar C, real scalar ell)
{
    real scalar p
    p = C.r_pool[ell]
    return(C.pool_entries[|C.pool_start[p] \ C.pool_start[p] + C.pool_size[p] - 1|]')
}

/* Units of comparison ell: target units, then control units, among the units observed
   at both dates. */
void _cdid_wsp_units(struct _cdid_wsp_data scalar d, struct cdid_cat scalar C, real scalar ell,
    real colvector idx_t, real colvector idx_c)
{
    real colvector obs, inK
    real rowvector K
    real scalar k
    obs = d.O[., C.r_s[ell]] :& d.O[., C.r_t[ell]]
    K = _cdid_wsp_pool(C, ell)
    inK = J(d.nA, 1, 0)
    for (k = 1; k <= cols(K); k++) inK = inK :| (d.cidx :== K[k])
    idx_t = selectindex(obs :& (d.cidx :== C.r_g[ell]))
    idx_c = selectindex(obs :& inK)
}

/* Unit loadings of comparison ell on the outcome change (R compute_pipeline_within_covariance):
   observation-rate weights by cohort, the target mean, and the IPW, OR or DR control
   loading. */
void _cdid_wsp_loadings(struct _cdid_wsp_data scalar d, struct cdid_cat scalar C,
    real scalar ell, string scalar method, real scalar eps, real scalar maxit,
    real colvector idx, real colvector om)
{
    struct cdid_logit_result scalar lf
    real colvector idx_t, idx_c, v_t, v_c, sel, b_c, w_ipw, om_c, piv, xbar_t, xatt, cx, D_u, cpos
    real matrix X_c, X_t, Zs, lv, Rk, H, Rf
    real rowvector K, tau, pv
    real scalar k, r_k, is_or, is_dr, rk, mx
    real colvector dd
    _cdid_wsp_units(d, C, ell, idx_t, idx_c)
    is_or = (method == "reg")
    is_dr = (method == "dr")
    v_t = d.w[idx_t] :/ max((sum(d.w[idx_t]) / d.mass[C.r_g[ell]], 1e-6))
    v_c = J(rows(idx_c), 1, 0)
    K = _cdid_wsp_pool(C, ell)
    for (k = 1; k <= cols(K); k++) {
        sel = selectindex(d.cidx[idx_c] :== K[k])
        if (rows(sel) == 0) continue
        r_k = max((sum(d.w[idx_c[sel]]) / d.mass[K[k]], 1e-6))
        v_c[sel] = d.w[idx_c[sel]] :/ r_k
    }
    if (!is_or) {
        Zs = d.Z[idx_t \ idx_c, .]
        if (cols(Zs) == 1 & all(Zs[., 1] :== Zs[1, 1]) & Zs[1, 1] != 0) b_c = v_c :* (sum(v_t) / sum(v_c))
        else {
            D_u = J(rows(idx_t), 1, 1) \ J(rows(idx_c), 1, 0)
            lf = cdid_logit_fit(Zs, D_u, v_t \ v_c, eps, maxit, 1)
            cpos = (rows(idx_t) + 1)::rows(D_u)
            lv = _cdid_ipw_logit_values(Zs[cpos, .] * lf.beta)
            b_c = v_c :* lv[., 2]
        }
        w_ipw = -b_c :/ sum(b_c)
    }
    if (!is_or & !is_dr) om_c = w_ipw
    else {
        X_c = d.Z[idx_c, .]
        X_t = d.Z[idx_t, .]
        hqrdp(X_c :* sqrt(v_c), H = ., tau = ., Rf = ., pv = .)
        dd = abs(diagonal(Rf))
        mx = max(dd)
        rk = sum(dd :> 1e-7 * mx)
        if (rk < 1) _cdid_ipw_fail("The outcome-regression control design has zero numerical rank on observed pair support.", 459)
        piv = pv[1..rk]'
        Rk = Rf[1..rk, 1..rk]
        xbar_t = (colsum(X_t[., piv'] :* v_t) :/ sum(v_t))'
        if (is_or) {
            cx = solveupper(Rk, solvelower(Rk', xbar_t))
            om_c = -v_c :* (X_c[., piv'] * cx)
        }
        else {
            xatt = (colsum(X_c[., piv'] :* b_c) :/ sum(b_c))'
            cx = solveupper(Rk, solvelower(Rk', xbar_t - xatt))
            om_c = w_ipw - v_c :* (X_c[., piv'] * cx)
        }
    }
    idx = idx_t \ idx_c
    om = (v_t :/ sum(v_t)) \ om_c
}

/* Within-unit levels covariance of comparisons rows (catalogue indices): entries
   sum_w * sum_i om_ai om_bi / w_i, signed by shared dates, for comparisons of the same
   target cohort. */
real matrix cdid_wsp_cov(struct _cdid_wsp_data scalar d, struct cdid_cat scalar C,
    real colvector rows_used, string scalar method, real scalar eps, real scalar maxit)
{
    struct _cdid_wsp_vec colvector lidx, lom
    real colvector mrk, inv
    real matrix V
    real scalar L, a, b, ea, eb, sg
    L = rows(rows_used)
    lidx = _cdid_wsp_vec(L, 1)
    lom = _cdid_wsp_vec(L, 1)
    for (a = 1; a <= L; a++) _cdid_wsp_loadings(d, C, rows_used[a], method, eps, maxit, lidx[a].v, lom[a].v)
    inv = (d.w :> 0) :/ (d.w + (d.w :<= 0))
    V = J(L, L, 0)
    mrk = J(d.nA, 1, 0)
    for (a = 1; a <= L; a++) {
        ea = rows_used[a]
        mrk[lidx[a].v] = lom[a].v :* inv[lidx[a].v]
        for (b = a; b <= L; b++) {
            eb = rows_used[b]
            if (C.r_g[ea] != C.r_g[eb]) continue
            sg = (C.r_t[ea] == C.r_t[eb]) + (C.r_s[ea] == C.r_s[eb]) -
                (C.r_t[ea] == C.r_s[eb]) - (C.r_s[ea] == C.r_t[eb])
            if (sg == 0) continue
            V[a, b] = sg * d.sumw * sum(mrk[lidx[b].v] :* lom[b].v)
            V[b, a] = V[a, b]
        }
        mrk[lidx[a].v] = J(rows(lidx[a].v), 1, 0)
    }
    return(0.5 :* (V + V'))
}


/* Exact-mask group of each comparison: target cohort, control pool and the exact set of
   analysis units observed at both dates (R mask keys of the within spanning forest). */
real colvector _cdid_wsp_key(struct _cdid_wsp_data scalar d, struct cdid_cat scalar C)
{
    struct _cdid_wsp_vec colvector mk
    real colvector idx_t, idx_c, u, key, ord, grp
    real matrix fp
    real scalar M, m, e, ng, j, hit
    M = rows(C.r_g)
    u = sin((1::d.nA) :* 1.618033988749895) :+ 2
    mk = _cdid_wsp_vec(M, 1)
    fp = J(M, 5, .)
    for (m = 1; m <= M; m++) {
        _cdid_wsp_units(d, C, m, idx_t, idx_c)
        mk[m].v = sort(idx_t \ idx_c, 1)
        fp[m, .] = (C.r_g[m], C.r_pool[m], rows(mk[m].v), sum(u[mk[m].v]), m)
    }
    ord = order(fp, (1, 2, 3, 4, 5))
    key = J(M, 1, 0)
    grp = J(0, 1, .)
    ng = 0
    for (e = 1; e <= M; e++) {
        m = ord[e]
        hit = 0
        if (e > 1) {
            if (all(fp[m, 1..4] :== fp[ord[e - 1], 1..4])) {
                for (j = 1; j <= rows(grp); j++) {
                    if (mk[grp[j]].v == mk[m].v) {
                        hit = key[grp[j]]
                        break
                    }
                }
            }
            else grp = J(0, 1, .)
        }
        if (!hit) {
            ng++
            hit = ng
            grp = grp \ m
        }
        key[m] = hit
    }
    return(key)
}

/* Covariance builder handed to the baseline IPW fit through its options. */
real matrix cdid_wsp_ipw_cov(struct cdid_cat scalar C, struct cdid_population scalar P, real matrix X,
    real colvector w, real colvector rows_used, real scalar eps, real scalar maxit)
{
    return(cdid_wsp_cov(_cdid_wsp_prepare(C, P, X, w), C, rows_used, "ipw", eps, maxit))
}

end
