*! Cross-pair chained DiD with selection-model observation weights (conditional missingness at random).
version 19.5
mata:

/* One observed-pair cell of the cross-pair catalogue: units of one cohort observed at both
   dates (analysis-unit and frame indices), selection weights and outcome changes. */
struct _cdid_xps_cell {
    real scalar ci, cs, ct, pc
    real colvector a, f, fw, dy
    real matrix ldt
    real colvector lda, gidx
    real scalar kpos
    real colvector rA, wA
    real scalar Mk, rraw, rc, ract, r
}

/* Outcome-regression block (one control cell) and its slope group. */
struct _cdid_xps_blk {
    real scalar cell, L, ws, alpha, my, grp, r
    real rowvector md
    real colvector trows
}

struct _cdid_xps_grp {
    real colvector beta
    real matrix M
    real colvector blocks
    real rowvector keep
}

struct cdid_xps_moments {
    real colvector mhat
    real matrix Phi, U
    real scalar n_trt_fits, n_or_groups, n_or_blocks
}

/* Reader and fitted selection models for the catalogue cells. */
struct cdid_xps_sel {
    struct cdid_sel_reader scalar rdr
    struct cdid_sel_fit_result scalar fit
    struct cdid_sel_up_map scalar um
    transmorphic pcl
}

struct cdid_xpn_pending {
    struct _cdid_xpn_rec colvector recs
    real scalar n, augmented, route
    string scalar method
    struct cdid_xps_sel scalar S
    real matrix linkcells
}

/* Rows of a recipe's keep mask that are nonzero (all rows when KP is empty). */
real colvector _cdid_xps_keep(real colvector rows, real matrix KP, real scalar r)
{
    if (rows(KP) == 0 | rows(rows) == 0) return(rows)
    return(select(rows, KP[rows, r] :!= 0))
}

real colvector _cdid_xps_cohort_rows(real colvector cohF, real colvector AF, real scalar g)
{
    if (g >= .) return(selectindex((cohF :>= .) :& (AF :== 1)))
    return(selectindex((cohF :== g) :& (AF :== 1)))
}

real scalar _cdid_xps_kpos(struct cdid_sel_reader scalar rdr, real scalar coh)
{
    real colvector h
    if (coh >= .) h = selectindex(rdr.cohorts :>= .)
    else h = selectindex(rdr.cohorts :== coh)
    if (rows(h) != 1) _cdid_ipw_fail("A cohort of the cross-pair catalogue is absent from the selection model.", 459)
    return(h[1])
}

/* Selection reader, cell map and fitted selection models over the target and control cells
   of every catalogue link. */
struct cdid_xps_sel scalar cdid_xps_sel_fit(struct cdid_public_bundle scalar B, struct cdid_xp_cat scalar C,
    real scalar eps, real scalar maxit, real scalar chunk, | real colvector links, real matrix Zr,
    real scalar icr, real colvector wr, real colvector keep)
{
    struct cdid_xps_sel scalar S
    struct cdid_sel_decl scalar decl
    struct cdid_population scalar pop
    real matrix specs
    real scalar N, univ_has_rho, int_col, l, i, hi
    real colvector lk
    string scalar univ
    pop = B.prepared.population
    N = pop.N
    univ = (strtrim(B.prepared.selection_universe) == "" ? "F" : strtrim(B.prepared.selection_universe))
    /* Optional recipe: the links it serves, its selection design and its weights. */
    hi = (args() >= 8 ? icr : B.prepared.sel_has_intercept)
    int_col = (hi ? 1 : 0)
    lk = (args() >= 6 ? links : (1::rows(C.g)))
    univ_has_rho = (B.prepared.sel_has_rho != 0 & B.prepared.sel_has_rho < . ? 1 :
        (univ != "F" | (rows(B.prepared.sel_rho) == N & any(B.prepared.sel_rho :!= 1))))
    decl = cdid_sel_norm_decl("conditional_mtar", 1, univ, univ_has_rho,
        (strtrim(B.prepared.selection_justification) != ""), B.prepared.selection_model,
        B.prepared.selection_pooling, B.prepared.selection_justification)
    if (args() >= 10) S.rdr = cdid_sel_reader_pop(pop, Zr, hi, int_col, wr, decl, B.prepared.sel_rho, chunk, keep)
    else if (args() >= 6) S.rdr = cdid_sel_reader_pop(pop, Zr, hi, int_col, wr, decl, B.prepared.sel_rho, chunk)
    else S.rdr = cdid_sel_reader_pop(pop, B.prepared.Z_sel, B.prepared.sel_has_intercept, int_col,
        B.prepared.weight.owner_weight, decl, B.prepared.sel_rho, chunk)
    specs = J(rows(lk), 5, .)
    for (i = 1; i <= rows(lk); i++) {
        l = lk[i]
        specs[i, ] = (_cdid_xps_kpos(S.rdr, C.coh[C.ci[l]]), C.cs[l], C.ct[l], C.cal[C.cs[l]], C.cal[C.ct[l]])
    }
    S.fit = cdid_sel_fit(S.rdr, cdid_sel_partition(S.rdr, specs), eps, maxit)
    S.um = cdid_sel_up_map_build(S.fit)
    S.pcl = asarray_create("real", 3)
    for (i = 1; i <= rows(S.fit.map.pair_cohorts); i++) {
        asarray(S.pcl, S.fit.map.pair_cohorts[i, 1..3], i)
    }
    return(S)
}

/* Selection-channel sensitivities of a cell: moment derivatives v (rows of the cell) are
   contracted with the cell's log-weight derivatives into G_theta and G_ak at column l. */
void _cdid_xps_sel_add(struct cdid_sel_moments scalar acc, struct _cdid_xps_cell scalar c,
    real colvector v, real scalar l, real scalar sel, real matrix Lr, real scalar ci)
{
    if (!sel) {
        if (c.ract) Lr[ci, l] = Lr[ci, l] - sum(v) / c.rc
        return
    }
    if (rows(c.gidx)) acc.G_theta[c.gidx, l] = acc.G_theta[c.gidx, l] + quadcross(c.ldt, v)
    acc.G_ak[c.kpos, l] = acc.G_ak[c.kpos, l] + quadcross(c.lda, v)
}

/* The same contraction into the moments of recipe r. */
void _cdid_xps_sel_add_r(struct cdid_sel_moments colvector acc, real scalar r, struct _cdid_xps_cell scalar c,
    real colvector v, real scalar l)
{
    if (rows(c.gidx)) acc[r].G_theta[c.gidx, l] = acc[r].G_theta[c.gidx, l] + quadcross(c.ldt, v)
    acc[r].G_ak[c.kpos, l] = acc[r].G_ak[c.kpos, l] + quadcross(c.lda, v)
}

/* Level loadings of a pair difference: + at the later date, - at the earlier date. */
void _cdid_xps_addU(real matrix U, real scalar nA, struct _cdid_xps_cell scalar c, real colvector x, real scalar l)
{
    U[(c.ct - 1) * nA :+ c.a, l] = U[(c.ct - 1) * nA :+ c.a, l] + x
    U[(c.cs - 1) * nA :+ c.a, l] = U[(c.cs - 1) * nA :+ c.a, l] - x
}

/* Within-unit covariance from level loadings (R within_cov_from_level_loadings): per target
   cohort, sum over dates of the scaled loading cross products, scale sqrt(sum w / w_i). */
real matrix cdid_xps_within_cov(real matrix U, real scalar nA, real colvector gvec, real colvector wA,
    | real matrix SC)
{
    real scalar L, T, tau, i
    real colvector sc, ug, cols_g
    real matrix V, Ut, S
    L = cols(U)
    T = rows(U) / nA
    if (args() < 5) sc = sqrt(sum(wA) :/ wA) :* (wA :> 0)
    V = J(L, L, 0)
    ug = uniqrows(gvec)
    for (i = 1; i <= rows(ug); i++) {
        /* Optional scales by target cohort (columns aligned with the sorted cohorts). */
        if (args() >= 5) sc = SC[, i]
        cols_g = selectindex(gvec :== ug[i])
        S = J(rows(cols_g), rows(cols_g), 0)
        for (tau = 1; tau <= T; tau++) {
            Ut = U[((tau - 1) * nA + 1)..(tau * nA), cols_g'] :* sc
            S = S + quadcross(Ut, Ut)
        }
        V[cols_g, cols_g] = 0.5 :* (S + S')
    }
    return(V)
}

void _cdid_xps_dated_fail()
{
    _cdid_ipw_fail("A covariate or sampling weight that varies over time is missing at a date read by the cross-pair comparisons (the last pre-treatment date of a target cohort for covariates and for fixweights(base_period), the first date for fixweights(first_period)); supply it at that date.", 459)
}

/* Rows of the design of recipe r: the columns of each recipe are stacked side by side
   (one recipe when covariates and weights are constant within units). */
real matrix _cdid_xps_xr(real matrix X, real scalar R, real scalar r, real colvector rows)
{
    real scalar k
    if (R == 1) return(X[rows, .])
    k = cols(X) / R
    if (k == 0) return(J(rows(rows), 0, .))
    return(X[rows, ((r - 1) * k + 1)..(r * k)])
}

/* Exact aliases of a rank-deficient outcome design: the dropped L-scaled columns must be
   reproduced by the kept ones (with a block constant when the design has a constant) on the
   fitting rows of every block and on the target rows where the regression is predicted. */
void cdid_xp_alias_check(real matrix Xc, real rowvector keep, real rowvector drop, real scalar icpt,
    struct _cdid_xps_blk colvector BK, real colvector blocks, real matrix Dfit, real colvector bfit,
    real matrix Xs, | real scalar R_in)
{
    real scalar R
    real matrix A, A0, Dt, H, Rf
    real rowvector tau, pv, c
    real colvector qcol
    real scalar j, b, err, tol, sc, q, rk
    R = (args() < 9 ? 1 : R_in)
    q = cols(Xs) / R
    hqrdp(Xc[, keep], H = ., tau = ., Rf = ., pv = .)
    rk = cols(keep)
    A0 = solveupper(Rf[1..rk, 1..rk], hqrdmultq1t(H, tau, Xc[, drop])[1..rk, ])
    A = A0
    A[pv', ] = A0
    sc = max((1, max(abs(Dfit))))
    err = 0
    for (j = 1; j <= rows(blocks); j++) {
        b = blocks[j]
        c = (icpt ? BK[b].md[drop] - BK[b].md[keep] * A : J(1, cols(drop), 0))
        qcol = selectindex(bfit :== j)
        err = max((err, max(abs((Dfit[qcol, drop] - Dfit[qcol, keep] * A) :- c))))
        if (rows(BK[b].trows)) {
            Dt = BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, BK[b].trows)
            sc = max((sc, max(abs(Dt))))
            err = max((err, max(abs((Dt[, drop] - Dt[, keep] * A) :- c))))
        }
    }
    tol = 64 * epsilon(1) * max((1, q)) * sc
    if (err >= . | err > tol)
        _cdid_ipw_fail("The outcome-regression design is rank deficient on the fitting rows of a control cohort without an exact alias that also holds on the target cohort; the counterfactual mean is not identified without extrapolation. Remove the redundant covariate.", 459)
}

/* Pooled or separate weighted least squares of the outcome change on L-scaled slopes with
   block intercepts (none without a constant). Exactly aliased columns are dropped. */
void _cdid_xps_or_fit(struct _cdid_xps_blk colvector BK, struct _cdid_xps_cell colvector CE,
    struct _cdid_xps_grp colvector GR, real matrix Xs, | real scalar icpt_in, real scalar R_in)
{
    real scalar g, j, b, q, rk, mx, icpt, k, R
    real matrix Xc, H, Rf, Dfit
    real colvector yc, qy, sw, bfit
    real rowvector tau, pv, keep, drop
    icpt = (args() < 5 ? 1 : icpt_in)
    R = (args() < 6 ? 1 : R_in)
    q = cols(Xs) / R
    for (g = 1; g <= rows(GR); g++) {
        if (q == 0) {
            GR[g].beta = J(0, 1, .)
            GR[g].M = J(0, 0, .)
            GR[g].keep = J(1, 0, .)
            for (j = 1; j <= rows(GR[g].blocks); j++) {
                b = GR[g].blocks[j]
                BK[b].alpha = BK[b].my
            }
            continue
        }
        Xc = J(0, q, .)
        yc = J(0, 1, .)
        Dfit = J(0, q, .)
        bfit = J(0, 1, .)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            sw = sqrt(CE[BK[b].cell].fw)
            Dfit = Dfit \ (BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, CE[BK[b].cell].f))
            bfit = bfit \ J(rows(sw), 1, j)
            Xc = Xc \ ((BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, CE[BK[b].cell].f)) :- BK[b].md) :* sw
            yc = yc \ (CE[BK[b].cell].dy :- BK[b].my) :* sw
        }
        hqrdp(Xc, H = ., tau = ., Rf = ., pv = .)
        mx = max(abs(diagonal(Rf)))
        rk = sum(abs(diagonal(Rf)) :> 1e-7 * mx)
        if (mx <= 0 | rk < 1) _cdid_ipw_fail("The outcome-regression design has no numerical rank on the observed pairs of a control cohort.", 459)
        keep = sort(pv[1..rk]', 1)'
        if (rk < q) {
            drop = sort(pv[(rk + 1)..q]', 1)'
            cdid_xp_alias_check(Xc, keep, drop, icpt, BK, GR[g].blocks, Dfit, bfit, Xs, R)
        }
        k = cols(keep)
        hqrdp(Xc[, keep], H = ., tau = ., Rf = ., pv = .)
        qy = hqrdmultq1t(H, tau, yc)
        GR[g].beta = J(q, 1, 0)
        GR[g].beta[keep[pv]'] = solveupper(Rf[1..k, 1..k], qy[1..k])
        GR[g].M = quadcross(Xc[, keep], Xc[, keep])
        GR[g].keep = keep
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            BK[b].alpha = (icpt ? BK[b].my - BK[b].md * GR[g].beta : 0)
        }
    }
}

/* Link moments and frame-unit influence functions (R evaluate_weighted_regimes_core with
   selection_est_method ipw on the cross-pair descriptors). Target links: selection-weighted
   mean outcome change; control links: control change reweighted toward the target cohort
   (selection-weighted on its universe rows) by IPW, OR or DR. */
struct cdid_xps_moments scalar cdid_xps_eval(struct cdid_xps_sel colvector S, real scalar sel,
    struct cdid_xp_cat scalar C, real matrix Y, real colvector aF, real matrix Xt, real matrix Xs,
    string scalar method, string scalar pooling, real scalar eps, real scalar maxit, real scalar want_U,
    real colvector wF, real colvector cohF, real colvector AF, | real scalar icpt_in,
    real scalar R_in, real colvector lrec_in, real matrix W_in, real matrix KP_in)
{
    struct cdid_xps_moments scalar out
    struct cdid_sel_moments colvector acc
    struct cdid_sel_moments scalar acc0
    struct _cdid_xps_cell colvector CE
    struct _cdid_xps_blk colvector BK
    struct _cdid_xps_grp colvector GR
    struct cdid_sel_pred scalar pr
    struct cdid_sel_log_deriv scalar ld
    transmorphic CK, BKK, TK
    real scalar N, L, q, l, c, b, g, kg, ag, Sg, Bv, Cv, Bol, Sm, isor, isdr, need_or, need_trt, ncel, nblk, ng, nc, j
    real colvector lc, lb, rows_g, fw_g, mh_g, beta_t, prob, odds_c, act_c, wm, om, eff, mh_c, h, sc, v, e_s, pos
    real colvector y_all, w_all, rows_all, a
    real matrix d_g, d_c, X_all, info, Lam_a, Pm, Wint, nu, dt, ldum, prb
    real rowvector dbp, dbm, Vg, lam
    real scalar key_g, nA, T
    real matrix Lr
    real scalar icpt, nreq, tix, kc
    real rowvector kp
    real scalar R, rl
    real colvector lrec
    real matrix W, KP
    struct _cdid_xpn_rec colvector NR
    struct _cdid_xpn_rec scalar NRr
    struct cdid_xpn_pending scalar NP
    transmorphic TKR
    icpt = (args() < 16 ? 1 : icpt_in)
    N = rows(wF)
    L = rows(C.g)
    /* Recipes: covariates read at the reference date of each link's target cohort and
       weights at their anchor date; one recipe when both are constant within units. */
    R = (args() < 17 ? 1 : R_in)
    lrec = (args() < 18 ? J(L, 1, 1) : lrec_in)
    W = (args() < 19 ? wF : W_in)
    /* Units each recipe reads under restricted reference eligibility (empty: all). */
    KP = (args() < 20 ? J(0, 0, .) : KP_in)
    if (sel & rows(S) != R) _cdid_ipw_fail("Each design block needs its own selection models.", 498)
    q = cols(Xs) / R
    isor = (method == "reg")
    isdr = (method == "dr")
    need_or = (isor | isdr)
    need_trt = !isor
    out.mhat = J(L, 1, .)
    out.Phi = J(N, L, 0)
    out.n_trt_fits = 0
    nA = rows(Y)
    T = rows(C.cal)
    out.U = (want_U ? J(T * nA, L, 0) : J(0, 0, .))
    if (sel) {
        acc = cdid_sel_moments(R, 1)
        for (j = 1; j <= R; j++) acc[j] = cdid_sel_moments_init(S[j].fit, L)
    }
    /* Cells. */
    CK = asarray_create("real", 4)
    lc = J(L, 1, 0)
    ncel = 0
    for (l = 1; l <= L; l++) {
        if (!asarray_contains(CK, (C.ci[l], C.cs[l], C.ct[l], lrec[l]))) {
            ncel++
            asarray(CK, (C.ci[l], C.cs[l], C.ct[l], lrec[l]), ncel)
        }
        lc[l] = asarray(CK, (C.ci[l], C.cs[l], C.ct[l], lrec[l]))
    }
    CE = _cdid_xps_cell(ncel, 1)
    for (l = 1; l <= L; l++) {
        c = lc[l]
        if (rows(CE[c].f)) continue
        CE[c].ci = C.ci[l]
        CE[c].cs = C.cs[l]
        CE[c].ct = C.ct[l]
        CE[c].r = lrec[l]
        CE[c].a = *C.obs[C.oi[l], C.mid[l]]
        CE[c].f = aF[CE[c].a]
        CE[c].dy = Y[CE[c].a, C.ct[l]] - Y[CE[c].a, C.cs[l]]
        if (sel) {
            rl = CE[c].r
            CE[c].pc = asarray(S[rl].pcl, (_cdid_xps_kpos(S[rl].rdr, C.coh[C.ci[l]]), C.cs[l], C.ct[l]))
            pr = cdid_sel_pred_pair(S[rl].fit, S[rl].rdr, CE[c].pc, CE[c].f)
            CE[c].fw = pr.sel_weight
            ld = cdid_sel_log_deriv_eval(S[rl].fit, S[rl].rdr, S[rl].um, CE[c].pc, CE[c].f)
            CE[c].ldt = ld.d_theta_local
            CE[c].lda = ld.d_ak
            CE[c].kpos = ld.cohort_pos
            CE[c].gidx = _cdid_sel_pc_global_idx(S[rl].fit, S[rl].um, CE[c].pc)
            if (R > 1 | args() >= 19) {
                if (missing(CE[c].fw)) _cdid_xps_dated_fail()
            }
        }
        else {
            /* Pair rate of the cohort's analysis units observed at both dates. */
            CE[c].rA = _cdid_xps_keep(_cdid_xps_cohort_rows(cohF, AF, C.coh[C.ci[l]]), KP, CE[c].r)
            CE[c].wA = W[CE[c].rA, CE[c].r]
            if (R > 1 | args() >= 19) {
                if (missing(CE[c].wA)) _cdid_xps_dated_fail()
                if (C.role[l] == 1 & cols(Xt)) {
                    if (missing(_cdid_xps_xr(Xt, R, CE[c].r, CE[c].f))) _cdid_xps_dated_fail()
                }
            }
            CE[c].Mk = sum(CE[c].wA)
            CE[c].rraw = (rows(CE[c].f) < rows(CE[c].rA) ? sum(W[CE[c].f, CE[c].r]) / CE[c].Mk : 1)
            CE[c].rc = max((CE[c].rraw, 1e-6))
            CE[c].ract = (rows(CE[c].f) < rows(CE[c].rA)) & (CE[c].rraw > 1e-6)
            CE[c].fw = W[CE[c].f, CE[c].r] :/ CE[c].rc
        }
    }
    Lr = (sel ? J(0, 0, .) : J(ncel, L, 0))
    ncel = rows(CE)
    /* Outcome-regression blocks, one per control cell, in control-link order. */
    lb = J(L, 1, 0)
    if (need_or) {
        if (!icpt & q == 0) _cdid_ipw_fail("The outcome design is empty; add covariates or a constant.", 459)
        BKK = asarray_create("real", 1)
        nblk = 0
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            if (!asarray_contains(BKK, lc[l])) {
                nblk++
                asarray(BKK, lc[l], nblk)
            }
            lb[l] = asarray(BKK, lc[l])
        }
        BK = _cdid_xps_blk(nblk, 1)
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            b = lb[l]
            if (BK[b].cell < .) continue
            BK[b].cell = lc[l]
            BK[b].r = lrec[l]
            BK[b].L = C.ct[l] - C.cs[l]
            BK[b].ws = sum(CE[lc[l]].fw)
            BK[b].my = (icpt ? sum(CE[lc[l]].fw :* CE[lc[l]].dy) / BK[b].ws : 0)
            if (q & (R > 1 | args() >= 19)) {
                if (missing(_cdid_xps_xr(Xs, R, BK[b].r, CE[lc[l]].f))) _cdid_xps_dated_fail()
            }
            BK[b].md = (q ? (icpt ? colsum(CE[lc[l]].fw :* (BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, CE[lc[l]].f))) :/ BK[b].ws : J(1, q, 0)) : J(1, 0, .))
            BK[b].trows = J(0, 1, .)
        }
        /* Target rows where each block's regression is predicted. */
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            if (sel) rows_g = S[lrec[l]].rdr.coh_u_idx[_cdid_xps_kpos(S[lrec[l]].rdr, C.g[l])].v
            else rows_g = _cdid_xps_keep(_cdid_xps_cohort_rows(cohF, AF, C.g[l]), KP, lrec[l])
            BK[lb[l]].trows = uniqrows(BK[lb[l]].trows \ rows_g)
        }
        if (pooling == "pooled") {
            GR = _cdid_xps_grp(1, 1)
            GR[1].blocks = (1::nblk)
            for (b = 1; b <= nblk; b++) BK[b].grp = 1
        }
        else {
            GR = _cdid_xps_grp(nblk, 1)
            for (b = 1; b <= nblk; b++) {
                GR[b].blocks = b
                BK[b].grp = b
            }
        }
        if (nblk) _cdid_xps_or_fit(BK, CE, GR, Xs, icpt, R)
        ng = rows(GR)
        Lam_a = J(nblk, L, 0)
        Pm = J(ng * q, L, 0)
    }
    else {
        nblk = 0
        ng = 0
    }
    out.n_or_groups = ng
    out.n_or_blocks = nblk
    TK = asarray_create("real", 3)
    TKR = asarray_create("real", 3)
    nreq = cdid_nuis_requested()
    NR = _cdid_xpn_rec(0, 1)
    for (l = 1; l <= L; l++) {
        c = lc[l]
        rl = lrec[l]
        if (C.role[l] == 1) {
            Sg = sum(CE[c].fw)
            out.mhat[l] = sum(CE[c].fw :* CE[c].dy) / Sg
            v = (CE[c].fw :/ Sg) :* (CE[c].dy :- out.mhat[l])
            out.Phi[CE[c].f, l] = out.Phi[CE[c].f, l] + N :* v
            if (sel) _cdid_xps_sel_add_r(acc, rl, CE[c], v, l)
            if (want_U) _cdid_xps_addU(out.U, nA, CE[c], CE[c].fw :/ Sg, l)
            continue
        }
        /* Target domain: universe rows of the target cohort with weights w a_g / rho. */
        if (sel) {
            kg = _cdid_xps_kpos(S[rl].rdr, C.g[l])
            rows_g = S[rl].rdr.coh_u_idx[kg].v
            ag = S[rl].fit.a_U[kg]
            if (rows(rows_g) == 0 | ag <= 0 | ag >= .) _cdid_ipw_fail("A target cohort has no repeat observers in the selection universe.", 459)
            fw_g = S[rl].rdr.w[rows_g] :* ag :/ S[rl].rdr.rho[rows_g]
            if (R > 1 | args() >= 19) {
                if (missing(fw_g)) _cdid_xps_dated_fail()
                if (cols(Xt)) {
                    if (missing(_cdid_xps_xr(Xt, R, rl, rows_g))) _cdid_xps_dated_fail()
                }
            }
        }
        else {
            kg = C.gi[l]
            rows_g = _cdid_xps_keep(_cdid_xps_cohort_rows(cohF, AF, C.g[l]), KP, rl)
            if (rows(rows_g) == 0) _cdid_ipw_fail("A target cohort has no analysis units.", 459)
            fw_g = W[rows_g, rl]
            if (R > 1 | args() >= 19) {
                if (missing(fw_g)) _cdid_xps_dated_fail()
                if (cols(Xt)) {
                    if (missing(_cdid_xps_xr(Xt, R, rl, rows_g))) _cdid_xps_dated_fail()
                }
            }
        }
        Sg = sum(fw_g)
        Bv = 0
        Cv = 0
        Bol = 0
        if (need_or) {
            b = lb[l]
            g = BK[b].grp
            if (q) {
                d_g = BK[b].L :* _cdid_xps_xr(Xs, R, rl, rows_g)
                mh_g = BK[b].alpha :+ d_g * GR[g].beta
                dbp = colsum(fw_g :* d_g) :/ Sg
            }
            else {
                mh_g = J(rows(rows_g), 1, BK[b].alpha)
                dbp = J(1, 0, .)
            }
            Bv = sum(fw_g :* mh_g) / Sg
            out.Phi[rows_g, l] = out.Phi[rows_g, l] + (N :* fw_g :/ Sg) :* (mh_g :- Bv)
        }
        if (need_trt) {
            nc = rows(CE[c].f)
            X_all = _cdid_xps_xr(Xt, R, rl, rows_g) \ _cdid_xps_xr(Xt, R, rl, CE[c].f)
            y_all = J(rows(rows_g), 1, 1) \ J(nc, 1, 0)
            w_all = fw_g \ CE[c].fw
            rows_all = rows_g \ CE[c].f
            /* One logit per target and control cell; without selection, the cells of a
               cohort observed completely share one logit. */
            kc = ((!sel & rows(CE[c].f) == rows(CE[c].rA)) ? -CE[c].ci : c)
            if (asarray_contains(TK, (kg, kc, rl))) {
                beta_t = asarray(TK, (kg, kc, rl))
                if (nreq) {
                    tix = asarray(TKR, (kg, kc, rl))
                    NR[tix].links = NR[tix].links \ l
                }
            }
            else {
                _cdid_sel_link_trt_logit(X_all, y_all, w_all, eps, maxit, beta_t = ., ldum = .)
                asarray(TK, (kg, kc, rl), beta_t)
                out.n_trt_fits = out.n_trt_fits + 1
            }
            prob = _cdid_ipw_logit_values(X_all * beta_t)
            info = quadcross(X_all, w_all :* prob[, 1] :* (1 :- prob[, 1]), X_all)
            if (nreq & !asarray_contains(TKR, (kg, kc, rl))) {
                NRr = cdid_xpn_wtd_trt(S[(sel ? rl : 1)], sel, CE[c], X_all, y_all, w_all, rows_all, rows(rows_g), beta_t, prob[, 1],
                    info, kg, (sel ? ag : 1), N, C.g[l], C.coh[C.ci[l]], C.cal[C.cs[l]], C.cal[C.ct[l]], (kc < 0), W[, rl], l)
                NR = NR \ NRr
                asarray(TKR, (kg, kc, rl), rows(NR))
            }
            pos = (rows(rows_g) + 1)::rows(X_all)
            odds_c = prob[pos, 2]
            act_c = prob[pos, 3]
            prob = prob[, 1]
            wm = CE[c].fw :* odds_c
            Sm = sum(wm)
            om = wm :/ Sm
            Cv = sum(om :* CE[c].dy)
            if (isdr) {
                if (q) {
                    d_c = BK[b].L :* _cdid_xps_xr(Xs, R, rl, CE[c].f)
                    mh_c = BK[b].alpha :+ d_c * GR[g].beta
                    dbm = colsum(om :* d_c)
                }
                else mh_c = J(nc, 1, BK[b].alpha)
                Bol = sum(om :* mh_c)
                eff = (CE[c].dy :- mh_c) :- (Cv - Bol)
            }
            else eff = CE[c].dy :- Cv
            out.Phi[CE[c].f, l] = out.Phi[CE[c].f, l] + (N :* wm :/ Sm) :* eff
            if (want_U) _cdid_xps_addU(out.U, nA, CE[c], om, l)
            Vg = colsum((wm :* act_c :* eff) :* _cdid_xps_xr(Xt, R, rl, CE[c].f)) :/ Sm
            h = N :* cholsolve(info, Vg')
            if (missing(h)) _cdid_ipw_fail("The treatment logit information is singular.", 459)
            sc = w_all :* (y_all :- prob) :* (X_all * h)
            out.Phi[rows_all, l] = out.Phi[rows_all, l] + sc
            if (sel) acc[rl].G_ak[kg, l] = acc[rl].G_ak[kg, l] + sum(sc[1..rows(rows_g)]) / (ag * N)
            v = (sc[pos] :+ (N :* wm :/ Sm) :* eff) :/ N
            if (sel) _cdid_xps_sel_add_r(acc, rl, CE[c], v, l)
            else _cdid_xps_sel_add(acc0, CE[c], v, l, sel, Lr, c)
        }
        if (isor) {
            out.mhat[l] = Bv
            if (icpt) Lam_a[b, l] = Lam_a[b, l] + 1
            if (q) Pm[((g - 1) * q + 1)..(g * q), l] = Pm[((g - 1) * q + 1)..(g * q), l] + dbp'
        }
        else if (isdr) {
            out.mhat[l] = Bv + Cv - Bol
            if (q) Pm[((g - 1) * q + 1)..(g * q), l] = Pm[((g - 1) * q + 1)..(g * q), l] + (dbp - dbm)'
        }
        else out.mhat[l] = Cv
    }
    /* Outcome-regression adjoints on the control cells. */
    if (nblk) {
        for (b = 1; b <= nblk; b++) {
            if (q & any(Lam_a[b, ] :!= 0)) {
                g = BK[b].grp
                Pm[((g - 1) * q + 1)..(g * q), ] = Pm[((g - 1) * q + 1)..(g * q), ] - BK[b].md' * Lam_a[b, ]
            }
        }
        for (b = 1; b <= nblk; b++) {
            g = BK[b].grp
            c = BK[b].cell
            if (q) e_s = CE[c].dy :- (BK[b].alpha :+ (BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, CE[c].f)) * GR[g].beta)
            else e_s = CE[c].dy :- BK[b].alpha
            lam = Lam_a[b, ]
            pos = selectindex(lam' :!= 0)
            if (rows(pos)) {
                Wint = (CE[c].fw :/ BK[b].ws) * lam[pos']
                for (j = 1; j <= rows(pos); j++) {
                    v = e_s :* Wint[, j]
                    out.Phi[CE[c].f, pos[j]] = out.Phi[CE[c].f, pos[j]] + N :* v
                    if (sel) _cdid_xps_sel_add_r(acc, CE[c].r, CE[c], v, pos[j])
                    else _cdid_xps_sel_add(acc0, CE[c], v, pos[j], sel, Lr, c)
                    if (want_U) _cdid_xps_addU(out.U, nA, CE[c], Wint[, j], pos[j])
                }
            }
            if (q) {
                kp = GR[g].keep
                nu = Pm[((g - 1) * q) :+ kp, ]
                pos = selectindex(colsum(abs(nu))' :> 0)
                if (rows(pos)) {
                    nu = N :* cholsolve(GR[g].M, nu[, pos'])
                    dt = (BK[b].L :* _cdid_xps_xr(Xs, R, BK[b].r, CE[c].f)[, kp]) :- BK[b].md[kp]
                    Wint = ((CE[c].fw :/ N) :* dt) * nu
                    for (j = 1; j <= rows(pos); j++) {
                        v = e_s :* Wint[, j]
                        out.Phi[CE[c].f, pos[j]] = out.Phi[CE[c].f, pos[j]] + N :* v
                        if (sel) _cdid_xps_sel_add_r(acc, CE[c].r, CE[c], v, pos[j])
                    else _cdid_xps_sel_add(acc0, CE[c], v, pos[j], sel, Lr, c)
                        if (want_U) _cdid_xps_addU(out.U, nA, CE[c], Wint[, j], pos[j])
                    }
                }
            }
        }
    }
    /* Reports of the outcome regressions and the capture of this call. */
    if (nreq) {
        for (g = 1; g <= ng; g++) {
            NRr = cdid_xpn_wtd_or(S[1], sel, BK, CE, GR[g], Xs, icpt, N, W[, 1], selectindex(_cdid_xpn_in(lb, GR[g].blocks) :& (lb :> 0)), C)
            NR = NR \ NRr
        }
        NP.recs = NR
        NP.n = rows(NR)
        NP.augmented = 0
        NP.route = (sel ? 2 : 1)
        NP.method = method
        if (sel) NP.S = S[1]
        NP.linkcells = J(L, 3, .)
        for (l = 1; l <= L; l++) NP.linkcells[l, ] = (C.coh[C.ci[l]], C.cs[l], C.ct[l])
        _cdid_xpn_put(NP)
    }
    /* Selection-model and repeat-rate channels. */
    if (sel) {
        for (j = 1; j <= R; j++) {
            if (max(abs(acc[j].G_theta)) > 0 | max(abs(acc[j].G_ak)) > 0)
                out.Phi = out.Phi + cdid_sel_eval_loading(S[j].fit, S[j].rdr, acc[j], (1::N))
        }
    }
    else {
        /* Pair-rate channels: r = W_obs / M_k on the cohort's analysis units. */
        for (c = 1; c <= ncel; c++) {
            if (!CE[c].ract) continue
            pos = selectindex(Lr[c, ]' :!= 0)
            if (rows(pos) == 0) continue
            out.Phi[CE[c].rA, pos'] = out.Phi[CE[c].rA, pos'] - ((CE[c].wA :* (N / CE[c].Mk) :* CE[c].rraw) * Lr[c, pos'])
            out.Phi[CE[c].f, pos'] = out.Phi[CE[c].f, pos'] + ((W[CE[c].f, CE[c].r] :* (N / CE[c].Mk)) * Lr[c, pos'])
        }
    }
    return(out)
}

/* Post-fit exact-mask spanning forest (R verify_post_fit_forest_reduction): within each
   (target, cohort, mask) group a chord is dropped only when its fitted moment and its
   unit influence function telescope along the tree path (relative tolerance 1e-9);
   kept links ascending, metric G = I + E E' on the verified chords. */
void cdid_xps_forest(struct cdid_xp_cat scalar C, real colvector mhat, real matrix Phi,
    real matrix U, real scalar nA, real colvector kept, real matrix G, real scalar maxdiff)
{
    real scalar M, ng, e, gr, m, nv, k, d, dd, u, v, ru, rv, x, y, i1, i2, nk, tol, tau, T
    real rowvector hasany
    real matrix Utr, Uca, Upr
    real colvector key, ord, idx, us, ut, dates, par, istree, tpos, dpos, depth, pv, ple, psg, queue
    real colvector tree_ids, cand_ids, mt, mc, mp, mscale, merr, ok, surv, nz, keptall, posmap, pp, kg
    real matrix gk, E, adj, Ptr, Pca, Ppr, Bg
    real rowvector pscale, perr
    pointer(real colvector) colvector bpos
    pointer(real matrix) colvector bmat
    tol = 1e-9
    M = rows(C.g)
    maxdiff = 0
    gk = (C.gi, C.ci, C.mid)
    key = J(M, 1, 0)
    ord = order((gk, (1::M)), (1, 2, 3, 4))
    ng = 0
    for (e = 1; e <= M; e++) {
        if (e == 1) ng = 1
        else if (any(gk[ord[e], ] :!= gk[ord[e - 1], ])) ng++
        key[ord[e]] = ng
    }
    keptall = J(0, 1, .)
    bpos = J(0, 1, NULL)
    bmat = J(0, 1, NULL)
    for (gr = 1; gr <= ng; gr++) {
        idx = selectindex(key :== gr)
        m = rows(idx)
        if (m == 1) {
            keptall = keptall \ idx
            continue
        }
        idx = idx[order((abs(C.cal[C.ct[idx]] - C.cal[C.cs[idx]]), C.cal[C.cs[idx]], C.cal[C.ct[idx]], idx), (1, 2, 3, 4))]
        dates = uniqrows(C.cs[idx] \ C.ct[idx])
        nv = rows(dates)
        us = J(m, 1, .)
        ut = J(m, 1, .)
        for (e = 1; e <= m; e++) {
            us[e] = selectindex(dates :== C.cs[idx[e]])
            ut[e] = selectindex(dates :== C.ct[idx[e]])
        }
        par = (1::nv)
        istree = J(m, 1, 0)
        for (e = 1; e <= m; e++) {
            ru = us[e]
            while (par[ru] != ru) {
                par[ru] = par[par[ru]]
                ru = par[ru]
            }
            rv = ut[e]
            while (par[rv] != rv) {
                par[rv] = par[par[rv]]
                rv = par[rv]
            }
            if (ru != rv) {
                par[rv] = ru
                istree[e] = 1
            }
        }
        tpos = selectindex(istree)
        dpos = selectindex(!istree)
        tree_ids = idx[tpos]
        k = rows(tpos)
        if (rows(dpos) == 0) {
            keptall = keptall \ tree_ids
            continue
        }
        adj = J(0, 4, .)
        for (e = 1; e <= k; e++) adj = adj \ (us[tpos[e]], ut[tpos[e]], e, 1) \ (ut[tpos[e]], us[tpos[e]], e, -1)
        depth = J(nv, 1, .)
        pv = J(nv, 1, 0)
        ple = J(nv, 1, 0)
        psg = J(nv, 1, 0)
        for (u = 1; u <= nv; u++) {
            if (depth[u] < . | !any(adj[, 1] :== u)) continue
            depth[u] = 0
            queue = u
            x = 1
            while (x <= rows(queue)) {
                v = queue[x]
                x++
                for (y = 1; y <= rows(adj); y++) {
                    if (adj[y, 1] != v) continue
                    if (depth[adj[y, 2]] < .) continue
                    depth[adj[y, 2]] = depth[v] + 1
                    pv[adj[y, 2]] = v
                    ple[adj[y, 2]] = adj[y, 3]
                    psg[adj[y, 2]] = adj[y, 4]
                    queue = queue \ adj[y, 2]
                }
            }
        }
        d = rows(dpos)
        E = J(k, d, 0)
        for (dd = 1; dd <= d; dd++) {
            u = us[dpos[dd]]
            v = ut[dpos[dd]]
            while (depth[u] > depth[v]) {
                E[ple[u], dd] = E[ple[u], dd] - psg[u]
                u = pv[u]
            }
            while (depth[v] > depth[u]) {
                E[ple[v], dd] = E[ple[v], dd] + psg[v]
                v = pv[v]
            }
            while (u != v) {
                E[ple[u], dd] = E[ple[u], dd] - psg[u]
                u = pv[u]
                E[ple[v], dd] = E[ple[v], dd] + psg[v]
                v = pv[v]
            }
        }
        cand_ids = idx[dpos]
        mt = mhat[tree_ids]
        mc = mhat[cand_ids]
        mp = E' * mt
        mscale = rowmax((J(d, 1, 1), abs(mc), colmax(abs(E) :* abs(mt))'))
        merr = abs(mc - mp)
        ok = (merr :< .) :& (merr :<= tol :* mscale)
        if (any(ok) & rows(Phi)) {
            surv = selectindex(ok)
            Ptr = Phi[, tree_ids]
            Pca = Phi[, cand_ids[surv]]
            nz = selectindex((rowsum(Ptr :!= 0) :> 0) :| (rowsum(Pca :!= 0) :> 0))
            if (rows(nz)) {
                Ppr = Ptr[nz, ] * E[, surv]
                Pca = Pca[nz, ]
                pscale = colmax((J(1, cols(Pca), 1) \ colmax(abs(Pca)) \ colmax(abs(Ppr))))
                perr = colmax(abs(Pca - Ppr))
                ok[surv] = ((perr :< .) :& (perr :<= tol :* pscale))'
            }
        }
        if (any(ok) & rows(U)) {
            T = rows(U) / nA
            for (tau = 1; tau <= T; tau++) {
                surv = selectindex(ok)
                if (length(surv) == 0) break
                Utr = U[((tau - 1) * nA + 1)..(tau * nA), tree_ids]
                Uca = U[((tau - 1) * nA + 1)..(tau * nA), cand_ids[surv]]
                nz = selectindex((rowsum(Utr :!= 0) :> 0) :| (rowsum(Uca :!= 0) :> 0))
                if (rows(nz) == 0) continue
                Upr = Utr[nz, ] * E[, surv]
                Uca = Uca[nz, ]
                hasany = (colsum(Uca :!= 0) :> 0) :| (colsum(Upr :!= 0) :> 0)
                pscale = colmax((J(1, cols(Uca), 1) \ colmax(abs(Uca)) \ colmax(abs(Upr))))
                perr = colmax(abs(Uca - Upr))
                ok[surv] = ok[surv] :& !(hasany :& ((perr :>= .) :| (perr :> tol :* pscale)))'
            }
        }
        if (any(ok)) maxdiff = max((maxdiff, max(select(merr, ok))))
        /* Chords that failed verification stay; a single candidate makes select() ambiguous
           in orientation, so the kept chords are indexed explicitly. */
        kg = tree_ids
        if (any(!ok)) kg = kg \ cand_ids[selectindex(!ok)]
        nk = rows(kg)
        Bg = I(nk)
        if (any(ok)) {
            Bg[1..k, 1..k] = I(k) + E[, selectindex(ok)] * E[, selectindex(ok)]'
            Bg[1..k, 1..k] = 0.5 :* (Bg[1..k, 1..k] + Bg[1..k, 1..k]')
        }
        keptall = keptall \ kg
        bpos = bpos \ &(kg[1..nk])
        bmat = bmat \ &(Bg[., .])
    }
    kept = sort(keptall, 1)
    posmap = J(M, 1, 0)
    posmap[kept] = (1::rows(kept))
    G = I(rows(kept))
    for (e = 1; e <= rows(bpos); e++) {
        pp = posmap[*bpos[e]]
        Bg = *bmat[e]
        for (i1 = 1; i1 <= rows(pp); i1++) {
            for (i2 = 1; i2 <= rows(pp); i2++) G[pp[i1], pp[i2]] = Bg[i1, i2]
        }
    }
    G = 0.5 :* (G + G')
}

end
