*! Cross-pair chained DiD with imputed outcome changes (selmethod(imputation)) and with the
*! augmented observation correction (selmethod(dr)) when covariates or weights vary over time.
version 19.5
mata:

/* Exactly aliased slopes of a counterfactual outcome design. Keeps the leading columns of the
   pivoted QR decomposition of the centered weighted design Xc and returns in A the coefficients
   that reproduce the dropped columns from the kept ones; drop is empty at full rank. */
void cdid_or_alias_keep(real matrix Xc, real scalar rank_tol, string scalar where,
    real rowvector keep, real rowvector drop, real matrix A)
{
    real matrix H, Rf, A0
    real rowvector tau, pv
    real scalar q, mx, rk
    q = cols(Xc)
    hqrdp(Xc, H = ., tau = ., Rf = ., pv = .)
    mx = max(abs(diagonal(Rf)))
    rk = sum(abs(diagonal(Rf)) :> rank_tol * mx)
    if (mx <= 0 | rk < 1 | mx >= .) _cdid_ipw_fail("The counterfactual outcome-regression design has no numerical rank on the frame rows of " + where + ".", 459)
    keep = sort(pv[1..rk]', 1)'
    if (rk == q) {
        drop = J(1, 0, .)
        A = J(rk, 0, .)
        return
    }
    drop = sort(pv[(rk + 1)..q]', 1)'
    hqrdp(Xc[, keep], H = ., tau = ., Rf = ., pv = .)
    A0 = solveupper(Rf[1..rk, 1..rk], hqrdmultq1t(H, tau, Xc[, drop])[1..rk, ])
    A = A0
    A[pv', ] = A0
}

/* Largest violation of an exact alias D[, drop] = D[, keep] A + c on the L-scaled rows D, and
   the largest absolute entry of D (scale of the tolerance), accumulated over calls. */
void cdid_or_alias_gap(real matrix D, real rowvector keep, real rowvector drop, real matrix A,
    real rowvector c, real scalar err, real scalar sc)
{
    if (rows(D) == 0) return
    sc = max((sc, max(abs(D))))
    err = max((err, max(abs((D[, drop] - D[, keep] * A) :- c))))
}

/* Weighted least squares of the outcome responses on the kept columns of Xc; dropped slopes are
   zero. */
void cdid_or_alias_solve(real matrix Xc, real colvector yc, real rowvector keep, real scalar q,
    real colvector beta, real matrix M)
{
    real matrix H, Rf
    real rowvector tau, pv
    real colvector qy
    real scalar k
    k = cols(keep)
    hqrdp(Xc[, keep], H = ., tau = ., Rf = ., pv = .)
    qy = hqrdmultq1t(H, tau, yc)
    beta = J(q, 1, 0)
    beta[keep[pv]'] = solveupper(Rf[1..k, 1..k], qy[1..k])
    M = quadcross(Xc[, keep], Xc[, keep])
}

/* Counterfactual outcome regressions of the imputed responses on the L-scaled slopes of each
   control cell's recipe over its cohort's frame rows (weights w a_k), with block intercepts;
   separate or pooled slopes (one slope shared by the blocks of every recipe when pooled).
   Exactly aliased columns are dropped when the alias, with one constant per block when the
   design has a constant, holds on every frame row of the cohorts each recipe reads (PR), where
   the regression is fitted and predicted; any other rank deficiency is refused. */
void _cdid_xprd_or_fit(struct _cdid_xpr_blk colvector BK, struct _cdid_xpr_cell colvector CE,
    real colvector crec, struct _cdid_xpr_grp colvector GR, real matrix XS, real scalar R,
    real scalar icpt, real scalar rank_tol, struct _cdid_xpr_mat colvector PR)
{
    real scalar g, j, b, q, c, err, sc
    real matrix Xc, D, A
    real colvector yc, sw, my
    real rowvector keep, drop, cst
    q = cols(XS) / R
    for (g = 1; g <= rows(GR); g++) {
        my = J(rows(GR[g].blocks), 1, 0)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            BK[b].wsum = sum(CE[c].w_syn)
            if (icpt) {
                my[j] = sum(CE[c].w_syn :* CE[c].Z) / BK[b].wsum
                BK[b].mean_d = (q ? colsum(CE[c].w_syn :* (BK[b].L :* _cdid_xps_xr(XS, R, crec[c], CE[c].rf))) :/ BK[b].wsum : J(1, 0, .))
            }
            else BK[b].mean_d = J(1, q, 0)
        }
        if (q == 0) {
            GR[g].beta = J(0, 1, .)
            GR[g].M = J(0, 0, .)
            GR[g].keep = J(1, 0, .)
            for (j = 1; j <= rows(GR[g].blocks); j++) BK[GR[g].blocks[j]].alpha = my[j]
            continue
        }
        Xc = J(0, q, .)
        yc = J(0, 1, .)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            sw = sqrt(CE[c].w_syn)
            D = BK[b].L :* _cdid_xps_xr(XS, R, crec[c], CE[c].rf)
            Xc = Xc \ (D :- BK[b].mean_d) :* sw
            yc = yc \ (CE[c].Z :- my[j]) :* sw
        }
        cdid_or_alias_keep(Xc, rank_tol, "a control cohort", keep = ., drop = ., A = .)
        if (cols(drop)) {
            err = sc = 0
            for (j = 1; j <= rows(GR[g].blocks); j++) {
                b = GR[g].blocks[j]
                c = BK[b].cell
                cst = (icpt ? BK[b].mean_d[drop] - BK[b].mean_d[keep] * A : J(1, cols(drop), 0))
                D = BK[b].L :* _cdid_xps_xr(XS, R, crec[c], PR[crec[c]].v)
                cdid_or_alias_gap(D, keep, drop, A, cst, err, sc)
            }
            if (err >= . | err > 64 * epsilon(1) * max((1, q)) * max((1, sc)))
                _cdid_ipw_fail("The counterfactual outcome-regression design is rank deficient on the frame rows of a control cohort without an exact alias that also holds where the regression is predicted; remove the redundant covariate.", 459)
        }
        cdid_or_alias_solve(Xc, yc, keep, q, GR[g].beta, GR[g].M)
        GR[g].keep = keep
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            BK[b].alpha = (icpt ? my[j] - BK[b].mean_d * GR[g].beta : 0)
        }
    }
}

/* Frame values of recipe r for the frame rows of the cohorts read by its links; missing values
   there are refused, other rows are filled with fill. */
real matrix _cdid_xprd_slice(real matrix X, real scalar R, real scalar r, real colvector inscope,
    real scalar fill, string scalar what)
{
    real matrix V
    real colvector k
    V = _cdid_xps_xr(X, R, r, (1::rows(X)))
    if (cols(V) == 0) return(V)
    k = selectindex(inscope)
    if (rows(k)) if (missing(V[k, ])) _cdid_ipw_fail(what + " must be nonmissing on every frame row of the cohorts compared at a reference date.", 459)
    k = selectindex(!inscope)
    if (rows(k)) V[k, ] = J(rows(k), cols(V), fill)
    return(V)
}

/* Link moments, frame-unit influence functions and level loadings of the cross-pair catalogue
   under imputed outcome changes (optionally augmented) when each link reads covariates and
   weights through its recipe (R evaluate_reconstruction_regimes_core on the recipe-cohort
   layout): every (cohort, recipe) pair is a cohort of its own with the recipe's weights,
   designs, repeat rates and masses; imputation pooling and observation models stay within a
   recipe; pooled counterfactual slopes are shared by every recipe; influence rows of a unit
   are summed over its recipes. */
struct cdid_xps_moments scalar cdid_xprd_eval(struct cdid_public_bundle scalar B, struct cdid_xp_cat scalar C,
    struct cdid_xps_sel colvector SS, real scalar augmented, string scalar method,
    string scalar pooling, string scalar imp_pooling, real scalar eps, real scalar maxit, real scalar want_U,
    real scalar rank_tol, real scalar chunk, real scalar max_cached, real scalar workspace, real scalar workcap,
    real scalar R, real colvector lrec, real matrix W, real matrix XT, real matrix XS, real scalar icpt,
    | real matrix KP)
{
    struct cdid_xps_moments scalar out
    struct cdid_population scalar pop
    struct cdid_imputation_model colvector mods
    struct cdid_imputation_pair rowvector pairs
    struct cdid_imputation_budget scalar budget
    struct cdid_sel_moments colvector acc
    struct cdid_sel_log_deriv scalar ld
    struct _cdid_xpr_cell colvector CE
    struct _cdid_xpr_blk colvector BK
    struct _cdid_xpr_grp colvector GR
    struct _cdid_xpr_trt colvector TR
    struct _cdid_xpr_mat colvector Limp_b, Lor_b, Ltrt, WR, ZR, XR, MZ, PR
    transmorphic CK, BKK, TK
    real scalar N, L, nA, T, q, isor, isdr, need_or, need_trt, l, c, ck, cg, b, g, j, k, pid, ci, kk, r
    real scalar ncel, nblk, ng, ntr, tid, dp, Tg, Bg, Cv, Bol, Sm, dm, gid, cid, ncoh, ngrp, gx
    real matrix pk, req, hit, Limp_a, Lor_a, La, Ptab, Q, Hm, nu, dor, lamM, Aco, Mone, MF
    real colvector cohF, A, lc, lb, crec, goff, mh, r_or, eff, wm, wj, cols_a, cols_b, used, inscope
    real colvector mhg, xi, rr, sc, pv, wr, a_coh
    real rowvector dbp, dbm, lam
    real matrix Zs, Xt
    pop = B.prepared.population
    N = pop.N
    nA = rows(pop.A_to_F)
    T = rows(C.cal)
    L = rows(C.g)
    cohF = pop.cohort
    A = pop.A
    q = cols(XS) / R
    isor = (method == "reg")
    isdr = (method == "dr")
    need_or = (isor | isdr)
    need_trt = !isor
    if (!icpt & q == 0) _cdid_ipw_fail("The outcome design is empty; add covariates or a constant.", 459)
    if (augmented & rows(SS) != R) _cdid_ipw_fail("Each design block needs its own observation models.", 498)
    /* Recipe frames: weights and designs on the frame rows of the cohorts each recipe reads. */
    WR = _cdid_xpr_mat(R, 1)
    ZR = _cdid_xpr_mat(R, 1)
    XR = _cdid_xpr_mat(R, 1)
    PR = _cdid_xpr_mat(R, 1)
    used = J(R, 1, 0)
    for (r = 1; r <= R; r++) {
        inscope = J(N, 1, 0)
        for (l = 1; l <= L; l++) {
            if (lrec[l] != r) continue
            used[r] = 1
            inscope = inscope :| (C.g[l] >= . ? cohF :>= . : cohF :== C.g[l])
            inscope = inscope :| (C.coh[C.ci[l]] >= . ? cohF :>= . : cohF :== C.coh[C.ci[l]])
        }
        /* Under a restricted reference eligibility context a recipe reads eligible units only. */
        if (args() >= 22) if (rows(KP)) inscope = inscope :& (KP[, r] :!= 0)
        PR[r].v = selectindex(inscope)
        WR[r].v = _cdid_xprd_slice(W, R, r, inscope, 1, "Sampling weights")
        if (any(inscope :& (WR[r].v :<= 0))) _cdid_ipw_fail("Sampling weights must be positive at the dates read by the cross-pair comparisons.", 459)
        ZR[r].v = _cdid_xprd_slice(XS, R, r, inscope, 0, "Covariates")
        if (need_trt) XR[r].v = _cdid_xprd_slice(XT, R, r, inscope, 0, "Covariates")
    }
    /* Stage 1: actual-change imputation cells of each recipe. */
    mods = cdid_imputation_model(R, 1)
    goff = J(R + 1, 1, 0)
    budget = cdid_imputation_default_budget(max((1, min((N, chunk)))), max_cached, workspace, workcap)
    for (r = 1; r <= R; r++) {
        goff[r + 1] = goff[r]
        if (!used[r]) continue
        pk = J(0, 2, .)
        req = J(0, 2, .)
        for (l = 1; l <= L; l++) {
            if (lrec[l] != r) continue
            pid = 0
            for (j = 1; j <= rows(pk); j++) {
                if (pk[j, 1] == C.cs[l] & pk[j, 2] == C.ct[l]) {
                    pid = j
                    break
                }
            }
            if (!pid) {
                pk = pk \ (C.cs[l], C.ct[l])
                pid = rows(pk)
            }
            req = req \ (editmissing(C.coh[C.ci[l]], 0), pid)
        }
        pairs = cdid_imputation_pair(1, rows(pk))
        for (j = 1; j <= rows(pk); j++) pairs[j] = cdid_imputation_make_pair(pk[j, 1], pk[j, 2], 1)
        if (args() >= 22 & rows(KP)) mods[r] = cdid_imputation_fit_static(pop, WR[r].v, ZR[r].v, pairs, req, imp_pooling,
            icpt, rank_tol, budget, J(1, 0, ""), KP[, r])
        else mods[r] = cdid_imputation_fit_static(pop, WR[r].v, ZR[r].v, pairs, req, imp_pooling, icpt, rank_tol, budget)
        goff[r + 1] = goff[r] + mods[r].n_groups
    }
    ngrp = goff[R + 1]
    ncoh = 0
    for (r = 1; r <= R; r++) if (used[r]) ncoh = max((ncoh, mods[r].n_cohorts))
    Aco = J(ncoh, R, 0)
    MF = J(ncoh, R, 0)
    Mone = J(ncoh, R, 0)
    MZ = _cdid_xpr_mat(R, 1)
    for (r = 1; r <= R; r++) {
        MZ[r].v = J(ncoh, q, 0)
        if (!used[r]) continue
        Aco[1..mods[r].n_cohorts, r] = mods[r].cohort_a_F[1, ]'
        MF[1..mods[r].n_cohorts, r] = mods[r].cohort_mass_F[1, ]'
        for (ci = 1; ci <= mods[r].n_cohorts; ci++) {
            rr = mods[r].coh_u[ci].rows
            if (rows(rr) == 0) continue
            Mone[ci, r] = sum(WR[r].v[rr] :* Aco[ci, r]) / N
            if (q) MZ[r].v[ci, ] = colsum((WR[r].v[rr] :* Aco[ci, r]) :* ZR[r].v[rr, ]) :/ N
        }
    }
    /* Cells (cohort, pair, recipe) in order of first use by the links. */
    CK = asarray_create("real", 4)
    lc = J(L, 1, 0)
    crec = J(0, 1, .)
    hit = J(0, 1, .)
    ncel = 0
    for (l = 1; l <= L; l++) {
        r = lrec[l]
        k = editmissing(C.coh[C.ci[l]], 0)
        if (!asarray_contains(CK, (k, C.cs[l], C.ct[l], r))) {
            cid = 0
            for (j = 1; j <= mods[r].n_cells; j++) {
                if (mods[r].cells[j].cohort == k & mods[r].cells[j].s_index == C.cs[l] & mods[r].cells[j].t_index == C.ct[l]) {
                    cid = j
                    break
                }
            }
            if (!cid) _cdid_ipw_fail("An imputation cell of the cross-pair catalogue was not fitted.", 498)
            ncel++
            asarray(CK, (k, C.cs[l], C.ct[l], r), ncel)
            hit = hit \ cid
            crec = crec \ r
        }
        lc[l] = asarray(CK, (k, C.cs[l], C.ct[l], r))
    }
    for (c = 1; c <= ncel; c++) {
        r = crec[c]
        ci = mods[r].cells[hit[c]].cohort_index
        if (Aco[ci, r] <= 0 | MF[ci, r] <= 0)
            _cdid_ipw_fail("Cohort " + mods[r].cohort_labels[ci] + " has no units observed at two or more dates; imputed changes need repeat observers.", 459)
    }
    if (augmented) {
        acc = cdid_sel_moments(R, 1)
        for (r = 1; r <= R; r++) if (used[r]) acc[r] = cdid_sel_moments_init(SS[r].fit, L)
    }
    CE = _cdid_xpr_cell(ncel, 1)
    for (c = 1; c <= ncel; c++) {
        r = crec[c]
        cid = hit[c]
        gid = mods[r].cells[cid].group_id
        ci = mods[r].cells[cid].cohort_index
        CE[c].cid = cid
        CE[c].ci = ci
        CE[c].gid = gid
        CE[c].L = mods[r].cells[cid].L
        CE[c].cs = mods[r].cells[cid].s_index
        CE[c].ct = mods[r].cells[cid].t_index
        CE[c].a = Aco[ci, r]
        CE[c].alpha = mods[r].cells[cid].alpha
        CE[c].wsum = mods[r].cells[cid].weight_sum
        CE[c].mean_d = mods[r].cells[cid].mean_d
        CE[c].rf = mods[r].coh_u[ci].rows
        CE[c].w_raw = WR[r].v[CE[c].rf]
        CE[c].w_syn = CE[c].w_raw :* CE[c].a
        k = cols(mods[r].groups[gid].kept)
        CE[c].d_imp = (k ? CE[c].L :* ZR[r].v[CE[c].rf, mods[r].groups[gid].kept] : J(rows(CE[c].rf), 0, .))
        CE[c].jhat = CE[c].alpha :+ (k ? CE[c].d_imp * mods[r].groups[gid].beta : J(rows(CE[c].rf), 1, 0))
        Q = _cdid_xpr_pairs(mods[r], mods[r].cells[cid].pair_id, CE[c].rf)
        CE[c].obs = selectindex(Q[, 1] :== 1)
        CE[c].rows_obs = CE[c].rf[CE[c].obs]
        CE[c].a_obs = pop.F_to_A[CE[c].rows_obs]
        CE[c].e_obs = Q[CE[c].obs, 2] - CE[c].jhat[CE[c].obs]
        CE[c].dtil_obs = (k ? CE[c].d_imp[CE[c].obs, ] :- (icpt ? CE[c].mean_d : J(1, k, 0)) : J(rows(CE[c].obs), 0, .))
        CE[c].Z = CE[c].jhat
        CE[c].omspi = J(rows(CE[c].rf), 1, 1)
        CE[c].pi_obs = J(0, 1, .)
        CE[c].T = J(0, 0, .)
        CE[c].th = J(0, 1, .)
        if (augmented & rows(CE[c].obs)) {
            CE[c].kpos = _cdid_xps_kpos(SS[r].rdr, (mods[r].cohorts[ci] == 0 ? . : mods[r].cohorts[ci]))
            CE[c].pc = asarray(SS[r].pcl, (CE[c].kpos, CE[c].cs, CE[c].ct))
            ld = cdid_sel_log_deriv_eval(SS[r].fit, SS[r].rdr, SS[r].um, CE[c].pc, CE[c].rows_obs)
            if (missing(ld.p_weight)) _cdid_ipw_fail("Nonfinite fitted observation probabilities on observed pairs.", 504)
            CE[c].pi_obs = 1 :/ ld.p_weight
            CE[c].T = ld.d_theta_local
            CE[c].th = _cdid_sel_pc_global_idx(SS[r].fit, SS[r].um, CE[c].pc)
            CE[c].Z[CE[c].obs] = CE[c].jhat[CE[c].obs] + CE[c].pi_obs :* CE[c].e_obs
            CE[c].omspi[CE[c].obs] = 1 :- CE[c].pi_obs
        }
    }
    Q = J(0, 0, .)
    /* Counterfactual outcome blocks, one per control cell, in control-link order. */
    lb = J(L, 1, 0)
    nblk = 0
    ng = 0
    if (need_or) {
        BKK = asarray_create("real", 1)
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            if (!asarray_contains(BKK, lc[l])) {
                nblk++
                asarray(BKK, lc[l], nblk)
            }
            lb[l] = asarray(BKK, lc[l])
        }
        BK = _cdid_xpr_blk(nblk, 1)
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            BK[lb[l]].cell = lc[l]
            BK[lb[l]].L = C.ct[l] - C.cs[l]
            BK[lb[l]].trows = J(0, 1, .)
        }
        if (pooling == "pooled" & nblk) {
            GR = _cdid_xpr_grp(1, 1)
            GR[1].blocks = (1::nblk)
            for (b = 1; b <= nblk; b++) BK[b].grp = 1
        }
        else {
            GR = _cdid_xpr_grp(nblk, 1)
            for (b = 1; b <= nblk; b++) {
                GR[b].blocks = b
                BK[b].grp = b
            }
        }
        if (nblk) _cdid_xprd_or_fit(BK, CE, crec, GR, XS, R, icpt, rank_tol, PR)
        ng = rows(GR)
    }
    out.n_or_groups = ng
    out.n_or_blocks = nblk
    out.mhat = J(L, 1, .)
    out.Phi = J(N, L, 0)
    out.U = (want_U ? J(T * nA, L, 0) : J(0, 0, .))
    Limp_a = J(ncel, L, 0)
    Limp_b = _cdid_xpr_mat(max((ngrp, 1)), 1)
    for (r = 1; r <= R; r++) {
        if (!used[r]) continue
        for (g = 1; g <= mods[r].n_groups; g++) Limp_b[goff[r] + g].v = J(cols(mods[r].groups[g].kept), L, 0)
    }
    Lor_a = J(nblk, L, 0)
    Lor_b = _cdid_xpr_mat(max((ng, 1)), 1)
    for (g = 1; g <= ng; g++) Lor_b[g].v = J(q, L, 0)
    La = J(ncoh * R, L, 0)
    TK = asarray_create("real", 3)
    TR = _cdid_xpr_trt(L, 1)
    Ltrt = _cdid_xpr_mat(L, 1)
    ntr = 0
    for (l = 1; l <= L; l++) {
        c = lc[l]
        r = lrec[l]
        wr = WR[r].v
        gx = goff[r] + CE[c].gid
        if (C.role[l] == 1) {
            dp = sum(CE[c].w_syn) / N
            Tg = sum(CE[c].w_syn :* CE[c].Z) / (N * dp)
            out.mhat[l] = Tg
            out.Phi[CE[c].rf, l] = out.Phi[CE[c].rf, l] + (CE[c].w_syn :/ dp) :* (CE[c].Z :- Tg)
            wj = CE[c].w_syn :* CE[c].omspi :/ (N * dp)
            if (icpt) Limp_a[c, l] = Limp_a[c, l] + sum(wj)
            if (cols(CE[c].d_imp)) Limp_b[gx].v[, l] = Limp_b[gx].v[, l] + quadcross(CE[c].d_imp, wj)
            if (augmented & rows(CE[c].obs)) {
                xi = (CE[c].w_syn[CE[c].obs] :/ (N * dp)) :* CE[c].pi_obs
                if (rows(CE[c].th)) acc[r].G_theta[CE[c].th, l] = acc[r].G_theta[CE[c].th, l] + quadcross(CE[c].T, xi :* CE[c].e_obs)
                if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, xi, l)
            }
            continue
        }
        /* Control link: target cohort g, control cell c of cohort k, recipe r. */
        cg = 0
        for (ci = 1; ci <= mods[r].n_cohorts; ci++) if (mods[r].cohorts[ci] == editmissing(C.g[l], 0)) cg = ci
        if (!cg) _cdid_ipw_fail("A target cohort is absent from the imputation fit.", 498)
        rr = mods[r].coh_u[cg].rows
        dp = Mone[cg, r]
        Bg = 0
        Cv = 0
        Bol = 0
        if (need_or) {
            b = lb[l]
            g = BK[b].grp
            dbp = (q ? (BK[b].L :* MZ[r].v[cg, ]) :/ dp : J(1, 0, .))
            Bg = BK[b].alpha + (q ? dbp * GR[g].beta : 0)
            mhg = BK[b].alpha :+ (q ? (BK[b].L :* ZR[r].v[rr, ]) * GR[g].beta : J(rows(rr), 1, 0))
            out.Phi[rr, l] = out.Phi[rr, l] + ((wr[rr] :* Aco[cg, r]) :/ dp) :* (mhg :- Bg)
        }
        if (need_trt) {
            ck = CE[c].ci
            if (asarray_contains(TK, (r, cg, ck))) tid = asarray(TK, (r, cg, ck))
            else {
                ntr++
                tid = ntr
                asarray(TK, (r, cg, ck), tid)
                TR[tid].cig = (r - 1) * ncoh + cg
                TR[tid].cik = (r - 1) * ncoh + ck
                TR[tid].rows = sort(rr \ mods[r].coh_u[ck].rows, 1)
                TR[tid].y = (mods[r].unit_cohort_index[TR[tid].rows] :== cg)
                TR[tid].w_syn = wr[TR[tid].rows] :* Aco[mods[r].unit_cohort_index[TR[tid].rows], r]
                TR[tid].X = XR[r].v[TR[tid].rows, ]
                _cdid_sel_link_trt_logit(TR[tid].X, TR[tid].y, TR[tid].w_syn, eps, maxit, pv = ., Hm = .)
                TR[tid].beta = pv
                Ptab = _cdid_ipw_logit_values(TR[tid].X * pv)
                TR[tid].prob = Ptab[, 1]
                TR[tid].odds = Ptab[, 2]
                TR[tid].act = Ptab[, 3]
                TR[tid].H = quadcross(TR[tid].X, TR[tid].w_syn :* TR[tid].prob :* (1 :- TR[tid].prob), TR[tid].X) :/ N
                TR[tid].pos_g = selectindex(TR[tid].y :== 1)
                TR[tid].pos_k = selectindex(TR[tid].y :== 0)
                TR[tid].C_a = (colsum((wr[TR[tid].rows] :* (TR[tid].y :- TR[tid].prob) :* TR[tid].X)[TR[tid].pos_g, ])' :/ N,
                    colsum((wr[TR[tid].rows] :* (TR[tid].y :- TR[tid].prob) :* TR[tid].X)[TR[tid].pos_k, ])' :/ N)
                Ltrt[tid].v = J(cols(TR[tid].X), L, 0)
            }
            wm = CE[c].w_syn :* TR[tid].odds[TR[tid].pos_k]
            Sm = sum(wm)
            dm = Sm / N
            if (dm <= 0 | dm >= .) _cdid_ipw_fail("The synthetic control normalizer is not positive.", 459)
            Cv = sum(wm :* CE[c].Z) / Sm
            if (isdr) {
                dor = (q ? BK[b].L :* ZR[r].v[CE[c].rf, ] : J(rows(CE[c].rf), 0, .))
                mh = BK[b].alpha :+ (q ? dor * GR[g].beta : J(rows(CE[c].rf), 1, 0))
                Bol = sum(wm :* mh) / Sm
                dbm = (q ? colsum(wm :* dor) :/ Sm : J(1, 0, .))
                eff = (CE[c].Z :- mh) :- (Cv - Bol)
            }
            else eff = CE[c].Z :- Cv
            out.Phi[CE[c].rf, l] = out.Phi[CE[c].rf, l] + (wm :/ dm) :* eff
            La[(r - 1) * ncoh + CE[c].ci, l] = La[(r - 1) * ncoh + CE[c].ci, l] + sum(wm :* eff) / (N * dm * CE[c].a)
            Ltrt[tid].v[, l] = Ltrt[tid].v[, l] + (colsum((wm :* TR[tid].act[TR[tid].pos_k] :* eff) :* TR[tid].X[TR[tid].pos_k, ]) :/ (N * dm))'
            wj = wm :* CE[c].omspi :/ (N * dm)
            if (icpt) Limp_a[c, l] = Limp_a[c, l] + sum(wj)
            if (cols(CE[c].d_imp)) Limp_b[gx].v[, l] = Limp_b[gx].v[, l] + quadcross(CE[c].d_imp, wj)
            if (augmented & rows(CE[c].obs)) {
                xi = (wm[CE[c].obs] :/ (N * dm)) :* CE[c].pi_obs
                if (rows(CE[c].th)) acc[r].G_theta[CE[c].th, l] = acc[r].G_theta[CE[c].th, l] + quadcross(CE[c].T, xi :* CE[c].e_obs)
                if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, xi, l)
            }
        }
        if (isor) {
            out.mhat[l] = Bg
            if (icpt) Lor_a[b, l] = Lor_a[b, l] + 1
            if (q) Lor_b[g].v[, l] = Lor_b[g].v[, l] + dbp'
        }
        else if (isdr) {
            out.mhat[l] = Bg + Cv - Bol
            if (q) Lor_b[g].v[, l] = Lor_b[g].v[, l] + (dbp - dbm)'
        }
        else out.mhat[l] = Cv
    }
    /* Treatment logits: score channel and repeat-rate channel. */
    for (tid = 1; tid <= ntr; tid++) {
        cols_b = selectindex(colsum(Ltrt[tid].v :!= 0)' :> 0)
        if (rows(cols_b) == 0) continue
        nu = cholsolve(TR[tid].H, Ltrt[tid].v[, cols_b'])
        if (missing(nu)) _cdid_ipw_fail("The treatment logit information is singular.", 459)
        sc = TR[tid].w_syn :* (TR[tid].y :- TR[tid].prob)
        out.Phi[TR[tid].rows, cols_b'] = out.Phi[TR[tid].rows, cols_b'] + sc :* (TR[tid].X * nu)
        La[TR[tid].cig, cols_b'] = La[TR[tid].cig, cols_b'] + TR[tid].C_a[, 1]' * nu
        La[TR[tid].cik, cols_b'] = La[TR[tid].cik, cols_b'] + TR[tid].C_a[, 2]' * nu
    }
    /* Counterfactual outcome regressions: intercept and slope adjoints on the control cells. */
    if (nblk) {
        if (icpt & q) {
            for (b = 1; b <= nblk; b++) {
                if (any(Lor_a[b, ] :!= 0)) Lor_b[BK[b].grp].v = Lor_b[BK[b].grp].v - BK[b].mean_d' * Lor_a[b, ]
            }
        }
        for (b = 1; b <= nblk; b++) {
            g = BK[b].grp
            c = BK[b].cell
            r = crec[c]
            gx = goff[r] + CE[c].gid
            lam = (icpt ? Lor_a[b, ] : J(1, L, 0))
            cols_a = selectindex(lam' :!= 0)
            cols_b = (q ? selectindex(colsum(Lor_b[g].v[GR[g].keep, ] :!= 0)' :> 0) : J(0, 1, .))
            if (rows(cols_a) + rows(cols_b) == 0) continue
            dor = (q ? BK[b].L :* ZR[r].v[CE[c].rf, ] : J(rows(CE[c].rf), 0, .))
            r_or = CE[c].Z :- (BK[b].alpha :+ (q ? dor * GR[g].beta : J(rows(CE[c].rf), 1, 0)))
            for (kk = 1; kk <= 2; kk++) {
                if (kk == 1) {
                    if (rows(cols_a) == 0) continue
                    Q = (CE[c].w_syn :/ BK[b].wsum) * lam[cols_a']
                    pv = cols_a
                }
                else {
                    if (rows(cols_b) == 0) continue
                    nu = N :* cholsolve(GR[g].M, Lor_b[g].v[GR[g].keep, cols_b'])
                    if (missing(nu)) _cdid_ipw_fail("The counterfactual outcome-regression Gram is singular.", 459)
                    Q = ((CE[c].w_syn :/ N) :* (dor[, GR[g].keep] :- (icpt ? BK[b].mean_d[GR[g].keep] : J(1, cols(GR[g].keep), 0)))) * nu
                    pv = cols_b
                }
                out.Phi[CE[c].rf, pv'] = out.Phi[CE[c].rf, pv'] + (N :* r_or) :* Q
                La[(r - 1) * ncoh + CE[c].ci, pv'] = La[(r - 1) * ncoh + CE[c].ci, pv'] + colsum((r_or :/ CE[c].a) :* Q)
                lamM = CE[c].omspi :* Q
                if (icpt) Limp_a[c, pv'] = Limp_a[c, pv'] + colsum(lamM)
                if (cols(CE[c].d_imp)) Limp_b[gx].v[, pv'] = Limp_b[gx].v[, pv'] + quadcross(CE[c].d_imp, lamM)
                if (augmented & rows(CE[c].obs)) {
                    if (rows(CE[c].th)) acc[r].G_theta[CE[c].th, pv'] = acc[r].G_theta[CE[c].th, pv'] +
                        quadcross(CE[c].T, (CE[c].pi_obs :* CE[c].e_obs) :* Q[CE[c].obs, ])
                    if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, CE[c].pi_obs :* Q[CE[c].obs, ], pv)
                }
            }
        }
    }
    /* Stage 1 adjoints on the observed pairs of each cell. */
    if (icpt) {
        for (c = 1; c <= ncel; c++) {
            gx = goff[crec[c]] + CE[c].gid
            if (cols(CE[c].d_imp) & any(Limp_a[c, ] :!= 0))
                Limp_b[gx].v = Limp_b[gx].v - CE[c].mean_d' * Limp_a[c, ]
        }
    }
    for (c = 1; c <= ncel; c++) {
        if (rows(CE[c].obs) == 0) continue
        r = crec[c]
        gid = CE[c].gid
        gx = goff[r] + gid
        lam = (icpt ? Limp_a[c, ] : J(1, L, 0))
        cols_a = selectindex(lam' :!= 0)
        cols_b = (cols(CE[c].d_imp) ? selectindex(colsum(Limp_b[gx].v :!= 0)' :> 0) : J(0, 1, .))
        for (kk = 1; kk <= 2; kk++) {
            if (kk == 1) {
                if (rows(cols_a) == 0) continue
                Q = (CE[c].w_raw[CE[c].obs] :/ CE[c].wsum) * lam[cols_a']
                pv = cols_a
            }
            else {
                if (rows(cols_b) == 0) continue
                nu = _cdid_nuisance_gram(mods[r].groups[gid].R, mods[r].groups[gid].pivot, N :* Limp_b[gx].v[, cols_b'])
                Q = ((CE[c].w_raw[CE[c].obs] :/ N) :* CE[c].dtil_obs) * nu
                pv = cols_b
            }
            out.Phi[CE[c].rows_obs, pv'] = out.Phi[CE[c].rows_obs, pv'] + (N :* CE[c].e_obs) :* Q
            if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, Q, pv)
        }
    }
    /* Repeat-rate channels of every (cohort, recipe) and observation-model channels. */
    for (r = 1; r <= R; r++) {
        if (!used[r]) continue
        for (ci = 1; ci <= mods[r].n_cohorts; ci++) {
            cols_a = selectindex(La[(r - 1) * ncoh + ci, ]' :!= 0)
            if (rows(cols_a) == 0) continue
            rr = mods[r].coh_u[ci].rows
            xi = (WR[r].v[rr] :/ MF[ci, r]) :* ((A[rr] :== 1) :- Aco[ci, r])
            out.Phi[rr, cols_a'] = out.Phi[rr, cols_a'] + xi * La[(r - 1) * ncoh + ci, cols_a']
        }
        if (augmented) {
            if (max(abs(acc[r].G_theta)) > 0) out.Phi = out.Phi + cdid_sel_eval_loading(SS[r].fit, SS[r].rdr, acc[r], (1::N))
        }
    }
    out.n_trt_fits = ntr
    return(out)
}

end
