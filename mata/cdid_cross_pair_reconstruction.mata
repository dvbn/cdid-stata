*! Cross-pair chained DiD with imputed outcome changes (selmethod(imputation)) and with the
*! augmented observation correction (selmethod(dr)).
version 19.5
mata:

/* One imputation cell (cohort, date pair) evaluated on every frame row of its cohort. */
struct _cdid_xpr_cell {
    real scalar cid, ci, gid, L, cs, ct, a, alpha, wsum, pc, kpos
    real colvector rf, w_raw, w_syn, jhat, Z, obs, rows_obs, a_obs, e_obs, pi_obs, omspi, th
    real matrix d_imp, dtil_obs, T
    real rowvector mean_d
}

/* Counterfactual outcome block (one control cell) fitted on its cohort's frame rows. */
struct _cdid_xpr_blk {
    real scalar cell, L, grp, alpha, wsum
    real rowvector mean_d
    real colvector trows
}

struct _cdid_xpr_grp {
    real colvector beta, blocks
    real matrix M
    real rowvector keep
}

/* Treatment logit of the target cohort against one control cohort on their frame rows. */
struct _cdid_xpr_trt {
    real scalar cig, cik
    real colvector rows, y, w_syn, prob, odds, act, pos_g, pos_k, beta
    real matrix X, H, C_a
}

struct _cdid_xpr_mat {
    real matrix v
}

/* Observed pairs (indicator, change) of frame rows for one model pair, read in chunks. */
real matrix _cdid_xpr_pairs(struct cdid_imputation_model scalar m, real scalar pid, real colvector rows)
{
    real matrix out
    real scalar lo, hi, B
    out = J(rows(rows), 2, 0)
    B = m.budget.max_chunk_rows
    for (lo = 1; lo <= rows(rows); lo = lo + B) {
        hi = min((rows(rows), lo + B - 1))
        out[lo..hi, ] = cdid_imputation_read_pair_chunk(m.reader, m.pop, m.pairs[pid], rows[lo..hi], 1)[, 1..2]
    }
    return(out)
}

/* Level loadings of a pair difference on the analysis units observed at both dates. */
void _cdid_xpr_addU(real matrix U, real scalar nA, real colvector a, real scalar cs, real scalar ct,
    real matrix x, real colvector cols)
{
    if (rows(a) == 0 | cols(x) == 0) return
    U[(ct - 1) * nA :+ a, cols'] = U[(ct - 1) * nA :+ a, cols'] + x
    U[(cs - 1) * nA :+ a, cols'] = U[(cs - 1) * nA :+ a, cols'] - x
}

/* Counterfactual outcome regressions of the imputed responses Z on the L-scaled slopes over
   the control cohorts' frame rows (weights w a_k), with block intercepts; separate or pooled
   slopes. Refuses a rank-deficient design. */
void _cdid_xpr_or_fit(struct _cdid_xpr_blk colvector BK, struct _cdid_xpr_cell colvector CE,
    struct _cdid_xpr_grp colvector GR, real matrix Xs, real scalar icpt, real scalar rank_tol)
{
    struct _cdid_xps_blk colvector BX
    real scalar g, j, b, q, rk, mx, c, k
    real matrix Xc, H, Rf, D, Dfit
    real colvector yc, qy, sw, my, bfit
    real rowvector tau, pv, keep, drop
    q = cols(Xs)
    for (g = 1; g <= rows(GR); g++) {
        my = J(rows(GR[g].blocks), 1, 0)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            BK[b].wsum = sum(CE[c].w_syn)
            if (icpt) {
                my[j] = sum(CE[c].w_syn :* CE[c].Z) / BK[b].wsum
                BK[b].mean_d = (q ? colsum(CE[c].w_syn :* (BK[b].L :* Xs[CE[c].rf, ])) :/ BK[b].wsum : J(1, 0, .))
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
        Dfit = J(0, q, .)
        bfit = J(0, 1, .)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            sw = sqrt(CE[c].w_syn)
            D = BK[b].L :* Xs[CE[c].rf, ]
            Dfit = Dfit \ D
            bfit = bfit \ J(rows(sw), 1, j)
            Xc = Xc \ (D :- BK[b].mean_d) :* sw
            yc = yc \ (CE[c].Z :- my[j]) :* sw
        }
        hqrdp(Xc, H = ., tau = ., Rf = ., pv = .)
        mx = max(abs(diagonal(Rf)))
        rk = sum(abs(diagonal(Rf)) :> rank_tol * mx)
        if (mx <= 0 | rk < 1) _cdid_ipw_fail("The counterfactual outcome-regression design has no numerical rank on the frame rows of a control cohort.", 459)
        keep = sort(pv[1..rk]', 1)'
        if (rk < q) {
            drop = sort(pv[(rk + 1)..q]', 1)'
            BX = _cdid_xps_blk(rows(BK), 1)
            for (j = 1; j <= rows(GR[g].blocks); j++) {
                b = GR[g].blocks[j]
                BX[b].L = BK[b].L
                BX[b].md = BK[b].mean_d
                BX[b].trows = BK[b].trows
            }
            cdid_xp_alias_check(Xc, keep, drop, icpt, BX, GR[g].blocks, Dfit, bfit, Xs)
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
            BK[b].alpha = (icpt ? my[j] - BK[b].mean_d * GR[g].beta : 0)
        }
    }
}

/* Link moments, frame-unit influence functions and level loadings of the cross-pair
   catalogue under imputed outcome changes, optionally augmented by inverse observation
   probabilities on the observed pairs (R evaluate_reconstruction_regimes_core on the
   cross-pair descriptors). Target links: mean of the cohort's response Z over its frame
   rows (weights w a_g); control links: control response reweighted toward the target
   cohort by a treatment logit (IPW), a counterfactual outcome regression (OR), or both. */
struct cdid_xps_moments scalar cdid_xpr_eval(struct cdid_public_bundle scalar B, struct cdid_xp_cat scalar C,
    struct cdid_xps_sel scalar S, real scalar augmented, real matrix Y, string scalar method,
    string scalar pooling, string scalar imp_pooling, real scalar eps, real scalar maxit, real scalar want_U,
    real scalar rank_tol, real scalar chunk, real scalar max_cached, real scalar workspace, real scalar workcap)
{
    struct cdid_xps_moments scalar out
    struct cdid_population scalar pop
    struct cdid_imputation_model scalar m
    struct cdid_imputation_pair rowvector pairs
    struct cdid_imputation_budget scalar budget
    struct cdid_sel_moments scalar acc
    struct cdid_sel_log_deriv scalar ld
    struct _cdid_xpr_cell colvector CE
    struct _cdid_xpr_blk colvector BK
    struct _cdid_xpr_grp colvector GR
    struct _cdid_xpr_trt colvector TR
    struct _cdid_xpr_mat colvector Limp_b, Lor_b, Ltrt
    transmorphic TK, BKK
    real scalar N, L, nA, T, p, q, icpt, isor, isdr, need_or, need_trt, l, c, ck, cg, b, g, j, k, pid, ci, kk
    real scalar ncel, nblk, ng, ntr, tid, dp, Tg, Bg, Cv, Bol, Sm, dm, nreq, gid, cid, lo
    real matrix pk, X, Zs, req, hit, Limp_a, Lor_a, La, Ptab, Q, Hm, nu, dor, lamM
    real colvector w, cohF, A, cellmap, lc, lb, lt, coh0, a_coh, mF_coh, M_one, mh, r_or, eff, wm, wj, cols_a, cols_b
    real colvector mhg, xi, rr, sc, pv
    real rowvector dbp, dbm, lam
    real matrix M_z
    pop = B.prepared.population
    N = pop.N
    nA = rows(pop.A_to_F)
    T = rows(C.cal)
    L = rows(C.g)
    w = B.prepared.weight.owner_weight
    cohF = pop.cohort
    A = pop.A
    icpt = B.prepared.x_has_intercept
    p = cols(B.prepared.X)
    X = B.prepared.X
    Zs = (icpt ? (p > 1 ? X[, 2..p] : J(N, 0, .)) : X)
    q = cols(Zs)
    isor = (method == "reg")
    isdr = (method == "dr")
    need_or = (isor | isdr)
    need_trt = !isor
    if (!icpt & q == 0) _cdid_ipw_fail("The outcome design is empty; add covariates or a constant.", 459)
    /* Stage 1: actual-change imputation cells requested by the links, in link order. */
    pk = J(0, 2, .)
    req = J(L, 2, .)
    for (l = 1; l <= L; l++) {
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
        req[l, ] = (editmissing(C.coh[C.ci[l]], 0), pid)
    }
    pairs = cdid_imputation_pair(1, rows(pk))
    for (j = 1; j <= rows(pk); j++) pairs[j] = cdid_imputation_make_pair(pk[j, 1], pk[j, 2], 1)
    budget = cdid_imputation_default_budget(max((1, min((N, chunk)))), max_cached, workspace, workcap)
    m = cdid_imputation_fit_static(pop, w, Zs, pairs, req, imp_pooling, icpt, rank_tol, budget)
    /* Cells in order of first use by the links. */
    cellmap = J(m.n_cells, 1, 0)
    lc = J(L, 1, 0)
    ncel = 0
    for (l = 1; l <= L; l++) {
        cid = 0
        for (j = 1; j <= m.n_cells; j++) {
            if (m.cells[j].cohort == req[l, 1] & m.cells[j].s_index == C.cs[l] & m.cells[j].t_index == C.ct[l]) {
                cid = j
                break
            }
        }
        if (!cid) _cdid_ipw_fail("An imputation cell of the cross-pair catalogue was not fitted.", 498)
        if (!cellmap[cid]) {
            ncel++
            cellmap[cid] = ncel
        }
        lc[l] = cellmap[cid]
    }
    a_coh = m.cohort_a_F[1, ]'
    mF_coh = m.cohort_mass_F[1, ]'
    for (ci = 1; ci <= m.n_cohorts; ci++) {
        if (any(m.unit_cohort_index :== ci) & (a_coh[ci] <= 0 | mF_coh[ci] <= 0)) {
            if (any(req[, 1] :== m.cohorts[ci]))
                _cdid_ipw_fail("Cohort " + m.cohort_labels[ci] + " has no units observed at two or more dates; imputed changes need repeat observers.", 459)
        }
    }
    M_one = J(m.n_cohorts, 1, 0)
    M_z = J(m.n_cohorts, q, 0)
    for (ci = 1; ci <= m.n_cohorts; ci++) {
        rr = m.coh_u[ci].rows
        if (rows(rr) == 0) continue
        M_one[ci] = sum(w[rr] :* a_coh[ci]) / N
        if (q) M_z[ci, ] = colsum((w[rr] :* a_coh[ci]) :* Zs[rr, ]) :/ N
    }
    if (augmented) acc = cdid_sel_moments_init(S.fit, L)
    CE = _cdid_xpr_cell(ncel, 1)
    for (cid = 1; cid <= m.n_cells; cid++) {
        c = cellmap[cid]
        if (!c) continue
        gid = m.cells[cid].group_id
        ci = m.cells[cid].cohort_index
        CE[c].cid = cid
        CE[c].ci = ci
        CE[c].gid = gid
        CE[c].L = m.cells[cid].L
        CE[c].cs = m.cells[cid].s_index
        CE[c].ct = m.cells[cid].t_index
        CE[c].a = a_coh[ci]
        CE[c].alpha = m.cells[cid].alpha
        CE[c].wsum = m.cells[cid].weight_sum
        CE[c].mean_d = m.cells[cid].mean_d
        CE[c].rf = m.coh_u[ci].rows
        CE[c].w_raw = w[CE[c].rf]
        CE[c].w_syn = CE[c].w_raw :* CE[c].a
        k = cols(m.groups[gid].kept)
        CE[c].d_imp = (k ? CE[c].L :* Zs[CE[c].rf, m.groups[gid].kept] : J(rows(CE[c].rf), 0, .))
        CE[c].jhat = CE[c].alpha :+ (k ? CE[c].d_imp * m.groups[gid].beta : J(rows(CE[c].rf), 1, 0))
        hit = _cdid_xpr_pairs(m, m.cells[cid].pair_id, CE[c].rf)
        CE[c].obs = selectindex(hit[, 1] :== 1)
        CE[c].rows_obs = CE[c].rf[CE[c].obs]
        CE[c].a_obs = pop.F_to_A[CE[c].rows_obs]
        CE[c].e_obs = hit[CE[c].obs, 2] - CE[c].jhat[CE[c].obs]
        CE[c].dtil_obs = (k ? CE[c].d_imp[CE[c].obs, ] :- (icpt ? CE[c].mean_d : J(1, k, 0)) : J(rows(CE[c].obs), 0, .))
        CE[c].Z = CE[c].jhat
        CE[c].omspi = J(rows(CE[c].rf), 1, 1)
        CE[c].pi_obs = J(0, 1, .)
        CE[c].T = J(0, 0, .)
        CE[c].th = J(0, 1, .)
        if (augmented & rows(CE[c].obs)) {
            CE[c].kpos = _cdid_xps_kpos(S.rdr, (m.cohorts[ci] == 0 ? . : m.cohorts[ci]))
            CE[c].pc = asarray(S.pcl, (CE[c].kpos, CE[c].cs, CE[c].ct))
            ld = cdid_sel_log_deriv_eval(S.fit, S.rdr, S.um, CE[c].pc, CE[c].rows_obs)
            if (missing(ld.p_weight)) _cdid_ipw_fail("Nonfinite fitted observation probabilities on observed pairs.", 504)
            CE[c].pi_obs = 1 :/ ld.p_weight
            CE[c].T = ld.d_theta_local
            CE[c].th = _cdid_sel_pc_global_idx(S.fit, S.um, CE[c].pc)
            CE[c].Z[CE[c].obs] = CE[c].jhat[CE[c].obs] + CE[c].pi_obs :* CE[c].e_obs
            CE[c].omspi[CE[c].obs] = 1 :- CE[c].pi_obs
        }
    }
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
        for (b = 1; b <= nblk; b++) BK[b].trows = J(0, 1, .)
        for (l = 1; l <= L; l++) {
            if (C.role[l] != 2) continue
            BK[lb[l]].cell = lc[l]
            BK[lb[l]].L = C.ct[l] - C.cs[l]
            for (ci = 1; ci <= m.n_cohorts; ci++) {
                if (m.cohorts[ci] == C.g[l]) BK[lb[l]].trows = uniqrows(BK[lb[l]].trows \ m.coh_u[ci].rows)
            }
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
        if (nblk) _cdid_xpr_or_fit(BK, CE, GR, Zs, icpt, rank_tol)
        ng = rows(GR)
    }
    out.n_or_groups = ng
    out.n_or_blocks = nblk
    out.mhat = J(L, 1, .)
    out.Phi = J(N, L, 0)
    out.U = (want_U ? J(T * nA, L, 0) : J(0, 0, .))
    Limp_a = J(ncel, L, 0)
    Limp_b = _cdid_xpr_mat(m.n_groups, 1)
    for (g = 1; g <= m.n_groups; g++) Limp_b[g].v = J(cols(m.groups[g].kept), L, 0)
    Lor_a = J(nblk, L, 0)
    Lor_b = _cdid_xpr_mat(ng, 1)
    for (g = 1; g <= ng; g++) Lor_b[g].v = J(q, L, 0)
    La = J(m.n_cohorts, L, 0)
    TK = asarray_create("real", 2)
    TR = _cdid_xpr_trt(L, 1)
    Ltrt = _cdid_xpr_mat(L, 1)
    ntr = 0
    for (l = 1; l <= L; l++) {
        c = lc[l]
        if (C.role[l] == 1) {
            dp = sum(CE[c].w_syn) / N
            Tg = sum(CE[c].w_syn :* CE[c].Z) / (N * dp)
            out.mhat[l] = Tg
            out.Phi[CE[c].rf, l] = out.Phi[CE[c].rf, l] + (CE[c].w_syn :/ dp) :* (CE[c].Z :- Tg)
            wj = CE[c].w_syn :* CE[c].omspi :/ (N * dp)
            if (icpt) Limp_a[c, l] = Limp_a[c, l] + sum(wj)
            if (cols(CE[c].d_imp)) Limp_b[CE[c].gid].v[, l] = Limp_b[CE[c].gid].v[, l] + quadcross(CE[c].d_imp, wj)
            if (augmented & rows(CE[c].obs)) {
                xi = (CE[c].w_syn[CE[c].obs] :/ (N * dp)) :* CE[c].pi_obs
                if (rows(CE[c].th)) acc.G_theta[CE[c].th, l] = acc.G_theta[CE[c].th, l] + quadcross(CE[c].T, xi :* CE[c].e_obs)
                if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, xi, l)
            }
            continue
        }
        /* Control link: target cohort g, control cell c of cohort k. */
        cg = 0
        for (ci = 1; ci <= m.n_cohorts; ci++) if (m.cohorts[ci] == C.g[l]) cg = ci
        if (!cg) _cdid_ipw_fail("A target cohort is absent from the imputation fit.", 498)
        rr = m.coh_u[cg].rows
        dp = M_one[cg]
        Bg = 0
        Cv = 0
        Bol = 0
        if (need_or) {
            b = lb[l]
            g = BK[b].grp
            dbp = (q ? (BK[b].L :* M_z[cg, ]) :/ dp : J(1, 0, .))
            Bg = BK[b].alpha + (q ? dbp * GR[g].beta : 0)
            mhg = BK[b].alpha :+ (q ? (BK[b].L :* Zs[rr, ]) * GR[g].beta : J(rows(rr), 1, 0))
            out.Phi[rr, l] = out.Phi[rr, l] + ((w[rr] :* a_coh[cg]) :/ dp) :* (mhg :- Bg)
        }
        if (need_trt) {
            ck = CE[c].ci
            if (asarray_contains(TK, (cg, ck))) tid = asarray(TK, (cg, ck))
            else {
                ntr++
                tid = ntr
                asarray(TK, (cg, ck), tid)
                TR[tid].cig = cg
                TR[tid].cik = ck
                TR[tid].rows = sort(rr \ m.coh_u[ck].rows, 1)
                TR[tid].y = (m.unit_cohort_index[TR[tid].rows] :== cg)
                TR[tid].w_syn = w[TR[tid].rows] :* a_coh[m.unit_cohort_index[TR[tid].rows]]
                TR[tid].X = X[TR[tid].rows, ]
                _cdid_sel_link_trt_logit(TR[tid].X, TR[tid].y, TR[tid].w_syn, eps, maxit, pv = ., Hm = .)
                TR[tid].beta = pv
                Ptab = _cdid_ipw_logit_values(TR[tid].X * pv)
                TR[tid].prob = Ptab[, 1]
                TR[tid].odds = Ptab[, 2]
                TR[tid].act = Ptab[, 3]
                TR[tid].H = quadcross(TR[tid].X, TR[tid].w_syn :* TR[tid].prob :* (1 :- TR[tid].prob), TR[tid].X) :/ N
                TR[tid].pos_g = selectindex(TR[tid].y :== 1)
                TR[tid].pos_k = selectindex(TR[tid].y :== 0)
                TR[tid].C_a = (colsum((w[TR[tid].rows] :* (TR[tid].y :- TR[tid].prob) :* TR[tid].X)[TR[tid].pos_g, ])' :/ N,
                    colsum((w[TR[tid].rows] :* (TR[tid].y :- TR[tid].prob) :* TR[tid].X)[TR[tid].pos_k, ])' :/ N)
                Ltrt[tid].v = J(cols(X), L, 0)
            }
            /* Control rows of cell c are the cohort rows in frame order, as pos_k. */
            wm = CE[c].w_syn :* TR[tid].odds[TR[tid].pos_k]
            Sm = sum(wm)
            dm = Sm / N
            if (dm <= 0 | dm >= .) _cdid_ipw_fail("The synthetic control normalizer is not positive.", 459)
            Cv = sum(wm :* CE[c].Z) / Sm
            if (isdr) {
                dor = (q ? BK[b].L :* Zs[CE[c].rf, ] : J(rows(CE[c].rf), 0, .))
                mh = BK[b].alpha :+ (q ? dor * GR[g].beta : J(rows(CE[c].rf), 1, 0))
                Bol = sum(wm :* mh) / Sm
                dbm = (q ? colsum(wm :* dor) :/ Sm : J(1, 0, .))
                eff = (CE[c].Z :- mh) :- (Cv - Bol)
            }
            else eff = CE[c].Z :- Cv
            out.Phi[CE[c].rf, l] = out.Phi[CE[c].rf, l] + (wm :/ dm) :* eff
            La[CE[c].ci, l] = La[CE[c].ci, l] + sum(wm :* eff) / (N * dm * CE[c].a)
            Ltrt[tid].v[, l] = Ltrt[tid].v[, l] + (colsum((wm :* TR[tid].act[TR[tid].pos_k] :* eff) :* TR[tid].X[TR[tid].pos_k, ]) :/ (N * dm))'
            wj = wm :* CE[c].omspi :/ (N * dm)
            if (icpt) Limp_a[c, l] = Limp_a[c, l] + sum(wj)
            if (cols(CE[c].d_imp)) Limp_b[CE[c].gid].v[, l] = Limp_b[CE[c].gid].v[, l] + quadcross(CE[c].d_imp, wj)
            if (augmented & rows(CE[c].obs)) {
                xi = (wm[CE[c].obs] :/ (N * dm)) :* CE[c].pi_obs
                if (rows(CE[c].th)) acc.G_theta[CE[c].th, l] = acc.G_theta[CE[c].th, l] + quadcross(CE[c].T, xi :* CE[c].e_obs)
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
            lam = (icpt ? Lor_a[b, ] : J(1, L, 0))
            cols_a = selectindex(lam' :!= 0)
            cols_b = (q ? selectindex(colsum(Lor_b[g].v[GR[g].keep, ] :!= 0)' :> 0) : J(0, 1, .))
            if (rows(cols_a) + rows(cols_b) == 0) continue
            dor = (q ? BK[b].L :* Zs[CE[c].rf, ] : J(rows(CE[c].rf), 0, .))
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
                La[CE[c].ci, pv'] = La[CE[c].ci, pv'] + colsum((r_or :/ CE[c].a) :* Q)
                lamM = CE[c].omspi :* Q
                if (icpt) Limp_a[c, pv'] = Limp_a[c, pv'] + colsum(lamM)
                if (cols(CE[c].d_imp)) Limp_b[CE[c].gid].v[, pv'] = Limp_b[CE[c].gid].v[, pv'] + quadcross(CE[c].d_imp, lamM)
                if (augmented & rows(CE[c].obs)) {
                    if (rows(CE[c].th)) acc.G_theta[CE[c].th, pv'] = acc.G_theta[CE[c].th, pv'] +
                        quadcross(CE[c].T, (CE[c].pi_obs :* CE[c].e_obs) :* Q[CE[c].obs, ])
                    if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, CE[c].pi_obs :* Q[CE[c].obs, ], pv)
                }
            }
        }
    }
    /* Stage 1 adjoints on the observed pairs of each cell. */
    if (icpt) {
        for (c = 1; c <= ncel; c++) {
            if (cols(CE[c].d_imp) & any(Limp_a[c, ] :!= 0))
                Limp_b[CE[c].gid].v = Limp_b[CE[c].gid].v - CE[c].mean_d' * Limp_a[c, ]
        }
    }
    for (c = 1; c <= ncel; c++) {
        if (rows(CE[c].obs) == 0) continue
        gid = CE[c].gid
        lam = (icpt ? Limp_a[c, ] : J(1, L, 0))
        cols_a = selectindex(lam' :!= 0)
        cols_b = (cols(CE[c].d_imp) ? selectindex(colsum(Limp_b[gid].v :!= 0)' :> 0) : J(0, 1, .))
        for (kk = 1; kk <= 2; kk++) {
            if (kk == 1) {
                if (rows(cols_a) == 0) continue
                Q = (CE[c].w_raw[CE[c].obs] :/ CE[c].wsum) * lam[cols_a']
                pv = cols_a
            }
            else {
                if (rows(cols_b) == 0) continue
                nu = _cdid_nuisance_gram(m.groups[gid].R, m.groups[gid].pivot, N :* Limp_b[gid].v[, cols_b'])
                Q = ((CE[c].w_raw[CE[c].obs] :/ N) :* CE[c].dtil_obs) * nu
                pv = cols_b
            }
            out.Phi[CE[c].rows_obs, pv'] = out.Phi[CE[c].rows_obs, pv'] + (N :* CE[c].e_obs) :* Q
            if (want_U) _cdid_xpr_addU(out.U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, Q, pv)
        }
    }
    /* Repeat-rate and observation-model channels. */
    for (ci = 1; ci <= m.n_cohorts; ci++) {
        cols_a = selectindex(La[ci, ]' :!= 0)
        if (rows(cols_a) == 0) continue
        rr = m.coh_u[ci].rows
        xi = (w[rr] :/ mF_coh[ci]) :* ((A[rr] :== 1) :- a_coh[ci])
        out.Phi[rr, cols_a'] = out.Phi[rr, cols_a'] + xi * La[ci, cols_a']
    }
    if (augmented) {
        if (max(abs(acc.G_theta)) > 0) out.Phi = out.Phi + cdid_sel_eval_loading(S.fit, S.rdr, acc, (1::N))
    }
    out.n_trt_fits = ntr
    if (cdid_nuis_requested()) cdid_xpn_recon(B, C, S, augmented, m, CE, TR, ntr, BK, GR, nblk, lc, lb, icpt, Zs, X,
        a_coh, mF_coh, method, pooling, imp_pooling)
    return(out)
}

end
