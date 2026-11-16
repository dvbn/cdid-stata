*! Within-unit GMM weighting for same-pair comparisons with imputed outcome changes
*! (selmethod(imputation)) and with the augmented observation correction (selmethod(dr)).
version 19.5
mata:

/* Counterfactual outcome block: the control cells of one comparison pool on one date pair. */
struct _cdid_wsr_blk {
    real scalar L, grp, alpha, wsum
    real rowvector mean_d, cells
}

/* Treatment logit of a target cohort against the cohorts of its control pool. */
struct _cdid_wsr_trt {
    real colvector rows, uci, odds
}

/* Counterfactual outcome regressions of the imputed responses on the L-scaled slopes over the
   frame rows of the control cells of each block (weights w a_k, slopes of each cell's recipe),
   with block intercepts; separate or pooled slopes. Exactly aliased columns are dropped when the
   alias, with one constant per block when the design has a constant, holds on every frame row of
   the cohorts each recipe reads (PR), where the regression is fitted and predicted; any other
   rank deficiency is refused. */
void _cdid_wsr_or_fit(struct _cdid_wsr_blk colvector BK, struct _cdid_xpr_cell colvector CE,
    real colvector crec, struct _cdid_xpr_grp colvector GR, struct _cdid_xpr_mat colvector ZR,
    real scalar q, real scalar icpt, real scalar rank_tol, struct _cdid_xpr_mat colvector PR)
{
    real scalar g, j, b, c, h, err, sc
    real matrix Xc, A
    real colvector yc, sw, my
    real rowvector keep, drop, cst
    for (g = 1; g <= rows(GR); g++) {
        my = J(rows(GR[g].blocks), 1, 0)
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            BK[b].wsum = 0
            BK[b].mean_d = J(1, q, 0)
            for (h = 1; h <= cols(BK[b].cells); h++) {
                c = BK[b].cells[h]
                BK[b].wsum = BK[b].wsum + sum(CE[c].w_syn)
                if (icpt) {
                    my[j] = my[j] + sum(CE[c].w_syn :* CE[c].Z)
                    if (q) BK[b].mean_d = BK[b].mean_d + colsum(CE[c].w_syn :* (BK[b].L :* ZR[crec[c]].v[CE[c].rf, ]))
                }
            }
            if (BK[b].wsum <= 0 | BK[b].wsum >= .) _cdid_ipw_fail("A control pool has no frame rows with positive weight.", 459)
            my[j] = my[j] / BK[b].wsum
            BK[b].mean_d = BK[b].mean_d :/ BK[b].wsum
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
            for (h = 1; h <= cols(BK[b].cells); h++) {
                c = BK[b].cells[h]
                sw = sqrt(CE[c].w_syn)
                Xc = Xc \ ((BK[b].L :* ZR[crec[c]].v[CE[c].rf, ]) :- BK[b].mean_d) :* sw
                yc = yc \ (CE[c].Z :- my[j]) :* sw
            }
        }
        cdid_or_alias_keep(Xc, rank_tol, "a control pool", keep = ., drop = ., A = .)
        if (cols(drop)) {
            err = sc = 0
            for (j = 1; j <= rows(GR[g].blocks); j++) {
                b = GR[g].blocks[j]
                cst = (icpt ? BK[b].mean_d[drop] - BK[b].mean_d[keep] * A : J(1, cols(drop), 0))
                for (h = 1; h <= cols(BK[b].cells); h++) {
                    c = BK[b].cells[h]
                    cdid_or_alias_gap(BK[b].L :* ZR[crec[c]].v[PR[crec[c]].v, ], keep, drop, A, cst, err, sc)
                }
            }
            if (err >= . | err > 64 * epsilon(1) * max((1, q)) * max((1, sc)))
                _cdid_ipw_fail("The counterfactual outcome-regression design is rank deficient on the frame rows of a control pool without an exact alias that also holds where the regression is predicted; remove the redundant covariate.", 459)
        }
        cdid_or_alias_solve(Xc, yc, keep, q, GR[g].beta, GR[g].M)
        GR[g].keep = keep
        for (j = 1; j <= rows(GR[g].blocks); j++) {
            b = GR[g].blocks[j]
            BK[b].alpha = (icpt ? my[j] - BK[b].mean_d * GR[g].beta : 0)
        }
    }
}

/* Observation-model reader of one recipe: selection design Z (constant first when hi) and
   weights w over the frame rows. */
struct cdid_xps_sel scalar _cdid_wsr_sel_reader(struct cdid_public_bundle scalar B, real matrix Z, real scalar hi,
    real colvector w, real scalar chunk, | real colvector keep)
{
    struct cdid_xps_sel scalar S
    struct cdid_sel_decl scalar decl
    real scalar N, univ_has_rho
    string scalar univ
    N = B.prepared.population.N
    univ = (strtrim(B.prepared.selection_universe) == "" ? "F" : strtrim(B.prepared.selection_universe))
    univ_has_rho = (B.prepared.sel_has_rho != 0 & B.prepared.sel_has_rho < . ? 1 :
        (univ != "F" | (rows(B.prepared.sel_rho) == N & any(B.prepared.sel_rho :!= 1))))
    decl = cdid_sel_norm_decl("conditional_mtar", 1, univ, univ_has_rho,
        (strtrim(B.prepared.selection_justification) != ""), B.prepared.selection_model,
        B.prepared.selection_pooling, B.prepared.selection_justification)
    if (args() >= 6) S.rdr = cdid_sel_reader_pop(B.prepared.population, Z, hi, (hi ? 1 : 0), w, decl, B.prepared.sel_rho, chunk, keep)
    else S.rdr = cdid_sel_reader_pop(B.prepared.population, Z, hi, (hi ? 1 : 0), w, decl, B.prepared.sel_rho, chunk)
    return(S)
}

/* Unit loadings on the outcome levels of the post comparisons of the same-pair catalogue
   under imputed outcome changes, optionally augmented by inverse observation probabilities
   on the observed pairs (R evaluate_reconstruction_regimes_core with compute_within on the
   same-pair descriptors). Each comparison is the target cohort's mean response minus the
   control response of its pool, reweighted toward the target cohort by a treatment logit
   (ipw), predicted by a counterfactual outcome regression (reg), or both (dr). With dated
   inputs (Sx), every (cohort, recipe) pair is a cohort of its own with the recipe's weights,
   designs and repeat rates; imputation and observation models stay within a recipe and pooled
   counterfactual slopes are shared by every recipe. Rows (date - 1) * n_A + analysis unit, one
   column per post comparison; gvec returns the target cohort and recv the recipe of each post
   comparison, mhat the fitted comparisons (post, then supported auxiliary). */
real matrix cdid_wsr_loadings(struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan,
    real scalar augmented, real scalar eps, real scalar maxit, real scalar chunk, real scalar rank_tol,
    real scalar max_cached, real scalar workspace, real scalar workcap, real colvector gvec,
    real colvector mhat, real colvector recv, | struct cdid_reference_inputs Sx)
{
    struct cdid_population scalar pop
    struct cdid_imputation_model colvector mods
    struct cdid_imputation_pair rowvector pairs
    struct cdid_imputation_budget scalar budget
    struct cdid_xps_sel colvector SS
    struct cdid_sel_log_deriv scalar ld
    struct _cdid_xpr_cell colvector CE
    struct _cdid_wsr_blk colvector BK
    struct _cdid_xpr_grp colvector GR
    struct _cdid_wsr_trt colvector TR
    struct _cdid_xpr_mat colvector Limp_b, Lor_b, WR, ZR, XR, ZS, MZ, PR
    transmorphic CK, BKK, TK
    real matrix X, Zs, pk, req, specs, U, Limp_a, Lor_a, KC, Ptab, Hm, Q, nu, dor, lamM, Xd, Aco, Mone
    real colvector w, calmap, lr, K, codes, kc, rr, wm, mh, pos, xi, wj, cols_a, cols_b, pv, rec, used, inscope, rows_in
    real colvector LG, LS, LT, LNK, bkey, tkey, tk, y, w_syn, beta, hit, crec, goff, cohc, ectx
    real matrix KP
    real rowvector dbp, dbm, lam, sdor
    real scalar N, T, nA, M, Rraw, NL, p, q, icpt, isor, isdr, need_or, need_trt, i, j, m, r, s, t
    real scalar pid, cid, c, ncel, nblk, ng, ntr, ci, gid, k, b, g, cg, ti, tid, sgn, dp, Tg, Bg
    real scalar Cv, Bol, Sm, Sz, Smm, dm, kk, h, L, maxk, isref, R, ncoh, ngrp, gx, dx, rtrt, hi
    string scalar imp_pooling, bk
    string colvector bcells, tcells
    pop = B.prepared.population
    N = pop.N
    T = rows(pop.calendar)
    nA = rows(pop.A_to_F)
    isref = (args() >= 14)
    isor = (B.prepared.method == "reg")
    isdr = (B.prepared.method == "dr")
    need_or = (isor | isdr)
    need_trt = !isor
    imp_pooling = (B.prepared.imputation_pooling == "" ? "separate" : B.prepared.imputation_pooling)
    if (isref) {
        icpt = Sx.roles[1].intercept
        q = cols(Sx.designs[Sx.roles[1].slot].value) - icpt
        rtrt = (Sx.roles[2].active ? 2 : 1)
    }
    else {
        icpt = B.prepared.x_has_intercept
        X = B.prepared.X
        p = cols(X)
        Zs = (icpt ? (p > 1 ? X[, 2..p] : J(N, 0, .)) : X)
        q = cols(Zs)
    }
    if (!icpt & q == 0) _cdid_ipw_fail("The outcome design is empty; add covariates or a constant.", 459)
    calmap = J(rows(plan.calendar), 1, .)
    for (j = 1; j <= rows(plan.calendar); j++) {
        k = selectindex(pop.calendar :== plan.calendar[j])
        if (rows(k) != 1) _cdid_pr_fail("cdid_native_reporting_input", "A catalogue date is absent from the panel calendar.")
        calmap[j] = k[1]
    }
    M = plan.M
    Rraw = plan.M + plan.Q
    lr = (1::M)
    if (Rraw > M) lr = lr \ (selectindex(plan.raw_supported[(M + 1)..Rraw]) :+ M)
    NL = rows(lr)
    rec = (isref ? Sx.view_recipe[lr] : J(NL, 1, 1))
    R = (isref ? max(rec) : 1)
    /* Comparisons: target cohort and sorted control cohorts (code 0 for never treated),
       calendar indices of the pair. */
    LG = LS = LT = LNK = J(NL, 1, .)
    maxk = 0
    for (m = 1; m <= NL; m++) {
        r = lr[m]
        maxk = max((maxk, plan.raw_K_ptr[r + 1] - plan.raw_K_ptr[r]))
    }
    KC = J(NL, max((1, maxk)), .)
    gvec = J(M, 1, .)
    for (m = 1; m <= NL; m++) {
        r = lr[m]
        LG[m] = editmissing(plan.cohorts[plan.raw[r, 1]], 0)
        if (m <= M) gvec[m] = _cdid_wsg_code(plan.cohorts[plan.raw[r, 1]])
        LS[m] = calmap[plan.raw[r, 2]]
        LT[m] = calmap[plan.raw[r, 3]]
        if (LT[m] <= LS[m]) _cdid_pr_fail("cdid_native_reporting_input", "A same-pair comparison is not ordered in time.")
        K = sort(editmissing(plan.cohorts[plan.raw_K_index[plan.raw_K_ptr[r]..(plan.raw_K_ptr[r + 1] - 1)]], 0), 1)
        LNK[m] = rows(K)
        if (rows(K) == 0) _cdid_pr_fail("cdid_native_reporting_input", "A same-pair comparison has no control cohort.")
        KC[m, 1..rows(K)] = K'
    }
    recv = rec[1..M]
    /* Recipe frames: weights and designs on the frame rows of the cohorts each recipe reads. */
    WR = _cdid_xpr_mat(R, 1)
    ZR = _cdid_xpr_mat(R, 1)
    XR = _cdid_xpr_mat(R, 1)
    ZS = _cdid_xpr_mat(R, 1)
    PR = _cdid_xpr_mat(R, 1)
    used = J(R, 1, 0)
    cohc = editmissing(pop.cohort, 0)
    hi = (B.prepared.sel_has_intercept ? 1 : 0)
    /* Units each recipe reads under restricted reference eligibility contexts. */
    KP = J(0, 0, .)
    if (isref) if (cols(Sx.recipes) >= 3 & rows(Sx.elig_index)) {
        KP = J(N, R, 1)
        for (r = 1; r <= R; r++) if (Sx.recipes[r, 3] > 0) {
            ectx = selectindex(Sx.elig_index :== Sx.recipes[r, 3])
            if (length(ectx) != 1) _cdid_pr_fail("cdid_native_reporting_input", "Unknown reference eligibility context.")
            KP[, r] = !Sx.elig_excluded[, ectx[1]]
        }
    }
    for (r = 1; r <= R; r++) {
        inscope = J(N, 1, 0)
        for (m = 1; m <= NL; m++) {
            if (rec[m] != r) continue
            used[r] = 1
            codes = LG[m] \ KC[m, 1..LNK[m]]'
            for (j = 1; j <= rows(codes); j++) inscope = inscope :| (cohc :== codes[j])
        }
        if (!used[r]) continue
        if (rows(KP)) inscope = inscope :& (KP[, r] :!= 0)
        PR[r].v = selectindex(inscope)
        if (!isref) {
            WR[r].v = B.prepared.weight.owner_weight
            ZR[r].v = Zs
            XR[r].v = X
            if (augmented) ZS[r].v = B.prepared.Z_sel
            continue
        }
        rows_in = selectindex(inscope)
        WR[r].v = J(N, 1, 1)
        WR[r].v[rows_in] = (Sx.weight_kind == "dated" ? cdid_ref_weight_rows(Sx, r, rows_in) : B.prepared.weight.owner_weight[rows_in])
        dx = (Sx.roles[1].dated ? Sx.recipes[r, 1] : 0)
        Xd = cdid_ref_design_rows(Sx, 1, dx * (N + 1) :+ rows_in)
        ZR[r].v = J(N, q, 0)
        if (q) ZR[r].v[rows_in, ] = (icpt ? Xd[, 2..cols(Xd)] : Xd)
        if (need_trt) {
            dx = (Sx.roles[rtrt].dated ? Sx.recipes[r, 1] : 0)
            Xd = cdid_ref_design_rows(Sx, rtrt, dx * (N + 1) :+ rows_in)
            XR[r].v = J(N, cols(Xd), 0)
            XR[r].v[rows_in, ] = Xd
        }
        if (augmented) {
            dx = (Sx.roles[3].dated ? Sx.recipes[r, 1] : 0)
            Xd = cdid_ref_design_rows(Sx, 3, dx * (N + 1) :+ rows_in)
            hi = Sx.designs[Sx.roles[3].slot].intercept
            ZS[r].v = J(N, cols(Xd), 0)
            if (hi) ZS[r].v[, 1] = J(N, 1, 1)
            ZS[r].v[rows_in, ] = Xd
        }
        if (missing(WR[r].v[rows_in]) | missing(ZR[r].v[rows_in, ]))
            _cdid_ipw_fail("Covariates and weights must be nonmissing on every frame row of the cohorts compared at a reference date.", 459)
    }
    /* Stage 1: actual-change imputation cells of every comparison, per recipe. */
    mods = cdid_imputation_model(R, 1)
    goff = J(R + 1, 1, 0)
    budget = cdid_imputation_default_budget(max((1, min((N, chunk)))), max_cached, workspace, workcap)
    SS = cdid_xps_sel(R, 1)
    for (r = 1; r <= R; r++) {
        goff[r + 1] = goff[r]
        if (!used[r]) continue
        pk = J(0, 2, .)
        req = J(0, 2, .)
        specs = J(0, 5, .)
        for (m = 1; m <= NL; m++) {
            if (rec[m] != r) continue
            pid = 0
            for (j = 1; j <= rows(pk); j++) {
                if (pk[j, 1] == LS[m] & pk[j, 2] == LT[m]) {
                    pid = j
                    break
                }
            }
            if (!pid) {
                pk = pk \ (LS[m], LT[m])
                pid = rows(pk)
            }
            codes = LG[m] \ KC[m, 1..LNK[m]]'
            for (j = 1; j <= rows(codes); j++) req = req \ (codes[j], pid)
        }
        pairs = cdid_imputation_pair(1, rows(pk))
        for (j = 1; j <= rows(pk); j++) pairs[j] = cdid_imputation_make_pair(pk[j, 1], pk[j, 2], 1)
        if (rows(KP)) mods[r] = cdid_imputation_fit_static(pop, WR[r].v, ZR[r].v, pairs, req, imp_pooling, icpt,
            rank_tol, budget, J(1, 0, ""), KP[, r])
        else mods[r] = cdid_imputation_fit_static(pop, WR[r].v, ZR[r].v, pairs, req, imp_pooling, icpt, rank_tol, budget)
        goff[r + 1] = goff[r] + mods[r].n_groups
        if (augmented) {
            if (rows(KP)) SS[r] = _cdid_wsr_sel_reader(B, ZS[r].v, hi, WR[r].v, chunk, KP[, r])
            else SS[r] = _cdid_wsr_sel_reader(B, ZS[r].v, hi, WR[r].v, chunk)
            for (m = 1; m <= NL; m++) {
                if (rec[m] != r) continue
                codes = LG[m] \ KC[m, 1..LNK[m]]'
                for (j = 1; j <= rows(codes); j++) {
                    specs = specs \ (_cdid_xps_kpos(SS[r].rdr, (codes[j] == 0 ? . : codes[j])), LS[m], LT[m],
                        pop.calendar[LS[m]], pop.calendar[LT[m]])
                }
            }
            _cdid_wsg_sel_fit(SS[r], specs, eps, maxit)
        }
    }
    ngrp = goff[R + 1]
    ncoh = 0
    for (r = 1; r <= R; r++) if (used[r]) ncoh = max((ncoh, mods[r].n_cohorts))
    Aco = J(ncoh, R, 0)
    Mone = J(ncoh, R, 0)
    MZ = _cdid_xpr_mat(R, 1)
    for (r = 1; r <= R; r++) {
        MZ[r].v = J(ncoh, q, 0)
        if (!used[r]) continue
        Aco[1..mods[r].n_cohorts, r] = mods[r].cohort_a_F[1, ]'
        for (ci = 1; ci <= mods[r].n_cohorts; ci++) {
            rr = mods[r].coh_u[ci].rows
            if (rows(rr) == 0) continue
            Mone[ci, r] = sum(WR[r].v[rr] :* Aco[ci, r]) / N
            if (q) MZ[r].v[ci, ] = colsum((WR[r].v[rr] :* Aco[ci, r]) :* ZR[r].v[rr, ]) :/ N
        }
    }
    /* Cells (cohort, pair, recipe) in order of first use. */
    CK = asarray_create("real", 4)
    ncel = 0
    hit = J(0, 1, .)
    crec = J(0, 1, .)
    for (m = 1; m <= NL; m++) {
        r = rec[m]
        codes = LG[m] \ KC[m, 1..LNK[m]]'
        for (j = 1; j <= rows(codes); j++) {
            if (asarray_contains(CK, (codes[j], LS[m], LT[m], r))) continue
            cid = 0
            for (i = 1; i <= mods[r].n_cells; i++) {
                if (mods[r].cells[i].cohort == codes[j] & mods[r].cells[i].s_index == LS[m] & mods[r].cells[i].t_index == LT[m]) {
                    cid = i
                    break
                }
            }
            if (!cid) _cdid_ipw_fail("An imputation cell of the same-pair catalogue was not fitted.", 498)
            ncel++
            asarray(CK, (codes[j], LS[m], LT[m], r), ncel)
            hit = hit \ cid
            crec = crec \ r
        }
    }
    CE = _cdid_xpr_cell(ncel, 1)
    for (c = 1; c <= ncel; c++) {
        r = crec[c]
        cid = hit[c]
        gid = mods[r].cells[cid].group_id
        ci = mods[r].cells[cid].cohort_index
        if (Aco[ci, r] <= 0 | mods[r].cohort_mass_F[1, ci] <= 0)
            _cdid_ipw_fail("Cohort " + mods[r].cohort_labels[ci] + " has no units observed at two or more dates; imputed changes need repeat observers.", 459)
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
        if (augmented & rows(CE[c].obs)) {
            CE[c].kpos = _cdid_xps_kpos(SS[r].rdr, (mods[r].cohorts[ci] == 0 ? . : mods[r].cohorts[ci]))
            CE[c].pc = asarray(SS[r].pcl, (CE[c].kpos, CE[c].cs, CE[c].ct))
            ld = cdid_sel_log_deriv_eval(SS[r].fit, SS[r].rdr, SS[r].um, CE[c].pc, CE[c].rows_obs)
            if (missing(ld.p_weight)) _cdid_ipw_fail("Nonfinite fitted observation probabilities on observed pairs.", 504)
            CE[c].pi_obs = 1 :/ ld.p_weight
            CE[c].Z[CE[c].obs] = CE[c].jhat[CE[c].obs] + CE[c].pi_obs :* CE[c].e_obs
            CE[c].omspi[CE[c].obs] = 1 :- CE[c].pi_obs
        }
    }
    Q = J(0, 0, .)
    /* Counterfactual outcome blocks keyed by the control cells of a comparison. */
    nblk = 0
    ng = 0
    bkey = J(NL, 1, 0)
    if (need_or) {
        BKK = asarray_create("string", 1)
        bcells = J(0, 1, "")
        for (m = 1; m <= NL; m++) {
            kc = J(LNK[m], 1, .)
            for (j = 1; j <= LNK[m]; j++) kc[j] = asarray(CK, (KC[m, j], LS[m], LT[m], rec[m]))
            bk = invtokens(strofreal(kc', "%12.0f"))
            if (!asarray_contains(BKK, bk)) {
                bcells = bcells \ bk
                asarray(BKK, bk, rows(bcells))
            }
            bkey[m] = asarray(BKK, bk)
        }
        nblk = rows(bcells)
        BK = _cdid_wsr_blk(nblk, 1)
        for (b = 1; b <= nblk; b++) {
            BK[b].cells = strtoreal(tokens(bcells[b]))
            BK[b].L = CE[BK[b].cells[1]].ct - CE[BK[b].cells[1]].cs
        }
        if (B.prepared.outcome_pooling == "pooled") {
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
        _cdid_wsr_or_fit(BK, CE, crec, GR, ZR, q, icpt, rank_tol, PR)
        ng = rows(GR)
    }
    /* Treatment logits of the target cohort against its control pool on their frame rows. */
    ntr = 0
    tkey = J(NL, 1, 0)
    if (need_trt) {
        TK = asarray_create("string", 1)
        tcells = J(0, 1, "")
        for (m = 1; m <= NL; m++) {
            r = rec[m]
            cg = asarray(CK, (LG[m], LS[m], LT[m], r))
            tk = J(LNK[m], 1, .)
            for (j = 1; j <= LNK[m]; j++) tk[j] = CE[asarray(CK, (KC[m, j], LS[m], LT[m], r))].ci
            bk = strofreal(r, "%12.0f") + " " + strofreal(CE[cg].ci, "%12.0f") + " " + invtokens(strofreal(sort(tk, 1)', "%12.0f"))
            if (!asarray_contains(TK, bk)) {
                tcells = tcells \ bk
                asarray(TK, bk, rows(tcells))
            }
            tkey[m] = asarray(TK, bk)
        }
        ntr = rows(tcells)
        TR = _cdid_wsr_trt(ntr, 1)
        for (tid = 1; tid <= ntr; tid++) {
            tk = strtoreal(tokens(tcells[tid]))'
            r = tk[1]
            tk = tk[2..rows(tk)]
            rr = J(0, 1, .)
            for (j = 1; j <= rows(tk); j++) rr = rr \ mods[r].coh_u[tk[j]].rows
            TR[tid].rows = sort(rr, 1)
            TR[tid].uci = mods[r].unit_cohort_index[TR[tid].rows]
            y = (TR[tid].uci :== tk[1])
            w_syn = WR[r].v[TR[tid].rows] :* Aco[TR[tid].uci, r]
            _cdid_sel_link_trt_logit(XR[r].v[TR[tid].rows, ], y, w_syn, eps, maxit, beta = ., Hm = .)
            Ptab = _cdid_ipw_logit_values(XR[r].v[TR[tid].rows, ] * beta)
            TR[tid].odds = Ptab[, 2]
        }
    }
    /* Comparison sensitivities. */
    mhat = J(NL, 1, .)
    U = J(T * nA, NL, 0)
    Limp_a = J(ncel, NL, 0)
    Limp_b = _cdid_xpr_mat(max((ngrp, 1)), 1)
    for (r = 1; r <= R; r++) {
        if (!used[r]) continue
        for (g = 1; g <= mods[r].n_groups; g++) Limp_b[goff[r] + g].v = J(cols(mods[r].groups[g].kept), NL, 0)
    }
    Lor_a = J(nblk, NL, 0)
    Lor_b = _cdid_xpr_mat(max((ng, 1)), 1)
    for (g = 1; g <= ng; g++) Lor_b[g].v = J(q, NL, 0)
    sgn = -1
    for (m = 1; m <= NL; m++) {
        r = rec[m]
        ti = asarray(CK, (LG[m], LS[m], LT[m], r))
        kc = J(LNK[m], 1, .)
        for (j = 1; j <= LNK[m]; j++) kc[j] = asarray(CK, (KC[m, j], LS[m], LT[m], r))
        L = LT[m] - LS[m]
        cg = CE[ti].ci
        /* Target arm. */
        dp = Mone[cg, r]
        Tg = sum(CE[ti].w_syn :* CE[ti].Z) / (N * dp)
        wj = CE[ti].w_syn :* CE[ti].omspi :/ (N * dp)
        gx = goff[r] + CE[ti].gid
        if (icpt) Limp_a[ti, m] = Limp_a[ti, m] + sum(wj)
        if (cols(CE[ti].d_imp)) Limp_b[gx].v[, m] = Limp_b[gx].v[, m] + quadcross(CE[ti].d_imp, wj)
        if (augmented & rows(CE[ti].obs)) {
            xi = (CE[ti].w_syn[CE[ti].obs] :/ (N * dp)) :* CE[ti].pi_obs
            _cdid_xpr_addU(U, nA, CE[ti].a_obs, CE[ti].cs, CE[ti].ct, xi, m)
        }
        /* Control arm. */
        Bg = Cv = Bol = 0
        if (need_or) {
            b = bkey[m]
            g = BK[b].grp
            dbp = (q ? (L :* MZ[r].v[cg, ]) :/ dp : J(1, 0, .))
            Bg = (icpt ? BK[b].alpha : 0) + (q ? dbp * GR[g].beta : 0)
        }
        if (need_trt) {
            tid = tkey[m]
            Sm = Sz = Smm = 0
            sdor = J(1, q, 0)
            for (j = 1; j <= LNK[m]; j++) {
                c = kc[j]
                pos = selectindex(TR[tid].uci :== CE[c].ci)
                if (rows(pos) != rows(CE[c].rf)) _cdid_ipw_fail("A control cohort is misaligned with its treatment logit.", 498)
                wm = CE[c].w_syn :* TR[tid].odds[pos]
                Sm = Sm + sum(wm)
                Sz = Sz + sum(wm :* CE[c].Z)
                if (isdr) {
                    dor = (q ? L :* ZR[r].v[CE[c].rf, ] : J(rows(CE[c].rf), 0, .))
                    mh = (icpt ? BK[b].alpha : 0) :+ (q ? dor * GR[g].beta : J(rows(CE[c].rf), 1, 0))
                    Smm = Smm + sum(wm :* mh)
                    if (q) sdor = sdor + colsum(wm :* dor)
                }
            }
            dm = Sm / N
            if (dm <= 0 | dm >= .) _cdid_ipw_fail("The synthetic control normalizer is not positive.", 459)
            Cv = Sz / Sm
            if (isdr) {
                Bol = Smm / Sm
                dbm = (q ? sdor :/ Sm : J(1, 0, .))
            }
            for (j = 1; j <= LNK[m]; j++) {
                c = kc[j]
                pos = selectindex(TR[tid].uci :== CE[c].ci)
                wm = CE[c].w_syn :* TR[tid].odds[pos]
                wj = sgn :* wm :* CE[c].omspi :/ (N * dm)
                gx = goff[r] + CE[c].gid
                if (icpt) Limp_a[c, m] = Limp_a[c, m] + sum(wj)
                if (cols(CE[c].d_imp)) Limp_b[gx].v[, m] = Limp_b[gx].v[, m] + quadcross(CE[c].d_imp, wj)
                if (augmented & rows(CE[c].obs)) {
                    xi = sgn :* (wm[CE[c].obs] :/ (N * dm)) :* CE[c].pi_obs
                    _cdid_xpr_addU(U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, xi, m)
                }
            }
        }
        if (isor) {
            mhat[m] = Tg + sgn * Bg
            if (icpt) Lor_a[b, m] = Lor_a[b, m] + sgn
            if (q) Lor_b[g].v[, m] = Lor_b[g].v[, m] + sgn :* dbp'
        }
        else if (isdr) {
            mhat[m] = Tg + sgn * (Bg + Cv - Bol)
            if (q) Lor_b[g].v[, m] = Lor_b[g].v[, m] + sgn :* (dbp - dbm)'
        }
        else mhat[m] = Tg + sgn * Cv
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
            lam = (icpt ? Lor_a[b, ] : J(1, NL, 0))
            cols_a = selectindex(lam' :!= 0)
            cols_b = (q ? selectindex(colsum(Lor_b[g].v[GR[g].keep, ] :!= 0)' :> 0) : J(0, 1, .))
            if (rows(cols_a) + rows(cols_b) == 0) continue
            if (rows(cols_b)) {
                nu = N :* cholsolve(GR[g].M, Lor_b[g].v[GR[g].keep, cols_b'])
                if (missing(nu)) _cdid_ipw_fail("The counterfactual outcome-regression Gram is singular.", 459)
            }
            for (h = 1; h <= cols(BK[b].cells); h++) {
                c = BK[b].cells[h]
                r = crec[c]
                gx = goff[r] + CE[c].gid
                for (kk = 1; kk <= 2; kk++) {
                    if (kk == 1) {
                        if (rows(cols_a) == 0) continue
                        Q = (CE[c].w_syn :/ BK[b].wsum) * lam[cols_a']
                        pv = cols_a
                    }
                    else {
                        if (rows(cols_b) == 0) continue
                        dor = BK[b].L :* ZR[r].v[CE[c].rf, ]
                        Q = ((CE[c].w_syn :/ N) :* (dor[, GR[g].keep] :- (icpt ? BK[b].mean_d[GR[g].keep] : J(1, cols(GR[g].keep), 0)))) * nu
                        pv = cols_b
                    }
                    lamM = CE[c].omspi :* Q
                    if (icpt) Limp_a[c, pv'] = Limp_a[c, pv'] + colsum(lamM)
                    if (cols(CE[c].d_imp)) Limp_b[gx].v[, pv'] = Limp_b[gx].v[, pv'] + quadcross(CE[c].d_imp, lamM)
                    if (augmented & rows(CE[c].obs))
                        _cdid_xpr_addU(U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, CE[c].pi_obs :* Q[CE[c].obs, ], pv)
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
        lam = (icpt ? Limp_a[c, ] : J(1, NL, 0))
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
            _cdid_xpr_addU(U, nA, CE[c].a_obs, CE[c].cs, CE[c].ct, Q, pv)
        }
    }
    return(U[., 1..M])
}

/* Within-unit weighting of the same-pair catalogue fitted with imputed outcome changes or the
   augmented observation correction (R compute_same_pair_within_bundle with
   evaluate_reconstruction_regimes_core): level loadings, the post-fit exact-mask spanning forest
   (recipe in the group key) and the within covariance of the kept comparisons. */
void cdid_wsr_within(struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan,
    struct cdid_cat scalar C, real colvector mhat, real scalar augmented, real scalar eps, real scalar maxit,
    real scalar chunk, real scalar rank_tol, real scalar max_cached, real scalar workspace, real scalar workcap,
    real colvector kept, real matrix G, real matrix V, | struct cdid_reference_inputs S_ref)
{
    struct cdid_reference_inputs scalar Sx
    real matrix U
    real colvector gvec, mfit, recv
    real scalar M, dev, isref
    isref = (args() >= 16)
    if (isref) {
        Sx = S_ref
        U = cdid_wsr_loadings(B, plan, augmented, eps, maxit, chunk, rank_tol, max_cached, workspace, workcap, gvec, mfit, recv, Sx)
    }
    else U = cdid_wsr_loadings(B, plan, augmented, eps, maxit, chunk, rank_tol, max_cached, workspace, workcap, gvec, mfit, recv)
    M = rows(mhat)
    if (M != cols(U) | M != rows(C.r_g)) _cdid_pr_fail("cdid_native_reporting_input", "The within loadings do not align with the catalogue.")
    dev = max(abs(mfit[1..M] - mhat)) / max((1, max(abs(mhat))))
    if (dev > 1e-8) _cdid_pr_fail("cdid_native_reporting_input", "The within loadings do not reproduce the fitted comparisons.")
    _cdid_wsg_within_finish(B, C, mhat, U, gvec, recv, isref, Sx, kept, G, V)
}

end
