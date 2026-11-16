*! Within-unit GMM weighting for same-pair comparisons fitted with selection-model observation
*! weights or with pooled outcome slopes.
version 19.5
mata:

/* Observed-pair cell of one cohort: analysis units observed at both dates, their frame rows
   and observation weights. */
struct _cdid_wsg_cell {
    real scalar cv, cs, ct, r
    real colvector a, f, fw
    real matrix xt, xs
}

/* Outcome-regression block: the control cells of one comparison pool on one date pair. */
struct _cdid_wsg_blk {
    real scalar L, ws, grp
    real rowvector md
    real rowvector cells
}

struct _cdid_wsg_grp {
    real colvector blocks
    real matrix M
}

real scalar _cdid_wsg_code(real scalar v)
{
    return(v >= . ? -8.98e307 : v)
}

/* Level loadings of a pair difference: + at the later date, - at the earlier date. */
void _cdid_wsg_addU(real matrix U, real scalar nA, struct _cdid_wsg_cell scalar c, real colvector x, real scalar l)
{
    if (rows(c.a) == 0) return
    U[(c.ct - 1) * nA :+ c.a, l] = U[(c.ct - 1) * nA :+ c.a, l] + x
    U[(c.cs - 1) * nA :+ c.a, l] = U[(c.cs - 1) * nA :+ c.a, l] - x
}

/* Selection reader of the bundle. */
struct cdid_xps_sel scalar _cdid_wsg_sel_reader(struct cdid_public_bundle scalar B, real scalar chunk)
{
    struct cdid_xps_sel scalar S
    struct cdid_sel_decl scalar decl
    struct cdid_population scalar pop
    real scalar N, univ_has_rho, int_col
    string scalar univ
    pop = B.prepared.population
    N = pop.N
    univ = (strtrim(B.prepared.selection_universe) == "" ? "F" : strtrim(B.prepared.selection_universe))
    int_col = (B.prepared.sel_has_intercept ? 1 : 0)
    univ_has_rho = (B.prepared.sel_has_rho != 0 & B.prepared.sel_has_rho < . ? 1 :
        (univ != "F" | (rows(B.prepared.sel_rho) == N & any(B.prepared.sel_rho :!= 1))))
    decl = cdid_sel_norm_decl("conditional_mtar", 1, univ, univ_has_rho,
        (strtrim(B.prepared.selection_justification) != ""), B.prepared.selection_model,
        B.prepared.selection_pooling, B.prepared.selection_justification)
    S.rdr = cdid_sel_reader_pop(pop, B.prepared.Z_sel, B.prepared.sel_has_intercept, int_col,
        B.prepared.weight.owner_weight, decl, B.prepared.sel_rho, chunk)
    return(S)
}

/* Selection models fitted on the observed-pair cells of the given specifications (cohort
   position in the selection reader, calendar indices and dates of the pair). */
void _cdid_wsg_sel_fit(struct cdid_xps_sel scalar S, real matrix specs, real scalar eps, real scalar maxit)
{
    real scalar i
    S.fit = cdid_sel_fit(S.rdr, cdid_sel_partition(S.rdr, specs), eps, maxit)
    S.um = cdid_sel_up_map_build(S.fit)
    S.pcl = asarray_create("real", 3)
    for (i = 1; i <= rows(S.fit.map.pair_cohorts); i++) {
        asarray(S.pcl, S.fit.map.pair_cohorts[i, 1..3], i)
    }
}

/* Unit weights of frame units for recipe r: the anchored weights of a recipe when the weights
   vary over time, the unit weights otherwise. */
real colvector _cdid_wsg_weights(struct cdid_public_bundle scalar B, real scalar r, real colvector units,
    real scalar isref, struct cdid_reference_inputs scalar S)
{
    if (rows(units) == 0) return(J(0, 1, .))
    if (isref) if (S.weight_kind == "dated") return(cdid_ref_weight_rows(S, r, units))
    return(B.prepared.weight.owner_weight[units])
}

/* Support of the fitted selection models of a recipe provider for cohort cv, the date pair
   (sv, tv) and reference recipe r. */
real scalar _cdid_wsg_selref_support(struct cdid_selref_provider scalar P, struct cdid_baseline_plan_state scalar plan,
    real scalar cv, real scalar sv, real scalar tv, real scalar r)
{
    real colvector k
    real scalar mr, j, pid
    if (r < 1 | r > rows(P.reference_recipe)) _cdid_pr_fail("cdid_native_reporting_input", "A comparison design block is absent from the selection models.")
    mr = P.reference_recipe[r]
    k = (cv >= . ? selectindex(plan.cohorts :>= .) : selectindex(plan.cohorts :== cv))
    if (length(k) != 1 | mr < 1) _cdid_pr_fail("cdid_native_reporting_input", "A comparison cell is absent from the selection models.")
    pid = 0
    for (j = 1; j <= rows(P.model.pairs); j++) {
        if (P.model.pairs[j].s_val == sv & P.model.pairs[j].t_val == tv) {
            pid = j
            break
        }
    }
    for (j = 1; j <= P.model.n_supports; j++) {
        if (P.model.supports[j].pair_id == pid & P.model.supports[j].cohort_pos == k[1] & P.model.supports[j].recipe_id == mr) return(j)
    }
    _cdid_pr_fail("cdid_native_reporting_input", "A comparison cell is absent from the selection models.")
    return(0)
}

/* Selection weights of the frame units f (observed at both dates) on support sid. */
real colvector _cdid_wsg_selref_weights(struct cdid_selref_provider scalar P, real scalar sid, real colvector f)
{
    struct cdid_selected_rows scalar R
    if (rows(f) == 0) return(J(0, 1, .))
    R = _cdid_selected_support_rows(P.model, sid, f)
    if (rows(R.unit) != rows(f)) _cdid_pr_fail("cdid_native_reporting_input", "The selection weights do not cover the observed pairs of a comparison cell.")
    if (any(R.unit :!= f)) _cdid_pr_fail("cdid_native_reporting_input", "The selection weights do not cover the observed pairs of a comparison cell.")
    return(R.q)
}

/* Unit loadings on the outcome levels of the post comparisons of the catalogue (R
   evaluate_weighted_regimes_core with compute_within on the same-pair descriptors): target
   change weights, control odds weights (ipw, dr) and the outcome-regression intercept and
   slope adjoints (reg, dr). Rows (date - 1) * n_A + analysis unit, one column per comparison;
   gvec returns the target cohort of each comparison. */
real matrix cdid_wsg_loadings(struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan,
    real scalar sel, real scalar eps, real scalar maxit, real scalar chunk, real colvector gvec,
    real colvector recv, | struct cdid_reference_inputs S_ref, pointer(struct cdid_selref_provider scalar) scalar SP)
{
    struct cdid_population scalar pop
    struct cdid_xps_sel scalar S
    struct _cdid_wsg_cell colvector CE
    struct cdid_reference_inputs scalar Sx
    real scalar isref, rtrt, rx, dx, dtt
    struct _cdid_wsg_blk colvector BK
    struct _cdid_wsg_grp colvector GR
    struct cdid_sel_pred scalar pr
    transmorphic CK, BKK
    real matrix ckm, grid, Y, O, U, Xt, Xs, specs, X_all, prob, Lam, Pm, ldum, Xc, H, Rf, nu, dt, Wint, d_c
    real colvector wA, lr, aF, cohA, wF, k, calmap, K, kc, rA, beta_t, y_all, w_all, fw_c, om_c, odds_c, pos, lam, bkey
    real rowvector dbp, dbm, tau, pv
    real scalar N, T, nA, M, NL, r, p, q, t, i, j, m, s, cv, ti, c, b, gr, nb, ng, need_or, need_trt, isor, isdr
    real scalar Sg, Rraw, L, mx, rk, Mk, rraw, rc, h, selref, sid
    real rowvector ectx
    string scalar pooling, bk
    string colvector bcells
    pop = B.prepared.population
    N = pop.N
    T = rows(pop.calendar)
    aF = pop.A_to_F
    nA = rows(aF)
    if (rows(pop.obs_grid) == N & cols(pop.obs_grid) == T) grid = pop.obs_grid[aF, .]
    else {
        grid = J(nA, T, 0)
        for (t = 1; t <= T; t++) grid[, t] = cdid_population_rows(pop, aF, J(nA, 1, pop.calendar[t]))
    }
    Y = J(nA, T, .)
    O = J(nA, T, 0)
    for (t = 1; t <= T; t++) {
        k = selectindex(grid[, t] :> 0)
        if (rows(k) == 0) continue
        Y[k, t] = pop.y[grid[k, t]]
        O[k, t] = (Y[k, t] :< .)
    }
    grid = J(0, 0, .)
    cohA = pop.cohort[aF]
    wF = B.prepared.weight.owner_weight
    p = cols(B.prepared.X)
    isref = (args() >= 9)
    if (isref) Sx = S_ref
    /* Selection weights of dated inputs come from the fitted recipe models of the provider. */
    selref = (sel & args() >= 10)
    if (sel & isref & !selref) _cdid_pr_fail("cdid_native_reporting_input", "Dated selection weights need the fitted selection models.")
    if (isref) {
        /* Dated inputs: outcome (role 1) and treatment (role 2, else role 1) designs at the
           reference date of each comparison's recipe; intercept in the first column. */
        rx = Sx.roles[1].slot
        rtrt = (Sx.roles[2].active ? 2 : 1)
        if (!Sx.roles[1].intercept) _cdid_pr_fail("cdid_chained_public_control_error",
            "gmmweighting(within) with covariates or weights that vary over time requires a constant in the covariate design.")
        q = cols(Sx.designs[rx].value) - 1
        Xt = Xs = J(0, 0, .)
    }
    else {
        Xt = (p ? B.prepared.X : J(N, 1, 1))
        Xs = (p > 1 ? B.prepared.X[, 2..p] : J(N, 0, .))
        q = cols(Xs)
    }
    isor = (B.prepared.method == "reg")
    isdr = (B.prepared.method == "dr")
    need_or = (isor | isdr)
    need_trt = !isor
    pooling = (B.prepared.outcome_pooling == "pooled" ? "pooled" : "separate")
    calmap = J(rows(plan.calendar), 1, .)
    for (j = 1; j <= rows(plan.calendar); j++) {
        k = selectindex(pop.calendar :== plan.calendar[j])
        if (rows(k) != 1) _cdid_pr_fail("cdid_native_reporting_input", "A catalogue date is absent from the panel calendar.")
        calmap[j] = k[1]
    }
    M = plan.M
    Rraw = plan.M + plan.Q
    /* Comparisons: the post catalogue, then the supported auxiliary comparisons, which enter
       the pooled outcome fit and the selection models. */
    lr = (1::M)
    if (Rraw > M) lr = lr \ (selectindex(plan.raw_supported[(M + 1)..Rraw]) :+ M)
    NL = rows(lr)
    recv = (isref ? Sx.view_recipe[lr] : J(NL, 1, 1))
    /* Selection models over the cells of every supported comparison (post and auxiliary). */
    if (sel & !selref) {
        S = _cdid_wsg_sel_reader(B, chunk)
        specs = J(0, 5, .)
        for (i = 1; i <= Rraw; i++) {
            if (i > M & !plan.raw_supported[i]) continue
            K = plan.cohorts[plan.raw[i, 1]] \ plan.cohorts[plan.raw_K_index[plan.raw_K_ptr[i]..(plan.raw_K_ptr[i + 1] - 1)]]
            for (j = 1; j <= rows(K); j++) {
                specs = specs \ (_cdid_xps_kpos(S.rdr, K[j]), calmap[plan.raw[i, 2]], calmap[plan.raw[i, 3]],
                    plan.calendar[plan.raw[i, 2]], plan.calendar[plan.raw[i, 3]])
            }
        }
        _cdid_wsg_sel_fit(S, specs, eps, maxit)
    }
    /* Cells of the comparisons: cohort, date pair and recipe. */
    CK = asarray_create("real", 4)
    ckm = J(0, 4, .)
    gvec = J(M, 1, .)
    for (m = 1; m <= NL; m++) {
        r = lr[m]
        K = plan.cohorts[plan.raw[r, 1]] \ plan.cohorts[plan.raw_K_index[plan.raw_K_ptr[r]..(plan.raw_K_ptr[r + 1] - 1)]]
        if (m <= M) gvec[m] = _cdid_wsg_code(K[1])
        s = calmap[plan.raw[r, 2]]
        t = calmap[plan.raw[r, 3]]
        for (j = 1; j <= rows(K); j++) {
            if (asarray_contains(CK, (_cdid_wsg_code(K[j]), s, t, recv[m]))) continue
            ckm = ckm \ (K[j], s, t, recv[m])
            asarray(CK, (_cdid_wsg_code(K[j]), s, t, recv[m]), rows(ckm))
        }
    }
    CE = _cdid_wsg_cell(rows(ckm), 1)
    for (c = 1; c <= rows(ckm); c++) {
        cv = ckm[c, 1]
        s = ckm[c, 2]
        t = ckm[c, 3]
        CE[c].cv = cv
        CE[c].cs = s
        CE[c].ct = t
        CE[c].r = ckm[c, 4]
        rA = (cv >= . ? selectindex(cohA :>= .) : selectindex(cohA :== cv))
        /* Restricted eligibility context of the recipe (missingreference(exclude)). */
        if (isref) if (cols(Sx.recipes) >= 3 & rows(rA)) if (Sx.recipes[CE[c].r, 3] > 0) {
            ectx = selectindex(Sx.elig_index :== Sx.recipes[CE[c].r, 3])
            if (length(ectx) != 1) _cdid_pr_fail("cdid_native_reporting_input", "Unknown reference eligibility context.")
            rA = select(rA, !Sx.elig_excluded[aF[rA], ectx[1]])
        }
        CE[c].a = (rows(rA) ? select(rA, O[rA, s] :& O[rA, t]) : J(0, 1, .))
        CE[c].f = aF[CE[c].a]
        if (isref) {
            dx = (Sx.roles[1].dated ? Sx.recipes[CE[c].r, 1] : 0)
            dtt = (Sx.roles[rtrt].dated ? Sx.recipes[CE[c].r, 1] : 0)
            CE[c].xs = cdid_ref_design_rows(Sx, 1, dx * (N + 1) :+ CE[c].f)
            CE[c].xs = (q ? CE[c].xs[., 2..(q + 1)] : J(rows(CE[c].f), 0, .))
            if (need_trt) CE[c].xt = cdid_ref_design_rows(Sx, rtrt, dtt * (N + 1) :+ CE[c].f)
        }
        else {
            CE[c].xs = Xs[CE[c].f, ]
            if (need_trt) CE[c].xt = Xt[CE[c].f, ]
        }
        if (selref) {
            sid = _cdid_wsg_selref_support(*SP, plan, cv, pop.calendar[s], pop.calendar[t], CE[c].r)
            CE[c].fw = _cdid_wsg_selref_weights(*SP, sid, CE[c].f)
        }
        else if (sel) {
            pr = cdid_sel_pred_pair(S.fit, S.rdr, asarray(S.pcl, (_cdid_xps_kpos(S.rdr, cv), s, t)), CE[c].f)
            CE[c].fw = pr.sel_weight
        }
        else {
            wA = _cdid_wsg_weights(B, CE[c].r, aF[rA], isref, Sx)
            Mk = sum(wA)
            rraw = (rows(CE[c].a) < rows(rA) ? sum(_cdid_wsg_weights(B, CE[c].r, CE[c].f, isref, Sx)) / Mk : 1)
            rc = max((rraw, 1e-6))
            CE[c].fw = _cdid_wsg_weights(B, CE[c].r, CE[c].f, isref, Sx) :/ rc
        }
    }
    U = J(T * nA, NL, 0)
    /* Outcome-regression blocks, keyed by control cells and date pair. */
    nb = 0
    bkey = J(NL, 1, 0)
    if (need_or) {
        BKK = asarray_create("string", 1)
        bcells = J(0, 1, "")
        for (m = 1; m <= NL; m++) {
            r = lr[m]
            K = plan.cohorts[plan.raw_K_index[plan.raw_K_ptr[r]..(plan.raw_K_ptr[r + 1] - 1)]]
            s = calmap[plan.raw[r, 2]]
            t = calmap[plan.raw[r, 3]]
            kc = J(rows(K), 1, .)
            for (j = 1; j <= rows(K); j++) kc[j] = asarray(CK, (_cdid_wsg_code(K[j]), s, t, recv[m]))
            bk = invtokens(strofreal(kc', "%12.0f"))
            if (!asarray_contains(BKK, bk)) {
                bcells = bcells \ bk
                asarray(BKK, bk, rows(bcells))
            }
            bkey[m] = asarray(BKK, bk)
        }
        BK = _cdid_wsg_blk(rows(bcells), 1)
        for (b = 1; b <= rows(bcells); b++) {
            kc = strtoreal(tokens(bcells[b]))'
            BK[b].cells = kc'
            BK[b].L = CE[kc[1]].ct - CE[kc[1]].cs
            BK[b].ws = 0
            BK[b].md = J(1, q, 0)
            for (j = 1; j <= rows(kc); j++) {
                BK[b].ws = BK[b].ws + sum(CE[kc[j]].fw)
                if (q) BK[b].md = BK[b].md + colsum(CE[kc[j]].fw :* (BK[b].L :* CE[kc[j]].xs))
            }
            if (BK[b].ws <= 0) _cdid_ipw_fail("A control pool has no observed outcome pairs.", 459)
            BK[b].md = BK[b].md :/ BK[b].ws
        }
        nb = rows(BK)
        if (pooling == "pooled") {
            GR = _cdid_wsg_grp(1, 1)
            GR[1].blocks = (1::nb)
            for (b = 1; b <= nb; b++) BK[b].grp = 1
        }
        else {
            GR = _cdid_wsg_grp(nb, 1)
            for (b = 1; b <= nb; b++) {
                GR[b].blocks = b
                BK[b].grp = b
            }
        }
        ng = rows(GR)
        if (q) {
            for (gr = 1; gr <= ng; gr++) {
                Xc = J(0, q, .)
                for (h = 1; h <= rows(GR[gr].blocks); h++) {
                    b = GR[gr].blocks[h]
                    for (j = 1; j <= cols(BK[b].cells); j++) {
                        c = BK[b].cells[j]
                        Xc = Xc \ ((BK[b].L :* CE[c].xs) :- BK[b].md) :* sqrt(CE[c].fw)
                    }
                }
                hqrdp(Xc, H = ., tau = ., Rf = ., pv = .)
                mx = max(abs(diagonal(Rf)))
                rk = sum(abs(diagonal(Rf)) :> 1e-7 * mx)
                if (rk < q | mx <= 0) _cdid_ipw_fail("The outcome-regression design is rank deficient on the observed pairs of a control pool; remove the redundant covariate.", 459)
                GR[gr].M = quadcross(Xc, Xc)
                Xc = H = J(0, 0, .)
            }
        }
        Lam = J(nb, NL, 0)
        Pm = J(ng * q, NL, 0)
    }
    /* Comparison loadings. */
    for (m = 1; m <= NL; m++) {
        r = lr[m]
        K = plan.cohorts[plan.raw_K_index[plan.raw_K_ptr[r]..(plan.raw_K_ptr[r + 1] - 1)]]
        s = calmap[plan.raw[r, 2]]
        t = calmap[plan.raw[r, 3]]
        L = t - s
        ti = asarray(CK, (_cdid_wsg_code(plan.cohorts[plan.raw[r, 1]]), s, t, recv[m]))
        kc = J(rows(K), 1, .)
        for (j = 1; j <= rows(K); j++) kc[j] = asarray(CK, (_cdid_wsg_code(K[j]), s, t, recv[m]))
        Sg = sum(CE[ti].fw)
        if (Sg <= 0) _cdid_ipw_fail("A target cohort has no observed outcome pair for a comparison.", 459)
        _cdid_wsg_addU(U, nA, CE[ti], CE[ti].fw :/ Sg, m)
        if (need_trt) {
            X_all = CE[ti].xt
            y_all = J(rows(CE[ti].f), 1, 1)
            w_all = CE[ti].fw
            for (j = 1; j <= rows(kc); j++) {
                X_all = X_all \ CE[kc[j]].xt
                y_all = y_all \ J(rows(CE[kc[j]].f), 1, 0)
                w_all = w_all \ CE[kc[j]].fw
            }
            _cdid_sel_link_trt_logit(X_all, y_all, w_all, eps, maxit, beta_t = ., ldum = .)
            prob = _cdid_ipw_logit_values(X_all * beta_t)
            pos = (rows(CE[ti].f) + 1)::rows(X_all)
            odds_c = prob[pos, 2]
            fw_c = w_all[pos]
            om_c = (fw_c :* odds_c)
            om_c = om_c :/ sum(om_c)
            h = 0
            for (j = 1; j <= rows(kc); j++) {
                c = kc[j]
                if (rows(CE[c].a)) _cdid_wsg_addU(U, nA, CE[c], -om_c[(h + 1)..(h + rows(CE[c].a))], m)
                h = h + rows(CE[c].a)
            }
        }
        if (need_or) {
            b = bkey[m]
            gr = BK[b].grp
            if (isor) Lam[b, m] = Lam[b, m] - 1
            if (q) {
                dbp = colsum(CE[ti].fw :* (L :* CE[ti].xs)) :/ Sg
                if (isdr) {
                    d_c = J(0, q, .)
                    for (j = 1; j <= rows(kc); j++) d_c = d_c \ (L :* CE[kc[j]].xs)
                    dbm = colsum(om_c :* d_c)
                    Pm[((gr - 1) * q + 1)..(gr * q), m] = Pm[((gr - 1) * q + 1)..(gr * q), m] - (dbp - dbm)'
                }
                else Pm[((gr - 1) * q + 1)..(gr * q), m] = Pm[((gr - 1) * q + 1)..(gr * q), m] - dbp'
            }
        }
    }
    /* Outcome-regression adjoints on the control cells. */
    if (need_or) {
        if (q) {
            for (b = 1; b <= nb; b++) {
                if (any(Lam[b, ] :!= 0)) {
                    gr = BK[b].grp
                    Pm[((gr - 1) * q + 1)..(gr * q), ] = Pm[((gr - 1) * q + 1)..(gr * q), ] - BK[b].md' * Lam[b, ]
                }
            }
        }
        for (b = 1; b <= nb; b++) {
            lam = Lam[b, ]'
            pos = selectindex(lam :!= 0)
            if (rows(pos)) {
                for (j = 1; j <= cols(BK[b].cells); j++) {
                    c = BK[b].cells[j]
                    Wint = (CE[c].fw :/ BK[b].ws) * lam[pos]'
                    for (h = 1; h <= rows(pos); h++) _cdid_wsg_addU(U, nA, CE[c], Wint[, h], pos[h])
                }
            }
        }
        if (q) {
            for (gr = 1; gr <= ng; gr++) {
                nu = Pm[((gr - 1) * q + 1)..(gr * q), ]
                pos = selectindex(colsum(abs(nu))' :> 0)
                if (rows(pos) == 0) continue
                nu = N :* cholsolve(GR[gr].M, nu[, pos'])
                if (missing(nu)) _cdid_ipw_fail("The outcome-regression Gram matrix is singular.", 459)
                for (h = 1; h <= rows(GR[gr].blocks); h++) {
                    b = GR[gr].blocks[h]
                    for (j = 1; j <= cols(BK[b].cells); j++) {
                        c = BK[b].cells[j]
                        dt = (BK[b].L :* CE[c].xs) :- BK[b].md
                        Wint = ((CE[c].fw :/ N) :* dt) * nu
                        for (i = 1; i <= rows(pos); i++) _cdid_wsg_addU(U, nA, CE[c], Wint[, i], pos[i])
                    }
                }
            }
        }
    }
    return(U[., 1..M])
}

/* Post-fit exact-mask spanning forest (R verify_post_fit_forest_reduction): within each
   group (C.gi, C.ci, C.mid) a chord is dropped only when its fitted moment and its
   unit influence function telescope along the tree path (relative tolerance 1e-9);
   kept links ascending, metric G = I + E E' on the verified chords. */
void _cdid_wsg_forest(struct cdid_xp_cat scalar C, real colvector mhat, real matrix Phi,
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
                if (rows(surv) == 0) break
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

/* Within scale of the analysis units for the comparisons of target cohort g (R
   build_within_scale_by_g): sqrt(sum w / w_i) with the unit weights, or, when the weights vary
   over time without a first-period reduction, sqrt(n_A / rel_i) with the weights at the anchor of
   the cohort's first comparison relative to their mean (other anchors fill units without one). */
real matrix _cdid_wsg_scale(struct cdid_public_bundle scalar B, real colvector gvec, real colvector recv,
    real scalar isref, struct cdid_reference_inputs scalar S)
{
    real colvector aF, w, ug, wi, rel, fb, ok, rg, miss, hit, wsup
    real matrix SC
    real scalar nA, i, j, idx, N1, shared, r
    aF = B.prepared.population.A_to_F
    nA = rows(aF)
    ug = uniqrows(gvec)
    SC = J(nA, rows(ug), 0)
    /* Restricted eligibility contexts (missingreference(exclude)): per target cohort, the
       eligible units' weights relative to their mean, with the common normalizer n_A. */
    if (isref) if (rows(S.elig_index)) {
        shared = (S.weight_kind != "dated")
        if (!shared) shared = (S.fix_weights == "first_period" & rows(S.V) == S.N & !missing(S.V) & all(S.V[aF] :> 0))
        N1 = S.N + 1
        for (i = 1; i <= rows(ug); i++) {
            r = recv[selectindex(gvec :== ug[i])[1]]
            if (shared) w = S.V[aF]
            else {
                hit = cdid_ref_lookup(S.recipes[r, 2] * N1 :+ aF, S.weight_keys)
                w = J(nA, 1, .)
                if (any(hit :> 0)) w[selectindex(hit :> 0)] = S.weights[hit[selectindex(hit :> 0)]]
            }
            ok = (w :< .) :& (w :> 0)
            if (cols(S.recipes) >= 3) if (S.recipes[r, 3] > 0) {
                idx = selectindex(S.elig_index :== S.recipes[r, 3])[1]
                ok = ok :& !S.elig_excluded[aF, idx]
            }
            if (any(ok)) {
                wsup = select(w, ok)
                SC[selectindex(ok), i] = sqrt(nA :* mean(wsup) :/ wsup)
            }
        }
        return(SC)
    }
    if (isref) if (S.weight_kind == "dated" & !(S.fix_weights == "first_period" & rows(S.V) == S.N & !missing(S.V) & all(S.V :> 0))) {
        wi = uniqrows(S.recipes[recv, 2])
        rel = J(nA, rows(wi), 0)
        fb = J(nA, 1, 0)
        N1 = S.N + 1
        for (j = 1; j <= rows(wi); j++) {
            hit = cdid_ref_lookup(wi[j] * N1 :+ aF, S.weight_keys)
            ok = (hit :> 0)
            wsup = J(nA, 1, 0)
            if (any(ok)) wsup[selectindex(ok)] = S.weights[hit[selectindex(ok)]]
            rel[., j] = (any(ok) ? wsup :/ mean(select(wsup, ok)) : wsup)
            rel[., j] = rel[., j] :* (rel[., j] :> 0)
            miss = (fb :<= 0) :& (rel[., j] :> 0)
            if (any(miss)) fb[selectindex(miss)] = rel[selectindex(miss), j]
        }
        for (i = 1; i <= rows(ug); i++) {
            idx = selectindex(wi :== S.recipes[recv[selectindex(gvec :== ug[i])[1]], 2])
            rg = rel[., idx]
            miss = (rg :<= 0) :& (fb :> 0)
            if (any(miss)) rg[selectindex(miss)] = fb[selectindex(miss)]
            SC[., i] = (rg :> 0) :* sqrt(nA :/ (rg + (rg :<= 0)))
        }
        return(SC)
    }
    if (isref) w = (S.weight_kind == "dated" ? S.V[aF] : B.prepared.weight.owner_weight[aF])
    else w = B.prepared.weight.owner_weight[aF]
    for (i = 1; i <= rows(ug); i++) SC[., i] = (w :> 0) :* sqrt(sum(w) :/ (w + (w :<= 0)))
    return(SC)
}

/* Within-unit covariance from scaled level loadings: per target cohort, the sum over dates of
   the cross products of the loadings times the cohort's unit scale. */
real matrix _cdid_wsg_cov(real matrix U, real scalar nA, real colvector gvec, real matrix SC)
{
    real scalar L, T, tau, i
    real colvector ug, cols_g
    real matrix V, Ut, Sg
    L = cols(U)
    T = rows(U) / nA
    V = J(L, L, 0)
    ug = uniqrows(gvec)
    for (i = 1; i <= rows(ug); i++) {
        cols_g = selectindex(gvec :== ug[i])
        Sg = J(rows(cols_g), rows(cols_g), 0)
        for (tau = 1; tau <= T; tau++) {
            Ut = U[((tau - 1) * nA + 1)..(tau * nA), cols_g'] :* SC[., i]
            Sg = Sg + quadcross(Ut, Ut)
        }
        V[cols_g, cols_g] = 0.5 :* (Sg + Sg')
    }
    return(V)
}

/* Within-unit weighting of the same-pair catalogue fitted with selection-model weights, pooled
   outcome slopes or dated inputs (R compute_same_pair_within_bundle): level loadings, the
   exact-mask spanning forest (target cohort, control pool, units observed at both dates, recipe)
   verified on the fitted moments and loadings, and the within covariance of the kept
   comparisons. */
void cdid_wsg_within(struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan,
    struct cdid_cat scalar C, real colvector mhat, real scalar sel, real scalar eps, real scalar maxit,
    real scalar chunk, real colvector kept, real matrix G, real matrix V,
    | struct cdid_reference_inputs S_ref, pointer(struct cdid_selref_provider scalar) scalar SP)
{
    struct cdid_reference_inputs scalar Sx
    real matrix U
    real colvector gvec, recv
    real scalar M, isref
    isref = (args() >= 12)
    if (isref) Sx = S_ref
    if (!isref & !B.prepared.x_has_intercept & cols(B.prepared.X) > 0)
        _cdid_pr_fail("cdid_chained_public_control_error",
            "gmmweighting(within) with selection() or outcomepooling(pooled) requires a constant in the covariate design.")
    if (args() >= 13) U = cdid_wsg_loadings(B, plan, sel, eps, maxit, chunk, gvec, recv, Sx, SP)
    else if (isref) U = cdid_wsg_loadings(B, plan, sel, eps, maxit, chunk, gvec, recv, Sx)
    else U = cdid_wsg_loadings(B, plan, sel, eps, maxit, chunk, gvec, recv)
    M = rows(mhat)
    if (M != cols(U) | M != rows(C.r_g)) _cdid_pr_fail("cdid_native_reporting_input", "The within loadings do not align with the catalogue.")
    recv = recv[1..M]
    _cdid_wsg_within_finish(B, C, mhat, U, gvec, recv, isref, Sx, kept, G, V)
}

/* Exact-mask spanning forest (target cohort, control pool, units observed at both dates,
   recipe) verified on the fitted moments and the level loadings U, and the within covariance
   of the kept comparisons. */
void _cdid_wsg_within_finish(struct cdid_public_bundle scalar B, struct cdid_cat scalar C,
    real colvector mhat, real matrix U, real colvector gvec, real colvector recv, real scalar isref,
    struct cdid_reference_inputs scalar Sx, real colvector kept, real matrix G, real matrix V)
{
    struct _cdid_wsp_data scalar d
    struct cdid_xp_cat scalar X
    real matrix SC, fp
    real colvector key, ord, kg
    real scalar nA, md, M, e, ng
    M = rows(mhat)
    nA = rows(B.prepared.population.A_to_F)
    d = _cdid_wsp_prepare(C, B.prepared.population, (isref ? J(B.prepared.population.N, 1, 1) : B.prepared.X), B.prepared.weight.owner_weight)
    key = _cdid_wsp_key(d, C)
    /* Recipe enters the group key. */
    fp = (key, recv)
    ord = order(fp, (1, 2))
    kg = J(M, 1, 0)
    ng = 0
    for (e = 1; e <= M; e++) {
        if (e == 1) ng = 1
        else if (any(fp[ord[e], ] :!= fp[ord[e - 1], ])) ng++
        kg[ord[e]] = ng
    }
    X.g = C.r_g
    X.gi = kg
    X.ci = J(M, 1, 0)
    X.mid = J(M, 1, 0)
    X.cs = C.r_s
    X.ct = C.r_t
    X.cal = C.calendar
    _cdid_wsg_forest(X, mhat, J(0, 0, .), U, nA, kept, G, md)
    SC = _cdid_wsg_scale(B, gvec, recv, isref, Sx)
    V = _cdid_wsg_cov(U[, kept], nA, gvec[kept], SC[., _cdid_wsg_gpos(gvec, gvec[kept])])
}

/* Positions of the target cohorts of the kept comparisons among the sorted cohorts of all. */
real rowvector _cdid_wsg_gpos(real colvector gall, real colvector gkept)
{
    real colvector ua, uk
    real rowvector pos
    real scalar i
    ua = uniqrows(gall)
    uk = uniqrows(gkept)
    pos = J(1, rows(uk), .)
    for (i = 1; i <= rows(uk); i++) pos[i] = selectindex(ua :== uk[i])[1]
    return(pos)
}

end
