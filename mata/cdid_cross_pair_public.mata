*! Cross-pair chained DiD (crosspair) from the prepared long data to the reported cells.
version 19.5
mata:

/* Request of the cross-pair route: parallel-trends assumption and its start date. */
void cdid_xp_request(real scalar on, string scalar pt, real scalar pt_start, | string scalar wname, string scalar cg)
{
    pointer() scalar p
    if (!on) {
        if (findexternal("__cdid_xp_req") != NULL) rmexternal("__cdid_xp_req")
        return
    }
    p = findexternal("__cdid_xp_req")
    if (p == NULL) p = crexternal("__cdid_xp_req")
    *p = (pt, strofreal(pt_start, "%21.17g"), (args() < 4 ? "" : wname), (args() < 5 ? "" : cg))
}

/* Control group of the cross-pair request (never-treated unless not-yet-treated was set). */
string scalar cdid_xp_control_group()
{
    pointer() scalar p
    p = findexternal("__cdid_xp_req")
    if (p == NULL) return("nevertreated")
    if (cols(*p) < 4) return("nevertreated")
    return((*p)[4] == "notyettreated" ? "notyettreated" : "nevertreated")
}

real scalar cdid_xp_requested()
{
    return(findexternal("__cdid_xp_req") != NULL)
}

/* Locate the engine cell of (g, t); 0 if none. */
real scalar _cdid_xp_cell(struct cdid_xp_fit scalar F, real scalar g, real scalar t)
{
    real colvector h
    h = selectindex((F.cells[, 1] :== g) :& (F.cells[, 2] :== t))
    return(rows(h) ? h[1] : 0)
}

/* Raw values of active variables at calendar index d (0: unit attributes) for the given
   frame units; missing where the unit has no row at that date. */
real matrix _cdid_xp_raw(struct cdid_reference_inputs scalar S, struct cdid_population scalar P,
    real scalar d, real colvector units, real rowvector vars)
{
    real matrix out
    real colvector rr, ok
    real scalar j, v
    out = J(rows(units), cols(vars), .)
    for (j = 1; j <= cols(vars); j++) {
        v = vars[j]
        if (S.kind[v] == 3) {
            if (d > 0) out[, j] = J(rows(units), 1, P.calendar[d])
        }
        else if (S.kind[v] == 1) out[, j] = S.attributes[units, v]
        else if (d > 0) {
            rr = cdid_population_rows(P, units, J(rows(units), 1, P.calendar[d]))
            ok = selectindex(rr :> 0)
            if (rows(ok)) out[ok, j] = S.data[rr[ok], v]
        }
    }
    return(out)
}

/* Design of a formula role at calendar index d (0: unit attributes) over all frame units:
   raw values for a plain formula, otherwise the encoded design of a role constant within
   units. Returns the full design (with its constant column first when it has one). */
real matrix _cdid_xp_role_design(struct cdid_public_bundle scalar B, struct cdid_population scalar P,
    real scalar role, real scalar d, real scalar icpt)
{
    real rowvector v
    real matrix raw
    real scalar sl, j, pl
    real colvector h
    string rowvector tk
    v = B.reference.roles[role].raw_index
    if (cols(v) == 0) {
        /* A role constant within units lists no raw variables; read them by name. */
        tk = tokens(B.reference.roles[role].formula)
        v = J(1, cols(tk), 0)
        for (j = 1; j <= cols(tk); j++) {
            h = selectindex(B.reference.raw_names' :== tk[j])
            if (length(h) == 1) v[j] = h[1]
        }
        if (cols(v)) {
            if (any(v :== 0)) v = J(1, 0, 0)
            else if (any(B.reference.kind[v] :!= 1)) v = J(1, 0, 0)
        }
    }
    if (cols(v)) pl = _cdid_ref_is_plain(B.reference.roles[role].formula, B.reference.raw_names[v])
    else pl = (strtrim(B.reference.roles[role].formula) == "")
    if (pl) {
        icpt = B.reference.roles[role].intercept
        raw = (cols(v) ? _cdid_xp_raw(B.reference, P, d, (1::P.N), v) : J(P.N, 0, .))
        return(icpt ? (J(P.N, 1, 1), raw) : raw)
    }
    if (B.reference.roles[role].dated)
        _cdid_pr_fail("cdid_chained_request_input_error", "crosspair with covariates that vary over time accepts plain numeric covariates; factor-variable notation is not available there.")
    sl = B.reference.roles[role].slot
    raw = J(P.N, cols(B.reference.designs[sl].value), .)
    raw[B.reference.designs[sl].keys, ] = B.reference.designs[sl].value
    icpt = B.reference.designs[sl].intercept
    return(raw)
}

/* Recipes of the cross-pair links with covariates or weights that vary over time: every link
   of target cohort g reads covariates at the reference date b_g, and weights at the first
   calendar date (fixweights(first_period)) or at b_g (fixweights(base_period)). Designs and
   weights of each recipe are stacked side by side over all frame units (ZR: selection
   designs when sel is set). */
void _cdid_xp_recipes(struct cdid_public_bundle scalar B, struct cdid_xp_cat scalar C,
    string scalar wname, real scalar R, real colvector lrec, real matrix W, real matrix XT,
    real matrix XS, real scalar icpt, real matrix rec, real colvector wcol, real scalar sel,
    real matrix ZR, real scalar sicpt, | real matrix exF, real colvector tctx, real colvector cidx,
    real matrix KP)
{
    struct cdid_population scalar P
    real scalar N, L, l, b, d, a, xd, wd, r, wi, restricted, cx
    real colvector hit
    real matrix X
    string scalar key
    transmorphic RK
    P = B.prepared.population
    N = P.N
    L = rows(C.g)
    xd = B.reference.roles[1].dated | (sel ? B.reference.roles[3].dated : 0)
    wd = (B.reference.weight_kind == "dated")
    wi = 0
    if (wd) {
        hit = selectindex(B.reference.raw_names' :== wname)
        if (length(hit) != 1) _cdid_pr_fail("cdid_chained_request_input_error", "The sampling-weight variable of the cross-pair fit could not be located.")
        wi = hit[1]
    }
    RK = asarray_create("string", 1)
    /* A restricted reference eligibility context (third column) is part of the recipe. */
    restricted = 0
    if (args() >= 15) restricted = (rows(tctx) == rows(C.tg) & cols(exF) > 0)
    rec = J(0, (restricted ? 3 : 2), .)
    lrec = J(L, 1, .)
    for (l = 1; l <= L; l++) {
        b = selectindex(P.calendar :== _cdid_pr_reference(P.calendar, C.g[l], B.request.anticipation))
        if (rows(b) != 1) _cdid_pr_fail("cdid_chained_request_input_error", "A target cohort of the cross-pair comparisons has no reference date in the calendar.")
        d = (xd ? b[1] : 0)
        a = (wd ? (B.reference.fix_weights == "first_period" ? 1 : b[1]) : 0)
        key = strofreal(d) + ":" + strofreal(a)
        cx = (restricted ? tctx[C.gi[l]] : 0)
        if (restricted) key = key + ":" + strofreal(cx > 0 ? cidx[cx] : 0)
        if (!asarray_contains(RK, key)) {
            rec = rec \ (restricted ? (d, a, cx) : (d, a))
            asarray(RK, key, rows(rec))
        }
        lrec[l] = asarray(RK, key)
    }
    R = rows(rec)
    W = J(N, R, .)
    XT = J(N, 0, .)
    XS = J(N, 0, .)
    ZR = J(N, 0, .)
    for (r = 1; r <= R; r++) {
        if (rec[r, 2] > 0) W[, r] = _cdid_xp_raw(B.reference, P, rec[r, 2], (1::N), wi)
        else W[, r] = B.prepared.weight.owner_weight
        X = _cdid_xp_role_design(B, P, 1, rec[r, 1], icpt)
        XT = XT, X
        XS = XS, (icpt ? (cols(X) > 1 ? X[, 2..cols(X)] : J(N, 0, .)) : X)
        if (sel) ZR = ZR, _cdid_xp_role_design(B, P, 3, rec[r, 1], sicpt)
    }
    if (any((W :< .) :& (W :<= 0))) _cdid_pr_fail("cdid_chained_request_input_error", "Sampling weights must be positive at the dates read by the cross-pair comparisons.")
    wcol = rec[, 2]
    /* Units each recipe reads: those eligible in its context (all units otherwise). */
    if (restricted) {
        KP = J(N, R, 1)
        for (r = 1; r <= R; r++) if (rec[r, 3] > 0) KP[, r] = (exF[, rec[r, 3]] :== 0)
    }
    else if (args() >= 18) KP = J(0, 0, .)
}

/* Within-unit scales under restricted reference eligibility (one column per sorted target
   cohort of gu): sqrt(n_A mean(w) / w_i) over the repeat observers eligible in the context
   of the cohort's recipe, with shared weights or the weights of that recipe. */
real matrix _cdid_xp_restricted_scales(struct cdid_xp_cat scalar C, real colvector gu,
    real colvector lrec, real matrix W, real colvector aF, real colvector wshared, real scalar shared,
    real matrix KP)
{
    real matrix SC
    real colvector ug, w, ix
    real scalar nA, i, r
    nA = rows(aF)
    ug = uniqrows(gu)
    SC = J(nA, rows(ug), 0)
    for (i = 1; i <= rows(ug); i++) {
        r = lrec[selectindex(C.g :== ug[i])[1]]
        w = (shared ? wshared : W[aF, r])
        ix = selectindex((w :< .) :& (w :> 0) :& (KP[aF, r] :!= 0))
        if (rows(ix)) SC[ix, i] = sqrt(nA * mean(w[ix]) :/ w[ix])
    }
    return(SC)
}

/* Within-unit scales by target cohort under weight recipes (R build_within_scale_by_g):
   sqrt(n_A / (w_i / mean w)) with the mean over the analysis units of the cohorts that the
   links of the same weight anchor read. Columns follow the sorted target cohorts of gu. */
real matrix _cdid_xp_recipe_scales(struct cdid_xp_cat scalar C, real colvector gu, real colvector lrec,
    real colvector wcol, real matrix W, real colvector aF, real colvector cohA, real colvector wstat)
{
    real matrix SC
    real colvector ug, inr, wa, coh, rel
    real scalar nA, i, l, a, j, k, m
    nA = rows(aF)
    ug = uniqrows(gu)
    SC = J(nA, rows(ug), 0)
    for (i = 1; i <= rows(ug); i++) {
        l = selectindex(C.g :== ug[i])[1]
        a = wcol[lrec[l]]
        if (a == 0) {
            SC[, i] = sqrt(sum(wstat) :/ wstat) :* (wstat :> 0)
            continue
        }
        inr = J(nA, 1, 0)
        for (j = 1; j <= rows(C.g); j++) {
            if (wcol[lrec[j]] != a) continue
            coh = C.g[j] \ C.coh[C.ci[j]]
            for (k = 1; k <= 2; k++) {
                if (coh[k] >= .) inr = inr :| (cohA :>= .)
                else inr = inr :| (cohA :== coh[k])
            }
        }
        wa = W[aF, lrec[l]]
        m = sum(select(wa, inr)) / sum(inr)
        rel = (inr :* (wa :< .)) :* wa :/ m
        rel = editmissing(rel, 0)
        SC[, i] = editmissing(sqrt(nA :/ rel), 0) :* (rel :> 0)
    }
    return(SC)
}

/* Names of the indicator columns of factor-variable terms: the generated columns from take the
   names of their terms to in the prepared covariate and selection designs. */
void cdid_xp_fv_names(string scalar bundle, string scalar from, string scalar to)
{
    pointer(struct cdid_public_bundle scalar) scalar p
    string rowvector f, t, xn, sn
    real scalar j, k
    p = findexternal(bundle)
    if (p == NULL) return
    f = tokens(from)
    t = tokens(to)
    if (cols(f) != cols(t)) return
    xn = p->prepared.x_names
    sn = p->prepared.sel_names
    for (j = 1; j <= cols(f); j++) {
        if (f[j] == t[j]) continue
        for (k = 1; k <= cols(xn); k++) if (xn[k] == f[j]) xn[k] = t[j]
        for (k = 1; k <= cols(sn); k++) if (sn[k] == f[j]) sn[k] = t[j]
    }
    p->prepared.x_names = xn
    p->prepared.sel_names = sn
}

/* 1 if any active variable takes two different nonmissing values within a unit. */
real scalar cdid_xp_varies(string scalar gid, string scalar vars, string scalar touse)
{
    string rowvector v
    real matrix D
    real colvector o, same
    real scalar j
    v = cdid_ref_formula_names(vars)
    for (j = 1; j <= cols(v); j++) {
        if (_st_varindex(v[j]) >= .) continue
        if (st_isstrvar(v[j])) continue
        D = st_data(., (gid, v[j]), touse)
        D = select(D, D[, 2] :< .)
        if (rows(D) < 2) continue
        o = order(D, (1, 2))
        D = D[o, ]
        same = (D[2..rows(D), 1] :== D[1..rows(D) - 1, 1])
        if (any(same :& (D[2..rows(D), 2] :!= D[1..rows(D) - 1, 2]))) return(1)
    }
    return(0)
}

struct cdid_public_fit_result scalar cdid_crosspair_public_fit(struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_public_report_buffer scalar buf
    struct cdid_xp_cat scalar C
    struct cdid_xp_fit scalar F
    struct cdid_xps_sel scalar S
    struct cdid_xps_moments scalar XM
    struct cdid_gmm_unit_gram scalar src
    pointer() scalar q
    real matrix Y, O, grid, G, postIF, auxIF, Z, Pc
    real colvector aF, w, coh, targets, used, rr, k, post_point, post_av, aux_point, aux_sup, post_ids, aux_ids, cix
    real colvector acell, asign
    string colvector post_status, post_reason, aux_reason
    string rowvector req
    real scalar N, T, nA, t, j, g, h, qd, sg, c, pt_start, Jt, Q, sel, md, p, gen, dated, Rr, icr
    real colvector lrec, wcol, wst
    real matrix Wr, XTr, XSr, recs, ZRr, exF, KP
    real colvector tctx, cidx
    real scalar sicr, ps, mref
    struct cdid_xps_sel colvector SS
    string scalar pt, wt, smeth
    sel = (B.prepared.selection_canonical != "none" & B.prepared.selection_canonical != "")
    smeth = (sel ? B.prepared.selection_est_method : "")
    if (sel & smeth != "ipw" & smeth != "imputation" & smeth != "dr")
        _cdid_pr_fail("cdid_chained_request_input_error", "crosspair with selection(conditional_mtar) requires selmethod(ipw), selmethod(imputation) or selmethod(dr).")
    dated = (B.prepared.design == "reference_recipes")
    cdid_xpn_clear()
    if (cdid_nuis_requested()) cdid_xpn_bundle(B)
    q = findexternal("__cdid_xp_req")
    req = *q
    pt = req[1]
    pt_start = strtoreal(req[2])
    if (pt == "ctrl" & pt_start >= .) pt_start = min(B.prepared.population.calendar)
    N = B.prepared.population.N
    T = rows(B.prepared.population.calendar)
    aF = B.prepared.population.A_to_F
    nA = rows(aF)
    if (nA < 2) _cdid_pr_fail("cdid_chained_request_input_error", "crosspair needs at least two units with two observed outcomes.")
    if (rows(B.prepared.population.obs_grid) == N & cols(B.prepared.population.obs_grid) == T)
        grid = B.prepared.population.obs_grid[aF, .]
    else {
        grid = J(nA, T, 0)
        for (t = 1; t <= T; t++) grid[, t] = cdid_population_rows(B.prepared.population, aF, J(nA, 1, B.prepared.population.calendar[t]))
    }
    Y = J(nA, T, .)
    O = J(nA, T, 0)
    for (t = 1; t <= T; t++) {
        k = selectindex(grid[, t] :> 0)
        if (rows(k) == 0) continue
        Y[k, t] = B.prepared.population.y[grid[k, t]]
        O[k, t] = (Y[k, t] :< .)
    }
    w = B.prepared.weight.owner_weight[aF]
    coh = B.prepared.population.cohort[aF]
    targets = uniqrows(select(coh, coh :< .))
    targets = select(targets, (max(B.prepared.population.calendar) + B.request.anticipation) :>= targets)
    /* Restricted reference eligibility contexts (missingreference(exclude)). */
    mref = 0
    tctx = cidx = J(0, 1, .)
    exF = KP = J(0, 0, .)
    if (dated) if (findexternal("__cdid_mref_state") != NULL) if (cdid_mrxp_active()) {
        exF = cdid_mrxp_excluded(N)
        if (cols(exF)) {
            mref = 1
            cidx = cdid_mrxp_index()
            tctx = cdid_mrxp_ctx_of(sort(uniqrows(targets), 1))
            cdid_mrxp_aux_clear()
        }
    }
    if (mref) C = cdid_xp_catalog(O, w, coh, B.prepared.population.calendar, targets, B.request.anticipation,
        B.request.control_group, exF[aF, ], tctx)
    else C = cdid_xp_catalog(O, w, coh, B.prepared.population.calendar, targets, B.request.anticipation, B.request.control_group)
    wt = (B.settings.gmm_weighting == "two_step" ? "optimal" : B.settings.gmm_weighting)
    if (rows(C.g) == 0) _cdid_pr_fail("cdid_chained_request_input_error", "The cross-pair comparisons are empty for the requested cohorts and control group.")
    gen = (sel | dated | ((B.prepared.method == "reg" | B.prepared.method == "dr") & B.prepared.outcome_pooling == "pooled"))
    if (mref) _cdid_xp_recipes(B, C, (cols(req) >= 3 ? req[3] : ""), Rr, lrec, Wr, XTr, XSr, icr, recs, wcol,
        sel, ZRr, sicr, exF, tctx, cidx, KP)
    else if (dated) _cdid_xp_recipes(B, C, (cols(req) >= 3 ? req[3] : ""), Rr, lrec, Wr, XTr, XSr, icr, recs, wcol,
        sel, ZRr, sicr)
    if (gen) {
        /* Selection-model weights or pooled outcome slopes: moments and frame-unit influence
           functions on the full catalogue, then the verified forest reduction. */
        if (dated & sel & smeth != "imputation") {
            /* Selection models of each recipe: its own design, weights and rates. */
            SS = cdid_xps_sel(Rr, 1)
            ps = cols(ZRr) / Rr
            for (j = 1; j <= Rr; j++) {
                if (mref) SS[j] = cdid_xps_sel_fit(B, C, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, ctrl.chunk,
                    selectindex(lrec :== j), (ps ? ZRr[, ((j - 1) * ps + 1)..(j * ps)] : J(N, 0, .)), sicr, Wr[, j], KP[, j])
                else SS[j] = cdid_xps_sel_fit(B, C, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, ctrl.chunk,
                    selectindex(lrec :== j), (ps ? ZRr[, ((j - 1) * ps + 1)..(j * ps)] : J(N, 0, .)), sicr, Wr[, j])
            }
        }
        else if (dated) SS = S
        else if (sel & smeth != "imputation") S = cdid_xps_sel_fit(B, C, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, ctrl.chunk)
        p = cols(B.prepared.X)
        if (sel & smeth != "ipw" & dated) {
            /* Imputed outcome changes on the recipe of each comparison. */
            XM = cdid_xprd_eval(B, C, SS, (smeth == "dr"), B.prepared.method,
                (B.prepared.outcome_pooling == "" ? "separate" : B.prepared.outcome_pooling),
                (B.prepared.imputation_pooling == "" ? "separate" : B.prepared.imputation_pooling),
                ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, (wt == "within"), ctrl.baseline.rank_tol,
                ctrl.chunk, ctrl.baseline.max_cached, ctrl.workspace_bytes, ctrl.workcap, Rr, lrec, Wr, XTr, XSr, icr, KP)
        }
        else if (sel & smeth != "ipw") {
            /* Imputed outcome changes, augmented by the observation model with selmethod(dr). */
            XM = cdid_xpr_eval(B, C, S, (smeth == "dr"), Y, B.prepared.method,
                (B.prepared.outcome_pooling == "" ? "separate" : B.prepared.outcome_pooling),
                (B.prepared.imputation_pooling == "" ? "separate" : B.prepared.imputation_pooling),
                ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, (wt == "within"), ctrl.baseline.rank_tol,
                ctrl.chunk, ctrl.baseline.max_cached, ctrl.workspace_bytes, ctrl.workcap)
        }
        else if (dated) XM = cdid_xps_eval(SS, sel, C, Y, aF, XTr, XSr, B.prepared.method,
            (B.prepared.outcome_pooling == "" ? "separate" : B.prepared.outcome_pooling),
            ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, (wt == "within"),
            B.prepared.weight.owner_weight, B.prepared.population.cohort, B.prepared.population.A,
            icr, Rr, lrec, Wr, KP)
        else XM = cdid_xps_eval(S, sel, C, Y, aF, (p ? B.prepared.X : J(N, 1, 1)),
            (B.prepared.x_has_intercept ? (p > 1 ? B.prepared.X[, 2..p] : J(N, 0, .)) : B.prepared.X),
            B.prepared.method, (B.prepared.outcome_pooling == "" ? "separate" : B.prepared.outcome_pooling),
            ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, (wt == "within"),
            B.prepared.weight.owner_weight, B.prepared.population.cohort, B.prepared.population.A,
            B.prepared.x_has_intercept)
        Pc = XM.Phi :- (colsum(XM.Phi) :/ N)
        if (wt == "identity") {
            used = (1::rows(C.g))
            G = I(rows(used))
        }
        else cdid_xps_forest(C, XM.mhat, (wt == "optimal" ? Pc : J(0, 0, .)), XM.U, nA, used, G, md)
        if (cdid_nuis_requested() & dated) cdid_xpn_dated(B, (sel & smeth != "imputation"))
        else if (cdid_nuis_requested()) cdid_xpn_finish(B, used, ctrl.inference.compute_inffunc, wt)
        F.mhat = XM.mhat[used]
        F.Phi = Pc[, used] :* (nA / N)
        Pc = J(0, 0, .)
        if (wt == "identity") F.V = J(rows(used), rows(used), 0)
        else if (wt == "within" & mref) {
            wst = (B.reference.weight_kind == "dated" & B.reference.fix_weights == "first_period" ?
                B.prepared.weight.V[aF] : B.prepared.weight.owner_weight[aF])
            F.V = cdid_xps_within_cov(XM.U[, used], nA, C.g[used], w,
                _cdid_xp_restricted_scales(C, C.g[used], lrec, Wr, aF, wst,
                    (B.reference.weight_kind != "dated" | (B.reference.fix_weights == "first_period" & !missing(wst) & all(wst :> 0))), KP))
        }
        else if (wt == "within" & dated) {
            /* Static weights keep their within-unit scale; first-period weights complete on the
               analysis units use it too; otherwise the scale follows each weight anchor. */
            wst = (B.reference.weight_kind == "dated" & B.reference.fix_weights == "first_period" ?
                B.prepared.weight.V[aF] : B.prepared.weight.owner_weight[aF])
            if (B.reference.weight_kind == "dated" & B.reference.fix_weights == "first_period" & !missing(wst) & all(wst :> 0))
                F.V = cdid_xps_within_cov(XM.U[, used], nA, C.g[used], wst)
            else F.V = cdid_xps_within_cov(XM.U[, used], nA, C.g[used], w,
                _cdid_xp_recipe_scales(C, C.g[used], lrec, (B.reference.weight_kind == "dated" ? wcol : J(rows(wcol), 1, 0)),
                    Wr, aF, coh, wst))
        }
        else if (wt == "within") F.V = cdid_xps_within_cov(XM.U[, used], nA, C.g[used], w)
        else {
            F.V = quadcross(F.Phi, F.Phi) :/ nA
            F.V = 0.5 :* (F.V + F.V')
        }
        src.n_units = N
        src.rows_consumed = N
        src.normalization_count = nA
        src.normalization_factor = 1 / nA
        src.sample_rank_upper = N - 1
        src.mean = colsum(F.Phi) :/ N
        src.retained_if = J(0, 0, .)
        src.covariance = F.V
        cdid_xp_targets(F, C, used, G, rows(C.g), B.request.anticipation, pt, pt_start,
            B.request.base_period, wt, ctrl.condition_trigger, ctrl.ridge_relative, src)
    }
    else {
        if (wt == "identity") {
            used = (1::rows(C.g))
            G = I(rows(used))
        }
        else {
            cdid_xp_forest(C)
            used = C.kept
            G = C.G
        }
        Z = (cols(B.prepared.X) ? B.prepared.X[aF, .] : J(nA, 1, 1))
        if (cols(Z) == 1 & B.prepared.x_has_intercept) Z = J(nA, 1, 1)
        if (cdid_nuis_requested()) cdid_xpn_base_init(B.prepared.method, (cols(Z) == 1 & all(Z :== Z[1, 1])))
        F = cdid_xp_solve(C, Y, w, used, G, rows(C.g), B.request.anticipation, pt, pt_start,
            B.request.base_period, wt, ctrl.condition_trigger, ctrl.ridge_relative, Z,
            B.prepared.method, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit)
        if (cdid_nuis_requested()) cdid_xpn_finish(B, used, ctrl.inference.compute_inffunc, wt)
    }
    cdid_chain_pending(cdid_chain_xp(C, F, used, B.request.anticipation, pt, pt_start, B.request.base_period,
        wt, B.request.report, B.request.post))
    /* Ledger quantities. */
    Jt = rows(B.request.post)
    Q = rows(B.request.aux)
    post_point = J(Jt, 1, .)
    post_av = J(Jt, 1, 0)
    post_status = J(Jt, 1, "disconnected_component")
    post_reason = J(Jt, 1, "the cross-pair comparisons do not connect this cell to the reference date")
    cix = J(Jt, 1, 0)
    for (h = 1; h <= Jt; h++) {
        c = _cdid_xp_cell(F, B.request.post[h, 1], B.request.post[h, 2])
        cix[h] = c
        if (c & F.avail[c]) {
            post_point[h] = F.att[c]
            post_av[h] = 1
            post_status[h] = "estimable"
            post_reason[h] = ""
        }
    }
    aux_point = J(Q, 1, .)
    aux_sup = J(Q, 1, 0)
    aux_reason = J(Q, 1, "the cross-pair comparisons do not connect this pre-treatment cell to the reference date")
    acell = J(Q, 1, 0)
    asign = J(Q, 1, 0)
    for (j = 1; j <= rows(B.request.report); j++) {
        if (B.request.report[j, 5] != 2) continue
        qd = B.request.report[j, 7]
        sg = B.request.report[j, 8]
        if (qd < 1 | qd > Q | qd >= .) continue
        c = _cdid_xp_cell(F, B.request.report[j, 1], B.request.report[j, 2])
        if (!c) continue
        if (!F.avail[c]) continue
        if (aux_sup[qd] & acell[qd] != c) _cdid_pr_fail("cdid_native_reporting_input", "Two cross-pair cells share one pretrend descriptor.")
        aux_point[qd] = sg * F.att[c]
        aux_sup[qd] = 1
        aux_reason[qd] = "supported"
        acell[qd] = c
        asign[qd] = sg
    }
    buf = cdid_public_report_start(B, post_point, post_av, post_status, post_reason, aux_point, aux_sup, aux_reason,
        ctrl.inference.compute_inffunc, ctrl.inference.limits)
    if (ctrl.inference.compute_inffunc) {
        post_ids = selectindex(post_av)
        aux_ids = selectindex(aux_sup)
        postIF = J(N, rows(post_ids), 0)
        auxIF = J(N, rows(aux_ids), 0)
        if (gen) {
            for (h = 1; h <= rows(post_ids); h++) postIF[, h] = (N / nA) :* F.IF[, cix[post_ids[h]]]
            for (h = 1; h <= rows(aux_ids); h++) auxIF[, h] = (asign[aux_ids[h]] * N / nA) :* F.IF[, acell[aux_ids[h]]]
        }
        else {
            for (h = 1; h <= rows(post_ids); h++) postIF[aF, h] = (N / nA) :* F.IF[, cix[post_ids[h]]]
            for (h = 1; h <= rows(aux_ids); h++) auxIF[aF, h] = (asign[aux_ids[h]] * N / nA) :* F.IF[, acell[aux_ids[h]]]
        }
        cdid_public_report_add(buf, (1::N), post_ids, postIF, aux_ids, auxIF)
    }
    out.payload = cdid_public_report_finish(buf, ctrl.inference.limits)
    out.payload.numeric_settings = out.payload.numeric_settings \ rows(C.g) \ rows(used) \ F.diag'
    out.payload.numeric_setting_names = out.payload.numeric_setting_names \ "cross_pair_links" \ "cross_pair_links_used" \
        ("GMM_lambda" \ "GMM_ridge_rho" \ "GMM_covariance_rank" \ "GMM_d_min" \ "GMM_d_max" \ "GMM_criterion")
    out.payload.text_settings = out.payload.text_settings \ ("cross_pair", "1") \ ("pt_assumption", pt) \
        ("GMM_trigger_reason", F.trigger_reason)
    out.inference = ctrl.inference
    cdid_state_validate(out.payload, ctrl.inference.limits)
    return(out)
}

/* Restricted reference eligibility contexts; called only while such a context is active. */
real scalar cdid_mrxp_active()
{
    return(cdid_mref_active())
}
real matrix cdid_mrxp_excluded(real scalar N)
{
    return(cdid_mref_excluded(N))
}
real colvector cdid_mrxp_index()
{
    return(cdid_mref_index())
}
real colvector cdid_mrxp_ctx_of(real colvector cohort_values)
{
    return(cdid_mref_ctx_of(cohort_values))
}
void cdid_mrxp_aux_clear()
{
    cdid_mref_aux_clear()
}

end
