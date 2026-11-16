*! Opt-in reports of the fitted nuisance models (keepnuisance); no refit, no default cost.
version 19.5
mata:

/* One captured fit. comp rows: selection, outcome, treatment, imputation;
   columns: name, status, code, reason.
   blocks columns: 1 n_cells 2 n_parameters 3 n_intercepts 4 n_slopes 5 n_fitting_rows
     6 n_units 7 n_fitting_clusters 8 unit_weight_sum 9 converged 10 iterations 11 rank
     12 n_pair_cohorts 13 n_comparisons 14 cohort (. if several)
   coef columns: 1 block 2 cell 3 estimate 4 std_error 5 z 6 p
   cells columns: 1 cell 2 block 3 type (1 pair, 2 date) 4 cohort 5 s 6 t 7 n_rows
     8 n_observed 9 observed_share_weighted 10 p_min 11 p_max 12 n_floor */
/* One role of the augmented DR+DR route (outcome, treatment, imputation).
   blocks columns: 1 n_parameters 2 n_intercepts 3 n_slopes 4 n_fitting_rows
     5 weight_sum 6 n_contributing_clusters 7 n_comparisons 8 backend_id
   block_txt columns: 1 id 2 label 3 equation 4 cohorts 5 kept slopes 6 aliased slopes
   coef columns as in cdid_nuis_cap; coef_txt columns: term, role. */
struct cdid_nuis_role {
    string scalar role, title, response, sample, weights, pooling
    real matrix blocks, coef
    string matrix block_txt, coef_txt
}

struct cdid_nuis_cap {
    string scalar token, mode, method, selmethod, selection, clustervar
    real scalar n_catalog, n_aux
    string matrix comp
    string scalar sel_model, sel_pooling, sel_universe
    string rowvector slope_names
    real scalar N, clustered, n_groups, has_intercept, floor, maxit
    real matrix blocks, coef, cells
    string matrix block_txt, coef_txt
    string colvector cell_key
    string matrix comparisons
    struct cdid_nuis_role rowvector roles
}

void cdid_nuis_begin(real scalar on)
{
    pointer() scalar p
    p = findexternal("__cdid_nuis_req")
    if (p == NULL) p = crexternal("__cdid_nuis_req")
    *p = (on != 0)
    if (findexternal("__cdid_nuis_pending") != NULL) rmexternal("__cdid_nuis_pending")
}

real scalar cdid_nuis_requested()
{
    pointer() scalar p
    p = findexternal("__cdid_nuis_req")
    if (p == NULL) return(0)
    return(*p == 1)
}

/* Captures by result token (the 64 most recent), so that a restored or
   loaded result finds its own reports. */
void _cdid_nuis_reg_put(string scalar token, struct cdid_nuis_cap scalar cap)
{
    pointer() scalar p, o, q
    string colvector ord
    p = findexternal("__cdid_nuis_registry")
    if (p == NULL) {
        p = crexternal("__cdid_nuis_registry")
        *p = asarray_create()
    }
    o = findexternal("__cdid_nuis_order")
    if (o == NULL) {
        o = crexternal("__cdid_nuis_order")
        *o = J(0, 1, "")
    }
    asarray(*p, token, cap)
    ord = *o
    if (rows(ord)) ord = select(ord, ord :!= token)
    ord = ord \ token
    while (rows(ord) > 64) {
        asarray_remove(*p, ord[1])
        ord = ord[2..rows(ord)]
    }
    *o = ord
    q = findexternal("__cdid_nuis_cap")
    if (q == NULL) q = crexternal("__cdid_nuis_cap")
    *q = cap
}

real scalar _cdid_nuis_reg_get(string scalar token, struct cdid_nuis_cap scalar cap)
{
    pointer() scalar p
    p = findexternal("__cdid_nuis_registry")
    if (p == NULL | token == "") return(0)
    if (!asarray_contains(*p, token)) return(0)
    cap = asarray(*p, token)
    return(1)
}

real scalar cdid_nuis_has(string scalar token)
{
    struct cdid_nuis_cap scalar cap
    return(_cdid_nuis_reg_get(token, cap))
}

/* Saved results: the reports follow the result state in the same file. */
void cdid_nuis_append(string scalar token, string scalar path)
{
    struct cdid_nuis_cap scalar cap
    real scalar fh, rc
    if (!_cdid_nuis_reg_get(token, cap)) return
    fh = _fopen(path, "a")
    if (fh < 0) _cdid_rs_fail("cannot add the nuisance reports to the saved file", -fh)
    _cdid_rs_write(fh, ("cdid-nuisance-reports", "1"))
    _cdid_rs_write(fh, cap)
    rc = _fclose(fh)
    if (rc < 0) _cdid_rs_fail("Native result-state close failed.", -rc)
}

void cdid_nuis_loaded(string scalar token, transmorphic matrix raw)
{
    struct cdid_nuis_cap scalar cap
    if (eltype(raw) != "struct") _cdid_rs_fail("Invalid saved nuisance reports.", 498)
    if (structname(raw) != "cdid_nuis_cap") _cdid_rs_fail("Invalid saved nuisance reports.", 498)
    cap = raw
    cap.token = token
    _cdid_nuis_reg_put(token, cap)
}

string scalar _cdid_nuis_num(real scalar x)
{
    if (x >= .) return("0")
    return(strofreal(x, "%21.0g"))
}

string scalar _cdid_nuis_g(real scalar g)
{
    return(g >= . ? "0" : strofreal(g, "%21.0g"))
}

/* Capture the observation logits of the static selected route (selmethod(ipw)
   with method(ipw|reg|dr)). Coefficient SEs use the estimator's own unit
   influence functions: Var(theta_j) = N^-2 sum_g (S_gj - Sbar_j)^2, with
   S_gj the sum of the unit influence functions of cluster g (or of unit g). */
struct cdid_nuis_cap scalar _cdid_nuis_sel_capture(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, real matrix pc_lookup, real colvector calendar,
    struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan,
    real scalar allow_legacy)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_sel_block_fit scalar blk
    struct cdid_sel_if_chunk scalar ch
    pointer() scalar p
    real scalar N, nb, b, c, j, k, cellj, d_s, has_int, b_pos, ncoef, row, lo, hi, chunk
    real scalar G, nl, l, pc, s_cal, t_cal, kk, nT, ncell, clus
    real matrix S, sum1, sum2, Zc, pr, pos_alpha, pos_beta
    real colvector var, se, ids, units, yc, wc, cl, hit, blks, coh_list
    string scalar lab, key, eqn, sl
    real colvector seen, ord, newid
    real scalar nseen, first_seen_on, legacy

    cap.mode = "chained"
    cap.comp = J(4, 4, "")
    N = fit.N
    nb = fit.n_blocks
    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept
    cap.N = N
    cap.sel_model = fit.selection_model
    cap.sel_pooling = fit.selection_pooling
    cap.sel_universe = fit.selection_universe
    cap.has_intercept = has_int
    cap.floor = 1e-6
    cap.maxit = fit.maxit
    cap.clustered = B.prepared.clustered
    cap.n_groups = (B.prepared.clustered ? B.prepared.C : N)
    cap.slope_names = J(1, d_s, "")
    for (j = 1; j <= d_s; j++) {
        k = fit.spec.slope_cols[j]
        if (k <= cols(B.prepared.sel_names)) cap.slope_names[j] = B.prepared.sel_names[k]
        else cap.slope_names[j] = "z" + strofreal(j)
    }

    /* Parameter positions in theta = (alpha over cells, beta stacked by block). */
    pos_alpha = J(fit.n_cells, 1, .)
    if (has_int) pos_alpha = (1::fit.n_cells)
    pos_beta = J(nb, 1, .)
    b_pos = fit.n_alpha + 1
    for (b = 1; b <= nb; b++) {
        pos_beta[b] = b_pos
        b_pos = b_pos + d_s
    }

    /* Influence-function moments, streamed over units. */
    cl = (B.prepared.clustered ? B.prepared.unit_cluster : (1::N))
    G = cap.n_groups
    S = J(G, fit.n_alpha + fit.n_beta, 0)
    chunk = max((1, min((N, floor(4e6 / max((1, fit.n_alpha + fit.n_beta)))))))
    for (lo = 1; lo <= N; lo = lo + chunk) {
        hi = min((N, lo + chunk - 1))
        ids = (lo::hi)
        ch = cdid_sel_eval_if_chunk(fit, rdr, ids, 0)
        if (B.prepared.clustered) {
            for (row = 1; row <= rows(ids); row++) {
                S[cl[ids[row]], ] = S[cl[ids[row]], ] + ch.theta_if[row, ]
            }
        }
        else S[ids, ] = ch.theta_if
    }
    sum1 = colsum(S)
    sum2 = colsum(S :* S)
    var = ((sum2 - (sum1 :* sum1) :/ G) :/ (N * N))'
    se = sqrt(var :* (var :> 0))

    /* Display order: the default separate pair logits of method(ipw) are listed
       in order of first use by the comparisons; those of method(reg) by pair and
       cohort; every other specification in the order of the fitted blocks. */
    legacy = allow_legacy & (fit.selection_model == "pair_logit" & fit.selection_pooling == "separate" & fit.spec.d > 0 &
        (B.prepared.method == "ipw" | B.prepared.method == "reg"))
    first_seen_on = (legacy & B.prepared.method == "ipw")
    seen = J(nb, 1, 0)
    nseen = 0
    /* Comparisons: every supported catalog descriptor, then every supported
       auxiliary descriptor, one entry each (an auxiliary comparison that repeats
       a catalog comparison is listed and counted again). */
    nT = rows(calendar)
    nl = rows(plan.raw)
    cap.comparisons = J(nl, 5, "")
    blks = J(nb, 1, 0)
    for (l = 1; l <= nl; l++) {
        if (!plan.raw_supported[l]) {
            cap.comparisons[l, 1] = "unsupported"
            continue
        }
        s_cal = plan.raw[l, 2]
        t_cal = plan.raw[l, 3]
        coh_list = plan.raw[l, 1] \ plan.raw_K_index[|plan.raw_K_ptr[l] \ plan.raw_K_ptr[l + 1] - 1|]
        cap.comparisons[l, ] = ((l <= plan.M ? "catalog" : "auxiliary"),
            _cdid_nuis_g(plan.cohorts[coh_list[1]]),
            _cdid_nuis_num(calendar[s_cal]), _cdid_nuis_num(calendar[t_cal]),
            _cdid_nuis_glist(plan.cohorts[coh_list[|2 \ rows(coh_list)|]], rows(coh_list) - 1))
        hit = J(nb, 1, 0)
        for (kk = 1; kk <= rows(coh_list); kk++) {
            pc = pc_lookup[(s_cal - 1) * nT + t_cal, coh_list[kk]]
            if (pc < 1) continue
            hit[fit.map.cells[fit.map.pair_cohorts[pc, 6], 7]] = 1
            hit[fit.map.cells[fit.map.pair_cohorts[pc, 7], 7]] = 1
        }
        blks = blks + hit
        if (first_seen_on) _cdid_nuis_seen(seen, nseen, coh_list, s_cal, t_cal, pc_lookup, nT, fit)
    }
    cap.n_catalog = sum(cap.comparisons[, 1] :== "catalog")

    /* Cells. */
    ncell = fit.n_cells
    cap.cells = J(ncell, 12, .)
    cap.cell_key = J(ncell, 1, "")
    for (c = 1; c <= ncell; c++) {
        k = fit.map.cells[c, 3]
        units = _cdid_sel_index((rdr.cohort_pos :== k) :& (rdr.U :== 1))
        cap.cells[c, 1] = c
        cap.cells[c, 2] = fit.cell_block[c]
        cap.cells[c, 3] = fit.map.cells[c, 2]
        cap.cells[c, 4] = fit.cohorts[k]
        if (fit.map.cells[c, 2] == 1) {
            cap.cells[c, 5] = fit.all_times[fit.map.cells[c, 5]]
            cap.cells[c, 6] = fit.all_times[fit.map.cells[c, 6]]
            cap.cell_key[c] = "g=" + _cdid_nuis_g(cap.cells[c, 4]) + ",(" + _cdid_nuis_num(cap.cells[c, 5]) + "," + _cdid_nuis_num(cap.cells[c, 6]) + ")"
        }
        else {
            cap.cells[c, 6] = fit.all_times[fit.map.cells[c, 4]]
            cap.cell_key[c] = "g=" + _cdid_nuis_g(cap.cells[c, 4]) + ",t=" + _cdid_nuis_num(cap.cells[c, 6])
        }
        cap.cells[c, 7] = rows(units)
        if (rows(units)) {
            yc = _cdid_sel_cell_y(rdr, fit.map, c, units)
            wc = rdr.w[units]
            pr = cdid_sel_pred_cell(fit, c, rdr.Z[units, ])
            cap.cells[c, 8] = sum(yc)
            cap.cells[c, 9] = sum(wc :* yc) / sum(wc)
            cap.cells[c, 10] = min(pr[, 2])
            cap.cells[c, 11] = max(pr[, 2])
            cap.cells[c, 12] = sum(pr[, 2] :< cap.floor)
        }
    }

    /* Blocks and coefficients. */
    cap.blocks = J(nb, 14, .)
    cap.block_txt = J(nb, 3, "")
    ncoef = 0
    for (b = 1; b <= nb; b++) ncoef = ncoef + (has_int ? fit.blocks[b].B_m : 0) + d_s
    cap.coef = J(ncoef, 6, .)
    cap.coef_txt = J(ncoef, 2, "")
    row = 0
    sl = ""
    for (j = 1; j <= d_s; j++) sl = sl + " + b" + strofreal(j) + "*" + cap.slope_names[j]
    for (b = 1; b <= nb; b++) {
        blk = fit.blocks[b]
        units = J(0, 1, .)
        for (j = 1; j <= blk.B_m; j++) {
            cellj = blk.cell_ids[j]
            units = units \ _cdid_sel_index((rdr.cohort_pos :== fit.map.cells[cellj, 3]) :& (rdr.U :== 1))
        }
        units = uniqrows(units)
        cap.blocks[b, 1] = blk.B_m
        cap.blocks[b, 3] = (has_int ? blk.B_m : 0)
        cap.blocks[b, 4] = d_s
        cap.blocks[b, 2] = cap.blocks[b, 3] + d_s
        cap.blocks[b, 5] = sum(cap.cells[blk.cell_ids, 7])
        cap.blocks[b, 6] = rows(units)
        if (B.prepared.clustered & rows(units)) cap.blocks[b, 7] = rows(uniqrows(cl[units]))
        cap.blocks[b, 8] = (rows(units) ? sum(rdr.w[units]) : 0)
        cap.blocks[b, 9] = 1
        cap.blocks[b, 10] = blk.iterations
        cap.blocks[b, 11] = blk.rank
        cap.blocks[b, 12] = sum(select(J(rows(fit.map.pair_cohorts), 1, 1),
            (fit.cell_block[fit.map.pair_cohorts[, 6]] :== b) :| (fit.cell_block[fit.map.pair_cohorts[, 7]] :== b)))
        cap.blocks[b, 13] = blks[b]
        cap.blocks[b, 14] = (rows(uniqrows(cap.cells[blk.cell_ids, 4])) == 1 ? cap.cells[blk.cell_ids[1], 4] : .)
        cap.block_txt[b, 1] = "S" + strofreal(b)
        key = cap.cell_key[blk.cell_ids[1]]
        if (cap.sel_pooling == "separate") lab = (cap.cells[blk.cell_ids[1], 3] == 1 ? "pair[" : "date[") + key + "]"
        else if (cap.sel_pooling == "within_cohort") lab = "within_cohort[g=" + _cdid_nuis_g(cap.blocks[b, 14]) + "]"
        else {
            coh_list = uniqrows(cap.cells[blk.cell_ids, 4])
            coh_list = editmissing(coh_list, 0)
            coh_list = sort(coh_list, 1)
            lab = "across_cohorts[g in {"
            for (j = 1; j <= rows(coh_list); j++) lab = lab + (j > 1 ? "," : "") + _cdid_nuis_num(coh_list[j])
            lab = lab + "}]"
        }
        cap.block_txt[b, 2] = lab
        if (cap.sel_model == "pair_logit") eqn = "logit Pr(S_i(s,t) = 1 | X_i, G_i, U_i = 1) = "
        else eqn = "logit Pr(R_it = 1 | X_i, G_i, U_i = 1) = "
        eqn = eqn + (has_int ? (blk.B_m > 1 ? "a[cell]" : "a") : "0") + sl
        cap.block_txt[b, 3] = eqn
        for (j = 1; j <= (has_int ? blk.B_m : 0); j++) {
            row++
            cellj = blk.cell_ids[j]
            cap.coef[row, 1] = b
            cap.coef[row, 2] = cellj
            cap.coef[row, 3] = fit.alpha[cellj]
            cap.coef[row, 4] = se[pos_alpha[cellj]]
            cap.coef_txt[row, ] = ("(Intercept)[" + cap.cell_key[cellj] + "]", "local_intercept")
        }
        for (j = 1; j <= d_s; j++) {
            row++
            cap.coef[row, 1] = b
            cap.coef[row, 3] = fit.beta[pos_beta[b] - fit.n_alpha + j - 1]
            cap.coef[row, 4] = se[pos_beta[b] + j - 1]
            cap.coef_txt[row, ] = (cap.slope_names[j], (blk.B_m > 1 ? "shared_slope" : "slope"))
        }
        if (B.prepared.clustered & cap.blocks[b, 7] < 2) {
            cap.coef[select((1::rows(cap.coef)), cap.coef[, 1] :== b), 4] = J(sum(cap.coef[, 1] :== b), 1, .)
        }
    }
    cap.coef[, 5] = cap.coef[, 3] :/ cap.coef[, 4]
    cap.coef[, 6] = 2 :* normal(-abs(cap.coef[, 5]))
    if (legacy) {
        if (first_seen_on) {
            ord = select(seen, seen :> 0)
            ord = ord \ select((1::nb), !_cdid_nuis_in((1::nb), ord))
        }
        else {
            ord = J(nb, 3, .)
            for (b = 1; b <= nb; b++) {
                c = fit.blocks[b].cell_ids[1]
                ord[b, ] = (fit.map.cells[c, 5], fit.map.cells[c, 6], fit.map.cells[c, 3])
            }
            ord = order(ord, (1, 2, 3))
        }
        _cdid_nuis_reorder(cap, ord)
    }
    cap.n_aux = sum(cap.comparisons[, 1] :== "auxiliary")
    cap.comparisons = select(cap.comparisons, cap.comparisons[, 1] :!= "unsupported")
    cap.comp[1, ] = ("selection", "available", "", "")
    return(cap)
}

void cdid_nuis_capture_selected(struct cdid_selected_model scalar m,
    struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan)
{
    struct cdid_nuis_cap scalar cap
    if (rows(m.contexts) != 1) {
        cap.mode = "chained"
        cap.comp = J(4, 4, "")
        cap.comp[1, ] = ("selection", "unavailable", "reference_recipe_route",
            "With covariates or weights that vary over time, the observation logits are fitted separately for each reference period and are not reported.")
        _cdid_nuis_store_pending(cap)
        return
    }
    cap = _cdid_nuis_sel_capture(m.contexts[1].sel_fit, m.reader, m.contexts[1].pc_lookup,
        m.calendar_values, B, plan, 1)
    _cdid_nuis_store_pending(cap)
}

string scalar _cdid_nuis_glist(real colvector g, real scalar n)
{
    string scalar out
    real scalar i
    out = ""
    for (i = 1; i <= n; i++) out = out + (i > 1 ? "," : "") + _cdid_nuis_g(g[i])
    return(out)
}

real colvector _cdid_nuis_in(real colvector x, real colvector set)
{
    real colvector r
    real scalar i
    r = J(rows(x), 1, 0)
    for (i = 1; i <= rows(set); i++) r = r :| (x :== set[i])
    return(r)
}

void _cdid_nuis_seen(real colvector seen, real scalar nseen, real colvector coh_list,
    real scalar s_cal, real scalar t_cal, real matrix pc_lookup, real scalar nT,
    struct cdid_sel_fit_result scalar fit)
{
    real scalar kk, pc, b
    for (kk = 1; kk <= rows(coh_list); kk++) {
        if (coh_list[kk] >= .) continue
        pc = pc_lookup[(s_cal - 1) * nT + t_cal, coh_list[kk]]
        if (pc < 1) continue
        b = fit.map.cells[fit.map.pair_cohorts[pc, 6], 7]
        if (!any(seen[|1 \ max((1, nseen))|] :== b) | nseen == 0) {
            nseen++
            seen[nseen] = b
        }
    }
}

void _cdid_nuis_reorder(struct cdid_nuis_cap scalar cap, real colvector ord)
{
    real colvector newid, rows_new
    real scalar b, nb
    string matrix bt
    nb = rows(ord)
    newid = J(nb, 1, .)
    newid[ord] = (1::nb)
    cap.blocks = cap.blocks[ord, ]
    bt = cap.block_txt[ord, ]
    for (b = 1; b <= nb; b++) bt[b, 1] = "S" + strofreal(b)
    cap.block_txt = bt
    cap.cells[, 2] = newid[cap.cells[, 2]]
    cap.coef[, 1] = newid[cap.coef[, 1]]
    rows_new = order((cap.coef[, 1], (1::rows(cap.coef))), (1, 2))
    cap.coef = cap.coef[rows_new, ]
    cap.coef_txt = cap.coef_txt[rows_new, ]
}

void _cdid_nuis_store_pending(struct cdid_nuis_cap scalar cap)
{
    pointer() scalar p
    p = findexternal("__cdid_nuis_pending")
    if (p == NULL) p = crexternal("__cdid_nuis_pending")
    *p = cap
}

/* Record the status of the observation logits on a route that does not
   report them. */
void cdid_nuis_mark_selection(string scalar status, string scalar code, string scalar reason)
{
    struct cdid_nuis_cap scalar cap
    cap.comp = J(4, 4, "")
    cap.comp[1, ] = ("selection", status, code, reason)
    _cdid_nuis_store_pending(cap)
}

/* Complete the capture after posting: fill the status of every component and
   bind the capture to the posted results. */
void cdid_nuis_finish(string scalar token, string scalar mode, string scalar method,
    string scalar selection, string scalar selmethod, | string scalar clustervar)
{
    struct cdid_nuis_cap scalar cap
    pointer() scalar p
    p = findexternal("__cdid_nuis_pending")
    if (p != NULL) {
        cap = *p
        rmexternal("__cdid_nuis_pending")
    }
    else cap.comp = J(4, 4, "")
    cap.token = token
    cap.clustervar = (args() >= 6 ? clustervar : "")
    cap.mode = mode
    cap.method = method
    cap.selection = selection
    cap.selmethod = (selection == "none" ? "" : selmethod)
    if (cdid_xp_requested() & cols(cap.roles) == 3 & rows(cap.comp) == 4) {
        if (cap.comp[2, 2] != "" & cap.comp[3, 2] != "" & cap.comp[4, 2] != "") cap.mode = "chained"
        else cap.mode = "crosspair"
    }
    else if (cdid_xp_requested()) cap.mode = "crosspair"
    if (cap.mode == "crosspair") {
        cap.mode = "crosspair"
        cap.comp = ("selection", "unavailable", "cross_pair_route",
            "Nuisance-model reports are not available for this cross-pair route.") \
            ("outcome", "unavailable", "cross_pair_route", "") \
            ("treatment", "unavailable", "cross_pair_route", "") \
            ("imputation", "unavailable", "cross_pair_route", "")
    }
    else if (mode == "ordinary") {
        cap.comp = ("selection", "unavailable", "ordinary_mode",
            "Ordinary mode is computed by csdid; cdid reports no nuisance models in ordinary mode.") \
            ("outcome", "unavailable", "ordinary_mode", "") \
            ("treatment", "unavailable", "ordinary_mode", "") \
            ("imputation", "unavailable", "ordinary_mode", "")
    }
    else {
        if (cap.comp[1, 2] == "") {
            if (selection == "none") cap.comp[1, ] = ("selection", "not_applicable", "no_selection_model",
                "selection(none) fits no observation model.")
            else if (method == "custom") cap.comp[1, ] = ("selection", "unavailable", "custom_callback_route",
                "A user-written estimator uses the selection fit through its own callback; its observation logits are not reported.")
            else if (selmethod == "imputation") cap.comp[1, ] = ("selection", "not_applicable", "pure_imputation",
                "selmethod(imputation) fits no observation model.")
            else cap.comp[1, ] = ("selection", "not_captured", "outside_current_slice",
                "Observation logits of this specification are not available in cdid_nuisance.")
        }
        if (cap.comp[2, 2] == "") {
            if (method == "ipw" | method == "custom") cap.comp[2, ] = ("outcome", "not_applicable", "no_outcome_regression",
                "method(" + method + ") fits no outcome regression.")
            else cap.comp[2, ] = ("outcome", "not_captured", "outside_current_slice",
                "Outcome-change regressions of method(" + method + ") are not available in cdid_nuisance.")
        }
        if (cap.comp[3, 2] == "") {
            if (method == "custom") cap.comp[3, ] = ("treatment", "unavailable", "custom_callback",
                "A user-written estimator supplies its own weights; cdid holds no fitted treatment model for it.")
            else if (method == "reg") cap.comp[3, ] = ("treatment", "not_applicable", "no_treatment_propensity",
                "method(reg) fits no treatment propensity.")
            else cap.comp[3, ] = ("treatment", "not_captured", "outside_current_slice",
                "Treatment propensity logits are not available in cdid_nuisance.")
        }
        if (cap.comp[4, 2] == "") {
            if (selection != "none" & (selmethod == "imputation" | selmethod == "dr"))
                cap.comp[4, ] = ("imputation", "not_captured", "outside_current_slice",
                    "Actual-change imputation regressions are not available in cdid_nuisance.")
            else cap.comp[4, ] = ("imputation", "not_applicable", "no_imputation",
                (selection == "none" ? "selection(none)" : "selmethod(" + selmethod + ")") + " fits no actual-change imputation.")
        }
    }
    _cdid_nuis_reg_put(token, cap)
    cdid_nuis_begin(0)
}

void _cdid_nuis_status_line(string rowvector r)
{
    printf("{txt}  %s: {res}%s{txt}", r[1], subinstr(r[2], "_", " "))
    if (r[4] != "") printf("{txt}. %s", r[4])
    printf("\n")
}

string scalar _cdid_nuis_smcl(string scalar s)
{
    return(subinstr(subinstr(subinstr(s, "{", char(1)), "}", "{c )-}"), char(1), "{c -(}"))
}

/* Display of the augmented roles (outcome, treatment, imputation). */
string scalar _cdid_nuis_role_title(string scalar role)
{
    if (role == "outcome") return("Counterfactual outcome regressions on the reconstructed response")
    if (role == "treatment") return("Treatment propensity logits on each risk set")
    return("Imputation regressions of the actual outcome changes")
}

void _cdid_nuis_role_header(struct cdid_nuis_role scalar R, string scalar wtext, string scalar se_line)
{
    if (R.title != "") {
        printf("\n{txt}%s, %s\n", R.title, (R.pooling == "" ? "" : R.pooling + (strpos(R.pooling, " ") ? "" : " pooling")))
        printf("{txt}  response:   %s\n", _cdid_nuis_smcl(R.response))
        printf("{txt}  sample:     %s\n", _cdid_nuis_smcl(R.sample))
        printf("{txt}  weights:    %s\n", _cdid_nuis_smcl(R.weights))
        printf("{txt}  %s\n", se_line)
        return
    }
    printf("\n{txt}%s (augmented correction), %s\n", _cdid_nuis_role_title(R.role),
        (R.pooling == "risk_set" ? "one logit per target cohort and control set" : R.pooling + " pooling"))
    if (R.role == "outcome") {
        printf("{txt}  response:   %s\n", _cdid_nuis_smcl("Z_i(k;s,t) = j_c + S_i(s,t) (Y_it - Y_is - j_c) / max(q_c, eps) if the pair is observed, j_c otherwise"))
        printf("{txt}  sample:     every unit of each control cohort k of the block\n")
        printf("{txt}  weights:    %s, times the repeat rate a_k of the control cohort\n", _cdid_nuis_smcl(wtext))
    }
    else if (R.role == "treatment") {
        printf("{txt}  response:   %s\n", _cdid_nuis_smcl("D_i = 1 if G_i = g, on the units with G_i in {g} or in a control cohort"))
        printf("{txt}  sample:     every unit of the target and control cohorts\n")
        printf("{txt}  weights:    %s, times the repeat rate of the unit's cohort\n", _cdid_nuis_smcl(wtext))
    }
    else {
        printf("{txt}  response:   Y_it - Y_is on the units with both outcomes observed (S_i(s,t) = 1)\n")
        printf("{txt}  sample:     units of the cell's cohort with both outcomes observed\n")
        printf("{txt}  weights:    %s\n", _cdid_nuis_smcl(wtext))
    }
    printf("{txt}  %s\n", se_line)
}

void _cdid_nuis_role_index(struct cdid_nuis_role scalar R)
{
    real scalar b
    printf("\n{txt}%s\n", (R.title != "" ? R.title : _cdid_nuis_role_title(R.role)))
    printf("{txt}  {lalign 7:Block} {lalign 34:Label} {ralign 6:Param} {ralign 8:Rows} {ralign 7:Compar}\n")
    printf("{txt}  {hline 66}\n")
    for (b = 1; b <= rows(R.blocks); b++) {
        printf("{txt}  {lalign 7:%s} {lalign 34:%s} {res}{ralign 6:%g} {ralign 8:%g} {ralign 7:%g}\n",
            R.block_txt[b, 1], abbrev(subinstr(subinstr(R.block_txt[b, 2], "{", "("), "}", ")"), 34),
            R.blocks[b, 1], R.blocks[b, 4], R.blocks[b, 7])
    }
}

void _cdid_nuis_role_block(struct cdid_nuis_role scalar R, real scalar b, string scalar se_line, real scalar post)
{
    real matrix T
    real colvector ix
    real scalar k, j
    string colvector terms, names
    ix = (rows(R.coef) ? select((1::rows(R.coef)), R.coef[, 1] :== b) : J(0, 1, .))
    k = rows(ix)
    printf("\n{txt}Block {res}%s{txt}: %s\n", R.block_txt[b, 1], _cdid_nuis_smcl(R.block_txt[b, 2]))
    printf("{txt}  %s\n", _cdid_nuis_smcl(R.block_txt[b, 3]))
    printf("{txt}  %g fitting rows, sum of weights %g; used by %g comparisons\n", R.blocks[b, 4], R.blocks[b, 5], R.blocks[b, 7])
    if (cols(R.block_txt) >= 6) if (R.block_txt[b, 6] != "") printf("{txt}  slopes dropped as exact aliases on the prediction rows: %s\n", R.block_txt[b, 6])
    printf("{txt}{hline 30}{c TT}{hline 46}\n")
    printf("{txt}{ralign 29:Term} {c |} {ralign 11:Coef.} {ralign 11:Std. err.} {ralign 9:z} {ralign 9:P>|z|}\n")
    printf("{txt}{hline 30}{c +}{hline 46}\n")
    for (j = 1; j <= k; j++) {
        printf("{txt}{ralign 29:%s} {c |} {res}{ralign 11:%9.0g} {ralign 11:%9.0g} {ralign 9:%7.2f} {ralign 9:%6.3f}\n",
            abbrev(subinstr(subinstr(R.coef_txt[ix[j], 1], "{", "("), "}", ")"), 29), R.coef[ix[j], 3], R.coef[ix[j], 4], R.coef[ix[j], 5], R.coef[ix[j], 6])
    }
    printf("{txt}{hline 30}{c BT}{hline 46}\n")
    printf("{txt}%s\n", se_line)
    if (!post | k == 0) return
    T = J(9, k, .)
    T[1, ] = R.coef[ix, 3]'
    T[2, ] = R.coef[ix, 4]'
    T[3, ] = R.coef[ix, 5]'
    T[4, ] = R.coef[ix, 6]'
    T[5, ] = (R.coef[ix, 3] - invnormal(.975) :* R.coef[ix, 4])'
    T[6, ] = (R.coef[ix, 3] + invnormal(.975) :* R.coef[ix, 4])'
    T[8, ] = J(1, k, invnormal(.975))
    T[9, ] = J(1, k, 0)
    names = J(k, 1, "")
    terms = R.coef_txt[ix, 1]
    for (j = 1; j <= k; j++) {
        if (R.coef_txt[ix[j], 2] == "local_intercept" | R.coef_txt[ix[j], 2] == "intercept")
            names[j] = "_cons" + (sum(R.coef_txt[ix, 2] :== "local_intercept") > 1 ? "_" + strofreal(j) : "")
        else names[j] = strtoname(R.coef_txt[ix[j], 1])
    }
    st_matrix("r(table)", T)
    st_matrixrowstripe("r(table)", (J(9, 1, ""), ("b" \ "se" \ "z" \ "pvalue" \ "ll" \ "ul" \ "df" \ "crit" \ "eform")))
    st_matrixcolstripe("r(table)", (J(k, 1, ""), names))
    st_global("r(terms)", invtokens(("`" :+ `"""' :+ terms :+ `"""' :+ "'")'))
    st_global("r(block)", R.block_txt[b, 1])
    st_global("r(label)", R.block_txt[b, 2])
    st_global("r(equation)", R.block_txt[b, 3])
    if (cols(R.block_txt) >= 6) {
        st_global("r(kept_slopes)", R.block_txt[b, 5])
        st_global("r(aliased_slopes)", R.block_txt[b, 6])
    }
    st_numscalar("r(n_fitting_rows)", R.blocks[b, 4])
    st_numscalar("r(n_comparisons)", R.blocks[b, 7])
}

void _cdid_nuis_role_show(struct cdid_nuis_role scalar R, string scalar block, string scalar se_line,
    string scalar wtext)
{
    real scalar nb, b, j
    nb = rows(R.blocks)
    if (block == "") {
        _cdid_nuis_role_header(R, wtext, se_line)
        for (b = 1; b <= nb; b++) _cdid_nuis_role_block(R, b, se_line, 0)
        st_matrix("r(blocks)", R.blocks[, 1..7])
        st_matrixrowstripe("r(blocks)", (J(nb, 1, ""), R.block_txt[, 1]))
        st_matrixcolstripe("r(blocks)", (J(7, 1, ""), ("n_parameters" \ "n_intercepts" \ "n_slopes" \
            "n_fitting_rows" \ "weight_sum" \ "n_contributing_clusters" \ "n_comparisons")))
        st_numscalar("r(n_blocks)", nb)
        return
    }
    b = .
    for (j = 1; j <= nb; j++) {
        if (strlower(R.block_txt[j, 1]) == strlower(strtrim(block)) | R.block_txt[j, 2] == strtrim(block)) b = j
    }
    if (b >= . & strtoreal(block) < . & strtoreal(block) >= 1 & strtoreal(block) <= nb & strtoreal(block) == floor(strtoreal(block)))
        b = strtoreal(block)
    if (b >= .) {
        errprintf("block(%s) not found; use a block identifier such as %s, its number, or its label\n", block, R.block_txt[1, 1])
        exit(198)
    }
    _cdid_nuis_role_block(R, b, se_line, 1)
}

/* Selection blocks proper: the roster probabilities Pr(U | cohort) of the fitting rosters A
   and B are kept in the last row with identifier roster and are reported in their own table. */
real scalar _cdid_nuis_nsel(struct cdid_nuis_cap scalar cap)
{
    real scalar nb
    nb = rows(cap.blocks)
    if (nb) if (cap.block_txt[nb, 1] == "roster") nb--
    return(nb)
}

/* Table of the roster probabilities, also returned in r(roster) (estimate and standard error
   by cohort). */
void _cdid_nuis_roster(struct cdid_nuis_cap scalar cap, string scalar se_line)
{
    real colvector ix
    real matrix T
    real scalar b, j, k
    string colvector names
    b = rows(cap.blocks)
    if (b == _cdid_nuis_nsel(cap)) return
    ix = (rows(cap.coef) ? select((1::rows(cap.coef)), cap.coef[, 1] :== b) : J(0, 1, .))
    k = rows(ix)
    printf("\n{txt}Roster probabilities: %s\n", _cdid_nuis_smcl(cap.block_txt[b, 3]))
    printf("{txt}  %g units in the cohorts of the frame; used to rescale the observation probabilities of the fitting roster\n", cap.blocks[b, 5])
    printf("{txt}{hline 30}{c TT}{hline 46}\n")
    printf("{txt}{ralign 29:Term} {c |} {ralign 11:Coef.} {ralign 11:Std. err.} {ralign 9:z} {ralign 9:P>|z|}\n")
    printf("{txt}{hline 30}{c +}{hline 46}\n")
    for (j = 1; j <= k; j++) {
        printf("{txt}{ralign 29:%s} {c |} {res}{ralign 11:%9.0g} {ralign 11:%9.0g} {ralign 9:%7.2f} {ralign 9:%6.3f}\n",
            abbrev(cap.coef_txt[ix[j], 1], 29), cap.coef[ix[j], 3], cap.coef[ix[j], 4], cap.coef[ix[j], 5], cap.coef[ix[j], 6])
    }
    printf("{txt}{hline 30}{c BT}{hline 46}\n")
    printf("{txt}%s\n", se_line)
    if (k == 0) return
    T = cap.coef[ix, (3, 4)]
    names = cap.coef_txt[ix, 1]
    st_matrix("r(roster)", T)
    st_matrixrowstripe("r(roster)", (J(k, 1, ""), strtoname(names)))
    st_matrixcolstripe("r(roster)", (J(2, 1, ""), ("b" \ "se")))
    st_global("r(roster_terms)", invtokens(("`" :+ `"""' :+ names :+ `"""' :+ "'")'))
}

/* Display and return the captured reports. */
void cdid_nuis_show(string scalar token, string scalar component, string scalar block,
    string scalar clustervar, string scalar wtext)
{
    struct cdid_nuis_cap scalar cap
    real scalar i, b, j, n, nb, avail
    real matrix sel, T
    string scalar comp, se_line
    string colvector terms

    st_rclear()
    comp = strlower(strtrim(component))
    if (comp != "" & comp != "selection" & comp != "outcome" & comp != "treatment" & comp != "imputation") {
        errprintf("component must be selection, outcome, treatment, or imputation\n")
        exit(198)
    }
    if (!_cdid_nuis_reg_get(token, cap)) {
        printf("{txt}No nuisance models were kept with the current results; refit with the {cmd:keepnuisance} option.\n")
        st_global("r(status)", "not_captured")
        return
    }
    if (cap.mode == "crosspair") {
        printf("{txt}Nuisance models are not reported for this {cmd:crosspair} route.\n")
        st_global("r(status)", "unavailable")
        st_global("r(reason)", cap.comp[1, 4])
        return
    }
    if (cap.mode == "ordinary") {
        printf("{txt}Nuisance models are not reported in ordinary mode, which {cmd:cdid} computes with {cmd:csdid}.\n")
        st_global("r(status)", "unavailable")
        st_global("r(reason)", cap.comp[1, 4])
        return
    }
    if (comp == "" & block != "") comp = "selection"
    st_global("r(method)", cap.method)
    st_global("r(selmethod)", cap.selmethod)
    for (i = 1; i <= 4; i++) {
        st_global("r(status_" + cap.comp[i, 1] + ")", cap.comp[i, 2])
    }
    avail = (cap.comp[1, 2] == "available")
    if (cap.clustered) se_line = "SE: cluster-robust by " + (cap.clustervar != "" ? cap.clustervar : clustervar) + " (" + strofreal(cap.n_groups) +
        " clusters): unit scores summed within cluster, centered, divided by N^2; no small-sample adjustment."
    else se_line = "SE: one score contribution per unit (N = " + strofreal(cap.N) +
        "), centered, divided by N^2; no small-sample adjustment."

    if (comp == "") {
        printf("\n{txt}Nuisance models of {cmd:cdid}, method({res}%s{txt})%s\n", cap.method,
            (cap.selmethod == "" ? "" : ", selmethod({res}" + cap.selmethod + "{txt})"))
        if (avail) {
            _cdid_nuis_header(cap, wtext, se_line)
            nb = _cdid_nuis_nsel(cap)
            printf("\n{txt}  {lalign 7:Block} {lalign 26:Label} {ralign 6:Cells} {ralign 6:Param} {ralign 8:Rows} {ralign 6:Units} {ralign 7:Compar}\n")
            printf("{txt}  {hline 72}\n")
            for (b = 1; b <= nb; b++) {
                printf("{txt}  {lalign 7:%s} {lalign 26:%s} {res}{ralign 6:%g} {ralign 6:%g} {ralign 8:%g} {ralign 6:%g} {ralign 7:%g}\n",
                    cap.block_txt[b, 1], abbrev(subinstr(subinstr(cap.block_txt[b, 2], "{", "("), "}", ")"), 26), cap.blocks[b, 1], cap.blocks[b, 2],
                    cap.blocks[b, 5], cap.blocks[b, 6], cap.blocks[b, 13])
            }
            T = cap.blocks[1..nb, (1, 2, 5, 6, 7, 8, 10, 11, 12, 13)]
            st_matrix("r(blocks)", T)
            st_matrixrowstripe("r(blocks)", (J(nb, 1, ""), cap.block_txt[1..nb, 1]))
            st_matrixcolstripe("r(blocks)", (J(10, 1, ""), ("n_cells" \ "n_parameters" \ "n_fitting_rows" \
                "n_units" \ "n_fitting_clusters" \ "unit_weight_sum" \ "iterations" \ "rank" \ "n_pair_cohorts" \ "n_comparisons")))
            st_numscalar("r(n_blocks)", nb)
            st_numscalar("r(N)", cap.N)
            if (nb < rows(cap.blocks)) printf("\n{txt}  Roster probabilities Pr(U | cohort): %g cohorts; {cmd:cdid_nuisance selection} reports them\n",
                cap.blocks[rows(cap.blocks), 2])
        }
        for (i = 1; i <= cols(cap.roles); i++) if (cap.comp[i + 1, 2] == "available") _cdid_nuis_role_index(cap.roles[i])
        printf("\n{txt}Components:\n")
        for (i = 1; i <= 4; i++) _cdid_nuis_status_line(cap.comp[i, ])
        if (avail & cols(cap.roles) == 3) printf("\n{txt}Coefficients: {cmd:cdid_nuisance} {it:component} [{cmd:, block(}{it:id}{cmd:)}], with {it:component} selection, outcome, treatment, or imputation\n")
        else if (avail) printf("\n{txt}Coefficients: {cmd:cdid_nuisance selection} or {cmd:cdid_nuisance selection, block(%s)}\n", cap.block_txt[1, 1])
        st_global("r(component)", "")
        return
    }

    i = (comp == "selection" ? 1 : (comp == "outcome" ? 2 : (comp == "treatment" ? 3 : 4)))
    st_global("r(component)", comp)
    st_global("r(status)", cap.comp[i, 2])
    if (cap.comp[i, 2] != "available") {
        printf("{txt}Nuisance component %s: {res}%s{txt}", comp, subinstr(cap.comp[i, 2], "_", " "))
        if (cap.comp[i, 4] != "") printf("{txt}. %s", cap.comp[i, 4])
        printf("\n")
        st_global("r(reason)", cap.comp[i, 4])
        return
    }
    if (i > 1) {
        if (cols(cap.roles) != 3) {
            errprintf("the %s report is not available\n", comp)
            exit(498)
        }
        _cdid_nuis_role_show(cap.roles[i - 1], block, se_line, wtext)
        return
    }
    nb = rows(cap.blocks)
    if (block == "") {
        _cdid_nuis_header(cap, wtext, se_line)
        for (b = 1; b <= _cdid_nuis_nsel(cap); b++) _cdid_nuis_block(cap, b, se_line, 0)
        _cdid_nuis_roster(cap, se_line)
        st_numscalar("r(n_blocks)", _cdid_nuis_nsel(cap))
        return
    }
    b = .
    for (j = 1; j <= nb; j++) {
        if (strlower(cap.block_txt[j, 1]) == strlower(strtrim(block)) | cap.block_txt[j, 2] == strtrim(block)) b = j
    }
    if (b >= . & strtoreal(block) < . & strtoreal(block) >= 1 & strtoreal(block) <= nb & strtoreal(block) == floor(strtoreal(block)))
        b = strtoreal(block)
    if (b >= .) {
        errprintf("block(%s) not found; use a block identifier such as %s, its number, or its label\n", block, cap.block_txt[1, 1])
        exit(198)
    }
    _cdid_nuis_block(cap, b, se_line, 1)
}

void _cdid_nuis_header(struct cdid_nuis_cap scalar cap, string scalar wtext, string scalar se_line)
{
    string scalar resp, templ, sl
    real scalar j
    sl = ""
    for (j = 1; j <= cols(cap.slope_names); j++) sl = sl + (j > 1 ? " " : "") + cap.slope_names[j]
    if (sl == "") sl = "none"
    if (cap.sel_model == "pair_logit") {
        resp = "S_i(s,t) = 1 if Y_is and Y_it are both observed"
        templ = "logit Pr(S_i(s,t) = 1 | X_i, G_i = g, U_i = 1)"
    }
    else {
        resp = "R_it = 1 if Y_it is observed; Pr(S_i(s,t) = 1) = Pr(R_is = 1) Pr(R_it = 1)"
        templ = "logit Pr(R_it = 1 | X_i, G_i = g, U_i = 1)"
    }
    if (cap.mode == "static") {
        printf("\n{txt}Observation logits: {res}%s{txt}, {res}%s{txt} pooling\n", cap.sel_model, cap.sel_pooling)
        if (cap.sel_model == "pair_logit")
            printf("{txt}  model:      logit Pr(both outcomes of the pair observed | role, X(t-1)) = intercepts + slopes\n")
        else printf("{txt}  model:      logit Pr(outcome observed at the date | role, X(t-1)) = intercepts + slopes; pair probability = product of the dates\n")
        printf("{txt}  covariates: %s\n", sl)
        printf("{txt}  sample:     switchers or stayers of each transition (U = %s); weights: unit weights\n", cap.sel_universe)
        printf("{txt}  %s\n", se_line)
        return
    }
    printf("\n{txt}Observation logits: {res}%s{txt}, {res}%s{txt} pooling\n", cap.sel_model, cap.sel_pooling)
    printf("{txt}  model:      %s = intercepts + slopes on the covariates\n", templ)
    printf("{txt}  response:   %s\n", resp)
    printf("{txt}  covariates: %s\n", sl)
    printf("{txt}  sample:     U = %s; weights: %s\n", cap.sel_universe, wtext)
    printf("{txt}  %s\n", se_line)
}

void _cdid_nuis_block(struct cdid_nuis_cap scalar cap, real scalar b, string scalar se_line, real scalar post)
{
    real matrix sel, T, C
    real colvector ix
    real scalar k, j
    string colvector terms, names
    ix = (rows(cap.coef) ? select((1::rows(cap.coef)), cap.coef[, 1] :== b) : J(0, 1, .))
    k = rows(ix)
    printf("\n{txt}Block {res}%s{txt}: %s\n", cap.block_txt[b, 1], _cdid_nuis_smcl(cap.block_txt[b, 2]))
    printf("{txt}  %s\n", _cdid_nuis_smcl(cap.block_txt[b, 3]))
    if (cap.block_txt[b, 1] == "roster") printf("{txt}  %g units in the cohorts of the frame\n", cap.blocks[b, 5])
    else printf("{txt}  %g units, %g fitting rows; converged in %g iterations, rank %g; used by %g comparisons\n",
        cap.blocks[b, 6], cap.blocks[b, 5], cap.blocks[b, 10], cap.blocks[b, 11], cap.blocks[b, 13])
    printf("{txt}{hline 30}{c TT}{hline 46}\n")
    printf("{txt}{ralign 29:Term} {c |} {ralign 11:Coef.} {ralign 11:Std. err.} {ralign 9:z} {ralign 9:P>|z|}\n")
    printf("{txt}{hline 30}{c +}{hline 46}\n")
    for (j = 1; j <= k; j++) {
        printf("{txt}{ralign 29:%s} {c |} {res}{ralign 11:%9.0g} {ralign 11:%9.0g} {ralign 9:%7.2f} {ralign 9:%6.3f}\n",
            abbrev(subinstr(subinstr(cap.coef_txt[ix[j], 1], "{", "("), "}", ")"), 29), cap.coef[ix[j], 3], cap.coef[ix[j], 4], cap.coef[ix[j], 5], cap.coef[ix[j], 6])
    }
    printf("{txt}{hline 30}{c BT}{hline 46}\n")
    printf("{txt}%s\n", se_line)
    if (!post | k == 0) return
    T = J(9, k, .)
    T[1, ] = cap.coef[ix, 3]'
    T[2, ] = cap.coef[ix, 4]'
    T[3, ] = cap.coef[ix, 5]'
    T[4, ] = cap.coef[ix, 6]'
    T[5, ] = (cap.coef[ix, 3] - invnormal(.975) :* cap.coef[ix, 4])'
    T[6, ] = (cap.coef[ix, 3] + invnormal(.975) :* cap.coef[ix, 4])'
    T[8, ] = J(1, k, invnormal(.975))
    T[9, ] = J(1, k, 0)
    names = J(k, 1, "")
    terms = cap.coef_txt[ix, 1]
    for (j = 1; j <= k; j++) {
        if (cap.coef_txt[ix[j], 2] == "local_intercept") names[j] = "_cons" + (sum(cap.coef_txt[ix, 2] :== "local_intercept") > 1 ? "_" + strofreal(j) : "")
        else names[j] = strtoname(cap.coef_txt[ix[j], 1])
    }
    st_matrix("r(table)", T)
    st_matrixrowstripe("r(table)", (J(9, 1, ""), ("b" \ "se" \ "z" \ "pvalue" \ "ll" \ "ul" \ "df" \ "crit" \ "eform")))
    st_matrixcolstripe("r(table)", (J(k, 1, ""), names))
    st_global("r(terms)", invtokens(("`" :+ `"""' :+ terms :+ `"""' :+ "'")'))
    st_global("r(block)", cap.block_txt[b, 1])
    st_global("r(label)", cap.block_txt[b, 2])
    st_global("r(equation)", cap.block_txt[b, 3])
    st_numscalar("r(n_units)", cap.blocks[b, 6])
    st_numscalar("r(n_fitting_rows)", cap.blocks[b, 5])
    st_numscalar("r(n_comparisons)", cap.blocks[b, 13])
}

end
