*! Public ordinary ATT custom-estimator adapter, current StataNow19.5.
*! GPL-3 adaptation of R cdid att_gt.R + did2.5.1 faster_mode=TRUE.
*! Uses cdid_co_call_cell() and existing ordinary inference/state/aggregation contracts.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_coatt(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_ord_custom_owner {
    real colvector rows, y, t, g, raw_w, w, iv, row_unit, tlist, glist
    real colvector unit_id, source_group, raw_unit_w, source_weight
    real colvector cluster_labels, cluster_positions
    string colvector scluster
    real matrix identity_map, raw_x, X
    string rowvector raw_xnames, xcols
    string scalar formula, control, base_period, fix_weights, cluster_name
    string scalar callback_name, payload_name
    real scalar panel, allow, eff_panel, eff_unbal, true_rcs, intercept
    real scalar anticipation, norm_factor, n_units, C, max_entries, cluster_is_str
}
void _cdid_oca_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_oca_ix(real matrix x)
{
    return(vec(selectindex(x)))
}
void cdid_ord_custom_intake(string scalar owner_name, string scalar sample,
    string scalar evalmask, string scalar y, string scalar time, string scalar group,
    string scalar id, string scalar weight, string scalar clusters, string scalar formula,
    real scalar intercept, real scalar panel, real scalar allow, string scalar control,
    string scalar base_period, real scalar anticipation, string scalar fixweights,
    string scalar callback_name, string scalar payload_name, real scalar workspace)
{
    struct cdid_ord_custom_owner scalar O
    pointer() scalar p
    string rowvector xnames, cn
    string colvector scl_vec, u_scl, slabels
    real colvector rr, y_vec, t_vec, g_vec, raw_w, w_vec, iv, cl_vec, ok, keep
    real colvector tl, gl, uq, ord, g_sort, pre_rowid, u_rows, u_cl, u_ix, distinct
    real matrix raw, info
    real scalar n, j, k, cap, latest, cutoff_t, first, is_bal, nT, first_r, lo, hi, mid
    if (findexternal(owner_name) != NULL) _cdid_oca_fail("Temporary ordinary custom object already exists.", 498)
    cdid_custom_validate_symbols(callback_name, payload_name)
    if (control != "nevertreated" & control != "notyettreated")
        _cdid_oca_fail("controlgroup() must be nevertreated or notyettreated.", 198)
    if (base_period != "varying" & base_period != "universal")
        _cdid_oca_fail("baseperiod() must be varying or universal.", 198)
    if (anticipation < 0 | anticipation >= . | anticipation != floor(anticipation))
        _cdid_oca_fail("anticipation() must be a nonnegative integer.", 198)
    if (fixweights != "" & fixweights != "varying" & fixweights != "base_period" & fixweights != "first_period")
        _cdid_oca_fail("fixweights() must be omitted, first_period, base_period, or varying.", 198)
    if (!panel & (fixweights == "base_period" | fixweights == "first_period"))
        _cdid_oca_fail("fixed unit anchors are not supported for true repeated cross sections.", 198)
    if (panel & fixweights == "varying")
        _cdid_oca_fail("fixweights(varying) is not supported with custom estimators on panel data.", 198)
    if (panel & id == "")
        _cdid_oca_fail("ordinary panel data require ivar().", 198)
    if (workspace < 2048 | workspace >= . | workspace != floor(workspace))
        _cdid_oca_fail("workspace() must be an explicit integer byte limit.", 198)
    cap = floor(workspace / 8)
    xnames = cdid_ref_formula_names(formula)
    cn = tokens(clusters)
    if (id != "" & cols(cn) > 0) cn = select(cn, cn :!= id)
    if (cols(cn) > 1)
        _cdid_oca_fail("At most one cluster variable beyond ivar() is supported.", 198)
    rr = _cdid_oca_ix(st_data(., sample) :!= 0); n = rows(rr)
    if (n == 0) _cdid_oca_fail("No requested observations.", 2000)
    if (2 * st_nobs() + n * (cols(xnames) + 14) > cap)
        _cdid_oca_fail("Raw intake exceeds the explicit workspace budget.", 908)
    cdid_ordinary_input_guard(id, time, group, sample)
    t_vec = st_data(rr, time); g_vec = st_data(rr, group); y_vec = st_data(rr, y)
    raw_w = (weight == "" ? J(n, 1, 1) : st_data(rr, weight))
    if (id != "") {
        if (!st_isnumvar(id)) _cdid_oca_fail("ivar() must be numeric in ordinary mode.", 198)
        iv = st_data(rr, id)
    }
    else iv = (1::n)
    O.cluster_is_str = 0
    cl_vec = J(n, 1, .); scl_vec = J(n, 1, "")
    if (cols(cn) == 1) {
        if (st_isstrvar(cn[1])) {
            O.cluster_is_str = 1
            scl_vec = st_sdata(rr, cn[1])
            if (id != "") {
                ok = (iv :< .) :& (scl_vec :!= "")
                if (any(ok)) {
                    keep = _cdid_oca_ix(ok)
                    ord = order(iv[keep], 1); keep = keep[ord]
                    for (j = 2; j <= rows(keep); j++)
                        if (iv[keep[j]] == iv[keep[j - 1]] & scl_vec[keep[j]] != scl_vec[keep[j - 1]])
                            _cdid_oca_fail("Time-varying ordinary cluster variables are not supported.", 459)
                }
            }
        }
        else {
            cl_vec = st_data(rr, cn[1])
            if (id != "") {
                ok = (iv :< .) :& (cl_vec :< .)
                if (any(ok)) {
                    keep = _cdid_oca_ix(ok)
                    if (rows(uniqrows((iv[keep], cl_vec[keep]))) != rows(uniqrows(iv[keep])))
                        _cdid_oca_fail("Time-varying ordinary cluster variables are not supported.", 459)
                }
            }
        }
    }
    ok = (y_vec :< .) :& (t_vec :< .) :& (g_vec :< .) :& (raw_w :< .) :& (iv :< .)
    if (cols(cn) == 1) ok = ok :& (O.cluster_is_str ? (scl_vec :!= "") : (cl_vec :< .))
    if (cols(xnames)) {
        raw = st_data(rr, xnames)
        ok = ok :& (rowmissing(raw) :== 0)
    }
    else raw = J(n, 0, .)
    ok = ok :& (st_data(rr, evalmask) :!= 0)
    rr = select(rr, ok); y_vec = select(y_vec, ok); t_vec = select(t_vec, ok)
    g_vec = select(g_vec, ok); raw_w = select(raw_w, ok); iv = select(iv, ok)
    cl_vec = select(cl_vec, ok); scl_vec = select(scl_vec, ok); raw = select(raw, ok)
    if (!rows(rr)) _cdid_oca_fail("All requested observations were missing or nonfinite.", 2000)
    if (any(raw_w :< 0) | mean(raw_w) <= 0)
        _cdid_oca_fail("Sampling weights must be nonnegative with a positive mean.", 498)
    /* did2.5.1 faster-mode standardization normalizes BEFORE support/balance settling. */
    O.norm_factor = 1 / mean(raw_w)
    w_vec = raw_w * O.norm_factor
    tl = uniqrows(sort(t_vec, 1))
    g_vec = (g_vec :> (max(tl) + anticipation)) :* 0 :+ (g_vec :<= (max(tl) + anticipation)) :* g_vec
    gl = uniqrows(sort(g_vec, 1)); latest = .
    if (!any(gl :== 0)) {
        latest = max(gl); cutoff_t = latest - anticipation
        keep = (t_vec :< cutoff_t)
        rr = select(rr, keep); y_vec = select(y_vec, keep); t_vec = select(t_vec, keep)
        g_vec = select(g_vec, keep); raw_w = select(raw_w, keep); w_vec = select(w_vec, keep)
        iv = select(iv, keep); cl_vec = select(cl_vec, keep); scl_vec = select(scl_vec, keep); raw = select(raw, keep)
        if (!rows(rr)) _cdid_oca_fail("No observations remain after the no-never-treated cutoff.", 2000)
        if (control == "nevertreated") g_vec = g_vec :* (g_vec :!= latest)
        tl = uniqrows(sort(t_vec, 1))
    }
    first = tl[1]
    keep = (g_vec :== 0) :| (g_vec :> (first + anticipation))
    rr = select(rr, keep); y_vec = select(y_vec, keep); t_vec = select(t_vec, keep)
    g_vec = select(g_vec, keep); raw_w = select(raw_w, keep); w_vec = select(w_vec, keep)
    iv = select(iv, keep); cl_vec = select(cl_vec, keep); scl_vec = select(scl_vec, keep); raw = select(raw, keep)
    if (!rows(rr)) _cdid_oca_fail("No valid groups remain after first-period treatment exclusions.", 498)
    tl = uniqrows(sort(t_vec, 1))
    gl = uniqrows(sort(select(g_vec, (g_vec :> 0) :& (g_vec :> (tl[1] + anticipation))), 1))
    if (control != "nevertreated" & !any(g_vec :== 0) & latest < .)
        gl = select(gl, gl :< latest)
    O.true_rcs = (!panel)
    O.eff_panel = panel
    O.eff_unbal = (panel & allow)
    if (O.eff_unbal) {
        is_bal = (rows(rr) == rows(uniqrows(iv)) * rows(tl))
        O.eff_unbal = !is_bal
        O.eff_panel = is_bal
        if (O.eff_unbal) O.true_rcs = 0
    }
    else if (O.eff_panel) {
        ord = order((iv, t_vec), (1, 2)); info = _cdid_ps_coatt(iv[ord], 1)
        keep = J(rows(rr), 1, 0)
        for (j = 1; j <= rows(info); j++) {
            u_ix = ord[|info[j, 1] \ info[j, 2]|]
            if (rows(u_ix) == rows(tl)) keep[u_ix] = J(rows(u_ix), 1, 1)
        }
        rr = select(rr, keep); y_vec = select(y_vec, keep); t_vec = select(t_vec, keep)
        g_vec = select(g_vec, keep); raw_w = select(raw_w, keep); w_vec = select(w_vec, keep)
        iv = select(iv, keep); cl_vec = select(cl_vec, keep); scl_vec = select(scl_vec, keep); raw = select(raw, keep)
        if (!rows(rr)) _cdid_oca_fail("All observations dropped while balancing the panel.", 498)
    }
    if (!rows(gl)) _cdid_oca_fail("No valid treated groups.", 498)
    if (rows(tl) < 2) _cdid_oca_fail("At least two time periods are required.", 498)
    uq = uniqrows(sort(g_vec, 1))
    for (j = 1; j <= rows(uq); j++) {
        if (sum(g_vec :== uq[j]) / rows(tl) < cols(xnames) + 5) {
            if (uq[j] == 0 & control == "nevertreated")
                _cdid_oca_fail("The never-treated group is too small; use notyettreated explicitly if appropriate.", 498)
        }
    }
    pre_rowid = (1::rows(rr))
    g_sort = (g_vec :== 0) :* (max(g_vec) + 1000) :+ g_vec
    if (O.eff_panel) ord = order((t_vec, g_sort, iv), (1, 2, 3))
    else ord = order((t_vec, g_sort, pre_rowid), (1, 2, 3))
    rr = rr[ord]; y_vec = y_vec[ord]; t_vec = t_vec[ord]; g_vec = g_vec[ord]
    raw_w = raw_w[ord]; w_vec = w_vec[ord]; iv = iv[ord]; cl_vec = cl_vec[ord]
    scl_vec = scl_vec[ord]; pre_rowid = pre_rowid[ord]; raw = raw[ord,]
    nT = rows(tl)
    if (O.eff_panel) {
        O.n_units = rows(rr) / nT
        u_rows = (1::O.n_units)
        O.unit_id = iv[u_rows]
        O.source_group = g_vec[u_rows]
        O.raw_unit_w = raw_w[u_rows]
        O.source_weight = w_vec[u_rows]
        O.identity_map = (O.unit_id, rr[u_rows], O.unit_id)
        u_cl = cl_vec[u_rows]; u_scl = scl_vec[u_rows]
        O.row_unit = J(rows(rr), 1, 0)
        for (k = 1; k <= nT; k++)
            O.row_unit[((k - 1) * O.n_units + 1)..(k * O.n_units)] = (1::O.n_units)
    }
    else if (O.eff_unbal) {
        O.row_unit = cdid_oc_first_positions(iv)
        O.n_units = max(O.row_unit)
        O.unit_id = O.source_group = O.raw_unit_w = O.source_weight = u_cl = J(O.n_units, 1, .)
        u_scl = J(O.n_units, 1, "")
        O.identity_map = J(O.n_units, 3, .)
        for (j = 1; j <= O.n_units; j++) {
            u_ix = _cdid_oca_ix(O.row_unit :== j)
            first_r = u_ix[1]
            O.unit_id[j] = iv[first_r]
            O.source_group[j] = g_vec[first_r]
            O.raw_unit_w[j] = mean(raw_w[u_ix])
            O.source_weight[j] = mean(w_vec[u_ix])
            O.identity_map[j,] = (O.unit_id[j], rr[first_r], O.unit_id[j])
            u_cl[j] = cl_vec[first_r]; u_scl[j] = scl_vec[first_r]
        }
    }
    else {
        O.n_units = rows(rr)
        O.row_unit = (1::O.n_units)
        O.unit_id = pre_rowid
        O.source_group = g_vec
        O.raw_unit_w = raw_w
        O.source_weight = w_vec
        O.identity_map = ((1::O.n_units), rr, (id == "" ? rr : iv))
        u_cl = cl_vec; u_scl = scl_vec
    }
    O.scluster = J(O.n_units, 1, "")
    if (cols(cn) == 0) {
        O.cluster_name = ""
        O.cluster_labels = (1::O.n_units)
        O.cluster_positions = (1::O.n_units)
        O.C = O.n_units
    }
    else {
        O.cluster_name = cn[1]
        if (O.cluster_is_str) {
            O.scluster = u_scl
            O.cluster_labels = (1::O.n_units)
            slabels = uniqrows(sort(u_scl, 1)); O.C = rows(slabels)
            O.cluster_positions = J(O.n_units, 1, .)
            for (j = 1; j <= O.n_units; j++)
                O.cluster_positions[j] = _cdid_oca_ix(slabels :== u_scl[j])[1]
        }
        else {
            O.cluster_labels = u_cl
            distinct = uniqrows(sort(u_cl, 1)); O.C = rows(distinct)
            O.cluster_positions = J(O.n_units, 1, .)
            for (j = 1; j <= O.n_units; j++) {
                lo = 1; hi = O.C
                while (lo <= hi) {
                    mid = floor((lo + hi) / 2)
                    if (distinct[mid] == u_cl[j]) break
                    if (distinct[mid] < u_cl[j]) lo = mid + 1
                    else hi = mid - 1
                }
                O.cluster_positions[j] = mid
            }
        }
    }
    O.rows = rr; O.y = y_vec; O.t = t_vec; O.g = g_vec; O.raw_w = raw_w; O.w = w_vec; O.iv = iv
    O.tlist = tl; O.glist = gl; O.raw_x = raw; O.raw_xnames = xnames
    O.formula = formula; O.intercept = intercept; O.panel = panel; O.allow = allow
    O.control = control; O.base_period = base_period; O.anticipation = anticipation
    O.fix_weights = fixweights; O.callback_name = callback_name; O.payload_name = payload_name
    O.max_entries = cap
    O.X = J(0, 0, .); O.xcols = J(1, 0, "")
    p = crexternal(owner_name)
    if (p == NULL) _cdid_oca_fail("Cannot create temporary ordinary custom object.", 498)
    *p = O
}
void cdid_ord_custom_stage(string scalar owner_name)
{
    pointer() scalar p
    struct cdid_ord_custom_owner scalar O
    real scalar j
    p = findexternal(owner_name)
    if (p == NULL) _cdid_oca_fail("Ordinary custom object is absent.", 498)
    O = *p
    if (st_nobs() != 0 | st_nvar() != 0) _cdid_oca_fail("Ordinary custom staging frame must be empty.", 498)
    st_addobs(rows(O.rows))
    for (j = 1; j <= cols(O.raw_x); j++) {
        (void) st_addvar("double", O.raw_xnames[j])
        st_store(., O.raw_xnames[j], O.raw_x[, j])
    }
}
void cdid_ord_custom_encode(string scalar owner_name, string scalar generated,
    string scalar expanded, real rowvector omitted)
{
    pointer() scalar p
    struct cdid_ord_custom_owner scalar O
    struct cdid_reference_inputs scalar S
    struct cdid_ref_design scalar D
    p = findexternal(owner_name)
    if (p == NULL) _cdid_oca_fail("Ordinary custom object is absent.", 498)
    O = *p
    D.keys = (1::st_nobs()); D.intercept = O.intercept; D.encoded = 0; D.formula = O.formula
    S.designs = D; S.n_slots = 1; S.max_entries = O.max_entries; S.retained_entries = 0
    cdid_ref_encoded(S, 1, generated, expanded, omitted)
    O.X = S.designs[1].value; O.xcols = S.designs[1].columns
    if (cols(O.X) < 1)
        _cdid_oca_fail("Custom ordinary estimation requires at least one covariate or intercept column.", 459)
    *p = O
}
string scalar cdid_ord_custom_fit_store(string scalar owner_name, string scalar limits_name,
    real scalar alpha, real scalar bstrap, real scalar cband, real scalar biters,
    real scalar compute_if, real scalar row_chunk, real scalar workspace)
{
    pointer() scalar p, plim, p_fn, p_pl
    struct cdid_ord_custom_owner scalar O
    struct cdid_state_limits scalar L
    struct cdid_infer_budget scalar budget
    struct cdid_ord_summary scalar S
    struct cdid_ord_boot scalar M
    struct cdid_callback_request scalar req
    struct cdid_co_response scalar cb_out
    struct cdid_ord_inference scalar inf
    struct cdid_result_state scalar state
    string matrix versions
    real matrix cells, F, X_cov
    real colvector is_C, valid_obs, D, y1, y0, y_v, post_v, ww, target_rows, fixed_w, na_w
    real colvector u_valid, inf_long
    real scalar nT, nG, tfac, K, g_idx, t_idx, pret_g, pret, g_val, t_val, base_val
    real scalar cell_idx, call_idx, w_idx, cov_idx, n1, role, n_t1, n_t0, n_c1, n_c0
    real scalar target_period, i, u_id, cap
    p = findexternal(owner_name)
    if (p == NULL) _cdid_oca_fail("Ordinary custom object disappeared.", 111)
    O = *p
    cdid_custom_validate_symbols(O.callback_name, O.payload_name)
    p_fn = findexternal(O.callback_name + "()")
    p_pl = (O.payload_name != "" ? findexternal(O.payload_name) : NULL)
    nT = rows(O.tlist); nG = rows(O.glist)
    tfac = (O.base_period == "universal" ? 0 : 1)
    if (nT == 2) cband = 0
    K = nG * (nT - tfac)
    if (O.n_units < 1 | K < 1 | biters < 1 | biters != floor(biters) | biters >= . |
        workspace < 1 | workspace >= . | row_chunk < 1 | row_chunk != floor(row_chunk))
        _cdid_oca_fail("Invalid ordinary custom inference dimensions or explicit execution controls.", 198)
    cap = floor(workspace / 8)
    L.payload_slots = cap; L.payload_string_bytes = floor(workspace)
    L.registry_slots = cap; L.registry_string_bytes = floor(workspace)
    if (O.n_units * K + rows(O.rows) * (cols(O.X) + 8) > cap)
        _cdid_oca_fail("Ordinary custom cell/IF matrices exceed the explicit workspace budget.", 908)
    budget = cdid_infer_budget_create(O.n_units, O.n_units, K, biters,
        min((O.n_units, row_chunk)), max((K, min((cap, row_chunk * K)))),
        min((biters, row_chunk)), max((O.n_units, min((cap, row_chunk * O.n_units)))),
        cap, floor(workspace))
    cells = J(K, 12, .)
    F = J(O.n_units, K, 0)
    cell_idx = 0; call_idx = 0
    for (g_idx = 1; g_idx <= nG; g_idx++) {
        g_val = O.glist[g_idx]
        pret_g = max(_cdid_oca_ix(O.tlist :< (g_val - O.anticipation)))
        for (t_idx = 1; t_idx <= nT - tfac; t_idx++) {
            cell_idx = cell_idx + 1
            pret = (O.base_period == "universal" | g_val <= O.tlist[t_idx + tfac] ? pret_g : t_idx)
            t_val = O.tlist[t_idx + tfac]
            base_val = O.tlist[pret]
            if (O.base_period == "universal" & base_val == t_val) {
                cells[cell_idx,] = (g_val, t_val, base_val, 1, 1, 0, ., t_val - g_val, 0, 0, 0, 0)
                continue
            }
            role = (t_val < g_val ? 2 : 3)
            if (O.eff_panel) {
                if (O.control == "nevertreated") is_C = (O.source_group :== 0)
                else is_C = (O.source_group :== 0) :| ((O.source_group :> (max((t_val, base_val)) + O.anticipation)) :& (O.source_group :!= g_val))
                valid_obs = _cdid_oca_ix((O.source_group :== g_val) :| is_C)
                n_t1 = sum(O.source_group[valid_obs] :== g_val); n_t0 = n_t1
                n_c1 = sum(O.source_group[valid_obs] :!= g_val); n_c0 = n_c1
                if (n_t1 == 0 | n_c1 == 0) {
                    cells[cell_idx,] = (g_val, t_val, base_val, role, 0, ., ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
                    F[, cell_idx] = J(O.n_units, 1, .)
                    continue
                }
                if (O.fix_weights == "") w_idx = min((pret, t_idx))
                else if (O.fix_weights == "base_period") w_idx = pret_g
                else w_idx = 1
                cov_idx = min((pret, t_idx))
                y1 = O.y[((t_idx + tfac - 1) * O.n_units) :+ valid_obs]
                y0 = O.y[((pret - 1) * O.n_units) :+ valid_obs]
                D = (O.source_group[valid_obs] :== g_val)
                ww = O.w[((w_idx - 1) * O.n_units) :+ valid_obs]
                X_cov = O.X[((cov_idx - 1) * O.n_units) :+ valid_obs,]
                call_idx = call_idx + 1
                req = cdid_callback_request_init()
                req.mode = "ordinary"; req.route = "panel"
                req.panel = 1; req.inffunc = compute_if; req.gradient_required = 0
                req.N = O.n_units; req.call_index = call_idx; req.raw_index = cell_idx; req.model_index = cell_idx
                req.g = g_val; req.s = base_val; req.t = t_val
                req.link_id = "ordinary_att_cell:" + _cdid_ord_key(cell_idx)
                req.rows_F = valid_obs; req.y0 = y0; req.y1 = y1
                req.D = D; req.i_weights = ww; req.covariates = X_cov; req.covariate_names = O.xcols
                cb_out = cdid_co_call_cell(p_fn, req, p_pl)
                if (!cb_out.available) {
                    errprintf("Warning: ordinary custom callback cell (g=%g, t=%g) is unavailable: %s\n", g_val, t_val, cb_out.reason)
                    cells[cell_idx,] = (g_val, t_val, base_val, role, 0, ., ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
                    F[, cell_idx] = J(O.n_units, 1, .)
                    continue
                }
                n1 = rows(valid_obs)
                if (compute_if) F[valid_obs, cell_idx] = (O.n_units / n1) :* cb_out.psi
                cells[cell_idx,] = (g_val, t_val, base_val, role, 1, cb_out.ATT, ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
            }
            else {
                if (O.control == "nevertreated") is_C = (O.g :== 0)
                else is_C = (O.g :== 0) :| ((O.g :> (max((t_val, base_val)) + O.anticipation)) :& (O.g :!= g_val))
                valid_obs = _cdid_oca_ix(((O.t :== t_val) :| (O.t :== base_val)) :& ((O.g :== g_val) :| is_C))
                if (O.eff_unbal & (O.fix_weights == "base_period" | O.fix_weights == "first_period")) {
                    target_period = (O.fix_weights == "base_period" ? O.tlist[pret_g] : O.tlist[1])
                    target_rows = _cdid_oca_ix(O.t :== target_period)
                    fixed_w = J(O.n_units, 1, .)
                    if (rows(target_rows)) fixed_w[O.row_unit[target_rows]] = O.w[target_rows]
                    na_w = (fixed_w[O.row_unit[valid_obs]] :>= .)
                    valid_obs = select(valid_obs, !na_w)
                    ww = fixed_w[O.row_unit[valid_obs]]
                }
                else ww = O.w[valid_obs]
                D = (O.g[valid_obs] :== g_val)
                post_v = (O.t[valid_obs] :== t_val)
                n_t1 = sum(D :* post_v); n_t0 = sum(D :* (1 :- post_v))
                n_c1 = sum((1 :- D) :* post_v); n_c0 = sum((1 :- D) :* (1 :- post_v))
                if (n_t1 == 0 | n_t0 == 0 | n_c1 == 0 | n_c0 == 0) {
                    cells[cell_idx,] = (g_val, t_val, base_val, role, 0, ., ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
                    F[, cell_idx] = J(O.n_units, 1, .)
                    continue
                }
                y_v = O.y[valid_obs]
                X_cov = O.X[valid_obs,]
                call_idx = call_idx + 1
                req = cdid_callback_request_init()
                req.mode = "ordinary"; req.route = (O.eff_unbal ? "unbalanced_panel" : "RCS")
                req.panel = 0; req.inffunc = compute_if; req.gradient_required = 0
                req.N = O.n_units; req.call_index = call_idx; req.raw_index = cell_idx; req.model_index = cell_idx
                req.g = g_val; req.s = base_val; req.t = t_val
                req.link_id = "ordinary_att_cell:" + _cdid_ord_key(cell_idx)
                req.rows_F = O.row_unit[valid_obs]; req.y = y_v; req.post = post_v
                req.D = D; req.i_weights = ww; req.covariates = X_cov; req.covariate_names = O.xcols
                cb_out = cdid_co_call_cell(p_fn, req, p_pl)
                if (!cb_out.available) {
                    errprintf("Warning: ordinary custom callback cell (g=%g, t=%g) is unavailable: %s\n", g_val, t_val, cb_out.reason)
                    cells[cell_idx,] = (g_val, t_val, base_val, role, 0, ., ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
                    F[, cell_idx] = J(O.n_units, 1, .)
                    continue
                }
                n1 = rows(valid_obs)
                if (compute_if) {
                    if (O.eff_unbal) {
                        inf_long = (O.n_units / n1) :* cb_out.psi
                        u_valid = O.row_unit[valid_obs]
                        for (i = 1; i <= rows(valid_obs); i++) {
                            u_id = u_valid[i]
                            F[u_id, cell_idx] = F[u_id, cell_idx] + inf_long[i]
                        }
                    }
                    else F[valid_obs, cell_idx] = (O.n_units / n1) :* cb_out.psi
                }
                cells[cell_idx,] = (g_val, t_val, base_val, role, 1, cb_out.ATT, ., t_val - g_val, n_t1, n_t0, n_c1, n_c0)
            }
        }
    }
    if (sum(cells[, 5] :& (cells[, 4] :!= 1)) == 0)
        _cdid_oca_fail("No valid (g,t) cells were estimated.", 459)
    M = cdid_ord_boot_none(K)
    if (compute_if) {
        S = cdid_ord_if_summary(F, O.cluster_positions, O.C, O.cluster_name != "", bstrap, row_chunk)
        cells[, 7] = S.source_se
        if (bstrap) M = cdid_ord_boot_stream(S, biters, alpha, 0, budget, "caller", J(0, 0, .))
    }
    inf = cdid_ord_report_cells(cells[, 6], cells[, 1], cells[, 2], (compute_if ? F : J(0, K, .)),
        O.cluster_positions, O.C, O.cluster_name != "", O.cluster_name != "",
        compute_if, bstrap, cband, alpha, row_chunk, M)
    versions = ("R_cdid", "1.3.0.0" \
        "did", "2.5.1" \
        "estimator", "user-supplied (method(custom)) ordinary estimator" \
        "rcs_dr_influence_function", "unaffected")
    /* Custom psi is user-returned: D1_exception is strictly 0. */
    state = cdid_ord_state_build(cells, (compute_if ? F : J(0, K, .)), O.identity_map,
        O.source_group, O.raw_unit_w, O.source_weight, O.norm_factor, O.tlist,
        O.cluster_labels, O.cluster_positions, O.C, O.cluster_name, 0, versions,
        inf, alpha, biters, L)
    if (O.cluster_name != "" & O.cluster_is_str) {
        state.cluster_kind = "string"
        state.cluster_string_id = O.scluster
        state.cluster_numeric_id = J(O.n_units, 1, .)
    }
    _cdid_ord_meta_matrix(state, "ord.source_sample_rows", O.rows)
    _cdid_ord_meta_text(state, "ord.panel_mode", (O.eff_panel ? "panel" : (O.eff_unbal ? "unbalanced_panel" : "rcs")))
    _cdid_ord_meta_text(state, "ord.compute_path", "native_custom_ordinary_att")
    _cdid_ord_meta_text(state, "ord.method", "custom")
    _cdid_ord_meta_text(state, "ord.callback", O.callback_name)
    _cdid_ord_meta_text(state, "ord.payload", O.payload_name)
    _cdid_ord_meta_number(state, "ord.n_callbacks", call_idx)
    _cdid_ord_meta_text(state, "ord.control_group", O.control)
    _cdid_ord_meta_text(state, "ord.base_period", O.base_period)
    _cdid_ord_meta_text(state, "ord.fix_weights", O.fix_weights)
    _cdid_ord_meta_number(state, "ord.anticipation", O.anticipation)
    _cdid_ord_meta_text(state, "ord.native_multiplier", "accepted native caller stream; increasing draw batches, column-major cluster order inside batch; distinct from R RNG")
    plim = crexternal(limits_name)
    if (plim == NULL) _cdid_oca_fail("Cannot create explicit state limits.", 498)
    *plim = L
    return(cdid_state_put(state, L))
}
/* Data rows retained by the custom ordinary intake (for e(sample)). */
real colvector cdid_ord_custom_rows(string scalar owner_name)
{
    pointer(struct cdid_ord_custom_owner scalar) scalar p
    p = findexternal(owner_name)
    if (p == NULL) _error(3499, "cdid_ord_custom_rows(): intake owner not found")
    return((*p).rows)
}
end
