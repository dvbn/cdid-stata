*! Public persistence/replay loader helpers for native cdid; current StataNow19.5.
*! Reuses cdid_result_state, cdid_result_bridge, cdid_ordinary_post,
*! cdid_chained_aggregate_public, cdid_frontend_display, and cdid_pretest_result.
*! Read-only restoration of saved state payloads: no data read, refit, or RNG draws.
version 19.5
mata:

real scalar _cdid_pload_has_text(struct cdid_result_state scalar P, string scalar key)
{
    if (rows(P.text_settings) == 0) return(0)
    return(rows(_cdid_rb_index(P.text_settings[, 1] :== key)) >= 1)
}

void _cdid_pload_sync_posting(string scalar token, string scalar posting, string scalar limits_name)
{
    struct cdid_state_registry scalar R
    struct cdid_result_state scalar P
    struct cdid_state_limits scalar L
    pointer() scalar p_lim, p_reg
    real scalar ix
    real colvector hit

    p_lim = findexternal(limits_name)
    if (p_lim == NULL) _cdid_rb_fail("Explicit state limits disappeared before posting.", 111)
    L = *p_lim
    R = _cdid_rs_registry()
    ix = _cdid_rs_index(R, token)
    P = *R.payloads[ix]
    if (rows(P.text_settings) == 0) return
    hit = _cdid_rb_index(P.text_settings[, 1] :== "cdid_bridge.posting_mode")
    if (rows(hit) == 1) {
        if (P.text_settings[hit[1], 2] != posting) {
            P.text_settings[hit[1], 2] = posting
            cdid_state_validate(P, L)
            *R.payloads[ix] = P
            R.costs[ix, ] = cdid_state_cost(P)
            p_reg = findexternal("CDID_PRIVATE_RESULT_REGISTRY")
            *p_reg = R
        }
    }
}

string scalar cdid_public_load_open(string scalar path, real scalar workspace, string scalar limits_name)
{
    struct cdid_state_limits scalar L
    struct cdid_result_state scalar P
    pointer() scalar p_lim
    string scalar tok, mode, kind

    if (strtrim(path) == "") {
        _cdid_rs_fail("Supply a saved result-state file path.", 198)
    }
    if (workspace < 8 | workspace >= . | workspace != floor(workspace) | workspace > 2^53 - 1) {
        _cdid_rs_fail("workspace() must be an explicit integer of at least 8 bytes.", 198)
    }
    L.payload_slots = max((1, floor(workspace / 8)))
    L.payload_string_bytes = floor(workspace)
    L.registry_slots = 2^53 - 1
    L.registry_string_bytes = 2^53 - 1

    tok = cdid_state_load(path, L)
    P = cdid_state_get(tok)
    mode = _cdid_ca_opt_text(P, "ord.mode", "")
    if (mode == "ordinary") {
        kind = _cdid_ca_opt_text(P, "ord.report.kind", "")
        if (kind != "ordinary_cells" & kind != "ordinary_simple" &
            kind != "ordinary_group" & kind != "ordinary_dynamic" &
            kind != "ordinary_calendar") {
            cdid_state_drop(tok)
            _cdid_rs_fail("Saved ordinary state lacks a valid ord.report.kind.", 498)
        }
        p_lim = crexternal(limits_name)
        if (p_lim == NULL) {
            cdid_state_drop(tok)
            _cdid_rs_fail("Cannot allocate temporary state limits.", 908)
        }
        *p_lim = L
        st_local("pload_mode", "ordinary")
        st_local("pload_kind", kind)
        st_local("pload_depvar", _cdid_ca_opt_text(P, "depvar", ""))
        return(tok)
    }

    kind = _cdid_ca_opt_text(P, "cdid_bridge.kind", "")
    if (kind != "cells" & kind != "aggregation") {
        cdid_state_drop(tok)
        _cdid_rs_fail("Saved chained state lacks a valid cdid_bridge.kind (cells or aggregation).", 498)
    }
    p_lim = crexternal(limits_name)
    if (p_lim == NULL) {
        cdid_state_drop(tok)
        _cdid_rs_fail("Cannot allocate temporary state limits.", 908)
    }
    *p_lim = L
    st_local("pload_mode", "chained")
    st_local("pload_kind", kind)
    st_local("pload_depvar", _cdid_ca_opt_text(P, "depvar", ""))
    return(tok)
}

void cdid_public_load_stage(
    string scalar token,
    string scalar limits_name,
    string scalar posting,
    string scalar bname,
    string scalar vname,
    string scalar oname,
    string scalar sename,
    string scalar aname,
    string scalar ciname,
    string scalar piname,
    string scalar mapname,
    real scalar matrix_limit)
{
    struct cdid_result_state scalar P
    real matrix b, V, Omega, se, analytic_se, ci, point_ci, coef_map
    real matrix pts, cell_av, se_vec, ase_vec, pci_mat, bci_mat, V_mat, Omega_mat, ov, ov_ase
    real colvector report_cells, report, av
    real scalar K, j, n_F, C, clustered, alpha, critical, normal_critical
    real scalar wald_df, wald_statistic, wald_p, wald_C, bstrap, cband, B, point_only
    real scalar ise, iase, icilo, icihi, iblo, ibhi
    real scalar has_rng_draws, draws_used, n_draws, ov_draws_used, ov_n_draws
    string scalar kind, status, inf_status, flag_law, warning, wald_status, wald_reason, rng_mode
    string matrix stripe
    string colvector retained_settings

    if (posting != "matrix" & posting != "state") {
        _cdid_rb_fail("Posting mode must be matrix or state.", 198)
    }
    if (matrix_limit < 1 | matrix_limit >= . | matrix_limit != floor(matrix_limit)) {
        _cdid_rb_fail("Current standard Stata matrix limit was not supplied.", 498)
    }

    P = cdid_state_get(token)
    kind = _cdid_rb_stored_text(P, "cdid_bridge.kind")
    if (kind != "cells" & kind != "aggregation") {
        _cdid_rb_fail("Saved chained state requires cdid_bridge.kind equal to cells or aggregation.", 498)
    }

    status = _cdid_rb_stored_text(P, "cdid_bridge.status")
    inf_status = _cdid_rb_stored_text(P, "cdid_bridge.inference_status")
    flag_law = _cdid_rb_stored_text(P, "cdid_bridge.flag_law")
    warning = _cdid_rb_stored_text(P, "cdid_bridge.warning")
    wald_status = _cdid_rb_stored_text(P, "cdid_bridge.wald_status")
    wald_reason = _cdid_rb_stored_text(P, "cdid_bridge.wald_reason")

    clustered = _cdid_rb_stored_number(P, "clustered")
    n_F = _cdid_rb_stored_number(P, "N")
    C = _cdid_rb_stored_number(P, "C")
    alpha = _cdid_rb_stored_number(P, "alpha")
    critical = _cdid_rb_stored_number(P, "critical")
    normal_critical = _cdid_rb_stored_number(P, "normal_critical")
    wald_df = _cdid_rb_stored_number(P, "wald_df")
    wald_statistic = _cdid_rb_stored_number(P, "wald_statistic")
    wald_p = _cdid_rb_stored_number(P, "wald_p")
    wald_C = _cdid_rb_stored_number(P, "wald_C")
    bstrap = _cdid_rb_stored_number(P, "bstrap")
    cband = _cdid_rb_stored_number(P, "cband")
    B = _cdid_rb_stored_number(P, "B")
    point_only = _cdid_rb_stored_number(P, "point_only")

    rng_mode = _cdid_ca_opt_text(P, "cdid_bridge.rng.mode", "")
    has_rng_draws = _cdid_ca_has_num(P, "cdid_bridge.rng.draws_used")
    draws_used = _cdid_ca_num_first(P, ("cdid_bridge.rng.draws_used"), 0)
    n_draws = _cdid_ca_num_first(P, ("cdid_bridge.rng.n_draws"), 0)

    b = se = analytic_se = J(1, 0, .)
    V = Omega = J(0, 0, .)
    point_ci = ci = J(2, 0, .)
    coef_map = J(0, 2, .)

    if (kind == "cells") {
        K = rows(P.point)
        ise = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.se")
        iase = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.analytic_se")
        icilo = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.ci_lower")
        icihi = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.ci_upper")
        iblo = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.band_lower")
        ibhi = _cdid_rb_column(P.ledger_columns, "cdid_bridge.cells.band_upper")
        report_cells = _cdid_ca_report_cells(P)
        report = (rows(report_cells) ? _cdid_rb_index(P.available[report_cells] :== 1) : J(0, 1, .))
        report = (rows(report) ? report_cells[report] : J(0, 1, .))
        if (rows(report)) {
            if (missing(P.point[report])) {
                _cdid_rb_fail("Available reporting points must be finite.", 498)
            }
        }
        if (posting == "matrix") {
            if (rows(report) > matrix_limit) {
                _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.", 908)
            }
            if (rows(report)) {
                b = P.point[report]'
                se = P.ledger[report, ise]'
                analytic_se = P.ledger[report, iase]'
                point_ci = P.ledger[report, (icilo, icihi)]'
                ci = P.ledger[report, (iblo, ibhi)]'
                coef_map = ((1::rows(report)), report)
                if (!point_only & rows(P.V_estimator) == K & cols(P.V_estimator) == K) {
                    V = P.V_estimator[report, report']
                    Omega = P.V_asymptotic[report, report']
                }
            }
        }
    }
    else {
        pts = _cdid_ca_read_matrix(P, "agg.points")
        cell_av = _cdid_ca_read_matrix(P, "agg.available")
        se_vec = _cdid_ca_read_matrix(P, "agg.se")
        ase_vec = _cdid_ca_read_matrix(P, "agg.analytic_se")
        pci_mat = _cdid_ca_read_matrix(P, "agg.point_ci")
        bci_mat = _cdid_ca_read_matrix(P, "agg.band_ci")
        V_mat = _cdid_ca_read_matrix(P, "agg.V_estimator")
        Omega_mat = _cdid_ca_read_matrix(P, "agg.V_asymptotic")
        ov = _cdid_ca_read_matrix(P, "agg.overall")
        ov_ase = _cdid_ca_read_matrix(P, "agg.overall_analytic_se")
        K = rows(pts)
        av = _cdid_rb_index(cell_av :== 1)
        if (posting == "matrix") {
            if (rows(av) > matrix_limit) {
                _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.", 908)
            }
            if (rows(av)) {
                b = pts[av]'
                se = se_vec[av]'
                analytic_se = ase_vec[av]'
                point_ci = pci_mat[av, ]'
                ci = bci_mat[av, ]'
                coef_map = ((1::rows(av)), av)
                if (!point_only & rows(V_mat) == K & cols(V_mat) == K) {
                    V = V_mat[av, av']
                    Omega = Omega_mat[av, av']
                }
            }
        }
        st_local("bridge_selector_warning", _cdid_ca_opt_text(P, "cdid_bridge.agg.selector_warning", ""))
        st_local("bridge_agg_type", _cdid_rb_stored_text(P, "cdid_bridge.agg.type"))
        st_local("bridge_target_def", _cdid_rb_stored_text(P, "cdid_bridge.agg.target_definition"))
        st_local("bridge_band_family", _cdid_rb_stored_text(P, "cdid_bridge.agg.band_family"))
        st_local("bridge_na_rm", _cdid_rb_key(_cdid_rb_stored_number(P, "agg.na_rm")))
        st_local("bridge_min_e", strofreal(_cdid_rb_stored_number(P, "agg.min_e"), "%24.17g"))
        st_local("bridge_max_e", strofreal(_cdid_rb_stored_number(P, "agg.max_e"), "%24.17g"))
        st_local("bridge_balance_e", strofreal(_cdid_rb_stored_number(P, "agg.balance_e"), "%24.17g"))
        st_local("bridge_overall_est", strofreal(ov[1, 1], "%24.17g"))
        st_local("bridge_overall_se", strofreal(ov[1, 2], "%24.17g"))
        st_local("bridge_overall_analytic_se", strofreal((rows(ov_ase) == 1 ? ov_ase[1, 1] : .), "%24.17g"))
        st_local("bridge_overall_ci_lo", strofreal(ov[1, 3], "%24.17g"))
        st_local("bridge_overall_ci_hi", strofreal(ov[1, 4], "%24.17g"))
        st_local("bridge_overall_status", _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.status"))
        st_local("bridge_overall_reason", _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.reason"))
        st_local("bridge_overall_inf_status", _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.inference_status"))
        st_local("bridge_overall_inf_kind", _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.inference_kind"))
        ov_draws_used = _cdid_ca_num_first(P, ("cdid_bridge.agg.overall.draws_used"), 0)
        ov_n_draws = _cdid_ca_num_first(P, ("cdid_bridge.agg.overall.n_draws"), 0)
        st_local("bridge_overall_draws_used", _cdid_rb_key(ov_draws_used))
        st_local("bridge_overall_n_draws", _cdid_rb_key(ov_n_draws))
    }

    _cdid_pload_sync_posting(token, posting, limits_name)

    stripe = J(cols(b), 2, "")
    for (j = 1; j <= cols(b); j++) {
        stripe[j, 2] = (kind == "cells" ? "cell_" : "agg_") + _cdid_rb_key(coef_map[j, 2])
    }
    if (cols(b)) {
        st_matrix(bname, b)
        st_matrixcolstripe(bname, stripe)
        st_matrix(sename, se)
        st_matrix(aname, analytic_se)
        st_matrix(ciname, ci)
        st_matrix(piname, point_ci)
        st_matrix(mapname, coef_map)
    }
    if (rows(V)) {
        st_matrix(vname, V)
        st_matrixcolstripe(vname, stripe)
        st_matrixrowstripe(vname, stripe)
        st_matrix(oname, Omega)
        st_matrixcolstripe(oname, stripe)
        st_matrixrowstripe(oname, stripe)
    }

    retained_settings = ("method" \ "gmm_weighting" \ "gmm_covariance" \ "base_period" \
        "control_group" \ "pt_assumption" \ "selection" \ "selection_canonical" \
        "selection_model" \ "selection_pooling" \ "selection_universe" \
        "selection_est_method" \ "imputation_pooling" \ "outcome_pooling" \ "selection_justification" \
        "inactive_arguments" \ "weight_source" \ "weight_reduction" \ "fix_weights_requested")
    for (j = 1; j <= rows(retained_settings); j++) {
        st_local("pl_has_" + retained_settings[j], _cdid_rb_key(_cdid_pload_has_text(P, retained_settings[j])))
        st_local("bridge_" + retained_settings[j], _cdid_ca_opt_text(P, retained_settings[j], ""))
    }
    st_local("bridge_has_anticipation", _cdid_rb_key(_cdid_ca_has_num(P, "anticipation")))
    if (_cdid_ca_has_num(P, "anticipation")) {
        st_local("bridge_anticipation", strofreal(_cdid_ca_num_first(P, ("anticipation"), 0), "%24.17g"))
    }
    st_local("bridge_has_rng", _cdid_rb_key(has_rng_draws))
    st_local("bridge_rng_mode", rng_mode)
    st_local("bridge_draws_used", _cdid_rb_key(draws_used))
    st_local("bridge_n_draws", _cdid_rb_key(n_draws))

    st_local("bridge_npost", _cdid_rb_key(cols(b)))
    st_local("bridge_hasV", _cdid_rb_key(rows(V) > 0))
    st_local("bridge_N_source", _cdid_rb_key(rows(P.point)))
    st_local("bridge_alpha", strofreal(alpha, "%24.17g"))
    st_local("bridge_B", strofreal(B, "%24.17g"))
    st_local("bridge_bstrap", strofreal(bstrap, "%24.17g"))
    st_local("bridge_cband", _cdid_rb_key(cband))
    st_local("bridge_clustered", _cdid_rb_key(clustered))
    st_local("bridge_N", _cdid_rb_key(n_F))
    st_local("bridge_C", _cdid_rb_key(C))
    st_local("bridge_K", _cdid_rb_key(K))
    st_local("bridge_wald_C", strofreal(wald_C, "%24.17g"))
    st_local("bridge_posting_mode", posting)
    st_local("bridge_kind", kind)
    st_local("bridge_status", status)
    st_local("bridge_inference_status", inf_status)
    st_local("bridge_flags", flag_law)
    st_local("bridge_warning", warning)
    st_local("bridge_wald_status", wald_status)
    st_local("bridge_wald_reason", wald_reason)
    st_local("bridge_critical", strofreal(critical, "%24.17g"))
    st_local("bridge_normal", strofreal(normal_critical, "%24.17g"))
    st_local("bridge_wald_df", strofreal(wald_df, "%24.17g"))
    st_local("bridge_wald_stat", strofreal(wald_statistic, "%24.17g"))
    st_local("bridge_wald_p", strofreal(wald_p, "%24.17g"))
    st_local("bridge_point_only", _cdid_rb_key(point_only))
}

void cdid_public_load_display(string scalar token)
{
    struct cdid_result_state scalar P
    P = cdid_state_get(token)
    if (_cdid_ca_opt_text(P, "cdid_bridge.kind", "") == "aggregation") cdid_chained_agg_replay(token)
    else cdid_frontend_display(token)
}

string scalar _cdid_pt_req_val(string matrix req, string scalar key, real scalar col, string scalar fallback)
{
    real colvector ix
    if (rows(req) == 0 | cols(req) < col) return(fallback)
    ix = _cdid_rb_index(req[, 1] :== key)
    if (rows(ix) >= 1) return(req[ix[rows(ix)], col])
    return(fallback)
}

real scalar _cdid_pt_req_has(string matrix req, string scalar key)
{
    if (rows(req) == 0 | cols(req) < 1) return(0)
    return(rows(_cdid_rb_index(req[, 1] :== key)) >= 1)
}

void cdid_pt_load_meta(string scalar name)
{
    struct cdid_pt_result scalar R
    pointer() scalar p
    string scalar panel, covars, noconstant, ivar, raw_allow, res_allow, outpool
    real scalar bs, cb, d1_exc, rf04_form, rf04_id, allow_sup

    p = findexternal(name)
    if (p == NULL) _cdid_pt_fail("Loaded conditional result object disappeared.", 111)
    R = *p

    bs = R.bstrap_requested
    cb = R.cband_requested
    panel = _cdid_pt_req_val(R.request, "panel", 3, "")
    covars = _cdid_pt_req_val(R.request, "covars", 3, R.formula)
    noconstant = _cdid_pt_req_val(R.request, "noconstant", 3, "")
    ivar = _cdid_pt_req_val(R.request, "ivar", 3, "")
    raw_allow = _cdid_pt_req_val(R.request, "allowunbalanced", 2, "")
    res_allow = _cdid_pt_req_val(R.request, "allowunbalanced", 3, "")
    outpool = _cdid_pt_req_val(R.request, "outcomepooling", 3, "")

    d1_exc = (R.mode != "ordinary" ? 0 : (panel == "" ? . : (panel == "0" & R.method == "dr")))
    if (_cdid_pt_req_has(R.request, "covars") | _cdid_pt_req_has(R.request, "noconstant")) {
        rf04_form = (R.mode == "ordinary" & covars == "" & noconstant == "")
    }
    else {
        rf04_form = (R.mode == "ordinary" ? . : 0)
    }
    if (_cdid_pt_req_has(R.request, "ivar")) {
        rf04_id = (R.mode != "ordinary" ? 0 : (panel == "" ? . : (panel == "0" & ivar == "")))
    }
    else {
        rf04_id = (R.mode == "ordinary" ? . : 0)
    }
    allow_sup = (_cdid_pt_req_has(R.request, "allowunbalanced") ? (raw_allow != "") : .)

    st_local("pt_bstrap_req", strofreal(bs, "%21.0f"))
    st_local("pt_cband_req", strofreal(cb, "%21.0f"))
    st_local("pt_d1_exc", strofreal(d1_exc, "%21.0f"))
    st_local("pt_rf04_form", strofreal(rf04_form, "%21.0f"))
    st_local("pt_rf04_id", strofreal(rf04_id, "%21.0f"))
    st_local("pt_allow_sup", strofreal(allow_sup, "%21.0f"))
    st_local("pt_allow_req", res_allow)
    st_local("pt_outcome_pooling", outpool)
    st_local("pt_panel_mode", (panel == "0" ? "rcs" : "panel"))
}

end
