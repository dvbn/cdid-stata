*! Private bounded public chained aggregation composition, current StataNow19.5.
*! Consumes cdid_result_state, cdid_inference, cdid_public_rng,
*! cdid_aggregation, and cdid_result_bridge modules without refitting or
*! retaining full BxC sign matrices in caller-stream mode.
version 19.5
mata:

// Persist only interface names on the just-created private fit staging token.
// The replacement token is registered before releasing that owned staging token.
string scalar cdid_chained_bind_names(string scalar token, string scalar request,
    string scalar id_name, string scalar cluster_name)
{
    struct cdid_result_state scalar P
    struct cdid_bridge_request scalar R
    pointer() scalar ptr
    string scalar next_token
    ptr = findexternal(request)
    if (ptr == NULL) _cdid_rb_fail("Fit inference request not found.", 111)
    R = *ptr
    P = cdid_state_get(token)
    P.text_settings = P.text_settings \ ("id_name", id_name) \
        ("cluster_name", (P.cluster_kind == "none" ? "" : cluster_name))
    next_token = cdid_state_put(P, R.limits)
    cdid_state_drop(token)
    return(next_token)
}

real scalar _cdid_ca_num_first(struct cdid_result_state scalar P,
    string colvector keys, real scalar fallback)
{
    real scalar k
    real colvector ix
    for (k = 1; k <= rows(keys); k++) {
        ix = _cdid_rb_index(P.numeric_setting_names :== keys[k])
        if (rows(ix) >= 1) {
            return(P.numeric_settings[ix[rows(ix)]])
        }
    }
    return(fallback)
}

string scalar _cdid_ca_opt_text(struct cdid_result_state scalar P,
    string scalar key, string scalar fallback)
{
    real colvector ix
    if (rows(P.text_settings) == 0) return(fallback)
    ix = _cdid_rb_index(P.text_settings[, 1] :== key)
    if (rows(ix) >= 1) {
        return(P.text_settings[ix[rows(ix)], 2])
    }
    return(fallback)
}

real scalar _cdid_ca_has_num(struct cdid_result_state scalar P, string scalar key)
{
    return(rows(_cdid_rb_index(P.numeric_setting_names :== key)) >= 1)
}

real matrix _cdid_ca_read_matrix(struct cdid_result_state scalar P, string scalar key)
{
    real scalar r, c, i, j
    real matrix X
    r = _cdid_rb_stored_number(P, key + ".rows")
    c = _cdid_rb_stored_number(P, key + ".cols")
    X = J(r, c, .)
    for (j = 1; j <= c; j++) {
        for (i = 1; i <= r; i++) {
            X[i, j] = _cdid_rb_stored_number(P, key + "." + _cdid_rb_key(i) + "." + _cdid_rb_key(j))
        }
    }
    return(X)
}

real colvector _cdid_ca_unit_cluster(struct cdid_result_state scalar P)
{
    real scalar N, i, g, found
    real colvector ix, pos, seen_num
    string colvector seen_str
    N = rows(P.cohort)
    if (cols(P.unit_columns) > 0) {
        ix = _cdid_rb_index(P.unit_columns :== "cdid_bridge.cells.cluster_position.1")
        if (rows(ix) == 1) return(P.unit_metadata[, ix[1]])
        ix = _cdid_rb_index(P.unit_columns :== "cdid_bridge.cluster_position.1")
        if (rows(ix) == 1) return(P.unit_metadata[, ix[1]])
        ix = _cdid_rb_index(P.unit_columns :== "cluster_position")
        if (rows(ix) == 1) return(P.unit_metadata[, ix[1]])
    }
    if (P.cluster_kind == "none") {
        return(1::N)
    }
    pos = J(N, 1, .)
    if (P.cluster_kind == "numeric") {
        seen_num = J(0, 1, .)
        for (i = 1; i <= N; i++) {
            found = 0
            for (g = 1; g <= rows(seen_num); g++) {
                if (seen_num[g] == P.cluster_numeric_id[i]) {
                    pos[i] = g
                    found = 1
                    break
                }
            }
            if (!found) {
                seen_num = seen_num \ P.cluster_numeric_id[i]
                pos[i] = rows(seen_num)
            }
        }
        return(pos)
    }
    seen_str = J(0, 1, "")
    for (i = 1; i <= N; i++) {
        found = 0
        for (g = 1; g <= rows(seen_str); g++) {
            if (seen_str[g] == P.cluster_string_id[i]) {
                pos[i] = g
                found = 1
                break
            }
        }
        if (!found) {
            seen_str = seen_str \ P.cluster_string_id[i]
            pos[i] = rows(seen_str)
        }
    }
    return(pos)
}

real colvector _cdid_ca_report_cells(struct cdid_result_state scalar P)
{
    real scalar K, j
    real colvector out
    string scalar role
    if (_cdid_ca_has_num(P, "cdid_bridge.cells.report_cells.rows")) {
        out = vec(_cdid_ca_read_matrix(P, "cells.report_cells"))
        return(rows(out) ? out : J(0, 1, .))
    }
    if (_cdid_ca_has_num(P, "cdid_bridge.report_cells.rows")) {
        out = vec(_cdid_ca_read_matrix(P, "report_cells"))
        return(rows(out) ? out : J(0, 1, .))
    }
    K = rows(P.point)
    out = J(0, 1, .)
    for (j = 1; j <= K; j++) {
        role = _cdid_ca_opt_text(P, "cell_role:" + _cdid_rb_key(j), "")
        if (role != "auxiliary") {
            out = out \ j
        }
    }
    return(out)
}

real colvector _cdid_ca_pre_cells(struct cdid_result_state scalar P)
{
    real scalar K, j
    real colvector out
    string scalar role
    if (_cdid_ca_has_num(P, "cdid_bridge.pre_cells.rows")) {
        out = vec(_cdid_ca_read_matrix(P, "pre_cells"))
        return(rows(out) ? out : J(0, 1, .))
    }
    K = rows(P.point)
    out = J(0, 1, .)
    for (j = 1; j <= K; j++) {
        role = _cdid_ca_opt_text(P, "cell_role:" + _cdid_rb_key(j), "")
        if (role == "pre") {
            out = out \ j
        }
    }
    return(out)
}

real colvector _cdid_ca_cohort_grid(struct cdid_result_state scalar P)
{
    real scalar ig
    real colvector raw, sub, cg
    raw = J(0, 1, .)
    if (_cdid_ca_has_num(P, "cdid_bridge.agg.request_cohort_grid.rows")) {
        cg = vec(_cdid_ca_read_matrix(P, "agg.request_cohort_grid"))
        if (rows(cg) > 0) {
            sub = select(cg, cg :< .)
            if (rows(sub) > 0) raw = raw \ sub
        }
    }
    if (rows(P.shares) > 0 & cols(P.shares) >= 1) {
        sub = select(P.shares[, 1], P.shares[, 1] :< .)
        if (rows(sub) > 0) raw = raw \ sub
    }
    if (rows(P.ledger) > 0 & cols(P.ledger_columns) > 0) {
        if (rows(_cdid_rb_index(P.ledger_columns :== "g")) == 1) {
            ig = _cdid_rb_column(P.ledger_columns, "g")
            sub = select(P.ledger[, ig], P.ledger[, ig] :< .)
            if (rows(sub) > 0) raw = raw \ sub
        }
    }
    if (rows(raw) == 0) return(J(0, 1, .))
    return(uniqrows(sort(raw, 1)))
}

void cdid_chained_agg_requests(
    string scalar token,
    string scalar type,
    real scalar min_e,
    real scalar max_e,
    real scalar balance_e,
    real scalar na_rm,
    string scalar bstrap_str,
    real scalar cband_req,
    string scalar biters_str,
    real scalar level,
    real scalar workspace,
    real scalar chunk,
    string scalar posting,
    string scalar requested_cluster,
    string scalar req_use_name,
    string scalar agg_use_name)
{
    struct cdid_result_state scalar P
    struct cdid_bridge_request scalar R
    struct cdid_agg_request scalar Q
    struct cdid_infer_flags scalar bin_flags
    struct cdid_agg_ledger scalar ledger
    real colvector post
    real scalar joint_targets
    real scalar n_F, src_clustered, clustered, C, bstrap, cband, alpha, B, B_inh, needs_draws
    real scalar bytes, entries, J_bound, B_bound, if_rows, if_cells, sign_rows, sign_cells
    real colvector unit_cluster
    string scalar endpoint, rc_trim, rc_low, stored_cl
    pointer() scalar p_req, p_agg

    if (strtrim(token) == "") {
        _cdid_rb_fail("Chained aggregation requires a retained state token in e(state_token).", 301)
    }
    P = cdid_state_get(token)
    if (_cdid_ca_opt_text(P, "ord.mode", "") == "ordinary") {
        _cdid_rb_fail("Chained aggregation requires a chained state payload.", 498)
    }
    if (_cdid_ca_num_first(P, ("cdid_bridge.point_only" \ "point_only"), 0) == 1) {
        _cdid_rb_fail("estat aggregations need influence functions; refit cdid without noinffunc.", 459)
    }
    if (_cdid_ca_num_first(P, ("compute_inffunc"), 1) == 0 | rows(P.IF) == 0) {
        _cdid_rb_fail("Point-only state cannot be aggregated: original influence functions are required.", 459)
    }
    if (type != "simple" & type != "group" & type != "dynamic" & type != "calendar") {
        _cdid_rb_fail("Aggregation type must be simple, group, dynamic, or calendar.", 498)
    }
    if (na_rm != 0 & na_rm != 1) {
        _cdid_rb_fail("na_rm must be 0 or 1.", 498)
    }
    if (posting != "matrix" & posting != "state") {
        _cdid_rb_fail("Posting mode must be matrix or state.", 498)
    }
    if (workspace < 8 | workspace >= . | workspace != floor(workspace) | workspace > 2^53-1 |
        chunk < 1 | chunk >= . | chunk != floor(chunk)) {
        _cdid_rb_fail("Workspace and chunk controls must be positive finite values.", 498)
    }

    if (bstrap_str == "-1" | bstrap_str == "" | bstrap_str == ".") {
        bstrap = _cdid_ca_num_first(P, ("cdid_bridge.cells.bstrap" \ "bstrap" \ "cdid_bridge.bstrap"), 0)
    }
    else {
        bstrap = strtoreal(bstrap_str)
        if (bstrap != 0 & bstrap != 1) {
            _cdid_rb_fail("bstrap() must be 0 or 1; omission or . inherits the fit setting.", 498)
        }
    }

    if (cband_req == -1) {
        cband = _cdid_ca_num_first(P, ("cdid_bridge.cells.cband" \ "cband" \ "cdid_bridge.cband"), 0)
    }
    else if (cband_req == 0 | cband_req == 1) {
        cband = cband_req
    }
    else {
        _cdid_rb_fail("cband() must be 0, 1, or omitted to inherit.", 498)
    }

    if (level >= .) {
        alpha = _cdid_ca_num_first(P, ("cdid_bridge.cells.alpha" \ "alpha" \ "alp" \ "cdid_bridge.alpha"), 0.05)
    }
    else {
        alpha = (level == 95 ? 0.05 : 1 - level / 100)
    }
    if (alpha <= 0 | alpha >= 1 | alpha >= .) {
        _cdid_rb_fail("Require explicit valid alpha in (0, 1).", 498)
    }

    endpoint = (type == "simple" ? "simple" : "aggte_bins")
    bin_flags = cdid_infer_resolve_flags(bstrap, cband, endpoint)
    needs_draws = (bin_flags.draws | bstrap == 1)

    if (biters_str == "-1" | biters_str == "") {
        B = _cdid_ca_num_first(P, ("biters" \ "cdid_bridge.cells.B" \ "cdid_bridge.B" \ "B"), 1000)
    }
    else if (biters_str == ".") {
        B = _cdid_ca_num_first(P, ("biters" \ "cdid_bridge.cells.B" \ "cdid_bridge.B" \ "B"), 1000)
    }
    else {
        B = strtoreal(biters_str)
    }
    if (B < 1 | B >= . | B != floor(B) | B > 2^31-1) {
        _cdid_rb_fail("biters() must be a positive whole number or omitted to inherit.", 498)
    }
    bin_flags = cdid_infer_resolve_flags(bstrap, cband, endpoint)

    n_F = rows(P.cohort)
    src_clustered = _cdid_ca_num_first(P, ("cdid_bridge.cells.clustered" \ "cdid_bridge.clustered"), (P.cluster_kind != "none"))
    if (src_clustered != 0) {
        unit_cluster = _cdid_ca_unit_cluster(P)
        C = _cdid_ca_num_first(P, ("cdid_bridge.cells.C" \ "cdid_bridge.C" \ "C"), max(unit_cluster))
    }
    else {
        unit_cluster = (1::n_F)
        C = n_F
    }

    rc_trim = strtrim(requested_cluster)
    rc_low = strlower(rc_trim)
    if (rc_trim == "" | rc_low == "inherit") {
        clustered = (src_clustered != 0)
    }
    else if (rc_low == "none" | rc_low == "off" | rc_low == "0" |
        rc_trim == _cdid_ca_opt_text(P, "id_name", "")) {
        if (src_clustered != 0) {
            stored_cl = _cdid_ca_opt_text(P, "cluster_name", _cdid_ca_opt_text(P, "cdid_bridge.cluster_name", ""))
            _cdid_rb_fail("cluster(" + rc_trim + ") is not available: the fit was clustered on " + stored_cl +
                ", and estat keeps the clustering of a chained fit, as R's aggte(); omit cluster() or specify cluster(" + stored_cl + ")", 498)
        }
        clustered = 0
    }
    else {
        if (src_clustered == 0) {
            _cdid_rb_fail("cluster(" + rc_trim + ") is not available: the fit was not clustered, and estat cannot add clustering after estimation; refit cdid with cluster(" + rc_trim + ")", 498)
        }
        stored_cl = _cdid_ca_opt_text(P, "cluster_name", _cdid_ca_opt_text(P, "cdid_bridge.cluster_name", ""))
        if (stored_cl == "" | rc_trim != stored_cl) {
            _cdid_rb_fail("cluster(" + rc_trim + ") is not available: the fit was clustered on " + stored_cl +
                ", and estat can use only that variable; refit cdid with cluster(" + rc_trim + ") to cluster on it", 498)
        }
        clustered = 1
    }

    bytes = floor(workspace)
    entries = max((1, floor(workspace / 8)))
    R.limits.payload_slots = entries
    R.limits.payload_string_bytes = bytes
    R.limits.registry_slots = 2^53-1
    R.limits.registry_string_bytes = 2^53-1

    // aggregation l.590 also summarizes the JOINT supported source
    // IF and cohort-share IF columns for the variance decomposition. That
    // dimension can exceed the number of source cells (5 + 2 = 7 in rb_payload).
    ledger = _cdid_rb_agg_ledger(P)
    post = _cdid_rb_index(ledger.time :>= ledger.g)
    joint_targets = sum(P.available :== 1)
    if (rows(post)) joint_targets = joint_targets + rows(uniqrows(sort(ledger.g[post], 1)))
    J_bound = max((1, rows(P.point), rows(P.calendar), rows(P.shares), joint_targets))
    B_bound = (B < . & B >= 1 ? B : 1)
    if_rows = min((n_F, max((1, floor(chunk)))))
    if_cells = max((J_bound, min((entries, if_rows * J_bound))))
    sign_rows = min((B_bound, max((1, floor(chunk)))))
    sign_cells = max((C, min((entries, sign_rows * C))))
    R.budget = cdid_infer_budget_create(n_F, C, J_bound, B_bound, if_rows, if_cells, sign_rows, sign_cells, entries, bytes)

    R.alpha = alpha
    R.compute_inffunc = 1
    R.bstrap = bstrap
    R.cband = cband
    R.B = B
    R.C = C
    R.clustered = clustered
    R.run_wald = 0
    R.unit_cluster = unit_cluster
    R.report_cells = _cdid_ca_report_cells(P)
    R.pre_cells = _cdid_ca_pre_cells(P)
    R.posting_mode = posting

    Q.type = type
    Q.na_rm = na_rm
    Q.balance_e = balance_e
    Q.min_e = min_e
    Q.max_e = max_e
    Q.alpha = alpha
    Q.bstrap = bstrap
    Q.cband = cband
    Q.B = B
    Q.calendar = P.calendar
    Q.cohort_grid = _cdid_ca_cohort_grid(P)

    p_req = crexternal(req_use_name)
    if (p_req == NULL) _cdid_rb_fail("Cannot allocate temporary bridge request.", 908)
    *p_req = R
    p_agg = crexternal(agg_use_name)
    if (p_agg == NULL) _cdid_rb_fail("Cannot allocate temporary aggregation request.", 908)
    *p_agg = Q
}

struct cdid_agg_inference scalar _cdid_pb_agg_apply_inference(
    struct cdid_agg_inference scalar base,
    real colvector estimates,
    real colvector cell_available,
    struct cdid_infer_flags scalar flags,
    real scalar B,
    real scalar alpha,
    struct cdid_infer_budget scalar budget,
    string scalar sign_mode,
    real matrix fixed_signs,
    string scalar err_prefix)
{
    struct cdid_agg_inference scalar out
    struct cdid_multiplier_state scalar ms
    struct cdid_multiplier_result scalar mult
    struct cdid_infer_result scalar rep
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar batch
    real colvector idx, positive
    real scalar entries, cap

    out = base
    if (out.n_available == 0) {
        return(out)
    }
    if (out.inference_status == "unavailable") {
        _cdid_rb_fail(err_prefix + out.reason, 459)
    }

    idx = _cdid_rb_index(cell_available :!= 0)
    mult = cdid_multiplier_none()
    positive = _cdid_rb_index(out.summary.analytic_se :> 0)
    if (flags.draws & rows(positive)) {
        entries = _cdid_infer_preflight(budget, out.summary.n_F, out.summary.C, out.summary.J, B)
        cap = min((budget.max_workspace_cells, floor(budget.max_workspace_bytes / 8)))
        if (sign_mode == "fixed") {
            entries = _cdid_infer_cells(rows(fixed_signs), cols(fixed_signs), cap, entries)
        }
        ms = cdid_multiplier_start(out.summary, positive, B, alpha, budget)
        stream = cdid_public_rng_start(B, out.summary.C, budget, sign_mode, fixed_signs)
        while (stream.cursor < stream.B) {
            batch = cdid_public_rng_next(stream)
            cdid_multiplier_add(ms, batch.draw_positions, batch.cluster_positions, batch.signs)
        }
        cdid_public_rng_complete(stream)
        mult = cdid_multiplier_finish(ms, 0)
    }
    rep = cdid_infer_report(estimates[idx], out.summary, flags, mult, alpha)
    if (rep.status != "computed") {
        _cdid_rb_fail(err_prefix + rep.reason, 459)
    }

    out.compact = rep
    out.inference_status = rep.status
    out.reason = rep.reason
    out.flag_law = rep.flag_law
    out.band_family = rep.band_family
    out.bootstrap = rep.bootstrap
    out.simultaneous = rep.simultaneous
    out.se[idx] = rep.se
    out.analytic_se[idx] = rep.analytic_se
    out.ci_lower[idx] = rep.ci_lower
    out.ci_upper[idx] = rep.ci_upper
    out.band_lower[idx] = rep.band_lower
    out.band_upper[idx] = rep.band_upper
    out.z_stat[idx] = rep.z_stat
    out.p_value[idx] = rep.p_value
    out.covariance[idx, idx'] = rep.covariance
    out.Omega[idx, idx'] = rep.Omega
    out.band_critical = rep.band_critical
    out.normal_critical = rep.normal_critical
    out.draws_used = rep.draws_used
    out.n_draws = rep.n_draws
    out.warning = rep.warning
    return(out)
}

struct cdid_bridge_result scalar _cdid_pb_agg_core(
    string scalar token,
    struct cdid_agg_request scalar Q,
    struct cdid_bridge_request scalar R,
    string scalar sign_mode,
    real matrix fixed_signs,
    real scalar matrix_limit)
{
    struct cdid_bridge_result scalar O
    struct cdid_result_state scalar P
    struct cdid_agg_ledger scalar L
    struct cdid_agg_frame scalar F
    struct cdid_agg_request scalar Q_base
    struct cdid_agg_result scalar A
    struct cdid_infer_flags scalar bin_flags, overall_flags
    real colvector av, point_flag
    real scalar j, reuse, overall_bstrap, draws_used
    string scalar endpoint

    if (sign_mode != "caller" & sign_mode != "fixed") {
        _cdid_rb_fail("Public sign source must be caller or private fixed signs.", 498)
    }
    P = cdid_state_get(token)
    _cdid_rb_request_check(P, R)
    point_flag = _cdid_rb_index(P.numeric_setting_names :== "cdid_bridge.point_only")
    if (rows(point_flag)) {
        if (_cdid_rb_stored_number(P, "point_only") == 1) {
            _cdid_rb_fail("estat aggregations need influence functions; refit cdid without noinffunc.", 459)
        }
    }
    if (!R.compute_inffunc | rows(P.IF) == 0) {
        _cdid_rb_fail("Point-only state cannot be aggregated: original influence functions are required.", 459)
    }
    if (R.run_wald) {
        _cdid_rb_fail("Aggregation preserves stored full-family pretests; request a new Wald through cell reporting.", 498)
    }
    if (R.alpha != Q.alpha | R.bstrap != Q.bstrap | R.cband != Q.cband | R.B != Q.B) {
        _cdid_rb_fail("Aggregation request and passed inference settings disagree.", 498)
    }
    if (rows(Q.calendar) != rows(P.calendar) | any(Q.calendar :!= P.calendar)) {
        _cdid_rb_fail("Aggregation calendar must be the stored original calendar.", 498)
    }

    endpoint = (Q.type == "simple" ? "simple" : "aggte_bins")
    bin_flags = cdid_infer_resolve_flags(Q.bstrap, Q.cband, endpoint)
    overall_bstrap = (Q.bstrap >= . ? 0 : Q.bstrap)
    overall_flags = cdid_infer_resolve_flags(overall_bstrap, 0, "cells")

    L = _cdid_rb_agg_ledger(P)
    F.n = rows(P.cohort)
    F.C = R.C
    F.clustered = R.clustered
    F.cohort = P.cohort
    F.cohort_known = P.cohort_known
    F.unit_cluster = R.unit_cluster
    F.frame_id = J(F.n, 1, "")
    for (j = 1; j <= F.n; j++) F.frame_id[j] = "F_row_" + _cdid_rb_key(j)
    if (P.id_kind == "string") F.frame_id = P.unit_string_id
    F.share_law = _cdid_rb_stored_text(P, "aggregation_share_law")
    F.w = P.weight
    if (F.share_law == "first_calendar_masses") {
        F.w = P.unit_metadata[, _cdid_rb_column(P.unit_columns, "first_calendar_weight")]
    }
    if (any(P.unit_columns :== "reference_share_excluded")) F.share_excluded = cdid_mrcp_state_excluded(P)

    // Compute exact aggregate point estimates, influence functions, cohort shares,
    // fixed-weight inference, variance decomposition, and IF summaries analytically
    // before any RNG draws.
    Q_base = Q
    Q_base.bstrap = 0
    Q_base.cband = 0
    Q_base.B = .
    A = cdid_aggregate(P.point, P.IF, P.available, L, F, Q_base, J(0, 0, .), R.budget)

    if (matrix_limit < . & R.posting_mode == "matrix") {
        if (sum(A.cell_available :== 1) > matrix_limit) {
            _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.", 908)
        }
    }

    A.inference = _cdid_pb_agg_apply_inference(
        A.inference,
        A.estimates,
        A.cell_available,
        bin_flags,
        Q.B,
        Q.alpha,
        R.budget,
        sign_mode,
        fixed_signs,
        "Requested aggregate inference refused: ")

    reuse = (Q.type == "simple") & !A.inference.simultaneous
    if (reuse) {
        A.overall.inference = A.inference
    }
    else {
        A.overall.inference.band_family = Q.type + "_overall"
        A.overall.inference = _cdid_pb_agg_apply_inference(
            A.overall.inference,
            J(1, 1, A.overall.estimate),
            J(1, 1, (A.overall.status == "available")),
            overall_flags,
            (overall_bstrap == 1 ? Q.B : .),
            Q.alpha,
            R.budget,
            sign_mode,
            fixed_signs,
            "Requested overall inference refused: ")
    }
    A.overall.inference_kind = (A.overall.inference.bootstrap ? "multiplier_bootstrap_iqr_se" : "analytic")
    A.overall.draws_reused_from = (reuse & A.inference.draws_used ? "inference" : "")
    if (A.overall.status == "available" & A.overall.inference.inference_status == "computed" & rows(A.overall.inference.se) == 1) {
        A.overall.se = A.overall.inference.se[1]
        A.overall.ci_lower = A.overall.inference.ci_lower[1]
        A.overall.ci_upper = A.overall.inference.ci_upper[1]
    }
    else {
        A.overall.se = .
        A.overall.ci_lower = .
        A.overall.ci_upper = .
    }

    O.payload = P
    _cdid_rb_clear_metadata(O.payload, 1)
    O.kind = "aggregation"
    O.n_F = F.n
    O.C = F.C
    O.K = A.M
    O.point_only = 0
    O.status = A.status
    O.inference_status = A.inference.inference_status
    O.flag_law = A.inference.flag_law
    O.warning = A.inference.warning
    O.normal_critical = A.inference.normal_critical
    O.critical = A.inference.band_critical
    _cdid_rb_preserve_wald(O, P)
    av = _cdid_rb_index(A.cell_available :== 1)
    if (rows(av)) {
        O.b = A.estimates[av]'
        O.se = A.inference.se[av]'
        O.analytic_se = A.inference.analytic_se[av]'
        O.V = A.inference.covariance[av, av']
        O.Omega = A.inference.Omega[av, av']
        O.point_ci = (A.inference.ci_lower[av], A.inference.ci_upper[av])'
        O.ci = (A.inference.band_lower[av], A.inference.band_upper[av])'
        O.coef_map = ((1::rows(av)), av)
    }
    else {
        // Mata's empty single subscript of a scalar is not an empty vector.
        // No supported bins means no posted coefficient or fabricated value.
        O.b = O.se = O.analytic_se = J(1, 0, .)
        O.V = O.Omega = J(0, 0, .)
        O.point_ci = O.ci = J(2, 0, .)
        O.coef_map = J(0, 2, .)
    }

    _cdid_rb_text(O.payload, "agg.selector_warning", A.selectors.warning)
    _cdid_rb_text(O.payload, "agg.type", A.type)
    _cdid_rb_text(O.payload, "agg.target_definition", A.target_definition)
    _cdid_rb_text(O.payload, "agg.overall.status", A.overall.status)
    _cdid_rb_text(O.payload, "agg.overall.reason", A.overall.reason)
    _cdid_rb_text(O.payload, "agg.overall.inference_status", A.overall.inference.inference_status)
    _cdid_rb_text(O.payload, "agg.overall.inference_kind", A.overall.inference_kind)
    _cdid_rb_text(O.payload, "agg.overall.draws_reused_from", A.overall.draws_reused_from)
    _cdid_rb_number(O.payload, "agg.overall.draws_used", A.overall.inference.draws_used)
    _cdid_rb_number(O.payload, "agg.overall.n_draws", A.overall.inference.n_draws)
    _cdid_rb_text(O.payload, "agg.inference_reason", A.inference.reason)
    _cdid_rb_text(O.payload, "agg.band_family", A.inference.band_family)
    _cdid_rb_number(O.payload, "agg.na_rm", Q.na_rm)
    _cdid_rb_number(O.payload, "agg.balance_e", Q.balance_e)
    _cdid_rb_number(O.payload, "agg.min_e", Q.min_e)
    _cdid_rb_number(O.payload, "agg.max_e", Q.max_e)
    _cdid_rb_matrix(O.payload, "agg.request_cohort_grid", Q.cohort_grid)
    _cdid_rb_matrix(O.payload, "agg.n_intended", A.n_intended)
    _cdid_rb_matrix(O.payload, "agg.n_used", A.n_used)
    _cdid_rb_matrix(O.payload, "agg.share_dependent", A.share_dependent)
    _cdid_rb_matrix(O.payload, "agg.share_denominator", A.share_denominator)
    _cdid_rb_matrix(O.payload, "agg.overall_weights", A.overall.weights)
    _cdid_rb_matrix(O.payload, "agg.overall_share_gradient", A.overall.share_gradient)
    _cdid_rb_matrix(O.payload, "agg.overall_intended", A.overall.intended)
    _cdid_rb_matrix(O.payload, "agg.overall_used", A.overall.used)
    _cdid_rb_matrix(O.payload, "agg.overall_excluded", A.overall.excluded)
    _cdid_rb_matrix(O.payload, "agg.overall_analytic_se", A.overall.inference.analytic_se)
    _cdid_rb_matrix(O.payload, "agg.overall_V_estimator", A.overall.inference.covariance)
    _cdid_rb_matrix(O.payload, "agg.overall_V_asymptotic", A.overall.inference.Omega)
    _cdid_rb_matrix(O.payload, "agg.fixed_se", A.fixed_inference.se)
    _cdid_rb_matrix(O.payload, "agg.Omega_fixed", A.decomposition.Omega_fixed)
    _cdid_rb_matrix(O.payload, "agg.Omega_share", A.decomposition.Omega_share)
    _cdid_rb_matrix(O.payload, "agg.Omega_cross", A.decomposition.Omega_cross)
    _cdid_rb_matrix(O.payload, "agg.Omega_total", A.decomposition.Omega_total)
    _cdid_rb_matrix(O.payload, "agg.Omega_theta_p", A.decomposition.Omega_theta_p)
    _cdid_rb_matrix(O.payload, "agg.points", A.estimates)
    _cdid_rb_matrix(O.payload, "agg.keys", A.keys)
    _cdid_rb_matrix(O.payload, "agg.available", A.cell_available)
    _cdid_rb_matrix(O.payload, "agg.restricted", A.restricted)
    _cdid_rb_matrix(O.payload, "agg.se", A.inference.se)
    _cdid_rb_matrix(O.payload, "agg.analytic_se", A.inference.analytic_se)
    _cdid_rb_matrix(O.payload, "agg.point_ci", (A.inference.ci_lower, A.inference.ci_upper))
    _cdid_rb_matrix(O.payload, "agg.band_ci", (A.inference.band_lower, A.inference.band_upper))
    _cdid_rb_matrix(O.payload, "agg.V_estimator", A.inference.covariance)
    _cdid_rb_matrix(O.payload, "agg.V_asymptotic", A.inference.Omega)
    _cdid_rb_matrix(O.payload, "agg.W", A.W)
    _cdid_rb_matrix(O.payload, "agg.D", A.D)
    _cdid_rb_matrix(O.payload, "agg.intended", A.intended)
    _cdid_rb_matrix(O.payload, "agg.used", A.used)
    _cdid_rb_matrix(O.payload, "agg.excluded", A.excluded)
    _cdid_rb_matrix(O.payload, "agg.overall", (A.overall.estimate, A.overall.se, A.overall.ci_lower, A.overall.ci_upper))
    _cdid_rb_matrix(O.payload, "agg.share_cohorts", A.shares.cohorts)
    _cdid_rb_matrix(O.payload, "agg.share_p", A.shares.p)
    _cdid_rb_unit_matrix(O.payload, "agg.IF_total", A.if_total)
    _cdid_rb_unit_matrix(O.payload, "agg.IF_fixed", A.if_fixed)
    _cdid_rb_unit_matrix(O.payload, "agg.IF_share", A.if_share)
    _cdid_rb_unit_matrix(O.payload, "agg.cohort_share_IF", A.shares.P)
    _cdid_rb_unit_matrix(O.payload, "agg.overall_IF", A.overall.if_full)
    for (j = 1; j <= A.M; j++) {
        _cdid_rb_text(O.payload, "agg.label." + _cdid_rb_key(j), A.labels[j])
        _cdid_rb_text(O.payload, "agg.reason." + _cdid_rb_key(j), A.cell_reason[j])
        _cdid_rb_text(O.payload, "agg.role." + _cdid_rb_key(j), A.bin_role[j])
        _cdid_rb_text(O.payload, "agg.support." + _cdid_rb_key(j), A.bin_support_status[j])
    }
    if (R.posting_mode == "state") {
        O.b = J(1, 0, .)
        O.V = O.Omega = J(0, 0, .)
    }
    draws_used = (A.inference.draws_used | A.overall.inference.draws_used)
    _cdid_rb_text(O.payload, "rng.mode", (draws_used ? (sign_mode == "caller" ? "caller_rng_stream_Rademacher" : "caller_fixed_signs") : "not_used"))
    _cdid_rb_number(O.payload, "rng.draws_used", draws_used)
    _cdid_rb_number(O.payload, "rng.n_draws", A.inference.n_draws)
    if (draws_used) {
        _cdid_rb_number(O.payload, "rng.batch_rows", min((Q.B, R.budget.max_sign_batch_rows, floor(R.budget.max_sign_batch_cells / R.C))))
    }
    _cdid_rb_common_metadata(O, R)
    cdid_state_validate(O.payload, R.limits)
    return(O)
}

struct cdid_bridge_result scalar cdid_public_agg_stream(
    string scalar token,
    struct cdid_agg_request scalar Q,
    struct cdid_bridge_request scalar R,
    string scalar sign_mode,
    real matrix fixed_signs)
{
    return(_cdid_pb_agg_core(token, Q, R, sign_mode, fixed_signs, .))
}

string scalar cdid_chained_agg_stage(
    string scalar token,
    string scalar type,
    real scalar min_e,
    real scalar max_e,
    real scalar balance_e,
    real scalar na_rm,
    string scalar bstrap_str,
    real scalar cband_req,
    string scalar biters_str,
    real scalar level,
    real scalar workspace,
    real scalar chunk,
    string scalar posting,
    string scalar requested_cluster,
    string scalar request_name,
    string scalar agg_name,
    string scalar req_use_name,
    string scalar agg_use_name,
    string scalar signs_name,
    string scalar sign_mode,
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
    struct cdid_bridge_request scalar R
    struct cdid_agg_request scalar Q
    struct cdid_infer_flags scalar bin_flags
    struct cdid_bridge_result scalar O
    real matrix signs, ov, ov_ase
    real scalar needs_signs, j, draws_used, n_draws, ov_draws_used, ov_n_draws
    string scalar endpoint, rng_mode
    string matrix stripe
    string colvector retained_settings
    pointer() scalar ptr

    if (matrix_limit < 1 | matrix_limit >= . | matrix_limit != floor(matrix_limit)) {
        _cdid_rb_fail("Current standard Stata matrix limit was not supplied.", 498)
    }
    if (sign_mode != "caller" & sign_mode != "fixed") {
        _cdid_rb_fail("Public sign source must be caller or private fixed signs.", 498)
    }

    if (request_name == "" & agg_name == "") {
        cdid_chained_agg_requests(token, type, min_e, max_e, balance_e, na_rm,
            bstrap_str, cband_req, biters_str, level, workspace, chunk, posting,
            requested_cluster, req_use_name, agg_use_name)
    }
    ptr = findexternal(req_use_name)
    if (ptr == NULL) _cdid_rb_fail("Explicit Mata bridge request not found.", 111)
    R = *ptr
    ptr = findexternal(agg_use_name)
    if (ptr == NULL) _cdid_rb_fail("Explicit aggregation request not found.", 111)
    Q = *ptr

    needs_signs = 0
    if (R.compute_inffunc) {
        endpoint = (Q.type == "simple" ? "simple" : "aggte_bins")
        bin_flags = cdid_infer_resolve_flags(Q.bstrap, Q.cband, endpoint)
        needs_signs = (bin_flags.draws | Q.bstrap == 1)
    }
    signs = J(0, 0, .)
    if (needs_signs & sign_mode == "fixed" & signs_name != "") {
        ptr = findexternal(signs_name)
        if (ptr == NULL) _cdid_rb_fail("Required supplied signs not found.", 111)
        signs = *ptr
    }

    O = _cdid_pb_agg_core(token, Q, R, sign_mode, signs, matrix_limit)
    if (cols(O.b) > matrix_limit) {
        _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.", 908)
    }

    stripe = J(cols(O.b), 2, "")
    for (j = 1; j <= cols(O.b); j++) {
        stripe[j, 2] = "agg_" + _cdid_rb_key(O.coef_map[j, 2])
    }
    if (cols(O.b)) {
        st_matrix(bname, O.b)
        st_matrixcolstripe(bname, stripe)
    }
    if (rows(O.V)) {
        st_matrix(vname, O.V)
        st_matrixcolstripe(vname, stripe)
        st_matrixrowstripe(vname, stripe)
        st_matrix(oname, O.Omega)
        st_matrixcolstripe(oname, stripe)
        st_matrixrowstripe(oname, stripe)
    }
    if (cols(O.b)) {
        st_matrix(sename, O.se)
        st_matrix(aname, O.analytic_se)
        st_matrix(ciname, O.ci)
        st_matrix(piname, O.point_ci)
        st_matrix(mapname, O.coef_map)
    }

    ov = _cdid_ca_read_matrix(O.payload, "agg.overall")
    ov_ase = _cdid_ca_read_matrix(O.payload, "agg.overall_analytic_se")
    rng_mode = _cdid_rb_stored_text(O.payload, "cdid_bridge.rng.mode")
    draws_used = _cdid_rb_stored_number(O.payload, "rng.draws_used")
    n_draws = _cdid_rb_stored_number(O.payload, "rng.n_draws")
    ov_draws_used = _cdid_rb_stored_number(O.payload, "agg.overall.draws_used")
    ov_n_draws = _cdid_rb_stored_number(O.payload, "agg.overall.n_draws")

    retained_settings = ("method" \ "gmm_weighting" \ "gmm_covariance" \ "base_period" \
        "control_group" \ "pt_assumption" \ "selection" \ "selection_canonical" \
        "selection_model" \ "selection_pooling" \ "selection_universe" \
        "selection_est_method" \ "imputation_pooling" \ "outcome_pooling" \ "weight_source" \
        "weight_reduction" \ "fix_weights_requested")
    for (j = 1; j <= rows(retained_settings); j++) {
        st_local("bridge_" + retained_settings[j], _cdid_ca_opt_text(O.payload, retained_settings[j], ""))
    }
    st_local("bridge_npost", _cdid_rb_key(cols(O.b)))
    st_local("bridge_hasV", _cdid_rb_key(rows(O.V) > 0))
    st_local("bridge_N_source", _cdid_rb_key(rows(O.payload.point)))
    st_local("bridge_alpha", strofreal(R.alpha, "%24.17g"))
    st_local("bridge_B", strofreal(R.B, "%24.17g"))
    st_local("bridge_bstrap", strofreal(R.bstrap, "%24.17g"))
    st_local("bridge_cband", _cdid_rb_key(R.cband))
    st_local("bridge_clustered", _cdid_rb_key(R.clustered))
    st_local("bridge_N", _cdid_rb_key(O.n_F))
    st_local("bridge_C", _cdid_rb_key(O.C))
    st_local("bridge_K", _cdid_rb_key(O.K))
    st_local("bridge_wald_C", strofreal(O.wald_C, "%24.17g"))
    st_local("bridge_posting_mode", R.posting_mode)
    st_local("bridge_kind", O.kind)
    st_local("bridge_status", O.status)
    st_local("bridge_inference_status", O.inference_status)
    st_local("bridge_flags", O.flag_law)
    st_local("bridge_warning", O.warning)
    st_local("bridge_wald_status", O.wald_status)
    st_local("bridge_wald_reason", O.wald_reason)
    st_local("bridge_critical", strofreal(O.critical, "%24.17g"))
    st_local("bridge_normal", strofreal(O.normal_critical, "%24.17g"))
    st_local("bridge_wald_df", strofreal(O.wald_df, "%24.17g"))
    st_local("bridge_wald_stat", strofreal(O.wald_statistic, "%24.17g"))
    st_local("bridge_wald_p", strofreal(O.wald_p, "%24.17g"))
    st_local("bridge_point_only", _cdid_rb_key(O.point_only))

    st_local("bridge_selector_warning", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.selector_warning"))
    st_local("bridge_agg_type", Q.type)
    st_local("bridge_target_def", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.target_definition"))
    st_local("bridge_band_family", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.band_family"))
    st_local("bridge_na_rm", _cdid_rb_key(Q.na_rm))
    st_local("bridge_min_e", strofreal(Q.min_e, "%24.17g"))
    st_local("bridge_max_e", strofreal(Q.max_e, "%24.17g"))
    st_local("bridge_balance_e", strofreal(Q.balance_e, "%24.17g"))
    st_local("bridge_overall_est", strofreal(ov[1, 1], "%24.17g"))
    st_local("bridge_overall_se", strofreal(ov[1, 2], "%24.17g"))
    st_local("bridge_overall_analytic_se", strofreal((rows(ov_ase) == 1 ? ov_ase[1, 1] : .), "%24.17g"))
    st_local("bridge_overall_ci_lo", strofreal(ov[1, 3], "%24.17g"))
    st_local("bridge_overall_ci_hi", strofreal(ov[1, 4], "%24.17g"))
    st_local("bridge_overall_status", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.overall.status"))
    st_local("bridge_overall_reason", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.overall.reason"))
    st_local("bridge_overall_inf_status", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.overall.inference_status"))
    st_local("bridge_overall_inf_kind", _cdid_rb_stored_text(O.payload, "cdid_bridge.agg.overall.inference_kind"))
    st_local("bridge_rng_mode", rng_mode)
    st_local("bridge_draws_used", _cdid_rb_key(draws_used))
    st_local("bridge_n_draws", _cdid_rb_key(n_draws))
    st_local("bridge_overall_draws_used", _cdid_rb_key(ov_draws_used))
    st_local("bridge_overall_n_draws", _cdid_rb_key(ov_n_draws))

    return(cdid_state_put(O.payload, R.limits))
}

void cdid_chained_agg_replay(string scalar token)
{
    struct cdid_result_state scalar P
    struct cdid_infer_flags scalar flags
    real scalar M, alpha, bstrap, cband, N, min_e, max_e, balance_e, na_rm
    real matrix ov, pts, av, se, band_ci
    string scalar agg_type, ov_status, ov_reason, wald_status, warning, endpoint

    if (strtrim(token) == "") {
        _cdid_rb_fail("Chained aggregation replay requires a state token.", 301)
    }
    P = cdid_state_get(token)
    if (_cdid_ca_opt_text(P, "cdid_bridge.kind", "") != "aggregation") {
        _cdid_rb_fail("State token does not contain a persisted chained aggregation result.", 301)
    }

    agg_type = _cdid_rb_stored_text(P, "cdid_bridge.agg.type")
    alpha = _cdid_rb_stored_number(P, "alpha")
    bstrap = _cdid_rb_stored_number(P, "bstrap")
    cband = _cdid_rb_stored_number(P, "cband")
    N = _cdid_rb_stored_number(P, "N")
    na_rm = _cdid_rb_stored_number(P, "agg.na_rm")
    balance_e = _cdid_rb_stored_number(P, "agg.balance_e")
    min_e = _cdid_rb_stored_number(P, "agg.min_e")
    max_e = _cdid_rb_stored_number(P, "agg.max_e")

    endpoint = (agg_type == "simple" ? "simple" : "aggte_bins")
    flags = cdid_infer_resolve_flags(bstrap, cband, endpoint)

    ov = _cdid_ca_read_matrix(P, "agg.overall")
    ov_status = _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.status")
    ov_reason = _cdid_rb_stored_text(P, "cdid_bridge.agg.overall.reason")
    pts = _cdid_ca_read_matrix(P, "agg.points")
    av = _cdid_ca_read_matrix(P, "agg.available")
    se = _cdid_ca_read_matrix(P, "agg.se")
    band_ci = _cdid_ca_read_matrix(P, "agg.band_ci")
    M = rows(pts)
    wald_status = _cdid_ca_opt_text(P, "cdid_bridge.wald_status", "not_computed")
    warning = _cdid_ca_opt_text(P, "cdid_bridge.warning", "")

    cdid_chained_agg_print(P, agg_type, alpha, flags, N, M, ov, ov_status, ov_reason, pts, av, se, band_ci,
        min_e, max_e, balance_e, na_rm, wald_status, warning)
}

/* Aggregation table shared by estat, replay and cdid_load. */
void cdid_chained_agg_print(struct cdid_result_state scalar P, string scalar agg_type, real scalar alpha,
    struct cdid_infer_flags scalar flags, real scalar N, real scalar M, real matrix ov, string scalar ov_status,
    string scalar ov_reason, real matrix pts, real matrix av, real matrix se, real matrix band_ci,
    real scalar min_e, real scalar max_e, real scalar balance_e, real scalar na_rm,
    string scalar wald_status, string scalar warning)
{
    real scalar j, level, k
    string scalar lbl, why, head, sel_warning, sel
    level = 100 * (1 - alpha)
    head = (agg_type == "group" ? "Group" : (agg_type == "dynamic" ? "Event time" : (agg_type == "calendar" ? "Time" : "")))
    _cdid_disp_chained_header(P, "Chained difference-in-differences", (agg_type == "simple" ? "" : "Number of " + (agg_type == "dynamic" ? "event times" : (agg_type == "group" ? "groups" : "periods"))), M, N)
    _cdid_disp_note("Aggregation: " + agg_type + (agg_type == "simple" ? " (average effect over all post-treatment cells)" : ""))
    sel = (min_e < . ? " mine(" + strofreal(min_e) + ")" : "") + (max_e < . ? " maxe(" + strofreal(max_e) + ")" : "") +
        (balance_e < . ? " balancee(" + strofreal(balance_e) + ")" : "") + (na_rm != 0 & na_rm < . ? " narm" : "")
    if (agg_type != "calendar" & sel != "") _cdid_disp_note("Event-time selection:" + sel)
    _cdid_disp_rule(78)
    printf("{txt}%-20s %10s %10s  %23s\n", head, "ATT", "Std. err.", _cdid_disp_interval_head(level, flags.simultaneous))
    _cdid_disp_rule(78)
    printf("{txt}%-20s {res}%10.7g %10.7g  %11.7g %11.7g\n", "Overall", ov[1, 1], ov[1, 2], ov[1, 3], ov[1, 4])
    if (ov_status != "available") _cdid_disp_unavailable(ov_reason)
    if (M > 0 & agg_type != "simple") {
        _cdid_disp_rule(78)
        for (j = 1; j <= M; j++) {
            lbl = _cdid_ca_opt_text(P, "cdid_bridge.agg.label." + _cdid_rb_key(j), "")
            k = strrpos(lbl, "_")
            if (k > 0) lbl = substr(lbl, k + 1, .)
            why = _cdid_ca_opt_text(P, "cdid_bridge.agg.reason." + _cdid_rb_key(j), "")
            if (why == "") why = _cdid_ca_opt_text(P, "cdid_bridge.agg.support." + _cdid_rb_key(j), "")
            printf("{txt}%20s {res}%10.7g %10.7g  %11.7g %11.7g\n", lbl, pts[j, 1], se[j, 1], band_ci[j, 1], band_ci[j, 2])
            if (!av[j, 1]) _cdid_disp_unavailable(why)
        }
    }
    _cdid_disp_rule(78)
    _cdid_disp_inference(0, flags.bootstrap, flags.simultaneous, level, "estimates")
    if (agg_type != "simple") _cdid_disp_note("The overall summary uses a pointwise interval.")
    _cdid_disp_wald(P, wald_status)
    if (strtrim(subinstr(warning, ";", "")) != "") _cdid_disp_note(warning)
    sel_warning = _cdid_ca_opt_text(P, "cdid_bridge.agg.selector_warning", "")
    if (sel_warning != "") _cdid_disp_note(sel_warning)
}


/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
transmorphic cdid_mrcp_state_excluded(transmorphic a1)
{
    return(cdid_mref_state_excluded(a1))
}

end
