*! Public custom-estimator chained binding and symbol validation, current StataNow19.5.
*! GPL-3 adaptation of R cdid chained_custom_ipw.R, chained_custom_selection.R,
*! chained_reference_routes.R, chained_selection_consumers.R.
*! Routes method(custom) callback(name) [payload(name)] across baseline and observation-weighted estimators.
version 19.5
mata:
real scalar _cdid_custom_is_ident(string scalar s)
{
    real scalar n, i, c
    n = strlen(s)
    if (n < 1 | n > 32) return(0)
    c = ascii(substr(s, 1, 1))
    if (!((c >= 65 & c <= 90) | (c >= 97 & c <= 122) | c == 95)) return(0)
    for (i = 2; i <= n; i++) {
        c = ascii(substr(s, i, 1))
        if (!((c >= 65 & c <= 90) | (c >= 97 & c <= 122) | (c >= 48 & c <= 57) | c == 95))
            return(0)
    }
    return(1)
}
void _cdid_custom_name_gate(string scalar callback_name, string scalar payload_name)
{
    if (!_cdid_custom_is_ident(callback_name)) {
        errprintf("method(custom) requires callback() naming a Mata function.\n")
        exit(198)
    }
    if (payload_name != "" & !_cdid_custom_is_ident(payload_name)) {
        errprintf("payload() must name a Mata global object.\n")
        exit(198)
    }
}
void cdid_custom_validate_symbols(string scalar callback_name, string scalar payload_name)
{
    pointer() scalar p_fn, p_pl
    if (!_cdid_custom_is_ident(callback_name)) {
        errprintf("method(custom) requires callback() naming a Mata function.\n")
        exit(198)
    }
    p_fn = findexternal(callback_name + "()")
    if (p_fn == NULL) {
        errprintf("Mata function %s() not found; define it before calling the command.\n", callback_name)
        exit(198)
    }
    if (payload_name != "") {
        if (!_cdid_custom_is_ident(payload_name)) {
            errprintf("payload() must name a Mata global object.\n")
            exit(198)
        }
        p_pl = findexternal(payload_name)
        if (p_pl == NULL) {
            errprintf("payload(%s) was not found as an active external Mata object.\n", payload_name)
            exit(198)
        }
    }
}
struct cdid_result_state scalar _cdid_custom_cb_report(
    struct cdid_public_bundle scalar B, struct cdid_cb_provider scalar fit,
    real colvector post_point, real colvector post_available, string colvector post_status,
    string colvector post_reason, real colvector post_est_ids, real matrix post_loadings,
    real scalar compute_if, struct cdid_state_limits scalar limits, real scalar chunk, real scalar cap)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    real colvector aux_point, aux_available, active, qids, postpos, auxpos, consumers, hit, uid
    real matrix values, postIF, auxIF
    real scalar Q, M, Jt, i, j, q, lo, hi, effective, perrow, retained
    Q = rows(B.request.aux); M = rows(fit.post_model); Jt = rows(post_est_ids)
    if (rows(post_loadings) != Jt | cols(post_loadings) != M | rows(fit.aux_model) != Q |
        chunk < 1 | chunk != floor(chunk))
        _cdid_pr_fail("cdid_native_reporting_input", "Custom baseline reporting dimensions are invalid.")
    aux_point = J(Q, 1, .); aux_available = (*fit.plan).aux_supported
    for (q = 1; q <= Q; q++) if (aux_available[q]) {
        if (fit.aux_model[q] < 1 | fit.aux_model[q] > fit.n_models)
            _cdid_pr_fail("cdid_native_reporting_input", "Supported auxiliary comparison lacks a custom model.")
        aux_point[q] = fit.points[fit.aux_model[q]]
    }
    buffer = cdid_public_report_start(B, post_point, post_available, post_status, post_reason,
        aux_point, aux_available, (*fit.plan).aux_reason, compute_if, limits)
    if (!compute_if) return(cdid_public_report_finish(buffer, limits))
    active = _cdid_rb_index(colsum(abs(post_loadings))' :> 0)
    if (Jt > 0 & any(fit.post_model[active] :< 1))
        _cdid_pr_fail("cdid_native_reporting_input", "Active GMM loading references an unsupported post link.")
    qids = J(0, 1, .)
    for (j = 1; j <= rows(B.request.report); j++) {
        if (B.request.report[j, 5] == 2 & buffer.payload.available[j])
            qids = qids \ B.request.report[j, 7]
    }
    if (rows(qids)) qids = uniqrows(sort(qids, 1))
    consumers = (rows(active) ? fit.post_model[active] : J(0, 1, .)) \
        (rows(qids) ? fit.aux_model[qids] : J(0, 1, .))
    if (rows(consumers)) consumers = uniqrows(sort(consumers, 1))
    postpos = J(rows(active), 1, .); auxpos = J(rows(qids), 1, .)
    for (i = 1; i <= rows(active); i++) { hit = _cdid_rb_index(consumers :== fit.post_model[active[i]]); postpos[i] = hit[1]; }
    for (i = 1; i <= rows(qids); i++) { hit = _cdid_rb_index(consumers :== fit.aux_model[qids[i]]); auxpos[i] = hit[1]; }
    retained = cdid_state_cost(buffer.payload)[1] + fit.workspace_entries + length(post_loadings)
    perrow = 16 * (rows(consumers) + 8) + 2 * (Jt + rows(qids) + rows(B.request.report))
    if (retained >= cap | perrow > cap - retained)
        _cdid_rs_fail("Custom baseline report IF replay exceeds explicit numeric workspace limits.", 908)
    effective = min((chunk, fit.chunk_rows, floor((cap - retained) / perrow)))
    if (effective < 1) _cdid_rs_fail("Custom baseline report IF replay admits no original-F row chunk.", 908)
    for (lo = 1; lo <= fit.n_F; lo = lo + effective) {
        hi = min((fit.n_F, lo + effective - 1)); uid = (lo::hi)
        values = (rows(consumers) ? cdid_cb_if_rows(fit, uid, consumers, "frame") : J(rows(uid), 0, .))
        postIF = (Jt == 0 ? J(rows(uid), 0, .) :
            (rows(active) ? values[, postpos'] * post_loadings[, active']' : J(rows(uid), Jt, 0)))
        auxIF = (rows(qids) ? values[, auxpos'] : J(rows(uid), 0, .))
        cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
    }
    out = cdid_public_report_finish(buffer, limits)
    return(out)
}
struct cdid_result_state scalar _cdid_custom_sc_report(
    struct cdid_public_bundle scalar B, struct cdid_sc_provider scalar fit,
    real colvector post_point, real colvector post_available, string colvector post_status,
    string colvector post_reason, real colvector post_est_ids, real matrix post_loadings,
    real scalar compute_if, struct cdid_state_limits scalar limits, real scalar chunk, real scalar cap)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    real colvector aux_point, aux_available, active, qids, postpos, auxpos, consumers, hit, uid
    real matrix values, postIF, auxIF
    real scalar Q, M, Jt, i, j, q, lo, hi, effective, perrow, retained
    Q = rows(B.request.aux); M = rows(fit.post_model); Jt = rows(post_est_ids)
    if (rows(post_loadings) != Jt | cols(post_loadings) != M | rows(fit.aux_model) != Q |
        chunk < 1 | chunk != floor(chunk))
        _cdid_pr_fail("cdid_native_reporting_input", "Custom selected reporting dimensions are invalid.")
    aux_point = J(Q, 1, .); aux_available = (*fit.plan).aux_supported
    for (q = 1; q <= Q; q++) if (aux_available[q]) {
        if (fit.aux_model[q] < 1 | fit.aux_model[q] > fit.n_models)
            _cdid_pr_fail("cdid_native_reporting_input", "Supported auxiliary comparison lacks a custom selected model.")
        aux_point[q] = fit.points[fit.aux_model[q]]
    }
    buffer = cdid_public_report_start(B, post_point, post_available, post_status, post_reason,
        aux_point, aux_available, (*fit.plan).aux_reason, compute_if, limits)
    if (!compute_if) return(cdid_public_report_finish(buffer, limits))
    active = _cdid_rb_index(colsum(abs(post_loadings))' :> 0)
    if (Jt > 0 & any(fit.post_model[active] :< 1))
        _cdid_pr_fail("cdid_native_reporting_input", "Active GMM loading references an unsupported post link.")
    qids = J(0, 1, .)
    for (j = 1; j <= rows(B.request.report); j++) {
        if (B.request.report[j, 5] == 2 & buffer.payload.available[j])
            qids = qids \ B.request.report[j, 7]
    }
    if (rows(qids)) qids = uniqrows(sort(qids, 1))
    consumers = (rows(active) ? fit.post_model[active] : J(0, 1, .)) \
        (rows(qids) ? fit.aux_model[qids] : J(0, 1, .))
    if (rows(consumers)) consumers = uniqrows(sort(consumers, 1))
    postpos = J(rows(active), 1, .); auxpos = J(rows(qids), 1, .)
    for (i = 1; i <= rows(active); i++) { hit = _cdid_rb_index(consumers :== fit.post_model[active[i]]); postpos[i] = hit[1]; }
    for (i = 1; i <= rows(qids); i++) { hit = _cdid_rb_index(consumers :== fit.aux_model[qids[i]]); auxpos[i] = hit[1]; }
    retained = cdid_state_cost(buffer.payload)[1] + fit.workspace_entries + length(post_loadings)
    perrow = 16 * (rows(consumers) + 8) + 2 * (Jt + rows(qids) + rows(B.request.report))
    if (retained >= cap | perrow > cap - retained)
        _cdid_rs_fail("Custom selected report IF replay exceeds explicit numeric workspace limits.", 908)
    effective = min((chunk, fit.chunk_rows, floor((cap - retained) / perrow)))
    if (effective < 1) _cdid_rs_fail("Custom selected report IF replay admits no original-F row chunk.", 908)
    for (lo = 1; lo <= fit.N; lo = lo + effective) {
        hi = min((fit.N, lo + effective - 1)); uid = (lo::hi)
        values = (rows(consumers) ? cdid_sc_if_rows(fit, uid, consumers, "frame") : J(rows(uid), 0, .))
        postIF = (Jt == 0 ? J(rows(uid), 0, .) :
            (rows(active) ? values[, postpos'] * post_loadings[, active']' : J(rows(uid), Jt, 0)))
        auxIF = (rows(qids) ? values[, auxpos'] : J(rows(uid), 0, .))
        cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
    }
    out = cdid_public_report_finish(buffer, limits)
    return(out)
}
struct cdid_result_state scalar _cdid_custom_scr_report(
    struct cdid_public_bundle scalar B, struct cdid_scr_provider scalar fit,
    real colvector post_point, real colvector post_available, string colvector post_status,
    string colvector post_reason, real colvector post_est_ids, real matrix post_loadings,
    real scalar compute_if, struct cdid_state_limits scalar limits, real scalar chunk, real scalar cap)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    real colvector aux_point, aux_available, active, qids, postpos, auxpos, consumers, hit, uid
    real matrix values, postIF, auxIF
    real scalar Q, M, Jt, i, j, q, lo, hi, effective, perrow, retained
    Q = rows(B.request.aux); M = rows(fit.post_model); Jt = rows(post_est_ids)
    if (rows(post_loadings) != Jt | cols(post_loadings) != M | rows(fit.aux_model) != Q |
        chunk < 1 | chunk != floor(chunk))
        _cdid_pr_fail("cdid_native_reporting_input", "Custom selected reference reporting dimensions are invalid.")
    aux_point = J(Q, 1, .); aux_available = (*fit.plan).aux_supported
    for (q = 1; q <= Q; q++) if (aux_available[q]) {
        if (fit.aux_model[q] < 1 | fit.aux_model[q] > fit.n_models)
            _cdid_pr_fail("cdid_native_reporting_input", "Supported auxiliary comparison lacks a custom selected reference model.")
        aux_point[q] = fit.points[fit.aux_model[q]]
    }
    buffer = cdid_public_report_start(B, post_point, post_available, post_status, post_reason,
        aux_point, aux_available, (*fit.plan).aux_reason, compute_if, limits)
    if (!compute_if) return(cdid_public_report_finish(buffer, limits))
    active = _cdid_rb_index(colsum(abs(post_loadings))' :> 0)
    if (Jt > 0 & any(fit.post_model[active] :< 1))
        _cdid_pr_fail("cdid_native_reporting_input", "Active GMM loading references an unsupported post link.")
    qids = J(0, 1, .)
    for (j = 1; j <= rows(B.request.report); j++) {
        if (B.request.report[j, 5] == 2 & buffer.payload.available[j])
            qids = qids \ B.request.report[j, 7]
    }
    if (rows(qids)) qids = uniqrows(sort(qids, 1))
    consumers = (rows(active) ? fit.post_model[active] : J(0, 1, .)) \
        (rows(qids) ? fit.aux_model[qids] : J(0, 1, .))
    if (rows(consumers)) consumers = uniqrows(sort(consumers, 1))
    postpos = J(rows(active), 1, .); auxpos = J(rows(qids), 1, .)
    for (i = 1; i <= rows(active); i++) { hit = _cdid_rb_index(consumers :== fit.post_model[active[i]]); postpos[i] = hit[1]; }
    for (i = 1; i <= rows(qids); i++) { hit = _cdid_rb_index(consumers :== fit.aux_model[qids[i]]); auxpos[i] = hit[1]; }
    retained = cdid_state_cost(buffer.payload)[1] + fit.workspace_entries + length(post_loadings)
    perrow = 16 * (rows(consumers) + 8) + 2 * (Jt + rows(qids) + rows(B.request.report))
    if (retained >= cap | perrow > cap - retained)
        _cdid_rs_fail("Custom selected reference report IF replay exceeds explicit numeric workspace limits.", 908)
    effective = min((chunk, fit.chunk_rows, floor((cap - retained) / perrow)))
    if (effective < 1) _cdid_rs_fail("Custom selected reference report IF replay admits no original-F row chunk.", 908)
    for (lo = 1; lo <= fit.N; lo = lo + effective) {
        hi = min((fit.N, lo + effective - 1)); uid = (lo::hi)
        values = (rows(consumers) ? cdid_scr_if_rows(fit, uid, consumers, "frame") : J(rows(uid), 0, .))
        postIF = (Jt == 0 ? J(rows(uid), 0, .) :
            (rows(active) ? values[, postpos'] * post_loadings[, active']' : J(rows(uid), Jt, 0)))
        auxIF = (rows(qids) ? values[, auxpos'] : J(rows(uid), 0, .))
        cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
    }
    out = cdid_public_report_finish(buffer, limits)
    return(out)
}
/* Records the identification chains and GMM comparison weights of a custom fit for
   cdid_chain. kept: comparisons used by GMM; points: estimates of the custom models. */
void _cdid_custom_chain_record(struct cdid_public_bundle scalar B, struct cdid_cat scalar C,
    struct cdid_design scalar D, real colvector kept, real colvector points, real colvector post_model,
    real colvector aux_model, real colvector aux_supported, real matrix loadings, real colvector est_ids)
{
    real colvector ap
    real scalar q
    ap = J(rows(B.request.aux), 1, .)
    for (q = 1; q <= rows(ap); q++) if (aux_supported[q]) ap[q] = points[aux_model[q]]
    cdid_chain_sp_record(B, C, D, kept, points[post_model], loadings, est_ids, ap, aux_supported)
}

struct cdid_result_state scalar _cdid_custom_fit_baseline(
    struct cdid_public_bundle scalar B, struct cdid_public_binding_control scalar ctrl,
    pointer(function) scalar cb_fn, pointer() scalar cb_payload)
{
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_population scalar owner
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar catalog
    struct cdid_cb_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_chained_ipw_postings scalar post
    struct cdid_chained_ipw_opts scalar opts
    struct cdid_cipw_row_basis scalar basis
    struct cdid_cb_weight_context scalar context
    struct cdid_gmm_result scalar sol
    struct cdid_result_state scalar out
    string scalar route
    real matrix G, loadings, postmeta, om
    real colvector point, diagnostic, reportpost, keep, ac, K
    string colvector diagnostic_names
    real scalar M, Jt, used, remaining, i, j, op, pid, q, sgn
    route = (B.prepared.design == "reference_recipes" ? "D60" : "static")
    co = cdid_catalog_default_opts(); co.mem_limit_bytes = ctrl.workspace_bytes; co.work_cap = ctrl.workcap
    co.conventions = "D10_D13_admissible_frame"
    pairs = cdid_cipw_pairs_opts(B.prepared.population, B.request.cohorts, co)
    C = cdid_catalog_build_opts(B.request.calendar, B.request.cohorts, B.request.target_cohorts,
        B.request.anticipation, pairs.occupancy, B.request.control_group, B.request.post, co)
    D = cdid_chained_ipw_build_design(C, B.request.post, B.request.pt_assumption, B.request.pt_start)
    B.request = cdid_public_request_tags(B.request, C)
    plan = cdid_baseline_plan(B.prepared.population, C, B.request.aux, B.request.aux_K_ptr, B.request.aux_K_index, ctrl.entries, ctrl.workcap)
    owner = B.prepared.population; reference = B.reference; reference.population = &owner
    cdid_ref_activate(reference)
    if (route == "static" & reference.weight_kind == "dated") {
        reference.weight_kind = "static"
        reference.recipes[, 2] = J(rows(reference.recipes), 1, 0)
        reference.weight_keys = (1::reference.N)
        reference.weights = reference.V
    }
    catalog = cdid_reference_catalog_make(plan, reference, ctrl.entries)
    fit = cdid_cb_prepare(&owner, &reference, &plan, &catalog, ctrl.baseline, route, cb_fn, cb_payload)
    sys = cdid_baseline_make_system(C, D, B.request.post, ctrl.entries, fit.workspace_entries)
    support = cdid_public_post_ledger(C, D, B.request, sys)
    M = sys.M; Jt = rows(B.request.post); point = J(Jt, 1, .); loadings = J(0, M, .)
    diagnostic = J(0, 1, .); diagnostic_names = J(0, 1, "")
    if (M > 0 & sys.p_retained > 0) {
        used = fit.workspace_entries + sys.workspace_entries
        if (B.settings.gmm_weighting == "identity") {
            used = _cdid_gmm_identity_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries
            sol = cdid_gmm_identity_map(fit.points[fit.post_model], sys.full_A, sys.H_est, remaining)
        }
        else {
            if (route == "D60") {
                keep = (1::M)
                G = I(M)
            }
            else {
                opts = cdid_chained_ipw_default_opts()
                opts.gmm_weighting = B.settings.gmm_weighting
                opts.reduce_structural_redundancy = 1
                opts.max_matrix_entries = ctrl.entries
                opts.work_cap = ctrl.workcap
                post = _cdid_cipw_build_postings(owner, C.calendar)
                ac = J(owner.n_A, 1, 0)
                for (j = 1; j <= owner.n_A; j++) ac[j] = _cdid_cohort_index(C.cohorts, C.nt_index, owner.cohort[owner.A_to_F[j]])
                om = J(C.n_cohorts, rows(C.op_g), 0)
                for (op = 1; op <= rows(C.op_g); op++) {
                    pid = C.op_pool[op]
                    K = C.pool_entries[C.pool_start[pid]..(C.pool_start[pid] + C.pool_size[pid] - 1)]
                    for (j = 1; j <= 1 + rows(K); j++) {
                        q = (j == 1 ? C.op_g[op] : K[j - 1])
                        om[q, op] = j
                    }
                }
                basis = _cdid_cipw_row_basis(C, D, post, ac, om, opts, 0, 0)
                keep = basis.retained
                G = (rows(basis.G) ? basis.G : I(rows(keep)))
            }
            used = _cdid_base_cells(2 * ctrl.chunk, rows(keep), ctrl.entries, used)
            used = _cdid_base_cells(3 * rows(keep), rows(keep), ctrl.entries, used)
            used = _cdid_gmm_workspace(rows(keep), sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries - rows(keep) * rows(keep)
            context = cdid_cb_weight_context_make(&fit, fit.post_model[keep], J(rows(keep), 1, 1))
            sol = cdid_gmm_fit_unit_norm(fit.points[fit.post_model[keep]], sys.full_A[keep,], G, M, sys.H_est,
                B.settings.gmm_weighting, ctrl.condition_trigger, ctrl.ridge_relative,
                fit.weight_rows, fit.normalization_count, ctrl.chunk, &cdid_cb_weight(), &context, remaining)
        }
        if (sys.n_estimable) point[sys.estimable_ids] = sol.targets
        if (B.settings.gmm_weighting == "identity" | route == "D60") loadings = sol.target_loadings
        else {
            loadings = J(rows(sol.target_loadings), M, 0)
            if (rows(keep)) loadings[, keep'] = sol.target_loadings
        }
        _cdid_custom_chain_record(B, C, D, ((B.settings.gmm_weighting == "identity" | route == "D60") ? (1::M) : keep),
            fit.points, fit.post_model, fit.aux_model, (*fit.plan).aux_supported, loadings, sys.estimable_ids)
        diagnostic = (sol.criterion \ sol.lambda \ sol.rho \ sol.trace_B \ sol.covariance_rank \
            sol.source_n_units \ sol.normalization_factor \ sol.sample_rank_upper)
        diagnostic_names = ("GMM_criterion" \ "GMM_lambda" \ "GMM_ridge_rho" \ "GMM_trace_B" \
            "GMM_covariance_rank" \ "GMM_source_replay_rows" \ "GMM_gram_normalization_factor" \ "GMM_sample_rank_upper")
        G = J(0, 0, .); sol.unit_gram_covariance = sol.coefficient_loadings = sol.target_loadings = J(0, 0, .)
        sol.coefficient_covariance = sol.covariance = J(0, 0, .)
    }
    else if (sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input", "An estimable post target has no GMM parameter system.")
    sys.full_A = sys.H_est = J(0, 0, .)
    out = _cdid_custom_cb_report(B, fit, point, sys.post_available, support.status, support.reason,
        sys.estimable_ids, loadings, ctrl.inference.compute_inffunc, ctrl.inference.limits, ctrl.chunk, ctrl.entries)
    postmeta = J(rows(B.request.report), cols(support.ledger), .)
    reportpost = _cdid_rb_index(B.request.report[, 5] :== 3)
    for (i = 1; i <= rows(reportpost); i++) { j = reportpost[i]; postmeta[j,] = support.ledger[B.request.report[j, 6],]; }
    out.ledger = out.ledger, postmeta
    out.ledger_columns = out.ledger_columns, "post_support." :+ support.columns
    out.numeric_settings = out.numeric_settings \ fit.n_models \ fit.n_callbacks \ diagnostic \
        ctrl.workspace_bytes \ ctrl.entries \ ctrl.chunk \ ctrl.workcap \ ctrl.condition_trigger \ ctrl.ridge_relative \
        ctrl.links_bound \ ctrl.J_bound
    out.numeric_setting_names = out.numeric_setting_names \ "n_custom_models" \ "n_callbacks" \ diagnostic_names \
        "numeric_workspace_bytes" \ "numeric_workspace_entries" \ "chunk_rows" \ "counted_work_limit" \
        "condition_trigger" \ "ridge_relative" \ "catalogue_links_bound" \ "inference_targets_bound"
    out.text_settings = out.text_settings \ ("custom_route", "baseline_" + route)
    out.numeric_settings = out.numeric_settings \ fit.weight_rows \ fit.normalization_count \ 0 \ 0
    out.numeric_setting_names = out.numeric_setting_names \ "GMM_replay_rows" \ "GMM_normalization_count" \ "n_selection_fits" \ "n_treatment_fits"
    out.text_settings = out.text_settings \ ("GMM_weight_row_domain", fit.weight_row_domain)
    out.source_versions = out.source_versions \ ("estimator", "user-supplied (method(custom)) baseline estimator; influence functions over all units")
    return(out)
}
struct cdid_result_state scalar _cdid_custom_fit_selected_static(
    struct cdid_public_bundle scalar B, struct cdid_public_binding_control scalar ctrl,
    pointer(function) scalar cb_fn, pointer() scalar cb_payload)
{
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_population scalar owner
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar catalog
    struct cdid_sel_decl scalar declaration
    struct cdid_sel_reader scalar reader
    struct cdid_selected_controls scalar ctl
    struct cdid_sc_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_sc_weight_context scalar context
    struct cdid_gmm_result scalar sol
    struct cdid_result_state scalar out
    real matrix X_A, G, loadings, postmeta
    real colvector point, diagnostic, reportpost, w_static
    string colvector diagnostic_names
    real scalar M, Jt, used, remaining, i, j, has_w
    co = cdid_catalog_default_opts(); co.mem_limit_bytes = ctrl.workspace_bytes; co.work_cap = ctrl.workcap
    co.conventions = "D10_D13_admissible_frame"
    pairs = cdid_cipw_pairs_opts(B.prepared.population, B.request.cohorts, co)
    C = cdid_catalog_build_opts(B.request.calendar, B.request.cohorts, B.request.target_cohorts,
        B.request.anticipation, pairs.occupancy, B.request.control_group, B.request.post, co)
    D = cdid_chained_ipw_build_design(C, B.request.post, B.request.pt_assumption, B.request.pt_start)
    B.request = cdid_public_request_tags(B.request, C)
    plan = cdid_baseline_plan(B.prepared.population, C, B.request.aux, B.request.aux_K_ptr, B.request.aux_K_index, ctrl.entries, ctrl.workcap)
    owner = B.prepared.population; reference = B.reference; reference.population = &owner
    cdid_ref_activate(reference)
    catalog = cdid_reference_catalog_make(plan, reference, ctrl.entries)
    has_w = (B.prepared.weight.source != "unity" & B.prepared.weight.source != "")
    w_static = (rows(B.prepared.weight.scientific_static) == owner.N ? B.prepared.weight.scientific_static : B.prepared.weight.V)
    declaration = cdid_sel_norm_decl(B.prepared.selection, !B.prepared.selection_formula_inherited,
        B.prepared.selection_universe, B.prepared.sel_has_rho, (strtrim(B.prepared.selection_justification) != ""), B.prepared.selection_model,
        B.prepared.selection_pooling, B.prepared.selection_justification)
    reader = cdid_sel_reader_pop(owner, B.prepared.Z_sel, B.prepared.sel_has_intercept,
        (B.prepared.sel_has_intercept ? 1 : 0), w_static, declaration, B.prepared.sel_rho, ctrl.chunk)
    ctl = cdid_selected_default_controls()
    ctl.est_method = "ipw"; ctl.eps = ctrl.baseline.glm_epsilon; ctl.maxit = ctrl.baseline.glm_maxit
    ctl.selection_model = B.prepared.selection_model; ctl.selection_pooling = B.prepared.selection_pooling
    ctl.rank_tol = ctrl.baseline.rank_tol; ctl.max_chunk_rows = ctrl.chunk
    ctl.max_work = (ctrl.workcap < . ? ctrl.workcap : 0); ctl.max_workspace_doubles = ctrl.entries
    X_A = B.prepared.X[owner.A_to_F,]
    fit = cdid_sc_prepare(&owner, &reader, &plan, &catalog, X_A, B.prepared.x_names, declaration, ctl, cb_fn, cb_payload)
    sys = cdid_baseline_make_system(C, D, B.request.post, ctrl.entries, fit.workspace_entries)
    support = cdid_public_post_ledger(C, D, B.request, sys)
    M = sys.M; Jt = rows(B.request.post); point = J(Jt, 1, .); loadings = J(0, M, .)
    diagnostic = J(0, 1, .); diagnostic_names = J(0, 1, "")
    if (M > 0 & sys.p_retained > 0) {
        used = fit.workspace_entries + sys.workspace_entries
        if (B.settings.gmm_weighting == "identity") {
            used = _cdid_gmm_identity_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries
            sol = cdid_gmm_identity_map(fit.points[fit.post_model], sys.full_A, sys.H_est, remaining)
        }
        else {
            used = _cdid_base_cells(2 * ctrl.chunk, M, ctrl.entries, used)
            used = _cdid_base_cells(3 * M, M, ctrl.entries, used)
            used = _cdid_gmm_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries - M * M
            G = I(M)
            context = cdid_sc_weight_context_make(&fit, fit.post_model, J(M, 1, 1))
            sol = cdid_gmm_fit_unit_norm(fit.points[fit.post_model], sys.full_A, G, M, sys.H_est,
                B.settings.gmm_weighting, ctrl.condition_trigger, ctrl.ridge_relative,
                fit.weight_rows, fit.normalization_count, ctrl.chunk, &cdid_sc_weight(), &context, remaining)
        }
        if (sys.n_estimable) point[sys.estimable_ids] = sol.targets
        loadings = sol.target_loadings
        _cdid_custom_chain_record(B, C, D, (1::M), fit.points, fit.post_model, fit.aux_model,
            (*fit.plan).aux_supported, loadings, sys.estimable_ids)
        diagnostic = (sol.criterion \ sol.lambda \ sol.rho \ sol.trace_B \ sol.covariance_rank \
            sol.source_n_units \ sol.normalization_factor \ sol.sample_rank_upper)
        diagnostic_names = ("GMM_criterion" \ "GMM_lambda" \ "GMM_ridge_rho" \ "GMM_trace_B" \
            "GMM_covariance_rank" \ "GMM_source_replay_rows" \ "GMM_gram_normalization_factor" \ "GMM_sample_rank_upper")
        G = J(0, 0, .); sol.unit_gram_covariance = sol.coefficient_loadings = sol.target_loadings = J(0, 0, .)
        sol.coefficient_covariance = sol.covariance = J(0, 0, .)
    }
    else if (sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input", "An estimable post target has no GMM parameter system.")
    sys.full_A = sys.H_est = J(0, 0, .)
    out = _cdid_custom_sc_report(B, fit, point, sys.post_available, support.status, support.reason,
        sys.estimable_ids, loadings, ctrl.inference.compute_inffunc, ctrl.inference.limits, ctrl.chunk, ctrl.entries)
    postmeta = J(rows(B.request.report), cols(support.ledger), .)
    reportpost = _cdid_rb_index(B.request.report[, 5] :== 3)
    for (i = 1; i <= rows(reportpost); i++) { j = reportpost[i]; postmeta[j,] = support.ledger[B.request.report[j, 6],]; }
    out.ledger = out.ledger, postmeta
    out.ledger_columns = out.ledger_columns, "post_support." :+ support.columns
    out.numeric_settings = out.numeric_settings \ fit.n_models \ fit.n_callbacks \ diagnostic \
        ctrl.workspace_bytes \ ctrl.entries \ ctrl.chunk \ ctrl.workcap \ ctrl.condition_trigger \ ctrl.ridge_relative \
        ctrl.links_bound \ ctrl.J_bound
    out.numeric_setting_names = out.numeric_setting_names \ "n_custom_models" \ "n_callbacks" \ diagnostic_names \
        "numeric_workspace_bytes" \ "numeric_workspace_entries" \ "chunk_rows" \ "counted_work_limit" \
        "condition_trigger" \ "ridge_relative" \ "catalogue_links_bound" \ "inference_targets_bound"
    out.text_settings = out.text_settings \ ("custom_route", "selected_static")
    out.numeric_settings = out.numeric_settings \ fit.weight_rows \ fit.normalization_count \ fit.n_selection_fits \ fit.n_treatment_fits
    out.numeric_setting_names = out.numeric_setting_names \ "GMM_replay_rows" \ "GMM_normalization_count" \ "n_selection_fits" \ "n_treatment_fits"
    out.text_settings = out.text_settings \ ("GMM_weight_row_domain", fit.weight_row_domain)
    out.source_versions = out.source_versions \ ("estimator", "user-supplied (method(custom)) estimator with observation weighting; influence functions over all units")
    return(out)
}
struct cdid_result_state scalar _cdid_custom_fit_selref(
    struct cdid_public_bundle scalar B, struct cdid_public_binding_control scalar ctrl,
    pointer(function) scalar cb_fn, pointer() scalar cb_payload)
{
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_population scalar owner
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar catalog
    struct cdid_selref_controls scalar selected_controls
    struct cdid_scr_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_scr_weight_context scalar context
    struct cdid_gmm_result scalar sol
    struct cdid_result_state scalar out
    real matrix G, loadings, postmeta
    real colvector point, diagnostic, reportpost
    string colvector diagnostic_names
    real scalar M, Jt, used, remaining, i, j
    co = cdid_catalog_default_opts(); co.mem_limit_bytes = ctrl.workspace_bytes; co.work_cap = ctrl.workcap
    co.conventions = "D10_D13_admissible_frame"
    pairs = cdid_cipw_pairs_opts(B.prepared.population, B.request.cohorts, co)
    C = cdid_catalog_build_opts(B.request.calendar, B.request.cohorts, B.request.target_cohorts,
        B.request.anticipation, pairs.occupancy, B.request.control_group, B.request.post, co)
    D = cdid_chained_ipw_build_design(C, B.request.post, B.request.pt_assumption, B.request.pt_start)
    B.request = cdid_public_request_tags(B.request, C)
    plan = cdid_baseline_plan(B.prepared.population, C, B.request.aux, B.request.aux_K_ptr, B.request.aux_K_index, ctrl.entries, ctrl.workcap)
    owner = B.prepared.population; reference = B.reference; reference.population = &owner
    cdid_ref_activate(reference)
    catalog = cdid_reference_catalog_make(plan, reference, ctrl.entries)
    selected_controls = cdid_selref_default_controls(owner.N)
    selected_controls.method = "ipw"; selected_controls.selection_model = B.prepared.selection_model
    selected_controls.selection_pooling = B.prepared.selection_pooling
    selected_controls.universe = B.prepared.selection_universe
    selected_controls.justification = B.prepared.selection_justification
    selected_controls.has_rho = B.prepared.sel_has_rho; selected_controls.rho = B.prepared.sel_rho
    selected_controls.eps = ctrl.baseline.glm_epsilon; selected_controls.maxit = ctrl.baseline.glm_maxit
    selected_controls.rank_tol = ctrl.baseline.rank_tol; selected_controls.chunk_rows = ctrl.chunk
    selected_controls.max_work = (ctrl.workcap < . ? ctrl.workcap : 0)
    selected_controls.max_entries = ctrl.entries
    selected_controls.points_only = 0
    fit = cdid_scr_prepare(&owner, &reference, &plan, &catalog, selected_controls, cb_fn, cb_payload)
    sys = cdid_baseline_make_system(C, D, B.request.post, ctrl.entries, fit.workspace_entries)
    support = cdid_public_post_ledger(C, D, B.request, sys)
    M = sys.M; Jt = rows(B.request.post); point = J(Jt, 1, .); loadings = J(0, M, .)
    diagnostic = J(0, 1, .); diagnostic_names = J(0, 1, "")
    if (M > 0 & sys.p_retained > 0) {
        used = fit.workspace_entries + sys.workspace_entries
        if (B.settings.gmm_weighting == "identity") {
            used = _cdid_gmm_identity_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries
            sol = cdid_gmm_identity_map(fit.points[fit.post_model], sys.full_A, sys.H_est, remaining)
        }
        else {
            used = _cdid_base_cells(2 * ctrl.chunk, M, ctrl.entries, used)
            used = _cdid_base_cells(3 * M, M, ctrl.entries, used)
            used = _cdid_gmm_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries - M * M
            G = I(M)
            context = cdid_scr_weight_context_make(&fit, fit.post_model, J(M, 1, 1))
            sol = cdid_gmm_fit_unit_norm(fit.points[fit.post_model], sys.full_A, G, M, sys.H_est,
                B.settings.gmm_weighting, ctrl.condition_trigger, ctrl.ridge_relative,
                fit.weight_rows, fit.normalization_count, ctrl.chunk, &cdid_scr_weight(), &context, remaining)
        }
        if (sys.n_estimable) point[sys.estimable_ids] = sol.targets
        loadings = sol.target_loadings
        _cdid_custom_chain_record(B, C, D, (1::M), fit.points, fit.post_model, fit.aux_model,
            (*fit.plan).aux_supported, loadings, sys.estimable_ids)
        diagnostic = (sol.criterion \ sol.lambda \ sol.rho \ sol.trace_B \ sol.covariance_rank \
            sol.source_n_units \ sol.normalization_factor \ sol.sample_rank_upper)
        diagnostic_names = ("GMM_criterion" \ "GMM_lambda" \ "GMM_ridge_rho" \ "GMM_trace_B" \
            "GMM_covariance_rank" \ "GMM_source_replay_rows" \ "GMM_gram_normalization_factor" \ "GMM_sample_rank_upper")
        G = J(0, 0, .); sol.unit_gram_covariance = sol.coefficient_loadings = sol.target_loadings = J(0, 0, .)
        sol.coefficient_covariance = sol.covariance = J(0, 0, .)
    }
    else if (sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input", "An estimable post target has no GMM parameter system.")
    sys.full_A = sys.H_est = J(0, 0, .)
    out = _cdid_custom_scr_report(B, fit, point, sys.post_available, support.status, support.reason,
        sys.estimable_ids, loadings, ctrl.inference.compute_inffunc, ctrl.inference.limits, ctrl.chunk, ctrl.entries)
    postmeta = J(rows(B.request.report), cols(support.ledger), .)
    reportpost = _cdid_rb_index(B.request.report[, 5] :== 3)
    for (i = 1; i <= rows(reportpost); i++) { j = reportpost[i]; postmeta[j,] = support.ledger[B.request.report[j, 6],]; }
    out.ledger = out.ledger, postmeta
    out.ledger_columns = out.ledger_columns, "post_support." :+ support.columns
    out.numeric_settings = out.numeric_settings \ fit.n_models \ fit.n_callbacks \ diagnostic \
        ctrl.workspace_bytes \ ctrl.entries \ ctrl.chunk \ ctrl.workcap \ ctrl.condition_trigger \ ctrl.ridge_relative \
        ctrl.links_bound \ ctrl.J_bound
    out.numeric_setting_names = out.numeric_setting_names \ "n_custom_models" \ "n_callbacks" \ diagnostic_names \
        "numeric_workspace_bytes" \ "numeric_workspace_entries" \ "chunk_rows" \ "counted_work_limit" \
        "condition_trigger" \ "ridge_relative" \ "catalogue_links_bound" \ "inference_targets_bound"
    out.text_settings = out.text_settings \ ("custom_route", "selected_reference")
    out.numeric_settings = out.numeric_settings \ fit.weight_rows \ fit.normalization_count \ fit.n_selection_fits \ fit.n_treatment_fits
    out.numeric_setting_names = out.numeric_setting_names \ "GMM_replay_rows" \ "GMM_normalization_count" \ "n_selection_fits" \ "n_treatment_fits"
    out.text_settings = out.text_settings \ ("GMM_weight_row_domain", fit.weight_row_domain)
    out.numeric_settings = out.numeric_settings \ fit.n_selection_blocks \ fit.n_selection_registry_fits
    out.numeric_setting_names = out.numeric_setting_names \ "n_selection_blocks" \ "n_selection_registry_fits"
    out.source_versions = out.source_versions \ ("estimator", "user-supplied (method(custom)) reference-date estimator with observation weighting; influence functions over all units")
    return(out)
}
void cdid_custom_public_store(string scalar bundle_name, string scalar request_name,
    string scalar callback_name, string scalar payload_name, string scalar gmm_weighting,
    real scalar include_intercept, real scalar workspace, real scalar chunk, real scalar workcap,
    real scalar condition, real scalar ridge, real scalar compute_if, real scalar bstrap,
    real scalar cband, real scalar biters, real scalar alpha, string scalar posting)
{
    pointer() scalar p, out, p_fn, p_pl
    struct cdid_public_bundle scalar B
    struct cdid_public_binding_control scalar ctrl
    struct cdid_result_state scalar state
    struct cdid_bridge_request scalar R
    string scalar token
    cdid_custom_validate_symbols(callback_name, payload_name)
    p_fn = findexternal(callback_name + "()")
    p_pl = (payload_name != "" ? findexternal(payload_name) : NULL)
    p = findexternal(bundle_name)
    if (p == NULL) _cdid_pr_fail("cdid_chained_request_input_error", "Prepared public bundle is absent.")
    if (findexternal(request_name) != NULL)
        _cdid_pr_fail("cdid_chained_request_input_error", "Requested bridge-request name is already in use.")
    B = *p
    if (posting != "matrix" & posting != "state")
        _cdid_pr_fail("cdid_chained_public_control_error", "posting must be matrix or state.")
    if (gmm_weighting != "optimal" & gmm_weighting != "two_step" & gmm_weighting != "identity")
        _cdid_pr_fail("cdid_chained_public_control_error", "gmmweighting() must be optimal, two_step, or identity.")
    B.settings.gmm_weighting = gmm_weighting
    ctrl = cdid_public_binding_controls(B, include_intercept, workspace, chunk, workcap,
        condition, ridge, compute_if, bstrap, cband, biters, alpha, posting)
    B.settings.alpha = ctrl.inference.alpha; B.settings.compute_if = ctrl.inference.compute_inffunc
    B.settings.bstrap = ctrl.inference.bstrap; B.settings.cband = ctrl.inference.cband; B.settings.biters = ctrl.inference.B
    B.prepared.method = "custom"; B.settings.method = "custom"
    if (B.prepared.selection_canonical == "none")
        state = _cdid_custom_fit_baseline(B, ctrl, p_fn, p_pl)
    else if (B.prepared.selection_canonical == "conditional_mtar" & B.prepared.design == "static_numeric")
        state = _cdid_custom_fit_selected_static(B, ctrl, p_fn, p_pl)
    else if (B.prepared.selection_canonical == "conditional_mtar" & B.prepared.design == "reference_recipes")
        state = _cdid_custom_fit_selref(B, ctrl, p_fn, p_pl)
    else _cdid_pr_fail("cdid_native_private_binding_scope", "Unsupported custom chained route.")
    state.text_settings = state.text_settings \ ("callback", callback_name) \ ("payload", payload_name)
    cdid_state_validate(state, ctrl.inference.limits)
    token = cdid_state_put(state, ctrl.inference.limits)
    R.budget = ctrl.inference.budget; R.limits = ctrl.inference.limits
    R.alpha = ctrl.inference.alpha; R.compute_inffunc = ctrl.inference.compute_inffunc
    R.bstrap = ctrl.inference.bstrap; R.cband = ctrl.inference.cband; R.B = ctrl.inference.B
    R.C = ctrl.inference.C; R.clustered = ctrl.inference.clustered; R.run_wald = ctrl.inference.run_wald
    R.unit_cluster = ctrl.inference.unit_cluster
    R.report_cells = ctrl.inference.report_cells; R.pre_cells = ctrl.inference.pre_cells
    R.posting_mode = posting
    out = crexternal(request_name)
    if (out == NULL) {
        cdid_state_drop(token)
        _cdid_pr_fail("cdid_chained_request_input_error", "Cannot create the explicit bridge request.")
    }
    *out = R
    st_local("fit_token", token)
}
void cdid_custom_pt_not_applicable(string scalar output, real scalar B, real scalar alpha,
    real scalar bs, real scalar cb, string scalar formula, string scalar clusters,
    string scalar signs, real scalar cap)
{
    struct cdid_pt_result scalar R
    R.mode = "chained"; R.method = "custom"
    R.status = "not_applicable"
    R.reason = "conditional pre-trends test is not supported with custom panel estimators"
    R.formula = formula; R.clusters = clusters
    R.rng = (signs == "" ? "caller" : "fixed")
    R.policy = "not available for method(custom) in chained mode, as in R"
    R.channel = "none"; R.membership = "none"
    R.stats = (., ., 0, ., B, alpha, ., ., ., 0, ., ., st_nobs())
    R.draws = J(0, 1, .)
    R.family_columns = ("cell", "group", "time", "from", "sign", "auxiliary", "unavailable",
        "design_index", "reference_date", "estimate", "scalar_process", "distinct_points", "coordinates")
    R.family = J(0, 13, .); R.family_reasons = J(0, 1, "")
    R.cell_ids = R.cell_status = J(0, 1, "")
    R.raw_grid_variables = cdid_ref_formula_names(formula)
    R.preprocess_columns = ("unit", "cohort", "cohort_known", "observed_rows", "A", "B", "frame_only")
    R.xcols = R.zcols = J(1, 0, "")
    R.id_labels = R.cluster_labels = R.grid_keys = J(0, 1, "")
    R.cluster_positions = R.cluster_values = R.retained_rows = R.grid_first = R.multiplicity = J(0, 1, .)
    R.rowmap = J(0, 4, .); R.preprocess = J(0, 7, .)
    R.missing_grid = R.unsupported_grid = J(0, 4, .)
    R.bstrap_requested = bs; R.cband_requested = cb; R.multipliers_used = 0; R.raw_CvM = .
    cdid_pt_put(output, R, cap)
}

/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
void cdid_mrcb_require_route(transmorphic a1)
{
    cdid_mref_require_route(a1)
}

end
