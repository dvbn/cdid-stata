*! Public static-numeric observation-weighted selection provider and fit glue; current StataNow 19.5.
*! Binds cdid_selected_prepare_registry() and cdid_selected_chunk() to common GMM/report/inference/state.
version 19.5
mata:

struct cdid_selected_profile_bundle {
    struct cdid_selected_registry scalar registry
    struct cdid_selected_recipe_spec colvector recipes
    real colvector consumer_recipe_id
    string scalar profile_kind
}

struct cdid_selected_provider {
    string scalar method, status, weight_row_domain
    real scalar n_F, n_A, normalization_count, weight_rows, chunk_rows, workspace_entries
    struct cdid_population scalar population
    struct cdid_baseline_plan_state scalar plan
    struct cdid_selected_profile_bundle scalar profiles
    struct cdid_selected_model scalar model
    real colvector points
    real colvector post_ids_sorted, post_map_sorted
    real scalar post_cache_ready
    real colvector eval_ids_sorted, eval_post_map
    real matrix cached_eval_phi
    real scalar cache_eval_active
}

struct cdid_sr_plan_cache {
    real scalar valid, N, R, n_A, entries_cap, anticipation
    string scalar control_group
    real colvector calendar, cohorts, target_cohorts, first_row, pop_cohort
    real matrix post, aux
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_baseline_plan_state scalar plan
}

real scalar _cdid_sel_try_cached_plan(
    struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl,
    struct cdid_pairs scalar pairs,
    struct cdid_cat scalar C,
    struct cdid_baseline_plan_state scalar plan)
{
    pointer(struct cdid_sr_plan_cache scalar) scalar p
    struct cdid_sr_plan_cache scalar c

    p = findexternal("_cdid_sr_cached_plan")
    if (p == NULL) return(0)
    c = *p
    rmexternal("_cdid_sr_cached_plan")
    if (!c.valid) return(0)
    if (ctrl.workcap < . | ctrl.entries != c.entries_cap | ctrl.workspace_bytes != 8 * c.entries_cap) return(0)
    if (c.N != B.prepared.population.N | c.R != B.prepared.population.R | c.n_A != B.prepared.population.n_A) return(0)
    if (c.anticipation != B.request.anticipation | c.control_group != B.request.control_group) return(0)
    if (c.calendar != B.request.calendar | c.cohorts != B.request.cohorts | c.target_cohorts != B.request.target_cohorts) return(0)
    if (c.post != B.request.post | c.aux != B.request.aux) return(0)
    if (c.first_row != B.prepared.population.first_row | c.pop_cohort != B.prepared.population.cohort) return(0)
    pairs = c.pairs
    C = c.C
    plan = c.plan
    return(1)
}

void _cdid_sel_pub_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

string scalar _cdid_sel_pub_func_key(struct cdid_baseline_plan_state scalar p, real scalar i)
{
    real colvector K
    real scalar j, g_val, s_val, t_val
    string scalar prefix, k_str

    prefix = (i <= p.M ? "cmpcat:CP04-post:" : "cmpcat:CP04-aux:")
    g_val = p.cohorts[p.raw[i, 1]]
    s_val = p.calendar[p.raw[i, 2]]
    t_val = p.calendar[p.raw[i, 3]]
    K = p.raw_K_index[p.raw_K_ptr[i]..(p.raw_K_ptr[i + 1] - 1)]
    k_str = ""
    for (j = 1; j <= rows(K); j++) {
        k_str = k_str + (j > 1 ? "," : "") + _cdid_pr_num(p.cohorts[K[j]])
    }
    return(prefix + "g=" + _cdid_pr_num(g_val) + "|s=" + _cdid_pr_num(s_val) +
        "|t=" + _cdid_pr_num(t_val) + "|K=" + k_str)
}

/* Clean profile-registry boundary: static-numeric profiles populate slot 1
   (calendar_index = 0) and recipe 1 ("static"). The subsequent dated-recipe
   adapter plugs multiple recipes and dated profile slots into this same struct
   without changing downstream provider preparation, GMM weighting, or report IF. */
struct cdid_selected_profile_bundle scalar cdid_selected_static_profiles(
    struct cdid_public_bundle scalar B,
    struct cdid_baseline_plan_state scalar plan)
{
    struct cdid_selected_profile_bundle scalar prof
    real matrix X_or
    real scalar N, p

    N = B.prepared.population.N
    p = cols(B.prepared.X)
    if (B.prepared.x_has_intercept & p < 1) {
        _cdid_sel_pub_fail("cdid_chained_request_covariate_error",
            "Declared treatment/outcome intercept requires column 1.", 198)
    }
    if (B.prepared.x_has_intercept & p > 1) X_or = B.prepared.X[, 2..p]
    else if (B.prepared.x_has_intercept) X_or = J(N, 0, .)
    else X_or = B.prepared.X

    prof.profile_kind = "static_numeric"
    prof.registry.sel_data = cdid_selected_matrix_slot(1, 1)
    prof.registry.sel_data[1].calendar_index = 0
    prof.registry.sel_data[1].value = B.prepared.Z_sel

    prof.registry.or_data = cdid_selected_matrix_slot(1, 1)
    prof.registry.or_data[1].calendar_index = 0
    prof.registry.or_data[1].value = X_or

    prof.registry.trt_data = cdid_selected_matrix_slot(1, 1)
    prof.registry.trt_data[1].calendar_index = 0
    prof.registry.trt_data[1].value = B.prepared.X

    prof.registry.weight_data = cdid_selected_matrix_slot(1, 1)
    prof.registry.weight_data[1].calendar_index = 0
    prof.registry.weight_data[1].value = B.prepared.weight.owner_weight

    prof.registry.rho_data = cdid_selected_matrix_slot(1, 1)
    prof.registry.rho_data[1].calendar_index = 0
    prof.registry.rho_data[1].value = B.prepared.sel_rho

    prof.recipes = cdid_selected_recipe_spec(1, 1)
    prof.recipes[1].recipe_id = 1
    prof.recipes[1].recipe_key = "static"
    prof.recipes[1].design_cal = 0
    prof.recipes[1].weight_cal = 0
    prof.recipes[1].sel_id = 1
    prof.recipes[1].or_id = 1
    prof.recipes[1].trt_id = 1
    prof.recipes[1].weight_id = 1
    prof.recipes[1].rho_id = 1
    prof.recipes[1].sel_has_intercept = B.prepared.sel_has_intercept
    prof.recipes[1].sel_int_col = (B.prepared.sel_has_intercept ? 1 : 0)
    prof.recipes[1].has_xformla = 1
    prof.recipes[1].has_rho = B.prepared.sel_has_rho
    prof.recipes[1].universe = B.prepared.selection_universe
    prof.recipes[1].justification = B.prepared.selection_justification

    prof.consumer_recipe_id = (plan.n_consumers > 0 ? J(plan.n_consumers, 1, 1) : J(0, 1, 1))
    return(prof)
}

/* Prepare the selected provider across all unique supported post + auxiliary
   consumers using the cdid_selected_prepare_registry() API. */
struct cdid_selected_provider scalar cdid_selected_provider_prepare(
    struct cdid_public_bundle scalar B,
    struct cdid_baseline_plan_state scalar plan,
    struct cdid_selected_profile_bundle scalar prof,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_selected_provider scalar prov
    struct cdid_selected_pair_in colvector pairs_in
    struct cdid_selected_link_spec colvector link_specs
    struct cdid_selected_controls scalar ctl
    real colvector cohort_idx, k_ctrl, hit, first_rep
    real scalar N, J, R_raw, pi, c, i, rep_i, s_cal, t_cal

    N = B.prepared.population.N
    J = rows(plan.cohorts)
    R_raw = plan.M + plan.Q
    if (B.prepared.selection_canonical == "conditional_mtar") {
        if (any(B.prepared.population.cohort_known :== 0)) {
            _cdid_sel_pub_fail("cdid_chained_request_cohort_error",
                "Corollary1 requires known cohorts on every F row, including frame-only rows.", 459)
        }
    }

    if (J > 1) {
        if (rows(uniqrows(sort(plan.cohorts, 1))) != J) {
            _cdid_sel_pub_fail("cdid_chained_request_cohort_error",
                "Known unit cohort is absent from catalogue cohort dictionary.", 498)
        }
    }
    cohort_idx = J(N, 1, 0)
    for (c = 1; c <= J; c++) {
        if (plan.cohorts[c] >= .) {
            hit = vec(selectindex((B.prepared.population.cohort_known :== 1) :& (B.prepared.population.cohort :>= .)))
        }
        else {
            hit = vec(selectindex((B.prepared.population.cohort_known :== 1) :& (B.prepared.population.cohort :== plan.cohorts[c])))
        }
        if (rows(hit) > 0) cohort_idx[hit] = J(rows(hit), 1, c)
    }
    if (any((B.prepared.population.cohort_known :== 1) :& (cohort_idx :== 0))) {
        _cdid_sel_pub_fail("cdid_chained_request_cohort_error",
            "Known unit cohort is absent from catalogue cohort dictionary.", 498)
    }

    prov.method = B.prepared.method
    prov.status = "ok"
    prov.weight_row_domain = "F on core scale psi=(n_A/N)*phi_F with 1/n_A normalization"
    prov.n_F = N
    prov.n_A = B.prepared.population.n_A
    prov.normalization_count = B.prepared.population.n_A
    prov.weight_rows = N
    prov.chunk_rows = ctrl.chunk
    prov.population = B.prepared.population
    prov.plan = plan
    prov.profiles = prof
    prov.points = J(plan.n_consumers, 1, .)
    prov.workspace_entries = plan.workspace_entries
    prov.post_ids_sorted = J(0, 1, 0)
    prov.post_map_sorted = J(0, 1, 0)
    prov.post_cache_ready = 0
    prov.eval_ids_sorted = J(0, 1, 0)
    prov.eval_post_map = J(0, 1, 0)
    prov.cached_eval_phi = J(0, 0, .)
    prov.cache_eval_active = 0

    if (plan.n_consumers == 0) {
        prov.model.prepared = 0
        prov.status = "empty_consumers"
        return(prov)
    }

    pairs_in = cdid_selected_pair_in(plan.n_pairs, 1)
    for (pi = 1; pi <= plan.n_pairs; pi++) {
        s_cal = plan.pair[pi, 1]
        t_cal = plan.pair[pi, 2]
        pairs_in[pi].s_cal = s_cal
        pairs_in[pi].t_cal = t_cal
        pairs_in[pi].s_val = plan.calendar[s_cal]
        pairs_in[pi].t_val = plan.calendar[t_cal]
    }

    first_rep = J(plan.n_consumers, 1, 0)
    for (i = 1; i <= R_raw; i++) {
        if (plan.raw_supported[i]) {
            c = plan.raw_consumer[i]
            if (c >= 1 & c <= plan.n_consumers) {
                if (first_rep[c] == 0) first_rep[c] = i
            }
        }
    }
    link_specs = cdid_selected_link_spec(plan.n_consumers, 1)
    for (c = 1; c <= plan.n_consumers; c++) {
        rep_i = first_rep[c]
        if (rep_i == 0) {
            _cdid_sel_pub_fail("cdid_selected_public_error",
                "Consumer index has no representative raw descriptor.", 498)
        }
        pi = plan.consumer[c, 2]
        s_cal = plan.pair[pi, 1]
        t_cal = plan.pair[pi, 2]
        k_ctrl = plan.raw_K_index[plan.raw_K_ptr[rep_i]..(plan.raw_K_ptr[rep_i + 1] - 1)]
        link_specs[c].link_id = (rep_i <= plan.M ?
            "cp04:catalog:" + _cdid_rb_key(rep_i) :
            "cp04:auxiliary:" + _cdid_rb_key(rep_i - plan.M))
        link_specs[c].functional_key = _cdid_sel_pub_func_key(plan, rep_i)
        link_specs[c].pair_id = pi
        link_specs[c].s_cal = s_cal
        link_specs[c].t_cal = t_cal
        link_specs[c].target_cohort_pos = plan.consumer[c, 1]
        link_specs[c].control_cohort_pos = k_ctrl
        link_specs[c].recipe_id = prof.consumer_recipe_id[c]
    }

    ctl.est_method = (B.prepared.method == "reg" ? "or" : B.prepared.method)
    ctl.outcome_pooling = ((B.prepared.method == "reg" | B.prepared.method == "dr") ?
        (B.prepared.outcome_pooling == "" ? "separate" : B.prepared.outcome_pooling) : "separate")
    ctl.selection = B.prepared.selection_canonical
    ctl.selection_model = B.prepared.selection_model
    ctl.selection_pooling = B.prepared.selection_pooling
    ctl.or_intercept = B.prepared.x_has_intercept
    ctl.eps = ctrl.baseline.glm_epsilon
    ctl.maxit = ctrl.baseline.glm_maxit
    ctl.rank_tol = ctrl.baseline.rank_tol
    /* When compute_inffunc == 0 (noinffunc), optimal/two_step GMM still requires
       core-scale moment IFs to form the optimal weighting matrix. */
    ctl.points_only = (ctrl.inference.compute_inffunc == 0 & B.settings.gmm_weighting == "identity")
    ctl.max_chunk_rows = ctrl.chunk
    ctl.max_work = (ctrl.workcap < . ? ctrl.workcap : 0)
    ctl.max_workspace_doubles = ctrl.entries

    prov.model = cdid_selected_prepare_registry(
        prof.recipes,
        prof.registry,
        pairs_in,
        cohort_idx,
        plan.cohorts,
        plan.calendar,
        B.prepared.population.A,
        B.prepared.weight.owner_weight,
        link_specs,
        ctl,
        &prov.population)

    prov.points = prov.model.estimates
    if (plan.M > 0) {
        if (missing(prov.points[plan.post_consumer])) {
            _cdid_sel_pub_fail("cdid_selected_public_error",
                "Selected provider produced nonfinite post catalog point estimate.", 459)
        }
        prov.post_ids_sorted = uniqrows(sort(plan.post_consumer, 1))
        prov.post_map_sorted = _cdid_sel_link_match_sorted(plan.post_consumer, prov.post_ids_sorted)
        prov.post_cache_ready = 1
    }
    else {
        prov.post_cache_ready = 1
    }
    for (i = 1; i <= plan.Q; i++) {
        if (prov.plan.aux_supported[i]) {
            if (prov.points[prov.plan.aux_consumer[i]] >= .) {
                prov.plan.aux_supported[i] = 0
                prov.plan.aux_reason[i] = "nonfinite_point"
            }
        }
    }
    prov.workspace_entries = _cdid_base_cells(1, prov.model.workspace_doubles, ctrl.entries, plan.workspace_entries)
    return(prov)
}

/* Callback for cdid_gmm_fit_unit_norm(): streams rows lo..hi on core scale
   (scale_core = 1, psi = (n_A/N)*phi_F) with 1/n_A normalizer. */
real matrix cdid_selected_post_core(pointer() scalar context, real scalar lo, real scalar hi)
{
    pointer(struct cdid_selected_provider scalar) scalar f
    struct cdid_selected_chunk_out scalar ch
    real colvector ids, map

    if (context == NULL) {
        _cdid_sel_pub_fail("cdid_selected_public_error", "Missing prepared selected provider context.", 498)
    }
    f = context
    if ((*f).plan.M == 0) return(J(hi - lo + 1, 0, 0))
    (*f).model.population_ref = &(*f).population
    if ((*f).cache_eval_active & lo == 1 & hi == (*f).n_F & rows((*f).eval_post_map) == (*f).plan.M) {
        ch = cdid_selected_chunk((*f).model, (lo::hi), (*f).eval_ids_sorted, 0, 0)
        (*f).cached_eval_phi = ch.phi
        return((*f).model.lambda :* ch.phi[, (*f).eval_post_map'])
    }
    if ((*f).post_cache_ready & rows((*f).post_map_sorted) == rows((*f).plan.post_consumer)) {
        if (all((*f).post_ids_sorted[(*f).post_map_sorted] :== (*f).plan.post_consumer)) {
            ch = cdid_selected_chunk((*f).model, (lo::hi), (*f).post_ids_sorted, 1, 0)
            return(ch.phi[, (*f).post_map_sorted'])
        }
    }
    ids = uniqrows(sort((*f).plan.post_consumer, 1))
    map = _cdid_sel_link_match_sorted((*f).plan.post_consumer, ids)
    ch = cdid_selected_chunk((*f).model, (lo::hi), ids, 1, 0)
    return(ch.phi[, map'])
}

/* Stream original-F frame-scale influence functions (scale_core = 0) from the
   prepared selected provider into the common report state. */
struct cdid_result_state scalar cdid_selected_provider_report(
    struct cdid_public_bundle scalar B,
    struct cdid_selected_provider scalar fit,
    real colvector post_point,
    real colvector post_available,
    string colvector post_status,
    string colvector post_reason,
    real colvector post_est_ids,
    real matrix post_loadings,
    real scalar compute_if,
    struct cdid_state_limits scalar limits,
    real scalar chunk,
    real scalar cap,
    | real matrix retained_if)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    struct cdid_selected_chunk_out scalar ch
    real colvector aux_point, aux_available, active, qids, postpos, auxpos, consumers, hit, uid
    real colvector aux_in_post, rem_consumers, rem_pos
    real matrix values, postIF, auxIF, ret_if
    real scalar Q, M, Jt, i, j, lo, hi, effective, perrow, retained, p_x, p_or, p_sel, use_ret, scale_f

    Q = rows(B.request.aux)
    M = rows(fit.plan.post_consumer)
    Jt = rows(post_est_ids)
    ret_if = (args() >= 13 ? retained_if : J(0, 0, .))
    if (rows(post_loadings) != Jt | cols(post_loadings) != M | rows(fit.plan.aux_consumer) != Q |
        chunk < 1 | chunk >= . | chunk != floor(chunk) | cap < 1 | cap >= . | cap != floor(cap)) {
        _cdid_pr_fail("cdid_native_reporting_input", "Provider/post loading identities or replay controls do not align.")
    }
    if (rows(post_est_ids)) {
        _cdid_rb_ids(post_est_ids, rows(B.request.post))
        if (any(post_available[post_est_ids] :!= 1) | rows(post_est_ids) != sum(post_available)) {
            _cdid_pr_fail("cdid_native_reporting_input", "Estimable post IDs must equal the full supported post set.")
        }
    }
    else if (sum(post_available)) {
        _cdid_pr_fail("cdid_native_reporting_input", "Supported post set has no passed target loading.")
    }
    aux_point = J(Q, 1, .)
    aux_available = fit.plan.aux_supported
    active = _cdid_rb_index(aux_available)
    if (rows(active)) aux_point[active] = fit.points[fit.plan.aux_consumer[active]]
    buffer = cdid_public_report_start(B, post_point, post_available, post_status, post_reason,
        aux_point, aux_available, fit.plan.aux_reason, compute_if, limits)

    if (compute_if) {
        fit.model.population_ref = &fit.population
        qids = J(0, 1, .)
        for (i = 1; i <= rows(B.request.report); i++) {
            if (B.request.report[i, 5] == 2 & buffer.payload.available[i]) {
                qids = qids \ B.request.report[i, 7]
            }
        }
        if (rows(qids)) qids = uniqrows(sort(qids, 1))
        consumers = (Jt ? fit.plan.post_consumer : J(0, 1, .))
        if (rows(qids)) consumers = consumers \ fit.plan.aux_consumer[qids]
        if (rows(consumers)) consumers = uniqrows(sort(consumers, 1))
        postpos = J((Jt ? M : 0), 1, .)
        for (j = 1; j <= rows(postpos); j++) {
            hit = _cdid_rb_index(consumers :== fit.plan.post_consumer[j])
            postpos[j] = hit[1]
        }
        auxpos = J(rows(qids), 1, .)
        for (j = 1; j <= rows(qids); j++) {
            hit = _cdid_rb_index(consumers :== fit.plan.aux_consumer[qids[j]])
            auxpos[j] = hit[1]
        }

        p_x = cols(B.prepared.X)
        p_or = cols(fit.profiles.registry.or_data[1].value)
        p_sel = cols(B.prepared.Z_sel)
        retained = cdid_state_cost(buffer.payload)[1] + fit.workspace_entries + length(post_loadings)
        perrow = 10 * (rows(consumers) + p_x + p_or + p_sel + 4) + 3 * (Jt + rows(qids) + rows(B.request.report))
        if (retained >= cap | perrow > cap - retained) {
            _cdid_rs_fail("Reported IF/provider replay exceeds explicit numeric workspace.", 908)
        }
        effective = min((chunk, fit.chunk_rows, floor((cap - retained) / perrow)))
        if (fit.cache_eval_active & rows(fit.cached_eval_phi) == fit.n_F &
            cols(fit.cached_eval_phi) == rows(consumers) & fit.eval_ids_sorted == consumers) {
            for (lo = 1; lo <= fit.n_F; lo = lo + effective) {
                hi = min((fit.n_F, lo + effective - 1))
                uid = (lo::hi)
                values = fit.cached_eval_phi[uid, ]
                postIF = (Jt ? values[, postpos'] * post_loadings' : J(hi - lo + 1, 0, .))
                auxIF = (rows(qids) ? values[, auxpos'] : J(hi - lo + 1, 0, .))
                cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
            }
        }
        else {
            use_ret = (rows(ret_if) == fit.n_F & cols(ret_if) == M & fit.model.lambda > 0 & fit.model.lambda < .)
            if (use_ret) {
                scale_f = 1 / fit.model.lambda
                aux_in_post = J(rows(qids), 1, 0)
                rem_consumers = J(0, 1, .)
                rem_pos = J(rows(qids), 1, 0)
                for (j = 1; j <= rows(qids); j++) {
                    hit = _cdid_rb_index(fit.plan.post_consumer :== fit.plan.aux_consumer[qids[j]])
                    if (rows(hit)) aux_in_post[j] = hit[1]
                    else rem_consumers = rem_consumers \ fit.plan.aux_consumer[qids[j]]
                }
                if (rows(rem_consumers)) {
                    rem_consumers = uniqrows(sort(rem_consumers, 1))
                    for (j = 1; j <= rows(qids); j++) if (!aux_in_post[j]) {
                        hit = _cdid_rb_index(rem_consumers :== fit.plan.aux_consumer[qids[j]])
                        rem_pos[j] = hit[1]
                    }
                }
                for (lo = 1; lo <= fit.n_F; lo = lo + effective) {
                    hi = min((fit.n_F, lo + effective - 1))
                    uid = (lo::hi)
                    postIF = (Jt ? (scale_f :* ret_if[uid, ]) * post_loadings' : J(hi - lo + 1, 0, .))
                    if (rows(qids)) {
                        auxIF = J(hi - lo + 1, rows(qids), .)
                        if (rows(rem_consumers) > 0) {
                            ch = cdid_selected_chunk(fit.model, uid, rem_consumers, 0, 0)
                            values = ch.phi
                        }
                        else values = J(hi - lo + 1, 0, .)
                        for (j = 1; j <= rows(qids); j++) {
                            if (aux_in_post[j]) auxIF[, j] = scale_f :* ret_if[uid, aux_in_post[j]]
                            else auxIF[, j] = values[, rem_pos[j]]
                        }
                    }
                    else auxIF = J(hi - lo + 1, 0, .)
                    cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
                }
            }
            else {
                for (lo = 1; lo <= fit.n_F; lo = lo + effective) {
                    hi = min((fit.n_F, lo + effective - 1))
                    uid = (lo::hi)
                    if (rows(consumers) > 0) {
                        ch = cdid_selected_chunk(fit.model, uid, consumers, 0, 0)
                        values = ch.phi
                    }
                    else {
                        values = J(hi - lo + 1, 0, .)
                    }
                    postIF = (Jt ? values[, postpos'] * post_loadings' : J(hi - lo + 1, 0, .))
                    auxIF = (rows(qids) ? values[, auxpos'] : J(hi - lo + 1, 0, .))
                    cdid_public_report_add(buffer, uid, post_est_ids, postIF, qids, auxIF)
                }
            }
        }
    }

    out = cdid_public_report_finish(buffer, limits)
    out.numeric_settings = out.numeric_settings \
        fit.weight_rows \ fit.normalization_count \
        fit.model.n_selection_blocks \ fit.model.n_selection_blocks \
        fit.model.n_treatment_fit_calls \ fit.model.n_or_block_fits \
        sum(B.prepared.U) \ cols(B.prepared.Z_sel) \
        B.prepared.sel_has_intercept \ B.prepared.selection_formula_inherited \
        B.prepared.sel_has_rho \ M \ Q
    out.numeric_setting_names = out.numeric_setting_names \
        "GMM_replay_rows" \ "GMM_normalization_count" \
        "n_selection_fits" \ "n_selection_blocks" \
        "n_treatment_fits" \ "n_pooled_OR_fits" \
        "N_U" \ "p_selection" \
        "sel_has_intercept" \ "selection_formula_inherited" \
        "sel_has_rho" \ "M_full" \ "Q_aux_union"
    out.text_settings = out.text_settings \
        ("GMM_weight_row_domain", fit.weight_row_domain) \
        ("provider_status", fit.status) \
        ("profile_kind", fit.profiles.profile_kind)
    out.source_versions = out.source_versions \
        ("estimator", "observation-weighted post-treatment and auxiliary comparisons; influence functions over all units")
    cdid_state_validate(out, limits)
    return(out)
}

/* Full selected public fit: builds comparison catalog, post system, static
   profile registry, selected provider, GMM solution, and original-F report state. */
struct cdid_public_fit_result scalar cdid_selected_public_fit(
    struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_public_binding_control scalar ctrl_fit
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_selected_profile_bundle scalar prof
    struct cdid_selected_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_gmm_result scalar sol
    real matrix G, loadings, postmeta, Vw, Lx
    real colvector point, diagnostic, reportpost, qids, eval_ids, kept, rk, ap, gq
    string colvector diagnostic_names
    real scalar M, Jt, used, remaining, i, j, N, p, def_chunk, eff_chunk

    if (B.prepared.design != "static_numeric" |
        (B.prepared.method != "ipw" & B.prepared.method != "reg" & B.prepared.method != "dr")) {
        _cdid_pr_fail("cdid_native_private_binding_scope",
            "Selected public binding requires static numeric ipw, reg, or dr.")
    }
    if (B.settings.gmm_weighting != "optimal" & B.settings.gmm_weighting != "two_step" &
        B.settings.gmm_weighting != "identity" & B.settings.gmm_weighting != "within") {
        _cdid_pr_fail("cdid_chained_public_control_error", "gmmweighting() must be optimal, two_step, identity, or within.")
    }

    B.settings.alpha = ctrl.inference.alpha
    B.settings.compute_if = ctrl.inference.compute_inffunc
    B.settings.bstrap = ctrl.inference.bstrap
    B.settings.cband = ctrl.inference.cband
    B.settings.biters = ctrl.inference.B

    if (!_cdid_sel_try_cached_plan(B, ctrl, pairs, C, plan)) {
        co = cdid_catalog_default_opts()
        co.mem_limit_bytes = ctrl.workspace_bytes
        co.work_cap = ctrl.workcap
        co.conventions = "D10_D13_admissible_frame"
        pairs = cdid_cipw_pairs_opts(B.prepared.population, B.request.cohorts, co)
        C = cdid_catalog_build_opts(B.request.calendar, B.request.cohorts, B.request.target_cohorts,
            B.request.anticipation, pairs.occupancy, B.request.control_group, B.request.post, co)
        plan = cdid_baseline_plan(B.prepared.population, C, B.request.aux,
            B.request.aux_K_ptr, B.request.aux_K_index, ctrl.entries, ctrl.workcap)
    }
    D = cdid_chained_ipw_build_design(C, B.request.post, B.request.pt_assumption, B.request.pt_start)
    B.request = cdid_public_request_tags(B.request, C)

    N = B.prepared.population.N
    p = max((1, cols(B.prepared.X), cols(B.prepared.Z_sel)))
    def_chunk = max((1, min((N, 16384, floor(ctrl.entries / max((1, 32 * (ctrl.links_bound + p + 1))))))))
    eff_chunk = ctrl.chunk
    if (ctrl.workcap >= . & ctrl.chunk == def_chunk &
        ctrl.entries >= 32 * (plan.M + p + 1) * min((N, 16384))) {
        eff_chunk = max((ctrl.chunk, min((N, 16384, floor(ctrl.entries / max((1, 32 * (plan.M + p + 1))))))))
    }
    ctrl_fit = ctrl
    ctrl_fit.chunk = eff_chunk

    prof = cdid_selected_static_profiles(B, plan)
    fit = cdid_selected_provider_prepare(B, plan, prof, ctrl_fit)
    fit.model.population_ref = &fit.population

    sys = cdid_baseline_make_system(C, D, B.request.post, ctrl.entries, fit.workspace_entries)
    support = cdid_public_post_ledger(C, D, B.request, sys)
    M = sys.M
    Jt = rows(B.request.post)
    point = J(Jt, 1, .)
    loadings = J(0, M, .)
    diagnostic = J(0, 1, .)
    diagnostic_names = J(0, 1, "")

    if (M > 0 & sys.p_retained > 0) {
        used = fit.workspace_entries + sys.workspace_entries
        used = _cdid_base_cells(8 * ctrl.chunk,
            M + cols(B.prepared.X) + cols(B.prepared.Z_sel) + 4, ctrl.entries, used)
        if (B.settings.gmm_weighting == "identity") {
            used = _cdid_gmm_identity_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries
            sol = cdid_gmm_identity_map(fit.points[plan.post_consumer], sys.full_A, sys.H_est, remaining)
        }
        else {
            used = _cdid_base_cells(2 * ctrl.chunk, M, ctrl.entries, used)
            used = _cdid_base_cells(3 * M, M, ctrl.entries, used)
            used = _cdid_gmm_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit.workspace_entries - sys.workspace_entries - M * M
            if (ctrl.inference.compute_inffunc & ctrl.workcap >= . & sys.n_estimable > 0 & eff_chunk >= fit.n_F) {
                qids = J(0, 1, .)
                for (i = 1; i <= rows(B.request.report); i++) {
                    if (B.request.report[i, 5] == 2 & B.request.report[i, 9] == 0) {
                        if (fit.plan.aux_supported[B.request.report[i, 7]]) {
                            qids = qids \ B.request.report[i, 7]
                        }
                    }
                }
                if (rows(qids)) {
                    qids = uniqrows(sort(qids, 1))
                    eval_ids = uniqrows(sort(plan.post_consumer \ fit.plan.aux_consumer[qids], 1))
                    if (rows(eval_ids) > rows(fit.post_ids_sorted) &
                        used + fit.n_F * rows(eval_ids) < remaining) {
                        fit.eval_ids_sorted = eval_ids
                        fit.eval_post_map = _cdid_sel_link_match_sorted(plan.post_consumer, eval_ids)
                        fit.cached_eval_phi = J(0, 0, .)
                        fit.cache_eval_active = 1
                    }
                }
            }
            if (B.settings.gmm_weighting == "within") {
                rk = fit.points[plan.post_consumer]
                cdid_wsg_within(B, plan, C, rk, 1, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit,
                    eff_chunk, kept, G, Vw)
                sol = cdid_gmm_fit_within_norm(rk[kept], sys.full_A[kept, ], Vw, G, M, sys.H_est,
                    ctrl.condition_trigger, ctrl.ridge_relative, fit.weight_rows, fit.normalization_count,
                    eff_chunk, &cdid_selected_post_core(), &fit, remaining, M, kept')
                Lx = J(rows(sol.target_loadings), M, 0)
                Lx[., kept'] = sol.target_loadings
                sol.target_loadings = Lx
                Lx = Vw = J(0, 0, .)
            }
            else {
                kept = (1::M)
                G = I(M)
                sol = cdid_gmm_fit_unit_norm(fit.points[plan.post_consumer], sys.full_A, G, M, sys.H_est,
                    B.settings.gmm_weighting, ctrl.condition_trigger, ctrl.ridge_relative,
                    fit.weight_rows, fit.normalization_count, eff_chunk,
                    &cdid_selected_post_core(), &fit, remaining)
            }
        }
        if (B.settings.gmm_weighting == "identity") kept = (1::M)
        if (sys.n_estimable) point[sys.estimable_ids] = sol.targets
        loadings = sol.target_loadings
        ap = J(rows(B.request.aux), 1, .)
        gq = selectindex(fit.plan.aux_supported)
        if (rows(gq)) ap[gq] = fit.points[fit.plan.aux_consumer[gq]]
        rk = fit.points[plan.post_consumer]
        cdid_chain_sp_record(B, C, D, kept, rk, loadings, sys.estimable_ids, ap, fit.plan.aux_supported)
        diagnostic = (sol.criterion \ sol.lambda \ sol.rho \ sol.trace_B \
            sol.covariance_rank \ sol.source_n_units \ sol.normalization_factor \ sol.sample_rank_upper)
        diagnostic_names = ("GMM_criterion" \ "GMM_lambda" \ "GMM_ridge_rho" \ "GMM_trace_B" \
            "GMM_covariance_rank" \ "GMM_source_replay_rows" \ "GMM_gram_normalization_factor" \ "GMM_sample_rank_upper")
        G = J(0, 0, .)
        sol.unit_gram_covariance = sol.coefficient_loadings = sol.target_loadings = J(0, 0, .)
        sol.coefficient_covariance = sol.covariance = J(0, 0, .)
    }
    else if (sys.n_estimable) {
        _cdid_pr_fail("cdid_native_reporting_input", "An estimable post target has no GMM parameter system.")
    }

    if (cdid_nuis_requested()) cdid_nuis_capture_selected(fit.model, B, plan)
    sys.full_A = sys.H_est = J(0, 0, .)
    out.payload = cdid_selected_provider_report(B, fit, point, sys.post_available,
        support.status, support.reason, sys.estimable_ids, loadings,
        ctrl.inference.compute_inffunc, ctrl.inference.limits, eff_chunk, ctrl.entries, sol.retained_if)
    fit.cached_eval_phi = J(0, 0, .)
    sol.retained_if = J(0, 0, .)
    postmeta = J(rows(B.request.report), cols(support.ledger), .)
    reportpost = _cdid_rb_index(B.request.report[, 5] :== 3)
    for (i = 1; i <= rows(reportpost); i++) {
        j = reportpost[i]
        postmeta[j, ] = support.ledger[B.request.report[j, 6], ]
    }
    out.payload.ledger = out.payload.ledger, postmeta
    out.payload.ledger_columns = out.payload.ledger_columns, "post_support." :+ support.columns
    out.payload.numeric_settings = out.payload.numeric_settings \ diagnostic \
        ctrl.workspace_bytes \ ctrl.entries \ ctrl.chunk \ ctrl.workcap \
        ctrl.condition_trigger \ ctrl.ridge_relative \ ctrl.links_bound \ ctrl.J_bound
    out.payload.numeric_setting_names = out.payload.numeric_setting_names \ diagnostic_names \
        "numeric_workspace_bytes" \ "numeric_workspace_entries" \ "chunk_rows" \
        "counted_work_limit" \ "condition_trigger" \ "ridge_relative" \
        "catalogue_links_bound" \ "inference_targets_bound"
    out.payload.text_settings = out.payload.text_settings \
        ("GMM_metric", "full original catalogue identity metric before covariance weighting") \
        ("registry_limit_meaning", "2^53-1 exact-count ceiling; not a session RAM/work budget")
    out.payload.source_versions = out.payload.source_versions \
        ("gmm", "GMM combination, normalized by the units with two observed outcomes; observation-weighted comparisons")
    out.inference = ctrl.inference
    cdid_state_validate(out.payload, ctrl.inference.limits)
    return(out)
}

end
