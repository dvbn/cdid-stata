*! Persistent nested owner record and full-union preparation for conditional imputation/augmented pretests.
*! Pinned R cdid 1.1.0.0 (commit ):
*! R/chained_imputation_public.R, R/chained_augmented_public.R, and R/chained_conditional_extract.R.
*! Retains live Stage-1 model, treatment/dated/selection readers, and pure/augmented backends
*! in a single external heap owner across CPTI preparation and piece-streaming replay.
*! Load after cdid_cpt_imputation, cdid_pretest_result, cdid_pretest_grid, cdid_imputation_public,
*! and cdid_augmented_public, and before cdid_pretest_chained.
version 19.5
mata:
mata set matastrict on

struct cdid_pt_imp_owner {
    string scalar route, design
    struct cdid_imputation_model scalar model
    real matrix Xtrt
    struct cdid_imp_dated_registry scalar registry
    struct cdid_imp_recipe_spec colvector recipes
    real matrix used
    struct cdid_imp_dated_reader scalar reader
    struct cdid_imp_backend scalar imp_backend
    struct cdid_sel_reader scalar qreader
    struct cdid_aug_backend scalar aug_backend
    real colvector design_index
    real scalar workspace_entries
}

/* Attach a canonical 1-slot static reference design record when static numeric
   intake (_cdid_request_prepare) left B.reference unpopulated. This allows
   _cdid_pretest_grid to evaluate the full-F outcome-X grid through the single
   shared fixed FV/numeric prediction path without mutating preparers. */
void _cdid_pt_ensure_static_ref(struct cdid_public_bundle scalar B,
    struct cdid_reference_inputs scalar raw)
{
    struct cdid_reference_inputs scalar S
    struct cdid_ref_design scalar D
    if (B.prepared.design == "reference_recipes") return
    if (B.reference.n_slots > 0 & B.reference.n_slots < .) return
    D.intercept = B.prepared.x_has_intercept
    D.dated = 0
    D.encoded = 1
    D.raw_names = raw.raw_names
    D.raw_index = (cols(raw.raw_names) ? (1..cols(raw.raw_names)) : J(1, 0, 0))
    D.expanded = raw.raw_names
    D.formula = invtokens(raw.raw_names)
    D.columns = (D.intercept ? ("_cons", raw.raw_names) : raw.raw_names)
    D.omitted = J(1, cols(D.expanded), 0)
    D.column_map = (cols(D.expanded) ? (D.intercept :+ (1..cols(D.expanded))) : J(1, 0, 0))
    D.keys = _cdid_ref_ix(B.prepared.required_X)
    D.value = B.prepared.X[D.keys, ]
    S = raw
    S.roles = cdid_ref_role(3, 1)
    S.roles[1].active = 1
    S.roles[1].intercept = D.intercept
    S.roles[1].dated = 0
    S.roles[1].slot = 1
    S.roles[1].formula = D.formula
    S.roles[1].raw_index = D.raw_index
    S.roles[1].keys = D.keys
    S.roles[1].first_view = J(rows(D.keys), 1, 1)
    S.roles[1].fit_keys = D.keys
    S.designs = D
    S.n_slots = 1
    S.retained_entries = rows(D.value) * cols(D.value)
    S.population = NULL
    B.reference = S
}

/* Fallback unavailable CPTI record when the complete M+Q union has zero
   supported consumer links. */
struct cdid_cpti scalar _cdid_pt_cpti_empty(struct cdid_public_bundle scalar B,
    string scalar kind, string scalar layout, real scalar workcap)
{
    struct cdid_cpti scalar ip
    struct cdid_cpti_cell_meta scalar cm
    real colvector sel_rows
    real scalar fam, n_cells, i, r
    ip.kind = kind
    ip.layout = layout
    ip.method = B.prepared.method
    ip.base_period = B.request.base_period
    ip.n_F = B.prepared.population.N
    ip.n_A = B.prepared.population.n_A
    ip.G = rows(B.request.cohorts)
    ip.n_recipes = 0
    ip.n_links = 0
    ip.n_primitive = 0
    ip.n_post_targets = rows(B.request.post)
    ip.n_aux = rows(B.request.aux)
    ip.chunk_rows = B.prepared.population.N
    ip.retained_bytes = 0
    ip.workspace_limit = 0
    ip.work_ops = 0
    ip.union_work = 0
    ip.work_limit = workcap
    ip.fingerprint = 0
    ip.ib = NULL
    ip.ab = NULL
    ip.orig_link = J(0, 1, 0)
    ip.aux_link = J(0, 1, 0)
    fam = (B.request.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_impb_ix((B.request.selection[, 1] :== fam) :& (B.request.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    ip.n_cells = n_cells
    ip.n_unavailable = n_cells
    ip.cells = cdid_cpti_cell_meta(n_cells, 1)
    if (n_cells == 0) {
        ip.status = "unavailable"
        ip.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + B.request.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        return(ip)
    }
    for (i = 1; i <= n_cells; i++) {
        r = sel_rows[i]
        cm = cdid_cpti_cell_meta()
        cm.g = B.request.selection[r, 2]
        cm.t = B.request.selection[r, 3]
        cm.from = B.request.selection[r, 5]
        cm.q = B.request.selection[r, 6]
        cm.sign = B.request.selection[r, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.link = 0
        cm.recipe = 0
        cm.supported = 0
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = 0
        cm.route = layout + " " + kind + " " + B.prepared.method
        cm.status = "unavailable"
        cm.reason = (B.request.selection[r, 8] != 0 ? B.request.selection_reason[r] : "Unsupported pre-treatment comparison.")
        ip.cells[i] = cm
    }
    ip.status = "unavailable"
    ip.reason = "The complete requested test is unavailable: " + _cdid_cpti_istr(ip.n_unavailable) + " declared cell(s) are unavailable. No available-only subset is tested."
    return(ip)
}

/* Prepare the persistent owner in place on the caller's external heap slot (*pOwn)
   and bind the CPTI producer over the complete M+Q union. Every borrowed
   pointer (&((*pOwn).model), &((*pOwn).Xtrt), &((*pOwn).registry), &((*pOwn).reader),
   &((*pOwn).imp_backend), &((*pOwn).qreader), &((*pOwn).aug_backend)) addresses
   the live external heap record directly and is never copied after binding. */
struct cdid_cpti scalar _cdid_pt_imp_owner_prepare(
    pointer(struct cdid_pt_imp_owner scalar) scalar pOwn,
    struct cdid_public_bundle scalar B,
    pointer(struct cdid_population scalar) scalar pp,
    struct cdid_baseline_plan_state scalar plan,
    pointer(struct cdid_reference_inputs scalar) scalar pS,
    struct cdid_reference_catalog scalar catalog,
    struct cdid_public_binding_control scalar ctl,
    real scalar workcap)
{
    struct cdid_cpti scalar cpti
    struct cdid_imputation_pair rowvector imp_pairs
    struct cdid_imputation_budget scalar imp_budget
    struct cdid_sel_decl scalar qdecl
    struct cdid_imp_matrix_slot colvector qdata
    real colvector recipe_qslot, link_pair_map, design_index
    real matrix links, controls, active_requests, Z
    string colvector link_ids
    string rowvector slope_names
    real scalar i, p, ext_entries, rec, aug
    string scalar kind, layout

    if (pOwn == NULL | pp == NULL)
        _cdid_pt_fail("Conditional imputation object requires external object and population pointers.", 498)
    *pOwn = cdid_pt_imp_owner()
    aug = (B.prepared.selection_est_method == "dr")
    kind = (aug ? "augmented" : "imputation")
    layout = (B.prepared.design == "reference_recipes" ? "dated" : "static")
    (*pOwn).route = kind
    (*pOwn).design = B.prepared.design

    if (B.prepared.design == "reference_recipes") {
        if (pS == NULL) _cdid_pt_fail("Conditional reference imputation object requires a reference design pointer.", 498)
        if (catalog.n_consumers < 1) {
            cpti = _cdid_pt_cpti_empty(B, kind, layout, workcap)
            (*pOwn).design_index = J(cpti.n_cells, 1, 0)
            (*pOwn).workspace_entries = _cdid_base_cells(1, (*pS).retained_entries, ctl.entries, plan.workspace_entries)
            return(cpti)
        }
        cdid_imp_reference_links(plan, catalog, links, controls, link_ids)
        cdid_imp_reference_registry(*pS, (*pOwn).registry, (*pOwn).recipes)
        (*pOwn).used = cdid_imp_used_from_links(*pp, (*pS).n_recipes, links, controls)
        ext_entries = plan.workspace_entries + (*pS).retained_entries + length(links) + length(controls) + length((*pOwn).used)
        for (i = 1; i <= rows((*pOwn).registry.z_data); i++) ext_entries = ext_entries + length((*pOwn).registry.z_data[i].value)
        for (i = 1; i <= rows((*pOwn).registry.x_data); i++) ext_entries = ext_entries + length((*pOwn).registry.x_data[i].value)
        for (i = 1; i <= rows((*pOwn).registry.w_data); i++) ext_entries = ext_entries + length((*pOwn).registry.w_data[i].value)
        if (8 * ext_entries >= ctl.workspace_bytes)
            _cdid_pt_fail("Reference imputation catalogue/design exceeds explicit numeric workspace.", 908)
        imp_budget = cdid_imputation_default_budget(ctl.chunk, ctl.baseline.max_cached,
            ctl.workspace_bytes - 8 * ext_entries, ctl.workcap)
        (*pOwn).reader = cdid_imp_dated_reader_make(pp, &((*pOwn).registry), (*pOwn).recipes, (*pOwn).used,
            B.prepared.method != "reg", imp_budget)
        (*pOwn).model = cdid_imp_fit_dated((*pOwn).reader, links, controls,
            B.prepared.imputation_pooling, ctl.baseline.rank_tol, imp_budget)
        if (aug) {
            qdecl = cdid_aug_route_gate(B.prepared.method, 0, B.prepared.selection,
                B.prepared.selection_universe, B.prepared.selection_model,
                B.prepared.selection_pooling, B.prepared.sel_has_rho, B.prepared.selection_justification)
            (*pOwn).imp_backend = cdid_imp_backend_prepare_d(&((*pOwn).model), &((*pOwn).reader),
                links, controls, link_ids, (B.prepared.method == "reg" ? "reg" : "ipw"),
                ctl.baseline.glm_epsilon, ctl.baseline.glm_maxit, ctl.chunk, imp_budget,
                "pretest_reference_augmented_factors", "", (B.prepared.method == "reg"))
            cdid_aug_reference_qslots(*pS, (*pOwn).recipes, qdata, recipe_qslot)
            (*pOwn).aug_backend = cdid_aug_backend_prepare_d(&((*pOwn).imp_backend), &((*pOwn).reader),
                qdata, recipe_qslot, (*pS).roles[3].intercept,
                ((*pS).roles[3].intercept ? 1 : 0), qdecl,
                ctl.chunk, imp_budget, "pretest_reference_augmented",
                B.prepared.method, B.prepared.outcome_pooling)
            qdata = cdid_imp_matrix_slot(0, 1)
            cpti = cdid_cpti_prepare_aug_dated(&((*pOwn).aug_backend), B.request,
                catalog.post_consumer, catalog.aux_consumer, workcap)
        }
        else {
            (*pOwn).imp_backend = cdid_imp_backend_prepare_d(&((*pOwn).model), &((*pOwn).reader),
                links, controls, link_ids, B.prepared.method,
                ctl.baseline.glm_epsilon, ctl.baseline.glm_maxit, ctl.chunk, imp_budget,
                "pretest_reference_imputation", B.prepared.outcome_pooling)
            cpti = cdid_cpti_prepare_dated(&((*pOwn).imp_backend), B.request,
                catalog.post_consumer, catalog.aux_consumer, workcap)
        }
        design_index = J(cpti.n_cells, 1, 0)
        for (i = 1; i <= cpti.n_cells; i++) {
            rec = cpti.cells[i].recipe
            if (cpti.cells[i].supported & rec >= 1 & rec <= rows((*pS).recipes))
                design_index[i] = (*pS).recipes[rec, 1]
        }
        (*pOwn).design_index = design_index
    }
    else {
        if (plan.n_consumers < 1) {
            cpti = _cdid_pt_cpti_empty(B, kind, layout, workcap)
            (*pOwn).design_index = J(cpti.n_cells, 1, 0)
            (*pOwn).workspace_entries = plan.workspace_entries
            return(cpti)
        }
        cdid_imputation_public_links(plan, links, controls, link_ids)
        cdid_imp_links_to_requests(B.request.calendar, links, controls, imp_pairs, active_requests, link_pair_map)
        p = cols(B.prepared.X)
        slope_names = (B.prepared.x_has_intercept ? (p > 1 ? B.prepared.x_names[2..p] : J(1, 0, "")) : B.prepared.x_names)
        Z = (B.prepared.x_has_intercept ? (p > 1 ? B.prepared.X[, 2..p] : J((*pp).N, 0, .)) : B.prepared.X)
        (*pOwn).Xtrt = (B.prepared.method == "reg" ? J((*pp).N, 0, .) : B.prepared.X)
        ext_entries = plan.workspace_entries + length(Z) + length((*pOwn).Xtrt) + length(links) +
            length(controls) + length(active_requests) + length(link_pair_map)
        if (8 * ext_entries >= ctl.workspace_bytes)
            _cdid_pt_fail("Imputation catalogue/design exceeds explicit numeric workspace.", 908)
        imp_budget = cdid_imputation_default_budget(ctl.chunk, ctl.baseline.max_cached,
            ctl.workspace_bytes - 8 * ext_entries, ctl.workcap)
        (*pOwn).model = cdid_imputation_fit_static(*pp, B.prepared.weight.owner_weight, Z,
            imp_pairs, active_requests, B.prepared.imputation_pooling, B.prepared.x_has_intercept,
            ctl.baseline.rank_tol, imp_budget, slope_names)
        if (aug) {
            qdecl = cdid_aug_route_gate(B.prepared.method, 0, B.prepared.selection,
                B.prepared.selection_universe, B.prepared.selection_model,
                B.prepared.selection_pooling, B.prepared.sel_has_rho, B.prepared.selection_justification)
            (*pOwn).imp_backend = cdid_imp_backend_prepare(&((*pOwn).model), &((*pOwn).Xtrt),
                links, controls, link_ids, (B.prepared.method == "reg" ? "reg" : "ipw"),
                ctl.baseline.glm_epsilon, ctl.baseline.glm_maxit, ctl.chunk, imp_budget,
                "pretest_static_augmented_factors", "", (B.prepared.method == "reg"))
            (*pOwn).qreader = cdid_sel_reader_pop(*pp, B.prepared.Z_sel, B.prepared.sel_has_intercept,
                (B.prepared.sel_has_intercept ? 1 : 0), B.prepared.weight.owner_weight, qdecl,
                J((*pp).N, 1, 1), ctl.chunk)
            (*pOwn).aug_backend = cdid_aug_backend_prepare(&((*pOwn).imp_backend), &((*pOwn).qreader),
                ctl.chunk, imp_budget, "pretest_static_augmented",
                B.prepared.method, B.prepared.outcome_pooling)
            cpti = cdid_cpti_prepare_aug_static(&((*pOwn).aug_backend), B.request, workcap)
        }
        else {
            (*pOwn).imp_backend = cdid_imp_backend_prepare(&((*pOwn).model), &((*pOwn).Xtrt),
                links, controls, link_ids, B.prepared.method,
                ctl.baseline.glm_epsilon, ctl.baseline.glm_maxit, ctl.chunk, imp_budget,
                "pretest_static_imputation", B.prepared.outcome_pooling)
            cpti = cdid_cpti_prepare_static(&((*pOwn).imp_backend), B.request, workcap)
        }
        (*pOwn).design_index = J(cpti.n_cells, 1, 0)
    }

    (*pOwn).workspace_entries = _cdid_base_cells(1, ext_entries + ceil((cpti.retained_bytes +
        (aug ? (*pOwn).aug_backend.scratch_bytes : (*pOwn).imp_backend.scratch_bytes)) / 8), ctl.entries, 0)
    return(cpti)
}

end
