*! Recipe-indexed (dated) observation-weighted selected provider bridge; private; StataNow 19.5.
*! GPL-3 adaptation of David Benatia, Christophe Bellego, Joel Cuerrier and Vincent
*! Dortet-Bernadet, cdid 1.1.0.0 /:
*! R/chained_reference_inputs.R (build_chained_reference_reader, chained_reference_x/_s,
*! chained_reference_weight, chained_reference_stat), R/chained_selection_models.R and
*! R/chained_selection_consumers.R (fit_selection_models, fit_cp04_selection_models and
*! fit_or_selection_models with `reference`). Conditional-MTAR observation weighting for
*! ipw/reg/dr only; not actual-change imputation or the augmented correction.
*! Load cdid_population, cdid_public_request, cdid_reference_inputs, cdid_baseline_provider,
*! cdid_reference_catalog, cdid_selection and cdid_selected_link (used_scope) first.
version 19.5
mata:

struct cdid_selref_controls {
    string scalar method, outcome_pooling, selection_model, selection_pooling, universe, justification
    real scalar has_rho, eps, maxit, rank_tol, chunk_rows, max_work, max_entries, points_only
    real colvector rho
}

/* Exact inputs handed to the cdid_selected_prepare_registry(..., used_scope=1). */
struct cdid_selref_inputs {
    struct cdid_selected_registry scalar registry
    struct cdid_selected_recipe_spec colvector recipes
    struct cdid_selected_pair_in colvector pairs
    struct cdid_selected_link_spec colvector links
    struct cdid_selected_controls scalar ctl
    real colvector cohort_idx, w_pop, recipe_reference, reference_recipe, consumer_recipe
    real colvector sel_cal, or_cal, trt_cal, w_cal
    real matrix scope, recipe_dates, recipe_keep
    string colvector recipe_keys
    real scalar planned_entries
}

/* Population, reference inputs, plan and catalogue are BORROWED: the caller keeps
   them alive and immutable for the provider's lifetime (O(1) checks on replay). */
struct cdid_selref_provider {
    string scalar method, status, weight_row_domain, profile_kind
    real scalar N, R_hist, n_F, n_A, M, Q, n_consumers, n_recipes, n_reference_recipes
    real scalar normalization_count, weight_rows, chunk_rows, workspace_entries
    pointer(struct cdid_population scalar) scalar population
    pointer(struct cdid_reference_inputs scalar) scalar reference
    pointer(struct cdid_baseline_plan_state scalar) scalar plan
    pointer(struct cdid_reference_catalog scalar) scalar catalog
    struct cdid_selected_model scalar model
    real colvector recipe_reference, reference_recipe, consumer_recipe
    real colvector sel_cal, or_cal, trt_cal, w_cal
    real matrix scope, recipe_dates
    string colvector recipe_keys
    real colvector points, post_consumer, aux_consumer, aux_supported
    string colvector aux_reason
    real colvector post_ids_sorted, post_map_sorted
    real scalar post_cache_ready
    real colvector eval_ids_sorted, eval_post_map
    real matrix cached_eval_phi
    real scalar cache_eval_active
}

void _cdid_selref_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_selref_ix(real matrix pred)
{
    if (rows(pred) * cols(pred) == 0) return(J(0, 1, 0))
    return(vec(selectindex(pred)))
}

string scalar _cdid_selref_num(real scalar x)
{
    return(x == . ? "Inf" : strtrim(strofreal(x, "%24.17g")))
}

real scalar _cdid_selref_pos(real colvector table, real scalar value)
{
    real colvector hit
    hit = _cdid_selref_ix(table :== value)
    if (rows(hit) != 1) _cdid_selref_fail("integrity", "Registry slot calendar index is absent or repeated.", 498)
    return(hit[1])
}

real scalar _cdid_selref_bit(real scalar mask, real scalar bit)
{
    return(mod(floor(mask / bit), 2) == 1)
}

struct cdid_selref_controls scalar cdid_selref_default_controls(real scalar N)
{
    struct cdid_selref_controls scalar c
    if (N < 1 | N >= . | N != floor(N)) _cdid_selref_fail("input", "N must be a positive integer.", 498)
    c.method = "dr"
    c.outcome_pooling = "separate"
    c.selection_model = "pair_logit"
    c.selection_pooling = "separate"
    c.universe = "F"
    c.justification = ""
    c.has_rho = 0
    c.eps = 1e-12
    c.maxit = 50
    c.rank_tol = 1e-7
    c.chunk_rows = 0
    c.max_work = 0
    c.max_entries = .
    c.points_only = 0
    c.rho = J(N, 1, 1)
    return(c)
}

/* One original-F slot of a reference design role at calendar index d (0 = dateless).
   Values exist ONLY at the role's actual required keys, read from its fixed union;
   every other row stays missing. drop_intercept removes the leading _cons column
   (outcome slopes; the OR intercept is the block intercept, as on the static path). */
real matrix _cdid_selref_design_slot(struct cdid_reference_inputs scalar S,
    real scalar role, real scalar d, real scalar drop_intercept)
{
    real colvector keys, units
    real matrix value, out
    real scalar slot, p0, drop, dd

    slot = S.roles[role].slot
    if (!S.roles[role].active | slot < 1 | slot > S.n_slots)
        _cdid_selref_fail("integrity", "Requested reference design role is inactive.", 498)
    if (!S.designs[slot].encoded)
        _cdid_selref_fail("integrity", "Reference design slot is not encoded.", 498)
    p0 = cols(S.designs[slot].value)
    drop = (drop_intercept & S.designs[slot].intercept)
    if (drop) {
        if (p0 < 1) _cdid_selref_fail("integrity", "Encoded intercept column is absent.", 498)
        if (S.designs[slot].columns[1] != "_cons")
            _cdid_selref_fail("integrity", "Encoded intercept is not the leading _cons column.", 498)
    }
    dd = (S.roles[role].dated ? d : 0)
    keys = S.roles[role].keys
    if (rows(keys)) keys = select(keys, floor(keys :/ (S.N + 1)) :== dd)
    out = J(S.N, p0 - drop, .)
    if (rows(keys) == 0 | cols(out) == 0) return(out)
    value = cdid_ref_design_rows(S, role, keys)
    if (drop) value = value[, 2..p0]
    units = keys :- dd * (S.N + 1)
    out[units, ] = value
    return(out)
}

/* One weight-anchor slot: anchored sampling weights at the actual required
   (anchor, unit) keys only; missing elsewhere. First-calendar V is not used. */
real colvector _cdid_selref_weight_slot(struct cdid_reference_inputs scalar S, real scalar a)
{
    real colvector pos, out
    out = J(S.N, 1, .)
    if (rows(S.weight_keys) == 0) return(out)
    pos = _cdid_selref_ix(floor(S.weight_keys :/ (S.N + 1)) :== a)
    if (rows(pos)) out[S.weight_keys[pos] :- a * (S.N + 1)] = S.weights[pos]
    return(out)
}

string scalar _cdid_selref_recipe_key(struct cdid_reference_inputs scalar S, real scalar r)
{
    real scalar d, a
    string scalar design, weights
    d = S.recipes[r, 1]
    a = S.recipes[r, 2]
    design = (d == 0 ? "unit" : "date=" + _cdid_selref_num(S.calendar[d]))
    weights = (a == 0 ? "dateless" : S.fix_weights + "@" + _cdid_selref_num(S.calendar[a]))
    if (cols(S.recipes) >= 3) if (S.recipes[r, 3] > 0)
        weights = weights + "|eligibility:b=" + _cdid_selref_num(S.calendar[S.recipes[r, 3]])
    return("formula:" + S.roles[1].formula + "|selection_formula:" + S.roles[3].formula +
        "|design:" + design + "|weights:" + weights)
}

/* Same comparison-catalogue functional key as the static selected provider. */
string scalar _cdid_selref_func_key(struct cdid_baseline_plan_state scalar p, real scalar i)
{
    real colvector K
    real scalar j
    string scalar k_str, prefix
    K = p.raw_K_index[|p.raw_K_ptr[i] \ p.raw_K_ptr[i + 1] - 1|]
    k_str = ""
    for (j = 1; j <= rows(K); j++) k_str = k_str + (j > 1 ? "," : "") + _cdid_selref_num(p.cohorts[K[j]])
    prefix = (i <= p.M ? "cmpcat:CP04-post:" : "cmpcat:CP04-aux:")
    return(prefix + "g=" + _cdid_selref_num(p.cohorts[p.raw[i, 1]]) +
        "|s=" + _cdid_selref_num(p.calendar[p.raw[i, 2]]) +
        "|t=" + _cdid_selref_num(p.calendar[p.raw[i, 3]]) + "|K=" + k_str)
}

real colvector _cdid_selref_K(struct cdid_reference_catalog scalar cat, real scalar c)
{
    if (cat.K_ptr[c + 1] <= cat.K_ptr[c])
        _cdid_selref_fail("integrity", "A catalogue consumer has no ordered control set.", 498)
    return(cat.K_index[|cat.K_ptr[c] \ cat.K_ptr[c + 1] - 1|])
}

/* Builds unique original-F registry slots per (role, design calendar index) and
   per weight anchor, model recipes = reference recipes actually used by
   catalogue consumers (key identity kept), one link per catalogue consumer.
   Validation reads only metadata and each slot's own required keys. */
struct cdid_selref_inputs scalar cdid_selref_build(
    pointer(struct cdid_population scalar) scalar population,
    pointer(struct cdid_reference_inputs scalar) scalar reference,
    pointer(struct cdid_baseline_plan_state scalar) scalar plan,
    pointer(struct cdid_reference_catalog scalar) scalar catalog,
    struct cdid_selref_controls scalar controls)
{
    struct cdid_selref_inputs scalar out
    real matrix usedco_mat
    real colvector U, hit, K, Kp, used_ref, cohorts_c
    real scalar N, J, n_raw, nc, nrec, c, i, j, k, r, m, d, a, mask, p_id, e
    real scalar d_sel, p_or, p_trt, lead, sel_dated, or_dated, trt_dated

    if (population == NULL | reference == NULL | plan == NULL | catalog == NULL)
        _cdid_selref_fail("owner", "Population, reference-input, full-plan and catalogue objects are required.", 498)
    if (controls.method != "ipw" & controls.method != "reg" & controls.method != "dr")
        _cdid_selref_fail("input", "method must be ipw, reg or dr.", 498)
    if (controls.outcome_pooling != "" & controls.outcome_pooling != "separate" & controls.outcome_pooling != "pooled")
        _cdid_selref_fail("input", "outcome_pooling must be separate or pooled.", 498)
    if (controls.points_only != 0 & controls.points_only != 1)
        _cdid_selref_fail("input", "points_only must be 0 or 1.", 498)
    if (controls.chunk_rows < 0 | controls.chunk_rows >= . | controls.chunk_rows != floor(controls.chunk_rows))
        _cdid_selref_fail("input", "chunk_rows must be 0 (all rows) or a positive integer.", 498)
    if (controls.max_entries < .) {
        if (controls.max_entries < 1 | controls.max_entries != floor(controls.max_entries))
            _cdid_selref_fail("input", "max_entries must be . or a positive integer.", 498)
    }

    /* One borrowed population owner; roster/calendar checked once here. */
    if ((*reference).population != population)
        _cdid_selref_fail("owner", "Reference inputs must borrow the population object supplied here.", 498)
    cdid_ref_activate(*reference)
    N = (*population).N
    if ((*reference).route != "selected" | (*reference).method != controls.method)
        _cdid_selref_fail("input", "Reference inputs were not prepared for route(selected) and this method.", 498)
    if (!(*reference).roles[3].active)
        _cdid_selref_fail("input", "The selected route requires an active selection design role.", 498)
    if (controls.method != "reg" & !(*reference).roles[2].active)
        _cdid_selref_fail("input", "ipw and dr require an active treatment design role.", 498)
    for (k = 1; k <= (*reference).n_slots; k++) {
        if (!(*reference).designs[k].encoded)
            _cdid_selref_fail("input", "Every reference design slot must be encoded before consumption.", 498)
    }

    if (controls.universe == "F") U = J(N, 1, 1)
    else if (controls.universe == "B") U = (*population).B
    else if (controls.universe == "A") U = (*population).A
    else _cdid_selref_fail("input", "Selection universe must be F, B or A.", 498)
    if (rows((*reference).unit_universe) != N)
        _cdid_selref_fail("input", "Reference selection roster has the wrong length.", 498)
    if (any((*reference).unit_universe :!= U))
        _cdid_selref_fail("input", "Reference selection roster differs from the declared F/B/A universe.", 498)
    if (rows(controls.rho) != N | cols(controls.rho) != 1)
        _cdid_selref_fail("input", "Known transport rho must be an original-F N x 1 vector.", 498)

    /* Same cohort rule and dictionary positions as the static selected provider. */
    if (any((*population).cohort_known :== 0))
        _cdid_selref_fail("cdid_chained_request_cohort_error",
            "Corollary1 requires known cohorts on every F row, including frame-only rows.", 459)
    J = rows((*plan).cohorts)
    if (J > 1) {
        if (rows(uniqrows(sort((*plan).cohorts, 1))) != J)
            _cdid_selref_fail("cdid_chained_request_cohort_error", "Known unit cohort is absent from catalogue cohort dictionary.", 498)
    }
    out.cohort_idx = J(N, 1, 0)
    for (k = 1; k <= J; k++) {
        hit = _cdid_selref_ix((*population).cohort :== (*plan).cohorts[k])
        if (rows(hit) > 0) out.cohort_idx[hit] = J(rows(hit), 1, k)
    }
    if (any(out.cohort_idx :== 0))
        _cdid_selref_fail("cdid_chained_request_cohort_error", "Known unit cohort is absent from catalogue cohort dictionary.", 498)

    n_raw = rows((*plan).raw)
    nc = (*catalog).n_consumers
    if ((*catalog).M != (*plan).M | (*catalog).Q != (*plan).Q | (*catalog).n_raw != n_raw |
        n_raw != (*plan).M + (*plan).Q | (*reference).n_views != n_raw |
        rows((*reference).view_recipe) != n_raw | nc < 1 | rows((*catalog).consumers) != nc |
        rows((*catalog).consumer_recipe) != nc | rows((*catalog).representative) != nc |
        rows((*catalog).post_consumer) != (*plan).M | rows((*catalog).aux_consumer) != (*plan).Q |
        rows((*catalog).K_ptr) != nc + 1)
        _cdid_selref_fail("input", "Plan, shared catalogue and reference views are not ONE original post plus both-auxiliary union.", 498)
    if (rows((*plan).calendar) != rows((*population).calendar))
        _cdid_selref_fail("input", "Plan and population calendars disagree.", 498)
    if (any((*plan).calendar :!= (*population).calendar))
        _cdid_selref_fail("input", "Plan and population calendars disagree.", 498)

    usedco_mat = J((*reference).n_recipes, J, -1)
    for (i = 1; i <= rows((*reference).used_cohorts); i++) {
        r = (*reference).used_cohorts[i, 1]
        k = (*reference).used_cohorts[i, 2]
        if (r >= 1 & r <= (*reference).n_recipes & k >= 1 & k <= J)
            usedco_mat[r, k] = (*reference).used_cohorts[i, 3]
    }
    used_ref = uniqrows(sort((*catalog).consumer_recipe, 1))
    nrec = rows(used_ref)
    if (any(used_ref :< 1) | any(used_ref :> (*reference).n_recipes) | any(used_ref :!= floor(used_ref)))
        _cdid_selref_fail("integrity", "Catalogue consumer design block is not a reference design block.", 498)
    out.recipe_reference = used_ref
    out.reference_recipe = J((*reference).n_recipes, 1, 0)
    out.reference_recipe[used_ref] = (1::nrec)
    out.consumer_recipe = out.reference_recipe[(*catalog).consumer_recipe]
    out.scope = J(nrec, J, 0)
    for (c = 1; c <= nc; c++) {
        i = (*catalog).representative[c]
        if (i < 1 | i > n_raw | i != floor(i))
            _cdid_selref_fail("integrity", "Catalogue representative is not an original raw view.", 498)
        if (!(*plan).raw_supported[i])
            _cdid_selref_fail("integrity", "Catalogue consumer represents an unsupported raw view.", 498)
        if (any((*catalog).consumers[c, 1..3] :!= (*plan).raw[i, ]))
            _cdid_selref_fail("integrity", "Catalogue consumer differs from its plan raw descriptor.", 498)
        r = (*catalog).consumer_recipe[c]
        if ((*catalog).consumers[c, 4] != r | (*reference).view_recipe[i] != r)
            _cdid_selref_fail("integrity", "Catalogue consumer design block differs from the reference view design block.", 498)
        K = _cdid_selref_K(*catalog, c)
        Kp = (*plan).raw_K_index[|(*plan).raw_K_ptr[i] \ (*plan).raw_K_ptr[i + 1] - 1|]
        if (rows(K) != rows(Kp))
            _cdid_selref_fail("integrity", "Catalogue consumer controls differ from the plan.", 498)
        if (any(K :!= Kp))
            _cdid_selref_fail("integrity", "Catalogue consumer controls differ from the plan.", 498)
        cohorts_c = (*plan).raw[i, 1] \ K
        m = out.consumer_recipe[c]
        for (k = 1; k <= rows(cohorts_c); k++) {
            mask = (cohorts_c[k] >= 1 & cohorts_c[k] <= J ? usedco_mat[r, cohorts_c[k]] : -1)
            if (mask < 0)
                _cdid_selref_fail("support", "A consumer cohort lies outside its reference design block's required supports.", 459)
            if (!_cdid_selref_bit(mask, 4) | !_cdid_selref_bit(mask, 8))
                _cdid_selref_fail("support", "Reference design block supports lack selection design or weight rows.", 459)
            if (controls.method != "ipw" & !_cdid_selref_bit(mask, 1))
                _cdid_selref_fail("support", "Reference design block supports lack outcome design rows.", 459)
            if (controls.method != "reg" & !_cdid_selref_bit(mask, 2))
                _cdid_selref_fail("support", "Reference design block supports lack treatment design rows.", 459)
            out.scope[m, cohorts_c[k]] = 1
        }
    }

    /* Restricted eligibility contexts (missingreference(exclude)): the units of
       each model recipe's context (all units when unrestricted). */
    out.recipe_keep = J(0, 0, .)
    if (cols((*reference).recipes) >= 3 & rows((*reference).elig_index)) {
        out.recipe_keep = J(N, nrec, 1)
        for (m = 1; m <= nrec; m++) {
            e = (*reference).recipes[used_ref[m], 3]
            if (e > 0) {
                hit = _cdid_selref_ix((*reference).elig_index :== e)
                if (rows(hit) != 1) _cdid_selref_fail("integrity", "Unknown reference eligibility context.", 498)
                out.recipe_keep[, m] = !(*reference).elig_excluded[, hit[1]]
            }
        }
    }
    out.recipe_dates = J(nrec, 4, .)
    out.recipe_keys = J(nrec, 1, "")
    for (m = 1; m <= nrec; m++) {
        r = used_ref[m]
        d = (*reference).recipes[r, 1]
        a = (*reference).recipes[r, 2]
        out.recipe_dates[m, 1] = d
        out.recipe_dates[m, 2] = a
        if (d > 0) out.recipe_dates[m, 3] = (*population).calendar[d]
        if (a > 0) out.recipe_dates[m, 4] = (*population).calendar[a]
        out.recipe_keys[m] = _cdid_selref_recipe_key(*reference, r)
    }

    sel_dated = (*reference).roles[3].dated
    or_dated = (*reference).roles[1].dated
    trt_dated = (controls.method != "reg" ? (*reference).roles[2].dated : 0)
    if ((sel_dated | or_dated | trt_dated) & any(out.recipe_dates[, 1] :== 0))
        _cdid_selref_fail("integrity", "A dated role has a design block without a design reference date.", 498)
    out.sel_cal = (sel_dated ? uniqrows(sort(out.recipe_dates[, 1], 1)) : 0)
    if (controls.method != "ipw") out.or_cal = (or_dated ? uniqrows(sort(out.recipe_dates[, 1], 1)) : 0)
    else out.or_cal = J(0, 1, .)
    if (controls.method != "reg") out.trt_cal = (trt_dated ? uniqrows(sort(out.recipe_dates[, 1], 1)) : 0)
    else out.trt_cal = J(0, 1, .)
    out.w_cal = uniqrows(sort(out.recipe_dates[, 2], 1))

    lead = (*reference).designs[(*reference).roles[3].slot].intercept
    d_sel = cols((*reference).designs[(*reference).roles[3].slot].value)
    p_or = 0
    if (controls.method != "ipw") {
        k = (*reference).roles[1].slot
        p_or = cols((*reference).designs[k].value) - (*reference).designs[k].intercept
    }
    p_trt = 0
    if (controls.method != "reg") p_trt = cols((*reference).designs[(*reference).roles[2].slot].value)
    /* Registry charge before any original-F slot is allocated. */
    out.planned_entries = _cdid_base_cells(N, rows(out.sel_cal) * d_sel + rows(out.or_cal) * p_or +
        rows(out.trt_cal) * p_trt + rows(out.w_cal) + 2, controls.max_entries, (*plan).workspace_entries)

    out.registry.sel_data = cdid_selected_matrix_slot(rows(out.sel_cal), 1)
    for (j = 1; j <= rows(out.sel_cal); j++) {
        out.registry.sel_data[j].calendar_index = out.sel_cal[j]
        out.registry.sel_data[j].value = _cdid_selref_design_slot(*reference, 3, out.sel_cal[j], 0)
    }
    out.registry.or_data = cdid_selected_matrix_slot(rows(out.or_cal), 1)
    for (j = 1; j <= rows(out.or_cal); j++) {
        out.registry.or_data[j].calendar_index = out.or_cal[j]
        out.registry.or_data[j].value = _cdid_selref_design_slot(*reference, 1, out.or_cal[j], 1)
    }
    out.registry.trt_data = cdid_selected_matrix_slot(rows(out.trt_cal), 1)
    for (j = 1; j <= rows(out.trt_cal); j++) {
        out.registry.trt_data[j].calendar_index = out.trt_cal[j]
        out.registry.trt_data[j].value = _cdid_selref_design_slot(*reference, 2, out.trt_cal[j], 0)
    }
    out.registry.weight_data = cdid_selected_matrix_slot(rows(out.w_cal), 1)
    for (j = 1; j <= rows(out.w_cal); j++) {
        out.registry.weight_data[j].calendar_index = out.w_cal[j]
        out.registry.weight_data[j].value = _cdid_selref_weight_slot(*reference, out.w_cal[j])
    }
    out.registry.rho_data = cdid_selected_matrix_slot(1, 1)
    out.registry.rho_data[1].calendar_index = 0
    out.registry.rho_data[1].value = controls.rho

    out.recipes = cdid_selected_recipe_spec(nrec, 1)
    for (m = 1; m <= nrec; m++) {
        d = out.recipe_dates[m, 1]
        out.recipes[m].recipe_id = m
        out.recipes[m].recipe_key = out.recipe_keys[m]
        out.recipes[m].design_cal = d
        out.recipes[m].weight_cal = out.recipe_dates[m, 2]
        out.recipes[m].sel_id = _cdid_selref_pos(out.sel_cal, (sel_dated ? d : 0))
        out.recipes[m].or_id = 0
        if (controls.method != "ipw") out.recipes[m].or_id = _cdid_selref_pos(out.or_cal, (or_dated ? d : 0))
        out.recipes[m].trt_id = 0
        if (controls.method != "reg") out.recipes[m].trt_id = _cdid_selref_pos(out.trt_cal, (trt_dated ? d : 0))
        out.recipes[m].weight_id = _cdid_selref_pos(out.w_cal, out.recipe_dates[m, 2])
        out.recipes[m].rho_id = 1
        out.recipes[m].sel_has_intercept = lead
        out.recipes[m].sel_int_col = (lead ? 1 : 0)
        out.recipes[m].has_xformla = 1
        out.recipes[m].has_rho = controls.has_rho
        out.recipes[m].universe = controls.universe
        out.recipes[m].justification = controls.justification
    }

    out.pairs = cdid_selected_pair_in((*plan).n_pairs, 1)
    for (p_id = 1; p_id <= (*plan).n_pairs; p_id++) {
        out.pairs[p_id].s_cal = (*plan).pair[p_id, 1]
        out.pairs[p_id].t_cal = (*plan).pair[p_id, 2]
        out.pairs[p_id].s_val = (*plan).calendar[(*plan).pair[p_id, 1]]
        out.pairs[p_id].t_val = (*plan).calendar[(*plan).pair[p_id, 2]]
    }

    out.links = cdid_selected_link_spec(nc, 1)
    for (c = 1; c <= nc; c++) {
        i = (*catalog).representative[c]
        m = out.consumer_recipe[c]
        p_id = _cdid_base_find2((*plan).pair, (*plan).raw[i, 2], (*plan).raw[i, 3])
        if (p_id < 1) _cdid_selref_fail("integrity", "Consumer pair is absent from the plan pair dictionary.", 498)
        if (i <= (*plan).M) out.links[c].link_id = "cp04:catalog:" + _cdid_selref_num(i)
        else out.links[c].link_id = "cp04:auxiliary:" + _cdid_selref_num(i - (*plan).M)
        out.links[c].link_id = out.links[c].link_id + "|recipe:" + _cdid_selref_num(m)
        out.links[c].functional_key = _cdid_selref_func_key(*plan, i) + "|reference:" + out.recipe_keys[m]
        out.links[c].pair_id = p_id
        out.links[c].s_cal = (*plan).raw[i, 2]
        out.links[c].t_cal = (*plan).raw[i, 3]
        out.links[c].target_cohort_pos = (*plan).raw[i, 1]
        out.links[c].control_cohort_pos = _cdid_selref_K(*catalog, c)
        out.links[c].recipe_id = m
    }

    out.ctl = cdid_selected_default_controls()
    out.ctl.est_method = (controls.method == "reg" ? "or" : controls.method)
    out.ctl.outcome_pooling = (controls.outcome_pooling == "" ? "separate" : controls.outcome_pooling)
    out.ctl.selection = "conditional_mtar"
    out.ctl.selection_model = controls.selection_model
    out.ctl.selection_pooling = controls.selection_pooling
    out.ctl.or_intercept = (*reference).roles[1].intercept
    out.ctl.eps = controls.eps
    out.ctl.maxit = controls.maxit
    out.ctl.rank_tol = controls.rank_tol
    out.ctl.max_chunk_rows = controls.chunk_rows
    out.ctl.max_work = controls.max_work
    out.ctl.max_workspace_doubles = controls.max_entries
    out.ctl.points_only = controls.points_only
    /* w_pop feeds only the static-only tau/eta detail channels (R: NULL under dated;
       chained_selection_consumers.R:875-889). It is the population owner weight of the
       reference_recipes reduction (cdid_population.mata:443). This provider never
       requests those details; V stays separate for aggregation. */
    out.w_pop = J(N, 1, 1)
    return(out)
}

void _cdid_selref_release(struct cdid_selref_inputs scalar b)
{
    b.registry.sel_data = cdid_selected_matrix_slot(0, 1)
    b.registry.or_data = cdid_selected_matrix_slot(0, 1)
    b.registry.trt_data = cdid_selected_matrix_slot(0, 1)
    b.registry.weight_data = cdid_selected_matrix_slot(0, 1)
    b.registry.rho_data = cdid_selected_matrix_slot(0, 1)
}

struct cdid_selref_provider scalar cdid_selref_prepare(
    pointer(struct cdid_population scalar) scalar population,
    pointer(struct cdid_reference_inputs scalar) scalar reference,
    pointer(struct cdid_baseline_plan_state scalar) scalar plan,
    pointer(struct cdid_reference_catalog scalar) scalar catalog,
    struct cdid_selref_controls scalar controls)
{
    struct cdid_selref_provider scalar prov
    struct cdid_selref_inputs scalar b
    real scalar i

    b = cdid_selref_build(population, reference, plan, catalog, controls)
    if (rows(b.recipe_keep)) prov.model = cdid_selected_prepare_registry(b.recipes, b.registry, b.pairs, b.cohort_idx,
        (*plan).cohorts, (*plan).calendar, (*population).A, b.w_pop, b.links, b.ctl, population, 1, b.recipe_keep)
    else prov.model = cdid_selected_prepare_registry(b.recipes, b.registry, b.pairs, b.cohort_idx,
        (*plan).cohorts, (*plan).calendar, (*population).A, b.w_pop, b.links, b.ctl, population, 1)
    _cdid_selref_release(b)

    prov.population = population
    prov.reference = reference
    prov.plan = plan
    prov.catalog = catalog
    prov.method = controls.method
    prov.status = "ok"
    prov.profile_kind = "reference_recipes"
    prov.weight_row_domain = "F on core scale psi=(n_A/N)*phi_F with 1/n_A normalization"
    prov.N = (*population).N
    prov.n_F = prov.N
    prov.R_hist = (*population).R
    prov.n_A = (*population).n_A
    prov.normalization_count = (*population).n_A
    prov.weight_rows = prov.N
    prov.chunk_rows = (controls.chunk_rows > 0 ? controls.chunk_rows : prov.N)
    prov.M = (*plan).M
    prov.Q = (*plan).Q
    prov.n_consumers = (*catalog).n_consumers
    prov.n_recipes = rows(b.recipe_reference)
    prov.n_reference_recipes = (*reference).n_recipes
    prov.recipe_reference = b.recipe_reference
    prov.reference_recipe = b.reference_recipe
    prov.consumer_recipe = b.consumer_recipe
    prov.sel_cal = b.sel_cal
    prov.or_cal = b.or_cal
    prov.trt_cal = b.trt_cal
    prov.w_cal = b.w_cal
    prov.scope = b.scope
    prov.recipe_dates = b.recipe_dates
    prov.recipe_keys = b.recipe_keys

    prov.points = prov.model.estimates
    prov.post_consumer = (*catalog).post_consumer
    prov.aux_consumer = (*catalog).aux_consumer
    prov.aux_supported = (*catalog).aux_supported
    prov.aux_reason = (*catalog).aux_reason
    if (prov.M > 0) {
        if (missing(prov.points[prov.post_consumer]))
            _cdid_selref_fail("cdid_selected_public_error", "Selected provider produced nonfinite post catalog point estimate.", 459)
        prov.post_ids_sorted = uniqrows(sort(prov.post_consumer, 1))
        prov.post_map_sorted = _cdid_sel_link_match_sorted(prov.post_consumer, prov.post_ids_sorted)
        prov.post_cache_ready = 1
    }
    else {
        prov.post_ids_sorted = J(0, 1, 0)
        prov.post_map_sorted = J(0, 1, 0)
        prov.post_cache_ready = 1
    }
    prov.eval_ids_sorted = J(0, 1, 0)
    prov.eval_post_map = J(0, 1, 0)
    prov.cached_eval_phi = J(0, 0, .)
    prov.cache_eval_active = 0
    for (i = 1; i <= prov.Q; i++) {
        if (prov.aux_supported[i]) {
            if (prov.points[prov.aux_consumer[i]] >= .) {
                prov.aux_supported[i] = 0
                prov.aux_reason[i] = "nonfinite_point"
            }
        }
    }
    prov.workspace_entries = _cdid_base_cells(1, prov.model.workspace_doubles +
        16 * prov.n_consumers + 8 * prov.n_recipes * (cols(prov.scope) + 4),
        controls.max_entries, (*plan).workspace_entries)
    return(prov)
}

/* O(1) borrowed-owner and lifetime checks; no roster, history or hash scan. */
void _cdid_selref_check(struct cdid_selref_provider scalar f)
{
    if (f.population == NULL | f.reference == NULL | f.plan == NULL | f.catalog == NULL)
        _cdid_selref_fail("owner", "A borrowed population/reference/plan/catalogue object no longer exists.", 498)
    if (!f.model.prepared | f.model.population_ref != f.population)
        _cdid_selref_fail("owner", "The selected model no longer borrows the original population object.", 498)
    if ((*f.population).N != f.N | (*f.population).R != f.R_hist | (*f.reference).population != f.population |
        (*f.reference).N != f.N | (*f.reference).n_recipes != f.n_reference_recipes |
        (*f.plan).M != f.M | (*f.plan).Q != f.Q | (*f.catalog).n_consumers != f.n_consumers)
        _cdid_selref_fail("owner", "Borrowed objects changed after selected-provider preparation.", 498)
}

/* GMM callback, same contract as cdid_selected_post_core: rows lo..hi on core scale
   (psi=(n_A/N)*phi_F), post columns in original plan order. The provider is used
   through its pointer; no per-chunk copy, refit or validation scan. */
real matrix cdid_selref_post_core(pointer() scalar context, real scalar lo, real scalar hi)
{
    pointer(struct cdid_selref_provider scalar) scalar f
    struct cdid_selected_chunk_out scalar ch
    real colvector ids, map

    if (context == NULL)
        _cdid_selref_fail("input", "Missing prepared selected reference provider context.", 498)
    f = context
    _cdid_selref_check(*f)
    if ((*f).M == 0) return(J(hi - lo + 1, 0, 0))
    if ((*f).cache_eval_active & lo == 1 & hi == (*f).n_F & rows((*f).eval_post_map) == (*f).M) {
        ch = cdid_selected_chunk((*f).model, (lo::hi), (*f).eval_ids_sorted, 0, 0)
        (*f).cached_eval_phi = ch.phi
        return((*f).model.lambda :* ch.phi[, (*f).eval_post_map'])
    }
    if ((*f).post_cache_ready & rows((*f).post_map_sorted) == rows((*f).post_consumer)) {
        if (all((*f).post_ids_sorted[(*f).post_map_sorted] :== (*f).post_consumer)) {
            ch = cdid_selected_chunk((*f).model, (lo::hi), (*f).post_ids_sorted, 1, 0)
            return(ch.phi[, (*f).post_map_sorted'])
        }
    }
    ids = uniqrows(sort((*f).post_consumer, 1))
    map = _cdid_sel_link_match_sorted((*f).post_consumer, ids)
    ch = cdid_selected_chunk((*f).model, (lo::hi), ids, 1, 0)
    return(ch.phi[, map'])
}

/* Original-F influence rows for requested catalogue consumers (frame scale when
   scale_core=0). Bounded in the requested chunk and consumer batch. */
real matrix cdid_selref_if_chunk(struct cdid_selref_provider scalar f, real colvector unit_ids,
    real colvector consumers, real scalar scale_core)
{
    struct cdid_selected_chunk_out scalar ch
    _cdid_selref_check(f)
    if (rows(consumers) < 1 | cols(consumers) != 1)
        _cdid_selref_fail("input", "Request a column of catalogue consumer IDs.", 498)
    ch = cdid_selected_chunk(f.model, unit_ids, consumers, scale_core, 0)
    return(ch.phi)
}

real colvector cdid_selref_post_points(struct cdid_selref_provider scalar f)
{
    _cdid_selref_check(f)
    return(f.M ? f.points[f.post_consumer] : J(0, 1, .))
}

real colvector cdid_selref_aux_points(struct cdid_selref_provider scalar f)
{
    real colvector out, act
    _cdid_selref_check(f)
    out = J(f.Q, 1, .)
    act = _cdid_selref_ix(f.aux_supported :== 1)
    if (rows(act)) out[act] = f.points[f.aux_consumer[act]]
    return(out)
}

/* First-calendar weights for aggregation shares; distinct from context weights and
   truthfully missing where the first period is absent. */
real colvector cdid_selref_first_calendar_V(struct cdid_selref_provider scalar f)
{
    _cdid_selref_check(f)
    return((*f.reference).V)
}

/* Diagnostic scope-checked reads: rows must be U rows of the recipe's own cohorts. */
void _cdid_selref_rows_in_scope(struct cdid_selref_provider scalar f, real scalar r,
    real colvector units)
{
    real colvector k, U
    real scalar code
    if (r < 1 | r > f.n_recipes | r != floor(r))
        _cdid_selref_fail("input", "Unknown model design block.", 498)
    if (cols(units) != 1 | rows(units) < 1)
        _cdid_selref_fail("input", "Rows must be a nonempty column of original-F positions.", 498)
    if (missing(units) | any(units :< 1) | any(units :> f.N) | any(units :!= floor(units)))
        _cdid_selref_fail("input", "Rows must be original-F positions.", 498)
    k = f.model.cohort_idx[units]
    if (any(k :< 1)) _cdid_selref_fail("support", "Requested row has no known cohort.", 459)
    if (any(f.scope[r, k'] :!= 1))
        _cdid_selref_fail("support", "Requested row lies outside the design block's used-cohort scope.", 459)
    code = f.model.contexts[r].universe_code
    if (f.model.reader_universe == code) U = f.model.reader.U[units]
    else U = f.model.universe_buffers[code].value[units]
    if (any(U :!= 1))
        _cdid_selref_fail("support", "Requested row lies outside the declared selection universe U.", 459)
}

real colvector cdid_selref_weight_rows(struct cdid_selref_provider scalar f, real scalar r,
    real colvector units)
{
    _cdid_selref_check(f)
    _cdid_selref_rows_in_scope(f, r, units)
    return(_cdid_sld_w_rows(f.model, r, units))
}

real matrix cdid_selref_selection_rows(struct cdid_selref_provider scalar f, real scalar r,
    real colvector units)
{
    real scalar id
    _cdid_selref_check(f)
    _cdid_selref_rows_in_scope(f, r, units)
    id = f.model.contexts[r].sel_id
    if (id == f.model.reader_sel_slot) return(f.model.reader.Z[units, ])
    return(f.model.profiles.sel_data[id].value[units, ])
}

real matrix cdid_selref_outcome_rows(struct cdid_selref_provider scalar f, real scalar r,
    real colvector units)
{
    real matrix value
    _cdid_selref_check(f)
    if (f.method == "ipw") _cdid_selref_fail("input", "Pure IPW reads no outcome design.", 498)
    _cdid_selref_rows_in_scope(f, r, units)
    value = _cdid_sld_or_rows(f.model, r, units)
    if (missing(value))
        _cdid_selref_fail("support", "Outcome design requested outside the design block's actual observed-pair support.", 459)
    return(value)
}

real matrix cdid_selref_treatment_rows(struct cdid_selref_provider scalar f, real scalar r,
    real colvector units)
{
    real matrix value
    _cdid_selref_check(f)
    if (f.method == "reg") _cdid_selref_fail("input", "Outcome regression reads no treatment design.", 498)
    _cdid_selref_rows_in_scope(f, r, units)
    value = _cdid_sld_trt_rows(f.model, r, units)
    if (missing(value))
        _cdid_selref_fail("support", "Treatment design requested outside the design block's actual observed-pair support.", 459)
    return(value)
}

/* Model records for qualification: counts, selection cells, recipe cohort law. */
string rowvector cdid_selref_count_names()
{
    return(("n_recipes", "n_reference_recipes", "n_consumers", "n_selection_fits",
        "n_selection_blocks", "n_treatment_fits", "n_pooled_OR_fits", "n_supports",
        "n_OR_blocks", "n_treatment_models", "n_selection_slots", "n_outcome_slots",
        "n_treatment_slots", "n_weight_slots"))
}

real rowvector cdid_selref_counts(struct cdid_selref_provider scalar f)
{
    return((f.n_recipes, f.n_reference_recipes, f.n_consumers, f.model.n_selection_fit_calls,
        f.model.n_selection_blocks, f.model.n_treatment_fit_calls, f.model.n_or_block_fits,
        f.model.n_supports, f.model.n_blocks, f.model.n_trt_models,
        rows(f.sel_cal), rows(f.or_cal), rows(f.trt_cal), rows(f.w_cal)))
}

/* Columns: recipe, reference recipe, pair cohort, cohort position, s_cal, t_cal,
   cell_s, cell_t, block_s, block_t. */
real matrix cdid_selref_cells(struct cdid_selref_provider scalar f)
{
    real matrix out, pc
    real scalar r, n_tot, at, nr
    n_tot = 0
    for (r = 1; r <= f.n_recipes; r++) n_tot = n_tot + rows(f.model.contexts[r].sel_fit.map.pair_cohorts)
    out = J(n_tot, 10, .)
    at = 0
    for (r = 1; r <= f.n_recipes; r++) {
        pc = f.model.contexts[r].sel_fit.map.pair_cohorts
        nr = rows(pc)
        if (nr > 0) {
            out[|(at + 1), 1 \ (at + nr), 10|] = (J(nr, 1, r), J(nr, 1, f.recipe_reference[r]),
                (1::nr), pc[, 1..3], pc[, 6..7],
                f.model.contexts[r].sel_fit.cell_block[pc[, 6]],
                f.model.contexts[r].sel_fit.cell_block[pc[, 7]])
            at = at + nr
        }
    }
    return(out)
}

real colvector cdid_selref_theta(struct cdid_selref_provider scalar f, real scalar r)
{
    if (r < 1 | r > f.n_recipes | r != floor(r)) _cdid_selref_fail("input", "Unknown model design block.", 498)
    return(f.model.contexts[r].sel_fit.theta)
}

/* Columns: recipe, reference recipe, cohort position, mass_U, repeat_mass, a_U on the
   recipe's own used cohorts (R stats$mass_U/repeat_mass/a). */
real matrix cdid_selref_rates(struct cdid_selref_provider scalar f)
{
    real matrix out
    real colvector hit
    real scalar r, n_tot, at, nr
    n_tot = sum(f.scope :!= 0)
    out = J(n_tot, 6, .)
    at = 0
    for (r = 1; r <= f.n_recipes; r++) {
        hit = _cdid_selref_ix(f.scope[r, ] :!= 0)
        nr = rows(hit)
        if (nr > 0) {
            out[|(at + 1), 1 \ (at + nr), 6|] = (J(nr, 1, r), J(nr, 1, f.recipe_reference[r]),
                hit, f.model.contexts[r].mass_U[hit], f.model.contexts[r].repeat_mass[hit],
                f.model.contexts[r].a_U[hit])
            at = at + nr
        }
    }
    return(out)
}

end
