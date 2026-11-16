*! dated reader and recipe-scoped Stage 1 actual-change imputation; StataNow 19.5.
*! Pinned R cdid 1.1.0.0 ():
*!   R/chained_reference_inputs.R (chained_reference_recipes, build_chained_reference_reader,
*!   chained_reference_x, chained_reference_weight, chained_reference_stat) and
*!   R/chained_imputation.R (imputation_build_catalog, imputation_fit_core Stage 1,
*!   imputation_evaluate_if_chunk_core and its nested adjoint closures).
*! Load cdid_population, cdid_nuisance, cdid_ipw_link, and cdid_imputation_model first.
*! Preserves Stage 1 cdid_imputation_model.mata unchanged; provides borrowed
*! unique-date/weight-anchor slot registry, used-cohort-scoped reader, and dated Stage 1
*! fit/predict/score/adjoint replay without N_layout duplication or full-F unused-cohort reads.
version 19.5
mata:
mata set matastrict on

/* One unique raw design or weight-anchor matrix slot on original F (N_F rows).
   calendar_index = 0 denotes a static (undated) slot; 1..T indexes pop.calendar. */
struct cdid_imp_matrix_slot {
    real scalar calendar_index, calendar_date
    real matrix value
}

/* Registry of unique canonical raw outcome-slope designs (z_data), treatment designs
   (x_data), and weight-anchor vectors (w_data) on original F. Formulas are encoded once
   on the canonical required union and fixed across all date/anchor slots. */
struct cdid_imp_dated_registry {
    real scalar N, p_raw, q, intercept, frozen_encoding
    string rowvector z_cols, x_cols
    string scalar formula_key
    struct cdid_imp_matrix_slot colvector z_data
    struct cdid_imp_matrix_slot colvector x_data
    struct cdid_imp_matrix_slot colvector w_data
}

/* One deduplicated dated recipe (design_date, weight_anchor_date) referencing shared slots
   in the borrowed cdid_imp_dated_registry without duplicating N_F matrices per recipe. */
struct cdid_imp_recipe_spec {
    real scalar recipe_id, design_cal, weight_cal, design_date, weight_date
    real scalar z_slot, x_slot, w_slot
    string scalar recipe_key
}

/* Borrowed used-cohort-scoped dated reader on original F.
   pop and reg are BORROWED pointers; no N_layout owner is ever allocated. */
struct cdid_imp_dated_reader {
    pointer(struct cdid_population scalar) scalar pop
    pointer(struct cdid_imp_dated_registry scalar) scalar reg
    real scalar N, n_A, lambda, G, n_recipes, p_raw, q, intercept, has_ps, is_dated
    real colvector cohorts, unit_cohort_index
    string colvector cohort_labels
    struct cdid_imp_recipe_spec colvector recipes
    real matrix used_cohorts, layout_memberships
    real scalar n_memberships, n_layout_rows
    real matrix cohort_mass_F, cohort_repeat_mass_F, cohort_a_F, cohort_colmax_F
    real matrix elig_keep
    real scalar borrow_sum, owned_entries, borrowed_entries, workspace_bytes
}

void _cdid_impd_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

struct cdid_imp_matrix_slot scalar cdid_imp_dated_slot_make(
    real scalar cal_idx,
    real scalar cal_date,
    real matrix val)
{
    struct cdid_imp_matrix_slot scalar s
    s.calendar_index = cal_idx
    s.calendar_date = cal_date
    s.value = val
    return(s)
}

struct cdid_imp_dated_registry scalar cdid_imp_dated_reg_init(
    real scalar N,
    real scalar p_raw,
    real scalar q,
    real scalar intercept,
    string rowvector z_cols,
    string rowvector x_cols,
    | string scalar formula_key)
{
    struct cdid_imp_dated_registry scalar reg
    real scalar j
    if (N < 1 | N >= . | N != floor(N) | p_raw < 0 | p_raw >= . | p_raw != floor(p_raw) |
        q < 0 | q >= . | q != floor(q) | (intercept != 0 & intercept != 1))
        _cdid_impd_fail("Invalid dimensions or intercept flag in cdid_imp_dated_reg_init.", 498)
    if (!intercept & p_raw == 0)
        _cdid_impd_fail("No-intercept outcome design is empty.", 459)
    reg.N = N
    reg.p_raw = p_raw
    reg.q = q
    reg.intercept = intercept
    reg.frozen_encoding = 1
    if (cols(z_cols) == p_raw) reg.z_cols = z_cols
    else {
        reg.z_cols = J(1, p_raw, "")
        for (j=1; j<=p_raw; j++) reg.z_cols[j] = "x" + strofreal(j, "%21.0g")
    }
    if (cols(x_cols) == q) reg.x_cols = x_cols
    else {
        reg.x_cols = J(1, q, "")
        for (j=1; j<=q; j++) reg.x_cols[j] = (intercept & j == 1 ? "(Intercept)" : "x" + strofreal(j - intercept, "%21.0g"))
    }
    if (cols(reg.z_cols) > 1) {
        if (cols(uniqrows(reg.z_cols')') != cols(reg.z_cols))
            _cdid_impd_fail("Duplicate column names in canonical outcome slope design.", 498)
    }
    if (cols(reg.x_cols) > 1) {
        if (cols(uniqrows(reg.x_cols')') != cols(reg.x_cols))
            _cdid_impd_fail("Duplicate column names in canonical treatment design.", 498)
    }
    reg.formula_key = (args() < 7 | formula_key == "" ? "canonical_frozen" : formula_key)
    reg.z_data = cdid_imp_matrix_slot(0, 1)
    reg.x_data = cdid_imp_matrix_slot(0, 1)
    reg.w_data = cdid_imp_matrix_slot(0, 1)
    return(reg)
}

real scalar _cdid_impd_find_slot(
    struct cdid_imp_matrix_slot colvector slots,
    real scalar cal_idx)
{
    real scalar i
    for (i=1; i<=rows(slots); i++) {
        if (slots[i].calendar_index == cal_idx) return(i)
    }
    return(0)
}

void cdid_imp_dated_reg_add_z(
    struct cdid_imp_dated_registry scalar reg,
    real scalar cal_idx,
    real scalar cal_date,
    real matrix Z)
{
    real scalar pos, n, j
    struct cdid_imp_matrix_slot colvector next
    if (cal_idx < 0 | cal_idx >= . | cal_idx != floor(cal_idx))
        _cdid_impd_fail("Invalid calendar_index for outcome design slot.", 498)
    if (rows(Z) != reg.N | cols(Z) != reg.p_raw)
        _cdid_impd_fail("Outcome design slot must have N_F rows and fixed p_raw columns.", 498)
    pos = _cdid_impd_find_slot(reg.z_data, cal_idx)
    if (pos > 0)
        _cdid_impd_fail("Duplicate outcome design slot for calendar_index.", 498)
    n = rows(reg.z_data)
    for (j=1; j<=n; j++) {
        if (reg.z_data[j].calendar_date == cal_date)
            _cdid_impd_fail("Duplicate outcome design slot for calendar_date.", 498)
    }
    next = cdid_imp_matrix_slot(n + 1, 1)
    if (n > 0) next[1..n] = reg.z_data
    next[n + 1] = cdid_imp_dated_slot_make(cal_idx, cal_date, Z)
    reg.z_data = next
}

void cdid_imp_dated_reg_add_x(
    struct cdid_imp_dated_registry scalar reg,
    real scalar cal_idx,
    real scalar cal_date,
    real matrix X)
{
    real scalar pos, n, j
    struct cdid_imp_matrix_slot colvector next
    if (cal_idx < 0 | cal_idx >= . | cal_idx != floor(cal_idx))
        _cdid_impd_fail("Invalid calendar_index for treatment design slot.", 498)
    if (rows(X) != reg.N | cols(X) != reg.q)
        _cdid_impd_fail("Treatment design slot must have N_F rows and fixed q columns.", 498)
    pos = _cdid_impd_find_slot(reg.x_data, cal_idx)
    if (pos > 0)
        _cdid_impd_fail("Duplicate treatment design slot for calendar_index.", 498)
    n = rows(reg.x_data)
    for (j=1; j<=n; j++) {
        if (reg.x_data[j].calendar_date == cal_date)
            _cdid_impd_fail("Duplicate treatment design slot for calendar_date.", 498)
    }
    next = cdid_imp_matrix_slot(n + 1, 1)
    if (n > 0) next[1..n] = reg.x_data
    next[n + 1] = cdid_imp_dated_slot_make(cal_idx, cal_date, X)
    reg.x_data = next
}

void cdid_imp_dated_reg_add_w(
    struct cdid_imp_dated_registry scalar reg,
    real scalar cal_idx,
    real scalar cal_date,
    real colvector W)
{
    real scalar pos, n, j
    struct cdid_imp_matrix_slot colvector next
    if (cal_idx < 0 | cal_idx >= . | cal_idx != floor(cal_idx))
        _cdid_impd_fail("Invalid calendar_index for weight-anchor slot.", 498)
    if (rows(W) != reg.N | cols(W) != 1)
        _cdid_impd_fail("Weight-anchor slot must be an N_F x 1 column vector.", 498)
    pos = _cdid_impd_find_slot(reg.w_data, cal_idx)
    if (pos > 0)
        _cdid_impd_fail("Duplicate weight-anchor slot for calendar_index.", 498)
    n = rows(reg.w_data)
    for (j=1; j<=n; j++) {
        if (reg.w_data[j].calendar_date == cal_date)
            _cdid_impd_fail("Duplicate weight-anchor slot for calendar_date.", 498)
    }
    next = cdid_imp_matrix_slot(n + 1, 1)
    if (n > 0) next[1..n] = reg.w_data
    next[n + 1] = cdid_imp_dated_slot_make(cal_idx, cal_date, W)
    reg.w_data = next
}

struct cdid_imp_recipe_spec scalar cdid_imp_recipe_spec_make(
    real scalar recipe_id,
    real scalar design_cal,
    real scalar weight_cal,
    real scalar design_date,
    real scalar weight_date,
    real scalar z_slot,
    real scalar x_slot,
    real scalar w_slot,
    | string scalar recipe_key)
{
    struct cdid_imp_recipe_spec scalar r
    r.recipe_id = recipe_id
    r.design_cal = design_cal
    r.weight_cal = weight_cal
    r.design_date = design_date
    r.weight_date = weight_date
    r.z_slot = z_slot
    r.x_slot = x_slot
    r.w_slot = w_slot
    if (args() >= 9 & recipe_key != "") r.recipe_key = recipe_key
    else {
        r.recipe_key = "design:" + strofreal(design_date, "%21.0g") +
            "|weights:" + strofreal(weight_date, "%21.0g")
    }
    return(r)
}

/* Derive canonical cohorts (never-treated 0 first, then ascending finite cohorts)
   and unit_cohort_index on original F from struct cdid_population. */
void _cdid_impd_cohorts(
    struct cdid_population scalar pop,
    real colvector cohorts_out,
    string colvector labels_out,
    real colvector unit_ci_out)
{
    real scalar N, G, k
    real colvector u_coh, sort_keys, coh_ord, ix
    N = pop.N
    if (N < 1 | rows(pop.cohort_known) != N | any(pop.cohort_known :!= 1))
        _cdid_impd_fail("Pure imputation requires full representative roster F with known cohort on every unit.", 459)
    u_coh = editmissing(pop.cohort, 0)
    if (any(u_coh :< 0)) _cdid_impd_fail("Cohort values must be nonnegative.", 459)
    cohorts_out = uniqrows(sort(u_coh, 1))
    G = rows(cohorts_out)
    sort_keys = J(G, 1, 0)
    for (k=1; k<=G; k++) sort_keys[k] = (cohorts_out[k] == 0 ? -1 : cohorts_out[k])
    coh_ord = order((sort_keys, (1::G)), (1, 2))
    cohorts_out = cohorts_out[coh_ord]
    labels_out = J(G, 1, "")
    unit_ci_out = J(N, 1, 0)
    for (k=1; k<=G; k++) {
        labels_out[k] = _cdid_imp_cohort_label(cohorts_out[k])
        ix = _cdid_imp_index(u_coh :== cohorts_out[k])
        if (rows(ix) > 0) unit_ci_out[ix] = J(rows(ix), 1, k)
    }
    if (any(unit_ci_out :== 0))
        _cdid_impd_fail("Failed to map unit cohort to canonical cohort dictionary.", 498)
}

/* Derive n_recipes x G binary used_cohorts matrix from link_spec (L x 3 or L x 4:
   s_date, t_date, target_g [, recipe_id]) and link_controls (M x 2: link_idx, ctrl_k). */
real matrix cdid_imp_used_from_links(
    struct cdid_population scalar pop,
    real scalar n_recipes,
    real matrix link_spec,
    real matrix link_controls)
{
    real colvector cohorts, unit_ci
    string colvector labels
    real matrix used
    real scalar L, ell, rec, ci, m, k
    if (n_recipes < 1 | n_recipes >= . | n_recipes != floor(n_recipes))
        _cdid_impd_fail("n_recipes must be a positive integer.", 498)
    _cdid_impd_cohorts(pop, cohorts=., labels=J(0,1,""), unit_ci=.)
    L = rows(link_spec)
    if (L < 1 | (cols(link_spec) != 3 & cols(link_spec) != 4) | missing(link_spec) |
        cols(link_controls) != 2 | missing(link_controls))
        _cdid_impd_fail("Invalid link_spec or link_controls in cdid_imp_used_from_links.", 498)
    used = J(n_recipes, rows(cohorts), 0)
    for (ell=1; ell<=L; ell++) {
        rec = (cols(link_spec) == 4 ? link_spec[ell, 4] : 1)
        if (rec < 1 | rec > n_recipes | rec != floor(rec))
            _cdid_impd_fail("link_spec design block index out of range 1..n_recipes.", 498)
        ci = _cdid_imp_find_exact(cohorts, link_spec[ell, 3])
        if (ci < 1 | cohorts[ci] == 0)
            _cdid_impd_fail("link_spec target cohort must be a known finite treated cohort.", 498)
        used[rec, ci] = 1
    }
    for (m=1; m<=rows(link_controls); m++) {
        ell = link_controls[m, 1]
        if (ell < 1 | ell > L | ell != floor(ell))
            _cdid_impd_fail("link_controls references an invalid link index.", 498)
        rec = (cols(link_spec) == 4 ? link_spec[ell, 4] : 1)
        k = _cdid_imp_find_exact(cohorts, link_controls[m, 2])
        if (k < 1)
            _cdid_impd_fail("link_controls references an unknown control cohort.", 498)
        used[rec, k] = 1
    }
    return(used)
}

/* Compute exact scalar borrow fingerprint over pop, reg, and used (recipe, cohort) slices.
   Unused cohorts are never read so NA values outside used support never pollute the sum. */
real scalar _cdid_impd_borrow_fp(
    struct cdid_population scalar pop,
    struct cdid_imp_dated_registry scalar reg,
    struct cdid_imp_recipe_spec colvector recipes,
    real matrix used_cohorts,
    real colvector unit_ci,
    real scalar has_ps,
    | real matrix keep)
{
    struct _cdid_imp_ivec rowvector coh_u
    real scalar fp, nr, G, r, ci, zs, xs, ws
    real colvector ix
    real matrix sum_w, sum_z, sum_x, done_w, done_z, done_x
    nr = rows(recipes)
    G = cols(used_cohorts)
    coh_u = _cdid_imp_ivec(1, G)
    for (ci=1; ci<=G; ci++) coh_u[ci].rows = _cdid_imp_index(unit_ci :== ci)
    sum_w = J(rows(reg.w_data), G, 0)
    done_w = J(rows(reg.w_data), G, 0)
    sum_z = J(rows(reg.z_data), G, 0)
    done_z = J(rows(reg.z_data), G, 0)
    sum_x = J(max((1, rows(reg.x_data))), G, 0)
    done_x = J(max((1, rows(reg.x_data))), G, 0)
    fp = pop.N + 3 * pop.n_A + 5 * rows(pop.calendar) + 7 * quadsum(pop.calendar) +
         11 * quadsum(pop.A) + 13 * reg.N + 17 * reg.p_raw + 19 * reg.q +
         23 * reg.intercept + 29 * reg.frozen_encoding + 31 * rows(reg.z_data) +
         37 * rows(reg.x_data) + 41 * rows(reg.w_data)
    for (r=1; r<=nr; r++) {
        zs = recipes[r].z_slot
        xs = recipes[r].x_slot
        ws = recipes[r].w_slot
        fp = fp + (r + 0.25) * (recipes[r].recipe_id + 2 * recipes[r].design_cal +
             3 * recipes[r].weight_cal + 5 * zs + 7 * xs + 11 * ws)
        for (ci=1; ci<=G; ci++) {
            if (!used_cohorts[r, ci]) continue
            ix = coh_u[ci].rows
            if (rows(ix) == 0) continue
            if (args() >= 7) if (rows(keep)) {
                ix = select(ix, keep[ix, r] :!= 0)
                if (rows(ix) == 0) continue
                fp = fp + (r * 10 + ci) * quadsum(reg.w_data[ws].value[ix])
                if (reg.p_raw > 0) fp = fp + (r * 10 + ci + 0.1) * quadsum(reg.z_data[zs].value[ix, ])
                if (has_ps & reg.q > 0 & xs > 0) fp = fp + (r * 10 + ci + 0.2) * quadsum(reg.x_data[xs].value[ix, ])
                continue
            }
            if (!done_w[ws, ci]) {
                sum_w[ws, ci] = quadsum(reg.w_data[ws].value[ix])
                done_w[ws, ci] = 1
            }
            fp = fp + (r * 10 + ci) * sum_w[ws, ci]
            if (reg.p_raw > 0) {
                if (!done_z[zs, ci]) {
                    sum_z[zs, ci] = quadsum(reg.z_data[zs].value[ix, ])
                    done_z[zs, ci] = 1
                }
                fp = fp + (r * 10 + ci + 0.1) * sum_z[zs, ci]
            }
            if (has_ps & reg.q > 0 & xs > 0) {
                if (!done_x[xs, ci]) {
                    sum_x[xs, ci] = quadsum(reg.x_data[xs].value[ix, ])
                    done_x[xs, ci] = 1
                }
                fp = fp + (r * 10 + ci + 0.2) * sum_x[xs, ci]
            }
        }
    }
    return(fp)
}

/* Construct borrowed used-cohort-scoped dated reader.
   Validates finite designs and positive finite weights ONLY on used (recipe, cohort)
   memberships; missing values on unused cohorts at a recipe's date/anchor are ignored. */
struct cdid_imp_dated_reader scalar cdid_imp_dated_reader_make(
    pointer(struct cdid_population scalar) scalar pop,
    pointer(struct cdid_imp_dated_registry scalar) scalar reg,
    struct cdid_imp_recipe_spec colvector recipes,
    real matrix used_cohorts,
    real scalar has_ps,
    | struct cdid_imputation_budget budget,
      real matrix keep)
{
    struct cdid_imp_dated_reader scalar rdr
    struct cdid_imputation_budget scalar b
    struct _cdid_imp_ivec rowvector coh_u, coh_u_A
    real scalar N, G, nr, p, q, r, r2, ci, zs, xs, ws, n_mem, n_lay, s, j
    real colvector ix, ix_A, w_sub
    real matrix z_sub, x_sub, mem

    if (pop == NULL | reg == NULL)
        _cdid_impd_fail("Dated reader requires non-NULL population and registry pointers.", 498)
    if (args() >= 6) {
        cdid_imputation_validate_budget(budget)
        b = budget
    }
    else b = cdid_imputation_default_budget(max((1, (*pop).N)), 16)

    if (!(*reg).frozen_encoding)
        _cdid_impd_fail("Dated registry must have a fixed canonical formula encoding.", 498)
    N = (*pop).N
    if ((*reg).N != N)
        _cdid_impd_fail("Registry row count N disagrees with population N_F.", 498)
    if ((*pop).n_A < 1 | (*pop).n_A > N)
        _cdid_impd_fail("Population repeat-observer count n_A must lie in 1..N_F.", 459)
    if (has_ps != 0 & has_ps != 1)
        _cdid_impd_fail("has_ps must be 0 or 1.", 498)

    _cdid_impd_cohorts(*pop, rdr.cohorts=., rdr.cohort_labels=J(0,1,""), rdr.unit_cohort_index=.)
    G = rows(rdr.cohorts)
    nr = rows(recipes)
    p = (*reg).p_raw
    q = (has_ps ? (*reg).q : 0)
    if (nr < 1)
        _cdid_impd_fail("Dated reader requires at least one design block.", 498)
    if (rows(used_cohorts) != nr | cols(used_cohorts) != G | missing(used_cohorts) |
        any((used_cohorts :!= 0) :& (used_cohorts :!= 1)))
        _cdid_impd_fail("used_cohorts must be an n_recipes x G binary matrix.", 498)
    /* Restricted reference eligibility contexts: recipe r reads only the units
       with keep[, r] != 0 (all units when keep is empty). */
    rdr.elig_keep = J(0, 0, .)
    if (args() >= 7) if (rows(keep)) {
        if (rows(keep) != N | cols(keep) != nr | missing(keep))
            _cdid_impd_fail("Recipe eligibility masks must be N_F x n_recipes.", 498)
        rdr.elig_keep = (keep :!= 0)
    }
    if (rows((*reg).z_data) < 1 | rows((*reg).w_data) < 1)
        _cdid_impd_fail("Dated registry requires at least one outcome-design slot and one weight slot.", 498)
    if (has_ps & (q < 1 | rows((*reg).x_data) < 1))
        _cdid_impd_fail("ipw and dr require positive treatment column count q and treatment design slots.", 498)

    /* Verify slot shapes and uniqueness in registry */
    for (s=1; s<=rows((*reg).z_data); s++) {
        if (rows((*reg).z_data[s].value) != N | cols((*reg).z_data[s].value) != p)
            _cdid_impd_fail("Outcome design slot dimensions disagree with N_F x p_raw.", 498)
    }
    for (s=1; s<=rows((*reg).x_data); s++) {
        if (rows((*reg).x_data[s].value) != N | cols((*reg).x_data[s].value) != (*reg).q)
            _cdid_impd_fail("Treatment design slot dimensions disagree with N_F x q.", 498)
    }
    for (s=1; s<=rows((*reg).w_data); s++) {
        if (rows((*reg).w_data[s].value) != N | cols((*reg).w_data[s].value) != 1)
            _cdid_impd_fail("Weight slot dimensions disagree with N_F x 1.", 498)
    }

    rdr.cohort_mass_F = J(nr, G, 0)
    rdr.cohort_repeat_mass_F = J(nr, G, 0)
    rdr.cohort_a_F = J(nr, G, .)
    rdr.cohort_colmax_F = J(nr * G, p, 0)
    mem = J(nr * G, 2, 0)
    n_mem = 0
    n_lay = 0
    coh_u = _cdid_imp_ivec(1, G)
    coh_u_A = _cdid_imp_ivec(1, G)
    for (ci=1; ci<=G; ci++) {
        coh_u[ci].rows = _cdid_imp_index(rdr.unit_cohort_index :== ci)
        coh_u_A[ci].rows = _cdid_imp_index((rdr.unit_cohort_index :== ci) :& ((*pop).A :== 1))
    }

    for (r=1; r<=nr; r++) {
        if (recipes[r].recipe_id != r)
            _cdid_impd_fail("Recipes must be indexed 1..n_recipes in canonical order.", 498)
        zs = recipes[r].z_slot
        xs = recipes[r].x_slot
        ws = recipes[r].w_slot
        if (zs < 1 | zs > rows((*reg).z_data) | zs != floor(zs) |
            ws < 1 | ws > rows((*reg).w_data) | ws != floor(ws))
            _cdid_impd_fail("Recipe references an invalid outcome-design or weight slot index.", 498)
        if (has_ps & (xs < 1 | xs > rows((*reg).x_data) | xs != floor(xs)))
            _cdid_impd_fail("Recipe references an invalid treatment-design slot index.", 498)
        if (recipes[r].design_cal > 0) {
            if (recipes[r].design_cal > rows((*pop).calendar) |
                (*pop).calendar[recipes[r].design_cal] != recipes[r].design_date)
                _cdid_impd_fail("Recipe design_cal/design_date disagrees with population calendar.", 498)
        }
        if (recipes[r].weight_cal > 0) {
            if (recipes[r].weight_cal > rows((*pop).calendar) |
                (*pop).calendar[recipes[r].weight_cal] != recipes[r].weight_date)
                _cdid_impd_fail("Recipe weight_cal/weight_date disagrees with population calendar.", 498)
        }
        for (r2=1; r2<r; r2++) {
            if (recipes[r2].design_cal == recipes[r].design_cal &
                recipes[r2].weight_cal == recipes[r].weight_cal &
                recipes[r2].recipe_key == recipes[r].recipe_key)
                _cdid_impd_fail("Duplicate design block specification in dated reader.", 498)
        }
        if (!any(used_cohorts[r, ] :== 1))
            _cdid_impd_fail("Every declared design block must have at least one used cohort.", 498)

        for (ci=1; ci<=G; ci++) {
            if (!used_cohorts[r, ci]) continue
            n_mem = n_mem + 1
            mem[n_mem, ] = (r, ci)
            ix = coh_u[ci].rows
            if (rows(rdr.elig_keep) & rows(ix)) ix = select(ix, rdr.elig_keep[ix, r])
            if (rows(ix) < 1)
                _cdid_impd_fail("Used cohort " + rdr.cohort_labels[ci] + " has zero units in F" +
                    (rows(rdr.elig_keep) ? " inside its reference eligibility context." : "."), 459)
            n_lay = n_lay + rows(ix)
            w_sub = (*reg).w_data[ws].value[ix]
            if (missing(w_sub) | any(w_sub :<= 0))
                _cdid_impd_fail("Missing, nonfinite, or nonpositive weights on used cohort " +
                    rdr.cohort_labels[ci] + " in recipe " + strofreal(r, "%21.0g") + ".", 459)
            if (p > 0) {
                z_sub = (*reg).z_data[zs].value[ix, ]
                if (missing(z_sub))
                    _cdid_impd_fail("Missing or nonfinite outcome covariates on used cohort " +
                        rdr.cohort_labels[ci] + " in recipe " + strofreal(r, "%21.0g") + ".", 459)
                rdr.cohort_colmax_F[(r - 1) * G + ci, ] = colmax(abs(z_sub))
            }
            if (has_ps & q > 0) {
                x_sub = (*reg).x_data[xs].value[ix, ]
                if (missing(x_sub))
                    _cdid_impd_fail("Missing or nonfinite treatment covariates on used cohort " +
                        rdr.cohort_labels[ci] + " in recipe " + strofreal(r, "%21.0g") + ".", 459)
                if ((*reg).intercept) {
                    if (!all(x_sub[, 1] :== 1))
                        _cdid_impd_fail("Treatment intercept column must equal 1 on all used cohort rows.", 459)
                }
            }
            ix_A = coh_u_A[ci].rows
            if (rows(rdr.elig_keep) & rows(ix_A)) ix_A = select(ix_A, rdr.elig_keep[ix_A, r])
            rdr.cohort_mass_F[r, ci] = quadsum(w_sub) / N
            rdr.cohort_repeat_mass_F[r, ci] = (rows(ix_A) ? quadsum((*reg).w_data[ws].value[ix_A]) / N : 0)
            if (!(rdr.cohort_mass_F[r, ci] > 0) | rdr.cohort_mass_F[r, ci] >= . |
                !(rdr.cohort_repeat_mass_F[r, ci] > 0) | rdr.cohort_repeat_mass_F[r, ci] >= .)
                _cdid_impd_fail("Used cohort " + rdr.cohort_labels[ci] + " in recipe " +
                    strofreal(r, "%21.0g") + " has zero full-F or repeat-observation mass (a_k <= 0).", 459)
            rdr.cohort_a_F[r, ci] = rdr.cohort_repeat_mass_F[r, ci] / rdr.cohort_mass_F[r, ci]
        }
    }

    rdr.pop = pop
    rdr.reg = reg
    rdr.N = N
    rdr.n_A = (*pop).n_A
    rdr.lambda = (*pop).n_A / N
    rdr.G = G
    rdr.n_recipes = nr
    rdr.p_raw = p
    rdr.q = q
    rdr.intercept = (*reg).intercept
    rdr.has_ps = has_ps
    rdr.is_dated = 1
    rdr.recipes = recipes
    rdr.used_cohorts = used_cohorts
    rdr.n_memberships = n_mem
    rdr.n_layout_rows = n_lay
    rdr.layout_memberships = (n_mem > 0 ? mem[1..n_mem, ] : J(0, 2, 0))
    rdr.owned_entries = N + G + nr * G * (4 + p) + 2 * n_mem + 9 * nr + 32
    rdr.borrowed_entries = N * (rows((*reg).z_data) * p + (has_ps ? rows((*reg).x_data) * q : 0) +
        rows((*reg).w_data) + 4)
    rdr.workspace_bytes = 8 * (rdr.owned_entries + rdr.borrowed_entries)
    _cdid_imp_check_budget(b, rdr.workspace_bytes, 8 * N * max((1, p + q)),
        16 * n_lay * max((1, p + q + 2)), "dated reader construction")
    rdr.borrow_sum = _cdid_impd_borrow_fp(*pop, *reg, recipes, used_cohorts, rdr.unit_cohort_index, has_ps, rdr.elig_keep)
    return(rdr)
}

void _cdid_impd_check_reader(struct cdid_imp_dated_reader scalar rdr)
{
    real scalar fp
    if (rdr.pop == NULL | rdr.reg == NULL)
        _cdid_impd_fail("Dated reader borrowed population or registry pointer is NULL; lifetime ended.", 498)
    if ((*(rdr.pop)).N != rdr.N | (*(rdr.pop)).n_A != rdr.n_A | (*(rdr.reg)).N != rdr.N |
        (*(rdr.reg)).p_raw != rdr.p_raw | (rdr.has_ps & (*(rdr.reg)).q != rdr.q))
        _cdid_impd_fail("Borrowed population or registry dimensions changed after dated reader creation.", 498)
    fp = _cdid_impd_borrow_fp(*(rdr.pop), *(rdr.reg), rdr.recipes, rdr.used_cohorts,
        rdr.unit_cohort_index, rdr.has_ps, rdr.elig_keep)
    if (fp != rdr.borrow_sum)
        _cdid_impd_fail("Borrowed population or dated registry mutated after reader creation.", 498)
}

void _cdid_impd_check_rows(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real colvector f_rows)
{
    real scalar B, i, ci
    real colvector ci_vec
    if (r < 1 | r > rdr.n_recipes | r != floor(r))
        _cdid_impd_fail("Recipe index out of range 1..n_recipes in dated reader.", 498)
    B = rows(f_rows)
    if (B < 1 | cols(f_rows) != 1 | missing(f_rows))
        _cdid_impd_fail("Requested original-F row vector must be a nonempty finite column vector.", 498)
    if (any(f_rows :< 1) | any(f_rows :> rdr.N) | any(f_rows :!= floor(f_rows)))
        _cdid_impd_fail("Requested row index out of original-F range 1..N_F.", 498)
    ci_vec = rdr.unit_cohort_index[f_rows]
    if (any(rdr.used_cohorts[r, ci_vec'] :== 0)) {
        for (i=1; i<=B; i++) {
            ci = ci_vec[i]
            if (!rdr.used_cohorts[r, ci])
                _cdid_impd_fail("Dated reader queried outside design block " + strofreal(r, "%21.0g") +
                    " used-cohort support (cohort " + rdr.cohort_labels[ci] + ").", 459)
        }
    }
}

real colvector _cdid_impd_w(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real colvector f_rows,
    | real scalar trusted)
{
    real scalar ws
    if (args() < 4 | trusted == 0 | trusted >= .) _cdid_impd_check_rows(rdr, r, f_rows)
    ws = rdr.recipes[r].w_slot
    return((*(rdr.reg)).w_data[ws].value[f_rows])
}

real matrix _cdid_impd_z(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real colvector f_rows,
    | real rowvector cc,
      real scalar trusted)
{
    real scalar zs
    if (args() < 5 | trusted == 0 | trusted >= .) _cdid_impd_check_rows(rdr, r, f_rows)
    zs = rdr.recipes[r].z_slot
    if (args() >= 4) {
        if (cols(cc) == 0) return(J(rows(f_rows), 0, 0))
        return((*(rdr.reg)).z_data[zs].value[f_rows, cc])
    }
    if (rdr.p_raw == 0) return(J(rows(f_rows), 0, 0))
    return((*(rdr.reg)).z_data[zs].value[f_rows, ])
}

real matrix _cdid_impd_z_all(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real colvector f_rows,
    | real scalar trusted)
{
    real scalar zs
    if (args() < 4 | trusted == 0 | trusted >= .) _cdid_impd_check_rows(rdr, r, f_rows)
    if (rdr.p_raw == 0) return(J(rows(f_rows), 0, 0))
    zs = rdr.recipes[r].z_slot
    return((*(rdr.reg)).z_data[zs].value[f_rows, ])
}

real matrix _cdid_impd_x(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real colvector f_rows,
    | real scalar trusted)
{
    real scalar xs
    if (!rdr.has_ps)
        _cdid_impd_fail("Dated reader treatment design queried when has_ps == 0.", 498)
    if (args() < 4 | trusted == 0 | trusted >= .) _cdid_impd_check_rows(rdr, r, f_rows)
    xs = rdr.recipes[r].x_slot
    return((*(rdr.reg)).x_data[xs].value[f_rows, ])
}

real scalar _cdid_impd_ak(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real scalar ci)
{
    if (r < 1 | r > rdr.n_recipes | r != floor(r) | ci < 1 | ci > rdr.G | ci != floor(ci))
        _cdid_impd_fail("Invalid design block or cohort index in _cdid_impd_ak.", 498)
    if (!rdr.used_cohorts[r, ci])
        _cdid_impd_fail("Repeat rate a_k queried for unused (design block, cohort) membership.", 459)
    return(rdr.cohort_a_F[r, ci])
}

real scalar _cdid_impd_mass(
    struct cdid_imp_dated_reader scalar rdr,
    real scalar r,
    real scalar ci)
{
    if (r < 1 | r > rdr.n_recipes | r != floor(r) | ci < 1 | ci > rdr.G | ci != floor(ci))
        _cdid_impd_fail("Invalid design block or cohort index in _cdid_impd_mass.", 498)
    if (!rdr.used_cohorts[r, ci])
        _cdid_impd_fail("Cohort mass queried for unused (design block, cohort) membership.", 459)
    return(rdr.cohort_mass_F[r, ci])
}

/* Verify alias reconstruction on observed-pair training rows AND used-cohort full-F
   prediction rows for a group using the scoped dated reader. */
void _cdid_impd_verify_group_alias(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar gid,
    real rowvector retained,
    real rowvector proposed,
    real matrix map,
    real matrix constants_Z,
    real rowvector colmax_prop,
    real scalar thr,
    real rowvector ctl_err_out,
    real rowvector trt_err_out,
    real scalar rows_checked_out)
{
    real scalar n_prop, B_r, loc, cid, ci, r, pid, n_coh, lo, hi, j, den, ratio
    real colvector cids, coh_rows, r_chunk, obs_ix, u, s0
    real matrix hit, Z_sub, u_mat, s0_mat
    n_prop = cols(proposed)
    ctl_err_out = J(1, n_prop, 0)
    trt_err_out = J(1, n_prop, 0)
    rows_checked_out = 0
    if (n_prop == 0) return
    cids = model.groups[gid].cell_ids
    B_r = rows(cids)
    for (loc=1; loc<=B_r; loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        pid = model.cells[cid].pair_id
        r = model.cells[cid].recipe
        coh_rows = _cdid_imp_rows_rc(model, r, ci)
        n_coh = rows(coh_rows)
        for (lo=1; lo<=n_coh; lo=lo+model.budget.max_chunk_rows) {
            hi = min((n_coh, lo + model.budget.max_chunk_rows - 1))
            r_chunk = coh_rows[lo..hi]
            hit = cdid_imputation_read_pair_chunk(model.reader, model.pop, model.pairs[pid], r_chunk, 0)
            obs_ix = _cdid_imp_index(hit[, 1] :== 1)
            if (rows(obs_ix) == 0) continue
            Z_sub = _cdid_impd_z_all(rdr, r, r_chunk[obs_ix], 1)
            u_mat = Z_sub[, proposed] - (cols(retained) ? Z_sub[, retained] * map : J(rows(Z_sub), n_prop, 0))
            s0_mat = abs(Z_sub[, proposed]) + (cols(retained) ? abs(Z_sub[, retained]) * abs(map) : J(rows(Z_sub), n_prop, 0))
            for (j=1; j<=n_prop; j++) {
                u = abs(u_mat[, j] :- constants_Z[loc, j])
                s0 = s0_mat[, j] :+ abs(constants_Z[loc, j])
                den = max((colmax_prop[j], 2.2250738585072014e-308))
                ratio = max(u :/ rowmax((s0, J(rows(u), 1, den))))
                if (ratio >= .) ratio = .
                if (ratio > ctl_err_out[j] | ratio >= .) ctl_err_out[j] = ratio
            }
            rows_checked_out = rows_checked_out + rows(Z_sub)
        }
    }
    for (loc=1; loc<=B_r; loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        r = model.cells[cid].recipe
        coh_rows = _cdid_imp_rows_rc(model, r, ci)
        n_coh = rows(coh_rows)
        for (lo=1; lo<=n_coh; lo=lo+model.budget.max_chunk_rows) {
            hi = min((n_coh, lo + model.budget.max_chunk_rows - 1))
            r_chunk = coh_rows[lo..hi]
            Z_sub = _cdid_impd_z_all(rdr, r, r_chunk, 1)
            u_mat = Z_sub[, proposed] - (cols(retained) ? Z_sub[, retained] * map : J(rows(Z_sub), n_prop, 0))
            s0_mat = abs(Z_sub[, proposed]) + (cols(retained) ? abs(Z_sub[, retained]) * abs(map) : J(rows(Z_sub), n_prop, 0))
            for (j=1; j<=n_prop; j++) {
                u = abs(u_mat[, j] :- constants_Z[loc, j])
                s0 = s0_mat[, j] :+ abs(constants_Z[loc, j])
                den = max((colmax_prop[j], 2.2250738585072014e-308))
                ratio = max(u :/ rowmax((s0, J(rows(u), 1, den))))
                if (ratio >= .) ratio = .
                if (ratio > trt_err_out[j] | ratio >= .) trt_err_out[j] = ratio
            }
            rows_checked_out = rows_checked_out + rows(Z_sub)
        }
    }
}

/* Recipe-scoped Stage 1 fitter on the borrowed cdid_imp_dated_reader.
   Reuses QR compression, FWL, rank/alias checks, cdid_pooled_wls_fit, and
   cdid_imputation_adjoint_prepare math while reading Z and W only on each recipe's
   used cohorts. Stores compact 0 x p placeholders in out.recipes[r] instead of
   duplicating N_F x p matrices per recipe. */
struct cdid_imputation_model scalar cdid_imp_fit_dated(
    struct cdid_imp_dated_reader scalar rdr,
    real matrix link_spec,
    real matrix link_controls,
    string scalar pooling,
    real scalar rank_tol,
    struct cdid_imputation_budget scalar budget)
{
    struct cdid_imputation_model scalar out
    struct cdid_imputation_pair rowvector raw_pairs
    struct cdid_wls_result scalar wfit
    real scalar N, G, nr, npair, p, r, pi, m, ci, pid, cid, gid, loc, j, k, rk
    real scalar n_coh, lo, hi, nn, scale, nPseudo, at, thr, orig_norm, res_norm
    real scalar exact_j, k_keep, f_wsum
    real colvector pair_ord, raw_to_canon_pair, link_pair_map
    real colvector coh_rows, r_chunk, obs_ix, ww, yy, cids, g_cohorts, rr, ix
    real colvector pseudo_y, pseudo_w, pseudo_block, v, qy, coef, jhat_chunk, yc
    real rowvector constant_cols, retained, proposed, aliased, aliased_pos, candidates
    real rowvector keep, colmax_prop, ctl_err, trt_err, tau, pivot, f_sumZ
    real matrix active_requests, pair_key_mat, cell_mat, hit, rawZ, D, F_c, Q, R1, V
    real matrix coef_rows, res_fac, A_plus, A_minus, pseudo_X, C, mean_D_raw, map
    real matrix constants_Z, constants_D, Z_keep, map_full, rc_Fw, rc_FZ

    _cdid_impd_check_reader(rdr)
    cdid_imputation_validate_budget(budget)
    if (pooling != "separate" & pooling != "within_cohort" & pooling != "across_cohorts")
        _cdid_imp_fail("imputation_pooling must be separate, within_cohort, or across_cohorts.", 498)
    if (rank_tol <= 0 | rank_tol >= 1 | rank_tol >= .)
        _cdid_imp_fail("Require rank_tol in (0, 1).", 498)

    N = rdr.N
    G = rdr.G
    nr = rdr.n_recipes
    p = rdr.p_raw

    out.cohorts = rdr.cohorts
    out.cohort_labels = rdr.cohort_labels
    out.unit_cohort_index = rdr.unit_cohort_index
    out.coh_u = _cdid_imp_ivec(1, G)
    for (ci=1; ci<=G; ci++) out.coh_u[ci].rows = _cdid_imp_index(out.unit_cohort_index :== ci)
    out.elig_keep = rdr.elig_keep
    out.cohort_mass_F = rdr.cohort_mass_F
    out.cohort_repeat_mass_F = rdr.cohort_repeat_mass_F
    out.cohort_a_F = rdr.cohort_a_F
    out.cohort_colmax_F = rdr.cohort_colmax_F
    out.cohort_tau_F = J(nr, G, 0)
    out.cohort_eta_F = J(nr, G, 0)
    out.mean_w_F = J(nr, 1, 0)
    out.mean_wA_F = J(nr, 1, 0)
    out.recipes = cdid_imputation_recipe(1, nr)
    for (r=1; r<=nr; r++) {
        out.recipes[r].recipe_id = rdr.recipes[r].recipe_id
        out.recipes[r].recipe_key = rdr.recipes[r].recipe_key
        out.recipes[r].weight = J(0, 1, 0)
        out.recipes[r].Z = J(0, p, 0)
        out.mean_w_F[r] = sum(out.cohort_mass_F[r, ])
        out.mean_wA_F[r] = sum(out.cohort_repeat_mass_F[r, ])
        if (out.mean_w_F[r] <= 0 | out.mean_wA_F[r] <= 0)
            _cdid_imp_fail("Total used-cohort F mass and repeat-observer A mass must be positive.", 459)
        out.cohort_tau_F[r, ] = out.cohort_mass_F[r, ] :/ out.mean_w_F[r]
        out.cohort_eta_F[r, ] = out.cohort_repeat_mass_F[r, ] :/ out.mean_wA_F[r]
    }

    cdid_imp_links_to_requests((*(rdr.pop)).calendar, link_spec, link_controls,
        raw_pairs=cdid_imputation_pair(1, 0), active_requests=., link_pair_map=.)

    npair = cols(raw_pairs)
    if (npair < 1) _cdid_imp_fail("Require at least one outcome pair.", 498)
    pair_key_mat = J(npair, 5, 0)
    for (pi=1; pi<=npair; pi++) {
        r = raw_pairs[pi].recipe
        if (r < 1 | r > nr | r != floor(r))
            _cdid_imp_fail("Pair design block index out of range 1..n_recipes.", 498)
        raw_pairs[pi].L = raw_pairs[pi].t_index - raw_pairs[pi].s_index
        if (raw_pairs[pi].pair_label == "")
            raw_pairs[pi].pair_label = "e" + strofreal(raw_pairs[pi].s_index, "%21.0g") +
                strofreal(raw_pairs[pi].t_index, "%21.0g")
        pair_key_mat[pi, ] = (raw_pairs[pi].L, raw_pairs[pi].s_index, raw_pairs[pi].t_index, r, pi)
    }
    pair_ord = order(pair_key_mat, (1, 2, 3, 4, 5))
    out.pairs = cdid_imputation_pair(1, npair)
    raw_to_canon_pair = J(npair, 1, 0)
    for (pi=1; pi<=npair; pi++) {
        out.pairs[pi] = raw_pairs[pair_ord[pi]]
        out.pairs[pi].source_pair_id = raw_pairs[pair_ord[pi]].pair_id
        if (out.pairs[pi].pair_key == "")
            out.pairs[pi].pair_key = "pair:" + strofreal(pair_ord[pi], "%21.0g")
        out.pairs[pi].pair_id = pi
        raw_to_canon_pair[pair_ord[pi]] = pi
    }

    cell_mat = J(G, npair, 0)
    for (m=1; m<=rows(active_requests); m++) {
        ci = _cdid_imp_find_exact(out.cohorts, active_requests[m, 1])
        pi = active_requests[m, 2]
        if (ci == 0 | pi < 1 | pi > npair | pi != floor(pi))
            _cdid_imp_fail("Active cell request references unknown cohort or pair slot.", 498)
        pid = raw_to_canon_pair[pi]
        r = out.pairs[pid].recipe
        if (!rdr.used_cohorts[r, ci])
            _cdid_imp_fail("Active cell request falls outside declared used_cohorts support.", 459)
        cell_mat[ci, pid] = 1
    }
    out.n_cells = sum(cell_mat)
    if (out.n_cells < 1) _cdid_imp_fail("No active imputation cells declared.", 498)
    out.cells = cdid_imputation_cell(1, out.n_cells)
    out.cell_lookup = J(G, npair, 0)
    out.active_cohorts = J(0, 1, 0)
    cid = 0
    /* Order cells by (cohort_index ci, recipe r, canonical pair_id pid) so multi-recipe
       catalogs group each (cohort, recipe) block's pairs in ascending (L, s, t) order. */
    for (ci=1; ci<=G; ci++) {
        if (any(cell_mat[ci, ] :== 1)) out.active_cohorts = out.active_cohorts \ ci
        for (r=1; r<=nr; r++) {
            for (pid=1; pid<=npair; pid++) {
                if (!cell_mat[ci, pid] | out.pairs[pid].recipe != r) continue
                cid = cid + 1
                out.cells[cid].cell_id = cid
                out.cells[cid].cohort_index = ci
                out.cells[cid].cohort = out.cohorts[ci]
                out.cells[cid].cohort_label = out.cohort_labels[ci]
                out.cells[cid].pair_id = pid
                out.cells[cid].recipe = r
                out.cells[cid].s_index = out.pairs[pid].s_index
                out.cells[cid].t_index = out.pairs[pid].t_index
                out.cells[cid].L = out.pairs[pid].L
                out.cells[cid].pair_label = out.pairs[pid].pair_label
                out.cells[cid].name = "k" + out.cohort_labels[ci] + "." + out.pairs[pid].pair_label +
                    ".r" + strofreal(out.recipes[r].recipe_id, "%21.0g")
                out.cell_lookup[ci, pid] = cid
            }
        }
    }

    if (pooling == "separate") {
        out.n_groups = out.n_cells
        out.groups = cdid_imputation_group(1, out.n_groups)
        for (cid=1; cid<=out.n_cells; cid++) {
            out.cells[cid].group_id = cid
            out.cells[cid].local_id = 1
            out.groups[cid].group_id = cid
            out.groups[cid].group_index0 = cid - 1
            out.groups[cid].recipe = out.cells[cid].recipe
            out.groups[cid].n_cells = 1
            out.groups[cid].name = "group" + strofreal(cid - 1, "%21.0g")
            out.groups[cid].pooling = "separate"
            out.groups[cid].cell_ids = cid
            out.groups[cid].cohorts = out.cells[cid].cohort_index
        }
    }
    else if (pooling == "within_cohort") {
        gid = 0
        for (m=1; m<=rows(out.active_cohorts); m++) {
            ci = out.active_cohorts[m]
            for (r=1; r<=nr; r++) {
                cids = J(0, 1, 0)
                for (cid=1; cid<=out.n_cells; cid++) {
                    if (out.cells[cid].cohort_index == ci & out.cells[cid].recipe == r)
                        cids = cids \ cid
                }
                if (rows(cids) > 0) gid = gid + 1
            }
        }
        out.n_groups = gid
        out.groups = cdid_imputation_group(1, out.n_groups)
        gid = 0
        for (m=1; m<=rows(out.active_cohorts); m++) {
            ci = out.active_cohorts[m]
            for (r=1; r<=nr; r++) {
                cids = J(0, 1, 0)
                for (cid=1; cid<=out.n_cells; cid++) {
                    if (out.cells[cid].cohort_index == ci & out.cells[cid].recipe == r)
                        cids = cids \ cid
                }
                if (rows(cids) == 0) continue
                gid = gid + 1
                for (loc=1; loc<=rows(cids); loc++) {
                    out.cells[cids[loc]].group_id = gid
                    out.cells[cids[loc]].local_id = loc
                }
                out.groups[gid].group_id = gid
                out.groups[gid].group_index0 = gid - 1
                out.groups[gid].recipe = r
                out.groups[gid].n_cells = rows(cids)
                out.groups[gid].name = "group" + strofreal(gid - 1, "%21.0g")
                out.groups[gid].pooling = "within_cohort"
                out.groups[gid].cell_ids = cids
                out.groups[gid].cohorts = ci
            }
        }
    }
    else {
        gid = 0
        for (r=1; r<=nr; r++) {
            cids = J(0, 1, 0)
            for (cid=1; cid<=out.n_cells; cid++) {
                if (out.cells[cid].recipe == r) cids = cids \ cid
            }
            if (rows(cids) > 0) gid = gid + 1
        }
        out.n_groups = gid
        out.groups = cdid_imputation_group(1, out.n_groups)
        gid = 0
        for (r=1; r<=nr; r++) {
            cids = J(0, 1, 0)
            g_cohorts = J(0, 1, 0)
            for (cid=1; cid<=out.n_cells; cid++) {
                if (out.cells[cid].recipe == r) {
                    cids = cids \ cid
                    if (!any(g_cohorts :== out.cells[cid].cohort_index))
                        g_cohorts = g_cohorts \ out.cells[cid].cohort_index
                }
            }
            if (rows(cids) == 0) continue
            gid = gid + 1
            for (loc=1; loc<=rows(cids); loc++) {
                out.cells[cids[loc]].group_id = gid
                out.cells[cids[loc]].local_id = loc
            }
            out.groups[gid].group_id = gid
            out.groups[gid].group_index0 = gid - 1
            out.groups[gid].recipe = r
            out.groups[gid].n_cells = rows(cids)
            out.groups[gid].name = "group" + strofreal(gid - 1, "%21.0g")
            out.groups[gid].pooling = "across_cohorts"
            out.groups[gid].cell_ids = cids
            out.groups[gid].cohorts = g_cohorts
        }
    }

    /* Preflight charges borrowed unique-date slots via rdr.workspace_bytes (not nr * N). */
    out.retained_bytes_bound = rdr.workspace_bytes + 8 * (out.n_cells * ((p + 2)^2 + 8 * p + 32) +
        out.n_groups * (4 * p^2 + 8 * p + 32) + budget.max_cached_pairs * 2 * max((1, rdr.n_A)) + 4096)
    out.scratch_bytes_bound = 8 * (min((N, budget.max_chunk_rows)) * (2 * p + 16) +
        2 * (p + 2) * max((1, out.n_cells)) * (p + 4) + 4096)
    out.work_ops_bound = 128 * N * (out.n_cells + 1) * (p + 2)^2 * (1 + ceil(ln(N + 1) / ln(2)))
    _cdid_imp_check_budget(budget, out.retained_bytes_bound, out.scratch_bytes_bound,
        out.work_ops_bound, "Stage 1 dated imputation fitting preflight")

    out.pop = *(rdr.pop)
    out.budget = budget
    out.reader = cdid_imputation_pair_reader_init(budget)
    out.N = N
    out.n_A = rdr.n_A
    out.lambda = rdr.lambda
    out.n_cohorts = G
    out.n_recipes = nr
    out.n_pairs = npair
    out.intercept = rdr.intercept
    out.p_raw = p
    out.rank_tol = rank_tol
    out.pooling = pooling
    out.slope_names = (*(rdr.reg)).z_cols

    /* Stage 1a: Stream observed pairs per cell (visited pair-major for cache locality) */
    rc_Fw = J(nr, G, .)
    rc_FZ = J(nr * G, p, 0)
    for (pid=1; pid<=npair; pid++) {
        for (ci=1; ci<=G; ci++) {
            cid = out.cell_lookup[ci, pid]
            if (cid == 0) continue
            r = out.cells[cid].recipe
            coh_rows = _cdid_imp_rows_rc(out, r, ci)
            n_coh = rows(coh_rows)
            F_c = J(0, p + 2, 0)
            out.cells[cid].n_pair = 0
            out.cells[cid].weight_sum = 0
            out.cells[cid].sumY = 0
            out.cells[cid].sumZ = J(1, p, 0)
            out.cells[cid].xmin = J(1, p, .)
            out.cells[cid].xmax = J(1, p, .)
            out.cells[cid].firstZ = J(1, p, .)
            out.cells[cid].colmax_obs = J(1, p, 0)
            out.cells[cid].F_n = n_coh
            if (rc_Fw[r, ci] >= .) {
                f_wsum = 0
                f_sumZ = J(1, p, 0)
                for (lo=1; lo<=n_coh; lo=lo+budget.max_chunk_rows) {
                    hi = min((n_coh, lo + budget.max_chunk_rows - 1))
                    r_chunk = coh_rows[lo..hi]
                    ww = _cdid_impd_w(rdr, r, r_chunk, 1)
                    f_wsum = f_wsum + sum(ww)
                    if (p > 0) f_sumZ = f_sumZ + quadcross(ww, _cdid_impd_z_all(rdr, r, r_chunk, 1))
                }
                rc_Fw[r, ci] = f_wsum
                if (p > 0) rc_FZ[(r - 1) * G + ci, ] = f_sumZ
            }
            out.cells[cid].F_weight_sum = rc_Fw[r, ci]
            out.cells[cid].F_sum_Z = (p > 0 ? rc_FZ[(r - 1) * G + ci, ] : J(1, 0, 0))
            for (lo=1; lo<=n_coh; lo=lo+budget.max_chunk_rows) {
                hi = min((n_coh, lo + budget.max_chunk_rows - 1))
                r_chunk = coh_rows[lo..hi]
                hit = cdid_imputation_read_pair_chunk(out.reader, out.pop, out.pairs[pid], r_chunk, 1)
                obs_ix = _cdid_imp_index(hit[, 1] :== 1)
                if (rows(obs_ix) == 0) continue
                ww = _cdid_impd_w(rdr, r, r_chunk[obs_ix], 1)
                yy = hit[obs_ix, 2]
                rawZ = _cdid_impd_z_all(rdr, r, r_chunk[obs_ix], 1)
                D = out.cells[cid].L :* rawZ
                if (missing(D)) _cdid_imp_fail("Calendar-scaled imputation design overflowed.", 504)
                if (p > 0) {
                    if (out.cells[cid].n_pair == 0) {
                        out.cells[cid].xmin = colmin(rawZ)
                        out.cells[cid].xmax = colmax(rawZ)
                        out.cells[cid].firstZ = rawZ[1, ]
                    }
                    else {
                        out.cells[cid].xmin = colmin(out.cells[cid].xmin \ rawZ)
                        out.cells[cid].xmax = colmax(out.cells[cid].xmax \ rawZ)
                    }
                    out.cells[cid].colmax_obs = colmax(out.cells[cid].colmax_obs \ abs(rawZ))
                    out.cells[cid].sumZ = out.cells[cid].sumZ + quadcross(ww, rawZ)
                }
                out.cells[cid].n_pair = out.cells[cid].n_pair + rows(obs_ix)
                out.cells[cid].weight_sum = out.cells[cid].weight_sum + sum(ww)
                out.cells[cid].sumY = out.cells[cid].sumY + sum(ww :* yy)
                F_c = _cdid_imp_compress(F_c \ ((J(rows(obs_ix), 1, 1), D, yy) :* sqrt(ww)))
            }
            if (out.cells[cid].n_pair < 1 | out.cells[cid].weight_sum <= 0 |
                out.cells[cid].weight_sum >= . | (out.cells[cid].weight_sum / N) <= 0)
                _cdid_imp_fail("Actual-change imputation cell `" + out.cells[cid].name + "` has zero observed-pair support.", 459)
            if (out.cells[cid].F_weight_sum <= 0 | out.cells[cid].F_weight_sum >= .)
                _cdid_imp_fail("Actual-change imputation cell `" + out.cells[cid].name + "` has zero full-F cohort mass.", 459)
            out.cells[cid].mass = out.cells[cid].weight_sum / N
            out.cells[cid].mean_y = out.cells[cid].sumY / out.cells[cid].weight_sum
            out.cells[cid].F_mass = out.cells[cid].F_weight_sum / N
            out.cells[cid].F_mean_Z = (p > 0 ? out.cells[cid].F_sum_Z :/ out.cells[cid].F_weight_sum : J(1, 0, 0))
            out.cells[cid].qr_factor = F_c
        }
    }

    /* Stage 1b: Fit each pooling group r in R_imp with used-cohort full-F alias verification */
    thr = 64 * epsilon(1) * max((1, p + out.intercept))
    for (gid=1; gid<=out.n_groups; gid++) {
        cids = out.groups[gid].cell_ids
        nPseudo = 2 * (p + 1) * rows(cids)
        pseudo_X = J(nPseudo, p, 0)
        pseudo_y = J(nPseudo, 1, 0)
        pseudo_w = J(nPseudo, 1, 0)
        pseudo_block = J(nPseudo, 1, 0)
        at = 0
        out.groups[gid].n_fit_rows = 0
        for (loc=1; loc<=rows(cids); loc++) {
            cid = cids[loc]
            out.groups[gid].n_fit_rows = out.groups[gid].n_fit_rows + out.cells[cid].n_pair
            F_c = out.cells[cid].qr_factor
            hqrd(F_c[, 1], Q=., tau=., R1=.)
            V = hqrdmultq(Q, tau, F_c[, 2..(p + 2)], 1)
            coef_rows = V[1, ] / R1[1, 1]
            res_fac = V[2..rows(V), ]
            if (out.intercept & p > 0) {
                constant_cols = (out.cells[cid].xmin :== out.cells[cid].firstZ) :&
                                (out.cells[cid].xmax :== out.cells[cid].firstZ)
                for (k=1; k<=p; k++) {
                    if (constant_cols[k]) {
                        coef_rows[k] = out.cells[cid].L * out.cells[cid].firstZ[k]
                        res_fac[, k] = J(rows(res_fac), 1, 0)
                    }
                }
            }
            nn = rows(res_fac)
            scale = sqrt(nn / out.cells[cid].weight_sum)
            A_plus = coef_rows :+ scale :* res_fac
            A_minus = coef_rows :- scale :* res_fac
            rr = (at + 1 .. at + 2 * nn)'
            if (p > 0) pseudo_X[rr, ] = A_plus[, 1..p] \ A_minus[, 1..p]
            pseudo_y[rr] = A_plus[, p + 1] \ A_minus[, p + 1]
            pseudo_w[rr] = J(2 * nn, 1, out.cells[cid].weight_sum / (2 * nn))
            pseudo_block[rr] = J(2 * nn, 1, loc)
            at = at + 2 * nn
        }

        C = pseudo_X
        mean_D_raw = J(rows(cids), p, 0)
        for (loc=1; loc<=rows(cids); loc++) {
            ix = _cdid_imp_index(pseudo_block :== loc)
            if (out.intercept & p > 0) {
                mean_D_raw[loc, ] = quadcross(pseudo_w[ix], pseudo_X[ix, ]) / sum(pseudo_w[ix])
                C[ix, ] = pseudo_X[ix, ] :- mean_D_raw[loc, ]
                for (j=1; j<=p; j++) {
                    if (all(pseudo_X[ix, j] :== pseudo_X[ix[1], j])) {
                        C[ix, j] = J(rows(ix), 1, 0)
                        mean_D_raw[loc, j] = pseudo_X[ix[1], j]
                    }
                }
            }
        }
        C = C :* sqrt(pseudo_w)
        retained = J(1, 0, 0)
        proposed = J(1, 0, 0)
        for (j=1; j<=p; j++) {
            v = C[, j]
            orig_norm = sqrt(quadcross(v, v))
            if (cols(retained) > 0) {
                hqrdp(C[, retained], Q=., tau=., R1=., pivot=.)
                qy = hqrdmultq1t(Q, tau, v)
                coef = solveupper(R1, qy, rk=., 0)
                coef = coef[invorder(pivot)]
                v = v - C[, retained] * coef
            }
            res_norm = sqrt(quadcross(v, v))
            if (orig_norm == 0 | res_norm <= rank_tol * orig_norm) proposed = (proposed, j)
            else retained = (retained, j)
        }

        aliased = J(1, 0, 0)
        aliased_pos = J(1, 0, 0)
        candidates = J(1, 0, 0)
        map = J(cols(retained), 0, 0)
        constants_Z = J(rows(cids), 0, 0)
        constants_D = J(rows(cids), 0, 0)
        ctl_err = J(1, 0, 0)
        trt_err = J(1, 0, 0)
        out.groups[gid].rows_checked = 0

        if (cols(proposed) > 0) {
            map = J(cols(retained), cols(proposed), 0)
            if (cols(retained) > 0) {
                hqrdp(C[, retained], Q=., tau=., R1=., pivot=.)
                qy = hqrdmultq1t(Q, tau, C[, proposed])
                map = solveupper(R1, qy, rk=., 0)
                map = map[invorder(pivot), ]
                if (!missing(map)) {
                    qy = hqrdmultq1t(Q, tau, C[, proposed] - C[, retained] * map)
                    V = solveupper(R1, qy, rk=., 0)
                    if (!missing(V)) map = map + V[invorder(pivot), ]
                }
            }
            constants_Z = J(rows(cids), cols(proposed), 0)
            if (out.intercept) {
                for (loc=1; loc<=rows(cids); loc++) {
                    cid = cids[loc]
                    constants_Z[loc, ] = (mean_D_raw[loc, proposed] -
                        (cols(retained) ? mean_D_raw[loc, retained] * map : J(1, cols(proposed), 0))) :/ out.cells[cid].L
                }
            }
            colmax_prop = J(1, cols(proposed), 0)
            for (loc=1; loc<=rows(cids); loc++) {
                cid = cids[loc]
                ci = out.cells[cid].cohort_index
                r = out.cells[cid].recipe
                colmax_prop = colmax(colmax_prop \ out.cells[cid].colmax_obs[proposed] \
                    out.cohort_colmax_F[(r - 1) * G + ci, proposed])
            }
            _cdid_impd_verify_group_alias(out, rdr, gid, retained, proposed, map, constants_Z,
                colmax_prop, thr, ctl_err=., trt_err=., out.groups[gid].rows_checked=.)
            for (j=1; j<=cols(proposed); j++) {
                exact_j = (ctl_err[j] <= thr & ctl_err[j] < .)
                if (exact_j & (trt_err[j] > thr | trt_err[j] >= .))
                    _cdid_imp_fail("Outcome alias does not preserve required full-F prediction support in " + out.groups[gid].name + ".", 459)
                if (exact_j & trt_err[j] <= thr) {
                    aliased = (aliased, proposed[j])
                    aliased_pos = (aliased_pos, j)
                }
                else {
                    candidates = (candidates, proposed[j])
                }
            }
            keep = (cols(retained) + cols(candidates) > 0 ? sort((retained, candidates)', 1)' : J(1, 0, 0))
            if (!out.intercept & cols(keep) == 0)
                _cdid_imp_fail("Outcome design has no retained columns on actual training and full-F rows.", 459)
            map_full = J(cols(keep), cols(aliased), 0)
            if (cols(retained) > 0 & cols(aliased) > 0) {
                for (j=1; j<=cols(aliased); j++) {
                    for (k=1; k<=cols(retained); k++) {
                        map_full[_cdid_imp_index(keep :== retained[k]), j] = map[k, aliased_pos[j]]
                    }
                }
            }
            map = map_full
            constants_Z = (cols(aliased) > 0 ? constants_Z[, aliased_pos] : J(rows(cids), 0, 0))
            constants_D = J(rows(cids), cols(aliased), 0)
            if (out.intercept & cols(aliased) > 0) {
                for (loc=1; loc<=rows(cids); loc++) {
                    constants_D[loc, ] = out.cells[cids[loc]].L :* constants_Z[loc, ]
                }
            }
        }
        else {
            keep = retained
        }

        k_keep = cols(keep)
        if (cols(candidates) == 0) {
            wfit = cdid_pooled_wls_fit(
                (k_keep > 0 ? pseudo_X[, keep] : J(nPseudo, 0, 0)),
                pseudo_y, pseudo_w, pseudo_block, J(nPseudo, 1, 1), N,
                out.intercept, J(0, k_keep, 0), J(0, 1, 0), rank_tol)
            if (cols(wfit.kept) != k_keep)
                _cdid_imp_fail("Pooled WLS final rank refused; no ridge or parameter rescue.", 459)
            out.groups[gid].beta = wfit.beta
            out.groups[gid].alpha = wfit.alpha
            out.groups[gid].R = wfit.R
            out.groups[gid].pivot = wfit.pivot
            out.groups[gid].rcond_R = wfit.rcond_R
            out.groups[gid].block_mass = wfit.mass
            out.groups[gid].mean_y = wfit.mean_y
            out.groups[gid].mean_D = (out.intercept & k_keep > 0 ? wfit.mean_X : J(rows(cids), k_keep, 0))
            out.groups[gid].rss = wfit.rss
        }
        else {
            yc = pseudo_y
            out.groups[gid].block_mass = J(rows(cids), 1, 0)
            out.groups[gid].mean_y = J(rows(cids), 1, 0)
            for (loc=1; loc<=rows(cids); loc++) {
                ix = _cdid_imp_index(pseudo_block :== loc)
                out.groups[gid].block_mass[loc] = sum(pseudo_w[ix])
                out.groups[gid].mean_y[loc] = sum(pseudo_w[ix] :* pseudo_y[ix]) / out.groups[gid].block_mass[loc]
                if (out.intercept) yc[ix] = pseudo_y[ix] :- out.groups[gid].mean_y[loc]
            }
            yc = yc :* sqrt(pseudo_w)
            Z_keep = C[, keep]
            hqrdp(Z_keep, Q=., tau=., R1=., pivot=.)
            if (rows(R1) != k_keep | missing(R1) | min(abs(diagonal(R1))) <= rank_tol * max(abs(diagonal(R1))))
                _cdid_imp_fail("Pooled WLS final rank refused; no ridge or parameter rescue.", 459)
            qy = hqrdmultq1t(Q, tau, yc)
            coef = solveupper(R1, qy, rk=., 0)
            if (rk != k_keep | missing(coef))
                _cdid_imp_fail("Nonfinite pooled WLS solve.", 504)
            out.groups[gid].beta = coef[invorder(pivot)]
            out.groups[gid].R = R1
            out.groups[gid].pivot = pivot
            out.groups[gid].rcond_R = _cdid_nuisance_rcond(R1)
            out.groups[gid].mean_D = (out.intercept ? mean_D_raw[, keep] : J(rows(cids), k_keep, 0))
            out.groups[gid].alpha = (out.intercept ? out.groups[gid].mean_y - out.groups[gid].mean_D * out.groups[gid].beta : J(rows(cids), 1, 0))
            v = pseudo_y - out.groups[gid].alpha[pseudo_block, 1] - pseudo_X[, keep] * out.groups[gid].beta
            out.groups[gid].rss = sum(pseudo_w :* v:^2)
        }

        out.groups[gid].p_raw = p
        out.groups[gid].p = k_keep
        out.groups[gid].intercept = out.intercept
        out.groups[gid].objective = out.groups[gid].rss / N
        out.groups[gid].rank_tol = rank_tol
        out.groups[gid].alias_threshold = thr
        out.groups[gid].kept = keep
        out.groups[gid].retained = retained
        out.groups[gid].proposed = proposed
        out.groups[gid].dropped = aliased
        out.groups[gid].candidates = candidates
        out.groups[gid].control_scaled = ctl_err
        out.groups[gid].treated_scaled = trt_err
        out.groups[gid].alias_map = map
        out.groups[gid].alias_constants_Z = constants_Z
        out.groups[gid].alias_constants_D = constants_D

        for (loc=1; loc<=rows(cids); loc++) {
            cid = cids[loc]
            ci = out.cells[cid].cohort_index
            r = out.cells[cid].recipe
            out.cells[cid].alpha = (out.intercept ? out.groups[gid].alpha[loc] : 0)
            out.cells[cid].mean_d = (out.intercept & k_keep > 0 ? out.groups[gid].mean_D[loc, ] : J(1, k_keep, 0))
            out.cells[cid].d_bar_j = (k_keep > 0 ? out.cells[cid].L :* out.cells[cid].F_mean_Z[keep] : J(1, 0, 0))

            coh_rows = _cdid_imp_rows_rc(out, r, ci)
            n_coh = rows(coh_rows)
            out.cells[cid].F_mean_j = 0
            out.cells[cid].F_sum_jZ = J(1, p, 0)
            for (lo=1; lo<=n_coh; lo=lo+budget.max_chunk_rows) {
                hi = min((n_coh, lo + budget.max_chunk_rows - 1))
                r_chunk = coh_rows[lo..hi]
                ww = _cdid_impd_w(rdr, r, r_chunk, 1)
                jhat_chunk = J(rows(r_chunk), 1, out.cells[cid].alpha)
                if (k_keep > 0) {
                    jhat_chunk = jhat_chunk + (out.cells[cid].L :* _cdid_impd_z(rdr, r, r_chunk, keep, 1)) * out.groups[gid].beta
                }
                out.cells[cid].F_mean_j = out.cells[cid].F_mean_j + sum(ww :* jhat_chunk)
                if (p > 0) {
                    out.cells[cid].F_sum_jZ = out.cells[cid].F_sum_jZ + quadcross(ww :* jhat_chunk, _cdid_impd_z_all(rdr, r, r_chunk, 1))
                }
            }
            out.cells[cid].F_mean_j = out.cells[cid].F_mean_j / out.cells[cid].F_weight_sum
            out.cells[cid].F_mean_jZ = (p > 0 ? out.cells[cid].F_sum_jZ :/ out.cells[cid].F_weight_sum : J(1, 0, 0))
            if (missing(out.cells[cid].alpha) | missing(out.cells[cid].F_mean_j) | missing(out.cells[cid].F_mean_jZ))
                _cdid_imp_fail("Nonfinite cell fitted parameters or full-F moments.", 504)
        }
    }
    return(out)
}

/* Predict j_hat_c(X_i) on a bounded chunk of F rows within recipe r(c)'s used-cohort support. */
real colvector cdid_imp_predict_cell_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar cid,
    real colvector f_rows)
{
    real scalar gid, r, B
    real colvector out
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imp_predict_cell_d.", 498)
    B = rows(f_rows)
    gid = model.cells[cid].group_id
    r = model.cells[cid].recipe
    _cdid_impd_check_rows(rdr, r, f_rows)
    out = J(B, 1, model.cells[cid].alpha)
    if (model.groups[gid].p > 0) {
        out = out + (model.cells[cid].L :* _cdid_impd_z(rdr, r, f_rows, model.groups[gid].kept, 1)) * model.groups[gid].beta
    }
    return(out)
}

/* Observed-pair cell residual score on a bounded chunk of F rows (dated reader).
   Reads Z and W only on rows belonging to cell cid's cohort (where S_e can be nonzero). */
real colvector cdid_imp_cell_score_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar cid,
    real colvector f_rows)
{
    real scalar B, ci, pid, r, gid
    real colvector out, ix_coh, obs_local, rr, jhat_rr
    real matrix hit
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imp_cell_score_d.", 498)
    B = rows(f_rows)
    out = J(B, 1, 0)
    ci = model.cells[cid].cohort_index
    ix_coh = _cdid_imp_index(model.unit_cohort_index[f_rows] :== ci)
    ix_coh = _cdid_imp_keep_ix(model, model.cells[cid].recipe, f_rows, ix_coh)
    if (rows(ix_coh) == 0) return(out)
    pid = model.cells[cid].pair_id
    r = model.cells[cid].recipe
    gid = model.cells[cid].group_id
    hit = cdid_imputation_read_pair_chunk(model.reader, model.pop, model.pairs[pid], f_rows[ix_coh], 1)
    obs_local = _cdid_imp_index(hit[, 1] :== 1)
    if (rows(obs_local) == 0) return(out)
    rr = ix_coh[obs_local]
    jhat_rr = J(rows(rr), 1, model.cells[cid].alpha)
    if (model.groups[gid].p > 0) {
        jhat_rr = jhat_rr + (model.cells[cid].L :* _cdid_impd_z(rdr, r, f_rows[rr], model.groups[gid].kept, 1)) * model.groups[gid].beta
    }
    out[rr] = _cdid_impd_w(rdr, r, f_rows[rr], 1) :* (hit[obs_local, 2] :- jhat_rr)
    return(out)
}

/* Group slope score h_{i, beta_r} on a bounded chunk of F rows (dated reader). */
real matrix cdid_imp_group_score_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar gid,
    real colvector f_rows)
{
    struct _cdid_imp_ivec rowvector c_ix
    real scalar B, k_keep, n_c, a_lo, a_hi, loc, cid0, cid_next, s0, t0, any_coh, cid, ci, r
    real colvector cids, ix_coh, obs_local, rr, jhat_rr, s_rr
    real matrix h_beta, hit, d_rr, d_tilde, act_keys
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (gid < 1 | gid > model.n_groups | gid != floor(gid))
        _cdid_imp_fail("Invalid group index in cdid_imp_group_score_d.", 498)
    B = rows(f_rows)
    k_keep = model.groups[gid].p
    h_beta = J(B, k_keep, 0)
    if (k_keep == 0) return(h_beta)
    c_ix = _cdid_imp_bucket_coh(model.unit_cohort_index[f_rows], model.n_cohorts)
    cids = model.groups[gid].cell_ids
    n_c = rows(cids)
    if (n_c > 1 & model.groups[gid].pooling == "across_cohorts") {
        act_keys = J(n_c, 4, 0)
        for (loc=1; loc<=n_c; loc++) {
            cid = cids[loc]
            act_keys[loc, ] = (model.cells[cid].L, model.cells[cid].s_index, model.cells[cid].t_index, cid)
        }
        cids = cids[order(act_keys, (1, 2, 3, 4))]
    }
    a_lo = 1
    while (a_lo <= n_c) {
        cid0 = cids[a_lo]
        s0 = model.cells[cid0].s_index
        t0 = model.cells[cid0].t_index
        a_hi = a_lo
        while (a_hi < n_c) {
            cid_next = cids[a_hi + 1]
            if (model.cells[cid_next].s_index != s0 | model.cells[cid_next].t_index != t0) break
            a_hi = a_hi + 1
        }
        any_coh = 0
        for (loc=a_lo; loc<=a_hi; loc++) {
            if (rows(c_ix[model.cells[cids[loc]].cohort_index].rows) > 0) {
                any_coh = 1
                break
            }
        }
        if (any_coh) {
            hit = cdid_imputation_read_pair_chunk(model.reader, model.pop,
                model.pairs[model.cells[cid0].pair_id], f_rows, 1)
            if (any(hit[, 1] :== 1)) {
                for (loc=a_lo; loc<=a_hi; loc++) {
                    cid = cids[loc]
                    ci = model.cells[cid].cohort_index
                    ix_coh = _cdid_imp_keep_ix(model, model.cells[cid].recipe, f_rows, c_ix[ci].rows)
                    if (rows(ix_coh) == 0) continue
                    obs_local = _cdid_imp_index(hit[ix_coh, 1] :== 1)
                    if (rows(obs_local) == 0) continue
                    rr = ix_coh[obs_local]
                    r = model.cells[cid].recipe
                    d_rr = model.cells[cid].L :* _cdid_impd_z(rdr, r, f_rows[rr], model.groups[gid].kept, 1)
                    jhat_rr = model.cells[cid].alpha :+ d_rr * model.groups[gid].beta
                    s_rr = _cdid_impd_w(rdr, r, f_rows[rr], 1) :* (hit[rr, 2] :- jhat_rr)
                    d_tilde = (model.intercept ? d_rr :- model.cells[cid].mean_d : d_rr)
                    h_beta[rr, ] = h_beta[rr, ] + s_rr :* d_tilde
                }
            }
        }
        a_lo = a_hi + 1
    }
    return(h_beta)
}

/* Forward chunk evaluation of group slope IF (dated reader). */
real matrix cdid_imp_beta_if_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar gid,
    real colvector f_rows,
    | string scalar scale)
{
    real matrix h_beta, xi_beta
    string scalar sc
    sc = (args() < 5 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    h_beta = cdid_imp_group_score_d(model, rdr, gid, f_rows)
    if (model.groups[gid].p == 0) return(h_beta)
    xi_beta = _cdid_nuisance_gram(model.groups[gid].R, model.groups[gid].pivot, (model.N :* h_beta)')'
    if (sc == "core") xi_beta = model.lambda :* xi_beta
    return(xi_beta)
}

/* Forward chunk evaluation of cell intercept IF (dated reader). */
real colvector cdid_imp_alpha_if_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    real scalar cid,
    real colvector f_rows,
    | string scalar scale)
{
    real scalar gid
    real colvector s_c, xi_alpha
    real matrix xi_beta
    string scalar sc
    sc = (args() < 5 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imp_alpha_if_d.", 498)
    if (!model.intercept) return(J(rows(f_rows), 1, 0))
    gid = model.cells[cid].group_id
    s_c = cdid_imp_cell_score_d(model, rdr, cid, f_rows)
    xi_alpha = s_c :/ model.cells[cid].mass
    if (model.groups[gid].p > 0) {
        xi_beta = cdid_imp_beta_if_d(model, rdr, gid, f_rows, "frame")
        xi_alpha = xi_alpha - xi_beta * model.cells[cid].mean_d'
    }
    if (sc == "core") xi_alpha = model.lambda :* xi_alpha
    return(xi_alpha)
}

/* Sparse adjoint pullback on a bounded chunk of F rows using the scoped dated reader.
   Sums a recurring original unit's contributions across all active cells/recipes into
   its single original-F row before any product. */
real matrix cdid_imp_adjoint_chunk_d(
    struct cdid_imputation_model scalar model,
    struct cdid_imp_dated_reader scalar rdr,
    struct cdid_imputation_adjoint scalar adj,
    real colvector f_rows,
    | string scalar scale)
{
    struct _cdid_imp_ivec rowvector c_ix
    real scalar B, M, n_act, a_lo, a_hi, a_idx, cid0, cid_next, s0, t0, any_coh, cid, gid, ci, r
    real colvector ix_coh, obs_local, rr, jhat_rr, s_rr
    real matrix out, hit, d_rr, d_tilde, val
    string scalar sc
    sc = (args() < 5 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    _cdid_imp_validate_chunk_rows(model, f_rows)
    B = rows(f_rows)
    M = adj.n_contrasts
    out = J(B, M, 0)
    n_act = rows(adj.active_cells)
    if (n_act == 0) return(out)
    c_ix = _cdid_imp_bucket_coh(model.unit_cohort_index[f_rows], model.n_cohorts)
    a_lo = 1
    while (a_lo <= n_act) {
        cid0 = adj.active_cells[a_lo]
        s0 = model.cells[cid0].s_index
        t0 = model.cells[cid0].t_index
        a_hi = a_lo
        while (a_hi < n_act) {
            cid_next = adj.active_cells[a_hi + 1]
            if (model.cells[cid_next].s_index != s0 | model.cells[cid_next].t_index != t0) break
            a_hi = a_hi + 1
        }
        any_coh = 0
        for (a_idx=a_lo; a_idx<=a_hi; a_idx++) {
            if (rows(c_ix[model.cells[adj.active_cells[a_idx]].cohort_index].rows) > 0) {
                any_coh = 1
                break
            }
        }
        if (any_coh) {
            hit = cdid_imputation_read_pair_chunk(model.reader, model.pop,
                model.pairs[model.cells[cid0].pair_id], f_rows, 1)
            if (any(hit[, 1] :== 1)) {
                for (a_idx=a_lo; a_idx<=a_hi; a_idx++) {
                    cid = adj.active_cells[a_idx]
                    ci = model.cells[cid].cohort_index
                    ix_coh = _cdid_imp_keep_ix(model, model.cells[cid].recipe, f_rows, c_ix[ci].rows)
                    if (rows(ix_coh) == 0) continue
                    obs_local = _cdid_imp_index(hit[ix_coh, 1] :== 1)
                    if (rows(obs_local) == 0) continue
                    rr = ix_coh[obs_local]
                    gid = model.cells[cid].group_id
                    r = model.cells[cid].recipe
                    if (model.groups[gid].p > 0) {
                        d_rr = model.cells[cid].L :* _cdid_impd_z(rdr, r, f_rows[rr], model.groups[gid].kept, 1)
                        jhat_rr = model.cells[cid].alpha :+ d_rr * model.groups[gid].beta
                        d_tilde = (model.intercept ? d_rr :- model.cells[cid].mean_d : d_rr)
                        val = (d_tilde * adj.groups[gid].nu) :+ adj.tilde_lambda[cid, ]
                        s_rr = _cdid_impd_w(rdr, r, f_rows[rr], 1) :* (hit[rr, 2] :- jhat_rr)
                        out[rr, ] = out[rr, ] + s_rr :* val
                    }
                    else {
                        s_rr = _cdid_impd_w(rdr, r, f_rows[rr], 1) :* (hit[rr, 2] :- model.cells[cid].alpha)
                        out[rr, ] = out[rr, ] + s_rr * adj.tilde_lambda[cid, ]
                    }
                }
            }
        }
        a_lo = a_hi + 1
    }
    if (sc == "core") out = model.lambda :* out
    return(out)
}

end
