*! Native selection-probability and model-IF primitives for cdid; StataNow 19.5.
*! Native adaptation of David Benatia et al., cdid 1.1.0.0 (GPL-3),
*! commit :
*! R/chained_selection_models.R, R/chained_cp04_corollary1.R,
*! R/chained_cp04_public.R, and R/chained_selection_consumers.R.
version 19.5
mata:

struct cdid_sel_decl {
    real scalar enabled, has_xformla, has_rho, has_just
    string scalar selection, canonical_assumption, selection_universe
    string scalar selection_model, selection_pooling, justification
}

struct cdid_sel_design_spec {
    real scalar N, d, d_s, has_intercept, int_col
    real rowvector slope_cols, col_order, inv_col_order
    string scalar kind
}

struct _cdid_sel_ivec {
    real colvector v
}

struct _cdid_sel_mvec {
    real matrix M
}

struct cdid_sel_reader {
    string scalar kind, selection_model, selection_pooling, selection_universe
    string scalar justification
    real scalar N, n_U, n_A, n_cohorts, n_times, max_chunk_rows, has_pop
    struct cdid_sel_design_spec scalar spec
    real matrix Z, date_responses, pair_responses
    real colvector pair_s_idx, pair_t_idx
    real colvector cohorts, all_times, cohort_pos, U, A, B, w, rho
    real colvector mass_F, mass_U, repeat_mass, a_U
    struct cdid_population scalar pop
    real matrix pop_obs_mat
    real scalar pop_obs_ready
    struct _cdid_sel_ivec colvector coh_u_idx
    real colvector coh_u_ref
}

struct cdid_sel_cell_map {
    string scalar selection_model, selection_pooling
    real scalar n_pair_cohorts, n_cells, n_blocks, has_ref
    real matrix pair_cohorts, cells, blocks
    real colvector block_cell_ids
}

struct cdid_sel_block_fit {
    real scalar block_id, d, d_s, B_m, p_m, has_intercept
    real scalar n_fit, n_obs_rows, iterations, deviance, rank, rank_tol, cutoff, rcond_R
    string scalar kind, solver
    real colvector cell_ids, cell_n, cell_mass, cell_pi, cell_r, m_b, alpha, beta, theta
    real matrix z_bar, R, H_chol
    real rowvector pivot
    real matrix H_inv
    struct _cdid_sel_ivec colvector cell_u_idx, cell_s_u, cell_eta_u, cell_p_u
    real scalar cell_cache_ready, dup_of_block
}

struct cdid_sel_fit_result {
    string scalar selection_model, selection_pooling, selection_universe
    real scalar N, n_U, n_A, n_cells, n_blocks, n_pair_cohorts
    real scalar n_alpha, n_beta, total_params, max_chunk_rows, eps, maxit
    struct cdid_sel_design_spec scalar spec
    struct cdid_sel_cell_map scalar map
    struct cdid_sel_block_fit vector blocks
    real colvector cohorts, all_times, mass_F, mass_U, repeat_mass, a_U
    real colvector alpha_start, beta_start, block_p, cell_block, cell_local
    real colvector alpha, beta, theta
    real scalar fit_w_sum
}

struct cdid_sel_budget {
    real scalar N, n_U, n_pair_cohorts, n_cells, n_blocks
    real scalar d, d_s, max_block_cells, max_block_params, total_params
    real scalar max_chunk_rows, dense_block_cols, compact_state_doubles
    real scalar stream_workspace_doubles, work_units, max_work, within_work_cap
}

struct cdid_sel_pred {
    real scalar cohort, a_k, mass_U_k, cell_id, cell_s, cell_t
    real colvector eta, eta_s, eta_t, p_raw, p_s, p_t, p_pair, p_weight
    real colvector active, w, rho, sel_weight, score_factor
}

struct cdid_sel_probes {
    real scalar n_targets, n_pair_cohorts, d_s
    real colvector target_pc, target_cohort, target_rho
    real colvector pair_probabilities, multipliers, active_mask
    real colvector eta_s, eta_t, p_s, p_t
    real matrix target_z_slope
}

struct cdid_sel_up_map {
    string scalar selection_model, selection_pooling
    real scalar n_pair_cohorts, n_cells, n_blocks, n_cohorts, total_params, total_channels
    real colvector cohorts, a_U, mass_U, alpha_start, beta_start, block_p
    real matrix pc_info
}

struct cdid_sel_log_deriv {
    real scalar pc, cohort, cohort_pos, a_k, mass_U_k
    real scalar block_id, local_b, block_s, local_s, block_t, local_t, p_pc
    real colvector active, p_pair, p_weight, p_s, p_t, d_ak
    real matrix d_theta_local
}

struct cdid_sel_moments {
    real scalar n_moments, total_params, n_cohorts
    real matrix G_theta, G_ak
}

struct cdid_sel_if_chunk {
    real colvector unit_ids
    real matrix xi_alpha, xi_beta, theta_if, cohort_a_if
}

struct cdid_sel_table_out {
    real colvector alpha, beta, coefficients, cohort_a, mass_U
    real colvector pair_probabilities, multipliers
    real matrix theta_if, cohort_a_if, pair_prob_if, multiplier_if
    real matrix moment_if, loading_theta, loading_ak, covariance
    struct cdid_sel_fit_result scalar fit
    struct cdid_sel_budget scalar budget
}

void _cdid_sel_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_sel_index(real matrix pred)
{
    return(vec(selectindex(pred)))
}

real scalar _cdid_sel_is_pos_int(real scalar x, real scalar min_val)
{
    if (x >= . | x != floor(x) | x < min_val) return(0)
    return(1)
}

/* Vectorized logistic link inverse matching cdid_logit_values[, 1] exactly. */
real colvector _cdid_sel_logistic(real colvector eta)
{
    real colvector p, hi, o
    real scalar n

    if (missing(eta)) _cdid_sel_fail("input", "Nonfinite logit predictor.", 498)
    n = rows(eta)
    if (n == 0) return(J(0, 1, 0))
    if (!any(eta :> 40)) {
        o = exp(eta)
        return(o :/ (1 :+ o))
    }
    hi = (eta :> 40)
    o = exp(eta :* (1 :- hi))
    p = o :/ (1 :+ o)
    p = p :* (1 :- hi) :+ (1 :/ (1 :+ exp(-eta :* hi))) :* hi
    return(p)
}

/* Binary search bounds [p1, p2] in strictly ascending u_idx for [lo, hi]. */
void _cdid_sel_range_bounds(real colvector u_idx, real scalar lo, real scalar hi,
    real scalar p1, real scalar p2)
{
    real scalar n, l, r, mid

    n = rows(u_idx)
    if (n == 0 | lo > hi) {
        p1 = 1
        p2 = 0
        return
    }
    if (u_idx[n] < lo | u_idx[1] > hi) {
        p1 = 1
        p2 = 0
        return
    }
    l = 1
    r = n
    while (l < r) {
        mid = floor((l + r) / 2)
        if (u_idx[mid] < lo) l = mid + 1
        else r = mid
    }
    p1 = l
    if (u_idx[p1] > hi) {
        p2 = p1 - 1
        return
    }
    l = p1
    r = n
    while (l < r) {
        mid = floor((l + r + 1) / 2)
        if (u_idx[mid] <= hi) l = mid
        else r = mid - 1
    }
    p2 = l
}

void _cdid_sel_ensure_pop_obs(struct cdid_sel_reader scalar rdr)
{
    real colvector obs_r, obs_u, obs_d, hit
    real scalar t

    if (!rdr.has_pop | rdr.pop_obs_ready == 1) return
    rdr.pop_obs_ready = 1
    if (rdr.N < 1 | rdr.n_times < 1 | rdr.N * rdr.n_times > 8000000) {
        rdr.pop_obs_mat = J(0, 0, 0)
        return
    }
    rdr.pop_obs_mat = J(rdr.N, rdr.n_times, 0)
    obs_r = rdr.pop.observed_order
    if (rows(obs_r) == 0) return
    obs_u = rdr.pop.ids.row_unit[obs_r]
    obs_d = rdr.pop.date[obs_r]
    for (t = 1; t <= rdr.n_times; t = t + 1) {
        hit = _cdid_sel_index(obs_d :== rdr.all_times[t])
        if (rows(hit) > 0) rdr.pop_obs_mat[obs_u[hit], t] = J(rows(hit), 1, 1)
    }
}

real colvector _cdid_sel_coh_u(struct cdid_sel_reader scalar rdr, real scalar k)
{
    real scalar j, need_build

    need_build = 0
    if (rows(rdr.coh_u_idx) != rdr.n_cohorts | rows(rdr.coh_u_ref) != rdr.N) need_build = 1
    else if (rdr.coh_u_ref != rdr.U) need_build = 1

    if (need_build) {
        rdr.coh_u_idx = _cdid_sel_ivec(rdr.n_cohorts)
        for (j = 1; j <= rdr.n_cohorts; j = j + 1) {
            rdr.coh_u_idx[j].v = _cdid_sel_index((rdr.cohort_pos :== j) :& (rdr.U :== 1))
        }
        rdr.coh_u_ref = rdr.U
    }
    return(rdr.coh_u_idx[k].v)
}

/* Selection probability denominator floor 1e-6 and strict derivative mask.
   Distinct from treatment propensity upper cap 1 - 1e-6 in cdid_nuisance. */
real colvector cdid_sel_clip_floor(real colvector p)
{
    if (missing(p)) _cdid_sel_fail("input", "Nonfinite selection probability.", 498)
    return(rowmax((p, J(rows(p), 1, 1e-6))))
}

real colvector cdid_sel_floor_act(real colvector p)
{
    if (missing(p)) _cdid_sel_fail("input", "Nonfinite selection probability.", 498)
    return(p :> 1e-6)
}

/* Normalize and validate selection declaration. When selection == "none",
   any active selection option is refused, and otherwise an inactive descriptor
   is returned without inspecting covariates, weights, or responses. */
struct cdid_sel_decl scalar cdid_sel_norm_decl(string scalar selection,
    real scalar has_xformla, string scalar selection_universe,
    real scalar has_rho, real scalar has_just,
    string scalar selection_model, string scalar selection_pooling,
    string scalar justification)
{
    struct cdid_sel_decl scalar d
    string scalar sel, univ, mod, pool, just

    sel = strtrim(selection)
    univ = strtrim(selection_universe)
    mod = strtrim(selection_model)
    pool = strtrim(selection_pooling)
    just = strtrim(justification)
    if (univ == "") univ = "F"
    if (mod == "") mod = "pair_logit"
    if (pool == "") pool = "separate"

    if (sel != "none" & sel != "conditional_mtar" & sel != "corollary1")
        _cdid_sel_fail("input", "selection must be 'none', 'conditional_mtar', or 'corollary1'.", 498)
    if (univ != "F" & univ != "B" & univ != "A")
        _cdid_sel_fail("input", "selection_universe must be 'F', 'B', or 'A'.", 498)
    if (mod != "pair_logit" & mod != "independent_period_logit")
        _cdid_sel_fail("input", "selection_model must be 'pair_logit' or 'independent_period_logit'.", 498)
    if (pool != "separate" & pool != "within_cohort" & pool != "across_cohorts")
        _cdid_sel_fail("input", "selpooling() must be separate, within_cohort, or across_cohorts.", 498)

    d.selection = sel
    d.selection_universe = univ
    d.selection_model = mod
    d.selection_pooling = pool
    d.justification = just
    d.has_xformla = (has_xformla != 0)
    d.has_rho = (has_rho != 0)
    d.has_just = (has_just != 0 | just != "")

    if (sel == "none") {
        if (d.has_xformla | univ != "F" | d.has_rho | d.has_just |
            mod != "pair_logit" | pool != "separate")
            _cdid_sel_fail("unsupported", "selmodel(), selpooling(), selcovars(), seluniverse(), selrho() and seljustification() require selection(conditional_mtar).", 459)
        d.enabled = 0
        d.canonical_assumption = "none"
        return(d)
    }

    if (sel == "corollary1" & mod != "pair_logit")
        _cdid_sel_fail("unsupported", "selection(corollary1) is allowed only with selmodel(pair_logit); specify selection(conditional_mtar).", 459)

    d.enabled = 1
    d.canonical_assumption = "conditional_mtar"
    return(d)
}

/* Validate universe U, repeat-observer A, ever-observed B, and transport rho
   against the normalized selection declaration. */
real scalar cdid_sel_check_mode(struct cdid_sel_decl scalar decl,
    real scalar N, real colvector U, real colvector A, real colvector B,
    real colvector rho)
{
    real colvector u_idx

    if (!_cdid_sel_is_pos_int(N, 1))
        _cdid_sel_fail("input", "Full-frame N must be a positive integer.", 498)
    if (rows(U) != N | cols(U) != 1 | missing(U) | any((U :!= 0) :& (U :!= 1)))
        _cdid_sel_fail("input", "U must be an N x 1 binary indicator vector.", 498)
    if (rows(A) != N | cols(A) != 1 | missing(A) | any((A :!= 0) :& (A :!= 1)))
        _cdid_sel_fail("input", "A must be an N x 1 binary indicator vector.", 498)
    if (rows(B) != N | cols(B) != 1 | missing(B) | any((B :!= 0) :& (B :!= 1)))
        _cdid_sel_fail("input", "B must be an N x 1 binary indicator vector.", 498)
    if (any(A :> B))
        _cdid_sel_fail("input", "Repeat-observer indicator A must be a subset of ever-observed B.", 498)

    if (!decl.enabled) return(1)

    if (decl.selection_universe == "F" & any(U :!= 1))
        _cdid_sel_fail("input", "selection_universe = 'F' requires U == 1 on all units.", 498)
    if (decl.selection_universe == "B" & any(U :!= B))
        _cdid_sel_fail("input", "selection_universe = 'B' requires U == B on all units.", 498)
    if (decl.selection_universe == "A" & any(U :!= A))
        _cdid_sel_fail("input", "selection_universe = 'A' requires U == A on all units.", 498)

    u_idx = _cdid_sel_index(U :== 1)
    if (rows(u_idx) == 0)
        _cdid_sel_fail("support", "Selection estimation universe U is empty.", 459)

    if (decl.selection_model == "independent_period_logit") {
        if (decl.selection_universe != "F" | any(U :!= 1))
            _cdid_sel_fail("unsupported", "independent_period_logit requires selection_universe = 'F' because conditioning on A or B breaks marginal period independence.", 459)
        if (decl.has_rho | decl.has_just | strtrim(decl.justification) != "")
            _cdid_sel_fail("unsupported", "independent_period_logit requires selection_universe = 'F' with unit transport rho = 1.", 459)
        if (rows(rho) != N | cols(rho) != 1 | missing(rho) | any(rho :!= 1))
            _cdid_sel_fail("unsupported", "independent_period_logit requires rho = 1 on all units.", 459)
        return(1)
    }

    if (decl.selection_universe == "F") {
        if (decl.has_rho)
            _cdid_sel_fail("unsupported", "selrho() is not allowed with seluniverse(F).", 459)
        if (rows(rho) != N | cols(rho) != 1 | missing(rho) | any(rho :!= 1))
            _cdid_sel_fail("unsupported", "selection_universe = 'F' requires rho = 1 on all units.", 459)
    }
    else {
        if (!decl.has_just | strtrim(decl.justification) == "")
            _cdid_sel_fail("unsupported", "seluniverse(B) and seluniverse(A) require seljustification().", 459)
        if (any((A :== 1) :& (U :!= 1)))
            _cdid_sel_fail("support", "Selection estimation universe U must contain all repeat-observer units A.", 459)
        if (rows(rho) != N | cols(rho) != 1)
            _cdid_sel_fail("input", "rho must have length N.", 498)
        if (missing(rho[u_idx]) | any(rho[u_idx] :<= 0))
            _cdid_sel_fail("support", "selection_rho must be finite and strictly positive on U.", 459)
    }
    return(1)
}

/* Extract and validate selection design matrix specification on U rows only.
   Rows outside U are never required to be finite. */
struct cdid_sel_design_spec scalar cdid_sel_design_classify(real matrix Z,
    real colvector U, real scalar has_intercept, real scalar int_col)
{
    struct cdid_sel_design_spec scalar s
    real colvector u_idx
    real scalar n, d, j

    n = rows(U)
    if (n == 0 | cols(U) != 1 | missing(U) | any((U :!= 0) :& (U :!= 1)))
        _cdid_sel_fail("input", "U must be a nonempty binary column vector.", 498)
    if (rows(Z) != n)
        _cdid_sel_fail("input", "Selection design matrix Z must have N rows.", 498)
    if (has_intercept != 0 & has_intercept != 1)
        _cdid_sel_fail("input", "has_intercept must be 0 or 1.", 498)

    u_idx = _cdid_sel_index(U :== 1)
    if (rows(u_idx) == 0)
        _cdid_sel_fail("support", "Selection estimation universe U has no active units.", 459)

    d = cols(Z)
    if (d > 0) {
        if (missing(Z[u_idx, ]))
            _cdid_sel_fail("input", "Selection design matrix Z must be finite on all active U rows.", 498)
    }

    s.N = n
    s.d = d
    s.has_intercept = has_intercept

    if (d == 0) {
        if (has_intercept == 1)
            _cdid_sel_fail("input", "has_intercept = 1 requires at least one column in Z.", 498)
        s.d_s = 0
        s.int_col = 0
        s.slope_cols = J(1, 0, 0)
        s.col_order = J(1, 0, 0)
        s.inv_col_order = J(1, 0, 0)
        s.kind = "zero_col"
        return(s)
    }

    if (has_intercept == 1) {
        if (!_cdid_sel_is_pos_int(int_col, 1) | int_col > d)
            _cdid_sel_fail("input", "Invalid intercept column index int_col.", 498)
        if (any(Z[u_idx, int_col] :!= 1))
            _cdid_sel_fail("input", "Declared intercept column in Z must equal 1 on all active U rows.", 498)
        s.int_col = int_col
        s.slope_cols = J(1, 0, 0)
        for (j = 1; j <= d; j = j + 1) {
            if (j != int_col) s.slope_cols = (s.slope_cols, j)
        }
        s.d_s = cols(s.slope_cols)
        s.col_order = (int_col, s.slope_cols)
        s.inv_col_order = invorder(s.col_order)
        s.kind = (s.d_s > 0 ? "intercept_slope" : "intercept_only")
        return(s)
    }

    s.int_col = 0
    s.slope_cols = 1..d
    s.d_s = d
    s.col_order = 1..d
    s.inv_col_order = 1..d
    s.kind = "no_intercept"
    return(s)
}

/* Sparse date observation indicator I_{i,t} on requested original-frame units.
   Queries all ever-observed units in B (including singletons in B \ A), while
   all-missing-Y rows and frame-only units in F \ B remain 0. */
real colvector cdid_sel_date_rows(struct cdid_population scalar pop,
    real colvector unit_ids, real scalar t_val)
{
    real colvector ans, cal_hit
    real scalar n_req, j, u, lo, hi, mid, r

    if (t_val >= .)
        _cdid_sel_fail("input", "Calendar date t_val must be finite.", 498)
    cal_hit = _cdid_sel_index(pop.calendar :== t_val)
    if (rows(cal_hit) != 1)
        _cdid_sel_fail("input", "Requested date t_val is not in population calendar.", 498)

    n_req = rows(unit_ids)
    if (cols(unit_ids) != 1 | (n_req > 0 & (missing(unit_ids) |
        any(unit_ids :< 1) | any(unit_ids :> pop.N) | any(unit_ids :!= floor(unit_ids)))))
        _cdid_sel_fail("input", "unit_ids must be valid 1..N original-frame indices.", 498)

    ans = J(n_req, 1, 0)
    for (j = 1; j <= n_req; j = j + 1) {
        u = unit_ids[j]
        if (pop.B[u] == 0 | pop.row_count[u] == 0) continue
        lo = pop.row_start[u]
        hi = lo + pop.row_count[u] - 1
        while (lo <= hi) {
            mid = floor((lo + hi) / 2)
            r = pop.row_order[mid]
            if (pop.date[r] == t_val) {
                if (pop.y[r] < .) ans[j] = 1
                break
            }
            if (pop.date[r] < t_val) lo = mid + 1
            else hi = mid - 1
        }
    }
    return(ans)
}

/* Sparse pair observation indicator I_{i,s}*I_{i,t} on requested original-frame
   units. Queries repeat-observer units in A with s_val < t_val. */
real colvector cdid_sel_pair_rows(struct cdid_population scalar pop,
    real colvector unit_ids, real scalar s_val, real scalar t_val)
{
    real colvector ans, s_hit, t_hit
    real scalar n_req, j, u, lo, hi, mid, r, has_s, has_t

    if (s_val >= . | t_val >= . | s_val >= t_val)
        _cdid_sel_fail("input", "Pair dates must be finite with s_val < t_val.", 498)
    s_hit = _cdid_sel_index(pop.calendar :== s_val)
    t_hit = _cdid_sel_index(pop.calendar :== t_val)
    if (rows(s_hit) != 1 | rows(t_hit) != 1)
        _cdid_sel_fail("input", "Pair dates must belong to population calendar.", 498)

    n_req = rows(unit_ids)
    if (cols(unit_ids) != 1 | (n_req > 0 & (missing(unit_ids) |
        any(unit_ids :< 1) | any(unit_ids :> pop.N) | any(unit_ids :!= floor(unit_ids)))))
        _cdid_sel_fail("input", "unit_ids must be valid 1..N original-frame indices.", 498)

    ans = J(n_req, 1, 0)
    for (j = 1; j <= n_req; j = j + 1) {
        u = unit_ids[j]
        if (pop.A[u] == 0 | pop.row_count[u] < 2) continue
        has_s = 0
        lo = pop.row_start[u]
        hi = lo + pop.row_count[u] - 1
        while (lo <= hi) {
            mid = floor((lo + hi) / 2)
            r = pop.row_order[mid]
            if (pop.date[r] == s_val) {
                if (pop.y[r] < .) has_s = 1
                break
            }
            if (pop.date[r] < s_val) lo = mid + 1
            else hi = mid - 1
        }
        if (!has_s) continue
        has_t = 0
        lo = pop.row_start[u]
        hi = lo + pop.row_count[u] - 1
        while (lo <= hi) {
            mid = floor((lo + hi) / 2)
            r = pop.row_order[mid]
            if (pop.date[r] == t_val) {
                if (pop.y[r] < .) has_t = 1
                break
            }
            if (pop.date[r] < t_val) lo = mid + 1
            else hi = mid - 1
        }
        if (has_t) ans[j] = 1
    }
    return(ans)
}

/* Exact label lookup for an explicitly label-valued probe vector. Partition
   inputs use cohort positions; probe mode is decided for the entire vector. */
real scalar _cdid_sel_cohort_pos(real scalar val, real colvector cohorts)
{
    real scalar k
    for(k=1;k<=rows(cohorts);k=k+1)if(cohorts[k]==val)return(k)
    _cdid_sel_fail("input", "Cohort label not found in reader cohorts.", 498)
    return(0)
}

/* Build a full-F selection reader backed by cdid_population. Zero-observer
   cohorts in F are preserved in cohorts and only refused if requested by an
   active selection cell. */
struct cdid_sel_reader scalar cdid_sel_reader_pop(struct cdid_population scalar pop,
    real matrix Z, real scalar has_intercept, real scalar int_col,
    real colvector w, struct cdid_sel_decl scalar decl,
    real colvector rho, real scalar max_chunk_rows, | real colvector keep)
{
    struct cdid_sel_reader scalar rdr
    real colvector U, known_idx, u_idx, fk_idx, uk_idx, ak_idx, K
    real scalar k

    if (!decl.enabled)
        _cdid_sel_fail("input", "Cannot build active selection reader when selection is 'none'.", 498)
    if (!_cdid_sel_is_pos_int(max_chunk_rows, 1))
        _cdid_sel_fail("input", "max_chunk_rows must be a positive integer.", 498)
    /* keep: units eligible in a restricted reference eligibility context; the reader
       reads only those units (weights and designs of the others are never used). */
    K = J(0, 1, .)
    if (args() >= 9) if (rows(keep) == pop.N) K = (keep :!= 0)
    if (rows(w) != pop.N | cols(w) != 1)
        _cdid_sel_fail("input", "Sampling weights w must be finite and strictly positive on all N units.", 498)
    if (rows(K)) {
        if (missing(select(w, K)) | any(select(w, K) :<= 0))
            _cdid_sel_fail("input", "Sampling weights w must be finite and strictly positive on all N units.", 498)
    }
    else if (missing(w) | any(w :<= 0))
        _cdid_sel_fail("input", "Sampling weights w must be finite and strictly positive on all N units.", 498)

    if (decl.selection_universe == "F") U = J(pop.N, 1, 1)
    else if (decl.selection_universe == "B") U = pop.B
    else U = pop.A

    (void) cdid_sel_check_mode(decl, pop.N, U, pop.A, pop.B, rho)
    if (rows(K)) U = U :* K
    rdr.spec = cdid_sel_design_classify(Z, U, has_intercept, int_col)

    u_idx = _cdid_sel_index(U :== 1)
    if (any(pop.cohort_known[u_idx] :== 0))
        _cdid_sel_fail("support", "Selection estimation universe U contains units with unknown cohort.", 459)

    known_idx = _cdid_sel_index(pop.cohort_known :== 1)
    if (rows(known_idx) == 0)
        _cdid_sel_fail("support", "No units with known cohort in population.", 459)

    rdr.cohorts = uniqrows(sort(pop.cohort[known_idx], 1))
    rdr.n_cohorts = rows(rdr.cohorts)
    rdr.cohort_pos = J(pop.N, 1, 0)
    for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
        fk_idx = _cdid_sel_index((pop.cohort_known :!= 0) :& (pop.cohort :== rdr.cohorts[k]))
        if (rows(fk_idx) > 0) rdr.cohort_pos[fk_idx] = J(rows(fk_idx), 1, k)
    }
    if (any(rdr.cohort_pos[u_idx] :== 0))
        _cdid_sel_fail("support", "Active U units could not be mapped to known cohorts.", 459)

    rdr.mass_F = J(rdr.n_cohorts, 1, 0)
    rdr.mass_U = J(rdr.n_cohorts, 1, 0)
    rdr.repeat_mass = J(rdr.n_cohorts, 1, 0)
    rdr.a_U = J(rdr.n_cohorts, 1, 0)
    rdr.coh_u_idx = _cdid_sel_ivec(rdr.n_cohorts)
    for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
        fk_idx = _cdid_sel_index(rdr.cohort_pos :== k)
        if (rows(K) & rows(fk_idx)) fk_idx = select(fk_idx, K[fk_idx])
        if (rows(fk_idx) > 0) rdr.mass_F[k] = sum(w[fk_idx]) / pop.N
        uk_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1))
        rdr.coh_u_idx[k].v = uk_idx
        if (rows(uk_idx) > 0) rdr.mass_U[k] = sum(w[uk_idx]) / pop.N
        ak_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1) :& (pop.A :== 1))
        if (rows(ak_idx) > 0) rdr.repeat_mass[k] = sum(w[ak_idx]) / pop.N
        rdr.a_U[k] = (rdr.mass_U[k] > 0 ? rdr.repeat_mass[k] / rdr.mass_U[k] : 0)
    }

    rdr.kind = "population"
    rdr.selection_model = decl.selection_model
    rdr.selection_pooling = decl.selection_pooling
    rdr.selection_universe = decl.selection_universe
    rdr.justification = decl.justification
    rdr.N = pop.N
    rdr.n_U = rows(u_idx)
    rdr.n_A = pop.n_A
    rdr.all_times = pop.calendar
    rdr.n_times = rows(pop.calendar)
    rdr.max_chunk_rows = max_chunk_rows
    rdr.has_pop = 1
    rdr.Z = Z
    rdr.date_responses = J(0, 0, 0)
    rdr.pair_responses = J(0, 0, 0)
    rdr.pair_s_idx = J(0, 1, 0)
    rdr.pair_t_idx = J(0, 1, 0)
    rdr.U = U
    rdr.A = pop.A
    rdr.B = pop.B
    rdr.w = w
    rdr.rho = rho
    rdr.pop = pop
    rdr.pop_obs_ready = 0
    rdr.coh_u_ref = U
    _cdid_sel_ensure_pop_obs(rdr)
    return(rdr)
}

/* Build a full-F selection reader from explicit matrices/tables (mirrors R's
   evaluate_selection_table_if reader setup). */
struct cdid_sel_reader scalar cdid_sel_reader_table(real matrix Z,
    real scalar has_intercept, real scalar int_col,
    real colvector cohort_index, real colvector cohorts_in,
    real colvector all_times_in, real colvector w,
    real colvector A_in, real colvector B_in, real colvector U_in,
    real colvector rho_in, struct cdid_sel_decl scalar decl,
    real matrix date_responses, real matrix pair_responses,
    real colvector pair_s_idx, real colvector pair_t_idx,
    real scalar max_chunk_rows)
{
    struct cdid_sel_reader scalar rdr
    real colvector U, A, B, rho, cohorts, all_times, u_idx, fk_idx, uk_idx, ak_idx
    real scalar n, k, u, n_t, found

    if (!decl.enabled)
        _cdid_sel_fail("input", "Cannot build active selection reader when selection is 'none'.", 498)
    if (!_cdid_sel_is_pos_int(max_chunk_rows, 1))
        _cdid_sel_fail("input", "max_chunk_rows must be a positive integer.", 498)

    n = rows(Z)
    if (n < 1)
        _cdid_sel_fail("input", "Z must have at least one row.", 498)
    if (rows(w) != n | cols(w) != 1 | missing(w) | any(w :<= 0))
        _cdid_sel_fail("input", "w must be a finite positive N x 1 vector.", 498)
    if (rows(cohort_index) != n | cols(cohort_index) != 1)
        _cdid_sel_fail("input", "cohort_index must be an N x 1 vector.", 498)

    if (rows(A_in) == 0) {
        if (rows(date_responses) == n & cols(date_responses) > 0)
            A = (rowsum(date_responses :!= 0) :>= 2)
        else A = J(n, 1, 1)
    }
    else A = A_in

    if (rows(B_in) == 0) {
        if (rows(date_responses) == n & cols(date_responses) > 0)
            B = (rowsum(date_responses :!= 0) :>= 1) :| A
        else B = J(n, 1, 1)
    }
    else B = B_in

    if (rows(U_in) == 0) {
        if (decl.selection_universe == "F") U = J(n, 1, 1)
        else if (decl.selection_universe == "B") U = B
        else U = A
    }
    else U = U_in

    rho = (rows(rho_in) == 0 ? J(n, 1, 1) : rho_in)

    (void) cdid_sel_check_mode(decl, n, U, A, B, rho)
    rdr.spec = cdid_sel_design_classify(Z, U, has_intercept, int_col)

    u_idx = _cdid_sel_index(U :== 1)
    if (rows(cohorts_in) == 0) cohorts = uniqrows(sort(cohort_index[u_idx], 1))
    else cohorts = cohorts_in
    if (rows(cohorts) < 1)
        _cdid_sel_fail("input", "cohorts must be nonempty.", 498)

    rdr.cohorts = cohorts
    rdr.n_cohorts = rows(cohorts)
    rdr.cohort_pos = J(n, 1, 0)
    for (u = 1; u <= n; u = u + 1) {
        found = 0
        for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
            if (cohort_index[u] == cohorts[k]) {
                rdr.cohort_pos[u] = k
                found = 1
                break
            }
        }
        if (!found & U[u] == 1)
            _cdid_sel_fail("input", "cohort_index on U contains values not in cohorts.", 498)
    }

    if (rows(all_times_in) == 0) {
        n_t = (cols(date_responses) > 0 ? cols(date_responses) :
              (rows(pair_t_idx) > 0 ? max(pair_t_idx) : 1))
        all_times = (1::n_t)
    }
    else all_times = all_times_in
    n_t = rows(all_times)
    if (n_t < 1 | missing(all_times))
        _cdid_sel_fail("input", "all_times must be a nonempty finite vector.", 498)

    if (rows(date_responses) > 0) {
        if (rows(date_responses) != n | cols(date_responses) != n_t | missing(date_responses))
            _cdid_sel_fail("input", "date_responses must be a finite N x length(all_times) matrix.", 498)
    }
    if (rows(pair_responses) > 0) {
        if (rows(pair_responses) != n | cols(pair_responses) != rows(pair_s_idx) |
            rows(pair_s_idx) != rows(pair_t_idx) | missing(pair_responses))
            _cdid_sel_fail("input", "pair_responses dimensions must match N and pair_s_idx/pair_t_idx.", 498)
    }

    rdr.mass_F = J(rdr.n_cohorts, 1, 0)
    rdr.mass_U = J(rdr.n_cohorts, 1, 0)
    rdr.repeat_mass = J(rdr.n_cohorts, 1, 0)
    rdr.a_U = J(rdr.n_cohorts, 1, 0)
    rdr.coh_u_idx = _cdid_sel_ivec(rdr.n_cohorts)
    for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
        fk_idx = _cdid_sel_index(rdr.cohort_pos :== k)
        if (rows(fk_idx) > 0) rdr.mass_F[k] = sum(w[fk_idx]) / n
        uk_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1))
        rdr.coh_u_idx[k].v = uk_idx
        if (rows(uk_idx) > 0) rdr.mass_U[k] = sum(w[uk_idx]) / n
        ak_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1) :& (A :== 1))
        if (rows(ak_idx) > 0) rdr.repeat_mass[k] = sum(w[ak_idx]) / n
        rdr.a_U[k] = (rdr.mass_U[k] > 0 ? rdr.repeat_mass[k] / rdr.mass_U[k] : 0)
    }

    rdr.kind = "table"
    rdr.selection_model = decl.selection_model
    rdr.selection_pooling = decl.selection_pooling
    rdr.selection_universe = decl.selection_universe
    rdr.justification = decl.justification
    rdr.N = n
    rdr.n_U = rows(u_idx)
    rdr.n_A = sum(A :> 0)
    rdr.all_times = all_times
    rdr.n_times = n_t
    rdr.max_chunk_rows = max_chunk_rows
    rdr.has_pop = 0
    rdr.Z = Z
    rdr.date_responses = date_responses
    rdr.pair_responses = pair_responses
    rdr.pair_s_idx = pair_s_idx
    rdr.pair_t_idx = pair_t_idx
    rdr.U = U
    rdr.A = A
    rdr.B = B
    rdr.w = w
    rdr.rho = rho
    rdr.pop_obs_mat = J(0, 0, 0)
    rdr.pop_obs_ready = 1
    rdr.coh_u_ref = U
    return(rdr)
}

/* Canonicalize pair-cohort specifications, build unique cells (pair cells or
   date cells), and partition cells into estimation blocks by pooling mode.
   pair_specs has cohort POSITION, calendar indices, optional labels/reference ID. */
struct cdid_sel_cell_map scalar cdid_sel_partition(struct cdid_sel_reader scalar rdr,
    real matrix pair_specs)
{
    struct cdid_sel_cell_map scalar m
    real matrix raw_pc, ord_pc, pc_dedup, raw_dates, ord_dates, dates_dedup
    real colvector o, keep_mask, b_cells, uniq_keys, c_match
    real scalar n_in, i, k, s_cal, t_cal, s_val, t_val, ref_id, has_ref
    real scalar n_pc, n_c, c, b_id, pos, key_c

    n_in = rows(pair_specs)
    if (n_in < 1 | (cols(pair_specs) != 3 & cols(pair_specs) != 4 &
        cols(pair_specs) != 5 & cols(pair_specs) != 6))
        _cdid_sel_fail("input", "pair_specs must be a nonempty matrix with 3, 4, 5, or 6 columns.", 498)

    has_ref = (cols(pair_specs) == 4 | cols(pair_specs) == 6)
    raw_pc = J(n_in, 8, 0)
    for (i = 1; i <= n_in; i = i + 1) {
        k = pair_specs[i, 1]
        if (!_cdid_sel_is_pos_int(k, 1) | k > rdr.n_cohorts)
            _cdid_sel_fail("input", "pair_specs requires cohort positions in 1..n_cohorts.", 498)
        s_cal = pair_specs[i, 2]
        t_cal = pair_specs[i, 3]
        if (!_cdid_sel_is_pos_int(s_cal, 1) | !_cdid_sel_is_pos_int(t_cal, 1) |
            s_cal >= t_cal | t_cal > rdr.n_times)
            _cdid_sel_fail("input", "pair_specs calendar indices must satisfy 1 <= s_cal < t_cal <= T.", 498)
        s_val = rdr.all_times[s_cal]
        t_val = rdr.all_times[t_cal]
        if (cols(pair_specs) >= 5) {
            if (pair_specs[i, 4] != s_val | pair_specs[i, 5] != t_val)
                _cdid_sel_fail("input", "pair_specs dates disagree with all_times[s_cal/t_cal].", 498)
        }
        ref_id = 0
        if (cols(pair_specs) == 4) ref_id = pair_specs[i, 4]
        else if (cols(pair_specs) == 6) ref_id = pair_specs[i, 6]
        if (has_ref & !_cdid_sel_is_pos_int(ref_id, 1))
            _cdid_sel_fail("input", "pair_specs reference key must be a positive integer.", 498)
        raw_pc[i, ] = (k, s_cal, t_cal, s_val, t_val, 0, 0, ref_id)
    }

    o = order(raw_pc, (1, 2, 3, 8))
    ord_pc = raw_pc[o, ]
    keep_mask = J(n_in, 1, 1)
    for (i = 2; i <= n_in; i = i + 1) {
        if (ord_pc[i, 1] == ord_pc[i - 1, 1] & ord_pc[i, 2] == ord_pc[i - 1, 2] &
            ord_pc[i, 3] == ord_pc[i - 1, 3] & ord_pc[i, 8] == ord_pc[i - 1, 8])
            keep_mask[i] = 0
    }
    pc_dedup = select(ord_pc, keep_mask :== 1)
    n_pc = rows(pc_dedup)

    if (rdr.selection_model == "pair_logit") {
        n_c = n_pc
        m.cells = J(n_c, 9, 0)
        for (c = 1; c <= n_c; c = c + 1) {
            pc_dedup[c, 6] = c
            pc_dedup[c, 7] = c
            m.cells[c, ] = (c, 1, pc_dedup[c, 1], 0, pc_dedup[c, 2], pc_dedup[c, 3], 0, 0, pc_dedup[c, 8])
        }
    }
    else {
        raw_dates = (pc_dedup[, (1, 2, 8)]) \ (pc_dedup[, (1, 3, 8)])
        o = order(raw_dates, (1, 2, 3))
        ord_dates = raw_dates[o, ]
        keep_mask = J(rows(ord_dates), 1, 1)
        for (i = 2; i <= rows(ord_dates); i = i + 1) {
            if (ord_dates[i, 1] == ord_dates[i - 1, 1] &
                ord_dates[i, 2] == ord_dates[i - 1, 2] &
                ord_dates[i, 3] == ord_dates[i - 1, 3])
                keep_mask[i] = 0
        }
        dates_dedup = select(ord_dates, keep_mask :== 1)
        n_c = rows(dates_dedup)
        m.cells = J(n_c, 9, 0)
        for (c = 1; c <= n_c; c = c + 1) {
            m.cells[c, ] = (c, 2, dates_dedup[c, 1], dates_dedup[c, 2], 0, 0, 0, 0, dates_dedup[c, 3])
        }
        for (i = 1; i <= n_pc; i = i + 1) {
            c_match = _cdid_sel_index((m.cells[, 3] :== pc_dedup[i, 1]) :&
                                      (m.cells[, 4] :== pc_dedup[i, 2]) :&
                                      (m.cells[, 9] :== pc_dedup[i, 8]))
            pc_dedup[i, 6] = c_match[1]
            c_match = _cdid_sel_index((m.cells[, 3] :== pc_dedup[i, 1]) :&
                                      (m.cells[, 4] :== pc_dedup[i, 3]) :&
                                      (m.cells[, 9] :== pc_dedup[i, 8]))
            pc_dedup[i, 7] = c_match[1]
        }
    }

    if (rdr.selection_pooling == "separate") {
        m.n_blocks = n_c
        m.blocks = J(n_c, 5, 0)
        m.block_cell_ids = (1::n_c)
        for (c = 1; c <= n_c; c = c + 1) {
            m.cells[c, 7] = c
            m.cells[c, 8] = 1
            m.blocks[c, ] = (c, c, 1, m.cells[c, 3], m.cells[c, 9])
        }
    }
    else if (rdr.selection_pooling == "within_cohort") {
        raw_dates = uniqrows(sort(m.cells[, (3, 9)], (1, 2)))
        m.n_blocks = rows(raw_dates)
        m.blocks = J(m.n_blocks, 5, 0)
        m.block_cell_ids = J(n_c, 1, 0)
        pos = 1
        for (b_id = 1; b_id <= m.n_blocks; b_id = b_id + 1) {
            b_cells = _cdid_sel_index((m.cells[, 3] :== raw_dates[b_id, 1]) :&
                                      (m.cells[, 9] :== raw_dates[b_id, 2]))
            m.blocks[b_id, ] = (b_id, pos, rows(b_cells), raw_dates[b_id, 1], raw_dates[b_id, 2])
            for (i = 1; i <= rows(b_cells); i = i + 1) {
                c = b_cells[i]
                m.cells[c, 7] = b_id
                m.cells[c, 8] = i
                m.block_cell_ids[pos] = c
                pos = pos + 1
            }
        }
    }
    else {
        if (!has_ref) {
            m.n_blocks = 1
            m.blocks = (1, 1, n_c, 0, 0)
            m.block_cell_ids = (1::n_c)
            for (c = 1; c <= n_c; c = c + 1) {
                m.cells[c, 7] = 1
                m.cells[c, 8] = c
            }
        }
        else {
            uniq_keys = uniqrows(sort(m.cells[, 9], 1))
            m.n_blocks = rows(uniq_keys)
            m.blocks = J(m.n_blocks, 5, 0)
            m.block_cell_ids = J(n_c, 1, 0)
            pos = 1
            for (b_id = 1; b_id <= m.n_blocks; b_id = b_id + 1) {
                key_c = uniq_keys[b_id]
                b_cells = _cdid_sel_index(m.cells[, 9] :== key_c)
                m.blocks[b_id, ] = (b_id, pos, rows(b_cells), 0, key_c)
                for (i = 1; i <= rows(b_cells); i = i + 1) {
                    c = b_cells[i]
                    m.cells[c, 7] = b_id
                    m.cells[c, 8] = i
                    m.block_cell_ids[pos] = c
                    pos = pos + 1
                }
            }
        }
    }

    m.selection_model = rdr.selection_model
    m.selection_pooling = rdr.selection_pooling
    m.n_pair_cohorts = n_pc
    m.n_cells = n_c
    m.has_ref = has_ref
    m.pair_cohorts = pc_dedup
    return(m)
}

/* Query cell binary response vector y on requested original-frame unit_ids. */
real colvector _cdid_sel_cell_y(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real scalar cell_id,
    real colvector unit_ids)
{
    real scalar kind_code, cal_idx, s_cal, t_cal
    real colvector hit

    kind_code = map.cells[cell_id, 2]
    cal_idx = map.cells[cell_id, 4]
    s_cal = map.cells[cell_id, 5]
    t_cal = map.cells[cell_id, 6]

    if (rdr.has_pop) {
        if (rows(unit_ids) == 0) return(J(0, 1, 0))
        _cdid_sel_ensure_pop_obs(rdr)
        if (rows(rdr.pop_obs_mat) == rdr.N) {
            if (kind_code == 2)
                return(rdr.pop_obs_mat[unit_ids, cal_idx] :& (rdr.pop.B[unit_ids] :!= 0) :& (rdr.pop.row_count[unit_ids] :> 0))
            return(rdr.pop_obs_mat[unit_ids, s_cal] :& rdr.pop_obs_mat[unit_ids, t_cal] :&
                (rdr.pop.A[unit_ids] :!= 0) :& (rdr.pop.row_count[unit_ids] :>= 2))
        }
        if (kind_code == 2)
            return(cdid_sel_date_rows(rdr.pop, unit_ids, rdr.all_times[cal_idx]))
        return(cdid_sel_pair_rows(rdr.pop, unit_ids, rdr.all_times[s_cal], rdr.all_times[t_cal]))
    }

    if (kind_code == 2) {
        if (cols(rdr.date_responses) < cal_idx)
            _cdid_sel_fail("input", "Table reader missing date_responses column for date cell.", 498)
        return(rdr.date_responses[unit_ids, cal_idx])
    }

    if (cols(rdr.pair_responses) > 0) {
        hit = _cdid_sel_index((rdr.pair_s_idx :== s_cal) :& (rdr.pair_t_idx :== t_cal))
        if (rows(hit) >= 1) return(rdr.pair_responses[unit_ids, hit[1]])
    }
    if (cols(rdr.date_responses) >= t_cal)
        return((rdr.date_responses[unit_ids, s_cal] :!= 0) :& (rdr.date_responses[unit_ids, t_cal] :!= 0))

    _cdid_sel_fail("input", "Table reader missing pair response for requested pair cell.", 498)
    return(J(0, 1, 0))
}

/* Evaluate work and memory budget for the selection cell map and reader.
   When max_work < ., enforces work_units <= max_work with r(459). */
struct cdid_sel_budget scalar cdid_sel_budget_check(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real scalar max_work)
{
    struct cdid_sel_budget scalar b
    real scalar m, b_m, p_m, d_s, has_int

    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept
    b.N = rdr.N
    b.n_U = rdr.n_U
    b.n_pair_cohorts = map.n_pair_cohorts
    b.n_cells = map.n_cells
    b.n_blocks = map.n_blocks
    b.d = rdr.spec.d
    b.d_s = d_s
    b.max_chunk_rows = rdr.max_chunk_rows
    b.max_block_cells = max(map.blocks[, 3])
    b.max_block_params = (has_int ? b.max_block_cells + d_s : d_s)

    b.total_params = 0
    b.compact_state_doubles = 0
    for (m = 1; m <= map.n_blocks; m = m + 1) {
        b_m = map.blocks[m, 3]
        p_m = (has_int ? b_m + d_s : d_s)
        b.total_params = b.total_params + p_m
        b.compact_state_doubles = b.compact_state_doubles + p_m + b_m * (1 + d_s) + d_s * d_s
    }
    b.dense_block_cols = d_s
    b.stream_workspace_doubles = rdr.max_chunk_rows * (d_s + 4) + (d_s + 1) * (d_s + 1)
    b.work_units = map.n_cells * rdr.n_U * max((1, rdr.spec.d))
    b.max_work = max_work
    b.within_work_cap = 1
    if (max_work < .) {
        if (!_cdid_sel_is_pos_int(max_work, 1))
            _cdid_sel_fail("input", "max_work must be a positive integer when supplied.", 498)
        b.within_work_cap = (b.work_units <= max_work)
        if (!b.within_work_cap)
            _cdid_sel_fail("bound", "Selection estimation work exceeds caller max_work cap.", 459)
    }
    return(b)
}

/* Compress a stacked matrix of chunk QR factors into a single upper-triangular
   matrix of at most cols(R_stack) rows. */
real matrix _cdid_sel_qr_compress(real matrix R_stack)
{
    real matrix Q, R
    real rowvector tau

    if (rows(R_stack) <= cols(R_stack)) return(R_stack)
    hqrd(R_stack, Q = ., tau = ., R = .)
    return(R)
}

/* Compute block deviance at current (alpha, beta) across all cells and chunks
   in block, tracking numerical boundary hits. */
real scalar _cdid_sel_block_dev(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real colvector cell_ids,
    real colvector alpha, real colvector beta, real scalar eps_prob,
    real scalar boundary_hit)
{
    real colvector u_idx, ch_idx, y_fit, w_fit, eta, mu
    real matrix D_ch
    real scalar B_m, d_s, has_int, dev, b, c, k, n_b, s_pos, e_pos, n_ch

    B_m = rows(cell_ids)
    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept
    dev = 0
    for (b = 1; b <= B_m; b = b + 1) {
        c = cell_ids[b]
        k = map.cells[c, 3]
        u_idx = _cdid_sel_coh_u(rdr, k)
        n_b = rows(u_idx)
        for (s_pos = 1; s_pos <= n_b; s_pos = s_pos + rdr.max_chunk_rows) {
            e_pos = min((n_b, s_pos + rdr.max_chunk_rows - 1))
            ch_idx = u_idx[|s_pos \ e_pos|]
            n_ch = rows(ch_idx)
            D_ch = (d_s > 0 ? rdr.Z[ch_idx, rdr.spec.slope_cols] : J(n_ch, 0, 0))
            y_fit = _cdid_sel_cell_y(rdr, map, c, ch_idx)
            w_fit = rdr.w[ch_idx]
            eta = (has_int ? (d_s > 0 ? alpha[b] :+ (D_ch * beta) : J(n_ch, 1, alpha[b])) :
                             (d_s > 0 ? D_ch * beta : J(n_ch, 1, 0)))
            if (missing(eta)) return(.)
            if (any(abs(eta) :> 30)) boundary_hit = 1
            mu = _cdid_sel_logistic(eta)
            if (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob)) {
                boundary_hit = 1
                mu = rowmin((rowmax((mu, J(n_ch, 1, eps_prob))), J(n_ch, 1, 1 - eps_prob)))
            }
            dev = dev + _cdid_logit_deviance(y_fit, mu, w_fit)
        }
    }
    return(dev)
}

/* Core IRLS and final-information solver for a pre-validated separate pair_logit
   block on cohort k with active U rows u_idx, design X_raw, and weights w_fit. */
struct cdid_sel_block_fit scalar _cdid_sel_fit_sep_pair_core(struct cdid_sel_reader scalar rdr,
    real scalar block_id, real colvector cell_ids, real scalar k,
    real colvector u_idx, real matrix X_raw, real colvector w_fit,
    real scalar ws_k, real colvector y_fit, real scalar eps, real scalar maxit)
{
    struct cdid_sel_block_fit scalar out
    real colvector wt, b_vec, bn, eta, p, sw, z, qy, signed_x
    real rowvector ord, tau, pivot, sx_min, sx_max
    real matrix D_ch, H, L_chol, A, Q, R, z_sol, signed_X
    real scalar d, d_s, has_int, n_fit, j, h, it, rk, dev, olddev, boundary, converged, qtol, wt_sum

    d = rdr.spec.d
    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept
    n_fit = rows(u_idx)
    if (missing(y_fit) | any((y_fit :!= 0) :& (y_fit :!= 1)))
        _cdid_sel_fail("input", "Invalid X, y, or w in separate pair selection block.", 498)
    if (all(y_fit :== 0) | all(y_fit :== 1))
        _cdid_sel_fail("rank", "No finite logit root: response is constant on the estimation subset.", 459)

    if (d > 0) {
        signed_X = (2 :* y_fit :- 1) :* X_raw
        sx_min = colmin(signed_X)
        sx_max = colmax(signed_X)
        if (any(((sx_min :>= 0) :& (sx_max :> 0)) :| ((sx_max :<= 0) :& (sx_min :< 0))))
            _cdid_sel_fail("rank", "No finite logit root: signed-coordinate separation witness.", 459)
    }

    b_vec = J(d, 1, 0)
    eta = J(n_fit, 1, 0)
    p = J(n_fit, 1, 0.5)
    olddev = _cdid_logit_deviance(y_fit, p, w_fit)
    converged = 0
    boundary = 0
    qtol = min((1e-7, eps / 1000))
    for (it = 1; it <= maxit; it = it + 1) {
        wt = w_fit :* p :* (1 :- p)
        H = cross(X_raw, wt, X_raw)
        L_chol = cholesky(H)
        bn = J(0, 1, .)
        if (!missing(L_chol)) {
            if (all(diagonal(L_chol) :> 100 * qtol :* sqrt(diagonal(H)))) {
                z = cross(X_raw, wt :* eta :+ w_fit :* (y_fit :- p))
                qy = solvelower(L_chol, z, rk = ., 0)
                if (rk == d & !missing(qy)) {
                    bn = solveupper(L_chol', qy, rk = ., 0)
                    if (rk != d | missing(bn)) bn = J(0, 1, .)
                }
            }
        }
        if (rows(bn) != d) {
            sw = sqrt(wt)
            A = X_raw :* sw
            hqrd(A, Q = ., tau = ., R = .)
            if (rows(R) != d | missing(R) | any(abs(diagonal(R)) :<= qtol :* sqrt(colsum(A :^ 2))'))
                _cdid_sel_fail("rank", "Covariates in the selection model are collinear on its estimation sample; remove redundant covariates.", 459)
            z = (eta :+ (y_fit :- p) :/ (p :* (1 :- p))) :* sw
            qy = hqrdmultq1t(Q, tau, z)
            bn = solveupper(R, qy, rk = ., 0)
            if (rk != d | missing(bn))
                _cdid_sel_fail("numerical", "Nonfinite logit IRLS step.", 504)
        }
        eta = X_raw * bn
        p = _cdid_sel_logistic(eta)
        dev = _cdid_logit_deviance(y_fit, p, w_fit)
        h = 0
        while (dev >= . | any(p :<= 0) | any(p :>= 1)) {
            boundary = 1
            h = h + 1
            if (h > maxit)
                _cdid_sel_fail("boundary", "Logit boundary step-halving failed.", 430)
            bn = (bn + b_vec) / 2
            eta = X_raw * bn
            p = _cdid_sel_logistic(eta)
            dev = _cdid_logit_deviance(y_fit, p, w_fit)
        }
        b_vec = bn
        if (abs(dev - olddev) / (0.1 + abs(dev)) < eps) {
            converged = 1
            break
        }
        olddev = dev
    }
    if (boundary)
        _cdid_sel_fail("boundary", "Logit encountered numerical boundary.", 430)
    if (!converged)
        _cdid_sel_fail("convergence", "Logit did not converge under supplied controls.", 430)
    if (any(abs(eta) :> 30) | any(p :<= 10 * epsilon(1)) | any(p :>= 1 - 10 * epsilon(1)))
        _cdid_sel_fail("boundary", "Logit fitted probabilities reached numerical boundary.", 430)

    wt = w_fit :* p :* (1 :- p)
    H = quadcross(X_raw, wt, X_raw) / rdr.N
    L_chol = cholesky(H)
    if (missing(L_chol) | min(diagonal(L_chol)) <= 0)
        _cdid_sel_fail("rank", "Selection logit information matrix is not positive definite.", 459)
    R = sqrt(rdr.N) :* L_chol'
    pivot = 1..d
    out.rcond_R = _cdid_nuisance_rcond(R)
    if (missing(H) | missing(R) | out.rcond_R <= 0)
        _cdid_sel_fail("rank", "Logit information is singular/nonfinite.", 459)

    ord = rdr.spec.col_order
    out.block_id = block_id
    out.d = d
    out.d_s = d_s
    out.B_m = 1
    out.p_m = d
    out.has_intercept = has_int
    out.kind = rdr.spec.kind
    out.cell_ids = cell_ids
    out.theta = b_vec[ord']
    if (has_int) {
        out.alpha = J(1, 1, b_vec[rdr.spec.int_col])
        out.beta = (d_s > 0 ? b_vec[rdr.spec.slope_cols'] : J(0, 1, 0))
    }
    else {
        out.alpha = J(0, 1, 0)
        out.beta = b_vec
    }
    out.H_chol = L_chol'
    z_sol = solvelower(L_chol, I(d), rk = ., 0)
    out.H_inv = (rk == d & !missing(z_sol) ? solveupper(out.H_chol, z_sol, rk = ., 0) : J(0, 0, 0))
    if (rk != d | missing(out.H_inv)) out.H_inv = J(0, 0, 0)
    out.R = R
    out.pivot = pivot
    wt_sum = sum(wt)
    out.cell_n = J(1, 1, n_fit)
    out.cell_mass = J(1, 1, ws_k / rdr.N)
    out.cell_pi = J(1, 1, sum(w_fit :* y_fit) / (rdr.mass_U[k] * rdr.N))
    out.cell_r = J(1, 1, out.cell_pi[1] / rdr.a_U[k])
    out.m_b = J(1, 1, wt_sum / rdr.N)
    out.z_bar = J(1, d_s, 0)
    if (has_int & d_s > 0) {
        D_ch = X_raw[, rdr.spec.slope_cols]
        out.z_bar = quadcross(wt, D_ch)' / wt_sum
    }
    out.n_fit = n_fit
    out.n_obs_rows = n_fit
    out.iterations = it
    out.deviance = dev
    out.rank = d
    out.rank_tol = qtol
    out.cutoff = 0
    out.solver = "cp04_separate_pair"
    out.cell_u_idx = _cdid_sel_ivec(1)
    out.cell_u_idx[1].v = u_idx
    out.cell_eta_u = _cdid_sel_ivec(1)
    out.cell_eta_u[1].v = eta
    out.cell_p_u = _cdid_sel_ivec(1)
    out.cell_p_u[1].v = p
    out.cell_s_u = _cdid_sel_ivec(1)
    out.cell_s_u[1].v = w_fit :* (y_fit :- p)
    out.cell_cache_ready = 1
    out.dup_of_block = 0
    return(out)
}

/* Separate positive-column pair_logit block fitter (cp04_logit / cdid_logit_fit). */
struct cdid_sel_block_fit scalar _cdid_sel_fit_sep_pair(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real scalar block_id,
    real colvector cell_ids, real scalar eps, real scalar maxit)
{
    real colvector u_idx, y_fit, w_fit
    real rowvector x_sum
    real matrix X_raw
    real scalar d, c, k, n_fit, j, h

    d = rdr.spec.d
    c = cell_ids[1]
    k = map.cells[c, 3]
    if (rdr.mass_U[k] <= 0 | rdr.a_U[k] <= 0)
        _cdid_sel_fail("support", "Selection block cohort has zero repeat-observer mass or U mass.", 459)
    u_idx = _cdid_sel_coh_u(rdr, k)
    n_fit = rows(u_idx)
    if (n_fit < d)
        _cdid_sel_fail("rank", "Selection block has fewer active rows than parameters.", 459)
    X_raw = rdr.Z[u_idx, ]
    w_fit = rdr.w[u_idx]
    y_fit = _cdid_sel_cell_y(rdr, map, c, u_idx)
    if (missing(X_raw) | missing(w_fit) | any(w_fit :<= 0) |
        missing(y_fit) | any((y_fit :!= 0) :& (y_fit :!= 1)))
        _cdid_sel_fail("input", "Invalid X, y, or w in separate pair selection block.", 498)
    if (all(y_fit :== 0) | all(y_fit :== 1))
        _cdid_sel_fail("rank", "No finite logit root: response is constant on the estimation subset.", 459)

    x_sum = colsum(X_raw)
    for (j = 1; j <= d; j = j + 1) {
        if (x_sum[j] == 0) {
            if (all(X_raw[, j] :== 0))
                _cdid_sel_fail("rank", "Logit zero predictor; model unchanged.", 459)
        }
        for (h = 1; h < j; h = h + 1) {
            if (x_sum[j] == x_sum[h]) {
                if (all(X_raw[, j] :== X_raw[, h]))
                    _cdid_sel_fail("rank", "A covariate in the selection model duplicates another column or the intercept (for example, a constant covariate); remove it.", 459)
            }
        }
    }
    return(_cdid_sel_fit_sep_pair_core(rdr, block_id, cell_ids, k, u_idx, X_raw, w_fit, sum(w_fit), y_fit, eps, maxit))
}

/* Build one IRLS or final-information WLS summary across all cells and chunks
   of a block using within-cell FWL and compact QR merging. */
real scalar _cdid_sel_block_wls(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real colvector cell_ids,
    real colvector alpha, real colvector beta, real scalar eps_prob,
    real scalar is_final, real scalar boundary_hit,
    real colvector m_b, real matrix z_bar, real colvector mean_y_vec,
    real matrix R_stack)
{
    real colvector u_idx, ch_idx, y_fit, w_fit, eta, mu, wt, z, z_c
    real colvector wd_b, d_min, d_max, mean_d
    real rowvector tau
    real matrix D_ch, D_c, A_ch, Q_ch, R_ch
    real rowvector d_row
    real scalar B_m, d_s, has_int, b, c, k, n_b, s_pos, e_pos, n_ch, ws_b, wy_b, mean_y, j

    B_m = rows(cell_ids)
    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept
    m_b = J(B_m, 1, 0)
    z_bar = J(B_m, d_s, 0)
    mean_y_vec = (has_int ? J(B_m, 1, 0) : J(0, 1, 0))
    R_stack = J(0, d_s + 1, 0)

    for (b = 1; b <= B_m; b = b + 1) {
        c = cell_ids[b]
        k = map.cells[c, 3]
        u_idx = _cdid_sel_coh_u(rdr, k)
        n_b = rows(u_idx)

        if (!has_int) {
            ws_b = 0
            for (s_pos = 1; s_pos <= n_b; s_pos = s_pos + rdr.max_chunk_rows) {
                e_pos = min((n_b, s_pos + rdr.max_chunk_rows - 1))
                ch_idx = u_idx[|s_pos \ e_pos|]
                n_ch = rows(ch_idx)
                D_ch = rdr.Z[ch_idx, rdr.spec.slope_cols]
                y_fit = _cdid_sel_cell_y(rdr, map, c, ch_idx)
                w_fit = rdr.w[ch_idx]
                eta = D_ch * beta
                if (missing(eta)) return(0)
                if (is_final) {
                    mu = _cdid_sel_logistic(eta)
                    if (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob) | any(mu :* (1 :- mu) :<= 0))
                        _cdid_sel_fail("boundary", "Selection block fitted probabilities reached boundary.", 430)
                }
                else {
                    if (any(abs(eta) :> 30)) boundary_hit = 1
                    mu = _cdid_sel_logistic(eta)
                    if (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob)) {
                        boundary_hit = 1
                        mu = rowmin((rowmax((mu, J(n_ch, 1, eps_prob))), J(n_ch, 1, 1 - eps_prob)))
                    }
                }
                wt = w_fit :* mu :* (1 :- mu)
                if (missing(wt) | any(wt :<= 0)) return(0)
                z = (is_final ? J(n_ch, 1, 0) : eta :+ (y_fit :- mu) :/ (mu :* (1 :- mu)))
                if (missing(z)) return(0)
                ws_b = ws_b + sum(wt)
                A_ch = (D_ch, z) :* sqrt(wt)
                hqrd(A_ch, Q_ch = ., tau = ., R_ch = .)
                R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
            }
            m_b[b] = ws_b / rdr.N
        }
        else {
            ws_b = 0
            wy_b = 0
            wd_b = J(d_s, 1, 0)
            d_min = J(d_s, 1, .)
            d_max = J(d_s, 1, .)
            for (s_pos = 1; s_pos <= n_b; s_pos = s_pos + rdr.max_chunk_rows) {
                e_pos = min((n_b, s_pos + rdr.max_chunk_rows - 1))
                ch_idx = u_idx[|s_pos \ e_pos|]
                n_ch = rows(ch_idx)
                D_ch = (d_s > 0 ? rdr.Z[ch_idx, rdr.spec.slope_cols] : J(n_ch, 0, 0))
                y_fit = _cdid_sel_cell_y(rdr, map, c, ch_idx)
                w_fit = rdr.w[ch_idx]
                eta = (d_s > 0 ? alpha[b] :+ (D_ch * beta) : J(n_ch, 1, alpha[b]))
                if (missing(eta)) return(0)
                if (is_final) {
                    mu = _cdid_sel_logistic(eta)
                    if (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob) | any(mu :* (1 :- mu) :<= 0))
                        _cdid_sel_fail("boundary", "Selection block fitted probabilities reached boundary.", 430)
                }
                else {
                    if (any(abs(eta) :> 30)) boundary_hit = 1
                    mu = _cdid_sel_logistic(eta)
                    if (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob)) {
                        boundary_hit = 1
                        mu = rowmin((rowmax((mu, J(n_ch, 1, eps_prob))), J(n_ch, 1, 1 - eps_prob)))
                    }
                }
                wt = w_fit :* mu :* (1 :- mu)
                if (missing(wt) | any(wt :<= 0)) return(0)
                z = (is_final ? J(n_ch, 1, 0) : eta :+ (y_fit :- mu) :/ (mu :* (1 :- mu)))
                if (missing(z)) return(0)
                ws_b = ws_b + sum(wt)
                wy_b = wy_b + sum(wt :* z)
                if (d_s > 0) {
                    wd_b = wd_b + quadcross(D_ch, wt)
                    for (j = 1; j <= d_s; j = j + 1) {
                        if (d_min[j] >= . | min(D_ch[, j]) < d_min[j]) d_min[j] = min(D_ch[, j])
                        if (d_max[j] >= . | max(D_ch[, j]) > d_max[j]) d_max[j] = max(D_ch[, j])
                    }
                }
            }
            if (ws_b <= 0 | ws_b >= .) return(0)
            m_b[b] = ws_b / rdr.N
            mean_y = wy_b / ws_b
            mean_y_vec[b] = mean_y
            mean_d = (d_s > 0 ? wd_b / ws_b : J(0, 1, 0))
            if (d_s > 0) {
                for (j = 1; j <= d_s; j = j + 1) {
                    if (d_min[j] == d_max[j]) mean_d[j] = d_min[j]
                }
                z_bar[b, ] = mean_d'
                d_row = mean_d'
                for (s_pos = 1; s_pos <= n_b; s_pos = s_pos + rdr.max_chunk_rows) {
                    e_pos = min((n_b, s_pos + rdr.max_chunk_rows - 1))
                    ch_idx = u_idx[|s_pos \ e_pos|]
                    n_ch = rows(ch_idx)
                    D_ch = rdr.Z[ch_idx, rdr.spec.slope_cols]
                    y_fit = _cdid_sel_cell_y(rdr, map, c, ch_idx)
                    w_fit = rdr.w[ch_idx]
                    eta = alpha[b] :+ (D_ch * beta)
                    mu = _cdid_sel_logistic(eta)
                    if (!is_final & (any(mu :<= eps_prob) | any(mu :>= 1 - eps_prob)))
                        mu = rowmin((rowmax((mu, J(n_ch, 1, eps_prob))), J(n_ch, 1, 1 - eps_prob)))
                    wt = w_fit :* mu :* (1 :- mu)
                    z = (is_final ? J(n_ch, 1, 0) : eta :+ (y_fit :- mu) :/ (mu :* (1 :- mu)))
                    D_c = D_ch :- d_row
                    for (j = 1; j <= d_s; j = j + 1) {
                        if (d_min[j] == d_max[j]) D_c[, j] = J(n_ch, 1, 0)
                    }
                    z_c = z :- mean_y
                    A_ch = (D_c, z_c) :* sqrt(wt)
                    hqrd(A_ch, Q_ch = ., tau = ., R_ch = .)
                    R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
                }
            }
        }
    }
    return(1)
}

/* Fit a single estimation block. Uses vectorized separate pair_logit solver for
   separate positive-column pair_logit, and within-cell FWL + compact QR merging
   for pooled blocks and independent_period_logit without any N x B_m dummy matrix. */
struct cdid_sel_block_fit scalar cdid_sel_fit_block(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real scalar block_id,
    real scalar eps, real scalar maxit)
{
    struct cdid_sel_block_fit scalar out
    struct _cdid_sel_ivec colvector u_arr, y_arr, w_arr
    struct _cdid_sel_mvec colvector D_arr
    real colvector cell_ids, u_idx, ch_idx, y_fit, w_fit, eta_b, mu_b, wt_b, z_b, z_c
    real colvector cell_n, cell_mass, cell_pi, cell_r, m_b, alpha, beta, alpha_new, beta_new
    real colvector sl_pos, sl_neg, s_coord, r_y, qy, d_abs, coef
    real rowvector tau, piv_step, mean_d_row
    real matrix D_ch, D_c, z_bar, R_stack, R_all, R_x, Q_x, R_step, A_all, A_b, R_ch, d_const_mat, d_min_mat
    real scalar start_pos, B_m, d, d_s, has_int, p_m, c, k, n_b, n_fit, n_obs_rows
    real scalar has_sep, b, s_pos, e_pos, c_pos, c_neg, j, iter, dev_old, dev
    real scalar boundary_hit, converged, rank_failed, rank_fail_iter, rank_tol, eps_prob
    real scalar ws_b, wy_b, mean_y, rk, max_d, cutoff, ok_wls, dummy_bh, use_fast, stack_direct, r_off

    if (!_cdid_sel_is_pos_int(block_id, 1) | block_id > map.n_blocks)
        _cdid_sel_fail("input", "Invalid block_id in cdid_sel_fit_block.", 498)
    if (eps <= 0 | eps >= . | !_cdid_sel_is_pos_int(maxit, 1))
        _cdid_sel_fail("input", "Invalid GLM controls eps/maxit.", 498)

    start_pos = map.blocks[block_id, 2]
    B_m = map.blocks[block_id, 3]
    cell_ids = map.block_cell_ids[|start_pos \ (start_pos + B_m - 1)|]
    d = rdr.spec.d
    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept
    p_m = (has_int ? B_m + d_s : d_s)

    if (B_m == 1 & rdr.selection_model == "pair_logit" & rdr.selection_pooling == "separate" & d > 0)
        return(_cdid_sel_fit_sep_pair(rdr, map, block_id, cell_ids, eps, maxit))

    out.block_id = block_id
    out.d = d
    out.d_s = d_s
    out.B_m = B_m
    out.p_m = p_m
    out.has_intercept = has_int
    out.kind = rdr.spec.kind
    out.cell_ids = cell_ids

    cell_n = J(B_m, 1, 0)
    cell_mass = J(B_m, 1, 0)
    cell_pi = J(B_m, 1, 0)
    cell_r = J(B_m, 1, 0)
    n_fit = 0
    for (b = 1; b <= B_m; b = b + 1) {
        c = cell_ids[b]
        k = map.cells[c, 3]
        if (rdr.mass_U[k] <= 0 | rdr.a_U[k] <= 0)
            _cdid_sel_fail("support", "Selection block cohort has zero repeat-observer mass or U mass.", 459)
        u_idx = _cdid_sel_coh_u(rdr, k)
        n_b = rows(u_idx)
        if (n_b < 1)
            _cdid_sel_fail("support", "Selection cell has zero estimation rows.", 459)
        cell_n[b] = n_b
        n_fit = n_fit + n_b
    }
    n_obs_rows = n_fit
    use_fast = (n_fit * max((1, d_s)) <= 2000000)

    has_sep = 0
    sl_pos = J(d_s, 1, 0)
    sl_neg = J(d_s, 1, 0)

    if (use_fast) {
        u_arr = _cdid_sel_ivec(B_m)
        y_arr = _cdid_sel_ivec(B_m)
        w_arr = _cdid_sel_ivec(B_m)
        D_arr = _cdid_sel_mvec(B_m)
        d_const_mat = J(B_m, d_s, 0)
        d_min_mat = J(B_m, d_s, 0)
        for (b = 1; b <= B_m; b = b + 1) {
            c = cell_ids[b]
            k = map.cells[c, 3]
            u_idx = _cdid_sel_coh_u(rdr, k)
            n_b = rows(u_idx)
            y_fit = _cdid_sel_cell_y(rdr, map, c, u_idx)
            w_fit = rdr.w[u_idx]
            if (missing(y_fit) | any((y_fit :!= 0) :& (y_fit :!= 1)))
                _cdid_sel_fail("input", "Selection cell response must be binary in {0, 1}.", 498)
            if (missing(w_fit) | any(w_fit :<= 0))
                _cdid_sel_fail("input", "Selection cell weights must be finite and positive.", 498)
            if (d_s > 0) {
                D_ch = rdr.Z[u_idx, rdr.spec.slope_cols]
                if (missing(D_ch))
                    _cdid_sel_fail("input", "Selection slope covariates must be finite on U.", 498)
                for (j = 1; j <= d_s; j = j + 1) {
                    s_coord = (2 :* y_fit :- 1) :* D_ch[, j]
                    if (any(s_coord :> 0)) sl_pos[j] = 1
                    if (any(s_coord :< 0)) sl_neg[j] = 1
                    d_min_mat[b, j] = min(D_ch[, j])
                    if (d_min_mat[b, j] == max(D_ch[, j])) d_const_mat[b, j] = 1
                }
            }
            else D_ch = J(n_b, 0, 0)
            c_pos = sum(y_fit :== 1)
            c_neg = n_b - c_pos
            ws_b = sum(w_fit)
            wy_b = sum(w_fit :* y_fit)
            cell_mass[b] = ws_b / rdr.N
            cell_pi[b] = wy_b / (rdr.mass_U[k] * rdr.N)
            cell_r[b] = cell_pi[b] / rdr.a_U[k]
            if (has_int & (c_pos == 0 | c_neg == 0)) has_sep = 1
            u_arr[b].v = u_idx
            y_arr[b].v = y_fit
            w_arr[b].v = w_fit
            D_arr[b].M = D_ch
        }
    }
    else {
        for (b = 1; b <= B_m; b = b + 1) {
            c = cell_ids[b]
            k = map.cells[c, 3]
            u_idx = _cdid_sel_coh_u(rdr, k)
            n_b = rows(u_idx)
            c_pos = 0
            c_neg = 0
            ws_b = 0
            wy_b = 0
            for (s_pos = 1; s_pos <= n_b; s_pos = s_pos + rdr.max_chunk_rows) {
                e_pos = min((n_b, s_pos + rdr.max_chunk_rows - 1))
                ch_idx = u_idx[|s_pos \ e_pos|]
                y_fit = _cdid_sel_cell_y(rdr, map, c, ch_idx)
                w_fit = rdr.w[ch_idx]
                if (missing(y_fit) | any((y_fit :!= 0) :& (y_fit :!= 1)))
                    _cdid_sel_fail("input", "Selection cell response must be binary in {0, 1}.", 498)
                if (missing(w_fit) | any(w_fit :<= 0))
                    _cdid_sel_fail("input", "Selection cell weights must be finite and positive.", 498)
                if (d_s > 0) {
                    D_ch = rdr.Z[ch_idx, rdr.spec.slope_cols]
                    if (missing(D_ch))
                        _cdid_sel_fail("input", "Selection slope covariates must be finite on U.", 498)
                    for (j = 1; j <= d_s; j = j + 1) {
                        s_coord = (2 :* y_fit :- 1) :* D_ch[, j]
                        if (any(s_coord :> 0)) sl_pos[j] = 1
                        if (any(s_coord :< 0)) sl_neg[j] = 1
                    }
                }
                c_pos = c_pos + sum(y_fit :== 1)
                c_neg = c_neg + sum(y_fit :== 0)
                ws_b = ws_b + sum(w_fit)
                wy_b = wy_b + sum(w_fit :* y_fit)
            }
            cell_mass[b] = ws_b / rdr.N
            cell_pi[b] = wy_b / (rdr.mass_U[k] * rdr.N)
            cell_r[b] = cell_pi[b] / rdr.a_U[k]
            if (has_int & (c_pos == 0 | c_neg == 0)) has_sep = 1
        }
    }

    if (n_fit < p_m)
        _cdid_sel_fail("rank", "Selection block has fewer active rows than parameters.", 459)
    if (!has_sep & d_s > 0) {
        for (j = 1; j <= d_s; j = j + 1) {
            if ((sl_pos[j] & !sl_neg[j]) | (!sl_pos[j] & sl_neg[j])) {
                has_sep = 1
                break
            }
        }
    }

    out.cell_n = cell_n
    out.cell_mass = cell_mass
    out.cell_pi = cell_pi
    out.cell_r = cell_r
    out.n_fit = n_fit
    out.n_obs_rows = n_obs_rows

    if (d == 0) {
        out.m_b = 0.25 :* cell_mass
        out.z_bar = J(B_m, 0, 0)
        out.alpha = J(0, 1, 0)
        out.beta = J(0, 1, 0)
        out.theta = J(0, 1, 0)
        out.R = J(0, 0, 0)
        out.H_chol = J(0, 0, 0)
        out.H_inv = J(0, 0, 0)
        out.pivot = J(1, 0, 0)
        out.iterations = 0
        out.deviance = sum(cell_mass) * rdr.N * (-2 * ln(0.5))
        out.rank = 0
        out.rank_tol = 1e-7
        out.cutoff = 0
        out.rcond_R = 1
        out.solver = "zero_col"
        out.cell_cache_ready = 0
        return(out)
    }

    rank_tol = max((1e-12, min((1e-7, eps / 1000))))
    eps_prob = 10 * epsilon(1)
    alpha = (has_int ? J(B_m, 1, 0) : J(0, 1, 0))
    beta = J(d_s, 1, 0)
    boundary_hit = 0
    dev_old = (use_fast ? sum(cell_mass) * rdr.N * (-2 * ln(0.5)) :
        _cdid_sel_block_dev(rdr, map, cell_ids, alpha, beta, eps_prob, boundary_hit))
    converged = 0
    rank_failed = 0
    rank_fail_iter = 0
    R_step = J(0, 0, 0)
    piv_step = J(1, 0, 0)
    cutoff = 0
    stack_direct = (use_fast & d_s > 0 & n_fit <= 16384)

    for (iter = 1; iter <= maxit; iter = iter + 1) {
        if (use_fast) {
            ok_wls = 1
            m_b = J(B_m, 1, 0)
            z_bar = J(B_m, d_s, 0)
            alpha_new = (has_int ? J(B_m, 1, 0) : J(0, 1, 0))
            if (stack_direct) A_all = J(n_fit, d_s + 1, 0)
            else R_stack = J(0, d_s + 1, 0)
            r_off = 1
            for (b = 1; b <= B_m; b = b + 1) {
                n_b = cell_n[b]
                eta_b = (has_int ? (d_s > 0 ? alpha[b] :+ (D_arr[b].M * beta) : J(n_b, 1, alpha[b])) :
                                   D_arr[b].M * beta)
                if (missing(eta_b)) {
                    ok_wls = 0
                    break
                }
                if (any(abs(eta_b) :> 30)) boundary_hit = 1
                mu_b = _cdid_sel_logistic(eta_b)
                if (any(mu_b :<= eps_prob) | any(mu_b :>= 1 - eps_prob)) {
                    boundary_hit = 1
                    mu_b = rowmin((rowmax((mu_b, J(n_b, 1, eps_prob))), J(n_b, 1, 1 - eps_prob)))
                }
                wt_b = w_arr[b].v :* mu_b :* (1 :- mu_b)
                if (missing(wt_b) | any(wt_b :<= 0)) {
                    ok_wls = 0
                    break
                }
                z_b = eta_b :+ (y_arr[b].v :- mu_b) :/ (mu_b :* (1 :- mu_b))
                if (missing(z_b)) {
                    ok_wls = 0
                    break
                }
                ws_b = sum(wt_b)
                if (ws_b <= 0 | ws_b >= .) {
                    ok_wls = 0
                    break
                }
                m_b[b] = ws_b / rdr.N
                if (!has_int) {
                    A_b = (D_arr[b].M, z_b) :* sqrt(wt_b)
                    if (stack_direct) A_all[|r_off, 1 \ (r_off + n_b - 1), (d_s + 1)|] = A_b
                    else {
                        hqrd(A_b, Q_x = ., tau = ., R_ch = .)
                        R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
                    }
                }
                else {
                    mean_y = sum(wt_b :* z_b) / ws_b
                    alpha_new[b] = mean_y
                    if (d_s > 0) {
                        mean_d_row = quadcross(wt_b, D_arr[b].M) / ws_b
                        for (j = 1; j <= d_s; j = j + 1) {
                            if (d_const_mat[b, j]) mean_d_row[j] = d_min_mat[b, j]
                        }
                        z_bar[b, ] = mean_d_row
                        D_c = D_arr[b].M :- mean_d_row
                        for (j = 1; j <= d_s; j = j + 1) {
                            if (d_const_mat[b, j]) D_c[, j] = J(n_b, 1, 0)
                        }
                        z_c = z_b :- mean_y
                        A_b = (D_c, z_c) :* sqrt(wt_b)
                        if (stack_direct) A_all[|r_off, 1 \ (r_off + n_b - 1), (d_s + 1)|] = A_b
                        else {
                            hqrd(A_b, Q_x = ., tau = ., R_ch = .)
                            R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
                        }
                    }
                }
                r_off = r_off + n_b
            }
            if (ok_wls & stack_direct) {
                hqrd(A_all, Q_x = ., tau = ., R_all = .)
            }
            else if (ok_wls & d_s > 0) R_all = _cdid_sel_qr_compress(R_stack)
        }
        else {
            ok_wls = _cdid_sel_block_wls(rdr, map, cell_ids, alpha, beta, eps_prob, 0,
                boundary_hit, m_b = ., z_bar = ., alpha_new = ., R_stack = .)
            if (ok_wls & d_s > 0) R_all = _cdid_sel_qr_compress(R_stack)
        }

        if (!ok_wls) {
            if (iter == 1) {
                rank_failed = 1
                rank_fail_iter = 1
            }
            else boundary_hit = 1
            break
        }

        if (d_s == 0) {
            beta_new = J(0, 1, 0)
            R_step = J(0, 0, 0)
            piv_step = J(1, 0, 0)
            cutoff = 0
        }
        else {
            if (rows(R_all) < d_s) {
                rank_failed = 1
                rank_fail_iter = iter
                break
            }
            R_x = R_all[, 1..d_s]
            r_y = R_all[, d_s + 1]
            hqrdp(R_x, Q_x = ., tau = ., R_step = ., piv_step = .)
            if (rows(R_step) != d_s | missing(R_step)) {
                rank_failed = 1
                rank_fail_iter = iter
                break
            }
            d_abs = abs(diagonal(R_step))
            max_d = max(d_abs)
            cutoff = rank_tol * max_d
            if (missing(d_abs) | max_d <= 0 | any(d_abs :<= cutoff)) {
                rank_failed = 1
                rank_fail_iter = iter
                break
            }
            qy = hqrdmultq1t(Q_x, tau, r_y)[|1 \ d_s|]
            coef = solveupper(R_step, qy, rk = ., 0)
            if (rk != d_s | missing(coef)) {
                rank_failed = 1
                rank_fail_iter = iter
                break
            }
            beta_new = coef[invorder(piv_step)]
            if (has_int) {
                alpha_new = alpha_new - z_bar * beta_new
            }
        }

        if (missing(alpha_new) | missing(beta_new)) {
            if (has_sep) break
            _cdid_sel_fail("numerical", "Nonfinite logit coefficients during IRLS.", 504)
        }
        alpha = alpha_new
        beta = beta_new
        if (use_fast) {
            dev = 0
            for (b = 1; b <= B_m; b = b + 1) {
                n_b = cell_n[b]
                eta_b = (has_int ? (d_s > 0 ? alpha[b] :+ (D_arr[b].M * beta) : J(n_b, 1, alpha[b])) :
                                   D_arr[b].M * beta)
                if (missing(eta_b)) {
                    dev = .
                    break
                }
                if (any(abs(eta_b) :> 30)) boundary_hit = 1
                mu_b = _cdid_sel_logistic(eta_b)
                if (any(mu_b :<= eps_prob) | any(mu_b :>= 1 - eps_prob)) {
                    boundary_hit = 1
                    mu_b = rowmin((rowmax((mu_b, J(n_b, 1, eps_prob))), J(n_b, 1, 1 - eps_prob)))
                }
                dev = dev + _cdid_logit_deviance(y_arr[b].v, mu_b, w_arr[b].v)
            }
        }
        else dev = _cdid_sel_block_dev(rdr, map, cell_ids, alpha, beta, eps_prob, boundary_hit)
        if (dev >= .) {
            boundary_hit = 1
            break
        }
        if (abs(dev - dev_old) / (0.1 + abs(dev)) < eps) {
            converged = 1
            break
        }
        dev_old = dev
    }

    if (rank_failed & (!has_sep | rank_fail_iter == 1))
        _cdid_sel_fail("rank", "Selection block design is rank deficient.", 459)
    if (has_sep)
        _cdid_sel_fail("rank", "No finite selection logit root: signed-coordinate separation witness.", 459)
    if (rank_failed)
        _cdid_sel_fail("rank", "Selection block design is rank deficient.", 459)
    if (boundary_hit)
        _cdid_sel_fail("boundary", "Selection block IRLS reached numerical boundary.", 430)
    if (!converged)
        _cdid_sel_fail("convergence", "Selection block IRLS failed to converge.", 430)

    if (use_fast) {
        ok_wls = 1
        m_b = J(B_m, 1, 0)
        z_bar = J(B_m, d_s, 0)
        out.cell_u_idx = u_arr
        out.cell_eta_u = _cdid_sel_ivec(B_m)
        out.cell_p_u = _cdid_sel_ivec(B_m)
        out.cell_s_u = _cdid_sel_ivec(B_m)
        if (stack_direct) A_all = J(n_fit, d_s, 0)
        else R_stack = J(0, d_s, 0)
        r_off = 1
        for (b = 1; b <= B_m; b = b + 1) {
            n_b = cell_n[b]
            eta_b = (has_int ? (d_s > 0 ? alpha[b] :+ (D_arr[b].M * beta) : J(n_b, 1, alpha[b])) :
                               D_arr[b].M * beta)
            if (missing(eta_b)) {
                ok_wls = 0
                break
            }
            mu_b = _cdid_sel_logistic(eta_b)
            if (any(mu_b :<= eps_prob) | any(mu_b :>= 1 - eps_prob) | any(mu_b :* (1 :- mu_b) :<= 0))
                _cdid_sel_fail("boundary", "Selection block fitted probabilities reached boundary.", 430)
            wt_b = w_arr[b].v :* mu_b :* (1 :- mu_b)
            if (missing(wt_b) | any(wt_b :<= 0)) {
                ok_wls = 0
                break
            }
            ws_b = sum(wt_b)
            if (ws_b <= 0 | ws_b >= .) {
                ok_wls = 0
                break
            }
            m_b[b] = ws_b / rdr.N
            out.cell_eta_u[b].v = eta_b
            out.cell_p_u[b].v = mu_b
            out.cell_s_u[b].v = w_arr[b].v :* (y_arr[b].v :- mu_b)
            if (!has_int) {
                A_b = D_arr[b].M :* sqrt(wt_b)
                if (stack_direct) A_all[|r_off, 1 \ (r_off + n_b - 1), d_s|] = A_b
                else {
                    hqrd(A_b, Q_x = ., tau = ., R_ch = .)
                    R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
                }
            }
            else if (d_s > 0) {
                mean_d_row = quadcross(wt_b, D_arr[b].M) / ws_b
                for (j = 1; j <= d_s; j = j + 1) {
                    if (d_const_mat[b, j]) mean_d_row[j] = d_min_mat[b, j]
                }
                z_bar[b, ] = mean_d_row
                D_c = D_arr[b].M :- mean_d_row
                for (j = 1; j <= d_s; j = j + 1) {
                    if (d_const_mat[b, j]) D_c[, j] = J(n_b, 1, 0)
                }
                A_b = D_c :* sqrt(wt_b)
                if (stack_direct) A_all[|r_off, 1 \ (r_off + n_b - 1), d_s|] = A_b
                else {
                    hqrd(A_b, Q_x = ., tau = ., R_ch = .)
                    R_stack = _cdid_sel_qr_compress(R_stack \ R_ch)
                }
            }
            r_off = r_off + n_b
        }
        if (!ok_wls)
            _cdid_sel_fail("rank", "Selection block final information is singular or nonfinite.", 459)
        if (d_s > 0) {
            if (stack_direct) hqrd(A_all, Q_x = ., tau = ., R_x = .)
            else R_x = _cdid_sel_qr_compress(R_stack)
        }
        out.cell_cache_ready = 1
    }
    else {
        dummy_bh = 0
        ok_wls = _cdid_sel_block_wls(rdr, map, cell_ids, alpha, beta, eps_prob, 1,
            dummy_bh, m_b = ., z_bar = ., alpha_new = ., R_stack = .)
        if (!ok_wls)
            _cdid_sel_fail("rank", "Selection block final information is singular or nonfinite.", 459)
        if (d_s > 0) {
            R_all = _cdid_sel_qr_compress(R_stack)
            R_x = (rows(R_all) > 0 ? R_all[, 1..d_s] : J(0, d_s, 0))
        }
        out.cell_cache_ready = 0
    }

    if (d_s > 0) {
        if (rows(R_x) < d_s)
            _cdid_sel_fail("rank", "Selection block final information is rank deficient.", 459)
        hqrdp(R_x, Q_x = ., tau = ., R_step = ., piv_step = .)
        if (rows(R_step) != d_s | missing(R_step))
            _cdid_sel_fail("rank", "Selection block final information is singular.", 459)
        d_abs = abs(diagonal(R_step))
        max_d = max(d_abs)
        cutoff = rank_tol * max_d
        if (missing(d_abs) | max_d <= 0 | any(d_abs :<= cutoff))
            _cdid_sel_fail("rank", "Selection block final information is singular or near-singular.", 459)
    }

    out.m_b = m_b
    out.z_bar = z_bar
    out.alpha = alpha
    out.beta = beta
    out.theta = (has_int ? (alpha \ beta) : beta)
    out.R = R_step
    out.H_chol = J(0, 0, 0)
    out.pivot = piv_step
    out.H_inv = (d_s > 0 ? rdr.N :* _cdid_nuisance_gram(R_step, piv_step, I(d_s)) : J(0, 0, 0))
    out.iterations = iter
    out.deviance = dev
    out.rank = p_m
    out.rank_tol = rank_tol
    out.cutoff = cutoff
    out.rcond_R = (d_s > 0 ? _cdid_nuisance_rcond(R_step) : 1)
    if (d_s > 0 & out.rcond_R <= 0)
        _cdid_sel_fail("rank", "Selection block final R factor has zero reciprocal condition.", 459)
    out.solver = "block_fwl_qr"
    return(out)
}

/* Fit all selection blocks in canonical block order and assemble compact
   selection model state. Global theta follows R's canonical stacking:
   all cell intercepts alpha_all (1..C if has_intercept) followed by all block
   slope vectors beta_all (blocks 1..M). */
struct cdid_sel_fit_result scalar cdid_sel_fit(struct cdid_sel_reader scalar rdr,
    struct cdid_sel_cell_map scalar map, real scalar eps, real scalar maxit)
{
    struct cdid_sel_fit_result scalar res
    struct _cdid_sel_ivec colvector seen_y
    real colvector cell_ids, u_idx, w_fit, y_fit, seen_m, seen_ysum
    real rowvector x_sum
    real matrix X_raw
    real scalar m, a_pos, b_pos, b_m, d, d_s, has_int, is_sep_pair
    real scalar start_pos, c, k, cur_k, n_fit, ws_k, j, h, n_seen_k, ysum, m0

    _cdid_sel_ensure_pop_obs(rdr)
    d = rdr.spec.d
    d_s = rdr.spec.d_s
    has_int = rdr.spec.has_intercept

    res.selection_model = rdr.selection_model
    res.selection_pooling = rdr.selection_pooling
    res.selection_universe = rdr.selection_universe
    res.N = rdr.N
    res.n_U = rdr.n_U
    res.n_A = rdr.n_A
    res.n_cells = map.n_cells
    res.n_blocks = map.n_blocks
    res.n_pair_cohorts = map.n_pair_cohorts
    res.max_chunk_rows = rdr.max_chunk_rows
    res.eps = eps
    res.maxit = maxit
    res.spec = rdr.spec
    res.map = map
    res.cohorts = rdr.cohorts
    res.all_times = rdr.all_times
    res.mass_F = rdr.mass_F
    res.mass_U = rdr.mass_U
    res.repeat_mass = rdr.repeat_mass
    res.a_U = rdr.a_U
    res.fit_w_sum = sum(rdr.w)

    res.blocks = cdid_sel_block_fit(map.n_blocks)
    res.alpha_start = J(map.n_blocks, 1, 0)
    res.beta_start = J(map.n_blocks, 1, 0)
    res.block_p = J(map.n_blocks, 1, 0)
    res.n_alpha = (has_int ? map.n_cells : 0)
    res.n_beta = map.n_blocks * d_s
    res.total_params = res.n_alpha + res.n_beta
    res.alpha = J(res.n_alpha, 1, 0)
    res.beta = J(res.n_beta, 1, 0)

    is_sep_pair = (rdr.selection_model == "pair_logit" & rdr.selection_pooling == "separate" & d > 0)
    if (is_sep_pair) {
        if (eps <= 0 | eps >= . | !_cdid_sel_is_pos_int(maxit, 1))
            _cdid_sel_fail("input", "Invalid GLM controls eps/maxit.", 498)
    }
    cur_k = 0
    n_seen_k = 0
    ws_k = 0

    a_pos = 1
    b_pos = 1
    for (m = 1; m <= map.n_blocks; m = m + 1) {
        if (is_sep_pair & map.blocks[m, 3] == 1) {
            start_pos = map.blocks[m, 2]
            cell_ids = map.block_cell_ids[|start_pos \ start_pos|]
            c = cell_ids[1]
            k = map.cells[c, 3]
            if (k != cur_k) {
                if (rdr.mass_U[k] <= 0 | rdr.a_U[k] <= 0)
                    _cdid_sel_fail("support", "Selection block cohort has zero repeat-observer mass or U mass.", 459)
                u_idx = _cdid_sel_coh_u(rdr, k)
                n_fit = rows(u_idx)
                if (n_fit < d)
                    _cdid_sel_fail("rank", "Selection block has fewer active rows than parameters.", 459)
                X_raw = rdr.Z[u_idx, ]
                w_fit = rdr.w[u_idx]
                if (missing(X_raw) | missing(w_fit) | any(w_fit :<= 0))
                    _cdid_sel_fail("input", "Invalid X, y, or w in separate pair selection block.", 498)
                x_sum = colsum(X_raw)
                for (j = 1; j <= d; j = j + 1) {
                    if (x_sum[j] == 0) {
                        if (all(X_raw[, j] :== 0))
                            _cdid_sel_fail("rank", "Logit zero predictor; model unchanged.", 459)
                    }
                    for (h = 1; h < j; h = h + 1) {
                        if (x_sum[j] == x_sum[h]) {
                            if (all(X_raw[, j] :== X_raw[, h]))
                                _cdid_sel_fail("rank", "A covariate in the selection model duplicates another column or the intercept (for example, a constant covariate); remove it.", 459)
                        }
                    }
                }
                ws_k = sum(w_fit)
                cur_k = k
                n_seen_k = 0
                seen_m = J(map.n_blocks, 1, 0)
                seen_ysum = J(map.n_blocks, 1, 0)
                seen_y = _cdid_sel_ivec(map.n_blocks)
            }
            y_fit = _cdid_sel_cell_y(rdr, map, c, u_idx)
            ysum = sum(y_fit)
            m0 = 0
            for (j = 1; j <= n_seen_k; j = j + 1) {
                if (seen_ysum[j] == ysum) {
                    if (all(y_fit :== seen_y[j].v)) {
                        m0 = seen_m[j]
                        break
                    }
                }
            }
            if (m0 > 0) {
                res.blocks[m] = res.blocks[m0]
                res.blocks[m].block_id = m
                res.blocks[m].cell_ids = cell_ids
                res.blocks[m].dup_of_block = m0
            }
            else {
                res.blocks[m] = _cdid_sel_fit_sep_pair_core(rdr, m, cell_ids, k, u_idx, X_raw, w_fit, ws_k, y_fit, eps, maxit)
                n_seen_k = n_seen_k + 1
                seen_m[n_seen_k] = m
                seen_ysum[n_seen_k] = ysum
                seen_y[n_seen_k].v = y_fit
            }
        }
        else {
            res.blocks[m] = cdid_sel_fit_block(rdr, map, m, eps, maxit)
        }
        b_m = res.blocks[m].B_m
        res.block_p[m] = res.blocks[m].p_m
        if (has_int) {
            res.alpha_start[m] = a_pos
            res.alpha[res.blocks[m].cell_ids] = res.blocks[m].alpha
            a_pos = a_pos + b_m
        }
        if (d_s > 0) {
            res.beta_start[m] = res.n_alpha + b_pos
            res.beta[|b_pos \ (b_pos + d_s - 1)|] = res.blocks[m].beta
            b_pos = b_pos + d_s
        }
    }
    res.theta = res.alpha \ res.beta
    res.cell_block = map.cells[, 7]
    res.cell_local = map.cells[, 8]
    return(res)
}

/* Predict raw cell linear predictor eta and raw logistic probability p_raw on
   supplied design rows Z_sub (column 1 = eta, column 2 = p_raw). */
real matrix cdid_sel_pred_cell(struct cdid_sel_fit_result scalar fit,
    real scalar cell_id, real matrix Z_sub)
{
    struct cdid_sel_block_fit scalar blk
    real colvector eta, p_raw
    real scalar n_sub, b_id, loc_b

    if (!_cdid_sel_is_pos_int(cell_id, 1) | cell_id > fit.n_cells)
        _cdid_sel_fail("input", "Invalid cell_id in cdid_sel_pred_cell.", 498)
    n_sub = rows(Z_sub)
    if (cols(Z_sub) != fit.spec.d | (n_sub > 0 & missing(Z_sub)))
        _cdid_sel_fail("input", "Z_sub must be finite with d columns.", 498)
    if (n_sub == 0) return(J(0, 2, 0))

    b_id = fit.cell_block[cell_id]
    loc_b = fit.cell_local[cell_id]
    blk = fit.blocks[b_id]

    if (blk.d == 0) eta = J(n_sub, 1, 0)
    else if (!blk.has_intercept) eta = Z_sub[, fit.spec.slope_cols] * blk.beta
    else if (blk.d_s == 0) eta = J(n_sub, 1, blk.alpha[loc_b])
    else eta = blk.alpha[loc_b] :+ (Z_sub[, fit.spec.slope_cols] * blk.beta)

    p_raw = _cdid_sel_logistic(eta)
    return((eta, p_raw))
}

/* Predict pair-cohort probabilities, floor activity, and selection weights on
   requested U unit_ids. For independent_period_logit, p_pair = p_s * p_t and
   the 1e-6 floor is applied only to the final pair product p_pair. */
struct cdid_sel_pred scalar cdid_sel_pred_pair(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, real scalar pc, real colvector unit_ids)
{
    struct cdid_sel_pred scalar pr
    real matrix Z_sub, pred_s, pred_t
    real scalar n_req, k, cell_s, cell_t

    if (!_cdid_sel_is_pos_int(pc, 1) | pc > fit.n_pair_cohorts)
        _cdid_sel_fail("input", "Invalid pair-cohort index pc.", 498)
    n_req = rows(unit_ids)
    if (cols(unit_ids) != 1 | (n_req > 0 & (missing(unit_ids) |
        any(unit_ids :< 1) | any(unit_ids :> fit.N) | any(unit_ids :!= floor(unit_ids)))))
        _cdid_sel_fail("input", "unit_ids must be valid 1..N integer indices.", 498)
    if (n_req > 0) {
        if (any(rdr.U[unit_ids] :!= 1))
            _cdid_sel_fail("input", "Cannot predict selection probabilities outside U.", 498)
    }

    k = fit.map.pair_cohorts[pc, 1]
    cell_s = fit.map.pair_cohorts[pc, 6]
    cell_t = fit.map.pair_cohorts[pc, 7]

    pr.cohort = k
    pr.a_k = fit.a_U[k]
    pr.mass_U_k = fit.mass_U[k]
    if (pr.a_k <= 0 | pr.mass_U_k <= 0)
        _cdid_sel_fail("support", "Pair-cohort has nonpositive repeat rate or U mass.", 459)

    pr.w = (n_req > 0 ? rdr.w[unit_ids] : J(0, 1, 0))
    pr.rho = (n_req > 0 ? rdr.rho[unit_ids] : J(0, 1, 0))

    if (fit.selection_model == "pair_logit") {
        Z_sub = (n_req > 0 ? rdr.Z[unit_ids, ] : J(0, fit.spec.d, 0))
        pred_s = cdid_sel_pred_cell(fit, cell_s, Z_sub)
        pr.cell_id = cell_s
        pr.cell_s = cell_s
        pr.cell_t = cell_s
        pr.eta = pred_s[, 1]
        pr.eta_s = pred_s[, 1]
        pr.eta_t = J(n_req, 1, 0)
        pr.p_raw = pred_s[, 2]
        pr.p_s = pred_s[, 2]
        pr.p_t = J(n_req, 1, 1)
        pr.p_pair = pr.p_raw
        pr.p_weight = cdid_sel_clip_floor(pr.p_pair)
        pr.active = cdid_sel_floor_act(pr.p_pair)
        pr.sel_weight = (pr.w :* pr.a_k) :/ (pr.rho :* pr.p_weight)
        pr.score_factor = pr.active :* (1 :- pr.p_pair) :* (pr.p_pair :/ pr.p_weight)
        return(pr)
    }

    Z_sub = (n_req > 0 ? rdr.Z[unit_ids, ] : J(0, fit.spec.d, 0))
    pred_s = cdid_sel_pred_cell(fit, cell_s, Z_sub)
    pred_t = cdid_sel_pred_cell(fit, cell_t, Z_sub)
    pr.cell_id = 0
    pr.cell_s = cell_s
    pr.cell_t = cell_t
    pr.eta = J(0, 1, 0)
    pr.eta_s = pred_s[, 1]
    pr.eta_t = pred_t[, 1]
    pr.p_raw = J(0, 1, 0)
    pr.p_s = pred_s[, 2]
    pr.p_t = pred_t[, 2]
    pr.p_pair = pr.p_s :* pr.p_t
    pr.p_weight = cdid_sel_clip_floor(pr.p_pair)
    pr.active = cdid_sel_floor_act(pr.p_pair)
    pr.sel_weight = (pr.w :* pr.a_k) :/ (pr.rho :* pr.p_weight)
    pr.score_factor = J(0, 1, 0)
    return(pr)
}

/* Build probe targets and evaluate pair probabilities and multipliers.
   When rows(probe_cohorts) == 0, builds default cohort U-weighted mean targets
   (one per pair-cohort). When rows(probe_cohorts) > 0, matches each probe row
   against all pair-cohorts of its cohort in canonical order, mirroring R's
   build_selection_probe_targets and evaluate_selection_table_if. */
struct cdid_sel_probes scalar cdid_sel_pred_probes(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, real colvector probe_cohorts,
    real matrix Z_probe, real colvector rho_probe)
{
    struct cdid_sel_probes scalar pb
    struct cdid_sel_block_fit scalar blk_s, blk_t
    real matrix Z_p_slope, target_z
    real colvector target_pc, target_k, target_rho, rk, match_q, q_clip, a_vec
    real rowvector z_eval
    real scalar n_p, n_pc, d_s, q, k, p_idx, j, n_t, c_s, c_t, b_s, l_s, b_t, l_t
    real scalar eta_1, eta_2, ps, pt, pp, rho_eval, all_labels

    n_p = rows(probe_cohorts)
    n_pc = fit.n_pair_cohorts
    d_s = fit.spec.d_s

    target_pc = J(0, 1, 0)
    target_k = J(0, 1, 0)
    target_rho = J(0, 1, 0)
    target_z = J(0, d_s, 0)

    if (n_p == 0) {
        for (q = 1; q <= n_pc; q = q + 1) {
            k = fit.map.pair_cohorts[q, 1]
            rk = _cdid_sel_coh_u(rdr, k)
            z_eval = (d_s > 0 ? quadcross(rdr.w[rk], rdr.Z[rk, fit.spec.slope_cols]) / sum(rdr.w[rk]) : J(1, 0, 0))
            rho_eval = sum(rdr.w[rk] :* rdr.rho[rk]) / sum(rdr.w[rk])
            target_pc = target_pc \ q
            target_k = target_k \ k
            target_rho = target_rho \ rho_eval
            if (d_s > 0) target_z = target_z \ z_eval
            else target_z = J(rows(target_pc), 0, 0)
        }
    }
    else {
        if (rows(Z_probe) != n_p | missing(Z_probe))
            _cdid_sel_fail("input", "Z_probe must be finite with rows(probe_cohorts) rows.", 498)
        if (cols(Z_probe) == d_s) Z_p_slope = Z_probe
        else if (cols(Z_probe) == fit.spec.d) Z_p_slope = (d_s > 0 ? Z_probe[, fit.spec.slope_cols] : J(n_p, 0, 0))
        else _cdid_sel_fail("input", "Z_probe column count must match selection design.", 498)

        if (rows(rho_probe) == 0) rho_eval = 1
        else if (rows(rho_probe) != n_p | cols(rho_probe) != 1 | missing(rho_probe) | any(rho_probe :<= 0))
            _cdid_sel_fail("input", "rho_probe must be finite and positive with length n_probes.", 498)

        all_labels=1
        for(p_idx=1;p_idx<=n_p;p_idx=p_idx+1)
            if(!any(fit.cohorts:==probe_cohorts[p_idx]))all_labels=0
        for (p_idx = 1; p_idx <= n_p; p_idx = p_idx + 1) {
            if(all_labels) k = _cdid_sel_cohort_pos(probe_cohorts[p_idx], fit.cohorts)
            else {
                k=probe_cohorts[p_idx]
                if(!_cdid_sel_is_pos_int(k,1) | k>rows(fit.cohorts))
                    _cdid_sel_fail("input", "Probe cohort positions must be in 1..n_cohorts.", 498)
            }
            rho_eval = (rows(rho_probe) == 0 ? 1 : rho_probe[p_idx])
            match_q = _cdid_sel_index(fit.map.pair_cohorts[, 1] :== k)
            for (j = 1; j <= rows(match_q); j = j + 1) {
                q = match_q[j]
                target_pc = target_pc \ q
                target_k = target_k \ k
                target_rho = target_rho \ rho_eval
                if (d_s > 0) target_z = target_z \ Z_p_slope[p_idx, ]
                else target_z = J(rows(target_pc), 0, 0)
            }
        }
        if (rows(target_pc) == 0)
            _cdid_sel_fail("input", "Supplied probes did not match any requested pair-cohorts.", 498)
    }

    n_t = rows(target_pc)
    pb.n_targets = n_t
    pb.n_pair_cohorts = n_pc
    pb.d_s = d_s
    pb.target_pc = target_pc
    pb.target_cohort = target_k
    pb.target_rho = target_rho
    pb.target_z_slope = target_z
    pb.pair_probabilities = J(n_t, 1, 0)
    pb.multipliers = J(n_t, 1, 0)
    pb.active_mask = J(n_t, 1, 0)
    pb.eta_s = J(n_t, 1, 0)
    pb.eta_t = J(n_t, 1, 0)
    pb.p_s = J(n_t, 1, 0)
    pb.p_t = J(n_t, 1, 0)

    for (j = 1; j <= n_t; j = j + 1) {
        q = target_pc[j]
        k = target_k[j]
        c_s = fit.map.pair_cohorts[q, 6]
        c_t = fit.map.pair_cohorts[q, 7]
        b_s = fit.cell_block[c_s]
        l_s = fit.cell_local[c_s]
        blk_s = fit.blocks[b_s]
        if (blk_s.d == 0) {
            eta_1 = 0
            ps = 0.5
        }
        else {
            eta_1 = (blk_s.has_intercept ? blk_s.alpha[l_s] : 0) +
                    (d_s > 0 ? sum(target_z[j, ] :* blk_s.beta') : 0)
            ps = _cdid_sel_logistic(J(1, 1, eta_1))[1, 1]
        }
        if (fit.selection_model == "pair_logit") {
            eta_2 = 0
            pt = 1
            pp = ps
        }
        else {
            b_t = fit.cell_block[c_t]
            l_t = fit.cell_local[c_t]
            blk_t = fit.blocks[b_t]
            if (blk_t.d == 0) {
                eta_2 = 0
                pt = 0.5
            }
            else {
                eta_2 = (blk_t.has_intercept ? blk_t.alpha[l_t] : 0) +
                        (d_s > 0 ? sum(target_z[j, ] :* blk_t.beta') : 0)
                pt = _cdid_sel_logistic(J(1, 1, eta_2))[1, 1]
            }
            pp = ps * pt
        }
        if (pp <= 0 | pp >= 1 | pp >= .)
            _cdid_sel_fail("boundary", "Probe pair selection probability reached boundary.", 430)
        pb.eta_s[j] = eta_1
        pb.eta_t[j] = eta_2
        pb.p_s[j] = ps
        pb.p_t[j] = pt
        pb.pair_probabilities[j] = pp
    }

    q_clip = cdid_sel_clip_floor(pb.pair_probabilities)
    pb.active_mask = cdid_sel_floor_act(pb.pair_probabilities)
    a_vec = fit.a_U[target_k]
    pb.multipliers = a_vec :/ (target_rho :* q_clip)
    return(pb)
}

/* Evaluate original-unit influence functions on a chunk unit_ids of 1..N for
   a single fitted block m. Sums cell scores by original unit i across all cells
   in block m before solving the block Fisher information system. */
void _cdid_sel_block_if_chunk(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, real scalar block_id,
    real colvector unit_ids, real matrix xi_alpha_out, real matrix xi_beta_out)
{
    struct cdid_sel_block_fit scalar blk
    real matrix Z_act, S_act, z_sol, xi_raw, xi_theta, s_mat, shift_mat
    real matrix D_act, score_beta, pred_b
    real colvector act_loc, act_u, y_act, w_act, s_b, s_sum, act_any, act_any_loc, act_any_u
    real rowvector ord
    real scalar n_req, c, k, b, rk, is_contig, lo, hi, p1, p2, use_cache

    n_req = rows(unit_ids)
    blk = fit.blocks[block_id]
    xi_alpha_out = (blk.has_intercept ? J(n_req, blk.B_m, 0) : J(n_req, 0, 0))
    xi_beta_out = J(n_req, blk.d_s, 0)
    if (blk.p_m == 0 | n_req == 0) return

    is_contig = 0
    lo = unit_ids[1]
    hi = unit_ids[n_req]
    if (hi - lo == n_req - 1) {
        if (unit_ids == (lo::hi)) is_contig = 1
    }
    use_cache = (is_contig & blk.cell_cache_ready == 1 & fit.N == rdr.N & fit.n_U == rdr.n_U)

    if (blk.solver == "cp04_separate_pair") {
        c = blk.cell_ids[1]
        k = fit.map.cells[c, 3]
        if (use_cache) {
            _cdid_sel_range_bounds(blk.cell_u_idx[1].v, lo, hi, p1 = ., p2 = .)
            if (p1 > p2) return
            act_u = blk.cell_u_idx[1].v[|p1 \ p2|]
            act_loc = act_u :- (lo - 1)
            s_b = blk.cell_s_u[1].v[|p1 \ p2|]
            S_act = rdr.Z[act_u, ] :* s_b
        }
        else {
            act_loc = _cdid_sel_index((rdr.cohort_pos[unit_ids] :== k) :& (rdr.U[unit_ids] :== 1))
            if (rows(act_loc) == 0) return
            act_u = unit_ids[act_loc]
            Z_act = rdr.Z[act_u, ]
            w_act = rdr.w[act_u]
            y_act = _cdid_sel_cell_y(rdr, fit.map, c, act_u)
            pred_b = cdid_sel_pred_cell(fit, c, Z_act)
            S_act = Z_act :* (w_act :* (y_act :- pred_b[, 2]))
        }
        if (rows(blk.H_inv) == blk.p_m) {
            xi_theta = (S_act * blk.H_inv')[, fit.spec.col_order]
            if (missing(xi_theta))
                _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate pair IF.", 504)
        }
        else {
            z_sol = solvelower(blk.H_chol', S_act', rk = ., 0)
            if (rk != blk.p_m | missing(z_sol))
                _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate pair IF.", 504)
            xi_raw = solveupper(blk.H_chol, z_sol, rk = ., 0)
            if (rk != blk.p_m | missing(xi_raw))
                _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate pair IF.", 504)
            ord = fit.spec.col_order
            xi_theta = (xi_raw')[, ord]
        }
        if (blk.has_intercept) {
            xi_alpha_out[act_loc, ] = xi_theta[, 1]
            if (blk.d_s > 0) xi_beta_out[act_loc, ] = xi_theta[, 2..blk.p_m]
        }
        else {
            xi_beta_out[act_loc, ] = xi_theta
        }
        return
    }

    s_mat = (blk.has_intercept ? J(n_req, blk.B_m, 0) : J(0, 0, 0))
    s_sum = J(n_req, 1, 0)
    shift_mat = (blk.has_intercept & blk.d_s > 0 ? J(n_req, blk.d_s, 0) : J(0, 0, 0))
    act_any = J(n_req, 1, 0)

    for (b = 1; b <= blk.B_m; b = b + 1) {
        c = blk.cell_ids[b]
        k = fit.map.cells[c, 3]
        if (use_cache) {
            _cdid_sel_range_bounds(blk.cell_u_idx[b].v, lo, hi, p1 = ., p2 = .)
            if (p1 > p2) continue
            act_u = blk.cell_u_idx[b].v[|p1 \ p2|]
            act_loc = act_u :- (lo - 1)
            s_b = blk.cell_s_u[b].v[|p1 \ p2|]
        }
        else {
            act_loc = _cdid_sel_index((rdr.cohort_pos[unit_ids] :== k) :& (rdr.U[unit_ids] :== 1))
            if (rows(act_loc) == 0) continue
            act_u = unit_ids[act_loc]
            Z_act = rdr.Z[act_u, ]
            w_act = rdr.w[act_u]
            y_act = _cdid_sel_cell_y(rdr, fit.map, c, act_u)
            pred_b = cdid_sel_pred_cell(fit, c, Z_act)
            s_b = w_act :* (y_act :- pred_b[, 2])
        }
        act_any[act_loc] = J(rows(act_loc), 1, 1)
        if (blk.has_intercept) {
            s_mat[act_loc, b] = s_b
            s_sum[act_loc] = s_sum[act_loc] + s_b
            if (blk.d_s > 0)
                shift_mat[act_loc, ] = shift_mat[act_loc, ] + s_b * blk.z_bar[b, ]
        }
        else {
            s_sum[act_loc] = s_sum[act_loc] + s_b
        }
    }

    act_any_loc = _cdid_sel_index(act_any :== 1)
    if (rows(act_any_loc) == 0) return
    act_any_u = unit_ids[act_any_loc]

    if (blk.d_s > 0) {
        D_act = rdr.Z[act_any_u, fit.spec.slope_cols]
        score_beta = (blk.has_intercept ? D_act :* s_sum[act_any_loc] - shift_mat[act_any_loc, ] :
                                          D_act :* s_sum[act_any_loc])
        if (rows(blk.H_inv) == blk.d_s) {
            xi_beta_out[act_any_loc, ] = score_beta * blk.H_inv'
        }
        else {
            xi_beta_out[act_any_loc, ] = (fit.N :* _cdid_nuisance_gram(blk.R, blk.pivot, score_beta'))'
        }
    }

    if (blk.has_intercept) {
        xi_alpha_out[act_any_loc, ] = s_mat[act_any_loc, ] :/ (blk.m_b')
        if (blk.d_s > 0)
            xi_alpha_out[act_any_loc, ] = xi_alpha_out[act_any_loc, ] -
                xi_beta_out[act_any_loc, ] * blk.z_bar'
    }
}

/* Evaluate original-unit influence functions for theta = (alpha_all, beta_all)
   and (optionally) a_k^U on an arbitrary chunk unit_ids of 1..N. */
struct cdid_sel_if_chunk scalar cdid_sel_eval_if_chunk(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, real colvector unit_ids,
    real scalar include_ak)
{
    struct cdid_sel_if_chunk scalar ch
    real matrix xa_m, xb_m
    real scalar n_req, m, d_s, has_int, b_pos, r, u, k

    n_req = rows(unit_ids)
    if (cols(unit_ids) != 1 | (n_req > 0 & (missing(unit_ids) |
        any(unit_ids :< 1) | any(unit_ids :> fit.N) | any(unit_ids :!= floor(unit_ids)))))
        _cdid_sel_fail("input", "unit_ids must be valid 1..N integer indices.", 498)

    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept
    ch.unit_ids = unit_ids
    ch.xi_alpha = J(n_req, fit.n_alpha, 0)
    ch.xi_beta = J(n_req, fit.n_beta, 0)

    b_pos = 1
    for (m = 1; m <= fit.n_blocks; m = m + 1) {
        _cdid_sel_block_if_chunk(fit, rdr, m, unit_ids, xa_m = ., xb_m = .)
        if (has_int & cols(xa_m) > 0)
            ch.xi_alpha[, fit.blocks[m].cell_ids'] = xa_m
        if (d_s > 0) {
            ch.xi_beta[, b_pos..(b_pos + d_s - 1)] = xb_m
            b_pos = b_pos + d_s
        }
    }
    ch.theta_if = (ch.xi_alpha, ch.xi_beta)

    if (!include_ak) {
        ch.cohort_a_if = J(0, 0, 0)
        return(ch)
    }

    ch.cohort_a_if = J(n_req, rows(fit.cohorts), 0)
    for (r = 1; r <= n_req; r = r + 1) {
        u = unit_ids[r]
        if (rdr.U[u] != 1) continue
        k = rdr.cohort_pos[u]
        if (k >= 1 & k <= rows(fit.cohorts)) {
            if (fit.mass_U[k] > 0)
                ch.cohort_a_if[r, k] = rdr.w[u] * (rdr.A[u] - fit.a_U[k]) / fit.mass_U[k]
        }
    }
    return(ch)
}

/* Build upstream channel map for pair-cohorts and cohort repeat-rate channels.
   pc_info columns (12):
   (cohort_pos, cell_id, cell_s, cell_t, block_id, local_b, block_s, local_s,
    block_t, local_t, p_pc, same_block). */
struct cdid_sel_up_map scalar cdid_sel_up_map_build(struct cdid_sel_fit_result scalar fit)
{
    struct cdid_sel_up_map scalar um
    real scalar n_pc, pc, k, c_s, c_t, b_s, l_s, b_t, l_t, d_s, has_int, p_pc, same_b

    n_pc = fit.n_pair_cohorts
    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept

    um.selection_model = fit.selection_model
    um.selection_pooling = fit.selection_pooling
    um.n_pair_cohorts = n_pc
    um.n_cells = fit.n_cells
    um.n_blocks = fit.n_blocks
    um.n_cohorts = rows(fit.cohorts)
    um.total_params = fit.total_params
    um.total_channels = fit.total_params + rows(fit.cohorts)
    um.cohorts = fit.cohorts
    um.a_U = fit.a_U
    um.mass_U = fit.mass_U
    um.alpha_start = fit.alpha_start
    um.beta_start = fit.beta_start
    um.block_p = fit.block_p
    um.pc_info = J(n_pc, 12, 0)

    for (pc = 1; pc <= n_pc; pc = pc + 1) {
        k = fit.map.pair_cohorts[pc, 1]
        c_s = fit.map.pair_cohorts[pc, 6]
        c_t = fit.map.pair_cohorts[pc, 7]
        b_s = fit.cell_block[c_s]
        l_s = fit.cell_local[c_s]
        b_t = fit.cell_block[c_t]
        l_t = fit.cell_local[c_t]
        if (fit.selection_model == "pair_logit") {
            p_pc = (has_int ? 1 + d_s : d_s)
            um.pc_info[pc, ] = (k, c_s, c_s, c_s, b_s, l_s, b_s, l_s, b_s, l_s, p_pc, 1)
        }
        else {
            same_b = (b_s == b_t)
            if (has_int) p_pc = (same_b ? 2 + d_s : 2 + 2 * d_s)
            else p_pc = (same_b ? d_s : 2 * d_s)
            um.pc_info[pc, ] = (k, 0, c_s, c_t, (same_b ? b_s : 0), 0, b_s, l_s, b_t, l_t, p_pc, same_b)
        }
    }
    return(um)
}

/* Compute rowwise derivative of log(w_sel) wrt local pair-cohort selection
   parameters (d_theta_local, gated by strict floor mask p_pair > 1e-6) and wrt
   a_k^U (d_ak = 1 / a_k^U). */
struct cdid_sel_log_deriv scalar cdid_sel_log_deriv_eval(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, struct cdid_sel_up_map scalar um,
    real scalar pc, real colvector unit_ids)
{
    struct cdid_sel_log_deriv scalar ld
    struct cdid_sel_pred scalar pr
    real matrix D_sub
    real colvector g_eta, g_s, g_t
    real scalar n_req, d_s, has_int, same_b

    pr = cdid_sel_pred_pair(fit, rdr, pc, unit_ids)
    n_req = rows(unit_ids)
    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept
    D_sub = (n_req > 0 & d_s > 0 ? rdr.Z[unit_ids, fit.spec.slope_cols] : J(n_req, 0, 0))

    ld.pc = pc
    ld.cohort = pr.cohort
    ld.cohort_pos = um.pc_info[pc, 1]
    ld.a_k = pr.a_k
    ld.mass_U_k = pr.mass_U_k
    ld.block_id = um.pc_info[pc, 5]
    ld.local_b = um.pc_info[pc, 6]
    ld.block_s = um.pc_info[pc, 7]
    ld.local_s = um.pc_info[pc, 8]
    ld.block_t = um.pc_info[pc, 9]
    ld.local_t = um.pc_info[pc, 10]
    ld.p_pc = um.pc_info[pc, 11]
    same_b = um.pc_info[pc, 12]
    ld.active = pr.active
    ld.p_pair = pr.p_pair
    ld.p_weight = pr.p_weight
    ld.p_s = pr.p_s
    ld.p_t = pr.p_t
    ld.d_ak = J(n_req, 1, 1 / pr.a_k)

    if (ld.p_pc == 0 | n_req == 0) {
        ld.d_theta_local = J(n_req, ld.p_pc, 0)
        return(ld)
    }

    if (fit.selection_model == "pair_logit") {
        g_eta = -pr.active :* (1 :- pr.p_pair)
        if (!has_int) ld.d_theta_local = D_sub :* g_eta
        else if (d_s == 0) ld.d_theta_local = g_eta
        else ld.d_theta_local = (g_eta, D_sub :* g_eta)
        return(ld)
    }

    g_s = -pr.active :* (1 :- pr.p_s)
    g_t = -pr.active :* (1 :- pr.p_t)
    if (!has_int) {
        if (same_b) ld.d_theta_local = D_sub :* (g_s + g_t)
        else ld.d_theta_local = (D_sub :* g_s, D_sub :* g_t)
        return(ld)
    }
    if (d_s == 0) {
        ld.d_theta_local = (g_s, g_t)
        return(ld)
    }
    if (same_b) ld.d_theta_local = (g_s, g_t, D_sub :* (g_s + g_t))
    else ld.d_theta_local = (g_s, g_t, D_sub :* g_s, D_sub :* g_t)
    return(ld)
}

/* Map local pair-cohort parameter positions to global theta = (alpha_all, beta_all)
   positions 1..total_params. */
real colvector _cdid_sel_pc_global_idx(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_up_map scalar um, real scalar pc)
{
    real scalar c_s, c_t, b_s, b_t, same_b, d_s, has_int, bst_s, bst_t
    real colvector idx

    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept
    if (fit.spec.d == 0) return(J(0, 1, 0))

    c_s = um.pc_info[pc, 3]
    c_t = um.pc_info[pc, 4]
    b_s = um.pc_info[pc, 7]
    b_t = um.pc_info[pc, 9]
    same_b = um.pc_info[pc, 12]
    bst_s = fit.beta_start[b_s]

    if (fit.selection_model == "pair_logit") {
        idx = J(0, 1, 0)
        if (has_int) idx = idx \ c_s
        if (d_s > 0) idx = idx \ (bst_s::(bst_s + d_s - 1))
        return(idx)
    }

    bst_t = fit.beta_start[b_t]
    idx = J(0, 1, 0)
    if (has_int) idx = idx \ c_s \ c_t
    if (d_s > 0) {
        idx = idx \ (bst_s::(bst_s + d_s - 1))
        if (!same_b) idx = idx \ (bst_t::(bst_t + d_s - 1))
    }
    return(idx)
}

/* Initialize consumer moment accumulator for n_moments columns. */
struct cdid_sel_moments scalar cdid_sel_moments_init(struct cdid_sel_fit_result scalar fit,
    real scalar n_moments)
{
    struct cdid_sel_moments scalar acc

    if (!_cdid_sel_is_pos_int(n_moments, 0))
        _cdid_sel_fail("input", "n_moments must be a nonnegative integer.", 498)
    acc.n_moments = n_moments
    acc.total_params = fit.total_params
    acc.n_cohorts = rows(fit.cohorts)
    acc.G_theta = J(fit.total_params, n_moments, 0)
    acc.G_ak = J(rows(fit.cohorts), n_moments, 0)
    return(acc)
}

// Accumulate consumer moment sensitivities wrt log(w_sel) on a chunk of units
// for pair-cohort pc. moment_deriv_w is n_chunk x n_moments representing
// observer_weights / N (matching R accumulate_selection_table_moments with
// sign = -1 already built into cdid_sel_log_deriv_eval).
void cdid_sel_add_moments(struct cdid_sel_moments scalar acc,
    struct cdid_sel_fit_result scalar fit, struct cdid_sel_reader scalar rdr,
    struct cdid_sel_up_map scalar um, real scalar pc,
    real colvector unit_ids, real matrix moment_deriv_w)
{
    struct cdid_sel_log_deriv scalar ld
    real colvector g_idx
    real scalar k

    if (rows(moment_deriv_w) != rows(unit_ids) | cols(moment_deriv_w) != acc.n_moments |
        missing(moment_deriv_w))
        _cdid_sel_fail("input", "moment_deriv_w dimensions must match (rows(unit_ids), n_moments).", 498)
    if (rows(unit_ids) == 0 | acc.n_moments == 0) return

    ld = cdid_sel_log_deriv_eval(fit, rdr, um, pc, unit_ids)
    g_idx = _cdid_sel_pc_global_idx(fit, um, pc)
    if (rows(g_idx) > 0)
        acc.G_theta[g_idx, ] = acc.G_theta[g_idx, ] + quadcross(ld.d_theta_local, moment_deriv_w)
    k = ld.cohort_pos
    acc.G_ak[k, ] = acc.G_ak[k, ] + quadcross(ld.d_ak, moment_deriv_w)
}

/* Evaluate consumer moment adjustment IF on a chunk unit_ids using the
   profiled contraction identity (mirrors R's evaluate_pair_selection_loading /
   evaluate_selection_table_loading) without allocating an n_chunk x B_m matrix. */
real matrix cdid_sel_eval_loading(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, struct cdid_sel_moments scalar acc,
    real colvector unit_ids)
{
    struct cdid_sel_block_fit scalar blk, blk_c
    struct _cdid_sel_mvec colvector t_arr
    real matrix out_if, out_k, G_alpha, G_beta, t_vec, z_sol, V_m, U_m, W_b
    real matrix Z_act, D_act, S_act, pred_b, t_beta_eff
    real colvector act_loc, act_u, y_act, w_act, s_b, s_sum, act_any, act_any_loc, act_any_u
    real colvector k_sub, m_sub, valid_k, scale_ak, has_t, nz_cols
    real rowvector ord
    real scalar n_req, m, m_end, m_cur, m0, n_run, j_run, j0, n_act_blk
    real scalar d_s, has_int, bst, b, c, k, rk, is_contig, lo, hi, p1, p2, use_cache

    n_req = rows(unit_ids)
    if (n_req == 0 | acc.n_moments == 0) return(J(n_req, acc.n_moments, 0))
    out_if = J(n_req, acc.n_moments, 0)
    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept

    if (max(abs(acc.G_ak)) > 0) {
        act_loc = _cdid_sel_index(rdr.U[unit_ids] :== 1)
        if (rows(act_loc) > 0) {
            act_u = unit_ids[act_loc]
            k_sub = rdr.cohort_pos[act_u]
            valid_k = (k_sub :>= 1) :& (k_sub :<= rows(fit.cohorts))
            if (!all(valid_k)) {
                act_loc = select(act_loc, valid_k)
                act_u = select(act_u, valid_k)
                k_sub = select(k_sub, valid_k)
            }
            if (rows(act_loc) > 0) {
                m_sub = fit.mass_U[k_sub]
                if (any(m_sub :<= 0)) {
                    valid_k = (m_sub :> 0)
                    act_loc = select(act_loc, valid_k)
                    act_u = select(act_u, valid_k)
                    k_sub = select(k_sub, valid_k)
                    m_sub = select(m_sub, valid_k)
                }
                if (rows(act_loc) > 0) {
                    scale_ak = rdr.w[act_u] :* (rdr.A[act_u] :- fit.a_U[k_sub]) :/ m_sub
                    out_if[act_loc, ] = out_if[act_loc, ] + acc.G_ak[k_sub, ] :* scale_ak
                }
            }
        }
    }

    if (fit.total_params == 0) return(out_if)
    if (max(abs(acc.G_theta)) == 0) return(out_if)

    is_contig = 0
    lo = unit_ids[1]
    hi = unit_ids[n_req]
    if (hi - lo == n_req - 1) {
        if (unit_ids == (lo::hi)) is_contig = 1
    }
    use_cache = (is_contig & fit.N == rdr.N & fit.n_U == rdr.n_U)

    for (m = 1; m <= fit.n_blocks; m = m + 1) {
        blk = fit.blocks[m]
        if (blk.p_m == 0) continue

        if (use_cache & blk.solver == "cp04_separate_pair" & blk.cell_cache_ready == 1) {
            c = blk.cell_ids[1]
            k = fit.map.cells[c, 3]
            m_end = m
            while (m_end < fit.n_blocks) {
                if (fit.blocks[m_end + 1].solver != "cp04_separate_pair" |
                    fit.blocks[m_end + 1].cell_cache_ready != 1 |
                    fit.blocks[m_end + 1].p_m != blk.p_m) break
                if (fit.map.cells[fit.blocks[m_end + 1].cell_ids[1], 3] != k) break
                m_end = m_end + 1
            }
            _cdid_sel_range_bounds(blk.cell_u_idx[1].v, lo, hi, p1 = ., p2 = .)
            if (p1 <= p2) {
                n_run = m_end - m + 1
                t_arr = _cdid_sel_mvec(n_run)
                has_t = J(n_run, 1, 0)
                ord = fit.spec.col_order
                for (j_run = 1; j_run <= n_run; j_run = j_run + 1) {
                    m_cur = m + j_run - 1
                    blk_c = fit.blocks[m_cur]
                    G_alpha = (has_int ? acc.G_theta[blk_c.cell_ids, ] : J(0, acc.n_moments, 0))
                    bst = fit.beta_start[m_cur]
                    G_beta = (d_s > 0 ? acc.G_theta[|bst, 1 \ (bst + d_s - 1), acc.n_moments|] : J(0, acc.n_moments, 0))
                    if (max(abs(G_alpha)) == 0 & max(abs(G_beta)) == 0) continue
                    t_vec = J(blk_c.p_m, acc.n_moments, 0)
                    t_vec[ord', ] = (has_int ? (G_alpha \ G_beta) : G_beta)
                    m0 = blk_c.dup_of_block
                    if (m0 >= m & m0 < m_cur) {
                        j0 = m0 - m + 1
                        if (has_t[j0]) t_arr[j0].M = t_arr[j0].M + t_vec
                        else {
                            t_arr[j0].M = t_vec
                            has_t[j0] = 1
                        }
                    }
                    else {
                        if (has_t[j_run]) t_arr[j_run].M = t_arr[j_run].M + t_vec
                        else {
                            t_arr[j_run].M = t_vec
                            has_t[j_run] = 1
                        }
                    }
                }
                n_act_blk = sum(has_t)
                if (n_act_blk > 0) {
                    act_u = blk.cell_u_idx[1].v[|p1 \ p2|]
                    act_loc = act_u :- (lo - 1)
                    Z_act = rdr.Z[act_u, ]
                    out_k = J(rows(act_u), acc.n_moments, 0)
                    for (j_run = 1; j_run <= n_run; j_run = j_run + 1) {
                        if (!has_t[j_run]) continue
                        nz_cols = _cdid_sel_index((colsum(abs(t_arr[j_run].M)) :> 0)')
                        if (rows(nz_cols) == 0) continue
                        blk_c = fit.blocks[m + j_run - 1]
                        s_b = blk_c.cell_s_u[1].v[|p1 \ p2|]
                        S_act = Z_act :* s_b
                        if (rows(blk_c.H_inv) == blk_c.p_m) {
                            V_m = blk_c.H_inv * t_arr[j_run].M[, nz_cols']
                            if (missing(V_m))
                                _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate moment loading.", 504)
                        }
                        else {
                            z_sol = solvelower(blk_c.H_chol', t_arr[j_run].M[, nz_cols'], rk = ., 0)
                            V_m = solveupper(blk_c.H_chol, z_sol, rk = ., 0)
                            if (rk != blk_c.p_m | missing(V_m))
                                _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate moment loading.", 504)
                        }
                        out_k[, nz_cols'] = out_k[, nz_cols'] + S_act * V_m
                    }
                    out_if[act_loc, ] = out_if[act_loc, ] + out_k
                }
            }
            m = m_end
            continue
        }

        G_alpha = (has_int ? acc.G_theta[blk.cell_ids, ] : J(0, acc.n_moments, 0))
        bst = fit.beta_start[m]
        G_beta = (d_s > 0 ? acc.G_theta[|bst, 1 \ (bst + d_s - 1), acc.n_moments|] : J(0, acc.n_moments, 0))
        if (max(abs(G_alpha)) == 0 & max(abs(G_beta)) == 0) continue

        if (blk.solver == "cp04_separate_pair") {
            c = blk.cell_ids[1]
            k = fit.map.cells[c, 3]
            if (use_cache & blk.cell_cache_ready == 1) {
                _cdid_sel_range_bounds(blk.cell_u_idx[1].v, lo, hi, p1 = ., p2 = .)
                if (p1 > p2) continue
                act_u = blk.cell_u_idx[1].v[|p1 \ p2|]
                act_loc = act_u :- (lo - 1)
                s_b = blk.cell_s_u[1].v[|p1 \ p2|]
                S_act = rdr.Z[act_u, ] :* s_b
            }
            else {
                act_loc = _cdid_sel_index((rdr.cohort_pos[unit_ids] :== k) :& (rdr.U[unit_ids] :== 1))
                if (rows(act_loc) == 0) continue
                act_u = unit_ids[act_loc]
                Z_act = rdr.Z[act_u, ]
                w_act = rdr.w[act_u]
                y_act = _cdid_sel_cell_y(rdr, fit.map, c, act_u)
                pred_b = cdid_sel_pred_cell(fit, c, Z_act)
                S_act = Z_act :* (w_act :* (y_act :- pred_b[, 2]))
            }
            ord = fit.spec.col_order
            t_vec = J(blk.p_m, acc.n_moments, 0)
            t_vec[ord', ] = (has_int ? (G_alpha \ G_beta) : G_beta)
            if (rows(blk.H_inv) == blk.p_m) {
                V_m = blk.H_inv * t_vec
                if (missing(V_m))
                    _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate moment loading.", 504)
            }
            else {
                z_sol = solvelower(blk.H_chol', t_vec, rk = ., 0)
                V_m = solveupper(blk.H_chol, z_sol, rk = ., 0)
                if (rk != blk.p_m | missing(V_m))
                    _cdid_sel_fail("numerical", "Nonfinite triangular solve in separate moment loading.", 504)
            }
            out_if[act_loc, ] = out_if[act_loc, ] + S_act * V_m
            continue
        }

        U_m = J(d_s, acc.n_moments, 0)
        if (d_s > 0) {
            t_beta_eff = (has_int ? G_beta - quadcross(blk.z_bar, G_alpha) : G_beta)
            if (max(abs(t_beta_eff)) > 0) {
                if (rows(blk.H_inv) == d_s) U_m = blk.H_inv * t_beta_eff
                else U_m = fit.N :* _cdid_nuisance_gram(blk.R, blk.pivot, t_beta_eff)
            }
        }

        s_sum = J(n_req, 1, 0)
        act_any = J(n_req, 1, 0)
        for (b = 1; b <= blk.B_m; b = b + 1) {
            c = blk.cell_ids[b]
            k = fit.map.cells[c, 3]
            if (use_cache & blk.cell_cache_ready == 1) {
                _cdid_sel_range_bounds(blk.cell_u_idx[b].v, lo, hi, p1 = ., p2 = .)
                if (p1 > p2) continue
                act_u = blk.cell_u_idx[b].v[|p1 \ p2|]
                act_loc = act_u :- (lo - 1)
                s_b = blk.cell_s_u[b].v[|p1 \ p2|]
            }
            else {
                act_loc = _cdid_sel_index((rdr.cohort_pos[unit_ids] :== k) :& (rdr.U[unit_ids] :== 1))
                if (rows(act_loc) == 0) continue
                act_u = unit_ids[act_loc]
                Z_act = rdr.Z[act_u, ]
                w_act = rdr.w[act_u]
                y_act = _cdid_sel_cell_y(rdr, fit.map, c, act_u)
                pred_b = cdid_sel_pred_cell(fit, c, Z_act)
                s_b = w_act :* (y_act :- pred_b[, 2])
            }
            act_any[act_loc] = J(rows(act_loc), 1, 1)
            s_sum[act_loc] = s_sum[act_loc] + s_b
            if (has_int) {
                W_b = G_alpha[b, ] :/ blk.m_b[b]
                if (d_s > 0) W_b = W_b - blk.z_bar[b, ] * U_m
                out_if[act_loc, ] = out_if[act_loc, ] + s_b * W_b
            }
        }
        if (d_s > 0 & max(abs(U_m)) > 0) {
            act_any_loc = _cdid_sel_index(act_any :== 1)
            if (rows(act_any_loc) > 0) {
                act_any_u = unit_ids[act_any_loc]
                D_act = rdr.Z[act_any_u, fit.spec.slope_cols]
                out_if[act_any_loc, ] = out_if[act_any_loc, ] +
                    (D_act :* s_sum[act_any_loc]) * U_m
            }
        }
    }
    return(out_if)
}

/* Evaluate probe pair-probability IF and multiplier IF on a chunk of units. */
void _cdid_sel_probe_if_chunk(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_probes scalar pb, struct cdid_sel_if_chunk scalar ch,
    real matrix pair_prob_if_ch, real matrix mult_if_ch)
{
    real colvector eta_s_if, eta_t_if, log_p_if, xp, xa, xm
    real scalar n_req, n_t, j, q, k, c_s, c_t, b_s, b_t, bst_s, bst_t, d_s, has_int
    real scalar pp, ps, pt, qc, act, rho_j, a_k

    n_req = rows(ch.unit_ids)
    n_t = pb.n_targets
    pair_prob_if_ch = J(n_req, n_t, 0)
    mult_if_ch = J(n_req, n_t, 0)
    if (n_req == 0 | n_t == 0) return

    d_s = fit.spec.d_s
    has_int = fit.spec.has_intercept

    for (j = 1; j <= n_t; j = j + 1) {
        q = pb.target_pc[j]
        k = pb.target_cohort[j]
        c_s = fit.map.pair_cohorts[q, 6]
        c_t = fit.map.pair_cohorts[q, 7]
        b_s = fit.cell_block[c_s]
        b_t = fit.cell_block[c_t]
        bst_s = fit.beta_start[b_s]
        bst_t = fit.beta_start[b_t]

        pp = pb.pair_probabilities[j]
        ps = pb.p_s[j]
        pt = pb.p_t[j]
        qc = max((pp, 1e-6))
        act = pb.active_mask[j]
        rho_j = pb.target_rho[j]
        a_k = fit.a_U[k]

        eta_s_if = J(n_req, 1, 0)
        if (fit.spec.d > 0) {
            if (has_int) eta_s_if = eta_s_if + ch.xi_alpha[, c_s]
            if (d_s > 0)
                eta_s_if = eta_s_if + ch.xi_beta[, (bst_s - fit.n_alpha)..(bst_s - fit.n_alpha + d_s - 1)] * pb.target_z_slope[j, ]'
        }

        if (fit.selection_model == "pair_logit") {
            log_p_if = (1 - pp) :* eta_s_if
        }
        else {
            eta_t_if = J(n_req, 1, 0)
            if (fit.spec.d > 0) {
                if (has_int) eta_t_if = eta_t_if + ch.xi_alpha[, c_t]
                if (d_s > 0)
                    eta_t_if = eta_t_if + ch.xi_beta[, (bst_t - fit.n_alpha)..(bst_t - fit.n_alpha + d_s - 1)] * pb.target_z_slope[j, ]'
            }
            log_p_if = (1 - ps) :* eta_s_if + (1 - pt) :* eta_t_if
        }

        xp = pp :* log_p_if
        xa = ch.cohort_a_if[, k]
        xm = xa :/ (rho_j * qc) - (a_k / (rho_j * (qc ^ 2))) :* (act :* xp)
        pair_prob_if_ch[, j] = xp
        mult_if_ch[, j] = xm
    }
}

/* Stream original-unit IF covariance across chunks of 1..N without allocating
   the full N x K_total IF matrix. */
real matrix cdid_sel_stream_cov(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_reader scalar rdr, struct cdid_sel_probes scalar pb,
    struct cdid_sel_moments scalar acc)
{
    struct cdid_sel_if_chunk scalar ch
    real matrix cov_mat, pp_ch, mult_ch, mom_ch, psi_ch
    real colvector ch_ids
    real scalar n_tot, s_pos, e_pos

    n_tot = fit.total_params + rows(fit.cohorts) + 2 * pb.n_targets + acc.n_moments
    cov_mat = J(n_tot, n_tot, 0)
    if (n_tot == 0) return(cov_mat)

    for (s_pos = 1; s_pos <= fit.N; s_pos = s_pos + rdr.max_chunk_rows) {
        e_pos = min((fit.N, s_pos + rdr.max_chunk_rows - 1))
        ch_ids = (s_pos::e_pos)
        ch = cdid_sel_eval_if_chunk(fit, rdr, ch_ids, 1)
        pp_ch = J(rows(ch_ids), 0, 0)
        mult_ch = J(rows(ch_ids), 0, 0)
        if (pb.n_targets > 0)
            _cdid_sel_probe_if_chunk(fit, pb, ch, pp_ch, mult_ch)
        mom_ch = (acc.n_moments > 0 ? cdid_sel_eval_loading(fit, rdr, acc, ch_ids) : J(rows(ch_ids), 0, 0))
        psi_ch = (ch.theta_if, ch.cohort_a_if, pp_ch, mult_ch, mom_ch)
        cov_mat = cov_mat + quadcross(psi_ch, psi_ch) / (fit.N ^ 2)
    }
    return(cov_mat)
}

/* Standalone explicit-table selection model fit + IF + probe + moment loading +
   streamed covariance driver (mirrors R's evaluate_selection_table_if). */
struct cdid_sel_table_out scalar cdid_sel_table_if(real matrix Z,
    real scalar has_intercept, real scalar int_col,
    real colvector cohort_index, real colvector cohorts_in,
    real colvector all_times_in, real colvector w,
    real colvector A_in, real colvector B_in, real colvector U_in,
    real colvector rho_in, struct cdid_sel_decl scalar decl,
    real matrix date_responses, real matrix pair_responses,
    real colvector pair_s_idx, real colvector pair_t_idx,
    real matrix pair_specs, real colvector probe_cohorts,
    real matrix Z_probe, real colvector rho_probe,
    real matrix moment_weights, real scalar eps, real scalar maxit,
    real scalar max_chunk_rows, real scalar max_work)
{
    struct cdid_sel_table_out scalar out
    struct cdid_sel_reader scalar rdr
    struct cdid_sel_cell_map scalar map
    struct cdid_sel_probes scalar pb
    struct cdid_sel_up_map scalar um
    struct cdid_sel_moments scalar acc
    struct cdid_sel_if_chunk scalar ch
    real colvector all_ids, u_idx
    real matrix mw_sub
    real scalar n_mom, pc, k

    rdr = cdid_sel_reader_table(Z, has_intercept, int_col, cohort_index, cohorts_in,
        all_times_in, w, A_in, B_in, U_in, rho_in, decl, date_responses,
        pair_responses, pair_s_idx, pair_t_idx, max_chunk_rows)
    map = cdid_sel_partition(rdr, pair_specs)
    out.budget = cdid_sel_budget_check(rdr, map, max_work)
    out.fit = cdid_sel_fit(rdr, map, eps, maxit)

    out.alpha = out.fit.alpha
    out.beta = out.fit.beta
    out.coefficients = out.fit.theta
    out.cohort_a = out.fit.a_U
    out.mass_U = out.fit.mass_U

    all_ids = (1::out.fit.N)
    ch = cdid_sel_eval_if_chunk(out.fit, rdr, all_ids, 1)
    out.theta_if = ch.theta_if
    out.cohort_a_if = ch.cohort_a_if

    pb = cdid_sel_pred_probes(out.fit, rdr, probe_cohorts, Z_probe, rho_probe)
    out.pair_probabilities = pb.pair_probabilities
    out.multipliers = pb.multipliers
    out.pair_prob_if = J(out.fit.N, pb.n_targets, 0)
    out.multiplier_if = J(out.fit.N, pb.n_targets, 0)
    if (pb.n_targets > 0)
        _cdid_sel_probe_if_chunk(out.fit, pb, ch, out.pair_prob_if, out.multiplier_if)

    n_mom = cols(moment_weights)
    acc = cdid_sel_moments_init(out.fit, n_mom)
    if (n_mom > 0) {
        if (rows(moment_weights) != out.fit.N)
            _cdid_sel_fail("input", "moment_weights must have N rows.", 498)
        um = cdid_sel_up_map_build(out.fit)
        for (pc = 1; pc <= out.fit.n_pair_cohorts; pc = pc + 1) {
            k = out.fit.map.pair_cohorts[pc, 1]
            u_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (rdr.U :== 1))
            if (rows(u_idx) == 0) continue
            mw_sub = moment_weights[u_idx, ] :/ out.fit.N
            cdid_sel_add_moments(acc, out.fit, rdr, um, pc, u_idx, mw_sub)
        }
    }
    out.loading_theta = acc.G_theta
    out.loading_ak = acc.G_ak
    out.moment_if = (n_mom > 0 ? cdid_sel_eval_loading(out.fit, rdr, acc, all_ids) : J(out.fit.N, 0, 0))
    out.covariance = cdid_sel_stream_cov(out.fit, rdr, pb, acc)
    return(out)
}

end
