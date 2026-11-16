*! Baseline static-input chained IPW assembly with source GMM; qualified on StataNow19.5.
*! Reference: R cdid 1.1.0.0, commit .
*! Consumes cdid_population.mata, cdid_catalog.mata, cdid_nuisance.mata,
*! cdid_ipw_link.mata, and cdid_gmm.mata.
version 14.2
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_cipw(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
mata set matastrict on

struct cdid_chained_ipw_opts {
    string scalar est_method, gmm_weighting, pt_assumption
    string scalar selection, fix_weights, conventions
    real scalar   pt_start, dated_design
    real scalar   glm_epsilon, glm_maxit
    real scalar   condition_trigger, ridge_relative
    real scalar   chunk_size, max_matrix_entries, work_cap
    real scalar   max_cached_nuisances, materialize_moment_covariance, reduce_structural_redundancy
    pointer() scalar within_cov
}

struct cdid_chained_ipw_postings {
    real scalar    n_A, n_dates, n_obs_a
    real colvector date_start, date_end, date_count
    real colvector post_a, post_y
    real matrix    obs_A, Y_A
}

struct cdid_chained_ipw_cache_entry {
    real scalar              active, op_index, stamp
    real colvector           pair_a
    struct cdid_ipw_nuisance scalar nu
}

struct cdid_chained_ipw_link_meta {
    real scalar    row_index, tp_index, op_index, pool_index
    real scalar    g_index, s_index, t_index, ref_index
    real scalar    g_label, s_label, t_label, ref_label
    real scalar    s_role, t_role, s_qualified, t_qualified
    real scalar    from_vertex, to_vertex, start_param, end_param
    real scalar    link_type, component_id, zero_loading
    real scalar    n_cohorts, p_covariates
    real scalar    n_pair_all, n_pair_risk, n_treated_obs, n_control_obs
    real scalar    estimate, treated_mean, control_mean
    real scalar    sum_a, sum_b, mean_a, mean_b
    real scalar    iterations, deviance, score_max, reused_nuisance
    string scalar  functional_key, row_fit_kind
    real colvector cohort_indices, group_a_count, group_a_mass
    real colvector pair_count, pair_mass, q, rates, rates_c, rates_active
    real colvector beta, G_beta, G_r, loading_beta, loading_rate
    real colvector rate_coef_s0, rate_coef_s1
}

/* Outcome-free optimal coordinates. Identity intentionally keeps all rows.
   E is implicit on retained rows and sparse signed tree paths on dropped rows.
   Every class has identical operator/static design/weights/GLM controls and
   exactly equal observed RISK-unit indices; no outcome enters this decision. */
struct cdid_cipw_row_basis {
    real scalar M, L, n_classes, n_paths, workspace_entries
    real colvector retained, dropped, class_id, class_rep, fit_order
    real colvector path_start, path_count, e_row, e_col, e_value
    real matrix G
}

struct cdid_chained_ipw_result {
    string scalar  status, global_reason, moment_covariance_status, gmm_status
    string scalar  control_group, pt_assumption, weighting
    real scalar    pt_all, pt_start, pt_start_index, anticipation
    real scalar    N_F, n_panel, n_A, n_B, n_singleton, n_unobserved, n_frame_only
    real scalar    n_dates, n_cohorts, p_covariates
    real scalar    M_links, n_operators, n_pools
    struct cdid_cipw_row_basis scalar row_basis
    real scalar    n_vertices, n_params_unreduced, n_params_retained, n_nullity, n_components
    real scalar    n_zero_loading_links
    real scalar    n_targets, n_estimable_targets, n_unavailable_targets
    real scalar    n_nuisance_fits, n_nuisance_reused
    real scalar    chunk_size, effective_chunk_size, n_batches
    real scalar    workspace_entries_bound, workspace_entries_limit, work_ops, work_cap
    real scalar    no_extra_calendar_window, no_extra_represented_restriction, no_extra_restrictions

    real colvector calendar, cohorts
    real colvector A_to_F, F_to_A, unit_cohort, unit_cohort_known, unit_weight

    real colvector moments, moment_if_mean, moment_fitted, moment_residual
    real matrix    moment_table
    struct cdid_chained_ipw_link_meta colvector links
    real matrix    Vr_estimator, Vr_analytical

    real colvector kept_param_indices, dropped_gauge_vertices, dropped_gauge_params
    real matrix    A_retained, H_estimable
    struct cdid_gmm_result scalar gmm
    real scalar    criterion

    real colvector target_g, target_t, target_g_index, target_t_index
    real colvector target_cat_supported, target_cat_ntouch, target_cat_reason
    real colvector target_in_design, target_vertex, target_param
    real colvector target_design_status, target_available
    string colvector target_cat_reason_name, target_design_status_name, target_unavailable_reason
    real colvector estimable_target_indices, unavailable_target_indices
    real colvector att, se
    real matrix    target_ledger

    real matrix    target_if
    real rowvector target_if_mean
    real matrix    V_estimator, V_analytical
    real matrix    V_estimator_estimable, V_analytical_estimable
}

/* Owned replay state: pointers remain local/stable throughout the solver call.
   No user supplied covariance/provenance boolean enters the internal Gram. */
struct cdid_cipw_replay_context {
    pointer() scalar population, catalog, postings, cohorts, operators, X, w, links
}

void _cdid_cipw_set_lasterr(string scalar cls)
{
    pointer() scalar p
    if ((p = findexternal("cdid_chained_ipw_lasterr")) == NULL) {
        p = crexternal("cdid_chained_ipw_lasterr")
    }
    if (p != NULL) *p = cls
}

string scalar cdid_chained_ipw_last_error()
{
    pointer() scalar p
    if ((p = findexternal("cdid_chained_ipw_lasterr")) == NULL) return("")
    if (eltype(*p) != "string" | rows(*p) != 1 | cols(*p) != 1) return("")
    return(*p)
}

void cdid_chained_ipw_clear_error()
{
    _cdid_cipw_set_lasterr("")
}

void _cdid_cipw_fail(string scalar cls, string scalar msg, real scalar rc)
{
    _cdid_cipw_set_lasterr(cls)
    errprintf("%s\n", msg)
    exit(rc)
}

real scalar _cdid_cipw_any(real matrix x)
{
    return(sum(x :!= 0) > 0)
}

real scalar _cdid_cipw_all(real matrix x)
{
    return(sum(x :!= 0) == rows(x) * cols(x))
}

real colvector _cdid_cipw_which(real colvector x)
{
    if (rows(x) == 0) return(J(0, 1, 0))
    return(vec(selectindex(x :!= 0)))
}

real colvector _cdid_cipw_seq(real scalar a, real scalar b)
{
    if (b < a) return(J(0, 1, 0))
    return((a::b))
}

real scalar _cdid_cipw_cells(real scalar a, real scalar b, real scalar cap, real scalar total)
{
    real scalar z
    if (a < 0 | b < 0 | a >= . | b >= . | a != floor(a) | b != floor(b))
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "Workspace dimensions must be finite nonnegative integers.", 498)
    if (a > 0 & b > (2^53 - 1) / a)
        _cdid_cipw_fail("cdid_chained_ipw_workspace_error", "Workspace dimension product exceeds exact integer range.", 908)
    z = a * b
    if (z > cap - total)
        _cdid_cipw_fail("cdid_chained_ipw_workspace_error", "Required workspace exceeds max_matrix_entries.", 908)
    return(total + z)
}

real scalar _cdid_cipw_work_add(real scalar ops, real scalar work_cap, real scalar total)
{
    if (ops < 0 | ops >= .)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "Work operation increment must be finite and nonnegative.", 498)
    if (ops > (2^53 - 1) - total)
        _cdid_cipw_fail("cdid_chained_ipw_work_cap_error", "Accumulated work operations exceed exact integer range.", 908)
    if (work_cap < . & total + ops > work_cap)
        _cdid_cipw_fail("cdid_chained_ipw_work_cap_error", "Required estimator work exceeds explicit work_cap.", 908)
    return(total + ops)
}

real scalar _cdid_cipw_lb(real colvector vec, real scalar left, real scalar right, real scalar target)
{
    real scalar lo, hi, mid, ans
    lo = left
    hi = right
    ans = right + 1
    while (lo <= hi) {
        mid = floor(lo + (hi - lo) / 2)
        if (vec[mid] >= target) {
            ans = mid
            hi = mid - 1
        }
        else lo = mid + 1
    }
    return(ans)
}

real scalar _cdid_cipw_ub(real colvector vec, real scalar left, real scalar right, real scalar target)
{
    real scalar lo, hi, mid, ans
    lo = left
    hi = right
    ans = left - 1
    while (lo <= hi) {
        mid = floor(lo + (hi - lo) / 2)
        if (vec[mid] <= target) {
            ans = mid
            lo = mid + 1
        }
        else hi = mid - 1
    }
    return(ans)
}

struct cdid_chained_ipw_opts scalar cdid_chained_ipw_default_opts()
{
    struct cdid_chained_ipw_opts scalar o
    o.est_method = "ipw"
    o.gmm_weighting = "identity"
    o.pt_assumption = "post"
    o.pt_start = .
    o.selection = "none"
    o.fix_weights = ""
    o.dated_design = 0
    o.conventions = "stata-port"
    o.glm_epsilon = 1e-12
    o.glm_maxit = 50
    o.condition_trigger = 1e-10
    o.ridge_relative = 1e-4
    o.chunk_size = 1024
    o.max_matrix_entries = 2^27
    o.work_cap = .
    o.max_cached_nuisances = 16
    o.within_cov = NULL
    o.reduce_structural_redundancy = 1
    o.materialize_moment_covariance = 0
    return(o)
}

void _cdid_cipw_check_opts(struct cdid_chained_ipw_opts scalar o)
{
    if (o.est_method != "ipw")
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_method", "Baseline static-input IPW estimator requires est_method == ipw.", 498)
    if (o.gmm_weighting != "identity" & o.gmm_weighting != "optimal" & o.gmm_weighting != "two_step" & o.gmm_weighting != "within")
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_weighting", "gmm_weighting must be identity, optimal, two_step or within.", 498)
    if (o.selection != "" & o.selection != "none")
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_selection", "Baseline static-input IPW estimator requires selection == none.", 498)
    if (o.dated_design != 0)
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_dated_design", "Baseline static-input IPW entry refuses dated covariate designs.", 498)
    if (o.fix_weights != "" & o.fix_weights != "static" & o.fix_weights != "first_period")
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_weights", "Baseline static-input IPW entry requires static or first_period fixed weights.", 498)
    if (o.pt_assumption != "post" & o.pt_assumption != "all" & o.pt_assumption != "ctrl")
        _cdid_cipw_fail("cdid_design_pt_assumption_error", "pt_assumption must be post, ctrl, or all.", 498)
    if (o.pt_assumption == "post" & o.pt_start < .)
        _cdid_cipw_fail("cdid_design_pt_start_error", "pt_start must be missing under pt_assumption post.", 498)
    if (o.pt_assumption == "all" & o.pt_start >= .)
        _cdid_cipw_fail("cdid_design_pt_start_error", "pt_start must be an explicit stored calendar date under pt_assumption all.", 498)
    if (o.conventions == "" | strpos(o.conventions, ",") | strpos(o.conventions, char(34)))
        _cdid_cipw_fail("cdid_catalog_input_error", "conventions must be a nonempty label without commas or quotes.", 498)
    if (o.glm_epsilon >= . | o.glm_epsilon <= 0 | o.glm_maxit >= . | o.glm_maxit < 1 |
        o.glm_maxit != floor(o.glm_maxit) | o.glm_maxit > 2^31 - 1)
        _cdid_cipw_fail("cdid_ipw_control_error", "Invalid GLM solver controls glm_epsilon / glm_maxit.", 498)
    if (o.condition_trigger >= . | o.condition_trigger <= 0 | o.condition_trigger >= 1 |
        o.ridge_relative >= . | o.ridge_relative <= 0)
        _cdid_cipw_fail("cdid_gmm_control_error", "Invalid GMM tuning controls condition_trigger / ridge_relative.", 498)
    if (o.chunk_size >= . | o.chunk_size < 1 | o.chunk_size != floor(o.chunk_size) | o.chunk_size > 2^53 - 1)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "chunk_size must be a positive finite integer.", 498)
    if (o.max_matrix_entries >= . | o.max_matrix_entries < 1 |
        o.max_matrix_entries != floor(o.max_matrix_entries) | o.max_matrix_entries > 2^53 - 1)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "max_matrix_entries must be a positive finite integer.", 498)
    if (o.work_cap < .) {
        if (o.work_cap > .)
            _cdid_cipw_fail("cdid_chained_ipw_bound_error", "work_cap cannot be an extended missing value.", 498)
        if (o.work_cap < 1 | o.work_cap != floor(o.work_cap) | o.work_cap > 2^53 - 1)
            _cdid_cipw_fail("cdid_chained_ipw_bound_error", "When supplied, work_cap must be a positive finite integer.", 498)
    }
    else if (o.work_cap > .)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "work_cap cannot be an extended missing value.", 498)
    if (o.reduce_structural_redundancy!=0 & o.reduce_structural_redundancy!=1)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error","reduce_structural_redundancy must be 0 or1.",498)
    if (o.materialize_moment_covariance!=0 & o.materialize_moment_covariance!=1)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error","materialize_moment_covariance must be 0 or1.",498)
    if (o.max_cached_nuisances >= . | o.max_cached_nuisances < 0 |
        o.max_cached_nuisances != floor(o.max_cached_nuisances) | o.max_cached_nuisances > 2^31 - 1)
        _cdid_cipw_fail("cdid_chained_ipw_bound_error", "max_cached_nuisances must be a nonnegative finite integer.", 498)
}

/* Extract the sorted known cohort dictionary from a prepared population P.
   Units with cohort_known == 0 (unknown frame-only units) never contribute
   a never-treated (.) entry; only cohort_known == 1 with cohort == . does. */
real colvector cdid_cipw_cohorts(struct cdid_population scalar P)
{
    real colvector known_idx, kc, fin_idx, out
    real scalar has_nt
    if (P.N < 1 | rows(P.cohort) != P.N | rows(P.cohort_known) != P.N)
        _cdid_cipw_fail("cdid_chained_ipw_population_error", "Invalid prepared population structure.", 498)
    known_idx = _cdid_cipw_which(P.cohort_known :!= 0)
    if (rows(known_idx) == 0)
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Population has zero units with known cohort.", 459)
    kc = P.cohort[known_idx]
    if (_cdid_cipw_any(kc :> .))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Extended missing values are invalid in cohorts.", 459)
    fin_idx = _cdid_cipw_which(kc :< .)
    has_nt = _cdid_cipw_any(kc :== .)
    out = J(0, 1, .)
    if (rows(fin_idx) > 0) out = uniqrows(sort(kc[fin_idx], 1))
    if (has_nt) out = out \ .
    return(out)
}

real colvector cdid_cipw_targets(struct cdid_population scalar P)
{
    real colvector coh, fin_idx
    coh = cdid_cipw_cohorts(P)
    fin_idx = _cdid_cipw_which(coh :< .)
    if (rows(fin_idx) == 0)
        _cdid_cipw_fail("cdid_chained_ipw_target_error", "At least one finite treated target cohort is required.", 459)
    return(coh[fin_idx])
}

real matrix cdid_cipw_post_targets(real colvector calendar, real colvector target_cohorts)
{
    real scalar nJ, nT, j, d, cnt, at
    real matrix out
    nJ = rows(target_cohorts)
    nT = rows(calendar)
    if (nJ < 1 | cols(target_cohorts) != 1 | missing(target_cohorts))
        _cdid_cipw_fail("cdid_chained_ipw_target_error", "target_cohorts must be a nonempty column of finite cohorts.", 498)
    if (nT < 1 | cols(calendar) != 1 | missing(calendar))
        _cdid_cipw_fail("cdid_chained_ipw_calendar_error", "calendar must be a nonempty column of finite dates.", 498)
    cnt = 0
    for (j = 1; j <= nJ; j++) {
        for (d = 1; d <= nT; d++) {
            if (calendar[d] >= target_cohorts[j]) cnt++
        }
    }
    if (cnt == 0) return(J(0, 2, .))
    out = J(cnt, 2, .)
    at = 0
    for (j = 1; j <= nJ; j++) {
        for (d = 1; d <= nT; d++) {
            if (calendar[d] >= target_cohorts[j]) {
                at++
                out[at, .] = (target_cohorts[j], calendar[d])
            }
        }
    }
    return(out)
}

struct cdid_pairs scalar cdid_cipw_pairs_opts(struct cdid_population scalar P,
    real colvector cohorts, struct cdid_catalog_opts scalar cat_opts)
{
    real scalar u, nt
    real colvector d_idx, c_idx, unit_c
    real matrix obs
    nt = _cdid_check_cohorts(cohorts, 0)
    if (P.R == 0) {
        obs = J(0, 4, .)
        return(cdid_observed_pairs_opts(obs, P.calendar, cohorts, cat_opts))
    }
    if (any(P.ids.row_unit :< 1) | any(P.ids.row_unit :> P.N))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Every panel row unit must have a known cohort in F.", 459)
    if (any(P.cohort_known[P.ids.row_unit] :!= 1))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Every panel row unit must have a known cohort in F.", 459)
    unit_c = _cdid_cohort_index_vec(cohorts, nt, P.cohort, P.cohort_known)
    d_idx = _cdid_population_lookup(P.date, P.calendar)
    c_idx = unit_c[P.ids.row_unit]
    if (any(d_idx :== 0) | any(c_idx :== 0))
        _cdid_cipw_fail("cdid_chained_ipw_dictionary_error", "Panel row date or cohort is absent from supplied dictionaries.", 459)
    obs = (P.ids.row_unit, d_idx, c_idx, P.y)
    return(cdid_observed_pairs_opts(obs, P.calendar, cohorts, cat_opts))
}

struct cdid_pairs scalar cdid_chained_ipw_observed_pairs(struct cdid_population scalar P,
    | real colvector cohorts)
{
    real colvector coh
    coh = (args() >= 2 ? cohorts : cdid_cipw_cohorts(P))
    return(cdid_cipw_pairs_opts(P, coh, cdid_catalog_default_opts()))
}

/* Authoritative R default qualification (R/chained_ipw_pipeline.R lines 5495..5508):
   qualify ALL conditional causal endpoints (anticipation_conditional == 3 or
   post_conditional == 5) present in retained comparison catalogue rows. */
real matrix cdid_cipw_qualified(struct cdid_cat scalar C)
{
    real scalar nr, ns, nt
    real colvector s_idx, t_idx, dq_gi, dq_ti
    real matrix pairs_idx, out
    nr = rows(C.r_tp)
    if (nr == 0) return(J(0, 2, .))
    s_idx = _cdid_cipw_which((C.r_srole :== 3) :| (C.r_srole :== 5))
    t_idx = _cdid_cipw_which((C.r_trole :== 3) :| (C.r_trole :== 5))
    ns = rows(s_idx)
    nt = rows(t_idx)
    if (ns + nt == 0) return(J(0, 2, .))
    dq_gi = J(0, 1, 0)
    dq_ti = J(0, 1, 0)
    if (ns > 0) {
        dq_gi = dq_gi \ C.r_g[s_idx]
        dq_ti = dq_ti \ C.r_s[s_idx]
    }
    if (nt > 0) {
        dq_gi = dq_gi \ C.r_g[t_idx]
        dq_ti = dq_ti \ C.r_t[t_idx]
    }
    pairs_idx = uniqrows(sort((dq_gi, dq_ti), (1, 2)))
    out = (C.cohorts[pairs_idx[., 1]], C.calendar[pairs_idx[., 2]])
    return(out)
}

/* Authoritative R design target filter (R/chained_ipw_pipeline.R lines 5514..5518):
   only row-bearing cohorts' requested targets enter the design target grid;
   every unsupported requested cell stays in the estimator target ledger. */
real matrix cdid_chained_ipw_design_targets(struct cdid_cat scalar C, real matrix requested_targets)
{
    real scalar q, i, gi, nt
    real colvector row_g, keep, keep_idx
    q = rows(requested_targets)
    if (q == 0 | rows(C.r_tp) == 0) return(J(0, 2, .))
    if (cols(requested_targets) != 2)
        _cdid_cipw_fail("cdid_chained_ipw_target_error", "requested_targets must have 2 columns (g, t).", 498)
    nt = C.nt_index
    row_g = uniqrows(sort(C.r_g, 1))
    keep = J(q, 1, 0)
    for (i = 1; i <= q; i++) {
        gi = _cdid_cohort_index(C.cohorts, nt, requested_targets[i, 1])
        if (gi > 0) keep[i] = (_cdid_bsearch(row_g, gi) > 0)
    }
    keep_idx = _cdid_cipw_which(keep)
    if (rows(keep_idx) == 0) return(J(0, 2, .))
    return(requested_targets[keep_idx, .])
}

struct cdid_design scalar cdid_chained_ipw_build_design(struct cdid_cat scalar C,
    real matrix requested_targets, string scalar pt_assumption, real scalar pt_start)
{
    real matrix req, dq, tpt
    req = requested_targets
    if (rows(req) == 0 & rows(C.rq_g) > 0) {
        req = (C.cohorts[C.rq_g], C.calendar[C.rq_t])
    }
    dq = cdid_cipw_qualified(C)
    tpt = cdid_chained_ipw_design_targets(C, req)
    return(cdid_catalog_design_ex(C, dq, tpt, pt_assumption, pt_start))
}

struct cdid_chained_ipw_postings scalar _cdid_cipw_build_postings(
    struct cdid_population scalar P, real colvector calendar)
{
    struct cdid_chained_ipw_postings scalar post
    real scalar nO, di, cnt_a, nT, k
    real colvector keep_r, raw_d, raw_a, raw_y, ord, a_obs, keep
    nT = rows(calendar)
    post.n_A = P.n_A
    post.n_dates = nT
    post.date_count = J(nT, 1, 0)
    post.date_start = J(nT, 1, 1)
    post.date_end = J(nT, 1, 0)
    post.post_a = J(0, 1, 0)
    post.post_y = J(0, 1, .)
    post.obs_A = J(0, 0, 0)
    post.Y_A = J(0, 0, .)
    post.n_obs_a = 0
    nO = rows(P.observed_order)
    if (nO == 0 | P.n_A == 0) return(post)

    a_obs = P.F_to_A[P.ids.row_unit[P.observed_order]]
    keep = _cdid_cipw_which(a_obs :> 0)
    cnt_a = rows(keep)
    if (cnt_a == 0) return(post)
    keep_r = P.observed_order[keep]
    raw_a = a_obs[keep]
    raw_d = _cdid_population_lookup(P.date[keep_r], calendar)
    if (any(raw_d :== 0))
        _cdid_cipw_fail("cdid_chained_ipw_calendar_error", "Observed A-unit date absent from calendar.", 459)
    raw_y = P.y[keep_r]
    for (di = 1; di <= nT; di++) post.date_count[di] = sum(raw_d :== di)
    ord = order((raw_d, raw_a), (1, 2))
    post.post_a = raw_a[ord]
    post.post_y = raw_y[ord]
    post.n_obs_a = cnt_a
    post.date_end = runningsum(post.date_count)
    post.date_start = post.date_end - post.date_count :+ 1
    if (cnt_a > 1) {
        if (any((raw_d[ord[|2 \ cnt_a|]] :== raw_d[ord[|1 \ cnt_a - 1|]]) :&
            (post.post_a[|2 \ cnt_a|] :<= post.post_a[|1 \ cnt_a - 1|])))
            _cdid_cipw_fail("cdid_chained_ipw_postings_error", "Calendar postings within date must be strictly increasing in A index.", 504)
    }
    if (P.n_A > 0 & nT > 0 & P.n_A * nT <= 2000000 & max(abs(post.post_y)) < 1e300) {
        post.obs_A = J(P.n_A, nT, 0)
        post.Y_A = J(P.n_A, nT, 0)
        for (di = 1; di <= nT; di++) {
            if (post.date_count[di] > 0) {
                keep = post.post_a[|post.date_start[di] \ post.date_end[di]|]
                post.obs_A[keep, di] = J(rows(keep), 1, 1)
                post.Y_A[keep, di] = post.post_y[|post.date_start[di] \ post.date_end[di]|]
            }
        }
    }
    return(post)
}

void _cdid_cipw_query_pair(struct cdid_chained_ipw_postings scalar post,
    real scalar s_idx, real scalar t_idx, real scalar a_lo, real scalar a_hi,
    real colvector pair_a, real colvector pair_dy)
{
    real scalar ls, hs, lt, ht, min_a, max_a, is_ptr, it_ptr, cap_hit, nhit, as_val, at_val
    real colvector tmp_a, tmp_dy
    pair_a = J(0, 1, 0)
    pair_dy = J(0, 1, .)
    if (a_lo > a_hi | s_idx < 1 | t_idx < 1 | s_idx > post.n_dates | t_idx > post.n_dates) return
    if (post.date_count[s_idx] == 0 | post.date_count[t_idx] == 0) return

    if (rows(post.obs_A) == post.n_A & a_lo >= 1 & a_hi <= post.n_A) {
        if (a_lo == 1 & a_hi == post.n_A) {
            tmp_a = vec(selectindex(post.obs_A[., s_idx] :& post.obs_A[., t_idx]))
        }
        else {
            tmp_a = vec(selectindex(post.obs_A[|a_lo, s_idx \ a_hi, s_idx|] :& post.obs_A[|a_lo, t_idx \ a_hi, t_idx|]))
            if (rows(tmp_a) > 0) tmp_a = tmp_a :+ (a_lo - 1)
        }
        if (rows(tmp_a) == 0) return
        pair_a = tmp_a
        pair_dy = post.Y_A[pair_a, t_idx] :- post.Y_A[pair_a, s_idx]
        if (missing(pair_dy))
            _cdid_cipw_fail("cdid_ipw_outcome_error", "Pair outcome difference overflow.", 504)
        return
    }

    ls = _cdid_cipw_lb(post.post_a, post.date_start[s_idx], post.date_end[s_idx], a_lo)
    hs = _cdid_cipw_ub(post.post_a, post.date_start[s_idx], post.date_end[s_idx], a_hi)
    if (ls > hs) return

    lt = _cdid_cipw_lb(post.post_a, post.date_start[t_idx], post.date_end[t_idx], a_lo)
    ht = _cdid_cipw_ub(post.post_a, post.date_start[t_idx], post.date_end[t_idx], a_hi)
    if (lt > ht) return

    min_a = max((post.post_a[ls], post.post_a[lt]))
    max_a = min((post.post_a[hs], post.post_a[ht]))
    if (min_a > max_a) return

    if (post.post_a[ls] < min_a) ls = _cdid_cipw_lb(post.post_a, ls, hs, min_a)
    if (post.post_a[hs] > max_a) hs = _cdid_cipw_ub(post.post_a, ls, hs, max_a)
    if (ls > hs) return

    if (post.post_a[lt] < min_a) lt = _cdid_cipw_lb(post.post_a, lt, ht, min_a)
    if (post.post_a[ht] > max_a) ht = _cdid_cipw_ub(post.post_a, lt, ht, max_a)
    if (lt > ht) return

    cap_hit = min((hs - ls + 1, ht - lt + 1))
    tmp_a = J(cap_hit, 1, 0)
    tmp_dy = J(cap_hit, 1, .)
    nhit = 0
    is_ptr = ls
    it_ptr = lt
    while (is_ptr <= hs & it_ptr <= ht) {
        as_val = post.post_a[is_ptr]
        at_val = post.post_a[it_ptr]
        if (as_val == at_val) {
            nhit++
            tmp_a[nhit] = as_val
            tmp_dy[nhit] = post.post_y[it_ptr] - post.post_y[is_ptr]
            is_ptr++
            it_ptr++
        }
        else if (as_val < at_val) is_ptr++
        else it_ptr++
    }
    if (nhit == 0) return
    pair_a = tmp_a[|1 \ nhit|]
    pair_dy = tmp_dy[|1 \ nhit|]
    if (missing(pair_dy))
        _cdid_cipw_fail("cdid_ipw_outcome_error", "Pair outcome difference overflow.", 504)
}

/* Only posting membership is read. In particular no subtraction of Y, logit,
   covariance eigenvalue or empirical rank can influence the structural map. */
real colvector _cdid_cipw_pair_indices(
    struct cdid_chained_ipw_postings scalar post, real scalar s, real scalar t)
{
    real scalar a, b, ae, be, n, cap
    real colvector hit
    if (rows(post.obs_A) == post.n_A) {
        if (post.date_count[s] == 0 | post.date_count[t] == 0) return(J(0, 1, 0))
        return(vec(selectindex(post.obs_A[., s] :& post.obs_A[., t])))
    }
    a=post.date_start[s]; ae=post.date_end[s]
    b=post.date_start[t]; be=post.date_end[t]
    cap=min((post.date_count[s],post.date_count[t]))
    hit=J(cap,1,0); n=0
    while(a<=ae & b<=be) {
        if(post.post_a[a]==post.post_a[b]) {
            n++;hit[n]=post.post_a[a];a++;b++
        }
        else if(post.post_a[a]<post.post_a[b])a++
        else b++
    }
    if(n==0)return(J(0,1,0))
    return(hit[|1\n|])
}

real colvector _cdid_cipw_risk_indices(
    struct cdid_chained_ipw_postings scalar post, struct cdid_cat scalar C,
    real colvector ac, real matrix om, real scalar ell)
{
    real colvector hit, keep
    if (rows(post.obs_A) == post.n_A) {
        if (post.date_count[C.r_s[ell]] == 0 | post.date_count[C.r_t[ell]] == 0) return(J(0, 1, 0))
        return(vec(selectindex(post.obs_A[., C.r_s[ell]] :& post.obs_A[., C.r_t[ell]] :& (om[ac, C.r_op[ell]] :> 0))))
    }
    hit=_cdid_cipw_pair_indices(post,C.r_s[ell],C.r_t[ell])
    if(rows(hit)==0)return(hit)
    keep=_cdid_cipw_which(om[ac[hit],C.r_op[ell]]:>0)
    if(rows(keep)==0)return(J(0,1,0))
    return(hit[keep])
}

/* Exact sparse E products; the dense M-by-L map is never required. */
real matrix _cdid_cipw_expand(struct cdid_cipw_row_basis scalar B,real matrix Z)
{
    real scalar h,k,i
    real matrix out
    if(rows(Z)!=B.L)_cdid_cipw_fail("cdid_ipw_reduction_error","Primitive expansion dimension mismatch.",3200)
    out=J(B.M,cols(Z),0)
    out[B.retained,.]=Z
    for(h=1;h<=rows(B.dropped);h++) {
        i=B.dropped[h]
        for(k=B.path_start[i];k<B.path_start[i]+B.path_count[i];k++)
            out[i,.]=out[i,.]+B.e_value[k]*Z[B.e_col[k],.]
    }
    return(out)
}

void _cdid_cipw_check_design(struct cdid_cipw_row_basis scalar B,real matrix A)
{
    real scalar h,i,k
    real rowvector row
    for(h=1;h<=rows(B.dropped);h++) {
        i=B.dropped[h];row=J(1,cols(A),0)
        for(k=B.path_start[i];k<B.path_start[i]+B.path_count[i];k++)
            row=row+B.e_value[k]*A[B.retained[B.e_col[k]],.]
        if(!all(row:==A[i,.]))
            _cdid_cipw_fail("cdid_ipw_reduction_error","Structural map does not reconstruct the complete retained design.",459)
    }
}

struct cdid_cipw_row_basis scalar _cdid_cipw_row_basis(
    struct cdid_cat scalar C,struct cdid_design scalar D,
    struct cdid_chained_ipw_postings scalar post,real colvector ac,real matrix om,
    struct cdid_chained_ipw_opts scalar opts,real scalar entries,real scalar work)
{
    struct cdid_cipw_row_basis scalar B
    struct _cdid_forest scalar F
    real scalar M,ell,c,k,hit,nc,nn,nv,i,u,v,d,at,off,top,npath,a,b,h,j,q,sgn
    real scalar ws,path_ws
    real colvector risk,reps,cl,mem,dates,from,to,vmap,keep,tree
    real colvector par,edge,sign,depth,rf,rt,stack,idx,ord,rowspec,full_to_ret
    real matrix paths
    pointer(real colvector) colvector rep_risk
    M=rows(C.r_tp)
    ws=_cdid_cipw_cells(6,M,opts.max_matrix_entries,entries)
    B.M=M;B.L=M;B.n_classes=0;B.n_paths=0
    B.retained=(1::M);B.dropped=J(0,1,0);B.class_id=J(M,1,0)
    B.class_rep=J(0,1,0);B.fit_order=(1::M)
    B.path_start=J(M,1,0);B.path_count=J(M,1,0)
    B.e_row=J(0,1,0);B.e_col=J(0,1,0);B.e_value=J(0,1,0);B.G=J(0,0,0)
    B.workspace_entries=ws-entries
    if(opts.gmm_weighting=="identity" | !opts.reduce_structural_redundancy)return(B)
    /* Linear proof/graph scratch reserved BEFORE its arrays and each path
       payload count established BEFORE any path/Gram allocation. */
    ws=_cdid_cipw_cells(78,M,opts.max_matrix_entries,ws)
    ws=_cdid_cipw_cells(5,post.n_A,opts.max_matrix_entries,ws)
    reps=J(M,1,0);cl=J(M,1,0);nc=0;rep_risk=J(M,1,NULL)
    for(ell=1;ell<=M;ell++) {
        risk=_cdid_cipw_risk_indices(post,C,ac,om,ell)
        work=_cdid_cipw_work_add(post.n_A,opts.work_cap,work)
        if(rows(risk)==0)_cdid_cipw_fail("cdid_ipw_pair_support_error","Structural row has zero risk observers.",459)
        hit=0
        for(c=1;c<=nc;c++) {
            if(C.r_op[ell]!=C.r_op[reps[c]])continue
            work=_cdid_cipw_work_add(post.n_A,opts.work_cap,work)
            if(rows(*rep_risk[c])!=rows(risk))continue
            if((C.r_s[ell]==C.r_s[reps[c]] & C.r_t[ell]==C.r_t[reps[c]]) |
                (post.date_count[C.r_s[ell]]==post.n_A & post.date_count[C.r_t[ell]]==post.n_A &
                 post.date_count[C.r_s[reps[c]]]==post.n_A & post.date_count[C.r_t[reps[c]]]==post.n_A)) {hit=c;break;}
            if(all((*rep_risk[c]):==risk)) {hit=c;break;}
        }
        if(hit==0) {nc++;reps[nc]=ell;rep_risk[nc]=&(risk:+0);hit=nc;}
        cl[ell]=hit
    }
    B.n_classes=nc;B.class_id=cl;B.class_rep=reps[|1\nc|]
    par=J(2*M,1,0);edge=J(2*M,1,0);sign=J(2*M,1,0);depth=J(2*M,1,-1)
    rf=J(M,1,0);rt=J(M,1,0);keep=J(M,1,0);off=0
    for(c=1;c<=nc;c++) {
        mem=_cdid_cipw_which(cl:==c)
        mem=mem[order((C.r_s[mem],C.r_t[mem],mem),(1,2,3))]
        nn=rows(mem)
        dates=uniqrows(sort(C.r_s[mem]\C.r_t[mem],1));nv=rows(dates)
        from=J(nn,1,0);to=J(nn,1,0);vmap=J(nv,1,0)
        for(i=1;i<=nn;i++) {
            ell=mem[i];a=_cdid_bsearch(dates,C.r_s[ell]);b=_cdid_bsearch(dates,C.r_t[ell])
            from[i]=a;to[i]=b
            if(vmap[a]!=0 & vmap[a]!=D.l_from[ell])
                _cdid_cipw_fail("cdid_ipw_class_vertex_error","One class date maps to different qualified vertices.",459)
            if(vmap[b]!=0 & vmap[b]!=D.l_to[ell])
                _cdid_cipw_fail("cdid_ipw_class_vertex_error","One class date maps to different qualified vertices.",459)
            vmap[a]=D.l_from[ell];vmap[b]=D.l_to[ell]
            if(i>1) {
                if(from[i]==from[i-1] & to[i]==to[i-1])
                    _cdid_cipw_fail("cdid_ipw_class_duplicate_error","Duplicate date pair within a structural class.",459)
            }
        }
        if(rows(uniqrows(sort(vmap,1)))!=nv)
            _cdid_cipw_fail("cdid_ipw_class_vertex_error","Class date-to-qualified-vertex map is not injective.",459)
        F=_cdid_bfs_forest(nv,from,to,J(0,1,0),(1::nn))
        idx=_cdid_cipw_which(F.pe:>0);tree=sort(F.pe[idx],1)
        keep[mem[tree]]=J(rows(tree),1,1)
        for(i=1;i<=nv;i++) {
            at=off+i
            if(F.par[i]>0)par[at]=off+F.par[i]
            if(F.pe[i]>0)edge[at]=mem[F.pe[i]]
            sign[at]=F.ps[i]
            if(F.par[i]==0)depth[at]=0
        }
        stack=J(nv,1,0)
        for(i=1;i<=nv;i++) {
            u=off+i;top=0
            while(depth[u]<0) {top++;stack[top]=u;u=par[u];}
            while(top>0) {v=stack[top];depth[v]=depth[par[v]]+1;top--;}
        }
        rf[mem]=off:+from;rt[mem]=off:+to;off=off+nv
    }
    B.retained=_cdid_cipw_which(keep);B.dropped=_cdid_cipw_which(!keep);B.L=rows(B.retained)
    full_to_ret=J(M,1,0);full_to_ret[B.retained]=(1::B.L)
    npath=0
    for(h=1;h<=rows(B.dropped);h++) {
        ell=B.dropped[h];u=rf[ell];v=rt[ell];nn=0
        while(u!=v) {
            if(depth[u]>=depth[v]) {u=par[u];nn++;}
            else {v=par[v];nn++;}
        }
        B.path_count[ell]=nn;B.path_start[ell]=npath+1;npath=npath+nn
    }
    path_ws=_cdid_cipw_cells(7,npath,opts.max_matrix_entries,ws)
    path_ws=_cdid_cipw_cells(B.L,B.L,opts.max_matrix_entries,path_ws)
    B.e_row=J(npath,1,0);B.e_col=J(npath,1,0);B.e_value=J(npath,1,0)
    at=0
    for(h=1;h<=rows(B.dropped);h++) {
        ell=B.dropped[h];u=rf[ell];v=rt[ell]
        while(u!=v) {
            at++;B.e_row[at]=ell
            if(depth[u]>=depth[v]) {
                B.e_col[at]=full_to_ret[edge[u]];B.e_value[at]=-sign[u];u=par[u]
            }
            else {
                B.e_col[at]=full_to_ret[edge[v]];B.e_value[at]=sign[v];v=par[v]
            }
        }
        a=B.path_start[ell];b=a+B.path_count[ell]-1
        ord=order(B.e_col[|a\b|],1)
        B.e_col[|a\b|]=B.e_col[a:-1:+ord]
        B.e_value[|a\b|]=B.e_value[a:-1:+ord]
    }
    /* Integer signed-path algebra: G=I+E_drop'E_drop, never G/M. */
    B.G=I(B.L)
    for(h=1;h<=rows(B.dropped);h++) {
        ell=B.dropped[h];a=B.path_start[ell];b=a+B.path_count[ell]-1
        for(i=a;i<=b;i++)for(j=a;j<=b;j++)
            B.G[B.e_col[i],B.e_col[j]]=B.G[B.e_col[i],B.e_col[j]]+B.e_value[i]*B.e_value[j]
    }
    B.n_paths=npath;B.workspace_entries=path_ws-entries
    B.fit_order=B.retained[order((cl[B.retained],B.retained),(1,2))]
    return(B)
}

/* Evaluate one full-frame chunk [first_full..last_full] across all M links
   from compact per-link parameters/loadings without refitting any logit. */
real matrix _cdid_cipw_eval_chunk(
    struct cdid_population scalar P,
    struct cdid_cat scalar C,
    struct cdid_chained_ipw_postings scalar post,
    real colvector a_coh_idx,
    real matrix op_coh_map,
    real matrix X_A,
    real colvector w_A,
    struct cdid_chained_ipw_link_meta colvector links,
    real scalar first_full,
    real scalar last_full)
{
    real scalar b, M, a_lo, a_hi, b_A, ell, op, scale
    real scalar nhit, nrk, p_cov, nC, has_fast_lv, s_i, t_i
    real colvector pair_a, pair_dy, hit_j, risk_pos, risk_a, risk_j, risk_dy
    real colvector eta_hit, z_load_hit, a_vec, c_vec, w_vec, row_vec, j_vec, in_m
    real colvector row_rk, w_rk, fitw_vec, val_vec, t_ix, c_ix, p_rk, o_rk
    real matrix Phi, X_hit, lv_hit, B_beta, B_load, C_s0, Eta_all, Zload_all, O_all, P_all
    real matrix obs_chunk, Y_chunk, J_chunk

    b = last_full - first_full + 1
    M = rows(links)
    Phi = J(b, M, 0)
    if (b < 1 | M < 1 | P.n_A < 1) return(Phi)

    a_lo = _cdid_cipw_lb(P.A_to_F, 1, P.n_A, first_full)
    a_hi = _cdid_cipw_ub(P.A_to_F, 1, P.n_A, last_full)
    if (a_lo > a_hi) return(Phi)
    b_A = a_hi - a_lo + 1
    scale = P.N / P.n_A
    a_vec = (a_lo::a_hi)
    c_vec = a_coh_idx[a_vec]
    w_vec = w_A[a_vec]
    row_vec = P.A_to_F[a_vec] :- first_full :+ 1

    if (rows(post.obs_A) == post.n_A & b_A * M <= 2000000) {
        p_cov = cols(X_A)
        nC = rows(op_coh_map)
        B_beta = J(p_cov, M, 0)
        B_load = J(p_cov, M, 0)
        C_s0 = J(nC, M, 0)
        for (ell = 1; ell <= M; ell++) {
            B_beta[., ell] = links[ell].beta
            B_load[., ell] = links[ell].loading_beta
            op = links[ell].op_index
            in_m = _cdid_cipw_which(op_coh_map[., op] :> 0)
            if (rows(in_m) > 0) C_s0[in_m, ell] = links[ell].rate_coef_s0[op_coh_map[in_m, op]]
        }
        Phi[row_vec, .] = (scale :* w_vec) :* C_s0[c_vec, .]
        X_hit = X_A[a_vec, .]
        Eta_all = X_hit * B_beta
        Zload_all = X_hit * B_load
        has_fast_lv = (!missing(Eta_all) & max(Eta_all) <= 13.815 & min(Eta_all) >= -700)
        if (has_fast_lv) {
            O_all = exp(Eta_all)
            P_all = O_all :/ (1 :+ O_all)
        }
        obs_chunk = post.obs_A[a_vec, .]
        Y_chunk = post.Y_A[a_vec, .]
        J_chunk = op_coh_map[c_vec, .]
        for (ell = 1; ell <= M; ell++) {
            op = links[ell].op_index
            s_i = links[ell].s_index
            t_i = links[ell].t_index
            if (post.date_count[s_i] == 0 | post.date_count[t_i] == 0) continue
            risk_pos = vec(selectindex(obs_chunk[., s_i] :& obs_chunk[., t_i] :& (J_chunk[., op] :> 0)))
            nrk = rows(risk_pos)
            if (nrk == 0) continue
            risk_j = J_chunk[risk_pos, op]
            risk_dy = Y_chunk[risk_pos, t_i] :- Y_chunk[risk_pos, s_i]
            if (missing(risk_dy))
                _cdid_cipw_fail("cdid_ipw_outcome_error", "Pair outcome difference overflow.", 504)
            row_rk = row_vec[risk_pos]
            w_rk = w_vec[risk_pos]
            fitw_vec = w_rk :/ links[ell].rates_c[risk_j]
            z_load_hit = Zload_all[risk_pos, ell]
            if (has_fast_lv) {
                p_rk = P_all[risk_pos, ell]
                o_rk = O_all[risk_pos, ell]
            }
            else {
                lv_hit = _cdid_ipw_logit_values(Eta_all[risk_pos, ell])
                p_rk = lv_hit[., 1]
                o_rk = lv_hit[., 2]
            }
            val_vec = scale :* (w_rk :* links[ell].rate_coef_s1[risk_j]) :+
                      scale :* (fitw_vec :* ((risk_j :== 1) :- p_rk) :* z_load_hit)
            t_ix = _cdid_cipw_which(risk_j :== 1)
            if (rows(t_ix) > 0) {
                val_vec[t_ix] = val_vec[t_ix] :+
                    scale :* (fitw_vec[t_ix] :/ links[ell].mean_a :* (risk_dy[t_ix] :- links[ell].treated_mean))
            }
            c_ix = _cdid_cipw_which(risk_j :> 1)
            if (rows(c_ix) > 0) {
                val_vec[c_ix] = val_vec[c_ix] :+
                    scale :* (-fitw_vec[c_ix] :* o_rk[c_ix] :/ links[ell].mean_b :* (risk_dy[c_ix] :- links[ell].control_mean))
            }
            Phi[row_rk, ell] = val_vec
        }
        if (missing(Phi))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite original-unit IF contribution during chunk replay.", 504)
        return(Phi)
    }

    for (ell = 1; ell <= M; ell++) {
        op = links[ell].op_index
        j_vec = op_coh_map[c_vec, op]
        in_m = _cdid_cipw_which(j_vec :> 0)
        if (rows(in_m) > 0) {
            Phi[row_vec[in_m], ell] = scale :* (w_vec[in_m] :* links[ell].rate_coef_s0[j_vec[in_m]])
        }
        _cdid_cipw_query_pair(post, links[ell].s_index, links[ell].t_index,
            a_lo, a_hi, pair_a = ., pair_dy = .)
        nhit = rows(pair_a)
        if (nhit == 0) continue

        hit_j = op_coh_map[a_coh_idx[pair_a], op]
        risk_pos = _cdid_cipw_which(hit_j :> 0)
        nrk = rows(risk_pos)
        if (nrk == 0) continue

        risk_a = pair_a[risk_pos]
        risk_j = hit_j[risk_pos]
        risk_dy = pair_dy[risk_pos]
        X_hit = X_A[risk_a, .]
        eta_hit = X_hit * links[ell].beta
        lv_hit = _cdid_ipw_logit_values(eta_hit)
        z_load_hit = X_hit * links[ell].loading_beta

        row_rk = P.A_to_F[risk_a] :- first_full :+ 1
        w_rk = w_A[risk_a]
        fitw_vec = w_rk :/ links[ell].rates_c[risk_j]
        val_vec = scale :* (w_rk :* links[ell].rate_coef_s1[risk_j]) :+
                  scale :* (fitw_vec :* ((risk_j :== 1) :- lv_hit[., 1]) :* z_load_hit)
        t_ix = _cdid_cipw_which(risk_j :== 1)
        if (rows(t_ix) > 0) {
            val_vec[t_ix] = val_vec[t_ix] :+
                scale :* (fitw_vec[t_ix] :/ links[ell].mean_a :* (risk_dy[t_ix] :- links[ell].treated_mean))
        }
        c_ix = _cdid_cipw_which(risk_j :> 1)
        if (rows(c_ix) > 0) {
            val_vec[c_ix] = val_vec[c_ix] :+
                scale :* (-fitw_vec[c_ix] :* lv_hit[c_ix, 2] :/ links[ell].mean_b :* (risk_dy[c_ix] :- links[ell].control_mean))
        }
        Phi[row_rk, ell] = val_vec
    }
    if (missing(Phi))
        _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite original-unit IF contribution during chunk replay.", 504)
    return(Phi)
}

real matrix _cdid_cipw_eval_a_chunk(
    struct cdid_population scalar P,
    struct cdid_cat scalar C,
    struct cdid_chained_ipw_postings scalar post,
    real colvector a_coh_idx,
    real matrix op_coh_map,
    real matrix X_A,
    real colvector w_A,
    struct cdid_chained_ipw_link_meta colvector links,
    real scalar first_full,
    real scalar last_full)
{
    real scalar b, M, a_lo, a_hi, b_A, ell, op, scale
    real scalar nhit, nrk, p_cov, nC, has_fast_lv, s_i, t_i
    real colvector pair_a, pair_dy, hit_j, risk_pos, risk_a, risk_j, risk_dy
    real colvector eta_hit, z_load_hit, a_vec, c_vec, w_vec, row_vec, j_vec, in_m
    real colvector row_rk, w_rk, fitw_vec, val_vec, t_ix, c_ix, p_rk, o_rk
    real matrix Phi, X_hit, lv_hit, B_beta, B_load, C_s0, Eta_all, Zload_all, O_all, P_all
    real matrix obs_chunk, Y_chunk, J_chunk

    b = last_full - first_full + 1
    M = rows(links)
    Phi = J(b, M, 0)
    if (b < 1 | M < 1 | P.n_A < 1) return(Phi)

    a_lo = first_full
    a_hi = last_full
    if (a_lo > a_hi) return(Phi)
    b_A = a_hi - a_lo + 1
    scale = 1
    a_vec = (a_lo::a_hi)
    c_vec = a_coh_idx[a_vec]
    w_vec = w_A[a_vec]
    row_vec = a_vec :- first_full :+ 1

    if (rows(post.obs_A) == post.n_A & b_A * M <= 2000000) {
        p_cov = cols(X_A)
        nC = rows(op_coh_map)
        B_beta = J(p_cov, M, 0)
        B_load = J(p_cov, M, 0)
        C_s0 = J(nC, M, 0)
        for (ell = 1; ell <= M; ell++) {
            B_beta[., ell] = links[ell].beta
            B_load[., ell] = links[ell].loading_beta
            op = links[ell].op_index
            in_m = _cdid_cipw_which(op_coh_map[., op] :> 0)
            if (rows(in_m) > 0) C_s0[in_m, ell] = links[ell].rate_coef_s0[op_coh_map[in_m, op]]
        }
        Phi = w_vec :* C_s0[c_vec, .]
        X_hit = X_A[a_vec, .]
        Eta_all = X_hit * B_beta
        Zload_all = X_hit * B_load
        has_fast_lv = (!missing(Eta_all) & max(Eta_all) <= 13.815 & min(Eta_all) >= -700)
        if (has_fast_lv) {
            O_all = exp(Eta_all)
            P_all = O_all :/ (1 :+ O_all)
        }
        obs_chunk = post.obs_A[a_vec, .]
        Y_chunk = post.Y_A[a_vec, .]
        J_chunk = op_coh_map[c_vec, .]
        for (ell = 1; ell <= M; ell++) {
            op = links[ell].op_index
            s_i = links[ell].s_index
            t_i = links[ell].t_index
            if (post.date_count[s_i] == 0 | post.date_count[t_i] == 0) continue
            risk_pos = vec(selectindex(obs_chunk[., s_i] :& obs_chunk[., t_i] :& (J_chunk[., op] :> 0)))
            nrk = rows(risk_pos)
            if (nrk == 0) continue
            risk_j = J_chunk[risk_pos, op]
            risk_dy = Y_chunk[risk_pos, t_i] :- Y_chunk[risk_pos, s_i]
            if (missing(risk_dy))
                _cdid_cipw_fail("cdid_ipw_outcome_error", "Pair outcome difference overflow.", 504)
            w_rk = w_vec[risk_pos]
            fitw_vec = w_rk :/ links[ell].rates_c[risk_j]
            z_load_hit = Zload_all[risk_pos, ell]
            if (has_fast_lv) {
                p_rk = P_all[risk_pos, ell]
                o_rk = O_all[risk_pos, ell]
            }
            else {
                lv_hit = _cdid_ipw_logit_values(Eta_all[risk_pos, ell])
                p_rk = lv_hit[., 1]
                o_rk = lv_hit[., 2]
            }
            val_vec = (w_rk :* links[ell].rate_coef_s1[risk_j]) :+
                      (fitw_vec :* ((risk_j :== 1) :- p_rk) :* z_load_hit)
            t_ix = _cdid_cipw_which(risk_j :== 1)
            if (rows(t_ix) > 0) {
                val_vec[t_ix] = val_vec[t_ix] :+
                    (fitw_vec[t_ix] :/ links[ell].mean_a :* (risk_dy[t_ix] :- links[ell].treated_mean))
            }
            c_ix = _cdid_cipw_which(risk_j :> 1)
            if (rows(c_ix) > 0) {
                val_vec[c_ix] = val_vec[c_ix] :+
                    (-fitw_vec[c_ix] :* o_rk[c_ix] :/ links[ell].mean_b :* (risk_dy[c_ix] :- links[ell].control_mean))
            }
            Phi[risk_pos, ell] = val_vec
        }
        if (missing(Phi))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite original-unit IF contribution during chunk replay.", 504)
        return(Phi)
    }

    for (ell = 1; ell <= M; ell++) {
        op = links[ell].op_index
        j_vec = op_coh_map[c_vec, op]
        in_m = _cdid_cipw_which(j_vec :> 0)
        if (rows(in_m) > 0) {
            Phi[row_vec[in_m], ell] = w_vec[in_m] :* links[ell].rate_coef_s0[j_vec[in_m]]
        }
        _cdid_cipw_query_pair(post, links[ell].s_index, links[ell].t_index,
            a_lo, a_hi, pair_a = ., pair_dy = .)
        nhit = rows(pair_a)
        if (nhit == 0) continue

        hit_j = op_coh_map[a_coh_idx[pair_a], op]
        risk_pos = _cdid_cipw_which(hit_j :> 0)
        nrk = rows(risk_pos)
        if (nrk == 0) continue

        risk_a = pair_a[risk_pos]
        risk_j = hit_j[risk_pos]
        risk_dy = pair_dy[risk_pos]
        X_hit = X_A[risk_a, .]
        eta_hit = X_hit * links[ell].beta
        lv_hit = _cdid_ipw_logit_values(eta_hit)
        z_load_hit = X_hit * links[ell].loading_beta

        row_rk = risk_a :- first_full :+ 1
        w_rk = w_A[risk_a]
        fitw_vec = w_rk :/ links[ell].rates_c[risk_j]
        val_vec = (w_rk :* links[ell].rate_coef_s1[risk_j]) :+
                  (fitw_vec :* ((risk_j :== 1) :- lv_hit[., 1]) :* z_load_hit)
        t_ix = _cdid_cipw_which(risk_j :== 1)
        if (rows(t_ix) > 0) {
            val_vec[t_ix] = val_vec[t_ix] :+
                (fitw_vec[t_ix] :/ links[ell].mean_a :* (risk_dy[t_ix] :- links[ell].treated_mean))
        }
        c_ix = _cdid_cipw_which(risk_j :> 1)
        if (rows(c_ix) > 0) {
            val_vec[c_ix] = val_vec[c_ix] :+
                (-fitw_vec[c_ix] :* lv_hit[c_ix, 2] :/ links[ell].mean_b :* (risk_dy[c_ix] :- links[ell].control_mean))
        }
        Phi[row_rk, ell] = val_vec
    }
    if (missing(Phi))
        _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite original-unit IF contribution during chunk replay.", 504)
    return(Phi)
}

/* Baseline source Gram consumes A-coordinate raw IF; reporting stays full F. */
real matrix _cdid_cipw_gram_provider(pointer() scalar state, real scalar lo, real scalar hi)
{
    struct cdid_cipw_replay_context scalar ctx
    ctx=*state
    return(_cdid_cipw_eval_a_chunk(*(ctx.population),*(ctx.catalog),*(ctx.postings),
        *(ctx.cohorts),*(ctx.operators),*(ctx.X),*(ctx.w),*(ctx.links),lo,hi))
}

void _cdid_cipw_init_empty_result(
    struct cdid_chained_ipw_result scalar res,
    struct cdid_population scalar P,
    struct cdid_cat scalar C,
    struct cdid_design scalar D,
    real scalar p_cov,
    real colvector w_F,
    struct cdid_chained_ipw_opts scalar opts)
{
    res.status = "ok"
    res.gmm_status="not_run"
    res.moment_covariance_status=(opts.materialize_moment_covariance ? "requested" : "not_requested")
    res.global_reason = ""
    res.control_group = C.control_group
    res.pt_assumption = (D.pt_all ? "all" : opts.pt_assumption)
    res.weighting = ((opts.gmm_weighting == "two_step" | opts.gmm_weighting == "within") ? "optimal" : opts.gmm_weighting)
    res.pt_all = (res.pt_assumption == "all")
    res.pt_start = (res.pt_all ? (D.pt_all ? D.pt_start : opts.pt_start) : .)
    res.pt_start_index = (res.pt_all ? _cdid_bsearch(C.calendar, res.pt_start) : 0)
    res.anticipation = C.anticipation

    res.N_F = P.N
    res.n_panel = P.n_panel
    res.n_A = P.n_A
    res.n_B = P.n_B
    res.n_singleton = P.n_B - P.n_A
    res.n_unobserved = P.n_panel - P.n_B
    res.n_frame_only = P.n_frame_only
    res.n_dates = C.n_dates
    res.n_cohorts = C.n_cohorts
    res.p_covariates = p_cov

    res.row_basis.M=rows(C.r_tp);res.row_basis.L=0;res.row_basis.n_classes=0
    res.row_basis.retained=J(0,1,0);res.row_basis.dropped=J(0,1,0)
    res.row_basis.G=J(0,0,0);res.row_basis.n_paths=0;res.row_basis.workspace_entries=0
    res.M_links = rows(C.r_tp)
    res.n_operators = rows(C.op_g)
    res.n_pools = C.n_pools
    res.n_vertices = D.n_vertices
    res.n_params_unreduced = D.n_params
    res.n_params_retained = D.n_rank
    res.n_nullity = D.nullity
    res.n_components = D.n_components
    res.n_zero_loading_links = rows(D.zero_loading_links)

    res.n_targets = 0
    res.n_estimable_targets = 0
    res.n_unavailable_targets = 0
    res.n_nuisance_fits = 0
    res.n_nuisance_reused = 0
    res.chunk_size = opts.chunk_size
    res.effective_chunk_size = min((opts.chunk_size, P.N))
    res.n_batches = ceil(P.N / res.effective_chunk_size)
    res.workspace_entries_bound = 0
    res.workspace_entries_limit = opts.max_matrix_entries
    res.work_ops = 0
    res.work_cap = opts.work_cap
    res.no_extra_calendar_window = (res.pt_all & D.built ? cdid_design_diag_value(D, "no_extra_calendar_window") : .)
    res.no_extra_represented_restriction = (res.pt_all & D.built ? cdid_design_diag_value(D, "no_extra_represented_restriction") : .)
    res.no_extra_restrictions = (res.pt_all & D.built ? cdid_design_diag_value(D, "no_extra_restrictions") : .)

    res.calendar = C.calendar
    res.cohorts = C.cohorts
    res.A_to_F = P.A_to_F
    res.F_to_A = P.F_to_A
    res.unit_cohort = P.cohort
    res.unit_cohort_known = P.cohort_known
    res.unit_weight = w_F

    res.moments = J(0, 1, .)
    res.moment_if_mean = J(0, 1, .)
    res.moment_fitted = J(0, 1, .)
    res.moment_residual = J(0, 1, .)
    res.moment_table = J(0, 14, .)
    res.links = cdid_chained_ipw_link_meta(0, 1)
    res.Vr_estimator = J(0, 0, .)
    res.Vr_analytical = J(0, 0, .)

    res.kept_param_indices = J(0, 1, 0)
    res.dropped_gauge_vertices = J(0, 1, 0)
    res.dropped_gauge_params = J(0, 1, 0)
    res.A_retained = J(0, 0, .)
    res.H_estimable = J(0, 0, .)
    res.criterion = .

    res.target_g = J(0, 1, .)
    res.target_t = J(0, 1, .)
    res.target_g_index = J(0, 1, 0)
    res.target_t_index = J(0, 1, 0)
    res.target_cat_supported = J(0, 1, 0)
    res.target_cat_ntouch = J(0, 1, 0)
    res.target_cat_reason = J(0, 1, 0)
    res.target_in_design = J(0, 1, 0)
    res.target_vertex = J(0, 1, 0)
    res.target_param = J(0, 1, 0)
    res.target_design_status = J(0, 1, 0)
    res.target_available = J(0, 1, 0)
    res.target_cat_reason_name = J(0, 1, "")
    res.target_design_status_name = J(0, 1, "")
    res.target_unavailable_reason = J(0, 1, "")
    res.estimable_target_indices = J(0, 1, 0)
    res.unavailable_target_indices = J(0, 1, 0)
    res.att = J(0, 1, .)
    res.se = J(0, 1, .)
    res.target_ledger = J(0, 14, .)

    res.target_if = J(P.N, 0, .)
    res.target_if_mean = J(1, 0, .)
    res.V_estimator = J(0, 0, .)
    res.V_analytical = J(0, 0, .)
    res.V_estimator_estimable = J(0, 0, .)
    res.V_analytical_estimable = J(0, 0, .)
}

void _cdid_cipw_build_ledger(
    struct cdid_chained_ipw_result scalar res,
    struct cdid_cat scalar C,
    struct cdid_design scalar D,
    real matrix req, | real scalar compute_target_if)
{
    real scalar K, i, gi, ti, v, st, M, keep_if
    real colvector gq, tq, hit
    real matrix sq, cl
    keep_if=(args()==5 ? compute_target_if : 1)
    if(keep_if!=0 & keep_if!=1)_cdid_cipw_fail("cdid_chained_ipw_option_error","compute_target_if must be zero or one.",498)
    K = rows(req)
    M = rows(C.r_tp)
    if (K > 0 & cols(req) != 2)
        _cdid_cipw_fail("cdid_chained_ipw_target_error", "requested_targets must have exactly 2 columns (g, t).", 498)
    gq = J(K, 1, 0)
    tq = J(K, 1, 0)
    for (i = 1; i <= K; i++) {
        if (req[i, 1] >= . | req[i, 2] >= .)
            _cdid_cipw_fail("cdid_chained_ipw_target_error", "Requested targets require finite numeric g and t.", 498)
        gi = _cdid_cohort_index(C.cohorts, C.nt_index, req[i, 1])
        if (gi == 0 | !_cdid_cipw_any(C.tg :== gi))
            _cdid_cipw_fail("cdid_chained_ipw_target_error", "Requested target cohort must belong to catalogue target cohorts.", 459)
        ti = _cdid_bsearch(C.calendar, req[i, 2])
        if (ti == 0)
            _cdid_cipw_fail("cdid_chained_ipw_target_error", "Requested target date must belong to the full calendar.", 459)
        if (req[i, 2] < req[i, 1])
            _cdid_cipw_fail("cdid_chained_ipw_target_error", "Requested post targets require t >= g.", 459)
        gq[i] = gi
        tq[i] = ti
    }
    if (K > 1) {
        sq = sort((gq, tq), (1, 2))
        if (_cdid_cipw_any((sq[|2, 1 \ K, 1|] :== sq[|1, 1 \ K - 1, 1|]) :&
            (sq[|2, 2 \ K, 2|] :== sq[|1, 2 \ K - 1, 2|])))
            _cdid_cipw_fail("cdid_chained_ipw_duplicate_target_error", "Duplicate requested (g, t) target.", 459)
    }

    cl = _cdid_cat_classify(C, gq, tq)
    res.n_targets = K
    if (K == 0) return

    res.target_g = req[., 1]
    res.target_t = req[., 2]
    res.target_g_index = gq
    res.target_t_index = tq
    res.target_cat_supported = cl[., 1]
    res.target_cat_ntouch = cl[., 2]
    res.target_cat_reason = cl[., 3]
    res.target_in_design = J(K, 1, 0)
    res.target_vertex = J(K, 1, 0)
    res.target_param = J(K, 1, 0)
    res.target_design_status = J(K, 1, 0)
    res.target_available = J(K, 1, 0)
    res.target_cat_reason_name = J(K, 1, "")
    res.target_design_status_name = J(K, 1, "")
    res.target_unavailable_reason = J(K, 1, "")

    for (i = 1; i <= K; i++) {
        res.target_cat_reason_name[i] = cdid_request_reason_name(res.target_cat_reason[i])
        if (D.built) {
            hit = _cdid_cipw_which((D.v_g :== gq[i]) :& (D.v_date :== tq[i]) :& (D.v_kind :== 6))
            if (rows(hit) == 1) {
                v = hit[1]
                st = (D.v_conn[v] ? 1 : (D.v_degree[v] == 0 ? 2 : 3))
                res.target_in_design[i] = 1
                res.target_vertex[i] = v
                res.target_param[i] = D.v_param[v]
                res.target_design_status[i] = st
                res.target_available[i] = (st == 1)
            }
        }
        res.target_design_status_name[i] = cdid_post_status_name(res.target_design_status[i])
        if (!res.target_available[i]) {
            if (M == 0) {
                res.target_unavailable_reason[i] = "empty_catalog:" + res.target_cat_reason_name[i]
            }
            else if (!res.target_cat_supported[i]) {
                res.target_unavailable_reason[i] = res.target_cat_reason_name[i]
            }
            else {
                res.target_unavailable_reason[i] = res.target_design_status_name[i]
            }
        }
    }

    res.estimable_target_indices = _cdid_cipw_which(res.target_available)
    res.unavailable_target_indices = _cdid_cipw_which(!res.target_available)
    res.n_estimable_targets = rows(res.estimable_target_indices)
    res.n_unavailable_targets = rows(res.unavailable_target_indices)
    res.att = J(K, 1, .)
    res.se = J(K, 1, .)
    res.target_ledger = (res.target_g, res.target_t, res.target_g_index, res.target_t_index,
        res.target_cat_supported, res.target_cat_ntouch, res.target_cat_reason,
        res.target_in_design, res.target_vertex, res.target_param,
        res.target_design_status, res.target_available, res.att, res.se)
    res.target_if = J((keep_if ? res.N_F : 0), K, .)
    res.target_if_mean = J(1, K, .)
    res.V_estimator = (keep_if ? J(K,K,.) : J(0,0,.))
    res.V_analytical = (keep_if ? J(K,K,.) : J(0,0,.))
}

void _cdid_cipw_link_row(struct cdid_chained_ipw_link_meta scalar meta,
    struct cdid_cat scalar C,struct cdid_design scalar D,real scalar ell)
{
    real scalar op,pid
    op=C.r_op[ell];pid=C.op_pool[op]
    meta.row_index = ell
    meta.tp_index = C.r_tp[ell]
    meta.op_index = op
    meta.pool_index = pid
    meta.g_index = C.r_g[ell]
    meta.s_index = C.r_s[ell]
    meta.t_index = C.r_t[ell]
    meta.ref_index = C.r_ref[ell]
    meta.g_label = C.cohorts[C.r_g[ell]]
    meta.s_label = C.calendar[C.r_s[ell]]
    meta.t_label = C.calendar[C.r_t[ell]]
    meta.ref_label = (C.r_ref[ell] > 0 ? C.calendar[C.r_ref[ell]] : .)
    meta.s_role = C.r_srole[ell]
    meta.t_role = C.r_trole[ell]
    meta.s_qualified = D.l_sq[ell]
    meta.t_qualified = D.l_tq[ell]
    meta.from_vertex = D.l_from[ell]
    meta.to_vertex = D.l_to[ell]
    meta.start_param = D.l_start[ell]
    meta.end_param = D.l_end[ell]
    meta.link_type = D.l_type[ell]
    meta.component_id = D.l_comp[ell]
    meta.zero_loading = (D.l_start[ell] == 0 & D.l_end[ell] == 0)
}

/* Core baseline static-input IPW assembly with full catalogue and source weighting. */
struct cdid_chained_ipw_result scalar cdid_chained_ipw_fit_opts(
    struct cdid_population scalar P,
    struct cdid_cat scalar C,
    struct cdid_design scalar D_in,
    real matrix X_F,
    real colvector w_F,
    real matrix requested_targets,
    struct cdid_chained_ipw_opts scalar opts,
    | real scalar compute_target_if,
      real scalar skip_target_cov,
      pointer(pointer(real colvector) colvector) scalar ps_co_out,
      pointer(struct cdid_ipw_nuisance colvector) scalar ps_nu_out)
{
    struct cdid_chained_ipw_result scalar res
    struct cdid_design scalar D
    struct cdid_chained_ipw_postings scalar post
    struct cdid_ipw_frame colvector op_frames
    struct cdid_ipw_nuisance scalar nu, nu_norm
    struct cdid_ipw_nuisance colvector pool_nu
    pointer(real colvector) colvector pool_co
    struct cdid_ipw_link scalar lk
    struct cdid_chained_ipw_cache_entry colvector cache
    struct cdid_cipw_replay_context scalar replay_context
    struct cdid_chained_ipw_link_meta colvector replay_links

    real scalar p_cov, u, a, ci, M, L, K, nop, op, pid, p_sz, j, ell, total_ws, total_work
    real scalar P_unred, P_ret, J_est, eff_chunk, b_first, b_last, ps, pe, rs, re
    real scalar it, h, ep, v, last_class, this_class, pair_all_n, keep_if, no_tcov, want_ps_pool, n_pool
    real scalar n_cache_slots, hit_slot, evict_slot, min_stamp, clock, c_i, v_root, p_root
    real colvector pop_cohorts, u_ci, a_coh_idx, a_coh_counts, a_known, w_A, pool_k, coh_idx_op, op_frame_built, ac_pos
    real colvector pair_a, pair_dy, risk_keep, keep_param_mask, full_to_ret, unanch_roots
    real rowvector sum_phi, mean_phi, T_mean
    real matrix primitive_A, primitive_m, metric_G, ac_spans
    real matrix req, X_A, op_coh_map, Phi_chunk, Phi_centered, cov_sum
    real matrix T_if, T_centered, T_cross, V_est_sub, V_asymp_sub, coef_cross, coef_if

    cdid_chained_ipw_clear_error()
    _cdid_cipw_check_opts(opts)
    keep_if=(args()>=8 & compute_target_if<. ? compute_target_if : 1)
    if(keep_if!=0 & keep_if!=1)_cdid_cipw_fail("cdid_chained_ipw_option_error","compute_target_if must be zero or one.",498)
    no_tcov=(args()>=9 & skip_target_cov==1)
    want_ps_pool=(args()>=11 & ps_co_out!=NULL & ps_nu_out!=NULL)
    if(want_ps_pool) {
        *ps_co_out=J(0,1,NULL)
        *ps_nu_out=cdid_ipw_nuisance(0,1)
    }

    if (P.N < 1 | rows(P.calendar) < 2)
        _cdid_cipw_fail("cdid_chained_ipw_population_error", "Prepared population requires N >= 1 and at least 2 calendar dates.", 498)
    if (rows(C.calendar) != rows(P.calendar) | _cdid_cipw_any(C.calendar :!= P.calendar))
        _cdid_cipw_fail("cdid_chained_ipw_dictionary_error", "Catalogue calendar disagrees with population calendar.", 498)

    pop_cohorts = cdid_cipw_cohorts(P)
    if (rows(C.cohorts) != rows(pop_cohorts) | _cdid_cipw_any(C.cohorts :!= pop_cohorts))
        _cdid_cipw_fail("cdid_chained_ipw_dictionary_error", "Catalogue cohorts dictionary disagrees with known population cohorts.", 498)
    if (any((P.cohort_known :!= 0) :& (P.cohort_known :!= 1)))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "cohort_known must be binary 0/1.", 498)
    if (any((!P.cohort_known) :& (P.A :!= 0)))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Unknown cohort is forbidden on pair-eligible A units.", 459)
    u_ci = _cdid_cohort_index_vec(C.cohorts, C.nt_index, P.cohort, P.cohort_known)
    if (any((P.cohort_known :!= 0) :& (u_ci :== 0)))
        _cdid_cipw_fail("cdid_chained_ipw_cohort_error", "Known unit cohort is absent from catalogue cohorts dictionary.", 459)

    p_cov = cols(X_F)
    if (rows(X_F) != P.N | p_cov < 1 | missing(X_F))
        _cdid_cipw_fail("cdid_ipw_covariate_error", "Static full-F design X_F must be a finite P.N x p matrix with p >= 1.", 498)
    if (rows(w_F) != P.N | cols(w_F) != 1 | missing(w_F) | _cdid_cipw_any(w_F :<= 0))
        _cdid_cipw_fail("cdid_ipw_weight_error", "Static full-F weights w_F must be a finite strictly positive P.N x 1 vector.", 498)

    req = requested_targets
    if (rows(req) == 0 & rows(C.rq_g) > 0) {
        req = (C.cohorts[C.rq_g], C.calendar[C.rq_t])
    }

    M = rows(C.r_tp)
    D = D_in
    if (M == 0) {
        if (D.built != 0)
            _cdid_cipw_fail("cdid_chained_ipw_design_error", "Catalogue has zero links but supplied design has built != 0.", 498)
        if (opts.pt_assumption == "all") {
            if (_cdid_bsearch(C.calendar, opts.pt_start) == 0)
                _cdid_cipw_fail("cdid_design_pt_start_error", "pt_start must belong to the supplied full calendar.", 459)
        }
    }
    else {
        if (!D.built) {
            D = cdid_chained_ipw_build_design(C, req, opts.pt_assumption, opts.pt_start)
        }
        if (D.n_links != M | rows(D.calendar) != rows(C.calendar) | _cdid_cipw_any(D.calendar :!= C.calendar) |
            rows(D.cohorts) != rows(C.cohorts) | _cdid_cipw_any(D.cohorts :!= C.cohorts) |
            D.anticipation != C.anticipation | _cdid_cipw_any(D.l_g :!= C.r_g) |
            _cdid_cipw_any(D.l_s :!= C.r_s) | _cdid_cipw_any(D.l_t :!= C.r_t) |
            _cdid_cipw_any(D.l_op :!= C.r_op))
            _cdid_cipw_fail("cdid_chained_ipw_design_error", "Supplied design D does not match catalogue C.", 498)
    }

    _cdid_cipw_init_empty_result(res, P, C, D, p_cov, w_F, opts)
    _cdid_cipw_build_ledger(res, C, D, req,keep_if)
    K = res.n_targets
    J_est = res.n_estimable_targets

    if (M == 0) {
        total_ws = 0
        if(keep_if) total_ws = _cdid_cipw_cells(P.N, K, opts.max_matrix_entries, total_ws)
        if(keep_if) total_ws = _cdid_cipw_cells(2 * K, K, opts.max_matrix_entries, total_ws)
        total_ws = _cdid_cipw_cells(14, K, opts.max_matrix_entries, total_ws)
        res.workspace_entries_bound = total_ws
        res.work_ops = _cdid_cipw_work_add(P.N + K + C.work, opts.work_cap, 0)
        res.gmm_status="not_run_empty_catalog"
        res.status = "empty_catalog"
        res.global_reason = "the comparison catalog retained zero links for the requested target cohorts and control group"
        return(res)
    }

    if (P.n_A < 1 | rows(P.A_to_F) != P.n_A)
        _cdid_cipw_fail("cdid_ipw_support_error", "At least one pair-eligible A-unit (n_obs >= 2) is required.", 459)

    nop = rows(C.op_g)
    P_unred = D.n_params
    P_ret = D.n_rank
    eff_chunk = res.effective_chunk_size
    n_cache_slots = min((opts.max_cached_nuisances, nop))

    // Checked workspace preflight before allocating estimator/GMM/replay arrays
    total_ws = 0
    total_ws = _cdid_cipw_cells(4, P.R + P.N + C.n_dates, opts.max_matrix_entries, total_ws)
    total_ws = _cdid_cipw_cells(C.n_cohorts, nop, opts.max_matrix_entries, total_ws)
    total_ws = _cdid_cipw_cells(M, 3 * p_cov + 12 * C.n_cohorts + 16, opts.max_matrix_entries, total_ws)
    total_ws = _cdid_cipw_cells(max((1, n_cache_slots)) * P.n_A, p_cov + 8, opts.max_matrix_entries, total_ws)
    if(keep_if) total_ws = _cdid_cipw_cells(eff_chunk, J_est, opts.max_matrix_entries, total_ws)
    if(opts.materialize_moment_covariance) {
        total_ws = _cdid_cipw_cells(4 * M, M, opts.max_matrix_entries, total_ws)
        total_ws = _cdid_cipw_cells(3 * eff_chunk, M, opts.max_matrix_entries, total_ws)
    }
    total_ws = _cdid_cipw_cells(M, max((1, P_ret)), opts.max_matrix_entries, total_ws)
    total_ws = _cdid_cipw_cells(J_est, max((1, P_ret)), opts.max_matrix_entries, total_ws)
    if(keep_if) total_ws = _cdid_cipw_cells(P.N, K, opts.max_matrix_entries, total_ws)
    if(keep_if) total_ws = _cdid_cipw_cells(2 * K, K, opts.max_matrix_entries, total_ws)
    total_ws = _cdid_cipw_cells(2 * J_est, J_est, opts.max_matrix_entries, total_ws)
    if (P_ret>0 & res.weighting=="identity") {
        total_ws=_cdid_gmm_identity_workspace(M,P_ret,J_est,opts.max_matrix_entries,total_ws)
        if(keep_if) {
            total_ws=_cdid_cipw_cells(2*P_ret,P_ret,opts.max_matrix_entries,total_ws)
            total_ws=_cdid_cipw_cells(eff_chunk,P_ret,opts.max_matrix_entries,total_ws)
        }
    }
    res.workspace_entries_bound = total_ws

    total_work = 0
    total_work = _cdid_cipw_work_add(C.work, opts.work_cap, total_work)
    total_work = _cdid_cipw_work_add(M * P.n_A * p_cov, opts.work_cap, total_work)
    total_work = _cdid_cipw_work_add(2 * P.N * M + P.N * J_est, opts.work_cap, total_work)

    // Map A-units to cohort dictionary indices
    a_coh_idx = u_ci[P.A_to_F]
    a_coh_counts = J(C.n_cohorts, 1, 0)
    ac_pos = select(a_coh_idx, a_coh_idx :> 0)
    if (rows(ac_pos) > 0) {
        ac_pos = sort(ac_pos, 1)
        ac_spans = _cdid_ps_cipw(ac_pos, 1)
        a_coh_counts[ac_pos[ac_spans[., 1]]] = ac_spans[., 2] :- ac_spans[., 1] :+ 1
    }
    a_known = P.cohort_known[P.A_to_F]
    w_A = w_F[P.A_to_F]
    X_A = X_F[P.A_to_F, .]
    if (sum(w_A) >= . | sum(w_A) <= 0)
        _cdid_cipw_fail("cdid_ipw_weight_error", "Total A-unit weight must be finite and positive.", 498)

    // Build compact operator cohort-to-group map (C.n_cohorts x n_operators)
    // Group 1 = target cohort g; groups 2..m = declared control cohorts in pool order; 0 = off-risk.
    op_coh_map = J(C.n_cohorts, nop, 0)
    for (op = 1; op <= nop; op++) {
        pid = C.op_pool[op]
        p_sz = C.pool_size[pid]
        pool_k = C.pool_entries[|C.pool_start[pid] \ C.pool_start[pid] + p_sz - 1|]
        op_coh_map[C.op_g[op], op] = 1
        for (j = 1; j <= p_sz; j++) {
            op_coh_map[pool_k[j], op] = j + 1
        }
    }

    // Build compact observed-outcome calendar postings over A-units
    post = _cdid_cipw_build_postings(P, C.calendar)

    res.row_basis=_cdid_cipw_row_basis(C,D,post,a_coh_idx,op_coh_map,opts,total_ws,total_work)
    L=res.row_basis.L
    total_ws=total_ws+res.row_basis.workspace_entries
    total_ws=_cdid_cipw_cells(eff_chunk,L,opts.max_matrix_entries,total_ws)
    if(P_ret>0 & res.weighting=="optimal") {
        total_ws=_cdid_cipw_cells(2*min((eff_chunk,P.n_A)),L,opts.max_matrix_entries,total_ws)
        total_ws=_cdid_cipw_cells(2*L,L,opts.max_matrix_entries,total_ws)
        total_ws=_cdid_cipw_cells(4,L,opts.max_matrix_entries,total_ws)
        total_ws=_cdid_gmm_workspace(L,P_ret,J_est,opts.max_matrix_entries,total_ws)
    }
    res.workspace_entries_bound=total_ws

    // Fit only primitive rows, class-contiguously when structurally reduced.
    // One current nuisance state per identical recipe; no fit during IF replay.

    res.moments = J(M, 1, .)
    res.moment_table = J(M, 14, .)
    res.links = cdid_chained_ipw_link_meta(M, 1)
    cache = cdid_chained_ipw_cache_entry(n_cache_slots, 1)
    for (c_i = 1; c_i <= n_cache_slots; c_i++) {
        cache[c_i].active = 0
        cache[c_i].op_index = 0
        cache[c_i].stamp = 0
        cache[c_i].pair_a = J(0, 1, 0)
    }
    clock = 0
    last_class=0
    op_frames = cdid_ipw_frame(nop, 1)
    op_frame_built = J(nop, 1, 0)
    n_pool = 0
    if (want_ps_pool) {
        pool_co = J(L, 1, NULL)
        pool_nu = cdid_ipw_nuisance(L, 1)
    }

    for (it = 1; it <= L; it++) {
        ell=res.row_basis.fit_order[it]
        op = C.r_op[ell]
        pid = C.op_pool[op]
        p_sz = C.pool_size[pid]
        pool_k = C.pool_entries[|C.pool_start[pid] \ C.pool_start[pid] + p_sz - 1|]
        coh_idx_op = C.op_g[op] \ pool_k

        _cdid_cipw_query_pair(post, C.r_s[ell], C.r_t[ell], 1, P.n_A, pair_a = ., pair_dy = .)
        if (rows(pair_a) < 1)
            _cdid_cipw_fail("cdid_ipw_pair_support_error", "Requested date pair has zero A-unit observations.", 459)

        pair_all_n=rows(pair_a)
        this_class=res.row_basis.class_id[ell]
        if(this_class>0) {
            risk_keep=_cdid_cipw_which(op_coh_map[a_coh_idx[pair_a],op]:>0)
            pair_a=pair_a[risk_keep];pair_dy=pair_dy[risk_keep]
        }
        hit_slot = 0
        if(this_class>0 & this_class==last_class)hit_slot=-1
        for (c_i = 1; c_i <= n_cache_slots & hit_slot==0; c_i++) {
            if (cache[c_i].active & cache[c_i].op_index == op & rows(cache[c_i].pair_a) == rows(pair_a)) {
                if (_cdid_cipw_all(cache[c_i].pair_a :== pair_a)) {
                    hit_slot = c_i
                    break
                }
            }
        }

        clock++
        if(hit_slot==-1)res.n_nuisance_reused=res.n_nuisance_reused+1
        else if (hit_slot > 0) {
            nu = cache[hit_slot].nu
            cache[hit_slot].stamp = clock
            res.n_nuisance_reused = res.n_nuisance_reused + 1
        }
        else {
            if (!op_frame_built[op]) {
                op_frames[op] = cdid_ipw_frame_create(op_coh_map[a_coh_idx, op], a_known, w_A,
                    P.A_to_F, P.N, 1 + p_sz)
                op_frame_built[op] = 1
            }
            nu = cdid_ipw_nuisance_prepare(op_frames[op], pair_a, X_A[pair_a, .],
                opts.glm_epsilon, opts.glm_maxit)
            res.n_nuisance_fits = res.n_nuisance_fits + 1
            if (want_ps_pool) {
                n_pool++
                pool_co[n_pool] = &(coh_idx_op :+ 0)
                nu_norm = nu
                nu_norm.n_pair_all = nu.n_pair_risk
                nu_norm.in_risk = (1::nu.n_pair_risk)
                pool_nu[n_pool] = nu_norm
            }
            if (n_cache_slots > 0) {
                evict_slot = 0
                min_stamp = clock + 1
                for (c_i = 1; c_i <= n_cache_slots; c_i++) {
                    if (!cache[c_i].active) {
                        evict_slot = c_i
                        break
                    }
                    if (cache[c_i].stamp < min_stamp) {
                        min_stamp = cache[c_i].stamp
                        evict_slot = c_i
                    }
                }
                cache[evict_slot].active = 1
                cache[evict_slot].op_index = op
                cache[evict_slot].stamp = clock
                cache[evict_slot].pair_a = pair_a
                cache[evict_slot].nu = nu
            }
        }

        last_class=this_class
        lk = cdid_ipw_link_outcome(nu, pair_dy)

        res.moments[ell] = lk.estimate
        _cdid_cipw_link_row(res.links[ell],C,D,ell)
        res.links[ell].row_fit_kind="fitted_primitive"
        res.links[ell].n_cohorts = nu.n_cohorts
        res.links[ell].p_covariates = nu.p
        res.links[ell].n_pair_all = pair_all_n
        res.links[ell].n_pair_risk = nu.n_pair_risk
        res.links[ell].n_treated_obs = nu.pair_count[1]
        res.links[ell].n_control_obs = nu.n_pair_risk - nu.pair_count[1]
        res.links[ell].estimate = lk.estimate
        res.links[ell].treated_mean = lk.treated_mean
        res.links[ell].control_mean = lk.control_mean
        res.links[ell].sum_a = nu.sum_a
        res.links[ell].sum_b = nu.sum_b
        res.links[ell].mean_a = nu.mean_a
        res.links[ell].mean_b = nu.mean_b
        res.links[ell].iterations = nu.iterations
        res.links[ell].deviance = nu.deviance
        res.links[ell].score_max = nu.score_max
        res.links[ell].reused_nuisance = (hit_slot != 0)
        res.links[ell].functional_key = D.op_keys[op]
        res.links[ell].cohort_indices = coh_idx_op
        res.links[ell].group_a_count = a_coh_counts[coh_idx_op]
        res.links[ell].group_a_mass = nu.q :* P.n_A
        res.links[ell].pair_count = nu.pair_count
        res.links[ell].pair_mass = nu.pair_mass
        res.links[ell].q = nu.q
        res.links[ell].rates = nu.rates
        res.links[ell].rates_c = nu.rates_c
        res.links[ell].rates_active = nu.rates_active
        res.links[ell].beta = lk.beta
        res.links[ell].G_beta = lk.G_beta
        res.links[ell].G_r = lk.G_r
        res.links[ell].loading_beta = lk.loading_beta
        res.links[ell].loading_rate = lk.loading_rate
        res.links[ell].rate_coef_s0 = lk.rate_coef_s0
        res.links[ell].rate_coef_s1 = lk.rate_coef_s1

        res.moment_table[ell, .] = (ell, res.links[ell].g_label, res.links[ell].s_label,
            res.links[ell].t_label, op, pid, lk.estimate, lk.treated_mean, lk.control_mean,
            pair_all_n, nu.n_pair_risk, res.links[ell].zero_loading,
            res.links[ell].reused_nuisance, nu.iterations)
    }
    if (want_ps_pool & n_pool > 0) {
        *ps_co_out = pool_co[|1\n_pool|]
        *ps_nu_out = pool_nu[|1\n_pool|]
    }

    for(h=1;h<=rows(res.row_basis.dropped);h++) {
        ell=res.row_basis.dropped[h]
        ep=res.row_basis.path_start[ell]
        v=res.row_basis.retained[res.row_basis.e_col[ep]]
        res.links[ell]=res.links[v]
        _cdid_cipw_link_row(res.links[ell],C,D,ell)
        res.links[ell].row_fit_kind="structural_E_image"
        res.links[ell].estimate=0;res.links[ell].treated_mean=0;res.links[ell].control_mean=0
        res.links[ell].G_beta=J(rows(res.links[v].G_beta),1,0)
        res.links[ell].G_r=J(rows(res.links[v].G_r),1,0)
        res.links[ell].loading_beta=J(rows(res.links[v].loading_beta),1,0)
        res.links[ell].loading_rate=J(rows(res.links[v].loading_rate),1,0)
        res.links[ell].rate_coef_s0=J(rows(res.links[v].rate_coef_s0),1,0)
        res.links[ell].rate_coef_s1=J(rows(res.links[v].rate_coef_s1),1,0)
        for(ep=res.row_basis.path_start[ell];ep<res.row_basis.path_start[ell]+res.row_basis.path_count[ell];ep++) {
            v=res.row_basis.retained[res.row_basis.e_col[ep]]
            a=res.row_basis.e_value[ep]
            res.links[ell].estimate=res.links[ell].estimate+a*res.links[v].estimate
            res.links[ell].treated_mean=res.links[ell].treated_mean+a*res.links[v].treated_mean
            res.links[ell].control_mean=res.links[ell].control_mean+a*res.links[v].control_mean
            res.links[ell].G_beta=res.links[ell].G_beta+a*res.links[v].G_beta
            res.links[ell].G_r=res.links[ell].G_r+a*res.links[v].G_r
            res.links[ell].loading_beta=res.links[ell].loading_beta+a*res.links[v].loading_beta
            res.links[ell].loading_rate=res.links[ell].loading_rate+a*res.links[v].loading_rate
            res.links[ell].rate_coef_s0=res.links[ell].rate_coef_s0+a*res.links[v].rate_coef_s0
            res.links[ell].rate_coef_s1=res.links[ell].rate_coef_s1+a*res.links[v].rate_coef_s1
        }
        res.links[ell].n_pair_all=rows(_cdid_cipw_pair_indices(post,C.r_s[ell],C.r_t[ell]))
        res.links[ell].reused_nuisance=1
        res.moments[ell]=res.links[ell].estimate
        res.moment_table[ell,.]=(ell,res.links[ell].g_label,res.links[ell].s_label,res.links[ell].t_label,
            res.links[ell].op_index,res.links[ell].pool_index,res.links[ell].estimate,
            res.links[ell].treated_mean,res.links[ell].control_mean,res.links[ell].n_pair_all,
            res.links[ell].n_pair_risk,res.links[ell].zero_loading,1,res.links[ell].iterations)
    }
    replay_links=res.links[res.row_basis.retained]

    // Release bounded nuisance cache before IF/covariance streaming
    cache = cdid_chained_ipw_cache_entry(0)

    // Point-only omits only the output/diagnostic full-F mean pass. Optimal
    // weighting still uses the complete source A-domain centered Gram below.
    mean_phi=J(1,L,.)
    if((keep_if | opts.materialize_moment_covariance) & !(P_ret>0 & res.weighting!="identity" & !opts.materialize_moment_covariance)) {
    // Pass 1: Stream full-frame unit chunks to compute full-F link IF column means
    sum_phi = J(1, L, 0)
    b_first = 1
    while (b_first <= P.N) {
        b_last = min((P.N, b_first + eff_chunk - 1))
        Phi_chunk = _cdid_cipw_eval_chunk(P, C, post, a_coh_idx, op_coh_map,
            X_A, w_A, replay_links, b_first, b_last)
        sum_phi = sum_phi + colsum(Phi_chunk)
        if (missing(sum_phi))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite IF sum during moment covariance streaming.", 504)
        b_first = b_last + 1
    }
    mean_phi = sum_phi :/ P.N
    res.moment_if_mean = _cdid_cipw_expand(res.row_basis,mean_phi')

    }
    else res.moment_if_mean=J(M,1,.)

    if(opts.materialize_moment_covariance) {
    // Pass 2: Stream full-frame unit chunks to accumulate centered M x M moment covariance
    cov_sum = J(M, M, 0)
    b_first = 1
    while (b_first <= P.N) {
        b_last = min((P.N, b_first + eff_chunk - 1))
        Phi_chunk = _cdid_cipw_eval_chunk(P, C, post, a_coh_idx, op_coh_map,
            X_A, w_A, replay_links, b_first, b_last)
        Phi_centered = _cdid_cipw_expand(res.row_basis,(Phi_chunk :- mean_phi)')'
        if (missing(Phi_centered))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite centered IF during moment covariance streaming.", 504)
        cov_sum = cov_sum + quadcross(Phi_centered, Phi_centered)
        if (missing(cov_sum))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite IF cross-product during moment covariance streaming.", 504)
        b_first = b_last + 1
    }
    for (j = 1; j <= M; j++) {
        for (ell = 1; ell < j; ell++) cov_sum[j, ell] = cov_sum[ell, j]
    }
    res.Vr_estimator = cov_sum :/ (P.N^2)
    res.Vr_analytical = cov_sum :/ P.N
    for (j = 1; j <= M; j++) {
        for (ell = 1; ell < j; ell++) {
            res.Vr_estimator[j, ell] = res.Vr_estimator[ell, j]
            res.Vr_analytical[j, ell] = res.Vr_analytical[ell, j]
        }
    }

        res.moment_covariance_status="materialized_full_F"
    }

    // Build full-column-rank retained design A_retained (M x P_ret) using forest gauge roots
    // for any unanchored connected components; never invent unavailable target levels.
    keep_param_mask = J(P_unred, 1, 1)
    if (D.nullity > 0) {
        unanch_roots = D.comp_roots[D.unanchored]
        res.dropped_gauge_vertices = unanch_roots
        res.dropped_gauge_params = J(D.nullity, 1, 0)
        for (j = 1; j <= D.nullity; j++) {
            v_root = unanch_roots[j]
            p_root = D.v_param[v_root]
            if (p_root < 1 | p_root > P_unred)
                _cdid_cipw_fail("cdid_chained_ipw_design_error", "Unanchored component root must be a non-anchored parameter vertex.", 504)
            keep_param_mask[p_root] = 0
            res.dropped_gauge_params[j] = p_root
        }
    }
    res.kept_param_indices = _cdid_cipw_which(keep_param_mask)
    if (rows(res.kept_param_indices) != P_ret)
        _cdid_cipw_fail("cdid_chained_ipw_design_error", "Retained parameter count does not match design rank.", 504)

    full_to_ret = J(P_unred, 1, 0)
    if (P_ret > 0) full_to_ret[res.kept_param_indices] = (1::P_ret)

    res.A_retained = J(M, P_ret, 0)
    for (ell = 1; ell <= M; ell++) {
        ps = D.l_start[ell]
        pe = D.l_end[ell]
        if (ps > 0) {
            rs = full_to_ret[ps]
            if (rs > 0) res.A_retained[ell, rs] = -1
        }
        if (pe > 0) {
            re = full_to_ret[pe]
            if (re > 0) res.A_retained[ell, re] = 1
        }
    }

    res.H_estimable = J(J_est, P_ret, 0)
    for (j = 1; j <= J_est; j++) {
        u = res.estimable_target_indices[j]
        ps = res.target_param[u]
        if (ps < 1 | ps > P_unred)
            _cdid_cipw_fail("cdid_chained_ipw_design_error", "Estimable target has invalid design parameter index.", 504)
        rs = full_to_ret[ps]
        if (rs < 1 | rs > P_ret)
            _cdid_cipw_fail("cdid_chained_ipw_design_error", "Estimable target vertex was unexpectedly dropped as unanchored gauge.", 504)
        res.H_estimable[j, rs] = 1
    }

    primitive_A=res.A_retained[res.row_basis.retained,.]
    _cdid_cipw_check_design(res.row_basis,res.A_retained)
    primitive_m=res.moments[res.row_basis.retained]
    // Same original catalogue, source primitive covariance, full metric and M.

    if (P_ret > 0) {
        if (res.weighting=="identity") {
            res.gmm = cdid_gmm_identity_map(res.moments,res.A_retained,res.H_estimable,opts.max_matrix_entries)
        }
        else {
            replay_context.population=&P; replay_context.catalog=&C; replay_context.postings=&post
            replay_context.cohorts=&a_coh_idx; replay_context.operators=&op_coh_map
            replay_context.X=&X_A; replay_context.w=&w_A; replay_context.links=&replay_links
            metric_G=(rows(res.row_basis.G) ? res.row_basis.G : I(L))
            if (opts.gmm_weighting == "within") {
                if (opts.within_cov == NULL)
                    _cdid_cipw_fail("cdid_chained_ipw_unsupported_weighting", "Within weighting requires its covariance builder.", 498)
                res.gmm = cdid_gmm_fit_within_norm(primitive_m, primitive_A,
                    (*opts.within_cov)(C, P, X_F, w_F, res.row_basis.retained, opts.glm_epsilon, opts.glm_maxit),
                    metric_G, M, res.H_estimable, opts.condition_trigger, opts.ridge_relative, P.n_A, P.n_A,
                    eff_chunk, &_cdid_cipw_gram_provider(), &replay_context, opts.max_matrix_entries)
            }
            else res.gmm = cdid_gmm_fit_unit_stream(primitive_m,primitive_A,metric_G,M,res.H_estimable,
                opts.gmm_weighting,opts.condition_trigger,opts.ridge_relative,P.n_A,
                eff_chunk,&_cdid_cipw_gram_provider(),&replay_context,opts.max_matrix_entries)
            if(keep_if & !opts.materialize_moment_covariance) {
                mean_phi = res.gmm.unit_gram_mean
                res.moment_if_mean = _cdid_cipw_expand(res.row_basis, mean_phi')
            }
        }
        res.moment_fitted = res.A_retained * res.gmm.coefficients
        res.moment_residual = res.moments-res.moment_fitted
        res.criterion = res.gmm.criterion
        res.gmm_status="solved"
    }
    else {
        // Under PT-all it is possible for all retained links to be zero-loading links
        // between anchored pre-treatment vertices (P_ret == 0, J_est == 0).
        res.moment_fitted = J(M, 1, 0)
        res.moment_residual = res.moments
        if(res.weighting=="identity") {
            res.criterion=quadcross(res.moments,res.moments)
            res.gmm_status="parameter_free_identity_algebra"
        }
        else {
            // Prepared-only K0 has no parameter solve. Its unsupported public
            // request is refused before GMM by the source request validator;
            // never mislabel the identity norm as an optimal criterion.
            res.criterion=.
            res.gmm_status="not_run_no_parameters"
        }
    }

    // Complete-F projected IF stream: target IF retained; coefficient
    // covariance streamed without constructing an M-by-M moment covariance.
    if(P_ret>0) {
        if(J_est>0) res.att[res.estimable_target_indices] = res.gmm.targets
        if(keep_if) {
        if(res.weighting!="identity" & rows(res.gmm.retained_if)==P.n_A & cols(res.gmm.retained_if)==L) {
            if(J_est>0) {
                if(P.n_A==P.N) {
                    res.target_if[.,res.estimable_target_indices]=res.gmm.retained_if*res.gmm.target_loadings'
                }
                else {
                    res.target_if[.,res.estimable_target_indices]=J(P.N,J_est,0)
                    res.target_if[P.A_to_F,res.estimable_target_indices]=((P.N/P.n_A):*res.gmm.retained_if)*res.gmm.target_loadings'
                }
            }
            res.gmm.retained_if=J(0,0,.)
        }
        else {
        coef_cross=(res.weighting=="identity" ? J(P_ret,P_ret,0) : J(0,0,0))
        b_first=1
        while(b_first<=P.N) {
            b_last=min((P.N,b_first+eff_chunk-1))
            Phi_chunk=_cdid_cipw_eval_chunk(P,C,post,a_coh_idx,op_coh_map,X_A,w_A,replay_links,b_first,b_last)
            if(J_est>0)res.target_if[b_first..b_last,res.estimable_target_indices]=Phi_chunk*res.gmm.target_loadings'
            if(res.weighting=="identity") {
                coef_if=(Phi_chunk:-mean_phi)*res.gmm.coefficient_loadings'
                coef_cross=coef_cross+quadcross(coef_if,coef_if)
            }
            b_first=b_last+1
        }
        if(res.weighting=="identity") {
            res.gmm.coefficient_covariance=coef_cross/(P.N^2)
            res.gmm.source_n_units=P.N;res.gmm.normalization_factor=1/(P.N^2)
            res.gmm.covariance_provenance="complete_original_F_projected_IF_sandwich"
        }
        }
        }
    }
    if(keep_if & J_est>0) {
        T_if = res.target_if[., res.estimable_target_indices]
        if (missing(T_if))
            _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite target IF values after GMM loading projection.", 504)
        T_mean = colsum(T_if) :/ P.N
        res.target_if_mean[res.estimable_target_indices] = T_mean
        if (!no_tcov) {
            T_centered = T_if :- T_mean
            T_cross = quadcross(T_centered, T_centered)
            if (missing(T_cross))
                _cdid_cipw_fail("cdid_ipw_numeric_error", "Nonfinite centered target IF cross-product.", 504)
            for (j = 1; j <= J_est; j++) {
                for (ell = 1; ell < j; ell++) T_cross[j, ell] = T_cross[ell, j]
            }
            V_est_sub = T_cross :/ (P.N^2)
            V_asymp_sub = T_cross :/ P.N
            for (j = 1; j <= J_est; j++) {
                for (ell = 1; ell < j; ell++) {
                    V_est_sub[j, ell] = V_est_sub[ell, j]
                    V_asymp_sub[j, ell] = V_asymp_sub[ell, j]
                }
            }
            if(res.weighting=="identity")res.gmm.covariance=V_est_sub
            res.V_estimator_estimable = V_est_sub
            res.V_analytical_estimable = V_asymp_sub
            res.se[res.estimable_target_indices] = sqrt(rowmax((diagonal(V_est_sub), J(J_est, 1, 0))))
            res.V_estimator[res.estimable_target_indices, res.estimable_target_indices] = V_est_sub
            res.V_analytical[res.estimable_target_indices, res.estimable_target_indices] = V_asymp_sub
        }
    }

    if (K > 0) {
        res.target_ledger[., 13] = res.att
        res.target_ledger[., 14] = res.se
    }

    if (K == 0 | J_est == K) {
        res.status = "ok"
        res.global_reason = ""
    }
    else if (J_est > 0) {
        res.status = "partial_unavailable"
        res.global_reason = "some requested post target cells are unavailable in the comparison catalog or disconnected in the PT design"
    }
    else {
        res.status = "all_unavailable"
        res.global_reason = "every requested post target cell is unavailable in the comparison catalog or disconnected in the PT design"
    }

    res.work_ops = total_work
    return(res)
}

struct cdid_chained_ipw_result scalar cdid_chained_ipw_fit(
    struct cdid_population scalar P,
    struct cdid_cat scalar C,
    struct cdid_design scalar D,
    real matrix X_F,
    real colvector w_F,
    real matrix requested_targets)
{
    return(cdid_chained_ipw_fit_opts(P, C, D, X_F, w_F, requested_targets,
        cdid_chained_ipw_default_opts()))
}

/* Complete population-to-estimator convenience wrapper (no stubs). */
struct cdid_chained_ipw_result scalar cdid_chained_ipw_from_population(
    struct cdid_population scalar P,
    real matrix X_F,
    real colvector w_F,
    string scalar control_group,
    real scalar anticipation,
    real matrix requested_targets,
    | struct cdid_chained_ipw_opts user_opts)
{
    struct cdid_chained_ipw_opts scalar opts
    struct cdid_catalog_opts scalar cat_opts
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    real colvector coh, tgt_coh, w_use
    real matrix X_use, req

    opts = (args() >= 7 ? user_opts : cdid_chained_ipw_default_opts())
    _cdid_cipw_check_opts(opts)

    coh = cdid_cipw_cohorts(P)
    tgt_coh = cdid_cipw_targets(P)
    req = requested_targets
    if (rows(req) == 0) {
        req = cdid_cipw_post_targets(P.calendar, tgt_coh)
    }
    X_use = (cols(X_F) == 0 & rows(X_F) == 0 ? J(P.N, 1, 1) : X_F)
    w_use = (rows(w_F) == 0 ? J(P.N, 1, 1) : w_F)

    cat_opts = cdid_catalog_default_opts()
    cat_opts.conventions = opts.conventions
    cat_opts.mem_limit_bytes = min((2^30, 8 * opts.max_matrix_entries))
    cat_opts.work_cap = opts.work_cap

    pairs = cdid_cipw_pairs_opts(P, coh, cat_opts)
    C = cdid_catalog_build_opts(P.calendar, coh, tgt_coh, anticipation,
        pairs.occupancy, control_group, req, cat_opts)
    D = cdid_chained_ipw_build_design(C, req, opts.pt_assumption, opts.pt_start)
    return(cdid_chained_ipw_fit_opts(P, C, D, X_use, w_use, req, opts))
}

/* Complete raw-data convenience wrapper (numeric IDs, static covariates/weights, no stubs). */
struct cdid_chained_ipw_result scalar cdid_chained_ipw_from_data(
    real colvector data_id,
    real colvector date,
    real colvector y,
    real colvector g,
    real matrix data_X,
    real colvector data_w,
    real colvector frame_id,
    real colvector frame_g,
    real matrix frame_X,
    real colvector frame_w,
    real scalar frame_supplied,
    real scalar frame_has_g,
    real scalar X_in_data,
    real scalar X_in_frame,
    real scalar w_in_data,
    real scalar w_in_frame,
    string scalar control_group,
    real scalar anticipation,
    real matrix requested_targets,
    | struct cdid_chained_ipw_opts user_opts)
{
    struct cdid_chained_ipw_opts scalar opts
    struct cdid_population_ids scalar ids
    struct cdid_population scalar P
    struct cdid_population_gate scalar gate
    struct cdid_population_weight scalar pw
    string scalar fix_w
    real scalar i, u, p_raw
    real matrix X_F, data_design, frame_design

    opts = (args() >= 20 ? user_opts : cdid_chained_ipw_default_opts())
    _cdid_cipw_check_opts(opts)

    ids = cdid_population_ids_numeric(data_id, frame_id, frame_supplied)
    P = cdid_population_prepare(ids, date, y, g, frame_g, frame_has_g)

    fix_w = (opts.fix_weights == "static" ? "" : opts.fix_weights)
    pw = cdid_population_weights(P, data_w, frame_w, w_in_data, w_in_frame, fix_w, 0)
    if (pw.reduction != "static" & pw.reduction != "first_period_static")
        _cdid_cipw_fail("cdid_chained_ipw_unsupported_weights", "Baseline static-input IPW refuses reference-dated weights.", 459)

    if (!X_in_data & !X_in_frame) {
        X_F = J(P.N, 1, 1)
    }
    else {
        if (X_in_data & rows(data_X) != P.R)
            _cdid_cipw_fail("cdid_chained_ipw_covariate_error", "data_X row count must equal panel row count P.R.", 459)
        if (frame_supplied & X_in_frame & rows(frame_X) != P.N)
            _cdid_cipw_fail("cdid_chained_ipw_covariate_error", "frame_X row count must equal full-frame count P.N.", 459)
        p_raw = (X_in_frame ? cols(frame_X) : cols(data_X))
        if (p_raw < 1)
            _cdid_cipw_fail("cdid_chained_ipw_covariate_error", "Supplied covariate matrix must have at least 1 column.", 459)
        if (X_in_data & X_in_frame & cols(data_X) != cols(frame_X))
            _cdid_cipw_fail("cdid_chained_ipw_covariate_error", "data_X and frame_X column counts disagree.", 459)

        if (frame_supplied & X_in_frame) {
            frame_design = (J(P.N, 1, 1), frame_X)
        }
        else {
            if (P.n_frame_only > 0)
                _cdid_cipw_fail("cdid_chained_ipw_covariate_error", "Frame-only units require frame_X covariates on the full-F static route.", 459)
            frame_design = J(P.N, 1 + p_raw, .)
            for (u = 1; u <= P.N; u++) {
                i = P.first_row[u]
                if (i > 0) frame_design[u, .] = (1, data_X[i, .])
            }
        }
        data_design = (X_in_data ? (J(P.R, 1, 1), data_X) : J(0, 0, .))
        gate = cdid_population_static_gate(P, frame_design, data_design, J(P.N, 1, 1),
            X_in_data, 0, 0, 1)
        if (!gate.accepted)
            _cdid_cipw_fail(gate.error_class, gate.message, 459)
        X_F = frame_design
    }

    return(cdid_chained_ipw_from_population(P, X_F, pw.owner_weight,
        control_group, anticipation, requested_targets, opts))
}

end
