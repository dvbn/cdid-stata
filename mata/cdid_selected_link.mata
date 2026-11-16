*! version 1.1.0  01nov2026
*! Selected IPW, Control-Change OR, and Treatment DR Link Composition
version 19.5

mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_slink(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}

/* --------------------------------------------------------------------------
   1. Struct Definitions
   -------------------------------------------------------------------------- */

struct cdid_selected_controls {
    string scalar est_method
    string scalar selection
    string scalar selection_model
    string scalar selection_pooling
    string scalar outcome_pooling
    real scalar   or_intercept
    real scalar   eps
    real scalar   maxit
    real scalar   rank_tol
    real scalar   max_chunk_rows
    real scalar   max_work
    real scalar   max_workspace_doubles
    real scalar   points_only
}

/* Production input: one matrix per unique modeled date/role, one vector
   per weight anchor, and one known transport vector per declared universe.
   Recipes contain IDs only. calendar_index=0 denotes a static profile. */
struct cdid_selected_matrix_slot {
    real scalar calendar_index
    real matrix value
}
struct cdid_selected_registry {
    struct cdid_selected_matrix_slot colvector sel_data, or_data, trt_data
    struct cdid_selected_matrix_slot colvector weight_data, rho_data
}
struct cdid_selected_recipe_spec {
    real scalar recipe_id, design_cal, weight_cal
    real scalar sel_id, or_id, trt_id, weight_id, rho_id
    real scalar sel_has_intercept, sel_int_col, has_xformla, has_rho
    string scalar recipe_key, universe, justification
}

/* Explicit-array diagnostic convenience only. Production callers use the
   registry/spec API and never construct full-F arrays for every recipe. */
struct cdid_selected_recipe_in {
    real scalar   recipe_id
    string scalar recipe_key
    real scalar   design_cal
    real scalar   weight_cal
    real matrix   Z_sel
    real scalar   sel_has_intercept
    real scalar   sel_int_col
    real scalar   has_xformla
    real scalar   has_rho
    real matrix   Z_or
    real matrix   V_trt
    real colvector w
    real colvector B
    real colvector U
    real colvector rho
    string scalar universe
    string scalar justification
}

struct cdid_selected_rows {
    real colvector loc, unit, dy, q, cohort_local
    real matrix score, R_ipw, R_dr
}

struct cdid_selected_pair_in {
    real scalar    s_cal
    real scalar    t_cal
    real scalar    s_val
    real scalar    t_val
}

struct cdid_selected_link_spec {
    string scalar  link_id
    real scalar    pair_id
    real scalar    s_cal
    real scalar    t_cal
    real scalar    target_cohort_pos
    real colvector control_cohort_pos
    real scalar    recipe_id
    string scalar  functional_key
}

// Sparse consumer derivatives in selection parameter coordinates.
struct cdid_selected_loaded_block {
    real colvector ids, col_pos
    real matrix values
}

struct cdid_selected_loading_piece {
    real colvector theta_ids, columns
    real scalar cohort
    real matrix G_theta
    real rowvector G_ak
}

struct cdid_selected_universe_buffer {
    real colvector value
}

struct cdid_selected_context {
    real scalar                 recipe_id
    string scalar               recipe_key
    real scalar                 design_cal
    real scalar                 weight_cal
    struct cdid_sel_design_spec scalar reader_spec
    real scalar reader_n_U, universe_code
    real scalar sel_id, or_id, trt_id, weight_id, rho_id
    real colvector mass_F, mass_U, repeat_mass, a_U
    string scalar universe, justification
    struct cdid_sel_cell_map scalar cell_map
    struct cdid_sel_fit_result scalar sel_fit
    struct cdid_sel_up_map scalar up_map
    struct cdid_selected_loading_piece colvector loading_pieces
    real scalar n_loading_pieces
    real scalar                 m_col_beta_start
    real scalar                 m_col_beta_end
    real scalar                 m_col_obj
    real scalar                 n_moments_total
    struct cdid_sel_moments scalar acc_full
    real scalar                 acc_full_ready
    real matrix                 pc_lookup
}

struct cdid_selected_support {
    real scalar    support_id
    real scalar    pair_id
    real scalar    cohort_pos
    real scalar    recipe_id
    real scalar    pc_id
    real scalar    a_k
    real scalar    mass_U_k
    real scalar    repeat_mass_k
    real scalar    rate_baseline
    real scalar    rate_c
    real scalar    rate_active
    real scalar    weight_mass
    real scalar    n_obs
    real scalar    sum_y
    real rowvector sum_z
    real rowvector xmin
    real rowvector xmax
    real rowvector first_z
    real rowvector colmax
    real matrix    F_qr
    struct cdid_sel_log_deriv scalar ld
    real colvector ld_g_idx
    real scalar    ld_ready
    real matrix    Z_or
}

struct cdid_selected_or_block {
    real scalar    block_id
    real scalar    pair_id
    real scalar    recipe_id
    real scalar    L
    real colvector control_cohort_pos
    real colvector support_ids
    real scalar    alpha
    real scalar    weight_sum
    real scalar    mass
    real scalar    mean_y
    real rowvector mean_d
    real scalar    p_kept
    real scalar    rss
    real scalar    objective
    real scalar    verif_scan
    real scalar    beta_col_start
    real scalar    beta_col_end
    real rowvector kept
    real rowvector dropped
    real matrix    alias_map
    real matrix    alias_constants
    real colvector beta
    real matrix    R
    real rowvector pivot
    real scalar    score_mean_alpha
    real rowvector score_mean_beta
    real scalar    m_col_alpha
    real scalar    m_col_beta_start
    real scalar    m_col_beta_end
    real colvector baseline_alpha_rate
    real matrix    baseline_beta_rate
    real colvector blk_obs_u
    real colvector blk_ha
    real matrix    blk_hb
    real matrix    blk_gram_inv
    real scalar    blk_cache_ready
}

struct cdid_sel_sep_blk_chunk {
    real colvector ha, xa
    real matrix    hb, xb
}

struct cdid_selected_or_shared {
    real scalar    intercept
    real scalar    p_raw
    real scalar    p_kept
    real scalar    rss
    real scalar    objective
    real scalar    verification_scan
    real rowvector kept
    real rowvector dropped
    real matrix    alias_map
    real matrix    alias_constants
    real colvector beta
    real matrix    R
    real rowvector pivot
    real matrix    beta_score_F
    real matrix    baseline_beta_rate_loading
    real colvector baseline_objective_rate_loading
    real matrix    shared_gram_inv
    real scalar    shared_gram_ready
}

struct cdid_selected_trt_model {
    real scalar    trt_id
    real scalar    pair_id
    real scalar    recipe_id
    real scalar    block_id
    real scalar    target_cohort_pos
    real scalar    n_risk
    real scalar    p_trt
    real colvector cohorts
    real colvector support_ids
    real colvector pc_ids
    real scalar    sum_wp
    real scalar    sum_wm
    real scalar    dp
    real scalar    dm
    real scalar    T_ipw
    real scalar    M_ipw
    real scalar    mu_ipw
    real scalar    T_dr
    real scalar    M_dr
    real scalar    mu_dr
    real scalar    mu_dr_decomp
    real colvector d_bar_plus
    real colvector d_bar_minus
    real colvector d_bar_diff
    real colvector beta_trt
    real matrix    H_chol
    real matrix    V_ipw
    real matrix    V_dr
    real scalar    m_col_C_start
    real scalar    m_col_C_end
    real scalar    m_col_U_ipw_start
    real scalar    m_col_U_ipw_end
    real scalar    m_col_U_dr_start
    real scalar    m_col_U_dr_end
    real matrix    baseline_C
    real matrix    baseline_U_ipw
    real matrix    baseline_U_dr
    real scalar    rep_tid
    real colvector fit_u_cache
    real colvector fit_q_cache
    real scalar    fit_q_sum
    real matrix    H_inv
    real colvector w_gamma_ipw
    real colvector w_gamma_dr
    real colvector phi_risk_ipw
    real colvector phi_risk_dr
    real scalar    m_col_fast
    real colvector baseline_fast_ipw
    real colvector baseline_fast_dr
    struct cdid_selected_rows scalar risk_cache
    real scalar    risk_cache_ready
}

struct cdid_selected_link_rec {
    real scalar    link_index
    string scalar  link_id
    string scalar  functional_key
    real scalar    pair_id
    real scalar    recipe_id
    real scalar    target_cohort_pos
    real scalar    target_support_id
    real scalar    block_id
    real scalar    trt_id
    real scalar    L
    real colvector control_cohort_pos
    real colvector control_support_ids
    real scalar    mu
    real scalar    mu_ipw
    real scalar    mu_or
    real scalar    mu_dr
    real scalar    d_target
    real scalar    num_target
    real rowvector c_target
    real scalar    m_col_or_mu
    real scalar    m_col_or_den
    real scalar    m_col_or_des_start
    real scalar    m_col_or_des_end
    real scalar    des_col_start
    real scalar    des_col_end
    real scalar    baseline_or_rate_loading
    real scalar    baseline_den_rate_loading
    real rowvector baseline_design_rate_loading
    real colvector target_mu_vec
    real scalar    target_mu_ready
}

struct cdid_selected_model {
    real scalar                          prepared
    real matrix                          recipe_keep
    string scalar                        outcome_pooling
    pointer() scalar                     population_ref
    struct cdid_sel_reader scalar         reader
    struct cdid_selected_registry scalar profiles
    real scalar reader_sel_slot, reader_w_slot, reader_rho_slot
    struct cdid_selected_universe_buffer colvector universe_buffers
    real scalar                          reader_ready, reader_active, reader_universe
    struct cdid_selected_controls scalar ctl
    real scalar                          N
    real scalar                          n_A
    real scalar                          lambda
    real scalar                          n_cohorts
    real scalar                          n_dates
    real scalar                          R
    real scalar                          n_pairs
    real scalar                          n_supports
    real scalar                          n_blocks
    real scalar                          n_trt_models
    real scalar                          n_links
    real scalar                          total_p_sep
    real scalar                          total_des_sep
    real colvector                       is_A
    real colvector                       unit_to_a
    real colvector                       a_to_unit
    real colvector                       cohort_idx
    real colvector                       cohort_values
    real colvector                       calendar_values
    real colvector                       w_pop
    real colvector                       tau
    real colvector                       eta
    real scalar                          mean_w
    real scalar                          mean_wA
    struct cdid_selected_context colvector contexts
    struct cdid_selected_pair_in colvector pairs
    struct cdid_selected_support colvector supports
    struct cdid_selected_or_block colvector or_blocks
    struct cdid_selected_or_shared scalar  or_shared
    struct cdid_selected_trt_model colvector trt_models
    struct cdid_selected_link_rec colvector links
    real colvector                       estimates
    real colvector                       mu_ipw_vec
    real colvector                       mu_or_vec
    real colvector                       mu_dr_vec
    real scalar                          n_selection_fit_calls
    real scalar                          n_selection_blocks
    real scalar                          n_treatment_fit_calls
    real scalar                          n_or_block_fits
    real scalar                          n_imputation_calls
    real scalar                          workspace_doubles
    real scalar                          work_ops
    struct cdid_selected_rows colvector  pair_full
    real scalar                          pair_full_ready
    struct cdid_selected_rows colvector  support_full
    real scalar                          support_full_ready
}

struct cdid_selected_chunk_out {
    real colvector unit_ids
    real colvector link_ids
    real colvector alpha_block_ids
    real scalar    scale_core
    real scalar    normalization_count
    real matrix    phi
    real matrix    psi_dir
    real matrix    psi_or
    real matrix    psi_ps
    real matrix    psi_sel
    real matrix    psi_ipw
    real matrix    psi_or_estimator
    real matrix    xi_beta
    real matrix    beta_score
    real matrix    xi_alpha
    real matrix    alpha_score
    real matrix    den_if
    real matrix    num_if
    real matrix    design_if
    real colvector objective_if
    real matrix    coord_ipw
    real matrix    coord_dr
    real matrix    xi_gamma
    real matrix    norm_A_if
    real matrix    tau_if
    real matrix    eta_if
    real colvector lambda_if
}

struct cdid_selected_link_out {
    real colvector estimates
    real colvector mu_ipw
    real colvector mu_or
    real colvector mu_dr
    real colvector or_beta
    real colvector or_alpha
    real scalar    or_objective
    real matrix    phi
    real matrix    phi_core
    real matrix    psi_dir
    real matrix    psi_or
    real matrix    psi_ps
    real matrix    psi_sel
    real matrix    psi_ipw
    real matrix    psi_or_estimator
    real matrix    xi_beta
    real matrix    beta_score
    real matrix    xi_alpha
    real matrix    alpha_score
    real matrix    den_if
    real matrix    num_if
    real matrix    design_if
    real colvector objective_if
    real matrix    coord_ipw
    real matrix    coord_dr
    real matrix    xi_gamma
    real matrix    norm_A_if
    real matrix    tau_if
    real matrix    eta_if
    real colvector lambda_if
    real matrix    V_unit
    real matrix    Omega_unit
    real matrix    V_cluster
    real matrix    Omega_cluster
    real scalar    n_selection_fit_calls
    real scalar    n_selection_blocks
    real scalar    n_treatment_fit_calls
    real scalar    n_or_block_fits
    real scalar    n_imputation_calls
}

/* --------------------------------------------------------------------------
   2. Error, Budget, and Utility Helpers
   -------------------------------------------------------------------------- */

void _cdid_sel_link_fail(string scalar kind, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_sel_link_index(real matrix mask)
{
    if (rows(mask)*cols(mask) == 0) return(J(0, 1, 0))
    return(vec(selectindex(mask :!= 0)))
}

void _cdid_sel_link_budget_check(real scalar ws_doubles, real scalar work_ops,
    real scalar max_work, real scalar max_workspace_doubles, string scalar stage)
{
    if (missing(ws_doubles) | ws_doubles < 0 | ws_doubles > 2^53-1) {
        _cdid_sel_link_fail("budget", "Workspace size is not a finite exact allocation count.", 460)
    }
    if (max_workspace_doubles != .) {
        if (max_workspace_doubles < 1 | max_workspace_doubles >= . |
            max_workspace_doubles != floor(max_workspace_doubles))
            _cdid_sel_link_fail("input", "Workspace cap must be an explicit positive integer or . for no cap.", 498)
        if (ws_doubles > max_workspace_doubles)
            _cdid_sel_link_fail("budget", "Workspace exceeds the caller's explicit allocation cap at "+stage+".", 460)
    }
    if (max_work < 0 | (max_work > 0 & max_work >= .)) {
        _cdid_sel_link_fail("input", "max_work must be 0 (unlimited) or a positive finite cap.", 498)
    }
    if (max_work > 0) {
        if (missing(work_ops) | work_ops < 0 | work_ops > max_work) {
            _cdid_sel_link_fail("budget",
                sprintf("Work operation cap (%g) exceeded at %s.", max_work, stage),
                460)
        }
    }
}

struct cdid_selected_controls scalar cdid_selected_default_controls()
{
    struct cdid_selected_controls scalar ctl

    ctl.est_method        = "dr"
    ctl.selection         = "conditional_mtar"
    ctl.selection_model   = "pair_logit"
    ctl.selection_pooling = "separate"
    ctl.outcome_pooling   = "separate"
    ctl.or_intercept      = 1
    ctl.eps               = 1e-12
    ctl.maxit             = 50
    ctl.rank_tol          = 1e-7
    ctl.max_chunk_rows    = 0
    ctl.max_work          = 0
    ctl.max_workspace_doubles = .
    ctl.points_only       = 0
    return(ctl)
}

struct cdid_selected_controls scalar _cdid_sel_link_norm_controls(
    struct cdid_selected_controls scalar raw)
{
    struct cdid_selected_controls scalar ctl

    ctl = raw
    if (ctl.est_method == "reg") ctl.est_method = "or"
    if (ctl.est_method != "ipw" & ctl.est_method != "or" & ctl.est_method != "dr") {
        _cdid_sel_link_fail("input",
            "est_method must be 'ipw', 'or' ('reg'), or 'dr'.", 459)
    }
    if (ctl.selection == "") ctl.selection = "conditional_mtar"
    if (ctl.selection == "mtar") ctl.selection = "conditional_mtar"
    if (ctl.selection != "conditional_mtar" &
        ctl.selection != "corollary1" &
        ctl.selection != "none") {
        _cdid_sel_link_fail("input",
            "selection must be 'conditional_mtar', 'corollary1', or 'none'.", 459)
    }
    if (ctl.selection == "none") {
        if ((ctl.selection_model != "" & ctl.selection_model != "pair_logit") |
            (ctl.selection_pooling != "" & ctl.selection_pooling != "separate")) {
            _cdid_sel_link_fail("unsupported",
                "selmodel(), selpooling(), selcovars(), seluniverse(), selrho() and seljustification() require selection(conditional_mtar).", 459)
        }
        ctl.selection_model   = "pair_logit"
        ctl.selection_pooling = "separate"
    }
    else {
        if (ctl.selection_model == "") ctl.selection_model = "pair_logit"
        if (ctl.selection_pooling == "") ctl.selection_pooling = "separate"
        if (ctl.selection_model != "pair_logit" &
            ctl.selection_model != "independent_period_logit") {
            _cdid_sel_link_fail("input",
                "selection_model must be 'pair_logit' or 'independent_period_logit'.", 459)
        }
        if (ctl.selection == "corollary1" & ctl.selection_model != "pair_logit") {
            _cdid_sel_link_fail("unsupported",
                "selection(corollary1) is allowed only with selmodel(pair_logit); specify selection(conditional_mtar).", 459)
        }
        if (ctl.selection == "corollary1") ctl.selection = "conditional_mtar"
        if (ctl.selection_pooling != "separate" &
            ctl.selection_pooling != "within_cohort" &
            ctl.selection_pooling != "across_cohorts") {
            _cdid_sel_link_fail("input",
                "selpooling() must be separate, within_cohort, or across_cohorts.", 459)
        }
    }
    if (ctl.or_intercept != 0 & ctl.or_intercept != 1) {
        _cdid_sel_link_fail("input", "or_intercept must be 0 or 1.", 459)
    }
    if (ctl.outcome_pooling == "") ctl.outcome_pooling = "separate"
    if (ctl.outcome_pooling != "separate" & ctl.outcome_pooling != "pooled") {
        _cdid_sel_link_fail("input",
            "outcome_pooling must be 'separate' or 'pooled'.", 459)
    }
    if (ctl.eps <= 0 | ctl.eps >= 1 | missing(ctl.eps)) {
        _cdid_sel_link_fail("input", "eps must be in (0, 1).", 459)
    }
    if (ctl.maxit < 1 | ctl.maxit >= . | ctl.maxit != floor(ctl.maxit)) {
        _cdid_sel_link_fail("input", "maxit must be a positive integer.", 459)
    }
    if (ctl.rank_tol <= 0 | ctl.rank_tol >= 1 | missing(ctl.rank_tol)) {
        _cdid_sel_link_fail("input", "rank_tol must be in (0, 1).", 459)
    }
    if (ctl.max_chunk_rows < 0 | ctl.max_chunk_rows >= . |
        ctl.max_chunk_rows != floor(ctl.max_chunk_rows)) {
        _cdid_sel_link_fail("input", "max_chunk_rows must be a nonnegative integer.", 459)
    }
    if (ctl.max_work < 0 | ctl.max_work >= .) {
        _cdid_sel_link_fail("input", "max_work must be nonnegative and finite.", 459)
    }
    if (ctl.points_only != 0 & ctl.points_only != 1) {
        _cdid_sel_link_fail("input", "points_only must be 0 or 1.", 459)
    }
    return(ctl)
}

/* Match sorted query vector against sorted table (strictly increasing);
   returns positions in table (0 where absent). */
real colvector _cdid_sel_link_match_sorted(real colvector query, real colvector table)
{
    real colvector out
    real scalar n_q, n_t, i, k

    n_q = rows(query)
    n_t = rows(table)
    out = J(n_q, 1, 0)
    if (n_q == 0 | n_t == 0) return(out)
    for (i = 1; i <= n_q; i = i + 1) {
        k = _cdid_ipw_lower_bound(table, query[i])
        if (k <= n_t) {
            if (table[k] == query[i]) out[i] = k
        }
    }
    return(out)
}

/* Table-convenience metadata only: identify unique calendar pairs referenced
   by links. Production supplies opaque pair positions directly. */
struct cdid_selected_pair_in colvector cdid_selected_pairs_from_matrix(
    real matrix Y_obs, real matrix obs_mat, real colvector calendar_values,
    struct cdid_selected_link_spec colvector link_specs)
{
    struct cdid_selected_pair_in colvector pairs
    real matrix pair_keys
    real scalar n_links, l, s_c, t_c, n_p, p, N, M

    N = rows(obs_mat)
    M = cols(obs_mat)
    n_links = rows(link_specs)
    if (rows(Y_obs) != N | cols(Y_obs) != M | rows(calendar_values) != M | n_links < 1) {
        _cdid_sel_link_fail("input", "Invalid Y_obs, obs_mat, or calendar_values dimensions.", 459)
    }
    pair_keys = J(n_links, 2, 0)
    for (l = 1; l <= n_links; l = l + 1) {
        s_c = link_specs[l].s_cal
        t_c = link_specs[l].t_cal
        if (s_c < 1 | t_c <= s_c | t_c > M | s_c != floor(s_c) | t_c != floor(t_c)) {
            _cdid_sel_link_fail("support", "Pair calendar indices must satisfy 1 <= s_cal < t_cal <= M.", 459)
        }
        pair_keys[l, ] = (s_c, t_c)
    }
    pair_keys = uniqrows(sort(pair_keys, (1, 2)))
    n_p = rows(pair_keys)
    pairs = cdid_selected_pair_in(n_p, 1)
    for (p = 1; p <= n_p; p = p + 1) {
        s_c = pair_keys[p, 1]
        t_c = pair_keys[p, 2]
        pairs[p].s_cal = s_c
        pairs[p].t_cal = t_c
        pairs[p].s_val = calendar_values[s_c]
        pairs[p].t_val = calendar_values[t_c]
    }
    return(pairs)
}

void _cdid_sld_slot_check(real scalar id, real scalar count)
{
    if(id<1 | id>count | id!=floor(id) | id>=.)
        _cdid_sel_link_fail("input","Active profile ID is outside the supplied unique-date registry.",459)
}
/* Diagnostic adapter deduplicates identical supplied arrays. Production
   supplies unique profile slots directly and does not build recipe arrays. */
real scalar _cdid_sld_slot_add(struct cdid_selected_matrix_slot colvector slots,
    real scalar used, real scalar date, real matrix value)
{
    real scalar j
    for(j=1;j<=used;j++) {
        if(slots[j].calendar_index!=date | rows(slots[j].value)!=rows(value) | cols(slots[j].value)!=cols(value)) continue
        if(all(vec(slots[j].value):==vec(value))) return(j)
    }
    used=used+1;slots[used].calendar_index=date;slots[used].value=value
    return(used)
}
/* Registry values are retained once; the sole reader leases only
   its active profile. No context owns an N-row design/weight/unit copy. */
real matrix _cdid_sld_sel_full(struct cdid_selected_model scalar m, real scalar id)
{
    return(id==m.reader_sel_slot ? m.reader.Z : m.profiles.sel_data[id].value)
}
real colvector _cdid_sld_w_full(struct cdid_selected_model scalar m, real scalar id)
{
    return(id==m.reader_w_slot ? m.reader.w : m.profiles.weight_data[id].value)
}
real colvector _cdid_sld_rho_full(struct cdid_selected_model scalar m, real scalar id)
{
    return(id==m.reader_rho_slot ? m.reader.rho : m.profiles.rho_data[id].value)
}
real matrix _cdid_sld_or_rows(struct cdid_selected_model scalar m, real scalar r,
    real colvector units, | real rowvector columns)
{
    real scalar id
    id=m.contexts[r].or_id
    if(args()==3) return(m.profiles.or_data[id].value[units,])
    return(m.profiles.or_data[id].value[units,columns])
}
real matrix _cdid_sld_trt_rows(struct cdid_selected_model scalar m, real scalar r, real colvector units)
{
    return(m.profiles.trt_data[m.contexts[r].trt_id].value[units,])
}
real colvector _cdid_sld_w_rows(struct cdid_selected_model scalar m, real scalar r, real colvector units)
{
    real scalar id
    id=m.contexts[r].weight_id
    if(id==m.reader_w_slot) return(m.reader.w[units])
    return(m.profiles.weight_data[id].value[units])
}
/* Restricted reference eligibility (missingreference(exclude)): 1 for the units
   inside recipe r's eligibility context; every unit when unrestricted. */
real colvector _cdid_sld_keep(struct cdid_selected_model scalar m, real scalar r, real colvector units)
{
    if (rows(m.recipe_keep) != m.N | cols(m.recipe_keep) < r | rows(units) == 0) return(J(rows(units), 1, 1))
    return(m.recipe_keep[units, r])
}
void _cdid_selected_activate_reader(struct cdid_selected_model scalar m, real scalar r)
{
    real scalar id,u
    if(m.reader_active==r) return
    id=m.contexts[r].sel_id
    if(id!=m.reader_sel_slot) {
        if(m.reader_sel_slot) swap(m.reader.Z,m.profiles.sel_data[m.reader_sel_slot].value)
        swap(m.reader.Z,m.profiles.sel_data[id].value);m.reader_sel_slot=id
    }
    id=m.contexts[r].weight_id
    if(id!=m.reader_w_slot) {
        if(m.reader_w_slot) swap(m.reader.w,m.profiles.weight_data[m.reader_w_slot].value)
        swap(m.reader.w,m.profiles.weight_data[id].value);m.reader_w_slot=id
    }
    id=m.contexts[r].rho_id
    if(id!=m.reader_rho_slot) {
        if(m.reader_rho_slot) swap(m.reader.rho,m.profiles.rho_data[m.reader_rho_slot].value)
        swap(m.reader.rho,m.profiles.rho_data[id].value);m.reader_rho_slot=id
    }
    u=m.contexts[r].universe_code
    if(m.reader_universe!=u) {
        if(m.reader_universe>0) swap(m.reader.U,m.universe_buffers[m.reader_universe].value)
        swap(m.reader.U,m.universe_buffers[u].value);m.reader_universe=u
    }
    m.reader.spec=m.contexts[r].reader_spec;m.reader.n_U=m.contexts[r].reader_n_U
    m.reader.selection_universe=m.contexts[r].universe;m.reader.justification=m.contexts[r].justification
    m.reader.mass_F=m.contexts[r].mass_F;m.reader.mass_U=m.contexts[r].mass_U
    m.reader.repeat_mass=m.contexts[r].repeat_mass;m.reader.a_U=m.contexts[r].a_U
    m.reader_active=r
}
/* A recipe restricted to a reference eligibility context reads its own roster through a
   dedicated universe buffer. */
void _cdid_selected_use_universe(struct cdid_selected_model scalar m, real scalar r, real colvector U_r)
{
    struct cdid_selected_universe_buffer colvector ubuf
    real scalar id
    if(m.reader_universe>0) swap(m.reader.U,m.universe_buffers[m.reader_universe].value)
    m.reader_universe=0;m.reader_active=0
    ubuf=cdid_selected_universe_buffer(rows(m.universe_buffers)+1,1)
    for(id=1;id<=rows(m.universe_buffers);id++) ubuf[id].value=m.universe_buffers[id].value
    ubuf[rows(ubuf)].value=U_r
    m.universe_buffers=ubuf
    m.contexts[r].universe_code=rows(m.universe_buffers)
    _cdid_selected_activate_reader(m,r)
}
void _cdid_selected_reader_store(struct cdid_selected_model scalar m, real scalar r,
    struct cdid_sel_reader scalar reader)
{
    if(!m.reader_ready) {
        m.reader=reader
        m.reader.Z=J(0,0,0);m.reader.w=m.reader.rho=m.reader.U=J(0,1,0)
        m.universe_buffers=cdid_selected_universe_buffer(3,1)
        m.universe_buffers[1].value=J(m.N,1,1)
        m.universe_buffers[2].value=reader.B;m.universe_buffers[3].value=reader.A
        m.reader_ready=1;m.reader_active=0;m.reader_universe=0
    }
    m.contexts[r].reader_spec=reader.spec;m.contexts[r].reader_n_U=reader.n_U
    m.contexts[r].mass_F=reader.mass_F;m.contexts[r].mass_U=reader.mass_U
    m.contexts[r].repeat_mass=reader.repeat_mass;m.contexts[r].a_U=reader.a_U
    m.contexts[r].universe=reader.selection_universe;m.contexts[r].justification=reader.justification
    m.contexts[r].universe_code=(reader.selection_universe=="F" ? 1 : (reader.selection_universe=="B" ? 2 : 3))
    _cdid_selected_activate_reader(m,r)
}

/* Used-scope reader for recipe-indexed (dated) registry slots. Additive twin of
   cdid_sel_reader_pop (unchanged): same struct, original F/U/A/B, cohort
   positions, declaration and known-rho checks. Z and w are validated, and the
   recipe law mass_U/repeat_mass/a_U is formed exactly as reader_pop does, only
   on the U rows of cohorts in cohort_scope (the recipe's own selection cells;
   R chained_reference_inputs.R:611-628). Other cohorts keep the kernel's
   absent-cohort encoding 0: their cells/pairs refuse and their a_k^U IF channel
   is skipped (R has no such rate coordinate). Out-of-scope rows are never read
   or altered. mass_F is a non-computational descriptor; dated defines it only
   when U=F. */
struct cdid_sel_reader scalar cdid_sel_reader_scope(struct cdid_population scalar pop,
    real matrix Z, real scalar has_intercept, real scalar int_col,
    real colvector w, struct cdid_sel_decl scalar decl,
    real colvector rho, real scalar max_chunk_rows, real colvector cohort_scope,
    | real colvector keep)
{
    struct cdid_sel_reader scalar rdr
    real colvector U, known_idx, u_idx, s_idx, uk_idx, ak_idx, scoped_U, mapped
    real scalar k, u

    if (!decl.enabled)
        _cdid_sel_fail("input", "Cannot build active selection reader when selection is 'none'.", 498)
    if (!_cdid_sel_is_pos_int(max_chunk_rows, 1))
        _cdid_sel_fail("input", "max_chunk_rows must be a positive integer.", 498)
    if (rows(w) != pop.N | cols(w) != 1 | rows(Z) != pop.N)
        _cdid_sel_fail("input", "Recipe selection design and weight slots must have original N rows.", 498)

    if (decl.selection_universe == "F") U = J(pop.N, 1, 1)
    else if (decl.selection_universe == "B") U = pop.B
    else U = pop.A

    (void) cdid_sel_check_mode(decl, pop.N, U, pop.A, pop.B, rho)
    /* A restricted eligibility context leaves its excluded units out of U. */
    if (args() >= 10) if (rows(keep) == pop.N) U = U :* (keep :!= 0)

    u_idx = _cdid_sel_index(U :== 1)
    if (any(pop.cohort_known[u_idx] :== 0))
        _cdid_sel_fail("support", "Selection estimation universe U contains units with unknown cohort.", 459)
    known_idx = _cdid_sel_index(pop.cohort_known :== 1)
    if (rows(known_idx) == 0)
        _cdid_sel_fail("support", "No units with known cohort in population.", 459)

    rdr.cohorts = uniqrows(sort(pop.cohort[known_idx], 1))
    rdr.n_cohorts = rows(rdr.cohorts)
    if (rows(cohort_scope) != rdr.n_cohorts | cols(cohort_scope) != 1 | missing(cohort_scope) |
        any((cohort_scope :!= 0) :& (cohort_scope :!= 1)) | !any(cohort_scope :== 1))
        _cdid_sel_fail("input", "cohort_scope must be a nonempty binary vector over the known cohorts.", 498)
    rdr.cohort_pos = J(pop.N, 1, 0)
    for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
        uk_idx = _cdid_sel_index((pop.cohort_known :== 1) :& (pop.cohort :== rdr.cohorts[k]))
        if (rows(uk_idx) > 0) rdr.cohort_pos[uk_idx] = J(rows(uk_idx), 1, k)
    }
    if (any(rdr.cohort_pos[u_idx] :== 0))
        _cdid_sel_fail("support", "Active U units could not be mapped to known cohorts.", 459)

    scoped_U = J(pop.N, 1, 0)
    mapped = _cdid_sel_index(rdr.cohort_pos :> 0)
    if (rows(mapped)) scoped_U[mapped] = cohort_scope[rdr.cohort_pos[mapped]] :* U[mapped]
    s_idx = _cdid_sel_index(scoped_U :== 1)
    rdr.spec = cdid_sel_design_classify(Z, scoped_U, has_intercept, int_col)
    if (missing(w[s_idx]) | any(w[s_idx] :<= 0))
        _cdid_sel_fail("support", "Recipe sampling weights must be finite and strictly positive on every used-cohort U row.", 459)

    rdr.mass_F = J(rdr.n_cohorts, 1, .)
    rdr.mass_U = J(rdr.n_cohorts, 1, 0)
    rdr.repeat_mass = J(rdr.n_cohorts, 1, 0)
    rdr.a_U = J(rdr.n_cohorts, 1, 0)
    for (k = 1; k <= rdr.n_cohorts; k = k + 1) {
        if (!cohort_scope[k]) continue
        uk_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1))
        if (rows(uk_idx) > 0) rdr.mass_U[k] = sum(w[uk_idx]) / pop.N
        ak_idx = _cdid_sel_index((rdr.cohort_pos :== k) :& (U :== 1) :& (pop.A :== 1))
        if (rows(ak_idx) > 0) rdr.repeat_mass[k] = sum(w[ak_idx]) / pop.N
        rdr.a_U[k] = (rdr.mass_U[k] > 0 ? rdr.repeat_mass[k] / rdr.mass_U[k] : 0)
        if (decl.selection_universe == "F") rdr.mass_F[k] = rdr.mass_U[k]
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
    return(rdr)
}

real colvector _cdid_sel_link_logistic(real colvector eta)
{
    real colvector p, hi_ix, lo_ix, o
    real scalar n

    if (missing(eta)) _cdid_nuisance_fail("Nonfinite logit predictor.", 498)
    n = rows(eta)
    if (n == 0) return(J(0, 1, 0))
    p = J(n, 1, 0)
    hi_ix = _cdid_sel_link_index(eta :> 40)
    if (rows(hi_ix) == n) return(1 :/ (1 :+ exp(-eta)))
    if (rows(hi_ix) == 0) {
        o = exp(eta)
        return(o :/ (1 :+ o))
    }
    p[hi_ix] = 1 :/ (1 :+ exp(-eta[hi_ix]))
    lo_ix = _cdid_sel_link_index(eta :<= 40)
    o = exp(eta[lo_ix])
    p[lo_ix] = o :/ (1 :+ o)
    return(p)
}

real matrix _cdid_sel_link_logit_vals(real colvector eta)
{
    real colvector p, active, odds, act_ix
    real scalar cap, cap_odds

    p = _cdid_sel_link_logistic(eta)
    active = cdid_ipw_odds_active(p)
    cap = 1 - 1e-6
    cap_odds = cap / (1 - cap)
    odds = J(rows(eta), 1, cap_odds)
    act_ix = _cdid_sel_link_index(active)
    if (rows(act_ix) > 0) odds[act_ix] = exp(eta[act_ix])
    if (missing(odds) | any(odds :<= 0)) {
        _cdid_ipw_fail("Propensity odds must be finite and strictly positive.", 504)
    }
    return((p, odds, active))
}

void _cdid_sel_link_trt_logit(real matrix X, real colvector D, real colvector w,
    real scalar eps, real scalar maxit, real colvector beta_out, real matrix info_out)
{
    real matrix Xa, Q, R, A, info_chol
    real colvector da, wa, ix, p, eta, z, sw, b, bn, qy, tau, signed_x
    real rowvector x_sum
    real scalar n, k, j, h, it, rk, dev, olddev, boundary, converged, qtol

    n = rows(X); k = cols(X)
    if (n == 0 | k == 0 | rows(D) != n | cols(D) != 1 | rows(w) != n | cols(w) != 1 |
        missing(X) | missing(D) | missing(w) | any(w :< 0) | any((D :!= 0) :& (D :!= 1)) |
        eps <= 0 | eps >= . | maxit < 1 | maxit != floor(maxit) | maxit >= .)
        _cdid_nuisance_fail("Require finite design, binary response, nonnegative weights and explicit GLM controls.", 498)
    ix = _cdid_sel_link_index(w :> 0)
    if (rows(ix) < k) _cdid_nuisance_fail("Active logit rows smaller than model dimension.", 459)
    Xa = (rows(ix) == n ? X : X[ix, ])
    da = (rows(ix) == n ? D : D[ix])
    wa = (rows(ix) == n ? w : w[ix])
    x_sum = colsum(Xa)
    for (j = 1; j <= k; j++) {
        if (x_sum[j] == 0) {
            if (all(Xa[, j] :== 0)) _cdid_nuisance_fail("Logit zero predictor; model unchanged.", 459)
        }
        for (h = 1; h < j; h++) {
            if (x_sum[j] == x_sum[h]) {
                if (all(Xa[, j] :== Xa[, h])) _cdid_nuisance_fail("A covariate in a logit nuisance model duplicates another column or the intercept (for example, a constant covariate); remove it.", 459)
            }
        }
        signed_x = (2 :* da :- 1) :* Xa[, j]
        if ((all(signed_x :>= 0) & any(signed_x :> 0)) | (all(signed_x :<= 0) & any(signed_x :< 0)))
            _cdid_nuisance_fail("No finite logit root: signed-coordinate separation witness.", 459)
    }
    b = J(k, 1, 0); eta = J(rows(ix), 1, 0); p = J(rows(ix), 1, .5)
    olddev = _cdid_logit_deviance(da, p, wa); converged = 0; boundary = 0
    qtol = min((1e-7, eps / 1000))
    for (it = 1; it <= maxit; it++) {
        sw = wa :* p :* (1 :- p)
        info_out = cross(Xa, sw, Xa)
        info_chol = cholesky(info_out)
        bn = J(0, 1, .)
        if (!missing(info_chol)) {
            if (all(diagonal(info_chol) :> 100 * qtol :* sqrt(diagonal(info_out)))) {
                z = cross(Xa, sw :* eta :+ wa :* (da :- p))
                qy = solvelower(info_chol, z, rk = ., 0)
                if (rk == k & !missing(qy)) {
                    bn = solveupper(info_chol', qy, rk = ., 0)
                    if (rk != k | missing(bn)) bn = J(0, 1, .)
                }
            }
        }
        if (rows(bn) != k) {
            sw = sqrt(sw); A = Xa :* sw
            hqrd(A, Q = ., tau = ., R = .)
            if (rows(R) != k | missing(R) | any(abs(diagonal(R)) :<= qtol :* sqrt(colsum(A :^ 2))'))
                _cdid_nuisance_fail("Covariates in a logit nuisance model are collinear on its estimation sample; remove redundant covariates.", 459)
            z = (eta :+ (da :- p) :/ (p :* (1 :- p))) :* sw
            qy = hqrdmultq1t(Q, tau, z)
            bn = solveupper(R, qy, rk = ., 0)
            if (rk != k | missing(bn)) _cdid_nuisance_fail("Nonfinite logit IRLS step.", 504)
        }
        eta = Xa * bn; p = _cdid_sel_link_logistic(eta)
        dev = _cdid_logit_deviance(da, p, wa)
        h = 0
        while (dev >= . | any(p :<= 0) | any(p :>= 1)) {
            boundary = 1; h = h + 1
            if (h > maxit) _cdid_nuisance_fail("Logit boundary step-halving failed.", 430)
            bn = (bn + b) / 2; eta = Xa * bn; p = _cdid_sel_link_logistic(eta)
            dev = _cdid_logit_deviance(da, p, wa)
        }
        b = bn
        if (abs(dev - olddev) / (.1 + abs(dev)) < eps) { converged = 1; break; }
        olddev = dev
    }
    if (boundary) _cdid_nuisance_fail("Logit encountered numerical boundary.", 430)
    if (!converged) _cdid_nuisance_fail("Logit did not converge under supplied controls.", 430)
    if (any(abs(eta) :> 30) | any(p :<= 10 * epsilon(1)) | any(p :>= 1 - 10 * epsilon(1)))
        _cdid_nuisance_fail("Logit fitted probabilities reached numerical boundary.", 430)
    if (rows(ix) < n) (void) _cdid_sel_link_logistic(X * b)
    info_out = quadcross(Xa, wa :* p :* (1 :- p), Xa)
    info_chol = cholesky(info_out)
    if (missing(info_chol) | min(diagonal(info_chol)) <= 0)
        _cdid_nuisance_fail("Logit information is not positive definite.", 459)
    R = info_chol'
    if (missing(info_out) | missing(R) | _cdid_nuisance_rcond(R) <= 0)
        _cdid_nuisance_fail("Logit information is singular/nonfinite.", 459)
    beta_out = b
}

struct cdid_selected_rows scalar _cdid_sel_link_slice_rows(
    struct cdid_selected_rows scalar full, real colvector units,
    real scalar N, real scalar want_eval, real scalar is_dr, real scalar p_trt)
{
    struct cdid_selected_rows scalar out
    real colvector u_sub, pos, hit, idx
    real scalar n_u, n_f, lo, hi, p1, p2

    n_u = rows(units)
    n_f = rows(full.unit)
    if (n_u == 0 | n_f == 0) {
        out.loc = out.unit = out.dy = out.q = out.cohort_local = J(0, 1, 0)
        if (want_eval) {
            out.score = J(0, p_trt, 0)
            out.R_ipw = J(0, 4, 0)
            out.R_dr  = (is_dr ? J(0, 4, 0) : J(0, 0, 0))
        }
        else out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
        return(out)
    }
    if (n_u == N) {
        if (units[1] == 1 & units[N] == N) {
            if (all(units :== (1::N))) {
                if (want_eval) return(full)
                out = full
                out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
                return(out)
            }
        }
    }
    lo = units[1]; hi = units[n_u]
    if (hi - lo + 1 == n_u) {
        if (n_u == 1 | all(units :== (lo::hi))) {
            p1 = _cdid_ipw_lower_bound(full.unit, lo)
            p2 = _cdid_ipw_lower_bound(full.unit, hi + 1) - 1
            if (p1 > p2) {
                out.loc = out.unit = out.dy = out.q = out.cohort_local = J(0, 1, 0)
                if (want_eval) {
                    out.score = J(0, p_trt, 0)
                    out.R_ipw = J(0, 4, 0)
                    out.R_dr  = (is_dr ? J(0, 4, 0) : J(0, 0, 0))
                }
                else out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
                return(out)
            }
            u_sub = full.unit[|p1 \ p2|]
            out.unit = u_sub
            out.loc  = u_sub :- (lo - 1)
            out.dy   = full.dy[|p1 \ p2|]
            out.q    = (rows(full.q) ? full.q[|p1 \ p2|] : J(0, 1, 0))
            out.cohort_local = (rows(full.cohort_local) ? full.cohort_local[|p1 \ p2|] : J(0, 1, 0))
            if (want_eval) {
                out.score = (p_trt > 0 ? full.score[|p1, 1 \ p2, p_trt|] : J(p2 - p1 + 1, 0, 0))
                out.R_ipw = full.R_ipw[|p1, 1 \ p2, 4|]
                out.R_dr  = (is_dr ? full.R_dr[|p1, 1 \ p2, 4|] : J(0, 0, 0))
            }
            else out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
            return(out)
        }
    }
    pos = _cdid_sel_link_match_sorted(units, full.unit)
    hit = _cdid_sel_link_index(pos :> 0)
    if (rows(hit) == 0) {
        out.loc = out.unit = out.dy = out.q = out.cohort_local = J(0, 1, 0)
        if (want_eval) {
            out.score = J(0, p_trt, 0)
            out.R_ipw = J(0, 4, 0)
            out.R_dr  = (is_dr ? J(0, 4, 0) : J(0, 0, 0))
        }
        else out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
        return(out)
    }
    idx = pos[hit]
    out.loc  = hit
    out.unit = units[hit]
    out.dy   = full.dy[idx]
    out.q    = (rows(full.q) ? full.q[idx] : J(0, 1, 0))
    out.cohort_local = (rows(full.cohort_local) ? full.cohort_local[idx] : J(0, 1, 0))
    if (want_eval) {
        out.score = full.score[idx, .]
        out.R_ipw = full.R_ipw[idx, .]
        out.R_dr  = (is_dr ? full.R_dr[idx, .] : J(0, 0, 0))
    }
    else out.score = out.R_ipw = out.R_dr = J(0, 0, 0)
    return(out)
}

/* The production model borrows one complete population owner.
   Pair reads return only the requested rows. */
struct cdid_selected_rows scalar _cdid_selected_pair_rows_pop(struct cdid_population scalar pop,
    real scalar s, real scalar t, real colvector units)
{
    struct cdid_selected_rows scalar out
    real colvector rs,rt,hit,keep
    rs=cdid_population_rows(pop,units,J(rows(units),1,s))
    rt=cdid_population_rows(pop,units,J(rows(units),1,t))
    hit=_cdid_sel_link_index((rs:>0):&(rt:>0))
    if(rows(hit)) {
        rs=rs[hit];rt=rt[hit]
        keep=_cdid_sel_link_index((pop.y[rs]:<.):&(pop.y[rt]:<.))
        hit=(rows(keep) ? hit[keep] : J(0,1,0))
        if(rows(keep)){rs=rs[keep];rt=rt[keep];}
    }
    out.loc=hit;out.unit=(rows(hit) ? units[hit] : J(0,1,0))
    out.dy=(rows(hit) ? pop.y[rt]-pop.y[rs] : J(0,1,0))
    out.q=J(0,1,0);out.cohort_local=J(0,1,0)
    out.score=out.R_ipw=out.R_dr=J(0,0,0)
    return(out)
}
struct cdid_selected_rows scalar _cdid_selected_pair_rows(struct cdid_selected_model scalar m,
    real scalar pair, real colvector units)
{
    if (m.pair_full_ready) {
        return(_cdid_sel_link_slice_rows(m.pair_full[pair], units, m.N, 0, 0, 0))
    }
    return(_cdid_selected_pair_rows_pop(*m.population_ref,m.pairs[pair].s_val,m.pairs[pair].t_val,units))
}
struct cdid_selected_rows scalar _cdid_sel_support_from_pair(struct cdid_selected_model scalar m,
    real scalar sid, struct cdid_selected_rows scalar pair_rows)
{
    struct cdid_selected_rows scalar out
    struct cdid_sel_pred scalar prediction
    real colvector take, pos
    real scalar r, k, p1, n_out, n_sup
    out.cohort_local=J(0,1,0);out.score=out.R_ipw=out.R_dr=J(0,0,0)
    if(rows(pair_rows.unit)==0) {
        out.loc=out.unit=out.dy=out.q=J(0,1,0);return(out)
    }
    r=m.supports[sid].recipe_id;k=m.supports[sid].cohort_pos
    take=_cdid_sel_link_index((m.cohort_idx[pair_rows.unit]:==k):&_cdid_sld_keep(m,r,pair_rows.unit))
    if(rows(take)==0) {
        out.loc=out.unit=out.dy=out.q=J(0,1,0);return(out)
    }
    out.loc=pair_rows.loc[take];out.unit=pair_rows.unit[take];out.dy=pair_rows.dy[take]
    if (m.support_full_ready) {
        n_out = rows(out.unit)
        n_sup = rows(m.support_full[sid].unit)
        if (n_out == n_sup) {
            if (out.unit[1] == m.support_full[sid].unit[1] & out.unit[n_out] == m.support_full[sid].unit[n_sup]) {
                if (all(out.unit :== m.support_full[sid].unit)) {
                    out.q = m.support_full[sid].q
                    return(out)
                }
            }
        }
        p1 = _cdid_ipw_lower_bound(m.support_full[sid].unit, out.unit[1])
        if (p1 + n_out - 1 <= n_sup) {
            if (all(m.support_full[sid].unit[|p1 \ (p1 + n_out - 1)|] :== out.unit)) {
                out.q = m.support_full[sid].q[|p1 \ (p1 + n_out - 1)|]
                return(out)
            }
        }
        pos = _cdid_sel_link_match_sorted(out.unit, m.support_full[sid].unit)
        out.q = m.support_full[sid].q[pos]
        return(out)
    }
    if(m.ctl.selection=="conditional_mtar") {
        _cdid_selected_activate_reader(m,r)
        prediction=cdid_sel_pred_pair(m.contexts[r].sel_fit,m.reader,m.supports[sid].pc_id,out.unit)
        out.q=prediction.sel_weight
    } else out.q=_cdid_sld_w_rows(m,r,out.unit)/m.supports[sid].rate_c
    return(out)
}
struct cdid_selected_rows scalar _cdid_selected_support_rows(struct cdid_selected_model scalar m,
    real scalar sid, real colvector units)
{
    if (m.support_full_ready) {
        return(_cdid_sel_link_slice_rows(m.support_full[sid], units, m.N, 0, 0, 0))
    }
    return(_cdid_sel_support_from_pair(m,sid,_cdid_selected_pair_rows(m,m.supports[sid].pair_id,units)))
}
struct cdid_selected_rows scalar _cdid_selected_treatment_rows(struct cdid_selected_model scalar m,
    real scalar tm, real colvector units, real scalar evaluate,
    | struct cdid_selected_rows colvector pair_cache)
{
    struct cdid_selected_rows scalar pair_rows,out
    struct cdid_sel_pred scalar prediction
    real scalar k,r,sid,b,L,p_b,pi,p1,n_sub,n_sup
    real colvector cohort_row,cmap,take,sub,D,prob,odds,active,wp,wm,residual,beta_b,pos
    real rowvector kept_b
    real matrix X,design,logit_values
    if (m.trt_models[tm].risk_cache_ready) {
        return(_cdid_sel_link_slice_rows(m.trt_models[tm].risk_cache, units, m.N,
            evaluate, m.ctl.est_method == "dr", m.trt_models[tm].p_trt))
    }
    r=m.trt_models[tm].recipe_id;pi=m.trt_models[tm].pair_id
    if(args()>=5 & rows(pair_cache)>=pi) pair_rows=pair_cache[pi]
    else pair_rows=_cdid_selected_pair_rows(m,pi,units)
    cohort_row=J(rows(pair_rows.unit),1,0)
    if(rows(pair_rows.unit)) {
        cmap=J(m.n_cohorts+1,1,0)
        cmap[m.trt_models[tm].cohorts:+1]=(1::rows(m.trt_models[tm].cohorts))
        cohort_row=cmap[m.cohort_idx[pair_rows.unit]:+1]:*_cdid_sld_keep(m,r,pair_rows.unit)
    }
    take=_cdid_sel_link_index(cohort_row:>0)
    if(rows(take)==0) {
        out.loc=out.unit=out.dy=out.q=out.cohort_local=J(0,1,0)
        out.score=J(0,m.trt_models[tm].p_trt,0);out.R_ipw=J(0,4,0)
        out.R_dr=(m.ctl.est_method=="dr" ? J(0,4,0) : J(0,0,0))
        return(out)
    }
    out.loc=pair_rows.loc[take];out.unit=pair_rows.unit[take];out.dy=pair_rows.dy[take]
    out.cohort_local=cohort_row[take]
    if (m.R == 1 & rows(pair_rows.q) == rows(pair_rows.unit)) {
        out.q = pair_rows.q[take]
    }
    else {
        out.q=J(rows(out.unit),1,0)
        if(m.ctl.selection=="conditional_mtar" & !m.support_full_ready) _cdid_selected_activate_reader(m,r)
        for(k=1;k<=rows(m.trt_models[tm].cohorts);k++) {
            sub=_cdid_sel_link_index(out.cohort_local:==k)
            if(rows(sub)==0) continue
            sid=m.trt_models[tm].support_ids[k]
            if (m.support_full_ready) {
                n_sub = rows(sub)
                n_sup = rows(m.support_full[sid].unit)
                if (n_sub == n_sup) {
                    if (out.unit[sub[1]] == m.support_full[sid].unit[1] & out.unit[sub[n_sub]] == m.support_full[sid].unit[n_sup]) {
                        if (all(out.unit[sub] :== m.support_full[sid].unit)) {
                            out.q[sub] = m.support_full[sid].q
                            continue
                        }
                    }
                }
                p1 = _cdid_ipw_lower_bound(m.support_full[sid].unit, out.unit[sub[1]])
                if (p1 + n_sub - 1 <= n_sup) {
                    if (all(m.support_full[sid].unit[|p1 \ (p1 + n_sub - 1)|] :== out.unit[sub])) {
                        out.q[sub] = m.support_full[sid].q[|p1 \ (p1 + n_sub - 1)|]
                        continue
                    }
                }
                pos = _cdid_sel_link_match_sorted(out.unit[sub], m.support_full[sid].unit)
                out.q[sub] = m.support_full[sid].q[pos]
            }
            else if(m.ctl.selection=="conditional_mtar") {
                prediction=cdid_sel_pred_pair(m.contexts[r].sel_fit,m.reader,m.supports[sid].pc_id,out.unit[sub])
                out.q[sub]=prediction.sel_weight
            } else out.q[sub]=_cdid_sld_w_rows(m,r,out.unit[sub])/m.supports[sid].rate_c
        }
    }
    if(!evaluate) return(out)
    if(rows(out.unit)==0) {
        out.score=J(0,m.trt_models[tm].p_trt,0);out.R_ipw=J(0,4,0)
        out.R_dr=(m.ctl.est_method=="dr" ? J(0,4,0) : J(0,0,0));return(out)
    }
    X=_cdid_sld_trt_rows(m,r,out.unit);D=(out.cohort_local:==1)
    logit_values=_cdid_sel_link_logit_vals(X*m.trt_models[tm].beta_trt)
    prob=logit_values[,1];odds=logit_values[,2];active=logit_values[,3]
    wp=out.q:*D;wm=out.q:*(1:-D):*odds
    out.score=X:*(out.q:*(D:-prob))
    out.R_ipw=(wp:*(out.dy:-m.trt_models[tm].T_ipw)/m.trt_models[tm].dp,
        wm:*(out.dy:-m.trt_models[tm].M_ipw)/m.trt_models[tm].dm,wp,wm)
    if(m.ctl.est_method=="dr") {
        b=m.trt_models[tm].block_id;L=m.or_blocks[b].L
        p_b=(m.outcome_pooling=="separate" ? m.or_blocks[b].p_kept : m.or_shared.p_kept)
        kept_b=(m.outcome_pooling=="separate" ? m.or_blocks[b].kept : m.or_shared.kept)
        beta_b=(m.outcome_pooling=="separate" ? m.or_blocks[b].beta : m.or_shared.beta)
        design=(p_b ? L*_cdid_sld_or_rows(m,r,out.unit,kept_b) : J(rows(out.unit),0,0))
        residual=out.dy:-m.or_blocks[b].alpha:-(p_b ? design*beta_b : J(rows(out.unit),1,0))
        out.R_dr=(wp:*(residual:-m.trt_models[tm].T_dr)/m.trt_models[tm].dp,
            wm:*(residual:-m.trt_models[tm].M_dr)/m.trt_models[tm].dm,wp,wm)
    }
    return(out)
}

/* Verify exact alias reconstruction on actual control and required target
   prediction support rows for every block b = 1..n_blocks. */
void _cdid_sel_link_verify_aliases(struct cdid_selected_model scalar m)
{
    struct cdid_selected_rows scalar read_rows
    real scalar n_drop, b, s, r, L_b, j, den, err, threshold
    real colvector used, obs_u, v, scale
    real rowvector col_max, drop_b, kept_b
    real matrix Z_s, map_b
    real scalar cm

    threshold = 64 * epsilon(1) * max((1, m.or_shared.p_raw + m.or_shared.intercept))
    if (m.outcome_pooling == "separate") {
        for (b = 1; b <= m.n_blocks; b = b + 1) {
            drop_b = m.or_blocks[b].dropped
            n_drop = cols(drop_b)
            if (n_drop == 0) continue
            kept_b = m.or_blocks[b].kept
            map_b  = m.or_blocks[b].alias_map
            used   = J(m.n_supports, 1, 0)
            used[m.or_blocks[b].support_ids] = J(rows(m.or_blocks[b].support_ids), 1, 1)
            for (j = 1; j <= m.n_links; j = j + 1) {
                if (m.links[j].block_id == b) used[m.links[j].target_support_id] = 1
            }
            col_max = J(1, n_drop, 0)
            for (s = 1; s <= m.n_supports; s = s + 1) {
                if (!used[s]) continue
                r = m.supports[s].recipe_id
                read_rows = _cdid_selected_support_rows(m, s, (1::m.N))
                obs_u = read_rows.unit
                if (rows(obs_u) > 0) {
                    col_max = colmax(col_max \ abs(_cdid_sld_or_rows(m, r, obs_u, drop_b)))
                }
            }
            L_b = m.or_blocks[b].L
            for (s = 1; s <= m.n_supports; s = s + 1) {
                if (!used[s]) continue
                r = m.supports[s].recipe_id
                read_rows = _cdid_selected_support_rows(m, s, (1::m.N))
                obs_u = read_rows.unit
                if (rows(obs_u) == 0) continue
                Z_s = _cdid_sld_or_rows(m, r, obs_u)
                for (j = 1; j <= n_drop; j = j + 1) {
                    cm = (m.or_shared.intercept ? m.or_blocks[b].alias_constants[1, j] / L_b : 0)
                    if (cols(kept_b) > 0) {
                        v = Z_s[, drop_b[j]] - Z_s[, kept_b] * map_b[, j] :- cm
                        scale = abs(Z_s[, drop_b[j]]) +
                            abs(Z_s[, kept_b]) * abs(map_b[, j]) :+ abs(cm)
                    }
                    else {
                        v = Z_s[, drop_b[j]] :- cm
                        scale = abs(Z_s[, drop_b[j]]) :+ abs(cm)
                    }
                    den = max((col_max[j], 2.2250738585072014e-308))
                    err = max(abs(v) :/ rowmax((scale, J(rows(v), 1, den))))
                    if (err > threshold | err >= .) {
                        _cdid_sel_link_fail("support",
                            "Pooled control design has an exact null direction that does not reconstruct required treated predictions.",
                            459)
                    }
                }
            }
            m.or_blocks[b].verif_scan = 1
        }
        return
    }

    n_drop = cols(m.or_shared.dropped)
    if (n_drop == 0) return

    col_max = J(1, n_drop, 0)
    for (s = 1; s <= m.n_supports; s = s + 1) {
        r = m.supports[s].recipe_id
        read_rows=_cdid_selected_support_rows(m,s,(1::m.N))
        obs_u=read_rows.unit
        if (rows(obs_u) > 0) {
            col_max = colmax(col_max \ abs(_cdid_sld_or_rows(m,r,obs_u,m.or_shared.dropped)))
        }
    }
    for (b = 1; b <= m.n_blocks; b = b + 1) {
        used = J(m.n_supports, 1, 0)
        used[m.or_blocks[b].support_ids] = J(rows(m.or_blocks[b].support_ids), 1, 1)
        for (j = 1; j <= m.n_links; j = j + 1) {
            if (m.links[j].block_id == b) used[m.links[j].target_support_id] = 1
        }
        L_b = m.or_blocks[b].L
        for (s = 1; s <= m.n_supports; s = s + 1) {
            if (!used[s]) continue
            r = m.supports[s].recipe_id
            read_rows=_cdid_selected_support_rows(m,s,(1::m.N))
        obs_u=read_rows.unit
            if (rows(obs_u) == 0) continue
            Z_s = _cdid_sld_or_rows(m,r,obs_u)
            for (j = 1; j <= n_drop; j = j + 1) {
                cm = (m.or_shared.intercept ? m.or_shared.alias_constants[b, j] / L_b : 0)
                if (cols(m.or_shared.kept) > 0) {
                    v = Z_s[, m.or_shared.dropped[j]] -
                        Z_s[, m.or_shared.kept] * m.or_shared.alias_map[, j] :- cm
                    scale = abs(Z_s[, m.or_shared.dropped[j]]) +
                        abs(Z_s[, m.or_shared.kept]) * abs(m.or_shared.alias_map[, j]) :+ abs(cm)
                }
                else {
                    v = Z_s[, m.or_shared.dropped[j]] :- cm
                    scale = abs(Z_s[, m.or_shared.dropped[j]]) :+ abs(cm)
                }
                den = max((col_max[j], 2.2250738585072014e-308))
                err = max(abs(v) :/ rowmax((scale, J(rows(v), 1, den))))
                if (err > threshold | err >= .) {
                    _cdid_sel_link_fail("support",
                        "Pooled control design has an exact null direction that does not reconstruct required treated predictions.",
                        459)
                }
            }
        }
    }
    m.or_shared.verification_scan = 1
}

/* Store only structurally used coordinates and moment columns.  No global
   parameter-by-all-moments matrix is retained. Derivatives already contain /N. */
void _cdid_selected_add_loading(struct cdid_selected_model scalar m, real scalar r,
    real scalar pc, real colvector units, real colvector columns, real matrix derivative)
{
    struct cdid_sel_log_deriv scalar ld
    struct cdid_selected_loading_piece scalar piece
    if(rows(units)==0 | rows(columns)==0) return
    if(rows(derivative)!=rows(units) | cols(derivative)!=rows(columns))
        _cdid_sel_link_fail("input", "Sparse loading derivative dimensions disagree.", 498)
    _cdid_selected_activate_reader(m,r)
    ld=cdid_sel_log_deriv_eval(m.contexts[r].sel_fit,m.reader,m.contexts[r].up_map,pc,units)
    piece.theta_ids=_cdid_sel_pc_global_idx(m.contexts[r].sel_fit,m.contexts[r].up_map,pc)
    piece.columns=columns; piece.cohort=ld.cohort_pos
    piece.G_theta=quadcross(ld.d_theta_local,derivative)
    piece.G_ak=quadcross(ld.d_ak,derivative)
    m.contexts[r].n_loading_pieces=m.contexts[r].n_loading_pieces+1
    if(m.contexts[r].n_loading_pieces>rows(m.contexts[r].loading_pieces))
        _cdid_sel_link_fail("integrity","Sparse selection loading capacity was underestimated.",498)
    m.contexts[r].loading_pieces[m.contexts[r].n_loading_pieces]=piece
}
void _cdid_selected_add_loading_sup(struct cdid_selected_model scalar m, real scalar sid,
    real colvector units, real colvector columns, real matrix derivative)
{
    struct cdid_selected_loading_piece scalar piece
    real scalar r, pc
    if (rows(units) == 0 | rows(columns) == 0) return
    r = m.supports[sid].recipe_id
    pc = m.supports[sid].pc_id
    if (m.supports[sid].ld_ready & rows(m.supports[sid].ld.d_ak) == rows(units)) {
        if (rows(derivative) != rows(units) | cols(derivative) != rows(columns))
            _cdid_sel_link_fail("input", "Sparse loading derivative dimensions disagree.", 498)
        piece.theta_ids = m.supports[sid].ld_g_idx
        piece.columns   = columns
        piece.cohort    = m.supports[sid].ld.cohort_pos
        piece.G_theta   = quadcross(m.supports[sid].ld.d_theta_local, derivative)
        piece.G_ak      = quadcross(m.supports[sid].ld.d_ak, derivative)
        m.contexts[r].n_loading_pieces = m.contexts[r].n_loading_pieces + 1
        if (m.contexts[r].n_loading_pieces > rows(m.contexts[r].loading_pieces))
            _cdid_sel_link_fail("integrity", "Sparse selection loading capacity was underestimated.", 498)
        m.contexts[r].loading_pieces[m.contexts[r].n_loading_pieces] = piece
        return
    }
    _cdid_selected_add_loading(m, r, pc, units, columns, derivative)
}
/* Inflate a accumulator from loading_pieces for the requested
   column batch, then contract shared model scores once for that batch. */
real matrix _cdid_selected_loading_pieces(struct cdid_selected_model scalar m, real scalar r,
    real colvector units, real colvector columns)
{
    struct cdid_sel_moments scalar acc
    struct cdid_selected_loading_piece scalar piece
    real colvector position, take, to
    real scalar j
    if(rows(columns)==0) return(J(rows(units),0,0))
    _cdid_selected_activate_reader(m,r)
    acc=cdid_sel_moments_init(m.contexts[r].sel_fit,rows(columns))
    position=J(m.contexts[r].n_moments_total,1,0);position[columns]=(1::rows(columns))
    for(j=1;j<=m.contexts[r].n_loading_pieces;j++) {
        piece=m.contexts[r].loading_pieces[j]
        take=_cdid_sel_link_index(position[piece.columns]:>0)
        if(rows(take)==0) continue
        to=position[piece.columns[take]]
        if(rows(piece.theta_ids))
            acc.G_theta[piece.theta_ids,to']=acc.G_theta[piece.theta_ids,to']+piece.G_theta[,take']
        acc.G_ak[piece.cohort,to']=acc.G_ak[piece.cohort,to']+piece.G_ak[,take']
    }
    return(cdid_sel_eval_loading(m.contexts[r].sel_fit,m.reader,acc,units))
}

real matrix _cdid_selected_loading(struct cdid_selected_model scalar m, real scalar r,
    real colvector units, real colvector columns)
{
    struct cdid_sel_moments scalar acc
    if(rows(columns)==0) return(J(rows(units),0,0))
    if (m.contexts[r].acc_full_ready) {
        _cdid_selected_activate_reader(m,r)
        acc.n_moments = rows(columns)
        acc.G_theta   = m.contexts[r].acc_full.G_theta[, columns']
        acc.G_ak      = m.contexts[r].acc_full.G_ak[, columns']
        return(cdid_sel_eval_loading(m.contexts[r].sel_fit,m.reader,acc,units))
    }
    return(_cdid_selected_loading_pieces(m, r, units, columns))
}

/* Build one bounded requested-family view. Global beta/objective forcing is
   summed across every recipe without retaining a unit-by-recipe output grid. */
struct cdid_selected_loaded_block colvector _cdid_selected_load_requests(
    struct cdid_selected_model scalar m, real colvector units, real colvector link_ids,
    real scalar details, real matrix beta_force, real colvector objective_force)
{
    struct cdid_selected_loaded_block colvector loaded
    real colvector ids,prefix,all_ids,trt_ids,blocks,beta_blocks,blk_mark,bblk_mark,trt_mark
    real matrix value
    real scalar r,j,l,b,tm,at,p,p_b,need_or,need_obj,is_sep,use_fast_acc,n_ids,pos_id
    loaded=cdid_selected_loaded_block(m.R,1)
    is_sep=(m.outcome_pooling=="separate")
    p=(is_sep ? 0 : m.or_shared.p_kept)
    beta_force=J(rows(units),p,0);objective_force=J(rows(units),1,0)
    need_or=(m.ctl.est_method!="ipw") & (m.ctl.est_method=="or" | details)
    need_obj=(m.ctl.est_method!="ipw") & details
    for(r=1;r<=m.R;r++) {
        use_fast_acc = (!details & m.contexts[r].acc_full_ready)
        prefix=J(0,1,0)
        if(!is_sep & m.ctl.est_method!="ipw" & p)
            prefix=(m.contexts[r].m_col_beta_start::m.contexts[r].m_col_beta_end)
        if(need_obj) prefix=prefix\m.contexts[r].m_col_obj
        blk_mark  = J(m.n_blocks, 1, 0)
        bblk_mark = J(m.n_blocks, 1, 0)
        trt_mark  = (m.n_trt_models > 0 ? J(m.n_trt_models, 1, 0) : J(0, 1, 0))
        n_ids = 0
        for(j=1;j<=rows(link_ids);j++) {
            l=link_ids[j]
            if(is_sep & m.ctl.est_method!="ipw") {
                b=m.links[l].block_id
                if(m.or_blocks[b].recipe_id==r) bblk_mark[b]=1
            }
            if(need_or & m.links[l].recipe_id==r) {
                b=m.links[l].block_id;blk_mark[b]=1
                if (details) {
                    p_b=(is_sep ? m.or_blocks[b].p_kept : p)
                    n_ids = n_ids + 2 + p_b
                }
                else n_ids = n_ids + 1
            }
            if(m.ctl.est_method!="or") {
                tm=m.links[l].trt_id
                if(m.trt_models[tm].recipe_id==r) trt_mark[tm]=1
            }
        }
        beta_blocks = _cdid_sel_link_index(bblk_mark)
        blocks      = _cdid_sel_link_index(blk_mark)
        trt_ids     = _cdid_sel_link_index(trt_mark)
        for(j=1;j<=rows(beta_blocks);j++) n_ids = n_ids + m.or_blocks[beta_blocks[j]].p_kept
        if(m.or_shared.intercept) n_ids = n_ids + rows(blocks)
        for(j=1;j<=rows(trt_ids);j++) {
            tm=trt_ids[j]
            if (!details & use_fast_acc & m.trt_models[tm].m_col_fast > 0) n_ids = n_ids + 1
            else if (!details) n_ids = n_ids + m.trt_models[tm].p_trt + 2
            else n_ids = n_ids + m.trt_models[tm].p_trt + (m.ctl.est_method=="dr" ? 8 : 4)
        }
        ids = J(n_ids, 1, 0)
        pos_id = 0
        for(j=1;j<=rows(link_ids);j++) {
            l=link_ids[j]
            if(need_or & m.links[l].recipe_id==r) {
                if (details) {
                    b=m.links[l].block_id
                    p_b=(is_sep ? m.or_blocks[b].p_kept : p)
                    ids[pos_id + 1] = m.links[l].m_col_or_mu
                    ids[pos_id + 2] = m.links[l].m_col_or_den
                    pos_id = pos_id + 2
                    if(p_b) {
                        ids[|(pos_id + 1) \ (pos_id + p_b)|] =
                            (m.links[l].m_col_or_des_start::m.links[l].m_col_or_des_end)
                        pos_id = pos_id + p_b
                    }
                }
                else {
                    pos_id = pos_id + 1
                    ids[pos_id] = m.links[l].m_col_or_mu
                }
            }
        }
        for(j=1;j<=rows(beta_blocks);j++) {
            b=beta_blocks[j]
            p_b=m.or_blocks[b].p_kept
            if(p_b>0) {
                ids[|(pos_id + 1) \ (pos_id + p_b)|] =
                    (m.or_blocks[b].m_col_beta_start::m.or_blocks[b].m_col_beta_end)
                pos_id = pos_id + p_b
            }
        }
        if(rows(blocks) & m.or_shared.intercept) {
            for(j=1;j<=rows(blocks);j++) {
                pos_id = pos_id + 1
                ids[pos_id] = m.or_blocks[blocks[j]].m_col_alpha
            }
        }
        for(j=1;j<=rows(trt_ids);j++) {
            tm=trt_ids[j]
            if (!details & use_fast_acc & m.trt_models[tm].m_col_fast > 0) {
                pos_id = pos_id + 1
                ids[pos_id] = m.trt_models[tm].m_col_fast
            }
            else if (!details) {
                p_b = m.trt_models[tm].p_trt
                ids[|(pos_id + 1) \ (pos_id + p_b)|] =
                    (m.trt_models[tm].m_col_C_start::m.trt_models[tm].m_col_C_end)
                pos_id = pos_id + p_b
                if (m.ctl.est_method == "dr") {
                    ids[pos_id + 1] = m.trt_models[tm].m_col_U_dr_start
                    ids[pos_id + 2] = m.trt_models[tm].m_col_U_dr_start + 1
                }
                else {
                    ids[pos_id + 1] = m.trt_models[tm].m_col_U_ipw_start
                    ids[pos_id + 2] = m.trt_models[tm].m_col_U_ipw_start + 1
                }
                pos_id = pos_id + 2
            }
            else {
                p_b = m.trt_models[tm].p_trt + 4
                ids[|(pos_id + 1) \ (pos_id + p_b)|] =
                    (m.trt_models[tm].m_col_C_start::m.trt_models[tm].m_col_U_ipw_end)
                pos_id = pos_id + p_b
                if(m.ctl.est_method=="dr") {
                    ids[|(pos_id + 1) \ (pos_id + 4)|] =
                        (m.trt_models[tm].m_col_U_dr_start::m.trt_models[tm].m_col_U_dr_end)
                    pos_id = pos_id + 4
                }
            }
        }
        if(rows(ids)) ids=uniqrows(sort(ids,1))
        all_ids=prefix\ids
        value = (details ? _cdid_selected_loading_pieces(m,r,units,all_ids) :
                           _cdid_selected_loading(m,r,units,all_ids))
        at=0
        if(!is_sep & m.ctl.est_method!="ipw" & p) { beta_force=beta_force+value[,1..p];at=p; }
        if(need_obj) { objective_force=objective_force+value[,at+1];at=at+1; }
        loaded[r].ids=ids
        loaded[r].col_pos = J((m.contexts[r].acc_full_ready ? cols(m.contexts[r].acc_full.G_theta) : m.contexts[r].n_moments_total), 1, 0)
        if (rows(ids)) loaded[r].col_pos[ids] = (1::rows(ids))
        loaded[r].values=(rows(ids) ? value[,(at+1)..cols(value)] : J(rows(units),0,0))
    }
    return(loaded)
}
real matrix _cdid_selected_loaded(struct cdid_selected_loaded_block colvector loaded,
    real scalar recipe, real colvector columns)
{
    real colvector positions
    if (rows(columns) == 0) return(J(rows(loaded[recipe].values), 0, 0))
    if (rows(loaded[recipe].col_pos) >= max(columns) & min(columns) >= 1) {
        positions = loaded[recipe].col_pos[columns]
    }
    else {
        positions=_cdid_sel_link_match_sorted(columns,loaded[recipe].ids)
    }
    if(any(positions:==0)) _cdid_sel_link_fail("replay","Requested loading was not prepared for this chunk family.",498)
    return(loaded[recipe].values[,positions'])
}

/* --------------------------------------------------------------------------
   3. Core Multi-Recipe Preparation Engine
   -------------------------------------------------------------------------- */

void _cdid_sel_fit_or_separate(struct cdid_selected_model scalar m)
{
    struct cdid_selected_rows scalar pair_rows, read_rows
    struct cdid_wls_result scalar wls_fit
    real scalar N, p_or, icpt, b, k, l, sid, r, L_b, mass_b, nn, scale_b, p_b
    real scalar total_p, total_des, sum_a, dres, can_cache_sep, tot_obs_b, at_b
    real colvector s_all, pseudo_y, pseudo_w, pseudo_block, u_k, dy_k, q_k, raw_res, beta_b
    real colvector u_all_b, ha_all_b, ord_b
    real rowvector constant_cols, rawZ_first, coef_rows, tau_q, sum_b, kept_b, c_l
    real matrix F_b, Q_mat, R_mat, V_mat, res_fac, A_ps, V_ps, pseudo_X, D_kept, centered, V_b, hb_all_b

    N = m.N
    p_or = m.or_shared.p_raw
    icpt = m.ctl.or_intercept
    m.n_or_block_fits = m.n_blocks
    m.or_shared.rss = 0
    total_p = 0
    total_des = 0

    for (b = 1; b <= m.n_blocks; b = b + 1) {
        L_b = m.or_blocks[b].L
        F_b = J(0, p_or + 2, 0)
        mass_b = 0
        constant_cols = J(1, p_or, 1)
        rawZ_first = J(1, p_or, 0)
        s_all = m.or_blocks[b].support_ids
        for (k = 1; k <= rows(s_all); k = k + 1) {
            sid = s_all[k]
            F_b = _cdid_or_compress(F_b \ m.supports[sid].F_qr)
            mass_b = mass_b + m.supports[sid].weight_mass
            if (k == 1 & p_or > 0) rawZ_first = m.supports[sid].first_z
            if (p_or > 0) {
                constant_cols = constant_cols :*
                    (m.supports[sid].xmin :== rawZ_first) :*
                    (m.supports[sid].xmax :== rawZ_first)
            }
        }
        hqrd(F_b[, 1], Q_mat = ., tau_q = ., R_mat = .)
        V_mat = hqrdmultq(Q_mat, tau_q, F_b[, (2..(p_or + 2))], 1)
        coef_rows = V_mat[1, ] / R_mat[1, 1]
        res_fac   = (rows(V_mat) > 1 ? V_mat[2..rows(V_mat), ] : J(1, p_or + 1, 0))
        if (icpt & p_or > 0) {
            for (k = 1; k <= p_or; k = k + 1) {
                if (constant_cols[k]) {
                    coef_rows[k] = L_b * rawZ_first[k]
                    res_fac[, k] = J(rows(res_fac), 1, 0)
                }
            }
        }
        nn = rows(res_fac)
        scale_b = sqrt(nn / mass_b)
        A_ps = coef_rows :+ scale_b :* res_fac
        V_ps = coef_rows :- scale_b :* res_fac
        pseudo_X = (p_or > 0 ? (A_ps[, 1..p_or] \ V_ps[, 1..p_or]) : J(2 * nn, 0, 0))
        pseudo_y = A_ps[, p_or + 1] \ V_ps[, p_or + 1]
        pseudo_w = J(2 * nn, 1, mass_b / (2 * nn))
        pseudo_block = J(2 * nn, 1, 1)

        wls_fit = cdid_pooled_wls_fit(pseudo_X, pseudo_y, pseudo_w, pseudo_block,
            J(2 * nn, 1, 1), 1, icpt, J(0, p_or, 0), J(0, 1, 0), m.ctl.rank_tol)
        p_b = rows(wls_fit.beta)
        m.or_blocks[b].p_kept           = p_b
        m.or_blocks[b].beta             = wls_fit.beta
        m.or_blocks[b].R                = wls_fit.R
        m.or_blocks[b].pivot            = wls_fit.pivot
        m.or_blocks[b].kept             = wls_fit.kept
        m.or_blocks[b].dropped          = wls_fit.dropped
        m.or_blocks[b].alias_map        = wls_fit.alias_map
        m.or_blocks[b].alias_constants  = wls_fit.alias_constants
        m.or_blocks[b].alpha            = (icpt ? wls_fit.alpha[1] : 0)
        m.or_blocks[b].weight_sum       = wls_fit.mass[1]
        m.or_blocks[b].mass             = wls_fit.mass[1] / N
        m.or_blocks[b].mean_y           = (icpt ? wls_fit.mean_y[1] : 0)
        m.or_blocks[b].mean_d           = (icpt & p_b > 0 ?
            wls_fit.mean_X[1, wls_fit.kept] : J(1, p_b, 0))
        m.or_blocks[b].verif_scan       = 0
        m.or_blocks[b].m_col_alpha      = 0
        m.or_blocks[b].m_col_beta_start = 0
        m.or_blocks[b].m_col_beta_end   = 0
        m.or_blocks[b].beta_col_start   = total_p + 1
        m.or_blocks[b].beta_col_end     = total_p + p_b
        total_p = total_p + p_b
    }
    m.total_p_sep = total_p

    _cdid_sel_link_verify_aliases(m)

    if (!m.ctl.points_only & m.ctl.selection == "none") {
        m.or_shared.baseline_objective_rate_loading = J(m.n_supports, 1, 0)
    }

    can_cache_sep = (!m.ctl.points_only &
        (m.ctl.max_workspace_doubles >= . |
         m.workspace_doubles + m.n_blocks * m.n_A * (p_or + 2) <= m.ctl.max_workspace_doubles))

    for (b = 1; b <= m.n_blocks; b = b + 1) {
        L_b   = m.or_blocks[b].L
        r     = m.or_blocks[b].recipe_id
        s_all = m.or_blocks[b].support_ids
        p_b   = m.or_blocks[b].p_kept
        kept_b = m.or_blocks[b].kept
        beta_b = m.or_blocks[b].beta
        m.or_blocks[b].rss = 0
        m.or_blocks[b].blk_cache_ready = 0
        sum_a = 0
        sum_b = J(1, p_b, 0)
        if (!m.ctl.points_only & m.ctl.selection == "none") {
            m.or_blocks[b].baseline_alpha_rate = J(rows(s_all), 1, 0)
            V_b = J(p_b, rows(s_all), 0)
        }
        if (can_cache_sep) {
            tot_obs_b = 0
            for (k = 1; k <= rows(s_all); k = k + 1) {
                tot_obs_b = tot_obs_b + m.supports[s_all[k]].n_obs
            }
            u_all_b  = J(tot_obs_b, 1, 0)
            ha_all_b = (icpt ? J(tot_obs_b, 1, 0) : J(0, 1, 0))
            hb_all_b = (p_b > 0 ? J(tot_obs_b, p_b, 0) : J(0, 0, 0))
            at_b = 0
        }
        if (!m.support_full_ready) {
            if (b == 1) pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
            else if (m.or_blocks[b].pair_id != m.or_blocks[b - 1].pair_id) {
                pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
            }
        }
        for (k = 1; k <= rows(s_all); k = k + 1) {
            sid = s_all[k]
            if (m.support_full_ready) {
                u_k = m.support_full[sid].unit; dy_k = m.support_full[sid].dy; q_k = m.support_full[sid].q
            }
            else {
                read_rows = _cdid_sel_support_from_pair(m, sid, pair_rows)
                u_k = read_rows.unit; dy_k = read_rows.dy; q_k = read_rows.q
            }
            if (p_b > 0) {
                D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                    L_b :* m.supports[sid].Z_or[, kept_b] :
                    L_b :* _cdid_sld_or_rows(m, r, u_k, kept_b))
            }
            else D_kept = J(rows(u_k), 0, 0)
            raw_res = dy_k :- m.or_blocks[b].alpha :-
                (p_b > 0 ? D_kept * beta_b : J(rows(u_k), 1, 0))
            m.or_blocks[b].rss = m.or_blocks[b].rss + sum(q_k :* (raw_res :^ 2))
            if (!m.ctl.points_only) {
                centered = (icpt & p_b > 0 ? D_kept :- m.or_blocks[b].mean_d : D_kept)
                if (icpt) sum_a = sum_a + sum(q_k :* raw_res)
                if (p_b > 0) sum_b = sum_b + quadcross(q_k :* raw_res, centered)
                if (can_cache_sep & rows(u_k) > 0) {
                    u_all_b[|(at_b + 1) \ (at_b + rows(u_k))|] = u_k
                    if (icpt) ha_all_b[|(at_b + 1) \ (at_b + rows(u_k))|] = q_k :* raw_res
                    if (p_b > 0) hb_all_b[|(at_b + 1), 1 \ (at_b + rows(u_k)), p_b|] = centered :* (q_k :* raw_res)
                    at_b = at_b + rows(u_k)
                }
                if (m.ctl.selection == "none" & m.supports[sid].rate_active) {
                    m.or_shared.baseline_objective_rate_loading[sid] =
                        m.or_shared.baseline_objective_rate_loading[sid] -
                        sum(q_k :* (raw_res :^ 2)) / (N * m.supports[sid].rate_c)
                    if (icpt) {
                        m.or_blocks[b].baseline_alpha_rate[k] =
                            -sum(q_k :* raw_res) / m.supports[sid].rate_c
                    }
                    if (p_b > 0) {
                        V_b[, k] = V_b[, k] -
                            (quadcross(q_k :* raw_res, centered))' / m.supports[sid].rate_c
                    }
                }
            }
        }
        m.or_blocks[b].objective        = m.or_blocks[b].rss / N
        m.or_blocks[b].score_mean_alpha = sum_a / N
        m.or_blocks[b].score_mean_beta  = sum_b / N
        m.or_shared.rss = m.or_shared.rss + m.or_blocks[b].rss
        if (can_cache_sep) {
            if (tot_obs_b > 1) {
                ord_b = order(u_all_b, 1)
                m.or_blocks[b].blk_obs_u = u_all_b[ord_b]
                m.or_blocks[b].blk_ha    = (icpt ? ha_all_b[ord_b] : J(0, 1, 0))
                m.or_blocks[b].blk_hb    = (p_b > 0 ? hb_all_b[ord_b, .] : J(0, 0, 0))
            }
            else {
                m.or_blocks[b].blk_obs_u = u_all_b
                m.or_blocks[b].blk_ha    = ha_all_b
                m.or_blocks[b].blk_hb    = hb_all_b
            }
            m.or_blocks[b].blk_gram_inv = (p_b > 0 ?
                (N :* _cdid_nuisance_gram(m.or_blocks[b].R, m.or_blocks[b].pivot, I(p_b)))' : J(0, 0, 0))
            m.or_blocks[b].blk_cache_ready = 1
        }
        if (!m.ctl.points_only & m.ctl.selection == "none") {
            m.or_blocks[b].baseline_beta_rate =
                (p_b > 0 ? _cdid_nuisance_gram(m.or_blocks[b].R, m.or_blocks[b].pivot, V_b) : J(0, rows(s_all), 0))
        }
    }
    m.or_shared.objective = m.or_shared.rss / N

    for (l = 1; l <= m.n_links; l = l + 1) {
        b   = m.links[l].block_id
        sid = m.links[l].target_support_id
        L_b = m.links[l].L
        mass_b = m.supports[sid].weight_mass
        p_b    = m.or_blocks[b].p_kept
        kept_b = m.or_blocks[b].kept
        beta_b = m.or_blocks[b].beta
        m.links[l].des_col_start = total_des + 1
        m.links[l].des_col_end   = total_des + p_b
        total_des = total_des + p_b
        c_l = (p_b > 0 ? L_b :* m.supports[sid].sum_z[kept_b] :/ mass_b : J(1, 0, 0))
        dres = m.supports[sid].sum_y -
            (icpt ? m.or_blocks[b].alpha * mass_b : 0) -
            (p_b > 0 ? L_b * (m.supports[sid].sum_z[kept_b] * beta_b) : 0)
        m.links[l].mu_or      = dres / mass_b
        m.links[l].c_target   = c_l
        m.links[l].d_target   = mass_b / N
        m.links[l].num_target = m.links[l].mu_or * m.links[l].d_target
        m.links[l].baseline_or_rate_loading = 0
        m.links[l].baseline_den_rate_loading = 0
        m.links[l].baseline_design_rate_loading = J(1, p_b, 0)
        m.links[l].target_mu_ready = 0
        if (m.ctl.selection == "none" & m.supports[sid].rate_active) {
            m.links[l].baseline_den_rate_loading = -m.links[l].d_target / m.supports[sid].rate_c
            if (p_b) m.links[l].baseline_design_rate_loading =
                -(L_b * m.supports[sid].sum_z[kept_b] - mass_b * c_l) /
                (N * m.supports[sid].rate_c)
            m.links[l].baseline_or_rate_loading =
                -(dres - m.links[l].mu_or * mass_b) / (N * m.supports[sid].rate_c * m.links[l].d_target)
        }
        if (can_cache_sep & m.ctl.est_method == "or" & m.support_full_ready) {
            r = m.links[l].recipe_id
            u_k = m.support_full[sid].unit
            dy_k = m.support_full[sid].dy
            q_k = m.support_full[sid].q
            if (p_b > 0) {
                D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                    L_b :* m.supports[sid].Z_or[, kept_b] :
                    L_b :* _cdid_sld_or_rows(m, r, u_k, kept_b))
            }
            else D_kept = J(rows(u_k), 0, 0)
            raw_res = dy_k :- (icpt ? m.or_blocks[b].alpha : 0) :-
                (p_b > 0 ? D_kept * beta_b : J(rows(u_k), 1, 0))
            m.links[l].target_mu_vec = q_k :* (raw_res :- m.links[l].mu_or) / m.links[l].d_target
            m.links[l].target_mu_ready = 1
        }
    }
    m.total_des_sep = total_des
}

struct cdid_selected_model scalar cdid_selected_prepare_registry(
    struct cdid_selected_recipe_spec colvector recipes_in,
    struct cdid_selected_registry scalar data_in,
    struct cdid_selected_pair_in colvector pairs_in,
    real colvector cohort_idx,
    real colvector cohort_values,
    real colvector calendar_values,
    real colvector is_A,
    real colvector w_pop,
    struct cdid_selected_link_spec colvector link_specs,
    struct cdid_selected_controls scalar raw_ctl,
    pointer() scalar population_ref,
    | real scalar used_scope, real matrix recipe_keep)
{
    struct cdid_selected_model scalar m
    struct cdid_population scalar population
    struct cdid_selected_rows scalar pair_rows, read_rows
    struct cdid_selected_controls scalar ctl
    struct cdid_sel_decl scalar decl
    struct cdid_sel_reader scalar temporary_reader
    struct cdid_sel_pred scalar pr
    struct cdid_wls_result scalar wls_fit
    struct cdid_logit_result scalar lfit
    real scalar N, n_A, J, M_dates, R, n_pairs, n_links, r, p, l, k, b, s, tm
    real scalar p_or, p_trt, d_sel_max, ws_est, work_est, chunk_r, max_w_sel, id
    real scalar int_col_r, has_x_r, has_rho_r, has_just_r
    real scalar s_c, t_c, pi, g_pos, sid, pc_id, L_b, mass_b, n_pseudo, at, nn, scale_b
    struct cdid_sel_moments scalar acc_r
    real scalar p_kept, p_b, b_c1, b_c2, dres, sum_wp, sum_wm, dp, dm, T_ipw, M_ipw, mu_ipw, T_dr, M_dr, mu_dr, q_sum
    real scalar m_col, k_idx, c_pos, can_cache, tot_sup, at_sup, n_tm_cur, prior_tm, rk_h1, rk_h2, tot_obs_b, at_b, m_fast
    string scalar univ_r
    string colvector link_ids_str
    real colvector unit_to_a, a_to_unit, k_ctrl, obs_u, dy_u, ix_k, u_k, dy_k, q_k, beta_b
    real colvector B_r, U_r, rho_r, pair_s_vec, pair_t_vec, unknown_units, sup_map, ord_blk, first_link_tm
    real colvector pseudo_y, pseudo_w, pseudo_block, raw_res, fw, wp, wm, wm_act, prob, odds, odds_act
    real colvector r_obs, m_hat, y_trt, coh_loc, s_all, pc_all, cohs_m, sub_k, rate_vec, rate_c, rate_act, moment_cols
    real colvector obs_ru, obs_dt, obs_y, ix_m, u_all_b, ha_all_b, ord_b
    real rowvector constant_cols, rawZ_first, coef_rows, c_l, kept_b
    real colvector db_plus, db_minus, db_diff, wg_ipw, wg_dr
    real matrix pair_specs_r, sup_keys, blk_keys, blk_ctrl_mat, trt_keys, trt_coh_mat, pair_resp_mat
    real matrix F_s, F_b, Q_mat, R_mat, V_mat, res_fac, A_ps, V_ps, pseudo_X, D_s, D_kept, centered
    real matrix X_trt, Z_or_m, D_obs, score_trt, R_ipw_m, R_dr_m, V_ipw_m, V_dr_m, H_mat, H_chol, H_inv_m, z_sol_h, Y_dense
    real rowvector tau_q
    real colvector known_units,known_values
    real scalar scoped
    real colvector scope_r
    struct cdid_selected_universe_buffer colvector ubuf

    ctl = _cdid_sel_link_norm_controls(raw_ctl)
    /* used_scope=1: recipe-indexed slots carry values only on each recipe's actual
       supports; static callers omit it and keep the full-F validation below. */
    scoped = (args() >= 12 ? used_scope : 0)
    if (scoped != 0 & scoped != 1) _cdid_sel_link_fail("input", "used_scope must be 0 or 1.", 498)
    if (scoped & ctl.selection != "conditional_mtar")
        _cdid_sel_link_fail("unsupported", "Recipe-scoped registry slots are consumed only by conditional_mtar selection.", 459)
    if(population_ref==NULL) _cdid_sel_link_fail("input","A borrowed compact population object is required.",498)
    population=*population_ref;m.population_ref=population_ref
    N = rows(cohort_idx)
    J = rows(cohort_values)
    M_dates = rows(calendar_values)
    R = rows(recipes_in)
    n_pairs = rows(pairs_in)
    n_links = rows(link_specs)

    if (N < 1 | J < 1 | M_dates < 2 | R < 1 | n_pairs < 1 | n_links < 1 |
        rows(is_A) != N | rows(w_pop) != N) {
        _cdid_sel_link_fail("input", "Invalid frame, calendar, design block, pair, or link dimensions.", 459)
    }
    if(population.N!=N | rows(population.calendar)!=M_dates | any(population.calendar:!=calendar_values) |
        any(population.A:!=is_A))
        _cdid_sel_link_fail("input","Borrowed population, modeled calendar and original-F/A map disagree.",498)
    if (missing(cohort_idx) | any(cohort_idx :< 0) | any(cohort_idx :> J) |
        any(cohort_idx :!= floor(cohort_idx)))
        _cdid_sel_link_fail("input", "cohort_idx must be known positions1..J or0 for explicitly unknown off-A metadata.",459)
    known_units=_cdid_sel_link_index(population.cohort_known:==1)
    known_values=uniqrows(sort(population.cohort[known_units],1))
    if(rows(known_values)!=J) _cdid_sel_link_fail("input","Known cohort registry differs from the borrowed population.",498)
    if(any(known_values:!=cohort_values))
        _cdid_sel_link_fail("input","Known cohort values/order differ from the borrowed population.",498)
    unknown_units=_cdid_sel_link_index(population.cohort_known:==0)
    if(rows(unknown_units)) {
        if(any(cohort_idx[unknown_units]:!=0) | any(is_A[unknown_units]:!=0))
            _cdid_sel_link_fail("input","Unknown cohort requires position0 strictly off A.",459)
    }
    if(rows(known_units)) {
        if(any(cohort_idx[known_units]:==0))
            _cdid_sel_link_fail("input","Known cohort has no registry position.",498)
        if(any(population.cohort[known_units]:!=cohort_values[cohort_idx[known_units]]))
            _cdid_sel_link_fail("input","Cohort positions disagree with borrowed known metadata.",498)
    }
    if (missing(is_A) | any((is_A :!= 0) :& (is_A :!= 1))) {
        _cdid_sel_link_fail("input", "is_A must be a binary 0/1 indicator on 1..N.", 459)
    }
    n_A = sum(is_A)
    if (n_A < 1 | n_A > N) {
        _cdid_sel_link_fail("support", "At least one repeat-observer unit in A is required.", 459)
    }
    if (missing(w_pop) | any(w_pop :<= 0)) {
        _cdid_sel_link_fail("input", "Population weights w_pop must be finite and strictly positive.", 459)
    }

    m.profiles=data_in
    m.reader_sel_slot=m.reader_w_slot=m.reader_rho_slot=0
    for(r=1;r<=R;r++) {
        _cdid_sld_slot_check(recipes_in[r].weight_id,rows(m.profiles.weight_data))
        if(ctl.selection=="conditional_mtar") {
            _cdid_sld_slot_check(recipes_in[r].sel_id,rows(m.profiles.sel_data))
            _cdid_sld_slot_check(recipes_in[r].rho_id,rows(m.profiles.rho_data))
        }
        if(ctl.est_method!="ipw") _cdid_sld_slot_check(recipes_in[r].or_id,rows(m.profiles.or_data))
        if(ctl.est_method!="or") _cdid_sld_slot_check(recipes_in[r].trt_id,rows(m.profiles.trt_data))
    }
    p_or  = (ctl.est_method == "or" | ctl.est_method == "dr" ? cols(m.profiles.or_data[recipes_in[1].or_id].value) : 0)
    p_trt = (ctl.est_method == "ipw" | ctl.est_method == "dr" ? cols(m.profiles.trt_data[recipes_in[1].trt_id].value) : 0)
    if ((ctl.est_method == "or" | ctl.est_method == "dr") & (!ctl.or_intercept & p_or == 0)) {
        _cdid_sel_link_fail("input", "Empty outcome formula refused: no intercept and no slope columns.", 459)
    }
    if ((ctl.est_method == "ipw" | ctl.est_method == "dr") & p_trt < 1) {
        _cdid_sel_link_fail("input", "Treatment propensity design requires at least one column.", 459)
    }

    d_sel_max = 0
    for (r = 1; r <= R; r = r + 1) {
        if (recipes_in[r].recipe_id != r | strtrim(recipes_in[r].recipe_key) == "") {
            _cdid_sel_link_fail("input", "Recipes must have canonical 1..R IDs and nonempty keys.", 459)
        }
        if (rows(m.profiles.weight_data[recipes_in[r].weight_id].value) != N | (!scoped & (missing(m.profiles.weight_data[recipes_in[r].weight_id].value) | any(m.profiles.weight_data[recipes_in[r].weight_id].value :<= 0)))) {
            _cdid_sel_link_fail("input", "Recipe sampling weights must be finite and strictly positive on 1..N.", 459)
        }
        if (ctl.est_method == "or" | ctl.est_method == "dr") {
            if (rows(m.profiles.or_data[recipes_in[r].or_id].value) != N | cols(m.profiles.or_data[recipes_in[r].or_id].value) != p_or) {
                _cdid_sel_link_fail("input", "Recipe Z_or must have N rows and consistent p_or columns.", 459)
            }
        }
        if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
            if (rows(m.profiles.trt_data[recipes_in[r].trt_id].value) != N | cols(m.profiles.trt_data[recipes_in[r].trt_id].value) != p_trt) {
                _cdid_sel_link_fail("input", "Recipe V_trt must have N rows and consistent p_trt columns.", 459)
            }
        }
        if (ctl.selection == "conditional_mtar") {
            if (rows(m.profiles.sel_data[recipes_in[r].sel_id].value) != N) {
                _cdid_sel_link_fail("input", "Recipe Z_sel must have N rows when selection is active.", 459)
            }
            d_sel_max = max((d_sel_max, cols(m.profiles.sel_data[recipes_in[r].sel_id].value)))
        }
        else {
            univ_r = (strtrim(recipes_in[r].universe) == "" ? "F" : strtrim(recipes_in[r].universe))
            has_x_r = (recipes_in[r].has_xformla != 0 & recipes_in[r].has_xformla < . ? 1 : 0)
            has_rho_r = (recipes_in[r].has_rho != 0 & recipes_in[r].has_rho < . ? 1 : 0)
            has_just_r = (strtrim(recipes_in[r].justification) != "")
            (void) cdid_sel_norm_decl("none", has_x_r, univ_r, has_rho_r, has_just_r,
                ctl.selection_model, ctl.selection_pooling, recipes_in[r].justification)
        }
    }

    /* Pre-allocation workspace and work budget check */
    /* Declared double-entry accounting for retained model designs/maps and
       bounded pair scratch. This is an allocation estimate, not measured RSS.
       There is no automatic workspace cap; the private cap is caller supplied. */
    ws_est=N*(p_or+16)+n_links*(p_or+p_trt+12)+n_pairs*4
    for(id=1;id<=rows(data_in.sel_data);id++) ws_est=ws_est+rows(data_in.sel_data[id].value)*cols(data_in.sel_data[id].value)
    for(id=1;id<=rows(data_in.or_data);id++) ws_est=ws_est+rows(data_in.or_data[id].value)*cols(data_in.or_data[id].value)
    for(id=1;id<=rows(data_in.trt_data);id++) ws_est=ws_est+rows(data_in.trt_data[id].value)*cols(data_in.trt_data[id].value)
    for(id=1;id<=rows(data_in.weight_data);id++) ws_est=ws_est+rows(data_in.weight_data[id].value)
    for(id=1;id<=rows(data_in.rho_data);id++) ws_est=ws_est+rows(data_in.rho_data[id].value)
    work_est = 128 * (ctl.maxit + 2) * R * N * (d_sel_max + p_or + p_trt + 4)^2 +
               128 * n_links * N * (p_or + p_trt + 4)
    _cdid_sel_link_budget_check(ws_est, work_est, ctl.max_work, ctl.max_workspace_doubles, "cdid_selected_prepare_registry")

    /* 1. Populate Original-F Map and Population Share Moments */
    m.prepared              = 0
    m.recipe_keep           = J(0, 0, .)
    if (args() >= 13) {
        if (rows(recipe_keep) == N & cols(recipe_keep) == R) m.recipe_keep = recipe_keep
        else if (rows(recipe_keep) | cols(recipe_keep))
            _cdid_sel_link_fail("input", "Design block eligibility masks must have N rows and one column per design block.", 498)
    }
    m.outcome_pooling       = ctl.outcome_pooling
    m.reader_ready=0;m.reader_active=0;m.reader_universe=0
    m.ctl                   = ctl
    m.total_p_sep           = 0
    m.total_des_sep         = 0
    m.N                     = N
    m.n_A                   = n_A
    m.lambda                = n_A / N
    m.n_cohorts             = J
    m.n_dates               = M_dates
    m.R                     = R
    m.n_pairs               = n_pairs
    m.n_links               = n_links
    m.is_A                  = is_A
    m.cohort_idx            = cohort_idx
    m.cohort_values         = cohort_values
    m.calendar_values       = calendar_values
    m.w_pop                 = w_pop
    m.n_selection_fit_calls = 0
    m.n_selection_blocks = 0
    m.n_treatment_fit_calls = 0
    m.n_or_block_fits       = 0
    m.n_imputation_calls    = 0
    m.workspace_doubles     = ws_est
    m.work_ops              = work_est

    a_to_unit = _cdid_sel_link_index(is_A :== 1)
    unit_to_a = J(N, 1, 0)
    unit_to_a[a_to_unit] = (1::n_A)
    m.a_to_unit = a_to_unit
    m.unit_to_a = unit_to_a

    m.mean_w  = sum(w_pop) / N
    m.mean_wA = sum(w_pop[a_to_unit]) / N
    m.tau     = J(J, 1, 0)
    m.eta     = J(J, 1, 0)
    for (k = 1; k <= J; k = k + 1) {
        ix_k = _cdid_sel_link_index(cohort_idx :== k)
        if (rows(ix_k) > 0) m.tau[k] = (sum(w_pop[ix_k]) / N) / m.mean_w
        ix_k = _cdid_sel_link_index((cohort_idx :== k) :& (is_A :== 1))
        if (rows(ix_k) > 0) m.eta[k] = (sum(w_pop[ix_k]) / N) / m.mean_wA
    }

    /* 2. Validate Pairs and Consumer Link Specifications */
    for (p = 1; p <= n_pairs; p = p + 1) {
        s_c = pairs_in[p].s_cal
        t_c = pairs_in[p].t_cal
        if (s_c < 1 | t_c <= s_c | t_c > M_dates |
            s_c != floor(s_c) | t_c != floor(t_c)) {
            _cdid_sel_link_fail("support", "Pair calendar indices must lie in 1..M with s_cal < t_cal.", 459)
        }
        if(pairs_in[p].s_val!=calendar_values[s_c] | pairs_in[p].t_val!=calendar_values[t_c])
            _cdid_sel_link_fail("input","Supplied pair values disagree with original calendar positions.",498)
    }
    m.pairs = pairs_in

    /* Validate link_specs: distinct nonempty link_ids, valid pairs, cohorts, recipes */
    link_ids_str = J(n_links, 1, "")
    for (l = 1; l <= n_links; l = l + 1) {
        if (strtrim(link_specs[l].link_id) == "") {
            _cdid_sel_link_fail("input", "Distinct nonempty link IDs required.", 459)
        }
        link_ids_str[l] = link_specs[l].link_id
        r = link_specs[l].recipe_id
        if (r < 1 | r > R | r != floor(r)) {
            _cdid_sel_link_fail("input", "Link recipe_id outside 1..R.", 459)
        }
        g_pos = link_specs[l].target_cohort_pos
        k_ctrl = link_specs[l].control_cohort_pos
        if (g_pos < 1 | g_pos > J | g_pos != floor(g_pos) | rows(k_ctrl) < 1 |
            missing(k_ctrl) | any(k_ctrl :< 1) | any(k_ctrl :> J) |
            any(k_ctrl :!= floor(k_ctrl)) | any(k_ctrl :== g_pos)) {
            _cdid_sel_link_fail("input", "Invalid target or control cohort positions in link_specs.", 459)
        }
    }
    if (n_links > 1) {
        link_ids_str = sort(link_ids_str, 1)
        if (any(link_ids_str[|2 \ n_links|] :== link_ids_str[|1 \ (n_links - 1)|])) {
            _cdid_sel_link_fail("input", "Duplicate link_id in link_specs.", 459)
        }
    }

    /* Precompute full-F observed pair outcome differences in one pass when within cap */
    m.pair_full_ready = 0
    m.support_full_ready = 0
    can_cache = (ctl.max_workspace_doubles >= . |
        m.workspace_doubles + N * M_dates + (3 * n_pairs + 4 * n_links) * n_A <= ctl.max_workspace_doubles)
    if (can_cache) {
        Y_dense = J(N, M_dates, .)
        if (rows(population.observed_order) > 0) {
            obs_ru = population.ids.row_unit[population.observed_order]
            obs_dt = population.date[population.observed_order]
            obs_y  = population.y[population.observed_order]
            for (k = 1; k <= M_dates; k = k + 1) {
                ix_m = _cdid_sel_link_index(obs_dt :== calendar_values[k])
                if (rows(ix_m) > 0) Y_dense[obs_ru[ix_m], k] = obs_y[ix_m]
            }
        }
        m.pair_full = cdid_selected_rows(n_pairs, 1)
        for (p = 1; p <= n_pairs; p = p + 1) {
            s_c = pairs_in[p].s_cal
            t_c = pairs_in[p].t_cal
            ix_k = _cdid_sel_link_index((Y_dense[, s_c] :< .) :& (Y_dense[, t_c] :< .))
            m.pair_full[p].loc = ix_k
            m.pair_full[p].unit = ix_k
            m.pair_full[p].dy = (rows(ix_k) > 0 ? Y_dense[ix_k, t_c] - Y_dense[ix_k, s_c] : J(0, 1, 0))
            m.pair_full[p].q = J(0, 1, 0)
            m.pair_full[p].cohort_local = J(0, 1, 0)
            m.pair_full[p].score = m.pair_full[p].R_ipw = m.pair_full[p].R_dr = J(0, 0, 0)
        }
        Y_dense = J(0, 0, 0)
        m.pair_full_ready = 1
    }

    /* 3. Build Canonical Support, OR Block, and Treatment Model Catalogues */
    /* Map each link to its pair_id */
    m.links = cdid_selected_link_rec(n_links, 1)
    sup_keys = J(n_links * (1 + J), 3, 0) /* (pair_id, cohort_pos, recipe_id) */
    at_sup = 0

    for (l = 1; l <= n_links; l = l + 1) {
        pi = link_specs[l].pair_id
        if (pi == . | pi == 0) {
            pi = 0
            for (p = 1; p <= n_pairs; p = p + 1) {
                if (pairs_in[p].s_cal == link_specs[l].s_cal &
                    pairs_in[p].t_cal == link_specs[l].t_cal) {
                    if (pi != 0) _cdid_sel_link_fail("input",
                        "Equal-endpoint supplied pairs require an explicit opaque pair_id.", 498)
                    pi = p
                }
            }
        }
        if (pi < 1 | pi >= . | pi > n_pairs | pi != floor(pi))
            _cdid_sel_link_fail("input", "Link pair_id is not a supplied pair position.", 498)
        if (pairs_in[pi].s_cal != link_specs[l].s_cal |
            pairs_in[pi].t_cal != link_specs[l].t_cal)
            _cdid_sel_link_fail("input", "Explicit pair_id and link calendar indices disagree.", 498)
        r = link_specs[l].recipe_id
        g_pos = link_specs[l].target_cohort_pos
        k_ctrl = uniqrows(sort(link_specs[l].control_cohort_pos, 1))

        m.links[l].link_index         = l
        m.links[l].link_id            = link_specs[l].link_id
        m.links[l].functional_key     = link_specs[l].functional_key
        m.links[l].pair_id            = pi
        m.links[l].recipe_id          = r
        m.links[l].target_cohort_pos  = g_pos
        m.links[l].control_cohort_pos = k_ctrl
        m.links[l].L                  = pairs_in[pi].t_cal - pairs_in[pi].s_cal
        m.links[l].target_mu_ready    = 0

        at_sup = at_sup + 1
        sup_keys[at_sup, ] = (pi, g_pos, r)
        for (k = 1; k <= rows(k_ctrl); k = k + 1) {
            at_sup = at_sup + 1
            sup_keys[at_sup, ] = (pi, k_ctrl[k], r)
        }
    }

    /* Canonicalize supports in (pair_id, cohort_pos, recipe_id) order */
    sup_keys = uniqrows(sort(sup_keys[1..at_sup, ], (1, 2, 3)))
    m.n_supports = rows(sup_keys)
    m.supports = cdid_selected_support(m.n_supports, 1)
    tot_sup = n_pairs * J * R
    sup_map = (tot_sup <= 500000 ? J(tot_sup, 1, 0) : J(0, 1, 0))
    for (s = 1; s <= m.n_supports; s = s + 1) {
        m.supports[s].support_id = s
        m.supports[s].pair_id    = sup_keys[s, 1]
        m.supports[s].cohort_pos = sup_keys[s, 2]
        m.supports[s].recipe_id  = sup_keys[s, 3]
        if (rows(sup_map) > 0) {
            sup_map[(sup_keys[s, 1] - 1) * J * R + (sup_keys[s, 2] - 1) * R + sup_keys[s, 3]] = s
        }
    }

    /* Resolve target_support_id and control_support_ids for each link */
    for (l = 1; l <= n_links; l = l + 1) {
        pi = m.links[l].pair_id
        r  = m.links[l].recipe_id
        g_pos = m.links[l].target_cohort_pos
        k_ctrl = m.links[l].control_cohort_pos
        if (rows(sup_map) > 0) {
            m.links[l].target_support_id = sup_map[(pi - 1) * J * R + (g_pos - 1) * R + r]
            m.links[l].control_support_ids = J(rows(k_ctrl), 1, 0)
            for (k = 1; k <= rows(k_ctrl); k = k + 1) {
                m.links[l].control_support_ids[k] = sup_map[(pi - 1) * J * R + (k_ctrl[k] - 1) * R + r]
            }
        }
        else {
            m.links[l].target_support_id = 0
            for (s = 1; s <= m.n_supports; s = s + 1) {
                if (m.supports[s].pair_id == pi &
                    m.supports[s].cohort_pos == g_pos &
                    m.supports[s].recipe_id == r) {
                    m.links[l].target_support_id = s
                    break
                }
            }
            m.links[l].control_support_ids = J(rows(k_ctrl), 1, 0)
            for (k = 1; k <= rows(k_ctrl); k = k + 1) {
                for (s = 1; s <= m.n_supports; s = s + 1) {
                    if (m.supports[s].pair_id == pi &
                        m.supports[s].cohort_pos == k_ctrl[k] &
                        m.supports[s].recipe_id == r) {
                        m.links[l].control_support_ids[k] = s
                        break
                    }
                }
            }
        }
    }

    /* Deduplicate OR Control Blocks by (pair_id, sorted_controls, recipe_id)
       in canonical pair-major, lexicographic control-cohort, recipe order */
    blk_keys = J(n_links, 2 + J, 0)
    for (l = 1; l <= n_links; l = l + 1) {
        blk_keys[l, 1] = m.links[l].pair_id
        blk_keys[l, 2 + J] = m.links[l].recipe_id
        k_ctrl = m.links[l].control_cohort_pos
        for (k = 1; k <= rows(k_ctrl); k = k + 1) {
            blk_keys[l, 1 + k] = k_ctrl[k]
        }
    }
    ord_blk = order(blk_keys, (1..(2 + J)))
    blk_ctrl_mat = uniqrows(blk_keys[ord_blk, ])
    m.n_blocks = rows(blk_ctrl_mat)
    m.or_blocks = cdid_selected_or_block(m.n_blocks, 1)
    for (b = 1; b <= m.n_blocks; b = b + 1) {
        pi = blk_ctrl_mat[b, 1]
        r  = blk_ctrl_mat[b, 2 + J]
        k_ctrl = select(blk_ctrl_mat[b, 2..(1 + J)]', blk_ctrl_mat[b, 2..(1 + J)]' :> 0)
        m.or_blocks[b].block_id           = b
        m.or_blocks[b].pair_id            = pi
        m.or_blocks[b].recipe_id          = r
        m.or_blocks[b].L                  = pairs_in[pi].t_cal - pairs_in[pi].s_cal
        m.or_blocks[b].control_cohort_pos = k_ctrl
        m.or_blocks[b].support_ids        = J(rows(k_ctrl), 1, 0)
        m.or_blocks[b].blk_cache_ready    = 0
        if (rows(sup_map) > 0) {
            for (k = 1; k <= rows(k_ctrl); k = k + 1) {
                m.or_blocks[b].support_ids[k] = sup_map[(pi - 1) * J * R + (k_ctrl[k] - 1) * R + r]
            }
        }
        else {
            for (k = 1; k <= rows(k_ctrl); k = k + 1) {
                for (s = 1; s <= m.n_supports; s = s + 1) {
                    if (m.supports[s].pair_id == pi &
                        m.supports[s].cohort_pos == k_ctrl[k] &
                        m.supports[s].recipe_id == r) {
                        m.or_blocks[b].support_ids[k] = s
                        break
                    }
                }
            }
        }
    }
    b = 1
    m.links[ord_blk[1]].block_id = 1
    for (l = 2; l <= n_links; l = l + 1) {
        if (any(blk_keys[ord_blk[l], ] :!= blk_keys[ord_blk[l - 1], ])) b = b + 1
        m.links[ord_blk[l]].block_id = b
    }

    /* Deduplicate Treatment Propensity Models by (pair_id, target_cohort, sorted_controls, recipe_id)
       preserving first-seen link traversal order (matching dr_build_treatment_catalog) */
    trt_keys = J(n_links, 3 + J, 0)
    first_link_tm = J(n_links, 1, 0)
    n_tm_cur = 0
    for (l = 1; l <= n_links; l = l + 1) {
        trt_coh_mat = J(1, 3 + J, 0)
        trt_coh_mat[1, 1] = m.links[l].pair_id
        trt_coh_mat[1, 2] = m.links[l].target_cohort_pos
        trt_coh_mat[1, 3 + J] = m.links[l].recipe_id
        k_ctrl = m.links[l].control_cohort_pos
        for (k = 1; k <= rows(k_ctrl); k = k + 1) {
            trt_coh_mat[1, 2 + k] = k_ctrl[k]
        }
        tm = 0
        for (k = 1; k <= n_tm_cur; k = k + 1) {
            if (all(trt_keys[k, ] :== trt_coh_mat)) {
                tm = k
                break
            }
        }
        if (tm == 0) {
            n_tm_cur = n_tm_cur + 1
            trt_keys[n_tm_cur, ] = trt_coh_mat
            first_link_tm[n_tm_cur] = l
            tm = n_tm_cur
        }
        m.links[l].trt_id = tm
    }
    m.n_trt_models = n_tm_cur
    m.trt_models = cdid_selected_trt_model(m.n_trt_models, 1)
    for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
        l = first_link_tm[tm]
        m.trt_models[tm].trt_id            = tm
        m.trt_models[tm].pair_id           = m.links[l].pair_id
        m.trt_models[tm].recipe_id         = m.links[l].recipe_id
        m.trt_models[tm].block_id          = m.links[l].block_id
        m.trt_models[tm].target_cohort_pos = m.links[l].target_cohort_pos
        m.trt_models[tm].cohorts           = m.links[l].target_cohort_pos \ m.links[l].control_cohort_pos
        m.trt_models[tm].support_ids       = m.links[l].target_support_id \ m.links[l].control_support_ids
        m.trt_models[tm].p_trt             = p_trt
        m.trt_models[tm].rep_tid           = tm
        m.trt_models[tm].m_col_fast        = 0
        m.trt_models[tm].risk_cache_ready  = 0
    }

    /* 4. Build Recipe Contexts and Fit Context-Shared Selection Models */
    m.contexts = cdid_selected_context(R, 1)
    for (r = 1; r <= R; r = r + 1) {
        m.contexts[r].recipe_id        = r
        m.contexts[r].recipe_key       = recipes_in[r].recipe_key
        m.contexts[r].design_cal       = recipes_in[r].design_cal
        m.contexts[r].weight_cal       = recipes_in[r].weight_cal
        m.contexts[r].sel_id=recipes_in[r].sel_id
        m.contexts[r].or_id=recipes_in[r].or_id
        m.contexts[r].trt_id=recipes_in[r].trt_id
        m.contexts[r].weight_id=recipes_in[r].weight_id
        m.contexts[r].rho_id=recipes_in[r].rho_id
        m.contexts[r].m_col_beta_start = 0
        m.contexts[r].m_col_beta_end   = 0
        m.contexts[r].m_col_obj        = 0
        m.contexts[r].n_moments_total  = 0
        m.contexts[r].n_loading_pieces = 0
        m.contexts[r].acc_full_ready   = 0

        if (ctl.selection == "conditional_mtar") {
            univ_r     = (strtrim(recipes_in[r].universe) == "" ? "F" : strtrim(recipes_in[r].universe))
            int_col_r  = (recipes_in[r].sel_has_intercept ?
                          (recipes_in[r].sel_int_col > 0 & recipes_in[r].sel_int_col < . ? recipes_in[r].sel_int_col : 1) : 0)
            has_x_r    = 1
            has_rho_r  = (recipes_in[r].has_rho != 0 & recipes_in[r].has_rho < . ? 1 :
                          (univ_r != "F" | (rows(_cdid_sld_rho_full(m,recipes_in[r].rho_id)) == N & any(_cdid_sld_rho_full(m,recipes_in[r].rho_id) :!= 1))))
            has_just_r = (strtrim(recipes_in[r].justification) != "")
            decl = cdid_sel_norm_decl("conditional_mtar", has_x_r, univ_r,
                has_rho_r, has_just_r, ctl.selection_model,
                ctl.selection_pooling, recipes_in[r].justification)

            B_r=population.B
            U_r=(univ_r=="F" ? J(N,1,1) : (univ_r=="B" ? B_r : is_A))
            rho_r=_cdid_sld_rho_full(m,recipes_in[r].rho_id)
            chunk_r = (ctl.max_chunk_rows > 0 ? ctl.max_chunk_rows : N)

            if(scoped) {
                /* The recipe's selection cells are exactly its links' target/control cohorts. */
                scope_r=J(J,1,0)
                for(l=1;l<=n_links;l++) if(m.links[l].recipe_id==r) {
                    scope_r[m.links[l].target_cohort_pos]=1
                    scope_r[m.links[l].control_cohort_pos]=J(rows(m.links[l].control_cohort_pos),1,1)
                }
                if(rows(m.recipe_keep)) temporary_reader=cdid_sel_reader_scope(population,_cdid_sld_sel_full(m,recipes_in[r].sel_id),
                    recipes_in[r].sel_has_intercept,int_col_r,_cdid_sld_w_full(m,recipes_in[r].weight_id),decl,rho_r,chunk_r,scope_r,
                    m.recipe_keep[,r])
                else temporary_reader=cdid_sel_reader_scope(population,_cdid_sld_sel_full(m,recipes_in[r].sel_id),
                    recipes_in[r].sel_has_intercept,int_col_r,_cdid_sld_w_full(m,recipes_in[r].weight_id),decl,rho_r,chunk_r,scope_r)
            }
            else temporary_reader=cdid_sel_reader_pop(population,_cdid_sld_sel_full(m,recipes_in[r].sel_id),
                recipes_in[r].sel_has_intercept,int_col_r,_cdid_sld_w_full(m,recipes_in[r].weight_id),decl,rho_r,chunk_r)
            _cdid_selected_reader_store(m,r,temporary_reader)
            /* A restricted recipe reads its own roster U E_g through a dedicated universe buffer. */
            if(rows(m.recipe_keep)) if(!all(m.recipe_keep[,r])) {
                U_r=U_r:*(m.recipe_keep[,r]:!=0)
                _cdid_selected_use_universe(m,r,U_r)
            }
            if(any(m.reader.U:!=U_r) | any(m.reader.B:!=B_r))
                _cdid_sel_link_fail("input","Declared F/B/A design block universe disagrees with compact population.",498)

            /* Collect active pair-cohorts for recipe r in link traversal order */
            at_sup = 0
            for (l = 1; l <= n_links; l = l + 1) {
                if (m.links[l].recipe_id == r) at_sup = at_sup + 1 + rows(m.links[l].control_cohort_pos)
            }
            pair_specs_r = J(at_sup, 5, 0)
            at_sup = 0
            for (l = 1; l <= n_links; l = l + 1) {
                if (m.links[l].recipe_id != r) continue
                pi = m.links[l].pair_id
                cohs_m = m.links[l].target_cohort_pos \ m.links[l].control_cohort_pos
                for (k = 1; k <= rows(cohs_m); k = k + 1) {
                    at_sup = at_sup + 1
                    pair_specs_r[at_sup, ] =
                        (cohs_m[k], pairs_in[pi].s_cal, pairs_in[pi].t_cal,
                         pairs_in[pi].s_val, pairs_in[pi].t_val)
                }
            }
            m.contexts[r].cell_map = cdid_sel_partition(m.reader, pair_specs_r)
            m.contexts[r].pc_lookup = J(M_dates * M_dates, J, 0)
            for (k = 1; k <= rows(m.contexts[r].cell_map.pair_cohorts); k = k + 1) {
                m.contexts[r].pc_lookup[(m.contexts[r].cell_map.pair_cohorts[k, 2] - 1) * M_dates +
                    m.contexts[r].cell_map.pair_cohorts[k, 3],
                    m.contexts[r].cell_map.pair_cohorts[k, 1]] = k
            }
            max_w_sel = (ctl.max_work > 0 & ctl.max_work < . ? ctl.max_work : .)
            (void) cdid_sel_budget_check(m.reader, m.contexts[r].cell_map, max_w_sel)
            m.contexts[r].sel_fit = cdid_sel_fit(m.reader, m.contexts[r].cell_map,
                ctl.eps, ctl.maxit)
            m.contexts[r].up_map = cdid_sel_up_map_build(m.contexts[r].sel_fit)
            m.n_selection_fit_calls=m.n_selection_fit_calls+1
            m.n_selection_blocks=m.n_selection_blocks+m.contexts[r].sel_fit.n_blocks
        }
    }

    /* 5. Evaluate Support Pair Observations, Selection Weights q, and Support Moments */
    if (can_cache) {
        m.support_full = cdid_selected_rows(m.n_supports, 1)
        if (R == 1) {
            for (p = 1; p <= n_pairs; p = p + 1) {
                m.pair_full[p].q = J(rows(m.pair_full[p].unit), 1, 0)
            }
        }
    }
    for (s = 1; s <= m.n_supports; s = s + 1) {
        pi = m.supports[s].pair_id
        c_pos = m.supports[s].cohort_pos
        r = m.supports[s].recipe_id
        m.supports[s].ld_ready = 0

        if (s == 1) {
            read_rows = (m.pair_full_ready ? m.pair_full[pi] : _cdid_selected_pair_rows(m, pi, (1::N)))
            obs_u = read_rows.unit; dy_u = read_rows.dy
        }
        else if (pi != m.supports[s - 1].pair_id) {
            read_rows = (m.pair_full_ready ? m.pair_full[pi] : _cdid_selected_pair_rows(m, pi, (1::N)))
            obs_u = read_rows.unit; dy_u = read_rows.dy
        }
        ix_k  = _cdid_sel_link_index((cohort_idx[obs_u] :== c_pos) :& _cdid_sld_keep(m, r, obs_u))
        if (rows(ix_k) < 1) {
            _cdid_sel_link_fail("support",
                "Positive actual-pair cohort support required on every declared support.", 459)
        }
        u_k  = obs_u[ix_k]
        dy_k = dy_u[ix_k]
        m.supports[s].n_obs     = rows(u_k)

        /* Check finite outcome/treatment designs on observed pair rows */
        if ((ctl.est_method == "or" | ctl.est_method == "dr") & p_or > 0) {
            Z_or_m = _cdid_sld_or_rows(m,r,u_k)
            if (missing(Z_or_m)) {
                _cdid_sel_link_fail("support",
                    "Missing or nonfinite outcome design values on observed pair support.", 459)
            }
        }
        if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
            if (missing(_cdid_sld_trt_rows(m,r,u_k))) {
                _cdid_sel_link_fail("support",
                    "Missing or nonfinite treatment design values on observed pair support.", 459)
            }
        }

        if (ctl.selection == "conditional_mtar") {
            /* Match pc_id in contexts[r].cell_map.pair_cohorts */
            pc_id = m.contexts[r].pc_lookup[(pairs_in[pi].s_cal - 1) * M_dates + pairs_in[pi].t_cal, c_pos]
            if (pc_id == 0) {
                _cdid_sel_link_fail("integrity", "Support pair-cohort missing from selection cell_map.", 498)
            }
            m.supports[s].pc_id = pc_id
            _cdid_selected_activate_reader(m,r)
            if (!ctl.points_only & can_cache) {
                m.supports[s].ld = cdid_sel_log_deriv_eval(m.contexts[r].sel_fit, m.reader, m.contexts[r].up_map, pc_id, u_k)
                m.supports[s].ld_g_idx = _cdid_sel_pc_global_idx(m.contexts[r].sel_fit, m.contexts[r].up_map, pc_id)
                m.supports[s].ld_ready = 1
                q_k = (m.reader.w[u_k] :* m.supports[s].ld.a_k) :/ (m.reader.rho[u_k] :* m.supports[s].ld.p_weight)
                if (missing(q_k) | any(q_k :<= 0)) {
                    _cdid_sel_link_fail("numerical", "Nonfinite or nonpositive fitted selection weights q.", 504)
                }
                m.supports[s].a_k           = m.supports[s].ld.a_k
                m.supports[s].mass_U_k      = m.supports[s].ld.mass_U_k
                m.supports[s].repeat_mass_k = m.supports[s].ld.a_k * m.supports[s].ld.mass_U_k * N
            }
            else {
                pr = cdid_sel_pred_pair(m.contexts[r].sel_fit, m.reader, pc_id, u_k)
                q_k = pr.sel_weight
                if (missing(q_k) | any(q_k :<= 0)) {
                    _cdid_sel_link_fail("numerical", "Nonfinite or nonpositive fitted selection weights q.", 504)
                }
                m.supports[s].a_k           = pr.a_k
                m.supports[s].mass_U_k      = pr.mass_U_k
                m.supports[s].repeat_mass_k = pr.a_k * pr.mass_U_k * N
            }
            m.supports[s].rate_baseline = 1
            m.supports[s].rate_c        = 1
            m.supports[s].rate_active   = 0
        }
        else {
            /* Baseline selection == "none" delegation */
            ix_k = _cdid_sel_link_index((cohort_idx :== c_pos) :& (is_A :== 1))
            m.supports[s].repeat_mass_k = sum(_cdid_sld_w_rows(m,r,ix_k))
            m.supports[s].mass_U_k      = m.supports[s].repeat_mass_k / N
            m.supports[s].a_k           = 1
            m.supports[s].pc_id         = 0
            m.supports[s].rate_baseline = sum(_cdid_sld_w_rows(m,r,u_k)) / m.supports[s].repeat_mass_k
            V_mat = cdid_ipw_rate_values(J(1, 1, m.supports[s].rate_baseline))
            m.supports[s].rate_c      = V_mat[1, 1]
            m.supports[s].rate_active = V_mat[1, 2]
            q_k = _cdid_sld_w_rows(m,r,u_k) :/ m.supports[s].rate_c
        }

        m.supports[s].weight_mass = sum(q_k)
        m.supports[s].sum_y       = sum(q_k :* dy_k)
        if (m.supports[s].weight_mass <= 0 | missing(m.supports[s].weight_mass)) {
            _cdid_sel_link_fail("support", "Nonpositive or nonfinite support weight mass.", 459)
        }
        if (can_cache) {
            m.support_full[s].loc  = u_k
            m.support_full[s].unit = u_k
            m.support_full[s].dy   = dy_k
            m.support_full[s].q    = q_k
            m.support_full[s].cohort_local = J(0, 1, 0)
            m.support_full[s].score = m.support_full[s].R_ipw = m.support_full[s].R_dr = J(0, 0, 0)
            if (R == 1) m.pair_full[pi].q[ix_k] = q_k
        }

        if (ctl.est_method == "or" | ctl.est_method == "dr") {
            L_b = pairs_in[pi].t_cal - pairs_in[pi].s_cal
            if (p_or > 0) {
                if (can_cache) m.supports[s].Z_or = Z_or_m
                D_s = L_b :* Z_or_m
                m.supports[s].sum_z   = quadcross(q_k, Z_or_m)
                m.supports[s].xmin    = colmin(Z_or_m)
                m.supports[s].xmax    = colmax(Z_or_m)
                m.supports[s].first_z = Z_or_m[1, ]
                m.supports[s].colmax  = colmax(abs(Z_or_m))
                m.supports[s].F_qr    = _cdid_or_compress((J(rows(u_k), 1, 1), D_s, dy_k) :* sqrt(q_k))
            }
            else {
                m.supports[s].sum_z   = J(1, 0, 0)
                m.supports[s].xmin    = J(1, 0, 0)
                m.supports[s].xmax    = J(1, 0, 0)
                m.supports[s].first_z = J(1, 0, 0)
                m.supports[s].colmax  = J(1, 0, 0)
                m.supports[s].F_qr    = _cdid_or_compress((J(rows(u_k), 1, 1), dy_k) :* sqrt(q_k))
            }
        }
    }
    if (can_cache) m.support_full_ready = 1

    /* Populate pc_ids on treatment models */
    for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
        s_all = m.trt_models[tm].support_ids
        pc_all = J(rows(s_all), 1, 0)
        for (k = 1; k <= rows(s_all); k = k + 1) {
            pc_all[k] = m.supports[s_all[k]].pc_id
        }
        m.trt_models[tm].pc_ids = pc_all
    }

    /* 6. Joint Untreated-Change OR Fit Across All Scientific Control Blocks & Recipes */
    m.or_shared.intercept         = ctl.or_intercept
    m.or_shared.p_raw             = p_or
    m.or_shared.p_kept            = 0
    m.or_shared.rss               = .
    m.or_shared.objective         = .
    m.or_shared.verification_scan = 0
    m.or_shared.kept              = J(1, 0, 0)
    m.or_shared.dropped           = J(1, 0, 0)
    m.or_shared.alias_map         = J(0, 0, 0)
    m.or_shared.alias_constants   = J(m.n_blocks, 0, 0)
    m.or_shared.beta              = J(0, 1, 0)
    m.or_shared.R                 = J(0, 0, 0)
    m.or_shared.pivot             = J(1, 0, 0)
    m.or_shared.beta_score_F      = J(0, 0, 0)
    m.or_shared.shared_gram_inv   = J(0, 0, 0)
    m.or_shared.shared_gram_ready = 0

    if (ctl.est_method == "or" | ctl.est_method == "dr") {
        if (ctl.outcome_pooling == "separate") {
            _cdid_sel_fit_or_separate(m)
        }
        else {
        n_pseudo = 2 * (p_or + 1) * m.n_blocks
        pseudo_X     = J(n_pseudo, p_or, 0)
        pseudo_y     = J(n_pseudo, 1, 0)
        pseudo_w     = J(n_pseudo, 1, 0)
        pseudo_block = J(n_pseudo, 1, 0)
        at = 0

        for (b = 1; b <= m.n_blocks; b = b + 1) {
            L_b = m.or_blocks[b].L
            F_b = J(0, p_or + 2, 0)
            mass_b = 0
            constant_cols = J(1, p_or, 1)
            rawZ_first = J(1, p_or, 0)
            s_all = m.or_blocks[b].support_ids
            for (k = 1; k <= rows(s_all); k = k + 1) {
                sid = s_all[k]
                F_b = _cdid_or_compress(F_b \ m.supports[sid].F_qr)
                mass_b = mass_b + m.supports[sid].weight_mass
                if (k == 1 & p_or > 0) rawZ_first = m.supports[sid].first_z
                if (p_or > 0) {
                    constant_cols = constant_cols :*
                        (m.supports[sid].xmin :== rawZ_first) :*
                        (m.supports[sid].xmax :== rawZ_first)
                }
            }
            hqrd(F_b[, 1], Q_mat = ., tau_q = ., R_mat = .)
            V_mat = hqrdmultq(Q_mat, tau_q, F_b[, (2..(p_or + 2))], 1)
            coef_rows = V_mat[1, ] / R_mat[1, 1]
            res_fac   = (rows(V_mat)>1 ? V_mat[2..rows(V_mat), ] : J(1,p_or+1,0))
            if (ctl.or_intercept & p_or > 0) {
                for (k = 1; k <= p_or; k = k + 1) {
                    if (constant_cols[k]) {
                        coef_rows[k] = L_b * rawZ_first[k]
                        res_fac[, k] = J(rows(res_fac), 1, 0)
                    }
                }
            }
            nn = rows(res_fac)
            scale_b = sqrt(nn / mass_b)
            A_ps = coef_rows :+ scale_b :* res_fac
            V_ps = coef_rows :- scale_b :* res_fac
            ix_k = ((at + 1)::(at + 2 * nn))
            if (p_or > 0) pseudo_X[ix_k, ] = A_ps[, 1..p_or] \ V_ps[, 1..p_or]
            pseudo_y[ix_k]     = A_ps[, p_or + 1] \ V_ps[, p_or + 1]
            pseudo_w[ix_k]     = J(2 * nn, 1, mass_b / (2 * nn))
            pseudo_block[ix_k] = J(2 * nn, 1, b)
            at = at + 2 * nn
        }

        /* Thin compressed factors can have fewer rows than the allocation
           upper bound. Only filled pseudo rows represent the WLS objective. */
        n_pseudo=at
        pseudo_X=pseudo_X[1..at,];pseudo_y=pseudo_y[1..at]
        pseudo_w=pseudo_w[1..at];pseudo_block=pseudo_block[1..at]
        wls_fit = cdid_pooled_wls_fit(pseudo_X, pseudo_y, pseudo_w, pseudo_block,
            J(n_pseudo, 1, 1), 1, ctl.or_intercept, J(0, p_or, 0), J(0, 1, 0), ctl.rank_tol)
        m.n_or_block_fits = 1

        m.or_shared.beta            = wls_fit.beta
        m.or_shared.R               = wls_fit.R
        m.or_shared.pivot           = wls_fit.pivot
        m.or_shared.kept            = wls_fit.kept
        m.or_shared.dropped         = wls_fit.dropped
        m.or_shared.alias_map       = wls_fit.alias_map
        m.or_shared.alias_constants = wls_fit.alias_constants
        p_kept = rows(wls_fit.beta)
        m.or_shared.p_kept          = p_kept
        if (!ctl.points_only) {
            m.or_shared.shared_gram_inv = (p_kept > 0 ?
                (N :* _cdid_nuisance_gram(m.or_shared.R, m.or_shared.pivot, I(p_kept)))' : J(0, 0, 0))
            m.or_shared.shared_gram_ready = 1
        }

        for (b = 1; b <= m.n_blocks; b = b + 1) {
            m.or_blocks[b].alpha      = (ctl.or_intercept ? wls_fit.alpha[b] : 0)
            m.or_blocks[b].weight_sum = wls_fit.mass[b]
            m.or_blocks[b].mass       = wls_fit.mass[b] / N
            m.or_blocks[b].mean_y     = (ctl.or_intercept ? wls_fit.mean_y[b] : 0)
            m.or_blocks[b].mean_d     = (ctl.or_intercept & p_kept > 0 ?
                wls_fit.mean_X[b, wls_fit.kept] : J(1, p_kept, 0))
            m.or_blocks[b].m_col_alpha = 0
            m.or_blocks[b].blk_cache_ready = 0
        }

        /* Verify aliases on actual observed control & required treated prediction rows */
        _cdid_sel_link_verify_aliases(m)

        /* Compute exact objective RSS/N and direct unit scores across all blocks */
        if (!ctl.points_only) {
            m.or_shared.beta_score_F  = J(N, p_kept, 0)
            if (ctl.selection == "none") {
                V_mat = J(p_kept, m.n_supports, 0)
                m.or_shared.baseline_objective_rate_loading=J(m.n_supports,1,0)
            }
        }
        m.or_shared.rss = 0
        for (b = 1; b <= m.n_blocks; b = b + 1) {
            L_b = m.or_blocks[b].L
            r   = m.or_blocks[b].recipe_id
            s_all = m.or_blocks[b].support_ids
            if (!ctl.points_only & ctl.selection == "none") {
                m.or_blocks[b].baseline_alpha_rate = J(rows(s_all), 1, 0)
            }
            if (!ctl.points_only & can_cache & ctl.or_intercept) {
                tot_obs_b = 0
                for (k = 1; k <= rows(s_all); k = k + 1) tot_obs_b = tot_obs_b + m.supports[s_all[k]].n_obs
                u_all_b  = J(tot_obs_b, 1, 0)
                ha_all_b = J(tot_obs_b, 1, 0)
                at_b = 0
            }
            if (!m.support_full_ready) {
                if (b == 1) pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
                else if (m.or_blocks[b].pair_id != m.or_blocks[b - 1].pair_id) {
                    pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
                }
            }
            for (k = 1; k <= rows(s_all); k = k + 1) {
                sid = s_all[k]
                if (m.support_full_ready) {
                    u_k = m.support_full[sid].unit; dy_k = m.support_full[sid].dy; q_k = m.support_full[sid].q
                }
                else {
                    read_rows = _cdid_sel_support_from_pair(m, sid, pair_rows)
                    u_k = read_rows.unit; dy_k = read_rows.dy; q_k = read_rows.q
                }
                if (p_kept > 0) {
                    D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                        L_b :* m.supports[sid].Z_or[, m.or_shared.kept] :
                        L_b :* _cdid_sld_or_rows(m,r,u_k,m.or_shared.kept))
                }
                else D_kept = J(rows(u_k), 0, 0)
                raw_res = dy_k :- m.or_blocks[b].alpha :-
                    (p_kept > 0 ? D_kept * m.or_shared.beta : J(rows(u_k), 1, 0))
                m.or_shared.rss = m.or_shared.rss + sum(q_k :* (raw_res :^ 2))

                if (!ctl.points_only) {
                    centered = (ctl.or_intercept & p_kept > 0 ? D_kept :- m.or_blocks[b].mean_d : D_kept)
                    if (p_kept > 0) {
                        m.or_shared.beta_score_F[u_k, ] = m.or_shared.beta_score_F[u_k, ] +
                            centered :* (q_k :* raw_res)
                    }
                    if (can_cache & ctl.or_intercept & rows(u_k) > 0) {
                        u_all_b[|(at_b + 1) \ (at_b + rows(u_k))|]  = u_k
                        ha_all_b[|(at_b + 1) \ (at_b + rows(u_k))|] = q_k :* raw_res
                        at_b = at_b + rows(u_k)
                    }
                    if (ctl.selection == "none" & m.supports[sid].rate_active) {
                        m.or_shared.baseline_objective_rate_loading[sid]=
                            m.or_shared.baseline_objective_rate_loading[sid]-
                            sum(q_k:*(raw_res:^2))/(N*m.supports[sid].rate_c)
                        if (ctl.or_intercept) {
                            m.or_blocks[b].baseline_alpha_rate[k] =
                                -sum(q_k :* raw_res) / m.supports[sid].rate_c
                        }
                        if (p_kept > 0) {
                            V_mat[, sid] = V_mat[, sid] -
                                (quadcross(q_k :* raw_res, centered))' / m.supports[sid].rate_c
                        }
                    }
                }
            }
            if (!ctl.points_only & can_cache & ctl.or_intercept) {
                if (tot_obs_b > 1) {
                    ord_b = order(u_all_b, 1)
                    m.or_blocks[b].blk_obs_u = u_all_b[ord_b]
                    m.or_blocks[b].blk_ha    = ha_all_b[ord_b]
                }
                else {
                    m.or_blocks[b].blk_obs_u = u_all_b
                    m.or_blocks[b].blk_ha    = ha_all_b
                }
                m.or_blocks[b].blk_hb = J(0, 0, 0)
                m.or_blocks[b].blk_gram_inv = J(0, 0, 0)
                m.or_blocks[b].blk_cache_ready = 1
            }
        }
        m.or_shared.objective = m.or_shared.rss / N
        if (!ctl.points_only & ctl.selection == "none") {
            m.or_shared.baseline_beta_rate_loading =
                _cdid_nuisance_gram(m.or_shared.R, m.or_shared.pivot, V_mat)
        }

        /* Compute OR target link moments (mu_or, d_target, num_target, c_target) */
        for (l = 1; l <= n_links; l = l + 1) {
            b   = m.links[l].block_id
            sid = m.links[l].target_support_id
            L_b = m.links[l].L
            mass_b = m.supports[sid].weight_mass
            c_l = (p_kept > 0 ? L_b :* m.supports[sid].sum_z[m.or_shared.kept] :/ mass_b : J(1, 0, 0))
            dres = m.supports[sid].sum_y -
                (ctl.or_intercept ? m.or_blocks[b].alpha * mass_b : 0) -
                (p_kept > 0 ? L_b * (m.supports[sid].sum_z[m.or_shared.kept] * m.or_shared.beta) : 0)
            m.links[l].mu_or      = dres / mass_b
            m.links[l].c_target   = c_l
            m.links[l].d_target   = mass_b / N
            m.links[l].num_target = m.links[l].mu_or * m.links[l].d_target
            m.links[l].baseline_or_rate_loading = 0
            m.links[l].baseline_den_rate_loading=0
            m.links[l].baseline_design_rate_loading=J(1,p_kept,0)
            m.links[l].target_mu_ready = 0
            if (ctl.selection == "none" & m.supports[sid].rate_active) {
                m.links[l].baseline_den_rate_loading=-m.links[l].d_target/m.supports[sid].rate_c
                if(p_kept) m.links[l].baseline_design_rate_loading=
                    -(L_b*m.supports[sid].sum_z[m.or_shared.kept]-mass_b*c_l)/
                    (N*m.supports[sid].rate_c)
                m.links[l].baseline_or_rate_loading =
                    -(dres - m.links[l].mu_or * mass_b) / (N * m.supports[sid].rate_c * m.links[l].d_target)
            }
            if (!ctl.points_only & can_cache & ctl.est_method == "or" & m.support_full_ready) {
                r = m.links[l].recipe_id
                u_k = m.support_full[sid].unit
                dy_k = m.support_full[sid].dy
                q_k = m.support_full[sid].q
                if (p_kept > 0) {
                    D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                        L_b :* m.supports[sid].Z_or[, m.or_shared.kept] :
                        L_b :* _cdid_sld_or_rows(m, r, u_k, m.or_shared.kept))
                }
                else D_kept = J(rows(u_k), 0, 0)
                raw_res = dy_k :- (ctl.or_intercept ? m.or_blocks[b].alpha : 0) :-
                    (p_kept > 0 ? D_kept * m.or_shared.beta : J(rows(u_k), 1, 0))
                m.links[l].target_mu_vec = q_k :* (raw_res :- m.links[l].mu_or) / m.links[l].d_target
                m.links[l].target_mu_ready = 1
            }
        }
        }
    }

    /* 7. Fit Deduplicated Risk-Set Treatment Propensity Models (IPW & DR) */
    if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
        p_kept = m.or_shared.p_kept
        for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
            pi     = m.trt_models[tm].pair_id
            r      = m.trt_models[tm].recipe_id
            b      = m.trt_models[tm].block_id
            cohs_m = m.trt_models[tm].cohorts
            s_all  = m.trt_models[tm].support_ids

            read_rows=(m.pair_full_ready ? _cdid_selected_treatment_rows(m,tm,J(0,1,0),0,m.pair_full) :
                _cdid_selected_treatment_rows(m,tm,(1::N),0))
            u_k=read_rows.unit;dy_k=read_rows.dy;coh_loc=read_rows.cohort_local
            y_trt=(coh_loc:==1);q_k=read_rows.q
            if(rows(u_k)<2) _cdid_sel_link_fail("support","Risk set has fewer than two observed units.",459)

            X_trt = _cdid_sld_trt_rows(m,r,u_k)
            q_sum = sum(q_k)
            prior_tm = 0
            if (can_cache) {
                for (k = 1; k < tm; k = k + 1) {
                    if (m.trt_models[k].fit_q_sum != q_sum) continue
                    if (m.trt_models[k].rep_tid != k) continue
                    if (m.trt_models[k].target_cohort_pos != m.trt_models[tm].target_cohort_pos) continue
                    if (m.contexts[m.trt_models[k].recipe_id].trt_id != m.contexts[r].trt_id) continue
                    if (rows(m.trt_models[k].cohorts) != rows(cohs_m)) continue
                    if (rows(m.trt_models[k].fit_q_cache) != rows(q_k)) continue
                    if (!all(m.trt_models[k].cohorts :== cohs_m)) continue
                    if (!all(m.trt_models[k].fit_u_cache :== u_k)) continue
                    if (!all(m.trt_models[k].fit_q_cache :== q_k)) continue
                    prior_tm = k
                    break
                }
            }

            m.n_treatment_fit_calls = m.n_treatment_fit_calls + 1
            if (prior_tm > 0) {
                m.trt_models[tm].rep_tid = prior_tm
                beta_b = m.trt_models[prior_tm].beta_trt
                H_chol = m.trt_models[prior_tm].H_chol
                H_inv_m = m.trt_models[prior_tm].H_inv
                V_mat = _cdid_sel_link_logit_vals(X_trt * beta_b)
                prob = V_mat[, 1]; odds = V_mat[, 2]; odds_act = V_mat[, 3]
                wp     = q_k :* y_trt
                wm     = q_k :* (1 :- y_trt) :* odds
                sum_wp = m.trt_models[prior_tm].sum_wp
                sum_wm = m.trt_models[prior_tm].sum_wm
                dp     = m.trt_models[prior_tm].dp
                dm     = m.trt_models[prior_tm].dm
                score_trt = (q_k :* (y_trt :- prob)) :* X_trt
                wm_act    = wm :* odds_act
                if (m.trt_models[prior_tm].pair_id == pi) {
                    T_ipw   = m.trt_models[prior_tm].T_ipw
                    M_ipw   = m.trt_models[prior_tm].M_ipw
                    mu_ipw  = m.trt_models[prior_tm].mu_ipw
                    R_ipw_m = (wp :* (dy_k :- T_ipw) / dp, wm :* (dy_k :- M_ipw) / dm, wp, wm)
                    V_ipw_m = m.trt_models[prior_tm].V_ipw
                    wg_ipw  = m.trt_models[prior_tm].w_gamma_ipw
                }
                else {
                    T_ipw  = sum(wp :* dy_k) / sum_wp
                    M_ipw  = sum(wm :* dy_k) / sum_wm
                    mu_ipw = T_ipw - M_ipw
                    if (missing(T_ipw) | missing(M_ipw) | missing(mu_ipw)) {
                        _cdid_sel_link_fail("numerical", "Nonfinite IPW means.", 504)
                    }
                    R_ipw_m = (wp :* (dy_k :- T_ipw) / dp, wm :* (dy_k :- M_ipw) / dm, wp, wm)
                    V_ipw_m = J(4, p_trt, 0)
                    V_ipw_m[2, ] = quadcross(wm_act :* (dy_k :- M_ipw), X_trt) / (N * dm)
                    V_ipw_m[4, ] = m.trt_models[prior_tm].V_ipw[4, ]
                    wg_ipw = (rows(H_inv_m) == p_trt ? H_inv_m * V_ipw_m[2, ]' : J(0, 1, 0))
                }
            }
            else {
                m.trt_models[tm].rep_tid = tm
                if (can_cache) {
                    m.trt_models[tm].fit_u_cache = u_k
                    m.trt_models[tm].fit_q_cache = q_k
                    m.trt_models[tm].fit_q_sum   = q_sum
                }
                _cdid_sel_link_trt_logit(X_trt, y_trt, q_k, ctl.eps, ctl.maxit, beta_b = ., H_mat = .)
                H_mat  = H_mat / N
                H_chol = cholesky(H_mat)'
                if (missing(H_chol)) {
                    _cdid_sel_link_fail("numerical",
                        "Treatment propensity Fisher information is not positive definite.", 459)
                }
                H_inv_m = J(0, 0, 0)
                if (!ctl.points_only) {
                    z_sol_h = solvelower(H_chol', I(p_trt), rk_h1 = ., 0)
                    H_inv_m = solveupper(H_chol, z_sol_h, rk_h2 = ., 0)
                    if (rk_h1 != p_trt | rk_h2 != p_trt | missing(H_inv_m)) H_inv_m = J(0, 0, 0)
                }

                V_mat = _cdid_sel_link_logit_vals(X_trt * beta_b)
                prob = V_mat[, 1]; odds = V_mat[, 2]; odds_act = V_mat[, 3]
                wp     = q_k :* y_trt
                wm     = q_k :* (1 :- y_trt) :* odds
                sum_wp = sum(wp)
                sum_wm = sum(wm)
                dp     = sum_wp / N
                dm     = sum_wm / N
                if (dp <= 0 | dm <= 0 | missing(dp) | missing(dm)) {
                    _cdid_sel_link_fail("support",
                        "Treatment propensity normalizers (dp, dm) must be finite and positive.", 459)
                }

                T_ipw  = sum(wp :* dy_k) / sum_wp
                M_ipw  = sum(wm :* dy_k) / sum_wm
                mu_ipw = T_ipw - M_ipw
                if (missing(T_ipw) | missing(M_ipw) | missing(mu_ipw)) {
                    _cdid_sel_link_fail("numerical", "Nonfinite IPW means.", 504)
                }

                score_trt = (q_k :* (y_trt :- prob)) :* X_trt
                R_ipw_m   = (wp :* (dy_k :- T_ipw) / dp, wm :* (dy_k :- M_ipw) / dm, wp, wm)
                wm_act    = wm :* odds_act
                V_ipw_m   = J(4, p_trt, 0)
                V_ipw_m[2, ] = quadcross(wm_act :* (dy_k :- M_ipw), X_trt) / (N * dm)
                V_ipw_m[4, ] = quadcross(wm_act, X_trt) / N
                wg_ipw = (rows(H_inv_m) == p_trt ? H_inv_m * V_ipw_m[2, ]' : J(0, 1, 0))
            }

            m.trt_models[tm].n_risk       = rows(u_k)
            m.trt_models[tm].sum_wp       = sum_wp
            m.trt_models[tm].sum_wm       = sum_wm
            m.trt_models[tm].dp           = dp
            m.trt_models[tm].dm           = dm
            m.trt_models[tm].T_ipw        = T_ipw
            m.trt_models[tm].M_ipw        = M_ipw
            m.trt_models[tm].mu_ipw       = mu_ipw
            m.trt_models[tm].beta_trt     = beta_b
            m.trt_models[tm].H_chol       = H_chol
            m.trt_models[tm].H_inv        = H_inv_m
            m.trt_models[tm].V_ipw        = V_ipw_m
            m.trt_models[tm].w_gamma_ipw  = wg_ipw
            m.trt_models[tm].phi_risk_ipw = (!ctl.points_only & can_cache & rows(wg_ipw) == p_trt ?
                (R_ipw_m[, 1] :- R_ipw_m[, 2]) :- (score_trt * wg_ipw) : J(0, 1, 0))

            if (ctl.selection == "none" & !ctl.points_only) {
                m.trt_models[tm].baseline_C     = J(p_trt, rows(cohs_m), 0)
                m.trt_models[tm].baseline_U_ipw = J(4, rows(cohs_m), 0)
                for (k = 1; k <= rows(cohs_m); k = k + 1) {
                    sid = s_all[k]
                    if (!m.supports[sid].rate_active) continue
                    sub_k = _cdid_sel_link_index(coh_loc :== k)
                    if (rows(sub_k) > 0) {
                        scale_b = -1 / (N * m.supports[sid].rate_c)
                        m.trt_models[tm].baseline_C[, k] =
                            (colsum(score_trt[sub_k, ]))' * scale_b
                        m.trt_models[tm].baseline_U_ipw[, k] =
                            (colsum(R_ipw_m[sub_k, ]))' * scale_b
                    }
                }
                m.trt_models[tm].baseline_fast_ipw = (rows(wg_ipw) == p_trt ?
                    (m.trt_models[tm].baseline_U_ipw[1, ]' :- m.trt_models[tm].baseline_U_ipw[2, ]') :-
                    (m.trt_models[tm].baseline_C' * wg_ipw) : J(0, 1, 0))
            }

            R_dr_m = J(0, 0, 0)
            if (ctl.est_method == "dr") {
                L_b = m.or_blocks[b].L
                p_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].p_kept : p_kept)
                kept_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].kept : m.or_shared.kept)
                beta_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].beta : m.or_shared.beta)
                D_obs = (p_b > 0 ? L_b :* _cdid_sld_or_rows(m,r,u_k,kept_b) : J(rows(u_k), 0, 0))
                m_hat = (ctl.or_intercept ? m.or_blocks[b].alpha : 0) :+
                    (p_b > 0 ? D_obs * beta_b : J(rows(u_k), 1, 0))
                r_obs = dy_k :- m_hat
                T_dr  = sum(wp :* r_obs) / sum_wp
                M_dr  = sum(wm :* r_obs) / sum_wm
                mu_dr = T_dr - M_dr
                db_plus  = (p_b > 0 ? (quadcross(wp, D_obs) / sum_wp)' : J(0, 1, 0))
                db_minus = (p_b > 0 ? (quadcross(wm, D_obs) / sum_wm)' : J(0, 1, 0))
                db_diff  = db_plus - db_minus

                R_dr_m = (wp :* (r_obs :- T_dr) / dp, wm :* (r_obs :- M_dr) / dm, wp, wm)
                V_dr_m = J(4, p_trt, 0)
                V_dr_m[2, ] = quadcross(wm_act :* (r_obs :- M_dr), X_trt) / (N * dm)
                V_dr_m[4, ] = V_ipw_m[4, ]
                wg_dr = (rows(H_inv_m) == p_trt ? H_inv_m * V_dr_m[2, ]' : J(0, 1, 0))

                m.trt_models[tm].T_dr         = T_dr
                m.trt_models[tm].M_dr         = M_dr
                m.trt_models[tm].mu_dr        = mu_dr
                m.trt_models[tm].mu_dr_decomp = mu_ipw - (p_b > 0 ? sum(db_diff :* beta_b) : 0)
                m.trt_models[tm].d_bar_plus   = db_plus
                m.trt_models[tm].d_bar_minus  = db_minus
                m.trt_models[tm].d_bar_diff   = db_diff
                m.trt_models[tm].V_dr         = V_dr_m
                m.trt_models[tm].w_gamma_dr   = wg_dr
                m.trt_models[tm].phi_risk_dr  = (!ctl.points_only & can_cache & rows(wg_dr) == p_trt ?
                    (R_dr_m[, 1] :- R_dr_m[, 2]) :- (score_trt * wg_dr) : J(0, 1, 0))

                if (ctl.selection == "none" & !ctl.points_only) {
                    m.trt_models[tm].baseline_U_dr = J(4, rows(cohs_m), 0)
                    for (k = 1; k <= rows(cohs_m); k = k + 1) {
                        sid = s_all[k]
                        if (!m.supports[sid].rate_active) continue
                        sub_k = _cdid_sel_link_index(coh_loc :== k)
                        if (rows(sub_k) > 0) {
                            scale_b = -1 / (N * m.supports[sid].rate_c)
                            m.trt_models[tm].baseline_U_dr[, k] =
                                (colsum(R_dr_m[sub_k, ]))' * scale_b
                        }
                    }
                    m.trt_models[tm].baseline_fast_dr = (rows(wg_dr) == p_trt ?
                        (m.trt_models[tm].baseline_U_dr[1, ]' :- m.trt_models[tm].baseline_U_dr[2, ]') :-
                        (m.trt_models[tm].baseline_C' * wg_dr) : J(0, 1, 0))
                }
            }
            if (!ctl.points_only & can_cache) {
                m.trt_models[tm].risk_cache.loc          = u_k
                m.trt_models[tm].risk_cache.unit         = u_k
                m.trt_models[tm].risk_cache.dy           = dy_k
                m.trt_models[tm].risk_cache.q            = q_k
                m.trt_models[tm].risk_cache.cohort_local = coh_loc
                m.trt_models[tm].risk_cache.score        = score_trt
                m.trt_models[tm].risk_cache.R_ipw        = R_ipw_m
                m.trt_models[tm].risk_cache.R_dr         = R_dr_m
                m.trt_models[tm].risk_cache_ready        = 1
            }
        }
        for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
            m.trt_models[tm].fit_u_cache = J(0, 1, 0)
            m.trt_models[tm].fit_q_cache = J(0, 1, 0)
        }
    }

    /* 8. Populate Link Point Estimates */
    m.estimates  = J(n_links, 1, 0)
    m.mu_ipw_vec = J(n_links, 1, .)
    m.mu_or_vec  = J(n_links, 1, .)
    m.mu_dr_vec  = J(n_links, 1, .)

    for (l = 1; l <= n_links; l = l + 1) {
        if (ctl.est_method == "ipw") {
            tm = m.links[l].trt_id
            m.links[l].mu     = m.trt_models[tm].mu_ipw
            m.links[l].mu_ipw = m.trt_models[tm].mu_ipw
            m.estimates[l]    = m.links[l].mu
            m.mu_ipw_vec[l]   = m.links[l].mu_ipw
        }
        else if (ctl.est_method == "or") {
            m.links[l].mu    = m.links[l].mu_or
            m.estimates[l]   = m.links[l].mu_or
            m.mu_or_vec[l]   = m.links[l].mu_or
        }
        else {
            tm = m.links[l].trt_id
            m.links[l].mu     = m.trt_models[tm].mu_dr
            m.links[l].mu_dr  = m.trt_models[tm].mu_dr
            m.links[l].mu_ipw = m.trt_models[tm].mu_ipw
            m.estimates[l]    = m.links[l].mu_dr
            m.mu_dr_vec[l]    = m.links[l].mu_dr
            m.mu_ipw_vec[l]   = m.links[l].mu_ipw
            m.mu_or_vec[l]    = m.links[l].mu_or
        }
    }

    /* 9. Accumulate Upstream Selection Moment Sensitivities (unless points_only) */
    if (!ctl.points_only & ctl.selection == "conditional_mtar") {
        p_kept = m.or_shared.p_kept
        for (r = 1; r <= R; r = r + 1) {
            m_col = 0
            if (ctl.est_method == "or" | ctl.est_method == "dr") {
                if (ctl.outcome_pooling == "separate") {
                    for (b = 1; b <= m.n_blocks; b = b + 1) {
                        if (m.or_blocks[b].recipe_id == r) {
                            p_b = m.or_blocks[b].p_kept
                            if (p_b > 0) {
                                m.or_blocks[b].m_col_beta_start = m_col + 1
                                m.or_blocks[b].m_col_beta_end   = m_col + p_b
                                m_col = m_col + p_b
                            }
                            if (ctl.or_intercept) {
                                m_col = m_col + 1
                                m.or_blocks[b].m_col_alpha = m_col
                            }
                        }
                    }
                }
                else {
                    if (p_kept > 0) {
                        m.contexts[r].m_col_beta_start = m_col + 1
                        m.contexts[r].m_col_beta_end   = m_col + p_kept
                        m_col = m_col + p_kept
                    }
                    if (ctl.or_intercept) {
                        for (b = 1; b <= m.n_blocks; b = b + 1) {
                            if (m.or_blocks[b].recipe_id == r) {
                                m_col = m_col + 1
                                m.or_blocks[b].m_col_alpha = m_col
                            }
                        }
                    }
                }
                m_col = m_col + 1
                m.contexts[r].m_col_obj = m_col

                for (l = 1; l <= n_links; l = l + 1) {
                    if (m.links[l].recipe_id == r) {
                        p_b = (ctl.outcome_pooling == "separate" ?
                            m.or_blocks[m.links[l].block_id].p_kept : p_kept)
                        m_col = m_col + 1
                        m.links[l].m_col_or_mu = m_col
                        m_col = m_col + 1
                        m.links[l].m_col_or_den = m_col
                        if (p_b > 0) {
                            m.links[l].m_col_or_des_start = m_col + 1
                            m.links[l].m_col_or_des_end   = m_col + p_b
                            m_col = m_col + p_b
                        }
                    }
                }
            }

            if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
                for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
                    if (m.trt_models[tm].recipe_id == r) {
                        m.trt_models[tm].m_col_C_start = m_col + 1
                        m.trt_models[tm].m_col_C_end   = m_col + p_trt
                        m_col = m_col + p_trt
                        m.trt_models[tm].m_col_U_ipw_start = m_col + 1
                        m.trt_models[tm].m_col_U_ipw_end   = m_col + 4
                        m_col = m_col + 4
                        if (ctl.est_method == "dr") {
                            m.trt_models[tm].m_col_U_dr_start = m_col + 1
                            m.trt_models[tm].m_col_U_dr_end   = m_col + 4
                            m_col = m_col + 4
                        }
                    }
                }
            }

            m.contexts[r].n_moments_total = m_col
            k_idx=0
            if(ctl.est_method!="ipw") {
                for(b=1;b<=m.n_blocks;b++) if(m.or_blocks[b].recipe_id==r)
                    k_idx=k_idx+rows(m.or_blocks[b].support_ids)
                for(l=1;l<=n_links;l++) if(m.links[l].recipe_id==r) k_idx=k_idx+1
            }
            if(ctl.est_method!="or") for(tm=1;tm<=m.n_trt_models;tm++)
                if(m.trt_models[tm].recipe_id==r) k_idx=k_idx+rows(m.trt_models[tm].cohorts)
            m.contexts[r].loading_pieces = cdid_selected_loading_piece(k_idx,1)
            m.contexts[r].n_loading_pieces=0
        }

        /* Add OR control block moments (F_beta, forcing_alpha, forcing_objective) */
        if (ctl.est_method == "or" | ctl.est_method == "dr") {
            for (b = 1; b <= m.n_blocks; b = b + 1) {
                r   = m.or_blocks[b].recipe_id
                L_b = m.or_blocks[b].L
                s_all = m.or_blocks[b].support_ids
                p_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].p_kept : p_kept)
                kept_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].kept : m.or_shared.kept)
                beta_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].beta : m.or_shared.beta)
                b_c1 = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].m_col_beta_start : m.contexts[r].m_col_beta_start)
                b_c2 = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].m_col_beta_end : m.contexts[r].m_col_beta_end)
                if (!m.support_full_ready) {
                    if (b == 1) pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
                    else if (m.or_blocks[b].pair_id != m.or_blocks[b - 1].pair_id) {
                        pair_rows = _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, (1::N))
                    }
                }
                for (k = 1; k <= rows(s_all); k = k + 1) {
                    sid   = s_all[k]
                    if (m.support_full_ready) {
                        u_k = m.support_full[sid].unit; dy_k = m.support_full[sid].dy; q_k = m.support_full[sid].q
                    }
                    else {
                        read_rows = _cdid_sel_support_from_pair(m, sid, pair_rows)
                        u_k = read_rows.unit; dy_k = read_rows.dy; q_k = read_rows.q
                    }
                    if (p_b > 0) {
                        D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                            L_b :* m.supports[sid].Z_or[, kept_b] :
                            L_b :* _cdid_sld_or_rows(m,r,u_k,kept_b))
                    }
                    else D_kept = J(rows(u_k), 0, 0)
                    raw_res = dy_k :- m.or_blocks[b].alpha :-
                        (p_b > 0 ? D_kept * beta_b : J(rows(u_k), 1, 0))
                    centered = (ctl.or_intercept & p_b > 0 ? D_kept :- m.or_blocks[b].mean_d : D_kept)

                    moment_cols=J(0,1,0);V_mat=J(rows(u_k),0,0)
                    if(p_b>0) {
                        moment_cols=moment_cols\(b_c1::b_c2)
                        V_mat=V_mat,(centered:*(q_k:*raw_res))/N
                    }
                    if(ctl.or_intercept) {
                        moment_cols=moment_cols\m.or_blocks[b].m_col_alpha
                        V_mat=V_mat,(q_k:*raw_res)/N
                    }
                    moment_cols=moment_cols\m.contexts[r].m_col_obj
                    V_mat=V_mat,(q_k:*(raw_res:^2))/N
                    _cdid_selected_add_loading_sup(m,sid,u_k,moment_cols,V_mat)
                }
            }

            /* Add OR target support moments (mu_or, den, design) */
            for (l = 1; l <= n_links; l = l + 1) {
                r     = m.links[l].recipe_id
                b     = m.links[l].block_id
                sid   = m.links[l].target_support_id
                L_b   = m.links[l].L
                p_b   = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].p_kept : p_kept)
                kept_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].kept : m.or_shared.kept)
                beta_b = (ctl.outcome_pooling == "separate" ? m.or_blocks[b].beta : m.or_shared.beta)
                if (m.support_full_ready) {
                    u_k = m.support_full[sid].unit; dy_k = m.support_full[sid].dy; q_k = m.support_full[sid].q
                }
                else {
                    read_rows=_cdid_selected_support_rows(m,sid,(1::N))
                    u_k=read_rows.unit;dy_k=read_rows.dy;q_k=read_rows.q
                }
                if (p_b > 0) {
                    D_kept = (rows(m.supports[sid].Z_or) == rows(u_k) ?
                        L_b :* m.supports[sid].Z_or[, kept_b] :
                        L_b :* _cdid_sld_or_rows(m,r,u_k,kept_b))
                }
                else D_kept = J(rows(u_k), 0, 0)
                raw_res = dy_k :- m.or_blocks[b].alpha :-
                    (p_b > 0 ? D_kept * beta_b : J(rows(u_k), 1, 0))

                moment_cols=(m.links[l].m_col_or_mu\m.links[l].m_col_or_den)
                V_mat=((q_k:*(raw_res:-m.links[l].mu_or))/(N*m.links[l].d_target),q_k/N)
                if(p_b>0) {
                    moment_cols=moment_cols\(m.links[l].m_col_or_des_start::m.links[l].m_col_or_des_end)
                    V_mat=V_mat,(q_k:*(D_kept:-m.links[l].c_target))/N
                }
                _cdid_selected_add_loading_sup(m,sid,u_k,moment_cols,V_mat)
            }
        }

        /* Add Treatment Propensity & IPW/DR risk-set moments (C, U_ipw, U_dr) */
        if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
            for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
                r      = m.trt_models[tm].recipe_id
                cohs_m = m.trt_models[tm].cohorts
                s_all  = m.trt_models[tm].support_ids
                read_rows=(m.trt_models[tm].risk_cache_ready ? m.trt_models[tm].risk_cache :
                    _cdid_selected_treatment_rows(m,tm,(1::N),1))
                u_k=read_rows.unit;coh_loc=read_rows.cohort_local

                for (k = 1; k <= rows(cohs_m); k = k + 1) {
                    sub_k = _cdid_sel_link_index(coh_loc :== k)
                    if (rows(sub_k) == 0) continue
                    moment_cols=(m.trt_models[tm].m_col_C_start::m.trt_models[tm].m_col_U_ipw_end)
                    V_mat=(read_rows.score[sub_k,],read_rows.R_ipw[sub_k,])/N
                    if(ctl.est_method=="dr") {
                        moment_cols=moment_cols\(m.trt_models[tm].m_col_U_dr_start::m.trt_models[tm].m_col_U_dr_end)
                        V_mat=V_mat,read_rows.R_dr[sub_k,]/N
                    }
                    _cdid_selected_add_loading_sup(m,s_all[k],u_k[sub_k],moment_cols,V_mat)
                }
            }
        }

        /* Pre-assemble dense per-recipe moment accumulators + combined fast columns */
        if (can_cache) {
            for (r = 1; r <= R; r = r + 1) {
                m_col = m.contexts[r].n_moments_total
                m_fast = m_col
                if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
                    for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
                        if (m.trt_models[tm].recipe_id != r) continue
                        if (ctl.est_method == "ipw" & rows(m.trt_models[tm].w_gamma_ipw) == p_trt) {
                            m_fast = m_fast + 1
                            m.trt_models[tm].m_col_fast = m_fast
                        }
                        else if (ctl.est_method == "dr" & rows(m.trt_models[tm].w_gamma_dr) == p_trt) {
                            m_fast = m_fast + 1
                            m.trt_models[tm].m_col_fast = m_fast
                        }
                    }
                }
                acc_r = cdid_sel_moments_init(m.contexts[r].sel_fit, m_fast)
                for (k = 1; k <= m.contexts[r].n_loading_pieces; k = k + 1) {
                    if (rows(m.contexts[r].loading_pieces[k].theta_ids)) {
                        acc_r.G_theta[m.contexts[r].loading_pieces[k].theta_ids, m.contexts[r].loading_pieces[k].columns'] =
                            acc_r.G_theta[m.contexts[r].loading_pieces[k].theta_ids, m.contexts[r].loading_pieces[k].columns'] +
                            m.contexts[r].loading_pieces[k].G_theta
                    }
                    acc_r.G_ak[m.contexts[r].loading_pieces[k].cohort, m.contexts[r].loading_pieces[k].columns'] =
                        acc_r.G_ak[m.contexts[r].loading_pieces[k].cohort, m.contexts[r].loading_pieces[k].columns'] +
                        m.contexts[r].loading_pieces[k].G_ak
                }
                if (ctl.est_method == "ipw" | ctl.est_method == "dr") {
                    for (tm = 1; tm <= m.n_trt_models; tm = tm + 1) {
                        if (m.trt_models[tm].recipe_id != r | m.trt_models[tm].m_col_fast == 0) continue
                        m_fast = m.trt_models[tm].m_col_fast
                        b_c1 = m.trt_models[tm].m_col_C_start
                        b_c2 = m.trt_models[tm].m_col_C_end
                        if (ctl.est_method == "ipw") {
                            k_idx = m.trt_models[tm].m_col_U_ipw_start
                            acc_r.G_theta[, m_fast] = acc_r.G_theta[, k_idx] - acc_r.G_theta[, k_idx + 1] -
                                acc_r.G_theta[, b_c1..b_c2] * m.trt_models[tm].w_gamma_ipw
                            acc_r.G_ak[, m_fast] = acc_r.G_ak[, k_idx] - acc_r.G_ak[, k_idx + 1] -
                                acc_r.G_ak[, b_c1..b_c2] * m.trt_models[tm].w_gamma_ipw
                        }
                        else {
                            k_idx = m.trt_models[tm].m_col_U_dr_start
                            acc_r.G_theta[, m_fast] = acc_r.G_theta[, k_idx] - acc_r.G_theta[, k_idx + 1] -
                                acc_r.G_theta[, b_c1..b_c2] * m.trt_models[tm].w_gamma_dr
                            acc_r.G_ak[, m_fast] = acc_r.G_ak[, k_idx] - acc_r.G_ak[, k_idx + 1] -
                                acc_r.G_ak[, b_c1..b_c2] * m.trt_models[tm].w_gamma_dr
                        }
                    }
                }
                m.contexts[r].acc_full       = acc_r
                m.contexts[r].acc_full_ready = 1
            }
        }
    }

    m.prepared = 1
    return(m)
}

void _cdid_sel_slice_obs_vec(
    real colvector obs_u,
    real colvector units,
    real scalar N,
    real scalar lo_c,
    real scalar hi_c,
    real colvector loc_out,
    real colvector idx_out)
{
    real scalar n_obs, n_u, lo, hi, p1, p2
    real colvector ord_u, pos, hit

    n_obs = rows(obs_u)
    n_u   = rows(units)
    if (n_obs == 0 | n_u == 0) {
        loc_out = idx_out = J(0, 1, 0)
        return
    }
    if (lo_c > 0) {
        if (lo_c == 1 & hi_c == N) {
            loc_out = obs_u
            idx_out = (1::n_obs)
            return
        }
        if (obs_u[n_obs] < lo_c | obs_u[1] > hi_c) {
            loc_out = idx_out = J(0, 1, 0)
            return
        }
        p1 = _cdid_ipw_lower_bound(obs_u, lo_c)
        if (p1 > n_obs) {
            loc_out = idx_out = J(0, 1, 0)
            return
        }
        if (obs_u[p1] > hi_c) {
            loc_out = idx_out = J(0, 1, 0)
            return
        }
        p2 = _cdid_ipw_lower_bound(obs_u, hi_c + 1) - 1
        if (p2 < p1) {
            loc_out = idx_out = J(0, 1, 0)
            return
        }
        idx_out = (p1::p2)
        loc_out = obs_u[idx_out] :- (lo_c - 1)
        return
    }
    if (units[n_u] - units[1] == n_u - 1) {
        if ((n_u == 1 ? 1 : all(units[|2 \ n_u|] - units[|1 \ (n_u - 1)|] :== 1))) {
            lo = units[1]; hi = units[n_u]
            if (lo == 1 & hi == N) {
                loc_out = obs_u; idx_out = (1::n_obs); return
            }
            if (obs_u[n_obs] < lo | obs_u[1] > hi) {
                loc_out = idx_out = J(0, 1, 0); return
            }
            p1 = _cdid_ipw_lower_bound(obs_u, lo)
            if (p1 > n_obs) {
                loc_out = idx_out = J(0, 1, 0)
                return
            }
            if (obs_u[p1] > hi) {
                loc_out = idx_out = J(0, 1, 0)
                return
            }
            p2 = _cdid_ipw_lower_bound(obs_u, hi + 1) - 1
            if (p2 < p1) {
                loc_out = idx_out = J(0, 1, 0)
                return
            }
            idx_out = (p1::p2)
            loc_out = obs_u[idx_out] :- (lo - 1)
            return
        }
    }
    ord_u = order(units, 1)
    pos = _cdid_sel_link_match_sorted(units[ord_u], obs_u)
    hit = _cdid_sel_link_index(pos :> 0)
    if (rows(hit) == 0) {
        loc_out = idx_out = J(0, 1, 0)
        return
    }
    loc_out = ord_u[hit]
    idx_out = pos[hit]
}

void _cdid_selected_or_local_scores(struct cdid_selected_model scalar m,
    real colvector units, real colvector block_ids, real scalar want_objective,
    real matrix alpha_score, real colvector objective_score,
    | struct cdid_selected_rows colvector pair_cache)
{
    struct cdid_selected_rows scalar pair_rows,read_rows
    real scalar b,k,sid,r,Lb,p_b,is_sep,prev_pi,has_cache,lo_c,hi_c,n_u
    real rowvector kept_b
    real colvector selected,row_positions,u,q,dy,residual,beta_b,blk_map,loc_b,idx_b
    real matrix design
    n_u=rows(units)
    alpha_score=J(n_u,rows(block_ids),0)
    objective_score=J(n_u,1,0)
    is_sep=(m.outcome_pooling=="separate");prev_pi=0
    has_cache=(args()>=7 & rows(pair_cache)>0)
    blk_map=J(m.n_blocks,1,0)
    if(rows(block_ids)>0) blk_map[block_ids]=(1::rows(block_ids))
    lo_c=0;hi_c=0
    if(n_u>0) {
        if(units[n_u]-units[1]==n_u-1) {
            if((n_u==1 ? 1 : all(units[|2\n_u|]-units[|1\(n_u-1)|]:==1))) {
                lo_c=units[1]
                hi_c=units[n_u]
            }
        }
    }
    if(!want_objective & m.or_shared.intercept) {
        for(k=1;k<=rows(block_ids);k++) {
            b=block_ids[k]
            if(m.or_blocks[b].blk_cache_ready) {
                _cdid_sel_slice_obs_vec(m.or_blocks[b].blk_obs_u,units,m.N,lo_c,hi_c,loc_b=.,idx_b=.)
                if(rows(loc_b)) alpha_score[loc_b,k]=m.or_blocks[b].blk_ha[idx_b]
            }
            else {
                pair_rows=(has_cache ? pair_cache[m.or_blocks[b].pair_id] :
                    _cdid_selected_pair_rows(m,m.or_blocks[b].pair_id,units))
                if(rows(pair_rows.unit)==0) continue
                r=m.or_blocks[b].recipe_id;Lb=m.or_blocks[b].L
                p_b=(is_sep ? m.or_blocks[b].p_kept : m.or_shared.p_kept)
                kept_b=(is_sep ? m.or_blocks[b].kept : m.or_shared.kept)
                beta_b=(is_sep ? m.or_blocks[b].beta : m.or_shared.beta)
                for(sid=1;sid<=rows(m.or_blocks[b].support_ids);sid++) {
                    read_rows=_cdid_sel_support_from_pair(m,m.or_blocks[b].support_ids[sid],pair_rows)
                    row_positions=read_rows.loc
                    if(rows(row_positions)==0) continue
                    u=read_rows.unit;q=read_rows.q;dy=read_rows.dy
                    design=(p_b ? Lb*_cdid_sld_or_rows(m,r,u,kept_b) : J(rows(u),0,0))
                    residual=dy:-m.or_blocks[b].alpha:-(p_b ? design*beta_b : J(rows(u),1,0))
                    alpha_score[row_positions,k]=alpha_score[row_positions,k]+q:*residual
                }
            }
        }
        return
    }
    for(b=1;b<=m.n_blocks;b++) {
        if(blk_map[b]==0 & !want_objective) continue
        if(m.or_blocks[b].pair_id!=prev_pi) {
            pair_rows=(has_cache ? pair_cache[m.or_blocks[b].pair_id] :
                _cdid_selected_pair_rows(m,m.or_blocks[b].pair_id,units))
            prev_pi=m.or_blocks[b].pair_id
        }
        if(rows(pair_rows.unit)==0) continue
        r=m.or_blocks[b].recipe_id;Lb=m.or_blocks[b].L
        p_b=(is_sep ? m.or_blocks[b].p_kept : m.or_shared.p_kept)
        kept_b=(is_sep ? m.or_blocks[b].kept : m.or_shared.kept)
        beta_b=(is_sep ? m.or_blocks[b].beta : m.or_shared.beta)
        for(k=1;k<=rows(m.or_blocks[b].support_ids);k++) {
            sid=m.or_blocks[b].support_ids[k]
            read_rows=_cdid_sel_support_from_pair(m,sid,pair_rows)
            row_positions=read_rows.loc
            if(rows(row_positions)==0) continue
            u=read_rows.unit;q=read_rows.q;dy=read_rows.dy
            design=(p_b ? Lb*_cdid_sld_or_rows(m,r,u,kept_b) : J(rows(u),0,0))
            residual=dy:-m.or_blocks[b].alpha:-(p_b ? design*beta_b : J(rows(u),1,0))
            if(m.or_shared.intercept & blk_map[b]>0)
                alpha_score[row_positions,blk_map[b]]=alpha_score[row_positions,blk_map[b]]+q:*residual
            if(want_objective) objective_score[row_positions]=objective_score[row_positions]+q:*(residual:^2)
        }
    }
}

/* Explicit-array diagnostic adapter. Both entry points use the SAME fit. */
struct cdid_selected_model scalar cdid_selected_prepare_multi(
    struct cdid_selected_recipe_in colvector recipes,
    struct cdid_selected_pair_in colvector pairs,
    real colvector cohort_idx, real colvector cohort_values,
    real colvector calendar_values, real colvector is_A, real colvector w_pop,
    real matrix S_date_all, struct cdid_selected_link_spec colvector links,
    struct cdid_selected_controls scalar ctl, pointer() scalar population_ref)
{
    struct cdid_selected_registry scalar data
    struct cdid_selected_recipe_spec colvector specs
    struct cdid_population scalar pop
    real scalar R,r,ns,no,nt,nw,nr,N
    real colvector U,rho
    if(population_ref==NULL) _cdid_sel_link_fail("input","A borrowed compact population is required.",498)
    pop=*population_ref;N=pop.N;R=rows(recipes);specs=cdid_selected_recipe_spec(R,1)
    data.sel_data=cdid_selected_matrix_slot(R,1);data.or_data=cdid_selected_matrix_slot(R,1)
    data.trt_data=cdid_selected_matrix_slot(R,1);data.weight_data=cdid_selected_matrix_slot(R,1)
    data.rho_data=cdid_selected_matrix_slot(R,1);ns=no=nt=nw=nr=0
    if(rows(S_date_all)) {
        if(rows(S_date_all)!=N | cols(S_date_all)!=rows(calendar_values))
            _cdid_sel_link_fail("input","Diagnostic observation mask dimensions disagree with the population.",498)
        if(any(((rowsum(S_date_all:!=0):>=1):|is_A):!=pop.B))
            _cdid_sel_link_fail("input","Diagnostic B mask disagrees with the population.",498)
    }
    for(r=1;r<=R;r++) {
        specs[r].recipe_id=recipes[r].recipe_id;specs[r].recipe_key=recipes[r].recipe_key
        specs[r].design_cal=recipes[r].design_cal;specs[r].weight_cal=recipes[r].weight_cal
        specs[r].sel_has_intercept=recipes[r].sel_has_intercept;specs[r].sel_int_col=recipes[r].sel_int_col
        specs[r].has_xformla=recipes[r].has_xformla;specs[r].has_rho=recipes[r].has_rho
        specs[r].universe=recipes[r].universe;specs[r].justification=recipes[r].justification
        specs[r].sel_id=specs[r].or_id=specs[r].trt_id=specs[r].rho_id=0
        specs[r].weight_id=_cdid_sld_slot_add(data.weight_data,nw,recipes[r].weight_cal,recipes[r].w)
        if(ctl.selection=="conditional_mtar") {
            U=(recipes[r].universe=="B" ? pop.B : (recipes[r].universe=="A" ? is_A : J(N,1,1)))
            if(rows(recipes[r].B)==N) {
                if(any(recipes[r].B:!=pop.B)) _cdid_sel_link_fail("input","Diagnostic B mask disagrees with the population.",498)
            }
            if(rows(recipes[r].U)==N) {
                if(any(recipes[r].U:!=U)) _cdid_sel_link_fail("input","Diagnostic universe disagrees with the population.",498)
            }
            rho=(rows(recipes[r].rho)==N ? recipes[r].rho : J(N,1,1))
            specs[r].rho_id=_cdid_sld_slot_add(data.rho_data,nr,0,rho)
            specs[r].sel_id=_cdid_sld_slot_add(data.sel_data,ns,recipes[r].design_cal,recipes[r].Z_sel)
        }
        if(ctl.est_method!="ipw") specs[r].or_id=_cdid_sld_slot_add(data.or_data,no,recipes[r].design_cal,recipes[r].Z_or)
        if(ctl.est_method!="or") specs[r].trt_id=_cdid_sld_slot_add(data.trt_data,nt,recipes[r].design_cal,recipes[r].V_trt)
    }
    if(ns) data.sel_data=data.sel_data[1::ns];else data.sel_data=cdid_selected_matrix_slot(0,1)
    if(no) data.or_data=data.or_data[1::no];else data.or_data=cdid_selected_matrix_slot(0,1)
    if(nt) data.trt_data=data.trt_data[1::nt];else data.trt_data=cdid_selected_matrix_slot(0,1)
    if(nw) data.weight_data=data.weight_data[1::nw];else data.weight_data=cdid_selected_matrix_slot(0,1)
    if(nr) data.rho_data=data.rho_data[1::nr];else data.rho_data=cdid_selected_matrix_slot(0,1)
    return(cdid_selected_prepare_registry(specs,data,pairs,cohort_idx,cohort_values,
        calendar_values,is_A,w_pop,links,ctl,population_ref))
}

/* Single-recipe convenience entry point */
struct cdid_selected_model scalar cdid_selected_prepare(
    real matrix Z_sel,
    real scalar sel_has_intercept,
    real matrix Z_or,
    real matrix V_trt,
    real matrix Y_obs,
    real matrix obs_mat,
    real colvector g,
    real colvector w,
    real colvector U,
    real colvector rho,
    string scalar universe,
    string scalar justification,
    real colvector calendar_values,
    struct cdid_selected_link_spec colvector link_specs,
    struct cdid_selected_controls scalar ctl,
    | struct cdid_selected_pair_in colvector supplied_pairs)
{
    struct cdid_population scalar population
    struct cdid_population_ids scalar population_ids
    real colvector frame_g,long_id,long_date,long_g,long_y,missing_rows,canonical_g
    struct cdid_selected_recipe_in colvector recipes_in
    struct cdid_selected_pair_in colvector pairs_in
    real colvector cohort_values, cohort_idx, is_A
    real scalar N, J, i, k

    N = rows(g)
    if (rows(obs_mat) != N | cols(obs_mat) != rows(calendar_values)) {
        _cdid_sel_link_fail("input", "obs_mat dimensions must match (rows(g), rows(calendar_values)).", 459)
    }
    is_A = (rowsum(obs_mat :== 1) :>= 2)
    /* This table convenience has known cohorts only:0 or . means never.
       Unknown off-U metadata uses prepare_multi with the population known flag. */
    if(any(g:>.) | any(g:<0)) _cdid_sel_link_fail("input","Table cohort values must be nonnegative or known-never .; use the flagged population API for unknown cohorts.",498)
    canonical_g=g
    missing_rows=_cdid_sel_link_index(canonical_g:==0)
    if(rows(missing_rows)) canonical_g[missing_rows]=J(rows(missing_rows),1,.)
    cohort_values = uniqrows(sort(canonical_g, 1))
    J = rows(cohort_values)
    cohort_idx = J(N, 1, 0)
    for (k = 1; k <= J; k = k + 1) {
        missing_rows = _cdid_sel_link_index(canonical_g :== cohort_values[k])
        if (rows(missing_rows)) cohort_idx[missing_rows] = J(rows(missing_rows), 1, k)
    }

    recipes_in = cdid_selected_recipe_in(1, 1)
    recipes_in[1].recipe_id         = 1
    recipes_in[1].recipe_key        = "static"
    recipes_in[1].design_cal        = 0
    recipes_in[1].weight_cal        = 0
    recipes_in[1].Z_sel             = Z_sel
    recipes_in[1].sel_has_intercept = sel_has_intercept
    recipes_in[1].sel_int_col       = (sel_has_intercept ? 1 : 0)
    recipes_in[1].has_xformla       = (ctl.selection == "none" ? 0 : 1)
    recipes_in[1].has_rho           = (ctl.selection == "none" ? 0 : (universe != "F" | (rows(rho) == N & any(rho :!= 1))))
    recipes_in[1].Z_or              = Z_or
    recipes_in[1].V_trt             = V_trt
    recipes_in[1].w                 = w
    recipes_in[1].B                 = (rowsum(obs_mat :!= 0) :>= 1) :| is_A
    recipes_in[1].U                 = U
    recipes_in[1].rho               = rho
    recipes_in[1].universe          = universe
    recipes_in[1].justification     = justification

    if (args() >= 16 & rows(supplied_pairs) > 0) {
        pairs_in = supplied_pairs
    }
    else {
        pairs_in = cdid_selected_pairs_from_matrix(Y_obs, obs_mat, calendar_values, link_specs)
    }

    long_y=vec(Y_obs')
    missing_rows=_cdid_sel_link_index(vec(obs_mat'):!=1)
    if(rows(missing_rows)) long_y[missing_rows]=J(rows(missing_rows),1,.)
    frame_g=canonical_g
    missing_rows=_cdid_sel_link_index(frame_g:>=.)
    if(rows(missing_rows)) frame_g[missing_rows]=J(rows(missing_rows),1,0)
    long_id=(1::N)#J(rows(calendar_values),1,1)
    long_date=J(N,1,1)#calendar_values
    long_g=frame_g#J(rows(calendar_values),1,1)
    population_ids=cdid_population_ids_numeric(long_id,(1::N),1)
    population=cdid_population_prepare(population_ids,long_date,long_y,long_g,frame_g,1)
    return(cdid_selected_prepare_multi(recipes_in,pairs_in,cohort_idx,
        cohort_values,calendar_values,is_A,w,J(0,0,0),link_specs,ctl,&population))
}

/* --------------------------------------------------------------------------
   4. Reusable Nuisance & Conditional-CvM Moment Evaluators
   -------------------------------------------------------------------------- */

/* One baseline rate channel at a time, in full original-F units. The
   contraction never allocates a unit-by-all-support rate IF grid. */
real colvector _cdid_selected_rate_one(struct cdid_selected_model scalar m,
    real colvector unit_ids, real scalar s,
    | struct cdid_selected_rows colvector pair_cache)
{
    struct cdid_selected_rows scalar read_rows
    real colvector value,observed,ix,u_ix
    real scalar r,c,pi
    value=J(rows(unit_ids),1,0);observed=J(rows(unit_ids),1,0)
    pi=m.supports[s].pair_id
    read_rows=(args()>=4 & rows(pair_cache)>=pi ? pair_cache[pi] : _cdid_selected_pair_rows(m,pi,unit_ids))
    if(rows(read_rows.loc)) observed[read_rows.loc]=J(rows(read_rows.loc),1,1)
    r=m.supports[s].recipe_id;c=m.supports[s].cohort_pos
    ix=_cdid_sel_link_index(m.is_A[unit_ids]:&(m.cohort_idx[unit_ids]:==c):&_cdid_sld_keep(m,r,unit_ids))
    if(rows(ix)) {
        u_ix=unit_ids[ix]
        value[ix]=m.N:*_cdid_sld_w_rows(m,r,u_ix):*(observed[ix]:-m.supports[s].rate_baseline)/
            m.supports[s].repeat_mass_k
    }
    return(value)
}
real matrix _cdid_selected_rate_contract(struct cdid_selected_model scalar m,
    real colvector unit_ids, real colvector support_ids, real matrix loading,
    | struct cdid_selected_rows colvector pair_cache)
{
    real matrix value
    real scalar j,has_cache
    if(rows(support_ids)!=rows(loading))
        _cdid_sel_link_fail("input","Baseline rate loading dimensions disagree.",498)
    value=J(rows(unit_ids),cols(loading),0)
    has_cache=(args()>=5 & rows(pair_cache)>0)
    for(j=1;j<=rows(support_ids);j++) if(any(loading[j,]:!=0)) {
        if(has_cache) value=value+_cdid_selected_rate_one(m,unit_ids,support_ids[j],pair_cache)*loading[j,]
        else value=value+_cdid_selected_rate_one(m,unit_ids,support_ids[j])*loading[j,]
    }
    return(value)
}

void _cdid_sel_sep_block_if(struct cdid_selected_model scalar m, real scalar b,
    real colvector units, real scalar use_loaded,
    struct cdid_selected_loaded_block colvector loaded, real scalar need_alpha,
    real colvector ha_out, real matrix hb_out,
    real colvector xa_out, real matrix xb_out,
    | struct cdid_selected_rows colvector pair_cache,
      real scalar lo_in,
      real scalar hi_in)
{
    struct cdid_selected_rows scalar pair_rows, read_rows
    real scalar n_req, r, Lb, p_b, icpt, k, sid, c1, c2, has_cache, lo_c, hi_c
    real rowvector kept_b
    real colvector beta_b, row_pos, u, q, dy, res, u_mu, loc_b, idx_b
    real matrix design, centered, eval_b, rhs_b

    n_req  = rows(units)
    r      = m.or_blocks[b].recipe_id
    Lb     = m.or_blocks[b].L
    p_b    = m.or_blocks[b].p_kept
    kept_b = m.or_blocks[b].kept
    beta_b = m.or_blocks[b].beta
    icpt   = m.or_shared.intercept
    has_cache = (args() >= 11 & rows(pair_cache) > 0)
    lo_c   = (args() >= 13 & !missing(lo_in) ? lo_in : 0)
    hi_c   = (args() >= 13 & !missing(hi_in) ? hi_in : 0)

    ha_out = (need_alpha & icpt ? J(n_req, 1, 0) : J(0, 1, 0))
    hb_out = J(n_req, p_b, 0)
    if (m.or_blocks[b].blk_cache_ready) {
        if ((need_alpha & icpt) | p_b > 0) {
            _cdid_sel_slice_obs_vec(m.or_blocks[b].blk_obs_u, units, m.N, lo_c, hi_c, loc_b = ., idx_b = .)
            if (rows(loc_b) > 0) {
                if (need_alpha & icpt) ha_out[loc_b] = m.or_blocks[b].blk_ha[idx_b]
                if (p_b > 0) hb_out[loc_b, ] = m.or_blocks[b].blk_hb[idx_b, .]
            }
        }
    }
    else {
        if (!(need_alpha & icpt)) ha_out = J(n_req, 1, 0)
        pair_rows = (has_cache ? pair_cache[m.or_blocks[b].pair_id] :
            _cdid_selected_pair_rows(m, m.or_blocks[b].pair_id, units))
        if (rows(pair_rows.unit) > 0) {
            for (k = 1; k <= rows(m.or_blocks[b].support_ids); k = k + 1) {
                sid = m.or_blocks[b].support_ids[k]
                read_rows = _cdid_sel_support_from_pair(m, sid, pair_rows)
                row_pos = read_rows.loc
                if (rows(row_pos) == 0) continue
                u  = read_rows.unit
                q  = read_rows.q
                dy = read_rows.dy
                design = (p_b > 0 ? Lb :* _cdid_sld_or_rows(m, r, u, kept_b) : J(rows(u), 0, 0))
                res = dy :- m.or_blocks[b].alpha :- (p_b > 0 ? design * beta_b : J(rows(u), 1, 0))
                if (need_alpha & icpt) ha_out[row_pos] = ha_out[row_pos] + q :* res
                if (p_b > 0) {
                    centered = (icpt ? design :- m.or_blocks[b].mean_d : design)
                    hb_out[row_pos, ] = hb_out[row_pos, ] + centered :* (q :* res)
                }
            }
        }
    }

    xb_out = J(n_req, p_b, 0)
    if (p_b > 0) {
        rhs_b = hb_out :- m.or_blocks[b].score_mean_beta
        if (m.ctl.selection == "conditional_mtar") {
            c1 = m.or_blocks[b].m_col_beta_start
            c2 = m.or_blocks[b].m_col_beta_end
            eval_b = (use_loaded ? _cdid_selected_loaded(loaded, r, (c1::c2)) :
                _cdid_selected_loading(m, r, units, (c1::c2)))
            if (m.or_blocks[b].blk_cache_ready & rows(m.or_blocks[b].blk_gram_inv) == p_b) {
                xb_out = (rhs_b + eval_b) * m.or_blocks[b].blk_gram_inv
            }
            else {
                xb_out = (m.N :* _cdid_nuisance_gram(m.or_blocks[b].R, m.or_blocks[b].pivot, (rhs_b + eval_b)'))'
            }
        }
        else {
            if (m.or_blocks[b].blk_cache_ready & rows(m.or_blocks[b].blk_gram_inv) == p_b) {
                xb_out = rhs_b * m.or_blocks[b].blk_gram_inv +
                    _cdid_selected_rate_contract(m, units, m.or_blocks[b].support_ids, m.or_blocks[b].baseline_beta_rate', pair_cache)
            }
            else {
                xb_out = (m.N :* _cdid_nuisance_gram(m.or_blocks[b].R, m.or_blocks[b].pivot, rhs_b'))' +
                    _cdid_selected_rate_contract(m, units, m.or_blocks[b].support_ids, m.or_blocks[b].baseline_beta_rate', pair_cache)
            }
        }
    }

    xa_out = J(n_req, 1, 0)
    if (need_alpha & icpt) {
        if (m.ctl.selection == "conditional_mtar") {
            u_mu = (use_loaded ? _cdid_selected_loaded(loaded, r, m.or_blocks[b].m_col_alpha) :
                _cdid_selected_loading(m, r, units, m.or_blocks[b].m_col_alpha))
        }
        else {
            u_mu = _cdid_selected_rate_contract(m, units, m.or_blocks[b].support_ids,
                m.or_blocks[b].baseline_alpha_rate, pair_cache) / m.N
        }
        xa_out = (ha_out :- m.or_blocks[b].score_mean_alpha + u_mu) / m.or_blocks[b].mass -
            (p_b > 0 ? xb_out * m.or_blocks[b].mean_d' : J(n_req, 1, 0))
    }
}

/* Evaluate untreated-change OR nuisance IFs (xi_beta, xi_alpha, objective_if)
   on unit_ids. Reusable by actual-change imputation, augmented double-DR, and
   conditional CvM consumers. */
void cdid_selected_or_eval_chunk(struct cdid_selected_model scalar m,
    real colvector unit_ids, real matrix xi_beta_out,
    real matrix xi_alpha_out, real colvector obj_if_out,
    | real colvector requested_blocks)
{
    struct cdid_selected_loaded_block colvector empty_loads
    real matrix h_beta,eval_beta,alpha_score,hb_tmp,xb_tmp
    real colvector blocks,objective_score,u_mu,ha_tmp,xa_tmp
    real scalar p_kept,r,b,j,tot_p,at
    if(!m.prepared | m.ctl.points_only)
        _cdid_sel_link_fail("replay","Model is not prepared for IF replay.",459)
    if(m.ctl.est_method=="ipw")
        _cdid_sel_link_fail("replay","OR IF requested from a pure IPW model.",459)
    blocks=(args()<6 ? (1::m.n_blocks) : requested_blocks)
    if(m.outcome_pooling=="separate") {
        _cdid_selected_or_local_scores(m,unit_ids,J(0,1,0),1,alpha_score=.,objective_score=.)
        obj_if_out=objective_score:-m.or_shared.objective
        if(m.ctl.selection=="conditional_mtar") {
            for(r=1;r<=m.R;r++)
                obj_if_out=obj_if_out+_cdid_selected_loading(m,r,unit_ids,m.contexts[r].m_col_obj)
        } else {
            obj_if_out=obj_if_out+_cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),
                m.or_shared.baseline_objective_rate_loading)
        }
        empty_loads=cdid_selected_loaded_block(0,1)
        tot_p=0
        for(j=1;j<=rows(blocks);j++) tot_p=tot_p+m.or_blocks[blocks[j]].p_kept
        xi_beta_out=J(rows(unit_ids),tot_p,0)
        xi_alpha_out=J(rows(unit_ids),rows(blocks),0)
        at=0
        for(j=1;j<=rows(blocks);j++) {
            b=blocks[j]
            _cdid_sel_sep_block_if(m,b,unit_ids,0,empty_loads,m.or_shared.intercept,
                ha_tmp=.,hb_tmp=.,xa_tmp=.,xb_tmp=.)
            if(m.or_blocks[b].p_kept>0) {
                xi_beta_out[,(at+1)..(at+m.or_blocks[b].p_kept)]=xb_tmp
                at=at+m.or_blocks[b].p_kept
            }
            if(m.or_shared.intercept) xi_alpha_out[,j]=xa_tmp
        }
        return
    }
    p_kept=m.or_shared.p_kept
    h_beta=(p_kept ? m.or_shared.beta_score_F[unit_ids,] : J(rows(unit_ids),0,0))
    _cdid_selected_or_local_scores(m,unit_ids,blocks,1,alpha_score,objective_score)
    xi_beta_out=J(rows(unit_ids),p_kept,0)
    xi_alpha_out=J(rows(unit_ids),rows(blocks),0)
    obj_if_out=objective_score:-m.or_shared.objective
    if(m.ctl.selection=="conditional_mtar") {
        eval_beta=J(rows(unit_ids),p_kept,0)
        for(r=1;r<=m.R;r++) {
            if(p_kept) eval_beta=eval_beta+_cdid_selected_loading(m,r,unit_ids,
                (m.contexts[r].m_col_beta_start::m.contexts[r].m_col_beta_end))
            obj_if_out=obj_if_out+_cdid_selected_loading(m,r,unit_ids,m.contexts[r].m_col_obj)
        }
        if(p_kept) {
            if (m.or_shared.shared_gram_ready & rows(m.or_shared.shared_gram_inv) == p_kept)
                xi_beta_out = (h_beta + eval_beta) * m.or_shared.shared_gram_inv
            else xi_beta_out=(m.N:*_cdid_nuisance_gram(m.or_shared.R,m.or_shared.pivot,(h_beta+eval_beta)'))'
        }
    } else {
        if(p_kept) {
            if (m.or_shared.shared_gram_ready & rows(m.or_shared.shared_gram_inv) == p_kept)
                xi_beta_out = h_beta * m.or_shared.shared_gram_inv +
                    _cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),m.or_shared.baseline_beta_rate_loading')
            else xi_beta_out=(m.N:*_cdid_nuisance_gram(m.or_shared.R,m.or_shared.pivot,h_beta'))'+
                _cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),m.or_shared.baseline_beta_rate_loading')
        }
        obj_if_out=obj_if_out+_cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),
            m.or_shared.baseline_objective_rate_loading)
    }
    if(m.or_shared.intercept) for(j=1;j<=rows(blocks);j++) {
        b=blocks[j];r=m.or_blocks[b].recipe_id
        if(m.ctl.selection=="conditional_mtar")
            u_mu=_cdid_selected_loading(m,r,unit_ids,m.or_blocks[b].m_col_alpha)
        else u_mu=_cdid_selected_rate_contract(m,unit_ids,m.or_blocks[b].support_ids,
            m.or_blocks[b].baseline_alpha_rate)/m.N
        xi_alpha_out[,j]=(alpha_score[,j]+u_mu)/m.or_blocks[b].mass
        if(p_kept) xi_alpha_out[,j]=xi_alpha_out[,j]-xi_beta_out*m.or_blocks[b].mean_d'
    }
}

/* Evaluate treatment propensity parameter IF xi_gamma (n_chunk x p_trt) on
   unit_ids for treatment model trt_id, Reusable by conditional CvM / custom consumers. */
real matrix cdid_selected_ps_eval_chunk(struct cdid_selected_model scalar m,
    real scalar trt_id, real colvector unit_ids)
{
    struct cdid_selected_rows scalar read_rows
    struct cdid_sel_moments scalar acc_c
    real matrix beta_score, C_xi_u, z_sol, xi_gamma
    real colvector pos_u, loc_u
    real scalar n_req, r, p_trt, c1, c2, rk

    if (!m.prepared | m.ctl.points_only) {
        _cdid_sel_link_fail("replay",
            "Model is not prepared for IF replay or was fitted in points_only mode.", 459)
    }
    if (trt_id < 1 | trt_id > m.n_trt_models | trt_id != floor(trt_id)) {
        _cdid_sel_link_fail("input", "Invalid trt_id in cdid_selected_ps_eval_chunk.", 459)
    }
    n_req = rows(unit_ids)
    p_trt = m.trt_models[trt_id].p_trt
    r     = m.trt_models[trt_id].recipe_id

    beta_score = J(n_req, p_trt, 0)
    read_rows=_cdid_selected_treatment_rows(m,trt_id,unit_ids,1)
    loc_u=read_rows.loc
    if(rows(loc_u)) beta_score[loc_u,]=read_rows.score


    if (m.ctl.selection == "conditional_mtar") {
        c1 = m.trt_models[trt_id].m_col_C_start
        c2 = m.trt_models[trt_id].m_col_C_end
        C_xi_u=_cdid_selected_loading(m,r,unit_ids,(c1::c2))
    }
    else {
        C_xi_u=_cdid_selected_rate_contract(m,unit_ids,m.trt_models[trt_id].support_ids,
            m.trt_models[trt_id].baseline_C')
    }

    if (rows(m.trt_models[trt_id].H_inv) == p_trt) {
        xi_gamma = (beta_score + C_xi_u) * m.trt_models[trt_id].H_inv
        if (missing(xi_gamma)) {
            _cdid_sel_link_fail("numerical", "Nonfinite triangular solve in cdid_selected_ps_eval_chunk.", 504)
        }
        return(xi_gamma)
    }
    z_sol = solvelower(m.trt_models[trt_id].H_chol', (beta_score + C_xi_u)', rk = ., 0)
    xi_gamma = (solveupper(m.trt_models[trt_id].H_chol, z_sol, rk = ., 0))'
    if (rk != p_trt | missing(xi_gamma)) {
        _cdid_sel_link_fail("numerical", "Nonfinite triangular solve in cdid_selected_ps_eval_chunk.", 504)
    }
    return(xi_gamma)
}

/* Generic moment IF evaluator on risk-set observers of treatment model trt_id
   (for conditional CvM / custom residual moment extractors):
   obs_moment_w is n_risk x n_m in canonical original-unit risk order (unnormalized q_i*m_i),
   d_gamma is n_m x p_trt derivative wrt treatment propensity slope gamma. */
real matrix cdid_selected_ps_moment_if_chunk(struct cdid_selected_model scalar m,
    real scalar trt_id, real matrix obs_moment_w, real matrix d_gamma,
    real colvector unit_ids)
{
    struct cdid_selected_rows scalar risk_rows
    struct cdid_sel_moments scalar acc_m
    real matrix direct_m, sel_m, xi_gamma
    real colvector pos_u, loc_u, cohs_m, pc_all, u_k, coh_loc, sub_k
    real rowvector mean_m
    real scalar n_req, n_m, p_trt, r, k, sid

    if (!m.prepared | m.ctl.points_only) {
        _cdid_sel_link_fail("replay",
            "Model is not prepared for IF replay or was fitted in points_only mode.", 459)
    }
    n_req = rows(unit_ids)
    n_m   = cols(obs_moment_w)
    p_trt = m.trt_models[trt_id].p_trt
    r     = m.trt_models[trt_id].recipe_id
    if (rows(obs_moment_w) != m.trt_models[trt_id].n_risk |
        rows(d_gamma) != n_m | cols(d_gamma) != p_trt) {
        _cdid_sel_link_fail("input", "Dimension mismatch in cdid_selected_ps_moment_if_chunk.", 459)
    }

    risk_rows=_cdid_selected_treatment_rows(m,trt_id,(1::m.N),0)
    mean_m = colsum(obs_moment_w) / m.N
    direct_m = J(n_req, 1, -1) * mean_m
    pos_u = _cdid_sel_link_match_sorted(unit_ids, risk_rows.unit)
    loc_u = _cdid_sel_link_index(pos_u :> 0)
    if (rows(loc_u) > 0) {
        direct_m[loc_u, ] = direct_m[loc_u, ] + obs_moment_w[pos_u[loc_u], ]
    }

    if (m.ctl.selection == "conditional_mtar") {
        _cdid_selected_activate_reader(m,r)
        acc_m = cdid_sel_moments_init(m.contexts[r].sel_fit, n_m)
        cohs_m  = m.trt_models[trt_id].cohorts
        pc_all  = m.trt_models[trt_id].pc_ids
        u_k     = risk_rows.unit
        coh_loc = risk_rows.cohort_local
        for (k = 1; k <= rows(cohs_m); k = k + 1) {
            sub_k = _cdid_sel_link_index(coh_loc :== k)
            if (rows(sub_k) > 0) {
                cdid_sel_add_moments(acc_m, m.contexts[r].sel_fit, m.reader,
                    m.contexts[r].up_map, pc_all[k], u_k[sub_k], obs_moment_w[sub_k, ] / m.N)
            }
        }
        sel_m = cdid_sel_eval_loading(m.contexts[r].sel_fit, m.reader, acc_m, unit_ids)
    }
    else {
        sel_m = J(n_req, n_m, 0)
        for(k=1;k<=rows(m.trt_models[trt_id].support_ids);k++) {
            sid=m.trt_models[trt_id].support_ids[k]
            if(!m.supports[sid].rate_active) continue
            sub_k=_cdid_sel_link_index(risk_rows.cohort_local:==k)
            if(rows(sub_k)) sel_m=sel_m+_cdid_selected_rate_one(m,unit_ids,sid)*
                (-colsum(obs_moment_w[sub_k,])/(m.N*m.supports[sid].rate_c))
        }
    }

    xi_gamma = cdid_selected_ps_eval_chunk(m, trt_id, unit_ids)
    return(direct_m + sel_m + xi_gamma * d_gamma')
}

/* --------------------------------------------------------------------------
   5. Complete Original-F Chunk IF Evaluator
   -------------------------------------------------------------------------- */

struct cdid_selected_chunk_out scalar cdid_selected_chunk(
    struct cdid_selected_model scalar m,
    real colvector unit_ids,
    | real colvector req_links,
      real scalar scale_core,
      real scalar want_details)
{
    struct cdid_selected_rows scalar read_rows
    struct cdid_selected_rows colvector pair_cache
    struct cdid_selected_chunk_out scalar out
    struct cdid_selected_loaded_block colvector selected_loads
    struct cdid_sel_sep_blk_chunk colvector sep_blks
    real matrix beta_force, alpha_score, hb_tmp, xb_tmp, dummy_a
    real colvector objective_force, objective_score, block_ids, b_position, ha_tmp, xa_tmp, beta_b
    real colvector link_ids, tids_req, pos_u, loc_u, u_mu, u_den, a_if, b_if, r_offset
    real colvector psi_dir_m, psi_or_m, psi_ps_m, psi_sel_m, psi_dr_m, psi_ipw_m
    real matrix h_beta, eval_beta, rhs_beta, xi_beta, xi_alpha
    real matrix direct_mu, direct_den, direct_des_j, u_des, des_if_j
    real matrix beta_score, C_xi_u, U_ipw_xi_u, U_dr_xi_u, direct_ipw, direct_dr
    real matrix z_sol, xi_gamma, coord_ipw_m, coord_dr_m, norm_A_m
    real matrix D_loc, res_loc, phi_trt, D_diff_mat, C_target_mat
    real colvector q_loc, dy_loc, eval_obj, blk_pos_map, trt_pos_map, req_trt_pos, idx_u
    real scalar n_req, n_l_req, n_t_req, p_kept, p_trt, r, b, j, l, tm, k, sid, L_l, d_l, mu_l
    real scalar c1, c2, rk, ws_chunk, work_chunk, is_sep, p_bound, p_b, need_a, tot_p, tot_des, at, des_at
    real scalar lo_c, hi_c, need_pair_cache, fast_trt
    real rowvector c_l, kept_b

    if (!m.prepared) {
        _cdid_sel_link_fail("replay", "Model is not prepared.", 459)
    }
    if (m.ctl.points_only) {
        _cdid_sel_link_fail("replay",
            "IF replay requested on points-only prepared model; never refit by chunk.", 459)
    }

    n_req = rows(unit_ids)
    if (n_req < 1 | cols(unit_ids) != 1 | missing(unit_ids) |
        any(unit_ids :< 1) | any(unit_ids :> m.N) | any(unit_ids :!= floor(unit_ids))) {
        _cdid_sel_link_fail("input", "unit_ids must be valid 1..N integer positions.", 459)
    }
    if (rows(uniqrows(unit_ids)) != n_req) {
        _cdid_sel_link_fail("input", "unit_ids must not contain duplicates.", 459)
    }
    if (m.ctl.max_chunk_rows > 0 & n_req > m.ctl.max_chunk_rows) {
        _cdid_sel_link_fail("input", "Chunk size exceeds max_chunk_rows.", 459)
    }

    if (args() < 3 | rows(req_links) == 0) {
        link_ids = (1::m.n_links)
    }
    else {
        link_ids = req_links
        if (cols(link_ids) != 1 | missing(link_ids) |
            any(link_ids :< 1) | any(link_ids :> m.n_links) |
            any(link_ids :!= floor(link_ids)) |
            rows(uniqrows(link_ids)) != rows(link_ids)) {
            _cdid_sel_link_fail("input", "link_ids must be distinct integers in 1..n_links.", 459)
        }
    }
    if (args() < 4 | missing(scale_core)) scale_core = 0
    if (args() < 5 | missing(want_details)) want_details = 1
    if (scale_core != 0 & scale_core != 1) {
        _cdid_sel_link_fail("input", "scale_core must be 0 or 1.", 459)
    }
    if (want_details != 0 & want_details != 1) {
        _cdid_sel_link_fail("input", "want_details must be 0 or 1.", 459)
    }

    n_l_req = rows(link_ids)
    is_sep  = (m.outcome_pooling == "separate")
    p_kept  = m.or_shared.p_kept
    p_bound = (is_sep ? m.or_shared.p_raw : p_kept)
    p_trt   = (m.n_trt_models > 0 ? m.trt_models[1].p_trt : 0)

    ws_chunk = m.workspace_doubles + n_req * (8 * n_l_req + 4 * p_bound + 4 * p_trt + 2 * m.n_cohorts + 32)
    work_chunk = m.work_ops + 64 * n_req * (n_l_req + m.R) * (p_bound + p_trt + 8)
    _cdid_sel_link_budget_check(ws_chunk, work_chunk, m.ctl.max_work, m.ctl.max_workspace_doubles, "cdid_selected_chunk")

    out.alpha_block_ids=J(0,1,0)
    out.unit_ids            = unit_ids
    out.link_ids            = link_ids
    out.scale_core          = scale_core
    out.normalization_count = (scale_core ? m.n_A : m.N)
    out.phi                 = J(n_req, n_l_req, 0)
    out.psi_dir             = J(0, 0, 0)
    out.psi_or              = J(0, 0, 0)
    out.psi_ps              = J(0, 0, 0)
    out.psi_sel             = J(0, 0, 0)
    out.psi_ipw             = J(0, 0, 0)
    out.psi_or_estimator    = J(0, 0, 0)
    out.xi_beta             = J(0, 0, 0)
    out.beta_score          = J(0, 0, 0)
    out.xi_alpha            = J(0, 0, 0)
    out.alpha_score         = J(0, 0, 0)
    out.den_if              = J(0, 0, 0)
    out.num_if              = J(0, 0, 0)
    out.design_if           = J(0, 0, 0)
    out.objective_if        = J(0, 1, 0)
    out.coord_ipw           = J(0, 0, 0)
    out.coord_dr            = J(0, 0, 0)
    out.xi_gamma            = J(0, 0, 0)
    out.norm_A_if           = J(0, 0, 0)
    out.tau_if              = J(0, 0, 0)
    out.eta_if              = J(0, 0, 0)
    out.lambda_if           = J(0, 1, 0)

    lo_c = 0; hi_c = 0
    if (unit_ids[n_req] - unit_ids[1] == n_req - 1) {
        if ((n_req == 1 ? 1 : all(unit_ids[|2 \ n_req|] - unit_ids[|1 \ (n_req - 1)|] :== 1))) {
            lo_c = unit_ids[1]; hi_c = unit_ids[n_req]
        }
    }

    selected_loads = cdid_selected_loaded_block(0, 1)
    if(m.ctl.selection=="conditional_mtar")
        selected_loads=_cdid_selected_load_requests(m,unit_ids,link_ids,want_details,beta_force,objective_force)
    need_pair_cache = (want_details == 1 | m.ctl.selection == "none" | !m.support_full_ready)
    if (need_pair_cache) {
        pair_cache = cdid_selected_rows(m.n_pairs, 1)
        for (k = 1; k <= m.n_pairs; k = k + 1) {
            pair_cache[k] = _cdid_selected_pair_rows(m, k, unit_ids)
        }
    }
    else {
        pair_cache = cdid_selected_rows(0, 1)
    }

    /* 2. Evaluate OR Nuisance and Link IFs (when OR or DR is active) */
    if (m.ctl.est_method == "or" | m.ctl.est_method == "dr") {
        if (is_sep) {
            block_ids = J(n_l_req, 1, 0)
            for (j = 1; j <= n_l_req; j = j + 1) block_ids[j] = m.links[link_ids[j]].block_id
            block_ids = uniqrows(sort(block_ids, 1))
            blk_pos_map = J(m.n_blocks, 1, 0)
            blk_pos_map[block_ids] = (1::rows(block_ids))
            out.alpha_block_ids = block_ids
            sep_blks = cdid_sel_sep_blk_chunk(rows(block_ids), 1)
            need_a = ((m.ctl.est_method == "or" | want_details) & m.or_shared.intercept ? 1 : 0)
            tot_p = 0
            for (k = 1; k <= rows(block_ids); k = k + 1) {
                b = block_ids[k]
                _cdid_sel_sep_block_if(m, b, unit_ids, 1, selected_loads, need_a,
                    ha_tmp = ., hb_tmp = ., xa_tmp = ., xb_tmp = ., pair_cache, lo_c, hi_c)
                sep_blks[k].ha = ha_tmp
                sep_blks[k].hb = hb_tmp
                sep_blks[k].xa = xa_tmp
                sep_blks[k].xb = xb_tmp
                tot_p = tot_p + m.or_blocks[b].p_kept
            }
            out.xi_beta = J(n_req, tot_p, 0)
            if (want_details) out.beta_score = J(n_req, tot_p, 0)
            at = 0
            for (k = 1; k <= rows(block_ids); k = k + 1) {
                p_b = m.or_blocks[block_ids[k]].p_kept
                if (p_b > 0) {
                    out.xi_beta[, (at + 1)..(at + p_b)] = sep_blks[k].xb
                    if (want_details) out.beta_score[, (at + 1)..(at + p_b)] = sep_blks[k].hb
                    at = at + p_b
                }
            }

            if (m.ctl.est_method == "or" & !want_details) {
                for (j = 1; j <= n_l_req; j = j + 1) {
                    l   = link_ids[j]
                    r   = m.links[l].recipe_id
                    b   = m.links[l].block_id
                    sid = m.links[l].target_support_id
                    c_l = m.links[l].c_target
                    p_b = m.or_blocks[b].p_kept
                    if (m.links[l].target_mu_ready & m.support_full_ready) {
                        _cdid_sel_slice_obs_vec(m.support_full[sid].unit, unit_ids, m.N, lo_c, hi_c, loc_u = ., idx_u = .)
                        if (rows(loc_u) > 0) out.phi[loc_u, j] = m.links[l].target_mu_vec[idx_u]
                    }
                    else {
                        L_l = m.links[l].L; d_l = m.links[l].d_target; mu_l = m.links[l].mu_or
                        read_rows = (rows(pair_cache) > 0 ?
                            _cdid_sel_support_from_pair(m, sid, pair_cache[m.supports[sid].pair_id]) :
                            _cdid_selected_support_rows(m, sid, unit_ids))
                        loc_u = read_rows.loc
                        if (rows(loc_u) > 0) {
                            D_loc = (p_b > 0 ? L_l :* _cdid_sld_or_rows(m, r, unit_ids[loc_u], m.or_blocks[b].kept) : J(rows(loc_u), 0, 0))
                            res_loc = read_rows.dy :- (m.or_shared.intercept ? m.or_blocks[b].alpha : 0) :-
                                (p_b > 0 ? D_loc * m.or_blocks[b].beta : J(rows(loc_u), 1, 0))
                            out.phi[loc_u, j] = read_rows.q :* (res_loc :- mu_l) / d_l
                        }
                    }
                    if (m.ctl.selection == "conditional_mtar") {
                        u_mu = _cdid_selected_loaded(selected_loads, r, m.links[l].m_col_or_mu)
                    }
                    else {
                        u_mu = _cdid_selected_rate_one(m, unit_ids, sid, pair_cache) :* m.links[l].baseline_or_rate_loading
                    }
                    k = blk_pos_map[b]
                    a_if = (m.or_shared.intercept ? sep_blks[k].xa : J(n_req, 1, 0))
                    b_if = (p_b > 0 ? sep_blks[k].xb * c_l' : J(n_req, 1, 0))
                    out.phi[, j] = out.phi[, j] - a_if - b_if + u_mu
                }
            }
            else if (m.ctl.est_method == "or" | want_details) {
                xi_alpha    = J(n_req, rows(block_ids), 0)
                alpha_score = J(n_req, rows(block_ids), 0)
                if (m.or_shared.intercept) {
                    for (k = 1; k <= rows(block_ids); k = k + 1) {
                        xi_alpha[, k]    = sep_blks[k].xa
                        alpha_score[, k] = sep_blks[k].ha
                    }
                }
                out.xi_alpha    = xi_alpha
                out.alpha_score = alpha_score
                _cdid_selected_or_local_scores(m, unit_ids, J(0, 1, 0), 1, dummy_a = ., objective_score = ., pair_cache)
                out.objective_if = objective_score :- m.or_shared.objective
                if (m.ctl.selection == "conditional_mtar") out.objective_if = out.objective_if + objective_force
                else out.objective_if = out.objective_if + _cdid_selected_rate_contract(m, unit_ids,
                    (1::m.n_supports), m.or_shared.baseline_objective_rate_loading, pair_cache)

                tot_des = 0
                for (j = 1; j <= n_l_req; j = j + 1) tot_des = tot_des + m.or_blocks[m.links[link_ids[j]].block_id].p_kept
                direct_mu  = J(n_req, n_l_req, 0)
                direct_den = J(n_req, n_l_req, 0)
                out.psi_or_estimator = J(n_req, n_l_req, 0)
                out.den_if           = J(n_req, n_l_req, 0)
                out.num_if           = J(n_req, n_l_req, 0)
                out.design_if        = J(n_req, tot_des, 0)
                des_at = 0

                for (j = 1; j <= n_l_req; j = j + 1) {
                    l    = link_ids[j]
                    r    = m.links[l].recipe_id
                    b    = m.links[l].block_id
                    sid  = m.links[l].target_support_id
                    L_l  = m.links[l].L
                    d_l  = m.links[l].d_target
                    mu_l = m.links[l].mu_or
                    c_l  = m.links[l].c_target
                    p_b    = m.or_blocks[b].p_kept
                    kept_b = m.or_blocks[b].kept
                    beta_b = m.or_blocks[b].beta

                    direct_den[, j] = J(n_req, 1, -d_l)
                    direct_des_j = J(n_req, p_b, 0)

                    read_rows = _cdid_sel_support_from_pair(m, sid, pair_cache[m.supports[sid].pair_id])
                    loc_u = read_rows.loc
                    if (rows(loc_u) > 0) {
                        q_loc  = read_rows.q
                        dy_loc = read_rows.dy
                        D_loc  = (p_b > 0 ? L_l :* _cdid_sld_or_rows(m, r, unit_ids[loc_u], kept_b) : J(rows(loc_u), 0, 0))
                        res_loc = dy_loc :- (m.or_shared.intercept ? m.or_blocks[b].alpha : 0) :-
                            (p_b > 0 ? D_loc * beta_b : J(rows(loc_u), 1, 0))
                        direct_mu[loc_u, j]  = q_loc :* (res_loc :- mu_l) / d_l
                        direct_den[loc_u, j] = direct_den[loc_u, j] + q_loc
                        if (p_b > 0) {
                            direct_des_j[loc_u, ] = q_loc :* (D_loc :- c_l)
                        }
                    }

                    if (m.ctl.selection == "conditional_mtar") {
                        u_mu  = _cdid_selected_loaded(selected_loads, r, m.links[l].m_col_or_mu)
                        u_den = _cdid_selected_loaded(selected_loads, r, m.links[l].m_col_or_den)
                        if (p_b > 0) {
                            c1 = m.links[l].m_col_or_des_start
                            c2 = m.links[l].m_col_or_des_end
                            u_des = _cdid_selected_loaded(selected_loads, r, (c1::c2))
                        }
                        else {
                            u_des = J(n_req, 0, 0)
                        }
                    }
                    else {
                        u_mu  = _cdid_selected_rate_one(m, unit_ids, sid, pair_cache)
                        u_den = u_mu :* m.links[l].baseline_den_rate_loading
                        u_des = u_mu * m.links[l].baseline_design_rate_loading
                        u_mu  = u_mu :* m.links[l].baseline_or_rate_loading
                    }

                    k = blk_pos_map[b]
                    a_if = (m.or_shared.intercept ? sep_blks[k].xa : J(n_req, 1, 0))
                    b_if = (p_b > 0 ? sep_blks[k].xb * c_l' : J(n_req, 1, 0))
                    out.psi_or_estimator[, j] = direct_mu[, j] - a_if - b_if + u_mu
                    out.den_if[, j]           = direct_den[, j] + u_den
                    out.num_if[, j]           = d_l :* (direct_mu[, j] + u_mu) +
                        mu_l :* out.den_if[, j] - d_l :* a_if - d_l :* b_if
                    if (p_b > 0) {
                        des_if_j = direct_des_j + u_des + out.den_if[, j] * c_l
                        out.design_if[, (des_at + 1)..(des_at + p_b)] = des_if_j
                        des_at = des_at + p_b
                    }
                }
                if (m.ctl.est_method == "or") {
                    out.phi = out.psi_or_estimator
                }
            }
        }
        else {
            h_beta = (p_kept > 0 ? m.or_shared.beta_score_F[unit_ids, ] : J(n_req, 0, 0))
            xi_beta = J(n_req, p_kept, 0)
            if (p_kept > 0) {
                if (m.ctl.selection == "conditional_mtar") {
                    eval_beta=beta_force
                    rhs_beta = h_beta + eval_beta
                    if (m.or_shared.shared_gram_ready & rows(m.or_shared.shared_gram_inv) == p_kept)
                        xi_beta = rhs_beta * m.or_shared.shared_gram_inv
                    else xi_beta = (m.N :* _cdid_nuisance_gram(m.or_shared.R, m.or_shared.pivot, rhs_beta'))'
                }
                else {
                    if (m.or_shared.shared_gram_ready & rows(m.or_shared.shared_gram_inv) == p_kept)
                        xi_beta = h_beta * m.or_shared.shared_gram_inv +
                            _cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),m.or_shared.baseline_beta_rate_loading',pair_cache)
                    else xi_beta = (m.N :* _cdid_nuisance_gram(m.or_shared.R, m.or_shared.pivot, h_beta'))' +
                        _cdid_selected_rate_contract(m,unit_ids,(1::m.n_supports),m.or_shared.baseline_beta_rate_loading',pair_cache)
                }
            }
            out.xi_beta = xi_beta
            if (want_details) out.beta_score = h_beta

            if (m.ctl.est_method == "or" & !want_details) {
                block_ids = J(n_l_req, 1, 0)
                for (j = 1; j <= n_l_req; j = j + 1) block_ids[j] = m.links[link_ids[j]].block_id
                block_ids = uniqrows(sort(block_ids, 1))
                blk_pos_map = J(m.n_blocks, 1, 0)
                blk_pos_map[block_ids] = (1::rows(block_ids))
                out.alpha_block_ids = block_ids
                xi_alpha = J(n_req, rows(block_ids), 0)
                if (m.or_shared.intercept) {
                    _cdid_selected_or_local_scores(m, unit_ids, block_ids, 0, alpha_score = ., objective_score = ., pair_cache)
                    for (k = 1; k <= rows(block_ids); k = k + 1) {
                        b = block_ids[k]; r = m.or_blocks[b].recipe_id
                        if (m.ctl.selection == "conditional_mtar")
                            u_mu = _cdid_selected_loaded(selected_loads, r, m.or_blocks[b].m_col_alpha)
                        else u_mu = _cdid_selected_rate_contract(m, unit_ids, m.or_blocks[b].support_ids,
                            m.or_blocks[b].baseline_alpha_rate, pair_cache) / m.N
                        xi_alpha[, k] = (alpha_score[, k] + u_mu) / m.or_blocks[b].mass
                        if (p_kept) xi_alpha[, k] = xi_alpha[, k] - xi_beta * m.or_blocks[b].mean_d'
                    }
                }
                C_target_mat = (p_kept > 0 ? J(p_kept, n_l_req, 0) : J(0, 0, 0))
                for (j = 1; j <= n_l_req; j = j + 1) {
                    l   = link_ids[j]
                    r   = m.links[l].recipe_id
                    b   = m.links[l].block_id
                    sid = m.links[l].target_support_id
                    if (m.links[l].target_mu_ready & m.support_full_ready) {
                        _cdid_sel_slice_obs_vec(m.support_full[sid].unit, unit_ids, m.N, lo_c, hi_c, loc_u = ., idx_u = .)
                        if (rows(loc_u) > 0) out.phi[loc_u, j] = m.links[l].target_mu_vec[idx_u]
                    }
                    else {
                        L_l = m.links[l].L; d_l = m.links[l].d_target; mu_l = m.links[l].mu_or
                        read_rows = (rows(pair_cache) > 0 ?
                            _cdid_sel_support_from_pair(m, sid, pair_cache[m.supports[sid].pair_id]) :
                            _cdid_selected_support_rows(m, sid, unit_ids))
                        loc_u = read_rows.loc
                        if (rows(loc_u) > 0) {
                            D_loc = (p_kept > 0 ? L_l :* _cdid_sld_or_rows(m, r, unit_ids[loc_u], m.or_shared.kept) : J(rows(loc_u), 0, 0))
                            res_loc = read_rows.dy :- (m.or_shared.intercept ? m.or_blocks[b].alpha : 0) :-
                                (p_kept > 0 ? D_loc * m.or_shared.beta : J(rows(loc_u), 1, 0))
                            out.phi[loc_u, j] = read_rows.q :* (res_loc :- mu_l) / d_l
                        }
                    }
                    if (m.ctl.selection == "conditional_mtar")
                        u_mu = _cdid_selected_loaded(selected_loads, r, m.links[l].m_col_or_mu)
                    else u_mu = _cdid_selected_rate_one(m, unit_ids, sid, pair_cache) :* m.links[l].baseline_or_rate_loading
                    if (m.or_shared.intercept) out.phi[, j] = out.phi[, j] - xi_alpha[, blk_pos_map[b]] + u_mu
                    else out.phi[, j] = out.phi[, j] + u_mu
                    if (p_kept > 0) C_target_mat[, j] = m.links[l].c_target'
                }
                if (p_kept > 0) out.phi = out.phi - xi_beta * C_target_mat
            }
            else if (m.ctl.est_method == "or" | want_details) {
                block_ids=J(n_l_req,1,0)
                for(j=1;j<=n_l_req;j++) block_ids[j]=m.links[link_ids[j]].block_id
                block_ids=uniqrows(sort(block_ids,1));out.alpha_block_ids=block_ids
                blk_pos_map = J(m.n_blocks, 1, 0)
                blk_pos_map[block_ids] = (1::rows(block_ids))
                _cdid_selected_or_local_scores(m,unit_ids,block_ids,1,alpha_score,objective_score,pair_cache)
                xi_alpha=J(n_req,rows(block_ids),0)
                if(m.or_shared.intercept) for(k=1;k<=rows(block_ids);k++) {
                    b=block_ids[k];r=m.or_blocks[b].recipe_id
                    if(m.ctl.selection=="conditional_mtar")
                        u_mu=_cdid_selected_loaded(selected_loads,r,m.or_blocks[b].m_col_alpha)
                    else u_mu=_cdid_selected_rate_contract(m,unit_ids,m.or_blocks[b].support_ids,
                m.or_blocks[b].baseline_alpha_rate,pair_cache)/m.N
                    xi_alpha[,k]=(alpha_score[,k]+u_mu)/m.or_blocks[b].mass
                    if(p_kept) xi_alpha[,k]=xi_alpha[,k]-xi_beta*m.or_blocks[b].mean_d'
                }
                out.xi_alpha=xi_alpha;out.alpha_score=alpha_score
                out.objective_if=objective_score:-m.or_shared.objective
                if(m.ctl.selection=="conditional_mtar") out.objective_if=out.objective_if+objective_force
                else out.objective_if=out.objective_if+_cdid_selected_rate_contract(m,unit_ids,
                    (1::m.n_supports),m.or_shared.baseline_objective_rate_loading,pair_cache)

                direct_mu  = J(n_req, n_l_req, 0)
                direct_den = J(n_req, n_l_req, 0)
                out.psi_or_estimator = J(n_req, n_l_req, 0)
                out.den_if           = J(n_req, n_l_req, 0)
                out.num_if           = J(n_req, n_l_req, 0)
                out.design_if        = J(n_req, n_l_req * p_kept, 0)

                for (j = 1; j <= n_l_req; j = j + 1) {
                    l   = link_ids[j]
                    r   = m.links[l].recipe_id
                    b   = m.links[l].block_id
                    sid = m.links[l].target_support_id
                    L_l = m.links[l].L
                    d_l = m.links[l].d_target
                    mu_l = m.links[l].mu_or
                    c_l  = m.links[l].c_target

                    direct_den[, j] = J(n_req, 1, -d_l)
                    direct_des_j = J(n_req, p_kept, 0)

                    read_rows=_cdid_sel_support_from_pair(m,sid,pair_cache[m.supports[sid].pair_id])
                    loc_u=read_rows.loc
                    if(rows(loc_u)>0) {
                        q_loc=read_rows.q;dy_loc=read_rows.dy
                        D_loc  = (p_kept > 0 ? L_l :* _cdid_sld_or_rows(m,r,unit_ids[loc_u],m.or_shared.kept) : J(rows(loc_u), 0, 0))
                        res_loc = dy_loc :- (m.or_shared.intercept ? m.or_blocks[b].alpha : 0) :-
                            (p_kept > 0 ? D_loc * m.or_shared.beta : J(rows(loc_u), 1, 0))
                        direct_mu[loc_u, j]  = q_loc :* (res_loc :- mu_l) / d_l
                        direct_den[loc_u, j] = direct_den[loc_u, j] + q_loc
                        if (p_kept > 0) {
                            direct_des_j[loc_u, ] = q_loc :* (D_loc :- c_l)
                        }
                    }

                    if (m.ctl.selection == "conditional_mtar") {
                        u_mu=_cdid_selected_loaded(selected_loads,r,m.links[l].m_col_or_mu)
                        u_den=_cdid_selected_loaded(selected_loads,r,m.links[l].m_col_or_den)
                        if (p_kept > 0) {
                            c1=m.links[l].m_col_or_des_start;c2=m.links[l].m_col_or_des_end
                            u_des=_cdid_selected_loaded(selected_loads,r,(c1::c2))
                        }
                        else {
                            u_des = J(n_req, 0, 0)
                        }
                    }
                    else {
                        u_mu=_cdid_selected_rate_one(m,unit_ids,sid,pair_cache)
                        u_den=u_mu:*m.links[l].baseline_den_rate_loading
                        u_des=u_mu*m.links[l].baseline_design_rate_loading
                        u_mu=u_mu:*m.links[l].baseline_or_rate_loading
                    }

                    k = blk_pos_map[b]
                    a_if=(m.or_shared.intercept ? xi_alpha[,k] : J(n_req,1,0))
                    b_if = (p_kept > 0 ? xi_beta * c_l' : J(n_req, 1, 0))
                    out.psi_or_estimator[, j] = direct_mu[, j] - a_if - b_if + u_mu
                    out.den_if[, j]           = direct_den[, j] + u_den
                    out.num_if[, j]           = d_l :* (direct_mu[, j] + u_mu) +
                        mu_l :* out.den_if[, j] - d_l :* a_if - d_l :* b_if
                    if (p_kept > 0) {
                        des_if_j = direct_des_j + u_des + out.den_if[, j] * c_l
                        out.design_if[, ((j - 1) * p_kept + 1)..(j * p_kept)] = des_if_j
                    }
                }
                if (m.ctl.est_method == "or") {
                    out.phi = out.psi_or_estimator
                }
            }
        }
    }

    /* 3. Evaluate Treatment Propensity, IPW, and DR IFs (when IPW or DR is active) */
    if (m.ctl.est_method == "ipw" | m.ctl.est_method == "dr") {
        req_trt_pos = J(n_l_req, 1, 0)
        for (j = 1; j <= n_l_req; j = j + 1) {
            req_trt_pos[j] = m.links[link_ids[j]].trt_id
        }
        tids_req = uniqrows(sort(req_trt_pos, 1))
        n_t_req  = rows(tids_req)
        trt_pos_map = J(m.n_trt_models, 1, 0)
        trt_pos_map[tids_req] = (1::n_t_req)
        req_trt_pos = trt_pos_map[req_trt_pos]

        fast_trt = (!want_details)
        if (fast_trt) {
            for (k = 1; k <= n_t_req; k = k + 1) {
                tm = tids_req[k]
                if (!m.trt_models[tm].risk_cache_ready) {
                    fast_trt = 0
                    break
                }
                if (m.ctl.est_method == "ipw" & rows(m.trt_models[tm].phi_risk_ipw) != m.trt_models[tm].n_risk) {
                    fast_trt = 0; break
                }
                if (m.ctl.est_method == "dr" & rows(m.trt_models[tm].phi_risk_dr) != m.trt_models[tm].n_risk) {
                    fast_trt = 0; break
                }
                if (m.ctl.selection == "conditional_mtar") {
                    r = m.trt_models[tm].recipe_id
                    if (!m.contexts[r].acc_full_ready | m.trt_models[tm].m_col_fast == 0) {
                        fast_trt = 0; break
                    }
                }
            }
        }

        if (fast_trt) {
            phi_trt = J(n_req, n_t_req, 0)
            D_diff_mat = (m.ctl.est_method == "dr" & !is_sep & p_kept > 0 ? J(p_kept, n_t_req, 0) : J(0, 0, 0))
            for (k = 1; k <= n_t_req; k = k + 1) {
                tm = tids_req[k]
                r  = m.trt_models[tm].recipe_id
                _cdid_sel_slice_obs_vec(m.trt_models[tm].risk_cache.unit, unit_ids, m.N, lo_c, hi_c, loc_u = ., idx_u = .)
                if (rows(loc_u) > 0) {
                    if (m.ctl.est_method == "ipw") phi_trt[loc_u, k] = m.trt_models[tm].phi_risk_ipw[idx_u]
                    else phi_trt[loc_u, k] = m.trt_models[tm].phi_risk_dr[idx_u]
                }
                if (m.ctl.selection == "conditional_mtar") {
                    phi_trt[, k] = phi_trt[, k] + _cdid_selected_loaded(selected_loads, r, m.trt_models[tm].m_col_fast)
                }
                else {
                    if (m.ctl.est_method == "ipw") {
                        phi_trt[, k] = phi_trt[, k] + _cdid_selected_rate_contract(m, unit_ids,
                            m.trt_models[tm].support_ids, m.trt_models[tm].baseline_fast_ipw, pair_cache)
                    }
                    else {
                        phi_trt[, k] = phi_trt[, k] + _cdid_selected_rate_contract(m, unit_ids,
                            m.trt_models[tm].support_ids, m.trt_models[tm].baseline_fast_dr, pair_cache)
                    }
                }
                if (m.ctl.est_method == "dr") {
                    if (is_sep) {
                        b = m.trt_models[tm].block_id
                        if (m.or_blocks[b].p_kept > 0) {
                            phi_trt[, k] = phi_trt[, k] - sep_blks[blk_pos_map[b]].xb * m.trt_models[tm].d_bar_diff
                        }
                    }
                    else if (p_kept > 0) {
                        D_diff_mat[, k] = m.trt_models[tm].d_bar_diff
                    }
                }
            }
            if (m.ctl.est_method == "dr" & !is_sep & p_kept > 0) {
                phi_trt = phi_trt - xi_beta * D_diff_mat
            }
            out.phi = phi_trt[, req_trt_pos']
        }
        else {
            if (want_details) {
                out.psi_ipw   = J(n_req, n_l_req, 0)
                out.coord_ipw = J(n_req, 4 * n_l_req, 0)
                out.xi_gamma  = J(n_req, p_trt * n_l_req, 0)
                out.norm_A_if = J(n_req, 2 * n_l_req, 0)
                if (m.ctl.est_method == "dr") {
                    out.psi_dir  = J(n_req, n_l_req, 0)
                    out.psi_or   = J(n_req, n_l_req, 0)
                    out.psi_ps   = J(n_req, n_l_req, 0)
                    out.psi_sel  = J(n_req, n_l_req, 0)
                    out.coord_dr = J(n_req, 4 * n_l_req, 0)
                }
            }

            for (k = 1; k <= n_t_req; k = k + 1) {
                tm = tids_req[k]
                r  = m.trt_models[tm].recipe_id

                beta_score = J(n_req, p_trt, 0)
                direct_ipw = J(n_req, 4, 0)
                direct_ipw[, 3] = J(n_req, 1, -m.trt_models[tm].dp)
                direct_ipw[, 4] = J(n_req, 1, -m.trt_models[tm].dm)
                if (m.ctl.est_method == "dr") {
                    direct_dr = J(n_req, 4, 0)
                    direct_dr[, 3] = J(n_req, 1, -m.trt_models[tm].dp)
                    direct_dr[, 4] = J(n_req, 1, -m.trt_models[tm].dm)
                }

                read_rows=(rows(pair_cache) > 0 ?
                    _cdid_selected_treatment_rows(m,tm,unit_ids,1,pair_cache) :
                    _cdid_selected_treatment_rows(m,tm,unit_ids,1))
                loc_u=read_rows.loc
                if(rows(loc_u)>0) {
                    beta_score[loc_u,]=read_rows.score
                    direct_ipw[loc_u,1..2]=read_rows.R_ipw[,1..2]
                    direct_ipw[loc_u,3..4]=direct_ipw[loc_u,3..4]+read_rows.R_ipw[,3..4]
                    if(m.ctl.est_method=="dr") {
                        direct_dr[loc_u,1..2]=read_rows.R_dr[,1..2]
                        direct_dr[loc_u,3..4]=direct_dr[loc_u,3..4]+read_rows.R_dr[,3..4]
                    }
                }

                if (m.ctl.selection == "conditional_mtar") {
                    c1 = m.trt_models[tm].m_col_C_start
                    c2 = m.trt_models[tm].m_col_C_end
                    C_xi_u = _cdid_selected_loaded(selected_loads,r,(c1::c2))
                    if (want_details) {
                        c1 = m.trt_models[tm].m_col_U_ipw_start
                        c2 = m.trt_models[tm].m_col_U_ipw_end
                        U_ipw_xi_u = _cdid_selected_loaded(selected_loads,r,(c1::c2))
                        if (m.ctl.est_method == "dr") {
                            c1 = m.trt_models[tm].m_col_U_dr_start
                            c2 = m.trt_models[tm].m_col_U_dr_end
                            U_dr_xi_u = _cdid_selected_loaded(selected_loads,r,(c1::c2))
                        }
                    }
                    else {
                        if (m.ctl.est_method == "dr") {
                            c1 = m.trt_models[tm].m_col_U_dr_start
                            U_dr_xi_u = (_cdid_selected_loaded(selected_loads,r,(c1 \ (c1 + 1))), J(n_req, 2, 0))
                            U_ipw_xi_u = J(n_req, 4, 0)
                        }
                        else {
                            c1 = m.trt_models[tm].m_col_U_ipw_start
                            U_ipw_xi_u = (_cdid_selected_loaded(selected_loads,r,(c1 \ (c1 + 1))), J(n_req, 2, 0))
                        }
                    }
                }
                else {
                    C_xi_u=_cdid_selected_rate_contract(m,unit_ids,m.trt_models[tm].support_ids,m.trt_models[tm].baseline_C',pair_cache)
                    U_ipw_xi_u=_cdid_selected_rate_contract(m,unit_ids,m.trt_models[tm].support_ids,m.trt_models[tm].baseline_U_ipw',pair_cache)
                    if (m.ctl.est_method == "dr") {
                        U_dr_xi_u=_cdid_selected_rate_contract(m,unit_ids,m.trt_models[tm].support_ids,m.trt_models[tm].baseline_U_dr',pair_cache)
                    }
                }

                if (rows(m.trt_models[tm].H_inv) == p_trt) {
                    xi_gamma = (beta_score + C_xi_u) * m.trt_models[tm].H_inv
                    rk = p_trt
                }
                else {
                    z_sol = solvelower(m.trt_models[tm].H_chol', (beta_score + C_xi_u)', rk = ., 0)
                    xi_gamma = (solveupper(m.trt_models[tm].H_chol, z_sol, rk = ., 0))'
                }
                if (rk != p_trt | missing(xi_gamma)) {
                    _cdid_sel_link_fail("numerical", "Nonfinite triangular solve in treatment IF replay.", 504)
                }

                coord_ipw_m = direct_ipw + U_ipw_xi_u + xi_gamma * m.trt_models[tm].V_ipw'
                psi_ipw_m   = coord_ipw_m[, 1] - coord_ipw_m[, 2]
                if (want_details) {
                    norm_A_m = (coord_ipw_m[, 3..4] :/ m.lambda) :-
                        ((m.is_A[unit_ids] :- m.lambda) * (m.trt_models[tm].dp, m.trt_models[tm].dm) :/ (m.lambda^2))
                }

                if (m.ctl.est_method == "dr") {
                    coord_dr_m = direct_dr + U_dr_xi_u + xi_gamma * m.trt_models[tm].V_dr'
                    psi_dir_m  = direct_dr[, 1] - direct_dr[, 2]
                    if (is_sep) {
                        b = m.trt_models[tm].block_id
                        p_b = m.or_blocks[b].p_kept
                        psi_or_m = (p_b > 0 ? -(sep_blks[blk_pos_map[b]].xb * m.trt_models[tm].d_bar_diff) : J(n_req, 1, 0))
                    }
                    else {
                        psi_or_m = (p_kept > 0 ? -(xi_beta * m.trt_models[tm].d_bar_diff) : J(n_req, 1, 0))
                    }
                    psi_ps_m   = -(xi_gamma * m.trt_models[tm].V_dr[2, ]')
                    psi_sel_m  = U_dr_xi_u[, 1] - U_dr_xi_u[, 2]
                    psi_dr_m   = psi_dir_m + psi_or_m + psi_ps_m + psi_sel_m
                    if (want_details) {
                        norm_A_m = (coord_dr_m[, 3..4] :/ m.lambda) :-
                            ((m.is_A[unit_ids] :- m.lambda) * (m.trt_models[tm].dp, m.trt_models[tm].dm) :/ (m.lambda^2))
                    }
                }

                for (j = 1; j <= n_l_req; j = j + 1) {
                    if (req_trt_pos[j] != k) continue
                    if (want_details) {
                        out.psi_ipw[, j] = psi_ipw_m
                        out.coord_ipw[, ((j - 1) * 4 + 1)..(j * 4)] = coord_ipw_m
                        out.xi_gamma[, ((j - 1) * p_trt + 1)..(j * p_trt)] = xi_gamma
                        out.norm_A_if[, ((j - 1) * 2 + 1)..(j * 2)] = norm_A_m
                    }
                    if (m.ctl.est_method == "ipw") {
                        out.phi[, j] = psi_ipw_m
                    }
                    else {
                        out.phi[, j] = psi_dr_m
                        if (want_details) {
                            out.psi_dir[, j] = psi_dir_m
                            out.psi_or[, j]  = psi_or_m
                            out.psi_ps[, j]  = psi_ps_m
                            out.psi_sel[, j] = psi_sel_m
                            out.coord_dr[, ((j - 1) * 4 + 1)..(j * 4)] = coord_dr_m
                        }
                    }
                }
            }
        }
    }

    /* 4. Population Share IF Channels (tau, eta, lambda) */
    if (want_details) {
        out.tau_if    = J(n_req, m.n_cohorts, 0)
        out.eta_if    = J(n_req, m.n_cohorts, 0)
        out.lambda_if = m.is_A[unit_ids] :- m.lambda
        for (k = 1; k <= m.n_cohorts; k = k + 1) {
            out.tau_if[, k] = m.w_pop[unit_ids] :*
                ((m.cohort_idx[unit_ids] :== k) :- m.tau[k]) :/ m.mean_w
            out.eta_if[, k] = m.w_pop[unit_ids] :* m.is_A[unit_ids] :*
                ((m.cohort_idx[unit_ids] :== k) :- m.eta[k]) :/ m.mean_wA
        }
    }

    if (missing(out.phi)) {
        _cdid_sel_link_fail("numerical", "Nonfinite values in complete link IF chunk.", 504)
    }
    if (scale_core) {
        out.phi = m.lambda :* out.phi
    }
    return(out)
}

/* --------------------------------------------------------------------------
   6. Streamed Covariance and Full-F Link IF Assembly
   -------------------------------------------------------------------------- */

/* Streamed unit and cluster covariance without materializing full N x L phi */
void cdid_selected_cov_stream(
    struct cdid_selected_model scalar m,
    real matrix V_unit,
    real matrix Omega_unit,
    | real scalar chunk_rows,
      real colvector cluster_id,
      real matrix V_cluster,
      real matrix Omega_cluster)
{
    struct cdid_selected_chunk_out scalar ch
    real colvector u_chunk, cl_vals, ix_c, c_uniq
    real matrix G_cl, info_cl
    real scalar step, lo, hi, n_cl

    if (!m.prepared | m.ctl.points_only) {
        _cdid_sel_link_fail("replay",
            "Model is not prepared for IF replay or was fitted in points_only mode.", 459)
    }
    step = (args() >= 4 & !missing(chunk_rows) & chunk_rows > 0 ? chunk_rows :
           (m.ctl.max_chunk_rows > 0 ? m.ctl.max_chunk_rows : m.N))
    if (step < 1 | step != floor(step)) {
        _cdid_sel_link_fail("input", "chunk_rows must be a positive integer.", 459)
    }

    Omega_unit = J(m.n_links, m.n_links, 0)
    n_cl = 0
    if (args() >= 5 & rows(cluster_id) == m.N) {
        ix_c = order(cluster_id, 1)
        cl_vals = cluster_id[ix_c]
        u_chunk = J(m.N, 1, 0)
        u_chunk[ix_c] = runningsum(m.N > 1 ? (1 \ (cl_vals[2..m.N] :!= cl_vals[1..(m.N - 1)])) : 1)
        cl_vals = u_chunk
        n_cl = cl_vals[ix_c[m.N]]
        G_cl = J(n_cl, m.n_links, 0)
    }

    for (lo = 1; lo <= m.N; lo = lo + step) {
        hi = min((m.N, lo + step - 1))
        u_chunk = (lo::hi)
        ch = cdid_selected_chunk(m, u_chunk, (1::m.n_links), 0, 0)
        Omega_unit = Omega_unit + quadcross(ch.phi, ch.phi)
        if (n_cl > 0) {
            ix_c = cl_vals[u_chunk]
            if (rows(u_chunk) > 1) {
                u_chunk = order((ix_c, (1::rows(u_chunk))), (1, 2))
                info_cl = _cdid_ps_slink(ix_c[u_chunk], 1)
                c_uniq  = ix_c[u_chunk[info_cl[, 1]]]
                G_cl[c_uniq, ] = G_cl[c_uniq, ] + panelsum(ch.phi[u_chunk, ], info_cl)
            }
            else {
                G_cl[ix_c[1], ] = G_cl[ix_c[1], ] + ch.phi
            }
        }
    }
    Omega_unit = Omega_unit / m.N
    V_unit     = Omega_unit / m.N
    if (n_cl > 0) {
        Omega_cluster = quadcross(G_cl, G_cl) / m.N
        V_cluster     = Omega_cluster / m.N
    }
    else {
        Omega_cluster = J(0, 0, 0)
        V_cluster     = J(0, 0, 0)
    }
}

/* Explicit full-output diagnostic/testing convenience. This intentionally
   materializes N-by-all-link outputs. Production consumers use chunk() with
   bounded requested families and never call this assembler. */
struct cdid_selected_link_out scalar cdid_selected_link_if(
    struct cdid_selected_model scalar m,
    | real scalar chunk_rows,
      real scalar batch_links,
      real colvector cluster_id)
{
    struct cdid_selected_link_out scalar out
    struct cdid_selected_chunk_out scalar ch
    real colvector u_chunk, l_batch, cl_vals, ix_c, c_uniq
    real matrix phi_u, G_cl, info_cl
    real scalar r_step, l_step, lo, hi, l_lo, l_hi, p_kept, p_trt, b
    real scalar is_sep, tot_p, tot_des, p_b, at, k, l, n_cl

    if (!m.prepared | m.ctl.points_only) {
        _cdid_sel_link_fail("replay",
            "Model is not prepared for IF replay or was fitted in points_only mode.", 459)
    }

    r_step = (args() >= 2 & !missing(chunk_rows) & chunk_rows > 0 ? chunk_rows :
             (m.ctl.max_chunk_rows > 0 ? m.ctl.max_chunk_rows : m.N))
    l_step = (args() >= 3 & !missing(batch_links) & batch_links > 0 ? batch_links : m.n_links)
    if (r_step < 1 | r_step != floor(r_step) | l_step < 1 | l_step != floor(l_step)) {
        _cdid_sel_link_fail("input", "chunk_rows and batch_links must be positive integers.", 459)
    }

    is_sep  = (m.outcome_pooling == "separate")
    p_kept  = m.or_shared.p_kept
    tot_p   = (is_sep ? m.total_p_sep : p_kept)
    tot_des = (is_sep ? m.total_des_sep : m.n_links * p_kept)
    p_trt   = (m.n_trt_models > 0 ? m.trt_models[1].p_trt : 0)

    out.estimates             = m.estimates
    out.mu_ipw                = m.mu_ipw_vec
    out.mu_or                 = m.mu_or_vec
    out.mu_dr                 = m.mu_dr_vec
    if (m.ctl.est_method == "ipw") {
        out.or_beta = J(0, 1, 0)
    }
    else if (is_sep) {
        out.or_beta = J(tot_p, 1, 0)
        for (b = 1; b <= m.n_blocks; b = b + 1) {
            if (m.or_blocks[b].p_kept > 0) {
                out.or_beta[m.or_blocks[b].beta_col_start..m.or_blocks[b].beta_col_end] = m.or_blocks[b].beta
            }
        }
    }
    else {
        out.or_beta = m.or_shared.beta
    }
    out.or_alpha              = (m.ctl.est_method=="ipw" ? J(0,1,0) : J(m.n_blocks,1,0))
    if(m.ctl.est_method!="ipw") for (b = 1; b <= m.n_blocks; b = b + 1) out.or_alpha[b] = m.or_blocks[b].alpha
    out.or_objective          = m.or_shared.objective
    out.n_selection_fit_calls = m.n_selection_fit_calls
    out.n_selection_blocks = m.n_selection_blocks
    out.n_treatment_fit_calls = m.n_treatment_fit_calls
    out.n_or_block_fits       = m.n_or_block_fits
    out.n_imputation_calls    = m.n_imputation_calls

    out.phi              = J(m.N, m.n_links, 0)
    out.phi_core         = J(m.N, m.n_links, 0)
    out.psi_dir          = (m.ctl.est_method == "dr" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.psi_or           = (m.ctl.est_method == "dr" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.psi_ps           = (m.ctl.est_method == "dr" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.psi_sel          = (m.ctl.est_method == "dr" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.psi_ipw          = (m.ctl.est_method != "or" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.psi_or_estimator = (m.ctl.est_method != "ipw" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.xi_beta          = (m.ctl.est_method != "ipw" ? J(m.N, tot_p, 0) : J(0, 0, 0))
    out.beta_score       = (m.ctl.est_method != "ipw" ? J(m.N, tot_p, 0) : J(0, 0, 0))
    out.xi_alpha         = (m.ctl.est_method != "ipw" ? J(m.N, m.n_blocks, 0) : J(0, 0, 0))
    out.alpha_score      = (m.ctl.est_method != "ipw" ? J(m.N, m.n_blocks, 0) : J(0, 0, 0))
    out.den_if           = (m.ctl.est_method != "ipw" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.num_if           = (m.ctl.est_method != "ipw" ? J(m.N, m.n_links, 0) : J(0, 0, 0))
    out.design_if        = (m.ctl.est_method != "ipw" ? J(m.N, tot_des, 0) : J(0, 0, 0))
    out.objective_if     = (m.ctl.est_method != "ipw" ? J(m.N, 1, 0) : J(0, 1, 0))
    out.coord_ipw        = (m.ctl.est_method != "or" ? J(m.N, 4 * m.n_links, 0) : J(0, 0, 0))
    out.coord_dr         = (m.ctl.est_method == "dr" ? J(m.N, 4 * m.n_links, 0) : J(0, 0, 0))
    out.xi_gamma         = (m.ctl.est_method != "or" ? J(m.N, p_trt * m.n_links, 0) : J(0, 0, 0))
    out.norm_A_if        = (m.ctl.est_method != "or" ? J(m.N, 2 * m.n_links, 0) : J(0, 0, 0))
    out.tau_if           = J(m.N, m.n_cohorts, 0)
    out.eta_if           = J(m.N, m.n_cohorts, 0)
    out.lambda_if        = J(m.N, 1, 0)

    for (lo = 1; lo <= m.N; lo = lo + r_step) {
        hi = min((m.N, lo + r_step - 1))
        u_chunk = (lo::hi)
        for (l_lo = 1; l_lo <= m.n_links; l_lo = l_lo + l_step) {
            l_hi = min((m.n_links, l_lo + l_step - 1))
            l_batch = (l_lo::l_hi)
            ch = cdid_selected_chunk(m, u_chunk, l_batch, 0, 1)

            out.phi[u_chunk, l_batch]      = ch.phi
            out.phi_core[u_chunk, l_batch] = m.lambda :* ch.phi

            if (l_lo == 1) {
                out.tau_if[u_chunk, ]   = ch.tau_if
                out.eta_if[u_chunk, ]   = ch.eta_if
                out.lambda_if[u_chunk]  = ch.lambda_if
                if (m.ctl.est_method != "ipw") {
                    if (!is_sep & p_kept > 0) {
                        out.xi_beta[u_chunk, ]    = ch.xi_beta
                        out.beta_score[u_chunk, ] = ch.beta_score
                    }
                    out.objective_if[u_chunk]  = ch.objective_if
                }
            }
            if (m.ctl.est_method != "ipw") {
                out.xi_alpha[u_chunk,ch.alpha_block_ids']=ch.xi_alpha
                out.alpha_score[u_chunk,ch.alpha_block_ids']=ch.alpha_score
                out.psi_or_estimator[u_chunk, l_batch] = ch.psi_or_estimator
                out.den_if[u_chunk, l_batch]           = ch.den_if
                out.num_if[u_chunk, l_batch]           = ch.num_if
                if (is_sep) {
                    at = 0
                    for (k = 1; k <= rows(ch.alpha_block_ids); k = k + 1) {
                        b = ch.alpha_block_ids[k]
                        p_b = m.or_blocks[b].p_kept
                        if (p_b > 0) {
                            out.xi_beta[u_chunk, m.or_blocks[b].beta_col_start..m.or_blocks[b].beta_col_end] =
                                ch.xi_beta[, (at + 1)..(at + p_b)]
                            out.beta_score[u_chunk, m.or_blocks[b].beta_col_start..m.or_blocks[b].beta_col_end] =
                                ch.beta_score[, (at + 1)..(at + p_b)]
                            at = at + p_b
                        }
                    }
                    at = 0
                    for (k = 1; k <= rows(l_batch); k = k + 1) {
                        l = l_batch[k]
                        p_b = m.or_blocks[m.links[l].block_id].p_kept
                        if (p_b > 0) {
                            out.design_if[u_chunk, m.links[l].des_col_start..m.links[l].des_col_end] =
                                ch.design_if[, (at + 1)..(at + p_b)]
                            at = at + p_b
                        }
                    }
                }
                else if (p_kept > 0) {
                    out.design_if[u_chunk, ((l_lo - 1) * p_kept + 1)..(l_hi * p_kept)] = ch.design_if
                }
            }
            if (m.ctl.est_method != "or") {
                out.psi_ipw[u_chunk, l_batch] = ch.psi_ipw
                out.coord_ipw[u_chunk, ((l_lo - 1) * 4 + 1)..(l_hi * 4)] = ch.coord_ipw
                out.xi_gamma[u_chunk, ((l_lo - 1) * p_trt + 1)..(l_hi * p_trt)] = ch.xi_gamma
                out.norm_A_if[u_chunk, ((l_lo - 1) * 2 + 1)..(l_hi * 2)] = ch.norm_A_if
            }
            if (m.ctl.est_method == "dr") {
                out.psi_dir[u_chunk, l_batch] = ch.psi_dir
                out.psi_or[u_chunk, l_batch]  = ch.psi_or
                out.psi_ps[u_chunk, l_batch]  = ch.psi_ps
                out.psi_sel[u_chunk, l_batch] = ch.psi_sel
                out.coord_dr[u_chunk, ((l_lo - 1) * 4 + 1)..(l_hi * 4)] = ch.coord_dr
            }
        }
    }

    out.Omega_unit = J(m.n_links, m.n_links, 0)
    n_cl = 0
    if (args() >= 4 & rows(cluster_id) == m.N) {
        ix_c = order(cluster_id, 1)
        cl_vals = cluster_id[ix_c]
        u_chunk = J(m.N, 1, 0)
        u_chunk[ix_c] = runningsum(m.N > 1 ? (1 \ (cl_vals[2..m.N] :!= cl_vals[1..(m.N - 1)])) : 1)
        cl_vals = u_chunk
        n_cl = cl_vals[ix_c[m.N]]
        G_cl = J(n_cl, m.n_links, 0)
    }

    for (lo = 1; lo <= m.N; lo = lo + r_step) {
        hi = min((m.N, lo + r_step - 1))
        u_chunk = (lo::hi)
        phi_u = out.phi[u_chunk, ]
        out.Omega_unit = out.Omega_unit + quadcross(phi_u, phi_u)
        if (n_cl > 0) {
            ix_c = cl_vals[u_chunk]
            if (rows(u_chunk) > 1) {
                u_chunk = order((ix_c, (1::rows(u_chunk))), (1, 2))
                info_cl = _cdid_ps_slink(ix_c[u_chunk], 1)
                c_uniq  = ix_c[u_chunk[info_cl[, 1]]]
                G_cl[c_uniq, ] = G_cl[c_uniq, ] + panelsum(phi_u[u_chunk, ], info_cl)
            }
            else {
                G_cl[ix_c[1], ] = G_cl[ix_c[1], ] + phi_u
            }
        }
    }
    out.Omega_unit = out.Omega_unit / m.N
    out.V_unit     = out.Omega_unit / m.N
    if (n_cl > 0) {
        out.Omega_cluster = quadcross(G_cl, G_cl) / m.N
        out.V_cluster     = out.Omega_cluster / m.N
    }
    else {
        out.Omega_cluster = J(0, 0, 0)
        out.V_cluster     = J(0, 0, 0)
    }

    return(out)
}

end
