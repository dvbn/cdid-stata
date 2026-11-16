*! Augmented observation-correction backend; StataNow 19.5.
*! selection_est_method "dr" with treatment est_method "dr" ONLY: conditional_mtar on the
*! representative roster F, rho = 1, known F cohorts, both q models, all q/j poolings.
*! Pinned R cdid 1.1.0.0 (), R/chained_augmented.R:
*!   augmented_cell_rows, fit_augmented_backend (steps 1 and 4-6), augmented_if_needed_sets,
*!   augmented_evaluate_if_chunk_core (non-layout path); R/chained_augmented_public.R gates.
*! Load the unchanged cdid_population, cdid_nuisance, cdid_ipw_link,
*! cdid_imputation_model, cdid_gmm, cdid_imputation_dated, cdid_imputation_backend and
*! cdid_selection modules first. No public command, option, default or tolerance is added.
version 19.5

mata:
mata set matastrict on

/* Augmented fields of one Stage 1 cell (source cells[[cid]] additions). Yaug_z is raw
   slope order (1 x p_raw); MW_d and MD_d are kept order of the cell's Stage 1 group
   (1 x p_g). theta_ids are NATIVE global selection theta positions (alpha_all, beta_all)
   in the local order of cdid_sel_log_deriv_eval; Rq_one is aligned with theta_ids.
   n_floor counts observed rows whose raw q <= 1e-6 (strict derivative mask off). */
struct _cdid_aug_cell {
    real scalar q_pc, n_obs, Yaug_one, MW_one, MD_one
    real scalar q_min, q_max, n_floor, w_max
    real rowvector Yaug_z, MW_d, MD_d
    real colvector theta_ids, Rq_one
}

/* One deduplicated block of the Z-response common OR (pair, recipe, ordered controls) with
   its upstream forcing. jf_alpha_beta is n_ctrl x p_raw and jf_beta_beta is
   p_or x (n_ctrl p_raw), raw slope order of each control cell group (source
   j_force_alpha_beta, j_force_beta_beta). q_alpha merges -Rq_one of the block's control
   cells by native theta ID. */
struct _cdid_aug_block {
    string scalar name
    real scalar block_id, pair_id, recipe, L, n_ctrl, weight_sum, mass, alpha, mean_y
    real scalar or_p, or_rss, or_rcond, or_adj_ready, score_mean_a
    real colvector ctrl_ci, ctrl_cell, links
    real rowvector mean_d, or_keep, or_pivot
    real colvector a_force_alpha, jf_alpha_alpha, q_alpha_ids, q_alpha, or_beta, score_mean_b, q_beta_ids
    real matrix a_force_beta, jf_beta_alpha, jf_alpha_beta, jf_beta_beta
    real matrix or_R, or_alias_map, or_alias_const, q_beta
    struct cdid_imputation_adjoint scalar or_adj
}

/* One computational link (R links[[ell]] record of step 6). nu_T is kept order of the target
   cell group; nu_C_k is n_ctrl x p_raw raw order. sel_* carry native theta IDs. */
struct _cdid_aug_link {
    string scalar link_id
    real scalar link_index, pair_id, recipe, L, target_ci, target_cell, or_block, trt_id, n_ctrl
    real colvector ctrl_ci, ctrl_cell
    real scalar mu, mu_ipw, mu_reg, mu_dr, d_plus, num_T, T_gl, num_Bg, B_gl, d_minus, num_C, C_ol, num_Bo, B_ol
    real scalar a_coeff_T_g, a_coeff_Bg_g, lam_T
    real rowvector nu_T, d_bar_plus_or, d_bar_minus_or, d_bar_diff_or
    real colvector a_coeff_C_k, a_coeff_Bo_k, m_minus_k, lam_C_k
    real matrix nu_C_k
    real rowvector denminus_gamma, g_C_gamma, g_Bo_gamma, g_dr_gamma
    real colvector sel_T_ids, sel_T, sel_C_ids, sel_C
}

/* Prepared augmented backend. ipw (est_method ipw backend, which itself borrows the
   Stage 1 fit and treatment design), model and qrdr (raw-W observation reader) are
   BORROWED: the caller owns them and keeps them alive and unmodified for the whole
   lifetime; every replay re-checks their fingerprints. The record OWNS the observation
   partition/fit/up-map, the compact augmented cell moments, the Z-response common OR
   factors, forcing and link loads. No N-length array is retained here. */
/* Compact recipe-local q state. Raw designs are unique-date slots, weights stay
   borrowed from the immutable dated reader, and only one original-F reader is live. */
struct _cdid_aug_qcontext {
    struct cdid_sel_fit_result scalar fit
    struct cdid_sel_up_map scalar up
    struct cdid_sel_design_spec scalar spec
    real scalar theta_offset, q_slot, w_slot, n_U
    real colvector mass_F, mass_U, repeat_mass, a_U, U
}
struct cdid_aug_backend {
    string scalar backend_id, method, selection_est_method, selection, selection_model
    string scalar selection_pooling, pooling, outcome_pooling, law, normalization
    pointer(struct cdid_imp_backend scalar) scalar ipw
    pointer(struct cdid_imputation_model scalar) scalar model
    pointer(struct cdid_sel_reader scalar) scalar qrdr
    struct cdid_imputation_budget scalar budget
    pointer(struct cdid_imp_dated_reader scalar) scalar qdated
    struct cdid_imp_matrix_slot colvector qdata
    struct _cdid_aug_qcontext colvector qcontexts
    struct cdid_sel_reader scalar qlease
    real scalar qactive, qslot_active, has_or, has_ps
    real scalar N, n_A, lambda, G, n_cells, n_groups, n_recipes, p_raw, intercept, q, d_sel
    real scalar rank_tol, glm_eps, glm_maxit, link_batch, chunk_rows, is_dated
    real scalar n_links, n_blocks, n_trt, n_q_params, n_q_pc
    real scalar n_or_fits, n_q_fits, n_ps_fits, n_imputation_fits
    real scalar ipw_sum, q_sum, retained_bytes, scratch_bytes, work_ops
    struct cdid_sel_fit_result scalar qfit
    struct cdid_sel_up_map scalar qum
    real colvector q_of_cell, sel_cohort_pos
    struct _cdid_aug_cell rowvector cells
    struct _cdid_aug_block rowvector blocks
    struct _cdid_aug_link rowvector links
    real scalar or_p, or_rss, or_rcond, or_rows_checked, or_threshold, or_adj_ready
    real rowvector or_keep, or_retained, or_proposed, or_dropped, or_candidates, or_pivot
    real rowvector or_ctl_scaled, or_trt_scaled
    real matrix or_R, or_alias_map, or_alias_const, q_beta
    real colvector or_beta, q_beta_ids
    struct cdid_imputation_adjoint scalar or_adj
}

/* Replay result: B rows (explicit original-F positions) by L requested links. comp is
   B x 10L and ref4 is B x 4L, link-major, in the source SCHEMA order: direct, rate,
   imputation, selection, treatment_direct, treatment_rate, OR_direct, OR_rate,
   OR_imputation, OR_selection; and j_score_pullback, q_score_pullback, all_rate_paths,
   remaining_direct_and_downstream_scores. Both are empty unless detail = 1. */
struct cdid_aug_backend_if {
    string scalar method, scale
    real colvector f_rows, links, id_num
    string colvector id_str
    real scalar normalization_count, N, n_A, lambda, detail
    real matrix phi, comp, ref4
}

/* Cached representative propensity group for batched augmented IF contraction. */
struct _cdid_aug_tp_grp {
    real scalar rtid
    real rowvector link_pos
    real matrix G_gamma
}

/* Pre-built augmented IF batch plan (detail = 0) shared across all row chunks. */
struct _cdid_aug_plan {
    real scalar ready, nl, so, nA, tot_po, sep
    real colvector lk, abids, bbids, bids, boff, zids, zpos, rep_tids, link_tp, or_bids, or_bpos
    real matrix or_coh_mask, A_load, Alpha_load, Beta_load, Raa_load, Rba_load
    struct _cdid_aug_tp_grp rowvector tp_grps
    struct cdid_imputation_adjoint scalar adj
    struct cdid_sel_moments scalar sm
}

/* GMM weight-provider context; the backend is borrowed with the same lifetime rule. */
struct cdid_aug_backend_wctx {
    pointer(struct cdid_aug_backend scalar) scalar backend
    real scalar n_moments, row_chunk
    real colvector need_links, col_pos, col_sign
    string colvector consumer_ids
    real scalar fast_ready
    struct _cdid_aug_plan scalar plan
}

/* Vectorized raw probability, capped odds, and strict derivative mask matching _cdid_ipw_logit_values. */
real matrix _cdid_aug_logit_vals(real colvector eta)
{
    real colvector p, odds, active, hi, ec
    real scalar cap, cap_odds
    if (missing(eta)) _cdid_nuisance_fail("Nonfinite logit predictor.", 498)
    if (rows(eta) == 0) return(J(0, 3, 0))
    hi = (eta :> 40)
    ec = exp(rowmin((eta, J(rows(eta), 1, 40))))
    p = (hi :* (1 :/ (1 :+ exp(-eta :* hi)))) :+ ((1 :- hi) :* (ec :/ (1 :+ ec)))
    active = cdid_ipw_odds_active(p)
    cap = 1 - 1e-6
    cap_odds = cap / (1 - cap)
    odds = active :* ec :+ (1 :- active) :* cap_odds
    if (missing(odds) | any(odds :<= 0))
        _cdid_ipw_fail("Propensity odds must be finite and strictly positive.", 504)
    return((p, odds, active))
}

struct _cdid_aug_bkacc {
    real scalar s_r
    real rowvector s_rd, s_vd
    real matrix s_vdd, Rqd
}

struct _cdid_aug_lacc {
    real scalar s_wm, s_wmZ, s_wmm, po, alpha_b
    real colvector s_wm_k, s_wodZ_k, s_wod_k, s_wodm_k, lamC, beta_b
    real rowvector s_wdZX, s_wdX, s_wdmX, s_wmdor, dbp, keep_b
    real matrix nuC
}

struct _cdid_aug_rcc {
    real scalar ready, has_x
    real colvector wr, v
    real matrix Zr, X
}

struct _cdid_aug_cev {
    real scalar ready
    real colvector zr, obs, piv, ev
    real matrix Tm
}

void _cdid_aug_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

/* Positions of native theta IDs inside a sorted unique ID support (all must be present). */
real colvector _cdid_aug_find(real colvector keys, real colvector ids)
{
    real colvector pos
    real scalar i
    pos = J(rows(ids), 1, 0)
    for (i=1; i<=rows(ids); i++) {
        pos[i] = _cdid_imp_find_exact(keys, ids[i])
        if (pos[i] < 1)
            _cdid_aug_fail("Native selection theta ID outside its merged loading support.", 498)
    }
    return(pos)
}

/* Sorted unique native theta IDs of a set of cells. */
real colvector _cdid_aug_theta_union(struct cdid_aug_backend scalar a, real colvector cids)
{
    real colvector ids
    real scalar i
    ids = J(0, 1, 0)
    for (i=1; i<=rows(cids); i++) ids = ids \ a.cells[cids[i]].theta_ids
    if (rows(ids) == 0) return(J(0, 1, 0))
    return(uniqrows(ids))
}

/* Order cells by calendar pair (s_index, t_index), recipe, cohort_index, and cid so pair
   chunks hit the LRU cache sequentially while keeping cohort order within each pair. */
real colvector _cdid_aug_order_cells(struct cdid_imputation_model scalar m, real colvector cids)
{
    real matrix keys
    real scalar n, i, cid
    n = rows(cids)
    if (n <= 1) return(cids)
    keys = J(n, 5, 0)
    for (i=1; i<=n; i++) {
        cid = cids[i]
        keys[i, ] = (m.cells[cid].s_index, m.cells[cid].t_index, m.cells[cid].recipe,
            m.cells[cid].cohort_index, cid)
    }
    return(cids[order(keys, (1, 2, 3, 4, 5))])
}

/* Private request gate in the source refusal order (chained_request_augmented_est_method,
   then chained_request_augmented_declaration): custom callback, treatment method, assumption,
   universe, then the selection declaration normalizer (model, pooling, legacy
   corollary1 alias, rho and justification). Refusals are explicit; nothing is substituted.
   The public binder calls this before any data are read. */
struct cdid_sel_decl scalar cdid_aug_route_gate(string scalar est_method, real scalar is_custom,
    string scalar selection, string scalar selection_universe, string scalar selection_model,
    string scalar selection_pooling, real scalar has_rho, string scalar justification)
{
    struct cdid_sel_decl scalar d
    if (is_custom != 0)
        _cdid_aug_fail("selection_est_method 'dr' (augmented observation correction) with a custom treatment callback is unsupported: callbacks have no specified augmented forcing or complete upstream score replay. Nothing is substituted.", 459)
    if (est_method != "dr" & est_method != "ipw" & est_method != "reg")
        _cdid_aug_fail("selection_est_method 'dr' (augmented observation correction) supports only est_method 'dr', 'ipw', or 'reg'; est_method '" + est_method + "' is refused before the data are evaluated and no other correction or treatment method is substituted.", 459)
    if (selection != "conditional_mtar" & selection != "corollary1")
        _cdid_aug_fail("selmethod(dr) requires selection(conditional_mtar).", 459)
    if (selection_universe != "F")
        _cdid_aug_fail("selection_est_method 'dr' is defined only on the representative roster selection_universe 'F'; B/A fitting rosters are unavailable for this method.", 459)
    d = cdid_sel_norm_decl(selection, 1, "F", has_rho, (strtrim(justification) != ""),
        selection_model, selection_pooling, justification)
    if (!d.enabled)
        _cdid_aug_fail("The augmented observation correction requires an active selection declaration.", 459)
    return(d)
}

/* Borrowed-owner fingerprints (exact, not a tolerance): the roster-propensity roots and
   repeat-rate forcing of the borrowed ipw backend, and the compact descriptors of the
   borrowed observation reader. Stage 1 itself is checked by _cdid_impb_check. */
real scalar _cdid_aug_ipw_sum(struct cdid_imp_backend scalar ib)
{
    real scalar s, tid
    s = ib.borrow_sum + 3 * ib.n_links + 5 * ib.n_trt + 7 * ib.q
    for (tid=1; tid<=ib.n_trt; tid++) {
        if (rows(ib.trts[tid].gamma) > 0)
            s = s + (tid + 0.5) * quadsum(ib.trts[tid].gamma :* (1::rows(ib.trts[tid].gamma)))
        s = s + tid * quadsum(vec(ib.trts[tid].a_forcing))
    }
    return(s)
}

real scalar _cdid_aug_q_sum(struct cdid_sel_reader scalar qr)
{
    real scalar G
    G = rows(qr.cohorts)
    if (G < 1) return(.)
    return(qr.N + 3 * qr.n_U + 5 * qr.n_A + 7 * rows(qr.Z) + 11 * cols(qr.Z) + 13 * qr.spec.d_s +
        quadsum(qr.mass_F :* (1::G)) + 17 * quadsum(qr.a_U :* (1::G)))
}

void _cdid_aug_check(struct cdid_aug_backend scalar a)
{
    pointer(struct cdid_imputation_model scalar) scalar pm
    if (a.ipw == NULL)
        _cdid_aug_fail("Borrowed ipw-factor backend pointer is NULL; its lifetime ended.", 498)
    if (!a.is_dated & a.qrdr == NULL)
        _cdid_aug_fail("Borrowed observation-probability reader pointer is NULL; its lifetime ended.", 498)
    _cdid_impb_need_model(*(a.ipw))
    pm = (*(a.ipw)).model
    if (pm != a.model)
        _cdid_aug_fail("Borrowed ipw-factor backend now points to a different first-step fit; prepare again.", 498)
    _cdid_impb_check(*(a.ipw), *pm)
    if (_cdid_aug_ipw_sum(*(a.ipw)) != a.ipw_sum)
        _cdid_aug_fail("Borrowed roster propensity factors changed after augmented preparation; prepare again.", 498)
    if (a.is_dated) {
        if (a.qdated == NULL) _cdid_aug_fail("Borrowed dated q reader object no longer exists.", 498)
        _cdid_impd_check_reader(*(a.qdated))
    }
    else if (_cdid_aug_q_sum(*(a.qrdr)) != a.q_sum)
        _cdid_aug_fail("Borrowed observation-probability reader changed after augmented preparation; prepare again.", 498)
}

/* Activate one recipe without retaining N_F-by-recipes readers/designs. The
   unique-date Q buffers are leased by swap, following the selected reader. */
void _cdid_aug_q_activate(struct cdid_aug_backend scalar a, real scalar r)
{
    real scalar slot
    pointer(struct cdid_imp_dated_reader scalar) scalar dr
    if (a.qactive == r) return
    if (r < 1 | r > a.n_recipes) _cdid_aug_fail("Invalid q design block.", 498)
    dr = a.qdated
    slot = a.qcontexts[r].q_slot
    if (a.qslot_active != slot) {
        if (a.qslot_active) swap(a.qlease.Z, a.qdata[a.qslot_active].value)
        swap(a.qlease.Z, a.qdata[slot].value)
        a.qslot_active = slot
    }
    a.qlease.w = (*(*dr).reg).w_data[a.qcontexts[r].w_slot].value
    a.qlease.spec = a.qcontexts[r].spec
    a.qlease.n_U = a.qcontexts[r].n_U
    a.qlease.mass_F = a.qcontexts[r].mass_F
    a.qlease.mass_U = a.qcontexts[r].mass_U
    a.qlease.repeat_mass = a.qcontexts[r].repeat_mass
    a.qlease.a_U = a.qcontexts[r].a_U
    if (rows(a.qcontexts[r].U)) a.qlease.U = a.qcontexts[r].U
    a.qactive = r
}

/* Contract only each recipe's theta rows, then sum on the ORIGINAL roster. All
   G_ak loads remain zero: repeat-rate paths are composed separately exactly once. */
real matrix _cdid_aug_q_loading(struct cdid_aug_backend scalar a,
    struct cdid_sel_moments scalar global_load, real colvector fr)
{
    struct cdid_sel_moments scalar local_load
    real scalar r, off, p
    real matrix out
    out = J(rows(fr), global_load.n_moments, 0)
    for (r=1; r<=a.n_recipes; r++) {
        p = a.qcontexts[r].fit.total_params
        if (p == 0) continue
        off = a.qcontexts[r].theta_offset
        local_load = cdid_sel_moments_init(a.qcontexts[r].fit, global_load.n_moments)
        local_load.G_theta = global_load.G_theta[(off+1)..(off+p), ]
        if (!any(local_load.G_theta :!= 0)) continue
        _cdid_aug_q_activate(a, r)
        out = out + cdid_sel_eval_loading(a.qcontexts[r].fit, a.qlease, local_load, fr)
    }
    return(out)
}

/* Source augmented_cell_rows on one bounded chunk rch of rows of the cell's cohort: the
   Stage 1 prediction j on every row and, on observed pairs only (S = 1), the RAW pair
   probability q (endpoint product for independent_period_logit), ONE floor
   pi = 1/max(q, 1e-6), e = Delta - j and Z = j + pi e. Elsewhere Z = j and neither Delta nor
   the inverse probability is read. need_T adds T_logit = -d log max(q, 1e-6)/d theta in the
   cell's native local theta order (cdid_sel_log_deriv_eval: already signed, zero where the
   strict q > 1e-6 mask is off); its d_ak rate column is NOT used. */
void _cdid_aug_cell_rows_z(struct cdid_aug_backend scalar a, struct cdid_imputation_model scalar m,
    struct cdid_sel_reader scalar qr, real scalar cid, real colvector rch, real scalar need_T,
    real colvector jh, real colvector zr, real colvector obs, real colvector qv,
    real colvector piv, real colvector ev, real matrix Tm)
{
    struct cdid_sel_log_deriv scalar ld
    struct cdid_sel_pred scalar pr
    real matrix hit
    real colvector qf
    real scalar pid, pc, dq
    pid = m.cells[cid].pair_id
    pc = a.q_of_cell[cid]
    if (a.is_dated) _cdid_aug_q_activate(a, m.cells[cid].recipe)
    dq = rows(a.cells[cid].theta_ids)
    hit = cdid_imputation_read_pair_chunk(m.reader, m.pop, m.pairs[pid], rch, 1)
    obs = _cdid_impb_ix(hit[, 1] :== 1)
    zr = jh
    qv = J(0, 1, 0)
    piv = J(0, 1, 0)
    ev = J(0, 1, 0)
    Tm = J(0, dq, 0)
    if (rows(obs) == 0) return
    if (need_T) {
        if (a.is_dated) ld = cdid_sel_log_deriv_eval(a.qcontexts[m.cells[cid].recipe].fit,
            a.qlease, a.qcontexts[m.cells[cid].recipe].up, pc, rch[obs])
        else ld = cdid_sel_log_deriv_eval(a.qfit, qr, a.qum, pc, rch[obs])
        qv = ld.p_pair
        qf = ld.p_weight
        Tm = ld.d_theta_local
        if (rows(Tm) != rows(obs) | cols(Tm) != dq)
            _cdid_aug_fail("Observation log-derivative rows disagree with the cell theta coordinates.", 498)
    }
    else {
        if (a.is_dated) pr = cdid_sel_pred_pair(a.qcontexts[m.cells[cid].recipe].fit,
            a.qlease, pc, rch[obs])
        else pr = cdid_sel_pred_pair(a.qfit, qr, pc, rch[obs])
        qv = pr.p_pair
        qf = pr.p_weight
    }
    if (rows(qv) != rows(obs) | rows(qf) != rows(obs) | missing(qv) | missing(qf))
        _cdid_aug_fail("Nonfinite fitted observation probabilities on observed pairs.", 504)
    piv = 1 :/ qf
    ev = hit[obs, 2] - jh[obs]
    zr[obs] = jh[obs] + piv :* ev
}

void _cdid_aug_cell_rows(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr, real scalar cid,
    real colvector rch, real scalar need_T, real colvector jh, real colvector zr,
    real colvector obs, real colvector qv, real colvector piv, real colvector ev, real matrix Tm)
{
    if (a.is_dated) _cdid_aug_q_activate(a, m.cells[cid].recipe)
    jh = _cdid_impb_pred_cell(ib, m, cid, rch)
    _cdid_aug_cell_rows_z(a, m, qr, cid, rch, need_T, jh, zr, obs, qv, piv, ev, Tm)
}

/* Step 0: original-F domain, borrowed-owner identity and source support gates. The observation
   reader and Stage 1 must share the original population (N_F, A, calendar), the RAW
   weight W and every unit cohort; cohorts are matched by VALUE because Stage 1 orders
   never-treated first and the selection reader orders it last. */
void _cdid_aug_init(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr)
{
    real scalar ci, k, val
    real colvector selpos
    if (a.has_ps) {
        if (ib.method != "ipw" | ib.has_ps != 1 | ib.has_or != 0 | ib.stage1_only)
            _cdid_aug_fail("The augmented ipw/dr backend borrows the est_method ipw preparation (propensity factors, no OR); a pure-imputation OR fitted on j is not the augmented Z-response OR.", 498)
    }
    else {
        if (ib.method != "reg" | !ib.stage1_only | ib.has_ps != 0 | ib.has_or != 0)
            _cdid_aug_fail("The augmented reg backend requires the first-step helper preparation (stage1_only reg, no PS or pure-imputation OR).", 498)
    }
    if (!a.is_dated & (ib.is_dated | m.n_recipes != 1))
        _cdid_aug_fail("This static augmented primitive accepts the static single-block first-step backend only; the dated route needs the scoped used-cohort observation reader.", 459)
    if (ib.glm_eps <= 0 | ib.glm_eps >= . | ib.glm_maxit < 1 | ib.glm_maxit >= . |
        ib.glm_maxit != floor(ib.glm_maxit))
        _cdid_aug_fail("Observation logits require an explicit positive GLM epsilon and integer maxit.", 498)
    if (m.N < 1 | rows(m.pop.cohort_known) != m.N | any(m.pop.cohort_known :!= 1))
        _cdid_aug_fail("The augmented observation correction requires known cohorts on every F row, including frame-only rows.", 459)
    if (!qr.has_pop | qr.kind != "population")
        _cdid_aug_fail("The observation-probability reader must be the population-backed cdid_sel_reader_pop reader.", 498)
    /* Under a restricted reference eligibility context the universe of a recipe is
       its eligible part of F (the reader passed here is that of the first recipe). */
    if (qr.selection_universe != "F" | rows(qr.U) != m.N)
        _cdid_aug_fail("The augmented observation correction requires the full representative roster universe F with U = 1 on every F row.", 459)
    if (rows(m.elig_keep)) {
        if (any(qr.U :!= m.elig_keep[, 1]))
            _cdid_aug_fail("The augmented observation correction requires the full representative roster universe F with U = 1 on every F row.", 459)
    }
    else {
        if (qr.n_U != m.N | any(qr.U :!= 1))
            _cdid_aug_fail("The augmented observation correction requires the full representative roster universe F with U = 1 on every F row.", 459)
    }
    if (rows(qr.rho) != m.N)
        _cdid_aug_fail("The augmented observation correction requires rho = 1 on every F row.", 459)
    if (any(qr.rho :!= 1))
        _cdid_aug_fail("The augmented observation correction requires rho = 1 on every F row.", 459)
    if (qr.selection_model == "independent_period_logit" & strtrim(qr.justification) != "")
        _cdid_aug_fail("The augmented observation correction on F uses no transport justification.", 459)
    if (qr.selection_model != "pair_logit" & qr.selection_model != "independent_period_logit")
        _cdid_aug_fail("Observation model must be pair_logit or independent_period_logit.", 498)
    if (qr.selection_pooling != "separate" & qr.selection_pooling != "within_cohort" &
        qr.selection_pooling != "across_cohorts")
        _cdid_aug_fail("Observation pooling must be separate, within_cohort or across_cohorts.", 498)
    if (qr.N != m.N | qr.n_A != m.n_A | rows(qr.A) != m.N | rows(qr.w) != m.N |
        rows(qr.all_times) != rows(m.pop.calendar) | rows(qr.cohort_pos) != m.N)
        _cdid_aug_fail("Observation reader and first-step fit do not share the original population.", 498)
    if (any(qr.A :!= m.pop.A) | any(qr.all_times :!= m.pop.calendar))
        _cdid_aug_fail("Observation reader and first-step fit disagree on repeat observers or the calendar.", 498)
    if (!a.is_dated) if (any(qr.w :!= m.recipes[1].weight))
        _cdid_aug_fail("Observation probabilities must be fitted with the same raw W as the first-step design block (never W a or W a/q).", 498)
    selpos = J(m.n_cohorts, 1, 0)
    for (ci=1; ci<=m.n_cohorts; ci++) {
        val = (m.cohorts[ci] == 0 ? . : m.cohorts[ci])
        for (k=1; k<=rows(qr.cohorts); k++) {
            if (qr.cohorts[k] == val) {
                selpos[ci] = k
                break
            }
        }
        if (selpos[ci] == 0)
            _cdid_aug_fail("First-step cohort " + m.cohort_labels[ci] + " is absent from the observation reader cohorts.", 498)
    }
    if (any(selpos[m.unit_cohort_index] :!= qr.cohort_pos))
        _cdid_aug_fail("Observation reader and first-step fit disagree on unit cohorts.", 498)
    a.sel_cohort_pos = selpos
    a.N = m.N
    a.n_A = m.n_A
    a.lambda = m.lambda
    a.G = m.n_cohorts
    a.n_cells = m.n_cells
    a.n_groups = m.n_groups
    a.n_recipes = m.n_recipes
    a.p_raw = m.p_raw
    a.intercept = m.intercept
    a.q = ib.q
    a.d_sel = cols(qr.Z)
    a.rank_tol = m.rank_tol
    a.glm_eps = ib.glm_eps
    a.glm_maxit = ib.glm_maxit
    a.pooling = m.pooling
    a.selection_model = qr.selection_model
    a.selection_pooling = qr.selection_pooling
    a.chunk_rows = min((a.budget.max_chunk_rows, ib.chunk_rows))
    a.n_trt = ib.n_trt
    a.n_imputation_fits = 0
    a.n_ps_fits = 0
    a.n_or_fits = 0
    a.n_q_fits = 0
    a.or_p = 0
    a.or_adj_ready = 0
    a.law = ("AUG-" + (a.method == "ipw" ? "IPW" : (a.method == "reg" ? "REG" : "DR")) + ":" +
        a.pooling + (a.has_or ? ":or-" + a.outcome_pooling : "") + ":" + a.selection_model + ":" +
        a.selection_pooling + ":" + (a.intercept ? "local-intercepts" : "no-intercept") +
        ":actual-pair:fixed-w:full-F")
    a.normalization = "frame IF on all original F rows; core = (n_A/N_F) frame after recipe summation; weight Gram centers all F core rows and divides by n_A"
    if (a.n_A < 1 | a.n_A > a.N | a.lambda != a.n_A / a.N)
        _cdid_aug_fail("Repeat-observer count n_A must lie in 1..N_F with lambda = n_A/N_F.", 459)
}

/* Link intake from the borrowed ipw registry (same Stage 1 pair/target/control cells
   and treatment keys), then source OR blocks keyed (pair_id [carries recipe], ordered controls)
   exactly as the _cdid_impb_registry. */
void _cdid_aug_registry(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m)
{
    real scalar L, ell, i, nb, bid, raw, pid, mlen, found
    real colvector kord, link_raw, ordb, rank_b
    real rowvector krow, kr
    real matrix bkeys
    L = ib.n_links
    if (L < 1) _cdid_aug_fail("The borrowed ipw backend has no links.", 498)
    a.n_links = L
    a.links = _cdid_aug_link(1, L)
    mlen = 0
    for (ell=1; ell<=L; ell++) {
        a.links[ell].link_id = ib.links[ell].link_id
        a.links[ell].link_index = ell
        a.links[ell].pair_id = ib.links[ell].pair_id
        a.links[ell].recipe = ib.links[ell].recipe
        a.links[ell].L = ib.links[ell].L
        a.links[ell].target_ci = ib.links[ell].target_ci
        a.links[ell].target_cell = ib.links[ell].target_cell
        a.links[ell].trt_id = ib.links[ell].trt_id
        a.links[ell].n_ctrl = ib.links[ell].n_ctrl
        a.links[ell].ctrl_ci = ib.links[ell].ctrl_ci
        a.links[ell].ctrl_cell = ib.links[ell].ctrl_cell
        a.links[ell].or_block = 0
        a.links[ell].mu = .
        a.links[ell].mu_ipw = .
        a.links[ell].mu_reg = .
        a.links[ell].mu_dr = .
        if (a.links[ell].n_ctrl < 1 | (a.has_ps & (a.links[ell].trt_id < 1 | a.links[ell].trt_id > ib.n_trt)))
            _cdid_aug_fail("Borrowed first-step link registry lacks a required treatment key or control set.", 498)
        if (a.links[ell].n_ctrl > mlen) mlen = a.links[ell].n_ctrl
    }
    if (a.has_or) {
        bkeys = J(0, 1 + mlen, 0)
        link_raw = J(L, 1, 0)
        for (ell=1; ell<=L; ell++) {
            krow = (a.links[ell].pair_id, a.links[ell].ctrl_ci', J(1, mlen - a.links[ell].n_ctrl, -1))
            found = 0
            for (i=1; i<=rows(bkeys); i++) {
                if (all(bkeys[i, ] :== krow)) {
                    found = i
                    break
                }
            }
            if (found == 0) {
                bkeys = bkeys \ krow
                found = rows(bkeys)
            }
            link_raw[ell] = found
        }
        nb = rows(bkeys)
        ordb = order(bkeys, _cdid_impb_seq(1, cols(bkeys)))
        rank_b = invorder(ordb)
        a.n_blocks = nb
        a.blocks = _cdid_aug_block(1, nb)
        for (bid=1; bid<=nb; bid++) {
            raw = ordb[bid]
            pid = bkeys[raw, 1]
            kr = bkeys[raw, 2..cols(bkeys)]
            kord = select(kr', kr' :> 0)
            a.blocks[bid].block_id = bid
            a.blocks[bid].pair_id = pid
            a.blocks[bid].recipe = m.pairs[pid].recipe
            a.blocks[bid].L = m.pairs[pid].L
            a.blocks[bid].n_ctrl = rows(kord)
            a.blocks[bid].ctrl_ci = kord
            a.blocks[bid].ctrl_cell = m.cell_lookup[kord, pid]
            a.blocks[bid].links = _cdid_impb_ix(rank_b[link_raw] :== bid)
            a.blocks[bid].name = (m.pairs[pid].pair_label + ":K={" +
                invtokens(m.cohort_labels[kord]', ",") + "}")
        }
        for (ell=1; ell<=L; ell++) a.links[ell].or_block = rank_b[link_raw[ell]]
    }
    else {
        a.n_blocks = 0
        a.blocks = _cdid_aug_block(1, 0)
    }
    a.cells = _cdid_aug_cell(1, a.n_cells)
}

/* catalog identity: the observation pair-cohort catalog IS C_imp. Every Stage 1 cell
   (cohort mapped by value, calendar pair) is one requested pair-cohort; the
   partition must return exactly one pair-cohort per cell, one-to-one. */
struct cdid_sel_cell_map scalar _cdid_aug_q_map(struct cdid_aug_backend scalar a,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr)
{
    struct cdid_sel_cell_map scalar qmap
    real matrix pspec
    real colvector seen
    real scalar cid, pc, hit
    pspec = J(a.n_cells, 3, 0)
    for (cid=1; cid<=a.n_cells; cid++) {
        if (m.cells[cid].recipe != 1)
            _cdid_aug_fail("Static augmented cells must belong to the single static design block.", 498)
        pspec[cid, ] = (a.sel_cohort_pos[m.cells[cid].cohort_index], m.cells[cid].s_index,
            m.cells[cid].t_index)
    }
    qmap = cdid_sel_partition(qr, pspec)
    if (qmap.n_pair_cohorts != a.n_cells)
        _cdid_aug_fail("The actual imputation cells and the observation-probability pair-cohort catalog are not one-to-one; no cell is reconstructed without its own fitted probability.", 498)
    a.q_of_cell = J(a.n_cells, 1, 0)
    seen = J(qmap.n_pair_cohorts, 1, 0)
    for (cid=1; cid<=a.n_cells; cid++) {
        hit = 0
        for (pc=1; pc<=qmap.n_pair_cohorts; pc++) {
            if (qmap.pair_cohorts[pc, 1] == pspec[cid, 1] & qmap.pair_cohorts[pc, 2] == pspec[cid, 2] &
                qmap.pair_cohorts[pc, 3] == pspec[cid, 3]) {
                if (hit > 0)
                    _cdid_aug_fail("Ambiguous observation pair-cohort for an imputation cell.", 498)
                hit = pc
            }
        }
        if (hit == 0)
            _cdid_aug_fail("The actual imputation cells and the observation-probability pair-cohort catalog are not one-to-one; no cell is reconstructed without its own fitted probability.", 498)
        if (seen[hit])
            _cdid_aug_fail("Two imputation cells resolve to one observation pair-cohort.", 498)
        seen[hit] = 1
        a.q_of_cell[cid] = hit
    }
    a.n_q_pc = qmap.n_pair_cohorts
    return(qmap)
}

/* Preflight of concurrently live retained bytes (the borrowed ipw backend already charges
   Stage 1 and the treatment design; the borrowed observation reader holds its design,
   vectors and population copy), the largest scan/fit/replay scratch and the opt-in
   work-ops estimate against the caller budget. */
void _cdid_aug_preflight(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    struct cdid_sel_cell_map scalar qmap)
{
    struct cdid_sel_budget scalar sb
    real scalar p, q, G, d, dq, ret, work, ell, bid, k, cmax, nlb, nbb, ntb, M_adj, M_sel
    real scalar sc_scan, sc_or, sc_rep, scr, rdr_bytes, cid, nq
    real colvector ncoh
    p = a.p_raw
    q = a.q
    G = a.G
    d = cols(qr.Z)
    dq = 2 + 2 * qr.spec.d_s
    cmax = a.chunk_rows
    sb = cdid_sel_budget_check(qr, qmap, .)
    /* Dated theta loads span ALL recipe contexts. The current reader/map is
       only one recipe; both retained p*nq and replay nq*M_sel use global nq. */
    nq = (a.is_dated ? a.n_q_params : sb.total_params)
    ncoh = J(G, 1, 0)
    for (k=1; k<=G; k++) ncoh[k] = rows(ib.coh_u[k].rows)
    ret = (128 + a.n_cells * (24 + 3 * p + 2 * dq) + sb.compact_state_doubles +
        9 * qmap.n_cells + 8 * qmap.n_pair_cohorts + 5 * qmap.n_blocks + 12 * qmap.n_pair_cohorts)
    work = 0
    for (cid=1; cid<=a.n_cells; cid++) {
        work = work + ncoh[m.cells[cid].cohort_index] * ((p + 2)^2 + d * (dq + 2) + 2 * p)
    }
    for (ell=1; ell<=a.n_links; ell++) {
        ret = ret + 64 + 5 * p + 4 * q + a.links[ell].n_ctrl * (8 + p) + dq * (1 + a.links[ell].n_ctrl)
        work = work + ncoh[a.links[ell].target_ci] * (p + 4)
        for (k=1; k<=a.links[ell].n_ctrl; k++) {
            work = work + ncoh[a.links[ell].ctrl_ci[k]] * (3 * p + 3 * q + d * (dq + 2) + 2 * dq + 12)
        }
    }
    for (bid=1; bid<=a.n_blocks; bid++) {
        ret = (ret + 32 + p + rows(a.blocks[bid].links) +
            a.blocks[bid].n_ctrl * (6 + 3 * p + p * p + 2 * dq))
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
            work = work + ncoh[a.blocks[bid].ctrl_ci[k]] * (4 * (p + 2)^2 + d * (dq + 2) + p * dq)
        }
    }
    ret = ret + 4 * p * p + 4 * p + a.n_cells * p * (p + 1) + a.n_groups * p * p + p * nq
    rdr_bytes = 8 * (a.N * (d + 7) + 2 * qr.pop.R + 18 * a.N + 4 * G)
    ret = 8 * ret + ib.retained_bytes + rdr_bytes
    sc_scan = cmax * (3 * p + q + d + 3 * dq + 24) + (p + 2) * (p + 2)
    sc_or = 2 * (p + 1) * a.n_blocks * (3 * p + 8) + a.n_cells * (p + 2) * (p + 4) + 2 * cmax * (p + 4)
    nlb = min((a.link_batch, a.n_links))
    nbb = min((nlb, a.n_blocks))
    ntb = min((nlb, a.n_trt))
    M_adj = nlb * (1 + a.has_ps) + (a.has_or ? nbb : 0)
    M_sel = nlb * (1 + a.has_ps) + (a.has_or ? nbb + (a.outcome_pooling == "separate" ? nbb * p : p) : 0)
    sc_rep = (cmax * (a.n_recipes * G + 4 * q * ntb + 13 * p + 6 * nbb + M_adj + M_sel + a.n_cells + 64) +
        a.n_cells * M_adj * (1 + p) + a.n_groups * p * M_adj + nq * M_sel)
    scr = 8 * max((sc_scan, sc_or, sc_rep, sb.stream_workspace_doubles)) + 4096
    work = 16 * work + 16 * sb.work_units * (a.glm_maxit + 2)
    _cdid_imp_check_budget(a.budget, ret, scr, work, "augmented backend preflight")
    a.retained_bytes = ret
    a.scratch_bytes = scr
    a.work_ops = work
}

/* fit: raw-W observation probabilities on F through the selection leaf
   with the same GLM controls as the roster logits (one source glm_control). */
void _cdid_aug_q_fit(struct cdid_aug_backend scalar a, struct cdid_sel_reader scalar qr,
    struct cdid_sel_cell_map scalar qmap)
{
    real scalar cid
    a.qfit = cdid_sel_fit(qr, qmap, a.glm_eps, a.glm_maxit)
    a.qum = cdid_sel_up_map_build(a.qfit)
    a.n_q_params = a.qfit.total_params
    a.n_q_fits = a.qfit.n_blocks
    for (cid=1; cid<=a.n_cells; cid++) {
        a.cells[cid].q_pc = a.q_of_cell[cid]
        a.cells[cid].theta_ids = _cdid_sel_pc_global_idx(a.qfit, a.qum, a.q_of_cell[cid])
    }
}

/* Step 4: augmentation scan of every cell on all full-F rows of its cohort (source cells_fitted
   additions): Yaug_one = H_N[v Z], Yaug_z = H_N[v Z z_raw], MW_one = H_N[v pi_full],
   MW_d = H_N[v pi_full d_imp], MD = M - MW (signed derivative loads, never fit weights)
   and Rq_one = -H_N[v pi_full e T_logit], with v = W a_k and pi_full = S/max(q, 1e-6). */
void _cdid_aug_cell_moments(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    struct _cdid_aug_rcc rowvector rcc, struct _cdid_aug_cev rowvector cev, real scalar use_cache)
{
    real scalar ci_ord, cid, ci, r, slot, ak, gid, pg, L, n_coh, lo, hi, p, dq, nobs
    real colvector c_ord, coh, rch, wr, v, vp, jh, zr, obs, qv, piv, ev
    real rowvector kept, allc
    real matrix Zr, Tm
    p = a.p_raw
    allc = _cdid_impb_seq(1, p)
    c_ord = _cdid_aug_order_cells(m, (1::a.n_cells))
    for (ci_ord=1; ci_ord<=rows(c_ord); ci_ord++) {
        cid = c_ord[ci_ord]
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        slot = (r - 1) * a.G + ci
        ak = _cdid_impb_ak(ib, m, r, ci)
        gid = m.cells[cid].group_id
        kept = m.groups[gid].kept
        pg = m.groups[gid].p
        L = m.cells[cid].L
        dq = rows(a.cells[cid].theta_ids)
        a.cells[cid].Yaug_one = 0
        a.cells[cid].Yaug_z = J(1, p, 0)
        a.cells[cid].MW_one = 0
        a.cells[cid].MW_d = J(1, pg, 0)
        a.cells[cid].Rq_one = J(dq, 1, 0)
        a.cells[cid].q_min = .
        a.cells[cid].q_max = .
        a.cells[cid].n_floor = 0
        a.cells[cid].w_max = 0
        nobs = 0
        coh = _cdid_impb_crows(ib, m, r, ci)
        n_coh = rows(coh)
        if (use_cache & n_coh <= a.chunk_rows) {
            if (!rcc[slot].ready) {
                wr = _cdid_impb_w(ib, m, r, coh, 1)
                rcc[slot].wr = wr
                rcc[slot].v = wr :* ak
                rcc[slot].Zr = (p > 0 ? _cdid_impb_z(ib, m, r, coh, allc, 1) : J(n_coh, 0, 0))
                rcc[slot].ready = 1
                rcc[slot].has_x = 0
            }
            Zr = rcc[slot].Zr
            v = rcc[slot].v
            jh = _cdid_impb_pred_z(m, cid, n_coh, Zr)
            _cdid_aug_cell_rows_z(a, m, qr, cid, coh, 1, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
            cev[cid].zr = zr
            cev[cid].obs = obs
            cev[cid].piv = piv
            cev[cid].ev = ev
            cev[cid].Tm = Tm
            cev[cid].ready = 1
            a.cells[cid].Yaug_one = quadsum(v :* zr)
            if (p > 0) a.cells[cid].Yaug_z = quadcross(v :* zr, Zr)
            if (rows(obs) > 0) {
                vp = v[obs] :* piv
                a.cells[cid].MW_one = quadsum(vp)
                if (pg > 0) a.cells[cid].MW_d = quadcross(vp, L :* Zr[obs, kept])
                if (dq > 0) a.cells[cid].Rq_one = -quadcross(Tm, vp :* ev)
                a.cells[cid].q_min = min(qv)
                a.cells[cid].q_max = max(qv)
                a.cells[cid].n_floor = sum(qv :<= 1e-6)
                a.cells[cid].w_max = max((0 \ vp))
                nobs = rows(obs)
            }
        }
        else {
            for (lo=1; lo<=n_coh; lo=lo+a.chunk_rows) {
                hi = min((n_coh, lo + a.chunk_rows - 1))
                rch = coh[lo..hi]
                Zr = (p > 0 ? _cdid_impb_z(ib, m, r, rch, allc, 1) : J(rows(rch), 0, 0))
                jh = _cdid_impb_pred_z(m, cid, rows(rch), Zr)
                _cdid_aug_cell_rows_z(a, m, qr, cid, rch, 1, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
                v = _cdid_impb_w(ib, m, r, rch, 1) :* ak
                a.cells[cid].Yaug_one = a.cells[cid].Yaug_one + quadsum(v :* zr)
                if (p > 0)
                    a.cells[cid].Yaug_z = a.cells[cid].Yaug_z + quadcross(v :* zr, Zr)
                if (rows(obs) == 0) continue
                vp = v[obs] :* piv
                a.cells[cid].MW_one = a.cells[cid].MW_one + quadsum(vp)
                if (pg > 0)
                    a.cells[cid].MW_d = a.cells[cid].MW_d + quadcross(vp, L :* Zr[obs, kept])
                if (dq > 0) a.cells[cid].Rq_one = a.cells[cid].Rq_one - quadcross(Tm, vp :* ev)
                a.cells[cid].q_min = (nobs == 0 ? min(qv) : min((a.cells[cid].q_min \ qv)))
                a.cells[cid].q_max = (nobs == 0 ? max(qv) : max((a.cells[cid].q_max \ qv)))
                a.cells[cid].n_floor = a.cells[cid].n_floor + sum(qv :<= 1e-6)
                a.cells[cid].w_max = max((a.cells[cid].w_max \ vp))
                nobs = nobs + rows(obs)
            }
        }
        if (nobs != m.cells[cid].n_pair)
            _cdid_aug_fail("Observed-pair support of cell " + m.cells[cid].name + " changed between the imputation fit and the augmentation scan.", 498)
        a.cells[cid].n_obs = nobs
        a.cells[cid].Yaug_one = a.cells[cid].Yaug_one / a.N
        a.cells[cid].Yaug_z = a.cells[cid].Yaug_z :/ a.N
        a.cells[cid].MW_one = a.cells[cid].MW_one / a.N
        a.cells[cid].MW_d = a.cells[cid].MW_d :/ a.N
        a.cells[cid].Rq_one = a.cells[cid].Rq_one :/ a.N
        a.cells[cid].MD_one = ib.cell_M_one[cid] - a.cells[cid].MW_one
        a.cells[cid].MD_d = (pg > 0 ? (L :* ib.cell_M_z[cid, kept]) - a.cells[cid].MW_d : J(1, 0, 0))
        if (missing(a.cells[cid].Yaug_one) | missing(a.cells[cid].Yaug_z) | missing(a.cells[cid].MW_one) |
            missing(a.cells[cid].MW_d) | missing(a.cells[cid].MD_one) | missing(a.cells[cid].MD_d) |
            missing(a.cells[cid].Rq_one))
            _cdid_aug_fail("Nonfinite augmented cell moments.", 504)
    }
}

/* Maximum full-F column scale of the proposed slope columns over every used support
   (control cells of all blocks and target cells of all links), as or_pooled_fit_blocks. */
real rowvector _cdid_aug_or_colmax(struct cdid_aug_backend scalar a,
    struct cdid_imputation_model scalar m, real rowvector proposed)
{
    real colvector used
    real rowvector out
    real scalar cid, bid, k, ell, row
    used = J(a.n_cells, 1, 0)
    for (bid=1; bid<=a.n_blocks; bid++) {
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) used[a.blocks[bid].ctrl_cell[k]] = 1
    }
    for (ell=1; ell<=a.n_links; ell++) used[a.links[ell].target_cell] = 1
    out = J(1, cols(proposed), 0)
    for (cid=1; cid<=a.n_cells; cid++) {
        if (!used[cid]) continue
        row = (m.cells[cid].recipe - 1) * m.n_cohorts + m.cells[cid].cohort_index
        out = colmax(out \ m.cohort_colmax_F[row, proposed])
    }
    return(out)
}

/* Strict reconstruction of proposed OR slope columns on actual original rows: every full-F
   cohort row of each block control cell (control role) and of each consuming link target
   cell (treated role), with the block constants (or_pooled_alias_check). Design-only. */
void _cdid_aug_or_alias(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, real rowvector retained, real rowvector proposed,
    real matrix map, real matrix constants_Z, real rowvector colmax_prop,
    real rowvector ctl_err, real rowvector trt_err)
{
    real scalar n_prop, cid, bid, ell, t, bb, ci, r, n_coh, lo, hi, j, den, ratio, is_c, is_t
    real colvector cb, tb, ub, coh, rch, u, s0
    real matrix Zs, u_mat, s0_mat
    n_prop = cols(proposed)
    ctl_err = J(1, n_prop, 0)
    trt_err = J(1, n_prop, 0)
    a.or_rows_checked = 0
    for (cid=1; cid<=a.n_cells; cid++) {
        cb = J(0, 1, 0)
        for (bid=1; bid<=a.n_blocks; bid++) {
            if (_cdid_impb_pos(a.blocks[bid].ctrl_cell, cid) > 0) cb = cb \ bid
        }
        tb = J(0, 1, 0)
        for (ell=1; ell<=a.n_links; ell++) {
            if (a.links[ell].target_cell == cid) tb = tb \ a.links[ell].or_block
        }
        if (rows(cb) + rows(tb) == 0) continue
        ub = uniqrows(cb \ tb)
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        coh = _cdid_impb_crows(ib, m, r, ci)
        n_coh = rows(coh)
        for (lo=1; lo<=n_coh; lo=lo+a.chunk_rows) {
            hi = min((n_coh, lo + a.chunk_rows - 1))
            rch = coh[lo..hi]
            Zs = _cdid_impb_z(ib, m, r, rch, _cdid_impb_seq(1, a.p_raw), 1)
            u_mat = Zs[, proposed] - (cols(retained) ? Zs[, retained] * map : J(rows(Zs), n_prop, 0))
            s0_mat = abs(Zs[, proposed]) + (cols(retained) ? abs(Zs[, retained]) * abs(map) : J(rows(Zs), n_prop, 0))
            for (t=1; t<=rows(ub); t++) {
                bb = ub[t]
                is_c = (_cdid_impb_pos(cb, bb) > 0)
                is_t = (_cdid_impb_pos(tb, bb) > 0)
                for (j=1; j<=n_prop; j++) {
                    u = abs(u_mat[, j] :- constants_Z[bb, j])
                    s0 = s0_mat[, j] :+ abs(constants_Z[bb, j])
                    den = max((colmax_prop[j], 2.2250738585072014e-308))
                    ratio = max(u :/ rowmax((s0, J(rows(u), 1, den))))
                    if (ratio >= .) ratio = .
                    if (is_c & (ratio > ctl_err[j] | ratio >= .)) ctl_err[j] = ratio
                    if (is_t & (ratio > trt_err[j] | ratio >= .)) trt_err[j] = ratio
                }
            }
            a.or_rows_checked = a.or_rows_checked + rows(Zs)
        }
    }
}

/* Block-scoped maximum full-F column scale of proposed slope columns for separate OR:
   strictly control cells of block bid and target cells of links consuming block bid. */
real rowvector _cdid_aug_or_colmax_blk(struct cdid_aug_backend scalar a,
    struct cdid_imputation_model scalar m, real scalar bid, real rowvector proposed)
{
    real colvector used, lk_b
    real rowvector out
    real scalar cid, k, i, row
    used = J(a.n_cells, 1, 0)
    for (k=1; k<=a.blocks[bid].n_ctrl; k++) used[a.blocks[bid].ctrl_cell[k]] = 1
    lk_b = a.blocks[bid].links
    for (i=1; i<=rows(lk_b); i++) used[a.links[lk_b[i]].target_cell] = 1
    out = J(1, cols(proposed), 0)
    for (cid=1; cid<=a.n_cells; cid++) {
        if (!used[cid]) continue
        row = (m.cells[cid].recipe - 1) * m.n_cohorts + m.cells[cid].cohort_index
        out = colmax(out \ m.cohort_colmax_F[row, proposed])
    }
    return(out)
}

/* Block-scoped alias verification for separate OR: strictly block bid control cells
   (control role) and consuming link target cells (treated role). */
void _cdid_aug_or_alias_blk(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, real scalar bid, real rowvector retained,
    real rowvector proposed, real matrix map, real rowvector const_row,
    real rowvector colmax_prop, real rowvector ctl_err, real rowvector trt_err)
{
    real scalar n_prop, cid, i, ci, r, n_coh, lo, hi, j, den, ratio, is_c, is_t
    real colvector lk_b, coh, rch, u, s0
    real matrix Zs, u_mat, s0_mat
    n_prop = cols(proposed)
    ctl_err = J(1, n_prop, 0)
    trt_err = J(1, n_prop, 0)
    lk_b = a.blocks[bid].links
    for (cid=1; cid<=a.n_cells; cid++) {
        is_c = (_cdid_impb_pos(a.blocks[bid].ctrl_cell, cid) > 0)
        is_t = 0
        for (i=1; i<=rows(lk_b); i++) {
            if (a.links[lk_b[i]].target_cell == cid) {
                is_t = 1
                break
            }
        }
        if (!is_c & !is_t) continue
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        coh = _cdid_impb_crows(ib, m, r, ci)
        n_coh = rows(coh)
        for (lo=1; lo<=n_coh; lo=lo+a.chunk_rows) {
            hi = min((n_coh, lo + a.chunk_rows - 1))
            rch = coh[lo..hi]
            Zs = _cdid_impb_z(ib, m, r, rch, _cdid_impb_seq(1, a.p_raw), 1)
            u_mat = Zs[, proposed] - (cols(retained) ? Zs[, retained] * map : J(rows(Zs), n_prop, 0))
            s0_mat = abs(Zs[, proposed]) + (cols(retained) ? abs(Zs[, retained]) * abs(map) : J(rows(Zs), n_prop, 0))
            for (j=1; j<=n_prop; j++) {
                u = abs(u_mat[, j] :- const_row[j])
                s0 = s0_mat[, j] :+ abs(const_row[j])
                den = max((colmax_prop[j], 2.2250738585072014e-308))
                ratio = max(u :/ rowmax((s0, J(rows(u), 1, den))))
                if (ratio >= .) ratio = .
                if (is_c & (ratio > ctl_err[j] | ratio >= .)) ctl_err[j] = ratio
                if (is_t & (ratio > trt_err[j] | ratio >= .)) trt_err[j] = ratio
            }
            a.or_rows_checked = a.or_rows_checked + rows(Zs)
        }
    }
}

/* Step 5: counterfactual OR on the reconstructed response Z over deduplicated blocks
   (pooled common slope or separate per-block): response Z of the (control cohort, block pair)
   cell on every full-F row of that cohort, POSITIVE weight v = W a_k (never v/q or signed v D),
   design d = L z. */
void _cdid_aug_fit_or(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    struct _cdid_aug_rcc rowvector rcc, struct _cdid_aug_cev rowvector cev)
{
    struct _cdid_impb_cacc rowvector acc
    struct cdid_wls_result scalar wfit
    real scalar p, nb, bid, k, ci_ord, cid, ci, r, slot, ak, n_coh, lo, hi, j, nn, scale, at, nPseudo
    real scalar thr, orig_norm, res_norm, rk, k_keep, exact_j, isfirst, wsum, L, n_b, tot_rss
    real colvector used_cids, c_ord, coh, rch, v, jh, zr, obs, qv, piv, ev, rr, ix, pseudo_y, pseudo_w, pseudo_block
    real colvector qy, coef, vv, yc, bvec, avec, mean_y, py_b, pw_b, blk_one
    real rowvector anchor, constb, coefrow, retained, proposed, aliased, aliased_pos
    real rowvector candidates, keep, colmax_prop, ctl_err, trt_err, allc, md_raw_b, const_row
    real matrix Zr, D, Fb, Q, R1, V, tau, pivot, res_fac, A_plus, A_minus, pseudo_X, C, Tm, px_b, Cb
    real matrix mean_D_raw, map, constants_Z, map_full, Z_keep, mean_D, Rfit

    p = a.p_raw
    nb = a.n_blocks
    allc = _cdid_impb_seq(1, p)
    acc = _cdid_impb_cacc(1, a.n_cells)
    for (cid=1; cid<=a.n_cells; cid++) {
        acc[cid].used = 0
        acc[cid].n = 0
        acc[cid].wsum = 0
        acc[cid].F = J(0, p + 2, 0)
        acc[cid].anchor = J(1, p, 0)
        acc[cid].constant = J(1, p, 1)
    }
    for (bid=1; bid<=nb; bid++) {
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) acc[a.blocks[bid].ctrl_cell[k]].used = 1
    }
    used_cids = J(0, 1, 0)
    for (cid=1; cid<=a.n_cells; cid++) {
        if (acc[cid].used) used_cids = used_cids \ cid
    }
    c_ord = _cdid_aug_order_cells(m, used_cids)

    /* 5a: bounded streaming of each used control cell: sqrt(v) [1, L z, Z] summaries. */
    for (ci_ord=1; ci_ord<=rows(c_ord); ci_ord++) {
        cid = c_ord[ci_ord]
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        slot = (r - 1) * a.G + ci
        L = m.cells[cid].L
        ak = _cdid_impb_ak(ib, m, r, ci)
        coh = _cdid_impb_crows(ib, m, r, ci)
        n_coh = rows(coh)
        if (cev[cid].ready & rcc[slot].ready & n_coh > 0) {
            v = rcc[slot].v
            Zr = rcc[slot].Zr
            zr = cev[cid].zr
            D = L :* Zr
            if (missing(D) | missing(zr) | missing(v))
                _cdid_aug_fail("Nonfinite augmented counterfactual OR support rows.", 504)
            if (p > 0) {
                acc[cid].anchor = D[1, ]
                acc[cid].constant = colmin(D :== acc[cid].anchor)
            }
            acc[cid].F = _cdid_imp_compress((J(n_coh, 1, 1), D, zr) :* sqrt(v))
            acc[cid].n = n_coh
            acc[cid].wsum = quadsum(v)
        }
        else {
            isfirst = 1
            for (lo=1; lo<=n_coh; lo=lo+a.chunk_rows) {
                hi = min((n_coh, lo + a.chunk_rows - 1))
                rch = coh[lo..hi]
                v = _cdid_impb_w(ib, m, r, rch, 1) :* ak
                Zr = (p > 0 ? _cdid_impb_z(ib, m, r, rch, allc, 1) : J(rows(rch), 0, 0))
                jh = _cdid_impb_pred_z(m, cid, rows(rch), Zr)
                _cdid_aug_cell_rows_z(a, m, qr, cid, rch, 0, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
                D = L :* Zr
                if (missing(D) | missing(zr) | missing(v))
                    _cdid_aug_fail("Nonfinite augmented counterfactual OR support rows.", 504)
                if (p > 0) {
                    if (isfirst) acc[cid].anchor = D[1, ]
                    acc[cid].constant = acc[cid].constant :& colmin(D :== acc[cid].anchor)
                }
                isfirst = 0
                acc[cid].F = _cdid_imp_compress(acc[cid].F \ ((J(rows(rch), 1, 1), D, zr) :* sqrt(v)))
                acc[cid].n = acc[cid].n + rows(rch)
                acc[cid].wsum = acc[cid].wsum + quadsum(v)
            }
        }
        if (acc[cid].n < 1 | !(acc[cid].wsum > 0) | acc[cid].wsum >= .)
            _cdid_aug_fail("Augmented OR control cell " + m.cells[cid].name + " has no positive synthetic full-F support.", 459)
    }

    /* 5b: merge supports per block; exact pseudo rows reproduce block means and crossproducts. */
    nPseudo = 2 * (p + 1) * nb
    pseudo_X = J(nPseudo, p, 0)
    pseudo_y = J(nPseudo, 1, 0)
    pseudo_w = J(nPseudo, 1, 0)
    pseudo_block = J(nPseudo, 1, 0)
    at = 0
    for (bid=1; bid<=nb; bid++) {
        Fb = J(0, p + 2, 0)
        wsum = 0
        constb = J(1, p, 1)
        anchor = J(1, p, 0)
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
            cid = a.blocks[bid].ctrl_cell[k]
            Fb = Fb \ acc[cid].F
            wsum = wsum + acc[cid].wsum
            if (p > 0) {
                if (k == 1) anchor = acc[cid].anchor
                constb = constb :& acc[cid].constant :& (acc[cid].anchor :== anchor)
            }
        }
        Fb = _cdid_imp_compress(Fb)
        a.blocks[bid].weight_sum = wsum
        a.blocks[bid].mass = wsum / a.N
        if (!(a.blocks[bid].mass > 0) | a.blocks[bid].mass >= .)
            _cdid_aug_fail("Augmented OR block " + a.blocks[bid].name + " has nonpositive synthetic control mass.", 459)
        hqrd(Fb[, 1], Q=., tau=., R1=.)
        V = hqrdmultq(Q, tau, Fb[, 2..(p + 2)], 1)
        coefrow = V[1, ] / R1[1, 1]
        res_fac = V[2..rows(V), ]
        if (a.intercept & p > 0) {
            for (j=1; j<=p; j++) {
                if (constb[j]) {
                    coefrow[j] = anchor[j]
                    res_fac[, j] = J(rows(res_fac), 1, 0)
                }
            }
        }
        nn = rows(res_fac)
        scale = sqrt(nn / wsum)
        A_plus = coefrow :+ scale :* res_fac
        A_minus = coefrow :- scale :* res_fac
        rr = (at + 1 .. at + 2 * nn)'
        if (p > 0) pseudo_X[rr, ] = A_plus[, 1..p] \ A_minus[, 1..p]
        pseudo_y[rr] = A_plus[, p + 1] \ A_minus[, p + 1]
        pseudo_w[rr] = J(2 * nn, 1, wsum / (2 * nn))
        pseudo_block[rr] = J(2 * nn, 1, bid)
        at = at + 2 * nn
    }
    acc = _cdid_impb_cacc(1, 1)

    if (a.outcome_pooling == "separate") {
        thr = 64 * epsilon(1) * max((1, p + a.intercept))
        a.or_rows_checked = 0
        a.or_threshold = thr
        tot_rss = 0
        for (bid=1; bid<=nb; bid++) {
            ix = _cdid_impb_ix(pseudo_block :== bid)
            n_b = rows(ix)
            px_b = pseudo_X[ix, ]
            py_b = pseudo_y[ix]
            pw_b = pseudo_w[ix]
            blk_one = J(n_b, 1, 1)
            Cb = px_b
            md_raw_b = J(1, p, 0)
            if (a.intercept & p > 0) {
                md_raw_b = quadcross(pw_b, px_b) / sum(pw_b)
                Cb = px_b :- md_raw_b
                for (j=1; j<=p; j++) {
                    if (all(px_b[, j] :== px_b[1, j])) {
                        Cb[, j] = J(n_b, 1, 0)
                        md_raw_b[j] = px_b[1, j]
                    }
                }
            }
            Cb = Cb :* sqrt(pw_b)
            retained = J(1, 0, 0)
            proposed = J(1, 0, 0)
            for (j=1; j<=p; j++) {
                vv = Cb[, j]
                orig_norm = sqrt(quadcross(vv, vv))
                if (cols(retained) > 0) {
                    hqrdp(Cb[, retained], Q=., tau=., R1=., pivot=.)
                    qy = hqrdmultq1t(Q, tau, vv)
                    coef = solveupper(R1, qy, rk=., 0)
                    coef = coef[invorder(pivot)]
                    vv = vv - Cb[, retained] * coef
                }
                res_norm = sqrt(quadcross(vv, vv))
                if (orig_norm == 0 | res_norm <= a.rank_tol * orig_norm) proposed = (proposed, j)
                else retained = (retained, j)
            }
            aliased = J(1, 0, 0)
            aliased_pos = J(1, 0, 0)
            candidates = J(1, 0, 0)
            map = J(cols(retained), 0, 0)
            const_row = J(1, 0, 0)
            ctl_err = J(1, 0, 0)
            trt_err = J(1, 0, 0)
            if (cols(proposed) > 0) {
                map = J(cols(retained), cols(proposed), 0)
                if (cols(retained) > 0) {
                    hqrdp(Cb[, retained], Q=., tau=., R1=., pivot=.)
                    qy = hqrdmultq1t(Q, tau, Cb[, proposed])
                    map = solveupper(R1, qy, rk=., 0)
                    map = map[invorder(pivot), ]
                    if (!missing(map)) {
                        qy = hqrdmultq1t(Q, tau, Cb[, proposed] - Cb[, retained] * map)
                        V = solveupper(R1, qy, rk=., 0)
                        if (!missing(V)) map = map + V[invorder(pivot), ]
                    }
                }
                const_row = J(1, cols(proposed), 0)
                if (a.intercept) {
                    const_row = (md_raw_b[proposed] -
                        (cols(retained) ? md_raw_b[retained] * map : J(1, cols(proposed), 0))) :/ a.blocks[bid].L
                }
                colmax_prop = _cdid_aug_or_colmax_blk(a, m, bid, proposed)
                _cdid_aug_or_alias_blk(a, ib, m, bid, retained, proposed, map, const_row, colmax_prop, ctl_err, trt_err)
                for (j=1; j<=cols(proposed); j++) {
                    exact_j = (ctl_err[j] <= thr & ctl_err[j] < .)
                    if (exact_j & (trt_err[j] > thr | trt_err[j] >= .))
                        _cdid_aug_fail("Control design of augmented counterfactual OR block " + a.blocks[bid].name + " has an exact null direction that does not reconstruct the required treated predictions; no ridge or singular-contrast solver is used; model unchanged. This computation failure does not establish population ATT nonidentification.", 459)
                    if (exact_j) {
                        aliased = (aliased, proposed[j])
                        aliased_pos = (aliased_pos, j)
                    }
                    else candidates = (candidates, proposed[j])
                }
                keep = (cols(retained) + cols(candidates) > 0 ? sort((retained, candidates)', 1)' : J(1, 0, 0))
                if (!a.intercept & cols(keep) == 0)
                    _cdid_aug_fail("Augmented OR block " + a.blocks[bid].name + " design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
                map_full = J(cols(keep), cols(aliased), 0)
                if (cols(retained) > 0 & cols(aliased) > 0) {
                    for (j=1; j<=cols(aliased); j++) {
                        for (k=1; k<=cols(retained); k++) {
                            map_full[_cdid_impb_pos(keep, retained[k]), j] = map[k, aliased_pos[j]]
                        }
                    }
                }
                map = map_full
                const_row = (cols(aliased) > 0 ? const_row[aliased_pos] : J(1, 0, 0))
            }
            else {
                keep = retained
                if (!a.intercept & cols(keep) == 0)
                    _cdid_aug_fail("Augmented OR block " + a.blocks[bid].name + " design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
            }
            k_keep = cols(keep)
            if (cols(candidates) == 0) {
                wfit = cdid_pooled_wls_fit((k_keep > 0 ? px_b[, keep] : J(n_b, 0, 0)),
                    py_b, pw_b, blk_one, blk_one, a.N,
                    a.intercept, J(0, k_keep, 0), J(0, 1, 0), a.rank_tol)
                if (cols(wfit.kept) != k_keep)
                    _cdid_aug_fail("Augmented OR block final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
                bvec = wfit.beta
                avec = wfit.alpha
                Rfit = wfit.R
                pivot = wfit.pivot
                a.blocks[bid].or_rcond = wfit.rcond_R
                mean_D = (a.intercept & k_keep > 0 ? wfit.mean_X : J(1, k_keep, 0))
                mean_y = wfit.mean_y
                a.blocks[bid].or_rss = wfit.rss
            }
            else {
                mean_y = J(1, 1, sum(pw_b :* py_b) / sum(pw_b))
                yc = (a.intercept ? py_b :- mean_y[1] : py_b) :* sqrt(pw_b)
                Z_keep = Cb[, keep]
                hqrdp(Z_keep, Q=., tau=., R1=., pivot=.)
                if (rows(R1) != k_keep | missing(R1) |
                    min(abs(diagonal(R1))) <= a.rank_tol * max(abs(diagonal(R1))))
                    _cdid_aug_fail("Augmented OR block final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
                qy = hqrdmultq1t(Q, tau, yc)
                coef = solveupper(R1, qy, rk=., 0)
                if (rk != k_keep | missing(coef)) _cdid_aug_fail("Nonfinite augmented OR block solve.", 504)
                bvec = coef[invorder(pivot)]
                Rfit = R1
                a.blocks[bid].or_rcond = _cdid_nuisance_rcond(R1)
                mean_D = (a.intercept ? md_raw_b[keep] : J(1, k_keep, 0))
                avec = (a.intercept ? mean_y - mean_D * bvec : J(1, 1, 0))
                vv = py_b - avec[blk_one] - px_b[, keep] * bvec
                a.blocks[bid].or_rss = sum(pw_b :* vv:^2)
            }
            if (missing(bvec) | missing(avec))
                _cdid_aug_fail("Nonfinite augmented counterfactual OR block coefficients.", 504)
            a.blocks[bid].or_p = k_keep
            a.blocks[bid].or_keep = keep
            a.blocks[bid].or_beta = bvec
            a.blocks[bid].or_R = Rfit
            a.blocks[bid].or_pivot = pivot
            a.blocks[bid].or_alias_map = map
            a.blocks[bid].or_alias_const = const_row
            a.blocks[bid].or_adj_ready = 0
            a.blocks[bid].alpha = (a.intercept ? avec[1] : 0)
            a.blocks[bid].mean_d = (a.intercept & k_keep > 0 ? mean_D[1, ] : J(1, k_keep, 0))
            a.blocks[bid].mean_y = mean_y[1]
            tot_rss = tot_rss + a.blocks[bid].or_rss
        }
        a.or_p = 0
        a.or_rss = tot_rss
        a.or_rcond = .
        a.or_keep = J(1, 0, 0)
        a.or_retained = J(1, 0, 0)
        a.or_proposed = J(1, 0, 0)
        a.or_dropped = J(1, 0, 0)
        a.or_candidates = J(1, 0, 0)
        a.or_ctl_scaled = J(1, 0, 0)
        a.or_trt_scaled = J(1, 0, 0)
        a.or_alias_map = J(0, 0, 0)
        a.or_alias_const = J(nb, 0, 0)
        a.or_beta = J(0, 1, 0)
        a.or_R = J(0, 0, 0)
        a.or_pivot = J(1, 0, 0)
        a.n_or_fits = nb
        return
    }

    /* 5c: original-order near-dependency proposal on the partialled weighted pooled factor. */
    C = pseudo_X
    mean_D_raw = J(nb, p, 0)
    for (bid=1; bid<=nb; bid++) {
        ix = _cdid_impb_ix(pseudo_block :== bid)
        if (a.intercept & p > 0) {
            mean_D_raw[bid, ] = quadcross(pseudo_w[ix], pseudo_X[ix, ]) / sum(pseudo_w[ix])
            C[ix, ] = pseudo_X[ix, ] :- mean_D_raw[bid, ]
            for (j=1; j<=p; j++) {
                if (all(pseudo_X[ix, j] :== pseudo_X[ix[1], j])) {
                    C[ix, j] = J(rows(ix), 1, 0)
                    mean_D_raw[bid, j] = pseudo_X[ix[1], j]
                }
            }
        }
    }
    C = C :* sqrt(pseudo_w)
    retained = J(1, 0, 0)
    proposed = J(1, 0, 0)
    for (j=1; j<=p; j++) {
        vv = C[, j]
        orig_norm = sqrt(quadcross(vv, vv))
        if (cols(retained) > 0) {
            hqrdp(C[, retained], Q=., tau=., R1=., pivot=.)
            qy = hqrdmultq1t(Q, tau, vv)
            coef = solveupper(R1, qy, rk=., 0)
            coef = coef[invorder(pivot)]
            vv = vv - C[, retained] * coef
        }
        res_norm = sqrt(quadcross(vv, vv))
        if (orig_norm == 0 | res_norm <= a.rank_tol * orig_norm) proposed = (proposed, j)
        else retained = (retained, j)
    }

    /* 5d: exact aliases must reconstruct control AND required treated full-F rows. */
    thr = 64 * epsilon(1) * max((1, p + a.intercept))
    aliased = J(1, 0, 0)
    aliased_pos = J(1, 0, 0)
    candidates = J(1, 0, 0)
    map = J(cols(retained), 0, 0)
    constants_Z = J(nb, 0, 0)
    ctl_err = J(1, 0, 0)
    trt_err = J(1, 0, 0)
    a.or_rows_checked = 0
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
        constants_Z = J(nb, cols(proposed), 0)
        if (a.intercept) {
            for (bid=1; bid<=nb; bid++) {
                constants_Z[bid, ] = (mean_D_raw[bid, proposed] -
                    (cols(retained) ? mean_D_raw[bid, retained] * map : J(1, cols(proposed), 0))) :/ a.blocks[bid].L
            }
        }
        colmax_prop = _cdid_aug_or_colmax(a, m, proposed)
        _cdid_aug_or_alias(a, ib, m, retained, proposed, map, constants_Z, colmax_prop, ctl_err, trt_err)
        for (j=1; j<=cols(proposed); j++) {
            exact_j = (ctl_err[j] <= thr & ctl_err[j] < .)
            if (exact_j & (trt_err[j] > thr | trt_err[j] >= .))
                _cdid_aug_fail("Pooled control design of the augmented counterfactual OR has an exact null direction that does not reconstruct the required treated predictions; no ridge or singular-contrast solver is used; model unchanged. This computation failure does not establish population ATT nonidentification.", 459)
            if (exact_j) {
                aliased = (aliased, proposed[j])
                aliased_pos = (aliased_pos, j)
            }
            else candidates = (candidates, proposed[j])
        }
        keep = (cols(retained) + cols(candidates) > 0 ? sort((retained, candidates)', 1)' : J(1, 0, 0))
        if (!a.intercept & cols(keep) == 0)
            _cdid_aug_fail("Augmented OR design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
        map_full = J(cols(keep), cols(aliased), 0)
        if (cols(retained) > 0 & cols(aliased) > 0) {
            for (j=1; j<=cols(aliased); j++) {
                for (k=1; k<=cols(retained); k++) {
                    map_full[_cdid_impb_pos(keep, retained[k]), j] = map[k, aliased_pos[j]]
                }
            }
        }
        map = map_full
        constants_Z = (cols(aliased) > 0 ? constants_Z[, aliased_pos] : J(nb, 0, 0))
    }
    else {
        keep = retained
    }

    /* 5e: one solve of the common slope on the kept columns (QR rank rule). */
    k_keep = cols(keep)
    if (cols(candidates) == 0) {
        wfit = cdid_pooled_wls_fit((k_keep > 0 ? pseudo_X[, keep] : J(nPseudo, 0, 0)),
            pseudo_y, pseudo_w, pseudo_block, J(nPseudo, 1, 1), a.N,
            a.intercept, J(0, k_keep, 0), J(0, 1, 0), a.rank_tol)
        if (cols(wfit.kept) != k_keep)
            _cdid_aug_fail("Augmented OR final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
        bvec = wfit.beta
        avec = wfit.alpha
        Rfit = wfit.R
        pivot = wfit.pivot
        a.or_rcond = wfit.rcond_R
        mean_D = (a.intercept & k_keep > 0 ? wfit.mean_X : J(nb, k_keep, 0))
        mean_y = wfit.mean_y
        a.or_rss = wfit.rss
    }
    else {
        yc = pseudo_y
        mean_y = J(nb, 1, 0)
        for (bid=1; bid<=nb; bid++) {
            ix = _cdid_impb_ix(pseudo_block :== bid)
            mean_y[bid] = sum(pseudo_w[ix] :* pseudo_y[ix]) / sum(pseudo_w[ix])
            if (a.intercept) yc[ix] = pseudo_y[ix] :- mean_y[bid]
        }
        yc = yc :* sqrt(pseudo_w)
        Z_keep = C[, keep]
        hqrdp(Z_keep, Q=., tau=., R1=., pivot=.)
        if (rows(R1) != k_keep | missing(R1) |
            min(abs(diagonal(R1))) <= a.rank_tol * max(abs(diagonal(R1))))
            _cdid_aug_fail("Augmented OR final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
        qy = hqrdmultq1t(Q, tau, yc)
        coef = solveupper(R1, qy, rk=., 0)
        if (rk != k_keep | missing(coef)) _cdid_aug_fail("Nonfinite augmented OR solve.", 504)
        bvec = coef[invorder(pivot)]
        Rfit = R1
        a.or_rcond = _cdid_nuisance_rcond(R1)
        mean_D = (a.intercept ? mean_D_raw[, keep] : J(nb, k_keep, 0))
        avec = (a.intercept ? mean_y - mean_D * bvec : J(nb, 1, 0))
        vv = pseudo_y - avec[pseudo_block] - pseudo_X[, keep] * bvec
        a.or_rss = sum(pseudo_w :* vv:^2)
    }
    if (missing(bvec) | missing(avec))
        _cdid_aug_fail("Nonfinite augmented counterfactual OR coefficients.", 504)
    a.or_p = k_keep
    a.or_keep = keep
    a.or_retained = retained
    a.or_proposed = proposed
    a.or_dropped = aliased
    a.or_candidates = candidates
    a.or_ctl_scaled = ctl_err
    a.or_trt_scaled = trt_err
    a.or_alias_map = map
    a.or_alias_const = constants_Z
    a.or_threshold = thr
    a.or_beta = bvec
    a.or_R = Rfit
    a.or_pivot = pivot
    for (bid=1; bid<=nb; bid++) {
        a.blocks[bid].alpha = (a.intercept ? avec[bid] : 0)
        a.blocks[bid].mean_d = (a.intercept & k_keep > 0 ? mean_D[bid, ] : J(1, k_keep, 0))
        a.blocks[bid].mean_y = mean_y[bid]
    }
    a.n_or_fits = 1
}

/* Step 5 forcing of every block constituent (control cohort k), with r = Z - m and
   D = 1 - pi_full (signed; D = 1 off observed pairs):
   a_force_alpha = H_N[W_r I_k r], a_force_beta = H_N[W_r I_k r dtilde];
   jf_alpha_alpha = MD_one (intercept), jf_alpha_beta = MD_d (kept -> raw positions);
   jf_beta_alpha = H_N[v D dtilde] (intercept; nonzero even for a singleton block, because
   centering uses v while this forcing uses v D), jf_beta_beta = H_N[v D dtilde d_imp'];
   q_alpha merges -Rq_one by native theta ID per block, and the GLOBAL common-slope q forcing
   q_beta = H_N[dtilde v pi_full e T_logit] merges ALL blocks before the slope solve. */
void _cdid_aug_or_forcing(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    struct _cdid_aug_rcc rowvector rcc, struct _cdid_aug_cev rowvector cev)
{
    struct _cdid_aug_bkacc rowvector bka
    real scalar p, po, bid, k, nk, cid, gid, pg, ci, r, slot, ak, L, n_coh, n_r, lo, hi, alpha_b, dq, sep
    real scalar tot_bk, idx, ci_ord, need_T, h, use_cev
    real colvector c_ord, cell_used, cell_need_T, hits, coh, rch, wr, v, mh, res, omp, wD
    real colvector jh, zr, obs, qv, piv, ev, used, pos, gpos, beta_b, ak_vec
    real rowvector kept_g, keep_b, allc
    real matrix bk_tab, Zr, dor, dt, dimp, afb, jba, jab, jbb, Tm
    p = a.p_raw
    allc = _cdid_impb_seq(1, p)
    sep = (a.outcome_pooling == "separate")
    if (sep) {
        a.q_beta_ids = J(0, 1, 0)
        a.q_beta = J(0, 0, 0)
    }
    else {
        used = J(0, 1, 0)
        for (bid=1; bid<=a.n_blocks; bid++) used = used \ a.blocks[bid].ctrl_cell
        a.q_beta_ids = _cdid_aug_theta_union(a, uniqrows(used))
        a.q_beta = J(a.or_p, rows(a.q_beta_ids), 0)
    }
    tot_bk = 0
    for (bid=1; bid<=a.n_blocks; bid++) tot_bk = tot_bk + a.blocks[bid].n_ctrl
    bka = _cdid_aug_bkacc(1, tot_bk)
    bk_tab = J(tot_bk, 3, 0)
    cell_used = J(a.n_cells, 1, 0)
    cell_need_T = J(a.n_cells, 1, 0)
    idx = 0
    for (bid=1; bid<=a.n_blocks; bid++) {
        po = (sep ? a.blocks[bid].or_p : a.or_p)
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
            idx++
            cid = a.blocks[bid].ctrl_cell[k]
            gid = m.cells[cid].group_id
            pg = m.groups[gid].p
            dq = rows(a.cells[cid].theta_ids)
            bk_tab[idx, ] = (bid, k, cid)
            bka[idx].s_r = 0
            bka[idx].s_rd = J(1, po, 0)
            bka[idx].s_vd = J(1, po, 0)
            bka[idx].s_vdd = J(po, pg, 0)
            bka[idx].Rqd = J(po, dq, 0)
            cell_used[cid] = 1
            if (po > 0) cell_need_T[cid] = 1
        }
    }
    c_ord = _cdid_aug_order_cells(m, _cdid_impb_ix(cell_used :== 1))
    for (ci_ord=1; ci_ord<=rows(c_ord); ci_ord++) {
        cid = c_ord[ci_ord]
        hits = _cdid_impb_ix(bk_tab[, 3] :== cid)
        need_T = cell_need_T[cid]
        gid = m.cells[cid].group_id
        kept_g = m.groups[gid].kept
        pg = m.groups[gid].p
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        slot = (r - 1) * a.G + ci
        L = m.cells[cid].L
        ak = _cdid_impb_ak(ib, m, r, ci)
        dq = rows(a.cells[cid].theta_ids)
        coh = _cdid_impb_crows(ib, m, r, ci)
        n_coh = rows(coh)
        use_cev = (cev[cid].ready & rcc[slot].ready)
        for (lo=1; lo<=n_coh; lo=lo+(use_cev ? n_coh : a.chunk_rows)) {
            if (use_cev) {
                n_r = n_coh
                Zr = rcc[slot].Zr
                wr = rcc[slot].wr
                v = rcc[slot].v
                zr = cev[cid].zr
                obs = cev[cid].obs
                piv = cev[cid].piv
                ev = cev[cid].ev
                Tm = cev[cid].Tm
            }
            else {
                hi = min((n_coh, lo + a.chunk_rows - 1))
                rch = coh[lo..hi]
                n_r = rows(rch)
                Zr = (p > 0 ? _cdid_impb_z(ib, m, r, rch, allc, 1) : J(n_r, 0, 0))
                jh = _cdid_impb_pred_z(m, cid, n_r, Zr)
                _cdid_aug_cell_rows_z(a, m, qr, cid, rch, need_T, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
                wr = _cdid_impb_w(ib, m, r, rch, 1)
                v = wr :* ak
            }
            if (need_T) {
                omp = J(n_r, 1, 1)
                if (rows(obs) > 0) omp[obs] = 1 :- piv
                wD = v :* omp
                dimp = (pg > 0 ? L :* Zr[, kept_g] : J(n_r, 0, 0))
            }
            for (h=1; h<=rows(hits); h++) {
                idx = hits[h]
                bid = bk_tab[idx, 1]
                po = (sep ? a.blocks[bid].or_p : a.or_p)
                alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
                if (po > 0) {
                    keep_b = (sep ? a.blocks[bid].or_keep : a.or_keep)
                    beta_b = (sep ? a.blocks[bid].or_beta : a.or_beta)
                    dor = L :* Zr[, keep_b]
                    mh = alpha_b :+ dor * beta_b
                }
                else mh = J(n_r, 1, alpha_b)
                res = zr - mh
                bka[idx].s_r = bka[idx].s_r + quadsum(wr :* res)
                if (po > 0) {
                    dt = (a.intercept ? dor :- a.blocks[bid].mean_d : dor)
                    bka[idx].s_rd = bka[idx].s_rd + quadcross(wr :* res, dt)
                    bka[idx].s_vd = bka[idx].s_vd + quadcross(wD, dt)
                    if (pg > 0) bka[idx].s_vdd = bka[idx].s_vdd + quadcross(wD :* dt, dimp)
                    if (rows(obs) > 0 & dq > 0)
                        bka[idx].Rqd = bka[idx].Rqd - quadcross(dt[obs, ] :* (v[obs] :* piv :* ev), Tm)
                }
            }
        }
    }
    idx = 0
    for (bid=1; bid<=a.n_blocks; bid++) {
        po = (sep ? a.blocks[bid].or_p : a.or_p)
        nk = a.blocks[bid].n_ctrl
        a.blocks[bid].a_force_alpha = J(nk, 1, 0)
        a.blocks[bid].jf_alpha_alpha = J(nk, 1, 0)
        a.blocks[bid].q_alpha_ids = _cdid_aug_theta_union(a, a.blocks[bid].ctrl_cell)
        a.blocks[bid].q_alpha = J(rows(a.blocks[bid].q_alpha_ids), 1, 0)
        if (sep) {
            a.blocks[bid].q_beta_ids = _cdid_aug_theta_union(a, a.blocks[bid].ctrl_cell)
            a.blocks[bid].q_beta = J(po, rows(a.blocks[bid].q_beta_ids), 0)
        }
        else {
            a.blocks[bid].q_beta_ids = J(0, 1, 0)
            a.blocks[bid].q_beta = J(0, 0, 0)
        }
        ak_vec = J(nk, 1, 0)
        afb = J(po, nk, 0)
        jba = J(po, nk, 0)
        jab = J(nk, p, 0)
        jbb = J(po, nk * p, 0)
        for (k=1; k<=nk; k++) {
            idx++
            cid = a.blocks[bid].ctrl_cell[k]
            gid = m.cells[cid].group_id
            kept_g = m.groups[gid].kept
            pg = m.groups[gid].p
            ci = m.cells[cid].cohort_index
            r = m.cells[cid].recipe
            ak_vec[k] = _cdid_impb_ak(ib, m, r, ci)
            dq = rows(a.cells[cid].theta_ids)
            a.blocks[bid].a_force_alpha[k] = bka[idx].s_r / a.N
            a.blocks[bid].jf_alpha_alpha[k] = (a.intercept ? a.cells[cid].MD_one : 0)
            if (pg > 0) jab[k, kept_g] = a.cells[cid].MD_d
            if (po > 0) {
                afb[, k] = bka[idx].s_rd' / a.N
                if (a.intercept) jba[, k] = bka[idx].s_vd' / a.N
                if (pg > 0) jbb[, (k - 1) * p :+ kept_g] = bka[idx].s_vdd / a.N
            }
            if (dq > 0) {
                pos = _cdid_aug_find(a.blocks[bid].q_alpha_ids, a.cells[cid].theta_ids)
                a.blocks[bid].q_alpha[pos] = a.blocks[bid].q_alpha[pos] - a.cells[cid].Rq_one
                if (po > 0) {
                    if (sep) {
                        gpos = _cdid_aug_find(a.blocks[bid].q_beta_ids, a.cells[cid].theta_ids)
                        a.blocks[bid].q_beta[, gpos] = a.blocks[bid].q_beta[, gpos] - bka[idx].Rqd :/ a.N
                    }
                    else {
                        gpos = _cdid_aug_find(a.q_beta_ids, a.cells[cid].theta_ids)
                        a.q_beta[, gpos] = a.q_beta[, gpos] - bka[idx].Rqd :/ a.N
                    }
                }
            }
        }
        a.blocks[bid].a_force_beta = afb
        a.blocks[bid].jf_beta_alpha = jba
        a.blocks[bid].jf_alpha_beta = jab
        a.blocks[bid].jf_beta_beta = jbb
        if (sep) {
            a.blocks[bid].score_mean_a = sum(ak_vec :* a.blocks[bid].a_force_alpha)
            a.blocks[bid].score_mean_b = (po > 0 ? afb * ak_vec : J(0, 1, 0))
            if (missing(a.blocks[bid].q_beta) | missing(a.blocks[bid].score_mean_a) |
                missing(a.blocks[bid].score_mean_b))
                _cdid_aug_fail("Nonfinite augmented block-local OR forcing.", 504)
        }
        else {
            a.blocks[bid].score_mean_a = 0
            a.blocks[bid].score_mean_b = J(0, 1, 0)
        }
        if (missing(a.blocks[bid].a_force_alpha) | missing(afb) | missing(jba) | missing(jab) |
            missing(jbb) | missing(a.blocks[bid].q_alpha))
            _cdid_aug_fail("Nonfinite augmented OR forcing.", 504)
    }
    if (!sep & missing(a.q_beta)) _cdid_aug_fail("Nonfinite augmented global shared-slope selection forcing.", 504)
}

/* Step 6: link ratios with the reconstructed Z numerators and q-independent normalizers,
   Z-fitted OR terms (when reg/dr), propensity control terms (when ipw/dr), signed MD
   adjoint weights and merged theta-only selection loadings. */
void _cdid_aug_link_sums(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    struct _cdid_aug_rcc rowvector rcc, struct _cdid_aug_cev rowvector cev)
{
    struct _cdid_aug_lacc rowvector la
    real scalar ell, g, cg, rg, ag, L, d_plus, num_T, T, bid, alpha_b, B_gl, num_Bg
    real scalar tid, rtid, nk, k, ck, ckc, rk, slot, ak, gk, pk, n_coh, lo, hi, p, po, q, dq, sep
    real scalar d_minus, Cl, B_ol, tot_lk, idx, ci_ord, n_hit, n_tp, h, tp, n_r, use_cev, prev_bid
    real colvector c_ord, u_ckc, u_rtid, hit_rows, coh, rch, wr, v, wm, wod, wdot
    real colvector omp, mh, pos, beta_b, jh, zr, obs, qv, piv, ev, s_wmD_h
    real rowvector keptk, dbp, keep_b, allc, s_wm_tp, s_wmZ_tp, s_wmD_tp, s_wodZ_tp, s_wod_tp
    real matrix lk_tab, hits, Zr, X, vals, wm_m, wod_m, wdot_m, wmD_m, dimp, dor, Tm, s_dimpD_h, Rqo_h
    real matrix s_wdZX_tp, s_wdX_tp, s_dimpD_tp, Rqo_tp
    p = a.p_raw
    q = a.q
    allc = _cdid_impb_seq(1, p)
    sep = (a.has_or & a.outcome_pooling == "separate")
    if (a.has_ps) la = _cdid_aug_lacc(1, a.n_links)
    tot_lk = 0
    for (ell=1; ell<=a.n_links; ell++) {
        g = a.links[ell].target_cell
        cg = a.links[ell].target_ci
        rg = m.cells[g].recipe
        ag = _cdid_impb_ak(ib, m, rg, cg)
        L = a.links[ell].L
        d_plus = ib.cell_M_one[g]
        if (!(d_plus > 0) | d_plus >= .)
            _cdid_aug_fail("Link " + a.links[ell].link_id + " has nonpositive treated denominator d_{l,+}.", 459)
        num_T = a.cells[g].Yaug_one
        T = num_T / d_plus
        a.links[ell].d_plus = d_plus
        a.links[ell].num_T = num_T
        a.links[ell].T_gl = T
        a.links[ell].a_coeff_T_g = ((a.cells[g].Yaug_one / ag) - T * (ib.cell_M_one[g] / ag)) / d_plus
        a.links[ell].lam_T = a.cells[g].MD_one / d_plus
        a.links[ell].nu_T = a.cells[g].MD_d :/ d_plus
        a.links[ell].sel_T_ids = a.cells[g].theta_ids
        a.links[ell].sel_T = -a.cells[g].Rq_one :/ d_plus
        if (missing(T) | missing(a.links[ell].a_coeff_T_g) | missing(a.links[ell].lam_T) |
            missing(a.links[ell].nu_T) | missing(a.links[ell].sel_T))
            _cdid_aug_fail("Nonfinite augmented link treated summary.", 504)

        po = 0
        keep_b = J(1, 0, 0)
        beta_b = J(0, 1, 0)
        alpha_b = 0
        dbp = J(1, 0, 0)
        B_gl = .
        if (a.has_or) {
            bid = a.links[ell].or_block
            po = (sep ? a.blocks[bid].or_p : a.or_p)
            keep_b = (sep ? a.blocks[bid].or_keep : a.or_keep)
            beta_b = (sep ? a.blocks[bid].or_beta : a.or_beta)
            alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
            dbp = (po > 0 ? (L :* ib.cell_M_z[g, keep_b]) :/ d_plus : J(1, 0, 0))
            B_gl = alpha_b + (po > 0 ? dbp * beta_b : 0)
            num_Bg = d_plus * B_gl
            a.links[ell].d_bar_plus_or = dbp
            a.links[ell].B_gl = B_gl
            a.links[ell].num_Bg = num_Bg
            a.links[ell].a_coeff_Bg_g = ((num_Bg / ag) - B_gl * (d_plus / ag)) / d_plus
            if (missing(dbp) | missing(B_gl) | missing(a.links[ell].a_coeff_Bg_g))
                _cdid_aug_fail("Nonfinite augmented link treated OR summary.", 504)
        }
        else {
            a.links[ell].d_bar_plus_or = J(1, 0, 0)
            a.links[ell].B_gl = .
            a.links[ell].num_Bg = .
            a.links[ell].a_coeff_Bg_g = .
        }

        if (a.has_ps) {
            nk = a.links[ell].n_ctrl
            tot_lk = tot_lk + nk
            la[ell].s_wm = 0
            la[ell].s_wmZ = 0
            la[ell].s_wmm = 0
            la[ell].po = po
            la[ell].alpha_b = alpha_b
            la[ell].keep_b = keep_b
            la[ell].beta_b = beta_b
            la[ell].dbp = dbp
            la[ell].s_wmdor = J(1, po, 0)
            la[ell].s_wm_k = J(nk, 1, 0)
            la[ell].s_wodZ_k = J(nk, 1, 0)
            la[ell].s_wod_k = J(nk, 1, 0)
            la[ell].s_wodm_k = J(nk, 1, 0)
            la[ell].s_wdZX = J(1, q, 0)
            la[ell].s_wdX = J(1, q, 0)
            la[ell].s_wdmX = J(1, q, 0)
            la[ell].lamC = J(nk, 1, 0)
            la[ell].nuC = J(nk, p, 0)
            a.links[ell].sel_C_ids = _cdid_aug_theta_union(a, a.links[ell].ctrl_cell)
            a.links[ell].sel_C = J(rows(a.links[ell].sel_C_ids), 1, 0)
        }
        else {
            a.links[ell].d_minus = .
            a.links[ell].num_C = .
            a.links[ell].C_ol = .
            a.links[ell].num_Bo = .
            a.links[ell].B_ol = .
            a.links[ell].d_bar_minus_or = J(1, 0, 0)
            a.links[ell].d_bar_diff_or = J(1, 0, 0)
            a.links[ell].m_minus_k = J(0, 1, .)
            a.links[ell].lam_C_k = J(0, 1, .)
            a.links[ell].nu_C_k = J(0, 0, .)
            a.links[ell].a_coeff_C_k = J(0, 1, .)
            a.links[ell].a_coeff_Bo_k = J(0, 1, .)
            a.links[ell].denminus_gamma = J(1, 0, .)
            a.links[ell].g_C_gamma = J(1, 0, .)
            a.links[ell].g_Bo_gamma = J(1, 0, .)
            a.links[ell].g_dr_gamma = J(1, 0, .)
            a.links[ell].sel_C_ids = J(0, 1, 0)
            a.links[ell].sel_C = J(0, 1, 0)
            a.links[ell].mu_ipw = .
            a.links[ell].mu_reg = (a.has_or ? T - B_gl : .)
            a.links[ell].mu_dr = .
            a.links[ell].mu = (a.method == "ipw" ? a.links[ell].mu_ipw :
                (a.method == "reg" ? a.links[ell].mu_reg : a.links[ell].mu_dr))
            if (missing(a.links[ell].mu))
                _cdid_aug_fail("Nonfinite augmented link point estimate.", 504)
        }
    }
    if (!a.has_ps) return

    lk_tab = J(tot_lk, 4, 0)
    idx = 0
    for (ell=1; ell<=a.n_links; ell++) {
        tid = a.links[ell].trt_id
        rtid = ib.trts[tid].rep_tid
        for (k=1; k<=a.links[ell].n_ctrl; k++) {
            idx++
            lk_tab[idx, ] = (ell, k, a.links[ell].ctrl_cell[k], rtid)
        }
    }
    u_ckc = uniqrows(lk_tab[, 3])
    c_ord = _cdid_aug_order_cells(m, u_ckc)
    for (ci_ord=1; ci_ord<=rows(c_ord); ci_ord++) {
        ckc = c_ord[ci_ord]
        hit_rows = _cdid_impb_ix(lk_tab[, 3] :== ckc)
        n_hit = rows(hit_rows)
        hits = lk_tab[hit_rows, ]
        u_rtid = uniqrows(hits[, 4])
        n_tp = rows(u_rtid)
        for (h=1; h<=n_hit; h++) hits[h, 4] = _cdid_impb_pos(u_rtid, hits[h, 4])
        ck = m.cells[ckc].cohort_index
        rk = m.cells[ckc].recipe
        slot = (rk - 1) * a.G + ck
        L = m.cells[ckc].L
        ak = _cdid_impb_ak(ib, m, rk, ck)
        gk = m.cells[ckc].group_id
        keptk = m.groups[gk].kept
        pk = m.groups[gk].p
        dq = rows(a.cells[ckc].theta_ids)
        s_wmD_h = J(n_hit, 1, 0)
        s_dimpD_h = J(n_hit, pk, 0)
        Rqo_h = J(dq, n_hit, 0)
        coh = _cdid_impb_crows(ib, m, rk, ck)
        n_coh = rows(coh)
        use_cev = (cev[ckc].ready & rcc[slot].ready)
        for (lo=1; lo<=n_coh; lo=lo+(use_cev ? n_coh : a.chunk_rows)) {
            if (use_cev) {
                n_r = n_coh
                wr = rcc[slot].wr
                v = rcc[slot].v
                Zr = rcc[slot].Zr
                if (!rcc[slot].has_x) {
                    rcc[slot].X = _cdid_impb_x(ib, rk, coh, 1)
                    rcc[slot].has_x = 1
                }
                X = rcc[slot].X
                zr = cev[ckc].zr
                obs = cev[ckc].obs
                piv = cev[ckc].piv
                ev = cev[ckc].ev
                Tm = cev[ckc].Tm
            }
            else {
                hi = min((n_coh, lo + a.chunk_rows - 1))
                rch = coh[lo..hi]
                n_r = rows(rch)
                wr = _cdid_impb_w(ib, m, rk, rch, 1)
                v = wr :* ak
                X = _cdid_impb_x(ib, rk, rch, 1)
                Zr = (p > 0 ? _cdid_impb_z(ib, m, rk, rch, allc, 1) : J(n_r, 0, 0))
                jh = _cdid_impb_pred_z(m, ckc, n_r, Zr)
                _cdid_aug_cell_rows_z(a, m, qr, ckc, rch, 1, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
            }
            omp = J(n_r, 1, 1)
            if (rows(obs) > 0) omp[obs] = 1 :- piv
            dimp = (pk > 0 ? L :* Zr[, keptk] : J(n_r, 0, 0))
            wm_m = J(n_r, n_tp, 0)
            wod_m = J(n_r, n_tp, 0)
            wdot_m = J(n_r, n_tp, 0)
            wmD_m = J(n_r, n_tp, 0)
            for (tp=1; tp<=n_tp; tp++) {
                rtid = u_rtid[tp]
                vals = _cdid_aug_logit_vals(X * ib.trts[rtid].gamma)
                wm_m[, tp] = v :* vals[, 2]
                wod_m[, tp] = wr :* vals[, 2]
                wdot_m[, tp] = wm_m[, tp] :* vals[, 3]
                wmD_m[, tp] = wm_m[, tp] :* omp
            }
            s_wm_tp = quadcolsum(wm_m)
            s_wmZ_tp = quadcross(zr, wm_m)
            s_wmD_tp = quadcolsum(wmD_m)
            s_wodZ_tp = quadcross(zr, wod_m)
            s_wod_tp = quadcolsum(wod_m)
            s_wdZX_tp = quadcross(wdot_m :* zr, X)
            s_wdX_tp = quadcross(wdot_m, X)
            s_dimpD_tp = (pk > 0 ? quadcross(wmD_m, dimp) : J(n_tp, 0, 0))
            Rqo_tp = (rows(obs) > 0 & dq > 0 ? quadcross(Tm, wm_m[obs, ] :* (piv :* ev)) : J(dq, n_tp, 0))
            prev_bid = 0
            mh = J(0, 1, 0)
            dor = J(0, 0, 0)
            for (h=1; h<=n_hit; h++) {
                ell = hits[h, 1]
                k = hits[h, 2]
                tp = hits[h, 4]
                la[ell].s_wm = la[ell].s_wm + s_wm_tp[tp]
                la[ell].s_wmZ = la[ell].s_wmZ + s_wmZ_tp[tp]
                la[ell].s_wm_k[k] = la[ell].s_wm_k[k] + s_wm_tp[tp]
                s_wmD_h[h] = s_wmD_h[h] + s_wmD_tp[tp]
                if (pk > 0) s_dimpD_h[h, ] = s_dimpD_h[h, ] + s_dimpD_tp[tp, ]
                la[ell].s_wodZ_k[k] = la[ell].s_wodZ_k[k] + s_wodZ_tp[tp]
                la[ell].s_wod_k[k] = la[ell].s_wod_k[k] + s_wod_tp[tp]
                la[ell].s_wdZX = la[ell].s_wdZX + s_wdZX_tp[tp, ]
                la[ell].s_wdX = la[ell].s_wdX + s_wdX_tp[tp, ]
                if (a.has_or) {
                    bid = a.links[ell].or_block
                    po = la[ell].po
                    if (bid != prev_bid) {
                        alpha_b = la[ell].alpha_b
                        if (po > 0) {
                            dor = L :* Zr[, la[ell].keep_b]
                            mh = alpha_b :+ dor * la[ell].beta_b
                        }
                        else mh = J(n_r, 1, alpha_b)
                        prev_bid = bid
                    }
                    wm = wm_m[, tp]
                    wod = wod_m[, tp]
                    wdot = wdot_m[, tp]
                    if (po > 0) la[ell].s_wmdor = la[ell].s_wmdor + quadcross(wm, dor)
                    la[ell].s_wmm = la[ell].s_wmm + quadsum(wm :* mh)
                    la[ell].s_wodm_k[k] = la[ell].s_wodm_k[k] + quadsum(wod :* mh)
                    la[ell].s_wdmX = la[ell].s_wdmX + quadcross(wdot :* mh, X)
                }
                if (rows(obs) > 0 & dq > 0) Rqo_h[, h] = Rqo_h[, h] - Rqo_tp[, tp]
            }
        }
        for (h=1; h<=n_hit; h++) {
            ell = hits[h, 1]
            k = hits[h, 2]
            la[ell].lamC[k] = s_wmD_h[h]
            if (pk > 0) la[ell].nuC[k, keptk] = s_dimpD_h[h, ]
            if (dq > 0) {
                pos = _cdid_aug_find(a.links[ell].sel_C_ids, a.cells[ckc].theta_ids)
                a.links[ell].sel_C[pos] = a.links[ell].sel_C[pos] - Rqo_h[, h]
            }
        }
    }

    for (ell=1; ell<=a.n_links; ell++) {
        T = a.links[ell].T_gl
        B_gl = a.links[ell].B_gl
        B_ol = .
        d_minus = la[ell].s_wm / a.N
        if (!(d_minus > 0) | d_minus >= .)
            _cdid_aug_fail("Link " + a.links[ell].link_id + " has nonpositive control denominator d_{l,-}.", 459)
        Cl = (la[ell].s_wmZ / a.N) / d_minus
        a.links[ell].d_minus = d_minus
        a.links[ell].num_C = la[ell].s_wmZ / a.N
        a.links[ell].C_ol = Cl
        a.links[ell].m_minus_k = (la[ell].s_wm_k :/ a.N) :/ d_minus
        a.links[ell].lam_C_k = (la[ell].lamC :/ a.N) :/ d_minus
        a.links[ell].nu_C_k = (la[ell].nuC :/ a.N) :/ d_minus
        a.links[ell].sel_C = (a.links[ell].sel_C :/ a.N) :/ d_minus
        a.links[ell].a_coeff_C_k = ((la[ell].s_wodZ_k :/ a.N) - Cl :* (la[ell].s_wod_k :/ a.N)) :/ d_minus
        a.links[ell].denminus_gamma = la[ell].s_wdX / a.N
        a.links[ell].g_C_gamma = ((la[ell].s_wdZX / a.N) - Cl * (la[ell].s_wdX / a.N)) / d_minus
        if (missing(Cl) | missing(a.links[ell].m_minus_k) | missing(a.links[ell].lam_C_k) |
            missing(a.links[ell].nu_C_k) | missing(a.links[ell].sel_C) |
            missing(a.links[ell].a_coeff_C_k) | missing(a.links[ell].g_C_gamma))
            _cdid_aug_fail("Nonfinite augmented link control summary.", 504)
        if (a.has_or) {
            po = la[ell].po
            a.links[ell].d_bar_minus_or = (po > 0 ? (la[ell].s_wmdor / a.N) / d_minus : J(1, 0, 0))
            a.links[ell].d_bar_diff_or = la[ell].dbp - a.links[ell].d_bar_minus_or
            a.links[ell].num_Bo = la[ell].s_wmm / a.N
            B_ol = la[ell].alpha_b + (po > 0 ? a.links[ell].d_bar_minus_or * la[ell].beta_b : 0)
            a.links[ell].B_ol = B_ol
            a.links[ell].a_coeff_Bo_k = ((la[ell].s_wodm_k :/ a.N) - B_ol :* (la[ell].s_wod_k :/ a.N)) :/ d_minus
            a.links[ell].g_Bo_gamma = ((la[ell].s_wdmX / a.N) - B_ol * (la[ell].s_wdX / a.N)) / d_minus
            a.links[ell].g_dr_gamma = a.links[ell].g_C_gamma - a.links[ell].g_Bo_gamma
            if (missing(B_ol) | missing(a.links[ell].d_bar_minus_or) |
                missing(a.links[ell].a_coeff_Bo_k) | missing(a.links[ell].g_Bo_gamma))
                _cdid_aug_fail("Nonfinite augmented link control OR summary.", 504)
        }
        else {
            a.links[ell].d_bar_minus_or = J(1, 0, 0)
            a.links[ell].d_bar_diff_or = J(1, 0, 0)
            a.links[ell].num_Bo = .
            a.links[ell].B_ol = .
            a.links[ell].a_coeff_Bo_k = J(0, 1, .)
            a.links[ell].g_Bo_gamma = J(1, 0, .)
            a.links[ell].g_dr_gamma = J(1, 0, .)
        }
        a.links[ell].mu_ipw = (a.has_ps ? T - Cl : .)
        a.links[ell].mu_reg = (a.has_or ? T - B_gl : .)
        a.links[ell].mu_dr = (a.has_ps & a.has_or ? (T - B_gl) - (Cl - B_ol) : .)
        a.links[ell].mu = (a.method == "ipw" ? a.links[ell].mu_ipw :
            (a.method == "reg" ? a.links[ell].mu_reg : a.links[ell].mu_dr))
        if (missing(a.links[ell].mu))
            _cdid_aug_fail("Nonfinite augmented link point estimate.", 504)
    }
}

/* Global or block-local OR slope j forcing adjoints (one shared adjoint in both modes). */
void _cdid_aug_or_adjoint(struct cdid_aug_backend scalar a, struct cdid_imputation_model scalar m)
{
    real scalar po, tot_po, at, p, bid, k, cl, j
    real matrix lam, nu
    p = a.p_raw
    a.or_adj_ready = 0
    if (!a.has_or) return
    if (a.outcome_pooling == "separate") {
        tot_po = 0
        for (bid=1; bid<=a.n_blocks; bid++) {
            po = a.blocks[bid].or_p
            a.blocks[bid].or_adj_ready = (po > 0)
            tot_po = tot_po + po
        }
        if (tot_po == 0) return
        lam = J(a.n_cells, tot_po, 0)
        nu = J(a.n_cells, p * tot_po, 0)
        at = 0
        for (bid=1; bid<=a.n_blocks; bid++) {
            po = a.blocks[bid].or_p
            if (po == 0) continue
            for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                cl = a.blocks[bid].ctrl_cell[k]
                for (j=1; j<=po; j++) {
                    if (a.intercept) lam[cl, at + j] = lam[cl, at + j] + a.blocks[bid].jf_beta_alpha[j, k]
                    nu[cl, ((at + j - 1) * p + 1)..((at + j) * p)] = (nu[cl, ((at + j - 1) * p + 1)..((at + j) * p)] +
                        a.blocks[bid].jf_beta_beta[j, ((k - 1) * p + 1)..(k * p)])
                }
            }
            at = at + po
        }
        a.or_adj = cdid_imputation_adjoint_prepare(m, lam, nu)
        a.or_adj_ready = 1
        at = 0
        for (bid=1; bid<=a.n_blocks; bid++) {
            po = a.blocks[bid].or_p
            if (po == 0) continue
            a.blocks[bid].or_adj = cdid_imputation_adjoint_slice(m, a.or_adj, at + 1, at + po)
            at = at + po
        }
        return
    }
    po = a.or_p
    if (po == 0) return
    lam = J(a.n_cells, po, 0)
    nu = J(a.n_cells, p * po, 0)
    for (bid=1; bid<=a.n_blocks; bid++) {
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
            cl = a.blocks[bid].ctrl_cell[k]
            for (j=1; j<=po; j++) {
                if (a.intercept) lam[cl, j] = lam[cl, j] + a.blocks[bid].jf_beta_alpha[j, k]
                nu[cl, ((j - 1) * p + 1)..(j * p)] = (nu[cl, ((j - 1) * p + 1)..(j * p)] +
                    a.blocks[bid].jf_beta_beta[j, ((k - 1) * p + 1)..(k * p)])
            }
        }
    }
    a.or_adj = cdid_imputation_adjoint_prepare(m, lam, nu)
    a.or_adj_ready = 1
}

/* Local OR intercept blocks needed by a link batch (source alpha_bids, non-details). */
real colvector _cdid_aug_abids(struct cdid_aug_backend scalar a, real colvector lk)
{
    real colvector out
    real scalar i
    if (!a.has_or | !a.intercept) return(J(0, 1, 0))
    out = J(rows(lk), 1, 0)
    for (i=1; i<=rows(lk); i++) out[i] = a.links[lk[i]].or_block
    return(uniqrows(out))
}

/* Unique OR blocks consumed by a link batch (used when outcome_pooling == "separate"). */
real colvector _cdid_aug_bbids(struct cdid_aug_backend scalar a, real colvector lk)
{
    real colvector out
    real scalar i
    if (!a.has_or) return(J(0, 1, 0))
    out = J(rows(lk), 1, 0)
    for (i=1; i<=rows(lk); i++) out[i] = a.links[lk[i]].or_block
    return(uniqrows(out))
}

/* Reconstructed cells needed by one link batch (source augmented_if_needed_sets, non-details):
   requested link target cells, link control cells when has_ps, plus control cells of every
   contributing OR block. */
real colvector _cdid_aug_zids(struct cdid_aug_backend scalar a, real colvector lk, real colvector bids)
{
    real colvector z
    real scalar i
    z = J(0, 1, 0)
    for (i=1; i<=rows(lk); i++) {
        z = z \ a.links[lk[i]].target_cell
        if (a.has_ps) z = z \ a.links[lk[i]].ctrl_cell
    }
    for (i=1; i<=rows(bids); i++) z = z \ a.blocks[bids[i]].ctrl_cell
    return(uniqrows(z))
}

/* One bounded Stage 1 adjoint per link batch. Columns: the signed
   MD treated loads of each batch link (lam_T on its target cell, nu_T on its group slopes),
   then (when has_ps) the signed control loads (lam_C_k, nu_C_k), then the OR alpha j loads
   (jf_alpha_alpha, jf_alpha_beta) of the batch blocks. */
struct cdid_imputation_adjoint scalar _cdid_aug_batch_adj(struct cdid_aug_backend scalar a,
    struct cdid_imputation_model scalar m, real colvector lk, real colvector abids)
{
    real scalar n_l, so, M, p, mi, ell, g, gid, k, cl, aa, bid, mc, ma
    real rowvector kept
    real matrix lam, nu
    n_l = rows(lk)
    so = n_l * (1 + a.has_ps)
    p = a.p_raw
    M = so + rows(abids)
    lam = J(a.n_cells, M, 0)
    nu = J(a.n_cells, p * M, 0)
    for (mi=1; mi<=n_l; mi++) {
        ell = lk[mi]
        g = a.links[ell].target_cell
        gid = m.cells[g].group_id
        kept = m.groups[gid].kept
        lam[g, mi] = lam[g, mi] + a.links[ell].lam_T
        if (cols(kept) > 0)
            nu[g, (mi - 1) * p :+ kept] = nu[g, (mi - 1) * p :+ kept] + a.links[ell].nu_T
        if (a.has_ps) {
            mc = n_l + mi
            for (k=1; k<=a.links[ell].n_ctrl; k++) {
                cl = a.links[ell].ctrl_cell[k]
                lam[cl, mc] = lam[cl, mc] + a.links[ell].lam_C_k[k]
                if (p > 0) {
                    nu[cl, ((mc - 1) * p + 1)..(mc * p)] = (nu[cl, ((mc - 1) * p + 1)..(mc * p)] +
                        a.links[ell].nu_C_k[k, ])
                }
            }
        }
    }
    for (aa=1; aa<=rows(abids); aa++) {
        bid = abids[aa]
        ma = so + aa
        for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
            cl = a.blocks[bid].ctrl_cell[k]
            lam[cl, ma] = lam[cl, ma] + a.blocks[bid].jf_alpha_alpha[k]
            if (p > 0) {
                nu[cl, ((ma - 1) * p + 1)..(ma * p)] = (nu[cl, ((ma - 1) * p + 1)..(ma * p)] +
                    a.blocks[bid].jf_alpha_beta[k, ])
            }
        }
    }
    return(cdid_imputation_adjoint_prepare(m, lam, nu))
}

/* Append OR slope j forcing adjoint columns onto the link-batch adjoint so every Stage 1 j
   pullback for a chunk executes in one combined calendar-pair-ordered sweep. */
struct cdid_imputation_adjoint scalar _cdid_aug_cat_or_adj(struct cdid_aug_backend scalar a,
    struct cdid_imputation_model scalar m, struct cdid_imputation_adjoint scalar base_adj,
    real colvector bbids)
{
    struct cdid_imputation_adjoint scalar out
    real scalar sep, M0, tot_po, at_all, M_tot, t, bid, pob, gid, at, loc, cid
    real colvector act_mask, ix_act, ac, boff_all, cids_g
    real rowvector cc, cc_src
    real matrix act_keys
    if (!a.has_or) return(base_adj)
    sep = (a.outcome_pooling == "separate")
    tot_po = 0
    at_all = 0
    boff_all = J(0, 1, 0)
    if (sep) {
        boff_all = J(a.n_blocks, 1, 0)
        for (bid=1; bid<=a.n_blocks; bid++) {
            boff_all[bid] = at_all
            at_all = at_all + a.blocks[bid].or_p
        }
        for (t=1; t<=rows(bbids); t++) {
            bid = bbids[t]
            pob = a.blocks[bid].or_p
            if (pob > 0 & (!a.blocks[bid].or_adj_ready | !a.or_adj_ready))
                _cdid_aug_fail("Block-local OR j adjoint is not prepared.", 498)
            tot_po = tot_po + pob
        }
    }
    else {
        tot_po = a.or_p
        if (tot_po > 0 & !a.or_adj_ready)
            _cdid_aug_fail("Global common-slope j adjoint is not prepared.", 498)
    }
    if (tot_po == 0) return(base_adj)
    M0 = base_adj.n_contrasts
    M_tot = M0 + tot_po
    out = base_adj
    out.n_contrasts = M_tot
    out.tilde_lambda = J(m.n_cells, M_tot, 0)
    if (M0 > 0) out.tilde_lambda[, 1..M0] = base_adj.tilde_lambda
    for (gid=1; gid<=m.n_groups; gid++) {
        out.groups[gid].nu = J(m.groups[gid].p, M_tot, 0)
        if (m.groups[gid].p > 0 & M0 > 0) out.groups[gid].nu[, 1..M0] = base_adj.groups[gid].nu
    }
    act_mask = J(m.n_cells, 1, 0)
    if (rows(base_adj.active_cells) > 0) act_mask[base_adj.active_cells] = J(rows(base_adj.active_cells), 1, 1)
    cc = (M0 + 1)..M_tot
    if (sep & rows(bbids) < a.n_blocks) {
        cc_src = J(1, tot_po, 0)
        at = 0
        for (t=1; t<=rows(bbids); t++) {
            bid = bbids[t]
            pob = a.blocks[bid].or_p
            if (pob == 0) continue
            cc_src[(at + 1)..(at + pob)] = (boff_all[bid] + 1)..(boff_all[bid] + pob)
            at = at + pob
        }
        out.tilde_lambda[, cc] = a.or_adj.tilde_lambda[, cc_src]
        if (m.intercept) {
            for (cid=1; cid<=m.n_cells; cid++) {
                if (any(out.tilde_lambda[cid, cc] :!= 0)) act_mask[cid] = 1
            }
        }
        for (gid=1; gid<=m.n_groups; gid++) {
            if (m.groups[gid].p > 0) {
                out.groups[gid].nu[, cc] = a.or_adj.groups[gid].nu[, cc_src]
                if (any(out.groups[gid].nu[, cc] :!= 0)) {
                    cids_g = m.groups[gid].cell_ids
                    act_mask[cids_g] = J(rows(cids_g), 1, 1)
                }
            }
        }
    }
    else {
        out.tilde_lambda[, cc] = a.or_adj.tilde_lambda
        for (gid=1; gid<=m.n_groups; gid++) {
            if (m.groups[gid].p > 0) out.groups[gid].nu[, cc] = a.or_adj.groups[gid].nu
        }
        ac = a.or_adj.active_cells
        if (rows(ac) > 0) act_mask[ac] = J(rows(ac), 1, 1)
    }
    ix_act = _cdid_impb_ix(act_mask :== 1)
    if (rows(ix_act) > 1) {
        act_keys = J(rows(ix_act), 4, 0)
        for (loc=1; loc<=rows(ix_act); loc++) {
            cid = ix_act[loc]
            act_keys[loc, ] = (m.cells[cid].L, m.cells[cid].s_index, m.cells[cid].t_index, cid)
        }
        ix_act = ix_act[order(act_keys, (1, 2, 3, 4))]
    }
    out.active_cells = ix_act
    return(out)
}

/* Theta-only observation loads of one link batch in ONE selection moment record:
   G_theta carries sel_T (and sel_C when has_ps), q_alpha (per batch block) and q_beta
   (global common slope or per-batch-block separate slopes) at native theta IDs; G_ak stays ZERO. */
struct cdid_sel_moments scalar _cdid_aug_batch_sel(struct cdid_aug_backend scalar a,
    real colvector lk, real colvector abids, | real colvector bbids)
{
    struct cdid_sel_moments scalar sm
    real scalar n_l, so, nA, po, tot_po, sep, mi, ell, aa, t, bid, j, c0, off
    real colvector ids, bb_in
    n_l = rows(lk)
    so = n_l * (1 + a.has_ps)
    nA = rows(abids)
    sep = (a.has_or & a.outcome_pooling == "separate")
    bb_in = (args() >= 4 ? bbids : (sep ? _cdid_aug_bbids(a, lk) : J(0, 1, 0)))
    tot_po = 0
    if (a.has_or) {
        if (sep) {
            for (t=1; t<=rows(bb_in); t++) tot_po = tot_po + a.blocks[bb_in[t]].or_p
        }
        else tot_po = a.or_p
    }
    if (a.is_dated) {
        sm.n_moments = so + nA + tot_po
        sm.total_params = a.n_q_params
        sm.n_cohorts = a.G
        sm.G_theta = J(a.n_q_params, sm.n_moments, 0)
        sm.G_ak = J(a.G, sm.n_moments, 0)
    }
    else sm = cdid_sel_moments_init(a.qfit, so + nA + tot_po)
    for (mi=1; mi<=n_l; mi++) {
        ell = lk[mi]
        ids = a.links[ell].sel_T_ids
        if (rows(ids) > 0) sm.G_theta[ids, mi] = sm.G_theta[ids, mi] + a.links[ell].sel_T
        if (a.has_ps) {
            ids = a.links[ell].sel_C_ids
            if (rows(ids) > 0) sm.G_theta[ids, n_l + mi] = sm.G_theta[ids, n_l + mi] + a.links[ell].sel_C
        }
    }
    for (aa=1; aa<=nA; aa++) {
        bid = abids[aa]
        ids = a.blocks[bid].q_alpha_ids
        c0 = so + aa
        if (rows(ids) > 0) sm.G_theta[ids, c0] = sm.G_theta[ids, c0] + a.blocks[bid].q_alpha
    }
    if (a.has_or) {
        if (sep) {
            off = so + nA
            for (t=1; t<=rows(bb_in); t++) {
                bid = bb_in[t]
                po = a.blocks[bid].or_p
                ids = a.blocks[bid].q_beta_ids
                if (po > 0 & rows(ids) > 0) {
                    sm.G_theta[ids, (off + 1)..(off + po)] = (sm.G_theta[ids, (off + 1)..(off + po)] +
                        a.blocks[bid].q_beta')
                }
                off = off + po
            }
        }
        else if (a.or_p > 0 & rows(a.q_beta_ids) > 0) {
            sm.G_theta[a.q_beta_ids, (so + nA + 1)..(so + nA + a.or_p)] = (
                sm.G_theta[a.q_beta_ids, (so + nA + 1)..(so + nA + a.or_p)] + a.q_beta')
        }
    }
    if (missing(sm.G_theta) | missing(sm.G_ak))
        _cdid_aug_fail("Nonfinite theta-only observation loading record.", 504)
    if (any(sm.G_ak :!= 0))
        _cdid_aug_fail("Observation loading record must not carry a repeat-rate load.", 498)
    return(sm)
}

/* Common OR slope IF rows: N (R'R)^-1 rhs with the raw Z-response OR factor (N exactly once). */
real matrix _cdid_aug_or_solve(struct cdid_aug_backend scalar a, real matrix rhs)
{
    if (a.or_p == 0) return(J(rows(rhs), 0, 0))
    return(_cdid_nuisance_gram(a.or_R, a.or_pivot, a.N :* rhs')')
}

/* Complete original-unit IF of one bounded row chunk and one link batch at frame scale
   (source augmented_evaluate_if_chunk_core, non-details). Supports ipw, reg and dr with
   pooled or separate Z-response OR. Channels in source SCHEMA order. */
void _cdid_aug_eval(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    real colvector fr, real colvector lk, real colvector abids, real colvector bbids,
    struct cdid_imputation_adjoint scalar adj, struct cdid_sel_moments scalar sm,
    real scalar detail, real matrix o_phi, real matrix o_comp, real matrix o_ref,
    | real matrix o_or_alpha, real matrix o_or_beta)
{
    struct _cdid_impb_ivec rowvector c_ix
    real scalar B, nl, G, nr, N, p, q, po, pob, tot_po, sep, nA, r, k, mi, ell, t, tid, rtid, tp
    real scalar bid, apos, bpos, L, alpha_b, cl, ck, cg, g, aa, rt, rb, rl, xcol, cid, zs, nz, so
    real colvector coh, wr, Av, rep_tids, ix, v, mh, rv, ixg, ixk, bids, zids, zpos
    real colvector jh, zr, obs, qv, piv, ev, boff, obeta, center_weight, rad_b, raa_b, raj_b, raq_b
    real colvector wp, Zg, mg, dT, aT, jT, sT, dBg, aBg, oBd, oBa, oBj, oBq
    real colvector wm, Zk, mk, aC, aBo, dC, jC, sC, tCd, tCa, dBo, tBd, tBa, oOd, oOa, oOj, oOq
    real rowvector cc, okeep, bcols, allc
    real matrix xa, va, Z_rc, X_rc, wm_cache, Jb, S, XGd, XGa, sd, sa, X, vals, dor, dt, Zc, Tm, comp
    real matrix rbd, rba, rbj, rbq, rad, raa, raj, raq
    real matrix rbd_b, rba_b, rbj_b, rbq_b, xbd_b, xba_b, xbj_b, xbq_b
    real matrix xbd, xba, xbj, xbq, xad, xaa, xaj, xaq

    B = rows(fr)
    nl = rows(lk)
    G = a.G
    nr = a.n_recipes
    N = a.N
    p = a.p_raw
    q = a.q
    po = a.or_p
    sep = (a.has_or & a.outcome_pooling == "separate")
    nA = rows(abids)
    so = nl * (1 + a.has_ps)
    allc = _cdid_impb_seq(1, p)
    coh = m.unit_cohort_index[fr]
    c_ix = _cdid_impb_bucket_coh(coh, G)
    Av = m.pop.A[fr]

    tot_po = 0
    boff = J(rows(bbids), 1, 0)
    if (a.has_or) {
        if (sep) {
            for (t=1; t<=rows(bbids); t++) {
                boff[t] = tot_po
                tot_po = tot_po + a.blocks[bbids[t]].or_p
            }
        }
        else tot_po = po
    }

    /* 1. Repeat-rate IF xi_a and cached (W a_k, Z, X) slices for every used (recipe r, cohort k) on F. */
    xa = J(B, nr * G, 0)
    va = J(B, nr * G, 0)
    Z_rc = J(B, nr * p, 0)
    X_rc = (a.has_ps ? J(B, nr * q, 0) : J(B, 0, 0))
    for (r=1; r<=nr; r++) {
        for (k=1; k<=G; k++) {
            if (a.is_dated) {
                if (!(*(a.qdated)).used_cohorts[r, k]) continue
            }
            if (!(ib.mass_F_mat[r, k] > 0)) continue
            ix = _cdid_impb_cix(c_ix, m, fr, r, k)
            if (rows(ix) == 0) continue
            wr = _cdid_impb_w(ib, m, r, fr[ix], 1)
            xa[ix, (r - 1) * G + k] = wr :* (Av[ix] :- ib.a_k_mat[r, k]) :/ ib.mass_F_mat[r, k]
            va[ix, (r - 1) * G + k] = wr :* ib.a_k_mat[r, k]
            if (p > 0) Z_rc[ix, ((r - 1) * p + 1)..(r * p)] = _cdid_impb_z(ib, m, r, fr[ix], allc, 1)
            if (a.has_ps & q > 0) X_rc[ix, ((r - 1) * q + 1)..(r * q)] = _cdid_impb_x(ib, r, fr[ix], 1)
        }
    }

    /* 2. Stage 1 j pullbacks of the signed MD loads; theta-only observation contractions. */
    Jb = _cdid_impb_adj_chunk(ib, m, adj, fr)
    if (a.is_dated) S = _cdid_aug_q_loading(a, sm, fr)
    else S = cdid_sel_eval_loading(a.qfit, qr, sm, fr)
    if (rows(Jb) != B | cols(Jb) != so + nA + tot_po | rows(S) != B | cols(S) != so + nA + tot_po)
        _cdid_aug_fail("Batch adjoint or observation contraction has an unexpected shape.", 498)

    /* 3. Reconstructed Z of the needed cells on this chunk, once each (bounded cell cache). */
    if (!a.has_or) bids = J(0, 1, 0)
    else if (sep) bids = bbids
    else bids = (po > 0 ? (1::a.n_blocks) : abids)
    zids = _cdid_aug_order_cells(m, _cdid_aug_zids(a, lk, bids))
    nz = rows(zids)
    zpos = J(a.n_cells, 1, 0)
    Zc = J(B, nz, 0)
    for (zs=1; zs<=nz; zs++) {
        cid = zids[zs]
        zpos[cid] = zs
        ck = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        ix = _cdid_impb_cix(c_ix, m, fr, r, ck)
        if (rows(ix) == 0) continue
        jh = _cdid_impb_pred_z(m, cid, rows(ix),
            (p > 0 ? Z_rc[ix, ((r - 1) * p + 1)..(r * p)] : J(rows(ix), 0, 0)))
        _cdid_aug_cell_rows_z(a, m, qr, cid, fr[ix], 0, jh, zr=., obs=., qv=., piv=., ev=., Tm=.)
        Zc[ix, zs] = zr
    }

    /* 4. Treatment gamma IF from the borrowed roots when has_ps (deduplicated by rep_tid). */
    rep_tids = J(0, 1, 0)
    if (a.has_ps) {
        rep_tids = J(nl, 1, 0)
        for (mi=1; mi<=nl; mi++) rep_tids[mi] = ib.trts[a.links[lk[mi]].trt_id].rep_tid
        rep_tids = uniqrows(rep_tids)
    }
    XGd = J(B, q * rows(rep_tids), 0)
    XGa = J(B, q * rows(rep_tids), 0)
    wm_cache = J(B, rows(rep_tids), 0)
    for (t=1; t<=rows(rep_tids); t++) {
        rtid = rep_tids[t]
        rt = ib.trts[rtid].recipe
        cc = ((t - 1) * q + 1)..(t * q)
        sd = J(B, q, 0)
        for (k=1; k<=ib.trts[rtid].n_coh; k++) {
            ck = ib.trts[rtid].all_ci[k]
            ix = _cdid_impb_cix(c_ix, m, fr, rt, ck)
            if (rows(ix) == 0) continue
            v = va[ix, (rt - 1) * G + ck]
            X = (q > 0 ? X_rc[ix, ((rt - 1) * q + 1)..(rt * q)] : J(rows(ix), 0, 0))
            vals = _cdid_aug_logit_vals(X * ib.trts[rtid].gamma)
            if (ck == ib.trts[rtid].target_ci) {
                sd[ix, ] = (v :* (1 :- vals[, 1])) :* X
            }
            else {
                sd[ix, ] = (v :* (0 :- vals[, 1])) :* X
                wm_cache[ix, t] = v :* vals[, 2]
            }
        }
        sa = xa[, (rt - 1) * G :+ (ib.trts[rtid].all_ci')] * ib.trts[rtid].a_forcing'
        XGd[, cc] = _cdid_nuisance_gram(ib.trts[rtid].R, ib.trts[rtid].pivot, N :* sd')'
        XGa[, cc] = _cdid_nuisance_gram(ib.trts[rtid].R, ib.trts[rtid].pivot, N :* sa')'
    }

    /* 5. Z-response OR IF: direct, rate, j (MD adjoints) and q (theta-only) channels. */
    xbd = J(B, tot_po, 0)
    xba = J(B, tot_po, 0)
    xbj = J(B, tot_po, 0)
    xbq = J(B, tot_po, 0)
    xad = J(B, nA, 0)
    xaa = J(B, nA, 0)
    xaj = J(B, nA, 0)
    xaq = J(B, nA, 0)
    if (a.has_or) {
        if (sep) {
            for (t=1; t<=rows(bbids); t++) {
                bid = bbids[t]
                pob = a.blocks[bid].or_p
                okeep = a.blocks[bid].or_keep
                obeta = a.blocks[bid].or_beta
                rb = a.blocks[bid].recipe
                apos = _cdid_impb_pos(abids, bid)
                L = a.blocks[bid].L
                alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
                rad_b = J(B, 1, 0)
                raa_b = J(B, 1, 0)
                rbd_b = J(B, pob, 0)
                rba_b = J(B, pob, 0)
                for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                    cl = a.blocks[bid].ctrl_cell[k]
                    ck = a.blocks[bid].ctrl_ci[k]
                    xcol = (rb - 1) * G + ck
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) > 0) {
                        v = va[ix, xcol]
                        if (pob > 0) {
                            dor = L :* Z_rc[ix, (rb - 1) * p :+ okeep]
                            mh = alpha_b :+ dor * obeta
                        }
                        else mh = J(rows(ix), 1, alpha_b)
                        rv = v :* (Zc[ix, zpos[cl]] - mh)
                        /* Center the direct scores on this OR block's weighted support. */
                        center_weight = v :/ a.blocks[bid].mass
                        if (apos > 0)
                            rad_b[ix] = rad_b[ix] + (rv :- center_weight :* a.blocks[bid].score_mean_a)
                        if (pob > 0) {
                            dt = (a.intercept ? dor :- a.blocks[bid].mean_d : dor)
                            rbd_b[ix, ] = rbd_b[ix, ] + (rv :* dt - center_weight * a.blocks[bid].score_mean_b')
                        }
                    }
                    if (apos > 0) raa_b = raa_b + a.blocks[bid].a_force_alpha[k] :* xa[, xcol]
                    if (pob > 0) rba_b = rba_b + xa[, xcol] * a.blocks[bid].a_force_beta[, k]'
                }
                if (pob > 0) {
                    bcols = (boff[t] + 1)..(boff[t] + pob)
                    rbj_b = Jb[, (so + nA) :+ bcols]
                    rbq_b = S[, (so + nA) :+ bcols]
                    xbd_b = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, N :* rbd_b')'
                    xba_b = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, N :* rba_b')'
                    xbj_b = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, N :* rbj_b')'
                    xbq_b = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, N :* rbq_b')'
                    xbd[, bcols] = xbd_b
                    xba[, bcols] = xba_b
                    xbj[, bcols] = xbj_b
                    xbq[, bcols] = xbq_b
                }
                if (apos > 0) {
                    raj_b = Jb[, so + apos]
                    raq_b = S[, so + apos]
                    xad[, apos] = rad_b :/ a.blocks[bid].mass
                    xaa[, apos] = raa_b :/ a.blocks[bid].mass
                    xaj[, apos] = raj_b :/ a.blocks[bid].mass
                    xaq[, apos] = raq_b :/ a.blocks[bid].mass
                    if (pob > 0) {
                        xad[, apos] = xad[, apos] - xbd_b * a.blocks[bid].mean_d'
                        xaa[, apos] = xaa[, apos] - xba_b * a.blocks[bid].mean_d'
                        xaj[, apos] = xaj[, apos] - xbj_b * a.blocks[bid].mean_d'
                        xaq[, apos] = xaq[, apos] - xbq_b * a.blocks[bid].mean_d'
                    }
                }
            }
        }
        else {
            rbd = J(B, po, 0)
            rba = J(B, po, 0)
            rbj = J(B, po, 0)
            rbq = J(B, po, 0)
            rad = J(B, nA, 0)
            raa = J(B, nA, 0)
            raj = J(B, nA, 0)
            raq = J(B, nA, 0)
            for (t=1; t<=rows(bids); t++) {
                bid = bids[t]
                rb = a.blocks[bid].recipe
                apos = _cdid_impb_pos(abids, bid)
                L = a.blocks[bid].L
                alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
                for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                    cl = a.blocks[bid].ctrl_cell[k]
                    ck = a.blocks[bid].ctrl_ci[k]
                    xcol = (rb - 1) * G + ck
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) > 0) {
                        v = va[ix, xcol]
                        if (po > 0) {
                            dor = L :* Z_rc[ix, (rb - 1) * p :+ a.or_keep]
                            mh = alpha_b :+ dor * a.or_beta
                        }
                        else mh = J(rows(ix), 1, alpha_b)
                        rv = v :* (Zc[ix, zpos[cl]] - mh)
                        if (apos > 0) rad[ix, apos] = rad[ix, apos] + rv
                        if (po > 0) {
                            dt = (a.intercept ? dor :- a.blocks[bid].mean_d : dor)
                            rbd[ix, ] = rbd[ix, ] + rv :* dt
                        }
                    }
                    if (apos > 0) raa[, apos] = raa[, apos] + a.blocks[bid].a_force_alpha[k] :* xa[, xcol]
                    if (po > 0) rba = rba + xa[, xcol] * a.blocks[bid].a_force_beta[, k]'
                }
            }
            if (po > 0) {
                rbj = Jb[, (so + nA + 1)..(so + nA + po)]
                rbq = S[, (so + nA + 1)..(so + nA + po)]
                xbd = _cdid_aug_or_solve(a, rbd)
                xba = _cdid_aug_or_solve(a, rba)
                xbj = _cdid_aug_or_solve(a, rbj)
                xbq = _cdid_aug_or_solve(a, rbq)
            }
            for (aa=1; aa<=nA; aa++) {
                bid = abids[aa]
                raj[, aa] = Jb[, so + aa]
                raq[, aa] = S[, so + aa]
                xad[, aa] = rad[, aa] :/ a.blocks[bid].mass
                xaa[, aa] = raa[, aa] :/ a.blocks[bid].mass
                xaj[, aa] = raj[, aa] :/ a.blocks[bid].mass
                xaq[, aa] = raq[, aa] :/ a.blocks[bid].mass
                if (po > 0) {
                    xad[, aa] = xad[, aa] - xbd * a.blocks[bid].mean_d'
                    xaa[, aa] = xaa[, aa] - xba * a.blocks[bid].mean_d'
                    xaj[, aa] = xaj[, aa] - xbj * a.blocks[bid].mean_d'
                    xaq[, aa] = xaq[, aa] - xbq * a.blocks[bid].mean_d'
                }
            }
        }
    }

    /* Optional: total coefficient IF of the Z-response OR (opt-in reports only). */
    if (args() >= 16) {
        o_or_alpha = xad + xaa + xaj + xaq
        o_or_beta = xbd + xba + xbj + xbq
    }

    /* 6. Link rows: ten source channels, summed on each original unit before any product. */
    o_phi = J(B, nl, 0)
    o_comp = J(0, 0, .)
    o_ref = J(0, 0, .)
    if (detail) {
        o_comp = J(B, 10 * nl, 0)
        o_ref = J(B, 4 * nl, 0)
    }
    pob = 0
    okeep = J(1, 0, 0)
    obeta = J(0, 1, 0)
    bcols = J(1, 0, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        rl = a.links[ell].recipe
        g = a.links[ell].target_cell
        cg = a.links[ell].target_ci
        L = a.links[ell].L
        apos = 0
        alpha_b = 0
        ixg = _cdid_impb_cix(c_ix, m, fr, rl, cg)
        wp = J(B, 1, 0)
        Zg = J(B, 1, 0)
        if (rows(ixg) > 0) {
            wp[ixg] = va[ixg, (rl - 1) * G + cg]
            Zg[ixg] = Zc[ixg, zpos[g]]
        }
        dT = wp :* (Zg :- a.links[ell].T_gl) :/ a.links[ell].d_plus
        aT = a.links[ell].a_coeff_T_g :* xa[, (rl - 1) * G + cg]
        jT = Jb[, mi]
        sT = S[, mi]

        if (a.has_or) {
            bid = a.links[ell].or_block
            apos = _cdid_impb_pos(abids, bid)
            pob = (sep ? a.blocks[bid].or_p : po)
            okeep = (sep ? a.blocks[bid].or_keep : a.or_keep)
            obeta = (sep ? a.blocks[bid].or_beta : a.or_beta)
            if (pob > 0) {
                if (sep) {
                    bpos = _cdid_impb_pos(bbids, bid)
                    bcols = (boff[bpos] + 1)..(boff[bpos] + pob)
                }
                else bcols = 1..pob
            }
            else bcols = J(1, 0, 0)
            alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
            mg = J(B, 1, 0)
            if (rows(ixg) > 0) {
                if (pob > 0) mg[ixg] = alpha_b :+ (L :* Z_rc[ixg, (rl - 1) * p :+ okeep]) * obeta
                else mg[ixg] = J(rows(ixg), 1, alpha_b)
            }
            dBg = wp :* (mg :- a.links[ell].B_gl) :/ a.links[ell].d_plus
            aBg = a.links[ell].a_coeff_Bg_g :* xa[, (rl - 1) * G + cg]
            oBd = J(B, 1, 0)
            oBa = J(B, 1, 0)
            oBj = J(B, 1, 0)
            oBq = J(B, 1, 0)
            if (apos > 0) {
                oBd = oBd + xad[, apos]
                oBa = oBa + xaa[, apos]
                oBj = oBj + xaj[, apos]
                oBq = oBq + xaq[, apos]
            }
            if (pob > 0) {
                oBd = oBd + xbd[, bcols] * a.links[ell].d_bar_plus_or'
                oBa = oBa + xba[, bcols] * a.links[ell].d_bar_plus_or'
                oBj = oBj + xbj[, bcols] * a.links[ell].d_bar_plus_or'
                oBq = oBq + xbq[, bcols] * a.links[ell].d_bar_plus_or'
            }
        }

        if (a.has_ps) {
            tid = a.links[ell].trt_id
            tp = _cdid_impb_pos(rep_tids, ib.trts[tid].rep_tid)
            cc = ((tp - 1) * q + 1)..(tp * q)
            wm = wm_cache[, tp]
            Zk = J(B, 1, 0)
            mk = J(B, 1, 0)
            aC = J(B, 1, 0)
            aBo = J(B, 1, 0)
            for (k=1; k<=a.links[ell].n_ctrl; k++) {
                ck = a.links[ell].ctrl_ci[k]
                cl = a.links[ell].ctrl_cell[k]
                xcol = (rl - 1) * G + ck
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) > 0) {
                    Zk[ixk] = Zc[ixk, zpos[cl]]
                    if (a.has_or) {
                        if (pob > 0) mk[ixk] = alpha_b :+ (L :* Z_rc[ixk, (rl - 1) * p :+ okeep]) * obeta
                        else mk[ixk] = J(rows(ixk), 1, alpha_b)
                    }
                }
                aC = aC + a.links[ell].a_coeff_C_k[k] :* xa[, xcol]
                if (a.has_or) aBo = aBo + a.links[ell].a_coeff_Bo_k[k] :* xa[, xcol]
            }
            dC = wm :* (Zk :- a.links[ell].C_ol) :/ a.links[ell].d_minus
            jC = Jb[, nl + mi]
            sC = S[, nl + mi]
            tCd = XGd[, cc] * (a.method == "dr" ? a.links[ell].g_dr_gamma' : a.links[ell].g_C_gamma')
            tCa = XGa[, cc] * (a.method == "dr" ? a.links[ell].g_dr_gamma' : a.links[ell].g_C_gamma')
            if (a.has_or) {
                dBo = wm :* (mk :- a.links[ell].B_ol) :/ a.links[ell].d_minus
                tBd = J(B, 1, 0)
                tBa = J(B, 1, 0)
                oOd = J(B, 1, 0)
                oOa = J(B, 1, 0)
                oOj = J(B, 1, 0)
                oOq = J(B, 1, 0)
                if (pob > 0) {
                    oBd = xbd[, bcols] * a.links[ell].d_bar_diff_or'
                    oBa = xba[, bcols] * a.links[ell].d_bar_diff_or'
                    oBj = xbj[, bcols] * a.links[ell].d_bar_diff_or'
                    oBq = xbq[, bcols] * a.links[ell].d_bar_diff_or'
                }
                else {
                    oBd = J(B, 1, 0)
                    oBa = J(B, 1, 0)
                    oBj = J(B, 1, 0)
                    oBq = J(B, 1, 0)
                }
            }
        }

        if (a.method == "ipw") {
            comp = (dT - dC, aT - aC, jT - jC, sT - sC, -tCd, -tCa, J(B, 4, 0))
        }
        else if (a.method == "reg") {
            comp = (dT - dBg, aT - aBg, jT, sT, J(B, 2, 0), -oBd, -oBa, -oBj, -oBq)
        }
        else {
            comp = ((dT - dBg) - (dC - dBo), (aT - aBg) - (aC - aBo), jT - jC, sT - sC,
                -(tCd - tBd), -(tCa - tBa), -oBd + oOd, -oBa + oOa, -oBj + oOj, -oBq + oOq)
        }
        o_phi[, mi] = quadrowsum(comp)
        if (detail) {
            o_comp[, ((mi - 1) * 10 + 1)..(mi * 10)] = comp
            o_ref[, ((mi - 1) * 4 + 1)..(mi * 4)] = (comp[, 3] + comp[, 9], comp[, 4] + comp[, 10],
                comp[, 2] + comp[, 8] + comp[, 6], comp[, 1] + comp[, 7] + comp[, 5])
        }
    }
    if (missing(o_phi)) _cdid_aug_fail("Nonfinite complete augmented link IF.", 504)
}

void _cdid_aug_set_ids(struct cdid_aug_backend_if scalar res,
    struct cdid_imputation_model scalar m, real colvector fr)
{
    if (m.pop.ids.string_id) {
        res.id_str = m.pop.ids.string_id_value[fr]
        res.id_num = J(0, 1, .)
    }
    else {
        res.id_num = m.pop.ids.numeric_id[fr]
        res.id_str = J(0, 1, "")
    }
}

/* Replay scratch (bounded chunk, link batch, Z cell cache, adjoint and theta contraction)
   plus requested output, charged with the concurrently live retained stores; the opt-in
   work cap applies to each replay call. */
void _cdid_aug_replay_budget(struct cdid_aug_backend scalar a, real scalar Bc, real scalar nl,
    real scalar B, real scalar detail)
{
    real scalar nlb, nbb, M_adj, M_sel, scr, work, xa_cols, p, bid, max_po
    p = a.p_raw
    nlb = min((nl, a.link_batch))
    nbb = min((nlb, a.n_blocks))
    max_po = a.or_p
    if (a.has_or & a.outcome_pooling == "separate") {
        max_po = 0
        for (bid=1; bid<=a.n_blocks; bid++) max_po = max((max_po, a.blocks[bid].or_p))
        max_po = nbb * max_po
    }
    M_adj = nlb * (1 + a.has_ps) + (a.has_or ? nbb : 0)
    M_sel = nlb * (1 + a.has_ps) + (a.has_or ? nbb + max_po : 0)
    xa_cols = a.n_recipes * a.G
    scr = (8 * (Bc * (xa_cols + 4 * a.q * nlb + 13 * p + 6 * nbb + M_adj + M_sel + a.n_cells + 64) +
        a.n_cells * M_adj * (1 + p) + a.n_groups * p * M_adj + a.n_q_params * M_sel +
        B * nl * (1 + 14 * detail)) + 4096)
    work = (16 * B * (nl * (xa_cols + 3) * (p + a.q + 4) + a.n_cells * (M_adj + 1) * (p + 1) +
        (a.n_blocks + nl) * (a.d_sel + 2) * (M_sel + 1)) * ceil(nl / nlb))
    _cdid_imp_check_budget(a.budget, a.retained_bytes, scr, work, "augmented backend replay")
}

/* True when evaluating all nl links in one detail=0 pass fits within max_workspace_bytes. */
real scalar _cdid_aug_can_full_batch(struct cdid_aug_backend scalar a, real scalar Bc,
    real scalar nl, real scalar B)
{
    real scalar nbb, M_adj, M_sel, scr, xa_cols, p, bid, max_po
    if (nl <= a.link_batch) return(1)
    p = a.p_raw
    nbb = min((nl, a.n_blocks))
    max_po = a.or_p
    if (a.has_or & a.outcome_pooling == "separate") {
        max_po = 0
        for (bid=1; bid<=a.n_blocks; bid++) max_po = max((max_po, a.blocks[bid].or_p))
        max_po = nbb * max_po
    }
    M_adj = nl * (1 + a.has_ps) + (a.has_or ? nbb : 0)
    M_sel = nl * (1 + a.has_ps) + (a.has_or ? nbb + max_po : 0)
    xa_cols = a.n_recipes * a.G
    scr = (8 * (Bc * (xa_cols + 4 * a.q * nl + 13 * p + 6 * nbb + M_adj + M_sel + a.n_cells + 64) +
        a.n_cells * M_adj * (1 + p) + a.n_groups * p * M_adj + a.n_q_params * M_sel +
        B * nl) + 4096)
    return(a.retained_bytes + scr <= a.budget.max_workspace_bytes)
}

/* Pre-build augmented IF (detail = 0) adjoint, theta loading, cell maps, and BLAS loading matrices. */
struct _cdid_aug_plan scalar _cdid_aug_make_plan(struct cdid_aug_backend scalar a,
    struct cdid_imp_backend scalar ib, struct cdid_imputation_model scalar m, real colvector lk)
{
    struct _cdid_aug_plan scalar plan
    struct cdid_imputation_adjoint scalar base_adj
    real scalar nl, G, nr, p, q, po, pob, tot_po, mi, ell, t, tb, bid, bpos, apos, k, ck, cg, rl, rb, xcol, n_rep, s
    real scalar c_tgt, c_ctl
    real colvector bids_z, or_mark, pos_v
    real rowvector bcols

    nl = rows(lk)
    G = a.G
    nr = a.n_recipes
    p = a.p_raw
    q = a.q
    po = a.or_p
    plan.lk = lk
    plan.nl = nl
    plan.so = nl * (1 + a.has_ps)
    plan.sep = (a.has_or & a.outcome_pooling == "separate")
    plan.abids = (a.has_or & a.intercept & a.method == "reg" ? _cdid_aug_abids(a, lk) : J(0, 1, 0))
    plan.nA = rows(plan.abids)
    plan.bbids = (plan.sep ? _cdid_aug_bbids(a, lk) : J(0, 1, 0))
    plan.bids = J(0, 1, 0)
    plan.boff = J(0, 1, 0)
    tot_po = 0
    if (a.has_or) {
        if (plan.sep) {
            plan.boff = J(rows(plan.bbids), 1, 0)
            for (t=1; t<=rows(plan.bbids); t++) {
                plan.boff[t] = tot_po
                tot_po = tot_po + a.blocks[plan.bbids[t]].or_p
            }
            bids_z = plan.bbids
        }
        else {
            plan.bids = (po > 0 ? (1::a.n_blocks) : plan.abids)
            tot_po = po
            bids_z = plan.bids
        }
    }
    else bids_z = J(0, 1, 0)
    plan.tot_po = tot_po

    base_adj = _cdid_aug_batch_adj(a, m, lk, plan.abids)
    plan.adj = _cdid_aug_cat_or_adj(a, m, base_adj, plan.bbids)
    plan.sm = _cdid_aug_batch_sel(a, lk, plan.abids, plan.bbids)

    plan.zids = _cdid_aug_order_cells(m, _cdid_aug_zids(a, lk, bids_z))
    plan.zpos = J(a.n_cells, 1, 0)
    for (t=1; t<=rows(plan.zids); t++) plan.zpos[plan.zids[t]] = t

    plan.rep_tids = J(0, 1, 0)
    plan.link_tp = J(nl, 1, 0)
    plan.tp_grps = _cdid_aug_tp_grp(1, 0)
    if (a.has_ps) {
        plan.rep_tids = J(nl, 1, 0)
        for (mi=1; mi<=nl; mi++) plan.rep_tids[mi] = ib.trts[a.links[lk[mi]].trt_id].rep_tid
        plan.rep_tids = uniqrows(plan.rep_tids)
        n_rep = rows(plan.rep_tids)
        for (mi=1; mi<=nl; mi++) {
            plan.link_tp[mi] = _cdid_impb_pos(plan.rep_tids, ib.trts[a.links[lk[mi]].trt_id].rep_tid)
        }
        plan.tp_grps = _cdid_aug_tp_grp(1, n_rep)
        for (t=1; t<=n_rep; t++) {
            pos_v = _cdid_impb_ix(plan.link_tp :== t)
            plan.tp_grps[t].rtid = plan.rep_tids[t]
            plan.tp_grps[t].link_pos = pos_v'
            plan.tp_grps[t].G_gamma = J(q, rows(pos_v), 0)
            for (s=1; s<=rows(pos_v); s++) {
                ell = lk[pos_v[s]]
                plan.tp_grps[t].G_gamma[, s] = (a.method == "ipw" ?
                    a.links[ell].g_C_gamma' : a.links[ell].g_dr_gamma')
            }
        }
    }

    plan.A_load = J(nr * G, nl, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        rl = a.links[ell].recipe
        cg = a.links[ell].target_ci
        c_tgt = a.links[ell].a_coeff_T_g - (a.has_or ? a.links[ell].a_coeff_Bg_g : 0)
        plan.A_load[(rl - 1) * G + cg, mi] = plan.A_load[(rl - 1) * G + cg, mi] + c_tgt
        if (a.has_ps) {
            for (k=1; k<=a.links[ell].n_ctrl; k++) {
                ck = a.links[ell].ctrl_ci[k]
                c_ctl = a.links[ell].a_coeff_C_k[k] - (a.has_or ? a.links[ell].a_coeff_Bo_k[k] : 0)
                plan.A_load[(rl - 1) * G + ck, mi] = plan.A_load[(rl - 1) * G + ck, mi] - c_ctl
            }
        }
    }

    plan.or_bids = J(0, 1, 0)
    plan.or_bpos = J(a.n_blocks, 1, 0)
    plan.or_coh_mask = J(0, G, 0)
    plan.Alpha_load = J(plan.nA, nl, 0)
    plan.Beta_load = J(tot_po, nl, 0)
    plan.Raa_load = J(nr * G, plan.nA, 0)
    plan.Rba_load = J(nr * G, po, 0)
    if (a.has_or) {
        or_mark = J(a.n_blocks, 1, 0)
        for (mi=1; mi<=nl; mi++) or_mark[a.links[lk[mi]].or_block] = 1
        if (!plan.sep) {
            for (t=1; t<=rows(plan.bids); t++) or_mark[plan.bids[t]] = 1
        }
        plan.or_bids = _cdid_impb_ix(or_mark :== 1)
        for (t=1; t<=rows(plan.or_bids); t++) plan.or_bpos[plan.or_bids[t]] = t
        plan.or_coh_mask = J(rows(plan.or_bids), G, 0)
        if (plan.sep) {
            for (t=1; t<=rows(plan.bbids); t++) {
                bid = plan.bbids[t]
                tb = plan.or_bpos[bid]
                if (a.blocks[bid].n_ctrl > 0)
                    plan.or_coh_mask[tb, a.blocks[bid].ctrl_ci'] = J(1, a.blocks[bid].n_ctrl, 1)
            }
        }
        else {
            for (t=1; t<=rows(plan.bids); t++) {
                bid = plan.bids[t]
                tb = plan.or_bpos[bid]
                if (a.blocks[bid].n_ctrl > 0)
                    plan.or_coh_mask[tb, a.blocks[bid].ctrl_ci'] = J(1, a.blocks[bid].n_ctrl, 1)
                rb = a.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                    xcol = (rb - 1) * G + a.blocks[bid].ctrl_ci[k]
                    if (apos > 0) plan.Raa_load[xcol, apos] = plan.Raa_load[xcol, apos] + a.blocks[bid].a_force_alpha[k]
                    if (po > 0) plan.Rba_load[xcol, ] = plan.Rba_load[xcol, ] + a.blocks[bid].a_force_beta[, k]'
                }
            }
        }
        for (mi=1; mi<=nl; mi++) {
            ell = lk[mi]
            bid = a.links[ell].or_block
            tb = plan.or_bpos[bid]
            plan.or_coh_mask[tb, a.links[ell].target_ci] = 1
            if (a.has_ps & a.links[ell].n_ctrl > 0)
                plan.or_coh_mask[tb, a.links[ell].ctrl_ci'] = J(1, a.links[ell].n_ctrl, 1)
            if (plan.nA > 0) {
                apos = _cdid_impb_pos(plan.abids, bid)
                if (apos > 0) plan.Alpha_load[apos, mi] = -1
            }
            pob = (plan.sep ? a.blocks[bid].or_p : po)
            if (pob > 0) {
                if (plan.sep) {
                    bpos = _cdid_impb_pos(plan.bbids, bid)
                    bcols = (plan.boff[bpos] + 1)..(plan.boff[bpos] + pob)
                }
                else bcols = 1..pob
                plan.Beta_load[bcols, mi] = -(a.method == "reg" ?
                    a.links[ell].d_bar_plus_or' : a.links[ell].d_bar_diff_or')
            }
        }
    }
    plan.ready = 1
    return(plan)
}

/* Batched augmented IF kernel for detail = 0 (frame scale). Groups pair reads across cohorts,
   combines j+q and direct+rate solves, and contracts loadings via BLAS matrix products. */
void _cdid_aug_eval_fast(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, real colvector fr,
    struct _cdid_aug_plan scalar plan, real matrix o_phi)
{
    struct _cdid_impb_ivec rowvector c_ix
    struct cdid_sel_pred scalar pr
    real scalar B, nl, G, nr, N, p, q, po, pob, tot_po, nA, so, r, k, mi, ell, t, tb, rtid, tp
    real scalar bid, apos, L, alpha_b, cl, ck, cg, g, aa, rt, rb, rl, cid, zs, nz, n_rep
    real scalar g0, g1, s_idx, t_idx, has_any, pc
    real colvector coh, wr, Av, ix, v, rv, ixg, ixk, jh, obs, rch_obs, piv, ev, obeta
    real colvector center_weight, rad_b, raa_b, rajq_b
    real rowvector cc, okeep, bcols, allc, lpos
    real matrix xa, va, Z_rc, X_rc, wm_cache, Jb, S, JS, XG_tot, sd, sa, X, vals, dor, dt
    real matrix Zc, Mc, hit_all, rbd, rba, rad, raa, rajq, rbd_b, rba_b
    real matrix xb_tot, xa_tot, xb_b, sol

    B = rows(fr)
    nl = plan.nl
    G = a.G
    nr = a.n_recipes
    N = a.N
    p = a.p_raw
    q = a.q
    po = a.or_p
    nA = plan.nA
    so = plan.so
    tot_po = plan.tot_po
    allc = _cdid_impb_seq(1, p)
    coh = m.unit_cohort_index[fr]
    c_ix = _cdid_impb_bucket_coh(coh, G)
    Av = m.pop.A[fr]

    xa = J(B, nr * G, 0)
    va = J(B, nr * G, 0)
    Z_rc = J(B, nr * p, 0)
    X_rc = (a.has_ps ? J(B, nr * q, 0) : J(B, 0, 0))
    for (r=1; r<=nr; r++) {
        for (k=1; k<=G; k++) {
            if (a.is_dated) {
                if (!(*(a.qdated)).used_cohorts[r, k]) continue
            }
            if (!(ib.mass_F_mat[r, k] > 0)) continue
            ix = _cdid_impb_cix(c_ix, m, fr, r, k)
            if (rows(ix) == 0) continue
            wr = _cdid_impb_w(ib, m, r, fr[ix], 1)
            xa[ix, (r - 1) * G + k] = wr :* (Av[ix] :- ib.a_k_mat[r, k]) :/ ib.mass_F_mat[r, k]
            va[ix, (r - 1) * G + k] = wr :* ib.a_k_mat[r, k]
            if (p > 0) Z_rc[ix, ((r - 1) * p + 1)..(r * p)] = _cdid_impb_z(ib, m, r, fr[ix], allc, 1)
            if (a.has_ps & q > 0) X_rc[ix, ((r - 1) * q + 1)..(r * q)] = _cdid_impb_x(ib, r, fr[ix], 1)
        }
    }

    Jb = _cdid_impb_adj_chunk(ib, m, plan.adj, fr)
    if (a.is_dated) S = _cdid_aug_q_loading(a, plan.sm, fr)
    else S = cdid_sel_eval_loading(a.qfit, *(a.qrdr), plan.sm, fr)
    if (rows(Jb) != B | cols(Jb) != so + nA + tot_po | rows(S) != B | cols(S) != so + nA + tot_po)
        _cdid_aug_fail("Batch adjoint or observation contraction has an unexpected shape.", 498)
    JS = Jb + S

    if (a.has_ps) o_phi = (JS[, 1..nl] - JS[, (nl + 1)..(2 * nl)]) + xa * plan.A_load
    else o_phi = JS[, 1..nl] + xa * plan.A_load

    nz = rows(plan.zids)
    Zc = J(B, nz, 0)
    g0 = 1
    while (g0 <= nz) {
        s_idx = m.cells[plan.zids[g0]].s_index
        t_idx = m.cells[plan.zids[g0]].t_index
        g1 = g0
        while (g1 < nz) {
            if (m.cells[plan.zids[g1 + 1]].s_index != s_idx |
                m.cells[plan.zids[g1 + 1]].t_index != t_idx) break
            g1++
        }
        has_any = 0
        for (zs=g0; zs<=g1; zs++) {
            if (rows(c_ix[m.cells[plan.zids[zs]].cohort_index].rows) > 0) {
                has_any = 1
                break
            }
        }
        if (has_any) {
            hit_all = cdid_imputation_read_pair_chunk(m.reader, m.pop,
                m.pairs[m.cells[plan.zids[g0]].pair_id], fr, 1)
            for (zs=g0; zs<=g1; zs++) {
                cid = plan.zids[zs]
                ck = m.cells[cid].cohort_index
                ix = _cdid_impb_cix(c_ix, m, fr, m.cells[cid].recipe, ck)
                if (rows(ix) == 0) continue
                r = m.cells[cid].recipe
                jh = _cdid_impb_pred_z(m, cid, rows(ix),
                    (p > 0 ? Z_rc[ix, ((r - 1) * p + 1)..(r * p)] : J(rows(ix), 0, 0)))
                obs = _cdid_impb_ix(hit_all[ix, 1] :== 1)
                if (rows(obs) > 0) {
                    pc = a.q_of_cell[cid]
                    rch_obs = (fr[ix])[obs]
                    if (a.is_dated) {
                        _cdid_aug_q_activate(a, r)
                        pr = cdid_sel_pred_pair(a.qcontexts[r].fit, a.qlease, pc, rch_obs)
                    }
                    else {
                        pr = cdid_sel_pred_pair(a.qfit, *(a.qrdr), pc, rch_obs)
                    }
                    if (rows(pr.p_pair) != rows(obs) | rows(pr.p_weight) != rows(obs) |
                        missing(pr.p_pair) | missing(pr.p_weight))
                        _cdid_aug_fail("Nonfinite fitted observation probabilities on observed pairs.", 504)
                    piv = 1 :/ pr.p_weight
                    ev = hit_all[ix[obs], 2] - jh[obs]
                    jh[obs] = jh[obs] + piv :* ev
                }
                Zc[ix, zs] = jh
            }
        }
        g0 = g1 + 1
    }

    Mc = J(B, rows(plan.or_bids), 0)
    if (a.has_or) {
        for (t=1; t<=rows(plan.or_bids); t++) {
            bid = plan.or_bids[t]
            rb = a.blocks[bid].recipe
            L = a.blocks[bid].L
            alpha_b = (a.intercept ? a.blocks[bid].alpha : 0)
            pob = (plan.sep ? a.blocks[bid].or_p : po)
            okeep = (plan.sep ? a.blocks[bid].or_keep : a.or_keep)
            obeta = (plan.sep ? a.blocks[bid].or_beta : a.or_beta)
            for (ck=1; ck<=G; ck++) {
                if (!plan.or_coh_mask[t, ck]) continue
                ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                if (rows(ix) == 0) continue
                if (pob > 0) Mc[ix, t] = alpha_b :+ (L :* Z_rc[ix, (rb - 1) * p :+ okeep]) * obeta
                else Mc[ix, t] = J(rows(ix), 1, alpha_b)
            }
        }
    }

    n_rep = rows(plan.rep_tids)
    wm_cache = J(B, n_rep, 0)
    if (a.has_ps) {
        XG_tot = J(B, q * n_rep, 0)
        for (t=1; t<=n_rep; t++) {
            rtid = plan.rep_tids[t]
            rt = ib.trts[rtid].recipe
            cc = ((t - 1) * q + 1)..(t * q)
            sd = J(B, q, 0)
            for (k=1; k<=ib.trts[rtid].n_coh; k++) {
                ck = ib.trts[rtid].all_ci[k]
                ix = _cdid_impb_cix(c_ix, m, fr, rt, ck)
                if (rows(ix) == 0) continue
                v = va[ix, (rt - 1) * G + ck]
                X = (q > 0 ? X_rc[ix, ((rt - 1) * q + 1)..(rt * q)] : J(rows(ix), 0, 0))
                vals = _cdid_aug_logit_vals(X * ib.trts[rtid].gamma)
                if (ck == ib.trts[rtid].target_ci) {
                    sd[ix, ] = (v :* (1 :- vals[, 1])) :* X
                }
                else {
                    sd[ix, ] = (v :* (0 :- vals[, 1])) :* X
                    wm_cache[ix, t] = v :* vals[, 2]
                }
            }
            sa = xa[, (rt - 1) * G :+ (ib.trts[rtid].all_ci')] * ib.trts[rtid].a_forcing'
            sol = _cdid_nuisance_gram(ib.trts[rtid].R, ib.trts[rtid].pivot, N :* (sd \ sa)')'
            XG_tot[, cc] = sol[1..B, ] + sol[(B + 1)..(2 * B), ]
            lpos = plan.tp_grps[t].link_pos
            if (cols(lpos) > 0) {
                o_phi[, lpos] = o_phi[, lpos] - XG_tot[, cc] * plan.tp_grps[t].G_gamma
            }
        }
    }

    if (a.has_or) {
        xa_tot = J(B, nA, 0)
        xb_tot = J(B, tot_po, 0)
        if (plan.sep) {
            for (t=1; t<=rows(plan.bbids); t++) {
                bid = plan.bbids[t]
                tb = plan.or_bpos[bid]
                pob = a.blocks[bid].or_p
                okeep = a.blocks[bid].or_keep
                rb = a.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                L = a.blocks[bid].L
                rad_b = J(B, 1, 0)
                rbd_b = J(B, pob, 0)
                for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                    cl = a.blocks[bid].ctrl_cell[k]
                    ck = a.blocks[bid].ctrl_ci[k]
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) == 0) continue
                    v = va[ix, (rb - 1) * G + ck]
                    rv = v :* (Zc[ix, plan.zpos[cl]] - Mc[ix, tb])
                    /* Center the direct scores on this OR block's weighted support. */
                    center_weight = v :/ a.blocks[bid].mass
                    if (apos > 0)
                        rad_b[ix] = rad_b[ix] + (rv :- center_weight :* a.blocks[bid].score_mean_a)
                    if (pob > 0) {
                        dor = L :* Z_rc[ix, (rb - 1) * p :+ okeep]
                        dt = (a.intercept ? dor :- a.blocks[bid].mean_d : dor)
                        rbd_b[ix, ] = rbd_b[ix, ] + (rv :* dt - center_weight * a.blocks[bid].score_mean_b')
                    }
                }
                raa_b = (apos > 0 ?
                    xa[, (rb - 1) * G :+ (a.blocks[bid].ctrl_ci')] * a.blocks[bid].a_force_alpha : J(B, 1, 0))
                rba_b = (pob > 0 ?
                    xa[, (rb - 1) * G :+ (a.blocks[bid].ctrl_ci')] * a.blocks[bid].a_force_beta' : J(B, 0, 0))
                xb_b = J(B, 0, 0)
                if (pob > 0) {
                    bcols = (plan.boff[t] + 1)..(plan.boff[t] + pob)
                    sol = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot,
                        N :* (rbd_b \ rba_b \ Jb[, (so + nA) :+ bcols] \ S[, (so + nA) :+ bcols])')'
                    xb_b = (sol[1..B, ] + sol[(B + 1)..(2 * B), ]) +
                        (sol[(2 * B + 1)..(3 * B), ] + sol[(3 * B + 1)..(4 * B), ])
                    xb_tot[, bcols] = xb_b
                }
                if (apos > 0) {
                    rajq_b = JS[, so + apos]
                    xa_tot[, apos] = (rad_b + raa_b + rajq_b) :/ a.blocks[bid].mass
                    if (pob > 0) xa_tot[, apos] = xa_tot[, apos] - xb_b * a.blocks[bid].mean_d'
                }
            }
        }
        else {
            rbd = J(B, po, 0)
            rad = J(B, nA, 0)
            for (t=1; t<=rows(plan.bids); t++) {
                bid = plan.bids[t]
                tb = plan.or_bpos[bid]
                rb = a.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                L = a.blocks[bid].L
                for (k=1; k<=a.blocks[bid].n_ctrl; k++) {
                    cl = a.blocks[bid].ctrl_cell[k]
                    ck = a.blocks[bid].ctrl_ci[k]
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) == 0) continue
                    v = va[ix, (rb - 1) * G + ck]
                    rv = v :* (Zc[ix, plan.zpos[cl]] - Mc[ix, tb])
                    if (apos > 0) rad[ix, apos] = rad[ix, apos] + rv
                    if (po > 0) {
                        dor = L :* Z_rc[ix, (rb - 1) * p :+ a.or_keep]
                        dt = (a.intercept ? dor :- a.blocks[bid].mean_d : dor)
                        rbd[ix, ] = rbd[ix, ] + rv :* dt
                    }
                }
            }
            raa = (nA > 0 ? xa * plan.Raa_load : J(B, 0, 0))
            rba = (po > 0 ? xa * plan.Rba_load : J(B, 0, 0))
            if (po > 0) {
                sol = _cdid_aug_or_solve(a,
                    rbd \ rba \ Jb[, (so + nA + 1)..(so + nA + po)] \ S[, (so + nA + 1)..(so + nA + po)])
                xb_tot = (sol[1..B, ] + sol[(B + 1)..(2 * B), ]) +
                    (sol[(2 * B + 1)..(3 * B), ] + sol[(3 * B + 1)..(4 * B), ])
            }
            if (nA > 0) {
                rajq = JS[, (so + 1)..(so + nA)]
                for (aa=1; aa<=nA; aa++) {
                    bid = plan.abids[aa]
                    xa_tot[, aa] = (rad[, aa] + raa[, aa] + rajq[, aa]) :/ a.blocks[bid].mass
                    if (po > 0) xa_tot[, aa] = xa_tot[, aa] - xb_tot * a.blocks[bid].mean_d'
                }
            }
        }
        if (nA > 0) o_phi = o_phi + xa_tot * plan.Alpha_load
        if (tot_po > 0) o_phi = o_phi + xb_tot * plan.Beta_load
    }

    for (mi=1; mi<=nl; mi++) {
        ell = plan.lk[mi]
        rl = a.links[ell].recipe
        g = a.links[ell].target_cell
        cg = a.links[ell].target_ci
        ixg = _cdid_impb_cix(c_ix, m, fr, rl, cg)
        if (!a.has_or) {
            if (rows(ixg) > 0) {
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ a.links[ell].d_plus) :*
                    (Zc[ixg, plan.zpos[g]] :- a.links[ell].T_gl))
            }
            tp = plan.link_tp[mi]
            for (k=1; k<=a.links[ell].n_ctrl; k++) {
                ck = a.links[ell].ctrl_ci[k]
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) == 0) continue
                cl = a.links[ell].ctrl_cell[k]
                o_phi[ixk, mi] = (o_phi[ixk, mi] -
                    (wm_cache[ixk, tp] :/ a.links[ell].d_minus) :*
                    (Zc[ixk, plan.zpos[cl]] :- a.links[ell].C_ol))
            }
        }
        else if (!a.has_ps) {
            if (rows(ixg) > 0) {
                tb = plan.or_bpos[a.links[ell].or_block]
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ a.links[ell].d_plus) :*
                    ((Zc[ixg, plan.zpos[g]] :- a.links[ell].T_gl) -
                     (Mc[ixg, tb] :- a.links[ell].B_gl)))
            }
        }
        else {
            tb = plan.or_bpos[a.links[ell].or_block]
            if (rows(ixg) > 0) {
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ a.links[ell].d_plus) :*
                    ((Zc[ixg, plan.zpos[g]] :- a.links[ell].T_gl) -
                     (Mc[ixg, tb] :- a.links[ell].B_gl)))
            }
            tp = plan.link_tp[mi]
            for (k=1; k<=a.links[ell].n_ctrl; k++) {
                ck = a.links[ell].ctrl_ci[k]
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) == 0) continue
                cl = a.links[ell].ctrl_cell[k]
                o_phi[ixk, mi] = (o_phi[ixk, mi] -
                    (wm_cache[ixk, tp] :/ a.links[ell].d_minus) :*
                    ((Zc[ixk, plan.zpos[cl]] :- a.links[ell].C_ol) -
                     (Mc[ixk, tb] :- a.links[ell].B_ol)))
            }
        }
    }
    if (missing(o_phi)) _cdid_aug_fail("Nonfinite complete augmented link IF.", 504)
}

/* Row replay. f_rows: unique original-F positions (any order, any count; bounded internally
   by chunk_rows). links: unique backend link indices, processed in batches of link_batch with
   ONE Stage 1 adjoint and ONE theta-only contraction record per batch. scale frame (phi_F,
   normalization N_F) or core ((n_A/N_F) phi_F after summation, normalization n_A). detail = 1
   adds the ten source channels and the four reference sums (same scale). */
struct cdid_aug_backend_if scalar cdid_aug_backend_rows(struct cdid_aug_backend scalar a,
    real colvector f_rows, real colvector links, | string scalar scale, real scalar detail)
{
    struct cdid_aug_backend_if scalar res
    struct _cdid_aug_plan scalar plan
    struct cdid_imputation_adjoint scalar adj
    struct cdid_sel_moments scalar sm
    pointer(struct cdid_imp_backend scalar) scalar pb
    pointer(struct cdid_imputation_model scalar) scalar pm
    pointer(struct cdid_sel_reader scalar) scalar pq
    real scalar B, nl, dmode, blo, bhi, lo, hi, step
    real colvector lk, rch, abids, bbids
    real matrix t_phi, t_comp, t_ref
    string scalar sc

    sc = (args() < 4 | scale == "" ? "frame" : scale)
    dmode = (args() < 5 ? 0 : detail)
    if (sc != "frame" & sc != "core") _cdid_aug_fail("IF scale must be frame or core.", 498)
    if (dmode != 0 & dmode != 1) _cdid_aug_fail("detail must be 0 or 1.", 498)
    _cdid_aug_check(a)
    pb = a.ipw
    pm = a.model
    pq = (a.is_dated ? &(a.qlease) : a.qrdr)
    B = rows(f_rows)
    if (B < 1 | missing(f_rows) | any(f_rows :< 1) | any(f_rows :> a.N) | any(f_rows :!= floor(f_rows)))
        _cdid_aug_fail("Replay rows must be original-F positions 1..N.", 498)
    if (B > 1) {
        if (!all(f_rows[2..B] :> f_rows[1..(B - 1)])) {
            if (rows(uniqrows(f_rows)) != B) _cdid_aug_fail("Replay rows must not repeat an original-F unit.", 498)
        }
    }
    nl = rows(links)
    if (nl < 1 | missing(links) | any(links :< 1) | any(links :> a.n_links) | any(links :!= floor(links)))
        _cdid_aug_fail("Invalid backend link indices.", 498)
    if (nl > 1) {
        if (!all(links[2..nl] :> links[1..(nl - 1)])) {
            if (rows(uniqrows(links)) != nl) _cdid_aug_fail("A link request must not repeat a link.", 498)
        }
    }
    _cdid_aug_replay_budget(a, min((B, a.chunk_rows)), nl, B, dmode)
    res.phi = J(B, nl, 0)
    res.comp = J(0, 0, .)
    res.ref4 = J(0, 0, .)
    if (!dmode) {
        step = (_cdid_aug_can_full_batch(a, min((B, a.chunk_rows)), nl, B) ? nl : a.link_batch)
        for (blo=1; blo<=nl; blo=blo+step) {
            bhi = min((nl, blo + step - 1))
            lk = links[blo..bhi]
            plan = _cdid_aug_make_plan(a, *pb, *pm, lk)
            for (lo=1; lo<=B; lo=lo+a.chunk_rows) {
                hi = min((B, lo + a.chunk_rows - 1))
                rch = f_rows[lo..hi]
                _cdid_aug_eval_fast(a, *pb, *pm, rch, plan, t_phi=.)
                res.phi[lo..hi, blo..bhi] = t_phi
            }
        }
    }
    else {
        res.comp = J(B, 10 * nl, .)
        res.ref4 = J(B, 4 * nl, .)
        for (blo=1; blo<=nl; blo=blo+a.link_batch) {
            bhi = min((nl, blo + a.link_batch - 1))
            lk = links[blo..bhi]
            abids = (a.has_or & a.intercept & a.method == "reg" ? _cdid_aug_abids(a, lk) : J(0, 1, 0))
            bbids = (a.has_or & a.outcome_pooling == "separate" ? _cdid_aug_bbids(a, lk) : J(0, 1, 0))
            adj = _cdid_aug_batch_adj(a, *pm, lk, abids)
            adj = _cdid_aug_cat_or_adj(a, *pm, adj, bbids)
            sm = _cdid_aug_batch_sel(a, lk, abids, bbids)
            for (lo=1; lo<=B; lo=lo+a.chunk_rows) {
                hi = min((B, lo + a.chunk_rows - 1))
                rch = f_rows[lo..hi]
                _cdid_aug_eval(a, *pb, *pm, *pq, rch, lk, abids, bbids, adj, sm, dmode, t_phi=., t_comp=., t_ref=.)
                res.phi[lo..hi, blo..bhi] = t_phi
                res.comp[lo..hi, ((blo - 1) * 10 + 1)..(bhi * 10)] = t_comp
                res.ref4[lo..hi, ((blo - 1) * 4 + 1)..(bhi * 4)] = t_ref
            }
        }
    }
    if (sc == "core") {
        res.phi = a.lambda :* res.phi
        if (dmode) {
            res.comp = a.lambda :* res.comp
            res.ref4 = a.lambda :* res.ref4
        }
    }
    res.method = a.method
    res.scale = sc
    res.f_rows = f_rows
    res.links = links
    res.N = a.N
    res.n_A = a.n_A
    res.lambda = a.lambda
    res.normalization_count = (sc == "core" ? a.n_A : a.N)
    res.detail = dmode
    _cdid_aug_set_ids(res, *pm, f_rows)
    return(res)
}

/* L x 1 augmented primary points mu in backend link order. */
real colvector cdid_aug_backend_points(struct cdid_aug_backend scalar a)
{
    real colvector out
    real scalar ell
    out = J(a.n_links, 1, .)
    for (ell=1; ell<=a.n_links; ell++) out[ell] = a.links[ell].mu
    return(out)
}

void cdid_aug_backend_describe(struct cdid_aug_backend scalar a)
{
    printf("{txt}augmented backend %s: selection_est_method=dr est_method=%s j-pooling=%s q-model=%s q-pooling=%s N_F=%g n_A=%g links=%g cells=%g\n",
        a.backend_id, a.method, a.pooling, a.selection_model, a.selection_pooling, a.N, a.n_A, a.n_links, a.n_cells)
    printf("{txt}  owned fits: observation logit blocks=%g (theta=%g) Z-response OR=%g (blocks=%g, p_OR=%g, pooling=%s); borrowed: Stage 1 groups=%g, PS keys=%g\n",
        a.n_q_fits, a.n_q_params, a.n_or_fits, a.n_blocks, a.or_p, a.outcome_pooling, a.n_groups, a.n_trt)
}

/* GMM weight-provider context. consumer_map is M x 2: (backend link index, or 0 for a
   reference-zero moment; sign +1/-1), in original scientific-catalog order; an empty map
   means identity over the backend links. The context borrows the backend pointer. */
struct cdid_aug_backend_wctx scalar cdid_aug_backend_wctx_make(
    pointer(struct cdid_aug_backend scalar) scalar backend,
    real matrix consumer_map,
    | string colvector consumer_ids)
{
    struct cdid_aug_backend_wctx scalar ctx
    real scalar M, j, nl
    real matrix cm
    if (backend == NULL) _cdid_aug_fail("Weight provider requires a prepared augmented backend pointer.", 498)
    nl = (*backend).n_links
    cm = consumer_map
    if (rows(cm) == 0) cm = ((1::nl), J(nl, 1, 1))
    M = rows(cm)
    if (cols(cm) != 2 | missing(cm))
        _cdid_aug_fail("consumer_map must be a finite M x 2 matrix (link index or 0, sign).", 498)
    if (any(cm[, 1] :< 0) | any(cm[, 1] :> nl) | any(cm[, 1] :!= floor(cm[, 1])) | any(abs(cm[, 2]) :!= 1))
        _cdid_aug_fail("consumer_map needs link indices in 0..L (0 = reference zero) and signs +1/-1.", 498)
    if (!any(cm[, 1] :> 0)) _cdid_aug_fail("consumer_map must consume at least one backend link.", 498)
    ctx.backend = backend
    ctx.n_moments = M
    ctx.need_links = uniqrows(select(cm[, 1], cm[, 1] :> 0))
    ctx.col_pos = J(M, 1, 0)
    for (j=1; j<=M; j++) {
        if (cm[j, 1] > 0) ctx.col_pos[j] = _cdid_impb_pos(ctx.need_links, cm[j, 1])
    }
    ctx.col_sign = cm[, 2]
    ctx.consumer_ids = J(0, 1, "")
    if (args() >= 3) ctx.consumer_ids = consumer_ids
    if (rows(ctx.consumer_ids) != 0 & rows(ctx.consumer_ids) != M)
        _cdid_aug_fail("consumer_ids must be empty or one per consumer moment.", 498)
    ctx.row_chunk = (*backend).chunk_rows
    ctx.fast_ready = 0
    if ((*backend).ipw != NULL & (*backend).model != NULL &
        _cdid_aug_can_full_batch(*backend, ctx.row_chunk, rows(ctx.need_links), ctx.row_chunk)) {
        ctx.plan = _cdid_aug_make_plan(*backend, *((*backend).ipw), *((*backend).model), ctx.need_links)
        ctx.fast_ready = 1
    }
    return(ctx)
}

/* GMM provider weight(context, lo, hi): (hi-lo+1) x M CORE IF rows lo..hi of the original F
   roster in consumer order (signs and reference zeros included). Every one of the N_F rows
   is replayable; normalization_count = n_A. Point-only GMM still consumes these rows for its
   optimal/two-step weight; only the final target IF/inference is suppressed downstream. */
real matrix cdid_aug_backend_weight(pointer() scalar context, real scalar lo, real scalar hi)
{
    struct cdid_aug_backend_wctx scalar ctx
    struct cdid_aug_backend_if scalar r
    pointer(struct cdid_aug_backend scalar) scalar pb
    real matrix out, phi, t_phi
    real scalar j, N, B, rlo, rhi
    if (context == NULL) _cdid_aug_fail("Missing prepared augmented weight-provider context.", 498)
    ctx = *context
    pb = ctx.backend
    if (pb == NULL) _cdid_aug_fail("Weight-provider backend pointer is NULL; its lifetime ended.", 498)
    N = (*pb).N
    if (lo < 1 | hi < lo | hi > N | lo != floor(lo) | hi != floor(hi))
        _cdid_aug_fail("Weight-provider rows must be a contiguous block of original F rows.", 498)
    B = hi - lo + 1
    if (ctx.fast_ready & _cdid_aug_can_full_batch(*pb, min((B, (*pb).chunk_rows)), rows(ctx.need_links), B)) {
        _cdid_aug_check(*pb)
        _cdid_aug_replay_budget(*pb, min((B, (*pb).chunk_rows)), rows(ctx.need_links), B, 0)
        phi = J(B, rows(ctx.need_links), 0)
        for (rlo=lo; rlo<=hi; rlo=rlo+(*pb).chunk_rows) {
            rhi = min((hi, rlo + (*pb).chunk_rows - 1))
            _cdid_aug_eval_fast(*pb, *((*pb).ipw), *((*pb).model), (rlo::rhi), ctx.plan, t_phi=.)
            phi[(rlo - lo + 1)..(rhi - lo + 1), ] = (*pb).lambda :* t_phi
        }
    }
    else {
        r = cdid_aug_backend_rows(*pb, (lo::hi), ctx.need_links, "core", 0)
        phi = r.phi
    }
    if (all(ctx.col_pos :> 0)) return(phi[, ctx.col_pos] :* (ctx.col_sign'))
    out = J(B, ctx.n_moments, 0)
    for (j=1; j<=ctx.n_moments; j++) {
        if (ctx.col_pos[j] > 0) out[, j] = ctx.col_sign[j] :* phi[, ctx.col_pos[j]]
    }
    return(out)
}

/* Complete-unit weight Gram through the GMM constructor: N_F replay rows, centered
   over F, divided by n_A (never a layout or off-A row count). */
struct cdid_gmm_unit_gram scalar cdid_aug_backend_gram(pointer() scalar context,
    real scalar chunk, real scalar cap)
{
    struct cdid_aug_backend_wctx scalar ctx
    pointer(struct cdid_aug_backend scalar) scalar pb
    if (context == NULL) _cdid_aug_fail("Missing prepared augmented weight-provider context.", 498)
    ctx = *context
    pb = ctx.backend
    if (pb == NULL) _cdid_aug_fail("Weight-provider backend pointer is NULL; its lifetime ended.", 498)
    return(_cdid_gmm_unit_gram_norm((*pb).N, ctx.n_moments, (*pb).n_A,
        chunk, &cdid_aug_backend_weight(), context, cap))
}

/* Prepare the static augmented backend.
   Concrete input seam (static): ipw = &IB with IB = cdid_imp_backend_prepare(&M, &XT, links,
   controls, ids, "ipw", glm_eps, glm_maxit, link_batch, budget) on a static Stage 1
   fit M = cdid_imputation_fit_static(pop, W, Zslopes, pairs, requests, j_pooling, intercept,
   rank_tol, budget) over the complete post + both-auxiliary link registry, XT the full-F
   xformla treatment design; qrdr = &QR with QR = cdid_sel_reader_pop(pop, Q, has_int, int_col,
   W, decl, J(N,1,1), chunk) for decl = cdid_aug_route_gate(...) and Q the full-F
   selection_xformla design (NULL inherits xformla). The same raw W, population, calendar and
   cohorts are verified exactly. GLM controls are the ipw backend's (one source glm_control).
   Recipe-indexed designs (refused above) would borrow
   cdid_imp_backend_prepare_d(..., "ipw", ...) -- the row accessors used here already
   dispatch to its dated reader, and Delta reads use the same Stage 1 pair reader -- and
   replace the single q context by one recipe-local q context per recipe on a
   USED-COHORT observation reader over original F (cdid_sel_reader_pop refuses missing
   unused-cohort slots by design; no fill, crop, fake universe or known cohort). */
struct cdid_aug_backend scalar cdid_aug_backend_prepare(
    pointer(struct cdid_imp_backend scalar) scalar ipw,
    pointer(struct cdid_sel_reader scalar) scalar qrdr,
    real scalar link_batch,
    struct cdid_imputation_budget scalar budget,
    | string scalar backend_id,
      string scalar method_in,
      string scalar outcome_pooling_in)
{
    struct cdid_aug_backend scalar a
    struct cdid_sel_cell_map scalar qmap
    struct _cdid_aug_rcc rowvector rcc
    struct _cdid_aug_cev rowvector cev
    pointer(struct cdid_imputation_model scalar) scalar pm
    real scalar ci, max_coh, c_bytes, use_cache
    string scalar mstr, op
    if (ipw == NULL)
        _cdid_aug_fail("Augmented preparation requires a first-step backend pointer.", 498)
    if (qrdr == NULL)
        _cdid_aug_fail("Augmented preparation requires an observation-probability reader pointer.", 498)
    if (link_batch < 1 | link_batch >= . | link_batch != floor(link_batch))
        _cdid_aug_fail("link_batch must be a positive integer.", 498)
    mstr = (args() >= 6 & strtrim(method_in) != "" ? strlower(strtrim(method_in)) : "dr")
    if (mstr != "ipw" & mstr != "reg" & mstr != "dr")
        _cdid_aug_fail("est_method must be ipw, reg, or dr.", 459)
    cdid_imputation_validate_budget(budget)
    _cdid_impb_need_model(*ipw)
    pm = (*ipw).model
    _cdid_impb_check(*ipw, *pm)
    a.is_dated = 0
    a.ipw = ipw
    a.model = pm
    a.qrdr = qrdr
    a.budget = budget
    a.link_batch = link_batch
    a.backend_id = ""
    if (args() >= 5) a.backend_id = backend_id
    a.method = mstr
    a.has_or = (mstr == "reg" | mstr == "dr")
    a.has_ps = (mstr == "ipw" | mstr == "dr")
    if (a.has_or) {
        op = (args() >= 7 ? strlower(strtrim(outcome_pooling_in)) : "")
        if (op == "") op = "separate"
        if (op != "separate" & op != "pooled")
            _cdid_aug_fail("outcome_pooling must be separate or pooled.", 198)
        a.outcome_pooling = op
    }
    else a.outcome_pooling = ""
    a.selection_est_method = "dr"
    a.selection = "conditional_mtar"
    _cdid_aug_init(a, *ipw, *pm, *qrdr)
    _cdid_aug_registry(a, *ipw, *pm)
    qmap = _cdid_aug_q_map(a, *pm, *qrdr)
    _cdid_aug_preflight(a, *ipw, *pm, *qrdr, qmap)
    _cdid_aug_q_fit(a, *qrdr, qmap)
    max_coh = 0
    for (ci=1; ci<=a.G; ci++) max_coh = max((max_coh, rows((*ipw).coh_u[ci].rows)))
    c_bytes = 8 * (a.n_recipes * a.N * (a.p_raw + (a.has_ps ? a.q : 0) + 2) +
        a.n_cells * max_coh * (2 * a.d_sel + 8))
    use_cache = (max_coh <= a.chunk_rows & a.retained_bytes + c_bytes <= a.budget.max_workspace_bytes)
    rcc = _cdid_aug_rcc(1, a.n_recipes * a.G)
    cev = _cdid_aug_cev(1, a.n_cells)
    for (ci=1; ci<=cols(rcc); ci++) {
        rcc[ci].ready = 0
        rcc[ci].has_x = 0
    }
    for (ci=1; ci<=a.n_cells; ci++) cev[ci].ready = 0
    _cdid_aug_cell_moments(a, *ipw, *pm, *qrdr, rcc, cev, use_cache)
    if (a.has_or) {
        _cdid_aug_fit_or(a, *ipw, *pm, *qrdr, rcc, cev)
        _cdid_aug_or_forcing(a, *ipw, *pm, *qrdr, rcc, cev)
    }
    _cdid_aug_link_sums(a, *ipw, *pm, *qrdr, rcc, cev)
    if (a.has_or) _cdid_aug_or_adjoint(a, *pm)
    a.ipw_sum = _cdid_aug_ipw_sum(*ipw)
    a.q_sum = _cdid_aug_q_sum(*qrdr)
    return(a)
}

/* Dated q leaf intake. qdata contains ONE matrix per canonical Q date, not one
   per recipe. The caller keeps dr/reg/model/ipw immutable. This backend owns qdata
   after value transfer and leases one slot at a time; there is no layout population. */
struct cdid_aug_backend scalar cdid_aug_backend_prepare_d(
    pointer(struct cdid_imp_backend scalar) scalar ipw,
    pointer(struct cdid_imp_dated_reader scalar) scalar dr,
    struct cdid_imp_matrix_slot colvector qdata, real colvector recipe_qslot,
    real scalar q_intercept, real scalar q_int_col, struct cdid_sel_decl scalar decl,
    real scalar link_batch, struct cdid_imputation_budget scalar budget,
    | string scalar backend_id,
      string scalar method_in,
      string scalar outcome_pooling_in)
{
    struct cdid_aug_backend scalar a
    struct cdid_sel_reader scalar qr
    struct cdid_sel_cell_map scalar qmap
    struct cdid_sel_budget scalar sb
    struct _cdid_aug_rcc rowvector rcc
    struct _cdid_aug_cev rowvector cev
    pointer(struct cdid_imputation_model scalar) scalar pm
    real colvector cohort_scope, cids, seen
    real matrix pspec
    real scalar r, j, cid, pc, hit, k, slot, ret, scratch, work, off, max_coh, c_bytes, use_cache
    string scalar mstr, op
    if (ipw == NULL | dr == NULL) _cdid_aug_fail("First-step and dated reader objects are required.", 498)
    mstr = (args() >= 11 & strtrim(method_in) != "" ? strlower(strtrim(method_in)) : "dr")
    if (mstr != "ipw" & mstr != "reg" & mstr != "dr")
        _cdid_aug_fail("est_method must be ipw, reg, or dr.", 459)
    cdid_imputation_validate_budget(budget)
    _cdid_impd_check_reader(*dr)
    _cdid_impb_need_model(*ipw)
    pm = (*ipw).model
    _cdid_impb_check(*ipw, *pm)
    if (!(*ipw).is_dated | (*ipw).rdr != dr)
        _cdid_aug_fail("Dated augmentation requires the same dated first-step reader.", 498)
    if (rows(recipe_qslot) != (*pm).n_recipes | cols(recipe_qslot) != 1 | missing(recipe_qslot) |
        any(recipe_qslot :< 1) | any(recipe_qslot :> rows(qdata)) | any(recipe_qslot :!= floor(recipe_qslot)))
        _cdid_aug_fail("Each dated design block requires one valid canonical Q slot.", 498)
    if (link_batch < 1 | link_batch >= . | link_batch != floor(link_batch))
        _cdid_aug_fail("link_batch must be a positive integer.", 498)
    a.is_dated=1; a.ipw=ipw; a.model=pm; a.qrdr=NULL; a.qdated=dr
    a.budget=budget; a.link_batch=link_batch
    a.backend_id=(args()<10 ? "" : backend_id)
    a.method=mstr; a.selection_est_method="dr"; a.selection="conditional_mtar"
    a.has_or=(mstr=="reg" | mstr=="dr")
    a.has_ps=(mstr=="ipw" | mstr=="dr")
    if (a.has_or) {
        op = (args() >= 12 ? strlower(strtrim(outcome_pooling_in)) : "")
        if (op == "") op = "separate"
        if (op != "separate" & op != "pooled")
            _cdid_aug_fail("outcome_pooling must be separate or pooled.", 198)
        a.outcome_pooling = op
    }
    else a.outcome_pooling = ""
    a.n_recipes=(*pm).n_recipes
    a.qcontexts=_cdid_aug_qcontext(a.n_recipes,1)
    a.qdata=qdata
    a.qactive=0;a.qslot_active=0
    a.q_of_cell=J((*pm).n_cells,1,0)
    a.n_q_params=0;a.n_q_pc=0;a.n_q_fits=0
    ret=0;scratch=0;work=0
    for(j=1;j<=rows(qdata);j++) {
        if(rows(qdata[j].value)!=(*pm).N) _cdid_aug_fail("Q slots require original-F rows.",498)
        ret=ret+length(qdata[j].value)
    }
    for(r=1;r<=a.n_recipes;r++) {
        slot=recipe_qslot[r]
        cohort_scope=J((*pm).n_cohorts,1,0)
        /* Selection leaf orders never-treated last; Stage1 orders it first. */
        for(k=1;k<=(*pm).n_cohorts;k++) {
            j=((*pm).cohorts[k]==0 ? (*pm).n_cohorts : sum(((*pm).cohorts:>0) :& ((*pm).cohorts:<(*pm).cohorts[k]))+1)
            cohort_scope[j]=(*dr).used_cohorts[r,k]
        }
        if(rows((*pm).elig_keep)) qr=cdid_sel_reader_scope(*(*dr).pop,qdata[slot].value,q_intercept,q_int_col,
            (*(*dr).reg).w_data[(*dr).recipes[r].w_slot].value,decl,J((*pm).N,1,1),budget.max_chunk_rows,cohort_scope,
            (*pm).elig_keep[,r])
        else qr=cdid_sel_reader_scope(*(*dr).pop,qdata[slot].value,q_intercept,q_int_col,
            (*(*dr).reg).w_data[(*dr).recipes[r].w_slot].value,decl,J((*pm).N,1,1),budget.max_chunk_rows,cohort_scope)
        if(r==1) {
            _cdid_aug_init(a,*ipw,*pm,qr)
            _cdid_aug_registry(a,*ipw,*pm)
            a.n_q_params=0;a.n_q_pc=0;a.n_q_fits=0
            a.qlease=qr;a.qlease.Z=J(0,0,0)
        }
        a.qcontexts[r].q_slot=slot
        a.qcontexts[r].w_slot=(*dr).recipes[r].w_slot
        a.qcontexts[r].spec=qr.spec;a.qcontexts[r].n_U=qr.n_U
        a.qcontexts[r].mass_F=qr.mass_F;a.qcontexts[r].mass_U=qr.mass_U
        a.qcontexts[r].repeat_mass=qr.repeat_mass;a.qcontexts[r].a_U=qr.a_U
        a.qcontexts[r].U=(rows((*pm).elig_keep) ? qr.U : J(0,1,.))
        cids=J(0,1,0)
        for(cid=1;cid<=a.n_cells;cid++)if((*pm).cells[cid].recipe==r)cids=cids\cid
        pspec=J(rows(cids),3,0)
        for(j=1;j<=rows(cids);j++) {
            cid=cids[j]
            pspec[j,]=(a.sel_cohort_pos[(*pm).cells[cid].cohort_index],(*pm).cells[cid].s_index,(*pm).cells[cid].t_index)
        }
        qmap=cdid_sel_partition(qr,pspec)
        if(qmap.n_pair_cohorts!=rows(cids)) _cdid_aug_fail("Recipe q and j cells must match one-to-one.",498)
        sb=cdid_sel_budget_check(qr,qmap,.)
        ret=ret+sb.compact_state_doubles+9*qmap.n_cells+8*qmap.n_pair_cohorts+5*qmap.n_blocks+12*qmap.n_pair_cohorts+6*a.G
        scratch=max((scratch,sb.stream_workspace_doubles))
        work=work+16*sb.work_units*(a.glm_maxit+2)
        /* Charge compact q state plus single live reader before fitting. The full
           all-block outcome/replay preflight below charges its scientific state. */
        _cdid_imp_check_budget(budget,(*ipw).retained_bytes+8*ret,
            8*(scratch+(*pm).N*(cols(qr.Z)+12))+4096,work,"dated augmented q preparation")
        a.qcontexts[r].fit=cdid_sel_fit(qr,qmap,a.glm_eps,a.glm_maxit)
        a.qcontexts[r].up=cdid_sel_up_map_build(a.qcontexts[r].fit)
        off=a.n_q_params;a.qcontexts[r].theta_offset=off
        a.n_q_params=a.n_q_params+a.qcontexts[r].fit.total_params
        a.n_q_fits=a.n_q_fits+a.qcontexts[r].fit.n_blocks
        a.n_q_pc=a.n_q_pc+qmap.n_pair_cohorts
        seen=J(qmap.n_pair_cohorts,1,0)
        for(j=1;j<=rows(cids);j++) {
            cid=cids[j];hit=0
            for(pc=1;pc<=qmap.n_pair_cohorts;pc++)if(all(qmap.pair_cohorts[pc,1..3]:==pspec[j,]))hit=pc
            if(hit==0) _cdid_aug_fail("A design block imputation cell has no q fit.",498)
            if(seen[hit]) _cdid_aug_fail("Two design block imputation cells resolve to one q fit.",498)
            seen[hit]=1;a.q_of_cell[cid]=hit;a.cells[cid].q_pc=hit
            a.cells[cid].theta_ids=off:+_cdid_sel_pc_global_idx(a.qcontexts[r].fit,a.qcontexts[r].up,hit)
        }
    }
    /* The static all-block estimate is conservative in scientific cell/
       link state. Add every recipe's compact q state and unique Q storage; keep
       its largest scan/replay scratch and explicit work refusal. */
    _cdid_aug_preflight(a,*ipw,*pm,qr,qmap)
    a.retained_bytes=a.retained_bytes+8*ret
    a.scratch_bytes=max((a.scratch_bytes,8*(scratch+(*pm).N*(cols(qr.Z)+12))+4096))
    a.work_ops=a.work_ops+work
    _cdid_imp_check_budget(budget,a.retained_bytes,a.scratch_bytes,a.work_ops,"dated augmented complete preparation")
    _cdid_aug_q_activate(a,1)
    max_coh=0
    for(k=1;k<=a.G;k++)max_coh=max((max_coh,rows((*ipw).coh_u[k].rows)))
    c_bytes=8*(a.n_recipes*a.N*(a.p_raw+(a.has_ps ? a.q : 0)+2)+a.n_cells*max_coh*(2*a.d_sel+8))
    use_cache=(max_coh<=a.chunk_rows & a.retained_bytes+c_bytes<=a.budget.max_workspace_bytes)
    rcc=_cdid_aug_rcc(1,a.n_recipes*a.G)
    cev=_cdid_aug_cev(1,a.n_cells)
    for(k=1;k<=cols(rcc);k++) {rcc[k].ready=0;rcc[k].has_x=0;}
    for(k=1;k<=a.n_cells;k++)cev[k].ready=0
    _cdid_aug_cell_moments(a,*ipw,*pm,a.qlease,rcc,cev,use_cache)
    if(a.has_or) {
        _cdid_aug_fit_or(a,*ipw,*pm,a.qlease,rcc,cev)
        _cdid_aug_or_forcing(a,*ipw,*pm,a.qlease,rcc,cev)
    }
    _cdid_aug_link_sums(a,*ipw,*pm,a.qlease,rcc,cev)
    if(a.has_or) _cdid_aug_or_adjoint(a,*pm)
    a.ipw_sum=_cdid_aug_ipw_sum(*ipw);a.q_sum=.
    return(a)
}

end
