*! Private selected (conditional_mtar) conditional-pretrend process piece producer; StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (commit ):
*! R/chained_conditional_extract.R (cpt_piece, cpt_two_arm_piece,
*! cpt_two_arm_rows, cpt_extract_two_arm, cpt_extract_or)
*! and R/chained_conditional_pretest.R (cpt_declared_cells).
*! Supports all 3 selected observation-weighted methods (ipw, reg, dr) across
*! selection_model (pair_logit, independent_period_logit), selection_pooling
*! (separate, within_cohort, across_cohorts), and selection_universe (F, B, A)
*! for both static (cdid_selected_model) and dated (cdid_selref_provider) owners.
*! Reuses the complete post + both-auxiliary union fit without refitting selection,
*! propensity, or pooled WLS nuisances by threshold or chunk.
version 19.5
mata:

struct cdid_cpt_sel_cell_meta {
    string scalar id, status, reason, route
    real scalar g, t, from, s_date, t2_date, s_index, t2_index
    real scalar q, link_idx, recipe, design_index, sign, supported
    real scalar estimate, estimate_scalar, n_chunks, chunk_rows
    real scalar n_up, p_trt, p_or, or_icpt
    real colvector up_channels, up_blocks
    struct cdid_cpt_spec scalar spec
}

struct cdid_cpt_selected {
    string scalar kind, method, selection_model, selection_pooling
    string scalar selection_universe, base_period, status, reason
    real scalar n_F, n_A, n_cells, n_unavailable, chunk_rows
    real scalar workspace_entries, workspace_limit, work_ops, work_limit
    real scalar n_sel_fits, n_trt_fits, n_or_fits
    pointer(struct cdid_selected_model scalar) scalar static_provider
    pointer(struct cdid_selref_provider scalar) scalar dated_provider
    struct cdid_cpt_sel_cell_meta colvector cells
}

void _cdid_cpt_sel_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_cpt_sel_ix(real matrix pred)
{
    if (rows(pred) * cols(pred) == 0) return(J(0, 1, 0))
    return(vec(selectindex(pred)))
}

/* Preserve first-seen order of integer IDs (matching R unique()). */
real colvector _cdid_cpt_sel_uniq(real colvector x)
{
    real colvector out
    real scalar n, i, k, seen
    n = rows(x)
    if (n == 0) return(J(0, 1, 0))
    out = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        seen = (k > 0 ? any(out[1..k] :== x[i]) : 0)
        if (!seen) {
            k++
            out[k] = x[i]
        }
    }
    return(k > 0 ? out[1..k] : J(0, 1, 0))
}

/* Match each element of query against table (1-based position; 0 if absent). */
real colvector _cdid_cpt_sel_match(real colvector query, real colvector table)
{
    real colvector out, hit
    real scalar n, i
    n = rows(query)
    out = J(n, 1, 0)
    if (n == 0 | rows(table) == 0) return(out)
    for (i = 1; i <= n; i++) {
        hit = _cdid_cpt_sel_ix(table :== query[i])
        if (rows(hit) > 0) out[i] = hit[1]
    }
    return(out)
}

/* Upstream channel IDs for one pair-cohort pc in recipe r:
   Channel 1..n_cohorts = cohort repeat-rate a_k^U;
   Channel n_cohorts + (1..total_params) = selection logit parameter theta_g.
   For pair-cohort pc, channel 1 is cohort k = um.pc_info[pc, 1] followed by
   n_cohorts :+ _cdid_sel_pc_global_idx(fit, um, pc), matching R
   selection_pair_cohort_log_deriv (ld.d_ak followed by ld.d_theta_local). */
real colvector _cdid_cpt_sel_pc_channels(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_up_map scalar um, real scalar pc)
{
    real scalar k
    real colvector g_idx, ch
    k = um.pc_info[pc, 1]
    g_idx = _cdid_sel_pc_global_idx(fit, um, pc)
    ch = k
    if (rows(g_idx) > 0) ch = ch \ (um.n_cohorts :+ g_idx)
    return(ch)
}

/* Active selection block IDs for a set of pair-cohort indices pc_ids. */
real colvector _cdid_cpt_sel_pc_blocks(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_up_map scalar um, real colvector pc_ids)
{
    real colvector raw_b
    real scalar i, pc, b_s, b_t
    raw_b = J(0, 1, 0)
    for (i = 1; i <= rows(pc_ids); i++) {
        pc = pc_ids[i]
        b_s = um.pc_info[pc, 7]
        b_t = um.pc_info[pc, 9]
        if (b_s >= 1 & b_s <= fit.n_blocks) {
            if (fit.blocks[b_s].p_m > 0) raw_b = raw_b \ b_s
        }
        if (fit.selection_model != "pair_logit" & b_t != b_s & b_t >= 1 & b_t <= fit.n_blocks) {
            if (fit.blocks[b_t].p_m > 0) raw_b = raw_b \ b_t
        }
    }
    return(_cdid_cpt_sel_uniq(raw_b))
}

/* Deduplicated upstream channel IDs across a treatment model's pair-cohorts. */
real colvector _cdid_cpt_sel_trt_channels(struct cdid_sel_fit_result scalar fit,
    struct cdid_sel_up_map scalar um, real colvector pc_ids)
{
    real colvector raw_ch
    real scalar k
    raw_ch = J(0, 1, 0)
    for (k = 1; k <= rows(pc_ids); k++) {
        raw_ch = raw_ch \ _cdid_cpt_sel_pc_channels(fit, um, pc_ids[k])
    }
    return(_cdid_cpt_sel_uniq(raw_ch))
}

/* Coordinate labels matching R chained_conditional_extract.R lines 255-256 and 297. */
string colvector _cdid_cpt_sel_coords_two_arm(real scalar is_dr,
    real scalar n_up, real scalar p_trt, real scalar or_icpt,
    real scalar reg_block, real scalar p_or, | real scalar is_sep)
{
    string colvector out
    string scalar b_pfx
    real scalar j, at, total, sep
    sep = (args() >= 7 & is_sep != 0 ? 1 : 0)
    total = n_up + p_trt + (is_dr & or_icpt ? 1 : 0) + (is_dr ? p_or : 0)
    out = J(total, 1, "")
    at = 0
    for (j = 1; j <= n_up; j++) {
        at++
        out[at] = "up:" + strtrim(strofreal(j, "%21.0f"))
    }
    for (j = 1; j <= p_trt; j++) {
        at++
        out[at] = "ps:" + strtrim(strofreal(j, "%21.0f"))
    }
    if (is_dr) {
        if (or_icpt) {
            at++
            out[at] = "or_alpha:" + strtrim(strofreal(reg_block, "%21.0f"))
        }
        b_pfx = (sep ? strtrim(strofreal(reg_block, "%21.0f")) + ":" : "")
        for (j = 1; j <= p_or; j++) {
            at++
            out[at] = "or_beta:" + b_pfx + strtrim(strofreal(j, "%21.0f"))
        }
    }
    return(out)
}

string colvector _cdid_cpt_sel_coords_or(real colvector up_channels,
    real scalar or_icpt, real scalar reg_block, real scalar p_or, | real scalar is_sep)
{
    string colvector out
    string scalar b_pfx
    real scalar n_ch, j, at, total, sep
    sep = (args() >= 5 & is_sep != 0 ? 1 : 0)
    n_ch = rows(up_channels)
    total = n_ch + (or_icpt ? 1 : 0) + p_or
    out = J(total, 1, "")
    at = 0
    for (j = 1; j <= n_ch; j++) {
        at++
        out[at] = "up:" + strtrim(strofreal(up_channels[j], "%21.0f"))
    }
    if (or_icpt) {
        at++
        out[at] = "or_alpha:" + strtrim(strofreal(reg_block, "%21.0f"))
    }
    b_pfx = (sep ? strtrim(strofreal(reg_block, "%21.0f")) + ":" : "")
    for (j = 1; j <= p_or; j++) {
        at++
        out[at] = "or_beta:" + b_pfx + strtrim(strofreal(j, "%21.0f"))
    }
    return(out)
}

/* Bounded evaluation of upstream selection IF columns (a_k^U and theta_g) on
   original-F rowsF for recipe r, evaluating only the selection blocks in up_blocks. */
real matrix _cdid_cpt_sel_up_if(struct cdid_selected_model scalar m,
    real scalar r, real colvector up_channels, real colvector up_blocks,
    real colvector rowsF)
{
    pointer(struct cdid_sel_fit_result scalar) scalar fit
    pointer(struct cdid_sel_up_map scalar) scalar um
    real matrix xi_up, xa_m, xb_m
    real colvector u_sub, c_sub, w_sub, a_sub, ix_k
    real scalar n_req, n_up, n_coh, d_s, has_int, j, ch, k, b_idx, blk_id, g, c, bst
    n_req = rows(rowsF)
    n_up = rows(up_channels)
    xi_up = J(n_req, n_up, 0)
    if (n_req == 0 | n_up == 0) return(xi_up)
    _cdid_selected_activate_reader(m, r)
    fit = &(m.contexts[r].sel_fit)
    um = &(m.contexts[r].up_map)
    n_coh = (*um).n_cohorts
    d_s = (*fit).spec.d_s
    has_int = (*fit).spec.has_intercept

    u_sub = m.reader.U[rowsF]
    c_sub = m.reader.cohort_pos[rowsF]
    w_sub = m.reader.w[rowsF]
    a_sub = m.reader.A[rowsF]
    for (j = 1; j <= n_up; j++) {
        ch = up_channels[j]
        if (ch >= 1 & ch <= n_coh) {
            k = ch
            if ((*fit).mass_U[k] > 0) {
                ix_k = _cdid_cpt_sel_ix((u_sub :== 1) :& (c_sub :== k))
                if (rows(ix_k) > 0) {
                    xi_up[ix_k, j] = w_sub[ix_k] :* (a_sub[ix_k] :- (*fit).a_U[k]) :/ (*fit).mass_U[k]
                }
            }
        }
    }

    for (b_idx = 1; b_idx <= rows(up_blocks); b_idx++) {
        blk_id = up_blocks[b_idx]
        xa_m = J(0, 0, 0); xb_m = J(0, 0, 0)
        _cdid_sel_block_if_chunk(*fit, m.reader, blk_id, rowsF, xa_m, xb_m)
        bst = (*fit).beta_start[blk_id]
        for (j = 1; j <= n_up; j++) {
            ch = up_channels[j]
            if (ch <= n_coh) continue
            g = ch - n_coh
            if (g <= (*fit).n_alpha) {
                c = g
                if (has_int & (*fit).cell_block[c] == blk_id & cols(xa_m) > 0) {
                    xi_up[, j] = xa_m[, (*fit).cell_local[c]]
                }
            }
            else if (d_s > 0 & g >= bst & g < bst + d_s & cols(xb_m) > 0) {
                xi_up[, j] = xb_m[, g - bst + 1]
            }
        }
    }
    return(xi_up)
}

/* Bounded evaluation of untreated-change OR nuisance IFs (xi_alpha for block b
   and local/shared xi_beta) on original-F rowsF without scanning other blocks for objective_score. */
void _cdid_cpt_sel_or_if(struct cdid_selected_model scalar m, real scalar b,
    real colvector rowsF, real matrix xi_alpha_out, real matrix xi_beta_out)
{
    struct cdid_selected_loaded_block colvector empty_loads
    real matrix h_beta, eval_beta, alpha_score, hb_tmp
    real colvector obj_unused, u_mu, ha_tmp, xa_tmp
    real scalar n_req, p_kept, r, r_b
    n_req = rows(rowsF)
    if (m.outcome_pooling == "separate") {
        p_kept = m.or_blocks[b].p_kept
        xi_beta_out = J(n_req, p_kept, 0)
        xi_alpha_out = J(n_req, (m.or_shared.intercept ? 1 : 0), 0)
        if (n_req == 0 | (!m.or_shared.intercept & p_kept == 0)) return
        empty_loads = cdid_selected_loaded_block(0, 1)
        _cdid_sel_sep_block_if(m, b, rowsF, 0, empty_loads, m.or_shared.intercept,
            ha_tmp = ., hb_tmp = ., xa_tmp = ., xi_beta_out)
        if (m.or_shared.intercept) xi_alpha_out = xa_tmp
        return
    }
    p_kept = m.or_shared.p_kept
    xi_beta_out = J(n_req, p_kept, 0)
    xi_alpha_out = J(n_req, (m.or_shared.intercept ? 1 : 0), 0)
    if (n_req == 0 | (!m.or_shared.intercept & p_kept == 0)) return
    h_beta = (p_kept > 0 ? m.or_shared.beta_score_F[rowsF, ] : J(n_req, 0, 0))
    alpha_score = J(0, 0, 0); obj_unused = J(0, 1, 0)
    if (m.or_shared.intercept) {
        _cdid_selected_or_local_scores(m, rowsF, J(1, 1, b), 0, alpha_score, obj_unused)
    }
    if (p_kept > 0) {
        eval_beta = J(n_req, p_kept, 0)
        for (r = 1; r <= m.R; r++) {
            eval_beta = eval_beta + _cdid_selected_loading(m, r, rowsF,
                (m.contexts[r].m_col_beta_start::m.contexts[r].m_col_beta_end))
        }
        xi_beta_out = (m.N :* _cdid_nuisance_gram(m.or_shared.R, m.or_shared.pivot, (h_beta + eval_beta)'))'
    }
    if (m.or_shared.intercept) {
        r_b = m.or_blocks[b].recipe_id
        u_mu = _cdid_selected_loading(m, r_b, rowsF, m.or_blocks[b].m_col_alpha)
        xi_alpha_out = (alpha_score[, 1] + u_mu) / m.or_blocks[b].mass
        if (p_kept > 0) {
            xi_alpha_out = xi_alpha_out - xi_beta_out * m.or_blocks[b].mean_d'
        }
    }
}

/* Verify that a prepared cdid_selected_model satisfies conditional_mtar IF replay invariants. */
void _cdid_cpt_sel_check_model(struct cdid_selected_model scalar m)
{
    if (!m.prepared)
        _cdid_cpt_sel_fail("owner", "Selected model is not prepared.", 498)
    if (m.population_ref == NULL)
        _cdid_cpt_sel_fail("owner", "Selected model must retain a population_ref pointer.", 498)
    if (m.ctl.points_only)
        _cdid_cpt_sel_fail("replay", "Selected model was prepared with points_only=1; IF replay requires complete nuisance state.", 459)
    if (m.ctl.selection != "conditional_mtar")
        _cdid_cpt_sel_fail("request", "cdid_cpt_selected requires ctl.selection == 'conditional_mtar'.", 498)
    if (m.ctl.est_method != "ipw" & m.ctl.est_method != "or" & m.ctl.est_method != "dr")
        _cdid_cpt_sel_fail("request", "Unsupported selected est_method: " + m.ctl.est_method, 498)
    if (m.ctl.selection_model != "pair_logit" & m.ctl.selection_model != "independent_period_logit")
        _cdid_cpt_sel_fail("request", "Unsupported selection_model: " + m.ctl.selection_model, 498)
    if (m.ctl.selection_pooling != "separate" & m.ctl.selection_pooling != "within_cohort" &
        m.ctl.selection_pooling != "across_cohorts")
        _cdid_cpt_sel_fail("request", "Unsupported selection_pooling: " + m.ctl.selection_pooling, 498)
}

/* Check whether a comparison (g_pos, s_cal, t_cal, K) has positive actual-pair support
   in population P on both target cohort g_pos and every control cohort in K. */
real scalar _cdid_cpt_sel_raw_supported(struct cdid_population scalar P,
    struct cdid_selected_model scalar m, real scalar g_pos,
    real scalar s_cal, real scalar t_cal, real colvector K,
    string scalar reason_out)
{
    struct cdid_selected_rows scalar pr
    real colvector obs_c
    real scalar k
    reason_out = ""
    if (rows(K) < 1) {
        reason_out = "no_clean_control_cohort"
        return(0)
    }
    pr = _cdid_selected_pair_rows_pop(P, m.calendar_values[s_cal], m.calendar_values[t_cal], (1::P.N))
    if (rows(pr.unit) < 1) {
        reason_out = "empty_pair_support"
        return(0)
    }
    obs_c = m.cohort_idx[pr.unit]
    if (!any(obs_c :== g_pos)) {
        reason_out = "target_pair_empty"
        return(0)
    }
    for (k = 1; k <= rows(K); k++) {
        if (!any(obs_c :== K[k])) {
            reason_out = "control_pair_empty"
            return(0)
        }
    }
    return(1)
}

/* Validate the complete original comparison + both-auxiliary union on static cdid_selected_model
   and map each auxiliary row q in 1..Q to its link index in m.links (0 if unsupported).
   Never fits or mutates any nuisance. */
void _cdid_cpt_sel_union_static(struct cdid_cpt_selected scalar S,
    struct cdid_population scalar P, struct cdid_public_request scalar req,
    struct cdid_selected_model scalar m, real colvector aux_link,
    real colvector aux_sup, string colvector aux_why)
{
    struct cdid_catalog_opts scalar opts
    struct cdid_cat scalar C
    real matrix occ
    real colvector K, link_matched, hit
    real scalar i, j, q, ell, pi, M, Q, g_pos, s_cal, t_cal, sup, found, used
    string scalar why
    opts = cdid_catalog_default_opts()
    opts.mem_limit_bytes = 8 * S.workspace_limit
    opts.work_cap = S.work_limit
    occ = _cdid_cpt_base_occ(P, req.cohorts, opts)
    C = cdid_catalog_build_opts(req.calendar, req.cohorts, req.target_cohorts,
        req.anticipation, occ, req.control_group, req.post, opts)
    M = rows(C.r_g); Q = rows(req.aux)
    S.work_ops = _cdid_base_work(C.work + M + Q + m.n_links, S.work_limit, S.work_ops)
    used = _cdid_base_cells(4, M + Q + m.n_links + 1, S.workspace_limit, S.workspace_entries)

    aux_link = J(Q, 1, 0)
    aux_sup = J(Q, 1, 0)
    aux_why = J(Q, 1, "")
    link_matched = J(m.n_links, 1, 0)

    for (i = 1; i <= M + Q; i++) {
        if (i <= M) {
            g_pos = C.r_g[i]
            s_cal = C.r_s[i]
            t_cal = C.r_t[i]
            j = C.r_pool[i]
            K = (C.pool_size[j] > 0 ?
                C.pool_entries[C.pool_start[j]..C.pool_start[j] + C.pool_size[j] - 1] :
                J(0, 1, 0))
        }
        else {
            q = i - M
            hit = _cdid_sel_link_match_sorted(J(1, 1, req.aux[q, 1]), req.cohorts)
            g_pos = hit[1]
            hit = _cdid_sel_link_match_sorted(J(1, 1, req.aux[q, 2]), req.calendar)
            s_cal = hit[1]
            hit = _cdid_sel_link_match_sorted(J(1, 1, req.aux[q, 3]), req.calendar)
            t_cal = hit[1]
            if (g_pos < 1 | s_cal < 1 | t_cal < 1)
                _cdid_cpt_sel_fail("identity", "Auxiliary comparison endpoints are absent from request dictionaries.", 498)
            K = (req.aux_K_ptr[q + 1] > req.aux_K_ptr[q] ?
                req.aux_K_index[req.aux_K_ptr[q]..req.aux_K_ptr[q + 1] - 1] :
                J(0, 1, 0))
        }
        S.work_ops = _cdid_base_work(rows(K), S.work_limit, S.work_ops)
        found = 0
        if (m.n_links == M + Q) {
            pi = m.links[i].pair_id
            if (m.links[i].recipe_id == 1 & m.links[i].target_cohort_pos == g_pos &
                m.pairs[pi].s_cal == s_cal & m.pairs[pi].t_cal == t_cal &
                rows(m.links[i].control_cohort_pos) == rows(K)) {
                if (all(m.links[i].control_cohort_pos :== K)) {
                    found = i
                    link_matched[i] = 1
                }
            }
        }
        if (!found) {
            for (ell = 1; ell <= m.n_links; ell++) {
                pi = m.links[ell].pair_id
                if (m.links[ell].recipe_id == 1 & m.links[ell].target_cohort_pos == g_pos &
                    m.pairs[pi].s_cal == s_cal & m.pairs[pi].t_cal == t_cal &
                    rows(m.links[ell].control_cohort_pos) == rows(K)) {
                    if (all(m.links[ell].control_cohort_pos :== K)) {
                        if (!found) found = ell
                        link_matched[ell] = 1
                    }
                }
            }
        }
        why = ""
        sup = _cdid_cpt_sel_raw_supported(P, m, g_pos, s_cal, t_cal, K, why)
        if (sup & !found) {
            _cdid_cpt_sel_fail("identity",
                "Static selected model does not retain the complete original comparison and both-auxiliary union.", 498)
        }
        if (!sup & found) {
            _cdid_cpt_sel_fail("identity",
                "Static selected model contains an unsupported comparison link.", 498)
        }
        if (i > M) {
            q = i - M
            aux_link[q] = found
            aux_sup[q] = (sup & found > 0 & m.estimates[found] < .)
            if (!aux_sup[q]) {
                aux_why[q] = (why != "" ? why : "nonfinite_point")
            }
        }
    }
    if (any(link_matched :== 0)) {
        _cdid_cpt_sel_fail("identity",
            "Static selected model contains links outside the request's original comparison and auxiliary union.", 498)
    }
}

/* Validate the complete original comparison + both-auxiliary union on dated cdid_selref_provider.
   Never fits or mutates any nuisance. */
void _cdid_cpt_sel_union_dated(struct cdid_cpt_selected scalar S,
    struct cdid_population scalar P, struct cdid_public_request scalar req,
    struct cdid_reference_catalog scalar cat, struct cdid_reference_inputs scalar ref)
{
    struct cdid_catalog_opts scalar opts
    struct cdid_cat scalar C
    real matrix raw, occ
    real colvector K, actualK
    real scalar i, j, M, Q, c, used
    opts = cdid_catalog_default_opts()
    opts.mem_limit_bytes = 8 * S.workspace_limit
    opts.work_cap = S.work_limit
    occ = _cdid_cpt_base_occ(P, req.cohorts, opts)
    C = cdid_catalog_build_opts(req.calendar, req.cohorts, req.target_cohorts,
        req.anticipation, occ, req.control_group, req.post, opts)
    M = rows(C.r_g); Q = rows(req.aux)
    raw = (C.cohorts[C.r_g], C.calendar[C.r_s], C.calendar[C.r_t]) \ req.aux
    S.work_ops = _cdid_base_work(C.work + M + Q, S.work_limit, S.work_ops)
    used = _cdid_base_cells(4, M + Q + 1, S.workspace_limit, S.workspace_entries)
    if (cat.M != M | cat.Q != Q | cat.n_raw != M + Q | rows(ref.view_meta) != M + Q)
        _cdid_cpt_sel_fail("identity", "Dated selected provider does not retain the complete original comparison and auxiliary union.", 498)
    if (!all(raw :== ref.view_meta[, 1..3]))
        _cdid_cpt_sel_fail("identity", "Dated original comparison or auxiliary endpoints differ from the request.", 498)
    for (i = 1; i <= M + Q; i++) {
        if (i <= M) {
            j = C.r_pool[i]
            K = C.pool_entries[C.pool_start[j]..C.pool_start[j] + C.pool_size[j] - 1]
        }
        else {
            j = i - M
            K = req.aux_K_index[req.aux_K_ptr[j]..req.aux_K_ptr[j + 1] - 1]
        }
        S.work_ops = _cdid_base_work(rows(K), S.work_limit, S.work_ops)
        c = cat.raw_consumer[i]
        if (!c) continue
        actualK = cat.K_index[cat.K_ptr[c]..cat.K_ptr[c + 1] - 1]
        if (rows(K) != rows(actualK))
            _cdid_cpt_sel_fail("identity", "A complete-union control set differs from the request.", 498)
        if (!all(K :== actualK))
            _cdid_cpt_sel_fail("identity", "A complete-union control set differs from the request.", 498)
    }
}

/* Prepare selected conditional-pretrend piece producer from a static cdid_selected_model
   that was fitted on the full original post + both-auxiliary union. */
struct cdid_cpt_selected scalar cdid_cpt_selected_prepare_static(
    pointer(struct cdid_selected_model scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    struct cdid_cpt_selected scalar S
    pointer(struct cdid_selected_model scalar) scalar m
    pointer(struct cdid_population scalar) scalar P
    struct cdid_cpt_sel_cell_meta scalar cm
    pointer(struct cdid_sel_fit_result scalar) scalar fit
    pointer(struct cdid_sel_up_map scalar) scalar um
    real colvector aux_link, aux_sup, sel_rows, s_vec, D_vec, pc_ids
    string colvector aux_why, coords
    real scalar fam, i, r_sel, q, ell, tm, b, sid, p_or, icpt, is_sep, n_cells, eff_workcap, ws_lim, ch_rows
    if (provider == NULL)
        _cdid_cpt_sel_fail("owner", "Static selected model pointer is required.", 498)
    m = provider
    _cdid_cpt_sel_check_model(*m)
    P = (*m).population_ref
    if ((*m).R != 1 | rows((*m).contexts) != 1 | (*m).contexts[1].design_cal != 0 | (*m).contexts[1].weight_cal != 0)
        _cdid_cpt_sel_fail("identity", "Static selected model must have a single dateless design block.", 498)
    if (workcap <. & (workcap < 1 | workcap != floor(workcap)))
        _cdid_cpt_sel_fail("bound", "workcap must be . or a positive integer.", 498)
    eff_workcap = ((*m).ctl.max_work > 0 & (*m).ctl.max_work <. ?
        (workcap <. ? min((workcap, (*m).ctl.max_work)) : (*m).ctl.max_work) : workcap)
    ws_lim = ((*m).ctl.max_workspace_doubles > 0 & (*m).ctl.max_workspace_doubles <. ?
        (*m).ctl.max_workspace_doubles : .)
    ch_rows = ((*m).ctl.max_chunk_rows > 0 & (*m).ctl.max_chunk_rows <. ?
        min(((*m).N, (*m).ctl.max_chunk_rows)) : (*m).N)
    if (rows((*m).calendar_values) != rows(req.calendar) | !all((*m).calendar_values :== req.calendar) |
        rows((*m).cohort_values) != rows(req.cohorts) | !all((*m).cohort_values :== req.cohorts) |
        (*P).N != (*m).N | (*P).n_A != (*m).n_A)
        _cdid_cpt_sel_fail("identity", "Static selected model must be fitted on the full original population and dictionaries of the request.", 498)
    if (req.base_period != "varying" & req.base_period != "universal")
        _cdid_cpt_sel_fail("request", "Invalid request base_period.", 498)

    S.kind = "static"
    S.method = ((*m).ctl.est_method == "or" ? "reg" : (*m).ctl.est_method)
    S.selection_model = (*m).ctl.selection_model
    S.selection_pooling = (*m).ctl.selection_pooling
    S.selection_universe = (*m).contexts[1].universe
    S.base_period = req.base_period
    S.n_F = (*m).N
    S.n_A = (*m).n_A
    S.chunk_rows = ch_rows
    S.workspace_entries = (*m).workspace_doubles
    S.workspace_limit = ws_lim
    S.work_limit = eff_workcap
    S.work_ops = _cdid_base_work(0, eff_workcap, (*m).work_ops)
    S.n_sel_fits = (*m).n_selection_fit_calls
    S.n_trt_fits = (*m).n_treatment_fit_calls
    S.n_or_fits = (*m).n_or_block_fits
    S.static_provider = provider
    S.dated_provider = NULL

    aux_link = J(0, 1, 0); aux_sup = J(0, 1, 0); aux_why = J(0, 1, "")
    _cdid_cpt_sel_union_static(S, *P, req, *m, aux_link, aux_sup, aux_why)

    fam = (req.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_cpt_sel_ix((req.selection[, 1] :== fam) :& (req.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    S.workspace_entries = _cdid_base_cells(24, n_cells + 1, S.workspace_limit, S.workspace_entries)
    S.work_ops = _cdid_base_work(rows(req.selection) + n_cells, S.work_limit, S.work_ops)
    S.n_cells = n_cells
    S.n_unavailable = 0
    S.cells = cdid_cpt_sel_cell_meta(n_cells, 1)
    if (n_cells == 0) {
        S.status = "unavailable"
        S.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + req.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        return(S)
    }

    is_sep = ((*m).outcome_pooling == "separate")
    p_or = (S.method == "ipw" ? 0 : (*m).or_shared.p_kept)
    icpt = (S.method == "ipw" ? 0 : (*m).or_shared.intercept)
    fit = &((*m).contexts[1].sel_fit)
    um = &((*m).contexts[1].up_map)

    for (i = 1; i <= n_cells; i++) {
        r_sel = sel_rows[i]
        cm.g = req.selection[r_sel, 2]
        cm.t = req.selection[r_sel, 3]
        cm.from = req.selection[r_sel, 5]
        cm.q = req.selection[r_sel, 6]
        cm.sign = req.selection[r_sel, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.s_index = 0
        cm.t2_index = 0
        cm.link_idx = 0
        cm.recipe = 0
        cm.design_index = 0
        cm.supported = 0
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = ch_rows
        cm.spec = cdid_cpt_spec()
        cm.n_up = 0
        cm.p_trt = 0
        cm.p_or = 0
        cm.or_icpt = 0
        cm.up_channels = J(0, 1, 0)
        cm.up_blocks = J(0, 1, 0)
        cm.route = "static selected " + S.method + " (" + S.selection_model + ", " +
            S.selection_pooling + ", universe=" + S.selection_universe + ")"
        if (req.selection[r_sel, 8] != 0) {
            cm.status = "unavailable"
            cm.reason = req.selection_reason[r_sel]
            S.n_unavailable = S.n_unavailable + 1
            S.cells[i] = cm
            continue
        }
        q = cm.q
        cm.s_date = req.aux[q, 2]
        cm.t2_date = req.aux[q, 3]
        if (!aux_sup[q] | aux_link[q] < 1) {
            cm.status = "unavailable"
            cm.reason = "unsupported auxiliary comparison: " + aux_why[q]
            S.n_unavailable = S.n_unavailable + 1
            S.cells[i] = cm
            continue
        }
        ell = aux_link[q]
        cm.link_idx = ell
        cm.recipe = 1
        cm.design_index = 0
        cm.s_index = (*m).pairs[(*m).links[ell].pair_id].s_cal
        cm.t2_index = (*m).pairs[(*m).links[ell].pair_id].t_cal
        cm.supported = 1
        cm.status = "available"
        cm.reason = ""
        cm.estimate = (*m).estimates[ell]
        cm.estimate_scalar = cm.sign * (*m).estimates[ell]
        cm.chunk_rows = min((S.n_F, ch_rows))
        cm.n_chunks = ceil(S.n_F / cm.chunk_rows)
        if (S.method == "reg") {
            sid = (*m).links[ell].target_support_id
            b = (*m).links[ell].block_id
            p_or = (is_sep ? (*m).or_blocks[b].p_kept : (*m).or_shared.p_kept)
            pc_ids = J(1, 1, (*m).supports[sid].pc_id)
            cm.up_channels = _cdid_cpt_sel_pc_channels(*fit, *um, pc_ids[1])
            cm.up_blocks = _cdid_cpt_sel_pc_blocks(*fit, *um, pc_ids)
            cm.n_up = rows(cm.up_channels)
            cm.p_trt = 0
            cm.p_or = p_or
            cm.or_icpt = icpt
            coords = _cdid_cpt_sel_coords_or(cm.up_channels, icpt, b, p_or, is_sep)
            s_vec = 1
            D_vec = (*m).links[ell].d_target
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        else {
            tm = (*m).links[ell].trt_id
            b = (*m).trt_models[tm].block_id
            p_or = (S.method == "dr" ? (is_sep ? (*m).or_blocks[b].p_kept : (*m).or_shared.p_kept) : 0)
            pc_ids = (*m).trt_models[tm].pc_ids
            cm.up_channels = _cdid_cpt_sel_trt_channels(*fit, *um, pc_ids)
            cm.up_blocks = _cdid_cpt_sel_pc_blocks(*fit, *um, pc_ids)
            cm.n_up = rows(cm.up_channels)
            cm.p_trt = (*m).trt_models[tm].p_trt
            cm.p_or = p_or
            cm.or_icpt = (S.method == "dr" ? icpt : 0)
            coords = _cdid_cpt_sel_coords_two_arm(S.method == "dr", cm.n_up, cm.p_trt, cm.or_icpt, b, cm.p_or, is_sep)
            s_vec = (1 \ -1)
            D_vec = ((*m).trt_models[tm].dp \ (*m).trt_models[tm].dm)
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        S.workspace_entries = _cdid_base_cells(1,
            cm.n_up + rows(cm.up_blocks) + 2 * cm.spec.R,
            S.workspace_limit, S.workspace_entries)
        S.cells[i] = cm
    }

    if (S.n_unavailable > 0) {
        S.status = "unavailable"
        S.reason = "The complete requested test is unavailable: " + strofreal(S.n_unavailable, "%21.0f") + " declared cell(s) are unavailable. No available-only subset is tested."
    }
    else {
        S.status = "available"
        S.reason = ""
    }
    return(S)
}

/* Prepare selected conditional-pretrend piece producer from a dated cdid_selref_provider
   that was fitted on the full original post + both-auxiliary union. */
struct cdid_cpt_selected scalar cdid_cpt_selected_prepare_dated(
    pointer(struct cdid_selref_provider scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    struct cdid_cpt_selected scalar S
    pointer(struct cdid_selref_provider scalar) scalar prov
    pointer(struct cdid_selected_model scalar) scalar m
    pointer(struct cdid_population scalar) scalar P
    pointer(struct cdid_reference_inputs scalar) scalar ref
    pointer(struct cdid_baseline_plan_state scalar) scalar plan
    pointer(struct cdid_reference_catalog scalar) scalar cat
    struct cdid_cpt_sel_cell_meta scalar cm
    pointer(struct cdid_sel_fit_result scalar) scalar fit
    pointer(struct cdid_sel_up_map scalar) scalar um
    real colvector sel_rows, s_vec, D_vec, pc_ids
    string colvector coords
    real scalar fam, i, r_sel, q, ell, rec, tm, b, sid, p_or, icpt, is_sep, n_cells, eff_workcap, ws_lim, ch_rows
    if (provider == NULL)
        _cdid_cpt_sel_fail("owner", "Dated selected provider pointer is required.", 498)
    prov = provider
    _cdid_selref_check(*prov)
    m = &((*prov).model)
    _cdid_cpt_sel_check_model(*m)
    P = (*prov).population
    ref = (*prov).reference
    plan = (*prov).plan
    cat = (*prov).catalog
    cdid_ref_activate(*ref)
    if (workcap <. & (workcap < 1 | workcap != floor(workcap)))
        _cdid_cpt_sel_fail("bound", "workcap must be . or a positive integer.", 498)
    eff_workcap = workcap
    if ((*m).ctl.max_work > 0 & (*m).ctl.max_work <.) eff_workcap = min((eff_workcap, (*m).ctl.max_work))
    ws_lim = ((*m).ctl.max_workspace_doubles > 0 & (*m).ctl.max_workspace_doubles <. ?
        min(((*plan).workspace_limit, (*m).ctl.max_workspace_doubles)) : (*plan).workspace_limit)
    ch_rows = min(((*prov).n_F, (*prov).chunk_rows))
    if ((*m).ctl.max_chunk_rows > 0 & (*m).ctl.max_chunk_rows <.) ch_rows = min((ch_rows, (*m).ctl.max_chunk_rows))
    if ((*cat).Q != rows(req.aux) |
        rows((*plan).calendar) != rows(req.calendar) | !all((*plan).calendar :== req.calendar) |
        rows((*plan).cohorts) != rows(req.cohorts) | !all((*plan).cohorts :== req.cohorts) |
        (*P).N != (*prov).n_F | (*P).n_A != (*prov).n_A | (*m).N != (*prov).n_F | (*m).n_A != (*prov).n_A)
        _cdid_cpt_sel_fail("identity", "Dated selected provider must be fitted on the full original post + both-auxiliary union of the request.", 498)
    if (req.base_period != "varying" & req.base_period != "universal")
        _cdid_cpt_sel_fail("request", "Invalid request base_period.", 498)

    S.kind = "dated"
    S.method = (*prov).method
    S.selection_model = (*m).ctl.selection_model
    S.selection_pooling = (*m).ctl.selection_pooling
    S.selection_universe = (rows((*m).contexts) >= 1 ? (*m).contexts[1].universe : "F")
    S.base_period = req.base_period
    S.n_F = (*prov).n_F
    S.n_A = (*prov).n_A
    S.chunk_rows = ch_rows
    S.workspace_entries = (*prov).workspace_entries
    S.workspace_limit = ws_lim
    S.work_limit = eff_workcap
    S.work_ops = _cdid_base_work(0, eff_workcap, (*plan).work_ops + (*m).work_ops)
    S.n_sel_fits = (*m).n_selection_fit_calls
    S.n_trt_fits = (*m).n_treatment_fit_calls
    S.n_or_fits = (*m).n_or_block_fits
    S.static_provider = NULL
    S.dated_provider = provider

    _cdid_cpt_sel_union_dated(S, *P, req, *cat, *ref)

    fam = (req.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_cpt_sel_ix((req.selection[, 1] :== fam) :& (req.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    S.workspace_entries = _cdid_base_cells(24, n_cells + 1, S.workspace_limit, S.workspace_entries)
    S.work_ops = _cdid_base_work(rows(req.selection) + n_cells, S.work_limit, S.work_ops)
    S.n_cells = n_cells
    S.n_unavailable = 0
    S.cells = cdid_cpt_sel_cell_meta(n_cells, 1)
    if (n_cells == 0) {
        S.status = "unavailable"
        S.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + req.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        return(S)
    }

    is_sep = ((*m).outcome_pooling == "separate")
    p_or = (S.method == "ipw" ? 0 : (*m).or_shared.p_kept)
    icpt = (S.method == "ipw" ? 0 : (*m).or_shared.intercept)

    for (i = 1; i <= n_cells; i++) {
        r_sel = sel_rows[i]
        cm.g = req.selection[r_sel, 2]
        cm.t = req.selection[r_sel, 3]
        cm.from = req.selection[r_sel, 5]
        cm.q = req.selection[r_sel, 6]
        cm.sign = req.selection[r_sel, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.s_index = 0
        cm.t2_index = 0
        cm.link_idx = 0
        cm.recipe = 0
        cm.design_index = 0
        cm.supported = 0
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = ch_rows
        cm.spec = cdid_cpt_spec()
        cm.n_up = 0
        cm.p_trt = 0
        cm.p_or = 0
        cm.or_icpt = 0
        cm.up_channels = J(0, 1, 0)
        cm.up_blocks = J(0, 1, 0)
        cm.route = "dated D60 selected " + S.method + " (" + S.selection_model + ", " +
            S.selection_pooling + ", universe=" + S.selection_universe + ", reference-date recipes)"
        if (req.selection[r_sel, 8] != 0) {
            cm.status = "unavailable"
            cm.reason = req.selection_reason[r_sel]
            S.n_unavailable = S.n_unavailable + 1
            S.cells[i] = cm
            continue
        }
        q = cm.q
        cm.s_date = req.aux[q, 2]
        cm.t2_date = req.aux[q, 3]
        if (!(*prov).aux_supported[q] | (*prov).aux_consumer[q] < 1) {
            cm.status = "unavailable"
            cm.reason = "unsupported auxiliary comparison: " + (*prov).aux_reason[q]
            S.n_unavailable = S.n_unavailable + 1
            S.cells[i] = cm
            continue
        }
        ell = (*prov).aux_consumer[q]
        rec = (*prov).consumer_recipe[ell]
        cm.link_idx = ell
        cm.recipe = rec
        cm.design_index = (*prov).recipe_dates[rec, 1]
        cm.s_index = (*m).pairs[(*m).links[ell].pair_id].s_cal
        cm.t2_index = (*m).pairs[(*m).links[ell].pair_id].t_cal
        cm.supported = 1
        cm.status = "available"
        cm.reason = ""
        cm.estimate = (*prov).points[ell]
        cm.estimate_scalar = cm.sign * (*prov).points[ell]
        cm.chunk_rows = min((S.n_F, ch_rows))
        cm.n_chunks = ceil(S.n_F / cm.chunk_rows)
        fit = &((*m).contexts[rec].sel_fit)
        um = &((*m).contexts[rec].up_map)
        if (S.method == "reg") {
            sid = (*m).links[ell].target_support_id
            b = (*m).links[ell].block_id
            p_or = (is_sep ? (*m).or_blocks[b].p_kept : (*m).or_shared.p_kept)
            pc_ids = J(1, 1, (*m).supports[sid].pc_id)
            cm.up_channels = _cdid_cpt_sel_pc_channels(*fit, *um, pc_ids[1])
            cm.up_blocks = _cdid_cpt_sel_pc_blocks(*fit, *um, pc_ids)
            cm.n_up = rows(cm.up_channels)
            cm.p_trt = 0
            cm.p_or = p_or
            cm.or_icpt = icpt
            coords = _cdid_cpt_sel_coords_or(cm.up_channels, icpt, b, p_or, is_sep)
            s_vec = 1
            D_vec = (*m).links[ell].d_target
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        else {
            tm = (*m).links[ell].trt_id
            b = (*m).trt_models[tm].block_id
            p_or = (S.method == "dr" ? (is_sep ? (*m).or_blocks[b].p_kept : (*m).or_shared.p_kept) : 0)
            pc_ids = (*m).trt_models[tm].pc_ids
            cm.up_channels = _cdid_cpt_sel_trt_channels(*fit, *um, pc_ids)
            cm.up_blocks = _cdid_cpt_sel_pc_blocks(*fit, *um, pc_ids)
            cm.n_up = rows(cm.up_channels)
            cm.p_trt = (*m).trt_models[tm].p_trt
            cm.p_or = p_or
            cm.or_icpt = (S.method == "dr" ? icpt : 0)
            coords = _cdid_cpt_sel_coords_two_arm(S.method == "dr", cm.n_up, cm.p_trt, cm.or_icpt, b, cm.p_or, is_sep)
            s_vec = (1 \ -1)
            D_vec = ((*m).trt_models[tm].dp \ (*m).trt_models[tm].dm)
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        S.workspace_entries = _cdid_base_cells(1,
            cm.n_up + rows(cm.up_blocks) + 2 * cm.spec.R,
            S.workspace_limit, S.workspace_entries)
        S.cells[i] = cm
    }

    if (S.n_unavailable > 0) {
        S.status = "unavailable"
        S.reason = "The complete requested test is unavailable: " + strofreal(S.n_unavailable, "%21.0f") + " declared cell(s) are unavailable. No available-only subset is tested."
    }
    else {
        S.status = "available"
        S.reason = ""
    }
    return(S)
}

/* Short aliases for prepare entry points. */
struct cdid_cpt_selected scalar cdid_cpt_sel_prepare_static(
    pointer(struct cdid_selected_model scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    return(cdid_cpt_selected_prepare_static(provider, req, workcap))
}

struct cdid_cpt_selected scalar cdid_cpt_sel_prepare_dated(
    pointer(struct cdid_selref_provider scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    return(cdid_cpt_selected_prepare_dated(provider, req, workcap))
}

struct cdid_cpt_spec scalar cdid_cpt_selected_cell_spec(
    struct cdid_cpt_selected scalar S, real scalar cell_idx)
{
    if (cell_idx < 1 | cell_idx > S.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_sel_fail("cell", "Invalid declared cell index.", 498)
    if (!S.cells[cell_idx].supported)
        _cdid_cpt_sel_fail("cell", "Requested cell is unavailable: " + S.cells[cell_idx].reason, 459)
    return(S.cells[cell_idx].spec)
}

real scalar cdid_cpt_selected_cell_chunks(
    struct cdid_cpt_selected scalar S, real scalar cell_idx)
{
    if (cell_idx < 1 | cell_idx > S.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_sel_fail("cell", "Invalid declared cell index.", 498)
    if (!S.cells[cell_idx].supported)
        _cdid_cpt_sel_fail("cell", "Requested cell is unavailable: " + S.cells[cell_idx].reason, 459)
    return(S.cells[cell_idx].n_chunks)
}

/* Internal piece evaluator on arbitrary unique original-F rows over cdid_selected_model
   (shared by static and recipe-indexed dated models). */
struct cdid_cpt_piece scalar _cdid_cpt_sel_piece_model_rows(
    struct cdid_cpt_selected scalar S,
    struct cdid_selected_model scalar m,
    struct cdid_cpt_sel_cell_meta scalar cm,
    real colvector rowsF,
    | transmorphic scalar cache)
{
    struct cdid_selected_rows scalar rd
    pointer(struct cdid_sel_fit_result scalar) scalar fit
    pointer(struct cdid_sel_up_map scalar) scalar um
    struct cdid_sel_log_deriv scalar ld
    real scalar n, ell, r, sid, tm, b, L, n_up, p_trt, p_or, icpt, Dl, n_obs, k, pc, is_sep, use_cache
    real rowvector kept_b
    real colvector pos, rows_rr, dy, q_w, kk, D_ind, odds, act, rho_obs, sub, ch_k, col_idx, beta_b
    real matrix xi_alpha, xi_up, xi_ps, xi_beta, xi, Z_or, d_or, nu, rh, lam, pim, X_trt, lv, lam_up
    string scalar tmkey, akey, bkey
    use_cache = (args() >= 5)
    n = rows(rowsF)
    ell = cm.link_idx
    r = cm.recipe
    n_up = cm.n_up
    p_trt = cm.p_trt
    p_or = cm.p_or
    icpt = cm.or_icpt
    Dl = cm.spec.Dl
    is_sep = (m.outcome_pooling == "separate")

    _cdid_selected_activate_reader(m, r)
    fit = &(m.contexts[r].sel_fit)
    um = &(m.contexts[r].up_map)
    xi_up = _cdid_cpt_sel_up_if(m, r, cm.up_channels, cm.up_blocks, rowsF)

    if (S.method == "reg") {
        sid = m.links[ell].target_support_id
        b = m.links[ell].block_id
        L = m.links[ell].L
        pc = m.supports[sid].pc_id
        kept_b = (is_sep ? m.or_blocks[b].kept : m.or_shared.kept)
        beta_b = (is_sep ? m.or_blocks[b].beta : m.or_shared.beta)
        xi_alpha = J(n, 0, 0); xi_beta = J(n, 0, 0)
        if (icpt | p_or > 0) {
            if (use_cache) {
                akey = "oa:" + strofreal(b); bkey = "ob:" + strofreal(b)
                if (asarray_contains(cache, akey)) {
                    xi_alpha = asarray(cache, akey); xi_beta = asarray(cache, bkey)
                }
                else {
                    _cdid_cpt_sel_or_if(m, b, rowsF, xi_alpha, xi_beta)
                    asarray(cache, akey, xi_alpha); asarray(cache, bkey, xi_beta)
                }
            }
            else _cdid_cpt_sel_or_if(m, b, rowsF, xi_alpha, xi_beta)
        }
        xi = xi_up
        if (icpt) xi = (xi, xi_alpha)
        if (p_or > 0) xi = (xi, xi_beta)

        rd = _cdid_selected_support_rows(m, sid, rowsF)
        n_obs = rows(rd.loc)
        nu = J(n, 1, 0); rh = J(n, 1, 0)
        lam = J(n, Dl, 0); pim = J(n, Dl, 0)
        if (n_obs > 0) {
            pos = rd.loc; rows_rr = rd.unit; dy = rd.dy; q_w = rd.q
            _cdid_selected_activate_reader(m, r)
            ld = cdid_sel_log_deriv_eval(*fit, m.reader, *um, pc, rows_rr)
            d_or = (p_or > 0 ? L :* _cdid_sld_or_rows(m, r, rows_rr, kept_b) : J(n_obs, 0, 0))
            rho_obs = dy :- (icpt ? m.or_blocks[b].alpha : 0) :-
                (p_or > 0 ? d_or * beta_b : J(n_obs, 1, 0))
            nu[pos, 1] = q_w
            rh[pos, 1] = rho_obs
            if (ld.p_pc > 0) lam[pos, 1..n_up] = (ld.d_ak, ld.d_theta_local)
            else lam[pos, 1..n_up] = ld.d_ak
            if (icpt) pim[pos, n_up + 1] = J(n_obs, 1, -1)
            if (p_or > 0) pim[pos, (n_up + icpt) :+ (1..p_or)] = -d_or
        }
        return(cdid_cpt_make_piece(rowsF, nu, rh, lam, J(0, 0, 0), pim, J(0, 0, 0),
            1 / cm.spec.D, xi, cm.spec))
    }

    tm = m.links[ell].trt_id
    b = m.trt_models[tm].block_id
    L = (S.method == "dr" ? m.or_blocks[b].L : 0)
    if (use_cache) {
        tmkey = "ps:" + strofreal(tm)
        if (asarray_contains(cache, tmkey)) xi_ps = asarray(cache, tmkey)
        else {
            xi_ps = cdid_selected_ps_eval_chunk(m, tm, rowsF)
            asarray(cache, tmkey, xi_ps)
        }
    }
    else xi_ps = cdid_selected_ps_eval_chunk(m, tm, rowsF)
    xi = (xi_up, xi_ps)
    if (S.method == "dr") {
        xi_alpha = J(n, 0, 0); xi_beta = J(n, 0, 0)
        if (icpt | p_or > 0) {
            if (use_cache) {
                akey = "oa:" + strofreal(b); bkey = "ob:" + strofreal(b)
                if (asarray_contains(cache, akey)) {
                    xi_alpha = asarray(cache, akey); xi_beta = asarray(cache, bkey)
                }
                else {
                    _cdid_cpt_sel_or_if(m, b, rowsF, xi_alpha, xi_beta)
                    asarray(cache, akey, xi_alpha); asarray(cache, bkey, xi_beta)
                }
            }
            else _cdid_cpt_sel_or_if(m, b, rowsF, xi_alpha, xi_beta)
        }
        if (icpt) xi = (xi, xi_alpha)
        if (p_or > 0) xi = (xi, xi_beta)
    }

    rd = _cdid_selected_treatment_rows(m, tm, rowsF, 0)
    n_obs = rows(rd.loc)
    pos = J(0, 1, 0); D_ind = J(0, 1, 0); q_w = J(0, 1, 0)
    X_trt = J(0, p_trt, 0); odds = J(0, 1, 0); act = J(0, 1, 0)
    rho_obs = J(0, 1, 0); lam_up = J(0, n_up, 0); d_or = J(0, p_or, 0)
    if (n_obs > 0) {
        pos = rd.loc; rows_rr = rd.unit; dy = rd.dy; kk = rd.cohort_local; q_w = rd.q
        D_ind = (kk :== 1)
        _cdid_selected_activate_reader(m, r)
        lam_up = J(n_obs, n_up, 0)
        for (k = 1; k <= rows(m.trt_models[tm].cohorts); k++) {
            sub = _cdid_cpt_sel_ix(kk :== k)
            if (rows(sub) == 0) continue
            pc = m.trt_models[tm].pc_ids[k]
            ld = cdid_sel_log_deriv_eval(*fit, m.reader, *um, pc, rows_rr[sub])
            ch_k = _cdid_cpt_sel_pc_channels(*fit, *um, pc)
            col_idx = _cdid_cpt_sel_match(ch_k, cm.up_channels)
            if (ld.p_pc > 0) lam_up[sub, col_idx'] = (ld.d_ak, ld.d_theta_local)
            else lam_up[sub, col_idx'] = ld.d_ak
        }
        X_trt = _cdid_sld_trt_rows(m, r, rows_rr)
        lv = _cdid_ipw_logit_values(X_trt * m.trt_models[tm].beta_trt)
        odds = lv[, 2]
        act = lv[, 3]
        if (S.method == "dr") {
            kept_b = (is_sep ? m.or_blocks[b].kept : m.or_shared.kept)
            beta_b = (is_sep ? m.or_blocks[b].beta : m.or_shared.beta)
            Z_or = (p_or > 0 ? _cdid_sld_or_rows(m, r, rows_rr, kept_b) : J(n_obs, 0, 0))
            d_or = L :* Z_or
            rho_obs = dy :- (icpt ? m.or_blocks[b].alpha : 0) :-
                (p_or > 0 ? d_or * beta_b : J(n_obs, 1, 0))
        }
        else {
            rho_obs = dy
        }
    }
    return(cdid_cpt_two_arm_piece(rowsF, pos, D_ind, q_w, X_trt, odds, act, rho_obs,
        lam_up, xi, m.trt_models[tm].dp, m.trt_models[tm].dm, n_up, p_trt, d_or, icpt, cm.spec))
}

/* Evaluate one bounded piece on caller-supplied unique original-F row positions. */
struct cdid_cpt_piece scalar cdid_cpt_selected_piece_rows(
    struct cdid_cpt_selected scalar S, real scalar cell_idx, real colvector rowsF)
{
    pointer(struct cdid_selected_model scalar) scalar ms
    pointer(struct cdid_selref_provider scalar) scalar fd
    struct cdid_cpt_sel_cell_meta scalar cm
    real scalar n, R, Dl, used, is_inc
    if (cell_idx < 1 | cell_idx > S.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_sel_fail("replay", "Invalid declared cell index.", 498)
    cm = S.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_sel_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    n = rows(rowsF)
    is_inc = (n <= 1 ? 1 : all(rowsF[(2..n)] :> rowsF[(1..(n - 1))]))
    if (n < 1 | cols(rowsF) != 1 | missing(rowsF) | any(rowsF :!= floor(rowsF)) |
        any(rowsF :< 1) | any(rowsF :> S.n_F) | (!is_inc & rows(uniqrows(sort(rowsF, 1))) != n) |
        n > S.chunk_rows)
        _cdid_cpt_sel_fail("replay", "rowsF must be 1..chunk_rows unique original-F positions in 1..N_F.", 498)
    R = cm.spec.R; Dl = cm.spec.Dl
    used = _cdid_base_cells(12 * n, 2 * R + 3 * Dl + 8, S.workspace_limit, S.workspace_entries)
    if (S.kind == "static") {
        if (S.static_provider == NULL)
            _cdid_cpt_sel_fail("owner", "Borrowed static selected model is absent.", 498)
        ms = S.static_provider
        _cdid_cpt_sel_check_model(*ms)
        if ((*ms).N != S.n_F | (*ms).n_A != S.n_A |
            (*ms).n_selection_fit_calls != S.n_sel_fits |
            (*ms).n_treatment_fit_calls != S.n_trt_fits |
            (*ms).n_or_block_fits != S.n_or_fits)
            _cdid_cpt_sel_fail("owner", "Borrowed static selected model dimensions or fit counters changed.", 498)
        return(_cdid_cpt_sel_piece_model_rows(S, *ms, cm, rowsF))
    }
    if (S.dated_provider == NULL)
        _cdid_cpt_sel_fail("owner", "Borrowed dated selected provider is absent.", 498)
    fd = S.dated_provider
    _cdid_selref_check(*fd)
    if ((*fd).n_F != S.n_F | (*fd).n_A != S.n_A |
        (*fd).model.n_selection_fit_calls != S.n_sel_fits |
        (*fd).model.n_treatment_fit_calls != S.n_trt_fits |
        (*fd).model.n_or_block_fits != S.n_or_fits)
        _cdid_cpt_sel_fail("owner", "Borrowed dated selected provider dimensions or fit counters changed.", 498)
    return(_cdid_cpt_sel_piece_model_rows(S, (*fd).model, cm, rowsF))
}

/* Replay chunk k (1..n_chunks) on original-F partition 1..N_F. */
struct cdid_cpt_piece scalar cdid_cpt_selected_piece_chunk(
    struct cdid_cpt_selected scalar S, real scalar cell_idx, real scalar chunk_k)
{
    struct cdid_cpt_sel_cell_meta scalar cm
    real scalar lo, hi
    real colvector rowsF
    if (cell_idx < 1 | cell_idx > S.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_sel_fail("replay", "Invalid declared cell index.", 498)
    cm = S.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_sel_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    if (chunk_k < 1 | chunk_k > cm.n_chunks | chunk_k != floor(chunk_k))
        _cdid_cpt_sel_fail("replay", "Invalid piece chunk index.", 498)
    lo = (chunk_k - 1) * cm.chunk_rows + 1
    hi = min((S.n_F, lo + cm.chunk_rows - 1))
    rowsF = (lo..hi)'
    return(cdid_cpt_selected_piece_rows(S, cell_idx, rowsF))
}

/* Stream all chunks of one cell into the cdid_conditional engine accumulator. */
struct cdid_cpt_cell scalar cdid_cpt_selected_build_cell(
    struct cdid_cpt_selected scalar S, real scalar cell_idx,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar k, n_chunks
    spec = cdid_cpt_selected_cell_spec(S, cell_idx)
    n_chunks = cdid_cpt_selected_cell_chunks(S, cell_idx)
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, S.n_F, keep_rows, budget)
    for (k = 1; k <= n_chunks; k++) {
        pc = cdid_cpt_selected_piece_chunk(S, cell_idx, k)
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Stream an arbitrary partition/permutation row_order of original-F positions in batches of
   piece_size into the cdid_conditional engine accumulator (for split/reorder verification). */
struct cdid_cpt_cell scalar cdid_cpt_sel_build_cell_rows(
    struct cdid_cpt_selected scalar S, real scalar cell_idx,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget,
    real colvector row_order, real scalar piece_size)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar st, en, n_rows, bsz
    spec = cdid_cpt_selected_cell_spec(S, cell_idx)
    n_rows = rows(row_order)
    if (n_rows != S.n_F | cols(row_order) != 1 | missing(row_order) |
        any(row_order :!= floor(row_order)) | any(row_order :< 1) | any(row_order :> S.n_F) |
        rows(uniqrows(sort(row_order, 1))) != S.n_F | piece_size < 1 | piece_size != floor(piece_size) |
        piece_size > S.chunk_rows)
        _cdid_cpt_sel_fail("replay", "Invalid custom row_order or piece_size.", 498)
    bsz = min((piece_size, budget.max_piece_rows))
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, S.n_F, keep_rows, budget)
    for (st = 1; st <= n_rows; st = st + bsz) {
        en = min((n_rows, st + bsz - 1))
        pc = cdid_cpt_selected_piece_rows(S, cell_idx, row_order[(st..en)'])
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Build all declared cells 1..S.n_cells from caller-supplied grids (either 1 shared grid or n_cells grids). */
struct cdid_cpt_cell colvector cdid_cpt_selected_build_all(
    struct cdid_cpt_selected scalar S, struct cdid_cpt_grid colvector grids,
    real colvector cmap, real scalar C, real scalar keep_rows,
    struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_cell colvector out
    struct cdid_cpt_acc colvector accs
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_sel_cell_meta scalar cm
    pointer(struct cdid_selected_model scalar) scalar ms
    pointer(struct cdid_selref_provider scalar) scalar fd
    transmorphic scalar cache
    real scalar i, k, jg, max_chunks, lo, hi, n, R, Dl, used
    real colvector rowsF
    if (S.status != "available" | S.n_cells < 1 | S.n_unavailable > 0)
        _cdid_cpt_sel_fail("unavailable", S.reason, 459)
    if (rows(grids) != 1 & rows(grids) != S.n_cells)
        _cdid_cpt_sel_fail("grid", "Supply either one common grid or one grid per declared cell.", 498)
    if (S.kind == "static") {
        if (S.static_provider == NULL)
            _cdid_cpt_sel_fail("owner", "Borrowed static selected model is absent.", 498)
        ms = S.static_provider
        _cdid_cpt_sel_check_model(*ms)
        if ((*ms).N != S.n_F | (*ms).n_A != S.n_A |
            (*ms).n_selection_fit_calls != S.n_sel_fits |
            (*ms).n_treatment_fit_calls != S.n_trt_fits |
            (*ms).n_or_block_fits != S.n_or_fits)
            _cdid_cpt_sel_fail("owner", "Borrowed static selected model dimensions or fit counters changed.", 498)
    }
    else {
        if (S.dated_provider == NULL)
            _cdid_cpt_sel_fail("owner", "Borrowed dated selected provider is absent.", 498)
        fd = S.dated_provider
        _cdid_selref_check(*fd)
        if ((*fd).n_F != S.n_F | (*fd).n_A != S.n_A |
            (*fd).model.n_selection_fit_calls != S.n_sel_fits |
            (*fd).model.n_treatment_fit_calls != S.n_trt_fits |
            (*fd).model.n_or_block_fits != S.n_or_fits)
            _cdid_cpt_sel_fail("owner", "Borrowed dated selected provider dimensions or fit counters changed.", 498)
        ms = &((*fd).model)
    }
    out = cdid_cpt_cell(S.n_cells, 1)
    accs = cdid_cpt_acc(S.n_cells, 1)
    max_chunks = 0
    for (i = 1; i <= S.n_cells; i++) {
        jg = (rows(grids) == 1 ? 1 : i)
        accs[i] = cdid_cpt_cell_start(cdid_cpt_selected_cell_spec(S, i), grids[jg], cmap, C, S.n_F, keep_rows, budget)
        if (S.cells[i].n_chunks > max_chunks) max_chunks = S.cells[i].n_chunks
    }
    for (k = 1; k <= max_chunks; k++) {
        cache = asarray_create("string", 1)
        for (i = 1; i <= S.n_cells; i++) {
            cm = S.cells[i]
            if (k > cm.n_chunks) continue
            lo = (k - 1) * cm.chunk_rows + 1
            hi = min((S.n_F, lo + cm.chunk_rows - 1))
            rowsF = (lo..hi)'
            n = rows(rowsF); R = cm.spec.R; Dl = cm.spec.Dl
            used = _cdid_base_cells(12 * n, 2 * R + 3 * Dl + 8, S.workspace_limit, S.workspace_entries)
            pc = _cdid_cpt_sel_piece_model_rows(S, *ms, cm, rowsF, cache)
            jg = (rows(grids) == 1 ? 1 : i)
            cdid_cpt_cell_add(accs[i], grids[jg], cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
        }
    }
    for (i = 1; i <= S.n_cells; i++) {
        jg = (rows(grids) == 1 ? 1 : i)
        out[i] = cdid_cpt_cell_finish(accs[i], grids[jg])
    }
    return(out)
}

/* Oriented scalar h=1 IF on frame scale (1/N_F normalization) from the underlying
   selected provider on caller-supplied rowsF (chunked by S.chunk_rows when rows(rowsF) > S.chunk_rows). */
real colvector cdid_cpt_sel_provider_if(
    struct cdid_cpt_selected scalar S, real scalar cell_idx, real colvector rowsF)
{
    struct cdid_cpt_sel_cell_meta scalar cm
    struct cdid_selected_chunk_out scalar ch
    real scalar n, st, en, bsz
    real colvector out, sub
    if (cell_idx < 1 | cell_idx > S.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_sel_fail("replay", "Invalid declared cell index.", 498)
    cm = S.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_sel_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    n = rows(rowsF)
    if (n < 1 | cols(rowsF) != 1 | missing(rowsF))
        _cdid_cpt_sel_fail("replay", "Invalid rowsF for provider IF replay.", 498)
    bsz = S.chunk_rows
    out = J(n, 1, 0)
    for (st = 1; st <= n; st = st + bsz) {
        en = min((n, st + bsz - 1))
        sub = rowsF[(st..en)']
        if (S.kind == "static") {
            if (S.static_provider == NULL)
                _cdid_cpt_sel_fail("owner", "Borrowed static selected model is absent.", 498)
            ch = cdid_selected_chunk(*S.static_provider, sub, J(1, 1, cm.link_idx), 0, 0)
            out[(st..en)'] = cm.sign :* ch.phi[, 1]
        }
        else {
            if (S.dated_provider == NULL)
                _cdid_cpt_sel_fail("owner", "Borrowed dated selected provider is absent.", 498)
            out[(st..en)'] = cm.sign :*
                cdid_selref_if_chunk(*S.dated_provider, sub, J(1, 1, cm.link_idx), 0)[, 1]
        }
    }
    return(out)
}

end
