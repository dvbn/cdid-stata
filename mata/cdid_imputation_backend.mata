*! Pure-imputation backend Stages 2-5 (static 1-recipe and dated multi-recipe); StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (), R/chained_imputation.R:
*!   imputation_build_catalog (OR block and treatment keys), imputation_fit_core Stages 2-5,
*!   imputation_if_needed_sets and imputation_evaluate_if_chunk_core; R/chained_or_pooled.R
*!   or_pooled_fit_blocks; R/chained_cp04_corollary1.R cp04_logit and cp04_solve;
*! R/chained_reference_inputs.R (dated recipe-indexed designs, weights and repeat rates).
*! Load cdid_population, cdid_nuisance, cdid_ipw_link, cdid_imputation_model (
*! Stage 1, unchanged), cdid_gmm, and cdid_imputation_dated before this file.
version 19.5


mata:
mata set matastrict on

/* One deduplicated computational link (source links_fitted record). Unavailable arms keep
   missing scalars and empty vectors. d_bar_plus_j is in the kept order of the target cell
   group; d_bar_minus_j_k is n_ctrl x p_raw in raw slope order (kept columns of each control
   cell group, zeros elsewhere). */
struct _cdid_impb_link {
    string scalar link_id
    real scalar link_index, pair_id, recipe, L, target_ci, target_cell, or_block, trt_id, n_ctrl
    real colvector ctrl_ci, ctrl_cell
    real scalar mu, mu_ipw, mu_reg, mu_dr
    real scalar d_plus, num_T, T_gl, num_Bg, B_gl, d_minus, num_C, C_ol, num_Bo, B_ol
    real scalar a_coeff_T_g, a_coeff_Bg_g
    real rowvector d_bar_plus_j, d_bar_plus_or, d_bar_minus_or, d_bar_diff_or
    real colvector a_coeff_C_k, a_coeff_Bo_k, m_minus_k
    real matrix d_bar_minus_j_k
    real rowvector denminus_gamma, g_C_gamma, g_Bo_gamma, g_dr_gamma
}

/* One shared counterfactual OR block (pair, recipe, ordered control cohorts) with upstream forcing.
   jf_alpha_beta is n_ctrl x p_raw and jf_beta_beta is p_or x (n_ctrl * p_raw), both in raw
   slope order of each control cell group (source j_force_alpha_beta, j_force_beta_beta). */
struct _cdid_impb_block {
    string scalar name
    real scalar block_id, pair_id, recipe, L, n_ctrl, weight_sum, mass, alpha, mean_y
    real colvector ctrl_ci, ctrl_cell, links
    real rowvector mean_d
    real colvector a_force_alpha, jf_alpha_alpha
    real matrix a_force_beta, jf_beta_alpha, jf_alpha_beta, jf_beta_beta
    real scalar or_p, or_rss, or_rcond, or_adj_ready, score_mean_a
    real rowvector or_keep, or_pivot
    real colvector or_beta, score_mean_b
    real matrix or_R, or_alias_map, or_alias_const
    struct cdid_imputation_adjoint scalar or_adj
}

/* One full-F treatment logit per source key (pair, recipe, target, ordered controls). R is the raw
   factor of sum_i v p (1 - p) X X' (no 1/N): gamma IF rows are N * Gram^-1 score, N once. */
struct _cdid_impb_trt {
    string scalar name
    real scalar trt_id, rep_tid, pair_id, recipe, target_ci, target_cell, n_coh, n_risk
    real scalar iterations, converged, deviance, score_max, rcond_R
    real scalar n_capped, n_capped_ctrl, p_min, p_max
    real colvector all_ci, links, gamma
    real matrix R, a_forcing
    real rowvector pivot
}

/* Per-control-cell synthetic OR support summary used only while fitting. */
struct _cdid_impb_cacc {
    real scalar used, n, wsum
    real rowvector anchor, constant
    real matrix F
}

/* Cached index bucket for fast cohort slicing without repeated N-row scans. */
struct _cdid_impb_ivec {
    real colvector rows
}

/* Per-link control accumulator used while streaming Stage 4 link summaries. */
struct _cdid_impb_lacc {
    real scalar s_wm, s_wmj, s_wmm, po, alpha_b
    real colvector s_wm_k, s_wodj_k, s_wod_k, s_wodm_k, obeta
    real rowvector s_wdjX, s_wdX, s_wdmX, s_wmdor, okeep, dbp_or
    real matrix s_wmdimp
}

/* Transient (recipe, cohort) and cell-prediction caches used only during backend prepare. */
struct _cdid_impb_rcc {
    real scalar ready, has_x
    real colvector wr, v
    real matrix Zr, X
}

struct _cdid_impb_cjc {
    real scalar ready
    real colvector jh
}

/* Prepared backend. model, xtrt (static), and rdr (dated) are BORROWED: the caller owns
   these objects and must keep them alive and unmodified for the whole backend lifetime;
   every replay entry re-checks the borrowed fingerprints. */
struct cdid_imp_backend {
    string scalar backend_id, method, pooling, outcome_pooling, law, selection, normalization
    pointer(struct cdid_imputation_model scalar) scalar model
    pointer(real matrix) scalar xtrt
    pointer(struct cdid_imp_dated_reader scalar) scalar rdr
    struct cdid_imputation_budget scalar budget
    struct _cdid_impb_ivec rowvector coh_u
    real scalar N, n_A, lambda, G, n_cells, n_groups, n_recipes, p_raw, intercept, q
    real scalar rank_tol, glm_eps, glm_maxit, link_batch, chunk_rows
    real scalar has_or, has_ps, stage1_only, uses_selection_models, is_dated
    real scalar n_links, n_blocks, n_trt
    real scalar n_imputation_fits, n_or_fits, n_ps_fits, n_obs_logits
    real scalar borrow_sum, retained_bytes, scratch_bytes, work_ops
    real colvector a_k, mass_F, cell_M_one, cell_Y_one
    real matrix a_k_mat, mass_F_mat, cell_M_z, cell_Y_z
    struct _cdid_impb_link rowvector links
    struct _cdid_impb_block rowvector blocks
    struct _cdid_impb_trt rowvector trts
    real scalar or_p, or_rss, or_rcond, or_rows_checked, or_threshold, or_adj_ready
    real rowvector or_keep, or_retained, or_proposed, or_dropped, or_candidates, or_pivot
    real rowvector or_ctl_scaled, or_trt_scaled
    real matrix or_R, or_alias_map, or_alias_const
    real colvector or_beta
    struct cdid_imputation_adjoint scalar or_adj
}

/* Replay result: B rows (explicit original-F positions) by L requested links. comp_* are
   B x 8L, link-major, in source channel order direct, a_direct, j_direct, treatment_direct,
   treatment_a, OR_direct, OR_a, OR_j. Unavailable methods stay missing. */
struct cdid_imp_backend_if {
    string scalar method, scale
    real colvector f_rows, links, id_num
    string colvector id_str
    real scalar normalization_count, N, n_A, lambda, detail
    real matrix phi, psi_ipw, psi_reg, psi_dr, comp_ipw, comp_reg, comp_dr
}

/* Cached representative propensity group for batched Stage 5 contraction. */
struct _cdid_impb_tp_grp {
    real scalar rtid
    real rowvector link_pos
    real matrix G_gamma
}

/* Pre-built Stage 5 batch plan (detail = 0) shared across all row chunks. */
struct _cdid_impb_plan {
    real scalar ready, nl, nA, tot_po, sep
    real colvector lk, abids, bbids, bids, boff, zids, zpos, rep_tids, link_tp, or_bids, or_bpos
    real matrix or_coh_mask, A_load, Alpha_load, Beta_load, Raa_load, Rba_load
    struct _cdid_impb_tp_grp rowvector tp_grps
    struct cdid_imputation_adjoint scalar adj
}

/* GMM weight-provider context; backend is borrowed with the same lifetime rule. */
struct cdid_imp_backend_wctx {
    pointer(struct cdid_imp_backend scalar) scalar backend
    real scalar n_moments, row_chunk
    real colvector need_links, col_pos, col_sign
    string colvector consumer_ids
    real scalar fast_ready
    struct _cdid_impb_plan scalar plan
}

/* Vectorized raw probability, capped odds, and strict derivative mask matching _cdid_ipw_logit_values. */
real matrix _cdid_impb_logit_vals(real colvector eta)
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

void _cdid_impb_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_impb_ix(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}

real rowvector _cdid_impb_seq(real scalar lo, real scalar hi)
{
    if (hi < lo) return(J(1, 0, 0))
    return(lo..hi)
}

real scalar _cdid_impb_pos(real vector v, real scalar x)
{
    real scalar i
    for (i=1; i<=length(v); i++) {
        if (v[i] == x) return(i)
    }
    return(0)
}

/* Row-reader boundary for Stages 2-5. Dispatches to borrowed dated reader in dated mode
   (enforcing recipe-indexed used-cohort support) or static Stage 1 arrays in static mode. */
real colvector _cdid_impb_w(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar r, real colvector fr,
    | real scalar trusted)
{
    real scalar tr
    tr = (args() >= 5 & trusted != 0 & trusted < . ? 1 : 0)
    if (b.is_dated) return(_cdid_impd_w(*(b.rdr), r, fr, tr))
    return(m.recipes[r].weight[fr])
}

real matrix _cdid_impb_z(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar r, real colvector fr, real rowvector cc,
    | real scalar trusted)
{
    real scalar tr
    if (cols(cc) == 0) return(J(rows(fr), 0, 0))
    tr = (args() >= 6 & trusted != 0 & trusted < . ? 1 : 0)
    if (b.is_dated) return(_cdid_impd_z(*(b.rdr), r, fr, cc, tr))
    return(m.recipes[r].Z[fr, cc])
}

real matrix _cdid_impb_x(struct cdid_imp_backend scalar b, real scalar r, real colvector fr,
    | real scalar trusted)
{
    pointer(real matrix) scalar px
    real scalar tr
    tr = (args() >= 4 & trusted != 0 & trusted < . ? 1 : 0)
    if (b.is_dated) return(_cdid_impd_x(*(b.rdr), r, fr, tr))
    px = b.xtrt
    return((*px)[fr, ])
}

real scalar _cdid_impb_ak(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar r, real scalar ci)
{
    if (b.is_dated) return(_cdid_impd_ak(*(b.rdr), r, ci))
    return(m.cohort_a_F[r, ci])
}

real colvector _cdid_impb_pred_cell(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar cid, real colvector fr)
{
    if (b.is_dated) return(cdid_imp_predict_cell_d(m, *(b.rdr), cid, fr))
    return(cdid_imp_predict_cell_chunk(m, cid, fr))
}

real colvector _cdid_impb_pred_z(struct cdid_imputation_model scalar m, real scalar cid,
    real scalar n_r, real matrix Zr)
{
    real scalar gid
    real colvector out
    gid = m.cells[cid].group_id
    out = J(n_r, 1, m.cells[cid].alpha)
    if (m.groups[gid].p > 0) {
        out = out + (m.cells[cid].L :* Zr[, m.groups[gid].kept]) * m.groups[gid].beta
    }
    return(out)
}

/* Rows of cohort ci read under recipe r (only eligible units when the recipe has a
   restricted reference eligibility context). */
real colvector _cdid_impb_crows(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar r, real scalar ci)
{
    real colvector rr
    rr = b.coh_u[ci].rows
    if (rows(m.elig_keep) == 0 | rows(rr) == 0) return(rr)
    return(select(rr, m.elig_keep[rr, r] :!= 0))
}

/* Positions of cohort k in the chunk fr whose units recipe r reads. */
real colvector _cdid_impb_cix(struct _cdid_impb_ivec rowvector c_ix,
    struct cdid_imputation_model scalar m, real colvector fr, real scalar r, real scalar k)
{
    real colvector ix
    ix = c_ix[k].rows
    if (rows(m.elig_keep) == 0 | rows(ix) == 0) return(ix)
    return(select(ix, m.elig_keep[fr[ix], r] :!= 0))
}

struct _cdid_impb_ivec rowvector _cdid_impb_bucket_coh(real colvector coh, real scalar G)
{
    struct _cdid_impb_ivec rowvector out
    real colvector uc
    real scalar B, c, j
    out = _cdid_impb_ivec(1, G)
    for (c=1; c<=G; c++) out[c].rows = J(0, 1, 0)
    B = rows(coh)
    if (B == 0) return(out)
    c = coh[1]
    if (all(coh :== c)) {
        out[c].rows = (1::B)
        return(out)
    }
    if (G <= 16) {
        for (c=1; c<=G; c++) out[c].rows = _cdid_impb_ix(coh :== c)
        return(out)
    }
    uc = uniqrows(coh)
    for (j=1; j<=rows(uc); j++) {
        c = uc[j]
        out[c].rows = _cdid_impb_ix(coh :== c)
    }
    return(out)
}

real matrix _cdid_impb_adj_chunk(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, struct cdid_imputation_adjoint scalar adj,
    real colvector fr)
{
    if (b.is_dated) return(cdid_imp_adjoint_chunk_d(m, *(b.rdr), adj, fr, "frame"))
    return(cdid_imputation_adjoint_chunk(m, adj, fr, "frame"))
}

real scalar _cdid_impb_N(struct cdid_imp_backend scalar b)
{
    return(b.N)
}

real scalar _cdid_impb_nA(struct cdid_imp_backend scalar b)
{
    return(b.n_A)
}

real scalar _cdid_impb_nlinks(struct cdid_imp_backend scalar b)
{
    return(b.n_links)
}

real scalar _cdid_impb_chunk(struct cdid_imp_backend scalar b)
{
    return(b.chunk_rows)
}

/* Identity fingerprint of the borrowed Stage 1 fit (exact, not a tolerance). */
real scalar _cdid_impb_borrow_sum(struct cdid_imputation_model scalar m)
{
    real scalar s, cid, gid
    s = m.N + 3 * m.n_cells + 7 * m.n_groups + 11 * m.p_raw + 13 * m.n_A + 17 * m.n_recipes
    for (cid=1; cid<=m.n_cells; cid++) s = s + cid * m.cells[cid].alpha
    for (gid=1; gid<=m.n_groups; gid++) {
        if (m.groups[gid].p > 0)
            s = s + (gid + 0.5) * quadsum(m.groups[gid].beta :* (1::m.groups[gid].p))
    }
    return(s)
}

void _cdid_impb_need_model(struct cdid_imp_backend scalar b)
{
    if (b.model == NULL)
        _cdid_impb_fail("Borrowed first-step fit pointer is NULL; its lifetime ended.", 498)
    if (b.is_dated & b.rdr == NULL)
        _cdid_impb_fail("Borrowed dated reader pointer is NULL; its lifetime ended.", 498)
}

void _cdid_impb_check(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m)
{
    pointer(real matrix) scalar px
    real scalar cur_sum
    if (b.is_dated) {
        if (b.rdr == NULL)
            _cdid_impb_fail("Borrowed dated reader pointer is NULL; its lifetime ended.", 498)
        _cdid_impd_check_reader(*(b.rdr))
        cur_sum = (*(b.rdr)).borrow_sum + _cdid_impb_borrow_sum(m)
    }
    else {
        cur_sum = _cdid_impb_borrow_sum(m)
    }
    if (m.N != b.N | m.n_A != b.n_A | m.n_cells != b.n_cells | m.n_groups != b.n_groups |
        m.p_raw != b.p_raw | m.n_recipes != b.n_recipes | cur_sum != b.borrow_sum)
        _cdid_impb_fail("Borrowed first-step fit changed after backend preparation; prepare the backend again.", 498)
    if (!b.is_dated & b.has_ps) {
        px = b.xtrt
        if (px == NULL)
            _cdid_impb_fail("Borrowed treatment design pointer is NULL; its lifetime ended.", 498)
        if (rows(*px) != b.N | cols(*px) != b.q)
            _cdid_impb_fail("Borrowed treatment design changed shape after backend preparation.", 498)
    }
}

/* Original-F domain record, inactive selection record and source support gates. */
void _cdid_impb_init(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m)
{
    pointer(real matrix) scalar px
    real scalar j, ci, r
    if (!b.is_dated & m.n_recipes != 1)
        _cdid_impb_fail("The static core accepts exactly one canonical design block; dated multi-block layouts require cdid_imp_backend_prepare_d.", 498)
    if (b.is_dated) {
        if (b.rdr == NULL)
            _cdid_impb_fail("Dated backend preparation requires a borrowed cdid_imp_dated_reader pointer.", 498)
        _cdid_impd_check_reader(*(b.rdr))
        if ((*(b.rdr)).N != m.N | (*(b.rdr)).n_A != m.n_A | (*(b.rdr)).G != m.n_cohorts |
            (*(b.rdr)).n_recipes != m.n_recipes | (*(b.rdr)).p_raw != m.p_raw |
            (*(b.rdr)).intercept != m.intercept)
            _cdid_impb_fail("Borrowed dated reader disagrees with first-step model dimensions or metadata.", 498)
    }
    if (m.N < 1 | m.n_cells < 1 | m.n_groups < 1 | rows(m.pop.cohort_known) != m.N |
        any(m.pop.cohort_known :!= 1))
        _cdid_impb_fail("Pure imputation requires a fitted first-step model on the full F roster with known cohorts.", 459)
    b.N = m.N
    b.n_A = m.n_A
    b.lambda = m.lambda
    b.G = m.n_cohorts
    b.coh_u = _cdid_impb_ivec(1, b.G)
    for (ci=1; ci<=b.G; ci++) b.coh_u[ci].rows = _cdid_impb_ix(m.unit_cohort_index :== ci)
    b.n_cells = m.n_cells
    b.n_groups = m.n_groups
    b.n_recipes = m.n_recipes
    b.p_raw = m.p_raw
    b.intercept = m.intercept
    b.rank_tol = m.rank_tol
    b.pooling = m.pooling
    b.chunk_rows = min((b.budget.max_chunk_rows, m.budget.max_chunk_rows))
    b.a_k = m.cohort_a_F[1, ]'
    b.mass_F = m.cohort_mass_F[1, ]'
    b.a_k_mat = m.cohort_a_F
    b.mass_F_mat = m.cohort_mass_F
    b.selection = "none"
    b.uses_selection_models = 0
    b.n_obs_logits = 0
    b.n_imputation_fits = m.n_groups
    b.n_or_fits = 0
    b.n_ps_fits = 0
    b.n_links = 0
    b.n_blocks = 0
    b.n_trt = 0
    b.or_p = 0
    b.or_rss = .
    b.or_rcond = .
    b.or_rows_checked = 0
    b.or_threshold = .
    b.or_adj_ready = 0
    b.or_keep = J(1, 0, 0)
    b.or_retained = J(1, 0, 0)
    b.or_proposed = J(1, 0, 0)
    b.or_dropped = J(1, 0, 0)
    b.or_candidates = J(1, 0, 0)
    b.or_pivot = J(1, 0, 0)
    b.or_ctl_scaled = J(1, 0, 0)
    b.or_trt_scaled = J(1, 0, 0)
    b.or_R = J(0, 0, 0)
    b.or_alias_map = J(0, 0, 0)
    b.or_alias_const = J(0, 0, 0)
    b.or_beta = J(0, 1, 0)
    b.law = ("IMPUTATION:" + b.method + ":" + b.pooling + ":" +
        (b.intercept ? "local-intercepts" : "no-intercept") +
        (b.is_dated ? ":actual-pair:W_g:full-F:recipe-cohort-layout" : ":actual-pair:fixed-w:full-F"))
    b.normalization = "frame IF on all original F rows; core = (n_A/N_F) frame; weight Gram centers all F core rows and divides by n_A"
    if (b.n_A < 1 | b.n_A > b.N)
        _cdid_impb_fail("Repeat-observer count n_A must lie in 1..N_F.", 459)
    if (b.lambda != b.n_A / b.N)
        _cdid_impb_fail("First-step lambda differs from n_A/N_F.", 504)
    if (b.is_dated) {
        for (r=1; r<=b.n_recipes; r++) {
            for (ci=1; ci<=b.G; ci++) {
                if (!(*(b.rdr)).used_cohorts[r, ci]) continue
                if (!(b.mass_F_mat[r, ci] > 0) | b.mass_F_mat[r, ci] >= . |
                    !(b.a_k_mat[r, ci] > 0) | b.a_k_mat[r, ci] >= .)
                    _cdid_impb_fail("Cohort " + m.cohort_labels[ci] + " in recipe " +
                        strofreal(r, "%21.0g") + " has zero full-roster or repeat-observation mass on F (a_k <= 0).", 459)
            }
        }
    }
    else {
        for (j=1; j<=rows(m.active_cohorts); j++) {
            ci = m.active_cohorts[j]
            if (!(b.mass_F[ci] > 0) | b.mass_F[ci] >= . | !(b.a_k[ci] > 0) | b.a_k[ci] >= .)
                _cdid_impb_fail("Cohort " + m.cohort_labels[ci] + " has zero full-roster or repeat-observation mass on F (a_k <= 0); synthetic imputation weights require positive repeat mass.", 459)
        }
    }
    b.q = 0
    if (b.has_ps) {
        if (b.is_dated) {
            if (!(*(b.rdr)).has_ps | (*(b.rdr)).q < 1)
                _cdid_impb_fail("ipw and dr require a dated reader prepared with has_ps == 1 and q >= 1.", 498)
            b.q = (*(b.rdr)).q
        }
        else {
            px = b.xtrt
            if (px == NULL)
                _cdid_impb_fail("ipw and dr require a borrowed full-F treatment design pointer.", 498)
            if (rows(*px) != b.N | cols(*px) < 1)
                _cdid_impb_fail("Treatment design needs one row per original F unit and at least one column.", 498)
            if (missing(*px))
                _cdid_impb_fail("Pure imputation requires a finite treatment design on every original F row.", 459)
            b.q = cols(*px)
        }
    }
}

/* Link intake against the Stage 1 registry, then source block and treatment keys:
   OR block = (pair_id [which carries recipe], ordered controls);
   treatment = (pair_id [which carries recipe], target, ordered controls). */
void _cdid_impb_registry(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real matrix link_spec, real matrix link_controls, string colvector link_ids)
{
    real scalar L, ell, s_idx, t_idx, rec, pid, i, k, ci, cell, mlen, nb, nt, found, bid, tid, raw
    real colvector cc, kord, cellv, link_raw, ordb, rank_b, tvec
    real rowvector krow, kr
    real matrix bkeys, tkeys

    L = rows(link_spec)
    if (L < 1 | (cols(link_spec) != 3 & cols(link_spec) != 4) | missing(link_spec) |
        cols(link_controls) != 2 | missing(link_controls))
        _cdid_impb_fail("Require a finite L x 3 (or L x 4) link_spec and a finite M x 2 link_controls.", 498)
    if (rows(link_ids) != 0 & rows(link_ids) != L)
        _cdid_impb_fail("link_ids must be empty or carry one opaque identifier per link.", 498)
    for (i=1; i<=rows(link_controls); i++) {
        if (link_controls[i, 1] < 1 | link_controls[i, 1] > L |
            link_controls[i, 1] != floor(link_controls[i, 1]))
            _cdid_impb_fail("link_controls references an unknown link index.", 498)
    }
    b.n_links = L
    b.links = _cdid_impb_link(1, L)
    mlen = 0
    for (ell=1; ell<=L; ell++) {
        s_idx = _cdid_imp_find_exact(m.pop.calendar, link_spec[ell, 1])
        t_idx = _cdid_imp_find_exact(m.pop.calendar, link_spec[ell, 2])
        rec = (cols(link_spec) == 4 ? link_spec[ell, 4] : 1)
        if (s_idx < 1 | t_idx <= s_idx)
            _cdid_impb_fail("Link dates must lie on the stored calendar with s before t.", 498)
        if (!b.is_dated & rec != 1)
            _cdid_impb_fail("The static core accepts design block 1 only; dated design blocks require cdid_imp_backend_prepare_d.", 498)
        if (rec < 1 | rec > m.n_recipes | rec != floor(rec))
            _cdid_impb_fail("Link design block index out of range 1..n_recipes.", 498)
        pid = 0
        for (i=1; i<=m.n_pairs; i++) {
            if (m.pairs[i].s_index == s_idx & m.pairs[i].t_index == t_idx & m.pairs[i].recipe == rec) {
                if (pid > 0) _cdid_impb_fail("Ambiguous first-step pair identity for a link.", 498)
                pid = i
            }
        }
        if (pid == 0) _cdid_impb_fail("Link pair is not a fitted first-step pair.", 498)
        ci = 0
        if (link_spec[ell, 3] > 0) ci = _cdid_imp_find_exact(m.cohorts, link_spec[ell, 3])
        if (ci < 1) _cdid_impb_fail("Link target must be a known treated cohort of the first-step fit.", 498)
        if (b.is_dated) {
            if (!(*(b.rdr)).used_cohorts[rec, ci])
                _cdid_impb_fail("Link target cohort is outside design block used_cohorts support.", 459)
        }
        cell = m.cell_lookup[ci, pid]
        if (cell < 1) _cdid_impb_fail("First-step fit lacks a link target cell.", 498)
        cc = J(0, 1, 0)
        for (i=1; i<=rows(link_controls); i++) {
            if (link_controls[i, 1] != ell) continue
            k = _cdid_imp_find_exact(m.cohorts, link_controls[i, 2])
            if (k < 1) _cdid_impb_fail("Link control cohort is not a known cohort of the first-step fit.", 498)
            if (k == ci) _cdid_impb_fail("Link control set contains its own target cohort.", 498)
            if (b.is_dated) {
                if (!(*(b.rdr)).used_cohorts[rec, k])
                    _cdid_impb_fail("Link control cohort is outside design block used_cohorts support.", 459)
            }
            cc = cc \ k
        }
        if (rows(cc) == 0) _cdid_impb_fail("Each link requires a nonempty control cohort set.", 498)
        kord = uniqrows(cc)
        if (rows(kord) != rows(cc)) _cdid_impb_fail("Duplicate control cohort within one link.", 498)
        cellv = J(rows(kord), 1, 0)
        for (k=1; k<=rows(kord); k++) {
            cellv[k] = m.cell_lookup[kord[k], pid]
            if (cellv[k] < 1) _cdid_impb_fail("First-step fit lacks a link control cell.", 498)
        }
        b.links[ell].link_index = ell
        b.links[ell].link_id = (rows(link_ids) ? link_ids[ell] : "link:" + strofreal(ell, "%21.0g"))
        b.links[ell].pair_id = pid
        b.links[ell].recipe = rec
        b.links[ell].L = m.pairs[pid].L
        b.links[ell].target_ci = ci
        b.links[ell].target_cell = cell
        b.links[ell].n_ctrl = rows(kord)
        b.links[ell].ctrl_ci = kord
        b.links[ell].ctrl_cell = cellv
        b.links[ell].or_block = 0
        b.links[ell].trt_id = 0
        b.links[ell].mu = .
        b.links[ell].mu_ipw = .
        b.links[ell].mu_reg = .
        b.links[ell].mu_dr = .
        if (rows(kord) > mlen) mlen = rows(kord)
    }

    if (b.has_or) {
        bkeys = J(0, 1 + mlen, 0)
        link_raw = J(L, 1, 0)
        for (ell=1; ell<=L; ell++) {
            krow = (b.links[ell].pair_id, b.links[ell].ctrl_ci', J(1, mlen - b.links[ell].n_ctrl, -1))
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
        b.n_blocks = nb
        b.blocks = _cdid_impb_block(1, nb)
        for (bid=1; bid<=nb; bid++) {
            raw = ordb[bid]
            pid = bkeys[raw, 1]
            kr = bkeys[raw, 2..cols(bkeys)]
            kord = select(kr', kr' :> 0)
            b.blocks[bid].block_id = bid
            b.blocks[bid].pair_id = pid
            b.blocks[bid].recipe = m.pairs[pid].recipe
            b.blocks[bid].L = m.pairs[pid].L
            b.blocks[bid].n_ctrl = rows(kord)
            b.blocks[bid].ctrl_ci = kord
            b.blocks[bid].ctrl_cell = m.cell_lookup[kord, pid]
            b.blocks[bid].links = _cdid_impb_ix(rank_b[link_raw] :== bid)
            b.blocks[bid].name = (m.pairs[pid].pair_label + ":K={" +
                invtokens(m.cohort_labels[kord]', ",") + "}" +
                (b.is_dated | m.n_recipes > 1 ? ":r" + strofreal(b.blocks[bid].recipe, "%21.0g") : ""))
        }
        for (ell=1; ell<=L; ell++) b.links[ell].or_block = rank_b[link_raw[ell]]
    }

    if (b.has_ps) {
        tkeys = J(0, 2 + mlen, 0)
        tvec = J(L, 1, 0)
        for (ell=1; ell<=L; ell++) {
            krow = (b.links[ell].pair_id, b.links[ell].target_ci, b.links[ell].ctrl_ci',
                J(1, mlen - b.links[ell].n_ctrl, -1))
            found = 0
            for (i=1; i<=rows(tkeys); i++) {
                if (all(tkeys[i, ] :== krow)) {
                    found = i
                    break
                }
            }
            if (found == 0) {
                tkeys = tkeys \ krow
                found = rows(tkeys)
            }
            tvec[ell] = found
            b.links[ell].trt_id = found
        }
        nt = rows(tkeys)
        b.n_trt = nt
        b.trts = _cdid_impb_trt(1, nt)
        for (tid=1; tid<=nt; tid++) {
            ell = _cdid_impb_pos(tvec, tid)
            b.trts[tid].trt_id = tid
            b.trts[tid].rep_tid = tid
            b.trts[tid].pair_id = b.links[ell].pair_id
            b.trts[tid].recipe = b.links[ell].recipe
            b.trts[tid].target_ci = b.links[ell].target_ci
            b.trts[tid].target_cell = b.links[ell].target_cell
            b.trts[tid].all_ci = b.links[ell].target_ci \ b.links[ell].ctrl_ci
            b.trts[tid].n_coh = rows(b.trts[tid].all_ci)
            b.trts[tid].links = _cdid_impb_ix(tvec :== tid)
            b.trts[tid].name = (m.pairs[b.links[ell].pair_id].pair_label + ":g=" +
                m.cohort_labels[b.links[ell].target_ci] + ":K={" +
                invtokens(m.cohort_labels[b.links[ell].ctrl_ci]', ",") + "}" +
                (b.is_dated | m.n_recipes > 1 ? ":r" + strofreal(b.trts[tid].recipe, "%21.0g") : ""))
        }
    }
}

/* Preflight of concurrently live retained bytes (including the borrowed Stage 1 state and the
   borrowed treatment design / dated reader), the largest fit/replay scratch and the opt-in
   work-ops estimate against the caller budget. */
void _cdid_impb_preflight(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m)
{
    real scalar p, q, G, ret, work, ell, bid, tid, k, nr, nrmax, cmax, nlb, nbb, ntb, M_b
    real scalar sc_ps, sc_or, sc_rep, scr, xa_cols
    real colvector ncoh
    p = b.p_raw
    q = b.q
    G = b.G
    xa_cols = b.n_recipes * G
    cmax = b.chunk_rows
    ncoh = J(G, 1, 0)
    for (k=1; k<=G; k++) ncoh[k] = rows(b.coh_u[k].rows)
    ret = 64 + b.n_cells * (2 + 2 * p)
    work = 0
    for (k=1; k<=b.n_cells; k++) work = work + ncoh[m.cells[k].cohort_index] * (p + 2)^2
    for (ell=1; ell<=b.n_links; ell++) {
        ret = ret + 48 + 3 * p + 4 * q + b.links[ell].n_ctrl * (4 + p)
        for (k=1; k<=b.links[ell].n_ctrl; k++) {
            work = work + ncoh[b.links[ell].ctrl_ci[k]] * (2 * p + 3 * q + 8)
        }
    }
    for (bid=1; bid<=b.n_blocks; bid++) {
        ret = ret + 32 + p + rows(b.blocks[bid].links) + b.blocks[bid].n_ctrl * (6 + 3 * p + p * p)
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
            work = work + 3 * ncoh[b.blocks[bid].ctrl_ci[k]] * (p + 2)^2
        }
    }
    nrmax = 0
    for (tid=1; tid<=b.n_trt; tid++) {
        nr = 0
        for (k=1; k<=b.trts[tid].n_coh; k++) nr = nr + ncoh[b.trts[tid].all_ci[k]]
        if (nr > nrmax) nrmax = nr
        ret = ret + 32 + 2 * q * q + q * (b.trts[tid].n_coh + 1) + b.trts[tid].n_coh
        work = work + nr * (q + 2)^2 * (b.glm_maxit + 2)
    }
    if (b.has_or) {
        if (b.outcome_pooling == "separate")
            ret = ret + b.n_blocks * (4 * p * p + 4 * p + b.n_cells * p * (p + 1) + b.n_groups * p * p)
        else
            ret = ret + 4 * p * p + 4 * p + b.n_cells * p * (p + 1) + b.n_groups * p * p
    }
    if (b.is_dated) ret = 8 * ret + m.retained_bytes_bound
    else ret = 8 * ret + m.retained_bytes_bound + 8 * b.N * q
    sc_ps = nrmax * (8 * q + 16)
    sc_or = 0
    if (b.has_or) sc_or = 2 * (p + 1) * b.n_blocks * (3 * p + 8) + b.n_cells * (p + 2) * (p + 4) + 2 * cmax * (p + 4)
    nlb = min((b.link_batch, b.n_links))
    nbb = min((nlb, b.n_blocks))
    ntb = min((nlb, b.n_trt))
    M_b = nlb * (1 + b.has_ps) + nbb
    sc_rep = (cmax * (xa_cols + 4 * q * ntb + 9 * p + 3 * nbb + M_b + 48) +
        b.n_cells * M_b * (1 + p) + b.n_groups * p * M_b)
    scr = 8 * max((sc_ps, sc_or, sc_rep)) + 4096
    work = 16 * work
    _cdid_imp_check_budget(b.budget, ret, scr, work, "imputation backend preflight")
    b.retained_bytes = ret
    b.scratch_bytes = scr
    b.work_ops = work
}

/* Full-F synthetic cell moments of every Stage 1 cell (source cells_fitted):
   M_one = H_N[W a_k I_k], Y_one = H_N[W a_k I_k j], M_z = H_N[W a_k I_k z],
   Y_z = H_N[W a_k I_k j z] over all raw slope columns. */
void _cdid_impb_cell_moments(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    struct _cdid_impb_rcc rowvector rcc, struct _cdid_impb_cjc rowvector cjc, real scalar use_cache)
{
    real scalar cid, ci, r, ak, n_coh, lo, hi, p, j, slot, s_v
    real colvector coh, rch, wr, v, jh, cids
    real rowvector allc, s_vz
    real matrix Zr
    p = b.p_raw
    allc = _cdid_impb_seq(1, p)
    b.cell_M_one = J(b.n_cells, 1, 0)
    b.cell_Y_one = J(b.n_cells, 1, 0)
    b.cell_M_z = J(b.n_cells, p, 0)
    b.cell_Y_z = J(b.n_cells, p, 0)
    for (r=1; r<=b.n_recipes; r++) {
        for (ci=1; ci<=b.G; ci++) {
            cids = J(0, 1, 0)
            for (cid=1; cid<=b.n_cells; cid++) {
                if (m.cells[cid].recipe == r & m.cells[cid].cohort_index == ci) cids = cids \ cid
            }
            if (rows(cids) == 0) continue
            ak = _cdid_impb_ak(b, m, r, ci)
            coh = _cdid_impb_crows(b, m, r, ci)
            n_coh = rows(coh)
            if (use_cache & n_coh <= b.chunk_rows) {
                slot = (r - 1) * b.G + ci
                if (n_coh > 0) {
                    wr = _cdid_impb_w(b, m, r, coh, 1)
                    v = wr :* ak
                    Zr = (p > 0 ? _cdid_impb_z(b, m, r, coh, allc, 1) : J(0, 0, 0))
                }
                else {
                    wr = J(0, 1, 0)
                    v = J(0, 1, 0)
                    Zr = J(0, p, 0)
                }
                rcc[slot].wr = wr
                rcc[slot].v = v
                rcc[slot].Zr = Zr
                rcc[slot].has_x = 0
                rcc[slot].ready = 1
                s_v = (n_coh > 0 ? quadsum(v) : 0)
                s_vz = (p > 0 & n_coh > 0 ? quadcross(v, Zr) : J(1, p, 0))
                for (j=1; j<=rows(cids); j++) {
                    cid = cids[j]
                    jh = (n_coh > 0 ? _cdid_impb_pred_z(m, cid, n_coh, Zr) : J(0, 1, 0))
                    cjc[cid].jh = jh
                    cjc[cid].ready = 1
                    b.cell_M_one[cid] = s_v
                    b.cell_Y_one[cid] = (n_coh > 0 ? quadsum(v :* jh) : 0)
                    if (p > 0) {
                        b.cell_M_z[cid, ] = s_vz
                        if (n_coh > 0) b.cell_Y_z[cid, ] = quadcross(v :* jh, Zr)
                    }
                }
                continue
            }
            for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
                hi = min((n_coh, lo + b.chunk_rows - 1))
                rch = coh[lo..hi]
                v = _cdid_impb_w(b, m, r, rch, 1) :* ak
                Zr = (p > 0 ? _cdid_impb_z(b, m, r, rch, allc, 1) : J(rows(rch), 0, 0))
                s_v = quadsum(v)
                s_vz = (p > 0 ? quadcross(v, Zr) : J(1, 0, 0))
                for (j=1; j<=rows(cids); j++) {
                    cid = cids[j]
                    jh = _cdid_impb_pred_z(m, cid, rows(rch), Zr)
                    b.cell_M_one[cid] = b.cell_M_one[cid] + s_v
                    b.cell_Y_one[cid] = b.cell_Y_one[cid] + quadsum(v :* jh)
                    if (p > 0) {
                        b.cell_M_z[cid, ] = b.cell_M_z[cid, ] + s_vz
                        b.cell_Y_z[cid, ] = b.cell_Y_z[cid, ] + quadcross(v :* jh, Zr)
                    }
                }
            }
        }
    }
    b.cell_M_one = b.cell_M_one :/ b.N
    b.cell_Y_one = b.cell_Y_one :/ b.N
    b.cell_M_z = b.cell_M_z :/ b.N
    b.cell_Y_z = b.cell_Y_z :/ b.N
    if (missing(b.cell_M_one) | missing(b.cell_Y_one) | missing(b.cell_M_z) | missing(b.cell_Y_z))
        _cdid_impb_fail("Nonfinite imputation full-F cell moments.", 504)
}

/* Maximum full-F column scale of the proposed slope columns over every used support
   (control cells of all blocks and target cells of all links), as or_pooled_fit_blocks. */
real rowvector _cdid_impb_or_colmax(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real rowvector proposed)
{
    real colvector used
    real rowvector out
    real scalar cid, bid, k, ell, row
    used = J(b.n_cells, 1, 0)
    for (bid=1; bid<=b.n_blocks; bid++) {
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) used[b.blocks[bid].ctrl_cell[k]] = 1
    }
    for (ell=1; ell<=b.n_links; ell++) used[b.links[ell].target_cell] = 1
    out = J(1, cols(proposed), 0)
    for (cid=1; cid<=b.n_cells; cid++) {
        if (!used[cid]) continue
        row = (m.cells[cid].recipe - 1) * m.n_cohorts + m.cells[cid].cohort_index
        out = colmax(out \ m.cohort_colmax_F[row, proposed])
    }
    return(out)
}

/* Strict reconstruction of proposed OR slope columns on actual original rows: every full-F
   cohort row of each block control cell (control role) and of each consuming link target
   cell (treated role), with the block constants (or_pooled_alias_check). */
void _cdid_impb_or_alias(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real rowvector retained, real rowvector proposed, real matrix map, real matrix constants_Z,
    real rowvector colmax_prop, real rowvector ctl_err, real rowvector trt_err)
{
    real scalar n_prop, cid, bid, ell, t, bb, ci, r, n_coh, lo, hi, j, den, ratio, is_c, is_t
    real colvector cb, tb, ub, coh, rch, u, s0
    real matrix Zs, u_mat, s0_mat
    n_prop = cols(proposed)
    ctl_err = J(1, n_prop, 0)
    trt_err = J(1, n_prop, 0)
    b.or_rows_checked = 0
    for (cid=1; cid<=b.n_cells; cid++) {
        cb = J(0, 1, 0)
        for (bid=1; bid<=b.n_blocks; bid++) {
            if (_cdid_impb_pos(b.blocks[bid].ctrl_cell, cid) > 0) cb = cb \ bid
        }
        tb = J(0, 1, 0)
        for (ell=1; ell<=b.n_links; ell++) {
            if (b.links[ell].target_cell == cid) tb = tb \ b.links[ell].or_block
        }
        if (rows(cb) + rows(tb) == 0) continue
        ub = uniqrows(cb \ tb)
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        coh = _cdid_impb_crows(b, m, r, ci)
        n_coh = rows(coh)
        for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
            hi = min((n_coh, lo + b.chunk_rows - 1))
            rch = coh[lo..hi]
            Zs = _cdid_impb_z(b, m, r, rch, _cdid_impb_seq(1, b.p_raw), 1)
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
            b.or_rows_checked = b.or_rows_checked + rows(Zs)
        }
    }
}

/* Block-scoped maximum full-F column scale of proposed slope columns for separate OR:
   strictly control cells of block bid and target cells of links consuming block bid. */
real rowvector _cdid_impb_or_colmax_blk(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar bid, real rowvector proposed)
{
    real colvector used, lk_b
    real rowvector out
    real scalar cid, k, i, row
    used = J(b.n_cells, 1, 0)
    for (k=1; k<=b.blocks[bid].n_ctrl; k++) used[b.blocks[bid].ctrl_cell[k]] = 1
    lk_b = b.blocks[bid].links
    for (i=1; i<=rows(lk_b); i++) used[b.links[lk_b[i]].target_cell] = 1
    out = J(1, cols(proposed), 0)
    for (cid=1; cid<=b.n_cells; cid++) {
        if (!used[cid]) continue
        row = (m.cells[cid].recipe - 1) * m.n_cohorts + m.cells[cid].cohort_index
        out = colmax(out \ m.cohort_colmax_F[row, proposed])
    }
    return(out)
}

/* Block-scoped alias verification for separate OR: strictly block bid control cells
   (control role) and consuming link target cells (treated role). */
void _cdid_impb_or_alias_blk(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real scalar bid, real rowvector retained, real rowvector proposed, real matrix map,
    real rowvector const_row, real rowvector colmax_prop, real rowvector ctl_err, real rowvector trt_err)
{
    real scalar n_prop, cid, i, ci, r, n_coh, lo, hi, j, den, ratio, is_c, is_t
    real colvector lk_b, coh, rch, u, s0
    real matrix Zs, u_mat, s0_mat
    n_prop = cols(proposed)
    ctl_err = J(1, n_prop, 0)
    trt_err = J(1, n_prop, 0)
    lk_b = b.blocks[bid].links
    for (cid=1; cid<=b.n_cells; cid++) {
        is_c = (_cdid_impb_pos(b.blocks[bid].ctrl_cell, cid) > 0)
        is_t = 0
        for (i=1; i<=rows(lk_b); i++) {
            if (b.links[lk_b[i]].target_cell == cid) {
                is_t = 1
                break
            }
        }
        if (!is_c & !is_t) continue
        ci = m.cells[cid].cohort_index
        r = m.cells[cid].recipe
        coh = _cdid_impb_crows(b, m, r, ci)
        n_coh = rows(coh)
        for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
            hi = min((n_coh, lo + b.chunk_rows - 1))
            rch = coh[lo..hi]
            Zs = _cdid_impb_z(b, m, r, rch, _cdid_impb_seq(1, b.p_raw), 1)
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
            b.or_rows_checked = b.or_rows_checked + rows(Zs)
        }
    }
}

/* Stage 2: counterfactual OR over deduplicated blocks (pooled common slope or separate per-block). */
void _cdid_impb_fit_or(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    struct _cdid_impb_rcc rowvector rcc, struct _cdid_impb_cjc rowvector cjc)
{
    struct _cdid_impb_cacc rowvector acc
    struct cdid_wls_result scalar wfit
    real scalar p, nb, bid, k, cid, ci, r, ak, n_coh, lo, hi, j, nn, scale, at, nPseudo, slot
    real scalar thr, orig_norm, res_norm, rk, k_keep, exact_j, isfirst, wsum, L, n_b, tot_rss
    real colvector coh, rch, v, jh, rr, ix, pseudo_y, pseudo_w, pseudo_block, qy, coef, vv
    real colvector yc, bvec, avec, mean_y, py_b, pw_b, blk_one
    real rowvector anchor, constb, coefrow, retained, proposed, aliased, aliased_pos
    real rowvector candidates, keep, colmax_prop, ctl_err, trt_err, allc, md_raw_b, const_row
    real matrix D, Fb, Q, R1, V, tau, pivot, res_fac, A_plus, A_minus, pseudo_X, C, px_b, Cb
    real matrix mean_D_raw, map, constants_Z, map_full, Z_keep, mean_D, Rfit

    p = b.p_raw
    nb = b.n_blocks
    allc = _cdid_impb_seq(1, p)
    acc = _cdid_impb_cacc(1, b.n_cells)
    for (cid=1; cid<=b.n_cells; cid++) {
        acc[cid].used = 0
        acc[cid].n = 0
        acc[cid].wsum = 0
        acc[cid].F = J(0, p + 2, 0)
        acc[cid].anchor = J(1, p, 0)
        acc[cid].constant = J(1, p, 1)
    }
    for (bid=1; bid<=nb; bid++) {
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) acc[b.blocks[bid].ctrl_cell[k]].used = 1
    }

    /* 2a: bounded streaming of each used control cell: sqrt(v) [1, L z, j] summaries. */
    for (r=1; r<=b.n_recipes; r++) {
        for (ci=1; ci<=b.G; ci++) {
            ix = J(0, 1, 0)
            for (cid=1; cid<=b.n_cells; cid++) {
                if (acc[cid].used & m.cells[cid].recipe == r & m.cells[cid].cohort_index == ci) ix = ix \ cid
            }
            if (rows(ix) == 0) continue
            ak = _cdid_impb_ak(b, m, r, ci)
            coh = _cdid_impb_crows(b, m, r, ci)
            n_coh = rows(coh)
            slot = (r - 1) * b.G + ci
            if (rcc[slot].ready & n_coh <= b.chunk_rows) {
                if (n_coh > 0) {
                    v = rcc[slot].v
                    Fb = rcc[slot].Zr
                    wsum = quadsum(v)
                    for (k=1; k<=rows(ix); k++) {
                        cid = ix[k]
                        L = m.cells[cid].L
                        jh = (cjc[cid].ready ? cjc[cid].jh : _cdid_impb_pred_z(m, cid, n_coh, Fb))
                        D = L :* Fb
                        if (missing(D) | missing(jh) | missing(v))
                            _cdid_impb_fail("Nonfinite synthetic counterfactual OR support rows.", 504)
                        if (p > 0) {
                            acc[cid].anchor = D[1, ]
                            acc[cid].constant = colmin(D :== acc[cid].anchor)
                        }
                        acc[cid].F = _cdid_imp_compress((J(n_coh, 1, 1), D, jh) :* sqrt(v))
                        acc[cid].n = n_coh
                        acc[cid].wsum = wsum
                    }
                }
                for (k=1; k<=rows(ix); k++) {
                    cid = ix[k]
                    if (acc[cid].n < 1 | !(acc[cid].wsum > 0) | acc[cid].wsum >= .)
                        _cdid_impb_fail("Shared OR control cell " + m.cells[cid].name + " has no positive synthetic full-F support.", 459)
                }
                continue
            }
            isfirst = 1
            for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
                hi = min((n_coh, lo + b.chunk_rows - 1))
                rch = coh[lo..hi]
                v = _cdid_impb_w(b, m, r, rch, 1) :* ak
                Fb = (p > 0 ? _cdid_impb_z(b, m, r, rch, allc, 1) : J(rows(rch), 0, 0))
                for (k=1; k<=rows(ix); k++) {
                    cid = ix[k]
                    L = m.cells[cid].L
                    jh = _cdid_impb_pred_z(m, cid, rows(rch), Fb)
                    D = L :* Fb
                    if (missing(D) | missing(jh) | missing(v))
                        _cdid_impb_fail("Nonfinite synthetic counterfactual OR support rows.", 504)
                    if (p > 0) {
                        if (isfirst) acc[cid].anchor = D[1, ]
                        acc[cid].constant = acc[cid].constant :& colmin(D :== acc[cid].anchor)
                    }
                    acc[cid].F = _cdid_imp_compress(acc[cid].F \ ((J(rows(rch), 1, 1), D, jh) :* sqrt(v)))
                    acc[cid].n = acc[cid].n + rows(rch)
                    acc[cid].wsum = acc[cid].wsum + quadsum(v)
                }
                isfirst = 0
            }
            for (k=1; k<=rows(ix); k++) {
                cid = ix[k]
                if (acc[cid].n < 1 | !(acc[cid].wsum > 0) | acc[cid].wsum >= .)
                    _cdid_impb_fail("Shared OR control cell " + m.cells[cid].name + " has no positive synthetic full-F support.", 459)
            }
        }
    }

    /* 2b: merge supports per block; exact pseudo rows reproduce block means and crossproducts. */
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
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
            cid = b.blocks[bid].ctrl_cell[k]
            Fb = Fb \ acc[cid].F
            wsum = wsum + acc[cid].wsum
            if (p > 0) {
                if (k == 1) anchor = acc[cid].anchor
                constb = constb :& acc[cid].constant :& (acc[cid].anchor :== anchor)
            }
        }
        Fb = _cdid_imp_compress(Fb)
        b.blocks[bid].weight_sum = wsum
        b.blocks[bid].mass = wsum / b.N
        if (!(b.blocks[bid].mass > 0) | b.blocks[bid].mass >= .)
            _cdid_impb_fail("Shared OR block " + b.blocks[bid].name + " has nonpositive synthetic control mass.", 459)
        hqrd(Fb[, 1], Q=., tau=., R1=.)
        V = hqrdmultq(Q, tau, Fb[, 2..(p + 2)], 1)
        coefrow = V[1, ] / R1[1, 1]
        res_fac = V[2..rows(V), ]
        if (b.intercept & p > 0) {
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

    if (b.outcome_pooling == "separate") {
        thr = 64 * epsilon(1) * max((1, p + b.intercept))
        b.or_rows_checked = 0
        b.or_threshold = thr
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
            if (b.intercept & p > 0) {
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
                if (orig_norm == 0 | res_norm <= b.rank_tol * orig_norm) proposed = (proposed, j)
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
                if (b.intercept) {
                    const_row = (md_raw_b[proposed] -
                        (cols(retained) ? md_raw_b[retained] * map : J(1, cols(proposed), 0))) :/ b.blocks[bid].L
                }
                colmax_prop = _cdid_impb_or_colmax_blk(b, m, bid, proposed)
                _cdid_impb_or_alias_blk(b, m, bid, retained, proposed, map, const_row, colmax_prop, ctl_err, trt_err)
                for (j=1; j<=cols(proposed); j++) {
                    exact_j = (ctl_err[j] <= thr & ctl_err[j] < .)
                    if (exact_j & (trt_err[j] > thr | trt_err[j] >= .))
                        _cdid_impb_fail("Control design of counterfactual OR block " + b.blocks[bid].name + " has an exact null direction that does not reconstruct the required treated predictions; no ridge or singular-contrast solver is used; model unchanged. This computation failure does not establish population ATT nonidentification.", 459)
                    if (exact_j) {
                        aliased = (aliased, proposed[j])
                        aliased_pos = (aliased_pos, j)
                    }
                    else candidates = (candidates, proposed[j])
                }
                keep = (cols(retained) + cols(candidates) > 0 ? sort((retained, candidates)', 1)' : J(1, 0, 0))
                if (!b.intercept & cols(keep) == 0)
                    _cdid_impb_fail("Counterfactual OR block " + b.blocks[bid].name + " design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
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
                if (!b.intercept & cols(keep) == 0)
                    _cdid_impb_fail("Counterfactual OR block " + b.blocks[bid].name + " design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
            }
            k_keep = cols(keep)
            if (cols(candidates) == 0) {
                wfit = cdid_pooled_wls_fit((k_keep > 0 ? px_b[, keep] : J(n_b, 0, 0)),
                    py_b, pw_b, blk_one, blk_one, b.N,
                    b.intercept, J(0, k_keep, 0), J(0, 1, 0), b.rank_tol)
                if (cols(wfit.kept) != k_keep)
                    _cdid_impb_fail("Counterfactual OR block final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
                bvec = wfit.beta
                avec = wfit.alpha
                Rfit = wfit.R
                pivot = wfit.pivot
                b.blocks[bid].or_rcond = wfit.rcond_R
                mean_D = (b.intercept & k_keep > 0 ? wfit.mean_X : J(1, k_keep, 0))
                mean_y = wfit.mean_y
                b.blocks[bid].or_rss = wfit.rss
            }
            else {
                mean_y = J(1, 1, sum(pw_b :* py_b) / sum(pw_b))
                yc = (b.intercept ? py_b :- mean_y[1] : py_b) :* sqrt(pw_b)
                Z_keep = Cb[, keep]
                hqrdp(Z_keep, Q=., tau=., R1=., pivot=.)
                if (rows(R1) != k_keep | missing(R1) |
                    min(abs(diagonal(R1))) <= b.rank_tol * max(abs(diagonal(R1))))
                    _cdid_impb_fail("Counterfactual OR block final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
                qy = hqrdmultq1t(Q, tau, yc)
                coef = solveupper(R1, qy, rk=., 0)
                if (rk != k_keep | missing(coef)) _cdid_impb_fail("Nonfinite counterfactual OR block solve.", 504)
                bvec = coef[invorder(pivot)]
                Rfit = R1
                b.blocks[bid].or_rcond = _cdid_nuisance_rcond(R1)
                mean_D = (b.intercept ? md_raw_b[keep] : J(1, k_keep, 0))
                avec = (b.intercept ? mean_y - mean_D * bvec : J(1, 1, 0))
                vv = py_b - avec[blk_one] - px_b[, keep] * bvec
                b.blocks[bid].or_rss = sum(pw_b :* vv:^2)
            }
            if (missing(bvec) | missing(avec))
                _cdid_impb_fail("Nonfinite counterfactual OR block coefficients.", 504)
            b.blocks[bid].or_p = k_keep
            b.blocks[bid].or_keep = keep
            b.blocks[bid].or_beta = bvec
            b.blocks[bid].or_R = Rfit
            b.blocks[bid].or_pivot = pivot
            b.blocks[bid].or_alias_map = map
            b.blocks[bid].or_alias_const = const_row
            b.blocks[bid].or_adj_ready = 0
            b.blocks[bid].alpha = (b.intercept ? avec[1] : 0)
            b.blocks[bid].mean_d = (b.intercept & k_keep > 0 ? mean_D[1, ] : J(1, k_keep, 0))
            b.blocks[bid].mean_y = mean_y[1]
            tot_rss = tot_rss + b.blocks[bid].or_rss
        }
        b.or_p = 0
        b.or_rss = tot_rss
        b.n_or_fits = nb
        return
    }

    /* 2c: original-order near-dependency proposal on the partialled weighted pooled factor. */
    C = pseudo_X
    mean_D_raw = J(nb, p, 0)
    for (bid=1; bid<=nb; bid++) {
        ix = _cdid_impb_ix(pseudo_block :== bid)
        if (b.intercept & p > 0) {
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
        if (orig_norm == 0 | res_norm <= b.rank_tol * orig_norm) proposed = (proposed, j)
        else retained = (retained, j)
    }

    /* 2d: exact aliases must reconstruct control AND required treated full-F rows. */
    thr = 64 * epsilon(1) * max((1, p + b.intercept))
    aliased = J(1, 0, 0)
    aliased_pos = J(1, 0, 0)
    candidates = J(1, 0, 0)
    map = J(cols(retained), 0, 0)
    constants_Z = J(nb, 0, 0)
    ctl_err = J(1, 0, 0)
    trt_err = J(1, 0, 0)
    b.or_rows_checked = 0
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
        if (b.intercept) {
            for (bid=1; bid<=nb; bid++) {
                constants_Z[bid, ] = (mean_D_raw[bid, proposed] -
                    (cols(retained) ? mean_D_raw[bid, retained] * map : J(1, cols(proposed), 0))) :/ b.blocks[bid].L
            }
        }
        colmax_prop = _cdid_impb_or_colmax(b, m, proposed)
        _cdid_impb_or_alias(b, m, retained, proposed, map, constants_Z, colmax_prop, ctl_err, trt_err)
        for (j=1; j<=cols(proposed); j++) {
            exact_j = (ctl_err[j] <= thr & ctl_err[j] < .)
            if (exact_j & (trt_err[j] > thr | trt_err[j] >= .))
                _cdid_impb_fail("Pooled control design of the shared counterfactual OR has an exact null direction that does not reconstruct the required treated predictions; no ridge or singular-contrast solver is used; model unchanged. This computation failure does not establish population ATT nonidentification.", 459)
            if (exact_j) {
                aliased = (aliased, proposed[j])
                aliased_pos = (aliased_pos, j)
            }
            else candidates = (candidates, proposed[j])
        }
        keep = (cols(retained) + cols(candidates) > 0 ? sort((retained, candidates)', 1)' : J(1, 0, 0))
        if (!b.intercept & cols(keep) == 0)
            _cdid_impb_fail("Shared OR design has no retained columns on actual control and treated rows; no reduced control regression can be computed.", 459)
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

    /* 2e: one solve of the common slope on the kept columns (QR rank rule). */
    k_keep = cols(keep)
    if (cols(candidates) == 0) {
        wfit = cdid_pooled_wls_fit((k_keep > 0 ? pseudo_X[, keep] : J(nPseudo, 0, 0)),
            pseudo_y, pseudo_w, pseudo_block, J(nPseudo, 1, 1), b.N,
            b.intercept, J(0, k_keep, 0), J(0, 1, 0), b.rank_tol)
        if (cols(wfit.kept) != k_keep)
            _cdid_impb_fail("Shared OR final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
        bvec = wfit.beta
        avec = wfit.alpha
        Rfit = wfit.R
        pivot = wfit.pivot
        b.or_rcond = wfit.rcond_R
        mean_D = (b.intercept & k_keep > 0 ? wfit.mean_X : J(nb, k_keep, 0))
        mean_y = wfit.mean_y
        b.or_rss = wfit.rss
    }
    else {
        yc = pseudo_y
        mean_y = J(nb, 1, 0)
        for (bid=1; bid<=nb; bid++) {
            ix = _cdid_impb_ix(pseudo_block :== bid)
            mean_y[bid] = sum(pseudo_w[ix] :* pseudo_y[ix]) / sum(pseudo_w[ix])
            if (b.intercept) yc[ix] = pseudo_y[ix] :- mean_y[bid]
        }
        yc = yc :* sqrt(pseudo_w)
        Z_keep = C[, keep]
        hqrdp(Z_keep, Q=., tau=., R1=., pivot=.)
        if (rows(R1) != k_keep | missing(R1) |
            min(abs(diagonal(R1))) <= b.rank_tol * max(abs(diagonal(R1))))
            _cdid_impb_fail("Shared OR final rank refused; no ridge or parameter rescue. This computation failure does not establish population ATT nonidentification.", 459)
        qy = hqrdmultq1t(Q, tau, yc)
        coef = solveupper(R1, qy, rk=., 0)
        if (rk != k_keep | missing(coef)) _cdid_impb_fail("Nonfinite shared OR solve.", 504)
        bvec = coef[invorder(pivot)]
        Rfit = R1
        b.or_rcond = _cdid_nuisance_rcond(R1)
        mean_D = (b.intercept ? mean_D_raw[, keep] : J(nb, k_keep, 0))
        avec = (b.intercept ? mean_y - mean_D * bvec : J(nb, 1, 0))
        vv = pseudo_y - avec[pseudo_block] - pseudo_X[, keep] * bvec
        b.or_rss = sum(pseudo_w :* vv:^2)
    }
    if (missing(bvec) | missing(avec))
        _cdid_impb_fail("Nonfinite shared counterfactual OR coefficients.", 504)
    b.or_p = k_keep
    b.or_keep = keep
    b.or_retained = retained
    b.or_proposed = proposed
    b.or_dropped = aliased
    b.or_candidates = candidates
    b.or_ctl_scaled = ctl_err
    b.or_trt_scaled = trt_err
    b.or_alias_map = map
    b.or_alias_const = constants_Z
    b.or_threshold = thr
    b.or_beta = bvec
    b.or_R = Rfit
    b.or_pivot = pivot
    for (bid=1; bid<=nb; bid++) {
        b.blocks[bid].alpha = (b.intercept ? avec[bid] : 0)
        b.blocks[bid].mean_d = (b.intercept & k_keep > 0 ? mean_D[bid, ] : J(1, k_keep, 0))
        b.blocks[bid].mean_y = mean_y[bid]
    }
    b.n_or_fits = 1
}

/* Stage 2 upstream forcing per (block, control cohort k), source or_block_fits:
   a_force_alpha = H_N[W_r I_k r], a_force_beta = H_N[W_r I_k r dtilde] with r = j - m;
   jf_alpha_alpha = M_one (intercept), jf_alpha_beta = M_d_kept = L M_z[kept];
   jf_beta_alpha = H_N[v I_k dtilde] (intercept), jf_beta_beta = H_N[v I_k dtilde d_imp']. */
void _cdid_impb_or_forcing(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    struct _cdid_impb_rcc rowvector rcc, struct _cdid_impb_cjc rowvector cjc)
{
    real scalar p, po, bid, k, nk, cid, gid, pg, ci, r, ak, L, n_coh, lo, hi, alpha_b, s_r, sep, slot
    real colvector coh, rch, wr, v, jh, mh, res, obeta, ak_vec
    real rowvector kept_g, okeep, s_rd, s_vd, allc
    real matrix Zr, dor, dt, dimp, s_vdd, afb, jba, jab, jbb
    p = b.p_raw
    allc = _cdid_impb_seq(1, p)
    sep = (b.outcome_pooling == "separate")
    for (bid=1; bid<=b.n_blocks; bid++) {
        po = (sep ? b.blocks[bid].or_p : b.or_p)
        okeep = (sep ? b.blocks[bid].or_keep : b.or_keep)
        obeta = (sep ? b.blocks[bid].or_beta : b.or_beta)
        nk = b.blocks[bid].n_ctrl
        L = b.blocks[bid].L
        alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
        b.blocks[bid].a_force_alpha = J(nk, 1, 0)
        b.blocks[bid].jf_alpha_alpha = J(nk, 1, 0)
        ak_vec = J(nk, 1, 0)
        afb = J(po, nk, 0)
        jba = J(po, nk, 0)
        jab = J(nk, p, 0)
        jbb = J(po, nk * p, 0)
        for (k=1; k<=nk; k++) {
            cid = b.blocks[bid].ctrl_cell[k]
            gid = m.cells[cid].group_id
            kept_g = m.groups[gid].kept
            pg = m.groups[gid].p
            ci = m.cells[cid].cohort_index
            r = m.cells[cid].recipe
            ak = _cdid_impb_ak(b, m, r, ci)
            ak_vec[k] = ak
            s_r = 0
            s_rd = J(1, po, 0)
            s_vd = J(1, po, 0)
            s_vdd = J(po, pg, 0)
            coh = _cdid_impb_crows(b, m, r, ci)
            n_coh = rows(coh)
            slot = (r - 1) * b.G + ci
            if (rcc[slot].ready & n_coh <= b.chunk_rows) {
                if (n_coh > 0) {
                    wr = rcc[slot].wr
                    v = rcc[slot].v
                    Zr = rcc[slot].Zr
                    jh = (cjc[cid].ready ? cjc[cid].jh : _cdid_impb_pred_z(m, cid, n_coh, Zr))
                    if (po > 0) {
                        dor = L :* Zr[, okeep]
                        mh = alpha_b :+ dor * obeta
                    }
                    else mh = J(n_coh, 1, alpha_b)
                    res = jh - mh
                    s_r = quadsum(wr :* res)
                    if (po > 0) {
                        dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                        s_rd = quadcross(wr :* res, dt)
                        s_vd = quadcross(v, dt)
                        if (pg > 0) {
                            dimp = L :* Zr[, kept_g]
                            s_vdd = quadcross(v :* dt, dimp)
                        }
                    }
                }
            }
            else {
                for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
                    hi = min((n_coh, lo + b.chunk_rows - 1))
                    rch = coh[lo..hi]
                    wr = _cdid_impb_w(b, m, r, rch, 1)
                    v = wr :* ak
                    Zr = (p > 0 ? _cdid_impb_z(b, m, r, rch, allc, 1) : J(rows(rch), 0, 0))
                    jh = _cdid_impb_pred_z(m, cid, rows(rch), Zr)
                    if (po > 0) {
                        dor = L :* Zr[, okeep]
                        mh = alpha_b :+ dor * obeta
                    }
                    else mh = J(rows(rch), 1, alpha_b)
                    res = jh - mh
                    s_r = s_r + quadsum(wr :* res)
                    if (po > 0) {
                        dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                        s_rd = s_rd + quadcross(wr :* res, dt)
                        s_vd = s_vd + quadcross(v, dt)
                        if (pg > 0) {
                            dimp = L :* Zr[, kept_g]
                            s_vdd = s_vdd + quadcross(v :* dt, dimp)
                        }
                    }
                }
            }
            b.blocks[bid].a_force_alpha[k] = s_r / b.N
            b.blocks[bid].jf_alpha_alpha[k] = (b.intercept ? b.cell_M_one[cid] : 0)
            if (pg > 0) jab[k, kept_g] = L :* b.cell_M_z[cid, kept_g]
            if (po > 0) {
                afb[, k] = s_rd' / b.N
                if (b.intercept) jba[, k] = s_vd' / b.N
                if (pg > 0) jbb[, (k - 1) * p :+ kept_g] = s_vdd / b.N
            }
        }
        b.blocks[bid].a_force_beta = afb
        b.blocks[bid].jf_beta_alpha = jba
        b.blocks[bid].jf_alpha_beta = jab
        b.blocks[bid].jf_beta_beta = jbb
        if (sep) {
            b.blocks[bid].score_mean_a = sum(ak_vec :* b.blocks[bid].a_force_alpha)
            b.blocks[bid].score_mean_b = (po > 0 ? afb * ak_vec : J(0, 1, 0))
        }
        else {
            b.blocks[bid].score_mean_a = 0
            b.blocks[bid].score_mean_b = J(0, 1, 0)
        }
        if (missing(b.blocks[bid].a_force_alpha) | missing(afb) | missing(jba) | missing(jab) | missing(jbb))
            _cdid_impb_fail("Nonfinite counterfactual OR upstream forcing.", 504)
    }
}

/* Stage 3: full-F treatment logit for one key on the target and ordered control cohorts
   under recipe rt = b.trts[tid].recipe: raw treatment columns X_rt, binary target indicator,
   weight v = W_rt a_{rt,k}. Identical (recipe, target, ordered controls) keys reuse prior fit. */
void _cdid_impb_fit_ps(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real scalar tid)
{
    struct cdid_logit_result scalar lf
    real scalar k, nc, rt, prior, n_r, at, n_ku, same_rw
    real colvector risk, cir, D, wr, v, pr, act, sel, akc, cu
    real matrix X, vals, score, af
    nc = b.trts[tid].n_coh
    rt = b.trts[tid].recipe
    for (prior=1; prior<=tid-1; prior++) {
        if (b.trts[prior].rep_tid != prior) continue
        same_rw = (b.trts[prior].recipe == rt)
        if (!same_rw & b.is_dated) {
            same_rw = ((*(b.rdr)).recipes[b.trts[prior].recipe].x_slot == (*(b.rdr)).recipes[rt].x_slot &
                (*(b.rdr)).recipes[b.trts[prior].recipe].w_slot == (*(b.rdr)).recipes[rt].w_slot)
        }
        if (!same_rw | b.trts[prior].target_ci != b.trts[tid].target_ci |
            b.trts[prior].n_coh != nc) continue
        if (!all(b.trts[prior].all_ci :== b.trts[tid].all_ci)) continue
        b.trts[tid].rep_tid = prior
        b.trts[tid].gamma = b.trts[prior].gamma
        b.trts[tid].R = b.trts[prior].R
        b.trts[tid].pivot = b.trts[prior].pivot
        b.trts[tid].iterations = b.trts[prior].iterations
        b.trts[tid].converged = b.trts[prior].converged
        b.trts[tid].deviance = b.trts[prior].deviance
        b.trts[tid].score_max = b.trts[prior].score_max
        b.trts[tid].rcond_R = b.trts[prior].rcond_R
        b.trts[tid].n_risk = b.trts[prior].n_risk
        b.trts[tid].a_forcing = b.trts[prior].a_forcing
        b.trts[tid].n_capped = b.trts[prior].n_capped
        b.trts[tid].n_capped_ctrl = b.trts[prior].n_capped_ctrl
        b.trts[tid].p_min = b.trts[prior].p_min
        b.trts[tid].p_max = b.trts[prior].p_max
        return
    }
    n_r = 0
    for (k=1; k<=nc; k++) n_r = n_r + rows(_cdid_impb_crows(b, m, rt, b.trts[tid].all_ci[k]))
    risk = J(n_r, 1, 0)
    at = 0
    for (k=1; k<=nc; k++) {
        cu = _cdid_impb_crows(b, m, rt, b.trts[tid].all_ci[k])
        n_ku = rows(cu)
        if (n_ku > 0) {
            risk[(at + 1)..(at + n_ku)] = cu
            at = at + n_ku
        }
    }
    if (n_r > 1) risk = sort(risk, 1)
    cir = m.unit_cohort_index[risk]
    D = (cir :== b.trts[tid].target_ci)
    wr = _cdid_impb_w(b, m, rt, risk, 1)
    akc = m.cohort_a_F[rt, ]'
    v = wr :* akc[cir]
    X = _cdid_impb_x(b, rt, risk, 1)
    lf = cdid_logit_fit(X, D, v, b.glm_eps, b.glm_maxit)
    b.trts[tid].gamma = lf.beta
    b.trts[tid].R = lf.R
    b.trts[tid].pivot = lf.pivot
    b.trts[tid].iterations = lf.iterations
    b.trts[tid].converged = lf.converged
    b.trts[tid].deviance = lf.deviance
    b.trts[tid].score_max = lf.score_max
    b.trts[tid].rcond_R = lf.rcond_R
    b.trts[tid].n_risk = rows(risk)
    vals = _cdid_impb_logit_vals(X * lf.beta)
    pr = vals[, 1]
    act = vals[, 3]
    score = (wr :* (D - pr)) :* X
    af = J(cols(X), nc, 0)
    for (k=1; k<=nc; k++) {
        sel = _cdid_impb_ix(cir :== b.trts[tid].all_ci[k])
        if (rows(sel) > 0) af[, k] = quadcolsum(score[sel, ])' / b.N
    }
    b.trts[tid].a_forcing = af
    b.trts[tid].n_capped = sum(act :== 0)
    b.trts[tid].n_capped_ctrl = sum((act :== 0) :& (D :== 0))
    b.trts[tid].p_min = min(pr)
    b.trts[tid].p_max = max(pr)
    if (missing(af) | missing(lf.beta))
        _cdid_impb_fail("Nonfinite imputation treatment propensity fit or repeat-rate forcing.", 504)
}

/* Stage 4: link points and partial derivatives (source links_fitted). Normalization uses
   the original N_F; source formulas are evaluated verbatim, including identically-zero
   a_coeff_T_g and a_coeff_Bg_g. Unavailable arms stay missing/empty. */
void _cdid_impb_link_sums(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    struct _cdid_impb_rcc rowvector rcc, struct _cdid_impb_cjc rowvector cjc)
{
    struct _cdid_impb_lacc rowvector la
    real scalar ell, g, cg, rg, ag, L, gid, pg, d_plus, num_T, T, bid, alpha_b, B_gl, num_Bg
    real scalar tid, rtid, nk, k, ck, ckc, rk, ak, gk, pk, n_coh, lo, hi, p, po, q, sep
    real scalar d_minus, Cl, num_C, B_ol, mu, n_hit, h, tp, cp, n_r, slot
    real colvector coh, rch, wr, v, wm, wdot, wod, jh, mh, obeta, u_rtid, u_ckc
    real rowvector kept, keptk, okeep, dbp_or, allc, s_wm_tp, s_wod_tp
    real matrix X, Zr, vals, dor, hits, odds_m, act_m, jh_m, wm_m, wod_m, wdot_m, s_wdX_tp, s_wmj_tc, s_wodj_tc
    p = b.p_raw
    allc = _cdid_impb_seq(1, p)
    po = 0
    okeep = J(1, 0, 0)
    obeta = J(0, 1, 0)
    sep = (b.outcome_pooling == "separate")
    q = b.q
    if (b.has_ps) la = _cdid_impb_lacc(1, b.n_links)
    for (ell=1; ell<=b.n_links; ell++) {
        g = b.links[ell].target_cell
        cg = b.links[ell].target_ci
        rg = m.cells[g].recipe
        ag = _cdid_impb_ak(b, m, rg, cg)
        L = b.links[ell].L
        gid = m.cells[g].group_id
        kept = m.groups[gid].kept
        pg = m.groups[gid].p
        d_plus = b.cell_M_one[g]
        if (!(d_plus > 0) | d_plus >= .)
            _cdid_impb_fail("Link " + b.links[ell].link_id + " has nonpositive treated denominator d_{l,+}.", 459)
        num_T = b.cell_Y_one[g]
        T = num_T / d_plus
        b.links[ell].d_plus = d_plus
        b.links[ell].num_T = num_T
        b.links[ell].T_gl = T
        b.links[ell].d_bar_plus_j = (pg > 0 ? (L :* b.cell_M_z[g, kept]) :/ d_plus : J(1, 0, 0))
        b.links[ell].a_coeff_T_g = ((num_T / ag) - T * (d_plus / ag)) / d_plus
        b.links[ell].num_Bg = .
        b.links[ell].B_gl = .
        b.links[ell].a_coeff_Bg_g = 0
        b.links[ell].d_bar_plus_or = J(1, 0, 0)
        b.links[ell].d_minus = .
        b.links[ell].num_C = .
        b.links[ell].C_ol = .
        b.links[ell].num_Bo = .
        b.links[ell].B_ol = .
        b.links[ell].d_bar_minus_or = J(1, 0, 0)
        b.links[ell].d_bar_diff_or = J(1, 0, 0)
        b.links[ell].a_coeff_C_k = J(0, 1, 0)
        b.links[ell].a_coeff_Bo_k = J(0, 1, 0)
        b.links[ell].m_minus_k = J(0, 1, 0)
        b.links[ell].d_bar_minus_j_k = J(0, p, 0)
        b.links[ell].denminus_gamma = J(1, 0, 0)
        b.links[ell].g_C_gamma = J(1, 0, 0)
        b.links[ell].g_Bo_gamma = J(1, 0, 0)
        b.links[ell].g_dr_gamma = J(1, 0, 0)
        alpha_b = 0
        B_gl = .
        dbp_or = J(1, 0, 0)
        if (b.has_or) {
            bid = b.links[ell].or_block
            po = (sep ? b.blocks[bid].or_p : b.or_p)
            okeep = (sep ? b.blocks[bid].or_keep : b.or_keep)
            obeta = (sep ? b.blocks[bid].or_beta : b.or_beta)
            alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
            dbp_or = (po > 0 ? (L :* b.cell_M_z[g, okeep]) :/ d_plus : J(1, 0, 0))
            B_gl = alpha_b + (po > 0 ? dbp_or * obeta : 0)
            num_Bg = d_plus * B_gl
            b.links[ell].d_bar_plus_or = dbp_or
            b.links[ell].B_gl = B_gl
            b.links[ell].num_Bg = num_Bg
            b.links[ell].mu_reg = T - B_gl
            b.links[ell].a_coeff_Bg_g = ((num_Bg / ag) - B_gl * (d_plus / ag)) / d_plus
        }
        if (!b.has_ps) {
            mu = b.links[ell].mu_reg
            if (mu >= .) _cdid_impb_fail("Nonfinite link primary pure-imputation estimate.", 504)
            b.links[ell].mu = mu
        }
        else {
            nk = b.links[ell].n_ctrl
            la[ell].s_wm = 0
            la[ell].s_wmj = 0
            la[ell].s_wmm = 0
            la[ell].po = po
            la[ell].alpha_b = alpha_b
            la[ell].okeep = okeep
            la[ell].obeta = obeta
            la[ell].dbp_or = dbp_or
            la[ell].s_wmdor = J(1, po, 0)
            la[ell].s_wm_k = J(nk, 1, 0)
            la[ell].s_wodj_k = J(nk, 1, 0)
            la[ell].s_wod_k = J(nk, 1, 0)
            la[ell].s_wodm_k = J(nk, 1, 0)
            la[ell].s_wmdimp = J(nk, p, 0)
            la[ell].s_wdjX = J(1, q, 0)
            la[ell].s_wdX = J(1, q, 0)
            la[ell].s_wdmX = J(1, q, 0)
        }
    }
    if (!b.has_ps) return

    for (rk=1; rk<=b.n_recipes; rk++) {
        for (ck=1; ck<=b.G; ck++) {
            n_hit = 0
            for (ell=1; ell<=b.n_links; ell++) {
                if (b.links[ell].recipe != rk) continue
                k = _cdid_impb_pos(b.links[ell].ctrl_ci, ck)
                if (k > 0) n_hit++
            }
            if (n_hit == 0) continue
            hits = J(n_hit, 4, 0)
            h = 0
            for (ell=1; ell<=b.n_links; ell++) {
                if (b.links[ell].recipe != rk) continue
                k = _cdid_impb_pos(b.links[ell].ctrl_ci, ck)
                if (k == 0) continue
                h++
                tid = b.links[ell].trt_id
                hits[h, ] = (ell, k, b.links[ell].ctrl_cell[k], b.trts[tid].rep_tid)
            }
            u_ckc = uniqrows(hits[, 3])
            u_rtid = uniqrows(hits[, 4])
            for (h=1; h<=n_hit; h++) {
                hits[h, 3] = _cdid_impb_pos(u_ckc, hits[h, 3])
                hits[h, 4] = _cdid_impb_pos(u_rtid, hits[h, 4])
            }
            ak = _cdid_impb_ak(b, m, rk, ck)
            coh = _cdid_impb_crows(b, m, rk, ck)
            n_coh = rows(coh)
            slot = (rk - 1) * b.G + ck
            for (lo=1; lo<=n_coh; lo=lo+b.chunk_rows) {
                hi = min((n_coh, lo + b.chunk_rows - 1))
                rch = coh[lo..hi]
                n_r = rows(rch)
                if (lo == 1 & hi == n_coh & rcc[slot].ready) {
                    wr = rcc[slot].wr
                    v = rcc[slot].v
                    Zr = rcc[slot].Zr
                }
                else {
                    wr = _cdid_impb_w(b, m, rk, rch, 1)
                    v = wr :* ak
                    Zr = (p > 0 ? _cdid_impb_z(b, m, rk, rch, allc, 1) : J(n_r, 0, 0))
                }
                X = _cdid_impb_x(b, rk, rch, 1)
                odds_m = J(n_r, rows(u_rtid), 0)
                act_m = J(n_r, rows(u_rtid), 0)
                for (tp=1; tp<=rows(u_rtid); tp++) {
                    rtid = u_rtid[tp]
                    vals = _cdid_impb_logit_vals(X * b.trts[rtid].gamma)
                    odds_m[, tp] = vals[, 2]
                    act_m[, tp] = vals[, 3]
                }
                jh_m = J(n_r, rows(u_ckc), 0)
                if (lo == 1 & hi == n_coh) {
                    for (cp=1; cp<=rows(u_ckc); cp++) {
                        ckc = u_ckc[cp]
                        jh_m[, cp] = (cjc[ckc].ready ? cjc[ckc].jh : _cdid_impb_pred_z(m, ckc, n_r, Zr))
                    }
                }
                else {
                    for (cp=1; cp<=rows(u_ckc); cp++) {
                        jh_m[, cp] = _cdid_impb_pred_z(m, u_ckc[cp], n_r, Zr)
                    }
                }
                wm_m = v :* odds_m
                wod_m = wr :* odds_m
                wdot_m = wm_m :* act_m
                s_wm_tp = quadcolsum(wm_m)
                s_wod_tp = quadcolsum(wod_m)
                s_wdX_tp = quadcross(wdot_m, X)
                s_wmj_tc = quadcross(wm_m, jh_m)
                s_wodj_tc = quadcross(wod_m, jh_m)
                for (h=1; h<=n_hit; h++) {
                    ell = hits[h, 1]
                    k = hits[h, 2]
                    cp = hits[h, 3]
                    tp = hits[h, 4]
                    ckc = b.links[ell].ctrl_cell[k]
                    gk = m.cells[ckc].group_id
                    keptk = m.groups[gk].kept
                    pk = m.groups[gk].p
                    L = b.links[ell].L
                    jh = jh_m[, cp]
                    wm = wm_m[, tp]
                    wdot = wdot_m[, tp]
                    la[ell].s_wm = la[ell].s_wm + s_wm_tp[tp]
                    la[ell].s_wmj = la[ell].s_wmj + s_wmj_tc[tp, cp]
                    la[ell].s_wm_k[k] = la[ell].s_wm_k[k] + s_wm_tp[tp]
                    if (pk > 0)
                        la[ell].s_wmdimp[k, keptk] = la[ell].s_wmdimp[k, keptk] + quadcross(wm, L :* Zr[, keptk])
                    la[ell].s_wodj_k[k] = la[ell].s_wodj_k[k] + s_wodj_tc[tp, cp]
                    la[ell].s_wod_k[k] = la[ell].s_wod_k[k] + s_wod_tp[tp]
                    la[ell].s_wdjX = la[ell].s_wdjX + quadcross(wdot :* jh, X)
                    la[ell].s_wdX = la[ell].s_wdX + s_wdX_tp[tp, ]
                    if (b.has_or) {
                        po = la[ell].po
                        alpha_b = la[ell].alpha_b
                        if (po > 0) {
                            dor = L :* Zr[, la[ell].okeep]
                            mh = alpha_b :+ dor * la[ell].obeta
                            la[ell].s_wmdor = la[ell].s_wmdor + quadcross(wm, dor)
                        }
                        else mh = J(n_r, 1, alpha_b)
                        wod = wod_m[, tp]
                        la[ell].s_wmm = la[ell].s_wmm + quadsum(wm :* mh)
                        la[ell].s_wodm_k[k] = la[ell].s_wodm_k[k] + quadsum(wod :* mh)
                        la[ell].s_wdmX = la[ell].s_wdmX + quadcross(wdot :* mh, X)
                    }
                }
            }
        }
    }

    for (ell=1; ell<=b.n_links; ell++) {
        T = b.links[ell].T_gl
        B_gl = b.links[ell].B_gl
        d_minus = la[ell].s_wm / b.N
        if (!(d_minus > 0) | d_minus >= .)
            _cdid_impb_fail("Link " + b.links[ell].link_id + " has nonpositive control IPW denominator d_{l,-}.", 459)
        num_C = la[ell].s_wmj / b.N
        Cl = num_C / d_minus
        b.links[ell].d_minus = d_minus
        b.links[ell].num_C = num_C
        b.links[ell].C_ol = Cl
        b.links[ell].mu_ipw = T - Cl
        b.links[ell].m_minus_k = (la[ell].s_wm_k / b.N) / d_minus
        b.links[ell].d_bar_minus_j_k = (la[ell].s_wmdimp / b.N) / d_minus
        b.links[ell].a_coeff_C_k = ((la[ell].s_wodj_k / b.N) - Cl * (la[ell].s_wod_k / b.N)) / d_minus
        b.links[ell].denminus_gamma = la[ell].s_wdX / b.N
        b.links[ell].g_C_gamma = ((la[ell].s_wdjX / b.N) - Cl * (la[ell].s_wdX / b.N)) / d_minus
        if (b.has_or) {
            po = la[ell].po
            alpha_b = la[ell].alpha_b
            b.links[ell].d_bar_minus_or = (po > 0 ? (la[ell].s_wmdor / b.N) / d_minus : J(1, 0, 0))
            b.links[ell].d_bar_diff_or = la[ell].dbp_or - b.links[ell].d_bar_minus_or
            b.links[ell].num_Bo = la[ell].s_wmm / b.N
            B_ol = alpha_b + (po > 0 ? b.links[ell].d_bar_minus_or * la[ell].obeta : 0)
            b.links[ell].B_ol = B_ol
            b.links[ell].mu_dr = (T - B_gl) - (Cl - B_ol)
            b.links[ell].a_coeff_Bo_k = ((la[ell].s_wodm_k / b.N) - B_ol * (la[ell].s_wod_k / b.N)) / d_minus
            b.links[ell].g_Bo_gamma = ((la[ell].s_wdmX / b.N) - B_ol * (la[ell].s_wdX / b.N)) / d_minus
            b.links[ell].g_dr_gamma = b.links[ell].g_C_gamma - b.links[ell].g_Bo_gamma
        }
        if (missing(b.links[ell].m_minus_k) | missing(b.links[ell].d_bar_minus_j_k) |
            missing(b.links[ell].a_coeff_C_k) | missing(b.links[ell].g_C_gamma) |
            missing(b.links[ell].a_coeff_Bo_k) | missing(b.links[ell].g_Bo_gamma))
            _cdid_impb_fail("Nonfinite link control IPW/DR summary.", 504)
        mu = (b.method == "ipw" ? b.links[ell].mu_ipw : (b.method == "reg" ? b.links[ell].mu_reg : b.links[ell].mu_dr))
        if (mu >= .) _cdid_impb_fail("Nonfinite link primary pure-imputation estimate.", 504)
        b.links[ell].mu = mu
    }
}

/* OR slope j forcing adjoint: ONE shared adjoint across blocks in both pooled and separate modes. */
void _cdid_impb_or_adjoint(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m)
{
    real scalar po, tot_po, at, p, bid, k, cl, j
    real matrix lam, nu
    p = b.p_raw
    b.or_adj_ready = 0
    if (!b.has_or) return
    if (b.outcome_pooling == "separate") {
        tot_po = 0
        for (bid=1; bid<=b.n_blocks; bid++) {
            po = b.blocks[bid].or_p
            b.blocks[bid].or_adj_ready = (po > 0)
            tot_po = tot_po + po
        }
        if (tot_po == 0) return
        lam = J(b.n_cells, tot_po, 0)
        nu = J(b.n_cells, p * tot_po, 0)
        at = 0
        for (bid=1; bid<=b.n_blocks; bid++) {
            po = b.blocks[bid].or_p
            if (po == 0) continue
            for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                cl = b.blocks[bid].ctrl_cell[k]
                for (j=1; j<=po; j++) {
                    if (b.intercept) lam[cl, at + j] = lam[cl, at + j] + b.blocks[bid].jf_beta_alpha[j, k]
                    nu[cl, ((at + j - 1) * p + 1)..((at + j) * p)] = (nu[cl, ((at + j - 1) * p + 1)..((at + j) * p)] +
                        b.blocks[bid].jf_beta_beta[j, ((k - 1) * p + 1)..(k * p)])
                }
            }
            at = at + po
        }
        b.or_adj = cdid_imputation_adjoint_prepare(m, lam, nu)
        b.or_adj_ready = 1
        at = 0
        for (bid=1; bid<=b.n_blocks; bid++) {
            po = b.blocks[bid].or_p
            if (po == 0) continue
            b.blocks[bid].or_adj = cdid_imputation_adjoint_slice(m, b.or_adj, at + 1, at + po)
            at = at + po
        }
        return
    }
    po = b.or_p
    if (po == 0) return
    lam = J(b.n_cells, po, 0)
    nu = J(b.n_cells, p * po, 0)
    for (bid=1; bid<=b.n_blocks; bid++) {
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
            cl = b.blocks[bid].ctrl_cell[k]
            for (j=1; j<=po; j++) {
                if (b.intercept) lam[cl, j] = lam[cl, j] + b.blocks[bid].jf_beta_alpha[j, k]
                nu[cl, ((j - 1) * p + 1)..(j * p)] = (nu[cl, ((j - 1) * p + 1)..(j * p)] +
                    b.blocks[bid].jf_beta_beta[j, ((k - 1) * p + 1)..(k * p)])
            }
        }
    }
    b.or_adj = cdid_imputation_adjoint_prepare(m, lam, nu)
    b.or_adj_ready = 1
}

/* Local OR intercept blocks needed by a link batch (source alpha_bids, non-details). */
real colvector _cdid_impb_abids(struct cdid_imp_backend scalar b, real colvector lk)
{
    real colvector a
    real scalar i
    if (!b.has_or | !b.intercept) return(J(0, 1, 0))
    a = J(rows(lk), 1, 0)
    for (i=1; i<=rows(lk); i++) a[i] = b.links[lk[i]].or_block
    return(uniqrows(a))
}

/* One bounded adjoint per link batch on the Stage 1 machinery. Columns: T loads of
   each batch link, then C loads (ipw/dr), then OR alpha j loads of the batch blocks. */
struct cdid_imputation_adjoint scalar _cdid_impb_batch_adj(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real colvector lk, real colvector abids)
{
    real scalar n_l, M, p, mi, ell, g, gid, k, cl, a, bid, mc, ma
    real rowvector kept
    real matrix lam, nu
    n_l = rows(lk)
    p = b.p_raw
    M = n_l * (1 + b.has_ps) + rows(abids)
    lam = J(b.n_cells, M, 0)
    nu = J(b.n_cells, p * M, 0)
    for (mi=1; mi<=n_l; mi++) {
        ell = lk[mi]
        g = b.links[ell].target_cell
        gid = m.cells[g].group_id
        kept = m.groups[gid].kept
        lam[g, mi] = 1
        if (cols(kept) > 0)
            nu[g, (mi - 1) * p :+ kept] = nu[g, (mi - 1) * p :+ kept] + b.links[ell].d_bar_plus_j
        if (b.has_ps) {
            mc = n_l + mi
            for (k=1; k<=b.links[ell].n_ctrl; k++) {
                cl = b.links[ell].ctrl_cell[k]
                lam[cl, mc] = lam[cl, mc] + b.links[ell].m_minus_k[k]
                if (p > 0) {
                    nu[cl, ((mc - 1) * p + 1)..(mc * p)] = (nu[cl, ((mc - 1) * p + 1)..(mc * p)] +
                        b.links[ell].d_bar_minus_j_k[k, ])
                }
            }
        }
    }
    for (a=1; a<=rows(abids); a++) {
        bid = abids[a]
        ma = n_l * (1 + b.has_ps) + a
        for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
            cl = b.blocks[bid].ctrl_cell[k]
            lam[cl, ma] = lam[cl, ma] + b.blocks[bid].jf_alpha_alpha[k]
            if (p > 0) {
                nu[cl, ((ma - 1) * p + 1)..(ma * p)] = (nu[cl, ((ma - 1) * p + 1)..(ma * p)] +
                    b.blocks[bid].jf_alpha_beta[k, ])
            }
        }
    }
    return(cdid_imputation_adjoint_prepare(m, lam, nu))
}

/* Append OR slope j forcing adjoint columns onto the link-batch adjoint so every Stage 1 j
   pullback for a chunk executes in one combined calendar-pair-ordered sweep. */
struct cdid_imputation_adjoint scalar _cdid_impb_cat_or_adj(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, struct cdid_imputation_adjoint scalar base_adj,
    real colvector bbids)
{
    struct cdid_imputation_adjoint scalar out
    real scalar sep, M0, tot_po, at_all, M_tot, t, bid, pob, gid, at, loc, cid
    real colvector act_mask, ix_act, ac, boff_all, cids_g
    real rowvector cc, cc_src
    real matrix act_keys
    if (!b.has_or) return(base_adj)
    sep = (b.outcome_pooling == "separate")
    tot_po = 0
    at_all = 0
    boff_all = J(0, 1, 0)
    if (sep) {
        boff_all = J(b.n_blocks, 1, 0)
        for (bid=1; bid<=b.n_blocks; bid++) {
            boff_all[bid] = at_all
            at_all = at_all + b.blocks[bid].or_p
        }
        for (t=1; t<=rows(bbids); t++) {
            bid = bbids[t]
            if (b.blocks[bid].or_adj_ready) tot_po = tot_po + b.blocks[bid].or_p
        }
    }
    else {
        if (b.or_adj_ready) tot_po = b.or_adj.n_contrasts
    }
    if (tot_po == 0) return(base_adj)
    M0 = base_adj.n_contrasts
    M_tot = M0 + tot_po
    out = base_adj
    out.n_contrasts = M_tot
    out.tilde_lambda = J(m.n_cells, M_tot, 0)
    out.tilde_lambda[, 1..M0] = base_adj.tilde_lambda
    for (gid=1; gid<=m.n_groups; gid++) {
        out.groups[gid].nu = J(m.groups[gid].p, M_tot, 0)
        if (m.groups[gid].p > 0) out.groups[gid].nu[, 1..M0] = base_adj.groups[gid].nu
    }
    act_mask = J(m.n_cells, 1, 0)
    if (rows(base_adj.active_cells) > 0) act_mask[base_adj.active_cells] = J(rows(base_adj.active_cells), 1, 1)
    cc = (M0 + 1)..M_tot
    if (sep & rows(bbids) < b.n_blocks) {
        cc_src = J(1, tot_po, 0)
        at = 0
        for (t=1; t<=rows(bbids); t++) {
            bid = bbids[t]
            if (!b.blocks[bid].or_adj_ready) continue
            pob = b.blocks[bid].or_p
            cc_src[(at + 1)..(at + pob)] = (boff_all[bid] + 1)..(boff_all[bid] + pob)
            at = at + pob
        }
        out.tilde_lambda[, cc] = b.or_adj.tilde_lambda[, cc_src]
        if (m.intercept) {
            for (cid=1; cid<=m.n_cells; cid++) {
                if (any(out.tilde_lambda[cid, cc] :!= 0)) act_mask[cid] = 1
            }
        }
        for (gid=1; gid<=m.n_groups; gid++) {
            if (m.groups[gid].p > 0) {
                out.groups[gid].nu[, cc] = b.or_adj.groups[gid].nu[, cc_src]
                if (any(out.groups[gid].nu[, cc] :!= 0)) {
                    cids_g = m.groups[gid].cell_ids
                    act_mask[cids_g] = J(rows(cids_g), 1, 1)
                }
            }
        }
    }
    else {
        out.tilde_lambda[, cc] = b.or_adj.tilde_lambda
        for (gid=1; gid<=m.n_groups; gid++) {
            if (m.groups[gid].p > 0) out.groups[gid].nu[, cc] = b.or_adj.groups[gid].nu
        }
        ac = b.or_adj.active_cells
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

/* Common OR slope IF rows: N * (R'R)^-1 rhs with the raw shared-OR factor (N exactly once). */
real matrix _cdid_impb_or_solve(struct cdid_imp_backend scalar b, real matrix rhs)
{
    if (b.or_p == 0) return(J(rows(rhs), 0, 0))
    return(_cdid_nuisance_gram(b.or_R, b.or_pivot, b.N :* rhs')')
}

/* Stage 5 kernel for one bounded row chunk and one link batch (frame scale). Recurring
   contributions of a unit across recipes/blocks/cells are summed into its single original-F
   row before any product. Only used (recipe, cohort) slices are read in dated mode. */
void _cdid_impb_eval(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real colvector fr, real colvector lk, real colvector abids,
    struct cdid_imputation_adjoint scalar adj, real scalar detail,
    real matrix o_phi, real matrix o_ipw, real matrix o_reg, real matrix o_dr,
    real matrix o_ci, real matrix o_cr, real matrix o_cd)
{
    struct _cdid_impb_ivec rowvector c_ix
    real scalar B, nl, G, nr, N, p, q, po, pob, tot_po, nA, r, k, mi, ell, t, tid, rtid, tp, bid, bpos, apos, L, alpha_b
    real scalar cl, ck, cg, g, a, joff, rt, rb, rl, xcol, sep, zs
    real colvector coh, wr, Av, rep_tids, ix, v, jh, mh, rv, ixg, ixk, bids, bbids, boff, obeta
    real colvector zmark, zids, zpos
    real colvector wp, jg, dT, aT, jT, mg, dBg, aBg, oBd, oBa, oBj, rad_b, raa_b, raj_b
    real colvector wm, jk, mk, aC, aBo, dC, jC, tCd, tCa, dBo, tBd, tBa, oOd, oOa, oOj
    real rowvector cc, okeep, bcols, allc
    real matrix xa, va, Z_rc, X_rc, Jc, wm_cache, XGd, XGa, sd, sa, X, vals, Jb, dor, dt
    real matrix rbd, rba, rbj, rad, raa, raj, xbd, xba, xbj, xad, xaa, xaj
    real matrix rbd_b, rba_b, rbj_b, xbd_b, xba_b, xbj_b
    real matrix ci_, cr_, cd_, prim

    B = rows(fr)
    nl = rows(lk)
    G = b.G
    nr = b.n_recipes
    N = b.N
    p = b.p_raw
    q = b.q
    po = b.or_p
    sep = (b.outcome_pooling == "separate")
    allc = _cdid_impb_seq(1, p)
    coh = m.unit_cohort_index[fr]
    c_ix = _cdid_impb_bucket_coh(coh, G)
    Av = m.pop.A[fr]

    /* Repeat-rate IF xi_a_{r, k} and cached (W a_k, Z, X) slices for every used (recipe r, cohort k) on F. */
    xa = J(B, nr * G, 0)
    va = J(B, nr * G, 0)
    Z_rc = J(B, nr * p, 0)
    X_rc = (b.has_ps ? J(B, nr * q, 0) : J(B, 0, 0))
    for (r=1; r<=nr; r++) {
        for (k=1; k<=G; k++) {
            if (b.is_dated) {
                if (!(*(b.rdr)).used_cohorts[r, k]) continue
            }
            if (!(b.mass_F_mat[r, k] > 0)) continue
            ix = _cdid_impb_cix(c_ix, m, fr, r, k)
            if (rows(ix) == 0) continue
            wr = _cdid_impb_w(b, m, r, fr[ix], 1)
            xa[ix, (r - 1) * G + k] = wr :* (Av[ix] :- b.a_k_mat[r, k]) :/ b.mass_F_mat[r, k]
            va[ix, (r - 1) * G + k] = wr :* b.a_k_mat[r, k]
            if (p > 0) Z_rc[ix, ((r - 1) * p + 1)..(r * p)] = _cdid_impb_z(b, m, r, fr[ix], allc, 1)
            if (b.has_ps & q > 0) X_rc[ix, ((r - 1) * q + 1)..(r * q)] = _cdid_impb_x(b, r, fr[ix], 1)
        }
    }

    /* Precompute Stage 1 cell predictions once per chunk across all needed link and OR cells. */
    bbids = J(0, 1, 0)
    bids = J(0, 1, 0)
    zmark = J(b.n_cells, 1, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        zmark[b.links[ell].target_cell] = 1
        if (b.has_ps) zmark[b.links[ell].ctrl_cell] = J(b.links[ell].n_ctrl, 1, 1)
    }
    if (b.has_or) {
        if (sep) {
            bbids = J(nl, 1, 0)
            for (mi=1; mi<=nl; mi++) bbids[mi] = b.links[lk[mi]].or_block
            bbids = uniqrows(bbids)
            for (t=1; t<=rows(bbids); t++) {
                bid = bbids[t]
                zmark[b.blocks[bid].ctrl_cell] = J(b.blocks[bid].n_ctrl, 1, 1)
            }
        }
        else {
            bids = (po > 0 ? (1::b.n_blocks) : abids)
            for (t=1; t<=rows(bids); t++) {
                bid = bids[t]
                zmark[b.blocks[bid].ctrl_cell] = J(b.blocks[bid].n_ctrl, 1, 1)
            }
        }
    }
    zids = _cdid_impb_ix(zmark :== 1)
    zpos = J(b.n_cells, 1, 0)
    Jc = J(B, rows(zids), 0)
    for (zs=1; zs<=rows(zids); zs++) {
        cl = zids[zs]
        zpos[cl] = zs
        ck = m.cells[cl].cohort_index
        r = m.cells[cl].recipe
        ix = _cdid_impb_cix(c_ix, m, fr, r, ck)
        if (rows(ix) > 0) {
            Jc[ix, zs] = _cdid_impb_pred_z(m, cl, rows(ix),
                (p > 0 ? Z_rc[ix, ((r - 1) * p + 1)..(r * p)] : J(rows(ix), 0, 0)))
        }
    }

    /* Stage 1 j pullbacks of the batch functionals (T, C, OR alpha j, and OR beta j loads). */
    Jb = _cdid_impb_adj_chunk(b, m, adj, fr)

    /* Treatment gamma IF: direct (v (D - p) X) and a (xi_a a_forcing') channels; deduplicated by rep_tid. */
    rep_tids = J(0, 1, 0)
    if (b.has_ps) {
        rep_tids = J(nl, 1, 0)
        for (mi=1; mi<=nl; mi++) rep_tids[mi] = b.trts[b.links[lk[mi]].trt_id].rep_tid
        rep_tids = uniqrows(rep_tids)
    }
    XGd = J(B, q * rows(rep_tids), 0)
    XGa = J(B, q * rows(rep_tids), 0)
    wm_cache = J(B, rows(rep_tids), 0)
    for (t=1; t<=rows(rep_tids); t++) {
        rtid = rep_tids[t]
        rt = b.trts[rtid].recipe
        cc = ((t - 1) * q + 1)..(t * q)
        sd = J(B, q, 0)
        for (k=1; k<=b.trts[rtid].n_coh; k++) {
            ck = b.trts[rtid].all_ci[k]
            ix = _cdid_impb_cix(c_ix, m, fr, rt, ck)
            if (rows(ix) == 0) continue
            v = va[ix, (rt - 1) * G + ck]
            X = X_rc[ix, ((rt - 1) * q + 1)..(rt * q)]
            vals = _cdid_impb_logit_vals(X * b.trts[rtid].gamma)
            if (ck == b.trts[rtid].target_ci) {
                sd[ix, ] = (v :* (1 :- vals[, 1])) :* X
            }
            else {
                sd[ix, ] = (v :* (0 :- vals[, 1])) :* X
                wm_cache[ix, t] = v :* vals[, 2]
            }
        }
        sa = xa[, (rt - 1) * G :+ (b.trts[rtid].all_ci')] * b.trts[rtid].a_forcing'
        XGd[, cc] = _cdid_nuisance_gram(b.trts[rtid].R, b.trts[rtid].pivot, N :* sd')'
        XGa[, cc] = _cdid_nuisance_gram(b.trts[rtid].R, b.trts[rtid].pivot, N :* sa')'
    }

    /* Counterfactual OR IF: direct, a and j channels (pooled common slope or separate per-block). */
    nA = rows(abids)
    boff = J(0, 1, 0)
    xad = J(B, nA, 0)
    xaa = J(B, nA, 0)
    xaj = J(B, nA, 0)
    xbd = J(B, po, 0)
    xba = J(B, po, 0)
    xbj = J(B, po, 0)
    if (b.has_or) {
        joff = nl * (1 + b.has_ps)
        if (sep) {
            boff = J(rows(bbids), 1, 0)
            tot_po = 0
            for (t=1; t<=rows(bbids); t++) {
                boff[t] = tot_po
                tot_po = tot_po + b.blocks[bbids[t]].or_p
            }
            xbd = J(B, tot_po, 0)
            xba = J(B, tot_po, 0)
            xbj = J(B, tot_po, 0)
            for (t=1; t<=rows(bbids); t++) {
                bid = bbids[t]
                pob = b.blocks[bid].or_p
                okeep = b.blocks[bid].or_keep
                obeta = b.blocks[bid].or_beta
                rb = b.blocks[bid].recipe
                apos = _cdid_impb_pos(abids, bid)
                L = b.blocks[bid].L
                alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
                rad_b = J(B, 1, 0)
                raa_b = J(B, 1, 0)
                rbd_b = J(B, pob, 0)
                rba_b = J(B, pob, 0)
                for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                    cl = b.blocks[bid].ctrl_cell[k]
                    ck = b.blocks[bid].ctrl_ci[k]
                    xcol = (rb - 1) * G + ck
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) > 0) {
                        v = va[ix, xcol]
                        jh = Jc[ix, zpos[cl]]
                        if (pob > 0) {
                            dor = L :* Z_rc[ix, (rb - 1) * p :+ okeep]
                            mh = alpha_b :+ dor * obeta
                        }
                        else mh = J(rows(ix), 1, alpha_b)
                        rv = v :* (jh - mh)
                        if (apos > 0) rad_b[ix] = rad_b[ix] + rv
                        if (pob > 0) {
                            dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                            rbd_b[ix, ] = rbd_b[ix, ] + rv :* dt
                        }
                    }
                    if (apos > 0) raa_b = raa_b + b.blocks[bid].a_force_alpha[k] :* xa[, xcol]
                    if (pob > 0) rba_b = rba_b + xa[, xcol] * b.blocks[bid].a_force_beta[, k]'
                }
                if (b.intercept & apos > 0) rad_b = rad_b :- b.blocks[bid].score_mean_a
                if (pob > 0) {
                    rbd_b = rbd_b :- (b.blocks[bid].score_mean_b')
                    bcols = (boff[t] + 1)..(boff[t] + pob)
                    rbj_b = Jb[, (joff + nA) :+ bcols]
                    xbd_b = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, N :* rbd_b')'
                    xba_b = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, N :* rba_b')'
                    xbj_b = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, N :* rbj_b')'
                    xbd[, bcols] = xbd_b
                    xba[, bcols] = xba_b
                    xbj[, bcols] = xbj_b
                }
                if (apos > 0) {
                    raj_b = Jb[, joff + apos]
                    xad[, apos] = rad_b :/ b.blocks[bid].mass
                    xaa[, apos] = raa_b :/ b.blocks[bid].mass
                    xaj[, apos] = raj_b :/ b.blocks[bid].mass
                    if (pob > 0) {
                        xad[, apos] = xad[, apos] - xbd_b * b.blocks[bid].mean_d'
                        xaa[, apos] = xaa[, apos] - xba_b * b.blocks[bid].mean_d'
                        xaj[, apos] = xaj[, apos] - xbj_b * b.blocks[bid].mean_d'
                    }
                }
            }
        }
        else {
            rbd = J(B, po, 0)
            rba = J(B, po, 0)
            rbj = J(B, po, 0)
            rad = J(B, nA, 0)
            raa = J(B, nA, 0)
            raj = J(B, nA, 0)
            for (t=1; t<=rows(bids); t++) {
                bid = bids[t]
                rb = b.blocks[bid].recipe
                apos = _cdid_impb_pos(abids, bid)
                L = b.blocks[bid].L
                alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
                for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                    cl = b.blocks[bid].ctrl_cell[k]
                    ck = b.blocks[bid].ctrl_ci[k]
                    xcol = (rb - 1) * G + ck
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) > 0) {
                        v = va[ix, xcol]
                        jh = Jc[ix, zpos[cl]]
                        if (po > 0) {
                            dor = L :* Z_rc[ix, (rb - 1) * p :+ b.or_keep]
                            mh = alpha_b :+ dor * b.or_beta
                        }
                        else mh = J(rows(ix), 1, alpha_b)
                        rv = v :* (jh - mh)
                        if (apos > 0) rad[ix, apos] = rad[ix, apos] + rv
                        if (po > 0) {
                            dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                            rbd[ix, ] = rbd[ix, ] + rv :* dt
                        }
                    }
                    if (apos > 0) raa[, apos] = raa[, apos] + b.blocks[bid].a_force_alpha[k] :* xa[, xcol]
                    if (po > 0) rba = rba + xa[, xcol] * b.blocks[bid].a_force_beta[, k]'
                }
            }
            if (po > 0) {
                rbj = Jb[, (joff + nA + 1)..(joff + nA + po)]
                xbd = _cdid_impb_or_solve(b, rbd)
                xba = _cdid_impb_or_solve(b, rba)
                xbj = _cdid_impb_or_solve(b, rbj)
            }
            for (a=1; a<=nA; a++) {
                bid = abids[a]
                raj[, a] = Jb[, joff + a]
                xad[, a] = rad[, a] :/ b.blocks[bid].mass
                xaa[, a] = raa[, a] :/ b.blocks[bid].mass
                xaj[, a] = raj[, a] :/ b.blocks[bid].mass
                if (po > 0) {
                    xad[, a] = xad[, a] - xbd * b.blocks[bid].mean_d'
                    xaa[, a] = xaa[, a] - xba * b.blocks[bid].mean_d'
                    xaj[, a] = xaj[, a] - xbj * b.blocks[bid].mean_d'
                }
            }
        }
    }

    o_phi = J(B, nl, 0)
    if (detail) {
        o_ipw = J(B, nl, .)
        o_reg = J(B, nl, .)
        o_dr = J(B, nl, .)
        o_ci = J(B, 8 * nl, .)
        o_cr = J(B, 8 * nl, .)
        o_cd = J(B, 8 * nl, .)
    }
    pob = 0
    okeep = J(1, 0, 0)
    obeta = J(0, 1, 0)
    bcols = J(1, 0, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        rl = b.links[ell].recipe
        g = b.links[ell].target_cell
        cg = b.links[ell].target_ci
        L = b.links[ell].L
        apos = 0
        alpha_b = 0
        ixg = _cdid_impb_cix(c_ix, m, fr, rl, cg)
        wp = J(B, 1, 0)
        jg = Jc[, zpos[g]]
        if (rows(ixg) > 0) wp[ixg] = va[ixg, (rl - 1) * G + cg]
        dT = wp :* (jg :- b.links[ell].T_gl) :/ b.links[ell].d_plus
        aT = b.links[ell].a_coeff_T_g :* xa[, (rl - 1) * G + cg]
        jT = Jb[, mi]
        if (b.has_or) {
            bid = b.links[ell].or_block
            apos = _cdid_impb_pos(abids, bid)
            pob = (sep ? b.blocks[bid].or_p : po)
            okeep = (sep ? b.blocks[bid].or_keep : b.or_keep)
            obeta = (sep ? b.blocks[bid].or_beta : b.or_beta)
            if (pob > 0) {
                if (sep) {
                    bpos = _cdid_impb_pos(bbids, bid)
                    bcols = (boff[bpos] + 1)..(boff[bpos] + pob)
                }
                else bcols = 1..pob
            }
            else bcols = J(1, 0, 0)
            alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
            mg = J(B, 1, 0)
            if (rows(ixg) > 0) {
                if (pob > 0) mg[ixg] = alpha_b :+ (L :* Z_rc[ixg, (rl - 1) * p :+ okeep]) * obeta
                else mg[ixg] = J(rows(ixg), 1, alpha_b)
            }
            dBg = wp :* (mg :- b.links[ell].B_gl) :/ b.links[ell].d_plus
            aBg = b.links[ell].a_coeff_Bg_g :* xa[, (rl - 1) * G + cg]
            oBd = J(B, 1, 0)
            oBa = J(B, 1, 0)
            oBj = J(B, 1, 0)
            if (apos > 0) {
                oBd = oBd + xad[, apos]
                oBa = oBa + xaa[, apos]
                oBj = oBj + xaj[, apos]
            }
            if (pob > 0) {
                oBd = oBd + xbd[, bcols] * b.links[ell].d_bar_plus_or'
                oBa = oBa + xba[, bcols] * b.links[ell].d_bar_plus_or'
                oBj = oBj + xbj[, bcols] * b.links[ell].d_bar_plus_or'
            }
        }
        if (b.has_ps) {
            tid = b.links[ell].trt_id
            tp = _cdid_impb_pos(rep_tids, b.trts[tid].rep_tid)
            cc = ((tp - 1) * q + 1)..(tp * q)
            wm = wm_cache[, tp]
            jk = J(B, 1, 0)
            mk = J(B, 1, 0)
            aC = J(B, 1, 0)
            aBo = J(B, 1, 0)
            for (k=1; k<=b.links[ell].n_ctrl; k++) {
                ck = b.links[ell].ctrl_ci[k]
                cl = b.links[ell].ctrl_cell[k]
                xcol = (rl - 1) * G + ck
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) > 0) {
                    jk[ixk] = Jc[ixk, zpos[cl]]
                    if (b.has_or) {
                        if (pob > 0) mk[ixk] = alpha_b :+ (L :* Z_rc[ixk, (rl - 1) * p :+ okeep]) * obeta
                        else mk[ixk] = J(rows(ixk), 1, alpha_b)
                    }
                }
                aC = aC + b.links[ell].a_coeff_C_k[k] :* xa[, xcol]
                if (b.has_or) aBo = aBo + b.links[ell].a_coeff_Bo_k[k] :* xa[, xcol]
            }
            dC = wm :* (jk :- b.links[ell].C_ol) :/ b.links[ell].d_minus
            jC = Jb[, nl + mi]
            tCd = XGd[, cc] * b.links[ell].g_C_gamma'
            tCa = XGa[, cc] * b.links[ell].g_C_gamma'
            if (b.has_or) {
                dBo = wm :* (mk :- b.links[ell].B_ol) :/ b.links[ell].d_minus
                tBd = XGd[, cc] * b.links[ell].g_Bo_gamma'
                tBa = XGa[, cc] * b.links[ell].g_Bo_gamma'
                oOd = J(B, 1, 0)
                oOa = J(B, 1, 0)
                oOj = J(B, 1, 0)
                if (apos > 0) {
                    oOd = oOd + xad[, apos]
                    oOa = oOa + xaa[, apos]
                    oOj = oOj + xaj[, apos]
                }
                if (pob > 0) {
                    oOd = oOd + xbd[, bcols] * b.links[ell].d_bar_minus_or'
                    oOa = oOa + xba[, bcols] * b.links[ell].d_bar_minus_or'
                    oOj = oOj + xbj[, bcols] * b.links[ell].d_bar_minus_or'
                }
            }
        }
        ci_ = J(B, 8, .)
        cr_ = J(B, 8, .)
        cd_ = J(B, 8, .)
        if (b.has_ps) ci_ = (dT - dC, aT - aC, jT - jC, -tCd, -tCa, J(B, 3, 0))
        if (b.has_or) cr_ = (dT - dBg, aT - aBg, jT, J(B, 2, 0), -oBd, -oBa, -oBj)
        if (b.has_ps & b.has_or)
            cd_ = ((dT - dBg) - (dC - dBo), (aT - aBg) - (aC - aBo), jT - jC,
                -(tCd - tBd), -(tCa - tBa), -oBd + oOd, -oBa + oOa, -oBj + oOj)
        prim = (b.method == "ipw" ? ci_ : (b.method == "reg" ? cr_ : cd_))
        o_phi[, mi] = quadrowsum(prim)
        if (detail) {
            if (b.has_ps) {
                o_ipw[, mi] = quadrowsum(ci_)
                o_ci[, ((mi - 1) * 8 + 1)..(mi * 8)] = ci_
            }
            if (b.has_or) {
                o_reg[, mi] = quadrowsum(cr_)
                o_cr[, ((mi - 1) * 8 + 1)..(mi * 8)] = cr_
            }
            if (b.has_ps & b.has_or) {
                o_dr[, mi] = quadrowsum(cd_)
                o_cd[, ((mi - 1) * 8 + 1)..(mi * 8)] = cd_
            }
        }
    }
    if (missing(o_phi)) _cdid_impb_fail("Nonfinite complete pure-imputation link IF.", 504)
}

void _cdid_impb_set_ids(struct cdid_imp_backend_if scalar res,
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

/* Replay scratch (bounded chunk and batch) plus requested output, charged with the
   concurrently live retained stores; the opt-in work cap applies to each replay call. */
void _cdid_impb_replay_budget(struct cdid_imp_backend scalar b, real scalar Bc,
    real scalar nl, real scalar B, real scalar detail)
{
    real scalar nlb, M_b, scr, work, xa_cols
    nlb = min((nl, b.link_batch))
    M_b = nlb * (1 + b.has_ps) + min((nlb, b.n_blocks))
    xa_cols = b.n_recipes * b.G
    scr = (8 * (Bc * (xa_cols + 4 * b.q * nlb + 9 * b.p_raw + 3 * nlb + M_b + 48) +
        b.n_cells * M_b * (1 + b.p_raw) + b.n_groups * b.p_raw * M_b +
        B * nl * (1 + 27 * detail)) + 4096)
    work = (16 * B * (nl * (xa_cols + 3) * (b.p_raw + b.q + 4) + b.n_cells * M_b * (b.p_raw + 1)) *
        ceil(nl / nlb))
    _cdid_imp_check_budget(b.budget, b.retained_bytes, scr, work, "imputation backend replay")
}

/* True when evaluating all nl links in one detail=0 pass fits within max_workspace_bytes. */
real scalar _cdid_impb_can_full_batch(struct cdid_imp_backend scalar b, real scalar Bc,
    real scalar nl, real scalar B)
{
    real scalar M_full, xa_cols, scr
    if (nl <= b.link_batch) return(1)
    M_full = nl * (1 + b.has_ps) + min((nl, b.n_blocks))
    xa_cols = b.n_recipes * b.G
    scr = (8 * (Bc * (xa_cols + 4 * b.q * nl + 9 * b.p_raw + 3 * nl + M_full + 48) +
        b.n_cells * M_full * (1 + b.p_raw) + b.n_groups * b.p_raw * M_full +
        B * nl) + 4096)
    return(b.retained_bytes + scr <= b.budget.max_workspace_bytes)
}

/* Pre-build Stage 5 (detail = 0) adjoint, cell/propensity index maps, and BLAS loading matrices. */
struct _cdid_impb_plan scalar _cdid_impb_make_plan(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real colvector lk)
{
    struct _cdid_impb_plan scalar plan
    struct cdid_imputation_adjoint scalar base_adj
    real scalar nl, G, nr, p, q, po, pob, tot_po, mi, ell, t, tb, bid, bpos, apos, k, ck, cg, rl, rb, xcol, n_rep, s
    real scalar c_tgt, c_ctl
    real colvector zmark, or_mark, pos_v
    real rowvector bcols

    nl = rows(lk)
    G = b.G
    nr = b.n_recipes
    p = b.p_raw
    q = b.q
    po = b.or_p
    plan.lk = lk
    plan.nl = nl
    plan.sep = (b.outcome_pooling == "separate")
    plan.abids = (b.has_or & b.intercept & b.method == "reg" ? _cdid_impb_abids(b, lk) : J(0, 1, 0))
    plan.nA = rows(plan.abids)
    plan.bbids = J(0, 1, 0)
    plan.bids = J(0, 1, 0)
    plan.boff = J(0, 1, 0)
    tot_po = 0
    if (b.has_or) {
        if (plan.sep) {
            plan.bbids = J(nl, 1, 0)
            for (mi=1; mi<=nl; mi++) plan.bbids[mi] = b.links[lk[mi]].or_block
            plan.bbids = uniqrows(plan.bbids)
            plan.boff = J(rows(plan.bbids), 1, 0)
            for (t=1; t<=rows(plan.bbids); t++) {
                plan.boff[t] = tot_po
                tot_po = tot_po + b.blocks[plan.bbids[t]].or_p
            }
        }
        else {
            plan.bids = (po > 0 ? (1::b.n_blocks) : plan.abids)
            tot_po = po
        }
    }
    plan.tot_po = tot_po
    base_adj = _cdid_impb_batch_adj(b, m, lk, plan.abids)
    plan.adj = _cdid_impb_cat_or_adj(b, m, base_adj, plan.bbids)

    zmark = J(b.n_cells, 1, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        zmark[b.links[ell].target_cell] = 1
        if (b.has_ps) zmark[b.links[ell].ctrl_cell] = J(b.links[ell].n_ctrl, 1, 1)
    }
    if (b.has_or) {
        if (plan.sep) {
            for (t=1; t<=rows(plan.bbids); t++) {
                bid = plan.bbids[t]
                zmark[b.blocks[bid].ctrl_cell] = J(b.blocks[bid].n_ctrl, 1, 1)
            }
        }
        else {
            for (t=1; t<=rows(plan.bids); t++) {
                bid = plan.bids[t]
                zmark[b.blocks[bid].ctrl_cell] = J(b.blocks[bid].n_ctrl, 1, 1)
            }
        }
    }
    plan.zids = _cdid_impb_ix(zmark :== 1)
    plan.zpos = J(b.n_cells, 1, 0)
    for (t=1; t<=rows(plan.zids); t++) plan.zpos[plan.zids[t]] = t

    plan.rep_tids = J(0, 1, 0)
    plan.link_tp = J(nl, 1, 0)
    plan.tp_grps = _cdid_impb_tp_grp(1, 0)
    if (b.has_ps) {
        plan.rep_tids = J(nl, 1, 0)
        for (mi=1; mi<=nl; mi++) plan.rep_tids[mi] = b.trts[b.links[lk[mi]].trt_id].rep_tid
        plan.rep_tids = uniqrows(plan.rep_tids)
        n_rep = rows(plan.rep_tids)
        for (mi=1; mi<=nl; mi++) {
            plan.link_tp[mi] = _cdid_impb_pos(plan.rep_tids, b.trts[b.links[lk[mi]].trt_id].rep_tid)
        }
        plan.tp_grps = _cdid_impb_tp_grp(1, n_rep)
        for (t=1; t<=n_rep; t++) {
            pos_v = _cdid_impb_ix(plan.link_tp :== t)
            plan.tp_grps[t].rtid = plan.rep_tids[t]
            plan.tp_grps[t].link_pos = pos_v'
            plan.tp_grps[t].G_gamma = J(q, rows(pos_v), 0)
            for (s=1; s<=rows(pos_v); s++) {
                ell = lk[pos_v[s]]
                plan.tp_grps[t].G_gamma[, s] = (b.method == "ipw" ?
                    b.links[ell].g_C_gamma' : b.links[ell].g_dr_gamma')
            }
        }
    }

    plan.A_load = J(nr * G, nl, 0)
    for (mi=1; mi<=nl; mi++) {
        ell = lk[mi]
        rl = b.links[ell].recipe
        cg = b.links[ell].target_ci
        c_tgt = b.links[ell].a_coeff_T_g - (b.has_or ? b.links[ell].a_coeff_Bg_g : 0)
        plan.A_load[(rl - 1) * G + cg, mi] = plan.A_load[(rl - 1) * G + cg, mi] + c_tgt
        if (b.has_ps) {
            for (k=1; k<=b.links[ell].n_ctrl; k++) {
                ck = b.links[ell].ctrl_ci[k]
                c_ctl = b.links[ell].a_coeff_C_k[k] - (b.has_or ? b.links[ell].a_coeff_Bo_k[k] : 0)
                plan.A_load[(rl - 1) * G + ck, mi] = plan.A_load[(rl - 1) * G + ck, mi] - c_ctl
            }
        }
    }

    plan.or_bids = J(0, 1, 0)
    plan.or_bpos = J(b.n_blocks, 1, 0)
    plan.or_coh_mask = J(0, G, 0)
    plan.Alpha_load = J(plan.nA, nl, 0)
    plan.Beta_load = J(tot_po, nl, 0)
    plan.Raa_load = J(nr * G, plan.nA, 0)
    plan.Rba_load = J(nr * G, po, 0)
    if (b.has_or) {
        or_mark = J(b.n_blocks, 1, 0)
        for (mi=1; mi<=nl; mi++) or_mark[b.links[lk[mi]].or_block] = 1
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
                if (b.blocks[bid].n_ctrl > 0)
                    plan.or_coh_mask[tb, b.blocks[bid].ctrl_ci'] = J(1, b.blocks[bid].n_ctrl, 1)
            }
        }
        else {
            for (t=1; t<=rows(plan.bids); t++) {
                bid = plan.bids[t]
                tb = plan.or_bpos[bid]
                if (b.blocks[bid].n_ctrl > 0)
                    plan.or_coh_mask[tb, b.blocks[bid].ctrl_ci'] = J(1, b.blocks[bid].n_ctrl, 1)
                rb = b.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                    xcol = (rb - 1) * G + b.blocks[bid].ctrl_ci[k]
                    if (apos > 0) plan.Raa_load[xcol, apos] = plan.Raa_load[xcol, apos] + b.blocks[bid].a_force_alpha[k]
                    if (po > 0) plan.Rba_load[xcol, ] = plan.Rba_load[xcol, ] + b.blocks[bid].a_force_beta[, k]'
                }
            }
        }
        for (mi=1; mi<=nl; mi++) {
            ell = lk[mi]
            bid = b.links[ell].or_block
            tb = plan.or_bpos[bid]
            plan.or_coh_mask[tb, b.links[ell].target_ci] = 1
            if (b.has_ps & b.links[ell].n_ctrl > 0)
                plan.or_coh_mask[tb, b.links[ell].ctrl_ci'] = J(1, b.links[ell].n_ctrl, 1)
            if (plan.nA > 0) {
                apos = _cdid_impb_pos(plan.abids, bid)
                if (apos > 0) plan.Alpha_load[apos, mi] = -1
            }
            pob = (plan.sep ? b.blocks[bid].or_p : po)
            if (pob > 0) {
                if (plan.sep) {
                    bpos = _cdid_impb_pos(plan.bbids, bid)
                    bcols = (plan.boff[bpos] + 1)..(plan.boff[bpos] + pob)
                }
                else bcols = 1..pob
                plan.Beta_load[bcols, mi] = -(b.method == "reg" ?
                    b.links[ell].d_bar_plus_or' : b.links[ell].d_bar_diff_or')
            }
        }
    }
    plan.ready = 1
    return(plan)
}

/* Batched Stage 5 kernel for detail = 0 (frame scale). Combines channel solves and contracts
   repeat-rate, propensity, Stage 1 j, and OR loadings via BLAS matrix products. */
void _cdid_impb_eval_fast(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real colvector fr, struct _cdid_impb_plan scalar plan, real matrix o_phi)
{
    struct _cdid_impb_ivec rowvector c_ix
    real scalar B, nl, G, nr, N, p, q, po, pob, tot_po, nA, r, k, mi, ell, t, tb, rtid, tp, bid, apos, L, alpha_b
    real scalar cl, ck, cg, g, a, joff, rt, rb, rl, xcol, zs, n_rep
    real colvector coh, wr, Av, ix, v, jh, rv, ixg, ixk, obeta, rad_b, raa_b, raj_b
    real rowvector cc, okeep, bcols, allc, lpos
    real matrix xa, va, Z_rc, X_rc, Jc, Mc, wm_cache, XG_tot, sd, sa, X, vals, Jb, dor, dt
    real matrix rbd, rba, rbj, rad, raa, raj, xb_tot, xa_tot, rbd_b, rba_b, rbj_b, xb_b

    B = rows(fr)
    nl = plan.nl
    G = b.G
    nr = b.n_recipes
    N = b.N
    p = b.p_raw
    q = b.q
    po = b.or_p
    nA = plan.nA
    tot_po = plan.tot_po
    allc = _cdid_impb_seq(1, p)
    coh = m.unit_cohort_index[fr]
    c_ix = _cdid_impb_bucket_coh(coh, G)
    Av = m.pop.A[fr]

    xa = J(B, nr * G, 0)
    va = J(B, nr * G, 0)
    Z_rc = J(B, nr * p, 0)
    X_rc = (b.has_ps ? J(B, nr * q, 0) : J(B, 0, 0))
    for (r=1; r<=nr; r++) {
        for (k=1; k<=G; k++) {
            if (b.is_dated) {
                if (!(*(b.rdr)).used_cohorts[r, k]) continue
            }
            if (!(b.mass_F_mat[r, k] > 0)) continue
            ix = _cdid_impb_cix(c_ix, m, fr, r, k)
            if (rows(ix) == 0) continue
            wr = _cdid_impb_w(b, m, r, fr[ix], 1)
            xa[ix, (r - 1) * G + k] = wr :* (Av[ix] :- b.a_k_mat[r, k]) :/ b.mass_F_mat[r, k]
            va[ix, (r - 1) * G + k] = wr :* b.a_k_mat[r, k]
            if (p > 0) Z_rc[ix, ((r - 1) * p + 1)..(r * p)] = _cdid_impb_z(b, m, r, fr[ix], allc, 1)
            if (b.has_ps & q > 0) X_rc[ix, ((r - 1) * q + 1)..(r * q)] = _cdid_impb_x(b, r, fr[ix], 1)
        }
    }

    Jc = J(B, rows(plan.zids), 0)
    for (zs=1; zs<=rows(plan.zids); zs++) {
        cl = plan.zids[zs]
        ck = m.cells[cl].cohort_index
        r = m.cells[cl].recipe
        ix = _cdid_impb_cix(c_ix, m, fr, r, ck)
        if (rows(ix) > 0) {
            Jc[ix, zs] = _cdid_impb_pred_z(m, cl, rows(ix),
                (p > 0 ? Z_rc[ix, ((r - 1) * p + 1)..(r * p)] : J(rows(ix), 0, 0)))
        }
    }

    Mc = J(B, rows(plan.or_bids), 0)
    if (b.has_or) {
        for (t=1; t<=rows(plan.or_bids); t++) {
            bid = plan.or_bids[t]
            rb = b.blocks[bid].recipe
            L = b.blocks[bid].L
            alpha_b = (b.intercept ? b.blocks[bid].alpha : 0)
            pob = (plan.sep ? b.blocks[bid].or_p : po)
            okeep = (plan.sep ? b.blocks[bid].or_keep : b.or_keep)
            obeta = (plan.sep ? b.blocks[bid].or_beta : b.or_beta)
            for (ck=1; ck<=G; ck++) {
                if (!plan.or_coh_mask[t, ck]) continue
                ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                if (rows(ix) == 0) continue
                if (pob > 0) Mc[ix, t] = alpha_b :+ (L :* Z_rc[ix, (rb - 1) * p :+ okeep]) * obeta
                else Mc[ix, t] = J(rows(ix), 1, alpha_b)
            }
        }
    }

    Jb = _cdid_impb_adj_chunk(b, m, plan.adj, fr)

    if (b.has_ps) o_phi = (Jb[, 1..nl] - Jb[, (nl + 1)..(2 * nl)]) + xa * plan.A_load
    else o_phi = Jb[, 1..nl] + xa * plan.A_load

    n_rep = rows(plan.rep_tids)
    wm_cache = J(B, n_rep, 0)
    if (b.has_ps) {
        XG_tot = J(B, q * n_rep, 0)
        for (t=1; t<=n_rep; t++) {
            rtid = plan.rep_tids[t]
            rt = b.trts[rtid].recipe
            cc = ((t - 1) * q + 1)..(t * q)
            sd = J(B, q, 0)
            for (k=1; k<=b.trts[rtid].n_coh; k++) {
                ck = b.trts[rtid].all_ci[k]
                ix = _cdid_impb_cix(c_ix, m, fr, rt, ck)
                if (rows(ix) == 0) continue
                v = va[ix, (rt - 1) * G + ck]
                X = X_rc[ix, ((rt - 1) * q + 1)..(rt * q)]
                vals = _cdid_impb_logit_vals(X * b.trts[rtid].gamma)
                if (ck == b.trts[rtid].target_ci) {
                    sd[ix, ] = (v :* (1 :- vals[, 1])) :* X
                }
                else {
                    sd[ix, ] = (v :* (0 :- vals[, 1])) :* X
                    wm_cache[ix, t] = v :* vals[, 2]
                }
            }
            sa = xa[, (rt - 1) * G :+ (b.trts[rtid].all_ci')] * b.trts[rtid].a_forcing'
            XG_tot[, cc] = _cdid_nuisance_gram(b.trts[rtid].R, b.trts[rtid].pivot, N :* (sd + sa)')'
            lpos = plan.tp_grps[t].link_pos
            if (cols(lpos) > 0) {
                o_phi[, lpos] = o_phi[, lpos] - XG_tot[, cc] * plan.tp_grps[t].G_gamma
            }
        }
    }

    if (b.has_or) {
        joff = nl * (1 + b.has_ps)
        xa_tot = J(B, nA, 0)
        xb_tot = J(B, tot_po, 0)
        if (plan.sep) {
            for (t=1; t<=rows(plan.bbids); t++) {
                bid = plan.bbids[t]
                tb = plan.or_bpos[bid]
                pob = b.blocks[bid].or_p
                okeep = b.blocks[bid].or_keep
                rb = b.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                L = b.blocks[bid].L
                rad_b = J(B, 1, 0)
                rbd_b = J(B, pob, 0)
                for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                    cl = b.blocks[bid].ctrl_cell[k]
                    ck = b.blocks[bid].ctrl_ci[k]
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) == 0) continue
                    v = va[ix, (rb - 1) * G + ck]
                    jh = Jc[ix, plan.zpos[cl]]
                    rv = v :* (jh - Mc[ix, tb])
                    if (apos > 0) rad_b[ix] = rad_b[ix] + rv
                    if (pob > 0) {
                        dor = L :* Z_rc[ix, (rb - 1) * p :+ okeep]
                        dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                        rbd_b[ix, ] = rbd_b[ix, ] + rv :* dt
                    }
                }
                raa_b = (apos > 0 ?
                    xa[, (rb - 1) * G :+ (b.blocks[bid].ctrl_ci')] * b.blocks[bid].a_force_alpha : J(B, 1, 0))
                rba_b = (pob > 0 ?
                    xa[, (rb - 1) * G :+ (b.blocks[bid].ctrl_ci')] * b.blocks[bid].a_force_beta' : J(B, 0, 0))
                if (b.intercept & apos > 0) rad_b = rad_b :- b.blocks[bid].score_mean_a
                xb_b = J(B, 0, 0)
                if (pob > 0) {
                    rbd_b = rbd_b :- (b.blocks[bid].score_mean_b')
                    bcols = (plan.boff[t] + 1)..(plan.boff[t] + pob)
                    rbj_b = Jb[, (joff + nA) :+ bcols]
                    xb_b = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot,
                        N :* (rbd_b + rba_b + rbj_b)')'
                    xb_tot[, bcols] = xb_b
                }
                if (apos > 0) {
                    raj_b = Jb[, joff + apos]
                    xa_tot[, apos] = (rad_b + raa_b + raj_b) :/ b.blocks[bid].mass
                    if (pob > 0) xa_tot[, apos] = xa_tot[, apos] - xb_b * b.blocks[bid].mean_d'
                }
            }
        }
        else {
            rbd = J(B, po, 0)
            rad = J(B, nA, 0)
            for (t=1; t<=rows(plan.bids); t++) {
                bid = plan.bids[t]
                tb = plan.or_bpos[bid]
                rb = b.blocks[bid].recipe
                apos = _cdid_impb_pos(plan.abids, bid)
                L = b.blocks[bid].L
                for (k=1; k<=b.blocks[bid].n_ctrl; k++) {
                    cl = b.blocks[bid].ctrl_cell[k]
                    ck = b.blocks[bid].ctrl_ci[k]
                    ix = _cdid_impb_cix(c_ix, m, fr, rb, ck)
                    if (rows(ix) == 0) continue
                    v = va[ix, (rb - 1) * G + ck]
                    jh = Jc[ix, plan.zpos[cl]]
                    rv = v :* (jh - Mc[ix, tb])
                    if (apos > 0) rad[ix, apos] = rad[ix, apos] + rv
                    if (po > 0) {
                        dor = L :* Z_rc[ix, (rb - 1) * p :+ b.or_keep]
                        dt = (b.intercept ? dor :- b.blocks[bid].mean_d : dor)
                        rbd[ix, ] = rbd[ix, ] + rv :* dt
                    }
                }
            }
            raa = (nA > 0 ? xa * plan.Raa_load : J(B, 0, 0))
            rba = (po > 0 ? xa * plan.Rba_load : J(B, 0, 0))
            if (po > 0) {
                rbj = Jb[, (joff + nA + 1)..(joff + nA + po)]
                xb_tot = _cdid_impb_or_solve(b, rbd + rba + rbj)
            }
            if (nA > 0) {
                raj = Jb[, (joff + 1)..(joff + nA)]
                for (a=1; a<=nA; a++) {
                    bid = plan.abids[a]
                    xa_tot[, a] = (rad[, a] + raa[, a] + raj[, a]) :/ b.blocks[bid].mass
                    if (po > 0) xa_tot[, a] = xa_tot[, a] - xb_tot * b.blocks[bid].mean_d'
                }
            }
        }
        if (nA > 0) o_phi = o_phi + xa_tot * plan.Alpha_load
        if (tot_po > 0) o_phi = o_phi + xb_tot * plan.Beta_load
    }

    for (mi=1; mi<=nl; mi++) {
        ell = plan.lk[mi]
        rl = b.links[ell].recipe
        g = b.links[ell].target_cell
        cg = b.links[ell].target_ci
        ixg = _cdid_impb_cix(c_ix, m, fr, rl, cg)
        if (!b.has_or) {
            if (rows(ixg) > 0) {
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ b.links[ell].d_plus) :*
                    (Jc[ixg, plan.zpos[g]] :- b.links[ell].T_gl))
            }
            tp = plan.link_tp[mi]
            for (k=1; k<=b.links[ell].n_ctrl; k++) {
                ck = b.links[ell].ctrl_ci[k]
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) == 0) continue
                cl = b.links[ell].ctrl_cell[k]
                o_phi[ixk, mi] = (o_phi[ixk, mi] -
                    (wm_cache[ixk, tp] :/ b.links[ell].d_minus) :*
                    (Jc[ixk, plan.zpos[cl]] :- b.links[ell].C_ol))
            }
        }
        else if (!b.has_ps) {
            if (rows(ixg) > 0) {
                tb = plan.or_bpos[b.links[ell].or_block]
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ b.links[ell].d_plus) :*
                    ((Jc[ixg, plan.zpos[g]] :- b.links[ell].T_gl) -
                     (Mc[ixg, tb] :- b.links[ell].B_gl)))
            }
        }
        else {
            tb = plan.or_bpos[b.links[ell].or_block]
            if (rows(ixg) > 0) {
                o_phi[ixg, mi] = (o_phi[ixg, mi] +
                    (va[ixg, (rl - 1) * G + cg] :/ b.links[ell].d_plus) :*
                    ((Jc[ixg, plan.zpos[g]] :- b.links[ell].T_gl) -
                     (Mc[ixg, tb] :- b.links[ell].B_gl)))
            }
            tp = plan.link_tp[mi]
            for (k=1; k<=b.links[ell].n_ctrl; k++) {
                ck = b.links[ell].ctrl_ci[k]
                ixk = _cdid_impb_cix(c_ix, m, fr, rl, ck)
                if (rows(ixk) == 0) continue
                cl = b.links[ell].ctrl_cell[k]
                o_phi[ixk, mi] = (o_phi[ixk, mi] -
                    (wm_cache[ixk, tp] :/ b.links[ell].d_minus) :*
                    ((Jc[ixk, plan.zpos[cl]] :- b.links[ell].C_ol) -
                     (Mc[ixk, tb] :- b.links[ell].B_ol)))
            }
        }
    }
    if (missing(o_phi)) _cdid_impb_fail("Nonfinite complete pure-imputation link IF.", 504)
}

/* Public row replay. f_rows: unique original-F positions (any order, any count; bounded
   internally by chunk_rows). links: unique backend link indices, processed in batches of
   link_batch with ONE adjoint per batch. scale frame (phi_F, normalization N) or core
   ((n_A/N) phi_F, normalization n_A). detail=1 adds psi_ipw/psi_reg/psi_dr and the 8 source
   channels for every method the backend supports; unsupported ones stay missing. */
struct cdid_imp_backend_if scalar cdid_imp_backend_rows(struct cdid_imp_backend scalar b,
    real colvector f_rows, real colvector links, | string scalar scale, real scalar detail)
{
    struct cdid_imp_backend_if scalar res
    struct _cdid_impb_plan scalar plan
    struct cdid_imputation_adjoint scalar adj
    pointer(struct cdid_imputation_model scalar) scalar pm
    real scalar B, nl, dmode, blo, bhi, lo, hi, c0, c1, mi, step
    real colvector lk, rch, abids, bbids
    real matrix t_phi, t_ipw, t_reg, t_dr, t_ci, t_cr, t_cd
    string scalar sc

    sc = (args() < 4 | scale == "" ? "frame" : scale)
    dmode = (args() < 5 ? 0 : detail)
    if (b.stage1_only) _cdid_impb_fail("First-step helper backend has no pure-imputation IF evaluation.", 498)
    if (sc != "frame" & sc != "core") _cdid_impb_fail("IF scale must be frame or core.", 498)
    if (dmode != 0 & dmode != 1) _cdid_impb_fail("detail must be 0 or 1.", 498)
    _cdid_impb_need_model(b)
    pm = b.model
    _cdid_impb_check(b, *pm)
    B = rows(f_rows)
    if (B < 1 | missing(f_rows) | any(f_rows :< 1) | any(f_rows :> b.N) | any(f_rows :!= floor(f_rows)))
        _cdid_impb_fail("Replay rows must be original-F positions 1..N.", 498)
    if (B > 1) {
        if (!all(f_rows[2..B] :> f_rows[1..(B - 1)])) {
            if (rows(uniqrows(f_rows)) != B) _cdid_impb_fail("Replay rows must not repeat an original-F unit.", 498)
        }
    }
    nl = rows(links)
    if (nl < 1 | missing(links) | any(links :< 1) | any(links :> b.n_links) | any(links :!= floor(links)))
        _cdid_impb_fail("Invalid backend link indices.", 498)
    if (nl > 1) {
        if (!all(links[2..nl] :> links[1..(nl - 1)])) {
            if (rows(uniqrows(links)) != nl) _cdid_impb_fail("A link request must not repeat a link.", 498)
        }
    }
    _cdid_impb_replay_budget(b, min((B, b.chunk_rows)), nl, B, dmode)
    res.phi = J(B, nl, 0)
    res.psi_ipw = J(0, 0, .)
    res.psi_reg = J(0, 0, .)
    res.psi_dr = J(0, 0, .)
    res.comp_ipw = J(0, 0, .)
    res.comp_reg = J(0, 0, .)
    res.comp_dr = J(0, 0, .)
    if (!dmode) {
        step = (_cdid_impb_can_full_batch(b, min((B, b.chunk_rows)), nl, B) ? nl : b.link_batch)
        for (blo=1; blo<=nl; blo=blo+step) {
            bhi = min((nl, blo + step - 1))
            lk = links[blo..bhi]
            plan = _cdid_impb_make_plan(b, *pm, lk)
            for (lo=1; lo<=B; lo=lo+b.chunk_rows) {
                hi = min((B, lo + b.chunk_rows - 1))
                rch = f_rows[lo..hi]
                _cdid_impb_eval_fast(b, *pm, rch, plan, t_phi=.)
                res.phi[lo..hi, blo..bhi] = t_phi
            }
        }
    }
    else {
        res.psi_ipw = J(B, nl, .)
        res.psi_reg = J(B, nl, .)
        res.psi_dr = J(B, nl, .)
        res.comp_ipw = J(B, 8 * nl, .)
        res.comp_reg = J(B, 8 * nl, .)
        res.comp_dr = J(B, 8 * nl, .)
        for (blo=1; blo<=nl; blo=blo+b.link_batch) {
            bhi = min((nl, blo + b.link_batch - 1))
            lk = links[blo..bhi]
            abids = _cdid_impb_abids(b, lk)
            bbids = J(0, 1, 0)
            if (b.has_or & b.outcome_pooling == "separate") {
                bbids = J(rows(lk), 1, 0)
                for (mi=1; mi<=rows(lk); mi++) bbids[mi] = b.links[lk[mi]].or_block
                bbids = uniqrows(bbids)
            }
            adj = _cdid_impb_batch_adj(b, *pm, lk, abids)
            adj = _cdid_impb_cat_or_adj(b, *pm, adj, bbids)
            c0 = (blo - 1) * 8 + 1
            c1 = bhi * 8
            for (lo=1; lo<=B; lo=lo+b.chunk_rows) {
                hi = min((B, lo + b.chunk_rows - 1))
                rch = f_rows[lo..hi]
                _cdid_impb_eval(b, *pm, rch, lk, abids, adj, dmode, t_phi=., t_ipw=., t_reg=., t_dr=.,
                    t_ci=., t_cr=., t_cd=.)
                res.phi[lo..hi, blo..bhi] = t_phi
                res.psi_ipw[lo..hi, blo..bhi] = t_ipw
                res.psi_reg[lo..hi, blo..bhi] = t_reg
                res.psi_dr[lo..hi, blo..bhi] = t_dr
                res.comp_ipw[lo..hi, c0..c1] = t_ci
                res.comp_reg[lo..hi, c0..c1] = t_cr
                res.comp_dr[lo..hi, c0..c1] = t_cd
            }
        }
    }
    if (sc == "core") {
        res.phi = b.lambda :* res.phi
        if (dmode) {
            res.psi_ipw = b.lambda :* res.psi_ipw
            res.psi_reg = b.lambda :* res.psi_reg
            res.psi_dr = b.lambda :* res.psi_dr
            res.comp_ipw = b.lambda :* res.comp_ipw
            res.comp_reg = b.lambda :* res.comp_reg
            res.comp_dr = b.lambda :* res.comp_dr
        }
    }
    res.method = b.method
    res.scale = sc
    res.f_rows = f_rows
    res.links = links
    res.N = b.N
    res.n_A = b.n_A
    res.lambda = b.lambda
    res.normalization_count = (sc == "core" ? b.n_A : b.N)
    res.detail = dmode
    _cdid_impb_set_ids(res, *pm, f_rows)
    return(res)
}

/* L x 4 point summary: primary mu, mu_ipw, mu_reg, mu_dr (missing when unavailable). */
real matrix cdid_imp_backend_points(struct cdid_imp_backend scalar b)
{
    real matrix out
    real scalar ell
    if (b.stage1_only) _cdid_impb_fail("First-step helper backend has no pure-imputation points.", 498)
    out = J(b.n_links, 4, .)
    for (ell=1; ell<=b.n_links; ell++) {
        out[ell, ] = (b.links[ell].mu, b.links[ell].mu_ipw, b.links[ell].mu_reg, b.links[ell].mu_dr)
    }
    return(out)
}

void cdid_imp_backend_describe(struct cdid_imp_backend scalar b)
{
    printf("{txt}imputation backend %s: method=%s pooling=%s dated=%g N_F=%g n_A=%g recipes=%g links=%g cells=%g groups=%g\n",
        b.backend_id, b.method, b.pooling, b.is_dated, b.N, b.n_A, b.n_recipes, b.n_links, b.n_cells, b.n_groups)
    printf("{txt}  fits: imputation groups=%g shared OR=%g (blocks=%g, p_OR=%g) treatment logits=%g observation logits=%g selection=%s\n",
        b.n_imputation_fits, b.n_or_fits, b.n_blocks, b.or_p, b.n_ps_fits, b.n_obs_logits, b.selection)
}

/* GMM weight-provider context. consumer_map is M x 2: (backend link index, or 0 for a
   reference-zero moment; sign +1/-1), in original scientific-catalog order; an empty map
   means identity over the backend links. The context borrows the backend pointer. */
struct cdid_imp_backend_wctx scalar cdid_imp_backend_wctx_make(
    pointer(struct cdid_imp_backend scalar) scalar backend,
    real matrix consumer_map,
    | string colvector consumer_ids)
{
    struct cdid_imp_backend_wctx scalar ctx
    real scalar M, j, nl
    real matrix cm
    if (backend == NULL) _cdid_impb_fail("Weight provider requires a prepared backend pointer.", 498)
    if ((*backend).stage1_only) _cdid_impb_fail("First-step helper backend cannot build a pure-imputation weight provider.", 498)
    nl = _cdid_impb_nlinks(*backend)
    cm = consumer_map
    if (rows(cm) == 0) cm = ((1::nl), J(nl, 1, 1))
    M = rows(cm)
    if (cols(cm) != 2 | missing(cm))
        _cdid_impb_fail("consumer_map must be a finite M x 2 matrix (link index or 0, sign).", 498)
    if (any(cm[, 1] :< 0) | any(cm[, 1] :> nl) | any(cm[, 1] :!= floor(cm[, 1])) | any(abs(cm[, 2]) :!= 1))
        _cdid_impb_fail("consumer_map needs link indices in 0..L (0 = reference zero) and signs +1/-1.", 498)
    if (!any(cm[, 1] :> 0)) _cdid_impb_fail("consumer_map must consume at least one backend link.", 498)
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
        _cdid_impb_fail("consumer_ids must be empty or one per consumer moment.", 498)
    ctx.row_chunk = _cdid_impb_chunk(*backend)
    ctx.fast_ready = 0
    if ((*backend).model != NULL & _cdid_impb_can_full_batch(*backend, ctx.row_chunk, rows(ctx.need_links), ctx.row_chunk)) {
        ctx.plan = _cdid_impb_make_plan(*backend, *((*backend).model), ctx.need_links)
        ctx.fast_ready = 1
    }
    return(ctx)
}

/* GMM provider weight(context, lo, hi): (hi-lo+1) x M CORE IF rows lo..hi of the original F
   roster in consumer order. Pure imputation replays ALL F rows (replay_rows = N_F) with
   normalization_count = n_A; this provider has no point-only switch. */
real matrix cdid_imp_backend_weight(pointer() scalar context, real scalar lo, real scalar hi)
{
    struct cdid_imp_backend_wctx scalar ctx
    struct cdid_imp_backend_if scalar r
    pointer(struct cdid_imp_backend scalar) scalar pb
    pointer(struct cdid_imputation_model scalar) scalar pm
    real matrix out, phi, t_phi
    real scalar j, N, B, rlo, rhi
    if (context == NULL) _cdid_impb_fail("Missing prepared imputation weight-provider context.", 498)
    ctx = *context
    pb = ctx.backend
    if (pb == NULL) _cdid_impb_fail("Weight-provider backend pointer is NULL; its lifetime ended.", 498)
    N = _cdid_impb_N(*pb)
    if (lo < 1 | hi < lo | hi > N | lo != floor(lo) | hi != floor(hi))
        _cdid_impb_fail("Weight-provider rows must be a contiguous block of original F rows.", 498)
    B = hi - lo + 1
    if (ctx.fast_ready & _cdid_impb_can_full_batch(*pb, min((B, (*pb).chunk_rows)), rows(ctx.need_links), B)) {
        _cdid_impb_need_model(*pb)
        pm = (*pb).model
        _cdid_impb_check(*pb, *pm)
        _cdid_impb_replay_budget(*pb, min((B, (*pb).chunk_rows)), rows(ctx.need_links), B, 0)
        phi = J(B, rows(ctx.need_links), 0)
        for (rlo=lo; rlo<=hi; rlo=rlo+(*pb).chunk_rows) {
            rhi = min((hi, rlo + (*pb).chunk_rows - 1))
            _cdid_impb_eval_fast(*pb, *pm, (rlo::rhi), ctx.plan, t_phi=.)
            phi[(rlo - lo + 1)..(rhi - lo + 1), ] = (*pb).lambda :* t_phi
        }
    }
    else {
        r = cdid_imp_backend_rows(*pb, (lo::hi), ctx.need_links, "core", 0)
        phi = r.phi
    }
    if (all(ctx.col_pos :> 0)) return(phi[, ctx.col_pos] :* (ctx.col_sign'))
    out = J(B, ctx.n_moments, 0)
    for (j=1; j<=ctx.n_moments; j++) {
        if (ctx.col_pos[j] > 0) out[, j] = ctx.col_sign[j] :* phi[, ctx.col_pos[j]]
    }
    return(out)
}

/* Complete-unit weight Gram through the GMM constructor: N_F replay rows,
   centered over F, divided by n_A. */
struct cdid_gmm_unit_gram scalar cdid_imp_backend_gram(pointer() scalar context,
    real scalar chunk, real scalar cap)
{
    struct cdid_imp_backend_wctx scalar ctx
    pointer(struct cdid_imp_backend scalar) scalar pb
    if (context == NULL) _cdid_impb_fail("Missing prepared imputation weight-provider context.", 498)
    ctx = *context
    pb = ctx.backend
    if (pb == NULL) _cdid_impb_fail("Weight-provider backend pointer is NULL; its lifetime ended.", 498)
    return(_cdid_gmm_unit_gram_norm(_cdid_impb_N(*pb), ctx.n_moments, _cdid_impb_nA(*pb),
        chunk, &cdid_imp_backend_weight(), context, cap))
}

/* Prepare the static pure-imputation backend on a borrowed Stage 1 fit (1 recipe). */
struct cdid_imp_backend scalar cdid_imp_backend_prepare(
    pointer(struct cdid_imputation_model scalar) scalar model,
    pointer(real matrix) scalar xtrt,
    real matrix link_spec,
    real matrix link_controls,
    string colvector link_ids,
    string scalar method,
    real scalar glm_eps,
    real scalar glm_maxit,
    real scalar link_batch,
    struct cdid_imputation_budget scalar budget,
    | string scalar backend_id,
      string scalar outcome_pooling_in,
      real scalar stage1_only)
{
    struct cdid_imp_backend scalar b
    struct _cdid_impb_rcc rowvector rcc
    struct _cdid_impb_cjc rowvector cjc
    pointer(struct cdid_imputation_model scalar) scalar pm
    real scalar tid, s1, ci, max_coh, c_bytes, use_cache
    string scalar op
    if (model == NULL)
        _cdid_impb_fail("Backend preparation requires a first-step fit pointer.", 498)
    if (method != "ipw" & method != "reg" & method != "dr")
        _cdid_impb_fail("est_method must be ipw, reg, or dr.", 498)
    if (link_batch < 1 | link_batch >= . | link_batch != floor(link_batch))
        _cdid_impb_fail("link_batch must be a positive integer.", 498)
    s1 = (args() >= 13 & stage1_only != 0 & stage1_only < . ? 1 : 0)
    if (s1 & method != "reg")
        _cdid_impb_fail("stage1_only helper mode requires method reg.", 498)
    cdid_imputation_validate_budget(budget)
    b.method = method
    b.stage1_only = s1
    b.has_or = (!s1 & (method == "reg" | method == "dr"))
    b.has_ps = (!s1 & (method == "ipw" | method == "dr"))
    if (b.has_or) {
        op = (args() >= 12 ? strlower(strtrim(outcome_pooling_in)) : "")
        if (op == "") op = "separate"
        if (op != "separate" & op != "pooled")
            _cdid_impb_fail("outcome_pooling must be separate or pooled.", 198)
        b.outcome_pooling = op
    }
    else b.outcome_pooling = ""
    b.is_dated = 0
    if (b.has_ps & (glm_eps <= 0 | glm_eps >= . | glm_maxit < 1 | glm_maxit >= . |
        glm_maxit != floor(glm_maxit)))
        _cdid_impb_fail("Treatment logits require an explicit positive GLM epsilon and integer maxit.", 498)
    pm = model
    b.model = model
    b.xtrt = NULL
    b.rdr = NULL
    if (b.has_ps) b.xtrt = xtrt
    b.budget = budget
    b.backend_id = ""
    if (args() >= 11) b.backend_id = backend_id
    b.link_batch = link_batch
    b.glm_eps = glm_eps
    b.glm_maxit = glm_maxit
    _cdid_impb_init(b, *pm)
    _cdid_impb_registry(b, *pm, link_spec, link_controls, link_ids)
    _cdid_impb_preflight(b, *pm)
    max_coh = 0
    for (ci=1; ci<=b.G; ci++) max_coh = max((max_coh, rows(b.coh_u[ci].rows)))
    c_bytes = 8 * (b.n_recipes * b.N * (b.p_raw + 2) + b.n_cells * max_coh)
    use_cache = (!s1 & (b.has_or | b.has_ps) & max_coh <= b.chunk_rows &
        b.retained_bytes + c_bytes <= b.budget.max_workspace_bytes)
    rcc = _cdid_impb_rcc(1, b.n_recipes * b.G)
    cjc = _cdid_impb_cjc(1, b.n_cells)
    for (ci=1; ci<=cols(rcc); ci++) rcc[ci].ready = 0
    for (ci=1; ci<=b.n_cells; ci++) cjc[ci].ready = 0
    _cdid_impb_cell_moments(b, *pm, rcc, cjc, use_cache)
    if (!s1) {
        if (b.has_or) {
            _cdid_impb_fit_or(b, *pm, rcc, cjc)
            _cdid_impb_or_forcing(b, *pm, rcc, cjc)
        }
        if (b.has_ps) {
            for (tid=1; tid<=b.n_trt; tid++) _cdid_impb_fit_ps(b, *pm, tid)
            b.n_ps_fits = b.n_trt
        }
        _cdid_impb_link_sums(b, *pm, rcc, cjc)
        _cdid_impb_or_adjoint(b, *pm)
    }
    b.borrow_sum = _cdid_impb_borrow_sum(*pm)
    return(b)
}

/* Prepare the dated multi-recipe pure-imputation backend on a borrowed
   Stage 1 dated fit (cdid_imp_fit_dated) and borrowed cdid_imp_dated_reader. */
struct cdid_imp_backend scalar cdid_imp_backend_prepare_d(
    pointer(struct cdid_imputation_model scalar) scalar model,
    pointer(struct cdid_imp_dated_reader scalar) scalar rdr,
    real matrix link_spec,
    real matrix link_controls,
    string colvector link_ids,
    string scalar method,
    real scalar glm_eps,
    real scalar glm_maxit,
    real scalar link_batch,
    struct cdid_imputation_budget scalar budget,
    | string scalar backend_id,
      string scalar outcome_pooling_in,
      real scalar stage1_only)
{
    struct cdid_imp_backend scalar b
    struct _cdid_impb_rcc rowvector rcc
    struct _cdid_impb_cjc rowvector cjc
    pointer(struct cdid_imputation_model scalar) scalar pm
    real scalar tid, s1, ci, max_coh, c_bytes, use_cache
    string scalar op
    if (model == NULL)
        _cdid_impb_fail("Backend preparation requires a first-step fit pointer.", 498)
    if (rdr == NULL)
        _cdid_impb_fail("Dated backend preparation requires a borrowed cdid_imp_dated_reader pointer.", 498)
    if (method != "ipw" & method != "reg" & method != "dr")
        _cdid_impb_fail("est_method must be ipw, reg, or dr.", 498)
    if (link_batch < 1 | link_batch >= . | link_batch != floor(link_batch))
        _cdid_impb_fail("link_batch must be a positive integer.", 498)
    s1 = (args() >= 13 & stage1_only != 0 & stage1_only < . ? 1 : 0)
    if (s1 & method != "reg")
        _cdid_impb_fail("stage1_only helper mode requires method reg.", 498)
    cdid_imputation_validate_budget(budget)
    b.method = method
    b.stage1_only = s1
    b.has_or = (!s1 & (method == "reg" | method == "dr"))
    b.has_ps = (!s1 & (method == "ipw" | method == "dr"))
    if (b.has_or) {
        op = (args() >= 12 ? strlower(strtrim(outcome_pooling_in)) : "")
        if (op == "") op = "separate"
        if (op != "separate" & op != "pooled")
            _cdid_impb_fail("outcome_pooling must be separate or pooled.", 198)
        b.outcome_pooling = op
    }
    else b.outcome_pooling = ""
    b.is_dated = 1
    if (b.has_ps & (glm_eps <= 0 | glm_eps >= . | glm_maxit < 1 | glm_maxit >= . |
        glm_maxit != floor(glm_maxit)))
        _cdid_impb_fail("Treatment logits require an explicit positive GLM epsilon and integer maxit.", 498)
    pm = model
    b.model = model
    b.xtrt = NULL
    b.rdr = rdr
    b.budget = budget
    b.backend_id = ""
    if (args() >= 11) b.backend_id = backend_id
    b.link_batch = link_batch
    b.glm_eps = glm_eps
    b.glm_maxit = glm_maxit
    _cdid_impb_init(b, *pm)
    _cdid_impb_registry(b, *pm, link_spec, link_controls, link_ids)
    _cdid_impb_preflight(b, *pm)
    max_coh = 0
    for (ci=1; ci<=b.G; ci++) max_coh = max((max_coh, rows(b.coh_u[ci].rows)))
    c_bytes = 8 * (b.n_recipes * b.N * (b.p_raw + 2) + b.n_cells * max_coh)
    use_cache = (!s1 & (b.has_or | b.has_ps) & max_coh <= b.chunk_rows &
        b.retained_bytes + c_bytes <= b.budget.max_workspace_bytes)
    rcc = _cdid_impb_rcc(1, b.n_recipes * b.G)
    cjc = _cdid_impb_cjc(1, b.n_cells)
    for (ci=1; ci<=cols(rcc); ci++) rcc[ci].ready = 0
    for (ci=1; ci<=b.n_cells; ci++) cjc[ci].ready = 0
    _cdid_impb_cell_moments(b, *pm, rcc, cjc, use_cache)
    if (!s1) {
        if (b.has_or) {
            _cdid_impb_fit_or(b, *pm, rcc, cjc)
            _cdid_impb_or_forcing(b, *pm, rcc, cjc)
        }
        if (b.has_ps) {
            for (tid=1; tid<=b.n_trt; tid++) _cdid_impb_fit_ps(b, *pm, tid)
            b.n_ps_fits = b.n_trt
        }
        _cdid_impb_link_sums(b, *pm, rcc, cjc)
        _cdid_impb_or_adjoint(b, *pm)
    }
    b.borrow_sum = (*rdr).borrow_sum + _cdid_impb_borrow_sum(*pm)
    return(b)
}

end
