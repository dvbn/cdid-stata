*! Stage 1 actual observed-outcome-change imputation model; StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (): R/chained_imputation.R.
*! Load cdid_population.mata, cdid_nuisance.mata, and cdid_ipw_link.mata first.
*! Stage 1 only: fixed-weight observed-pair WLS j_hat_{k,e}(X), bounded LRU pair reader,
*! canonical cell/group/recipe registry, full-F alias checks, and original-unit IF pullback.
version 19.0
mata:

struct cdid_imputation_budget {
    real scalar max_chunk_rows
    real scalar max_cached_pairs
    real scalar max_workspace_bytes
    real scalar max_work_ops
    real scalar enforce_work_ops
}

struct cdid_imputation_recipe {
    real scalar recipe_id
    string scalar recipe_key
    real colvector weight
    real matrix Z
}

struct cdid_imputation_pair {
    real scalar pair_id, source_pair_id, s_index, t_index, recipe, L
    string scalar pair_label, pair_key
}

struct _cdid_imp_cache_slot {
    real scalar pair_id, s_index, t_index, recipe, stamp, valid
    string scalar pair_key
    real colvector obs_f, dy
}

struct _cdid_imp_date_post {
    real colvector u, y
}

struct cdid_imputation_pair_reader {
    real scalar max_cached_pairs, chunk_rows, clock
    real scalar n_queries, n_hits, n_misses, n_evictions, current_cached_pairs
    struct _cdid_imp_cache_slot rowvector slots
    struct _cdid_imp_date_post rowvector date_posts
    real scalar date_ready, date_N, date_R, date_T
}

struct cdid_imputation_cell {
    real scalar cell_id, cohort_index, cohort, pair_id, recipe, s_index, t_index, L
    real scalar group_id, local_id, n_pair, weight_sum, mass, alpha, mean_y, sumY
    real scalar F_n, F_weight_sum, F_mass, F_mean_j
    string scalar name, cohort_label, pair_label
    real rowvector mean_d, sumZ, F_sum_Z, F_mean_Z, F_sum_jZ, F_mean_jZ, d_bar_j
    real rowvector xmin, xmax, firstZ, colmax_obs
    real matrix qr_factor
}

struct cdid_imputation_group {
    real scalar group_id, group_index0, recipe, n_cells, n_fit_rows, p_raw, p, intercept
    real scalar rss, objective, rank_tol, alias_threshold, rcond_R, rows_checked
    string scalar name, pooling
    real colvector cell_ids, cohorts, alpha, beta, block_mass, mean_y
    real rowvector kept, retained, proposed, dropped, candidates, pivot
    real rowvector control_scaled, treated_scaled
    real matrix mean_D, R, alias_map, alias_constants_Z, alias_constants_D
}

struct _cdid_imp_group_nu {
    real matrix nu
}

struct cdid_imputation_adjoint {
    real scalar n_cells, n_groups, n_contrasts, intercept
    real colvector active_cells
    real matrix tilde_lambda
    struct _cdid_imp_group_nu rowvector groups
}

struct cdid_imputation_contrast_if {
    real colvector f_rows, estimate
    real matrix phi, direct, j_pullback
    string scalar scale
    real scalar normalization_count
}

struct _cdid_imp_ivec {
    real colvector rows
}

struct cdid_imputation_model {
    struct cdid_population scalar pop
    struct cdid_imputation_budget scalar budget
    struct cdid_imputation_pair_reader scalar reader
    struct cdid_imputation_recipe rowvector recipes
    struct cdid_imputation_pair rowvector pairs
    struct cdid_imputation_cell rowvector cells
    struct cdid_imputation_group rowvector groups
    struct _cdid_imp_ivec rowvector coh_u
    real scalar N, n_A, lambda, n_cohorts, n_recipes, n_pairs, n_cells, n_groups
    real scalar intercept, p_raw, rank_tol
    string scalar pooling
    string rowvector slope_names
    real colvector cohorts, unit_cohort_index, active_cohorts
    string colvector cohort_labels
    real matrix cohort_mass_F, cohort_repeat_mass_F, cohort_a_F
    real matrix cohort_tau_F, cohort_eta_F
    real colvector mean_w_F, mean_wA_F
    real matrix cohort_colmax_F, cell_lookup
    real matrix elig_keep
    real scalar retained_bytes_bound, scratch_bytes_bound, work_ops_bound
}

void _cdid_imp_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_imp_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}

real colvector _cdid_imp_coh_rows(struct cdid_imputation_model scalar model, real scalar ci)
{
    if (cols(model.coh_u) == model.n_cohorts) return(model.coh_u[ci].rows)
    return(_cdid_imp_index(model.unit_cohort_index :== ci))
}

/* Rows of cohort ci read under recipe r. A recipe with a restricted reference
   eligibility context reads only the units that are eligible in that context. */
real colvector _cdid_imp_rows_rc(struct cdid_imputation_model scalar model, real scalar r, real scalar ci)
{
    real colvector rr
    rr = _cdid_imp_coh_rows(model, ci)
    if (rows(model.elig_keep) == 0 | rows(rr) == 0) return(rr)
    return(select(rr, model.elig_keep[rr, r] :!= 0))
}

/* Positions ix into the chunk f_rows whose units recipe r reads. */
real colvector _cdid_imp_keep_ix(struct cdid_imputation_model scalar model, real scalar r,
    real colvector f_rows, real colvector ix)
{
    if (rows(model.elig_keep) == 0 | rows(ix) == 0) return(ix)
    return(select(ix, model.elig_keep[f_rows[ix], r] :!= 0))
}

struct _cdid_imp_ivec rowvector _cdid_imp_bucket_coh(real colvector coh, real scalar G)
{
    struct _cdid_imp_ivec rowvector out
    real colvector uc
    real scalar B, c, j
    out = _cdid_imp_ivec(1, G)
    for (c=1; c<=G; c++) out[c].rows = J(0, 1, 0)
    B = rows(coh)
    if (B == 0) return(out)
    c = coh[1]
    if (all(coh :== c)) {
        out[c].rows = (1::B)
        return(out)
    }
    if (G <= 16) {
        for (c=1; c<=G; c++) out[c].rows = _cdid_imp_index(coh :== c)
        return(out)
    }
    uc = uniqrows(coh)
    for (j=1; j<=rows(uc); j++) {
        c = uc[j]
        out[c].rows = _cdid_imp_index(coh :== c)
    }
    return(out)
}

real scalar _cdid_imp_lex(real colvector a, real colvector b)
{
    real scalar i
    for (i=1; i<=min((rows(a), rows(b))); i++) {
        if (a[i] < b[i]) return(-1)
        if (a[i] > b[i]) return(1)
    }
    return(sign(rows(a) - rows(b)))
}

real scalar _cdid_imp_lower_bound(real colvector sorted_vec, real scalar target)
{
    real scalar lo, hi, mid
    lo = 1
    hi = rows(sorted_vec)
    while (lo <= hi) {
        mid = floor((lo + hi) / 2)
        if (sorted_vec[mid] < target) lo = mid + 1
        else hi = mid - 1
    }
    return(lo)
}

real scalar _cdid_imp_find_exact(real colvector sorted_vec, real scalar target)
{
    real scalar pos
    if (rows(sorted_vec) == 0) return(0)
    pos = _cdid_imp_lower_bound(sorted_vec, target)
    if (pos > rows(sorted_vec)) return(0)
    if (sorted_vec[pos] != target) return(0)
    return(pos)
}

string scalar _cdid_imp_cohort_label(real scalar k)
{
    if (k == 0 | k >= .) return("0")
    return(strofreal(k, "%21.0g"))
}

/* Bounded augmented QR compression in original column order, including response.
   Reuses the _cdid_or_compress convention from cdid_or_dr_link.mata. */
real matrix _cdid_imp_compress(real matrix A)
{
    real matrix Q, R
    real rowvector tau, pivot
    real scalar k
    k = cols(A)
    if (rows(A) == 0) return(J(0, k, 0))
    if (rows(A) < k) A = A \ J(k - rows(A), k, 0)
    if (missing(A)) _cdid_imp_fail("Nonfinite weighted cell rows or QR factor.", 504)
    hqrdp(A, Q=., tau=., R=., pivot=.)
    if (missing(R)) _cdid_imp_fail("Nonfinite augmented cell QR factor.", 504)
    return(R[, invorder(pivot)])
}

/* Caller-controlled workspace and chunk budget.
   Default 1 GiB workspace is caller policy, not a hardcoded fixed kernel refusal.
   Work-ops cap is optional and enforced only when explicitly supplied (< .). */
struct cdid_imputation_budget scalar cdid_imputation_default_budget(
    real scalar max_chunk_rows,
    real scalar max_cached_pairs,
    | real scalar max_workspace_bytes,
      real scalar max_work_ops)
{
    struct cdid_imputation_budget scalar b
    b.max_chunk_rows = max_chunk_rows
    b.max_cached_pairs = max_cached_pairs
    b.max_workspace_bytes = (args() < 3 | max_workspace_bytes >= . ? 1073741824 : max_workspace_bytes)
    b.max_work_ops = (args() < 4 ? . : max_work_ops)
    b.enforce_work_ops = (b.max_work_ops < .)
    cdid_imputation_validate_budget(b)
    return(b)
}

void cdid_imputation_validate_budget(struct cdid_imputation_budget scalar b)
{
    if (b.max_chunk_rows < 1 | b.max_chunk_rows >= . | b.max_chunk_rows != floor(b.max_chunk_rows) |
        b.max_cached_pairs < 1 | b.max_cached_pairs >= . | b.max_cached_pairs != floor(b.max_cached_pairs) |
        b.max_workspace_bytes <= 0 | b.max_workspace_bytes >= . |
        (b.enforce_work_ops != 0 & b.enforce_work_ops != 1))
        _cdid_imp_fail("Invalid explicit caller budget (require integer chunk/cache bounds and positive finite workspace bytes).", 498)
    if (b.enforce_work_ops & (b.max_work_ops <= 0 | b.max_work_ops >= .))
        _cdid_imp_fail("Explicit max_work_ops budget must be positive and finite when enabled.", 498)
}

void _cdid_imp_check_budget(
    struct cdid_imputation_budget scalar b,
    real scalar retained_bytes,
    real scalar scratch_bytes,
    real scalar work_ops,
    string scalar context)
{
    cdid_imputation_validate_budget(b)
    if (retained_bytes < 0 | scratch_bytes < 0 | work_ops < 0 |
        retained_bytes >= . | scratch_bytes >= . | work_ops >= .)
        _cdid_imp_fail("Nonfinite workspace or work-ops estimate in " + context + ".", 504)
    if (retained_bytes + scratch_bytes > b.max_workspace_bytes)
        _cdid_imp_fail("Caller workspace byte budget exceeded in " + context + ".", 498)
    if (b.enforce_work_ops & work_ops > b.max_work_ops)
        _cdid_imp_fail("Caller work-ops budget exceeded in " + context + ".", 498)
}

struct cdid_imputation_pair scalar cdid_imputation_make_pair(
    real scalar s_index,
    real scalar t_index,
    | real scalar recipe,
      string scalar pair_key)
{
    struct cdid_imputation_pair scalar pr
    pr.pair_id = 0
    pr.source_pair_id = 0
    pr.s_index = s_index
    pr.t_index = t_index
    pr.recipe = (args() < 3 | recipe >= . | recipe == 0 ? 1 : recipe)
    pr.L = t_index - s_index
    pr.pair_label = "e" + strofreal(s_index, "%21.0g") + strofreal(t_index, "%21.0g")
    pr.pair_key = (args() < 4 ? "" : pair_key)
    return(pr)
}

/* Initialize bounded LRU pair reader. Never materializes N x all_pairs or N x T x T. */
struct cdid_imputation_pair_reader scalar cdid_imputation_pair_reader_init(
    struct cdid_imputation_budget scalar budget)
{
    struct cdid_imputation_pair_reader scalar rdr
    real scalar i
    cdid_imputation_validate_budget(budget)
    rdr.max_cached_pairs = budget.max_cached_pairs
    rdr.chunk_rows = budget.max_chunk_rows
    rdr.clock = 0
    rdr.n_queries = 0
    rdr.n_hits = 0
    rdr.n_misses = 0
    rdr.n_evictions = 0
    rdr.current_cached_pairs = 0
    rdr.slots = _cdid_imp_cache_slot(1, rdr.max_cached_pairs)
    for (i=1; i<=rdr.max_cached_pairs; i++) {
        rdr.slots[i].pair_id = 0
        rdr.slots[i].s_index = 0
        rdr.slots[i].t_index = 0
        rdr.slots[i].recipe = 0
        rdr.slots[i].pair_key = ""
        rdr.slots[i].stamp = 0
        rdr.slots[i].valid = 0
        rdr.slots[i].obs_f = J(0, 1, 0)
        rdr.slots[i].dy = J(0, 1, 0)
    }
    rdr.date_posts = _cdid_imp_date_post(1, 0)
    rdr.date_ready = 0
    rdr.date_N = 0
    rdr.date_R = 0
    rdr.date_T = 0
    return(rdr)
}

/* Lazily index observed panel rows per calendar date in O(R_obs) memory. */
void _cdid_imp_ensure_date_posts(
    struct cdid_imputation_pair_reader scalar reader,
    struct cdid_population scalar pop)
{
    real scalar T, n_ord, d
    real colvector ord, d_vec, y_vec, u_vec, valid_y, ix
    T = rows(pop.calendar)
    n_ord = rows(pop.observed_order)
    if (reader.date_ready == 1 & reader.date_N == pop.N &
        reader.date_R == n_ord & reader.date_T == T & cols(reader.date_posts) == T) {
        return
    }
    reader.date_posts = _cdid_imp_date_post(1, T)
    if (n_ord > 0) {
        ord = pop.observed_order
        d_vec = pop.date[ord]
        y_vec = pop.y[ord]
        u_vec = pop.ids.row_unit[ord]
        valid_y = (y_vec :< .)
        for (d=1; d<=T; d++) {
            ix = _cdid_imp_index((d_vec :== pop.calendar[d]) :& valid_y)
            reader.date_posts[d].u = (rows(ix) ? u_vec[ix] : J(0, 1, 0))
            reader.date_posts[d].y = (rows(ix) ? y_vec[ix] : J(0, 1, 0))
        }
    }
    else {
        for (d=1; d<=T; d++) {
            reader.date_posts[d].u = J(0, 1, 0)
            reader.date_posts[d].y = J(0, 1, 0)
        }
    }
    reader.date_ready = 1
    reader.date_N = pop.N
    reader.date_R = n_ord
    reader.date_T = T
}

/* Populate one compact observed-pair posting slot from sparse population postings
   in bounded chunks over pop.observed_order. */
void _cdid_imp_populate_pair_slot(
    struct _cdid_imp_cache_slot scalar slot,
    struct cdid_population scalar pop,
    struct cdid_imputation_pair scalar pr,
    real scalar chunk_rows)
{
    real scalar n_ord, step, lo, hi, s_date, t_date, nobs, has_carry, carry_u, carry_ys, m, k
    real colvector ord_ch, d_ch, st_ix, r_st, u_st, d_st, y_st, pair_hit, obs_f, obs_dy
    s_date = pop.calendar[pr.s_index]
    t_date = pop.calendar[pr.t_index]
    obs_f = J(rows(pop.A_to_F), 1, 0)
    obs_dy = J(rows(pop.A_to_F), 1, 0)
    nobs = 0
    n_ord = rows(pop.observed_order)
    step = max((1, 2 * chunk_rows))
    has_carry = 0
    carry_u = 0
    carry_ys = 0
    for (lo=1; lo<=n_ord; lo=lo+step) {
        hi = min((n_ord, lo + step - 1))
        ord_ch = pop.observed_order[lo..hi]
        d_ch = pop.date[ord_ch]
        st_ix = _cdid_imp_index(((d_ch :== s_date) :| (d_ch :== t_date)) :& (pop.y[ord_ch] :< .))
        if (rows(st_ix) == 0) continue
        r_st = ord_ch[st_ix]
        u_st = pop.ids.row_unit[r_st]
        d_st = d_ch[st_ix]
        y_st = pop.y[r_st]
        if (has_carry) {
            u_st = carry_u \ u_st
            d_st = s_date \ d_st
            y_st = carry_ys \ y_st
        }
        m = rows(u_st)
        if (m >= 2) {
            pair_hit = _cdid_imp_index((u_st[1..(m - 1)] :== u_st[2..m]) :&
                (d_st[1..(m - 1)] :== s_date) :& (d_st[2..m] :== t_date))
            k = rows(pair_hit)
            if (k > 0) {
                obs_f[(nobs + 1)..(nobs + k)] = u_st[pair_hit]
                obs_dy[(nobs + 1)..(nobs + k)] = y_st[pair_hit :+ 1] :- y_st[pair_hit]
                nobs = nobs + k
            }
        }
        if (d_st[m] == s_date) {
            has_carry = 1
            carry_u = u_st[m]
            carry_ys = y_st[m]
        }
        else {
            has_carry = 0
        }
    }
    obs_f = (nobs ? obs_f[1..nobs] : J(0, 1, 0))
    obs_dy = (nobs ? obs_dy[1..nobs] : J(0, 1, 0))
    if (missing(obs_dy))
        _cdid_imp_fail("Nonfinite observed outcome difference for " + pr.pair_label + ".", 504)
    slot.pair_id = pr.pair_id
    slot.s_index = pr.s_index
    slot.t_index = pr.t_index
    slot.recipe = pr.recipe
    slot.pair_key = pr.pair_key
    slot.valid = 1
    slot.obs_f = obs_f
    slot.dy = obs_dy
}

void _cdid_imp_populate_slot_fast(
    struct cdid_imputation_pair_reader scalar reader,
    struct _cdid_imp_cache_slot scalar slot,
    struct cdid_population scalar pop,
    struct cdid_imputation_pair scalar pr)
{
    real colvector us, ys, ut, yt, pos_s, hit_t, obs_f, obs_dy
    if (pop.N > 262144) {
        _cdid_imp_populate_pair_slot(slot, pop, pr, reader.chunk_rows)
        return
    }
    _cdid_imp_ensure_date_posts(reader, pop)
    us = reader.date_posts[pr.s_index].u
    ys = reader.date_posts[pr.s_index].y
    ut = reader.date_posts[pr.t_index].u
    yt = reader.date_posts[pr.t_index].y
    if (rows(us) == 0 | rows(ut) == 0) {
        obs_f = J(0, 1, 0)
        obs_dy = J(0, 1, 0)
    }
    else {
        pos_s = J(pop.N, 1, 0)
        pos_s[us] = (1::rows(us))
        hit_t = _cdid_imp_index(pos_s[ut] :> 0)
        if (rows(hit_t) > 0) {
            obs_f = ut[hit_t]
            obs_dy = yt[hit_t] :- ys[pos_s[obs_f]]
        }
        else {
            obs_f = J(0, 1, 0)
            obs_dy = J(0, 1, 0)
        }
    }
    if (missing(obs_dy))
        _cdid_imp_fail("Nonfinite observed outcome difference for " + pr.pair_label + ".", 504)
    slot.pair_id = pr.pair_id
    slot.s_index = pr.s_index
    slot.t_index = pr.t_index
    slot.recipe = pr.recipe
    slot.pair_key = pr.pair_key
    slot.valid = 1
    slot.obs_f = obs_f
    slot.dy = obs_dy
}

/* Read observed pair indicator S and outcome change dy = Y_t - Y_s for a bounded
   chunk of original-F rows using the LRU cache. Returns B x 2 matrix (S, dy). */
real matrix cdid_imputation_read_pair_chunk(
    struct cdid_imputation_pair_reader scalar reader,
    struct cdid_population scalar pop,
    struct cdid_imputation_pair scalar pr,
    real colvector f_rows,
    real scalar read_dy)
{
    real scalar B, i, slot_idx, empty_idx, lru_idx, min_stamp, pos, n_obs, p_lo, p_hi
    real scalar f_min, f_max, span, f0
    real colvector sub_idx, pos_map, hit_q, ix_hit
    real matrix out
    B = rows(f_rows)
    if (B < 1 | B > reader.chunk_rows | cols(f_rows) != 1 | missing(f_rows) |
        any(f_rows :< 1) | any(f_rows :> pop.N) | any(f_rows :!= floor(f_rows)))
        _cdid_imp_fail("Invalid or oversized original-F row chunk in pair reader.", 498)
    reader.n_queries = reader.n_queries + 1
    reader.clock = reader.clock + 1
    slot_idx = 0
    empty_idx = 0
    lru_idx = 1
    min_stamp = .
    for (i=1; i<=reader.max_cached_pairs; i++) {
        if (reader.slots[i].valid) {
            if (reader.slots[i].s_index == pr.s_index &
                reader.slots[i].t_index == pr.t_index) {
                slot_idx = i
                break
            }
            if (reader.slots[i].stamp < min_stamp) {
                min_stamp = reader.slots[i].stamp
                lru_idx = i
            }
        }
        else if (empty_idx == 0) {
            empty_idx = i
        }
    }
    if (slot_idx > 0) {
        reader.n_hits = reader.n_hits + 1
        reader.slots[slot_idx].stamp = reader.clock
    }
    else {
        reader.n_misses = reader.n_misses + 1
        if (empty_idx > 0) {
            slot_idx = empty_idx
            reader.current_cached_pairs = reader.current_cached_pairs + 1
        }
        else {
            slot_idx = lru_idx
            reader.n_evictions = reader.n_evictions + 1
        }
        _cdid_imp_populate_slot_fast(reader, reader.slots[slot_idx], pop, pr)
        reader.slots[slot_idx].stamp = reader.clock
    }
    out = J(B, 2, 0)
    n_obs = rows(reader.slots[slot_idx].obs_f)
    if (n_obs == 0) return(out)
    if (B == 1) {
        pos = _cdid_imp_find_exact(reader.slots[slot_idx].obs_f, f_rows[1])
        if (pos > 0) {
            out[1, 1] = 1
            if (read_dy) out[1, 2] = reader.slots[slot_idx].dy[pos]
        }
        return(out)
    }
    if (all(f_rows[2..B] :> f_rows[1..(B - 1)])) {
        p_lo = _cdid_imp_lower_bound(reader.slots[slot_idx].obs_f, f_rows[1])
        if (p_lo > n_obs) return(out)
        p_hi = _cdid_imp_lower_bound(reader.slots[slot_idx].obs_f, f_rows[B] + 1) - 1
        if (p_lo > p_hi) return(out)
        span = f_rows[B] - f_rows[1] + 1
        if (span == B) {
            sub_idx = reader.slots[slot_idx].obs_f[p_lo..p_hi] :- (f_rows[1] - 1)
            out[sub_idx, 1] = J(rows(sub_idx), 1, 1)
            if (read_dy) out[sub_idx, 2] = reader.slots[slot_idx].dy[p_lo..p_hi]
            return(out)
        }
        if (span <= 262144) {
            f0 = f_rows[1] - 1
            pos_map = J(span, 1, 0)
            pos_map[reader.slots[slot_idx].obs_f[p_lo..p_hi] :- f0] = (p_lo::p_hi)
            hit_q = pos_map[f_rows :- f0]
            ix_hit = _cdid_imp_index(hit_q :> 0)
            if (rows(ix_hit) > 0) {
                out[ix_hit, 1] = J(rows(ix_hit), 1, 1)
                if (read_dy) out[ix_hit, 2] = reader.slots[slot_idx].dy[hit_q[ix_hit]]
            }
            return(out)
        }
        pos = p_lo
        for (i=1; i<=B; i++) {
            while (pos <= p_hi) {
                if (reader.slots[slot_idx].obs_f[pos] >= f_rows[i]) break
                pos = pos + 1
            }
            if (pos > p_hi) break
            if (reader.slots[slot_idx].obs_f[pos] == f_rows[i]) {
                out[i, 1] = 1
                if (read_dy) out[i, 2] = reader.slots[slot_idx].dy[pos]
                pos = pos + 1
            }
        }
        return(out)
    }
    f_min = min(f_rows)
    f_max = max(f_rows)
    p_lo = _cdid_imp_lower_bound(reader.slots[slot_idx].obs_f, f_min)
    if (p_lo <= n_obs) {
        p_hi = _cdid_imp_lower_bound(reader.slots[slot_idx].obs_f, f_max + 1) - 1
        span = f_max - f_min + 1
        if (p_lo <= p_hi & span <= 262144) {
            f0 = f_min - 1
            pos_map = J(span, 1, 0)
            pos_map[reader.slots[slot_idx].obs_f[p_lo..p_hi] :- f0] = (p_lo::p_hi)
            hit_q = pos_map[f_rows :- f0]
            ix_hit = _cdid_imp_index(hit_q :> 0)
            if (rows(ix_hit) > 0) {
                out[ix_hit, 1] = J(rows(ix_hit), 1, 1)
                if (read_dy) out[ix_hit, 2] = reader.slots[slot_idx].dy[hit_q[ix_hit]]
            }
            return(out)
        }
    }
    for (i=1; i<=B; i++) {
        pos = _cdid_imp_find_exact(reader.slots[slot_idx].obs_f, f_rows[i])
        if (pos > 0) {
            out[i, 1] = 1
            if (read_dy) out[i, 2] = reader.slots[slot_idx].dy[pos]
        }
    }
    return(out)
}

/* Build canonical deduplicated pairs and active (cohort_value, raw_pair_slot)
   requests from link specifications:
   link_spec is L x 3 or L x 4: (s_date, t_date, target_cohort_g [, recipe_id]).
   control_cohorts_list is M x 2: (link_index, control_cohort_k). */
void cdid_imp_links_to_requests(
    real colvector calendar,
    real matrix link_spec,
    real matrix control_cohorts_list,
    struct cdid_imputation_pair rowvector pairs_out,
    real matrix active_requests_out,
    real colvector link_pair_map_out)
{
    real scalar L_links, ell, s_idx, t_idx, rec, found, np, m, l_idx
    struct cdid_imputation_pair rowvector raw_p
    real matrix req
    L_links = rows(link_spec)
    if (L_links < 1 | (cols(link_spec) != 3 & cols(link_spec) != 4) | missing(link_spec) |
        rows(calendar) < 2 | missing(calendar) | cols(control_cohorts_list) != 2 |
        missing(control_cohorts_list))
        _cdid_imp_fail("Invalid link_spec or control_cohorts_list in cdid_imp_links_to_requests.", 498)
    raw_p = cdid_imputation_pair(1, L_links)
    np = 0
    link_pair_map_out = J(L_links, 1, 0)
    req = J(0, 2, 0)
    for (ell=1; ell<=L_links; ell++) {
        s_idx = _cdid_imp_find_exact(calendar, link_spec[ell, 1])
        t_idx = _cdid_imp_find_exact(calendar, link_spec[ell, 2])
        rec = (cols(link_spec) == 4 ? link_spec[ell, 4] : 1)
        if (s_idx < 1 | t_idx <= s_idx | rec < 1 | rec != floor(rec))
            _cdid_imp_fail("Link dates must belong to calendar with s < t and positive integer design block.", 498)
        found = 0
        for (m=1; m<=np; m++) {
            if (raw_p[m].s_index == s_idx & raw_p[m].t_index == t_idx & raw_p[m].recipe == rec) {
                found = m
                break
            }
        }
        if (found == 0) {
            np = np + 1
            found = np
            raw_p[np] = cdid_imputation_make_pair(s_idx, t_idx, rec)
        }
        link_pair_map_out[ell] = found
        req = req \ (link_spec[ell, 3], found)
    }
    for (m=1; m<=rows(control_cohorts_list); m++) {
        l_idx = control_cohorts_list[m, 1]
        if (l_idx < 1 | l_idx > L_links | l_idx != floor(l_idx))
            _cdid_imp_fail("Invalid link index in control_cohorts_list.", 498)
        req = req \ (control_cohorts_list[m, 2], link_pair_map_out[l_idx])
    }
    pairs_out = cdid_imputation_pair(1, np)
    for (m=1; m<=np; m++) pairs_out[m] = raw_p[m]
    active_requests_out = req
}

/* Verify alias reconstruction on actual observed-pair training rows AND
   all full-F cohort prediction rows for a group, matching or_pooled_alias_check
   and or_pooled_fit_blocks in R/chained_or_pooled.R. */
void _cdid_imp_verify_group_alias(
    struct cdid_imputation_model scalar model,
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
    /* First pass: training supports (observed-pair rows of each cell in group) */
    for (loc=1; loc<=B_r; loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        pid = model.cells[cid].pair_id
        r = model.cells[cid].recipe
        coh_rows = _cdid_imp_coh_rows(model, ci)
        n_coh = rows(coh_rows)
        for (lo=1; lo<=n_coh; lo=lo+model.budget.max_chunk_rows) {
            hi = min((n_coh, lo + model.budget.max_chunk_rows - 1))
            r_chunk = coh_rows[lo..hi]
            hit = cdid_imputation_read_pair_chunk(model.reader, model.pop, model.pairs[pid], r_chunk, 0)
            obs_ix = _cdid_imp_index(hit[, 1] :== 1)
            if (rows(obs_ix) == 0) continue
            Z_sub = model.recipes[r].Z[r_chunk[obs_ix], ]
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
    /* Second pass: required prediction supports (ALL full-F rows of each cell's cohort) */
    for (loc=1; loc<=B_r; loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        r = model.cells[cid].recipe
        coh_rows = _cdid_imp_coh_rows(model, ci)
        n_coh = rows(coh_rows)
        for (lo=1; lo<=n_coh; lo=lo+model.budget.max_chunk_rows) {
            hi = min((n_coh, lo + model.budget.max_chunk_rows - 1))
            r_chunk = coh_rows[lo..hi]
            Z_sub = model.recipes[r].Z[r_chunk, ]
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

/* Main Stage 1 fitter:
   - Canonicalizes cohorts, pairs (by L, s_index, t_index, recipe), active cells C_imp,
     and slope-pooling partition R_imp ("separate", "within_cohort", "across_cohorts" per recipe).
   - Streams observed pairs via bounded LRU pair reader and compresses per-cell QR factors.
   - Fits each group r in R_imp using cdid_pooled_wls_fit with full-F cohort alias verification. */
struct cdid_imputation_model scalar cdid_imputation_fit(
    struct cdid_population scalar pop,
    struct cdid_imputation_recipe rowvector recipes,
    struct cdid_imputation_pair rowvector raw_pairs,
    real matrix active_cell_requests,
    string scalar pooling,
    real scalar intercept,
    real scalar rank_tol,
    struct cdid_imputation_budget scalar budget,
    | string rowvector slope_names,
      real colvector elig)
{
    struct cdid_imputation_model scalar out
    struct cdid_wls_result scalar wfit
    real scalar N, G, nr, npair, p, r, k, pi, m, ci, pid, cid, gid, loc, j, rk
    real scalar n_coh, lo, hi, nn, scale, nPseudo, at, thr, orig_norm, res_norm
    real scalar exact_j, k_keep
    real colvector u_coh, sort_keys, coh_ord, pair_ord, raw_to_canon_pair, ix, ix_A
    real colvector coh_rows, r_chunk, obs_ix, ww, yy, cids, g_cohorts, rr
    real colvector pseudo_y, pseudo_w, pseudo_block, v, qy, coef, jhat_chunk, yc
    real rowvector constant_cols, retained, proposed, aliased, aliased_pos, candidates
    real rowvector keep, colmax_prop, ctl_err, trt_err, tau, pivot
    real matrix pair_key_mat, cell_mat, hit, rawZ, D, F_c, Q, R1, V, coef_rows, res_fac
    real matrix A_plus, A_minus, pseudo_X, C, mean_D_raw, map, constants_Z, constants_D
    real matrix Z_keep, map_full, rc_Fw, rc_FZ
    real scalar f_wsum
    real rowvector f_sumZ

    cdid_imputation_validate_budget(budget)
    if (pooling != "separate" & pooling != "within_cohort" & pooling != "across_cohorts")
        _cdid_imp_fail("imputation_pooling must be separate, within_cohort, or across_cohorts.", 498)
    if ((intercept != 0 & intercept != 1) | rank_tol <= 0 | rank_tol >= 1 | rank_tol >= .)
        _cdid_imp_fail("Require binary intercept flag and rank_tol in (0, 1).", 498)
    N = pop.N
    if (N < 1 | rows(pop.cohort_known) != N | any(pop.cohort_known :!= 1))
        _cdid_imp_fail("Imputation requires full representative roster F with known cohort on every unit.", 459)

    /* Canonical cohort ordering: never-treated (0 / .) first, then finite cohorts ascending */
    u_coh = editmissing(pop.cohort, 0)
    if (any(u_coh :< 0)) _cdid_imp_fail("Cohort values must be nonnegative.", 459)
    out.cohorts = uniqrows(sort(u_coh, 1))
    G = rows(out.cohorts)
    sort_keys = J(G, 1, 0)
    for (k=1; k<=G; k++) sort_keys[k] = (out.cohorts[k] == 0 ? -1 : out.cohorts[k])
    coh_ord = order((sort_keys, (1::G)), (1, 2))
    out.cohorts = out.cohorts[coh_ord]
    out.cohort_labels = J(G, 1, "")
    for (k=1; k<=G; k++) out.cohort_labels[k] = _cdid_imp_cohort_label(out.cohorts[k])
    out.unit_cohort_index = J(N, 1, 0)
    out.coh_u = _cdid_imp_ivec(1, G)
    for (k=1; k<=G; k++) {
        ix = _cdid_imp_index(u_coh :== out.cohorts[k])
        out.coh_u[k].rows = ix
        if (rows(ix) > 0) out.unit_cohort_index[ix] = J(rows(ix), 1, k)
    }
    if (any(out.unit_cohort_index :== 0))
        _cdid_imp_fail("Failed to map unit cohort to canonical cohort dictionary.", 498)
    /* A restricted reference eligibility context: every cohort reads its eligible units. */
    out.elig_keep = J(0, 0, .)
    if (args() >= 10) if (rows(elig) == N & cols(elig) == 1) {
        out.elig_keep = (elig :!= 0)
        for (k=1; k<=G; k++) {
            ix = out.coh_u[k].rows
            if (rows(ix)) out.coh_u[k].rows = select(ix, out.elig_keep[ix])
        }
    }

    /* Validate recipes and compute per-recipe cohort masses and repeat rates on F */
    nr = cols(recipes)
    if (nr < 1) _cdid_imp_fail("Require at least one imputation design block.", 498)
    p = cols(recipes[1].Z)
    if (!intercept & p == 0) _cdid_imp_fail("No-intercept outcome design is empty.", 459)
    out.cohort_mass_F = J(nr, G, 0)
    out.cohort_repeat_mass_F = J(nr, G, 0)
    out.cohort_a_F = J(nr, G, 0)
    out.cohort_tau_F = J(nr, G, 0)
    out.cohort_eta_F = J(nr, G, 0)
    out.mean_w_F = J(nr, 1, 0)
    out.mean_wA_F = J(nr, 1, 0)
    for (r=1; r<=nr; r++) {
        if (rows(recipes[r].weight) != N | cols(recipes[r].weight) != 1 |
            missing(recipes[r].weight) | any(recipes[r].weight :<= 0) |
            rows(recipes[r].Z) != N | cols(recipes[r].Z) != p | missing(recipes[r].Z))
            _cdid_imp_fail("Each design block requires positive finite N x 1 weights and finite N x p full-F design Z.", 498)
        if (recipes[r].recipe_id == 0 | recipes[r].recipe_id >= .) recipes[r].recipe_id = r
        if (recipes[r].recipe_key == "") recipes[r].recipe_key = "recipe:" + strofreal(recipes[r].recipe_id, "%21.0g")
        for (k=1; k<=G; k++) {
            ix = out.coh_u[k].rows
            ix_A = (rows(ix) ? select(ix, pop.A[ix] :== 1) : J(0, 1, 0))
            out.cohort_mass_F[r, k] = sum(recipes[r].weight[ix]) / N
            out.cohort_repeat_mass_F[r, k] = (rows(ix_A) ? sum(recipes[r].weight[ix_A]) / N : 0)
            out.cohort_a_F[r, k] = (out.cohort_mass_F[r, k] > 0 ? out.cohort_repeat_mass_F[r, k] / out.cohort_mass_F[r, k] : 0)
        }
        out.mean_w_F[r] = sum(out.cohort_mass_F[r, ])
        out.mean_wA_F[r] = sum(out.cohort_repeat_mass_F[r, ])
        if (out.mean_w_F[r] <= 0 | out.mean_wA_F[r] <= 0)
            _cdid_imp_fail("Total F mass and repeat-observer A mass must be positive.", 459)
        out.cohort_tau_F[r, ] = out.cohort_mass_F[r, ] :/ out.mean_w_F[r]
        out.cohort_eta_F[r, ] = out.cohort_repeat_mass_F[r, ] :/ out.mean_wA_F[r]
    }

    /* Validate and canonically order pairs by (L, s_index, t_index, recipe, slot) */
    npair = cols(raw_pairs)
    if (npair < 1) _cdid_imp_fail("Require at least one outcome pair.", 498)
    pair_key_mat = J(npair, 5, 0)
    for (pi=1; pi<=npair; pi++) {
        if (raw_pairs[pi].recipe == 0 | raw_pairs[pi].recipe >= .) raw_pairs[pi].recipe = 1
        r = raw_pairs[pi].recipe
        if (r < 1 | r > nr | r != floor(r) |
            raw_pairs[pi].s_index < 1 | raw_pairs[pi].t_index <= raw_pairs[pi].s_index |
            raw_pairs[pi].t_index > rows(pop.calendar) |
            raw_pairs[pi].s_index != floor(raw_pairs[pi].s_index) |
            raw_pairs[pi].t_index != floor(raw_pairs[pi].t_index))
            _cdid_imp_fail("Invalid pair calendar indices or design block index.", 498)
        raw_pairs[pi].L = raw_pairs[pi].t_index - raw_pairs[pi].s_index
        if (raw_pairs[pi].pair_label == "")
            raw_pairs[pi].pair_label = "e" + strofreal(raw_pairs[pi].s_index, "%21.0g") + strofreal(raw_pairs[pi].t_index, "%21.0g")
        pair_key_mat[pi, ] = (raw_pairs[pi].L, raw_pairs[pi].s_index, raw_pairs[pi].t_index, r, pi)
    }
    pair_ord = order(pair_key_mat, (1, 2, 3, 4, 5))
    out.pairs = cdid_imputation_pair(1, npair)
    raw_to_canon_pair = J(npair, 1, 0)
    for (pi=1; pi<=npair; pi++) {
        out.pairs[pi] = raw_pairs[pair_ord[pi]]
        out.pairs[pi].source_pair_id = raw_pairs[pair_ord[pi]].pair_id
        if (out.pairs[pi].pair_key == "") out.pairs[pi].pair_key = "pair:" + strofreal(pair_ord[pi], "%21.0g")
        out.pairs[pi].pair_id = pi
        raw_to_canon_pair[pair_ord[pi]] = pi
    }

    /* Build active cells C_imp from active_cell_requests: M x 2 (cohort_value, raw_pair_slot) */
    if (rows(active_cell_requests) < 1 | cols(active_cell_requests) != 2 | missing(active_cell_requests))
        _cdid_imp_fail("Require nonempty M x 2 active_cell_requests (cohort_value, raw_pair_slot).", 498)
    cell_mat = J(G, npair, 0)
    for (m=1; m<=rows(active_cell_requests); m++) {
        ci = _cdid_imp_find_exact(out.cohorts, active_cell_requests[m, 1])
        pi = active_cell_requests[m, 2]
        if (ci == 0 | pi < 1 | pi > npair | pi != floor(pi))
            _cdid_imp_fail("Active cell request references unknown cohort or pair slot.", 498)
        pid = raw_to_canon_pair[pi]
        cell_mat[ci, pid] = 1
    }
    out.n_cells = sum(cell_mat)
    if (out.n_cells < 1) _cdid_imp_fail("No active imputation cells declared.", 498)
    out.cells = cdid_imputation_cell(1, out.n_cells)
    out.cell_lookup = J(G, npair, 0)
    out.active_cohorts = J(0, 1, 0)
    cid = 0
    for (ci=1; ci<=G; ci++) {
        if (any(cell_mat[ci, ] :== 1)) out.active_cohorts = out.active_cohorts \ ci
        for (pid=1; pid<=npair; pid++) {
            if (!cell_mat[ci, pid]) continue
            cid = cid + 1
            r = out.pairs[pid].recipe
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
                (nr > 1 ? ".r" + strofreal(recipes[r].recipe_id, "%21.0g") : "")
            out.cell_lookup[ci, pid] = cid
        }
    }

    /* Build slope-pooling partition R_imp inside each recipe */
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
        /* One group per active (cohort_index, recipe) in canonical order */
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
        /* across_cohorts: one group per active recipe (global group0 when nr == 1) */
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

    /* Preflight workspace and work-ops bounds against caller budget */
    out.retained_bytes_bound = 8 * (nr * N * (p + 2) + out.n_cells * ((p + 2)^2 + 8 * p + 32) +
        out.n_groups * (4 * p^2 + 8 * p + 32) + budget.max_cached_pairs * 2 * max((1, pop.n_A)) + 4096)
    out.scratch_bytes_bound = 8 * (min((N, budget.max_chunk_rows)) * (2 * p + 16) +
        2 * (p + 2) * max((1, out.n_cells)) * (p + 4) + 4096)
    out.work_ops_bound = 128 * N * (out.n_cells + 1) * (p + 2)^2 * (1 + ceil(ln(N + 1) / ln(2)))
    _cdid_imp_check_budget(budget, out.retained_bytes_bound, out.scratch_bytes_bound,
        out.work_ops_bound, "Stage 1 imputation fitting preflight")

    out.pop = pop
    out.budget = budget
    out.reader = cdid_imputation_pair_reader_init(budget)
    out.recipes = recipes
    out.N = N
    out.n_A = pop.n_A
    out.lambda = pop.n_A / N
    out.n_cohorts = G
    out.n_recipes = nr
    out.n_pairs = npair
    out.intercept = intercept
    out.p_raw = p
    out.rank_tol = rank_tol
    out.pooling = pooling
    if (args() >= 9 & cols(slope_names) == p) out.slope_names = slope_names
    else {
        out.slope_names = J(1, p, "")
        for (j=1; j<=p; j++) out.slope_names[j] = "x" + strofreal(j, "%21.0g")
    }

    /* Precompute full-F cohort colmax per (recipe, cohort_index) in bounded chunks */
    out.cohort_colmax_F = J(nr * G, p, 0)
    if (p > 0) {
        for (r=1; r<=nr; r++) {
            for (m=1; m<=rows(out.active_cohorts); m++) {
                ci = out.active_cohorts[m]
                coh_rows = out.coh_u[ci].rows
                n_coh = rows(coh_rows)
                for (lo=1; lo<=n_coh; lo=lo+budget.max_chunk_rows) {
                    hi = min((n_coh, lo + budget.max_chunk_rows - 1))
                    r_chunk = coh_rows[lo..hi]
                    out.cohort_colmax_F[(r - 1) * G + ci, ] = colmax(
                        out.cohort_colmax_F[(r - 1) * G + ci, ] \ abs(recipes[r].Z[r_chunk, ]))
                }
            }
        }
    }

    /* Stage 1a: Stream observed pairs per cell (visited pair-major for cache locality)
       and accumulate O((p+2)^2) QR summaries */
    rc_Fw = J(nr, G, .)
    rc_FZ = J(nr * G, p, 0)
    for (pid=1; pid<=npair; pid++) {
        for (ci=1; ci<=G; ci++) {
            cid = out.cell_lookup[ci, pid]
            if (cid == 0) continue
            r = out.cells[cid].recipe
            coh_rows = out.coh_u[ci].rows
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
                    ww = recipes[r].weight[r_chunk]
                    f_wsum = f_wsum + sum(ww)
                    if (p > 0) f_sumZ = f_sumZ + quadcross(ww, recipes[r].Z[r_chunk, ])
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
                ww = recipes[r].weight[r_chunk[obs_ix]]
                yy = hit[obs_ix, 2]
                rawZ = recipes[r].Z[r_chunk[obs_ix], ]
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

    /* Stage 1b: Fit each pooling group r in R_imp with full-F cohort alias verification */
    thr = 64 * epsilon(1) * max((1, p + intercept))
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
            if (intercept & p > 0) {
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

        /* Original-order near-dependency proposal on the partialled weighted group factor */
        C = pseudo_X
        mean_D_raw = J(rows(cids), p, 0)
        for (loc=1; loc<=rows(cids); loc++) {
            ix = _cdid_imp_index(pseudo_block :== loc)
            if (intercept & p > 0) {
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
            if (intercept) {
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
            _cdid_imp_verify_group_alias(out, gid, retained, proposed, map, constants_Z,
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
            if (!intercept & cols(keep) == 0)
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
            if (intercept & cols(aliased) > 0) {
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
            /* Reuse cdid_pooled_wls_fit on the kept columns */
            wfit = cdid_pooled_wls_fit(
                (k_keep > 0 ? pseudo_X[, keep] : J(nPseudo, 0, 0)),
                pseudo_y, pseudo_w, pseudo_block, J(nPseudo, 1, 1), N,
                intercept, J(0, k_keep, 0), J(0, 1, 0), rank_tol)
            if (cols(wfit.kept) != k_keep)
                _cdid_imp_fail("Pooled WLS final rank refused; no ridge or parameter rescue.", 459)
            out.groups[gid].beta = wfit.beta
            out.groups[gid].alpha = wfit.alpha
            out.groups[gid].R = wfit.R
            out.groups[gid].pivot = wfit.pivot
            out.groups[gid].rcond_R = wfit.rcond_R
            out.groups[gid].block_mass = wfit.mass
            out.groups[gid].mean_y = wfit.mean_y
            out.groups[gid].mean_D = (intercept & k_keep > 0 ? wfit.mean_X : J(rows(cids), k_keep, 0))
            out.groups[gid].rss = wfit.rss
        }
        else {
            /* Non-exact candidates retained by actual-row alias check: final QR rank check on keep */
            yc = pseudo_y
            out.groups[gid].block_mass = J(rows(cids), 1, 0)
            out.groups[gid].mean_y = J(rows(cids), 1, 0)
            for (loc=1; loc<=rows(cids); loc++) {
                ix = _cdid_imp_index(pseudo_block :== loc)
                out.groups[gid].block_mass[loc] = sum(pseudo_w[ix])
                out.groups[gid].mean_y[loc] = sum(pseudo_w[ix] :* pseudo_y[ix]) / out.groups[gid].block_mass[loc]
                if (intercept) yc[ix] = pseudo_y[ix] :- out.groups[gid].mean_y[loc]
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
            out.groups[gid].mean_D = (intercept ? mean_D_raw[, keep] : J(rows(cids), k_keep, 0))
            out.groups[gid].alpha = (intercept ? out.groups[gid].mean_y - out.groups[gid].mean_D * out.groups[gid].beta : J(rows(cids), 1, 0))
            v = pseudo_y - out.groups[gid].alpha[pseudo_block, 1] - pseudo_X[, keep] * out.groups[gid].beta
            out.groups[gid].rss = sum(pseudo_w :* v:^2)
        }

        out.groups[gid].p_raw = p
        out.groups[gid].p = k_keep
        out.groups[gid].intercept = intercept
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

        /* Populate cell fitted parameters and stream full-F cell prediction moments */
        for (loc=1; loc<=rows(cids); loc++) {
            cid = cids[loc]
            ci = out.cells[cid].cohort_index
            r = out.cells[cid].recipe
            out.cells[cid].alpha = (intercept ? out.groups[gid].alpha[loc] : 0)
            out.cells[cid].mean_d = (intercept & k_keep > 0 ? out.groups[gid].mean_D[loc, ] : J(1, k_keep, 0))
            out.cells[cid].d_bar_j = (k_keep > 0 ? out.cells[cid].L :* out.cells[cid].F_mean_Z[keep] : J(1, 0, 0))

            coh_rows = out.coh_u[ci].rows
            n_coh = rows(coh_rows)
            out.cells[cid].F_mean_j = 0
            out.cells[cid].F_sum_jZ = J(1, p, 0)
            for (lo=1; lo<=n_coh; lo=lo+budget.max_chunk_rows) {
                hi = min((n_coh, lo + budget.max_chunk_rows - 1))
                r_chunk = coh_rows[lo..hi]
                ww = recipes[r].weight[r_chunk]
                jhat_chunk = J(rows(r_chunk), 1, out.cells[cid].alpha)
                if (k_keep > 0) {
                    jhat_chunk = jhat_chunk + (out.cells[cid].L :* recipes[r].Z[r_chunk, keep]) * out.groups[gid].beta
                }
                out.cells[cid].F_mean_j = out.cells[cid].F_mean_j + sum(ww :* jhat_chunk)
                if (p > 0) {
                    out.cells[cid].F_sum_jZ = out.cells[cid].F_sum_jZ + quadcross(ww :* jhat_chunk, recipes[r].Z[r_chunk, ])
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

/* Static single-recipe convenience wrapper for cdid_imputation_fit */
struct cdid_imputation_model scalar cdid_imputation_fit_static(
    struct cdid_population scalar pop,
    real colvector weight_F,
    real matrix Z_F,
    struct cdid_imputation_pair rowvector raw_pairs,
    real matrix active_cell_requests,
    string scalar pooling,
    real scalar intercept,
    real scalar rank_tol,
    struct cdid_imputation_budget scalar budget,
    | string rowvector slope_names,
      real colvector keep)
{
    struct cdid_imputation_recipe rowvector rec
    rec = cdid_imputation_recipe(1, 1)
    rec[1].recipe_id = 1
    rec[1].weight = weight_F
    rec[1].Z = Z_F
    if (args() >= 11) {
        return(cdid_imputation_fit(pop, rec, raw_pairs, active_cell_requests,
            pooling, intercept, rank_tol, budget, slope_names, keep))
    }
    if (args() >= 10) {
        return(cdid_imputation_fit(pop, rec, raw_pairs, active_cell_requests,
            pooling, intercept, rank_tol, budget, slope_names))
    }
    return(cdid_imputation_fit(pop, rec, raw_pairs, active_cell_requests,
        pooling, intercept, rank_tol, budget))
}

void _cdid_imp_validate_chunk_rows(struct cdid_imputation_model scalar model, real colvector f_rows)
{
    real scalar B
    B = rows(f_rows)
    if (B < 1 | B > model.budget.max_chunk_rows | cols(f_rows) != 1 | missing(f_rows) |
        any(f_rows :< 1) | any(f_rows :> model.N) | any(f_rows :!= floor(f_rows)))
        _cdid_imp_fail("Invalid or oversized original-F row chunk.", 498)
    if (B > 1) {
        if (!all(f_rows[2..B] :> f_rows[1..(B - 1)])) {
            if (rows(uniqrows(f_rows)) != B)
                _cdid_imp_fail("Original-F row chunk must not contain duplicate row indices.", 498)
        }
    }
}

/* Predict j_hat_{c}(X_i) = alpha_c + L_c * Z_{r, i, kept_r} * beta_r on a bounded chunk of F rows */
real colvector cdid_imp_predict_cell_chunk(
    struct cdid_imputation_model scalar model,
    real scalar cid,
    real colvector f_rows)
{
    real scalar gid, r, B
    real colvector out
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imp_predict_cell_chunk.", 498)
    B = rows(f_rows)
    gid = model.cells[cid].group_id
    r = model.cells[cid].recipe
    out = J(B, 1, model.cells[cid].alpha)
    if (model.groups[gid].p > 0) {
        out = out + (model.cells[cid].L :* model.recipes[r].Z[f_rows, model.groups[gid].kept]) * model.groups[gid].beta
    }
    return(out)
}

/* Predict j_hat_{c}(X_i) across all active cells on a bounded chunk of F rows */
real matrix cdid_imp_predict_cells_chunk(
    struct cdid_imputation_model scalar model,
    real colvector f_rows)
{
    real scalar cid, B
    real matrix out
    _cdid_imp_validate_chunk_rows(model, f_rows)
    B = rows(f_rows)
    out = J(B, model.n_cells, 0)
    for (cid=1; cid<=model.n_cells; cid++) {
        out[, cid] = cdid_imp_predict_cell_chunk(model, cid, f_rows)
    }
    return(out)
}

/* Evaluate observed-pair cell residual score:
   s_{i,c} = w_{r,i} * I(g_i == k_c) * S_{i, e_c} * (dy_{i, e_c} - j_hat_c(X_i))
   Strictly zero off cohort k_c or when S_{i, e_c} == 0 (including off-A and A=1, S_e=0 units). */
real colvector cdid_imputation_cell_score_chunk(
    struct cdid_imputation_model scalar model,
    real scalar cid,
    real colvector f_rows)
{
    real scalar B, ci, pid, r, gid
    real colvector out, ix_coh, obs_local, rr, jhat_rr
    real matrix hit
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imputation_cell_score_chunk.", 498)
    B = rows(f_rows)
    out = J(B, 1, 0)
    ci = model.cells[cid].cohort_index
    ix_coh = _cdid_imp_index(model.unit_cohort_index[f_rows] :== ci)
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
        jhat_rr = jhat_rr + (model.cells[cid].L :* model.recipes[r].Z[f_rows[rr], model.groups[gid].kept]) * model.groups[gid].beta
    }
    out[rr] = model.recipes[r].weight[f_rows[rr]] :* (hit[obs_local, 2] :- jhat_rr)
    return(out)
}

/* Evaluate group slope score h_{i, beta_r} = sum_{c in C_r} s_{i,c} * tilde_d_{i,c} (B x p_r).
   Sums a recurring unit's contributions across all member cells c in C_r BEFORE any solve/product. */
real matrix cdid_imp_group_score_chunk(
    struct cdid_imputation_model scalar model,
    real scalar gid,
    real colvector f_rows)
{
    struct _cdid_imp_ivec rowvector c_ix
    real scalar B, k_keep, loc, cid, ci, pid, r
    real colvector cids, ix_coh, obs_local, rr, jhat_rr, s_rr
    real matrix h_beta, hit, d_rr, d_tilde, act_keys
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (gid < 1 | gid > model.n_groups | gid != floor(gid))
        _cdid_imp_fail("Invalid group index in cdid_imp_group_score_chunk.", 498)
    B = rows(f_rows)
    k_keep = model.groups[gid].p
    h_beta = J(B, k_keep, 0)
    if (k_keep == 0) return(h_beta)
    c_ix = _cdid_imp_bucket_coh(model.unit_cohort_index[f_rows], model.n_cohorts)
    cids = model.groups[gid].cell_ids
    if (rows(cids) > 1 & model.groups[gid].pooling == "across_cohorts") {
        act_keys = J(rows(cids), 4, 0)
        for (loc=1; loc<=rows(cids); loc++) {
            cid = cids[loc]
            act_keys[loc, ] = (model.cells[cid].L, model.cells[cid].s_index, model.cells[cid].t_index, cid)
        }
        cids = cids[order(act_keys, (1, 2, 3, 4))]
    }
    for (loc=1; loc<=rows(cids); loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        ix_coh = c_ix[ci].rows
        if (rows(ix_coh) == 0) continue
        pid = model.cells[cid].pair_id
        r = model.cells[cid].recipe
        hit = cdid_imputation_read_pair_chunk(model.reader, model.pop, model.pairs[pid], f_rows[ix_coh], 1)
        obs_local = _cdid_imp_index(hit[, 1] :== 1)
        if (rows(obs_local) == 0) continue
        rr = ix_coh[obs_local]
        d_rr = model.cells[cid].L :* model.recipes[r].Z[f_rows[rr], model.groups[gid].kept]
        jhat_rr = model.cells[cid].alpha :+ d_rr * model.groups[gid].beta
        s_rr = model.recipes[r].weight[f_rows[rr]] :* (hit[obs_local, 2] :- jhat_rr)
        d_tilde = (model.intercept ? d_rr :- model.cells[cid].mean_d : d_rr)
        h_beta[rr, ] = h_beta[rr, ] + s_rr :* d_tilde
    }
    return(h_beta)
}

/* Forward chunk evaluation of group slope IF:
   xi_{i, beta_r}^imp = N * (R_r' R_r)^{-1} h_{i, beta_r} on "frame" scale
   (or multiplied by lambda = n_A / N on "core" scale). */
real matrix cdid_imputation_beta_if_chunk(
    struct cdid_imputation_model scalar model,
    real scalar gid,
    real colvector f_rows,
    | string scalar scale)
{
    real matrix h_beta, xi_beta
    string scalar sc
    sc = (args() < 4 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    h_beta = cdid_imp_group_score_chunk(model, gid, f_rows)
    if (model.groups[gid].p == 0) return(h_beta)
    xi_beta = _cdid_nuisance_gram(model.groups[gid].R, model.groups[gid].pivot, (model.N :* h_beta)')'
    if (sc == "core") xi_beta = model.lambda :* xi_beta
    return(xi_beta)
}

/* Forward chunk evaluation of cell intercept IF:
   xi_{i, alpha_c}^imp = (s_{i,c} / mass_c) - xi_{i, beta_r}^imp * mean_d_c' */
real colvector cdid_imputation_alpha_if_chunk(
    struct cdid_imputation_model scalar model,
    real scalar cid,
    real colvector f_rows,
    | string scalar scale)
{
    real scalar gid
    real colvector s_c, xi_alpha
    real matrix xi_beta
    string scalar sc
    sc = (args() < 4 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (cid < 1 | cid > model.n_cells | cid != floor(cid))
        _cdid_imp_fail("Invalid cell index in cdid_imputation_alpha_if_chunk.", 498)
    if (!model.intercept) return(J(rows(f_rows), 1, 0))
    gid = model.cells[cid].group_id
    s_c = cdid_imputation_cell_score_chunk(model, cid, f_rows)
    xi_alpha = s_c :/ model.cells[cid].mass
    if (model.groups[gid].p > 0) {
        xi_beta = cdid_imputation_beta_if_chunk(model, gid, f_rows, "frame")
        xi_alpha = xi_alpha - xi_beta * model.cells[cid].mean_d'
    }
    if (sc == "core") xi_alpha = model.lambda :* xi_alpha
    return(xi_alpha)
}

/* Precompute sparse adjoint loadings for M linear functionals of Stage 1 parameters:
   lambda_mat is n_cells x M (loadings on alpha_c).
   nu_raw_mat is n_cells x (p_raw * M) (loadings on beta_{r(c)} in raw slope column order;
   for each contrast m=1..M, columns ((m-1)*p_raw + 1)..(m*p_raw) hold the p_raw slope loadings). */
struct cdid_imputation_adjoint scalar cdid_imputation_adjoint_prepare(
    struct cdid_imputation_model scalar model,
    real matrix lambda_mat,
    real matrix nu_raw_mat)
{
    struct cdid_imputation_adjoint scalar adj
    real scalar M, p, gid, k_keep, loc, cid, has_lam, has_nu
    real colvector cids, active_mask, g_active, ix_act
    real rowvector keep, lam_c, nu_row
    real matrix rhs, nu_c_kept, act_keys
    M = cols(lambda_mat)
    p = model.p_raw
    if (rows(lambda_mat) != model.n_cells | M < 1 |
        rows(nu_raw_mat) != model.n_cells | cols(nu_raw_mat) != p * M |
        missing(lambda_mat) | missing(nu_raw_mat))
        _cdid_imp_fail("Invalid lambda_mat or nu_raw_mat dimensions in cdid_imputation_adjoint_prepare.", 498)
    adj.n_cells = model.n_cells
    adj.n_groups = model.n_groups
    adj.n_contrasts = M
    adj.intercept = model.intercept
    adj.tilde_lambda = J(model.n_cells, M, 0)
    adj.groups = _cdid_imp_group_nu(1, model.n_groups)
    active_mask = J(model.n_cells, 1, 0)
    g_active = J(model.n_groups, 1, 0)
    for (gid=1; gid<=model.n_groups; gid++) {
        cids = model.groups[gid].cell_ids
        k_keep = model.groups[gid].p
        keep = model.groups[gid].kept
        rhs = J(k_keep, M, 0)
        for (loc=1; loc<=rows(cids); loc++) {
            cid = cids[loc]
            has_lam = 0
            if (model.intercept) {
                lam_c = lambda_mat[cid, ]
                has_lam = any(lam_c :!= 0)
                if (has_lam) {
                    adj.tilde_lambda[cid, ] = lam_c :/ model.cells[cid].mass
                    active_mask[cid] = 1
                }
            }
            if (k_keep > 0) {
                nu_row = nu_raw_mat[cid, ]
                has_nu = any(nu_row :!= 0)
                if (has_nu) {
                    nu_c_kept = (M == 1 ? nu_row[keep]' : colshape(nu_row, p)[, keep]')
                    rhs = rhs + nu_c_kept
                }
                if (has_lam) {
                    rhs = rhs - model.cells[cid].mean_d' * lam_c
                }
            }
        }
        if (k_keep > 0 & any(rhs :!= 0)) {
            adj.groups[gid].nu = _cdid_nuisance_gram(model.groups[gid].R, model.groups[gid].pivot, model.N :* rhs)
            g_active[gid] = 1
            for (loc=1; loc<=rows(cids); loc++) active_mask[cids[loc]] = 1
        }
        else {
            adj.groups[gid].nu = J(k_keep, M, 0)
        }
    }
    ix_act = _cdid_imp_index(active_mask :== 1)
    if (rows(ix_act) > 1) {
        act_keys = J(rows(ix_act), 4, 0)
        for (loc=1; loc<=rows(ix_act); loc++) {
            cid = ix_act[loc]
            act_keys[loc, ] = (model.cells[cid].L, model.cells[cid].s_index, model.cells[cid].t_index, cid)
        }
        ix_act = ix_act[order(act_keys, (1, 2, 3, 4))]
    }
    adj.active_cells = ix_act
    return(adj)
}

struct cdid_imputation_adjoint scalar cdid_imputation_adjoint_slice(
    struct cdid_imputation_model scalar model,
    struct cdid_imputation_adjoint scalar full_adj,
    real scalar c_lo,
    real scalar c_hi)
{
    struct cdid_imputation_adjoint scalar sub
    real scalar M_sub, gid, n_a, loc, cid
    real rowvector cols
    real colvector g_nz, keep_a
    cols = c_lo..c_hi
    M_sub = cols(cols)
    sub.n_cells = full_adj.n_cells
    sub.n_groups = full_adj.n_groups
    sub.n_contrasts = M_sub
    sub.intercept = full_adj.intercept
    sub.tilde_lambda = full_adj.tilde_lambda[, cols]
    sub.groups = _cdid_imp_group_nu(1, full_adj.n_groups)
    g_nz = J(full_adj.n_groups, 1, 0)
    for (gid=1; gid<=full_adj.n_groups; gid++) {
        if (model.groups[gid].p > 0) {
            sub.groups[gid].nu = full_adj.groups[gid].nu[, cols]
            if (any(sub.groups[gid].nu :!= 0)) g_nz[gid] = 1
        }
        else {
            sub.groups[gid].nu = J(0, M_sub, 0)
        }
    }
    n_a = rows(full_adj.active_cells)
    if (n_a == 0) {
        sub.active_cells = J(0, 1, 0)
        return(sub)
    }
    keep_a = J(n_a, 1, 0)
    for (loc=1; loc<=n_a; loc++) {
        cid = full_adj.active_cells[loc]
        if (g_nz[model.cells[cid].group_id]) keep_a[loc] = 1
        else if (full_adj.intercept) {
            if (any(sub.tilde_lambda[cid, ] :!= 0)) keep_a[loc] = 1
        }
    }
    sub.active_cells = select(full_adj.active_cells, keep_a)
    return(sub)
}

/* Evaluate sparse adjoint pullback on a bounded chunk of F rows:
   out[i, m] = sum_{c in C_imp} s_{i,c} * (tilde_lambda_{c,m} + tilde_d_{i,c} * tilde_nu_{r(c), m})
   Never materializes row-by-cell or row-by-group parameter IF matrices. */
real matrix cdid_imputation_adjoint_chunk(
    struct cdid_imputation_model scalar model,
    struct cdid_imputation_adjoint scalar adj,
    real colvector f_rows,
    | string scalar scale)
{
    struct _cdid_imp_ivec rowvector c_ix
    real scalar B, M, n_act, a_lo, a_hi, a_idx, cid0, cid_next, s0, t0, any_coh, cid, gid, ci, r
    real colvector ix_coh, obs_local, rr, jhat_rr, s_rr, g_done
    real rowvector mcols
    real matrix out, hit, d_rr, d_tilde, val, g_mask
    string scalar sc
    sc = (args() < 4 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    _cdid_imp_validate_chunk_rows(model, f_rows)
    B = rows(f_rows)
    M = adj.n_contrasts
    out = J(B, M, 0)
    n_act = rows(adj.active_cells)
    if (n_act == 0) return(out)
    g_done = J(model.n_groups, 1, 0)
    g_mask = J(model.n_groups, M, 0)
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
                    ix_coh = c_ix[ci].rows
                    if (rows(ix_coh) == 0) continue
                    obs_local = _cdid_imp_index(hit[ix_coh, 1] :== 1)
                    if (rows(obs_local) == 0) continue
                    rr = ix_coh[obs_local]
                    gid = model.cells[cid].group_id
                    r = model.cells[cid].recipe
                    /* Only contrasts with a nonzero loading on this cell can change out;
                       the remaining columns would receive exact zeros. */
                    if (!g_done[gid]) {
                        g_mask[gid, ] = (model.groups[gid].p > 0 ? (colsum(adj.groups[gid].nu :!= 0) :> 0) : J(1, M, 0))
                        g_done[gid] = 1
                    }
                    mcols = selectindex(g_mask[gid, ] :| (adj.tilde_lambda[cid, ] :!= 0))
                    if (cols(mcols) == 0) continue
                    if (model.groups[gid].p > 0) {
                        d_rr = model.cells[cid].L :* model.recipes[r].Z[f_rows[rr], model.groups[gid].kept]
                        jhat_rr = model.cells[cid].alpha :+ d_rr * model.groups[gid].beta
                        d_tilde = (model.intercept ? d_rr :- model.cells[cid].mean_d : d_rr)
                        s_rr = model.recipes[r].weight[f_rows[rr]] :* (hit[rr, 2] :- jhat_rr)
                        if (cols(mcols) == M) {
                            val = (d_tilde * adj.groups[gid].nu) :+ adj.tilde_lambda[cid, ]
                            out[rr, ] = out[rr, ] + s_rr :* val
                        }
                        else {
                            val = (d_tilde * adj.groups[gid].nu[, mcols]) :+ adj.tilde_lambda[cid, mcols]
                            out[rr, mcols] = out[rr, mcols] + s_rr :* val
                        }
                    }
                    else {
                        s_rr = model.recipes[r].weight[f_rows[rr]] :* (hit[rr, 2] :- model.cells[cid].alpha)
                        if (cols(mcols) == M) out[rr, ] = out[rr, ] + s_rr * adj.tilde_lambda[cid, ]
                        else out[rr, mcols] = out[rr, mcols] + s_rr * adj.tilde_lambda[cid, mcols]
                    }
                }
            }
        }
        a_lo = a_hi + 1
    }
    if (sc == "core") out = model.lambda :* out
    return(out)
}

/* Complete original-unit IF evaluation (phi = direct + j_pullback) for K requested
   cell prediction contrasts on F:
   contrast_weights is n_cells x K, where column m specifies weights omega_{c,m} on
   cell mean predictions mu_c = E_{F, w}[j_hat_c(X) | G = k_c].
   For a single target cell c_g with weight 1, estimate is T_gl and phi is psi_T.
   Every recurring unit's contributions across cells and recipes are summed into one
   row per original F unit before any outer product or covariance step. */
struct cdid_imputation_contrast_if scalar cdid_imp_contrast_if_chunk(
    struct cdid_imputation_model scalar model,
    real matrix contrast_weights,
    real colvector f_rows,
    | string scalar scale)
{
    struct cdid_imputation_contrast_if scalar res
    struct cdid_imputation_adjoint scalar adj
    real scalar B, K, p, cid, gid, ci, r, m, c_lo, omega
    real colvector cell_means, ix_coh, rr, jhat_rr, direct_c
    real matrix lam_mat, nu_raw_mat
    string scalar sc
    sc = (args() < 4 | scale == "" ? "frame" : scale)
    if (sc != "frame" & sc != "core")
        _cdid_imp_fail("IF scale must be frame or core.", 498)
    _cdid_imp_validate_chunk_rows(model, f_rows)
    if (rows(contrast_weights) != model.n_cells | cols(contrast_weights) < 1 | missing(contrast_weights))
        _cdid_imp_fail("contrast_weights must be an n_cells x K finite matrix.", 498)
    B = rows(f_rows)
    K = cols(contrast_weights)
    p = model.p_raw
    cell_means = J(model.n_cells, 1, 0)
    lam_mat = J(model.n_cells, K, 0)
    nu_raw_mat = J(model.n_cells, p * K, 0)
    for (cid=1; cid<=model.n_cells; cid++) {
        cell_means[cid] = model.cells[cid].F_mean_j
        gid = model.cells[cid].group_id
        for (m=1; m<=K; m++) {
            omega = contrast_weights[cid, m]
            if (omega == 0) continue
            if (model.intercept) lam_mat[cid, m] = lam_mat[cid, m] + omega
            if (model.groups[gid].p > 0) {
                c_lo = (m - 1) * p
                nu_raw_mat[cid, c_lo :+ model.groups[gid].kept] =
                    nu_raw_mat[cid, c_lo :+ model.groups[gid].kept] + omega :* model.cells[cid].d_bar_j
            }
        }
    }
    adj = cdid_imputation_adjoint_prepare(model, lam_mat, nu_raw_mat)
    res.f_rows = f_rows
    res.estimate = (contrast_weights' * cell_means)
    res.direct = J(B, K, 0)
    for (cid=1; cid<=model.n_cells; cid++) {
        if (all(contrast_weights[cid, ] :== 0)) continue
        ci = model.cells[cid].cohort_index
        ix_coh = _cdid_imp_index(model.unit_cohort_index[f_rows] :== ci)
        if (rows(ix_coh) == 0) continue
        r = model.cells[cid].recipe
        jhat_rr = cdid_imp_predict_cell_chunk(model, cid, f_rows[ix_coh])
        direct_c = model.recipes[r].weight[f_rows[ix_coh]] :*
            (jhat_rr :- model.cells[cid].F_mean_j) :/ model.cells[cid].F_mass
        res.direct[ix_coh, ] = res.direct[ix_coh, ] + direct_c * contrast_weights[cid, ]
    }
    res.j_pullback = cdid_imputation_adjoint_chunk(model, adj, f_rows, "frame")
    if (sc == "core") {
        res.direct = model.lambda :* res.direct
        res.j_pullback = model.lambda :* res.j_pullback
    }
    res.phi = res.direct + res.j_pullback
    res.scale = sc
    res.normalization_count = (sc == "core" ? model.n_A : model.N)
    return(res)
}

/* Complete fixed-weight log-weight directional derivative d(beta_r, alpha_{C_r}) / d(log w)
   for caller-supplied N x D direction columns on F. Streams in bounded chunks. */
real matrix cdid_imp_logweight_direction(
    struct cdid_imputation_model scalar model,
    real scalar gid,
    real matrix f_directions)
{
    real scalar D_cols, k_keep, B_r, loc, cid, ci, pid, r, n_coh, lo, hi
    real colvector cids, coh_rows, r_chunk, obs_ix, s_rr
    real matrix force_beta, force_alpha, hit, d_rr, d_tilde, db, da
    if (gid < 1 | gid > model.n_groups | gid != floor(gid) |
        rows(f_directions) != model.N | cols(f_directions) < 1 | missing(f_directions))
        _cdid_imp_fail("Invalid group index or f_directions matrix in cdid_imp_logweight_direction.", 498)
    D_cols = cols(f_directions)
    k_keep = model.groups[gid].p
    cids = model.groups[gid].cell_ids
    B_r = rows(cids)
    force_beta = J(k_keep, D_cols, 0)
    force_alpha = J(B_r, D_cols, 0)
    for (loc=1; loc<=B_r; loc++) {
        cid = cids[loc]
        ci = model.cells[cid].cohort_index
        pid = model.cells[cid].pair_id
        r = model.cells[cid].recipe
        coh_rows = _cdid_imp_coh_rows(model, ci)
        n_coh = rows(coh_rows)
        for (lo=1; lo<=n_coh; lo=lo+model.budget.max_chunk_rows) {
            hi = min((n_coh, lo + model.budget.max_chunk_rows - 1))
            r_chunk = coh_rows[lo..hi]
            hit = cdid_imputation_read_pair_chunk(model.reader, model.pop, model.pairs[pid], r_chunk, 1)
            obs_ix = _cdid_imp_index(hit[, 1] :== 1)
            if (rows(obs_ix) == 0) continue
            if (k_keep > 0) {
                d_rr = model.cells[cid].L :* model.recipes[r].Z[r_chunk[obs_ix], model.groups[gid].kept]
                s_rr = model.recipes[r].weight[r_chunk[obs_ix]] :*
                    (hit[obs_ix, 2] :- model.cells[cid].alpha :- d_rr * model.groups[gid].beta)
                d_tilde = (model.intercept ? d_rr :- model.cells[cid].mean_d : d_rr)
                force_beta = force_beta + quadcross(s_rr :* d_tilde, f_directions[r_chunk[obs_ix], ])
            }
            else {
                s_rr = model.recipes[r].weight[r_chunk[obs_ix]] :* (hit[obs_ix, 2] :- model.cells[cid].alpha)
            }
            if (model.intercept) {
                force_alpha[loc, ] = force_alpha[loc, ] +
                    quadcross(s_rr, f_directions[r_chunk[obs_ix], ]) :/ model.cells[cid].weight_sum
            }
        }
    }
    db = (k_keep > 0 ? _cdid_nuisance_gram(model.groups[gid].R, model.groups[gid].pivot, force_beta) : J(0, D_cols, 0))
    da = J(B_r, D_cols, 0)
    if (model.intercept) {
        for (loc=1; loc<=B_r; loc++) {
            cid = cids[loc]
            da[loc, ] = force_alpha[loc, ] - (k_keep > 0 ? model.cells[cid].mean_d * db : J(1, D_cols, 0))
        }
    }
    return(db \ da)
}

end
