*! cdid_catalog.mata -- native comparison catalogue and PT design
*! Mirrors R cdid 1.1.0.0: construct_observed_pairs / materialize_observed_pairs
*! (occupancy only), build_comparison_catalog, catalog_design_handoff and
*! chained_design_core on the calendar path. Stata 14.2 Mata; no frames.
*! Defines structures and functions: run once per Mata session (or after
*! -mata clear-). Standalone: no other cdid Mata file is required.
version 14.2
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_cat(transmorphic matrix V, real scalar idcol)
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

/* ======================================================================
CONVENTIONS (interface contract for the later engine)

Dictionaries
  calendar  T x 1, strictly increasing, nonmissing (|x| < 2^1023), T >= 2
            for the catalogue (T >= 1 for observed pairs).
  cohorts   C x 1, strictly increasing finite values, optionally followed by
            exactly one system missing value (.) = the never-treated cohort
            (R +Inf). It sorts last in Mata exactly as +Inf sorts last in R, so cohort
            index C is the never-treated index when present. Extended missing
            values .a-.z are refused everywhere. 0 is an ordinary finite label.
  The caller sorts dictionaries; positions in them are the "indices" below.

Indices and "none"
  Every index column is a 1-based dictionary/table position. Where R reports
  NA for an index (no reference, no horizon, anchored endpoint, no operator),
  this module stores 0. Index columns are never missing. Values (labels) are
  obtained exactly as calendar[index] / cohorts[index].

Occupancy (catalogue input)
  n x 4 real matrix (cohort index, s index, t index, count) in DICTIONARY
  INDICES, s < t, count a whole number in [0, 2^53-1]; one row per cell.
  cdid_occupancy_from_labels() converts R-style label rows (cohort, s, t,
  count; cohort . = never treated) by exact equality, refusing off-grid
  values; cdid_occupancy_to_labels() is the exact inverse.

Labels in other inputs
  targets (J x 1) and requested_post / declared / design targets (q x 2:
  g, t) are LABELS matched exactly (==) to the dictionaries; never strings.

Numeric comparisons
  Cutoffs u + anticipation < k use an exact TwoSum sign (R exact_cmp_sum).
  On Stata-representable inputs R never refuses; if a + b exceeds Stata's
  double range (Mata returns missing) the true sum exceeds every finite
  cohort, so the sign is +1. Plain comparisons (t >= g, region t < g,
  pt_start windows) are plain double comparisons as in R.

Codes (numeric; name helpers below return the R strings)
  catalogue endpoint role  1 reference 2 pre_private 3 anticipation_conditional
                           4 anticipation_private 5 post_conditional 6 post_private
  coordinate kind          1 reference 2 shared_conditional 3 private
  target_pair status       1 retained 2 empty_pool 3 no_admissible_reference
  target status            1 reference_available 2 no_admissible_reference
  request reason           0 supported 1 no_admissible_reference
                           2 no_target_pair_occupancy 3 no_nonempty_control_pool
  design vertex role       1 reference 2 pre 3 anticipation_nuisance
                           4 anticipation 5 post_nuisance 6 post
  link type                1 pre_only 2 straddling 3 anticipation_only
                           4 anticipation_to_post 5 post_only
  post-target status       1 estimable 2 isolated_unsupported 3 disconnected_component
                           (0 = requested target not in the maintained design)

Refusals
  Invalid input prints its message and exits with return code 3498; the R
  condition class is also stored and returned by cdid_catalog_last_error().
  Unavailable targets are NOT refusals: they are reported per target.

Resources
  struct cdid_catalog_opts: optional metadata bounds (max_candidates,
  max_rows, max_pool_entries; missing = none), the public opt-in work cap
  (work_cap; missing = none) and the numeric-memory policy (mem_limit_bytes,
  default 2^30 = 1 GiB source profile). Exact counts are computed before the
  corresponding allocation; exceeding a bound refuses, nothing is dropped.
====================================================================== */

// ------------------------------------------------------------ refusals
void _cdid_set_lasterr(string scalar cls)
{
    pointer() scalar p
    if ((p = findexternal("cdid_catalog_lasterr")) == NULL) p = crexternal("cdid_catalog_lasterr")
    if (p != NULL) *p = cls
}

string scalar cdid_catalog_last_error()
{
    pointer() scalar p
    if ((p = findexternal("cdid_catalog_lasterr")) == NULL) return("")
    if (eltype(*p) != "string" | rows(*p) != 1 | cols(*p) != 1) return("")
    return(*p)
}

void cdid_catalog_clear_error()
{
    _cdid_set_lasterr("")
}

void _cdid_fail(string scalar cls, string scalar msg)
{
    _cdid_set_lasterr(cls)
    errprintf("%s\n", msg)
    exit(3498)
}

// ------------------------------------------------------------- helpers
real scalar _cdid_any(real matrix x)
{
    return(sum(x :!= 0) > 0)
}

real scalar _cdid_all(real matrix x)
{
    return(sum(x :!= 0) == rows(x) * cols(x))
}

real colvector _cdid_which(real colvector x)
{
    if (rows(x) == 0) return(J(0, 1, .))
    if (!_cdid_any(x)) return(J(0, 1, .))
    return(select((1::rows(x)), x :!= 0))
}

real colvector _cdid_seq(real scalar a, real scalar b)
{
    if (b < a) return(J(0, 1, .))
    return((a::b))
}

real scalar _cdid_is_whole(real matrix x)
{
    if (rows(x) * cols(x) == 0) return(1)
    if (hasmissing(x)) return(0)
    return(_cdid_all(x :== floor(x)))
}

string scalar _cdid_istr(real scalar x)
{
    return(strtrim(strofreal(x, "%21.0f")))
}

// Exact position of x in the strictly increasing nonmissing vector d (0 = absent).
real scalar _cdid_bsearch(real colvector d, real scalar x)
{
    real scalar lo, hi, mid
    if (x >= .) return(0)
    lo = 1
    hi = rows(d)
    while (lo <= hi) {
        mid = floor((lo + hi) / 2)
        if (d[mid] == x) return(mid)
        if (d[mid] < x) lo = mid + 1
        else hi = mid - 1
    }
    return(0)
}

// Exact cohort index of a label; system missing maps to the never-treated index.
real scalar _cdid_cohort_index(real colvector coh, real scalar nt, real scalar x)
{
    if (x == .) return(nt)
    if (x > .) return(0)
    return(_cdid_bsearch(coh[_cdid_seq(1, rows(coh) - (nt > 0))], x))
}

// Exact sign of (a + b) - c; c == . is the never-treated cohort (+Inf).
real scalar _cdid_cmp_sum(real scalar a, real scalar b, real scalar c)
{
    real scalar s, aa, bb, err
    if (c == .) return(-1)
    s = a + b
    if (s >= .) return(1)
    if (s < c) return(-1)
    if (s > c) return(1)
    if (abs(a) >= abs(b)) {
        bb = s - a
        err = (a - (s - bb)) + (b - bb)
    }
    else {
        aa = s - b
        err = (b - (s - aa)) + (a - aa)
    }
    if (err >= .) _cdid_fail("cdid_design_arithmetic_error", "inconclusive floating-point cutoff comparison")
    return(err < 0 ? -1 : (err > 0 ? 1 : 0))
}

// --------------------------------------------------------- dictionaries
void _cdid_check_calendar(real matrix cal, real scalar minlen)
{
    real scalar n
    n = rows(cal)
    if (cols(cal) != 1 | n < minlen)
        _cdid_fail("cdid_design_calendar_error", "calendar must be a column of at least " + _cdid_istr(minlen) + " dates")
    if (hasmissing(cal)) _cdid_fail("cdid_design_calendar_error", "calendar dates must be nonmissing")
    if (n > 1) {
        if (!_cdid_all(cal[|2 \ n|] :> cal[|1 \ n - 1|]))
            _cdid_fail("cdid_design_calendar_error", "calendar must be strictly increasing (sorted, unique)")
    }
}

// Returns the never-treated index (0 if none).
real scalar _cdid_check_cohorts(real matrix coh, real scalar need_finite)
{
    real scalar n, nt, nf
    n = rows(coh)
    if (cols(coh) != 1 | n < 1) _cdid_fail("cdid_design_cohort_error", "cohorts must be a nonempty column")
    nt = (coh[n] == . ? n : 0)
    nf = n - (nt > 0)
    if (nf > 0) {
        if (hasmissing(coh[|1 \ nf|]))
            _cdid_fail("cdid_design_cohort_error", "only one trailing system missing (never treated) is allowed")
        if (nf > 1) {
            if (!_cdid_all(coh[|2 \ nf|] :> coh[|1 \ nf - 1|]))
                _cdid_fail("cdid_design_cohort_error", "finite cohorts must be strictly increasing (sorted, unique)")
        }
    }
    if (coh[n] > .) _cdid_fail("cdid_design_cohort_error", "extended missing values are not cohorts")
    if (need_finite & nf == 0) _cdid_fail("cdid_design_cohort_error", "at least one finite cohort is required")
    return(nt)
}

void _cdid_check_anticipation(real matrix a)
{
    if (rows(a) != 1 | cols(a) != 1)
        _cdid_fail("cdid_design_anticipation_error", "anticipation must be a scalar")
    if (a[1, 1] >= . | a[1, 1] < 0 | a[1, 1] != floor(a[1, 1]))
        _cdid_fail("cdid_design_anticipation_error", "anticipation must be a finite nonnegative whole number")
}

// --------------------------------------------------- resource controls
struct cdid_catalog_opts {
    string scalar conventions
    real scalar   mem_limit_bytes
    real scalar   max_candidates, max_rows, max_pool_entries
    real scalar   work_cap
}

struct cdid_catalog_opts scalar cdid_catalog_default_opts()
{
    struct cdid_catalog_opts scalar o
    o.conventions = "stata-port"
    o.mem_limit_bytes = 2^30
    o.max_candidates = .
    o.max_rows = .
    o.max_pool_entries = .
    o.work_cap = .
    return(o)
}

void _cdid_check_opts(struct cdid_catalog_opts scalar o)
{
    if (o.conventions == "" | strpos(o.conventions, ",") | strpos(o.conventions, char(34)))
        _cdid_fail("cdid_catalog_input_error", "conventions must be a nonempty label without commas or quotes")
    if (o.mem_limit_bytes >= . | o.mem_limit_bytes <= 0)
        _cdid_fail("cdid_catalog_bound_error", "mem_limit_bytes must be a positive finite number")
}

// Refuse before allocating ENTRIES doubles beyond the numeric-memory policy.
void _cdid_mem(struct cdid_catalog_opts scalar o, real scalar entries, string scalar what)
{
    if (entries >= . | 8 * entries > o.mem_limit_bytes)
        _cdid_fail("cdid_catalog_bound_error", what + " needs " + _cdid_istr(8 * entries) +
            " bytes, above the numeric-memory policy of " + _cdid_istr(o.mem_limit_bytes) + " bytes")
}

void _cdid_bound(real scalar value, real scalar bound, string scalar what)
{
    if (bound < . & value > bound)
        _cdid_fail("cdid_catalog_bound_error", "exact " + what + " " + _cdid_istr(value) +
            " exceeds the explicit bound " + _cdid_istr(bound))
}

// ----------------------------------------------------------- structures
// Observed-outcome postings and pair occupancy (R construct_observed_pairs +
// cohort_pair_cells). Units are full-frame indices 1..n_units; units with 0
// or 1 observed outcome stay in the unit table but contribute no pairs.
struct cdid_pairs {
    real colvector calendar, cohorts
    real scalar    nt_index
    real scalar    n_input_rows, n_observed_rows, n_na_rows
    real scalar    n_units, n_absent_units, n_unobserved_units
    real scalar    n_singleton_units, n_pair_eligible_units, n_pairs
    real scalar    n_patterns, was_canonical
    real colvector row_order
    real colvector unit_cohort, unit_nrows, unit_nobs, unit_npairs
    real colvector unit_first, unit_last, unit_single_date, unit_single_y
    real colvector unit_obs_start, unit_obs_end
    real colvector obs_unit, obs_date, obs_y
    real matrix    occupancy
}

// Comparison catalogue (R build_comparison_catalog). Index-valued members
// use 0 for "none". Row r pools: pool_entries[pool_start[p]..+pool_size[p]-1]
// with p = r_pool[r]; pools are typed sorted cohort-index sets.
struct cdid_cat {
    real colvector calendar, cohorts
    real scalar    anticipation, nt_index, nyt, n_dates, n_cohorts, mem_limit_bytes
    string scalar  control_group, conventions
    real colvector cstar
    real scalar    n_cutoff_cmp
    real colvector tg, tg_ref, tg_hfirst, tg_status
    real colvector tg_ntp, tg_ncand, tg_nne, tg_nrows
    real scalar    n_occ, n_pos
    real colvector pos_k, pos_s, pos_t, pos_count, pos_pair
    real colvector pair_s, pair_t, pair_start, pair_end, pair_inf
    real colvector tp_j, tp_g, tp_pair, tp_cell, tp_s, tp_t, tp_count, tp_nother
    real colvector tp_ref, tp_h0, tp_nh, tp_ncand, tp_inf, tp_finlo, tp_finhi
    real colvector tp_nne, tp_orig_ne, tp_nrows, tp_status
    real scalar    n_candidates, n_nonempty, work
    real colvector r_tp, r_g, r_s, r_t, r_ref, r_m, r_inf, r_kmin, r_psize
    real colvector r_orig, r_nh, r_h1, r_h2, r_op, r_pool
    real colvector r_srole, r_trole, r_scoord, r_tcoord
    real scalar    n_pools, pool_entries_preintern
    real colvector pool_size, pool_start, pool_nt, pool_entries
    real colvector op_g, op_pool
    real colvector co_g, co_date, co_kind, co_op
    real colvector rq_g, rq_t, rq_supported, rq_ntouch, rq_reason
    real scalar    n_ctx
    real colvector pos_ctx, pair_ctx, tg_ctx
}

// ------------------------------------------------------- observed pairs
struct cdid_pairs scalar cdid_observed_pairs(real matrix obs, real matrix calendar,
    real matrix cohorts)
{
    return(cdid_observed_pairs_opts(obs, calendar, cohorts, cdid_catalog_default_opts()))
}

// Exact (cohort, s, t, count) occupancy over observed pairs, sorted by (cohort, s, t).
real matrix _cdid_pairs_occupancy_dense(struct cdid_pairs scalar P, real colvector elig,
    real scalar nU, real scalar nT, real scalar nC)
{
    real colvector v, uc, S, T, vals, keep, rows_k
    real matrix Ind, Ok, out
    real scalar k, s, at, np
    v = J(nU * nT, 1, 0)
    v[(P.obs_unit :- 1) :* nT :+ P.obs_date] = J(rows(P.obs_unit), 1, 1)
    Ind = colshape(v, nT)
    v = J(0, 1, 0)
    uc = P.unit_cohort[elig]
    np = nT * (nT - 1) / 2
    S = J(np, 1, 0)
    T = J(np, 1, 0)
    at = 0
    for (s = 1; s < nT; s++) {
        S[|at + 1 \ at + nT - s|] = J(nT - s, 1, s)
        T[|at + 1 \ at + nT - s|] = ((s + 1)::nT)
        at = at + nT - s
    }
    out = J(0, 4, .)
    for (k = 1; k <= nC; k++) {
        rows_k = elig[_cdid_which(uc :== k)]
        if (rows(rows_k) == 0) continue
        Ok = cross(Ind[rows_k, .], Ind[rows_k, .])
        vals = vec(Ok)[(T :- 1) :* nT :+ S]
        keep = _cdid_which(vals :> 0)
        if (rows(keep)) out = out \ (J(rows(keep), 1, k), S[keep], T[keep], vals[keep])
    }
    return(out)
}

// obs: n x 4 (full-frame unit index, calendar index, cohort index, y); a
// missing y (any of . .a-.z) is an unobserved outcome, never an endpoint.
struct cdid_pairs scalar cdid_observed_pairs_opts(real matrix obs, real matrix calendar,
    real matrix cohorts, struct cdid_catalog_opts scalar opt)
{
    struct cdid_pairs scalar P
    real scalar n, nT, nC, nU, nO, ne, W, i, a, u, w, ng, E, at, len, k, cnt, s0, e0, no, same, newg
    real colvector unit, dte, coh, y, ord, dup, oi, kk, single, elig, gfirst, glast, pow2, obs_pow
    real matrix info, oinfo, sig, cells, out, pinfo, cinfo

    _cdid_check_opts(opt)
    _cdid_check_calendar(calendar, 1)
    P.nt_index = _cdid_check_cohorts(cohorts, 0)
    P.calendar = calendar
    P.cohorts = cohorts
    nT = rows(calendar)
    nC = rows(cohorts)
    if (cols(obs) != 4 & rows(obs) > 0)
        _cdid_fail("cdid_pairs_data_error", "observations need 4 columns: unit, date index, cohort index, y")
    n = rows(obs)
    P.n_input_rows = n
    P.occupancy = J(0, 4, .)
    P.row_order = J(0, 1, .)
    P.was_canonical = 1
    P.n_patterns = 0
    if (n == 0) {
        P.n_observed_rows = 0
        P.n_na_rows = 0
        P.n_units = 0
        P.n_absent_units = 0
        P.n_unobserved_units = 0
        P.n_singleton_units = 0
        P.n_pair_eligible_units = 0
        P.n_pairs = 0
        return(P)
    }
    unit = obs[., 1]
    dte = obs[., 2]
    coh = obs[., 3]
    y = obs[., 4]
    if (!_cdid_is_whole(unit) | _cdid_any(unit :< 1))
        _cdid_fail("cdid_pairs_id_error", "unit index must be a nonmissing positive whole number")
    if (!_cdid_is_whole(dte) | _cdid_any((dte :< 1) :| (dte :> nT)))
        _cdid_fail("cdid_pairs_time_error", "date index must be a nonmissing calendar position")
    if (!_cdid_is_whole(coh) | _cdid_any((coh :< 1) :| (coh :> nC)))
        _cdid_fail("cdid_pairs_cohort_error", "cohort index must be a nonmissing cohort position")
    if (n > 1) {
        P.was_canonical = _cdid_all((unit[|2 \ n|] :> unit[|1 \ n - 1|]) :|
            ((unit[|2 \ n|] :== unit[|1 \ n - 1|]) :& (dte[|2 \ n|] :> dte[|1 \ n - 1|])))
    }
    ord = (1::n)
    if (!P.was_canonical) {
        _cdid_mem(opt, 8 * n, "canonical row sort")
        ord = order((unit, dte, (1::n)), (1, 2, 3))
        unit = unit[ord]
        dte = dte[ord]
        coh = coh[ord]
        y = y[ord]
    }
    P.row_order = ord
    if (n > 1) {
        dup = unit[|2 \ n|] :== unit[|1 \ n - 1|]
        if (_cdid_any(dup :& (dte[|2 \ n|] :== dte[|1 \ n - 1|])))
            _cdid_fail("cdid_pairs_duplicate_error", "duplicate (unit, date) rows are invalid")
        if (_cdid_any(dup :& (coh[|2 \ n|] :!= coh[|1 \ n - 1|])))
            _cdid_fail("cdid_pairs_cohort_change_error", "within-unit cohort changes are invalid")
    }
    nU = max(unit)
    _cdid_mem(opt, 12 * nU + 4 * n, "unit table and postings")
    info = _cdid_ps_cat(unit, 1)
    kk = unit[info[., 1]]
    P.unit_nrows = J(nU, 1, 0)
    P.unit_cohort = J(nU, 1, 0)
    P.unit_nrows[kk] = info[., 2] - info[., 1] :+ 1
    P.unit_cohort[kk] = coh[info[., 1]]

    oi = _cdid_which(y :< .)
    nO = rows(oi)
    P.obs_unit = unit[oi]
    P.obs_date = dte[oi]
    P.obs_y = y[oi]
    P.unit_nobs = J(nU, 1, 0)
    P.unit_obs_start = J(nU, 1, 0)
    P.unit_obs_end = J(nU, 1, 0)
    P.unit_first = J(nU, 1, 0)
    P.unit_last = J(nU, 1, 0)
    P.unit_single_date = J(nU, 1, 0)
    P.unit_single_y = J(nU, 1, .)
    if (nO > 0) {
        oinfo = _cdid_ps_cat(P.obs_unit, 1)
        kk = P.obs_unit[oinfo[., 1]]
        P.unit_nobs[kk] = oinfo[., 2] - oinfo[., 1] :+ 1
        P.unit_obs_start[kk] = oinfo[., 1]
        P.unit_obs_end[kk] = oinfo[., 2]
        P.unit_first[kk] = P.obs_date[oinfo[., 1]]
        P.unit_last[kk] = P.obs_date[oinfo[., 2]]
    }
    single = _cdid_which(P.unit_nobs :== 1)
    if (rows(single) > 0) {
        P.unit_single_date[single] = P.unit_first[single]
        P.unit_single_y[single] = P.obs_y[P.unit_obs_start[single]]
    }
    P.unit_npairs = P.unit_nobs :* (P.unit_nobs :- 1) :/ 2
    P.n_pairs = sum(P.unit_npairs)
    if (P.n_pairs > 2^53 - 1)
        _cdid_fail("cdid_pairs_count_error", "total pair count exceeds exact integer arithmetic")
    P.n_observed_rows = nO
    P.n_na_rows = n - nO
    P.n_units = nU
    P.n_absent_units = sum(P.unit_nrows :== 0)
    P.n_unobserved_units = sum((P.unit_nrows :> 0) :& (P.unit_nobs :== 0))
    P.n_singleton_units = rows(single)
    elig = _cdid_which(P.unit_nobs :>= 2)
    ne = rows(elig)
    P.n_pair_eligible_units = ne
    if (ne == 0) return(P)
    if (max(P.obs_y) - min(P.obs_y) >= .) {
        for (i = 1; i <= ne; i++) {
            u = elig[i]
            if (max(P.obs_y[|P.unit_obs_start[u] \ P.unit_obs_end[u]|]) -
                min(P.obs_y[|P.unit_obs_start[u] \ P.unit_obs_end[u]|]) >= .)
                _cdid_fail("cdid_pairs_outcome_error", "an observed outcome pair difference exceeds the Stata double range")
        }
    }

    // Units with the same cohort and observed-date set share one pair set:
    // group them by an exact bitset signature (52 dates per exact word).
    W = ceil(nT / 52)
    _cdid_mem(opt, 2 * ne * (W + 2), "observation-pattern signatures")
    sig = J(ne, W + 2, 0)
    sig[., 1] = P.unit_cohort[elig]
    sig[., W + 2] = elig
    if (W == 1 & _cdid_all(P.unit_nobs[elig] :== nT)) {
        sig[., 2] = J(ne, 1, 2^nT - 1)
    }
    else if (W == 1) {
        pow2 = 2 :^ ((0::(nT - 1)))
        obs_pow = pow2[P.obs_date]
        for (i = 1; i <= ne; i++) {
            u = elig[i]
            if (P.unit_nobs[u] == nT) sig[i, 2] = 2^nT - 1
            else sig[i, 2] = sum(obs_pow[|P.unit_obs_start[u] \ P.unit_obs_end[u]|])
        }
    }
    else {
        for (i = 1; i <= ne; i++) {
            u = elig[i]
            for (a = P.unit_obs_start[u]; a <= P.unit_obs_end[u]; a++) {
                w = floor((P.obs_date[a] - 1) / 52) + 2
                sig[i, w] = sig[i, w] + 2^mod(P.obs_date[a] - 1, 52)
            }
        }
    }
    _sort(sig, (1..(W + 2)))
    if (ne == 1) {
        gfirst = 1
        glast = 1
    }
    else if (W == 1) {
        gfirst = 1 \ (1 :+ _cdid_which((sig[|2, 1 \ ne, 1|] :!= sig[|1, 1 \ ne - 1, 1|]) :|
            (sig[|2, 2 \ ne, 2|] :!= sig[|1, 2 \ ne - 1, 2|])))
        glast = (rows(gfirst) > 1 ? (gfirst[|2 \ rows(gfirst)|] :- 1) \ ne : ne)
    }
    else {
        gfirst = 1 \ (1 :+ _cdid_which(rowsum(sig[|2, 1 \ ne, W + 1|] :!= sig[|1, 1 \ ne - 1, W + 1|]) :> 0))
        glast = (rows(gfirst) > 1 ? (gfirst[|2 \ rows(gfirst)|] :- 1) \ ne : ne)
    }
    ng = rows(gfirst)
    P.n_patterns = ng
    E = 0
    for (i = 1; i <= ng; i++) E = E + P.unit_npairs[sig[gfirst[i], W + 2]]
    _cdid_mem(opt, 8 * E, "pattern pair cells")
    /* Pattern occupancy counts are exact integer sums, so the (cohort, s, t) table is
       formed from per-cohort 0/1 observation matrices when they are small. */
    if (E > 1 & nU * nT <= 4000000 & nT >= 2) {
        P.occupancy = _cdid_pairs_occupancy_dense(P, elig, nU, nT, nC)
        return(P)
    }
    cells = J(E, 4, .)
    at = 0
    for (i = 1; i <= ng; i++) {
        u = sig[gfirst[i], W + 2]
        k = sig[gfirst[i], 1]
        cnt = glast[i] - gfirst[i] + 1
        s0 = P.unit_obs_start[u]
        e0 = P.unit_obs_end[u]
        for (a = s0; a < e0; a++) {
            len = e0 - a
            cells[|at + 1, 1 \ at + len, 4|] = (J(len, 1, k), J(len, 1, P.obs_date[a]),
                P.obs_date[|a + 1 \ e0|], J(len, 1, cnt))
            at = at + len
        }
    }
    _sort(cells, (1, 2, 3))
    if (E == 1) {
        P.occupancy = cells
        return(P)
    }
    gfirst = 1 \ (1 :+ _cdid_which((cells[|2, 1 \ E, 1|] :!= cells[|1, 1 \ E - 1, 1|]) :|
        (cells[|2, 2 \ E, 2|] :!= cells[|1, 2 \ E - 1, 2|]) :|
        (cells[|2, 3 \ E, 3|] :!= cells[|1, 3 \ E - 1, 3|])))
    glast = (rows(gfirst) > 1 ? (gfirst[|2 \ rows(gfirst)|] :- 1) \ E : E)
    no = rows(gfirst)
    if (no == E) {
        P.occupancy = cells
        return(P)
    }
    out = cells[gfirst, .]
    for (i = 1; i <= no; i++) {
        if (glast[i] > gfirst[i]) {
            out[i, 4] = sum(cells[|gfirst[i], 4 \ glast[i], 4|])
        }
    }
    P.occupancy = out
    return(P)
}

// Lazy per-unit pairs (s index, t index, y_s, y_t, y_t - y_s) for one unit.
real matrix cdid_pairs_unit_pairs(struct cdid_pairs scalar P, real scalar u)
{
    real scalar s0, e0, a, at, len
    real matrix out
    if (u >= . | u < 1 | u > P.n_units | u != floor(u)) _cdid_fail("cdid_pairs_id_error", "unit index out of range")
    if (P.unit_nobs[u] < 2) return(J(0, 5, .))
    s0 = P.unit_obs_start[u]
    e0 = P.unit_obs_end[u]
    out = J(P.unit_npairs[u], 5, .)
    at = 0
    for (a = s0; a < e0; a++) {
        len = e0 - a
        out[|at + 1, 1 \ at + len, 5|] = (J(len, 1, P.obs_date[a]), P.obs_date[|a + 1 \ e0|],
            J(len, 1, P.obs_y[a]), P.obs_y[|a + 1 \ e0|], P.obs_y[|a + 1 \ e0|] :- P.obs_y[a])
        at = at + len
    }
    return(out)
}

// ------------------------------------------------------------ catalogue
struct cdid_cat scalar cdid_catalog_build(real matrix calendar, real matrix cohorts,
    real matrix targets, real matrix anticipation, real matrix occupancy,
    string scalar control_group, real matrix requested_post)
{
    return(cdid_catalog_build_opts(calendar, cohorts, targets, anticipation, occupancy,
        control_group, requested_post, cdid_catalog_default_opts()))
}

struct cdid_cat scalar cdid_catalog_build_opts(real matrix calendar, real matrix cohorts,
    real matrix targets, real matrix anticipation, real matrix occupancy,
    string scalar control_group, real matrix requested_post, struct cdid_catalog_opts scalar opt)
{
    struct cdid_cat scalar C
    real colvector tg_ctx
    real matrix sub
    _cdid_cat_spec(C, calendar, cohorts, targets, anticipation, control_group, opt)
    _cdid_cat_cells(C, occupancy, opt)
    // Restricted eligibility contexts (missingreference(exclude)): each context
    // reads the owner occupancy minus the pair counts of its excluded units.
    if (findexternal("__cdid_mref_state") != NULL) {
        if (cdid_mrca_catalog(C.calendar, C.cohorts, C.nt_index, C.tg, tg_ctx, sub))
            _cdid_cat_cells_ctx(C, occupancy, tg_ctx, sub, opt)
    }
    _cdid_cat_requests(C, requested_post)
    _cdid_cat_cutoffs(C)
    _cdid_cat_target_pairs(C)
    _cdid_cat_rows(C, opt)
    _cdid_cat_pools(C, opt)
    _cdid_cat_roles(C)
    _cdid_cat_tables(C)
    _cdid_cat_classify_requests(C)
    return(C)
}

void _cdid_cat_spec(struct cdid_cat scalar C, real matrix calendar, real matrix cohorts,
    real matrix targets, real matrix anticipation, string scalar control_group,
    struct cdid_catalog_opts scalar opt)
{
    real scalar j, k, nJ
    _cdid_check_opts(opt)
    _cdid_check_calendar(calendar, 2)
    C.nt_index = _cdid_check_cohorts(cohorts, 1)
    _cdid_check_anticipation(anticipation)
    C.calendar = calendar
    C.cohorts = cohorts
    C.anticipation = anticipation[1, 1]
    C.n_dates = rows(calendar)
    C.n_cohorts = rows(cohorts)
    if (control_group != "nevertreated" & control_group != "notyettreated")
        _cdid_fail("cdid_catalog_input_error", "control_group must be nevertreated or notyettreated")
    C.control_group = control_group
    C.nyt = (control_group == "notyettreated")
    C.conventions = opt.conventions
    C.mem_limit_bytes = opt.mem_limit_bytes
    nJ = rows(targets)
    if (cols(targets) != 1 | nJ < 1 | hasmissing(targets))
        _cdid_fail("cdid_catalog_target_error", "targets must be a nonempty column of finite cohorts")
    C.tg = J(nJ, 1, 0)
    for (j = 1; j <= nJ; j++) {
        k = _cdid_cohort_index(C.cohorts, C.nt_index, targets[j, 1])
        if (k == 0) _cdid_fail("cdid_catalog_target_error", "every target must belong to the cohort universe")
        C.tg[j] = k
    }
    C.tg = sort(C.tg, 1)
    if (nJ > 1) {
        if (_cdid_any(C.tg[|2 \ nJ|] :== C.tg[|1 \ nJ - 1|]))
            _cdid_fail("cdid_catalog_target_error", "targets must be unique")
    }
}

void _cdid_cat_cells(struct cdid_cat scalar C, real matrix occ, struct cdid_catalog_opts scalar opt)
{
    real scalar n, np, npair
    real matrix so, pos
    real colvector newp

    n = rows(occ)
    if (n > 0 & cols(occ) != 4)
        _cdid_fail("cdid_catalog_occupancy_error", "occupancy must have exactly the columns cohort, s, t, count")
    C.n_occ = n
    C.n_ctx = 0
    C.pos_ctx = C.pair_ctx = C.tg_ctx = J(0, 1, .)
    pos = J(0, 4, .)
    if (n > 0) {
        if (!_cdid_is_whole(occ))
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy entries must be nonmissing whole numbers")
        if (_cdid_any((occ[., 1] :< 1) :| (occ[., 1] :> C.n_cohorts)))
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy cohort outside the cohort universe")
        if (_cdid_any((occ[., 2] :< 1) :| (occ[., 2] :> C.n_dates) :| (occ[., 3] :< 1) :| (occ[., 3] :> C.n_dates)))
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy dates must be stored dates of the full calendar")
        if (_cdid_any((occ[., 4] :< 0) :| (occ[., 4] :> 2^53 - 1)))
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy counts must be whole numbers in [0, 2^53 - 1]")
        if (_cdid_any(occ[., 2] :>= occ[., 3]))
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy pairs require s < t")
        _cdid_mem(opt, 12 * n, "occupancy sort")
        so = sort(occ, (2, 3, 1))
        if (n > 1) {
            if (_cdid_any((so[|2, 1 \ n, 1|] :== so[|1, 1 \ n - 1, 1|]) :& (so[|2, 2 \ n, 2|] :== so[|1, 2 \ n - 1, 2|]) :&
                (so[|2, 3 \ n, 3|] :== so[|1, 3 \ n - 1, 3|])))
                _cdid_fail("cdid_catalog_duplicate_error", "duplicate occupancy (cohort, s, t); supply one count per cell")
        }
        pos = select(so, so[., 4] :> 0)
    }
    np = rows(pos)
    C.n_pos = np
    C.pos_k = pos[., 1]
    C.pos_s = pos[., 2]
    C.pos_t = pos[., 3]
    C.pos_count = pos[., 4]
    C.pos_pair = J(0, 1, .)
    C.pair_s = J(0, 1, .)
    C.pair_t = J(0, 1, .)
    C.pair_start = J(0, 1, .)
    C.pair_end = J(0, 1, .)
    C.pair_inf = J(0, 1, .)
    if (np == 0) return
    newp = J(np, 1, 1)
    if (np > 1) newp[|2 \ np|] = (pos[|2, 2 \ np, 2|] :!= pos[|1, 2 \ np - 1, 2|]) :| (pos[|2, 3 \ np, 3|] :!= pos[|1, 3 \ np - 1, 3|])
    C.pos_pair = runningsum(newp)
    C.pair_start = _cdid_which(newp)
    npair = rows(C.pair_start)
    C.pair_end = J(npair, 1, np)
    if (npair > 1) C.pair_end[|1 \ npair - 1|] = C.pair_start[|2 \ npair|] :- 1
    C.pair_s = C.pos_s[C.pair_start]
    C.pair_t = C.pos_t[C.pair_start]
    C.pair_inf = (C.pos_k[C.pair_end] :== C.nt_index)
}

// Restricted eligibility contexts: every context receives every owner cell
// with its excluded units' pair counts removed; positive cells are sorted by
// (context, s, t, cohort), context 0 being the owner occupancy, and a pair is
// (context, s, t). A target reads only the cells of its own context.
void _cdid_cat_cells_ctx(struct cdid_cat scalar C, real matrix occ, real colvector tg_ctx,
    real matrix sub, struct cdid_catalog_opts scalar opt)
{
    real scalar n, nctx, c, i, np, npair
    real matrix so, all, pos
    real colvector hit, newp
    n = rows(occ)
    nctx = max((0 \ tg_ctx \ (rows(sub) ? sub[., 1] : J(0, 1, .))))
    C.n_ctx = nctx
    C.tg_ctx = tg_ctx
    if (nctx == 0) {
        C.pos_ctx = J(C.n_pos, 1, 0)
        C.pair_ctx = J(rows(C.pair_s), 1, 0)
        return
    }
    _cdid_mem(opt, 6 * n * (nctx + 1), "eligibility occupancy cells")
    so = (n ? sort(occ, (2, 3, 1)) : J(0, 4, .))
    all = (J(n, 1, 0), so)
    for (c = 1; c <= nctx; c++) all = all \ (J(n, 1, c), so)
    for (i = 1; i <= rows(sub); i++) {
        if (sub[i, 1] < 1 | sub[i, 1] > nctx)
            _cdid_fail("cdid_catalog_occupancy_error", "eligibility occupancy context outside the restricted contexts")
        hit = _cdid_which((so[., 1] :== sub[i, 2]) :& (so[., 2] :== sub[i, 3]) :& (so[., 3] :== sub[i, 4]))
        if (rows(hit) != 1)
            _cdid_fail("cdid_catalog_occupancy_error", "eligibility occupancy subtracts a pair cell absent from the object occupancy")
        all[sub[i, 1] * n + hit[1], 5] = all[sub[i, 1] * n + hit[1], 5] - sub[i, 5]
    }
    if (_cdid_any(all[., 5] :< 0))
        _cdid_fail("cdid_catalog_occupancy_error", "eligibility occupancy subtracts more units than the object occupancy holds")
    pos = select(all, all[., 5] :> 0)
    pos = (rows(pos) ? sort(pos, (1, 3, 4, 2)) : J(0, 5, .))
    np = rows(pos)
    C.n_pos = np
    C.pos_ctx = pos[., 1]
    C.pos_k = pos[., 2]
    C.pos_s = pos[., 3]
    C.pos_t = pos[., 4]
    C.pos_count = pos[., 5]
    C.pos_pair = J(0, 1, .)
    C.pair_s = J(0, 1, .)
    C.pair_t = J(0, 1, .)
    C.pair_start = J(0, 1, .)
    C.pair_end = J(0, 1, .)
    C.pair_inf = J(0, 1, .)
    C.pair_ctx = J(0, 1, .)
    if (np == 0) return
    newp = J(np, 1, 1)
    if (np > 1) newp[|2 \ np|] = (pos[|2, 1 \ np, 1|] :!= pos[|1, 1 \ np - 1, 1|]) :|
        (pos[|2, 3 \ np, 3|] :!= pos[|1, 3 \ np - 1, 3|]) :| (pos[|2, 4 \ np, 4|] :!= pos[|1, 4 \ np - 1, 4|])
    C.pos_pair = runningsum(newp)
    C.pair_start = _cdid_which(newp)
    npair = rows(C.pair_start)
    C.pair_end = J(npair, 1, np)
    if (npair > 1) C.pair_end[|1 \ npair - 1|] = C.pair_start[|2 \ npair|] :- 1
    C.pair_s = C.pos_s[C.pair_start]
    C.pair_t = C.pos_t[C.pair_start]
    C.pair_ctx = C.pos_ctx[C.pair_start]
    C.pair_inf = (C.pos_k[C.pair_end] :== C.nt_index)
}

// cstar[d] = first cohort index k with calendar[d] + anticipation < cohorts[k]
// (n_cohorts + 1 if none): exact two-pointer scan as in R. Reference r_g is the
// last d with cstar[d] <= index(g); hfirst is the first stored date >= g.
void _cdid_cat_cutoffs(struct cdid_cat scalar C)
{
    real scalar d, ci, kk, below, j, nJ
    C.cstar = J(C.n_dates, 1, 0)
    ci = 1
    C.n_cutoff_cmp = 0
    for (d = 1; d <= C.n_dates; d++) {
        while (1) {
            C.n_cutoff_cmp = C.n_cutoff_cmp + 1
            kk = (ci < C.n_cohorts ? ci : C.n_cohorts)
            below = (_cdid_cmp_sum(C.calendar[d], C.anticipation, C.cohorts[kk]) < 0)
            if (ci > C.n_cohorts | below) break
            ci++
        }
        C.cstar[d] = ci
    }
    nJ = rows(C.tg)
    C.tg_ref = J(nJ, 1, 0)
    C.tg_hfirst = J(nJ, 1, 0)
    for (j = 1; j <= nJ; j++) {
        C.tg_ref[j] = sum(C.cstar :<= C.tg[j])
        C.tg_hfirst[j] = sum(C.calendar :< C.cohorts[C.tg[j]]) + 1
    }
    C.tg_status = 1 :+ (C.tg_ref :== 0)
    if (C.n_ctx > 0) C.tg_status = C.tg_status + 2 :* ((C.tg_ref :> 0) :& (C.tg_ctx :> 0))
}

void _cdid_cat_target_pairs(struct cdid_cat scalar C)
{
    real scalar ntp
    real colvector tpos, tj, idx
    real matrix M

    tpos = J(C.n_cohorts, 1, 0)
    tpos[C.tg] = (1::rows(C.tg))
    idx = J(0, 1, .)
    if (C.n_pos > 0) {
        tj = tpos[C.pos_k]
        if (C.n_ctx > 0) idx = _cdid_which((tj :> 0) :& (C.pos_ctx :== C.tg_ctx[rowmax((tj, J(rows(tj), 1, 1)))]))
        else idx = _cdid_which(tj :> 0)
    }
    ntp = rows(idx)
    if (ntp == 0) {
        C.tp_j = J(0, 1, .)
        C.tp_g = J(0, 1, .)
        C.tp_pair = J(0, 1, .)
        C.tp_cell = J(0, 1, .)
        C.tp_s = J(0, 1, .)
        C.tp_t = J(0, 1, .)
        C.tp_count = J(0, 1, .)
        C.tp_nother = J(0, 1, .)
        C.tp_ref = J(0, 1, .)
        C.tp_h0 = J(0, 1, .)
        C.tp_nh = J(0, 1, .)
        C.tp_ncand = J(0, 1, .)
        C.tp_inf = J(0, 1, .)
        C.tp_finlo = J(0, 1, .)
        C.tp_finhi = J(0, 1, .)
        return
    }
    M = (tj[idx], C.pos_pair[idx], idx)
    M = M[order(M, (1, 2)), .]
    C.tp_j = M[., 1]
    C.tp_pair = M[., 2]
    C.tp_cell = M[., 3]
    C.tp_g = C.tg[C.tp_j, 1]
    C.tp_s = C.pair_s[C.tp_pair, 1]
    C.tp_t = C.pair_t[C.tp_pair, 1]
    C.tp_count = C.pos_count[C.tp_cell, 1]
    C.tp_nother = C.pair_end[C.tp_pair, 1] - C.pair_start[C.tp_pair, 1]
    C.tp_ref = C.tg_ref[C.tp_j, 1]
    C.tp_h0 = rowmax((C.tg_hfirst[C.tp_j, 1], C.tp_t))
    C.tp_nh = C.nyt :* rowmax((J(ntp, 1, 0), C.n_dates :- C.tp_h0 :+ 1))
    C.tp_ncand = (C.tp_ref :> 0) :* (1 :+ C.tp_nh)
    C.tp_inf = C.pair_inf[C.tp_pair, 1]
    // Later finite cohorts of a target pair are the cells after the target's
    // own cell in its pair (cells are sorted by cohort), minus a final NT cell.
    C.tp_finlo = C.tp_cell :+ 1
    C.tp_finhi = C.pair_end[C.tp_pair, 1] - C.tp_inf
}

// Candidates: per target pair the original NYT/NT pool (h = 0) and, under NYT,
// the fixed-horizon reevaluations h = max(first date >= g, t), ..., T. A pool is
// the suffix of size m of the later finite cohorts plus NT when present at the
// pair. Rows = distinct nonempty (target pair, m), ordered by (pair, -m, h).
void _cdid_cat_rows(struct cdid_cat scalar C, struct cdid_catalog_opts scalar opt)
{
    real scalar ntp, nc, q, n, at, f, i, cc, k, nne, nr, fpos
    real colvector ctp, ch, cm, fin, lo, ne, idx, newr, kpos, w

    ntp = rows(C.tp_j)
    nc = sum(C.tp_ncand)
    C.n_candidates = nc
    _cdid_bound(nc, opt.max_candidates, "candidate count")
    _cdid_bound(nc, opt.work_cap, "catalogue work (candidates)")
    _cdid_mem(opt, 8 * nc, "candidate arrays")
    C.tp_nne = J(ntp, 1, 0)
    C.tp_orig_ne = J(ntp, 1, 0)
    C.tp_nrows = J(ntp, 1, 0)
    C.n_nonempty = 0
    ctp = J(nc, 1, 0)
    ch = J(nc, 1, 0)
    cm = J(nc, 1, 0)
    at = 0
    for (q = 1; q <= ntp; q++) {
        n = C.tp_ncand[q]
        if (n == 0) continue
        ctp[|at + 1 \ at + n|] = J(n, 1, q)
        if (n > 1) ch[|at + 2 \ at + n|] = (C.tp_h0[q]::C.n_dates)
        if (C.nyt) {
            f = C.tp_finhi[q] - C.tp_finlo[q] + 1
            if (f > 0) {
                fin = C.pos_k[|C.tp_finlo[q] \ C.tp_finhi[q]|]
                lo = C.cstar[C.tp_t[q]]
                if (n > 1) lo = lo \ C.cstar[|C.tp_h0[q] \ C.n_dates|]
                // Both lists are sorted. Count the upper suffix exactly as
                // R findInterval(), without a cohort-by-horizon matrix.
                fpos = 1
                for (i = 1; i <= n; i++) {
                    while (fpos <= f) {
                        if (fin[fpos] >= lo[i]) break
                        fpos++
                    }
                    cm[at + i] = f - fpos + 1
                }
            }
        }
        at = at + n
    }
    C.r_tp = J(0, 1, .)
    C.r_m = J(0, 1, .)
    C.r_orig = J(0, 1, .)
    C.r_nh = J(0, 1, .)
    C.r_h1 = J(0, 1, .)
    C.r_h2 = J(0, 1, .)
    if (nc > 0) {
        ne = (cm :> 0) :| C.tp_inf[ctp, 1]
        for (i = 1; i <= nc; i++) {
            if (ne[i]) {
                C.tp_nne[ctp[i]] = C.tp_nne[ctp[i]] + 1
                if (ch[i] == 0) C.tp_orig_ne[ctp[i]] = 1
            }
        }
        idx = _cdid_which(ne)
        nne = rows(idx)
        C.n_nonempty = nne
        if (nne > 0) {
            idx = idx[order((ctp[idx], -cm[idx], ch[idx]), (1, 2, 3))]
            newr = J(nne, 1, 1)
            if (nne > 1) newr[|2 \ nne|] = ((ctp[idx[|2 \ nne|]] :!= ctp[idx[|1 \ nne - 1|]]) :|
                (cm[idx[|2 \ nne|]] :!= cm[idx[|1 \ nne - 1|]]))
            nr = sum(newr)
            _cdid_bound(nr, opt.max_rows, "retained row count")
            C.r_tp = J(nr, 1, 0)
            C.r_m = J(nr, 1, 0)
            C.r_orig = J(nr, 1, 0)
            C.r_nh = J(nr, 1, 0)
            C.r_h1 = J(nr, 1, 0)
            C.r_h2 = J(nr, 1, 0)
            k = 0
            for (i = 1; i <= nne; i++) {
                cc = idx[i]
                if (newr[i]) {
                    k++
                    C.r_tp[k] = ctp[cc]
                    C.r_m[k] = cm[cc]
                }
                if (ch[cc] == 0) C.r_orig[k] = 1
                else {
                    C.r_nh[k] = C.r_nh[k] + 1
                    if (C.r_h1[k] == 0) C.r_h1[k] = ch[cc]
                    C.r_h2[k] = ch[cc]
                }
            }
        }
    }
    nr = rows(C.r_tp)
    C.r_g = C.tp_g[C.r_tp, 1]
    C.r_s = C.tp_s[C.r_tp, 1]
    C.r_t = C.tp_t[C.r_tp, 1]
    C.r_ref = C.tp_ref[C.r_tp, 1]
    C.r_inf = C.tp_inf[C.r_tp, 1]
    C.r_psize = C.r_m + C.r_inf
    C.r_kmin = J(nr, 1, 0)
    w = _cdid_which(C.r_m :> 0)
    if (rows(w) > 0) {
        kpos = C.tp_finhi[C.r_tp[w], 1] - C.r_m[w] :+ 1
        C.r_kmin[w] = C.pos_k[kpos]
    }
    for (i = 1; i <= nr; i++) C.tp_nrows[C.r_tp[i]] = C.tp_nrows[C.r_tp[i]] + 1
    C.tp_status = (C.tp_ref :== 0) :* 3 + (C.tp_ref :> 0) :* (2 :- C.tp_orig_ne)
}

// Lexicographic order of sorted index sequences with a proper prefix first:
// exactly R's radix sort of zero-padded comma-joined index keys.
real scalar _cdid_lexless(real colvector e, real scalar sa, real scalar la,
    real scalar sb, real scalar lb)
{
    real scalar i, n
    n = (la < lb ? la : lb)
    for (i = 0; i < n; i++) {
        if (e[sa + i] != e[sb + i]) return(e[sa + i] < e[sb + i])
    }
    return(la < lb)
}

// Stable bottom-up merge sort of sequences (start st, length ln) in e.
real colvector _cdid_lex_order(real colvector e, real colvector st, real colvector ln)
{
    real scalar n, w, lo, mid, hi, i, j, k
    real colvector a, b
    n = rows(st)
    if (n == 0) return(J(0, 1, .))
    a = (1::n)
    b = J(n, 1, 0)
    for (w = 1; w < n; w = 2 * w) {
        for (lo = 1; lo <= n; lo = lo + 2 * w) {
            mid = (lo + w - 1 < n ? lo + w - 1 : n)
            hi = (lo + 2 * w - 1 < n ? lo + 2 * w - 1 : n)
            i = lo
            j = mid + 1
            k = lo
            while (i <= mid & j <= hi) {
                if (_cdid_lexless(e, st[a[j]], ln[a[j]], st[a[i]], ln[a[i]])) {
                    b[k] = a[j]
                    j++
                }
                else {
                    b[k] = a[i]
                    i++
                }
                k++
            }
            while (i <= mid) {
                b[k] = a[i]
                i++
                k++
            }
            while (j <= hi) {
                b[k] = a[j]
                j++
                k++
            }
        }
        a = b
    }
    return(a)
}

// Intern pools: rows with identical cohort sets (exact 52-bit bitset words)
// share one pool; pool ids follow R's lexicographic key order.
void _cdid_cat_pools(struct cdid_cat scalar C, struct cdid_catalog_opts scalar opt)
{
    real scalar nr, W, r, a, k, w, i, ng, newg, at, hi, p, tot
    real colvector grp, rep, ulen, ustart, uent, po, gpid
    real matrix sig

    nr = rows(C.r_tp)
    C.pool_entries_preintern = sum(C.r_psize)
    _cdid_bound(C.pool_entries_preintern, opt.max_pool_entries, "pool entries")
    C.work = C.n_candidates + C.n_nonempty + 2 * nr + C.pool_entries_preintern
    _cdid_bound(C.work, opt.work_cap, "catalogue work")
    C.n_pools = 0
    C.pool_size = J(0, 1, .)
    C.pool_start = J(0, 1, .)
    C.pool_nt = J(0, 1, .)
    C.pool_entries = J(0, 1, .)
    C.r_pool = J(0, 1, .)
    if (nr == 0) return
    W = ceil(C.n_cohorts / 52)
    _cdid_mem(opt, 2 * nr * (W + 1) + 3 * C.pool_entries_preintern, "pool interning")
    sig = J(nr, W + 1, 0)
    sig[., W + 1] = (1::nr)
    for (r = 1; r <= nr; r++) {
        if (C.r_m[r] > 0) {
            for (a = C.tp_finhi[C.r_tp[r]] - C.r_m[r] + 1; a <= C.tp_finhi[C.r_tp[r]]; a++) {
                k = C.pos_k[a]
                w = floor((k - 1) / 52) + 1
                sig[r, w] = sig[r, w] + 2^mod(k - 1, 52)
            }
        }
        if (C.r_inf[r]) {
            k = C.nt_index
            w = floor((k - 1) / 52) + 1
            sig[r, w] = sig[r, w] + 2^mod(k - 1, 52)
        }
    }
    _sort(sig, (1..(W + 1)))
    grp = J(nr, 1, 0)
    rep = J(nr, 1, 0)
    ng = 0
    for (i = 1; i <= nr; i++) {
        newg = 1
        if (i > 1) newg = (sig[i, (1..W)] != sig[i - 1, (1..W)])
        if (newg) {
            ng++
            rep[ng] = sig[i, W + 1]
        }
        grp[sig[i, W + 1]] = ng
    }
    rep = rep[|1 \ ng|]
    ulen = C.r_psize[rep]
    ustart = runningsum(ulen) - ulen :+ 1
    tot = sum(ulen)
    uent = J(tot, 1, 0)
    for (i = 1; i <= ng; i++) {
        r = rep[i]
        at = ustart[i] - 1
        if (C.r_m[r] > 0) {
            hi = C.tp_finhi[C.r_tp[r]]
            uent[|at + 1 \ at + C.r_m[r]|] = C.pos_k[|hi - C.r_m[r] + 1 \ hi|]
            at = at + C.r_m[r]
        }
        if (C.r_inf[r]) uent[at + 1] = C.nt_index
    }
    po = _cdid_lex_order(uent, ustart, ulen)
    gpid = J(ng, 1, 0)
    gpid[po] = (1::ng)
    C.r_pool = gpid[grp, 1]
    C.n_pools = ng
    C.pool_size = ulen[po]
    C.pool_start = runningsum(C.pool_size) - C.pool_size :+ 1
    C.pool_entries = J(tot, 1, 0)
    for (p = 1; p <= ng; p++) {
        a = ustart[po[p]]
        hi = a + ulen[po[p]] - 1
        C.pool_entries[|C.pool_start[p] \ C.pool_start[p] + C.pool_size[p] - 1|] = uent[|a \ hi|]
    }
    C.pool_nt = (C.pool_entries[C.pool_start + C.pool_size :- 1] :== C.nt_index)
}

real colvector _cdid_role_kind(real colvector role)
{
    return((role :== 1) + 2 :* ((role :== 3) :| (role :== 5)) + 3 :* ((role :== 2) :| (role :== 4) :| (role :== 6)))
}

// Conditional endpoint role of date u (per row): reference / pre_private by
// the exact reference; otherwise eligible iff the pool lies in E_u(anticipation),
// i.e. its smallest finite cohort index is >= cstar[u] (NT always clean).
real colvector _cdid_cat_role(struct cdid_cat scalar C, real colvector u)
{
    real colvector pre, isref, elig, ant
    pre = (u :< C.r_ref)
    isref = (u :== C.r_ref)
    elig = (C.r_kmin :== 0) :| (C.r_kmin :>= C.cstar[u])
    ant = (C.calendar[u] :< C.cohorts[C.r_g])
    return(isref + 2 :* pre + ((!pre) :& (!isref)) :* (ant :* (3 :+ (!elig)) + (!ant) :* (5 :+ (!elig))))
}

void _cdid_cat_roles(struct cdid_cat scalar C)
{
    real scalar nr, i, r, nop, newop, nref, nE, e, nco, newc
    real colvector o, opg, opp, ks, kt, id, rr
    real matrix E

    nr = rows(C.r_tp)
    C.r_op = J(nr, 1, 0)
    C.op_g = J(0, 1, .)
    C.op_pool = J(0, 1, .)
    C.r_srole = J(0, 1, .)
    C.r_trole = J(0, 1, .)
    C.r_scoord = J(0, 1, .)
    C.r_tcoord = J(0, 1, .)
    if (nr > 0) {
        o = order((C.r_g, C.r_pool, (1::nr)), (1, 2, 3))
        opg = J(nr, 1, 0)
        opp = J(nr, 1, 0)
        nop = 0
        for (i = 1; i <= nr; i++) {
            r = o[i]
            newop = 1
            if (i > 1) newop = (C.r_g[r] != C.r_g[o[i - 1]] | C.r_pool[r] != C.r_pool[o[i - 1]])
            if (newop) {
                nop++
                opg[nop] = C.r_g[r]
                opp[nop] = C.r_pool[r]
            }
            C.r_op[r] = nop
        }
        C.op_g = opg[|1 \ nop|]
        C.op_pool = opp[|1 \ nop|]
        C.r_srole = _cdid_cat_role(C, C.r_s)
        C.r_trole = _cdid_cat_role(C, C.r_t)
    }
    // Coordinates: shared reference/conditional vertices per (g, u); private
    // vertices per (g, operator, u); every referenced target is anchored.
    rr = _cdid_which(C.tg_ref :> 0)
    nref = rows(rr)
    E = J(0, 4, .)
    if (nref > 0) E = (C.tg[rr], C.tg_ref[rr], J(nref, 1, 1), J(nref, 1, 0))
    if (nr > 0) {
        ks = _cdid_role_kind(C.r_srole)
        kt = _cdid_role_kind(C.r_trole)
        E = E \ (C.r_g, C.r_s, ks, (ks :== 3) :* C.r_op) \ (C.r_g, C.r_t, kt, (kt :== 3) :* C.r_op)
    }
    nE = rows(E)
    C.co_g = J(0, 1, .)
    C.co_date = J(0, 1, .)
    C.co_kind = J(0, 1, .)
    C.co_op = J(0, 1, .)
    if (nE == 0) return
    o = order((E, (1::nE)), (1, 2, 3, 4, 5))
    id = J(nE, 1, 0)
    C.co_g = J(nE, 1, 0)
    C.co_date = J(nE, 1, 0)
    C.co_kind = J(nE, 1, 0)
    C.co_op = J(nE, 1, 0)
    nco = 0
    for (i = 1; i <= nE; i++) {
        e = o[i]
        newc = 1
        if (i > 1) newc = (E[e, .] != E[o[i - 1], .])
        if (newc) {
            nco++
            C.co_g[nco] = E[e, 1]
            C.co_date[nco] = E[e, 2]
            C.co_kind[nco] = E[e, 3]
            C.co_op[nco] = E[e, 4]
        }
        id[e] = nco
    }
    C.co_g = C.co_g[|1 \ nco|]
    C.co_date = C.co_date[|1 \ nco|]
    C.co_kind = C.co_kind[|1 \ nco|]
    C.co_op = C.co_op[|1 \ nco|]
    if (nr > 0) {
        C.r_scoord = id[|nref + 1 \ nref + nr|]
        C.r_tcoord = id[|nref + nr + 1 \ nref + 2 * nr|]
    }
}

void _cdid_cat_tables(struct cdid_cat scalar C)
{
    real scalar j, nJ
    real colvector mj
    nJ = rows(C.tg)
    C.tg_ntp = J(nJ, 1, 0)
    C.tg_ncand = J(nJ, 1, 0)
    C.tg_nne = J(nJ, 1, 0)
    C.tg_nrows = J(nJ, 1, 0)
    for (j = 1; j <= nJ; j++) {
        mj = (C.tp_j :== j)
        C.tg_ntp[j] = sum(mj)
        C.tg_ncand[j] = sum(mj :* C.tp_ncand)
        C.tg_nne[j] = sum(mj :* C.tp_nne)
        C.tg_nrows[j] = sum(mj :* C.tp_nrows)
    }
}

void _cdid_cat_requests(struct cdid_cat scalar C, real matrix req)
{
    real scalar q, i, g, t
    real matrix sq
    q = rows(req)
    if (q > 0 & cols(req) != 2)
        _cdid_fail("cdid_catalog_request_error", "requested_post must have exactly the two columns g and t")
    C.rq_g = J(q, 1, 0)
    C.rq_t = J(q, 1, 0)
    for (i = 1; i <= q; i++) {
        if (req[i, 1] >= . | req[i, 2] >= .)
            _cdid_fail("cdid_catalog_request_error", "requested post targets require finite numeric g and t")
        g = _cdid_cohort_index(C.cohorts, C.nt_index, req[i, 1])
        if (g == 0 | !_cdid_any(C.tg :== g))
            _cdid_fail("cdid_catalog_request_error", "a requested post target cohort must be one of the targets")
        t = _cdid_bsearch(C.calendar, req[i, 2])
        if (t == 0)
            _cdid_fail("cdid_catalog_request_error", "a requested post date must be a stored date of the full calendar")
        if (req[i, 2] < req[i, 1])
            _cdid_fail("cdid_catalog_request_error", "requested post targets require t >= g")
        C.rq_g[i] = g
        C.rq_t[i] = t
    }
    if (q > 1) {
        sq = sort((C.rq_g, C.rq_t), (1, 2))
        if (_cdid_any((sq[|2, 1 \ q, 1|] :== sq[|1, 1 \ q - 1, 1|]) :& (sq[|2, 2 \ q, 2|] :== sq[|1, 2 \ q - 1, 2|])))
            _cdid_fail("cdid_catalog_duplicate_error", "duplicate requested post target")
    }
}

// Catalogue support of (g, t) cohort/date indices: columns supported (0/1),
// rows touching (all post_conditional), reason code (0 if supported).
real matrix _cdid_cat_classify(struct cdid_cat scalar C, real colvector gq, real colvector tq)
{
    real scalar q, i, g, t, touch, cond, hit, reason
    real colvector ms, mt
    real matrix out
    q = rows(gq)
    out = J(q, 3, 0)
    for (i = 1; i <= q; i++) {
        g = gq[i]
        t = tq[i]
        ms = (C.r_g :== g) :& (C.r_s :== t)
        mt = (C.r_g :== g) :& (C.r_t :== t)
        touch = sum(ms) + sum(mt)
        cond = sum(ms :& (C.r_srole :== 5)) + sum(mt :& (C.r_trole :== 5))
        if (touch != cond)
            _cdid_fail("cdid_catalog_internal_error", "internal role invariant failed: a retained row touches a post date privately")
        if (cond > 0) {
            out[i, .] = (1, cond, 0)
            continue
        }
        hit = _cdid_any((C.tp_g :== g) :& ((C.tp_s :== t) :| (C.tp_t :== t)))
        reason = (hit ? 3 : 2)
        if (reason == 2 & C.n_ctx > 0) {
            if (sum((C.tg :== g) :* (C.tg_ctx :> 0)) > 0) reason = 4
        }
        if (sum((C.tg :== g) :* C.tg_ref) == 0) reason = 1
        out[i, .] = (0, 0, reason)
    }
    return(out)
}

void _cdid_cat_classify_requests(struct cdid_cat scalar C)
{
    real matrix cl
    cl = _cdid_cat_classify(C, C.rq_g, C.rq_t)
    C.rq_supported = cl[., 1]
    C.rq_ntouch = cl[., 2]
    C.rq_reason = cl[., 3]
}

// -------------------------------------------------- occupancy conversion
// Label rows (cohort, s, t, count) -> index rows; cohort . = never treated.
real matrix cdid_occupancy_from_labels(real matrix occ, real matrix calendar, real matrix cohorts)
{
    real scalar n, i, nt
    real matrix out
    _cdid_check_calendar(calendar, 1)
    nt = _cdid_check_cohorts(cohorts, 0)
    n = rows(occ)
    if (n > 0 & cols(occ) != 4)
        _cdid_fail("cdid_catalog_occupancy_error", "occupancy must have exactly the columns cohort, s, t, count")
    out = J(n, 4, .)
    for (i = 1; i <= n; i++) {
        out[i, 1] = _cdid_cohort_index(cohorts, nt, occ[i, 1])
        if (out[i, 1] == 0)
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy cohort outside the supplied cohort universe")
        out[i, 2] = _cdid_bsearch(calendar, occ[i, 2])
        out[i, 3] = _cdid_bsearch(calendar, occ[i, 3])
        if (out[i, 2] == 0 | out[i, 3] == 0)
            _cdid_fail("cdid_catalog_occupancy_error", "occupancy dates must be exactly stored dates of the full calendar")
        out[i, 4] = occ[i, 4]
    }
    return(out)
}

real matrix cdid_occupancy_to_labels(real matrix occ, real matrix calendar, real matrix cohorts)
{
    real scalar n, nt
    _cdid_check_calendar(calendar, 1)
    nt = _cdid_check_cohorts(cohorts, 0)
    if (nt > 0 & nt != rows(cohorts))
        _cdid_fail("cdid_catalog_internal_error", "never-treated cohort must be the last dictionary entry")
    n = rows(occ)
    if (n == 0) return(J(0, 4, .))
    if (cols(occ) != 4 | !_cdid_is_whole(occ[., (1, 2, 3)]))
        _cdid_fail("cdid_catalog_occupancy_error", "index occupancy must be whole numbers in four columns")
    if (_cdid_any((occ[., 1] :< 1) :| (occ[., 1] :> rows(cohorts)) :| (occ[., 2] :< 1) :| (occ[., 2] :> rows(calendar)) :|
        (occ[., 3] :< 1) :| (occ[., 3] :> rows(calendar))))
        _cdid_fail("cdid_catalog_occupancy_error", "occupancy index outside the dictionaries")
    return((cohorts[occ[., 1]], calendar[occ[., 2]], calendar[occ[., 3]], occ[., 4]))
}

// --------------------------------------------------------------- design
// Spanning forest of R incidence_forest(): CSR adjacency sorted by (vertex,
// neighbour, canonical link rank); BFS seeds are the anchors (in the given
// order) then all vertices; if a component holds two or more anchors, parents
// are recomputed by multi-source BFS from that component's anchors.
struct _cdid_forest {
    real scalar    nc
    real colvector comp, par, pe, ps, deg, roots
}

struct _cdid_forest scalar _cdid_bfs_forest(real scalar nv, real colvector frm,
    real colvector to, real colvector anchors, real colvector eorder)
{
    struct _cdid_forest scalar F
    real scalar m, i, root, nc, head, tail, a, b, kk, na, multi, cc, p
    real colvector start, queue, vis, seeds, ca, oa, asorted, csorted, sc
    real matrix A

    m = rows(frm)
    F.nc = 0
    F.deg = J(nv, 1, 0)
    if (nv == 0) {
        F.comp = J(0, 1, .)
        F.par = J(0, 1, .)
        F.pe = J(0, 1, .)
        F.ps = J(0, 1, .)
        F.roots = J(0, 1, .)
        return(F)
    }
    A = J(0, 5, .)
    if (m > 0) {
        A = (frm, to, (1::m), J(m, 1, 1), eorder) \ (to, frm, (1::m), J(m, 1, -1), eorder)
        A = A[order(A, (1, 2, 5)), .]
        for (i = 1; i <= 2 * m; i++) F.deg[A[i, 1]] = F.deg[A[i, 1]] + 1
    }
    start = runningsum(F.deg) - F.deg
    F.comp = J(nv, 1, 0)
    F.par = J(nv, 1, 0)
    F.pe = J(nv, 1, 0)
    F.ps = J(nv, 1, 0)
    F.roots = J(nv, 1, 0)
    queue = J(nv, 1, 0)
    nc = 0
    seeds = anchors \ (1::nv)
    for (i = 1; i <= rows(seeds); i++) {
        root = seeds[i]
        if (F.comp[root] != 0) continue
        nc++
        F.roots[nc] = root
        F.comp[root] = nc
        head = 1
        tail = 1
        queue[1] = root
        while (head <= tail) {
            a = queue[head]
            head++
            for (kk = start[a] + 1; kk <= start[a] + F.deg[a]; kk++) {
                b = A[kk, 2]
                if (F.comp[b] == 0) {
                    F.comp[b] = nc
                    F.par[b] = a
                    F.pe[b] = A[kk, 3]
                    F.ps[b] = A[kk, 4]
                    tail++
                    queue[tail] = b
                }
            }
        }
    }
    F.nc = nc
    if (nc > 0) F.roots = F.roots[|1 \ nc|]
    else F.roots = J(0, 1, .)
    na = rows(anchors)
    multi = 0
    if (na > 1) {
        ca = F.comp[anchors]
        sc = sort(ca, 1)
        multi = _cdid_any(sc[|2 \ na|] :== sc[|1 \ na - 1|])
    }
    if (multi) {
        F.par = J(nv, 1, 0)
        F.pe = J(nv, 1, 0)
        F.ps = J(nv, 1, 0)
        vis = J(nv, 1, 0)
        oa = order((ca, (1::na)), (1, 2))
        asorted = anchors[oa]
        csorted = ca[oa]
        p = 1
        for (cc = 1; cc <= nc; cc++) {
            tail = 0
            while (p <= na) {
                if (csorted[p] != cc) break
                tail++
                queue[tail] = asorted[p]
                vis[asorted[p]] = 1
                p++
            }
            if (tail == 0) {
                tail = 1
                queue[1] = F.roots[cc]
                vis[F.roots[cc]] = 1
            }
            head = 1
            while (head <= tail) {
                a = queue[head]
                head++
                for (kk = start[a] + 1; kk <= start[a] + F.deg[a]; kk++) {
                    b = A[kk, 2]
                    if (!vis[b]) {
                        vis[b] = 1
                        F.par[b] = a
                        F.pe[b] = A[kk, 3]
                        F.ps[b] = A[kk, 4]
                        tail++
                        queue[tail] = b
                    }
                }
            }
        }
    }
    return(F)
}

// Maintained chained design of a catalogue (R build_comparison_catalog_design
// on catalog_design_handoff links). Links = every catalogue row, in row order.
struct cdid_design {
    real scalar      built, pt_all, pt_ctrl, pt_start, pt_start_index, anticipation
    real colvector   calendar, cohorts
    real scalar      n_links, n_vertices, n_params, n_cohorts, n_components, n_rank, nullity
    real colvector   v_g, v_date, v_kind, v_op, v_fi, v_param, v_anchor, v_degree
    real colvector   v_comp, v_conn, v_tconn, v_par, v_pe, v_ps
    real colvector   l_row, l_g, l_s, l_t, l_ref, l_op, l_sq, l_tq, l_from, l_to
    real colvector   l_start, l_end, l_type, l_comp, l_erank
    real colvector   p_vertex, pt_vertex, pt_status
    real colvector   comp_roots, comp_anchored, unanchored
    real colvector   dc_g, dc_ref
    real matrix      by_cohort
    string rowvector by_cohort_names
    real scalar      unrestricted_n_parameters
    real colvector   free_post_param_indices, restricted_post_param_indices, zero_loading_links
    real matrix      avail
    string colvector op_keys
}

real scalar _cdid_tref(struct cdid_cat scalar C, real scalar g)
{
    return(sum((C.tg :== g) :* C.tg_ref))
}

// 0/1 membership of each q in the sorted key vector s (exact).
real colvector _cdid_member(real colvector s, real colvector q)
{
    real scalar i
    real colvector out
    out = J(rows(q), 1, 0)
    if (rows(s) == 0) return(out)
    for (i = 1; i <= rows(q); i++) out[i] = (_cdid_bsearch(s, q[i]) > 0)
    return(out)
}

string scalar cdid_catalog_functional_key(struct cdid_cat scalar C, real scalar op)
{
    if (op >= . | op < 1 | op > rows(C.op_g) | op != floor(op))
        _cdid_fail("cdid_catalog_input_error", "operator id out of range")
    return("cmpcat:" + C.conventions + ":g" + _cdid_istr(C.op_g[op]) + ":K" + _cdid_istr(C.op_pool[op]))
}

// Public pipeline default (chained_ipw_pipeline.R): qualify every clean
// causal endpoint present in retained rows, including anticipation dates.
// Requested post targets remain the reporting grid; rowless cohorts retain
// availability records without entering the maintained design.
struct cdid_design scalar cdid_catalog_design(struct cdid_cat scalar C, string scalar pt_assumption,
    real scalar pt_start)
{
    real scalar q, i
    real colvector ix
    real matrix lab, declared
    q = rows(C.rq_g)
    lab = J(q, 2, .)
    for (i = 1; i <= q; i++) lab[i, .] = (C.cohorts[C.rq_g[i]], C.calendar[C.rq_t[i]])
    declared = J(0, 2, .)
    ix = _cdid_which((C.r_srole :== 3) :| (C.r_srole :== 5))
    if (rows(ix) > 0) declared = declared \ (C.cohorts[C.r_g[ix, 1], 1], C.calendar[C.r_s[ix, 1], 1])
    ix = _cdid_which((C.r_trole :== 3) :| (C.r_trole :== 5))
    if (rows(ix) > 0) declared = declared \ (C.cohorts[C.r_g[ix, 1], 1], C.calendar[C.r_t[ix, 1], 1])
    if (rows(declared) > 1) declared = uniqrows(sort(declared, (1, 2)))
    return(cdid_catalog_design_ex(C, declared, lab, pt_assumption, pt_start))
}

// declared, targets: q x 2 LABELS (g, t). pt_start: . under "post".
struct cdid_design scalar cdid_catalog_design_ex(struct cdid_cat scalar C, real matrix declared,
    real matrix targets, string scalar pt_assumption, real scalar pt_start)
{
    struct cdid_design scalar D
    real scalar nd, nq, i, g, u
    real colvector dg, du, qg, qt
    real matrix sq, cl

    if (pt_assumption == "post") {
        if (pt_start < .)
            _cdid_fail("cdid_design_pt_start_error", "pt_start must be missing (not supplied) under pt_assumption post")
        D.pt_all = 0
        D.pt_ctrl = 0
        D.pt_start = .
        D.pt_start_index = 0
    }
    else if (pt_assumption == "all") {
        if (pt_start >= .)
            _cdid_fail("cdid_design_pt_start_error", "under pt_assumption all, pt_start must be an explicit stored calendar date")
        D.pt_start_index = _cdid_bsearch(C.calendar, pt_start)
        if (D.pt_start_index == 0)
            _cdid_fail("cdid_design_pt_start_error", "pt_start must belong to the supplied full calendar")
        D.pt_all = 1
        D.pt_ctrl = 0
        D.pt_start = pt_start
    }
    else if (pt_assumption == "ctrl") {
        if (pt_start >= .) pt_start = C.calendar[1]
        D.pt_start_index = _cdid_bsearch(C.calendar, pt_start)
        if (D.pt_start_index == 0)
            _cdid_fail("cdid_design_pt_start_error", "pt_start must belong to the supplied full calendar")
        D.pt_all = 0
        D.pt_ctrl = 1
        D.pt_start = pt_start
    }
    else _cdid_fail("cdid_design_pt_assumption_error", "pt_assumption must be post, ctrl, or all")

    nd = rows(declared)
    if (nd > 0 & cols(declared) != 2)
        _cdid_fail("cdid_catalog_qualification_error", "declared qualifications need the two columns g and t")
    dg = J(nd, 1, 0)
    du = J(nd, 1, 0)
    for (i = 1; i <= nd; i++) {
        if (declared[i, 1] >= . | declared[i, 2] >= .)
            _cdid_fail("cdid_catalog_qualification_error", "declared qualifications need finite g and t")
        g = _cdid_cohort_index(C.cohorts, C.nt_index, declared[i, 1])
        u = _cdid_bsearch(C.calendar, declared[i, 2])
        if (g == 0 | u == 0)
            _cdid_fail("cdid_catalog_qualification_error", "declared qualifications need a referenced catalogue target and a stored calendar date")
        if (_cdid_tref(C, g) == 0)
            _cdid_fail("cdid_catalog_qualification_error", "declared qualifications need a referenced catalogue target and a stored calendar date")
        if (_cdid_cmp_sum(C.calendar[u], C.anticipation, C.cohorts[g]) < 0)
            _cdid_fail("cdid_catalog_qualification_error", "only causal endpoints (u + anticipation >= g) can be declared qualified")
        dg[i] = g
        du[i] = u
    }
    if (nd > 1) {
        sq = sort((dg, du), (1, 2))
        if (_cdid_any((sq[|2, 1 \ nd, 1|] :== sq[|1, 1 \ nd - 1, 1|]) :& (sq[|2, 2 \ nd, 2|] :== sq[|1, 2 \ nd - 1, 2|])))
            _cdid_fail("cdid_catalog_duplicate_error", "duplicate declared qualification")
    }
    nq = rows(targets)
    if (nq > 0 & cols(targets) != 2)
        _cdid_fail("cdid_design_target_error", "design targets need the two columns g and t")
    qg = J(nq, 1, 0)
    qt = J(nq, 1, 0)
    for (i = 1; i <= nq; i++) {
        if (targets[i, 1] >= . | targets[i, 2] >= .)
            _cdid_fail("cdid_design_target_error", "design targets need finite g and t")
        qg[i] = _cdid_cohort_index(C.cohorts, C.nt_index, targets[i, 1])
        qt[i] = _cdid_bsearch(C.calendar, targets[i, 2])
        if (qg[i] == 0 | qt[i] == 0)
            _cdid_fail("cdid_design_target_error", "design targets need a catalogue target cohort and a stored calendar date")
        if (!_cdid_any(C.tg :== qg[i]) | targets[i, 2] < targets[i, 1])
            _cdid_fail("cdid_design_target_error", "invalid cohort or post target date (t >= g required)")
    }
    D.calendar = C.calendar
    D.cohorts = C.cohorts
    D.anticipation = C.anticipation
    D.op_keys = J(rows(C.op_g), 1, "")
    for (i = 1; i <= rows(C.op_g); i++) D.op_keys[i] = cdid_catalog_functional_key(C, i)
    _cdid_design_core(D, C, dg, du, qg, qt)
    cl = _cdid_cat_classify(C, qg, qt)
    _cdid_design_avail(D, qg, qt, cl)
    return(D)
}

real colvector _cdid_vkind(real colvector role, real colvector q)
{
    return((role :== 1) + 2 :* (role :== 2) + ((role :== 3) :| (role :== 4)) :* (3 :+ q) +
        ((role :== 5) :| (role :== 6)) :* (5 :+ q))
}

real colvector _cdid_isnuis(real colvector kind)
{
    return((kind :== 2) :| (kind :== 3) :| (kind :== 5))
}

string rowvector _cdid_by_cohort_names(real scalar pt_all)
{
    string rowvector nm
    nm = ("cohort_index", "ref_index", "n_links", "n_vertices", "n_parameters", "n_pre_parameters",
        "n_post_parameters", "n_components", "rank", "nullity", "cycle_rank", "n_estimable_post",
        "n_isolated_post", "n_disconnected_post", "n_post_nuisance_parameters",
        "n_anticipation_nuisance_parameters", "n_anticipation_parameters", "n_estimable_anticipation")
    if (pt_all) nm = (nm, "pt_start_index", "n_extra_calendar_dates", "n_extra_represented_dates",
        "n_restricted_pre_parameters", "n_zero_loading_links", "no_extra_calendar_window",
        "no_extra_represented_restriction", "no_extra_restrictions")
    return(nm)
}

void _cdid_design_empty(struct cdid_design scalar D)
{
    D.built = 0
    D.n_links = 0
    D.n_vertices = 0
    D.n_params = 0
    D.n_cohorts = 0
    D.n_components = 0
    D.n_rank = 0
    D.nullity = 0
    D.v_g = J(0, 1, .)
    D.v_date = D.v_g
    D.v_kind = D.v_g
    D.v_op = D.v_g
    D.v_fi = D.v_g
    D.v_param = D.v_g
    D.v_anchor = D.v_g
    D.v_degree = D.v_g
    D.v_comp = D.v_g
    D.v_conn = D.v_g
    D.v_tconn = D.v_g
    D.v_par = D.v_g
    D.v_pe = D.v_g
    D.v_ps = D.v_g
    D.l_row = D.v_g
    D.l_g = D.v_g
    D.l_s = D.v_g
    D.l_t = D.v_g
    D.l_ref = D.v_g
    D.l_op = D.v_g
    D.l_sq = D.v_g
    D.l_tq = D.v_g
    D.l_from = D.v_g
    D.l_to = D.v_g
    D.l_start = D.v_g
    D.l_end = D.v_g
    D.l_type = D.v_g
    D.l_comp = D.v_g
    D.l_erank = D.v_g
    D.p_vertex = D.v_g
    D.pt_vertex = D.v_g
    D.pt_status = D.v_g
    D.comp_roots = D.v_g
    D.comp_anchored = D.v_g
    D.unanchored = D.v_g
    D.dc_g = D.v_g
    D.dc_ref = D.v_g
    D.by_cohort_names = _cdid_by_cohort_names(D.pt_all)
    D.by_cohort = J(0, cols(D.by_cohort_names), .)
    D.unrestricted_n_parameters = 0
    D.free_post_param_indices = D.v_g
    D.restricted_post_param_indices = D.v_g
    D.zero_loading_links = D.v_g
}

void _cdid_design_core(struct cdid_design scalar D, struct cdid_cat scalar C,
    real colvector dg, real colvector du, real colvector qg, real colvector qt)
{
    struct _cdid_forest scalar F, TF
    struct cdid_catalog_opts scalar mo
    real scalar m, nT, ndc, ngr, nop, i, k, p10, maxd, nE, e, nv, newv, ntk, h, g, nl, nul, np, ext, rep, rst
    real colvector dkey, dcg, keepq, ndig, pad, fo, fi, ks, kt, fvs, fvt, o, vid, from, to, w
    real colvector anchors, eo, erank, tk, tid, tanch, nonref, posnr, zero, compg, mv, mp, rmask
    real matrix grid, V, BC

    _cdid_design_empty(D)
    m = rows(C.r_tp)
    if (m == 0) return
    // Numeric-memory preflight with the catalogue's policy: vertex keys, CSR
    // adjacency of both forests, sort keys and per-vertex/per-link outputs are
    // fewer than 128 doubles per link, target or cohort.
    mo = cdid_catalog_default_opts()
    if (C.mem_limit_bytes < .) mo.mem_limit_bytes = C.mem_limit_bytes
    _cdid_mem(mo, 128 * (m + rows(qg) + rows(C.tg)), "design links, vertices and forests")
    D.built = 1
    nT = C.n_dates
    D.n_links = m
    D.l_row = (1::m)
    D.l_g = C.r_g
    D.l_s = C.r_s
    D.l_t = C.r_t
    D.l_ref = C.r_ref
    D.l_op = C.r_op
    // Qualification (R catalog_design_handoff): the reference is always
    // qualified; a conditional causal endpoint is qualified iff declared.
    dkey = J(0, 1, .)
    if (rows(dg) > 0) dkey = sort((dg :- 1) :* nT + du, 1)
    D.l_sq = (C.r_srole :== 1) :| (((C.r_srole :== 3) :| (C.r_srole :== 5)) :& _cdid_member(dkey, (C.r_g :- 1) :* nT + C.r_s))
    D.l_tq = (C.r_trole :== 1) :| (((C.r_trole :== 3) :| (C.r_trole :== 5)) :& _cdid_member(dkey, (C.r_g :- 1) :* nT + C.r_t))

    dcg = uniqrows(C.r_g)
    ndc = rows(dcg)
    D.n_cohorts = ndc
    D.dc_g = dcg
    D.dc_ref = J(ndc, 1, 0)
    for (i = 1; i <= ndc; i++) D.dc_ref[i] = _cdid_tref(C, dcg[i])
    grid = J(0, 2, .)
    if (rows(qg) > 0) {
        keepq = _cdid_member(dcg, qg)
        if (_cdid_any(keepq)) grid = uniqrows(select((qg, qt), keepq))
    }
    ngr = rows(grid)

    // Functional-key rank inside a cohort. R sorts "cmpcat:<conv>:g<gi>:K<pool>"
    // bytewise; within one cohort that is the decimal-string order of the pool
    // id (K10 < K2): right-pad the digits, then fewer digits first.
    nop = rows(C.op_g)
    ndig = J(nop, 1, 1)
    for (i = 1; i <= nop; i++) {
        p10 = 10
        while (C.op_pool[i] >= p10) {
            ndig[i] = ndig[i] + 1
            p10 = p10 * 10
        }
    }
    maxd = max(ndig)
    pad = C.op_pool
    for (i = 1; i <= nop; i++) {
        for (k = ndig[i]; k < maxd; k++) pad[i] = pad[i] * 10
    }
    fo = order((C.op_g, pad, ndig, (1::nop)), (1, 2, 3, 4))
    fi = J(nop, 1, 0)
    fi[fo] = (1::nop)

    // Vertex keys (cohort, kind, date, functional rank) sorted as R does by
    // (cohort, role kind, time, functional rank). Kinds 2, 3, 5 are private
    // per operator; 1 (reference), 4 and 6 (qualified causal) are shared.
    ks = _cdid_vkind(C.r_srole, D.l_sq)
    kt = _cdid_vkind(C.r_trole, D.l_tq)
    fvs = _cdid_isnuis(ks) :* fi[C.r_op, 1]
    fvt = _cdid_isnuis(kt) :* fi[C.r_op, 1]
    // Under shared control pre-trends, pre-reference differentials from pt_start
    // onward are shared across control pools.
    if (D.pt_ctrl) {
        fvs = fvs :* !((ks :== 2) :& (C.calendar[C.r_s] :>= D.pt_start))
        fvt = fvt :* !((kt :== 2) :& (C.calendar[C.r_t] :>= D.pt_start))
    }
    V = (dcg, J(ndc, 1, 1), D.dc_ref, J(ndc, 1, 0)) \ (C.r_g, ks, C.r_s, fvs) \ (C.r_g, kt, C.r_t, fvt)
    if (ngr > 0) V = V \ (grid[., 1], J(ngr, 1, 6), grid[., 2], J(ngr, 1, 0))
    nE = rows(V)
    o = order((V, (1::nE)), (1, 2, 3, 4, 5))
    vid = J(nE, 1, 0)
    D.v_g = J(nE, 1, 0)
    D.v_kind = J(nE, 1, 0)
    D.v_date = J(nE, 1, 0)
    D.v_fi = J(nE, 1, 0)
    nv = 0
    for (i = 1; i <= nE; i++) {
        e = o[i]
        newv = 1
        if (i > 1) newv = (V[e, .] != V[o[i - 1], .])
        if (newv) {
            nv++
            D.v_g[nv] = V[e, 1]
            D.v_kind[nv] = V[e, 2]
            D.v_date[nv] = V[e, 3]
            D.v_fi[nv] = V[e, 4]
        }
        vid[e] = nv
    }
    D.v_g = D.v_g[|1 \ nv|]
    D.v_kind = D.v_kind[|1 \ nv|]
    D.v_date = D.v_date[|1 \ nv|]
    D.v_fi = D.v_fi[|1 \ nv|]
    D.n_vertices = nv
    D.v_op = J(nv, 1, 0)
    w = _cdid_which(D.v_fi :> 0)
    if (rows(w) > 0) D.v_op[w] = fo[D.v_fi[w], 1]
    D.v_anchor = (D.v_kind :== 1) :| (D.pt_all :* (D.v_kind :== 2) :* (C.calendar[D.v_date] :>= D.pt_start))
    D.v_param = runningsum(!D.v_anchor) :* (!D.v_anchor)
    D.p_vertex = _cdid_which(!D.v_anchor)
    D.n_params = rows(D.p_vertex)
    from = vid[|ndc + 1 \ ndc + m|]
    to = vid[|ndc + m + 1 \ ndc + 2 * m|]
    D.l_from = from
    D.l_to = to
    D.l_start = D.v_param[from]
    D.l_end = D.v_param[to]

    eo = order((C.r_g, C.r_s, C.r_t, fi[C.r_op, 1], (1::m)), (1, 2, 3, 4, 5))
    erank = J(m, 1, 0)
    erank[eo] = (1::m)
    D.l_erank = erank
    anchors = _cdid_which(D.v_anchor)
    F = _cdid_bfs_forest(nv, from, to, anchors, erank)
    D.n_components = F.nc
    D.v_degree = F.deg
    D.v_comp = F.comp
    D.v_par = F.par
    D.v_pe = F.pe
    D.v_ps = F.ps
    D.comp_roots = F.roots
    D.comp_anchored = J(F.nc, 1, 0)
    D.comp_anchored[F.comp[anchors], 1] = J(rows(anchors), 1, 1)
    D.v_conn = D.comp_anchored[F.comp, 1]
    D.unanchored = _cdid_which(!D.comp_anchored)
    D.nullity = rows(D.unanchored)
    D.n_rank = D.n_params - D.nullity
    D.l_comp = F.comp[from]

    // Naive (cohort, date) time graph: exposes false collapsed connectivity.
    tk = order((D.v_g, D.v_date, (1::nv)), (1, 2, 3))
    tid = J(nv, 1, 0)
    ntk = 0
    for (i = 1; i <= nv; i++) {
        e = tk[i]
        newv = 1
        if (i > 1) newv = (D.v_g[e] != D.v_g[tk[i - 1]] | D.v_date[e] != D.v_date[tk[i - 1]])
        if (newv) ntk++
        tid[e] = ntk
    }
    TF = _cdid_bfs_forest(ntk, tid[from], tid[to], uniqrows(tid[anchors]), erank)
    tanch = J(TF.nc, 1, 0)
    tanch[TF.comp[tid[anchors]], 1] = J(rows(anchors), 1, 1)
    D.v_tconn = tanch[TF.comp[tid], 1]

    D.l_type = J(m, 1, 5)
    for (i = 1; i <= m; i++) {
        if (C.r_t[i] <= C.r_ref[i]) D.l_type[i] = 1
        else if (C.r_s[i] < C.r_ref[i]) D.l_type[i] = 2
        else if (C.calendar[C.r_t[i]] < C.cohorts[C.r_g[i]]) D.l_type[i] = 3
        else if (C.calendar[C.r_s[i]] < C.cohorts[C.r_g[i]] & C.r_s[i] != C.r_ref[i]) D.l_type[i] = 4
    }
    D.pt_vertex = _cdid_which(D.v_kind :== 6)
    D.pt_status = J(rows(D.pt_vertex), 1, 0)
    for (i = 1; i <= rows(D.pt_vertex); i++) {
        k = D.pt_vertex[i]
        D.pt_status[i] = (D.v_conn[k] ? 1 : (D.v_degree[k] == 0 ? 2 : 3))
    }

    zero = (D.l_start :== 0) :& (D.l_end :== 0)
    D.zero_loading_links = _cdid_which(zero)
    compg = D.v_g[F.roots]
    BC = J(ndc, cols(D.by_cohort_names), 0)
    for (h = 1; h <= ndc; h++) {
        g = dcg[h]
        mv = (D.v_g :== g)
        mp = mv :& (!D.v_anchor)
        nl = sum(C.r_g :== g)
        nul = sum((compg :== g) :& (!D.comp_anchored))
        np = sum(mp)
        BC[h, (1..18)] = (g, D.dc_ref[h], nl, sum(mv), np, sum(mp :& (D.v_kind :== 2)),
            sum(mp :& (D.v_kind :== 6)), sum(compg :== g), np - nul, nul, nl - (np - nul),
            sum(mp :& (D.v_kind :== 6) :& D.v_conn), sum(mp :& (D.v_kind :== 6) :& (D.v_degree :== 0)),
            sum(mp :& (D.v_kind :== 6) :& (!D.v_conn) :& (D.v_degree :> 0)), sum(mp :& (D.v_kind :== 5)),
            sum(mp :& (D.v_kind :== 3)), sum(mp :& (D.v_kind :== 4)), sum(mp :& (D.v_kind :== 4) :& D.v_conn))
        if (D.pt_all) {
            ext = sum((C.calendar :>= D.pt_start) :& (C.calendar :< C.calendar[D.dc_ref[h]]))
            rmask = mv :& (D.v_kind :== 2) :& (C.calendar[D.v_date] :>= D.pt_start)
            rep = 0
            if (_cdid_any(rmask)) rep = rows(uniqrows(select(D.v_date, rmask)))
            rst = sum(mv :& (D.v_kind :== 2) :& D.v_anchor)
            BC[h, (19..26)] = (D.pt_start_index, ext, rep, rst, sum(zero :& (C.r_g :== g)), (ext == 0), (rst == 0), (rst == 0))
        }
    }
    D.by_cohort = BC
    nonref = _cdid_which(D.v_kind :!= 1)
    D.unrestricted_n_parameters = rows(nonref)
    posnr = J(nv, 1, 0)
    if (rows(nonref) > 0) posnr[nonref] = (1::rows(nonref))
    D.free_post_param_indices = posnr[D.p_vertex]
    D.restricted_post_param_indices = posnr[_cdid_which((D.v_kind :!= 1) :& D.v_anchor)]
}

// Per requested target: catalogue support then maintained-design status.
// Columns: g index, t index, catalogue supported, rows touching, catalogue
// reason, in design, parameter index, design status (0 = not in design),
// available (estimable level in the maintained design).
void _cdid_design_avail(struct cdid_design scalar D, real colvector qg, real colvector qt, real matrix cl)
{
    real scalar q, i, v, st
    real colvector hit
    q = rows(qg)
    D.avail = J(q, 9, 0)
    for (i = 1; i <= q; i++) {
        D.avail[i, (1..5)] = (qg[i], qt[i], cl[i, 1], cl[i, 2], cl[i, 3])
        if (D.built) {
            hit = _cdid_which((D.v_g :== qg[i]) :& (D.v_date :== qt[i]) :& (D.v_kind :== 6))
            if (rows(hit) == 1) {
                v = hit[1]
                st = (D.v_conn[v] ? 1 : (D.v_degree[v] == 0 ? 2 : 3))
                D.avail[i, (6..9)] = (1, D.v_param[v], st, (st == 1))
            }
        }
    }
}

// ------------------------------------------------ incidence operations
// W x (links x k): x[end] - x[start] with anchored endpoints grounded at 0.
// No dense links-by-parameters matrix is formed.
real matrix cdid_design_mult(struct cdid_design scalar D, real matrix x)
{
    real matrix xx, out
    if (rows(x) != D.n_params | hasmissing(x))
        _cdid_fail("cdid_design_input_error", "invalid incidence product input (rows must equal the parameter count; no missing)")
    xx = J(1, cols(x), 0) \ x
    out = xx[D.l_end :+ 1, .] - xx[D.l_start :+ 1, .]
    if (hasmissing(out)) _cdid_fail("cdid_design_product_error", "incidence product overflow")
    return(out)
}

// W' y (parameters x k), accumulated in link order per side as R rowsum():
// minus the start-side sums, then plus the end-side sums.
real matrix cdid_design_tmult(struct cdid_design scalar D, real matrix y)
{
    real scalar l
    real matrix st, en, out
    if (rows(y) != D.n_links | hasmissing(y))
        _cdid_fail("cdid_design_input_error", "invalid incidence product input (rows must equal the link count; no missing)")
    st = J(D.n_params, cols(y), 0)
    en = J(D.n_params, cols(y), 0)
    for (l = 1; l <= D.n_links; l++) {
        if (D.l_start[l] > 0) st[D.l_start[l], .] = st[D.l_start[l], .] + y[l, .]
        if (D.l_end[l] > 0) en[D.l_end[l], .] = en[D.l_end[l], .] + y[l, .]
    }
    out = (-st) + en
    if (hasmissing(out)) _cdid_fail("cdid_design_product_error", "incidence product overflow")
    return(out)
}

void _cdid_entry_bound(real scalar need, real scalar max_entries)
{
    if (max_entries >= . | max_entries < 0 | max_entries != floor(max_entries))
        _cdid_fail("cdid_design_bound_error", "supply an explicit finite whole entry bound")
    if (need > max_entries)
        _cdid_fail("cdid_design_bound_error", "dense reference of " + _cdid_istr(need) +
            " entries exceeds the explicit entry bound")
}

// Test/debug only: dense W (links x parameters), refused beyond max_entries.
real matrix cdid_design_materialize(struct cdid_design scalar D, real scalar max_entries)
{
    real scalar l
    real matrix W
    _cdid_entry_bound(D.n_links * D.n_params, max_entries)
    W = J(D.n_links, D.n_params, 0)
    for (l = 1; l <= D.n_links; l++) {
        if (D.l_start[l] > 0) W[l, D.l_start[l]] = -1
        if (D.l_end[l] > 0) W[l, D.l_end[l]] = 1
    }
    return(W)
}

// Null-space basis: one indicator column per unanchored component.
real matrix cdid_design_nullspace(struct cdid_design scalar D, real scalar max_entries)
{
    real scalar cc
    real matrix B
    _cdid_entry_bound(D.n_params * D.nullity, max_entries)
    B = J(D.n_params, D.nullity, 0)
    for (cc = 1; cc <= D.nullity; cc++) {
        B[., cc] = select(D.v_comp, !D.v_anchor) :== D.unanchored[cc]
    }
    return(B)
}

// Signed root path of a vertex (link index, sign), sorted by link index.
real matrix cdid_design_root_path(struct cdid_design scalar D, real scalar v)
{
    real scalar k, cur
    real matrix out
    if (v >= . | v < 1 | v > D.n_vertices | v != floor(v)) _cdid_fail("cdid_design_input_error", "vertex index out of range")
    k = 0
    cur = v
    while (D.v_par[cur] != 0) {
        k++
        cur = D.v_par[cur]
    }
    out = J(k, 2, 0)
    cur = v
    for (k = 1; k <= rows(out); k++) {
        out[k, .] = (D.v_pe[cur], D.v_ps[cur])
        cur = D.v_par[cur]
    }
    if (rows(out) > 1) out = sort(out, (1, 2))
    return(out)
}

// Level witness w (link, sign) with W'w = e_param, or J(0,2,.) if not estimable.
real matrix cdid_design_witness(struct cdid_design scalar D, real scalar param)
{
    if (param >= . | param < 1 | param > D.n_params | param != floor(param))
        _cdid_fail("cdid_design_input_error", "unknown or noninteger parameter index")
    if (!D.v_conn[D.p_vertex[param]]) return(J(0, 2, .))
    return(cdid_design_root_path(D, D.p_vertex[param]))
}

// theta_b - theta_a is estimable iff both lie in one component or both
// levels are estimable (R check_contrast_estimable).
real scalar cdid_design_contrast_estimable(struct cdid_design scalar D, real scalar pa, real scalar pb)
{
    real scalar va, vb
    if (pa >= . | pb >= . | pa < 1 | pb < 1 | pa > D.n_params | pb > D.n_params | pa != floor(pa) | pb != floor(pb))
        _cdid_fail("cdid_design_input_error", "unknown or noninteger parameter index")
    va = D.p_vertex[pa]
    vb = D.p_vertex[pb]
    return(D.v_comp[va] == D.v_comp[vb] | (D.v_conn[va] & D.v_conn[vb]))
}

// Contrast witness (link, sign) with W'w = e_b - e_a; J(0,2,.) when the
// contrast is not estimable or trivially zero (use the function above).
real matrix cdid_design_contrast_witness(struct cdid_design scalar D, real scalar pa, real scalar pb)
{
    real scalar va, vb, i, j
    real matrix A, B, out
    if (!cdid_design_contrast_estimable(D, pa, pb)) return(J(0, 2, .))
    va = D.p_vertex[pa]
    vb = D.p_vertex[pb]
    A = cdid_design_root_path(D, va)
    B = cdid_design_root_path(D, vb)
    out = J(0, 2, .)
    i = 1
    j = 1
    while (i <= rows(A) | j <= rows(B)) {
        if (j > rows(B)) {
            out = out \ (A[i, 1], -A[i, 2])
            i++
        }
        else if (i > rows(A)) {
            out = out \ B[j, .]
            j++
        }
        else if (A[i, 1] < B[j, 1]) {
            out = out \ (A[i, 1], -A[i, 2])
            i++
        }
        else if (B[j, 1] < A[i, 1]) {
            out = out \ B[j, .]
            j++
        }
        else {
            if (B[j, 2] - A[i, 2] != 0) out = out \ (A[i, 1], B[j, 2] - A[i, 2])
            i++
            j++
        }
    }
    return(out)
}

// ------------------------------------------------------- names/accessors
string scalar _cdid_pick(string rowvector nm, real scalar code)
{
    if (code >= . | code < 1 | code > cols(nm) | code != floor(code)) return("")
    return(nm[code])
}

string scalar cdid_catalog_role_name(real scalar code)
{
    return(_cdid_pick(("reference", "pre_private", "anticipation_conditional", "anticipation_private",
        "post_conditional", "post_private"), code))
}

string scalar cdid_coordinate_kind_name(real scalar code)
{
    return(_cdid_pick(("reference", "shared_conditional", "private"), code))
}

string scalar cdid_target_pair_status_name(real scalar code)
{
    return(_cdid_pick(("retained", "empty_pool", "no_admissible_reference"), code))
}

string scalar cdid_target_status_name(real scalar code)
{
    return(_cdid_pick(("reference_available", "no_admissible_reference", "reference_available_eligibility_restricted"), code))
}

string scalar cdid_request_reason_name(real scalar code)
{
    if (code == 0) return("supported")
    return(_cdid_pick(("no_admissible_reference", "no_target_pair_occupancy", "no_nonempty_control_pool",
        "no_target_pair_occupancy_in_eligibility_context"), code))
}

string scalar cdid_design_role_name(real scalar code)
{
    return(_cdid_pick(("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post"), code))
}

string scalar cdid_link_type_name(real scalar code)
{
    return(_cdid_pick(("pre_only", "straddling", "anticipation_only", "anticipation_to_post", "post_only"), code))
}

string scalar cdid_post_status_name(real scalar code)
{
    if (code == 0) return("not_in_design")
    return(_cdid_pick(("estimable", "isolated_unsupported", "disconnected_component"), code))
}

// Catalogue resource counts under the R resources names.
real scalar cdid_catalog_count(struct cdid_cat scalar C, string scalar key)
{
    if (key == "occupancy_rows") return(C.n_occ)
    if (key == "requested_rows") return(rows(C.rq_g))
    if (key == "positive_cells") return(C.n_pos)
    if (key == "target_pairs") return(rows(C.tp_j))
    if (key == "candidates") return(C.n_candidates)
    if (key == "nonempty_candidates") return(C.n_nonempty)
    if (key == "rows") return(rows(C.r_tp))
    if (key == "pool_entries_preintern") return(C.pool_entries_preintern)
    if (key == "pool_entries_interned") return(rows(C.pool_entries))
    if (key == "pools") return(C.n_pools)
    if (key == "operators") return(rows(C.op_g))
    if (key == "coordinates") return(rows(C.co_g))
    if (key == "cutoff_comparisons") return(C.n_cutoff_cmp)
    if (key == "work") return(C.work)
    _cdid_fail("cdid_catalog_input_error", "unknown catalogue count " + key)
    return(.)
}

real scalar _cdid_bc_sum(struct cdid_design scalar D, string scalar name)
{
    real scalar j
    for (j = 1; j <= cols(D.by_cohort_names); j++) {
        if (D.by_cohort_names[j] == name) return(sum(D.by_cohort[., j]))
    }
    _cdid_fail("cdid_design_input_error", "by-cohort column " + name + " is not defined for this design")
    return(.)
}

// Design diagnostics under the R rank_diagnostics / pt_all_spec names.
real scalar cdid_design_diag_value(struct cdid_design scalar D, string scalar key)
{
    real colvector k6
    k6 = (D.v_kind :== 6) :& (!D.v_anchor)
    if (key == "n_links") return(D.n_links)
    if (key == "n_vertices") return(D.n_vertices)
    if (key == "n_parameters" | key == "restricted_n_parameters") return(D.n_params)
    if (key == "n_cohorts") return(D.n_cohorts)
    if (key == "n_components") return(D.n_components)
    if (key == "rank") return(D.n_rank)
    if (key == "nullity") return(D.nullity)
    if (key == "cycle_rank") return(D.n_links - D.n_rank)
    if (key == "has_full_column_rank") return(D.nullity == 0)
    if (key == "n_post_targets") return(rows(D.pt_vertex))
    if (key == "n_estimable_post_targets") return(sum(D.pt_status :== 1))
    if (key == "n_isolated_post_targets") return(sum(D.pt_status :== 2))
    if (key == "n_disconnected_post_targets") return(sum(D.pt_status :== 3))
    if (key == "n_pre_parameters") return(sum((D.v_kind :== 2) :& (!D.v_anchor)))
    if (key == "n_post_nuisance_parameters") return(sum((D.v_kind :== 5) :& (!D.v_anchor)))
    if (key == "n_post_parameters") return(sum(k6))
    if (key == "n_anticipation_nuisance_parameters") return(sum((D.v_kind :== 3) :& (!D.v_anchor)))
    if (key == "n_anticipation_parameters") return(sum((D.v_kind :== 4) :& (!D.v_anchor)))
    if (key == "n_estimable_anticipation") return(sum((D.v_kind :== 4) :& (!D.v_anchor) :& D.v_conn))
    if (key == "unrestricted_n_parameters") return(D.unrestricted_n_parameters)
    if (key == "pt_start_index") return(D.pt_start_index)
    if (!D.pt_all) _cdid_fail("cdid_design_input_error", "diagnostic " + key + " is defined only under pt_assumption all")
    if (key == "n_extra_calendar_dates") return(_cdid_bc_sum(D, "n_extra_calendar_dates"))
    if (key == "n_extra_represented_dates") return(_cdid_bc_sum(D, "n_extra_represented_dates"))
    if (key == "n_restricted_pre_parameters") return(_cdid_bc_sum(D, "n_restricted_pre_parameters"))
    if (key == "n_zero_loading_links") return(rows(D.zero_loading_links))
    if (key == "no_extra_calendar_window") return(_cdid_bc_sum(D, "n_extra_calendar_dates") == 0)
    if (key == "no_extra_represented_restriction" | key == "no_extra_restrictions" | key == "restriction_is_identity")
        return(_cdid_bc_sum(D, "n_restricted_pre_parameters") == 0)
    _cdid_fail("cdid_design_input_error", "unknown design diagnostic " + key)
    return(.)
}


/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
transmorphic cdid_mrca_catalog(transmorphic a1, transmorphic a2, transmorphic a3, transmorphic a4, transmorphic a5, transmorphic a6)
{
    return(cdid_mref_catalog(a1, a2, a3, a4, a5, a6))
}

end
