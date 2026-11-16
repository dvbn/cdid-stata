*! Private baseline conditional-pretrend process piece producer; StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (commit ):
*! R/chained_conditional_extract.R (cpt_piece, cpt_two_arm_piece,
*! cpt_extract_legacy_ipw, cpt_two_arm_rows, cpt_extract_two_arm, cpt_extract_or)
*! and R/chained_conditional_pretest.R (cpt_declared_cells).
*! Supports all 3 baseline methods (ipw, reg, dr) for both static
*! (cdid_baseline_provider) and dated (cdid_db_provider) owners.
*! Reuses the complete post + both-auxiliary union fit without refitting logits
*! or WLS by threshold or chunk.
version 19.5
mata:

struct cdid_cpt_base_ps_factor {
    real scalar p, m, sum_a, sum_b, D_T, D_M
    real matrix H_chol
}

struct cdid_cpt_base_cell_meta {
    string scalar id, status, reason, route
    real scalar g, t, from, s_date, t2_date, s_index, t2_index
    real scalar q, consumer, recipe, design_index, sign, supported, is_legacy_ipw
    real scalar estimate, estimate_scalar, n_chunks, chunk_rows
    struct cdid_cpt_spec scalar spec
}

struct cdid_cpt_piece {
    real colvector rows
    real matrix nu, rho, P, xi, blam
    struct cdid_cpt_spec scalar spec
}

struct cdid_cpt_baseline {
    string scalar kind, method, base_period, status, reason
    real scalar n_F, n_A, n_cells, n_unavailable, chunk_rows
    real scalar workspace_entries, workspace_limit, work_ops, work_limit
    pointer(struct cdid_baseline_provider scalar) scalar static_provider
    pointer(struct cdid_db_provider scalar) scalar dated_provider
    struct cdid_cpt_base_ps_factor colvector ps_factor
    struct cdid_cpt_base_cell_meta colvector cells
}

void _cdid_cpt_base_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

real colvector _cdid_cpt_base_ix(real matrix pred)
{
    return(vec(selectindex(pred)))
}

/* Solve (H_chol * H_chol') * X' = B' for X (n x p) using lower-triangular H_chol. */
real matrix _cdid_cpt_base_chol_solve(real matrix H_chol, real matrix B)
{
    real matrix Z, X
    real scalar rk, p
    p = rows(H_chol)
    if (p == 0 | rows(B) == 0) return(J(rows(B), p, 0))
    Z = solvelower(H_chol, B', rk=., 0)
    if (rk != p | missing(Z))
        _cdid_cpt_base_fail("numeric", "Propensity Hessian triangular solve failed.", 504)
    X = solveupper(H_chol', Z, rk=., 0)'
    if (rk != p | missing(X))
        _cdid_cpt_base_fail("numeric", "Propensity Hessian upper triangular solve failed.", 504)
    return(X)
}

/* Exact R cpt_piece(rows, nu, rho, lam, pim, sD, xi) for 1-term or 2-term baseline cells. */
struct cdid_cpt_piece scalar cdid_cpt_make_piece(real colvector rows,
    real matrix nu, real matrix rho,
    real matrix lam1, real matrix lam2,
    real matrix pim1, real matrix pim2,
    real colvector sD, real matrix xi,
    struct cdid_cpt_spec scalar spec)
{
    struct cdid_cpt_piece scalar pc
    real scalar n, R, Dl
    n = rows(rows); R = spec.R; Dl = spec.Dl
    if (R < 1 | R > 2 | cols(rows) != 1 | rows(nu) != n | cols(nu) != R | rows(rho) != n | cols(rho) != R |
        rows(xi) != n | cols(xi) != Dl | rows(sD) != R | cols(sD) != 1)
        _cdid_cpt_base_fail("piece", "Invalid baseline process piece dimensions.", 498)
    if ((length(lam1) & (rows(lam1)!=n | cols(lam1)!=Dl)) |
        (length(lam2) & (rows(lam2)!=n | cols(lam2)!=Dl)) |
        (length(pim1) & (rows(pim1)!=n | cols(pim1)!=Dl)) |
        (length(pim2) & (rows(pim2)!=n | cols(pim2)!=Dl)))
        _cdid_cpt_base_fail("piece", "A supplied loading matrix has invalid dimensions.", 498)
    pc.rows = rows
    pc.nu = nu
    pc.rho = rho
    pc.xi = xi
    pc.spec = spec
    pc.P = J(n, Dl, 0)
    pc.blam = J(Dl, R, 0)
    if (Dl == 0 | n == 0) return(pc)
    if (rows(lam1) == n & cols(lam1) == Dl) {
        pc.P = pc.P + ((sD[1] :* nu[, 1] :* rho[, 1]) :* lam1)
        pc.blam[, 1] = colsum(nu[, 1] :* lam1)'
    }
    if (rows(pim1) == n & cols(pim1) == Dl) {
        pc.P = pc.P + ((sD[1] :* nu[, 1]) :* pim1)
    }
    if (R >= 2) {
        if (rows(lam2) == n & cols(lam2) == Dl) {
            pc.P = pc.P + ((sD[2] :* nu[, 2] :* rho[, 2]) :* lam2)
            pc.blam[, 2] = colsum(nu[, 2] :* lam2)'
        }
        if (rows(pim2) == n & cols(pim2) == Dl) {
            pc.P = pc.P + ((sD[2] :* nu[, 2]) :* pim2)
        }
    }
    if (missing(pc.nu) | missing(pc.rho) | missing(pc.P) | missing(pc.xi) | missing(pc.blam))
        _cdid_cpt_base_fail("numeric", "Nonfinite baseline process piece for cell " + spec.id + ".", 504)
    return(pc)
}

/* Exact R cpt_two_arm_piece(rows, rr, D, fw, X, odds, active, rho, lam_up, xi,
   D_T, D_M, n_up, p_trt, d_or, or_icpt). */
struct cdid_cpt_piece scalar cdid_cpt_two_arm_piece(real colvector rows,
    real colvector rr, real colvector D, real colvector fw, real matrix X,
    real colvector odds, real colvector active, real colvector rho,
    real matrix lam_up, real matrix xi, real scalar D_T, real scalar D_M,
    real scalar n_up, real scalar p_trt, real matrix d_or, real scalar or_icpt,
    struct cdid_cpt_spec scalar spec)
{
    real scalar n, Dl, n_obs, off, p_or
    real matrix nu, rh, lamT, lamM, pim
    n = rows(rows); Dl = spec.Dl; n_obs = rows(rr)
    nu = J(n, 2, 0); rh = J(n, 2, 0)
    lamT = J(n, Dl, 0); lamM = J(n, Dl, 0)
    if (n_obs) {
        nu[rr, 1] = fw :* D
        nu[rr, 2] = fw :* (1 :- D) :* odds
        rh[rr, 1] = rho
        rh[rr, 2] = rho
        if (n_up) {
            lamT[rr, (1..n_up)] = lam_up
            lamM[rr, (1..n_up)] = lam_up
        }
        if (p_trt) {
            lamM[rr, n_up :+ (1..p_trt)] = active :* X
        }
    }
    pim = J(0, 0, 0)
    p_or = cols(d_or)
    if (or_icpt | rows(d_or) > 0 | p_or > 0) {
        if (Dl > n_up + p_trt) {
            pim = J(n, Dl, 0)
            off = n_up + p_trt
            if (or_icpt) {
                if (n_obs) pim[rr, off + 1] = J(n_obs, 1, -1)
                off = off + 1
            }
            if (p_or & n_obs) {
                pim[rr, off :+ (1..p_or)] = -d_or
            }
        }
    }
    return(cdid_cpt_make_piece(rows, nu, rh, lamT, lamM, pim, pim, (1 / D_T \ -1 / D_M), xi, spec))
}

/* Coordinate label builders matching R chained_conditional_extract.R lines 130, 250-251, 289. */
string colvector _cdid_cpt_base_coords_two_arm(real scalar is_legacy_ipw, real scalar is_dr,
    real scalar m, real scalar p_trt, real scalar or_icpt, real scalar reg_block, real scalar p_or,
    | real scalar separate_in)
{
    string colvector out
    string scalar prefix, b_tag
    real scalar j, at, total, sep
    sep = (args() >= 8 & separate_in != 0 & separate_in < .)
    total = m + p_trt + (is_dr & or_icpt ? 1 : 0) + (is_dr ? p_or : 0)
    out = J(total, 1, "")
    prefix = (is_legacy_ipw ? "rate:" : "up:")
    at = 0
    for (j = 1; j <= m; j++) {
        at++
        out[at] = prefix + strtrim(strofreal(j, "%21.0f"))
    }
    for (j = 1; j <= p_trt; j++) {
        at++
        out[at] = "ps:" + strtrim(strofreal(j, "%21.0f"))
    }
    if (is_dr) {
        b_tag = strtrim(strofreal(reg_block, "%21.0f"))
        if (or_icpt) {
            at++
            out[at] = "or_alpha:" + b_tag
        }
        for (j = 1; j <= p_or; j++) {
            at++
            out[at] = (sep ? "or_beta:" + b_tag + ":" + strtrim(strofreal(j, "%21.0f")) :
                "or_beta:" + strtrim(strofreal(j, "%21.0f")))
        }
    }
    return(out)
}

string colvector _cdid_cpt_base_coords_or(real scalar support_id, real scalar or_icpt,
    real scalar reg_block, real scalar p_or, | real scalar separate_in)
{
    string colvector out
    string scalar b_tag
    real scalar j, at, total, sep
    sep = (args() >= 5 & separate_in != 0 & separate_in < .)
    total = 1 + (or_icpt ? 1 : 0) + p_or
    out = J(total, 1, "")
    out[1] = "up:" + strtrim(strofreal(support_id, "%21.0f"))
    at = 1
    b_tag = strtrim(strofreal(reg_block, "%21.0f"))
    if (or_icpt) {
        at++
        out[at] = "or_alpha:" + b_tag
    }
    for (j = 1; j <= p_or; j++) {
        at++
        out[at] = (sep ? "or_beta:" + b_tag + ":" + strtrim(strofreal(j, "%21.0f")) :
            "or_beta:" + strtrim(strofreal(j, "%21.0f")))
    }
    return(out)
}

/* Recover bounded p x p PS Cholesky factor H_chol and exact sum_a, sum_b from the SAME
   fitted ps.beta and fixed risk design once per unique PS model in static provider. */
void _cdid_cpt_base_ps_static(struct cdid_cpt_baseline scalar B,
    struct cdid_baseline_provider scalar f)
{
    struct cdid_base_ps scalar ps
    real scalar nclass, c, ell, rep, j, n_risk
    real colvector risk, rg, ix, D, fw, p_raw
    real matrix X_risk, lv, info, H, H_chol
    nclass = rows(f.propensity)
    B.ps_factor = cdid_cpt_base_ps_factor(nclass, 1)
    if (f.method == "reg" | nclass == 0) return
    for (c = 1; c <= nclass; c++) {
        ps = f.propensity[c]
        rep = 0
        for (ell = 1; ell <= f.plan.n_consumers; ell++) {
            if (f.point[ell].ps_id == c) {
                rep = ell
                break
            }
        }
        if (rep == 0)
            _cdid_cpt_base_fail("ps_owner", "Fitted propensity class has no consumer mapping.", 498)
        B.workspace_entries = _cdid_base_cells(ps.p, ps.p + 4, B.workspace_limit, B.workspace_entries)
        B.work_ops = _cdid_base_work(f.n_A, B.work_limit, B.work_ops)
        risk = _cdid_base_risk_obs(f.population, f.plan, rep)
        n_risk = rows(risk)
        if (n_risk < 1 | rows(ps.beta) != ps.p | cols(f.X_A) != ps.p)
            _cdid_cpt_base_fail("ps_owner", "Risk observations or treatment columns disagree with fitted propensity.", 498)
        rg = J(n_risk, 1, 0)
        for (j = 1; j <= ps.n_cohorts; j++) {
            ix = _cdid_cpt_base_ix(f.plan.a_cohort[risk] :== ps.cohorts[j])
            if (rows(ix) != ps.pair_count[j])
                _cdid_cpt_base_fail("ps_owner", "Risk cohort count changed after provider preparation.", 498)
            rg[ix] = J(rows(ix), 1, j)
        }
        if (any(rg :== 0))
            _cdid_cpt_base_fail("ps_owner", "Risk observer outside propensity cohort set.", 498)
        D = (rg :== 1)
        fw = f.w_A[risk] :/ ps.rates_c[rg]
        X_risk = f.X_A[risk, .]
        lv = _cdid_ipw_logit_values(X_risk * ps.beta)
        p_raw = lv[, 1]
        info = quadcross(X_risk, fw :* p_raw :* (1 :- p_raw), X_risk)
        H = info / f.n_A
        H_chol = cholesky(H)
        if (missing(H_chol) | any(diagonal(H_chol) :<= 0))
            _cdid_cpt_base_fail("rank", "Static propensity Hessian is not positive definite.", 459)
        B.ps_factor[c].p = ps.p
        B.ps_factor[c].m = ps.n_cohorts
        ix = _cdid_cpt_base_ix(D)
        B.ps_factor[c].sum_a = sum(fw[ix])
        ix = _cdid_cpt_base_ix(!D)
        B.ps_factor[c].sum_b = sum(fw[ix] :* lv[ix, 2])
        B.ps_factor[c].D_T = B.ps_factor[c].sum_a / f.n_F
        B.ps_factor[c].D_M = B.ps_factor[c].sum_b / f.n_F
        B.ps_factor[c].H_chol = H_chol
    }
}

/* Recover bounded p x p PS Cholesky factor H_chol and exact sum_a, sum_b from the SAME
   fitted ps.beta and fixed dated risk design once per consumer PS model in dated provider. */
void _cdid_cpt_base_ps_dated(struct cdid_cpt_baseline scalar B,
    struct cdid_db_provider scalar f)
{
    struct cdid_base_ps scalar ps
    struct cdid_population scalar P
    struct cdid_reference_inputs scalar S
    real scalar nclass, c, pi, j, n_risk, unit_work
    real colvector risk, rg, ix, D, rw, fw, p_raw
    real matrix rr, X_risk, lv, info, H, H_chol
    nclass = rows(f.propensity)
    B.ps_factor = cdid_cpt_base_ps_factor(nclass, 1)
    if (f.method == "reg" | nclass == 0) return
    P = *f.population
    S = *f.reference
    unit_work = 2 * (1 + ceil(ln(rows(f.plan.calendar) + 1) / ln(2)))
    for (c = 1; c <= nclass; c++) {
        ps = f.propensity[c]
        pi = f.plan.consumer[c, 2]
        B.workspace_entries = _cdid_base_cells(ps.p, ps.p + 4, B.workspace_limit, B.workspace_entries)
        B.work_ops = _cdid_base_work(f.n_A * (unit_work + ps.n_cohorts), B.work_limit, B.work_ops)
        rr = _cdid_base_pair_rows(P, f.plan.pair[pi, 1], f.plan.pair[pi, 2], 1, f.n_A)
        rg = J(rows(rr), 1, 0)
        for (j = 1; j <= ps.n_cohorts; j++) {
            ix = _cdid_cpt_base_ix(f.plan.a_cohort[rr[, 1]] :== ps.cohorts[j])
            if (rows(ix)) rg[ix] = J(rows(ix), 1, j)
        }
        ix = _cdid_cpt_base_ix(rg :> 0)
        risk = rr[ix, 1]
        rg = rg[ix]
        n_risk = rows(risk)
        if (n_risk < 1)
            _cdid_cpt_base_fail("ps_owner", "Dated risk support has no pair observers.", 459)
        for (j = 1; j <= ps.n_cohorts; j++) {
            if (sum(rg :== j) != ps.pair_count[j])
                _cdid_cpt_base_fail("ps_owner", "Dated risk cohort count disagrees with fitted propensity.", 498)
        }
        rw = _cdid_db_w(S, f.plan, pi, risk)
        X_risk = _cdid_db_X(S, f.plan, pi, risk, 2)
        if (cols(X_risk) != ps.p | rows(ps.beta) != ps.p)
            _cdid_cpt_base_fail("ps_owner", "Dated treatment design width disagrees with fitted propensity.", 498)
        D = (rg :== 1)
        fw = rw :/ ps.rates_c[rg]
        lv = _cdid_ipw_logit_values(X_risk * ps.beta)
        p_raw = lv[, 1]
        info = quadcross(X_risk, fw :* p_raw :* (1 :- p_raw), X_risk)
        H = info / f.n_A
        H_chol = cholesky(H)
        if (missing(H_chol) | any(diagonal(H_chol) :<= 0))
            _cdid_cpt_base_fail("rank", "Dated propensity Hessian is not positive definite.", 459)
        B.ps_factor[c].p = ps.p
        B.ps_factor[c].m = ps.n_cohorts
        ix = _cdid_cpt_base_ix(D)
        B.ps_factor[c].sum_a = sum(fw[ix])
        ix = _cdid_cpt_base_ix(!D)
        B.ps_factor[c].sum_b = sum(fw[ix] :* lv[ix, 2])
        B.ps_factor[c].D_T = B.ps_factor[c].sum_a / f.n_F
        B.ps_factor[c].D_M = B.ps_factor[c].sum_b / f.n_F
        B.ps_factor[c].H_chol = H_chol
    }
}

real matrix _cdid_cpt_base_occ(struct cdid_population scalar P,
    real colvector cohorts, struct cdid_catalog_opts scalar opts)
{
    struct cdid_pairs scalar pairs
    real scalar nT, nC, n_A, nt, k, s, t, cnt, at
    real colvector a_f, ac, a_ci, y_aug, ik, fin_coh
    real matrix r_grid, Obs, Obs_k, M_k, out
    nT = rows(P.calendar); nC = rows(cohorts); n_A = P.n_A
    if (P.R == 0 | n_A < 1 | nT < 2 | nC < 1 | rows(P.obs_grid) != P.N |
        cols(P.obs_grid) != nT | rows(P.A_to_F) != n_A | max(P.y) - min(P.y) >= . |
        (opts.mem_limit_bytes < . & 8 * n_A * (nT + 4) > opts.mem_limit_bytes)) {
        pairs = cdid_cipw_pairs_opts(P, cohorts, opts)
        return(pairs.occupancy)
    }
    nt = _cdid_check_cohorts(cohorts, 0)
    a_f = P.A_to_F
    if (any(P.cohort_known[a_f] :!= 1)) {
        pairs = cdid_cipw_pairs_opts(P, cohorts, opts)
        return(pairs.occupancy)
    }
    ac = P.cohort[a_f]
    fin_coh = (nt > 0 ? (nC > 1 ? cohorts[1..(nC - 1)] : J(0, 1, .)) : cohorts)
    a_ci = _cdid_population_lookup(ac, fin_coh)
    if (nt > 0 & any(ac :== .)) a_ci[_cdid_cpt_base_ix(ac :== .)] = J(sum(ac :== .), 1, nt)
    if (any(a_ci :< 1) | any(a_ci :> nC)) {
        pairs = cdid_cipw_pairs_opts(P, cohorts, opts)
        return(pairs.occupancy)
    }
    r_grid = P.obs_grid[a_f, .]
    y_aug = 0 \ (P.y :< .)
    Obs = colshape(y_aug[colshape(r_grid, 1) :+ 1], nT)
    out = J(nC * nT * (nT - 1) / 2, 4, 0)
    at = 0
    for (k = 1; k <= nC; k++) {
        ik = _cdid_cpt_base_ix(a_ci :== k)
        if (rows(ik) == 0) continue
        Obs_k = Obs[ik, .]
        M_k = cross(Obs_k, Obs_k)
        for (s = 1; s <= nT - 1; s++) {
            for (t = s + 1; t <= nT; t++) {
                cnt = M_k[s, t]
                if (cnt > 0) {
                    at++
                    out[at, 1] = k
                    out[at, 2] = s
                    out[at, 3] = t
                    out[at, 4] = cnt
                }
            }
        }
    }
    if (at == 0) return(J(0, 4, .))
    return(out[1..at, .])
}

/* Validate the actual retained original catalogue plus BOTH auxiliary families.
   A public post target is not a primitive comparison row. This metadata pass uses
   catalogue construction once; it never fits or changes any nuisance. */
void _cdid_cpt_base_union_check(struct cdid_cpt_baseline scalar B,
    struct cdid_population scalar P,struct cdid_public_request scalar req,
    struct cdid_baseline_plan_state scalar original,
    struct cdid_reference_catalog scalar cat,struct cdid_reference_inputs scalar S)
{
    struct cdid_catalog_opts scalar opts
    struct cdid_cat scalar C
    real matrix raw, occ
    real colvector ptr,K,actualK
    real scalar i,j,at,M,Q,c,nk,used
    opts=cdid_catalog_default_opts()
    opts.mem_limit_bytes=8*B.workspace_limit
    opts.work_cap=B.work_limit
    occ=_cdid_cpt_base_occ(P,req.cohorts,opts)
    C=cdid_catalog_build_opts(req.calendar,req.cohorts,req.target_cohorts,
        req.anticipation,occ,req.control_group,req.post,opts)
    M=rows(C.r_g);Q=rows(req.aux)
    raw=(C.cohorts[C.r_g],C.calendar[C.r_s],C.calendar[C.r_t])\req.aux
    B.work_ops=_cdid_base_work(C.work+M+Q,B.work_limit,B.work_ops)
    used=_cdid_base_cells(4,M+Q+1,B.workspace_limit,B.workspace_entries)
    if(B.kind=="static") {
        if(original.M!=M | original.Q!=Q | rows(original.raw)!=M+Q)
            _cdid_cpt_base_fail("identity","Provider does not retain the complete original comparison and auxiliary union.",498)
        if(!all(raw:==(original.cohorts[original.raw[,1]],
            original.calendar[original.raw[,2]],original.calendar[original.raw[,3]])))
            _cdid_cpt_base_fail("identity","Original comparison or auxiliary endpoints differ from the request.",498)
    }
    else {
        if(cat.M!=M | cat.Q!=Q | cat.n_raw!=M+Q | rows(S.view_meta)!=M+Q)
            _cdid_cpt_base_fail("identity","Dated provider does not retain the complete original comparison and auxiliary union.",498)
        if(!all(raw:==S.view_meta[,1..3]))
            _cdid_cpt_base_fail("identity","Dated original comparison or auxiliary endpoints differ from the request.",498)
    }
    for(i=1;i<=M+Q;i++) {
        if(i<=M) {
            j=C.r_pool[i]
            K=C.pool_entries[C.pool_start[j]..C.pool_start[j]+C.pool_size[j]-1]
        }
        else {
            j=i-M
            K=req.aux_K_index[req.aux_K_ptr[j]..req.aux_K_ptr[j+1]-1]
        }
        B.work_ops=_cdid_base_work(rows(K),B.work_limit,B.work_ops)
        if(B.kind=="static")actualK=original.raw_K_index[original.raw_K_ptr[i]..original.raw_K_ptr[i+1]-1]
        else {
            c=cat.raw_consumer[i]
            if(!c)continue
            actualK=cat.K_index[cat.K_ptr[c]..cat.K_ptr[c+1]-1]
        }
        if(rows(K)!=rows(actualK))
            _cdid_cpt_base_fail("identity","A complete-union control set differs from the request.",498)
        if(!all(K:==actualK))
            _cdid_cpt_base_fail("identity","A complete-union control set differs from the request.",498)
    }
}

/* Prepare baseline conditional-pretrend piece producer from a static cdid_baseline_provider
   that was fitted on the full original post + both-auxiliary union. */
struct cdid_cpt_baseline scalar cdid_cpt_baseline_prepare_static(
    pointer(struct cdid_baseline_provider scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    struct cdid_cpt_baseline scalar B
    struct cdid_baseline_provider scalar f
    struct cdid_cpt_base_cell_meta scalar cm
    struct cdid_reference_catalog scalar no_cat
    struct cdid_reference_inputs scalar no_ref
    real scalar fam, i, r, q, ell, c, b, s, m, p_trt, p_or, icpt, n_cells, sep
    real colvector sel_rows, s_vec, D_vec
    string colvector coords
    if (provider == NULL)
        _cdid_cpt_base_fail("owner", "Static baseline provider pointer is required.", 498)
    f = *provider
    if (workcap <. & (workcap < 1 | workcap != floor(workcap)))
        _cdid_cpt_base_fail("bound", "workcap must be . or a positive integer.", 498)
    if (f.plan.Q != rows(req.aux) |
        rows(f.plan.calendar) != rows(req.calendar) | !all(f.plan.calendar :== req.calendar) |
        rows(f.plan.cohorts) != rows(req.cohorts) | !all(f.plan.cohorts :== req.cohorts) |
        f.population.N != f.n_F | f.population.n_A != f.n_A)
        _cdid_cpt_base_fail("identity", "Static provider must be fitted on the full original post + both-auxiliary union of the request.", 498)
    if (req.base_period != "varying" & req.base_period != "universal")
        _cdid_cpt_base_fail("request", "Invalid request base_period.", 498)
    B.kind = "static"
    B.method = f.method
    B.base_period = req.base_period
    B.n_F = f.n_F
    B.n_A = f.n_A
    B.chunk_rows = f.chunk_rows
    B.workspace_entries = f.workspace_entries
    B.workspace_limit = f.workspace_limit
    B.work_limit = workcap
    B.work_ops = _cdid_base_work(0, workcap, f.plan.work_ops)
    B.static_provider = provider
    B.dated_provider = NULL
    _cdid_cpt_base_union_check(B,f.population,req,f.plan,no_cat,no_ref)
    _cdid_cpt_base_ps_static(B, f)

    fam = (req.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_cpt_base_ix((req.selection[, 1] :== fam) :& (req.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    B.workspace_entries = _cdid_base_cells(16, n_cells + 1, B.workspace_limit, B.workspace_entries)
    B.work_ops = _cdid_base_work(rows(req.selection) + n_cells, B.work_limit, B.work_ops)
    B.n_cells = n_cells
    B.n_unavailable = 0
    B.cells = cdid_cpt_base_cell_meta(n_cells, 1)
    if (n_cells == 0) {
        B.status = "unavailable"
        B.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + req.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        return(B)
    }

    sep = (f.method != "ipw" & f.outcome.outcome_pooling == "separate")
    icpt = (f.method == "ipw" ? 0 : f.outcome.intercept)
    for (i = 1; i <= n_cells; i++) {
        r = sel_rows[i]
        cm.g = req.selection[r, 2]
        cm.t = req.selection[r, 3]
        cm.from = req.selection[r, 5]
        cm.q = req.selection[r, 6]
        cm.sign = req.selection[r, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.s_index = 0
        cm.t2_index = 0
        cm.consumer = 0
        cm.recipe = 0
        cm.design_index = 0
        cm.supported = 0
        cm.is_legacy_ipw = (f.method == "ipw")
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = f.chunk_rows
        cm.route = (f.method == "ipw" ?
            "CP02 risk-set IPW (A frame, observed-pair rates and one propensity logit per auxiliary link)" :
            "static baseline " + f.method)
        if (req.selection[r, 8] != 0) {
            cm.status = "unavailable"
            cm.reason = req.selection_reason[r]
            B.n_unavailable = B.n_unavailable + 1
            B.cells[i] = cm
            continue
        }
        q = cm.q
        cm.s_date = req.aux[q, 2]
        cm.t2_date = req.aux[q, 3]
        cm.s_index = f.plan.raw[f.plan.M + q, 2]
        cm.t2_index = f.plan.raw[f.plan.M + q, 3]
        if (!f.plan.aux_supported[q] | f.plan.aux_consumer[q] < 1) {
            cm.status = "unavailable"
            cm.reason = "unsupported auxiliary comparison: " + f.plan.aux_reason[q]
            B.n_unavailable = B.n_unavailable + 1
            B.cells[i] = cm
            continue
        }
        ell = f.plan.aux_consumer[q]
        cm.consumer = ell
        cm.supported = 1
        cm.status = "available"
        cm.reason = ""
        cm.estimate = f.points[ell]
        cm.estimate_scalar = cm.sign * f.points[ell]
        if (f.method == "reg") {
            b = f.plan.consumer[ell, 3]
            s = f.plan.consumer[ell, 4]
            p_or = (sep ? rows(f.outcome.blocks[b].beta) : rows(f.outcome.beta))
            coords = _cdid_cpt_base_coords_or(s, icpt, b, p_or, sep)
            s_vec = 1
            D_vec = f.outcome.link_d[ell]
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
            cm.chunk_rows = min((f.n_F, f.chunk_rows))
            cm.n_chunks = ceil(f.n_F / cm.chunk_rows)
        }
        else {
            c = f.point[ell].ps_id
            b = f.plan.consumer[ell, 3]
            m = f.propensity[c].n_cohorts
            p_trt = f.propensity[c].p
            p_or = (f.method == "dr" ? (sep ? rows(f.outcome.blocks[b].beta) : rows(f.outcome.beta)) : 0)
            coords = _cdid_cpt_base_coords_two_arm(cm.is_legacy_ipw, f.method == "dr", m, p_trt, icpt, b, p_or, sep)
            s_vec = (1 \ -1)
            D_vec = (B.ps_factor[c].D_T \ B.ps_factor[c].D_M)
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
            if (cm.is_legacy_ipw) {
                cm.chunk_rows = min((f.n_A, f.chunk_rows))
                cm.n_chunks = ceil(f.n_A / cm.chunk_rows)
            }
            else {
                cm.chunk_rows = min((f.n_F, f.chunk_rows))
                cm.n_chunks = ceil(f.n_F / cm.chunk_rows)
            }
        }
        B.cells[i] = cm
    }

    if (B.n_unavailable > 0) {
        B.status = "unavailable"
        B.reason = "The complete requested test is unavailable: " + strofreal(B.n_unavailable, "%21.0f") + " declared cell(s) are unavailable. No available-only subset is tested."
    }
    else {
        B.status = "available"
        B.reason = ""
    }
    return(B)
}

/* Prepare baseline conditional-pretrend piece producer from a dated cdid_db_provider
   that was fitted on the full original post + both-auxiliary union. */
struct cdid_cpt_baseline scalar cdid_cpt_baseline_prepare_dated(
    pointer(struct cdid_db_provider scalar) scalar provider,
    struct cdid_public_request scalar req,
    real scalar workcap)
{
    struct cdid_cpt_baseline scalar B
    struct cdid_db_provider scalar f
    struct cdid_population scalar P
    struct cdid_reference_inputs scalar S
    struct cdid_cpt_base_cell_meta scalar cm
    struct cdid_baseline_plan_state scalar no_original
    real scalar fam, i, r, q, ell, c, b, s, m, p_trt, p_or, icpt, n_cells, eff_workcap, sep
    real colvector sel_rows, s_vec, D_vec
    string colvector coords
    if (provider == NULL)
        _cdid_cpt_base_fail("owner", "Dated baseline provider pointer is required.", 498)
    f = *provider
    if (f.population == NULL | f.reference == NULL)
        _cdid_cpt_base_fail("owner", "Dated baseline provider must retain population and reference pointers.", 498)
    P = *f.population
    S = *f.reference
    cdid_ref_activate(S)
    if (workcap <. & (workcap < 1 | workcap != floor(workcap)))
        _cdid_cpt_base_fail("bound", "workcap must be . or a positive integer.", 498)
    eff_workcap = (workcap <. ? min((workcap, f.plan.work_limit)) : f.plan.work_limit)
    if (f.catalog.Q != rows(req.aux) |
        rows(f.plan.calendar) != rows(req.calendar) | !all(f.plan.calendar :== req.calendar) |
        rows(f.plan.cohorts) != rows(req.cohorts) | !all(f.plan.cohorts :== req.cohorts) |
        P.N != f.n_F | P.n_A != f.n_A)
        _cdid_cpt_base_fail("identity", "Dated provider must be fitted on the full original post + both-auxiliary union of the request.", 498)
    if (req.base_period != "varying" & req.base_period != "universal")
        _cdid_cpt_base_fail("request", "Invalid request base_period.", 498)
    B.kind = "dated"
    B.method = f.method
    B.base_period = req.base_period
    B.n_F = f.n_F
    B.n_A = f.n_A
    B.chunk_rows = f.chunk_rows
    B.workspace_entries = f.workspace_entries
    B.workspace_limit = f.workspace_limit
    B.work_limit = eff_workcap
    B.work_ops = _cdid_base_work(0, eff_workcap, f.work_ops)
    B.static_provider = NULL
    B.dated_provider = provider
    _cdid_cpt_base_union_check(B,P,req,no_original,f.catalog,S)
    _cdid_cpt_base_ps_dated(B, f)

    fam = (req.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_cpt_base_ix((req.selection[, 1] :== fam) :& (req.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    B.workspace_entries = _cdid_base_cells(16, n_cells + 1, B.workspace_limit, B.workspace_entries)
    B.work_ops = _cdid_base_work(rows(req.selection) + n_cells, B.work_limit, B.work_ops)
    B.n_cells = n_cells
    B.n_unavailable = 0
    B.cells = cdid_cpt_base_cell_meta(n_cells, 1)
    if (n_cells == 0) {
        B.status = "unavailable"
        B.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + req.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        return(B)
    }

    sep = (f.method != "ipw" & f.outcome.outcome_pooling == "separate")
    icpt = (f.method == "ipw" ? 0 : f.outcome.intercept)
    for (i = 1; i <= n_cells; i++) {
        r = sel_rows[i]
        cm.g = req.selection[r, 2]
        cm.t = req.selection[r, 3]
        cm.from = req.selection[r, 5]
        cm.q = req.selection[r, 6]
        cm.sign = req.selection[r, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.s_index = 0
        cm.t2_index = 0
        cm.consumer = 0
        cm.recipe = 0
        cm.design_index = 0
        cm.supported = 0
        cm.is_legacy_ipw = 0
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = f.chunk_rows
        cm.route = "dated D60 baseline " + f.method + " (reference-date recipes)"
        if (req.selection[r, 8] != 0) {
            cm.status = "unavailable"
            cm.reason = req.selection_reason[r]
            B.n_unavailable = B.n_unavailable + 1
            B.cells[i] = cm
            continue
        }
        q = cm.q
        cm.s_date = req.aux[q, 2]
        cm.t2_date = req.aux[q, 3]
        if (!f.catalog.aux_supported[q] | f.catalog.aux_consumer[q] < 1) {
            cm.status = "unavailable"
            cm.reason = "unsupported auxiliary comparison: " + f.catalog.aux_reason[q]
            B.n_unavailable = B.n_unavailable + 1
            B.cells[i] = cm
            continue
        }
        ell = f.catalog.aux_consumer[q]
        cm.consumer = ell
        cm.recipe = f.catalog.consumer_recipe[ell]
        cm.design_index = S.recipes[cm.recipe, 1]
        cm.s_index = f.plan.pair[f.plan.consumer[ell, 2], 1]
        cm.t2_index = f.plan.pair[f.plan.consumer[ell, 2], 2]
        cm.supported = 1
        cm.status = "available"
        cm.reason = ""
        cm.estimate = f.points[ell]
        cm.estimate_scalar = cm.sign * f.points[ell]
        cm.chunk_rows = min((f.n_F, f.chunk_rows))
        cm.n_chunks = ceil(f.n_F / cm.chunk_rows)
        if (f.method == "reg") {
            b = f.plan.consumer[ell, 3]
            s = f.plan.consumer[ell, 4]
            p_or = (sep ? rows(f.outcome.blocks[b].beta) : rows(f.outcome.beta))
            coords = _cdid_cpt_base_coords_or(s, icpt, b, p_or, sep)
            s_vec = 1
            D_vec = f.outcome.link_d[ell]
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        else {
            c = f.point[ell].ps_id
            b = f.plan.consumer[ell, 3]
            m = f.propensity[c].n_cohorts
            p_trt = f.propensity[c].p
            p_or = (f.method == "dr" ? (sep ? rows(f.outcome.blocks[b].beta) : rows(f.outcome.beta)) : 0)
            coords = _cdid_cpt_base_coords_two_arm(0, f.method == "dr", m, p_trt, icpt, b, p_or, sep)
            s_vec = (1 \ -1)
            D_vec = (B.ps_factor[c].D_T \ B.ps_factor[c].D_M)
            cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, s_vec, D_vec, coords)
        }
        B.cells[i] = cm
    }

    if (B.n_unavailable > 0) {
        B.status = "unavailable"
        B.reason = "The complete requested test is unavailable: " + strofreal(B.n_unavailable, "%21.0f") + " declared cell(s) are unavailable. No available-only subset is tested."
    }
    else {
        B.status = "available"
        B.reason = ""
    }
    return(B)
}

struct cdid_cpt_spec scalar cdid_cpt_baseline_cell_spec(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx)
{
    if (cell_idx < 1 | cell_idx > B.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_base_fail("cell", "Invalid declared cell index.", 498)
    if (!B.cells[cell_idx].supported)
        _cdid_cpt_base_fail("cell", "Requested cell is unavailable: " + B.cells[cell_idx].reason, 459)
    return(B.cells[cell_idx].spec)
}

real scalar cdid_cpt_baseline_cell_chunks(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx)
{
    if (cell_idx < 1 | cell_idx > B.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_base_fail("cell", "Invalid declared cell index.", 498)
    if (!B.cells[cell_idx].supported)
        _cdid_cpt_base_fail("cell", "Requested cell is unavailable: " + B.cells[cell_idx].reason, 459)
    return(B.cells[cell_idx].n_chunks)
}

/* Internal static piece evaluator on arbitrary unique original-F rows. */
struct cdid_cpt_piece scalar _cdid_cpt_base_piece_static_rows(
    struct cdid_cpt_baseline scalar B,
    struct cdid_baseline_provider scalar f,
    struct cdid_cpt_base_cell_meta scalar cm,
    real colvector rowsF,
    | transmorphic scalar cache)
{
    struct cdid_base_ps scalar ps
    struct cdid_cpt_base_ps_factor scalar pf
    real scalar n, ell, c, pi, b, s, L, m, p_trt, p_or, icpt, Dl, j, n_obs, scale, target_coh, sep, use_cache
    real colvector aa, gi, pos, ai_obs, dy, gj, D_ind, fw, odds, act, rho_obs, q_w, beta_or, ix_a, ac_sub, ix_j
    real rowvector kept
    real matrix rr, xr, score, Z_trt, lv, lam_up, xb, xi, betaIF, alphaIF, Z_or, d_or, nu, rh, lam, pim
    string scalar bkey, akey
    use_cache = (args() >= 5)
    n = rows(rowsF)
    ell = cm.consumer
    if (f.method == "reg") {
        b = f.plan.consumer[ell, 3]
        s = f.plan.consumer[ell, 4]
        pi = f.plan.support[s, 1]
        L = f.plan.pair[pi, 2] - f.plan.pair[pi, 1]
        icpt = f.outcome.intercept
        sep = (f.outcome.outcome_pooling == "separate")
        kept = (sep ? f.outcome.blocks[b].kept : f.outcome.kept)
        beta_or = (sep ? f.outcome.blocks[b].beta : f.outcome.beta)
        p_or = rows(beta_or)
        Dl = cm.spec.Dl
        xi = _cdid_base_rate_if(f, s, rowsF)
        if (p_or | icpt) {
            if (use_cache) {
                bkey = (sep ? "b:" + strofreal(b) : "b:0")
                if (asarray_contains(cache, bkey)) betaIF = asarray(cache, bkey)
                else {
                    betaIF = (sep ? _cdid_base_blk_beta_if(f, b, rowsF) : _cdid_base_beta_if(f, rowsF))
                    asarray(cache, bkey, betaIF)
                }
            }
            else betaIF = (sep ? _cdid_base_blk_beta_if(f, b, rowsF) : _cdid_base_beta_if(f, rowsF))
        }
        else betaIF = J(n, 0, 0)
        if (icpt) {
            if (use_cache) {
                akey = "a:" + strofreal(b)
                if (asarray_contains(cache, akey)) alphaIF = asarray(cache, akey)
                else {
                    alphaIF = _cdid_base_alpha_if(f, b, rowsF, betaIF)
                    asarray(cache, akey, alphaIF)
                }
            }
            else alphaIF = _cdid_base_alpha_if(f, b, rowsF, betaIF)
            xi = (xi, alphaIF)
        }
        if (p_or) xi = (xi, betaIF)
        target_coh = f.plan.active_cohorts[f.plan.support[s, 2]]
        rr = _cdid_base_frame_pair(f, pi, rowsF, target_coh, 1)
        n_obs = rows(rr)
        nu = J(n, 1, 0); rh = J(n, 1, 0)
        lam = J(n, Dl, 0); pim = J(n, Dl, 0)
        if (n_obs) {
            pos = rr[, 1]; ai_obs = rr[, 2]; dy = rr[, 3]
            q_w = f.w_A[ai_obs] / f.outcome.rates_c[s]
            d_or = (p_or ? L :* f.Z_A[ai_obs, kept] : J(n_obs, 0, 0))
            rho_obs = dy :- (icpt ? f.outcome.alpha[b] : 0) :- (p_or ? d_or * beta_or : J(n_obs, 1, 0))
            nu[pos, 1] = q_w
            rh[pos, 1] = rho_obs
            lam[pos, 1] = J(n_obs, 1, (f.outcome.rates_active[s] ? -1 / f.outcome.rates[s] : 0))
            if (icpt) pim[pos, 2] = J(n_obs, 1, -1)
            if (p_or) pim[pos, (1 + icpt) :+ (1..p_or)] = -d_or
        }
        return(cdid_cpt_make_piece(rowsF, nu, rh, lam, J(0, 0, 0), pim, J(0, 0, 0), 1 / cm.spec.D, xi, cm.spec))
    }

    c = f.point[ell].ps_id
    ps = f.propensity[c]
    pf = B.ps_factor[c]
    pi = f.plan.consumer[ell, 2]
    b = f.plan.consumer[ell, 3]
    L = f.plan.pair[pi, 2] - f.plan.pair[pi, 1]
    m = ps.n_cohorts
    p_trt = ps.p
    icpt = (f.method == "dr" ? f.outcome.intercept : 0)
    sep = (f.method == "dr" & f.outcome.outcome_pooling == "separate")
    kept = (f.method == "dr" ? (sep ? f.outcome.blocks[b].kept : f.outcome.kept) : J(1, 0, 0))
    beta_or = (f.method == "dr" ? (sep ? f.outcome.blocks[b].beta : f.outcome.beta) : J(0, 1, 0))
    p_or = rows(beta_or)
    scale = f.n_F / f.n_A

    aa = f.plan.F_to_A[rowsF]
    gi = J(n, 1, 0)
    xr = J(n, m, 0)
    ix_a = _cdid_cpt_base_ix(aa :> 0)
    if (rows(ix_a)) {
        ac_sub = f.plan.a_cohort[aa[ix_a]]
        for (j = 1; j <= m; j++) {
            ix_j = _cdid_cpt_base_ix(ac_sub :== ps.cohorts[j])
            if (rows(ix_j)) {
                gi[ix_a[ix_j]] = J(rows(ix_j), 1, j)
                xr[ix_a[ix_j], j] = (-ps.rates[j] / ps.q[j]) :* f.w_A[aa[ix_a[ix_j]]]
            }
        }
    }
    rr = _cdid_base_frame_pair(f, pi, rowsF, ps.cohorts, 1)
    n_obs = rows(rr)
    score = J(n, p_trt, 0)
    pos = J(0, 1, 0); D_ind = J(0, 1, 0); fw = J(0, 1, 0)
    Z_trt = J(0, p_trt, 0); odds = J(0, 1, 0); act = J(0, 1, 0)
    rho_obs = J(0, 1, 0); lam_up = J(0, m, 0); d_or = J(0, p_or, 0)
    if (n_obs) {
        pos = rr[, 1]; ai_obs = rr[, 2]; dy = rr[, 3]
        gj = gi[pos]
        lam_up = J(n_obs, m, 0)
        for (j = 1; j <= m; j++) {
            ix_j = _cdid_cpt_base_ix(gj :== j)
            if (rows(ix_j)) {
                xr[pos[ix_j], j] = ((1 - ps.rates[j]) / ps.q[j]) :* f.w_A[ai_obs[ix_j]]
                if (ps.rates_active[j]) lam_up[ix_j, j] = J(rows(ix_j), 1, -1 / ps.rates[j])
            }
        }
        Z_trt = f.X_A[ai_obs, .]
        lv = _cdid_ipw_logit_values(Z_trt * ps.beta)
        D_ind = (gj :== 1)
        fw = f.w_A[ai_obs] :/ ps.rates_c[gj]
        score[pos, .] = (fw :* (D_ind :- lv[, 1])) :* Z_trt
        odds = lv[, 2]
        act = lv[, 3]
        if (f.method == "dr") {
            Z_or = (p_or ? f.Z_A[ai_obs, kept] : J(n_obs, 0, 0))
            d_or = L :* Z_or
            rho_obs = dy :- (icpt ? f.outcome.alpha[b] : 0) :- (p_or ? d_or * beta_or : J(n_obs, 1, 0))
        }
        else {
            rho_obs = dy
        }
    }
    xb = _cdid_cpt_base_chol_solve(pf.H_chol, score + xr * (ps.A_beta_r'))
    xi = scale :* (xr, xb)
    if (f.method == "dr") {
        if (p_or | icpt) {
            if (use_cache) {
                bkey = (sep ? "b:" + strofreal(b) : "b:0")
                if (asarray_contains(cache, bkey)) betaIF = asarray(cache, bkey)
                else {
                    betaIF = (sep ? _cdid_base_blk_beta_if(f, b, rowsF) : _cdid_base_beta_if(f, rowsF))
                    asarray(cache, bkey, betaIF)
                }
            }
            else betaIF = (sep ? _cdid_base_blk_beta_if(f, b, rowsF) : _cdid_base_beta_if(f, rowsF))
        }
        else betaIF = J(n, 0, 0)
        if (icpt) {
            if (use_cache) {
                akey = "a:" + strofreal(b)
                if (asarray_contains(cache, akey)) alphaIF = asarray(cache, akey)
                else {
                    alphaIF = _cdid_base_alpha_if(f, b, rowsF, betaIF)
                    asarray(cache, akey, alphaIF)
                }
            }
            else alphaIF = _cdid_base_alpha_if(f, b, rowsF, betaIF)
            xi = (xi, alphaIF)
        }
        if (p_or) xi = (xi, betaIF)
    }
    return(cdid_cpt_two_arm_piece(rowsF, pos, D_ind, fw, Z_trt, odds, act, rho_obs,
        lam_up, xi, pf.D_T, pf.D_M, m, p_trt, d_or, icpt, cm.spec))
}

/* Internal dated piece evaluator on arbitrary unique original-F rows. */
struct cdid_cpt_piece scalar _cdid_cpt_base_piece_dated_rows(
    struct cdid_cpt_baseline scalar B,
    struct cdid_db_provider scalar f,
    struct cdid_cpt_base_cell_meta scalar cm,
    real colvector rowsF,
    | transmorphic scalar cache)
{
    struct cdid_base_ps scalar ps
    struct cdid_cpt_base_ps_factor scalar pf
    pointer(struct cdid_reference_inputs scalar) scalar pS
    real scalar n, ell, c, pi, b, s, L, m, p_trt, p_or, icpt, Dl, j, n_obs, scale, target_coh, sep, use_cache
    real colvector aa, gi, pos, ai_obs, dy, gj, D_ind, rw, fw, odds, act, rho_obs, q_w, w_gi, beta_or, ix_a, ac_sub, ix_gi, ix_j
    real rowvector kept
    real matrix rr, xr, score, Z_trt, lv, lam_up, xb, xi, betaIF, alphaIF, Z_or, d_or, nu, rh, lam, pim
    string scalar bkey, akey
    use_cache = (args() >= 5)
    pS = f.reference
    n = rows(rowsF)
    ell = cm.consumer
    if (f.method == "reg") {
        b = f.plan.consumer[ell, 3]
        s = f.plan.consumer[ell, 4]
        pi = f.plan.support[s, 1]
        L = f.plan.pair[pi, 2] - f.plan.pair[pi, 1]
        icpt = f.outcome.intercept
        sep = (f.outcome.outcome_pooling == "separate")
        kept = (sep ? f.outcome.blocks[b].kept : f.outcome.kept)
        beta_or = (sep ? f.outcome.blocks[b].beta : f.outcome.beta)
        p_or = rows(beta_or)
        Dl = cm.spec.Dl
        xi = _cdid_db_rate_if(f, s, rowsF)
        if (p_or | icpt) {
            if (use_cache) {
                bkey = (sep ? "b:" + strofreal(b) : "b:0")
                if (asarray_contains(cache, bkey)) betaIF = asarray(cache, bkey)
                else {
                    betaIF = (sep ? _cdid_db_blk_beta_if(f, b, rowsF) : _cdid_db_beta_if(f, rowsF))
                    asarray(cache, bkey, betaIF)
                }
            }
            else betaIF = (sep ? _cdid_db_blk_beta_if(f, b, rowsF) : _cdid_db_beta_if(f, rowsF))
        }
        else betaIF = J(n, 0, 0)
        if (icpt) {
            if (use_cache) {
                akey = "a:" + strofreal(b)
                if (asarray_contains(cache, akey)) alphaIF = asarray(cache, akey)
                else {
                    alphaIF = _cdid_db_alpha_if(f, b, rowsF, betaIF)
                    asarray(cache, akey, alphaIF)
                }
            }
            else alphaIF = _cdid_db_alpha_if(f, b, rowsF, betaIF)
            xi = (xi, alphaIF)
        }
        if (p_or) xi = (xi, betaIF)
        target_coh = f.plan.active_cohorts[f.plan.support[s, 2]]
        rr = _cdid_db_frame_pair(f, pi, rowsF, target_coh, 1)
        n_obs = rows(rr)
        nu = J(n, 1, 0); rh = J(n, 1, 0)
        lam = J(n, Dl, 0); pim = J(n, Dl, 0)
        if (n_obs) {
            pos = rr[, 1]; ai_obs = rr[, 2]; dy = rr[, 3]
            q_w = _cdid_db_w(*pS, f.plan, pi, ai_obs) / f.outcome.rates_c[s]
            d_or = (p_or ? L :* _cdid_db_Z(*pS, f.plan, pi, ai_obs)[, kept] : J(n_obs, 0, 0))
            rho_obs = dy :- (icpt ? f.outcome.alpha[b] : 0) :- (p_or ? d_or * beta_or : J(n_obs, 1, 0))
            nu[pos, 1] = q_w
            rh[pos, 1] = rho_obs
            lam[pos, 1] = J(n_obs, 1, (f.outcome.rates_active[s] ? -1 / f.outcome.rates[s] : 0))
            if (icpt) pim[pos, 2] = J(n_obs, 1, -1)
            if (p_or) pim[pos, (1 + icpt) :+ (1..p_or)] = -d_or
        }
        return(cdid_cpt_make_piece(rowsF, nu, rh, lam, J(0, 0, 0), pim, J(0, 0, 0), 1 / cm.spec.D, xi, cm.spec))
    }

    c = f.point[ell].ps_id
    ps = f.propensity[c]
    pf = B.ps_factor[c]
    pi = f.plan.consumer[ell, 2]
    b = f.plan.consumer[ell, 3]
    L = f.plan.pair[pi, 2] - f.plan.pair[pi, 1]
    m = ps.n_cohorts
    p_trt = ps.p
    icpt = (f.method == "dr" ? f.outcome.intercept : 0)
    sep = (f.method == "dr" & f.outcome.outcome_pooling == "separate")
    kept = (f.method == "dr" ? (sep ? f.outcome.blocks[b].kept : f.outcome.kept) : J(1, 0, 0))
    beta_or = (f.method == "dr" ? (sep ? f.outcome.blocks[b].beta : f.outcome.beta) : J(0, 1, 0))
    p_or = rows(beta_or)
    scale = f.n_F / f.n_A

    aa = f.plan.F_to_A[rowsF]
    gi = J(n, 1, 0)
    ix_a = _cdid_cpt_base_ix(aa :> 0)
    if (rows(ix_a)) {
        ac_sub = f.plan.a_cohort[aa[ix_a]]
        for (j = 1; j <= m; j++) {
            ix_j = _cdid_cpt_base_ix(ac_sub :== ps.cohorts[j])
            if (rows(ix_j)) gi[ix_a[ix_j]] = J(rows(ix_j), 1, j)
        }
    }
    rr = _cdid_db_frame_pair(f, pi, rowsF, ps.cohorts, 1)
    n_obs = rows(rr)
    xr = J(n, m, 0)
    ix_gi = _cdid_cpt_base_ix(gi :> 0)
    if (rows(ix_gi)) {
        w_gi = _cdid_db_w(*pS, f.plan, pi, aa[ix_gi])
        for (j = 1; j <= m; j++) {
            ix_j = _cdid_cpt_base_ix(gi[ix_gi] :== j)
            if (rows(ix_j)) xr[ix_gi[ix_j], j] = (-ps.rates[j] / ps.q[j]) :* w_gi[ix_j]
        }
    }
    score = J(n, p_trt, 0)
    pos = J(0, 1, 0); D_ind = J(0, 1, 0); fw = J(0, 1, 0)
    Z_trt = J(0, p_trt, 0); odds = J(0, 1, 0); act = J(0, 1, 0)
    rho_obs = J(0, 1, 0); lam_up = J(0, m, 0); d_or = J(0, p_or, 0)
    if (n_obs) {
        pos = rr[, 1]; ai_obs = rr[, 2]; dy = rr[, 3]
        gj = gi[pos]
        rw = _cdid_db_w(*pS, f.plan, pi, ai_obs)
        lam_up = J(n_obs, m, 0)
        for (j = 1; j <= m; j++) {
            ix_j = _cdid_cpt_base_ix(gj :== j)
            if (rows(ix_j)) {
                xr[pos[ix_j], j] = ((1 - ps.rates[j]) / ps.q[j]) :* rw[ix_j]
                if (ps.rates_active[j]) lam_up[ix_j, j] = J(rows(ix_j), 1, -1 / ps.rates[j])
            }
        }
        Z_trt = _cdid_db_X(*pS, f.plan, pi, ai_obs, 2)
        lv = _cdid_ipw_logit_values(Z_trt * ps.beta)
        D_ind = (gj :== 1)
        fw = rw :/ ps.rates_c[gj]
        score[pos, .] = (fw :* (D_ind :- lv[, 1])) :* Z_trt
        odds = lv[, 2]
        act = lv[, 3]
        if (f.method == "dr") {
            Z_or = (p_or ? _cdid_db_Z(*pS, f.plan, pi, ai_obs)[, kept] : J(n_obs, 0, 0))
            d_or = L :* Z_or
            rho_obs = dy :- (icpt ? f.outcome.alpha[b] : 0) :- (p_or ? d_or * beta_or : J(n_obs, 1, 0))
        }
        else {
            rho_obs = dy
        }
    }
    xb = _cdid_cpt_base_chol_solve(pf.H_chol, score + xr * (ps.A_beta_r'))
    xi = scale :* (xr, xb)
    if (f.method == "dr") {
        if (p_or | icpt) {
            if (use_cache) {
                bkey = (sep ? "b:" + strofreal(b) : "b:0")
                if (asarray_contains(cache, bkey)) betaIF = asarray(cache, bkey)
                else {
                    betaIF = (sep ? _cdid_db_blk_beta_if(f, b, rowsF) : _cdid_db_beta_if(f, rowsF))
                    asarray(cache, bkey, betaIF)
                }
            }
            else betaIF = (sep ? _cdid_db_blk_beta_if(f, b, rowsF) : _cdid_db_beta_if(f, rowsF))
        }
        else betaIF = J(n, 0, 0)
        if (icpt) {
            if (use_cache) {
                akey = "a:" + strofreal(b)
                if (asarray_contains(cache, akey)) alphaIF = asarray(cache, akey)
                else {
                    alphaIF = _cdid_db_alpha_if(f, b, rowsF, betaIF)
                    asarray(cache, akey, alphaIF)
                }
            }
            else alphaIF = _cdid_db_alpha_if(f, b, rowsF, betaIF)
            xi = (xi, alphaIF)
        }
        if (p_or) xi = (xi, betaIF)
    }
    return(cdid_cpt_two_arm_piece(rowsF, pos, D_ind, fw, Z_trt, odds, act, rho_obs,
        lam_up, xi, pf.D_T, pf.D_M, m, p_trt, d_or, icpt, cm.spec))
}

/* Evaluate one bounded piece on caller-supplied unique original-F row positions. */
struct cdid_cpt_piece scalar cdid_cpt_baseline_piece_rows(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx, real colvector rowsF)
{
    pointer(struct cdid_baseline_provider scalar) scalar pfs
    pointer(struct cdid_db_provider scalar) scalar pfd
    struct cdid_cpt_base_cell_meta scalar cm
    real scalar n, R, Dl, used, is_inc
    if (cell_idx < 1 | cell_idx > B.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_base_fail("replay", "Invalid declared cell index.", 498)
    cm = B.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_base_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    n = rows(rowsF)
    is_inc = (n <= 1 ? 1 : all(rowsF[(2..n)] :> rowsF[(1..(n - 1))]))
    if (n < 1 | cols(rowsF) != 1 | missing(rowsF) | any(rowsF :!= floor(rowsF)) |
        any(rowsF :< 1) | any(rowsF :> B.n_F) | (!is_inc & rows(uniqrows(sort(rowsF, 1))) != n) |
        n > B.chunk_rows)
        _cdid_cpt_base_fail("replay", "rowsF must be 1..chunk_rows unique original-F positions in 1..N_F.", 498)
    R = cm.spec.R; Dl = cm.spec.Dl
    used = _cdid_base_cells(12 * n, 2 * R + 3 * Dl + 8, B.workspace_limit, B.workspace_entries)
    if (B.kind == "static") {
        if (B.static_provider == NULL)
            _cdid_cpt_base_fail("owner", "Borrowed static baseline provider is absent.", 498)
        pfs = B.static_provider
        if ((*pfs).n_F != B.n_F | (*pfs).n_A != B.n_A)
            _cdid_cpt_base_fail("owner", "Borrowed static baseline provider dimensions changed.", 498)
        return(_cdid_cpt_base_piece_static_rows(B, *pfs, cm, rowsF))
    }
    if (B.dated_provider == NULL)
        _cdid_cpt_base_fail("owner", "Borrowed dated baseline provider is absent.", 498)
    pfd = B.dated_provider
    if ((*pfd).population == NULL | (*pfd).reference == NULL | (*pfd).n_F != B.n_F | (*pfd).n_A != B.n_A)
        _cdid_cpt_base_fail("owner", "Borrowed dated baseline provider or its objects changed.", 498)
    cdid_ref_activate(*((*pfd).reference))
    return(_cdid_cpt_base_piece_dated_rows(B, *pfd, cm, rowsF))
}

/* Replay chunk k (1..n_chunks) using the exact R chunk partition:
   - Static legacy IPW: chunks 1..n_A mapped to original F via A_to_F
   - Static REG/DR and all dated IPW/REG/DR: chunks 1..n_F. */
struct cdid_cpt_piece scalar cdid_cpt_baseline_piece_chunk(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx, real scalar chunk_k)
{
    struct cdid_cpt_base_cell_meta scalar cm
    real scalar lo, hi
    real colvector rowsF
    if (cell_idx < 1 | cell_idx > B.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_base_fail("replay", "Invalid declared cell index.", 498)
    cm = B.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_base_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    if (chunk_k < 1 | chunk_k > cm.n_chunks | chunk_k != floor(chunk_k))
        _cdid_cpt_base_fail("replay", "Invalid piece chunk index.", 498)
    lo = (chunk_k - 1) * cm.chunk_rows + 1
    if (cm.is_legacy_ipw) {
        if (B.static_provider == NULL)
            _cdid_cpt_base_fail("owner", "Borrowed static baseline provider is absent.", 498)
        hi = min((B.n_A, lo + cm.chunk_rows - 1))
        rowsF = (*B.static_provider).plan.A_to_F[(lo..hi)']
    }
    else {
        hi = min((B.n_F, lo + cm.chunk_rows - 1))
        rowsF = (lo..hi)'
    }
    return(cdid_cpt_baseline_piece_rows(B, cell_idx, rowsF))
}

/* Stream all chunks of one cell into the cdid_conditional engine accumulator. */
struct cdid_cpt_cell scalar cdid_cpt_baseline_build_cell(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar k, n_chunks
    spec = cdid_cpt_baseline_cell_spec(B, cell_idx)
    n_chunks = cdid_cpt_baseline_cell_chunks(B, cell_idx)
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, B.n_F, keep_rows, budget)
    for (k = 1; k <= n_chunks; k++) {
        pc = cdid_cpt_baseline_piece_chunk(B, cell_idx, k)
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Stream an arbitrary partition/permutation row_order of original-F positions in batches of
   piece_size into the cdid_conditional engine accumulator (for split/reorder verification). */
struct cdid_cpt_cell scalar cdid_cpt_base_build_cell_rows(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget,
    real colvector row_order, real scalar piece_size)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar st, en, n_rows, bsz
    spec = cdid_cpt_baseline_cell_spec(B, cell_idx)
    n_rows = rows(row_order)
    if (n_rows != B.n_F | cols(row_order) != 1 | missing(row_order) |
        any(row_order :!= floor(row_order)) | any(row_order :< 1) | any(row_order :> B.n_F) |
        rows(uniqrows(sort(row_order,1))) != B.n_F | piece_size < 1 | piece_size != floor(piece_size) | piece_size > B.chunk_rows)
        _cdid_cpt_base_fail("replay", "Invalid custom row_order or piece_size.", 498)
    bsz = min((piece_size, budget.max_piece_rows))
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, B.n_F, keep_rows, budget)
    for (st = 1; st <= n_rows; st = st + bsz) {
        en = min((n_rows, st + bsz - 1))
        pc = cdid_cpt_baseline_piece_rows(B, cell_idx, row_order[(st..en)'])
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Build all declared cells 1..B.n_cells from caller-supplied grids (either 1 shared grid or n_cells grids). */
struct cdid_cpt_cell colvector cdid_cpt_baseline_build_all(
    struct cdid_cpt_baseline scalar B, struct cdid_cpt_grid colvector grids,
    real colvector cmap, real scalar C, real scalar keep_rows,
    struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_cell colvector out
    struct cdid_cpt_acc colvector accs
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_base_cell_meta scalar cm
    pointer(struct cdid_baseline_provider scalar) scalar pfs
    pointer(struct cdid_db_provider scalar) scalar pfd
    transmorphic scalar cache
    real scalar i, k, jg, max_chunks, lo, hi, n, R, Dl, used
    real colvector rowsF
    if (B.status != "available" | B.n_cells < 1 | B.n_unavailable > 0)
        _cdid_cpt_base_fail("unavailable", B.reason, 459)
    if (rows(grids) != 1 & rows(grids) != B.n_cells)
        _cdid_cpt_base_fail("grid", "Supply either one common grid or one grid per declared cell.", 498)
    if (B.kind == "static") {
        if (B.static_provider == NULL)
            _cdid_cpt_base_fail("owner", "Borrowed static baseline provider is absent.", 498)
        pfs = B.static_provider
        if ((*pfs).n_F != B.n_F | (*pfs).n_A != B.n_A)
            _cdid_cpt_base_fail("owner", "Borrowed static baseline provider dimensions changed.", 498)
    }
    else {
        if (B.dated_provider == NULL)
            _cdid_cpt_base_fail("owner", "Borrowed dated baseline provider is absent.", 498)
        pfd = B.dated_provider
        if ((*pfd).population == NULL | (*pfd).reference == NULL | (*pfd).n_F != B.n_F | (*pfd).n_A != B.n_A)
            _cdid_cpt_base_fail("owner", "Borrowed dated baseline provider or its objects changed.", 498)
        cdid_ref_activate(*((*pfd).reference))
    }
    out = cdid_cpt_cell(B.n_cells, 1)
    accs = cdid_cpt_acc(B.n_cells, 1)
    max_chunks = 0
    for (i = 1; i <= B.n_cells; i++) {
        jg = (rows(grids) == 1 ? 1 : i)
        accs[i] = cdid_cpt_cell_start(cdid_cpt_baseline_cell_spec(B, i), grids[jg], cmap, C, B.n_F, keep_rows, budget)
        if (B.cells[i].n_chunks > max_chunks) max_chunks = B.cells[i].n_chunks
    }
    for (k = 1; k <= max_chunks; k++) {
        cache = asarray_create("string", 1)
        for (i = 1; i <= B.n_cells; i++) {
            cm = B.cells[i]
            if (k > cm.n_chunks) continue
            lo = (k - 1) * cm.chunk_rows + 1
            if (cm.is_legacy_ipw) {
                hi = min((B.n_A, lo + cm.chunk_rows - 1))
                rowsF = (*pfs).plan.A_to_F[(lo..hi)']
            }
            else {
                hi = min((B.n_F, lo + cm.chunk_rows - 1))
                rowsF = (lo..hi)'
            }
            n = rows(rowsF); R = cm.spec.R; Dl = cm.spec.Dl
            used = _cdid_base_cells(12 * n, 2 * R + 3 * Dl + 8, B.workspace_limit, B.workspace_entries)
            if (B.kind == "static") pc = _cdid_cpt_base_piece_static_rows(B, *pfs, cm, rowsF, cache)
            else pc = _cdid_cpt_base_piece_dated_rows(B, *pfd, cm, rowsF, cache)
            jg = (rows(grids) == 1 ? 1 : i)
            cdid_cpt_cell_add(accs[i], grids[jg], cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
        }
    }
    for (i = 1; i <= B.n_cells; i++) {
        jg = (rows(grids) == 1 ? 1 : i)
        out[i] = cdid_cpt_cell_finish(accs[i], grids[jg])
    }
    return(out)
}

/* Oriented scalar h=1 IF from the underlying provider on caller-supplied rowsF
   (chunked by B.chunk_rows when rows(rowsF) > B.chunk_rows). */
real colvector cdid_cpt_base_provider_if(
    struct cdid_cpt_baseline scalar B, real scalar cell_idx, real colvector rowsF)
{
    struct cdid_cpt_base_cell_meta scalar cm
    real scalar n, st, en, bsz
    real colvector out, sub
    if (cell_idx < 1 | cell_idx > B.n_cells | cell_idx != floor(cell_idx))
        _cdid_cpt_base_fail("replay", "Invalid declared cell index.", 498)
    cm = B.cells[cell_idx]
    if (!cm.supported)
        _cdid_cpt_base_fail("replay", "Requested cell is unavailable: " + cm.reason, 459)
    n = rows(rowsF)
    if (n < 1 | cols(rowsF) != 1 | missing(rowsF))
        _cdid_cpt_base_fail("replay", "Invalid rowsF for provider IF replay.", 498)
    bsz = B.chunk_rows
    out = J(n, 1, 0)
    for (st = 1; st <= n; st = st + bsz) {
        en = min((n, st + bsz - 1))
        sub = rowsF[(st..en)']
        if (B.kind == "static") {
            if (B.static_provider == NULL)
                _cdid_cpt_base_fail("owner", "Borrowed static baseline provider is absent.", 498)
            out[(st..en)'] = cm.sign :* cdid_baseline_if_rows(*B.static_provider, sub, cm.consumer, "frame")
        }
        else {
            if (B.dated_provider == NULL)
                _cdid_cpt_base_fail("owner", "Borrowed dated baseline provider is absent.", 498)
            out[(st..en)'] = cm.sign :* cdid_db_if_rows(*B.dated_provider, sub, cm.consumer, "frame")
        }
    }
    return(out)
}

end
