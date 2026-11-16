*! Private native ORDINARY conditional-pretrend producer and reducer, current StataNow 19.5.
*! Built-in ipw/reg/dr h-weighted moments and full original-frame influence processes, bounded
*! grid/row/draw streaming, ordinary source CvM reduction. No public parser here.
*! Pinned semantics (reference only, nothing copied): cdid 1.1.0.0 @
*! ordinary mode -> did 2.5.1 conditional_did_pretest.R, pre_process_did.R, compute.att_gt.R,
*! utility_functions.R and DRDID 1.3.0 std_ipw_did_panel.R, reg_did_panel.R, drdid_panel.R,
*! std_ipw_did_rc.R, reg_did_rc.R, drdid_rc.R (fastglm_fit.R controls).
*! Requires, unchanged: cdid_conditional.mata (cdid_cpt_budget/grid/indicator,
*! _cdid_cpt_rowsum_into, R-4.5.2 _cdid_cpt_type1_quantile), cdid_inference.mata +
*! cdid_public_rng.mata (bounded caller/fixed Rademacher stream) and MIT csdid 2.0.2
*! csdid__logit_fit/csdid__invlogit (fastglmPure binomial method-3 IRLS, R linkinv clamps).
*!
*! One declared pre cell, one distinct grid point u, full original frame of N units:
*!   r_j(u)   = h_u(j) y_j          panel: y_j = Y_time - Y_base of unit j; RCS: Y of row j
*!   T(u)     = sum_j h_u(j) y_j Gam_j                   L sufficient statistics (plain sums)
*!   J(u)     = T(u)'g                                    source ATT(g,t) on the h-weighted outcome
*!   Phi_j(u) = alpha_j r_j(u) + Bm_j T(u) for pair rows j; exactly 0 off the pair
*! alpha = (N/n1) lc and Bm = (N/n1) NL (compute.att_gt n/n1 embedding). Gam, lc, NL and g are
*! y-free and built ONCE per cell from the y-invariant nuisances (PS fit/clip/trim and Hessian,
*! weighted OR QR factors and Gram inverses, weight masses) by transcribing the DRDID
*! influence-function algebra as linear forms in T. Outcome-regression coefficients
*! bOR(u) = Kb * sum_j h_u(j) y_j sqrt(w_j) Q1_j (pivoted Householder QR, as fastglm method 0)
*! therefore change with u, and so do all their IF terms. h*IF(h=1) is never used.
*! CvM = (N/M) sum_cells sum_m J(u_m)^2 over the UNWEIGHTED empirical grid (M = N rows, exact
*! binary64 duplicates kept as multiplicities, zero sampling weights included).
*! Jb(u) = N^-1 sum_c s_bc sum_{j in c} Phi_j(u): raw uncentered cluster SUMS, one common sign
*! row per draw for every cell and grid point; CvMb = (N/M) sum sum Jb^2. Critical = type 1,
*! p = 1-ecdf(CvMb)(CvM) = mean(CvMb > CvM), no plus one. Finite all-zero J/Phi gives CvM 0,
*! zero draws, critical 0 and p 0 (no chained unavailable gate). Any NA cell refuses the family.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_ocond(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
/* Explicit old-preprocessed input. The public adapter (next section) must build it with
   old pre_process_did semantics; see the R documentation for the required metadata. */
struct cdid_oc_recipe {
    real scalar panel, intercept_only
    string scalar control, method, d1_policy
    real colvector row_unit, row_period, row_y, row_w, row_g
    real colvector key, cmap, tlist
    real matrix X, Z, cells
}
struct cdid_oc_problem {
    struct cdid_cpt_budget scalar budget
    struct cdid_cpt_grid scalar grid
    real scalar panel, N, M, R, nT, k, p, K, C, nyt, intercept_only
    real scalar d1_corrected, trim_level, identity_mode, n_atoms
    string scalar method, control, d1_policy, membership_law
    real colvector tlist, row_y, row_w, row_g, row_tp, row_unit
    real colvector unit_g, cmap, atom, inc_atom, inc_point, key
    real matrix X, cells, unit_row
}
struct cdid_oc_model {
    real scalar cell, group, time, base, status, n1, L, k, embed
    real scalar gridded, provided, n_events
    string scalar method, reason, channel
    real colvector prow, y, alpha, g, J, ev_c, ev_a, ev_v
    real matrix Gam, Bm, T, Bc
}
struct cdid_oc_mult {
    real scalar N, M, C, B, K, cursor, finished
    real colvector CvMb
    real matrix contrib
}
struct cdid_oc_result {
    real scalar N, M, K, C, B, alpha, CvM, CvMcval, CvMpval, n_distinct
    real colvector CvMb, mult
    real matrix J, contrib, cells
    string scalar status, reason, method, channel, membership_law, reducer
}
struct cdid_oc_xproc {
    real scalar N, Md, M, K, C
    real colvector mult, cmap
    real matrix J, S, jset
}
void _cdid_oc_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_oc_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
real scalar _cdid_oc_whole(real scalar x, real scalar lo)
{
    return(x<. & x==floor(x) & x>=lo & x<=2^53-1)
}
real scalar _cdid_oc_cap(struct cdid_cpt_budget scalar b)
{
    return(min((b.max_workspace_cells, floor(b.max_workspace_bytes/8))))
}
real scalar _cdid_oc_charge(real scalar a, real scalar b, real scalar cap, real scalar total)
{
    if (a>0 & b>(2^53-1)/a)
        _cdid_oc_fail("Ordinary conditional allocation product exceeds exact integer range.", 908)
    if (a*b>cap-total)
        _cdid_oc_fail("Ordinary conditional numeric storage exceeds the explicit workspace budget.", 908)
    return(total+a*b)
}
real scalar _cdid_oc_tindex(real colvector tlist, real scalar v)
{
    real colvector ix
    ix = _cdid_oc_index(tlist:==v)
    return(rows(ix)==1 ? ix[1] : 0)
}
/* Positions 1..G of x in first-appearance order of its exact values (-0 = +0); this is
   R unique()/match() order, used for test.mboot clusters and RCS identity groups. */
real colvector cdid_oc_first_positions(real colvector x)
{
    real scalar n, G
    real colvector orig, grp, firstg, ord, rnk, out
    real matrix W, info
    n = rows(x)
    if (n==0) return(J(0,1,.))
    if (cols(x)!=1 | missing(x)) _cdid_oc_fail("Labels must be a finite numeric column.", 498)
    W = sort((x:+0, (1::n)), (1,2))
    orig = W[.,2]
    if (n>1) grp = 1 \ (1:+runningsum(W[(2::n),1]:!=W[(1::(n-1)),1]))
    else grp = 1
    info = _cdid_ps_ocond(grp, 1)
    G = rows(info)
    firstg = orig[info[.,1]]
    ord = order(firstg, 1)
    rnk = J(G,1,0)
    rnk[ord] = (1::G)
    out = J(n,1,0)
    out[orig] = vec(rnk[grp])
    return(out)
}
/* Declared ordinary pre family (compute.att_gt varying base, pretreatment cells only):
   group-major, then time; cell (g, tlist[t+1], tlist[t]) whenever g > tlist[t+1]. */
real matrix cdid_oc_pre_cells(real colvector glist, real colvector tlist)
{
    real matrix out
    real scalar g, t, nT
    nT = rows(tlist)
    if (cols(glist)!=1 | cols(tlist)!=1 | nT<2 | missing(glist) | missing(tlist))
        _cdid_oc_fail("Pre cells need finite glist and at least two periods.", 498)
    if (any(tlist[(2::nT)]:<=tlist[(1::(nT-1))]))
        _cdid_oc_fail("tlist must be strictly increasing.", 498)
    out = J(0,3,.)
    for (g=1; g<=rows(glist); g++) {
        for (t=1; t<nT; t++) {
            if (glist[g]>tlist[t+1]) out = out \ (glist[g], tlist[t+1], tlist[t])
        }
    }
    return(out)
}
struct cdid_oc_problem scalar cdid_oc_create(struct cdid_oc_recipe scalar rc,
    struct cdid_cpt_budget scalar budget)
{
    struct cdid_oc_problem scalar P
    real scalar R, nT, N, j, c, tp, cap, total
    real colvector lin, ur, grp, chk
    real matrix pairs
    grp = J(0,1,.)
    pairs = J(0,2,.)
    if (rc.panel!=0 & rc.panel!=1)
        _cdid_oc_fail("panel must be 1 (old-preprocessed balanced panel) or 0 (true RCS).", 498)
    if (rc.method!="ipw" & rc.method!="reg" & rc.method!="dr")
        _cdid_oc_fail("Built-in method must be ipw, reg or dr; a custom callback needs its own route.", 498)
    if (rc.control!="nevertreated" & rc.control!="notyettreated")
        _cdid_oc_fail("control must be nevertreated or notyettreated.", 498)
    if (rc.intercept_only!=0 & rc.intercept_only!=1)
        _cdid_oc_fail("intercept_only is the explicit 0/1 source flag (xformla == ~1).", 498)
    R = rows(rc.row_y)
    if (R<1 | cols(rc.row_y)!=1 | rows(rc.row_period)!=R | cols(rc.row_period)!=1 |
        rows(rc.row_w)!=R | cols(rc.row_w)!=1 | rows(rc.row_g)!=R | cols(rc.row_g)!=1 |
        rows(rc.X)!=R)
        _cdid_oc_fail("Periods, outcomes, raw weights, cohorts and nuisance design must share the retained rows.", 498)
    if (missing(rc.row_y) | missing(rc.row_period) | missing(rc.row_w) | missing(rc.row_g) | missing(rc.X))
        _cdid_oc_fail("Retained rows must be finite (old pre_process_did complete finite screening).", 498)
    if (any(rc.row_w:<0)) _cdid_oc_fail("Raw sampling weights must be nonnegative; zero weights stay in the frame.", 498)
    if (any(rc.row_g:<0)) _cdid_oc_fail("Cohorts are 0 (never treated after source recoding) or positive.", 498)
    nT = rows(rc.tlist)
    if (nT<2 | cols(rc.tlist)!=1 | missing(rc.tlist)) _cdid_oc_fail("tlist needs at least two finite periods.", 498)
    if (any(rc.tlist[(2::nT)]:<=rc.tlist[(1::(nT-1))])) _cdid_oc_fail("tlist must be strictly increasing.", 498)
    P.row_tp = J(R,1,0)
    for (j=1; j<=nT; j++) P.row_tp = P.row_tp + j:*(rc.row_period:==rc.tlist[j])
    if (any(P.row_tp:<1)) _cdid_oc_fail("Every retained row period must belong to tlist.", 498)
    if (rc.panel) {
        if (rows(rc.row_unit)!=R | cols(rc.row_unit)!=1 | missing(rc.row_unit) |
            any(rc.row_unit:!=floor(rc.row_unit)) | any(rc.row_unit:<1))
            _cdid_oc_fail("Panel rows need integer unit positions 1..N in sorted original-ID order.", 498)
        N = max(rc.row_unit)
        if (N*nT!=R) _cdid_oc_fail("The old-preprocessed panel must be balanced: N units times T periods.", 498)
        lin = (rc.row_unit:-1):*nT :+ P.row_tp
        if (rows(uniqrows(lin))!=R) _cdid_oc_fail("Each panel unit must appear exactly once per period.", 498)
        ur = J(R,1,0)
        ur[lin] = (1::R)
        P.unit_row = rowshape(ur, N)
        P.unit_g = rc.row_g[P.unit_row[.,1]]
        if (any(rc.row_g:!=P.unit_g[rc.row_unit]))
            _cdid_oc_fail("The cohort must be time invariant within each unit.", 498)
        if (rows(rc.key)!=0) _cdid_oc_fail("Panel identity is the unit position; no RCS key is accepted.", 498)
        P.row_unit = rc.row_unit
    }
    else {
        N = R
        P.unit_row = J(0,0,.)
        P.unit_g = rc.row_g
        P.row_unit = (1::R)
    }
    if (rows(rc.Z)!=N | missing(rc.Z))
        _cdid_oc_fail("Grid design Z needs one finite row per original unit (panel first calendar rows; RCS every row).", 498)
    if (rows(rc.cmap)!=N | cols(rc.cmap)!=1 | missing(rc.cmap) | any(rc.cmap:!=floor(rc.cmap)) | any(rc.cmap:<1))
        _cdid_oc_fail("Cluster positions need one positive integer per original unit.", 498)
    chk = cdid_oc_first_positions(rc.cmap)
    if (any(chk:!=rc.cmap))
        _cdid_oc_fail("Cluster positions must be 1..C in first appearance order of the source labels (test.mboot unique()).", 498)
    if (cols(rc.cells)!=3 | rows(rc.cells)<1 | missing(rc.cells))
        _cdid_oc_fail("At least one declared pre cell (group,time,base) is required; the source refuses an empty pre family.", 498)
    for (c=1; c<=rows(rc.cells); c++) {
        tp = _cdid_oc_tindex(rc.tlist, rc.cells[c,2])
        if (tp<2) _cdid_oc_fail("Cell time must be a period after the first.", 498)
        if (rc.cells[c,3]!=rc.tlist[tp-1]) _cdid_oc_fail("Ordinary pretest base is the adjacent earlier period (varying base, anticipation 0).", 498)
        if (!(rc.cells[c,1]>rc.cells[c,2])) _cdid_oc_fail("Declared cells must be pre cells: group > time.", 498)
        if (!any(P.unit_g:==rc.cells[c,1])) _cdid_oc_fail("Every declared cohort must exist in the retained frame.", 498)
    }
    if (rows(uniqrows(rc.cells[.,(1,2)]))!=rows(rc.cells)) _cdid_oc_fail("Declared pre cells must be distinct.", 498)
    if (rc.intercept_only) {
        if (cols(rc.X)!=1 | any(rc.X:!=1)) _cdid_oc_fail("intercept_only requires the single unit nuisance column of ~1.", 498)
    }
    if (!rc.panel & rc.method=="dr") {
        if (rc.d1_policy=="" | rc.d1_policy=="corrected") P.d1_corrected = 1
        else if (rc.d1_policy=="literal") P.d1_corrected = 0
        else _cdid_oc_fail("RCS DR IF channel is corrected (default) or literal.", 498)
        P.d1_policy = (P.d1_corrected ? "corrected" : "literal")
    }
    else {
        if (rc.d1_policy!="" & rc.d1_policy!="inactive" & rc.d1_policy!="corrected" & rc.d1_policy!="literal")
            _cdid_oc_fail("Unknown D1 policy.", 498)
        P.d1_corrected = 0
        P.d1_policy = "inactive"
    }
    P.C = max(rc.cmap)
    if (N>budget.max_units | R>budget.max_units | P.C>budget.max_clusters | rows(rc.cells)>budget.max_cells |
        cols(rc.X)>budget.max_coords)
        _cdid_oc_fail("Ordinary conditional dimensions exceed explicit caller limits.", 908)
    cap = _cdid_oc_cap(budget)
    total = _cdid_oc_charge(R, cols(rc.X)+8, cap, 0)
    total = _cdid_oc_charge(N, nT+cols(rc.Z)+8, cap, total)
    P.budget = budget
    /* The engine's nonempty one-distinct-point branch indexes a scalar
       rank as a row. Preserve its bytes: use its equivalent zero-column singleton
       construction, charge the original design envelope, then retain Z metadata. */
    if (cols(rc.Z)>0) {
        if (all(rc.Z:==rc.Z[1,.])) {
            (void) _cdid_oc_charge(N, 3*cols(rc.Z)+9, cap, 0)
            P.grid = cdid_cpt_grid_create(J(N,0,.), budget)
            P.grid.p = cols(rc.Z)
            P.grid.zd = rc.Z[1,.]
        }
        else P.grid = cdid_cpt_grid_create(rc.Z, budget)
    }
    else P.grid = cdid_cpt_grid_create(rc.Z, budget)
    if (P.grid.n_distinct>budget.max_distinct_points)
        _cdid_oc_fail("Distinct grid points exceed max_distinct_points.", 908)
    P.identity_mode = 0
    P.inc_atom = J(0,1,.)
    P.inc_point = J(0,1,.)
    if (rc.panel) {
        P.key = (1::N)
        P.membership_law = "panel units: covariate values of the first period"
    }
    else if (rows(rc.key)==0) {
        P.key = (1::N)
        P.membership_law = "repeated cross-sections: each observation is its own unit"
    }
    else {
        if (rows(rc.key)!=N | cols(rc.key)!=1 | missing(rc.key))
            _cdid_oc_fail("A supplied RCS identity needs one finite value per retained row.", 498)
        if (rows(uniqrows((rc.key, P.row_tp)))!=N)
            _cdid_oc_fail("With repeated cross-sections, ivar() must be unique within each period.", 498)
        if (rows(uniqrows((rc.key, rc.row_g)))!=rows(uniqrows(rc.key)))
            _cdid_oc_fail("With repeated cross-sections, gvar() must be constant within each value of ivar().", 498)
        P.key = rc.key
        grp = cdid_oc_first_positions(rc.key)
        if (max(grp)==N) P.membership_law = "repeated cross-sections: each identifier is one observation"
        else {
            P.identity_mode = 1
            P.membership_law = "repeated cross-sections: all observations sharing an identifier"
        }
    }
    if (P.identity_mode) {
        P.atom = grp
        P.n_atoms = max(grp)
        pairs = uniqrows((grp, P.grid.map))
        P.inc_atom = pairs[.,1]
        P.inc_point = pairs[.,2]
        if (P.n_atoms>budget.max_distinct_points)
            _cdid_oc_fail("Identity groups exceed max_distinct_points.", 908)
    }
    else {
        P.atom = P.grid.map
        P.n_atoms = P.grid.n_distinct
    }
    P.panel = rc.panel
    P.N = N
    P.M = N
    P.R = R
    P.nT = nT
    P.k = cols(rc.X)
    P.p = cols(rc.Z)
    P.K = rows(rc.cells)
    P.nyt = (rc.control=="notyettreated")
    P.intercept_only = rc.intercept_only
    P.trim_level = 0.995
    P.method = rc.method
    P.control = rc.control
    P.tlist = rc.tlist
    P.row_y = rc.row_y
    P.row_w = rc.row_w
    P.row_g = rc.row_g
    P.cmap = rc.cmap
    P.X = rc.X
    P.cells = rc.cells
    return(P)
}
/* Exact 1-norm reciprocal condition number, 0 for singular/invalid input (R rcond()). */
real scalar _cdid_oc_rcond(real matrix A)
{
    real matrix Ai
    real scalar a1, i1
    if (rows(A)==0 | rows(A)!=cols(A) | missing(A)) return(0)
    a1 = max(colsum(abs(A)))
    if (a1==0 | a1>=.) return(0)
    Ai = luinv(A, 0)
    if (missing(Ai)) return(0)
    i1 = max(colsum(abs(Ai)))
    if (i1==0 | i1>=.) return(0)
    return(1/(a1*i1))
}
/* fastglmPure binomial IRLS (method 3, tol 1e-8, maxit 100) for ONE non-unit column. The
   csdid replica's one-column shortcut is the intercept closed form and is valid only for ~1. */
real colvector _cdid_oc_logit_one(real colvector d, real matrix x, real colvector w)
{
    real scalar iter, h, b, dev, dev_old
    real colvector mu, omu, eta, z, ww
    mu = (w:*d:+0.5):/(w:+1)
    omu = 1:-mu
    eta = ln(mu:/omu)
    dev_old = 1e308
    b = .
    for (iter=1; iter<=100; iter++) {
        ww = w:*mu:*omu
        z = eta + (d-mu):/(mu:*omu)
        h = quadcross(x, ww, x)
        if (!(h>0 & h<.)) return(.)
        b = quadcross(x, ww:*z)/h
        if (b>=.) return(.)
        eta = x*b
        mu = csdid__invlogit(eta)
        omu = 1:-mu
        dev = -2*quadsum(w:*(d:*ln(mu)+(1:-d):*ln(omu)))
        if (abs(dev-dev_old)/(abs(dev)+0.1)<1e-8) break
        dev_old = dev
    }
    return(b)
}
real colvector _cdid_oc_logit(real colvector d, real matrix x, real colvector w)
{
    if (cols(x)==1) {
        if (!all(x:==1)) return(_cdid_oc_logit_one(d, x, w))
    }
    return(csdid__logit_fit(d, x, w))
}
/* did utility_functions.R:87-94: UNWEIGHTED preliminary logit on every pair row (zero
   sampling weights included); intercept-only closed form outside the 1e-6 knife edge. */
real scalar _cdid_oc_overlap_fail(real matrix X, real colvector D, real scalar intercept_only)
{
    real scalar pbar
    real colvector b, ps
    if (rows(D)==0) return(1)
    if (intercept_only) {
        pbar = mean(D)
        if (abs(pbar-0.999)>1e-6) return(pbar>=0.999)
    }
    b = _cdid_oc_logit(D, X, J(rows(D),1,1))
    if (missing(b)) return(1)
    ps = csdid__invlogit(X*b)
    if (missing(ps)) return(1)
    return(max(ps)>=0.999)
}
/* did utility_functions.R:107-113: unweighted control crossprod, double precision. */
real scalar _cdid_oc_rcond_fail(real matrix X, real colvector D, real scalar intercept_only)
{
    real matrix Xc
    if (intercept_only) return(!any(D:==0))
    Xc = select(X, D:==0)
    if (rows(Xc)==0) return(1)
    return(_cdid_oc_rcond(cross(Xc, Xc))<epsilon(1))
}
/* DRDID PS block: fastglm_fit(int.cov, D, binomial, i.weights, method=3); pmin(ps, 1-1e-6);
   rcond(XtWX.ps) guard; Hessian.ps = chol2inv(chol(XtWX.ps))*n; asy.lin.rep.ps rows. */
real scalar _cdid_oc_pscore(real matrix X, real colvector D, real colvector iw,
    real colvector ps, real matrix asy)
{
    real scalar n1
    real colvector b, Wv
    real matrix XtWX, Hps
    n1 = rows(X)
    b = _cdid_oc_logit(D, X, iw)
    if (missing(b)) return(1)
    ps = csdid__invlogit(X*b)
    if (missing(ps)) return(1)
    ps = rowmin((ps, J(n1,1,1-1e-6)))
    Wv = ps:*(1:-ps):*iw
    XtWX = quadcross(X, Wv, X)
    if (_cdid_oc_rcond(XtWX)<epsilon(1)) return(2)
    Hps = cholinv(XtWX)
    if (missing(Hps)) return(2)
    asy = ((iw:*(D:-ps)):*X)*(Hps:*n1)
    return(0)
}
string scalar _cdid_oc_psreason(real scalar st)
{
    if (st==1) return("propensity score fit failed (NA coefficients)")
    return("propensity score design singular: rcond(XtWX.ps) < eps or no Cholesky factor")
}
/* Weighted OR on one block (fastglm gaussian method 0 = pivoted QR of sqrt(w)X) with the
   DRDID rcond(XpX) guard. Returns Kb (bOR = Kb * Q1'(sqrt(w) r)), Ginv = (X'WX)^-1 and the
   block's statistic coefficients Qs = sqrt(w) Q1 rows. Fewer positive weights than columns
   is rank deficient (fastglm NA coefficient -> DRDID stop). */
real scalar _cdid_oc_ols(real matrix Xs, real colvector ws, real scalar n1,
    real matrix Kb, real matrix Ginv, real matrix Qs)
{
    real scalar k, rk
    real matrix G, A, H, R1, R1inv, tau, piv
    k = cols(Xs)
    if (k<1 | rows(Xs)<1) return(1)
    if (sum(ws:>0)<k) return(1)
    G = quadcross(Xs, ws, Xs)
    if (_cdid_oc_rcond(G:/n1)<epsilon(1)) return(1)
    A = Xs:*sqrt(ws)
    H = J(0,0,.)
    tau = J(0,0,.)
    R1 = J(0,0,.)
    piv = .
    hqrdp(A, H, tau, R1, piv)
    if (rows(R1)!=k | cols(R1)!=k | missing(R1)) return(1)
    if (any(diagonal(R1):==0)) return(1)
    rk = .
    R1inv = solveupper(R1, I(k), rk, 0)
    if (missing(R1inv) | rk!=k) return(1)
    Kb = R1inv[invorder(piv), .]
    Ginv = Kb*Kb'
    Qs = hqrdq1(H, tau):*sqrt(ws)
    if (missing(Kb) | missing(Qs) | rows(Qs)!=rows(Xs) | cols(Qs)!=k) return(1)
    return(0)
}
real rowvector _cdid_oc_e(real scalar L, real scalar j)
{
    real rowvector v
    v = J(1, L, 0)
    v[j] = 1
    return(v)
}
/* k x L linear form selecting the vector statistic stored at columns j0..j0+k-1. */
real matrix _cdid_oc_eblock(real scalar L, real scalar j0, real scalar k)
{
    real matrix A
    A = J(k, L, 0)
    A[., (j0..(j0+k-1))] = I(k)
    return(A)
}
/* k x L linear form of OR coefficients bOR = Kb * (QR statistic at j0..j0+k-1). */
real matrix _cdid_oc_kblock(real scalar L, real scalar j0, real matrix Kb)
{
    real matrix A
    real scalar k
    k = rows(Kb)
    A = J(k, L, 0)
    A[., (j0..(j0+k-1))] = Kb
    return(A)
}
real scalar _cdid_oc_badmass(real rowvector S)
{
    return(any(S:<=0) | missing(S))
}
/* sgn * [w (r - F)/mw - w att/mw]: a DRDID leading ratio term (F empty means F = 0). */
void _cdid_oc_add_ratio(real colvector lc, real matrix NL, real scalar sgn,
    real colvector w, real scalar mw, real matrix F, real rowvector att)
{
    real colvector q
    q = w:/mw
    lc = lc + sgn:*q
    if (rows(F)) NL = NL - sgn:*(q:*F)
    NL = NL - sgn:*(q*att)
}
/* sgn * wo (r - Xb) (X v): an OLS asymptotic-linear-representation term times a fixed
   derivative vector v = XpX.inv * M, with the residual's bOR(u) kept as a linear form. */
void _cdid_oc_add_resid(real colvector lc, real matrix NL, real scalar sgn,
    real colvector wo, real matrix X, real colvector v, real matrix Xb)
{
    real colvector q
    q = wo:*(X*v)
    lc = lc + sgn:*q
    NL = NL - sgn:*(q:*Xb)
}
void _cdid_oc_cellfail(struct cdid_oc_model scalar m, string scalar why)
{
    m.status = 1
    m.reason = why
    m.L = 0
    m.alpha = J(0,1,.)
    m.g = J(0,1,.)
    m.Gam = J(0,0,.)
    m.Bm = J(0,0,.)
}
void _cdid_oc_setcell(struct cdid_oc_model scalar m, real matrix Gam, real colvector lc,
    real matrix NL, real colvector g)
{
    if (missing(Gam) | missing(lc) | missing(NL) | missing(g)) {
        _cdid_oc_cellfail(m, "nonfinite estimator coefficients (source ATT/IF NA or NaN)")
        return
    }
    m.L = cols(Gam)
    m.Gam = Gam
    m.g = g
    m.alpha = m.embed:*lc
    m.Bm = m.embed:*NL
}
/* DRDID std_ipw_did_panel.R:87-150 on r = h*deltaY. */
void _cdid_oc_ipw_panel(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector iw, real scalar trim_level)
{
    real scalar n1, k, L, st, S1, S0
    real colvector ps, trim, w1, w0, lc, sX0
    real matrix asy, Gam, NL, M2
    real rowvector eta1, eta0
    n1 = rows(X)
    k = cols(X)
    ps = J(0,1,.)
    asy = J(0,0,.)
    st = _cdid_oc_pscore(X, D, iw, ps, asy)
    if (st) {
        _cdid_oc_cellfail(m, _cdid_oc_psreason(st))
        return
    }
    trim = (D:==1) :| (ps:<trim_level)
    w1 = trim:*iw:*D
    w0 = trim:*iw:*ps:*(1:-D):/(1:-ps)
    S1 = quadsum(w1)
    S0 = quadsum(w0)
    if (_cdid_oc_badmass((S1, S0))) {
        _cdid_oc_cellfail(m, "zero treated or comparison weight mass (source ratio NaN)")
        return
    }
    L = 2+k
    Gam = (w1, w0, w0:*X)
    eta1 = _cdid_oc_e(L,1):/S1
    eta0 = _cdid_oc_e(L,2):/S0
    sX0 = quadcross(X, w0)
    M2 = (_cdid_oc_eblock(L,3,k) - sX0*eta0):/n1
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    _cdid_oc_add_ratio(lc, NL, 1, w1, S1/n1, J(0,0,.), eta1)
    _cdid_oc_add_ratio(lc, NL, -1, w0, S0/n1, J(0,0,.), eta0)
    NL = NL - (asy*M2):/(S0/n1)
    _cdid_oc_setcell(m, Gam, lc, NL, (eta1-eta0)')
}
/* DRDID reg_did_panel.R:104-165 on r = h*deltaY. */
void _cdid_oc_reg_panel(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector iw)
{
    real scalar n1, k, L, S1, mw1
    real colvector w1, wols, ctrl, lc, sX1
    real matrix Kb, Ginv, Qs, Gam, bOR, Xb, NL
    real rowvector et, ec
    n1 = rows(X)
    k = cols(X)
    ctrl = _cdid_oc_index(D:==0)
    Kb = J(0,0,.)
    Ginv = J(0,0,.)
    Qs = J(0,0,.)
    if (_cdid_oc_ols(X[ctrl,.], iw[ctrl], n1, Kb, Ginv, Qs)) {
        _cdid_oc_cellfail(m, "control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    w1 = iw:*D
    wols = iw:*(1:-D)
    S1 = quadsum(w1)
    if (_cdid_oc_badmass(S1)) {
        _cdid_oc_cellfail(m, "zero treated weight mass (source ratio NaN)")
        return
    }
    mw1 = S1/n1
    L = 1+k
    Gam = J(n1,L,0)
    Gam[.,1] = w1
    Gam[ctrl,(2..(1+k))] = Qs
    bOR = _cdid_oc_kblock(L, 2, Kb)
    Xb = X*bOR
    sX1 = quadcross(X, w1)
    et = _cdid_oc_e(L,1):/S1
    ec = (sX1'*bOR):/S1
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    _cdid_oc_add_ratio(lc, NL, 1, w1, mw1, J(0,0,.), et)
    NL = NL - ((w1:*Xb) - w1*ec):/mw1
    _cdid_oc_add_resid(lc, NL, -1/mw1, wols, X, (n1:*Ginv)*(sX1:/n1), Xb)
    _cdid_oc_setcell(m, Gam, lc, NL, (et-ec)')
}
/* DRDID drdid_panel.R:103-201 on r = h*deltaY. */
void _cdid_oc_dr_panel(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector iw, real scalar trim_level)
{
    real scalar n1, k, L, st, S1, S0, mw1, mw0
    real colvector ps, trim, w1, w0, wols, ctrl, lc, sX1, sX0
    real matrix asy, Kb, Ginv, Qs, Gam, bOR, Xb, NL, G0, M2, XpXinv
    real rowvector et, ec
    n1 = rows(X)
    k = cols(X)
    ps = J(0,1,.)
    asy = J(0,0,.)
    st = _cdid_oc_pscore(X, D, iw, ps, asy)
    if (st) {
        _cdid_oc_cellfail(m, _cdid_oc_psreason(st))
        return
    }
    ctrl = _cdid_oc_index(D:==0)
    Kb = J(0,0,.)
    Ginv = J(0,0,.)
    Qs = J(0,0,.)
    if (_cdid_oc_ols(X[ctrl,.], iw[ctrl], n1, Kb, Ginv, Qs)) {
        _cdid_oc_cellfail(m, "control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    trim = (D:==1) :| (ps:<trim_level)
    w1 = trim:*iw:*D
    w0 = trim:*iw:*ps:*(1:-D):/(1:-ps)
    wols = iw:*(1:-D)
    S1 = quadsum(w1)
    S0 = quadsum(w0)
    if (_cdid_oc_badmass((S1, S0))) {
        _cdid_oc_cellfail(m, "zero treated or comparison weight mass (source ratio NaN)")
        return
    }
    mw1 = S1/n1
    mw0 = S0/n1
    L = 2+2*k
    Gam = J(n1,L,0)
    Gam[.,1] = w1
    Gam[.,2] = w0
    Gam[ctrl,(3..(2+k))] = Qs
    Gam[.,((3+k)..(2+2*k))] = w0:*X
    bOR = _cdid_oc_kblock(L, 3, Kb)
    Xb = X*bOR
    sX1 = quadcross(X, w1)
    sX0 = quadcross(X, w0)
    et = (_cdid_oc_e(L,1) - sX1'*bOR):/S1
    ec = (_cdid_oc_e(L,2) - sX0'*bOR):/S0
    XpXinv = n1:*Ginv
    G0 = quadcross(X, w0, X)
    M2 = (_cdid_oc_eblock(L,3+k,k) - G0*bOR - sX0*ec):/n1
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    _cdid_oc_add_ratio(lc, NL, 1, w1, mw1, Xb, et)
    _cdid_oc_add_resid(lc, NL, -1/mw1, wols, X, XpXinv*(sX1:/n1), Xb)
    _cdid_oc_add_ratio(lc, NL, -1, w0, mw0, Xb, ec)
    NL = NL - (asy*M2):/mw0
    _cdid_oc_add_resid(lc, NL, 1/mw0, wols, X, XpXinv*(sX0:/n1), Xb)
    _cdid_oc_setcell(m, Gam, lc, NL, (et-ec)')
}
/* DRDID std_ipw_did_rc.R:81-160 on r = h*Y. */
void _cdid_oc_ipw_rc(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector post, real colvector iw, real scalar trim_level)
{
    real scalar n1, k, L, st
    real colvector ps, trim, wtp, wtq, wcp, wcq, lc
    real matrix asy, Gam, NL, M2p, M2q
    real rowvector Sm, atp, atq, acp, acq
    n1 = rows(X)
    k = cols(X)
    ps = J(0,1,.)
    asy = J(0,0,.)
    st = _cdid_oc_pscore(X, D, iw, ps, asy)
    if (st) {
        _cdid_oc_cellfail(m, _cdid_oc_psreason(st))
        return
    }
    trim = (D:==1) :| (ps:<trim_level)
    wtp = trim:*iw:*D:*(1:-post)
    wtq = trim:*iw:*D:*post
    wcp = trim:*iw:*ps:*(1:-D):*(1:-post):/(1:-ps)
    wcq = trim:*iw:*ps:*(1:-D):*post:/(1:-ps)
    Sm = (quadsum(wtp), quadsum(wtq), quadsum(wcp), quadsum(wcq))
    if (_cdid_oc_badmass(Sm)) {
        _cdid_oc_cellfail(m, "zero treated or comparison weight mass in a period (source ratio NaN)")
        return
    }
    L = 4+2*k
    Gam = (wtp, wtq, wcp, wcq, wcp:*X, wcq:*X)
    atp = _cdid_oc_e(L,1):/Sm[1]
    atq = _cdid_oc_e(L,2):/Sm[2]
    acp = _cdid_oc_e(L,3):/Sm[3]
    acq = _cdid_oc_e(L,4):/Sm[4]
    M2p = (_cdid_oc_eblock(L,5,k) - quadcross(X, wcp)*acp):/Sm[3]
    M2q = (_cdid_oc_eblock(L,5+k,k) - quadcross(X, wcq)*acq):/Sm[4]
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    _cdid_oc_add_ratio(lc, NL, 1, wtq, Sm[2]/n1, J(0,0,.), atq)
    _cdid_oc_add_ratio(lc, NL, -1, wtp, Sm[1]/n1, J(0,0,.), atp)
    _cdid_oc_add_ratio(lc, NL, -1, wcq, Sm[4]/n1, J(0,0,.), acq)
    _cdid_oc_add_ratio(lc, NL, 1, wcp, Sm[3]/n1, J(0,0,.), acp)
    NL = NL - asy*(M2q - M2p)
    _cdid_oc_setcell(m, Gam, lc, NL, ((atq-atp) - (acq-acp))')
}
/* DRDID reg_did_rc.R:86-184 on r = h*Y. */
void _cdid_oc_reg_rc(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector post, real colvector iw)
{
    real scalar n1, k, L, mc
    real colvector ipre, ipost, wtp, wtq, wc, wop, woq, lc, sXc
    real matrix Kp, Gp, Qp, Kq, Gq, Qq, Gam, bp, bq, Xbp, Xbq, NL
    real rowvector Sm, etp, etq, ecn
    n1 = rows(X)
    k = cols(X)
    ipre = _cdid_oc_index((D:==0) :& (post:==0))
    ipost = _cdid_oc_index((D:==0) :& (post:==1))
    Kp = J(0,0,.)
    Gp = J(0,0,.)
    Qp = J(0,0,.)
    Kq = J(0,0,.)
    Gq = J(0,0,.)
    Qq = J(0,0,.)
    if (_cdid_oc_ols(X[ipre,.], iw[ipre], n1, Kp, Gp, Qp)) {
        _cdid_oc_cellfail(m, "pre-treatment control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    if (_cdid_oc_ols(X[ipost,.], iw[ipost], n1, Kq, Gq, Qq)) {
        _cdid_oc_cellfail(m, "post-treatment control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    wtp = iw:*D:*(1:-post)
    wtq = iw:*D:*post
    wc = iw:*D
    wop = iw:*(1:-D):*(1:-post)
    woq = iw:*(1:-D):*post
    Sm = (quadsum(wtp), quadsum(wtq), quadsum(wc))
    if (_cdid_oc_badmass(Sm)) {
        _cdid_oc_cellfail(m, "zero treated weight mass in a period (source ratio NaN)")
        return
    }
    L = 2+2*k
    Gam = J(n1,L,0)
    Gam[.,1] = wtp
    Gam[.,2] = wtq
    Gam[ipre,(3..(2+k))] = Qp
    Gam[ipost,((3+k)..(2+2*k))] = Qq
    bp = _cdid_oc_kblock(L, 3, Kp)
    bq = _cdid_oc_kblock(L, 3+k, Kq)
    Xbp = X*bp
    Xbq = X*bq
    sXc = quadcross(X, wc)
    etp = _cdid_oc_e(L,1):/Sm[1]
    etq = _cdid_oc_e(L,2):/Sm[2]
    ecn = (sXc'*(bq - bp)):/Sm[3]
    mc = Sm[3]/n1
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    _cdid_oc_add_ratio(lc, NL, 1, wtq, Sm[2]/n1, J(0,0,.), etq)
    _cdid_oc_add_ratio(lc, NL, -1, wtp, Sm[1]/n1, J(0,0,.), etp)
    NL = NL - ((wc:*(Xbq - Xbp)) - wc*ecn):/mc
    _cdid_oc_add_resid(lc, NL, -1/mc, woq, X, (n1:*Gq)*(sXc:/n1), Xbq)
    _cdid_oc_add_resid(lc, NL, 1/mc, wop, X, (n1:*Gp)*(sXc:/n1), Xbp)
    _cdid_oc_setcell(m, Gam, lc, NL, ((etq-etp) - ecn)')
}
/* DRDID drdid_rc.R:87-346 on r = h*Y. sigma = +1 literal installed IF (pre-OR terms added);
   sigma = -1 corrected influence function (inf.treat.or/inf.cont.or use post - pre). Points equal. */
void _cdid_oc_dr_rc(struct cdid_oc_model scalar m, real matrix X, real colvector D,
    real colvector post, real colvector iw, real scalar trim_level, real scalar sigma)
{
    real scalar n1, k, L, st, o1, o2, o3, o4, o5, o6
    real colvector ps, trim, i1, i2, i3, i4, wtp, wtq, wcp, wcq, wd, wd1, wd0
    real colvector wo1, wo2, wo3, wo4, lc, sXtp, sXtq, sXcp, sXcq, sXd, sXd1, sXd0
    real colvector M1q, M1p, M3q, M3p, momq, momp
    real matrix asy, K1, G1, Q1, K2, G2, Q2, K3, G3, Q3, K4, G4, Q4, Gam
    real matrix b1, b2, b3, b4, Xb1, Xb2, Xb3, Xb4, OYC, dq, dp, NL, M2p, M2q
    real rowvector Sm, atp, atq, acp, acq, adq, ad1q, adp, ad0p
    n1 = rows(X)
    k = cols(X)
    ps = J(0,1,.)
    asy = J(0,0,.)
    st = _cdid_oc_pscore(X, D, iw, ps, asy)
    if (st) {
        _cdid_oc_cellfail(m, _cdid_oc_psreason(st))
        return
    }
    i1 = _cdid_oc_index((D:==0) :& (post:==0))
    i2 = _cdid_oc_index((D:==0) :& (post:==1))
    i3 = _cdid_oc_index((D:==1) :& (post:==0))
    i4 = _cdid_oc_index((D:==1) :& (post:==1))
    K1 = J(0,0,.)
    G1 = J(0,0,.)
    Q1 = J(0,0,.)
    K2 = J(0,0,.)
    G2 = J(0,0,.)
    Q2 = J(0,0,.)
    K3 = J(0,0,.)
    G3 = J(0,0,.)
    Q3 = J(0,0,.)
    K4 = J(0,0,.)
    G4 = J(0,0,.)
    Q4 = J(0,0,.)
    if (_cdid_oc_ols(X[i1,.], iw[i1], n1, K1, G1, Q1)) {
        _cdid_oc_cellfail(m, "pre-treatment control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    if (_cdid_oc_ols(X[i2,.], iw[i2], n1, K2, G2, Q2)) {
        _cdid_oc_cellfail(m, "post-treatment control outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    if (_cdid_oc_ols(X[i3,.], iw[i3], n1, K3, G3, Q3)) {
        _cdid_oc_cellfail(m, "pre-treatment treated outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    if (_cdid_oc_ols(X[i4,.], iw[i4], n1, K4, G4, Q4)) {
        _cdid_oc_cellfail(m, "post-treatment treated outcome regression singular or rank deficient (DRDID stop)")
        return
    }
    trim = (D:==1) :| (ps:<trim_level)
    wtp = trim:*iw:*D:*(1:-post)
    wtq = trim:*iw:*D:*post
    wcp = trim:*iw:*ps:*(1:-D):*(1:-post):/(1:-ps)
    wcq = trim:*iw:*ps:*(1:-D):*post:/(1:-ps)
    wd = trim:*iw:*D
    wd1 = trim:*iw:*D:*post
    wd0 = trim:*iw:*D:*(1:-post)
    Sm = (quadsum(wtp), quadsum(wtq), quadsum(wcp), quadsum(wcq), quadsum(wd), quadsum(wd1), quadsum(wd0))
    if (_cdid_oc_badmass(Sm)) {
        _cdid_oc_cellfail(m, "zero treated or comparison weight mass in a period (source ratio NaN)")
        return
    }
    L = 4+6*k
    o1 = 5
    o2 = 5+k
    o3 = 5+2*k
    o4 = 5+3*k
    o5 = 5+4*k
    o6 = 5+5*k
    Gam = J(n1,L,0)
    Gam[.,(1..4)] = (wtp, wtq, wcp, wcq)
    Gam[i1,(o1..(o1+k-1))] = Q1
    Gam[i2,(o2..(o2+k-1))] = Q2
    Gam[i3,(o3..(o3+k-1))] = Q3
    Gam[i4,(o4..(o4+k-1))] = Q4
    Gam[.,(o5..(o5+k-1))] = wcp:*X
    Gam[.,(o6..(o6+k-1))] = wcq:*X
    b1 = _cdid_oc_kblock(L, o1, K1)
    b2 = _cdid_oc_kblock(L, o2, K2)
    b3 = _cdid_oc_kblock(L, o3, K3)
    b4 = _cdid_oc_kblock(L, o4, K4)
    Xb1 = X*b1
    Xb2 = X*b2
    Xb3 = X*b3
    Xb4 = X*b4
    OYC = post:*Xb2 + (1:-post):*Xb1
    atp = (_cdid_oc_e(L,1) - quadcross(wtp, OYC)):/Sm[1]
    atq = (_cdid_oc_e(L,2) - quadcross(wtq, OYC)):/Sm[2]
    acp = (_cdid_oc_e(L,3) - quadcross(wcp, OYC)):/Sm[3]
    acq = (_cdid_oc_e(L,4) - quadcross(wcq, OYC)):/Sm[4]
    dq = Xb4 - Xb2
    dp = Xb3 - Xb1
    adq = quadcross(wd, dq):/Sm[5]
    ad1q = quadcross(wd1, dq):/Sm[6]
    adp = quadcross(wd, dp):/Sm[5]
    ad0p = quadcross(wd0, dp):/Sm[7]
    wo1 = iw:*(1:-D):*(1:-post)
    wo2 = iw:*(1:-D):*post
    wo3 = iw:*D:*(1:-post)
    wo4 = iw:*D:*post
    sXtp = quadcross(X, wtp)
    sXtq = quadcross(X, wtq)
    sXcp = quadcross(X, wcp)
    sXcq = quadcross(X, wcq)
    sXd = quadcross(X, wd)
    sXd1 = quadcross(X, wd1)
    sXd0 = quadcross(X, wd0)
    M1q = -sXtq:/Sm[2]
    M1p = -sXtp:/Sm[1]
    M3q = -sXcq:/Sm[4]
    M3p = -sXcp:/Sm[3]
    momq = sXd:/Sm[5] - sXd1:/Sm[6]
    momp = sXd:/Sm[5] - sXd0:/Sm[7]
    M2p = (_cdid_oc_eblock(L,o5,k) - quadcross(X, wcp, OYC) - sXcp*acp):/Sm[3]
    M2q = (_cdid_oc_eblock(L,o6,k) - quadcross(X, wcq, OYC) - sXcq*acq):/Sm[4]
    lc = J(n1,1,0)
    NL = J(n1,L,0)
    /* IF = inf.treat - inf.cont + inf.eff + inf.or (drdid_rc.R:331-346) */
    _cdid_oc_add_ratio(lc, NL, 1, wtq, Sm[2]/n1, OYC, atq)
    _cdid_oc_add_ratio(lc, NL, -1, wtp, Sm[1]/n1, OYC, atp)
    _cdid_oc_add_resid(lc, NL, 1, wo2, X, (n1:*G2)*M1q, Xb2)
    _cdid_oc_add_resid(lc, NL, sigma, wo1, X, (n1:*G1)*M1p, Xb1)
    _cdid_oc_add_ratio(lc, NL, -1, wcq, Sm[4]/n1, OYC, acq)
    _cdid_oc_add_ratio(lc, NL, 1, wcp, Sm[3]/n1, OYC, acp)
    NL = NL - asy*(M2q - M2p)
    _cdid_oc_add_resid(lc, NL, -1, wo2, X, (n1:*G2)*M3q, Xb2)
    _cdid_oc_add_resid(lc, NL, -sigma, wo1, X, (n1:*G1)*M3p, Xb1)
    NL = NL + (wd:/(Sm[5]/n1)):*dq - (wd:/(Sm[5]/n1))*adq
    NL = NL - (wd1:/(Sm[6]/n1)):*dq + (wd1:/(Sm[6]/n1))*ad1q
    NL = NL - (wd:/(Sm[5]/n1)):*dp + (wd:/(Sm[5]/n1))*adp
    NL = NL + (wd0:/(Sm[7]/n1)):*dp - (wd0:/(Sm[7]/n1))*ad0p
    _cdid_oc_add_resid(lc, NL, 1, wo4, X, (n1:*G4)*momq, Xb4)
    _cdid_oc_add_resid(lc, NL, -1, wo2, X, (n1:*G2)*momq, Xb2)
    _cdid_oc_add_resid(lc, NL, -1, wo3, X, (n1:*G3)*momp, Xb3)
    _cdid_oc_add_resid(lc, NL, 1, wo1, X, (n1:*G1)*momp, Xb1)
    _cdid_oc_setcell(m, Gam, lc, NL, ((atq-atp) - (acq-acp) + (adq-ad1q) - (adp-ad0p))')
}
/* Response-independent preparation of one declared cell: source pair rows and controls,
   guards, nuisances and the cell's linear process representation. A source-NA cell is
   returned with status 1 and its reason (never dropped); no grid or draw work happens. */
struct cdid_oc_model scalar cdid_oc_model_prepare(struct cdid_oc_problem scalar P, real scalar c)
{
    struct cdid_oc_model scalar m
    real scalar tpp, tpb, grp, mw
    real colvector Gf, Cf, D, w, post, rr, gr, sel, iw
    real matrix X
    if (!_cdid_oc_whole(c,1) | c>P.K) _cdid_oc_fail("Cell index must name a declared pre cell.", 498)
    m.cell = c
    m.group = P.cells[c,1]
    m.time = P.cells[c,2]
    m.base = P.cells[c,3]
    m.method = P.method
    m.status = 0
    m.reason = ""
    m.gridded = 0
    m.provided = 0
    m.n_events = 0
    m.k = P.k
    m.L = 0
    if (P.d1_policy=="corrected") m.channel = "corrected doubly robust influence function for repeated cross-sections"
    else if (P.d1_policy=="literal") m.channel = "influence function of the doubly robust estimator for repeated cross-sections as computed by DRDID"
    else m.channel = "influence function of the ordinary estimator"
    tpp = _cdid_oc_tindex(P.tlist, m.time)
    tpb = _cdid_oc_tindex(P.tlist, m.base)
    grp = m.group
    if (P.panel) {
        Gf = (P.unit_g:==grp)
        if (P.nyt) Cf = (P.unit_g:==0) :| ((P.unit_g:>m.time) :& (P.unit_g:!=grp))
        else Cf = (P.unit_g:==0)
        m.prow = _cdid_oc_index(Gf :| Cf)
        D = Gf[m.prow]
        m.y = P.row_y[P.unit_row[m.prow, tpp]] - P.row_y[P.unit_row[m.prow, tpb]]
        w = P.row_w[P.unit_row[m.prow, tpb]]
        X = P.X[P.unit_row[m.prow, tpb], .]
        post = J(rows(m.prow),1,0)
    }
    else {
        rr = _cdid_oc_index((P.row_tp:==tpp) :| (P.row_tp:==tpb))
        gr = P.row_g[rr]
        if (P.nyt) sel = _cdid_oc_index((gr:==grp) :| (gr:==0) :| ((gr:>m.time) :& (gr:!=grp)))
        else sel = _cdid_oc_index((gr:==grp) :| (gr:==0))
        m.prow = rr[sel]
        D = (P.row_g[m.prow]:==grp)
        post = (P.row_tp[m.prow]:==tpp)
        m.y = P.row_y[m.prow]
        w = P.row_w[m.prow]
        X = P.X[m.prow, .]
    }
    m.n1 = rows(m.prow)
    if (m.n1==0) {
        _cdid_oc_cellfail(m, "empty 2x2 comparison")
        return(m)
    }
    m.embed = P.N/m.n1
    if (!P.panel) {
        if (sum(D:*post)==0) {
            _cdid_oc_cellfail(m, "source skip: no units in group "+strofreal(grp)+" in time period "+strofreal(m.time))
            return(m)
        }
        if (sum(D:*(1:-post))==0) {
            _cdid_oc_cellfail(m, "source skip: no units in group "+strofreal(grp)+" in time period "+strofreal(m.base))
            return(m)
        }
        if (sum((1:-D):*post)==0) {
            _cdid_oc_cellfail(m, "source skip: no available control units in time period "+strofreal(m.time))
            return(m)
        }
        if (sum((1:-D):*(1:-post))==0) {
            _cdid_oc_cellfail(m, "source skip: no available control units in time period "+strofreal(m.base))
            return(m)
        }
    }
    if (P.k==0) {
        _cdid_oc_cellfail(m, "zero-column nuisance design (explicit ~0): no saved source reference; refused")
        return(m)
    }
    if (P.method=="ipw" | P.method=="dr") {
        if (_cdid_oc_overlap_fail(X, D, P.intercept_only)) {
            _cdid_oc_cellfail(m, "overlap condition violated (did overlap_check_fail, unweighted max ps >= 0.999)")
            return(m)
        }
    }
    if (P.method=="reg" | P.method=="dr") {
        if (_cdid_oc_rcond_fail(X, D, P.intercept_only)) {
            _cdid_oc_cellfail(m, "control covariate matrix singular or ill-conditioned (did rcond_check_fail)")
            return(m)
        }
    }
    mw = mean(w)
    if (!(mw>0 & mw<.)) {
        _cdid_oc_cellfail(m, "nonpositive pair mean of raw weights (DRDID i.weights/mean NaN)")
        return(m)
    }
    iw = w:/mw
    if (P.panel) {
        if (P.method=="ipw") _cdid_oc_ipw_panel(m, X, D, iw, P.trim_level)
        else if (P.method=="reg") _cdid_oc_reg_panel(m, X, D, iw)
        else _cdid_oc_dr_panel(m, X, D, iw, P.trim_level)
    }
    else {
        if (P.method=="ipw") _cdid_oc_ipw_rc(m, X, D, post, iw, P.trim_level)
        else if (P.method=="reg") _cdid_oc_reg_rc(m, X, D, post, iw)
        else _cdid_oc_dr_rc(m, X, D, post, iw, P.trim_level, (P.d1_corrected ? -1 : 1))
    }
    if (!m.status) {
        if (m.L>P.budget.max_terms) _cdid_oc_fail("Cell sufficient statistics exceed max_terms.", 908)
    }
    return(m)
}
struct cdid_oc_model colvector cdid_oc_prepare_all(struct cdid_oc_problem scalar P)
{
    struct cdid_oc_model colvector models
    struct cdid_oc_model scalar m
    real scalar c
    models = cdid_oc_model(P.K, 1)
    for (c=1; c<=P.K; c++) {
        m = cdid_oc_model_prepare(P, c)
        models[c] = m
    }
    return(models)
}
/* Complete declared family: cell, group, time, base, status (0 finite, 1 source NA). */
real matrix cdid_oc_family(struct cdid_oc_model colvector models)
{
    real matrix out
    real scalar c
    out = J(rows(models), 5, .)
    for (c=1; c<=rows(models); c++) {
        out[c,.] = (models[c].cell, models[c].group, models[c].time, models[c].base, models[c].status)
    }
    return(out)
}
string colvector cdid_oc_family_reasons(struct cdid_oc_model colvector models)
{
    string colvector out
    real scalar c
    out = J(rows(models), 1, "")
    for (c=1; c<=rows(models); c++) out[c] = models[c].reason
    return(out)
}
/* Source policy: one NA cell makes CvM and every CvMb NA and quantile() stops. The native
   test refuses the whole family before any grid, draw or RNG work; nothing is subset. */
void cdid_oc_family_check(struct cdid_oc_model colvector models)
{
    real scalar c, bad
    bad = 0
    for (c=1; c<=rows(models); c++) {
        if (models[c].status) {
            bad = bad+1
            errprintf("cdid_ordinary_conditional: cell %g (group %g, time %g, base %g) is NA: %s\n",
                models[c].cell, models[c].group, models[c].time, models[c].base, models[c].reason)
        }
    }
    if (bad) _cdid_oc_fail(strofreal(bad)+" of "+strofreal(rows(models))+" declared pre cells are NA; the source test refuses the full family (no subset, no Wald substitute).", 459)
}
/* Atom memberships at distinct grid points pts: own-row atoms are distinct grid points;
   repeated RCS identities take the union of their rows' orthant indicators. */
real matrix _cdid_oc_member(struct cdid_oc_problem scalar P, real colvector pts)
{
    real matrix Hd, Ha
    real scalar Md
    Md = P.grid.n_distinct
    if (P.grid.pH == NULL) P.grid.pH = &(J(0,0,.))
    if (rows(*P.grid.pH)!=Md | cols(*P.grid.pH)!=Md) {
        if (Md>=1 & Md<=4096 & Md*Md<=_cdid_oc_cap(P.budget)/4) {
            *P.grid.pH = cdid_cpt_indicator(P.grid.zd, P.grid.zd)
        }
    }
    if (rows(*P.grid.pH)==Md & cols(*P.grid.pH)==Md) {
        Hd = (rows(pts)==Md ? *P.grid.pH : (*P.grid.pH)[., pts])
    }
    else {
        Hd = cdid_cpt_indicator(P.grid.zd, P.grid.zd[pts,.])
    }
    if (!P.identity_mode) return(Hd)
    Ha = J(P.n_atoms, rows(pts), 0)
    Ha = _cdid_cpt_rowsum_into(Ha, Hd[P.inc_point,.], P.inc_atom)
    return(Ha:>0)
}
/* Statistics T and J at every distinct grid point in bounded grid chunks. */
void cdid_oc_model_grid(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m)
{
    real scalar Md, gc, st, en, cap, total
    real colvector pts
    real matrix Ha, CA
    if (m.status) _cdid_oc_fail("A source-NA cell has no process; the full family refuses.", 459)
    if (m.gridded) _cdid_oc_fail("Cell grid statistics are already computed.", 498)
    Md = P.grid.n_distinct
    gc = min((P.budget.grid_chunk, Md))
    cap = _cdid_oc_cap(P.budget)
    total = _cdid_oc_charge(P.n_atoms, m.L+gc, cap, 0)
    total = _cdid_oc_charge(Md+rows(P.inc_atom), gc, cap, total)
    total = _cdid_oc_charge(Md, m.L+2, cap, total)
    total = _cdid_oc_charge(m.n1, m.L, cap, total)
    CA = J(P.n_atoms, m.L, 0)
    CA = _cdid_cpt_rowsum_into(CA, m.Gam:*m.y, P.atom[m.prow])
    m.T = J(Md, m.L, 0)
    for (st=1; st<=Md; st=st+gc) {
        en = min((st+gc-1, Md))
        pts = (st::en)
        Ha = _cdid_oc_member(P, pts)
        m.T[pts,.] = quadcross(Ha, CA)
    }
    m.J = m.T*m.g
    if (missing(m.T) | missing(m.J))
        _cdid_oc_fail("Nonfinite sufficient statistics for cell "+strofreal(m.cell)+".", 504)
    m.gridded = 1
}
/* Bounded explicit process block Phi[units, pts] on the full original frame (0 off pair). */
real matrix cdid_oc_process_block(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m,
    real colvector units, real colvector pts)
{
    real colvector pos, sel, uu, jj
    real matrix out, Ha
    real scalar cap, total
    if (!m.gridded) _cdid_oc_fail("Process rows need the cell grid statistics.", 498)
    if (cols(units)!=1 | cols(pts)!=1 | missing(units) | missing(pts) | any(units:!=floor(units)) |
        any(pts:!=floor(pts)) | any(units:<1) | any(units:>P.N) | any(pts:<1) | any(pts:>P.grid.n_distinct))
        _cdid_oc_fail("Process rows are original units 1..N and columns distinct grid points.", 498)
    if (rows(units)>P.budget.max_piece_rows) _cdid_oc_fail("Process block rows exceed max_piece_rows.", 908)
    cap = _cdid_oc_cap(P.budget)
    total = _cdid_oc_charge(rows(units)+P.n_atoms+P.grid.n_distinct+rows(P.inc_atom), rows(pts)+1, cap, 0)
    total = _cdid_oc_charge(P.N, 1, cap, total)
    out = J(rows(units), rows(pts), 0)
    pos = J(P.N,1,0)
    pos[m.prow] = (1::m.n1)
    sel = _cdid_oc_index(pos[units]:>0)
    if (rows(sel)) {
        Ha = _cdid_oc_member(P, pts)
        uu = units[sel]
        jj = pos[uu]
        out[sel,.] = Ha[P.atom[uu],.]:*(m.alpha[jj]:*m.y[jj]) + m.Bm[jj,.]*(m.T[pts,.]')
    }
    return(out)
}
/* Compact multiplier provider: cluster x atom sums of the local terms (sparse events,
   sorted by atom) and raw cluster sums of the nonlocal loadings (C x L). */
void cdid_oc_provider(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m)
{
    real colvector v, key, o, skey, sval, keep, cc, aa
    real matrix info, Bc
    real scalar cap, total
    if (m.status | !m.gridded) _cdid_oc_fail("The provider needs a finite gridded cell.", 498)
    if (m.provided) _cdid_oc_fail("Cell provider state is already built.", 498)
    if (P.C>(2^53-1)/P.n_atoms) _cdid_oc_fail("Clusters times atoms exceed the exact event-key range.", 908)
    cap = _cdid_oc_cap(P.budget)
    total = _cdid_oc_charge(m.n1, 6, cap, 0)
    total = _cdid_oc_charge(P.C, m.L, cap, total)
    Bc = J(P.C, m.L, 0)
    Bc = _cdid_cpt_rowsum_into(Bc, m.Bm, P.cmap[m.prow])
    m.Bc = Bc
    v = m.alpha:*m.y
    key = (P.cmap[m.prow]:-1):*P.n_atoms :+ P.atom[m.prow]
    o = order(key, 1)
    key = key[o]
    v = v[o]
    info = _cdid_ps_ocond(key, 1)
    sval = panelsum(v, info)
    skey = key[info[.,1]]
    keep = _cdid_oc_index(sval:!=0)
    if (rows(keep)) {
        skey = skey[keep]
        sval = sval[keep]
        cc = floor((skey:-1):/P.n_atoms):+1
        aa = skey - (cc:-1):*P.n_atoms
        o = order((aa, cc), (1,2))
        m.ev_c = cc[o]
        m.ev_a = aa[o]
        m.ev_v = sval[o]
    }
    else {
        m.ev_c = J(0,1,.)
        m.ev_a = J(0,1,.)
        m.ev_v = J(0,1,.)
    }
    m.n_events = rows(m.ev_v)
    if (m.n_events>P.budget.max_events) _cdid_oc_fail("Cluster-by-atom events exceed max_events.", 908)
    if (missing(m.Bc) | missing(m.ev_v)) _cdid_oc_fail("Nonfinite provider state for cell "+strofreal(m.cell)+".", 504)
    m.provided = 1
}
/* Diagnostic raw cluster sums S[c,u] = sum_{j in c} Phi_j(u) (C x rows(pts)). */
real matrix cdid_oc_cluster_sums(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m,
    real colvector pts)
{
    real matrix Ha, S
    if (!m.provided) _cdid_oc_fail("Cluster sums need the cell provider.", 498)
    Ha = _cdid_oc_member(P, pts)
    S = J(P.C, rows(pts), 0)
    if (m.n_events) S = _cdid_cpt_rowsum_into(S, Ha[m.ev_a,.]:*m.ev_v, m.ev_c)
    return(S + m.Bc*(m.T[pts,.]'))
}
/* Diagnostic draw process Jb = E S / N for supplied signs E (b x C). */
real matrix cdid_oc_draw_process(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m,
    real matrix E, real colvector pts)
{
    if (cols(E)!=P.C | missing(E) | any((E:!=1):&(E:!=-1))) _cdid_oc_fail("Signs must be b x C values -1 or +1.", 498)
    return((E*cdid_oc_cluster_sums(P, m, pts)):/P.N)
}
real colvector cdid_oc_unit_membership(struct cdid_oc_problem scalar P, real scalar d)
{
    real matrix Ha
    if (!_cdid_oc_whole(d,1) | d>P.grid.n_distinct) _cdid_oc_fail("Unknown distinct grid point.", 498)
    Ha = _cdid_oc_member(P, d)
    return(Ha[P.atom, 1])
}
/* Source h-weighted outcome on every retained row at distinct grid point d. */
real colvector cdid_oc_weighted_outcome(struct cdid_oc_problem scalar P, real scalar d)
{
    real colvector h
    h = cdid_oc_unit_membership(P, d)
    return(P.row_y:*h[P.row_unit])
}
/* CvM = (N/M) sum_cells sum_m J^2 = n*sum(apply(J^2,2,mean)) over all M grid rows. */
real scalar cdid_oc_statistic(struct cdid_oc_problem scalar P, struct cdid_oc_model colvector models)
{
    real scalar c, out
    if (rows(models)!=P.K) _cdid_oc_fail("The statistic needs the complete declared family.", 498)
    out = 0
    for (c=1; c<=P.K; c++) {
        if (models[c].status | !models[c].gridded) _cdid_oc_fail("Every declared cell must be finite and gridded.", 459)
        out = out + quadsum(P.grid.mult:*(models[c].J:^2))
    }
    return((P.N/P.M)*out)
}
/* CvMb contribution of one sign batch E (b x C) for one cell, bounded grid/event chunks. */
real colvector _cdid_oc_batch_cell(struct cdid_oc_problem scalar P, struct cdid_oc_model scalar m,
    real matrix E)
{
    real scalar b, ne, st, en, Md, gc, ec
    real colvector out, ii, dsub, o, dd, pts
    real matrix Q, Em, Emt, Ssum, info, EB, Ha, Jb
    b = rows(E)
    Md = P.grid.n_distinct
    gc = min((P.budget.grid_chunk, Md))
    ec = max((1, P.budget.event_chunk))
    out = J(b,1,0)
    Q = J(b, P.n_atoms, 0)
    ne = m.n_events
    for (st=1; st<=ne; st=st+ec) {
        en = min((ne, st+ec-1))
        ii = (st::en)
        Em = E[., m.ev_c[ii]]:*(m.ev_v[ii]')
        dsub = m.ev_a[ii]
        o = order(dsub, 1)
        dd = dsub[o]
        info = _cdid_ps_ocond(dd, 1)
        Emt = Em'
        Ssum = panelsum(Emt[o,.], info)
        dd = dd[info[.,1]]
        Q[., dd] = Q[., dd] + Ssum'
    }
    EB = E*m.Bc
    for (st=1; st<=Md; st=st+gc) {
        en = min((st+gc-1, Md))
        pts = (st::en)
        Ha = _cdid_oc_member(P, pts)
        Jb = (Q*Ha + EB*(m.T[pts,.]')):/P.N
        out = out + (P.N/P.M):*((Jb:*Jb)*P.grid.mult[pts])
    }
    return(out)
}
struct cdid_oc_mult scalar cdid_oc_mult_start(struct cdid_oc_problem scalar P,
    struct cdid_oc_model colvector models, real scalar B)
{
    struct cdid_oc_mult scalar S
    real scalar c, br, gc, ec, cap, total, one, worst
    if (!_cdid_oc_whole(B,1) | B>2^31-1) _cdid_oc_fail("The ordinary test always draws B>=1 multipliers.", 498)
    if (rows(models)!=P.K | cols(models)!=1) _cdid_oc_fail("The complete declared pre family is required (no subset).", 498)
    cdid_oc_family_check(models)
    for (c=1; c<=P.K; c++) {
        if (models[c].cell!=c) _cdid_oc_fail("Cells must follow the declared family order.", 498)
        if (!models[c].gridded | !models[c].provided) _cdid_oc_fail("Every cell needs grid statistics and a provider before drawing.", 498)
    }
    if (B>P.budget.max_draws) _cdid_oc_fail("Draws exceed max_draws.", 908)
    br = min((P.budget.max_sign_batch_rows, floor(P.budget.max_sign_batch_cells/P.C), B))
    if (br<1) _cdid_oc_fail("The sign batch budget admits no complete cluster row.", 908)
    gc = min((P.budget.grid_chunk, P.grid.n_distinct))
    cap = _cdid_oc_cap(P.budget)
    total = _cdid_oc_charge(B, P.K+1, cap, 0)
    worst = 0
    for (c=1; c<=P.K; c++) {
        ec = max((1, min((models[c].n_events, P.budget.event_chunk))))
        one = br*(P.C+P.n_atoms+models[c].L+2*gc+1) + (P.n_atoms+P.grid.n_distinct+rows(P.inc_atom))*gc + ec*(2*br+3)
        worst = max((worst, one))
    }
    total = _cdid_oc_charge(worst, 1, cap, total)
    S.N = P.N
    S.M = P.M
    S.C = P.C
    S.B = B
    S.K = P.K
    S.cursor = 0
    S.finished = 0
    S.CvMb = J(B,1,.)
    S.contrib = J(B,P.K,.)
    return(S)
}
void cdid_oc_mult_add(struct cdid_oc_problem scalar P, struct cdid_oc_model colvector models,
    struct cdid_oc_mult scalar S, struct cdid_public_sign_batch scalar E)
{
    real scalar b, c
    real colvector expected
    b = rows(E.signs)
    if (S.finished | b<1 | cols(E.signs)!=S.C | rows(E.draw_positions)!=b | cols(E.draw_positions)!=1 |
        b>S.B-S.cursor | cols(E.cluster_positions)!=S.C | rows(E.cluster_positions)!=1 | rows(models)!=S.K)
        _cdid_oc_fail("Invalid or completed ordinary multiplier sign stream.", 498)
    if (b>P.budget.max_sign_batch_rows | b*S.C>P.budget.max_sign_batch_cells)
        _cdid_oc_fail("Sign batch exceeds the explicit row/cell budget.", 908)
    expected = ((S.cursor+1)::(S.cursor+b))
    if (any(E.draw_positions:!=expected) | any(E.cluster_positions:!=(1..S.C)))
        _cdid_oc_fail("Signs must keep declared draw order and first-appearance cluster order.", 498)
    if (missing(E.signs) | any((E.signs:!=1):&(E.signs:!=-1)))
        _cdid_oc_fail("Multiplier signs must be exactly -1 or +1.", 498)
    for (c=1; c<=S.K; c++) S.contrib[expected, c] = _cdid_oc_batch_cell(P, models[c], E.signs)
    S.CvMb[expected] = quadrowsum(S.contrib[expected,.])
    if (missing(S.CvMb[expected])) _cdid_oc_fail("Nonfinite ordinary multiplier draw.", 504)
    S.cursor = S.cursor+b
}
/* did test.mboot/MP.TEST reduction on a complete draw vector. No zero-process gate. */
struct cdid_oc_result scalar cdid_oc_reduce(real scalar CvM, real colvector CvMb, real scalar alpha)
{
    struct cdid_oc_result scalar r
    if (!(alpha>0 & alpha<1)) _cdid_oc_fail("alpha must lie strictly inside (0,1).", 498)
    if (rows(CvMb)<1 | cols(CvMb)!=1) _cdid_oc_fail("The ordinary test always has B>=1 draws.", 498)
    if (CvM>=. | missing(CvMb))
        _cdid_oc_fail("Nonfinite statistic or draw: the source type-1 quantile refuses (no subset or Wald fallback).", 504)
    r.CvM = CvM
    r.CvMb = CvMb
    r.B = rows(CvMb)
    r.alpha = alpha
    r.CvMcval = _cdid_cpt_type1_quantile(CvMb, 1-alpha)
    r.CvMpval = 1 - quadsum(CvMb:<=CvM)/rows(CvMb)
    r.status = "returned"
    r.reason = ""
    r.reducer = "type-1 critical; p = 1-ecdf(CvMb)(CvM) = mean(CvMb > CvM), strict, no plus one; finite zero kept"
    return(r)
}
struct cdid_oc_result scalar cdid_oc_finish(struct cdid_oc_problem scalar P,
    struct cdid_oc_model colvector models, struct cdid_oc_mult scalar S, real scalar alpha)
{
    struct cdid_oc_result scalar r
    real scalar c
    if (S.finished | S.cursor!=S.B | rows(models)!=S.K)
        _cdid_oc_fail("Every declared draw must be consumed exactly once before finalization.", 498)
    r = cdid_oc_reduce(cdid_oc_statistic(P, models), S.CvMb, alpha)
    r.N = P.N
    r.M = P.M
    r.K = P.K
    r.C = P.C
    r.n_distinct = P.grid.n_distinct
    r.mult = P.grid.mult
    r.cells = P.cells
    r.contrib = S.contrib
    r.J = J(P.grid.n_distinct, P.K, .)
    for (c=1; c<=P.K; c++) r.J[.,c] = models[c].J
    r.method = P.method
    r.channel = models[1].channel
    r.membership_law = P.membership_law
    S.finished = 1
    return(r)
}
/* Full built-in route: prepare, refuse a NA family before touching RNG, grid, provider,
   bounded caller/fixed sign stream, ordinary reducer. sign_mode "caller" consumes the
   caller's Stata RNG in documented batches; "fixed" consumes none. */
struct cdid_oc_result scalar cdid_oc_run(struct cdid_oc_problem scalar P, real scalar B,
    string scalar sign_mode, real matrix fixed, real scalar alpha)
{
    struct cdid_oc_model colvector models
    struct cdid_oc_mult scalar S
    struct cdid_infer_budget scalar ib
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar E
    real scalar c
    if (!(alpha>0 & alpha<1)) _cdid_oc_fail("alpha must lie strictly inside (0,1).", 498)
    if (!_cdid_oc_whole(B,1)) _cdid_oc_fail("The ordinary test always draws B>=1 multipliers.", 498)
    models = cdid_oc_prepare_all(P)
    cdid_oc_family_check(models)
    for (c=1; c<=P.K; c++) {
        cdid_oc_model_grid(P, models[c])
        cdid_oc_provider(P, models[c])
    }
    S = cdid_oc_mult_start(P, models, B)
    ib = cdid_infer_budget_create(P.budget.max_units, P.budget.max_clusters, P.budget.max_cells,
        P.budget.max_draws, P.budget.max_sign_batch_rows, P.budget.max_sign_batch_cells,
        P.budget.max_sign_batch_rows, P.budget.max_sign_batch_cells,
        P.budget.max_workspace_cells, P.budget.max_workspace_bytes)
    stream = cdid_public_rng_start(B, P.C, ib, sign_mode, fixed)
    while (stream.cursor<stream.B) {
        E = cdid_public_rng_next(stream)
        cdid_oc_mult_add(P, models, S, E)
    }
    cdid_public_rng_complete(stream)
    return(cdid_oc_finish(P, models, S, alpha))
}
/* Bounded explicit-process ingestion (supplied Phi blocks; used for the source helper-policy
   arrays and as the documented consumer for a future native callback route). Storage is
   C x Md x K raw cluster sums and is refused beyond the explicit workspace budget. */
struct cdid_oc_xproc scalar cdid_oc_xp_start(real scalar N, real colvector mult, real scalar K,
    real colvector cmap, struct cdid_cpt_budget scalar budget)
{
    struct cdid_oc_xproc scalar X
    real scalar Md, cap, total, C
    real colvector chk
    Md = rows(mult)
    if (!_cdid_oc_whole(N,1) | !_cdid_oc_whole(K,1) | Md<1 | cols(mult)!=1 | missing(mult) |
        any(mult:!=floor(mult)) | any(mult:<1))
        _cdid_oc_fail("Explicit process needs N>=1, K>=1 and positive integer grid multiplicities.", 498)
    if (rows(cmap)!=N | cols(cmap)!=1 | missing(cmap) | any(cmap:!=floor(cmap)) | any(cmap:<1))
        _cdid_oc_fail("Cluster positions need one positive integer per original unit.", 498)
    chk = cdid_oc_first_positions(cmap)
    if (any(chk:!=cmap)) _cdid_oc_fail("Cluster positions must be in first-appearance order.", 498)
    C = max(cmap)
    if (N>budget.max_units | C>budget.max_clusters | K>budget.max_cells | Md>budget.max_distinct_points)
        _cdid_oc_fail("Explicit process dimensions exceed explicit caller limits.", 908)
    cap = _cdid_oc_cap(budget)
    total = _cdid_oc_charge(C, Md*K, cap, 0)
    total = _cdid_oc_charge(Md, 2*K+1, cap, total)
    X.N = N
    X.Md = Md
    X.M = quadsum(mult)
    X.K = K
    X.C = C
    X.mult = mult
    X.cmap = cmap
    X.J = J(Md, K, .)
    X.jset = J(Md, K, 0)
    X.S = J(C, Md*K, 0)
    return(X)
}
void cdid_oc_xp_add(struct cdid_oc_xproc scalar X, real scalar cell, real colvector units,
    real colvector pts, real matrix Phi)
{
    real rowvector blk
    real matrix Sb
    if (!_cdid_oc_whole(cell,1) | cell>X.K | cols(units)!=1 | cols(pts)!=1 | missing(units) | missing(pts) |
        any(units:!=floor(units)) | any(pts:!=floor(pts)) | any(units:<1) | any(units:>X.N) |
        any(pts:<1) | any(pts:>X.Md) | rows(Phi)!=rows(units) | cols(Phi)!=rows(pts))
        _cdid_oc_fail("Explicit process block must be units x distinct points of a declared cell.", 498)
    if (rows(uniqrows(pts))!=rows(pts)) _cdid_oc_fail("Explicit process columns must be distinct.", 498)
    if (missing(Phi))
        _cdid_oc_fail("Nonfinite explicit process value: the source quantile refuses NA (OC21); nothing is dropped.", 504)
    blk = ((cell-1)*X.Md):+pts'
    Sb = X.S[., blk]
    Sb = _cdid_cpt_rowsum_into(Sb, Phi, X.cmap[units])
    X.S[., blk] = Sb
}
void cdid_oc_xp_set_j(struct cdid_oc_xproc scalar X, real scalar cell, real colvector pts,
    real colvector Jv)
{
    if (!_cdid_oc_whole(cell,1) | cell>X.K | cols(pts)!=1 | rows(Jv)!=rows(pts) | cols(Jv)!=1 |
        missing(pts) | any(pts:!=floor(pts)) | any(pts:<1) | any(pts:>X.Md))
        _cdid_oc_fail("Observed J must be supplied by declared cell and distinct point.", 498)
    if (missing(Jv)) _cdid_oc_fail("Nonfinite observed J: the source statistic is NA and the test refuses.", 504)
    X.J[pts, cell] = Jv
    X.jset[pts, cell] = J(rows(pts),1,1)
}
real scalar cdid_oc_xp_statistic(struct cdid_oc_xproc scalar X)
{
    if (any(X.jset:==0)) _cdid_oc_fail("Every declared cell and grid point needs its observed J.", 498)
    return((X.N/X.M)*quadsum(X.mult:*quadrowsum(X.J:^2)))
}
real matrix _cdid_oc_xp_contrib(struct cdid_oc_xproc scalar X, real matrix E)
{
    real matrix out, Jb
    real scalar c
    out = J(rows(E), X.K, 0)
    for (c=1; c<=X.K; c++) {
        Jb = (E*X.S[., ((((c-1)*X.Md)+1)..(c*X.Md))]):/X.N
        out[.,c] = (X.N/X.M):*((Jb:*Jb)*X.mult)
    }
    return(out)
}
struct cdid_oc_mult scalar cdid_oc_xp_mult_start(struct cdid_oc_xproc scalar X, real scalar B,
    struct cdid_cpt_budget scalar budget)
{
    struct cdid_oc_mult scalar S
    real scalar cap, total
    if (!_cdid_oc_whole(B,1) | B>2^31-1) _cdid_oc_fail("The ordinary test always draws B>=1 multipliers.", 498)
    if (B>budget.max_draws) _cdid_oc_fail("Draws exceed max_draws.", 908)
    cap = _cdid_oc_cap(budget)
    total = _cdid_oc_charge(B, X.K+1, cap, 0)
    S.N = X.N
    S.M = X.M
    S.C = X.C
    S.B = B
    S.K = X.K
    S.cursor = 0
    S.finished = 0
    S.CvMb = J(B,1,.)
    S.contrib = J(B,X.K,.)
    return(S)
}
void cdid_oc_xp_mult_add(struct cdid_oc_xproc scalar X, struct cdid_oc_mult scalar S,
    struct cdid_public_sign_batch scalar E)
{
    real scalar b
    real colvector expected
    b = rows(E.signs)
    if (S.finished | b<1 | cols(E.signs)!=S.C | rows(E.draw_positions)!=b | b>S.B-S.cursor |
        cols(E.cluster_positions)!=S.C)
        _cdid_oc_fail("Invalid or completed explicit multiplier sign stream.", 498)
    expected = ((S.cursor+1)::(S.cursor+b))
    if (any(E.draw_positions:!=expected) | any(E.cluster_positions:!=(1..S.C)))
        _cdid_oc_fail("Signs must keep declared draw order and first-appearance cluster order.", 498)
    if (missing(E.signs) | any((E.signs:!=1):&(E.signs:!=-1)))
        _cdid_oc_fail("Multiplier signs must be exactly -1 or +1.", 498)
    S.contrib[expected,.] = _cdid_oc_xp_contrib(X, E.signs)
    S.CvMb[expected] = quadrowsum(S.contrib[expected,.])
    S.cursor = S.cursor+b
}
struct cdid_oc_result scalar cdid_oc_xp_finish(struct cdid_oc_xproc scalar X,
    struct cdid_oc_mult scalar S, real scalar alpha)
{
    struct cdid_oc_result scalar r
    if (S.finished | S.cursor!=S.B) _cdid_oc_fail("Every declared draw must be consumed exactly once.", 498)
    r = cdid_oc_reduce(cdid_oc_xp_statistic(X), S.CvMb, alpha)
    r.N = X.N
    r.M = X.M
    r.K = X.K
    r.C = X.C
    r.n_distinct = X.Md
    r.mult = X.mult
    r.J = X.J
    r.contrib = S.contrib
    r.method = "explicit"
    r.channel = "supplied explicit process"
    r.membership_law = "supplied"
    S.finished = 1
    return(r)
}
end
