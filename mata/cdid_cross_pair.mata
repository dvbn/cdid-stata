*! Cross-pair chained DiD (crosspair): within-unit change moments of target and clean control cohorts.
version 19.5
mata:

/* Links of the cross-pair catalogue. role 1 target, 2 control; ci cohort index;
   mid observation-mask group within the cohort; gi target index. */
struct cdid_xp_cat {
    real colvector coh, cal, tg, tref, ucoh
    real colvector g, gi, ci, cs, ct, ref, role, mid, oi
    real colvector kept, dropped
    real matrix G
    pointer(real colvector) matrix obs
}

/* One reported fit. role: 2 outcome, 3 treatment, 4 imputation. links: catalogue links that
   read the fit; est, se: coefficients and standard errors; cell: local cell number of an
   intercept; term, crole: term names and roles. */
struct _cdid_xpn_rec {
    real scalar role, nrows, wsum, nclus
    real colvector links, est, se, cell
    string colvector term, crole
    string scalar label, eqn, cohorts, kept, aliased, skey
    string colvector pairs
    real scalar g, k, mask
}

struct _cdid_xpn_mat {
    real matrix v
    real colvector cols
}

struct _cdid_xp_pairs {
    real colvector cs, ct, mid
}

struct _cdid_xp_vec {
    real colvector v
}

/* Observation-mask groups of one cohort: every date pair observed on at least one unit,
   ordered by (s, t); pairs with the same set of observed units share a group. */
void _cdid_xp_masks(real matrix O, real colvector units, real colvector cs, real colvector ct,
    real colvector mid, pointer(real colvector) rowvector obs)
{
    real scalar T, a, b, n
    real colvector m
    transmorphic A
    string scalar key
    T = cols(O)
    cs = J(0, 1, .)
    ct = J(0, 1, .)
    mid = J(0, 1, .)
    obs = J(1, 0, NULL)
    n = 0
    A = asarray_create()
    for (a = 1; a < T; a++) {
        for (b = a + 1; b <= T; b++) {
            m = O[units, a] :& O[units, b]
            if (!any(m)) continue
            key = invtokens(strofreal(m'), "")
            if (!asarray_contains(A, key)) {
                n++
                asarray(A, key, n)
                obs = obs, &(units[selectindex(m)])
            }
            cs = cs \ a
            ct = ct \ b
            mid = mid \ asarray(A, key)
        }
    }
}

/* Catalogue: for each target cohort g with a reference date, its own pairs, then the pairs
   of each clean control cohort (k > g with not-yet-treated controls, never-treated only
   otherwise) whose later date is untreated (t + anticipation < k). With restricted
   reference eligibility contexts (exA: analysis units excluded in each context, tctx: the
   context of each sorted target cohort), the pairs and observation masks of a restricted
   target and of its controls use only the units eligible in its context; oi gives the row
   of obs that holds the observed units of a link. */
struct cdid_xp_cat scalar cdid_xp_catalog(real matrix O, real colvector w, real colvector ucoh,
    real colvector cal, real colvector targets, real scalar delta, string scalar control_group,
    | real matrix exA, real colvector tctx)
{
    struct cdid_xp_cat scalar C
    real scalar K, k, j, g, gk, r, i, nk, restricted, cx, slot
    real scalar q
    real colvector units, cs, ct, mid, ctrl, keep, sl
    real matrix slotof
    pointer(real colvector) rowvector ob
    pointer(real colvector) matrix OBS
    struct _cdid_xp_pairs colvector PP
    real colvector np
    C.cal = cal
    C.ucoh = ucoh
    C.coh = uniqrows(ucoh)
    K = rows(C.coh)
    C.tg = sort(uniqrows(targets), 1)
    restricted = 0
    if (args() >= 9) restricted = (rows(tctx) == rows(C.tg) & cols(exA) > 0 & rows(exA) == rows(w))
    /* Rows K*cx+k hold the pairs of cohort k among the units eligible in context cx. */
    slotof = J(K, (restricted ? cols(exA) : 0), 0)
    nk = K * (restricted ? 1 + cols(exA) : 1)
    np = J(nk, 1, 0)
    PP = _cdid_xp_pairs(nk, 1)
    OBS = J(nk, 1, NULL)
    for (k = 1; k <= K; k++) {
        units = selectindex((ucoh :== C.coh[k]) :& (w :> 0))
        if (rows(units) == 0) continue
        _cdid_xp_masks(O, units, cs, ct, mid, ob)
        np[k] = rows(cs)
        PP[k].cs = cs
        PP[k].ct = ct
        PP[k].mid = mid
        if (cols(ob) > cols(OBS)) OBS = OBS, J(rows(OBS), cols(ob) - cols(OBS), NULL)
        for (q = 1; q <= cols(ob); q++) OBS[k, q] = ob[q]
    }
    C.tref = J(rows(C.tg), 1, .)
    C.g = J(0, 1, .)
    C.gi = J(0, 1, .)
    C.ci = J(0, 1, .)
    C.cs = J(0, 1, .)
    C.ct = J(0, 1, .)
    C.ref = J(0, 1, .)
    C.role = J(0, 1, .)
    C.mid = J(0, 1, .)
    C.oi = J(0, 1, .)
    for (j = 1; j <= rows(C.tg); j++) {
        g = C.tg[j]
        r = select(cal, (cal :+ delta) :< g)
        if (rows(r) == 0) continue
        C.tref[j] = max(r)
        gk = selectindex(C.coh :== g)
        if (rows(gk) != 1) continue
        gk = gk[1]
        cx = (restricted ? tctx[j] : 0)
        if (cx > 0) {
            for (i = 1; i <= K; i++) {
                if (slotof[i, cx] > 0) continue
                units = selectindex((ucoh :== C.coh[i]) :& (w :> 0) :& (exA[, cx] :== 0))
                slot = K * cx + i
                slotof[i, cx] = slot
                if (rows(units) == 0) continue
                _cdid_xp_masks(O, units, cs, ct, mid, ob)
                np[slot] = rows(cs)
                PP[slot].cs = cs
                PP[slot].ct = ct
                PP[slot].mid = mid
                if (cols(ob) > cols(OBS)) OBS = OBS, J(rows(OBS), cols(ob) - cols(OBS), NULL)
                for (q = 1; q <= cols(ob); q++) OBS[slot, q] = ob[q]
            }
            sl = slotof[, cx]
        }
        else sl = (1::K)
        if (np[sl[gk]] == 0) continue
        C.g = C.g \ J(np[sl[gk]], 1, g)
        C.gi = C.gi \ J(np[sl[gk]], 1, j)
        C.ci = C.ci \ J(np[sl[gk]], 1, gk)
        C.oi = C.oi \ J(np[sl[gk]], 1, sl[gk])
        C.cs = C.cs \ PP[sl[gk]].cs
        C.ct = C.ct \ PP[sl[gk]].ct
        C.ref = C.ref \ J(np[sl[gk]], 1, C.tref[j])
        C.role = C.role \ J(np[sl[gk]], 1, 1)
        C.mid = C.mid \ PP[sl[gk]].mid
        if (control_group == "nevertreated") ctrl = selectindex(C.coh :>= .)
        else ctrl = selectindex(C.coh :> g)
        for (i = 1; i <= rows(ctrl); i++) {
            k = ctrl[i]
            if (np[sl[k]] == 0) continue
            if (C.coh[k] >= .) keep = J(np[sl[k]], 1, 1)
            else keep = (cal[PP[sl[k]].ct] :+ delta) :< C.coh[k]
            nk = sum(keep)
            if (nk == 0) continue
            C.g = C.g \ J(nk, 1, g)
            C.gi = C.gi \ J(nk, 1, j)
            C.ci = C.ci \ J(nk, 1, k)
            C.oi = C.oi \ J(nk, 1, sl[k])
            C.cs = C.cs \ select(PP[sl[k]].cs, keep)
            C.ct = C.ct \ select(PP[sl[k]].ct, keep)
            C.ref = C.ref \ J(nk, 1, C.tref[j])
            C.role = C.role \ J(nk, 1, 2)
            C.mid = C.mid \ select(PP[sl[k]].mid, keep)
        }
    }
    C.obs = OBS
    return(C)
}

/* Spanning-forest reduction within each (target, cohort, mask) group: edges in order of
   (t - s, s, t); a dropped edge expands along its tree path; G = I + E'E per group. */
void cdid_xp_forest(struct cdid_xp_cat scalar C)
{
    real scalar M, gr, ng, e, m, nv, u, v, ru, rv, k, d, dd, x, y
    real colvector key, ord, idx, us, ut, dates, par, istree, tpos, dpos, depth, pv, ple, psg, queue
    real colvector keptall, Gpos, posmap, pp
    real matrix gk, E, EE, adj
    pointer(real colvector) colvector bpos
    pointer(real matrix) colvector bmat
    real scalar nbk, i1, i2
    M = rows(C.g)
    gk = (C.gi, C.ci, C.mid)
    key = J(M, 1, 0)
    ord = order((gk, (1::M)), (1, 2, 3, 4))
    ng = 0
    for (e = 1; e <= M; e++) {
        if (e == 1) ng = 1
        else if (any(gk[ord[e], ] :!= gk[ord[e - 1], ])) ng++
        key[ord[e]] = ng
    }
    keptall = J(0, 1, .)
    bpos = J(0, 1, NULL)
    bmat = J(0, 1, NULL)
    for (gr = 1; gr <= ng; gr++) {
        idx = selectindex(key :== gr)
        m = rows(idx)
        if (m == 1) {
            keptall = keptall \ idx
            continue
        }
        idx = idx[order((C.cal[C.ct[idx]] - C.cal[C.cs[idx]], C.cal[C.cs[idx]], C.cal[C.ct[idx]]), (1, 2, 3))]
        dates = uniqrows(C.cs[idx] \ C.ct[idx])
        nv = rows(dates)
        us = J(m, 1, .)
        ut = J(m, 1, .)
        for (e = 1; e <= m; e++) {
            us[e] = selectindex(dates :== C.cs[idx[e]])
            ut[e] = selectindex(dates :== C.ct[idx[e]])
        }
        par = (1::nv)
        istree = J(m, 1, 0)
        for (e = 1; e <= m; e++) {
            ru = us[e]
            while (par[ru] != ru) {
                par[ru] = par[par[ru]]
                ru = par[ru]
            }
            rv = ut[e]
            while (par[rv] != rv) {
                par[rv] = par[par[rv]]
                rv = par[rv]
            }
            if (ru != rv) {
                par[rv] = ru
                istree[e] = 1
            }
        }
        tpos = selectindex(istree)
        dpos = selectindex(!istree)
        keptall = keptall \ idx[tpos]
        if (rows(dpos) == 0) continue
        k = rows(tpos)
        /* tree adjacency: (neighbour, edge, sign) */
        adj = J(0, 4, .)
        for (e = 1; e <= k; e++) {
            adj = adj \ (us[tpos[e]], ut[tpos[e]], e, 1) \ (ut[tpos[e]], us[tpos[e]], e, -1)
        }
        depth = J(nv, 1, .)
        pv = J(nv, 1, 0)
    ple = J(nv, 1, 0)
    psg = J(nv, 1, 0)
        for (u = 1; u <= nv; u++) {
            if (depth[u] < . | !any(adj[, 1] :== u)) continue
            depth[u] = 0
            queue = u
            x = 1
            while (x <= rows(queue)) {
                v = queue[x]
                x++
                for (y = 1; y <= rows(adj); y++) {
                    if (adj[y, 1] != v) continue
                    if (depth[adj[y, 2]] < .) continue
                    depth[adj[y, 2]] = depth[v] + 1
                    pv[adj[y, 2]] = v
                    ple[adj[y, 2]] = adj[y, 3]
                    psg[adj[y, 2]] = adj[y, 4]
                    queue = queue \ adj[y, 2]
                }
            }
        }
        d = rows(dpos)
        E = J(d, k, 0)
        for (dd = 1; dd <= d; dd++) {
            u = us[dpos[dd]]
            v = ut[dpos[dd]]
            while (depth[u] > depth[v]) {
                E[dd, ple[u]] = E[dd, ple[u]] - psg[u]
                u = pv[u]
            }
            while (depth[v] > depth[u]) {
                E[dd, ple[v]] = E[dd, ple[v]] + psg[v]
                v = pv[v]
            }
            while (u != v) {
                E[dd, ple[u]] = E[dd, ple[u]] - psg[u]
                u = pv[u]
                E[dd, ple[v]] = E[dd, ple[v]] + psg[v]
                v = pv[v]
            }
        }
        bpos = bpos \ &(idx[tpos])
        bmat = bmat \ &(cross(E, E))
    }
    C.kept = sort(keptall, 1)
    posmap = J(M, 1, 0)
    posmap[C.kept] = (1::rows(C.kept))
    C.dropped = selectindex(posmap :== 0)
    C.G = I(rows(C.kept))
    for (nbk = 1; nbk <= rows(bpos); nbk++) {
        pp = posmap[*bpos[nbk]]
        EE = *bmat[nbk]
        for (i1 = 1; i1 <= rows(pp); i1++) {
            for (i2 = 1; i2 <= rows(pp); i2++) C.G[pp[i1], pp[i2]] = C.G[pp[i1], pp[i2]] + EE[i1, i2]
        }
    }
    C.G = 0.5 :* (C.G + C.G')
}

/* Endpoint vertex of link l at calendar index u: (role, key) with roles 1 reference,
   2 pre_target, 3 pre_shared, 4 pre (key = control cohort index), 5 anticipation_target,
   6 trend_post, 7 post_target. */
real rowvector _cdid_xp_endpoint(struct cdid_xp_cat scalar C, real scalar l, real scalar u,
    real scalar delta, string scalar pt, real scalar pt_start)
{
    real scalar g, tv
    g = C.g[l]
    tv = C.cal[u]
    if (tv == C.ref[l]) return((1, 0))
    if (tv + delta < g) {
        if (pt == "all" & tv >= pt_start) return((3, 0))
        if (C.role[l] == 1) return((2, 0))
        if (pt == "ctrl" & tv >= pt_start) return((3, 0))
        return((4, C.ci[l]))
    }
    if (C.role[l] == 1) return(((tv < g ? 5 : 7), 0))
    return((6, 0))
}

/* Design on the used links: vertices (target index, date index, role, key), edges,
   anchored components and retained parameters (one gauge vertex dropped per
   unanchored component). */
struct cdid_xp_design {
    real matrix VX, A
    real colvector from, to, anch, vpar
    real scalar P
}

struct cdid_xp_design scalar _cdid_xp_design(struct cdid_xp_cat scalar C, real colvector used,
    real scalar delta, string scalar pt, real scalar pt_start)
{
    struct cdid_xp_design scalar D
    real scalar L, l, i, ru, rv, nv, u
    real colvector comp, keepv
    real rowvector ep
    transmorphic VA
    string scalar key
    L = rows(used)
    VA = asarray_create()
    D.VX = J(0, 4, .)
    D.from = J(L, 1, .)
    D.to = J(L, 1, .)
    for (l = 1; l <= L; l++) {
        for (i = 1; i <= 2; i++) {
            u = (i == 1 ? C.cs[used[l]] : C.ct[used[l]])
            ep = _cdid_xp_endpoint(C, used[l], u, delta, pt, pt_start)
            key = strofreal(C.gi[used[l]]) + "|" + strofreal(u) + "|" + strofreal(ep[1]) + "|" + strofreal(ep[2])
            if (!asarray_contains(VA, key)) {
                D.VX = D.VX \ (C.gi[used[l]], u, ep)
                asarray(VA, key, rows(D.VX))
            }
            if (i == 1) D.from[l] = asarray(VA, key)
            else D.to[l] = asarray(VA, key)
        }
    }
    nv = rows(D.VX)
    comp = (1::nv)
    for (l = 1; l <= L; l++) {
        ru = D.from[l]
        while (comp[ru] != ru) ru = comp[ru]
        rv = D.to[l]
        while (comp[rv] != rv) rv = comp[rv]
        if (ru != rv) comp[rv] = ru
    }
    for (i = 1; i <= nv; i++) {
        ru = i
        while (comp[ru] != ru) ru = comp[ru]
        comp[i] = ru
    }
    D.anch = J(nv, 1, 0)
    for (i = 1; i <= nv; i++) if (D.VX[i, 3] == 1) D.anch = D.anch :| (comp :== comp[i])
    keepv = (D.VX[, 3] :!= 1)
    for (i = 1; i <= nv; i++) if (!D.anch[i] & comp[i] == i) keepv[i] = 0
    D.P = sum(keepv)
    D.vpar = J(nv, 1, 0)
    if (D.P) D.vpar[selectindex(keepv)] = (1::D.P)
    D.A = J(L, D.P, 0)
    for (l = 1; l <= L; l++) {
        if (D.vpar[D.to[l]]) D.A[l, D.vpar[D.to[l]]] = D.A[l, D.vpar[D.to[l]]] + 1
        if (D.vpar[D.from[l]]) D.A[l, D.vpar[D.from[l]]] = D.A[l, D.vpar[D.from[l]]] - 1
    }
    return(D)
}

/* Retained parameter of the anchored vertex (target j, date u) with one of the roles. */
real scalar _cdid_xp_param(struct cdid_xp_design scalar D, real scalar j, real scalar u, real rowvector roles)
{
    real colvector h
    real scalar r
    for (r = 1; r <= cols(roles); r++) {
        h = selectindex((D.VX[, 1] :== j) :& (D.VX[, 2] :== u) :& (D.VX[, 3] :== roles[r]))
        if (rows(h) == 1) {
            if (D.anch[h[1]] & D.vpar[h[1]]) return(D.vpar[h[1]])
            return(0)
        }
    }
    return(0)
}


/* Covariate-adjusted moment of a control link (IPW, OR or DR): moment, unit IF added to
   Phi[, l], and the unit loadings (index, weight) for the within covariance. */
struct _cdid_xp_lg {
    real scalar ready, uncond, B_c, r_raw, M_k
    real colvector b_c, act_c, p_u, w_u, D_u, idx_u, idx_all, A_r, beta
    real matrix Z_u, Z_c, H_chol
}

void _cdid_xp_logit(struct _cdid_xp_lg scalar R, real matrix Ztr, real colvector w, real colvector idx_t,
    real colvector idx_c, real colvector idx_all, real scalar nA, real scalar eps, real scalar maxit)
{
    struct cdid_logit_result scalar lf
    real scalar M_k, r_c, sub
    real colvector w_c, v_c, cpos
    real matrix lv
    M_k = sum(w[idx_all])
    w_c = w[idx_c]
    sub = rows(idx_c) < rows(idx_all)
    R.r_raw = (sub ? sum(w_c) / M_k : 1)
    r_c = max((R.r_raw, 1e-6))
    v_c = w_c :/ r_c
    R.M_k = M_k
    R.idx_all = idx_all
    R.uncond = (cols(Ztr) == 1)
    if (R.uncond) {
        R.b_c = v_c :* (sum(w[idx_t]) / sum(v_c))
        R.B_c = sum(R.b_c)
        R.ready = 1
        return
    }
    R.idx_u = idx_t \ idx_c
    R.D_u = J(rows(idx_t), 1, 1) \ J(rows(idx_c), 1, 0)
    R.w_u = w[idx_t] \ v_c
    R.Z_u = Ztr[R.idx_u, .]
    R.Z_c = Ztr[idx_c, .]
    lf = cdid_logit_fit(R.Z_u, R.D_u, R.w_u, eps, maxit, 1)
    R.beta = lf.beta
    lv = _cdid_ipw_logit_values(R.Z_u * lf.beta)
    R.p_u = lv[, 1]
    cpos = (rows(idx_t) + 1)::rows(R.idx_u)
    R.act_c = lv[cpos, 3]
    R.b_c = v_c :* lv[cpos, 2]
    R.B_c = sum(R.b_c)
    R.H_chol = cholesky(cross(R.Z_u, R.Z_u :* (R.w_u :* R.p_u :* (1 :- R.p_u))) :/ nA)
    if (missing(R.H_chol)) _cdid_ipw_fail("Weighted propensity Hessian is not positive definite.", 459)
    R.A_r = J(0, 1, .)
    if (sub & R.r_raw > 1e-6) R.A_r = (colsum(R.Z_c :* (w_c :* R.p_u[cpos])) :/ (nA * r_c * r_c))'
    R.ready = 1
}

void _cdid_xp_logit_if(struct _cdid_xp_lg scalar R, real colvector res, real colvector w, real scalar nA,
    real colvector idx_c, real matrix Phi, real scalar l)
{
    real colvector G, h
    real scalar lam
    if (R.uncond) return
    G = (colsum(R.Z_c :* (R.b_c :* res :* R.act_c)) :/ R.B_c)'
    h = solveupper(R.H_chol', solvelower(R.H_chol, G))
    Phi[R.idx_u, l] = Phi[R.idx_u, l] + R.w_u :* (R.D_u - R.p_u) :* (R.Z_u * h)
    if (rows(R.A_r)) {
        lam = R.A_r' * h
        Phi[R.idx_all, l] = Phi[R.idx_all, l] - (w[R.idx_all] :* (nA / R.M_k)) :* (R.r_raw * lam)
        Phi[idx_c, l] = Phi[idx_c, l] + (w[idx_c] :* (nA / R.M_k)) :* lam
    }
}

/* Weighted least squares with column pivoting (rank tolerance 1e-7): retained columns and
   their triangular factor. A rank-deficient design is accepted when its dropped columns are
   exact linear aliases of the retained ones on the fitting rows and on the target rows Xt. */
void _cdid_xp_wls(real matrix X, real colvector w, real colvector y, real colvector piv,
    real matrix Rk, real colvector gam, | real matrix Xt)
{
    real matrix Xw, H, Rf, A
    real rowvector tau, pv
    real colvector d, qy, drop
    real scalar rk, mx, err, tol
    Xw = X :* sqrt(w)
    hqrdp(Xw, H = ., tau = ., Rf = ., pv = .)
    d = abs(diagonal(Rf))
    mx = max(d)
    rk = sum(d :> 1e-7 * mx)
    if (rk < 1 | mx <= 0) _cdid_ipw_fail("The outcome-regression design has no numerical rank on a control cohort's observed pairs.", 459)
    piv = pv[1..rk]'
    Rk = Rf[1..rk, 1..rk]
    if (rk < cols(X)) {
        drop = pv[(rk + 1)..cols(X)]'
        A = solveupper(Rk, hqrdmultq1t(H, tau, Xw[, drop'])[1..rk, ])
        err = max(abs(X[, drop'] - X[, piv'] * A))
        tol = 64 * epsilon(1) * cols(X) * max((1, max(abs(X))))
        if (args() >= 7) {
            err = max((err, max(abs(Xt[, drop'] - Xt[, piv'] * A))))
            tol = 64 * epsilon(1) * cols(X) * max((1, max(abs(X)), max(abs(Xt))))
        }
        if (err >= . | err > tol)
            _cdid_ipw_fail("The outcome-regression design is rank deficient on a control cohort's observed pairs without an exact alias that also holds on the target cohort; the counterfactual mean is not identified without extrapolation. Remove the redundant covariate.", 459)
    }
    qy = hqrdmultq1t(H, tau, y :* sqrt(w))
    gam = solveupper(Rk, qy[1..rk])
}

/* Cross-pair fit: every cell (g, t) of the targets, post cells from the design under the
   requested trends, pre cells from it when their control level is shared there and
   otherwise from the same moments under shared control pre-trends. att and the unit IF
   (analysis units) are missing for unavailable cells; the reference cell is zero. */
struct cdid_xp_fit {
    real matrix cells, IF
    real colvector att, avail
    real colvector mhat, cm, cp
    real matrix Phi, V, Lam
    real rowvector diag
    string scalar trigger_reason
}

/* Moments with covariates (R evaluate_cross_pair_moments, adjusted path): target links are
   weighted means; control links are adjusted by IPW, OR or DR toward the target cohort. */
void _cdid_xp_cov_moments(struct cdid_xp_fit scalar F, struct cdid_xp_cat scalar C, real matrix Y,
    real colvector w, real colvector used, real matrix Z, string scalar method, real scalar eps,
    real scalar maxit, string scalar weighting, struct cdid_gmm_unit_gram scalar src)
{
    struct _cdid_xp_lg scalar LG
    real scalar L, nA, l, gk, W_t, B_c, mu_bar, e_bar, a, b, sg, sumw, isor, isdr
    real colvector idx, idx_t, idx_c, idx_all, dy, om, w_t, w_c, piv, gam, mu_t, mu_c, e_c, xbar_t, cx, wadj, xatt, mrk
    real matrix Rk, Ztr
    struct _cdid_xp_vec colvector lidx, lom
    transmorphic LC, OK, TKR
    string scalar key
    real scalar nreq, rix
    pointer() scalar pnr
    L = rows(used)
    nreq = cdid_nuis_requested()
    pnr = findexternal("__cdid_xpn_base")
    if (pnr == NULL) nreq = 0
    OK = asarray_create("real", 3)
    TKR = asarray_create()
    nA = rows(Y)
    isor = (method == "reg")
    isdr = (method == "dr")
    Ztr = Z
    if (!isor & !isdr & cols(Z) > 1) {
        if (all(Z[, 1] :== Z[1, 1]) & all(Z[, 2..cols(Z)] :== Z[1, 2..cols(Z)])) Ztr = Z[, 1]
    }
    F.mhat = J(L, 1, .)
    F.Phi = J(nA, L, 0)
    lidx = _cdid_xp_vec(L, 1)
    lom = _cdid_xp_vec(L, 1)
    LC = asarray_create()
    for (l = 1; l <= L; l++) {
        idx = *C.obs[C.oi[used[l]], C.mid[used[l]]]
        dy = Y[idx, C.ct[used[l]]] - Y[idx, C.cs[used[l]]]
        if (C.role[used[l]] == 1) {
            om = w[idx] :/ sum(w[idx])
            F.mhat[l] = sum(om :* dy)
            F.Phi[idx, l] = nA :* om :* (dy :- F.mhat[l])
            lidx[l].v = idx
            lom[l].v = om
            continue
        }
        gk = selectindex(C.coh :== C.g[used[l]])
        idx_t = selectindex((C.ucoh :== C.coh[gk[1]]) :& (w :> 0))
        idx_all = selectindex((C.ucoh :== C.coh[C.ci[used[l]]]) :& (w :> 0))
        w_t = w[idx_t]
        W_t = sum(w_t)
        idx_c = idx
        w_c = w[idx_c]
        if (!isor) {
            key = strofreal(C.gi[used[l]]) + ":" + strofreal(C.ci[used[l]]) + ":" + strofreal(C.mid[used[l]])
            if (asarray_contains(LC, key)) {
                LG = asarray(LC, key)
                if (nreq & asarray_contains(TKR, key)) cdid_xpn_base_link(asarray(TKR, key), used[l], C.cal[C.cs[used[l]]], C.cal[C.ct[used[l]]])
            }
            else {
                LG = _cdid_xp_lg()
                _cdid_xp_logit(LG, (isdr ? Z : Ztr), w, idx_t, idx_c, idx_all, nA, eps, maxit)
                asarray(LC, key, LG)
                if (nreq & !LG.uncond) asarray(TKR, key, cdid_xpn_base_trt(LG, w, nA, C, used[l], Y))
            }
        }
        if (!isor & !isdr) {
            om = LG.b_c :/ LG.B_c
            F.mhat[l] = sum(om :* dy)
            F.Phi[idx_c, l] = F.Phi[idx_c, l] + nA :* om :* (dy :- F.mhat[l])
            _cdid_xp_logit_if(LG, dy :- F.mhat[l], w, nA, idx_c, F.Phi, l)
        }
        else {
            _cdid_xp_wls(Z[idx_c, .], w_c, dy, piv, Rk, gam, Z[idx_t, .])
            mu_t = Z[idx_t, piv'] * gam
            mu_c = Z[idx_c, piv'] * gam
            e_c = dy - mu_c
            if (nreq) {
                if (asarray_contains(OK, (C.ci[used[l]], C.cs[used[l]], C.ct[used[l]])))
                    cdid_xpn_base_link(asarray(OK, (C.ci[used[l]], C.cs[used[l]], C.ct[used[l]])), used[l], C.cal[C.cs[used[l]]], C.cal[C.ct[used[l]]])
                else asarray(OK, (C.ci[used[l]], C.cs[used[l]], C.ct[used[l]]),
                    cdid_xpn_base_or(Z[idx_c, piv'], w_c, e_c, piv, Rk, gam, idx_c, nA, C, used[l], cols(Z)))
            }
            xbar_t = (colsum(Z[idx_t, piv'] :* w_t) :/ W_t)'
            mu_bar = sum(w_t :* mu_t) / W_t
            if (isor) {
                F.mhat[l] = mu_bar
                cx = solveupper(Rk, solvelower(Rk', xbar_t))
                om = w_c :* (Z[idx_c, piv'] * cx)
                F.Phi[idx_t, l] = F.Phi[idx_t, l] + (nA :* w_t :/ W_t) :* (mu_t :- mu_bar)
                F.Phi[idx_c, l] = F.Phi[idx_c, l] + nA :* om :* e_c
            }
            else {
                B_c = LG.B_c
                e_bar = sum(LG.b_c :* e_c) / B_c
                F.mhat[l] = mu_bar + e_bar
                xatt = (colsum(Z[idx_c, piv'] :* LG.b_c) :/ B_c)'
                cx = solveupper(Rk, solvelower(Rk', xbar_t - xatt))
                wadj = w_c :* (Z[idx_c, piv'] * cx)
                om = (LG.b_c :/ B_c) + wadj
                F.Phi[idx_t, l] = F.Phi[idx_t, l] + (nA :* w_t :/ W_t) :* (mu_t :- mu_bar)
                F.Phi[idx_c, l] = F.Phi[idx_c, l] + (nA :* LG.b_c :/ B_c) :* (e_c :- e_bar) + nA :* wadj :* e_c
                _cdid_xp_logit_if(LG, e_c :- e_bar, w, nA, idx_c, F.Phi, l)
            }
        }
        lidx[l].v = idx_c
        lom[l].v = om
    }
    src.mean = colsum(F.Phi) :/ nA
    F.Phi = F.Phi :- src.mean
    if (weighting == "within") {
        sumw = sum(w)
        F.V = J(L, L, 0)
        mrk = J(nA, 1, 0)
        for (a = 1; a <= L; a++) {
            mrk[lidx[a].v] = lom[a].v :/ w[lidx[a].v]
            for (b = a; b <= L; b++) {
                if (C.gi[used[a]] != C.gi[used[b]]) continue
                sg = (C.ct[used[a]] == C.ct[used[b]]) + (C.cs[used[a]] == C.cs[used[b]]) -
                    (C.ct[used[a]] == C.cs[used[b]]) - (C.cs[used[a]] == C.ct[used[b]])
                if (sg == 0) continue
                F.V[a, b] = sg * sumw * sum(mrk[lidx[b].v] :* lom[b].v)
                F.V[b, a] = F.V[a, b]
            }
            mrk[lidx[a].v] = J(rows(lidx[a].v), 1, 0)
        }
        F.V = 0.5 :* (F.V + F.V')
    }
    else {
        F.V = quadcross(F.Phi, F.Phi) :/ nA
        F.V = 0.5 :* (F.V + F.V')
    }
}


real rowvector _cdid_xp_spec(struct cdid_xp_cat scalar C, struct cdid_xp_design scalar D, real scalar j,
    real scalar u, real scalar from, real scalar refi, real scalar post)
{
    real rowvector h
    real scalar pg, pl
    h = J(1, D.P, 0)
    if (post) {
        pg = _cdid_xp_param(D, j, u, 7)
        pl = _cdid_xp_param(D, j, u, 6)
        if (!pg | !pl) return(J(1, 0, .))
        h[pg] = 1
        h[pl] = -1
        return(h)
    }
    if (u != refi) {
        pg = _cdid_xp_param(D, j, u, (2, 5))
        pl = _cdid_xp_param(D, j, u, (3, 6))
        if (!pg | !pl) return(J(1, 0, .))
        h[pg] = h[pg] + 1
        h[pl] = h[pl] - 1
    }
    if (from != refi) {
        pg = _cdid_xp_param(D, j, from, (2, 5))
        pl = _cdid_xp_param(D, j, from, (3, 6))
        if (!pg | !pl) return(J(1, 0, .))
        h[pg] = h[pg] - 1
        h[pl] = h[pl] + 1
    }
    return(h)
}

void _cdid_xp_gmm(struct cdid_xp_fit scalar F, struct cdid_xp_design scalar D, real matrix H,
    real colvector cix, real matrix G, real scalar M_all, string scalar weighting, real scalar trigger,
    real scalar ridge, struct cdid_gmm_unit_gram scalar src, real scalar main)
{
    struct cdid_gmm_result scalar sol
    if (rows(H) == 0) return
    if (weighting == "identity") sol = cdid_gmm_identity_map(F.mhat, D.A, H, 2^53 - 1)
    else sol = _cdid_gmm_fit_core(F.mhat, D.A, F.V, G, M_all, H, "optimal", trigger, ridge,
        2^53 - 1, src, weighting != "within")
    F.att[cix] = sol.targets
    F.avail[cix] = J(rows(cix), 1, 1)
    F.Lam[., cix] = sol.target_loadings'
    F.IF[, cix] = F.Phi * sol.target_loadings'
    if (main) {
        F.diag = (sol.lambda, sol.rho, sol.covariance_rank, sol.d_min, sol.d_max, sol.criterion)
        F.trigger_reason = sol.trigger_reason
    }
}

struct cdid_xp_fit scalar cdid_xp_solve(struct cdid_xp_cat scalar C, real matrix Y, real colvector w,
    real colvector used, real matrix G, real scalar M_all, real scalar delta, string scalar pt,
    real scalar pt_start, string scalar base, string scalar weighting, real scalar trigger, real scalar ridge,
    | real matrix Z, string scalar method, real scalar eps, real scalar maxit)
{
    struct cdid_xp_fit scalar F
    struct cdid_xp_design scalar D, Dp
    struct cdid_gmm_unit_gram scalar src
    real scalar L, nA, l, j, a, b, sg, ov, sumw, u, from, refi, n, post, k
    real colvector idx, dy, om, Wl, mrk, cm, cp
    real matrix H, Hp, spec
    real rowvector h
    L = rows(used)
    nA = rows(Y)
    nA = rows(Y)
    src.n_units = nA
    src.rows_consumed = nA
    src.normalization_count = nA
    src.normalization_factor = 1 / nA
    src.sample_rank_upper = nA - 1
    src.retained_if = J(0, 0, .)
    if (args() < 14) Z = J(rows(Y), 1, 1)
    if (args() < 15) method = "dr"
    if (args() < 16) eps = 1e-8
    if (args() < 17) maxit = 25
    if (cols(Z) > 1 | !all(Z :== Z[1, 1])) {
        _cdid_xp_cov_moments(F, C, Y, w, used, Z, method, eps, maxit, weighting, src)
    }
    else {
        /* Moments and unit influence functions. */
        F.mhat = J(L, 1, .)
        F.Phi = J(nA, L, 0)
        Wl = J(L, 1, .)
        for (l = 1; l <= L; l++) {
            idx = *C.obs[C.oi[used[l]], C.mid[used[l]]]
            dy = Y[idx, C.ct[used[l]]] - Y[idx, C.cs[used[l]]]
            Wl[l] = sum(w[idx])
            om = w[idx] :/ Wl[l]
            F.mhat[l] = sum(om :* dy)
            F.Phi[idx, l] = nA :* om :* (dy :- F.mhat[l])
        }
        src.n_units = nA
        src.rows_consumed = nA
        src.normalization_count = nA
        src.normalization_factor = 1 / nA
        src.sample_rank_upper = nA - 1
        src.mean = colsum(F.Phi) :/ nA
        src.retained_if = J(0, 0, .)
        F.Phi = F.Phi :- src.mean
        if (weighting == "within") {
            sumw = sum(w)
            F.V = J(L, L, 0)
            mrk = J(nA, 1, 0)
            for (a = 1; a <= L; a++) {
                for (b = a; b <= L; b++) {
                    if (C.gi[used[a]] != C.gi[used[b]] | C.ci[used[a]] != C.ci[used[b]]) continue
                    sg = (C.ct[used[a]] == C.ct[used[b]]) + (C.cs[used[a]] == C.cs[used[b]]) -
                        (C.ct[used[a]] == C.cs[used[b]]) - (C.cs[used[a]] == C.ct[used[b]])
                    if (sg == 0) continue
                    if (C.mid[used[a]] == C.mid[used[b]]) ov = Wl[a]
                    else {
                        idx = *C.obs[C.oi[used[a]], C.mid[used[a]]]
                        mrk[idx] = J(rows(idx), 1, 1)
                        idx = *C.obs[C.oi[used[b]], C.mid[used[b]]]
                        ov = sum(w[idx] :* mrk[idx])
                        idx = *C.obs[C.oi[used[a]], C.mid[used[a]]]
                        mrk[idx] = J(rows(idx), 1, 0)
                    }
                    F.V[a, b] = sg * sumw * ov / (Wl[a] * Wl[b])
                    F.V[b, a] = F.V[a, b]
                }
            }
            F.V = 0.5 :* (F.V + F.V')
        }
        else {
            F.V = quadcross(F.Phi, F.Phi) :/ nA
            F.V = 0.5 :* (F.V + F.V')
        }
    }
    src.covariance = F.V
    cdid_xp_targets(F, C, used, G, M_all, delta, pt, pt_start, base, weighting, trigger, ridge, src)
    return(F)
}

/* Reported cells from the fitted moments F.mhat, their unit influence functions F.Phi and
   weight covariance F.V on the used links. */
void cdid_xp_targets(struct cdid_xp_fit scalar F, struct cdid_xp_cat scalar C, real colvector used,
    real matrix G, real scalar M_all, real scalar delta, string scalar pt, real scalar pt_start,
    string scalar base, string scalar weighting, real scalar trigger, real scalar ridge,
    struct cdid_gmm_unit_gram scalar src)
{
    struct cdid_xp_design scalar D, Dp
    real scalar nA, j, u, from, refi, n, post, L
    real colvector cm, cp
    real matrix H, Hp
    real rowvector h
    nA = rows(F.Phi)
    L = rows(used)
    /* Cells: every target and calendar date. */
    n = 0
    for (j = 1; j <= rows(C.tg); j++) n = n + rows(C.cal)
    F.cells = J(n, 3, .)
    n = 0
    for (j = 1; j <= rows(C.tg); j++) {
        for (u = 1; u <= rows(C.cal); u++) {
            n++
            F.cells[n, ] = (C.tg[j], C.cal[u], (C.cal[u] + delta >= C.tg[j]))
        }
    }
    F.att = J(n, 1, .)
    F.avail = J(n, 1, 0)
    F.IF = J(nA, n, .)
    F.diag = J(1, 6, .)
    F.trigger_reason = "none"
    D = _cdid_xp_design(C, used, delta, pt, pt_start)
    Dp = _cdid_xp_design(C, used, delta, "ctrl", min(C.cal))
    H = J(0, D.P, 0)
    Hp = J(0, Dp.P, 0)
    cm = J(0, 1, .)
    cp = J(0, 1, .)
    n = 0
    for (j = 1; j <= rows(C.tg); j++) {
        refi = (C.tref[j] < . ? selectindex(C.cal :== C.tref[j]) : .)
        for (u = 1; u <= rows(C.cal); u++) {
            n++
            if (refi >= .) continue
            post = F.cells[n, 3] & (C.cal[u] >= C.tg[j])
            if (!post & u == refi & base == "universal") {
                F.att[n] = 0
                F.avail[n] = 1
                F.IF[, n] = J(nA, 1, 0)
                continue
            }
            from = (post ? . : (base == "universal" ? refi : u - 1))
            if (!post & from < 1) continue
            h = _cdid_xp_spec(C, D, j, u, from, refi, post)
            if (cols(h)) {
                H = H \ h
                cm = cm \ n
                continue
            }
            if (post) continue
            h = _cdid_xp_spec(C, Dp, j, u, from, refi, 0)
            if (cols(h)) {
                Hp = Hp \ h
                cp = cp \ n
            }
        }
    }
    F.cm = cm
    F.cp = cp
    F.Lam = J(L, rows(F.cells), 0)
    _cdid_xp_gmm(F, D, H, cm, G, M_all, weighting, trigger, ridge, src, 1)
    _cdid_xp_gmm(F, Dp, Hp, cp, G, M_all, weighting, trigger, ridge, src, 0)
}

/* Exact-mask spanning forest of the catalogue (R within_mask_spanning_forest): comparisons
   grouped by target cohort, control pool and the exact set of units observed at both
   dates; kept rows ascending and their full-catalogue metric G = I + E'E. */
void cdid_wsp_forest(struct _cdid_wsp_data scalar d, struct cdid_cat scalar C,
    real colvector kept, real matrix G)
{
    struct cdid_xp_cat scalar X
    struct _cdid_xp_vec colvector mk
    real colvector idx_t, idx_c, u, key, ord, grp
    real matrix fp
    real scalar M, m, e, ng, j, hit
    M = rows(C.r_g)
    u = sin((1::d.nA) :* 1.618033988749895) :+ 2
    mk = _cdid_xp_vec(M, 1)
    fp = J(M, 5, .)
    for (m = 1; m <= M; m++) {
        _cdid_wsp_units(d, C, m, idx_t, idx_c)
        mk[m].v = sort(idx_t \ idx_c, 1)
        fp[m, .] = (C.r_g[m], C.r_pool[m], rows(mk[m].v), sum(u[mk[m].v]), m)
    }
    ord = order(fp, (1, 2, 3, 4, 5))
    key = J(M, 1, 0)
    grp = J(0, 1, .)
    ng = 0
    for (e = 1; e <= M; e++) {
        m = ord[e]
        hit = 0
        if (e > 1) {
            if (all(fp[m, 1..4] :== fp[ord[e - 1], 1..4])) {
                for (j = 1; j <= rows(grp); j++) {
                    if (mk[grp[j]].v == mk[m].v) {
                        hit = key[grp[j]]
                        break
                    }
                }
            }
            else grp = J(0, 1, .)
        }
        if (!hit) {
            ng++
            hit = ng
            grp = grp \ m
        }
        key[m] = hit
    }
    X.g = C.r_g
    X.gi = key
    X.ci = J(M, 1, 0)
    X.mid = J(M, 1, 0)
    X.cs = C.r_s
    X.ct = C.r_t
    X.cal = C.calendar
    cdid_xp_forest(X)
    kept = X.kept
    G = X.G
}

end
