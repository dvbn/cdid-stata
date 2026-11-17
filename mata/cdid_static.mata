*! Static binary treatment with entry, exit and re-entry (dvar()): complete-panel estimator.
version 19.5
mata:

/* Transition cells c = (t, d): switchers G (D_{t-1} = d, D_t = 1 - d) and stayers C
   (D_{t-1} = D_t = d) among units observed at t - 1 and t; sign sigma = 1 - 2d. The cell
   contrast theta_c compares the mean outcome change of switchers with that of stayers,
   adjusted by IPW (odds of a switcher logit on the covariates at t - 1, over switchers and
   stayers), by a stayer outcome regression with a cell intercept and slopes shared within
   a pooling block (reg), or both (dr). tau_c = sigma theta_c. The aggregate weights the
   signed cell contrasts of the targets by their switcher shares, with the influence of
   the random shares. */
struct _cdid_st_cell {
    real matrix Z, ell, Xt, chiCov, S, chib
    real colvector act, odds, chimC, zbar, cov, b, chia, gam
    real scalar mC, piC, alpha, orblock
}

struct _cdid_st_nf {
    real scalar ncell, nblock, ds, P, N
    real colvector cgrp, cs, cu, cblock, calpha, theta, U, iter, nfit, nobs
    real matrix bbeta, IF
}

struct cdid_st_fit {
    real matrix cells
    string colvector cid
    real colvector theta, tau
    real matrix IF
    real scalar agg, agg_avail, total_sw
    real colvector agg_if, weights
    string colvector status
    real scalar has_rep, rep_q, rep_j
    struct _cdid_st_cell colvector repW
    struct _cdid_st_nf scalar rsel, rimp
    real matrix rgrp
}

struct _cdid_st_or {
    real colvector alpha, cells
    real matrix b, S, mz, ell_score
    real scalar p
}

/* Unweighted logit of switchers against stayers on the risk set; returns odds, active
   flags, probabilities (all units, zero off the risk set) and the Cholesky factor of
   X'diag(p(1-p))X / n. */
void _cdid_st_logit(real matrix Xt, real colvector G, real colvector R, real scalar n, real scalar eps,
    real scalar maxit, real colvector odds, real colvector act, real colvector prob, real matrix Hc,
    | real colvector beta)
{
    struct cdid_logit_result scalar lf
    real colvector ri
    real matrix lv
    ri = selectindex(R)
    lf = cdid_logit_fit(Xt[ri, .], G[ri], J(rows(ri), 1, 1), eps, maxit, 1)
    lv = _cdid_ipw_logit_values(Xt[ri, .] * lf.beta)
    odds = act = prob = J(n, 1, 0)
    prob[ri] = lv[., 1]
    odds[ri] = lv[., 2]
    act[ri] = lv[., 3]
    Hc = cholesky(cross(Xt[ri, .], lv[., 1] :* (1 :- lv[., 1]), Xt[ri, .]) :/ n)
    if (missing(Hc)) _cdid_st_fail("The switcher logit Hessian is not positive definite.")
    if (args() >= 11) beta = lf.beta
}

void _cdid_st_fail(string scalar msg)
{
    errprintf("%s\n", msg)
    exit(459)
}

void _cdid_st_targets(real matrix D, real rowvector dirs, real rowvector tsel, real colvector tval,
    struct cdid_st_fit scalar F, real matrix Gm)
{
    real scalar n, T, nc, k, t, d, j, J
    real colvector G, C, tgt
    n = rows(D)
    T = cols(D)
    /* Candidate cells (t = 2..T, d = 0 then 1). */
    nc = 2 * (T - 1)
    F.cells = J(nc, 9, .)
    F.cid = J(nc, 1, "")
    Gm = J(n, nc, 0)
    k = 0
    for (t = 2; t <= T; t++) for (d = 0; d <= 1; d++) {
        k++
        G = (D[., t - 1] :== d) :& (D[., t] :== 1 - d)
        C = (D[., t - 1] :== d) :& (D[., t] :== d)
        Gm[., k] = G
        F.cells[k, .] = (tval[t], tval[t - 1], d, 1 - 2 * d, sum(G :| C), sum(G), sum(C), 0, k)
        F.cid[k] = "t" + strofreal(tval[t], "%21.0g") + "_" + (d == 0 ? "entry" : "exit")
    }
    /* Declared targets: cells with switchers in the requested directions, or the listed dates. */
    tgt = J(nc, 1, 0)
    if (cols(tsel) == 0) {
        for (k = 1; k <= nc; k++) tgt[k] = any(dirs :== F.cells[k, 3]) & F.cells[k, 6] >= 1
        if (!any(tgt)) for (k = 1; k <= nc; k++) tgt[k] = any(dirs :== F.cells[k, 3])
    }
    else for (k = 1; k <= nc; k++) tgt[k] = any(dirs :== F.cells[k, 3]) & any(tsel :== F.cells[k, 1])
    tgt = selectindex(tgt)
    F.cells = F.cells[tgt, .]
    F.cid = F.cid[tgt]
    Gm = Gm[., tgt']
    J = rows(tgt)
    F.status = J(J, 1, "ok")
}

/* Switcher-share aggregate over the declared targets, only when all are available. */
void _cdid_st_aggregate(struct cdid_st_fit scalar F, real matrix Gm, real scalar n, real scalar compute_if)
{
    real scalar Pi
    real colvector sth
    real matrix GA
    F.agg = .
    F.agg_avail = (rows(F.cells) > 0 & all(F.cells[., 8]))
    F.total_sw = sum(F.cells[., 6])
    F.weights = F.cells[., 6] :/ F.total_sw
    F.agg_if = J(0, 1, .)
    if (F.agg_avail) {
        sth = F.tau
        F.agg = sum(F.weights :* sth)
        if (compute_if) {
            Pi = F.total_sw / n
            GA = Gm :- (F.cells[., 6] :/ n)'
            F.agg_if = F.IF * F.weights + GA * ((sth :- F.agg) :/ Pi)
        }
    }
}

/* Outcome-regression blocks of the declared targets: one per target (separate), one for all
   (pooled), or a user partition "c c ... | c c ..." of cells written entry:<date> or exit:<date>.
   blk gives the block of each target; extra lists the training cells that are not targets
   (t index, d, block). */
void _cdid_st_or_blocks(string scalar pooling, struct cdid_st_fit scalar F, real colvector tval,
    real matrix D, real colvector blk, real matrix extra)
{
    string rowvector parts, toks
    string scalar tk, dir, dt, miss
    real scalar J, b, i, d, tv, ti, j, nb
    real matrix seen
    J = rows(F.cells)
    extra = J(0, 3, .)
    if (pooling == "separate" | pooling == "") {
        blk = (1::J)
        return
    }
    if (pooling == "pooled") {
        blk = J(J, 1, 1)
        return
    }
    blk = J(J, 1, 0)
    seen = J(0, 2, .)
    parts = tokens(pooling, "|")
    nb = 0
    for (i = 1; i <= cols(parts); i++) {
        if (parts[i] == "|") continue
        toks = tokens(parts[i])
        if (cols(toks) == 0) continue
        nb++
        for (b = 1; b <= cols(toks); b++) {
            tk = strlower(toks[b])
            if (strpos(tk, ":") == 0) _cdid_st_fail("outcomepooling(): cells are written entry:date or exit:date; " + tk + " is not.")
            dir = substr(tk, 1, strpos(tk, ":") - 1)
            dt = substr(tk, strpos(tk, ":") + 1, .)
            if (substr(dt, 1, 1) == "t") dt = substr(dt, 2, .)
            tv = strtoreal(dt)
            if ((dir != "entry" & dir != "exit") | tv >= .) _cdid_st_fail("outcomepooling(): cells are written entry:date or exit:date; " + toks[b] + " is not.")
            d = (dir == "exit")
            ti = selectindex(tval :== tv)
            if (rows(ti) != 1) _cdid_st_fail("outcomepooling(): " + toks[b] + " is not a transition of the panel (its date must be observed and follow the first date).")
            if (ti[1] < 2) _cdid_st_fail("outcomepooling(): " + toks[b] + " is not a transition of the panel (its date must be observed and follow the first date).")
            if (rows(seen)) if (any((seen[, 1] :== ti[1]) :& (seen[, 2] :== d))) _cdid_st_fail("outcomepooling(): a cell cannot belong to two blocks (" + toks[b] + ").")
            seen = seen \ (ti[1], d)
            j = selectindex((F.cells[, 1] :== tv) :& (F.cells[, 3] :== d))
            if (rows(j)) blk[j[1]] = nb
            else {
                if (!any((D[., ti[1] - 1] :== d) :& (D[., ti[1]] :== d))) _cdid_st_fail("outcomepooling(): training cell " + toks[b] + " has no stayers.")
                extra = extra \ (ti[1], d, nb)
            }
        }
    }
    if (any(blk :== 0)) {
        miss = ""
        for (j = 1; j <= J; j++) if (!blk[j]) miss = miss + (miss == "" ? "" : ", ") + (F.cells[j, 3] ? "exit:" : "entry:") + strofreal(F.cells[j, 1], "%21.0g")
        _cdid_st_fail("outcomepooling() omits the requested cell(s) " + miss + ".")
    }
}

struct cdid_st_fit scalar cdid_static_complete(real matrix Y, real matrix D, real matrix Xall, real scalar p,
    string scalar method, real rowvector dirs, real rowvector tsel, string scalar pooling, real scalar eps,
    real scalar maxit, real scalar compute_if, real colvector tval, | real scalar rep)
{
    struct cdid_st_fit scalar F
    struct _cdid_st_or colvector OR
    real colvector gam
    real scalar n, T, nc, k, t, d, j, J, c, use_logit, use_or, nb, b, mu1, mu0, Q, pi, q0, Pi
    real colvector G, C, R, dy, odds, act, prob, e, ell, ga, lam, kz, zb1, zb0, tgt, av, blocks, cl, s, om, sth
    real colvector blk, ex, ea
    real matrix Xt, Z, Hc, Gm, Zc, Sb, phi, GA, extra, em
    string scalar m
    n = rows(Y)
    T = cols(Y)
    m = (method == "or" ? "reg" : method)
    use_logit = (m == "ipw" | m == "dr")
    use_or = (m == "reg" | m == "dr")
    if (args() < 13) rep = 0
    _cdid_st_targets(D, dirs, tsel, tval, F, Gm)
    J = rows(F.cells)
    F.has_rep = rep
    F.rep_q = F.rep_j = 0
    if (rep) F.repW = _cdid_st_cell(J, 1)
    for (j = 1; j <= J; j++) {
        if (F.cells[j, 6] < 1) F.status[j] = "empty_switchers"
        else if (F.cells[j, 7] < 1) F.status[j] = "empty_stayers"
        F.cells[j, 8] = (F.status[j] == "ok")
    }
    av = selectindex(F.cells[., 8])
    F.theta = F.tau = J(J, 1, .)
    F.IF = J(n, J, .)
    if (compute_if & rows(av)) F.IF[., av'] = J(n, rows(av), 0)
    /* Outcome regressions by pooling block over the available targets. */
    if (use_or & rows(av)) {
        _cdid_st_or_blocks(pooling, F, tval, D, blk, extra)
        blocks = blk[av]
        nb = max(blk)
        OR = _cdid_st_or(nb, 1)
        for (b = 1; b <= nb; b++) {
            cl = select(av, blocks :== b)
            if (rows(cl) == 0) continue
            OR[b].cells = cl
            OR[b].p = p
            OR[b].alpha = J(rows(cl), 1, .)
            OR[b].mz = J(rows(cl), p, .)
            Sb = J(p, p, 0)
            kz = J(p, 1, 0)
            for (j = 1; j <= rows(cl); j++) {
                c = cl[j]
                k = F.cells[c, 9]
                t = 1 + ceil(k / 2)
                C = (D[., t - 1] :== F.cells[c, 3]) :& (D[., t] :== F.cells[c, 3])
                s = selectindex(C)
                dy = Y[s, t] - Y[s, t - 1]
                if (p) {
                    Z = Xall[s, ((t - 2) * p + 1)..((t - 1) * p)]
                    if (missing(Z)) _cdid_st_fail("Covariates must be nonmissing at date t - 1 for every unit at risk of switching at t.")
                    OR[b].mz[j, .] = mean(Z)
                    Zc = Z :- OR[b].mz[j, .]
                    Sb = Sb + cross(Zc, Zc)
                    kz = kz + cross(Zc, dy :- mean(dy))
                }
                OR[b].alpha[j] = mean(dy)
            }
            /* Training cells that are not targets. */
            ex = (rows(extra) ? selectindex(extra[, 3] :== b) : J(0, 1, .))
            ea = J(rows(ex), 1, .)
            em = J(rows(ex), p, .)
            for (j = 1; j <= rows(ex); j++) {
                t = extra[ex[j], 1]
                C = (D[., t - 1] :== extra[ex[j], 2]) :& (D[., t] :== extra[ex[j], 2])
                s = selectindex(C)
                dy = Y[s, t] - Y[s, t - 1]
                if (p) {
                    Z = Xall[s, ((t - 2) * p + 1)..((t - 1) * p)]
                    if (missing(Z)) _cdid_st_fail("Covariates must be nonmissing at date t - 1 for every unit at risk of switching at t.")
                    em[j, .] = mean(Z)
                    Zc = Z :- em[j, .]
                    Sb = Sb + cross(Zc, Zc)
                    kz = kz + cross(Zc, dy :- mean(dy))
                }
                ea[j] = mean(dy)
            }
            OR[b].S = Sb
            OR[b].b = (p ? cholsolve(Sb, kz) : J(0, 1, .))
            if (p) if (missing(OR[b].b)) _cdid_st_fail("The stayer outcome regression is rank deficient.")
            if (p) for (j = 1; j <= rows(cl); j++) OR[b].alpha[j] = OR[b].alpha[j] - OR[b].mz[j, .] * OR[b].b
            if (p) for (j = 1; j <= rows(ex); j++) ea[j] = ea[j] - em[j, .] * OR[b].b
            /* Slope scores summed over the block's cells. */
            OR[b].ell_score = J(n, p, 0)
            if (p & compute_if) {
                for (j = 1; j <= rows(cl); j++) {
                    c = cl[j]
                    k = F.cells[c, 9]
                    t = 1 + ceil(k / 2)
                    C = (D[., t - 1] :== F.cells[c, 3]) :& (D[., t] :== F.cells[c, 3])
                    Z = editmissing(Xall[., ((t - 2) * p + 1)..((t - 1) * p)], 0) :* C
                    e = ((Y[., t] - Y[., t - 1]) - (OR[b].alpha[j] :+ Z * OR[b].b)) :* C
                    OR[b].ell_score = OR[b].ell_score + (C :* e) :* (Z :- OR[b].mz[j, .])
                }
                for (j = 1; j <= rows(ex); j++) {
                    t = extra[ex[j], 1]
                    C = (D[., t - 1] :== extra[ex[j], 2]) :& (D[., t] :== extra[ex[j], 2])
                    Z = editmissing(Xall[., ((t - 2) * p + 1)..((t - 1) * p)], 0) :* C
                    e = ((Y[., t] - Y[., t - 1]) - (ea[j] :+ Z * OR[b].b)) :* C
                    OR[b].ell_score = OR[b].ell_score + (C :* e) :* (Z :- em[j, .])
                }
            }
        }
    }
    /* Cell contrasts and influence functions. */
    for (j = 1; j <= rows(av); j++) {
        c = av[j]
        k = F.cells[c, 9]
        t = 1 + ceil(k / 2)
        d = F.cells[c, 3]
        G = (D[., t - 1] :== d) :& (D[., t] :== 1 - d)
        C = (D[., t - 1] :== d) :& (D[., t] :== d)
        R = G :| C
        dy = Y[., t] - Y[., t - 1]
        pi = sum(G) / n
        q0 = sum(C) / n
        Z = (p ? Xall[., ((t - 2) * p + 1)..((t - 1) * p)] : J(n, 0, .))
        if (p) if (missing(select(Z, R)))
            _cdid_st_fail("Covariates must be nonmissing at date t - 1 for every unit at risk of switching at t.")
        Z = editmissing(Z, 0) :* R
        Xt = J(n, 1, 1), Z
        Xt = Xt :* R
        Xt[., 1] = J(n, 1, 1)
        e = dy
        if (use_or) {
            b = blocks[j]
            ga = selectindex(OR[b].cells :== c)
            e = dy - (OR[b].alpha[ga] :+ (p ? Z * OR[b].b : J(n, 1, 0)))
        }
        e = e :* R
        if (rep & use_or) {
            F.repW[c].orblock = b
            F.repW[c].alpha = OR[b].alpha[ga]
            F.repW[c].b = OR[b].b
            F.repW[c].zbar = (p ? OR[b].mz[ga, .]' : J(0, 1, .))
            if (compute_if) {
                F.repW[c].chib = (p ? cholsolve(OR[b].S :/ n, OR[b].ell_score')' : J(n, 0, .))
                F.repW[c].chia = (C :* e) :/ q0 - (p ? F.repW[c].chib * F.repW[c].zbar : J(n, 1, 0))
            }
        }
        if (use_logit) {
            _cdid_st_logit(Xt, G, R, n, eps, maxit, odds, act, prob, Hc, gam)
            Q = sum(C :* odds) / n
            if (rep) {
                F.repW[c].gam = gam
                if (compute_if) F.repW[c].ell = solveupper(Hc', solvelower(Hc, ((R :* (G - prob)) :* Xt)'))'
            }
        }
        mu1 = sum(G :* e) / sum(G)
        mu0 = (use_logit ? sum(C :* odds :* e) / sum(C :* odds) : 0)
        F.theta[c] = mu1 - mu0
        F.tau[c] = F.cells[c, 4] * F.theta[c]
        if (!compute_if) continue
        phi = (G :/ pi) :* (e :- mu1)
        if (use_logit) {
            phi = phi - (C :* odds :/ Q) :* (e :- mu0)
            ga = colsum(Xt :* (C :* act :* odds :* (e :- mu0)))' :/ sum(C :* odds)
            lam = solveupper(Hc', solvelower(Hc, ga))
            phi = phi - (R :* (G - prob)) :* (Xt * lam)
        }
        if (use_or) {
            if (m == "reg") phi = phi - (C :* e) :/ q0
            if (p) {
                zb1 = (colsum(Z :* G) :/ sum(G))'
                if (m == "reg") {
                    ga = selectindex(OR[b].cells :== c)
                    zb0 = OR[b].mz[ga, .]'
                }
                else zb0 = (colsum(Z :* (C :* odds)) :/ sum(C :* odds))'
                lam = cholsolve(OR[b].S :/ n, zb1 - zb0)
                phi = phi - OR[b].ell_score * lam
            }
        }
        if (missing(phi) | missing(F.theta[c])) _cdid_st_fail("Nonfinite cell estimate or influence function.")
        F.IF[., c] = F.cells[c, 4] :* phi
    }
    _cdid_st_aggregate(F, Gm, n, compute_if)
    return(F)
}

/* Missing outcomes (no selection model). For each requested cell and role (switchers G,
   stayers C), the role graph has a vertex for every date and an edge (s, u), s < u, whenever
   a role unit has observed outcomes at both dates; the edge moment is the (odds-weighted for
   tilted stayers) mean change of those units. Role levels relative to t - 1 are estimated by
   one joint GMM over all requested cells and roles. Outcome regressions reconstruct the stayer
   mean change and the covariance with the covariates through the least-squares flow of the
   stayer graph; centered links subtract the predicted switch-step change on edges that cross
   from before t to t or later. */
struct _cdid_stg {
    real matrix E
    real colvector dates, comp
    real scalar base_v, targ_v, connected, ncomp
}

real scalar _cdid_st_find(real colvector par, real scalar x)
{
    real scalar r, y
    r = x
    while (par[r] != r) r = par[r]
    while (par[x] != r) {
        y = par[x]
        par[x] = r
        x = y
    }
    return(r)
}

struct _cdid_stg scalar _cdid_st_graph(real matrix O, real colvector A, real scalar tb, real scalar tt)
{
    struct _cdid_stg scalar g
    real scalar T, s, u, k, j, ra, rb, nc
    real matrix E
    real colvector par, Oa, lab, map
    T = cols(O)
    E = J(T * (T - 1) / 2, 2, .)
    k = 0
    for (s = 1; s < T; s++) {
        Oa = A :& O[., s]
        if (!any(Oa)) continue
        for (u = s + 1; u <= T; u++) {
            if (any(Oa :& O[., u])) {
                k++
                E[k, .] = (s, u)
            }
        }
    }
    g.base_v = g.targ_v = g.connected = g.ncomp = 0
    if (!k) {
        g.E = J(0, 2, .)
        g.dates = g.comp = J(0, 1, .)
        return(g)
    }
    g.E = E[|1, 1 \ k, 2|]
    g.dates = uniqrows(vec(g.E))
    par = (1::T)
    for (j = 1; j <= k; j++) {
        ra = _cdid_st_find(par, g.E[j, 1])
        rb = _cdid_st_find(par, g.E[j, 2])
        if (ra != rb) par[max((ra, rb))] = min((ra, rb))
    }
    lab = J(rows(g.dates), 1, .)
    for (j = 1; j <= rows(g.dates); j++) lab[j] = _cdid_st_find(par, g.dates[j])
    map = J(T, 1, 0)
    nc = 0
    g.comp = J(rows(g.dates), 1, .)
    for (j = 1; j <= rows(g.dates); j++) {
        if (!map[lab[j]]) map[lab[j]] = ++nc
        g.comp[j] = map[lab[j]]
    }
    g.ncomp = nc
    for (j = 1; j <= rows(g.dates); j++) {
        if (g.dates[j] == tb) g.base_v = j
        if (g.dates[j] == tt) g.targ_v = j
    }
    if (g.base_v & g.targ_v) g.connected = (g.comp[g.base_v] == g.comp[g.targ_v])
    return(g)
}

/* Least-squares flow of the stayer graph from t - 1 to t (zero outside the endpoint component). */
real colvector _cdid_st_flow(struct _cdid_stg scalar g)
{
    real scalar ce, k, nv, j, a, b
    real colvector vin, pcol, ein, ell, u, varpi
    real matrix Dm, Dg
    ce = g.comp[g.base_v]
    vin = selectindex(g.comp :== ce)
    nv = rows(vin)
    pcol = J(rows(g.dates), 1, 0)
    k = 0
    for (j = 1; j <= nv; j++) if (vin[j] != g.base_v) pcol[vin[j]] = ++k
    ein = J(0, 1, .)
    for (j = 1; j <= rows(g.E); j++) {
        a = selectindex(g.dates :== g.E[j, 1])
        if (g.comp[a] == ce) ein = ein \ j
    }
    Dm = J(rows(ein), nv - 1, 0)
    for (j = 1; j <= rows(ein); j++) {
        a = pcol[selectindex(g.dates :== g.E[ein[j], 1])]
        b = pcol[selectindex(g.dates :== g.E[ein[j], 2])]
        if (a) Dm[j, a] = Dm[j, a] - 1
        if (b) Dm[j, b] = Dm[j, b] + 1
    }
    ell = J(nv - 1, 1, 0)
    ell[pcol[g.targ_v]] = 1
    Dg = cross(Dm, Dm)
    u = cholsolve(Dg, ell)
    if (missing(u)) _cdid_st_fail("The stayer graph flow is not identified.")
    varpi = J(rows(g.E), 1, 0)
    varpi[ein] = Dm * u
    return(varpi)
}

/* Observation models for missing outcomes. Nuisance cells are (role group, date pair) for the
   pair logit and the change imputation, or (role group, date) for the date logit; each cell has
   a local intercept and shares its slopes within its block (separate: the cell; within_group: the
   role group; pooled: all cells). */
struct _cdid_st_grp {
    real colvector A
    real matrix Xs, Xi, E
    real colvector pcell, icell, dcell
}

struct _cdid_st_obs {
    string scalar route, selmodel
    real scalar uses_q, uses_j, per_period
    real colvector U, rho
    struct _cdid_st_nf scalar sel, imp
    struct _cdid_st_grp colvector grp
}

/* Assigns blocks to the nuisance cells under the pooling rule. */
void _cdid_st_nf_blocks(struct _cdid_st_nf scalar f, string scalar pooling)
{
    real scalar k, b, P
    real colvector gmap
    f.cblock = J(f.ncell, 1, .)
    if (pooling == "separate") f.cblock = (1::f.ncell)
    else if (pooling == "pooled") f.cblock = J(f.ncell, 1, 1)
    else {
        gmap = J(max(f.cgrp), 1, 0)
        b = 0
        for (k = 1; k <= f.ncell; k++) {
            if (!gmap[f.cgrp[k]]) gmap[f.cgrp[k]] = ++b
            f.cblock[k] = gmap[f.cgrp[k]]
        }
    }
    f.nblock = max(f.cblock)
    f.calpha = J(f.ncell, 1, .)
    f.bbeta = J(f.nblock, f.ds, .)
    P = 0
    for (b = 1; b <= f.nblock; b++) {
        for (k = 1; k <= f.ncell; k++) if (f.cblock[k] == b) f.calpha[k] = ++P
        if (f.ds) {
            f.bbeta[b, .] = (P + 1)..(P + f.ds)
            P = P + f.ds
        }
    }
    f.P = P
    f.theta = J(P, 1, 0)
}

/* Fitting rows, responses and design of nuisance cell k. */
void _cdid_st_nf_rows(struct _cdid_st_nf scalar f, struct _cdid_st_grp colvector grp, real matrix Y,
    real matrix O, real scalar k, real scalar wls, real scalar per_period, real colvector rows,
    real colvector y, real matrix X)
{
    real colvector A, ob
    A = grp[f.cgrp[k]].A
    if (wls) {
        ob = A :& O[., f.cs[k]] :& O[., f.cu[k]]
        rows = selectindex(ob)
        y = Y[rows, f.cu[k]] - Y[rows, f.cs[k]]
        X = grp[f.cgrp[k]].Xi[rows, .]
        if (per_period) X = X :* (f.cu[k] - f.cs[k])
    }
    else {
        rows = selectindex(A :& f.U)
        if (f.cu[k]) y = O[rows, f.cs[k]] :& O[rows, f.cu[k]]
        else y = O[rows, f.cs[k]]
        X = grp[f.cgrp[k]].Xs[rows, .]
    }
}

real colvector _cdid_st_nf_eta(struct _cdid_st_nf scalar f, real scalar k, real matrix X)
{
    if (f.ds) return(f.theta[f.calpha[k]] :+ X * f.theta[f.bbeta[f.cblock[k], .]'])
    return(J(rows(X), 1, f.theta[f.calpha[k]]))
}

/* Block-wise fit (logit by IRLS with local intercepts, or least squares) and the original-unit
   influence function of all parameters (Fisher or Gram partialled on the cell means). */
void _cdid_st_nf_fit(struct _cdid_st_nf scalar f, struct _cdid_st_grp colvector grp, real matrix Y,
    real matrix O, real scalar wls, real scalar per_period, real scalar eps, real scalar maxit,
    real scalar need_if)
{
    real scalar b, k, it, dev, dev_old, conv, N, ds, nk, m_k
    real colvector cells, rows, y, eta, mu, w, z, wsum, yb, beta, sc, alpha
    real matrix X, G, xb, gvec, Xc, st, xib, Gc
    real rowvector can_increase, can_decrease, varies, min_unobs, max_unobs, min_obs, max_obs
    N = f.N
    ds = f.ds
    f.iter = J(f.nblock, 1, 0)
    f.nfit = f.nobs = J(f.ncell, 1, 0)
    if (need_if) f.IF = J(N, f.P, 0)
    for (b = 1; b <= f.nblock; b++) {
        cells = selectindex(f.cblock :== b)
        if (wls) {
            G = J(ds, ds, 0)
            gvec = J(ds, 1, 0)
            xb = J(rows(cells), ds, .)
            yb = wsum = J(rows(cells), 1, .)
            for (k = 1; k <= rows(cells); k++) {
                _cdid_st_nf_rows(f, grp, Y, O, cells[k], 1, per_period, rows, y, X)
                nk = rows(rows)
                if (nk < 1) _cdid_st_fail("An imputation cell has no unit observed at both dates.")
                f.nfit[cells[k]] = f.nobs[cells[k]] = nk
                wsum[k] = nk
                yb[k] = mean(y)
                if (ds) {
                    xb[k, .] = mean(X)
                    Xc = X :- xb[k, .]
                    G = G + cross(Xc, Xc)
                    gvec = gvec + cross(Xc, y :- yb[k])
                }
            }
            beta = J(0, 1, .)
            if (ds) {
                beta = cholsolve(G, gvec)
                if (missing(beta)) _cdid_st_fail("The change imputation regression is rank deficient.")
                f.theta[f.bbeta[b, .]'] = beta
            }
            for (k = 1; k <= rows(cells); k++)
                f.theta[f.calpha[cells[k]]] = yb[k] - (ds ? xb[k, .] * beta : 0)
        }
        else {
            dev_old = 0
            can_increase = can_decrease = J(1, ds, 1)
            varies = J(1, ds, 0)
            for (k = 1; k <= rows(cells); k++) {
                _cdid_st_nf_rows(f, grp, Y, O, cells[k], 0, 0, rows, y, X)
                if (rows(rows) < 1) _cdid_st_fail("An observation logit cell has no unit.")
                f.nfit[cells[k]] = rows(rows)
                f.nobs[cells[k]] = sum(y)
                if (f.nobs[cells[k]] == 0 | f.nobs[cells[k]] == f.nfit[cells[k]])
                    _cdid_st_fail("The observation logit has no finite root: an entire cell is observed or unobserved.")
                if (ds) {
                    /* With a local intercept, each cell can set its own cutoff.
                       A shared slope must separate in the same direction in
                       every cell, with a strict inequality on at least one row. */
                    Xc = X[selectindex(y :== 0), .]
                    min_unobs = colmin(Xc)
                    max_unobs = colmax(Xc)
                    Xc = X[selectindex(y :== 1), .]
                    min_obs = colmin(Xc)
                    max_obs = colmax(Xc)
                    can_increase = can_increase :& (max_unobs :<= min_obs)
                    can_decrease = can_decrease :& (max_obs :<= min_unobs)
                    varies = varies :| (min_unobs :!= max_unobs) :|
                        (min_obs :!= max_obs) :| (min_unobs :!= min_obs)
                }
                dev_old = dev_old + _cdid_st_bdev(y, J(rows(y), 1, .5))
            }
            if (any((can_increase :| can_decrease) :& varies))
                _cdid_st_fail("The observation logit has no finite root: separation by a selection covariate.")
            conv = 0
            for (it = 1; it <= maxit; it++) {
                G = J(ds, ds, 0)
                gvec = J(ds, 1, 0)
                xb = J(rows(cells), ds, .)
                yb = wsum = J(rows(cells), 1, .)
                for (k = 1; k <= rows(cells); k++) {
                    _cdid_st_nf_rows(f, grp, Y, O, cells[k], 0, 0, rows, y, X)
                    eta = _cdid_st_nf_eta(f, cells[k], X)
                    mu = _cdid_st_mu(eta)
                    w = mu :* (1 :- mu)
                    if (any(w :<= 0)) _cdid_st_fail("An observation logit reached probabilities of 0 or 1; the observation model is not identified.")
                    z = eta + (y - mu) :/ w
                    wsum[k] = sum(w)
                    yb[k] = sum(w :* z) / wsum[k]
                    if (ds) {
                        xb[k, .] = colsum(w :* X) :/ wsum[k]
                        Xc = X :- xb[k, .]
                        G = G + cross(Xc, w, Xc)
                        gvec = gvec + cross(Xc, w, z :- yb[k])
                    }
                }
                if (ds) {
                    beta = cholsolve(G, gvec)
                    if (missing(beta)) _cdid_st_fail("The observation logit has a rank-deficient design.")
                    f.theta[f.bbeta[b, .]'] = beta
                }
                for (k = 1; k <= rows(cells); k++)
                    f.theta[f.calpha[cells[k]]] = yb[k] - (ds ? xb[k, .] * beta : 0)
                dev = 0
                for (k = 1; k <= rows(cells); k++) {
                    _cdid_st_nf_rows(f, grp, Y, O, cells[k], 0, 0, rows, y, X)
                    dev = dev + _cdid_st_bdev(y, _cdid_st_mu(_cdid_st_nf_eta(f, cells[k], X)))
                }
                if (dev >= .) _cdid_st_fail("The observation logit did not converge.")
                if (abs(dev - dev_old) / (.1 + abs(dev)) < eps) {
                    conv = 1
                    f.iter[b] = it
                    break
                }
                dev_old = dev
            }
            if (!conv) _cdid_st_fail("The observation logit did not converge; the observed and unobserved units may be separated by the covariates.")
        }
        if (!need_if) continue
        /* Final information: cell masses, weighted covariate means, slope Gram; scores summed by unit. */
        G = J(ds, ds, 0)
        xb = J(rows(cells), ds, .)
        wsum = J(rows(cells), 1, .)
        st = J(N, ds, 0)
        sc = J(N, rows(cells), 0)
        for (k = 1; k <= rows(cells); k++) {
            _cdid_st_nf_rows(f, grp, Y, O, cells[k], wls, per_period, rows, y, X)
            eta = _cdid_st_nf_eta(f, cells[k], X)
            if (wls) {
                w = J(rows(rows), 1, 1)
                sc[rows, k] = y - eta
            }
            else {
                mu = cdid_logit_values(eta)[., 1]
                if (any(mu :<= 10 * epsilon(1)) | any(mu :>= 1 - 10 * epsilon(1)))
                    _cdid_st_fail("An observation logit reached probabilities of 0 or 1; the observation model is not identified.")
                w = mu :* (1 :- mu)
                sc[rows, k] = y - mu
            }
            wsum[k] = sum(w)
            if (ds) {
                xb[k, .] = colsum(w :* X) :/ wsum[k]
                Xc = X :- xb[k, .]
                G = G + cross(Xc, w, Xc)
                st[rows, .] = st[rows, .] + sc[rows, k] :* Xc
            }
        }
        xib = J(N, ds, 0)
        if (ds) {
            xib = cholsolve(G :/ N, st')'
            if (missing(xib)) _cdid_st_fail("Singular observation-model information.")
            f.IF[., f.bbeta[b, .]] = xib
        }
        for (k = 1; k <= rows(cells); k++) {
            m_k = wsum[k] / N
            f.IF[., f.calpha[cells[k]]] = sc[., k] :/ m_k - (ds ? xib * xb[k, .]' : J(N, 1, 0))
        }
    }
}

real colvector _cdid_st_mu(real colvector eta)
{
    real colvector mu, e
    real scalar epsp
    if (missing(eta)) _cdid_st_fail("Nonfinite observation logit index.")
    e = exp(eta)
    mu = (eta :>= 0) :* (1 :/ (1 :+ exp(-eta))) + (eta :< 0) :* (e :/ (1 :+ e))
    epsp = 10 * epsilon(1)
    mu = rowmin((rowmax((mu, J(rows(mu), 1, epsp))), J(rows(mu), 1, 1 - epsp)))
    return(mu)
}

real scalar _cdid_st_bdev(real colvector y, real colvector mu)
{
    return(2 * sum(y :* ln(1 :/ mu) + (1 :- y) :* ln(1 :/ (1 :- mu))))
}

/* Pair observation probability on all units (zero outside the role) and its channel
   contribution n x m for the coefficient matrix coef (rows of the role), scaled by 1/S. */
void _cdid_st_qpred(struct _cdid_st_obs scalar ob, real scalar g, real scalar e, real colvector q,
    real colvector k2, real matrix Fq)
{
    real colvector A
    real scalar k
    A = ob.grp[g].A
    if (ob.selmodel == "pair_logit") {
        k = ob.grp[g].pcell[e]
        q = cdid_logit_values(_cdid_st_nf_eta(ob.sel, k, ob.grp[g].Xs))[., 1] :* A
        k2 = k
        Fq = q :* (1 :- q)
    }
    else {
        k2 = (ob.grp[g].dcell[ob.grp[g].E[e, 1]] \ ob.grp[g].dcell[ob.grp[g].E[e, 2]])
        Fq = cdid_logit_values(_cdid_st_nf_eta(ob.sel, k2[1], ob.grp[g].Xs))[., 1],
             cdid_logit_values(_cdid_st_nf_eta(ob.sel, k2[2], ob.grp[g].Xs))[., 1]
        q = Fq[., 1] :* Fq[., 2] :* A
        Fq = (Fq[., 2] :* Fq[., 1] :* (1 :- Fq[., 1]), Fq[., 1] :* Fq[., 2] :* (1 :- Fq[., 2])) :* A
    }
}

real matrix _cdid_st_chan(struct _cdid_st_nf scalar f, real colvector ks, real matrix Fd, real matrix X,
    real matrix coef, real scalar S)
{
    real matrix out
    real scalar j, k
    real colvector v
    out = J(rows(coef), cols(coef), 0)
    for (j = 1; j <= rows(ks); j++) {
        k = ks[j]
        out = out + f.IF[., f.calpha[k]] * (colsum(coef :* Fd[., j]) :/ S)
        if (f.ds) out = out + f.IF[., f.bbeta[f.cblock[k], .]] * ((X' * (coef :* Fd[., j])) :/ S)
    }
    return(out)
}

/* Edge moment of role group g on its edge e: mean (kind 1, with tilt h) or covariance with Z
   (kind 2), with the direct influence and the tilt, observation and imputation channels. */
void _cdid_st_edge(struct _cdid_st_obs scalar ob, real scalar g, real scalar e, real matrix Y, real matrix O,
    real scalar kind, real colvector h, real scalar tilt, struct _cdid_st_cell scalar Wc, real matrix Z,
    real scalar need_if, real rowvector est, real matrix IF)
{
    real scalar n, s, u, S, M, k, gap
    real colvector A, obs, dy, Rw, q, qf, phi, act, dphi, j, w, Yv, ks
    real matrix Fw, Gy, Zbar, Zc, coef, Fq, Xj
    string scalar r
    n = rows(Y)
    r = ob.route
    s = ob.grp[g].E[e, 1]
    u = ob.grp[g].E[e, 2]
    A = ob.grp[g].A
    obs = A :& O[., s] :& O[., u]
    dy = editmissing(Y[., u] - Y[., s], 0) :* obs
    Rw = ((r == "baseline" | r == "weighting") ? obs : A)
    if (ob.uses_q) {
        _cdid_st_qpred(ob, g, e, q, ks, Fq)
        qf = rowmax((q, J(n, 1, 1e-6)))
        act = obs :* (q :> 1e-6)
        phi = obs :/ (ob.rho :* qf)
        dphi = -act :/ (ob.rho :* qf :* qf)
    }
    if (ob.uses_j) {
        k = ob.grp[g].icell[e]
        gap = (ob.per_period ? u - s : 1)
        Xj = ob.grp[g].Xi :* gap
        j = _cdid_st_nf_eta(ob.imp, k, Xj) :* A
    }
    if (r == "baseline") {
        w = h :* obs
        Yv = dy
    }
    else if (r == "weighting") {
        w = h :* phi
        Yv = dy
    }
    else if (r == "imputation") {
        w = h :* A
        Yv = j
    }
    else {
        w = h :* A
        Yv = (j + phi :* (dy - j)) :* A
    }
    S = sum(w)
    if (S <= 0 | S >= .) _cdid_st_fail("An observed-pair moment has a zero or nonfinite weight sum.")
    M = sum(w :* Yv) / S
    if (kind == 1) {
        est = M
        Fw = (Yv :- M) :* Rw
        Gy = Rw
    }
    else {
        Zbar = colsum(w :* Z) :/ S
        Zc = (Z :- Zbar) :* Rw
        est = (Zc' * (w :* Yv))' :/ S
        Fw = Zc :* (Yv :- M) - Rw * est
        Gy = Zc
    }
    if (!need_if) {
        IF = J(n, cols(est), 0)
        return
    }
    IF = n :* (w :* Fw) :/ S
    if (kind == 1 & tilt) {
        coef = (r == "weighting" ? phi : Rw) :* Fw
        IF = IF + Wc.ell * ((Wc.Xt' * ((Wc.act :* Wc.odds) :* coef)) :/ S)
    }
    if (ob.uses_q) {
        if (r == "weighting") coef = (h :* dphi) :* Fw
        else coef = (h :* (dy - j) :* dphi) :* Gy
        IF = IF + _cdid_st_chan(ob.sel, ks, Fq, ob.grp[g].Xs, coef, S)
    }
    if (ob.uses_j) {
        coef = (h :* (r == "imputation" ? A : (A - phi))) :* Gy
        IF = IF + _cdid_st_chan(ob.imp, k, A, Xj, coef, S)
    }
}

struct cdid_st_fit scalar cdid_static_missing(real matrix Y, real matrix D, real matrix Xall, real scalar p,
    string scalar method, real rowvector dirs, real rowvector tsel, string scalar pooling, real scalar eps,
    real scalar maxit, real scalar compute_if, real colvector tval, string scalar weighting,
    real scalar trigger, real scalar ridge, real scalar cap, string scalar route, string scalar selmodel,
    string scalar selpool, string scalar imppool, real scalar per_period, real matrix Xsel, real scalar ps,
    real colvector Usel, real colvector rho, | real scalar rep)
{
    real colvector gam
    struct _cdid_st_obs scalar ob
    real scalar g, e, kk
    real rowvector est
    real matrix Zs, IFe
    real colvector ones
    struct cdid_st_fit scalar F
    struct _cdid_stg colvector gG, gC
    struct _cdid_st_cell colvector W
    struct cdid_gmm_result scalar gm
    struct cdid_gmm_unit_gram scalar src
    real scalar n, T, J, j, c, k, t, d, use_logit, use_or, nb, b, ja, L, Lc, col, nG, nCr, Sk, need_if
    real matrix O, Gm, Z, Xt, Hc, IFk, Phi, A, H, ZC, chiCov, ZcC, chiAb, chiB, chib, Ap
    real matrix chi_zG, chi_zCr, Jg, PhiG, PhiC, tmp
    real colvector blk
    real matrix extra
    real colvector G, C, R, odds, act, prob, Mk, varpi, obs, dy, zb, cov, blocks, cl, r, step
    real colvector chi_alpha, chi_cG, chi_cC, rG, rC, av, Bp, bP, zG, zCr, hC
    n = rows(Y)
    T = cols(Y)
    O = (Y :< .)
    use_logit = (method == "ipw" | method == "dr")
    use_or = (method == "reg" | method == "dr")
    need_if = (compute_if | weighting != "identity")
    if (args() < 26) rep = 0
    F.has_rep = rep
    F.rep_q = F.rep_j = 0
    _cdid_st_targets(D, dirs, tsel, tval, F, Gm)
    J = rows(F.cells)
    gG = gC = _cdid_stg(J, 1)
    for (j = 1; j <= J; j++) {
        t = 1 + ceil(F.cells[j, 9] / 2)
        d = F.cells[j, 3]
        G = (D[., t - 1] :== d) :& (D[., t] :== 1 - d)
        C = (D[., t - 1] :== d) :& (D[., t] :== d)
        if (F.cells[j, 6] < 1) F.status[j] = "empty_switchers"
        else if (F.cells[j, 7] < 1) F.status[j] = "empty_stayers"
        else {
            gG[j] = _cdid_st_graph(O, G, t - 1, t)
            gC[j] = _cdid_st_graph(O, C, t - 1, t)
            if (!gG[j].connected) F.status[j] = "disconnected_switcher_graph"
            else if (!gC[j].connected) F.status[j] = "disconnected_stayer_graph"
        }
        F.cells[j, 8] = (F.status[j] == "ok")
    }
    av = selectindex(F.cells[., 8])
    F.theta = F.tau = J(J, 1, .)
    F.IF = J(n, J, .)
    if (rows(av)) {
        W = _cdid_st_cell(J, 1)
        /* Covariates at t - 1 and treatment logits on the full risk sets. */
        for (ja = 1; ja <= rows(av); ja++) {
            c = av[ja]
            t = 1 + ceil(F.cells[c, 9] / 2)
            d = F.cells[c, 3]
            R = (D[., t - 1] :== d)
            G = R :& (D[., t] :== 1 - d)
            Z = (p ? Xall[., ((t - 2) * p + 1)..((t - 1) * p)] : J(n, 0, .))
            if (p) if (missing(select(Z, R)))
                _cdid_st_fail("Covariates must be nonmissing at date t - 1 for every unit at risk of switching at t.")
            W[c].Z = editmissing(Z, 0) :* R
            if (use_logit) {
                Xt = (J(n, 1, 1), W[c].Z) :* R
                _cdid_st_logit(Xt, G, R, n, eps, maxit, odds, act, prob, Hc, gam)
                W[c].gam = gam
                W[c].Xt = Xt
                W[c].odds = odds
                W[c].act = act
                if (need_if) W[c].ell = solveupper(Hc', solvelower(Hc, ((R :* (G - prob)) :* Xt)'))'
            }
        }
        /* Observation groups (switchers and stayers of every available cell) and models. */
        ones = J(n, 1, 1)
        ob.route = route
        ob.selmodel = selmodel
        ob.uses_q = (route == "weighting" | route == "augmented")
        ob.uses_j = (route == "imputation" | route == "augmented")
        ob.per_period = per_period
        ob.U = Usel
        ob.rho = rho
        ob.grp = _cdid_st_grp(2 * rows(av), 1)
        for (ja = 1; ja <= rows(av); ja++) {
            c = av[ja]
            t = 1 + ceil(F.cells[c, 9] / 2)
            d = F.cells[c, 3]
            for (kk = 1; kk <= 2; kk++) {
                g = 2 * (ja - 1) + kk
                ob.grp[g].A = (D[., t - 1] :== d) :& (D[., t] :== (kk == 1 ? 1 - d : d))
                ob.grp[g].E = (kk == 1 ? gG[c].E : gC[c].E)
                ob.grp[g].Xi = W[c].Z
                if (ob.uses_q) {
                    Zs = (ps ? Xsel[., ((t - 2) * ps + 1)..((t - 1) * ps)] : J(n, 0, .))
                    if (ps) if (missing(select(Zs, ob.grp[g].A :& Usel)))
                        _cdid_st_fail("Observation-model covariates must be nonmissing at date t - 1 for the switchers and stayers.")
                    ob.grp[g].Xs = editmissing(Zs, 0) :* (ob.grp[g].A :& Usel)
                }
            }
        }
        if (ob.uses_q) {
            ob.sel.N = n
            ob.sel.ds = ps
            ob.sel.U = Usel
            ob.sel.cgrp = ob.sel.cs = ob.sel.cu = J(0, 1, .)
            for (g = 1; g <= rows(ob.grp); g++) {
                ob.grp[g].pcell = J(rows(ob.grp[g].E), 1, 0)
                ob.grp[g].dcell = J(T, 1, 0)
                for (e = 1; e <= rows(ob.grp[g].E); e++) {
                    if (selmodel == "pair_logit") {
                        ob.sel.cgrp = ob.sel.cgrp \ g
                        ob.sel.cs = ob.sel.cs \ ob.grp[g].E[e, 1]
                        ob.sel.cu = ob.sel.cu \ ob.grp[g].E[e, 2]
                        ob.grp[g].pcell[e] = rows(ob.sel.cgrp)
                    }
                    else {
                        for (kk = 1; kk <= 2; kk++) {
                            if (ob.grp[g].dcell[ob.grp[g].E[e, kk]]) continue
                            ob.sel.cgrp = ob.sel.cgrp \ g
                            ob.sel.cs = ob.sel.cs \ ob.grp[g].E[e, kk]
                            ob.sel.cu = ob.sel.cu \ 0
                            ob.grp[g].dcell[ob.grp[g].E[e, kk]] = rows(ob.sel.cgrp)
                        }
                    }
                }
            }
            ob.sel.ncell = rows(ob.sel.cgrp)
            _cdid_st_nf_blocks(ob.sel, selpool)
            _cdid_st_nf_fit(ob.sel, ob.grp, Y, O, 0, 0, eps, maxit, need_if)
        }
        if (ob.uses_j) {
            ob.imp.N = n
            ob.imp.ds = p
            ob.imp.U = J(n, 1, 1)
            ob.imp.cgrp = ob.imp.cs = ob.imp.cu = J(0, 1, .)
            for (g = 1; g <= rows(ob.grp); g++) {
                ob.grp[g].icell = J(rows(ob.grp[g].E), 1, 0)
                for (e = 1; e <= rows(ob.grp[g].E); e++) {
                    ob.imp.cgrp = ob.imp.cgrp \ g
                    ob.imp.cs = ob.imp.cs \ ob.grp[g].E[e, 1]
                    ob.imp.cu = ob.imp.cu \ ob.grp[g].E[e, 2]
                    ob.grp[g].icell[e] = rows(ob.imp.cgrp)
                }
            }
            ob.imp.ncell = rows(ob.imp.cgrp)
            _cdid_st_nf_blocks(ob.imp, imppool)
            _cdid_st_nf_fit(ob.imp, ob.grp, Y, O, 1, per_period, eps, maxit, need_if)
        }
        /* Stayer outcome regressions reconstructed from the stayer graphs. */
        if (use_or) {
            for (ja = 1; ja <= rows(av); ja++) {
                c = av[ja]
                t = 1 + ceil(F.cells[c, 9] / 2)
                d = F.cells[c, 3]
                C = (D[., t - 1] :== d) :& (D[., t] :== d)
                Z = W[c].Z
                varpi = _cdid_st_flow(gC[c])
                g = 2 * ja
                Mk = J(rows(gC[c].E), 1, 0)
                IFk = J(n, rows(gC[c].E), 0)
                cov = J(p, 1, 0)
                chiCov = J(n, p, 0)
                for (k = 1; k <= rows(gC[c].E); k++) {
                    if (varpi[k] == 0) continue
                    _cdid_st_edge(ob, g, k, Y, O, 1, ones, 0, W[c], J(n, 0, .), need_if, est, IFe)
                    Mk[k] = est
                    IFk[., k] = IFe
                    if (p) {
                        _cdid_st_edge(ob, g, k, Y, O, 2, ones, 0, W[c], Z, need_if, est, IFe)
                        cov = cov + varpi[k] :* est'
                        if (need_if) chiCov = chiCov + varpi[k] :* IFe
                    }
                }
                W[c].mC = varpi' * Mk
                W[c].chimC = IFk * varpi
                W[c].cov = cov
                W[c].chiCov = chiCov
                W[c].piC = sum(C) / n
                W[c].zbar = (p ? (colsum(Z :* C) :/ sum(C))' : J(0, 1, .))
                ZcC = (p ? (Z :- W[c].zbar') :* C : J(n, 0, .))
                W[c].S = cross(ZcC, ZcC) :/ sum(C)
            }
            _cdid_st_or_blocks(pooling, F, tval, D, blk, extra)
            if (rows(extra)) _cdid_st_fail("outcomepooling(): with missing outcomes, every cell of a block must be a requested cell; training cells outside the requested cells are not yet available.")
            blocks = blk[av]
            nb = max(blk)
            for (b = 1; b <= nb; b++) {
                cl = select(av, blocks :== b)
                if (rows(cl) == 0) continue
                Ap = J(p, p, 0)
                Bp = J(p, 1, 0)
                for (j = 1; j <= rows(cl); j++) {
                    Ap = Ap + W[cl[j]].piC :* W[cl[j]].S
                    Bp = Bp + W[cl[j]].piC :* W[cl[j]].cov
                }
                bP = J(0, 1, .)
                if (p) {
                    bP = cholsolve(Ap, Bp)
                    if (missing(bP)) _cdid_st_fail("The stayer outcome regression is rank deficient.")
                }
                chib = J(n, p, 0)
                if (p & need_if) {
                    chiAb = chiB = J(n, p, 0)
                    for (j = 1; j <= rows(cl); j++) {
                        c = cl[j]
                        t = 1 + ceil(F.cells[c, 9] / 2)
                        d = F.cells[c, 3]
                        C = (D[., t - 1] :== d) :& (D[., t] :== d)
                        ZcC = (W[c].Z :- W[c].zbar') :* C
                        chiAb = chiAb + ZcC :* (ZcC * bP) :- (W[c].piC :* W[c].S * bP)'
                        chiB = chiB + (C :- W[c].piC) * W[c].cov' + W[c].piC :* W[c].chiCov
                    }
                    chib = cholsolve(Ap, (chiB - chiAb)')'
                }
                for (j = 1; j <= rows(cl); j++) {
                    c = cl[j]
                    W[c].orblock = b
                    W[c].alpha = W[c].mC - (p ? W[c].zbar' * bP : 0)
                    W[c].b = bP
                    W[c].chib = chib
                    if (need_if) {
                        chi_alpha = W[c].chimC
                        if (p) {
                            t = 1 + ceil(F.cells[c, 9] / 2)
                            d = F.cells[c, 3]
                            C = (D[., t - 1] :== d) :& (D[., t] :== d)
                            ZcC = (W[c].Z :- W[c].zbar') :* C
                            chi_alpha = chi_alpha - (ZcC :/ W[c].piC) * bP - chib * W[c].zbar
                        }
                        W[c].chia = chi_alpha
                    }
                }
            }
        }
        /* Links of every available cell: switchers, then stayers (IPW, DR). */
        L = 0
        for (ja = 1; ja <= rows(av); ja++) {
            c = av[ja]
            L = L + rows(gG[c].E) + (use_logit ? rows(gC[c].E) : 0)
        }
        r = J(L, 1, .)
        Phi = J(n, L, 0)
        col = 0
        for (ja = 1; ja <= rows(av); ja++) {
            c = av[ja]
            t = 1 + ceil(F.cells[c, 9] / 2)
            d = F.cells[c, 3]
            G = (D[., t - 1] :== d) :& (D[., t] :== 1 - d)
            C = (D[., t - 1] :== d) :& (D[., t] :== d)
            Z = W[c].Z
            g = 2 * ja - 1
            rG = J(rows(gG[c].E), 1, .)
            PhiG = J(n, rows(gG[c].E), 0)
            for (k = 1; k <= rows(gG[c].E); k++) {
                _cdid_st_edge(ob, g, k, Y, O, 1, ones, 0, W[c], J(n, 0, .), need_if, est, IFe)
                rG[k] = est
                PhiG[., k] = IFe
            }
            if (use_logit) {
                g = 2 * ja
                rC = J(rows(gC[c].E), 1, .)
                PhiC = J(n, rows(gC[c].E), 0)
                for (k = 1; k <= rows(gC[c].E); k++) {
                    _cdid_st_edge(ob, g, k, Y, O, 1, W[c].odds, 1, W[c], J(n, 0, .), need_if, est, IFe)
                    rC[k] = est
                    PhiC[., k] = IFe
                }
            }
            if (use_or) {
                bP = W[c].b
                nG = sum(G)
                zG = (p ? (colsum(Z :* G) :/ nG)' : J(0, 1, .))
                step = (gG[c].E[., 2] :>= t) - (gG[c].E[., 1] :>= t)
                rG = rG - step :* (W[c].alpha + (p ? zG' * bP : 0))
                if (need_if) {
                    chi_cG = W[c].chia
                    if (p) {
                        chi_zG = ((Z :- zG') :* G) :/ (nG / n)
                        chi_cG = chi_cG + W[c].chib * zG + chi_zG * bP
                    }
                    PhiG = PhiG - chi_cG * step'
                }
                if (method == "dr") {
                    hC = W[c].odds
                    nCr = sum(C :* hC)
                    zCr = (p ? (colsum(Z :* (C :* hC)) :/ nCr)' : J(0, 1, .))
                    step = (gC[c].E[., 2] :>= t) - (gC[c].E[., 1] :>= t)
                    rC = rC - step :* (W[c].alpha + (p ? zCr' * bP : 0))
                    if (need_if) {
                        chi_cC = W[c].chia
                        if (p) {
                            Jg = ((Z :- zCr') :* C)' * (W[c].Xt :* (W[c].act :* hC)) :/ nCr
                            chi_zCr = ((Z :- zCr') :* (C :* hC)) :/ (nCr / n) + W[c].ell * Jg'
                            chi_cC = chi_cC + W[c].chib * zCr + chi_zCr * bP
                        }
                        PhiC = PhiC - chi_cC * step'
                    }
                }
            }
            Lc = rows(rG)
            r[|col + 1 \ col + Lc|] = rG
            if (need_if) Phi[|1, col + 1 \ n, col + Lc|] = PhiG
            col = col + Lc
            if (use_logit) {
                Lc = rows(rC)
                r[|col + 1 \ col + Lc|] = rC
                if (need_if) Phi[|1, col + 1 \ n, col + Lc|] = PhiC
                col = col + Lc
            }
        }
        /* Joint design: level of each role at every active date relative to t - 1; one gauge
           level dropped in every component that does not contain t - 1. */
        _cdid_st_design(gG, gC, av, use_logit, A, H)
        if (missing(r)) _cdid_st_fail("Nonfinite observed-pair moment.")
        if (weighting == "identity") gm = cdid_gmm_identity_map(r, A, H, cap)
        else {
            src.n_units = n
            src.rows_consumed = n
            src.normalization_count = n
            src.normalization_factor = 1 / n
            src.sample_rank_upper = n - 1
            src.mean = colsum(Phi) :/ n
            src.covariance = makesymmetric(cross(Phi, Phi) :/ n)
            src.retained_if = J(0, 0, .)
            gm = _cdid_gmm_fit_core(r, A, src.covariance, I(L), L, H, "optimal", trigger, ridge, cap, src)
        }
        for (ja = 1; ja <= rows(av); ja++) {
            c = av[ja]
            F.theta[c] = gm.targets[ja]
            F.tau[c] = F.cells[c, 4] * F.theta[c]
        }
        if (compute_if) F.IF[., av'] = (Phi * gm.target_loadings') :* F.cells[av, 4]'
        if (rep) {
            F.repW = W
            F.rep_q = ob.uses_q
            F.rep_j = ob.uses_j
            if (ob.uses_q) F.rsel = ob.sel
            if (ob.uses_j) F.rimp = ob.imp
            F.rgrp = J(rows(ob.grp), 2, .)
            for (ja = 1; ja <= rows(av); ja++) F.rgrp[|2 * ja - 1, 1 \ 2 * ja, 2|] = (av[ja], 1 \ av[ja], 2)
        }
    }
    _cdid_st_aggregate(F, Gm, n, compute_if)
    return(F)
}

void _cdid_st_design(struct _cdid_stg colvector gG, struct _cdid_stg colvector gC, real colvector av,
    real scalar use_logit, real matrix A, real matrix H)
{
    struct _cdid_stg scalar g
    real scalar ja, role, L, P, row, j, a, b, ce, k, nroles
    real colvector pc, tgt
    real matrix Ab
    nroles = (use_logit ? 2 : 1)
    L = P = 0
    for (ja = 1; ja <= rows(av); ja++) for (role = 1; role <= nroles; role++) {
        g = (role == 1 ? gG[av[ja]] : gC[av[ja]])
        L = L + rows(g.E)
        P = P + rows(g.dates) - g.ncomp
    }
    A = J(L, P, 0)
    H = J(rows(av), P, 0)
    row = 0
    k = 0
    for (ja = 1; ja <= rows(av); ja++) for (role = 1; role <= nroles; role++) {
        g = (role == 1 ? gG[av[ja]] : gC[av[ja]])
        ce = g.comp[g.base_v]
        pc = J(rows(g.dates), 1, 0)
        for (j = 1; j <= rows(g.dates); j++) {
            if (j == g.base_v) continue
            if (g.comp[j] != ce) if (j == min(selectindex(g.comp :== g.comp[j]))) continue
            pc[j] = ++k
        }
        for (j = 1; j <= rows(g.E); j++) {
            a = pc[selectindex(g.dates :== g.E[j, 1])]
            b = pc[selectindex(g.dates :== g.E[j, 2])]
            row++
            if (a) A[row, a] = A[row, a] - 1
            if (b) A[row, b] = A[row, b] + 1
        }
        H[ja, pc[g.targ_v]] = (role == 1 ? 1 : -1)
    }
}

/* Reports of the fitted nuisance models (keepnuisance): coefficients with standard errors from
   their original-unit influence functions (cluster sums, centered, divided by N^2). */
real rowvector _cdid_st_nse(real matrix IF, real colvector cl, real scalar C, real scalar n)
{
    real matrix S
    real rowvector s1, s2, v
    real scalar i
    if (cols(IF) == 0) return(J(1, 0, .))
    if (C == n) S = IF
    else {
        S = J(C, cols(IF), 0)
        for (i = 1; i <= n; i++) S[cl[i], .] = S[cl[i], .] + IF[i, .]
    }
    s1 = colsum(S)
    s2 = colsum(S :* S)
    v = (s2 - s1 :* s1 :/ C) :/ (n * n)
    return(sqrt(v :* (v :> 0)))
}

string scalar _cdid_st_num(real scalar x)
{
    return(strofreal(x, "%21.0g"))
}

string scalar _cdid_st_tlab(struct cdid_st_fit scalar F, real scalar c)
{
    return("t=" + _cdid_st_num(F.cells[c, 1]) + ", prior D=" + _cdid_st_num(F.cells[c, 3]))
}

void _cdid_st_coefrows(real matrix coef, string matrix txt, real scalar b, real colvector est,
    real rowvector se, string colvector terms, string colvector roles)
{
    real matrix add
    add = J(rows(est), 1, b), J(rows(est), 1, .), est, (cols(se) ? se' : J(rows(est), 1, .)), J(rows(est), 2, .)
    coef = coef \ add
    txt = txt \ (terms, roles)
}

void _cdid_st_zp(real matrix coef)
{
    if (!rows(coef)) return
    coef[., 5] = coef[., 3] :/ coef[., 4]
    coef[., 6] = 2 :* normal(-abs(coef[., 5]))
}

struct cdid_nuis_role scalar _cdid_st_role_imp(struct _cdid_st_nf scalar f, real colvector cl,
    real scalar C, real scalar n, real scalar compute_if, string rowvector names, string scalar pooling,
    real scalar per_period, string colvector ckey)
{
    struct cdid_nuis_role scalar R
    real scalar b, j, k, nc
    real colvector cells, est
    real rowvector se, idx
    string colvector terms, roles
    string scalar sl, lab
    R.role = "imputation"
    R.pooling = pooling
    R.weights = "unit weights"
    R.title = "Imputation regressions of the observed outcome changes"
    R.response = "Y(u) - Y(s) on the units of the role observed at both dates"
    R.sample = "switchers or stayers of each transition with both outcomes observed"
    R.blocks = J(f.nblock, 8, .)
    R.block_txt = J(f.nblock, 6, "")
    R.coef = J(0, 6, .)
    R.coef_txt = J(0, 2, "")
    sl = ""
    for (j = 1; j <= cols(names); j++) sl = sl + " + b" + strofreal(j) + "*" + (per_period ? "(u-s)*" : "") + names[j]
    for (b = 1; b <= f.nblock; b++) {
        cells = selectindex(f.cblock :== b)
        nc = rows(cells)
        idx = f.calpha[cells]'
        if (f.ds) idx = idx, f.bbeta[b, .]
        est = f.theta[idx']
        se = (compute_if ? _cdid_st_nse(f.IF[., idx], cl, C, n) : J(1, 0, .))
        terms = roles = J(0, 1, "")
        for (k = 1; k <= nc; k++) {
            terms = terms \ (nc == 1 ? "(Intercept)" : "(Intercept)[" + ckey[cells[k]] + "]")
            roles = roles \ "local_intercept"
        }
        for (j = 1; j <= f.ds; j++) {
            terms = terms \ names[j]
            roles = roles \ (nc == 1 ? "slope" : "shared_slope")
        }
        _cdid_st_coefrows(R.coef, R.coef_txt, b, est, se, terms, roles)
        R.blocks[b, .] = (nc + f.ds, nc, f.ds, sum(f.nfit[cells]), sum(f.nfit[cells]), ., nc, .)
        lab = (nc == 1 ? ckey[cells[1]] : (pooling == "pooled" ? "all transitions and roles" :
            regexr(ckey[cells[1]], ", pair .*$", "")))
        R.block_txt[b, .] = ("I" + strofreal(b), lab, "E[Y(u) - Y(s) | role, X(t-1)] = a[pair]" + sl, "", invtokens(names), "")
    }
    _cdid_st_zp(R.coef)
    return(R)
}

struct cdid_nuis_cap scalar _cdid_st_cap(struct cdid_st_fit scalar F, real scalar n, real colvector cl,
    real scalar C, real scalar clustered, string scalar method, string scalar selection, string scalar selmethod,
    string rowvector xn, string rowvector sn, string scalar selmodel, string scalar selpool,
    string scalar imppool, string scalar universe, real scalar per_period, real scalar compute_if,
    real colvector tval, string scalar outpool, real scalar missing_route)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar R
    struct _cdid_st_nf scalar f
    real colvector av, cells, est, tr
    real rowvector se, idx
    real matrix IFm
    real scalar j, c, k, b, nb, nc, g, p, role, ds
    string colvector terms, roles, ckey
    string scalar sl, lab
    av = selectindex(F.cells[., 8])
    p = cols(xn)
    cap.mode = "static"
    cap.method = method
    cap.selection = selection
    cap.selmethod = (selection == "none" ? "" : selmethod)
    cap.N = n
    cap.clustered = clustered
    cap.n_groups = C
    cap.has_intercept = 1
    cap.floor = 1e-6
    cap.maxit = 50
    cap.sel_model = selmodel
    cap.sel_pooling = selpool
    cap.sel_universe = universe
    cap.slope_names = sn
    cap.comparisons = J(0, 5, "")
    cap.n_catalog = cap.n_aux = 0
    cap.comp = ("selection", "not_applicable", "no_selection_logit", "This route fits no observation-selection logit.") \
        ("outcome", "not_applicable", "no_outcome_regression", "method(ipw) fits no outcome regression.") \
        ("treatment", "not_applicable", "no_treatment_propensity", "method(reg) fits no treatment propensity.") \
        ("imputation", "not_applicable", "no_imputation", "This route fits no actual-change imputation model.")
    cap.roles = cdid_nuis_role(1, 3)
    cap.roles[1].role = "outcome"
    cap.roles[2].role = "treatment"
    cap.roles[3].role = "imputation"
    if (!rows(av)) {
        for (k = 1; k <= 4; k++) {
            if ((k == 1 & F.rep_q) | (k == 2 & method != "ipw") | (k == 3 & method != "reg") | (k == 4 & F.rep_j))
                cap.comp[k, 2..4] = ("unavailable", "no_fitted_targets", "The requested transition cells have no available fitted target.")
        }
        return(cap)
    }
    /* Treatment logits. */
    if (method != "reg") {
        R.role = "treatment"
        R.title = "Switcher propensity logits on each transition risk set"
        R.response = "1 if the unit switches at t, among the units with D(t-1) = d"
        R.sample = "switchers and same-prior-state stayers of the transition"
        R.weights = "unit weights"
        R.pooling = "one logit per transition"
        R.blocks = J(rows(av), 8, .)
        R.block_txt = J(rows(av), 6, "")
        R.coef = J(0, 6, .)
        R.coef_txt = J(0, 2, "")
        sl = ""
        for (j = 1; j <= p; j++) sl = sl + " + b" + strofreal(j) + "*" + xn[j]
        for (j = 1; j <= rows(av); j++) {
            c = av[j]
            est = F.repW[c].gam
            se = (compute_if ? _cdid_st_nse(F.repW[c].ell, cl, C, n) : J(1, 0, .))
            terms = "(Intercept)" \ (p ? xn' : J(0, 1, ""))
            roles = "intercept" \ J(p, 1, "slope")
            _cdid_st_coefrows(R.coef, R.coef_txt, j, est, se, terms, roles)
            R.blocks[j, .] = (p + 1, 1, p, F.cells[c, 5], F.cells[c, 5], ., 1, .)
            R.block_txt[j, .] = ("T" + strofreal(j), _cdid_st_tlab(F, c), "logit Pr(switch at t | D(t-1) = d, X(t-1)) = a" + sl, "", invtokens(xn), "")
        }
        _cdid_st_zp(R.coef)
        cap.roles[2] = R
        cap.comp[3, 2..4] = ("available", "", "")
    }
    /* Outcome regressions. */
    if (method != "ipw") {
        R = cdid_nuis_role()
        R.role = "outcome"
        R.title = (missing_route ? "Stayer outcome-change regressions reconstructed from the observed date pairs" :
            "Stayer outcome-change regressions")
        R.response = "Y(t) - Y(t-1) of the stayers"
        R.sample = "same-prior-state stayers of the transitions of the block"
        R.weights = "unit weights"
        R.pooling = outpool
        nb = 0
        for (j = 1; j <= rows(av); j++) nb = max((nb, F.repW[av[j]].orblock))
        R.blocks = J(nb, 8, .)
        R.block_txt = J(nb, 6, "")
        R.coef = J(0, 6, .)
        R.coef_txt = J(0, 2, "")
        sl = ""
        for (j = 1; j <= p; j++) sl = sl + " + b" + strofreal(j) + "*" + xn[j]
        for (b = 1; b <= nb; b++) {
            cells = J(0, 1, .)
            for (j = 1; j <= rows(av); j++) if (F.repW[av[j]].orblock == b) cells = cells \ av[j]
            nc = rows(cells)
            est = J(nc, 1, .)
            IFm = J(n, 0, .)
            terms = roles = J(0, 1, "")
            for (k = 1; k <= nc; k++) {
                est[k] = F.repW[cells[k]].alpha
                if (compute_if) IFm = IFm, F.repW[cells[k]].chia
                terms = terms \ (nc == 1 ? "(Intercept)" : "(Intercept)[" + _cdid_st_tlab(F, cells[k]) + "]")
                roles = roles \ "local_intercept"
            }
            if (p) {
                est = est \ F.repW[cells[1]].b
                if (compute_if) IFm = IFm, F.repW[cells[1]].chib
                terms = terms \ xn'
                roles = roles \ J(p, 1, (nc == 1 ? "slope" : "shared_slope"))
            }
            se = (compute_if ? _cdid_st_nse(IFm, cl, C, n) : J(1, 0, .))
            _cdid_st_coefrows(R.coef, R.coef_txt, b, est, se, terms, roles)
            R.blocks[b, .] = (nc + p, nc, p, sum(F.cells[cells, 7]), sum(F.cells[cells, 7]), ., nc, .)
            lab = (nc == 1 ? _cdid_st_tlab(F, cells[1]) : "pooled transitions")
            R.block_txt[b, .] = ("O" + strofreal(b), lab, "E[Y(t) - Y(t-1) | stayer, X(t-1)] = a[transition]" + sl, "", invtokens(xn), "")
        }
        _cdid_st_zp(R.coef)
        cap.roles[1] = R
        cap.comp[2, 2..4] = ("available", "", "")
    }
    /* Observation logits. */
    if (F.rep_q) {
        f = F.rsel
        ds = f.ds
        ckey = J(f.ncell, 1, "")
        cap.cells = J(f.ncell, 12, .)
        for (k = 1; k <= f.ncell; k++) {
            g = f.cgrp[k]
            c = F.rgrp[g, 1]
            role = F.rgrp[g, 2]
            ckey[k] = _cdid_st_tlab(F, c) + ", " + (role == 1 ? "switchers" : "stayers") + ", " +
                (f.cu[k] ? "pair (" + _cdid_st_num(tval[f.cs[k]]) + "," + _cdid_st_num(tval[f.cu[k]]) + ")" :
                "date " + _cdid_st_num(tval[f.cs[k]]))
            cap.cells[k, .] = (k, f.cblock[k], (f.cu[k] ? 1 : 2), F.cells[c, 1], (f.cu[k] ? tval[f.cs[k]] : .),
                tval[(f.cu[k] ? f.cu[k] : f.cs[k])], f.nfit[k], f.nobs[k], f.nobs[k] / f.nfit[k], ., ., .)
        }
        cap.cell_key = ckey
        cap.blocks = J(f.nblock, 14, .)
        cap.block_txt = J(f.nblock, 3, "")
        cap.coef = J(0, 6, .)
        cap.coef_txt = J(0, 2, "")
        sl = ""
        for (j = 1; j <= ds; j++) sl = sl + " + b" + strofreal(j) + "*" + sn[j]
        for (b = 1; b <= f.nblock; b++) {
            cells = selectindex(f.cblock :== b)
            nc = rows(cells)
            idx = f.calpha[cells]'
            if (ds) idx = idx, f.bbeta[b, .]
            est = f.theta[idx']
            se = (compute_if ? _cdid_st_nse(f.IF[., idx], cl, C, n) : J(1, 0, .))
            terms = roles = J(0, 1, "")
            for (k = 1; k <= nc; k++) {
                terms = terms \ (nc == 1 ? "(Intercept)" : "(Intercept)[" + ckey[cells[k]] + "]")
                roles = roles \ "local_intercept"
            }
            for (j = 1; j <= ds; j++) {
                terms = terms \ sn[j]
                roles = roles \ (nc == 1 ? "slope" : "shared_slope")
            }
            _cdid_st_coefrows(cap.coef, cap.coef_txt, b, est, se, terms, roles)
            tr = uniqrows(cap.cells[cells, 4])
            cap.blocks[b, .] = (nc, nc + ds, nc, ds, sum(f.nfit[cells]), (nc == 1 ? f.nfit[cells[1]] : .), .,
                sum(f.nfit[cells]), 1, f.iter[b], nc + ds, ., nc, (rows(tr) == 1 ? tr[1] : .))
            lab = (nc == 1 ? ckey[cells[1]] : (selpool == "pooled" ? "all transitions and roles" :
                regexr(ckey[cells[1]], ", (pair|date) .*$", "")))
            cap.block_txt[b, .] = ("S" + strofreal(b), lab, (selmodel == "pair_logit" ?
                "logit Pr(both outcomes observed | role, X(t-1)) = " :
                "logit Pr(outcome observed at the date | role, X(t-1)) = ") + (nc > 1 ? "a[cell]" : "a") + sl)
        }
        _cdid_st_zp(cap.coef)
        cap.comp[1, 2..4] = ("available", "", "")
    }
    /* Imputation regressions. */
    if (F.rep_j) {
        f = F.rimp
        ckey = J(f.ncell, 1, "")
        for (k = 1; k <= f.ncell; k++) {
            g = f.cgrp[k]
            c = F.rgrp[g, 1]
            role = F.rgrp[g, 2]
            ckey[k] = _cdid_st_tlab(F, c) + ", " + (role == 1 ? "switchers" : "stayers") + ", pair (" +
                _cdid_st_num(tval[f.cs[k]]) + "," + _cdid_st_num(tval[f.cu[k]]) + ")"
        }
        cap.roles[3] = _cdid_st_role_imp(f, cl, C, n, compute_if, xn, imppool, per_period, ckey)
        cap.comp[4, 2..4] = ("available", "", "")
    }
    return(cap)
}

real scalar cdid_static_seq()
{
    pointer() scalar q
    q = findexternal("__cdid_static_seq")
    if (q == NULL) {
        q = crexternal("__cdid_static_seq")
        *q = 0
    }
    *q = *q + 1
    return(*q)
}

void cdid_static_nuis_finish(string scalar token, string scalar clustervar)
{
    struct cdid_nuis_cap scalar cap
    pointer() scalar p
    p = findexternal("__cdid_nuis_pending")
    if (p == NULL) return
    cap = *p
    rmexternal("__cdid_nuis_pending")
    cap.token = token
    cap.clustervar = clustervar
    _cdid_nuis_reg_put(token, cap)
}

/* Inference for the available cells and the aggregate. Cells use the shared IF summary,
   analytical covariance and the multiplier family; aggregate draws reuse the same signs
   but do not enter the simultaneous band. */
struct cdid_st_inf {
    real colvector av, se, analytic_se, bootstrap_se, lo, hi
    real matrix V
    real scalar crit, band_crit, agg_se, agg_analytic_se, agg_bootstrap_se, C, draws_used
    string scalar band_status
}

struct cdid_st_inf scalar cdid_static_infer(struct cdid_st_fit scalar F, real colvector unit_cluster,
    real scalar C, real scalar clustered, real scalar bstrap, real scalar cband, real scalar B,
    real scalar alpha, real scalar workspace)
{
    struct cdid_st_inf scalar O
    struct cdid_infer_budget scalar budget
    struct cdid_if_state scalar st
    struct cdid_if_summary scalar S
    struct cdid_infer_flags scalar flags, analytic
    struct cdid_multiplier_state scalar ms
    struct cdid_multiplier_result scalar mult
    struct cdid_infer_result scalar rep
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar batch
    real scalar n, Ja, cells, sr, iq
    real colvector positive, aggsum, aggdraw, ord
    n = rows(F.IF)
    O.av = selectindex(F.cells[., 8])
    Ja = rows(O.av)
    O.C = C
    O.crit = invnormal(1 - alpha / 2)
    O.band_crit = .
    O.band_status = "not_requested"
    O.agg_se = O.agg_analytic_se = O.agg_bootstrap_se = .
    O.draws_used = 0
    O.se = O.analytic_se = O.bootstrap_se = O.lo = O.hi = J(rows(F.cells), 1, .)
    O.V = J(0, 0, .)
    if (!Ja) return(O)
    cells = floor(workspace / 8)
    sr = max((1, min((B, floor(cells / max((1, 8 * C)))))))
    budget = cdid_infer_budget_create(n, C, Ja, B, n, n * Ja, sr, sr * C, cells, workspace)
    st = cdid_if_start(n, unit_cluster, C, Ja, clustered, budget)
    cdid_if_add(st, (1::n), F.IF[., O.av'])
    S = cdid_if_finish(st)
    O.V = S.covariance
    O.analytic_se[O.av] = S.analytic_se
    flags = cdid_infer_resolve_flags(bstrap, cband, "cells")
    analytic = cdid_infer_resolve_flags(0, 0, "cells")
    mult = cdid_multiplier_none()
    rep = cdid_infer_report(F.tau[O.av], S, analytic, mult, alpha)
    if (rep.status != "computed") _cdid_st_fail("Requested inference refused: " + subinstr(rep.reason, "_", " "))
    if (F.agg_avail) {
        aggsum = _cdid_st_clsum(F.agg_if, unit_cluster, C)
        aggsum = aggsum :- mean(aggsum)
        O.agg_analytic_se = sqrt(sum(aggsum :^ 2)) / n
        if (O.agg_analytic_se <= 0 | O.agg_analytic_se >= .) O.agg_analytic_se = .
    }
    positive = selectindex(S.analytic_se :> 0)
    aggdraw = J(0, 1, .)
    if (flags.draws & rows(positive)) {
        ms = cdid_multiplier_start(S, positive, B, alpha, budget)
        stream = cdid_public_rng_start(B, C, budget, "caller", J(0, 0, .))
        if (F.agg_avail) aggdraw = J(B, 1, .)
        while (stream.cursor < stream.B) {
            batch = cdid_public_rng_next(stream)
            cdid_multiplier_add(ms, batch.draw_positions, batch.cluster_positions, batch.signs)
            if (F.agg_avail) aggdraw[batch.draw_positions] = (batch.signs * aggsum) :/ sqrt(C)
        }
        cdid_public_rng_complete(stream)
        mult = cdid_multiplier_finish(ms, 0)
        O.draws_used = 1
        O.bootstrap_se[O.av[positive]] = mult.candidate_se :/ mult.target_accepted
        if (rows(aggdraw) & O.agg_analytic_se < . & !missing(aggdraw)) {
            ord = sort(aggdraw, 1)
            iq = ord[cdid_type1_index(B, .75)] - ord[cdid_type1_index(B, .25)]
            O.agg_bootstrap_se = iq / (invnormal(.75) - invnormal(.25)) * sqrt(C) / n
            if (O.agg_bootstrap_se <= 0 | O.agg_bootstrap_se >= .) O.agg_bootstrap_se = .
        }
        if (flags.simultaneous) {
            if (!all(F.cells[., 8])) O.band_status = "unavailable_target_cell"
            else if (!mult.accepted | mult.critical_value >= .) O.band_status = "unavailable_degenerate_family"
            else {
                O.band_status = "ok"
                O.band_crit = mult.critical_value
            }
        }
    }
    O.se = (flags.bootstrap ? O.bootstrap_se : O.analytic_se)
    O.agg_se = (flags.bootstrap ? O.agg_bootstrap_se : O.agg_analytic_se)
    if (O.band_status == "ok") {
        O.lo = F.tau - O.band_crit :* O.se
        O.hi = F.tau + O.band_crit :* O.se
    }
    else {
        O.lo = F.tau - O.crit :* O.se
        O.hi = F.tau + O.crit :* O.se
    }
    return(O)
}

real colvector _cdid_st_clsum(real colvector x, real colvector cl, real scalar C)
{
    real colvector out
    real scalar i
    out = J(C, 1, 0)
    for (i = 1; i <= rows(x); i++) out[cl[i]] = out[cl[i]] + x[i]
    return(out)
}

string scalar _cdid_st_tlabel(real scalar t)
{
    string scalar s
    s = strofreal(t, "%21.0g")
    s = subinstr(subinstr(s, "-", "m"), ".", "_")
    return("t" + s)
}

/* Reads the sorted balanced long panel, fits, and posts e(). */
void cdid_static_run(string scalar yv, string scalar dv, string scalar xv, string scalar timev,
    string scalar clv, string scalar touse, string scalar method, string scalar sw,
    string scalar stimes, string scalar pooling, real scalar compute_if, real scalar bstrap,
    real scalar cband, real scalar B, real scalar level, real scalar workspace, real scalar clustered,
    string scalar weighting, real scalar trigger, real scalar ridge, string scalar route, string scalar selmodel,
    string scalar selpool, string scalar imppool, real scalar per_period, string scalar sv,
    string scalar universe, string scalar rhov, real scalar keepnuis, string scalar xnames,
    string scalar snames)
{
    struct cdid_st_fit scalar F
    struct cdid_st_inf scalar I
    struct cdid_nuis_cap scalar cap
    real matrix Yl, Dl, Xl, Y, D, X, cells, b, V, Xs, Rr
    real colvector tl, tval, cl, av, status, Usel, rho
    real rowvector dirs, tsel
    real scalar n, T, p, alpha, j, k, ps
    string matrix stripe
    string colvector xn
    Yl = st_data(., yv, touse)
    Dl = st_data(., dv, touse)
    tl = st_data(., timev, touse)
    xn = tokens(xv)'
    p = rows(xn)
    Xl = (p ? st_data(., xv, touse) : J(rows(Yl), 0, .))
    ps = cols(tokens(sv))
    tval = uniqrows(tl)
    T = rows(tval)
    n = rows(Yl) / T
    if (n != floor(n) | T < 2) _cdid_st_fail("dvar() requires at least two dates and one row for every unit at every date.")
    if (any(colshape(tl, T) :!= J(n, 1, tval'))) _cdid_st_fail("dvar() requires one row for every unit at every date.")
    Y = colshape(Yl, T)
    D = colshape(Dl, T)
    X = (p ? colshape(Xl, T * p) : J(n, 0, .))
    Xs = (ps ? colshape(st_data(., sv, touse), T * ps) : J(n, 0, .))
    if (missing(D) | any((D :!= 0) :& (D :!= 1))) _cdid_st_fail("dvar() must be a binary 0/1 treatment indicator with no missing values.")
    Usel = J(n, 1, 1)
    rho = J(n, 1, 1)
    if (universe != "F") {
        Usel = rowsum(Y :< .) :>= (universe == "A" ? 2 : 1)
        if (rhov != "") {
            Rr = colshape(st_data(., rhov, touse), T)
            for (j = 1; j <= n; j++) {
                if (!Usel[j]) continue
                if (missing(Rr[j, .]) | any(Rr[j, .] :<= 0) | any(Rr[j, .] :!= Rr[j, 1]))
                    _cdid_st_fail("selrho() must be finite, positive and constant within each unit of the fitting universe.")
                rho[j] = Rr[j, 1]
            }
        }
    }
    dirs = (sw == "entry" ? 0 : (sw == "exit" ? 1 : (0, 1)))
    tsel = (stimes == "" ? J(1, 0, .) : strtoreal(tokens(stimes)))
    if (cols(tsel)) {
        for (j = 1; j <= cols(tsel); j++) {
            if (!any(tval[|2 \ T|] :== tsel[j])) {
                if (tsel[j] == tval[1]) _cdid_st_fail("switchtimes() cannot include the first date; transitions require a prior date.")
                _cdid_st_fail("switchtimes() value " + strofreal(tsel[j], "%21.0g") + " is not a transition date of the panel.")
            }
        }
        tsel = uniqrows(tsel')'
    }
    cl = (clustered ? colshape(st_data(., clv, touse), T)[., 1] : (1::n))
    alpha = 1 - level / 100
    if (missing(Y) | route != "baseline") F = cdid_static_missing(Y, D, X, p, method, dirs, tsel, pooling, 1e-12, 50,
        compute_if, tval, weighting, trigger, ridge, floor(workspace / 8), route, selmodel, selpool, imppool,
        per_period, Xs, ps, Usel, rho, keepnuis)
    else F = cdid_static_complete(Y, D, X, p, method, dirs, tsel, pooling, 1e-12, 50, compute_if, tval, keepnuis)
    if (keepnuis) {
        cap = _cdid_st_cap(F, n, cl, (clustered ? max(cl) : n), clustered, method,
            (route == "baseline" ? "none" : "conditional_mtar"),
            (route == "weighting" ? "ipw" : (route == "imputation" ? "imputation" : "dr")),
            tokens(xnames), tokens(snames), selmodel, selpool, imppool, universe, per_period, compute_if, tval,
            pooling, (missing(Y) | route != "baseline"))
        _cdid_nuis_store_pending(cap)
    }
    if (compute_if) I = cdid_static_infer(F, cl, max(cl), clustered, bstrap, cband, B, alpha, workspace)
    else {
        I.av = selectindex(F.cells[., 8])
        I.se = I.analytic_se = I.bootstrap_se = I.lo = I.hi = J(rows(F.cells), 1, .)
        I.agg_se = I.agg_analytic_se = I.agg_bootstrap_se = I.band_crit = .
        I.crit = invnormal(1 - alpha / 2)
        I.C = (clustered ? max(cl) : n)
        I.band_status = "not_requested"
        I.draws_used = 0
    }
    av = I.av
    status = J(rows(F.cells), 1, 0)
    for (k = 1; k <= rows(F.cells); k++) {
        status[k] = (F.status[k] == "empty_switchers" ? 1 : (F.status[k] == "empty_stayers" ? 2 :
            (F.status[k] == "disconnected_switcher_graph" ? 3 : (F.status[k] == "disconnected_stayer_graph" ? 4 : 0))))
    }
    cells = F.cells[., (1, 2, 3, 4, 8, 6, 7)], F.tau, F.theta, I.se, I.lo, I.hi, I.analytic_se, I.bootstrap_se, status
    stripe = J(rows(cells), 1, ""), F.cid
    st_matrix("__cdid_st_cells", cells)
    st_matrixcolstripe("__cdid_st_cells", (J(15, 1, ""), ("time" \ "tprev" \ "d" \ "sign" \ "available" \
        "switchers" \ "stayers" \ "att" \ "theta" \ "se" \ "lower" \ "upper" \ "analytic_se" \ "bootstrap_se" \ "status")))
    st_matrixrowstripe("__cdid_st_cells", stripe)
    if (rows(av)) {
        b = F.tau[av]'
        V = (compute_if ? I.V : J(rows(av), rows(av), 0))
        stripe = J(rows(av), 2, "")
        for (k = 1; k <= rows(av); k++) {
            stripe[k, 1] = (F.cells[av[k], 3] == 0 ? "entry" : "exit")
            stripe[k, 2] = _cdid_st_tlabel(F.cells[av[k], 1])
        }
        st_matrix("__cdid_st_b", b)
        st_matrixcolstripe("__cdid_st_b", stripe)
        st_matrix("__cdid_st_V", V)
        st_matrixcolstripe("__cdid_st_V", stripe)
        st_matrixrowstripe("__cdid_st_V", stripe)
    }
    st_numscalar("__cdid_st_navail", rows(av))
    st_numscalar("__cdid_st_n", n)
    st_numscalar("__cdid_st_T", T)
    st_numscalar("__cdid_st_C", I.C)
    st_numscalar("__cdid_st_agg", F.agg)
    st_numscalar("__cdid_st_agg_se", I.agg_se)
    st_numscalar("__cdid_st_agg_ase", I.agg_analytic_se)
    st_numscalar("__cdid_st_agg_bse", I.agg_bootstrap_se)
    st_numscalar("__cdid_st_crit", I.crit)
    st_numscalar("__cdid_st_bcrit", I.band_crit)
    st_numscalar("__cdid_st_draws", I.draws_used)
    st_numscalar("__cdid_st_switchers", F.total_sw)
    st_local("band_status", I.band_status)
    st_numscalar("__cdid_st_nmiss", sum(missing(Y)))
}

void cdid_static_display()
{
    real matrix M
    real scalar j, lev, agg, ase, crit, bootstrap, band, nmiss
    string scalar m, reason
    M = st_matrix("e(cells)")
    lev = st_numscalar("e(level)")
    m = st_global("e(method)")
    m = (m == "reg" ? "outcome regression (reg)" : (m == "ipw" ? "inverse probability weighting (ipw)" : "doubly robust (dr)"))
    bootstrap = (st_global("e(vce)") == "bootstrap")
    band = (st_global("e(band)") == "simultaneous")
    nmiss = st_numscalar("e(N_missing)")
    printf("\n")
    _cdid_st_line("Difference-in-differences, switching", "Number of units", st_numscalar("e(N)"))
    _cdid_st_line("Method: " + m, "Number of periods", st_numscalar("e(T)"))
    if (st_global("e(clustvar)") != "") _cdid_st_line("Switchers: " + st_global("e(switchers)"), "Number of clusters", st_numscalar("e(N_clust)"))
    else _cdid_st_line("Switchers: " + st_global("e(switchers)"), "", .)
    _cdid_st_line("Outcome pooling: " + st_global("e(outcome_pooling)"), (nmiss > 0 ? "Missing outcomes" : ""), nmiss)
    if (nmiss > 0 | st_global("e(selection)") == "conditional_mtar") {
        m = st_global("e(selection_est_method)")
        printf("{txt}Observation model: %s\n", (st_global("e(selection)") != "conditional_mtar" ? "none" :
            (m == "ipw" ? "inverse probability weighting" : (m == "imputation" ? "change imputation" :
            "augmented (weighting and imputation)"))))
    }
    if (st_global("e(clustvar)") != "") printf("{txt}(Std. err. adjusted for %s clusters in %s)\n", strofreal(st_numscalar("e(N_clust)"), "%12.0fc"), st_global("e(clustvar)"))
    printf("{txt}%s\n", 78 * "-")
    printf("{txt}%6s %-6s %9s %7s %10s %10s  %23s\n", "Time", "Switch", "Switchers", "Stayers", "ATT", "Std. err.",
        (st_global("e(vce)") == "none" ? "" : sprintf("[%g%% %s]", lev, (band ? "simult. band" : "conf. interval"))))
    printf("{txt}%s\n", 78 * "-")
    for (j = 1; j <= rows(M); j++) {
        printf("{txt}%6.0g %-6s {res}%9.0f %7.0f %10.7g %10.7g  %11.7g %11.7g\n", M[j, 1], (M[j, 3] == 0 ? "entry" : "exit"),
            M[j, 6], M[j, 7], M[j, 8], M[j, 10], M[j, 11], M[j, 12])
        if (!M[j, 5]) {
            reason = (M[j, 15] == 1 ? "no switchers" : (M[j, 15] == 2 ? "no stayers" :
                (M[j, 15] == 3 ? "the switchers' observed dates do not connect t-1 and t" :
                "the stayers' observed dates do not connect t-1 and t")))
            printf("{p 4 6 2}{txt}Unavailable: %s{p_end}\n", reason)
        }
    }
    printf("{txt}%s\n", 78 * "-")
    agg = st_numscalar("e(att)")
    ase = st_numscalar("e(se_att)")
    crit = invnormal(1 - st_numscalar("e(alpha)") / 2)
    if (agg < .) printf("{txt}%-32s {res}%10.7g %10.7g  %11.7g %11.7g\n", "Switcher-weighted average", agg, ase, agg - crit * ase, agg + crit * ase)
    else printf("{txt}Switcher-weighted average: unavailable (one or more requested cells are unavailable)\n")
    printf("{txt}%s\n", 78 * "-")
    if (st_global("e(vce)") == "none") printf("{p 0 6 2}{txt}Point estimates only (noinffunc).{p_end}\n")
    else printf("{p 0 6 2}{txt}Std. err.: %s; %s.%s The switcher-weighted average uses a pointwise interval.{p_end}\n",
        (bootstrap ? "multiplier bootstrap (interquartile range of " + strofreal(st_numscalar("e(biters)")) + " draws)" : "analytic"),
        (band ? sprintf("simultaneous %g%% band over the cells", lev) : sprintf("pointwise %g%% intervals", lev)),
        (bootstrap ? " e(V) holds the analytic covariance." : ""))
}

void _cdid_st_line(string scalar left, string scalar label, real scalar value)
{
    if (label == "") printf("{txt}%s\n", left)
    else printf("{txt}%-44s %-22s= {res}%9s\n", left, label, (value >= . ? "." : strofreal(value, "%9.0fc")))
}
end
