*! Observation-corrected cohort-by-date estimators (preset(cohort_time_fe) with selection())
*! and the nuisance-model reports of preset(cohort_time_fe).
version 19.5
mata:

void _cdid_ewr_set(struct _cdid_ewr_env scalar E)
{
    pointer() scalar p
    p = findexternal("__cdid_ewr_env")
    if (p == NULL) p = crexternal("__cdid_ewr_env")
    *p = E
}

struct _cdid_ewr_env scalar _cdid_ewr_get()
{
    pointer() scalar p
    struct _cdid_ewr_env scalar E
    p = findexternal("__cdid_ewr_env")
    if (p == NULL) {
        E.on = 0
        return(E)
    }
    E = *p
    return(E)
}

real scalar _cdid_ewr_on()
{
    pointer() scalar p
    struct _cdid_ewr_env scalar E
    p = findexternal("__cdid_ewr_env")
    if (p == NULL) return(0)
    E = *p
    return(E.on == 1)
}

void _cdid_ewr_clear()
{
    if (findexternal("__cdid_ewr_env") != NULL) rmexternal("__cdid_ewr_env")
}

/* ------------------------------------------------------------------------------------------
   Observation models of the corrected estimators. Cells c = (k, s, t) are cohort-date pairs
   with an observed change; on every unit of cohort k the corrected pair response is
     ipw         Z = S (Y_t - Y_s) / q_safe
     imputation  Z = j
     dr          Z = j + S (Y_t - Y_s - j) / q_safe
   with S = 1 when both outcomes are observed, q the fitted observation probability (pair
   logit, or the product of two date logits), q_safe = max(q, 1e-6) on the representative
   frame and Pr(U | k) rho max(q, 1e-6) on the rosters A and B, and j = a_c + (t - s) z'b the
   fitted actual change. Parameters, in order: selection blocks (local intercepts, then
   slopes), Pr(U | k) by cohort (rosters A and B), imputation groups (local intercepts, then
   slopes). IF holds N-scaled original-unit influence functions of all parameters.
   ------------------------------------------------------------------------------------------ */
struct _cdid_ewo {
    real scalar N, K, C, P, need_q, need_j, roster, datemodel, has_int, ds, icpt, p
    real colvector coh, U, rho, cell_k, cell_s, cell_t, cell_npair
    real matrix Y, O, Zs, Xs
    real scalar nsc, nsb
    real colvector sc_k, sc_s, sc_t, sc_block, sc_alpha, sc_nU, sc_pi, sb_iter, sb_rank
    real colvector cell_sc1, cell_sc2
    real matrix sb_beta
    real colvector coh_vals, n_k, nU, pi_val, pi_index
    real scalar nig
    real colvector ig_order, cell_group, cell_alpha, cell_mass, ig_nfit
    real matrix ig_beta, ig_keep, cell_meand
    real colvector theta, owner
    real matrix IF
}

void _cdid_ewo_fail(string scalar msg)
{
    errprintf("%s\n", msg)
    exit(459)
}

/* Fitting rows and binary response of selection cell c (units of its cohort in U). */
void _cdid_ewo_scrows(struct _cdid_ewo scalar f, real scalar c, real colvector rows, real colvector y)
{
    rows = selectindex((f.coh :== f.sc_k[c]) :& f.U)
    if (f.sc_s[c]) y = f.O[rows, f.sc_s[c]] :& f.O[rows, f.sc_t[c]]
    else y = f.O[rows, f.sc_t[c]]
}

real colvector _cdid_ewo_eta(struct _cdid_ewo scalar f, real scalar c, real matrix D)
{
    real colvector eta
    eta = J(rows(D), 1, 0)
    if (f.has_int) eta = eta :+ f.theta[f.sc_alpha[c]]
    if (f.ds) eta = eta + D * f.theta[f.sb_beta[f.sc_block[c], .]']
    return(eta)
}

/* Solve of a symmetric positive definite system after Jacobi scaling; missing when the
   scaled system is numerically singular. */
real matrix _cdid_ewo_solve(real matrix G, real matrix b)
{
    real colvector sc
    real matrix Gs, L
    real scalar j
    if (rows(G) == 0) return(J(0, cols(b), .))
    sc = sqrt(diagonal(G))
    if (missing(sc) | any(sc :<= 0)) return(J(rows(b), cols(b), .))
    Gs = makesymmetric(G :/ (sc * sc'))
    L = cholesky(Gs)
    if (missing(L)) return(J(rows(b), cols(b), .))
    for (j = 1; j <= rows(L); j++) if (L[j, j]^2 < 1e-14) return(J(rows(b), cols(b), .))
    return(cholsolve(Gs, b :/ sc) :/ sc)
}

/* Block IRLS of the selection logit with local intercepts and shared slopes (start at zero,
   relative deviance criterion), then the original-unit influence functions. */
void _cdid_ewo_logit(struct _cdid_ewo scalar f, real scalar b, real scalar eps, real scalar maxit,
    real scalar need_if)
{
    real colvector cells, rows, y, eta, mu, w, z, wsum, yb, beta, sc, xa
    real matrix D, G, xb, Xc, gv, sgn, st, scm, xib
    real rowvector pos, neg
    real scalar k, nc, it, dev, dev_old, conv, nfit, pm, N
    N = f.N
    cells = selectindex(f.sc_block :== b)
    nc = rows(cells)
    pm = (f.has_int ? nc : 0) + f.ds
    f.sb_rank[b] = pm
    nfit = 0
    dev_old = 0
    pos = neg = J(1, f.ds, 0)
    for (k = 1; k <= nc; k++) {
        _cdid_ewo_scrows(f, cells[k], rows, y)
        if (!rows(rows)) _cdid_ewo_fail("An observation logit cell has no fitting unit.")
        nfit = nfit + rows(rows)
        f.sc_nU[cells[k]] = rows(rows)
        f.sc_pi[cells[k]] = sum(y) / rows(rows)
        if (f.ds) {
            sgn = (2 :* y :- 1) :* f.Zs[rows, .]
            pos = pos :| (colsum(sgn :> 0) :> 0)
            neg = neg :| (colsum(sgn :< 0) :> 0)
        }
        dev_old = dev_old + _cdid_st_bdev(y, J(rows(y), 1, .5))
    }
    if (nfit < pm) _cdid_ewo_fail("An observation logit has fewer fitting rows than parameters.")
    if (!f.has_int & !f.ds) {
        f.sb_iter[b] = 0
        return
    }
    for (k = 1; k <= nc; k++) {
        _cdid_ewo_scrows(f, cells[k], rows, y)
        if (f.has_int & (all(y :== 1) | all(y :== 0)))
            _cdid_ewo_fail("An observation logit cell has every outcome pair observed or none; the logit has no finite solution. Use a coarser selpooling() or fewer selection covariates.")
    }
    if (f.ds) if (any(pos :!= neg))
        _cdid_ewo_fail("An observation logit is separated by a selection covariate; the logit has no finite solution.")
    conv = 0
    for (it = 1; it <= maxit; it++) {
        G = J(f.ds, f.ds, 0)
        gv = J(f.ds, 1, 0)
        xb = J(nc, f.ds, .)
        yb = wsum = J(nc, 1, .)
        for (k = 1; k <= nc; k++) {
            _cdid_ewo_scrows(f, cells[k], rows, y)
            D = f.Zs[rows, .]
            eta = _cdid_ewo_eta(f, cells[k], D)
            if (any(abs(eta) :> 30)) _cdid_ewo_fail("An observation logit reached fitted probabilities of 0 or 1; the observation model is not identified.")
            mu = _cdid_st_mu(eta)
            w = mu :* (1 :- mu)
            z = eta + (y - mu) :/ w
            wsum[k] = sum(w)
            yb[k] = sum(w :* z) / wsum[k]
            if (f.ds) {
                if (f.has_int) {
                    xb[k, .] = colsum(w :* D) :/ wsum[k]
                    Xc = D :- xb[k, .]
                    G = G + cross(Xc, w, Xc)
                    gv = gv + cross(Xc, w, z :- yb[k])
                }
                else {
                    G = G + cross(D, w, D)
                    gv = gv + cross(D, w, z)
                }
            }
        }
        if (f.ds) {
            beta = _cdid_ewo_solve(G, gv)
            if (missing(beta)) _cdid_ewo_fail("An observation logit has a rank-deficient design.")
            f.theta[f.sb_beta[b, .]'] = beta
        }
        if (f.has_int) for (k = 1; k <= nc; k++) f.theta[f.sc_alpha[cells[k]]] = yb[k] - (f.ds ? xb[k, .] * beta : 0)
        dev = 0
        for (k = 1; k <= nc; k++) {
            _cdid_ewo_scrows(f, cells[k], rows, y)
            eta = _cdid_ewo_eta(f, cells[k], f.Zs[rows, .])
            if (any(abs(eta) :> 30)) _cdid_ewo_fail("An observation logit reached fitted probabilities of 0 or 1; the observation model is not identified.")
            dev = dev + _cdid_st_bdev(y, _cdid_st_mu(eta))
        }
        if (dev >= .) _cdid_ewo_fail("An observation logit did not converge.")
        if (abs(dev - dev_old) / (.1 + abs(dev)) < eps) {
            conv = 1
            f.sb_iter[b] = it
            break
        }
        dev_old = dev
    }
    if (!conv) _cdid_ewo_fail("An observation logit did not converge.")
    /* Final information: cell masses, Fisher-weighted covariate means and the slope Gram. */
    G = J(f.ds, f.ds, 0)
    xb = J(nc, f.ds, .)
    wsum = J(nc, 1, .)
    st = J(N, f.ds, 0)
    scm = J(N, nc, 0)
    for (k = 1; k <= nc; k++) {
        _cdid_ewo_scrows(f, cells[k], rows, y)
        D = f.Zs[rows, .]
        mu = cdid_logit_values(_cdid_ewo_eta(f, cells[k], D))[., 1]
        if (any(mu :<= 10 * epsilon(1)) | any(mu :>= 1 - 10 * epsilon(1)))
            _cdid_ewo_fail("An observation logit reached fitted probabilities of 0 or 1; the observation model is not identified.")
        w = mu :* (1 :- mu)
        sc = y - mu
        scm[rows, k] = sc
        wsum[k] = sum(w)
        if (f.ds) {
            if (f.has_int) {
                xb[k, .] = colsum(w :* D) :/ wsum[k]
                Xc = D :- xb[k, .]
                G = G + cross(Xc, w, Xc)
                st[rows, .] = st[rows, .] + sc :* Xc
            }
            else {
                G = G + cross(D, w, D)
                st[rows, .] = st[rows, .] + sc :* D
            }
        }
    }
    if (!need_if) return
    xib = J(N, f.ds, 0)
    if (f.ds) {
        xib = (_cdid_ewo_solve(G :/ N, st'))'
        if (missing(xib)) _cdid_ewo_fail("The information of an observation logit is singular.")
        f.IF[., f.sb_beta[b, .]] = xib
    }
    if (f.has_int) {
        for (k = 1; k <= nc; k++) {
            xa = scm[., k] :/ (wsum[k] / N)
            if (f.ds) xa = xa - xib * xb[k, .]'
            f.IF[., f.sc_alpha[cells[k]]] = xa
        }
    }
}

/* Slope columns of imputation group g that are kept (R or_pooled_fit_blocks alias
   convention): a pivoted QR of the cell-demeaned design d = (t - s) z on the observed pairs
   (tolerance 1e-7) proposes the dependent columns; a proposed column is dropped only when
   z_P = c_c + z_R' A reconstructs it to roundoff on the observed pairs of every cell and on
   every unit of the cell's cohort, refused when it reconstructs on the observed pairs only,
   and kept otherwise. */
real rowvector _cdid_ewo_imp_keep(struct _cdid_ewo scalar f, real scalar g)
{
    real colvector cells, rows, ctl, trt, cc, ret, prop, pr
    real matrix X, z, md, u, s0, map, A, Q, R, rr
    real rowvector keep, cm, zb, ratio, nrm
    real scalar k, nc, c, L, j, thr, r0, nrj
    cells = select(f.ig_order, f.cell_group[f.ig_order] :== g)
    nc = rows(cells)
    keep = J(1, f.p, 1)
    if (f.p == 0) return(keep)
    X = J(0, f.p, .)
    md = J(nc, f.p, 0)
    cm = J(1, f.p, 0)
    for (k = 1; k <= nc; k++) {
        c = cells[k]
        rows = selectindex((f.coh :== f.cell_k[c]) :& f.O[., f.cell_s[c]] :& f.O[., f.cell_t[c]])
        L = f.cell_t[c] - f.cell_s[c]
        z = f.Xs[rows, .]
        if (f.icpt) {
            md[k, .] = mean(z)
            X = X \ (L :* (z :- md[k, .]))
        }
        else X = X \ (L :* z)
        cm = colmax((cm \ colmax(abs(f.Xs[selectindex(f.coh :== f.cell_k[c]), .]))))
    }
    /* Left-to-right rank proposal: a column is dependent when its residual on the retained
       columns has norm below 1e-7 times its own norm. */
    ret = J(0, 1, .)
    prop = J(0, 1, .)
    nrm = sqrt(colsum(X :^ 2))
    for (j = 1; j <= f.p; j++) {
        if (rows(ret)) {
            A = X[., ret']
            rr = X[., j] - A * qrsolve(A, X[., j])
            nrj = sqrt(sum(rr :^ 2))
        }
        else nrj = nrm[j]
        if (nrj >= 1e-7 * (nrm[j] > 0 ? nrm[j] : 1) & nrm[j] > 0) ret = ret \ j
        else prop = prop \ j
    }
    if (rows(prop) == 0) return(keep)
    A = X[., ret']
    if (rows(ret)) {
        map = qrsolve(A, X[., prop'])
        map = map + qrsolve(A, X[., prop'] - A * map)
    }
    else map = J(0, rows(prop), 0)
    thr = 64 * epsilon(1) * max((1, f.p + f.icpt))
    ctl = trt = J(rows(prop), 1, 0)
    for (k = 1; k <= nc; k++) {
        c = cells[k]
        zb = (f.icpt ? md[k, .] : J(1, f.p, 0))
        cc = (f.icpt ? (zb[prop'] - (rows(ret) ? zb[ret'] * map : J(1, rows(prop), 0)))' : J(rows(prop), 1, 0))
        for (r0 = 1; r0 <= 2; r0++) {
            if (r0 == 1) rows = selectindex((f.coh :== f.cell_k[c]) :& f.O[., f.cell_s[c]] :& f.O[., f.cell_t[c]])
            else rows = selectindex(f.coh :== f.cell_k[c])
            z = f.Xs[rows, .]
            u = z[., prop'] - (rows(ret) ? z[., ret'] * map : J(rows(z), rows(prop), 0))
            s0 = abs(z[., prop']) + (rows(ret) ? abs(z[., ret']) * abs(map) : J(rows(z), rows(prop), 0))
            ratio = colmax(abs(u :- cc') :/ _cdid_ewo_floor(s0 :+ abs(cc'), cm[prop'], smallestdouble()))
            if (missing(ratio)) ratio = editmissing(ratio, maxdouble())
            if (r0 == 1) ctl = rowmax((ctl, ratio'))
            else trt = rowmax((trt, ratio'))
        }
    }
    pr = (ctl :<= thr) :& (trt :> thr)
    if (any(pr)) _cdid_ewo_fail("The actual-change imputation design has an exact null direction on the observed pairs that does not reconstruct the predictions on the cohort; no estimates are returned.")
    for (j = 1; j <= rows(prop); j++) if (ctl[j] <= thr & trt[j] <= thr) keep[prop[j]] = 0
    return(keep)
}

/* Elementwise maximum of a matrix, a row of column floors, and a scalar floor. */
real matrix _cdid_ewo_floor(real matrix a, real rowvector b, real scalar c)
{
    real matrix m
    m = a
    m = m :* (m :>= b) + (J(rows(a), 1, 1) * b) :* (m :< b)
    m = m :* (m :>= c) + c :* (m :< c)
    return(m)
}

/* Actual-change imputation of one group: least squares of Y_t - Y_s on (t - s) z over the
   observed pairs, with local intercepts and slopes shared within the group. */
void _cdid_ewo_imp(struct _cdid_ewo scalar f, real scalar g, real scalar need_if)
{
    real colvector cells, rows, y, beta, r, yb
    real matrix d, G, gv, md, hb, ha, xib
    real scalar k, nc, c, N, L, pk
    real rowvector kc
    N = f.N
    kc = selectindex(f.ig_keep[g, .])
    pk = cols(kc)
    cells = select(f.ig_order, f.cell_group[f.ig_order] :== g)
    nc = rows(cells)
    G = J(pk, pk, 0)
    gv = J(pk, 1, 0)
    md = J(nc, pk, 0)
    yb = J(nc, 1, 0)
    f.ig_nfit[g] = 0
    for (k = 1; k <= nc; k++) {
        c = cells[k]
        rows = selectindex((f.coh :== f.cell_k[c]) :& f.O[., f.cell_s[c]] :& f.O[., f.cell_t[c]])
        if (!rows(rows)) _cdid_ewo_fail("An actual-change imputation cell has no observed pair.")
        L = f.cell_t[c] - f.cell_s[c]
        y = f.Y[rows, f.cell_t[c]] - f.Y[rows, f.cell_s[c]]
        d = L :* f.Xs[rows, kc]
        f.ig_nfit[g] = f.ig_nfit[g] + rows(rows)
        f.cell_mass[c] = rows(rows) / N
        if (f.icpt) {
            yb[k] = mean(y)
            if (pk) {
                md[k, .] = mean(d)
                G = G + cross(d :- md[k, .], d :- md[k, .])
                gv = gv + cross(d :- md[k, .], y :- yb[k])
            }
        }
        else {
            G = G + cross(d, d)
            gv = gv + cross(d, y)
        }
        if (pk) f.cell_meand[c, kc] = md[k, .]
    }
    beta = J(0, 1, .)
    if (pk) {
        beta = _cdid_ewo_solve(G, gv)
        if (missing(beta)) _cdid_ewo_fail("The actual-change imputation regression is rank deficient on the observed pairs; use a coarser imppooling() or fewer covariates.")
        f.theta[f.ig_beta[g, kc]'] = beta
    }
    if (f.icpt) for (k = 1; k <= nc; k++) f.theta[f.cell_alpha[cells[k]]] = yb[k] - (pk ? md[k, .] * beta : 0)
    if (!need_if) return
    hb = J(N, pk, 0)
    ha = J(N, nc, 0)
    for (k = 1; k <= nc; k++) {
        c = cells[k]
        rows = selectindex((f.coh :== f.cell_k[c]) :& f.O[., f.cell_s[c]] :& f.O[., f.cell_t[c]])
        L = f.cell_t[c] - f.cell_s[c]
        y = f.Y[rows, f.cell_t[c]] - f.Y[rows, f.cell_s[c]]
        d = L :* f.Xs[rows, kc]
        r = y :- (f.icpt ? f.theta[f.cell_alpha[c]] : 0)
        if (pk) r = r - d * beta
        ha[rows, k] = r
        if (pk) hb[rows, .] = hb[rows, .] + r :* (f.icpt ? d :- md[k, .] : d)
    }
    xib = J(N, pk, 0)
    if (pk) {
        xib = N :* (_cdid_ewo_solve(G, hb'))'
        f.IF[., f.ig_beta[g, kc]] = xib
    }
    if (f.icpt) for (k = 1; k <= nc; k++) {
        f.IF[., f.cell_alpha[cells[k]]] = ha[., k] :/ f.cell_mass[cells[k]] - (pk ? xib * md[k, .]' : J(N, 1, 0))
    }
}

/* cells: C x 3 (cohort value, s, t) in the order of the fitted cohort-date pairs. Zs and Xs
   hold the slope columns of the selection and imputation designs. */
struct _cdid_ewo scalar _cdid_ewo_fit(real matrix Y, real matrix O, real colvector coh, real matrix cells,
    string scalar method, string scalar selmodel, string scalar selpool, string scalar imppool,
    real matrix Zs, real scalar sel_int, real colvector U, real colvector rho, string scalar universe,
    real matrix Xs, real scalar imp_int, real scalar eps, real scalar maxit, real scalar need_if)
{
    struct _cdid_ewo scalar f
    real matrix raw, key
    real colvector ord, kk, nobs
    real scalar i, b, P, c, nb, m
    f.N = rows(Y)
    f.K = cols(Y)
    f.Y = Y
    f.O = O
    f.coh = coh
    f.C = rows(cells)
    f.cell_k = cells[., 1]
    f.cell_s = cells[., 2]
    f.cell_t = cells[., 3]
    f.need_q = (method == "ipw" | method == "dr")
    f.need_j = (method == "imputation" | method == "dr")
    f.datemodel = (selmodel == "independent_period_logit")
    f.coh_vals = uniqrows(f.cell_k)
    f.n_k = J(rows(f.coh_vals), 1, 0)
    for (i = 1; i <= rows(f.coh_vals); i++) f.n_k[i] = sum(coh :== f.coh_vals[i])
    if (any(f.n_k :< 1)) _cdid_ewo_fail("A cohort of a fitted pair has no unit.")
    f.cell_npair = J(f.C, 1, 0)
    for (i = 1; i <= f.C; i++) {
        f.cell_npair[i] = sum((coh :== f.cell_k[i]) :& O[., f.cell_s[i]] :& O[., f.cell_t[i]])
        if (f.cell_npair[i] < 1) _cdid_ewo_fail("A fitted cohort-date pair has no observed change.")
    }
    f.U = J(f.N, 1, 1)
    f.rho = J(f.N, 1, 1)
    f.roster = 0
    f.has_int = f.ds = f.icpt = f.p = 0
    f.nsc = f.nsb = f.nig = 0
    P = 0
    f.theta = J(0, 1, .)
    f.owner = J(0, 1, .)
    f.pi_index = J(0, 1, .)
    f.pi_val = J(0, 1, .)
    f.nU = J(rows(f.coh_vals), 1, 0)
    if (f.need_q) {
        f.U = (U :!= 0)
        f.rho = rho
        if (universe == "F") {
            if (!all(f.U) | any(rho :!= 1)) _cdid_ewo_fail("The representative roster U = F requires every unit in U and rho = 1.")
        }
        else {
            if (f.datemodel) _cdid_ewo_fail("selmodel(independent_period_logit) is available only with seluniverse(F); conditioning on A or B removes the independence of the observation dates.")
            if (method != "ipw") _cdid_ewo_fail("Fitting rosters A and B are available only with the pair-logit IPW observation correction, selmethod(ipw).")
            nobs = rowsum(O)
            if (any((nobs :>= 2) :& !f.U)) _cdid_ewo_fail("The fitting roster U must contain every unit with at least two observed outcomes.")
            if (any(missing(select(rho, f.U))) | any(select(rho, f.U) :<= 0)) _cdid_ewo_fail("The known ratio rho must be finite and positive on U.")
            f.roster = 1
        }
        f.has_int = (sel_int != 0)
        f.ds = cols(Zs)
        f.Zs = Zs
        if (f.ds) {
            f.Zs = editmissing(Zs, 0) :* f.U
            if (missing(select(Zs, f.U))) _cdid_ewo_fail("The selection covariates must be nonmissing on the fitting roster U.")
        }
        /* Selection cells: pairs (k, s, t) or dates (k, d), sorted by cohort then dates. */
        if (!f.datemodel) {
            ord = order(cells, (1, 2, 3))
            f.nsc = f.C
            f.sc_k = cells[ord, 1]
            f.sc_s = cells[ord, 2]
            f.sc_t = cells[ord, 3]
            f.cell_sc1 = J(f.C, 1, 0)
            f.cell_sc1[ord] = (1::f.C)
            f.cell_sc2 = J(f.C, 1, 0)
        }
        else {
            raw = (cells[., 1], cells[., 2]) \ (cells[., 1], cells[., 3])
            key = uniqrows(raw)
            f.nsc = rows(key)
            f.sc_k = key[., 1]
            f.sc_s = J(f.nsc, 1, 0)
            f.sc_t = key[., 2]
            f.cell_sc1 = f.cell_sc2 = J(f.C, 1, 0)
            for (i = 1; i <= f.C; i++) {
                f.cell_sc1[i] = selectindex((key[., 1] :== cells[i, 1]) :& (key[., 2] :== cells[i, 2]))
                f.cell_sc2[i] = selectindex((key[., 1] :== cells[i, 1]) :& (key[., 2] :== cells[i, 3]))
            }
        }
        if (selpool == "separate") f.sc_block = (1::f.nsc)
        else if (selpool == "within_cohort") {
            f.sc_block = J(f.nsc, 1, 0)
            for (i = 1; i <= f.nsc; i++) f.sc_block[i] = selectindex(f.coh_vals :== f.sc_k[i])
            kk = uniqrows(f.sc_block)
            for (i = 1; i <= f.nsc; i++) f.sc_block[i] = selectindex(kk :== f.sc_block[i])
        }
        else f.sc_block = J(f.nsc, 1, 1)
        f.nsb = max(f.sc_block)
        f.sc_alpha = J(f.nsc, 1, 0)
        f.sb_beta = J(f.nsb, f.ds, 0)
        f.owner = J(0, 1, .)
        for (b = 1; b <= f.nsb; b++) {
            if (f.has_int) for (i = 1; i <= f.nsc; i++) if (f.sc_block[i] == b) {
                f.sc_alpha[i] = ++P
                f.owner = f.owner \ b
            }
            for (i = 1; i <= f.ds; i++) {
                f.sb_beta[b, i] = ++P
                f.owner = f.owner \ b
            }
        }
        f.sc_nU = f.sc_pi = J(f.nsc, 1, .)
        f.sb_iter = f.sb_rank = J(f.nsb, 1, .)
        if (f.roster) {
            for (i = 1; i <= rows(f.coh_vals); i++) f.nU[i] = sum((coh :== f.coh_vals[i]) :& f.U)
            f.pi_val = f.nU :/ f.n_k
            f.pi_index = J(rows(f.coh_vals), 1, 0)
            nb = f.nsb
            for (i = 1; i <= rows(f.coh_vals); i++) {
                f.pi_index[i] = ++P
                f.owner = f.owner \ (++nb)
            }
        }
        else {
            for (i = 1; i <= rows(f.coh_vals); i++) f.nU[i] = sum(coh :== f.coh_vals[i])
        }
    }
    if (f.need_j) {
        f.icpt = (imp_int != 0)
        f.p = cols(Xs)
        f.Xs = Xs
        if (f.p) if (missing(Xs)) _cdid_ewo_fail("The covariates of the actual-change imputation must be nonmissing for every unit.")
        if (!f.icpt & !f.p) _cdid_ewo_fail("The actual-change imputation requires an intercept or at least one covariate.")
        ord = order((f.cell_k, f.cell_t - f.cell_s, f.cell_s, f.cell_t), (1, 2, 3, 4))
        f.ig_order = ord
        f.cell_group = J(f.C, 1, 0)
        if (imppool == "separate") f.cell_group[ord] = (1::f.C)
        else if (imppool == "within_cohort") {
            for (i = 1; i <= f.C; i++) f.cell_group[i] = selectindex(f.coh_vals :== f.cell_k[i])
        }
        else f.cell_group = J(f.C, 1, 1)
        f.nig = max(f.cell_group)
        f.cell_alpha = J(f.C, 1, 0)
        f.ig_beta = J(f.nig, f.p, 0)
        f.ig_keep = J(f.nig, f.p, 1)
        for (b = 1; b <= f.nig; b++) f.ig_keep[b, .] = _cdid_ewo_imp_keep(f, b)
        f.cell_mass = J(f.C, 1, .)
        f.cell_meand = J(f.C, f.p, 0)
        f.ig_nfit = J(f.nig, 1, 0)
        nb = (rows(f.owner) ? max(f.owner) : 0)
        for (b = 1; b <= f.nig; b++) {
            nb++
            if (f.icpt) for (i = 1; i <= f.C; i++) {
                c = ord[i]
                if (f.cell_group[c] == b) {
                    f.cell_alpha[c] = ++P
                    f.owner = f.owner \ nb
                }
            }
            for (i = 1; i <= f.p; i++) {
                if (!f.ig_keep[b, i]) continue
                f.ig_beta[b, i] = ++P
                f.owner = f.owner \ nb
            }
        }
    }
    f.P = P
    f.theta = J(P, 1, 0)
    f.IF = (need_if ? J(f.N, P, 0) : J(f.N, 0, .))
    if (f.need_q) {
        for (b = 1; b <= f.nsb; b++) _cdid_ewo_logit(f, b, eps, maxit, need_if)
        if (f.roster) {
            for (i = 1; i <= rows(f.coh_vals); i++) {
                f.theta[f.pi_index[i]] = f.pi_val[i]
                if (need_if) {
                    m = f.n_k[i]
                    f.IF[., f.pi_index[i]] = (coh :== f.coh_vals[i]) :* ((f.N / m) :* (f.U :- f.pi_val[i]))
                }
            }
        }
    }
    if (f.need_j) for (b = 1; b <= f.nig; b++) _cdid_ewo_imp(f, b, need_if)
    if (missing(f.theta)) _cdid_ewo_fail("The observation models returned nonfinite coefficients.")
    return(f)
}

/* Corrected pair response of fitted cell c on every unit of its cohort (rows), with its
   gradient in all observation parameters when requested. */
void _cdid_ewo_read(struct _cdid_ewo scalar f, real scalar c, real colvector rows, real colvector val,
    real matrix Gr, real scalar grad)
{
    real colvector S, obs, ro, dy, j, q, qf, scl, den, e, act, fac, qa, oth, eta
    real matrix Dz, Bj, dq, D, Breg
    real scalar n, s, t, L, kix, a, sc, idx
    real rowvector kc
    rows = selectindex(f.coh :== f.cell_k[c])
    n = rows(rows)
    s = f.cell_s[c]
    t = f.cell_t[c]
    L = t - s
    S = f.O[rows, s] :& f.O[rows, t]
    obs = selectindex(S)
    dy = J(n, 1, 0)
    if (rows(obs)) dy[obs] = f.Y[rows[obs], t] - f.Y[rows[obs], s]
    Gr = (grad ? J(n, f.P, 0) : J(n, 0, .))
    j = J(n, 1, 0)
    if (f.need_j) {
        Dz = (f.p ? L :* f.Xs[rows, .] : J(n, 0, .))
        j = J(n, 1, (f.icpt ? f.theta[f.cell_alpha[c]] : 0))
        if (f.p) {
            kc = selectindex(f.ig_keep[f.cell_group[c], .])
            if (cols(kc)) j = j + Dz[., kc] * f.theta[f.ig_beta[f.cell_group[c], kc]']
        }
        Bj = (f.icpt ? (J(n, 1, 1), Dz) : Dz)
    }
    val = j
    den = J(0, 1, .)
    if (f.need_q & rows(obs)) {
        ro = rows[obs]
        D = f.Zs[ro, .]
        Breg = (f.has_int ? (J(rows(ro), 1, 1), D) : D)
        if (!f.datemodel) {
            sc = f.cell_sc1[c]
            q = (f.has_int | f.ds ? cdid_logit_values(_cdid_ewo_eta(f, sc, D))[., 1] : J(rows(ro), 1, .5))
            if (grad) {
                dq = J(rows(ro), f.P, 0)
                idx = 0
                if (f.has_int) dq[., f.sc_alpha[sc]] = q :* (1 :- q)
                if (f.ds) dq[., f.sb_beta[f.sc_block[sc], .]] = (q :* (1 :- q)) :* D
            }
        }
        else {
            qa = (f.has_int | f.ds ? cdid_logit_values(_cdid_ewo_eta(f, f.cell_sc1[c], D))[., 1] : J(rows(ro), 1, .5))
            oth = (f.has_int | f.ds ? cdid_logit_values(_cdid_ewo_eta(f, f.cell_sc2[c], D))[., 1] : J(rows(ro), 1, .5))
            q = qa :* oth
            if (grad) {
                dq = J(rows(ro), f.P, 0)
                for (a = 1; a <= 2; a++) {
                    sc = (a == 1 ? f.cell_sc1[c] : f.cell_sc2[c])
                    eta = (a == 1 ? qa : oth)
                    e = (a == 1 ? oth : qa) :* eta :* (1 :- eta)
                    if (f.has_int) dq[., f.sc_alpha[sc]] = dq[., f.sc_alpha[sc]] + e
                    if (f.ds) dq[., f.sb_beta[f.sc_block[sc], .]] = dq[., f.sb_beta[f.sc_block[sc], .]] + e :* D
                }
            }
        }
        qf = rowmax((q, J(rows(q), 1, 1e-6)))
        kix = selectindex(f.coh_vals :== f.cell_k[c])
        scl = (f.roster ? f.pi_val[kix] :* f.rho[ro] : J(rows(ro), 1, 1))
        den = scl :* qf
        e = dy[obs] - j[obs]
        val[obs] = j[obs] + e :/ den
        if (grad) {
            act = (q :> 1e-6)
            Gr[obs, .] = Gr[obs, .] + ((-e :* act) :/ (scl :* qf :* qf)) :* dq
        }
    }
    if (missing(val)) _cdid_ewo_fail("The corrected pair responses are not finite.")
    if (!grad) return
    if (f.roster) {
        kix = selectindex(f.coh_vals :== f.cell_k[c])
        Gr[., f.pi_index[kix]] = -val :/ f.pi_val[kix]
    }
    if (f.need_j) {
        fac = J(n, 1, 1)
        if (rows(den)) fac[obs] = 1 :- 1 :/ den
        if (f.icpt) Gr[., f.cell_alpha[c]] = Gr[., f.cell_alpha[c]] + fac
        if (f.p) {
            kc = selectindex(f.ig_keep[f.cell_group[c], .])
            if (cols(kc)) Gr[., f.ig_beta[f.cell_group[c], kc]] = Gr[., f.ig_beta[f.cell_group[c], kc]] + fac :* Dz[., kc]
        }
    }
}

/* ------------------------------------------------------------------------------------------
   Pair catalog and the corrected least-squares projection on linked slots.
   ------------------------------------------------------------------------------------------ */
struct _cdid_ewc_cat {
    real colvector groups, nunits
    real matrix cells, observed
    real scalar K, omitted
}

/* Every cohort-date pair (s < t) observed for at least one unit of the cohort, by cohort,
   then t, then s. cells: group index, s, t, number of units observing both dates. */
struct _cdid_ewc_cat scalar _cdid_ewc_catalog(real matrix O, real colvector coh)
{
    struct _cdid_ewc_cat scalar C
    real matrix cnt
    real colvector u
    real scalar gi, s, t, K
    K = cols(O)
    C.K = K
    C.groups = uniqrows(coh)
    C.cells = J(0, 4, .)
    C.observed = J(K, rows(C.groups), 0)
    C.nunits = J(rows(C.groups), 1, 0)
    for (gi = 1; gi <= rows(C.groups); gi++) {
        u = selectindex(coh :== C.groups[gi])
        C.nunits[gi] = rows(u)
        cnt = cross(O[u, .], O[u, .])
        C.observed[., gi] = (diagonal(cnt) :> 0)
        for (t = 2; t <= K; t++) for (s = 1; s < t; s++) if (cnt[s, t] > 0) C.cells = C.cells \ (gi, s, t, cnt[s, t])
    }
    C.omitted = rows(C.groups) * K * (K - 1) / 2 - rows(C.cells)
    return(C)
}

struct _cdid_ewc_geo {
    real scalar A
    real colvector use, left, right, keep_slots, estimable, comp
    real matrix slot
}

/* Slots linked by the admitted pairs; one gauge (lowest slot) per component that carries a
   requested contrast. */
struct _cdid_ewc_geo scalar _cdid_ewc_geometry(struct _cdid_ewc_cat scalar C, real matrix slot,
    real colvector use, real matrix req)
{
    struct _cdid_ewc_geo scalar Gm
    real matrix link
    real colvector left, right, used, act, first
    real scalar i, A
    A = max(slot)
    Gm.A = A
    Gm.slot = slot
    left = right = J(rows(use), 1, .)
    for (i = 1; i <= rows(use); i++) {
        left[i] = slot[C.cells[use[i], 2], C.cells[use[i], 1]]
        right[i] = slot[C.cells[use[i], 3], C.cells[use[i], 1]]
    }
    link = J(A, A, 0)
    for (i = 1; i <= rows(use); i++) link[left[i], right[i]] = link[right[i], left[i]] = 1
    Gm.comp = _cdid_ew_components(link)
    Gm.estimable = J(rows(req), 1, 0)
    for (i = 1; i <= rows(req); i++) Gm.estimable[i] = (Gm.comp[req[i, 1]] == Gm.comp[req[i, 2]])
    used = (rows(req) ? uniqrows(select(Gm.comp[req[., 1]], Gm.estimable)) : J(0, 1, .))
    act = J(rows(use), 1, 0)
    for (i = 1; i <= rows(use); i++) act[i] = (any(used :== Gm.comp[left[i]]) & left[i] != right[i])
    act = selectindex(act)
    Gm.use = use[act]
    Gm.left = left[act]
    Gm.right = right[act]
    first = J(A, 1, 0)
    for (i = 1; i <= A; i++) first[i] = (i == 1 ? 1 : !any(Gm.comp[1..i - 1] :== Gm.comp[i]))
    Gm.keep_slots = J(0, 1, .)
    for (i = 1; i <= A; i++) if (!first[i] & any(used :== Gm.comp[i])) Gm.keep_slots = Gm.keep_slots \ i
    return(Gm)
}

struct _cdid_ewc_proj {
    real matrix beta, score, Gs
    real colvector keep, sc
    real scalar p, q, A
}

/* Corrected least squares of the pair responses on W (b[slot t] - b[slot s]) with the full
   cohort feature Gram on every admitted pair; scores add the generated-response channel. */
struct _cdid_ewc_proj scalar _cdid_ewc_project(struct _cdid_ewc_cat scalar C, real matrix W,
    real colvector coh, struct _cdid_ewc_geo scalar Gm, struct _cdid_ewo scalar ob, real scalar scores)
{
    struct _cdid_ewc_proj scalar R
    real matrix gram, rhs, load, FG, Gr, coords, vals, adj, out, Wg
    real colvector rows, val, sv, tv, e, rr, gu, gi
    real scalar n, q, A, p, i, g, cell, k, a
    n = rows(W)
    q = cols(W)
    A = Gm.A
    p = A * q
    R.p = p
    R.q = q
    R.A = A
    coords = J(A, q, .)
    for (a = 1; a <= q; a++) coords[., a] = (1::A) :+ (a - 1) * A
    R.keep = (rows(Gm.keep_slots) ? vec(coords[Gm.keep_slots, .]) : J(0, 1, .))
    R.beta = J(A, q, 0)
    R.score = (scores ? J(n, p, 0) : J(n, 0, .))
    if (!rows(R.keep)) return(R)
    gram = J(p, p, 0)
    rhs = J(p, 1, 0)
    load = (scores ? J(p, ob.P, 0) : J(p, 0, .))
    for (i = 1; i <= rows(Gm.use); i++) {
        cell = Gm.use[i]
        g = C.cells[cell, 1]
        sv = coords[Gm.left[i], .]'
        tv = coords[Gm.right[i], .]'
        _cdid_ewo_read(ob, _cdid_ewc_fitpos(ob, C, cell), rows, val, Gr, scores)
        Wg = W[rows, .]
        FG = cross(Wg, Wg)
        gram[tv, tv] = gram[tv, tv] + FG
        gram[sv, sv] = gram[sv, sv] + FG
        gram[tv, sv] = gram[tv, sv] - FG
        gram[sv, tv] = gram[sv, tv] - FG
        rhs[tv] = rhs[tv] + cross(Wg, val)
        rhs[sv] = rhs[sv] - cross(Wg, val)
        if (scores & ob.P) {
            load[tv, .] = load[tv, .] + cross(Wg, Gr)
            load[sv, .] = load[sv, .] - cross(Wg, Gr)
        }
    }
    R.Gs = gram[R.keep, R.keep]
    R.sc = sqrt(diagonal(R.Gs))
    if (missing(R.sc) | any(R.sc :<= 0)) _cdid_ewc_rankfail()
    R.Gs = makesymmetric(R.Gs :/ (R.sc * R.sc'))
    if (missing(cholesky(R.Gs))) _cdid_ewc_rankfail()
    if (1 / cond(R.Gs) <= epsilon(1)) _cdid_ewc_rankfail()
    rr = J(p, 1, 0)
    rr[R.keep] = cholsolve(R.Gs, rhs[R.keep] :/ R.sc) :/ R.sc
    R.beta = colshape(rr, A)'
    if (!scores) return(R)
    adj = (ob.P ? cholsolve(R.Gs, load[R.keep, .] :/ R.sc) :/ R.sc : J(rows(R.keep), 0, .))
    vals = J(n, p, 0)
    for (i = 1; i <= rows(Gm.use); i++) {
        cell = Gm.use[i]
        sv = coords[Gm.left[i], .]'
        tv = coords[Gm.right[i], .]'
        _cdid_ewo_read(ob, _cdid_ewc_fitpos(ob, C, cell), rows, val, Gr, 0)
        Wg = W[rows, .]
        e = n :* (val - Wg * (R.beta[Gm.right[i], .] - R.beta[Gm.left[i], .])')
        vals[rows, tv'] = vals[rows, tv'] + Wg :* e
        vals[rows, sv'] = vals[rows, sv'] - Wg :* e
    }
    out = (cholsolve(R.Gs, vals[., R.keep]' :/ R.sc) :/ R.sc)'
    if (ob.P) out = out + ob.IF * adj'
    R.score[., R.keep'] = out
    return(R)
}

void _cdid_ewc_rankfail()
{
    errprintf("the observation-corrected cohort-by-date projection is rank deficient after fixing one level per connected set of slots; no estimates are returned\n")
    exit(459)
}

/* Position of catalog cell `cell` among the fitted cells of the observation models. */
real scalar _cdid_ewc_fitpos(struct _cdid_ewo scalar ob, struct _cdid_ewc_cat scalar C, real scalar cell)
{
    real colvector h
    h = selectindex((ob.cell_k :== C.groups[C.cells[cell, 1]]) :& (ob.cell_s :== C.cells[cell, 2]) :& (ob.cell_t :== C.cells[cell, 3]))
    if (rows(h) != 1) _cdid_ewo_fail("A fitted cohort-date pair has no observation model.")
    return(h[1])
}

/* Features: intercept and covariates centered at the frame mean and divided by their largest
   absolute deviation. */
real matrix _cdid_ewc_features(real matrix Z, real colvector center, real colvector amp)
{
    real matrix cen
    if (!cols(Z)) {
        center = amp = J(0, 1, .)
        return(J(rows(Z), 1, 1))
    }
    center = mean(Z)'
    cen = Z :- center'
    amp = colmax(abs(cen))'
    amp = amp + (amp :== 0)
    return((J(rows(Z), 1, 1), cen :/ amp'))
}

/* Corrected regression of the cohort-by-date fixed-effects model. */
struct cdid_etwfe_fit scalar cdid_etwfe_corrected_reg(real matrix Y, real matrix O, real colvector coh,
    real colvector cal, real scalar never, real matrix Zc, struct _cdid_ewc_set scalar S,
    real scalar compute_if)
{
    struct cdid_etwfe_fit scalar out
    struct _cdid_ewc_cat scalar C
    struct _cdid_ewc_geo scalar Gm
    struct _cdid_ewo scalar ob
    struct _cdid_ewc_proj scalar Pj
    struct _cdid_ewr_env scalar E
    real matrix slot, nodes, W, tx, diffm, Cm, V, cellsfit
    real colvector tdate, tgroup, center, amp, ids, att, status, u, v, fitc, cohg
    real scalar K, n, nG, gi, t, g, i, a, k, q, keepnuis
    K = cols(O)
    n = rows(O)
    E = _cdid_ewr_get()
    keepnuis = E.on
    C = _cdid_ewc_catalog(O, coh)
    nG = rows(C.groups)
    slot = J(K, nG, .)
    tdate = tgroup = J(0, 1, .)
    for (gi = 1; gi <= nG; gi++) {
        g = C.groups[gi]
        for (t = 1; t <= K; t++) {
            if (g == 0) slot[t, gi] = t
            else if (never ? t != g - 1 : t >= g) {
                slot[t, gi] = (gi - 1) * K + t
                tdate = tdate \ t
                tgroup = tgroup \ gi
            }
            else slot[t, gi] = t
        }
    }
    nodes = J(rows(tdate), 2, .)
    for (i = 1; i <= rows(tdate); i++) nodes[i, .] = (slot[tdate[i], tgroup[i]], tdate[i])
    Gm = _cdid_ewc_geometry(C, slot, (1::rows(C.cells)), nodes)
    W = _cdid_ewc_features(Zc, center, amp)
    q = cols(W)
    fitc = Gm.use
    if (rows(fitc)) {
        fitc = sort(fitc, 1)
        cellsfit = (C.groups[C.cells[fitc, 1]], C.cells[fitc, 2..3])
        ob = _cdid_ewo_fit(Y, O, coh, cellsfit, S.selmethod, S.selmodel, S.selpool, S.imppool, S.Zs, S.sel_int,
            S.U, S.rho, S.universe, S.Xs, S.imp_int, S.eps, S.maxit, (compute_if | keepnuis))
    }
    else {
        ob.P = 0
        ob.IF = J(n, 0, .)
    }
    Pj = _cdid_ewc_project(C, W, coh, Gm, ob, (compute_if | keepnuis))
    tx = J(nG, q, .)
    for (gi = 1; gi <= nG; gi++) tx[gi, .] = mean(W[selectindex(coh :== C.groups[gi]), .])
    ids = selectindex(Gm.estimable)
    att = J(rows(tdate), 1, .)
    diffm = J(0, q, .)
    if (rows(ids)) {
        diffm = Pj.beta[nodes[ids, 1], .] - Pj.beta[nodes[ids, 2], .]
        att[ids] = rowsum(diffm :* tx[tgroup[ids], .])
    }
    status = J(rows(tdate), 1, 0)
    for (i = 1; i <= rows(tdate); i++) status[i] = (Gm.estimable[i] ? 2 : C.observed[tdate[i], tgroup[i]])
    out.targets = (cal[C.groups[tgroup]], cal[tdate], tdate :< C.groups[tgroup], att, status)
    out.IF = J(n, rows(tdate), .)
    if (keepnuis) _cdid_ewr_obs_capture(ob, E)
    if (keepnuis) _cdid_ewr_proj_capture(C, Gm, Pj, "corrected_levels", rows(Gm.use), E)
    if (!compute_if | !rows(ids)) return(out)
    Cm = J(Pj.p, rows(ids), 0)
    for (a = 1; a <= q; a++) {
        for (k = 1; k <= rows(ids); k++) {
            Cm[nodes[ids[k], 1] + (a - 1) * Pj.A, k] = tx[tgroup[ids[k]], a]
            Cm[nodes[ids[k], 2] + (a - 1) * Pj.A, k] = Cm[nodes[ids[k], 2] + (a - 1) * Pj.A, k] - tx[tgroup[ids[k]], a]
        }
    }
    V = Pj.score * Cm
    for (gi = 2; gi <= nG; gi++) {
        v = selectindex(tgroup[ids] :== gi)
        if (!rows(v)) continue
        u = selectindex(coh :== C.groups[gi])
        V[u, v'] = V[u, v'] + (n / rows(u)) :* ((W[u, .] :- tx[gi, .]) * diffm[v, .]')
    }
    out.IF[., ids'] = V
    return(out)
}

/* ------------------------------------------------------------------------------------------
   Corrected treatment IPW and DR: links of each treated cohort g (its own pairs, and the pairs
   of the donor cohorts untreated at t), with responses normalized on the full cohort frame.
   ------------------------------------------------------------------------------------------ */
struct cdid_etwfe_fit scalar cdid_etwfe_corrected_tm(real matrix Y, real matrix O, real colvector coh,
    real colvector cal, real scalar never, real scalar dr, real matrix Zc, struct _cdid_ewc_set scalar S,
    string scalar weighting, real scalar compute_if, real scalar trigger, real scalar ridge)
{
    struct _cdid_ewc_cat scalar C
    struct _cdid_ewc_geo scalar Gc
    struct _cdid_ewo scalar ob
    struct _cdid_ewc_proj scalar Pj
    struct _cdid_ew_tm scalar T
    struct _cdid_ewr_env scalar E
    struct cdid_xp_design scalar D
    struct _cdid_ew_odds colvector fits
    struct cdid_logit_result scalar lf
    real matrix lay, W, X, mu, delta, req, cellsfit, Gr, lv, gam, imb, Hm, IFl, slot, models, grid
    real colvector treated, ok, keep, comp, lcomp, rcomp, rows, val, tilt, act, e, fitc, center, amp, rk, donor, tgt, cl
    real colvector mhat, den, rmean, mix, gx, v
    real rowvector load
    real scalar K, n, nG, L, j, gi, g, c, s, t, i, a, m, nm, q, scored, keepnuis, nyt, p1, p2, k, gg
    K = cols(O)
    n = rows(O)
    E = _cdid_ewr_get()
    keepnuis = E.on
    nyt = !never
    scored = (compute_if | weighting == "optimal")
    C = _cdid_ewc_catalog(O, coh)
    nG = rows(C.groups)
    treated = select(C.groups, C.groups :> 0)
    /* Layout: (g, cell, donor) sorted by g, s, t, target first, donor. */
    lay = J(0, 5, .)
    for (i = 1; i <= rows(treated); i++) {
        g = treated[i]
        for (c = 1; c <= rows(C.cells); c++) {
            gg = C.groups[C.cells[c, 1]]
            if (gg == g) lay = lay \ (g, c, ., C.cells[c, 2], C.cells[c, 3])
            else if (gg == 0 | (nyt & C.cells[c, 3] < gg)) lay = lay \ (g, c, gg, C.cells[c, 2], C.cells[c, 3])
        }
    }
    if (rows(lay)) lay = lay[order((lay[., (1, 4, 5)], lay[., 3] :< ., editmissing(lay[., 3], -1)), (1, 2, 3, 4, 5)), .]
    /* Keep the links in incidence components of observed, requested contrasts. */
    if (rows(lay)) {
        D = _cdid_ew_incidence(lay[., (1, 4, 5, 3)], treated, K, never)
        comp = _cdid_ew_inc_comp(D)
        keep = J(0, 1, .)
        for (i = 1; i <= rows(treated); i++) {
            gi = selectindex(C.groups :== treated[i])
            for (t = 1; t <= K; t++) {
                if (!C.observed[t, gi]) continue
                if (!(t >= treated[i] | (never & t < treated[i] - 1))) continue
                if (t < treated[i]) {
                    p1 = _cdid_ew_vertex(D, i, t, 2)
                    p2 = _cdid_ew_vertex(D, i, t, 4)
                }
                else {
                    p1 = _cdid_ew_vertex(D, i, t, 7)
                    p2 = _cdid_ew_vertex(D, i, t, 6)
                }
                if (p1 == 0 | p2 == 0) continue
                if (comp[p1] == comp[p2]) keep = keep \ comp[p1]
            }
        }
        lcomp = comp[D.from]
        ok = J(rows(lay), 1, 0)
        for (i = 1; i <= rows(lay); i++) ok[i] = any(keep :== lcomp[i])
        lay = select(lay, ok)
    }
    L = rows(lay)
    tgt = (L ? (lay[., 3] :>= .) : J(0, 1, .))
    /* Calendar projection on the clean pairs (DR with donor links). */
    W = _cdid_ewc_features(Zc, center, amp)
    q = cols(W)
    Gc.use = J(0, 1, .)
    if (dr & L) if (any(!tgt)) {
        slot = J(K, nG, .)
        for (gi = 1; gi <= nG; gi++) slot[., gi] = (1::K)
        rk = J(0, 1, .)
        for (c = 1; c <= rows(C.cells); c++) {
            gg = C.groups[C.cells[c, 1]]
            if (gg == 0 | (nyt & C.cells[c, 3] < gg)) rk = rk \ c
        }
        req = select(lay[., (4, 5)], !tgt)
        Gc = _cdid_ewc_geometry(C, slot, rk, req)
    }
    fitc = (L ? uniqrows(lay[., 2] \ Gc.use) : J(0, 1, .))
    if (rows(fitc)) {
        cellsfit = (C.groups[C.cells[fitc, 1]], C.cells[fitc, 2..3])
        ob = _cdid_ewo_fit(Y, O, coh, cellsfit, S.selmethod, S.selmodel, S.selpool, S.imppool, S.Zs, S.sel_int,
            S.U, S.rho, S.universe, S.Xs, S.imp_int, S.eps, S.maxit, (scored | keepnuis))
    }
    else {
        ob.P = 0
        ob.IF = J(n, 0, .)
    }
    if (keepnuis & rows(fitc)) _cdid_ewr_obs_capture(ob, E)
    T.links = (L ? lay[., (1, 4, 5, 3)] : J(0, 4, .))
    T.mhat = J(L, 1, .)
    T.Phi = J(0, 0, .)
    if (!L) {
        if (keepnuis) _cdid_ewr_tm_empty(E)
        return(_cdid_ew_tm_solve(T, O, coh, never, weighting, compute_if, trigger, ridge))
    }
    X = J(n, 1, 1), Zc
    mu = J(nG, q, .)
    for (gi = 1; gi <= nG; gi++) mu[gi, .] = mean(W[selectindex(coh :== C.groups[gi]), .])
    /* Treatment logits, in order of first use. */
    models = J(0, 2, .)
    for (j = 1; j <= L; j++) {
        if (tgt[j]) continue
        if (!rows(models)) models = lay[j, (1, 3)]
        else if (!any((models[., 1] :== lay[j, 1]) :& (models[., 2] :== lay[j, 3]))) models = models \ lay[j, (1, 3)]
    }
    nm = rows(models)
    fits = _cdid_ew_odds(max((nm, 1)), 1)
    for (m = 1; m <= nm; m++) {
        fits[m].g = models[m, 1]
        fits[m].k = models[m, 2]
        fits[m].cg = selectindex(C.groups :== models[m, 1])
        fits[m].ids = selectindex((coh :== models[m, 1]) :| (coh :== models[m, 2]))
        fits[m].cols = selectindex((lay[., 1] :== models[m, 1]) :& (lay[., 3] :== models[m, 2]))
        lf = cdid_logit_fit(X[fits[m].ids, .], (coh[fits[m].ids] :== models[m, 1]), J(rows(fits[m].ids), 1, 1), S.eps, S.maxit, 1)
        fits[m].beta = lf.beta
        lv = _cdid_ipw_logit_values(X[fits[m].ids, .] * lf.beta)
        Hm = cross(X[fits[m].ids, .], lv[., 1] :* (1 :- lv[., 1]), X[fits[m].ids, .]) :/ n
        fits[m].Hc = cholesky(Hm)
        if (missing(fits[m].Hc)) _cdid_ew_tm_fail()
    }
    delta = J(L, q, 0)
    Pj.p = 0
    if (dr & any(!tgt)) {
        Pj = _cdid_ewc_project(C, W, coh, Gc, ob, (scored | keepnuis))
        for (j = 1; j <= L; j++) if (!tgt[j]) delta[j, .] = Pj.beta[lay[j, 5], .] - Pj.beta[lay[j, 4], .]
    }
    den = rmean = J(L, 1, .)
    gam = imb = J(L, q, 0)
    if (scored) gam = J(L, cols(X), 0)
    load = J(1, ob.P, 0)
    IFl = (scored ? J(n, L, 0) : J(n, 0, .))
    for (j = 1; j <= L; j++) {
        c = _cdid_ewc_fitpos(ob, C, lay[j, 2])
        _cdid_ewo_read(ob, c, rows, val, Gr, scored)
        if (tgt[j]) {
            tilt = J(rows(rows), 1, 1)
            act = J(rows(rows), 1, 0)
            e = val
        }
        else {
            m = selectindex((models[., 1] :== lay[j, 1]) :& (models[., 2] :== lay[j, 3]))
            lv = _cdid_ipw_logit_values(X[rows, .] * fits[m].beta)
            tilt = lv[., 2]
            act = lv[., 3]
            e = (dr ? val - W[rows, .] * delta[j, .]' : val)
        }
        den[j] = sum(tilt)
        rmean[j] = sum(tilt :* e) / den[j]
        gi = selectindex(C.groups :== lay[j, 1])
        T.mhat[j] = rmean[j] + mu[gi, .] * delta[j, .]'
        if (!scored) continue
        /* Scores of link j. */
        if (ob.P) IFl[., j] = ob.IF * (cross(Gr, tilt) :/ den[j])
        if (tgt[j]) IFl[rows, j] = IFl[rows, j] + n :* (val :- rmean[j]) :/ den[j]
        else {
            gx = cross(X[rows, .], tilt :* act :* (e :- rmean[j])) :/ den[j]
            IFl[rows, j] = IFl[rows, j] + n :* tilt :* (e :- rmean[j]) :/ den[j]
            v = fits[m].ids
            lv = _cdid_ipw_logit_values(X[v, .] * fits[m].beta)
            Hm = X[v, .] :* ((coh[v] :== fits[m].g) - lv[., 1])
            IFl[v, j] = IFl[v, j] + solveupper(fits[m].Hc', solvelower(fits[m].Hc, Hm'))' * gx
            if (dr) {
                rk = selectindex(coh :== lay[j, 1])
                IFl[rk, j] = IFl[rk, j] + (n / rows(rk)) :* ((W[rk, .] :- mu[gi, .]) * delta[j, .]')
                mix = mu[gi, .]' - cross(W[rows, .], tilt) :/ den[j]
                for (a = 1; a <= q; a++) IFl[., j] = IFl[., j] + (Pj.score[., lay[j, 5] + (a - 1) * K] - Pj.score[., lay[j, 4] + (a - 1) * K]) :* mix[a]
            }
        }
    }
    if (missing(T.mhat) | any(den :<= 0)) _cdid_ew_tm_fail()
    if (scored) {
        if (missing(IFl)) _cdid_ew_tm_fail()
        T.Phi = IFl
    }
    if (keepnuis) {
        _cdid_ewr_tm_capture(fits, nm, X, coh, E)
        if (dr & any(!tgt)) _cdid_ewr_proj_capture(C, Gc, Pj, "corrected_calendar", sum(!tgt), E)
    }
    return(_cdid_ew_tm_solve(T, O, coh, never, weighting, compute_if, trigger, ridge))
}

/* ------------------------------------------------------------------------------------------
   Reports of the fitted nuisance models (keepnuisance). Standard errors use the N-scaled
   original-unit influence functions: sums within clusters (or units), centered, divided by
   N^2, missing for a block with fewer than two contributing clusters.
   ------------------------------------------------------------------------------------------ */
real rowvector _cdid_ewr_var(real matrix IF, real colvector owner, struct _cdid_ewr_env scalar E)
{
    real matrix Sg, hit
    real rowvector s1, s2, v
    real colvector bl, cols
    real scalar n, G, i, b
    n = rows(IF)
    if (!cols(IF)) return(J(1, 0, .))
    if (!E.compute_if) return(J(1, cols(IF), .))
    if (E.clustered) {
        G = E.C
        Sg = J(G, cols(IF), 0)
        for (i = 1; i <= n; i++) Sg[E.cl[i], .] = Sg[E.cl[i], .] + IF[i, .]
    }
    else {
        G = n
        Sg = IF
    }
    if (G < 2) return(J(1, cols(IF), .))
    s1 = colsum(Sg)
    s2 = colsum(Sg :* Sg)
    v = (s2 - s1 :* s1 :/ G) :/ (n * n)
    v = v :* (v :> 0)
    if (E.clustered) {
        bl = uniqrows(owner)
        for (i = 1; i <= rows(bl); i++) {
            cols = selectindex(owner :== bl[i])
            hit = J(G, 1, 0)
            for (b = 1; b <= n; b++) if (any(IF[b, cols'] :!= 0)) hit[E.cl[b]] = 1
            if (sum(hit) < 2) v[cols'] = J(1, rows(cols), .)
        }
    }
    return(v)
}

real scalar _cdid_ewr_nclus(real matrix IF, struct _cdid_ewr_env scalar E)
{
    real colvector hit
    real scalar b
    if (!E.clustered) return(.)
    hit = J(E.C, 1, 0)
    for (b = 1; b <= rows(IF); b++) if (any(IF[b, .] :!= 0)) hit[E.cl[b]] = 1
    return(sum(hit))
}

string scalar _cdid_ewr_num(real scalar x)
{
    return(strofreal(x, "%21.0g"))
}

real scalar _cdid_ewr_calv(struct _cdid_ewr_env scalar E, real scalar k)
{
    return(k == 0 ? 0 : E.cal[k])
}

void _cdid_ewr_zp(real matrix coef)
{
    if (!rows(coef)) return
    coef[., 5] = coef[., 3] :/ coef[., 4]
    coef[., 6] = 2 :* normal(-abs(coef[., 5]))
}

struct cdid_nuis_cap scalar _cdid_ewr_pending()
{
    pointer() scalar p
    struct cdid_nuis_cap scalar cap
    struct _cdid_ewr_env scalar E
    p = findexternal("__cdid_nuis_pending")
    if (p != NULL) {
        cap = *p
        if (rows(cap.comp) == 4 & cols(cap.roles) == 3) return(cap)
    }
    E = _cdid_ewr_get()
    cap.mode = "chained"
    cap.comp = ("selection", "", "", "") \ ("outcome", "", "", "") \ ("treatment", "", "", "") \ ("imputation", "", "", "")
    cap.roles = cdid_nuis_role(1, 3)
    cap.roles[1].role = "outcome"
    cap.roles[2].role = "treatment"
    cap.roles[3].role = "imputation"
    cap.N = E.N
    cap.clustered = E.clustered
    cap.n_groups = (E.clustered ? E.C : E.N)
    cap.has_intercept = 1
    cap.floor = 1e-6
    cap.maxit = 50
    cap.comparisons = J(0, 5, "")
    cap.n_catalog = cap.n_aux = 0
    return(cap)
}

/* Outcome regression on linked slots. beta: all coefficients; keep: fitted coordinates;
   cidx: slot by feature coordinates; sd, sg: calendar date and cohort (0 never) of each slot;
   gauge: fixed slots; IFk: N-scaled influence functions of the fitted coordinates. */
void _cdid_ewr_outcome(real colvector beta, real colvector keep, real matrix cidx, real colvector sd,
    real colvector sg, real colvector gauge, real matrix IFk, real scalar nfit, real scalar ncomp,
    string scalar model, string scalar sample, real colvector center, real colvector amp,
    struct _cdid_ewr_env scalar E)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar R
    real colvector all, kept, pos
    real rowvector v
    real scalar i, a, s, q, k, row
    string rowvector feat
    string scalar eq, nm, coh_txt
    cap = _cdid_ewr_pending()
    if (!rows(keep)) {
        cap.comp[2, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted outcome regression in this call.")
        _cdid_nuis_store_pending(cap)
        return
    }
    q = cols(cidx)
    feat = "(Intercept)"
    if (q > 1) feat = feat, E.xnames[1..q - 1]
    kept = J(0, 1, .)
    for (s = 1; s <= rows(cidx); s++) if (any(keep :== cidx[s, 1])) kept = kept \ s
    all = uniqrows(selectindex(gauge) \ kept)
    v = _cdid_ewr_var(IFk, J(rows(keep), 1, 1), E)
    R.role = "outcome"
    R.title = "Calendar outcome regression, one linked slot system"
    R.response = (model == "levels" ? "observed outcome levels Y_it" : (model == "calendar_changes" ?
        "observed untreated pair changes" : "observation-corrected pair response Z_i(k;s,t)"))
    R.sample = sample
    R.weights = "unit weights"
    R.pooling = "one linked slot system"
    R.coef = J(rows(all) * q, 6, .)
    R.coef_txt = J(rows(all) * q, 2, "")
    row = 0
    coh_txt = ""
    for (i = 1; i <= rows(all); i++) {
        s = all[i]
        for (a = 1; a <= q; a++) {
            row++
            k = cidx[s, a]
            pos = selectindex(keep :== k)
            R.coef[row, .] = (1, i, beta[k], (rows(pos) ? sqrt(v[pos[1]]) : .), ., .)
            nm = (sg[s] == 0 ? "date[" + _cdid_ewr_num(sd[s]) + "]" : "cohort_date[" + _cdid_ewr_num(sg[s]) + "," + _cdid_ewr_num(sd[s]) + "]")
            R.coef_txt[row, .] = (nm + ":" + feat[a], (gauge[s] ? "fixed_gauge" : (a == 1 ? "level" : "scaled_slope")))
        }
        if (!strpos("," + coh_txt + ",", "," + _cdid_ewr_num(sg[s]) + ","))
            coh_txt = coh_txt + (coh_txt == "" ? "" : ",") + _cdid_ewr_num(sg[s])
    }
    _cdid_ewr_zp(R.coef)
    eq = (model == "levels" ? "Y_it = alpha_i + W_i' b[slot(G_i,t)]" : (model == "calendar_changes" ?
        "Y_it - Y_is = W_i' (b[t] - b[s]) + e_i(s,t)" : (model == "corrected_levels" ?
        "Z_i(k;s,t) = W_i' (b[slot(k,t)] - b[slot(k,s)]) + e_i(k;s,t)" : "Z_i(k;s,t) = W_i' (b[t] - b[s]) + e_i(k;s,t)")))
    eq = eq + (q > 1 ? " with W_i = (1, (Z_i - center) / amplitude)" : " with W_i = 1") +
        ". The lowest slot of each fitted connected component is fixed at zero."
    if (model == "levels") eq = eq + " Unit effects are profiled out."
    R.blocks = (rows(keep), 0, rows(keep), nfit, nfit, _cdid_ewr_nclus(IFk, E), ncomp, .)
    R.block_txt = ("calendar", "Anchored calendar coefficients", eq, coh_txt, (q > 1 ? invtokens(E.xnames[1..q - 1], ",") : ""), "")
    cap.roles[1] = R
    cap.comp[2, 2..4] = ("available", "", "")
    _cdid_nuis_store_pending(cap)
}

void _cdid_ewr_proj_capture(struct _cdid_ewc_cat scalar C, struct _cdid_ewc_geo scalar Gm,
    struct _cdid_ewc_proj scalar Pj, string scalar model, real scalar ncomp, struct _cdid_ewr_env scalar E)
{
    real matrix cidx
    real colvector sd, sg, gauge, ends, center, amp
    real scalar s, a, K, nfit, i
    K = C.K
    if (!rows(Pj.keep)) {
        _cdid_ewr_outcome(J(0, 1, .), J(0, 1, .), J(0, 1, .), J(0, 1, .), J(0, 1, .), J(0, 1, .), J(rows(Pj.score), 0, .), 0, 0, model, "", J(0, 1, .), J(0, 1, .), E)
        return
    }
    cidx = J(Pj.A, Pj.q, .)
    for (a = 1; a <= Pj.q; a++) cidx[., a] = (1::Pj.A) :+ (a - 1) * Pj.A
    sd = sg = J(Pj.A, 1, .)
    for (s = 1; s <= Pj.A; s++) {
        sd[s] = E.cal[mod(s - 1, K) + 1]
        sg[s] = (model == "corrected_calendar" ? 0 : _cdid_ewr_calv(E, C.groups[floor((s - 1) / K) + 1]))
    }
    ends = uniqrows(Gm.left \ Gm.right)
    gauge = J(Pj.A, 1, 0)
    for (i = 1; i <= rows(ends); i++) if (!any(Gm.keep_slots :== ends[i])) gauge[ends[i]] = 1
    nfit = 0
    for (i = 1; i <= rows(Gm.use); i++) nfit = nfit + C.nunits[C.cells[Gm.use[i], 1]]
    _cdid_ewr_outcome(vec(Pj.beta), Pj.keep, cidx, sd, sg, gauge, Pj.score[., Pj.keep'], nfit, ncomp, model,
        "full cohort-frame units on supported pair equations", J(0, 1, .), J(0, 1, .), E)
}

/* Treatment propensity logits of each target and donor cohort. */
void _cdid_ewr_tm_capture(struct _cdid_ew_odds colvector fits, real scalar nm, real matrix X,
    real colvector coh, struct _cdid_ewr_env scalar E)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar R
    real matrix IFm, lv, u
    real rowvector v
    real scalar m, j, q, row, n, nk
    string rowvector terms
    string scalar sl
    cap = _cdid_ewr_pending()
    if (!nm) {
        cap.comp[3, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted treatment propensity in this call.")
        _cdid_nuis_store_pending(cap)
        return
    }
    n = rows(X)
    q = cols(X)
    terms = "(Intercept)"
    if (q > 1) terms = terms, E.xnames[1..q - 1]
    sl = ""
    for (j = 2; j <= q; j++) sl = sl + " + b[" + terms[j] + "] * " + terms[j]
    R.role = "treatment"
    R.title = "Two-cohort treatment propensity logits"
    R.response = "D_i = 1{G_i = g}"
    R.sample = "all frame units in each target-donor risk set"
    R.weights = "unit weights"
    R.pooling = "two-cohort risk set"
    R.blocks = J(nm, 8, .)
    R.block_txt = J(nm, 6, "")
    R.coef = J(nm * q, 6, .)
    R.coef_txt = J(nm * q, 2, "")
    row = 0
    for (m = 1; m <= nm; m++) {
        IFm = J(n, q, 0)
        lv = _cdid_ipw_logit_values(X[fits[m].ids, .] * fits[m].beta)
        u = X[fits[m].ids, .] :* ((coh[fits[m].ids] :== fits[m].g) - lv[., 1])
        IFm[fits[m].ids, .] = solveupper(fits[m].Hc', solvelower(fits[m].Hc, u'))'
        v = _cdid_ewr_var(IFm, J(q, 1, 1), E)
        for (j = 1; j <= q; j++) {
            row++
            R.coef[row, .] = (m, ., fits[m].beta[j], sqrt(v[j]), ., .)
            R.coef_txt[row, .] = (terms[j], (j == 1 ? "intercept" : "slope"))
        }
        nk = rows(fits[m].ids)
        R.blocks[m, .] = (q, 1, q - 1, nk, nk, _cdid_ewr_nclus(IFm, E), rows(fits[m].cols), .)
        R.block_txt[m, .] = ("treatment_" + strofreal(m), "target " + _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].g)) + " vs donor " +
            _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].k)), "logit Pr(G_i = " + _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].g)) +
            " | Z_i, G_i in {" + _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].g)) + "," + _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].k)) +
            "}) = a" + sl, _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].g)) + "," + _cdid_ewr_num(_cdid_ewr_calv(E, fits[m].k)),
            (q > 1 ? invtokens(terms[2..q], ",") : ""), "")
    }
    _cdid_ewr_zp(R.coef)
    cap.roles[2] = R
    cap.comp[3, 2..4] = ("available", "", "")
    _cdid_nuis_store_pending(cap)
}

void _cdid_ewr_tm_empty(struct _cdid_ewr_env scalar E)
{
    struct cdid_nuis_cap scalar cap
    cap = _cdid_ewr_pending()
    cap.comp[3, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted treatment propensity in this call.")
    _cdid_nuis_store_pending(cap)
}

/* Selection logits, roster probabilities and actual-change imputations of a corrected fit. */
void _cdid_ewr_obs_capture(struct _cdid_ewo scalar ob, struct _cdid_ewr_env scalar E)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar R
    real rowvector v, kc, ac
    real colvector cells, idx
    real scalar b, k, c, nc, row, nb, j, g, nrow, ncl
    string colvector key
    string scalar sl, lab, eqn, gtxt
    cap = _cdid_ewr_pending()
    v = _cdid_ewr_var(ob.IF, ob.owner, E)
    if (ob.need_q) {
        cap.sel_model = E.selmodel
        cap.sel_pooling = E.selpool
        cap.sel_universe = E.universe
        cap.has_intercept = ob.has_int
        cap.slope_names = (ob.ds ? E.snames[1..ob.ds] : J(1, 0, ""))
        key = J(ob.nsc, 1, "")
        cap.cells = J(ob.nsc, 12, .)
        for (c = 1; c <= ob.nsc; c++) {
            if (ob.sc_s[c]) key[c] = "g=" + _cdid_ewr_num(_cdid_ewr_calv(E, ob.sc_k[c])) + ",(" + _cdid_ewr_num(E.cal[ob.sc_s[c]]) + "," + _cdid_ewr_num(E.cal[ob.sc_t[c]]) + ")"
            else key[c] = "g=" + _cdid_ewr_num(_cdid_ewr_calv(E, ob.sc_k[c])) + ",t=" + _cdid_ewr_num(E.cal[ob.sc_t[c]])
            cap.cells[c, .] = (c, ob.sc_block[c], (ob.sc_s[c] ? 1 : 2), _cdid_ewr_calv(E, ob.sc_k[c]),
                (ob.sc_s[c] ? E.cal[ob.sc_s[c]] : .), E.cal[ob.sc_t[c]], ob.sc_nU[c], round(ob.sc_nU[c] * ob.sc_pi[c]), ob.sc_pi[c], ., ., .)
        }
        cap.cell_key = key
        nb = ob.nsb + (ob.roster ? 1 : 0)
        cap.blocks = J(nb, 14, .)
        cap.block_txt = J(nb, 3, "")
        cap.coef = J(0, 6, .)
        cap.coef_txt = J(0, 2, "")
        sl = ""
        for (j = 1; j <= ob.ds; j++) sl = sl + " + b" + strofreal(j) + "*" + cap.slope_names[j]
        for (b = 1; b <= ob.nsb; b++) {
            cells = selectindex(ob.sc_block :== b)
            nc = rows(cells)
            if (ob.has_int) for (k = 1; k <= nc; k++) {
                c = cells[k]
                cap.coef = cap.coef \ (b, c, ob.theta[ob.sc_alpha[c]], sqrt(v[ob.sc_alpha[c]]), ., .)
                cap.coef_txt = cap.coef_txt \ ("a[" + key[c] + "]", "local_intercept")
            }
            for (j = 1; j <= ob.ds; j++) {
                cap.coef = cap.coef \ (b, ., ob.theta[ob.sb_beta[b, j]], sqrt(v[ob.sb_beta[b, j]]), ., .)
                cap.coef_txt = cap.coef_txt \ (cap.slope_names[j], (nc > 1 ? "shared_slope" : "slope"))
            }
            idx = (ob.has_int ? ob.sc_alpha[cells] : J(0, 1, .)) \ (ob.ds ? ob.sb_beta[b, .]' : J(0, 1, .))
            ncl = (rows(idx) ? _cdid_ewr_nclus(ob.IF[., idx'], E) : .)
            nrow = sum(ob.sc_nU[cells])
            gtxt = _cdid_ewr_glist(E, uniqrows(ob.sc_k[cells]))
            cap.blocks[b, .] = (nc, (ob.has_int ? nc : 0) + ob.ds, (ob.has_int ? nc : 0), ob.ds, nrow,
                _cdid_ewr_units(ob, cells), ncl, nrow, 1, ob.sb_iter[b], ob.sb_rank[b], ., _cdid_ewr_ncomp(ob, b),
                (rows(uniqrows(ob.sc_k[cells])) == 1 ? _cdid_ewr_calv(E, ob.sc_k[cells[1]]) : .))
            if (E.selpool == "separate") lab = (ob.sc_s[cells[1]] ? "pair[" : "date[") + key[cells[1]] + "]"
            else if (E.selpool == "within_cohort") lab = "within_cohort[g=" + _cdid_ewr_num(_cdid_ewr_calv(E, ob.sc_k[cells[1]])) + "]"
            else lab = "across_cohorts[g in {" + gtxt + "}]"
            eqn = (ob.sc_s[cells[1]] ? "logit Pr(S_i(s,t) = 1 | X_i, G_i, U_i = 1) = " : "logit Pr(R_it = 1 | X_i, G_i) = ")
            eqn = eqn + (ob.has_int ? (nc > 1 ? "a[cell]" : "a") : (ob.ds ? "" : "0")) + (ob.has_int ? sl : substr(sl, 4, .))
            cap.block_txt[b, .] = ("selection_" + strofreal(b), lab, eqn)
        }
        if (ob.roster) {
            b = ob.nsb + 1
            for (k = 1; k <= rows(ob.coh_vals); k++) {
                cap.coef = cap.coef \ (b, ., ob.pi_val[k], sqrt(v[ob.pi_index[k]]), ., .)
                cap.coef_txt = cap.coef_txt \ ("Pr(U | G=" + _cdid_ewr_num(_cdid_ewr_calv(E, ob.coh_vals[k])) + ")", "roster_probability")
            }
            cap.blocks[b, .] = (0, rows(ob.coh_vals), 0, 0, sum(ob.n_k), sum(ob.n_k), _cdid_ewr_nclus(ob.IF[., ob.pi_index'], E), sum(ob.n_k), ., ., ., ., ., .)
            cap.block_txt[b, .] = ("roster", "Pr(U | cohort)", "Pr(U_i = 1 | G_i = g), estimated by the share of cohort units in U")
        }
        _cdid_ewr_zp(cap.coef)
        cap.comp[1, 2..4] = ("available", "", "")
    }
    if (ob.need_j) {
        R.role = "imputation"
        R.title = "Actual-change imputation regressions"
        R.response = "observed outcome change Y_it - Y_is on a supported pair"
        R.sample = "observed pairs within each cohort"
        R.weights = "unit weights"
        R.pooling = E.imppool
        R.blocks = J(ob.nig, 8, .)
        R.block_txt = J(ob.nig, 6, "")
        R.coef = J(0, 6, .)
        R.coef_txt = J(0, 2, "")
        key = J(ob.C, 1, "")
        for (c = 1; c <= ob.C; c++) key[c] = "g=" + _cdid_ewr_num(_cdid_ewr_calv(E, ob.cell_k[c])) + "," + _cdid_ewr_num(E.cal[ob.cell_s[c]]) + ":" + _cdid_ewr_num(E.cal[ob.cell_t[c]])
        for (g = 1; g <= ob.nig; g++) {
            kc = selectindex(ob.ig_keep[g, .])
            ac = selectindex(!ob.ig_keep[g, .])
            sl = ""
            if (cols(kc)) {
                for (j = 1; j <= cols(kc); j++) sl = sl + (j > 1 ? " + " : "") + "b[" + E.xnames[kc[j]] + "] * " + E.xnames[kc[j]]
                sl = "(t-s) * (" + sl + ")"
            }
            cells = select(ob.ig_order, ob.cell_group[ob.ig_order] :== g)
            nc = rows(cells)
            lab = ""
            if (ob.icpt) for (k = 1; k <= nc; k++) {
                c = cells[k]
                R.coef = R.coef \ (g, c, ob.theta[ob.cell_alpha[c]], sqrt(v[ob.cell_alpha[c]]), ., .)
                R.coef_txt = R.coef_txt \ ("a[" + key[c] + "]", "local_intercept")
            }
            for (k = 1; k <= nc; k++) lab = lab + (k > 1 ? ", " : "") + key[cells[k]]
            for (j = 1; j <= cols(kc); j++) {
                R.coef = R.coef \ (g, ., ob.theta[ob.ig_beta[g, kc[j]]], sqrt(v[ob.ig_beta[g, kc[j]]]), ., .)
                R.coef_txt = R.coef_txt \ (E.xnames[kc[j]], "slope")
            }
            idx = (ob.icpt ? ob.cell_alpha[cells] : J(0, 1, .)) \ (cols(kc) ? ob.ig_beta[g, kc]' : J(0, 1, .))
            R.blocks[g, .] = (rows(idx), (ob.icpt ? nc : 0), cols(kc), ob.ig_nfit[g], ob.ig_nfit[g], _cdid_ewr_nclus(ob.IF[., idx'], E), nc, .)
            R.block_txt[g, .] = ("imputation_" + strofreal(g), lab, "E[Y_it - Y_is | G_i = k, X_i, S_i(s,t) = 1] = " +
                (ob.icpt ? "a[k,s,t]" + (cols(kc) ? " + " : "") : "") + sl, _cdid_ewr_glist(E, uniqrows(ob.cell_k[cells])),
                (cols(kc) ? invtokens(E.xnames[kc], ",") : ""), (cols(ac) ? invtokens(E.xnames[ac], ",") : ""))
        }
        _cdid_ewr_zp(R.coef)
        cap.roles[3] = R
        cap.comp[4, 2..4] = ("available", "", "")
    }
    _cdid_nuis_store_pending(cap)
}

string scalar _cdid_ewr_glist(struct _cdid_ewr_env scalar E, real colvector k)
{
    string scalar s
    real scalar i
    s = ""
    for (i = 1; i <= rows(k); i++) s = s + (i > 1 ? "," : "") + _cdid_ewr_num(_cdid_ewr_calv(E, k[i]))
    return(s)
}

real scalar _cdid_ewr_units(struct _cdid_ewo scalar ob, real colvector cells)
{
    real colvector k
    real scalar i, s
    k = uniqrows(ob.sc_k[cells])
    s = 0
    for (i = 1; i <= rows(k); i++) s = s + ob.sc_nU[selectindex(ob.sc_k :== k[i])[1]]
    return(s)
}

/* Number of fitted pairs whose observation probability uses block b. */
real scalar _cdid_ewr_ncomp(struct _cdid_ewo scalar ob, real scalar b)
{
    real scalar i, n
    n = 0
    for (i = 1; i <= ob.C; i++) {
        if (ob.sc_block[ob.cell_sc1[i]] == b) n++
        else if (ob.cell_sc2[i]) if (ob.sc_block[ob.cell_sc2[i]] == b) n++
    }
    return(n)
}

/* Final statuses of the four components of a cohort_time_fe capture. */
void _cdid_ewr_finish(string scalar method, string scalar selection, string scalar selmethod)
{
    struct cdid_nuis_cap scalar cap
    cap = _cdid_ewr_pending()
    if (cap.comp[1, 2] == "") {
        if (selection == "none") cap.comp[1, 2..4] = ("not_applicable", "no_selection_model", "selection(none) fits no observation model.")
        else if (selmethod == "imputation") cap.comp[1, 2..4] = ("not_applicable", "pure_imputation", "selmethod(imputation) fits no observation logit.")
        else cap.comp[1, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted observation logit in this call.")
    }
    if (cap.comp[2, 2] == "") {
        if (method == "ipw") cap.comp[2, 2..4] = ("not_applicable", "no_outcome_regression", "method(ipw) fits no outcome regression.")
        else cap.comp[2, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted outcome regression in this call.")
    }
    if (cap.comp[3, 2] == "") {
        if (method == "reg") cap.comp[3, 2..4] = ("not_applicable", "no_treatment_propensity", "method(reg) fits no treatment propensity.")
        else cap.comp[3, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted treatment propensity in this call.")
    }
    if (cap.comp[4, 2] == "") {
        if (selection == "none" | selmethod == "ipw") cap.comp[4, 2..4] = ("not_applicable", "no_imputation", "This observation-correction choice fits no actual-change imputation model.")
        else cap.comp[4, 2..4] = ("unavailable", "no_supported_model", "No supported equation required a fitted actual-change imputation model in this call.")
    }
    _cdid_nuis_store_pending(cap)
}

/* Within least squares of the uncorrected regression (levels). The report solves every
   component without covariates and the components with a reported contrast otherwise. */
void _cdid_ewr_levels(real matrix Y, real matrix O, real colvector group, real matrix sbd,
    real colvector seen, real colvector tc, real matrix feat, real matrix gram, real matrix rhs,
    real colvector comp, real colvector kept, real matrix cidx, real colvector nobs)
{
    struct _cdid_ewr_env scalar E
    real colvector raw, sd, sg, first, keepS, gauge, keep, used, d, s, y, fit, e, bvec, grp
    real matrix U, IFk
    real scalar nslot, nF, K, k, i, a, n, nfit
    E = _cdid_ewr_get()
    n = rows(Y)
    K = cols(Y)
    nslot = rows(cidx)
    nF = cols(feat)
    raw = selectindex(seen)
    sd = sg = J(nslot, 1, .)
    for (k = 1; k <= nslot; k++) {
        sd[k] = E.cal[mod(raw[k] - 1, K) + 1]
        i = floor((raw[k] - 1) / K) + 1
        sg[k] = (i == 1 ? 0 : E.cal[tc[i - 1]])
    }
    first = J(nslot, 1, 0)
    for (k = 1; k <= nslot; k++) first[k] = (k == 1 ? 1 : !any(comp[1..k - 1] :== comp[k]))
    if (nF == 1) {
        keepS = !first
        gauge = first
    }
    else {
        keepS = kept
        used = uniqrows(select(comp, kept))
        gauge = J(nslot, 1, 0)
        for (k = 1; k <= nslot; k++) gauge[k] = (first[k] & any(used :== comp[k]))
    }
    keep = (any(keepS) ? vec(cidx[selectindex(keepS), .]) : J(0, 1, .))
    if (!rows(keep)) {
        _cdid_ewr_outcome(J(0, 1, .), keep, cidx, sd, sg, gauge, J(n, 0, .), 0, 0, "levels", "", J(0, 1, .), J(0, 1, .), E)
        return
    }
    bvec = J(rows(cidx) * nF, 1, 0)
    bvec[keep] = _cdid_ewo_solve(gram[keep, keep], vec(rhs)[keep])
    if (missing(bvec)) _cdid_ew_fail()
    U = J(n, rows(bvec), 0)
    nfit = 0
    for (i = 1; i <= n; i++) {
        if (nobs[i] < 2) continue
        d = selectindex(O[i, .]')
        s = sbd[d, group[i]]
        if (nF == 1) nfit = nfit + rows(d)
        else if (gauge[s[1]] | keepS[s[1]]) nfit = nfit + rows(d)
        y = Y[i, d']'
        fit = J(rows(d), 1, 0)
        for (a = 1; a <= nF; a++) fit = fit + bvec[cidx[s, a]] :* feat[i, a]
        e = (y :- mean(y)) - (fit :- mean(fit))
        for (a = 1; a <= nF; a++) U[i, cidx[s, a]'] = U[i, cidx[s, a]'] + (e :* feat[i, a])'
    }
    IFk = n :* (_cdid_ewo_solve(gram[keep, keep], U[., keep']'))'
    _cdid_ewr_outcome(bvec, keep, cidx, sd, sg, gauge, IFk, nfit, 0, "levels",
        "observed rows on repeat-observed units", J(0, 1, .), J(0, 1, .), E)
}

/* Common untreated calendar regression of the uncorrected DR fit. */
void _cdid_ewr_calendar(real matrix Y, real matrix O, real colvector coh, real scalar never,
    real matrix W, real matrix B, real colvector keepc, real colvector comp, real colvector used_d,
    real scalar ndonor)
{
    struct _cdid_ewr_env scalar E
    real matrix Cm, e, uu, sc, IFk, pc, cidx, hitm
    real colvector mcol, gauge, sd, sg, usedc, hit
    real scalar N, K, q, rr, d, a, k, d1, d2, nfit
    E = _cdid_ewr_get()
    N = rows(Y)
    K = cols(Y)
    q = cols(W)
    cidx = J(K, q, .)
    for (d = 1; d <= K; d++) for (a = 1; a <= q; a++) cidx[d, a] = (d - 1) * q + a
    usedc = uniqrows(comp[used_d])
    gauge = J(K, 1, 0)
    for (d = 1; d <= K; d++) gauge[d] = ((d == 1 ? 1 : !any(comp[1..d - 1] :== comp[d])) & any(usedc :== comp[d]))
    sd = E.cal
    sg = J(K, 1, 0)
    Cm = O :* ((coh :== 0) :| (!never :& ((J(N, 1, 1) * (1..K)) :< coh)))
    mcol = rowsum(Cm)
    hit = selectindex(mcol :> 1)
    pc = cross(Cm[hit, .], Cm[hit, .])
    nfit = 0
    for (d1 = 1; d1 < K; d1++) for (d2 = d1 + 1; d2 <= K; d2++) if (any(usedc :== comp[d1]) & any(usedc :== comp[d2])) nfit = nfit + pc[d1, d2]
    e = Cm :* (Y - W * B')
    uu = Cm :* (mcol :* e :- rowsum(e))
    sc = J(N, rows(keepc), 0)
    for (rr = 1; rr <= rows(keepc); rr++) {
        d = ceil(keepc[rr] / q)
        a = keepc[rr] - (d - 1) * q
        sc[., rr] = uu[., d] :* W[., a]
    }
    IFk = N :* (_cdid_ewo_solve(_cdid_ewr_calgram(Cm, W, keepc, q), sc'))'
    _cdid_ewr_outcome(vec(B'), keepc, cidx, sd, sg, gauge, IFk, nfit, ndonor, "calendar_changes",
        "clean observed pairs", J(0, 1, .), J(0, 1, .), E)
}

/* Gram of the clean-pair calendar regression on the fitted coordinates. */
real matrix _cdid_ewr_calgram(real matrix Cm0, real matrix W, real colvector keepc, real scalar q)
{
    real matrix Cm, Kx, gram, u
    real colvector hit, mcol
    real scalar K, d, a
    K = cols(Cm0)
    hit = selectindex(rowsum(Cm0) :> 1)
    Cm = Cm0[hit, .]
    u = W[hit, .]
    mcol = rowsum(Cm)
    Kx = J(rows(Cm), K * q, 0)
    for (d = 1; d <= K; d++) for (a = 1; a <= q; a++) Kx[., (d - 1) * q + a] = Cm[., d] :* u[., a]
    gram = -cross(Kx, Kx)
    for (d = 1; d <= K; d++) gram[((d - 1) * q + 1)..(d * q), ((d - 1) * q + 1)..(d * q)] =
        gram[((d - 1) * q + 1)..(d * q), ((d - 1) * q + 1)..(d * q)] + cross(u, mcol :* Cm[., d], u)
    return(gram[keepc, keepc])
}

end
