*! Linear cohort-by-date fixed-effects regression (preset(cohort_time_fe)).
version 19.5
mata:

/* Fit of y_it = alpha_i + mu_s(i,t) + x_i'b_s(i,t) on observed rows, where the slot s is
   the common date t (never-treated rows and rows outside the free cells) or the free
   cohort-date cell (g, t): with never-treated controls every date except g - 1, with
   not-yet-treated controls dates t >= g. Unit effects are eliminated within units; one
   gauge slot per connected component of the slots linked by observed rows is fixed at
   zero. Contrasts are evaluated at the cohort mean of the covariates. targets: cohort,
   period, pre flag, estimate, status (0 no support, 1 not identified, 2 identified);
   IF: N-scaled original-unit influence functions (missing for unidentified targets). */
struct cdid_etwfe_fit {
    real matrix targets, IF
}

/* Settings shared by the report builders of one fit: unit clusters, names and calendar. */
struct _cdid_ewr_env {
    real colvector cl, cal
    real scalar C, clustered, on, N, compute_if
    string rowvector xnames, snames
    string scalar selmodel, selpool, imppool, universe, method, selmethod
}

/* Observation-model settings of one corrected fit. */
struct _cdid_ewc_set {
    string scalar selmethod, selmodel, selpool, imppool, universe
    real matrix Zs, Xs
    real scalar sel_int, imp_int, eps, maxit
    real colvector U, rho
}

struct _cdid_ew_pat {
    real colvector ids
    real rowvector m
}

/* Units of one group partitioned by observed-date pattern. */
struct _cdid_ew_pat colvector _cdid_ew_patterns(real matrix O, real colvector ids)
{
    struct _cdid_ew_pat colvector P
    real matrix M
    real colvector ord
    real scalar i, n, k
    n = rows(ids)
    if (n == 0) return(_cdid_ew_pat(0, 1))
    M = O[ids, .]
    ord = order((M, ids), (1..(cols(M) + 1)))
    P = _cdid_ew_pat(n, 1)
    k = 0
    for (i = 1; i <= n; i++) {
        if (i == 1 | any(M[ord[i], .] :!= M[ord[max((i - 1, 1))], .])) {
            k++
            P[k].m = M[ord[i], .]
            P[k].ids = J(0, 1, .)
        }
        P[k].ids = P[k].ids \ ids[ord[i]]
    }
    return(P[1..k])
}

real colvector _cdid_ew_components(real matrix link)
{
    real colvector comp, front, reach
    real scalar p, s, nc
    p = rows(link)
    comp = J(p, 1, 0)
    nc = 0
    for (s = 1; s <= p; s++) {
        if (comp[s]) continue
        nc++
        comp[s] = nc
        front = s
        while (rows(front)) {
            reach = (colsum(link[front, .]) :> 0)'
            front = selectindex(reach :& (comp :== 0))
            if (rows(front)) comp[front] = J(rows(front), 1, nc)
        }
    }
    return(comp)
}

struct cdid_etwfe_fit scalar cdid_etwfe_core(real matrix Y, real matrix O, real colvector q,
    real colvector cal, real colvector tc, real scalar never, real matrix Zc, real scalar compute_if)
{
    struct cdid_etwfe_fit scalar out
    struct _cdid_ew_pat colvector Pt
    real scalar n, K, nG, nF, nfull, nslot, ncoef, j, t, gg, a, b, k, u, nid, rc
    real matrix freem, raw, sbd, feat, cmean, gram, rhs, link, cidx, beta, sc, mix, G, L, C, tr, W, Yw, SW, blk, fit, res, score, cen
    real colvector group, seen, lookup, tdate, tgroup, cell, dslot, nobs, rep, ids, s, used, comp, ident, kept, keep, scale, att, center, amp, roster, cols, status, ss, first
    real rowvector m
    n = rows(Y)
    K = cols(Y)
    nG = rows(tc)
    group = q :+ 1
    freem = J(K, nG, 0)
    for (j = 1; j <= nG; j++) freem[., j] = (never ? ((1::K) :!= tc[j] - 1) : ((1::K) :>= tc[j]))
    raw = J(K, nG + 1, .)
    raw[., 1] = (1::K)
    for (j = 1; j <= nG; j++) raw[., j + 1] = freem[., j] :* ((1::K) :+ j * K) + (!freem[., j]) :* (1::K)
    nfull = K * (nG + 1)
    seen = J(nfull, 1, 0)
    for (gg = 1; gg <= nG + 1; gg++) {
        ids = selectindex(group :== gg)
        if (rows(ids) == 0) continue
        m = (colsum(O[ids, .]) :> 0)
        for (t = 1; t <= K; t++) if (m[t]) seen[raw[t, gg]] = 1
    }
    nslot = sum(seen)
    lookup = J(nfull, 1, 0)
    lookup[selectindex(seen)] = (1::nslot)
    sbd = J(K, nG + 1, 0)
    for (gg = 1; gg <= nG + 1; gg++) sbd[., gg] = lookup[raw[., gg]]
    tdate = J(0, 1, .)
    tgroup = J(0, 1, .)
    for (j = 1; j <= nG; j++) {
        for (t = 1; t <= K; t++) if (freem[t, j]) {
            tdate = tdate \ t
            tgroup = tgroup \ (j + 1)
        }
    }
    cell = J(rows(tdate), 1, 0)
    dslot = J(rows(tdate), 1, 0)
    for (k = 1; k <= rows(tdate); k++) {
        cell[k] = sbd[tdate[k], tgroup[k]]
        dslot[k] = sbd[tdate[k], 1]
    }
    /* Features: intercept and centered, scaled covariates. */
    if (cols(Zc)) {
        center = mean(Zc)'
        cen = Zc :- center'
        amp = colmax(abs(cen))'
        amp = amp + (amp :== 0)
        feat = J(n, 1, 1), cen :/ amp'
    }
    else feat = J(n, 1, 1)
    nF = cols(feat)
    cmean = J(nG + 1, nF, .)
    for (gg = 1; gg <= nG + 1; gg++) {
        roster = selectindex(group :== gg)
        if (rows(roster)) cmean[gg, .] = mean(feat[roster, .])
    }
    nobs = rowsum(O)
    rep = (nobs :> 1)
    ncoef = nslot * nF
    cidx = J(nslot, nF, .)
    for (a = 1; a <= nF; a++) cidx[., a] = (1::nslot) :+ (a - 1) * nslot
    gram = J(ncoef, ncoef, 0)
    rhs = J(nslot, nF, 0)
    link = J(nslot, nslot, 0)
    for (gg = 1; gg <= nG + 1; gg++) {
        ids = selectindex((group :== gg) :& rep)
        if (rows(ids) == 0) continue
        Pt = _cdid_ew_patterns(O, ids)
        for (k = 1; k <= rows(Pt); k++) {
            used = selectindex(Pt[k].m')
            s = sbd[used, gg]
            W = feat[Pt[k].ids, .]
            Yw = Y[Pt[k].ids, used']
            Yw = Yw :- rowsum(Yw) :/ rows(used)
            rhs[s, .] = rhs[s, .] + cross(Yw, W)
            link[s, s] = J(rows(s), rows(s), 1)
            SW = cross(W, W)
            blk = I(rows(used)) - J(rows(used), rows(used), 1 / rows(used))
            for (a = 1; a <= nF; a++) {
                for (b = 1; b <= nF; b++) gram[cidx[s, a], cidx[s, b]] = gram[cidx[s, a], cidx[s, b]] + SW[a, b] :* blk
            }
        }
    }
    comp = _cdid_ew_components(link)
    ident = (cell :> 0) :& (dslot :> 0)
    for (k = 1; k <= rows(ident); k++) if (ident[k]) ident[k] = (comp[cell[k]] == comp[dslot[k]])
    first = J(nslot, 1, 0)
    for (k = 1; k <= nslot; k++) first[k] = (k == 1 ? 1 : !any(comp[1..k - 1] :== comp[k]))
    kept = J(nslot, 1, 0)
    for (k = 1; k <= nslot; k++) {
        if (first[k]) continue
        for (j = 1; j <= rows(ident); j++) if (ident[j]) if (comp[cell[j]] == comp[k]) {
            kept[k] = 1
            break
        }
    }
    keep = vec(cidx[selectindex(kept), .])
    status = (cell :> 0) + ident
    beta = J(ncoef, 1, 0)
    if (rows(keep)) {
        G = gram[keep, keep]
        if (missing(G) | any(diagonal(G) :<= 0)) _cdid_ew_fail()
        scale = sqrt(diagonal(G))
        G = G :/ (scale * scale')
        G = makesymmetric(G)
        L = cholesky(G)
        if (missing(L)) _cdid_ew_fail()
        rc = 1 / cond(G)
        if (rc <= epsilon(1)) _cdid_ew_fail()
        beta[keep] = (cholsolve(G, vec(rhs)[keep] :/ scale)) :/ scale
    }
    beta = colshape(beta, nslot)'
    if (_cdid_ewr_on()) _cdid_ewr_levels(Y, O, group, sbd, seen, tc, feat, gram, rhs, comp, kept, cidx, nobs)
    nid = sum(ident)
    ids = selectindex(ident)
    sc = J(nid, nF, .)
    mix = cmean[tgroup, .]
    att = J(rows(tdate), 1, .)
    if (nid) {
        sc = beta[cell[ids], .] - beta[dslot[ids], .]
        att[ids] = rowsum(sc :* mix[ids, .])
    }
    out.targets = (cal[tc[tgroup :- 1]], cal[tdate], tdate :< tc[tgroup :- 1], att, status)
    out.IF = J(n, rows(tdate), .)
    if (!compute_if) return(out)
    out.IF[., ids] = J(n, nid, 0)
    if (!nid) return(out)
    C = J(ncoef, nid, 0)
    for (a = 1; a <= nF; a++) {
        for (k = 1; k <= nid; k++) {
            C[cidx[cell[ids[k]], a], k] = mix[ids[k], a]
            C[cidx[dslot[ids[k]], a], k] = C[cidx[dslot[ids[k]], a], k] - mix[ids[k], a]
        }
    }
    tr = J(ncoef, nid, 0)
    if (rows(keep)) tr[keep, .] = cholsolve(G, C[keep, .] :/ scale) :/ scale
    for (gg = 1; gg <= nG + 1; gg++) {
        ids = selectindex((group :== gg) :& rep)
        if (rows(ids) == 0) continue
        Pt = _cdid_ew_patterns(O, ids)
        for (k = 1; k <= rows(Pt); k++) {
            used = selectindex(Pt[k].m')
            s = sbd[used, gg]
            W = feat[Pt[k].ids, .]
            Yw = Y[Pt[k].ids, used']
            Yw = Yw :- rowsum(Yw) :/ rows(used)
            fit = W * beta[s, .]'
            res = Yw - (fit :- rowsum(fit) :/ rows(used))
            score = J(rows(Pt[k].ids), nid, 0)
            for (a = 1; a <= nF; a++) score = score + (res * tr[cidx[s, a], .]) :* W[., a]
            out.IF[Pt[k].ids, selectindex(ident)'] = n :* score
        }
    }
    if (nF > 1) {
        ids = selectindex(ident)
        for (gg = 2; gg <= nG + 1; gg++) {
            cols = selectindex(tgroup[ids] :== gg)
            if (rows(cols) == 0) continue
            roster = selectindex(group :== gg)
            out.IF[roster, ids[cols]'] = out.IF[roster, ids[cols]'] +
                (n / rows(roster)) :* ((feat[roster, .] :- cmean[gg, .]) * sc[cols, .]')
        }
    }
    return(out)
}

void _cdid_ew_fail()
{
    errprintf("the cohort-by-date fixed-effects system is numerically singular after fixing one level per connected set of periods; no estimates are returned\n")
    exit(459)
}


/* Treatment IPW and DR on cohort-standardized calendar paths (R etwfe_treatment_moments and
   etwfe_treatment_core). Links: for each treated cohort g, its own changes (s, t) and the
   changes of each donor cohort k untreated at t (never-treated, or with not-yet-treated
   controls any other cohort). A donor change is reweighted toward cohort g by the odds of
   the logit of g against k on (1, x); with DR, a common untreated calendar regression of
   changes on (1, x) fitted on the clean pairs is evaluated at the mean of cohort g and
   the odds-weighted residual change is added. Scores are N-scaled original-unit IFs. */
struct _cdid_ew_tm {
    real matrix links, Phi
    real colvector mhat
}

struct _cdid_ew_odds {
    real colvector ids, beta, cols
    real matrix Hc
    real scalar g, k, cg
}

struct _cdid_ew_tm scalar _cdid_ew_tm_fit(real matrix Y, real matrix O, real colvector coh,
    real matrix Z, real scalar never, real scalar dr, real scalar scores, real scalar eps, real scalar maxit)
{
    struct _cdid_ew_tm scalar T
    struct _cdid_ew_odds colvector fits
    struct cdid_logit_result scalar lf
    real scalar N, K, q, nl, c, g, k, s, t, L, m, a, nm, ci, ck, rc, d, nb, rr
    real colvector level, nroster, nobs, ids, mem, clean_ids, used_d, comp, keepc, scale, beta, model, gg, ss, tt, dd, tgt, cols, hit, act, mass, rmean, center, amp, cc, rgs, rks
    real matrix W, X, Wbar, B, G, Lc, gram, rhs, link, Cm, Kx, yc, E, Ow, M0, Den, M1, M2, M3, OX, lv, gamma, imb, P, u, phi, e, uu, H, sc, F, PN, PDY, pairs
    real colvector mcol
    N = rows(Y)
    K = cols(Y)
    level = uniqrows(coh)
    nl = rows(level)
    nroster = J(nl, 1, 0)
    for (c = 1; c <= nl; c++) nroster[c] = sum(coh :== level[c])
    nobs = rowsum(O)
    q = cols(Z) + 1
    center = (cols(Z) ? mean(Z)' : J(0, 1, .))
    W = J(N, q, 1)
    amp = J(q - 1, 1, 1)
    for (a = 1; a <= q - 1; a++) {
        cc = Z[., a] :- center[a]
        if (max(abs(cc)) > 0) amp[a] = max(abs(cc))
        W[., a + 1] = cc :/ amp[a]
    }
    X = J(N, 1, 1), Z
    PN = PDY = J(K, K * nl, 0)
    for (c = 1; c <= nl; c++) {
        mem = selectindex((coh :== level[c]) :& (nobs :> 1))
        if (rows(mem)) {
            PN[., ((c - 1) * K + 1)..(c * K)] = cross(O[mem, .], O[mem, .])
            P = cross(O[mem, .], Y[mem, .])
            PDY[., ((c - 1) * K + 1)..(c * K)] = P - P'
        }
    }
    /* Links in (g, s, t, target first, donor) order. */
    gg = ss = tt = dd = J(0, 1, .)
    for (ci = 1; ci <= nl; ci++) {
        g = level[ci]
        if (g == 0) continue
        for (s = 1; s < K; s++) for (t = s + 1; t <= K; t++) {
            if (PN[s, (ci - 1) * K + t] > 0) {
                gg = gg \ g
                ss = ss \ s
                tt = tt \ t
                dd = dd \ .
            }
            for (ck = 1; ck <= nl; ck++) {
                k = level[ck]
                if (k == g) continue
                if (never & k != 0) continue
                if (PN[s, (ck - 1) * K + t] > 0 & (k == 0 | t < k)) {
                    gg = gg \ g
                    ss = ss \ s
                    tt = tt \ t
                    dd = dd \ k
                }
            }
        }
    }
    L = rows(gg)
    T.links = (gg, ss, tt, dd)
    T.mhat = J(L, 1, .)
    mass = J(L, 1, .)
    rmean = J(L, 1, 0)
    tgt = (dd :>= .)
    for (m = 1; m <= L; m++) if (tgt[m]) {
        ci = selectindex(level :== gg[m])
        mass[m] = PN[ss[m], (ci - 1) * K + tt[m]]
        T.mhat[m] = PDY[ss[m], (ci - 1) * K + tt[m]] / mass[m]
    }
    /* Odds models in role order. */
    model = J(L, 1, 0)
    pairs = J(0, 3, .)
    for (ci = 1; ci <= nl; ci++) {
        g = level[ci]
        if (g == 0) continue
        for (ck = 1; ck <= nl; ck++) {
            k = level[ck]
            if (k == g | (never & k != 0)) continue
            if (any((gg :== g) :& (dd :== k))) pairs = pairs \ (g, k, ci)
        }
    }
    nm = rows(pairs)
    fits = _cdid_ew_odds(max((nm, 1)), 1)
    for (m = 1; m <= nm; m++) {
        hit = selectindex((gg :== pairs[m, 1]) :& (dd :== pairs[m, 2]))
        fits[m].g = pairs[m, 1]
        fits[m].k = pairs[m, 2]
        fits[m].cg = pairs[m, 3]
        fits[m].cols = hit
        model[hit] = J(rows(hit), 1, m)
    }
    /* Common untreated calendar regression (DR). */
    B = J(K, q, 0)
    if (dr & nm) {
        clean_ids = selectindex((nobs :> 1) :& ((coh :== 0) :| (!never :& (coh :> 2))))
        gram = J(K * q, K * q, 0)
        rhs = J(1, K * q, 0)
        link = J(K, K, 0)
        if (rows(clean_ids)) {
            Cm = O[clean_ids, .] :* ((coh[clean_ids] :== 0) :| (!never :& ((J(rows(clean_ids), 1, 1) * (1..K)) :< coh[clean_ids])))
            hit = selectindex(rowsum(Cm) :> 1)
            Cm = Cm[hit, .]
            u = W[clean_ids[hit], .]
            yc = Y[clean_ids[hit], .] :* Cm
            mcol = rowsum(Cm)
            link = (cross(Cm, Cm) :> 0)
            Kx = J(rows(Cm), K * q, 0)
            for (d = 1; d <= K; d++) for (a = 1; a <= q; a++) Kx[., (d - 1) * q + a] = Cm[., d] :* u[., a]
            gram = -cross(Kx, Kx)
            for (d = 1; d <= K; d++) gram[((d - 1) * q + 1)..(d * q), ((d - 1) * q + 1)..(d * q)] =
                gram[((d - 1) * q + 1)..(d * q), ((d - 1) * q + 1)..(d * q)] + cross(u, mcol :* Cm[., d], u)
            e = mcol :* yc :- rowsum(yc)
            for (d = 1; d <= K; d++) for (a = 1; a <= q; a++) rhs[(d - 1) * q + a] = sum(Kx[., (d - 1) * q + a] :* e[., d])
        }
        comp = _cdid_ew_components(link)
        used_d = select(ss \ tt, (dd \ dd) :< .)
        keepc = J(K * q, 1, 0)
        for (d = 1; d <= K; d++) {
            if (d > 1) if (!any(comp[1..d - 1] :== comp[d])) continue
            if (d == 1) continue
            if (!any(comp[used_d] :== comp[d])) continue
            keepc[((d - 1) * q + 1)..(d * q)] = J(q, 1, 1)
        }
        keepc = selectindex(keepc)
        if (!rows(keepc)) _cdid_ew_tm_fail()
        G = gram[keepc, keepc]
        if (missing(G) | any(diagonal(G) :<= 0)) _cdid_ew_tm_fail()
        scale = sqrt(diagonal(G))
        G = makesymmetric(G :/ (scale * scale'))
        Lc = cholesky(G)
        if (missing(Lc)) _cdid_ew_tm_fail()
        if (1 / cond(G) <= epsilon(1)) _cdid_ew_tm_fail()
        beta = J(K * q, 1, 0)
        beta[keepc] = cholsolve(G, rhs[keepc]' :/ scale) :/ scale
        B = rowshape(beta, K)
    }
    /* Odds and donor moments. */
    Wbar = J(nl, q, .)
    for (c = 1; c <= nl; c++) Wbar[c, .] = mean(W[selectindex(coh :== level[c]), .])
    if (scores) {
        gamma = J(q, L, 0)
        imb = J(q, L, 0)
    }
    for (m = 1; m <= nm; m++) {
        ci = fits[m].cg
        ck = selectindex(level :== fits[m].k)
        ids = selectindex(coh :== fits[m].g) \ selectindex(coh :== fits[m].k)
        fits[m].ids = ids
        lf = cdid_logit_fit(X[ids, .], (coh[ids] :== fits[m].g), J(rows(ids), 1, 1), eps, maxit, 1)
        fits[m].beta = lf.beta
        lv = _cdid_ipw_logit_values(X[ids, .] * lf.beta)
        H = cross(X[ids, .], lv[., 1] :* (1 :- lv[., 1]), X[ids, .]) :/ N
        fits[m].Hc = cholesky(H)
        if (missing(fits[m].Hc)) _cdid_ew_tm_fail()
        mem = selectindex((coh :== fits[m].k) :& (nobs :> 1))
        cols = fits[m].cols
        M0 = Den = J(K, K, 0)
        if (scores) M1 = M2 = M3 = J(K * q, K, 0)
        if (rows(mem)) {
            lv = _cdid_ipw_logit_values(X[mem, .] * lf.beta)
            E = (dr ? O[mem, .] :* (Y[mem, .] - W[mem, .] * B') : Y[mem, .])
            Ow = O[mem, .] :* lv[., 2]
            M0 = cross(Ow, E)
            Den = cross(Ow, O[mem, .])
            if (scores) {
                OX = J(rows(mem), K * q, 0)
                for (a = 1; a <= q; a++) OX[., ((a - 1) * K + 1)..(a * K)] = Ow :* (X[mem, a] :* lv[., 3])
                M1 = cross(OX, E)
                M2 = cross(OX, O[mem, .])
                if (dr) {
                    for (a = 1; a <= q; a++) OX[., ((a - 1) * K + 1)..(a * K)] = Ow :* W[mem, a]
                    M3 = cross(OX, O[mem, .])
                }
            }
        }
        for (rr = 1; rr <= rows(cols); rr++) {
            c = cols[rr]
            mass[c] = Den[ss[c], tt[c]]
            rmean[c] = (M0[ss[c], tt[c]] - M0[tt[c], ss[c]]) / mass[c]
            T.mhat[c] = (B[tt[c], .] - B[ss[c], .]) * Wbar[ci, .]' + rmean[c]
            if (scores) {
                for (a = 1; a <= q; a++) {
                    gamma[a, c] = (M1[(a - 1) * K + ss[c], tt[c]] - M1[(a - 1) * K + tt[c], ss[c]] - rmean[c] * M2[(a - 1) * K + ss[c], tt[c]]) / mass[c]
                    if (dr) imb[a, c] = Wbar[ci, a] - M3[(a - 1) * K + ss[c], tt[c]] / mass[c]
                }
            }
        }
    }
    if (missing(T.mhat)) _cdid_ew_tm_fail()
    if (_cdid_ewr_on()) {
        _cdid_ewr_tm_capture(fits, nm, X, coh, _cdid_ewr_get())
        if (dr & nm) _cdid_ewr_calendar(Y, O, coh, never, W, B, keepc, comp, used_d, sum(!tgt))
    }
    T.Phi = J(0, 0, .)
    if (!scores) return(T)
    /* Scores. */
    T.Phi = J(N, L, 0)
    for (ci = 1; ci <= nl; ci++) {
        if (level[ci] == 0) continue
        cols = selectindex(tgt :& (gg :== level[ci]))
        ids = selectindex(coh :== level[ci])
        if (!rows(cols) | !rows(ids)) continue
        for (rr = 1; rr <= rows(cols); rr++) {
            c = cols[rr]
            T.Phi[ids, c] = N :* O[ids, ss[c]] :* O[ids, tt[c]] :* (Y[ids, tt[c]] - Y[ids, ss[c]] :- T.mhat[c]) :/ mass[c]
        }
    }
    P = W * B'
    for (m = 1; m <= nm; m++) {
        ids = fits[m].ids
        cols = fits[m].cols
        lv = _cdid_ipw_logit_values(X[ids, .] * fits[m].beta)
        u = X[ids, .] :* ((coh[ids] :== fits[m].g) - lv[., 1])
        sc = solveupper(fits[m].Hc', solvelower(fits[m].Hc, u'))'
        T.Phi[ids, cols'] = T.Phi[ids, cols'] + sc * gamma[., cols']
        rks = selectindex(coh[ids] :== fits[m].k)
        if (rows(rks)) {
            E = (dr ? O[ids[rks], .] :* (Y[ids[rks], .] - P[ids[rks], .]) : Y[ids[rks], .])
            for (rr = 1; rr <= rows(cols); rr++) {
                c = cols[rr]
                T.Phi[ids[rks], c] = T.Phi[ids[rks], c] + N :* lv[rks, 2] :* O[ids[rks], ss[c]] :* O[ids[rks], tt[c]] :*
                    (E[., tt[c]] - E[., ss[c]] :- rmean[c]) :/ mass[c]
            }
        }
        rgs = selectindex(coh[ids] :== fits[m].g)
        if (dr & rows(rgs)) {
            T.Phi[ids[rgs], cols'] = T.Phi[ids[rgs], cols'] + (N / nroster[fits[m].cg]) :*
                ((W[ids[rgs], .] :- Wbar[fits[m].cg, .]) * (B[tt[cols], .] - B[ss[cols], .])')
        }
    }
    if (dr & nm) {
        if (rows(keepc)) {
            Cm = O :* ((coh :== 0) :| (!never :& ((J(N, 1, 1) * (1..K)) :< coh)))
            mcol = rowsum(Cm)
            e = Cm :* (Y - P)
            uu = Cm :* (mcol :* e :- rowsum(e))
            sc = J(N, rows(keepc), 0)
            for (rr = 1; rr <= rows(keepc); rr++) {
                d = ceil(keepc[rr] / q)
                a = keepc[rr] - (d - 1) * q
                sc[., rr] = uu[., d] :* W[., a]
            }
            phi = J(N, K * q, 0)
            phi[., keepc'] = N :* (cholsolve(G, sc' :/ scale) :/ scale)'
            cols = selectindex(!tgt)
            for (rr = 1; rr <= rows(cols); rr++) {
                c = cols[rr]
                for (a = 1; a <= q; a++) T.Phi[., c] = T.Phi[., c] + (phi[., (tt[c] - 1) * q + a] - phi[., (ss[c] - 1) * q + a]) :* imb[a, c]
            }
        }
    }
    if (missing(T.Phi)) _cdid_ew_tm_fail()
    return(T)
}

void _cdid_ew_tm_fail()
{
    errprintf("the cohort-by-date treatment moments could not be computed: a logit or the common calendar regression is singular or not finite\n")
    exit(459)
}

/* GMM fit of the treatment moments on the cross-pair incidence with shared calendar levels
   for the donors. */
struct cdid_etwfe_fit scalar cdid_etwfe_treatment(real matrix Y, real matrix O, real colvector coh,
    real matrix Z, real scalar never, real scalar dr, string scalar weighting, real scalar compute_if,
    real scalar eps, real scalar maxit, real scalar trigger, real scalar ridge)
{
    struct _cdid_ew_tm scalar T
    real scalar scored
    scored = (compute_if | weighting == "optimal")
    T = _cdid_ew_tm_fit(Y, O, coh, Z, never, dr, scored, eps, maxit)
    return(_cdid_ew_tm_solve(T, O, coh, never, weighting, compute_if, trigger, ridge))
}

/* Incidence of the treatment links (g, s, t, donor; donor missing for target links): control
   links share the calendar vertices (key 0). */
struct cdid_xp_design scalar _cdid_ew_incidence(real matrix links, real colvector treated,
    real scalar K, real scalar never)
{
    struct cdid_xp_cat scalar C
    real scalar L, i
    L = rows(links)
    C.cal = (1::K)
    C.tg = treated
    C.tref = treated :- 1
    C.g = links[., 1]
    C.cs = links[., 2]
    C.ct = links[., 3]
    C.role = 1 :+ (links[., 4] :< .)
    C.ci = J(L, 1, 0)
    C.gi = J(L, 1, 0)
    for (i = 1; i <= L; i++) C.gi[i] = selectindex(treated :== C.g[i])
    C.ref = C.g :- 1
    return(_cdid_xp_design(C, (1::L), 0, (never ? "post" : "all"), (never ? . : 1)))
}

/* Connected components (root vertex) of the incidence. */
real colvector _cdid_ew_inc_comp(struct cdid_xp_design scalar D)
{
    real colvector comp
    real scalar i, ru, rv
    comp = (1::rows(D.VX))
    for (i = 1; i <= rows(D.from); i++) {
        ru = D.from[i]
        while (comp[ru] != ru) ru = comp[ru]
        rv = D.to[i]
        while (comp[rv] != rv) rv = comp[rv]
        if (ru != rv) comp[rv] = ru
    }
    for (i = 1; i <= rows(comp); i++) {
        ru = i
        while (comp[ru] != ru) ru = comp[ru]
        comp[i] = ru
    }
    return(comp)
}

/* GMM fit of the treatment moments T on the cross-pair incidence. */
struct cdid_etwfe_fit scalar _cdid_ew_tm_solve(struct _cdid_ew_tm scalar T, real matrix O,
    real colvector coh, real scalar never, string scalar weighting, real scalar compute_if,
    real scalar trigger, real scalar ridge)
{
    struct cdid_etwfe_fit scalar out
    struct cdid_xp_design scalar D
    struct cdid_gmm_unit_gram scalar src
    struct cdid_gmm_result scalar sol
    real scalar N, K, L, j, t, nt, ru, rv, i, k, ok
    real colvector treated, tq, tp, status, comp, v, cf, hrow, est, obsc
    real matrix H, V
    real rowvector h
    N = rows(O)
    K = cols(O)
    treated = uniqrows(select(coh, coh :> 0))
    tq = tp = J(0, 1, .)
    for (j = 1; j <= rows(treated); j++) for (t = 1; t <= K; t++) {
        if (never ? t != treated[j] - 1 : t >= treated[j]) {
            tq = tq \ j
            tp = tp \ t
        }
    }
    nt = rows(tq)
    status = J(nt, 1, 0)
    for (i = 1; i <= nt; i++) status[i] = (sum(O[selectindex(coh :== treated[tq[i]]), tp[i]]) > 0)
    out.targets = (treated[tq], tp, tp :< treated[tq], J(nt, 1, .), status)
    out.IF = J(N, nt, .)
    L = rows(T.mhat)
    if (!L) return(out)
    D = _cdid_ew_incidence(T.links, treated, K, never)
    comp = _cdid_ew_inc_comp(D)
    H = J(0, D.P, 0)
    est = J(0, 1, .)
    for (i = 1; i <= nt; i++) {
        if (!status[i]) continue
        if (tp[i] < treated[tq[i]]) {
            v = _cdid_ew_vertex(D, tq[i], tp[i], 2) \ _cdid_ew_vertex(D, tq[i], tp[i], 4)
        }
        else v = _cdid_ew_vertex(D, tq[i], tp[i], 7) \ _cdid_ew_vertex(D, tq[i], tp[i], 6)
        if (any(v :== 0)) continue
        cf = (1 \ -1)
        ok = 1
        for (k = 1; k <= 2; k++) {
            if (D.anch[v[k]]) continue
            if (sum(cf :* (comp[v] :== comp[v[k]])) != 0) ok = 0
        }
        if (!ok) continue
        h = J(1, D.P, 0)
        for (k = 1; k <= 2; k++) if (D.vpar[v[k]]) h[D.vpar[v[k]]] = h[D.vpar[v[k]]] + cf[k]
        H = H \ h
        est = est \ i
    }
    if (!rows(est) | !D.P) return(out)
    if (weighting == "identity") sol = cdid_gmm_identity_map(T.mhat, D.A, H, 2^53 - 1)
    else {
        V = cross(T.Phi, T.Phi) :/ N
        V = makesymmetric(V)
        src.n_units = N
        src.rows_consumed = N
        src.normalization_count = N
        src.normalization_factor = 1 / N
        src.sample_rank_upper = N - 1
        src.mean = colsum(T.Phi) :/ N
        src.retained_if = J(0, 0, .)
        src.covariance = V
        sol = _cdid_gmm_fit_core(T.mhat, D.A, V, I(L), L, H, "optimal", trigger, ridge, 2^53 - 1, src)
    }
    out.targets[est, 4] = sol.targets
    out.targets[est, 5] = J(rows(est), 1, 2)
    if (compute_if) out.IF[., est'] = T.Phi * sol.target_loadings'
    return(out)
}

real scalar _cdid_ew_vertex(struct cdid_xp_design scalar D, real scalar j, real scalar u, real scalar role)
{
    real colvector h
    h = selectindex((D.VX[., 1] :== j) :& (D.VX[., 2] :== u) :& (D.VX[., 3] :== role))
    return(rows(h) == 1 ? h[1] : 0)
}

/* Public fit of preset(cohort_time_fe), regression branch. Reported cells: with
   never-treated controls every cohort-date cell except the reference date g - 1 (pre cells
   through their direct descriptors, post cells through the post ledger); with
   not-yet-treated controls the post cells only. */
struct cdid_public_fit_result scalar cdid_etwfe_public_fit(struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_public_report_buffer scalar buf
    struct cdid_public_bundle scalar B2
    struct cdid_etwfe_fit scalar F
    struct _cdid_ewr_env scalar E
    struct _cdid_ewc_set scalar S
    real matrix Y, O, grid, Zc, postIF, auxIF
    real colvector cal, coh, q, tc, keep, post_point, post_av, aux_point, aux_sup, rr, k, post_ids, aux_ids, hit, cellof, tgt, r, cohc
    string colvector post_status, post_reason, aux_reason
    string rowvector xn
    real scalar N, K, t, j, never, h, qd, sg, c, p, corrected
    corrected = (B.prepared.selection_canonical != "none" & B.prepared.selection_canonical != "")
    if (B.prepared.design == "reference_recipes")
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires covariates that are constant within units and no weights.")
    if (B.prepared.method != "reg" & B.prepared.method != "ipw" & B.prepared.method != "dr")
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires method(reg), method(ipw), or method(dr).")
    N = B.prepared.population.N
    cal = B.prepared.population.calendar
    K = rows(cal)
    if (K < 2) _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires consecutive positive integer study dates.")
    if (any(cal :<= 0) | any(cal :!= floor(cal)) | any(cal[2..K] - cal[1..K - 1] :!= 1))
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires consecutive positive integer study dates.")
    never = (B.request.control_group == "nevertreated")
    if (!never & B.request.pt_start != cal[1])
        _cdid_pr_fail("cdid_chained_request_input_error", "With not-yet-treated controls, preset(cohort_time_fe) requires ptassumption(all) and ptstart() at the first study date.")
    coh = B.prepared.population.cohort
    if (any(B.prepared.population.cohort_known :== 0))
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires known cohorts.")
    tc = uniqrows(select(coh, coh :< .))
    if (!rows(tc) | !any(coh :>= .))
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires a never-treated cohort and treated cohorts with adoption dates after the first study date.")
    for (j = 1; j <= rows(tc); j++) {
        hit = selectindex(cal[2..K] :== tc[j])
        if (rows(hit) != 1) _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires a never-treated cohort and treated cohorts with adoption dates after the first study date.")
        tc[j] = hit[1] + 1
    }
    q = J(N, 1, 0)
    for (j = 1; j <= rows(tc); j++) q = q + j :* (coh :== cal[tc[j]])
    if (rows(B.prepared.population.obs_grid) == N & cols(B.prepared.population.obs_grid) == K)
        grid = B.prepared.population.obs_grid
    else {
        grid = J(N, K, 0)
        for (t = 1; t <= K; t++) grid[., t] = cdid_population_rows(B.prepared.population, (1::N), J(N, 1, cal[t]))
    }
    Y = J(N, K, 0)
    O = J(N, K, 0)
    for (t = 1; t <= K; t++) {
        k = selectindex(grid[., t] :> 0)
        if (rows(k) == 0) continue
        Y[k, t] = B.prepared.population.y[grid[k, t]]
        O[k, t] = (Y[k, t] :< .)
        Y[k, t] = Y[k, t] :* O[k, t]
        Y[k, t] = editmissing(Y[k, t], 0)
    }
    p = cols(B.prepared.X)
    Zc = (B.prepared.x_has_intercept ? (p > 1 ? B.prepared.X[., 2..p] : J(N, 0, .)) : B.prepared.X)
    xn = (B.prepared.x_has_intercept ? (p > 1 ? B.prepared.x_names[2..p] : J(1, 0, "")) : B.prepared.x_names)
    /* Report settings (keepnuisance). */
    E.on = cdid_nuis_requested()
    E.cal = cal
    E.N = N
    E.clustered = B.prepared.clustered
    E.C = (B.prepared.clustered ? B.prepared.C : N)
    E.cl = (B.prepared.clustered ? B.prepared.unit_cluster : (1::N))
    E.compute_if = ctrl.inference.compute_inffunc
    E.xnames = xn
    E.method = B.prepared.method
    E.selmethod = B.prepared.selection_est_method
    E.selmodel = B.prepared.selection_model
    E.selpool = B.prepared.selection_pooling
    E.imppool = B.prepared.imputation_pooling
    E.universe = B.prepared.selection_universe
    E.snames = J(1, 0, "")
    if (corrected & cols(B.prepared.Z_sel)) E.snames = (B.prepared.sel_has_intercept ? (cols(B.prepared.sel_names) > 1 ? B.prepared.sel_names[2..cols(B.prepared.sel_names)] : J(1, 0, "")) : B.prepared.sel_names)
    _cdid_ewr_set(E)
    if (findexternal("__cdid_nuis_pending") != NULL) rmexternal("__cdid_nuis_pending")
    cohc = (q :> 0) :* tc[q + (q :== 0)]
    if (!corrected & cols(Zc)) if (missing(Zc))
        _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) requires covariates that are nonmissing for every unit.")
    if (corrected) {
        Zc = J(N, 0, .)
        for (j = 1; j <= rows(B.prepared.columns); j++) Zc = Zc, B.prepared.columns[j].unit_value
        if (missing(Zc)) _cdid_pr_fail("cdid_chained_request_input_error", "preset(cohort_time_fe) with selection() requires covariates that are nonmissing for every unit.")
        S.selmethod = B.prepared.selection_est_method
        S.selmodel = B.prepared.selection_model
        S.selpool = B.prepared.selection_pooling
        S.imppool = B.prepared.imputation_pooling
        S.universe = B.prepared.selection_universe
        S.sel_int = B.prepared.sel_has_intercept
        S.Zs = (B.prepared.sel_has_intercept ? (cols(B.prepared.Z_sel) > 1 ? B.prepared.Z_sel[., 2..cols(B.prepared.Z_sel)] : J(N, 0, .)) : B.prepared.Z_sel)
        if (S.selmethod == "imputation") {
            S.sel_int = 0
            S.Zs = J(N, 0, .)
        }
        S.Xs = Zc
        S.imp_int = B.prepared.x_has_intercept
        S.U = B.prepared.U
        S.rho = B.prepared.sel_rho
        S.eps = ctrl.baseline.glm_epsilon
        S.maxit = ctrl.baseline.glm_maxit
        if (B.prepared.method == "reg") F = cdid_etwfe_corrected_reg(Y, O, cohc, cal, never, Zc, S, ctrl.inference.compute_inffunc)
        else {
            F = cdid_etwfe_corrected_tm(Y, O, cohc, cal, never, (B.prepared.method == "dr"), Zc, S,
                ((B.settings.gmm_weighting == "identity") ? "identity" : "optimal"), ctrl.inference.compute_inffunc,
                ctrl.condition_trigger, ctrl.ridge_relative)
            F.targets[., 1] = cal[F.targets[., 1]]
            F.targets[., 2] = cal[F.targets[., 2]]
        }
    }
    else if (B.prepared.method == "reg") F = cdid_etwfe_core(Y, O, q, cal, tc, never, Zc, ctrl.inference.compute_inffunc)
    else {
        F = cdid_etwfe_treatment(Y, O, (q :> 0) :* tc[q + (q :== 0)], Zc, never, (B.prepared.method == "dr"),
            ((B.settings.gmm_weighting == "identity") ? "identity" : "optimal"), ctrl.inference.compute_inffunc,
            ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit, ctrl.condition_trigger, ctrl.ridge_relative)
        F.targets[., 1] = cal[F.targets[., 1]]
        F.targets[., 2] = cal[F.targets[., 2]]
    }
    if (E.on) _cdid_ewr_finish(B.prepared.method, (corrected ? "conditional_mtar" : "none"), B.prepared.selection_est_method)
    _cdid_ewr_clear()
    /* Reported cells. */
    B2 = B
    keep = J(rows(B.request.report), 1, 0)
    cellof = J(rows(B.request.report), 1, 0)
    for (j = 1; j <= rows(B.request.report); j++) {
        r = selectindex((F.targets[., 1] :== B.request.report[j, 1]) :& (F.targets[., 2] :== B.request.report[j, 2]))
        if (rows(r) != 1) continue
        if (B.request.report[j, 5] == 3 | (never & B.request.report[j, 5] == 2 & B.request.report[j, 7] > 0)) {
            keep[j] = 1
            cellof[j] = r[1]
        }
    }
    keep = selectindex(keep)
    cellof = cellof[keep]
    B2.request.report = B.request.report[keep, .]
    B2.request.report_reason = B.request.report_reason[keep]
    B2.request.pt_window_tag = B.request.pt_window_tag[keep]
    B2.request.pt_operator_covered = B.request.pt_operator_covered[keep]
    post_point = J(rows(B.request.post), 1, .)
    post_av = J(rows(B.request.post), 1, 0)
    post_status = J(rows(B.request.post), 1, "no_support")
    post_reason = J(rows(B.request.post), 1, "There is no observed outcome in this cohort-date cell.")
    aux_point = J(rows(B.request.aux), 1, .)
    aux_sup = J(rows(B.request.aux), 1, 0)
    aux_reason = J(rows(B.request.aux), 1, "There is no observed outcome in this cohort-date cell.")
    tgt = J(rows(B2.request.report), 1, 0)
    for (j = 1; j <= rows(B2.request.report); j++) {
        c = cellof[j]
        if (B2.request.report[j, 5] == 3) {
            h = B2.request.report[j, 6]
            if (F.targets[c, 5] == 2) {
                post_point[h] = F.targets[c, 4]
                post_av[h] = 1
                post_status[h] = "identified"
                post_reason[h] = ""
            }
            else if (F.targets[c, 5] == 1) {
                post_status[h] = "unidentified"
                post_reason[h] = "The cohort-date and common date levels are not linked by observed differences."
            }
        }
        else {
            qd = B2.request.report[j, 7]
            sg = B2.request.report[j, 8]
            if (F.targets[c, 5] == 2) {
                aux_point[qd] = sg * F.targets[c, 4]
                aux_sup[qd] = 1
                aux_reason[qd] = "supported"
            }
            else if (F.targets[c, 5] == 1) aux_reason[qd] = "The cohort-date and common date levels are not linked by observed differences."
        }
    }
    buf = cdid_public_report_start(B2, post_point, post_av, post_status, post_reason, aux_point, aux_sup, aux_reason,
        ctrl.inference.compute_inffunc, ctrl.inference.limits)
    if (ctrl.inference.compute_inffunc) {
        post_ids = selectindex(post_av)
        aux_ids = selectindex(aux_sup)
        postIF = J(N, rows(post_ids), 0)
        auxIF = J(N, rows(aux_ids), 0)
        for (j = 1; j <= rows(B2.request.report); j++) {
            c = cellof[j]
            if (F.targets[c, 5] != 2) continue
            if (B2.request.report[j, 5] == 3) {
                h = selectindex(post_ids :== B2.request.report[j, 6])
                postIF[., h] = F.IF[., c]
            }
            else {
                h = selectindex(aux_ids :== B2.request.report[j, 7])
                auxIF[., h] = B2.request.report[j, 8] :* F.IF[., c]
            }
        }
        cdid_public_report_add(buf, (1::N), post_ids, postIF, aux_ids, auxIF)
    }
    out.payload = cdid_public_report_finish(buf, ctrl.inference.limits)
    for (j = 1; j <= rows(out.payload.reason); j++) {
        c = cellof[j]
        if (out.payload.available[j]) out.payload.reason[j] = ""
        else if (F.targets[c, 5] == 1) out.payload.reason[j] = "The cohort-date and common date levels are not linked by observed differences."
        else out.payload.reason[j] = "There is no observed outcome in this cohort-date cell."
    }
    out.payload.text_settings = out.payload.text_settings \ ("preset", "cohort_time_fe") \
        ("cohort_time_fe_model", "linear unit and date fixed effects with cohort-by-date effects") \
        ("covariance", "original-unit score sandwich, no degrees-of-freedom multiplier")
    out.inference = ctrl.inference
    out.inference.report_cells = (rows(B2.request.report) ? (1::rows(B2.request.report)) : J(0, 1, .))
    out.inference.pre_cells = selectindex(B2.request.report[., 5] :== 2)
    cdid_state_validate(out.payload, ctrl.inference.limits)
    return(out)
}

void cdid_etwfe_request(real scalar on)
{
    pointer() scalar p
    p = findexternal("__cdid_etwfe_req")
    if (p == NULL) p = crexternal("__cdid_etwfe_req")
    *p = (on != 0)
}

real scalar cdid_etwfe_requested()
{
    pointer() scalar p
    p = findexternal("__cdid_etwfe_req")
    if (p == NULL) return(0)
    return(*p == 1)
}

end
