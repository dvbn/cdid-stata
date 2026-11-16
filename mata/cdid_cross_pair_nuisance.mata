*! Nuisance-model reports (keepnuisance) of cross-pair fits.
version 19.5
mata:

/* Column names of the covariate design, the constant written (Intercept). */
string rowvector _cdid_xpn_xnames(struct cdid_public_bundle scalar B)
{
    string rowvector nm
    real scalar j, p
    p = cols(B.prepared.X)
    nm = J(1, p, "")
    for (j = 1; j <= p; j++) {
        if (j <= cols(B.prepared.x_names)) nm[j] = B.prepared.x_names[j]
        else nm[j] = "x" + strofreal(j)
        if (nm[j] == "_cons") nm[j] = "(Intercept)"
    }
    return(nm)
}

string scalar _cdid_xpn_gtxt(real scalar g)
{
    return((g >= . | g == 0) ? "0" : strofreal(g, "%21.0g"))
}

/* Sorted cohort list, never-treated last, written with 0 for never treated. */
string scalar _cdid_xpn_cohorts(real colvector g)
{
    real colvector u
    string scalar s
    real scalar i
    u = uniqrows(editvalue(g, 0, .))
    s = ""
    for (i = 1; i <= rows(u); i++) s = s + (i > 1 ? "," : "") + _cdid_xpn_gtxt(u[i])
    return(s)
}

/* Right-hand side of a regression equation: intercept symbol, then b1*x1 + ... scaled by L. */
string scalar _cdid_xpn_rhs(string rowvector kept, string scalar isym, string scalar L)
{
    string scalar sl, out
    real scalar j
    sl = ""
    for (j = 1; j <= cols(kept); j++) sl = sl + (j > 1 ? " + " : "") + "b" + strofreal(j) + "*" + kept[j]
    if (sl != "" & L != "1") sl = L + "*(" + sl + ")"
    out = isym
    if (sl != "") out = out + (out != "" ? " + " : "") + sl
    return(out == "" ? "0" : out)
}

string scalar _cdid_xpn_logit_rhs(string rowvector terms)
{
    string scalar out
    real scalar j, n
    out = ""
    n = cols(terms)
    for (j = 1; j <= n; j++) {
        if (n > 7 & j == 7) {
            out = out + " + ... (" + strofreal(n - 6) + " more terms)"
            break
        }
        out = out + (j > 1 ? " + " : "") + "c" + strofreal(j - 1) + (terms[j] == "(Intercept)" ? "" : "*" + terms[j])
    }
    return(out == "" ? "0" : out)
}

/* Diagonal variance of coefficient influence functions xi (N-scaled, frame units): centered
   sums over units, or over clusters, divided by N^2. */
real rowvector _cdid_xpn_var(real matrix xi, struct cdid_public_bundle scalar B, real scalar nclus)
{
    real matrix S
    real colvector cl, hit
    real scalar N, G
    N = rows(xi)
    if (B.prepared.clustered) {
        cl = B.prepared.unit_cluster
        G = B.prepared.C
        S = J(G, cols(xi), 0)
        hit = selectindex(rowsum(xi :!= 0) :> 0)
        nclus = rows(uniqrows(cl[hit]))
        S = _cdid_xpn_csum(xi, cl, G)
        if (nclus < 2) return(J(1, cols(xi), .))
    }
    else {
        G = N
        S = xi
        nclus = .
    }
    return((colsum(S :* S) - (colsum(S) :^ 2) :/ G) :/ (N * N))
}

real matrix _cdid_xpn_csum(real matrix xi, real colvector cl, real scalar G)
{
    real matrix S, xs
    real colvector o, cs
    real scalar i, lo
    S = J(G, cols(xi), 0)
    o = order(cl, 1)
    cs = cl[o]
    xs = xi[o, ]
    lo = 1
    for (i = 1; i <= rows(cs); i++) {
        if (i == rows(cs) | cs[min((i + 1, rows(cs)))] != cs[i]) {
            S[cs[i], ] = colsum(xs[lo..i, ])
            lo = i + 1
        }
    }
    return(S)
}

void _cdid_xpn_put(struct cdid_xpn_pending scalar P)
{
    pointer() scalar p
    p = findexternal("__cdid_xpn_pending")
    if (p == NULL) p = crexternal("__cdid_xpn_pending")
    *p = P
}

void cdid_xpn_clear()
{
    if (findexternal("__cdid_xpn_pending") != NULL) rmexternal("__cdid_xpn_pending")
    if (findexternal("__cdid_xpn_base") != NULL) rmexternal("__cdid_xpn_base")
    if (findexternal("__cdid_xpn_base_info") != NULL) rmexternal("__cdid_xpn_base_info")
}

/* Capture of the imputation, treatment and counterfactual outcome fits of the reconstruction
   route (selmethod(imputation|dr)); coefficient influence functions follow every upstream
   channel the estimator propagates (repeat rates, imputation fits, observation logits). */
void cdid_xpn_recon(struct cdid_public_bundle scalar B, struct cdid_xp_cat scalar C,
    struct cdid_xps_sel scalar S, real scalar augmented, struct cdid_imputation_model scalar m,
    struct _cdid_xpr_cell colvector CE, struct _cdid_xpr_trt colvector TR, real scalar ntr,
    struct _cdid_xpr_blk colvector BK, struct _cdid_xpr_grp colvector GR, real scalar nblk,
    real colvector lc, real colvector lb, real scalar icpt, real matrix Zs, real matrix X,
    real colvector a_coh, real colvector mF_coh, string scalar method, string scalar pooling,
    string scalar imp_pooling)
{
    struct cdid_xpn_pending scalar P
    struct _cdid_xpn_rec scalar R
    struct _cdid_xpn_mat colvector IX, Aib, Aia
    struct cdid_sel_moments scalar acc
    real scalar N, L, nc, c, j, g, gid, k, b, ci, p, pk, t, nci, l, col, nr, ng
    real colvector w, A, gorder, cells, rr, xi_a, keys, rf, e, fw, f, cohs
    real matrix xi, Sb, Minv, dt, hb, H, Ath, Atha, Aab, Aaa, d, xa, cm
    real rowvector v, kp, ha
    string rowvector xn, sn, kept
    string scalar L_txt, lab
    string colvector tk
    N = rows(X)
    L = rows(C.g)
    w = B.prepared.weight.owner_weight
    A = B.prepared.population.A
    xn = _cdid_xpn_xnames(B)
    sn = (icpt ? (cols(xn) > 1 ? xn[2..cols(xn)] : J(1, 0, "")) : xn)
    P.recs = _cdid_xpn_rec(0, 1)
    nc = rows(CE)
    /* Imputation groups in order of first use by the cells. */
    gorder = J(0, 1, .)
    for (c = 1; c <= nc; c++) if (!any(gorder :== CE[c].gid)) gorder = gorder \ CE[c].gid
    IX = _cdid_xpn_mat(m.n_groups, 1)
    for (j = 1; j <= rows(gorder); j++) {
        gid = gorder[j]
        cells = J(0, 1, .)
        for (c = 1; c <= nc; c++) if (CE[c].gid == gid) cells = cells \ c
        kp = m.groups[gid].kept
        p = cols(kp)
        nci = (icpt ? rows(cells) : 0)
        Sb = J(N, p, 0)
        xi = J(N, nci + p, 0)
        for (k = 1; k <= rows(cells); k++) {
            c = cells[k]
            if (rows(CE[c].obs) == 0) continue
            e = CE[c].w_raw[CE[c].obs] :* CE[c].e_obs
            if (p) Sb[CE[c].rows_obs, ] = Sb[CE[c].rows_obs, ] + e :* CE[c].dtil_obs
            if (icpt) xi[CE[c].rows_obs, k] = xi[CE[c].rows_obs, k] + N :* e :/ CE[c].wsum
        }
        if (p) {
            Minv = _cdid_nuisance_gram(m.groups[gid].R, m.groups[gid].pivot, N :* I(p))
            xi[, (nci + 1)..(nci + p)] = Sb * Minv
            if (icpt) for (k = 1; k <= rows(cells); k++) xi[, k] = xi[, k] - xi[, (nci + 1)..(nci + p)] * CE[cells[k]].mean_d'
        }
        IX[gid].v = xi
        IX[gid].cols = cells
        R = _cdid_xpn_rec()
        R.role = 4
        R.est = J(0, 1, .)
        R.term = J(0, 1, "")
        R.crole = J(0, 1, "")
        R.cell = J(0, 1, .)
        lab = ""
        kept = (p ? sn[kp] : J(1, 0, ""))
        cohs = J(0, 1, .)
        for (k = 1; k <= rows(cells); k++) {
            c = cells[k]
            lab = lab + (k > 1 ? "; " : "") + "g=" + _cdid_xpn_gtxt(m.cohorts[CE[c].ci]) + ",(" +
                strofreal(C.cal[CE[c].cs], "%21.0g") + "," + strofreal(C.cal[CE[c].ct], "%21.0g") + ")"
            cohs = cohs \ m.cohorts[CE[c].ci]
            if (icpt) {
                R.est = R.est \ CE[c].alpha
                R.term = R.term \ "(Intercept)[g=" + _cdid_xpn_gtxt(m.cohorts[CE[c].ci]) + ",(" +
                    strofreal(C.cal[CE[c].cs], "%21.0g") + "," + strofreal(C.cal[CE[c].ct], "%21.0g") + ")]"
                R.crole = R.crole \ "local_intercept"
                R.cell = R.cell \ c
            }
        }
        for (k = 1; k <= p; k++) {
            R.est = R.est \ m.groups[gid].beta[k]
            R.term = R.term \ kept[k]
            R.crole = R.crole \ (rows(cells) > 1 ? "shared_slope" : "slope")
            R.cell = R.cell \ .
        }
        R.label = imp_pooling + "[" + lab + "]"
        L_txt = (rows(uniqrows(CE[cells[1]].L \ _cdid_xpn_Ls(CE, cells))) == 1 ? strofreal(CE[cells[1]].L) : "L(s,t)")
        if (rows(cells) == 1) {
            c = cells[1]
            R.eqn = "E[Y_i" + strofreal(C.cal[CE[c].ct], "%21.0g") + " - Y_i" + strofreal(C.cal[CE[c].cs], "%21.0g") +
                " | X_i, G_i = " + _cdid_xpn_gtxt(m.cohorts[CE[c].ci]) + ", S_i(" + strofreal(C.cal[CE[c].cs], "%21.0g") + "," +
                strofreal(C.cal[CE[c].ct], "%21.0g") + ") = 1] = " + _cdid_xpn_rhs(kept, (icpt ? "a" : ""), L_txt)
        }
        else R.eqn = "E[Y_it - Y_is | X_i, G_i = " + (rows(uniqrows(cohs)) == 1 ? _cdid_xpn_gtxt(cohs[1]) : "k") +
            ", S_i(s,t) = 1] = " + _cdid_xpn_rhs(kept, (icpt ? "a[c]" : ""), L_txt)
        R.cohorts = _cdid_xpn_cohorts(cohs)
        R.kept = invtokens(kept, ",")
        R.aliased = invtokens(_cdid_xpn_setdiff(sn, kept), ",")
        R.nrows = 0
        R.wsum = 0
        for (k = 1; k <= rows(cells); k++) {
            R.nrows = R.nrows + rows(CE[cells[k]].obs)
            R.wsum = R.wsum + CE[cells[k]].wsum
        }
        R.links = selectindex(_cdid_xpn_in(lc, cells))
        v = _cdid_xpn_var(xi, B, R.nclus)
        R.se = sqrt(v')
        P.recs = P.recs \ R
    }
    /* Treatment logits, in the order of their (target, control) cohort keys. */
    if (ntr) {
        tk = J(ntr, 1, "")
        for (t = 1; t <= ntr; t++) tk[t] = strofreal(TR[t].cig) + "|" + strofreal(TR[t].cik)
        keys = order(tk, 1)
        for (j = 1; j <= ntr; j++) {
            t = keys[j]
            p = cols(TR[t].X)
            xi = J(N, p, 0)
            xi[TR[t].rows, ] = (TR[t].w_syn :* (TR[t].y :- TR[t].prob)) :* TR[t].X
            for (k = 1; k <= 2; k++) {
                ci = (k == 1 ? TR[t].cig : TR[t].cik)
                rf = m.coh_u[ci].rows
                xi_a = (w[rf] :/ mF_coh[ci]) :* ((A[rf] :== 1) :- a_coh[ci])
                xi[rf, ] = xi[rf, ] + xi_a * TR[t].C_a[, k]'
            }
            xi = xi * invsym(TR[t].H)
            R = _cdid_xpn_rec()
            R.role = 3
            R.est = TR[t].beta
            R.term = xn'
            R.crole = J(p, 1, "slope")
            for (k = 1; k <= p; k++) if (xn[k] == "(Intercept)") R.crole[k] = "intercept"
            R.cell = J(p, 1, .)
            R.label = "risk_set[g=" + _cdid_xpn_gtxt(m.cohorts[TR[t].cig]) + " vs K={" + _cdid_xpn_gtxt(m.cohorts[TR[t].cik]) + "}]"
            R.eqn = "logit Pr(G_i = " + _cdid_xpn_gtxt(m.cohorts[TR[t].cig]) + " | X_i, G_i in {" +
                _cdid_xpn_gtxt(m.cohorts[TR[t].cig]) + "," + _cdid_xpn_gtxt(m.cohorts[TR[t].cik]) + "}) = " + _cdid_xpn_logit_rhs(xn)
            R.cohorts = _cdid_xpn_cohorts(m.cohorts[TR[t].cig] \ m.cohorts[TR[t].cik])
            R.kept = invtokens(select(xn, xn :!= "(Intercept)"), ",")
            R.aliased = ""
            R.nrows = rows(TR[t].rows)
            R.wsum = sum(TR[t].w_syn)
            R.links = J(0, 1, .)
            for (l = 1; l <= L; l++) {
                if (C.role[l] != 2) continue
                if (CE[lc[l]].ci == TR[t].cik & m.cohorts[TR[t].cig] == C.g[l]) R.links = R.links \ l
            }
            v = _cdid_xpn_var(xi, B, R.nclus)
            R.se = sqrt(v')
            P.recs = P.recs \ R
        }
    }
    /* Counterfactual outcome regressions on the reconstructed response. */
    ng = rows(GR)
    if (nblk) for (g = 1; g <= ng; g++) {
        kp = GR[g].keep
        pk = cols(kp)
        nr = rows(GR[g].blocks)
        Sb = J(N, pk, 0)
        xa = J(N, (icpt ? nr : 0), 0)
        Aab = J(m.n_cohorts, pk, 0)
        Aaa = J(m.n_cohorts, (icpt ? nr : 0), 0)
        Aib = _cdid_xpn_mat(m.n_groups, 1)
        Aia = _cdid_xpn_mat(m.n_groups, 1)
        for (gid = 1; gid <= m.n_groups; gid++) {
            Aib[gid].v = J(pk, cols(IX[gid].v), 0)
            Aia[gid].v = J((icpt ? nr : 0), cols(IX[gid].v), 0)
        }
        if (augmented) {
            acc = cdid_sel_moments_init(S.fit, pk + (icpt ? nr : 0))
            Ath = acc.G_theta
        }
        cohs = J(0, 1, .)
        for (j = 1; j <= nr; j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            ci = CE[c].ci
            cohs = cohs \ m.cohorts[ci]
            rf = CE[c].rf
            fw = CE[c].w_syn
            d = (pk ? BK[b].L :* Zs[rf, kp] : J(rows(rf), 0, .))
            e = CE[c].Z :- (BK[b].alpha :+ (pk ? d * GR[g].beta[kp] : J(rows(rf), 1, 0)))
            dt = (pk ? d :- (icpt ? BK[b].mean_d[kp] : J(1, pk, 0)) : J(rows(rf), 0, .))
            gid = CE[c].gid
            col = selectindex(IX[gid].cols :== c)
            if (pk) {
                Sb[rf, ] = Sb[rf, ] + (fw :* e) :* dt
                Aab[ci, ] = Aab[ci, ] + colsum((fw :* e) :* dt) :/ (N * CE[c].a)
                hb = (CE[c].omspi :* fw :/ N) :* dt
                if (icpt) Aib[gid].v[, col] = Aib[gid].v[, col] + colsum(hb)'
                if (cols(CE[c].d_imp)) Aib[gid].v[, (cols(IX[gid].v) - cols(CE[c].d_imp) + 1)..cols(IX[gid].v)] =
                    Aib[gid].v[, (cols(IX[gid].v) - cols(CE[c].d_imp) + 1)..cols(IX[gid].v)] + quadcross(hb, CE[c].d_imp)
            }
            if (icpt) {
                xa[rf, j] = xa[rf, j] + N :* fw :* e :/ BK[b].wsum
                Aaa[ci, j] = Aaa[ci, j] + sum(N :* fw :* e :/ BK[b].wsum) / (N * CE[c].a)
                ha = (CE[c].omspi :* fw :/ BK[b].wsum)'
                Aia[gid].v[j, col] = Aia[gid].v[j, col] + sum(ha)
                if (cols(CE[c].d_imp)) Aia[gid].v[j, (cols(IX[gid].v) - cols(CE[c].d_imp) + 1)..cols(IX[gid].v)] =
                    Aia[gid].v[j, (cols(IX[gid].v) - cols(CE[c].d_imp) + 1)..cols(IX[gid].v)] + ha * CE[c].d_imp
            }
            if (augmented & rows(CE[c].obs) & rows(CE[c].th)) {
                f = CE[c].pi_obs :* CE[c].e_obs
                if (pk) Ath[CE[c].th, 1..pk] = Ath[CE[c].th, 1..pk] + quadcross(CE[c].T, (f :* fw[CE[c].obs] :/ N) :* dt[CE[c].obs, ])
                if (icpt) Ath[CE[c].th, pk + j] = Ath[CE[c].th, pk + j] + quadcross(CE[c].T, f :* fw[CE[c].obs] :/ BK[b].wsum)
            }
        }
        /* Channels: repeat rates, imputation fits, observation logits. */
        cm = (Sb, xa)
        for (ci = 1; ci <= m.n_cohorts; ci++) {
            if (!any((Aab[ci, ], Aaa[ci, ]) :!= 0)) continue
            rf = m.coh_u[ci].rows
            xi_a = (w[rf] :/ mF_coh[ci]) :* ((A[rf] :== 1) :- a_coh[ci])
            cm[rf, ] = cm[rf, ] + xi_a * (Aab[ci, ], Aaa[ci, ])
        }
        for (gid = 1; gid <= m.n_groups; gid++) {
            if (!any((Aib[gid].v \ Aia[gid].v) :!= 0)) continue
            cm = cm + IX[gid].v * (Aib[gid].v \ Aia[gid].v)'
        }
        if (augmented) {
            if (any(Ath :!= 0)) {
                acc.G_theta = Ath
                cm = cm + cdid_sel_eval_loading(S.fit, S.rdr, acc, (1::N))
            }
        }
        xi = J(N, (icpt ? nr : 0) + pk, 0)
        if (pk) {
            Minv = N :* invsym(GR[g].M)
            xi[, ((icpt ? nr : 0) + 1)..((icpt ? nr : 0) + pk)] = cm[, 1..pk] * Minv
        }
        if (icpt) for (j = 1; j <= nr; j++) {
            xi[, j] = cm[, pk + j]
            if (pk) xi[, j] = xi[, j] - xi[, (nr + 1)..(nr + pk)] * BK[GR[g].blocks[j]].mean_d[kp]'
        }
        R = _cdid_xpn_rec()
        R.role = 2
        R.est = J(0, 1, .)
        R.term = J(0, 1, "")
        R.crole = J(0, 1, "")
        R.cell = J(0, 1, .)
        kept = (pk ? sn[kp] : J(1, 0, ""))
        lab = ""
        for (j = 1; j <= nr; j++) {
            b = GR[g].blocks[j]
            c = BK[b].cell
            if (icpt) {
                R.est = R.est \ BK[b].alpha
                R.term = R.term \ "(Intercept)[" + _cdid_xpn_bkey(CE[c], m, C) + "]"
                R.crole = R.crole \ "local_intercept"
                R.cell = R.cell \ j
            }
        }
        for (k = 1; k <= pk; k++) {
            R.est = R.est \ GR[g].beta[kp[k]]
            R.term = R.term \ kept[k]
            R.crole = R.crole \ (nr > 1 ? "shared_slope" : "slope")
            R.cell = R.cell \ .
        }
        if (nr > 1) R.label = "pooled[" + strofreal(nr) + " control block(s), common slope]"
        else R.label = "block[" + _cdid_xpn_bkey(CE[BK[GR[g].blocks[1]].cell], m, C) + "]"
        L_txt = (rows(uniqrows(_cdid_xpn_bL(BK, GR[g].blocks))) == 1 ? strofreal(BK[GR[g].blocks[1]].L) : "L(s,t)")
        R.eqn = "E[Z_i(k;s,t) | X_i, G_i = k] = " + _cdid_xpn_rhs(kept, (icpt ? (nr > 1 ? "a[b]" : "a") : ""), L_txt) + ",  " +
            (nr > 1 ? "k in K_b for each control block b = ((s,t), K_b)" : "block " + _cdid_xpn_bkey(CE[BK[GR[g].blocks[1]].cell], m, C))
        R.cohorts = _cdid_xpn_cohorts(cohs)
        R.kept = invtokens(kept, ",")
        R.aliased = invtokens(_cdid_xpn_setdiff(sn, kept), ",")
        R.nrows = 0
        R.wsum = 0
        for (j = 1; j <= nr; j++) {
            R.nrows = R.nrows + rows(CE[BK[GR[g].blocks[j]].cell].rf)
            R.wsum = R.wsum + BK[GR[g].blocks[j]].wsum
        }
        R.links = selectindex(_cdid_xpn_in(lb, GR[g].blocks) :& (lb :> 0))
        v = _cdid_xpn_var(xi, B, R.nclus)
        R.se = sqrt(v')
        P.recs = P.recs \ R
    }
    P.n = rows(P.recs)
    P.augmented = augmented
    P.route = 3
    P.method = method
    if (augmented) P.S = S
    P.linkcells = J(L, 3, .)
    for (l = 1; l <= L; l++) P.linkcells[l, ] = (editmissing(C.coh[C.ci[l]], .), C.cs[l], C.ct[l])
    _cdid_xpn_put(P)
}

real colvector _cdid_xpn_Ls(struct _cdid_xpr_cell colvector CE, real colvector cells)
{
    real colvector out
    real scalar k
    out = J(rows(cells), 1, .)
    for (k = 1; k <= rows(cells); k++) out[k] = CE[cells[k]].L
    return(out)
}

real colvector _cdid_xpn_bL(struct _cdid_xpr_blk colvector BK, real colvector blocks)
{
    real colvector out
    real scalar k
    out = J(rows(blocks), 1, .)
    for (k = 1; k <= rows(blocks); k++) out[k] = BK[blocks[k]].L
    return(out)
}

string scalar _cdid_xpn_bkey(struct _cdid_xpr_cell scalar c, struct cdid_imputation_model scalar m,
    struct cdid_xp_cat scalar C)
{
    return("k=" + _cdid_xpn_gtxt(m.cohorts[c.ci]) + ",(" + strofreal(C.cal[c.cs], "%21.0g") + "," +
        strofreal(C.cal[c.ct], "%21.0g") + ")")
}

real colvector _cdid_xpn_in(real colvector x, real colvector set)
{
    real colvector r
    real scalar i
    r = J(rows(x), 1, 0)
    for (i = 1; i <= rows(set); i++) r = r :| (x :== set[i])
    return(r)
}

string rowvector _cdid_xpn_setdiff(string rowvector a, string rowvector b)
{
    string rowvector out
    real scalar i
    out = J(1, 0, "")
    for (i = 1; i <= cols(a); i++) if (!any(b :== a[i])) out = out, a[i]
    return(out)
}


/* Upstream channel of coefficient scores Sc (rows of a cell): selection-model and repeat-rate
   coordinates through the cell's log-weight derivatives, or the cell's observation rate. */
void _cdid_xpn_cell_chan(real matrix cm, struct cdid_xps_sel scalar S, real scalar sel,
    struct _cdid_xps_cell scalar c, real matrix Sc, struct cdid_sel_moments scalar acc,
    real scalar N, real colvector wF)
{
    real rowvector a
    if (sel) {
        if (rows(c.gidx)) acc.G_theta[c.gidx, ] = acc.G_theta[c.gidx, ] + quadcross(c.ldt, Sc) :/ N
        acc.G_ak[c.kpos, ] = acc.G_ak[c.kpos, ] + quadcross(c.lda, Sc) :/ N
        return
    }
    if (!c.ract) return
    a = -colsum(Sc) :/ (N * c.rc)
    cm[c.rA, ] = cm[c.rA, ] - (c.wA :* (N / c.Mk) :* c.rraw) * a
    cm[c.f, ] = cm[c.f, ] + (wF[c.f] :* (N / c.Mk)) * a
}

/* Treatment logit of a control link on the weighted route: coefficient influence functions
   with the selection-model and repeat-rate channels (selection) or the observation-rate
   channel of the control cell (no selection model). */
struct _cdid_xpn_rec scalar cdid_xpn_wtd_trt(struct cdid_xps_sel scalar S, real scalar sel,
    struct _cdid_xps_cell scalar c, real matrix X, real colvector y, real colvector w,
    real colvector rows_all, real scalar ng, real colvector beta, real colvector prob, real matrix info,
    real scalar kg, real scalar ag, real scalar N, real scalar g, real scalar k, real scalar s,
    real scalar t, real scalar complete, real colvector wF, real scalar l)
{
    struct _cdid_xpn_rec scalar R
    struct cdid_sel_moments scalar acc
    struct cdid_public_bundle scalar B
    real matrix Sc, cm, xi
    real colvector pc
    real scalar p, j
    string rowvector xn
    pointer() scalar pb
    pb = findexternal("__cdid_xpn_bundle")
    B = *pb
    p = cols(X)
    Sc = (w :* (y :- prob)) :* X
    cm = J(N, p, 0)
    cm[rows_all, ] = cm[rows_all, ] + Sc
    pc = (ng + 1)::rows(X)
    if (sel) {
        acc = cdid_sel_moments_init(S.fit, p)
        acc.G_ak[kg, ] = acc.G_ak[kg, ] + colsum(Sc[1..ng, ]) :/ (ag * N)
        _cdid_xpn_cell_chan(cm, S, sel, c, Sc[pc, ], acc, N, wF)
        cm = cm + cdid_sel_eval_loading(S.fit, S.rdr, acc, (1::N))
    }
    else _cdid_xpn_cell_chan(cm, S, sel, c, Sc[pc, ], acc, N, wF)
    xi = cm * invsym(info :/ N)
    xn = _cdid_xpn_xnames(B)
    R.role = 3
    R.est = beta
    R.term = xn'
    R.crole = J(p, 1, "slope")
    for (j = 1; j <= p; j++) if (xn[j] == "(Intercept)") R.crole[j] = "intercept"
    R.cell = J(p, 1, .)
    R.label = "risk_set[g=" + _cdid_xpn_gtxt(g) + " vs K={" + _cdid_xpn_gtxt(k) + "}" +
        (complete ? ", complete A rows]" : ", observed (" + strofreal(s, "%21.0g") + "," + strofreal(t, "%21.0g") + ")]")
    R.eqn = "logit Pr(G_i = " + _cdid_xpn_gtxt(g) + " | X_i, (G_i = " + _cdid_xpn_gtxt(g) + (sel ? ", U_i = 1)" : ", A_i = 1)") +
        " or (G_i = " + _cdid_xpn_gtxt(k) + (complete ? ", A_i = 1)" : ", S_i(" + strofreal(s, "%21.0g") + "," + strofreal(t, "%21.0g") + ") = 1)") +
        ") = " + _cdid_xpn_logit_rhs(xn)
    R.cohorts = _cdid_xpn_cohorts(g \ k)
    R.kept = invtokens(select(xn, xn :!= "(Intercept)"), ",")
    R.aliased = ""
    R.nrows = rows(X)
    R.wsum = sum(w)
    R.links = l
    R.se = sqrt(_cdid_xpn_var(xi, B, R.nclus)')
    return(R)
}

/* Outcome-change regression group of the weighted route (separate block or pooled blocks). */
struct _cdid_xpn_rec scalar cdid_xpn_wtd_or(struct cdid_xps_sel scalar S, real scalar sel,
    struct _cdid_xps_blk colvector BK, struct _cdid_xps_cell colvector CE, struct _cdid_xps_grp scalar G,
    real matrix Xs, real scalar icpt, real scalar N, real colvector wF, real colvector links,
    struct cdid_xp_cat scalar C)
{
    struct _cdid_xpn_rec scalar R
    struct cdid_sel_moments scalar acc
    struct cdid_public_bundle scalar B
    real matrix cm, xi, d, dt, Sc
    real colvector e, fw, cohs
    real rowvector kp
    real scalar pk, nr, j, b, c, k, na
    string rowvector xn, sn, kept
    string scalar L_txt, key
    pointer() scalar pb
    pb = findexternal("__cdid_xpn_bundle")
    B = *pb
    kp = G.keep
    pk = cols(kp)
    nr = rows(G.blocks)
    na = (icpt ? nr : 0)
    cm = J(N, pk + na, 0)
    if (sel) acc = cdid_sel_moments_init(S.fit, pk + na)
    cohs = J(0, 1, .)
    for (j = 1; j <= nr; j++) {
        b = G.blocks[j]
        c = BK[b].cell
        cohs = cohs \ C.coh[CE[c].ci]
        fw = CE[c].fw
        d = (pk ? BK[b].L :* Xs[CE[c].f, kp] : J(rows(fw), 0, .))
        e = CE[c].dy :- (BK[b].alpha :+ (pk ? d * G.beta[kp] : J(rows(fw), 1, 0)))
        dt = (pk ? d :- (icpt ? BK[b].md[kp] : J(1, pk, 0)) : J(rows(fw), 0, .))
        Sc = J(rows(fw), pk + na, 0)
        if (pk) Sc[, 1..pk] = (fw :* e) :* dt
        if (icpt) Sc[, pk + j] = N :* fw :* e :/ BK[b].ws
        cm[CE[c].f, ] = cm[CE[c].f, ] + Sc
        _cdid_xpn_cell_chan(cm, S, sel, CE[c], Sc, acc, N, wF)
    }
    if (sel) cm = cm + cdid_sel_eval_loading(S.fit, S.rdr, acc, (1::N))
    xi = J(N, na + pk, 0)
    if (pk) xi[, (na + 1)..(na + pk)] = cm[, 1..pk] * (N :* invsym(G.M))
    if (icpt) for (j = 1; j <= nr; j++) {
        xi[, j] = cm[, pk + j]
        if (pk) xi[, j] = xi[, j] - xi[, (na + 1)..(na + pk)] * BK[G.blocks[j]].md[kp]'
    }
    xn = _cdid_xpn_xnames(B)
    sn = (icpt ? (cols(xn) > 1 ? xn[2..cols(xn)] : J(1, 0, "")) : xn)
    kept = (pk ? sn[kp] : J(1, 0, ""))
    R.role = 2
    R.est = J(0, 1, .)
    R.term = J(0, 1, "")
    R.crole = J(0, 1, "")
    R.cell = J(0, 1, .)
    for (j = 1; j <= nr; j++) {
        b = G.blocks[j]
        c = BK[b].cell
        key = "k=" + _cdid_xpn_gtxt(C.coh[CE[c].ci]) + ",(" + strofreal(C.cal[CE[c].cs], "%21.0g") + "," + strofreal(C.cal[CE[c].ct], "%21.0g") + ")"
        if (j == 1) R.label = key
        if (icpt) {
            R.est = R.est \ BK[b].alpha
            R.term = R.term \ "(Intercept)[" + key + "]"
            R.crole = R.crole \ "local_intercept"
            R.cell = R.cell \ j
        }
    }
    for (k = 1; k <= pk; k++) {
        R.est = R.est \ G.beta[kp[k]]
        R.term = R.term \ kept[k]
        R.crole = R.crole \ (nr > 1 ? "shared_slope" : "slope")
        R.cell = R.cell \ .
    }
    key = R.label
    R.label = (nr > 1 ? "pooled[" + strofreal(nr) + " control block(s), common slope]" : "block[" + key + "]")
    L_txt = (rows(uniqrows(_cdid_xpn_xbL(BK, G.blocks))) == 1 ? strofreal(BK[G.blocks[1]].L) : "L(s,t)")
    R.eqn = "E[Y_it - Y_is | X_i, G_i = k, S_i(s,t) = 1] = " + _cdid_xpn_rhs(kept, (icpt ? (nr > 1 ? "a[b]" : "a") : ""), L_txt) +
        ",  " + (nr > 1 ? "for each control block b = (k,(s,t))" : "block " + key)
    R.cohorts = _cdid_xpn_cohorts(cohs)
    R.kept = invtokens(kept, ",")
    R.aliased = invtokens(_cdid_xpn_setdiff(sn, kept), ",")
    R.nrows = 0
    R.wsum = 0
    for (j = 1; j <= nr; j++) {
        R.nrows = R.nrows + rows(CE[BK[G.blocks[j]].cell].f)
        R.wsum = R.wsum + BK[G.blocks[j]].ws
    }
    R.links = links
    R.se = sqrt(_cdid_xpn_var(xi, B, R.nclus)')
    return(R)
}

real colvector _cdid_xpn_xbL(struct _cdid_xps_blk colvector BK, real colvector blocks)
{
    real colvector out
    real scalar k
    out = J(rows(blocks), 1, .)
    for (k = 1; k <= rows(blocks); k++) out[k] = BK[blocks[k]].L
    return(out)
}

/* The bundle of the current fit, read by the report builders. */
void cdid_xpn_bundle(struct cdid_public_bundle scalar B)
{
    pointer() scalar p
    p = findexternal("__cdid_xpn_bundle")
    if (p == NULL) p = crexternal("__cdid_xpn_bundle")
    *p = B
}


/* Keep only the selection blocks read by the retained links (and their cells), renumbered in
   their fitted order, as R lists them after the spanning-forest reduction. */
void _cdid_xpn_sel_used(struct cdid_nuis_cap scalar cap)
{
    real colvector keep, newb, keepc, newc, kr, cc
    real scalar nb, b
    nb = rows(cap.blocks)
    if (nb == 0) return
    keep = selectindex(cap.blocks[, 13] :> 0)
    if (rows(keep) == nb) return
    newb = J(nb, 1, 0)
    if (rows(keep)) newb[keep] = (1::rows(keep))
    cap.blocks = cap.blocks[keep, ]
    cap.block_txt = cap.block_txt[keep, ]
    for (b = 1; b <= rows(keep); b++) cap.block_txt[b, 1] = "S" + strofreal(b)
    keepc = selectindex(newb[cap.cells[, 2]] :> 0)
    newc = J(rows(cap.cells), 1, 0)
    if (rows(keepc)) newc[keepc] = (1::rows(keepc))
    cap.cells = cap.cells[keepc, ]
    cap.cell_key = cap.cell_key[keepc]
    if (rows(keepc)) {
        cap.cells[, 1] = (1::rows(keepc))
        cap.cells[, 2] = newb[cap.cells[, 2]]
    }
    kr = selectindex(newb[cap.coef[, 1]] :> 0)
    cap.coef = cap.coef[kr, ]
    cap.coef_txt = cap.coef_txt[kr, ]
    if (rows(kr)) {
        cap.coef[, 1] = newb[cap.coef[, 1]]
        cc = selectindex(cap.coef[, 2] :< .)
        if (rows(cc)) cap.coef[cc, 2] = newc[cap.coef[cc, 2]]
    }
}

/* Baseline route (no selection model; analysis units): records of the outcome regressions
   and treatment logits, kept in a session external until the retained links are known. */
void cdid_xpn_base_init(string scalar method, real scalar uncond)
{
    pointer() scalar p, q
    struct _cdid_xpn_rec colvector R
    p = findexternal("__cdid_xpn_base")
    if (p == NULL) p = crexternal("__cdid_xpn_base")
    R = _cdid_xpn_rec(0, 1)
    *p = R
    q = findexternal("__cdid_xpn_base_info")
    if (q == NULL) q = crexternal("__cdid_xpn_base_info")
    *q = (method, strofreal(uncond))
}

real scalar _cdid_xpn_base_add(struct _cdid_xpn_rec scalar R)
{
    pointer() scalar p
    struct _cdid_xpn_rec colvector V
    p = findexternal("__cdid_xpn_base")
    V = *p
    V = V \ R
    *p = V
    return(rows(V))
}

void cdid_xpn_base_link(real scalar ix, real scalar l, real scalar s, real scalar t)
{
    pointer() scalar p
    struct _cdid_xpn_rec colvector V
    string scalar pr
    p = findexternal("__cdid_xpn_base")
    V = *p
    V[ix].links = V[ix].links \ l
    pr = "(" + strofreal(s, "%21.0g") + "," + strofreal(t, "%21.0g") + ")"
    if (!any(V[ix].pairs :== pr)) V[ix].pairs = V[ix].pairs \ pr
    *p = V
}

/* Variance on the analysis units (normalizer n_A); clusters of the analysis units. */
real rowvector _cdid_xpn_var_a(real matrix xi, real scalar nclus)
{
    struct cdid_public_bundle scalar B
    pointer() scalar pb
    real matrix S
    real colvector cl, hit
    real scalar n, G
    pb = findexternal("__cdid_xpn_bundle")
    B = *pb
    n = rows(xi)
    if (B.prepared.clustered) {
        cl = B.prepared.unit_cluster[B.prepared.population.A_to_F]
        G = B.prepared.C
        hit = selectindex(rowsum(xi :!= 0) :> 0)
        nclus = rows(uniqrows(cl[hit]))
        if (nclus < 2) return(J(1, cols(xi), .))
        S = _cdid_xpn_csum(xi, cl, G)
    }
    else {
        G = n
        S = xi
        nclus = .
    }
    return((colsum(S :* S) - (colsum(S) :^ 2) :/ G) :/ (n * n))
}

/* Weighted least squares of the observed control change on its observed pairs (fixed weights). */
real scalar cdid_xpn_base_or(real matrix Xk, real colvector w, real colvector e, real colvector piv,
    real matrix Rk, real colvector gam, real colvector idx, real scalar nA, struct cdid_xp_cat scalar C,
    real scalar l, real scalar p_all)
{
    struct _cdid_xpn_rec scalar R
    struct cdid_public_bundle scalar B
    pointer() scalar pb
    real matrix xi
    real colvector ord
    real scalar k, icpt
    string rowvector xn, kept, sl
    pb = findexternal("__cdid_xpn_bundle")
    B = *pb
    xn = _cdid_xpn_xnames(B)
    xi = J(nA, cols(Xk), 0)
    xi[idx, ] = nA :* (((w :* e) :* Xk) * invsym(Rk' * Rk))
    ord = order(piv, 1)
    R.role = 2
    R.est = gam[ord]
    xi = xi[, ord']
    R.term = xn[piv[ord]']'
    R.crole = J(rows(R.term), 1, "slope")
    icpt = 0
    for (k = 1; k <= rows(R.term); k++) if (R.term[k] == "(Intercept)") {
        R.crole[k] = "intercept"
        icpt = 1
    }
    R.cell = J(rows(R.term), 1, .)
    kept = select(R.term', R.term' :!= "(Intercept)")
    sl = _cdid_xpn_setdiff(xn, R.term')
    R.label = "pair[k=" + _cdid_xpn_gtxt(C.coh[C.ci[l]]) + ",(" + strofreal(C.cal[C.cs[l]], "%21.0g") + "," + strofreal(C.cal[C.ct[l]], "%21.0g") + ")]"
    R.eqn = "E[Y_i" + strofreal(C.cal[C.ct[l]], "%21.0g") + " - Y_i" + strofreal(C.cal[C.cs[l]], "%21.0g") + " | X_i, G_i = " +
        _cdid_xpn_gtxt(C.coh[C.ci[l]]) + ", S_i(" + strofreal(C.cal[C.cs[l]], "%21.0g") + "," + strofreal(C.cal[C.ct[l]], "%21.0g") +
        ") = 1] = " + _cdid_xpn_rhs(kept, (icpt ? "a" : ""), "1")
    R.cohorts = _cdid_xpn_gtxt(C.coh[C.ci[l]])
    R.kept = invtokens(kept, ",")
    R.aliased = invtokens(sl, ",")
    R.nrows = rows(idx)
    R.wsum = sum(w)
    R.links = l
    R.pairs = J(0, 1, "")
    R.se = sqrt(_cdid_xpn_var_a(xi, R.nclus)')
    return(_cdid_xpn_base_add(R))
}

/* %.12g text of a positive number, as the C library writes it. */
string scalar _cdid_xpn_g12(real scalar x)
{
    real scalar e, d
    string scalar s
    if (x == 0) return("0")
    e = floor(log10(abs(x)))
    if (e < -4 | e >= 12) return(strtrim(strofreal(x, "%21.11e")))
    d = 11 - e
    s = strtrim(strofreal(x, "%40." + strofreal(d) + "f"))
    if (strpos(s, ".")) {
        while (substr(s, strlen(s), 1) == "0") s = substr(s, 1, strlen(s) - 1)
        if (substr(s, strlen(s), 1) == ".") s = substr(s, 1, strlen(s) - 1)
    }
    return(s)
}

/* Observation-mask number of the pair (a, b) of a cohort in the numbering of the reference
   implementation (units grouped by fingerprint keys of their observed pair sets). */
real scalar _cdid_xpn_rmask(real matrix O, real scalar a, real scalar b)
{
    real scalar n, T, i, j, np, cur, h, r, hit
    real matrix cnt, fp, P, reps
    real colvector u, ord, ids
    string colvector keys, lev
    n = rows(O)
    T = cols(O)
    if (n <= 1) return(1)
    cnt = quadcross(O, O)
    u = sin((1::n) :* 1.618033988749895) :+ 2
    fp = quadcross(O, O :* u)
    P = J(0, 2, .)
    for (i = 1; i < T; i++) for (j = i + 1; j <= T; j++) if (cnt[i, j] > 0.5) P = P \ (i, j)
    np = rows(P)
    hit = 0
    for (i = 1; i <= np; i++) if (cnt[P[i, 1], P[i, 2]] != n) hit = 1
    if (!hit) return(1)
    keys = J(np, 1, "")
    for (i = 1; i <= np; i++) keys[i] = strofreal(round(cnt[P[i, 1], P[i, 2]]), "%21.0f") + ":" + _cdid_xpn_g12(fp[P[i, 1], P[i, 2]])
    lev = uniqrows(keys)
    lev = sort(lev, 1)
    ids = J(np, 1, 0)
    cur = 0
    for (h = 1; h <= rows(lev); h++) {
        ord = selectindex(keys :== lev[h])
        reps = J(n, 0, .)
        r = J(0, 1, .)
        for (i = 1; i <= rows(ord); i++) {
            hit = 0
            for (j = 1; j <= cols(reps); j++) {
                if (all(reps[, j] :== (O[, P[ord[i], 1]] :& O[, P[ord[i], 2]]))) {
                    hit = r[j]
                    break
                }
            }
            if (!hit) {
                cur++
                hit = cur
                reps = reps, (O[, P[ord[i], 1]] :& O[, P[ord[i], 2]])
                r = r \ cur
            }
            ids[ord[i]] = hit
        }
    }
    for (i = 1; i <= np; i++) if (P[i, 1] == a & P[i, 2] == b) return(ids[i])
    return(.)
}

/* Treatment logit of a control link on the baseline route: target units against the control
   cohort's units observed on the link's mask (weights w / r), with the rate channel. */
real scalar cdid_xpn_base_trt(struct _cdid_xp_lg scalar LG, real colvector w, real scalar nA,
    struct cdid_xp_cat scalar C, real scalar l, real matrix Y)
{
    struct _cdid_xpn_rec scalar R
    struct cdid_public_bundle scalar B
    pointer() scalar pb
    real matrix xi, H
    real colvector f, units, idx_c
    real scalar p, k, gi, mk
    string rowvector xn
    pb = findexternal("__cdid_xpn_bundle")
    B = *pb
    p = cols(LG.Z_u)
    xi = J(nA, p, 0)
    xi[LG.idx_u, ] = (LG.w_u :* (LG.D_u - LG.p_u)) :* LG.Z_u
    if (rows(LG.A_r)) {
        idx_c = select(LG.idx_u, LG.D_u :== 0)
        f = J(nA, 1, 0)
        f[LG.idx_all] = f[LG.idx_all] - w[LG.idx_all] :* (nA / LG.M_k) :* LG.r_raw
        f[idx_c] = f[idx_c] + w[idx_c] :* (nA / LG.M_k)
        xi = xi + f * LG.A_r'
    }
    H = LG.H_chol * LG.H_chol'
    xi = xi * invsym(H)
    xn = _cdid_xpn_xnames(B)
    if (cols(xn) != p) xn = xn[1..p]
    R.role = 3
    R.est = LG.beta
    R.term = xn'
    R.crole = J(p, 1, "slope")
    for (k = 1; k <= p; k++) if (xn[k] == "(Intercept)") R.crole[k] = "intercept"
    R.cell = J(p, 1, .)
    units = selectindex((C.ucoh :== C.coh[C.ci[l]]) :& (w :> 0))
    mk = _cdid_xpn_rmask(Y[units, ] :< ., C.cs[l], C.ct[l])
    gi = C.gi[l]
    R.g = C.g[l]
    R.k = C.coh[C.ci[l]]
    R.mask = mk
    R.skey = strofreal(gi) + ":" + strofreal(C.ci[l]) + ":" + strofreal(mk)
    R.pairs = "(" + strofreal(C.cal[C.cs[l]], "%21.0g") + "," + strofreal(C.cal[C.ct[l]], "%21.0g") + ")"
    R.cohorts = _cdid_xpn_cohorts(C.g[l] \ C.coh[C.ci[l]])
    R.kept = invtokens(select(xn, xn :!= "(Intercept)"), ",")
    R.aliased = ""
    R.nrows = rows(LG.idx_u)
    R.wsum = sum(LG.w_u)
    R.links = l
    R.se = sqrt(_cdid_xpn_var_a(xi, R.nclus)')
    return(_cdid_xpn_base_add(R))
}

/* Reports of the fits read by the retained links (used), as the capture of this call. */
void cdid_xpn_finish(struct cdid_public_bundle scalar B, real colvector used, real scalar compute_if,
    string scalar wt)
{
    struct cdid_xpn_pending scalar P
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar Rl
    struct cdid_baseline_plan_state scalar plan
    struct _cdid_xpn_rec scalar R
    pointer() scalar p
    real scalar i, role, nb, row, k, n_used, nT, ncoh, l, kpos, pc
    real colvector use, ordr
    real matrix pcl
    string scalar pre
    string rowvector info
    string colvector tk
    real scalar uncond
    p = findexternal("__cdid_xpn_pending")
    if (p == NULL) {
        /* Baseline route. */
        p = findexternal("__cdid_xpn_base")
        if (p == NULL) return
        P.recs = *p
        rmexternal("__cdid_xpn_base")
        p = findexternal("__cdid_xpn_base_info")
        info = *p
        rmexternal("__cdid_xpn_base_info")
        P.n = rows(P.recs)
        P.route = 0
        P.augmented = 0
        P.method = info[1]
        uncond = strtoreal(info[2])
        /* Treatment logits in the order of their (target, cohort, mask) keys. */
        tk = J(P.n, 1, "")
        for (i = 1; i <= P.n; i++) tk[i] = (P.recs[i].role == 3 ? P.recs[i].skey : "")
        ordr = selectindex(tk :== "") \ select(order(tk, 1), tk[order(tk, 1)] :!= "")
        P.recs = P.recs[ordr]
        for (i = 1; i <= P.n; i++) {
            if (P.recs[i].role != 3) continue
            P.recs[i].label = "risk_set[g=" + _cdid_xpn_gtxt(P.recs[i].g) + " vs k=" + _cdid_xpn_gtxt(P.recs[i].k) +
                " observed on mask " + strofreal(P.recs[i].mask) + ": pairs " + invtokens(P.recs[i].pairs', ",") + "]"
            P.recs[i].eqn = "logit Pr(G_i = " + _cdid_xpn_gtxt(P.recs[i].g) + " | X_i, G_i = " + _cdid_xpn_gtxt(P.recs[i].g) +
                " or (G_i = " + _cdid_xpn_gtxt(P.recs[i].k) + ", S_i" + P.recs[i].pairs[1] + " = 1)) = " +
                _cdid_xpn_logit_rhs(P.recs[i].term')
        }
    }
    else {
        P = *p
        rmexternal("__cdid_xpn_pending")
        uncond = 0
    }
    if (P.augmented | P.route == 2) {
        /* Observation logits, on the retained links. */
        nT = rows(B.prepared.population.calendar)
        ncoh = rows(P.S.rdr.cohorts)
        pcl = J(nT * nT, ncoh, 0)
        plan.raw = J(rows(used), 3, .)
        for (i = 1; i <= rows(used); i++) {
            l = used[i]
            kpos = _cdid_xps_kpos(P.S.rdr, P.linkcells[l, 1])
            pc = asarray(P.S.pcl, (kpos, P.linkcells[l, 2], P.linkcells[l, 3]))
            pcl[(P.linkcells[l, 2] - 1) * nT + P.linkcells[l, 3], kpos] = pc
            plan.raw[i, ] = (kpos, P.linkcells[l, 2], P.linkcells[l, 3])
        }
        plan.M = rows(used)
        plan.raw_supported = J(rows(used), 1, 1)
        plan.raw_K_ptr = (1::rows(used) + 1)
        plan.raw_K_index = plan.raw[, 1]
        plan.cohorts = P.S.rdr.cohorts
        cap = _cdid_nuis_sel_capture(P.S.fit, P.S.rdr, pcl, B.prepared.population.calendar, B, plan, 0)
        _cdid_xpn_sel_used(cap)
    }
    else {
        cap.comp = J(4, 4, "")
        cap.comp[1, ] = ("selection", "not_applicable", "no_selection_model",
            (P.route == 3 ? "This cross-pair route fits no observation-selection logit (selmethod(imputation))." :
            (P.route == 0 ? "This cross-pair route fits no observation-selection logit (baseline cross-pair moments)." :
            "This cross-pair route fits no observation-selection logit (selection(none)).")))
        cap.N = B.prepared.population.N
        cap.clustered = B.prepared.clustered
        cap.n_groups = (B.prepared.clustered ? B.prepared.C : cap.N)
        cap.comparisons = J(0, 5, "")
        cap.n_catalog = cap.n_aux = 0
    }
    cap.mode = "chained"
    cap.roles = cdid_nuis_role(1, 3)
    for (role = 2; role <= 4; role++) {
        Rl = cdid_nuis_role()
        Rl.role = (role == 2 ? "outcome" : (role == 3 ? "treatment" : "imputation"))
        Rl.blocks = J(0, 8, .)
        Rl.block_txt = J(0, 6, "")
        Rl.coef = J(0, 6, .)
        Rl.coef_txt = J(0, 2, "")
        if (P.route == 3) {
            if (role == 2) {
                Rl.title = "Counterfactual control outcome regressions on the reconstructed response Z_c"
                Rl.response = (P.augmented ? "Z_c = j_c + S_i(s,t) (Y_it - Y_is - j_c) / max(q_c, eps) on observed pairs, j_c otherwise" :
                    "Z_c = j_c, the fitted imputation, on every frame unit of the control cell")
                Rl.sample = "every frame unit of each control cell's cohort"
                Rl.weights = "unit weights times the repeat rate a_k of the control cohort"
                Rl.pooling = B.prepared.outcome_pooling
            }
            else if (role == 3) {
                Rl.title = "Treatment propensity logits of the cross-pair control links"
                Rl.response = "D_i = 1{G_i = g} on the units of the target and control cohorts"
                Rl.sample = "every frame unit of the target and control cohorts"
                Rl.weights = "unit weights times the repeat rate of the unit's cohort"
                Rl.pooling = "one logit per target and control cohort"
            }
            else {
                Rl.title = "Actual-change imputation regressions"
                Rl.response = "Y_it - Y_is on the units with both outcomes observed (S_i(s,t) = 1)"
                Rl.sample = "units of the cell's cohort with both outcomes observed"
                Rl.weights = "unit weights"
                Rl.pooling = B.prepared.imputation_pooling
            }
        }
        else {
            if (role == 2 & P.route == 0) {
                Rl.title = "Control outcome-change regressions on the observed actual change"
                Rl.response = "Delta_i(s,t) = Y_it - Y_is on control units observed at both dates"
                Rl.sample = "control-cohort units observed at both pair dates; one regression per control cohort, mask and pair"
                Rl.weights = "unit weights (fixed)"
                Rl.pooling = "separate"
            }
            else if (role == 3 & P.route == 0) {
                Rl.title = "Treatment propensity logits of the cross-pair control links"
                Rl.response = "D_i = 1{G_i = g} on the fit's risk set"
                Rl.sample = "target cohort g units and the control cohort's units observed on the link's mask"
                Rl.weights = "target: unit weights; controls: unit weights divided by the observed share r of the cohort's weight on the mask"
                Rl.pooling = "one logit per target cohort, control cohort and mask"
            }
            else if (role == 2) {
                Rl.title = "Control outcome-change regressions on the observed actual change"
                Rl.response = "Delta_i(s,t) = Y_it - Y_is on control units observed at both dates"
                Rl.sample = "control-cohort units observed at both pair dates"
                Rl.weights = (P.route == 2 ? "unit weights times a_k / (rho_i max(p_i(s,t), eps)) on observed pairs" :
                    "unit weights, divided by the observed share r of the cohort's weight on a sub-sampled pair")
                Rl.pooling = B.prepared.outcome_pooling
            }
            else if (role == 3) {
                Rl.title = "Treatment propensity logits of the cross-pair control links"
                Rl.response = "D_i = 1{G_i = g} on the fit's risk set"
                Rl.sample = "target cohort g units and the control cohort's units observed at both pair dates"
                Rl.weights = (P.route == 2 ? "target: unit weights times a_g / rho_i; controls: unit weights times a_k / (rho_i max(p_i, eps))" :
                    "unit weights, divided by the observed share r on a sub-sampled pair")
                Rl.pooling = "one logit per target cohort, control cohort and pair"
            }
        }
        if (Rl.pooling == "") Rl.pooling = "separate"
        nb = 0
        pre = (role == 2 ? "O" : (role == 3 ? "T" : "I"))
        for (i = 1; i <= P.n; i++) {
            R = P.recs[i]
            if (R.role != role) continue
            use = (rows(R.links) ? select(R.links, _cdid_xpn_in(R.links, used)) : J(0, 1, .))
            n_used = rows(use)
            if (n_used == 0) continue
            nb++
            if (!compute_if & wt != "optimal") R.se = J(rows(R.est), 1, .)
            Rl.blocks = Rl.blocks \ (rows(R.est), sum(R.crole :== "intercept" :| R.crole :== "local_intercept"),
                sum(R.crole :== "slope" :| R.crole :== "shared_slope"), R.nrows, R.wsum, R.nclus, n_used, .)
            Rl.block_txt = Rl.block_txt \ (pre + strofreal(nb), R.label, R.eqn, R.cohorts, R.kept, R.aliased)
            for (k = 1; k <= rows(R.est); k++) {
                Rl.coef = Rl.coef \ (nb, R.cell[k], R.est[k], R.se[k], R.est[k] / R.se[k], 2 * normal(-abs(R.est[k] / R.se[k])))
                Rl.coef_txt = Rl.coef_txt \ (R.term[k], R.crole[k])
            }
        }
        cap.roles[role - 1] = Rl
        if (nb) cap.comp[role, ] = (Rl.role, "available", "", "")
        else if (role == 4) cap.comp[role, ] = ("imputation", "not_applicable", "no_imputation",
            "This cross-pair route fits no actual-change imputation.")
        else if (role == 3 & P.method == "reg") cap.comp[role, ] = ("treatment", "not_applicable", "no_treatment_propensity",
            "method(reg) fits no treatment propensity.")
        else if (role == 2 & P.method == "ipw") cap.comp[role, ] = ("outcome", "not_applicable", "no_outcome_regression",
            "method(ipw) fits no outcome regression.")
        else if (uncond) cap.comp[role, ] = (Rl.role, "not_applicable", "closed_form_unconditional",
            "The design is a single constant column, so the cross-pair moments use closed-form weighted cell means; no " +
            (role == 3 ? "logit" : "regression") + " was fitted.")
        else cap.comp[role, ] = (Rl.role, "not_applicable", "no_fit_used",
            "No " + Rl.role + " fit of this route is read by a retained cross-pair link.")
    }
    _cdid_nuis_store_pending(cap)
}

/* Fits with covariates or weights that vary over time are evaluated per reference date; their
   nuisance models are not reported (as R, status unavailable, code reference_recipe_route). */
void cdid_xpn_dated(struct cdid_public_bundle scalar B, real scalar sel)
{
    struct cdid_nuis_cap scalar cap
    struct cdid_nuis_role scalar Rl
    real scalar role
    string scalar why
    if (findexternal("__cdid_xpn_pending") != NULL) rmexternal("__cdid_xpn_pending")
    if (findexternal("__cdid_xpn_base") != NULL) rmexternal("__cdid_xpn_base")
    if (findexternal("__cdid_xpn_base_info") != NULL) rmexternal("__cdid_xpn_base_info")
    why = "With covariates or weights that vary over time, the cross-pair nuisance models are fitted separately for each reference date and are not reported."
    cap.comp = J(4, 4, "")
    if (sel) cap.comp[1, ] = ("selection", "unavailable", "reference_recipe_route",
        "With covariates or weights that vary over time, the selection logits are fitted separately for each reference date and are not reported.")
    else cap.comp[1, ] = ("selection", "not_applicable", "no_selection_model",
        (B.prepared.selection_est_method == "imputation" & B.prepared.selection_canonical != "none" & B.prepared.selection_canonical != "" ?
        "This cross-pair route fits no observation-selection logit (selmethod(imputation))." :
        "This cross-pair route fits no observation-selection logit (selection(none))."))
    cap.N = B.prepared.population.N
    cap.clustered = B.prepared.clustered
    cap.n_groups = (B.prepared.clustered ? B.prepared.C : cap.N)
    cap.comparisons = J(0, 5, "")
    cap.n_catalog = cap.n_aux = 0
    cap.mode = "chained"
    cap.roles = cdid_nuis_role(1, 3)
    for (role = 2; role <= 4; role++) {
        Rl = cdid_nuis_role()
        Rl.role = (role == 2 ? "outcome" : (role == 3 ? "treatment" : "imputation"))
        Rl.blocks = J(0, 8, .)
        Rl.block_txt = J(0, 6, "")
        Rl.coef = J(0, 6, .)
        Rl.coef_txt = J(0, 2, "")
        cap.roles[role - 1] = Rl
        cap.comp[role, ] = (Rl.role, "unavailable", "reference_recipe_route", why)
    }
    _cdid_nuis_store_pending(cap)
}

end
