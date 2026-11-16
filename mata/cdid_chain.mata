*! Identification chains and GMM comparison weights of chained fits (cdid_chain).
version 19.5
mata:

/* Chain information of one fit: the comparisons used by GMM (g, s, t, role 1 target /
   2 control / 3 paired comparison, control cohort, estimate), and for every reported cell
   (g, time, ref, base, role 1 reference / 2 pre / 3 post, available, requires chaining,
   column of Lam) the GMM weights Lam on the comparisons and the witness path wit
   (cell, comparison, weight). */
struct cdid_chain_info {
    string scalar token, kind, pt, weighting
    real scalar pt_start, n_moments, n_params
    real matrix links, cells, Lam, wit
    real colvector mhat
    string colvector pool
}

string scalar _cdid_chain_lbl(real scalar x)
{
    if (x >= .) return("Inf")
    return(strtrim(strofreal(x, "%21.0g")))
}

void cdid_chain_pending(struct cdid_chain_info scalar I)
{
    pointer() scalar p
    p = findexternal("__cdid_chain_pending")
    if (p == NULL) p = crexternal("__cdid_chain_pending")
    *p = I
}

void cdid_chain_clear()
{
    if (findexternal("__cdid_chain_pending") != NULL) rmexternal("__cdid_chain_pending")
}

/* Results by token (the 64 most recent), like the nuisance reports. */
void _cdid_chain_put(string scalar token, struct cdid_chain_info scalar I)
{
    pointer() scalar p, o
    string colvector ord
    p = findexternal("__cdid_chain_registry")
    if (p == NULL) {
        p = crexternal("__cdid_chain_registry")
        *p = asarray_create()
    }
    o = findexternal("__cdid_chain_order")
    if (o == NULL) {
        o = crexternal("__cdid_chain_order")
        *o = J(0, 1, "")
    }
    I.token = token
    asarray(*p, token, I)
    ord = *o
    if (rows(ord)) ord = select(ord, ord :!= token)
    ord = ord \ token
    while (rows(ord) > 64) {
        asarray_remove(*p, ord[1])
        ord = ord[2..rows(ord)]
    }
    *o = ord
}

real scalar _cdid_chain_get(string scalar token, struct cdid_chain_info scalar I)
{
    pointer() scalar p
    p = findexternal("__cdid_chain_registry")
    if (p == NULL | token == "") return(0)
    if (!asarray_contains(*p, token)) return(0)
    I = asarray(*p, token)
    return(1)
}

real scalar cdid_chain_has(string scalar token)
{
    struct cdid_chain_info scalar I
    return(_cdid_chain_get(token, I))
}

/* Attach the pending chain information of the fit just stored to its token. */
void cdid_chain_finish(string scalar token)
{
    pointer() scalar p
    p = findexternal("__cdid_chain_pending")
    if (p == NULL) return
    _cdid_chain_put(token, *p)
    rmexternal("__cdid_chain_pending")
}

/* Chain information of a fit also describes the aggregations computed from it. */
void cdid_chain_copy(string scalar from, string scalar to)
{
    struct cdid_chain_info scalar I
    if (_cdid_chain_get(from, I)) _cdid_chain_put(to, I)
}

void cdid_chain_append(string scalar token, string scalar path)
{
    struct cdid_chain_info scalar I
    real scalar fh, rc
    if (!_cdid_chain_get(token, I)) return
    fh = _fopen(path, "a")
    if (fh < 0) _cdid_rs_fail("cannot add the chain information to the saved file", -fh)
    _cdid_rs_write(fh, ("cdid-chain-info", "1"))
    _cdid_rs_write(fh, I)
    rc = _fclose(fh)
    if (rc < 0) _cdid_rs_fail("Native result-state close failed.", -rc)
}

void cdid_chain_loaded(string scalar token, transmorphic matrix raw)
{
    if (eltype(raw) != "struct") _cdid_rs_fail("Invalid saved chain information.", 498)
    if (structname(raw) != "cdid_chain_info") _cdid_rs_fail("Invalid saved chain information.", 498)
    _cdid_chain_put(token, raw)
}

/* Cross-pair design graph in the vertex order of the R package cdid: vertices
   (target index, date index, role, control cohort) sorted by target, role, date and the
   control label; reference vertices anchor the spanning forest. */
struct _cdid_chain_graph {
    real matrix V
    real colvector from, to, anchors
    struct _cdid_forest scalar F
}

struct _cdid_chain_graph scalar _cdid_chain_xp_graph(struct cdid_xp_cat scalar C, real colvector used,
    real scalar delta, string scalar pt, real scalar pt_start, real matrix ledger)
{
    struct _cdid_chain_graph scalar Gr
    real matrix cand
    string colvector fk, ufk
    real colvector rowg, fkr, ord, keep, vix
    real rowvector ep
    real scalar L, l, i, j, k, u, gi, nv, n
    string scalar key
    transmorphic A
    L = rows(used)
    rowg = uniqrows(C.gi[used])
    cand = J(0, 4, .)
    fk = J(0, 1, "")
    for (j = 1; j <= rows(rowg); j++) {
        gi = rowg[j]
        cand = cand \ (gi, selectindex(C.cal :== C.tref[gi]), 1, 0)
        fk = fk \ ""
    }
    for (i = 1; i <= 2; i++) {
        for (l = 1; l <= L; l++) {
            u = (i == 1 ? C.cs[used[l]] : C.ct[used[l]])
            ep = _cdid_xp_endpoint(C, used[l], u, delta, pt, pt_start)
            cand = cand \ (C.gi[used[l]], u, ep)
            fk = fk \ (ep[1] == 4 ? "ctrl:" + _cdid_chain_lbl(C.coh[ep[2]]) : ((ep[1] == 2 | ep[1] == 5 | ep[1] == 7) ? "target" : ""))
        }
    }
    for (i = 1; i <= 2; i++) {
        for (j = 1; j <= rows(ledger); j++) {
            k = selectindex(C.tg :== ledger[j, 1])
            if (rows(k) != 1) continue
            if (!any(rowg :== k)) continue
            u = selectindex(C.cal :== ledger[j, 2])
            if (rows(u) != 1) continue
            cand = cand \ (k, u, (i == 1 ? 7 : 6), 0)
            fk = fk \ (i == 1 ? "target" : "")
        }
    }
    n = rows(cand)
    A = asarray_create()
    keep = J(n, 1, 0)
    vix = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        key = strofreal(cand[i, 1]) + "|" + strofreal(cand[i, 2]) + "|" + strofreal(cand[i, 3]) + "|" + fk[i]
        if (!asarray_contains(A, key)) {
            asarray(A, key, i)
            keep[i] = 1
        }
        vix[i] = asarray(A, key)
    }
    ufk = uniqrows(fk)
    fkr = J(n, 1, 0)
    for (i = 1; i <= n; i++) fkr[i] = selectindex(ufk :== fk[i])
    keep = selectindex(keep)
    ord = keep[order((cand[keep, 1], cand[keep, 3], C.cal[cand[keep, 2]], fkr[keep], keep), (1, 2, 3, 4, 5))]
    nv = rows(ord)
    Gr.V = cand[ord, .]
    k = J(n, 1, 0)
    k[ord] = (1::nv)
    vix = k[vix]
    Gr.from = vix[(rows(rowg) + 1)::(rows(rowg) + L)]
    Gr.to = vix[(rows(rowg) + L + 1)::(rows(rowg) + 2 * L)]
    Gr.anchors = selectindex(Gr.V[, 3] :== 1)
    Gr.F = _cdid_bfs_forest(nv, Gr.from, Gr.to, Gr.anchors, (1::L))
    return(Gr)
}

real scalar _cdid_chain_vertex(struct _cdid_chain_graph scalar Gr, real scalar gi, real scalar u,
    real rowvector roles)
{
    real colvector h
    real scalar r
    h = J(0, 1, .)
    for (r = 1; r <= cols(roles); r++) h = h \ selectindex((Gr.V[, 1] :== gi) :& (Gr.V[, 2] :== u) :& (Gr.V[, 3] :== roles[r]))
    return(rows(h) == 1 ? h[1] : 0)
}

/* Witness of a contrast (vertices, coefficients): the sum of the coefficients times the
   root paths of the vertices in the spanning forest; empty when not identified. */
real matrix _cdid_chain_witness(struct _cdid_chain_graph scalar Gr, real colvector v, real colvector c, real scalar L)
{
    real colvector w, unanch
    real scalar i, cur, comp
    if (any(v :== 0)) return(J(0, 2, .))
    unanch = J(Gr.F.nc, 1, 1)
    for (i = 1; i <= rows(Gr.anchors); i++) unanch[Gr.F.comp[Gr.anchors[i]]] = 0
    for (comp = 1; comp <= Gr.F.nc; comp++) {
        if (!unanch[comp]) continue
        if (sum(c :* (Gr.F.comp[v] :== comp)) != 0) return(J(0, 2, .))
    }
    w = J(L, 1, 0)
    for (i = 1; i <= rows(v); i++) {
        cur = v[i]
        while (Gr.F.par[cur] != 0) {
            w[Gr.F.pe[cur]] = w[Gr.F.pe[cur]] + c[i] * Gr.F.ps[cur]
            cur = Gr.F.par[cur]
        }
    }
    i = selectindex(w :!= 0)
    if (rows(i) == 0) return(J(0, 2, .))
    return((i, w[i]))
}

/* Contrast of cell (target index j, date u, base from, reference refi) in the graph. */
void _cdid_chain_xp_contrast(struct _cdid_chain_graph scalar Gr, real scalar j, real scalar u,
    real scalar from, real scalar refi, real scalar post, real colvector v, real colvector c)
{
    if (post) {
        v = _cdid_chain_vertex(Gr, j, u, 7) \ _cdid_chain_vertex(Gr, j, u, 6)
        c = (1 \ -1)
        return
    }
    v = J(0, 1, .)
    c = J(0, 1, .)
    if (u != refi) {
        v = v \ _cdid_chain_vertex(Gr, j, u, (2, 5)) \ _cdid_chain_vertex(Gr, j, u, (3, 6))
        c = c \ 1 \ -1
    }
    if (from != refi) {
        v = v \ _cdid_chain_vertex(Gr, j, from, (2, 5)) \ _cdid_chain_vertex(Gr, j, from, (3, 6))
        c = c \ -1 \ 1
    }
}

/* Chain information of a cross-pair fit for the reported cells. */
struct cdid_chain_info scalar cdid_chain_xp(struct cdid_xp_cat scalar C, struct cdid_xp_fit scalar F,
    real colvector used, real scalar delta, string scalar pt, real scalar pt_start, string scalar base,
    string scalar weighting, real matrix report, real matrix ledger)
{
    struct cdid_chain_info scalar I
    struct _cdid_chain_graph scalar Gm, Ga
    real scalar L, r, c, j, u, refi, from, post, lo, hi, direct, k, aux
    real colvector v, cf
    real matrix W
    L = rows(used)
    I.kind = "crosspair"
    I.pt = pt
    I.pt_start = pt_start
    I.weighting = weighting
    I.n_moments = L
    I.links = (C.g[used], C.cal[C.cs[used]], C.cal[C.ct[used]], C.role[used], C.coh[C.ci[used]])
    I.mhat = F.mhat
    I.pool = J(L, 1, "")
    for (k = 1; k <= L; k++) {
        I.pool[k] = (C.role[used[k]] == 1 ? "target:" + _cdid_chain_lbl(C.g[used[k]]) : "{" + _cdid_chain_lbl(C.coh[C.ci[used[k]]]) + "}")
    }
    Gm = _cdid_chain_xp_graph(C, used, delta, pt, pt_start, ledger)
    Ga = _cdid_chain_xp_graph(C, used, delta, "ctrl", min(C.cal), ledger)
    I.n_params = 0
    for (k = 1; k <= Gm.F.nc; k++) {
        if (any(Gm.F.comp[Gm.anchors] :== k)) I.n_params = I.n_params + sum((Gm.F.comp :== k) :& (Gm.V[., 3] :!= 1))
        else I.n_params = I.n_params + sum(Gm.F.comp :== k) - 1
    }
    I.cells = J(rows(report), 8, 0)
    I.Lam = J(L, 0, .)
    I.wit = J(0, 3, .)
    for (r = 1; r <= rows(report); r++) {
        I.cells[r, 1..5] = report[r, 1..5]
        c = _cdid_xp_cell(F, report[r, 1], report[r, 2])
        if (!c) continue
        I.cells[r, 6] = F.avail[c]
        if (!F.avail[c] | report[r, 5] == 1) continue
        j = selectindex(C.tg :== report[r, 1])
        u = selectindex(C.cal :== report[r, 2])
        if (rows(j) != 1 | rows(u) != 1) continue
        refi = selectindex(C.cal :== C.tref[j])
        if (rows(refi) != 1) continue
        post = F.cells[c, 3] & (C.cal[u] >= C.tg[j])
        from = (post ? . : (base == "universal" ? refi : u - 1))
        aux = any(F.cp :== c)
        if (aux | any(F.cm :== c)) {
            I.Lam = I.Lam, F.Lam[., c]
            I.cells[r, 8] = cols(I.Lam)
            if (aux) _cdid_chain_xp_contrast(Ga, j, u, from, refi, post, v, cf)
            else _cdid_chain_xp_contrast(Gm, j, u, from, refi, post, v, cf)
            W = (aux ? _cdid_chain_witness(Ga, v, cf, L) : _cdid_chain_witness(Gm, v, cf, L))
            if (rows(W)) I.wit = I.wit \ (J(rows(W), 1, r), W)
        }
        lo = min((report[r, 4], report[r, 2]))
        hi = max((report[r, 4], report[r, 2]))
        if (report[r, 4] >= . | lo == hi) direct = 1
        else direct = any((C.role :== 1) :& (C.g :== C.tg[j]) :& (C.cal[C.cs] :== lo) :& (C.cal[C.ct] :== hi)) &
            any((C.role :== 2) :& (C.g :== C.tg[j]) :& (C.cal[C.cs] :== lo) :& (C.cal[C.ct] :== hi))
        I.cells[r, 7] = !direct
    }
    return(I)
}


/* Catalogue restricted to rows used (ascending), for the design of the comparisons kept by
   an exact reduction. */
struct cdid_cat scalar _cdid_chain_subcat(struct cdid_cat scalar C, real colvector used)
{
    struct cdid_cat scalar S
    S = C
    S.r_g = C.r_g[used]
    S.r_s = C.r_s[used]
    S.r_t = C.r_t[used]
    S.r_op = C.r_op[used]
    S.r_ref = C.r_ref[used]
    S.r_srole = C.r_srole[used]
    S.r_trole = C.r_trole[used]
    S.r_tp = C.r_tp[used]
    S.r_pool = C.r_pool[used]
    return(S)
}

string scalar _cdid_chain_poollbl(struct cdid_cat scalar C, real colvector K)
{
    string scalar out
    real scalar k
    out = "{"
    for (k = 1; k <= rows(K); k++) out = out + (k > 1 ? "," : "") + _cdid_chain_lbl(C.cohorts[K[k]])
    return(out + "}")
}

/* Chain information of a same-pair fit. used: catalogue rows of the GMM comparisons
   (ascending); Du: their design; Lam: GMM weights of the estimable post targets (post
   rows est_ids) on used; mhat: estimates of used. Pre-treatment cells are single
   auxiliary comparisons (aux descriptors with values aux_point). */
struct cdid_chain_info scalar cdid_chain_sp(struct cdid_cat scalar C, struct cdid_design scalar Dall,
    real colvector used, struct cdid_design scalar Du, real colvector mhat, real matrix Lam,
    real colvector est_ids, real matrix report, real matrix aux, real colvector auxK_ptr,
    real colvector auxK_idx, real colvector aux_point, real colvector aux_sup, string scalar pt,
    real scalar pt_start, string scalar weighting)
{
    struct cdid_chain_info scalar I
    real scalar L, r, h, q, k, col, lo, hi, gi, nl
    real matrix W
    real colvector K
    L = rows(used)
    I.kind = "samepair"
    I.pt = pt
    I.pt_start = pt_start
    I.weighting = weighting
    I.n_moments = L
    I.n_params = Du.n_rank
    I.links = (C.cohorts[C.r_g[used]], C.calendar[C.r_s[used]], C.calendar[C.r_t[used]], J(L, 1, 3), J(L, 1, .))
    I.mhat = mhat
    I.pool = J(L, 1, "")
    for (k = 1; k <= L; k++) {
        K = C.pool_entries[|C.pool_start[C.r_pool[used[k]]] \ C.pool_start[C.r_pool[used[k]]] + C.pool_size[C.r_pool[used[k]]] - 1|]
        I.pool[k] = _cdid_chain_poollbl(C, K)
    }
    I.cells = J(rows(report), 8, 0)
    I.Lam = J(L, 0, .)
    I.wit = J(0, 3, .)
    for (r = 1; r <= rows(report); r++) {
        I.cells[r, 1..5] = report[r, 1..5]
        if (report[r, 5] == 1) {
            I.cells[r, 6] = 1
            continue
        }
        if (report[r, 5] == 2) {
            q = report[r, 7]
            if (q < 1 | q >= . | q > rows(aux_sup)) continue
            if (!aux_sup[q]) continue
            I.cells[r, 6] = 1
            K = auxK_idx[|auxK_ptr[q] \ auxK_ptr[q + 1] - 1|]
            I.links = I.links \ (aux[q, 1], aux[q, 2], aux[q, 3], 4, .)
            I.mhat = I.mhat \ aux_point[q]
            I.pool = I.pool \ _cdid_chain_poollbl(C, K)
            I.wit = I.wit \ (r, rows(I.links), report[r, 8])
            continue
        }
        h = report[r, 6]
        if (h < 1 | h >= .) continue
        col = selectindex(est_ids :== h)
        if (rows(col) != 1) continue
        I.cells[r, 6] = 1
        I.Lam = I.Lam, Lam[col, .]'
        I.cells[r, 8] = cols(I.Lam)
        if (Du.avail[h, 7] >= 1 & Du.avail[h, 7] < .) {
            W = cdid_design_witness(Du, Du.avail[h, 7])
            if (rows(W)) I.wit = I.wit \ (J(rows(W), 1, r), W)
        }
        lo = min((report[r, 4], report[r, 2]))
        hi = max((report[r, 4], report[r, 2]))
        if (report[r, 4] < . & lo != hi) {
            gi = _cdid_bsearch(C.cohorts, report[r, 1])
            nl = _cdid_bsearch(C.calendar, lo)
            k = _cdid_bsearch(C.calendar, hi)
            I.cells[r, 7] = !any(Dall.l_tq :& (Dall.l_g :== gi) :& (Dall.l_s :== nl) :& (Dall.l_t :== k))
        }
    }
    return(I)
}

/* Records the chain information of a same-pair fit for cdid_chain. kept: catalogue rows of
   the GMM comparisons (ascending); Lam: GMM weights of the estimable post targets on the full
   catalogue (zero outside kept); mhat: estimates of the full catalogue. */
void cdid_chain_sp_record(struct cdid_public_bundle scalar B, struct cdid_cat scalar C,
    struct cdid_design scalar D, real colvector kept, real colvector mhat, real matrix Lam,
    real colvector est_ids, real colvector aux_point, real colvector aux_sup)
{
    struct cdid_design scalar Dk
    real scalar M
    M = rows(C.r_g)
    if (rows(kept) < M) {
        Dk = cdid_chained_ipw_build_design(_cdid_chain_subcat(C, kept), B.request.post,
            B.request.pt_assumption, B.request.pt_start)
        cdid_chain_pending(cdid_chain_sp(C, D, kept, Dk, mhat[kept], Lam[., kept'], est_ids,
            B.request.report, B.request.aux, B.request.aux_K_ptr, B.request.aux_K_index, aux_point,
            aux_sup, B.request.pt_assumption, B.request.pt_start, B.settings.gmm_weighting))
    }
    else cdid_chain_pending(cdid_chain_sp(C, D, (1::M), D, mhat, Lam, est_ids, B.request.report,
        B.request.aux, B.request.aux_K_ptr, B.request.aux_K_index, aux_point, aux_sup,
        B.request.pt_assumption, B.request.pt_start, B.settings.gmm_weighting))
}

/* Values of the supported auxiliary comparisons (missing elsewhere). */
real colvector cdid_chain_auxpts(real colvector points, real colvector aux_consumer, real colvector aux_supported)
{
    real colvector ap, gq
    ap = J(rows(aux_consumer), 1, .)
    gq = selectindex(aux_supported :!= 0)
    if (rows(gq)) ap[gq] = points[aux_consumer[gq]]
    return(ap)
}

end
