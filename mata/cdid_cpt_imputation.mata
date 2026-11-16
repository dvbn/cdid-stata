*! Private pure-imputation and augmented conditional-pretrend process piece producer; StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (commit ):
*! R/chained_conditional_extract.R cpt_piece (lines 16-27) and cpt_extract_imputation
*! (lines 347-509, kind "imputation" and kind "augmented"); the coordinate IF rows follow
*! R/chained_imputation.R imputation_evaluate_if_chunk_core details mode (xi_a, imputation
*! alpha/beta, treatment gamma, shared OR alpha/beta) and the native augmented
*! backend (Z-response OR with its q channel, theta-only selection contractions).
*! Borrows a live cdid_imp_backend (static or dated; ipw, reg, dr) or a live
*! cdid_aug_backend (dr, conditional_mtar, F, rho = 1), already fitted on the
*! complete original post + both-auxiliary union. Nothing is refit by threshold or chunk;
*! no h * final ATT IF is substituted; no N-by-all-moment or parameter IF state is kept.
*! Native conditional pretrend adapter for imputation and augmented correction.
*! Load after cdid_conditional, cdid_cpt_baseline, cdid_public_request,
*! cdid_imputation_model, cdid_imputation_dated, cdid_imputation_backend and
*! cdid_augmented_backend (and their own prerequisites).
version 19.5
mata:
mata set matastrict on

/* One declared pretrend cell of the requested family. terms are source term codes
   1 = T, 2 = Bg, 3 = C, 4 = Bo in source order. a_ci/cids are the cohort and Stage 1 cell
   coordinates (target first, then link controls for ipw/dr); gids are the cells' p > 0
   Stage 1 groups (first appearance) with column offsets goff; cell_gpos maps each cid to
   its gids position (0 when the group has no kept slope); th holds sorted native theta IDs
   (augmented only). adj_alpha is the bounded Stage 1 adjoint of the cell's OR block alpha
   j forcing (jf_alpha_alpha, jf_alpha_beta). */
struct cdid_cpti_cell_meta {
    string scalar id, status, reason, route
    real scalar g, t, from, q, sign, link, recipe, supported
    real scalar s_date, t2_date, estimate, estimate_scalar, n_chunks, chunk_rows
    real scalar n_terms, Dl, icpt, use_or, sep_or, p_trt, p_or, tid, bid
    real scalar off_ia, off_ib, off_ps, off_oa, off_ob, off_th, work_row
    real colvector terms, a_ci, cids, gids, goff, cell_gpos, th
    struct cdid_cpt_spec scalar spec
    struct cdid_imputation_adjoint scalar adj_alpha
}

/* Piece producer. ib/ab are BORROWED live owners: the caller keeps them (and
   everything they borrow) alive and unmodified; every replay re-runs the owner
   checks plus this producer's fingerprint. kind is "imputation" or "augmented"; layout
   is "static" or "dated". orig_link (M x 1) and aux_link (Q x 1) are the verified backend
   link of every primitive original comparison and every auxiliary comparison. */
struct cdid_cpti {
    string scalar kind, layout, method, base_period, status, reason
    real scalar n_F, n_A, G, n_recipes, n_links, n_primitive, n_post_targets, n_aux
    real scalar n_cells, n_unavailable, chunk_rows
    real scalar retained_bytes, workspace_limit, work_ops, union_work, work_limit, fingerprint
    pointer(struct cdid_imp_backend scalar) scalar ib
    pointer(struct cdid_aug_backend scalar) scalar ab
    real colvector orig_link, aux_link
    struct cdid_cpti_cell_meta colvector cells
}

void _cdid_cpti_fail(string scalar cls, string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}

string scalar _cdid_cpti_istr(real scalar x)
{
    return(strtrim(strofreal(x, "%21.0f")))
}

/* Request never-treated cohort (.) as the Stage 1 value 0. */
real colvector _cdid_cpti_cval(real colvector v)
{
    real colvector out
    out = v
    if (any(out :>= .)) out[selectindex(out :>= .)] = J(sum(out :>= .), 1, 0)
    return(out)
}

/* Exact R cpt_piece(rows, nu, rho, lam, pim, sD, xi) for R >= 1 terms. lam and pim are
   rows x (R Dl), term-major; every term carries both loadings (possibly zero). */
struct cdid_cpt_piece scalar _cdid_cpti_piece_make(real colvector fr, real matrix nu,
    real matrix rho, real matrix lam, real matrix pim, real colvector sD, real matrix xi,
    struct cdid_cpt_spec scalar spec)
{
    struct cdid_cpt_piece scalar pc
    real scalar n, R, Dl, r
    real rowvector cc
    n = rows(fr); R = spec.R; Dl = spec.Dl
    if (R < 1 | cols(fr) != 1 | rows(nu) != n | cols(nu) != R | rows(rho) != n | cols(rho) != R |
        rows(xi) != n | cols(xi) != Dl | rows(lam) != n | cols(lam) != R * Dl |
        rows(pim) != n | cols(pim) != R * Dl | rows(sD) != R | cols(sD) != 1)
        _cdid_cpti_fail("piece", "Invalid imputation process piece dimensions.", 498)
    pc.rows = fr
    pc.nu = nu
    pc.rho = rho
    pc.xi = xi
    pc.spec = spec
    pc.P = J(n, Dl, 0)
    pc.blam = J(Dl, R, 0)
    if (Dl == 0 | n == 0) return(pc)
    for (r=1; r<=R; r++) {
        cc = ((r - 1) * Dl + 1)..(r * Dl)
        pc.P = pc.P + ((sD[r] :* nu[, r] :* rho[, r]) :* lam[, cc])
        pc.blam[, r] = colsum(nu[, r] :* lam[, cc])'
        pc.P = pc.P + ((sD[r] :* nu[, r]) :* pim[, cc])
    }
    if (missing(pc.nu) | missing(pc.rho) | missing(pc.P) | missing(pc.xi) | missing(pc.blam))
        _cdid_cpti_fail("numeric", "Nonfinite imputation process piece for cell " + spec.id + ".", 504)
    return(pc)
}

/* Observation reader handed to _cdid_aug_cell_rows: the static borrowed reader, or the
   dated lease (the dated route reads a.qlease internally either way). */
pointer(struct cdid_sel_reader scalar) scalar _cdid_cpti_aug_qr(struct cdid_aug_backend scalar a)
{
    if (a.is_dated) return(&(a.qlease))
    return(a.qrdr)
}

/* Exact fingerprint of the borrowed backend record (not a tolerance). */
real scalar _cdid_cpti_fp(struct cdid_cpti scalar ip)
{
    pointer(struct cdid_imp_backend scalar) scalar pb
    pointer(struct cdid_aug_backend scalar) scalar pa
    real scalar s, ell, bid
    if (ip.kind == "augmented") {
        pa = ip.ab
        pb = (*pa).ipw
    }
    else pb = ip.ib
    s = (*pb).N + 3 * (*pb).n_links + 5 * (*pb).n_blocks + 7 * (*pb).n_trt + 11 * (*pb).or_p + (*pb).borrow_sum
    if ((*pb).has_or & (*pb).outcome_pooling == "separate") {
        for (bid=1; bid<=(*pb).n_blocks; bid++) s = s + (11 + bid) * (*pb).blocks[bid].or_p
    }
    for (ell=1; ell<=(*pb).n_links; ell++) s = s + ell * (*pb).links[ell].mu
    if (ip.kind == "augmented") {
        s = s + (*pa).ipw_sum + 13 * (*pa).or_p + 17 * (*pa).n_blocks + 19 * (*pa).n_q_params
        if ((*pa).has_or & (*pa).outcome_pooling == "separate") {
            for (bid=1; bid<=(*pa).n_blocks; bid++) s = s + (13 + bid) * (*pa).blocks[bid].or_p
        }
        for (ell=1; ell<=(*pa).n_links; ell++) s = s + (ell + 0.5) * (*pa).links[ell].mu
    }
    return(s)
}

/* borrowed-owner checks, then this producer's record fingerprint. */
void _cdid_cpti_check(struct cdid_cpti scalar ip)
{
    pointer(struct cdid_imputation_model scalar) scalar pm
    if (ip.kind == "augmented") {
        if (ip.ab == NULL)
            _cdid_cpti_fail("owner", "Borrowed augmented backend pointer is NULL; its lifetime ended.", 498)
        _cdid_aug_check(*(ip.ab))
    }
    else {
        if (ip.ib == NULL)
            _cdid_cpti_fail("owner", "Borrowed imputation backend pointer is NULL; its lifetime ended.", 498)
        _cdid_impb_need_model(*(ip.ib))
        pm = (*(ip.ib)).model
        _cdid_impb_check(*(ip.ib), *pm)
    }
    if (_cdid_cpti_fp(ip) != ip.fingerprint)
        _cdid_cpti_fail("owner", "Borrowed backend record changed after conditional preparation; prepare again.", 498)
}

/* Complete fitted union: rebuild the original comparison catalogue once (metadata
   only) and match every primitive original row and every auxiliary row to a backend link
   with the same calendar endpoints, target cohort and ordered control set; every backend
   link must be consumed. Supplied maps are verified; empty maps are derived and must be
   unique (dated multi-recipe unions need the caller's explicit recipe-aware maps). */
void _cdid_cpti_union(struct cdid_cpti scalar ip, struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, struct cdid_public_request scalar req,
    real colvector orig_link, real colvector aux_link)
{
    struct cdid_catalog_opts scalar opts
    struct cdid_cat scalar C
    real matrix occ
    real scalar M, Q, i, j, ell, n_hit, hit, gv, sv, tv, pid, given
    real colvector K, kv, used, lg, ls, lt, map
    opts = cdid_catalog_default_opts()
    opts.mem_limit_bytes = ip.workspace_limit
    opts.work_cap = ip.work_limit
    occ = _cdid_cpt_base_occ(m.pop, req.cohorts, opts)
    C = cdid_catalog_build_opts(req.calendar, req.cohorts, req.target_cohorts,
        req.anticipation, occ, req.control_group, req.post, opts)
    M = rows(C.r_g)
    Q = rows(req.aux)
    if (rows(C.cohorts) != rows(req.cohorts) | !all(C.cohorts :== req.cohorts))
        _cdid_cpti_fail("identity", "Catalogue and request cohort dictionaries disagree.", 498)
    ip.n_primitive = M
    ip.n_aux = Q
    ip.n_post_targets = rows(req.post)
    ip.work_ops = ip.work_ops + C.work + (M + Q) * (b.n_links + 4)
    ip.union_work = ip.work_ops
    if (ip.work_limit < . & ip.work_ops > ip.work_limit)
        _cdid_cpti_fail("bound", "Union identity check exceeds the explicit counted-work cap.", 908)
    if (rows(orig_link) != 0 & (rows(orig_link) != M | cols(orig_link) != 1))
        _cdid_cpti_fail("identity", "orig_link must be empty or one backend link per primitive original comparison.", 498)
    if (rows(aux_link) != 0 & (rows(aux_link) != Q | cols(aux_link) != 1))
        _cdid_cpti_fail("identity", "aux_link must be empty or one backend link per auxiliary comparison.", 498)
    lg = J(b.n_links, 1, 0); ls = J(b.n_links, 1, 0); lt = J(b.n_links, 1, 0)
    for (ell=1; ell<=b.n_links; ell++) {
        pid = b.links[ell].pair_id
        lg[ell] = m.cohorts[b.links[ell].target_ci]
        ls[ell] = m.pop.calendar[m.pairs[pid].s_index]
        lt[ell] = m.pop.calendar[m.pairs[pid].t_index]
    }
    map = J(M + Q, 1, 0)
    used = J(b.n_links, 1, 0)
    for (i=1; i<=M + Q; i++) {
        if (i <= M) {
            gv = C.cohorts[C.r_g[i]]; sv = C.calendar[C.r_s[i]]; tv = C.calendar[C.r_t[i]]
            j = C.r_pool[i]
            if (C.pool_size[j] < 1)
                _cdid_cpti_fail("identity", "A primitive original comparison has an empty control pool.", 498)
            K = C.pool_entries[C.pool_start[j]..C.pool_start[j] + C.pool_size[j] - 1]
            given = .
            if (rows(orig_link)) given = orig_link[i]
        }
        else {
            j = i - M
            gv = req.aux[j, 1]; sv = req.aux[j, 2]; tv = req.aux[j, 3]
            if (req.aux_K_ptr[j + 1] <= req.aux_K_ptr[j])
                _cdid_cpti_fail("identity", "An auxiliary comparison has an empty control pool.", 498)
            K = req.aux_K_index[req.aux_K_ptr[j]..req.aux_K_ptr[j + 1] - 1]
            given = .
            if (rows(aux_link)) given = aux_link[j]
        }
        kv = sort(_cdid_cpti_cval(req.cohorts[K]), 1)
        n_hit = 0
        hit = 0
        if (given < .) {
            if (given >= 1 & given <= b.n_links & given == floor(given)) {
                ell = given
                if (lg[ell] == (gv >= . ? 0 : gv) & ls[ell] == sv & lt[ell] == tv &
                    b.links[ell].n_ctrl == rows(kv)) {
                    if (all(sort(m.cohorts[b.links[ell].ctrl_ci], 1) :== kv)) hit = ell
                }
            }
            if (hit == 0)
                _cdid_cpti_fail("identity", "A supplied union map entry does not name a backend link with the same endpoints, target cohort and control set.", 498)
        }
        else {
            for (ell=1; ell<=b.n_links; ell++) {
                if (lg[ell] != (gv >= . ? 0 : gv) | ls[ell] != sv | lt[ell] != tv) continue
                if (b.links[ell].n_ctrl != rows(kv)) continue
                if (!all(sort(m.cohorts[b.links[ell].ctrl_ci], 1) :== kv)) continue
                n_hit++
                if (hit == 0) hit = ell
            }
            if (n_hit == 0)
                _cdid_cpti_fail("identity", "The backend is not fitted on the complete original post + both-auxiliary union of the request.", 498)
            if (n_hit > 1)
                _cdid_cpti_fail("identity", "Several backend links (design blocks) match one comparison; supply explicit block-specific orig_link/aux_link maps.", 498)
        }
        map[i] = hit
        used[hit] = 1
    }
    if (!all(used))
        _cdid_cpti_fail("identity", "The backend carries links outside the request's original post + both-auxiliary union; its shared fits are not this request's union fit.", 498)
    ip.orig_link = J(0, 1, 0)
    ip.aux_link = J(0, 1, 0)
    if (M > 0) ip.orig_link = map[1..M]
    if (Q > 0) ip.aux_link = map[(M + 1)..(M + Q)]
}

/* Coordinates, terms, spec, OR-alpha adjoint, chunking and counted work of one available
   cell, read from the live backend (source cpt_extract_imputation lines 361-396). */
void _cdid_cpti_layout(struct cdid_cpti scalar ip, struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, struct cdid_cpti_cell_meta scalar cm)
{
    pointer(struct cdid_aug_backend scalar) scalar pa
    real scalar ell, aug, ipwdr, cp, gid, gp, j, at, nk, rl, nth, R, Dl, per_row, fixed, Bmax
    real scalar d_plus, d_minus, k, cl, p, nblk
    real colvector svec, Dvec, ids, s_all, D_all
    real matrix lamA, nuA
    string colvector coords
    aug = (ip.kind == "augmented")
    if (aug) pa = ip.ab
    ell = cm.link
    rl = b.links[ell].recipe
    nk = b.links[ell].n_ctrl
    ipwdr = (ip.method == "ipw" | ip.method == "dr")
    cm.use_or = (ip.method == "reg" | ip.method == "dr")
    cm.icpt = m.intercept
    cm.recipe = rl
    cm.terms = 1
    if (cm.use_or) cm.terms = cm.terms \ 2
    if (ipwdr) cm.terms = cm.terms \ 3
    if (ip.method == "dr") cm.terms = cm.terms \ 4
    cm.n_terms = rows(cm.terms)
    cm.a_ci = b.links[ell].target_ci
    cm.cids = b.links[ell].target_cell
    if (ipwdr) {
        cm.a_ci = cm.a_ci \ b.links[ell].ctrl_ci
        cm.cids = cm.cids \ b.links[ell].ctrl_cell
    }
    cm.gids = J(0, 1, 0)
    cm.goff = J(0, 1, 0)
    cm.cell_gpos = J(rows(cm.cids), 1, 0)
    at = 0
    for (cp=1; cp<=rows(cm.cids); cp++) {
        gid = m.cells[cm.cids[cp]].group_id
        if (m.groups[gid].p < 1) continue
        gp = _cdid_impb_pos(cm.gids, gid)
        if (gp == 0) {
            cm.gids = cm.gids \ gid
            cm.goff = cm.goff \ at
            at = at + m.groups[gid].p
            gp = rows(cm.gids)
        }
        cm.cell_gpos[cp] = gp
    }
    cm.tid = (ipwdr ? b.links[ell].trt_id : 0)
    cm.p_trt = (ipwdr ? b.q : 0)
    cm.bid = 0
    cm.sep_or = 0
    cm.p_or = 0
    if (cm.use_or & aug) {
        cm.bid = (*pa).links[ell].or_block
        cm.sep_or = ((*pa).outcome_pooling == "separate")
        cm.p_or = (cm.sep_or ? (*pa).blocks[cm.bid].or_p : (*pa).or_p)
    }
    else if (cm.use_or) {
        cm.bid = b.links[ell].or_block
        cm.sep_or = (b.outcome_pooling == "separate")
        cm.p_or = (cm.sep_or ? b.blocks[cm.bid].or_p : b.or_p)
    }
    cm.th = J(0, 1, 0)
    if (aug) {
        ids = J(0, 1, 0)
        for (cp=1; cp<=rows(cm.cids); cp++) ids = ids \ (*pa).cells[cm.cids[cp]].theta_ids
        if (rows(ids)) cm.th = uniqrows(ids)
    }
    nth = rows(cm.th)
    cm.off_ia = rows(cm.a_ci)
    cm.off_ib = cm.off_ia + (cm.icpt ? rows(cm.cids) : 0)
    cm.off_ps = cm.off_ib + at
    cm.off_oa = cm.off_ps + cm.p_trt
    cm.off_ob = cm.off_oa + (cm.use_or & cm.icpt)
    cm.off_th = cm.off_ob + cm.p_or
    Dl = cm.off_th + nth
    cm.Dl = Dl
    coords = J(Dl, 1, "")
    for (j=1; j<=rows(cm.a_ci); j++) coords[j] = "a:" + _cdid_cpti_istr((rl - 1) * ip.G + cm.a_ci[j])
    if (cm.icpt) {
        for (cp=1; cp<=rows(cm.cids); cp++) coords[cm.off_ia + cp] = "ia:" + _cdid_cpti_istr(cm.cids[cp])
    }
    for (gp=1; gp<=rows(cm.gids); gp++) {
        for (j=1; j<=m.groups[cm.gids[gp]].p; j++) {
            coords[cm.off_ib + cm.goff[gp] + j] = "ib:" + _cdid_cpti_istr(cm.gids[gp]) + ":" + _cdid_cpti_istr(j)
        }
    }
    for (j=1; j<=cm.p_trt; j++) coords[cm.off_ps + j] = "ps:" + _cdid_cpti_istr(j)
    if (cm.use_or & cm.icpt) coords[cm.off_oa + 1] = "oa:" + _cdid_cpti_istr(cm.bid)
    for (j=1; j<=cm.p_or; j++) {
        if (cm.sep_or) coords[cm.off_ob + j] = "ob:" + _cdid_cpti_istr(cm.bid) + ":" + _cdid_cpti_istr(j)
        else coords[cm.off_ob + j] = "ob:" + _cdid_cpti_istr(j)
    }
    for (j=1; j<=nth; j++) coords[cm.off_th + j] = "th:" + _cdid_cpti_istr(cm.th[j])
    if (aug) {
        d_plus = (*pa).links[ell].d_plus
        d_minus = (*pa).links[ell].d_minus
        cm.estimate = (*pa).links[ell].mu
    }
    else {
        d_plus = b.links[ell].d_plus
        d_minus = b.links[ell].d_minus
        cm.estimate = b.links[ell].mu
    }
    s_all = (1 \ -1 \ -1 \ 1)
    D_all = (d_plus \ d_plus \ d_minus \ d_minus)
    svec = s_all[cm.terms]
    Dvec = D_all[cm.terms]
    cm.spec = cdid_cpt_spec_create(cm.id, cm.sign, svec, Dvec, coords)
    cm.estimate_scalar = cm.sign * cm.estimate
    if (cm.estimate >= .)
        _cdid_cpti_fail("numeric", "Nonfinite backend link estimate for cell " + cm.id + ".", 504)
    /* Bounded alpha j adjoint of the cell's OR block (pure: j loads; augmented: MD loads). */
    p = b.p_raw
    if (cm.use_or & cm.icpt) {
        lamA = J(b.n_cells, 1, 0)
        nuA = J(b.n_cells, p, 0)
        if (aug) {
            for (k=1; k<=(*pa).blocks[cm.bid].n_ctrl; k++) {
                cl = (*pa).blocks[cm.bid].ctrl_cell[k]
                lamA[cl] = lamA[cl] + (*pa).blocks[cm.bid].jf_alpha_alpha[k]
                if (p > 0) nuA[cl, ] = nuA[cl, ] + (*pa).blocks[cm.bid].jf_alpha_beta[k, ]
            }
        }
        else {
            for (k=1; k<=b.blocks[cm.bid].n_ctrl; k++) {
                cl = b.blocks[cm.bid].ctrl_cell[k]
                lamA[cl] = lamA[cl] + b.blocks[cm.bid].jf_alpha_alpha[k]
                if (p > 0) nuA[cl, ] = nuA[cl, ] + b.blocks[cm.bid].jf_alpha_beta[k, ]
            }
        }
        cm.adj_alpha = cdid_imputation_adjoint_prepare(m, lamA, nuA)
    }
    /* Conservative numeric envelope of one piece: lam/pim (2 R Dl), xi/P/scratch (4 Dl),
       nu/rho, design/score reads and theta contraction; fixed loads and factor scratch. */
    R = cm.n_terms
    per_row = 2 * R * Dl + 4 * Dl + 2 * R + 4 * (cm.p_trt + b.p_raw + cm.p_or + 8) + 2 * nth +
        3 * (nk + 1) + 32
    fixed = Dl * R + 2 * Dl + (b.p_raw + 2)^2 + (cm.p_or + 2)^2 + (cm.p_trt + 2)^2 +
        b.n_cells + b.n_groups * (b.p_raw + 1)
    if (aug) fixed = fixed + (*pa).n_q_params * (nth + 1 + cm.p_or) + (*pa).G * (nth + 1 + cm.p_or)
    Bmax = floor((ip.workspace_limit - ip.retained_bytes - 8 * fixed - 4096) / (8 * per_row))
    if (Bmax < 1)
        _cdid_cpti_fail("bound", "No admissible conditional piece row batch within the explicit workspace limit.", 908)
    cm.chunk_rows = min((ip.n_F, ip.chunk_rows, Bmax))
    cm.n_chunks = ceil(ip.n_F / cm.chunk_rows)
    nblk = b.n_blocks
    if (aug) nblk = (*pa).n_blocks
    if (cm.sep_or) nblk = 1
    cm.work_row = 16 * (2 * R * Dl + Dl * (cm.p_trt + b.p_raw + 4) +
        (nblk + 1) * (cm.p_or + 2) * (nk + 2) + rows(cm.cids) * (b.p_raw + 2)^2 + nth * (cm.p_or + 2))
}

/* Shared preparation of all four public entries. */
struct cdid_cpti scalar _cdid_cpti_prepare(string scalar kind, string scalar layout,
    pointer(struct cdid_imp_backend scalar) scalar pib,
    pointer(struct cdid_aug_backend scalar) scalar pab,
    struct cdid_public_request scalar req, real colvector orig_link, real colvector aux_link,
    real scalar workcap)
{
    struct cdid_cpti scalar ip
    struct cdid_cpti_cell_meta scalar cm
    pointer(struct cdid_imp_backend scalar) scalar pb
    pointer(struct cdid_imputation_model scalar) scalar pm
    real scalar fam, i, r, q, n_cells, eff, meta, total_work, bid
    real colvector sel_rows
    if (workcap < . & (workcap < 1 | workcap != floor(workcap)))
        _cdid_cpti_fail("bound", "workcap must be . or a positive integer.", 498)
    ip.kind = kind
    ip.layout = layout
    ip.ib = NULL
    ip.ab = NULL
    if (kind == "augmented") {
        if (pab == NULL)
            _cdid_cpti_fail("owner", "An accepted augmented backend pointer is required.", 498)
        _cdid_aug_check(*pab)
        if (((*pab).method != "ipw" & (*pab).method != "reg" & (*pab).method != "dr") |
            (*pab).selection_est_method != "dr" | (*pab).selection != "conditional_mtar")
            _cdid_cpti_fail("route", "The augmented conditional process requires est_method ipw, reg or dr with selection_est_method dr under conditional_mtar.", 498)
        if ((*pab).has_or) {
            if ((*pab).outcome_pooling == "separate") {
                for (bid=1; bid<=(*pab).n_blocks; bid++) {
                    if ((*pab).blocks[bid].or_p > 0 & !(*pab).blocks[bid].or_adj_ready)
                        _cdid_cpti_fail("owner", "Augmented block-local OR j adjoint is not prepared.", 498)
                }
            }
            else if ((*pab).or_p > 0 & !(*pab).or_adj_ready)
                _cdid_cpti_fail("owner", "Augmented common-slope j adjoint is not prepared.", 498)
        }
        ip.ab = pab
        pb = (*pab).ipw
        ip.method = (*pab).method
        ip.chunk_rows = (*pab).chunk_rows
        ip.retained_bytes = (*pab).retained_bytes
        ip.workspace_limit = (*pab).budget.max_workspace_bytes
        eff = ((*pab).budget.enforce_work_ops ? (*pab).budget.max_work_ops : .)
    }
    else if (kind == "imputation") {
        if (pib == NULL)
            _cdid_cpti_fail("owner", "An accepted pure-imputation backend pointer is required.", 498)
        _cdid_impb_need_model(*pib)
        pm = (*pib).model
        _cdid_impb_check(*pib, *pm)
        if ((*pib).method != "ipw" & (*pib).method != "reg" & (*pib).method != "dr")
            _cdid_cpti_fail("route", "Pure imputation supports est_method ipw, reg or dr only.", 498)
        if ((*pib).has_or) {
            if ((*pib).outcome_pooling == "separate") {
                for (bid=1; bid<=(*pib).n_blocks; bid++) {
                    if ((*pib).blocks[bid].or_p > 0 & !(*pib).blocks[bid].or_adj_ready)
                        _cdid_cpti_fail("owner", "Pure-imputation block-local OR j adjoint is not prepared.", 498)
                }
            }
            else if ((*pib).or_p > 0 & !(*pib).or_adj_ready)
                _cdid_cpti_fail("owner", "Pure-imputation common-slope j adjoint is not prepared.", 498)
        }
        ip.ib = pib
        pb = pib
        ip.method = (*pib).method
        ip.chunk_rows = (*pib).chunk_rows
        ip.retained_bytes = (*pib).retained_bytes
        ip.workspace_limit = (*pib).budget.max_workspace_bytes
        eff = ((*pib).budget.enforce_work_ops ? (*pib).budget.max_work_ops : .)
    }
    else _cdid_cpti_fail("route", "Unknown conditional imputation kind.", 498)
    _cdid_impb_need_model(*pb)
    pm = (*pb).model
    if ((layout == "dated") != ((*pb).is_dated == 1))
        _cdid_cpti_fail("route", "Use the " + ((*pb).is_dated ? "dated" : "static") + " preparation entry for this backend.", 498)
    ip.work_limit = (workcap < . ? (eff < . ? min((workcap, eff)) : workcap) : eff)
    ip.work_ops = 0
    ip.n_F = (*pb).N
    ip.n_A = (*pb).n_A
    ip.G = (*pb).G
    ip.n_recipes = (*pb).n_recipes
    ip.n_links = (*pb).n_links
    if (req.base_period != "varying" & req.base_period != "universal")
        _cdid_cpti_fail("request", "Invalid request base_period.", 498)
    if (rows(req.calendar) != rows((*pm).pop.calendar) | (*pm).pop.N != ip.n_F)
        _cdid_cpti_fail("identity", "Backend population and request calendar disagree.", 498)
    if (!all(req.calendar :== (*pm).pop.calendar))
        _cdid_cpti_fail("identity", "Backend population and request calendar disagree.", 498)
    if (cols(req.selection) != 9)
        _cdid_cpti_fail("request", "Request selection table must have 9 columns.", 498)
    ip.base_period = req.base_period
    _cdid_cpti_union(ip, *pb, *pm, req, orig_link, aux_link)

    fam = (req.base_period == "varying" ? 1 : 2)
    sel_rows = _cdid_impb_ix((req.selection[, 1] :== fam) :& (req.selection[, 8] :!= 1))
    n_cells = rows(sel_rows)
    ip.n_cells = n_cells
    ip.n_unavailable = 0
    ip.cells = cdid_cpti_cell_meta(n_cells, 1)
    meta = 64 * (n_cells + 1)
    if (ip.retained_bytes + 8 * meta > ip.workspace_limit)
        _cdid_cpti_fail("bound", "Conditional cell metadata exceeds the explicit workspace limit.", 908)
    if (n_cells == 0) {
        ip.status = "unavailable"
        ip.reason = "No declared pre-treatment cell exists in the base_period = " + char(34) + req.base_period + char(34) + " family (reference zeros excluded): there is no conditional pretrend restriction to test."
        ip.fingerprint = _cdid_cpti_fp(ip)
        return(ip)
    }
    total_work = 0
    for (i=1; i<=n_cells; i++) {
        r = sel_rows[i]
        cm = cdid_cpti_cell_meta()
        cm.g = req.selection[r, 2]
        cm.t = req.selection[r, 3]
        cm.from = req.selection[r, 5]
        cm.q = req.selection[r, 6]
        cm.sign = req.selection[r, 7]
        cm.id = "pre(g=" + _cdid_pr_num(cm.g) + ",t=" + _cdid_pr_num(cm.t) + ")"
        cm.s_date = .
        cm.t2_date = .
        cm.link = 0
        cm.recipe = 0
        cm.supported = 0
        cm.estimate = .
        cm.estimate_scalar = .
        cm.n_chunks = 0
        cm.chunk_rows = 0
        cm.route = ip.layout + " " + ip.kind + " " + ip.method
        if (req.selection[r, 8] != 0) {
            cm.status = "unavailable"
            cm.reason = req.selection_reason[r]
            ip.n_unavailable = ip.n_unavailable + 1
            ip.cells[i] = cm
            continue
        }
        q = cm.q
        if (q < 1 | q > ip.n_aux | q != floor(q) | (cm.sign != 1 & cm.sign != -1))
            _cdid_cpti_fail("request", "Available declared cell lacks an auxiliary row or a +/-1 sign.", 498)
        cm.s_date = req.aux[q, 2]
        cm.t2_date = req.aux[q, 3]
        cm.link = ip.aux_link[q]
        cm.supported = 1
        cm.status = "available"
        cm.reason = ""
        _cdid_cpti_layout(ip, *pb, *pm, cm)
        total_work = total_work + ip.n_F * cm.work_row
        ip.cells[i] = cm
    }
    ip.work_ops = ip.work_ops + total_work
    if (ip.work_limit < . & ip.work_ops > ip.work_limit)
        _cdid_cpti_fail("bound", "Counted conditional piece work exceeds the explicit work cap.", 908)
    if (ip.n_unavailable > 0) {
        ip.status = "unavailable"
        ip.reason = "The complete requested test is unavailable: " + _cdid_cpti_istr(ip.n_unavailable) + " declared cell(s) are unavailable. No available-only subset is tested."
    }
    else {
        ip.status = "available"
        ip.reason = ""
    }
    ip.fingerprint = _cdid_cpti_fp(ip)
    return(ip)
}

/* Public private entries. The static ones derive the (unique) union map; the dated ones
   accept the caller's recipe-aware maps (J(0,1,.) derives when unique). */
struct cdid_cpti scalar cdid_cpti_prepare_static(
    pointer(struct cdid_imp_backend scalar) scalar ib,
    struct cdid_public_request scalar req, real scalar workcap)
{
    return(_cdid_cpti_prepare("imputation", "static", ib, NULL, req, J(0, 1, .), J(0, 1, .), workcap))
}

struct cdid_cpti scalar cdid_cpti_prepare_dated(
    pointer(struct cdid_imp_backend scalar) scalar ib,
    struct cdid_public_request scalar req, real colvector orig_link, real colvector aux_link,
    real scalar workcap)
{
    return(_cdid_cpti_prepare("imputation", "dated", ib, NULL, req, orig_link, aux_link, workcap))
}

struct cdid_cpti scalar cdid_cpti_prepare_aug_static(
    pointer(struct cdid_aug_backend scalar) scalar ab,
    struct cdid_public_request scalar req, real scalar workcap)
{
    return(_cdid_cpti_prepare("augmented", "static", NULL, ab, req, J(0, 1, .), J(0, 1, .), workcap))
}

struct cdid_cpti scalar cdid_cpti_prepare_aug_dated(
    pointer(struct cdid_aug_backend scalar) scalar ab,
    struct cdid_public_request scalar req, real colvector orig_link, real colvector aux_link,
    real scalar workcap)
{
    return(_cdid_cpti_prepare("augmented", "dated", NULL, ab, req, orig_link, aux_link, workcap))
}

struct cdid_cpt_spec scalar cdid_cpti_cell_spec(struct cdid_cpti scalar ip, real scalar ci)
{
    if (ci < 1 | ci > ip.n_cells | ci != floor(ci))
        _cdid_cpti_fail("cell", "Invalid declared cell index.", 498)
    if (!ip.cells[ci].supported)
        _cdid_cpti_fail("cell", "Requested cell is unavailable: " + ip.cells[ci].reason, 459)
    return(ip.cells[ci].spec)
}

real scalar cdid_cpti_cell_chunks(struct cdid_cpti scalar ip, real scalar ci)
{
    if (ci < 1 | ci > ip.n_cells | ci != floor(ci))
        _cdid_cpti_fail("cell", "Invalid declared cell index.", 498)
    if (!ip.cells[ci].supported)
        _cdid_cpti_fail("cell", "Requested cell is unavailable: " + ip.cells[ci].reason, 459)
    return(ip.cells[ci].n_chunks)
}

/* Repeat-rate IF xi_{a_{r,k}} = W_r I_k (A - a_{r,k}) / mass_{r,k} on the chunk (frame
   scale, off-A rows included), exactly as the backend replay kernel. */
real colvector _cdid_cpti_xa(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real colvector fr, real colvector coh, real colvector Av, real scalar r, real scalar k)
{
    real colvector out, ix, wr
    out = J(rows(fr), 1, 0)
    if (b.is_dated) {
        if (!(*(b.rdr)).used_cohorts[r, k]) return(out)
    }
    if (!(b.mass_F_mat[r, k] > 0)) return(out)
    ix = _cdid_impb_ix(coh :== k)
    if (rows(ix) == 0) return(out)
    wr = _cdid_impb_w(b, m, r, fr[ix])
    out[ix] = wr :* (Av[ix] :- b.a_k_mat[r, k]) :/ b.mass_F_mat[r, k]
    return(out)
}

/* Treatment gamma IF of key tid: direct v (D - p) X and repeat-rate a_forcing channels
   (source trt_if[[tid]]$gamma = direct + a), raw Gram factor with N exactly once. */
real matrix _cdid_cpti_gamma(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real colvector fr, real colvector coh, real colvector Av, real scalar tid)
{
    real scalar B, q, rt, k, nc
    real colvector mask, ix, cix, yv, wr, v
    real matrix sd, sa, xa, X, vals
    B = rows(fr)
    q = b.q
    rt = b.trts[tid].recipe
    nc = b.trts[tid].n_coh
    sd = J(B, q, 0)
    mask = J(B, 1, 0)
    for (k=1; k<=nc; k++) mask = mask :| (coh :== b.trts[tid].all_ci[k])
    ix = _cdid_impb_ix(mask)
    if (rows(ix) > 0) {
        cix = coh[ix]
        yv = (cix :== b.trts[tid].target_ci)
        wr = _cdid_impb_w(b, m, rt, fr[ix])
        v = wr :* (b.a_k_mat[rt, cix]')
        X = _cdid_impb_x(b, rt, fr[ix])
        vals = _cdid_ipw_logit_values(X * b.trts[tid].gamma)
        sd[ix, ] = (v :* (yv - vals[, 1])) :* X
    }
    xa = J(B, nc, 0)
    for (k=1; k<=nc; k++) xa[, k] = _cdid_cpti_xa(b, m, fr, coh, Av, rt, b.trts[tid].all_ci[k])
    sa = xa * b.trts[tid].a_forcing'
    return(_cdid_nuisance_gram(b.trts[tid].R, b.trts[tid].pivot, b.N :* sd')' +
        _cdid_nuisance_gram(b.trts[tid].R, b.trts[tid].pivot, b.N :* sa')')
}

/* Stage 1 coordinate IF rows (frame scale) through the static or dated APIs. */
real colvector _cdid_cpti_alpha_if(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar cid, real colvector fr)
{
    if (b.is_dated) return(cdid_imp_alpha_if_d(m, *(b.rdr), cid, fr, "frame"))
    return(cdid_imputation_alpha_if_chunk(m, cid, fr, "frame"))
}

real matrix _cdid_cpti_beta_if(struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar gid, real colvector fr)
{
    if (b.is_dated) return(cdid_imp_beta_if_d(m, *(b.rdr), gid, fr, "frame"))
    return(cdid_imputation_beta_if_chunk(m, gid, fr, "frame"))
}

/* Pure-imputation shared OR IF of block bid (alpha, when intercepts exist) and of the
   common slope: direct + a + j channels, each solved separately and summed (source
   or_if$alpha / or_if$beta; j by the bounded adjoints). */
void _cdid_cpti_or_imp(struct cdid_imp_backend scalar b, struct cdid_imputation_model scalar m,
    real colvector fr, real colvector coh, real colvector Av, real scalar bid,
    struct cdid_imputation_adjoint scalar adja, real colvector oa, real matrix ob)
{
    real scalar B, po, t, bb, rb, L, alpha_b, k, cl, ck, isb, mass, sep
    real colvector bids, ix, v, jh, mh, rv, xc, rad, raa, raj, obeta
    real rowvector okeep
    real matrix rbd, rba, rbj, dor, dt, xbd, xba, xbj, Jaj
    B = rows(fr)
    sep = (b.outcome_pooling == "separate")
    po = (sep ? b.blocks[bid].or_p : b.or_p)
    okeep = (sep ? b.blocks[bid].or_keep : b.or_keep)
    obeta = (sep ? b.blocks[bid].or_beta : b.or_beta)
    rbd = J(B, po, 0)
    rba = J(B, po, 0)
    rad = J(B, 1, 0)
    raa = J(B, 1, 0)
    bids = (sep ? J(1, 1, bid) : (po > 0 ? (1::b.n_blocks) : J(1, 1, bid)))
    for (t=1; t<=rows(bids); t++) {
        bb = bids[t]
        isb = (bb == bid & b.intercept)
        if (!isb & po == 0) continue
        rb = b.blocks[bb].recipe
        L = b.blocks[bb].L
        alpha_b = (b.intercept ? b.blocks[bb].alpha : 0)
        for (k=1; k<=b.blocks[bb].n_ctrl; k++) {
            cl = b.blocks[bb].ctrl_cell[k]
            ck = b.blocks[bb].ctrl_ci[k]
            ix = _cdid_impb_ix(coh :== ck)
            if (rows(ix) > 0) {
                v = _cdid_impb_w(b, m, rb, fr[ix]) :* b.a_k_mat[rb, ck]
                jh = _cdid_impb_pred_cell(b, m, cl, fr[ix])
                if (po > 0) {
                    dor = L :* _cdid_impb_z(b, m, rb, fr[ix], okeep)
                    mh = alpha_b :+ dor * obeta
                }
                else mh = J(rows(ix), 1, alpha_b)
                rv = v :* (jh - mh)
                if (isb) rad[ix] = rad[ix] + rv
                if (po > 0) {
                    dt = (b.intercept ? dor :- b.blocks[bb].mean_d : dor)
                    rbd[ix, ] = rbd[ix, ] + rv :* dt
                }
            }
            xc = _cdid_cpti_xa(b, m, fr, coh, Av, rb, ck)
            if (isb) raa = raa + b.blocks[bb].a_force_alpha[k] :* xc
            if (po > 0) rba = rba + xc * b.blocks[bb].a_force_beta[, k]'
        }
    }
    if (sep) {
        if (b.intercept) rad = rad :- b.blocks[bid].score_mean_a
        if (po > 0) rbd = rbd :- (b.blocks[bid].score_mean_b')
    }
    xbd = J(B, po, 0)
    xba = J(B, po, 0)
    xbj = J(B, po, 0)
    ob = J(B, po, 0)
    if (po > 0) {
        if (sep) {
            if (!b.blocks[bid].or_adj_ready)
                _cdid_cpti_fail("owner", "Block-local OR j adjoint is not prepared.", 498)
            rbj = _cdid_impb_adj_chunk(b, m, b.blocks[bid].or_adj, fr)
            xbd = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, b.N :* rbd')'
            xba = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, b.N :* rba')'
            xbj = _cdid_nuisance_gram(b.blocks[bid].or_R, b.blocks[bid].or_pivot, b.N :* rbj')'
        }
        else {
            if (!b.or_adj_ready) _cdid_cpti_fail("owner", "Common-slope j adjoint is not prepared.", 498)
            rbj = _cdid_impb_adj_chunk(b, m, b.or_adj, fr)
            xbd = _cdid_impb_or_solve(b, rbd)
            xba = _cdid_impb_or_solve(b, rba)
            xbj = _cdid_impb_or_solve(b, rbj)
        }
        ob = xbd + xba + xbj
    }
    oa = J(B, 1, 0)
    if (b.intercept) {
        Jaj = _cdid_impb_adj_chunk(b, m, adja, fr)
        raj = Jaj[, 1]
        mass = b.blocks[bid].mass
        if (po > 0) {
            oa = ((rad :/ mass) - xbd * b.blocks[bid].mean_d') +
                ((raa :/ mass) - xba * b.blocks[bid].mean_d') +
                ((raj :/ mass) - xbj * b.blocks[bid].mean_d')
        }
        else oa = (rad :/ mass) + (raa :/ mass) + (raj :/ mass)
    }
    if (missing(oa) | missing(ob)) _cdid_cpti_fail("numeric", "Nonfinite shared OR coordinate IF.", 504)
}

/* Theta-only selection contraction of one augmented piece in ONE moment record:
   unit loads at the cell's theta IDs (theta IF columns), then the block q_alpha load
   (when intercepts exist) and the OR q_beta loads (p_OR columns). G_ak stays zero. */
real matrix _cdid_cpti_sel(struct cdid_aug_backend scalar a, struct cdid_sel_reader scalar qr,
    real colvector th, real scalar bid, real scalar want_alpha, real colvector fr)
{
    struct cdid_sel_moments scalar sm
    real scalar nth, po, nmom, j, c0
    real colvector ids
    real matrix S
    nth = rows(th)
    po = 0
    if (bid > 0 & a.has_or) {
        po = (a.outcome_pooling == "separate" ? a.blocks[bid].or_p : a.or_p)
    }
    else want_alpha = 0
    nmom = nth + want_alpha + po
    if (nmom == 0) return(J(rows(fr), 0, 0))
    if (a.is_dated) {
        sm.n_moments = nmom
        sm.total_params = a.n_q_params
        sm.n_cohorts = a.G
        sm.G_theta = J(a.n_q_params, nmom, 0)
        sm.G_ak = J(a.G, nmom, 0)
    }
    else sm = cdid_sel_moments_init(a.qfit, nmom)
    for (j=1; j<=nth; j++) sm.G_theta[th[j], j] = 1
    if (want_alpha & bid > 0) {
        ids = a.blocks[bid].q_alpha_ids
        c0 = nth + 1
        if (rows(ids) > 0) sm.G_theta[ids, c0] = sm.G_theta[ids, c0] + a.blocks[bid].q_alpha
    }
    if (po > 0 & bid > 0) {
        if (a.outcome_pooling == "separate") {
            ids = a.blocks[bid].q_beta_ids
            if (rows(ids) > 0) {
                for (j=1; j<=po; j++) {
                    c0 = nth + want_alpha + j
                    sm.G_theta[ids, c0] = sm.G_theta[ids, c0] + a.blocks[bid].q_beta[j, ]'
                }
            }
        }
        else if (rows(a.q_beta_ids) > 0) {
            for (j=1; j<=po; j++) {
                c0 = nth + want_alpha + j
                sm.G_theta[a.q_beta_ids, c0] = sm.G_theta[a.q_beta_ids, c0] + a.q_beta[j, ]'
            }
        }
    }
    if (missing(sm.G_theta) | any(sm.G_ak :!= 0))
        _cdid_cpti_fail("numeric", "Invalid theta-only observation loading record.", 504)
    if (a.is_dated) S = _cdid_aug_q_loading(a, sm, fr)
    else S = cdid_sel_eval_loading(a.qfit, qr, sm, fr)
    if (rows(S) != rows(fr) | cols(S) != nmom)
        _cdid_cpti_fail("integrity", "Observation contraction has an unexpected shape.", 498)
    return(S)
}

/* Augmented Z-response OR IF of block bid and its slope: direct (v (Z - m) on
   control rows), rate, j (MD-loaded adjoints) and q (theta-only) channels, each solved
   separately and summed. sqa/sqb are the q rhs columns from _cdid_cpti_sel. */
void _cdid_cpti_or_aug(struct cdid_aug_backend scalar a, struct cdid_imp_backend scalar ib,
    struct cdid_imputation_model scalar m, struct cdid_sel_reader scalar qr,
    real colvector fr, real colvector coh, real colvector Av, real scalar bid,
    struct cdid_imputation_adjoint scalar adja, real colvector sqa, real matrix sqb,
    real colvector oa, real matrix ob)
{
    real scalar B, po, t, bb, rb, L, alpha_b, k, cl, ck, isb, mass, sep
    real colvector bids, ix, v, mh, rv, center_weight, xc, rad, raa, raj, jh, zr, obs, qv, piv, ev, obeta
    real rowvector okeep
    real matrix rbd, rba, rbj, dor, dt, direct_beta, xbd, xba, xbj, xbq, Tm, Jaj
    B = rows(fr)
    sep = (a.outcome_pooling == "separate")
    po = (sep ? a.blocks[bid].or_p : a.or_p)
    okeep = (sep ? a.blocks[bid].or_keep : a.or_keep)
    obeta = (sep ? a.blocks[bid].or_beta : a.or_beta)
    rbd = J(B, po, 0)
    rba = J(B, po, 0)
    rad = J(B, 1, 0)
    raa = J(B, 1, 0)
    bids = (sep ? J(1, 1, bid) : (po > 0 ? (1::a.n_blocks) : J(1, 1, bid)))
    for (t=1; t<=rows(bids); t++) {
        bb = bids[t]
        isb = (bb == bid & a.intercept)
        if (!isb & po == 0) continue
        rb = a.blocks[bb].recipe
        L = a.blocks[bb].L
        alpha_b = (a.intercept ? a.blocks[bb].alpha : 0)
        for (k=1; k<=a.blocks[bb].n_ctrl; k++) {
            cl = a.blocks[bb].ctrl_cell[k]
            ck = a.blocks[bb].ctrl_ci[k]
            ix = _cdid_impb_ix(coh :== ck)
            if (rows(ix) > 0) {
                v = _cdid_impb_w(ib, m, rb, fr[ix]) :* ib.a_k_mat[rb, ck]
                _cdid_aug_cell_rows(a, ib, m, qr, cl, fr[ix], 0, jh=., zr=., obs=., qv=., piv=., ev=., Tm=.)
                if (po > 0) {
                    dor = L :* _cdid_impb_z(ib, m, rb, fr[ix], okeep)
                    mh = alpha_b :+ dor * obeta
                }
                else mh = J(rows(ix), 1, alpha_b)
                rv = v :* (zr - mh)
                /* Center the direct scores on this OR block's weighted support. */
                if (sep) center_weight = v :/ a.blocks[bb].mass
                if (isb) {
                    rad[ix] = rad[ix] + (sep ? rv :- center_weight :* a.blocks[bb].score_mean_a : rv)
                }
                if (po > 0) {
                    dt = (a.intercept ? dor :- a.blocks[bb].mean_d : dor)
                    direct_beta = rv :* dt
                    if (sep) direct_beta = direct_beta - center_weight * a.blocks[bb].score_mean_b'
                    rbd[ix, ] = rbd[ix, ] + direct_beta
                }
            }
            xc = _cdid_cpti_xa(ib, m, fr, coh, Av, rb, ck)
            if (isb) raa = raa + a.blocks[bb].a_force_alpha[k] :* xc
            if (po > 0) rba = rba + xc * a.blocks[bb].a_force_beta[, k]'
        }
    }
    xbd = J(B, po, 0)
    xba = J(B, po, 0)
    xbj = J(B, po, 0)
    xbq = J(B, po, 0)
    ob = J(B, po, 0)
    if (po > 0) {
        if (sep) {
            if (!a.blocks[bid].or_adj_ready)
                _cdid_cpti_fail("owner", "Augmented block-local OR j adjoint is not prepared.", 498)
            rbj = _cdid_impb_adj_chunk(ib, m, a.blocks[bid].or_adj, fr)
            xbd = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, a.N :* rbd')'
            xba = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, a.N :* rba')'
            xbj = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, a.N :* rbj')'
            xbq = _cdid_nuisance_gram(a.blocks[bid].or_R, a.blocks[bid].or_pivot, a.N :* sqb')'
        }
        else {
            if (!a.or_adj_ready) _cdid_cpti_fail("owner", "Augmented common-slope j adjoint is not prepared.", 498)
            rbj = _cdid_impb_adj_chunk(ib, m, a.or_adj, fr)
            xbd = _cdid_aug_or_solve(a, rbd)
            xba = _cdid_aug_or_solve(a, rba)
            xbj = _cdid_aug_or_solve(a, rbj)
            xbq = _cdid_aug_or_solve(a, sqb)
        }
        ob = xbd + xba + xbj + xbq
    }
    oa = J(B, 1, 0)
    if (a.intercept) {
        Jaj = _cdid_impb_adj_chunk(ib, m, adja, fr)
        raj = Jaj[, 1]
        mass = a.blocks[bid].mass
        if (po > 0) {
            oa = ((rad :/ mass) - xbd * a.blocks[bid].mean_d') +
                ((raa :/ mass) - xba * a.blocks[bid].mean_d') +
                ((raj :/ mass) - xbj * a.blocks[bid].mean_d') +
                ((sqa :/ mass) - xbq * a.blocks[bid].mean_d')
        }
        else oa = (rad :/ mass) + (raa :/ mass) + (raj :/ mass) + (sqa :/ mass)
    }
    if (missing(oa) | missing(ob)) _cdid_cpti_fail("numeric", "Nonfinite augmented OR coordinate IF.", 504)
}

/* Source cell_resp: pure imputation response j with unit derivative factor; augmented
   response Z = j + S (Delta - j) / max(q, 1e-6) with factor 1 - S/max(q, 1e-6) and the
   explicit theta derivative e pi T_logit on observed pairs. d = L z_kept of the cell group. */
void _cdid_cpti_resp(struct cdid_cpti scalar ip, struct cdid_imp_backend scalar b,
    struct cdid_imputation_model scalar m, real scalar cid, real colvector rch,
    real colvector rho, real colvector fac, real matrix d, real matrix pt)
{
    pointer(struct cdid_aug_backend scalar) scalar pa
    pointer(struct cdid_sel_reader scalar) scalar pq
    real scalar n, gid, r, dq
    real colvector jh, zr, obs, qv, piv, ev
    real matrix Tm
    n = rows(rch)
    gid = m.cells[cid].group_id
    r = m.cells[cid].recipe
    if (m.groups[gid].p > 0) d = m.cells[cid].L :* _cdid_impb_z(b, m, r, rch, m.groups[gid].kept)
    else d = J(n, 0, 0)
    fac = J(n, 1, 1)
    if (ip.kind != "augmented") {
        rho = _cdid_impb_pred_cell(b, m, cid, rch)
        pt = J(n, 0, 0)
        return
    }
    pa = ip.ab
    pq = _cdid_cpti_aug_qr(*pa)
    dq = rows((*pa).cells[cid].theta_ids)
    _cdid_aug_cell_rows(*pa, b, m, *pq, cid, rch, 1, jh=., zr=., obs=., qv=., piv=., ev=., Tm=.)
    rho = zr
    pt = J(n, dq, 0)
    if (rows(obs) > 0) {
        fac[obs] = 1 :- piv
        if (dq > 0) pt[obs, ] = (ev :* piv) :* Tm
    }
}

/* Source put_resp into the term block starting at column tb + 1 for cell position cp. */
void _cdid_cpti_put(struct cdid_cpti scalar ip, struct cdid_cpti_cell_meta scalar cm,
    real matrix pim, real scalar tb, real colvector ii, real scalar cp,
    real colvector fac, real matrix d, real matrix pt)
{
    real scalar gp
    real colvector thpos
    if (cm.icpt) pim[ii, tb + cm.off_ia + cp] = fac
    gp = cm.cell_gpos[cp]
    if (gp > 0 & cols(d) > 0) pim[ii, tb + cm.off_ib + cm.goff[gp] :+ (1..cols(d))] = fac :* d
    if (ip.kind == "augmented" & cols(pt) > 0) {
        thpos = _cdid_aug_find(cm.th, (*(ip.ab)).cells[cm.cids[cp]].theta_ids)
        pim[ii, tb + cm.off_th :+ thpos'] = pt
    }
}

/* Source rows_piece for one available cell on unique original-F rows fr. */
struct cdid_cpt_piece scalar _cdid_cpti_eval(struct cdid_cpti scalar ip, real scalar ci,
    real colvector fr, | transmorphic scalar cache)
{
    struct cdid_cpti_cell_meta scalar cm
    pointer(struct cdid_imp_backend scalar) scalar pb
    pointer(struct cdid_imputation_model scalar) scalar pm
    pointer(struct cdid_aug_backend scalar) scalar pa
    pointer(struct cdid_sel_reader scalar) scalar pq
    real scalar aug, B, R, Dl, ell, rl, cg, g, L, rT, rBg, rC, rBo, tb, ag, akc, alpha_b
    real scalar po, q, nk, k, ck, cl, j, gp, nth, n, use_cache
    real colvector coh, Av, ig, ik, nug, nuk, resp, fac, sD, oa, odds, act, sqa, mo
    real matrix nu, rho, lam, pim, xi, d, pt, z, X, vals, ob, S, sqb, xb
    real rowvector okeep
    real colvector obeta
    string scalar ckey, ckey2
    use_cache = (args() >= 4)
    cm = ip.cells[ci]
    aug = (ip.kind == "augmented")
    if (aug) {
        pa = ip.ab
        pb = (*pa).ipw
        pq = _cdid_cpti_aug_qr(*pa)
    }
    else pb = ip.ib
    pm = (*pb).model
    B = rows(fr)
    R = cm.n_terms
    Dl = cm.Dl
    ell = cm.link
    rl = (*pb).links[ell].recipe
    cg = (*pb).links[ell].target_ci
    g = (*pb).links[ell].target_cell
    L = (*pb).links[ell].L
    nk = (*pb).links[ell].n_ctrl
    po = cm.p_or
    q = cm.p_trt
    nth = rows(cm.th)
    rT = 1
    rBg = _cdid_impb_pos(cm.terms, 2)
    rC = _cdid_impb_pos(cm.terms, 3)
    rBo = _cdid_impb_pos(cm.terms, 4)
    okeep = J(1, 0, 0)
    obeta = J(0, 1, 0)
    alpha_b = 0
    if (cm.use_or) {
        if (aug) {
            okeep = (cm.sep_or ? (*pa).blocks[cm.bid].or_keep : (*pa).or_keep)
            obeta = (cm.sep_or ? (*pa).blocks[cm.bid].or_beta : (*pa).or_beta)
            alpha_b = ((*pa).intercept ? (*pa).blocks[cm.bid].alpha : 0)
        }
        else {
            okeep = (cm.sep_or ? (*pb).blocks[cm.bid].or_keep : (*pb).or_keep)
            obeta = (cm.sep_or ? (*pb).blocks[cm.bid].or_beta : (*pb).or_beta)
            alpha_b = ((*pb).intercept ? (*pb).blocks[cm.bid].alpha : 0)
        }
    }
    coh = (*pm).unit_cohort_index[fr]
    Av = (*pm).pop.A[fr]
    nu = J(B, R, 0)
    rho = J(B, R, 0)
    lam = J(B, R * Dl, 0)
    pim = J(B, R * Dl, 0)
    xi = J(B, Dl, 0)

    /* Treated target rows: T (and Bg). */
    ig = _cdid_impb_ix(coh :== cg)
    n = rows(ig)
    if (n > 0) {
        ag = (*pb).a_k_mat[rl, cg]
        nug = _cdid_impb_w(*pb, *pm, rl, fr[ig]) :* ag
        _cdid_cpti_resp(ip, *pb, *pm, g, fr[ig], resp=., fac=., d=., pt=.)
        tb = (rT - 1) * Dl
        nu[ig, rT] = nug
        rho[ig, rT] = resp
        lam[ig, tb + 1] = J(n, 1, 1 / ag)
        _cdid_cpti_put(ip, cm, pim, tb, ig, 1, fac, d, pt)
        if (rBg > 0) {
            z = (po > 0 ? L :* _cdid_impb_z(*pb, *pm, rl, fr[ig], okeep) : J(n, 0, 0))
            mo = J(n, 1, (cm.icpt ? alpha_b : 0))
            if (po > 0) mo = mo + z * obeta
            tb = (rBg - 1) * Dl
            nu[ig, rBg] = nug
            rho[ig, rBg] = mo
            lam[ig, tb + 1] = J(n, 1, 1 / ag)
            if (cm.icpt) pim[ig, tb + cm.off_oa + 1] = J(n, 1, 1)
            if (po > 0) pim[ig, tb + cm.off_ob :+ (1..po)] = z
        }
    }

    /* Control cohorts: C (and Bo) with propensity odds. */
    if (rC > 0) {
        for (k=1; k<=nk; k++) {
            ck = (*pb).links[ell].ctrl_ci[k]
            cl = (*pb).links[ell].ctrl_cell[k]
            ik = _cdid_impb_ix(coh :== ck)
            n = rows(ik)
            if (n == 0) continue
            X = _cdid_impb_x(*pb, rl, fr[ik])
            vals = _cdid_ipw_logit_values(X * (*pb).trts[cm.tid].gamma)
            odds = vals[, 2]
            act = vals[, 3]
            akc = (*pb).a_k_mat[rl, ck]
            nuk = (_cdid_impb_w(*pb, *pm, rl, fr[ik]) :* akc) :* odds
            _cdid_cpti_resp(ip, *pb, *pm, cl, fr[ik], resp=., fac=., d=., pt=.)
            tb = (rC - 1) * Dl
            nu[ik, rC] = nuk
            rho[ik, rC] = resp
            lam[ik, tb + 1 + k] = J(n, 1, 1 / akc)
            lam[ik, tb + cm.off_ps :+ (1..q)] = act :* X
            _cdid_cpti_put(ip, cm, pim, tb, ik, 1 + k, fac, d, pt)
            if (rBo > 0) {
                z = (po > 0 ? L :* _cdid_impb_z(*pb, *pm, rl, fr[ik], okeep) : J(n, 0, 0))
                mo = J(n, 1, (cm.icpt ? alpha_b : 0))
                if (po > 0) mo = mo + z * obeta
                tb = (rBo - 1) * Dl
                nu[ik, rBo] = nuk
                rho[ik, rBo] = mo
                lam[ik, tb + 1 + k] = J(n, 1, 1 / akc)
                lam[ik, tb + cm.off_ps :+ (1..q)] = act :* X
                if (cm.icpt) pim[ik, tb + cm.off_oa + 1] = J(n, 1, 1)
                if (po > 0) pim[ik, tb + cm.off_ob :+ (1..po)] = z
            }
        }
    }

    /* Coordinate IF rows of exactly the coordinates this cell reads (frame scale). */
    for (j=1; j<=rows(cm.a_ci); j++) {
        if (use_cache) {
            ckey = "xa:" + strofreal(rl) + ":" + strofreal(cm.a_ci[j])
            if (asarray_contains(cache, ckey)) xi[, j] = asarray(cache, ckey)
            else {
                xi[, j] = _cdid_cpti_xa(*pb, *pm, fr, coh, Av, rl, cm.a_ci[j])
                asarray(cache, ckey, xi[, j])
            }
        }
        else xi[, j] = _cdid_cpti_xa(*pb, *pm, fr, coh, Av, rl, cm.a_ci[j])
    }
    if (cm.icpt) {
        for (j=1; j<=rows(cm.cids); j++) {
            if (use_cache) {
                ckey = "ia:" + strofreal(cm.cids[j])
                if (asarray_contains(cache, ckey)) xi[, cm.off_ia + j] = asarray(cache, ckey)
                else {
                    xi[, cm.off_ia + j] = _cdid_cpti_alpha_if(*pb, *pm, cm.cids[j], fr)
                    asarray(cache, ckey, xi[, cm.off_ia + j])
                }
            }
            else xi[, cm.off_ia + j] = _cdid_cpti_alpha_if(*pb, *pm, cm.cids[j], fr)
        }
    }
    for (gp=1; gp<=rows(cm.gids); gp++) {
        if (use_cache) {
            ckey = "ib:" + strofreal(cm.gids[gp])
            if (asarray_contains(cache, ckey)) xb = asarray(cache, ckey)
            else {
                xb = _cdid_cpti_beta_if(*pb, *pm, cm.gids[gp], fr)
                asarray(cache, ckey, xb)
            }
        }
        else xb = _cdid_cpti_beta_if(*pb, *pm, cm.gids[gp], fr)
        if (rows(xb) != B | cols(xb) != (*pm).groups[cm.gids[gp]].p)
            _cdid_cpti_fail("integrity", "First-step group slope IF width differs from its kept columns.", 498)
        xi[, cm.off_ib + cm.goff[gp] :+ (1..cols(xb))] = xb
    }
    if (q > 0) {
        if (use_cache) {
            ckey = "ps:" + strofreal(cm.tid)
            if (asarray_contains(cache, ckey)) xi[, cm.off_ps :+ (1..q)] = asarray(cache, ckey)
            else {
                xi[, cm.off_ps :+ (1..q)] = _cdid_cpti_gamma(*pb, *pm, fr, coh, Av, cm.tid)
                asarray(cache, ckey, xi[, cm.off_ps :+ (1..q)])
            }
        }
        else xi[, cm.off_ps :+ (1..q)] = _cdid_cpti_gamma(*pb, *pm, fr, coh, Av, cm.tid)
    }
    if (cm.use_or) {
        if (aug) {
            S = _cdid_cpti_sel(*pa, *pq, cm.th, cm.bid, cm.icpt, fr)
            sqa = J(B, 1, 0)
            sqb = J(B, 0, 0)
            if (cm.icpt) sqa = S[, nth + 1]
            if (po > 0) sqb = S[, nth + cm.icpt :+ (1..po)]
            _cdid_cpti_or_aug(*pa, *pb, *pm, *pq, fr, coh, Av, cm.bid, cm.adj_alpha, sqa, sqb, oa=., ob=.)
            if (nth > 0) xi[, cm.off_th :+ (1..nth)] = S[, 1..nth]
        }
        else if (use_cache) {
            ckey = "oa:" + strofreal(cm.bid)
            ckey2 = "ob:" + strofreal(cm.bid)
            if (asarray_contains(cache, ckey)) {
                oa = asarray(cache, ckey)
                ob = asarray(cache, ckey2)
            }
            else {
                _cdid_cpti_or_imp(*pb, *pm, fr, coh, Av, cm.bid, cm.adj_alpha, oa=., ob=.)
                asarray(cache, ckey, oa)
                asarray(cache, ckey2, ob)
            }
        }
        else _cdid_cpti_or_imp(*pb, *pm, fr, coh, Av, cm.bid, cm.adj_alpha, oa=., ob=.)
        if (cm.icpt) xi[, cm.off_oa + 1] = oa
        if (po > 0) xi[, cm.off_ob :+ (1..po)] = ob
    }
    else if (aug & nth > 0) {
        S = _cdid_cpti_sel(*pa, *pq, cm.th, 0, 0, fr)
        xi[, cm.off_th :+ (1..nth)] = S[, 1..nth]
    }
    sD = cm.spec.s :/ cm.spec.D
    return(_cdid_cpti_piece_make(fr, nu, rho, lam, pim, sD, xi, cm.spec))
}

/* One bounded piece on caller-supplied unique original-F rows (1..chunk_rows of them). */
struct cdid_cpt_piece scalar cdid_cpti_piece_rows(struct cdid_cpti scalar ip, real scalar ci,
    real colvector rowsF)
{
    real scalar n
    if (ci < 1 | ci > ip.n_cells | ci != floor(ci))
        _cdid_cpti_fail("replay", "Invalid declared cell index.", 498)
    if (!ip.cells[ci].supported)
        _cdid_cpti_fail("replay", "Requested cell is unavailable: " + ip.cells[ci].reason, 459)
    n = rows(rowsF)
    if (n < 1 | cols(rowsF) != 1 | missing(rowsF))
        _cdid_cpti_fail("replay", "rowsF must be 1..chunk_rows unique original-F positions in 1..N_F.", 498)
    if (any(rowsF :!= floor(rowsF)) | any(rowsF :< 1) | any(rowsF :> ip.n_F) |
        n > ip.cells[ci].chunk_rows)
        _cdid_cpti_fail("replay", "rowsF must be 1..chunk_rows unique original-F positions in 1..N_F.", 498)
    if (!(n <= 1 ? 1 : all(rowsF[(2..n)] :> rowsF[(1..(n - 1))]))) {
        if (rows(uniqrows(rowsF)) != n)
            _cdid_cpti_fail("replay", "rowsF must not repeat an original-F unit.", 498)
    }
    _cdid_cpti_check(ip)
    return(_cdid_cpti_eval(ip, ci, rowsF))
}

/* Chunk k of the full-F partition 1..N_F (source pieces(k) on original rows). */
struct cdid_cpt_piece scalar cdid_cpti_piece_chunk(struct cdid_cpti scalar ip, real scalar ci,
    real scalar chunk_k)
{
    real scalar lo, hi
    if (ci < 1 | ci > ip.n_cells | ci != floor(ci))
        _cdid_cpti_fail("replay", "Invalid declared cell index.", 498)
    if (!ip.cells[ci].supported)
        _cdid_cpti_fail("replay", "Requested cell is unavailable: " + ip.cells[ci].reason, 459)
    if (chunk_k < 1 | chunk_k > ip.cells[ci].n_chunks | chunk_k != floor(chunk_k))
        _cdid_cpti_fail("replay", "Invalid piece chunk index.", 498)
    lo = (chunk_k - 1) * ip.cells[ci].chunk_rows + 1
    hi = min((ip.n_F, lo + ip.cells[ci].chunk_rows - 1))
    return(cdid_cpti_piece_rows(ip, ci, (lo::hi)))
}

/* Stream every original-F row of one cell into the engine accumulator. */
struct cdid_cpt_cell scalar cdid_cpti_build_cell(struct cdid_cpti scalar ip, real scalar ci,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar st, en, bsz
    spec = cdid_cpti_cell_spec(ip, ci)
    _cdid_cpti_check(ip)
    bsz = min((ip.cells[ci].chunk_rows, budget.max_piece_rows))
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, ip.n_F, keep_rows, budget)
    for (st=1; st<=ip.n_F; st=st+bsz) {
        en = min((ip.n_F, st + bsz - 1))
        pc = _cdid_cpti_eval(ip, ci, (st::en))
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Stream a caller permutation row_order of 1..N_F in batches of piece_size (split and
   reorder verification); the scientific fits are never touched. */
struct cdid_cpt_cell scalar cdid_cpti_build_cell_rows(struct cdid_cpti scalar ip, real scalar ci,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget,
    real colvector row_order, real scalar piece_size)
{
    struct cdid_cpt_acc scalar acc
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar st, en, n_rows, bsz
    spec = cdid_cpti_cell_spec(ip, ci)
    n_rows = rows(row_order)
    if (n_rows != ip.n_F | cols(row_order) != 1 | missing(row_order))
        _cdid_cpti_fail("replay", "Invalid custom row_order or piece_size.", 498)
    if (any(row_order :!= floor(row_order)) | any(row_order :< 1) | any(row_order :> ip.n_F) |
        piece_size < 1 | piece_size != floor(piece_size) | piece_size > ip.cells[ci].chunk_rows)
        _cdid_cpti_fail("replay", "Invalid custom row_order or piece_size.", 498)
    if (rows(uniqrows(row_order)) != ip.n_F)
        _cdid_cpti_fail("replay", "row_order must be a permutation of 1..N_F.", 498)
    bsz = min((piece_size, budget.max_piece_rows))
    acc = cdid_cpt_cell_start(spec, grid, cmap, C, ip.n_F, keep_rows, budget)
    for (st=1; st<=n_rows; st=st+bsz) {
        en = min((n_rows, st + bsz - 1))
        pc = cdid_cpti_piece_rows(ip, ci, row_order[st..en])
        cdid_cpt_cell_add(acc, grid, cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
    }
    return(cdid_cpt_cell_finish(acc, grid))
}

/* Every declared cell of the requested family, or a 459 refusal of the WHOLE test when any
   declared cell is unavailable (never an available-only subset). One common grid or one per
   declared cell. */
struct cdid_cpt_cell colvector cdid_cpti_build_all(struct cdid_cpti scalar ip,
    struct cdid_cpt_grid colvector grids, real colvector cmap, real scalar C,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_cell colvector out
    struct cdid_cpt_acc colvector accs
    struct cdid_cpt_piece scalar pc
    struct cdid_cpt_spec scalar spec
    real scalar i, st, en, bsz, same_bsz
    real colvector rowsF
    transmorphic scalar cache
    if (ip.status != "available" | ip.n_cells < 1 | ip.n_unavailable > 0)
        _cdid_cpti_fail("unavailable", ip.reason, 459)
    if (rows(grids) != 1 & rows(grids) != ip.n_cells)
        _cdid_cpti_fail("grid", "Supply either one common grid or one grid per declared cell.", 498)
    _cdid_cpti_check(ip)
    bsz = min((ip.cells[1].chunk_rows, budget.max_piece_rows))
    same_bsz = 1
    for (i=2; i<=ip.n_cells; i++) {
        if (min((ip.cells[i].chunk_rows, budget.max_piece_rows)) != bsz) same_bsz = 0
    }
    out = cdid_cpt_cell(ip.n_cells, 1)
    if (!same_bsz) {
        for (i=1; i<=ip.n_cells; i++) {
            if (rows(grids) == 1) out[i] = cdid_cpti_build_cell(ip, i, grids[1], cmap, C, keep_rows, budget)
            else out[i] = cdid_cpti_build_cell(ip, i, grids[i], cmap, C, keep_rows, budget)
        }
        return(out)
    }
    accs = cdid_cpt_acc(ip.n_cells, 1)
    for (i=1; i<=ip.n_cells; i++) {
        spec = cdid_cpti_cell_spec(ip, i)
        if (rows(grids) == 1) accs[i] = cdid_cpt_cell_start(spec, grids[1], cmap, C, ip.n_F, keep_rows, budget)
        else accs[i] = cdid_cpt_cell_start(spec, grids[i], cmap, C, ip.n_F, keep_rows, budget)
    }
    for (st=1; st<=ip.n_F; st=st+bsz) {
        en = min((ip.n_F, st + bsz - 1))
        rowsF = (st::en)
        cache = asarray_create("string", 1)
        for (i=1; i<=ip.n_cells; i++) {
            pc = _cdid_cpti_eval(ip, i, rowsF, cache)
            if (rows(grids) == 1) cdid_cpt_cell_add(accs[i], grids[1], cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
            else cdid_cpt_cell_add(accs[i], grids[i], cmap, pc.rows, pc.nu, pc.rho, pc.P, pc.xi, pc.blam)
        }
    }
    for (i=1; i<=ip.n_cells; i++) {
        if (rows(grids) == 1) out[i] = cdid_cpt_cell_finish(accs[i], grids[1])
        else out[i] = cdid_cpt_cell_finish(accs[i], grids[i])
    }
    return(out)
}

/* Oriented scalar h = 1 IF of the backend link on rowsF (frame scale): an
   identity check of the process at h = 1, never a substitute for it. */
real colvector cdid_cpti_provider_if(struct cdid_cpti scalar ip, real scalar ci, real colvector rowsF)
{
    struct cdid_imp_backend_if scalar r1
    struct cdid_aug_backend_if scalar r2
    if (ci < 1 | ci > ip.n_cells | ci != floor(ci))
        _cdid_cpti_fail("replay", "Invalid declared cell index.", 498)
    if (!ip.cells[ci].supported)
        _cdid_cpti_fail("replay", "Requested cell is unavailable: " + ip.cells[ci].reason, 459)
    _cdid_cpti_check(ip)
    if (ip.kind == "augmented") {
        r2 = cdid_aug_backend_rows(*(ip.ab), rowsF, ip.cells[ci].link, "frame", 0)
        return(ip.cells[ci].sign :* r2.phi[, 1])
    }
    r1 = cdid_imp_backend_rows(*(ip.ib), rowsF, ip.cells[ci].link, "frame", 0)
    return(ip.cells[ci].sign :* r1.phi[, 1])
}

void cdid_cpti_describe(struct cdid_cpti scalar ip)
{
    real scalar i
    printf("{txt}conditional %s %s %s: status=%s N_F=%g n_A=%g links=%g primitive=%g post_targets=%g aux=%g cells=%g unavailable=%g\n",
        ip.layout, ip.kind, ip.method, ip.status, ip.n_F, ip.n_A, ip.n_links, ip.n_primitive,
        ip.n_post_targets, ip.n_aux, ip.n_cells, ip.n_unavailable)
    for (i=1; i<=ip.n_cells; i++) {
        printf("{txt}  %s q=%g link=%g sign=%g status=%s terms=%g coords=%g chunks=%g x %g\n",
            ip.cells[i].id, ip.cells[i].q, ip.cells[i].link, ip.cells[i].sign, ip.cells[i].status,
            ip.cells[i].n_terms, ip.cells[i].Dl, ip.cells[i].n_chunks, ip.cells[i].chunk_rows)
    }
}

end
