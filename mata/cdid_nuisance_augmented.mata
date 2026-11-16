*! Opt-in reports of the augmented doubly robust nuisance fits (keepnuisance); no refit.
version 19.5
mata:

/* Pair-cohort lookup of a selection fit in the coordinates of the comparison
   plan: row (s - 1) * T + t, column the plan's cohort index. Matching is by
   calendar and cohort values; never-treated cohorts are coded 0 on both sides. */
real matrix _cdid_nuis_pc_lookup(struct cdid_sel_fit_result scalar fit,
    struct cdid_baseline_plan_state scalar plan)
{
    real scalar nT, pc, g
    real colvector s, t, k, pg
    real matrix L
    nT = rows(plan.calendar)
    L = J(nT * nT, rows(plan.cohorts), 0)
    pg = editmissing(plan.cohorts, 0)
    for (pc = 1; pc <= rows(fit.map.pair_cohorts); pc++) {
        g = editmissing(fit.cohorts[fit.map.pair_cohorts[pc, 1]], 0)
        s = selectindex(plan.calendar :== fit.all_times[fit.map.pair_cohorts[pc, 2]])
        t = selectindex(plan.calendar :== fit.all_times[fit.map.pair_cohorts[pc, 3]])
        k = selectindex(pg :== g)
        if (rows(s) == 1 & rows(t) == 1 & rows(k) == 1) L[(s - 1) * nT + t, k] = pc
    }
    return(L)
}

string scalar _cdid_nuis_pair(real colvector cal, real scalar s, real scalar t)
{
    return("(" + _cdid_nuis_num(cal[s]) + "," + _cdid_nuis_num(cal[t]) + ")")
}

string scalar _cdid_nuis_slopes(string rowvector names, real rowvector kept, real scalar L)
{
    string scalar out
    real scalar j
    if (cols(kept) == 0) return("")
    out = ""
    for (j = 1; j <= cols(kept); j++) out = out + (j > 1 ? " + " : "") + "b" + strofreal(j) + "*" + names[kept[j]]
    return((L == 1 ? "" : _cdid_nuis_num(L) + "*") + "(" + out + ")")
}

/* Capture on the augmented route (selmethod(dr)): the observation logits, and
   with method(dr) the counterfactual outcome regressions, the treatment
   propensity logits and the imputation regressions. */
void cdid_nuis_capture_augmented(struct cdid_aug_backend scalar a,
    struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan)
{
    struct cdid_nuis_cap scalar cap
    string scalar why
    real scalar i
    if (a.is_dated | (a.qfit.selection_universe != "F" & a.qfit.selection_universe != "")) {
        if (a.is_dated) why = "With covariates or weights that vary over time, the augmented correction fits its models separately for each reference period; they are not reported."
        else why = "The augmented correction is reported when its observation logits are fitted on every unit of the study population."
        cap.mode = "chained"
        cap.comp = J(4, 4, "")
        cap.comp[1, ] = ("selection", "unavailable", "reference_recipe_route", why)
        if (a.method == "dr") {
            for (i = 2; i <= 4; i++) cap.comp[i, ] = ((i == 2 ? "outcome" : (i == 3 ? "treatment" : "imputation")),
                "unavailable", "reference_recipe_route", why)
        }
        _cdid_nuis_store_pending(cap)
        return
    }
    cap = _cdid_nuis_sel_capture(a.qfit, *(a.qrdr), _cdid_nuis_pc_lookup(a.qfit, plan),
        plan.calendar, B, plan, 0)
    if (a.method == "dr" & a.has_or & a.has_ps) _cdid_nuis_aug_roles(cap, a, B, plan)
    _cdid_nuis_store_pending(cap)
}

void _cdid_nuis_aug_roles(struct cdid_nuis_cap scalar cap, struct cdid_aug_backend scalar a,
    struct cdid_public_bundle scalar B, struct cdid_baseline_plan_state scalar plan)
{
    pointer(struct cdid_imp_backend scalar) scalar pib
    pointer(struct cdid_imputation_model scalar) scalar pm
    struct cdid_imputation_adjoint scalar adj
    struct cdid_sel_moments scalar sm
    struct _cdid_impb_ivec rowvector c_ix
    struct cdid_nuis_role scalar R
    real scalar N, G, nb, sep, icpt, q, nO, nT, nI, nf, W, f, b, j, k, kk, ck, cid, gid, rtid, rt
    real scalar lo, hi, Bc, chunk, G_cl, clus, row, at, nr, ell, pid, ncoef, slope_role, i
    real colvector lk, abids, bbids, fr, coh, Av, wr, ix, reps, fsize, fstart, boff, cl, rows_k
    real colvector sum1, sum2, var, se, ncl, tmp, fx
    real matrix xi, ora, orb, phi, comp, ref4, X, vals, sd, XA, sa, xb, S, seenc, hitm
    real rowvector kept
    string rowvector snames, tnames
    string scalar ctrl, lab

    pib = a.ipw
    pm = a.model
    N = a.N
    G = a.G
    nb = a.n_blocks
    sep = (a.outcome_pooling == "separate")
    icpt = a.intercept
    q = a.q

    /* Design names: the treatment design and the raw slope columns. */
    tnames = J(1, q, "")
    for (j = 1; j <= q; j++) {
        if (j <= cols(B.prepared.x_names)) tnames[j] = B.prepared.x_names[j]
        else tnames[j] = "x" + strofreal(j)
    }
    if (B.prepared.x_has_intercept & q >= 1) tnames[1] = "(Intercept)"
    snames = J(1, a.p_raw, "")
    for (j = 1; j <= a.p_raw; j++) {
        k = j + (B.prepared.x_has_intercept ? 1 : 0)
        snames[j] = (k <= cols(B.prepared.x_names) ? B.prepared.x_names[k] : "x" + strofreal(j))
    }
    rows_k = J(G, 1, 0)
    for (k = 1; k <= G; k++) rows_k[k] = sum((*pm).unit_cohort_index :== k)

    /* Fit layout: outcome fits, treatment fits (one per distinct risk set), imputation groups. */
    nO = (sep ? nb : 1)
    reps = J(0, 1, .)
    for (rtid = 1; rtid <= (*pib).n_trt; rtid++) if ((*pib).trts[rtid].rep_tid == rtid) reps = reps \ rtid
    nT = rows(reps)
    nI = (*pm).n_groups
    nf = nO + nT + nI
    fsize = J(nf, 1, 0)
    boff = J(nb, 1, 0)
    if (sep) {
        at = 0
        for (b = 1; b <= nb; b++) {
            boff[b] = at
            at = at + a.blocks[b].or_p
            fsize[b] = (icpt ? 1 : 0) + a.blocks[b].or_p
        }
    }
    else fsize[1] = (icpt ? nb : 0) + a.or_p
    for (f = 1; f <= nT; f++) fsize[nO + f] = q
    for (gid = 1; gid <= nI; gid++) fsize[nO + nT + gid] = (icpt ? (*pm).groups[gid].n_cells : 0) + (*pm).groups[gid].p
    fstart = runningsum(fsize) - fsize :+ 1
    W = sum(fsize)

    /* One pass over the units: coefficient influence functions of every fit. */
    clus = B.prepared.clustered
    G_cl = (clus ? B.prepared.C : N)
    cl = (clus ? B.prepared.unit_cluster : J(0, 1, .))
    sum1 = J(1, W, 0)
    sum2 = J(1, W, 0)
    S = (clus ? J(G_cl, W, 0) : J(0, 0, .))
    seenc = (clus ? J(G_cl, nf, 0) : J(0, 0, .))
    lk = 1
    abids = (icpt ? (1::nb) : J(0, 1, 0))
    bbids = (sep ? (1::nb) : J(0, 1, 0))
    adj = _cdid_aug_batch_adj(a, *pm, lk, abids)
    adj = _cdid_aug_cat_or_adj(a, *pm, adj, bbids)
    sm = _cdid_aug_batch_sel(a, lk, abids, bbids)
    chunk = max((1, a.chunk_rows))
    for (lo = 1; lo <= N; lo = lo + chunk) {
        hi = min((N, lo + chunk - 1))
        fr = (lo::hi)
        Bc = rows(fr)
        _cdid_aug_eval(a, *pib, *pm, *(a.qrdr), fr, lk, abids, bbids, adj, sm, 0,
            phi = ., comp = ., ref4 = ., ora = ., orb = .)
        xi = J(Bc, W, 0)
        /* Outcome. */
        if (sep) {
            for (b = 1; b <= nb; b++) {
                at = fstart[b]
                if (icpt) {
                    xi[, at] = ora[, b]
                    at++
                }
                if (a.blocks[b].or_p > 0) xi[, at..(at + a.blocks[b].or_p - 1)] = orb[, (boff[b] + 1)..(boff[b] + a.blocks[b].or_p)]
            }
        }
        else {
            at = fstart[1]
            if (icpt) {
                xi[, at..(at + nb - 1)] = ora[, 1..nb]
                at = at + nb
            }
            if (a.or_p > 0) xi[, at..(at + a.or_p - 1)] = orb[, 1..a.or_p]
        }
        /* Treatment: weighted logit score plus the repeat-rate channel. */
        coh = (*pm).unit_cohort_index[fr]
        c_ix = _cdid_impb_bucket_coh(coh, G)
        Av = (*pm).pop.A[fr]
        if (q > 0) for (f = 1; f <= nT; f++) {
            rtid = reps[f]
            rt = (*pib).trts[rtid].recipe
            sd = J(Bc, q, 0)
            XA = J(Bc, (*pib).trts[rtid].n_coh, 0)
            for (k = 1; k <= (*pib).trts[rtid].n_coh; k++) {
                ck = (*pib).trts[rtid].all_ci[k]
                if (!((*pib).mass_F_mat[rt, ck] > 0)) continue
                ix = c_ix[ck].rows
                if (rows(ix) == 0) continue
                wr = _cdid_impb_w(*pib, *pm, rt, fr[ix], 1)
                XA[ix, k] = wr :* (Av[ix] :- (*pib).a_k_mat[rt, ck]) :/ (*pib).mass_F_mat[rt, ck]
                X = _cdid_impb_x(*pib, rt, fr[ix], 1)
                vals = _cdid_aug_logit_vals(X * (*pib).trts[rtid].gamma)
                sd[ix, ] = ((wr :* (*pib).a_k_mat[rt, ck]) :* ((ck == (*pib).trts[rtid].target_ci) :- vals[, 1])) :* X
            }
            sa = XA * (*pib).trts[rtid].a_forcing'
            at = fstart[nO + f]
            xi[, at..(at + q - 1)] = _cdid_nuisance_gram((*pib).trts[rtid].R, (*pib).trts[rtid].pivot, N :* sd')' +
                _cdid_nuisance_gram((*pib).trts[rtid].R, (*pib).trts[rtid].pivot, N :* sa')'
        }
        /* Imputation: observed-pair weighted least-squares scores. */
        for (gid = 1; gid <= nI; gid++) {
            at = fstart[nO + nT + gid]
            xb = ((*pm).groups[gid].p > 0 ? cdid_imputation_beta_if_chunk(*pm, gid, fr, "frame") : J(Bc, 0, 0))
            if (icpt) {
                for (kk = 1; kk <= (*pm).groups[gid].n_cells; kk++) {
                    cid = (*pm).groups[gid].cell_ids[kk]
                    xi[, at] = cdid_imputation_cell_score_chunk(*pm, cid, fr) :/ (*pm).cells[cid].mass
                    if ((*pm).groups[gid].p > 0) xi[, at] = xi[, at] - xb * (*pm).cells[cid].mean_d'
                    at++
                }
            }
            if ((*pm).groups[gid].p > 0) xi[, at..(at + (*pm).groups[gid].p - 1)] = xb
        }
        if (clus) {
            hitm = J(Bc, nf, 0)
            for (f = 1; f <= nf; f++) if (fsize[f] > 0) hitm[, f] = rowsum(xi[, fstart[f]..(fstart[f] + fsize[f] - 1)] :!= 0) :> 0
            for (row = 1; row <= Bc; row++) {
                S[cl[fr[row]], ] = S[cl[fr[row]], ] + xi[row, ]
                seenc[cl[fr[row]], ] = seenc[cl[fr[row]], ] :| hitm[row, ]
            }
        }
        else {
            sum1 = sum1 + colsum(xi)
            sum2 = sum2 + colsum(xi :* xi)
        }
    }
    if (clus) {
        sum1 = colsum(S)
        sum2 = colsum(S :* S)
        ncl = colsum(seenc)'
    }
    else ncl = J(nf, 1, .)
    var = ((sum2 - (sum1 :* sum1) :/ G_cl) :/ (N * N))'
    se = sqrt(var :* (var :> 0))
    if (clus) {
        for (f = 1; f <= nf; f++) if (ncl[f] < 2 & fsize[f] > 0) se[fstart[f]..(fstart[f] + fsize[f] - 1)] = J(fsize[f], 1, .)
    }

    cap.roles = cdid_nuis_role(1, 3)

    /* Outcome tables. */
    R.role = "outcome"
    R.pooling = (sep ? "separate" : "pooled")
    R.blocks = J(nO, 8, .)
    R.block_txt = J(nO, 6, "")
    R.coef = J(sum(fsize[1..nO]), 6, .)
    R.coef_txt = J(rows(R.coef), 2, "")
    row = 0
    for (f = 1; f <= nO; f++) {
        R.block_txt[f, 1] = "O" + strofreal(f)
        if (sep) {
            b = f
            pid = a.blocks[b].pair_id
            ctrl = _cdid_nuis_glist((*pm).cohorts[a.blocks[b].ctrl_ci], a.blocks[b].n_ctrl)
            R.block_txt[f, 2] = "block[" + a.blocks[b].name + "]"
            R.block_txt[f, 3] = "E[Z_i(k;" + _cdid_nuis_num(plan.calendar[(*pm).pairs[pid].s_index]) + "," +
                _cdid_nuis_num(plan.calendar[(*pm).pairs[pid].t_index]) + ") | X_i, G_i = k] = " +
                (icpt ? "a" : "0") + (a.blocks[b].or_p > 0 ? " + " + _cdid_nuis_slopes(snames, a.blocks[b].or_keep, a.blocks[b].L) : "") +
                ",  k in {" + ctrl + "}"
            R.block_txt[f, 4] = ctrl
            R.block_txt[f, 5] = (a.blocks[b].or_p > 0 ? invtokens(snames[a.blocks[b].or_keep], ",") : "")
            R.blocks[f, ] = (fsize[f], icpt, a.blocks[b].or_p, sum(rows_k[a.blocks[b].ctrl_ci]), a.blocks[b].weight_sum,
                ncl[f], 0, b)
            if (icpt) {
                row++
                R.coef[row, (1, 2, 3)] = (f, b, a.blocks[b].alpha)
                R.coef_txt[row, ] = ("(Intercept)[" + _cdid_nuis_pair(plan.calendar, (*pm).pairs[pid].s_index, (*pm).pairs[pid].t_index) +
                    ",K={" + ctrl + "}]", "local_intercept")
            }
            for (j = 1; j <= a.blocks[b].or_p; j++) {
                row++
                R.coef[row, (1, 3)] = (f, a.blocks[b].or_beta[j])
                R.coef_txt[row, ] = (snames[a.blocks[b].or_keep[j]], "slope")
            }
        }
        else {
            R.block_txt[f, 2] = "pooled[" + strofreal(nb) + " control block(s), common slope]"
            R.block_txt[f, 3] = "E[Z_i(k;s,t) | X_i, G_i = k] = " + (icpt ? "a[b]" : "0") +
                (a.or_p > 0 ? " + " + subinstr(_cdid_nuis_slopes(snames, a.or_keep, 1), "(", "L(s,t)*(", 1) : "") +
                ",  k in K_b for each control block b"
            tmp = J(0, 1, .)
            for (b = 1; b <= nb; b++) tmp = tmp \ a.blocks[b].ctrl_ci
            tmp = uniqrows(tmp)
            R.block_txt[f, 4] = _cdid_nuis_glist((*pm).cohorts[tmp], rows(tmp))
            R.block_txt[f, 5] = (a.or_p > 0 ? invtokens(snames[a.or_keep], ",") : "")
            R.blocks[f, ] = (fsize[f], (icpt ? nb : 0), a.or_p, 0, 0, ncl[f], sum(plan.raw_supported :!= 0), 1)
            for (b = 1; b <= nb; b++) {
                R.blocks[f, 4] = R.blocks[f, 4] + sum(rows_k[a.blocks[b].ctrl_ci])
                R.blocks[f, 5] = R.blocks[f, 5] + a.blocks[b].weight_sum
                if (icpt) {
                    pid = a.blocks[b].pair_id
                    row++
                    R.coef[row, (1, 2, 3)] = (f, b, a.blocks[b].alpha)
                    R.coef_txt[row, ] = ("(Intercept)[" + _cdid_nuis_pair(plan.calendar, (*pm).pairs[pid].s_index, (*pm).pairs[pid].t_index) +
                        ",K={" + _cdid_nuis_glist((*pm).cohorts[a.blocks[b].ctrl_ci], a.blocks[b].n_ctrl) + "}]", "local_intercept")
                }
            }
            for (j = 1; j <= a.or_p; j++) {
                row++
                R.coef[row, (1, 3)] = (f, a.or_beta[j])
                R.coef_txt[row, ] = (snames[a.or_keep[j]], "shared_slope")
            }
        }
    }
    if (sep) for (i = 1; i <= rows(plan.raw); i++) if (plan.raw_supported[i]) {
        ell = plan.raw_consumer[i]
        R.blocks[a.links[ell].or_block, 7] = R.blocks[a.links[ell].or_block, 7] + 1
    }
    R.coef[, 4] = se[1..rows(R.coef)]
    cap.roles[1] = R

    /* Treatment tables. */
    R.role = "treatment"
    R.pooling = "risk_set"
    R.blocks = J(nT, 8, .)
    R.block_txt = J(nT, 6, "")
    R.coef = J(nT * q, 6, .)
    R.coef_txt = J(nT * q, 2, "")
    row = 0
    for (f = 1; f <= nT; f++) {
        rtid = reps[f]
        rt = (*pib).trts[rtid].recipe
        tmp = select((*pib).trts[rtid].all_ci, (*pib).trts[rtid].all_ci :!= (*pib).trts[rtid].target_ci)
        ctrl = _cdid_nuis_glist((*pm).cohorts[tmp], rows(tmp))
        R.block_txt[f, 1] = "T" + strofreal(f)
        R.block_txt[f, 2] = "risk_set[g=" + _cdid_nuis_g((*pm).cohorts[(*pib).trts[rtid].target_ci]) + " vs K={" + ctrl + "}]"
        lab = ""
        for (j = 1; j <= q; j++) lab = lab + (j > 1 ? " + " : "") + "c" + strofreal(j - 1) + (tnames[j] == "(Intercept)" ? "" : "*" + tnames[j])
        R.block_txt[f, 3] = "logit Pr(G_i = " + _cdid_nuis_g((*pm).cohorts[(*pib).trts[rtid].target_ci]) + " | X_i, G_i in {" +
            _cdid_nuis_g((*pm).cohorts[(*pib).trts[rtid].target_ci]) + (ctrl != "" ? "," + ctrl : "") + "}) = " + lab
        R.block_txt[f, 4] = _cdid_nuis_glist((*pm).cohorts[(*pib).trts[rtid].all_ci], (*pib).trts[rtid].n_coh)
        R.block_txt[f, 5] = invtokens(select(tnames, tnames :!= "(Intercept)"), ",")
        R.blocks[f, ] = (q, sum(tnames :== "(Intercept)"), q - sum(tnames :== "(Intercept)"), (*pib).trts[rtid].n_risk,
            sum((*pib).a_k_mat[rt, (*pib).trts[rtid].all_ci'] :* (*pib).mass_F_mat[rt, (*pib).trts[rtid].all_ci']) * N,
            ncl[nO + f], 0, rtid)
        for (j = 1; j <= q; j++) {
            row++
            R.coef[row, (1, 3, 4)] = (f, (*pib).trts[rtid].gamma[j], se[fstart[nO + f] + j - 1])
            R.coef_txt[row, ] = (tnames[j], (tnames[j] == "(Intercept)" ? "intercept" : "slope"))
        }
    }
    for (i = 1; i <= rows(plan.raw); i++) {
        if (!plan.raw_supported[i]) continue
        ell = plan.raw_consumer[i]
        rtid = (*pib).trts[a.links[ell].trt_id].rep_tid
        fx = selectindex(reps :== rtid)
        if (rows(fx) == 1) R.blocks[fx[1], 7] = R.blocks[fx[1], 7] + 1
    }
    cap.roles[2] = R

    /* Imputation tables. */
    R.role = "imputation"
    R.pooling = (*pm).pooling
    R.blocks = J(nI, 8, .)
    R.block_txt = J(nI, 6, "")
    ncoef = sum(fsize[(nO + nT + 1)..nf])
    R.coef = J(ncoef, 6, .)
    R.coef_txt = J(ncoef, 2, "")
    row = 0
    for (gid = 1; gid <= nI; gid++) {
        kept = (*pm).groups[gid].kept
        nr = (*pm).groups[gid].n_cells
        lab = ""
        tmp = J(nr, 1, .)
        for (kk = 1; kk <= nr; kk++) {
            cid = (*pm).groups[gid].cell_ids[kk]
            pid = (*pm).cells[cid].pair_id
            tmp[kk] = (*pm).cells[cid].cohort_index
            lab = lab + (kk > 1 ? "; " : "") + "g=" + _cdid_nuis_g((*pm).cells[cid].cohort) + "," +
                _cdid_nuis_pair(plan.calendar, (*pm).pairs[pid].s_index, (*pm).pairs[pid].t_index)
        }
        tmp = uniqrows(tmp)
        R.block_txt[gid, 1] = "I" + strofreal(gid)
        R.block_txt[gid, 2] = (*pm).groups[gid].pooling + "[" + (*pm).groups[gid].name + ": " + lab + "]"
        cid = (*pm).groups[gid].cell_ids[1]
        pid = (*pm).cells[cid].pair_id
        if (nr == 1) R.block_txt[gid, 3] = "E[Y_i" + _cdid_nuis_num(plan.calendar[(*pm).pairs[pid].t_index]) + " - Y_i" +
            _cdid_nuis_num(plan.calendar[(*pm).pairs[pid].s_index]) + " | X_i, G_i = " + _cdid_nuis_g((*pm).cells[cid].cohort) +
            ", S_i" + _cdid_nuis_pair(plan.calendar, (*pm).pairs[pid].s_index, (*pm).pairs[pid].t_index) + " = 1] = " +
            (icpt ? "a" : "0") + ((*pm).groups[gid].p > 0 ? " + " + _cdid_nuis_slopes(snames, kept, (*pm).cells[cid].L) : "")
        else R.block_txt[gid, 3] = "E[Y_it - Y_is | X_i, G_i = g, S_i(s,t) = 1] = " +
            (icpt ? (rows(tmp) == 1 ? "a[(s,t)]" : "a[g,(s,t)]") : "0") +
            ((*pm).groups[gid].p > 0 ? " + " + subinstr(_cdid_nuis_slopes(snames, kept, 1), "(", "L(s,t)*(", 1) : "")
        R.block_txt[gid, 4] = _cdid_nuis_glist((*pm).cohorts[tmp], rows(tmp))
        R.block_txt[gid, 5] = ((*pm).groups[gid].p > 0 ? invtokens(snames[kept], ",") : "")
        R.blocks[gid, ] = (fsize[nO + nT + gid], (icpt ? nr : 0), (*pm).groups[gid].p, 0, 0, ncl[nO + nT + gid], 0, gid)
        for (kk = 1; kk <= nr; kk++) {
            cid = (*pm).groups[gid].cell_ids[kk]
            R.blocks[gid, 4] = R.blocks[gid, 4] + (*pm).cells[cid].n_pair
            R.blocks[gid, 5] = R.blocks[gid, 5] + (*pm).cells[cid].weight_sum
        }
        at = fstart[nO + nT + gid]
        slope_role = (nr > 1)
        if (icpt) {
            for (kk = 1; kk <= nr; kk++) {
                cid = (*pm).groups[gid].cell_ids[kk]
                pid = (*pm).cells[cid].pair_id
                row++
                R.coef[row, ] = (gid, cid, (*pm).cells[cid].alpha, se[at], ., .)
                R.coef_txt[row, ] = ("(Intercept)[g=" + _cdid_nuis_g((*pm).cells[cid].cohort) + "," +
                    _cdid_nuis_pair(plan.calendar, (*pm).pairs[pid].s_index, (*pm).pairs[pid].t_index) + "]", "local_intercept")
                at++
            }
        }
        for (j = 1; j <= (*pm).groups[gid].p; j++) {
            row++
            R.coef[row, ] = (gid, ., (*pm).groups[gid].beta[j], se[at], ., .)
            R.coef_txt[row, ] = (snames[kept[j]], (slope_role ? "shared_slope" : "slope"))
            at++
        }
    }
    for (i = 1; i <= rows(plan.raw); i++) {
        if (!plan.raw_supported[i]) continue
        ell = plan.raw_consumer[i]
        tmp = (*pm).cells[a.links[ell].target_cell].group_id
        for (k = 1; k <= a.links[ell].n_ctrl; k++) tmp = tmp \ (*pm).cells[a.links[ell].ctrl_cell[k]].group_id
        tmp = uniqrows(tmp)
        for (k = 1; k <= rows(tmp); k++) R.blocks[tmp[k], 7] = R.blocks[tmp[k], 7] + 1
    }
    cap.roles[3] = R

    for (f = 1; f <= 3; f++) {
        cap.roles[f].coef[, 5] = cap.roles[f].coef[, 3] :/ cap.roles[f].coef[, 4]
        cap.roles[f].coef[, 6] = 2 :* normal(-abs(cap.roles[f].coef[, 5]))
    }
    cap.comp[2, ] = ("outcome", "available", "", "")
    cap.comp[3, ] = ("treatment", "available", "", "")
    cap.comp[4, ] = ("imputation", "available", "", "")
}

end
