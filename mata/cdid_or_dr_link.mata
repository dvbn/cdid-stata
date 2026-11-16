*! Baseline pooled untreated-change OR and treatment DR, current Stata.
*! Pinned R R cdid: chained_or_baseline/pooled.R, chained_dr_backend.R.
*! Load cdid_nuisance.mata and cdid_ipw_link.mata first. No selection fitting.
version 19.0
mata:
struct cdid_or_recipe {
    struct cdid_ipw_frame scalar frame
    real matrix Z
}
struct cdid_or_pair {
    real scalar s_index, t_index, recipe
    real colvector a, dy, pos_a
}
struct _cdid_or_summary {
    real matrix F
}
struct cdid_or_block_fit {
    real rowvector kept, dropped, pivot
    real colvector beta
    real matrix R, alias_map, alias_constants, mean_D, beta_rate_loading
    real rowvector score_mean_beta
    real scalar score_mean_alpha
}
struct cdid_or_backend {
    string scalar outcome_pooling
    struct cdid_or_recipe rowvector recipes
    struct cdid_or_pair rowvector pairs
    struct cdid_or_block_fit rowvector blocks
    real scalar n_A, n_F, n_supports, n_blocks, n_links, intercept, rank_tol, chunk_rows
    real scalar n_fit_rows, objective
    real matrix support, links, sumZ, mean_D, R, alias_map, alias_constants
    real matrix beta_score_A, beta_rate_loading
    real colvector block_pair, block_start, block_support, rates, rates_c, rates_active
    real colvector pair_mass, sumY, block_mass, alpha, beta, alpha_rate
    real colvector link_mu, link_d, link_numerator, link_rate_loading
    real matrix link_c
    real rowvector kept, dropped, pivot
}
struct cdid_or_if_state {
    string scalar outcome_pooling
    real scalar first_full, last_full, n_A, n_F, n_supports, n_blocks
    real matrix beta
    real colvector fitted_beta
}
struct cdid_or_if_chunk {
    real colvector full_index, phi
    real matrix components
}
struct cdid_dr_link {
    struct cdid_ipw_link scalar residual_ipw
    real scalar link_id, estimate, estimate_or, estimate_ipw, denominator_plus, denominator_minus
    real colvector d_bar_plus, d_bar_minus, d_bar_diff
}
void _cdid_or_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}
real colvector _cdid_or_index(real matrix x)
{
    return(vec(selectindex(x)))
}
/* Bounded augmented QR factor in original column order, including response.
   Padding zero rows changes no cross-product and handles tiny support cells. */
real matrix _cdid_or_compress(real matrix A)
{
    real matrix Q, R
    real rowvector tau, pivot
    real scalar k
    k=cols(A)
    if(rows(A)==0) return(J(0,k,0))
    if(rows(A)<k) A=A\J(k-rows(A),k,0)
    if(missing(A)) _cdid_or_fail("Nonfinite weighted support rows or factor.",504)
    hqrdp(A,Q=.,tau=.,R=.,pivot=.)
    if(missing(R)) _cdid_or_fail("Nonfinite augmented support QR.",504)
    return(R[,invorder(pivot)])
}
real scalar _cdid_or_pair_pos(struct cdid_or_pair scalar pair, real scalar ai)
{
    real scalar k
    if(rows(pair.pos_a)>=ai & ai>=1) return(pair.pos_a[ai])
    k=_cdid_ipw_lower_bound(pair.a,ai)
    if(k>rows(pair.a)) return(0)
    if(pair.a[k]!=ai) return(0)
    return(k)
}
real colvector _cdid_or_pair_pos_vec(struct cdid_or_pair scalar pair, real colvector ai, real scalar n_A)
{
    real colvector pos
    real scalar i
    if(rows(pair.pos_a)==n_A) return(pair.pos_a[ai])
    pos=J(rows(ai),1,0)
    for(i=1;i<=rows(ai);i++) pos[i]=_cdid_or_pair_pos(pair,ai[i])
    return(pos)
}
real scalar _cdid_or_lex(real colvector a, real colvector b)
{
    real scalar i
    for(i=1;i<=min((rows(a),rows(b)));i++) {
        if(a[i]<b[i]) return(-1)
        if(a[i]>b[i]) return(1)
    }
    return(sign(rows(a)-rows(b)))
}
/* One stored posting per caller-supplied canonical pair/recipe identity,
   not one observer array per link. Array indices carry opaque pair identity;
   equal numeric endpoints can denote distinct source pair objects.
   support is S x2: pair-id, global cohort-id; sorted by those numeric keys.
   CSR blocks: block_pair[b], members[start[b]..start[b+1]-1]. Each member is
   a support-id. Blocks are pair-major, then lexicographic sorted cohort sets;
   duplicates refuse so every scientific block has objective multiplier one.
   links is L x2: block-id,target-support-id; repeated consumers are permitted.
   Recipes have the same original F/A map and global cohort IDs, but may carry
   different fixed designs/weights. Z excludes the explicit intercept column.
   Pair endpoints are CALENDAR INDICES; L=t_index-s_index, never label distance.
   rank_tol and fit-time chunk_rows are explicit, with no hidden work budget. */
struct cdid_or_backend scalar cdid_or_prepare(struct cdid_or_recipe rowvector recipes,
    struct cdid_or_pair rowvector pairs, real matrix support,
    real colvector block_pair, real colvector block_start, real colvector block_support,
    real matrix links, real scalar intercept, real scalar rank_tol, real scalar chunk_rows,
    | string scalar outcome_pooling_in)
{
    struct cdid_or_backend scalar out
    struct cdid_wls_result scalar fit
    struct _cdid_or_summary rowvector summaries
    struct cdid_ipw_frame scalar fr
    string scalar outcome_pooling
    real scalar nr, npair, p, k, kb, max_k, S, B, M, Jl, r, pi, s, b, m, j, lo, hi, nn, sid, L
    real scalar mass, first, nPseudo, at, dres, scale, sum_sa, n_ctrl_b
    real matrix A, F, Q, R, V, residual_factor, pseudo_X, coef_rows, force_beta, rawZ, D, centered
    real matrix xmin, xmax, firstZ
    real rowvector tau, constant_columns, sum_sb
    real colvector ix, ai, ww, yy, pseudo_y, pseudo_w, pseudo_block, oldpool, newpool
    real colvector rr, raw_res, fw, sa
    outcome_pooling=(args()<11 ? "" : strtrim(outcome_pooling_in))
    if(outcome_pooling=="") outcome_pooling="separate"
    if(outcome_pooling!="pooled" & outcome_pooling!="separate")
        _cdid_or_fail("Invalid outcome_pooling control.",498)
    nr=cols(recipes); npair=cols(pairs); S=rows(support); B=rows(block_pair); M=rows(block_support); Jl=rows(links)
    if(nr<1 | npair<1 | S<1 | B<1 | Jl<1 | cols(support)!=2 | cols(links)!=2 |
        cols(block_pair)!=1 | cols(block_start)!=1 | rows(block_start)!=B+1 | cols(block_support)!=1 |
        missing(support) | missing(links) | missing(block_pair) | missing(block_start) | missing(block_support) |
        any(support:!=floor(support)) | any(links:!=floor(links)) | any(block_pair:!=floor(block_pair)) |
        any(block_start:!=floor(block_start)) | any(block_support:!=floor(block_support)) |
        any(block_pair:<1) | any(block_pair:>npair) | any(block_support:<1) | any(block_support:>S) |
        block_start[1]!=1 | block_start[B+1]!=M+1 | any(block_start[2..B+1]:<=block_start[1..B]) |
        any(support[,1]:<1) | any(support[,1]:>npair) | any(support[,2]:<1) |
        any(links[,1]:<1) | any(links[,1]:>B) | any(links[,2]:<1) | any(links[,2]:>S) |
        (intercept!=0 & intercept!=1) | rank_tol<=0 | rank_tol>=1 |
        chunk_rows<1 | chunk_rows>=. | chunk_rows!=floor(chunk_rows))
        _cdid_or_fail("Invalid explicit canonical support/block/link inputs or execution controls.",498)
    fr=recipes[1].frame; p=cols(recipes[1].Z)
    if(!intercept & p==0) _cdid_or_fail("No-intercept outcome design is empty.",459)
    for(r=1;r<=nr;r++) {
        if(recipes[r].frame.n_A!=fr.n_A | recipes[r].frame.n_F!=fr.n_F |
            recipes[r].frame.n_cohorts!=fr.n_cohorts | any(recipes[r].frame.group:!=fr.group) |
            any(recipes[r].frame.a_to_full:!=fr.a_to_full) | rows(recipes[r].Z)!=fr.n_A |
            cols(recipes[r].Z)!=p | missing(recipes[r].Z) | any(fr.group:<1))
            _cdid_or_fail("Recipes require one common known-cohort F/A map and matching outcome columns.",498)
    }
    for(pi=1;pi<=npair;pi++) {
        r=pairs[pi].recipe
        if(r<1 | r>nr | r!=floor(r) | pairs[pi].s_index<1 |
            pairs[pi].t_index<=pairs[pi].s_index | pairs[pi].t_index>=. |
            pairs[pi].s_index!=floor(pairs[pi].s_index) | pairs[pi].t_index!=floor(pairs[pi].t_index) |
            rows(pairs[pi].a)<1 | cols(pairs[pi].a)!=1 | cols(pairs[pi].dy)!=1 |
            rows(pairs[pi].dy)!=rows(pairs[pi].a) | missing(pairs[pi].a) | missing(pairs[pi].dy) |
            any(pairs[pi].a:!=floor(pairs[pi].a)) | any(pairs[pi].a:<1) | any(pairs[pi].a:>fr.n_A))
            _cdid_or_fail("Invalid compact actual-pair postings or calendar indices.",498)
        nn=rows(pairs[pi].a)
        if(nn>1) if(any(pairs[pi].a[2..nn]:<=pairs[pi].a[1..nn-1]))
            _cdid_or_fail("Pair A indices must be strictly increasing.",498)
        if(pi>1) {
            if(_cdid_or_lex((pairs[pi-1].s_index\pairs[pi-1].t_index\pairs[pi-1].recipe),
                (pairs[pi].s_index\pairs[pi].t_index\pairs[pi].recipe))>0)
                _cdid_or_fail("Caller-supplied pair/design block identities must follow nondecreasing calendar order.",498)
        }
    }
    for(s=1;s<=S;s++) {
        if(support[s,2]>fr.n_cohorts) _cdid_or_fail("Support cohort outside common dictionary.",498)
        if(s>1) if(_cdid_or_lex(support[s-1,]',support[s,]')>=0)
            _cdid_or_fail("Support identities must be unique and canonical.",498)
    }
    oldpool=J(0,1,0)
    for(b=1;b<=B;b++) {
        sa=block_support[block_start[b]..block_start[b+1]-1]
        if(any(support[sa,1]:!=block_pair[b])) _cdid_or_fail("Control support belongs to a different pair/design block.",498)
        newpool=support[sa,2]
        if(rows(sa)>1) if(any(newpool[2..rows(sa)]:<=newpool[1..rows(sa)-1]))
            _cdid_or_fail("Control cohorts must be sorted and unique within each block.",498)
        if(b>1) {
            if(block_pair[b]<block_pair[b-1]) _cdid_or_fail("Control blocks must be pair-major.",498)
            if(block_pair[b]==block_pair[b-1]) if(_cdid_or_lex(oldpool,newpool)>=0)
                _cdid_or_fail("Duplicate or noncanonical scientific control block; no objective replication is allowed.",498)
        }
        oldpool=newpool
    }
    for(j=1;j<=Jl;j++) {
        b=links[j,1]; sid=links[j,2]
        if(support[sid,1]!=block_pair[b] | any(block_support[block_start[b]..block_start[b+1]-1]:==sid))
            _cdid_or_fail("Target support must match its block pair and lie outside the control pool.",498)
    }
    out.outcome_pooling=outcome_pooling
    out.recipes=recipes; out.pairs=pairs; out.support=support; out.links=links
    out.block_pair=block_pair; out.block_start=block_start; out.block_support=block_support
    out.n_A=fr.n_A; out.n_F=fr.n_F; out.n_supports=S; out.n_blocks=B; out.n_links=Jl
    for(pi=1;pi<=npair;pi++) {
        out.pairs[pi].pos_a=J(out.n_A,1,0)
        out.pairs[pi].pos_a[out.pairs[pi].a]=(1::rows(out.pairs[pi].a))
    }
    out.intercept=intercept; out.rank_tol=rank_tol; out.chunk_rows=chunk_rows
    out.pair_mass=J(S,1,0); out.sumY=J(S,1,0); out.sumZ=J(S,p,0)
    out.rates=J(S,1,0); xmin=J(S,p,.); xmax=xmin; firstZ=xmin
    summaries=_cdid_or_summary(1,S)
    for(s=1;s<=S;s++) {
        pi=support[s,1]; r=pairs[pi].recipe; fr=recipes[r].frame; L=pairs[pi].t_index-pairs[pi].s_index
        /* Canonical mass uses the same whole filtered, increasing-A vector
           sum as IPW. QR chunk boundaries must not change rates. */
        ix=_cdid_or_index(fr.group[pairs[pi].a,1]:==support[s,2])
        ai=pairs[pi].a[ix,1]
        out.pair_mass[s]=sum(fr.weight[ai,1])
        F=J(0,p+2,0); first=1; nn=rows(pairs[pi].a)
        for(lo=1;lo<=nn;lo=lo+chunk_rows) {
            hi=min((nn,lo+chunk_rows-1)); rr=(lo..hi)'
            ix=_cdid_or_index(fr.group[pairs[pi].a[rr],1]:==support[s,2])
            if(rows(ix)==0) continue
            rr=rr[ix,1]; ai=pairs[pi].a[rr,1]; ww=fr.weight[ai,1]; yy=pairs[pi].dy[rr,1]
            rawZ=recipes[r].Z[ai,]; D=L:*rawZ
            if(missing(D)) _cdid_or_fail("Calendar-scaled outcome design overflowed.",504)
            if(p) {
                if(first) { xmin[s,]=colmin(rawZ); xmax[s,]=colmax(rawZ); firstZ[s,]=rawZ[1,]; }
                else { xmin[s,]=colmin(xmin[s,]\rawZ); xmax[s,]=colmax(xmax[s,]\rawZ); }
            }
            first=0; out.sumY[s]=out.sumY[s]+sum(ww:*yy)
            if(p) out.sumZ[s,]=out.sumZ[s,]+quadcross(ww,rawZ)
            F=_cdid_or_compress(F\((J(rows(ai),1,1),D,yy):*sqrt(ww)))
        }
        if(out.pair_mass[s]<=0 | out.pair_mass[s]>=.)
            _cdid_or_fail("Every declared pair/cohort support requires positive observed mass; none is dropped.",459)
        out.rates[s]=out.pair_mass[s]/fr.group_mass[support[s,2]]
        summaries[s].F=F
    }
    A=cdid_ipw_rate_values(out.rates); out.rates_c=A[,1]; out.rates_active=A[,2]
    /* A bounded artificial row representation of each exact augmented QR:
       paired +/- centered rows preserve mass, means and centered Gram/response.
       Its size is 2(p+1) per block, independent of original-unit count. */
    nPseudo=2*(p+1)*B
    pseudo_X=J(nPseudo,p,0); pseudo_y=J(nPseudo,1,0); pseudo_w=J(nPseudo,1,0); pseudo_block=J(nPseudo,1,0)
    at=0
    for(b=1;b<=B;b++) {
        pi=block_pair[b]; L=pairs[pi].t_index-pairs[pi].s_index
        F=J(0,p+2,0); mass=0; constant_columns=J(1,p,1); rawZ=J(1,p,0)
        for(m=block_start[b];m<block_start[b+1];m++) {
            sid=block_support[m]; F=_cdid_or_compress(F\(summaries[sid].F/sqrt(out.rates_c[sid])))
            mass=mass+out.pair_mass[sid]/out.rates_c[sid]
            if(m==block_start[b]) rawZ=firstZ[sid,]
            if(p) constant_columns=constant_columns:*(xmin[sid,]:==rawZ):*(xmax[sid,]:==rawZ)
        }
        hqrd(F[,1],Q=.,tau=.,R=.)
        V=hqrdmultq(Q,tau,F[,(2..p+2)],1)
        coef_rows=V[1,]/R[1,1]
        residual_factor=V[2..rows(V),]
        if(intercept & p) for(k=1;k<=p;k++) if(constant_columns[k]) {
            coef_rows[k]=L*rawZ[k]; residual_factor[,k]=J(rows(residual_factor),1,0)
        }
        nn=rows(residual_factor); scale=sqrt(nn/mass)
        A=coef_rows:+scale:*residual_factor
        V=coef_rows:-scale:*residual_factor
        rr=(at+1..at+2*nn)'
        if(p) pseudo_X[rr,]=(A[,1..p]\V[,1..p])
        pseudo_y[rr]=(A[,p+1]\V[,p+1]); pseudo_w[rr]=J(2*nn,1,mass/(2*nn)); pseudo_block[rr]=J(2*nn,1,b)
        at=at+2*nn
    }
    if(out.outcome_pooling=="separate") {
        out.blocks=cdid_or_block_fit(1,B)
        out.alpha=J(B,1,0); out.block_mass=J(B,1,0); out.beta=J(0,1,0)
        out.R=J(0,0,0); out.pivot=J(1,0,0); out.kept=J(1,0,0); out.dropped=J(1,0,0)
        out.alias_map=J(0,0,0); out.alias_constants=J(0,0,0); out.mean_D=J(B,0,0)
        out.beta_score_A=J(out.n_A,0,0); out.beta_rate_loading=J(0,S,0)
        at=0; max_k=0
        for(b=1;b<=B;b++) {
            nn=p+1; rr=(at+1..at+2*nn)'; at=at+2*nn
            fit=cdid_pooled_wls_fit((p ? pseudo_X[rr,] : J(2*nn,0,0)),pseudo_y[rr],pseudo_w[rr],
                J(2*nn,1,1),J(2*nn,1,1),1,intercept,J(0,p,0),J(0,1,0),rank_tol,1)
            out.blocks[b].beta=fit.beta; out.blocks[b].R=fit.R; out.blocks[b].pivot=fit.pivot
            out.blocks[b].kept=fit.kept; out.blocks[b].dropped=fit.dropped
            out.blocks[b].alias_map=fit.alias_map; out.blocks[b].alias_constants=fit.alias_constants
            out.blocks[b].mean_D=(intercept ? fit.mean_X[,fit.kept] : J(1,rows(fit.beta),0))
            out.alpha[b]=fit.alpha[1]; out.block_mass[b]=fit.mass[1]
            if(rows(fit.beta)>max_k) max_k=rows(fit.beta)
        }
        _cdid_or_alias_rows(out)
        out.alpha_rate=J(M,1,0); out.n_fit_rows=0; out.objective=0
        for(b=1;b<=B;b++) {
            pi=block_pair[b]; r=pairs[pi].recipe; fr=recipes[r].frame; L=pairs[pi].t_index-pairs[pi].s_index
            kb=rows(out.blocks[b].beta); n_ctrl_b=block_start[b+1]-block_start[b]
            force_beta=J(kb,n_ctrl_b,0); sum_sa=0; sum_sb=J(1,kb,0)
            for(m=block_start[b];m<block_start[b+1];m++) {
                sid=block_support[m]; ix=_cdid_or_index(fr.group[pairs[pi].a,1]:==support[sid,2])
                ai=pairs[pi].a[ix,1]; fw=fr.weight[ai,1]/out.rates_c[sid]
                D=L:*recipes[r].Z[ai,out.blocks[b].kept]
                raw_res=pairs[pi].dy[ix,1]:-out.alpha[b]:-(kb ? D*out.blocks[b].beta : 0)
                centered=(intercept ? D:-out.blocks[b].mean_D : D)
                coef_rows=quadcross(fw:*raw_res,centered)
                if(intercept) sum_sa=sum_sa+sum(fw:*raw_res)
                if(kb) sum_sb=sum_sb+coef_rows
                if(out.rates_active[sid]) {
                    if(intercept) out.alpha_rate[m]=-sum(fw:*raw_res)/out.rates_c[sid]
                    if(kb) force_beta[,m-block_start[b]+1]=force_beta[,m-block_start[b]+1]-coef_rows'/out.rates_c[sid]
                }
                out.n_fit_rows=out.n_fit_rows+rows(ai); out.objective=out.objective+sum(fw:*raw_res:^2)/out.n_F
            }
            out.blocks[b].score_mean_alpha=(intercept ? sum_sa/out.n_F : 0)
            out.blocks[b].score_mean_beta=(kb ? sum_sb/out.n_F : J(1,0,0))
            out.blocks[b].beta_rate_loading=_cdid_nuisance_gram(out.blocks[b].R,out.blocks[b].pivot,force_beta)
        }
        out.link_mu=J(Jl,1,0); out.link_c=J(Jl,max_k,0); out.link_d=J(Jl,1,0)
        out.link_numerator=J(Jl,1,0); out.link_rate_loading=J(Jl,1,0)
        for(j=1;j<=Jl;j++) {
            b=links[j,1]; sid=links[j,2]; pi=support[sid,1]; L=pairs[pi].t_index-pairs[pi].s_index
            mass=out.pair_mass[sid]; kb=rows(out.blocks[b].beta)
            dres=out.sumY[sid]-out.alpha[b]*mass-(kb ? L*(out.sumZ[sid,out.blocks[b].kept]*out.blocks[b].beta) : 0)
            out.link_mu[j]=dres/mass
            if(kb) out.link_c[j,1..kb]=L:*out.sumZ[sid,out.blocks[b].kept]/mass
            out.link_d[j]=mass/(out.n_F*out.rates_c[sid]); out.link_numerator[j]=out.link_mu[j]*out.link_d[j]
            if(out.rates_active[sid]) out.link_rate_loading[j]=-(dres-out.link_mu[j]*mass)/
                (out.n_F*out.rates_c[sid]^2*out.link_d[j])
        }
        return(out)
    }
    out.blocks=cdid_or_block_fit(1,0)
    fit=cdid_pooled_wls_fit(pseudo_X,pseudo_y,pseudo_w,pseudo_block,J(nPseudo,1,1),1,
        intercept,J(0,p,0),J(0,1,0),rank_tol,1)
    out.beta=fit.beta; out.alpha=fit.alpha; out.R=fit.R; out.pivot=fit.pivot
    out.kept=fit.kept; out.dropped=fit.dropped; out.alias_map=fit.alias_map; out.alias_constants=fit.alias_constants
    out.mean_D=(intercept ? fit.mean_X[,fit.kept] : J(B,rows(fit.beta),0)); out.block_mass=fit.mass
    _cdid_or_alias_rows(out)
    k=rows(out.beta); out.beta_score_A=J(out.n_A,k,0); force_beta=J(k,S,0); out.alpha_rate=J(M,1,0)
    out.n_fit_rows=0; out.objective=0
    for(b=1;b<=B;b++) {
        pi=block_pair[b]; r=pairs[pi].recipe; fr=recipes[r].frame; L=pairs[pi].t_index-pairs[pi].s_index
        for(m=block_start[b];m<block_start[b+1];m++) {
            sid=block_support[m]; ix=_cdid_or_index(fr.group[pairs[pi].a,1]:==support[sid,2])
            ai=pairs[pi].a[ix,1]; fw=fr.weight[ai,1]/out.rates_c[sid]
            D=L:*recipes[r].Z[ai,out.kept]
            raw_res=pairs[pi].dy[ix,1]:-out.alpha[b]:-D*out.beta
            centered=(intercept ? D:-out.mean_D[b,] : D)
            coef_rows=quadcross(fw:*raw_res,centered)
            if(k) out.beta_score_A[ai,]=out.beta_score_A[ai,]+centered:*(fw:*raw_res)
            if(out.rates_active[sid]) {
                if(intercept) out.alpha_rate[m]=-sum(fw:*raw_res)/out.rates_c[sid]
                if(k) force_beta[,sid]=force_beta[,sid]-coef_rows'/out.rates_c[sid]
            }
            out.n_fit_rows=out.n_fit_rows+rows(ai); out.objective=out.objective+sum(fw:*raw_res:^2)/out.n_F
        }
    }
    out.beta_rate_loading=_cdid_nuisance_gram(out.R,out.pivot,force_beta)
    out.link_mu=J(Jl,1,0); out.link_c=J(Jl,k,0); out.link_d=J(Jl,1,0)
    out.link_numerator=J(Jl,1,0); out.link_rate_loading=J(Jl,1,0)
    for(j=1;j<=Jl;j++) {
        b=links[j,1]; sid=links[j,2]; pi=support[sid,1]; L=pairs[pi].t_index-pairs[pi].s_index
        mass=out.pair_mass[sid]
        dres=out.sumY[sid]-out.alpha[b]*mass-L*(out.sumZ[sid,out.kept]*out.beta)
        out.link_mu[j]=dres/mass; out.link_c[j,]=L:*out.sumZ[sid,out.kept]/mass
        out.link_d[j]=mass/(out.n_F*out.rates_c[sid]); out.link_numerator[j]=out.link_mu[j]*out.link_d[j]
        if(out.rates_active[sid]) out.link_rate_loading[j]=-(dres-out.link_mu[j]*mass)/
            (out.n_F*out.rates_c[sid]^2*out.link_d[j])
    }
    return(out)
}
/* Verify exact aliases on ACTUAL fitting and required prediction rows, not
   merely on compressed rows. Each (block,support) is checked once. */
void _cdid_or_alias_rows(struct cdid_or_backend scalar model)
{
    real scalar b, s, pi, r, j, L, den, err, threshold
    real colvector used, ai, ix, v, scale, cm
    real rowvector column_max, dropped, kept
    real matrix Z
    threshold=64*epsilon(1)*max((1,cols(model.recipes[1].Z)+model.intercept))
    if(model.outcome_pooling=="separate") {
        for(b=1;b<=model.n_blocks;b++) {
            dropped=model.blocks[b].dropped; kept=model.blocks[b].kept
            if(cols(dropped)==0) continue
            used=J(model.n_supports,1,0)
            used[model.block_support[model.block_start[b]..model.block_start[b+1]-1]]=J(model.block_start[b+1]-model.block_start[b],1,1)
            ix=_cdid_or_index(model.links[,1]:==b)
            if(rows(ix)) used[model.links[ix,2]]=J(rows(ix),1,1)
            column_max=J(1,cols(dropped),0)
            for(s=1;s<=model.n_supports;s++) if(used[s]) {
                pi=model.support[s,1]; r=model.pairs[pi].recipe
                ai=model.pairs[pi].a
                ix=_cdid_or_index(model.recipes[r].frame.group[ai,1]:==model.support[s,2]); ai=ai[ix,1]
                if(rows(ai)) column_max=colmax(column_max\abs(model.recipes[r].Z[ai,dropped]))
            }
            for(s=1;s<=model.n_supports;s++) if(used[s]) {
                pi=model.support[s,1]; r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
                ai=model.pairs[pi].a; ix=_cdid_or_index(model.recipes[r].frame.group[ai,1]:==model.support[s,2]); ai=ai[ix,1]
                Z=model.recipes[r].Z[ai,]
                for(j=1;j<=cols(dropped);j++) {
                    cm=(model.intercept ? model.blocks[b].alias_constants[1,j]/L : 0)
                    v=Z[,dropped[j]]-(cols(kept) ? Z[,kept]*model.blocks[b].alias_map[,j] : J(rows(Z),1,0)):-cm
                    scale=abs(Z[,dropped[j]])+(cols(kept) ? abs(Z[,kept])*abs(model.blocks[b].alias_map[,j]) : J(rows(Z),1,0)):+abs(cm)
                    den=max((column_max[j],2.2250738585072014e-308))
                    err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                    if(err>threshold | err>=.) _cdid_or_fail("Outcome alias does not reconstruct actual control/required prediction support.",459)
                }
            }
        }
        return
    }
    if(cols(model.dropped)==0) return
    /* The source guard uses maximum original-row column magnitudes. */
    column_max=J(1,cols(model.dropped),0)
    for(s=1;s<=model.n_supports;s++) {
        pi=model.support[s,1]; r=model.pairs[pi].recipe
        ai=model.pairs[pi].a
        ix=_cdid_or_index(model.recipes[r].frame.group[ai,1]:==model.support[s,2]); ai=ai[ix,1]
        if(rows(ai)) column_max=colmax(column_max\abs(model.recipes[r].Z[ai,model.dropped]))
    }
    for(b=1;b<=model.n_blocks;b++) {
        used=J(model.n_supports,1,0)
        used[model.block_support[model.block_start[b]..model.block_start[b+1]-1]]=J(model.block_start[b+1]-model.block_start[b],1,1)
        ix=_cdid_or_index(model.links[,1]:==b)
        if(rows(ix)) used[model.links[ix,2]]=J(rows(ix),1,1)
        for(s=1;s<=model.n_supports;s++) if(used[s]) {
            pi=model.support[s,1]; r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
            ai=model.pairs[pi].a; ix=_cdid_or_index(model.recipes[r].frame.group[ai,1]:==model.support[s,2]); ai=ai[ix,1]
            Z=model.recipes[r].Z[ai,]
            for(j=1;j<=cols(model.dropped);j++) {
                cm=(model.intercept ? model.alias_constants[b,j]/L : 0)
                v=Z[,model.dropped[j]]-Z[,model.kept]*model.alias_map[,j]:-cm
                scale=abs(Z[,model.dropped[j]])+abs(Z[,model.kept])*abs(model.alias_map[,j]):+abs(cm)
                den=max((column_max[j],2.2250738585072014e-308))
                err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                if(err>threshold | err>=.) _cdid_or_fail("Outcome alias does not reconstruct actual control/required prediction support.",459)
            }
        }
    }
}
real colvector cdid_or_rate_if_chunk(struct cdid_or_backend scalar model,
    real scalar sid, real scalar first_full, real scalar last_full)
{
    real colvector out, ai_all, ix_a, ai, ix_g, row_g, ai_g, pos
    real scalar pi, r, coh, gmass
    if(sid<1 | sid>model.n_supports | sid!=floor(sid) | first_full<1 | last_full>model.n_F |
        first_full>last_full | first_full!=floor(first_full) | last_full!=floor(last_full))
        _cdid_or_fail("Invalid rate support or original-F row chunk.",498)
    pi=model.support[sid,1]; r=model.pairs[pi].recipe; coh=model.support[sid,2]
    out=J(last_full-first_full+1,1,0)
    ai_all=model.recipes[r].frame.unit_to_a[(first_full..last_full)']
    ix_a=_cdid_or_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]
        ix_g=_cdid_or_index(model.recipes[r].frame.group[ai]:==coh)
        if(rows(ix_g)) {
            row_g=ix_a[ix_g]; ai_g=ai[ix_g]
            pos=_cdid_or_pair_pos_vec(model.pairs[pi],ai_g,model.n_A)
            gmass=model.recipes[r].frame.group_mass[coh]
            out[row_g]=model.n_F:*model.recipes[r].frame.weight[ai_g]:*((pos:>0):-model.rates[sid]):/gmass
        }
    }
    return(out)
}
real matrix _cdid_or_blk_beta_if_chunk(struct cdid_or_backend scalar model,
    real scalar b, real scalar first_full, real scalar last_full)
{
    real scalar n, kb, pi, r, L, m, sid, coh
    real matrix score, out, D, centered
    real colvector up, ai_all, ix_a, ai, pos_all, ix_p, row_p, ai_p, pos_p, ix_s, row_s, ai_s, pos_s, fw, res
    n=last_full-first_full+1; kb=rows(model.blocks[b].beta)
    if(kb==0) return(J(n,0,0))
    score=J(n,kb,0)
    pi=model.block_pair[b]; r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
    ai_all=model.recipes[r].frame.unit_to_a[(first_full..last_full)']
    ix_a=_cdid_or_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]
        pos_all=_cdid_or_pair_pos_vec(model.pairs[pi],ai,model.n_A)
        ix_p=_cdid_or_index(pos_all:>0)
        if(rows(ix_p)) {
            row_p=ix_a[ix_p]; ai_p=ai[ix_p]; pos_p=pos_all[ix_p]
            for(m=model.block_start[b];m<model.block_start[b+1];m++) {
                sid=model.block_support[m]; coh=model.support[sid,2]
                ix_s=_cdid_or_index(model.recipes[r].frame.group[ai_p]:==coh)
                if(rows(ix_s)) {
                    row_s=row_p[ix_s]; ai_s=ai_p[ix_s]; pos_s=pos_p[ix_s]
                    fw=(model.n_F/model.rates_c[sid]):*model.recipes[r].frame.weight[ai_s]
                    D=L:*model.recipes[r].Z[ai_s,model.blocks[b].kept]
                    res=model.pairs[pi].dy[pos_s]:-model.alpha[b]:-D*model.blocks[b].beta
                    centered=(model.intercept ? D:-model.blocks[b].mean_D : D)
                    score[row_s,]=score[row_s,]+centered:*(fw:*res)
                }
            }
        }
    }
    score=score:-(model.n_F:*model.blocks[b].score_mean_beta)
    out=_cdid_nuisance_gram(model.blocks[b].R,model.blocks[b].pivot,score')'
    for(m=model.block_start[b];m<model.block_start[b+1];m++) {
        if(!any(model.blocks[b].beta_rate_loading[,m-model.block_start[b]+1])) continue
        sid=model.block_support[m]
        up=cdid_or_rate_if_chunk(model,sid,first_full,last_full)
        out=out+up*model.blocks[b].beta_rate_loading[,m-model.block_start[b]+1]'
    }
    return(out)
}
struct cdid_or_if_state scalar cdid_or_if_state_chunk(struct cdid_or_backend scalar model,
    real scalar first_full, real scalar last_full)
{
    struct cdid_or_if_state scalar out
    real scalar s, n
    real matrix score
    real colvector up, ai_all, ix_a
    if(first_full<1 | last_full>model.n_F | first_full>last_full |
        first_full!=floor(first_full) | last_full!=floor(last_full))
        _cdid_or_fail("Invalid original-F IF chunk.",498)
    n=last_full-first_full+1
    out.outcome_pooling=model.outcome_pooling
    out.first_full=first_full; out.last_full=last_full; out.n_A=model.n_A; out.n_F=model.n_F
    out.n_supports=model.n_supports; out.n_blocks=model.n_blocks; out.fitted_beta=model.beta
    if(model.outcome_pooling=="separate") {
        out.beta=J(n,0,0)
        return(out)
    }
    score=J(n,rows(model.beta),0)
    if(rows(model.beta)) {
        ai_all=model.recipes[1].frame.unit_to_a[(first_full..last_full)']
        ix_a=_cdid_or_index(ai_all:>0)
        if(rows(ix_a)) score[ix_a,]=model.n_F:*model.beta_score_A[ai_all[ix_a],]
    }
    out.beta=_cdid_nuisance_gram(model.R,model.pivot,score')'
    if(rows(model.beta)) {
        for(s=1;s<=model.n_supports;s++) {
            if(!any(model.beta_rate_loading[,s])) continue
            up=cdid_or_rate_if_chunk(model,s,first_full,last_full)
            out.beta=out.beta+up*model.beta_rate_loading[,s]'
        }
    }
    return(out)
}
void _cdid_or_check_state(struct cdid_or_backend scalar model, struct cdid_or_if_state scalar state)
{
    if(state.outcome_pooling!=model.outcome_pooling | state.n_A!=model.n_A | state.n_F!=model.n_F |
        state.n_supports!=model.n_supports | state.n_blocks!=model.n_blocks |
        any(state.fitted_beta:!=model.beta) |
        rows(state.beta)!=state.last_full-state.first_full+1 | cols(state.beta)!=rows(model.beta))
        _cdid_or_fail("IF state must come from the same immutable prepared backend.",498)
}
real colvector _cdid_or_alpha_if_with_beta(struct cdid_or_backend scalar model, real scalar b,
    struct cdid_or_if_state scalar state, real matrix betaIF)
{
    real colvector out, up, beta, ai_all, ix_a, ai, pos_all, ix_p, row_p, ai_p, pos_p, ix_s, row_s, ai_s, pos_s, res
    real rowvector kept, mean_D
    real scalar m, sid, pi, r, L, kb, sep, coh
    out=J(rows(state.beta),1,0)
    if(!model.intercept) return(out)
    sep=(model.outcome_pooling=="separate")
    kept=(sep ? model.blocks[b].kept : model.kept)
    beta=(sep ? model.blocks[b].beta : model.beta)
    mean_D=(sep ? model.blocks[b].mean_D : model.mean_D[b,])
    kb=rows(beta)
    pi=model.block_pair[b]; r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
    ai_all=model.recipes[r].frame.unit_to_a[(state.first_full..state.last_full)']
    ix_a=_cdid_or_index(ai_all:>0)
    ai=(rows(ix_a) ? ai_all[ix_a] : J(0,1,0))
    pos_all=(rows(ai) ? _cdid_or_pair_pos_vec(model.pairs[pi],ai,model.n_A) : J(0,1,0))
    ix_p=(rows(pos_all) ? _cdid_or_index(pos_all:>0) : J(0,1,0))
    row_p=(rows(ix_p) ? ix_a[ix_p] : J(0,1,0))
    ai_p=(rows(ix_p) ? ai[ix_p] : J(0,1,0))
    pos_p=(rows(ix_p) ? pos_all[ix_p] : J(0,1,0))
    for(m=model.block_start[b];m<model.block_start[b+1];m++) {
        sid=model.block_support[m]; coh=model.support[sid,2]
        if(rows(ix_p)) {
            ix_s=_cdid_or_index(model.recipes[r].frame.group[ai_p]:==coh)
            if(rows(ix_s)) {
                row_s=row_p[ix_s]; ai_s=ai_p[ix_s]; pos_s=pos_p[ix_s]
                res=model.pairs[pi].dy[pos_s]:-model.alpha[b]:-(kb ? L:*(model.recipes[r].Z[ai_s,kept]*beta) : 0)
                out[row_s]=out[row_s]+(model.n_F/(model.rates_c[sid]*model.block_mass[b])):*model.recipes[r].frame.weight[ai_s]:*res
            }
        }
        if(model.alpha_rate[m]!=0) {
            up=cdid_or_rate_if_chunk(model,sid,state.first_full,state.last_full)
            out=out+up*(model.alpha_rate[m]/model.block_mass[b])
        }
    }
    if(sep) out=out:-(model.n_F*model.blocks[b].score_mean_alpha/model.block_mass[b])
    if(kb) out=out-betaIF*mean_D'
    return(out)
}
real colvector cdid_or_alpha_if_chunk(struct cdid_or_backend scalar model, real scalar b,
    struct cdid_or_if_state scalar state)
{
    real matrix betaIF
    _cdid_or_check_state(model,state)
    if(b<1 | b>model.n_blocks | b!=floor(b)) _cdid_or_fail("Unknown OR intercept block.",498)
    if(!model.intercept) return(J(rows(state.beta),1,0))
    betaIF=(model.outcome_pooling=="separate" ?
        _cdid_or_blk_beta_if_chunk(model,b,state.first_full,state.last_full) : state.beta)
    return(_cdid_or_alpha_if_with_beta(model,b,state,betaIF))
}
struct cdid_or_if_chunk scalar cdid_or_link_if_chunk(struct cdid_or_backend scalar model,
    real scalar ell, struct cdid_or_if_state scalar state)
{
    struct cdid_or_if_chunk scalar out
    real scalar b, sid, pi, r, L, kb, sep, coh
    real rowvector kept
    real colvector beta, ai_all, ix_a, ai, ix_g, row_g, ai_g, pos, ix_p, row_p, ai_p, pos_p, res
    real matrix betaIF
    _cdid_or_check_state(model,state)
    if(ell<1 | ell>model.n_links | ell!=floor(ell)) _cdid_or_fail("Unknown OR consumer.",498)
    b=model.links[ell,1]; sid=model.links[ell,2]; pi=model.support[sid,1]; coh=model.support[sid,2]
    r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
    sep=(model.outcome_pooling=="separate")
    kept=(sep ? model.blocks[b].kept : model.kept)
    beta=(sep ? model.blocks[b].beta : model.beta)
    kb=rows(beta)
    betaIF=(sep ? _cdid_or_blk_beta_if_chunk(model,b,state.first_full,state.last_full) : state.beta)
    out.full_index=(state.first_full..state.last_full)'; out.components=J(rows(state.beta),4,0)
    ai_all=model.recipes[r].frame.unit_to_a[out.full_index]
    ix_a=_cdid_or_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]
        ix_g=_cdid_or_index(model.recipes[r].frame.group[ai]:==coh)
        if(rows(ix_g)) {
            row_g=ix_a[ix_g]; ai_g=ai[ix_g]
            pos=_cdid_or_pair_pos_vec(model.pairs[pi],ai_g,model.n_A)
            ix_p=_cdid_or_index(pos:>0)
            if(rows(ix_p)) {
                row_p=row_g[ix_p]; ai_p=ai_g[ix_p]; pos_p=pos[ix_p]
                res=model.pairs[pi].dy[pos_p]:-model.alpha[b]:-(kb ? L:*(model.recipes[r].Z[ai_p,kept]*beta) : 0):-model.link_mu[ell]
                out.components[row_p,1]=(model.n_F/model.pair_mass[sid]):*model.recipes[r].frame.weight[ai_p]:*res
            }
        }
    }
    out.components[,2]=-_cdid_or_alpha_if_with_beta(model,b,state,betaIF)
    if(kb) out.components[,3]=-betaIF*model.link_c[ell,1..kb]'
    if(model.link_rate_loading[ell]!=0) {
        out.components[,4]=cdid_or_rate_if_chunk(model,sid,state.first_full,state.last_full)*model.link_rate_loading[ell]
    }
    out.phi=rowsum(out.components)
    if(missing(out.phi)) _cdid_or_fail("Nonfinite original-unit OR influence function.",504)
    return(out)
}
/* The caller owns and reuses the outcome-free treatment nuisance. No logit
   fit occurs here. Ordered target/control groups must match this OR consumer. */
struct cdid_dr_link scalar cdid_dr_link_prepare(struct cdid_or_backend scalar model,
    real scalar ell, struct cdid_ipw_nuisance scalar nu)
{
    struct cdid_dr_link scalar out
    struct cdid_ipw_link scalar raw
    real scalar b, sid, pi, r, L, k, kb, sep
    real colvector co, sids, expected, residual_y, ix, fw, beta
    real rowvector kept
    real matrix D
    if(ell<1 | ell>model.n_links | ell!=floor(ell)) _cdid_or_fail("Unknown DR consumer.",498)
    b=model.links[ell,1]; sid=model.links[ell,2]; pi=model.support[sid,1]; r=model.pairs[pi].recipe
    L=model.pairs[pi].t_index-model.pairs[pi].s_index
    sids=sid\model.block_support[model.block_start[b]..model.block_start[b+1]-1]
    co=model.support[sids,2]
    if(nu.n_A!=model.n_A | nu.n_F!=model.n_F | nu.n_cohorts!=rows(sids) |
        nu.n_pair_all!=rows(model.pairs[pi].a) | any(nu.rates:!=model.rates[sids,1]))
        _cdid_or_fail("Prepared treatment nuisance has a different support/design block.",498)
    expected=J(rows(nu.risk_a),1,0)
    for(k=1;k<=rows(co);k++) {
        ix=_cdid_or_index(model.recipes[r].frame.group[nu.risk_a,1]:==co[k]); expected[ix]=J(rows(ix),1,k)
    }
    if(any(expected:!=nu.risk_group) | any(model.pairs[pi].a[nu.in_risk,1]:!=nu.risk_a) |
        any(nu.risk_weight:!=model.recipes[r].frame.weight[nu.risk_a,1]))
        _cdid_or_fail("Treatment nuisance ordered rows/weights do not match the OR design block.",498)
    sep=(model.outcome_pooling=="separate")
    kept=(sep ? model.blocks[b].kept : model.kept)
    beta=(sep ? model.blocks[b].beta : model.beta)
    kb=rows(beta)
    residual_y=model.pairs[pi].dy
    D=L:*model.recipes[r].Z[nu.risk_a,kept]
    residual_y[nu.in_risk]=residual_y[nu.in_risk,1]:-model.alpha[b]:-(kb ? D*beta : 0)
    out.residual_ipw=cdid_ipw_link_outcome(nu,residual_y)
    raw=cdid_ipw_link_outcome(nu,model.pairs[pi].dy)
    ix=_cdid_or_index(nu.risk_group:==1)
    out.d_bar_plus=(kb ? quadcross(nu.fit_weight[ix,1],D[ix,])'/nu.sum_a : J(0,1,0))
    ix=_cdid_or_index(nu.risk_group:>1)
    fw=nu.fit_weight[ix,1]:*nu.odds[ix,1]
    out.d_bar_minus=(kb ? quadcross(fw,D[ix,])'/nu.sum_b : J(0,1,0))
    out.d_bar_diff=out.d_bar_plus-out.d_bar_minus
    out.link_id=ell; out.estimate=out.residual_ipw.estimate; out.estimate_or=model.link_mu[ell]
    out.estimate_ipw=raw.estimate; out.denominator_plus=nu.sum_a/model.n_F; out.denominator_minus=nu.sum_b/model.n_F
    return(out)
}
/* Components: treated direct, control direct, propensity (raw score plus ALL
   rate-to-propensity channels), direct rate, shared outcome-beta. Alpha
   cancels from the two normalized residual means; it is not omitted in OR. */
struct cdid_or_if_chunk scalar cdid_dr_link_if_chunk(struct cdid_or_backend scalar model,
    struct cdid_dr_link scalar link, struct cdid_ipw_frame scalar frame,
    struct cdid_ipw_nuisance scalar nu, struct cdid_or_if_state scalar state)
{
    struct cdid_or_if_chunk scalar out
    real scalar ell, b, sid, pi, r, L, n, factor, kb, sep
    real rowvector kept
    real colvector beta, r2ps, ai_all, ix_a, ai, jg, ix_g, row_g, ai_g, j_g, pos_pair, up
    real colvector ix_r, row_r, ai_r, j_r, pos_r, residual, ix_t, ix_c, pt, pc
    real matrix betaIF
    _cdid_or_check_state(model,state)
    ell=link.link_id; b=model.links[ell,1]; sid=model.links[ell,2]; pi=model.support[sid,1]
    r=model.pairs[pi].recipe; L=model.pairs[pi].t_index-model.pairs[pi].s_index
    if(frame.n_A!=nu.n_A | frame.n_F!=nu.n_F | any(frame.q:!=nu.q) |
        any(link.residual_ipw.beta:!=nu.beta) | any(link.residual_ipw.rates:!=nu.rates))
        _cdid_or_fail("DR treatment frame/nuisance/link do not match.",498)
    sep=(model.outcome_pooling=="separate")
    kept=(sep ? model.blocks[b].kept : model.kept)
    beta=(sep ? model.blocks[b].beta : model.beta)
    kb=rows(beta)
    betaIF=(sep ? _cdid_or_blk_beta_if_chunk(model,b,state.first_full,state.last_full) : state.beta)
    n=rows(state.beta); out.full_index=(state.first_full..state.last_full)'; out.components=J(n,5,0)
    factor=model.n_F/model.n_A
    r2ps=nu.A_beta_r'*link.residual_ipw.loading_beta
    ai_all=frame.unit_to_a[out.full_index]
    ix_a=_cdid_or_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]; jg=frame.group[ai]
        ix_g=_cdid_or_index(jg:>0)
        if(rows(ix_g)) {
            row_g=ix_a[ix_g]; ai_g=ai[ix_g]; j_g=jg[ix_g]
            pos_pair=_cdid_or_pair_pos_vec(model.pairs[pi],ai_g,model.n_A)
            up=factor:*frame.weight[ai_g]:*((pos_pair:>0):-nu.rates[j_g]):/nu.q[j_g]
            out.components[row_g,3]=up:*r2ps[j_g]
            out.components[row_g,4]=up:*link.residual_ipw.G_r[j_g]
            ix_r=_cdid_or_index(pos_pair:>0)
            if(rows(ix_r)) {
                row_r=row_g[ix_r]; ai_r=ai_g[ix_r]; j_r=j_g[ix_r]
                pos_r=_cdid_ipw_risk_pos(nu,ai_r)
                residual=model.pairs[pi].dy[nu.in_risk[pos_r]]:-model.alpha[b]:-(kb ? L:*(model.recipes[r].Z[ai_r,kept]*beta) : 0)
                out.components[row_r,3]=out.components[row_r,3]+factor:*nu.fit_weight[pos_r]:*((j_r:==1):-nu.raw_probability[pos_r]):*(nu.X[pos_r,]*link.residual_ipw.loading_beta)
                ix_t=_cdid_or_index(j_r:==1)
                if(rows(ix_t)) {
                    pt=pos_r[ix_t]
                    out.components[row_r[ix_t],1]=(factor/nu.mean_a):*nu.fit_weight[pt]:*(residual[ix_t]:-link.residual_ipw.treated_mean)
                }
                ix_c=_cdid_or_index(j_r:!=1)
                if(rows(ix_c)) {
                    pc=pos_r[ix_c]
                    out.components[row_r[ix_c],2]=-(factor/nu.mean_b):*nu.fit_weight[pc]:*nu.odds[pc]:*(residual[ix_c]:-link.residual_ipw.control_mean)
                }
            }
        }
    }
    if(kb) out.components[,5]=-betaIF*link.d_bar_diff
    out.phi=rowsum(out.components)
    if(missing(out.phi)) _cdid_or_fail("Nonfinite original-unit DR influence function.",504)
    return(out)
}
end
