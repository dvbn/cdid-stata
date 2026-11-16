*! Common baseline IPW/reg/DR scientific plan; independent of any GMM solver.
*! Pinned R cdid.
*! Load population/catalogue/nuisance/IPW/OR-DR modules first.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_base(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_baseline_plan_state {
    real scalar n_F,n_A,M,Q,n_consumers,n_pairs,n_supports,n_blocks
    real scalar workspace_entries,workspace_limit,work_ops
    real colvector calendar,cohorts,A_to_F,F_to_A,a_cohort,active_cohorts,cohort_to_active,cohort_A_count
    real matrix pair,support,consumer,raw
    real colvector raw_K_ptr,raw_K_index,raw_consumer,raw_supported,raw_pair_counts
    real colvector raw_count_ptr,raw_group_count,post_consumer,aux_consumer,aux_supported
    string colvector aux_reason
    real colvector block_pair,block_ptr,block_support
    real matrix row_AT,pair_cum
    real colvector pair_a,pair_rs,pair_rt
    real colvector raw_mref_skip
}
void _cdid_base_fail(string scalar cls,string scalar msg,real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}
real colvector _cdid_base_index(real matrix x)
{
    return(vec(selectindex(x)))
}
real scalar _cdid_base_cells(real scalar a,real scalar b,real scalar cap,real scalar used)
{
    if(a>=. | b>=. | a<0 | b<0 | a!=floor(a) | b!=floor(b) | used>=. | used<0 |
        (a>0 & b>(2^53-1)/a))_cdid_base_fail("bound","Invalid/overflowing numeric workspace dimensions.",908)
    if(a*b>cap-used)_cdid_base_fail("workspace","Required numeric workspace exceeds explicit max_entries.",908)
    return(used+a*b)
}
real scalar _cdid_base_work(real scalar x,real scalar cap,real scalar used)
{
    if(x>=. | x<0 | x>(2^53-1)-used | (cap<. & x>cap-used))
        _cdid_base_fail("work","Required work exceeds explicit work cap or exact range.",908)
    return(used+x)
}
real scalar _cdid_base_find2(real matrix x,real scalar a,real scalar b)
{
    real scalar lo,hi,m
    lo=1;hi=rows(x)
    while(lo<=hi) {
        m=floor((lo+hi)/2)
        if(x[m,1]==a & x[m,2]==b)return(m)
        if(x[m,1]<a | (x[m,1]==a & x[m,2]<b))lo=m+1
        else hi=m-1
    }
    return(0)
}
/* One actual pair at a time. Only original row addresses and missing-Y
   membership are read here; changes are formed by the fit reader later. */
real matrix _cdid_base_pair_rows(struct cdid_population scalar P,real scalar s,real scalar t,
    real scalar alo,real scalar ahi)
{
    real colvector ai,rs,rt,ix
    if(ahi<alo)return(J(0,3,0))
    ai=(alo..ahi)'
    rs=cdid_population_rows(P,P.A_to_F[ai],J(rows(ai),1,P.calendar[s]))
    rt=cdid_population_rows(P,P.A_to_F[ai],J(rows(ai),1,P.calendar[t]))
    ix=_cdid_base_index((rs:>0):&(rt:>0))
    if(rows(ix)==0)return(J(0,3,0))
    ai=ai[ix];rs=rs[ix];rt=rt[ix]
    ix=_cdid_base_index((P.y[rs]:<.):&(P.y[rt]:<.))
    if(rows(ix)==0)return(J(0,3,0))
    return((ai[ix],rs[ix],rt[ix]))
}
real matrix _cdid_base_build_row_AT(struct cdid_population scalar P)
{
    real matrix row_AT
    real colvector aa,ix,d_sub,hit,ti_vec
    real scalar T,ti,k,j
    T=rows(P.calendar)
    if(P.n_A<=0 | T<=0 | P.n_A*T>8000000) return(J(0,0,0))
    if(rows(P.obs_grid)==P.N & cols(P.obs_grid)==T) return(P.obs_grid[P.A_to_F,.])
    row_AT=J(P.n_A,T,0)
    if(P.R<=0) return(row_AT)
    aa=P.F_to_A[P.ids.row_unit]
    ix=_cdid_base_index(aa:>0)
    if(rows(ix)==0) return(row_AT)
    aa=aa[ix];d_sub=P.date[ix]
    if(T<=64) {
        for(ti=1;ti<=T;ti++) {
            hit=_cdid_base_index(d_sub:==P.calendar[ti])
            if(rows(hit))row_AT[aa[hit],ti]=ix[hit]
        }
    }
    else {
        ti_vec=_cdid_population_lookup(d_sub,P.calendar)
        hit=_cdid_base_index(ti_vec:>0)
        for(k=1;k<=rows(hit);k++) {
            j=hit[k];row_AT[aa[j],ti_vec[j]]=ix[j]
        }
    }
    return(row_AT)
}
void _cdid_base_build_pair_cache(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar p)
{
    real scalar pi,s,t,tot,n_pi
    real colvector rs,rt,ix,ix2,ai,has_a,yobs
    real matrix ok_AT,rAT
    p.row_AT=_cdid_base_build_row_AT(P)
    p.pair_cum=J(0,0,0);p.pair_a=J(0,1,0);p.pair_rs=J(0,1,0);p.pair_rt=J(0,1,0)
    if(rows(p.row_AT)!=P.n_A | p.n_pairs<=0 | P.n_A<=0 | P.n_A*p.n_pairs>4000000) return
    p.pair_cum=J(P.n_A+1,p.n_pairs,0)
    yobs=J(rows(P.y)+1,1,0)
    yobs[|1\rows(P.y)|]=(P.y:<.)
    /* Row 0 (no actual row) maps to the trailing zero, so one product marks pair support. */
    rAT=editvalue(p.row_AT,0,rows(P.y)+1)
    ok_AT=colshape(yobs[vec(rAT)],P.n_A)'
    tot=0
    for(pi=1;pi<=p.n_pairs;pi++) {
        has_a=ok_AT[,p.pair[pi,1]]:*ok_AT[,p.pair[pi,2]]
        p.pair_cum[,pi]=tot:+(0\runningsum(has_a))
        tot=p.pair_cum[P.n_A+1,pi]
    }
    if(tot<=0) return
    p.pair_a=J(tot,1,0);p.pair_rs=J(tot,1,0);p.pair_rt=J(tot,1,0)
    for(pi=1;pi<=p.n_pairs;pi++) {
        n_pi=p.pair_cum[P.n_A+1,pi]-p.pair_cum[1,pi]
        if(n_pi<=0)continue
        s=p.pair[pi,1];t=p.pair[pi,2]
        ai=_cdid_base_index(ok_AT[,s]:*ok_AT[,t])
        tot=p.pair_cum[1,pi]
        p.pair_a[|tot+1\tot+n_pi|]=ai
        p.pair_rs[|tot+1\tot+n_pi|]=p.row_AT[ai,s]
        p.pair_rt[|tot+1\tot+n_pi|]=p.row_AT[ai,t]
    }
}
real matrix _cdid_base_pair_slice(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar p,
    real scalar pi,real scalar alo,real scalar ahi)
{
    real scalar klo,khi,s,t
    real colvector ai,rs,rt,ix
    if(ahi<alo)return(J(0,3,0))
    if(rows(p.pair_cum)==p.n_A+1 & cols(p.pair_cum)==p.n_pairs & alo>=1 & ahi<=p.n_A) {
        klo=p.pair_cum[alo,pi]+1;khi=p.pair_cum[ahi+1,pi]
        if(khi<klo)return(J(0,3,0))
        return((p.pair_a[|klo\khi|],p.pair_rs[|klo\khi|],p.pair_rt[|klo\khi|]))
    }
    if(rows(p.row_AT)==p.n_A & alo>=1 & ahi<=p.n_A) {
        s=p.pair[pi,1];t=p.pair[pi,2]
        ai=(alo..ahi)';rs=p.row_AT[ai,s];rt=p.row_AT[ai,t]
        ix=_cdid_base_index((rs:>0):&(rt:>0))
        if(rows(ix)==0)return(J(0,3,0))
        ai=ai[ix];rs=rs[ix];rt=rt[ix]
        ix=_cdid_base_index((P.y[rs]:<.):&(P.y[rt]:<.))
        if(rows(ix)==0)return(J(0,3,0))
        return((ai[ix],rs[ix],rt[ix]))
    }
    return(_cdid_base_pair_rows(P,p.pair[pi,1],p.pair[pi,2],alo,ahi))
}
real scalar _cdid_base_same_K(struct cdid_baseline_plan_state scalar p,real scalar i,real scalar j)
{
    real scalar ni,nj
    ni=p.raw_K_ptr[i+1]-p.raw_K_ptr[i];nj=p.raw_K_ptr[j+1]-p.raw_K_ptr[j]
    if(ni!=nj)return(0)
    return(all(p.raw_K_index[p.raw_K_ptr[i]..p.raw_K_ptr[i+1]-1]:==
        p.raw_K_index[p.raw_K_ptr[j]..p.raw_K_ptr[j+1]-1]))
}
real scalar _cdid_base_K_lex(struct cdid_baseline_plan_state scalar p,real scalar i,real scalar j)
{
    real scalar a,b,ni,nj,k
    ni=p.raw_K_ptr[i+1]-p.raw_K_ptr[i];nj=p.raw_K_ptr[j+1]-p.raw_K_ptr[j]
    for(k=0;k<min((ni,nj));k++) {
        a=p.raw_K_index[p.raw_K_ptr[i]+k];b=p.raw_K_index[p.raw_K_ptr[j]+k]
        if(a!=b)return(sign(a-b))
    }
    return(sign(ni-nj))
}
/* Source pipeline descriptors comprise ALL original post comparisons followed
   by the pure request layer's deduplicated BOTH-family auxiliary union.
   Pair supports are decided from counts; no outcome value selects a block. */
struct cdid_baseline_plan_state scalar cdid_baseline_plan(struct cdid_population scalar P,
    struct cdid_cat scalar C,real matrix aux,real colvector aux_K_ptr,real colvector aux_K_index,
    real scalar cap,real scalar workcap)
{
    struct cdid_baseline_plan_state scalar p
    real scalar R,M,Q,Ktotal,used,work,i,j,k,u,a,g,s,t,op,pool,nK,at,pi,nc,nb,hit,b,old,ordj,S
    real scalar rr,sz,total_support,ci,ai
    real colvector K,rawpair,counts,members,crep,brep,cb,bo,bmap,ix,active,raworder,c_ord
    real matrix pairrows,needed,consumer,c_info,coh_ind,all_counts
    real scalar has_cache,n_pr,mref_on,ectx,sup_e
    real colvector view_ctx,aux_miss,cnt_e,pr_e,lost
    real matrix elig_A,prows_e
    string scalar kpool,lostc
    if(cap>=. | cap<1 | cap!=floor(cap) | cap>2^53-1 | workcap>. |
        (workcap<. & (workcap<1 | workcap!=floor(workcap))))
        _cdid_base_fail("bound","Supply explicit positive numeric-entry cap and missing/positive work cap.",498)
    M=rows(C.r_g);Q=rows(aux);R=M+Q
    if(cols(aux)!=3 | missing(aux) | rows(aux_K_ptr)!=Q+1 | cols(aux_K_ptr)!=1 |
        cols(aux_K_index)!=1 | missing(aux_K_ptr) | missing(aux_K_index) |
        any(aux_K_ptr:!=floor(aux_K_ptr)) | any(aux_K_index:!=floor(aux_K_index)) |
        aux_K_ptr[1]!=1 | aux_K_ptr[Q+1]!=rows(aux_K_index)+1 |
        any(aux_K_index:<1) | any(aux_K_index:>C.n_cohorts))
        _cdid_base_fail("request","Invalid explicit auxiliary descriptor/control CSR.",498)
    if(rows(P.calendar)!=rows(C.calendar) | !all(P.calendar:==C.calendar))
        _cdid_base_fail("dictionary","Population and catalogue calendars disagree.",498)
    Ktotal=rows(aux_K_index)
    for(i=1;i<=M;i++)Ktotal=_cdid_base_cells(1,C.pool_size[C.r_pool[i]],cap,Ktotal)
    used=_cdid_base_cells(96,R,cap,0)
    used=_cdid_base_cells(20,Ktotal,cap,used)
    used=_cdid_base_cells(16,P.N+P.n_A+C.n_cohorts,cap,used)
    p.n_F=P.N;p.n_A=P.n_A;p.M=M;p.Q=Q;p.calendar=C.calendar;p.cohorts=C.cohorts
    p.workspace_limit=cap;p.A_to_F=P.A_to_F;p.F_to_A=P.F_to_A
    p.raw=J(R,3,0);p.raw_K_ptr=J(R+1,1,1);p.raw_K_index=J(Ktotal,1,0)
    p.raw_consumer=J(R,1,0);p.raw_supported=J(R,1,0);p.raw_pair_counts=J(R,1,0)
    p.raw_count_ptr=J(R+1,1,1);p.raw_group_count=J(R+Ktotal,1,0)
    p.post_consumer=J(M,1,0);p.aux_consumer=J(Q,1,0);p.aux_supported=J(Q,1,0)
    p.aux_reason=J(Q,1,"");p.a_cohort=J(P.n_A,1,0);p.cohort_A_count=J(C.n_cohorts,1,0)
    p.cohort_to_active=J(C.n_cohorts,1,0)
    if(P.n_A>0) {
        if(any(P.cohort_known[P.A_to_F]:!=1))_cdid_base_fail("cohort","A-unit cohort must be known.",459)
        p.a_cohort=_cdid_cohort_index_vec(C.cohorts,C.nt_index,P.cohort[P.A_to_F],J(P.n_A,1,1))
        if(any(p.a_cohort:==0))_cdid_base_fail("dictionary","A-unit cohort absent from catalogue dictionary.",498)
        c_ord=sort(p.a_cohort,1)
        c_info=_cdid_ps_base(c_ord,1)
        p.cohort_A_count[c_ord[c_info[,1]]]=c_info[,2]:-c_info[,1]:+1
    }
    p.active_cohorts=_cdid_base_index(p.cohort_A_count:>0)
    if(rows(p.active_cohorts))p.cohort_to_active[p.active_cohorts]=(1..rows(p.active_cohorts))'
    at=0
    for(i=1;i<=R;i++) {
        if(i<=M) {
            pool=C.r_pool[i];nK=C.pool_size[pool]
            K=C.pool_entries[C.pool_start[pool]..C.pool_start[pool]+nK-1]
            p.raw[i,]=(C.r_g[i],C.r_s[i],C.r_t[i])
        }
        else {
            j=i-M;nK=aux_K_ptr[j+1]-aux_K_ptr[j]
            if(nK<1)_cdid_base_fail("request","An auxiliary descriptor requires a nonempty declared control set.",498)
            K=aux_K_index[aux_K_ptr[j]..aux_K_ptr[j+1]-1]
            g=_cdid_cohort_index(C.cohorts,C.nt_index,aux[j,1])
            s=_cdid_bsearch(C.calendar,aux[j,2]);t=_cdid_bsearch(C.calendar,aux[j,3])
            if(g==0 | aux[j,1]>=. | s==0 | t==0 | s>=t)
                _cdid_base_fail("request","Auxiliary g/s/t absent or unordered in shared dictionaries.",498)
            p.raw[i,]=(g,s,t)
        }
        if(nK<1 | any(K:==p.raw[i,1]))_cdid_base_fail("control","Empty control set or target included in controls.",498)
        if(nK>1)if(any(K[2..nK]:<=K[1..nK-1]))_cdid_base_fail("control","Controls must be ordered unique cohort indices.",498)
        p.raw_K_ptr[i]=at+1;p.raw_K_index[at+1..at+nK]=K;at=at+nK
        p.raw_count_ptr[i+1]=p.raw_count_ptr[i]+nK+1
    }
    p.raw_K_ptr[R+1]=at+1
    p.pair=(R ? uniqrows(sort(p.raw[,2..3],(1,2))) : J(0,2,0));p.n_pairs=rows(p.pair)
    _cdid_base_build_pair_cache(P,p)
    rawpair=J(R,1,0)
    for(i=1;i<=R;i++)rawpair[i]=_cdid_base_find2(p.pair,p.raw[i,2],p.raw[i,3])
    work=0
    /* Per-pair cohort counts are exact integer sums over the cached pair membership. */
    has_cache=(rows(p.pair_cum)==P.n_A+1 & cols(p.pair_cum)==p.n_pairs & P.n_A>0 & p.n_pairs>0)
    if(has_cache) {
        coh_ind=J(P.n_A*C.n_cohorts,1,0)
        coh_ind[(1::P.n_A)+(p.a_cohort:-1):*P.n_A]=J(P.n_A,1,1)
        coh_ind=colshape(coh_ind,P.n_A)'
        all_counts=cross(p.pair_cum[|2,1\P.n_A+1,p.n_pairs|]-p.pair_cum[|1,1\P.n_A,p.n_pairs|],coh_ind)
    }
    /* Restricted eligibility contexts (missingreference(exclude)): views of a
       restricted target count only its eligible units; auxiliary comparisons
       without eligible support, or whose own earlier-date covariate is missing,
       are unavailable and demand nothing. */
    p.raw_mref_skip=J(R,1,0)
    mref_on=0;aux_miss=J(Q,1,0)
    if(findexternal("__cdid_mref_state")!=NULL) {
        mref_on=cdid_mrbp_active()
        if(Q) aux_miss=cdid_mrbp_aux_missing(Q)
    }
    view_ctx=J(R,1,0);elig_A=J(P.n_A,0,.)
    if(mref_on & R) {
        view_ctx=cdid_mrbp_ctx_of(C.cohorts[p.raw[,1]])
        elig_A=cdid_mrbp_excluded(P.N)
        elig_A=(cols(elig_A) & P.n_A ? !elig_A[P.A_to_F,] : J(P.n_A,cols(elig_A),1))
    }
    for(pi=1;pi<=p.n_pairs;pi++) {
        if(has_cache) n_pr=p.pair_cum[P.n_A+1,pi]-p.pair_cum[1,pi]
        else {
            pairrows=_cdid_base_pair_slice(P,p,pi,1,P.n_A)
            n_pr=rows(pairrows)
        }
        work=_cdid_base_work(2*P.n_A*(1+ceil(ln(rows(P.calendar)+1)/ln(2))),workcap,work)
        if(has_cache) counts=all_counts[pi,]'
        else {
            counts=J(C.n_cohorts,1,0)
            if(n_pr) {
                c_ord=sort(p.a_cohort[pairrows[,1]],1)
                c_info=_cdid_ps_base(c_ord,1)
                counts[c_ord[c_info[,1]]]=c_info[,2]-c_info[,1]:+1
            }
        }
        members=_cdid_base_index(rawpair:==pi)
        for(j=1;j<=rows(members);j++) {
            i=members[j];g=p.raw[i,1];K=p.raw_K_index[p.raw_K_ptr[i]..p.raw_K_ptr[i+1]-1]
            p.raw_pair_counts[i]=n_pr
            p.raw_group_count[p.raw_count_ptr[i]..p.raw_count_ptr[i+1]-1]=counts[g\K]
            p.raw_supported[i]=all(p.cohort_A_count[g\K]:>0) & all(counts[g\K]:>0)
            if(!p.raw_supported[i]) {
                if(i<=M)_cdid_base_fail("post_support","Original catalogue row lacks declared target/control pair support; none is dropped.",459)
                p.aux_reason[i-M]="unsupported_pair"
            }
            else if(view_ctx[i]>0) {
                ectx=view_ctx[i]
                prows_e=_cdid_base_pair_slice(P,p,pi,1,P.n_A);cnt_e=J(C.n_cohorts,1,0)
                if(rows(prows_e)) {
                    pr_e=select(prows_e[,1],elig_A[prows_e[,1],ectx])
                    for(u=1;u<=rows(pr_e);u++) cnt_e[p.a_cohort[pr_e[u]]]=cnt_e[p.a_cohort[pr_e[u]]]+1
                }
                p.raw_pair_counts[i]=sum(cnt_e[g\K])
                p.raw_group_count[p.raw_count_ptr[i]..p.raw_count_ptr[i+1]-1]=cnt_e[g\K]
                sup_e=all(cnt_e[g\K]:>0)
                if(!sup_e) {
                    if(i<=M)_cdid_base_fail("post_support","A main comparison lost its support inside its eligibility context; the restricted catalogue and the comparisons disagree.",498)
                    lost=g\K;lost=select(lost,cnt_e[lost]:<1);lostc=""
                    for(u=1;u<=rows(lost);u++) lostc=lostc+(u>1 ? ", " : "")+_cdid_base_num(C.cohorts[lost[u]])
                    p.raw_supported[i]=0;p.raw_mref_skip[i]=1
                    p.aux_reason[i-M]="Unavailable in the restricted eligibility context of target "+_cdid_base_num(C.cohorts[g])+
                        " (missingreference(exclude)): no eligible unit of cohort(s) "+lostc+" is observed at both dates of ("+
                        _cdid_base_num(C.calendar[p.raw[i,2]])+", "+_cdid_base_num(C.calendar[p.raw[i,3]])+")"
                }
            }
            if(i>M & p.raw_supported[i]) if(aux_miss[i-M]) {
                p.raw_supported[i]=0;p.raw_mref_skip[i]=1
                p.aux_reason[i-M]=cdid_mrbp_aux_missing_reason(i-M)
            }
            if(i>M & p.raw_mref_skip[i]) {
                kpool=""
                for(u=1;u<=rows(K);u++) kpool=kpool+(u>1 ? "," : "")+_cdid_base_num(C.cohorts[K[u]])
                cdid_mrbp_aux_note(i-M,C.cohorts[g],C.calendar[p.raw[i,2]],C.calendar[p.raw[i,3]],kpool,p.aux_reason[i-M])
                p.aux_reason[i-M]=p.aux_reason[i-M]+" [g="+_cdid_base_num(C.cohorts[g])+", s="+_cdid_base_num(C.calendar[p.raw[i,2]])+
                    ", t="+_cdid_base_num(C.calendar[p.raw[i,3]])+", K={"+kpool+"}]"
            }
        }
    }
    /* Consumer identity is exact(g,pair,ordered controls); functional labels
       never duplicate a scientific control block in the pooled objective. */
    crep=J(R,1,0);brep=J(R,1,0);cb=J(R,1,0);nc=0;nb=0
    for(i=1;i<=R;i++) {
        if(!p.raw_supported[i])continue
        hit=0
        for(j=1;j<=nc;j++) {
            k=crep[j]
            if(p.raw[i,1]==p.raw[k,1] & rawpair[i]==rawpair[k])if(_cdid_base_same_K(p,i,k)){hit=j;break;}
        }
        if(hit==0) {
            nc++;crep[nc]=i;hit=nc;b=0
            for(j=1;j<=nb;j++)if(rawpair[i]==rawpair[brep[j]])if(_cdid_base_same_K(p,i,brep[j])){b=j;break;}
            if(b==0){nb++;brep[nb]=i;b=nb;}
            cb[nc]=b
        }
        p.raw_consumer[i]=hit
    }
    p.n_consumers=nc;p.n_blocks=nb;p.consumer=J(nc,4,0)
    bo=(nb ? (1..nb)' : J(0,1,0))
    for(i=2;i<=nb;i++) {
        old=bo[i];j=i-1
        while(j>=1) {
            k=bo[j]
            if(rawpair[brep[k]]<rawpair[brep[old]])break
            if(rawpair[brep[k]]==rawpair[brep[old]])if(_cdid_base_K_lex(p,brep[k],brep[old])<=0)break
            bo[j+1]=bo[j];j--
        }
        bo[j+1]=old
    }
    bmap=J(nb,1,0);if(nb)bmap[bo]=(1..nb)'
    total_support=0
    for(j=1;j<=nc;j++){i=crep[j];total_support=total_support+1+p.raw_K_ptr[i+1]-p.raw_K_ptr[i];}
    used=_cdid_base_cells(4,total_support,cap,used)
    needed=J(total_support,2,0);at=0
    for(j=1;j<=nc;j++) {
        i=crep[j];K=p.raw[i,1]\p.raw_K_index[p.raw_K_ptr[i]..p.raw_K_ptr[i+1]-1]
        for(k=1;k<=rows(K);k++){at++;needed[at,]=(rawpair[i],p.cohort_to_active[K[k]]);}
    }
    p.support=(total_support ? uniqrows(sort(needed,(1,2))) : J(0,2,0));p.n_supports=rows(p.support)
    p.block_pair=J(nb,1,0);p.block_ptr=J(nb+1,1,1);at=0
    for(b=1;b<=nb;b++){i=brep[bo[b]];p.block_ptr[b+1]=p.block_ptr[b]+p.raw_K_ptr[i+1]-p.raw_K_ptr[i];}
    used=_cdid_base_cells(2,p.block_ptr[nb+1]-1,cap,used)
    p.block_support=J(p.block_ptr[nb+1]-1,1,0)
    for(b=1;b<=nb;b++) {
        i=brep[bo[b]];p.block_pair[b]=rawpair[i]
        K=p.raw_K_index[p.raw_K_ptr[i]..p.raw_K_ptr[i+1]-1]
        for(k=1;k<=rows(K);k++){at++;p.block_support[at]=_cdid_base_find2(p.support,rawpair[i],p.cohort_to_active[K[k]]);}
    }
    for(j=1;j<=nc;j++){i=crep[j];p.consumer[j,]=(p.raw[i,1],rawpair[i],bmap[cb[j]],_cdid_base_find2(p.support,rawpair[i],p.cohort_to_active[p.raw[i,1]]));}
    if(M)p.post_consumer=p.raw_consumer[1..M]
    if(Q){p.aux_consumer=p.raw_consumer[M+1..R];p.aux_supported=p.raw_supported[M+1..R];}
    p.workspace_entries=used;p.work_ops=work
    return(p)
}
struct cdid_baseline_controls {
    string scalar method,outcome_pooling
    real scalar intercept,glm_epsilon,glm_maxit,rank_tol,chunk_rows,max_entries,workcap,max_cached
}
struct cdid_base_support_qr {
    real matrix F
}
struct cdid_base_or_block {
    real rowvector kept,dropped,pivot
    real colvector beta
    real matrix R,alias_map,alias_constants,mean_D,beta_rate_loading
    real rowvector score_mean_beta
    real scalar score_mean_alpha
}
struct cdid_base_or_fit {
    string scalar outcome_pooling
    struct cdid_base_or_block colvector blocks
    real scalar n_F,n_A,n_supports,n_blocks,n_consumers,intercept,p,n_fit_rows,objective,retained_entries,peak_entries
    real colvector pair_mass,sumY,rates,rates_c,rates_active,alpha,beta,block_mass,alpha_rate
    real matrix sumZ,mean_D,R,alias_map,alias_constants,beta_score_A,beta_rate_loading,link_c
    real colvector link_mu,link_d,link_rate_loading
    real rowvector kept,dropped,pivot
}
struct cdid_base_ps {
    real scalar n_A,n_F,n_cohorts,p,iterations,deviance,score_max,mean_a,mean_b
    real colvector cohorts,q,pair_count,pair_mass,rates,rates_c,rates_active,beta
    real matrix A_beta_r
}
struct cdid_base_point {
    real scalar estimate,ps_id,or_id
    struct cdid_ipw_link scalar residual_ipw
    real colvector d_bar_diff
}
struct cdid_base_ps_cache {
    real scalar rep_ell,pi
    real colvector risk_a,risk_pos_a,raw_probability,odds,cmap
}
struct cdid_baseline_provider {
    string scalar method,status,weight_row_domain
    real scalar n_F,n_A,core_to_frame,normalization_count,weight_rows,sample_rank_upper
    real scalar n_ps_fits,n_ps_reuses,n_or_fits,workspace_entries,workspace_limit,chunk_rows
    struct cdid_population scalar population
    struct cdid_baseline_plan_state scalar plan
    real matrix X_A,Z_A
    real colvector w_A,points,group_mass
    struct cdid_base_or_fit scalar outcome
    struct cdid_base_ps colvector propensity
    struct cdid_base_point colvector point
    struct cdid_base_ps_cache colvector ps_cache
}
/* The caller records these source settings explicitly. No environment/read of
   c(level), RNG, estimator default or solver selection belongs here. */
struct cdid_baseline_controls scalar cdid_baseline_control_make(string scalar method,
    real scalar intercept,real scalar epsilon,real scalar maxit,real scalar rank_tol,
    real scalar chunk_rows,real scalar max_entries,real scalar workcap,real scalar max_cached,
    | string scalar outcome_pooling_in)
{
    struct cdid_baseline_controls scalar out
    string scalar outcome_pooling
    outcome_pooling=(args()<10 ? "" : strtrim(outcome_pooling_in))
    if(outcome_pooling=="")outcome_pooling="separate"
    if((method!="ipw" & method!="reg" & method!="dr") | (intercept!=0 & intercept!=1) |
        (outcome_pooling!="separate" & outcome_pooling!="pooled") |
        epsilon>=. | epsilon<=0 | maxit>=. | maxit<1 | maxit!=floor(maxit) |
        rank_tol>=. | rank_tol<=0 | rank_tol>=1 | chunk_rows>=. | chunk_rows<1 |
        chunk_rows!=floor(chunk_rows) | max_entries>=. | max_entries<1 |
        max_entries!=floor(max_entries) | max_entries>2^53-1 | workcap>. |
        (workcap<. & (workcap<1 | workcap!=floor(workcap))) | max_cached>=. |
        max_cached<0 | max_cached!=floor(max_cached))
        _cdid_base_fail("control","Invalid explicit baseline method/solver/resource controls.",498)
    out.method=method;out.outcome_pooling=outcome_pooling
    out.intercept=intercept;out.glm_epsilon=epsilon;out.glm_maxit=maxit
    out.rank_tol=rank_tol;out.chunk_rows=chunk_rows;out.max_entries=max_entries
    out.workcap=workcap;out.max_cached=max_cached
    return(out)
}
/* Native sparse pair reader. At most one requested pair's current row chunk;
   no all-pair observed-index or outcome-change array is retained. */
real matrix _cdid_base_read(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar p,
    real scalar pi,real scalar alo,real scalar ahi,real colvector groups)
{
    real matrix r
    real colvector dy,keep,cmap
    r=_cdid_base_pair_slice(P,p,pi,alo,ahi)
    if(rows(r)==0)return(J(0,2,0))
    if(rows(groups)==1)keep=_cdid_base_index(p.a_cohort[r[,1]]:==groups[1])
    else {
        cmap=J(rows(p.cohorts),1,0);cmap[groups]=J(rows(groups),1,1)
        keep=_cdid_base_index(cmap[p.a_cohort[r[,1]]])
    }
    if(rows(keep)==0)return(J(0,2,0))
    r=r[keep,]
    dy=P.y[r[,3]]-P.y[r[,2]]
    if(missing(dy))_cdid_base_fail("outcome","Pair outcome change overflowed.",504)
    return((r[,1],dy))
}
/* Build transient CSR support buckets across s = 1..n_supports for fast OR fit. */
real scalar _cdid_base_build_sup_csr(struct cdid_population scalar P,
    struct cdid_baseline_plan_state scalar p,
    real colvector sup_lo,real colvector sup_hi,
    real colvector sup_a,real colvector sup_rs,real colvector sup_rt)
{
    real scalar S,n_act,pi,tot,at,m
    real matrix sup_of_pa,rr,s_info
    real colvector sid,keep,all_sid,all_a,all_rs,all_rt,ord
    S=p.n_supports;n_act=rows(p.active_cohorts)
    sup_lo=J(S,1,1);sup_hi=J(S,1,0)
    sup_a=J(0,1,0);sup_rs=J(0,1,0);sup_rt=J(0,1,0)
    if(S<=0 | p.n_pairs<=0 | n_act<=0 | p.n_A*p.n_pairs>1000000) return(0)
    sup_of_pa=J(n_act,p.n_pairs,0)
    for(m=1;m<=S;m++) sup_of_pa[p.support[m,2],p.support[m,1]]=m
    tot=0
    for(pi=1;pi<=p.n_pairs;pi++) {
        rr=_cdid_base_pair_slice(P,p,pi,1,p.n_A)
        tot=tot+rows(rr)
    }
    if(tot<=0) return(1)
    all_sid=J(tot,1,0);all_a=J(tot,1,0);all_rs=J(tot,1,0);all_rt=J(tot,1,0)
    at=0
    for(pi=1;pi<=p.n_pairs;pi++) {
        rr=_cdid_base_pair_slice(P,p,pi,1,p.n_A)
        if(rows(rr)==0)continue
        sid=sup_of_pa[p.cohort_to_active[p.a_cohort[rr[,1]]],pi]
        keep=_cdid_base_index(sid:>0)
        m=rows(keep)
        if(m==0)continue
        all_sid[|at+1\at+m|]=sid[keep]
        all_a[|at+1\at+m|]=rr[keep,1]
        all_rs[|at+1\at+m|]=rr[keep,2]
        all_rt[|at+1\at+m|]=rr[keep,3]
        at=at+m
    }
    if(at<=0) return(1)
    if(at<tot) {
        all_sid=all_sid[|1\at|];all_a=all_a[|1\at|]
        all_rs=all_rs[|1\at|];all_rt=all_rt[|1\at|]
    }
    ord=order((all_sid,all_a),(1,2))
    all_sid=all_sid[ord];sup_a=all_a[ord];sup_rs=all_rs[ord];sup_rt=all_rt[ord]
    s_info=_cdid_ps_base(all_sid,1)
    if(rows(s_info)) {
        sid=all_sid[s_info[,1]]
        sup_lo[sid]=s_info[,1];sup_hi[sid]=s_info[,2]
    }
    return(1)
}
/* Check the alias map on ACTUAL fitting and required prediction
   rows. Neither compressed rows nor a declared alias can waive this check. */
void _cdid_base_or_alias(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar p,
    real matrix Z,struct cdid_base_or_fit scalar f,real scalar chunk)
{
    real scalar b,s,pi,L,j,lo,hi,ai,den,err,threshold,cm
    real colvector used,ix,a,v,scale
    real matrix rr,z
    real rowvector maxima,dropped,kept
    threshold=64*epsilon(1)*max((1,cols(Z)+f.intercept))
    if(f.outcome_pooling=="separate") {
        for(b=1;b<=p.n_blocks;b++) {
            dropped=f.blocks[b].dropped;kept=f.blocks[b].kept
            if(cols(dropped)==0)continue
            used=J(p.n_supports,1,0)
            used[p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1]]=J(p.block_ptr[b+1]-p.block_ptr[b],1,1)
            ix=_cdid_base_index(p.consumer[,3]:==b)
            if(rows(ix))used[p.consumer[ix,4]]=J(rows(ix),1,1)
            maxima=J(1,cols(dropped),0)
            for(s=1;s<=p.n_supports;s++)if(used[s]) {
                pi=p.support[s,1]
                for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                    rr=_cdid_base_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                    if(rows(rr)==0)continue
                    ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                    if(rows(ix))maxima=colmax(maxima\abs(Z[rr[ix,1],dropped]))
                }
            }
            for(s=1;s<=p.n_supports;s++)if(used[s]) {
                pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
                for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                    rr=_cdid_base_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                    if(rows(rr)==0)continue
                    ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                    if(rows(ix)==0)continue
                    z=Z[rr[ix,1],]
                    for(j=1;j<=cols(dropped);j++) {
                        cm=(f.intercept ? f.blocks[b].alias_constants[1,j]/L : 0)
                        v=z[,dropped[j]]-(cols(kept) ? z[,kept]*f.blocks[b].alias_map[,j] : J(rows(z),1,0)):-cm
                        scale=abs(z[,dropped[j]])+(cols(kept) ? abs(z[,kept])*abs(f.blocks[b].alias_map[,j]) : J(rows(z),1,0)):+abs(cm)
                        den=max((maxima[j],2.2250738585072014e-308))
                        err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                        if(err>=. | err>threshold)_cdid_base_fail("alias","Outcome alias does not preserve actual fitting/required prediction support.",459)
                    }
                }
            }
        }
        return
    }
    if(cols(f.dropped)==0)return
    maxima=J(1,cols(f.dropped),0)
    for(s=1;s<=p.n_supports;s++) {
        pi=p.support[s,1]
        for(lo=1;lo<=p.n_A;lo=lo+chunk) {
            rr=_cdid_base_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
            if(rows(rr)==0)continue
            ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
            if(rows(ix))maxima=colmax(maxima\abs(Z[rr[ix,1],f.dropped]))
        }
    }
    for(b=1;b<=p.n_blocks;b++) {
        used=J(p.n_supports,1,0)
        used[p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1]]=J(p.block_ptr[b+1]-p.block_ptr[b],1,1)
        ix=_cdid_base_index(p.consumer[,3]:==b)
        if(rows(ix))used[p.consumer[ix,4]]=J(rows(ix),1,1)
        for(s=1;s<=p.n_supports;s++)if(used[s]) {
            pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
            for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                rr=_cdid_base_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                if(rows(rr)==0)continue
                ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                if(rows(ix)==0)continue
                z=Z[rr[ix,1],]
                for(j=1;j<=cols(f.dropped);j++) {
                    cm=(f.intercept ? f.alias_constants[b,j]/L : 0)
                    v=z[,f.dropped[j]]-z[,f.kept]*f.alias_map[,j]:-cm
                    scale=abs(z[,f.dropped[j]])+abs(z[,f.kept])*abs(f.alias_map[,j]):+abs(cm)
                    den=max((maxima[j],2.2250738585072014e-308))
                    err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                    if(err>=. | err>threshold)_cdid_base_fail("alias","Outcome alias does not preserve actual fitting/required prediction support.",459)
                }
            }
        }
    }
}

/* Same augmented-QR/FWL algebra, with an actual sparse pair reader
   instead of a resident array of every pair's N observations. */
struct cdid_base_or_fit scalar _cdid_base_or_fit(struct cdid_population scalar P,
    struct cdid_baseline_plan_state scalar p,real matrix Z,real colvector w,
    struct cdid_baseline_controls scalar ctrl,real scalar retained_base)
{
    struct cdid_base_or_fit scalar out
    struct cdid_base_support_qr colvector summary
    struct cdid_wls_result scalar fit
    real scalar S,B,Jl,d,s,b,m,j,pi,L,lo,hi,first,nn,mass,k,kb,max_k,at,nPseudo,scale,den,used,sum_sa,n_ctrl_b,has_csr
    real colvector groupmass,ai,ww,yy,rr,ix,fw,rawres,sup_lo,sup_hi,sup_a,sup_rs,sup_rt
    real rowvector constant,baseZ,sum_sb
    real matrix rows0,rawZ,D,F,V,Q,R,A,coefrows,resfactor,xmin,xmax,firstZ,values
    real matrix pseudoX,pseudoY,pseudoW,pseudoB,force,centered
    real rowvector tau
    S=p.n_supports;B=p.n_blocks;Jl=p.n_consumers;d=cols(Z)
    if(S<1 | B<1 | Jl<1 | rows(Z)!=p.n_A | rows(w)!=p.n_A |
        cols(w)!=1 | missing(Z) | missing(w) | any(w:<=0))
        _cdid_base_fail("outcome_design","Require finite A-supported outcome design and positive A weights.",498)
    used=_cdid_base_cells(S,d+2,ctrl.max_entries,0)
    used=_cdid_base_cells(12*used,d+2,ctrl.max_entries,retained_base)
    den=_cdid_base_cells(B,d+2,ctrl.max_entries,0)
    used=_cdid_base_cells(40*den,d+2,ctrl.max_entries,used)
    used=_cdid_base_cells(8*p.n_A,d+2,ctrl.max_entries,used)
    used=_cdid_base_cells(12*S+12*B,Jl>0 ? 1 : 0,ctrl.max_entries,used)
    used=_cdid_base_cells(4*Jl,d+4,ctrl.max_entries,used)
    out.peak_entries=used
    out.outcome_pooling=(ctrl.outcome_pooling=="" ? "separate" : ctrl.outcome_pooling)
    out.n_F=p.n_F;out.n_A=p.n_A;out.n_supports=S;out.n_blocks=B;out.n_consumers=Jl
    out.intercept=ctrl.intercept;out.p=d
    out.pair_mass=J(S,1,0);out.sumY=J(S,1,0);out.sumZ=J(S,d,0);out.rates=J(S,1,0)
    xmin=J(S,d,.);xmax=xmin;firstZ=xmin;summary=cdid_base_support_qr(S,1)
    groupmass=J(rows(p.active_cohorts),1,0)
    for(j=1;j<=rows(groupmass);j++) {
        ix=_cdid_base_index(p.a_cohort:==p.active_cohorts[j]);groupmass[j]=sum(w[ix])
    }
    has_csr=_cdid_base_build_sup_csr(P,p,sup_lo=.,sup_hi=.,sup_a=.,sup_rs=.,sup_rt=.)
    for(s=1;s<=S;s++) {
        pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
        if(has_csr) {
            if(sup_hi[s]>=sup_lo[s]) {
                ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                out.pair_mass[s]=sum(w[ai])
                yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                if(missing(yy))_cdid_base_fail("outcome","Pair outcome change overflowed.",504)
                ww=w[ai];rawZ=Z[ai,];D=L:*rawZ
                if(missing(D))_cdid_base_fail("outcome_design","Calendar-scaled outcome design overflowed.",504)
                if(d) {
                    xmin[s,]=colmin(rawZ);xmax[s,]=colmax(rawZ);firstZ[s,]=rawZ[1,]
                    out.sumZ[s,]=quadcross(ww,rawZ)
                }
                out.sumY[s]=sum(ww:*yy)
                F=_cdid_or_compress(((J(rows(ai),1,1),D,yy):*sqrt(ww)))
            }
            else F=J(0,d+2,0)
        }
        else {
            rows0=_cdid_base_pair_slice(P,p,pi,1,p.n_A)
            ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rows0[,1]]]:==p.support[s,2])
            ai=rows0[ix,1];out.pair_mass[s]=sum(w[ai])
            F=J(0,d+2,0);first=1
            for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                hi=min((p.n_A,lo+ctrl.chunk_rows-1))
                rows0=_cdid_base_read(P,p,pi,lo,hi,p.active_cohorts[p.support[s,2]])
                if(rows(rows0)==0)continue
                ai=rows0[,1];ww=w[ai];yy=rows0[,2];rawZ=Z[ai,];D=L:*rawZ
                if(missing(D))_cdid_base_fail("outcome_design","Calendar-scaled outcome design overflowed.",504)
                if(d) {
                    if(first){xmin[s,]=colmin(rawZ);xmax[s,]=colmax(rawZ);firstZ[s,]=rawZ[1,];}
                    else{xmin[s,]=colmin(xmin[s,]\rawZ);xmax[s,]=colmax(xmax[s,]\rawZ);}
                }
                first=0;out.sumY[s]=out.sumY[s]+sum(ww:*yy)
                if(d)out.sumZ[s,]=out.sumZ[s,]+quadcross(ww,rawZ)
                F=_cdid_or_compress(F\((J(rows(ai),1,1),D,yy):*sqrt(ww)))
            }
        }
        if(out.pair_mass[s]<=0 | out.pair_mass[s]>=.)
            _cdid_base_fail("support","Declared support has nonpositive/nonfinite observed weight mass; none is dropped.",459)
        out.rates[s]=out.pair_mass[s]/groupmass[p.support[s,2]];summary[s].F=F
    }
    values=cdid_ipw_rate_values(out.rates);out.rates_c=values[,1];out.rates_active=values[,2]
    nPseudo=2*(d+1)*B
    pseudoX=J(nPseudo,d,0);pseudoY=J(nPseudo,1,0);pseudoW=pseudoY;pseudoB=pseudoY;at=0
    for(b=1;b<=B;b++) {
        pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
        F=J(0,d+2,0);mass=0;constant=J(1,d,1);baseZ=J(1,d,0)
        for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
            s=p.block_support[m]
            F=_cdid_or_compress(F\(summary[s].F/sqrt(out.rates_c[s])))
            mass=mass+out.pair_mass[s]/out.rates_c[s]
            if(m==p.block_ptr[b])baseZ=firstZ[s,]
            if(d)constant=constant:*(xmin[s,]:==baseZ):*(xmax[s,]:==baseZ)
        }
        hqrd(F[,1],Q=.,tau=.,R=.)
        V=hqrdmultq(Q,tau,F[,(2..d+2)],1);coefrows=V[1,]/R[1,1]
        resfactor=V[2..rows(V),]
        if(ctrl.intercept & d)for(k=1;k<=d;k++)if(constant[k]) {
            coefrows[k]=L*baseZ[k];resfactor[,k]=J(rows(resfactor),1,0)
        }
        nn=rows(resfactor);scale=sqrt(nn/mass)
        A=coefrows:+scale:*resfactor;V=coefrows:-scale:*resfactor
        rr=(at+1..at+2*nn)'
        if(d)pseudoX[rr,]=(A[,1..d]\V[,1..d])
        pseudoY[rr]=(A[,d+1]\V[,d+1]);pseudoW[rr]=J(2*nn,1,mass/(2*nn));pseudoB[rr]=J(2*nn,1,b)
        at=at+2*nn
    }
    if(at!=nPseudo)_cdid_base_fail("summary","Unexpected compressed support factor dimension.",3200)
    if(out.outcome_pooling=="separate") {
        out.blocks=cdid_base_or_block(B,1)
        out.alpha=J(B,1,0);out.block_mass=J(B,1,0);out.beta=J(0,1,0)
        out.R=J(0,0,0);out.pivot=J(1,0,0);out.kept=J(1,0,0);out.dropped=J(1,0,0)
        out.alias_map=J(0,0,0);out.alias_constants=J(0,0,0);out.mean_D=J(B,0,0)
        out.beta_score_A=J(0,0,0);out.beta_rate_loading=J(0,S,0)
        at=0;max_k=0
        for(b=1;b<=B;b++) {
            nn=d+1;rr=(at+1..at+2*nn)';at=at+2*nn
            fit=cdid_pooled_wls_fit((d ? pseudoX[rr,] : J(2*nn,0,0)),pseudoY[rr],pseudoW[rr],
                J(2*nn,1,1),J(2*nn,1,1),1,ctrl.intercept,J(0,d,0),J(0,1,0),ctrl.rank_tol,1)
            out.blocks[b].beta=fit.beta;out.blocks[b].R=fit.R;out.blocks[b].pivot=fit.pivot
            out.blocks[b].kept=fit.kept;out.blocks[b].dropped=fit.dropped
            out.blocks[b].alias_map=fit.alias_map;out.blocks[b].alias_constants=fit.alias_constants
            out.blocks[b].mean_D=(ctrl.intercept ? fit.mean_X[,fit.kept] : J(1,rows(fit.beta),0))
            out.alpha[b]=fit.alpha[1];out.block_mass[b]=fit.mass[1]
            if(rows(fit.beta)>max_k)max_k=rows(fit.beta)
        }
        _cdid_base_or_alias(P,p,Z,out,ctrl.chunk_rows)
        out.alpha_rate=J(rows(p.block_support),1,0);out.n_fit_rows=0;out.objective=0
        for(b=1;b<=B;b++) {
            pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
            kb=rows(out.blocks[b].beta);n_ctrl_b=p.block_ptr[b+1]-p.block_ptr[b]
            force=J(kb,n_ctrl_b,0);sum_sa=0;sum_sb=J(1,kb,0)
            for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
                s=p.block_support[m]
                if(has_csr) {
                    if(sup_hi[s]<sup_lo[s])continue
                    ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                    yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                    fw=w[ai]/out.rates_c[s];D=L:*Z[ai,out.blocks[b].kept]
                    rawres=yy:-out.alpha[b]:-(kb ? D*out.blocks[b].beta : 0)
                    centered=(ctrl.intercept ? D:-out.blocks[b].mean_D : D)
                    coefrows=quadcross(fw:*rawres,centered)
                    if(ctrl.intercept)sum_sa=sum_sa+sum(fw:*rawres)
                    if(kb)sum_sb=sum_sb+coefrows
                    if(out.rates_active[s]) {
                        if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                        if(kb)force[,m-p.block_ptr[b]+1]=force[,m-p.block_ptr[b]+1]-coefrows'/out.rates_c[s]
                    }
                    out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                }
                else {
                    for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                        rows0=_cdid_base_read(P,p,pi,lo,min((p.n_A,lo+ctrl.chunk_rows-1)),p.active_cohorts[p.support[s,2]])
                        if(rows(rows0)==0)continue
                        ai=rows0[,1];fw=w[ai]/out.rates_c[s];D=L:*Z[ai,out.blocks[b].kept]
                        rawres=rows0[,2]:-out.alpha[b]:-(kb ? D*out.blocks[b].beta : 0)
                        centered=(ctrl.intercept ? D:-out.blocks[b].mean_D : D)
                        coefrows=quadcross(fw:*rawres,centered)
                        if(ctrl.intercept)sum_sa=sum_sa+sum(fw:*rawres)
                        if(kb)sum_sb=sum_sb+coefrows
                        if(out.rates_active[s]) {
                            if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                            if(kb)force[,m-p.block_ptr[b]+1]=force[,m-p.block_ptr[b]+1]-coefrows'/out.rates_c[s]
                        }
                        out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                    }
                }
            }
            out.blocks[b].score_mean_alpha=(ctrl.intercept ? sum_sa/p.n_F : 0)
            out.blocks[b].score_mean_beta=(kb ? sum_sb/p.n_F : J(1,0,0))
            out.blocks[b].beta_rate_loading=_cdid_nuisance_gram(out.blocks[b].R,out.blocks[b].pivot,force)
        }
        out.link_mu=J(Jl,1,0);out.link_c=J(Jl,max_k,0);out.link_d=J(Jl,1,0);out.link_rate_loading=J(Jl,1,0)
        for(j=1;j<=Jl;j++) {
            b=p.consumer[j,3];s=p.consumer[j,4];pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
            mass=out.pair_mass[s];kb=rows(out.blocks[b].beta)
            den=out.sumY[s]-out.alpha[b]*mass-(kb ? L*(out.sumZ[s,out.blocks[b].kept]*out.blocks[b].beta) : 0)
            out.link_mu[j]=den/mass
            if(kb)out.link_c[j,1..kb]=L:*out.sumZ[s,out.blocks[b].kept]/mass
            out.link_d[j]=mass/(p.n_F*out.rates_c[s])
            if(out.rates_active[s])out.link_rate_loading[j]=-(den-out.link_mu[j]*mass)/(p.n_F*out.rates_c[s]^2*out.link_d[j])
        }
        out.retained_entries=11+5*S+2*B+rows(p.block_support)+S*d+B*(max_k*max_k+5*max_k+1)+
            max_k*rows(p.block_support)+Jl*max_k+3*Jl
        return(out)
    }
    out.blocks=cdid_base_or_block(0,1)
    fit=cdid_pooled_wls_fit(pseudoX,pseudoY,pseudoW,pseudoB,J(nPseudo,1,1),1,
        ctrl.intercept,J(0,d,0),J(0,1,0),ctrl.rank_tol,1)
    out.beta=fit.beta;out.alpha=fit.alpha;out.R=fit.R;out.pivot=fit.pivot
    out.kept=fit.kept;out.dropped=fit.dropped;out.alias_map=fit.alias_map;out.alias_constants=fit.alias_constants
    out.mean_D=(ctrl.intercept ? fit.mean_X[,fit.kept] : J(B,rows(fit.beta),0));out.block_mass=fit.mass
    _cdid_base_or_alias(P,p,Z,out,ctrl.chunk_rows)
    k=rows(out.beta);out.beta_score_A=J(p.n_A,k,0);force=J(k,S,0)
    out.alpha_rate=J(rows(p.block_support),1,0);out.n_fit_rows=0;out.objective=0
    for(b=1;b<=B;b++) {
        pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
        for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
            s=p.block_support[m]
            if(has_csr) {
                if(sup_hi[s]<sup_lo[s])continue
                ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                fw=w[ai]/out.rates_c[s];D=L:*Z[ai,out.kept]
                rawres=yy:-out.alpha[b]:-D*out.beta
                centered=(ctrl.intercept ? D:-out.mean_D[b,] : D)
                coefrows=quadcross(fw:*rawres,centered)
                if(k)out.beta_score_A[ai,]=out.beta_score_A[ai,]+centered:*(fw:*rawres)
                if(out.rates_active[s]) {
                    if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                    if(k)force[,s]=force[,s]-coefrows'/out.rates_c[s]
                }
                out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
            }
            else {
                for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                    rows0=_cdid_base_read(P,p,pi,lo,min((p.n_A,lo+ctrl.chunk_rows-1)),p.active_cohorts[p.support[s,2]])
                    if(rows(rows0)==0)continue
                    ai=rows0[,1];fw=w[ai]/out.rates_c[s];D=L:*Z[ai,out.kept]
                    rawres=rows0[,2]:-out.alpha[b]:-D*out.beta
                    centered=(ctrl.intercept ? D:-out.mean_D[b,] : D)
                    coefrows=quadcross(fw:*rawres,centered)
                    if(k)out.beta_score_A[ai,]=out.beta_score_A[ai,]+centered:*(fw:*rawres)
                    if(out.rates_active[s]) {
                        if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                        if(k)force[,s]=force[,s]-coefrows'/out.rates_c[s]
                    }
                    out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                }
            }
        }
    }
    out.beta_rate_loading=_cdid_nuisance_gram(out.R,out.pivot,force)
    out.link_mu=J(Jl,1,0);out.link_c=J(Jl,k,0);out.link_d=J(Jl,1,0);out.link_rate_loading=J(Jl,1,0)
    for(j=1;j<=Jl;j++) {
        b=p.consumer[j,3];s=p.consumer[j,4];pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
        mass=out.pair_mass[s];den=out.sumY[s]-out.alpha[b]*mass-L*(out.sumZ[s,out.kept]*out.beta)
        out.link_mu[j]=den/mass;out.link_c[j,]=L:*out.sumZ[s,out.kept]/mass
        out.link_d[j]=mass/(p.n_F*out.rates_c[s])
        if(out.rates_active[s])out.link_rate_loading[j]=-(den-out.link_mu[j]*mass)/(p.n_F*out.rates_c[s]^2*out.link_d[j])
    }
    out.retained_entries=11+5*S+2*B+3*k+rows(p.block_support)+S*d+B*k+k*k+
        (k+B)*cols(out.dropped)+p.n_A*k+k*S+Jl*k+3*Jl+cols(out.dropped)
    return(out)
}

/* Membership and optional changes on arbitrary original-F rows, filtered to
   this actual consumer's cohort supports BEFORE reading outcome differences. */
real matrix _cdid_base_frame_pair(struct cdid_baseline_provider scalar f,
    real scalar pi,real colvector rowsF,real colvector groups,real scalar readY)
{
    real colvector aa,ix,keep,rs,rt,dy,cmap,ai
    real scalar s,t
    aa=f.plan.F_to_A[rowsF];ix=_cdid_base_index(aa:>0)
    if(rows(ix)==0)return(J(0,3,0))
    if(rows(groups)==1)keep=_cdid_base_index(f.plan.a_cohort[aa[ix]]:==groups[1])
    else {
        cmap=J(rows(f.plan.cohorts),1,0);cmap[groups]=J(rows(groups),1,1)
        keep=_cdid_base_index(cmap[f.plan.a_cohort[aa[ix]]])
    }
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=aa[ix]
    if(rows(f.plan.row_AT)==f.n_A) {
        s=f.plan.pair[pi,1];t=f.plan.pair[pi,2]
        rs=f.plan.row_AT[ai,s];rt=f.plan.row_AT[ai,t]
    }
    else {
        rs=cdid_population_rows(f.population,rowsF[ix],J(rows(ix),1,f.plan.calendar[f.plan.pair[pi,1]]))
        rt=cdid_population_rows(f.population,rowsF[ix],J(rows(ix),1,f.plan.calendar[f.plan.pair[pi,2]]))
    }
    keep=_cdid_base_index((rs:>0):&(rt:>0))
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=ai[keep];rs=rs[keep];rt=rt[keep]
    keep=_cdid_base_index((f.population.y[rs]:<.):&(f.population.y[rt]:<.))
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=ai[keep];rs=rs[keep];rt=rt[keep];dy=J(rows(ix),1,0)
    if(readY)dy=f.population.y[rt]-f.population.y[rs]
    if(missing(dy))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
    return((ix,ai,dy))
}
real colvector _cdid_base_rate_if_rr(struct cdid_baseline_provider scalar f,
    real scalar sid,real colvector rowsF,real matrix rr)
{
    real colvector out,aa,ix,keep
    real scalar cohort,mass
    out=J(rows(rowsF),1,0);aa=f.plan.F_to_A[rowsF]
    cohort=f.plan.active_cohorts[f.plan.support[sid,2]]
    mass=f.group_mass[f.plan.support[sid,2]]
    ix=_cdid_base_index(aa:>0)
    if(rows(ix)) {
        keep=_cdid_base_index(f.plan.a_cohort[aa[ix]]:==cohort)
        if(rows(keep)) {
            ix=ix[keep]
            out[ix]=(-f.n_F:*f.w_A[aa[ix]]:*f.outcome.rates[sid]):/mass
        }
    }
    if(rows(rr))out[rr[,1]]=out[rr[,1]]:+f.n_F:*f.w_A[rr[,2]]/mass
    return(out)
}
real colvector _cdid_base_rate_if(struct cdid_baseline_provider scalar f,
    real scalar sid,real colvector rowsF)
{
    real matrix rr
    real scalar cohort,pi
    cohort=f.plan.active_cohorts[f.plan.support[sid,2]];pi=f.plan.support[sid,1]
    rr=_cdid_base_frame_pair(f,pi,rowsF,cohort,0)
    return(_cdid_base_rate_if_rr(f,sid,rowsF,rr))
}
/* Every requested consumer shares this ONE beta IF chunk in pooled mode. */
real matrix _cdid_base_beta_if(struct cdid_baseline_provider scalar f,real colvector rowsF)
{
    real matrix score,out
    real colvector aa,ix,up
    real scalar k,s
    k=rows(f.outcome.beta);score=J(rows(rowsF),k,0)
    if(k==0)return(score)
    aa=f.plan.F_to_A[rowsF];ix=_cdid_base_index(aa:>0)
    if(rows(ix))score[ix,]=f.n_F:*f.outcome.beta_score_A[aa[ix],]
    out=_cdid_nuisance_gram(f.outcome.R,f.outcome.pivot,score')'
    for(s=1;s<=f.plan.n_supports;s++) {
        if(!any(f.outcome.beta_rate_loading[,s]))continue
        up=_cdid_base_rate_if(f,s,rowsF)
        out=out+up*f.outcome.beta_rate_loading[,s]'
    }
    return(out)
}
/* Per-block beta IF chunk on-the-fly for outcome_pooling == "separate". */
real matrix _cdid_base_blk_beta_if(struct cdid_baseline_provider scalar f,
    real scalar b,real colvector rowsF)
{
    real scalar n,kb,pi,L,m,s,c_s
    real matrix score,out,rr_all,rr,D_all,D,centered
    real colvector ai,fw,res,up,b_co,ix_s
    n=rows(rowsF);kb=rows(f.outcome.blocks[b].beta)
    if(kb==0)return(J(n,0,0))
    score=J(n,kb,0)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_base_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all)) D_all=L:*f.Z_A[rr_all[,2],f.outcome.blocks[b].kept]
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        if(rows(ix_s)) {
            rr=rr_all[ix_s,];ai=rr[,2];D=D_all[ix_s,]
            res=rr[,3]:-f.outcome.alpha[b]:-D*f.outcome.blocks[b].beta
            fw=f.n_F:*f.w_A[ai]/f.outcome.rates_c[s]
            centered=(f.outcome.intercept ? D:-f.outcome.blocks[b].mean_D : D)
            score[rr[,1],]=score[rr[,1],]+centered:*(fw:*res)
        }
    }
    score=score:-(f.n_F:*f.outcome.blocks[b].score_mean_beta)
    out=_cdid_nuisance_gram(f.outcome.blocks[b].R,f.outcome.blocks[b].pivot,score')'
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        if(!any(f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]))continue
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]'
    }
    return(out)
}
real colvector _cdid_base_alpha_if(struct cdid_baseline_provider scalar f,
    real scalar b,real colvector rowsF,real matrix betaIF)
{
    real colvector out,up,ai,res,fw,beta,b_co,ix_s
    real rowvector kept,mean_D
    real scalar m,s,pi,L,kb,sep,c_s
    real matrix rr_all,rr,D_all,D
    out=J(rows(rowsF),1,0)
    if(!f.outcome.intercept)return(out)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    mean_D=(sep ? f.outcome.blocks[b].mean_D : f.outcome.mean_D[b,])
    kb=rows(beta)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_base_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all) & kb) D_all=L:*f.Z_A[rr_all[,2],kept]
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        if(rows(rr)) {
            ai=rr[,2];D=(kb ? D_all[ix_s,] : J(rows(rr),0,0))
            res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0)
            fw=f.n_F:*f.w_A[ai]/f.outcome.rates_c[s]/f.outcome.block_mass[b]
            out[rr[,1]]=out[rr[,1]]+fw:*res
        }
        if(f.outcome.alpha_rate[m]!=0) {
            up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
            out=out+up*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(sep)out=out:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
    if(kb)out=out-betaIF*mean_D'
    return(out)
}
/* Combined separate-block betaIF and alphaIF sharing one pair query and D*beta. */
void _cdid_base_blk_beta_alpha_if(struct cdid_baseline_provider scalar f,
    real scalar b,real colvector rowsF,real matrix betaIF,real colvector alphaIF)
{
    real scalar n,kb,pi,L,m,s,c_s,do_alpha
    real matrix score,rr_all,rr,D_all,D,centered
    real colvector ai,fw,res_all,res,up,b_co,ix_s
    n=rows(rowsF);kb=rows(f.outcome.blocks[b].beta);do_alpha=(f.outcome.intercept!=0)
    betaIF=J(n,kb,0);alphaIF=J(n,1,0)
    if(kb==0 & !do_alpha) return
    score=J(n,kb,0)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_base_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all)) {
        if(kb) {
            D_all=L:*f.Z_A[rr_all[,2],f.outcome.blocks[b].kept]
            res_all=rr_all[,3]:-f.outcome.alpha[b]:-D_all*f.outcome.blocks[b].beta
        }
        else res_all=rr_all[,3]:-f.outcome.alpha[b]
    }
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        if(rows(rr)) {
            ai=rr[,2];res=res_all[ix_s]
            fw=f.n_F:*f.w_A[ai]/f.outcome.rates_c[s]
            if(kb) {
                D=D_all[ix_s,]
                centered=(f.outcome.intercept ? D:-f.outcome.blocks[b].mean_D : D)
                score[rr[,1],]=score[rr[,1],]+centered:*(fw:*res)
            }
            if(do_alpha) alphaIF[rr[,1]]=alphaIF[rr[,1]]+(fw/f.outcome.block_mass[b]):*res
        }
        if(do_alpha & f.outcome.alpha_rate[m]!=0) {
            up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
            alphaIF=alphaIF+up*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(kb) {
        score=score:-(f.n_F:*f.outcome.blocks[b].score_mean_beta)
        betaIF=_cdid_nuisance_gram(f.outcome.blocks[b].R,f.outcome.blocks[b].pivot,score')'
        for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
            if(!any(f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]))continue
            s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
            ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
            rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
            up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
            betaIF=betaIF+up*f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]'
        }
    }
    if(do_alpha) {
        alphaIF=alphaIF:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
        if(kb)alphaIF=alphaIF-betaIF*f.outcome.blocks[b].mean_D'
    }
}
real colvector _cdid_base_or_if(struct cdid_baseline_provider scalar f,
    real scalar ell,real colvector rowsF,real matrix betaIF,| real colvector alphaIF_in)
{
    real scalar b,s,pi,L,kb,sep
    real colvector out,ai,res,up,beta,alphaIF
    real rowvector kept
    real matrix rr,D
    b=f.plan.consumer[ell,3];s=f.plan.consumer[ell,4];pi=f.plan.support[s,1]
    L=f.plan.pair[pi,2]-f.plan.pair[pi,1];out=J(rows(rowsF),1,0)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    kb=rows(beta)
    rr=_cdid_base_frame_pair(f,pi,rowsF,f.plan.active_cohorts[f.plan.support[s,2]],1)
    if(rows(rr)) {
        ai=rr[,2];D=L:*f.Z_A[ai,kept]
        res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0):-f.outcome.link_mu[ell]
        out[rr[,1]]=f.n_F:*f.w_A[ai]/f.outcome.pair_mass[s]:*res
    }
    alphaIF=(args()>=5 & rows(alphaIF_in)==rows(rowsF) ? alphaIF_in : _cdid_base_alpha_if(f,b,rowsF,betaIF))
    out=out-alphaIF
    if(f.outcome.link_rate_loading[ell]!=0) {
        up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.link_rate_loading[ell]
    }
    if(kb)out=out-betaIF*f.outcome.link_c[ell,1..kb]'
    return(out)
}

real colvector _cdid_base_cohorts(struct cdid_baseline_plan_state scalar p,real scalar ell)
{
    real scalar b
    real colvector s
    b=p.consumer[ell,3];s=p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1]
    return(p.consumer[ell,1]\p.active_cohorts[p.support[s,2]])
}
real colvector _cdid_base_risk_obs(struct cdid_population scalar P,
    struct cdid_baseline_plan_state scalar p,real scalar ell)
{
    real matrix rr
    real colvector co,keep,cmap
    real scalar pi
    pi=p.consumer[ell,2];co=_cdid_base_cohorts(p,ell)
    rr=_cdid_base_pair_slice(P,p,pi,1,p.n_A)
    if(rows(rr)==0)return(J(0,1,0))
    cmap=J(rows(p.cohorts),1,0);cmap[co]=J(rows(co),1,1)
    keep=_cdid_base_index(cmap[p.a_cohort[rr[,1]]])
    return(rows(keep) ? rr[keep,1] : J(0,1,0))
}
/* One compact summary per exact outcome-free treatment recipe; no retained
   N-vector per model. A current full fitted nuisance is consumed only during
   this class's point fits, then released before the next class/replay. */
struct cdid_baseline_provider scalar cdid_baseline_prepare(struct cdid_baseline_plan_state scalar p,
    struct cdid_population scalar P,real matrix X_F,real matrix Z_F,real colvector w_F,
    struct cdid_baseline_controls scalar ctrl,
    | pointer(real colvector) colvector ext_co,struct cdid_ipw_nuisance colvector ext_nu)
{
    struct cdid_baseline_provider scalar out
    struct cdid_ipw_frame scalar frame
    struct cdid_ipw_nuisance scalar nu
    struct cdid_ipw_link scalar link
    transmorphic scalar co_buckets,frame_cache,map_cache,ext_buckets
    pointer(real colvector) colvector rep_co,rep_risk
    string scalar key_co
    real scalar N,A,nlink,px,pz,used,work,i,j,k,c,nclass,ell,rep,pi,b,L,hit,has_ext,ext_hit,j_cand,use_row_AT
    real colvector ids,classid,reps,risk,other,co,map,cls,keep,dy,ix,fw,beta,c_list,t_ix,c_ix,ext_list
    real rowvector kept
    real matrix rr,D,c_info
    N=p.n_F;A=p.n_A;nlink=p.n_consumers;px=(ctrl.method=="reg" ? 0 : cols(X_F));pz=(ctrl.method=="ipw" ? 0 : cols(Z_F))
    if(P.N!=N | P.n_A!=A | rows(w_F)!=N | cols(w_F)!=1 |
        rows(p.A_to_F)!=A | !all(P.A_to_F:==p.A_to_F))
        _cdid_base_fail("input","Prepared population/plan/weight dimensions disagree.",498)
    if(A)if(missing(w_F[P.A_to_F]) | any(w_F[P.A_to_F]:<=0))
        _cdid_base_fail("weight","Weights must be finite positive on A.",498)
    if(ctrl.method!="reg")if(rows(X_F)!=N | px<1)
        _cdid_base_fail("treatment_design","Supply explicit full-row treatment matrix with at least one column.",498)
    if(ctrl.method!="ipw")if(rows(Z_F)!=N | (!ctrl.intercept & pz==0))
        _cdid_base_fail("outcome_design","Supply explicit outcome columns/intercept with original-F row shape.",498)
    used=_cdid_base_cells(8,P.R+N+A,ctrl.max_entries,p.workspace_entries)
    used=_cdid_base_cells(8*nlink,px+2*rows(p.cohorts)+16,ctrl.max_entries,used)
    if(ctrl.method!="reg") {
        used=_cdid_base_cells(10*A,px+4,ctrl.max_entries,used)
        used=_cdid_base_cells(8*px,px,ctrl.max_entries,used)
        k=_cdid_base_cells(nlink,px+1,ctrl.max_entries,0)
        used=_cdid_base_cells(4*k,rows(p.cohorts)+4,ctrl.max_entries,used)
    }
    if(ctrl.method!="ipw")used=_cdid_base_cells(8*A,pz+4,ctrl.max_entries,used)
    out.method=ctrl.method;out.status=(nlink ? "prepared" : "no_supported_consumers")
    out.n_F=N;out.n_A=A;out.core_to_frame=(A ? N/A : .);out.normalization_count=A
    out.weight_row_domain=(ctrl.method=="ipw" ? "A" : "F")
    out.weight_rows=(ctrl.method=="ipw" ? A : N);out.sample_rank_upper=out.weight_rows-1
    out.n_or_fits=0;out.n_ps_fits=0;out.n_ps_reuses=0
    out.population=P;out.plan=p
    if(!(rows(out.plan.row_AT)==A & (out.plan.n_pairs<=0 | A<=0 | A*out.plan.n_pairs>4000000 |
        (rows(out.plan.pair_cum)==A+1 & cols(out.plan.pair_cum)==out.plan.n_pairs)))) {
        _cdid_base_build_pair_cache(P,out.plan)
    }
    out.w_A=(A ? w_F[P.A_to_F] : J(0,1,0))
    out.X_A=J(A,0,0);out.Z_A=J(A,0,0);out.points=J(nlink,1,.)
    out.group_mass=J(rows(out.plan.active_cohorts),1,0)
    for(j=1;j<=rows(out.plan.active_cohorts);j++) {
        ix=_cdid_base_index(out.plan.a_cohort:==out.plan.active_cohorts[j]);out.group_mass[j]=sum(out.w_A[ix])
    }
    out.point=cdid_base_point(nlink,1);out.propensity=cdid_base_ps(0,1);out.ps_cache=cdid_base_ps_cache(0,1)
    out.outcome.p=0;out.outcome.beta=J(0,1,0);out.outcome.alpha=J(0,1,0)
    out.workspace_entries=used;out.workspace_limit=ctrl.max_entries;out.chunk_rows=ctrl.chunk_rows
    if(nlink==0)return(out)
    if(ctrl.method!="reg") {
        out.X_A=X_F[P.A_to_F,]
        if(missing(out.X_A))_cdid_base_fail("treatment_design","Treatment design must be finite on A.",498)
        if(ctrl.method=="ipw" & missing(X_F))
            _cdid_base_fail("treatment_design","Pure public baseline IPW requires finite full-F treatment design.",498)
    }
    if(ctrl.method!="ipw") {
        out.Z_A=Z_F[P.A_to_F,]
        if(missing(out.Z_A))_cdid_base_fail("outcome_design","Outcome design must be finite on A.",498)
        out.outcome=_cdid_base_or_fit(P,out.plan,out.Z_A,out.w_A,ctrl,used)
        out.n_or_fits=(out.outcome.outcome_pooling=="separate" ? out.plan.n_blocks : 1)
        out.workspace_entries=_cdid_base_cells(1,out.outcome.retained_entries,ctrl.max_entries,used)
        if(ctrl.method=="reg") {
            out.points=out.outcome.link_mu
            for(ell=1;ell<=nlink;ell++){out.point[ell].estimate=out.points[ell];out.point[ell].or_id=ell;out.point[ell].ps_id=0;}
            return(out)
        }
    }
    /* Same model ingredients = ordered g/K, exact RISK observers, fixed X/w
       and the explicit GLM controls. No outcome/residual selects a class. */
    classid=J(nlink,1,0);reps=J(nlink,1,0);nclass=0;work=out.plan.work_ops
    co_buckets=asarray_create("string",1)
    rep_co=J(nlink,1,NULL);rep_risk=J(nlink,1,NULL)
    for(ell=1;ell<=nlink;ell++) {
        co=_cdid_base_cohorts(out.plan,ell);risk=_cdid_base_risk_obs(P,out.plan,ell);hit=0
        work=_cdid_base_work(A,ctrl.workcap,work)
        key_co=invtokens(strofreal(co)')
        if(asarray_contains(co_buckets,key_co)) {
            c_list=asarray(co_buckets,key_co)
            for(j=1;j<=rows(c_list);j++) {
                c=c_list[j];rep=reps[c]
                work=_cdid_base_work(A,ctrl.workcap,work)
                if(out.plan.consumer[ell,2]==out.plan.consumer[rep,2]) {hit=c;break;}
                other=*rep_risk[c]
                if(rows(risk)!=rows(other))continue
                if(all(risk:==other)){hit=c;break;}
            }
        }
        if(hit==0) {
            nclass++;reps[nclass]=ell;hit=nclass
            rep_co[nclass]=&(co:+0);rep_risk[nclass]=&(risk:+0)
            if(asarray_contains(co_buckets,key_co)) asarray(co_buckets,key_co,(asarray(co_buckets,key_co)\nclass))
            else asarray(co_buckets,key_co,nclass)
        }
        classid[ell]=hit
    }
    out.propensity=cdid_base_ps(nclass,1);out.ps_cache=cdid_base_ps_cache(nclass,1)
    ids=order((classid,(1::nlink)),(1,2))
    c_info=_cdid_ps_base(classid[ids],1)
    frame_cache=asarray_create("string",1)
    map_cache=asarray_create("string",1)
    has_ext=(args()>=8 & rows(ext_co)==rows(ext_nu) & rows(ext_co)>0)
    ext_buckets=asarray_create("string",1)
    if(has_ext) {
        for(j=1;j<=rows(ext_co);j++) {
            if(ext_co[j]==NULL) continue
            key_co=invtokens(strofreal(*ext_co[j])')
            if(asarray_contains(ext_buckets,key_co)) asarray(ext_buckets,key_co,(asarray(ext_buckets,key_co)\j))
            else asarray(ext_buckets,key_co,j)
        }
    }
    use_row_AT=(rows(out.plan.row_AT)==A)
    for(c=1;c<=nclass;c++) {
        rep=reps[c];co=*rep_co[c];risk=*rep_risk[c]
        key_co=invtokens(strofreal(co)')
        if(asarray_contains(map_cache,key_co)) map=asarray(map_cache,key_co)
        else {
            map=J(rows(out.plan.cohorts),1,0);map[co]=(1..rows(co))'
            asarray(map_cache,key_co,map)
        }
        ext_hit=0
        if(has_ext & asarray_contains(ext_buckets,key_co)) {
            ext_list=asarray(ext_buckets,key_co)
            for(j=1;j<=rows(ext_list);j++) {
                j_cand=ext_list[j]
                if(ext_nu[j_cand].n_A==A & ext_nu[j_cand].n_F==N & ext_nu[j_cand].p==px &
                    ext_nu[j_cand].epsilon==ctrl.glm_epsilon & ext_nu[j_cand].maxit==ctrl.glm_maxit &
                    rows(ext_nu[j_cand].risk_a)==rows(risk)) {
                    if(all(ext_nu[j_cand].risk_a:==risk)) {
                        nu=ext_nu[j_cand];ext_hit=1;break
                    }
                }
            }
        }
        if(!ext_hit) {
            if(asarray_contains(frame_cache,key_co)) frame=asarray(frame_cache,key_co)
            else {
                frame=cdid_ipw_frame_create(map[out.plan.a_cohort],J(A,1,1),out.w_A,out.plan.A_to_F,N,rows(co))
                asarray(frame_cache,key_co,frame)
            }
            nu=cdid_ipw_nuisance_prepare(frame,risk,out.X_A[risk,],ctrl.glm_epsilon,ctrl.glm_maxit)
        }
        out.n_ps_fits=out.n_ps_fits+1
        out.propensity[c].n_A=A;out.propensity[c].n_F=N;out.propensity[c].n_cohorts=rows(co)
        out.propensity[c].p=nu.p;out.propensity[c].cohorts=co;out.propensity[c].q=nu.q
        out.propensity[c].pair_count=nu.pair_count;out.propensity[c].pair_mass=nu.pair_mass
        out.propensity[c].rates=nu.rates;out.propensity[c].rates_c=nu.rates_c;out.propensity[c].rates_active=nu.rates_active
        out.propensity[c].beta=nu.beta;out.propensity[c].A_beta_r=nu.A_beta_r
        out.propensity[c].mean_a=nu.mean_a;out.propensity[c].mean_b=nu.mean_b
        out.propensity[c].iterations=nu.iterations;out.propensity[c].deviance=nu.deviance;out.propensity[c].score_max=nu.score_max
        out.ps_cache[c].rep_ell=rep;out.ps_cache[c].pi=out.plan.consumer[rep,2]
        out.ps_cache[c].risk_a=nu.risk_a;out.ps_cache[c].risk_pos_a=nu.risk_pos_a
        out.ps_cache[c].raw_probability=nu.raw_probability;out.ps_cache[c].odds=nu.odds
        out.ps_cache[c].cmap=map
        t_ix=_cdid_base_index(nu.risk_group:==1);c_ix=_cdid_base_index(nu.risk_group:>1)
        fw=nu.fit_weight[c_ix]:*nu.odds[c_ix]
        cls=ids[|c_info[c,1]\c_info[c,2]|]
        for(i=1;i<=rows(cls);i++) {
            ell=cls[i];pi=out.plan.consumer[ell,2];b=out.plan.consumer[ell,3];L=out.plan.pair[pi,2]-out.plan.pair[pi,1]
            if(use_row_AT) {
                if(rows(risk)) {
                    dy=P.y[out.plan.row_AT[risk,out.plan.pair[pi,2]]]-P.y[out.plan.row_AT[risk,out.plan.pair[pi,1]]]
                    if(missing(dy))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
                }
                else dy=J(0,1,0)
            }
            else {
                rr=_cdid_base_read(P,out.plan,pi,1,A,co)
                if(rows(rr)!=rows(risk))_cdid_base_fail("recipe","Observed risk rows changed during prepared replay.",498)
                if(!all(rr[,1]:==risk))_cdid_base_fail("recipe","Observed risk identity differs from shared propensity design block.",498)
                dy=rr[,2]
            }
            out.point[ell].ps_id=c;out.point[ell].or_id=(ctrl.method=="dr" ? ell : 0)
            out.point[ell].d_bar_diff=J(0,1,0)
            if(ctrl.method=="dr") {
                kept=(out.outcome.outcome_pooling=="separate" ? out.outcome.blocks[b].kept : out.outcome.kept)
                beta=(out.outcome.outcome_pooling=="separate" ? out.outcome.blocks[b].beta : out.outcome.beta)
                D=L:*out.Z_A[risk,kept]
                dy=dy:-out.outcome.alpha[b]:-(rows(beta) ? D*beta : 0)
                if(rows(beta)) {
                    out.point[ell].d_bar_diff=quadcross(nu.fit_weight[t_ix],D[t_ix,])'/nu.sum_a-quadcross(fw,D[c_ix,])'/nu.sum_b
                }
            }
            link=cdid_ipw_link_outcome(nu,dy);out.point[ell].residual_ipw=link
            out.points[ell]=link.estimate;out.point[ell].estimate=link.estimate
            if(i>1)out.n_ps_reuses=out.n_ps_reuses+1
        }
    }
    return(out)
}
/* Compact IPW component equations, with optional DR residualization
   and the global shared-beta channel. Complete original-F rows, zero off A. */
real colvector _cdid_base_ipw_if(struct cdid_baseline_provider scalar f,
    real scalar ell,real colvector rowsF,real matrix betaIF)
{
    struct cdid_base_ps scalar ps
    struct cdid_ipw_link scalar lk
    real scalar c,pi,b,L,factor
    real colvector out,aa,ix,cj,keep,cmap,zload,eta,dy,beta,fitw,t_ix,c_ix
    real rowvector kept
    real matrix rr,lv
    c=f.point[ell].ps_id;ps=f.propensity[c];lk=f.point[ell].residual_ipw
    pi=f.plan.consumer[ell,2];b=f.plan.consumer[ell,3];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    aa=f.plan.F_to_A[rowsF];out=J(rows(rowsF),1,0);factor=f.n_F/f.n_A
    cmap=J(rows(f.plan.cohorts),1,0);cmap[ps.cohorts]=(1::rows(ps.cohorts))
    ix=_cdid_base_index(aa:>0)
    if(rows(ix)) {
        cj=cmap[f.plan.a_cohort[aa[ix]]]
        keep=_cdid_base_index(cj:>0)
        if(rows(keep)) {
            ix=ix[keep];cj=cj[keep]
            out[ix]=(factor:*f.w_A[aa[ix]]):*lk.rate_coef_s0[cj]
        }
    }
    rr=_cdid_base_frame_pair(f,pi,rowsF,ps.cohorts,1)
    if(rows(rr)) {
        eta=f.X_A[rr[,2],]*ps.beta;lv=_cdid_ipw_logit_values(eta)
        zload=f.X_A[rr[,2],]*lk.loading_beta;dy=rr[,3]
        if(f.method=="dr") {
            kept=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].kept : f.outcome.kept)
            beta=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].beta : f.outcome.beta)
            dy=dy:-f.outcome.alpha[b]:-(rows(beta) ? L:*f.Z_A[rr[,2],kept]*beta : 0)
        }
        cj=cmap[f.plan.a_cohort[rr[,2]]]
        fitw=f.w_A[rr[,2]]:/ps.rates_c[cj]
        out[rr[,1]]=factor:*(f.w_A[rr[,2]]:*lk.rate_coef_s1[cj]+fitw:*((cj:==1):-lv[,1]):*zload)
        t_ix=_cdid_base_index(cj:==1)
        if(rows(t_ix))out[rr[t_ix,1]]=out[rr[t_ix,1]]+factor:*fitw[t_ix]:/ps.mean_a:*(dy[t_ix]:-lk.treated_mean)
        c_ix=_cdid_base_index(cj:>1)
        if(rows(c_ix))out[rr[c_ix,1]]=out[rr[c_ix,1]]-factor:*fitw[c_ix]:*lv[c_ix,2]:/ps.mean_b:*(dy[c_ix]:-lk.control_mean)
    }
    if(f.method=="dr" & rows(f.point[ell].d_bar_diff))out=out-betaIF*f.point[ell].d_bar_diff
    return(out)
}
/* Batched propensity-class IF evaluation across all requested consumers in class c
   (without the DR outcome-beta channel, which is subtracted per block / pooled). */
real matrix _cdid_base_ps_class_if(struct cdid_baseline_provider scalar f,
    real scalar c,real colvector rowsF,real colvector aa,real colvector ix_pos_a,
    real colvector ells)
{
    struct cdid_base_ps scalar ps
    real scalar m_c,B,factor,p_trt,n_coh,q,ell,pi,b,L,s_t,t_t,g_idx
    real colvector cmap,cj_a,keep,ix_a,ai_a,pos_v,ix_r,ai_r,pos_r,cj_r,w_r,fitw,raw_p,odds_r
    real colvector t_ix,c_ix,dy_one,beta,grp_key,grp_ord,sub_q
    real rowvector t_mean,c_mean,kept
    real matrix out_c,Rate0,Rate1,Lambda,Z_load,DY,OnPair,grp_info,rr_fb
    B=rows(rowsF);m_c=rows(ells);out_c=J(B,m_c,0)
    if(m_c==0 | B==0) return(out_c)
    if(rows(f.ps_cache)<c | rows(f.plan.row_AT)!=f.n_A) {
        for(q=1;q<=m_c;q++) out_c[,q]=_cdid_base_ipw_if(f,ells[q],rowsF,J(B,0,0))
        return(out_c)
    }
    ps=f.propensity[c];cmap=f.ps_cache[c].cmap
    n_coh=ps.n_cohorts;p_trt=ps.p;factor=f.n_F/f.n_A
    if(rows(ix_pos_a)==0) return(out_c)
    Rate0=J(n_coh,m_c,0);Rate1=J(n_coh,m_c,0);Lambda=J(p_trt,m_c,0)
    t_mean=J(1,m_c,0);c_mean=J(1,m_c,0)
    for(q=1;q<=m_c;q++) {
        ell=ells[q]
        Rate0[,q]=f.point[ell].residual_ipw.rate_coef_s0
        Rate1[,q]=f.point[ell].residual_ipw.rate_coef_s1
        Lambda[,q]=f.point[ell].residual_ipw.loading_beta
        t_mean[q]=f.point[ell].residual_ipw.treated_mean
        c_mean[q]=f.point[ell].residual_ipw.control_mean
    }
    ai_a=aa[ix_pos_a];cj_a=cmap[f.plan.a_cohort[ai_a]]
    keep=_cdid_base_index(cj_a:>0)
    if(rows(keep)) {
        ix_a=ix_pos_a[keep];cj_a=cj_a[keep]
        out_c[ix_a,]=(factor:*f.w_A[aa[ix_a]]):*Rate0[cj_a,]
    }
    if(rows(f.ps_cache[c].risk_pos_a)==f.n_A) {
        pos_v=f.ps_cache[c].risk_pos_a[aa[ix_pos_a]]
        keep=_cdid_base_index(pos_v:>0)
        if(rows(keep)==0) return(out_c)
        ix_r=ix_pos_a[keep];ai_r=aa[ix_r];pos_r=pos_v[keep]
    }
    else {
        rr_fb=_cdid_base_frame_pair(f,f.ps_cache[c].pi,rowsF,ps.cohorts,0)
        if(rows(rr_fb)==0) return(out_c)
        ix_r=rr_fb[,1];ai_r=rr_fb[,2]
        pos_r=J(rows(ai_r),1,0)
        for(q=1;q<=rows(ai_r);q++) pos_r[q]=_cdid_ipw_lower_bound(f.ps_cache[c].risk_a,ai_r[q])
    }
    raw_p=f.ps_cache[c].raw_probability[pos_r]
    odds_r=f.ps_cache[c].odds[pos_r]
    cj_r=cmap[f.plan.a_cohort[ai_r]]
    w_r=f.w_A[ai_r];fitw=w_r:/ps.rates_c[cj_r]
    Z_load=f.X_A[ai_r,]*Lambda
    OnPair=factor:*(w_r:*Rate1[cj_r,]+(fitw:*((cj_r:==1):-raw_p)):*Z_load)
    DY=J(rows(ai_r),m_c,0)
    grp_key=(f.method=="dr" ? f.plan.consumer[ells,3] : f.plan.consumer[ells,2])
    grp_ord=order((grp_key,(1::m_c)),(1,2))
    grp_info=_cdid_ps_base(grp_key[grp_ord],1)
    for(g_idx=1;g_idx<=rows(grp_info);g_idx++) {
        sub_q=grp_ord[|grp_info[g_idx,1]\grp_info[g_idx,2]|]
        if(f.method=="dr") {
            b=grp_key[sub_q[1]];pi=f.plan.block_pair[b]
        }
        else pi=grp_key[sub_q[1]]
        s_t=f.plan.pair[pi,1];t_t=f.plan.pair[pi,2]
        dy_one=f.population.y[f.plan.row_AT[ai_r,t_t]]-f.population.y[f.plan.row_AT[ai_r,s_t]]
        if(missing(dy_one))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
        if(f.method=="dr") {
            L=t_t-s_t
            kept=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].kept : f.outcome.kept)
            beta=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].beta : f.outcome.beta)
            dy_one=dy_one:-f.outcome.alpha[b]:-(rows(beta) ? L:*f.Z_A[ai_r,kept]*beta : 0)
        }
        if(rows(sub_q)==1) DY[,sub_q[1]]=dy_one
        else DY[,sub_q]=dy_one*J(1,rows(sub_q),1)
    }
    t_ix=_cdid_base_index(cj_r:==1)
    if(rows(t_ix)) OnPair[t_ix,]=OnPair[t_ix,]+(factor:*fitw[t_ix]:/ps.mean_a):*(DY[t_ix,]:-t_mean)
    c_ix=_cdid_base_index(cj_r:>1)
    if(rows(c_ix)) OnPair[c_ix,]=OnPair[c_ix,]-(factor:*fitw[c_ix]:*odds_r[c_ix]:/ps.mean_b):*(DY[c_ix,]:-c_mean)
    out_c[ix_r,]=OnPair
    return(out_c)
}
void _cdid_base_coh_split(struct cdid_baseline_provider scalar f,
    real colvector aa,real colvector ix_pos_a,
    pointer(real colvector) colvector coh_ix_f,pointer(real colvector) colvector coh_ai)
{
    real scalar n_ac,ac,g,n_pos
    real colvector empty_v,ai_all,ac_all,ord_ac,ac_sorted,sub
    real matrix ac_spans
    n_ac=rows(f.plan.active_cohorts)
    empty_v=J(0,1,0)
    coh_ix_f=J(n_ac,1,&empty_v);coh_ai=J(n_ac,1,&empty_v)
    n_pos=rows(ix_pos_a)
    if(n_pos==0 | n_ac==0) return
    ai_all=aa[ix_pos_a]
    ac_all=f.plan.cohort_to_active[f.plan.a_cohort[ai_all]]
    ord_ac=(n_pos==1 ? 1 : (all(ac_all[|2\n_pos|]:>=ac_all[|1\n_pos-1|]) ? (1::n_pos) : order((ac_all,(1::n_pos)),(1,2))))
    ac_sorted=ac_all[ord_ac]
    ac_spans=_cdid_ps_base(ac_sorted,1)
    for(g=1;g<=rows(ac_spans);g++) {
        ac=ac_sorted[ac_spans[g,1]]
        if(ac>0) {
            sub=ord_ac[|ac_spans[g,1]\ac_spans[g,2]|]
            coh_ix_f[ac]=&(ix_pos_a[sub]:+0)
            coh_ai[ac]=&(ai_all[sub]:+0)
        }
    }
}
void _cdid_base_ac_pair(struct cdid_baseline_provider scalar f,real scalar pi,
    real colvector ix_ac,real colvector ai_ac,real scalar readY,
    real colvector keep_obs,real colvector ix_obs,real colvector ai_obs,real colvector dy_obs)
{
    real scalar n_ac_r,s_t,t_t
    real colvector rs,rt
    n_ac_r=rows(ai_ac)
    if(n_ac_r==0) {
        keep_obs=ix_obs=ai_obs=J(0,1,0);dy_obs=J(0,1,0);return
    }
    s_t=f.plan.pair[pi,1];t_t=f.plan.pair[pi,2]
    rs=f.plan.row_AT[ai_ac,s_t];rt=f.plan.row_AT[ai_ac,t_t]
    if(all((rs:>0):&(rt:>0))) {
        if(readY) {
            dy_obs=f.population.y[rt]-f.population.y[rs]
            if(!missing(dy_obs)) {
                keep_obs=(1::n_ac_r);ix_obs=ix_ac;ai_obs=ai_ac;return
            }
        }
        else if(all((f.population.y[rs]:<.):&(f.population.y[rt]:<.))) {
            keep_obs=(1::n_ac_r);ix_obs=ix_ac;ai_obs=ai_ac;dy_obs=J(n_ac_r,1,0);return
        }
    }
    keep_obs=_cdid_base_index((rs:>0):&(rt:>0))
    if(rows(keep_obs)==0) {
        ix_obs=ai_obs=J(0,1,0);dy_obs=J(0,1,0);return
    }
    rs=rs[keep_obs];rt=rt[keep_obs]
    keep_obs=keep_obs[_cdid_base_index((f.population.y[rs]:<.):&(f.population.y[rt]:<.))]
    if(rows(keep_obs)==0) {
        ix_obs=ai_obs=J(0,1,0);dy_obs=J(0,1,0);return
    }
    ix_obs=ix_ac[keep_obs];ai_obs=ai_ac[keep_obs]
    if(readY) {
        dy_obs=f.population.y[f.plan.row_AT[ai_obs,t_t]]-f.population.y[f.plan.row_AT[ai_obs,s_t]]
        if(missing(dy_obs))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
    }
    else dy_obs=J(rows(keep_obs),1,0)
}
real colvector _cdid_base_rate_if_ac(struct cdid_baseline_provider scalar f,
    real scalar sid,real colvector ai_ac,real colvector keep_obs,real colvector ai_obs)
{
    real colvector up_ac
    real scalar mass
    if(rows(ai_ac)==0) return(J(0,1,0))
    mass=f.group_mass[f.plan.support[sid,2]]
    up_ac=(-f.n_F:*f.w_A[ai_ac]:*f.outcome.rates[sid]):/mass
    if(rows(keep_obs)) up_ac[keep_obs]=up_ac[keep_obs]:+f.n_F:*f.w_A[ai_obs]/mass
    return(up_ac)
}
real matrix _cdid_base_beta_if_fast(struct cdid_baseline_provider scalar f,
    real scalar B,real colvector aa,real colvector ix_pos_a,
    pointer(real colvector) colvector coh_ix_f,pointer(real colvector) colvector coh_ai)
{
    real matrix score,out
    real colvector ix_ac,ai_ac,keep_obs,ix_obs,ai_obs,dy_obs,up_ac
    real scalar k,s,ac
    k=rows(f.outcome.beta);score=J(B,k,0)
    if(k==0)return(score)
    if(rows(ix_pos_a))score[ix_pos_a,]=f.n_F:*f.outcome.beta_score_A[aa[ix_pos_a],]
    out=_cdid_nuisance_gram(f.outcome.R,f.outcome.pivot,score')'
    for(s=1;s<=f.plan.n_supports;s++) {
        if(!any(f.outcome.beta_rate_loading[,s]))continue
        ac=f.plan.support[s,2];ix_ac=*coh_ix_f[ac];ai_ac=*coh_ai[ac]
        if(rows(ix_ac)==0)continue
        _cdid_base_ac_pair(f,f.plan.support[s,1],ix_ac,ai_ac,0,keep_obs=.,ix_obs=.,ai_obs=.,dy_obs=.)
        up_ac=_cdid_base_rate_if_ac(f,s,ai_ac,keep_obs,ai_obs)
        out[ix_ac,]=out[ix_ac,]+up_ac*f.outcome.beta_rate_loading[,s]'
    }
    return(out)
}
real colvector _cdid_base_alpha_if_fast(struct cdid_baseline_provider scalar f,
    real scalar b,real scalar B,
    pointer(real colvector) colvector coh_ix_f,pointer(real colvector) colvector coh_ai,
    real matrix betaIF)
{
    real colvector out,up_ac,res,fw,beta,ix_ac,ai_ac,keep_obs,ix_obs,ai_obs,dy_obs
    real rowvector kept,mean_D
    real scalar m,s,ac,pi,L,kb,sep
    real matrix D
    out=J(B,1,0)
    if(!f.outcome.intercept)return(out)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    mean_D=(sep ? f.outcome.blocks[b].mean_D : f.outcome.mean_D[b,])
    kb=rows(beta)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];ac=f.plan.support[s,2]
        ix_ac=*coh_ix_f[ac];ai_ac=*coh_ai[ac]
        if(rows(ix_ac)==0)continue
        _cdid_base_ac_pair(f,pi,ix_ac,ai_ac,1,keep_obs=.,ix_obs=.,ai_obs=.,dy_obs=.)
        if(rows(ix_obs)) {
            D=(kb ? L:*f.Z_A[ai_obs,kept] : J(rows(ix_obs),0,0))
            res=dy_obs:-f.outcome.alpha[b]:-(kb ? D*beta : 0)
            fw=f.n_F:*f.w_A[ai_obs]/f.outcome.rates_c[s]/f.outcome.block_mass[b]
            out[ix_obs]=out[ix_obs]+fw:*res
        }
        if(f.outcome.alpha_rate[m]!=0) {
            up_ac=_cdid_base_rate_if_ac(f,s,ai_ac,keep_obs,ai_obs)
            out[ix_ac]=out[ix_ac]+up_ac*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(sep)out=out:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
    if(kb)out=out-betaIF*mean_D'
    return(out)
}
void _cdid_base_blk_ba_if_fast(struct cdid_baseline_provider scalar f,
    real scalar b,real scalar B,
    pointer(real colvector) colvector coh_ix_f,pointer(real colvector) colvector coh_ai,
    real scalar want_alpha,real matrix betaIF,real colvector alphaIF)
{
    real scalar kb,pi,L,m,s,ac,do_alpha
    real matrix score,D,centered
    real colvector fw,res,up_ac,ix_ac,ai_ac,keep_obs,ix_obs,ai_obs,dy_obs
    kb=rows(f.outcome.blocks[b].beta);do_alpha=(want_alpha!=0 & f.outcome.intercept!=0)
    betaIF=J(B,kb,0);alphaIF=J(B,1,0)
    if(kb==0 & !do_alpha) return
    score=J(B,kb,0)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];ac=f.plan.support[s,2]
        ix_ac=*coh_ix_f[ac];ai_ac=*coh_ai[ac]
        if(rows(ix_ac)==0)continue
        _cdid_base_ac_pair(f,pi,ix_ac,ai_ac,1,keep_obs=.,ix_obs=.,ai_obs=.,dy_obs=.)
        if(rows(ix_obs)) {
            if(kb) {
                D=L:*f.Z_A[ai_obs,f.outcome.blocks[b].kept]
                res=dy_obs:-f.outcome.alpha[b]:-D*f.outcome.blocks[b].beta
            }
            else res=dy_obs:-f.outcome.alpha[b]
            fw=f.n_F:*f.w_A[ai_obs]/f.outcome.rates_c[s]
            if(kb) {
                centered=(f.outcome.intercept ? D:-f.outcome.blocks[b].mean_D : D)
                score[ix_obs,]=score[ix_obs,]+centered:*(fw:*res)
            }
            if(do_alpha) alphaIF[ix_obs]=alphaIF[ix_obs]+(fw/f.outcome.block_mass[b]):*res
        }
        if(do_alpha & f.outcome.alpha_rate[m]!=0) {
            up_ac=_cdid_base_rate_if_ac(f,s,ai_ac,keep_obs,ai_obs)
            alphaIF[ix_ac]=alphaIF[ix_ac]+up_ac*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(kb) {
        score=score:-(f.n_F:*f.outcome.blocks[b].score_mean_beta)
        betaIF=_cdid_nuisance_gram(f.outcome.blocks[b].R,f.outcome.blocks[b].pivot,score')'
        for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
            if(!any(f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]))continue
            s=f.plan.block_support[m];ac=f.plan.support[s,2]
            ix_ac=*coh_ix_f[ac];ai_ac=*coh_ai[ac]
            if(rows(ix_ac)==0)continue
            _cdid_base_ac_pair(f,pi,ix_ac,ai_ac,0,keep_obs=.,ix_obs=.,ai_obs=.,dy_obs=.)
            up_ac=_cdid_base_rate_if_ac(f,s,ai_ac,keep_obs,ai_obs)
            betaIF[ix_ac,]=betaIF[ix_ac,]+up_ac*f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]'
        }
    }
    if(do_alpha) {
        alphaIF=alphaIF:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
        if(kb)alphaIF=alphaIF-betaIF*f.outcome.blocks[b].mean_D'
    }
}
real matrix _cdid_base_or_block_fast(struct cdid_baseline_provider scalar f,
    real scalar b,real colvector ells,real scalar B,
    pointer(real colvector) colvector coh_ix_f,pointer(real colvector) colvector coh_ai,
    real colvector alphaIF)
{
    real scalar m_b,i,ell,s,ac,pi,L,kb,sep
    real colvector beta,res,up_ac,ix_ac,ai_ac,keep_obs,ix_obs,ai_obs,dy_obs
    real rowvector kept
    real matrix out_b,D
    m_b=rows(ells)
    out_b=(0:-alphaIF)*J(1,m_b,1)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    kb=rows(beta)
    for(i=1;i<=m_b;i++) {
        ell=ells[i];s=f.plan.consumer[ell,4];pi=f.plan.support[s,1];ac=f.plan.support[s,2]
        ix_ac=*coh_ix_f[ac];ai_ac=*coh_ai[ac]
        if(rows(ix_ac)==0)continue
        L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
        _cdid_base_ac_pair(f,pi,ix_ac,ai_ac,1,keep_obs=.,ix_obs=.,ai_obs=.,dy_obs=.)
        if(rows(ix_obs)) {
            D=(kb ? L:*f.Z_A[ai_obs,kept] : J(rows(ix_obs),0,0))
            res=dy_obs:-f.outcome.alpha[b]:-(kb ? D*beta : 0):-f.outcome.link_mu[ell]
            out_b[ix_obs,i]=out_b[ix_obs,i]+f.n_F:*f.w_A[ai_obs]/f.outcome.pair_mass[s]:*res
        }
        if(f.outcome.link_rate_loading[ell]!=0) {
            up_ac=_cdid_base_rate_if_ac(f,s,ai_ac,keep_obs,ai_obs)
            out_b[ix_ac,i]=out_b[ix_ac,i]+up_ac*f.outcome.link_rate_loading[ell]
        }
    }
    return(out_b)
}
real colvector _cdid_base_or_if_nobeta(struct cdid_baseline_provider scalar f,
    real scalar ell,real colvector rowsF,real colvector alphaIF)
{
    real scalar b,s,pi,L,kb,sep
    real colvector out,ai,res,up,beta
    real rowvector kept
    real matrix rr,D
    b=f.plan.consumer[ell,3];s=f.plan.consumer[ell,4];pi=f.plan.support[s,1]
    L=f.plan.pair[pi,2]-f.plan.pair[pi,1];out=J(rows(rowsF),1,0)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    kb=rows(beta)
    rr=_cdid_base_frame_pair(f,pi,rowsF,f.plan.active_cohorts[f.plan.support[s,2]],1)
    if(rows(rr)) {
        ai=rr[,2];D=L:*f.Z_A[ai,kept]
        res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0):-f.outcome.link_mu[ell]
        out[rr[,1]]=f.n_F:*f.w_A[ai]/f.outcome.pair_mass[s]:*res
    }
    out=out-alphaIF
    if(f.outcome.link_rate_loading[ell]!=0) {
        up=_cdid_base_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.link_rate_loading[ell]
    }
    return(out)
}
real matrix cdid_baseline_if_rows(struct cdid_baseline_provider scalar f,
    real colvector rowsF,real colvector consumers,string scalar scale)
{
    real scalar B,nout,i,j,b,c,kb,k,used,use_ps_batch,use_fast_or
    real colvector ub,alphaIF,aa,ix_pos_a,ps_ids,ps_ord,idx_c,b_ids,b_ord,idx_b
    pointer(real colvector) colvector coh_ix_f,coh_ai
    real matrix betaIF,out,ps_info,b_info,Dbar,out_b
    B=rows(rowsF);nout=rows(consumers)
    if(cols(rowsF)!=1 | cols(consumers)!=1 | missing(rowsF) | missing(consumers) |
        any(rowsF:!=floor(rowsF)) | any(rowsF:<1) | any(rowsF:>f.n_F) |
        any(consumers:!=floor(consumers)) | any(consumers:<1) | any(consumers:>f.plan.n_consumers) |
        (B>1 ? (!all(rowsF[|2\B|]:>rowsF[|1\B-1|]) ? rows(uniqrows(sort(rowsF,1)))!=B : 0) : 0) |
        (nout>1 ? (!all(consumers[|2\nout|]:>consumers[|1\nout-1|]) ? rows(uniqrows(sort(consumers,1)))!=nout : 0) : 0) |
        B>f.chunk_rows | (scale!="frame" & scale!="core"))
        _cdid_base_fail("replay","Invalid original-F row/consumer batch or explicit IF scale.",498)
    used=_cdid_base_cells(8*B,nout+cols(f.X_A)+cols(f.Z_A)+4,f.workspace_limit,f.workspace_entries)
    out=J(B,nout,0);if(B==0 | nout==0)return(out)
    use_fast_or=(rows(f.plan.row_AT)==f.n_A)
    if(f.method=="reg") {
        b_ids=f.plan.consumer[consumers,3]
        b_ord=order((b_ids,(1::nout)),(1,2))
        b_info=_cdid_ps_base(b_ids[b_ord],1)
        if(use_fast_or) {
            aa=f.plan.F_to_A[rowsF];ix_pos_a=_cdid_base_index(aa:>0)
            _cdid_base_coh_split(f,aa,ix_pos_a,coh_ix_f,coh_ai)
            if(f.outcome.outcome_pooling=="separate") {
                for(j=1;j<=rows(b_info);j++) {
                    idx_b=b_ord[|b_info[j,1]\b_info[j,2]|];b=b_ids[idx_b[1]]
                    _cdid_base_blk_ba_if_fast(f,b,B,coh_ix_f,coh_ai,1,betaIF=.,alphaIF=.)
                    kb=rows(f.outcome.blocks[b].beta)
                    out_b=_cdid_base_or_block_fast(f,b,consumers[idx_b],B,coh_ix_f,coh_ai,alphaIF)
                    if(kb) out_b=out_b-betaIF*f.outcome.link_c[consumers[idx_b],1..kb]'
                    out[,idx_b]=out_b
                }
            }
            else {
                betaIF=_cdid_base_beta_if_fast(f,B,aa,ix_pos_a,coh_ix_f,coh_ai);k=rows(f.outcome.beta)
                for(j=1;j<=rows(b_info);j++) {
                    idx_b=b_ord[|b_info[j,1]\b_info[j,2]|];b=b_ids[idx_b[1]]
                    alphaIF=_cdid_base_alpha_if_fast(f,b,B,coh_ix_f,coh_ai,betaIF)
                    out[,idx_b]=_cdid_base_or_block_fast(f,b,consumers[idx_b],B,coh_ix_f,coh_ai,alphaIF)
                }
                if(k) out=out-betaIF*f.outcome.link_c[consumers,1..k]'
            }
        }
        else if(f.outcome.outcome_pooling=="separate") {
            for(j=1;j<=rows(b_info);j++) {
                idx_b=b_ord[|b_info[j,1]\b_info[j,2]|];b=b_ids[idx_b[1]]
                _cdid_base_blk_beta_alpha_if(f,b,rowsF,betaIF=.,alphaIF=.)
                kb=rows(f.outcome.blocks[b].beta)
                for(i=1;i<=rows(idx_b);i++) out[,idx_b[i]]=_cdid_base_or_if_nobeta(f,consumers[idx_b[i]],rowsF,alphaIF)
                if(kb) out[,idx_b]=out[,idx_b]-betaIF*f.outcome.link_c[consumers[idx_b],1..kb]'
            }
        }
        else {
            betaIF=_cdid_base_beta_if(f,rowsF);k=rows(f.outcome.beta)
            for(j=1;j<=rows(b_info);j++) {
                idx_b=b_ord[|b_info[j,1]\b_info[j,2]|];b=b_ids[idx_b[1]]
                alphaIF=_cdid_base_alpha_if(f,b,rowsF,betaIF)
                for(i=1;i<=rows(idx_b);i++) out[,idx_b[i]]=_cdid_base_or_if_nobeta(f,consumers[idx_b[i]],rowsF,alphaIF)
            }
            if(k) out=out-betaIF*f.outcome.link_c[consumers,1..k]'
        }
    }
    else {
        aa=f.plan.F_to_A[rowsF];ix_pos_a=_cdid_base_index(aa:>0)
        use_ps_batch=(rows(f.ps_cache)==rows(f.propensity) & rows(f.plan.row_AT)==f.n_A)
        if(use_ps_batch) {
            ps_ids=J(nout,1,0)
            for(i=1;i<=nout;i++) ps_ids[i]=f.point[consumers[i]].ps_id
            ps_ord=order((ps_ids,(1::nout)),(1,2))
            ps_info=_cdid_ps_base(ps_ids[ps_ord],1)
            for(j=1;j<=rows(ps_info);j++) {
                idx_c=ps_ord[|ps_info[j,1]\ps_info[j,2]|];c=ps_ids[idx_c[1]]
                out[,idx_c]=_cdid_base_ps_class_if(f,c,rowsF,aa,ix_pos_a,consumers[idx_c])
            }
            if(f.method=="dr") {
                _cdid_base_coh_split(f,aa,ix_pos_a,coh_ix_f,coh_ai)
                if(f.outcome.outcome_pooling=="separate") {
                    b_ids=f.plan.consumer[consumers,3]
                    b_ord=order((b_ids,(1::nout)),(1,2))
                    b_info=_cdid_ps_base(b_ids[b_ord],1)
                    for(j=1;j<=rows(b_info);j++) {
                        idx_b=b_ord[|b_info[j,1]\b_info[j,2]|];b=b_ids[idx_b[1]]
                        kb=rows(f.outcome.blocks[b].beta)
                        if(kb==0)continue
                        _cdid_base_blk_ba_if_fast(f,b,B,coh_ix_f,coh_ai,0,betaIF=.,alphaIF=.)
                        Dbar=J(kb,rows(idx_b),0)
                        for(i=1;i<=rows(idx_b);i++) {
                            if(rows(f.point[consumers[idx_b[i]]].d_bar_diff)==kb) Dbar[,i]=f.point[consumers[idx_b[i]]].d_bar_diff
                        }
                        out[,idx_b]=out[,idx_b]-betaIF*Dbar
                    }
                }
                else {
                    k=rows(f.outcome.beta)
                    if(k>0) {
                        betaIF=_cdid_base_beta_if_fast(f,B,aa,ix_pos_a,coh_ix_f,coh_ai)
                        Dbar=J(k,nout,0)
                        for(i=1;i<=nout;i++) {
                            if(rows(f.point[consumers[i]].d_bar_diff)==k) Dbar[,i]=f.point[consumers[i]].d_bar_diff
                        }
                        out=out-betaIF*Dbar
                    }
                }
            }
        }
        else {
            if(f.method=="dr" & f.outcome.outcome_pooling=="separate") {
                ub=uniqrows(sort(f.plan.consumer[consumers,3],1))
                for(j=1;j<=rows(ub);j++) {
                    b=ub[j];betaIF=_cdid_base_blk_beta_if(f,b,rowsF)
                    for(i=1;i<=nout;i++)if(f.plan.consumer[consumers[i],3]==b)
                        out[,i]=_cdid_base_ipw_if(f,consumers[i],rowsF,betaIF)
                }
            }
            else {
                betaIF=(f.method=="dr" ? _cdid_base_beta_if(f,rowsF) : J(B,0,0))
                for(i=1;i<=nout;i++) out[,i]=_cdid_base_ipw_if(f,consumers[i],rowsF,betaIF)
            }
        }
    }
    if(scale=="core")out=out/f.core_to_frame
    if(missing(out))_cdid_base_fail("numeric","Nonfinite complete original-unit IF.",504)
    return(out)
}
real matrix cdid_baseline_if_chunk(struct cdid_baseline_provider scalar f,
    real scalar firstF,real scalar lastF,real colvector consumers,string scalar scale)
{
    if(firstF>=. | lastF>=. | firstF<1 | lastF>f.n_F | firstF>lastF |
        firstF!=floor(firstF) | lastF!=floor(lastF))_cdid_base_fail("replay","Invalid original-F chunk endpoints.",498)
    return(cdid_baseline_if_rows(f,(firstF..lastF)',consumers,scale))
}
/* Source weighting row domain and normalizer are separate from IF scale. */
real matrix cdid_baseline_weight_chunk(struct cdid_baseline_provider scalar f,
    real scalar first,real scalar last,real colvector consumers)
{
    real colvector rowsF
    if(first>=. | last>=. | first<1 | last>f.weight_rows | first>last |
        first!=floor(first) | last!=floor(last))_cdid_base_fail("replay","Invalid weighting replay row endpoints.",498)
    rowsF=(first..last)'
    if(f.weight_row_domain=="A")rowsF=f.plan.A_to_F[rowsF]
    return(cdid_baseline_if_rows(f,rowsF,consumers,"core"))
}


/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
transmorphic cdid_mrbp_active()
{
    return(cdid_mref_active())
}

transmorphic cdid_mrbp_aux_missing(transmorphic a1)
{
    return(cdid_mref_aux_missing(a1))
}

transmorphic cdid_mrbp_aux_missing_reason(transmorphic a1)
{
    return(cdid_mref_aux_missing_reason(a1))
}

void cdid_mrbp_aux_note(transmorphic a1, transmorphic a2, transmorphic a3, transmorphic a4, transmorphic a5, transmorphic a6)
{
    cdid_mref_aux_note(a1, a2, a3, a4, a5, a6)
}

transmorphic cdid_mrbp_ctx_of(transmorphic a1)
{
    return(cdid_mref_ctx_of(a1))
}

transmorphic cdid_mrbp_excluded(transmorphic a1)
{
    return(cdid_mref_excluded(a1))
}

string scalar _cdid_base_num(real scalar x)
{
    return(x==. ? "Inf" : strtrim(strofreal(x,"%24.17g")))
}

end
