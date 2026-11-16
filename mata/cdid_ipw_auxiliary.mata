*! Explicit auxiliary-only IPW plan; currentStataNow19.5.
*! Native adaptation of the cdid_baseline_plan implementation.
*! R cdid (GPL-3).
*! No catalogue cropping: the interface declares only auxiliary descriptors.
*! Full population/shared dictionaries and a COMPLETE provider plan are retained.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_ipwaux(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_baseline_plan_state scalar cdid_ipw_auxiliary_plan(struct cdid_population scalar P,
    real colvector cohorts,real matrix aux,real colvector aux_K_ptr,real colvector aux_K_index,
    real scalar cap,real scalar workcap)
{
    struct cdid_baseline_plan_state scalar p
    real scalar R,M,Q,Ktotal,used,work,i,j,k,u,a,g,s,t,op,pool,nK,at,pi,nc,nb,hit,b,old,ordj,S
    real scalar rr,sz,total_support,ci,ai,nt_index
    real colvector canonical
    real colvector K,rawpair,counts,members,crep,brep,cb,bo,bmap,ix,active,raworder,c_ord
    real matrix pairrows,needed,consumer,c_info
    if(cap>=. | cap<1 | cap!=floor(cap) | cap>2^53-1 | workcap>. |
        (workcap<. & (workcap<1 | workcap!=floor(workcap))))
        _cdid_base_fail("bound","Supply explicit positive numeric-entry cap and missing/positive work cap.",498)
    M=0;Q=rows(aux);R=Q
    canonical=cdid_cipw_cohorts(P)
    if(rows(cohorts)!=rows(canonical))_cdid_base_fail("dictionary","Auxiliary cohort dictionary differs from the full population.",498)
    if(rows(cohorts))if(!all(cohorts:==canonical))_cdid_base_fail("dictionary","Auxiliary cohort identities differ from the full population.",498)
    nt_index=_cdid_check_cohorts(cohorts,0)
    if(cols(aux)!=3 | missing(aux) | rows(aux_K_ptr)!=Q+1 | cols(aux_K_ptr)!=1 |
        cols(aux_K_index)!=1 | missing(aux_K_ptr) | missing(aux_K_index) |
        any(aux_K_ptr:!=floor(aux_K_ptr)) | any(aux_K_index:!=floor(aux_K_index)) |
        aux_K_ptr[1]!=1 | aux_K_ptr[Q+1]!=rows(aux_K_index)+1 |
        any(aux_K_index:<1) | any(aux_K_index:>rows(cohorts)))
        _cdid_base_fail("request","Invalid explicit auxiliary descriptor/control CSR.",498)
    Ktotal=rows(aux_K_index)
    used=_cdid_base_cells(96,R,cap,0)
    used=_cdid_base_cells(20,Ktotal,cap,used)
    used=_cdid_base_cells(16,P.N+P.n_A+rows(cohorts),cap,used)
    p.n_F=P.N;p.n_A=P.n_A;p.M=M;p.Q=Q;p.calendar=P.calendar;p.cohorts=cohorts
    p.workspace_limit=cap;p.A_to_F=P.A_to_F;p.F_to_A=P.F_to_A
    p.raw=J(R,3,0);p.raw_K_ptr=J(R+1,1,1);p.raw_K_index=J(Ktotal,1,0)
    p.raw_consumer=J(R,1,0);p.raw_supported=J(R,1,0);p.raw_pair_counts=J(R,1,0)
    p.raw_count_ptr=J(R+1,1,1);p.raw_group_count=J(R+Ktotal,1,0)
    p.post_consumer=J(M,1,0);p.aux_consumer=J(Q,1,0);p.aux_supported=J(Q,1,0)
    p.aux_reason=J(Q,1,"");p.a_cohort=J(P.n_A,1,0);p.cohort_A_count=J(rows(cohorts),1,0)
    p.cohort_to_active=J(rows(cohorts),1,0)
    if(P.n_A>0) {
        if(any(P.cohort_known[P.A_to_F]:!=1))_cdid_base_fail("cohort","A-unit cohort must be known.",459)
        p.a_cohort=_cdid_cohort_index_vec(cohorts,nt_index,P.cohort[P.A_to_F],J(P.n_A,1,1))
        if(any(p.a_cohort:==0))_cdid_base_fail("dictionary","A-unit cohort absent from catalogue dictionary.",498)
        c_ord=sort(p.a_cohort,1)
        c_info=_cdid_ps_ipwaux(c_ord,1)
        p.cohort_A_count[c_ord[c_info[,1]]]=c_info[,2]:-c_info[,1]:+1
    }
    p.active_cohorts=_cdid_base_index(p.cohort_A_count:>0)
    if(rows(p.active_cohorts))p.cohort_to_active[p.active_cohorts]=(1..rows(p.active_cohorts))'
    at=0
    for(i=1;i<=R;i++) {
        j=i;nK=aux_K_ptr[j+1]-aux_K_ptr[j]
        if(nK<1)_cdid_base_fail("request","An auxiliary descriptor requires a nonempty declared control set.",498)
        K=aux_K_index[aux_K_ptr[j]..aux_K_ptr[j+1]-1]
        g=_cdid_cohort_index(cohorts,nt_index,aux[j,1])
        s=_cdid_bsearch(P.calendar,aux[j,2]);t=_cdid_bsearch(P.calendar,aux[j,3])
        if(g==0 | aux[j,1]>=. | s==0 | t==0 | s>=t)
            _cdid_base_fail("request","Auxiliary g/s/t absent or unordered in shared dictionaries.",498)
        p.raw[i,]=(g,s,t)
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
    for(pi=1;pi<=p.n_pairs;pi++) {
        pairrows=_cdid_base_pair_slice(P,p,pi,1,P.n_A)
        work=_cdid_base_work(2*P.n_A*(1+ceil(ln(rows(P.calendar)+1)/ln(2))),workcap,work)
        counts=J(rows(cohorts),1,0)
        if(rows(pairrows)>0) {
            c_ord=sort(p.a_cohort[pairrows[,1]],1)
            c_info=_cdid_ps_ipwaux(c_ord,1)
            counts[c_ord[c_info[,1]]]=c_info[,2]-c_info[,1]:+1
        }
        members=_cdid_base_index(rawpair:==pi)
        for(j=1;j<=rows(members);j++) {
            i=members[j];g=p.raw[i,1];K=p.raw_K_index[p.raw_K_ptr[i]..p.raw_K_ptr[i+1]-1]
            p.raw_pair_counts[i]=rows(pairrows)
            p.raw_group_count[p.raw_count_ptr[i]..p.raw_count_ptr[i+1]-1]=counts[g\K]
            p.raw_supported[i]=all(p.cohort_A_count[g\K]:>0) & all(counts[g\K]:>0)
            if(!p.raw_supported[i]) {
                p.aux_reason[i]="unsupported_pair"
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
struct cdid_baseline_provider scalar cdid_ipw_auxiliary_prepare(
    struct cdid_population scalar P,real colvector cohorts,real matrix aux,
    real colvector Kptr,real colvector Kindex,real matrix X_F,real colvector w_F,
    struct cdid_baseline_controls scalar controls,
    | pointer(real colvector) colvector ext_co,struct cdid_ipw_nuisance colvector ext_nu)
{
    struct cdid_baseline_plan_state scalar plan
    if(controls.method!="ipw")_cdid_base_fail("method","Auxiliary-only preparation requires explicit IPW controls; no OR model is substituted.",498)
    plan=cdid_ipw_auxiliary_plan(P,cohorts,aux,Kptr,Kindex,controls.max_entries,controls.workcap)
    if(args()>=10) return(cdid_baseline_prepare(plan,P,X_F,J(P.N,0,.),w_F,controls,ext_co,ext_nu))
    return(cdid_baseline_prepare(plan,P,X_F,J(P.N,0,.),w_F,controls))
}
end
