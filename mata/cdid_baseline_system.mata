*! Common baseline post-system handoff; no nuisance fits, no IF rebase.
*! Reuses catalogue gauge/target predicates and complete-unit GMM.
version 19.5
mata:
struct cdid_baseline_system {
    real scalar M,p_retained,n_estimable,workspace_entries
    real matrix full_A,H_est
    real colvector kept_params,dropped_gauge_params,dropped_gauge_vertices
    real colvector estimable_ids,post_available,post_status,post_cat_supported,post_cat_reason
    real colvector post_param,post_vertex
    string colvector post_status_name,post_reason
}
struct cdid_baseline_system scalar cdid_baseline_make_system(struct cdid_cat scalar C,
    struct cdid_design scalar D,real matrix req,real scalar cap,real scalar retained_entries)
{
    struct cdid_baseline_system scalar o
    real scalar M,P,K,i,g,t,v,j,ps,pe,rs,re,used
    real colvector gq,tq,hit,keep,full_to_ret
    real matrix cl,sq
    M=rows(C.r_g);P=D.n_params;K=rows(req)
    if(cap<1 | cap>=. | cap!=floor(cap) | cap>2^53-1 | retained_entries<0 | retained_entries>=. |
        retained_entries!=floor(retained_entries) | cols(req)!=2 | missing(req) |
        (M>0 & D.n_links!=M))_cdid_base_fail("system","Invalid explicit post-system dimensions/budget.",498)
    /* Count before dense matrices or target/gauge arrays are allocated. */
    used=_cdid_base_cells(24,K,cap,retained_entries)
    used=_cdid_base_cells(8,P+D.n_vertices,cap,used)
    used=_cdid_base_cells(M,D.n_rank,cap,used)
    used=_cdid_base_cells(K,D.n_rank,cap,used)
    o.workspace_entries=used-retained_entries;o.M=M;o.p_retained=D.n_rank
    o.post_available=J(K,1,0);o.post_status=J(K,1,0);o.post_param=J(K,1,0);o.post_vertex=J(K,1,0)
    o.post_status_name=J(K,1,"");o.post_reason=J(K,1,"");gq=J(K,1,0);tq=J(K,1,0)
    for(i=1;i<=K;i++) {
        g=_cdid_cohort_index(C.cohorts,C.nt_index,req[i,1]);t=_cdid_bsearch(C.calendar,req[i,2])
        if(g==0 | t==0 | req[i,2]<req[i,1])_cdid_base_fail("system","Post target outside catalogue/calendar.",459)
        if(!any(C.tg:==g))_cdid_base_fail("system","Post cohort outside catalogue targets.",459)
        gq[i]=g;tq[i]=t
    }
    if(K>1) {
        sq=sort((gq,tq),(1,2))
        if(any((sq[2..K,1]:==sq[1..K-1,1]):&(sq[2..K,2]:==sq[1..K-1,2])))
            _cdid_base_fail("system","Duplicate requested post target.",459)
    }
    cl=_cdid_cat_classify(C,gq,tq);o.post_cat_supported=cl[,1];o.post_cat_reason=cl[,3]
    for(i=1;i<=K;i++) {
        if(D.built) {
            hit=_cdid_base_index((D.v_g:==gq[i]):&(D.v_date:==tq[i]):&(D.v_kind:==6))
            if(rows(hit)==1) {
                v=hit[1];o.post_vertex[i]=v;o.post_param[i]=D.v_param[v]
                o.post_status[i]=(D.v_conn[v] ? 1 : (D.v_degree[v]==0 ? 2 : 3))
                o.post_available[i]=(o.post_status[i]==1)
            }
        }
        o.post_status_name[i]=cdid_post_status_name(o.post_status[i])
        if(!o.post_available[i]) {
            if(M==0)o.post_reason[i]="empty_catalog:"+cdid_request_reason_name(o.post_cat_reason[i])
            else if(!o.post_cat_supported[i])o.post_reason[i]=cdid_request_reason_name(o.post_cat_reason[i])
            else o.post_reason[i]=o.post_status_name[i]
        }
    }
    o.estimable_ids=_cdid_base_index(o.post_available);o.n_estimable=rows(o.estimable_ids)
    keep=J(P,1,1);o.dropped_gauge_vertices=J(0,1,0);o.dropped_gauge_params=J(0,1,0)
    if(D.nullity>0) {
        o.dropped_gauge_vertices=D.comp_roots[D.unanchored]
        o.dropped_gauge_params=D.v_param[o.dropped_gauge_vertices]
        if(any(o.dropped_gauge_params:<1) | any(o.dropped_gauge_params:>P))
            _cdid_base_fail("system","Invalid accepted forest gauge root.",504)
        keep[o.dropped_gauge_params]=J(D.nullity,1,0)
    }
    o.kept_params=_cdid_base_index(keep)
    if(rows(o.kept_params)!=o.p_retained)_cdid_base_fail("system","Gauge count differs from accepted design rank.",504)
    full_to_ret=J(P,1,0)
    if(o.p_retained)full_to_ret[o.kept_params]=(1..o.p_retained)'
    o.full_A=J(M,o.p_retained,0);o.H_est=J(o.n_estimable,o.p_retained,0)
    for(i=1;i<=M;i++) {
        ps=D.l_start[i];pe=D.l_end[i]
        if(ps>0){rs=full_to_ret[ps];if(rs>0)o.full_A[i,rs]=-1;}
        if(pe>0){re=full_to_ret[pe];if(re>0)o.full_A[i,re]=1;}
    }
    for(j=1;j<=o.n_estimable;j++) {
        ps=o.post_param[o.estimable_ids[j]]
        if(ps<1 | ps>P)_cdid_base_fail("system","Invalid estimable parameter.",504)
        rs=full_to_ret[ps]
        if(rs<1)_cdid_base_fail("system","Estimable level was dropped as a gauge.",504)
        o.H_est[j,rs]=1
    }
    return(o)
}
/* Pointer lifetime belongs to the calling fit. Complete original C row order,
   even when identical raw consumers were deduplicated by the provider. */
real matrix cdid_baseline_post_core(pointer() scalar context,real scalar lo,real scalar hi)
{
    struct cdid_baseline_provider scalar f
    real colvector ids,map
    real matrix z
    real scalar j,k
    if(context==NULL)_cdid_base_fail("system","Missing prepared provider context.",498)
    f=*context;ids=uniqrows(sort(f.plan.post_consumer,1));map=J(f.plan.M,1,0)
    for(j=1;j<=f.plan.M;j++)for(k=1;k<=rows(ids);k++)if(f.plan.post_consumer[j]==ids[k])map[j]=k
    z=cdid_baseline_weight_chunk(f,lo,hi,ids)
    return(z[,map])
}
end
