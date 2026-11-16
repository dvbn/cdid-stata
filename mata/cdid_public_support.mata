*! Private reporting-only source status/reason translation.
*! R cdid1039120 chained_ipw_pipeline.R3320-3399; no design/fit changes.
version 19.5
mata:
struct cdid_public_post_support {
    string colvector status,reason
    real matrix ledger
    string rowvector columns
}
struct cdid_public_post_support scalar cdid_frontend_post_support(
    struct cdid_cat scalar C,struct cdid_design scalar D,struct cdid_public_request scalar R,
    real colvector available,real colvector status,real colvector cat_supported,
    real colvector cat_reason,real colvector vertex,real colvector param)
{
    struct cdid_public_post_support scalar out
    real scalar Jt,i,j,g,t,hasrows,touch,qualified,component,v
    real colvector hit,K
    string scalar pool,pooltxt,targetstatus,catreason
    Jt=rows(R.post)
    if(rows(available)!=Jt | rows(status)!=Jt | rows(cat_supported)!=Jt | rows(cat_reason)!=Jt | rows(vertex)!=Jt | rows(param)!=Jt)
        _cdid_pr_fail("cdid_native_reporting_input","Passed post support vectors do not align with the declared table of cells.")
    out.status=out.reason=J(Jt,1,"")
    out.ledger=J(Jt,8,.)
    out.columns=("catalogue_rows_touching","qualified_post_rows_touching","catalogue_supported","catalogue_reason","design_status","design_vertex","design_parameter","component")
    for(i=1;i<=Jt;i++) {
        g=_cdid_cohort_index(C.cohorts,C.nt_index,R.post[i,1]);t=_cdid_bsearch(C.calendar,R.post[i,2])
        hasrows=touch=qualified=0;component=.;targetstatus=catreason=""
        for(j=1;j<=rows(C.r_g);j++) if(C.r_g[j]==g) {
            hasrows=1;touch=touch+(C.r_s[j]==t | C.r_t[j]==t)
            qualified=qualified+((C.r_s[j]==t & C.r_srole[j]==5) | (C.r_t[j]==t & C.r_trole[j]==5))
        }
        K=(R.post_K_ptr[i]<R.post_K_ptr[i+1] ? R.post_K_index[|R.post_K_ptr[i]\R.post_K_ptr[i+1]-1|] : J(0,1,.));pool=""
        for(j=1;j<=rows(K);j++)pool=pool+(j>1 ? "," : "")+_cdid_pr_num(R.cohorts[K[j]])
        pooltxt="clean control pool at t="+_cdid_pr_num(R.post[i,2])+" under the declared control group: {"+pool+"}"
        hit=_cdid_rb_index(C.tg:==g)
        if(rows(hit))targetstatus=cdid_target_status_name(C.tg_status[hit[1]])
        if(cat_reason[i]>0)catreason=cdid_request_reason_name(cat_reason[i])
        v=vertex[i];if(v>0)component=D.v_comp[v]
        if(!hasrows) {
            out.status[i]="rowless_cohort"
            out.reason[i]="no retained comparison catalog row exists for cohort g="+_cdid_pr_num(R.post[i,1])+
                (targetstatus!="" ? " (catalog target status: "+targetstatus+")" : "")+"; "+pooltxt
        }
        else if(status[i]==0) {
            out.status[i]="not_represented";out.reason[i]="cell is not represented in the structural design; "+pooltxt
        }
        else if(available[i]) {
            out.status[i]="estimable"
            out.reason[i]=(D.pt_all ? "level anchored to the reference or a represented PT-all pre-date through retained comparison rows (prop:graph)" : "level anchored to the cohort reference through retained comparison rows (prop:graph)")
        }
        else if(status[i]==2) {
            out.status[i]="isolated_unsupported"
            out.reason[i]="no retained comparison row has (g="+_cdid_pr_num(R.post[i,1])+",t="+_cdid_pr_num(R.post[i,2])+
                ") as a qualified post endpoint; catalog rows touching this date: "+_cdid_rb_key(touch)+
                ", qualified post rows touching it: "+_cdid_rb_key(qualified)+"; "+pooltxt+
                (catreason!="" ? "; catalog request reason: "+catreason : "")
        }
        else {
            out.status[i]="disconnected_component"
            out.reason[i]="cell lies in design component "+_cdid_rb_key(component)+
                ", which no retained comparison row connects to the cohort reference; its level is a gauge coordinate, not an ATT estimate (prop:def); estimable contrasts inside the component are retained in the post-contrast basis; "+pooltxt
        }
        out.ledger[i,]=(touch,qualified,cat_supported[i],cat_reason[i],status[i],v,param[i],component)
    }
    return(out)
}
struct cdid_public_post_support scalar cdid_public_post_ledger(
    struct cdid_cat scalar C,struct cdid_design scalar D,struct cdid_public_request scalar R,
    struct cdid_baseline_system scalar sys)
{
    return(cdid_frontend_post_support(C,D,R,sys.post_available,sys.post_status,
        sys.post_cat_supported,sys.post_cat_reason,sys.post_vertex,sys.post_param))
}
end
