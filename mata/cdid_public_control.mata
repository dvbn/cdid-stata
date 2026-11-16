*! Private native translation of R chained_public_control.R31,57,210-286.
*! Numerical workspace1GiB; optional work cap; source-derived chunk and draw budgets.
version 19.5
mata:
struct cdid_public_binding_control {
    struct cdid_baseline_controls scalar baseline
    struct cdid_bridge_request scalar inference
    real scalar workspace_bytes,entries,chunk,workcap,links_bound,J_bound
    real scalar condition_trigger,ridge_relative
}
struct cdid_public_binding_control scalar cdid_public_binding_controls(
    struct cdid_public_bundle scalar B,real scalar intercept,real scalar workspace_bytes,
    real scalar chunk_override,real scalar workcap,real scalar condition_trigger,real scalar ridge_relative,
    real scalar compute_if,real scalar bstrap,real scalar cband,real scalar biters,real scalar alpha,
    string scalar posting_mode)
{
    struct cdid_public_binding_control scalar out
    struct cdid_bridge_request scalar Q
    struct cdid_infer_flags scalar flags
    real scalar N,C_eff,T,G,Jg,p,pairs,auxbound,parambound,br,ir,ifcells,signcells,perunit,exact
    exact=2^53-1;N=B.prepared.population.N;T=rows(B.request.calendar);G=rows(B.request.cohorts)
    Jg=rows(B.request.target_cohorts);p=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    C_eff=(B.prepared.clustered ? B.prepared.C : N)
    if(biters<1 | biters>=. | biters!=floor(biters) | biters>2^31-1)
        _cdid_pr_fail("cdid_chained_public_control_error","biters() must be a positive integer.")
    if(workspace_bytes>=. | workspace_bytes<8 | workspace_bytes!=floor(workspace_bytes) | workspace_bytes>exact)
        _cdid_pr_fail("cdid_chained_public_control_error","workspace() must be a positive integer number of bytes.")
    if(chunk_override>. | (chunk_override<. & (chunk_override<1 | chunk_override!=floor(chunk_override))))
        _cdid_pr_fail("cdid_chained_public_control_error","chunk() must be a positive integer.")
    if(workcap>. | (workcap<. & (workcap<1 | workcap!=floor(workcap) | workcap>exact)))
        _cdid_pr_fail("cdid_chained_public_control_error","workcap() must be a positive integer.")
    if(condition_trigger<=0 | condition_trigger>=1)
        _cdid_pr_fail("cdid_chained_public_control_error","condition() must be strictly between 0 and 1.")
    if(ridge_relative<0 | ridge_relative>=.)
        _cdid_pr_fail("cdid_chained_public_control_error","ridge() must be nonnegative.")
    if(workspace_bytes>=. | workspace_bytes<8 | workspace_bytes!=floor(workspace_bytes) | workspace_bytes>exact |
        chunk_override>. | (chunk_override<. & (chunk_override<1 | chunk_override!=floor(chunk_override))) |
        workcap>. | (workcap<. & (workcap<1 | workcap!=floor(workcap) | workcap>exact)) |
        condition_trigger<=0 | condition_trigger>=1 | ridge_relative<0 | ridge_relative>=. |
        biters<1 | biters>=. | biters!=floor(biters) | biters>2^31-1 |
        alpha<=0 | alpha>=1 | compute_if<0 | compute_if>1 | compute_if!=floor(compute_if) |
        (posting_mode!="matrix" & posting_mode!="state"))
        _cdid_pr_fail("cdid_chained_public_control_error","invalid computational or inference option")
    // Independent bootstrap/band laws are inherited unchanged from inference.
    flags=cdid_infer_resolve_flags(bstrap,cband,"cells")
    out.workspace_bytes=workspace_bytes;out.entries=floor(workspace_bytes/8)
    pairs=_cdid_gmm_cells(T,T-1,exact,0)/2
    out.links_bound=_cdid_gmm_cells(_cdid_gmm_cells(Jg,pairs,exact,0),G,exact,0)
    auxbound=_cdid_gmm_cells(2*Jg,T,exact,0)
    parambound=_cdid_gmm_cells(2,out.links_bound,exact,Jg)
    out.J_bound=parambound+auxbound+2*rows(B.request.post)+Jg
    if(out.J_bound>exact) _cdid_rs_fail("Source target-bound sum exceeds exact count arithmetic.",908)
    perunit=_cdid_gmm_cells(32,out.links_bound+p+1,exact,0)
    out.chunk=(chunk_override<. ? max((1,min((N,chunk_override)))) : max((1,min((N,16384,floor(out.entries/max((1,perunit))))))))
    out.workcap=workcap;out.condition_trigger=condition_trigger;out.ridge_relative=ridge_relative
    out.baseline=cdid_baseline_control_make(B.prepared.method,intercept,1e-12,50,1e-7,out.chunk,out.entries,workcap,16)
    out.baseline.outcome_pooling=(B.prepared.outcome_pooling=="" ? "separate" : B.prepared.outcome_pooling)
    br=(chunk_override<. ? min((biters,out.chunk)) : min((biters,out.chunk,max((1,floor(out.entries/max((1,16*C_eff))))))))
    ir=(chunk_override<. ? out.chunk : min((out.chunk,max((1,floor(out.entries/max((1,8*out.J_bound))))))))
    ifcells=_cdid_gmm_cells(ir,out.J_bound,exact,0);signcells=_cdid_gmm_cells(br,C_eff,exact,0)
    Q.budget=cdid_infer_budget_create(N,C_eff,out.J_bound,biters,ir,ifcells,br,signcells,out.entries,workspace_bytes)
    Q.limits.payload_slots=out.entries;Q.limits.payload_string_bytes=workspace_bytes
    // This is a representability ceiling, not a session RAM/work budget.
    Q.limits.registry_slots=Q.limits.registry_string_bytes=exact
    Q.alpha=alpha;Q.compute_inffunc=compute_if;Q.bstrap=bstrap;Q.cband=cband;Q.B=biters
    Q.C=C_eff;Q.unit_cluster=(B.prepared.clustered ? B.prepared.unit_cluster : (1::N));Q.clustered=B.prepared.clustered;Q.run_wald=1
    Q.report_cells=(rows(B.request.report) ? (1::rows(B.request.report)) : J(0,1,.))
    Q.pre_cells=_cdid_rb_index(B.request.report[,5]:==2);Q.posting_mode=posting_mode
    out.inference=Q
    return(out)
}
end
