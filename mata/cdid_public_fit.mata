*! Private actual long-data default fit assembly, StataNow 19.5.
*! Native reg/DR full-catalog route; provider/system/GMM components.
*! No statistical default fallback, user data access, or generated multiplier matrix.
version 19.5
mata:
struct cdid_public_fit_result scalar cdid_public_fit_bundle(struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_baseline_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_gmm_result scalar sol
    struct cdid_baseline_controls scalar base_ctrl
    real matrix Z,G,loadings,postmeta,Vw,Lx
    real colvector kept,rk,ap,gq
    struct _cdid_wsp_data scalar wsd
    real colvector point,diagnostic,reportpost
    string colvector diagnostic_names
    real scalar M,Jt,used,remaining,i,j,p,keepif,p_ctrl,def_pu,def_ch,eff_pu,eff_chunk,max_rem
    if(cdid_xp_requested())return(cdid_crosspair_public_fit(B,ctrl))
    if(cdid_etwfe_requested())return(cdid_etwfe_public_fit(B,ctrl))
    if(B.settings.gmm_weighting=="within") {
    }
    if(B.prepared.selection_est_method=="dr")return(cdid_augmented_public_fit(B,ctrl))
    if(B.prepared.selection_est_method=="imputation")return(cdid_imputation_public_fit(B,ctrl))
    if(B.prepared.selection_canonical=="conditional_mtar") {
        if(B.prepared.design=="reference_recipes") {
            if(cdid_nuis_requested())cdid_nuis_mark_selection("unavailable","reference_recipe_route",
                "With covariates or weights that vary over time, the observation logits are fitted separately for each reference period and are not reported.")
            return(cdid_selref_public_fit(B,ctrl))
        }
        return(cdid_selected_public_fit(B,ctrl))
    }
    if(B.prepared.design=="reference_recipes")return(cdid_baseline_reference_fit(B,ctrl))
    if(B.prepared.method=="ipw")return(cdid_frontend_ipw_fit(B,ctrl))
    if(B.prepared.design!="static_numeric" | (B.prepared.method!="dr" & B.prepared.method!="reg"))
        _cdid_pr_fail("cdid_native_private_binding_scope","This private full-catalog binding requires static numeric reg or dr; pure-IPW optimal must use the accepted structural basis, never a full-row substitute.")
    if(B.settings.gmm_weighting!="optimal" & B.settings.gmm_weighting!="two_step" & B.settings.gmm_weighting!="identity" & B.settings.gmm_weighting!="within")
        _cdid_pr_fail("cdid_chained_public_control_error","gmmweighting() must be optimal, two_step, identity, or within.")
    B.settings.alpha=ctrl.inference.alpha;B.settings.compute_if=ctrl.inference.compute_inffunc
    B.settings.bstrap=ctrl.inference.bstrap;B.settings.cband=ctrl.inference.cband;B.settings.biters=ctrl.inference.B
    co=cdid_catalog_default_opts();co.mem_limit_bytes=ctrl.workspace_bytes;co.work_cap=ctrl.workcap
    co.conventions="D10_D13_admissible_frame"
    pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
    C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
        B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
    D=cdid_chained_ipw_build_design(C,B.request.post,B.request.pt_assumption,B.request.pt_start)
    B.request=cdid_public_request_tags(B.request,C)
    // ALL original post rows and BOTH reporting families precede ONE shared OR.
    plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,ctrl.entries,ctrl.workcap)
    p_ctrl=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    def_pu=_cdid_gmm_cells(32,ctrl.links_bound+p_ctrl+1,2^53-1,0)
    def_ch=max((1,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,def_pu)))))))
    eff_chunk=ctrl.chunk
    if(ctrl.chunk==def_ch & ctrl.baseline.chunk_rows==ctrl.chunk) {
        eff_pu=_cdid_gmm_cells(32,rows(C.r_tp)+p_ctrl+1,2^53-1,0)
        eff_chunk=max((ctrl.chunk,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,eff_pu)))))))
    }
    base_ctrl=ctrl.baseline;base_ctrl.chunk_rows=eff_chunk
    p=cols(B.prepared.X)
    Z=(ctrl.baseline.intercept ? (p>1 ? B.prepared.X[,2..p] : J(B.prepared.population.N,0,.)) : B.prepared.X)
    fit=cdid_baseline_prepare(plan,B.prepared.population,B.prepared.X,Z,B.prepared.weight.owner_weight,base_ctrl)
    sys=cdid_baseline_make_system(C,D,B.request.post,ctrl.entries,fit.workspace_entries)
    support=cdid_public_post_ledger(C,D,B.request,sys)
    M=sys.M;Jt=rows(B.request.post);point=J(Jt,1,.);loadings=J(0,M,.)
    diagnostic=J(0,1,.);diagnostic_names=J(0,1,"")
    if(M>0 & sys.p_retained>0) {
        if(eff_chunk>ctrl.chunk) {
            max_rem=floor((ctrl.entries-fit.workspace_entries-sys.workspace_entries-3*M*M)/max((1,10*(M+cols(fit.X_A)+cols(fit.Z_A)+4))))
            eff_chunk=(max_rem>=ctrl.chunk ? min((eff_chunk,max_rem)) : ctrl.chunk)
        }
        used=fit.workspace_entries+sys.workspace_entries
        used=_cdid_base_cells(8*eff_chunk,M+cols(fit.X_A)+cols(fit.Z_A)+4,ctrl.entries,used)
        if(B.settings.gmm_weighting=="identity") {
            used=_cdid_gmm_identity_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries
            sol=cdid_gmm_identity_map(fit.points[plan.post_consumer],sys.full_A,sys.H_est,remaining)
        }
        else {
            used=_cdid_base_cells(2*eff_chunk,M,ctrl.entries,used)
            used=_cdid_base_cells(3*M,M,ctrl.entries,used)
            used=_cdid_gmm_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries-M*M
            if(B.settings.gmm_weighting=="within") {
                rk=fit.points[plan.post_consumer]
                if(B.prepared.outcome_pooling=="pooled")
                    cdid_wsg_within(B,plan,C,rk,0,ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,eff_chunk,kept,G,Vw)
                else {
                    wsd=_cdid_wsp_prepare(C,B.prepared.population,B.prepared.X,B.prepared.weight.owner_weight)
                    cdid_wsp_forest(wsd,C,kept,G)
                    Vw=cdid_wsp_cov(wsd,C,kept,B.prepared.method,ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit)
                }
                sol=cdid_gmm_fit_within_norm(rk[kept],sys.full_A[kept,],Vw,G,M,sys.H_est,
                    ctrl.condition_trigger,ctrl.ridge_relative,fit.weight_rows,fit.normalization_count,
                    eff_chunk,&cdid_baseline_post_core(),&fit,remaining,M,kept')
                Lx=J(rows(sol.target_loadings),M,0);Lx[.,kept']=sol.target_loadings;sol.target_loadings=Lx
                Lx=J(0,0,.);Vw=J(0,0,.)
            }
            else {
                G=I(M)
                sol=cdid_gmm_fit_unit_norm(fit.points[plan.post_consumer],sys.full_A,G,M,sys.H_est,
                    B.settings.gmm_weighting,ctrl.condition_trigger,ctrl.ridge_relative,
                    fit.weight_rows,fit.normalization_count,eff_chunk,&cdid_baseline_post_core(),&fit,remaining)
            }
        }
        if(sys.n_estimable) point[sys.estimable_ids]=sol.targets
        loadings=sol.target_loadings
        ap=J(rows(B.request.aux),1,.);gq=selectindex(fit.plan.aux_supported)
        if(rows(gq)) ap[gq]=fit.points[fit.plan.aux_consumer[gq]]
        rk=fit.points[plan.post_consumer]
        if(B.settings.gmm_weighting=="within") {
            Lx=loadings[.,kept']
            cdid_chain_pending(cdid_chain_sp(C,D,kept,cdid_chained_ipw_build_design(_cdid_chain_subcat(C,kept),
                B.request.post,B.request.pt_assumption,B.request.pt_start),rk[kept],Lx,sys.estimable_ids,B.request.report,
                B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,ap,fit.plan.aux_supported,
                B.request.pt_assumption,B.request.pt_start,B.settings.gmm_weighting))
            Lx=J(0,0,.)
        }
        else cdid_chain_pending(cdid_chain_sp(C,D,(1::M),D,rk,loadings,sys.estimable_ids,B.request.report,
            B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,ap,fit.plan.aux_supported,
            B.request.pt_assumption,B.request.pt_start,B.settings.gmm_weighting))
        diagnostic=(sol.criterion\sol.lambda\sol.rho\sol.trace_B\sol.covariance_rank\sol.source_n_units\sol.normalization_factor\sol.sample_rank_upper)
        diagnostic_names=("GMM_criterion"\"GMM_lambda"\"GMM_ridge_rho"\"GMM_trace_B"\"GMM_covariance_rank"\"GMM_source_replay_rows"\"GMM_gram_normalization_factor"\"GMM_sample_rank_upper")
        // Only final reported IF is retained for durable no-refit results.
        G=J(0,0,.);sol.unit_gram_covariance=sol.coefficient_loadings=sol.target_loadings=J(0,0,.)
        sol.coefficient_covariance=sol.covariance=J(0,0,.)
    }
    else if(sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input","An estimable post target has no GMM parameter system.")
    sys.full_A=sys.H_est=J(0,0,.)
    out.payload=cdid_public_provider_report(B,fit,point,sys.post_available,support.status,support.reason,
        sys.estimable_ids,loadings,ctrl.inference.compute_inffunc,ctrl.inference.limits,eff_chunk,ctrl.entries,sol.retained_if)
    sol.retained_if=J(0,0,.)
    postmeta=J(rows(B.request.report),cols(support.ledger),.)
    reportpost=_cdid_rb_index(B.request.report[,5]:==3)
    for(i=1;i<=rows(reportpost);i++){j=reportpost[i];postmeta[j,]=support.ledger[B.request.report[j,6],];}
    out.payload.ledger=out.payload.ledger,postmeta
    out.payload.ledger_columns=out.payload.ledger_columns,"post_support.":+support.columns
    out.payload.numeric_settings=out.payload.numeric_settings\diagnostic\ctrl.workspace_bytes\ctrl.entries\ctrl.chunk\ctrl.workcap\ctrl.condition_trigger\ctrl.ridge_relative\ctrl.links_bound\ctrl.J_bound
    out.payload.numeric_setting_names=out.payload.numeric_setting_names\diagnostic_names\"numeric_workspace_bytes"\"numeric_workspace_entries"\"chunk_rows"\"counted_work_limit"\"condition_trigger"\"ridge_relative"\"catalogue_links_bound"\"inference_targets_bound"
    out.payload.text_settings=out.payload.text_settings\("GMM_metric","full original catalogue identity metric before covariance weighting")\("registry_limit_meaning","2^53-1 exact-count ceiling; not a session RAM/work budget")
    out.payload.source_versions=out.payload.source_versions\("gmm","GMM combination, normalized by the units with two observed outcomes; outcome-regression and doubly robust comparisons")
    out.inference=ctrl.inference
    cdid_state_validate(out.payload,ctrl.inference.limits)
    return(out)
}
void cdid_public_default_store(string scalar bundle_name,string scalar request_name,
    string scalar weighting,real scalar intercept,real scalar workspace_bytes,
    real scalar chunk,real scalar workcap,real scalar trigger,real scalar ridge,
    real scalar compute_if,real scalar bstrap,real scalar cband,real scalar biters,
    real scalar alpha,string scalar posting_mode)
{
    struct cdid_public_bundle scalar B
    struct cdid_public_binding_control scalar ctrl
    struct cdid_public_fit_result scalar fit
    pointer() scalar p
    string scalar token
    p=findexternal(bundle_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Prepared long-data input is absent.")
    if(findexternal(request_name)!=NULL) _cdid_pr_fail("cdid_chained_request_input_error","Output request name already exists; no user object overwritten.")
    B=*p;B.settings.gmm_weighting=weighting
    ctrl=cdid_public_binding_controls(B,intercept,workspace_bytes,chunk,workcap,trigger,ridge,compute_if,bstrap,cband,biters,alpha,posting_mode)
    fit=cdid_public_fit_bundle(B,ctrl)
    token=cdid_state_put(fit.payload,ctrl.inference.limits);st_local("fit_token",token)
    p=crexternal(request_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Cannot allocate the explicit inference request.")
    *p=fit.inference
}
end
