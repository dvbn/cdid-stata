*! Public baseline IPW composition: exact post-treatment estimator plus explicit aux-only IPW.
*! Full original catalogue is preserved; source fullM/G/ridge branch is not replaced.
version 19.5
mata:
struct cdid_public_fit_result scalar cdid_frontend_ipw_fit(
    struct cdid_public_bundle scalar B,struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_catalog_opts scalar co
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D, Du
    struct cdid_chained_ipw_opts scalar options
    struct cdid_chained_ipw_result scalar post
    struct cdid_baseline_provider scalar aux
    struct cdid_baseline_controls scalar base_ctrl
    struct cdid_public_post_support scalar support
    struct cdid_public_report_buffer scalar buffer
    pointer(real colvector) colvector ps_co
    struct cdid_ipw_nuisance colvector ps_nu
    real colvector aux_point,good,qids,consumers,map,hit,reportpost,cu
    real matrix auxIF,values,postmeta,postIF
    real scalar i,j,Q,lo,hi,used,perrow,chunk,p_ctrl,def_pu,def_ch,eff_pu,eff_chunk
    if(B.prepared.method!="ipw" | B.prepared.design!="static_numeric")
        _cdid_pr_fail("cdid_chained_request_input_error","Baseline IPW requires the explicit static numeric treatment design.")
    B.settings.alpha=ctrl.inference.alpha;B.settings.compute_if=ctrl.inference.compute_inffunc
    B.settings.bstrap=ctrl.inference.bstrap;B.settings.cband=ctrl.inference.cband;B.settings.biters=ctrl.inference.B
    co=cdid_catalog_default_opts();co.mem_limit_bytes=ctrl.workspace_bytes;co.work_cap=ctrl.workcap
    co.conventions="D10_D13_admissible_frame"
    pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
    C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
        B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
    D=cdid_chained_ipw_build_design(C,B.request.post,B.request.pt_assumption,B.request.pt_start)
    B.request=cdid_public_request_tags(B.request,C)
    p_ctrl=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    def_pu=_cdid_gmm_cells(32,ctrl.links_bound+p_ctrl+1,2^53-1,0)
    def_ch=max((1,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,def_pu)))))))
    eff_chunk=ctrl.chunk
    if(ctrl.chunk==def_ch & ctrl.baseline.chunk_rows==ctrl.chunk) {
        eff_pu=_cdid_gmm_cells(32,rows(C.r_tp)+p_ctrl+1,2^53-1,0)
        eff_chunk=max((ctrl.chunk,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,eff_pu)))))))
    }
    base_ctrl=ctrl.baseline;base_ctrl.chunk_rows=eff_chunk
    options=cdid_chained_ipw_default_opts()
    options.est_method="ipw";options.gmm_weighting=B.settings.gmm_weighting
    options.pt_assumption=B.request.pt_assumption;options.pt_start=B.request.pt_start
    options.selection="none";options.fix_weights=B.prepared.weight.fix_weights
    options.dated_design=0;options.glm_epsilon=ctrl.baseline.glm_epsilon;options.glm_maxit=ctrl.baseline.glm_maxit
    options.condition_trigger=ctrl.condition_trigger;options.ridge_relative=ctrl.ridge_relative
    options.chunk_size=eff_chunk;options.max_matrix_entries=ctrl.entries;options.work_cap=ctrl.workcap
    options.max_cached_nuisances=ctrl.baseline.max_cached
    options.materialize_moment_covariance=0;options.reduce_structural_redundancy=1
    if(B.settings.gmm_weighting=="within")options.within_cov=&cdid_wsp_ipw_cov()
    // This is the exact full-catalog estimator. Identity keeps all
    // rows; optimal follows retained basis/G=E'E/fullM and the source handler.
    post=cdid_chained_ipw_fit_opts(B.prepared.population,C,D,B.prepared.X,
        B.prepared.weight.owner_weight,B.request.post,options,ctrl.inference.compute_inffunc,1,&ps_co,&ps_nu)
    support=cdid_frontend_post_support(C,D,B.request,post.target_available,post.target_design_status,
        post.target_cat_supported,post.target_cat_reason,post.target_vertex,post.target_param)
    // No post catalogue is passed or cropped here. Only BOTH auxiliary families
    // are planned/fitted; matching post propensity nuisances are reused directly.
    aux=cdid_ipw_auxiliary_prepare(B.prepared.population,B.request.cohorts,B.request.aux,
        B.request.aux_K_ptr,B.request.aux_K_index,B.prepared.X,B.prepared.weight.owner_weight,base_ctrl,ps_co,ps_nu)
    Q=rows(B.request.aux);aux_point=J(Q,1,.);good=_cdid_rb_index(aux.plan.aux_supported)
    if(rows(good))aux_point[good]=aux.points[aux.plan.aux_consumer[good]]
    buffer=cdid_public_report_start(B,post.att,post.target_available,support.status,support.reason,
        aux_point,aux.plan.aux_supported,aux.plan.aux_reason,ctrl.inference.compute_inffunc,ctrl.inference.limits)
    if(ctrl.inference.compute_inffunc) {
        qids=J(0,1,.)
        for(i=1;i<=rows(B.request.report);i++)if(B.request.report[i,5]==2 & buffer.payload.available[i])qids=qids\B.request.report[i,7]
        if(rows(qids))qids=uniqrows(sort(qids,1))
        consumers=(rows(qids) ? uniqrows(sort(aux.plan.aux_consumer[qids],1)) : J(0,1,.))
        map=J(rows(qids),1,.)
        for(i=1;i<=rows(qids);i++){hit=_cdid_rb_index(consumers:==aux.plan.aux_consumer[qids[i]]);map[i]=hit[1];}
        used=cdid_state_cost(buffer.payload)[1]+aux.workspace_entries+post.workspace_entries_bound
        perrow=10*(rows(consumers)+cols(aux.X_A)+4)+3*(rows(qids)+rows(B.request.post)+rows(B.request.report))
        if(used>=ctrl.entries | perrow>ctrl.entries-used)_cdid_rs_fail("Combined post/aux/report workspace exceeds explicit numeric limits.",908)
        chunk=min((eff_chunk,aux.chunk_rows,floor((ctrl.entries-used)/perrow)))
        for(lo=1;lo<=post.N_F;lo=lo+chunk) {
            hi=min((post.N_F,lo+chunk-1))
            values=(rows(consumers) ? cdid_baseline_if_chunk(aux,lo,hi,consumers,"frame") : J(hi-lo+1,0,.))
            auxIF=(rows(qids) ? values[,map'] : J(hi-lo+1,0,.))
            postIF=(rows(post.estimable_target_indices) ? post.target_if[lo..hi,post.estimable_target_indices'] : J(hi-lo+1,0,.))
            cdid_public_report_add(buffer,(lo::hi),post.estimable_target_indices,postIF,qids,auxIF)
        }
    }
    if(post.gmm_status=="solved" & rows(post.estimable_target_indices)) {
        cu=(cols(post.gmm.target_loadings)==post.M_links ? (1::post.M_links) : post.row_basis.retained)
        if(rows(cu)==post.M_links) Du=D
        else Du=cdid_chained_ipw_build_design(_cdid_chain_subcat(C,cu),B.request.post,B.request.pt_assumption,B.request.pt_start)
        cdid_chain_pending(cdid_chain_sp(C,D,cu,Du,post.moments[cu],post.gmm.target_loadings,post.estimable_target_indices,
            B.request.report,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,aux_point,aux.plan.aux_supported,
            B.request.pt_assumption,B.request.pt_start,B.settings.gmm_weighting))
    }
    out.payload=cdid_public_report_finish(buffer,ctrl.inference.limits)
    reportpost=_cdid_rb_index(B.request.report[,5]:==3);postmeta=J(rows(B.request.report),cols(support.ledger),.)
    for(i=1;i<=rows(reportpost);i++){j=reportpost[i];postmeta[j,]=support.ledger[B.request.report[j,6],];}
    out.payload.ledger=out.payload.ledger,postmeta;out.payload.ledger_columns=out.payload.ledger_columns,"post_support.":+support.columns
    out.payload.numeric_settings=out.payload.numeric_settings\post.M_links\post.row_basis.L\post.n_nuisance_fits\aux.n_ps_fits\aux.plan.M\aux.n_or_fits\post.gmm.source_n_units\
        post.n_A\post.gmm.lambda\post.gmm.rho\post.gmm.trace_B\post.criterion\ctrl.workspace_bytes\ctrl.entries\ctrl.chunk\ctrl.workcap\
        ctrl.condition_trigger\ctrl.ridge_relative\ctrl.links_bound\ctrl.J_bound
    out.payload.numeric_setting_names=out.payload.numeric_setting_names\"M_full"\"L_retained"\"n_post_nuisance_fits"\"n_aux_nuisance_fits"\"aux_post_rows_planned"\"n_pooled_OR_fits"\"GMM_replay_rows"\
        "GMM_normalization_count"\"GMM_lambda"\"GMM_ridge_rho"\"GMM_trace_B"\"GMM_criterion"\"numeric_workspace_bytes"\"numeric_workspace_entries"\"chunk_rows"\"counted_work_limit"\
        "condition_trigger"\"ridge_relative"\"catalogue_links_bound"\"inference_targets_bound"
    out.payload.text_settings=out.payload.text_settings\("GMM_metric","exact S3 retained basis, full-catalogue E'E and fullM")\("GMM_status",post.gmm_status)\("GMM_weight_row_domain","A for source raw weight IF; identity map has no weight Gram")\
        ("GMM_covariance_handler",post.gmm.covariance_handler)\("GMM_trigger_reason",post.gmm.trigger_reason)\
        ("auxiliary_fit_scope","explicit aux-only IPW domain; zero original post rows planned")\("registry_limit_meaning","exact-count ceiling; not session RAM/work limit")
    out.payload.source_versions=out.payload.source_versions\("estimator","inverse probability weighting with optional influence functions")\("auxiliary","inverse probability weighting for auxiliary comparisons")
    out.inference=ctrl.inference;cdid_state_validate(out.payload,ctrl.inference.limits)
    return(out)
}
end
