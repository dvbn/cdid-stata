*! Private complete reporting scatter for native long-data binding, StataNow 19.5.
*! GPL-3 adaptation of R cdid chained_ipw_pipeline.R3785-3890.
*! Consumes fitted points/IF chunks; no nuisance/GMM fitting or statistical RNG.
version 19.5
mata:
struct cdid_public_report_buffer {
    struct cdid_result_state scalar payload
    real colvector seen
    real scalar compute_if
}
struct cdid_public_report_buffer scalar cdid_public_report_start(
    struct cdid_public_bundle scalar B,real colvector post_point,real colvector post_available,
    string colvector post_status,string colvector post_reason,real colvector aux_point,
    real colvector aux_supported,string colvector aux_reason,real scalar compute_if,
    struct cdid_state_limits scalar limits)
{
    struct cdid_public_report_buffer scalar out
    struct cdid_result_state scalar S
    struct cdid_population scalar P
    struct cdid_public_request scalar R
    struct cdid_agg_frame scalar F
    struct cdid_agg_shares scalar sh
    real scalar N,K,Jt,Q,j,q,h,u,slots,at_ts
    real colvector first_source,kpool,u_has
    string matrix cell_ts
    string scalar role,tid,support,why,pool,sel_str,sel_canon,jkey
    P=B.prepared.population;R=B.request;N=P.N;K=rows(R.report);Jt=rows(R.post);Q=rows(R.aux)
    if(rows(post_point)!=Jt | rows(post_available)!=Jt | rows(post_status)!=Jt | rows(post_reason)!=Jt |
        rows(aux_point)!=Q | rows(aux_supported)!=Q | rows(aux_reason)!=Q |
        (compute_if!=0 & compute_if!=1)) _cdid_pr_fail("cdid_native_reporting_input","Passed post/auxiliary table of cells quantities do not align.")
    if(missing(post_available) | any((post_available:!=0):&(post_available:!=1)) |
        missing(aux_supported) | any((aux_supported:!=0):&(aux_supported:!=1)))
        _cdid_pr_fail("cdid_native_reporting_input","Passed support flags must be explicit zero/one.")
    S=cdid_state_empty();S.id_kind=(P.ids.string_id ? "string" : "numeric")
    S.cluster_kind=(B.prepared.clustered ? B.prepared.cluster_kind : "none")
    S.ledger=R.report;S.ledger_columns=("g","time","ref","base_t","role_code","post_index","auxiliary_index","sign","request_status","selection_index")
    S.point=J(K,1,.);S.available=J(K,1,0);S.role=R.report[,5];S.reason=R.report_reason
    S.unit_numeric_id=(P.ids.string_id ? (1::N) : P.ids.numeric_id)
    S.unit_string_id=(P.ids.string_id ? P.ids.string_id_value : J(N,1,""))
    S.cohort=P.cohort;S.cohort_known=P.cohort_known;S.weight=B.prepared.weight.V
    first_source=J(N,1,.);u_has=_cdid_rb_index(P.first_row:>0)
    if(rows(u_has)) first_source[u_has]=B.prepared.source_row[P.first_row[u_has]]
    S.unit_metadata=(P.row_count,P.observed_count,P.A,P.B,P.frame_only,P.F_to_A,P.F_to_B,first_source,B.prepared.weight.owner_weight,B.prepared.weight.V)
    S.unit_columns=("supplied_rows","observed_rows","A","B","frame_only","F_to_A","F_to_B","first_source_row","owner_weight","first_calendar_weight")
    S.cluster_numeric_id=(B.prepared.clustered ? B.prepared.cluster_numeric_id : J(N,1,.))
    S.cluster_string_id=(B.prepared.clustered ? B.prepared.cluster_string_id : J(N,1,""))
    F.n=N;F.C=(B.prepared.clustered ? B.prepared.C : N);F.clustered=B.prepared.clustered
    F.cohort=P.cohort;F.cohort_known=P.cohort_known;F.w=B.prepared.weight.V
    F.unit_cluster=(B.prepared.clustered ? B.prepared.unit_cluster : (1::N))
    F.frame_id=(missing(B.prepared.weight.V) ? (P.ids.string_id ? P.ids.string_id_value : strofreal(P.ids.numeric_id,"%24.17g")) : J(0,1,""))
    F.share_law="first_calendar_masses"
    /* Restricted eligibility contexts: treated units outside their own cohort's
       context leave that cohort's mass (aggregation then describes the
       eligible treated population). */
    if(findexternal("__cdid_mref_state")!=NULL) if(cdid_mrpr_active()) {
        F.share_excluded=cdid_mrpr_share_excluded(P)
        S.unit_metadata=S.unit_metadata,F.share_excluded
        S.unit_columns=S.unit_columns,"reference_share_excluded"
    }
    sh=_cdid_agg_shares(F,R.target_cohorts)
    S.shares=(sh.cohorts,sh.p,sh.mass_complete,J(rows(R.target_cohorts),1,sh.shares_known))
    S.share_columns=("cohort","share","mass_complete","shares_known")
    S.share_IF=(compute_if & sh.shares_known ? sh.P : J(0,rows(R.target_cohorts),.))
    S.share_IF_columns=strofreal(R.target_cohorts',"%24.17g");S.calendar=P.calendar
    S.numeric_settings=(R.anticipation\R.pt_start\B.settings.alpha\B.settings.bstrap\B.settings.cband\B.settings.biters\compute_if\N\P.n_A\P.n_B)
    S.numeric_setting_names=("anticipation"\"pt_start"\"alpha"\"bstrap"\"cband"\"biters"\"compute_inffunc"\"N_F"\"N_A"\"N_B")
    sel_str=(B.prepared.selection=="" ? "none" : B.prepared.selection)
    sel_canon=(B.prepared.selection_canonical=="" ? "none" : B.prepared.selection_canonical)
    S.text_settings=("method",B.prepared.method\"gmm_weighting",B.settings.gmm_weighting\"base_period",R.base_period\
        "control_group",R.control_group\"pt_assumption",R.pt_assumption\"selection",sel_str\
        "selection_canonical",sel_canon\
        "selection_model",(B.prepared.selection_model=="" ? "pair_logit" : B.prepared.selection_model)\
        "selection_pooling",(B.prepared.selection_pooling=="" ? "separate" : B.prepared.selection_pooling)\
        "selection_universe",(B.prepared.selection_universe=="" ? "F" : B.prepared.selection_universe)\
        "selection_est_method",(B.prepared.selection_est_method=="" ? "ipw" : B.prepared.selection_est_method)\
        "imputation_pooling",B.prepared.imputation_pooling\
        "outcome_pooling",((B.prepared.method=="reg" | B.prepared.method=="dr") ? (B.prepared.outcome_pooling=="" ? "separate" : B.prepared.outcome_pooling) : "")\
        "selection_justification",B.prepared.selection_justification\
        "inactive_arguments",(B.prepared.selection_est_method=="dr" ? B.prepared.inactive_arguments :
            (B.prepared.inactive_arguments=="" ? ((B.prepared.method=="ipw" | B.prepared.method=="custom") ? "outcome_pooling imputation_pooling" : "imputation_pooling") : B.prepared.inactive_arguments))\
        "gmm_covariance",(B.prepared.gmm_covariance=="" ? "tikhonov" : B.prepared.gmm_covariance)\
        "design",B.prepared.design\
        "covariates",invtokens(B.prepared.x_names)\
        "selection_covariates",invtokens(B.prepared.sel_names)\
        "weight_source",B.prepared.weight.source\
        "weight_reduction",B.prepared.weight.reduction\"fix_weights_requested",B.prepared.weight.fix_weights\
        "aggregation_share_law","first_calendar_masses"\"share_status",(sh.shares_known ? "known" : sh.share_reason)\
        "IF_scale","original full-F frame influence; no additional rebase"\
        "role_convention","1=reference zero;2=direct unrestricted pre;3=post GMM target")
    if(B.selected_static_status!="")S.text_settings=S.text_settings\
        ("selected_static_gate",B.selected_static_status)\("selected_static_gate_reason",B.selected_static_reason)
    S.source_versions=("R_cdid","1.3.0.0"\
        "population","all units of the unit table"\"report","post-treatment, pre-treatment and reference cells")
    if(B.baseline_static_status!="")S.text_settings=S.text_settings\
        ("baseline_static_gate",B.baseline_static_status)\("baseline_static_gate_reason",B.baseline_static_reason)
    cell_ts=J(K*(3+(R.pt_tags_enabled!=0)),2,"");at_ts=0
    for(j=1;j<=K;j++) {
        role=(S.role[j]==3 ? "post" : (S.role[j]==1 ? "reference" : "pre"))
        tid=(role=="post" ? "ATT" : role)+"(g="+_cdid_pr_num(R.report[j,1])+",t="+_cdid_pr_num(R.report[j,2])+")"
        support="unavailable_pretrend";why=S.reason[j]
        if(S.role[j]==3) {
            h=R.report[j,6];S.available[j]=post_available[h];support=post_status[h];why=post_reason[h]
            if(S.available[j]) S.point[j]=post_point[h]
        }
        else if(S.role[j]==1) {
            S.available[j]=1;S.point[j]=0;support="deterministic_reference_zero"
            why="universal base-period reference cell: exact zero with zero influence function"
        }
        else if(R.report[j,9]==0) {
            q=R.report[j,7]
            if(q<1) _cdid_pr_fail("cdid_native_reporting_input","Pending direct pretrend lacks its auxiliary descriptor.")
            if(aux_supported[q]) {
                S.available[j]=1;S.point[j]=R.report[j,8]*aux_point[q];support="supported_direct_pretrend"
                kpool=R.aux_K_index[|R.aux_K_ptr[q]\R.aux_K_ptr[q+1]-1|];pool=""
                for(u=1;u<=rows(kpool);u++) pool=pool+(u>1 ? "," : "")+_cdid_pr_num(R.cohorts[kpool[u]])
                why="direct auxiliary comparison "+_cdid_pr_num(R.report[j,4])+"->"+_cdid_pr_num(R.report[j,2])+" with controls {"+pool+"}"
            }
            else why="unsupported_pair: "+aux_reason[q]
        }
        S.reason[j]=why
        jkey=_cdid_rb_key(j)
        at_ts++;cell_ts[at_ts,]=("cell_target_id:"+jkey,tid)
        at_ts++;cell_ts[at_ts,]=("cell_role:"+jkey,role)
        at_ts++;cell_ts[at_ts,]=("cell_support:"+jkey,support)
        if(R.pt_tags_enabled) {
            at_ts++;cell_ts[at_ts,]=("cell_pt_window_tag:"+jkey,R.pt_window_tag[j])
        }
    }
    if(rows(cell_ts)) S.text_settings=S.text_settings\cell_ts
    if(findexternal("__cdid_mref_state")!=NULL) cdid_mrpr_state_record(S)
    if(R.pt_tags_enabled) {
        S.ledger=S.ledger,R.pt_operator_covered;S.ledger_columns=S.ledger_columns,"pt_operator_covered"
    }
    if(any(S.available)) if(missing(S.point[_cdid_rb_index(S.available)]))
        _cdid_pr_fail("cdid_native_reporting_input","A supported fitted reporting point is not finite.")
    S.IF=J(0,K,.)
    slots=cdid_state_cost(S)[1]
    if(compute_if) {
        if(K>0) if(N>floor((limits.payload_slots-slots)/K)) _cdid_rs_fail("Full reported IF exceeds explicit payload limit.",908)
        S.IF=J(N,K,.)
    }
    cdid_state_validate(S,limits)
    out.payload=S;out.compute_if=compute_if;out.seen=J(N,1,0)
    return(out)
}
void cdid_public_report_add(struct cdid_public_report_buffer scalar out,real colvector units,
    real colvector post_ids,real matrix post_IF,real colvector aux_ids,real matrix aux_IF)
{
    real scalar K,j,h,q,p_max,a_max,pos
    real colvector post_pos,aux_pos
    real matrix block
    K=rows(out.payload.point)
    if(!out.compute_if) _cdid_pr_fail("cdid_native_reporting_input","Point-only reporting has no IF replay.")
    _cdid_rb_ids(units,rows(out.seen))
    if(any(out.seen[units])) _cdid_pr_fail("cdid_native_reporting_input","Reported full-F IF rows were already supplied.")
    if(rows(post_IF)!=rows(units) | cols(post_IF)!=rows(post_ids) | rows(aux_IF)!=rows(units) | cols(aux_IF)!=rows(aux_ids))
        _cdid_pr_fail("cdid_native_reporting_input","Passed IF chunks and descriptor identities do not align.")
    p_max=(rows(post_ids) ? max(post_ids) : 0);a_max=(rows(aux_ids) ? max(aux_ids) : 0)
    if(rows(post_ids)) _cdid_rb_ids(post_ids,p_max)
    if(rows(aux_ids)) _cdid_rb_ids(aux_ids,a_max)
    post_pos=J(p_max,1,0);if(rows(post_ids)) post_pos[post_ids]=(1::rows(post_ids))
    aux_pos=J(a_max,1,0);if(rows(aux_ids)) aux_pos[aux_ids]=(1::rows(aux_ids))
    block=J(rows(units),K,.)
    for(j=1;j<=K;j++) if(out.payload.available[j]) {
        if(out.payload.role[j]==1) block[,j]=J(rows(units),1,0)
        else {
            if(out.payload.role[j]==3) {
                h=out.payload.ledger[j,6]
                pos=(h>=1 & h<=p_max ? post_pos[h] : 0)
                if(pos<1) _cdid_pr_fail("cdid_native_reporting_input","A supported post target is missing from the replay chunk.")
                if(missing(post_IF[,pos])) _cdid_pr_fail("cdid_native_reporting_input","Supported reporting IF chunk contains missing values.")
                block[,j]=post_IF[,pos]
            }
            else {
                q=out.payload.ledger[j,7]
                pos=(q>=1 & q<=a_max ? aux_pos[q] : 0)
                if(pos<1) _cdid_pr_fail("cdid_native_reporting_input","A supported direct pretrend is missing from the replay chunk.")
                if(missing(aux_IF[,pos])) _cdid_pr_fail("cdid_native_reporting_input","Supported reporting IF chunk contains missing values.")
                block[,j]=out.payload.ledger[j,8]*aux_IF[,pos]
            }
        }
    }
    out.payload.IF[units,]=block;out.seen[units]=J(rows(units),1,1)
}
struct cdid_result_state scalar cdid_public_report_finish(struct cdid_public_report_buffer scalar out,
    struct cdid_state_limits scalar limits)
{
    if(out.compute_if) if(any(out.seen:!=1)) _cdid_pr_fail("cdid_native_reporting_input","Full original-F reporting IF replay is incomplete.")
    cdid_state_validate(out.payload,limits)
    return(out.payload)
}

/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
transmorphic cdid_mrpr_active()
{
    return(cdid_mref_active())
}

transmorphic cdid_mrpr_share_excluded(transmorphic a1)
{
    return(cdid_mref_share_excluded(a1))
}

void cdid_mrpr_state_record(transmorphic a1)
{
    cdid_mref_state_record(a1)
}

end
