*! Private bounded public cells inference composition, current StataNow19.5.
*! Consumes fixed bridge helpers/result type and inference/state kernels.
*! Fixed fixed-sign bridge entry stays unchanged; generated BxC signs are never retained.
version 19.5
mata:
struct cdid_infer_result scalar _cdid_pb_inference(real colvector theta,
    struct cdid_if_summary scalar S,struct cdid_bridge_request scalar R,
    string scalar sign_mode,real matrix fixed_signs)
{
    struct cdid_infer_flags scalar flags,analytic
    struct cdid_multiplier_state scalar ms
    struct cdid_multiplier_result scalar mult
    struct cdid_infer_result scalar out,check
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar batch
    real colvector positive
    real scalar entries,cap
    flags=cdid_infer_resolve_flags(R.bstrap,R.cband,"cells")
    mult=cdid_multiplier_none()
    // Source validates all marginal scales before any RNG use. A zero-scale
    // nonconstant member cannot be removed from the multiplier family.
    analytic=cdid_infer_resolve_flags(0,0,"cells")
    check=cdid_infer_report(theta,S,analytic,mult,R.alpha)
    if(check.status!="computed") _cdid_rb_fail("Requested inference refused: "+check.reason,459)
    positive=_cdid_rb_index(S.analytic_se:>0)
    if(flags.draws & rows(positive)) {
        // preflight includes bounded IF/sign/cluster/rawdraw scratch.
        // A private supplied fixed matrix additionally remains caller-retained.
        entries=_cdid_infer_preflight(R.budget,S.n_F,S.C,S.J,R.B)
        cap=min((R.budget.max_workspace_cells,floor(R.budget.max_workspace_bytes/8)))
        if(sign_mode=="fixed") entries=_cdid_infer_cells(rows(fixed_signs),cols(fixed_signs),cap,entries)
        ms=cdid_multiplier_start(S,positive,R.B,R.alpha,R.budget)
        stream=cdid_public_rng_start(R.B,R.C,R.budget,sign_mode,fixed_signs)
        while(stream.cursor<stream.B) {
            batch=cdid_public_rng_next(stream)
            cdid_multiplier_add(ms,batch.draw_positions,batch.cluster_positions,batch.signs)
        }
        cdid_public_rng_complete(stream)
        mult=cdid_multiplier_finish(ms,0)
    }
    out=cdid_infer_report(theta,S,flags,mult,R.alpha)
    if(out.status!="computed") _cdid_rb_fail("Requested inference refused: "+out.reason,459)
    return(out)
}
struct cdid_bridge_result scalar _cdid_public_bridge_cells_state(struct cdid_result_state scalar P,
    struct cdid_bridge_request scalar R, string scalar sign_mode, real matrix fixed_signs)
{
    struct cdid_bridge_result scalar O
    struct cdid_if_summary scalar S, sub
    struct cdid_infer_result scalar I
    struct cdid_wald_result scalar W
    real colvector av, report, pos, av_pos
    real matrix fields
    string rowvector names
    real scalar K, retained
    O.payload=P; _cdid_rb_clear_metadata(O.payload,0)
    O.n_F=rows(P.cohort); O.C=R.C; K=rows(P.point); O.K=K
    O.kind="cells"; O.point_only=!R.compute_inffunc; O.status="inference_not_computed"
    O.flag_law="compute_inffunc_false_point_estimates_only"; O.warning=""
    I=cdid_infer_point_only((rows(R.report_cells) ? P.point[R.report_cells] : J(0,1,.)),O.n_F,O.C,R.alpha)
    O.normal_critical=O.critical=I.normal_critical
    O.pre_cells=R.pre_cells; O.wald_C=R.C
    O.wald_df=rows(R.pre_cells); O.wald_statistic=O.wald_p=.
    O.wald_status="not_computed"; O.wald_reason="influence functions not requested"
    report=(rows(R.report_cells) ? _cdid_rb_index(P.available[R.report_cells]:==1) : J(0,1,.))
    report=(rows(report) ? R.report_cells[report] : J(0,1,.))
    if(rows(report)) if(missing(P.point[report])) _cdid_rb_fail("Available reporting points must be finite.",498)
    fields=J(K,8,.)
    O.payload.V_estimator=O.payload.V_asymptotic=J(0,0,.)
    if(R.compute_inffunc) {
        retained=cdid_state_cost(P)[1]-length(P.V_estimator)-length(P.V_asymptotic)
        if(K) if(K>floor((R.limits.payload_slots-retained)/(2*K)))
            _cdid_rb_fail("Full-table covariances exceed explicit payload storage limits.",908)
        av=_cdid_rb_index(P.available:==1)
        O.payload.V_estimator=O.payload.V_asymptotic=J(K,K,.)
        if(rows(av)) {
            S=_cdid_rb_summary(P,R,av)
            O.payload.V_estimator[av,av']=S.covariance; O.payload.V_asymptotic[av,av']=S.Omega
        }
        if(rows(report)) {
            av_pos=J(K,1,0); av_pos[av]=(1::rows(av)); pos=av_pos[report]
            sub=_cdid_rb_subsummary(S,pos); I=_cdid_pb_inference(P.point[report],sub,R,sign_mode,fixed_signs)
            fields[report,]=(I.se,I.analytic_se,I.ci_lower,I.ci_upper,I.band_lower,I.band_upper,I.z_stat,I.p_value)
            O.status=I.status; O.flag_law=I.flag_law; O.warning=I.warning
            O.critical=I.band_critical; O.normal_critical=I.normal_critical
        } else { O.status="no_supported_reporting_targets"; O.flag_law="no_supported_reporting_targets"; }
        if(R.run_wald) {
            if(rows(R.pre_cells)) W=cdid_wald_full(P.point[R.pre_cells],O.payload.V_estimator[R.pre_cells,R.pre_cells'],
                P.available[R.pre_cells],R.C,R.budget.max_workspace_cells)
            else W=cdid_wald_full(J(0,1,.),J(0,0,.),J(0,1,.),R.C,R.budget.max_workspace_cells)
            O.wald_df=W.df; O.wald_statistic=W.statistic; O.wald_p=W.p_value
            O.wald_status=W.status; O.wald_reason=W.reason
        }
    }
    names=("se","analytic_se","ci_lower","ci_upper","band_lower","band_upper","z","p")
    O.payload.ledger=O.payload.ledger,fields; O.payload.ledger_columns=O.payload.ledger_columns,"cdid_bridge.cells.":+names
    _cdid_rb_text(O.payload,"cells.status",O.status); _cdid_rb_text(O.payload,"cells.flag_law",O.flag_law)
    _cdid_rb_number(O.payload,"cells.C",R.C); _cdid_rb_number(O.payload,"cells.clustered",R.clustered)
    _cdid_rb_number(O.payload,"cells.alpha",R.alpha); _cdid_rb_number(O.payload,"cells.B",R.B)
    _cdid_rb_number(O.payload,"cells.bstrap",R.bstrap); _cdid_rb_number(O.payload,"cells.cband",R.cband)
    _cdid_rb_number(O.payload,"cells.critical",O.critical); _cdid_rb_number(O.payload,"cells.normal_critical",O.normal_critical)
    _cdid_rb_unit_matrix(O.payload,"cells.cluster_position",R.unit_cluster)
    _cdid_rb_matrix(O.payload,"cells.report_cells",R.report_cells)
    O.b=(rows(report) ? P.point[report]' : J(1,0,.)); O.se=(rows(report) ? fields[report,1]' : J(1,0,.)); O.analytic_se=(rows(report) ? fields[report,2]' : J(1,0,.))
    O.point_ci=(rows(report) ? fields[report,(3,4)]' : J(2,0,.)); O.ci=(rows(report) ? fields[report,(5,6)]' : J(2,0,.))
    O.V=O.Omega=J(0,0,.)
    if(R.compute_inffunc & rows(report)) { O.V=O.payload.V_estimator[report,report']; O.Omega=O.payload.V_asymptotic[report,report']; }
    O.coef_map=(rows(report) ? ((1::rows(report)),report) : J(0,2,.))
    if(R.posting_mode=="state") { O.b=J(1,0,.); O.V=O.Omega=J(0,0,.); }
    // The preceding clear removed prior bridge-derived metadata; append actual RNG use once.
    _cdid_rb_text(O.payload,"rng.mode",(I.draws_used ? (sign_mode=="caller" ? "caller_rng_stream_Rademacher" : "caller_fixed_signs") : "not_used"))
    _cdid_rb_number(O.payload,"rng.draws_used",I.draws_used)
    _cdid_rb_number(O.payload,"rng.n_draws",I.n_draws)
    if(I.draws_used) _cdid_rb_number(O.payload,"rng.batch_rows",min((R.B,R.budget.max_sign_batch_rows,floor(R.budget.max_sign_batch_cells/R.C))))
    O.inference_status=O.status
    _cdid_rb_common_metadata(O,R); cdid_state_validate(O.payload,R.limits)
    return(O)
}
struct cdid_bridge_result scalar cdid_public_bridge_cells_stream(string scalar token,
    struct cdid_bridge_request scalar R, string scalar sign_mode, real matrix fixed_signs)
{
    struct cdid_result_state scalar P
    P=cdid_state_get(token); _cdid_rb_request_check(P,R)
    return(_cdid_public_bridge_cells_state(P,R,sign_mode,fixed_signs))
}
string scalar cdid_public_bridge_stage(string scalar token, string scalar request_name, string scalar signs_name,
    string scalar sign_mode, string scalar bname, string scalar vname, string scalar oname,
    string scalar sename, string scalar aname, string scalar ciname, string scalar piname, string scalar mapname,
    real scalar matrix_limit)
{
    struct cdid_bridge_request scalar R
    struct cdid_infer_flags scalar flags
    struct cdid_bridge_result scalar O
    struct cdid_result_state scalar source
    real matrix signs
    real scalar needs_signs, j
    string matrix stripe
    pointer() scalar ptr
    ptr=findexternal(request_name); if(ptr==NULL) _cdid_rb_fail("Explicit Mata bridge request not found.",111)
    R=*ptr
    // Validate Stata posting limits before any generated statistical draw.
    if(matrix_limit<1 | matrix_limit>=. | matrix_limit!=floor(matrix_limit))
        _cdid_rb_fail("Current standard Stata matrix limit was not supplied.",498)
    if(sign_mode!="caller" & sign_mode!="fixed")
        _cdid_rb_fail("Public sign source must be caller or private fixed signs.",498)
    source=cdid_state_get(token)
    _cdid_rb_request_check(source,R)
    if(R.posting_mode=="matrix") if((rows(R.report_cells) ? sum(source.available[R.report_cells]:==1) : 0)>matrix_limit)
        _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.",908)
    needs_signs=0
    if(R.compute_inffunc) {
        flags=cdid_infer_resolve_flags(R.bstrap,R.cband,"cells")
        needs_signs=flags.draws
    }
    signs=J(0,0,.)
    if(needs_signs & sign_mode=="fixed" & signs_name!="") {
        ptr=findexternal(signs_name)
        if(ptr!=NULL) signs=*ptr
    }
    O=_cdid_public_bridge_cells_state(source,R,sign_mode,signs)
    if(matrix_limit<1 | matrix_limit>=. | matrix_limit!=floor(matrix_limit))
        _cdid_rb_fail("Current standard Stata matrix limit was not supplied.",498)
    if(cols(O.b)>matrix_limit)
        _cdid_rb_fail("Mapped coefficient vector exceeds standard Stata matrix-expression limit; use explicit state posting; full scientific results remain in state.",908)
    stripe=J(cols(O.b),2,"")
    for(j=1;j<=cols(O.b);j++) stripe[j,2]=(O.kind=="cells" ? "cell_" : "agg_")+_cdid_rb_key(O.coef_map[j,2])
    if(cols(O.b)) { st_matrix(bname,O.b); st_matrixcolstripe(bname,stripe); }
    if(rows(O.V)) {
        st_matrix(vname,O.V); st_matrixcolstripe(vname,stripe); st_matrixrowstripe(vname,stripe)
        st_matrix(oname,O.Omega); st_matrixcolstripe(oname,stripe); st_matrixrowstripe(oname,stripe)
    }
    if(cols(O.b)) {
        st_matrix(sename,O.se); st_matrix(aname,O.analytic_se); st_matrix(ciname,O.ci)
        st_matrix(piname,O.point_ci); st_matrix(mapname,O.coef_map)
    }
    st_local("bridge_npost",_cdid_rb_key(cols(O.b))); st_local("bridge_hasV",_cdid_rb_key(rows(O.V)>0))
    st_local("bridge_N_source",_cdid_rb_key(rows(O.payload.point)))
    st_local("bridge_alpha",strofreal(R.alpha,"%24.17g")); st_local("bridge_B",strofreal(R.B,"%24.17g"))
    st_local("bridge_bstrap",strofreal(R.bstrap,"%24.17g")); st_local("bridge_cband",_cdid_rb_key(R.cband))
    st_local("bridge_clustered",_cdid_rb_key(R.clustered))
    st_local("bridge_N",_cdid_rb_key(O.n_F)); st_local("bridge_C",_cdid_rb_key(O.C)); st_local("bridge_K",_cdid_rb_key(O.K))
    st_local("bridge_wald_C",strofreal(O.wald_C,"%24.17g"))
    st_local("bridge_posting_mode",R.posting_mode); st_local("bridge_kind",O.kind); st_local("bridge_status",O.status); st_local("bridge_inference_status",O.inference_status); st_local("bridge_flags",O.flag_law)
    st_local("bridge_warning",O.warning); st_local("bridge_wald_status",O.wald_status); st_local("bridge_wald_reason",O.wald_reason)
    st_local("bridge_critical",strofreal(O.critical,"%24.17g")); st_local("bridge_normal",strofreal(O.normal_critical,"%24.17g"))
    st_local("bridge_wald_df",strofreal(O.wald_df,"%24.17g")); st_local("bridge_wald_stat",strofreal(O.wald_statistic,"%24.17g"))
    st_local("bridge_wald_p",strofreal(O.wald_p,"%24.17g")); st_local("bridge_point_only",_cdid_rb_key(O.point_only))
    return(cdid_state_put(O.payload,R.limits))
}
end
