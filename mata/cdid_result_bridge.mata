*! Private native result/inference/aggregation composition; current StataNow19.5.
*! Requires cdid_result_state, cdid_inference, cdid_aggregation modules.
version 19.5
mata:
struct cdid_bridge_request {
    struct cdid_infer_budget scalar budget
    struct cdid_state_limits scalar limits
    real scalar alpha, compute_inffunc, bstrap, cband, B, C, clustered, run_wald
    real colvector unit_cluster, report_cells, pre_cells
    string scalar posting_mode
}
struct cdid_bridge_result {
    struct cdid_result_state scalar payload
    real matrix b, V, Omega, se, analytic_se, ci, point_ci, coef_map
    real scalar n_F, C, K, point_only, critical, normal_critical
    real scalar wald_df, wald_statistic, wald_p, wald_C
    real colvector pre_cells
    string scalar kind, status, inference_status, flag_law, warning, wald_status, wald_reason
}
void _cdid_rb_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_rb_index(real matrix x)
{
    return(vec(selectindex(x)))
}
string scalar _cdid_rb_key(real scalar i)
{
    return(strtrim(strofreal(i,"%21.0f")))
}
void _cdid_rb_ids(real colvector ids, real scalar K)
{
    real scalar n
    n=rows(ids)
    if(cols(ids)!=1 | missing(ids) | any(ids:<1) | any(ids:>K) | any(ids:!=floor(ids)))
        _cdid_rb_fail("Requested cells must be explicit unique table of cells positions.",498)
    if(n>1) {
        if(!all(ids[|2\n|]:>ids[|1\n-1|])) {
            if(rows(uniqrows(sort(ids,1)))!=n)
                _cdid_rb_fail("Requested cell positions must not repeat.",498)
        }
    }
}
// cdid_bridge.* is the private bridge's owned result-metadata namespace.
real colvector _cdid_rb_keep_metadata(string colvector names, real scalar keep_cells)
{
    if(rows(names)==0) return(J(0,1,.))
    return(_cdid_rb_index((substr(names,1,12):!="cdid_bridge.") :|
        (keep_cells :& (substr(names,1,18):=="cdid_bridge.cells."))))
}
void _cdid_rb_clear_metadata(struct cdid_result_state scalar P, real scalar keep_cells)
{
    real colvector ix
    ix=_cdid_rb_keep_metadata(P.numeric_setting_names,keep_cells)
    P.numeric_setting_names=P.numeric_setting_names[ix]; P.numeric_settings=P.numeric_settings[ix]
    ix=_cdid_rb_keep_metadata(P.text_settings[,1],keep_cells)
    P.text_settings=P.text_settings[ix,]
    ix=_cdid_rb_keep_metadata(P.unit_columns',keep_cells)
    P.unit_columns=P.unit_columns[ix']; P.unit_metadata=P.unit_metadata[,ix']
    ix=_cdid_rb_keep_metadata(P.ledger_columns',keep_cells)
    P.ledger_columns=P.ledger_columns[ix']; P.ledger=P.ledger[,ix']
}
void _cdid_rb_text(struct cdid_result_state scalar P, string scalar key, string scalar value)
{
    P.text_settings=P.text_settings\("cdid_bridge."+key,value)
}
void _cdid_rb_number(struct cdid_result_state scalar P, string scalar key, real scalar value)
{
    P.numeric_setting_names=P.numeric_setting_names\("cdid_bridge."+key)
    P.numeric_settings=P.numeric_settings\value
}
// Fixed bridge outputs use labeled numeric settings, not a new file format.
void _cdid_rb_matrix(struct cdid_result_state scalar P, string scalar key, real matrix X)
{
    real scalar i, j, at
    string colvector names
    real colvector values
    names=J(length(X)+2,1,""); values=J(length(X)+2,1,.)
    names[1]="cdid_bridge."+key+".rows"; values[1]=rows(X)
    names[2]="cdid_bridge."+key+".cols"; values[2]=cols(X); at=2
    for(j=1;j<=cols(X);j++) for(i=1;i<=rows(X);i++) {
        at=at+1; names[at]="cdid_bridge."+key+"."+_cdid_rb_key(i)+"."+_cdid_rb_key(j)
        values[at]=X[i,j]
    }
    P.numeric_setting_names=P.numeric_setting_names\names
    P.numeric_settings=P.numeric_settings\values
}
void _cdid_rb_unit_matrix(struct cdid_result_state scalar P, string scalar key, real matrix X)
{
    real scalar j
    string rowvector names
    if(rows(X)==0) { _cdid_rb_matrix(P,key+".absent",X); return; }
    if(rows(X)!=rows(P.cohort)) _cdid_rb_fail("Derived unit rows lost original-F alignment.",498)
    names=J(1,cols(X),"")
    for(j=1;j<=cols(X);j++) names[j]="cdid_bridge."+key+"."+_cdid_rb_key(j)
    P.unit_metadata=P.unit_metadata,X; P.unit_columns=P.unit_columns,names
}
real scalar _cdid_rb_column(string rowvector names, string scalar wanted)
{
    real colvector ix
    ix=_cdid_rb_index(names:==wanted)
    if(rows(ix)!=1) _cdid_rb_fail("Stored payload requires exactly one named column: "+wanted,498)
    return(ix[1])
}
string scalar _cdid_rb_stored_text(struct cdid_result_state scalar P, string scalar key)
{
    real colvector ix
    ix=_cdid_rb_index(P.text_settings[,1]:==key)
    if(rows(ix)!=1) _cdid_rb_fail("Stored payload requires exactly one text setting: "+key,498)
    return(P.text_settings[ix[1],2])
}
real scalar _cdid_rb_stored_number(struct cdid_result_state scalar P, string scalar key)
{
    real colvector ix
    ix=_cdid_rb_index(P.numeric_setting_names:=="cdid_bridge."+key)
    if(rows(ix)!=1) _cdid_rb_fail("Stored bridge result lacks a unique numeric field: "+key,498)
    return(P.numeric_settings[ix[1]])
}
void _cdid_rb_preserve_wald(struct cdid_bridge_result scalar O, struct cdid_result_state scalar P)
{
    real colvector ix
    real scalar i, n
    O.wald_df=O.wald_statistic=O.wald_p=O.wald_C=.
    O.pre_cells=J(0,1,.)
    O.wald_status="not_computed"; O.wald_reason="source state has no stored full-family pretest"
    ix=_cdid_rb_index(P.text_settings[,1]:=="cdid_bridge.wald_status")
    if(rows(ix)==0) return
    O.wald_status=_cdid_rb_stored_text(P,"cdid_bridge.wald_status")
    O.wald_reason=_cdid_rb_stored_text(P,"cdid_bridge.wald_reason")
    O.wald_df=_cdid_rb_stored_number(P,"wald_df"); O.wald_statistic=_cdid_rb_stored_number(P,"wald_statistic")
    O.wald_p=_cdid_rb_stored_number(P,"wald_p"); O.wald_C=_cdid_rb_stored_number(P,"wald_C")
    n=_cdid_rb_stored_number(P,"pre_cells.rows"); O.pre_cells=J(n,1,.)
    for(i=1;i<=n;i++) O.pre_cells[i]=_cdid_rb_stored_number(P,"pre_cells."+_cdid_rb_key(i)+".1")
}
void _cdid_rb_request_check(struct cdid_result_state scalar P, struct cdid_bridge_request scalar R)
{
    real scalar N, i
    real colvector o
    struct cdid_infer_flags scalar checked_flags
    N=rows(P.cohort)
    cdid_state_validate(P,R.limits)
    if(R.alpha<=0 | R.alpha>=1 | R.alpha>=. | (R.compute_inffunc!=0 & R.compute_inffunc!=1) |
      (R.clustered!=0 & R.clustered!=1) | (R.run_wald!=0 & R.run_wald!=1) |
      (R.posting_mode!="matrix" & R.posting_mode!="state"))
        _cdid_rb_fail("Require explicit valid alpha, IF/Wald/cluster flags and posting mode.",498)
    if(R.compute_inffunc) checked_flags=cdid_infer_resolve_flags(R.bstrap,R.cband,"cells")
    _cdid_rb_ids(R.report_cells,rows(P.point)); _cdid_rb_ids(R.pre_cells,rows(P.point))
    if(R.C<1 | R.C>=. | R.C!=floor(R.C) | rows(R.unit_cluster)!=N | cols(R.unit_cluster)!=1 |
      missing(R.unit_cluster) | any(R.unit_cluster:<1) | any(R.unit_cluster:>R.C) |
      any(R.unit_cluster:!=floor(R.unit_cluster)))
        _cdid_rb_fail("Supply declared cluster positions for every original-F unit.",498)
    if(!R.clustered) {
        if(R.C!=N | any(R.unit_cluster:!=(1::N)))
            _cdid_rb_fail("Unclustered reporting requires each original unit in its own position.",498)
    }
    else {
        if(P.cluster_kind=="none") _cdid_rb_fail("Clustered reporting requires stored full-unit cluster IDs.",498)
        o=(P.cluster_kind=="numeric" ? order(P.cluster_numeric_id,1) : order(P.cluster_string_id,1))
        for(i=2;i<=N;i++) {
            if(P.cluster_kind=="numeric") {
                if(P.cluster_numeric_id[o[i]]==P.cluster_numeric_id[o[i-1]] & R.unit_cluster[o[i]]!=R.unit_cluster[o[i-1]])
                    _cdid_rb_fail("One stored cluster ID was split across declared positions.",498)
            } else if(P.cluster_string_id[o[i]]==P.cluster_string_id[o[i-1]] & R.unit_cluster[o[i]]!=R.unit_cluster[o[i-1]])
                _cdid_rb_fail("One stored cluster ID was split across declared positions.",498)
        }
        o=order(R.unit_cluster,1)
        for(i=2;i<=N;i++) if(R.unit_cluster[o[i]]==R.unit_cluster[o[i-1]]) {
            if(P.cluster_kind=="numeric") {
                if(P.cluster_numeric_id[o[i]]!=P.cluster_numeric_id[o[i-1]])
                    _cdid_rb_fail("Different stored cluster IDs were merged into one position.",498)
            } else if(P.cluster_string_id[o[i]]!=P.cluster_string_id[o[i-1]])
                _cdid_rb_fail("Different stored cluster IDs were merged into one position.",498)
        }
    }
    if(R.compute_inffunc & rows(P.IF)!=N) _cdid_rb_fail("Requested inference needs stored original-F influence functions.",459)
}
struct cdid_if_summary scalar _cdid_rb_summary(struct cdid_result_state scalar P,
    struct cdid_bridge_request scalar R, real colvector av)
{
    struct cdid_if_state scalar s
    real scalar step, lo, hi
    s=cdid_if_start(rows(P.cohort),R.unit_cluster,R.C,rows(av),R.clustered,R.budget)
    step=min((R.budget.max_if_batch_rows,floor(R.budget.max_if_batch_cells/rows(av))))
    if(step<1) _cdid_rb_fail("IF batch budget admits no complete original-unit row.",908)
    for(lo=1;lo<=rows(P.cohort);lo=lo+step) {
        hi=min((rows(P.cohort),lo+step-1))
        cdid_if_add(s,(lo::hi),P.IF[(lo..hi)',av'])
    }
    return(cdid_if_finish(s))
}
struct cdid_if_summary scalar _cdid_rb_subsummary(struct cdid_if_summary scalar S, real colvector ix)
{
    struct cdid_if_summary scalar T
    T=S; T.J=rows(ix)
    T.analytic_se=S.analytic_se[ix]; T.absolute_if_sum=S.absolute_if_sum[ix']
    T.cluster_mean=S.cluster_mean[ix']; T.raw_sums=S.raw_sums[,ix']
    T.centered_sums=S.centered_sums[,ix']; T.covariance=S.covariance[ix,ix']
    T.Omega=S.Omega[ix,ix']
    return(T)
}
struct cdid_infer_result scalar _cdid_rb_inference(real colvector theta,
    struct cdid_if_summary scalar S, struct cdid_bridge_request scalar R, real matrix signs)
{
    struct cdid_infer_flags scalar flags
    struct cdid_multiplier_state scalar ms
    struct cdid_multiplier_result scalar mult
    struct cdid_infer_result scalar out
    real colvector positive
    real scalar lo, hi, step
    flags=cdid_infer_resolve_flags(R.bstrap,R.cband,"cells")
    positive=_cdid_rb_index(S.analytic_se:>0); mult=cdid_multiplier_none()
    if(flags.draws & rows(positive)) {
        if(rows(signs)!=R.B | cols(signs)!=R.C) _cdid_rb_fail("Supply the complete declared draw-by-cluster sign matrix.",498)
        ms=cdid_multiplier_start(S,positive,R.B,R.alpha,R.budget)
        step=min((R.budget.max_sign_batch_rows,floor(R.budget.max_sign_batch_cells/R.C)))
        if(step<1) _cdid_rb_fail("Sign batch budget admits no draw row.",908)
        for(lo=1;lo<=R.B;lo=lo+step) {
            hi=min((R.B,lo+step-1)); cdid_multiplier_add(ms,(lo::hi),(1..R.C),signs[|lo,1\hi,R.C|])
        }
        mult=cdid_multiplier_finish(ms,0)
    }
    out=cdid_infer_report(theta,S,flags,mult,R.alpha)
    if(out.status!="computed") _cdid_rb_fail("Requested inference refused: "+out.reason,459)
    return(out)
}
void _cdid_rb_common_metadata(struct cdid_bridge_result scalar O, struct cdid_bridge_request scalar R)
{
    _cdid_rb_text(O.payload,"kind",O.kind); _cdid_rb_text(O.payload,"status",O.status)
    _cdid_rb_text(O.payload,"inference_status",O.inference_status)
    _cdid_rb_text(O.payload,"flag_law",O.flag_law); _cdid_rb_text(O.payload,"warning",O.warning)
    _cdid_rb_text(O.payload,"posting_mode",R.posting_mode)
    _cdid_rb_text(O.payload,"wald_status",O.wald_status); _cdid_rb_text(O.payload,"wald_reason",O.wald_reason)
    _cdid_rb_number(O.payload,"clustered",R.clustered); _cdid_rb_matrix(O.payload,"report_cells",R.report_cells)
    _cdid_rb_number(O.payload,"N",O.n_F); _cdid_rb_number(O.payload,"C",O.C)
    _cdid_rb_number(O.payload,"alpha",R.alpha); _cdid_rb_number(O.payload,"critical",O.critical)
    _cdid_rb_number(O.payload,"normal_critical",O.normal_critical)
    _cdid_rb_number(O.payload,"wald_df",O.wald_df); _cdid_rb_number(O.payload,"wald_statistic",O.wald_statistic)
    _cdid_rb_number(O.payload,"wald_C",O.wald_C); _cdid_rb_matrix(O.payload,"pre_cells",O.pre_cells)
    _cdid_rb_number(O.payload,"bstrap",R.bstrap); _cdid_rb_number(O.payload,"cband",R.cband); _cdid_rb_number(O.payload,"B",R.B)
    _cdid_rb_number(O.payload,"wald_p",O.wald_p); _cdid_rb_number(O.payload,"point_only",O.point_only)
    _cdid_rb_matrix(O.payload,"coef_map",O.coef_map)
    _cdid_rb_unit_matrix(O.payload,"cluster_position",R.unit_cluster)
    if(!any(O.payload.source_versions[,1]:=="cdid_private_result_bridge"))
        O.payload.source_versions=O.payload.source_versions\("cdid_private_result_bridge","1")
}
struct cdid_bridge_result scalar cdid_bridge_cells(string scalar token,
    struct cdid_bridge_request scalar R, real matrix signs)
{
    struct cdid_bridge_result scalar O
    struct cdid_result_state scalar P
    struct cdid_if_summary scalar S, sub
    struct cdid_infer_result scalar I
    struct cdid_wald_result scalar W
    real colvector av, report, pos, hit
    real matrix fields
    string rowvector names
    real scalar j, K, retained
    P=cdid_state_get(token); _cdid_rb_request_check(P,R)
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
            hit=J(K,1,0); hit[av]=(1::rows(av)); pos=hit[report]
            sub=_cdid_rb_subsummary(S,pos); I=_cdid_rb_inference(P.point[report],sub,R,signs)
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
    O.inference_status=O.status
    _cdid_rb_common_metadata(O,R); cdid_state_validate(O.payload,R.limits)
    return(O)
}
struct cdid_agg_ledger scalar _cdid_rb_agg_ledger(struct cdid_result_state scalar P)
{
    struct cdid_agg_ledger scalar L
    real scalar j, K
    K=rows(P.point)
    L.g=P.ledger[,_cdid_rb_column(P.ledger_columns,"g")]
    L.time=P.ledger[,_cdid_rb_column(P.ledger_columns,"time")]
    L.ref=P.ledger[,_cdid_rb_column(P.ledger_columns,"ref")]
    L.base_t=P.ledger[,_cdid_rb_column(P.ledger_columns,"base_t")]
    L.target_id=L.role=L.support_status=J(K,1,"")
    for(j=1;j<=K;j++) {
        L.target_id[j]=_cdid_rb_stored_text(P,"cell_target_id:"+_cdid_rb_key(j))
        L.role[j]=_cdid_rb_stored_text(P,"cell_role:"+_cdid_rb_key(j))
        L.support_status[j]=_cdid_rb_stored_text(P,"cell_support:"+_cdid_rb_key(j))
    }
    return(L)
}
struct cdid_bridge_result scalar cdid_bridge_aggregate(string scalar token,
    struct cdid_agg_request scalar Q, struct cdid_bridge_request scalar R, real matrix signs)
{
    struct cdid_bridge_result scalar O
    struct cdid_result_state scalar P
    struct cdid_agg_ledger scalar L
    struct cdid_agg_frame scalar F
    struct cdid_agg_result scalar A
    real colvector av, point_flag
    real scalar j
    P=cdid_state_get(token); _cdid_rb_request_check(P,R)
    point_flag=_cdid_rb_index(P.numeric_setting_names:=="cdid_bridge.point_only")
    if(rows(point_flag)) if(_cdid_rb_stored_number(P,"point_only")==1)
        _cdid_rb_fail("estat aggregations need influence functions; refit cdid without noinffunc.",459)
    if(!R.compute_inffunc | rows(P.IF)==0) _cdid_rb_fail("Point-only state cannot be aggregated: original influence functions are required.",459)
    if(R.run_wald) _cdid_rb_fail("Aggregation preserves stored full-family pretests; request a new Wald through cell reporting.",498)
    if(R.alpha!=Q.alpha | R.bstrap!=Q.bstrap | R.cband!=Q.cband | R.B!=Q.B)
        _cdid_rb_fail("Aggregation request and passed inference settings disagree.",498)
    if(rows(Q.calendar)!=rows(P.calendar) | any(Q.calendar:!=P.calendar))
        _cdid_rb_fail("Aggregation calendar must be the stored original calendar.",498)
    L=_cdid_rb_agg_ledger(P); F.n=rows(P.cohort); F.C=R.C; F.clustered=R.clustered
    F.cohort=P.cohort; F.cohort_known=P.cohort_known; F.unit_cluster=R.unit_cluster
    F.frame_id=J(F.n,1,""); for(j=1;j<=F.n;j++) F.frame_id[j]="F_row_"+_cdid_rb_key(j)
    if(P.id_kind=="string") F.frame_id=P.unit_string_id
    F.share_law=_cdid_rb_stored_text(P,"aggregation_share_law"); F.w=P.weight
    if(F.share_law=="first_calendar_masses") F.w=P.unit_metadata[,_cdid_rb_column(P.unit_columns,"first_calendar_weight")]
    if(any(P.unit_columns:=="reference_share_excluded")) F.share_excluded=cdid_mrrb_state_excluded(P)
    A=cdid_aggregate(P.point,P.IF,P.available,L,F,Q,signs,R.budget)
    if(A.inference.inference_status=="unavailable" & A.inference.n_available>0)
        _cdid_rb_fail("Requested aggregate inference refused: "+A.inference.reason,459)
    if(A.overall.inference.inference_status=="unavailable" & A.overall.inference.n_available>0)
        _cdid_rb_fail("Requested overall inference refused: "+A.overall.inference.reason,459)
    O.payload=P; _cdid_rb_clear_metadata(O.payload,1)
    O.kind="aggregation"; O.n_F=F.n; O.C=F.C; O.K=A.M; O.point_only=0
    O.status=A.status; O.inference_status=A.inference.inference_status; O.flag_law=A.inference.flag_law; O.warning=A.inference.warning
    O.normal_critical=A.inference.normal_critical; O.critical=A.inference.band_critical
    _cdid_rb_preserve_wald(O,P)
    av=_cdid_rb_index(A.cell_available:==1); O.b=A.estimates[av]'; O.se=A.inference.se[av]'
    O.analytic_se=A.inference.analytic_se[av]'; O.V=A.inference.covariance[av,av']
    O.Omega=A.inference.Omega[av,av']; O.point_ci=(A.inference.ci_lower[av],A.inference.ci_upper[av])'
    O.ci=(A.inference.band_lower[av],A.inference.band_upper[av])'; O.coef_map=(rows(av) ? ((1::rows(av)),av) : J(0,2,.))
    _cdid_rb_text(O.payload,"agg.type",A.type); _cdid_rb_text(O.payload,"agg.target_definition",A.target_definition)
    _cdid_rb_text(O.payload,"agg.overall.status",A.overall.status); _cdid_rb_text(O.payload,"agg.overall.reason",A.overall.reason)
    _cdid_rb_text(O.payload,"agg.overall.inference_status",A.overall.inference.inference_status)
    _cdid_rb_text(O.payload,"agg.overall.inference_kind",A.overall.inference_kind)
    _cdid_rb_text(O.payload,"agg.inference_reason",A.inference.reason)
    _cdid_rb_text(O.payload,"agg.band_family",A.inference.band_family)
    _cdid_rb_number(O.payload,"agg.na_rm",Q.na_rm); _cdid_rb_number(O.payload,"agg.balance_e",Q.balance_e)
    _cdid_rb_number(O.payload,"agg.min_e",Q.min_e); _cdid_rb_number(O.payload,"agg.max_e",Q.max_e)
    _cdid_rb_matrix(O.payload,"agg.request_cohort_grid",Q.cohort_grid)
    _cdid_rb_matrix(O.payload,"agg.n_intended",A.n_intended); _cdid_rb_matrix(O.payload,"agg.n_used",A.n_used)
    _cdid_rb_matrix(O.payload,"agg.share_dependent",A.share_dependent); _cdid_rb_matrix(O.payload,"agg.share_denominator",A.share_denominator)
    _cdid_rb_matrix(O.payload,"agg.overall_weights",A.overall.weights)
    _cdid_rb_matrix(O.payload,"agg.overall_share_gradient",A.overall.share_gradient)
    _cdid_rb_matrix(O.payload,"agg.overall_intended",A.overall.intended)
    _cdid_rb_matrix(O.payload,"agg.overall_used",A.overall.used); _cdid_rb_matrix(O.payload,"agg.overall_excluded",A.overall.excluded)
    _cdid_rb_matrix(O.payload,"agg.overall_analytic_se",A.overall.inference.analytic_se)
    _cdid_rb_matrix(O.payload,"agg.overall_V_estimator",A.overall.inference.covariance)
    _cdid_rb_matrix(O.payload,"agg.overall_V_asymptotic",A.overall.inference.Omega)
    _cdid_rb_matrix(O.payload,"agg.fixed_se",A.fixed_inference.se)
    _cdid_rb_matrix(O.payload,"agg.Omega_fixed",A.decomposition.Omega_fixed)
    _cdid_rb_matrix(O.payload,"agg.Omega_share",A.decomposition.Omega_share)
    _cdid_rb_matrix(O.payload,"agg.Omega_cross",A.decomposition.Omega_cross)
    _cdid_rb_matrix(O.payload,"agg.Omega_total",A.decomposition.Omega_total)
    _cdid_rb_matrix(O.payload,"agg.Omega_theta_p",A.decomposition.Omega_theta_p)
    _cdid_rb_matrix(O.payload,"agg.points",A.estimates); _cdid_rb_matrix(O.payload,"agg.keys",A.keys)
    _cdid_rb_matrix(O.payload,"agg.available",A.cell_available); _cdid_rb_matrix(O.payload,"agg.restricted",A.restricted)
    _cdid_rb_matrix(O.payload,"agg.se",A.inference.se); _cdid_rb_matrix(O.payload,"agg.analytic_se",A.inference.analytic_se)
    _cdid_rb_matrix(O.payload,"agg.point_ci",(A.inference.ci_lower,A.inference.ci_upper))
    _cdid_rb_matrix(O.payload,"agg.band_ci",(A.inference.band_lower,A.inference.band_upper))
    _cdid_rb_matrix(O.payload,"agg.V_estimator",A.inference.covariance); _cdid_rb_matrix(O.payload,"agg.V_asymptotic",A.inference.Omega)
    _cdid_rb_matrix(O.payload,"agg.W",A.W); _cdid_rb_matrix(O.payload,"agg.D",A.D)
    _cdid_rb_matrix(O.payload,"agg.intended",A.intended); _cdid_rb_matrix(O.payload,"agg.used",A.used)
    _cdid_rb_matrix(O.payload,"agg.excluded",A.excluded)
    _cdid_rb_matrix(O.payload,"agg.overall",(A.overall.estimate,A.overall.se,A.overall.ci_lower,A.overall.ci_upper))
    _cdid_rb_matrix(O.payload,"agg.share_cohorts",A.shares.cohorts); _cdid_rb_matrix(O.payload,"agg.share_p",A.shares.p)
    _cdid_rb_unit_matrix(O.payload,"agg.IF_total",A.if_total); _cdid_rb_unit_matrix(O.payload,"agg.IF_fixed",A.if_fixed)
    _cdid_rb_unit_matrix(O.payload,"agg.IF_share",A.if_share); _cdid_rb_unit_matrix(O.payload,"agg.cohort_share_IF",A.shares.P)
    _cdid_rb_unit_matrix(O.payload,"agg.overall_IF",A.overall.if_full)
    for(j=1;j<=A.M;j++) {
        _cdid_rb_text(O.payload,"agg.label."+_cdid_rb_key(j),A.labels[j])
        _cdid_rb_text(O.payload,"agg.reason."+_cdid_rb_key(j),A.cell_reason[j])
        _cdid_rb_text(O.payload,"agg.role."+_cdid_rb_key(j),A.bin_role[j])
        _cdid_rb_text(O.payload,"agg.support."+_cdid_rb_key(j),A.bin_support_status[j])
    }
    if(R.posting_mode=="state") { O.b=J(1,0,.); O.V=O.Omega=J(0,0,.); }
    _cdid_rb_common_metadata(O,R); cdid_state_validate(O.payload,R.limits)
    return(O)
}
string scalar cdid_bridge_stage(string scalar token, string scalar request_name, string scalar signs_name,
    string scalar agg_name, string scalar bname, string scalar vname, string scalar oname,
    string scalar sename, string scalar aname, string scalar ciname, string scalar piname, string scalar mapname,
    real scalar matrix_limit)
{
    struct cdid_bridge_request scalar R
    struct cdid_agg_request scalar Q
    struct cdid_infer_flags scalar flags
    struct cdid_bridge_result scalar O
    real matrix signs
    real scalar needs_signs, j
    string matrix stripe
    pointer() scalar ptr
    ptr=findexternal(request_name); if(ptr==NULL) _cdid_rb_fail("Explicit Mata bridge request not found.",111)
    R=*ptr
    needs_signs=0
    if(R.compute_inffunc) {
        if(agg_name=="") flags=cdid_infer_resolve_flags(R.bstrap,R.cband,"cells")
        else {
            ptr=findexternal(agg_name); if(ptr==NULL) _cdid_rb_fail("Explicit aggregation request not found.",111)
            Q=*ptr; flags=cdid_infer_resolve_flags(Q.bstrap,Q.cband,(Q.type=="simple" ? "simple" : "aggte_bins"))
        }
        needs_signs=flags.draws
    }
    signs=J(0,0,.)
    if(needs_signs & signs_name!="") {
        ptr=findexternal(signs_name); if(ptr==NULL) _cdid_rb_fail("Required supplied signs not found.",111)
        signs=*ptr
    }
    O=(agg_name=="" ? cdid_bridge_cells(token,R,signs) : cdid_bridge_aggregate(token,Q,R,signs))
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
void cdid_bridge_save(string scalar token, string scalar path, string scalar request_name)
{
    struct cdid_bridge_request scalar R
    pointer() scalar ptr
    ptr=findexternal(request_name); if(ptr==NULL) _cdid_rb_fail("Bridge request disappeared before save.",111)
    R=*ptr; cdid_state_save(token,path,R.limits)
}

/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
transmorphic cdid_mrrb_state_excluded(transmorphic a1)
{
    return(cdid_mref_state_excluded(a1))
}

end
