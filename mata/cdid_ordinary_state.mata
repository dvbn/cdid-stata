*! Ordinary payload and no-refit consumer adapter for result registry.
*! Source IF/unit/weight ordering and normalization must be explicit caller inputs.
version 19.5
mata:
string scalar _cdid_ord_key(real scalar i)
{
    return(strtrim(strofreal(i,"%21.0f")))
}
void _cdid_ord_meta_number(struct cdid_result_state scalar P,string scalar key,real scalar value)
{
    P.numeric_setting_names=P.numeric_setting_names\key;P.numeric_settings=P.numeric_settings\value
}
void _cdid_ord_meta_text(struct cdid_result_state scalar P,string scalar key,string scalar value)
{
    P.text_settings=P.text_settings\(key,value)
}
void _cdid_ord_meta_matrix(struct cdid_result_state scalar P,string scalar key,real matrix X)
{
    real scalar i,j
    _cdid_ord_meta_number(P,key+".rows",rows(X));_cdid_ord_meta_number(P,key+".cols",cols(X))
    for(j=1;j<=cols(X);j++) for(i=1;i<=rows(X);i++) _cdid_ord_meta_number(P,key+"."+_cdid_ord_key(i)+"."+_cdid_ord_key(j),X[i,j])
}
real scalar _cdid_ord_number(struct cdid_result_state scalar P,string scalar key)
{
    real colvector ix
    ix=_cdid_ord_ix(P.numeric_setting_names:==key)
    if(rows(ix)!=1) _cdid_ord_fail("Ordinary state requires one numeric field: "+key,498)
    return(P.numeric_settings[ix[1]])
}
string scalar _cdid_ord_text(struct cdid_result_state scalar P,string scalar key)
{
    real colvector ix
    ix=_cdid_ord_ix(P.text_settings[,1]:==key)
    if(rows(ix)!=1) _cdid_ord_fail("Ordinary state requires one text field: "+key,498)
    return(P.text_settings[ix[1],2])
}
/* Read the fit's cluster name without copying its retained influence functions. */
string scalar _cdid_ord_fit_cluster(string scalar token)
{
    struct cdid_state_registry scalar registry
    pointer(struct cdid_result_state scalar) scalar fit

    registry=_cdid_rs_registry()
    fit=registry.payloads[_cdid_rs_index(registry,token)]
    return(_cdid_ord_text(*fit,"ord.cluster_name"))
}
real matrix _cdid_ord_matrix(struct cdid_result_state scalar P,string scalar key)
{
    real matrix X
    real scalar i,j
    X=J(_cdid_ord_number(P,key+".rows"),_cdid_ord_number(P,key+".cols"),.)
    for(j=1;j<=cols(X);j++) for(i=1;i<=rows(X);i++) X[i,j]=_cdid_ord_number(P,key+"."+_cdid_ord_key(i)+"."+_cdid_ord_key(j))
    return(X)
}
real scalar _cdid_ord_column(string rowvector names,string scalar key)
{
    real colvector ix
    ix=_cdid_ord_ix(names:==key)
    if(rows(ix)!=1) _cdid_ord_fail("Ordinary state requires one unit column: "+key,498)
    return(ix[1])
}
void _cdid_ord_remove(struct cdid_result_state scalar P,string scalar prefix)
{
    real colvector keep
    keep=_cdid_ord_ix(substr(P.numeric_setting_names,1,strlen(prefix)):!=prefix)
    P.numeric_setting_names=P.numeric_setting_names[keep];P.numeric_settings=P.numeric_settings[keep]
    keep=_cdid_ord_ix(substr(P.text_settings[,1],1,strlen(prefix)):!=prefix);P.text_settings=P.text_settings[keep,]
    keep=_cdid_ord_ix(substr(P.unit_columns',1,strlen(prefix)):!=prefix)
    P.unit_columns=P.unit_columns[keep'];P.unit_metadata=P.unit_metadata[,keep']
}
void _cdid_ord_unit_matrix(struct cdid_result_state scalar P,string scalar key,real matrix X)
{
    real scalar j
    string rowvector names
    if(rows(X)!=rows(P.cohort)) _cdid_ord_fail("Ordinary derived IF lost original-unit alignment.",498)
    names=J(1,cols(X),"");for(j=1;j<=cols(X);j++) names[j]=key+"."+_cdid_ord_key(j)
    P.unit_metadata=P.unit_metadata,X;P.unit_columns=P.unit_columns,names
}
void cdid_ord_state_attach_cells(struct cdid_result_state scalar P,
    struct cdid_ord_inference scalar O,real scalar alpha,real scalar B)
{
    _cdid_ord_remove(P,"ord.report.");_cdid_ord_remove(P,"ord.cells.")
    _cdid_ord_meta_text(P,"ord.report.kind","ordinary_cells")
    _cdid_ord_meta_text(P,"ord.report.warning",O.warning)
    _cdid_ord_meta_text(P,"ord.report.covariance_origin","native raw Gram; source MP covariance separately labeled")
    _cdid_ord_meta_number(P,"ord.report.point_only",!O.compute_inffunc)
    _cdid_ord_meta_number(P,"ord.report.alpha",alpha);_cdid_ord_meta_number(P,"ord.report.B",B)
    _cdid_ord_meta_number(P,"ord.report.bstrap",O.bstrap);_cdid_ord_meta_number(P,"ord.report.cband",O.cband)
    _cdid_ord_meta_number(P,"ord.report.critical",O.critical);_cdid_ord_meta_number(P,"ord.report.normal",O.normal_critical)
    _cdid_ord_meta_number(P,"ord.report.calls",O.calls)
    _cdid_ord_meta_matrix(P,"ord.report.point",P.point);_cdid_ord_meta_matrix(P,"ord.report.labels",(1::rows(P.point)))
    _cdid_ord_meta_matrix(P,"ord.report.se",O.se);_cdid_ord_meta_matrix(P,"ord.report.ci",O.ci)
    _cdid_ord_meta_matrix(P,"ord.report.point_ci",O.point_ci)
    if(O.compute_inffunc) {
        P.V_estimator=O.summary.V;P.V_asymptotic=O.summary.Omega
        _cdid_ord_meta_matrix(P,"ord.report.analytic_se",O.summary.analytic_se)
        _cdid_ord_meta_matrix(P,"ord.report.V",O.summary.V);_cdid_ord_meta_matrix(P,"ord.report.Omega",O.summary.Omega)
        _cdid_ord_meta_matrix(P,"ord.cells.source_Omega",O.summary.source_Omega)
        _cdid_ord_meta_text(P,"ord.cells.source_covariance_method",O.summary.source_covariance_method)
        _cdid_ord_meta_matrix(P,"ord.cells.source_analytic_excluded",O.summary.source_excluded)
    }
    else {
        P.IF=J(0,rows(P.point),.);P.V_estimator=P.V_asymptotic=J(0,0,.)
        _cdid_ord_meta_matrix(P,"ord.report.analytic_se",J(rows(P.point),1,.))
        _cdid_ord_meta_matrix(P,"ord.report.V",J(0,0,.));_cdid_ord_meta_matrix(P,"ord.report.Omega",J(0,0,.))
        _cdid_ord_meta_matrix(P,"ord.cells.source_Omega",J(0,0,.))
    }
    _cdid_ord_meta_number(P,"ord.cells.compute_inffunc",O.compute_inffunc)
    _cdid_ord_meta_number(P,"ord.cells.alpha",alpha);_cdid_ord_meta_number(P,"ord.cells.B",B)
    _cdid_ord_meta_number(P,"ord.cells.bstrap",O.bstrap);_cdid_ord_meta_number(P,"ord.cells.cband",O.cband)
    _cdid_ord_meta_matrix(P,"ord.cells.se",O.se);_cdid_ord_meta_matrix(P,"ord.cells.ci",O.ci)
    _cdid_ord_meta_matrix(P,"ord.cells.bootstrap_helper_V",O.bootstrap.V)
    _cdid_ord_meta_matrix(P,"ord.cells.bootstrap_retained",O.bootstrap.retained)
    _cdid_ord_meta_text(P,"ord.cells.bootstrap_V_scale","cov(raw bres), B-1 denominator, retained targets only; absent when return_V=FALSE")
    _cdid_ord_meta_text(P,"ord.cells.wald_status",O.wald.status);_cdid_ord_meta_text(P,"ord.cells.wald_reason",O.wald.reason)
    _cdid_ord_meta_number(P,"ord.cells.wald_statistic",O.wald.statistic);_cdid_ord_meta_number(P,"ord.cells.wald_p",O.wald.p)
    _cdid_ord_meta_number(P,"ord.cells.wald_p_unrounded",O.wald.p_unrounded);_cdid_ord_meta_number(P,"ord.cells.wald_df",O.wald.df)
    _cdid_ord_meta_matrix(P,"ord.cells.wald_declared",O.wald.declared);_cdid_ord_meta_matrix(P,"ord.cells.wald_used",O.wald.used)
    _cdid_ord_meta_matrix(P,"ord.cells.wald_excluded",O.wald.excluded)
    _cdid_ord_meta_text(P,"ord.conditional_test","separate conditional CvM producer required; ordinary Wald is not conditional CvM")
}
struct cdid_result_state scalar cdid_ord_state_build(real matrix cells,real matrix F,
    real matrix identity_map,real colvector source_group,real colvector raw_weight,
    real colvector source_weight,real scalar normalization_factor,real colvector calendar,
    real colvector cluster_labels,real colvector cluster_positions,real scalar C,
    string scalar cluster_name,real scalar D1_exception,string matrix source_versions,
    struct cdid_ord_inference scalar O,real scalar alpha,real scalar B,
    struct cdid_state_limits scalar limits)
{
    struct cdid_result_state scalar P
    real scalar n,k,j
    real colvector groups
    n=rows(source_group);k=rows(cells)
    if(cols(cells)!=12 | cols(F)!=k | rows(identity_map)!=n | cols(identity_map)!=3 |
        rows(source_weight)!=n | rows(raw_weight)!=n | rows(cluster_positions)!=n | rows(cluster_labels)!=n |
        missing(source_weight) | any(source_weight:<0) | sum(source_weight)<=0 | missing(raw_weight) | any(raw_weight:<0) |
        normalization_factor<=0 | normalization_factor>=. | missing(source_group) |
        max(abs(source_weight-raw_weight*normalization_factor):/(1:+abs(source_weight)))>1e-10 |
        (D1_exception!=0 & D1_exception!=1))
        _cdid_ord_fail("Explicit ordinary source identity/weight normalization metadata is invalid.",498)
    P=cdid_state_empty();P.id_kind="numeric";P.cluster_kind=(cluster_name=="" ? "none" : "numeric")
    P.ledger=cells;P.ledger_columns=("group","time","base","role","available","ATT","source_analytic_SE","event","N","N_treated","N_control","N_other")
    P.point=cells[,6];P.available=cells[,5];P.role=cells[,4];P.IF=F
    P.reason=J(k,1,"");for(j=1;j<=k;j++) if(!P.available[j]) P.reason[j]="source ATT unavailable"
    P.unit_numeric_id=identity_map[,1];P.unit_string_id=J(n,1,"")
    P.cohort=source_group;P.cohort_known=J(n,1,1);P.weight=source_weight
    P.unit_metadata=(identity_map[,2..3],raw_weight,source_weight,cluster_positions)
    P.unit_columns=("source_row","input_id","ordinary_raw_unit_weight","ordinary_source_weight","ordinary_cluster_position")
    P.cluster_numeric_id=cluster_labels;P.cluster_string_id=J(n,1,"");P.calendar=calendar;P.source_versions=source_versions
    groups=uniqrows(sort(cells[,1],1));P.shares=(groups,J(rows(groups),1,.));P.share_columns=("cohort","source_mass")
    P.share_IF=J(n,rows(groups),.);P.share_IF_columns=J(1,rows(groups),"")
    for(j=1;j<=rows(groups);j++) {
        P.shares[j,2]=mean(source_weight:*(source_group:==groups[j]))
        P.share_IF[,j]=source_weight:*(source_group:==groups[j]):-P.shares[j,2]
        P.share_IF_columns[j]="mass_IF_"+_cdid_ord_key(j)
    }
    _cdid_ord_meta_text(P,"ord.mode","ordinary")
    _cdid_ord_meta_text(P,"ord.source_weight_law","source first-period/time-invariant or mean-per-internal-ID weights, source IF order")
    _cdid_ord_meta_text(P,"ord.cluster_name",cluster_name);_cdid_ord_meta_number(P,"ord.C",C)
    _cdid_ord_meta_number(P,"ord.weight_normalization_factor",normalization_factor)
    _cdid_ord_meta_number(P,"ord.D1_exception",D1_exception)
    cdid_ord_state_attach_cells(P,O,alpha,B);cdid_state_validate(P,limits)
    return(P)
}
/* Original cells and IF always remain the authoritative base payload. Active
   aggregate outputs live in an ordinary-owned namespace, preserving re-aggregation. */
struct cdid_result_state scalar cdid_ord_state_aggregate(struct cdid_result_state scalar P,
    struct cdid_ord_agg_request scalar R,string scalar requested_cluster,
    real scalar row_chunk,transmorphic scalar calls,struct cdid_state_limits scalar limits)
{
    struct cdid_ord_aggregate_result scalar A
    struct cdid_ord_summary scalar summary
    real colvector cluster,point,se,labels
    real scalar C,clustered,j,k,normal
    string scalar stored,warning
    real matrix allIF,V,ci,pointci
    cdid_state_validate(P,limits)
    if(_cdid_ord_text(P,"ord.mode")!="ordinary") _cdid_ord_fail("Ordinary adapter requires an ordinary payload.",498)
    if(rows(P.IF)==0 | !_cdid_ord_number(P,"ord.cells.compute_inffunc"))
        _cdid_ord_fail("Ordinary point-only state cannot be aggregated.",459)
    if(R.bstrap==.) R.bstrap=_cdid_ord_number(P,"ord.cells.bstrap")
    if(R.cband==.) R.cband=_cdid_ord_number(P,"ord.cells.cband")
    if(R.B==.) R.B=_cdid_ord_number(P,"ord.cells.B")
    if(R.alpha==.) R.alpha=_cdid_ord_number(P,"ord.cells.alpha")
    stored=_cdid_ord_text(P,"ord.cluster_name");warning=""
    if(requested_cluster=="inherit") requested_cluster=stored
    clustered=(requested_cluster!="" & requested_cluster==stored)
    if(requested_cluster!="" & !clustered) {
        warning="Warning: cluster("+requested_cluster+") is unavailable after this fit. "
        if(stored=="") warning=warning+"The fit did not retain a cluster variable. "
        else warning=warning+"The fit retained clustering on "+stored+". "
        warning=warning+"Reporting unit-level standard errors. Refit cdid with cluster("+requested_cluster+") for clustered inference."
    }
    cluster=(clustered ? P.unit_metadata[,_cdid_ord_column(P.unit_columns,"ordinary_cluster_position")] : (1::rows(P.cohort)))
    C=(clustered ? _cdid_ord_number(P,"ord.C") : rows(P.cohort))
    A=cdid_ord_aggregate_points(P.IF,P.point,P.ledger[,1],P.ledger[,2],P.cohort,P.weight,P.calendar,R)
    if(A.status!="computed") _cdid_ord_fail(A.reason,459)
    A=cdid_ord_aggregate_report(A,R,cluster,C,clustered,row_chunk,calls)
    _cdid_ord_remove(P,"ord.report.");_cdid_ord_remove(P,"ord.agg.")
    k=rows(A.point);normal=invnormal(1-R.alpha/2)
    point=A.point\A.overall_point;se=A.se\A.overall_se
    labels=A.egt\.;allIF=(A.IF,A.overall_IF)
    /* Full joint bin+overall covariance is native IF-derived, not source AGGTE. */
    summary=cdid_ord_if_summary(allIF,cluster,C,clustered,R.bstrap,row_chunk);V=summary.V
    ci=(point:-A.critical*se,point:+A.critical*se);pointci=(point:-normal*se,point:+normal*se)
    if(k) ci[k+1,]=pointci[k+1,]
    if(warning!="" & A.warning!="" & substr(A.warning,1,1)!=";") warning=warning+"; "
    _cdid_ord_meta_text(P,"ord.report.kind","ordinary_"+A.type);_cdid_ord_meta_text(P,"ord.report.warning",warning+A.warning)
    _cdid_ord_meta_text(P,"ord.report.covariance_origin",A.covariance_origin)
    _cdid_ord_meta_number(P,"ord.report.point_only",0);_cdid_ord_meta_number(P,"ord.report.alpha",R.alpha)
    _cdid_ord_meta_number(P,"ord.report.B",R.B);_cdid_ord_meta_number(P,"ord.report.bstrap",A.bstrap);_cdid_ord_meta_number(P,"ord.report.cband",A.cband)
    _cdid_ord_meta_number(P,"ord.report.critical",A.critical);_cdid_ord_meta_number(P,"ord.report.normal",normal)
    _cdid_ord_meta_number(P,"ord.report.calls",A.calls);_cdid_ord_meta_number(P,"ord.report.C",C);_cdid_ord_meta_number(P,"ord.report.clustered",clustered)
    _cdid_ord_meta_matrix(P,"ord.report.point",point);_cdid_ord_meta_matrix(P,"ord.report.labels",labels)
    _cdid_ord_meta_matrix(P,"ord.report.se",se);_cdid_ord_meta_matrix(P,"ord.report.analytic_se",sqrt(diagonal(V)))
    _cdid_ord_meta_matrix(P,"ord.report.ci",ci);_cdid_ord_meta_matrix(P,"ord.report.point_ci",pointci)
    _cdid_ord_meta_matrix(P,"ord.report.V",V);_cdid_ord_meta_matrix(P,"ord.report.Omega",rows(P.cohort)*V)
    _cdid_ord_unit_matrix(P,"ord.report.IF",allIF)
    _cdid_ord_meta_matrix(P,"ord.agg.used",A.used);_cdid_ord_meta_matrix(P,"ord.agg.excluded",A.excluded)
    _cdid_ord_meta_matrix(P,"ord.agg.coefficients",A.coefficients);_cdid_ord_meta_matrix(P,"ord.agg.overall_coefficients",A.overall_coefficients)
    _cdid_ord_meta_matrix(P,"ord.agg.groups",A.groups);_cdid_ord_meta_matrix(P,"ord.agg.pg",A.pg)
    _cdid_ord_meta_number(P,"ord.agg.na_rm",R.na_rm);_cdid_ord_meta_number(P,"ord.agg.min_e",R.min_e)
    _cdid_ord_meta_number(P,"ord.agg.max_e",R.max_e);_cdid_ord_meta_number(P,"ord.agg.balance_e",R.balance_e)
    _cdid_ord_meta_text(P,"ord.agg.cluster_requested",requested_cluster);_cdid_ord_meta_text(P,"ord.agg.cluster_honored",clustered ? stored : "")
    for(j=1;j<=rows(A.call_kinds);j++) _cdid_ord_meta_text(P,"ord.agg.call."+_cdid_ord_key(j),A.call_kinds[j])
    cdid_state_validate(P,limits)
    return(P)
}
/* Read-only replay: no fitting, bootstrap generation, sorting or RNG calls. */
void cdid_ord_replay(string scalar token)
{
    struct cdid_result_state scalar P
    real colvector point,se,labels,order
    real matrix ci
    real scalar i,j,alpha,bstrap,level,simultaneous,rcs,C
    string scalar kind,type,head,cl,method
    P=cdid_state_get(token);point=_cdid_ord_matrix(P,"ord.report.point");se=_cdid_ord_matrix(P,"ord.report.se")
    labels=_cdid_ord_matrix(P,"ord.report.labels")
    kind=_cdid_ord_text(P,"ord.report.kind")
    ci=_cdid_ord_matrix(P,"ord.report.ci")
    if(rows(ci)!=rows(point) | cols(ci)!=2) ci=J(rows(point),2,.)
    alpha=_cdid_ord_number(P,"ord.report.alpha");bstrap=_cdid_ord_number(P,"ord.report.bstrap")
    level=100*(1-alpha)
    type=subinstr(kind,"ordinary_","")
    simultaneous=(st_global("e(cmd)")=="cdid" & st_global("e(band)")=="simultaneous")
    rcs=anyof(("repeated-cross-section","rcs"),_cdid_ca_opt_text(P,"ord.panel_mode","panel"))
    cl=st_global("e(clustvar)")
    C=_cdid_disp_e("N_clust")
    method=_cdid_ca_opt_text(P,"ord.method",st_global("e(method)"))
    printf("\n")
    _cdid_disp_line("Difference-in-differences (ordinary mode)",(rcs ? "Number of obs" : "Number of units"),rows(P.cohort))
    if(type=="cells") _cdid_disp_line((method!="" ? "Method: "+_cdid_disp_label(method) : ""),"Number of cells",rows(point))
    else _cdid_disp_line((method!="" ? "Method: "+_cdid_disp_label(method) : ""),
        (type=="simple" ? "" : "Number of "+(type=="dynamic" ? "event times" : (type=="group" ? "groups" : "periods"))),
        rows(point)-(type!="simple"))
    if(cl!="" & C<.) _cdid_disp_line("","Number of clusters",C)
    if(type!="cells") _cdid_disp_note("Aggregation: "+type+(type=="simple" ? " (average effect over all post-treatment cells)" : ""))
    _cdid_disp_rule(78)
    if(type=="cells" & rows(P.ledger)==rows(point) & cols(P.ledger)>=2) {
        printf("{txt}%8s %8s  %10s %10s  %23s\n","Group","Time","ATT","Std. err.",_cdid_disp_interval_head(level,simultaneous))
        _cdid_disp_rule(78)
        for(i=1;i<=rows(point);i++) printf("{txt}%8.0g %8.0g  {res}%10.7g %10.7g  %11.7g %11.7g\n",P.ledger[i,1],P.ledger[i,2],point[i],se[i],ci[i,1],ci[i,2])
    }
    else {
        head=(type=="group" ? "Group" : (type=="dynamic" ? "Event time" : (type=="calendar" ? "Time" : "")))
        printf("{txt}%-20s %10s %10s  %23s\n",head,"ATT","Std. err.",_cdid_disp_interval_head(level,simultaneous))
        _cdid_disp_rule(78)
        order=J(0,1,.)
        for(i=1;i<=rows(point);i++) if(labels[i]>=.) order=order\i
        for(i=1;i<=rows(point);i++) if(labels[i]<.) order=order\i
        for(j=1;j<=rows(order);j++) {
            i=order[j]
            if(labels[i]>=.) printf("{txt}%-20s {res}%10.7g %10.7g  %11.7g %11.7g\n","Overall",point[i],se[i],ci[i,1],ci[i,2])
            else {
                if(j==2 & labels[order[1]]>=.) _cdid_disp_rule(78)
                printf("{txt}%20s {res}%10.7g %10.7g  %11.7g %11.7g\n",strofreal(labels[i]),point[i],se[i],ci[i,1],ci[i,2])
            }
        }
    }
    _cdid_disp_rule(78)
    _cdid_disp_inference(0,(bstrap==1),simultaneous,level,(type=="cells" ? "cells" : "estimates"))
    if(type!="cells" & type!="simple") _cdid_disp_note("The overall summary uses a pointwise interval.")
    _cdid_disp_note(ustrregexra(_cdid_ca_opt_text(P,"ord.report.warning",""),"^[; ]+",""))
    if(st_global("e(cmd)")=="cdid" & st_global("e(mode)")=="ordinary")
        _cdid_disp_wald_line(st_global("e(wald_status)"),_cdid_disp_e("wald_statistic"),_cdid_disp_e("wald_df"),
            _cdid_disp_e("wald_p"),st_global("e(wald_reason)"))
}
end
