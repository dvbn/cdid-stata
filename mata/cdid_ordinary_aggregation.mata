*! Ordinary aggregation: independently written ratio-derivative operators.
*! Source definitions: did2.5.1 compute.aggte, wif, get_agg_inf_func, getSE.
*! Requires cdid_ordinary_inference. No chained weighting/inference policy used.
version 19.5
mata:
struct cdid_ord_agg_request {
    string scalar type
    real scalar min_e, max_e, balance_e, na_rm, bstrap, cband, B, alpha
}
struct cdid_ord_aggregate_result {
    string scalar type, status, reason, warning, covariance_origin
    real colvector egt, point, overall_point, intended, used, excluded, groups, pg
    real matrix IF, overall_IF, coefficients, forcing, overall_coefficients
    real colvector se, overall_se
    real matrix V, Omega, overall_V, overall_Omega, mass_IF
    real scalar critical, bstrap, cband, calls
    string colvector call_kinds
}
real colvector _cdid_ord_recode(real colvector values, real colvector axis)
{
    real scalar i
    real colvector out, pos
    out=J(rows(values),1,.)
    for(i=1;i<=rows(values);i++) {
        if(values[i]==0) { out[i]=0; continue; }
        pos=_cdid_ord_ix(axis:==values[i])
        if(rows(pos)) out[i]=pos[1]
    }
    return(out)
}
real scalar _cdid_ord_le(real scalar x, real scalar y)
{
    if(x==.b) return(y<. | y==.a | y==.b)
    if(x==.a) return(y==.a)
    if(y==.a) return(x<. | x==.a)
    if(y==.b) return(x==.b)
    return(x<=y & x<.)
}
real colvector _cdid_ord_window(real colvector x, real scalar low, real scalar high)
{
    real scalar i
    real colvector keep
    keep=J(rows(x),1,0)
    for(i=1;i<=rows(x);i++) keep[i]=_cdid_ord_le(low,x[i]) & _cdid_ord_le(x[i],high)
    return(keep)
}
/* Composition of normalized cohort masses and their derivative. The input
   source weights may contain zeros. No comparison-specific renormalization. */
real colvector _cdid_ord_ratio_IF(real matrix F, real colvector theta,
    real colvector positions, real colvector cell_mass, real colvector cell_group,
    real colvector source_group, real colvector source_weight, real scalar estimated,
    real colvector coefficient, real colvector forcing)
{
    real scalar j, total, target, k
    real colvector out, mass_term
    if(!rows(positions)) _cdid_ord_fail("No valid estimates for this aggregation.",459)
    total=sum(cell_mass[positions]); coefficient=J(rows(theta),1,0)
    if(total<=0 | total>=.) _cdid_ord_fail("Selected aggregation cohort mass is not positive and finite.",459)
    coefficient[positions]=cell_mass[positions]/total
    target=quadcross(theta[positions],coefficient[positions])
    out=F[,positions']*coefficient[positions]; forcing=J(rows(F),1,0)
    if(estimated) for(j=1;j<=rows(positions);j++) {
        k=positions[j]
        mass_term=source_weight:*(source_group:==cell_group[k]):-cell_mass[k]
        forcing=forcing:+mass_term*((theta[k]-target)/total)
    }
    return(out:+forcing)
}
struct cdid_ord_aggregate_result scalar cdid_ord_aggregate_points(real matrix full_IF,
    real colvector full_point, real colvector full_group, real colvector full_time,
    real colvector source_group, real colvector source_weight, real colvector calendar,
    struct cdid_ord_agg_request scalar R)
{
    struct cdid_ord_aggregate_result scalar O
    real scalar n,k,j,h,g,lo,hi,maxdate,minaxis,q,theta
    real colvector keep,glist,axis,group,time,point,origgroup,origtime,pg,positions
    real colvector gc,tc,cell_mass,good,balance,labels,coefficient,forcing,overall_keep
    real matrix F
    n=rows(full_IF);k=rows(full_point)
    O.type=R.type;O.status="refused";O.reason="";O.warning=""
    O.covariance_origin="IF-derived raw analytical Gram; AGGTE returns no joint covariance"
    O.egt=O.point=O.overall_point=O.se=O.overall_se=O.intended=O.used=O.excluded=O.groups=O.pg=J(0,1,.)
    O.IF=O.overall_IF=O.coefficients=O.forcing=O.overall_coefficients=O.V=O.Omega=O.overall_V=O.overall_Omega=O.mass_IF=J(0,0,.)
    O.critical=invnormal(1-R.alpha/2);O.bstrap=R.bstrap;O.cband=R.cband;O.calls=0;O.call_kinds=J(0,1,"")
    if(!any(R.type:==("simple","group","dynamic","calendar")) | (R.na_rm!=0 & R.na_rm!=1) |
        (R.bstrap!=0 & R.bstrap!=1) | (R.cband!=0 & R.cband!=1) | R.alpha<=0 | R.alpha>=1 | R.alpha>=.)
        _cdid_ord_fail("Invalid explicit ordinary aggregation request.",198)
    if(R.balance_e!=. & (R.balance_e>=. | R.balance_e<0 | R.balance_e!=floor(R.balance_e))) _cdid_ord_fail("balance_e must be absent or a nonnegative integer.",198)
    if(n<1) { O.reason="compute_inffunc=FALSE: ordinary aggregation requires stored influence functions"; return(O); }
    if(cols(full_IF)!=k | rows(full_group)!=k | rows(full_time)!=k | rows(source_group)!=n |
        rows(source_weight)!=n | missing(source_weight) | any(source_weight:<0) | sum(source_weight)<=0 |
        missing(full_group) | missing(full_time) | missing(calendar) | !rows(calendar))
        _cdid_ord_fail("Ordinary full-grid/source-unit dimensions or source weights are invalid.",498)
    O.intended=(1::k);keep=_cdid_ord_ix(full_point:<.)
    if(!R.na_rm & rows(keep)!=k) { O.reason="missing ATT values: set na.rm=TRUE to remove"; return(O); }
    if(!rows(keep)) { O.reason="all ATT values are NA"; return(O); }
    /* Source group na.rm prefilter is on ORIGINAL labels, before recoding. */
    if(R.na_rm & R.type=="group") {
        glist=uniqrows(sort(full_group[keep],1));good=J(rows(keep),1,0)
        for(j=1;j<=rows(glist);j++) {
            positions=_cdid_ord_ix((full_group[keep]:==glist[j]) :& (full_time[keep]:>=glist[j]) :&
                _cdid_ord_window(full_time[keep]:-glist[j],.b,R.max_e))
            if(rows(positions)) good=good:|(full_group[keep]:==glist[j])
        }
        keep=select(keep,good)
        if(!rows(keep)) { O.reason="no groups have nonmissing post-treatment ATT in window"; return(O); }
    }
    O.used=keep;good=J(k,1,1);good[keep]=J(rows(keep),1,0);O.excluded=_cdid_ord_ix(good)
    point=full_point[keep];origgroup=full_group[keep];origtime=full_time[keep];F=full_IF[,keep']
    glist=uniqrows(sort(origgroup,1));axis=uniqrows(sort((calendar\glist),1))
    group=_cdid_ord_recode(origgroup,axis);time=_cdid_ord_recode(origtime,axis)
    gc=_cdid_ord_recode(glist,axis);tc=uniqrows(sort(time,1));maxdate=axis[max(time)];minaxis=axis[1]
    O.groups=glist;O.pg=J(rows(glist),1,0);O.mass_IF=J(n,rows(glist),0)
    for(j=1;j<=rows(glist);j++) {
        O.pg[j]=mean(source_weight:*(source_group:==glist[j]))
        O.mass_IF[,j]=source_weight:*(source_group:==glist[j]):-O.pg[j]
    }
    cell_mass=J(rows(point),1,0)
    for(j=1;j<=rows(point);j++) cell_mass[j]=O.pg[_cdid_ord_ix(glist:==origgroup[j])[1]]
    if(R.type=="simple") {
        positions=_cdid_ord_ix((group:<=time):&_cdid_ord_window(time:-group,.b,R.max_e))
        if(!rows(positions)) { O.reason="no valid post-treatment estimates for simple aggregation"; return(O); }
        O.overall_IF=_cdid_ord_ratio_IF(F,point,positions,cell_mass,origgroup,source_group,source_weight,1,coefficient,forcing)
        O.overall_point=quadcross(point[positions],coefficient[positions])
        O.overall_coefficients=J(1,k,0);O.overall_coefficients[1,keep']=coefficient'
        O.IF=J(n,0,.);O.forcing=forcing;O.status="computed";return(O)
    }
    balance=J(rows(point),1,1)
    if(R.type=="group") labels=glist
    else if(R.type=="calendar") {
        labels=J(0,1,.)
        for(j=1;j<=rows(tc);j++) if(any((time:==tc[j]):&(group:<=time))) labels=labels\axis[tc[j]]
        if(R.min_e!=.b | R.max_e!=.a | R.balance_e<.) O.warning="mine(), maxe() and balancee() are ignored by estat calendar, as by R's aggte()"
    }
    else {
        if(R.balance_e<.) {
            if(R.balance_e<0 | R.balance_e!=floor(R.balance_e)) _cdid_ord_fail("balance_e must be absent or a nonnegative integer.",198)
            balance=(maxdate:-origgroup):>=R.balance_e
        }
        labels=uniqrows(sort(select(origtime:-origgroup,balance),1))
        if(R.balance_e<.) labels=select(labels,_cdid_ord_window(labels,R.balance_e-maxdate+minaxis,R.balance_e))
        labels=select(labels,_cdid_ord_window(labels,R.min_e,R.max_e))
    }
    if(!rows(labels)) { O.reason="no valid bins in requested aggregation window";return(O); }
    O.egt=labels;O.point=J(rows(labels),1,.);O.IF=J(n,rows(labels),.)
    O.coefficients=J(rows(labels),k,0);O.forcing=J(n,rows(labels),0)
    for(j=1;j<=rows(labels);j++) {
        if(R.type=="group") positions=_cdid_ord_ix((origgroup:==labels[j]):&(group:<=time):&_cdid_ord_window(time:-group,.b,R.max_e))
        else if(R.type=="calendar") positions=_cdid_ord_ix((origtime:==labels[j]):&(group:<=time))
        else positions=_cdid_ord_ix((origtime:-origgroup:==labels[j]):&balance)
        if(!rows(positions)) { O.reason="no valid estimates for a requested aggregation bin";return(O); }
        O.IF[,j]=_cdid_ord_ratio_IF(F,point,positions,cell_mass,origgroup,source_group,source_weight,R.type!="group",coefficient,forcing)
        O.coefficients[j,keep']=coefficient';O.forcing[,j]=forcing
        O.point[j]=quadcross(point[positions],coefficient[positions])
    }
    if(R.type=="group") {
        positions=(1::rows(labels))
        O.overall_IF=_cdid_ord_ratio_IF(O.IF,O.point,positions,O.pg,labels,source_group,source_weight,1,coefficient,forcing)
        O.overall_point=quadcross(O.point,coefficient)
    }
    else {
        overall_keep=(R.type=="dynamic" ? labels:>=0 : J(rows(labels),1,1))
        positions=_cdid_ord_ix(overall_keep)
        if(!rows(positions)) { O.reason="no nonnegative event bins for ordinary overall dynamic aggregation";return(O); }
        coefficient=J(rows(labels),1,0);coefficient[positions]=J(rows(positions),1,1/rows(positions))
        O.overall_point=quadcross(O.point[positions],coefficient[positions])
        O.overall_IF=O.IF[,positions']*coefficient[positions]
    }
    O.overall_coefficients=coefficient'*O.coefficients
    O.status="computed"
    return(O)
}
/* Inference operates on already composed full-unit IFs. Bootstrap objects are
   supplied per actual source call, never silently reused between stages. */
struct cdid_ord_aggregate_result scalar cdid_ord_aggregate_report(
    struct cdid_ord_aggregate_result scalar O, struct cdid_ord_agg_request scalar R,
    real colvector cluster, real scalar C, real scalar clustered, real scalar row_chunk,
    transmorphic scalar boot_calls)
{
    struct cdid_ord_summary scalar S
    struct cdid_ord_boot scalar M
    real scalar j, call, normal, k
    string scalar key
    if(O.status!="computed") return(O)
    call=0;normal=invnormal(1-R.alpha/2);O.critical=normal;O.call_kinds=J(0,1,"")
    k=cols(O.IF);O.se=J(k,1,.)
    if(k) {
        S=cdid_ord_if_summary(O.IF,cluster,C,clustered,R.bstrap,row_chunk)
        O.V=S.V;O.Omega=S.Omega
        for(j=1;j<=k;j++) {
            if(R.bstrap) {
                call=call+1;key=strtrim(strofreal(call,"%21.0f"))
                if(!asarray_contains(boot_calls,key)) _cdid_ord_fail("Missing separate ordinary per-bin bootstrap call.",498)
                M=asarray(boot_calls,key)
                if(M.K!=1 | M.n!=rows(O.IF) | M.C!=C | M.B!=R.B) _cdid_ord_fail("Per-bin bootstrap call shape/provenance mismatch.",498)
                O.se[j]=M.se[1];O.call_kinds=O.call_kinds\("bin:"+strtrim(strofreal(j,"%21.0f")))
            }
            else O.se[j]=S.analytic_se[j]
        }
        if(R.cband) {
            call=call+1;key=strtrim(strofreal(call,"%21.0f"))
            if(!asarray_contains(boot_calls,key)) _cdid_ord_fail("Missing distinct ordinary family-band bootstrap call.",498)
            M=asarray(boot_calls,key)
            if(M.K!=k | M.n!=rows(O.IF) | M.C!=C | M.B!=R.B) _cdid_ord_fail("Family-band bootstrap call shape/provenance mismatch.",498)
            O.critical=M.critical;O.call_kinds=O.call_kinds\"family_band"
            if(!R.bstrap) O.warning=O.warning+"; with bstrap(0) cband(1), the simultaneous band uses multiplier draws while the standard errors stay analytic, as in R's aggte()"
            if(O.critical>=. | O.critical<normal) {
                O.critical=normal;O.cband=0;O.warning=O.warning+"; the simultaneous critical value is not usable and pointwise intervals are reported, as in R's aggte()"
            }
            if(O.critical>=7) O.warning=O.warning+"; the simultaneous critical value is very large"
        }
    }
    S=cdid_ord_if_summary(O.overall_IF,cluster,C,clustered,R.bstrap,row_chunk)
    O.overall_V=S.V;O.overall_Omega=S.Omega
    if(R.bstrap) {
        call=call+1;key=strtrim(strofreal(call,"%21.0f"))
        if(!asarray_contains(boot_calls,key)) _cdid_ord_fail("Missing separate ordinary overall bootstrap call.",498)
        M=asarray(boot_calls,key)
        if(M.K!=1 | M.n!=rows(O.overall_IF) | M.C!=C | M.B!=R.B) _cdid_ord_fail("Overall bootstrap call shape/provenance mismatch.",498)
        O.overall_se=M.se[1];O.call_kinds=O.call_kinds\"overall"
    }
    else O.overall_se=S.analytic_se
    for(j=1;j<=rows(O.se);j++) if(O.se[j]<=_cdid_ord_tiny()) O.se[j]=.
    if(O.overall_se[1]<=_cdid_ord_tiny()) O.overall_se[1]=.
    O.calls=call
    return(O)
}
end
