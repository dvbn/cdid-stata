*! Read-only user result utilities for retained cdid results, current StataNow19.5.
*! Consumes cdid_result_state registry payloads in place (no copy of the
*! full-unit influence functions), the cdid_result_bridge, chained aggregation and
*! ordinary-state metadata layouts, and the current e(). No estimation, inference,
*! multiplier draws, user-data access, sorting, or state creation/mutation.
version 19.5
mata:

struct cdid_ur_view {
    pointer(struct cdid_result_state scalar) scalar p
    string scalar token, mode, kind, agg_type, ord_kind
    transmorphic scalar text_ix, num_ix, unit_ix, ledger_ix
}

struct cdid_ur_cells {
    real colvector g, time, ref, base, event, role_code, available, estimate
    real colvector se, analytic_se, ci_lower, ci_upper, lower, upper, zstat, pvalue
    string colvector role, target_id, support, reason, pt_tag
    real matrix extra
    string rowvector extra_names
    string scalar interval_kind, se_kind, status, note
    real scalar point_only, critical, normal_critical, alpha
}

struct cdid_ur_agg {
    string scalar type, interval_kind, se_kind, overall_se_kind, status
    string scalar inference_status, inference_reason, overall_status, overall_reason
    string scalar target_definition, band_family, selector_warning
    real colvector key, estimate, se, analytic_se, ci_lower, ci_upper, lower, upper
    real colvector available, restricted, n_intended, n_used, share_dependent
    real colvector share_denominator, fixed_se, overall_weight, share_cohorts
    string colvector label, reason, role, support
    real rowvector overall, ov_coef
    real scalar overall_analytic_se, overall_ci_lower, overall_ci_upper
    real scalar na_rm, min_e, max_e, balance_e, alpha, critical, normal_critical, M
    real matrix intended, used, excluded, W, coef
    real colvector ov_intended, ov_used, ov_excluded, na_excluded
}

struct cdid_ur_report {
    string colvector section, label, rname, text
    real colvector value, isnum
}

// ---------------------------------------------------------------- basics
void _cdid_ur_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}

real colvector _cdid_ur_seq(real scalar n)
{
    return(n>=1 ? (1::n) : J(0,1,.))
}

transmorphic scalar _cdid_ur_index(string colvector keys)
{
    transmorphic scalar A
    real scalar i
    A=asarray_create("string",1)
    asarray_notfound(A,0)
    for(i=1;i<=rows(keys);i++) asarray(A,keys[i],i)
    return(A)
}

// SMCL-safe display of retained free text (reasons can contain braces).
string scalar _cdid_ur_smcl(string scalar s)
{
    string scalar out
    out=subinstr(s,"{",char(1))
    out=subinstr(out,"}","{c )-}")
    return(subinstr(out,char(1),"{c -(}"))
}

string scalar _cdid_ur_num2s(real scalar x)
{
    if(x>=.) return(strofreal(x))
    return(strtrim(strofreal(x,"%21.10g")))
}

// ---------------------------------------------------------------- retained readers
string scalar _cdid_ur_text(struct cdid_ur_view scalar V, string scalar key, string scalar fallback)
{
    pointer(struct cdid_result_state scalar) scalar p
    real scalar i
    p=V.p; i=asarray(V.text_ix,key)
    if(i==0) return(fallback)
    return(p->text_settings[i,2])
}

real scalar _cdid_ur_hasnum(struct cdid_ur_view scalar V, string scalar key)
{
    return(asarray(V.num_ix,key)>0)
}

real scalar _cdid_ur_num(struct cdid_ur_view scalar V, string scalar key, real scalar fallback)
{
    pointer(struct cdid_result_state scalar) scalar p
    real scalar i
    p=V.p; i=asarray(V.num_ix,key)
    if(i==0) return(fallback)
    return(p->numeric_settings[i])
}

// Labeled matrices use the layout key.rows, key.cols, key.i.j.
real matrix _cdid_ur_matrix(struct cdid_ur_view scalar V, string scalar key)
{
    pointer(struct cdid_result_state scalar) scalar p
    real scalar r, c, i, j, k
    real matrix X
    p=V.p
    r=_cdid_ur_num(V,key+".rows",.); c=_cdid_ur_num(V,key+".cols",.)
    if(r>=. | c>=. | r<0 | c<0 | r!=floor(r) | c!=floor(c))
        _cdid_ur_fail("the retained result lacks the stored table "+key,498)
    X=J(r,c,.)
    for(j=1;j<=c;j++) for(i=1;i<=r;i++) {
        k=asarray(V.num_ix,key+"."+_cdid_rb_key(i)+"."+_cdid_rb_key(j))
        if(k==0) _cdid_ur_fail("the retained table "+key+" is incomplete",498)
        X[i,j]=p->numeric_settings[k]
    }
    return(X)
}

real matrix _cdid_ur_shape(real matrix X, real scalar r, real scalar c, string scalar key)
{
    if(r==0) return(J(0,c,.))
    if(rows(X)!=r | cols(X)!=c)
        _cdid_ur_fail("the retained table "+key+" does not have its documented dimensions",498)
    return(X)
}

real colvector _cdid_ur_cvec(struct cdid_ur_view scalar V, string scalar key, real scalar r)
{
    return(vec(_cdid_ur_shape(_cdid_ur_matrix(V,key),r,1,key)))
}

real colvector _cdid_ur_lcol(struct cdid_ur_view scalar V, string scalar name)
{
    pointer(struct cdid_result_state scalar) scalar p
    real scalar j
    p=V.p; j=asarray(V.ledger_ix,name)
    if(j==0) return(J(rows(p->point),1,.))
    return(p->ledger[,j])
}

real scalar _cdid_ur_hasu(struct cdid_ur_view scalar V, string scalar name)
{
    return(asarray(V.unit_ix,name)>0)
}

real colvector _cdid_ur_ucol(struct cdid_ur_view scalar V, string scalar name)
{
    pointer(struct cdid_result_state scalar) scalar p
    real scalar j
    p=V.p; j=asarray(V.unit_ix,name)
    if(j==0) return(J(rows(p->cohort),1,.))
    return(p->unit_metadata[,j])
}

string scalar _cdid_ur_role_name(real scalar code)
{
    return(code==3 ? "post" : (code==2 ? "pre" : (code==1 ? "reference" : "")))
}

// Chained: known never-treated has cohort . with cohort_known 1; unknown has 0.
// Ordinary: retained source cohort 0 denotes never treated.
string colvector _cdid_ur_cohort_status(struct cdid_ur_view scalar V)
{
    pointer(struct cdid_result_state scalar) scalar p
    string colvector s
    real scalar u, c
    p=V.p; s=J(rows(p->cohort),1,"")
    for(u=1;u<=rows(s);u++) {
        c=p->cohort[u]
        if(V.mode=="ordinary") s[u]=(c>=. ? "unknown" : (c==0 ? "never_treated" : "treated"))
        else s[u]=(p->cohort_known[u]!=1 ? "unknown" : (c>=. ? "never_treated" : "treated"))
    }
    return(s)
}

// ---------------------------------------------------------------- open/validate
// The current e() must be cdid results whose retained state is in this session.
// No registry is created when none exists; nothing is reconstructed from data.
struct cdid_ur_view scalar cdid_ur_open()
{
    struct cdid_ur_view scalar V
    struct cdid_state_registry scalar R
    pointer(struct cdid_result_state scalar) scalar p
    real colvector hit
    real matrix nf
    string scalar cmd, mode, token, rk, at, bk
    cmd=st_global("e(cmd)"); mode=st_global("e(mode)"); token=st_global("e(state_token)")
    if(cmd!="cdid") _cdid_ur_fail("no cdid results found; run cdid or load results with cdid_load"+(cmd=="" ? "" : " (the current results are from "+cmd+")"),301)
    if(mode!="chained" & mode!="ordinary")
        _cdid_ur_fail("the current cdid results are incomplete; rerun cdid",301)
    if(strtrim(token)=="") _cdid_ur_fail("the current cdid results are incomplete; rerun cdid",301)
    if(findexternal("CDID_PRIVATE_RESULT_REGISTRY")==NULL)
        _cdid_ur_fail("the results of cdid are no longer in memory; rerun cdid or load them with cdid_load",111)
    R=_cdid_rs_registry()
    hit=vec(selectindex(R.tokens:==token))
    if(rows(hit)!=1)
        _cdid_ur_fail("the results of cdid are no longer in memory; rerun cdid or load them with cdid_load",111)
    p=R.payloads[hit[1]]
    V.p=p; V.token=token; V.mode=mode
    V.text_ix=_cdid_ur_index(p->text_settings[,1])
    V.num_ix=_cdid_ur_index(p->numeric_setting_names)
    V.unit_ix=_cdid_ur_index(p->unit_columns')
    V.ledger_ix=_cdid_ur_index(p->ledger_columns')
    if((_cdid_ur_text(V,"ord.mode","")=="ordinary")!=(mode=="ordinary"))
        _cdid_ur_fail("e(mode) is "+mode+" but the retained state belongs to the other mode; the current results are inconsistent",498)
    if(mode=="ordinary") {
        V.ord_kind=_cdid_ur_text(V,"ord.report.kind","")
        if(substr(V.ord_kind,1,9)!="ordinary_")
            _cdid_ur_fail("the retained ordinary state records no reported result kind",498)
        V.kind=(V.ord_kind=="ordinary_cells" ? "cells" : "aggregation")
        V.agg_type=(V.kind=="aggregation" ? substr(V.ord_kind,10,.) : "")
    }
    else {
        V.ord_kind=""
        bk=_cdid_ur_text(V,"cdid_bridge.kind","")
        V.kind=(bk=="aggregation" ? "aggregation" : "cells")
        V.agg_type=(V.kind=="aggregation" ? _cdid_ur_text(V,"cdid_bridge.agg.type","") : "")
        if(V.kind=="aggregation" & V.agg_type=="")
            _cdid_ur_fail("the retained chained aggregation records no aggregation type",498)
    }
    rk=st_global("e(result_kind)")
    if(rk!="") {
        if(mode=="ordinary") {
            if(rk!=V.ord_kind) _cdid_ur_fail("e(result_kind) does not describe the retained state; the current results are inconsistent",498)
        }
        else if(rk!=V.kind)
            _cdid_ur_fail("e(result_kind) does not describe the retained state; the current results are inconsistent",498)
    }
    at=st_global("e(agg_type)")
    if(at!="" & at!=V.agg_type)
        _cdid_ur_fail("e(agg_type) does not describe the retained state; the current results are inconsistent",498)
    nf=st_numscalar("e(N_full)")
    if(rows(nf)==1 & cols(nf)==1) if(nf[1,1]<. & nf[1,1]!=rows(p->cohort))
        _cdid_ur_fail("e(N_full) does not match the retained full-unit frame; the current results are inconsistent",498)
    return(V)
}

// ---------------------------------------------------------------- cell ledger view
struct cdid_ur_cells scalar cdid_ur_cells_view(struct cdid_ur_view scalar V)
{
    struct cdid_ur_cells scalar C
    struct cdid_infer_flags scalar flags
    pointer(struct cdid_result_state scalar) scalar p
    real scalar K, j, bs, cb
    real matrix ci, pci
    real colvector keep
    string rowvector core
    string scalar key
    p=V.p; K=rows(p->point)
    C.role_code=p->role; C.available=p->available; C.estimate=p->point; C.reason=p->reason
    C.role=C.target_id=C.support=C.pt_tag=J(K,1,"")
    for(j=1;j<=K;j++) {
        key=_cdid_rb_key(j)
        C.role[j]=_cdid_ur_text(V,"cell_role:"+key,_cdid_ur_role_name(C.role_code[j]))
        C.target_id[j]=_cdid_ur_text(V,"cell_target_id:"+key,"")
        C.support[j]=_cdid_ur_text(V,"cell_support:"+key,"")
        C.pt_tag[j]=_cdid_ur_text(V,"cell_pt_window_tag:"+key,"")
    }
    C.zstat=C.pvalue=C.ci_lower=C.ci_upper=J(K,1,.)
    C.note=""; C.critical=C.normal_critical=C.alpha=.
    if(V.mode=="ordinary") {
        C.g=_cdid_ur_lcol(V,"group"); C.time=_cdid_ur_lcol(V,"time"); C.base=_cdid_ur_lcol(V,"base")
        C.ref=J(K,1,.); C.event=_cdid_ur_lcol(V,"event")
        core=("group","time","base","role","available","ATT","event")
        C.point_only=(_cdid_ur_num(V,"ord.cells.compute_inffunc",.)==0)
        C.alpha=_cdid_ur_num(V,"ord.cells.alpha",.)
        if(V.ord_kind=="ordinary_cells") {
            C.se=vec(_cdid_ur_matrix(V,"ord.report.se"))
            C.analytic_se=vec(_cdid_ur_matrix(V,"ord.report.analytic_se"))
            ci=_cdid_ur_matrix(V,"ord.report.ci"); pci=_cdid_ur_matrix(V,"ord.report.point_ci")
            C.critical=_cdid_ur_num(V,"ord.report.critical",.)
            C.normal_critical=_cdid_ur_num(V,"ord.report.normal",.)
        }
        else {
            // An ordinary aggregation replaces ord.report.*; ord.cells.* keeps the
            // retained cell SE and active interval. e(analytic_se) is defined as
            // sqrt(diagonal(e(V))) of the retained estimator covariance.
            C.se=vec(_cdid_ur_matrix(V,"ord.cells.se")); ci=_cdid_ur_matrix(V,"ord.cells.ci"); pci=J(K,2,.)
            if(K>0 & rows(p->V_estimator)==K) C.analytic_se=sqrt(diagonal(p->V_estimator))
            else C.analytic_se=J(K,1,.)
            C.note="cell pointwise intervals are not retained after an ordinary aggregation; lower/upper are the retained cell intervals"
        }
        if(K>0) {
            if(rows(C.se)!=K | rows(C.analytic_se)!=K | rows(ci)!=K | cols(ci)!=2 | rows(pci)!=K | cols(pci)!=2)
                _cdid_ur_fail("retained ordinary cell inference does not align with the table of cells",498)
            C.lower=ci[,1]; C.upper=ci[,2]; C.ci_lower=pci[,1]; C.ci_upper=pci[,2]
        }
        else C.se=C.analytic_se=C.lower=C.upper=J(0,1,.)
        bs=_cdid_ur_num(V,"ord.cells.bstrap",0); cb=_cdid_ur_num(V,"ord.cells.cband",0)
        C.status=(C.point_only ? "point_estimates_only" : "computed")
        C.interval_kind=(C.point_only ? "none" : ((bs==1 & cb==1) ? "simultaneous_band" : "pointwise"))
        C.se_kind=(C.point_only ? "none" : (bs==1 ? "multiplier_bootstrap_iqr" : "source_analytic"))
    }
    else {
        C.g=_cdid_ur_lcol(V,"g"); C.time=_cdid_ur_lcol(V,"time")
        C.ref=_cdid_ur_lcol(V,"ref"); C.base=_cdid_ur_lcol(V,"base_t")
        C.event=C.time:-C.g
        core=("g","time","ref","base_t"),"cdid_bridge.cells.":+("se","analytic_se","ci_lower","ci_upper","band_lower","band_upper","z","p")
        C.se=_cdid_ur_lcol(V,"cdid_bridge.cells.se"); C.analytic_se=_cdid_ur_lcol(V,"cdid_bridge.cells.analytic_se")
        C.ci_lower=_cdid_ur_lcol(V,"cdid_bridge.cells.ci_lower"); C.ci_upper=_cdid_ur_lcol(V,"cdid_bridge.cells.ci_upper")
        C.lower=_cdid_ur_lcol(V,"cdid_bridge.cells.band_lower"); C.upper=_cdid_ur_lcol(V,"cdid_bridge.cells.band_upper")
        C.zstat=_cdid_ur_lcol(V,"cdid_bridge.cells.z"); C.pvalue=_cdid_ur_lcol(V,"cdid_bridge.cells.p")
        if(asarray(V.ledger_ix,"cdid_bridge.cells.se")==0) {
            C.status="inference_not_attached"; C.point_only=.
            C.interval_kind=C.se_kind="none"
            C.note="the retained state carries no cell inference fields"
        }
        else {
            C.status=_cdid_ur_text(V,"cdid_bridge.cells.status","")
            C.point_only=(C.status=="inference_not_computed")
            C.alpha=_cdid_ur_num(V,"cdid_bridge.cells.alpha",.)
            C.critical=_cdid_ur_num(V,"cdid_bridge.cells.critical",.)
            C.normal_critical=_cdid_ur_num(V,"cdid_bridge.cells.normal_critical",.)
            if(C.point_only) C.interval_kind=C.se_kind="none"
            else {
                bs=_cdid_ur_num(V,"cdid_bridge.cells.bstrap",0); cb=_cdid_ur_num(V,"cdid_bridge.cells.cband",0)
                flags=cdid_infer_resolve_flags(bs,cb,"cells")
                C.interval_kind=(flags.simultaneous ? "simultaneous_band" : "pointwise")
                C.se_kind=(flags.bootstrap ? "multiplier_bootstrap_iqr" : "analytic")
            }
        }
    }
    keep=J(0,1,.)
    for(j=1;j<=cols(p->ledger_columns);j++) if(!any(core:==p->ledger_columns[j])) keep=keep\j
    if(rows(keep)) {
        C.extra=p->ledger[,keep']; C.extra_names=p->ledger_columns[keep']
    }
    else {
        C.extra=J(K,0,.); C.extra_names=J(1,0,"")
    }
    return(C)
}

// ---------------------------------------------------------------- aggregation views
string scalar _cdid_ur_fmtkey(real scalar x)
{
    return(strtrim(strofreal(x,"%21.0g")))
}

struct cdid_ur_agg scalar _cdid_ur_agg_ord(struct cdid_ur_view scalar V)
{
    struct cdid_ur_agg scalar A
    pointer(struct cdid_result_state scalar) scalar p
    real colvector pt, lab, se, ase
    real matrix ci, pci, cf
    real scalar n, M, K, m, bs, cb
    string scalar prefix
    p=V.p; K=rows(p->point); A.type=V.agg_type
    pt=vec(_cdid_ur_matrix(V,"ord.report.point")); n=rows(pt); M=n-1
    if(M<0) _cdid_ur_fail("the retained ordinary aggregation has no overall row",498)
    lab=vec(_cdid_ur_matrix(V,"ord.report.labels")); se=vec(_cdid_ur_matrix(V,"ord.report.se"))
    ase=vec(_cdid_ur_matrix(V,"ord.report.analytic_se"))
    ci=_cdid_ur_matrix(V,"ord.report.ci"); pci=_cdid_ur_matrix(V,"ord.report.point_ci")
    if(rows(lab)!=n | rows(se)!=n | rows(ase)!=n | rows(ci)!=n | cols(ci)!=2 | rows(pci)!=n | cols(pci)!=2)
        _cdid_ur_fail("the retained ordinary aggregation tables do not align",498)
    A.M=M
    if(M) {
        A.key=lab[|1\M|]; A.estimate=pt[|1\M|]; A.se=se[|1\M|]; A.analytic_se=ase[|1\M|]
        A.lower=ci[|1,1\M,1|]; A.upper=ci[|1,2\M,2|]; A.ci_lower=pci[|1,1\M,1|]; A.ci_upper=pci[|1,2\M,2|]
    }
    else A.key=A.estimate=A.se=A.analytic_se=A.lower=A.upper=A.ci_lower=A.ci_upper=J(0,1,.)
    A.available=(A.estimate:<.)
    A.restricted=A.n_intended=A.share_dependent=A.share_denominator=A.fixed_se=A.overall_weight=J(M,1,.)
    prefix=(A.type=="group" ? "group_" : (A.type=="calendar" ? "calendar_" : "event_"))
    A.label=A.reason=A.support=J(M,1,""); A.role=J(M,1,"aggregate")
    for(m=1;m<=M;m++) A.label[m]=prefix+_cdid_ur_fmtkey(A.key[m])
    cf=_cdid_ur_matrix(V,"ord.agg.coefficients")
    A.coef=_cdid_ur_shape(cf,M,K,"ord.agg.coefficients")
    A.n_used=J(M,1,.)
    for(m=1;m<=M;m++) A.n_used[m]=sum((A.coef[m,]:!=0):&(A.coef[m,]:<.))
    A.ov_coef=_cdid_ur_shape(_cdid_ur_matrix(V,"ord.agg.overall_coefficients"),1,K,"ord.agg.overall_coefficients")
    A.na_excluded=vec(_cdid_ur_matrix(V,"ord.agg.excluded"))
    A.share_cohorts=vec(_cdid_ur_matrix(V,"ord.agg.groups"))
    A.intended=A.used=A.excluded=A.W=J(0,0,.)
    A.ov_intended=A.ov_used=A.ov_excluded=J(0,1,.)
    A.overall=(pt[n],se[n],ci[n,1],ci[n,2]); A.overall_ci_lower=pci[n,1]; A.overall_ci_upper=pci[n,2]
    A.overall_analytic_se=ase[n]
    A.na_rm=_cdid_ur_num(V,"ord.agg.na_rm",.); A.min_e=_cdid_ur_num(V,"ord.agg.min_e",.)
    A.max_e=_cdid_ur_num(V,"ord.agg.max_e",.); A.balance_e=_cdid_ur_num(V,"ord.agg.balance_e",.)
    A.alpha=_cdid_ur_num(V,"ord.report.alpha",.); A.critical=_cdid_ur_num(V,"ord.report.critical",.)
    A.normal_critical=_cdid_ur_num(V,"ord.report.normal",.)
    bs=_cdid_ur_num(V,"ord.report.bstrap",0); cb=_cdid_ur_num(V,"ord.report.cband",0)
    A.interval_kind=(cb==1 ? "simultaneous_band" : "pointwise")
    A.se_kind=(bs==1 ? "multiplier_bootstrap_iqr" : "analytic"); A.overall_se_kind=A.se_kind
    A.status="computed"; A.inference_status="computed"
    A.inference_reason=_cdid_ur_text(V,"ord.report.covariance_origin","")
    A.overall_status=(pt[n]<. ? "available" : "unavailable"); A.overall_reason=""
    A.target_definition=""; A.band_family=""
    A.selector_warning=_cdid_ur_text(V,"ord.report.warning","")
    return(A)
}

struct cdid_ur_agg scalar cdid_ur_agg_view(struct cdid_ur_view scalar V)
{
    struct cdid_ur_agg scalar A
    struct cdid_infer_flags scalar flags
    pointer(struct cdid_result_state scalar) scalar p
    real matrix pci, bci, ov, oase
    real scalar M, R, j, bs, cb
    string scalar pre, key
    if(V.kind!="aggregation")
        _cdid_ur_fail("the current results are not an aggregation; run estat simple, group, dynamic or calendar first",498)
    if(V.mode=="ordinary") return(_cdid_ur_agg_ord(V))
    p=V.p; R=rows(p->point); pre="cdid_bridge.agg."
    A.type=V.agg_type
    A.estimate=vec(_cdid_ur_matrix(V,pre+"points")); M=rows(A.estimate); A.M=M
    A.key=_cdid_ur_cvec(V,pre+"keys",M)
    A.available=_cdid_ur_cvec(V,pre+"available",M)
    A.restricted=_cdid_ur_cvec(V,pre+"restricted",M)
    A.se=_cdid_ur_cvec(V,pre+"se",M)
    A.analytic_se=_cdid_ur_cvec(V,pre+"analytic_se",M)
    A.fixed_se=_cdid_ur_cvec(V,pre+"fixed_se",M)
    A.n_intended=_cdid_ur_cvec(V,pre+"n_intended",M)
    A.n_used=_cdid_ur_cvec(V,pre+"n_used",M)
    A.share_dependent=_cdid_ur_cvec(V,pre+"share_dependent",M)
    A.share_denominator=_cdid_ur_cvec(V,pre+"share_denominator",M)
    A.overall_weight=_cdid_ur_cvec(V,pre+"overall_weights",M)
    pci=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"point_ci"),M,2,pre+"point_ci")
    bci=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"band_ci"),M,2,pre+"band_ci")
    if(M) {
        A.ci_lower=pci[,1]; A.ci_upper=pci[,2]; A.lower=bci[,1]; A.upper=bci[,2]
    }
    else A.ci_lower=A.ci_upper=A.lower=A.upper=J(0,1,.)
    A.intended=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"intended"),M,R,pre+"intended")
    A.used=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"used"),M,R,pre+"used")
    A.excluded=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"excluded"),M,R,pre+"excluded")
    A.W=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"W"),M,R,pre+"W")
    ov=_cdid_ur_shape(_cdid_ur_matrix(V,pre+"overall"),1,4,pre+"overall")
    A.overall=ov; A.overall_ci_lower=ov[1,3]; A.overall_ci_upper=ov[1,4]
    oase=_cdid_ur_matrix(V,pre+"overall_analytic_se")
    A.overall_analytic_se=(length(oase) ? oase[1,1] : .)
    A.ov_intended=vec(_cdid_ur_matrix(V,pre+"overall_intended"))
    A.ov_used=vec(_cdid_ur_matrix(V,pre+"overall_used"))
    A.ov_excluded=vec(_cdid_ur_matrix(V,pre+"overall_excluded"))
    A.share_cohorts=vec(_cdid_ur_matrix(V,pre+"share_cohorts"))
    A.coef=J(0,0,.); A.ov_coef=J(1,0,.); A.na_excluded=J(0,1,.)
    A.label=A.reason=A.role=A.support=J(M,1,"")
    for(j=1;j<=M;j++) {
        key=_cdid_rb_key(j)
        A.label[j]=_cdid_ur_text(V,pre+"label."+key,"")
        A.reason[j]=_cdid_ur_text(V,pre+"reason."+key,"")
        A.role[j]=_cdid_ur_text(V,pre+"role."+key,"")
        A.support[j]=_cdid_ur_text(V,pre+"support."+key,"")
    }
    A.na_rm=_cdid_ur_num(V,pre+"na_rm",.); A.min_e=_cdid_ur_num(V,pre+"min_e",.)
    A.max_e=_cdid_ur_num(V,pre+"max_e",.); A.balance_e=_cdid_ur_num(V,pre+"balance_e",.)
    A.alpha=_cdid_ur_num(V,"cdid_bridge.alpha",.); A.critical=_cdid_ur_num(V,"cdid_bridge.critical",.)
    A.normal_critical=_cdid_ur_num(V,"cdid_bridge.normal_critical",.)
    bs=_cdid_ur_num(V,"cdid_bridge.bstrap",0); cb=_cdid_ur_num(V,"cdid_bridge.cband",0)
    flags=cdid_infer_resolve_flags(bs,cb,(A.type=="simple" ? "simple" : "aggte_bins"))
    A.interval_kind=(flags.simultaneous ? "simultaneous_band" : "pointwise")
    A.se_kind=(flags.bootstrap ? "multiplier_bootstrap_iqr" : "analytic")
    A.overall_se_kind=_cdid_ur_text(V,pre+"overall.inference_kind","")
    A.status=_cdid_ur_text(V,"cdid_bridge.status","")
    A.inference_status=_cdid_ur_text(V,"cdid_bridge.inference_status","")
    A.inference_reason=_cdid_ur_text(V,pre+"inference_reason","")
    A.overall_status=_cdid_ur_text(V,pre+"overall.status","")
    A.overall_reason=_cdid_ur_text(V,pre+"overall.reason","")
    A.target_definition=_cdid_ur_text(V,pre+"target_definition","")
    A.band_family=_cdid_ur_text(V,pre+"band_family","")
    A.selector_warning=_cdid_ur_text(V,pre+"selector_warning","")
    return(A)
}

// ---------------------------------------------------------------- dataset writers
void _cdid_ur_newobs(real scalar n)
{
    if(st_nobs()!=0 | st_nvar()!=0) _cdid_ur_fail("internal: the output frame is not empty",498)
    if(n>0) st_addobs(n)
}

void _cdid_ur_nv(string scalar name, string scalar type, real colvector x, string scalar label)
{
    (void) st_addvar(type,name)
    if(rows(x)) st_store(.,name,x)
    st_varlabel(name,substr(label,1,80))
}

void _cdid_ur_sv(string scalar name, string colvector s, string scalar label)
{
    real scalar L
    L=(rows(s) ? max(strlen(s)) : 0)
    (void) st_addvar((L>2045 ? "strL" : "str"+strofreal(max((1,L)))),name)
    if(rows(s)) st_sstore(.,name,s)
    st_varlabel(name,substr(label,1,80))
}

// Retained column names become valid, unique Stata names; the exact source
// name is kept in the variable characteristic [cdid_source] and its label.
string scalar _cdid_ur_varname(string scalar raw, string scalar prefix, transmorphic scalar used)
{
    string scalar s, base, name, suf
    real rowvector a, ok
    real scalar k
    s=raw
    if(substr(s,1,12)=="cdid_bridge.") s=substr(s,13,.)
    if(s=="") s="column"
    a=ascii(s)
    ok=((a:>=48):&(a:<=57)):|((a:>=65):&(a:<=90)):|((a:>=97):&(a:<=122)):|(a:==95)
    a=a:*ok+95:*(1:-ok)
    base=substr(prefix+char(a),1,32); name=base; k=1
    while(asarray_contains(used,name)) {
        k=k+1; suf="_"+_cdid_rb_key(k)
        name=substr(base,1,32-strlen(suf))+suf
    }
    asarray(used,name,1)
    return(name)
}

void _cdid_ur_srcvar(string scalar raw, real colvector x, transmorphic scalar used)
{
    string scalar vn
    vn=_cdid_ur_varname(raw,"src_",used)
    _cdid_ur_nv(vn,"double",x,raw)
    st_global(vn+"[cdid_source]",raw)
}

void _cdid_ur_idvar(struct cdid_ur_view scalar V)
{
    pointer(struct cdid_result_state scalar) scalar p
    p=V.p
    if(p->id_kind=="string") _cdid_ur_sv("id",p->unit_string_id,"Original unit identifier (string)")
    else _cdid_ur_nv("id","double",p->unit_numeric_id,"Original unit identifier (numeric)")
}

void _cdid_ur_chars(struct cdid_ur_view scalar V, string scalar type)
{
    pointer(struct cdid_result_state scalar) scalar p
    string scalar sv
    real scalar i
    p=V.p; sv=""
    for(i=1;i<=rows(p->source_versions);i++)
        sv=sv+(i>1 ? "; " : "")+p->source_versions[i,1]+"="+p->source_versions[i,2]
    st_global("_dta[cdid_export_type]",type)
    st_global("_dta[cdid_mode]",V.mode)
    st_global("_dta[cdid_result_kind]",V.kind)
    st_global("_dta[cdid_agg_type]",V.agg_type)
    st_global("_dta[cdid_N_full]",_cdid_rb_key(rows(p->cohort)))
    st_global("_dta[cdid_K_cells]",_cdid_rb_key(rows(p->point)))
    st_global("_dta[cdid_source_versions]",sv)
    st_global("_dta[cdid_cmdline]",st_global("e(cmdline)"))
}

void cdid_ur_export_cells(struct cdid_ur_view scalar V)
{
    struct cdid_ur_cells scalar C
    transmorphic scalar used
    real scalar K, j
    C=cdid_ur_cells_view(V); K=rows(C.estimate)
    _cdid_ur_newobs(K)
    used=_cdid_ur_index(("cell","g","time","ref","base","event","role_code","role","available",
        "estimate","se","analytic_se","lower","upper","ci_lower","ci_upper","zstat","pvalue",
        "interval_kind","se_kind","target_id","support","reason","pt_window_tag")')
    _cdid_ur_nv("cell","long",_cdid_ur_seq(K),"Position in the retained original cell ledger")
    _cdid_ur_nv("g","double",C.g,"Cohort (first-treatment date)")
    _cdid_ur_nv("time","double",C.time,"Calendar date")
    _cdid_ur_nv("ref","double",C.ref,"Reference date (chained ledger)")
    _cdid_ur_nv("base","double",C.base,"Base date")
    _cdid_ur_nv("event","double",C.event,"Event time t - g in date units")
    _cdid_ur_nv("role_code","byte",C.role_code,"Role code: 1 reference, 2 pre-trend, 3 post target")
    _cdid_ur_sv("role",C.role,"Retained role label")
    _cdid_ur_nv("available","byte",C.available,"Retained availability (1 = available)")
    _cdid_ur_nv("estimate","double",C.estimate,"Retained point estimate")
    _cdid_ur_nv("se","double",C.se,"Reported standard error (retained)")
    _cdid_ur_nv("analytic_se","double",C.analytic_se,"Analytical SE from retained estimator covariance")
    _cdid_ur_nv("lower","double",C.lower,"Reported interval lower bound (retained)")
    _cdid_ur_nv("upper","double",C.upper,"Reported interval upper bound (retained)")
    _cdid_ur_nv("ci_lower","double",C.ci_lower,"Pointwise interval lower bound (retained)")
    _cdid_ur_nv("ci_upper","double",C.ci_upper,"Pointwise interval upper bound (retained)")
    _cdid_ur_nv("zstat","double",C.zstat,"z statistic (retained)")
    _cdid_ur_nv("pvalue","double",C.pvalue,"p-value (retained)")
    _cdid_ur_sv("interval_kind",J(K,1,C.interval_kind),"Reported interval: simultaneous_band, pointwise or none")
    _cdid_ur_sv("se_kind",J(K,1,C.se_kind),"Reported SE kind")
    _cdid_ur_sv("target_id",C.target_id,"Retained target identifier")
    _cdid_ur_sv("support",C.support,"Retained support status")
    _cdid_ur_sv("reason",cdid_text_unavailable_v(C.reason),"Availability reason")
    _cdid_ur_sv("pt_window_tag",C.pt_tag,"Retained parallel-trends window tag")
    for(j=1;j<=cols(C.extra_names);j++) _cdid_ur_srcvar(C.extra_names[j],C.extra[,j],used)
    _cdid_ur_chars(V,"cells")
    st_global("_dta[cdid_inference_status]",C.status)
    st_global("_dta[cdid_alpha]",_cdid_ur_num2s(C.alpha))
    st_global("_dta[cdid_critical]",_cdid_ur_num2s(C.critical))
    st_global("_dta[cdid_normal_critical]",_cdid_ur_num2s(C.normal_critical))
    st_global("_dta[cdid_note]",C.note)
}

void cdid_ur_export_units(struct cdid_ur_view scalar V)
{
    pointer(struct cdid_result_state scalar) scalar p
    transmorphic scalar used
    string scalar nm
    real scalar N, j
    p=V.p; N=rows(p->cohort)
    _cdid_ur_newobs(N)
    used=_cdid_ur_index(("unit","id","cohort","cohort_known","cohort_status","weight","cluster")')
    _cdid_ur_nv("unit","long",_cdid_ur_seq(N),"Original unit position in the retained full-unit frame")
    _cdid_ur_idvar(V)
    _cdid_ur_nv("cohort","double",p->cohort,"Retained cohort")
    _cdid_ur_nv("cohort_known","byte",p->cohort_known,"Retained cohort-known flag")
    _cdid_ur_sv("cohort_status",_cdid_ur_cohort_status(V),"treated, never_treated or unknown")
    _cdid_ur_nv("weight","double",p->weight,"Retained unit weight")
    if(p->cluster_kind=="numeric") _cdid_ur_nv("cluster","double",p->cluster_numeric_id,"Retained cluster identifier")
    else if(p->cluster_kind=="string") _cdid_ur_sv("cluster",p->cluster_string_id,"Retained cluster identifier")
    for(j=1;j<=cols(p->unit_columns);j++) {
        nm=p->unit_columns[j]
        if(substr(nm,1,16)=="cdid_bridge.agg." | substr(nm,1,14)=="ord.report.IF.") continue
        _cdid_ur_srcvar(nm,p->unit_metadata[,j],used)
    }
    _cdid_ur_chars(V,"units")
    st_global("_dta[cdid_id_kind]",p->id_kind)
    st_global("_dta[cdid_cluster_kind]",p->cluster_kind)
}

// Long layout: one row per original unit and retained influence target.
// Families: 1 cell, 2 fit_share, 3 aggregate, 4 aggregate_fixed,
// 5 aggregate_share, 6 overall, 7 aggregate_cohort_share.
void cdid_ur_export_influence(struct cdid_ur_view scalar V)
{
    pointer(struct cdid_result_state scalar) scalar p
    struct cdid_ur_cells scalar C
    struct cdid_ur_agg scalar A
    real matrix plan, X
    real colvector val, obs, units
    real rowvector ix
    real scalar N, K, t, row, total, cap, cc, j, col, idv
    string scalar fams, unitkey, ovkey
    string colvector famnames
    p=V.p; N=rows(p->cohort); K=rows(p->point)
    if(rows(p->IF)==0)
        _cdid_ur_fail("the current results have no influence functions (noinffunc); rerun cdid without noinffunc",459)
    C=cdid_ur_cells_view(V)
    // plan columns: family, target, source(0 IF,1 share_IF,2 unit_metadata), column, g, time, key, available
    plan=J(K,8,.)
    for(j=1;j<=K;j++) plan[j,]=(1,j,0,j,C.g[j],C.time[j],.,C.available[j])
    if(rows(p->share_IF)==N & cols(p->share_IF)>0) {
        cc=asarray(_cdid_ur_index(p->share_columns'),"cohort")
        for(j=1;j<=cols(p->share_IF);j++)
            plan=plan\(2,j,1,j,((cc>0 & rows(p->shares)>=j) ? p->shares[j,cc] : .),.,.,1)
    }
    if(V.kind=="aggregation") {
        A=cdid_ur_agg_view(V)
        if(V.mode=="ordinary") {
            unitkey="ord.report.IF."
            for(j=1;j<=A.M;j++) {
                col=asarray(V.unit_ix,unitkey+_cdid_rb_key(j))
                if(col==0) _cdid_ur_fail("the retained ordinary aggregate influence functions are incomplete",498)
                plan=plan\(3,j,2,col,.,.,A.key[j],A.available[j])
            }
            col=asarray(V.unit_ix,unitkey+_cdid_rb_key(A.M+1))
            if(col==0) _cdid_ur_fail("the retained ordinary overall influence function is absent",498)
            plan=plan\(6,1,2,col,.,.,.,(A.overall_status=="available"))
        }
        else {
            famnames=("IF_total"\"IF_fixed"\"IF_share")
            for(t=1;t<=3;t++) for(j=1;j<=A.M;j++) {
                col=asarray(V.unit_ix,"cdid_bridge.agg."+famnames[t]+"."+_cdid_rb_key(j))
                if(col) plan=plan\(2+t,j,2,col,.,.,A.key[j],A.available[j])
            }
            ovkey="cdid_bridge.agg.overall_IF.1"
            col=asarray(V.unit_ix,ovkey)
            if(col) plan=plan\(6,1,2,col,.,.,.,(A.overall_status=="available"))
            for(j=1;j<=rows(A.share_cohorts);j++) {
                col=asarray(V.unit_ix,"cdid_bridge.agg.cohort_share_IF."+_cdid_rb_key(j))
                if(col) plan=plan\(7,j,2,col,A.share_cohorts[j],.,.,1)
            }
        }
    }
    total=N*rows(plan)
    cap=c("max_N_theory")
    if(cap<. & total>cap)
        _cdid_ur_fail("the influence-function dataset needs "+_cdid_rb_key(total)+" observations, more than this Stata allows",908)
    _cdid_ur_newobs(total)
    ix=st_addvar(("long","byte","long","double","double","double","byte","double"),
        ("unit","family","target","g","time","key","available","value"))
    idv=(p->id_kind=="string")
    if(idv) (void) st_addvar((max(strlen(p->unit_string_id))>2045 ? "strL" : "str"+strofreal(max((1,max(strlen(p->unit_string_id)))))),"id")
    else (void) st_addvar("double","id")
    units=_cdid_ur_seq(N); row=0
    for(t=1;t<=rows(plan);t++) {
        obs=_cdid_ur_seq(N):+row
        if(plan[t,3]==0) val=p->IF[,plan[t,4]]
        else if(plan[t,3]==1) val=p->share_IF[,plan[t,4]]
        else val=p->unit_metadata[,plan[t,4]]
        X=(units,J(N,1,plan[t,1]),J(N,1,plan[t,2]),J(N,1,plan[t,5]),J(N,1,plan[t,6]),J(N,1,plan[t,7]),J(N,1,plan[t,8]),val)
        st_store(obs,ix,X)
        if(idv) st_sstore(obs,"id",p->unit_string_id)
        else st_store(obs,"id",p->unit_numeric_id)
        row=row+N
    }
    st_varlabel("unit","Original unit position in the retained full-unit frame")
    st_varlabel("id","Original unit identifier")
    st_varlabel("family","Influence-function family")
    st_varlabel("target","Target position within its family")
    st_varlabel("g","Cohort (cell and cohort-share families)")
    st_varlabel("time","Calendar date (cell family)")
    st_varlabel("key","Aggregation bin key (aggregate families)")
    st_varlabel("available","Retained target availability")
    st_varlabel("value","Retained influence-function value")
    famnames=("cell"\"fit_share"\"aggregate"\"aggregate_fixed"\"aggregate_share"\"overall"\"aggregate_cohort_share")
    st_vlmodify("cdid_if_family",(1::7),famnames)
    st_varvaluelabel("family","cdid_if_family")
    fams=""
    for(j=1;j<=7;j++) if(any(plan[,1]:==j)) fams=fams+(fams=="" ? "" : " ")+famnames[j]
    _cdid_ur_chars(V,"influence")
    st_global("_dta[cdid_if_families]",fams)
    st_global("_dta[cdid_if_scale]",_cdid_ur_text(V,"IF_scale","not recorded in the retained state"))
    st_global("_dta[cdid_covariance_scale]",(st_global("e(covariance_scale)")=="" ? "not recorded in the current e()" : st_global("e(covariance_scale)")))
}

void cdid_ur_export_agg(struct cdid_ur_view scalar V)
{
    struct cdid_ur_agg scalar A
    struct cdid_ur_cells scalar C
    real colvector mb, mc, ob, oc, idx, bin, key, est, se, ase, lo, hi, cil, ciu, fse
    real colvector av, rs, nint, nuse, sdep, sden, ovw, cell, g, tt, intd, usd, exd, wt
    string colvector rec, lab, st, role, sup, why, ik, sk
    real scalar M, R, n, m, r, i, ord, nx, nm, no
    A=cdid_ur_agg_view(V); C=cdid_ur_cells_view(V)
    M=A.M; R=rows(C.estimate); ord=(V.mode=="ordinary")
    mb=mc=ob=oc=J(0,1,.)
    for(m=1;m<=M;m++) {
        if(ord) idx=vec(selectindex((A.coef[m,]:!=0):&(A.coef[m,]:<.)))
        else idx=vec(selectindex((A.intended[m,]:==1):|(A.used[m,]:==1):|(A.excluded[m,]:==1):|((A.W[m,]:!=0):&(A.W[m,]:<.))))
        if(rows(idx)) {
            mb=mb\J(rows(idx),1,m); mc=mc\idx
        }
    }
    if(ord) {
        oc=vec(selectindex((A.ov_coef:!=0):&(A.ov_coef:<.))); ob=J(rows(oc),1,.)
    }
    else {
        ob=A.ov_intended; oc=J(rows(ob),1,.)
    }
    nx=rows(A.na_excluded); nm=rows(mb); no=rows(ob)
    n=M+1+nm+no+nx
    rec=lab=st=role=sup=why=ik=sk=J(n,1,"")
    bin=key=est=se=ase=lo=hi=cil=ciu=fse=av=rs=nint=nuse=sdep=sden=ovw=cell=g=tt=intd=usd=exd=wt=J(n,1,.)
    for(m=1;m<=M;m++) {
        rec[m]="bin"; bin[m]=m; key[m]=A.key[m]; lab[m]=A.label[m]
        est[m]=A.estimate[m]; se[m]=A.se[m]; ase[m]=A.analytic_se[m]; lo[m]=A.lower[m]; hi[m]=A.upper[m]
        cil[m]=A.ci_lower[m]; ciu[m]=A.ci_upper[m]; fse[m]=A.fixed_se[m]; av[m]=A.available[m]
        rs[m]=A.restricted[m]; nint[m]=A.n_intended[m]; nuse[m]=A.n_used[m]; sdep[m]=A.share_dependent[m]
        sden[m]=A.share_denominator[m]; ovw[m]=A.overall_weight[m]
        st[m]=(A.available[m]==1 ? "available" : "unavailable"); role[m]=A.role[m]; sup[m]=A.support[m]; why[m]=A.reason[m]
        ik[m]=A.interval_kind; sk[m]=A.se_kind
    }
    i=M+1
    rec[i]="overall"; lab[i]="overall"; est[i]=A.overall[1]; se[i]=A.overall[2]; lo[i]=A.overall[3]; hi[i]=A.overall[4]
    cil[i]=A.overall_ci_lower; ciu[i]=A.overall_ci_upper; ase[i]=A.overall_analytic_se
    av[i]=(A.overall_status=="available"); st[i]=A.overall_status; why[i]=A.overall_reason
    ik[i]="pointwise"; sk[i]=A.overall_se_kind
    for(r=1;r<=nm;r++) {
        i=M+1+r; m=mb[r]; rec[i]="member"; bin[i]=m; key[i]=A.key[m]; lab[i]=A.label[m]
        cell[i]=mc[r]; g[i]=C.g[mc[r]]; tt[i]=C.time[mc[r]]
        if(ord) {
            usd[i]=1; wt[i]=A.coef[m,mc[r]]
        }
        else {
            intd[i]=A.intended[m,mc[r]]; usd[i]=A.used[m,mc[r]]; exd[i]=A.excluded[m,mc[r]]; wt[i]=A.W[m,mc[r]]
        }
    }
    for(r=1;r<=no;r++) {
        i=M+1+nm+r; rec[i]="overall_member"
        if(ord) {
            cell[i]=oc[r]; g[i]=C.g[oc[r]]; tt[i]=C.time[oc[r]]; usd[i]=1; wt[i]=A.ov_coef[oc[r]]
        }
        else {
            m=ob[r]; bin[i]=m; key[i]=A.key[m]; lab[i]=A.label[m]; intd[i]=1
            usd[i]=any(A.ov_used:==m); exd[i]=any(A.ov_excluded:==m); wt[i]=A.overall_weight[m]
        }
    }
    for(r=1;r<=nx;r++) {
        i=M+1+nm+no+r; rec[i]="na_rm_excluded"; cell[i]=A.na_excluded[r]
        g[i]=C.g[A.na_excluded[r]]; tt[i]=C.time[A.na_excluded[r]]; exd[i]=1
    }
    _cdid_ur_newobs(n)
    _cdid_ur_sv("record",rec,"bin, overall, member, overall_member or na_rm_excluded")
    _cdid_ur_sv("agg_type",J(n,1,A.type),"Aggregation type")
    _cdid_ur_nv("bin","long",bin,"Aggregation bin position")
    _cdid_ur_nv("key","double",key,"Bin key: cohort, calendar date or event time")
    _cdid_ur_sv("label",lab,"Retained bin label")
    _cdid_ur_nv("estimate","double",est,"Retained aggregate estimate")
    _cdid_ur_nv("se","double",se,"Reported standard error (retained)")
    _cdid_ur_nv("analytic_se","double",ase,"Analytical SE (retained)")
    _cdid_ur_nv("lower","double",lo,"Reported interval lower bound (retained)")
    _cdid_ur_nv("upper","double",hi,"Reported interval upper bound (retained)")
    _cdid_ur_nv("ci_lower","double",cil,"Pointwise interval lower bound (retained)")
    _cdid_ur_nv("ci_upper","double",ciu,"Pointwise interval upper bound (retained)")
    _cdid_ur_nv("fixed_se","double",fse,"Fixed-share analytical SE (retained diagnostic)")
    _cdid_ur_sv("interval_kind",ik,"Reported interval kind")
    _cdid_ur_sv("se_kind",sk,"Reported SE kind")
    _cdid_ur_nv("available","byte",av,"Retained availability")
    _cdid_ur_sv("status",st,"Retained bin or overall status")
    _cdid_ur_nv("restricted","byte",rs,"Restricted target under na.rm")
    _cdid_ur_nv("n_intended","long",nint,"Intended cells")
    _cdid_ur_nv("n_used","long",nuse,"Used cells")
    _cdid_ur_nv("share_dependent","byte",sdep,"Share-dependent bin")
    _cdid_ur_nv("share_denominator","double",sden,"Share denominator")
    _cdid_ur_nv("overall_weight","double",ovw,"Weight of the bin in the overall summary")
    _cdid_ur_sv("role",role,"Retained bin role")
    _cdid_ur_sv("support",sup,"Retained bin support status")
    _cdid_ur_sv("reason",cdid_text_unavailable_v(why),"Availability reason")
    _cdid_ur_nv("cell","long",cell,"Original cell ledger position (membership rows)")
    _cdid_ur_nv("g","double",g,"Cohort of the member cell")
    _cdid_ur_nv("time","double",tt,"Calendar date of the member cell")
    _cdid_ur_nv("intended","byte",intd,"Member intended")
    _cdid_ur_nv("used","byte",usd,"Member used")
    _cdid_ur_nv("excluded","byte",exd,"Member excluded")
    _cdid_ur_nv("weight","double",wt,"Retained aggregation weight of the member")
    _cdid_ur_chars(V,"aggregate")
    st_global("_dta[cdid_target_definition]",A.target_definition)
    st_global("_dta[cdid_band_family]",A.band_family)
    st_global("_dta[cdid_inference_status]",A.inference_status)
    st_global("_dta[cdid_inference_reason]",A.inference_reason)
    st_global("_dta[cdid_selector_warning]",A.selector_warning)
    st_global("_dta[cdid_alpha]",_cdid_ur_num2s(A.alpha))
    st_global("_dta[cdid_critical]",_cdid_ur_num2s(A.critical))
    st_global("_dta[cdid_normal_critical]",_cdid_ur_num2s(A.normal_critical))
    st_global("_dta[cdid_na_rm]",_cdid_ur_num2s(A.na_rm))
    st_global("_dta[cdid_min_e]",_cdid_ur_num2s(A.min_e))
    st_global("_dta[cdid_max_e]",_cdid_ur_num2s(A.max_e))
    st_global("_dta[cdid_balance_e]",_cdid_ur_num2s(A.balance_e))
}

// ---------------------------------------------------------------- diagnostics table
void _cdid_ur_addn(struct cdid_ur_report scalar T, string scalar sec, string scalar lab,
    string scalar rname, real scalar v, string scalar why)
{
    T.section=T.section\sec; T.label=T.label\lab; T.rname=T.rname\rname
    T.value=T.value\v; T.isnum=T.isnum\1; T.text=T.text\(v>=. ? why : "")
}

void _cdid_ur_addt(struct cdid_ur_report scalar T, string scalar sec, string scalar lab,
    string scalar rname, string scalar s)
{
    T.section=T.section\sec; T.label=T.label\lab; T.rname=T.rname\rname
    T.value=T.value\.; T.isnum=T.isnum\0; T.text=T.text\s
}

void _cdid_ur_sett(struct cdid_ur_report scalar T, struct cdid_ur_view scalar V, string scalar sec,
    string scalar lab, string scalar key, string scalar rname)
{
    _cdid_ur_addt(T,sec,lab,rname,_cdid_ur_text(V,key,"not retained"))
}

void _cdid_ur_setn(struct cdid_ur_report scalar T, struct cdid_ur_view scalar V, string scalar sec,
    string scalar lab, string scalar key, string scalar rname)
{
    _cdid_ur_addn(T,sec,lab,rname,_cdid_ur_num(V,key,.),(_cdid_ur_hasnum(V,key) ? "retained as missing" : "not retained"))
}

string scalar _cdid_ur_tally(string colvector s)
{
    string colvector u
    string scalar out
    real scalar i
    if(rows(s)==0) return("none")
    u=uniqrows(s); out=""
    for(i=1;i<=rows(u);i++)
        out=out+(i>1 ? "; " : "")+(u[i]=="" ? "(unlabeled)" : u[i])+"="+_cdid_rb_key(sum(s:==u[i]))
    return(out)
}

real scalar _cdid_ur_count(real colvector x)
{
    return(missing(x) ? . : sum((x:==1)))
}

// Plain values shared with e(): standard errors, intervals, random numbers, inactive options.
string scalar _cdid_ur_vce(string scalar s)
{
    if(s=="analytic" | s=="source_analytic") return("analytic")
    if(strpos(s,"bootstrap")) return("bootstrap")
    if(strpos(s,"analytic")) return("analytic")
    return(s)
}

string scalar _cdid_ur_band(string scalar s)
{
    if(s=="simultaneous_band") return("simultaneous")
    return(s)
}

string scalar _cdid_ur_rng(string scalar s)
{
    if(s=="caller") return("Stata random-number generator (set seed); draws only for the multiplier bootstrap or simultaneous bands")
    if(s=="fixed") return("fixed multiplier signs supplied with the request")
    if(s=="loaded") return("results loaded from a file; no new multiplier draws")
    if(s=="not_used") return("not used: no multiplier draws")
    return(subinstr(s,"_"," "))
}

string scalar _cdid_ur_options(string scalar s)
{
    string rowvector a
    string scalar out, o
    real scalar i
    a=tokens(s); out=""
    for(i=1;i<=cols(a);i++) {
        o=a[i]
        if(o=="imputation_pooling") o="imppooling()"
        else if(o=="selection_model") o="selmodel()"
        else if(o=="selection_pooling") o="selpooling()"
        else if(o=="selection_xformla") o="selcovars()"
        else if(o=="outcome_pooling") o="outcomepooling()"
        out=out+(i>1 ? " " : "")+o
    }
    return(out)
}

struct cdid_ur_report scalar cdid_ur_report_build(struct cdid_ur_view scalar V)
{
    struct cdid_ur_report scalar T
    struct cdid_ur_cells scalar C
    struct cdid_ur_agg scalar A
    pointer(struct cdid_result_state scalar) scalar p
    real colvector post, pre, ref, av, x
    string colvector cs, hits
    string scalar sec, s
    real scalar N, K, ord, v, i, clu
    p=V.p; N=rows(p->cohort); K=rows(p->point); ord=(V.mode=="ordinary")
    T.section=T.label=T.rname=T.text=J(0,1,""); T.value=T.isnum=J(0,1,.)
    C=cdid_ur_cells_view(V)
    sec="Result"
    _cdid_ur_addt(T,sec,"Mode","mode",V.mode)
    _cdid_ur_addt(T,sec,"Current result","result_kind",V.kind)
    _cdid_ur_addt(T,sec,"Aggregation type","agg_type",(V.kind=="aggregation" ? V.agg_type : "none"))
    _cdid_ur_addt(T,sec,"Stored results format","results_format",p->schema)
    _cdid_ur_addt(T,sec,"Unit identifier type","id_type",p->id_kind)
    _cdid_ur_addt(T,sec,"Cluster identifier type","cluster_id_type",p->cluster_kind)

    sec="Units"
    cs=_cdid_ur_cohort_status(V)
    _cdid_ur_addn(T,sec,"Units","N_units",N,"")
    _cdid_ur_addn(T,sec,"Treated-cohort units","N_treated",sum(cs:=="treated"),"")
    _cdid_ur_addn(T,sec,"Never-treated units","N_never_treated",sum(cs:=="never_treated"),"")
    _cdid_ur_addn(T,sec,"Unknown-cohort units","N_unknown_cohort",sum(cs:=="unknown"),"")
    if(_cdid_ur_hasu(V,"supplied_rows")) {
        x=_cdid_ur_ucol(V,"supplied_rows"); v=(missing(x) ? . : sum(x:>0))
    }
    else v=.
    _cdid_ur_addn(T,sec,"Units with data rows","N_with_rows",v,"not retained")
    v=(_cdid_ur_hasu(V,"frame_only") ? _cdid_ur_count(_cdid_ur_ucol(V,"frame_only")) : .)
    _cdid_ur_addn(T,sec,"Units only in frame()","N_frame_only",v,"not retained")
    v=(_cdid_ur_hasnum(V,"N_B") ? _cdid_ur_num(V,"N_B",.) : (_cdid_ur_hasu(V,"B") ? _cdid_ur_count(_cdid_ur_ucol(V,"B")) : .))
    _cdid_ur_addn(T,sec,"Units with an observed outcome","N_with_outcome",v,"not retained")
    v=(_cdid_ur_hasnum(V,"N_A") ? _cdid_ur_num(V,"N_A",.) : (_cdid_ur_hasu(V,"A") ? _cdid_ur_count(_cdid_ur_ucol(V,"A")) : .))
    _cdid_ur_addn(T,sec,"Units with two observed outcomes","N_with_two_outcomes",v,"not retained")
    if(p->cluster_kind=="numeric") v=rows(uniqrows(p->cluster_numeric_id))
    else if(p->cluster_kind=="string") v=rows(uniqrows(p->cluster_string_id))
    else v=.
    _cdid_ur_addn(T,sec,"Distinct cluster identifiers","N_cluster_ids",v,"unclustered")
    // First-calendar masses may legitimately be absent. Mata sum()
    // omits missings, so guard the full-population summaries explicitly.
    v=(N>0 & !missing(p->weight))
    _cdid_ur_addn(T,sec,"Unit weight sum","weight_sum",(v ? sum(p->weight) : .),(N ? "retained weights include missing values" : "no units"))
    _cdid_ur_addn(T,sec,"Unit weight minimum","weight_min",(v ? min(p->weight) : .),(N ? "retained weights include missing values" : "no units"))
    _cdid_ur_addn(T,sec,"Unit weight maximum","weight_max",(v ? max(p->weight) : .),(N ? "retained weights include missing values" : "no units"))

    sec="ATT(g,t) cells"
    post=(p->role:==3); pre=(p->role:==2); ref=(p->role:==1); av=(p->available:==1)
    _cdid_ur_addn(T,sec,"Cells","N_cells",K,"")
    _cdid_ur_addn(T,sec,"Post-treatment cells","N_post",sum(post),"")
    _cdid_ur_addn(T,sec,"Pre-treatment cells","N_pre",sum(pre),"")
    _cdid_ur_addn(T,sec,"Reference cells","N_reference",sum(ref),"")
    _cdid_ur_addn(T,sec,"Cells with another role","N_other_role",K-sum(post:|pre:|ref),"")
    _cdid_ur_addn(T,sec,"Available cells","N_available",sum(av),"")
    _cdid_ur_addn(T,sec,"Available post-treatment cells","N_post_available",sum(post:&av),"")
    _cdid_ur_addn(T,sec,"Available pre-treatment cells","N_pre_available",sum(pre:&av),"")
    _cdid_ur_addn(T,sec,"Available reference cells","N_reference_available",sum(ref:&av),"")
    _cdid_ur_addn(T,sec,"Unavailable cells","N_unavailable",K-sum(av),"")
    _cdid_ur_addn(T,sec,"Cells with a standard error","N_with_se",sum(C.se:<.),"")

    sec="Identification and GMM"
    if(!ord) {
        _cdid_ur_addt(T,sec,"Post-treatment cell support","post_support",_cdid_ur_tally(select(C.support,post)))
        _cdid_ur_addt(T,sec,"Pre-treatment cell support","pre_support",_cdid_ur_tally(select(C.support,pre)))
        _cdid_ur_sett(T,V,sec,"GMM weighting","gmm_weighting","gmm_weighting")
        _cdid_ur_sett(T,V,sec,"GMM covariance","gmm_covariance","gmm_covariance")
        _cdid_ur_sett(T,V,sec,"GMM fit status","GMM_status","gmm_status")
        _cdid_ur_sett(T,V,sec,"GMM covariance regularization","GMM_covariance_handler","gmm_regularization")
        _cdid_ur_sett(T,V,sec,"GMM regularization reason","GMM_trigger_reason","gmm_regularization_reason")
        _cdid_ur_setn(T,V,sec,"GMM covariance rank","GMM_covariance_rank","gmm_covariance_rank")
        _cdid_ur_setn(T,V,sec,"GMM empirical IF covariance rank bound","GMM_sample_rank_upper","gmm_rank_bound")
        _cdid_ur_setn(T,V,sec,"GMM criterion","GMM_criterion","gmm_criterion")
        _cdid_ur_setn(T,V,sec,"GMM Tikhonov lambda","GMM_lambda","gmm_lambda")
        _cdid_ur_setn(T,V,sec,"GMM eigenvalue ratio (min/max)","GMM_ridge_rho","gmm_eigenvalue_ratio")
        _cdid_ur_setn(T,V,sec,"GMM covariance trace","GMM_trace_B","gmm_covariance_trace")
        _cdid_ur_setn(T,V,sec,"Eigenvalue-ratio threshold, condition()","condition_trigger","condition")
        _cdid_ur_setn(T,V,sec,"Relative ridge, ridge()","ridge_relative","ridge")
        _cdid_ur_addt(T,sec,"","","GMM covariance rank and conditioning are distinct from the identification of the cells, which the support statuses report.")
    }
    else _cdid_ur_addt(T,sec,"GMM fit status","gmm_status","not applicable to ordinary results")

    sec="Model and selection settings"
    if(!ord) {
        _cdid_ur_sett(T,V,sec,"Method","method","method")
        _cdid_ur_sett(T,V,sec,"Covariate design","design","covariate_design")
        _cdid_ur_sett(T,V,sec,"Covariate design columns","covariates","covariate_columns")
        _cdid_ur_sett(T,V,sec,"Selection","selection","selection")
        _cdid_ur_sett(T,V,sec,"Selection assumption","selection_canonical","selection_canonical")
        _cdid_ur_sett(T,V,sec,"Selection model","selection_model","selection_model")
        _cdid_ur_sett(T,V,sec,"Selection pooling","selection_pooling","selection_pooling")
        _cdid_ur_sett(T,V,sec,"Selection universe","selection_universe","selection_universe")
        _cdid_ur_sett(T,V,sec,"Selection estimation method","selection_est_method","selection_est_method")
        _cdid_ur_sett(T,V,sec,"Imputation pooling","imputation_pooling","imputation_pooling")
        _cdid_ur_sett(T,V,sec,"Outcome pooling","outcome_pooling","outcome_pooling")
        _cdid_ur_sett(T,V,sec,"Selection covariates","selection_covariates","selection_covariates")
        s=_cdid_ur_text(V,"inactive_arguments","not retained")
        _cdid_ur_addt(T,sec,"Inactive options","inactive_options",(s=="not retained" ? s : _cdid_ur_options(s)))
        _cdid_ur_sett(T,V,sec,"Nuisance-model fit status","provider_status","nuisance_fit_status")
        _cdid_ur_sett(T,V,sec,"Covariate profile","profile_kind","covariate_profile")
        _cdid_ur_sett(T,V,sec,"Observation-weighted static check","selected_static_gate","obs_weighted_check")
        _cdid_ur_sett(T,V,sec,"Observation-weighted check reason","selected_static_gate_reason","obs_weighted_check_reason")
        _cdid_ur_sett(T,V,sec,"Baseline static check","baseline_static_gate","baseline_check")
        _cdid_ur_sett(T,V,sec,"Baseline check reason","baseline_static_gate_reason","baseline_check_reason")
        _cdid_ur_setn(T,V,sec,"Treatment fits","n_treatment_fits","n_treatment_fits")
        _cdid_ur_setn(T,V,sec,"Pooled outcome-regression fits","n_pooled_OR_fits","n_pooled_outcome_fits")
        _cdid_ur_setn(T,V,sec,"Selection fits","n_selection_fits","n_selection_fits")
        _cdid_ur_setn(T,V,sec,"Imputation fits","n_imputation_fits","n_imputation_fits")
    }
    else {
        _cdid_ur_sett(T,V,sec,"Method","ord.method","method")
        _cdid_ur_sett(T,V,sec,"Data layout","ord.panel_mode","data_layout")
        _cdid_ur_sett(T,V,sec,"csdid computation","ord.compute_path","csdid_computation")
        _cdid_ur_sett(T,V,sec,"Cluster variable","ord.cluster_name","clustvar")
    }

    sec="Reference dates and weights"
    if(!ord) {
        _cdid_ur_sett(T,V,sec,"Base period","base_period","base_period")
        _cdid_ur_sett(T,V,sec,"Control group","control_group","control_group")
        _cdid_ur_sett(T,V,sec,"Parallel-trends assumption","pt_assumption","pt_assumption")
        _cdid_ur_setn(T,V,sec,"Parallel-trends start","pt_start","pt_start")
        _cdid_ur_setn(T,V,sec,"Anticipation","anticipation","anticipation")
        _cdid_ur_sett(T,V,sec,"Weight source","weight_source","weight_source")
        _cdid_ur_sett(T,V,sec,"Weight reduction","weight_reduction","weight_reduction")
        _cdid_ur_sett(T,V,sec,"Fixed weights requested","fix_weights_requested","fix_weights_requested")
        _cdid_ur_sett(T,V,sec,"Cohort shares for aggregation","aggregation_share_law","cohort_shares")
        _cdid_ur_sett(T,V,sec,"Cohort share status","share_status","cohort_share_status")
    }
    else {
        _cdid_ur_sett(T,V,sec,"Base period","ord.base_period","base_period")
        _cdid_ur_sett(T,V,sec,"Control group","ord.control_group","control_group")
        _cdid_ur_setn(T,V,sec,"Anticipation","ord.anticipation","anticipation")
        _cdid_ur_sett(T,V,sec,"Fixed weights requested","ord.fix_weights","fix_weights_requested")
        _cdid_ur_sett(T,V,sec,"Weight normalization","ord.source_weight_law","weight_norm")
        _cdid_ur_setn(T,V,sec,"Weight normalization factor","ord.weight_normalization_factor","weight_norm_factor")
    }

    sec="Numerical safeguards and resources"
    hits=J(0,1,"")
    for(i=1;i<=rows(p->text_settings);i++) if(strpos(strlower(p->text_settings[i,1]),"clip"))
        hits=hits\(p->text_settings[i,1]+"="+p->text_settings[i,2])
    for(i=1;i<=rows(p->numeric_setting_names);i++) if(strpos(strlower(p->numeric_setting_names[i]),"clip"))
        hits=hits\(p->numeric_setting_names[i]+"="+_cdid_ur_num2s(p->numeric_settings[i]))
    _cdid_ur_addt(T,sec,"Probability clipping","clipping",(rows(hits) ? invtokens(hits',"; ") : "not retained in this result"))
    if(!ord) {
        _cdid_ur_setn(T,V,sec,"Workspace in bytes, workspace()","numeric_workspace_bytes","workspace")
        _cdid_ur_setn(T,V,sec,"Rows per chunk, chunk()","chunk_rows","chunk")
        _cdid_ur_setn(T,V,sec,"Work limit, workcap()","counted_work_limit","workcap")
    }

    sec="Pretrends"
    s=(ord ? "ord.cells." : "cdid_bridge.")
    _cdid_ur_sett(T,V,sec,"Pretrend Wald status",s+"wald_status","wald_status")
    _cdid_ur_addt(T,sec,"Pretrend Wald reason","wald_reason",cdid_text_reason(_cdid_ur_text(V,s+"wald_reason","")))
    _cdid_ur_setn(T,V,sec,"Wald degrees of freedom",s+"wald_df","wald_df")
    _cdid_ur_setn(T,V,sec,"Wald statistic",s+"wald_statistic","wald_statistic")
    _cdid_ur_setn(T,V,sec,"Wald p-value",s+"wald_p","wald_p")
    _cdid_ur_addt(T,sec,"Conditional pretrend test","conditional_test","not retained in this result")

    sec="Cell inference"
    _cdid_ur_addt(T,sec,"Inference status","inference_status",C.status)
    _cdid_ur_addt(T,sec,"Intervals","band",_cdid_ur_band(C.interval_kind))
    _cdid_ur_addt(T,sec,"Standard errors","vce",_cdid_ur_vce(C.se_kind))
    _cdid_ur_addn(T,sec,"Point estimates only","point_only",C.point_only,"not recorded")
    _cdid_ur_addn(T,sec,"Alpha","alpha",C.alpha,"not retained")
    s=(ord ? "ord.cells." : "cdid_bridge.cells.")
    _cdid_ur_setn(T,V,sec,"bstrap()",s+"bstrap","bstrap")
    _cdid_ur_setn(T,V,sec,"cband()",s+"cband","cband")
    _cdid_ur_setn(T,V,sec,"Multiplier draws requested",s+"B","B_requested")
    if(!ord) {
        clu=_cdid_ur_num(V,"cdid_bridge.cells.clustered",.)
        _cdid_ur_addn(T,sec,"Clustered inference","clustered",clu,(_cdid_ur_hasnum(V,"cdid_bridge.cells.clustered") ? "retained as missing" : "not retained"))
        v=_cdid_ur_num(V,"cdid_bridge.cells.C",.)
    }
    else {
        clu=(p->cluster_kind!="none")
        _cdid_ur_addn(T,sec,"Clustered inference","clustered",clu,"")
        v=_cdid_ur_num(V,"ord.C",.)
    }
    _cdid_ur_addn(T,sec,"Number of clusters","N_clust",(clu==1 ? v : .),(clu==1 ? "not retained" : "unclustered"))
    if(!ord) _cdid_ur_addt(T,sec,"Random-number source","rng",_cdid_ur_rng(_cdid_ur_text(V,"cdid_bridge.rng.mode","not retained")))

    if(V.kind=="aggregation") {
        A=cdid_ur_agg_view(V)
        sec="Current aggregation"
        _cdid_ur_addn(T,sec,"Aggregates","N_aggregates",A.M,"")
        _cdid_ur_addn(T,sec,"Available aggregates","N_aggregates_available",_cdid_ur_count(A.available),"")
        _cdid_ur_addn(T,sec,"Aggregates restricted by narm","N_aggregates_narm",(ord ? . : _cdid_ur_count(A.restricted)),"not recorded")
        _cdid_ur_addt(T,sec,"Overall status","overall_status",A.overall_status)
        _cdid_ur_addt(T,sec,"Overall reason","overall_reason",cdid_text_unavailable(A.overall_reason))
        _cdid_ur_addt(T,sec,"Intervals of the aggregates","agg_band",_cdid_ur_band(A.interval_kind))
        _cdid_ur_addt(T,sec,"Standard errors of the aggregates","agg_vce",_cdid_ur_vce(A.se_kind))
        _cdid_ur_addt(T,sec,"Standard error of the overall summary","overall_vce",_cdid_ur_vce(A.overall_se_kind))
        _cdid_ur_addn(T,sec,"narm","na_rm",A.na_rm,"not retained")
        _cdid_ur_addn(T,sec,"mine()","min_e",A.min_e,"not set")
        _cdid_ur_addn(T,sec,"maxe()","max_e",A.max_e,"not set")
        _cdid_ur_addn(T,sec,"balancee()","balance_e",A.balance_e,"not set")
        _cdid_ur_addt(T,sec,"Event-time selection note","selector_warning",A.selector_warning)
    }
    return(T)
}

void _cdid_ur_show(struct cdid_ur_report scalar T)
{
    real scalar i, w
    string scalar last, v
    last=""; w=0
    for(i=1;i<=rows(T.label);i++) w=max((w,strlen(T.label[i])))
    w=max((44,w+5))
    for(i=1;i<=rows(T.label);i++) {
        if(T.section[i]!=last) {
            printf("\n{txt}%s\n",_cdid_ur_smcl(T.section[i])); last=T.section[i]
        }
        if(T.label[i]=="") {
            printf("{txt}{p 2 2 2}%s{p_end}\n",_cdid_ur_smcl(T.text[i]))
            continue
        }
        if(T.isnum[i]) v=((T.value[i]<. | T.text[i]=="") ? _cdid_ur_num2s(T.value[i]) : T.text[i])
        else v=(T.text[i]=="" ? "(empty)" : T.text[i])
        if(T.rname[i]=="cohort_shares" & v=="first_calendar_masses") v="cohort weights at the first calendar date"
        printf("{txt}  %s{col %g}{res}%s\n",_cdid_ur_smcl(T.label[i]),w,_cdid_ur_smcl(v))
    }
}

void _cdid_ur_post(struct cdid_ur_report scalar T)
{
    real scalar i
    st_rclear()
    for(i=1;i<=rows(T.label);i++) {
        if(T.rname[i]=="") continue
        if(T.isnum[i]) {
            st_numscalar("r("+T.rname[i]+")",T.value[i])
            if(T.value[i]>=. & T.text[i]!="") st_global("r("+T.rname[i]+"_status)",T.text[i])
        }
        else st_global("r("+T.rname[i]+")",T.text[i])
    }
}

void _cdid_ur_detail(struct cdid_ur_view scalar V)
{
    struct cdid_ur_cells scalar C
    struct cdid_ur_agg scalar A
    real scalar j, any_un
    C=cdid_ur_cells_view(V); any_un=0
    printf("\n{txt}Unavailable ATT(g,t) cells\n")
    for(j=1;j<=rows(C.estimate);j++) if(C.available[j]!=1) {
        any_un=1
        printf("{txt}  cell %g (cohort %s, date %s, %s): {res}%s\n",j,_cdid_ur_num2s(C.g[j]),_cdid_ur_num2s(C.time[j]),
            _cdid_ur_smcl(cdid_text_role(C.role[j])),_cdid_ur_smcl(C.reason[j]=="" ? "(no reason stored)" : cdid_text_unavailable(C.reason[j])))
    }
    if(!any_un) printf("{txt}  none\n")
    if(V.kind=="aggregation") {
        A=cdid_ur_agg_view(V); any_un=0
        printf("\n{txt}Unavailable aggregates\n")
        for(j=1;j<=A.M;j++) if(A.available[j]!=1) {
            any_un=1
            printf("{txt}  %s: {res}%s\n",_cdid_ur_smcl(cdid_text_bin(A.label[j])),_cdid_ur_smcl(A.reason[j]=="" ? "(no reason stored)" : cdid_text_unavailable(A.reason[j])))
        }
        if(!any_un) printf("{txt}  none\n")
    }
}

void cdid_ur_diagnostics(real scalar detail, real scalar matdim)
{
    struct cdid_ur_view scalar V
    struct cdid_ur_report scalar T
    struct cdid_ur_cells scalar C
    struct cdid_ur_agg scalar A
    real scalar K, M
    real matrix X
    V=cdid_ur_open(); T=cdid_ur_report_build(V); C=cdid_ur_cells_view(V)
    printf("\n{txt}Diagnostics of the current cdid results (from the stored results; no data are read)\n")
    _cdid_ur_show(T)
    if(detail) _cdid_ur_detail(V)
    _cdid_ur_post(T)
    K=rows(C.estimate)
    if(K>0 & K<=matdim) {
        X=(C.g,C.time,C.role_code,C.available,C.estimate,C.se,C.lower,C.upper)
        st_matrix("r(cells)",X)
        st_matrixcolstripe("r(cells)",(J(8,1,""),("g"\"time"\"role_code"\"available"\"estimate"\"se"\"lower"\"upper")))
    }
    else if(K>matdim) st_global("r(cells_matrix)","not returned: the number of cells exceeds c(max_matdim); use cdid_export")
    if(V.kind=="aggregation") {
        A=cdid_ur_agg_view(V); M=A.M
        if(M>0 & M<=matdim) {
            X=(A.key,A.available,A.estimate,A.se,A.lower,A.upper)
            st_matrix("r(aggregate)",X)
            st_matrixcolstripe("r(aggregate)",(J(6,1,""),("key"\"available"\"estimate"\"se"\"lower"\"upper")))
        }
        else if(M>matdim) st_global("r(aggregate_matrix)","not returned: the bins exceed c(max_matdim); use cdid_export")
    }
}

// ---------------------------------------------------------------- summary export
void cdid_ur_export_summary(struct cdid_ur_view scalar V)
{
    pointer(struct cdid_result_state scalar) scalar p
    struct cdid_ur_report scalar T
    string colvector cat, key, txt, nm
    real colvector val, coh
    real scalar i, j, cc, n
    real matrix em
    p=V.p; T=cdid_ur_report_build(V)
    cat=J(rows(T.label),1,"diagnostic"); key=T.rname
    for(i=1;i<=rows(key);i++) if(key[i]=="") key[i]=(T.label[i]=="" ? "note" : T.label[i])
    txt=T.text; val=T.value; coh=J(rows(T.label),1,.)
    cat=cat\J(3,1,"state"); key=key\("schema"\"id_kind"\"cluster_kind")
    txt=txt\(p->schema\p->id_kind\p->cluster_kind); val=val\J(3,1,.); coh=coh\J(3,1,.)
    n=rows(p->text_settings)
    if(n) {
        cat=cat\J(n,1,"text_setting"); key=key\p->text_settings[,1]; txt=txt\p->text_settings[,2]
        val=val\J(n,1,.); coh=coh\J(n,1,.)
    }
    n=rows(p->numeric_setting_names)
    if(n) {
        cat=cat\J(n,1,"numeric_setting"); key=key\p->numeric_setting_names; txt=txt\J(n,1,"")
        val=val\p->numeric_settings; coh=coh\J(n,1,.)
    }
    n=rows(p->source_versions)
    if(n) {
        cat=cat\J(n,1,"source_version"); key=key\p->source_versions[,1]; txt=txt\p->source_versions[,2]
        val=val\J(n,1,.); coh=coh\J(n,1,.)
    }
    cc=asarray(_cdid_ur_index(p->share_columns'),"cohort")
    for(i=1;i<=rows(p->shares);i++) for(j=1;j<=cols(p->shares);j++) {
        cat=cat\"cohort_share"; key=key\p->share_columns[j]; txt=txt\""
        val=val\p->shares[i,j]; coh=coh\(cc ? p->shares[i,cc] : .)
    }
    nm=st_dir("e()","macro","*")
    for(i=1;i<=rows(nm);i++) {
        cat=cat\"e_macro"; key=key\nm[i]; txt=txt\st_global("e("+nm[i]+")"); val=val\.; coh=coh\.
    }
    nm=st_dir("e()","numscalar","*")
    for(i=1;i<=rows(nm);i++) {
        cat=cat\"e_scalar"; key=key\nm[i]; txt=txt\""; val=val\st_numscalar("e("+nm[i]+")"); coh=coh\.
    }
    nm=st_dir("e()","matrix","*")
    for(i=1;i<=rows(nm);i++) {
        em=st_matrix("e("+nm[i]+")")
        cat=cat\"e_matrix"; key=key\nm[i]; txt=txt\(_cdid_rb_key(rows(em))+" x "+_cdid_rb_key(cols(em))); val=val\.; coh=coh\.
    }
    _cdid_ur_newobs(rows(cat))
    _cdid_ur_sv("category",cat,"Record category")
    _cdid_ur_sv("key",key,"Retained key or r() name")
    _cdid_ur_sv("text",txt,"Text value")
    _cdid_ur_nv("value","double",val,"Numeric value")
    _cdid_ur_nv("cohort","double",coh,"Cohort (cohort_share rows)")
    _cdid_ur_chars(V,"summary")
}

// ---------------------------------------------------------------- public entries
void cdid_ur_check(string scalar purpose)
{
    struct cdid_ur_view scalar V
    pointer(struct cdid_result_state scalar) scalar p
    V=cdid_ur_open(); p=V.p
    if(purpose=="aggregate" & V.kind!="aggregation")
        _cdid_ur_fail("type(aggregate) requires an aggregation; run estat simple, group, dynamic or calendar first (type(cells) saves the group-time cells)",498)
    if(purpose=="influence" & rows(p->IF)==0)
        _cdid_ur_fail("the current results have no influence functions (noinffunc); rerun cdid without noinffunc",459)
    st_local("ur_mode",V.mode); st_local("ur_kind",V.kind); st_local("ur_agg_type",V.agg_type)
}

void cdid_ur_export(string scalar type)
{
    struct cdid_ur_view scalar V
    V=cdid_ur_open()
    if(type=="cells") cdid_ur_export_cells(V)
    else if(type=="aggregate") cdid_ur_export_agg(V)
    else if(type=="units") cdid_ur_export_units(V)
    else if(type=="influence") cdid_ur_export_influence(V)
    else if(type=="summary") cdid_ur_export_summary(V)
    else _cdid_ur_fail("type() must be cells, aggregate, units, influence or summary",198)
    st_local("ur_rows",_cdid_rb_key(st_nobs()))
}

// ---------------------------------------------------------------- plot data
// Option-name classifier for native graph options passed through by cdid_plot.
string scalar cdid_ur_optname(string scalar token)
{
    string scalar s
    string colvector full
    real colvector minlen
    real scalar k, at
    s=strlower(strtrim(token)); at=strpos(s,"(")
    if(at) s=strtrim(substr(s,1,at-1))
    full=("title"\"subtitle"\"note"\"caption"\"xtitle"\"ytitle"\"legend"\"xlabel"\"by"\"name")
    minlen=(2\3\4\3\3\3\3\3\2\4)
    for(k=1;k<=rows(full);k++)
        if(strlen(s)>=minlen[k] & s==substr(full[k],1,strlen(s))) return(full[k])
    return("")
}

// Split "file[, options]" respecting simple and compound double quotes.
void cdid_ur_split(string scalar specname, string scalar fileloc, string scalar optloc)
{
    string scalar s, f, rest
    real scalar at, i, q
    s=strtrim(st_local(specname)); f=""; rest=""
    if(substr(s,1,2)==char(96)+char(34)) {
        at=strpos(s,char(34)+char(39))
        if(at==0) _cdid_ur_fail("unbalanced compound quotes in a file specification",198)
        f=substr(s,3,at-3); rest=strtrim(substr(s,at+2,.))
    }
    else if(substr(s,1,1)==char(34)) {
        at=strpos(substr(s,2,.),char(34))
        if(at==0) _cdid_ur_fail("unbalanced quotes in a file specification",198)
        f=substr(s,2,at-1); rest=strtrim(substr(s,at+2,.))
    }
    else {
        at=strpos(s,",")
        if(at) {
            f=strtrim(substr(s,1,at-1)); rest=strtrim(substr(s,at,.))
        }
        else f=s
    }
    if(rest!="") {
        if(substr(rest,1,1)!=",") _cdid_ur_fail("expected a comma after the file name in: "+s,198)
        rest=strtrim(substr(rest,2,.))
    }
    q=0
    for(i=1;i<=strlen(f);i++) if(substr(f,i,1)==char(34)) q=1
    if(q) _cdid_ur_fail("file names may not contain double quotes",198)
    st_local(fileloc,f); st_local(optloc,rest)
}

void cdid_ur_plot_data(real scalar want_cells, string scalar groups)
{
    struct cdid_ur_view scalar V
    struct cdid_ur_cells scalar C
    struct cdid_ur_agg scalar A
    real colvector req, sel, x, grp, tt, ev, est, lo, hi, av, keepg, present, drawn, lv, capped
    string colvector series, lab, why, svals
    string scalar ptype, xt, yt, ikind, skind, note, gnote, byvar
    real scalar n, i, alpha, usecells, nomit, nocap
    V=cdid_ur_open()
    usecells=(V.kind=="cells" | want_cells)
    if(!usecells & V.agg_type=="simple")
        _cdid_ur_fail("estat simple gives a single overall effect; specify cells to graph the group-time estimates, or run estat group, dynamic or calendar",498)
    req=(strtrim(groups)=="" ? J(0,1,.) : strtoreal(tokens(groups))')
    gnote=""
    if(usecells) {
        C=cdid_ur_cells_view(V); n=rows(C.estimate)
        x=C.time; grp=C.g; tt=C.time; ev=C.event; est=C.estimate; lo=C.lower; hi=C.upper; av=C.available
        series=J(n,1,"")
        for(i=1;i<=n;i++) {
            if(C.role_code[i]==2) series[i]="Pre"
            else if(C.role_code[i]==1) series[i]="Reference"
            else if(C.role_code[i]==3) series[i]="Post"
            else series[i]=(C.time[i]<C.g[i] ? "Pre" : "Post")
        }
        lab=C.target_id; why=C.reason; ikind=C.interval_kind; skind=C.se_kind; alpha=C.alpha
        ptype="cells"; xt="Time"; yt="ATT(g,t)"; byvar="group"; keepg=C.g
    }
    else {
        A=cdid_ur_agg_view(V); n=A.M
        x=A.key; est=A.estimate; lo=A.lower; hi=A.upper; av=A.available; lab=A.label; why=A.reason
        grp=tt=ev=J(n,1,.)
        if(A.type=="group") grp=A.key
        else if(A.type=="calendar") tt=A.key
        else ev=A.key
        series=J(n,1,"Post")
        for(i=1;i<=n;i++) {
            if(A.role[i]=="reference_aggregate") series[i]="Reference"
            else if(A.type=="dynamic" & A.key[i]<0) series[i]="Pre"
        }
        ikind=A.interval_kind; skind=A.se_kind; alpha=A.alpha
        ptype="aggregate_"+A.type; yt="ATT"; byvar=""
        xt=(A.type=="dynamic" ? "Event time (t - g)" : (A.type=="group" ? "Cohort (first treated)" : "Time"))
        keepg=(A.type=="group" ? A.key : J(n,1,.))
        if(rows(req) & A.type!="group") {
            gnote="group() applies to cell and group-aggregation plots; ignored for the "+A.type+" aggregation"
            req=J(0,1,.)
        }
    }
    sel=J(n,1,1)
    if(rows(req)) {
        present=J(rows(req),1,0)
        for(i=1;i<=rows(req);i++) present[i]=any(keepg:==req[i])
        if(all(present)) {
            sel=J(n,1,0)
            for(i=1;i<=rows(req);i++) sel=sel:|(keepg:==req[i])
        }
        else gnote="some cohorts in group() are not in the results; all cohorts are shown"
    }
    sel=vec(selectindex(sel)); n=rows(sel)
    if(n==0) _cdid_ur_fail("nothing to graph: the results contain no estimates",498)
    x=x[sel]; grp=grp[sel]; tt=tt[sel]; ev=ev[sel]; est=est[sel]; lo=lo[sel]; hi=hi[sel]; av=av[sel]
    series=series[sel]; lab=lab[sel]; why=why[sel]
    drawn=(av:==1):&(est:<.); capped=drawn:&(lo:<.):&(hi:<.)
    if(!any(drawn)) {
        errprintf("nothing to plot: every selected target is unavailable\n")
        for(i=1;i<=min((n,10));i++) errprintf("  %s: %s\n",(lab[i]=="" ? "target "+_cdid_rb_key(sel[i]) : lab[i]),(why[i]=="" ? "(no retained reason)" : cdid_text_unavailable(why[i])))
        exit(498)
    }
    _cdid_ur_newobs(n)
    _cdid_ur_sv("plot_type",J(n,1,ptype),"Plot type")
    _cdid_ur_sv("series",series,"Pre, Post or Reference")
    _cdid_ur_nv("x","double",x,"Horizontal axis value")
    _cdid_ur_sv("x_label",strtrim(strofreal(x,"%21.0g")),"Horizontal axis label")
    _cdid_ur_nv("estimate","double",est,"Retained point estimate")
    _cdid_ur_nv("ci_low","double",lo,"Retained reported interval lower bound")
    _cdid_ur_nv("ci_high","double",hi,"Retained reported interval upper bound")
    _cdid_ur_nv("group","double",grp,"Cohort")
    _cdid_ur_nv("time","double",tt,"Calendar date")
    _cdid_ur_nv("event_time","double",ev,"Event time t - g")
    _cdid_ur_nv("available","byte",av,"Retained availability")
    _cdid_ur_nv("drawn","byte",drawn,"Drawn on the graph")
    _cdid_ur_nv("target","long",sel,"Position in the retained ledger or bins")
    _cdid_ur_sv("label",lab,"Retained target label")
    _cdid_ur_sv("interval_kind",J(n,1,ikind),"Reported interval kind")
    _cdid_ur_sv("reason",cdid_text_unavailable_v(why),"Availability reason")
    svals=("Pre"\"Post"\"Reference")
    for(i=1;i<=3;i++) {
        st_local("ur_n_"+svals[i],_cdid_rb_key(sum(drawn:&(series:==svals[i]))))
        st_local("ur_cap_"+svals[i],_cdid_rb_key(sum(capped:&(series:==svals[i]))))
    }
    nomit=sum(!drawn); nocap=sum(drawn:&!capped)
    if(ikind=="none" | alpha>=.) note="Point estimates only; no intervals retained."
    else note="Bars: "+strtrim(strofreal(100*(1-alpha),"%9.0g"))+"% "+
        (ikind=="simultaneous_band" ? "simultaneous band" : "pointwise intervals")+
        (skind=="multiplier_bootstrap_iqr" ? " (multiplier bootstrap IQR standard errors)." : ".")
    if(nomit) note=note+" Not drawn (unavailable): "+_cdid_rb_key(nomit)+"."
    if(nocap & ikind!="none") note=note+" Drawn without a retained interval: "+_cdid_rb_key(nocap)+"."
    lv=uniqrows(select(x,drawn))
    st_local("ur_xlevels",(rows(lv)<=20 ? invtokens(strtrim(strofreal(lv',"%21.0g"))) : ""))
    st_local("ur_ptype",ptype); st_local("ur_xtitle",xt); st_local("ur_ytitle",yt)
    st_local("ur_byvar",byvar); st_local("ur_note",note); st_local("ur_gnote",gnote)
    st_local("ur_n_drawn",_cdid_rb_key(sum(drawn))); st_local("ur_n_omitted",_cdid_rb_key(nomit))
    st_local("ur_rows",_cdid_rb_key(n))
}

end
