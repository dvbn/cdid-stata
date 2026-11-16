*! Ordinary result-container translation only. No chained inference dispatch.
version 19.5
mata:
string scalar cdid_ord_post_stage(string scalar token,string scalar limits_name,
    string scalar posting,string scalar bname,string scalar vname,string scalar oname,
    string scalar sourcev,string scalar sourceo,string scalar sename,string scalar aname,
    string scalar ciname,string scalar piname,string scalar mapname,real scalar matrix_limit)
{
    struct cdid_result_state scalar P
    struct cdid_state_limits scalar limits
    pointer() scalar ptr
    real colvector point,se,analytic,av
    real matrix V,Omega,sourceOmega,ci,pi,map
    string matrix stripe
    string scalar kind,prefix
    real scalar j,hasV,pointonly,C,clustered
    ptr=findexternal(limits_name)
    if(ptr==NULL) _cdid_ord_fail("Explicit state limits disappeared before posting.",111)
    limits=*ptr;P=cdid_state_get(token);cdid_state_validate(P,limits)
    if(_cdid_ord_text(P,"ord.mode")!="ordinary") _cdid_ord_fail("Ordinary posting requires an ordinary payload.",498)
    if(posting!="matrix" & posting!="state") _cdid_ord_fail("Posting must explicitly request matrix or state.",198)
    point=_cdid_ord_matrix(P,"ord.report.point");se=_cdid_ord_matrix(P,"ord.report.se")
    analytic=_cdid_ord_matrix(P,"ord.report.analytic_se");ci=_cdid_ord_matrix(P,"ord.report.ci")
    pi=_cdid_ord_matrix(P,"ord.report.point_ci");V=_cdid_ord_matrix(P,"ord.report.V");Omega=_cdid_ord_matrix(P,"ord.report.Omega")
    kind=_cdid_ord_text(P,"ord.report.kind");pointonly=_cdid_ord_number(P,"ord.report.point_only")
    av=(posting=="matrix" ? _cdid_ord_ix(point:<.) : J(0,1,.))
    if(kind=="ordinary_cells") av=(posting=="matrix" ? _cdid_ord_ix((point:<.):&P.available) : J(0,1,.))
    if(matrix_limit<1 | matrix_limit>=. | matrix_limit!=floor(matrix_limit)) _cdid_ord_fail("Current Stata matrix limit must be explicit.",498)
    if(rows(av)>matrix_limit) _cdid_ord_fail("Mapped ordinary coefficient vector exceeds standard matrix size; explicit state posting retains the complete result.",908)
    hasV=(!pointonly & rows(V)>0 & rows(av)>0)
    map=(rows(av) ? ((1::rows(av)),av) : J(0,2,.));stripe=J(rows(av),2,"")
    prefix=(kind=="ordinary_cells" ? "cell_" : "agg_")
    for(j=1;j<=rows(av);j++) stripe[j,2]=prefix+_cdid_ord_key(av[j])
    if(rows(av)) {
        st_matrix(bname,point[av]');st_matrixcolstripe(bname,stripe)
        st_matrix(sename,se[av]');st_matrix(aname,analytic[av]')
        st_matrix(ciname,ci[av,]');st_matrix(piname,pi[av,]');st_matrix(mapname,map)
    }
    if(hasV) {
        st_matrix(vname,V[av,av']);st_matrix(oname,Omega[av,av'])
        st_matrixcolstripe(vname,stripe);st_matrixrowstripe(vname,stripe)
        st_matrixcolstripe(oname,stripe);st_matrixrowstripe(oname,stripe)
    }
    sourceOmega=_cdid_ord_matrix(P,"ord.cells.source_Omega")
    st_local("ord_has_sourceV",_cdid_ord_key(kind=="ordinary_cells" & rows(sourceOmega)>0 & rows(av)>0))
    if(kind=="ordinary_cells" & rows(sourceOmega)>0 & rows(av)>0) {
        st_matrix(sourceo,sourceOmega[av,av']);st_matrix(sourcev,sourceOmega[av,av']/rows(P.cohort))
        st_matrixcolstripe(sourceo,stripe);st_matrixrowstripe(sourceo,stripe)
        st_matrixcolstripe(sourcev,stripe);st_matrixrowstripe(sourcev,stripe)
    }
    if(kind=="ordinary_cells") {C=_cdid_ord_number(P,"ord.C");clustered=P.cluster_kind!="none";}
    else {C=_cdid_ord_number(P,"ord.report.C");clustered=_cdid_ord_number(P,"ord.report.clustered");}
    st_local("ord_N",_cdid_ord_key(rows(P.cohort)));st_local("ord_C",_cdid_ord_key(C));st_local("ord_clustered",_cdid_ord_key(clustered))
    st_local("ord_K",_cdid_ord_key(rows(point)));st_local("ord_sourceK",_cdid_ord_key(rows(P.point)))
    st_local("ord_npost",_cdid_ord_key(rows(av)));st_local("ord_hasV",_cdid_ord_key(hasV));st_local("ord_point_only",_cdid_ord_key(pointonly))
    st_local("ord_kind",kind);st_local("ord_warning",_cdid_ord_text(P,"ord.report.warning"))
    st_local("ord_covariance_origin",_cdid_ord_text(P,"ord.report.covariance_origin"))
    st_local("ord_wald_status",_cdid_ord_text(P,"ord.cells.wald_status"));st_local("ord_wald_reason",_cdid_ord_text(P,"ord.cells.wald_reason"))
    st_local("ord_wald_df",strofreal(_cdid_ord_number(P,"ord.cells.wald_df"),"%24.17g"))
    st_local("ord_wald_statistic",strofreal(_cdid_ord_number(P,"ord.cells.wald_statistic"),"%24.17g"))
    st_local("ord_wald_p",strofreal(_cdid_ord_number(P,"ord.cells.wald_p"),"%24.17g"))
    st_local("ord_wald_p_unrounded",strofreal(_cdid_ord_number(P,"ord.cells.wald_p_unrounded"),"%24.17g"))
    st_local("ord_alpha",strofreal(_cdid_ord_number(P,"ord.report.alpha"),"%24.17g"))
    st_local("ord_critical",strofreal(_cdid_ord_number(P,"ord.report.critical"),"%24.17g"))
    st_local("ord_normal",strofreal(_cdid_ord_number(P,"ord.report.normal"),"%24.17g"))
    st_local("ord_bstrap",_cdid_ord_key(_cdid_ord_number(P,"ord.report.bstrap")))
    st_local("ord_cband",_cdid_ord_key(_cdid_ord_number(P,"ord.report.cband")))
    st_local("ord_B",_cdid_ord_key(_cdid_ord_number(P,"ord.report.B")))
    return(cdid_state_put(P,limits))
}
void cdid_ord_save(string scalar token,string scalar path,string scalar limits_name)
{
    pointer() scalar ptr
    struct cdid_state_limits scalar limits
    ptr=findexternal(limits_name);if(ptr==NULL) _cdid_ord_fail("Explicit state limits disappeared before save.",111)
    limits=*ptr;cdid_state_save(token,path,limits)
}
end
