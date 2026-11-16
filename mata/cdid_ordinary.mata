*! Internal ordinary-mode extraction for cdid, current StataNow target.
*! Backend: unchanged MIT csdid2.0.2, be5d272debb1c98e18c7fba79f51b2b3effd3fc4.
*! Doubly robust repeated cross-sections use the corrected influence function.
version 19.5
mata:
void _cdid_ordinary_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_ordinary_index(real matrix x)
{
    return(vec(selectindex(x)))
}
/* did2.5.1 validate_args runs these checks BEFORE complete-case screening,
   including true RCS with a supplied idname. RCS IDs need not be globally
   unique, because did later creates distinct observation identities. */
void cdid_ordinary_input_guard(string scalar idname, string scalar tname,
    string scalar gname, string scalar requested)
{
    real colvector use, id, tt, gg, ix, o
    real scalar n
    use=st_data(.,requested):!=0; gg=st_data(.,gname); tt=st_data(.,tname)
    if (any(select(gg,use):<0)) _cdid_ordinary_fail("gvar must be nonnegative; zero denotes never treated.",459)
    if (idname=="") return
    id=st_data(.,idname)
    ix=_cdid_ordinary_index(use :& (id:<.) :& (gg:<.))
    n=rows(ix)
    if (n>1) {
        o=order((id[ix],gg[ix]),(1,2)); ix=ix[o]
        if (any((id[ix[|2\n|]]:==id[ix[|1\n-1|]]) :&
            (gg[ix[|2\n|]]:!=gg[ix[|1\n-1|]])))
            _cdid_ordinary_fail("gvar must be constant within each supplied ID, including RCS IDs.",459)
    }
    ix=_cdid_ordinary_index(use :& (id:<.) :& (tt:<.))
    n=rows(ix)
    if (n>1) {
        o=order((id[ix],tt[ix]),(1,2)); ix=ix[o]
        if (any((id[ix[|2\n|]]:==id[ix[|1\n-1|]]) :&
            (tt[ix[|2\n|]]:==tt[ix[|1\n-1|]])))
            _cdid_ordinary_fail("Duplicate supplied (ID,time) rows are not allowed, including missing-outcome rows.",459)
    }
}
/* A common affine translation preserves calendar distances, never-treated 0,
   raw IDs and all outcome/design data. Refuse nonrepresentable translations:
   neither loss of distinct dates nor rounding away anticipation is allowed. */
real scalar cdid_ordinary_calendar(string scalar tname, string scalar gname,
    string scalar requested, string scalar outtime, string scalar outgroup,
    real scalar anticipation)
{
    real colvector tt, gg, use, it, ig, axis, shifted, edge
    real scalar delta, n
    tt=st_data(.,tname); gg=st_data(.,gname); use=st_data(.,requested):!=0
    it=_cdid_ordinary_index(use :& (tt:<.))
    ig=_cdid_ordinary_index(use :& (gg:>0) :& (gg:<.))
    if (!rows(it)) return(0)
    delta=max((0,1-min(tt[it])))
    if (delta==0) return(0)
    if (delta>=.) _cdid_ordinary_fail("Calendar translation is not finite.",459)
    axis=uniqrows(sort((tt[it]\gg[ig]),1)); n=rows(axis)
    shifted=axis:+delta
    if (missing(shifted) | any(shifted:<1) | any((shifted:-delta):!=axis))
        _cdid_ordinary_fail("Calendar translation cannot preserve supplied dates exactly.",459)
    if (n>1) {
        if (any(shifted[|2\n|]:<=shifted[|1\n-1|]) |
            any((shifted[|2\n|]-shifted[|1\n-1|]):!=(axis[|2\n|]-axis[|1\n-1|])))
            _cdid_ordinary_fail("Calendar translation would collide or change date gaps.",459)
    }
    if (rows(ig)) {
        edge=gg[ig]:-anticipation
        if (missing(edge) | missing(edge:+delta) |
            any(((gg[ig]:+delta):-anticipation):!=(edge:+delta)) |
            any(((edge:+delta):-delta):!=edge))
            _cdid_ordinary_fail("Calendar translation would change anticipation boundaries.",459)
    }
    tt[it]=tt[it]:+delta; gg[ig]=gg[ig]:+delta
    st_store(.,outtime,tt); st_store(.,outgroup,gg)
    return(delta)
}
/* Full-grid covariance is embedded from the backend's reported analytical V.
   Available universal reference rows are structural zeros; unavailable rows
   and columns remain missing. IF is NOT recentered or rescaled here. */
void cdid_ordinary_collect(string scalar idname, real scalar rcs, real scalar calendar_shift,
    string scalar sourcevar, string scalar requested, string scalar sample,
    string scalar cellname, string scalar ifname, string scalar vname,
    string scalar asyname, string scalar unitname, string scalar mapname,
    string scalar groupname, string scalar samplename, string scalar excludedname)
{
    real matrix a, inf, u, bv, V, cells, mapping, gp
    real colvector available, structural, keep, aa, kind, rowsample, excluded, src, wanted, used
    real colvector id, o, look, hit
    real scalar N, K, j, lo, hi, mid, k
    a=st_matrix("e(attgt)"); inf=st_matrix("e(inffunc)"); u=st_matrix("e(unit_group)")
    bv=st_matrix("e(V)"); gp=st_matrix("e(group_prob)")
    N=st_numscalar("e(N_units)"); K=rows(a)
    if (st_global("e(version)")!="2.0.2" | K==0 | cols(a)!=10 |
        rows(inf)!=N | cols(inf)!=K | rows(u)!=N | cols(u)<3)
        _cdid_ordinary_fail("Qualified csdid2.0.2 full-storage results are required.",498)
    if (calendar_shift!=0) {
        a[,1]=a[,1]:-calendar_shift
        a[,2]=a[,2]:-calendar_shift; a[,10]=a[,10]:-calendar_shift
        u[,2]=u[,2]:-calendar_shift:*((u[,2]:>0):&(u[,2]:<.))
        if (cols(u)>=4) u[,4]=u[,4]:-calendar_shift
        gp[,1]=gp[,1]:-calendar_shift:*((gp[,1]:>0):&(gp[,1]:<.))
    }
    available=a[,4]:<.; structural=(a[,2]:==a[,10]); aa=_cdid_ordinary_index(available)
    if (rows(aa)) if (missing(inf[,aa])) _cdid_ordinary_fail("An available backend cell has a nonfinite influence function.",498)
    keep=_cdid_ordinary_index(available :& !structural)
    if (rows(bv)!=rows(keep) | cols(bv)!=rows(keep))
        _cdid_ordinary_fail("Backend analytical covariance order/dimension does not match its nonreference cells.",498)
    V=J(K,K,.)
    if (rows(aa)) V[aa,aa]=J(rows(aa),rows(aa),0)
    if (rows(keep)) V[keep,keep]=bv
    kind=2:*(a[,2]:<a[,1]) :+ 3:*(a[,2]:>=a[,1])
    kind=kind:*(1:-structural):+structural
    /* group,time,base,role(1reference/2pre/3post),available,ATT,SE,event,counts */
    cells=(a[,1],a[,2],a[,10],kind,available,a[,4],a[,5],a[,3],a[,6..9])
    src=st_data(.,sourcevar); wanted=st_data(.,requested):!=0; used=st_data(.,sample):!=0
    rowsample=_cdid_ordinary_index(used); excluded=_cdid_ordinary_index(wanted :& !used)
    if (rows(rowsample)!=st_numscalar("e(N)")) _cdid_ordinary_fail("Backend sample marker/count mismatch.",498)
    mapping=J(N,3,.)
    if (rcs) {
        if (any(u[,1]:<1) | any(u[,1]:>rows(src)) | any(u[,1]:!=floor(u[,1])))
            _cdid_ordinary_fail("Invalid backend RCS observation identity.",498)
        mapping[,1]=(1::N); mapping[,2]=src[u[,1]]
        mapping[,3]=idname=="" ? mapping[,2] : st_data(u[,1],idname)
        u[,1]=mapping[,1]
    } else {
        id=st_data(.,idname); mapping[,1]=u[,1]; mapping[,3]=u[,1]
        o=order(u[,1],1); look=u[o,1]
        for (j=1;j<=rows(rowsample);j++) {
            k=rowsample[j]; lo=1; hi=N
            while (lo<=hi) {
                mid=floor((lo+hi)/2)
                if (look[mid]==id[k]) {
                    hit=o[mid]
                    if (mapping[hit,2]>=.) mapping[hit,2]=src[k]
                    break
                }
                if (look[mid]<id[k]) lo=mid+1
                else hi=mid-1
            }
        }
        if (missing(mapping)) _cdid_ordinary_fail("Panel unit/source-row map is incomplete.",498)
    }
    if (cols(gp)!=3 | any(gp[,3]:!=N)) _cdid_ordinary_fail("Backend cohort-share denominator does not match IF units.",498)
    st_matrix(cellname,cells); st_matrix(ifname,inf); st_matrix(vname,V); st_matrix(asyname,N*V)
    st_matrix(unitname,u); st_matrix(mapname,mapping); st_matrix(groupname,gp)
    st_matrix(samplename,src[rowsample])
    if (rows(excluded)) st_matrix(excludedname,src[excluded])
}
end
