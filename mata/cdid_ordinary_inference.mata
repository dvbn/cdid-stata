*! Ordinary raw-IF inference. Independent native implementation from statistical definitions.
*! did2.5.1 mboot/att_gt are reference only. No GPL-2 code reused.
*! Tagged values: . = NA/NaN; .a = +Inf; .b = -Inf.
version 19.5
mata:
struct cdid_ord_summary {
    real scalar n, C, K, clustered
    real matrix sums, V, Omega, source_Omega
    real colvector analytic_se, source_se, source_excluded
    string scalar covariance_method, source_covariance_method
}
struct cdid_ord_boot {
    real scalar n, C, B, K, critical, return_V
    real matrix bres, V
    real colvector retained, se, scale, maxima
    string scalar axis_policy
}
struct cdid_ord_wald {
    real scalar statistic, p, p_unrounded, df
    real colvector declared, used, excluded
    string scalar status, reason
}
struct cdid_ord_inference {
    struct cdid_ord_summary scalar summary
    struct cdid_ord_boot scalar bootstrap
    struct cdid_ord_wald scalar wald
    real colvector se
    real matrix ci, point_ci
    real scalar critical, normal_critical, bstrap, cband, compute_inffunc, calls
    string scalar warning, kind
}
void _cdid_ord_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_ord_ix(real matrix x)
{
    return(vec(selectindex(x)))
}
real scalar _cdid_ord_tiny()
{
    return(10*sqrt(epsilon(1)))
}
real scalar _cdid_ord_sum(real matrix x)
{
    if (any(x:>=.)) return(.)
    return(sum(x))
}
/* R type-1 quantile, including explicit infinity tags and na.rm=TRUE. */
real scalar cdid_ord_quantile(real matrix x, real scalar p)
{
    real colvector keep, tag, o
    real scalar k
    if (p<0 | p>1 | p>=.) _cdid_ord_fail("Quantile probability must be finite in [0,1].",198)
    x=vec(x)
    keep=_cdid_ord_ix((x:<.) :| (x:==.a) :| (x:==.b))
    if (!rows(keep)) return(.)
    x=x[keep]; tag=(x:==.a):-(x:==.b)
    o=order((tag,x),(1,2)); x=x[o]
    k=max((1,ceil(rows(x)*p)))
    return(x[k])
}
/* Explicit partition supplied in source order. No sorting or IF centering. */
struct cdid_ord_summary scalar cdid_ord_if_summary(real matrix F,
    real colvector cluster, real scalar C, real scalar clustered, real scalar bstrap,
    real scalar row_chunk)
{
    struct cdid_ord_summary scalar S
    real scalar i, j, k, lo, hi, n
    real matrix iid
    n=rows(F); k=cols(F)
    if (n<1 | k<1 | C<1 | C>n | C!=floor(C) | C>=. | rows(cluster)!=n |
        missing(cluster) | any(cluster:<1) | any(cluster:>C) | any(cluster:!=floor(cluster)) |
        row_chunk<1 | row_chunk!=floor(row_chunk) | row_chunk>=. |
        (clustered!=0 & clustered!=1) | (bstrap!=0 & bstrap!=1))
        _cdid_ord_fail("Invalid full-unit IF partition or explicit row chunk.",498)
    if (!clustered & (C!=n | any(cluster:!=(1::n))))
        _cdid_ord_fail("Unit inference requires the complete source identity partition.",498)
    if (rows(uniqrows(sort(cluster,1)))!=C) _cdid_ord_fail("Every declared cluster must contain a unit.",498)
    S.n=n; S.C=C; S.K=k; S.clustered=clustered; S.sums=J(C,k,0)
    for(lo=1;lo<=n;lo=lo+row_chunk) {
        hi=min((n,lo+row_chunk-1))
        for(i=lo;i<=hi;i++) S.sums[cluster[i],]=S.sums[cluster[i],]+F[i,]
    }
    S.V=J(k,k,.); iid=J(k,k,.)
    for(i=1;i<=k;i++) for(j=i;j<=k;j++) {
        S.V[i,j]=S.V[j,i]=_cdid_ord_sum(S.sums[,i]:*S.sums[,j])/(n*n)
        iid[i,j]=iid[j,i]=_cdid_ord_sum(F[,i]:*F[,j])/n
    }
    S.Omega=n*S.V
    S.source_Omega=(clustered & !bstrap ? S.Omega : iid)
    S.analytic_se=sqrt(diagonal(S.V)); S.source_se=sqrt(diagonal(S.source_Omega)/n)
    S.source_excluded=_cdid_ord_ix((S.source_se:>=.) :| (S.source_se:<=_cdid_ord_tiny()))
    if(rows(S.source_excluded)) S.source_se[S.source_excluded]=J(rows(S.source_excluded),1,.)
    S.covariance_method=(clustered ? "raw_cluster_sum_Gram/n^2" : "raw_unit_IF_Gram/n^2")
    S.source_covariance_method=(clustered & !bstrap ? "MP.V_analytical=raw_cluster_sum_Gram/n" : "MP.V_analytical=raw_unit_IF_Gram/n")
    return(S)
}
/* Linear draw stage supports bounded row/cluster chunks; no B-by-C requirement. */
void cdid_ord_draw_add(real matrix Z, real matrix sums, real colvector draws,
    real rowvector clusters, real matrix signs, real scalar C)
{
    if (rows(signs)!=rows(draws) | cols(signs)!=cols(clusters) | C!=rows(sums) |
        missing(signs) | any((signs:!=1):&(signs:!=-1)) | missing(draws) | missing(clusters) |
        any(draws:<1) | any(draws:>rows(Z)) | any(draws:!=floor(draws)) |
        any(clusters:<1) | any(clusters:>C) | any(clusters:!=floor(clusters)) |
        cols(Z)!=cols(sums)) _cdid_ord_fail("Invalid explicit sign/cluster chunk.",498)
    Z[draws,]=Z[draws,]+signs*sums[clusters,]/sqrt(C)
}
struct cdid_ord_boot scalar cdid_ord_boot_none(real scalar K)
{
    struct cdid_ord_boot scalar M
    M.n=M.C=M.B=0; M.K=K; M.critical=.; M.return_V=0
    M.bres=M.V=J(0,0,.); M.retained=M.scale=M.maxima=J(0,1,.)
    M.se=J(K,1,.); M.axis_policy="draw-by-target; approved R-F03 correction"
    return(M)
}
/* Reduce actual source bres or native raw draws. original_positions maps columns
   back to the full target grid. One draw remains one row, including K>1. */
struct cdid_ord_boot scalar cdid_ord_bres_finish(real matrix Z,
    real colvector original_positions, real scalar K, real scalar n, real scalar C,
    real scalar alpha, real scalar return_V)
{
    struct cdid_ord_boot scalar M
    real scalar j, i, b, s, q, z, energy, v
    real colvector keep, scales, maxima, finite
    real matrix centered
    b=rows(Z)
    if (b<1 | n<1 | C<1 | C>n | alpha<=0 | alpha>=1 | alpha>=. |
        rows(original_positions)!=cols(Z) | cols(original_positions)!=1 |
        missing(original_positions) | any(original_positions:<1) | any(original_positions:>K) |
        any(original_positions:!=floor(original_positions)) | (return_V!=0 & return_V!=1))
        _cdid_ord_fail("Invalid ordinary bootstrap dimensions, column map or alpha.",498)
    if(rows(uniqrows(sort(original_positions,1)))!=rows(original_positions))
        _cdid_ord_fail("Bootstrap target map must not repeat.",498)
    keep=J(0,1,.)
    for(j=1;j<=cols(Z);j++) {
        energy=_cdid_ord_sum(Z[,j]:^2)
        if (_cdid_ord_sum(Z[,j])<. & energy<. & energy>_cdid_ord_tiny()) keep=keep\j
    }
    M=cdid_ord_boot_none(K); M.n=n; M.C=C; M.B=b; M.return_V=return_V
    M.bres=Z[,keep']; M.retained=original_positions[keep]; M.scale=J(rows(keep),1,.)
    for(j=1;j<=rows(keep);j++) {
        M.scale[j]=(cdid_ord_quantile(M.bres[,j],.75)-cdid_ord_quantile(M.bres[,j],.25))/(invnormal(.75)-invnormal(.25))
        M.se[M.retained[j]]=M.scale[j]*sqrt(C)/n
    }
    M.V=J(0,0,.)
    if(return_V) {
        M.V=J(rows(keep),rows(keep),.)
        if(b>1 & rows(keep)) {
            centered=M.bres:-colsum(M.bres)/b
            M.V=quadcross(centered,centered)/(b-1)
        }
    }
    maxima=J(b,1,.b)
    for(i=1;i<=b;i++) for(j=1;j<=rows(keep);j++) {
        z=M.bres[i,j]; s=M.scale[j]
        if(z>=. | s>=.) continue
        q=(s==0 ? (z==0 ? . : .a) : abs(z/s))
        if(q==.a) maxima[i]=.a
        else if(q<. & maxima[i]!=.a) maxima[i]=(maxima[i]==.b ? q : max((maxima[i],q)))
    }
    M.maxima=maxima
    finite=_cdid_ord_ix(maxima:<.)
    M.critical=cdid_ord_quantile(maxima[finite],1-alpha)
    return(M)
}
real scalar cdid_ord_cell_band(struct cdid_ord_boot scalar M, real scalar alpha)
{
    real scalar i, j, s, q
    real colvector maxima
    maxima=J(M.B,1,.b)
    for(j=1;j<=cols(M.bres);j++) {
        s=M.scale[j]
        if(s>=. | s<=_cdid_ord_tiny()) continue
        for(i=1;i<=M.B;i++) {
            q=abs(M.bres[i,j]/s)
            if(q<.) maxima[i]=(maxima[i]==.b ? q : max((maxima[i],q)))
        }
    }
    return(cdid_ord_quantile(maxima,1-alpha))
}
/* Exact reciprocal one-norm condition, using documented LU tolerance zero.
   Same numerical convention as the MIT ordinary backend. */
real scalar cdid_ord_rcond1(real matrix V)
{
    real matrix inverse
    real scalar a,b
    a=norm(V,1)
    if(a<=0 | a>=.) return(0)
    inverse=luinv(V,0)
    if(missing(inverse)) return(0)
    b=norm(inverse,1)
    if(b<=0 | b>=.) return(0)
    return(1/(a*b))
}
struct cdid_ord_wald scalar cdid_ord_wald_compute(real colvector point,
    real colvector group, real colvector time, struct cdid_ord_summary scalar S,
    real scalar extra_cluster_requested)
{
    struct cdid_ord_wald scalar W
    real colvector keep
    real matrix V
    real scalar j
    W.statistic=W.p=W.p_unrounded=.; W.df=0; W.status="unavailable"; W.reason=""
    W.declared=_cdid_ord_ix(group:>time); W.used=W.excluded=J(0,1,.)
    if(extra_cluster_requested & S.source_covariance_method!="MP.V_analytical=raw_cluster_sum_Gram/n") {
        W.reason="source omits Wald: requested extra clustering absent from source analytical covariance"
        return(W)
    }
    keep=J(rows(W.declared),1,1)
    for(j=1;j<=rows(W.declared);j++) if(any(S.source_excluded:==W.declared[j])) { keep[j]=0; W.excluded=W.excluded\W.declared[j]; }
    W.used=select(W.declared,keep); W.df=rows(W.used)
    if(!W.df) { W.reason="no pre-treatment cells remain after source analytic variance filter"; return(W); }
    V=S.source_Omega[W.used,W.used']
    if(missing(V)) { W.reason="source pre-treatment covariance contains NA"; return(W); }
    if(cdid_ord_rcond1(V)<=epsilon(1)) { W.reason="source pre-treatment covariance is singular"; return(W); }
    W.statistic=S.n*quadcross(point[W.used],lusolve(V,point[W.used],0))
    W.p_unrounded=chi2tail(W.df,W.statistic)
    W.p=round(W.p_unrounded,1e-5); W.status="computed"; W.reason="source ordinary filtered Wald"
    return(W)
}
struct cdid_ord_inference scalar cdid_ord_report_cells(real colvector point,
    real colvector group, real colvector time, real matrix F, real colvector cluster,
    real scalar C, real scalar clustered, real scalar extra_cluster_requested,
    real scalar compute_inffunc, real scalar bstrap, real scalar cband,
    real scalar alpha, real scalar row_chunk, struct cdid_ord_boot scalar M)
{
    struct cdid_ord_inference scalar O
    real colvector valid
    real scalar j, k
    k=rows(point)
    if(rows(group)!=k | rows(time)!=k | cols(F)!=k | (compute_inffunc!=0 & compute_inffunc!=1) |
        (bstrap!=0 & bstrap!=1) | (cband!=0 & cband!=1) | alpha<=0 | alpha>=1 | alpha>=.)
        _cdid_ord_fail("Invalid ordinary cell reporting request.",498)
    O.compute_inffunc=compute_inffunc; O.bstrap=bstrap*compute_inffunc
    O.cband=cband*compute_inffunc; O.calls=0; O.warning=""; O.kind="ordinary_cells"
    O.normal_critical=invnormal(1-alpha/2); O.critical=O.normal_critical
    O.se=J(k,1,.); O.bootstrap=cdid_ord_boot_none(k)
    O.wald.statistic=O.wald.p=O.wald.p_unrounded=.; O.wald.df=0
    O.wald.declared=O.wald.used=O.wald.excluded=J(0,1,.)
    O.wald.status="unavailable"; O.wald.reason="compute_inffunc=FALSE"
    if(compute_inffunc) {
        O.summary=cdid_ord_if_summary(F,cluster,C,clustered,bstrap,row_chunk)
        O.se=O.summary.source_se
        if(bstrap) {
            if(M.K!=k | M.n!=rows(F) | M.C!=C | M.B<1) _cdid_ord_fail("Bootstrap result does not match the complete ordinary cell grid.",498)
            O.bootstrap=M; O.calls=1
            valid=_cdid_ord_ix(O.summary.source_se:<.)
            O.se[valid]=M.se[valid]
            if(cband) O.critical=cdid_ord_cell_band(M,alpha)
        }
        valid=_cdid_ord_ix(O.se:<=_cdid_ord_tiny())
        if(rows(valid)) O.se[valid]=J(rows(valid),1,.)
        O.wald=cdid_ord_wald_compute(point,group,time,O.summary,extra_cluster_requested)
        if(O.critical<. & O.critical>=7) O.warning="source cell simultaneous critical value is very large"
    }
    O.point_ci=(point:-O.normal_critical*O.se,point:+O.normal_critical*O.se)
    O.ci=(point:-O.critical*O.se,point:+O.critical*O.se)
    /* Tagged infinite criticals require an explicit absent interval, not .a arithmetic. */
    if(O.critical>=.) O.ci=J(k,2,.)
    return(O)
}
end
