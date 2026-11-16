*! Native retained-design GMM algebra for cdid; Stata 14.2 or later.
*! Reference: R cdid 1.1.0.0, commit .
version 14.2
mata:

struct cdid_gmm_result {
    real colvector coefficients, targets, residual
    real matrix coefficient_loadings, target_loadings
    real matrix coefficient_covariance, covariance
    real scalar M_full, L_retained, p_parameters, n_targets
    real scalar workspace_entries_bound, workspace_entries_limit
    real scalar lambda, covariance_scale, trace_B, rho, d_min, d_max
    real scalar covariance_rank, zero_covariance, handler_invoked
    real scalar eigen_residual, eigen_bound, eigen_refined
    real scalar factor_residual, whitening_residual, qr_residual
    real scalar rcond_factor_1, rcond_qr_1, sigma_relative
    real scalar stationarity_residual, loading_residual, criterion
    string scalar weighting_requested, weighting, trigger_reason
    string scalar covariance_handler
    real scalar source_n_units, sample_rank_upper, normalization_factor, normalization_count
    real matrix unit_gram_covariance
    real rowvector unit_gram_mean
    string scalar covariance_provenance
    real matrix retained_if
}
/* Created by complete replay, never by a user PSD flag. */
struct cdid_gmm_unit_gram {
    real scalar n_units, rows_consumed, normalization_factor, sample_rank_upper, normalization_count
    real matrix covariance
    real rowvector mean
    real matrix retained_if
}

void _cdid_gmm_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}

real scalar _cdid_gmm_norm(real matrix X)
{
    real scalar a
    if (rows(X)==0 | cols(X)==0) return(0)
    if (missing(X)) return(.)
    a = max(abs(X))
    if (a==0) return(0)
    return(a*sqrt(sum((X:/a):^2)))
}

real scalar _cdid_gmm_ratio(real matrix X, real scalar scale)
{
    if (rows(X)==0 | cols(X)==0) return(0)
    if (scale>0 & scale<.) return(_cdid_gmm_norm(X)/scale)
    if (scale==0 & !missing(X) & all(X:==0)) return(0)
    return(.)
}

/* No inverse is stored: obtain the inverse 1-norm one RHS at a time.
   This is a computed norm diagnostic, whereas R uses LAPACK dtrcon's estimate. */
real scalar _cdid_gmm_rcond1(real matrix T, real scalar upper)
{
    real scalar j, n, z, maxcol, rank_T, norm_T
    real colvector b, x
    real matrix Xinv
    n = rows(T)
    if (n==0) return(0)
    norm_T = max(colsum(abs(T)))
    if (norm_T==0 | norm_T>=.) return(0)
    if (n<=1024) {
        if (upper) Xinv = solveupper(T,I(n),rank_T=.,0)
        else Xinv = solvelower(T,I(n),rank_T=.,0)
        if (rank_T!=n | missing(Xinv)) return(0)
        maxcol = max(colsum(abs(Xinv)))
        if (maxcol==0 | maxcol>=.) return(0)
        return((1/norm_T)/maxcol)
    }
    maxcol = 0
    for (j=1; j<=n; j++) {
        b = J(n,1,0)
        b[j] = 1
        if (upper) x = solveupper(T,b,rank_T=.,0)
        else x = solvelower(T,b,rank_T=.,0)
        if (rank_T!=n | missing(x)) return(0)
        z = sum(abs(x))
        if (z>=.) return(0)
        maxcol = max((maxcol,z))
    }
    if (maxcol==0) return(0)
    return((1/norm_T)/maxcol)
}

real scalar _cdid_gmm_cells(real scalar a, real scalar b,
    real scalar cap, real scalar total)
{
    real scalar z
    if (a>0 & b>(2^53-1)/a)
        _cdid_gmm_fail("Workspace dimension product exceeds exact integer range.",908)
    z = a*b
    if (z>cap-total)
        _cdid_gmm_fail("Required GMM workspace exceeds max_matrix_entries.",908)
    return(total+z)
}

real scalar _cdid_gmm_workspace(real scalar L, real scalar P, real scalar Jt,
    real scalar cap, real scalar total)
{
    total=_cdid_gmm_cells(30*L,L,cap,total)
    total=_cdid_gmm_cells(6*L,P,cap,total)
    total=_cdid_gmm_cells(8*L,P+Jt,cap,total)
    total=_cdid_gmm_cells(6*P,P+Jt,cap,total)
    total=_cdid_gmm_cells(4*P,P,cap,total)
    total=_cdid_gmm_cells(4*Jt,Jt,cap,total)
    total=_cdid_gmm_cells(56,L,cap,total)
    total=_cdid_gmm_cells(16,2*P+Jt,cap,total)
    return(total)
}

real scalar _cdid_gmm_identity_workspace(real scalar L, real scalar P, real scalar Jt,
    real scalar cap, real scalar total)
{
    total=_cdid_gmm_cells(8*L,P,cap,total)
    total=_cdid_gmm_cells(4*L,P+Jt,cap,total)
    total=_cdid_gmm_cells(6*P,P+Jt,cap,total)
    total=_cdid_gmm_cells(4*P,P,cap,total)
    total=_cdid_gmm_cells(4*Jt,Jt,cap,total)
    total=_cdid_gmm_cells(16,L+P+Jt,cap,total)
    return(total)
}
/* The QR equations, shared by all whitening branches and the
   covariance-free identity entry. No covariance/identity-metric allocation. */
struct cdid_gmm_result scalar _cdid_gmm_qr_map(real matrix X, real matrix y, real matrix H)
{
    struct cdid_gmm_result scalar out
    real scalar L, P, Jt, rank_T, denom
    real matrix Qh, Rq, gramA, Hboth, Z, Lambda, ep
    real rowvector tau, piv, da
    real colvector qy, beta_p, e
    L=rows(X); P=cols(X); Jt=rows(H)
    /* Full-rank domain: never truncate, ridge, drop, or select parameters. */
    hqrdp(X,Qh=.,tau=.,Rq=.,piv=.)
    if (missing(Rq) | any(diagonal(Rq):==0))
        _cdid_gmm_fail("Retained design rank/QR failure; covariance ridge cannot identify parameters.",459)
    out.rcond_qr_1=_cdid_gmm_rcond1(Rq,1)
    gramA=quadcross(X,X)
    if (missing(gramA)) _cdid_gmm_fail("Whitened design diagnostic Gram overflowed.",459)
    da=symeigenvalues(gramA)
    out.sigma_relative=sqrt(max((0,min(da)))/max(da))
    if (missing(da) | min(da)<=0 | out.sigma_relative>=. | out.sigma_relative<1e-10 |
        out.rcond_qr_1<1e-10)
        _cdid_gmm_fail("Retained design rank/conditioning refused; no parameter ridge or rank repair.",459)
    qy=hqrdmultq1t(Qh,tau,y)
    beta_p=solveupper(Rq,qy,rank_T=.,0)
    out.coefficients=beta_p[invorder(piv)]
    out.targets=H*out.coefficients
    Hboth=I(P)\H
    /* Reorder target coefficient rows, not target columns. */
    Z=solvelower(Rq',Hboth[,piv]',rank_T=.,0)
    Lambda=hqrdmultq(Qh,tau,Z\J(L-P,P+Jt,0),0)
    out.coefficient_loadings=Lambda[,1..P]'
    out.target_loadings=(Jt ? Lambda[,P+1..P+Jt]' : J(0,L,0))
    out.coefficient_covariance=J(0,0,0); out.covariance=J(0,0,0)
    out.retained_if=J(0,0,.)
    out.residual=y-X*out.coefficients
    ep=X[,piv]-hqrdq1(Qh,tau)*Rq
    out.qr_residual=_cdid_gmm_ratio(ep,_cdid_gmm_norm(X))
    e=y-X*out.coefficients
    out.criterion=quadcross(e,e)
    denom=_cdid_gmm_norm(X)*(_cdid_gmm_norm(y)+_cdid_gmm_norm(X)*_cdid_gmm_norm(out.coefficients))
    out.stationarity_residual=_cdid_gmm_ratio(quadcross(X,e),denom)
    out.loading_residual=_cdid_gmm_ratio(out.target_loadings*X-H,_cdid_gmm_norm(H))
    if(missing(out.coefficients) | missing(out.targets) | missing(Lambda) |
        out.qr_residual>=. | out.criterion>=. | out.stationarity_residual>=. | out.loading_residual>=.)
        _cdid_gmm_fail("Nonfinite GMM map or diagnostics.",504)
    return(out)
}

/* r: L x 1 primitive moments; A: L x P retained full-column-rank design.
   Vr: L x L exactly symmetric input covariance, in caller-chosen units.
   G: L x L exactly symmetric SPD full-catalogue metric E'E; M: full count.
   H: J x P targets, possibly J=0. No catalogue support or gauge selection here.
   weighting: identity, optimal, or two_step (alias for optimal).
   Tikhonov only; source defaults supplied explicitly are 1e-10 and 1e-4.
   Optional cap is a caller-declared matrix-entry workspace, not a work limit.
   Coefficient and target loadings multiply primitive moments/IF columns.
   Returned covariance = loadings * Vr * loadings'; no N scaling is inferred. */
struct cdid_gmm_result scalar _cdid_gmm_fit_core(real matrix r, real matrix A,
    real matrix Vr, real matrix G, real scalar M, real matrix H,
    string scalar weighting, real scalar condition_trigger,
    real scalar ridge_relative, real scalar max_matrix_entries,
    struct cdid_gmm_unit_gram scalar source)
{
    struct cdid_gmm_result scalar out, qr_map
    real scalar L, P, Jt, i, j, total, cap, rank_T, trigger, zero
    real scalar factor_floor, nB, rc, denom, internal_gram
    real matrix C, B, U, Qh, Rq, F, X, y, D, z, gramA
    real matrix Z, Lambda, Hboth, ep, U2, Bfit
    real rowvector d, piv, tau, da, order_d
    real colvector qy, beta_p, e
    L = rows(A); P = cols(A); Jt = rows(H)
    if (L<1 | P<1 | P>L | rows(r)!=L | cols(r)!=1 |
        rows(Vr)!=L | cols(Vr)!=L | rows(G)!=L | cols(G)!=L | cols(H)!=P)
        _cdid_gmm_fail("Nonconformable primitive/design/covariance/target dimensions.",3200)
    if (missing(r) | missing(A) | missing(Vr) | missing(G) | missing(H))
        _cdid_gmm_fail("All matrix inputs must be finite and nonmissing.",498)
    if (M>=. | M<L | M!=floor(M) | M>2^31-1)
        _cdid_gmm_fail("M_full must be an integer between L and 2^31-1.",498)
    if (condition_trigger>=. | condition_trigger<=0 | condition_trigger>=1 |
        ridge_relative>=. | ridge_relative<=0)
        _cdid_gmm_fail("Positive finite tuning required; condition_trigger must be below 1.",498)
    if (weighting!="identity" & weighting!="optimal" & weighting!="two_step")
        _cdid_gmm_fail("weighting must be identity, optimal, or two_step.",498)
    if (!all(Vr:==Vr') | !all(G:==G'))
        _cdid_gmm_fail("Vr and G must be exactly symmetric.",498)
    cap = 2^53-1
    if (max_matrix_entries>=. | max_matrix_entries<1 |
        max_matrix_entries!=floor(max_matrix_entries) | max_matrix_entries>cap)
        _cdid_gmm_fail("max_matrix_entries must be a positive exact finite integer.",498)
    cap = max_matrix_entries
    internal_gram=source.n_units>0 & source.n_units<.
    out.covariance_provenance="arbitrary_covariance"
    out.source_n_units=.; out.sample_rank_upper=.; out.normalization_factor=.; out.normalization_count=.
    out.unit_gram_covariance=J(0,0,0); out.unit_gram_mean=J(1,0,0); out.retained_if=J(0,0,.)
    if (internal_gram) {
        if (source.rows_consumed!=source.n_units | source.n_units!=floor(source.n_units) |
            source.normalization_count<1 | source.normalization_count>source.n_units |
            source.normalization_count!=floor(source.normalization_count) |
            source.normalization_factor!=1/source.normalization_count | source.sample_rank_upper!=source.n_units-1 |
            rows(source.covariance)!=L | cols(source.covariance)!=L)
            _cdid_gmm_fail("Internal complete-unit Gram metadata is inconsistent.",506)
        if (!all(Vr:==source.covariance))
            _cdid_gmm_fail("Internal complete-unit Gram does not match weight covariance.",506)
        out.covariance_provenance="internally_constructed_complete_unit_if_gram"
        out.source_n_units=source.n_units; out.sample_rank_upper=source.sample_rank_upper
        out.normalization_factor=source.normalization_factor; out.normalization_count=source.normalization_count
        out.unit_gram_covariance=source.covariance; out.unit_gram_mean=source.mean; out.retained_if=source.retained_if
    }
    /* Same handler allocation formula; also composed by engine. */
    total=_cdid_gmm_workspace(L,P,Jt,cap,0)
    out.workspace_entries_bound=total; out.workspace_entries_limit=cap
    C = cholesky(G)'
    if (missing(C) | any(diagonal(C):<=0))
        _cdid_gmm_fail("Full-catalogue metric G is not positive definite.",506)
    out.weighting_requested=weighting
    out.weighting=(weighting=="two_step" ? "optimal" : weighting)
    out.covariance_handler="tikhonov"
    out.M_full=M; out.L_retained=L; out.p_parameters=P; out.n_targets=Jt
    out.lambda=0; out.handler_invoked=0; out.zero_covariance=0
    out.trace_B=.; out.covariance_scale=.; out.rho=.; out.d_min=.; out.d_max=.
    out.covariance_rank=.; out.eigen_residual=.; out.eigen_bound=.; out.eigen_refined=0
    out.rcond_factor_1=.; out.factor_residual=0; out.whitening_residual=0
    out.trigger_reason="not_applicable_identity"
    trigger=0; zero=0
    D=C*A; z=C*r
    if (missing(D) | missing(z))
        _cdid_gmm_fail("Full-catalogue transformation overflowed.",504)
    if (out.weighting=="optimal") {
        B=C*Vr*C'
        if (missing(B)) _cdid_gmm_fail("Transformed covariance overflowed.",506)
        /* R copies the internally computed upper triangle to the lower. */
        for (j=1;j<=L;j++) for (i=1;i<j;i++) B[j,i]=B[i,j]
        zero=all(Vr:==0) & all(B:==0)
        symeigensystem(B,U=.,d=.)
        if (missing(U) | missing(d))
            _cdid_gmm_fail("Transformed covariance eigendecomposition failed.",506)
        out.d_min=min(d); out.d_max=max(d)
        nB=_cdid_gmm_norm(B)
        out.eigen_bound=100*L*epsilon(1)*nB
        out.eigen_residual=(zero ? 0 : _cdid_gmm_norm(B-U*diag(d)*U'))
        /* The one source-authorized eigensystem refinement; B is unchanged. */
        if (!zero & out.eigen_residual>out.eigen_bound & out.eigen_residual<. & out.eigen_bound<.) {
            hqrd(U,Qh=.,tau=.,Rq=.)
            U2=hqrdq1(Qh,tau)
            da=colsum(U2:*(B*U2))
            order_d=order(da',-1)'
            da=da[order_d]; U2=U2[,order_d]
            Bfit=U2*diag(da)*U2'
            rc=_cdid_gmm_norm(B-Bfit)
            if (!missing(da) & rc<=out.eigen_bound) {
                U=U2; d=da; out.eigen_residual=rc; out.eigen_refined=1
                out.d_min=min(d); out.d_max=max(d)
            }
        }
        if (!zero) {
            if (out.eigen_bound>=. | out.eigen_residual>=. | out.eigen_residual>out.eigen_bound)
                _cdid_gmm_fail("Covariance eigen residual exceeds the source roundoff bound.",506)
            if (out.d_min < -(out.eigen_bound+out.eigen_residual) |
                any(diagonal(B):< -(out.eigen_bound+out.eigen_residual)))
                _cdid_gmm_fail("Weight covariance is materially indefinite.",506)
            if (out.d_min<0 & !internal_gram)
                _cdid_gmm_fail("Negative covariance eigenvalue requires authenticated internal IF-Gram provenance.",506)
            out.trace_B=trace(B); out.covariance_scale=out.trace_B/M
            if (out.trace_B<=0 | out.trace_B>=. | out.covariance_scale<=0 |
                out.covariance_scale>=. | out.d_max<=0)
                _cdid_gmm_fail("Nonzero covariance has nonpositive or nonfinite scale.",506)
        }
        else { out.trace_B=0; out.covariance_scale=0; }
        out.rho=(zero | out.d_max<=0 ? 0 : max((out.d_min,0))/out.d_max)
        out.covariance_rank=sum(d:>(out.eigen_bound+out.eigen_residual))
        trigger=zero | out.d_min<=0 | out.rho<condition_trigger
        if (zero) out.trigger_reason="zero_covariance"
        else if (out.d_min<0) out.trigger_reason="roundoff_negative_psd"
        else if (out.d_min==0) out.trigger_reason="exact_zero_eigenvalue"
        else if (out.rho<condition_trigger) out.trigger_reason="spectral_condition_below_trigger"
        else out.trigger_reason="well_conditioned_pd"
        if (internal_gram & !trigger) {
            if (L>source.sample_rank_upper) {
                trigger=1; out.trigger_reason="known_unit_covariance_rank_bound"
            }
        }
        if (!trigger) {
            /* Match R's ordinary branch, whose primitive whitening cancels G. */
            F=cholesky(Vr)
            if (missing(F) | any(diagonal(F):<=0)) {
                trigger=1
                out.trigger_reason="cholesky_failure_on_valid_psd"
            }
        }
        if (trigger) {
            out.handler_invoked=1; out.zero_covariance=zero
            if (!zero) {
                out.lambda=ridge_relative*out.covariance_scale
                if (out.lambda<=0 | out.lambda>=.)
                    _cdid_gmm_fail("Covariance ridge lambda is nonpositive or nonfinite.",506)
                B=B+out.lambda*I(L)
                F=cholesky(B)
                if (missing(F) | any(diagonal(F):<=0))
                    _cdid_gmm_fail("Regularized covariance Cholesky failed; no fallback.",506)
                X=solvelower(F,D,rank_T=.,0)
                y=solvelower(F,z,rank_T=.,0)
                out.factor_residual=_cdid_gmm_ratio(F*F'-B,_cdid_gmm_norm(B))
                out.whitening_residual=_cdid_gmm_ratio(F*X-D,_cdid_gmm_norm(D))
            }
            else {
                X=D; y=z; F=C'
                out.factor_residual=_cdid_gmm_ratio(C'*C-G,_cdid_gmm_norm(G))
            }
        }
        else {
            X=solvelower(F,A,rank_T=.,0)
            y=solvelower(F,r,rank_T=.,0)
            out.factor_residual=_cdid_gmm_ratio(F*F'-Vr,_cdid_gmm_norm(Vr))
            out.whitening_residual=_cdid_gmm_ratio(F*X-A,_cdid_gmm_norm(A))
        }
        out.rcond_factor_1=_cdid_gmm_rcond1(F',1)
        factor_floor=min((1e-10,sqrt(condition_trigger)))
        if (!zero & (out.rcond_factor_1<factor_floor | out.factor_residual>=. |
            out.factor_residual>1e-10 | out.whitening_residual>=. | out.whitening_residual>1e-10))
            _cdid_gmm_fail("Covariance factor conditioning or factor/whitening residual refused.",506)
    }
    else { X=D; y=z; }
    if (missing(X) | missing(y))
        _cdid_gmm_fail("Whitened system contains nonfinite values.",504)
    qr_map=_cdid_gmm_qr_map(X,y,H)
    out.coefficients=qr_map.coefficients; out.targets=qr_map.targets
    Lambda=(qr_map.coefficient_loadings\qr_map.target_loadings)'
    if (out.weighting=="identity" | zero) Lambda=C'*Lambda
    else if (trigger) Lambda=C'*solveupper(F',Lambda,rank_T=.,0)
    else Lambda=solveupper(F',Lambda,rank_T=.,0)
    out.coefficient_loadings=Lambda[,1..P]'
    out.target_loadings=(Jt ? Lambda[,P+1..P+Jt]' : J(0,L,0))
    out.coefficient_covariance=out.coefficient_loadings*Vr*out.coefficient_loadings'
    out.covariance=out.target_loadings*Vr*out.target_loadings'
    out.residual=r-A*out.coefficients
    out.rcond_qr_1=qr_map.rcond_qr_1; out.sigma_relative=qr_map.sigma_relative
    out.qr_residual=qr_map.qr_residual; out.criterion=qr_map.criterion
    out.stationarity_residual=qr_map.stationarity_residual
    out.loading_residual=_cdid_gmm_ratio(out.target_loadings*A-H,_cdid_gmm_norm(H))
    if (missing(out.coefficients) | missing(out.targets) | missing(Lambda) |
        missing(out.coefficient_covariance) | missing(out.covariance) |
        out.qr_residual>=. | out.criterion>=. | out.stationarity_residual>=. |
        out.loading_residual>=.)
        _cdid_gmm_fail("Nonfinite GMM solution, covariance, or diagnostics.",504)
    return(out)
}

/* Full-row identity map. Covariances remain explicitly absent here and are
   supplied by the engine's complete-F projected IF stream, never fake zeros. */
struct cdid_gmm_result scalar cdid_gmm_identity_map(real matrix r,real matrix A,real matrix H,real scalar cap)
{
    struct cdid_gmm_result scalar out
    real scalar L, P, Jt, entries
    rmexternal("__cdid_gmm_base_ret")
    L=rows(A); P=cols(A); Jt=rows(H)
    if(L<1 | P<1 | P>L | rows(r)!=L | cols(r)!=1 | cols(H)!=P)
        _cdid_gmm_fail("Nonconformable full-row identity map.",3200)
    if(missing(r) | missing(A) | missing(H))_cdid_gmm_fail("Identity map inputs must be finite.",498)
    if(cap<1 | cap>=. | cap!=floor(cap) | cap>2^53-1)_cdid_gmm_fail("Invalid identity map workspace cap.",498)
    entries=_cdid_gmm_identity_workspace(L,P,Jt,cap,0)
    out=_cdid_gmm_qr_map(A,r,H)
    out.M_full=L;out.L_retained=L;out.p_parameters=P;out.n_targets=Jt
    out.workspace_entries_bound=entries;out.workspace_entries_limit=cap
    out.weighting_requested="identity";out.weighting="identity";out.covariance_handler="tikhonov"
    out.trigger_reason="not_applicable_identity";out.lambda=0;out.handler_invoked=0;out.zero_covariance=0
    out.covariance_scale=.;out.trace_B=.;out.rho=.;out.d_min=.;out.d_max=.;out.covariance_rank=.
    out.eigen_residual=.;out.eigen_bound=.;out.eigen_refined=0
    out.factor_residual=0;out.whitening_residual=0;out.rcond_factor_1=.
    out.source_n_units=.;out.sample_rank_upper=.;out.normalization_factor=.;out.normalization_count=.
    out.unit_gram_covariance=J(0,0,0);out.unit_gram_mean=J(1,0,0);out.retained_if=J(0,0,.)
    out.covariance_provenance="identity_map_covariance_from_projected_IF"
    return(out)
}

/* The arbitrary-covariance API remains strict. No provenance override. */
struct cdid_gmm_result scalar cdid_gmm_fit(real matrix r, real matrix A,
    real matrix Vr, real matrix G, real scalar M, real matrix H,
    string scalar weighting, real scalar condition_trigger,
    real scalar ridge_relative, | real scalar max_matrix_entries)
{
    struct cdid_gmm_unit_gram scalar source
    real scalar cap
    cap=(args()==10 ? max_matrix_entries : 2^53-1)
    source.n_units=0; source.covariance=J(0,0,0); source.retained_if=J(0,0,.)
    return(_cdid_gmm_fit_core(r,A,Vr,G,M,H,weighting,condition_trigger,ridge_relative,cap,source))
}
/* Complete unit-IF Gram is FORMED here, with fixed complete-row coverage.
   Provider(context,lo,hi) is owned by the calling engine; it returns raw IF
   rows lo..hi, with no covariance/PSD boolean from callers. The stream
   centers over n replay rows, normalizes by its explicit count, and carries rank n-1.
   No N-by-L IF matrix, security-seal or interval-certificate framework. */
struct cdid_gmm_unit_gram scalar _cdid_gmm_unit_gram_norm(real scalar n, real scalar L,
    real scalar normalization_count, real scalar chunk, pointer() scalar provider, pointer() scalar context, real scalar cap)
{
    struct cdid_gmm_unit_gram scalar source
    real scalar lo, hi, entries, i, j, batch
    real matrix Z, Zc, Z_full
    real rowvector sum_phi
    if(n<1 | n>=. | n!=floor(n) | n>2^31-1 | L<1 | L>=. | L!=floor(L) |
        normalization_count<1 | normalization_count>=. | normalization_count>n |
        normalization_count!=floor(normalization_count) |
        chunk<1 | chunk>=. | chunk!=floor(chunk) | provider==NULL | context==NULL)
        _cdid_gmm_fail("Invalid internally owned complete-unit IF replay.",498)
    if(cap<1 | cap>=. | cap!=floor(cap) | cap>2^53-1)
        _cdid_gmm_fail("Invalid complete-unit Gram workspace cap.",498)
    batch=min((n,chunk)); entries=0
    entries=_cdid_gmm_cells(2*batch,L,cap,entries)
    entries=_cdid_gmm_cells(2*L,L,cap,entries)
    entries=_cdid_gmm_cells(4,L,cap,entries)
    sum_phi=J(1,L,0)
    if (entries + n * L <= cap) {
        Z_full = J(n, L, .)
        for(lo=1;lo<=n;lo=lo+batch) {
            hi=min((n,lo+batch-1)); Z=(*provider)(context,lo,hi)
            if(rows(Z)!=hi-lo+1 | cols(Z)!=L | missing(Z))
                _cdid_gmm_fail("Owned IF replay does not cover its complete finite row batch.",504)
            Z_full[|lo,1 \ hi,L|] = Z
            sum_phi=sum_phi+colsum(Z)
        }
        source.n_units=n; source.rows_consumed=n; source.normalization_count=normalization_count
        source.normalization_factor=1/normalization_count
        source.sample_rank_upper=n-1; source.mean=sum_phi/n
        Zc=Z_full:-source.mean
        // Double-precision BLAS Gram (differs from quadcross by rounding only, about 1e-16 relative).
        source.covariance=cross(Zc,Zc)/normalization_count
        source.covariance=makesymmetric(source.covariance')
        if(missing(source.covariance) | missing(source.mean))
            _cdid_gmm_fail("Nonfinite internally constructed unit-IF Gram.",504)
        source.retained_if=Z_full
        return(source)
    }
    source.retained_if=J(0,0,.)
    for(lo=1;lo<=n;lo=lo+batch) {
        hi=min((n,lo+batch-1)); Z=(*provider)(context,lo,hi)
        if(rows(Z)!=hi-lo+1 | cols(Z)!=L | missing(Z))
            _cdid_gmm_fail("Owned IF replay does not cover its complete finite row batch.",504)
        sum_phi=sum_phi+colsum(Z)
    }
    source.n_units=n; source.rows_consumed=0; source.normalization_count=normalization_count
    source.normalization_factor=1/normalization_count
    source.sample_rank_upper=n-1; source.mean=sum_phi/n
    source.covariance=J(L,L,0)
    for(lo=1;lo<=n;lo=lo+batch) {
        hi=min((n,lo+batch-1)); Z=(*provider)(context,lo,hi)
        if(rows(Z)!=hi-lo+1 | cols(Z)!=L | missing(Z))
            _cdid_gmm_fail("Owned IF replay does not cover its complete finite Gram batch.",504)
        Zc=Z:-source.mean
        source.covariance=source.covariance+cross(Zc,Zc)
        source.rows_consumed=source.rows_consumed+rows(Z)
    }
    source.covariance=source.covariance/normalization_count
    source.covariance=makesymmetric(source.covariance')
    if(missing(source.covariance) | missing(source.mean))
        _cdid_gmm_fail("Nonfinite internally constructed unit-IF Gram.",504)
    return(source)
}
/* Old entry preserves normalization_count=n exactly. */
struct cdid_gmm_unit_gram scalar _cdid_gmm_unit_gram(real scalar n,real scalar L,
    real scalar chunk,pointer() scalar provider,pointer() scalar context,real scalar cap)
{
    return(_cdid_gmm_unit_gram_norm(n,L,n,chunk,provider,context,cap))
}
/* R replays complete F rows, centers over F, normalizes by n_A.
   The constructor forms the actual Gram; no covariance/provenance override. */
struct cdid_gmm_result scalar cdid_gmm_fit_unit_norm(real matrix r,real matrix A,
    real matrix G,real scalar M,real matrix H,string scalar weighting,
    real scalar trigger,real scalar ridge,real scalar replay_rows,real scalar normalization_count,
    real scalar chunk,pointer() scalar provider,pointer() scalar context,real scalar cap)
{
    struct cdid_gmm_unit_gram scalar source
    source=_cdid_gmm_unit_gram_norm(replay_rows,rows(A),normalization_count,chunk,provider,context,cap)
    return(_cdid_gmm_fit_core(r,A,source.covariance,G,M,H,weighting,trigger,ridge,cap,source))
}
struct cdid_gmm_result scalar cdid_gmm_fit_unit_stream(real matrix r, real matrix A,
    real matrix G, real scalar M, real matrix H, string scalar weighting,
    real scalar condition_trigger, real scalar ridge_relative, real scalar n_units,
    real scalar chunk, pointer() scalar provider, pointer() scalar context, real scalar cap)
{
    struct cdid_gmm_unit_gram scalar source
    source=_cdid_gmm_unit_gram(n_units,rows(A),chunk,provider,context,cap)
    return(_cdid_gmm_fit_core(r,A,source.covariance,G,M,H,weighting,condition_trigger,ridge_relative,cap,source))
}
/* GMM with an explicit within-unit weight covariance Vw. The complete unit-IF Gram of the
   moments is still formed from the owned replay, so the reported covariance and the
   retained IF are those of the estimator, not of the weight. With Lprov and cols, the
   provider returns Lprov columns of which cols are the GMM moments; the retained IF keeps
   all Lprov columns. */
struct cdid_gmm_result scalar cdid_gmm_fit_within_norm(real matrix r, real matrix A,
    real matrix Vw, real matrix G, real scalar M, real matrix H, real scalar trigger,
    real scalar ridge, real scalar replay_rows, real scalar normalization_count,
    real scalar chunk, pointer() scalar provider, pointer() scalar context, real scalar cap,
    | real scalar Lprov, real rowvector cols)
{
    struct cdid_gmm_unit_gram scalar source
    struct cdid_gmm_result scalar out
    real matrix S
    real rowvector mu
    if (args() < 16) {
        Lprov = rows(A)
        cols = (1..rows(A))
    }
    source = _cdid_gmm_unit_gram_norm(replay_rows, Lprov, normalization_count, chunk, provider, context, cap)
    S = source.covariance[cols, cols]
    mu = source.mean[cols]
    source.covariance = Vw
    out = _cdid_gmm_fit_core(r, A, Vw, G, M, H, "optimal", trigger, ridge, cap, source)
    out.weighting_requested = "within"
    out.covariance_provenance = "analytical_within_unit_levels_weight"
    out.unit_gram_covariance = S
    out.unit_gram_mean = mu
    out.coefficient_covariance = out.coefficient_loadings * S * out.coefficient_loadings'
    out.covariance = out.target_loadings * S * out.target_loadings'
    return(out)
}

end
