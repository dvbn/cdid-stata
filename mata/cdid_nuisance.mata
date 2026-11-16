*! Native nuisance primitives for cdid; Stata 14.2 or later.
*! R reference: cdid 1.1.0.0, .
version 14.2
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_nuis(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_logit_result {
    real colvector beta, rawp, clippedp, clip_active
    real matrix row_score, row_dlogweight, information, R
    real rowvector pivot
    real scalar n_rows, n_active, iterations, converged, deviance, score_max, rcond_R
}
struct cdid_wls_result {
    real colvector beta, alpha, block_ids, mass, mean_y, fitted, residual, prediction
    real colvector row_alpha_direct, row_block
    real matrix mean_X, row_score, row_dlogweight_beta, R, alias_map, alias_constants
    real rowvector kept, dropped, pivot
    real scalar n_rows, n_blocks, fullN, intercept, rss, rank_tol, alias_threshold, rcond_R
}
void _cdid_nuisance_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}
/* selectindex() preserves row orientation for a scalar/row input. Force
   a column even for empty scalar predicates before checking row counts. */
real colvector _cdid_nuisance_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
real scalar _cdid_nuisance_max(real matrix A)
{
    if (rows(A)*cols(A)==0) return(0)
    return(max(abs(A)))
}
real scalar _cdid_nuisance_rcond(real matrix R)
{
    real scalar j, z, n, rk, norm_R
    real colvector v, b
    real matrix Xinv
    n=rows(R); z=0
    if (n==0) return(1)
    norm_R=max(colsum(abs(R)))
    if (norm_R==0 | norm_R>=.) return(0)
    if (n<=1024) {
        Xinv=solveupper(R,I(n),rk=.,0)
        if (rk!=n | missing(Xinv)) return(0)
        z=max(colsum(abs(Xinv)))
        if (z==0 | z>=.) return(0)
        return((1/norm_R)/z)
    }
    for (j=1;j<=n;j++) {
        b=J(n,1,0); b[j]=1
        v=solveupper(R,b,rk=.,0)
        if (rk!=n | missing(v)) return(0)
        z=max((z,sum(abs(v))))
    }
    if (z==0 | z>=.) return(0)
    return((1/norm_R)/z)
}
/* Solve (A'A)x=b using the retained QR, never a Gram inverse. */
real matrix _cdid_nuisance_gram(real matrix R, real rowvector pivot, real matrix b)
{
    real matrix x
    real scalar rk
    if (rows(R)==0) return(J(0,cols(b),0))
    x=solvelower(R',b[pivot,],rk=.,0)
    x=solveupper(R,x,rk=.,0)
    if (missing(x) | rk!=rows(R)) _cdid_nuisance_fail("Nonfinite triangular solve.",504)
    return(x[invorder(pivot),])
}
/* Raw, treatment-capped p and strict derivative mask. No selection floor. */
real matrix cdid_logit_values(real colvector eta)
{
    real colvector p, a, c, ix
    p=_cdid_logit_p(eta)
    a=p:<(1-1e-6)
    c=p
    ix=_cdid_nuisance_index(!a)
    if (rows(ix)) c[ix]=J(rows(ix),1,1-1e-6)
    return((p,c,a))
}
/* Raw probability column of cdid_logit_values(), elementwise identical. */
real colvector _cdid_logit_p(real colvector eta)
{
    real colvector p, hi, lo, o
    if (missing(eta)) _cdid_nuisance_fail("Nonfinite logit predictor.",498)
    if (!any(eta:>40)) {
        o=exp(eta)
        return(o:/(1:+o))
    }
    p=J(rows(eta),1,0)
    hi=_cdid_nuisance_index(eta:>40)
    lo=_cdid_nuisance_index(eta:<=40)
    if (rows(hi)) p[hi]=1:/(1:+exp(-eta[hi]))
    if (rows(lo)) { o=exp(eta[lo]); p[lo]=o:/(1:+o); }
    return(p)
}
real scalar _cdid_logit_deviance(real colvector D, real colvector p, real colvector w)
{
    return(-2*sum(w:*(D:*ln(p)+(1:-D):*ln(1:-p))))
}
/* X is the complete explicit model matrix; add a unit column for an intercept.
   Zero weights are inactive for rank, separation and boundary checks.
   row_dlogweight = [sum w p(1-p) XX']^-1 w X(D-p), with POSITIVE sign.
   Its rows sum to the logweight derivative; multiply unit sums by full-frame
   N to obtain that frame's coefficient IF. No estimated-weight channel is hidden. */
struct cdid_logit_result scalar cdid_logit_fit(real matrix X, real colvector D,
    real colvector w, real scalar eps, real scalar maxit, | real scalar skip_row_if)
{
    struct cdid_logit_result scalar out
    real matrix Xa, Q, R, A, vals, information_chol, SXa, H_step, L_step
    real colvector da, wa, ix, p, eta, z, sw, b, bn, qy, tau, w_irls, z_raw, rhs_step, rhs_ref, d_bn, h_diag
    real rowvector pivot, cmin, cmax, csum, smin, smax
    real scalar n, na, k, j, h, it, rk, dev, olddev, boundary, converged, qtol, no_row_if, chol_ok
    no_row_if=(args()>=6 & skip_row_if==1)
    n=rows(X); k=cols(X)
    if (n==0 | k==0 | rows(D)!=n | cols(D)!=1 | rows(w)!=n | cols(w)!=1 |
        missing(X) | missing(D) | missing(w) | any(w:<0) | any((D:!=0):&(D:!=1)) |
        eps<=0 | eps>=. | maxit<1 | maxit!=floor(maxit) | maxit>=.)
        _cdid_nuisance_fail("Require finite design, binary response, nonnegative weights and explicit GLM controls.",498)
    ix=_cdid_nuisance_index(w:>0); na=rows(ix); out.n_active=na; out.n_rows=n
    if (na<k) _cdid_nuisance_fail("Active logit rows smaller than model dimension.",459)
    Xa=(na==n ? X : X[ix,]); da=(na==n ? D : D[ix]); wa=(na==n ? w : w[ix])
    cmin=colmin(Xa); cmax=colmax(Xa); csum=colsum(Xa)
    SXa=Xa:*(2:*da:-1); smin=colmin(SXa); smax=colmax(SXa)
    for (j=1;j<=k;j++) {
        if (cmin[j]==0 & cmax[j]==0) _cdid_nuisance_fail("Logit zero predictor; model unchanged.",459)
        for (h=1;h<j;h++) {
            if (csum[j]==csum[h] & cmin[j]==cmin[h] & cmax[j]==cmax[h] & Xa[1,j]==Xa[1,h] & Xa[na,j]==Xa[na,h]) {
                if (all(Xa[,j]:==Xa[,h]))
                    _cdid_nuisance_fail("A covariate in the propensity-score model duplicates another column or the intercept (for example, a constant covariate); remove it.",459)
            }
        }
        if ((smin[j]>=0 & smax[j]>0) | (smax[j]<=0 & smin[j]<0))
            _cdid_nuisance_fail("No finite logit root: signed-coordinate separation witness.",459)
    }
    b=J(k,1,0); eta=J(na,1,0); p=J(na,1,.5)
    olddev=_cdid_logit_deviance(da,p,wa); converged=0; boundary=0
    qtol=min((1e-7,eps/1000)); pivot=1..k
    for (it=1;it<=maxit;it++) {
        w_irls=wa:*p:*(1:-p)
        z_raw=eta+(da:-p):/(p:*(1:-p))
        H_step=cross(Xa,w_irls,Xa)
        chol_ok=0
        if (!missing(H_step)) {
            h_diag=diagonal(H_step)
            if (all(h_diag:>0)) {
                L_step=cholesky(H_step)
                if (!missing(L_step)) {
                    if (all(diagonal(L_step):>1e-4:*sqrt(h_diag))) {
                        rhs_step=cross(Xa,w_irls:*z_raw)
                        bn=solvelower(L_step,rhs_step,rk=.,0)
                        if (rk==k & !missing(bn)) bn=solveupper(L_step',bn,rk=.,0)
                        if (rk==k & !missing(bn)) {
                            rhs_ref=quadcross(Xa,w_irls:*(z_raw-Xa*bn))
                            d_bn=solvelower(L_step,rhs_ref,rk=.,0)
                            if (rk==k & !missing(d_bn)) d_bn=solveupper(L_step',d_bn,rk=.,0)
                            if (rk==k & !missing(d_bn)) { bn=bn+d_bn; chol_ok=1; }
                        }
                    }
                }
            }
        }
        if (!chol_ok) {
            sw=sqrt(w_irls); A=Xa:*sw
            /* For a full-rank glm.fit, dqrdc2 retains input order.
               Its rank tolerance is relative to each original column norm. */
            hqrd(A,Q=.,tau=.,R=.)
            if (rows(R)!=k | missing(R) | any(abs(diagonal(R)):<=qtol:*sqrt(colsum(A:^2))'))
                _cdid_nuisance_fail("Covariates in the propensity-score model are collinear on its estimation sample; remove redundant covariates.",459)
            z=z_raw:*sw
            qy=hqrdmultq1t(Q,tau,z)
            bn=solveupper(R,qy,rk=.,0)
            if (rk!=k | missing(bn)) _cdid_nuisance_fail("Nonfinite logit IRLS step.",504)
        }
        eta=Xa*bn; p=_cdid_logit_p(eta)
        dev=_cdid_logit_deviance(da,p,wa)
        h=0
        while (dev>=. | any(p:<=0) | any(p:>=1)) {
            boundary=1; h=h+1
            if (h>maxit) _cdid_nuisance_fail("Logit boundary step-halving failed.",430)
            bn=(bn+b)/2; eta=Xa*bn; p=_cdid_logit_p(eta)
            dev=_cdid_logit_deviance(da,p,wa)
        }
        b=bn
        if (abs(dev-olddev)/(.1+abs(dev))<eps) { converged=1; break; }
        olddev=dev
    }
    if (boundary) _cdid_nuisance_fail("Logit encountered numerical boundary.",430)
    if (!converged) _cdid_nuisance_fail("Logit did not converge under supplied controls.",430)
    /* R 4.5.2 family.c logit_linkinv maps active |eta|>30 to epsilon tails. */
    if (any(abs(eta):>30) | any(p:<=10*epsilon(1)) | any(p:>=1-10*epsilon(1)))
        _cdid_nuisance_fail("Logit fitted probabilities reached numerical boundary.",430)
    out.beta=b
    if (na<n) vals=cdid_logit_values(X*b)
    else vals=cdid_logit_values(eta)
    out.rawp=vals[,1]; out.clippedp=vals[,2]; out.clip_active=vals[,3]
    out.row_score=X:*(w:*(D:-out.rawp))
    out.information=quadcross(Xa,wa:*p:*(1:-p),Xa)
    information_chol=cholesky(out.information)
    if (missing(information_chol) | min(diagonal(information_chol))<=0)
        _cdid_nuisance_fail("Logit information is not positive definite.",459)
    if (no_row_if) {
        out.R=information_chol'; out.pivot=pivot; out.rcond_R=_cdid_nuisance_rcond(out.R)
        if (missing(out.information) | missing(out.R) | out.rcond_R<=0)
            _cdid_nuisance_fail("Logit information is singular/nonfinite.",459)
        out.row_dlogweight=J(0,0,0)
    }
    else {
        hqrd(Xa:*sqrt(wa:*p:*(1:-p)),Q=.,tau=.,R=.)
        out.R=R; out.pivot=pivot; out.rcond_R=_cdid_nuisance_rcond(R)
        if (missing(out.information) | missing(R) | out.rcond_R<=0)
            _cdid_nuisance_fail("Logit information is singular/nonfinite.",459)
        out.row_dlogweight=_cdid_nuisance_gram(R,pivot,out.row_score')'
    }
    out.iterations=it; out.converged=converged; out.deviance=dev
    out.score_max=_cdid_nuisance_max(colsum(out.row_score))
    return(out)
}
real matrix _cdid_nuisance_units(real matrix rowscore, real colvector unit, real scalar fullN)
{
    real matrix ans, info
    real colvector seen, ord
    real scalar n
    n=rows(unit)
    if (n!=rows(rowscore) | cols(unit)!=1 | missing(unit) | fullN<1 |
        fullN>=. | fullN!=floor(fullN) | any(unit:<1) | any(unit:>fullN) | any(unit:!=floor(unit)))
        _cdid_nuisance_fail("Unit IDs must index the explicit original frame.",498)
    ans=J(fullN,cols(rowscore),0)
    if (cols(rowscore)==0 | n==0) return(ans)
    if (n<=fullN) {
        seen=J(fullN,1,0); seen[unit]=J(n,1,1)
        if (sum(seen)==n) {
            ans[unit,]=rowscore
            return(ans)
        }
    }
    ord=order((unit,(1::n)),(1,2))
    info=_cdid_ps_nuis(unit[ord],1)
    ans[unit[ord[info[,1]]],]=panelsum(rowscore[ord,],info)
    return(ans)
}
real matrix cdid_logit_unit_if(struct cdid_logit_result scalar f, real colvector unit, real scalar fullN)
{
    return(fullN:*_cdid_nuisance_units(f.row_dlogweight,unit,fullN))
}
/* Already calendar-step-scaled X, with original-unit IDs and full-frame N.
   Arbitrary integer block IDs are allowed. Prediction rows must name a fitted
   block. FWL avoids an N by block-dummy design. Exact aliases require strict
   reconstruction on fit AND prediction rows. beta uses kept-column order. */
struct cdid_wls_result scalar cdid_pooled_wls_fit(real matrix X, real colvector y,
    real colvector w, real colvector block, real colvector unit, real scalar fullN,
    real scalar intercept, real matrix prediction_X, real colvector prediction_block,
    real scalar rank_tol, | real scalar skip_row_if)
{
    struct cdid_wls_result scalar out
    real scalar n, p, B, b, j, k, rk, err, den, cap, original_norm, residual_norm, n_pred, no_row_if, sorted_b
    real colvector ix, bi, pi, yc, v, qy, tau, coef, tmp, ord, b_sorted
    real matrix C, Q, R, Z, V, map, constants, scales, spans, X_b
    real rowvector keep, drop, pivot, candidate, newdrop, newkeep, cmin_b, cmax_b, col_norms
    no_row_if=(args()>=11 & skip_row_if==1)
    n=rows(X); p=cols(X); n_pred=rows(prediction_X)
    if (n==0 | rows(y)!=n | rows(w)!=n | rows(block)!=n | rows(unit)!=n |
        cols(y)!=1 | cols(w)!=1 | cols(block)!=1 | missing(X) | missing(y) |
        missing(w) | missing(block) | any(w:<=0) | any(block:!=floor(block)) |
        (intercept!=0 & intercept!=1) | cols(prediction_X)!=p |
        n_pred!=rows(prediction_block) | cols(prediction_block)!=1 |
        missing(prediction_X) | missing(prediction_block) | rank_tol<=0 | rank_tol>=1)
        _cdid_nuisance_fail("Invalid pooled WLS inputs or explicit rank tolerance.",498)
    if (!intercept & p==0) _cdid_nuisance_fail("No-intercept outcome design is empty.",459)
    tmp=_cdid_nuisance_units(J(n,0,0),unit,fullN)
    sorted_b=(n==1 ? 1 : all(block[|2\n|]:>=block[|1\n-1|]))
    if (sorted_b) { ord=(1::n); b_sorted=block; }
    else { ord=order((block,(1::n)),(1,2)); b_sorted=block[ord]; }
    spans=_cdid_ps_nuis(b_sorted,1)
    out.block_ids=b_sorted[spans[,1]]; B=rows(spans); bi=J(n,1,0)
    pi=J(n_pred,1,0)
    out.mass=J(B,1,0); out.mean_X=J(B,p,0); out.mean_y=J(B,1,0)
    C=X; yc=y
    for (b=1;b<=B;b++) {
        ix=(sorted_b ? (spans[b,1]..spans[b,2])' : ord[spans[b,1]..spans[b,2]])
        bi[ix]=J(rows(ix),1,b)
        out.mass[b]=sum(w[ix]); out.mean_y[b]=sum(w[ix]:*y[ix])/out.mass[b]
        if (p) {
            X_b=X[ix,]
            out.mean_X[b,]=quadcross(w[ix],X_b)/out.mass[b]
        }
        if (intercept) {
            if (p) C[ix,]=X_b:-out.mean_X[b,]
            yc[ix]=y[ix]:-out.mean_y[b]
            /* The R kernel sets actual block-constant centered columns to
               exact zero; retain them globally if another block varies. */
            if (p) {
                cmin_b=colmin(X_b); cmax_b=colmax(X_b)
                for (j=1;j<=p;j++) if (cmin_b[j]==cmax_b[j]) {
                    C[ix,j]=J(rows(ix),1,0); out.mean_X[b,j]=cmin_b[j]
                }
            }
        }
        if (n_pred) {
            if (B==1 & all(prediction_block:==out.block_ids[1])) pi=J(n_pred,1,1)
            else {
                ix=_cdid_nuisance_index(prediction_block:==out.block_ids[b])
                if (rows(ix)) pi[ix]=J(rows(ix),1,b)
            }
        }
    }
    if (n_pred & any(pi:==0)) _cdid_nuisance_fail("Prediction requires an unfitted block.",459)
    C=C:*sqrt(w); yc=yc:*sqrt(w); keep=J(1,0,0); drop=J(1,0,0)
    cap=64*epsilon(1)*max((1,p+intercept)); out.alias_threshold=cap
    if (p>=2) {
        col_norms=sqrt(colsum(C:^2))
        if (all(col_norms:>0)) {
            hqrd(C,Q=.,tau=.,R=.)
            if (rows(R)==p & !missing(R)) {
                if (all(abs(diagonal(R)):>(16*p*rank_tol):*col_norms')) keep=(1..p)
            }
        }
    }
    /* Original-order near-dependency proposal, followed by the source's
       arithmetic-roundoff check. Non-exact candidates remain in the model. */
    if (cols(keep)==0 & p>0) {
        for (j=1;j<=p;j++) {
            v=C[,j]; original_norm=sqrt(quadcross(v,v))
            if (cols(keep)) {
                hqrdp(C[,keep],Q=.,tau=.,R=.,pivot=.)
                qy=hqrdmultq1t(Q,tau,v); coef=solveupper(R,qy,rk=.,0)
                coef=coef[invorder(pivot)]; v=v-C[,keep]*coef
            }
            residual_norm=sqrt(quadcross(v,v))
            if (original_norm==0 | residual_norm<=rank_tol*original_norm) drop=(drop,j)
            else keep=(keep,j)
        }
    }
    map=J(cols(keep),cols(drop),0); constants=J(B,cols(drop),0)
    if (cols(drop)) {
        if (cols(keep)) {
            hqrdp(C[,keep],Q=.,tau=.,R=.,pivot=.)
            qy=hqrdmultq1t(Q,tau,C[,drop]); map=solveupper(R,qy,rk=.,0)
            map=map[invorder(pivot),]
            qy=hqrdmultq1t(Q,tau,C[,drop]-C[,keep]*map)
            V=solveupper(R,qy,rk=.,0); map=map+V[invorder(pivot),]
        }
        if (intercept) constants=out.mean_X[,drop]-out.mean_X[,keep]*map
        candidate=J(1,0,0); newdrop=J(1,0,0)
        for (j=1;j<=cols(drop);j++) {
            v=X[,drop[j]]-X[,keep]*map[,j]-constants[bi,j]
            scales=abs(X[,drop[j]])+abs(X[,keep])*abs(map[,j])+abs(constants[bi,j])
            den=max((_cdid_nuisance_max(X[,drop[j]]),_cdid_nuisance_max(prediction_X[,drop[j]]),2.2250738585072014e-308))
            err=max(abs(v):/rowmax((scales,J(n,1,den))))
            if (err>cap) { candidate=(candidate,drop[j]); continue; }
            if (n_pred) {
                v=prediction_X[,drop[j]]-prediction_X[,keep]*map[,j]-constants[pi,j]
                scales=abs(prediction_X[,drop[j]])+abs(prediction_X[,keep])*abs(map[,j])+abs(constants[pi,j])
                err=max(abs(v):/rowmax((scales,J(rows(v),1,den))))
                if (err>cap) _cdid_nuisance_fail("Outcome alias does not preserve required prediction support.",459)
            }
            newdrop=(newdrop,drop[j])
        }
        if (cols(candidate)) {
            /* Preserve previously checked alias maps with zero loadings on
               added candidates, as in the source preprocessing. */
            newkeep=sort((keep,candidate)',1)'
            Z=J(cols(newkeep),cols(newdrop),0)
            if (cols(keep)) for (j=1;j<=cols(newdrop);j++) {
                ix=_cdid_nuisance_index(drop:==newdrop[j])
                for (b=1;b<=cols(keep);b++) Z[_cdid_nuisance_index(newkeep:==keep[b]),j]=map[b,ix]
            }
            keep=newkeep; map=Z; drop=newdrop
            constants=J(B,cols(drop),0)
            if (intercept) constants=out.mean_X[,drop]-out.mean_X[,keep]*map
        }
    }
    k=cols(keep)
    if (!intercept & k==0) _cdid_nuisance_fail("Outcome design has no retained columns.",459)
    out.beta=J(0,1,0); out.R=J(0,0,0); out.pivot=J(1,0,0); out.rcond_R=1
    if (k) {
        Z=C[,keep]; hqrdp(Z,Q=.,tau=.,R=.,pivot=.)
        if (rows(R)!=k | missing(R) | min(abs(diagonal(R)))<=rank_tol*max(abs(diagonal(R))))
            _cdid_nuisance_fail("Pooled WLS final rank refused; no ridge or parameter rescue.",459)
        qy=hqrdmultq1t(Q,tau,yc); coef=solveupper(R,qy,rk=.,0)
        out.beta=coef[invorder(pivot)]; out.R=R; out.pivot=pivot
        out.rcond_R=(no_row_if ? 1 : _cdid_nuisance_rcond(R))
    }
    out.alpha=(intercept ? out.mean_y-out.mean_X[,keep]*out.beta : J(B,1,0))
    out.fitted=out.alpha[bi,1]+X[,keep]*out.beta; out.residual=y-out.fitted
    out.prediction=out.alpha[pi,1]+prediction_X[,keep]*out.beta
    if (no_row_if) {
        out.row_score=J(0,0,0)
        out.row_dlogweight_beta=J(0,0,0)
        out.row_alpha_direct=J(0,1,0)
    }
    else {
        out.row_score=(intercept ? X[,keep]-out.mean_X[bi,keep] : X[,keep]):*(w:*out.residual)
        out.row_dlogweight_beta=_cdid_nuisance_gram(out.R,out.pivot,out.row_score')'
        out.row_alpha_direct=(intercept ? w:*out.residual:/out.mass[bi,1] : J(n,1,0))
    }
    out.kept=keep; out.dropped=drop; out.alias_map=map; out.alias_constants=constants
    out.n_rows=n; out.n_blocks=B; out.fullN=fullN; out.intercept=intercept; out.row_block=block
    out.rss=sum(w:*out.residual:^2); out.rank_tol=rank_tol
    if (missing(out.beta) | missing(out.alpha) | missing(out.residual) | missing(out.prediction) |
        (!no_row_if & missing(out.row_dlogweight_beta))) _cdid_nuisance_fail("Nonfinite WLS fit or derivative.",504)
    return(out)
}
real matrix cdid_wls_unit_beta_if(struct cdid_wls_result scalar f, real colvector unit)
{
    return(f.fullN:*_cdid_nuisance_units(f.row_dlogweight_beta,unit,f.fullN))
}
/* One intercept block at a time, with the common-beta contribution. */
real colvector cdid_wls_unit_alpha_if(struct cdid_wls_result scalar f,
    real colvector unit, real colvector block, real scalar block_id)
{
    real colvector ix, direct
    if (rows(block)!=f.n_rows | cols(block)!=1 | missing(block) | any(block:!=f.row_block))
        _cdid_nuisance_fail("Block rows do not match fit.",498)
    ix=_cdid_nuisance_index(f.block_ids:==block_id)
    if (rows(ix)!=1) _cdid_nuisance_fail("Unknown intercept block.",498)
    if (!f.intercept) return(J(f.fullN,1,0))
    direct=f.row_alpha_direct:*(block:==block_id)
    return(f.fullN:*_cdid_nuisance_units(direct,unit,f.fullN)-
        cdid_wls_unit_beta_if(f,unit)*f.mean_X[ix,f.kept]')
}
/* Complete fixed-weight derivative for caller-supplied logweight directions.
   Output rows are kept beta then block alpha; columns are directions.
   Caller composes upstream estimated-weight channels using these loadings. */
real matrix cdid_wls_logweight_direction(struct cdid_wls_result scalar f,
    real colvector block, real matrix directions)
{
    real matrix db, da, spans
    real colvector ord
    if (rows(directions)!=f.n_rows | rows(block)!=f.n_rows | cols(block)!=1 |
        missing(block) | any(block:!=f.row_block) | missing(directions))
        _cdid_nuisance_fail("Invalid logweight direction rows.",498)
    db=quadcross(f.row_dlogweight_beta,directions); da=J(f.n_blocks,cols(directions),0)
    if (f.intercept) {
        ord=order((block,(1::f.n_rows)),(1,2))
        spans=_cdid_ps_nuis(block[ord],1)
        da=panelsum((f.row_alpha_direct:*directions)[ord,],spans)-f.mean_X[,f.kept]*db
    }
    return(db\da)
}
end
