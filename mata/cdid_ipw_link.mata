*! Baseline native weighted-pair IPW link and full-unit IF; Stata 14.2+.
*! Reference: cdid 1.1.0.0, commit .
*! Load cdid_nuisance.mata before this file. No population/catalogue decisions here.
version 14.2
mata:
struct cdid_ipw_frame {
    real scalar n_A, n_F, n_cohorts
    real colvector group, weight, a_to_full, unit_to_a
    real colvector group_count, group_mass, q
}
struct cdid_ipw_nuisance {
    real scalar n_A, n_F, n_cohorts, n_pair_all, n_pair_risk, p
    real scalar sum_a, sum_b, mean_a, mean_b, epsilon, maxit
    real scalar iterations, deviance, score_max
    real colvector risk_a, in_risk, risk_group, risk_weight
    real colvector pair_count, pair_mass, rates, rates_c, rates_active, q
    real colvector beta, raw_probability, odds, odds_active, fit_weight
    real matrix X, information, H, H_chol, A_beta_r, row_score
    real colvector risk_pos_a
}
struct cdid_ipw_link {
    real scalar estimate, treated_mean, control_mean, n_A, n_F, n_cohorts, p
    real colvector G_beta, G_r, loading_beta, loading_rate
    real colvector rate_coef_s0, rate_coef_s1, beta, rates, q
}
struct cdid_ipw_if_chunk {
    real colvector full_index, phi
    real matrix components
    real scalar scale_A_to_F
}
void _cdid_ipw_fail(string scalar msg, real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}
real colvector _cdid_ipw_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
/* Consumed rates and strict floor derivative mask; equality is inactive. */
real matrix cdid_ipw_rate_values(real colvector rates)
{
    if (missing(rates) | any(rates:<=0) | any(rates:>1))
        _cdid_ipw_fail("Pair observation rates must be finite and in (0,1].",498)
    return((rowmax((rates,J(rows(rates),1,1e-6))),rates:>1e-6))
}
real colvector cdid_ipw_odds_active(real colvector raw_probability)
{
    if (missing(raw_probability) | any(raw_probability:<0) | any(raw_probability:>1))
        _cdid_ipw_fail("Raw probabilities must be finite and in [0,1].",498)
    return(raw_probability:<(1-1e-6))
}
/* Exact source arithmetic for raw probabilities and consumed treatment odds. */
real matrix _cdid_ipw_logit_values(real colvector eta)
{
    real matrix values
    real colvector odds, active, ix
    real scalar cap, cap_odds
    values=cdid_logit_values(eta)
    active=cdid_ipw_odds_active(values[,1])
    cap=1-1e-6; cap_odds=cap/(1-cap)
    odds=J(rows(eta),1,cap_odds)
    ix=_cdid_ipw_index(active)
    if(rows(ix)) odds[ix]=exp(eta[ix])
    if(missing(odds) | any(odds:<=0))
        _cdid_ipw_fail("Propensity odds must be finite and strictly positive.",504)
    return((values[,1],odds,active))
}
/* Ordered group map: 1=target, 2..m=declared controls, 0=known off-risk.
   The population layer owns actual cohort labels, including never-treated.
   Unknown cohort is forbidden on A. Unknown off-A metadata is not consumed.
   The one-to-one a_to_full map preserves singleton and frame-only zero rows. */
struct cdid_ipw_frame scalar cdid_ipw_frame_create(real colvector a_group,
    real colvector a_cohort_known, real colvector a_weight,
    real colvector a_to_full, real scalar n_F, real scalar n_cohorts)
{
    struct cdid_ipw_frame scalar out
    real scalar n, j
    real colvector ix
    n=rows(a_group)
    if(n<1 | cols(a_group)!=1 | rows(a_cohort_known)!=n | cols(a_cohort_known)!=1 |
        rows(a_weight)!=n | cols(a_weight)!=1 | rows(a_to_full)!=n | cols(a_to_full)!=1 |
        missing(a_group) | missing(a_cohort_known) | missing(a_weight) | missing(a_to_full) |
        any(a_cohort_known:!=1) | any(a_weight:<=0) |
        n_F>=. | n_F<n | n_F!=floor(n_F) | n_F>2^31-1 |
        n_cohorts>=. | n_cohorts<2 | n_cohorts!=floor(n_cohorts) |
        any(a_group:!=floor(a_group)) | any(a_group:<0) | any(a_group:>n_cohorts) |
        any(a_to_full:!=floor(a_to_full)) | any(a_to_full:<1) | any(a_to_full:>n_F))
        _cdid_ipw_fail("Require known A cohorts, ordered groups, positive A weights and an explicit original-F index map.",498)
    if(n_cohorts>n)
        _cdid_ipw_fail("Every declared cohort requires positive A count; declared count exceeds A size.",459)
    if(sum(a_weight)>=. | sum(a_weight)<=0)
        _cdid_ipw_fail("Total A-unit weight must be finite and positive.",498)
    if(n>1) {
        if(!all(a_to_full[|2\n|]:>a_to_full[|1\n-1|])) {
            if(rows(uniqrows(sort(a_to_full,1)))!=n)
                _cdid_ipw_fail("a_to_full must be one-to-one.",498)
        }
    }
    out.n_A=n; out.n_F=n_F; out.n_cohorts=n_cohorts
    out.group=a_group; out.weight=a_weight; out.a_to_full=a_to_full
    out.unit_to_a=J(n_F,1,0); out.unit_to_a[a_to_full]=(1..n)'
    out.group_count=J(n_cohorts,1,0); out.group_mass=J(n_cohorts,1,0)
    for(j=1;j<=n_cohorts;j++) {
        ix=_cdid_ipw_index(a_group:==j)
        out.group_count[j]=rows(ix); out.group_mass[j]=sum(a_weight[ix])
    }
    if(any(out.group_count:<1) | missing(out.group_mass) | any(out.group_mass:<=0))
        _cdid_ipw_fail("Every declared target/control cohort requires positive A count and weight mass; none is dropped.",459)
    out.q=out.group_mass/n
    if(missing(out.q) | any(out.q:<=0))
        _cdid_ipw_fail("Invalid A-normalized cohort weight masses.",504)
    return(out)
}
/* Outcome-free propensity/rate stage. pair_a is the strictly increasing compact
   A index of ALL supplied pair observers; pair_X rows have the same order.
   Off-risk pair observers are retained in n_pair_all and never enter the fit.
   epsilon and maxit are explicit; source public settings are 1e-12 and 50.
   Reuse this result with cdid_ipw_link_outcome for the same nuisance recipe. */
struct cdid_ipw_nuisance scalar cdid_ipw_nuisance_prepare(struct cdid_ipw_frame scalar frame,
    real colvector pair_a, real matrix pair_X, real scalar epsilon, real scalar maxit)
{
    struct cdid_ipw_nuisance scalar out
    struct cdid_logit_result scalar fitted
    real scalar b, j, n, p
    real colvector gi, ix, D, weights, residual
    real matrix rv, lv
    n=rows(pair_a); p=cols(pair_X)
    if(n<1) _cdid_ipw_fail("Requested date pair has zero A-unit observations.",459)
    if(cols(pair_a)!=1 | rows(pair_X)!=n | p<1 | missing(pair_a) | missing(pair_X) |
        any(pair_a:!=floor(pair_a)) | any(pair_a:<1) | any(pair_a:>frame.n_A) |
        epsilon>=. | epsilon<=0 | maxit>=. | maxit<1 | maxit!=floor(maxit) | maxit>2^31-1)
        _cdid_ipw_fail("Invalid compact pair data or explicit GLM controls.",498)
    if(n>1) if(any(pair_a[2..n]:<=pair_a[1..n-1]))
        _cdid_ipw_fail("Compact pair A indices must be strictly increasing and unique.",498)
    gi=frame.group[pair_a]
    out.in_risk=_cdid_ipw_index(gi:>0)
    b=rows(out.in_risk)
    if(b<1) _cdid_ipw_fail("Requested pair has zero declared risk-set observations.",459)
    out.n_A=frame.n_A; out.n_F=frame.n_F; out.n_cohorts=frame.n_cohorts
    out.n_pair_all=n; out.n_pair_risk=b; out.p=p
    out.risk_a=pair_a[out.in_risk]; out.risk_group=gi[out.in_risk]
    out.risk_pos_a=J(out.n_A,1,0); out.risk_pos_a[out.risk_a]=(1::b)
    out.risk_weight=frame.weight[out.risk_a]; out.X=pair_X[out.in_risk,]
    out.pair_count=J(out.n_cohorts,1,0); out.pair_mass=J(out.n_cohorts,1,0)
    for(j=1;j<=out.n_cohorts;j++) {
        ix=_cdid_ipw_index(out.risk_group:==j)
        out.pair_count[j]=rows(ix); out.pair_mass[j]=sum(out.risk_weight[ix])
    }
    if(any(out.pair_count:<1) | missing(out.pair_mass) | any(out.pair_mass:<=0))
        _cdid_ipw_fail("Every declared cohort requires positive observed-pair count and mass; none is dropped.",459)
    out.rates=out.pair_mass:/frame.group_mass
    rv=cdid_ipw_rate_values(out.rates)
    out.rates_c=rv[,1]; out.rates_active=rv[,2]; out.q=frame.q
    D=out.risk_group:==1
    out.fit_weight=out.risk_weight:/out.rates_c[out.risk_group]
    if(missing(out.fit_weight) | any(out.fit_weight:<=0))
        _cdid_ipw_fail("Inverse-rate propensity weights must be finite and positive.",504)
    fitted=cdid_logit_fit(out.X,D,out.fit_weight,epsilon,maxit,1)
    out.beta=fitted.beta; out.iterations=fitted.iterations
    out.deviance=fitted.deviance; out.score_max=fitted.score_max
    out.epsilon=epsilon; out.maxit=maxit
    lv=_cdid_ipw_logit_values(out.X*out.beta)
    out.raw_probability=lv[,1]; out.odds=lv[,2]; out.odds_active=lv[,3]
    ix=_cdid_ipw_index(D)
    out.sum_a=sum(out.fit_weight[ix])
    ix=_cdid_ipw_index(!D)
    out.sum_b=sum(out.fit_weight[ix]:*out.odds[ix])
    out.mean_a=out.sum_a/out.n_A; out.mean_b=out.sum_b/out.n_A
    if(out.sum_a<=0 | out.sum_a>=. | out.sum_b<=0 | out.sum_b>=. |
        out.mean_a<=0 | out.mean_a>=. | out.mean_b<=0 | out.mean_b>=.)
        _cdid_ipw_fail("Invalid treated/control outcome normalization denominator.",504)
    residual=D:-out.raw_probability
    /* Borrow the nuisance result instead of repeating its O(np^2) cross-product. */
    out.row_score=fitted.row_score
    out.information=fitted.information
    out.H=out.information/out.n_A
    if(missing(out.H)) _cdid_ipw_fail("Nonfinite weighted propensity Hessian.",504)
    out.H_chol=cholesky(out.H)
    if(missing(out.H_chol) | any(diagonal(out.H_chol):<=0))
        _cdid_ipw_fail("Weighted propensity Hessian is not positive definite.",459)
    out.A_beta_r=J(p,out.n_cohorts,0)
    for(j=1;j<=out.n_cohorts;j++) if(out.rates_active[j]) {
        ix=_cdid_ipw_index(out.risk_group:==j)
        weights=out.risk_weight[ix]:*residual[ix]:/(out.rates_c[j]^2)
        out.A_beta_r[,j]=-quadcross(out.X[ix,],weights)/out.n_A
    }
    if(missing(out.A_beta_r)) _cdid_ipw_fail("Nonfinite ALL-cohort score-rate Jacobian.",504)
    return(out)
}
/* Outcome-only step; no propensity refit. pair_dy uses ALL pair observer order. */
struct cdid_ipw_link scalar cdid_ipw_link_outcome(struct cdid_ipw_nuisance scalar nu,
    real colvector pair_dy)
{
    struct cdid_ipw_link scalar out
    real scalar j, rk
    real colvector dy, ix_t, ix_c, ix, wt, res, Zsolve
    real matrix Xc
    if(rows(pair_dy)!=nu.n_pair_all | cols(pair_dy)!=1 | missing(pair_dy))
        _cdid_ipw_fail("Finite outcome changes must match all compact pair observers.",498)
    dy=pair_dy[nu.in_risk]
    ix_t=_cdid_ipw_index(nu.risk_group:==1); ix_c=_cdid_ipw_index(nu.risk_group:!=1)
    out.n_A=nu.n_A; out.n_F=nu.n_F; out.n_cohorts=nu.n_cohorts; out.p=nu.p
    out.treated_mean=sum(nu.fit_weight[ix_t]:*dy[ix_t])/nu.sum_a
    wt=nu.fit_weight[ix_c]:*nu.odds[ix_c]
    out.control_mean=sum(wt:*dy[ix_c])/nu.sum_b
    out.estimate=out.treated_mean-out.control_mean
    res=dy[ix_c]:-out.control_mean; Xc=nu.X[ix_c,]
    out.G_beta=-quadcross(Xc,nu.odds_active[ix_c]:*wt:*res)/nu.sum_b
    out.G_r=J(nu.n_cohorts,1,0)
    for(j=2;j<=nu.n_cohorts;j++) if(nu.rates_active[j]) {
        ix=_cdid_ipw_index(nu.risk_group:==j)
        out.G_r[j]=sum(nu.risk_weight[ix]:*nu.odds[ix]:*(dy[ix]:-out.control_mean))/
            (nu.rates_c[j]^2*nu.sum_b)
    }
    Zsolve=solvelower(nu.H_chol,out.G_beta,rk=.,0)
    out.loading_beta=solveupper(nu.H_chol',Zsolve,rk=.,0)
    if(rk!=nu.p) _cdid_ipw_fail("Propensity Hessian loading solve refused.",459)
    out.loading_rate=out.G_r+nu.A_beta_r'*out.loading_beta
    out.rate_coef_s0=-(nu.rates:/nu.q):*out.loading_rate
    out.rate_coef_s1=((1:-nu.rates):/nu.q):*out.loading_rate
    out.beta=nu.beta; out.rates=nu.rates; out.q=nu.q
    if(out.estimate>=. | out.treated_mean>=. | out.control_mean>=. |
        missing(out.G_beta) | missing(out.G_r) | missing(out.loading_beta) |
        missing(out.loading_rate) | missing(out.rate_coef_s0) | missing(out.rate_coef_s1))
        _cdid_ipw_fail("Nonfinite link means, partials, or replay loadings.",504)
    return(out)
}
real scalar _cdid_ipw_lower_bound(real colvector indices, real scalar value)
{
    real scalar lo, hi, mid
    lo=1; hi=rows(indices)
    while(lo<=hi) {
        mid=floor(lo+(hi-lo)/2)
        if(indices[mid]<value) lo=mid+1
        else hi=mid-1
    }
    return(lo)
}
real colvector _cdid_ipw_risk_pos(struct cdid_ipw_nuisance scalar nu, real colvector ai)
{
    real colvector pos
    real scalar i, k, nobs
    if(rows(nu.risk_pos_a)==nu.n_A) return(nu.risk_pos_a[ai])
    pos=J(rows(ai),1,0); nobs=rows(nu.risk_a)
    for(i=1;i<=rows(ai);i++) {
        k=_cdid_ipw_lower_bound(nu.risk_a,ai[i])
        if(k<=nobs) if(nu.risk_a[k]==ai[i]) pos[i]=k
    }
    return(pos)
}
/* Emit only one requested original-F chunk, including exact off-A zero rows.
   Components are rate, raw beta_score, treated, control, in that order.
   The same compact pair_dy must accompany point and replay calls. */
struct cdid_ipw_if_chunk scalar cdid_ipw_link_if_chunk(struct cdid_ipw_frame scalar frame,
    struct cdid_ipw_nuisance scalar nu, struct cdid_ipw_link scalar link,
    real colvector pair_dy, real scalar first_full, real scalar last_full)
{
    struct cdid_ipw_if_chunk scalar out
    real scalar b
    real colvector ai_all, ix_a, ai, jg, ix_g, row_g, ai_g, j_g, k_g, ix_r
    real colvector row_r, ai_r, j_r, k_r, dy_r, ix_t, ix_c, kt, kc
    if(first_full>=. | last_full>=. | first_full!=floor(first_full) |
        last_full!=floor(last_full) | first_full<1 | last_full>frame.n_F | first_full>last_full |
        rows(pair_dy)!=nu.n_pair_all | cols(pair_dy)!=1 | missing(pair_dy))
        _cdid_ipw_fail("Invalid full-frame chunk or compact outcome changes.",498)
    if(frame.n_A!=nu.n_A | frame.n_F!=nu.n_F | frame.n_cohorts!=nu.n_cohorts |
        link.n_A!=nu.n_A | link.n_F!=nu.n_F | link.n_cohorts!=nu.n_cohorts |
        link.p!=nu.p | any(link.beta:!=nu.beta) | any(link.rates:!=nu.rates) | any(link.q:!=nu.q) |
        any(frame.q:!=nu.q))
        _cdid_ipw_fail("Frame, nuisance and link design block do not match.",498)
    b=last_full-first_full+1; out.full_index=(first_full..last_full)'
    out.components=J(b,4,0); out.scale_A_to_F=frame.n_F/frame.n_A
    ai_all=frame.unit_to_a[out.full_index]
    ix_a=_cdid_ipw_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]; jg=frame.group[ai]
        ix_g=_cdid_ipw_index(jg:>0)
        if(rows(ix_g)) {
            row_g=ix_a[ix_g]; ai_g=ai[ix_g]; j_g=jg[ix_g]
            out.components[row_g,1]=frame.weight[ai_g]:*link.rate_coef_s0[j_g]
            k_g=_cdid_ipw_risk_pos(nu,ai_g)
            ix_r=_cdid_ipw_index(k_g:>0)
            if(rows(ix_r)) {
                row_r=row_g[ix_r]; ai_r=ai_g[ix_r]; j_r=j_g[ix_r]; k_r=k_g[ix_r]
                dy_r=pair_dy[nu.in_risk[k_r]]
                out.components[row_r,1]=frame.weight[ai_r]:*link.rate_coef_s1[j_r]
                out.components[row_r,2]=nu.fit_weight[k_r]:*((j_r:==1):-nu.raw_probability[k_r]):*(nu.X[k_r,]*link.loading_beta)
                ix_t=_cdid_ipw_index(j_r:==1)
                if(rows(ix_t)) {
                    kt=k_r[ix_t]
                    out.components[row_r[ix_t],3]=(nu.fit_weight[kt]:/nu.mean_a):*(dy_r[ix_t]:-link.treated_mean)
                }
                ix_c=_cdid_ipw_index(j_r:!=1)
                if(rows(ix_c)) {
                    kc=k_r[ix_c]
                    out.components[row_r[ix_c],4]=-(nu.fit_weight[kc]:*nu.odds[kc]:/nu.mean_b):*(dy_r[ix_c]:-link.control_mean)
                }
            }
        }
    }
    out.components=out.scale_A_to_F:*out.components
    out.phi=rowsum(out.components)
    if(missing(out.components) | missing(out.phi))
        _cdid_ipw_fail("Nonfinite original-unit IF chunk.",504)
    return(out)
}
/* Optional nuisance IF chunk, with rate columns then beta columns.
   Useful for downstream estimated-weight channels and independent derivative
   checks; no N by all-links allocation. The raw score is never clipped. */
real matrix cdid_ipw_nuisance_if_chunk(struct cdid_ipw_frame scalar frame,
    struct cdid_ipw_nuisance scalar nu, real scalar first_full, real scalar last_full)
{
    real scalar b, j, rk
    real colvector ai_all, ix_a, ai, jg, ix_g, row_g, ai_g, j_g, k_g, ix_r, row_r, ai_r, j_r, k_r, ix_j
    real matrix rate_if, forcing, beta_if
    if(first_full>=. | last_full>=. | first_full!=floor(first_full) |
        last_full!=floor(last_full) | first_full<1 | last_full>frame.n_F | first_full>last_full |
        frame.n_A!=nu.n_A | frame.n_F!=nu.n_F | frame.n_cohorts!=nu.n_cohorts | any(frame.q:!=nu.q))
        _cdid_ipw_fail("Invalid nuisance IF frame or chunk.",498)
    b=last_full-first_full+1
    rate_if=J(b,nu.n_cohorts,0); forcing=J(b,nu.p,0)
    ai_all=frame.unit_to_a[(first_full..last_full)']
    ix_a=_cdid_ipw_index(ai_all:>0)
    if(rows(ix_a)) {
        ai=ai_all[ix_a]; jg=frame.group[ai]
        ix_g=_cdid_ipw_index(jg:>0)
        if(rows(ix_g)) {
            row_g=ix_a[ix_g]; ai_g=ai[ix_g]; j_g=jg[ix_g]
            for(j=1;j<=nu.n_cohorts;j++) {
                ix_j=_cdid_ipw_index(j_g:==j)
                if(rows(ix_j)) rate_if[row_g[ix_j],j]=-frame.weight[ai_g[ix_j]]:*(nu.rates[j]/nu.q[j])
            }
            k_g=_cdid_ipw_risk_pos(nu,ai_g)
            ix_r=_cdid_ipw_index(k_g:>0)
            if(rows(ix_r)) {
                row_r=row_g[ix_r]; ai_r=ai_g[ix_r]; j_r=j_g[ix_r]; k_r=k_g[ix_r]
                for(j=1;j<=nu.n_cohorts;j++) {
                    ix_j=_cdid_ipw_index(j_r:==j)
                    if(rows(ix_j)) rate_if[row_r[ix_j],j]=frame.weight[ai_r[ix_j]]:*((1-nu.rates[j])/nu.q[j])
                }
                forcing[row_r,]=nu.row_score[k_r,]
            }
        }
    }
    forcing=forcing+rate_if*nu.A_beta_r'
    beta_if=solvelower(nu.H_chol,forcing',rk=.,0)
    beta_if=solveupper(nu.H_chol',beta_if,rk=.,0)'
    if(rk!=nu.p | missing(beta_if) | missing(rate_if))
        _cdid_ipw_fail("Nonfinite nuisance IF loading solve.",504)
    return((frame.n_F/frame.n_A):*(rate_if,beta_if))
}
end
