*! Native cohort-share aggregation primitive for cdid; current StataNow 19.5.
*! Port of pinned R cdid 1.1.0.0, commit ,
*! R/chained_ipw_pipeline.R: ipw_pipeline_aggregation_selection (l.4379-4451),
*! aggregate_chained_ipw_targets (l.4453-4834) and the supported-subvector
*! inference wiring of ipw_pipeline_supported_inference (l.4257-4328).
*! Conventional policy only (R `enclosure = NULL`): no interval certificates.
*! Covariance, multiplier and flag laws come from the cdid_inference.mata,
*! which must be compiled first (struct cdid_infer_* are referenced here).
*! No refit, no RNG, no N-by-all-moments IF and no draw tensor are created.
version 19.5
mata:
/* ---------------------------------------------------------------------------
   Inputs. Ledger rows are the FULL declared reporting ledger (post, pre and
   reference rows) in the frontend's order; frame vectors are in original-F
   unit order including frame-only/zero-IF units. cohort=. with cohort_known=1
   is known never-treated (R +Inf); cohort_known=0 is an unknown
   frame-only cohort (R NA). w may be missing only under share_law
   "first_calendar_masses" (R: V = w_t1 absent where no share reads it).
   Selector scalars: balance_e=. is R NULL; min_e=. is -Inf; max_e=. is +Inf.
   bstrap=. is the private R NULL flag law (inference translation).
   --------------------------------------------------------------------------- */
struct cdid_agg_ledger {
    real colvector g, time, ref, base_t
    string colvector target_id, role, support_status
}
struct cdid_agg_frame {
    real scalar n, C, clustered
    real colvector cohort, cohort_known, w, unit_cluster
    string colvector frame_id
    string scalar share_law
    real colvector share_excluded
}
struct cdid_agg_request {
    string scalar type
    real scalar na_rm, balance_e, min_e, max_e, alpha, bstrap, cband, B
    real colvector calendar, cohort_grid
}
struct cdid_agg_selection {
    string scalar type, convention, scope, warning
    real scalar balance_e, min_e, max_e, nondefault, applied, ignored, M, R
    real colvector keys, grid, horizon, selected, include
    real colvector eligible_cohorts, excluded_cohorts, event_window, declared_event_times
    real matrix sets
}
struct cdid_agg_shares {
    real scalar K, shares_known, first_masses, incomplete_masses, n_unknown
    real colvector cohorts, p, mass_complete
    string colvector mass_missing
    string scalar share_reason
    real matrix P
}
struct cdid_agg_inference {
    string scalar status, inference_status, reason, flag_law, band_family, warning
    real scalar M, n_available, n_F, C, bootstrap, simultaneous, draws_used, n_draws
    real scalar normal_critical, band_critical
    real colvector available, estimates, se, analytic_se, z_stat, p_value
    real colvector ci_lower, ci_upper, band_lower, band_upper
    real matrix covariance, Omega
    struct cdid_infer_result scalar compact
    struct cdid_if_summary scalar summary
}
struct cdid_agg_overall {
    string scalar status, reason, inference_kind, draws_reused_from
    real scalar estimate, se, ci_lower, ci_upper, share_dependent, share_denominator
    real colvector if_full, weights, share_gradient, intended, used, excluded
    struct cdid_agg_inference scalar inference
}
struct cdid_agg_decomposition {
    real scalar computed, shares_known
    real matrix Omega_total, Omega_fixed, Omega_share, Omega_cross, Omega_theta_p
}
struct cdid_agg_result {
    string scalar type, status, target_definition, flag_law, released_warning
    real scalar na_rm, n, R, M, K
    struct cdid_agg_selection scalar selectors
    struct cdid_agg_shares scalar shares
    real colvector keys
    string colvector labels, cell_reason, bin_role, bin_support_status
    real matrix intended, used, excluded
    real colvector n_intended, n_used, cell_available, restricted
    real colvector share_dependent, share_denominator, estimates
    real matrix W, D, if_total, if_fixed, if_share
    struct cdid_agg_inference scalar inference, fixed_inference
    struct cdid_agg_decomposition scalar decomposition
    struct cdid_agg_overall scalar overall
}
void _cdid_agg_fail(string scalar cls, string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_agg_index(real matrix predicate)
{
    if (length(predicate)==0) return(J(0,1,.))
    return(vec(selectindex(predicate)))
}
real colvector _cdid_agg_subset(real colvector x, real colvector mask)
{
    real colvector ix
    if (rows(x)==0) return(J(0,1,.))
    ix=_cdid_agg_index(mask)
    if (rows(ix)==0) return(J(0,1,.))
    return(x[ix])
}
/* Sorted unique finite values (R sort(unique(x)) after dropping NA/Inf). */
real colvector _cdid_agg_uniq(real colvector x)
{
    real colvector z
    z=_cdid_agg_subset(x,x:<.)
    if (rows(z)==0) return(J(0,1,.))
    return(uniqrows(sort(z,1)))
}
real colvector _cdid_agg_lookup(real colvector x, real colvector table)
{
    real colvector ans
    real scalar j, lo, hi, mid
    ans=J(rows(x),1,0)
    for (j=1;j<=rows(x);j++) {
        lo=1; hi=rows(table)
        while (lo<=hi) {
            mid=floor((lo+hi)/2)
            if (table[mid]==x[j]) { ans[j]=mid; break; }
            if (table[mid]<x[j]) lo=mid+1
            else hi=mid-1
        }
    }
    return(ans)
}
string scalar _cdid_agg_join(string colvector s, string scalar sep)
{
    string scalar out
    real scalar i
    out=""
    for (i=1;i<=rows(s);i++) out=out+(i>1 ? sep : "")+s[i]
    return(out)
}
string colvector _cdid_agg_uniqstr(string colvector s)
{
    string colvector out
    real scalar i
    out=J(0,1,"")
    for (i=1;i<=rows(s);i++) if (!anyof(out,s[i])) out=out\s[i]
    return(out)
}
/* R as.character() of a double for integers and ordinary decimals
   ("2004", "-1", "0.5"); Stata's bare leading "." is restored to "0.". */
string scalar _cdid_agg_fmt(real scalar x)
{
    string scalar s
    if (x>=.) return("NA")
    s=strtrim(sprintf("%21.15g",x))
    if (substr(s,1,1)==".") s="0"+s
    else if (substr(s,1,2)=="-.") s="-0"+substr(s,2,.)
    return(s)
}
string colvector _cdid_agg_fmtv(real colvector x)
{
    string colvector s
    real scalar i
    s=J(rows(x),1,"")
    for (i=1;i<=rows(x);i++) s[i]=_cdid_agg_fmt(x[i])
    return(s)
}
/* ---------------------------------------------------------------------------
   Released did 2.5.1 selectors bound to the full declared ledger
   (R l.4379-4451). Intended sets are formed BEFORE any missing-cell handling.
   simple/group: g<=t and recoded-grid horizon pos(t)-pos(g)<=max_e;
   dynamic: e=t-g in ORIGINAL date units, balance_e window, then [min_e,max_e];
   calendar: selectors ignored with the released warning text retained.
   --------------------------------------------------------------------------- */
struct cdid_agg_selection scalar _cdid_agg_selection_plan(string scalar type, real colvector g,
    real colvector tt, real scalar balance_e, real scalar min_e, real scalar max_e,
    real colvector calendar, real colvector cohort_grid)
{
    struct cdid_agg_selection scalar s
    real colvector post, keep, e_all, eseq, pos_t, pos_g, hit
    real scalar R, r, max_t, min_grid
    R=rows(g)
    if (type!="simple" & type!="group" & type!="dynamic" & type!="calendar")
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","type must be simple, group, dynamic or calendar.",498)
    if (R<1 | cols(g)!=1 | rows(tt)!=R | cols(tt)!=1 | missing(g) | missing(tt))
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","The declared table of cells needs one finite cohort and date per row.",498)
    if (balance_e<. & (balance_e<0 | balance_e!=floor(balance_e)))
        _cdid_agg_fail("cdid_ipw_pipeline_selector_error","balance_e must be missing (NULL) or a nonnegative whole number.",498)
    s.type=type; s.balance_e=balance_e; s.min_e=min_e; s.max_e=max_e; s.R=R
    s.nondefault=(balance_e<. | min_e<. | max_e<.)
    s.applied=s.nondefault & type!="calendar"; s.ignored=s.nondefault & type=="calendar"
    s.warning=""
    s.grid=_cdid_agg_uniq(calendar\cohort_grid\g\tt)
    s.scope="selectors choose among already reported ledger cells; the estimation catalog, moments, weights and IFs are unchanged"
    s.horizon=J(R,1,.); s.include=J(R,1,1); s.selected=J(R,1,0)
    s.event_window=J(0,1,.); s.declared_event_times=J(0,1,.); s.excluded_cohorts=J(0,1,.)
    post=_cdid_agg_index(tt:>=g)
    if (type=="calendar") {
        s.convention="released did 2.5.1: min_e, max_e and balance_e are ignored for calendar aggregation (warning)"
        if (s.nondefault) s.warning="mine(), maxe() and balancee() are ignored by estat calendar, as by R's aggte(); the calendar-time effects are not restricted."
        s.keys=(rows(post) ? _cdid_agg_uniq(tt[post]) : J(0,1,.))
        if (rows(post)) s.selected[post]=J(rows(post),1,1)
        s.eligible_cohorts=(rows(post) ? _cdid_agg_uniq(g[post]) : J(0,1,.))
    }
    else if (type=="simple" | type=="group") {
        s.convention=(type=="simple" ?
            "released did 2.5.1: cells with g <= t and recoded-grid horizon pos(t) - pos(g) <= max_e; min_e and balance_e have no effect" :
            "released did 2.5.1: per cohort, cells with g <= t and recoded-grid horizon pos(t) - pos(g) <= max_e; min_e and balance_e have no effect")
        pos_t=_cdid_agg_lookup(tt,s.grid); pos_g=_cdid_agg_lookup(g,s.grid)
        s.horizon=pos_t-pos_g
        keep=_cdid_agg_index((tt:>=g):&(max_e>=. ? J(R,1,1) : s.horizon:<=max_e))
        if (rows(keep)) s.selected[keep]=J(rows(keep),1,1)
        s.eligible_cohorts=(rows(post) ? _cdid_agg_uniq(g[post]) : J(0,1,.))
        if (type=="simple") s.keys=J(1,1,.)
        else s.keys=s.eligible_cohorts
    }
    else {
        s.convention="released did 2.5.1: e = t - g in original date units; balance_e keeps cohorts with max(t) - g >= balance_e and horizons in [balance_e - max(t) + min(grid), balance_e]; then [min_e, max_e]"
        e_all=tt-g; max_t=max(tt); min_grid=min(s.grid)
        s.declared_event_times=_cdid_agg_uniq(e_all)
        if (balance_e<.) {
            s.include=(max_t:-g):>=balance_e
            eseq=_cdid_agg_uniq(_cdid_agg_subset(e_all,s.include))
            eseq=_cdid_agg_subset(eseq,(eseq:<=balance_e):&(eseq:>=balance_e-max_t+min_grid))
        }
        else eseq=s.declared_event_times
        if (min_e<.) eseq=_cdid_agg_subset(eseq,eseq:>=min_e)
        if (max_e<.) eseq=_cdid_agg_subset(eseq,eseq:<=max_e)
        if (rows(eseq)==0)
            _cdid_agg_fail("cdid_ipw_pipeline_selector_error","no event time falls within the window set by mine(), maxe() and balancee(); the group-time results are unchanged.",498)
        if (s.nondefault & !any(eseq:>=0))
            _cdid_agg_fail("cdid_ipw_pipeline_selector_error","No valid att_gt() estimates found for this aggregation: the explicit event window contains no nonnegative horizon (released did 2.5.1 selector refusal).",498)
        s.keys=eseq; s.event_window=eseq
        for (r=1;r<=R;r++) if (s.include[r]) {
            hit=_cdid_agg_lookup(e_all[r],eseq)
            if (hit[1]>0) s.selected[r]=1
        }
        s.eligible_cohorts=_cdid_agg_uniq(_cdid_agg_subset(g,s.include))
        s.excluded_cohorts=_cdid_agg_uniq(_cdid_agg_subset(g,!s.include))
    }
    s.M=rows(s.keys)
    s.sets=J(0,0,.)
    return(s)
}
/* Materialize after actual dimensions are known and budgeted by aggregate. */
void _cdid_agg_membership(struct cdid_agg_selection scalar s, real colvector g, real colvector tt)
{
    real scalar r
    real colvector key
    s.sets=J(s.M,s.R,0)
    for(r=1;r<=s.R;r++) if(s.selected[r]) {
        if(s.type=="simple") s.sets[1,r]=1
        else {
            key=_cdid_agg_lookup((s.type=="group" ? g[r] : (s.type=="calendar" ? tt[r] : tt[r]-g[r])),s.keys)
            s.sets[key[1],r]=1
        }
    }
}
/* Standalone selection preserves its original unbudgeted API. */
struct cdid_agg_selection scalar cdid_agg_select(string scalar type, real colvector g,
    real colvector tt, real scalar balance_e, real scalar min_e, real scalar max_e,
    real colvector calendar, real colvector cohort_grid)
{
    struct cdid_agg_selection scalar s
    s=_cdid_agg_selection_plan(type,g,tt,balance_e,min_e,max_e,calendar,cohort_grid)
    _cdid_agg_membership(s,g,tt)
    return(s)
}
/* ---------------------------------------------------------------------------
   Full-frame weighted cohort shares and their ratio IF (R l.4507-4557).
   Ordinary law: p_k = sum(w 1{g=k})/sum(w), P[,k] = w (1{g=k} - p_k)/mean(w).
   First-calendar masses with gaps: p_k = sum(u)/n, P[,k] = u - p_k with
   u = w on cohort k, 0 elsewhere; a cohort lacking V is flagged incomplete and
   refused only by aggregates that read its mass. Never-treated V is not demanded.
   No zero-times-missing product is formed.
   --------------------------------------------------------------------------- */
struct cdid_agg_shares scalar _cdid_agg_shares(struct cdid_agg_frame scalar F, real colvector cohorts)
{
    struct cdid_agg_shares scalar sh
    real colvector ind, gap, u, hit
    real scalar k, n, totalw, meanw, restricted
    string scalar ids
    n=F.n; sh.K=rows(cohorts); sh.cohorts=cohorts
    sh.p=J(sh.K,1,.); sh.P=J(n,sh.K,0); sh.mass_complete=J(sh.K,1,1); sh.mass_missing=J(sh.K,1,"")
    sh.first_masses=(F.share_law=="first_calendar_masses")
    sh.shares_known=all(F.cohort_known); sh.n_unknown=sum(!F.cohort_known)
    sh.incomplete_masses=sh.first_masses & missing(F.w)>0
    /* Restricted eligibility contexts (missingreference(exclude)): the mass of a
       restricted cohort counts only its eligible treated units,
       m_g = mean(V 1[G=g] E_g), with contribution V 1[G=g] E_g - m_g. */
    restricted=0
    if (rows(F.share_excluded)==n & n>0) restricted=any(F.share_excluded:!=0)
    if (restricted & !sh.first_masses)
        _cdid_agg_fail("cdid_ipw_pipeline_aggregation_error","Restricted eligibility masses need first-calendar masses.",498)
    if (sh.shares_known & (sh.incomplete_masses | restricted)) {
        for (k=1;k<=sh.K;k++) {
            ind=(F.cohort:==cohorts[k])
            if (restricted) ind=ind:&(F.share_excluded:==0)
            gap=_cdid_agg_index(ind:&(F.w:>=.))
            if (rows(gap)) {
                sh.mass_complete[k]=0
                ids=""
                if (rows(F.frame_id)==n) ids=" (ids "+_cdid_agg_join(F.frame_id[gap[|1\min((5,rows(gap)))|]],", ")+
                    (rows(gap)>5 ? ", ..." : "")+")"
                sh.mass_missing[k]="cohort "+_cdid_agg_fmt(cohorts[k])+" lacks w at the first calendar date for "+
                    strofreal(rows(gap))+" unit(s)"+ids
                continue
            }
            hit=_cdid_agg_index(ind); u=J(n,1,0)
            if (rows(hit)) u[hit]=F.w[hit]
            sh.p[k]=sum(u)/n; sh.P[,k]=u:-sh.p[k]
        }
    }
    else if (sh.shares_known) {
        totalw=sum(F.w); meanw=totalw/n
        for (k=1;k<=sh.K;k++) {
            ind=(F.cohort:==cohorts[k])
            sh.p[k]=sum(F.w:*ind)/totalw
            sh.P[,k]=F.w:*(ind:-sh.p[k])/meanw
        }
    }
    sh.share_reason="missing_frame_cohort_share: "+strofreal(sh.n_unknown)+
        " frame unit(s) have unknown cohort, so population shares of cohorts {"+
        _cdid_agg_join(_cdid_agg_fmtv(cohorts),",")+
        "} are not identified by the supplied frame; no share is inferred from eligible records"
    return(sh)
}
/* Full-frame IF rows stream into the declared-cluster covariance. */
struct cdid_if_summary scalar _cdid_agg_if_summary(real matrix Z, struct cdid_agg_frame scalar F,
    struct cdid_infer_budget scalar budget)
{
    struct cdid_if_state scalar st
    real scalar Jt, ir, lo, hi
    Jt=cols(Z)
    st=cdid_if_start(F.n,F.unit_cluster,F.C,Jt,F.clustered,budget)
    ir=min((F.n,budget.max_if_batch_rows,floor(budget.max_if_batch_cells/Jt)))
    if (ir<1) _cdid_agg_fail("cdid_ipw_pipeline_bound_error","IF batch budget admits no row.",908)
    for (lo=1;lo<=F.n;lo=lo+ir) {
        hi=min((lo+ir-1,F.n))
        cdid_if_add(st,(lo::hi),Z[|lo,1\hi,Jt|])
    }
    return(cdid_if_finish(st))
}
struct cdid_agg_inference scalar _cdid_agg_empty_inference(real scalar M,
    struct cdid_agg_frame scalar F, string scalar reason)
{
    struct cdid_agg_inference scalar out
    out.M=M; out.n_available=0; out.n_F=F.n; out.C=F.C
    out.status="empty_supported_selection"; out.inference_status="empty_supported_selection"
    out.reason=reason; out.flag_law="empty_selection_no_draws"; out.band_family="none"; out.warning=""
    out.bootstrap=0; out.simultaneous=0; out.draws_used=0; out.n_draws=.
    out.normal_critical=.; out.band_critical=.
    out.available=J(M,1,0)
    out.estimates=out.se=out.analytic_se=out.z_stat=out.p_value=J(M,1,.)
    out.ci_lower=out.ci_upper=out.band_lower=out.band_upper=J(M,1,.)
    out.covariance=out.Omega=J(M,M,.)
    return(out)
}
/* supported-subvector inference (R l.4257-4328): estimates and IF columns
   are subset to the available positions BEFORE any cluster sum, centering,
   scale inversion or multiplier; results scatter into missing-initialized
   positions. Signs are caller-supplied B x C; the kernel draws nothing. */
struct cdid_agg_inference scalar _cdid_agg_supported_inference(real colvector est,
    real matrix Psi, real colvector available, struct cdid_agg_frame scalar F, real scalar alpha,
    struct cdid_infer_flags scalar flags, real scalar B, real matrix signs,
    struct cdid_infer_budget scalar budget, string scalar reason)
{
    struct cdid_agg_inference scalar out
    struct cdid_if_summary scalar cov
    struct cdid_multiplier_state scalar ms
    struct cdid_multiplier_result scalar mult
    struct cdid_infer_result scalar sub
    real colvector av, positive
    real scalar M, j, lo, hi, br
    M=rows(est); av=_cdid_agg_index(available)
    out=_cdid_agg_empty_inference(M,F,reason); out.available=available
    if (rows(av)==0) return(out)
    for (j=1;j<=rows(av);j++) if (est[av[j]]>=. | missing(Psi[,av[j]]))
        _cdid_agg_fail("cdid_ipw_pipeline_integrity_error","An available target must have a finite estimate and finite influence-function column.",498)
    cov=_cdid_agg_if_summary(Psi[,av],F,budget)
    positive=_cdid_agg_index(cov.analytic_se:>0)
    mult=cdid_multiplier_none()
    if (flags.draws & rows(positive)) {
        if (B<1 | B>=. | B!=floor(B) | B>2^31-1)
            _cdid_agg_fail("cdid_ipw_pipeline_bootstrap_error","The resolved flag law needs draws: supply an explicit positive integer draw count.",498)
        if (rows(signs)!=B | cols(signs)!=F.C)
            _cdid_agg_fail("cdid_ipw_pipeline_bootstrap_error","The resolved flag law needs draws: supply an explicit B x C matrix of +/-1 signs; the aggregation kernel generates none.",498)
        ms=cdid_multiplier_start(cov,positive,B,alpha,budget)
        br=min((B,budget.max_sign_batch_rows,floor(budget.max_sign_batch_cells/F.C)))
        if (br<1) _cdid_agg_fail("cdid_ipw_pipeline_bound_error","Sign batch budget admits no draw row.",908)
        for (lo=1;lo<=B;lo=lo+br) {
            hi=min((lo+br-1,B))
            cdid_multiplier_add(ms,(lo::hi),(1..F.C),signs[|lo,1\hi,F.C|])
        }
        mult=cdid_multiplier_finish(ms,0)
    }
    sub=cdid_infer_report(est[av],cov,flags,mult,alpha)
    /* summary keeps the declared-cluster raw/centered sums (R cluster_if,
       raw_cluster_sums) for the supported columns, in av order. */
    out.compact=sub; out.summary=cov; out.n_available=rows(av)
    out.status=(rows(av)==M ? "all_requested_supported" : "supported_subvector")
    out.inference_status=sub.status; out.reason=sub.reason
    out.flag_law=sub.flag_law; out.band_family=sub.band_family; out.warning=sub.warning
    out.bootstrap=sub.bootstrap; out.simultaneous=sub.simultaneous
    out.draws_used=sub.draws_used; out.n_draws=sub.n_draws
    out.normal_critical=sub.normal_critical; out.band_critical=sub.band_critical
    out.estimates[av]=sub.estimates; out.se[av]=sub.se; out.analytic_se[av]=sub.analytic_se
    out.z_stat[av]=sub.z_stat; out.p_value[av]=sub.p_value
    out.ci_lower[av]=sub.ci_lower; out.ci_upper[av]=sub.ci_upper
    out.band_lower[av]=sub.band_lower; out.band_upper[av]=sub.band_upper
    out.covariance[av,av]=sub.covariance; out.Omega[av,av]=sub.Omega
    return(out)
}
void _cdid_agg_check_inputs(real colvector theta, real matrix IF, real colvector available,
    struct cdid_agg_ledger scalar L, struct cdid_agg_frame scalar F, struct cdid_agg_request scalar Q)
{
    real scalar R, n, j
    real colvector av
    R=rows(theta); n=F.n
    if (Q.type!="simple" & Q.type!="group" & Q.type!="dynamic" & Q.type!="calendar")
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","type must be simple, group, dynamic or calendar.",498)
    if (Q.na_rm!=0 & Q.na_rm!=1) _cdid_agg_fail("cdid_ipw_pipeline_input_error","na.rm must be a scalar logical (0/1).",498)
    if (Q.alpha<=0 | Q.alpha>=1 | Q.alpha>=.) _cdid_agg_fail("cdid_ipw_pipeline_input_error","alpha must lie strictly inside (0,1).",498)
    if (n<1 | n>=. | n!=floor(n) | R<1 | cols(theta)!=1 | rows(IF)!=n | cols(IF)!=R |
        rows(available)!=R | cols(available)!=1 | missing(available) | any((available:!=0):&(available:!=1)))
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","available must be a complete 0/1 selector over the R cells and IF must be n x R.",498)
    if (rows(L.g)!=R | rows(L.time)!=R | rows(L.ref)!=R | rows(L.base_t)!=R | rows(L.target_id)!=R |
        rows(L.role)!=R | rows(L.support_status)!=R | missing(L.g) | missing(L.time))
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","Reference metadata must cover the declared target rows before support selection.",498)
    av=_cdid_agg_index(available)
    for (j=1;j<=rows(av);j++) if (theta[av[j]]>=. | missing(IF[,av[j]]))
        _cdid_agg_fail("cdid_ipw_pipeline_input_error","available cells need finite estimates and IF columns.",498)
    if (rows(F.cohort)!=n | rows(F.cohort_known)!=n | rows(F.w)!=n | rows(F.unit_cluster)!=n |
        missing(F.cohort_known) | any((F.cohort_known:!=0):&(F.cohort_known:!=1)) |
        any((F.cohort:<.):&(!F.cohort_known)) | (rows(F.frame_id)!=0 & rows(F.frame_id)!=n) |
        (F.share_law!="" & F.share_law!="first_calendar_masses"))
        _cdid_agg_fail("cdid_ipw_pipeline_aggregation_error","Aggregation needs same-frame cohort labels (unknown only for genuinely unknown frame-only cohorts), a 0/1 known flag and a share law of empty or first_calendar_masses.",498)
    if (F.share_law=="first_calendar_masses") {
        if (any((F.w:<.):&(F.w:<=0)))
            _cdid_agg_fail("cdid_ipw_pipeline_aggregation_error","Present first-calendar masses must be positive.",498)
    }
    else if (missing(F.w) | any(F.w:<=0))
        _cdid_agg_fail("cdid_ipw_pipeline_aggregation_error","Aggregation needs positive fixed weights on every frame unit.",498)
}
/* Serial whole-buffer envelope of this kernel's own allocations, measured
   against the caller's explicit budget; the inference kernels
   preflight their own summaries, draws and batches. No default cap. */
void _cdid_agg_preflight(struct cdid_infer_budget scalar b, real scalar n, real scalar R,
    real scalar M, real scalar K, real scalar na, real scalar plan_entries)
{
    real scalar cap, total
    cap=min((b.max_workspace_cells,floor(b.max_workspace_bytes/8)))
    total=plan_entries
    total=_cdid_infer_cells(3*n,M,cap,total)
    total=_cdid_infer_cells(n,K,cap,total)
    total=_cdid_infer_cells(n,na+K,cap,total)
    total=_cdid_infer_cells(2*n,1,cap,total)
    total=_cdid_infer_cells(5*M,R+K,cap,total)
    total=_cdid_infer_cells(12*M,M,cap,total)
    total=_cdid_infer_cells(na+K,na+K,cap,total)
}
/* ---------------------------------------------------------------------------
   Population cohort-share aggregation (R aggregate_chained_ipw_targets).
   theta: R x 1 cell estimates (missing allowed where unavailable);
   IF: n x R full-F target IF (columns of unavailable cells never read);
   available: R x 1 explicit selector from the reporting map;
   signs: B x C +/-1 matrix when the resolved law draws, else may be empty.
   --------------------------------------------------------------------------- */
struct cdid_agg_result scalar cdid_aggregate(real colvector theta, real matrix IF,
    real colvector available, struct cdid_agg_ledger scalar L, struct cdid_agg_frame scalar F,
    struct cdid_agg_request scalar Q, real matrix signs, struct cdid_infer_budget scalar budget)
{
    struct cdid_agg_result scalar out
    struct cdid_infer_flags scalar flags, fixed_flags, overall_flags
    struct cdid_if_summary scalar joint
    real scalar n, R, M, K, m, j, na, nam, den, S_g, reads_mass, reuse, n_coh, ok, plan_entries, plan_cap
    real colvector post, cohorts, rr, miss, used, kk, w, av, am, known_ref, outer, outer_missing, surv
    string colvector reasons
    string scalar bin_endpoint, cell_reason, prefix
    real matrix Ott, Otp, Opp, Wa, Da, Of, Os, Oc
    _cdid_agg_check_inputs(theta,IF,available,L,F,Q)
    n=F.n; R=rows(theta)
    out.type=Q.type; out.na_rm=Q.na_rm; out.n=n; out.R=R
    /* Flag law is validated before any selection (R l.4466-4470). */
    bin_endpoint=(Q.type=="simple" ? "simple" : "aggte_bins")
    flags=cdid_infer_resolve_flags(Q.bstrap,Q.cband,bin_endpoint)
    out.flag_law=flags.law
    // Linear key-planning numeric payload before metadata/sort temporaries.
    // No R-squared substitution. Excludes allocator/string/BLAS internals/RSS.
    plan_cap=min((budget.max_workspace_cells,floor(budget.max_workspace_bytes/8)))
    plan_entries=_cdid_infer_cells(12,R,plan_cap,0)
    plan_entries=_cdid_infer_cells(6,rows(Q.calendar)+rows(Q.cohort_grid)+2*R,plan_cap,plan_entries)
    post=_cdid_agg_index(L.time:>=L.g)
    cohorts=(rows(post) ? _cdid_agg_uniq(L.g[post]) : J(0,1,.)); K=rows(cohorts)
    out.K=K
    out.selectors=_cdid_agg_selection_plan(Q.type,L.g,L.time,Q.balance_e,Q.min_e,Q.max_e,Q.calendar,Q.cohort_grid)
    M=out.selectors.M; out.M=M; out.keys=out.selectors.keys
    av=_cdid_agg_index(available); na=rows(av)
    _cdid_agg_preflight(budget,n,R,M,K,(K ? na : 0),plan_entries)
    _cdid_agg_membership(out.selectors,L.g,L.time)
    prefix=(Q.type=="group" ? "group_" : (Q.type=="calendar" ? "calendar_" : "event_"))
    out.labels=J(M,1,"ATT_simple")
    if (Q.type!="simple") for (m=1;m<=M;m++) out.labels[m]=prefix+_cdid_agg_fmt(out.keys[m])
    out.intended=out.selectors.sets; out.used=J(M,R,0); out.excluded=J(M,R,0)
    out.n_intended=(M ? rowsum(out.intended) : J(0,1,.)); out.n_used=J(M,1,0)
    out.cell_available=J(M,1,0); out.restricted=J(M,1,0); out.share_dependent=J(M,1,0)
    out.share_denominator=J(M,1,.); out.estimates=J(M,1,.)
    out.W=J(M,R,0); out.D=J(M,K,0)
    out.if_total=out.if_fixed=out.if_share=J(n,M,.)
    out.bin_role=J(M,1,"aggregate"); out.bin_support_status=J(M,1,"aggregate")
    reasons=J(M,1,"")
    out.released_warning=""
    out.shares=_cdid_agg_shares(F,cohorts)
    out.overall.estimate=.; out.overall.se=.; out.overall.ci_lower=.; out.overall.ci_upper=.
    out.overall.if_full=J(n,1,.); out.overall.weights=J(M,1,0); out.overall.share_gradient=J(K,1,0)
    out.overall.share_dependent=0; out.overall.share_denominator=.
    out.overall.intended=out.overall.used=out.overall.excluded=J(0,1,.)
    out.overall.inference_kind="analytic"; out.overall.draws_reused_from=""
    out.decomposition.computed=0; out.decomposition.shares_known=out.shares.shares_known
    out.decomposition.Omega_total=out.decomposition.Omega_fixed=J(M,M,.)
    out.decomposition.Omega_share=out.decomposition.Omega_cross=J(M,M,.)
    out.decomposition.Omega_theta_p=J(0,0,.)
    if (K==0) {
        cell_reason="no_declared_post_cells: the reporting ledger declares no post cell"
        out.cell_reason=J(M,1,cell_reason); out.status="unavailable"
        out.inference=_cdid_agg_empty_inference(M,F,cell_reason)
        out.fixed_inference=_cdid_agg_empty_inference(M,F,cell_reason)
        out.overall.status="unavailable"; out.overall.reason=cell_reason
        out.overall.inference=_cdid_agg_empty_inference(1,F,cell_reason)
        out.target_definition=""
        return(out)
    }
    /* Authoritative universal-reference rows from source metadata only; a
       zero point or cancelled IF never establishes this role (R l.3953-3959). */
    known_ref=(L.role:=="reference"):&(L.support_status:=="deterministic_reference_zero"):&
        (L.time:==L.ref):&(L.time:==L.base_t)
    /* Per-aggregate intended sets, restricted sets and local weights (R l.4566-4618). */
    for (m=1;m<=M;m++) {
        rr=_cdid_agg_index(out.intended[m,]')
        miss=_cdid_agg_subset(rr,!available[rr])
        if (rows(miss)) out.excluded[m,miss]=J(1,rows(miss),1)
        if (rows(rr)==0) { reasons[m]="empty_intended_set"; continue; }
        if (!Q.na_rm & rows(miss)) {
            reasons[m]="missing_targets (na.rm = FALSE): "+
                _cdid_agg_join(L.target_id[miss]:+" [":+L.support_status[miss]:+"]","; ")
            continue
        }
        used=_cdid_agg_subset(rr,available[rr])
        if (rows(used)==0) {
            reasons[m]="empty_restricted_set (na.rm = TRUE): every intended cell is unavailable: "+_cdid_agg_join(L.target_id[miss],"; ")
            continue
        }
        out.restricted[m]=rows(miss)>0
        out.used[m,used]=J(1,rows(used),1); out.n_used[m]=rows(used)
        kk=_cdid_agg_lookup(L.g[used],cohorts)
        if (any(kk:==0)) _cdid_agg_fail("cdid_ipw_pipeline_integrity_error","A used cell's cohort has no post cell in the declared table of cells; the cohort share set is formed from post cells only.",498)
        n_coh=rows(uniqrows(sort(kk,1)))
        out.share_dependent[m]=(Q.type!="group") & n_coh>=2
        if (out.share_dependent[m] & !out.shares.shares_known) { reasons[m]=out.shares.share_reason; continue; }
        if (out.share_dependent[m] & !all(out.shares.mass_complete[kk])) {
            reasons[m]="missing_first_period_weight (D61): "+
                _cdid_agg_join(_cdid_agg_uniqstr(_cdid_agg_strsubset(out.shares.mass_missing[kk],!out.shares.mass_complete[kk])),"; ")+
                "; first-calendar masses are never replaced by first-observed, mean or anchor weights"
            continue
        }
        if (Q.type!="group" & out.shares.shares_known & (!out.shares.incomplete_masses | out.share_dependent[m])) {
            out.share_denominator[m]=sum(out.shares.p[kk])
            if (out.share_denominator[m]>=. | out.share_denominator[m]<=0) { reasons[m]="nonpositive_share_mass"; continue; }
        }
        if (Q.type=="group" | n_coh==1) {
            /* Fixed uniform weights; the one-cohort ratio map is degree-0 with exact zero share derivative. */
            w=J(rows(used),1,1/rows(used)); den=.
        }
        else { den=out.share_denominator[m]; w=out.shares.p[kk]:/den; }
        out.estimates[m]=sum(w:*theta[used]); out.W[m,used]=w'; out.cell_available[m]=1
        if (out.share_dependent[m]) for (j=1;j<=rows(used);j++)
            out.D[m,kk[j]]=out.D[m,kk[j]]+(theta[used[j]]-out.estimates[m])/den
        if (all(known_ref[used])) { out.bin_role[m]="reference_aggregate"; out.bin_support_status[m]="deterministic_reference_zero"; }
    }
    am=_cdid_agg_index(out.cell_available); nam=rows(am); na=rows(av)
    /* Local products on supported columns only (R l.4621-4627). */
    if (nam) {
        out.if_fixed[,am]=IF[,av]*out.W[am,av]'
        out.if_share[,am]=(out.shares.shares_known ? out.shares.P*out.D[am,]' : J(n,nam,0))
        out.if_total[,am]=out.if_fixed[,am]+out.if_share[,am]
    }
    cell_reason=(nam<M ? "unavailable aggregates: "+_cdid_agg_join(_cdid_agg_strsubset(out.labels,!out.cell_available),", ") :
        "all intended aggregates available")
    out.inference=_cdid_agg_supported_inference(out.estimates,out.if_total,out.cell_available,F,Q.alpha,flags,Q.B,signs,budget,cell_reason)
    /* Released compute.aggte warning exactly when band draws ran under bstrap=FALSE, cband=TRUE. */
    if (out.inference.draws_used & out.inference.simultaneous & !out.inference.bootstrap &
        out.inference.flag_law=="released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band")
        out.released_warning=out.inference.warning
    /* Fixed-share-only SEs: analytic diagnostic, never draws (R l.4682-4685). */
    fixed_flags=cdid_infer_resolve_flags(0,0,"cells")
    out.fixed_inference=_cdid_agg_supported_inference(out.estimates,out.if_fixed,out.cell_available,F,Q.alpha,fixed_flags,Q.B,J(0,0,.),budget,cell_reason)
    /* Cluster-aware fixed/share/cross decomposition through the
       declared-cluster operator on the joint (theta IF, share IF) stream (R l.4706-4730). */
    if (nam) {
        joint=_cdid_agg_if_summary((out.shares.shares_known ? (IF[,av],out.shares.P) : IF[,av]),F,budget)
        Ott=joint.Omega[|1,1\na,na|]; Wa=out.W[am,av]; Of=Wa*Ott*Wa'
        if (out.shares.shares_known) {
            Otp=joint.Omega[|1,na+1\na,na+K|]; Opp=joint.Omega[|na+1,na+1\na+K,na+K|]
            Da=out.D[am,]; Os=Da*Opp*Da'; Oc=Wa*Otp*Da'+Da*Otp'*Wa'
            out.decomposition.Omega_theta_p=J(R,K,.); out.decomposition.Omega_theta_p[av,]=Otp
        }
        else { Os=Oc=J(nam,nam,0); }
        out.decomposition.Omega_total[am,am]=Of+Os+Oc; out.decomposition.Omega_fixed[am,am]=Of
        out.decomposition.Omega_share[am,am]=Os; out.decomposition.Omega_cross[am,am]=Oc
        out.decomposition.computed=1
    }
    /* Outer summary: intended bins (na.rm=FALSE) or surviving bins renormalized (na.rm=TRUE). */
    outer=(Q.type=="dynamic" ? _cdid_agg_index(out.keys:>=0) : (1::M))
    outer_missing=_cdid_agg_subset(outer,!out.cell_available[outer])
    surv=(Q.na_rm ? _cdid_agg_subset(outer,out.cell_available[outer]) : (rows(outer_missing) ? J(0,1,.) : outer))
    out.overall.status="available"; out.overall.reason=""
    if (rows(outer)==0) { out.overall.status="unavailable"; out.overall.reason="no_outer_bins: the declared grid has no nonnegative outer bin"; }
    else if (!Q.na_rm & rows(outer_missing)) {
        out.overall.status="unavailable"
        out.overall.reason="missing_inner_aggregates (na.rm = FALSE): "+_cdid_agg_join(out.labels[outer_missing],", ")
    }
    else if (rows(surv)==0) {
        out.overall.status="unavailable"
        out.overall.reason="empty_restricted_outer_set (na.rm = TRUE): every intended bin is unavailable: "+_cdid_agg_join(out.labels[outer],", ")
    }
    else if (Q.type=="group" & rows(surv)>=2) {
        // Mata & is not R &&: access cohort masses only inside the group branch.
        if (!out.shares.shares_known) {
            out.overall.status="unavailable"; out.overall.reason=out.shares.share_reason
        }
        else if (!all(out.shares.mass_complete[surv])) {
            out.overall.status="unavailable"
            out.overall.reason="missing_first_period_weight (D61): "+
                _cdid_agg_join(_cdid_agg_strsubset(out.shares.mass_missing[surv],!out.shares.mass_complete[surv]),"; ")
        }
    }
    /* A one-bin group overall reads no mass (exact degree-0 identity). */
    reads_mass=(!out.shares.incomplete_masses | rows(surv)>=2)
    if (Q.type=="group" & out.shares.shares_known & rows(surv) & reads_mass) {
        // R sum propagates absent masses; Mata sum would omit missing values.
        if (missing(out.shares.p[surv])) out.overall.share_denominator=.
        else out.overall.share_denominator=sum(out.shares.p[surv])
    }
    if (out.overall.status=="available" & Q.type=="group" & out.shares.shares_known & reads_mass &
        (out.overall.share_denominator>=. | out.overall.share_denominator<=0)) {
        out.overall.status="unavailable"; out.overall.reason="nonpositive_share_mass"
    }
    out.overall.intended=outer; out.overall.used=surv
    out.overall.excluded=_cdid_agg_subset(outer,!_cdid_agg_member(outer,surv))
    if (out.overall.status=="available") {
        if (Q.type=="group") {
            if (rows(surv)>=2) {
                S_g=sum(out.shares.p[surv]); out.overall.weights[surv]=out.shares.p[surv]/S_g
                out.overall.estimate=sum(out.overall.weights[surv]:*out.estimates[surv])
                out.overall.share_gradient[surv]=(out.estimates[surv]:-out.overall.estimate)/S_g
                out.overall.share_dependent=1
            }
            else { out.overall.weights[surv]=1; out.overall.estimate=out.estimates[surv[1]]; }
        }
        else {
            out.overall.weights[surv]=J(rows(surv),1,1/rows(surv))
            out.overall.estimate=sum(out.overall.weights[surv]:*out.estimates[surv])
        }
        out.overall.if_full=out.if_total[,surv]*out.overall.weights[surv]
        if (out.overall.share_dependent) out.overall.if_full=out.overall.if_full+out.shares.P*out.overall.share_gradient
    }
    /* Overall: pointwise by upstream convention; multiplier IQR SE only when
       bstrap is explicitly TRUE, analytic otherwise (R l.4780-4787). */
    reuse=(Q.type=="simple") & !out.inference.simultaneous
    if (reuse) out.overall.inference=out.inference
    else {
        overall_flags=cdid_infer_resolve_flags((Q.bstrap>=. ? 0 : Q.bstrap),0,"cells")
        out.overall.inference=_cdid_agg_supported_inference(out.overall.estimate,out.overall.if_full,
            (out.overall.status=="available"),F,Q.alpha,overall_flags,Q.B,signs,budget,
            (out.overall.reason=="" ? "overall unavailable" : out.overall.reason))
    }
    out.overall.se=out.overall.inference.se[1]; out.overall.ci_lower=out.overall.inference.ci_lower[1]
    out.overall.ci_upper=out.overall.inference.ci_upper[1]
    out.overall.inference_kind=(out.overall.inference.bootstrap ? "multiplier_bootstrap_iqr_se" : "analytic")
    out.overall.draws_reused_from=(reuse & out.inference.draws_used ? "inference" : "")
    out.cell_reason=J(M,1,"")
    for (m=1;m<=M;m++) {
        if (!out.cell_available[m]) out.cell_reason[m]=reasons[m]
        else if (Q.na_rm & out.restricted[m]) out.cell_reason[m]="restricted_target (na.rm = TRUE): weights renormalized over the supported subset"
        else out.cell_reason[m]="intended_target"
    }
    out.target_definition=(Q.na_rm ?
        "explicit restricted targets: each aggregate averages its SUPPORTED intended cells with renormalized weights; outer summaries omit empty bins and renormalize over surviving bins (did 2.5.1 na.rm = TRUE convention); the intended sets, exclusions and changed targets are displayed in cell_status and overall" :
        "intended targets: each aggregate requires every intended cell; an aggregate with an unavailable intended cell is NA with the missing IDs and reasons; unrelated aggregates survive")
    ok=(nam==M & out.overall.status=="available")
    out.status=(ok ? "available" : (nam ? "partially_available" : "unavailable"))
    return(out)
}
string colvector _cdid_agg_strsubset(string colvector s, real colvector mask)
{
    real colvector ix
    ix=_cdid_agg_index(mask)
    if (rows(ix)==0) return(J(0,1,""))
    return(s[ix])
}
real colvector _cdid_agg_member(real colvector x, real colvector set)
{
    real colvector out
    real scalar i
    out=J(rows(x),1,0)
    for (i=1;i<=rows(x);i++) out[i]=anyof(set,x[i])
    return(out)
}
end
