*! Native full-frame inference kernels, current StataNow 19.5.
*! Pinned R cdid 1.1.0.0: chained_multiplier.R, chained_pretest.R, pipeline flags.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_infer(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_infer_budget {
    real scalar max_units, max_clusters, max_targets, max_draws
    real scalar max_if_batch_rows, max_if_batch_cells
    real scalar max_sign_batch_rows, max_sign_batch_cells
    real scalar max_workspace_cells, max_workspace_bytes
}
struct cdid_if_state {
    struct cdid_infer_budget scalar budget
    real scalar n_F, C, J, clustered, n_seen, finished
    real colvector unit_cluster, seen, cluster_count
    real rowvector absolute_if_sum
    real matrix raw_sums
}
struct cdid_if_summary {
    real scalar n_F, C, J, clustered
    real colvector cluster_count, analytic_se
    real rowvector absolute_if_sum, cluster_mean
    real matrix raw_sums, centered_sums, covariance, Omega
}
struct cdid_multiplier_state {
    struct cdid_infer_budget scalar budget
    real scalar n_F, C, J, source_J, B, alpha, cursor, finished
    real colvector target_positions, analytic_variance
    real matrix centered_sums, raw_draws
}
struct cdid_multiplier_result {
    real scalar n_F, C, J, source_J, B, alpha, accepted, critical_value
    real scalar q25_index, q75_index, critical_index, iqr_normalizer
    real colvector target_positions, analytic_se, q25, q75, raw_iqr, raw_scale
    real colvector candidate_se, bootstrap_scale, bootstrap_se, sup_t, target_accepted
    real matrix raw_draws
    string scalar status
    string colvector analytic_status, bootstrap_status, target_status
}
struct cdid_infer_flags {
    real scalar bootstrap, simultaneous, draws, released_false_true
    string scalar law, endpoint
}
struct cdid_infer_result {
    real scalar n_F, C, J, bootstrap, simultaneous, draws_used, n_draws
    real scalar normal_critical, band_critical
    real colvector estimates, se, analytic_se, z_stat, p_value
    real colvector ci_lower, ci_upper, band_lower, band_upper
    real matrix covariance, Omega
    string scalar status, reason, flag_law, band_family, warning
}
struct cdid_wald_result {
    real scalar df, statistic, p_value, rcond_upper_1, covariance_rank_upper
    real colvector solve_residual
    string scalar status, reason
}
void _cdid_infer_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_infer_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
struct cdid_infer_budget scalar cdid_infer_budget_create(real scalar max_units,
    real scalar max_clusters, real scalar max_targets, real scalar max_draws,
    real scalar max_if_batch_rows, real scalar max_if_batch_cells,
    real scalar max_sign_batch_rows, real scalar max_sign_batch_cells,
    real scalar max_workspace_cells, real scalar max_workspace_bytes)
{
    struct cdid_infer_budget scalar b
    real rowvector values
    values=(max_units,max_clusters,max_targets,max_draws,max_if_batch_rows,
        max_if_batch_cells,max_sign_batch_rows,max_sign_batch_cells,
        max_workspace_cells,max_workspace_bytes)
    if(missing(values) | any(values:<1) | any(values:!=floor(values)) | any(values:>2^53-1))
        _cdid_infer_fail("Every inference budget is an explicit positive finite integer <=2^53-1.",498)
    b.max_units=max_units; b.max_clusters=max_clusters; b.max_targets=max_targets; b.max_draws=max_draws
    b.max_if_batch_rows=max_if_batch_rows; b.max_if_batch_cells=max_if_batch_cells
    b.max_sign_batch_rows=max_sign_batch_rows; b.max_sign_batch_cells=max_sign_batch_cells
    b.max_workspace_cells=max_workspace_cells; b.max_workspace_bytes=max_workspace_bytes
    return(b)
}
real scalar _cdid_infer_cells(real scalar a, real scalar b, real scalar cap, real scalar total)
{
    if(a>0 & b>(2^53-1)/a)
        _cdid_infer_fail("Inference allocation product exceeds exact integer range.",908)
    if(a*b>cap-total)
        _cdid_infer_fail("Inference storage exceeds the explicit workspace budget.",908)
    return(total+a*b)
}
/* Serial whole-buffer envelope; numeric payload, not peak RSS or CPU bound.
   Includes retained summaries, BxJ draws/product scratch and one sign/IF batch.
   No tight default work cap and no N-by-B sign allocation. */
real scalar _cdid_infer_preflight(struct cdid_infer_budget scalar b,
    real scalar n, real scalar C, real scalar J, real scalar B)
{
    real scalar cap, total, ir, br
    if(n>b.max_units | C>b.max_clusters | J>b.max_targets | B>b.max_draws)
        _cdid_infer_fail("Inference dimensions exceed explicit caller limits.",908)
    cap=min((b.max_workspace_cells,floor(b.max_workspace_bytes/8)))
    total=0; ir=min((n,b.max_if_batch_rows)); br=min((B,b.max_sign_batch_rows))
    total=_cdid_infer_cells(4*C,J,cap,total)
    total=_cdid_infer_cells(8*J,J,cap,total)
    total=_cdid_infer_cells(2*B,J,cap,total)
    total=_cdid_infer_cells(6,n,cap,total)
    total=_cdid_infer_cells(4,C,cap,total)
    total=_cdid_infer_cells(8,B,cap,total)
    total=_cdid_infer_cells(40,J,cap,total)
    total=_cdid_infer_cells(2*ir,J,cap,total)
    total=_cdid_infer_cells(4*br,C,cap,total)
    return(total)
}
/* Unit IDs are original-F integer positions, with a reversible external ID
   dictionary owned by population. Cluster columns use declared positions1..C.
   Every F unit is required once, including zero IF/frame-only units. */
struct cdid_if_state scalar cdid_if_start(real scalar n_F, real colvector unit_cluster,
    real scalar C, real scalar Jt, real scalar clustered, struct cdid_infer_budget scalar budget)
{
    struct cdid_if_state scalar s
    real scalar entries
    if(n_F<1 | n_F>=. | n_F!=floor(n_F) | n_F>2^31-1 |
        C<1 | C>=. | C!=floor(C) | C>2^31-1 | Jt<1 | Jt>=. | Jt!=floor(Jt) | Jt>2^31-1 |
        (clustered!=0 & clustered!=1) | rows(unit_cluster)!=n_F | cols(unit_cluster)!=1 |
        missing(unit_cluster) | any(unit_cluster:<1) | any(unit_cluster:>C) |
        any(unit_cluster:!=floor(unit_cluster)))
        _cdid_infer_fail("Require full-F cluster positions and finite positive inference dimensions.",498)
    if(clustered & C<2) _cdid_infer_fail("At least two declared clusters are required for clustered inference.",459)
    if(!clustered) if(C!=n_F | any(unit_cluster:!=(1..n_F)'))
        _cdid_infer_fail("Unclustered inference uses one declared cluster for each original unit.",498)
    entries=_cdid_infer_preflight(budget,n_F,C,Jt,0)
    s.budget=budget; s.n_F=n_F; s.C=C; s.J=Jt; s.clustered=clustered
    s.n_seen=0; s.finished=0; s.unit_cluster=unit_cluster
    s.seen=J(n_F,1,0); s.cluster_count=J(C,1,0)
    s.absolute_if_sum=J(1,Jt,0); s.raw_sums=J(C,Jt,0)
    return(s)
}
void cdid_if_add(struct cdid_if_state scalar s, real colvector unit_positions, real matrix IF)
{
    real scalar b, sorted_c
    real colvector cl_b, ord_c, cl_s, u_cl
    real matrix sp_c
    b=rows(IF)
    if(s.finished | b<1 | cols(IF)!=s.J | rows(unit_positions)!=b | cols(unit_positions)!=1 |
        b>s.n_F-s.n_seen | missing(unit_positions) | any(unit_positions:!=floor(unit_positions)) |
        any(unit_positions:<1) | any(unit_positions:>s.n_F))
        _cdid_infer_fail("Invalid or completed full-unit IF stream.",498)
    if(b>s.budget.max_if_batch_rows | b>floor(s.budget.max_if_batch_cells/s.J))
        _cdid_infer_fail("IF batch exceeds explicit row/cell budget.",908)
    if(b>1) {
        if(!all(unit_positions[|2\b|]:>unit_positions[|1\b-1|])) {
            if(rows(uniqrows(sort(unit_positions,1)))!=b)
                _cdid_infer_fail("A full-F unit IF row was supplied more than once.",498)
        }
    }
    if(any(s.seen[unit_positions]))
        _cdid_infer_fail("A full-F unit IF row was supplied more than once.",498)
    if(missing(IF)) _cdid_infer_fail("Complete IF rows must be finite and already weighted.",498)
    s.absolute_if_sum=s.absolute_if_sum+colsum(abs(IF))
    if(!s.clustered) {
        s.raw_sums[unit_positions,]=s.raw_sums[unit_positions,]+IF
        s.cluster_count[unit_positions]=s.cluster_count[unit_positions]:+1
    }
    else {
        cl_b=s.unit_cluster[unit_positions]
        sorted_c=(b==1 ? 1 : all(cl_b[|2\b|]:>=cl_b[|1\b-1|]))
        if(sorted_c) {
            sp_c=_cdid_ps_infer(cl_b,1)
            u_cl=cl_b[sp_c[,1]]
            s.raw_sums[u_cl,]=s.raw_sums[u_cl,]+panelsum(IF,sp_c)
        }
        else {
            ord_c=order((cl_b,(1::b)),(1,2))
            cl_s=cl_b[ord_c]
            sp_c=_cdid_ps_infer(cl_s,1)
            u_cl=cl_s[sp_c[,1]]
            s.raw_sums[u_cl,]=s.raw_sums[u_cl,]+panelsum(IF[ord_c,],sp_c)
        }
        s.cluster_count[u_cl]=s.cluster_count[u_cl]+(sp_c[,2]:-sp_c[,1]:+1)
    }
    s.seen[unit_positions]=J(b,1,1); s.n_seen=s.n_seen+b
}
struct cdid_if_summary scalar cdid_if_finish(struct cdid_if_state scalar s)
{
    struct cdid_if_summary scalar out
    real matrix cross
    if(s.finished | s.n_seen!=s.n_F | !all(s.seen))
        _cdid_infer_fail("Every original-F unit must occur exactly once before IF finalization.",498)
    out.n_F=s.n_F; out.C=s.C; out.J=s.J; out.clustered=s.clustered
    out.cluster_count=s.cluster_count; out.absolute_if_sum=s.absolute_if_sum
    out.raw_sums=s.raw_sums; out.cluster_mean=colsum(s.raw_sums)/s.C
    /* Center across DECLARED cluster sums, never before unit aggregation. */
    out.centered_sums=s.raw_sums:-out.cluster_mean
    cross=quadcross(out.centered_sums,out.centered_sums)
    cross=.5:*(cross+cross')
    out.covariance=cross/(s.n_F*s.n_F); out.Omega=cross/s.n_F
    out.analytic_se=sqrt(rowmax((diagonal(out.covariance),J(s.J,1,0))))
    if(missing(out.raw_sums) | missing(out.centered_sums) | missing(out.covariance) |
        missing(out.Omega) | missing(out.absolute_if_sum))
        _cdid_infer_fail("Nonfinite full-frame covariance or declared cluster sums.",504)
    s.finished=1
    return(out)
}
real scalar cdid_type1_index(real scalar B, real scalar probability)
{
    if(B<1 | B>=. | B!=floor(B) | B>2^31-1 | probability<0 | probability>1 | probability>=.)
        _cdid_infer_fail("Type1 quantiles require a positive draw count and probability in [0,1].",498)
    return(max((1,min((B,ceil(B*probability))))))
}
real scalar cdid_type1_quantile(real colvector values, real scalar probability)
{
    real colvector ordered
    real scalar k
    if(missing(values) | cols(values)!=1) _cdid_infer_fail("Type1 quantile values must be a finite column.",498)
    k=cdid_type1_index(rows(values),probability)
    ordered=sort(values,1)
    return(ordered[k])
}
/* target_positions fixes the complete declared multiplier family. The report
   helper additionally checks that it is exactly the positive-scale target set,
   so a refused member cannot silently be removed from the requested family. */
struct cdid_multiplier_state scalar cdid_multiplier_start(struct cdid_if_summary scalar cov,
    real colvector target_positions, real scalar B, real scalar alpha, struct cdid_infer_budget scalar budget)
{
    struct cdid_multiplier_state scalar s
    real scalar entries, Jm
    Jm=rows(target_positions)
    if(Jm<1 | cols(target_positions)!=1 | missing(target_positions) |
        any(target_positions:!=floor(target_positions)) | any(target_positions:<1) |
        any(target_positions:>cov.J) | rows(uniqrows(sort(target_positions,1)))!=Jm |
        B<1 | B>=. | B!=floor(B) | B>2^31-1 | alpha<=0 | alpha>=1 | alpha>=.)
        _cdid_infer_fail("Invalid declared multiplier family, draw count or alpha.",498)
    entries=_cdid_infer_preflight(budget,cov.n_F,cov.C,cov.J,B)
    s.budget=budget; s.n_F=cov.n_F; s.C=cov.C; s.J=Jm; s.source_J=cov.J
    s.B=B; s.alpha=alpha; s.cursor=0; s.finished=0
    s.target_positions=target_positions; s.analytic_variance=diagonal(cov.covariance)[target_positions]
    s.centered_sums=cov.centered_sums[,target_positions]; s.raw_draws=J(B,Jm,.)
    return(s)
}
void cdid_multiplier_add(struct cdid_multiplier_state scalar s,
    real colvector draw_positions, real rowvector cluster_positions, real matrix signs)
{
    real scalar b
    real colvector expected
    b=rows(signs)
    if(s.finished | b<1 | cols(signs)!=s.C | rows(draw_positions)!=b | cols(draw_positions)!=1 |
        b>s.B-s.cursor | cols(cluster_positions)!=s.C | rows(cluster_positions)!=1 |
        missing(draw_positions) | missing(cluster_positions))
        _cdid_infer_fail("Invalid or completed multiplier sign stream.",498)
    if(b>s.budget.max_sign_batch_rows | b>floor(s.budget.max_sign_batch_cells/s.C))
        _cdid_infer_fail("Sign batch exceeds explicit row/cell budget.",908)
    expected=(s.cursor+1..s.cursor+b)'
    if(any(draw_positions:!=expected) | any(cluster_positions:!=(1..s.C)))
        _cdid_infer_fail("Supplied signs must preserve exact declared draw and cluster order.",498)
    if(missing(signs) | any((signs:!=1):&(signs:!=-1)))
        _cdid_infer_fail("Multiplier signs must be exactly -1 or+1.",498)
    s.raw_draws[expected,]=(signs*s.centered_sums)/sqrt(s.C)
    s.cursor=s.cursor+b
}
struct cdid_multiplier_result scalar cdid_multiplier_none()
{
    struct cdid_multiplier_result scalar out
    out.n_F=0; out.C=0; out.J=0; out.source_J=0; out.B=0; out.alpha=.
    out.accepted=0; out.status="not_requested"; out.critical_value=.
    out.target_positions=J(0,1,0); out.raw_draws=J(0,0,0)
    return(out)
}
struct cdid_multiplier_result scalar cdid_multiplier_finish(struct cdid_multiplier_state scalar s,
    real scalar return_raw_draws)
{
    struct cdid_multiplier_result scalar out
    real scalar j, B, sc, iq, se, valid_cov
    real colvector ordered, maximum, standardized
    if(s.finished | s.cursor!=s.B | (return_raw_draws!=0 & return_raw_draws!=1))
        _cdid_infer_fail("Multiplier stream must contain all draws exactly once before finalization.",498)
    out.n_F=s.n_F; out.C=s.C; out.J=s.J; out.source_J=s.source_J
    out.B=s.B; out.alpha=s.alpha; out.target_positions=s.target_positions
    out.q25_index=cdid_type1_index(s.B,.25); out.q75_index=cdid_type1_index(s.B,.75)
    out.critical_index=cdid_type1_index(s.B,1-s.alpha)
    out.iqr_normalizer=invnormal(.75)-invnormal(.25)
    out.q25=out.q75=out.raw_iqr=out.raw_scale=out.candidate_se=J(s.J,1,.)
    out.bootstrap_scale=out.bootstrap_se=out.analytic_se=J(s.J,1,.)
    out.analytic_status=out.bootstrap_status=out.target_status=J(s.J,1,"")
    out.target_accepted=J(s.J,1,0); out.sup_t=J(s.B,1,.)
    out.critical_value=.; out.accepted=0; out.status="refused"
    valid_cov=!missing(s.analytic_variance)
    for(j=1;j<=s.J;j++) {
        if(!valid_cov | s.analytic_variance[j]>=.) out.analytic_status[j]="nonfinite_analytical_variance"
        else if(s.analytic_variance[j]<=0) out.analytic_status[j]="zero_analytical_variance"
        else { out.analytic_status[j]="ok"; out.analytic_se[j]=sqrt(s.analytic_variance[j]); }
        if(missing(s.raw_draws[,j])) out.bootstrap_status[j]="nonfinite_bootstrap_scale"
        else {
            ordered=sort(s.raw_draws[,j],1)
            out.q25[j]=ordered[out.q25_index]; out.q75[j]=ordered[out.q75_index]
            iq=out.q75[j]-out.q25[j]; sc=iq/out.iqr_normalizer; se=sc*sqrt(s.C)/s.n_F
            out.raw_iqr[j]=iq; out.raw_scale[j]=sc
            if(iq>=. | sc>=. | se>=.) out.bootstrap_status[j]="nonfinite_bootstrap_scale"
            else if(iq<=0 | sc<=0 | se<=0) {
                out.bootstrap_status[j]="zero_bootstrap_scale"; out.candidate_se[j]=se
            }
            else { out.bootstrap_status[j]="ok"; out.candidate_se[j]=se; }
        }
        out.target_accepted[j]=out.analytic_status[j]=="ok" & out.bootstrap_status[j]=="ok"
        if(out.target_accepted[j]) out.target_status[j]="ok"
        else if(out.analytic_status[j]!="ok" & out.bootstrap_status[j]!="ok") {
            if(out.analytic_status[j]=="zero_analytical_variance" & out.bootstrap_status[j]=="zero_bootstrap_scale")
                out.target_status[j]="zero_analytical_and_bootstrap_scale"
            else out.target_status[j]=out.analytic_status[j]+"__"+out.bootstrap_status[j]
        }
        else if(out.analytic_status[j]!="ok") out.target_status[j]=out.analytic_status[j]
        else out.target_status[j]=out.bootstrap_status[j]
    }
    if(all(out.target_accepted)) {
        maximum=J(s.B,1,0)
        for(j=1;j<=s.J;j++) {
            standardized=abs(s.raw_draws[,j]):/out.raw_scale[j]
            if(missing(standardized)) {
                out.target_accepted[j]=0; out.bootstrap_status[j]="nonfinite_standardized_draw"
                out.target_status[j]="nonfinite_standardized_draw"
            }
            else maximum=rowmax((maximum,standardized))
        }
        if(all(out.target_accepted)) {
            out.accepted=1; out.status="ok"; out.bootstrap_scale=out.raw_scale
            out.bootstrap_se=out.candidate_se; out.sup_t=maximum
            out.critical_value=cdid_type1_quantile(maximum,1-s.alpha)
        }
    }
    out.raw_draws=(return_raw_draws ? s.raw_draws : J(0,s.J,0))
    s.finished=1
    return(out)
}
/* bstrap=. is an explicit translation of the private R NULL flag law.
   Public boolean options pass0/1; no omitted default is read here. */
struct cdid_infer_flags scalar cdid_infer_resolve_flags(real scalar bstrap,
    real scalar cband, string scalar endpoint)
{
    struct cdid_infer_flags scalar f
    if((bstrap!=0 & bstrap!=1 & bstrap!=.) | (cband!=0 & cband!=1) |
        (endpoint!="cells" & endpoint!="aggte_bins" & endpoint!="simple"))
        _cdid_infer_fail("Invalid explicit bootstrap/band flags or inference endpoint.",498)
    f.endpoint=endpoint; f.released_false_true=0
    if(bstrap==.) {
        f.bootstrap=cband; f.simultaneous=cband; f.draws=cband
        f.law="private_bootstrap_iff_cband"
    }
    else if(bstrap | !cband) {
        f.bootstrap=bstrap; f.simultaneous=cband & endpoint!="simple"; f.draws=bstrap
        f.law="explicit_bstrap_and_cband"
    }
    else {
        f.bootstrap=0; f.released_false_true=1
        if(endpoint=="aggte_bins") {
            f.simultaneous=1; f.draws=1
            f.law="released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band"
        }
        else {
            f.simultaneous=0; f.draws=0
            f.law="released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws"
        }
    }
    return(f)
}
struct cdid_infer_result scalar cdid_infer_point_only(real colvector estimates,
    real scalar n_F, real scalar C, real scalar alpha)
{
    struct cdid_infer_result scalar out
    real scalar Jt
    Jt=rows(estimates)
    if(n_F<1 | n_F>=. | n_F!=floor(n_F) | C<1 | C>=. | C!=floor(C) |
        alpha<=0 | alpha>=1 | alpha>=.) _cdid_infer_fail("Invalid point-only frame counts or alpha.",498)
    out.n_F=n_F; out.C=C; out.J=Jt; out.estimates=estimates
    out.se=out.analytic_se=out.z_stat=out.p_value=out.ci_lower=out.ci_upper=
        out.band_lower=out.band_upper=J(Jt,1,.)
    out.bootstrap=0; out.simultaneous=0; out.draws_used=0; out.n_draws=.
    out.normal_critical=invnormal(1-alpha/2); out.band_critical=out.normal_critical
    out.covariance=J(0,0,0); out.Omega=J(0,0,0)
    out.status="inference_not_computed"; out.reason="point estimates only; no IF or multiplier work"
    out.flag_law="compute_inffunc_false_point_estimates_only"; out.band_family="none"; out.warning=""
    return(out)
}
struct cdid_infer_result scalar cdid_infer_report(real colvector estimates,
    struct cdid_if_summary scalar cov, struct cdid_infer_flags scalar flags,
    struct cdid_multiplier_result scalar mult, real scalar alpha)
{
    struct cdid_infer_result scalar out
    real scalar Jt, deterministic_all
    real colvector deterministic, positive
    Jt=cov.J
    if(rows(estimates)!=Jt | cols(estimates)!=1 | missing(estimates) | alpha<=0 | alpha>=1 | alpha>=.)
        _cdid_infer_fail("Finite inference estimates must align with the complete target set and alpha.",498)
    out=cdid_infer_point_only(estimates,cov.n_F,cov.C,alpha)
    out.covariance=cov.covariance; out.Omega=cov.Omega
    out.analytic_se=cov.analytic_se; out.se=cov.analytic_se
    out.flag_law=flags.law; out.status="unavailable"; out.reason=""
    deterministic=(cov.absolute_if_sum':==0):&(estimates:==0)
    positive=_cdid_infer_index(cov.analytic_se:>0)
    if(any((cov.analytic_se:<=0):&!deterministic)) {
        out.reason="nonconstant_targets_require_positive_marginal_inference_scales"
        out.se=J(Jt,1,.)
        return(out)
    }
    deterministic_all=all(deterministic)
    out.normal_critical=(deterministic_all ? 0 : invnormal(1-alpha/2))
    out.band_critical=out.normal_critical
    if(flags.draws & rows(positive)) {
        if(!mult.accepted | mult.n_F!=cov.n_F | mult.C!=cov.C | mult.source_J!=cov.J |
            mult.alpha!=alpha | rows(mult.target_positions)!=rows(positive)) {
            out.reason="complete_multiplier_family_refused_or_missing"; out.se=J(Jt,1,.)
            return(out)
        }
        if(any(mult.target_positions:!=positive)) {
            out.reason="complete_multiplier_family_order_or_membership_mismatch"; out.se=J(Jt,1,.)
            return(out)
        }
        if(flags.bootstrap) out.se[positive]=mult.bootstrap_se
        if(flags.simultaneous) out.band_critical=mult.critical_value
        out.draws_used=1; out.n_draws=mult.B
        out.bootstrap=flags.bootstrap; out.simultaneous=flags.simultaneous
        if(flags.released_false_true & flags.endpoint=="aggte_bins")
            out.warning="Used bootstrap procedure to compute simultaneous confidence band"
    }
    out.ci_lower=estimates-out.normal_critical*out.se; out.ci_upper=estimates+out.normal_critical*out.se
    out.band_lower=estimates-out.band_critical*out.se; out.band_upper=estimates+out.band_critical*out.se
    out.z_stat=J(Jt,1,.); out.p_value=J(Jt,1,.)
    if(rows(positive)) {
        out.z_stat[positive]=estimates[positive]:/out.se[positive]
        out.p_value[positive]=2:*normal(-abs(out.z_stat[positive]))
    }
    out.status="computed"
    out.band_family=(deterministic_all ? "deterministic_exact_zero_over_supported_ids" :
        (out.simultaneous ? "simultaneous_multiplier_sup_t_over_supported_ids" : "pointwise_normal_over_supported_ids"))
    if(missing(out.se) | missing(out.ci_lower) | missing(out.ci_upper) |
        missing(out.band_lower) | missing(out.band_upper)) {
        out.status="unavailable"; out.reason="nonfinite_inference_intervals"
    }
    return(out)
}
/* Computed exact inverse1-norm via one RHS at a time, no inverse matrix.
   R reports LAPACK dtrcon's estimate. Use U=F' for Mata's lower Cholesky F. */
real scalar _cdid_infer_rcond_upper(real matrix U)
{
    real scalar j, n, rk, maxcol, z, norm_U
    real colvector b, x
    real matrix Xinv
    n=rows(U); maxcol=0; norm_U=max(colsum(abs(U)))
    if(norm_U<=0 | norm_U>=.) return(0)
    if(n>=1 & n<=1024) {
        Xinv=solveupper(U,I(n),rk=.,0)
        if(rk!=n | missing(Xinv)) return(0)
        maxcol=max(colsum(abs(Xinv)))
        if(maxcol<=0 | maxcol>=.) return(0)
        return((1/norm_U)/maxcol)
    }
    for(j=1;j<=n;j++) {
        b=J(n,1,0); b[j]=1
        x=solveupper(U,b,rk=.,0)
        if(rk!=n | missing(x)) return(0)
        z=sum(abs(x)); if(z>=.) return(0)
        maxcol=max((maxcol,z))
    }
    if(maxcol<=0) return(0)
    return((1/norm_U)/maxcol)
}
/* theta/V/available are the FULL declared nonreference pretrend set.
   Structural reference identities are excluded by the reporting caller.
   No missing restriction is dropped; no MP/Tikhonov or interval route. */
struct cdid_wald_result scalar cdid_wald_full(real colvector theta, real matrix V,
    real colvector available, real scalar independent_rows, real scalar max_matrix_entries)
{
    struct cdid_wald_result scalar out
    real scalar K, entries, rk
    real matrix F
    real colvector x
    K=rows(theta); out.df=K; out.statistic=.; out.p_value=.; out.rcond_upper_1=.
    out.solve_residual=J(0,1,0); out.covariance_rank_upper=.
    out.status="unavailable"; out.reason=""
    if(K==0) {
        out.status="not_applicable"; out.reason="no_declared_nonreference_pretrend_comparison"
        return(out)
    }
    if(rows(available)!=K | cols(available)!=1 | missing(available) |
        any((available:!=0):&(available:!=1)) | independent_rows<1 | independent_rows>=. |
        independent_rows!=floor(independent_rows) | max_matrix_entries<1 | max_matrix_entries>=. |
        max_matrix_entries!=floor(max_matrix_entries) | max_matrix_entries>2^53-1)
        _cdid_infer_fail("Wald availability, independent row count and workspace must be explicit valid controls.",498)
    if(any(!available)) {
        out.reason="declared_pretrend_unavailable_no_available_only_subset"
        return(out)
    }
    out.covariance_rank_upper=independent_rows-1
    if(K>out.covariance_rank_upper) {
        out.reason="declared_restrictions_exceed_centered_covariance_rank_bound"
        return(out)
    }
    if(cols(theta)!=1 | rows(V)!=K | cols(V)!=K | missing(theta) | missing(V)) {
        out.reason="declared_pretrend_estimates_or_joint_covariance_not_finite"
        return(out)
    }
    entries=_cdid_infer_cells(4*K,K,max_matrix_entries,0)
    if(!all(V:==V')) {
        out.reason="joint_pretrend_covariance_not_symmetric"
        return(out)
    }
    F=cholesky(V)
    if(missing(F) | any(diagonal(F):<=0)) {
        out.reason="full_joint_covariance_has_no_usable_Cholesky_no_pseudoinverse"
        return(out)
    }
    out.rcond_upper_1=_cdid_infer_rcond_upper(F')
    if(out.rcond_upper_1<1e-10 | out.rcond_upper_1>=.) {
        out.reason="full_joint_covariance_fails_existing_Cholesky_conditioning_policy"
        return(out)
    }
    x=solvelower(F,theta,rk=.,0); x=solveupper(F',x,rk=.,0)
    if(rk!=K | missing(x)) { out.reason="full_joint_covariance_solve_unusable"; return(out); }
    out.statistic=quadcross(theta,x); out.solve_residual=V*x-theta
    out.p_value=chi2tail(K,out.statistic)
    if(out.statistic>=. | out.statistic<0 | out.p_value>=. | missing(out.solve_residual)) {
        out.statistic=.; out.p_value=.; out.reason="full_set_Wald_not_finite_and_nonnegative"
        return(out)
    }
    out.status="computed"
    return(out)
}
end
