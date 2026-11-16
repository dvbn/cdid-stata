*! Private fitted-provider to full reporting IF composition; no nuisance refits.
*! Requires the separately qualified baseline provider and public report scatter.
version 19.5
mata:
struct cdid_result_state scalar cdid_public_provider_report(
    struct cdid_public_bundle scalar B,struct cdid_baseline_provider scalar fit,
    real colvector post_point,real colvector post_available,string colvector post_status,
    string colvector post_reason,real colvector post_est_ids,real matrix post_loadings,
    real scalar compute_if,struct cdid_state_limits scalar limits,real scalar chunk,real scalar cap,
    | real matrix retained_if)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    real colvector aux_point,aux_available,active,qids,postpos,auxpos,consumers,hit,uid,aux_in_post,rem_consumers,rem_pos
    real matrix values,postIF,auxIF,ret_if
    real scalar Q,M,Jt,i,j,q,lo,hi,effective,perrow,retained,use_ret
    Q=rows(B.request.aux);M=rows(fit.plan.post_consumer);Jt=rows(post_est_ids)
    ret_if=(args()>=13 ? retained_if : J(0,0,.))
    if(rows(post_loadings)!=Jt | cols(post_loadings)!=M | rows(fit.plan.aux_consumer)!=Q |
        chunk<1 | chunk>=. | chunk!=floor(chunk) | cap<1 | cap>=. | cap!=floor(cap))
        _cdid_pr_fail("cdid_native_reporting_input","Provider/post loading identities or replay controls do not align.")
    if(rows(post_est_ids)) {
        _cdid_rb_ids(post_est_ids,rows(B.request.post))
        if(any(post_available[post_est_ids]:!=1) | rows(post_est_ids)!=sum(post_available))
            _cdid_pr_fail("cdid_native_reporting_input","Estimable post IDs must equal the full supported post set.")
    }
    else if(sum(post_available)) _cdid_pr_fail("cdid_native_reporting_input","Supported post set has no passed target loading.")
    aux_point=J(Q,1,.);aux_available=fit.plan.aux_supported
    active=_cdid_rb_index(aux_available)
    if(rows(active)) aux_point[active]=fit.points[fit.plan.aux_consumer[active]]
    buffer=cdid_public_report_start(B,post_point,post_available,post_status,post_reason,
        aux_point,aux_available,fit.plan.aux_reason,compute_if,limits)
    if(compute_if) {
        qids=J(0,1,.)
        for(i=1;i<=rows(B.request.report);i++) if(B.request.report[i,5]==2 & buffer.payload.available[i])
            qids=qids\B.request.report[i,7]
        if(rows(qids)) qids=uniqrows(sort(qids,1))
        consumers=(Jt ? fit.plan.post_consumer : J(0,1,.))
        if(rows(qids)) consumers=consumers\fit.plan.aux_consumer[qids]
        if(rows(consumers)) consumers=uniqrows(sort(consumers,1))
        postpos=J((Jt ? M : 0),1,.)
        for(j=1;j<=rows(postpos);j++) {
            hit=_cdid_rb_index(consumers:==fit.plan.post_consumer[j]);postpos[j]=hit[1]
        }
        auxpos=J(rows(qids),1,.)
        for(j=1;j<=rows(qids);j++) {
            hit=_cdid_rb_index(consumers:==fit.plan.aux_consumer[qids[j]]);auxpos[j]=hit[1]
        }
        // Retained result + provider + target loading, plus conservative bounded
        // provider/scatter scratch. Never construct originalF-by-allmoments IF.
        retained=cdid_state_cost(buffer.payload)[1]+fit.workspace_entries+length(post_loadings)
        perrow=10*(rows(consumers)+cols(fit.X_A)+cols(fit.Z_A)+4)+3*(Jt+rows(qids)+rows(B.request.report))
        if(retained>=cap | perrow>cap-retained) _cdid_rs_fail("Reported IF/provider replay exceeds explicit numeric workspace.",908)
        effective=min((chunk,fit.chunk_rows,floor((cap-retained)/perrow)))
        use_ret=(rows(ret_if)==fit.n_F & cols(ret_if)==M & fit.weight_row_domain=="F")
        if(use_ret) {
            aux_in_post=J(rows(qids),1,0);rem_consumers=J(0,1,.);rem_pos=J(rows(qids),1,0)
            for(j=1;j<=rows(qids);j++) {
                hit=_cdid_rb_index(fit.plan.post_consumer:==fit.plan.aux_consumer[qids[j]])
                if(rows(hit)) aux_in_post[j]=hit[1]
                else rem_consumers=rem_consumers\fit.plan.aux_consumer[qids[j]]
            }
            if(rows(rem_consumers)) {
                rem_consumers=uniqrows(sort(rem_consumers,1))
                for(j=1;j<=rows(qids);j++) if(!aux_in_post[j]) {
                    hit=_cdid_rb_index(rem_consumers:==fit.plan.aux_consumer[qids[j]])
                    rem_pos[j]=hit[1]
                }
            }
            for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                postIF=(Jt ? (fit.core_to_frame:*ret_if[uid,])*post_loadings' : J(hi-lo+1,0,.))
                if(rows(qids)) {
                    auxIF=J(hi-lo+1,rows(qids),.)
                    values=(rows(rem_consumers) ? cdid_baseline_if_chunk(fit,lo,hi,rem_consumers,"frame") : J(hi-lo+1,0,.))
                    for(j=1;j<=rows(qids);j++) {
                        if(aux_in_post[j]) auxIF[,j]=fit.core_to_frame:*ret_if[uid,aux_in_post[j]]
                        else auxIF[,j]=values[,rem_pos[j]]
                    }
                }
                else auxIF=J(hi-lo+1,0,.)
                cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
            }
        }
        else {
            for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                values=(rows(consumers) ? cdid_baseline_if_chunk(fit,lo,hi,consumers,"frame") : J(hi-lo+1,0,.))
                postIF=(Jt ? values[,postpos']*post_loadings' : J(hi-lo+1,0,.))
                auxIF=(rows(qids) ? values[,auxpos'] : J(hi-lo+1,0,.))
                cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
            }
        }
    }
    out=cdid_public_report_finish(buffer,limits)
    out.numeric_settings=out.numeric_settings\fit.weight_rows\fit.normalization_count\fit.n_ps_fits\fit.n_or_fits\M\Q
    out.numeric_setting_names=out.numeric_setting_names\"GMM_replay_rows"\"GMM_normalization_count"\"n_treatment_fits"\"n_pooled_OR_fits"\"M_full"\"Q_aux_union"
    out.text_settings=out.text_settings\("GMM_weight_row_domain",fit.weight_row_domain)\("provider_status",fit.status)
    out.source_versions=out.source_versions\("estimator","post-treatment and auxiliary comparisons; influence functions over all units")
    cdid_state_validate(out,limits)
    return(out)
}
end
