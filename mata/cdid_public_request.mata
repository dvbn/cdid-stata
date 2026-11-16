*! Private public-request/auxiliary plan, native current StataNow19.5.
*! GPL-3 adaptation of Benatia, Bellego, Cuerrier, Dortet-Bernadet R cdid,
*! : chained_did_request.R901 and
*! chained_ipw_pipeline.R9(ipw_auxiliary_spec),3241(request ledger).
*! No nuisance/model/solver/IF fitting, availability filtering or public command.
version 19.5
mata:
struct cdid_public_request {
    string scalar control_group, base_period, pt_assumption
    real scalar anticipation, pt_start, pt_tags_enabled
    real colvector calendar, cohorts, panel_cohorts, target_cohorts, cohorts_without_post, target_ref
    real matrix post, aux, selection, report
    real colvector post_ref, post_K_ptr, post_K_index, aux_K_ptr, aux_K_index
    string colvector selection_reason, report_reason, pt_window_tag
    real colvector pt_operator_covered
}
struct cdid_public_settings {
    string scalar method, gmm_weighting, control_group, base_period, pt_assumption
    real scalar anticipation, pt_start, alpha, bstrap, cband, biters, compute_if
    string scalar selection, selection_canonical, selection_model, selection_pooling
    string scalar selection_universe, selection_est_method, selection_justification
    string scalar inactive_arguments, gmm_covariance, imputation_pooling, outcome_pooling
}
void _cdid_pr_fail(string scalar cls, string scalar message)
{
    errprintf("%s\n", message)
    exit(459)
}
string scalar _cdid_pr_num(real scalar x)
{
    return(x==. ? "Inf" : strtrim(strofreal(x,"%24.17g")))
}
void _cdid_pr_budget(real scalar used, real scalar max_entries)
{
    if(max_entries<1 | max_entries>=. | max_entries!=floor(max_entries))
        _cdid_pr_fail("cdid_chained_request_input_error","max_entries must be an explicit positive integer.")
    if(used>=. | used>max_entries) _cdid_pr_fail("cdid_ipw_pipeline_bound_error","Declared request/auxiliary numeric payload exceeds caller max_entries.")
}
struct cdid_public_settings scalar cdid_public_defaults()
{
    struct cdid_public_settings scalar S
    S.method="dr"; S.gmm_weighting="optimal"; S.control_group="nevertreated"
    S.base_period="varying"; S.pt_assumption="post"; S.anticipation=0; S.pt_start=.
    S.alpha=.05; S.bstrap=S.cband=S.compute_if=1; S.biters=1000
    S.selection="none"; S.selection_canonical="none"; S.selection_model="pair_logit"
    S.selection_pooling="separate"; S.selection_universe="F"; S.selection_est_method="ipw"
    S.imputation_pooling=""; S.outcome_pooling="separate"; S.selection_justification=""
    S.inactive_arguments="imputation_pooling"; S.gmm_covariance="tikhonov"
    return(S)
}
real scalar _cdid_pr_reference(real colvector cal, real scalar g, real scalar ant)
{
    real scalar j, ref
    ref=.
    for(j=1;j<=rows(cal);j++) if(_cdid_cmp_sum(cal[j],ant,g)<0) ref=cal[j]
    return(ref)
}
real colvector _cdid_pr_pool(real colvector coh, real scalar g, real scalar t,
    real scalar ant, string scalar control)
{
    real colvector out
    real scalar k, n
    out=J(rows(coh),1,0);n=0
    for(k=1;k<=rows(coh);k++) {
        if(control=="nevertreated") {
            if(coh[k]==.) {
                n++
                out[n]=k
            }
        }
        else if(coh[k]!=g & _cdid_cmp_sum(t,ant,coh[k])<0) {
            n++
            out[n]=k
        }
    }
    return(n ? out[|1\n|] : J(0,1,0))
}
string scalar _cdid_pr_key(real scalar gi, real scalar si, real scalar ti, real colvector K)
{
    real scalar j
    string scalar key
    key=_cdid_istr(gi)+":"+_cdid_istr(si)+":"+_cdid_istr(ti)
    for(j=1;j<=rows(K);j++) key=key+":"+_cdid_istr(K[j])
    return(key)
}
struct cdid_public_request scalar cdid_public_request_build(struct cdid_population scalar P,
    real colvector target_override, string scalar control, real scalar ant,
    string scalar base, string scalar pt, real scalar ptstart, real scalar max_entries)
{
    struct cdid_public_request scalar R
    real colvector ix, finite, withpost, K, fromix, pick, orderix
    real scalar Jg,T,S,Q,i,j,f,at,fi,g,t,ref,from,si,ti,gi,status,sgn,q,totalK,postK,npost,rr,base_entries
    real matrix aux,selection,report
    string colvector reason,rreason
    string scalar key,why,universe
    transmorphic scalar seen,pools,postpools
    _cdid_pr_budget(0,max_entries)
    if(P.N<1 | P.R<1 | rows(P.calendar)<1) _cdid_pr_fail("cdid_chained_request_input_error","data must be nonempty.")
    if(control!="nevertreated" & control!="notyettreated") _cdid_pr_fail("cdid_chained_request_input_error","controlgroup() must be nevertreated or notyettreated.")
    if(base!="varying" & base!="universal") _cdid_pr_fail("cdid_chained_request_input_error","baseperiod() must be varying or universal.")
    if(ant<0 | ant>=. | ant!=floor(ant)) _cdid_pr_fail("cdid_chained_request_input_error","anticipation() must be a nonnegative integer.")
    if(pt!="post" & pt!="all" & pt!="ctrl") _cdid_pr_fail("cdid_chained_request_input_error","pt_assumption must be post, ctrl, or all.")
    if(pt=="ctrl" & ptstart<.) if(!any(P.calendar:==ptstart)) _cdid_pr_fail("cdid_chained_request_pt_start","ptstart() must be one of the periods in the data.")
    if((pt=="post" & ptstart<.) | (pt=="all" & (ptstart>=. | !any(P.calendar:==ptstart))))
        _cdid_pr_fail("cdid_chained_request_pt_start","ptassumption(all) requires ptstart(#), which must be one of the periods in the data; ptstart() is not allowed with ptassumption(post).")
    R.control_group=control;R.base_period=base;R.anticipation=ant;R.pt_assumption=pt;R.pt_start=ptstart
    R.calendar=P.calendar;R.cohorts=cdid_cipw_cohorts(P);T=rows(P.calendar)
    ix=vec(selectindex((P.row_count:>0):&(P.cohort_known:!=0)))
    R.panel_cohorts=uniqrows(sort(P.cohort[ix],1))
    finite=select(R.panel_cohorts,R.panel_cohorts:<.)
    withpost=select(finite,finite:<=P.calendar[T]);R.cohorts_without_post=select(finite,finite:>P.calendar[T])
    if(rows(target_override)) {
        if(cols(target_override)!=1 | missing(target_override)) _cdid_pr_fail("cdid_chained_request_input_error","Target cohorts must be finite numeric DATA cohorts with post dates.")
        for(j=1;j<=rows(target_override);j++) if(!any(withpost:==target_override[j]))
            _cdid_pr_fail("cdid_chained_request_input_error","targetcohorts() must list treated cohorts present in the data and observed at or after treatment.")
        R.target_cohorts=uniqrows(sort(target_override,1))
    }
    else R.target_cohorts=withpost
    Jg=rows(R.target_cohorts)
    if(!Jg) _cdid_pr_fail("cdid_chained_request_input_error","No treated cohort has a calendar date at or after its treatment date; nothing can be estimated.")
    if(findexternal("__cdid_mref_state")!=NULL) cdid_mrpa_set_request(R.target_cohorts,P.calendar,ant)
    R.target_ref=J(Jg,1,.);S=0
    for(j=1;j<=Jg;j++) {
        R.target_ref[j]=_cdid_pr_reference(P.calendar,R.target_cohorts[j],ant)
        for(i=1;i<=T;i++) if(P.calendar[i]<R.target_cohorts[j]) S=S+1+(i>1)
    }
    R.post=cdid_cipw_post_targets(P.calendar,R.target_cohorts);npost=rows(R.post)
    base_entries=32*S+16*npost+8*(T+rows(R.cohorts)+Jg)
    _cdid_pr_budget(base_entries,max_entries)
    R.post_ref=J(npost,1,.);R.post_K_ptr=J(npost+1,1,1)
    postpools=asarray_create("real",1);postK=0
    for(i=1;i<=npost;i++) {
        R.post_ref[i]=_cdid_pr_reference(P.calendar,R.post[i,1],ant)
        K=_cdid_pr_pool(R.cohorts,R.post[i,1],R.post[i,2],ant,control)
        postK=postK+rows(K);_cdid_pr_budget(base_entries+2*postK,max_entries)
        asarray(postpools,i,K);R.post_K_ptr[i+1]=postK+1
    }
    R.post_K_index=J(postK,1,0)
    for(i=1;i<=npost;i++) if(R.post_K_ptr[i]<R.post_K_ptr[i+1]) R.post_K_index[|R.post_K_ptr[i]\R.post_K_ptr[i+1]-1|]=asarray(postpools,i)
    aux=J(S,3,.);selection=J(S,9,.);reason=J(S,1,"")
    seen=asarray_create("string",1);pools=asarray_create("real",1)
    Q=at=totalK=0;universe=""
    for(j=1;j<=rows(R.cohorts);j++) universe=universe+(j>1 ? "," : "")+_cdid_pr_num(R.cohorts[j])
    for(f=1;f<=2;f++) {
        fi=0
        for(j=1;j<=Jg;j++) {
            g=R.target_cohorts[j];ref=R.target_ref[j]
            ix=vec(selectindex(R.cohorts:==g));gi=ix[1]
            for(i=1+(f==1);i<=T;i++) {
                t=P.calendar[i]
                if(t>=g) continue
                at++;fi++;from=.;q=0;sgn=.;status=0;why=""
                if(ref>=.) {
                    status=2
                    why="no_admissible_reference: no stored calendar date t < g - anticipation for cohort g="+_cdid_pr_num(g)
                }
                else {
                    from=(f==1 ? P.calendar[i-1] : ref)
                    if(f==2 & t==ref) {
                        status=1
                        sgn=0
                    }
                    else {
                        K=_cdid_pr_pool(R.cohorts,g,max((t,from)),ant,control)
                        if(!rows(K)) {
                            status=3
                            why="no_clean_control_cohort: no cohort of the universe {"+universe+"} other than g="+_cdid_pr_num(g)+" is untreated at date "+_cdid_pr_num(max((t,from)))+" with anticipation "+_cdid_pr_num(ant)+" (comparison dates "+_cdid_pr_num(min((t,from)))+"->"+_cdid_pr_num(max((t,from)))+")"
                        }
                        else {
                            fromix=vec(selectindex(P.calendar:==from));si=min((i,fromix[1]));ti=max((i,fromix[1]))
                            key=_cdid_pr_key(gi,si,ti,K)
                            if(asarray_contains(seen,key)) q=asarray(seen,key)
                            else {
                                Q++;q=Q;aux[Q,]=(g,P.calendar[si],P.calendar[ti]);asarray(seen,key,Q)
                                totalK=totalK+rows(K);_cdid_pr_budget(base_entries+2*(postK+totalK),max_entries)
                                asarray(pools,Q,K)
                            }
                            sgn=(t>from ? 1 : -1)
                        }
                    }
                }
                selection[at,]=(f,g,t,ref,from,q,sgn,status,fi);reason[at]=why
            }
        }
    }
    R.aux=(Q ? aux[|1,1\Q,3|] : J(0,3,.));R.selection=selection;R.selection_reason=reason
    R.aux_K_ptr=J(Q+1,1,1);R.aux_K_index=J(totalK,1,0);rr=1
    for(q=1;q<=Q;q++) {
        K=asarray(pools,q);R.aux_K_index[|rr\rr+rows(K)-1|]=K;rr=rr+rows(K);R.aux_K_ptr[q+1]=rr
    }
    pick=vec(selectindex(selection[,1]:==(base=="varying" ? 1 : 2)))
    report=J(rows(pick)+npost,10,.);rreason=J(rows(report),1,"")
    for(i=1;i<=rows(pick);i++) {
        q=pick[i];report[i,]=(selection[q,2..5],(selection[q,8]==1 ? 1 : 2),0,selection[q,6..8],q);rreason[i]=reason[q]
    }
    for(i=1;i<=npost;i++) report[rows(pick)+i,]=(R.post[i,],R.post_ref[i],R.post_ref[i],3,i,0,1,0,0)
    orderix=order(report,(1,2));R.report=report[orderix,];R.report_reason=rreason[orderix]
    R.pt_tags_enabled=(pt=="all");R.pt_window_tag=J(rows(report),1,"");R.pt_operator_covered=J(rows(report),1,.)
    if(pt=="all") for(i=1;i<=rows(report);i++) {
        if(R.report[i,5]==3) R.pt_window_tag[i]="causal"
        else if(R.report[i,5]==1) R.pt_window_tag[i]="reference"
        else {
            from=R.report[i,4];t=R.report[i,2];g=R.report[i,1];q=R.report[i,7]
            R.pt_operator_covered[i]=0
            if(from>=.) R.pt_window_tag[i]="unavailable_reference"
            else if(_cdid_cmp_sum(max((t,from)),ant,g)>=0) R.pt_window_tag[i]="anticipation"
            else R.pt_window_tag[i]=(min((t,from))>=ptstart ? "in_pt_window" : "pre_pt_start")
            if(q>0) {
                K=R.aux_K_index[|R.aux_K_ptr[q]\R.aux_K_ptr[q+1]-1|]
                if(rows(K)==1) if(R.cohorts[K[1]]==.) R.pt_operator_covered[i]=1
            }
        }
    }
    return(R)
}
struct cdid_public_request scalar cdid_public_request_tags(struct cdid_public_request scalar R,struct cdid_cat scalar C)
{
    real scalar i,q,r,p
    real colvector K,other
    if(rows(C.cohorts)!=rows(R.cohorts) | !all(C.cohorts:==R.cohorts)) _cdid_pr_fail("cdid_chained_request_input_error","Request/catalogue cohort dictionaries disagree.")
    if(R.pt_assumption!="all") return(R)
    for(i=1;i<=rows(R.report);i++) if(R.report[i,5]==2 & !R.pt_operator_covered[i]) {
        q=R.report[i,7]
        if(q<1) continue
        K=R.aux_K_index[|R.aux_K_ptr[q]\R.aux_K_ptr[q+1]-1|]
        for(r=1;r<=rows(C.r_g);r++) if(C.cohorts[C.r_g[r]]==R.report[i,1]) {
            p=C.r_pool[r];other=C.pool_entries[|C.pool_start[p]\C.pool_start[p]+C.pool_size[p]-1|]
            if(rows(K)==rows(other)) if(all(K:==other)) {
                R.pt_operator_covered[i]=1
                break
            }
        }
    }
    return(R)
}

/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
void cdid_mrpa_set_request(transmorphic a1, transmorphic a2, transmorphic a3)
{
    cdid_mref_set_request(a1, a2, a3)
}

end
