*! Public static numeric actual-change imputation composition; StataNow19.5.
*! R cdid,
*! chained_imputation_public.R and chained_ipw_pipeline.R. Scientific kernels unchanged.
version 19.5
mata:

/* One link per exact functional/recipe identity. Original raw post+both-family
   auxiliary descriptors remain in plan.raw and are not removed from GMM/reporting. */
void cdid_imputation_public_links(struct cdid_baseline_plan_state scalar plan,
    real matrix links,real matrix controls,string colvector ids)
{
    real scalar c,i,rep,k,g,n_ctrl,at
    real colvector K,reps
    links=J(plan.n_consumers,3,.)
    ids=J(plan.n_consumers,1,"")
    reps=J(plan.n_consumers,1,0)
    n_ctrl=0
    for(c=1;c<=plan.n_consumers;c++) {
        rep=0
        for(i=1;i<=plan.M+plan.Q;i++) if(plan.raw_supported[i] & plan.raw_consumer[i]==c) {
            rep=i
            break
        }
        if(!rep) _cdid_pr_fail("cdid_native_reporting_input","Imputation consumer has no complete-catalogue descriptor.")
        reps[c]=rep
        n_ctrl=n_ctrl+(plan.raw_K_ptr[rep+1]-plan.raw_K_ptr[rep])
    }
    controls=J(n_ctrl,2,.)
    at=0
    for(c=1;c<=plan.n_consumers;c++) {
        rep=reps[c]
        g=plan.cohorts[plan.raw[rep,1]]
        links[c,]=(plan.calendar[plan.raw[rep,2]],plan.calendar[plan.raw[rep,3]],g)
        K=plan.raw_K_index[plan.raw_K_ptr[rep]..plan.raw_K_ptr[rep+1]-1]
        for(k=1;k<=rows(K);k++) {
            at++
            controls[at,]=(c,(plan.cohorts[K[k]]>=. ? 0 : plan.cohorts[K[k]]))
        }
        ids[c]=(rep<=plan.M ? "cp04:catalog:"+_cdid_rb_key(rep) : "cp04:auxiliary:"+_cdid_rb_key(rep-plan.M))
    }
}


/* ONE shared catalogue owns scientific identities. This adapter only translates
   its qualified descriptors to the already dated link/reader API. */
void cdid_imp_reference_links(struct cdid_baseline_plan_state scalar plan,
    struct cdid_reference_catalog scalar catalog,real matrix links,
    real matrix controls,string colvector ids)
{
    real scalar c,i,k,g,n_ctrl,at
    real colvector K
    n_ctrl=(catalog.n_consumers ? catalog.K_ptr[catalog.n_consumers+1]-1 : 0)
    links=J(catalog.n_consumers,4,.);controls=J(n_ctrl,2,.);ids=J(catalog.n_consumers,1,"")
    at=0
    for(c=1;c<=catalog.n_consumers;c++) {
        i=catalog.representative[c];g=plan.cohorts[catalog.consumers[c,1]]
        links[c,]=(plan.calendar[catalog.consumers[c,2]],plan.calendar[catalog.consumers[c,3]],g,catalog.consumer_recipe[c])
        K=catalog.K_index[catalog.K_ptr[c]..catalog.K_ptr[c+1]-1]
        for(k=1;k<=rows(K);k++) {
            at++
            controls[at,]=(c,(plan.cohorts[K[k]]>=. ? 0 : plan.cohorts[K[k]]))
        }
        ids[c]=(i<=plan.M ? "cp04:catalog:"+_cdid_rb_key(i) : "cp04:auxiliary:"+_cdid_rb_key(i-plan.M))
    }
}
/* Recipe eligibility masks (N x n_recipes): recipe r reads the units eligible in
   its reference eligibility context; empty when no recipe is restricted. */
real matrix cdid_imp_reference_keep(struct cdid_reference_inputs scalar S)
{
    real matrix keep
    real colvector hit
    real scalar r,e
    if(cols(S.recipes)<3 | rows(S.elig_index)==0) return(J(0,0,.))
    keep=J(S.N,S.n_recipes,1)
    for(r=1;r<=S.n_recipes;r++) {
        e=S.recipes[r,3]
        if(e>0 & e<.) {
            hit=vec(selectindex(S.elig_index:==e))
            if(rows(hit)!=1) _cdid_imp_fail("Unknown reference eligibility context of a design block.",498)
            keep[,r]=!S.elig_excluded[,hit[1]]
        }
    }
    return(keep)
}
/* Recipes that no catalogue comparison uses (an auxiliary comparison made
   unavailable in its eligibility context) are dropped; the rest are renumbered. */
void cdid_imp_reference_compact(real matrix links,struct cdid_imp_recipe_spec colvector recipes,real matrix keep)
{
    struct cdid_imp_recipe_spec colvector out
    real colvector ru,map
    real scalar j
    if(rows(links)==0) return
    ru=uniqrows(links[,4])
    if(rows(ru)==rows(recipes)) return
    map=J(rows(recipes),1,0);map[ru]=(1::rows(ru))
    links[,4]=map[links[,4]]
    out=cdid_imp_recipe_spec(rows(ru),1)
    for(j=1;j<=rows(ru);j++) {out[j]=recipes[ru[j]];out[j].recipe_id=j;}
    recipes=out
    if(rows(keep)) keep=keep[,ru]
}
void cdid_imp_reference_registry(struct cdid_reference_inputs scalar S,
    struct cdid_imp_dated_registry scalar registry,struct cdid_imp_recipe_spec colvector recipes)
{
    real scalar role,slot,j,k,cal,r,z,x,w,p,q,intercept
    real colvector dates,hit,keys,units
    real matrix value
    string rowvector zn,xn
    cdid_ref_activate(S)
    intercept=S.roles[1].intercept;slot=S.roles[1].slot
    p=cols(S.designs[slot].value)-intercept
    zn=(p ? S.designs[slot].columns[intercept+1..intercept+p] : J(1,0,""))
    q=0;xn=J(1,0,"")
    if(S.roles[2].active) {q=cols(S.designs[S.roles[2].slot].value);xn=S.designs[S.roles[2].slot].columns;}
    registry=cdid_imp_dated_reg_init(S.N,p,q,intercept,zn,xn,"shared_reference_frozen_union")
    /* Each unique role/date slot is materialized ONCE; unused cohort rows keep
       their raw missing support. Recipes borrow these slots, never dense clones. */
    for(role=1;role<=2;role++)if(S.roles[role].active) {
        slot=S.roles[role].slot;keys=S.designs[slot].keys
        dates=uniqrows(sort(floor(keys:/(S.N+1)),1))
        for(j=1;j<=rows(dates);j++) {
            cal=dates[j];hit=_cdid_ref_ix(floor(keys:/(S.N+1)):==cal)
            units=keys[hit]:-cal*(S.N+1)
            value=J(S.N,cols(S.designs[slot].value),.)
            value[units,]=cdid_ref_design_rows(S,role,keys[hit])
            if(role==1) {
                if(intercept)value=(p ? value[,2..p+1] : J(S.N,0,.))
                cdid_imp_dated_reg_add_z(registry,cal,(cal ? S.calendar[cal] : .),value)
            }
            else cdid_imp_dated_reg_add_x(registry,cal,(cal ? S.calendar[cal] : .),value)
        }
    }
    dates=uniqrows(sort(floor(S.weight_keys:/(S.N+1)),1))
    for(j=1;j<=rows(dates);j++) {
        cal=dates[j];hit=_cdid_ref_ix(floor(S.weight_keys:/(S.N+1)):==cal)
        units=S.weight_keys[hit]:-cal*(S.N+1);value=J(S.N,1,.)
        value[units]=S.weights[hit]
        cdid_imp_dated_reg_add_w(registry,cal,(cal ? S.calendar[cal] : .),value)
    }
    recipes=cdid_imp_recipe_spec(S.n_recipes,1)
    for(r=1;r<=S.n_recipes;r++) {
        cal=S.recipes[r,1]
        z=_cdid_impd_find_slot(registry.z_data,(S.roles[1].dated ? cal : 0))
        x=(S.roles[2].active ? _cdid_impd_find_slot(registry.x_data,(S.roles[2].dated ? cal : 0)) : 0)
        w=_cdid_impd_find_slot(registry.w_data,S.recipes[r,2])
        recipes[r]=cdid_imp_recipe_spec_make(r,cal,S.recipes[r,2],(cal ? S.calendar[cal] : .),
            (S.recipes[r,2] ? S.calendar[S.recipes[r,2]] : .),z,x,w,"shared:"+_cdid_rb_key(r))
    }
}

/* Durable compact scientific provenance. Raw descriptors, controls, support,
   link aliases and points persist without primitive full-F-by-all-moment IF. */
void cdid_imputation_public_metadata(struct cdid_result_state scalar out,
    struct cdid_baseline_plan_state scalar plan,pointer(struct cdid_imp_backend scalar) scalar backend,
    struct cdid_cat scalar C)
{
    real scalar i,c,k,ni,no,np,nb,ncells,ngroups,nlocal,p_or,pscols,link_point,n_raw,b6,b5
    real colvector K,add_num
    string colvector coh_str,add_names
    string matrix add_txt
    string scalar key,control_key,id,source_key,actual_key,reason,act_prefix
    ni=no=np=nb=ncells=ngroups=nlocal=p_or=pscols=0
    if(backend!=NULL) {
        ni=(*backend).n_imputation_fits;no=(*backend).n_or_fits;np=(*backend).n_ps_fits
        nb=((*backend).has_or ? (*backend).n_blocks : 0)
        ncells=(*backend).n_cells;ngroups=(*backend).n_groups
        nlocal=(*backend).intercept*ncells;p_or=(*backend).or_p
        pscols=((*backend).has_ps ? (*backend).q : 0)
    }
    act_prefix="cmpcat:IMP-"+strupper(_cdid_rb_stored_text(out,"method"))+":"+_cdid_rb_stored_text(out,"imputation_pooling")+":"
    coh_str=J(rows(plan.cohorts),1,"")
    for(k=1;k<=rows(plan.cohorts);k++) coh_str[k]=_cdid_pr_num(plan.cohorts[k])
    n_raw=plan.M+plan.Q
    add_num=J(6*n_raw,1,.)
    add_names=J(6*n_raw,1,"")
    add_txt=J(5*n_raw,2,"")
    for(i=1;i<=n_raw;i++) {
        key=_cdid_rb_key(i);c=plan.raw_consumer[i]
        id=(i<=plan.M ? "cp04:catalog:"+key : "cp04:auxiliary:"+_cdid_rb_key(i-plan.M))
        K=plan.raw_K_index[plan.raw_K_ptr[i]..plan.raw_K_ptr[i+1]-1]
        if(rows(K)==1) control_key=coh_str[K[1]]
        else {
            control_key=""
            for(k=1;k<=rows(K);k++) control_key=control_key+(k>1 ? "," : "")+coh_str[K[k]]
        }
        link_point=(plan.raw_supported[i] ? (*backend).links[c].mu : .)
        b6=6*(i-1);b5=5*(i-1)
        add_num[b6+1]=link_point
        add_num[b6+2]=plan.cohorts[plan.raw[i,1]]
        add_num[b6+3]=plan.calendar[plan.raw[i,2]]
        add_num[b6+4]=plan.calendar[plan.raw[i,3]]
        add_num[b6+5]=c
        add_num[b6+6]=plan.raw_supported[i]
        add_names[b6+1]="source_link_point:"+key
        add_names[b6+2]="source_link_g:"+key
        add_names[b6+3]="source_link_s:"+key
        add_names[b6+4]="source_link_t:"+key
        add_names[b6+5]="source_link_consumer:"+key
        add_names[b6+6]="source_link_supported:"+key
        if(i<=plan.M) source_key=cdid_catalog_functional_key(C,C.r_op[i])
        else source_key="cmpcat:CP04-aux:g"+coh_str[plan.raw[i,1]]+":K"+control_key
        actual_key=act_prefix+source_key
        reason=(plan.raw_supported[i] ? "supported" : plan.aux_reason[i-plan.M])
        add_txt[b5+1,]=("source_link_id:"+key,id)
        add_txt[b5+2,]=("source_link_controls:"+key,control_key)
        add_txt[b5+3,]=("source_link_operator:"+key,source_key)
        add_txt[b5+4,]=("source_link_functional:"+key,actual_key)
        add_txt[b5+5,]=("source_link_reason:"+key,reason)
    }
    out.numeric_settings=out.numeric_settings\plan.n_F\plan.n_A\0\0\np\no\ni\ncells\ngroups\nlocal\nb\p_or\pscols\plan.M\plan.Q\plan.n_consumers\add_num
    out.numeric_setting_names=out.numeric_setting_names\"GMM_replay_rows"\"GMM_normalization_count"\"n_selection_fits"\"observation_logit_models"\
        "n_treatment_fits"\"n_pooled_OR_fits"\"n_imputation_fits"\"n_imputation_cells"\"n_imputation_groups"\"n_imputation_local_intercepts"\
        "n_OR_blocks"\"p_OR_retained"\"p_treatment"\"M_full"\"Q_aux_union"\"n_unique_consumers"\add_names
    out.text_settings=out.text_settings\("GMM_weight_row_domain","all original F; core psi=(n_A/N_F)*phi_F; covariance normalization n_A")\
        ("profile_kind","static_numeric")\("provider_status","ok")\
        ("imputation_prediction_universe","every required cohort row in original F")\
        ("observation_model","none: pure actual-change imputation")\
        ("duration_convention","stored-calendar index difference")\add_txt
    out.source_versions=out.source_versions\("estimator","change imputation for missing outcomes; one fit over all comparisons; influence functions over all units")
}

struct cdid_result_state scalar cdid_imputation_public_report(
    struct cdid_public_bundle scalar B,struct cdid_baseline_plan_state scalar plan,
    pointer(struct cdid_imp_backend scalar) scalar backend,real colvector points,
    real scalar fit_entries,
    real colvector post_point,real colvector post_available,string colvector post_status,
    string colvector post_reason,real colvector post_est_ids,real matrix post_loadings,
    real scalar compute_if,struct cdid_state_limits scalar limits,real scalar chunk,real scalar cap,
    | real matrix retained_if)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    struct cdid_imp_backend_if scalar replay
    real colvector aux_point,aux_available,active,qids,postpos,auxpos,consumers,hit,uid,aux_in_post,rem_consumers,rem_pos
    real matrix values,postIF,auxIF,ret_if
    real scalar Q,M,Jt,i,j,q,lo,hi,effective,perrow,retained,use_ret,scale_f
    Q=rows(B.request.aux);M=rows(plan.post_consumer);Jt=rows(post_est_ids)
    ret_if=(args()>=15 ? retained_if : J(0,0,.))
    if(rows(post_loadings)!=Jt | cols(post_loadings)!=M | rows(plan.aux_consumer)!=Q |
        chunk<1 | chunk>=. | chunk!=floor(chunk) | cap<1 | cap>=. | cap!=floor(cap))
        _cdid_pr_fail("cdid_native_reporting_input","Provider/post loading identities or replay controls do not align.")
    if(rows(post_est_ids)) {
        _cdid_rb_ids(post_est_ids,rows(B.request.post))
        if(any(post_available[post_est_ids]:!=1) | rows(post_est_ids)!=sum(post_available))
            _cdid_pr_fail("cdid_native_reporting_input","Estimable post IDs must equal the full supported post set.")
    }
    else if(sum(post_available)) _cdid_pr_fail("cdid_native_reporting_input","Supported post set has no passed target loading.")
    aux_point=J(Q,1,.);aux_available=plan.aux_supported
    active=_cdid_rb_index(aux_available)
    if(rows(active)) aux_point[active]=points[plan.aux_consumer[active]]
    buffer=cdid_public_report_start(B,post_point,post_available,post_status,post_reason,
        aux_point,aux_available,plan.aux_reason,compute_if,limits)
    if(compute_if) {
        qids=J(0,1,.)
        for(i=1;i<=rows(B.request.report);i++) if(B.request.report[i,5]==2 & buffer.payload.available[i])
            qids=qids\B.request.report[i,7]
        if(rows(qids)) qids=uniqrows(sort(qids,1))
        consumers=(Jt ? plan.post_consumer : J(0,1,.))
        if(rows(qids)) consumers=consumers\plan.aux_consumer[qids]
        if(rows(consumers)) consumers=uniqrows(sort(consumers,1))
        postpos=J((Jt ? M : 0),1,.)
        for(j=1;j<=rows(postpos);j++) {
            hit=_cdid_rb_index(consumers:==plan.post_consumer[j]);postpos[j]=hit[1]
        }
        auxpos=J(rows(qids),1,.)
        for(j=1;j<=rows(qids);j++) {
            hit=_cdid_rb_index(consumers:==plan.aux_consumer[qids[j]]);auxpos[j]=hit[1]
        }
        // Retained result + provider + target loading, plus conservative bounded
        // provider/scatter scratch. Never construct originalF-by-allmoments IF.
        retained=cdid_state_cost(buffer.payload)[1]+fit_entries+length(post_loadings)
        perrow=10*(rows(consumers)+2*cols(B.prepared.X)+4)+3*(Jt+rows(qids)+rows(B.request.report))
        if(retained>=cap | perrow>cap-retained) _cdid_rs_fail("Reported IF/provider replay exceeds explicit numeric workspace.",908)
        effective=min((chunk,floor((cap-retained)/perrow)))
        use_ret=(backend!=NULL & rows(ret_if)==B.prepared.population.N & cols(ret_if)==M)
        if(use_ret) use_ret=((*backend).lambda>0 & (*backend).lambda<.)
        if(use_ret) {
            scale_f=1/(*backend).lambda
            aux_in_post=J(rows(qids),1,0);rem_consumers=J(0,1,.);rem_pos=J(rows(qids),1,0)
            for(j=1;j<=rows(qids);j++) {
                hit=_cdid_rb_index(plan.post_consumer:==plan.aux_consumer[qids[j]])
                if(rows(hit)) aux_in_post[j]=hit[1]
                else rem_consumers=rem_consumers\plan.aux_consumer[qids[j]]
            }
            if(rows(rem_consumers)) {
                rem_consumers=uniqrows(sort(rem_consumers,1))
                for(j=1;j<=rows(qids);j++) if(!aux_in_post[j]) {
                    hit=_cdid_rb_index(rem_consumers:==plan.aux_consumer[qids[j]])
                    rem_pos[j]=hit[1]
                }
            }
            for(lo=1;lo<=B.prepared.population.N;lo=lo+effective) {
                hi=min((B.prepared.population.N,lo+effective-1));uid=(lo::hi)
                postIF=(Jt ? (scale_f:*ret_if[uid,])*post_loadings' : J(hi-lo+1,0,.))
                if(rows(qids)) {
                    auxIF=J(hi-lo+1,rows(qids),.)
                    values=J(hi-lo+1,0,.)
                    if(rows(rem_consumers)) {
                        replay=cdid_imp_backend_rows(*backend,uid,rem_consumers,"frame",0)
                        values=replay.phi
                    }
                    for(j=1;j<=rows(qids);j++) {
                        if(aux_in_post[j]) auxIF[,j]=scale_f:*ret_if[uid,aux_in_post[j]]
                        else auxIF[,j]=values[,rem_pos[j]]
                    }
                }
                else auxIF=J(hi-lo+1,0,.)
                cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
            }
        }
        else {
            for(lo=1;lo<=B.prepared.population.N;lo=lo+effective) {
                hi=min((B.prepared.population.N,lo+effective-1));uid=(lo::hi)
                values=J(hi-lo+1,0,.)
                if(rows(consumers)) {
                    replay=cdid_imp_backend_rows(*backend,uid,consumers,"frame",0)
                    values=replay.phi
                }
                postIF=(Jt ? values[,postpos']*post_loadings' : J(hi-lo+1,0,.))
                auxIF=(rows(qids) ? values[,auxpos'] : J(hi-lo+1,0,.))
                cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
            }
        }
    }
    out=cdid_public_report_finish(buffer,limits)
    return(out)
}

struct cdid_public_fit_result scalar cdid_imputation_public_fit(
    struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    /* All borrowed objects are declared here and survive preparation, GMM and report.
       Never return a backend whose model/X/reader belong to a nested local scope. */
    struct cdid_imputation_model scalar model
    struct cdid_population scalar reference_population
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar reference_catalog
    struct cdid_imp_dated_registry scalar reference_registry
    struct cdid_imp_recipe_spec colvector reference_recipes
    struct cdid_imp_dated_reader scalar reference_reader
    real matrix reference_used, reference_keep
    struct cdid_imputation_pair rowvector imp_pairs
    struct cdid_imputation_budget scalar budget
    struct cdid_imp_backend scalar backend
    struct cdid_imp_backend_wctx scalar weight_context
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_gmm_result scalar sol
    real matrix Vw, Lx
    real colvector rk, kept
    real matrix G, loadings, postmeta, Z, Xtrt, links, controls, active_requests, point_table
    real colvector point, diagnostic, reportpost, link_pair_map, points
    string colvector diagnostic_names, link_ids, rec_names
    real scalar M, Jt, used, remaining, i, j, p, fit_entries, external_entries
    real scalar N_pop, p_ctrl, def_perunit, def_chunk, is_def_chunk, avail_entries, eff_chunk, gmm_chunk, gmm_avail, gmm_per_u
    pointer(struct cdid_imp_backend scalar) scalar backend_ref

    if ((B.prepared.design != "static_numeric" & B.prepared.design != "reference_recipes") |
        (B.prepared.method != "ipw" & B.prepared.method != "reg" & B.prepared.method != "dr")) {
        _cdid_pr_fail("cdid_native_private_binding_scope",
            "Pure-imputation public binding requires qualified static/reference ipw, reg, or dr.")
    }
    if (B.settings.gmm_weighting != "optimal" & B.settings.gmm_weighting != "two_step" &
        B.settings.gmm_weighting != "identity" & B.settings.gmm_weighting != "within") {
        _cdid_pr_fail("cdid_chained_public_control_error", "gmmweighting() must be optimal, two_step, identity, or within.")
    }

    B.settings.alpha = ctrl.inference.alpha
    B.settings.compute_if = ctrl.inference.compute_inffunc
    B.settings.bstrap = ctrl.inference.bstrap
    B.settings.cband = ctrl.inference.cband
    B.settings.biters = ctrl.inference.B

    co = cdid_catalog_default_opts()
    co.mem_limit_bytes = ctrl.workspace_bytes
    co.work_cap = ctrl.workcap
    co.conventions = "D10_D13_admissible_frame|IMP-"+B.prepared.method+"-"+B.prepared.imputation_pooling+"-F-actual-change-fixed-w"
    pairs = cdid_cipw_pairs_opts(B.prepared.population, B.request.cohorts, co)
    C = cdid_catalog_build_opts(B.request.calendar, B.request.cohorts, B.request.target_cohorts,
        B.request.anticipation, pairs.occupancy, B.request.control_group, B.request.post, co)
    D = cdid_chained_ipw_build_design(C, B.request.post, B.request.pt_assumption, B.request.pt_start)
    B.request = cdid_public_request_tags(B.request, C)

    plan = cdid_baseline_plan(B.prepared.population, C, B.request.aux,
        B.request.aux_K_ptr, B.request.aux_K_index, ctrl.entries, ctrl.workcap)
    if(B.prepared.design=="reference_recipes") {
        reference_population=B.prepared.population
        reference=B.reference;reference.population=&reference_population
        cdid_ref_activate(reference)
        reference_catalog=cdid_reference_catalog_make(plan,reference,ctrl.entries)
        plan.raw_consumer=reference_catalog.raw_consumer
        plan.post_consumer=reference_catalog.post_consumer;plan.aux_consumer=reference_catalog.aux_consumer
        plan.n_consumers=reference_catalog.n_consumers
    }
    N_pop=B.prepared.population.N
    p_ctrl=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    def_perunit=_cdid_gmm_cells(32,ctrl.links_bound+p_ctrl+1,2^53-1,0)
    def_chunk=max((1,min((N_pop,16384,floor(ctrl.entries/max((1,def_perunit)))))))
    is_def_chunk=(ctrl.chunk==def_chunk)
    backend_ref=NULL
    points=J(plan.n_consumers,1,.)
    fit_entries=plan.workspace_entries
    if(plan.n_consumers) {
        if(B.prepared.design=="reference_recipes") {
            cdid_imp_reference_links(plan,reference_catalog,links,controls,link_ids)
            cdid_imp_reference_registry(reference,reference_registry,reference_recipes)
            reference_keep=cdid_imp_reference_keep(reference)
            cdid_imp_reference_compact(links,reference_recipes,reference_keep)
            reference_used=cdid_imp_used_from_links(reference_population,rows(reference_recipes),links,controls)
            external_entries=plan.workspace_entries+reference.retained_entries+length(links)+length(controls)+length(reference_used)
            for(i=1;i<=rows(reference_registry.z_data);i++)external_entries=external_entries+length(reference_registry.z_data[i].value)
            for(i=1;i<=rows(reference_registry.x_data);i++)external_entries=external_entries+length(reference_registry.x_data[i].value)
            for(i=1;i<=rows(reference_registry.w_data);i++)external_entries=external_entries+length(reference_registry.w_data[i].value)
            if(8*external_entries>=ctrl.workspace_bytes)_cdid_rs_fail("Reference imputation catalogue/design exceeds explicit numeric workspace.",908)
            eff_chunk=ctrl.chunk
            if(is_def_chunk) {
                avail_entries=max((0,floor((ctrl.workspace_bytes-8*external_entries)/8)-8192))
                eff_chunk=max((ctrl.chunk,min((N_pop,16384,floor(avail_entries/max((1,64*(rows(links)+reference_registry.p_raw+reference_registry.q+8))))))))
            }
            budget=cdid_imputation_default_budget(eff_chunk,ctrl.baseline.max_cached,ctrl.workspace_bytes-8*external_entries,ctrl.workcap)
            if(rows(reference_keep)) reference_reader=cdid_imp_dated_reader_make(&reference_population,&reference_registry,reference_recipes,reference_used,B.prepared.method!="reg",budget,reference_keep)
            else reference_reader=cdid_imp_dated_reader_make(&reference_population,&reference_registry,reference_recipes,reference_used,B.prepared.method!="reg",budget)
            model=cdid_imp_fit_dated(reference_reader,links,controls,B.prepared.imputation_pooling,ctrl.baseline.rank_tol,budget)
            backend=cdid_imp_backend_prepare_d(&model,&reference_reader,links,controls,link_ids,B.prepared.method,
                ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,ctrl.chunk,budget,"public_reference_imputation",
                B.prepared.outcome_pooling)
        }
        else {
        cdid_imputation_public_links(plan,links,controls,link_ids)
        cdid_imp_links_to_requests(B.request.calendar,links,controls,imp_pairs,active_requests,link_pair_map)
        p=cols(B.prepared.X)
        Z=(B.prepared.x_has_intercept ? (p>1 ? B.prepared.X[,2..p] : J(B.prepared.population.N,0,.)) : B.prepared.X)
        Xtrt=(B.prepared.method=="reg" ? J(B.prepared.population.N,0,.) : B.prepared.X)
        external_entries=plan.workspace_entries+length(Z)+length(Xtrt)+length(links)+length(controls)+length(active_requests)+length(link_pair_map)
        if(8*external_entries>=ctrl.workspace_bytes)
            _cdid_rs_fail("Imputation catalogue/design exceeds explicit numeric workspace.",908)
        eff_chunk=ctrl.chunk
        if(is_def_chunk) {
            avail_entries=max((0,floor((ctrl.workspace_bytes-8*external_entries)/8)-8192))
            eff_chunk=max((ctrl.chunk,min((N_pop,16384,floor(avail_entries/max((1,64*(rows(links)+p_ctrl+8))))))))
        }
        budget=cdid_imputation_default_budget(eff_chunk,ctrl.baseline.max_cached,
            ctrl.workspace_bytes-8*external_entries,ctrl.workcap)
        model=cdid_imputation_fit_static(B.prepared.population,B.prepared.weight.owner_weight,Z,
            imp_pairs,active_requests,B.prepared.imputation_pooling,B.prepared.x_has_intercept,
            ctrl.baseline.rank_tol,budget)
        backend=cdid_imp_backend_prepare(&model,&Xtrt,links,controls,link_ids,B.prepared.method,
            ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,ctrl.chunk,budget,"public_static_imputation",
            B.prepared.outcome_pooling)
        }
        backend_ref=&backend
        point_table=cdid_imp_backend_points(backend)
        points=point_table[,1]
        fit_entries=external_entries+ceil((backend.retained_bytes+backend.scratch_bytes)/8)+length(points)
        if(plan.M) if(missing(points[plan.post_consumer]))
            _cdid_pr_fail("cdid_native_reporting_input","Imputation post catalogue has a nonfinite fitted point.")
        for(i=1;i<=plan.Q;i++) if(plan.aux_supported[i]) {
            if(points[plan.aux_consumer[i]]>=.) {
                plan.aux_supported[i]=0
                plan.aux_reason[i]="nonfinite_point"
            }
        }
    }
    sys = cdid_baseline_make_system(C, D, B.request.post, ctrl.entries, fit_entries)
    support = cdid_public_post_ledger(C, D, B.request, sys)
    M = sys.M
    Jt = rows(B.request.post)
    point = J(Jt, 1, .)
    loadings = J(0, M, .)
    diagnostic = J(0, 1, .)
    diagnostic_names = J(0, 1, "")

    gmm_chunk = ctrl.chunk
    if (M > 0 & sys.p_retained > 0) {
        if (is_def_chunk) {
            gmm_avail = (ctrl.entries - fit_entries - sys.workspace_entries - 4 * M * M -
                8 * (M + sys.p_retained + sys.n_estimable + 16)^2)
            gmm_per_u = 16 * M + 8 * (cols(B.prepared.X) + cols(B.prepared.Z_sel) + 4)
            if (gmm_avail > 0) {
                gmm_chunk = max((ctrl.chunk, min((N_pop, 16384, floor(gmm_avail / max((1, gmm_per_u)))))))
            }
        }
        used = fit_entries + sys.workspace_entries
        used = _cdid_base_cells(8 * gmm_chunk,
            M + cols(B.prepared.X) + cols(B.prepared.Z_sel) + 4, ctrl.entries, used)
        if (B.settings.gmm_weighting == "identity") {
            used = _cdid_gmm_identity_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit_entries - sys.workspace_entries
            sol = cdid_gmm_identity_map(points[plan.post_consumer], sys.full_A, sys.H_est, remaining)
        }
        else {
            used = _cdid_base_cells(2 * gmm_chunk, M, ctrl.entries, used)
            used = _cdid_base_cells(3 * M, M, ctrl.entries, used)
            used = _cdid_gmm_workspace(M, sys.p_retained, sys.n_estimable, ctrl.entries, used)
            remaining = ctrl.entries - fit_entries - sys.workspace_entries - M * M
            /* Original scientific moment order; full-M metric and ridge scale.
               R cp04_complete_reduction returns NULL for cdid_imputation_backend. */
            weight_context=cdid_imp_backend_wctx_make(&backend,(plan.post_consumer,J(M,1,1)))
            if (B.settings.gmm_weighting == "within") {
                rk = points[plan.post_consumer]
                if (B.prepared.design == "reference_recipes")
                    cdid_wsr_within(B, plan, C, rk, 0, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit,
                        gmm_chunk, ctrl.baseline.rank_tol, ctrl.baseline.max_cached, ctrl.workspace_bytes,
                        ctrl.workcap, kept, G, Vw, reference)
                else cdid_wsr_within(B, plan, C, rk, 0, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit,
                    gmm_chunk, ctrl.baseline.rank_tol, ctrl.baseline.max_cached, ctrl.workspace_bytes,
                    ctrl.workcap, kept, G, Vw)
                sol = cdid_gmm_fit_within_norm(rk[kept], sys.full_A[kept, ], Vw, G, M, sys.H_est,
                    ctrl.condition_trigger, ctrl.ridge_relative,
                    B.prepared.population.N, B.prepared.population.n_A, gmm_chunk,
                    &cdid_imp_backend_weight(), &weight_context, remaining, M, kept')
                Lx = J(rows(sol.target_loadings), M, 0)
                Lx[., kept'] = sol.target_loadings
                sol.target_loadings = Lx
                Lx = Vw = J(0, 0, .)
            }
            else {
                kept = (1::M)
                G = I(M)
                sol = cdid_gmm_fit_unit_norm(points[plan.post_consumer], sys.full_A, G, M, sys.H_est,
                    B.settings.gmm_weighting, ctrl.condition_trigger, ctrl.ridge_relative,
                    B.prepared.population.N, B.prepared.population.n_A, gmm_chunk,
                    &cdid_imp_backend_weight(), &weight_context, remaining)
            }
        }
        if (B.settings.gmm_weighting == "identity") kept = (1::M)
        if (sys.n_estimable) point[sys.estimable_ids] = sol.targets
        loadings = sol.target_loadings
        cdid_chain_sp_record(B, C, D, kept, points[plan.post_consumer], loadings, sys.estimable_ids,
            cdid_chain_auxpts(points, plan.aux_consumer, plan.aux_supported), plan.aux_supported)
        diagnostic = (sol.criterion \ sol.lambda \ sol.rho \ sol.trace_B \
            sol.covariance_rank \ sol.source_n_units \ sol.normalization_factor \ sol.sample_rank_upper)
        diagnostic_names = ("GMM_criterion" \ "GMM_lambda" \ "GMM_ridge_rho" \ "GMM_trace_B" \
            "GMM_covariance_rank" \ "GMM_source_replay_rows" \ "GMM_gram_normalization_factor" \ "GMM_sample_rank_upper")
        G = J(0, 0, .)
        sol.unit_gram_covariance = sol.coefficient_loadings = sol.target_loadings = J(0, 0, .)
        sol.coefficient_covariance = sol.covariance = J(0, 0, .)
    }
    else if (sys.n_estimable) {
        _cdid_pr_fail("cdid_native_reporting_input", "An estimable post target has no GMM parameter system.")
    }

    sys.full_A = sys.H_est = J(0, 0, .)
    out.payload = cdid_imputation_public_report(B, plan, backend_ref, points, fit_entries, point, sys.post_available,
        support.status, support.reason, sys.estimable_ids, loadings,
        ctrl.inference.compute_inffunc, ctrl.inference.limits, gmm_chunk, ctrl.entries, sol.retained_if)
    sol.retained_if = J(0, 0, .)
    postmeta = J(rows(B.request.report), cols(support.ledger), .)
    reportpost = _cdid_rb_index(B.request.report[, 5] :== 3)
    for (i = 1; i <= rows(reportpost); i++) {
        j = reportpost[i]
        postmeta[j, ] = support.ledger[B.request.report[j, 6], ]
    }
    out.payload.ledger = out.payload.ledger, postmeta
    out.payload.ledger_columns = out.payload.ledger_columns, "post_support." :+ support.columns
    out.payload.numeric_settings = out.payload.numeric_settings \ diagnostic \
        ctrl.workspace_bytes \ ctrl.entries \ ctrl.chunk \ ctrl.workcap \
        ctrl.condition_trigger \ ctrl.ridge_relative \ ctrl.links_bound \ ctrl.J_bound
    out.payload.numeric_setting_names = out.payload.numeric_setting_names \ diagnostic_names \
        "numeric_workspace_bytes" \ "numeric_workspace_entries" \ "chunk_rows" \
        "counted_work_limit" \ "condition_trigger" \ "ridge_relative" \
        "catalogue_links_bound" \ "inference_targets_bound"
    cdid_imputation_public_metadata(out.payload,plan,backend_ref,C)
    if(B.prepared.design=="reference_recipes") {
        for(i=1;i<=rows(out.payload.text_settings);i++)if(out.payload.text_settings[i,1]=="profile_kind")out.payload.text_settings[i,2]="reference_recipes"
        rec_names=J(rows(plan.raw),1,"")
        for(i=1;i<=rows(plan.raw);i++) rec_names[i]="source_link_recipe:"+_cdid_rb_key(i)
        out.payload.numeric_settings=out.payload.numeric_settings\reference.n_recipes\reference_catalog.n_blocks\reference.view_recipe
        out.payload.numeric_setting_names=out.payload.numeric_setting_names\"n_reference_recipes"\"n_reference_OR_blocks"\rec_names
    }
    out.payload.text_settings = out.payload.text_settings \
        ("GMM_metric", "full original catalogue identity metric before covariance weighting") \
        ("registry_limit_meaning", "2^53-1 exact-count ceiling; not a session RAM/work budget")
    out.payload.source_versions = out.payload.source_versions \
        ("gmm", "GMM combination, normalized by the units with two observed outcomes; change-imputation comparisons")
    out.inference = ctrl.inference
    cdid_state_validate(out.payload, ctrl.inference.limits)
    return(out)
}

real scalar _cdid_mat_varies_by_unit(real colvector row_unit, real scalar N, real matrix M)
{
    real scalar R, K, j, i, u
    real colvector first_u, idx_map, u_val, col_j
    R = rows(M)
    K = cols(M)
    if (R <= 1 | K == 0 | N < 1) return(0)
    if (!missing(M)) {
        first_u = J(N, 1, 0)
        for (i=R; i>=1; i--) first_u[row_unit[i]] = i
        idx_map = first_u[row_unit]
        for (j=1; j<=K; j++) {
            col_j = M[, j]
            if (any(col_j :!= col_j[idx_map])) return(1)
        }
        return(0)
    }
    for (j=1; j<=K; j++) {
        u_val = J(N, 1, .)
        for (i=1; i<=R; i++) {
            if (M[i, j] >= .) continue
            u = row_unit[i]
            if (u_val[u] >= .) u_val[u] = M[i, j]
            else if (M[i, j] != u_val[u]) return(1)
        }
    }
    return(0)
}

void _cdid_try_static_prepare_store(string scalar input_name, string scalar output_name,
    string scalar method, real scalar noconstant, string scalar fixweights,
    real colvector targets, string scalar control, real scalar ant, string scalar base,
    string scalar pt, real scalar ptstart, real scalar max_entries,
    string scalar selection_in, string scalar selection_model_in,
    string scalar selection_pooling_in, string scalar selection_universe_in,
    string scalar selection_justification_in, string scalar selection_est_method_in,
    string scalar imputation_pooling_in, string scalar gmm_covariance_in,
    real scalar selnoconstant_in, real scalar selconstant_in, real scalar selzero_in,
    string scalar outcome_pooling_in)
{
    struct cdid_ci_input scalar I
    struct cdid_public_bundle scalar B
    struct cdid_population_ids scalar ids
    pointer() scalar p_in, p_out
    real matrix M_chk
    real colvector same_id
    real scalar R, K, N_est, needs_ref, j, q_x
    p_in = findexternal(input_name)
    if (p_in == NULL) _cdid_pr_fail("cdid_chained_request_input_error", "Transient input is absent.")
    if (findexternal(output_name) != NULL)
        _cdid_pr_fail("cdid_chained_request_input_error", "Requested output name already exists; no user object is overwritten.")
    I = *p_in
    R = rows(I.date)
    q_x = (cols(I.frame_X) > 0 ? cols(I.frame_X) : cols(I.X))
    K = cols(I.X) + cols(I.Q) + (rows(I.w) > 0 ? 1 : 0)
    needs_ref = 0
    if (max_entries < 1048576 | (q_x == 0 & noconstant)) {
        needs_ref = 1
    }
    else if (I.string_id) {
        if (any(I.string_ids :== "") | (I.frame_supplied ? any(I.frame_string_ids :== "") : 0)) {
            needs_ref = 1
        }
    }
    else {
        if (missing(I.id) | (I.frame_supplied ? missing(I.frame_id) : 0)) {
            needs_ref = 1
        }
    }
    if (!needs_ref & I.frame_supplied) {
        if ((cols(I.X) > 0 & cols(I.frame_X) == 0) |
            (cols(I.Q) > 0 & cols(I.frame_Q) == 0) |
            (rows(I.w) > 0 & rows(I.frame_w) == 0)) {
            needs_ref = 1
        }
        else if (cols(I.X) > 0 ? (cols(I.X) != cols(I.frame_X) ? 1 : any(I.x_names :!= I.frame_x_names)) : 0) {
            needs_ref = 1
        }
        else if (cols(I.Q) > 0 ? (cols(I.Q) != cols(I.frame_Q) ? 1 : any(I.q_names :!= I.frame_q_names)) : 0) {
            needs_ref = 1
        }
        else if (missing(I.frame_X) | missing(I.frame_Q) |
            (rows(I.frame_w) > 0 ? (missing(I.frame_w) | any(I.frame_w :<= 0)) : 0)) {
            needs_ref = 1
        }
    }
    if (!needs_ref & K > 0 & R >= 1) {
        M_chk = J(R, 0, .)
        if (cols(I.X) > 0) M_chk = (M_chk, I.X)
        if (cols(I.Q) > 0) M_chk = (M_chk, I.Q)
        if (rows(I.w) > 0) M_chk = (M_chk, I.w)
        if (missing(M_chk) | (rows(I.w) > 0 ? any(I.w :<= 0) : 0)) {
            needs_ref = 1
        }
        else if (R >= 2) {
            if (!I.string_id) {
                if (all(I.id[2..R] :>= I.id[1..(R - 1)])) {
                    same_id = (I.id[2..R] :== I.id[1..(R - 1)])
                    for (j=1; j<=K; j++) {
                        if (any(same_id :& (M_chk[2..R, j] :!= M_chk[1..(R - 1), j]))) {
                            needs_ref = 1
                            break
                        }
                    }
                }
                else {
                    ids = cdid_population_ids_numeric(I.id, J(0, 1, .), 0)
                    N_est = rows(ids.numeric_id)
                    needs_ref = _cdid_mat_varies_by_unit(ids.row_unit, N_est, M_chk)
                }
            }
            else {
                ids = cdid_population_ids_string(I.string_ids, J(0, 1, ""), 0)
                N_est = rows(ids.string_id_value)
                needs_ref = _cdid_mat_varies_by_unit(ids.row_unit, N_est, M_chk)
            }
        }
    }
    if (needs_ref) {
        I.X = J(R, 0, .)
        I.x_names = J(1, 0, "")
        I.w = J(0, 1, .)
        I.Q = J(0, 0, .)
        I.q_names = J(1, 0, "")
        I.frame_X = J(0, 0, .)
        I.frame_x_names = J(1, 0, "")
        I.frame_w = J(0, 1, .)
        I.frame_Q = J(0, 0, .)
        I.frame_q_names = J(1, 0, "")
        *p_in = I
        _cdid_public_prepare_store(input_name, output_name, method, noconstant, fixweights,
            targets, control, ant, base, pt, ptstart, max_entries,
            selection_in, selection_model_in, selection_pooling_in, selection_universe_in,
            selection_justification_in, selection_est_method_in, imputation_pooling_in,
            gmm_covariance_in, selnoconstant_in, selconstant_in, selzero_in, outcome_pooling_in)
        st_local("fast_static_done", "0")
        return
    }
    B.prepared = cdid_public_static_prepare(I, method, noconstant, fixweights,
        selection_in, selection_model_in, selection_pooling_in, selection_universe_in,
        selection_justification_in, selection_est_method_in, imputation_pooling_in,
        gmm_covariance_in, selnoconstant_in, selconstant_in, selzero_in, outcome_pooling_in, 0)
    B.request = cdid_public_request_build(B.prepared.population, targets, control, ant, base, pt, ptstart, max_entries)
    B.settings = cdid_public_defaults()
    B.settings.method = method
    B.settings.control_group = control
    B.settings.anticipation = ant
    B.settings.base_period = base
    B.settings.pt_assumption = pt
    B.settings.pt_start = ptstart
    B.settings.selection = B.prepared.selection
    B.settings.selection_canonical = B.prepared.selection_canonical
    B.settings.selection_model = B.prepared.selection_model
    B.settings.selection_pooling = B.prepared.selection_pooling
    B.settings.selection_universe = B.prepared.selection_universe
    B.settings.selection_est_method = B.prepared.selection_est_method
    B.settings.selection_justification = B.prepared.selection_justification
    B.settings.imputation_pooling = B.prepared.imputation_pooling
    B.settings.outcome_pooling = B.prepared.outcome_pooling
    B.settings.inactive_arguments = B.prepared.inactive_arguments
    B.settings.gmm_covariance = B.prepared.gmm_covariance
    p_out = crexternal(output_name)
    if (p_out == NULL) _cdid_pr_fail("cdid_chained_request_input_error", "Cannot create requested output object.")
    *p_out = B
    st_local("pr_N", strofreal(B.prepared.population.N))
    st_local("pr_A", strofreal(B.prepared.population.n_A))
    st_local("pr_post", strofreal(rows(B.request.post)))
    st_local("pr_aux", strofreal(rows(B.request.aux)))
    st_local("pr_report", strofreal(rows(B.request.report)))
    st_local("fast_static_done", "1")
}

end
