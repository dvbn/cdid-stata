*! Native public augmented observation correction, static and dated recipes; StataNow19.5.
*! GPL-3 adaptation of R cdid chained_augmented_public.R.
version 19.5
mata:
void _cdid_aug_public_gate(string scalar method,string scalar selection,string scalar universe,
    string scalar model,string scalar pooling,real scalar has_rho,string scalar justification)
{
    (void)cdid_aug_route_gate(method,0,selection,universe,model,pooling,has_rho,justification)
}
void _cdid_aug_reference_map(string scalar name,string scalar xf,string scalar fxf,
    string scalar qf,string scalar fqf,string scalar w,string scalar fw)
{
    string rowvector x,f
    string matrix map
    real scalar role,j,k
    pointer() scalar p
    map=J(0,2,"")
    for(role=1;role<=2;role++) {
        x=cdid_ref_formula_names(role==1 ? xf : qf)
        f=cdid_ref_formula_names(role==1 ? fxf : fqf)
        if(cols(f)) {
            if(cols(f)!=cols(x))_cdid_ref_fail("Data/frame formula mappings disagree in dimension.",459)
            for(j=1;j<=cols(x);j++) {
                k=0
                if(rows(map)) {
                    if(any(map[,1]:==x[j])) {
                        k=selectindex(map[,1]:==x[j])[1]
                        if(map[k,2]!=f[j])_cdid_ref_fail("Formula roles disagree on a raw frame mapping.",459)
                    }
                }
                if(!k)map=map\(x[j],f[j])
            }
        }
    }
    if(fw!="")map=map\(w,fw)
    p=crexternal(name);if(p==NULL)_cdid_ref_fail("Transient frame map name is in use.",498)
    *p=map
}

/* Shared canonical Q date slots. One per unique role/date, not per recipe;
   retain raw missing values outside the shared required-key union. */
void cdid_aug_reference_qslots(struct cdid_reference_inputs scalar S,
    struct cdid_imp_recipe_spec colvector recipes,
    struct cdid_imp_matrix_slot colvector qdata, real colvector recipe_qslot)
{
    real colvector keys,dates,hit,units
    real matrix value
    real scalar slot,j,d,r
    slot=S.roles[3].slot
    if(!S.roles[3].active | slot<1 | !S.designs[slot].encoded)
        _cdid_ref_fail("Augmentation requires the active encoded Q design.",498)
    keys=S.designs[slot].keys
    dates=uniqrows(sort(floor(keys:/(S.N+1)),1))
    qdata=cdid_imp_matrix_slot(rows(dates),1)
    for(j=1;j<=rows(dates);j++) {
        d=dates[j];hit=_cdid_ref_ix(floor(keys:/(S.N+1)):==d)
        units=keys[hit]:-d*(S.N+1)
        value=J(S.N,cols(S.designs[slot].value),.)
        if(rows(units))value[units,]=cdid_ref_design_rows(S,3,keys[hit])
        qdata[j]=cdid_imp_dated_slot_make(d,(d ? S.calendar[d] : .),value)
    }
    recipe_qslot=J(rows(recipes),1,0)
    for(r=1;r<=rows(recipes);r++) {
        d=(S.roles[3].dated ? recipes[r].design_cal : 0)
        hit=_cdid_ref_ix(dates:==d)
        if(rows(hit)!=1)_cdid_ref_fail("A used augmented design block has no canonical Q date slot.",498)
        recipe_qslot[r]=hit[1]
    }
}

void cdid_augmented_public_metadata(struct cdid_result_state scalar out,
    struct cdid_baseline_plan_state scalar plan,pointer(struct cdid_aug_backend scalar) scalar backend,
    struct cdid_cat scalar C)
{
    real scalar i,c,k,ni,no,np,nb,nq,ncells,ngroups,nlocal,p_or,pscols,link_point,n_raw,b6,b5
    real colvector K,add_num
    string colvector coh_str,add_names
    string matrix add_txt
    string scalar key,control_key,id,source_key,actual_key,reason,act_prefix
    nq=ni=no=np=nb=ncells=ngroups=nlocal=p_or=pscols=0
    if(backend!=NULL) {
        ni=(*backend).n_imputation_fits;no=(*backend).n_or_fits;np=(*backend).n_ps_fits
        nb=(*backend).n_blocks;nq=(*backend).n_q_fits
        ncells=(*backend).n_cells;ngroups=(*backend).n_groups
        nlocal=(*backend).intercept*ncells;p_or=(*backend).or_p
        pscols=(*backend).q
    }
    act_prefix="cmpcat:AUG-"+strupper(_cdid_rb_stored_text(out,"method"))+":"+_cdid_rb_stored_text(out,"imputation_pooling")+":"
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
    out.numeric_settings=out.numeric_settings\plan.n_F\plan.n_A\nq\nq\np\no\ni\ncells\ngroups\nlocal\nb\p_or\pscols\plan.M\plan.Q\plan.n_consumers\add_num
    out.numeric_setting_names=out.numeric_setting_names\"GMM_replay_rows"\"GMM_normalization_count"\"n_selection_fits"\"observation_logit_models"\
        "n_treatment_fits"\"n_pooled_OR_fits"\"n_imputation_fits"\"n_imputation_cells"\"n_imputation_groups"\"n_imputation_local_intercepts"\
        "n_OR_blocks"\"p_OR_retained"\"p_treatment"\"M_full"\"Q_aux_union"\"n_unique_consumers"\add_names
    out.text_settings=out.text_settings\("GMM_weight_row_domain","all original F; core psi=(n_A/N_F)*phi_F; covariance normalization n_A")\
        ("profile_kind","static_numeric")\("provider_status","ok")\
        ("imputation_prediction_universe","every required cohort row in original F")\
        ("observation_model","Z=j+S*(Delta-j)/max(q,1e-6); raw-W q, W*a PS/OR")\
        ("duration_convention","stored-calendar index difference")\add_txt
    out.source_versions=out.source_versions\("estimator","augmented correction for missing outcomes: change imputation, observation probabilities, treatment propensities and outcome regression; influence functions over all units")
}

struct cdid_result_state scalar cdid_augmented_public_report(
    struct cdid_public_bundle scalar B,struct cdid_baseline_plan_state scalar plan,
    pointer(struct cdid_aug_backend scalar) scalar backend,real colvector points,
    real scalar fit_entries,
    real colvector post_point,real colvector post_available,string colvector post_status,
    string colvector post_reason,real colvector post_est_ids,real matrix post_loadings,
    real scalar compute_if,struct cdid_state_limits scalar limits,real scalar chunk,real scalar cap,
    | real matrix retained_if)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    struct cdid_aug_backend_if scalar replay
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
                        replay=cdid_aug_backend_rows(*backend,uid,rem_consumers,"frame",0)
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
                    replay=cdid_aug_backend_rows(*backend,uid,consumers,"frame",0)
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

struct cdid_public_fit_result scalar cdid_augmented_public_fit(
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
    struct cdid_imp_backend scalar ipw_backend
    struct cdid_aug_backend scalar backend
    struct cdid_sel_reader scalar qreader
    struct cdid_sel_decl scalar qdecl
    struct cdid_imp_matrix_slot colvector qdata
    real colvector recipe_qslot
    struct cdid_aug_backend_wctx scalar weight_context
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_gmm_result scalar sol
    real matrix Vw, Lx
    real colvector rk, kept
    real matrix G, loadings, postmeta, Z, Xtrt, links, controls, active_requests, point_table
    real colvector point, diagnostic, reportpost, link_pair_map, points
    string colvector diagnostic_names, link_ids, rec_names
    real scalar M, Jt, used, remaining, i, j, p, fit_entries, external_entries
    real scalar N_pop, p_ctrl, def_perunit, def_chunk, is_def_chunk, avail_entries, eff_chunk, gmm_chunk, gmm_avail, gmm_per_u, q_cols
    pointer(struct cdid_aug_backend scalar) scalar backend_ref

    if ((B.prepared.design != "static_numeric" & B.prepared.design != "reference_recipes") |
        (B.prepared.method != "ipw" & B.prepared.method != "reg" & B.prepared.method != "dr") |
        B.prepared.selection_est_method != "dr") {
        _cdid_pr_fail("cdid_native_private_binding_scope",
            "Augmented public binding requires method ipw, reg, or dr with selection_est_method dr.")
    }
    if (B.settings.gmm_weighting != "optimal" & B.settings.gmm_weighting != "two_step" &
        B.settings.gmm_weighting != "identity" & B.settings.gmm_weighting != "within") {
        _cdid_pr_fail("cdid_chained_public_control_error", "gmmweighting() must be optimal, two_step, identity, or within.")
    }

    qdecl=cdid_aug_route_gate(B.prepared.method,0,B.prepared.selection,B.prepared.selection_universe,
        B.prepared.selection_model,B.prepared.selection_pooling,B.prepared.sel_has_rho,B.prepared.selection_justification)
    B.settings.alpha = ctrl.inference.alpha
    B.settings.compute_if = ctrl.inference.compute_inffunc
    B.settings.bstrap = ctrl.inference.bstrap
    B.settings.cband = ctrl.inference.cband
    B.settings.biters = ctrl.inference.B

    co = cdid_catalog_default_opts()
    co.mem_limit_bytes = ctrl.workspace_bytes
    co.work_cap = ctrl.workcap
    co.conventions = "D10_D13_admissible_frame|AUG-"+B.prepared.method+"-"+B.prepared.imputation_pooling+"-F-actual-change-fixed-w"
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
            if(8*external_entries>=ctrl.workspace_bytes)_cdid_rs_fail("Reference augmented catalogue/design exceeds explicit numeric workspace.",908)
            eff_chunk=ctrl.chunk
            if(is_def_chunk) {
                q_cols=(reference.roles[3].active & reference.roles[3].slot>=1 ? cols(reference.designs[reference.roles[3].slot].value) : 0)
                avail_entries=max((0,floor((ctrl.workspace_bytes-8*external_entries)/8)-8192))
                eff_chunk=max((ctrl.chunk,min((N_pop,16384,floor(avail_entries/max((1,64*(rows(links)+reference_registry.p_raw+reference_registry.q+q_cols+8))))))))
            }
            budget=cdid_imputation_default_budget(eff_chunk,ctrl.baseline.max_cached,ctrl.workspace_bytes-8*external_entries,ctrl.workcap)
            if(rows(reference_keep)) reference_reader=cdid_imp_dated_reader_make(&reference_population,&reference_registry,reference_recipes,reference_used,B.prepared.method!="reg",budget,reference_keep)
            else reference_reader=cdid_imp_dated_reader_make(&reference_population,&reference_registry,reference_recipes,reference_used,B.prepared.method!="reg",budget)
            model=cdid_imp_fit_dated(reference_reader,links,controls,B.prepared.imputation_pooling,ctrl.baseline.rank_tol,budget)
            ipw_backend=cdid_imp_backend_prepare_d(&model,&reference_reader,links,controls,link_ids,
                (B.prepared.method=="reg" ? "reg" : "ipw"),
                ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,ctrl.chunk,budget,
                "public_reference_augmented_factors","",(B.prepared.method=="reg"))
            cdid_aug_reference_qslots(reference,reference_recipes,qdata,recipe_qslot)
            backend=cdid_aug_backend_prepare_d(&ipw_backend,&reference_reader,qdata,recipe_qslot,
                reference.roles[3].intercept,(reference.roles[3].intercept ? 1 : 0),qdecl,
                ctrl.chunk,budget,"public_reference_augmented",
                B.prepared.method,B.prepared.outcome_pooling)
            qdata=cdid_imp_matrix_slot(0,1)
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
            eff_chunk=max((ctrl.chunk,min((N_pop,16384,floor(avail_entries/max((1,64*(rows(links)+p_ctrl+cols(B.prepared.Z_sel)+8))))))))
        }
        budget=cdid_imputation_default_budget(eff_chunk,ctrl.baseline.max_cached,
            ctrl.workspace_bytes-8*external_entries,ctrl.workcap)
        model=cdid_imputation_fit_static(B.prepared.population,B.prepared.weight.owner_weight,Z,
            imp_pairs,active_requests,B.prepared.imputation_pooling,B.prepared.x_has_intercept,
            ctrl.baseline.rank_tol,budget)
        ipw_backend=cdid_imp_backend_prepare(&model,&Xtrt,links,controls,link_ids,
            (B.prepared.method=="reg" ? "reg" : "ipw"),
            ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,ctrl.chunk,budget,
            "public_static_augmented_factors","",(B.prepared.method=="reg"))
        qreader=cdid_sel_reader_pop(B.prepared.population,B.prepared.Z_sel,B.prepared.sel_has_intercept,
            (B.prepared.sel_has_intercept ? 1 : 0),B.prepared.weight.owner_weight,qdecl,
            J(B.prepared.population.N,1,1),eff_chunk)
        backend=cdid_aug_backend_prepare(&ipw_backend,&qreader,ctrl.chunk,budget,
            "public_static_augmented",B.prepared.method,B.prepared.outcome_pooling)
        }
        backend_ref=&backend
        if(cdid_nuis_requested())cdid_nuis_capture_augmented(backend,B,plan)
        points=cdid_aug_backend_points(backend)
        fit_entries=external_entries+ceil((backend.retained_bytes+backend.scratch_bytes)/8)+length(points)
        if(plan.M) if(missing(points[plan.post_consumer]))
            _cdid_pr_fail("cdid_native_reporting_input","Augmented post catalogue has a nonfinite fitted point.")
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
            weight_context=cdid_aug_backend_wctx_make(&backend,(plan.post_consumer,J(M,1,1)))
            if (B.settings.gmm_weighting == "within") {
                rk = points[plan.post_consumer]
                if (B.prepared.design == "reference_recipes")
                    cdid_wsr_within(B, plan, C, rk, 1, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit,
                        gmm_chunk, ctrl.baseline.rank_tol, ctrl.baseline.max_cached, ctrl.workspace_bytes,
                        ctrl.workcap, kept, G, Vw, reference)
                else cdid_wsr_within(B, plan, C, rk, 1, ctrl.baseline.glm_epsilon, ctrl.baseline.glm_maxit,
                    gmm_chunk, ctrl.baseline.rank_tol, ctrl.baseline.max_cached, ctrl.workspace_bytes,
                    ctrl.workcap, kept, G, Vw)
                sol = cdid_gmm_fit_within_norm(rk[kept], sys.full_A[kept, ], Vw, G, M, sys.H_est,
                    ctrl.condition_trigger, ctrl.ridge_relative,
                    B.prepared.population.N, B.prepared.population.n_A, gmm_chunk,
                    &cdid_aug_backend_weight(), &weight_context, remaining, M, kept')
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
                    &cdid_aug_backend_weight(), &weight_context, remaining)
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
    out.payload = cdid_augmented_public_report(B, plan, backend_ref, points, fit_entries, point, sys.post_available,
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
    cdid_augmented_public_metadata(out.payload,plan,backend_ref,C)
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
        ("gmm", "GMM combination, normalized by the units with two observed outcomes; augmented doubly robust comparisons")
    out.inference = ctrl.inference
    cdid_state_validate(out.payload, ctrl.inference.limits)
    return(out)
}

end
