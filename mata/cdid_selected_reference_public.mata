*! Public conditional-MTAR observation-weighted reference consumer; StataNow19.5.
version 19.5
mata:
void _cdid_sr_prepare_store(string scalar input,string scalar output,string scalar method,
 real scalar nc,string scalar fixweights,real colvector targets,string scalar control,
 real scalar ant,string scalar base,string scalar pt,real scalar start,real scalar cap,
 string scalar selection,string scalar model,string scalar pooling,string scalar universe,
 string scalar justification,string scalar covariance,real scalar si,real scalar inherited,
 | string scalar outcome_pooling_in)
{
 struct cdid_ci_input scalar I
 struct cdid_public_bundle scalar B
 struct cdid_sel_decl scalar decl
 struct cdid_population_ids scalar ids
 pointer() scalar p
 string scalar outpool
 real scalar out_active
 p=findexternal(input);if(p==NULL)_cdid_pr_fail("cdid_chained_request_input_error","Transient input is absent.")
 if(findexternal(output)!=NULL)_cdid_pr_fail("cdid_chained_request_input_error","Requested output name already exists; no user object is overwritten.")
 I=*p
 if(method!="ipw" & method!="reg" & method!="dr")_cdid_pr_fail("cdid_chained_request_input_error","method() must be dr, reg, ipw, or custom.")
 out_active=(method=="reg" | method=="dr")
 if(out_active) {
  outpool=(args()>=21 ? strlower(strtrim(outcome_pooling_in)) : "")
  if(outpool=="")outpool="separate"
  if(outpool!="separate" & outpool!="pooled")_cdid_prep_role_fail("cdid_chained_request_outcome_pooling","outcomepooling() must be separate or pooled.",198)
 }
 else outpool=""
 covariance=strlower(strtrim(covariance));if(covariance=="")covariance="tikhonov"
 if(covariance!="tikhonov")_cdid_prep_role_fail("cdid_chained_request_input_error","gmmcovariance() must be tikhonov.",198)
 decl=cdid_sel_norm_decl(selection,1,universe,I.frame_has_rho|I.data_has_rho,strtrim(justification)!="",model,pooling,justification)
 if(decl.canonical_assumption!="conditional_mtar")_cdid_pr_fail("cdid_chained_request_input_error","Selected reference preparation requires conditional_mtar.")
 if(I.string_id)ids=cdid_population_ids_string(I.string_ids,I.frame_string_ids,I.frame_supplied)
 else ids=cdid_population_ids_numeric(I.id,I.frame_id,I.frame_supplied)
 B.prepared.population=cdid_population_prepare(ids,I.date,I.y,I.g,I.frame_g,I.frame_has_g)
 if(any(B.prepared.population.cohort_known:==0))_cdid_pr_fail("cdid_chained_request_cohort_error","Corollary1 requires known cohorts on every F row, including frame-only rows.")
 B.prepared.method=method;B.prepared.design="unencoded";B.prepared.source_row=I.source_row
 B.prepared.required_X=B.prepared.population.A;B.prepared.x_has_intercept=!nc
 B.prepared.selection=selection;B.prepared.selection_canonical=decl.canonical_assumption
 B.prepared.selection_model=decl.selection_model;B.prepared.selection_pooling=decl.selection_pooling
 B.prepared.selection_universe=decl.selection_universe;B.prepared.selection_justification=decl.justification
 B.prepared.selection_est_method="ipw";B.prepared.inactive_arguments=(out_active ? "imputation_pooling" : "outcome_pooling imputation_pooling")
 B.prepared.gmm_covariance=covariance;B.prepared.imputation_pooling="";B.prepared.outcome_pooling=outpool
 B.prepared.U=(decl.selection_universe=="F" ? J(B.prepared.population.N,1,1) :
  (decl.selection_universe=="B" ? B.prepared.population.B : B.prepared.population.A))
 if(sum(B.prepared.U)==0)_cdid_pr_fail("cdid_chained_request_input_error","Selected universe U has zero units.")
 B.prepared.sel_has_intercept=si;B.prepared.selection_formula_inherited=inherited
 B.prepared.sel_has_rho=I.frame_has_rho|I.data_has_rho;B.prepared.sel_rho=J(B.prepared.population.N,1,.)
 B.prepared.sel_rho[_cdid_ref_ix(B.prepared.U)]=J(sum(B.prepared.U),1,1)
 if(B.prepared.sel_has_rho)B.prepared.sel_rho=cdid_ci_transport(B.prepared.population,I,B.prepared.U)
 B.prepared.X=B.prepared.Z_sel=J(0,0,.);B.prepared.x_names=B.prepared.sel_names=J(1,0,"")
 B.prepared.columns=B.prepared.sel_columns=J(0,1,cdid_population_column())
 B.prepared.weight.fix_weights=fixweights
 B.prepared.clustered=cdid_ci_clusters(B.prepared.population,I,B.prepared.cluster_kind,
  B.prepared.cluster_numeric_id,B.prepared.cluster_string_id,B.prepared.unit_cluster,B.prepared.C)
 B.request=cdid_public_request_build(B.prepared.population,targets,control,ant,base,pt,start,cap)
 B.settings=cdid_public_defaults();B.settings.method=method;B.settings.control_group=control
 B.settings.anticipation=ant;B.settings.base_period=base;B.settings.pt_assumption=pt;B.settings.pt_start=start
 B.settings.selection=selection;B.settings.selection_canonical=decl.canonical_assumption
 B.settings.selection_model=decl.selection_model;B.settings.selection_pooling=decl.selection_pooling
 B.settings.selection_universe=decl.selection_universe;B.settings.selection_est_method="ipw"
 B.settings.selection_justification=decl.justification;B.settings.inactive_arguments=B.prepared.inactive_arguments
 B.settings.gmm_covariance=covariance;B.settings.imputation_pooling="";B.settings.outcome_pooling=outpool
 p=crexternal(output);if(p==NULL)_cdid_ref_fail("Selected output name already exists.",498)
 *p=B
 st_local("pr_N",strofreal(B.prepared.population.N));st_local("pr_A",strofreal(B.prepared.population.n_A))
 st_local("pr_post",strofreal(rows(B.request.post)));st_local("pr_aux",strofreal(rows(B.request.aux)))
 st_local("pr_report",strofreal(rows(B.request.report)))
}
void _cdid_sr_map(string scalar name,string scalar xf,string scalar xframe,
 string scalar sf,string scalar sframe,string scalar w,string scalar wframe)
{
 string rowvector x,f
 string matrix map
 real scalar role,j,k,found
 pointer() scalar p
 map=J(0,2,"")
 for(role=1;role<=3;role++) {
  x=(role==1 ? cdid_ref_formula_names(xf) : (role==2 ? cdid_ref_formula_names(sf) : tokens(w)))
  f=(role==1 ? cdid_ref_formula_names(xframe) : (role==2 ? cdid_ref_formula_names(sframe) : tokens(wframe)))
  if(cols(f)) {
   if(cols(x)!=cols(f))_cdid_ref_fail("Data/frame raw formula mappings disagree in dimension.",459)
   for(j=1;j<=cols(x);j++) {
    found=0
    for(k=1;k<=rows(map);k++)if(map[k,1]==x[j]) {
     if(map[k,2]!=f[j])_cdid_ref_fail("Shared formula variable has conflicting frame mappings.",459)
     found=1
    }
    if(!found)map=map\(x[j],f[j])
   }
  }
 }
 p=crexternal(name);if(p==NULL)_cdid_ref_fail("Transient frame map is in use.",498)
 *p=map
}
void _cdid_sr_reference_begin(string scalar bundle,string scalar popname,string scalar planname,
 string scalar rowname,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_catalog_opts scalar co
 struct cdid_pairs scalar pairs
 struct cdid_cat scalar C
 struct cdid_baseline_plan_state scalar plan
 struct cdid_sr_plan_cache scalar cache
 pointer() scalar p
 B=*findexternal(bundle)
 co=cdid_catalog_default_opts();co.mem_limit_bytes=8*cap;co.conventions="D10_D13_admissible_frame"
 pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
 C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
  B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
 plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,cap,co.work_cap)
 p=crexternal(popname);if(p==NULL)_cdid_ref_fail("Transient population name is in use.",498)
 *p=B.prepared.population
 p=crexternal(planname);if(p==NULL)_cdid_ref_fail("Transient plan name is in use.",498)
 *p=plan
 p=crexternal(rowname);if(p==NULL)_cdid_ref_fail("Transient source-row name is in use.",498)
 *p=B.prepared.source_row
 cache.valid=1
 cache.N=B.prepared.population.N;cache.R=B.prepared.population.R;cache.n_A=B.prepared.population.n_A
 cache.entries_cap=cap;cache.anticipation=B.request.anticipation;cache.control_group=B.request.control_group
 cache.calendar=B.request.calendar;cache.cohorts=B.request.cohorts;cache.target_cohorts=B.request.target_cohorts
 cache.first_row=B.prepared.population.first_row;cache.pop_cohort=B.prepared.population.cohort
 cache.post=B.request.post;cache.aux=B.request.aux
 cache.pairs=pairs;cache.C=C;cache.plan=plan
 p=findexternal("_cdid_sr_cached_plan")
 if(p==NULL) p=crexternal("_cdid_sr_cached_plan")
 if(p!=NULL) *p=cache
}
/* Reuse the baseline real-row staging and ONE FV encoder. Only the fitting
   roster/intercept differ: all selected X on A, selection design on U. */
void _cdid_sr_probe_start(string scalar bundle,string scalar raw,string scalar probe,
 string scalar encoded,string scalar formula,real scalar selection,real scalar cap)
{
 struct cdid_public_bundle scalar B,T
 B=*findexternal(bundle);T=B
 T.prepared.required_X=(selection ? B.prepared.U : B.prepared.population.A)
 T.prepared.x_has_intercept=(selection ? B.prepared.sel_has_intercept : B.prepared.x_has_intercept)
 *findexternal(bundle)=T
 _cdid_bs_start(bundle,raw,probe,encoded,formula,cap)
 *findexternal(bundle)=B
}
void _cdid_sr_probe_share(string scalar bundle,string scalar sprobe,string scalar xprobe,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar s,x
 struct cdid_population scalar P
 real colvector u_A,in_U
 real matrix X
 pointer() scalar p
 real scalar n_A,n_hist_A,intercept,shared
 B=*findexternal(bundle);s=*findexternal(sprobe);P=*s.raw.population
 u_A=_cdid_ref_ix(B.prepared.population.A)
 intercept=B.prepared.x_has_intercept
 shared=0
 if(s.accepted & rows(u_A)>0 & (cols(s.unit_raw)+intercept)>0 &
    _cdid_ref_is_plain(s.formula,s.names)) {
  x.raw=s.raw;x.raw.max_entries=cap
  x.names=s.names;x.variables=s.variables
  x.units=u_A;x.history=J(0,1,0)
  x.formula=s.formula;x.ready=1;x.accepted=1;x.reason=""
  if(rows(u_A)==rows(s.units)) x.unit_raw=s.unit_raw
  else {
   in_U=_cdid_ref_ix(B.prepared.population.A[s.units])
   x.unit_raw=s.unit_raw[in_U,]
  }
  x.history_raw=J(0,cols(x.names),.)
  n_A=rows(u_A)
  n_hist_A=(n_A==P.N ? P.R : sum(P.row_count[u_A]))
  x.retained=_cdid_ref_cost(P.R,cols(x.raw.raw_names)+3,cap,0)
  x.retained=_cdid_ref_cost(P.N,2*cols(x.raw.raw_names)+3,cap,x.retained)
  x.retained=_cdid_ref_cost(n_A+n_hist_A,cols(x.names)+3,cap,x.retained)
  x.retained=_cdid_ref_cost(n_A+n_hist_A,cols(x.names)+intercept,cap,x.retained)
  X=(intercept ? (J(n_A,1,1),x.unit_raw) : x.unit_raw)
  x.design.keys=u_A;x.design.value=X
  x.design.raw_index=J(1,0,0);x.design.raw_names=J(1,0,"")
  x.design.formula=x.formula;x.design.intercept=intercept
  x.design.encoded=1;x.design.dated=0
  x.design.columns=(intercept ? ("_cons",x.names) : x.names)
  x.design.expanded=x.names
  x.design.omitted=J(1,cols(x.names),0)
  x.design.column_map=(cols(x.names) ? (intercept+1..intercept+cols(x.names)) : J(1,0,0))
  p=crexternal(xprobe);if(p==NULL)_cdid_ref_fail("Static probe object already exists.",498)
  *p=x
  shared=1
 }
 st_local("sr_shared",strofreal(shared))
}
void _cdid_sr_probe_finish(string scalar name,string scalar encoded,real scalar selection,real scalar rc)
{
 struct cdid_bs_probe scalar q
 struct cdid_reference_inputs scalar E
 struct cdid_population scalar P
 real matrix X,H
 real colvector hit,u_pos
 real scalar j
 q=*findexternal(name);P=*q.raw.population
 if(rc!=0 & rc!=198 & rc!=459)exit(rc)
 if(q.ready & rc==0) {
  E=*findexternal(encoded)
  if(rows(E.designs[1].value)==0 & E.designs[1].encoded) {
   if(!selection & cols(E.designs[1].columns)==0)q.reason="Covariates must be supplied on the fitting roster."
   else {
    q.accepted=1
    if(rows(q.history) & cols(q.unit_raw)>0) {
     if(rows(q.units)==P.N) {
      hit=(rows(q.history)==P.R ? P.ids.row_unit : P.ids.row_unit[q.history])
     }
     else {u_pos=J(P.N,1,0);u_pos[q.units]=(1::rows(q.units));hit=u_pos[P.ids.row_unit[q.history]];}
     for(j=1;j<=cols(q.unit_raw);j++) {
      if(!all(q.history_raw[,j]:==q.unit_raw[hit,j])) {
       q.accepted=0;q.reason="Evaluated design varies within unit or disagrees with the unit frame."
       break
      }
     }
    }
    if(q.accepted) {
     X=(E.designs[1].intercept ? (J(rows(q.units),1,1),q.unit_raw) : q.unit_raw)
     q.design=E.designs[1];q.design.keys=q.units;q.design.value=X;q.design.raw_index=J(1,0,0);q.design.raw_names=J(1,0,"");q.retained=E.retained_entries
    }
   }
  }
  else {
   X=E.designs[1].value[(1::rows(q.units)),]
   if(!selection & cols(X)==0)q.reason="Covariates must be supplied on the fitting roster."
   else {
    q.accepted=1
    if(rows(q.history)) {
     H=E.designs[1].value[rows(q.units)+1..rows(E.designs[1].value),]
     if(rows(q.units)==P.N)hit=P.ids.row_unit[q.history]
     else {u_pos=J(P.N,1,0);u_pos[q.units]=(1::rows(q.units));hit=u_pos[P.ids.row_unit[q.history]];}
     if(!all(H:==X[hit,])){q.accepted=0;q.reason="Evaluated design varies within unit or disagrees with the unit frame.";}
    }
    if(q.accepted){q.design=E.designs[1];q.design.keys=q.units;q.design.value=X;q.design.raw_index=J(1,0,0);q.design.raw_names=J(1,0,"");q.retained=E.retained_entries;}
   }
  }
 }
 else if(rc!=0)q.reason="Evaluated static FV design is missing, nonfinite, or invalid (native code "+strofreal(rc)+")."
 *findexternal(name)=q
}
/* Classify both formulas only after a failed evaluated gate, as source's
   covariate_gate does. A failed static role is re-raised unless its own formula
   qualifies as dated; a passing role retains its source-specific static cache. */
void _cdid_sr_reference(string scalar bundle,string scalar raw,string scalar xprobe,
 string scalar sprobe,string scalar planname,string scalar output,string scalar time,
 string scalar weight,string scalar fixweights,real scalar anticipation,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar x,s
 struct cdid_reference_inputs scalar S
 struct cdid_baseline_plan_state scalar plan
 string rowvector names
 string scalar xr,sr,why
 real rowvector ix,wi
 real scalar marker,xd,sd,role,slot,j
 pointer() scalar p
 B=*findexternal(bundle);x=*findexternal(xprobe);s=*findexternal(sprobe);plan=*findexternal(planname)
 marker=!x.accepted|!s.accepted;xd=sd=0
 /* Cross-pair fits read the raw values of static covariates on the whole selection universe. */
 names=((marker | cdid_xp_requested()) ? x.raw.raw_names : tokens(weight));ix=_cdid_ref_names(x.raw.raw_names,invtokens(names))
 S=cdid_ref_source(x.raw.population,x.raw.source_rows,names,x.raw.data[,ix],x.raw.attributes[,ix],
  x.raw.has_data[ix],x.raw.has_frame[ix],time,cap)
 if(marker) {
  ix=_cdid_ref_names(S.raw_names,invtokens(x.names));if(cols(ix))xd=any(S.kind[ix]:!=1)
  ix=_cdid_ref_names(S.raw_names,invtokens(s.names));if(cols(ix))sd=any(S.kind[ix]:!=1)
  if(!x.accepted & !xd)_cdid_pr_fail("cdid_chained_request_covariate_error",x.reason)
  if(!s.accepted & !sd)_cdid_pr_fail("cdid_chained_request_covariate_error",s.reason)
 }
 xr=(xd ? invtokens(x.names) : "");sr=(sd ? invtokens(s.names) : "")
 wi=(weight=="" ? J(1,0,0) : _cdid_ref_names(S.raw_names,weight))
 if(cols(wi))if(S.has_data[wi[1]])if(any(select(S.data[,wi[1]],S.data[,wi[1]]:<.):<=0))
  _cdid_ref_fail("Sampling weights must be finite and positive where supplied.",459)
 cdid_ref_require(S,plan,anticipation,"selected",B.prepared.method,B.prepared.U,fixweights,
  x.formula,(B.prepared.method=="reg" ? "" : x.formula),s.formula,
  xr,(B.prepared.method=="reg" ? "" : xr),sr,B.prepared.x_has_intercept,
  (B.prepared.method=="reg" ? 0 : B.prepared.x_has_intercept),B.prepared.sel_has_intercept,
  (cols(wi) ? wi[1] : 0))
 cdid_ref_designs(S)
 for(role=1;role<=3;role++)if(S.roles[role].active & !S.roles[role].dated) {
  slot=S.roles[role].slot
  if(role==3)S.designs[slot]=s.design
  else S.designs[slot]=x.design
  S.retained_entries=_cdid_ref_cost(rows(S.designs[slot].value),cols(S.designs[slot].value),cap,S.retained_entries)
 }
 S.retained_entries=_cdid_ref_cost(1,x.retained+s.retained,cap,S.retained_entries)
 B.selected_static_status=(marker ? "failed_dated_fallback" : "passed")
 why=(!x.accepted ? "outcome: "+x.reason : "")+(!s.accepted ? " selection: "+s.reason : "")
 B.selected_static_reason=why
 /* This helper owns the exact static / first_period_static /
    reference_recipes reduction and canonical population owner descriptor. */
 if(weight=="")B.prepared.weight=cdid_population_weights(B.prepared.population,J(0,1,.),J(0,1,.),0,0,fixweights,marker)
 else {
  j=wi[1]
  B.prepared.weight=cdid_population_weights(B.prepared.population,
   (S.has_data[j] ? S.data[,j] : J(0,1,.)),(S.has_frame[j] ? S.attributes[,j] : J(0,1,.)),
   S.has_data[j],S.has_frame[j],fixweights,marker)
 }
 if(B.prepared.weight.reduction!="static" & B.prepared.weight.reduction!="first_period_static" & B.prepared.weight.reduction!="reference_recipes")
  _cdid_ref_fail("Unknown accepted weight reduction.",498)
 B.prepared.weight.V=S.V
 *findexternal(bundle)=B
 p=crexternal(output);if(p==NULL)_cdid_ref_fail("Reference output name is in use.",498)
 *p=S
}
void _cdid_sr_slot_meta(string scalar name,real scalar slot)
{
 struct cdid_reference_inputs scalar S
 S=*findexternal(name)
 if(!S.designs[slot].encoded & _cdid_ref_is_plain(S.designs[slot].formula,S.designs[slot].raw_names)) {
  _cdid_ref_fast_encode(S,slot)
  *findexternal(name)=S
 }
 st_local("sr_encoded",strofreal(S.designs[slot].encoded))
 _cdid_ref_meta(name,slot)
}
void _cdid_sr_attach(string scalar bundle,string scalar reference)
{
 struct cdid_public_bundle scalar B
 struct cdid_reference_inputs scalar S
 real scalar active,r,slot
 real matrix V
 B=*findexternal(bundle);S=*findexternal(reference)
 active=B.prepared.weight.reduction=="reference_recipes"
 B.prepared.design=(active ? "reference_recipes" : "static_numeric")
 for(r=1;r<=3;r=r+2) {
  slot=S.roles[r].slot
  V=J((active ? 0 : S.N),cols(S.designs[slot].value),.)
  if(!active)V[S.designs[slot].keys,]=S.designs[slot].value
  if(r==1){B.prepared.X=V;B.prepared.x_names=S.designs[slot].columns;}
  else{B.prepared.Z_sel=V;B.prepared.sel_names=S.designs[slot].columns;}
 }
 S.population=NULL;B.reference=S;*findexternal(bundle)=B
}

void cdid_selected_reference_metadata(struct cdid_result_state scalar out,
 struct cdid_selref_provider scalar f,struct cdid_baseline_plan_state scalar plan,
 struct cdid_reference_inputs scalar S,struct cdid_reference_catalog scalar catalog)
{
 real scalar i,c,r,d,a
 real colvector K
 string scalar key,id
 out.numeric_settings=out.numeric_settings\f.n_reference_recipes\f.n_recipes\f.n_consumers
 out.numeric_setting_names=out.numeric_setting_names\"n_reference_recipes"\"n_selected_recipes"\"n_unique_consumers"
 out.text_settings=out.text_settings\("reference_route","selected")\("reference_weight_kind",S.weight_kind)
 for(i=1;i<=rows(plan.raw);i++) {
  key=_cdid_rb_key(i);c=catalog.raw_consumer[i];r=S.view_recipe[i]
  id=(i<=plan.M ? "cp04:catalog:"+_cdid_rb_key(i) : "cp04:auxiliary:"+_cdid_rb_key(i-plan.M))
  out.numeric_settings=out.numeric_settings\(c>0 ? f.points[c] : .)\plan.cohorts[plan.raw[i,1]]\
   plan.calendar[plan.raw[i,2]]\plan.calendar[plan.raw[i,3]]\r\plan.raw_supported[i]\c
  out.numeric_setting_names=out.numeric_setting_names\("source_link_point:"+key)\("source_link_g:"+key)\
   ("source_link_s:"+key)\("source_link_t:"+key)\("source_link_recipe:"+key)\("source_link_supported:"+key)\("source_link_consumer:"+key)
  K=plan.raw_K_index[plan.raw_K_ptr[i]..plan.raw_K_ptr[i+1]-1]
  d=S.recipes[r,1];a=S.recipes[r,2]
  out.numeric_settings=out.numeric_settings\(d ? S.calendar[d] : .)\(a ? S.calendar[a] : .)\sum(2:^(K:-1))
  out.numeric_setting_names=out.numeric_setting_names\("source_link_design:"+key)\("source_link_anchor:"+key)\("source_link_controls:"+key)
  out.text_settings=out.text_settings\("source_link_id:"+key,id)
 }
 out.source_versions=out.source_versions\("estimator","observation-weighted reference-date estimator; influence functions over all units")
}

struct cdid_result_state scalar cdid_selected_reference_report(
    struct cdid_public_bundle scalar B,struct cdid_selref_provider scalar fit,
    real colvector post_point,real colvector post_available,string colvector post_status,
    string colvector post_reason,real colvector post_est_ids,real matrix post_loadings,
    real scalar compute_if,struct cdid_state_limits scalar limits,real scalar chunk,real scalar cap,
    | real matrix retained_if)
{
    struct cdid_public_report_buffer scalar buffer
    struct cdid_result_state scalar out
    real colvector aux_point,aux_available,active,qids,postpos,auxpos,consumers,hit,uid,aux_in_post,rem_consumers,rem_pos
    real matrix values,postIF,auxIF,ret_if
    real scalar Q,M,Jt,i,j,q,lo,hi,effective,perrow,retained,use_ret,scale_f
    Q=rows(B.request.aux);M=rows(fit.post_consumer);Jt=rows(post_est_ids)
    ret_if=(args()>=13 ? retained_if : J(0,0,.))
    if(rows(post_loadings)!=Jt | cols(post_loadings)!=M | rows(fit.aux_consumer)!=Q |
        chunk<1 | chunk>=. | chunk!=floor(chunk) | cap<1 | cap>=. | cap!=floor(cap))
        _cdid_pr_fail("cdid_native_reporting_input","Provider/post loading identities or replay controls do not align.")
    if(rows(post_est_ids)) {
        _cdid_rb_ids(post_est_ids,rows(B.request.post))
        if(any(post_available[post_est_ids]:!=1) | rows(post_est_ids)!=sum(post_available))
            _cdid_pr_fail("cdid_native_reporting_input","Estimable post IDs must equal the full supported post set.")
    }
    else if(sum(post_available)) _cdid_pr_fail("cdid_native_reporting_input","Supported post set has no passed target loading.")
    aux_point=J(Q,1,.);aux_available=fit.aux_supported
    active=_cdid_rb_index(aux_available)
    if(rows(active)) aux_point[active]=fit.points[fit.aux_consumer[active]]
    buffer=cdid_public_report_start(B,post_point,post_available,post_status,post_reason,
        aux_point,aux_available,fit.aux_reason,compute_if,limits)
    if(compute_if) {
        qids=J(0,1,.)
        for(i=1;i<=rows(B.request.report);i++) if(B.request.report[i,5]==2 & buffer.payload.available[i])
            qids=qids\B.request.report[i,7]
        if(rows(qids)) qids=uniqrows(sort(qids,1))
        consumers=(Jt ? fit.post_consumer : J(0,1,.))
        if(rows(qids)) consumers=consumers\fit.aux_consumer[qids]
        if(rows(consumers)) consumers=uniqrows(sort(consumers,1))
        postpos=J((Jt ? M : 0),1,.)
        for(j=1;j<=rows(postpos);j++) {
            hit=_cdid_rb_index(consumers:==fit.post_consumer[j]);postpos[j]=hit[1]
        }
        auxpos=J(rows(qids),1,.)
        for(j=1;j<=rows(qids);j++) {
            hit=_cdid_rb_index(consumers:==fit.aux_consumer[qids[j]]);auxpos[j]=hit[1]
        }
        // Retained result + provider + target loading, plus conservative bounded
        // provider/scatter scratch. Never construct originalF-by-allmoments IF.
        retained=cdid_state_cost(buffer.payload)[1]+fit.workspace_entries+length(post_loadings)
        perrow=10*(rows(consumers)+2*cols(B.prepared.X)+cols(B.prepared.Z_sel)+4)+3*(Jt+rows(qids)+rows(B.request.report))
        if(retained>=cap | perrow>cap-retained) _cdid_rs_fail("Reported IF/provider replay exceeds explicit numeric workspace.",908)
        effective=min((chunk,fit.chunk_rows,floor((cap-retained)/perrow)))
        if(fit.cache_eval_active & rows(fit.cached_eval_phi)==fit.n_F &
            cols(fit.cached_eval_phi)==rows(consumers) & fit.eval_ids_sorted==consumers) {
            for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                values=fit.cached_eval_phi[uid,]
                postIF=(Jt ? values[,postpos']*post_loadings' : J(hi-lo+1,0,.))
                auxIF=(rows(qids) ? values[,auxpos'] : J(hi-lo+1,0,.))
                cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
            }
        }
        else {
            use_ret=(rows(ret_if)==fit.n_F & cols(ret_if)==M & fit.model.lambda>0 & fit.model.lambda<.)
            if(use_ret) {
                scale_f=1/fit.model.lambda
                aux_in_post=J(rows(qids),1,0);rem_consumers=J(0,1,.);rem_pos=J(rows(qids),1,0)
                for(j=1;j<=rows(qids);j++) {
                    hit=_cdid_rb_index(fit.post_consumer:==fit.aux_consumer[qids[j]])
                    if(rows(hit)) aux_in_post[j]=hit[1]
                    else rem_consumers=rem_consumers\fit.aux_consumer[qids[j]]
                }
                if(rows(rem_consumers)) {
                    rem_consumers=uniqrows(sort(rem_consumers,1))
                    for(j=1;j<=rows(qids);j++) if(!aux_in_post[j]) {
                        hit=_cdid_rb_index(rem_consumers:==fit.aux_consumer[qids[j]])
                        rem_pos[j]=hit[1]
                    }
                }
                for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                    hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                    postIF=(Jt ? (scale_f:*ret_if[uid,])*post_loadings' : J(hi-lo+1,0,.))
                    if(rows(qids)) {
                        auxIF=J(hi-lo+1,rows(qids),.)
                        values=(rows(rem_consumers) ? cdid_selref_if_chunk(fit,uid,rem_consumers,0) : J(hi-lo+1,0,.))
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
                for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                    hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                    values=(rows(consumers) ? cdid_selref_if_chunk(fit,uid,consumers,0) : J(hi-lo+1,0,.))
                    postIF=(Jt ? values[,postpos']*post_loadings' : J(hi-lo+1,0,.))
                    auxIF=(rows(qids) ? values[,auxpos'] : J(hi-lo+1,0,.))
                    cdid_public_report_add(buffer,uid,post_est_ids,postIF,qids,auxIF)
                }
            }
        }
    }
    out=cdid_public_report_finish(buffer,limits)
    out.numeric_settings=out.numeric_settings\fit.weight_rows\fit.normalization_count\fit.model.n_treatment_fit_calls\fit.model.n_or_block_fits\M\Q
    out.numeric_setting_names=out.numeric_setting_names\"GMM_replay_rows"\"GMM_normalization_count"\"n_treatment_fits"\"n_pooled_OR_fits"\"M_full"\"Q_aux_union"
    out.numeric_settings=out.numeric_settings\fit.model.n_selection_blocks\fit.model.n_selection_blocks\sum(B.prepared.U)\cols(B.prepared.Z_sel)\B.prepared.sel_has_intercept\B.prepared.selection_formula_inherited\B.prepared.sel_has_rho
    out.numeric_setting_names=out.numeric_setting_names\"n_selection_fits"\"n_selection_blocks"\"N_U"\"p_selection"\"sel_has_intercept"\"selection_formula_inherited"\"sel_has_rho"
    out.text_settings=out.text_settings\("GMM_weight_row_domain",fit.weight_row_domain)\("provider_status",fit.status)\("profile_kind",fit.profile_kind)
    out.source_versions=out.source_versions\("estimator","post-treatment and auxiliary comparisons; influence functions over all units")
    cdid_state_validate(out,limits)
    return(out)
}

struct cdid_public_fit_result scalar cdid_selref_public_fit(struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_selref_provider scalar fit
    struct cdid_selref_controls scalar selected_controls
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_population scalar owner
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar catalog
    struct cdid_gmm_result scalar sol
    real matrix Z,G,loadings,postmeta,Vw,Lx
    real colvector point,diagnostic,reportpost,qids,eval_ids,rk,kept
    string colvector diagnostic_names
    real scalar M,Jt,used,remaining,i,j,p,keepif,N,def_chunk,eff_chunk
    if(B.prepared.design!="reference_recipes" | B.prepared.selection_canonical!="conditional_mtar")
        _cdid_pr_fail("cdid_native_private_binding_scope","Selected binding requires the active conditional-MTAR reference route.")
    if(B.settings.gmm_weighting!="optimal" & B.settings.gmm_weighting!="two_step" & B.settings.gmm_weighting!="identity" & B.settings.gmm_weighting!="within")
        _cdid_pr_fail("cdid_chained_public_control_error","gmmweighting() must be optimal, two_step, identity, or within.")
    B.settings.alpha=ctrl.inference.alpha;B.settings.compute_if=ctrl.inference.compute_inffunc
    B.settings.bstrap=ctrl.inference.bstrap;B.settings.cband=ctrl.inference.cband;B.settings.biters=ctrl.inference.B
    if(!_cdid_sel_try_cached_plan(B,ctrl,pairs,C,plan)) {
        co=cdid_catalog_default_opts();co.mem_limit_bytes=ctrl.workspace_bytes;co.work_cap=ctrl.workcap
        co.conventions="D10_D13_admissible_frame"
        pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
        C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
            B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
        plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,ctrl.entries,ctrl.workcap)
    }
    D=cdid_chained_ipw_build_design(C,B.request.post,B.request.pt_assumption,B.request.pt_start)
    B.request=cdid_public_request_tags(B.request,C)
    owner=B.prepared.population;reference=B.reference;reference.population=&owner
    cdid_ref_activate(reference)
    catalog=cdid_reference_catalog_make(plan,reference,ctrl.entries)
    N=owner.N
    p=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    def_chunk=max((1,min((N,16384,floor(ctrl.entries/max((1,32*(ctrl.links_bound+p+1))))))))
    eff_chunk=ctrl.chunk
    if(ctrl.workcap>=. & ctrl.chunk==def_chunk & ctrl.entries>=32*(plan.M+p+1)*min((N,16384))) {
        eff_chunk=max((ctrl.chunk,min((N,16384,floor(ctrl.entries/max((1,32*(plan.M+p+1))))))))
    }
    selected_controls=cdid_selref_default_controls(owner.N)
    selected_controls.method=B.prepared.method;selected_controls.selection_model=B.prepared.selection_model
    selected_controls.outcome_pooling=(B.prepared.outcome_pooling=="" ? "separate" : B.prepared.outcome_pooling)
    selected_controls.selection_pooling=B.prepared.selection_pooling;selected_controls.universe=B.prepared.selection_universe
    selected_controls.justification=B.prepared.selection_justification;selected_controls.has_rho=B.prepared.sel_has_rho
    selected_controls.rho=B.prepared.sel_rho;selected_controls.eps=ctrl.baseline.glm_epsilon
    selected_controls.maxit=ctrl.baseline.glm_maxit;selected_controls.rank_tol=ctrl.baseline.rank_tol
    selected_controls.chunk_rows=eff_chunk;selected_controls.max_work=(ctrl.workcap<. ? ctrl.workcap : 0)
    selected_controls.max_entries=ctrl.entries
    selected_controls.points_only=(!ctrl.inference.compute_inffunc & B.settings.gmm_weighting=="identity")
    fit=cdid_selref_prepare(&owner,&reference,&plan,&catalog,selected_controls)
    fit.workspace_entries=_cdid_base_cells(1,reference.retained_entries,ctrl.entries,fit.workspace_entries)
    sys=cdid_baseline_make_system(C,D,B.request.post,ctrl.entries,fit.workspace_entries)
    support=cdid_public_post_ledger(C,D,B.request,sys)
    M=sys.M;Jt=rows(B.request.post);point=J(Jt,1,.);loadings=J(0,M,.)
    diagnostic=J(0,1,.);diagnostic_names=J(0,1,"")
    if(M>0 & sys.p_retained>0) {
        used=fit.workspace_entries+sys.workspace_entries
        used=_cdid_base_cells(8*ctrl.chunk,M+2*cols(B.prepared.X)+4,ctrl.entries,used)
        if(B.settings.gmm_weighting=="identity") {
            used=_cdid_gmm_identity_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries
            sol=cdid_gmm_identity_map(fit.points[catalog.post_consumer],sys.full_A,sys.H_est,remaining)
        }
        else {
            used=_cdid_base_cells(2*ctrl.chunk,M,ctrl.entries,used)
            used=_cdid_base_cells(3*M,M,ctrl.entries,used)
            used=_cdid_gmm_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries-M*M
            if(ctrl.inference.compute_inffunc & ctrl.workcap>=. & sys.n_estimable>0 & eff_chunk>=fit.n_F) {
                qids=J(0,1,.)
                for(i=1;i<=rows(B.request.report);i++) {
                    if(B.request.report[i,5]==2 & B.request.report[i,9]==0) {
                        if(fit.aux_supported[B.request.report[i,7]])
                            qids=qids\B.request.report[i,7]
                    }
                }
                if(rows(qids)) {
                    qids=uniqrows(sort(qids,1))
                    eval_ids=uniqrows(sort(catalog.post_consumer\fit.aux_consumer[qids],1))
                    if(rows(eval_ids)>rows(fit.post_ids_sorted) & used+fit.n_F*rows(eval_ids)<remaining) {
                        fit.eval_ids_sorted=eval_ids
                        fit.eval_post_map=_cdid_sel_link_match_sorted(catalog.post_consumer,eval_ids)
                        fit.cached_eval_phi=J(0,0,.)
                        fit.cache_eval_active=1
                    }
                }
            }
            if(B.settings.gmm_weighting=="within") {
                rk=fit.points[catalog.post_consumer]
                cdid_wsg_within(B,plan,C,rk,1,ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,eff_chunk,kept,G,Vw,reference,&fit)
                sol=cdid_gmm_fit_within_norm(rk[kept],sys.full_A[kept,],Vw,G,M,sys.H_est,ctrl.condition_trigger,ctrl.ridge_relative,
                    fit.weight_rows,fit.normalization_count,eff_chunk,&cdid_selref_post_core(),&fit,remaining,M,kept')
                Lx=J(rows(sol.target_loadings),M,0);Lx[.,kept']=sol.target_loadings;sol.target_loadings=Lx
                Lx=Vw=J(0,0,.)
            }
            else {
                kept=(1::M)
                G=I(M)
                sol=cdid_gmm_fit_unit_norm(fit.points[catalog.post_consumer],sys.full_A,G,M,sys.H_est,
                    B.settings.gmm_weighting,ctrl.condition_trigger,ctrl.ridge_relative,
                    fit.weight_rows,fit.normalization_count,eff_chunk,&cdid_selref_post_core(),&fit,remaining)
            }
        }
        if(B.settings.gmm_weighting=="identity") kept=(1::M)
        if(sys.n_estimable) point[sys.estimable_ids]=sol.targets
        loadings=sol.target_loadings
        cdid_chain_sp_record(B,C,D,kept,fit.points[catalog.post_consumer],loadings,sys.estimable_ids,
            cdid_chain_auxpts(fit.points,fit.aux_consumer,fit.aux_supported),fit.aux_supported)
        diagnostic=(sol.criterion\sol.lambda\sol.rho\sol.trace_B\sol.covariance_rank\sol.source_n_units\sol.normalization_factor\sol.sample_rank_upper)
        diagnostic_names=("GMM_criterion"\"GMM_lambda"\"GMM_ridge_rho"\"GMM_trace_B"\"GMM_covariance_rank"\"GMM_source_replay_rows"\"GMM_gram_normalization_factor"\"GMM_sample_rank_upper")
        // Only final reported IF is retained for durable no-refit results.
        G=J(0,0,.);sol.unit_gram_covariance=sol.coefficient_loadings=sol.target_loadings=J(0,0,.)
        sol.coefficient_covariance=sol.covariance=J(0,0,.)
    }
    else if(sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input","An estimable post target has no GMM parameter system.")
    sys.full_A=sys.H_est=J(0,0,.)
    out.payload=cdid_selected_reference_report(B,fit,point,sys.post_available,support.status,support.reason,
        sys.estimable_ids,loadings,ctrl.inference.compute_inffunc,ctrl.inference.limits,eff_chunk,ctrl.entries,sol.retained_if)
    fit.cached_eval_phi=J(0,0,.)
    sol.retained_if=J(0,0,.)
    postmeta=J(rows(B.request.report),cols(support.ledger),.)
    reportpost=_cdid_rb_index(B.request.report[,5]:==3)
    for(i=1;i<=rows(reportpost);i++){j=reportpost[i];postmeta[j,]=support.ledger[B.request.report[j,6],];}
    out.payload.ledger=out.payload.ledger,postmeta
    out.payload.ledger_columns=out.payload.ledger_columns,"post_support.":+support.columns
    out.payload.numeric_settings=out.payload.numeric_settings\diagnostic\ctrl.workspace_bytes\ctrl.entries\ctrl.chunk\ctrl.workcap\ctrl.condition_trigger\ctrl.ridge_relative\ctrl.links_bound\ctrl.J_bound
    out.payload.numeric_setting_names=out.payload.numeric_setting_names\diagnostic_names\"numeric_workspace_bytes"\"numeric_workspace_entries"\"chunk_rows"\"counted_work_limit"\"condition_trigger"\"ridge_relative"\"catalogue_links_bound"\"inference_targets_bound"
    out.payload.text_settings=out.payload.text_settings\("GMM_metric","full original catalogue identity metric before covariance weighting")\("registry_limit_meaning","2^53-1 exact-count ceiling; not a session RAM/work budget")
    out.payload.source_versions=out.payload.source_versions\("gmm","GMM combination of all post-treatment comparisons, normalized by the units with two observed outcomes")
    cdid_selected_reference_metadata(out.payload,fit,plan,reference,catalog)
    out.inference=ctrl.inference
    cdid_state_validate(out.payload,ctrl.inference.limits)
    return(out)
}

end
