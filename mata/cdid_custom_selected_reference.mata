*! Private selected custom dated callback provider.
*! Fixed R/chained_custom_selection.R32-167/192-226.
*! Reuses shared reference/profile/scoped-reader primitives; no builtin treatment fit.
version 19.5
mata:
struct cdid_scr_model {
 real scalar representative,recipe,context,loading_col,target,s,t,m,estimate,callback_if_mean
 string scalar link_id
 real colvector cohorts,selection_pc,obs_A,psi,gradient
}
struct cdid_scr_recipe {
 struct cdid_sel_cell_map scalar map
 struct cdid_sel_fit_result scalar fit
 struct cdid_sel_up_map scalar upstream
 struct cdid_sel_moments scalar loading
 real colvector models
}
struct cdid_scr_provider {
 pointer(struct cdid_population scalar) scalar population
 pointer(struct cdid_reference_inputs scalar) scalar reference
 pointer(struct cdid_baseline_plan_state scalar) scalar plan
 pointer(struct cdid_reference_catalog scalar) scalar catalog
 /* Only profile leasing/reader fields are populated. No builtin scientific preparation. */
 struct cdid_selected_model scalar reader_owner
 struct cdid_scr_recipe colvector recipes
 struct cdid_scr_model colvector models
 real colvector points,raw_model,post_model,aux_model,recipe_reference,consumer_context
 real matrix keep
 real scalar N,n_A,n_models,n_recipes,n_callbacks,n_selection_fits,n_selection_registry_fits,n_selection_blocks,n_treatment_fits,n_endpoint_reads
 real scalar weight_rows,normalization_count,sample_rank_upper,core_to_frame,chunk_rows
 real scalar workspace_entries,workspace_peak_entries,workspace_limit,work_ops,max_selection_params
 string scalar weight_row_domain,status
}
real scalar _cdid_scr_pc(struct cdid_sel_cell_map scalar map,real scalar k,real scalar s,real scalar t)
{
 real colvector ix
 ix=vec(selectindex((map.pair_cohorts[,1]:==k):&(map.pair_cohorts[,2]:==s):&(map.pair_cohorts[,3]:==t)))
 if(rows(ix)!=1)_cdid_callback_fail("selection_identity","Required reference pair/cohort cell is not uniquely mapped.",498)
 return(ix[1])
}
void _cdid_scr_owner(struct cdid_scr_provider scalar f)
{
 if(f.population==NULL | f.reference==NULL | f.plan==NULL | f.catalog==NULL)
  _cdid_callback_fail("owner","Original scientific objects are required for selected reference replay.",498)
 _cdid_ref_check_population(*f.reference)
 if((*f.reference).population!=f.population | (*f.population).N!=f.N | (*f.population).n_A!=f.n_A |
  (*f.plan).n_F!=f.N | (*f.plan).n_A!=f.n_A | (*f.catalog).n_consumers!=f.n_models)
  _cdid_callback_fail("owner","Selected reference scientific object identity changed.",498)
}
real matrix _cdid_scr_observers(struct cdid_scr_provider scalar f,real scalar j,real scalar lo,real scalar hi)
{
 real matrix rr
 real colvector take
 real scalar i
 rr=_cdid_base_pair_rows(*f.population,f.models[j].s,f.models[j].t,lo,hi);take=J(rows(rr),1,0)
 for(i=1;i<=rows(rr);i++)take[i]=any(f.models[j].cohorts:==(*f.plan).a_cohort[rr[i,1]])
 /* A recipe restricted to a reference eligibility context reads its eligible units only. */
 if(rows(f.keep) & rows(rr))take=take:*(f.keep[(*f.plan).A_to_F[rr[,1]],f.models[j].context]:!=0)
 take=vec(selectindex(take));return(rows(take) ? rr[take,] : J(0,3,.))
}
void _cdid_scr_callback(struct cdid_scr_provider scalar f,real scalar j,pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_callback_request scalar q
 struct cdid_callback_result scalar val
 struct cdid_sel_pred scalar pred
 struct cdid_sel_moments scalar one
 real matrix rr
 real colvector gi,ix,co,sub
 real scalar r,m,h,slot,d,lo,hi,n,at,k,c
 r=f.models[j].context;_cdid_selected_activate_reader(f.reader_owner,r)
 q=cdid_callback_request_init();m=f.models[j].m;h=rows(f.models[j].cohorts)
 q.mode="chained";q.route="conditional_mtar_dated";q.link_id=f.models[j].link_id
 q.panel=q.inffunc=q.gradient_required=1;q.N=f.N;q.n_A=f.n_A
 q.call_index=q.model_index=j;q.raw_index=f.models[j].representative;q.recipe=f.models[j].recipe
 q.g=(*f.plan).cohorts[f.models[j].target];q.K=(*f.plan).cohorts[f.models[j].cohorts[2..h]]
 q.s=(*f.plan).calendar[f.models[j].s];q.t=(*f.plan).calendar[f.models[j].t]
 q.rows_A=q.rows_F=q.y0=q.y1=q.D=q.i_weights=J(m,1,.)
 slot=(*f.reference).roles[2].slot;q.covariate_names=(*f.reference).designs[slot].columns
 q.covariates=J(m,cols(q.covariate_names),.);gi=J(m,1,.)
 d=((*f.reference).roles[2].dated ? (*f.reference).recipes[q.recipe,1] : 0);at=0
 for(lo=1;lo<=f.n_A;lo=lo+f.chunk_rows) {
  hi=min((f.n_A,lo+f.chunk_rows-1));rr=_cdid_scr_observers(f,j,lo,hi);n=rows(rr)
  if(n==0)continue
  if(at+n>m)_cdid_callback_fail("identity","Actual selected reference support exceeds its admitted count.",498)
  ix=((at+1)::(at+n));q.rows_A[ix]=rr[,1];q.rows_F[ix]=(*f.plan).A_to_F[rr[,1]]
  q.y0[ix]=(*f.population).y[rr[,2]];q.y1[ix]=(*f.population).y[rr[,3]]
  q.covariates[ix,]=cdid_ref_design_rows(*f.reference,2,d*(f.N+1):+q.rows_F[ix])
  co=(*f.plan).a_cohort[rr[,1]]
  for(k=1;k<=h;k++) {
   sub=vec(selectindex(co:==f.models[j].cohorts[k]))
   if(rows(sub)) {
    gi[ix[sub]]=J(rows(sub),1,k)
    pred=cdid_sel_pred_pair(f.recipes[r].fit,f.reader_owner.reader,f.models[j].selection_pc[k],q.rows_F[ix[sub]])
    q.i_weights[ix[sub]]=pred.sel_weight
   }
  }
  at=at+n
 }
 if(at!=m | missing(q.y0) | missing(q.y1) | missing(q.y1-q.y0) | missing(q.covariates) |
  missing(q.i_weights) | any(q.i_weights:<=0) | missing(gi))
  _cdid_callback_fail("selected_input","Actual reference callback support/endpoints/design/weights are invalid.",504)
 q.D=(gi:==1)
 val=cdid_callback_invoke_chained(callback,q,payload);f.n_callbacks=f.n_callbacks+1
 f.points[j]=f.models[j].estimate=val.ATT[1,1];f.models[j].callback_if_mean=mean(val.psi)
 f.models[j].obs_A=q.rows_A;f.models[j].psi=val.psi;f.models[j].gradient=val.log_weight_gradient
 one=cdid_sel_moments_init(f.recipes[r].fit,1)
 for(k=1;k<=h;k++) {
  sub=vec(selectindex(gi:==k));if(rows(sub)==0)continue
  /* Total gradient enters LD'g literally, with no /m or /N. */
  cdid_sel_add_moments(one,f.recipes[r].fit,f.reader_owner.reader,f.recipes[r].upstream,
   f.models[j].selection_pc[k],q.rows_F[sub],val.log_weight_gradient[sub])
 }
 c=f.models[j].loading_col;f.recipes[r].loading.G_theta[,c]=one.G_theta;f.recipes[r].loading.G_ak[,c]=one.G_ak
 f.n_endpoint_reads=f.n_endpoint_reads+2*m
}
struct cdid_scr_provider scalar cdid_scr_prepare(
 pointer(struct cdid_population scalar) scalar population,pointer(struct cdid_reference_inputs scalar) scalar reference,
 pointer(struct cdid_baseline_plan_state scalar) scalar plan,pointer(struct cdid_reference_catalog scalar) scalar catalog,
 struct cdid_selref_controls scalar controls,pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_scr_provider scalar f
 struct cdid_selref_inputs scalar b
 struct cdid_sel_reader scalar temporary
 struct cdid_sel_decl scalar declaration
 struct cdid_sel_budget scalar budget
 real matrix specs
 real colvector count,cohorts,reader_cohort
 real scalar i,j,k,r,c,rep,m,pcrows,at,used,maxm,maxh,p,work,cap
 if(callback==NULL)_cdid_callback_fail("owner","A scientific callback is required.",498)
 if(controls.method!="ipw")_cdid_callback_fail("scope","Selected custom callbacks require selection IPW.",498)
 if(controls.selection_model=="independent_period_logit" &
  (controls.universe!="F" | controls.has_rho | strtrim(controls.justification)!=""))
  _cdid_callback_fail("scope","Date-product custom selection requires F with no rho/justification.",459)
 if(controls.chunk_rows<1 | controls.chunk_rows>=. | controls.chunk_rows!=floor(controls.chunk_rows) |
  controls.max_entries<1 | controls.max_entries>=. | controls.max_entries!=floor(controls.max_entries))
  _cdid_callback_fail("control","Explicit existing chunk/workspace controls are required.",498)
 b=cdid_selref_build(population,reference,plan,catalog,controls)
 f.population=population;f.reference=reference;f.plan=plan;f.catalog=catalog
 f.N=(*population).N;f.n_A=(*population).n_A;f.n_models=(*catalog).n_consumers;f.n_recipes=rows(b.recipes)
 f.recipe_reference=b.recipe_reference;f.consumer_context=b.consumer_recipe
 f.keep=b.recipe_keep
 f.raw_model=(*catalog).raw_consumer;f.post_model=(*catalog).post_consumer;f.aux_model=(*catalog).aux_consumer
 f.weight_rows=f.N;f.normalization_count=f.n_A;f.sample_rank_upper=f.N-1;f.core_to_frame=f.N/f.n_A
 f.weight_row_domain="F";f.chunk_rows=controls.chunk_rows;f.workspace_limit=controls.max_entries
 f.n_callbacks=f.n_selection_fits=f.n_selection_registry_fits=f.n_selection_blocks=f.n_treatment_fits=f.n_endpoint_reads=0;f.status="preflight"
 cap=(controls.max_work>0 ? controls.max_work : .);p=cols((*reference).designs[(*reference).roles[2].slot].value)
 if(p<1)_cdid_callback_fail("callback_design","The actual chained callback design must be nonempty.",459)
 /* Retain unique selection/anchor/rho profiles, not copied per recipe or model.
    Callback X stays in the single borrowed encoded reference union. */
 f.reader_owner.N=f.N;f.reader_owner.reader_ready=f.reader_owner.reader_active=f.reader_owner.reader_universe=0
 f.reader_owner.reader_sel_slot=f.reader_owner.reader_w_slot=f.reader_owner.reader_rho_slot=0
 f.reader_owner.profiles=b.registry;f.reader_owner.profiles.trt_data=cdid_selected_matrix_slot(0,1)
 f.reader_owner.contexts=cdid_selected_context(f.n_recipes,1)
 f.recipes=cdid_scr_recipe(f.n_recipes,1);f.models=cdid_scr_model(f.n_models,1);f.points=J(f.n_models,1,.)
 used=_cdid_base_cells(1,2*b.planned_entries+2*_cdid_cb_population_entries(*population)+
  _cdid_cb_catalog_entries(*catalog)+(*reference).retained_entries+(*plan).workspace_entries+64*f.n_models+64*f.n_recipes+24*f.N,
  f.workspace_limit,0)
 maxm=maxh=f.max_selection_params=0;work=(*plan).work_ops
 for(j=1;j<=f.n_models;j++) {
  rep=(*catalog).representative[j];r=f.consumer_context[j]
  f.models[j].representative=rep;f.models[j].recipe=(*catalog).consumer_recipe[j];f.models[j].context=r
  f.models[j].target=(*catalog).consumers[j,1];f.models[j].s=(*catalog).consumers[j,2];f.models[j].t=(*catalog).consumers[j,3]
  f.models[j].link_id=b.links[j].link_id
  f.models[j].cohorts=f.models[j].target\(*catalog).K_index[(*catalog).K_ptr[j]..((*catalog).K_ptr[j+1]-1)]
  count=(*plan).raw_group_count[(*plan).raw_count_ptr[rep]..((*plan).raw_count_ptr[rep+1]-1)]
  m=sum(count);if(m<2 | any(count:<1))_cdid_callback_fail("support","Positive actual-pair callback cohort support is required.",459)
  f.models[j].m=m;used=_cdid_base_cells(1,3*m+2*rows(count)+24,f.workspace_limit,used)
  maxm=max((maxm,m));maxh=max((maxh,rows(count)))
 }
 /* Admit every recipe and the ENTIRE retained scientific registry before calls. */
 for(r=1;r<=f.n_recipes;r++) {
  f.reader_owner.contexts[r].sel_id=b.recipes[r].sel_id
  f.reader_owner.contexts[r].weight_id=b.recipes[r].weight_id
  f.reader_owner.contexts[r].rho_id=b.recipes[r].rho_id
  declaration=cdid_sel_norm_decl("conditional_mtar",1,controls.universe,controls.has_rho,
   strtrim(controls.justification)!="",controls.selection_model,controls.selection_pooling,controls.justification)
  if(rows(f.keep)) temporary=cdid_sel_reader_scope(*population,_cdid_sld_sel_full(f.reader_owner,b.recipes[r].sel_id),
   b.recipes[r].sel_has_intercept,b.recipes[r].sel_int_col,_cdid_sld_w_full(f.reader_owner,b.recipes[r].weight_id),
   declaration,_cdid_sld_rho_full(f.reader_owner,b.recipes[r].rho_id),f.chunk_rows,b.scope[r,]',f.keep[,r])
  else temporary=cdid_sel_reader_scope(*population,_cdid_sld_sel_full(f.reader_owner,b.recipes[r].sel_id),
   b.recipes[r].sel_has_intercept,b.recipes[r].sel_int_col,_cdid_sld_w_full(f.reader_owner,b.recipes[r].weight_id),
   declaration,_cdid_sld_rho_full(f.reader_owner,b.recipes[r].rho_id),f.chunk_rows,b.scope[r,]')
  _cdid_selected_reader_store(f.reader_owner,r,temporary)
  if(rows(f.keep)) if(!all(f.keep[,r])) _cdid_selected_use_universe(f.reader_owner,r,temporary.U)
  reader_cohort=cdid_ref_lookup((*plan).cohorts,f.reader_owner.reader.cohorts)
  if(any(reader_cohort:<1))_cdid_callback_fail("identity","Scientific/selection cohort dictionaries differ.",498)
  f.recipes[r].models=vec(selectindex(f.consumer_context:==r));pcrows=0
  for(i=1;i<=rows(f.recipes[r].models);i++)pcrows=pcrows+rows(f.models[f.recipes[r].models[i]].cohorts)
  specs=J(pcrows,3,.);at=0
  for(i=1;i<=rows(f.recipes[r].models);i++) {
   j=f.recipes[r].models[i];f.models[j].loading_col=i;cohorts=f.models[j].cohorts
   for(k=1;k<=rows(cohorts);k++){at=at+1;specs[at,]=(reader_cohort[cohorts[k]],f.models[j].s,f.models[j].t);}
  }
  f.recipes[r].map=cdid_sel_partition(f.reader_owner.reader,specs)
  budget=cdid_sel_budget_check(f.reader_owner.reader,f.recipes[r].map,cap)
  f.max_selection_params=max((f.max_selection_params,budget.total_params))
  used=_cdid_base_cells(1,8*pcrows+4*budget.compact_state_doubles+budget.stream_workspace_doubles+
   2*(budget.total_params+rows(reader_cohort))*rows(f.recipes[r].models),f.workspace_limit,used)
  work=_cdid_base_work(budget.work_units+pcrows*f.n_A,cap,work)
 }
 _cdid_selref_release(b)
 f.workspace_entries=used
 f.workspace_peak_entries=_cdid_base_cells(1,max((4*maxm*(p+f.n_models+16)+64*maxh,
  16*min((f.N,f.chunk_rows))*(f.n_models+f.max_selection_params+rows((*plan).cohorts)+8)))+256,f.workspace_limit,used)
 f.work_ops=_cdid_base_work(maxm*f.n_models*(p+maxh+1),cap,work)
 /* Complete recipe-keyed selection fitting on post/BOTH-aux union first. */
 for(r=1;r<=f.n_recipes;r++) {
  _cdid_selected_activate_reader(f.reader_owner,r)
  f.recipes[r].fit=cdid_sel_fit(f.reader_owner.reader,f.recipes[r].map,controls.eps,controls.maxit)
  f.n_selection_registry_fits=f.n_selection_registry_fits+1
  f.n_selection_fits=f.n_selection_fits+f.recipes[r].fit.n_blocks
  f.n_selection_blocks=f.n_selection_blocks+f.recipes[r].fit.n_blocks
  f.recipes[r].upstream=cdid_sel_up_map_build(f.recipes[r].fit)
  f.recipes[r].loading=cdid_sel_moments_init(f.recipes[r].fit,rows(f.recipes[r].models))
  for(i=1;i<=rows(f.recipes[r].models);i++) {
   j=f.recipes[r].models[i];cohorts=f.models[j].cohorts;f.models[j].selection_pc=J(rows(cohorts),1,.)
   for(k=1;k<=rows(cohorts);k++)f.models[j].selection_pc[k]=_cdid_scr_pc(f.recipes[r].map,reader_cohort[cohorts[k]],f.models[j].s,f.models[j].t)
  }
 }
 /* Preserve global first-occurrence registry order, not recipe-major callback order. */
 for(j=1;j<=f.n_models;j++)_cdid_scr_callback(f,j,callback,payload)
 f.status="prepared";return(f)
}
real matrix cdid_scr_if_rows(struct cdid_scr_provider scalar f,real colvector rowsF,real colvector models,string scalar scale)
{
 real matrix out,up
 real colvector aa,loc,ix,matchobs,cols_r
 real scalar r,i,j,v,used
 _cdid_scr_owner(f)
 if(f.status!="prepared" | rows(rowsF)<1 | rows(rowsF)>f.chunk_rows | cols(rowsF)!=1 | missing(rowsF) |
  any(rowsF:<1) | any(rowsF:>f.N) | any(rowsF:!=floor(rowsF)) | rows(uniqrows(sort(rowsF,1)))!=rows(rowsF) |
  rows(models)<1 | cols(models)!=1 | missing(models) | any(models:<1) | any(models:>f.n_models) |
  any(models:!=floor(models)) | rows(uniqrows(sort(models,1)))!=rows(models) | (scale!="frame" & scale!="core"))
  _cdid_callback_fail("batch","Invalid selected reference bounded original-F/model request.",498)
 used=_cdid_base_cells(1,16*rows(rowsF)*(f.n_models+f.max_selection_params+rows((*f.plan).cohorts)+8)+256,f.workspace_limit,f.workspace_entries)
 out=J(rows(rowsF),rows(models),0);aa=(*f.plan).F_to_A[rowsF]
 for(r=1;r<=f.n_recipes;r++) {
  loc=vec(selectindex(f.consumer_context[models]:==r));if(rows(loc)==0)continue
  _cdid_selected_activate_reader(f.reader_owner,r)
  up=cdid_sel_eval_loading(f.recipes[r].fit,f.reader_owner.reader,f.recipes[r].loading,rowsF)
  cols_r=J(rows(loc),1,.);for(i=1;i<=rows(loc);i++)cols_r[i]=f.models[models[loc[i]]].loading_col
  out[,loc]=up[,cols_r]
 }
 for(j=1;j<=rows(models);j++) {
  v=models[j];matchobs=cdid_ref_lookup(aa,f.models[v].obs_A);ix=vec(selectindex(matchobs:>0))
  if(rows(ix))out[ix,j]=out[ix,j]+f.N/f.models[v].m:*f.models[v].psi[matchobs[ix]]
 }
 if(scale=="core")out=out:/f.core_to_frame
 if(missing(out))_cdid_callback_fail("numeric","Nonfinite selected reference full original-F influence.",504)
 return(out)
}
struct cdid_scr_weight_context {
 pointer(struct cdid_scr_provider scalar) scalar provider
 real colvector models,signs
}
struct cdid_scr_weight_context scalar cdid_scr_weight_context_make(pointer(struct cdid_scr_provider scalar) scalar provider,
 real colvector models,real colvector signs)
{
 struct cdid_scr_weight_context scalar q
 if(provider==NULL)_cdid_callback_fail("owner","A selected reference provider is required.",498)
 if((*provider).status!="prepared" | rows(models)<1 | rows(models)!=rows(signs) | cols(models)!=1 | cols(signs)!=1 |
  missing(models) | any(models:<1) | any(models:>(*provider).n_models) | any(models:!=floor(models)) | missing(signs) |
  any((signs:!=1):&(signs:!=-1)))_cdid_callback_fail("input","Invalid selected reference signed context.",498)
 q.provider=provider;q.models=models;q.signs=signs;return(q)
}
real matrix cdid_scr_weight(pointer() scalar payload,real scalar lo,real scalar hi)
{
 struct cdid_scr_weight_context scalar q
 real colvector ids,map
 real matrix out
 if(payload==NULL)_cdid_callback_fail("owner","A selected reference weight context is required.",498)
 q=*payload;if(q.provider==NULL)_cdid_callback_fail("owner","A selected reference provider is required.",498)
 if(lo<1 | hi<lo | hi>(*q.provider).N | hi-lo+1>(*q.provider).chunk_rows | lo!=floor(lo) | hi!=floor(hi))
  _cdid_callback_fail("batch","Invalid selected reference weight row range.",498)
 ids=uniqrows(sort(q.models,1));map=cdid_ref_lookup(q.models,ids)
 out=cdid_scr_if_rows(*q.provider,(lo::hi),ids,"core")
 return(out[,map]:*q.signs')
}
end
