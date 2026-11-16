*! Private static conditional_mtar custom-IPW provider.
*! Fixed R cdid/, chained_custom_selection.R32-226.
*! Same fixed cdid_callback SDK; no built-in treatment/OR constructor.
version 19.5
mata:
struct cdid_sc_model {
 real scalar representative,target,s,t,m,estimate,callback_if_mean
 real colvector cohorts,selection_pc,obs_A,psi,gradient
}
struct cdid_sc_provider {
 pointer(struct cdid_population scalar) scalar population
 pointer(struct cdid_sel_reader scalar) scalar reader
 pointer(struct cdid_baseline_plan_state scalar) scalar plan
 pointer(struct cdid_reference_catalog scalar) scalar catalog
 struct cdid_sel_cell_map scalar selection_map
 struct cdid_sel_fit_result scalar selection_fit
 struct cdid_sel_up_map scalar upstream_map
 struct cdid_sel_moments scalar loading
 struct cdid_sc_model colvector models
 real matrix callback_X
 string rowvector callback_columns
 real colvector points,raw_model,post_model,aux_model,reader_cohort
 real scalar N,n_A,n_models,n_callbacks,n_selection_fits,n_treatment_fits,n_endpoint_reads
 real scalar weight_rows,normalization_count,sample_rank_upper,core_to_frame,chunk_rows
 real scalar workspace_entries,workspace_peak_entries,workspace_limit,work_ops
 string scalar weight_row_domain,status
}
real scalar _cdid_sc_pc(struct cdid_sel_cell_map scalar map,real scalar k,real scalar s,real scalar t)
{
 real colvector ix
 ix=vec(selectindex((map.pair_cohorts[,1]:==k):&(map.pair_cohorts[,2]:==s):&(map.pair_cohorts[,3]:==t)))
 if(rows(ix)!=1)_cdid_callback_fail("selection_identity","Required actual pair/cohort selection cell is not uniquely mapped.",498)
 return(ix[1])
}
real matrix _cdid_sc_observers(struct cdid_sc_provider scalar f,real scalar j,real scalar lo,real scalar hi)
{
 real matrix rr
 real colvector take
 real scalar i
 rr=_cdid_base_pair_rows(*f.population,f.models[j].s,f.models[j].t,lo,hi);take=J(rows(rr),1,0)
 for(i=1;i<=rows(rr);i++)take[i]=any(f.models[j].cohorts:==(*f.plan).a_cohort[rr[i,1]])
 take=vec(selectindex(take))
 return(rows(take) ? rr[take,] : J(0,3,.))
}
void _cdid_sc_callback(struct cdid_sc_provider scalar f,real scalar j,
 pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_callback_request scalar q
 struct cdid_callback_result scalar r
 struct cdid_sel_pred scalar pred
 struct cdid_sel_moments scalar one
 real matrix rr
 real colvector gi,co,ix,sub,units
 real scalar lo,hi,at,n,k,m,h
 q=cdid_callback_request_init();m=f.models[j].m;h=rows(f.models[j].cohorts)
 q.mode="chained";q.route="conditional_mtar_static";q.link_id="model:"+strofreal(j)
 q.panel=q.inffunc=q.gradient_required=1;q.N=f.N;q.n_A=f.n_A
 q.call_index=q.model_index=j;q.raw_index=f.models[j].representative;q.recipe=1
 q.g=(*f.plan).cohorts[f.models[j].target];q.K=(*f.plan).cohorts[f.models[j].cohorts[2..h]]
 q.s=(*f.plan).calendar[f.models[j].s];q.t=(*f.plan).calendar[f.models[j].t]
 q.rows_A=q.rows_F=q.y0=q.y1=q.D=q.i_weights=J(m,1,.)
 q.covariates=J(m,cols(f.callback_X),.);q.covariate_names=f.callback_columns;gi=J(m,1,.)
 at=0
 for(lo=1;lo<=f.n_A;lo=lo+f.chunk_rows) {
  hi=min((f.n_A,lo+f.chunk_rows-1));rr=_cdid_sc_observers(f,j,lo,hi);n=rows(rr)
  if(n==0)continue
  if(at+n>m)_cdid_callback_fail("identity","Actual selected observer count exceeds admitted support.",498)
  ix=((at+1)::(at+n));q.rows_A[ix]=rr[,1];q.rows_F[ix]=(*f.plan).A_to_F[rr[,1]]
  q.y0[ix]=(*f.population).y[rr[,2]];q.y1[ix]=(*f.population).y[rr[,3]]
  q.covariates[ix,]=f.callback_X[rr[,1],];co=(*f.plan).a_cohort[rr[,1]]
  for(k=1;k<=h;k++) {
   sub=vec(selectindex(co:==f.models[j].cohorts[k]))
   if(rows(sub)) {
    gi[ix[sub]]=J(rows(sub),1,k)
    pred=cdid_sel_pred_pair(f.selection_fit,*f.reader,f.models[j].selection_pc[k],q.rows_F[ix[sub]])
    q.i_weights[ix[sub]]=pred.sel_weight
   }
  }
  at=at+n
 }
 if(at!=m | missing(q.y0) | missing(q.y1) | missing(q.y1-q.y0) | missing(q.covariates) |
  missing(q.i_weights) | any(q.i_weights:<=0) | missing(gi))
  _cdid_callback_fail("selected_input","Actual selected callback support/endpoints/X/weights are invalid.",504)
 q.D=(gi:==1)
 /* The SDK validates all supplied output first, THEN requires explicit g.
    This refusal is after selection fitting and actual callback. */
 r=cdid_callback_invoke_chained(callback,q,payload);f.n_callbacks=f.n_callbacks+1
 f.points[j]=r.ATT[1,1];f.models[j].estimate=r.ATT[1,1]
 f.models[j].callback_if_mean=mean(r.psi);f.models[j].obs_A=q.rows_A
 f.models[j].psi=r.psi;f.models[j].gradient=r.log_weight_gradient
 one=cdid_sel_moments_init(f.selection_fit,1)
 for(k=1;k<=h;k++) {
  sub=vec(selectindex(gi:==k));if(rows(sub)==0)continue
  /* Literal TOTAL g, no /m or /N; primitive computes LD' g. */
  cdid_sel_add_moments(one,f.selection_fit,*f.reader,f.upstream_map,
   f.models[j].selection_pc[k],q.rows_F[sub],r.log_weight_gradient[sub])
 }
 f.loading.G_theta[,j]=one.G_theta;f.loading.G_ak[,j]=one.G_ak
 f.n_endpoint_reads=f.n_endpoint_reads+2*m
}
struct cdid_sc_provider scalar cdid_sc_prepare(
 pointer(struct cdid_population scalar) scalar population,
 pointer(struct cdid_sel_reader scalar) scalar reader,
 pointer(struct cdid_baseline_plan_state scalar) scalar plan,
 pointer(struct cdid_reference_catalog scalar) scalar catalog,
 real matrix X_A,string rowvector columns,struct cdid_sel_decl scalar declaration,
 struct cdid_selected_controls scalar control,pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_sc_provider scalar f
 struct cdid_sel_budget scalar budget
 real scalar i,j,k,r,nc,m,maxm,maxh,used,retained,work,pcrows,at,np
 real colvector cohorts,count
 real matrix specs
 /* Existing public-source scope; do not force inactive data to decide this. */
 if(!declaration.enabled | declaration.canonical_assumption!="conditional_mtar" | control.est_method!="ipw")
  _cdid_callback_fail("scope","This custom adapter requires conditional_mtar with selection IPW.",498)
 if(declaration.selection_model=="independent_period_logit") {
  if(declaration.selection_universe!="F" | declaration.has_rho | declaration.has_just)
   _cdid_callback_fail("scope","Date-product custom selection requires F with no rho/justification.",459)
 }
 if(callback==NULL | population==NULL | reader==NULL | plan==NULL | catalog==NULL)
  _cdid_callback_fail("owner","Callback/original population/selection reader/plan/catalog required.",498)
 if((*reader).selection_model!=declaration.selection_model | (*reader).selection_pooling!=declaration.selection_pooling |
  (*reader).selection_universe!=declaration.selection_universe | (*reader).N!=(*population).N |
  (*plan).n_F!=(*population).N | (*plan).n_A!=(*population).n_A |
  (*catalog).n_raw!=rows((*plan).raw) | any((*catalog).consumer_recipe:!=1))
  _cdid_callback_fail("identity","This static selected-model adapter requires matching original objects and one actual static design block.",498)
 if(rows(X_A)!=(*population).n_A | cols(X_A)<1 | cols(columns)!=cols(X_A) | missing(X_A))
  _cdid_callback_fail("callback_design","Finite actual callback X is required on A only; no global rank gate is imposed.",459)
 if(!(*reader).has_pop | (*reader).kind!="population" | (*reader).pop.R!=(*population).R |
  rows((*reader).cohorts)<1 | rows((*plan).A_to_F)!=(*population).n_A |
  any((*reader).pop.ids.row_unit:!=(*population).ids.row_unit) |
  any((*reader).pop.ids.numeric_id:!=(*population).ids.numeric_id) |
  any((*reader).pop.ids.string_id_value:!=(*population).ids.string_id_value) |
  any((*reader).pop.date:!=(*population).date) | any((*reader).pop.y:!=(*population).y) |
  any((*reader).A:!=(*population).A) | any((*reader).B:!=(*population).B) |
  any((*reader).pop.cohort:!=(*population).cohort) |
  any((*reader).pop.cohort_known:!=(*population).cohort_known) |
  any((*plan).A_to_F:!=(*population).A_to_F) | any((*plan).F_to_A:!=(*population).F_to_A) |
  rows((*reader).all_times)!=rows((*plan).calendar))
  _cdid_callback_fail("identity","Selection reader, population and plan must describe the same original identities/endpoints.",498)
 if(any((*reader).all_times:!=(*plan).calendar))
  _cdid_callback_fail("identity","Selection and scientific calendars must match.",498)
 if(control.selection_model!=declaration.selection_model | control.selection_pooling!=declaration.selection_pooling |
  control.selection!="conditional_mtar" | control.eps<=0 | control.eps>=1 | control.eps>=. |
  control.maxit<1 | control.maxit>=. | control.maxit!=floor(control.maxit) |
  control.max_chunk_rows<1 | control.max_chunk_rows>=. | control.max_chunk_rows!=floor(control.max_chunk_rows) |
  control.max_workspace_doubles<1 | control.max_workspace_doubles>=. |
  control.max_workspace_doubles!=floor(control.max_workspace_doubles))
  _cdid_callback_fail("control","Existing explicit selected chunk/workspace controls are required.",498)
 f.population=population;f.reader=reader;f.plan=plan;f.catalog=catalog
 f.N=(*population).N;f.n_A=(*population).n_A;f.n_models=(*catalog).n_consumers
 f.weight_rows=f.N;f.normalization_count=f.n_A;f.sample_rank_upper=f.N-1;f.core_to_frame=f.N/f.n_A
 f.weight_row_domain="F";f.chunk_rows=control.max_chunk_rows;f.workspace_limit=control.max_workspace_doubles
 f.n_callbacks=f.n_treatment_fits=f.n_endpoint_reads=0;f.status="preflight"
 f.raw_model=(*catalog).raw_consumer;f.post_model=(*catalog).post_consumer;f.aux_model=(*catalog).aux_consumer
 f.callback_X=X_A;f.callback_columns=columns
 f.reader_cohort=cdid_ref_lookup((*plan).cohorts,(*reader).cohorts)
 if(any(f.reader_cohort:<1))
  _cdid_callback_fail("identity","Every scientific catalogue cohort must map to the selection dictionary.",498)
 /* Borrowed original scientific owners, selection reader/profile arrays and
    entire source registry remain charged. Opaque callback workspace is not RSS. */
 used=_cdid_base_cells(1,(*plan).workspace_entries+_cdid_cb_population_entries(*population)+_cdid_cb_catalog_entries(*catalog),f.workspace_limit,0)
 used=_cdid_base_cells(1,_cdid_cb_population_entries((*reader).pop)+length((*reader).Z)+18*f.N+length((*reader).date_responses)+length((*reader).pair_responses)+length(X_A)+256,f.workspace_limit,used)
 pcrows=0;maxm=maxh=0
 for(i=1;i<=f.n_models;i++)pcrows=pcrows+1+(*catalog).K_ptr[i+1]-(*catalog).K_ptr[i]
 used=_cdid_base_cells(1,8*pcrows+12*rows((*plan).raw)+32*f.n_models,f.workspace_limit,used)
 specs=J(pcrows,3,.);at=0
 for(i=1;i<=f.n_models;i++) {
  cohorts=(*catalog).consumers[i,1]\(*catalog).K_index[(*catalog).K_ptr[i]..((*catalog).K_ptr[i+1]-1)]
  for(k=1;k<=rows(cohorts);k++){at=at+1;specs[at,]=(f.reader_cohort[cohorts[k]],(*catalog).consumers[i,2],(*catalog).consumers[i,3]);}
 }
 f.selection_map=cdid_sel_partition(*reader,specs)
 budget=cdid_sel_budget_check(*reader,f.selection_map,(control.max_work>0 ? control.max_work : .))
 /* Actual selection fits occur on post+BOTH-aux union BEFORE callbacks. */
 used=_cdid_base_cells(1,4*budget.compact_state_doubles+budget.stream_workspace_doubles+2*(budget.total_params+rows((*reader).cohorts))*f.n_models,f.workspace_limit,used)
 for(i=1;i<=f.n_models;i++) {
  r=(*catalog).representative[i]
  count=(*plan).raw_group_count[(*plan).raw_count_ptr[r]..((*plan).raw_count_ptr[r+1]-1)]
  m=sum(count);nc=rows(count)
  if(m<2 | any(count:<1))_cdid_callback_fail("support","Every selected callback cohort requires actual pair support.",459)
  used=_cdid_base_cells(1,3*m+2*nc+24,f.workspace_limit,used);maxm=max((maxm,m));maxh=max((maxh,nc))
 }
 f.workspace_entries=used
 f.workspace_peak_entries=_cdid_base_cells(1,max((4*maxm*(cols(X_A)+f.n_models+16)+64*maxh,
  16*min((f.N,f.chunk_rows))*(f.n_models+budget.total_params+rows((*reader).cohorts)+8)))+256,f.workspace_limit,used)
 f.work_ops=_cdid_base_work(budget.work_units+pcrows*f.n_A+maxm*f.n_models*(cols(X_A)+maxh+1),(control.max_work>0 ? control.max_work : .),(*plan).work_ops)
 f.selection_fit=cdid_sel_fit(*reader,f.selection_map,control.eps,control.maxit)
 f.n_selection_fits=f.selection_fit.n_blocks
 f.upstream_map=cdid_sel_up_map_build(f.selection_fit)
 f.loading=cdid_sel_moments_init(f.selection_fit,f.n_models)
 f.models=cdid_sc_model(f.n_models,1);f.points=J(f.n_models,1,.)
 for(i=1;i<=f.n_models;i++) {
  r=(*catalog).representative[i];f.models[i].representative=r
  f.models[i].target=(*catalog).consumers[i,1];f.models[i].s=(*catalog).consumers[i,2];f.models[i].t=(*catalog).consumers[i,3]
  f.models[i].cohorts=f.models[i].target\(*catalog).K_index[(*catalog).K_ptr[i]..((*catalog).K_ptr[i+1]-1)]
  count=(*plan).raw_group_count[(*plan).raw_count_ptr[r]..((*plan).raw_count_ptr[r+1]-1)]
  f.models[i].m=sum(count);f.models[i].selection_pc=J(rows(count),1,.)
  for(k=1;k<=rows(count);k++)f.models[i].selection_pc[k]=_cdid_sc_pc(f.selection_map,f.reader_cohort[f.models[i].cohorts[k]],f.models[i].s,f.models[i].t)
 }
 for(i=1;i<=f.n_models;i++)_cdid_sc_callback(f,i,callback,payload)
 f.status="prepared"
 return(f)
}
real matrix cdid_sc_if_rows(struct cdid_sc_provider scalar f,real colvector rowsF,real colvector models,string scalar scale)
{
 real matrix out,all_up
 real colvector aa,matchobs,ix
 real scalar j,v,used
 if(f.status!="prepared" | f.population==NULL | f.reader==NULL | f.plan==NULL)
  _cdid_callback_fail("owner","Prepared selected callback objects are required.",498)
 if(rows(rowsF)<1 | rows(rowsF)>f.chunk_rows | cols(rowsF)!=1 | missing(rowsF) | any(rowsF:<1) | any(rowsF:>f.N) |
  any(rowsF:!=floor(rowsF)) | rows(uniqrows(sort(rowsF,1)))!=rows(rowsF) | rows(models)<1 | cols(models)!=1 |
  missing(models) | any(models:<1) | any(models:>f.n_models) | any(models:!=floor(models)) |
  rows(uniqrows(sort(models,1)))!=rows(models) | (scale!="frame" & scale!="core"))
  _cdid_callback_fail("batch","Invalid bounded original-F/model batch or scale.",498)
 used=_cdid_base_cells(1,16*rows(rowsF)*(f.n_models+f.selection_fit.total_params+rows(f.selection_fit.cohorts)+8)+256,f.workspace_limit,f.workspace_entries)
 all_up=cdid_sel_eval_loading(f.selection_fit,*f.reader,f.loading,rowsF)
 out=all_up[,models];aa=(*f.plan).F_to_A[rowsF]
 for(j=1;j<=rows(models);j++) {
  v=models[j];matchobs=cdid_ref_lookup(aa,f.models[v].obs_A);ix=vec(selectindex(matchobs:>0))
  if(rows(ix))out[ix,j]=out[ix,j]+f.N/f.models[v].m:*f.models[v].psi[matchobs[ix]]
 }
 if(scale=="core")out=out:/f.core_to_frame
 if(missing(out))_cdid_callback_fail("numeric","Nonfinite full original-F selected callback IF.",504)
 return(out)
}
struct cdid_sc_weight_context {
 pointer(struct cdid_sc_provider scalar) scalar provider
 real colvector models,signs
}
struct cdid_sc_weight_context scalar cdid_sc_weight_context_make(pointer(struct cdid_sc_provider scalar) scalar provider,
 real colvector models,real colvector signs)
{
 struct cdid_sc_weight_context scalar q
 real scalar used
 if(provider==NULL)_cdid_callback_fail("owner","Selected callback provider is required.",498)
 if((*provider).status!="prepared" | rows(models)<1 | rows(models)!=rows(signs) | cols(models)!=1 | cols(signs)!=1 | missing(signs) |
  missing(models) | any(models:<1) | any(models:>(*provider).n_models) | any(models:!=floor(models)) |
  any((signs:!=1):&(signs:!=-1)))_cdid_callback_fail("input","Invalid selected callback signed context.",498)
 used=_cdid_base_cells(2,rows(models),(*provider).workspace_limit,(*provider).workspace_entries)
 q.provider=provider;q.models=models;q.signs=signs
 return(q)
}
real matrix cdid_sc_weight(pointer() scalar context,real scalar lo,real scalar hi)
{
 struct cdid_sc_weight_context scalar q
 real colvector ids,map
 real matrix z
 real scalar used
 if(context==NULL)_cdid_callback_fail("owner","Selected callback weight context required.",498)
 q=*context
 if(q.provider==NULL)_cdid_callback_fail("owner","Selected provider object is required.",498)
 if(lo<1 | hi<lo | hi>(*q.provider).weight_rows | lo>=. | hi>=. | lo!=floor(lo) | hi!=floor(hi) |
  hi-lo+1>(*q.provider).chunk_rows)_cdid_callback_fail("batch","Invalid selected weight-provider range.",498)
 used=_cdid_base_cells(1,2*rows(q.models)+16*(hi-lo+1)*(rows(q.models)+8)+128,
  (*q.provider).workspace_limit,(*q.provider).workspace_entries)
 ids=uniqrows(sort(q.models,1));map=cdid_ref_lookup(q.models,ids)
 z=cdid_sc_if_rows(*q.provider,(lo::hi),ids,"core")
 return(z[,map]:*q.signs')
}
end
