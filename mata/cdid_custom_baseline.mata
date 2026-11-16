*! Private selection=none STATIC/dated custom IPW provider, StataNow19.5.
*! Source chained_custom_ipw.R, chained_reference_routes.R, chained_ipw_pipeline.R.
*! Borrows ONE original population/plan/reference/catalog; no nuisance fit or encoder.
version 19.5
mata:
struct cdid_cb_model {
 real scalar raw_index,recipe,g,s,t,m,estimate,callback_if_mean,has_gradient,endpoints_read,n_columns
 string scalar derivative_source,link_id
 real colvector cohorts,obs_A,psi,repeat_mass,pair_mass,pair_count,rates,rates_use,rates_active,b_rate
}
struct cdid_cb_provider {
 string scalar route,status,weight_row_domain
 real scalar n_F,n_A,n_models,n_callbacks,n_endpoint_reads,normalization_count,weight_rows,sample_rank_upper,core_to_frame
 real scalar chunk_rows,workspace_limit,work_limit,borrowed_entries,owned_entries,retained_observer_entries
 real scalar workspace_entries,workspace_peak_entries,callback_entries,replay_entries,work_ops
 pointer(struct cdid_population scalar) scalar population
 pointer(struct cdid_reference_inputs scalar) scalar reference
 pointer(struct cdid_baseline_plan_state scalar) scalar plan
 pointer(struct cdid_reference_catalog scalar) scalar catalog
 real colvector raw_model,post_model,aux_model,points,call_order
 struct cdid_cb_model colvector models
}
real scalar _cdid_cb_catalog_entries(struct cdid_reference_catalog scalar C)
{
 return(length(C.consumers)+length(C.blocks)+length(C.representative)+length(C.consumer_recipe)+
 length(C.consumer_block)+length(C.K_ptr)+length(C.K_index)+length(C.block_K_ptr)+length(C.block_K_index)+
 length(C.raw_consumer)+length(C.post_consumer)+length(C.aux_consumer)+length(C.raw_supported)+length(C.aux_supported)+5)
}
real scalar _cdid_cb_population_entries(struct cdid_population scalar P)
{
 return(length(P.ids.numeric_id)+length(P.ids.row_unit)+length(P.date)+length(P.y)+length(P.calendar)+
 length(P.cohort)+length(P.cohort_known)+length(P.first_row)+length(P.row_count)+length(P.observed_count)+
 length(P.A)+length(P.B)+length(P.frame_only)+length(P.row_order)+length(P.row_start)+length(P.observed_order)+
 length(P.A_to_F)+length(P.F_to_A)+length(P.B_to_F)+length(P.F_to_B)+10)
}
void _cdid_cb_check_owner(struct cdid_cb_provider scalar f)
{
 if(f.population==NULL | f.reference==NULL | f.plan==NULL | f.catalog==NULL)
  _cdid_callback_fail("owner","All original scientific objects must remain available through replay.",498)
 _cdid_ref_check_population(*f.reference)
 if((*f.reference).population!=f.population | (*f.population).N!=f.n_F | (*f.population).n_A!=f.n_A |
  (*f.plan).n_F!=f.n_F | (*f.plan).n_A!=f.n_A | rows(f.raw_model)!=rows((*f.plan).raw))
  _cdid_callback_fail("owner","Borrowed scientific object identity/dimensions changed.",498)
}
/* Shared sparse actual-pair reader returns addresses only; risk membership is
   filtered BEFORE either endpoint value is copied into a callback request. */
real matrix _cdid_cb_observers(struct cdid_population scalar P,
 struct cdid_baseline_plan_state scalar plan,struct cdid_cb_model scalar model,
 real scalar lo,real scalar hi)
{
 real matrix rr
 real colvector keep
 real scalar j
 rr=_cdid_base_pair_rows(P,model.s,model.t,lo,hi)
 keep=J(rows(rr),1,0)
 for(j=1;j<=rows(rr);j++)keep[j]=any(model.cohorts:==plan.a_cohort[rr[j,1]])
 keep=vec(selectindex(keep))
 return(rows(keep) ? rr[keep,] : J(0,3,.))
}
/* Units eligible in the reference eligibility context of a recipe (empty when the
   recipe is unrestricted). */
real colvector _cdid_cb_elig(struct cdid_reference_inputs scalar S,real scalar recipe)
{
 real colvector hit
 real scalar e
 if(cols(S.recipes)<3 | rows(S.elig_index)==0)return(J(0,1,.))
 e=S.recipes[recipe,3]
 if(!(e>0 & e<.))return(J(0,1,.))
 hit=vec(selectindex(S.elig_index:==e))
 if(rows(hit)!=1)_cdid_callback_fail("identity","Unknown reference eligibility context of a design block.",498)
 return(!S.elig_excluded[,hit[1]])
}
real matrix _cdid_cb_X(struct cdid_reference_inputs scalar S,
 struct cdid_baseline_plan_state scalar plan,struct cdid_cb_model scalar model,real colvector ai)
{
 real scalar d
 d=(S.roles[2].dated ? S.recipes[model.recipe,1] : 0)
 return(cdid_ref_design_rows(S,2,d*(S.N+1):+plan.A_to_F[ai]))
}
/* This phase allocates only compact model metadata. Admission of the WHOLE
   catalogue, largest callback handoff and largest bounded replay precedes calls. */
struct cdid_cb_provider scalar cdid_cb_preflight(
 pointer(struct cdid_population scalar) scalar population,
 pointer(struct cdid_reference_inputs scalar) scalar reference,
 pointer(struct cdid_baseline_plan_state scalar) scalar plan,
 pointer(struct cdid_reference_catalog scalar) scalar catalog,
 struct cdid_baseline_controls scalar ctrl,string scalar route)
{
 struct cdid_cb_provider scalar f
 real scalar i,j,k,h,m,rep,hit,slot,p,used,models,nraw,maxh,maxm,work,unit
 real colvector representatives,counts
 if(population==NULL | reference==NULL | plan==NULL | catalog==NULL)
  _cdid_callback_fail("owner","Original population/reference/plan/catalog objects are required.",498)
 if(route!="static" & route!="D60")_cdid_callback_fail("scope","This provider supports static and dated baseline only.",498)
 cdid_ref_activate(*reference)
 if((*reference).population!=population | (*reference).route!="baseline" | (*reference).method!="ipw" |
  ctrl.method!="ipw" | (*plan).n_F!=(*population).N | (*plan).n_A!=(*population).n_A |
  (*catalog).n_raw!=rows((*plan).raw) | (*catalog).M!=(*plan).M |
  rows((*reference).view_recipe)!=rows((*plan).raw))
  _cdid_callback_fail("identity","Require matching original baseline-IPW objects and complete raw catalogue.",498)
 if(ctrl.max_entries<1 | ctrl.max_entries>=. | ctrl.max_entries!=floor(ctrl.max_entries) |
  ctrl.chunk_rows<1 | ctrl.chunk_rows>=. | ctrl.chunk_rows!=floor(ctrl.chunk_rows) |
  ctrl.workcap>. | (ctrl.workcap<. & (ctrl.workcap<1 | ctrl.workcap!=floor(ctrl.workcap))))
  _cdid_callback_fail("control","Require existing explicit workspace, chunk and work controls.",498)
 if((*population).n_A<1)_cdid_callback_fail("support","An actual A population is required.",459)
 if(!(*reference).roles[2].active)_cdid_callback_fail("design","The actual callback design role is inactive.",498)
 slot=(*reference).roles[2].slot
 if(slot<1 | slot>(*reference).n_slots)_cdid_callback_fail("design","Invalid callback design slot.",498)
 if(!(*reference).designs[slot].encoded)_cdid_callback_fail("design","Callback design must be genuinely encoded before preparation.",498)
 p=cols((*reference).designs[slot].value)
 if(p<1)_cdid_callback_fail("design","The chained callback source requires a nonempty actual design.",459)
 if(route=="static") {
  if((*reference).roles[2].dated | (*reference).weight_kind!="static")
   _cdid_callback_fail("scope","Dated designs/fixed-anchor weights require the dated estimator.",498)
  /* Actual static source requires finite F design, without any rank gate. */
  if(rows((*reference).designs[slot].value)!=(*population).N | missing((*reference).designs[slot].value))
   _cdid_callback_fail("design","Static callbacks require their finite original-F design.",459)
 }
 f.population=population;f.reference=reference;f.plan=plan;f.catalog=catalog;f.route=route
 f.n_F=(*population).N;f.n_A=(*population).n_A
 f.weight_rows=(route=="static" ? f.n_A : f.n_F);f.normalization_count=f.n_A
 f.sample_rank_upper=f.weight_rows-1;f.core_to_frame=f.n_F/f.n_A
 f.weight_row_domain=(route=="static" ? "A" : "F")
 f.workspace_limit=ctrl.max_entries;f.work_limit=ctrl.workcap;f.chunk_rows=ctrl.chunk_rows
 f.n_callbacks=f.n_endpoint_reads=0;f.status="preflight"
 nraw=rows((*plan).raw)
 used=_cdid_base_cells(1,_cdid_cb_population_entries(*population),ctrl.max_entries,0)
 used=_cdid_base_cells(1,(*plan).workspace_entries,ctrl.max_entries,used)
 used=_cdid_base_cells(1,(*reference).retained_entries,ctrl.max_entries,used)
 used=_cdid_base_cells(1,_cdid_cb_catalog_entries(*catalog),ctrl.max_entries,used)
 f.borrowed_entries=used
 /* Reserve registry/map/scalar metadata before allocating model maps. */
 used=_cdid_base_cells(1,12*nraw+64,ctrl.max_entries,used)
 f.raw_model=J(nraw,1,0);representatives=J(nraw,1,0);models=0
 if(route=="D60") {
  models=(*catalog).n_consumers
  f.raw_model=(*catalog).raw_consumer
  if(models)representatives[1..models]=(*catalog).representative
 }
 else for(i=1;i<=nraw;i++) {
  if(!(*plan).raw_supported[i])continue
  hit=0
  /* Source static auxiliary searches POST points only and uses first match.
     Each original post occurrence always invokes, even with identical keys. */
  if(i>(*plan).M)for(j=1;j<=(*plan).M;j++) {
   if(all((*plan).raw[i,]:==(*plan).raw[j,]))if(_cdid_base_same_K(*plan,i,j)) {hit=f.raw_model[j];break;}
  }
  if(hit==0){models++;representatives[models]=i;hit=models;}
  f.raw_model[i]=hit
 }
 f.n_models=models
 maxh=maxm=0;f.retained_observer_entries=0
 work=_cdid_base_work(0,ctrl.workcap,(*plan).work_ops)
 unit=2*(1+ceil(ln(rows((*plan).calendar)+1)/ln(2)))
 /* Admit every model before allocating retained model vectors or invoking any
    callback. Support counts are already outcome-free plan metadata. */
 for(i=1;i<=models;i++) {
  rep=representatives[i];h=1+(*plan).raw_K_ptr[rep+1]-(*plan).raw_K_ptr[rep]
  counts=(*plan).raw_group_count[(*plan).raw_count_ptr[rep]..(*plan).raw_count_ptr[rep+1]-1]
  m=sum(counts)
  if(m<2 | any(counts:<1))_cdid_callback_fail("support","Every callback cohort requires actual observed pair support.",459)
  used=_cdid_base_cells(1,2*m+9*h+24,ctrl.max_entries,used)
  f.retained_observer_entries=f.retained_observer_entries+2*m
  maxh=max((maxh,h));maxm=max((maxm,m))
  work=_cdid_base_work(4*f.n_A*(unit+h)+8*m*(p+h+1),ctrl.workcap,work)
 }
 f.callback_entries=4*maxm*(p+16)+64*maxh+128
 f.replay_entries=16*min((ctrl.chunk_rows,f.n_F))*(max((1,nraw))+8)+64*maxh+128
 f.workspace_peak_entries=_cdid_base_cells(1,max((f.callback_entries,f.replay_entries)),ctrl.max_entries,used)
 f.models=cdid_cb_model(models,1);f.points=J(models,1,.)
 f.call_order=(models ? representatives[1..models] : J(0,1,.))
 f.post_model=((*plan).M ? f.raw_model[1..(*plan).M] : J(0,1,.))
 f.aux_model=((*plan).Q ? f.raw_model[(*plan).M+1..nraw] : J(0,1,.))
 for(i=1;i<=models;i++) {
  rep=representatives[i]
  f.models[i].raw_index=rep;f.models[i].recipe=(*reference).view_recipe[rep]
  f.models[i].g=(*plan).raw[rep,1];f.models[i].s=(*plan).raw[rep,2];f.models[i].t=(*plan).raw[rep,3]
  f.models[i].cohorts=f.models[i].g\(*plan).raw_K_index[(*plan).raw_K_ptr[rep]..(*plan).raw_K_ptr[rep+1]-1]
  counts=(*plan).raw_group_count[(*plan).raw_count_ptr[rep]..(*plan).raw_count_ptr[rep+1]-1]
  f.models[i].m=sum(counts);f.models[i].n_columns=p;f.models[i].pair_count=counts
  f.models[i].link_id=(rep<=(*plan).M ? "post:" : "reporting:")+strofreal(rep<=(*plan).M ? rep : rep-(*plan).M)
 }
 f.owned_entries=used-f.borrowed_entries;f.workspace_entries=used
 f.work_ops=work
 return(f)
}
void _cdid_cb_fit_model(struct cdid_cb_provider scalar f,real scalar j,
 pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_callback_request scalar q
 struct cdid_callback_result scalar result
 real matrix rr,clip
 real colvector ix,gi,ww,ai,co,allmass,pm,counts,riskw,el
 real scalar lo,hi,k,m,h,at,n,slot
 q=cdid_callback_request_init();m=f.models[j].m;h=rows(f.models[j].cohorts)
 q.mode="chained";q.route="baseline_"+(f.route=="D60" ? "dated" : f.route);q.link_id=f.models[j].link_id
 q.panel=q.inffunc=1;q.gradient_required=0;q.N=f.n_F;q.n_A=f.n_A
 q.call_index=j;q.raw_index=f.models[j].raw_index;q.model_index=j
 q.g=(*f.plan).cohorts[f.models[j].g];q.K=(*f.plan).cohorts[f.models[j].cohorts[2..h]]
 q.s=(*f.plan).calendar[f.models[j].s];q.t=(*f.plan).calendar[f.models[j].t];q.recipe=f.models[j].recipe
 q.rows_A=q.rows_F=q.y0=q.y1=q.D=q.i_weights=J(m,1,.)
 slot=(*f.reference).roles[2].slot;q.covariate_names=(*f.reference).designs[slot].columns
 q.covariates=J(m,f.models[j].n_columns,.)
 gi=riskw=J(m,1,.);allmass=pm=counts=J(h,1,0)
 /* A restricted recipe reads only the units eligible in its context. */
 el=_cdid_cb_elig(*f.reference,f.models[j].recipe)
 /* Denominator uses actual fixed-anchor W on all A rows of each cohort. */
 for(k=1;k<=h;k++) {
  ai=vec(selectindex((*f.plan).a_cohort:==f.models[j].cohorts[k]))
  if(rows(el) & rows(ai))ai=select(ai,el[(*f.plan).A_to_F[ai]])
  ww=cdid_ref_weight_rows(*f.reference,f.models[j].recipe,(*f.plan).A_to_F[ai])
  allmass[k]=sum(ww)
 }
 at=0
 for(lo=1;lo<=f.n_A;lo=lo+f.chunk_rows) {
  hi=min((f.n_A,lo+f.chunk_rows-1))
  rr=_cdid_cb_observers(*f.population,*f.plan,f.models[j],lo,hi)
  if(rows(el) & rows(rr))rr=select(rr,el[(*f.plan).A_to_F[rr[,1]]])
  n=rows(rr)
  if(n==0)continue
  if(at+n>m)_cdid_callback_fail("identity","Actual callback observer count exceeds admitted support.",498)
  ix=((at+1)::(at+n));ai=rr[,1]
  q.rows_A[ix]=ai;q.rows_F[ix]=(*f.plan).A_to_F[ai]
  /* Read BOTH original endpoints, not an artificial dy/zero pair. */
  q.y0[ix]=(*f.population).y[rr[,2]];q.y1[ix]=(*f.population).y[rr[,3]]
  q.covariates[ix,]=_cdid_cb_X(*f.reference,*f.plan,f.models[j],ai)
  riskw[ix]=cdid_ref_weight_rows(*f.reference,f.models[j].recipe,(*f.plan).A_to_F[ai])
  co=(*f.plan).a_cohort[ai]
  for(k=1;k<=h;k++) {
   ww=vec(selectindex(co:==f.models[j].cohorts[k]))
   if(rows(ww))gi[ix[ww]]=J(rows(ww),1,k)
  }
  at=at+n
 }
 if(at!=m | missing(q.y0) | missing(q.y1) | missing(q.y1-q.y0) | missing(q.covariates) |
  missing(riskw) | any(riskw:<=0) | missing(gi))
  _cdid_callback_fail("input_numeric","Callback input support, endpoints, design or weights are invalid.",504)
 for(k=1;k<=h;k++) {
  ix=vec(selectindex(gi:==k));counts[k]=rows(ix);pm[k]=sum(riskw[ix])
 }
 if(any(counts:!=f.models[j].pair_count) | missing(allmass) | any(allmass:<=0) | missing(pm) | any(pm:<=0))
  _cdid_callback_fail("support","Actual callback cohort support disagrees with the original plan.",459)
 f.models[j].repeat_mass=allmass;f.models[j].pair_mass=pm;f.models[j].rates=pm:/allmass
 clip=cdid_ipw_rate_values(f.models[j].rates)
 f.models[j].rates_use=clip[,1];f.models[j].rates_active=clip[,2]
 q.D=(gi:==1);q.i_weights=riskw:/clip[gi,1]
 if(missing(q.i_weights) | any(q.i_weights:<=0))_cdid_callback_fail("input_numeric","Invalid actual inverse-rate callback weights.",504)
 /* Caller owns payload; no dereference/copy/evaluation before actual invocation. */
 result=cdid_callback_invoke_chained(callback,q,payload)
 f.n_callbacks=f.n_callbacks+1
 f.models[j].estimate=result.ATT[1,1];f.points[j]=f.models[j].estimate
 f.models[j].callback_if_mean=mean(result.psi);f.models[j].has_gradient=result.has_gradient
 f.models[j].derivative_source=(result.has_gradient ? "explicit_log_weight_gradient" : "identity_psi_over_m (existing baseline rule)")
 f.models[j].obs_A=q.rows_A;f.models[j].psi=result.psi
 f.models[j].b_rate=cdid_callback_rate_loading(result,gi,f.models[j].rates)
 f.models[j].endpoints_read=2*m;f.n_endpoint_reads=f.n_endpoint_reads+2*m
 /* q, optional g and observer scratch die here. Only raw psi/indices and compact
    cohort/loadings persist. No callback request stores a private owner pointer. */
}
struct cdid_cb_provider scalar cdid_cb_prepare(
 pointer(struct cdid_population scalar) scalar population,
 pointer(struct cdid_reference_inputs scalar) scalar reference,
 pointer(struct cdid_baseline_plan_state scalar) scalar plan,
 pointer(struct cdid_reference_catalog scalar) scalar catalog,
 struct cdid_baseline_controls scalar ctrl,string scalar route,
 pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_cb_provider scalar f
 real scalar j
 if(callback==NULL)_cdid_callback_fail("input","A scientific callback is required before preparation.",498)
 f=cdid_cb_preflight(population,reference,plan,catalog,ctrl,route)
 for(j=1;j<=f.n_models;j++)_cdid_cb_fit_model(f,j,callback,payload)
 f.status=(f.n_models ? "prepared" : "no_supported_consumers")
 return(f)
}
/* An original-F batch in caller order. All replay uses retained psi/indices;
   there are no callback calls, nuisance fits or per-chunk primitive refits. */
real matrix cdid_cb_if_rows(struct cdid_cb_provider scalar f,real colvector rowsF,
 real colvector models,string scalar scale)
{
 real matrix out,rr
 real colvector aa,matchobs,ix,co,ww,S,obspos,el
 real scalar j,k,v,m,i,n,B
 _cdid_cb_check_owner(f);B=rows(rowsF)
 if(cols(rowsF)!=1 | B<1 | B>f.chunk_rows | missing(rowsF) | any(rowsF:<1) | any(rowsF:>f.n_F) |
  any(rowsF:!=floor(rowsF)) | rows(uniqrows(sort(rowsF,1)))!=B | cols(models)!=1 | rows(models)<1 |
  missing(models) | any(models:<1) | any(models:>f.n_models) | any(models:!=floor(models)) |
  rows(uniqrows(sort(models,1)))!=rows(models) | (scale!="frame" & scale!="core"))
  _cdid_callback_fail("batch","Require unique explicit original-F/model batches and frame/core scale.",498)
 n=_cdid_base_cells(1,16*B*(rows(models)+8)+128,f.workspace_limit,f.workspace_entries)
 out=J(B,rows(models),0);aa=(*f.plan).F_to_A[rowsF]
 for(j=1;j<=rows(models);j++) {
  v=models[j];m=f.models[v].m
  el=_cdid_cb_elig(*f.reference,f.models[v].recipe)
  /* Binary matching references retained observed indices directly, no psi copy. */
  matchobs=cdid_ref_lookup(aa,f.models[v].obs_A)
  ix=vec(selectindex(matchobs:>0))
  if(rows(ix))out[ix,j]=(f.n_F/m):*f.models[v].psi[matchobs[ix]]
  for(k=1;k<=rows(f.models[v].cohorts);k++) {
   if(f.models[v].b_rate[k]==0)continue
   co=J(B,1,0)
   for(i=1;i<=B;i++)if(aa[i]>0)co[i]=((*f.plan).a_cohort[aa[i]]==f.models[v].cohorts[k])
   if(rows(el))co=co:*el[rowsF]
   ix=vec(selectindex(co));if(rows(ix)==0)continue
   ww=cdid_ref_weight_rows(*f.reference,f.models[v].recipe,rowsF[ix])
   S=(matchobs[ix]:>0)
   out[ix,j]=out[ix,j]+f.n_F:*ww:*(S:-f.models[v].rates[k]):/f.models[v].repeat_mass[k]:*f.models[v].b_rate[k]
  }
 }
 if(scale=="core")out=out:/f.core_to_frame
 if(missing(out))_cdid_callback_fail("numeric","Nonfinite complete callback IF replay.",504)
 return(out)
}
/* Source GMM row domain is route-specific: static centers literal A/core
   rows; dated centers complete original-F/core rows with n_A normalization.
   Full-F IF pullback is a SEPARATE explicit consumer below, never synthetic
   zero padding inserted into static GMM's centering population. */
real matrix cdid_cb_weight_chunk(struct cdid_cb_provider scalar f,real scalar lo,real scalar hi,
 real colvector models)
{
 real colvector rowsF
 if(lo<1 | hi<lo | hi>f.weight_rows | lo!=floor(lo) | hi!=floor(hi))
  _cdid_callback_fail("batch","Invalid source weight-domain row range.",498)
 rowsF=(lo::hi)
 if(f.route=="static")rowsF=(*f.plan).A_to_F[rowsF]
 return(cdid_cb_if_rows(f,rowsF,models,"core"))
}
real matrix cdid_cb_fullF_core_chunk(struct cdid_cb_provider scalar f,real scalar lo,real scalar hi,
 real colvector models)
{
 if(lo<1 | hi<lo | hi>f.n_F | lo!=floor(lo) | hi!=floor(hi))
  _cdid_callback_fail("batch","Invalid explicit original-F IF row range.",498)
 return(cdid_cb_if_rows(f,(lo::hi),models,"core"))
}
struct cdid_cb_weight_context {
 pointer(struct cdid_cb_provider scalar) scalar provider
 real colvector models,signs
}
struct cdid_cb_weight_context scalar cdid_cb_weight_context_make(
 pointer(struct cdid_cb_provider scalar) scalar provider,real colvector models,real colvector signs)
{
 struct cdid_cb_weight_context scalar q
 real scalar j
 if(provider==NULL)_cdid_callback_fail("owner","Callback provider is required.",498)
 if(rows(models)!=rows(signs) | rows(models)<1 | cols(models)!=1 | cols(signs)!=1 | missing(signs) |
  any((signs:!=1):&(signs:!=-1)) | missing(models) | any(models:<1) | any(models:>(*provider).n_models) |
  any(models:!=floor(models)))_cdid_callback_fail("input","Invalid signed model composition.",498)
 j=_cdid_base_cells(2,rows(models),(*provider).workspace_limit,(*provider).workspace_entries)
 q.provider=provider;q.models=models;q.signs=signs
 return(q)
}
real matrix cdid_cb_weight(pointer() scalar context,real scalar lo,real scalar hi)
{
 struct cdid_cb_weight_context scalar q
 real colvector ids,map
 real matrix z
 real scalar j
 if(context==NULL)_cdid_callback_fail("owner","Callback weight context is required.",498)
 q=*context
 j=_cdid_base_cells(1,2*rows(q.models)+16*(hi-lo+1)*(rows(q.models)+8)+128,(*q.provider).workspace_limit,(*q.provider).workspace_entries)
 ids=uniqrows(sort(q.models,1));map=cdid_ref_lookup(q.models,ids)
 z=cdid_cb_weight_chunk(*q.provider,lo,hi,ids)
 return(z[,map]:*q.signs')
}
real matrix cdid_cb_post_core(pointer() scalar context,real scalar lo,real scalar hi)
{
 struct cdid_cb_provider scalar f
 real colvector ids,map
 real matrix z
 if(context==NULL)_cdid_callback_fail("owner","Callback provider is required.",498)
 f=*context;ids=uniqrows(sort(f.post_model,1));map=cdid_ref_lookup(f.post_model,ids)
 z=cdid_cb_weight_chunk(f,lo,hi,ids)
 return(z[,map])
}
end
