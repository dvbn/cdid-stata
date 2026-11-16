*! Pre-GMM baseline/selected/imputation/augmented conditional producer, complete-family executor and shared signs.
version 19.5
mata:
struct cdid_pt_chain {
 struct cdid_public_bundle scalar bundle
 struct cdid_public_binding_control scalar control
 struct cdid_cpt_baseline scalar process
 struct cdid_cpt_selected scalar selected
 struct cdid_cpti scalar cpti
 real colvector cpti_design_index
 real scalar selected_active,cpti_active,n_cells,workspace_entries,work_ops
 string scalar status,reason
 struct cdid_cpt_budget scalar budget
 real colvector cmap
 string colvector cluster_labels
 real scalar C,cap,bytes,draws,alpha,input_rows
}
/* Literal public declaration predicates precede data/formula evaluation:
   R/chained_selection_consumers.R normalize_selection_declaration and
   R/chained_cp04_public.R cp04_validate_declaration. Lower mode validation
   still checks actual U/rho; it is not a substitute for this source-order gate. */
void _cdid_pt_selection_declaration(string scalar selection,real scalar hasq,
 string scalar universe,real scalar hasrho,real scalar hasjust,string scalar model,
 string scalar pooling,string scalar justification)
{
 struct cdid_sel_decl scalar d
 d=cdid_sel_norm_decl(selection,hasq,universe,hasrho,hasjust,model,pooling,justification)
 if(!d.enabled)return
 if(d.selection_model=="independent_period_logit") {
  if(d.selection_universe!="F" | d.has_rho | d.has_just)
   _cdid_pt_fail("Independent-period selection requires F without transport declarations.",459)
 }
 else {
  if(d.selection_universe=="F" & d.has_rho)
   _cdid_pt_fail("F uses unit transport; a transport column is supported only for A/B.",459)
  if(d.selection_universe!="F" & (!d.has_just | d.justification==""))
   _cdid_pt_fail("A/B selection requires an explicit transport justification.",459)
 }
}
void _cdid_pt_reference_begin(string scalar bundle,string scalar popname,string scalar planname,
 string scalar rowname,real scalar cap)
{
 pointer(struct cdid_public_bundle scalar) scalar pB
 struct cdid_catalog_opts scalar co
 struct cdid_cat scalar C
 struct cdid_baseline_plan_state scalar plan
 real matrix occ
 pointer() scalar p
 pB=findexternal(bundle)
 co=cdid_catalog_default_opts();co.mem_limit_bytes=8*cap;co.conventions="D10_D13_admissible_frame"
 occ=_cdid_cpt_base_occ(pB->prepared.population,pB->request.cohorts,co)
 C=cdid_catalog_build_opts(pB->request.calendar,pB->request.cohorts,pB->request.target_cohorts,
  pB->request.anticipation,occ,pB->request.control_group,pB->request.post,co)
 plan=cdid_baseline_plan(pB->prepared.population,C,pB->request.aux,pB->request.aux_K_ptr,pB->request.aux_K_index,cap,co.work_cap)
 p=crexternal(popname);if(p==NULL)_cdid_ref_fail("Transient population name is in use.",498)
 *p=pB->prepared.population
 p=crexternal(planname);if(p==NULL)_cdid_ref_fail("Transient plan name is in use.",498)
 *p=plan
 p=crexternal(rowname);if(p==NULL)_cdid_ref_fail("Transient source-row name is in use.",498)
 *p=pB->prepared.source_row
}
void _cdid_pt_subset_raw(string scalar srcname,string scalar dstname,string scalar raw_names)
{
 pointer(struct cdid_reference_inputs scalar) scalar ps,pd
 struct cdid_reference_inputs scalar S
 string rowvector names
 real rowvector idx
 real scalar is_same
 ps=findexternal(srcname)
 if(ps==NULL | findexternal(dstname)!=NULL)_cdid_ref_fail("Source raw staging is absent or output name already exists.",498)
 names=tokens(raw_names)
 is_same=(cols(names)==cols(ps->raw_names))
 if(is_same & cols(names)>0)is_same=all(names:==ps->raw_names)
 if(is_same) {
  S=*ps
 }
 else {
  idx=_cdid_ref_names(ps->raw_names,raw_names)
  S.population=ps->population;S.N=ps->N;S.R=ps->R;S.source_rows=ps->source_rows;S.raw_names=names
  S.data=ps->data[,idx];S.attributes=ps->attributes[,idx]
  S.has_data=ps->has_data[idx];S.has_frame=ps->has_frame[idx]
 }
 pd=crexternal(dstname);if(pd==NULL)_cdid_ref_fail("Cannot create the owned raw staging object.",498)
 *pd=S
}
void _cdid_pt_chain_begin(string scalar bundle,string scalar population,string scalar planname,
 string scalar reference,string scalar provider,string scalar process,string scalar grids,
 string scalar budgetname,string scalar rawname,string scalar catalogname,real scalar draws,real scalar alpha,real scalar bytes,
 real scalar chunk,real scalar workcap)
{
 struct cdid_pt_chain scalar T
 struct cdid_public_bundle scalar B
 struct cdid_population scalar P
 struct cdid_baseline_plan_state scalar plan
 struct cdid_baseline_provider scalar staticfit
 struct cdid_db_provider scalar datedfit
 struct cdid_selected_provider scalar selectedfit
 struct cdid_selected_profile_bundle scalar profiles
 struct cdid_selref_provider scalar selectedref
 struct cdid_selref_controls scalar sc
 struct cdid_reference_inputs scalar S,raw
 struct cdid_reference_catalog scalar catalog
 struct cdid_public_binding_control scalar ctl
 real matrix slopes
 real scalar p,k,u,i,used,max_terms,max_coords,piece_rows,s_rows,s_cells,g_chunk,e_chunk
 pointer() scalar out,pS,pfit,pCat
 pointer(struct cdid_pt_imp_owner scalar) scalar pImp
 B=*findexternal(bundle);P=*findexternal(population);plan=*findexternal(planname);raw=*findexternal(rawname)
 if(B.prepared.selection_canonical!="none" & B.prepared.selection_canonical!="conditional_mtar")
  _cdid_pt_fail("Conditional tests for this observation correction are not available in this build.",498)
 if(B.prepared.selection_est_method!="ipw" & B.prepared.selection_est_method!="imputation" & B.prepared.selection_est_method!="dr")
  _cdid_pt_fail("Conditional tests for this observation correction are not available in this build.",498)
 if(B.prepared.selection_est_method!="ipw" & B.prepared.selection_canonical!="conditional_mtar")
  _cdid_pt_fail("Imputation and augmented conditional pretests require active conditional_mtar.",459)
 if(B.prepared.selection_est_method=="dr" & B.prepared.method!="ipw" & B.prepared.method!="reg" & B.prepared.method!="dr")
  _cdid_pt_fail("Augmented conditional pretests require treatment method ipw, reg, or dr.",459)
 T.selected_active=(B.prepared.selection_canonical=="conditional_mtar" & B.prepared.selection_est_method=="ipw")
 T.cpti_active=(B.prepared.selection_canonical=="conditional_mtar" & (B.prepared.selection_est_method=="imputation" | B.prepared.selection_est_method=="dr"))
 if(B.prepared.design!="reference_recipes" & (B.reference.n_slots==0 | B.reference.n_slots>=.)) {
  _cdid_pt_ensure_static_ref(B,raw)
  *findexternal(bundle)=B
 }
 /* Only nuisance/resource controls are resolved here. Conditional
    inference owns its common sign stream and never runs final post GMM. */
 ctl=cdid_public_binding_controls(B,B.prepared.x_has_intercept,bytes,chunk,workcap,1e-10,1e-4,1,0,0,draws,alpha,"state")
 if(T.cpti_active) {
  pImp=crexternal(provider);if(pImp==NULL)_cdid_pt_fail("Conditional provider object already exists.",498)
 }
 else {
  pfit=crexternal(provider);if(pfit==NULL)_cdid_pt_fail("Conditional provider object already exists.",498)
 }
 if(B.prepared.design=="reference_recipes") {
  S=B.reference;S.population=findexternal(population);cdid_ref_activate(S)
  pS=crexternal(reference);if(pS==NULL)_cdid_pt_fail("Conditional reference object already exists.",498)
  *pS=S;catalog=cdid_reference_catalog_make(plan,S,ctl.entries)
  if(T.cpti_active) {
   T.cpti=_cdid_pt_imp_owner_prepare(pImp,B,findexternal(population),plan,pS,catalog,ctl,workcap)
   T.cpti_design_index=(*pImp).design_index
  }
  else if(T.selected_active) {
   /* Selected reference providers borrow all four owners for replay. */
   pCat=crexternal(catalogname);if(pCat==NULL)_cdid_pt_fail("Conditional catalogue object already exists.",498)
   *pCat=catalog
   sc=cdid_selref_default_controls(P.N)
   sc.method=B.prepared.method;sc.selection_model=B.prepared.selection_model
   sc.outcome_pooling=(B.prepared.outcome_pooling=="" ? "separate" : B.prepared.outcome_pooling)
   sc.selection_pooling=B.prepared.selection_pooling;sc.universe=B.prepared.selection_universe
   sc.justification=B.prepared.selection_justification;sc.has_rho=B.prepared.sel_has_rho
   sc.rho=B.prepared.sel_rho;sc.eps=ctl.baseline.glm_epsilon
   sc.maxit=ctl.baseline.glm_maxit;sc.rank_tol=ctl.baseline.rank_tol
   sc.chunk_rows=ctl.chunk;sc.max_work=(workcap<. ? workcap : 0)
   sc.max_entries=ctl.entries;sc.points_only=0
   selectedref=cdid_selref_prepare(findexternal(population),pS,findexternal(planname),pCat,sc)
   *pfit=selectedref
   T.selected=cdid_cpt_selected_prepare_dated(pfit,B.request,workcap)
  }
  else {
   datedfit=cdid_db_prepare(findexternal(population),pS,plan,catalog,ctl.baseline)
   *pfit=datedfit
   T.process=cdid_cpt_baseline_prepare_dated(pfit,B.request,workcap)
  }
 }
 else if(T.cpti_active) {
  T.cpti=_cdid_pt_imp_owner_prepare(pImp,B,findexternal(population),plan,NULL,catalog,ctl,workcap)
  T.cpti_design_index=(*pImp).design_index
 }
 else if(T.selected_active) {
  profiles=cdid_selected_static_profiles(B,plan)
  selectedfit=cdid_selected_provider_prepare(B,plan,profiles,ctl)
  /* Rebind the returned model to the immutable caller-owned population before
     retaining it; no pointer into a local returned provider survives. */
  selectedfit.model.population_ref=findexternal(population)
  *pfit=selectedfit.model
  T.selected=cdid_cpt_selected_prepare_static(pfit,B.request,workcap)
 }
 else {
  p=cols(B.prepared.X)
  slopes=(B.prepared.x_has_intercept ? (p>1 ? B.prepared.X[,2..p] : J(P.N,0,.)) : B.prepared.X)
  staticfit=cdid_baseline_prepare(plan,P,B.prepared.X,slopes,B.prepared.weight.owner_weight,ctl.baseline)
  *pfit=staticfit
  T.process=cdid_cpt_baseline_prepare_static(pfit,B.request,workcap)
 }
 if(T.cpti_active) {
  T.n_cells=T.cpti.n_cells;T.status=T.cpti.status;T.reason=T.cpti.reason
  T.workspace_entries=(*pImp).workspace_entries;T.work_ops=T.cpti.work_ops
 }
 else if(T.selected_active) {
  T.n_cells=T.selected.n_cells;T.status=T.selected.status;T.reason=T.selected.reason
  T.workspace_entries=T.selected.workspace_entries;T.work_ops=T.selected.work_ops
  if(B.prepared.design!="reference_recipes")T.workspace_entries=_cdid_base_cells(1,plan.workspace_entries,ctl.entries,T.workspace_entries)
 }
 else {
  T.n_cells=T.process.n_cells;T.status=T.process.status;T.reason=T.process.reason
  T.workspace_entries=T.process.workspace_entries;T.work_ops=T.process.work_ops
 }
 T.bundle=B;T.control=ctl;T.cap=ctl.entries;T.bytes=bytes;T.draws=draws;T.alpha=alpha
 T.input_rows=st_nobs();T.C=(B.prepared.clustered ? B.prepared.C : P.N)
 T.cmap=(B.prepared.clustered ? B.prepared.unit_cluster : (1::P.N))
 if(!B.prepared.clustered) {
  T.cluster_labels=(P.ids.string_id ? P.ids.string_id_value : strofreal(P.ids.numeric_id,"%24.17g"))
 }
 else {
  T.cluster_labels=J(T.C,1,"")
  for(u=1;u<=P.N;u++) {
   k=T.cmap[u]
   T.cluster_labels[k]=(B.prepared.cluster_kind=="string" ?
    B.prepared.cluster_string_id[u] : strofreal(B.prepared.cluster_numeric_id[u],"%24.17g"))
  }
 }
 if(T.cpti_active) {
  max_terms=4;max_coords=1
  for(i=1;i<=T.n_cells;i++)if(T.cpti.cells[i].supported) {
   if(T.cpti.cells[i].spec.R>max_terms)max_terms=T.cpti.cells[i].spec.R
   if(T.cpti.cells[i].spec.Dl>max_coords)max_coords=T.cpti.cells[i].spec.Dl
  }
  piece_rows=max((ctl.chunk,T.cpti.chunk_rows))
 }
 else {
  max_terms=32767;max_coords=32767;piece_rows=ctl.chunk
 }
 if(chunk>=. & draws*(T.C+2*P.N)*4<=T.cap) {
  s_rows=min((draws,max((ctl.chunk,512))))
  s_cells=max((T.C,s_rows*T.C))
  g_chunk=max((ctl.chunk,min((P.N,(P.N*P.N*8<=T.cap ? 2048 : 512)))))
  e_chunk=max((ctl.chunk,min((P.N,2048))))
 }
 else {
  s_rows=min((draws,ctl.chunk))
  s_cells=max((T.C,ctl.chunk*T.C))
  g_chunk=ctl.chunk
  e_chunk=ctl.chunk
 }
 T.budget=cdid_cpt_budget_create(P.N,T.C,max((1,T.n_cells)),draws,P.N,max_terms,max_coords,
  piece_rows,P.N,s_rows,s_cells,
  g_chunk,e_chunk,T.cap,bytes,workcap)
 out=crexternal(process);if(out==NULL)_cdid_pt_fail("Conditional process object already exists.",498)
 *out=T
 out=crexternal(budgetname);if(out==NULL)_cdid_pt_fail("Conditional budget object already exists.",498)
 *out=T.budget
 used=_cdid_base_cells(raw.R,cols(raw.raw_names)+2,T.cap,T.workspace_entries)
 used=_cdid_base_cells(raw.N,2*cols(raw.raw_names),T.cap,used)
 used=_cdid_base_cells(1,B.reference.retained_entries,T.cap,used)
 _cdid_pt_grids_start(grids,T.n_cells,P,T.cap,used)
 st_local("pt_cells",strofreal(T.n_cells));st_local("pt_process_available",strofreal(T.status=="available"))
}
void _cdid_pt_chain_cell(string scalar process,real scalar cell)
{
 pointer(struct cdid_pt_chain scalar) scalar pT
 pT=findexternal(process)
 if(pT->cpti_active) {
  st_local("pt_date",strofreal(pT->cpti.cells[cell].s_date,"%24.17g"))
  st_local("pt_design",strofreal(pT->cpti_design_index[cell],"%21.0f"))
 }
 else if(pT->selected_active) {
  st_local("pt_date",strofreal(pT->selected.cells[cell].s_date,"%24.17g"))
  st_local("pt_design",strofreal(pT->selected.cells[cell].design_index,"%21.0f"))
 }
 else {
  st_local("pt_date",strofreal(pT->process.cells[cell].s_date,"%24.17g"))
  st_local("pt_design",strofreal(pT->process.cells[cell].design_index,"%21.0f"))
 }
}
void _cdid_pt_chain_run(string scalar process,string scalar gridname,string scalar output,
 string scalar signs,string scalar posting,real scalar matsize,real scalar bs,real scalar cb,
 string scalar formula,string scalar clusters)
{
 pointer(struct cdid_pt_chain scalar) scalar pT
 pointer(struct cdid_pt_grid_set scalar) scalar pG
 pointer(struct cdid_population scalar) scalar pP
 struct cdid_pt_result scalar R
 struct cdid_cpt_grid colvector grids
 struct cdid_cpt_cell colvector cells
 struct cdid_cpt_mult_state scalar acc
 struct cdid_cpt_result scalar result
 struct cdid_public_sign_stream scalar stream
 struct cdid_public_sign_batch scalar batch
 struct cdid_infer_budget scalar ib
 real matrix fixed,E_all
 real scalar i,j,K,nd,used,work,rng_rows
 string scalar reason
 pT=findexternal(process);pG=findexternal(gridname);pP=&(pT->bundle.prepared.population);K=pT->n_cells
 R.mode="chained";R.method=pT->bundle.prepared.method;R.status="unavailable";R.reason=pT->reason
 R.formula=formula;R.clusters=clusters;R.rng=(signs=="" ? "caller" : "fixed")
 R.policy="all pre-treatment comparisons of the chained model before GMM combination; covariate grid of all units; unit weights from the first period"
 if(pT->selected_active)R.policy=R.policy+"; missing outcomes weighted by observation probabilities, universe "+pT->bundle.prepared.selection_universe+"; grid of the outcome covariates"
 else if(pT->cpti_active & pT->cpti.kind=="imputation")R.policy=R.policy+"; missing outcome changes imputed ("+pT->bundle.prepared.imputation_pooling+" pooling); grid of the outcome covariates"
 else if(pT->cpti_active & pT->cpti.kind=="augmented")R.policy=R.policy+"; augmented correction for missing outcomes ("+pT->bundle.prepared.selection_model+", "+pT->bundle.prepared.selection_pooling+" pooling, imputation "+pT->bundle.prepared.imputation_pooling+" pooling); grid of the outcome covariates"
 if(pT->cpti_active & pT->cpti.kind=="augmented" & pT->bundle.prepared.selection_justification!="")
  R.policy=R.policy+"; selection justification: "+pT->bundle.prepared.selection_justification
 R.channel="influence functions of the conditional pretrend moments";R.membership="distinct covariate values of all units, each counted with its unweighted multiplicity"
 R.stats=(pP->N,pP->N,K,pT->C,pT->draws,pT->alpha,.,.,.,.,pP->R,pP->R,pT->input_rows)
 R.draws=J(0,1,.);R.family=J(K,13,.);R.family_reasons=J(K,1,"")
 R.cell_ids=R.cell_status=J(K,1,"")
 R.raw_grid_variables=cdid_ref_formula_names(formula)
 R.preprocess_columns=("unit","cohort","cohort_known","observed_rows","A","B","frame_only")
 R.family_columns=("cell","group","time","from","sign","auxiliary","unavailable","design_index","reference_date","estimate","scalar_process","distinct_points","coordinates")
 R.xcols=pT->bundle.prepared.x_names;R.zcols=pT->bundle.prepared.x_names
 R.id_labels=pG->id_labels;R.cluster_positions=pT->cmap;R.cluster_labels=pT->cluster_labels
 R.cluster_values=(pT->bundle.prepared.cluster_kind=="numeric" ? pT->bundle.prepared.cluster_numeric_id : J(pP->N,1,.))
 R.retained_rows=pT->bundle.prepared.source_row
 R.rowmap=(pT->bundle.prepared.source_row,pP->ids.row_unit,pP->date,(pP->y:<.))
 R.preprocess=((1::pP->N),pP->cohort,pP->cohort_known,pP->observed_count,pP->A,pP->B,pP->frame_only)
 R.missing_grid=R.unsupported_grid=J(0,4,.)
 R.grid_first=R.multiplicity=J(0,1,.);R.grid_keys=J(K,1,"")
 R.bstrap_requested=bs;R.cband_requested=cb;R.multipliers_used=0;R.raw_CvM=.
 for(i=1;i<=K;i++) {
  if(pT->cpti_active) {
   R.family[i,]=(i,pT->cpti.cells[i].g,pT->cpti.cells[i].t,pT->cpti.cells[i].from,pT->cpti.cells[i].sign,
    pT->cpti.cells[i].q,!pT->cpti.cells[i].supported,pT->cpti_design_index[i],pT->cpti.cells[i].s_date,
    pT->cpti.cells[i].estimate_scalar,.,.,pT->cpti.cells[i].spec.Dl)
   R.family_reasons[i]=pT->cpti.cells[i].reason
   R.cell_ids[i]=pT->cpti.cells[i].id;R.cell_status[i]=pT->cpti.cells[i].status
  }
  else if(pT->selected_active) {
   R.family[i,]=(i,pT->selected.cells[i].g,pT->selected.cells[i].t,pT->selected.cells[i].from,pT->selected.cells[i].sign,
    pT->selected.cells[i].q,!pT->selected.cells[i].supported,pT->selected.cells[i].design_index,pT->selected.cells[i].s_date,
    pT->selected.cells[i].estimate,.,.,pT->selected.cells[i].spec.Dl)
   R.family_reasons[i]=pT->selected.cells[i].reason
   R.cell_ids[i]=pT->selected.cells[i].id;R.cell_status[i]=pT->selected.cells[i].status
  }
  else {
   R.family[i,]=(i,pT->process.cells[i].g,pT->process.cells[i].t,pT->process.cells[i].from,pT->process.cells[i].sign,
    pT->process.cells[i].q,!pT->process.cells[i].supported,pT->process.cells[i].design_index,pT->process.cells[i].s_date,
    pT->process.cells[i].estimate,.,.,pT->process.cells[i].spec.Dl)
   R.family_reasons[i]=pT->process.cells[i].reason
   R.cell_ids[i]=pT->process.cells[i].id;R.cell_status[i]=pT->process.cells[i].status
  }
 }
 if(pT->status!="available"){cdid_pt_put(output,R,pT->cap);return;}
 reason="";nd=0
 for(i=1;i<=rows(pG->records);i++) {
  if(pG->records[i].status!="available") {
   reason=reason+" "+pG->records[i].reason
   if(rows(R.missing_grid)==0)R.missing_grid=pG->records[i].missing
   if(rows(R.unsupported_grid)==0)R.unsupported_grid=pG->records[i].unsupported
  }
  else nd=nd+pG->records[i].grid.n_distinct
 }
 R.stats[10]=nd
 if(reason!="") {R.reason="Complete requested test unavailable: "+reason+" Estimated ATTs are unaffected; no unit is deleted.";cdid_pt_put(output,R,pT->cap);return;}
 if(pT->bundle.prepared.clustered & pT->C<2) {
  R.reason="Clustered multiplier inference needs at least two clusters; the conditional test is unavailable."
  cdid_pt_put(output,R,pT->cap);return
 }
 if(posting=="matrix" & max((pT->draws,K,nd))>matsize)
  _cdid_pt_fail("Requested conditional matrix results exceed matsize; choose posting(state) before drawing.",908)
 for(i=1;i<=K;i++) {
  j=pG->cell_record[i]
  if(j<1 | j>rows(pG->records))_cdid_pt_fail("A declared conditional cell lacks its fixed grid.",498)
  R.grid_keys[i]=pG->records[j].key
  R.family[i,12]=pG->records[j].grid.n_distinct
 }
 if(rows(pG->records)==1) {
  grids=pG->records[1].grid
 }
 else {
  grids=cdid_cpt_grid(K,1)
  for(i=1;i<=K;i++)grids[i]=pG->records[pG->cell_record[i]].grid
 }
 used=pG->retained
 used=_cdid_base_cells(pT->draws,2,pT->cap,used)
 /* complete-family builder refuses any unavailable declared cell.
    Its pieces carry residuals and original-F nuisance coordinate IFs. */
 if(pT->control.workcap<.) {
  work=pT->control.workcap-pT->work_ops-pG->work_ops
  if(work<1)_cdid_pt_fail("Fitted conditional registry exhausts the explicit counted-work limit.",908)
  pT->budget.max_work_ops=work
 }
 used=_cdid_base_cells(1,cdid_pt_cost(R),pT->cap,used)
 if(used>=pT->cap)_cdid_pt_fail("Retained conditional objects exhaust the numeric workspace before process construction.",908)
 pT->budget.max_workspace_cells=pT->cap-used;pT->budget.max_workspace_bytes=8*(pT->cap-used)
 if(pT->cpti_active)cells=cdid_cpti_build_all(pT->cpti,grids,pT->cmap,pT->C,0,pT->budget)
 else if(pT->selected_active)cells=cdid_cpt_selected_build_all(pT->selected,grids,pT->cmap,pT->C,0,pT->budget)
 else cells=cdid_cpt_baseline_build_all(pT->process,grids,pT->cmap,pT->C,0,pT->budget)
 for(j=1;j<=rows(pG->records);j++) {
  R.grid_first=R.grid_first\pG->records[j].grid.first
  R.multiplicity=R.multiplicity\pG->records[j].grid.mult
 }
 used=_cdid_base_cells(rows(R.grid_first),2,pT->cap,used)
 used=_cdid_base_cells(pT->draws,1,pT->cap,used)
 acc=cdid_cpt_mult_start(cells,pP->N,pT->C,pT->draws,pT->budget)
 fixed=(signs=="" ? J(0,0,.) : st_matrix(signs))
 rng_rows=min((pT->draws,pT->control.chunk,floor(max((pT->C,pT->control.chunk*pT->C))/pT->C)))
 ib=cdid_infer_budget_create(pP->N,pT->C,max((1,K)),pT->draws,pT->control.chunk,max((1,pT->cap)),
  rng_rows,max((pT->C,rng_rows*pT->C)),pT->cap,pT->bytes)
 stream=cdid_public_rng_start(pT->draws,pT->C,ib,R.rng,fixed)
 if(acc.batch_rows>=pT->draws & pT->draws>0 & acc.batch_rows>rng_rows) {
  E_all=J(pT->draws,pT->C,0)
  while(stream.cursor<stream.B) {
   batch=cdid_public_rng_next(stream)
   E_all[batch.draw_positions,]=batch.signs
  }
  cdid_cpt_mult_add(acc,cells,(1::pT->draws),(1..pT->C),E_all)
  E_all=J(0,0,.)
 }
 else {
  while(stream.cursor<stream.B) {
   batch=cdid_public_rng_next(stream)
   cdid_cpt_mult_add(acc,cells,batch.draw_positions,batch.cluster_positions,batch.signs)
  }
 }
 cdid_public_rng_complete(stream);R.multipliers_used=1
 result=cdid_cpt_finish(acc,cells,pT->alpha);R.raw_CvM=result.CvM
 R.status=result.status;R.reason=result.reason
 if(result.status=="available") {R.stats[7..9]=(result.CvM,result.CvMcval,result.CvMpval);R.draws=result.CvMb;}
 for(i=1;i<=K;i++){R.family[i,11]=cells[i].sigma*cells[i].mu1;R.family[i,7]=0;R.cell_status[i]="tested";}
 cdid_pt_put(output,R,pT->cap)
}
end
