*! Public baseline reference/static-FV assembly; one shared raw reader and encoder.
version 19.5
mata:
struct cdid_bs_probe {
 struct cdid_reference_inputs scalar raw
 struct cdid_ref_design scalar design
 real colvector units,history
 real matrix unit_raw,history_raw
 string rowvector names
 real rowvector variables
 real scalar ready,accepted,retained
 string scalar reason,formula
}
void _cdid_bs_start(string scalar bundle,string scalar rawname,string scalar name,
 string scalar encoded,string scalar formula,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar q
 struct cdid_reference_inputs scalar E
 struct cdid_ref_design scalar D
 struct cdid_population scalar P
 pointer() scalar p
 real scalar n,plain
 B=*findexternal(bundle);q.raw=*findexternal(rawname);q.raw.max_entries=cap;P=*q.raw.population
 q.names=cdid_ref_formula_names(formula);q.variables=_cdid_ref_names(q.raw.raw_names,invtokens(q.names))
 q.units=_cdid_ref_ix(B.prepared.required_X);q.history=J(0,1,0)
 q.formula=formula;q.ready=1;q.accepted=0;q.reason="";plain=0
 q.unit_raw=J(rows(q.units),cols(q.names),.);q.history_raw=J(0,cols(q.names),.)
 if(rows(q.units)==0){q.ready=0;q.reason="Covariates must be supplied on the fitting roster.";}
 if(cols(q.names)) {
  if(P.ids.frame_supplied) {
   if(!all(q.raw.has_frame[q.variables])) {
    q.ready=0;q.reason="Covariates must all be supplied on the unit frame for the static gate."
   }
   else q.unit_raw=q.raw.attributes[q.units,q.variables]
  }
  else {
   q.unit_raw=q.raw.data[P.first_row[q.units],q.variables]
  }
  if(all(q.raw.has_data[q.variables])) {
   q.history=_cdid_ref_ix(B.prepared.required_X[P.ids.row_unit])
   q.history_raw=q.raw.data[q.history,q.variables]
  }
  else if(any(q.raw.has_data[q.variables])) {
   q.ready=0;q.reason="Covariates are split between data and frame."
  }
  if(q.ready)if(missing(q.unit_raw)|missing(q.history_raw)) {
   q.ready=0;q.reason="Covariates contain missing or nonfinite values for units used in estimation; drop those units or impute the covariates before calling cdid."
  }
 }
 q.retained=_cdid_ref_cost(P.R,cols(q.raw.raw_names)+3,cap,0)
 q.retained=_cdid_ref_cost(P.N,2*cols(q.raw.raw_names)+3,cap,q.retained)
 q.retained=_cdid_ref_cost(rows(q.unit_raw)+rows(q.history_raw),cols(q.names)+3,cap,q.retained)
 if(q.ready) {
  n=rows(q.unit_raw)+rows(q.history_raw)
  D.keys=(1::n);D.raw_names=q.names;D.raw_index=(cols(q.names) ? (1..cols(q.names)) : J(1,0,0))
  D.formula=formula;D.intercept=B.prepared.x_has_intercept;D.encoded=0;D.dated=0
  E.designs=D;E.n_slots=1;E.max_entries=cap;E.retained_entries=q.retained
  plain=_cdid_ref_is_plain(formula,q.names)
  if(plain) {
   cdid_ref_encoding_budget(E,1,cols(q.names))
   E.retained_entries=_cdid_ref_cost(n,cols(q.names)+D.intercept,cap,E.retained_entries)
   E.designs[1].columns=(D.intercept ? ("_cons",q.names) : q.names)
   E.designs[1].expanded=q.names
   E.designs[1].omitted=J(1,cols(q.names),0)
   E.designs[1].column_map=(cols(q.names) ? (D.intercept+1..D.intercept+cols(q.names)) : J(1,0,0))
   E.designs[1].encoded=1
  }
  p=crexternal(encoded);if(p==NULL)_cdid_ref_fail("Static encoding object already exists.",498)
  *p=E
 }
 p=crexternal(name);if(p==NULL)_cdid_ref_fail("Static probe object already exists.",498)
 *p=q;st_local("bs_ready",strofreal(q.ready));st_local("bs_plain",strofreal(plain))
}
void _cdid_bs_stage(string scalar name)
{
 struct cdid_bs_probe scalar q
 real matrix raw
 real scalar j
 q=*findexternal(name);raw=q.unit_raw\q.history_raw
 if(st_nobs()!=0|st_nvar()!=0)_cdid_ref_fail("Static FV probe needs its empty owned frame.",498)
 st_addobs(rows(raw))
 for(j=1;j<=cols(raw);j++){(void)st_addvar("double",q.names[j]);st_store(.,q.names[j],raw[,j]);}
}
void _cdid_bs_finish(string scalar bundle,string scalar name,string scalar encoded,
 string scalar time,real scalar encode_rc)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar q
 struct cdid_reference_inputs scalar E,C
 struct cdid_population scalar P
 real colvector hit,u_pos
 real matrix history_X,unit_X
 B=*findexternal(bundle);q=*findexternal(name);P=*q.raw.population
 if(encode_rc!=0 & encode_rc!=198 & encode_rc!=459)exit(encode_rc)
 if(q.ready & encode_rc==0) {
  E=*findexternal(encoded)
  if(rows(E.designs[1].value)==0 & E.designs[1].encoded) {
   if(cols(E.designs[1].columns)==0) {
    q.reason="Covariates must be supplied on the fitting roster."
   }
   else {
    q.accepted=1
    if(rows(q.history) & cols(q.unit_raw)>0) {
     if(rows(q.units)==P.N)hit=P.ids.row_unit[q.history]
     else {u_pos=J(P.N,1,0);u_pos[q.units]=(1::rows(q.units));hit=u_pos[P.ids.row_unit[q.history]];}
     if(!all(q.history_raw:==q.unit_raw[hit,])) {
      q.accepted=0;q.reason="Covariates vary within unit across time or disagree with the unit frame."
     }
    }
    if(q.accepted) {
     unit_X=(E.designs[1].intercept ? (J(rows(q.units),1,1),q.unit_raw) : q.unit_raw)
     q.design=E.designs[1];q.design.keys=q.units;q.design.value=unit_X
     q.design.raw_index=J(1,0,0);q.design.raw_names=J(1,0,"")
     q.retained=E.retained_entries
    }
   }
  }
  else {
   unit_X=E.designs[1].value[(1::rows(q.units)),]
   if(cols(E.designs[1].value)==0) {
    q.reason="Covariates must be supplied on the fitting roster."
   }
   else {
    q.accepted=1
    if(rows(q.history)) {
     history_X=E.designs[1].value[rows(q.units)+1..rows(E.designs[1].value),]
     if(rows(q.units)==P.N)hit=P.ids.row_unit[q.history]
     else {u_pos=J(P.N,1,0);u_pos[q.units]=(1::rows(q.units));hit=u_pos[P.ids.row_unit[q.history]];}
     if(!all(history_X:==unit_X[hit,])) {
      q.accepted=0;q.reason="Covariates vary within unit across time or disagree with the unit frame."
     }
    }
    if(q.accepted) {
     q.design=E.designs[1];q.design.keys=q.units;q.design.value=unit_X
     q.design.raw_index=J(1,0,0);q.design.raw_names=J(1,0,"")
     q.retained=E.retained_entries
    }
   }
  }
 }
 else if(encode_rc!=0)q.reason="The evaluated static FV design is missing, nonfinite, or invalid (native code "+strofreal(encode_rc)+")."
 B.baseline_static_reason=q.reason
 B.baseline_static_status=(q.accepted ? "passed" : "failed_dated_fallback")
 *findexternal(bundle)=B;*findexternal(name)=q
 if(!q.accepted) {
  C=cdid_ref_source(q.raw.population,q.raw.source_rows,q.names,q.raw.data[,q.variables],
   q.raw.attributes[,q.variables],q.raw.has_data[q.variables],q.raw.has_frame[q.variables],time,q.raw.max_entries)
  if(!any((C.kind:==2):|(C.kind:==3)))_cdid_pr_fail("cdid_chained_request_covariate_error",q.reason)
 }
 st_local("bs_accepted",strofreal(q.accepted))
}
void _cdid_bs_reference(string scalar bundle,string scalar name,string scalar planname,
 string scalar output,string scalar weight,string scalar fixweights,string scalar time,
 real scalar anticipation,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar q
 struct cdid_reference_inputs scalar S
 struct cdid_baseline_plan_state scalar plan
 real rowvector wi,hasdata,hasframe
 real matrix data,frame
 string rowvector names
 real scalar slot
 pointer() scalar p
 B=*findexternal(bundle);q=*findexternal(name);plan=*findexternal(planname)
 wi=(weight=="" ? J(1,0,0) : _cdid_ref_names(q.raw.raw_names,weight))
 names=(weight=="" ? J(1,0,"") : weight)
 data=q.raw.data[,wi];frame=q.raw.attributes[,wi]
 hasdata=q.raw.has_data[wi];hasframe=q.raw.has_frame[wi]
 S=cdid_ref_source(q.raw.population,q.raw.source_rows,names,data,frame,hasdata,hasframe,time,cap)
 if(cols(wi))if(hasdata[1]) {
  if(any(select(data[,1],data[,1]:<.):<=0))_cdid_ref_fail("Sampling weights must be finite and positive where supplied.",459)
 }
 cdid_ref_require(S,plan,anticipation,"baseline",B.prepared.method,J(S.N,1,1),fixweights,
  q.formula,(B.prepared.method=="reg" ? "" : q.formula),"","","","",
  B.prepared.x_has_intercept,(B.prepared.method=="reg" ? 0 : B.prepared.x_has_intercept),0,(cols(wi) ? 1 : 0))
 cdid_ref_designs(S)
 for(slot=1;slot<=S.n_slots;slot++) {
  S.designs[slot]=q.design
  S.retained_entries=_cdid_ref_cost(rows(q.design.value),cols(q.design.value),cap,S.retained_entries)
 }
 S.retained_entries=_cdid_ref_cost(1,q.retained,cap,S.retained_entries)
 p=crexternal(output);if(p==NULL)_cdid_ref_fail("Reference output object already exists.",498)
 *p=S
}
/* Build actual population/settings before reference classification, without a
   temporary fitted design or fabricated weights. The shared reader attaches
   the only scientific X and weight values after its real support validation. */
struct cdid_public_prepared scalar cdid_baseline_input_prepare(struct cdid_ci_input scalar I,
 string scalar method,real scalar noconstant,string scalar fixweights,string scalar covariance,
 | string scalar outcome_pooling_in)
{
 struct cdid_public_prepared scalar B
 struct cdid_population_ids scalar ids
 struct cdid_sel_decl scalar decl
 string scalar gcov,outpool
 real scalar out_active
 if(method!="ipw" & method!="reg" & method!="dr")_cdid_pr_fail("cdid_chained_request_input_error","method() must be dr, reg, ipw, or custom.")
 if(noconstant!=0 & noconstant!=1)_cdid_pr_fail("cdid_chained_request_input_error","noconstant must be binary.")
 out_active=(method=="reg" | method=="dr")
 if(out_active) {
  outpool=(args()>=6 ? strlower(strtrim(outcome_pooling_in)) : "")
  if(outpool=="")outpool="separate"
  if(outpool!="separate" & outpool!="pooled")_cdid_prep_role_fail("cdid_chained_request_outcome_pooling","outcomepooling() must be separate or pooled.",198)
 }
 else outpool=""
 gcov=strlower(strtrim(covariance));if(gcov=="")gcov="tikhonov"
 if(gcov!="tikhonov")_cdid_prep_role_fail("cdid_chained_request_input_error","gmmcovariance() must be tikhonov.",198)
 decl=cdid_sel_norm_decl("none",0,"F",I.frame_has_rho|I.data_has_rho,0,"pair_logit","separate","")
 if(I.string_id)ids=cdid_population_ids_string(I.string_ids,I.frame_string_ids,I.frame_supplied)
 else ids=cdid_population_ids_numeric(I.id,I.frame_id,I.frame_supplied)
 B.population=cdid_population_prepare(ids,I.date,I.y,I.g,I.frame_g,I.frame_has_g)
 if(method=="ipw" & any(B.population.cohort_known:==0))
  _cdid_pr_fail("cdid_chained_request_cohort_error","method(ipw) requires a known gvar() value for every unit, including population-frame units.")
 B.method=method;B.design="unencoded";B.source_row=I.source_row
 B.columns=J(0,1,cdid_population_column());B.X=J(0,0,.);B.x_names=J(1,0,"")
 B.x_has_intercept=(noconstant==0)
 B.required_X=((method=="ipw" | cdid_etwfe_requested()) ? J(B.population.N,1,1) : B.population.A)
 B.selection=decl.selection;B.selection_canonical=decl.canonical_assumption
 B.selection_model=decl.selection_model;B.selection_pooling=decl.selection_pooling
 B.selection_universe=decl.selection_universe;B.selection_est_method="ipw"
 B.selection_justification=decl.justification;B.imputation_pooling="";B.outcome_pooling=outpool
 B.inactive_arguments=(out_active ? "imputation_pooling" : "outcome_pooling imputation_pooling");B.gmm_covariance=gcov
 B.sel_columns=J(0,1,cdid_population_column());B.Z_sel=J(B.population.N,0,.)
 B.sel_names=J(1,0,"");B.U=J(B.population.N,1,1);B.sel_rho=J(B.population.N,1,1)
 B.sel_has_rho=0;B.sel_has_intercept=B.x_has_intercept;B.selection_formula_inherited=1
 B.weight.fix_weights=fixweights
 B.clustered=cdid_ci_clusters(B.population,I,B.cluster_kind,B.cluster_numeric_id,
  B.cluster_string_id,B.unit_cluster,B.C)
 return(B)
}
void _cdid_baseline_prepare_store(string scalar input_name,string scalar output_name,
    string scalar method,real scalar noconstant,string scalar fixweights,
    real colvector targets,string scalar control,real scalar ant,string scalar base,
    string scalar pt,real scalar ptstart,real scalar max_entries,
    | string scalar selection_in,string scalar selection_model_in,
      string scalar selection_pooling_in,string scalar selection_universe_in,
      string scalar selection_justification_in,string scalar selection_est_method_in,
      string scalar imputation_pooling_in,string scalar gmm_covariance_in,
      real scalar selnoconstant_in,real scalar selconstant_in,real scalar selzero_in,
      string scalar outcome_pooling_in)
{
    struct cdid_ci_input scalar I
    struct cdid_public_bundle scalar B
    pointer() scalar p
    p=findexternal(input_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Transient input is absent.")
    if(findexternal(output_name)!=NULL) _cdid_pr_fail("cdid_chained_request_input_error","Requested output name already exists; no user object is overwritten.")
    I=*p
    B.prepared=cdid_baseline_input_prepare(I,method,noconstant,fixweights,
        (args()>=20 ? gmm_covariance_in : "tikhonov"),
        (args()>=24 ? outcome_pooling_in : ""))
    B.request=cdid_public_request_build(B.prepared.population,targets,control,ant,base,pt,ptstart,max_entries)
    B.settings=cdid_public_defaults();B.settings.method=method;B.settings.control_group=control
    B.settings.anticipation=ant;B.settings.base_period=base;B.settings.pt_assumption=pt;B.settings.pt_start=ptstart
    B.settings.selection=B.prepared.selection
    B.settings.selection_canonical=B.prepared.selection_canonical
    B.settings.selection_model=B.prepared.selection_model
    B.settings.selection_pooling=B.prepared.selection_pooling
    B.settings.selection_universe=B.prepared.selection_universe
    B.settings.selection_est_method=B.prepared.selection_est_method
    B.settings.selection_justification=B.prepared.selection_justification
    B.settings.imputation_pooling=B.prepared.imputation_pooling
    B.settings.outcome_pooling=B.prepared.outcome_pooling
    B.settings.inactive_arguments=B.prepared.inactive_arguments
    B.settings.gmm_covariance=B.prepared.gmm_covariance
    p=crexternal(output_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Cannot create requested output object.")
    *p=B
    st_local("pr_N",strofreal(B.prepared.population.N));st_local("pr_A",strofreal(B.prepared.population.n_A))
    st_local("pr_post",strofreal(rows(B.request.post)));st_local("pr_aux",strofreal(rows(B.request.aux)))
    st_local("pr_report",strofreal(rows(B.request.report)))
}
void _cdid_bs_pop_export(string scalar bundle,string scalar popname,string scalar rowname)
{
 struct cdid_public_bundle scalar B
 pointer() scalar p
 B=*findexternal(bundle)
 p=crexternal(popname);if(p==NULL)_cdid_ref_fail("Transient population name is in use.",498)
 *p=B.prepared.population
 p=crexternal(rowname);if(p==NULL)_cdid_ref_fail("Transient source-row name is in use.",498)
 *p=B.prepared.source_row
}
void _cdid_bs_begin_plan(string scalar bundle,string scalar planname,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_catalog_opts scalar co
 struct cdid_pairs scalar pairs
 struct cdid_cat scalar C
 struct cdid_baseline_plan_state scalar plan
 pointer() scalar p
 B=*findexternal(bundle)
 co=cdid_catalog_default_opts();co.mem_limit_bytes=8*cap;co.conventions="D10_D13_admissible_frame"
 pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
 C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
  B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
 plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,cap,co.work_cap)
 p=crexternal(planname);if(p==NULL)_cdid_ref_fail("Transient plan name is in use.",498)
 *p=plan
}
void _cdid_bs_try_static_attach(string scalar bundle,string scalar name,
 string scalar weight,string scalar fixweights,string scalar time,real scalar cap)
{
 struct cdid_public_bundle scalar B
 struct cdid_bs_probe scalar q
 struct cdid_population scalar P
 struct cdid_reference_inputs scalar S
 real rowvector wi
 real colvector wd,wf,u_act,q_w,r1
 real scalar hd,hf,ok_w,used
 st_local("bs_static_done","0")
 B=*findexternal(bundle);q=*findexternal(name);P=*q.raw.population
 if(!q.accepted) return
 if(fixweights!="" & fixweights!="first_period" & fixweights!="base_period") return
 ok_w=0;hd=0;hf=0;wd=J(0,1,.);wf=J(0,1,.)
 if(weight=="") ok_w=1
 else if(weight!=time) {
  wi=_cdid_ref_names(q.raw.raw_names,weight)
  if(cols(wi)==1) {
   hd=q.raw.has_data[wi[1]];hf=q.raw.has_frame[wi[1]]
   if((hd | hf) & (hf | P.n_frame_only==0)) {
    wd=(hd ? q.raw.data[,wi[1]] : J(0,1,.))
    wf=(hf ? q.raw.attributes[,wi[1]] : J(0,1,.))
    if((!hd | (!missing(wd) & all(wd:>0))) & (!hf | (!missing(wf) & all(wf:>0)))) {
     if(hd & hf) {
      if(all(wd:==wf[P.ids.row_unit])) ok_w=1
     }
     else if(!hd & hf) ok_w=1
     else {
      u_act=_cdid_ref_ix(!P.frame_only);q_w=J(P.N,1,1)
      if(rows(u_act)) q_w[u_act]=wd[P.first_row[u_act]]
      if(all(wd:==q_w[P.ids.row_unit])) ok_w=1
      else if(fixweights=="first_period") {
       r1=cdid_population_rows(P,(1::P.N),J(P.N,1,P.calendar[1]))
       if(all(r1:>0)) ok_w=1
      }
     }
    }
   }
  }
 }
 if(!ok_w) return
 B.prepared.weight=cdid_population_weights(P,wd,wf,hd,hf,fixweights,0)
 if(B.prepared.weight.reduction!="static" & B.prepared.weight.reduction!="first_period_static") return
 if(weight!="") B.prepared.weight.source=(hf ? "frame" : "data")
 B.prepared.design="static_numeric"
 B.prepared.x_names=q.design.columns
 B.prepared.X=J(P.N,cols(q.design.value),.)
 B.prepared.X[q.design.keys,]=q.design.value
 B.prepared.required_X=((B.prepared.method=="ipw" | cdid_etwfe_requested()) ? J(P.N,1,1) : P.A)
 S.original_ids.string_id=P.ids.string_id;S.original_ids.frame_supplied=P.ids.frame_supplied
 S.original_ids.numeric_id=P.ids.numeric_id;S.original_ids.string_id_value=P.ids.string_id_value
 S.original_ids.row_unit=J(0,1,0);S.calendar=P.calendar
 S.N=P.N;S.R=P.R;S.max_entries=cap
 S.roles=cdid_ref_role(3,1)
 S.roles[1].active=1;S.roles[1].intercept=q.design.intercept;S.roles[1].dated=0;S.roles[1].slot=1
 S.roles[1].formula=invtokens(tokens(q.formula));S.roles[1].keys=q.design.keys
 S.roles[1].first_view=J(rows(q.design.keys),1,1);S.roles[1].fit_keys=q.design.keys
 if(B.prepared.method!="reg") {
  S.roles[2]=S.roles[1];S.roles[2].slot=2
  S.designs=(q.design\q.design);S.n_slots=2
 }
 else {
  S.designs=q.design;S.n_slots=1
 }
 used=_cdid_ref_cost(rows(q.design.value),cols(q.design.value)*S.n_slots,cap,0)
 S.retained_entries=_cdid_ref_cost(1,q.retained,cap,used)
 S.weight_kind=(B.prepared.weight.kind=="varying" ? "dated" : "static")
 S.fix_weights=fixweights;S.V=B.prepared.weight.V;S.population=NULL
 B.reference=S
 *findexternal(bundle)=B
 st_local("reference_active","0")
 st_local("bs_static_done","1")
}
void _cdid_baseline_reference_begin(string scalar bundle,string scalar popname,string scalar planname,
 string scalar rowname,real scalar cap)
{
 _cdid_bs_pop_export(bundle,popname,rowname)
 _cdid_bs_begin_plan(bundle,planname,cap)
}
void _cdid_baseline_reference_attach(string scalar bundle,string scalar refname,string scalar wname)
{
 struct cdid_public_bundle scalar B
 struct cdid_reference_inputs scalar S
 real scalar active,j,slot,first_static
 real colvector keys
 B=*findexternal(bundle);S=*findexternal(refname)
 /* FV spelling alone never activates dated. Its static A/F basis is fixed by
    the common reader, then passed to the unchanged static scientific engine. */
 /* Source first_period_static is the existing static engine with w_t1.
    Complete first-calendar W on original F is required, not just A. */
 first_static=(B.baseline_static_status=="passed" & S.weight_kind=="dated" &
  S.fix_weights=="first_period" & !missing(S.V))
 active=(B.baseline_static_status!="passed" | (S.weight_kind=="dated" & !first_static));slot=S.roles[1].slot
 st_local("reference_active",strofreal(active))
 B.prepared.x_names=S.designs[slot].columns
 if(active) {
  B.prepared.design="reference_recipes"
  B.prepared.X=J(0,cols(S.designs[slot].value),.)
 }
 else {
  B.prepared.design="static_numeric"
  B.prepared.X=J(S.N,cols(S.designs[slot].value),.)
  keys=S.designs[slot].keys
  B.prepared.X[keys,]=S.designs[slot].value
 }
 B.prepared.required_X=((B.prepared.method=="ipw" | cdid_etwfe_requested()) ? J(S.N,1,1) : B.prepared.population.A)
 if(S.weight_kind=="dated" & !first_static) {
  B.prepared.weight.kind="varying";B.prepared.weight.fix_weights=S.fix_weights
  B.prepared.weight.reduction="reference_recipes"
  B.prepared.weight.V=S.V;B.prepared.weight.owner_weight=S.V
  B.prepared.weight.scientific_static=J(0,1,.)
 }
 else {
  if(wname=="")B.prepared.weight=cdid_population_weights(B.prepared.population,J(0,1,.),J(0,1,.),0,0,S.fix_weights,0)
  else {
   j=selectindex(S.raw_names:==wname)[1]
   B.prepared.weight=cdid_population_weights(B.prepared.population,
    (S.has_data[j] ? S.data[,j] : J(0,1,.)),(S.has_frame[j] ? S.attributes[,j] : J(0,1,.)),
    S.has_data[j],S.has_frame[j],S.fix_weights,0)
  }
  if(B.prepared.weight.reduction!="static" & B.prepared.weight.reduction!="first_period_static")
   _cdid_pr_fail("cdid_reference_support_error","The static route requires the accepted static or complete first-period weight reduction.")
  if(first_static & B.prepared.weight.reduction!="first_period_static")
   _cdid_pr_fail("cdid_reference_support_error","Complete first-period W must retain its accepted static weight reduction.")
 }
 if(active) {
  B.prepared.weight.reduction="reference_recipes"
  B.prepared.weight.V=S.V
 }
 if(wname!="") {
  j=selectindex(S.raw_names:==wname)[1]
  B.prepared.weight.source=(S.has_frame[j] ? "frame" : "data")
 }
 S.population=NULL;B.reference=S
 *findexternal(bundle)=B
}

real matrix cdid_baseline_reference_core(pointer() scalar context,real scalar lo,real scalar hi)
{
 struct cdid_db_provider scalar f
 real colvector ids,map
 real matrix values
 real scalar j,k
 f=*context;ids=uniqrows(sort(f.catalog.post_consumer,1));map=J(f.catalog.M,1,0)
 for(j=1;j<=f.catalog.M;j++)for(k=1;k<=rows(ids);k++)if(f.catalog.post_consumer[j]==ids[k])map[j]=k
 values=cdid_db_weight_chunk(f,lo,hi,ids)
 return(values[,map])
}
void cdid_baseline_reference_metadata(struct cdid_result_state scalar out,
 struct cdid_db_provider scalar f,struct cdid_baseline_plan_state scalar plan,
 struct cdid_reference_inputs scalar S)
{
 real scalar i,c,r
 string scalar key,id
 real colvector K
 out.numeric_settings=out.numeric_settings\S.n_recipes\f.plan.n_blocks\f.plan.n_consumers\f.work_ops\f.workspace_peak_entries
 out.numeric_setting_names=out.numeric_setting_names\"n_reference_recipes"\"n_reference_OR_blocks"\"n_unique_consumers"\"provider_counted_work"\"provider_peak_numeric_entries"
 out.text_settings=out.text_settings\("profile_kind","reference_recipes")\("reference_route","baseline")\
  ("reference_weight_kind",S.weight_kind)\("static_S3","inactive: actual D60 baseline reference route")
 for(i=1;i<=rows(plan.raw);i++) {
  key=_cdid_rb_key(i);c=f.catalog.raw_consumer[i];r=S.view_recipe[i]
  id=(i<=plan.M ? "cp04:catalog:"+_cdid_rb_key(i) : "cp04:auxiliary:"+_cdid_rb_key(i-plan.M))
  out.numeric_settings=out.numeric_settings\(c>0 ? f.points[c] : .)\plan.cohorts[plan.raw[i,1]]\
   plan.calendar[plan.raw[i,2]]\plan.calendar[plan.raw[i,3]]\r\plan.raw_supported[i]\c
  out.numeric_setting_names=out.numeric_setting_names\("source_link_point:"+key)\("source_link_g:"+key)\
   ("source_link_s:"+key)\("source_link_t:"+key)\("source_link_recipe:"+key)\("source_link_supported:"+key)\("source_link_consumer:"+key)
  out.text_settings=out.text_settings\("source_link_id:"+key,id)
 }
 out.source_versions=out.source_versions\("estimator","reference-date baseline estimator; common outcome regression; influence functions over all units")
}
struct cdid_result_state scalar cdid_baseline_reference_report(
    struct cdid_public_bundle scalar B,struct cdid_db_provider scalar fit,
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
    Q=rows(B.request.aux);M=rows(fit.catalog.post_consumer);Jt=rows(post_est_ids)
    ret_if=(args()>=13 ? retained_if : J(0,0,.))
    if(rows(post_loadings)!=Jt | cols(post_loadings)!=M | rows(fit.catalog.aux_consumer)!=Q |
        chunk<1 | chunk>=. | chunk!=floor(chunk) | cap<1 | cap>=. | cap!=floor(cap))
        _cdid_pr_fail("cdid_native_reporting_input","Provider/post loading identities or replay controls do not align.")
    if(rows(post_est_ids)) {
        _cdid_rb_ids(post_est_ids,rows(B.request.post))
        if(any(post_available[post_est_ids]:!=1) | rows(post_est_ids)!=sum(post_available))
            _cdid_pr_fail("cdid_native_reporting_input","Estimable post IDs must equal the full supported post set.")
    }
    else if(sum(post_available)) _cdid_pr_fail("cdid_native_reporting_input","Supported post set has no passed target loading.")
    aux_point=J(Q,1,.);aux_available=fit.catalog.aux_supported
    active=_cdid_rb_index(aux_available)
    if(rows(active)) aux_point[active]=fit.points[fit.catalog.aux_consumer[active]]
    buffer=cdid_public_report_start(B,post_point,post_available,post_status,post_reason,
        aux_point,aux_available,fit.catalog.aux_reason,compute_if,limits)
    if(compute_if) {
        qids=J(0,1,.)
        for(i=1;i<=rows(B.request.report);i++) if(B.request.report[i,5]==2 & buffer.payload.available[i])
            qids=qids\B.request.report[i,7]
        if(rows(qids)) qids=uniqrows(sort(qids,1))
        consumers=(Jt ? fit.catalog.post_consumer : J(0,1,.))
        if(rows(qids)) consumers=consumers\fit.catalog.aux_consumer[qids]
        if(rows(consumers)) consumers=uniqrows(sort(consumers,1))
        postpos=J((Jt ? M : 0),1,.)
        for(j=1;j<=rows(postpos);j++) {
            hit=_cdid_rb_index(consumers:==fit.catalog.post_consumer[j]);postpos[j]=hit[1]
        }
        auxpos=J(rows(qids),1,.)
        for(j=1;j<=rows(qids);j++) {
            hit=_cdid_rb_index(consumers:==fit.catalog.aux_consumer[qids[j]]);auxpos[j]=hit[1]
        }
        // Retained result + provider + target loading, plus conservative bounded
        // provider/scatter scratch. Never construct originalF-by-allmoments IF.
        retained=cdid_state_cost(buffer.payload)[1]+fit.workspace_entries+length(post_loadings)
        perrow=10*(rows(consumers)+2*cols(B.prepared.X)+4)+3*(Jt+rows(qids)+rows(B.request.report))
        if(retained>=cap | perrow>cap-retained) _cdid_rs_fail("Reported IF/provider replay exceeds explicit numeric workspace.",908)
        effective=min((chunk,fit.chunk_rows,floor((cap-retained)/perrow)))
        use_ret=(rows(ret_if)==fit.n_F & cols(ret_if)==M & fit.weight_row_domain=="F")
        if(use_ret) {
            aux_in_post=J(rows(qids),1,0);rem_consumers=J(0,1,.);rem_pos=J(rows(qids),1,0)
            for(j=1;j<=rows(qids);j++) {
                hit=_cdid_rb_index(fit.catalog.post_consumer:==fit.catalog.aux_consumer[qids[j]])
                if(rows(hit)) aux_in_post[j]=hit[1]
                else rem_consumers=rem_consumers\fit.catalog.aux_consumer[qids[j]]
            }
            if(rows(rem_consumers)) {
                rem_consumers=uniqrows(sort(rem_consumers,1))
                for(j=1;j<=rows(qids);j++) if(!aux_in_post[j]) {
                    hit=_cdid_rb_index(rem_consumers:==fit.catalog.aux_consumer[qids[j]])
                    rem_pos[j]=hit[1]
                }
            }
            for(lo=1;lo<=fit.n_F;lo=lo+effective) {
                hi=min((fit.n_F,lo+effective-1));uid=(lo::hi)
                postIF=(Jt ? (fit.core_to_frame:*ret_if[uid,])*post_loadings' : J(hi-lo+1,0,.))
                if(rows(qids)) {
                    auxIF=J(hi-lo+1,rows(qids),.)
                    values=(rows(rem_consumers) ? cdid_db_if_rows(fit,uid,rem_consumers,"frame") : J(hi-lo+1,0,.))
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
                values=(rows(consumers) ? cdid_db_if_rows(fit,uid,consumers,"frame") : J(hi-lo+1,0,.))
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
struct cdid_public_fit_result scalar cdid_baseline_reference_fit(struct cdid_public_bundle scalar B,
    struct cdid_public_binding_control scalar ctrl)
{
    struct cdid_public_fit_result scalar out
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_catalog_opts scalar co
    struct cdid_baseline_plan_state scalar plan
    struct cdid_db_provider scalar fit
    struct cdid_baseline_system scalar sys
    struct cdid_public_post_support scalar support
    struct cdid_population scalar owner
    struct cdid_reference_inputs scalar reference
    struct cdid_reference_catalog scalar catalog
    struct cdid_gmm_result scalar sol
    struct cdid_baseline_controls scalar base_ctrl
    real matrix Z,G,loadings,postmeta,Vw,Lx
    real colvector point,diagnostic,reportpost,kept,rk
    string colvector diagnostic_names
    real scalar M,Jt,used,remaining,i,j,p,keepif,p_ctrl,def_pu,def_ch,eff_pu,eff_chunk,max_rem
    if(B.prepared.design!="reference_recipes" | B.prepared.selection_canonical!="none")
        _cdid_pr_fail("cdid_native_private_binding_scope","Reference baseline binding requires the active baseline reference route.")
    if(B.settings.gmm_weighting!="optimal" & B.settings.gmm_weighting!="two_step" & B.settings.gmm_weighting!="identity" & B.settings.gmm_weighting!="within")
        _cdid_pr_fail("cdid_chained_public_control_error","gmmweighting() must be optimal, two_step, identity, or within.")
    B.settings.alpha=ctrl.inference.alpha;B.settings.compute_if=ctrl.inference.compute_inffunc
    B.settings.bstrap=ctrl.inference.bstrap;B.settings.cband=ctrl.inference.cband;B.settings.biters=ctrl.inference.B
    co=cdid_catalog_default_opts();co.mem_limit_bytes=ctrl.workspace_bytes;co.work_cap=ctrl.workcap
    co.conventions="D10_D13_admissible_frame"
    pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
    C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
        B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
    D=cdid_chained_ipw_build_design(C,B.request.post,B.request.pt_assumption,B.request.pt_start)
    B.request=cdid_public_request_tags(B.request,C)
    // ALL original post rows and BOTH reporting families precede ONE shared OR.
    plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,B.request.aux_K_ptr,B.request.aux_K_index,ctrl.entries,ctrl.workcap)
    p_ctrl=max((1,cols(B.prepared.X),cols(B.prepared.Z_sel)))
    def_pu=_cdid_gmm_cells(32,ctrl.links_bound+p_ctrl+1,2^53-1,0)
    def_ch=max((1,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,def_pu)))))))
    eff_chunk=ctrl.chunk
    if(ctrl.chunk==def_ch & ctrl.baseline.chunk_rows==ctrl.chunk) {
        eff_pu=_cdid_gmm_cells(32,rows(C.r_tp)+p_ctrl+1,2^53-1,0)
        eff_chunk=max((ctrl.chunk,min((B.prepared.population.N,16384,floor(ctrl.entries/max((1,eff_pu)))))))
    }
    base_ctrl=ctrl.baseline;base_ctrl.chunk_rows=eff_chunk
    owner=B.prepared.population;reference=B.reference;reference.population=&owner
    cdid_ref_activate(reference)
    catalog=cdid_reference_catalog_make(plan,reference,ctrl.entries)
    fit=cdid_db_prepare(&owner,&reference,plan,catalog,base_ctrl)
    sys=cdid_baseline_make_system(C,D,B.request.post,ctrl.entries,fit.workspace_entries)
    support=cdid_public_post_ledger(C,D,B.request,sys)
    M=sys.M;Jt=rows(B.request.post);point=J(Jt,1,.);loadings=J(0,M,.)
    diagnostic=J(0,1,.);diagnostic_names=J(0,1,"")
    if(M>0 & sys.p_retained>0) {
        if(eff_chunk>ctrl.chunk) {
            max_rem=floor((ctrl.entries-fit.workspace_entries-sys.workspace_entries-3*M*M)/max((1,10*(M+2*cols(B.prepared.X)+4))))
            eff_chunk=(max_rem>=ctrl.chunk ? min((eff_chunk,max_rem)) : ctrl.chunk)
        }
        used=fit.workspace_entries+sys.workspace_entries
        used=_cdid_base_cells(8*eff_chunk,M+2*cols(B.prepared.X)+4,ctrl.entries,used)
        if(B.settings.gmm_weighting=="identity") {
            used=_cdid_gmm_identity_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries
            sol=cdid_gmm_identity_map(fit.points[catalog.post_consumer],sys.full_A,sys.H_est,remaining)
        }
        else {
            used=_cdid_base_cells(2*eff_chunk,M,ctrl.entries,used)
            used=_cdid_base_cells(3*M,M,ctrl.entries,used)
            used=_cdid_gmm_workspace(M,sys.p_retained,sys.n_estimable,ctrl.entries,used)
            remaining=ctrl.entries-fit.workspace_entries-sys.workspace_entries-M*M
            if(B.settings.gmm_weighting=="within") {
                rk=fit.points[catalog.post_consumer]
                cdid_wsg_within(B,plan,C,rk,0,ctrl.baseline.glm_epsilon,ctrl.baseline.glm_maxit,eff_chunk,kept,G,Vw,reference)
                sol=cdid_gmm_fit_within_norm(rk[kept],sys.full_A[kept,],Vw,G,M,sys.H_est,
                    ctrl.condition_trigger,ctrl.ridge_relative,fit.weight_rows,fit.normalization_count,
                    eff_chunk,&cdid_baseline_reference_core(),&fit,remaining,M,kept')
                Lx=J(rows(sol.target_loadings),M,0);Lx[.,kept']=sol.target_loadings;sol.target_loadings=Lx
                Lx=Vw=J(0,0,.)
            }
            else {
                kept=(1::M)
                G=I(M)
                sol=cdid_gmm_fit_unit_norm(fit.points[catalog.post_consumer],sys.full_A,G,M,sys.H_est,
                    B.settings.gmm_weighting,ctrl.condition_trigger,ctrl.ridge_relative,
                    fit.weight_rows,fit.normalization_count,eff_chunk,&cdid_baseline_reference_core(),&fit,remaining)
            }
        }
        if(B.settings.gmm_weighting=="identity") kept=(1::M)
        if(sys.n_estimable) point[sys.estimable_ids]=sol.targets
        loadings=sol.target_loadings
        cdid_chain_sp_record(B,C,D,kept,fit.points[catalog.post_consumer],loadings,sys.estimable_ids,
            cdid_chain_auxpts(fit.points,fit.catalog.aux_consumer,fit.catalog.aux_supported),fit.catalog.aux_supported)
        diagnostic=(sol.criterion\sol.lambda\sol.rho\sol.trace_B\sol.covariance_rank\sol.source_n_units\sol.normalization_factor\sol.sample_rank_upper)
        diagnostic_names=("GMM_criterion"\"GMM_lambda"\"GMM_ridge_rho"\"GMM_trace_B"\"GMM_covariance_rank"\"GMM_source_replay_rows"\"GMM_gram_normalization_factor"\"GMM_sample_rank_upper")
        // Only final reported IF is retained for durable no-refit results.
        G=J(0,0,.);sol.unit_gram_covariance=sol.coefficient_loadings=sol.target_loadings=J(0,0,.)
        sol.coefficient_covariance=sol.covariance=J(0,0,.)
    }
    else if(sys.n_estimable) _cdid_pr_fail("cdid_native_reporting_input","An estimable post target has no GMM parameter system.")
    sys.full_A=sys.H_est=J(0,0,.)
    out.payload=cdid_baseline_reference_report(B,fit,point,sys.post_available,support.status,support.reason,
        sys.estimable_ids,loadings,ctrl.inference.compute_inffunc,ctrl.inference.limits,eff_chunk,ctrl.entries,sol.retained_if)
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
    cdid_baseline_reference_metadata(out.payload,fit,plan,reference)
    out.inference=ctrl.inference
    cdid_state_validate(out.payload,ctrl.inference.limits)
    return(out)
}

end
