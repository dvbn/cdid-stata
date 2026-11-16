*! Full-original-F conditional grid prediction through the ONE fixed native FV kernel.
version 19.5
mata:
struct cdid_pt_grid_record {
 string scalar key,status,reason,source
 real scalar date,design_index,reference_rows
 real matrix raw,fit_raw,Z,missing,unsupported
 string rowvector names,columns,factors
 struct cdid_ref_design scalar frozen
 struct cdid_cpt_grid scalar grid
}
struct cdid_pt_grid_set {
 struct cdid_pt_grid_record colvector records
 real colvector cell_record
 string colvector keys,id_labels
 real scalar retained,cap,work_ops
}
void _cdid_pt_grids_start(string scalar name,real scalar n,struct cdid_population scalar P,real scalar cap,real scalar retained)
{
 struct cdid_pt_grid_set scalar G
 pointer() scalar p
 G.records=cdid_pt_grid_record(0,1);G.keys=J(0,1,"");G.cell_record=J(n,1,0)
 G.id_labels=(P.ids.string_id ? P.ids.string_id_value : strofreal(P.ids.numeric_id,"%24.17g"))
 G.retained=_cdid_base_cells(n,2,cap,retained);G.cap=cap;G.work_ops=0
 p=crexternal(name);if(p==NULL)_cdid_pt_fail("Conditional grid collection object already exists.",498)
 *p=G
}
void _cdid_pt_grid_request(string scalar collection,string scalar bundle,string scalar rawname,
 string scalar pending,string scalar encoded,real scalar cell,real scalar date,real scalar design)
{
 pointer(struct cdid_pt_grid_set scalar) scalar pG
 pointer(struct cdid_public_bundle scalar) scalar pB
 pointer(struct cdid_reference_inputs scalar) scalar praw
 pointer(struct cdid_population scalar) scalar pP
 struct cdid_pt_grid_record scalar q
 struct cdid_reference_inputs scalar S,E
 struct cdid_ref_design scalar D
 real colvector first,rr,ix,bad,units,keys
 real rowvector vars
 real scalar i,j,u,slot,found,datum,admitted,unit_static
 pointer() scalar p
 pG=findexternal(collection);pB=findexternal(bundle);praw=findexternal(rawname);pP=praw->population
 slot=pB->reference.roles[1].slot;D=pB->reference.designs[slot]
 q.source=(pB->reference.roles[1].dated ? "dated" : "static")
 q.names=cdid_ref_formula_names(D.formula)
 vars=_cdid_ref_names(praw->raw_names,invtokens(q.names))
 unit_static=0
 if(q.source=="static") {
  unit_static=(cols(vars)==0 | all(pB->prepared.required_X:!=0) |
   (pP->ids.frame_supplied & all(praw->has_frame[vars])))
 }
 q.key=(q.source=="dated" ? "dated:"+strofreal(design,"%21.0f") :
  (unit_static ? "static:unit" : "static:"+strofreal(date,"%21x")))
 found=0
 for(i=1;i<=rows(pG->keys);i++)if(pG->keys[i]==q.key)found=i
 st_local("pt_grid_cached",strofreal(found>0))
 if(found){pG->cell_record[cell]=found;return;}
 S=pB->reference;S.population=pP
 q.frozen=D;q.date=date;q.design_index=design;q.status="available";q.reason=""
 q.columns=D.columns;q.factors=J(1,0,"")
 q.missing=q.unsupported=J(0,4,.)
 pG->work_ops=pG->work_ops+64*(pP->N+rows(D.keys))*(cols(q.names)+cols(D.value)+4)*(1+ceil(log(pP->R+pP->N+2)/log(2)))
 admitted=_cdid_base_cells(pP->N,2*cols(q.names)+cols(D.value)+4,pG->cap,pG->retained)
 admitted=_cdid_base_cells(rows(D.keys),cols(q.names),pG->cap,admitted)
 q.raw=J(pP->N,cols(q.names),.)
 if(q.source=="dated") {
  vars=_cdid_ref_names(S.raw_names,invtokens(q.names));keys=design*(S.N+1):+(1::pP->N)
  q.raw=cdid_ref_raw_query(S,keys,vars)
  q.fit_raw=cdid_ref_raw_query(S,D.keys,vars);q.reference_rows=rows(D.keys)
 }
 else {
  first=pP->first_row
  if(unit_static)ix=first
  else {
   rr=cdid_population_rows(*pP,(1::pP->N),J(pP->N,1,date))
   ix=(pB->prepared.required_X:!=0):*first:+(pB->prepared.required_X:==0):*rr
  }
  units=_cdid_ref_ix(ix:>0)
  for(j=1;j<=cols(q.names);j++) {
   if(pP->ids.frame_supplied & praw->has_frame[vars[j]])q.raw[,j]=praw->attributes[,vars[j]]
   else if(praw->has_data[vars[j]] & rows(units))q.raw[units,j]=praw->data[ix[units],vars[j]]
  }
  q.fit_raw=q.raw[_cdid_ref_ix(pB->prepared.required_X),];q.reference_rows=rows(q.fit_raw)
 }
 for(j=1;j<=cols(q.names);j++) {
  bad=_cdid_ref_ix(q.raw[,j]:>=.)
  for(i=1;i<=rows(bad);i++) {
   u=bad[i];datum=.
   if(q.source=="dated") {
    vars=_cdid_ref_names(S.raw_names,q.names[j]);if(S.kind[vars[1]]!=1)datum=date
   }
   else {
    vars=_cdid_ref_names(praw->raw_names,q.names[j])
    if(!(pP->ids.frame_supplied & praw->has_frame[vars[1]]) & !pB->prepared.required_X[u])datum=date
   }
   q.missing=q.missing\(u,datum,j,.)
  }
 }
 if(rows(q.missing)){q.status="unavailable";q.reason="Missing covariate values required by the full-F conditional grid; no carry-forward, substitution or deletion.";}
 q.frozen.keys=(1::pP->N);q.frozen.encoded=0;q.frozen.raw_names=q.names
 q.frozen.raw_index=(cols(q.names) ? (1..cols(q.names)) : J(1,0,0))
 q.frozen.value=J(0,0,.)
 E.designs=q.frozen;E.n_slots=1;E.max_entries=pG->cap
 E.retained_entries=_cdid_base_cells(pP->N,cols(q.names)+cols(D.value)+4,pG->cap,pG->retained)
 E.retained_entries=_cdid_base_cells(rows(q.fit_raw),cols(q.names),(*pG).cap,E.retained_entries)
 p=crexternal(encoded);if(p==NULL)_cdid_pt_fail("Grid encoding object already exists.",498)
 *p=E
 p=crexternal(pending);if(p==NULL)_cdid_pt_fail("Pending grid object already exists.",498)
 *p=q
 st_local("pt_grid_expanded",invtokens(D.expanded));st_local("pt_grid_ready",strofreal(q.status=="available"))
}
void _cdid_pt_grid_factor(string scalar pending,string scalar variable)
{
 struct cdid_pt_grid_record scalar q
 real rowvector var
 real colvector allowed,bad
 real scalar j,u
 q=*findexternal(pending)
 if(any(q.factors:==variable))return
 var=_cdid_ref_names(q.names,variable);j=var[1]
 q.factors=q.factors,variable
 allowed=uniqrows(sort(q.fit_raw[,j],1));bad=_cdid_ref_ix(cdid_ref_lookup(q.raw[,j],allowed):==0)
 for(u=1;u<=rows(bad);u++)if(q.raw[bad[u],j]<.)
  q.unsupported=q.unsupported\(bad[u],(q.source=="dated" ? q.date : .),j,q.raw[bad[u],j])
 if(rows(q.unsupported)){q.status="unavailable";q.reason="Grid covariate levels absent from the frozen fitting design; no level is added.";}
 *findexternal(pending)=q;st_local("pt_grid_ready",strofreal(q.status=="available"))
}
void _cdid_pt_grid_stage(string scalar pending)
{
 struct cdid_pt_grid_record scalar q
 real scalar j
 q=*findexternal(pending)
 if(st_nobs()!=0 | st_nvar()!=0)_cdid_pt_fail("Conditional prediction needs its empty owned frame.",498)
 st_addobs(rows(q.raw))
 for(j=1;j<=cols(q.raw);j++){(void)st_addvar("double",q.names[j]);st_store(.,q.names[j],q.raw[,j]);}
}
struct cdid_cpt_grid scalar cdid_pt_grid_create(real matrix Z,struct cdid_cpt_budget scalar budget)
{
 struct cdid_cpt_grid scalar g
 real scalar used
 used=_cdid_base_cells(rows(Z),cols(Z)+4,budget.max_workspace_cells,0)
 if(cols(Z)>0)if(all(Z:==Z[1,])) {
  /* Exact qualified constant-grid adapter: charge the real design envelope,
     use the singleton zero-column grid, then restore p/zd. */
  g=cdid_cpt_grid_create(J(rows(Z),0,.),budget);g.p=cols(Z);g.zd=Z[1,]
  return(g)
 }
 return(cdid_cpt_grid_create(Z,budget))
}
void _cdid_pt_grid_finish(string scalar collection,string scalar bundle,string scalar pending,
 string scalar encoded,real scalar cell,real scalar design,struct cdid_cpt_budget scalar budget)
{
 struct cdid_pt_grid_set scalar G
 struct cdid_pt_grid_record scalar q
 struct cdid_reference_inputs scalar E,S
 struct cdid_public_bundle scalar B
 real colvector hit,on,units
 real scalar j,slot
 G=*findexternal(collection);q=*findexternal(pending);B=*findexternal(bundle);S=B.reference
 slot=S.roles[1].slot
 if(q.status=="available") {
  E=*findexternal(encoded);q.Z=E.designs[1].value
  if(cols(q.Z)!=cols(S.designs[slot].value) | any(E.designs[1].columns:!=q.columns) | missing(q.Z)) {
   q.status="unavailable";q.reason="Frozen grid prediction is nonfinite or changes the fitted column map."
  }
  else if(q.source=="static") {
   units=S.designs[slot].keys
   if(!all(q.Z[units,]:==S.designs[slot].value)) {
    q.status="unavailable";q.reason="Static grid prediction does not reproduce the frozen fitting design."
   }
  }
  else {
   hit=cdid_ref_lookup(design*(S.N+1):+(1::S.N),S.designs[slot].keys);on=_cdid_ref_ix(hit:>0)
   if(rows(on))if(!all(q.Z[on,]:==S.designs[slot].value[hit[on],])) {
    q.status="unavailable";q.reason="Dated grid prediction does not reproduce canonical fitted rows."
   }
  }
  if(q.status=="available")q.grid=cdid_pt_grid_create(q.Z,budget)
 }
 if(q.status=="available") {
  G.retained=_cdid_base_cells(rows(q.Z),cols(q.Z)+1,G.cap,G.retained)
  G.retained=_cdid_base_cells(q.grid.n_distinct,cols(q.Z)+3,G.cap,G.retained)
 }
 q.raw=q.fit_raw=J(0,0,.)
 G.records=G.records\q;G.keys=G.keys\q.key;G.cell_record[cell]=rows(G.keys)
 *findexternal(collection)=G
}
end
