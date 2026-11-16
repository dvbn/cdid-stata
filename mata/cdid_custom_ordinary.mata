*! Private ordinary custom pretest producer.
*! Actual hY callback at EVERY original grid-row occurrence; no linear IF shortcut.
*! Fixed did2.5.1 compute.att_gt/conditional_did_pretest.
version 19.5
mata:
struct cdid_co_call {
 pointer(function) scalar callback
 pointer() scalar payload
 struct cdid_callback_request scalar request
 struct cdid_callback_result scalar result
 real scalar returned
}
struct cdid_co_response {
 real scalar available,rc,ATT
 real colvector psi
 string scalar reason
}
struct cdid_co_process {
 struct cdid_oc_xproc scalar process
 real scalar N,M,K,n_callbacks,available,workspace_entries,work_ops
 real matrix cells,cell_rc
 string matrix cell_reason
 string scalar membership_law,status
}
/* Stata's documented non-aborting _stata() allows an ordinary cell callback
   exception to be recorded without applying the strict chained validator.
   The only external object is this private temporary call owner; no scientific
   population pointer is exposed. Always release it after the protected call. */
void _cdid_co_call_store(string scalar name)
{
 pointer(struct cdid_co_call scalar) scalar p
 p=findexternal(name)
 if(p==NULL)_cdid_callback_fail("owner","Ordinary temporary callback object is absent.",498)
 (*p).result=cdid_callback_call((*p).callback,(*p).request,(*p).payload)
 (*p).returned=1
}
struct cdid_co_response scalar cdid_co_call_cell(pointer(function) scalar callback,
 struct cdid_callback_request scalar q,pointer() scalar payload)
{
 struct cdid_co_call scalar holder
 struct cdid_co_response scalar out
 pointer(struct cdid_co_call scalar) scalar p
 string scalar name,rhold
 real scalar rc,m
 name=st_tempname();rhold=st_tempname()
 if(findexternal(name)!=NULL)_cdid_callback_fail("owner","Ordinary temporary name already exists.",498)
 holder.callback=callback;holder.payload=payload;holder.request=q;holder.returned=0
 p=crexternal(name);*p=holder
 stata("_return hold "+rhold,1,1)
 rc=_stata("mata: _cdid_co_call_store("+char(34)+name+char(34)+")",1,1)
 stata("_return restore "+rhold,1,1)
 holder=*p;rmexternal(name)
 out.available=0;out.rc=rc;out.ATT=.;out.psi=J(0,1,.);out.reason="callback_error"
 if(rc!=0 | !holder.returned)return(out)
 /* Ordinary ignores optional gradient output. ATT/IF contract failures and
    NaN/missing ATT become cell-unavailable, never a chained-global validator.
    The full grid is still evaluated before family disposition. */
 if(holder.result.has_ATT!=1 | rows(holder.result.ATT)!=1 | cols(holder.result.ATT)!=1) {
  out.rc=503;out.reason="ordinary_ATT_shape";return(out)
 }
 if(missing(holder.result.ATT)) {out.rc=504;out.reason="ordinary_missing_ATT";return(out);}
 out.ATT=holder.result.ATT[1,1];m=rows(q.D)
 if(q.inffunc) {
  if(holder.result.has_psi!=1 | rows(holder.result.psi)!=m | cols(holder.result.psi)!=1) {
   out.rc=503;out.reason="ordinary_IF_shape";return(out)
  }
  if(missing(holder.result.psi)){out.rc=504;out.reason="ordinary_missing_IF";return(out);}
  out.psi=holder.result.psi
 }
 out.available=1;out.rc=0;out.reason=""
 return(out)
}
/* Input is the actual old-preprocessed ordinary recipe, including evaluated
   X and first-calendar panel/all-row RCS grid, complete raw-X screening and
   exact supplied/generated identities. This constructor performs no built-in
   method, overlap/rank, treatment or outcome fit. Public intake stays separate. */
struct cdid_co_process scalar cdid_co_prepare(struct cdid_oc_recipe scalar rc,
 struct cdid_cpt_budget scalar budget,string rowvector columns,pointer(function) scalar callback,pointer() scalar payload)
{
 struct cdid_co_process scalar out
 struct cdid_callback_request scalar q
 struct cdid_co_response scalar ans
 real scalar N,R,T,M,K,C,m,c,p,j,tp,bp,g,at,cap,used,n1,nt
 real colvector unit_g,row_tp,lin,ur,group,mark,h,take,rows0,rows1,units,Dv,Cv,idx,cluster
 real matrix unit_row,X0
 if(callback==NULL)_cdid_callback_fail("input","Ordinary custom callback function pointer is required.",498)
 R=rows(rc.row_y);T=rows(rc.tlist);K=rows(rc.cells)
 if((rc.panel!=0 & rc.panel!=1) | R<1 | T<2 | K<1 | cols(rc.cells)!=3 |
  cols(rc.row_y)!=1 | cols(rc.row_period)!=1 | cols(rc.row_g)!=1 | cols(rc.row_w)!=1 | cols(rc.tlist)!=1 |
  rows(rc.row_period)!=R | rows(rc.row_g)!=R | rows(rc.row_w)!=R | rows(rc.X)!=R | cols(columns)!=cols(rc.X) |
  missing(rc.row_y) | missing(rc.row_period) | missing(rc.row_g) | missing(rc.row_w) | missing(rc.X) |
  any(rc.row_w:<0) | mean(rc.row_w)<=0 | missing(rc.cells) | missing(rc.tlist))
  _cdid_callback_fail("ordinary_intake","Require finite actual old-preprocessed rows, raw nonnegative weights and declared cells.",498)
 if(any(rc.row_g:<0) | any(rc.tlist[(2::T)]:<=rc.tlist[(1::(T-1))]))
  _cdid_callback_fail("ordinary_intake","Source cohorts are zero/positive and the source calendar must strictly increase.",498)
 /* Validate the complete declared family before invoking any callback. */
 for(c=1;c<=K;c++) {
  tp=_cdid_oc_tindex(rc.tlist,rc.cells[c,2]);bp=_cdid_oc_tindex(rc.tlist,rc.cells[c,3])
  if(tp<2 | bp!=tp-1 | rc.cells[c,1]<=rc.cells[c,2])
   _cdid_callback_fail("ordinary_cell","Require actual adjacent varying-base pretreatment cells.",498)
 }
 if(rc.control!="nevertreated" & rc.control!="notyettreated")
  _cdid_callback_fail("ordinary_control","Invalid ordinary comparison control.",498)
 row_tp=J(R,1,0)
 for(tp=1;tp<=T;tp++)row_tp=row_tp+tp:*(rc.row_period:==rc.tlist[tp])
 if(any(row_tp:<1))_cdid_callback_fail("ordinary_time","Retained rows must use the actual source calendar.",498)
 if(rc.panel) {
  if(rows(rc.row_unit)!=R | missing(rc.row_unit) | any(rc.row_unit:<1) | any(rc.row_unit:!=floor(rc.row_unit)))
   _cdid_callback_fail("ordinary_identity","Panel rows require original unit positions.",498)
  N=max(rc.row_unit);lin=(rc.row_unit:-1):*T:+row_tp
  if(N*T!=R | rows(uniqrows(sort(lin,1)))!=R)_cdid_callback_fail("ordinary_panel","The old-preprocessed panel must be balanced.",498)
  ur=J(R,1,0);ur[lin]=(1::R);unit_row=rowshape(ur,N);unit_g=rc.row_g[unit_row[,1]]
  if(any(rc.row_g:!=unit_g[rc.row_unit]))_cdid_callback_fail("ordinary_cohort","Panel cohorts must be time invariant.",498)
  group=(1::N)
 }
 else {
  N=R;unit_g=rc.row_g;unit_row=J(0,0,.)
  if(rows(rc.key)) {
   if(rows(rc.key)!=N | missing(rc.key) | rows(uniqrows(sort((rc.key,row_tp),(1,2))))!=N |
    rows(uniqrows(sort((rc.key,rc.row_g),(1,2))))!=rows(uniqrows(sort(rc.key,1))))
    _cdid_callback_fail("ordinary_identity","Supplied RCS identities must preserve source cohort and unique period rows.",498)
   group=cdid_oc_first_positions(rc.key)
  }
  else group=(1::N)
 }
 M=rows(rc.Z);C=max(rc.cmap)
 if(M!=N | missing(rc.Z) | rows(rc.cmap)!=N | missing(rc.cmap) | any(cdid_oc_first_positions(rc.cmap):!=rc.cmap))
  _cdid_callback_fail("ordinary_grid","Actual first-calendar panel/all-RCS grid and source cluster order are required.",498)
 if(N>budget.max_units | C>budget.max_clusters | K>budget.max_cells | M>budget.max_distinct_points | cols(rc.X)>budget.max_coords)
  _cdid_callback_fail("ordinary_bound","Full source grid occurrences/cells exceed existing explicit process limits.",908)
 cap=_cdid_oc_cap(budget);used=_cdid_oc_charge(R,cols(rc.X)+8,cap,0)
 used=_cdid_oc_charge(N,T+cols(rc.Z)+8,cap,used)
 used=_cdid_oc_charge(C,M*K,cap,used);used=_cdid_oc_charge(M,2*K+1,cap,used)
 used=_cdid_oc_charge(4*N,cols(rc.X)+8,cap,used)
 out.workspace_entries=used;out.work_ops=N*M*(cols(rc.Z)+K*(cols(rc.X)+8))
 if(budget.max_work_ops<.)if(out.work_ops>budget.max_work_ops)_cdid_callback_fail("ordinary_work","Actual full-grid callback work exceeds existing explicit cap.",908)
 /* Multiplicity1 for EACH source occurrence, including equal thresholds.
    A stateful nonlinear callback is never deduplicated on geometric X equality. */
 out.process=cdid_oc_xp_start(N,J(M,1,1),K,rc.cmap,budget)
 out.N=N;out.M=M;out.K=K;out.cells=rc.cells;out.n_callbacks=0;out.available=1
 out.cell_rc=J(M,K,0);out.cell_reason=J(M,K,"")
 out.membership_law=(rc.panel ? "panel first-calendar ID membership" : "RCS union over actual supplied/generated ID grid rows")
 for(m=1;m<=M;m++) {
  h=J(N,1,1)
  for(p=1;p<=cols(rc.Z);p++)h=h:*(rc.Z[,p]:<=rc.Z[m,p])
  if(!rc.panel) {
   mark=J(max(group),1,0)
   for(j=1;j<=N;j++)if(h[j])mark[group[j]]=1
   h=mark[group]
  }
  for(c=1;c<=K;c++) {
   g=rc.cells[c,1];tp=_cdid_oc_tindex(rc.tlist,rc.cells[c,2]);bp=_cdid_oc_tindex(rc.tlist,rc.cells[c,3])
   if(tp<2 | bp!=tp-1 | g<=rc.cells[c,2])_cdid_callback_fail("ordinary_cell","Require actual adjacent varying-base pretreatment cells.",498)
   Dv=(unit_g:==g);Cv=(unit_g:==0)
   if(rc.control=="notyettreated")Cv=Cv:|((unit_g:>rc.cells[c,2]):&(unit_g:!=g))
   q=cdid_callback_request_init();q.mode="ordinary";q.route=(rc.panel ? "panel" : "RCS")
   q.panel=rc.panel;q.inffunc=1;q.gradient_required=0;q.N=N;q.call_index=out.n_callbacks+1
   q.grid_index=m;q.model_index=c;q.raw_index=m;q.g=g;q.s=rc.cells[c,3];q.t=rc.cells[c,2]
   q.link_id="ordinary_cell:"+strofreal(c)+":grid:"+strofreal(m)
   if(rc.panel) {
    units=vec(selectindex(Dv:|Cv));rows0=unit_row[units,bp];rows1=unit_row[units,tp]
    q.rows_F=units;q.D=Dv[units];q.y0=rc.row_y[rows0]:*h[units];q.y1=rc.row_y[rows1]:*h[units]
    q.covariate_names=columns;q.covariates=rc.X[rows0,];q.i_weights=rc.row_w[rows0]
   }
   else {
    units=vec(selectindex((Dv:|Cv):&((row_tp:==tp):|(row_tp:==bp))))
    q.rows_F=units;q.D=Dv[units];q.post=(row_tp[units]:==tp);q.y=rc.row_y[units]:*h[units]
    q.covariate_names=columns;q.covariates=rc.X[units,];q.i_weights=rc.row_w[units]
    if(sum(q.D:*q.post)==0 | sum(q.D:*(1:-q.post))==0 | sum((1:-q.D):*q.post)==0 | sum((1:-q.D):*(1:-q.post))==0) {
     out.available=0;out.cell_rc[m,c]=459;out.cell_reason[m,c]="source_empty_group_period";continue
    }
   }
   n1=rows(q.D);out.n_callbacks=out.n_callbacks+1
   ans=cdid_co_call_cell(callback,q,payload)
   if(!ans.available) {
    out.available=0;out.cell_rc[m,c]=ans.rc;out.cell_reason[m,c]=ans.reason
    errprintf("Ordinary custom callback cell (%g,%g), grid %g is unavailable: %s\n",g,q.t,m,ans.reason)
    continue
   }
   cdid_oc_xp_set_j(out.process,c,m,ans.ATT)
   cdid_oc_xp_add(out.process,c,units,m,(N/n1):*ans.psi)
  }
 }
 out.status=(out.available ? "prepared" : "ordinary_family_unavailable")
 return(out)
}
end
