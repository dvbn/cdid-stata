*! recipe-indexed baseline provider.
*! Source R cdid/; read reference inputs; no encoding/main dispatch.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_dbase(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_db_plan {
 real scalar n_F,n_A,n_consumers,n_pairs,n_supports,n_blocks,workspace_entries,work_ops,work_limit
 real colvector calendar,cohorts,A_to_F,F_to_A,a_cohort,active_cohorts,cohort_to_active
 real matrix pair,support,consumer
 real colvector pair_recipe,block_pair,block_ptr,block_support
 real matrix row_AT,pair_cum,w_hit_A,x1_hit_A,x2_hit_A
 real colvector pair_a,pair_rs,pair_rt,pair_w_col,pair_x1_col,pair_x2_col
 real scalar ref_n_wkeys,ref_slot1_nkeys,ref_slot2_nkeys
 real colvector pair_ctx
 real matrix elig_A
}
struct cdid_db_or_fit {
    string scalar outcome_pooling
    struct cdid_base_or_block colvector blocks
    real scalar n_F,n_A,n_supports,n_blocks,n_consumers,intercept,p,n_fit_rows,objective,retained_entries,peak_entries
    real colvector repeat_mass,pair_mass,sumY,rates,rates_c,rates_active,alpha,beta,block_mass,alpha_rate
    real matrix sumZ,mean_D,R,alias_map,alias_constants,beta_score_A,beta_rate_loading,link_c
    real colvector link_mu,link_d,link_rate_loading
    real rowvector kept,dropped,pivot
}
struct cdid_db_ps_cache {
    real scalar rep_ell,pi
    real colvector risk_a,risk_pos_a,raw_probability,odds,cmap
}
struct cdid_db_nu_entry {
    real scalar rep_ell,pi_rep
    struct cdid_ipw_nuisance scalar nu
    real colvector t_ix,c_ix,fw_ctrl,cmap
}
struct cdid_db_provider {
 string scalar method,status,weight_row_domain
 real scalar n_F,n_A,core_to_frame,normalization_count,weight_rows,sample_rank_upper
 real scalar n_ps_fits,n_ps_reuses,n_or_fits,workspace_entries,workspace_peak_entries,workspace_limit,chunk_rows,work_ops
 pointer(struct cdid_population scalar) scalar population
 pointer(struct cdid_reference_inputs scalar) scalar reference
 struct cdid_db_plan scalar plan
 struct cdid_reference_catalog scalar catalog
 real colvector points,group_mass
 struct cdid_db_or_fit scalar outcome
 struct cdid_base_ps colvector propensity
 struct cdid_base_point colvector point
 struct cdid_db_ps_cache colvector ps_cache
 real colvector ps_class_id
}
real scalar _cdid_db_pos(real colvector v,real scalar x)
{
 real colvector hit
 hit=_cdid_base_index(v:==x)
 if(rows(hit)!=1)_cdid_base_fail("dated_identity","A required cohort is not uniquely mapped.",498)
 return(hit[1])
}
/* Eligibility context of a pair's recipe (missingreference(exclude)): A rows
   of the pair's eligible units; all rows when the recipe is unrestricted. */
real colvector _cdid_db_keep_mask(struct cdid_db_plan scalar p,real scalar pi,real colvector ai)
{
 if(rows(p.pair_ctx)!=p.n_pairs | pi<1 | pi>p.n_pairs) return(J(rows(ai),1,1))
 if(p.pair_ctx[pi]==0 | rows(ai)==0) return(J(rows(ai),1,1))
 return(p.elig_A[ai,p.pair_ctx[pi]])
}
real colvector _cdid_db_keep(struct cdid_db_plan scalar p,real scalar pi,real colvector ai)
{
 if(rows(p.pair_ctx)!=p.n_pairs | pi<1 | pi>p.n_pairs) return(ai)
 if(p.pair_ctx[pi]==0 | rows(ai)==0) return(ai)
 return(select(ai,p.elig_A[ai,p.pair_ctx[pi]]))
}
struct cdid_db_plan scalar _cdid_db_plan_make(struct cdid_baseline_plan_state scalar original,
 struct cdid_reference_catalog scalar C,real scalar cap,real scalar retained,real scalar workcap)
{
 struct cdid_db_plan scalar p
 transmorphic scalar pairs,supports
 real scalar i,j,b,pi,s,k,g,used,np,ns,nbs,n_coh
 real colvector K,rec_buf,bsup_buf
 real matrix pair_buf,sup_buf
 string scalar key,sk
 /* Preflight the numeric registry envelope before dynamic concatenations.
    Shared source/original plan/catalog owners are charged by the factory. */
 used=_cdid_base_cells(8,3*C.n_consumers+rows(C.K_index)+C.n_blocks+rows(C.block_K_index),cap,retained)
 used=_cdid_base_cells(1,4*original.n_A+original.n_F+3*rows(original.cohorts)+rows(original.calendar)+20,cap,used)
 p.workspace_entries=used-retained;p.work_limit=workcap
 p.work_ops=_cdid_base_work(0,workcap,original.work_ops)
 p.n_F=original.n_F;p.n_A=original.n_A;p.calendar=original.calendar;p.cohorts=original.cohorts
 p.A_to_F=original.A_to_F;p.F_to_A=original.F_to_A;p.a_cohort=original.a_cohort
 p.active_cohorts=original.active_cohorts;p.cohort_to_active=original.cohort_to_active
 p.n_consumers=C.n_consumers;p.n_blocks=C.n_blocks
 n_coh=rows(p.cohort_to_active)
 pair_buf=J(C.n_consumers+C.n_blocks,2,.);rec_buf=J(C.n_consumers+C.n_blocks,1,.)
 sup_buf=J(C.n_consumers+rows(C.K_index),2,.)
 p.consumer=J(C.n_consumers,4,.)
 np=0;ns=0
 pairs=asarray_create("string",1);supports=asarray_create("string",1)
 for(i=1;i<=C.n_consumers;i++) {
  p.work_ops=_cdid_base_work(1+C.K_ptr[i+1]-C.K_ptr[i],p.work_limit,p.work_ops)
  key=strofreal(C.consumers[i,2])+":"+strofreal(C.consumers[i,3])+":"+strofreal(C.consumers[i,4])
  if(!asarray_contains(pairs,key)) {
   np++;asarray(pairs,key,np);pair_buf[np,]=C.consumers[i,(2,3)];rec_buf[np]=C.consumers[i,4]
  }
  pi=asarray(pairs,key);K=C.consumers[i,1]\C.K_index[C.K_ptr[i]..C.K_ptr[i+1]-1]
  for(j=1;j<=rows(K);j++) {
   k=K[j]
   if(k>=1 & k<=n_coh) {
    g=p.cohort_to_active[k]
    if(g<=0) g=_cdid_db_pos(p.active_cohorts,k)
   }
   else g=_cdid_db_pos(p.active_cohorts,k)
   sk=key+":"+strofreal(k)
   if(!asarray_contains(supports,sk)) {
    ns++;asarray(supports,sk,ns);sup_buf[ns,]=(pi,g)
   }
  }
  p.consumer[i,]=(C.consumers[i,1],pi,C.consumer_block[i],asarray(supports,key+":"+strofreal(C.consumers[i,1])))
 }
 p.pair=(np ? pair_buf[|1,1\np,2|] : J(0,2,.))
 p.pair_recipe=(np ? rec_buf[|1\np|] : J(0,1,.))
 p.support=(ns ? sup_buf[|1,1\ns,2|] : J(0,2,.))
 p.block_pair=J(p.n_blocks,1,.);p.block_ptr=J(p.n_blocks+1,1,1)
 bsup_buf=J(rows(C.block_K_index),1,.);nbs=0
 for(b=1;b<=p.n_blocks;b++) {
  p.work_ops=_cdid_base_work(1+C.block_K_ptr[b+1]-C.block_K_ptr[b],p.work_limit,p.work_ops)
  key=strofreal(C.blocks[b,1])+":"+strofreal(C.blocks[b,2])+":"+strofreal(C.blocks[b,3]);p.block_pair[b]=asarray(pairs,key)
  K=C.block_K_index[C.block_K_ptr[b]..C.block_K_ptr[b+1]-1]
  for(j=1;j<=rows(K);j++) {nbs++;bsup_buf[nbs]=asarray(supports,key+":"+strofreal(K[j]));}
  p.block_ptr[b+1]=nbs+1
 }
 p.block_support=(nbs ? bsup_buf[|1\nbs|] : J(0,1,.))
 p.n_pairs=np;p.n_supports=ns
 p.ref_n_wkeys=-1;p.ref_slot1_nkeys=-1;p.ref_slot2_nkeys=-1
 return(p)
}
real colvector _cdid_db_sorted_lookup(real colvector keys,real colvector table)
{
 real colvector hit
 real scalar n,nt,i,j,lo,hi,mid,k1
 n=rows(keys);nt=rows(table);hit=J(n,1,0)
 if(n==0 | nt==0) return(hit)
 k1=keys[1];lo=1;hi=nt
 while(lo<=hi) {
  mid=floor((lo+hi)/2)
  if(table[mid]<k1) lo=mid+1
  else hi=mid-1
 }
 j=lo
 for(i=1;i<=n;i++) {
  k1=keys[i]
  while(j<=nt) {
   if(table[j]>=k1) break
   j++
  }
  if(j>nt) break
  if(table[j]==k1) {hit[i]=j;j++;}
 }
 return(hit)
}
void _cdid_db_build_caches(struct cdid_population scalar P,struct cdid_reference_inputs scalar S,
 struct cdid_db_plan scalar p,real matrix orig_row_AT)
{
 real scalar pi,s,t,tot,n_pi,T,c,slot1,slot2
 real colvector rs,rt,ix,ix2,ai,has_a,anc_pi,u_anc,anc_map,d_pi,u_d,d_map
 T=rows(P.calendar)
 p.row_AT=(rows(orig_row_AT)==P.n_A & cols(orig_row_AT)==T ? orig_row_AT : _cdid_base_build_row_AT(P))
 p.pair_cum=J(0,0,0);p.pair_a=J(0,1,0);p.pair_rs=J(0,1,0);p.pair_rt=J(0,1,0)
 p.w_hit_A=J(0,0,0);p.x1_hit_A=J(0,0,0);p.x2_hit_A=J(0,0,0)
 p.pair_w_col=J(0,1,0);p.pair_x1_col=J(0,1,0);p.pair_x2_col=J(0,1,0)
 p.ref_n_wkeys=-1;p.ref_slot1_nkeys=-1;p.ref_slot2_nkeys=-1
 if(p.n_pairs<=0 | P.n_A<=0) return
 if(rows(p.row_AT)==P.n_A & P.n_A*p.n_pairs<=4000000) {
  p.pair_cum=J(P.n_A+1,p.n_pairs,0);tot=0
  for(pi=1;pi<=p.n_pairs;pi++) {
   s=p.pair[pi,1];t=p.pair[pi,2]
   rs=p.row_AT[,s];rt=p.row_AT[,t]
   ix=_cdid_base_index((rs:>0):&(rt:>0));has_a=J(P.n_A,1,0)
   if(rows(ix)) {
    ix2=_cdid_base_index((P.y[rs[ix]]:<.):&(P.y[rt[ix]]:<.))
    if(rows(ix2)) {ai=_cdid_db_keep(p,pi,ix[ix2]);if(rows(ai)) has_a[ai]=J(rows(ai),1,1);}
   }
   p.pair_cum[,pi]=tot:+(0\runningsum(has_a));tot=p.pair_cum[P.n_A+1,pi]
  }
  if(tot>0) {
   p.pair_a=J(tot,1,0);p.pair_rs=J(tot,1,0);p.pair_rt=J(tot,1,0)
   for(pi=1;pi<=p.n_pairs;pi++) {
    n_pi=p.pair_cum[P.n_A+1,pi]-p.pair_cum[1,pi]
    if(n_pi<=0) continue
    s=p.pair[pi,1];t=p.pair[pi,2]
    rs=p.row_AT[,s];rt=p.row_AT[,t]
    ix=_cdid_base_index((rs:>0):&(rt:>0))
    ix2=_cdid_base_index((P.y[rs[ix]]:<.):&(P.y[rt[ix]]:<.))
    ai=_cdid_db_keep(p,pi,ix[ix2]);tot=p.pair_cum[1,pi]
    p.pair_a[|tot+1\tot+n_pi|]=ai
    p.pair_rs[|tot+1\tot+n_pi|]=rs[ai]
    p.pair_rt[|tot+1\tot+n_pi|]=rt[ai]
   }
  }
 }
 if(P.n_A*(T+1)>8000000) return
 anc_pi=S.recipes[p.pair_recipe,2]
 u_anc=uniqrows(sort(anc_pi,1))
 anc_map=J(T+1,1,0);anc_map[u_anc:+1]=(1::rows(u_anc))
 p.pair_w_col=anc_map[anc_pi:+1]
 p.w_hit_A=J(P.n_A,rows(u_anc),0)
 for(c=1;c<=rows(u_anc);c++) p.w_hit_A[,c]=_cdid_db_sorted_lookup(u_anc[c]*(S.N+1):+p.A_to_F,S.weight_keys)
 p.ref_n_wkeys=rows(S.weight_keys)
 if(rows(S.roles)>=1) {
  if(S.roles[1].active) {
   slot1=S.roles[1].slot
   if(slot1>=1 & slot1<=S.n_slots) {
    if(S.designs[slot1].encoded) {
     d_pi=(S.roles[1].dated ? S.recipes[p.pair_recipe,1] : J(p.n_pairs,1,0))
     u_d=uniqrows(sort(d_pi,1))
     d_map=J(T+1,1,0);d_map[u_d:+1]=(1::rows(u_d))
     p.pair_x1_col=d_map[d_pi:+1]
     p.x1_hit_A=J(P.n_A,rows(u_d),0)
     for(c=1;c<=rows(u_d);c++) p.x1_hit_A[,c]=_cdid_db_sorted_lookup(u_d[c]*(S.N+1):+p.A_to_F,S.designs[slot1].keys)
     p.ref_slot1_nkeys=rows(S.designs[slot1].keys)
    }
   }
  }
 }
 if(rows(S.roles)>=2) {
  if(S.roles[2].active) {
   slot2=S.roles[2].slot
   if(slot2>=1 & slot2<=S.n_slots) {
    if(S.designs[slot2].encoded) {
     c=0
     if(rows(S.roles)>=1) {
      if(S.roles[1].active & slot2==S.roles[1].slot & S.roles[2].dated==S.roles[1].dated & rows(p.x1_hit_A)==P.n_A) c=1
     }
     if(c) {
      p.x2_hit_A=p.x1_hit_A;p.pair_x2_col=p.pair_x1_col;p.ref_slot2_nkeys=p.ref_slot1_nkeys
     }
     else {
      d_pi=(S.roles[2].dated ? S.recipes[p.pair_recipe,1] : J(p.n_pairs,1,0))
      u_d=uniqrows(sort(d_pi,1))
      d_map=J(T+1,1,0);d_map[u_d:+1]=(1::rows(u_d))
      p.pair_x2_col=d_map[d_pi:+1]
      p.x2_hit_A=J(P.n_A,rows(u_d),0)
      for(c=1;c<=rows(u_d);c++) p.x2_hit_A[,c]=_cdid_db_sorted_lookup(u_d[c]*(S.N+1):+p.A_to_F,S.designs[slot2].keys)
      p.ref_slot2_nkeys=rows(S.designs[slot2].keys)
     }
    }
   }
  }
 }
}
real matrix _cdid_db_pair_slice(struct cdid_population scalar P,struct cdid_db_plan scalar p,
 real scalar pi,real scalar alo,real scalar ahi)
{
 real scalar klo,khi,s,t
 real colvector ai,rs,rt,ix
 real matrix rr3
 if(ahi<alo) return(J(0,3,0))
 if(rows(p.pair_cum)==p.n_A+1 & cols(p.pair_cum)==p.n_pairs & alo>=1 & ahi<=p.n_A) {
  klo=p.pair_cum[alo,pi]+1;khi=p.pair_cum[ahi+1,pi]
  if(khi<klo) return(J(0,3,0))
  return((p.pair_a[|klo\khi|],p.pair_rs[|klo\khi|],p.pair_rt[|klo\khi|]))
 }
 if(rows(p.row_AT)==p.n_A & alo>=1 & ahi<=p.n_A) {
  s=p.pair[pi,1];t=p.pair[pi,2]
  ai=(alo..ahi)';rs=p.row_AT[ai,s];rt=p.row_AT[ai,t]
  ix=_cdid_base_index((rs:>0):&(rt:>0))
  if(rows(ix)==0) return(J(0,3,0))
  ai=ai[ix];rs=rs[ix];rt=rt[ix]
  ix=_cdid_base_index((P.y[rs]:<.):&(P.y[rt]:<.))
  if(rows(ix)==0) return(J(0,3,0))
  ai=ai[ix];rs=rs[ix];rt=rt[ix]
  ix=_cdid_base_index(_cdid_db_keep_mask(p,pi,ai))
  if(rows(ix)==0) return(J(0,3,0))
  return((ai[ix],rs[ix],rt[ix]))
 }
 rr3=_cdid_base_pair_rows(P,p.pair[pi,1],p.pair[pi,2],alo,ahi)
 if(rows(rr3)==0) return(rr3)
 ix=_cdid_base_index(_cdid_db_keep_mask(p,pi,rr3[,1]))
 return(rows(ix) ? rr3[ix,] : J(0,3,0))
}
/* Count the explicit preparation scans using the native work helper.
   These are row/address/group scan units, not CPU time or a FLOP certificate.
   IF replay is immutable and never consumes a hidden cumulative work budget. */
void _cdid_db_pair_work(struct cdid_db_plan scalar p,real scalar n,real scalar groups)
{
 real scalar unit
 unit=2*(1+ceil(ln(rows(p.calendar)+1)/ln(2)))+groups
 p.work_ops=_cdid_base_work(n*unit,p.work_limit,p.work_ops)
}
/* Return exactly the required recipe/date rows; do not fill/read unused values. */
real matrix _cdid_db_X(struct cdid_reference_inputs scalar S,struct cdid_db_plan scalar p,
 real scalar pi,real colvector ai,real scalar role)
{
 real scalar d,slot
 real colvector units,hit
 if(role==1 & pi>=1 & pi<=p.n_pairs & rows(p.x1_hit_A)==p.n_A & rows(S.roles)>=1) {
  if(S.roles[1].active) {
   slot=S.roles[1].slot
   if(slot>=1 & slot<=S.n_slots) {
    if(S.designs[slot].encoded & rows(S.designs[slot].keys)==p.ref_slot1_nkeys & cols(ai)==1 & !missing(ai)) {
     _cdid_ref_check_population(S)
     if(rows(ai)==0) return(J(0,cols(S.designs[slot].value),.))
     if(min(ai)>=1 & max(ai)<=p.n_A) {
      hit=p.x1_hit_A[ai,p.pair_x1_col[pi]]
      if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
      return(S.designs[slot].value[hit,.])
     }
    }
   }
  }
 }
 else if(role==2 & pi>=1 & pi<=p.n_pairs & rows(p.x2_hit_A)==p.n_A & rows(S.roles)>=2) {
  if(S.roles[2].active) {
   slot=S.roles[2].slot
   if(slot>=1 & slot<=S.n_slots) {
    if(S.designs[slot].encoded & rows(S.designs[slot].keys)==p.ref_slot2_nkeys & cols(ai)==1 & !missing(ai)) {
     _cdid_ref_check_population(S)
     if(rows(ai)==0) return(J(0,cols(S.designs[slot].value),.))
     if(min(ai)>=1 & max(ai)<=p.n_A) {
      hit=p.x2_hit_A[ai,p.pair_x2_col[pi]]
      if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
      return(S.designs[slot].value[hit,.])
     }
    }
   }
  }
 }
 units=p.A_to_F[ai];d=(S.roles[role].dated ? S.recipes[p.pair_recipe[pi],1] : 0)
 return(cdid_ref_design_rows(S,role,d*(p.n_F+1):+units))
}
real matrix _cdid_db_Z(struct cdid_reference_inputs scalar S,struct cdid_db_plan scalar p,
 real scalar pi,real colvector ai)
{
 real matrix X
 real scalar slot,nc
 real colvector hit
 if(pi>=1 & pi<=p.n_pairs & rows(p.x1_hit_A)==p.n_A & rows(S.roles)>=1) {
  if(S.roles[1].active) {
   slot=S.roles[1].slot
   if(slot>=1 & slot<=S.n_slots) {
    if(S.designs[slot].encoded & rows(S.designs[slot].keys)==p.ref_slot1_nkeys & cols(ai)==1 & !missing(ai)) {
     _cdid_ref_check_population(S)
     nc=cols(S.designs[slot].value)
     if(rows(ai)==0) return(J(0,(S.roles[1].intercept ? max((0,nc-1)) : nc),.))
     if(min(ai)>=1 & max(ai)<=p.n_A) {
      hit=p.x1_hit_A[ai,p.pair_x1_col[pi]]
      if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
      if(S.roles[1].intercept) return(nc>1 ? S.designs[slot].value[hit,2..nc] : J(rows(ai),0,.))
      return(S.designs[slot].value[hit,.])
     }
    }
   }
  }
 }
 X=_cdid_db_X(S,p,pi,ai,1)
 if(S.roles[1].intercept)return(cols(X)>1 ? X[,2..cols(X)] : J(rows(X),0,.))
 return(X)
}
real colvector _cdid_db_w(struct cdid_reference_inputs scalar S,struct cdid_db_plan scalar p,
 real scalar pi,real colvector ai)
{
 real colvector hit
 if(pi>=1 & pi<=p.n_pairs & rows(p.w_hit_A)==p.n_A & rows(S.weight_keys)==p.ref_n_wkeys & cols(ai)==1 & !missing(ai)) {
  _cdid_ref_check_population(S)
  if(rows(ai)==0) return(J(0,1,0))
  if(min(ai)>=1 & max(ai)<=p.n_A) {
   hit=p.w_hit_A[ai,p.pair_w_col[pi]]
   if(any(hit:==0)) _cdid_ref_fail("A weight query is outside its actual required support.",459)
   return(vec(S.weights[hit]))
  }
 }
 return(cdid_ref_weight_rows(S,p.pair_recipe[pi],p.A_to_F[ai]))
}
/* Pointer-taking internal readers avoiding *f.reference struct copies on hot paths. */
real colvector _cdid_db_w_ptr(pointer(struct cdid_reference_inputs scalar) scalar ref,
 struct cdid_db_plan scalar p,real scalar pi,real colvector ai)
{
 real colvector hit
 if(rows(ai)==0) return(J(0,1,0))
 if(rows(p.w_hit_A)==p.n_A & rows((*ref).weight_keys)==p.ref_n_wkeys) {
  hit=p.w_hit_A[ai,p.pair_w_col[pi]]
  if(any(hit:==0)) _cdid_ref_fail("A weight query is outside its actual required support.",459)
  return(vec((*ref).weights[hit]))
 }
 return(_cdid_db_w(*ref,p,pi,ai))
}
real matrix _cdid_db_X2_ptr(pointer(struct cdid_reference_inputs scalar) scalar ref,
 struct cdid_db_plan scalar p,real scalar pi,real colvector ai)
{
 real scalar slot
 real colvector hit
 if(rows(p.x2_hit_A)==p.n_A & rows((*ref).roles)>=2) {
  if((*ref).roles[2].active) {
   slot=(*ref).roles[2].slot
   if(slot>=1 & slot<=(*ref).n_slots) {
    if((*ref).designs[slot].encoded & rows((*ref).designs[slot].keys)==p.ref_slot2_nkeys) {
     if(rows(ai)==0) return(J(0,cols((*ref).designs[slot].value),.))
     hit=p.x2_hit_A[ai,p.pair_x2_col[pi]]
     if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
     return((*ref).designs[slot].value[hit,.])
    }
   }
  }
 }
 return(_cdid_db_X(*ref,p,pi,ai,2))
}
real matrix _cdid_db_Z_kept_ptr(pointer(struct cdid_reference_inputs scalar) scalar ref,
 struct cdid_db_plan scalar p,real scalar pi,real colvector ai,real rowvector kept)
{
 real scalar slot,off
 real colvector hit
 if(rows(p.x1_hit_A)==p.n_A & rows((*ref).roles)>=1) {
  if((*ref).roles[1].active) {
   slot=(*ref).roles[1].slot
   if(slot>=1 & slot<=(*ref).n_slots) {
    if((*ref).designs[slot].encoded & rows((*ref).designs[slot].keys)==p.ref_slot1_nkeys) {
     if(rows(ai)==0 | cols(kept)==0) {
      if(rows(ai)>0) {
       hit=p.x1_hit_A[ai,p.pair_x1_col[pi]]
       if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
      }
      return(J(rows(ai),cols(kept),0))
     }
     hit=p.x1_hit_A[ai,p.pair_x1_col[pi]]
     if(any(hit:==0)) _cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
     off=(*ref).roles[1].intercept
     return((*ref).designs[slot].value[hit,off:+kept])
    }
   }
  }
 }
 return(_cdid_db_Z(*ref,p,pi,ai)[,kept])
}
real matrix _cdid_db_read(struct cdid_population scalar P,struct cdid_db_plan scalar p,
    real scalar pi,real scalar alo,real scalar ahi,real colvector groups)
{
    real matrix r
    real colvector dy,keep,cmap
    _cdid_db_pair_work(p,max((0,ahi-alo+1)),rows(groups))
    r=_cdid_db_pair_slice(P,p,pi,alo,ahi)
    if(rows(r)==0)return(J(0,2,0))
    if(rows(groups)==1)keep=_cdid_base_index(p.a_cohort[r[,1]]:==groups[1])
    else {
        cmap=J(rows(p.cohorts),1,0);cmap[groups]=J(rows(groups),1,1)
        keep=_cdid_base_index(cmap[p.a_cohort[r[,1]]])
    }
    if(rows(keep)==0)return(J(0,2,0))
    r=r[keep,]
    dy=P.y[r[,3]]-P.y[r[,2]]
    if(missing(dy))_cdid_base_fail("outcome","Pair outcome change overflowed.",504)
    return((r[,1],dy))
}
real scalar _cdid_db_build_sup_csr(struct cdid_population scalar P,struct cdid_db_plan scalar p,
 real colvector sup_lo,real colvector sup_hi,
 real colvector sup_a,real colvector sup_rs,real colvector sup_rt)
{
 real scalar S,n_act,pi,tot,at,m
 real matrix sup_of_pa,rr,s_info
 real colvector sid,keep,all_sid,all_a,all_rs,all_rt,ord
 S=p.n_supports;n_act=rows(p.active_cohorts)
 sup_lo=J(S,1,1);sup_hi=J(S,1,0)
 sup_a=J(0,1,0);sup_rs=J(0,1,0);sup_rt=J(0,1,0)
 if(S<=0 | p.n_pairs<=0 | n_act<=0 | p.n_A*p.n_pairs>1000000) return(0)
 sup_of_pa=J(n_act,p.n_pairs,0)
 for(m=1;m<=S;m++) sup_of_pa[p.support[m,2],p.support[m,1]]=m
 tot=0
 for(pi=1;pi<=p.n_pairs;pi++) {
  rr=_cdid_db_pair_slice(P,p,pi,1,p.n_A);tot=tot+rows(rr)
 }
 if(tot<=0) return(1)
 all_sid=J(tot,1,0);all_a=J(tot,1,0);all_rs=J(tot,1,0);all_rt=J(tot,1,0)
 at=0
 for(pi=1;pi<=p.n_pairs;pi++) {
  rr=_cdid_db_pair_slice(P,p,pi,1,p.n_A)
  if(rows(rr)==0) continue
  sid=sup_of_pa[p.cohort_to_active[p.a_cohort[rr[,1]]],pi]
  keep=_cdid_base_index(sid:>0);m=rows(keep)
  if(m==0) continue
  all_sid[|at+1\at+m|]=sid[keep]
  all_a[|at+1\at+m|]=rr[keep,1]
  all_rs[|at+1\at+m|]=rr[keep,2]
  all_rt[|at+1\at+m|]=rr[keep,3]
  at=at+m
 }
 if(at<=0) return(1)
 if(at<tot) {
  all_sid=all_sid[|1\at|];all_a=all_a[|1\at|]
  all_rs=all_rs[|1\at|];all_rt=all_rt[|1\at|]
 }
 ord=order((all_sid,all_a),(1,2))
 all_sid=all_sid[ord];sup_a=all_a[ord];sup_rs=all_rs[ord];sup_rt=all_rt[ord]
 s_info=_cdid_ps_dbase(all_sid,1)
 if(rows(s_info)) {
  sid=all_sid[s_info[,1]]
  sup_lo[sid]=s_info[,1];sup_hi[sid]=s_info[,2]
 }
 return(1)
}

void _cdid_db_or_alias(struct cdid_population scalar P,struct cdid_db_plan scalar p,
    struct cdid_reference_inputs scalar S,struct cdid_db_or_fit scalar f,real scalar chunk)
{
    real scalar b,s,pi,L,j,lo,hi,ai,den,err,threshold,cm
    real colvector used,ix,a,v,scale
    real matrix rr,z
    real rowvector maxima,dropped,kept
    if(f.outcome_pooling=="separate") {
        for(b=1;b<=p.n_blocks;b++) {
            dropped=f.blocks[b].dropped;kept=f.blocks[b].kept
            if(cols(dropped)==0)continue
            threshold=64*epsilon(1)*max((1,cols(kept)+cols(dropped)+f.intercept))
            used=J(p.n_supports,1,0)
            used[p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1]]=J(p.block_ptr[b+1]-p.block_ptr[b],1,1)
            ix=_cdid_base_index(p.consumer[,3]:==b)
            if(rows(ix))used[p.consumer[ix,4]]=J(rows(ix),1,1)
            maxima=J(1,cols(dropped),0)
            for(s=1;s<=p.n_supports;s++)if(used[s]) {
                pi=p.support[s,1]
                for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                    _cdid_db_pair_work(p,min((p.n_A,lo+chunk-1))-lo+1,1)
                    rr=_cdid_db_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                    if(rows(rr)==0)continue
                    ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                    if(rows(ix))maxima=colmax(maxima\abs(_cdid_db_Z(S,p,pi,rr[ix,1])[,dropped]))
                }
            }
            for(s=1;s<=p.n_supports;s++)if(used[s]) {
                pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
                for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                    _cdid_db_pair_work(p,min((p.n_A,lo+chunk-1))-lo+1,1)
                    rr=_cdid_db_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                    if(rows(rr)==0)continue
                    ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                    if(rows(ix)==0)continue
                    z=_cdid_db_Z(S,p,pi,rr[ix,1])
                    for(j=1;j<=cols(dropped);j++) {
                        cm=(f.intercept ? f.blocks[b].alias_constants[1,j]/L : 0)
                        v=z[,dropped[j]]-(cols(kept) ? z[,kept]*f.blocks[b].alias_map[,j] : J(rows(z),1,0)):-cm
                        scale=abs(z[,dropped[j]])+(cols(kept) ? abs(z[,kept])*abs(f.blocks[b].alias_map[,j]) : J(rows(z),1,0)):+abs(cm)
                        den=max((maxima[j],2.2250738585072014e-308))
                        err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                        if(err>=. | err>threshold)_cdid_base_fail("alias","Outcome alias does not preserve actual fitting/required prediction support.",459)
                    }
                }
            }
        }
        return
    }
    if(cols(f.dropped)==0)return
    maxima=J(1,cols(f.dropped),0)
    for(s=1;s<=p.n_supports;s++) {
        pi=p.support[s,1]
        for(lo=1;lo<=p.n_A;lo=lo+chunk) {
            _cdid_db_pair_work(p,min((p.n_A,lo+chunk-1))-lo+1,1)
            rr=_cdid_db_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
            if(rows(rr)==0)continue
            ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
            if(rows(ix))maxima=colmax(maxima\abs(_cdid_db_Z(S,p,pi,rr[ix,1])[,f.dropped]))
        }
    }
    threshold=64*epsilon(1)*max((1,cols(f.kept)+cols(f.dropped)+f.intercept))
    for(b=1;b<=p.n_blocks;b++) {
        used=J(p.n_supports,1,0)
        used[p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1]]=J(p.block_ptr[b+1]-p.block_ptr[b],1,1)
        ix=_cdid_base_index(p.consumer[,3]:==b)
        if(rows(ix))used[p.consumer[ix,4]]=J(rows(ix),1,1)
        for(s=1;s<=p.n_supports;s++)if(used[s]) {
            pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
            for(lo=1;lo<=p.n_A;lo=lo+chunk) {
                _cdid_db_pair_work(p,min((p.n_A,lo+chunk-1))-lo+1,1)
                rr=_cdid_db_pair_slice(P,p,pi,lo,min((p.n_A,lo+chunk-1)))
                if(rows(rr)==0)continue
                ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rr[,1]]]:==p.support[s,2])
                if(rows(ix)==0)continue
                z=_cdid_db_Z(S,p,pi,rr[ix,1])
                for(j=1;j<=cols(f.dropped);j++) {
                    cm=(f.intercept ? f.alias_constants[b,j]/L : 0)
                    v=z[,f.dropped[j]]-z[,f.kept]*f.alias_map[,j]:-cm
                    scale=abs(z[,f.dropped[j]])+abs(z[,f.kept])*abs(f.alias_map[,j]):+abs(cm)
                    den=max((maxima[j],2.2250738585072014e-308))
                    err=max(abs(v):/rowmax((scale,J(rows(v),1,den))))
                    if(err>=. | err>threshold)_cdid_base_fail("alias","Outcome alias does not preserve actual fitting/required prediction support.",459)
                }
            }
        }
    }
}

struct cdid_db_or_fit scalar _cdid_db_or_fit(struct cdid_population scalar P,
    struct cdid_db_plan scalar p,struct cdid_reference_inputs scalar Ref,
    struct cdid_baseline_controls scalar ctrl,real scalar retained_base)
{
    struct cdid_db_or_fit scalar out
    struct cdid_base_support_qr colvector summary
    struct cdid_wls_result scalar fit
    pointer(struct cdid_reference_inputs scalar) scalar pref
    pointer(real colvector) colvector coh_a
    real scalar S,B,Jl,d,s,b,m,j,pi,L,lo,hi,first,nn,mass,k,kb,max_k,at,nPseudo,scale,den,used,sum_sa,n_ctrl_b,has_csr
    real colvector groupmass,ai,ww,yy,rr,ix,fw,rawres,sup_lo,sup_hi,sup_a,sup_rs,sup_rt
    real rowvector constant,baseZ,sum_sb
    real matrix rows0,rawZ,D,F,V,Q,R,A,coefrows,resfactor,xmin,xmax,firstZ,values
    real matrix pseudoX,pseudoY,pseudoW,pseudoB,force,centered
    real rowvector tau
    pref=&Ref
    S=p.n_supports;B=p.n_blocks;Jl=p.n_consumers;d=cols(Ref.designs[Ref.roles[1].slot].value)-ctrl.intercept
    if(S<1 | B<1 | Jl<1 | d<0 | (!ctrl.intercept & d==0))
        _cdid_base_fail("outcome_design","Require an encoded nonempty outcome design and supported blocks.",498)
    used=_cdid_base_cells(S,d+2,ctrl.max_entries,0)
    used=_cdid_base_cells(12*used,d+2,ctrl.max_entries,retained_base)
    den=_cdid_base_cells(B,d+2,ctrl.max_entries,0)
    used=_cdid_base_cells(40*den,d+2,ctrl.max_entries,used)
    used=_cdid_base_cells(8*p.n_A,d+2,ctrl.max_entries,used)
    used=_cdid_base_cells(12*S+12*B,Jl>0 ? 1 : 0,ctrl.max_entries,used)
    used=_cdid_base_cells(4*Jl,d+4,ctrl.max_entries,used)
    out.peak_entries=used
    out.outcome_pooling=(ctrl.outcome_pooling=="" ? "separate" : ctrl.outcome_pooling)
    out.n_F=p.n_F;out.n_A=p.n_A;out.n_supports=S;out.n_blocks=B;out.n_consumers=Jl
    out.intercept=ctrl.intercept;out.p=d
    out.pair_mass=J(S,1,0);out.sumY=J(S,1,0);out.sumZ=J(S,d,0);out.rates=J(S,1,0)
    xmin=J(S,d,.);xmax=xmin;firstZ=xmin;summary=cdid_base_support_qr(S,1)
    coh_a=J(rows(p.active_cohorts),1,NULL)
    for(j=1;j<=rows(p.active_cohorts);j++) {
        ix=_cdid_base_index(p.a_cohort:==p.active_cohorts[j])
        coh_a[j]=&(ix:+0)
    }
    groupmass=J(S,1,0)
    for(j=1;j<=S;j++) {
        p.work_ops=_cdid_base_work(p.n_A,p.work_limit,p.work_ops)
        groupmass[j]=sum(_cdid_db_w_ptr(pref,p,p.support[j,1],_cdid_db_keep(p,p.support[j,1],*coh_a[p.support[j,2]])))
    }
    out.repeat_mass=groupmass
    has_csr=_cdid_db_build_sup_csr(P,p,sup_lo=.,sup_hi=.,sup_a=.,sup_rs=.,sup_rt=.)
    for(s=1;s<=S;s++) {
        pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
        if(has_csr) {
            _cdid_db_pair_work(p,p.n_A,1)
            if(sup_hi[s]>=sup_lo[s]) {
                ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                ww=_cdid_db_w_ptr(pref,p,pi,ai)
                out.pair_mass[s]=sum(ww)
                _cdid_db_pair_work(p,p.n_A,1)
                yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                if(missing(yy))_cdid_base_fail("outcome","Pair outcome change overflowed.",504)
                rawZ=_cdid_db_Z(Ref,p,pi,ai);D=L:*rawZ
                if(missing(D))_cdid_base_fail("outcome_design","Calendar-scaled outcome design overflowed.",504)
                if(d) {
                    xmin[s,]=colmin(rawZ);xmax[s,]=colmax(rawZ);firstZ[s,]=rawZ[1,]
                    out.sumZ[s,]=quadcross(ww,rawZ)
                }
                out.sumY[s]=sum(ww:*yy)
                F=_cdid_or_compress(((J(rows(ai),1,1),D,yy):*sqrt(ww)))
            }
            else {
                _cdid_db_pair_work(p,p.n_A,1)
                F=J(0,d+2,0)
            }
        }
        else {
            _cdid_db_pair_work(p,p.n_A,1)
            rows0=_cdid_db_pair_slice(P,p,pi,1,p.n_A)
            ix=_cdid_base_index(p.cohort_to_active[p.a_cohort[rows0[,1]]]:==p.support[s,2])
            ai=rows0[ix,1];out.pair_mass[s]=sum(_cdid_db_w_ptr(pref,p,pi,ai))
            F=J(0,d+2,0);first=1
            for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                hi=min((p.n_A,lo+ctrl.chunk_rows-1))
                rows0=_cdid_db_read(P,p,pi,lo,hi,p.active_cohorts[p.support[s,2]])
                if(rows(rows0)==0)continue
                ai=rows0[,1];ww=_cdid_db_w_ptr(pref,p,pi,ai);yy=rows0[,2];rawZ=_cdid_db_Z(Ref,p,pi,ai);D=L:*rawZ
                if(missing(D))_cdid_base_fail("outcome_design","Calendar-scaled outcome design overflowed.",504)
                if(d) {
                    if(first){xmin[s,]=colmin(rawZ);xmax[s,]=colmax(rawZ);firstZ[s,]=rawZ[1,];}
                    else{xmin[s,]=colmin(xmin[s,]\rawZ);xmax[s,]=colmax(xmax[s,]\rawZ);}
                }
                first=0;out.sumY[s]=out.sumY[s]+sum(ww:*yy)
                if(d)out.sumZ[s,]=out.sumZ[s,]+quadcross(ww,rawZ)
                F=_cdid_or_compress(F\((J(rows(ai),1,1),D,yy):*sqrt(ww)))
            }
        }
        if(out.pair_mass[s]<=0 | out.pair_mass[s]>=.)
            _cdid_base_fail("support","Declared support has nonpositive/nonfinite observed weight mass; none is dropped.",459)
        out.rates[s]=out.pair_mass[s]/groupmass[s];summary[s].F=F
    }
    values=cdid_ipw_rate_values(out.rates);out.rates_c=values[,1];out.rates_active=values[,2]
    nPseudo=2*(d+1)*B
    pseudoX=J(nPseudo,d,0);pseudoY=J(nPseudo,1,0);pseudoW=pseudoY;pseudoB=pseudoY;at=0
    for(b=1;b<=B;b++) {
        pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
        F=J(0,d+2,0);mass=0;constant=J(1,d,1);baseZ=J(1,d,0)
        for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
            s=p.block_support[m]
            F=_cdid_or_compress(F\(summary[s].F/sqrt(out.rates_c[s])))
            mass=mass+out.pair_mass[s]/out.rates_c[s]
            if(m==p.block_ptr[b])baseZ=firstZ[s,]
            if(d)constant=constant:*(xmin[s,]:==baseZ):*(xmax[s,]:==baseZ)
        }
        hqrd(F[,1],Q=.,tau=.,R=.)
        V=hqrdmultq(Q,tau,F[,(2..d+2)],1);coefrows=V[1,]/R[1,1]
        resfactor=V[2..rows(V),]
        if(ctrl.intercept & d)for(k=1;k<=d;k++)if(constant[k]) {
            coefrows[k]=L*baseZ[k];resfactor[,k]=J(rows(resfactor),1,0)
        }
        nn=rows(resfactor);scale=sqrt(nn/mass)
        A=coefrows:+scale:*resfactor;V=coefrows:-scale:*resfactor
        rr=(at+1..at+2*nn)'
        if(d)pseudoX[rr,]=(A[,1..d]\V[,1..d])
        pseudoY[rr]=(A[,d+1]\V[,d+1]);pseudoW[rr]=J(2*nn,1,mass/(2*nn));pseudoB[rr]=J(2*nn,1,b)
        at=at+2*nn
    }
    if(at!=nPseudo)_cdid_base_fail("summary","Unexpected compressed support factor dimension.",3200)
    if(out.outcome_pooling=="separate") {
        out.blocks=cdid_base_or_block(B,1)
        out.alpha=J(B,1,0);out.block_mass=J(B,1,0);out.beta=J(0,1,0)
        out.R=J(0,0,0);out.pivot=J(1,0,0);out.kept=J(1,0,0);out.dropped=J(1,0,0)
        out.alias_map=J(0,0,0);out.alias_constants=J(0,0,0);out.mean_D=J(B,0,0)
        out.beta_score_A=J(0,0,0);out.beta_rate_loading=J(0,S,0)
        at=0;max_k=0
        for(b=1;b<=B;b++) {
            nn=d+1;rr=(at+1..at+2*nn)';at=at+2*nn
            fit=cdid_pooled_wls_fit((d ? pseudoX[rr,] : J(2*nn,0,0)),pseudoY[rr],pseudoW[rr],
                J(2*nn,1,1),J(2*nn,1,1),1,ctrl.intercept,J(0,d,0),J(0,1,0),ctrl.rank_tol,1)
            out.blocks[b].beta=fit.beta;out.blocks[b].R=fit.R;out.blocks[b].pivot=fit.pivot
            out.blocks[b].kept=fit.kept;out.blocks[b].dropped=fit.dropped
            out.blocks[b].alias_map=fit.alias_map;out.blocks[b].alias_constants=fit.alias_constants
            out.blocks[b].mean_D=(ctrl.intercept ? fit.mean_X[,fit.kept] : J(1,rows(fit.beta),0))
            out.alpha[b]=fit.alpha[1];out.block_mass[b]=fit.mass[1]
            if(rows(fit.beta)>max_k)max_k=rows(fit.beta)
        }
        _cdid_db_or_alias(P,p,Ref,out,ctrl.chunk_rows)
        out.alpha_rate=J(rows(p.block_support),1,0);out.n_fit_rows=0;out.objective=0
        for(b=1;b<=B;b++) {
            pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
            kb=rows(out.blocks[b].beta);n_ctrl_b=p.block_ptr[b+1]-p.block_ptr[b]
            force=J(kb,n_ctrl_b,0);sum_sa=0;sum_sb=J(1,kb,0)
            for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
                s=p.block_support[m]
                if(has_csr) {
                    _cdid_db_pair_work(p,p.n_A,1)
                    if(sup_hi[s]<sup_lo[s])continue
                    ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                    yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                    fw=_cdid_db_w_ptr(pref,p,pi,ai)/out.rates_c[s]
                    D=L:*_cdid_db_Z_kept_ptr(pref,p,pi,ai,out.blocks[b].kept)
                    rawres=yy:-out.alpha[b]:-(kb ? D*out.blocks[b].beta : 0)
                    centered=(ctrl.intercept ? D:-out.blocks[b].mean_D : D)
                    coefrows=quadcross(fw:*rawres,centered)
                    if(ctrl.intercept)sum_sa=sum_sa+sum(fw:*rawres)
                    if(kb)sum_sb=sum_sb+coefrows
                    if(out.rates_active[s]) {
                        if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                        if(kb)force[,m-p.block_ptr[b]+1]=force[,m-p.block_ptr[b]+1]-coefrows'/out.rates_c[s]
                    }
                    out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                }
                else {
                    for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                        rows0=_cdid_db_read(P,p,pi,lo,min((p.n_A,lo+ctrl.chunk_rows-1)),p.active_cohorts[p.support[s,2]])
                        if(rows(rows0)==0)continue
                        ai=rows0[,1];fw=_cdid_db_w_ptr(pref,p,pi,ai)/out.rates_c[s]
                        D=L:*_cdid_db_Z_kept_ptr(pref,p,pi,ai,out.blocks[b].kept)
                        rawres=rows0[,2]:-out.alpha[b]:-(kb ? D*out.blocks[b].beta : 0)
                        centered=(ctrl.intercept ? D:-out.blocks[b].mean_D : D)
                        coefrows=quadcross(fw:*rawres,centered)
                        if(ctrl.intercept)sum_sa=sum_sa+sum(fw:*rawres)
                        if(kb)sum_sb=sum_sb+coefrows
                        if(out.rates_active[s]) {
                            if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                            if(kb)force[,m-p.block_ptr[b]+1]=force[,m-p.block_ptr[b]+1]-coefrows'/out.rates_c[s]
                        }
                        out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                    }
                }
            }
            out.blocks[b].score_mean_alpha=(ctrl.intercept ? sum_sa/p.n_F : 0)
            out.blocks[b].score_mean_beta=(kb ? sum_sb/p.n_F : J(1,0,0))
            out.blocks[b].beta_rate_loading=_cdid_nuisance_gram(out.blocks[b].R,out.blocks[b].pivot,force)
        }
        out.link_mu=J(Jl,1,0);out.link_c=J(Jl,max_k,0);out.link_d=J(Jl,1,0);out.link_rate_loading=J(Jl,1,0)
        for(j=1;j<=Jl;j++) {
            b=p.consumer[j,3];s=p.consumer[j,4];pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
            mass=out.pair_mass[s];kb=rows(out.blocks[b].beta)
            den=out.sumY[s]-out.alpha[b]*mass-(kb ? L*(out.sumZ[s,out.blocks[b].kept]*out.blocks[b].beta) : 0)
            out.link_mu[j]=den/mass
            if(kb)out.link_c[j,1..kb]=L:*out.sumZ[s,out.blocks[b].kept]/mass
            out.link_d[j]=mass/(p.n_F*out.rates_c[s])
            if(out.rates_active[s])out.link_rate_loading[j]=-(den-out.link_mu[j]*mass)/(p.n_F*out.rates_c[s]^2*out.link_d[j])
        }
        out.retained_entries=11+5*S+2*B+rows(p.block_support)+S*d+B*(max_k*max_k+5*max_k+1)+
            max_k*rows(p.block_support)+Jl*max_k+3*Jl
        return(out)
    }
    out.blocks=cdid_base_or_block(0,1)
    fit=cdid_pooled_wls_fit(pseudoX,pseudoY,pseudoW,pseudoB,J(nPseudo,1,1),1,
        ctrl.intercept,J(0,d,0),J(0,1,0),ctrl.rank_tol,1)
    out.beta=fit.beta;out.alpha=fit.alpha;out.R=fit.R;out.pivot=fit.pivot
    out.kept=fit.kept;out.dropped=fit.dropped;out.alias_map=fit.alias_map;out.alias_constants=fit.alias_constants
    out.mean_D=(ctrl.intercept ? fit.mean_X[,fit.kept] : J(B,rows(fit.beta),0));out.block_mass=fit.mass
    _cdid_db_or_alias(P,p,Ref,out,ctrl.chunk_rows)
    k=rows(out.beta);out.beta_score_A=J(p.n_A,k,0);force=J(k,S,0)
    out.alpha_rate=J(rows(p.block_support),1,0);out.n_fit_rows=0;out.objective=0
    for(b=1;b<=B;b++) {
        pi=p.block_pair[b];L=p.pair[pi,2]-p.pair[pi,1]
        for(m=p.block_ptr[b];m<p.block_ptr[b+1];m++) {
            s=p.block_support[m]
            if(has_csr) {
                _cdid_db_pair_work(p,p.n_A,1)
                if(sup_hi[s]<sup_lo[s])continue
                ai=sup_a[|sup_lo[s]\sup_hi[s]|]
                yy=P.y[sup_rt[|sup_lo[s]\sup_hi[s]|]]-P.y[sup_rs[|sup_lo[s]\sup_hi[s]|]]
                fw=_cdid_db_w_ptr(pref,p,pi,ai)/out.rates_c[s]
                D=L:*_cdid_db_Z_kept_ptr(pref,p,pi,ai,out.kept)
                rawres=yy:-out.alpha[b]:-D*out.beta
                centered=(ctrl.intercept ? D:-out.mean_D[b,] : D)
                coefrows=quadcross(fw:*rawres,centered)
                if(k)out.beta_score_A[ai,]=out.beta_score_A[ai,]+centered:*(fw:*rawres)
                if(out.rates_active[s]) {
                    if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                    if(k)force[,s]=force[,s]-coefrows'/out.rates_c[s]
                }
                out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
            }
            else {
                for(lo=1;lo<=p.n_A;lo=lo+ctrl.chunk_rows) {
                    rows0=_cdid_db_read(P,p,pi,lo,min((p.n_A,lo+ctrl.chunk_rows-1)),p.active_cohorts[p.support[s,2]])
                    if(rows(rows0)==0)continue
                    ai=rows0[,1];fw=_cdid_db_w_ptr(pref,p,pi,ai)/out.rates_c[s]
                    D=L:*_cdid_db_Z_kept_ptr(pref,p,pi,ai,out.kept)
                    rawres=rows0[,2]:-out.alpha[b]:-D*out.beta
                    centered=(ctrl.intercept ? D:-out.mean_D[b,] : D)
                    coefrows=quadcross(fw:*rawres,centered)
                    if(k)out.beta_score_A[ai,]=out.beta_score_A[ai,]+centered:*(fw:*rawres)
                    if(out.rates_active[s]) {
                        if(ctrl.intercept)out.alpha_rate[m]=out.alpha_rate[m]-sum(fw:*rawres)/out.rates_c[s]
                        if(k)force[,s]=force[,s]-coefrows'/out.rates_c[s]
                    }
                    out.n_fit_rows=out.n_fit_rows+rows(ai);out.objective=out.objective+sum(fw:*rawres:^2)/p.n_F
                }
            }
        }
    }
    out.beta_rate_loading=_cdid_nuisance_gram(out.R,out.pivot,force)
    out.link_mu=J(Jl,1,0);out.link_c=J(Jl,k,0);out.link_d=J(Jl,1,0);out.link_rate_loading=J(Jl,1,0)
    for(j=1;j<=Jl;j++) {
        b=p.consumer[j,3];s=p.consumer[j,4];pi=p.support[s,1];L=p.pair[pi,2]-p.pair[pi,1]
        mass=out.pair_mass[s];den=out.sumY[s]-out.alpha[b]*mass-L*(out.sumZ[s,out.kept]*out.beta)
        out.link_mu[j]=den/mass;out.link_c[j,]=L:*out.sumZ[s,out.kept]/mass
        out.link_d[j]=mass/(p.n_F*out.rates_c[s])
        if(out.rates_active[s])out.link_rate_loading[j]=-(den-out.link_mu[j]*mass)/(p.n_F*out.rates_c[s]^2*out.link_d[j])
    }
    out.retained_entries=11+5*S+2*B+3*k+rows(p.block_support)+S*d+B*k+k*k+
        (k+B)*cols(out.dropped)+p.n_A*k+k*S+Jl*k+3*Jl+cols(out.dropped)
    return(out)
}

real matrix _cdid_db_frame_pair(struct cdid_db_provider scalar f,
    real scalar pi,real colvector rowsF,real colvector groups,real scalar readY)
{
    real colvector aa,ix,keep,rs,rt,dy,cmap,ai
    real scalar s,t
    aa=f.plan.F_to_A[rowsF];ix=_cdid_base_index(aa:>0)
    if(rows(ix)==0)return(J(0,3,0))
    if(rows(groups)==1)keep=_cdid_base_index(f.plan.a_cohort[aa[ix]]:==groups[1])
    else {
        cmap=J(rows(f.plan.cohorts),1,0);cmap[groups]=J(rows(groups),1,1)
        keep=_cdid_base_index(cmap[f.plan.a_cohort[aa[ix]]])
    }
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=aa[ix]
    keep=_cdid_base_index(_cdid_db_keep_mask(f.plan,pi,ai))
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=ai[keep]
    if(rows(f.plan.row_AT)==f.n_A) {
        s=f.plan.pair[pi,1];t=f.plan.pair[pi,2]
        rs=f.plan.row_AT[ai,s];rt=f.plan.row_AT[ai,t]
    }
    else {
        rs=cdid_population_rows(*f.population,rowsF[ix],J(rows(ix),1,f.plan.calendar[f.plan.pair[pi,1]]))
        rt=cdid_population_rows(*f.population,rowsF[ix],J(rows(ix),1,f.plan.calendar[f.plan.pair[pi,2]]))
    }
    keep=_cdid_base_index((rs:>0):&(rt:>0))
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=ai[keep];rs=rs[keep];rt=rt[keep]
    keep=_cdid_base_index(((*f.population).y[rs]:<.):&((*f.population).y[rt]:<.))
    if(rows(keep)==0)return(J(0,3,0))
    ix=ix[keep];ai=ai[keep];rs=rs[keep];rt=rt[keep];dy=J(rows(ix),1,0)
    if(readY)dy=(*f.population).y[rt]-(*f.population).y[rs]
    if(missing(dy))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
    return((ix,ai,dy))
}

real colvector _cdid_db_rate_if_rr(struct cdid_db_provider scalar f,
    real scalar sid,real colvector rowsF,real matrix rr)
{
    real colvector out,aa,ix,keep
    real scalar cohort,pi,mass
    out=J(rows(rowsF),1,0);aa=f.plan.F_to_A[rowsF]
    cohort=f.plan.active_cohorts[f.plan.support[sid,2]];pi=f.plan.support[sid,1]
    mass=f.group_mass[sid]
    ix=_cdid_base_index(aa:>0)
    if(rows(ix)) {
        keep=_cdid_base_index(f.plan.a_cohort[aa[ix]]:==cohort)
        if(rows(keep)) {
            ix=ix[keep]
            ix=ix[_cdid_base_index(_cdid_db_keep_mask(f.plan,pi,aa[ix]))]
            if(rows(ix)) out[ix]=(-f.n_F:*_cdid_db_w_ptr(f.reference,f.plan,pi,aa[ix]):*f.outcome.rates[sid]):/mass
        }
    }
    if(rows(rr))out[rr[,1]]=out[rr[,1]]:+f.n_F:*_cdid_db_w_ptr(f.reference,f.plan,pi,rr[,2])/mass
    return(out)
}

real colvector _cdid_db_rate_if(struct cdid_db_provider scalar f,
    real scalar sid,real colvector rowsF)
{
    real matrix rr
    real scalar cohort,pi
    cohort=f.plan.active_cohorts[f.plan.support[sid,2]];pi=f.plan.support[sid,1]
    rr=_cdid_db_frame_pair(f,pi,rowsF,cohort,0)
    return(_cdid_db_rate_if_rr(f,sid,rowsF,rr))
}

real matrix _cdid_db_beta_if(struct cdid_db_provider scalar f,real colvector rowsF)
{
    real matrix score,out
    real colvector aa,ix,up
    real scalar k,s
    k=rows(f.outcome.beta);score=J(rows(rowsF),k,0)
    if(k==0)return(score)
    aa=f.plan.F_to_A[rowsF];ix=_cdid_base_index(aa:>0)
    if(rows(ix))score[ix,]=f.n_F:*f.outcome.beta_score_A[aa[ix],]
    out=_cdid_nuisance_gram(f.outcome.R,f.outcome.pivot,score')'
    for(s=1;s<=f.plan.n_supports;s++) {
        if(!any(f.outcome.beta_rate_loading[,s]))continue
        up=_cdid_db_rate_if(f,s,rowsF)
        out=out+up*f.outcome.beta_rate_loading[,s]'
    }
    return(out)
}

real matrix _cdid_db_blk_beta_if(struct cdid_db_provider scalar f,
    real scalar b,real colvector rowsF)
{
    real scalar n,kb,pi,L,m,s,c_s
    real matrix score,out,rr_all,rr,Z_all,D,centered
    real colvector ai,w_all,fw,res,up,b_co,ix_s
    n=rows(rowsF);kb=rows(f.outcome.blocks[b].beta)
    if(kb==0)return(J(n,0,0))
    score=J(n,kb,0)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_db_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all)) {
        Z_all=L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,rr_all[,2],f.outcome.blocks[b].kept)
        w_all=_cdid_db_w_ptr(f.reference,f.plan,pi,rr_all[,2])
        for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
            s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
            ix_s=_cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s)
            if(rows(ix_s)) {
                rr=rr_all[ix_s,];D=Z_all[ix_s,]
                res=rr[,3]:-f.outcome.alpha[b]:-D*f.outcome.blocks[b].beta
                fw=f.n_F:*w_all[ix_s]/f.outcome.rates_c[s]
                centered=(f.outcome.intercept ? D:-f.outcome.blocks[b].mean_D : D)
                score[rr[,1],]=score[rr[,1],]+centered:*(fw:*res)
            }
        }
    }
    score=score:-(f.n_F:*f.outcome.blocks[b].score_mean_beta)
    out=_cdid_nuisance_gram(f.outcome.blocks[b].R,f.outcome.blocks[b].pivot,score')'
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        if(!any(f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]))continue
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]'
    }
    return(out)
}

real colvector _cdid_db_alpha_if(struct cdid_db_provider scalar f,
    real scalar b,real colvector rowsF,real matrix betaIF)
{
    real colvector out,up,w_all,res,fw,beta,b_co,ix_s
    real rowvector kept,mean_D
    real scalar m,s,pi,L,kb,sep,c_s
    real matrix rr_all,rr,Z_all,D
    out=J(rows(rowsF),1,0)
    if(!f.outcome.intercept)return(out)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    mean_D=(sep ? f.outcome.blocks[b].mean_D : f.outcome.mean_D[b,])
    kb=rows(beta)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_db_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all)) {
        Z_all=(kb ? L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,rr_all[,2],kept) : J(rows(rr_all),0,0))
        w_all=_cdid_db_w_ptr(f.reference,f.plan,pi,rr_all[,2])
    }
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        if(rows(rr)) {
            D=(kb ? Z_all[ix_s,] : J(rows(rr),0,0))
            res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0)
            fw=f.n_F:*w_all[ix_s]/f.outcome.rates_c[s]/f.outcome.block_mass[b]
            out[rr[,1]]=out[rr[,1]]+fw:*res
        }
        if(f.outcome.alpha_rate[m]!=0) {
            up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
            out=out+up*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(sep)out=out:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
    if(kb)out=out-betaIF*mean_D'
    return(out)
}
void _cdid_db_blk_beta_alpha_if(struct cdid_db_provider scalar f,
    real scalar b,real colvector rowsF,real matrix betaIF,real colvector alphaIF)
{
    real scalar n,kb,pi,L,m,s,c_s,do_alpha
    real matrix score,rr_all,rr,Z_all,D,centered
    real colvector w_all,fw,res_all,res,up,b_co,ix_s
    n=rows(rowsF);kb=rows(f.outcome.blocks[b].beta);do_alpha=(f.outcome.intercept!=0)
    betaIF=J(n,kb,0);alphaIF=J(n,1,0)
    if(kb==0 & !do_alpha) return
    score=J(n,kb,0)
    pi=f.plan.block_pair[b];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
    b_co=f.plan.active_cohorts[f.plan.support[f.plan.block_support[f.plan.block_ptr[b]..f.plan.block_ptr[b+1]-1],2]]
    rr_all=_cdid_db_frame_pair(f,pi,rowsF,b_co,1)
    if(rows(rr_all)) {
        w_all=_cdid_db_w_ptr(f.reference,f.plan,pi,rr_all[,2])
        if(kb) {
            Z_all=L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,rr_all[,2],f.outcome.blocks[b].kept)
            res_all=rr_all[,3]:-f.outcome.alpha[b]:-Z_all*f.outcome.blocks[b].beta
        }
        else res_all=rr_all[,3]:-f.outcome.alpha[b]
    }
    for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
        s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
        ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
        rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
        if(rows(rr)) {
            res=res_all[ix_s]
            fw=f.n_F:*w_all[ix_s]/f.outcome.rates_c[s]
            if(kb) {
                D=Z_all[ix_s,]
                centered=(f.outcome.intercept ? D:-f.outcome.blocks[b].mean_D : D)
                score[rr[,1],]=score[rr[,1],]+centered:*(fw:*res)
            }
            if(do_alpha) alphaIF[rr[,1]]=alphaIF[rr[,1]]+(fw/f.outcome.block_mass[b]):*res
        }
        if(do_alpha & f.outcome.alpha_rate[m]!=0) {
            up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
            alphaIF=alphaIF+up*(f.outcome.alpha_rate[m]/f.outcome.block_mass[b])
        }
    }
    if(kb) {
        score=score:-(f.n_F:*f.outcome.blocks[b].score_mean_beta)
        betaIF=_cdid_nuisance_gram(f.outcome.blocks[b].R,f.outcome.blocks[b].pivot,score')'
        for(m=f.plan.block_ptr[b];m<f.plan.block_ptr[b+1];m++) {
            if(!any(f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]))continue
            s=f.plan.block_support[m];c_s=f.plan.active_cohorts[f.plan.support[s,2]]
            ix_s=(rows(rr_all) ? _cdid_base_index(f.plan.a_cohort[rr_all[,2]]:==c_s) : J(0,1,0))
            rr=(rows(ix_s) ? rr_all[ix_s,] : J(0,3,0))
            up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
            betaIF=betaIF+up*f.outcome.blocks[b].beta_rate_loading[,m-f.plan.block_ptr[b]+1]'
        }
    }
    if(do_alpha) {
        alphaIF=alphaIF:-(f.n_F*f.outcome.blocks[b].score_mean_alpha/f.outcome.block_mass[b])
        if(kb)alphaIF=alphaIF-betaIF*f.outcome.blocks[b].mean_D'
    }
}

real colvector _cdid_db_or_if(struct cdid_db_provider scalar f,
    real scalar ell,real colvector rowsF,real matrix betaIF,| real colvector alphaIF_in)
{
    real scalar b,s,pi,L,kb,sep
    real colvector out,ai,res,up,beta,alphaIF
    real rowvector kept
    real matrix rr,D
    b=f.plan.consumer[ell,3];s=f.plan.consumer[ell,4];pi=f.plan.support[s,1]
    L=f.plan.pair[pi,2]-f.plan.pair[pi,1];out=J(rows(rowsF),1,0)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    kb=rows(beta)
    rr=_cdid_db_frame_pair(f,pi,rowsF,f.plan.active_cohorts[f.plan.support[s,2]],1)
    if(rows(rr)) {
        ai=rr[,2];D=(kb ? L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,ai,kept) : J(rows(ai),0,0))
        res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0):-f.outcome.link_mu[ell]
        out[rr[,1]]=f.n_F:*_cdid_db_w_ptr(f.reference,f.plan,pi,ai)/f.outcome.pair_mass[s]:*res
    }
    alphaIF=(args()>=5 & rows(alphaIF_in)==rows(rowsF) ? alphaIF_in : _cdid_db_alpha_if(f,b,rowsF,betaIF))
    out=out-alphaIF
    if(f.outcome.link_rate_loading[ell]!=0) {
        up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.link_rate_loading[ell]
    }
    if(kb)out=out-betaIF*f.outcome.link_c[ell,1..kb]'
    return(out)
}
real colvector _cdid_db_or_if_nobeta(struct cdid_db_provider scalar f,
    real scalar ell,real colvector rowsF,real colvector alphaIF)
{
    real scalar b,s,pi,L,kb,sep
    real colvector out,ai,res,up,beta
    real rowvector kept
    real matrix rr,D
    b=f.plan.consumer[ell,3];s=f.plan.consumer[ell,4];pi=f.plan.support[s,1]
    L=f.plan.pair[pi,2]-f.plan.pair[pi,1];out=J(rows(rowsF),1,0)
    sep=(f.outcome.outcome_pooling=="separate")
    kept=(sep ? f.outcome.blocks[b].kept : f.outcome.kept)
    beta=(sep ? f.outcome.blocks[b].beta : f.outcome.beta)
    kb=rows(beta)
    rr=_cdid_db_frame_pair(f,pi,rowsF,f.plan.active_cohorts[f.plan.support[s,2]],1)
    if(rows(rr)) {
        ai=rr[,2];D=(kb ? L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,ai,kept) : J(rows(ai),0,0))
        res=rr[,3]:-f.outcome.alpha[b]:-(kb ? D*beta : 0):-f.outcome.link_mu[ell]
        out[rr[,1]]=f.n_F:*_cdid_db_w_ptr(f.reference,f.plan,pi,ai)/f.outcome.pair_mass[s]:*res
    }
    out=out-alphaIF
    if(f.outcome.link_rate_loading[ell]!=0) {
        up=_cdid_db_rate_if_rr(f,s,rowsF,rr)
        out=out+up*f.outcome.link_rate_loading[ell]
    }
    return(out)
}
/* scoped factory: all risk rows and every used-cohort A mass come from
   original A/F positions. No fake all-A weight vector, zero fill, or cropped frame. */
struct cdid_ipw_nuisance scalar _cdid_db_nuisance_from_risk(struct cdid_reference_inputs scalar S,
 struct cdid_db_plan scalar p,real scalar pi,real colvector co,real colvector risk,
 pointer(real colvector) colvector coh_allA,struct cdid_baseline_controls scalar ctrl)
{
 struct cdid_ipw_nuisance scalar out
 struct cdid_logit_result scalar fitted
 pointer(struct cdid_reference_inputs scalar) scalar pref
 real scalar j,n
 real colvector ix,allA,groupmass,D,weights,residual
 real matrix rv,lv
 pref=&S;n=rows(risk)
 if(n<1)_cdid_base_fail("support","Requested dated risk support has no actual pair observers.",459)
 out.n_A=p.n_A;out.n_F=p.n_F;out.n_cohorts=rows(co);out.n_pair_all=n;out.n_pair_risk=n
 out.in_risk=(1::n);out.risk_a=risk
 out.risk_pos_a=J(p.n_A,1,0);out.risk_pos_a[risk]=(1::n)
 out.risk_weight=_cdid_db_w_ptr(pref,p,pi,risk)
 out.X=_cdid_db_X2_ptr(pref,p,pi,risk);out.p=cols(out.X);out.risk_group=J(n,1,0)
 out.pair_count=out.pair_mass=groupmass=J(rows(co),1,0)
 for(j=1;j<=rows(co);j++) {
  p.work_ops=_cdid_base_work(n+p.n_A,p.work_limit,p.work_ops)
  ix=_cdid_base_index(p.a_cohort[risk]:==co[j]);out.risk_group[ix]=J(rows(ix),1,j)
  out.pair_count[j]=rows(ix);out.pair_mass[j]=sum(out.risk_weight[ix])
  allA=J(0,1,0)
  if(rows(coh_allA)>=co[j]) {
   if(coh_allA[co[j]]!=NULL) allA=*coh_allA[co[j]]
   else allA=_cdid_base_index(p.a_cohort:==co[j])
  }
  else allA=_cdid_base_index(p.a_cohort:==co[j])
  groupmass[j]=sum(_cdid_db_w_ptr(pref,p,pi,_cdid_db_keep(p,pi,allA)))
 }
 if(missing(out.X)|missing(out.risk_weight)|any(out.risk_weight:<=0)|any(groupmass:<=0)|any(out.pair_mass:<=0))
  _cdid_base_fail("support","Finite used risk design and positive used A/pair masses required.",459)
 out.rates=out.pair_mass:/groupmass;rv=cdid_ipw_rate_values(out.rates)
 out.rates_c=rv[,1];out.rates_active=rv[,2];out.q=groupmass/p.n_A
 D=(out.risk_group:==1);out.fit_weight=out.risk_weight:/out.rates_c[out.risk_group]
 fitted=cdid_logit_fit(out.X,D,out.fit_weight,ctrl.glm_epsilon,ctrl.glm_maxit,1)
 out.beta=fitted.beta;out.iterations=fitted.iterations;out.deviance=fitted.deviance;out.score_max=fitted.score_max
 out.epsilon=ctrl.glm_epsilon;out.maxit=ctrl.glm_maxit
 lv=_cdid_ipw_logit_values(out.X*out.beta);out.raw_probability=lv[,1];out.odds=lv[,2];out.odds_active=lv[,3]
 ix=_cdid_base_index(D);out.sum_a=sum(out.fit_weight[ix]);ix=_cdid_base_index(!D)
 out.sum_b=sum(out.fit_weight[ix]:*out.odds[ix]);out.mean_a=out.sum_a/p.n_A;out.mean_b=out.sum_b/p.n_A
 if(out.mean_a<=0|out.mean_b<=0|out.mean_a>=.|out.mean_b>=.)_cdid_base_fail("denominator","Nonpositive/nonfinite dated ratio normalizer.",504)
 out.row_score=fitted.row_score;out.information=fitted.information;out.H=out.information/p.n_A
 out.H_chol=cholesky(out.H)
 if(missing(out.H_chol)|any(diagonal(out.H_chol):<=0))_cdid_base_fail("rank","Dated PS information is not positive definite.",459)
 residual=D:-out.raw_probability;out.A_beta_r=J(out.p,rows(co),0)
 for(j=1;j<=rows(co);j++)if(out.rates_active[j]) {
  ix=_cdid_base_index(out.risk_group:==j)
  weights=out.risk_weight[ix]:*residual[ix]:/(out.rates_c[j]^2)
  out.A_beta_r[,j]=-quadcross(out.X[ix,],weights)/p.n_A
 }
 return(out)
}
struct cdid_ipw_nuisance scalar _cdid_db_nuisance(struct cdid_population scalar P,
 struct cdid_reference_inputs scalar S,struct cdid_db_plan scalar p,
 real scalar ell,struct cdid_baseline_controls scalar ctrl)
{
 pointer(real colvector) colvector empty_coh
 real scalar pi,b
 real colvector co
 real matrix rr
 pi=p.consumer[ell,2];b=p.consumer[ell,3]
 co=p.consumer[ell,1]\p.active_cohorts[p.support[p.block_support[p.block_ptr[b]..p.block_ptr[b+1]-1],2]]
 rr=_cdid_db_read(P,p,pi,1,p.n_A,co)
 empty_coh=J(0,1,NULL)
 return(_cdid_db_nuisance_from_risk(S,p,pi,co,rr[,1],empty_coh,ctrl))
}
struct cdid_db_provider scalar cdid_db_prepare(pointer(struct cdid_population scalar) scalar population,
 pointer(struct cdid_reference_inputs scalar) scalar reference,
 struct cdid_baseline_plan_state scalar original,struct cdid_reference_catalog scalar catalog,
 struct cdid_baseline_controls scalar ctrl)
{
 struct cdid_db_provider scalar out
 struct cdid_population scalar P
 struct cdid_reference_inputs scalar S
 struct cdid_ipw_nuisance scalar nu
 struct cdid_ipw_link scalar link
 struct cdid_db_nu_entry colvector nu_cache
 transmorphic scalar ps_buckets
 pointer(real colvector) colvector coh_allA
 pointer(struct cdid_reference_inputs scalar) scalar pref
 string scalar key_ps
 real scalar ell,pi,b,L,j,px,pz,used,registry,r_rec,x_date,w_anc,n_unique_ps,hit,c_cand,rep_ell
 real colvector co,risk,dy,ix,fw,beta,c_list,cmap
 real rowvector kept
 real matrix rr,D
 if(population==NULL|reference==NULL)_cdid_base_fail("owner","Original population and reference objects required.",498)
 P=*population;S=*reference;cdid_ref_activate(S);pref=&S
 if(S.population!=population|S.route!="baseline"|S.method!=ctrl.method|P.N!=original.n_F|P.n_A!=original.n_A)
  _cdid_base_fail("identity","Original population/reference/route/method binding disagrees.",498)
 if(ctrl.intercept!=S.roles[1].intercept)_cdid_base_fail("design","Outcome intercept metadata disagrees.",498)
 /* All live scientific owners remain part of the charged numeric workspace. */
 registry=length(catalog.consumers)+length(catalog.blocks)+length(catalog.representative)+
  length(catalog.consumer_recipe)+length(catalog.consumer_block)+length(catalog.K_ptr)+
  length(catalog.K_index)+length(catalog.block_K_ptr)+length(catalog.block_K_index)+
  length(catalog.raw_consumer)+length(catalog.post_consumer)+length(catalog.aux_consumer)+
  length(catalog.raw_supported)+length(catalog.aux_supported)+5
 used=_cdid_base_cells(1,original.workspace_entries,ctrl.max_entries,0)
 used=_cdid_base_cells(1,S.retained_entries,ctrl.max_entries,used)
 used=_cdid_base_cells(1,registry,ctrl.max_entries,used)
 out.plan=_cdid_db_plan_make(original,catalog,ctrl.max_entries,used,ctrl.workcap)
 used=_cdid_base_cells(1,out.plan.workspace_entries,ctrl.max_entries,used)
 out.catalog=catalog;out.population=population;out.reference=reference
 out.method=ctrl.method;out.n_F=P.N;out.n_A=P.n_A;out.core_to_frame=(P.n_A ? P.N/P.n_A : .)
 out.normalization_count=P.n_A;out.weight_row_domain="F"
 out.weight_rows=P.N;out.sample_rank_upper=out.weight_rows-1
 out.chunk_rows=ctrl.chunk_rows;out.workspace_limit=ctrl.max_entries
 out.n_ps_fits=out.n_ps_reuses=out.n_or_fits=0
 out.status=(catalog.n_consumers ? "prepared" : "no_supported_consumers")
 px=(ctrl.method=="reg" ? 0 : cols(S.designs[S.roles[2].slot].value))
 pz=(ctrl.method=="ipw" ? 0 : cols(S.designs[S.roles[1].slot].value)-ctrl.intercept)
 used=_cdid_base_cells(12*catalog.n_consumers,px+2*rows(original.cohorts)+16,ctrl.max_entries,used)
 used=_cdid_base_cells(10*P.n_A,px+pz+4,ctrl.max_entries,used)
 used=_cdid_base_cells(8*px,px,ctrl.max_entries,used)
 out.workspace_entries=used;out.workspace_peak_entries=used;out.work_ops=out.plan.work_ops;out.points=J(catalog.n_consumers,1,.)
 out.point=cdid_base_point(catalog.n_consumers,1);out.propensity=cdid_base_ps(catalog.n_consumers,1)
 out.ps_cache=cdid_db_ps_cache(0,1);out.ps_class_id=J(catalog.n_consumers,1,0)
 out.group_mass=J(out.plan.n_supports,1,0);out.outcome.p=0;out.outcome.beta=J(0,1,0);out.outcome.alpha=J(0,1,0)
 out.plan.pair_ctx=J(out.plan.n_pairs,1,0);out.plan.elig_A=J(P.n_A,0,.)
 if(rows(S.elig_index)) {
  if(cols(S.recipes)<3 | rows(S.elig_excluded)!=P.N | cols(S.elig_excluded)!=rows(S.elig_index))
   _cdid_base_fail("eligibility","Reference eligibility contexts and design blocks disagree.",498)
  for(pi=1;pi<=out.plan.n_pairs;pi++) {
   r_rec=S.recipes[out.plan.pair_recipe[pi],3]
   if(r_rec>0) {
    ix=_cdid_base_index(S.elig_index:==r_rec)
    if(rows(ix)!=1) _cdid_base_fail("eligibility","Unknown reference eligibility context.",498)
    out.plan.pair_ctx[pi]=ix[1]
   }
  }
  out.plan.elig_A=(P.n_A ? !S.elig_excluded[P.A_to_F,] : J(0,rows(S.elig_index),.))
 }
 if(catalog.n_consumers==0)return(out)
 _cdid_db_build_caches(P,S,out.plan,original.row_AT)
 if(ctrl.method!="ipw") {
  out.outcome=_cdid_db_or_fit(P,out.plan,S,ctrl,used)
  out.n_or_fits=(out.outcome.outcome_pooling=="separate" ? out.plan.n_blocks : 1)
  out.group_mass=out.outcome.repeat_mass
  out.workspace_entries=_cdid_base_cells(1,out.outcome.retained_entries,ctrl.max_entries,used)
  out.workspace_peak_entries=max((out.workspace_entries,out.outcome.peak_entries))
  out.work_ops=out.plan.work_ops
  if(ctrl.method=="reg") {
   out.points=out.outcome.link_mu
   for(ell=1;ell<=catalog.n_consumers;ell++){out.point[ell].estimate=out.points[ell];out.point[ell].or_id=ell;out.point[ell].ps_id=0;}
   out.propensity=cdid_base_ps(0,1);return(out)
  }
 }
 /* Exact qualified consumer keys own PS models; deduplicate identical recipes. */
 coh_allA=J(rows(out.plan.cohorts),1,NULL)
 for(j=1;j<=rows(out.plan.active_cohorts);j++) {
  r_rec=out.plan.active_cohorts[j]
  ix=_cdid_base_index(out.plan.a_cohort:==r_rec)
  coh_allA[r_rec]=&(ix:+0)
 }
 nu_cache=cdid_db_nu_entry(catalog.n_consumers,1)
 ps_buckets=asarray_create("string",1);n_unique_ps=0
 for(ell=1;ell<=catalog.n_consumers;ell++) {
  pi=out.plan.consumer[ell,2];b=out.plan.consumer[ell,3];L=out.plan.pair[pi,2]-out.plan.pair[pi,1]
  co=out.plan.consumer[ell,1]\out.plan.active_cohorts[out.plan.support[out.plan.block_support[out.plan.block_ptr[b]..out.plan.block_ptr[b+1]-1],2]]
  rr=_cdid_db_read(P,out.plan,pi,1,P.n_A,co)
  risk=rr[,1];dy=rr[,2]
  if(rows(risk)<1)_cdid_base_fail("support","Requested dated risk support has no actual pair observers.",459)
  r_rec=out.plan.pair_recipe[pi]
  x_date=(S.roles[2].dated ? S.recipes[r_rec,1] : 0);w_anc=S.recipes[r_rec,2]
  key_ps=strofreal(x_date)+":"+strofreal(w_anc)+":"+invtokens(strofreal(co)')
  if(cols(S.recipes)>=3) key_ps=key_ps+":e"+strofreal(S.recipes[r_rec,3])
  hit=0
  if(asarray_contains(ps_buckets,key_ps)) {
   c_list=asarray(ps_buckets,key_ps)
   for(j=1;j<=rows(c_list);j++) {
    c_cand=c_list[j]
    if(pi==nu_cache[c_cand].pi_rep) {hit=c_cand;break;}
    if(rows(risk)==rows(nu_cache[c_cand].nu.risk_a)) {
     if(all(risk:==nu_cache[c_cand].nu.risk_a)) {hit=c_cand;break;}
    }
   }
  }
  if(hit==0) {
   nu=_cdid_db_nuisance_from_risk(S,out.plan,pi,co,risk,coh_allA,ctrl)
   n_unique_ps++;hit=n_unique_ps
   nu_cache[hit].rep_ell=ell;nu_cache[hit].pi_rep=pi;nu_cache[hit].nu=nu
   nu_cache[hit].t_ix=_cdid_base_index(nu.risk_group:==1)
   nu_cache[hit].c_ix=_cdid_base_index(nu.risk_group:>1)
   nu_cache[hit].fw_ctrl=nu.fit_weight[nu_cache[hit].c_ix]:*nu.odds[nu_cache[hit].c_ix]
   cmap=J(rows(out.plan.cohorts),1,0);cmap[co]=(1::rows(co))
   nu_cache[hit].cmap=cmap
   if(asarray_contains(ps_buckets,key_ps)) asarray(ps_buckets,key_ps,(asarray(ps_buckets,key_ps)\hit))
   else asarray(ps_buckets,key_ps,hit)
   out.propensity[ell].n_A=P.n_A;out.propensity[ell].n_F=P.N;out.propensity[ell].n_cohorts=nu.n_cohorts
   out.propensity[ell].cohorts=co;out.propensity[ell].p=nu.p;out.propensity[ell].q=nu.q
   out.propensity[ell].pair_count=nu.pair_count;out.propensity[ell].pair_mass=nu.pair_mass
   out.propensity[ell].rates=nu.rates;out.propensity[ell].rates_c=nu.rates_c;out.propensity[ell].rates_active=nu.rates_active
   out.propensity[ell].beta=nu.beta;out.propensity[ell].A_beta_r=nu.A_beta_r
   out.propensity[ell].mean_a=nu.mean_a;out.propensity[ell].mean_b=nu.mean_b
   out.propensity[ell].iterations=nu.iterations;out.propensity[ell].deviance=nu.deviance;out.propensity[ell].score_max=nu.score_max
  }
  else {
   nu=nu_cache[hit].nu;rep_ell=nu_cache[hit].rep_ell
   out.plan.work_ops=_cdid_base_work(rows(co)*(nu.n_pair_risk+P.n_A),out.plan.work_limit,out.plan.work_ops)
   out.propensity[ell]=out.propensity[rep_ell]
   out.n_ps_reuses=out.n_ps_reuses+1
  }
  out.n_ps_fits=out.n_ps_fits+1
  _cdid_db_pair_work(out.plan,P.n_A,rows(co))
  out.ps_class_id[ell]=hit
  out.point[ell].ps_id=ell;out.point[ell].or_id=(ctrl.method=="dr" ? ell : 0);out.point[ell].d_bar_diff=J(0,1,0)
  if(ctrl.method=="dr") {
   kept=(out.outcome.outcome_pooling=="separate" ? out.outcome.blocks[b].kept : out.outcome.kept)
   beta=(out.outcome.outcome_pooling=="separate" ? out.outcome.blocks[b].beta : out.outcome.beta)
   if(rows(beta)) {
    D=L:*_cdid_db_Z_kept_ptr(pref,out.plan,pi,nu.risk_a,kept)
    dy=dy:-out.outcome.alpha[b]:-D*beta
    out.point[ell].d_bar_diff=quadcross(nu.fit_weight[nu_cache[hit].t_ix],D[nu_cache[hit].t_ix,])'/nu.sum_a-quadcross(nu_cache[hit].fw_ctrl,D[nu_cache[hit].c_ix,])'/nu.sum_b
   }
   else dy=dy:-out.outcome.alpha[b]
  }
  link=cdid_ipw_link_outcome(nu,dy);out.point[ell].residual_ipw=link;out.points[ell]=link.estimate;out.point[ell].estimate=link.estimate
 }
 out.ps_cache=cdid_db_ps_cache(n_unique_ps,1)
 for(j=1;j<=n_unique_ps;j++) {
  out.ps_cache[j].rep_ell=nu_cache[j].rep_ell;out.ps_cache[j].pi=nu_cache[j].pi_rep
  out.ps_cache[j].risk_a=nu_cache[j].nu.risk_a;out.ps_cache[j].risk_pos_a=nu_cache[j].nu.risk_pos_a
  out.ps_cache[j].raw_probability=nu_cache[j].nu.raw_probability;out.ps_cache[j].odds=nu_cache[j].nu.odds
  out.ps_cache[j].cmap=nu_cache[j].cmap
 }
 out.work_ops=out.plan.work_ops
 return(out)
}
real colvector _cdid_db_ipw_if(struct cdid_db_provider scalar f,
 real scalar ell,real colvector rowsF,real matrix betaIF)
{
 struct cdid_base_ps scalar ps
 struct cdid_ipw_link scalar lk
 real scalar c,pi,b,L,factor
 real colvector out,aa,ix,cj,keep,cmap,zload,eta,dy,ww,beta,fitw,t_ix,c_ix
 real rowvector kept
 real matrix rr,lv,X_rr
 c=f.point[ell].ps_id;ps=f.propensity[c];lk=f.point[ell].residual_ipw
 pi=f.plan.consumer[ell,2];b=f.plan.consumer[ell,3];L=f.plan.pair[pi,2]-f.plan.pair[pi,1]
 aa=f.plan.F_to_A[rowsF];out=J(rows(rowsF),1,0);factor=f.n_F/f.n_A
 cmap=J(rows(f.plan.cohorts),1,0);cmap[ps.cohorts]=(1::rows(ps.cohorts))
 ix=_cdid_base_index(aa:>0)
 if(rows(ix)) {
  cj=cmap[f.plan.a_cohort[aa[ix]]]
  keep=_cdid_base_index(cj:>0)
  if(rows(keep)) {
   ix=ix[keep];cj=cj[keep]
   keep=_cdid_base_index(_cdid_db_keep_mask(f.plan,pi,aa[ix]))
   ix=(rows(keep) ? ix[keep] : J(0,1,.));cj=(rows(keep) ? cj[keep] : J(0,1,.))
   if(rows(ix)) out[ix]=(factor:*_cdid_db_w_ptr(f.reference,f.plan,pi,aa[ix])):*lk.rate_coef_s0[cj]
  }
 }
 rr=_cdid_db_frame_pair(f,pi,rowsF,ps.cohorts,1)
 if(rows(rr)) {
  X_rr=_cdid_db_X2_ptr(f.reference,f.plan,pi,rr[,2])
  eta=X_rr*ps.beta;lv=_cdid_ipw_logit_values(eta)
  zload=X_rr*lk.loading_beta;dy=rr[,3]
  ww=_cdid_db_w_ptr(f.reference,f.plan,pi,rr[,2])
  if(f.method=="dr") {
   kept=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].kept : f.outcome.kept)
   beta=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].beta : f.outcome.beta)
   dy=dy:-f.outcome.alpha[b]:-(rows(beta) ? L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,rr[,2],kept)*beta : 0)
  }
  cj=cmap[f.plan.a_cohort[rr[,2]]]
  fitw=ww:/ps.rates_c[cj]
  out[rr[,1]]=factor:*(ww:*lk.rate_coef_s1[cj]+fitw:*((cj:==1):-lv[,1]):*zload)
  t_ix=_cdid_base_index(cj:==1)
  if(rows(t_ix))out[rr[t_ix,1]]=out[rr[t_ix,1]]+factor:*fitw[t_ix]:/ps.mean_a:*(dy[t_ix]:-lk.treated_mean)
  c_ix=_cdid_base_index(cj:>1)
  if(rows(c_ix))out[rr[c_ix,1]]=out[rr[c_ix,1]]-factor:*fitw[c_ix]:*lv[c_ix,2]:/ps.mean_b:*(dy[c_ix]:-lk.control_mean)
 }
 if(f.method=="dr" & rows(f.point[ell].d_bar_diff))out=out-betaIF*f.point[ell].d_bar_diff
 return(out)
}
real matrix _cdid_db_ps_class_if(struct cdid_db_provider scalar f,
 real scalar c_ps,real colvector rowsF,real colvector aa,real colvector ix_pos_a,
 real colvector ells)
{
 struct cdid_base_ps scalar ps
 real scalar m_c,B,factor,p_trt,n_coh,q,ell,pi_rep,pi,b,L,s_t,t_t,g_idx
 real colvector cmap,cj_a,keep,ix_a,ai_a,pos_v,ix_r,ai_r,pos_r,cj_r,w_r,fitw,raw_p,odds_r
 real colvector t_ix,c_ix,dy_one,beta,grp_key,grp_ord,sub_q
 real rowvector t_mean,c_mean,kept
 real matrix out_c,Rate0,Rate1,Lambda,X_rr,Z_load,DY,OnPair,grp_info,rr_fb
 B=rows(rowsF);m_c=rows(ells);out_c=J(B,m_c,0)
 if(m_c==0 | B==0) return(out_c)
 if(rows(f.ps_cache)<c_ps | rows(f.plan.row_AT)!=f.n_A) {
  for(q=1;q<=m_c;q++) out_c[,q]=_cdid_db_ipw_if(f,ells[q],rowsF,J(B,0,0))
  return(out_c)
 }
 ps=f.propensity[f.ps_cache[c_ps].rep_ell];pi_rep=f.ps_cache[c_ps].pi;cmap=f.ps_cache[c_ps].cmap
 n_coh=ps.n_cohorts;p_trt=ps.p;factor=f.n_F/f.n_A
 if(rows(ix_pos_a)==0) return(out_c)
 Rate0=J(n_coh,m_c,0);Rate1=J(n_coh,m_c,0);Lambda=J(p_trt,m_c,0)
 t_mean=J(1,m_c,0);c_mean=J(1,m_c,0)
 for(q=1;q<=m_c;q++) {
  ell=ells[q]
  Rate0[,q]=f.point[ell].residual_ipw.rate_coef_s0
  Rate1[,q]=f.point[ell].residual_ipw.rate_coef_s1
  Lambda[,q]=f.point[ell].residual_ipw.loading_beta
  t_mean[q]=f.point[ell].residual_ipw.treated_mean
  c_mean[q]=f.point[ell].residual_ipw.control_mean
 }
 ai_a=aa[ix_pos_a];cj_a=cmap[f.plan.a_cohort[ai_a]]
 keep=_cdid_base_index(cj_a:>0)
 if(rows(keep)) {
  ix_a=ix_pos_a[keep];cj_a=cj_a[keep]
  keep=_cdid_base_index(_cdid_db_keep_mask(f.plan,pi_rep,aa[ix_a]))
  if(rows(keep)) {
   ix_a=ix_a[keep];cj_a=cj_a[keep]
   out_c[ix_a,]=(factor:*_cdid_db_w_ptr(f.reference,f.plan,pi_rep,aa[ix_a])):*Rate0[cj_a,]
  }
 }
 if(rows(f.ps_cache[c_ps].risk_pos_a)==f.n_A) {
  pos_v=f.ps_cache[c_ps].risk_pos_a[aa[ix_pos_a]]
  keep=_cdid_base_index(pos_v:>0)
  if(rows(keep)==0) return(out_c)
  ix_r=ix_pos_a[keep];ai_r=aa[ix_r];pos_r=pos_v[keep]
 }
 else {
  rr_fb=_cdid_db_frame_pair(f,pi_rep,rowsF,ps.cohorts,0)
  if(rows(rr_fb)==0) return(out_c)
  ix_r=rr_fb[,1];ai_r=rr_fb[,2]
  pos_r=J(rows(ai_r),1,0)
  for(q=1;q<=rows(ai_r);q++) pos_r[q]=_cdid_ipw_lower_bound(f.ps_cache[c_ps].risk_a,ai_r[q])
 }
 raw_p=f.ps_cache[c_ps].raw_probability[pos_r]
 odds_r=f.ps_cache[c_ps].odds[pos_r]
 cj_r=cmap[f.plan.a_cohort[ai_r]]
 w_r=_cdid_db_w_ptr(f.reference,f.plan,pi_rep,ai_r)
 fitw=w_r:/ps.rates_c[cj_r]
 X_rr=_cdid_db_X2_ptr(f.reference,f.plan,pi_rep,ai_r)
 Z_load=X_rr*Lambda
 OnPair=factor:*(w_r:*Rate1[cj_r,]+(fitw:*((cj_r:==1):-raw_p)):*Z_load)
 DY=J(rows(ai_r),m_c,0)
 grp_key=(f.method=="dr" ? f.plan.consumer[ells,3] : f.plan.consumer[ells,2])
 grp_ord=order((grp_key,(1::m_c)),(1,2))
 grp_info=_cdid_ps_dbase(grp_key[grp_ord],1)
 for(g_idx=1;g_idx<=rows(grp_info);g_idx++) {
  sub_q=grp_ord[|grp_info[g_idx,1]\grp_info[g_idx,2]|]
  if(f.method=="dr") {
   b=grp_key[sub_q[1]];pi=f.plan.block_pair[b]
  }
  else pi=grp_key[sub_q[1]]
  s_t=f.plan.pair[pi,1];t_t=f.plan.pair[pi,2]
  dy_one=(*f.population).y[f.plan.row_AT[ai_r,t_t]]-(*f.population).y[f.plan.row_AT[ai_r,s_t]]
  if(missing(dy_one))_cdid_base_fail("outcome","Requested pair outcome change overflowed.",504)
  if(f.method=="dr") {
   L=t_t-s_t
   kept=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].kept : f.outcome.kept)
   beta=(f.outcome.outcome_pooling=="separate" ? f.outcome.blocks[b].beta : f.outcome.beta)
   dy_one=dy_one:-f.outcome.alpha[b]:-(rows(beta) ? L:*_cdid_db_Z_kept_ptr(f.reference,f.plan,pi,ai_r,kept)*beta : 0)
  }
  if(rows(sub_q)==1) DY[,sub_q[1]]=dy_one
  else DY[,sub_q]=dy_one*J(1,rows(sub_q),1)
 }
 t_ix=_cdid_base_index(cj_r:==1)
 if(rows(t_ix)) OnPair[t_ix,]=OnPair[t_ix,]+(factor:*fitw[t_ix]:/ps.mean_a):*(DY[t_ix,]:-t_mean)
 c_ix=_cdid_base_index(cj_r:>1)
 if(rows(c_ix)) OnPair[c_ix,]=OnPair[c_ix,]-(factor:*fitw[c_ix]:*odds_r[c_ix]:/ps.mean_b):*(DY[c_ix,]:-c_mean)
 out_c[ix_r,]=OnPair
 return(out_c)
}
real matrix cdid_db_if_rows(struct cdid_db_provider scalar f,real colvector rowsF,
 real colvector consumers,string scalar scale)
{
 real scalar B,nout,i,j,b,c_ps,kb,k,used,use_ps_batch
 real colvector ub,alphaIF,aa,ix_pos_a,ps_ids,ps_ord,idx_c,b_ids,b_ord,idx_b
 real matrix out,betaIF,ps_info,b_info,Dbar
 if(f.reference==NULL|f.population==NULL)_cdid_base_fail("owner","Borrowed original objects are absent.",498)
 B=rows(rowsF);nout=rows(consumers);cdid_ref_activate(*f.reference)
 if(cols(rowsF)!=1|cols(consumers)!=1|missing(rowsF)|missing(consumers)|B>f.chunk_rows|
  any(rowsF:<1)|any(rowsF:>f.n_F)|any(rowsF:!=floor(rowsF))|
  any(consumers:<1)|any(consumers:>f.plan.n_consumers)|any(consumers:!=floor(consumers))|
  (B>1 ? (!all(rowsF[|2\B|]:>rowsF[|1\B-1|]) ? rows(uniqrows(sort(rowsF,1)))!=B : 0) : 0)|
  (nout>1 ? (!all(consumers[|2\nout|]:>consumers[|1\nout-1|]) ? rows(uniqrows(sort(consumers,1)))!=nout : 0) : 0)|
  (scale!="frame"&scale!="core"))_cdid_base_fail("replay","Invalid originalF row/link batch or scale.",498)
 used=_cdid_base_cells(10*B,nout+(f.method=="ipw" ? 0 : (f.outcome.outcome_pooling=="separate" ? f.outcome.p : rows(f.outcome.beta)))+4,f.workspace_limit,f.workspace_entries)
 out=J(B,nout,0);if(B==0|nout==0)return(out)
 if(f.method=="reg") {
  b_ids=f.plan.consumer[consumers,3]
  b_ord=order((b_ids,(1::nout)),(1,2))
  b_info=_cdid_ps_dbase(b_ids[b_ord],1)
  if(f.outcome.outcome_pooling=="separate") {
   for(i=1;i<=rows(b_info);i++) {
    idx_b=b_ord[|b_info[i,1]\b_info[i,2]|];b=b_ids[idx_b[1]]
    _cdid_db_blk_beta_alpha_if(f,b,rowsF,betaIF=.,alphaIF=.)
    kb=rows(f.outcome.blocks[b].beta)
    for(j=1;j<=rows(idx_b);j++) out[,idx_b[j]]=_cdid_db_or_if_nobeta(f,consumers[idx_b[j]],rowsF,alphaIF)
    if(kb) out[,idx_b]=out[,idx_b]-betaIF*f.outcome.link_c[consumers[idx_b],1..kb]'
   }
  }
  else {
   betaIF=_cdid_db_beta_if(f,rowsF);k=rows(f.outcome.beta)
   for(i=1;i<=rows(b_info);i++) {
    idx_b=b_ord[|b_info[i,1]\b_info[i,2]|];b=b_ids[idx_b[1]]
    alphaIF=_cdid_db_alpha_if(f,b,rowsF,betaIF)
    for(j=1;j<=rows(idx_b);j++) out[,idx_b[j]]=_cdid_db_or_if_nobeta(f,consumers[idx_b[j]],rowsF,alphaIF)
   }
   if(k) out=out-betaIF*f.outcome.link_c[consumers,1..k]'
  }
 }
 else {
  aa=f.plan.F_to_A[rowsF];ix_pos_a=_cdid_base_index(aa:>0)
  use_ps_batch=(rows(f.ps_class_id)==f.plan.n_consumers & rows(f.ps_cache)>0 & rows(f.plan.row_AT)==f.n_A)
  if(use_ps_batch) {
   ps_ids=f.ps_class_id[consumers]
   ps_ord=order((ps_ids,(1::nout)),(1,2))
   ps_info=_cdid_ps_dbase(ps_ids[ps_ord],1)
   for(i=1;i<=rows(ps_info);i++) {
    idx_c=ps_ord[|ps_info[i,1]\ps_info[i,2]|];c_ps=ps_ids[idx_c[1]]
    out[,idx_c]=_cdid_db_ps_class_if(f,c_ps,rowsF,aa,ix_pos_a,consumers[idx_c])
   }
   if(f.method=="dr") {
    if(f.outcome.outcome_pooling=="separate") {
     b_ids=f.plan.consumer[consumers,3]
     b_ord=order((b_ids,(1::nout)),(1,2))
     b_info=_cdid_ps_dbase(b_ids[b_ord],1)
     for(i=1;i<=rows(b_info);i++) {
      idx_b=b_ord[|b_info[i,1]\b_info[i,2]|];b=b_ids[idx_b[1]]
      kb=rows(f.outcome.blocks[b].beta)
      if(kb==0)continue
      betaIF=_cdid_db_blk_beta_if(f,b,rowsF)
      Dbar=J(kb,rows(idx_b),0)
      for(j=1;j<=rows(idx_b);j++) {
       if(rows(f.point[consumers[idx_b[j]]].d_bar_diff)==kb) Dbar[,j]=f.point[consumers[idx_b[j]]].d_bar_diff
      }
      out[,idx_b]=out[,idx_b]-betaIF*Dbar
     }
    }
    else {
     k=rows(f.outcome.beta)
     if(k>0) {
      betaIF=_cdid_db_beta_if(f,rowsF)
      Dbar=J(k,nout,0)
      for(j=1;j<=nout;j++) {
       if(rows(f.point[consumers[j]].d_bar_diff)==k) Dbar[,j]=f.point[consumers[j]].d_bar_diff
      }
      out=out-betaIF*Dbar
     }
    }
   }
  }
  else {
   if(f.method=="dr" & f.outcome.outcome_pooling=="separate") {
    ub=uniqrows(sort(f.plan.consumer[consumers,3],1))
    for(i=1;i<=rows(ub);i++) {
     b=ub[i];betaIF=_cdid_db_blk_beta_if(f,b,rowsF)
     for(j=1;j<=nout;j++)if(f.plan.consumer[consumers[j],3]==b)
      out[,j]=_cdid_db_ipw_if(f,consumers[j],rowsF,betaIF)
    }
   }
   else {
    betaIF=(f.method=="dr" ? _cdid_db_beta_if(f,rowsF) : J(B,0,0))
    for(j=1;j<=nout;j++) out[,j]=_cdid_db_ipw_if(f,consumers[j],rowsF,betaIF)
   }
  }
 }
 if(scale=="core")out=out/f.core_to_frame
 if(missing(out))_cdid_base_fail("numeric","Nonfinite complete dated baseline IF.",504)
 return(out)
}
real matrix cdid_db_weight_chunk(struct cdid_db_provider scalar f,real scalar lo,real scalar hi,
 real colvector consumers)
{
 real colvector rowsF
 if(lo<1|hi<lo|hi>f.weight_rows|lo!=floor(lo)|hi!=floor(hi))_cdid_base_fail("replay","Invalid weight replay interval.",498)
 rowsF=(lo::hi);if(f.weight_row_domain=="A")rowsF=f.plan.A_to_F[rowsF]
 return(cdid_db_if_rows(f,rowsF,consumers,"core"))
}
struct cdid_db_weight_context {
 pointer(struct cdid_db_provider scalar) scalar provider
 real colvector consumers,sign
}
struct cdid_db_weight_context scalar cdid_db_weight_context_make(
 pointer(struct cdid_db_provider scalar) scalar provider,real colvector consumers,real colvector sign)
{
 struct cdid_db_weight_context scalar ctx
 if(provider==NULL|rows(consumers)!=rows(sign)|cols(consumers)!=1|cols(sign)!=1|missing(sign)|any(abs(sign):!=1))
  _cdid_base_fail("context","Provider, explicit consumers and signs required.",498)
 ctx.provider=provider;ctx.consumers=consumers;ctx.sign=sign;return(ctx)
}
real matrix cdid_db_weight(pointer() scalar context,real scalar lo,real scalar hi)
{
 struct cdid_db_weight_context scalar ctx
 if(context==NULL)_cdid_base_fail("context","Weight context is absent.",498)
 ctx=*context
 return(cdid_db_weight_chunk(*ctx.provider,lo,hi,ctx.consumers):*ctx.sign')
}
end
