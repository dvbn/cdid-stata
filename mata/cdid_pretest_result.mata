*! Standalone conditional-test result record; no ATT covariance or synthetic points.
version 19.5
mata:
struct cdid_pt_result {
 string scalar mode,method,status,reason,formula,clusters,rng,policy,channel,membership
 real rowvector stats
 real colvector draws,retained_rows,grid_first,multiplicity,cluster_positions,cluster_values
 real matrix family,rowmap,preprocess,missing_grid,unsupported_grid
 string matrix request
 string rowvector family_columns,xcols,zcols,raw_grid_variables,preprocess_columns
 string colvector family_reasons,cluster_labels,id_labels,grid_keys,cell_ids,cell_status
 real scalar bstrap_requested,cband_requested,multipliers_used,raw_CvM
}
void _cdid_pt_fail(string scalar msg,real scalar rc)
{
    errprintf("%s\n", msg)
    exit(rc)
}
void cdid_pt_assert_free(string scalar name)
{
 if(name!="")if(findexternal(name)!=NULL)_cdid_pt_fail("Requested conditional result object already exists.",498)
}
real scalar cdid_pt_cost(struct cdid_pt_result scalar R)
{
 real scalar textbytes
 textbytes=sum(strlen((R.mode,R.method,R.status,R.reason,R.formula,R.clusters,R.rng,R.policy,R.channel,R.membership)))+
  sum(strlen(R.family_columns))+sum(strlen(R.xcols))+sum(strlen(R.zcols))+sum(strlen(R.family_reasons))+
  sum(strlen(R.cluster_labels))+sum(strlen(R.id_labels))+sum(strlen(R.grid_keys))+sum(strlen(R.raw_grid_variables))+
  sum(strlen(R.preprocess_columns))+sum(strlen(R.cell_ids))+sum(strlen(R.cell_status))+sum(strlen(R.request))
 return(ceil(textbytes/8)+length(R.stats)+length(R.draws)+length(R.retained_rows)+length(R.grid_first)+length(R.multiplicity)+
 length(R.cluster_positions)+length(R.cluster_values)+length(R.family)+length(R.rowmap)+length(R.preprocess)+
 length(R.missing_grid)+length(R.unsupported_grid))
}
void cdid_pt_validate(struct cdid_pt_result scalar R,real scalar cap)
{
 if((R.mode!="ordinary" & R.mode!="chained") | (R.status!="available" & R.status!="unavailable" & R.status!="not_applicable") |
  rows(R.stats)!=1 | cols(R.stats)!=13 | rows(R.draws)>R.stats[5] | cols(R.draws)!=1 |
  rows(R.family)!=rows(R.family_reasons) | cols(R.family)!=cols(R.family_columns) |
  cdid_pt_cost(R)>cap)_cdid_pt_fail("Invalid or over-budget conditional result record.",498)
 if(R.status=="available") {
  if(rows(R.draws)!=R.stats[5] | missing(R.draws) | missing(R.stats[7..9]))
   _cdid_pt_fail("Available conditional result lacks finite statistics or complete B-by-one draws.",498)
 }
 else if(rows(R.draws) | any(R.stats[7..9]:<.))_cdid_pt_fail("Unavailable conditional public fields must be missing with absent draws.",498)
}
void cdid_pt_put(string scalar name,struct cdid_pt_result scalar R,real scalar cap)
{
 pointer() scalar p
 cdid_pt_validate(R,cap)
 p=crexternal(name);if(p==NULL)_cdid_pt_fail("Requested conditional result object already exists.",498)
 *p=R
}
void cdid_pt_drop(string scalar name)
{
 struct cdid_pt_result scalar R
 pointer() scalar p
 p=findexternal(name);if(p==NULL)return
 if(eltype(*p)!="struct")_cdid_pt_fail("Requested object is not a conditional result.",498)
 R=*p
 if(structname(R)!="cdid_pt_result")_cdid_pt_fail("Requested object has a different result kind.",498)
 rmexternal(name)
}
void cdid_pt_from_ordinary(string scalar source,string scalar output,real scalar bs,real scalar cb,real scalar cap)
{
 struct cdid_op_state scalar O
 struct cdid_pt_result scalar R
 O=*findexternal(source)
 R.mode="ordinary";R.method=O.method;R.status=O.status;R.reason=O.reason;R.formula=O.formula
 R.rng=O.rng;R.policy=O.policy;R.channel=O.channel;R.membership=O.membership
 R.stats=O.stats;R.draws=O.draws;R.family=O.family;R.family_reasons=O.family_reasons
 R.family_columns=("cell","group","time","base","unavailable")
 R.preprocess_columns=("source_row","supplied_id","time","cohort","raw_weight","row_unit")
 R.raw_grid_variables=cdid_ref_formula_names(O.formula)
 R.retained_rows=O.retained_rows;R.grid_first=O.grid_first;R.multiplicity=O.multiplicity
 R.cluster_positions=O.cluster_positions;R.cluster_values=O.cluster_values;R.cluster_labels=O.cluster_labels
 R.rowmap=O.rowmap;R.preprocess=O.preprocess;R.xcols=O.xcols;R.zcols=O.zcols
 R.bstrap_requested=bs;R.cband_requested=cb;R.multipliers_used=(O.status=="available")
 R.raw_CvM=O.stats[7]
 R.missing_grid=R.unsupported_grid=J(0,4,.)
 cdid_pt_put(output,R,cap)
}
void cdid_pt_save(string scalar name,string scalar path,real scalar cap)
{
 struct cdid_pt_result scalar R
 real scalar fh,rc
 R=*findexternal(name);cdid_pt_validate(R,cap)
 fh=_fopen(path,"w")
 if(fh<0)_cdid_pt_fail("Cannot create conditional result file; existing files are never overwritten.",-fh)
 _cdid_rs_write(fh,("cdid-native-conditional-result","1"))
 _cdid_rs_write(fh,cdid_pt_cost(R));_cdid_rs_write(fh,R)
 rc=_fclose(fh);if(rc<0)_cdid_pt_fail("Conditional result close failed; partial file retained.",-rc)
}
void cdid_pt_load(string scalar path,string scalar output,real scalar cap)
{
 real scalar fh,rc
 transmorphic matrix head,cost,raw,tail
 struct cdid_pt_result scalar R
 if(findexternal(output)!=NULL)_cdid_pt_fail("Requested conditional result object already exists.",498)
 fh=_fopen(path,"r");if(fh<0)_cdid_pt_fail("Cannot open conditional result file.",-fh)
 head=_cdid_rs_read(fh)
 if(eltype(head)!="string" | rows(head)!=1 | cols(head)!=2){(void)_fclose(fh);_cdid_pt_fail("Unknown conditional result header.",498);}
 if(head[1]!="cdid-native-conditional-result" | head[2]!="1"){(void)_fclose(fh);_cdid_pt_fail("The file was not saved by cdid_pretest, or by a version of cdid that cdid_pretest_load cannot read.",498);}
 cost=_cdid_rs_read(fh)
 if(eltype(cost)!="real" | rows(cost)!=1 | cols(cost)!=1){(void)_fclose(fh);_cdid_pt_fail("Invalid conditional result size record.",498);}
 if(cost>=. | cost<0 | cost>cap | cost!=floor(cost)){(void)_fclose(fh);_cdid_pt_fail("Conditional result exceeds explicit numeric storage.",908);}
 raw=_cdid_rs_read(fh)
 tail=_fgetmatrix(fh,1);rc=fstatus(fh);(void)_fclose(fh)
 if(rc!=-1 | eltype(raw)!="struct")_cdid_pt_fail("Invalid conditional result payload or end marker.",498)
 R=raw;if(structname(R)!="cdid_pt_result")_cdid_pt_fail("Conditional result record has a different type.",498)
 cdid_pt_validate(R,cap);if(cdid_pt_cost(R)!=cost)_cdid_pt_fail("Conditional result size record disagrees with payload.",498)
 cdid_pt_put(output,R,cap)
}
void cdid_pt_request_meta(string scalar name)
{
 pointer(struct cdid_pt_result scalar) scalar pR
 string rowvector fields
 string matrix req
 real scalar i
 pR=findexternal(name);pR->clusters=st_local("clustervars")
 fields=tokens("request time gvar ivar covars noconstant if in weight exp mode panel rcs allowunbalanced method controlgroup wvar cluster clustervars bstrap cband biters level pl cores anticipation baseperiod fixweights frame frameivar framegvar framecovars framewvar framecluster selection selmodel selpooling selcovars frameselcovars selnoconstant selconstant selzero seluniverse selrho seljustification selmethod imppooling outcomepooling workspace chunk workcap posting store saving")
 req=J(cols(fields)+3,3,"")
 for(i=1;i<=cols(fields);i++)req[i,]=(fields[i],st_local("raw_"+fields[i]),st_local(fields[i]))
 req[cols(fields)+1,]=("caller_rng_before",st_local("rng0"),st_local("rng0"))
 req[cols(fields)+2,]=("caller_rng_after",st_local("rng_after"),st_local("rng_after"))
 req[cols(fields)+3,]=("sorting_rng_returned",st_local("sort0"),st_local("sort0"))
 pR->request=req
}
void cdid_pt_export(string scalar name,string scalar stats,string scalar draws,string scalar family,
 string scalar grid,string scalar posting,real scalar matsize,real scalar cap)
{
 pointer(struct cdid_pt_result scalar) scalar pR
 pR=findexternal(name);cdid_pt_validate(*pR,cap)
 if(posting!="matrix" & posting!="state")_cdid_pt_fail("posting() must be matrix or state.",198)
 if(posting=="matrix") {
  if(max((rows(pR->draws),rows(pR->family),cols(pR->family),rows(pR->grid_first)))>matsize)
   _cdid_pt_fail("Conditional matrix posting exceeds matsize; choose posting(state) before drawing.",908)
  if(rows(pR->draws))st_matrix(draws,pR->draws)
  if(rows(pR->family))st_matrix(family,pR->family)
  if(rows(pR->grid_first))st_matrix(grid,(pR->grid_first,pR->multiplicity))
 }
 st_matrix(stats,pR->stats)
 st_local("pt_status",pR->status);st_local("pt_reason",pR->reason);st_local("pt_mode",pR->mode)
 st_local("pt_method",pR->method);st_local("pt_formula",pR->formula);st_local("pt_rng",pR->rng)
 st_local("pt_policy",pR->policy);st_local("pt_channel",pR->channel);st_local("pt_membership",pR->membership)
 st_local("pt_xcols",invtokens(pR->xcols));st_local("pt_zcols",invtokens(pR->zcols))
 st_local("pt_family_columns",invtokens(pR->family_columns));st_local("pt_cluster",pR->clusters)
 st_local("pt_draws",strofreal(rows(pR->draws)));st_local("pt_family",strofreal(rows(pR->family)))
 st_local("pt_grid",strofreal(rows(pR->grid_first)));st_local("pt_used",strofreal(pR->multipliers_used))
}
end
