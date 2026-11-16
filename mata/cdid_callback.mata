*! Private callback protocol; no public callback/payload spelling is implied.
*! Fixed source: R cdid/, chained_custom_ipw.R45-138.
version 19.5
mata:
struct cdid_callback_request {
 string scalar mode,route,link_id
 real scalar panel,inffunc,gradient_required,N,n_A,call_index,raw_index,model_index,g,s,t,recipe,grid_index
 real colvector rows_F,rows_A,y1,y0,y,D,post,i_weights,K
 real matrix covariates
 string rowvector covariate_names
}
struct cdid_callback_result {
 real scalar has_ATT,has_psi,has_gradient
 real matrix ATT,psi,log_weight_gradient
}
void _cdid_callback_fail(string scalar cls,string scalar message,real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
struct cdid_callback_request scalar cdid_callback_request_init()
{
 struct cdid_callback_request scalar q
 q.mode=q.route=q.link_id=""
 q.panel=q.inffunc=q.gradient_required=0
 q.N=q.n_A=q.call_index=q.raw_index=q.model_index=q.g=q.s=q.t=q.recipe=q.grid_index=.
 q.rows_F=q.rows_A=q.y1=q.y0=q.y=q.D=q.post=q.i_weights=q.K=J(0,1,.)
 q.covariates=J(0,0,.);q.covariate_names=J(1,0,"")
 return(q)
}
/* Real matrices intentionally permit malformed dimensions for truthful validation.
   The private typed protocol replaces R named fields by explicit presence flags. */
struct cdid_callback_result scalar cdid_callback_result_make(real matrix att,
 real matrix psi,real matrix gradient)
{
 struct cdid_callback_result scalar out
 out.has_ATT=out.has_psi=1;out.has_gradient=(length(gradient)>0)
 out.ATT=att;out.psi=psi;out.log_weight_gradient=gradient
 return(out)
}
struct cdid_callback_result scalar cdid_callback_call(pointer(function) scalar callback,
 struct cdid_callback_request scalar request,pointer() scalar payload)
{
 if(callback==NULL)_cdid_callback_fail("input","A scientific callback function pointer is required.",498)
 return((*callback)(request,payload))
}
/* Strict CHAINED result law only. Ordinary cell warnings/NA are a separate adapter.
   Validate optional g before requiring it, as the selected source actually does. */
void cdid_callback_validate_chained(struct cdid_callback_result scalar out,
 real scalar m,real scalar require_gradient)
{
 if(m<2 | m>=. | m!=floor(m) | (require_gradient!=0 & require_gradient!=1))
  _cdid_callback_fail("input","Invalid observed count or gradient requirement.",498)
 if(out.has_ATT!=1 | out.has_psi!=1)
  _cdid_callback_fail("output","Callback result must contain ATT and att.inf.func.",498)
 if(rows(out.ATT)!=1 | cols(out.ATT)!=1 | missing(out.ATT))
  _cdid_callback_fail("output_numeric","Callback ATT must be a finite numeric scalar.",504)
 if(rows(out.psi)!=m | cols(out.psi)!=1)
  _cdid_callback_fail("output_dimension","Callback att.inf.func must have one entry per passed row.",503)
 if(missing(out.psi))_cdid_callback_fail("output_numeric","Callback att.inf.func contains nonfinite entries.",504)
 if(out.has_gradient!=0 & out.has_gradient!=1)
  _cdid_callback_fail("output","Invalid optional log_weight_gradient presence flag.",498)
 if(out.has_gradient) {
  if(rows(out.log_weight_gradient)!=m | cols(out.log_weight_gradient)!=1)
   _cdid_callback_fail("output_dimension","Callback log_weight_gradient must have one entry per passed row.",503)
  if(missing(out.log_weight_gradient))
   _cdid_callback_fail("output_numeric","Callback log_weight_gradient contains nonfinite entries.",504)
 }
 else if(length(out.log_weight_gradient))
  _cdid_callback_fail("output","Absent log_weight_gradient must have an empty private field.",498)
 if(require_gradient & !out.has_gradient)
  _cdid_callback_fail("derivative","The selected route requires explicit total log_weight_gradient after output validation.",498)
}
struct cdid_callback_result scalar cdid_callback_invoke_chained(pointer(function) scalar callback,
 struct cdid_callback_request scalar request,pointer() scalar payload)
{
 struct cdid_callback_result scalar out
 out=cdid_callback_call(callback,request,payload)
 cdid_callback_validate_chained(out,rows(request.D),request.gradient_required)
 return(out)
}
/* g is a TOTAL derivative in ATT units. Never divide supplied g by m or N.
   Only selection=none uses the released omitted-g identity psi/m rule. */
real colvector cdid_callback_rate_loading(struct cdid_callback_result scalar out,
 real colvector group,real colvector rates)
{
 real scalar k,m
 real colvector b,ix,value
 real matrix clipped
 m=rows(group);cdid_callback_validate_chained(out,m,0)
 if(cols(group)!=1 | missing(group) | any(group:<1) | any(group:>rows(rates)) | any(group:!=floor(group)))
  _cdid_callback_fail("input","Invalid passed-row cohort positions.",498)
 clipped=cdid_ipw_rate_values(rates);b=J(rows(rates),1,0)
 value=(out.has_gradient ? out.log_weight_gradient : out.psi:/m)
 for(k=1;k<=rows(rates);k++)if(clipped[k,2]) {
  ix=vec(selectindex(group:==k));if(rows(ix))b[k]=-sum(value[ix])/rates[k]
 }
 if(missing(b))_cdid_callback_fail("numeric","Nonfinite callback cohort-rate sensitivity.",504)
 return(b)
}
end
