*! Ordinary policy adapter for bounded native caller/fixed sign stream.
*! Requires cdid_inference budgets and cdid_public_rng; no chained reducer used.
version 19.5
mata:
struct cdid_ord_boot scalar cdid_ord_boot_stream(struct cdid_ord_summary scalar S,
    real scalar B, real scalar alpha, real scalar return_V,
    struct cdid_infer_budget scalar budget, string scalar mode, real matrix fixed)
{
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar batch
    real matrix Z
    real scalar storage, cap
    /* Reuse arithmetic budget checks only, not chained IF/degeneracy policies. */
    storage=_cdid_infer_preflight(budget,S.n,S.C,S.K,B)
    cap=min((budget.max_workspace_cells,floor(budget.max_workspace_bytes/8)))
    if(B*S.K>cap-storage) _cdid_ord_fail("Ordinary raw draw storage exceeds explicit workspace.",908)
    stream=cdid_public_rng_start(B,S.C,budget,mode,fixed)
    Z=J(B,S.K,0)
    while(stream.cursor<B) {
        batch=cdid_public_rng_next(stream)
        cdid_ord_draw_add(Z,S.sums,batch.draw_positions,batch.cluster_positions,batch.signs,S.C)
    }
    cdid_public_rng_complete(stream)
    return(cdid_ord_bres_finish(Z,(1::S.K),S.K,S.n,S.C,alpha,return_V))
}
/* Schedule is bin getSE calls, one family call if requested, overall getSE.
   Each fixed-call matrix has its own explicit map entry. */
transmorphic scalar cdid_ord_aggregate_draws(struct cdid_ord_aggregate_result scalar O,
    struct cdid_ord_agg_request scalar R, real colvector cluster, real scalar C,
    real scalar clustered, real scalar row_chunk, struct cdid_infer_budget scalar budget,
    string scalar mode, transmorphic scalar fixed_calls)
{
    transmorphic scalar calls
    struct cdid_ord_summary scalar S
    struct cdid_ord_boot scalar M
    real scalar j, call
    string scalar key
    real matrix fixed
    calls=asarray_create("string",1);call=0
    if(O.status!="computed") _cdid_ord_fail("Refused aggregation cannot draw multipliers.",459)
    for(j=1;j<=cols(O.IF);j++) if(R.bstrap) {
        call=call+1;key=strtrim(strofreal(call,"%21.0f"));fixed=J(0,0,.)
        if(mode=="fixed") {
            if(!asarray_contains(fixed_calls,key)) _cdid_ord_fail("Missing declared per-bin fixed signs.",498)
            fixed=asarray(fixed_calls,key)
        }
        S=cdid_ord_if_summary(O.IF[,j],cluster,C,clustered,R.bstrap,row_chunk)
        M=cdid_ord_boot_stream(S,R.B,R.alpha,0,budget,mode,fixed);asarray(calls,key,M)
    }
    if(cols(O.IF) & R.cband) {
        call=call+1;key=strtrim(strofreal(call,"%21.0f"));fixed=J(0,0,.)
        if(mode=="fixed") {
            if(!asarray_contains(fixed_calls,key)) _cdid_ord_fail("Missing separate fixed family signs.",498)
            fixed=asarray(fixed_calls,key)
        }
        S=cdid_ord_if_summary(O.IF,cluster,C,clustered,R.bstrap,row_chunk)
        M=cdid_ord_boot_stream(S,R.B,R.alpha,0,budget,mode,fixed);asarray(calls,key,M)
    }
    if(R.bstrap) {
        call=call+1;key=strtrim(strofreal(call,"%21.0f"));fixed=J(0,0,.)
        if(mode=="fixed") {
            if(!asarray_contains(fixed_calls,key)) _cdid_ord_fail("Missing separate fixed overall signs.",498)
            fixed=asarray(fixed_calls,key)
        }
        S=cdid_ord_if_summary(O.overall_IF,cluster,C,clustered,R.bstrap,row_chunk)
        M=cdid_ord_boot_stream(S,R.B,R.alpha,0,budget,mode,fixed);asarray(calls,key,M)
    }
    return(calls)
}
end
