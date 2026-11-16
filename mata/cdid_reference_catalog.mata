*! Shared scientific consumer registry qualified by reference recipe; StataNow19.5.
*! One raw map for baseline, selected and imputation consumers. No encoding or fit.
version 19.5
mata:
struct cdid_reference_catalog {
    real scalar M,Q,n_raw,n_consumers,n_blocks
    real matrix consumers,blocks
    real colvector representative,consumer_recipe,consumer_block
    real colvector K_ptr,K_index,block_K_ptr,block_K_index
    real colvector raw_consumer,post_consumer,aux_consumer,raw_supported,aux_supported
    string colvector keys,block_keys,raw_reason,aux_reason
}
string scalar _cdid_refcat_key(struct cdid_baseline_plan_state scalar p,
    real scalar i,real scalar recipe,real scalar include_g)
{
    string scalar key
    real scalar k
    key=(include_g ? strofreal(p.raw[i,1],"%21.0f")+":" : "")+
        strofreal(p.raw[i,2],"%21.0f")+":"+strofreal(p.raw[i,3],"%21.0f")+":"+strofreal(recipe,"%21.0f")+":"
    for(k=p.raw_K_ptr[i];k<p.raw_K_ptr[i+1];k++)key=key+strofreal(p.raw_K_index[k],"%21.0f")+"," 
    return(key)
}
struct cdid_reference_catalog scalar cdid_reference_catalog_make(
    struct cdid_baseline_plan_state scalar plan,
    struct cdid_reference_inputs scalar reference,real scalar cap)
{
    struct cdid_reference_catalog scalar out
    transmorphic scalar seen,blocks
    real scalar i,r,c,b,used,nk
    real colvector K
    string scalar key,bkey
    out.M=plan.M;out.Q=plan.Q;out.n_raw=rows(plan.raw)
    if(out.n_raw!=out.M+out.Q | rows(reference.view_recipe)!=out.n_raw)
        _cdid_ref_fail("Reference design blocks must align with ALL original post and both-auxiliary raw views.",498)
    out.raw_supported=plan.raw_supported;out.aux_supported=plan.aux_supported
    out.raw_reason=J(out.n_raw,1,"");out.aux_reason=plan.aux_reason
    out.raw_consumer=J(out.n_raw,1,0)
    out.consumers=J(0,4,.);out.blocks=J(0,3,.)
    out.representative=out.consumer_recipe=out.consumer_block=J(0,1,0)
    out.K_ptr=out.block_K_ptr=1;out.K_index=out.block_K_index=J(0,1,0)
    out.keys=out.block_keys=J(0,1,"")
    seen=asarray_create("string",1);blocks=asarray_create("string",1)
    used=_cdid_ref_cost(out.n_raw,3,cap,0)
    for(i=1;i<=out.n_raw;i++) {
        if(!plan.raw_supported[i]) {
            if(i>plan.M)out.raw_reason[i]=plan.aux_reason[i-plan.M]
            continue
        }
        r=reference.view_recipe[i]
        if(r<1 | r>reference.n_recipes | r!=floor(r))_cdid_ref_fail("Invalid raw-view reference design block.",498)
        key=_cdid_refcat_key(plan,i,r,1);bkey=_cdid_refcat_key(plan,i,r,0)
        if(!asarray_contains(seen,key)) {
            K=plan.raw_K_index[plan.raw_K_ptr[i]..plan.raw_K_ptr[i+1]-1];nk=rows(K)
            used=_cdid_ref_cost(1,7+nk,cap,used)
            c=rows(out.consumers)+1
            out.consumers=out.consumers\(plan.raw[i,],r)
            out.representative=out.representative\i;out.consumer_recipe=out.consumer_recipe\r
            out.K_index=out.K_index\K;out.K_ptr=out.K_ptr\(rows(out.K_index)+1)
            out.keys=out.keys\key
            if(!asarray_contains(blocks,bkey)) {
                used=_cdid_ref_cost(1,4+nk,cap,used)
                b=rows(out.blocks)+1;asarray(blocks,bkey,b)
                out.blocks=out.blocks\(plan.raw[i,2],plan.raw[i,3],r)
                out.block_keys=out.block_keys\bkey
                out.block_K_index=out.block_K_index\K;out.block_K_ptr=out.block_K_ptr\(rows(out.block_K_index)+1)
            }
            b=asarray(blocks,bkey);out.consumer_block=out.consumer_block\b
            asarray(seen,key,c)
        }
        out.raw_consumer[i]=asarray(seen,key)
    }
    out.n_consumers=rows(out.consumers);out.n_blocks=rows(out.blocks)
    out.post_consumer=(out.M ? out.raw_consumer[(1::out.M)] : J(0,1,0))
    out.aux_consumer=(out.Q ? out.raw_consumer[((out.M+1)::out.n_raw)] : J(0,1,0))
    return(out)
}
end
