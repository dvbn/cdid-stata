*! Private bounded caller-stream Rademacher adapter, current StataNow19.5.
*! R cdid chained_ipw_pipeline.R4090-4122: caller RNG advances, no seed reset;
*! batches follow draw order and column-major cluster order within each batch.
*! Stata and R RNG sequences are distinct. Fixed signs consume no RNG.
version 19.5
mata:
struct cdid_public_sign_stream {
    real scalar B,C,batch_rows,cursor,finished
    string scalar mode,ordering
    real matrix fixed
}
struct cdid_public_sign_batch {
    real colvector draw_positions
    real rowvector cluster_positions
    real matrix signs
}
void _cdid_pr_rng_fail(string scalar message,real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
struct cdid_public_sign_stream scalar cdid_public_rng_start(real scalar B,real scalar C,
    struct cdid_infer_budget scalar budget,string scalar mode,real matrix fixed)
{
    struct cdid_public_sign_stream scalar S
    real scalar cap
    if(B<1 | C<1 | B>=. | C>=. | B!=floor(B) | C!=floor(C) | B>2^31-1 | C>2^31-1)
        _cdid_pr_rng_fail("Explicit positive integer B and C are required.",498)
    if(B>budget.max_draws | C>budget.max_clusters) _cdid_pr_rng_fail("Draw/cluster dimensions exceed the explicit inference budget.",908)
    S.batch_rows=min((B,budget.max_sign_batch_rows,floor(budget.max_sign_batch_cells/C)))
    if(S.batch_rows<1) _cdid_pr_rng_fail("The sign budget cannot hold a complete declared cluster row.",908)
    if(mode!="caller" & mode!="fixed") _cdid_pr_rng_fail("Sign mode must be caller or fixed; no seed or alternative distribution is substituted.",498)
    if(mode=="fixed") {
        if(rows(fixed)!=B | cols(fixed)!=C | missing(fixed) | any((fixed:!=1):&(fixed:!=-1)))
            _cdid_pr_rng_fail("Fixed signs must have exact B by C shape with values -1 or +1.",498)
    }
    else if(rows(fixed)*cols(fixed)) _cdid_pr_rng_fail("Caller RNG mode cannot also supply fixed signs.",498)
    cap=min((budget.max_workspace_cells,floor(budget.max_workspace_bytes/8)))
    if(2*S.batch_rows*C+(mode=="fixed")*B*C>cap) _cdid_pr_rng_fail("Sign storage exceeds the explicit numeric workspace budget.",908)
    S.B=B;S.C=C;S.mode=mode;S.cursor=0;S.finished=0
    S.ordering="increasing draw batches; column-major within each batch"
    S.fixed=(mode=="fixed" ? fixed : J(0,0,.))
    return(S)
}
struct cdid_public_sign_batch scalar cdid_public_rng_next(struct cdid_public_sign_stream scalar S)
{
    struct cdid_public_sign_batch scalar E
    real scalar b
    real colvector u
    if(S.finished | S.cursor>=S.B) _cdid_pr_rng_fail("The declared sign stream is exhausted or finalized.",498)
    b=min((S.batch_rows,S.B-S.cursor))
    E.draw_positions=((S.cursor+1)::(S.cursor+b));E.cluster_positions=1..S.C
    if(S.mode=="fixed") E.signs=S.fixed[E.draw_positions,]
    else {
        // All validation/allocation sizes are known before touching caller RNG.
        u=runiform(b*S.C,1)
        E.signs=2:*(colshape(u,b)':>=.5):-1
    }
    S.cursor=S.cursor+b
    return(E)
}
void cdid_public_rng_complete(struct cdid_public_sign_stream scalar S)
{
    if(S.finished | S.cursor!=S.B) _cdid_pr_rng_fail("Every declared draw must be consumed exactly once before finalization.",498)
    S.finished=1
}
end
