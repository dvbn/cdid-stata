*! Private old-did preprocessing and ordinary conditional result consumer, StataNow19.5.
*! Reference: did2.5.1 pre_process_did.R and conditional_did_pretest.R; no R code copied.
*! Requires cdid_oc_*, cdid_co_*, cdid_public_rng_* and common cdid_ref_encoded/formula_names.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_opre(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_op_input {
    struct cdid_oc_recipe scalar recipe
    real scalar inputN, requestedN, panel, allow, intercept, maxentries, varying_w
    real matrix rowmap
    real colvector rows, ids, glist, gridrows, cluster_raw
    string colvector cluster_labels
    string rowvector xcols, zcols
    string scalar formula, cluster, idname, callback_name, payload_name
}
struct cdid_op_state {
    real rowvector stats
    real colvector draws, multiplicity, grid_first, retained_rows, cluster_positions, cluster_values
    real matrix family, rowmap, preprocess
    string colvector family_reasons, cluster_labels
    string rowvector xcols, zcols
    string scalar status, reason, method, formula, membership, channel, rng, policy
    string scalar callback_name, payload_name
    real scalar n_callbacks
}
void _cdid_op_fail(string scalar message,real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real scalar _cdid_op_is_ident(string scalar s)
{
    real scalar n,i,c
    n=strlen(s)
    if(n<1 | n>32)return(0)
    c=ascii(substr(s,1,1))
    if(!((c>=65 & c<=90) | (c>=97 & c<=122) | c==95))return(0)
    for(i=2;i<=n;i++) {
        c=ascii(substr(s,i,1))
        if(!((c>=65 & c<=90) | (c>=97 & c<=122) | (c>=48 & c<=57) | c==95))return(0)
    }
    return(1)
}
void _cdid_op_validate_symbols(string scalar callback_name,string scalar payload_name)
{
    if(!_cdid_op_is_ident(callback_name)) {
        errprintf("method(custom) requires callback() naming a Mata function.\n")
        exit(198)
    }
    if(findexternal(callback_name+"()")==NULL) {
        errprintf("Mata function %s() not found; define it before calling the command.\n",callback_name)
        exit(198)
    }
    if(payload_name!="") {
        if(!_cdid_op_is_ident(payload_name)) {
            errprintf("payload() must name a Mata global object.\n")
            exit(198)
        }
        if(findexternal(payload_name)==NULL) {
            errprintf("payload(%s) was not found as an active external Mata object.\n",payload_name)
            exit(198)
        }
    }
}
real colvector _cdid_op_ix(real matrix p)
{
    return(vec(selectindex(p)))
}
/* First appearance labels, including string cluster labels. No inference recoding is
   permitted to reorder labels. Numeric keys use their exact binary64 hex representation. */
real colvector _cdid_op_cmap(string colvector labels,string colvector distinct)
{
    transmorphic scalar A
    real colvector out
    real scalar i,j
    A=asarray_create("string",1);out=J(rows(labels),1,0);distinct=J(0,1,"")
    for(i=1;i<=rows(labels);i++) {
        if(!asarray_contains(A,labels[i])) {j=rows(distinct)+1;distinct=distinct\labels[i];asarray(A,labels[i],j);}
        out[i]=asarray(A,labels[i])
    }
    return(out)
}
/* Raw source order is held until the final (id,time) order; true RCS .rowid is
   created after screening/support, so that its final order is the retained raw order. */
void cdid_op_intake(string scalar name,string scalar sample,string scalar evalmask,string scalar y,
    string scalar time,string scalar group,string scalar id,string scalar weight,
    string scalar clusters,string scalar formula,real scalar panel,real scalar allow,
    real scalar intercept,string scalar method,string scalar control,real scalar cap,
    string scalar ordervar,string scalar groupvar,|string scalar callback_name,string scalar payload_name)
{
    struct cdid_op_input scalar O
    pointer scalar p
    real colvector rr,t,g,w,iv,ok,unit,ord,drop,tl,gl,ix,uq,cnt,gfirst,keep,cl
    real matrix raw,info,trip
    string rowvector xnames,cn
    string colvector labs,discard
    real scalar n,j,latest,first,nraw
    if(findexternal(name)!=NULL)_cdid_op_fail("Temporary intake object already exists.",498)
    if(method=="custom")_cdid_op_validate_symbols(callback_name,payload_name)
    xnames=cdid_ref_formula_names(formula);cn=tokens(clusters)
    if(cols(cn)>2 | (cols(cn)==2 & (id=="" | !any(cn:==id))))
        _cdid_op_fail("At most one cluster variable beyond ivar() is supported.",198)
    rr=_cdid_op_ix(st_data(.,sample):!=0);n=rows(rr)
    if(n==0)_cdid_op_fail("No requested observations.",2000)
    if(2*st_nobs()+n*(cols(xnames)+10)>cap)_cdid_op_fail("Raw intake exceeds the explicit workspace budget.",908)
    t=st_data(rr,time);g=st_data(rr,group)
    if(any(g:<0))_cdid_op_fail("Cohorts must be zero or positive treatment dates.",498)
    raw=st_data(rr,y);w=(weight=="" ? J(n,1,1) : st_data(rr,weight))
    ok=(raw:<.):&(t:<.):&(g:<.):&(w:<.)
    if(cols(xnames)) {raw=st_data(rr,xnames);ok=ok:*(rowmissing(raw):==0);}
    if(id!="") {
        if(st_isnumvar(id)) {iv=st_data(rr,id);ok=ok:*(iv:<.);}
        else {labs=st_sdata(rr,id);ok=ok:*(labs:!="");iv=J(n,1,.);}
    }
    else iv=(1::n)
    for(j=1;j<=cols(cn);j++) {
        if(st_isnumvar(cn[j]))ok=ok:*(st_data(rr,cn[j]):<.)
        else ok=ok:*(st_sdata(rr,cn[j]):!="")
    }
    O.rowmap=((1::st_nobs()),J(st_nobs(),1,0))
    O.rowmap[rr,2]=J(n,1,1)
    O.rowmap[select(rr,ok:==0),2]=J(sum(ok:==0),1,2)
    keep=ok:&(st_data(rr,evalmask):==0)
    O.rowmap[select(rr,keep),2]=J(sum(keep),1,6)
    ok=ok:&(st_data(rr,evalmask):!=0)
    rr=select(rr,ok);t=select(t,ok);g=select(g,ok);w=select(w,ok);iv=select(iv,ok)
    if(!rows(rr))_cdid_op_fail("All requested observations were missing or nonfinite.",2000)
    if(any(w:<0) | mean(w)<=0)_cdid_op_fail("Raw sampling weights must be nonnegative with a positive mean.",498)
    nraw=rows(rr);O.varying_w=0
    if(panel & id!="") {
        if(st_isnumvar(id)) {
        ord=order((iv,t),(1,2));info=_cdid_ps_opre(iv[ord],1)
        for(j=1;j<=rows(info);j++) {
            ix=ord[|info[j,1]\info[j,2]|]
            if(max(w[ix])-min(w[ix])>sqrt(epsilon(1)))O.varying_w=1
        }
        }
    }
    tl=uniqrows(sort(t,1))
    g= (g:>max(tl)):*0 + (g:<=max(tl)):*g
    gl=uniqrows(sort(g,1));latest=.
    if(!any(gl:==0)) {
        latest=max(gl);keep=t:<latest
        O.rowmap[select(rr,keep:==0),2]=J(sum(keep:==0),1,3)
        rr=select(rr,keep);t=select(t,keep);g=select(g,keep);w=select(w,keep);iv=select(iv,keep)
        if(control=="nevertreated")g=g:*(g:!=latest)
        printf("No never-treated group: %s; periods at or after %g excluded.\n",(control=="nevertreated" ? "latest cohort coerced to never-treated" : "latest cohort retained only as control"),latest)
    }
    if(!rows(rr))_cdid_op_fail("No observations remain after the no-never-treated cutoff.",2000)
    tl=uniqrows(sort(t,1));first=tl[1];keep=(g:==0):|(g:>first)
    O.rowmap[select(rr,keep:==0),2]=J(sum(keep:==0),1,4)
    rr=select(rr,keep);t=select(t,keep);g=select(g,keep);w=select(w,keep);iv=select(iv,keep)
    if(!rows(rr))_cdid_op_fail("No valid groups remain after first-period treatment exclusions.",498)
    tl=uniqrows(sort(t,1));gl=uniqrows(sort(g,1));gl=select(gl,(gl:>0):&(gl:>tl[1]))
    if(control=="notyettreated" & !any(g:==0) & latest<.)gl=select(gl,gl:<latest)
    if(id!="") {
        if(!st_isnumvar(id))_cdid_op_fail("A supplied ordinary identity must be numeric after finite-row screening.",498)
        ord=order((iv,t),(1,2));info=_cdid_ps_opre(iv[ord],1)
        if(rows(uniqrows((iv,g)))!=rows(info))_cdid_op_fail("The cohort must be constant within supplied identity.",498)
        if(rows(uniqrows((iv,t)))!=rows(rr))_cdid_op_fail("A supplied identity must be unique by time.",498)
        for(j=1;j<=cols(cn);j++)if(cn[j]!=id) {
            labs=(st_isnumvar(cn[j]) ? strofreal(st_data(rr,cn[j]):+0,"%21x") : st_sdata(rr,cn[j]))
            for(n=1;n<=rows(info);n++) {
                ix=ord[|info[n,1]\info[n,2]|]
                if(any(labs[ix]:!=labs[ix[1]]))_cdid_op_fail("Time-varying cluster variables are not supported.",498)
            }
        }
    }
    if(panel & !allow) {
        keep=J(rows(rr),1,0)
        for(j=1;j<=rows(info);j++) {
            ix=ord[|info[j,1]\info[j,2]|]
            if(rows(ix)==rows(tl))keep[ix]=J(rows(ix),1,1)
        }
        O.rowmap[select(rr,keep:==0),2]=J(sum(keep:==0),1,5)
        rr=select(rr,keep);t=select(t,keep);g=select(g,keep);w=select(w,keep);iv=select(iv,keep)
        if(!rows(rr))_cdid_op_fail("All observations dropped while balancing the panel.",498)
    }
    if(!rows(gl))_cdid_op_fail("No valid treated groups.",498)
    uq=uniqrows(sort(g,1))
    for(j=1;j<=rows(uq);j++)if(sum(g:==uq[j])/rows(tl)<cols(xnames)+5) {
        printf("Small ordinary group %g: %g observations per calendar period.\n",uq[j],sum(g:==uq[j])/rows(tl))
        if(uq[j]==0 & control=="nevertreated")_cdid_op_fail("The never-treated group is too small; use notyettreated explicitly if appropriate.",498)
    }
    if(rows(tl)<2)_cdid_op_fail("At least two periods are required.",498)
    if(max(gl)<=tl[2])_cdid_op_fail("There are no pre-treatment periods available for the conditional test.",498)
    if(allow)_cdid_op_fail("The conditional test does not support allowunbalanced.",498)
    ord=(panel ? order((iv,t),(1,2)) : (1::rows(rr)))
    rr=rr[ord];t=t[ord];g=g[ord];w=w[ord];iv=iv[ord]
    if(!panel & id=="")iv=(1::rows(rr))
    O.rows=rr;O.ids=iv;O.glist=gl;O.panel=panel;O.allow=allow;O.intercept=intercept
    O.inputN=st_nobs();O.requestedN=sum(O.rowmap[,2]:!=0);O.maxentries=cap
    O.formula=formula;O.cluster=clusters;O.idname=id
    O.callback_name=callback_name;O.payload_name=payload_name
    O.recipe.panel=panel;O.recipe.intercept_only=(formula=="" & intercept)
    O.recipe.method=method;O.recipe.control=control;O.recipe.d1_policy=(method=="custom" ? "unaffected" : "corrected")
    O.recipe.row_period=t;O.recipe.row_y=st_data(rr,y);O.recipe.row_w=w;O.recipe.row_g=g
    O.recipe.tlist=tl;O.recipe.cells=cdid_oc_pre_cells(gl,tl)
    O.recipe.key=(panel | id=="" ? J(0,1,.) : iv)
    if(panel) {
        O.recipe.row_unit=cdid_oc_first_positions(iv)
        O.gridrows=_cdid_op_ix(t:==tl[1])
    }
    else {O.recipe.row_unit=J(0,1,.);O.gridrows=(1::rows(rr));}
    if(panel & id!="")cn=select(cn,cn:!=id)
    if(cols(cn)>1)_cdid_op_fail("The source true-RCS multiplier supports at most one requested cluster variable after .rowid replacement.",498)
    if(!cols(cn)) {
        O.recipe.cmap=(1::rows(O.gridrows));O.cluster_raw=O.recipe.cmap
        O.cluster_labels=strofreal(O.cluster_raw:+0,"%21x")
    }
    else {
        cl=rr[O.gridrows]
        if(st_isnumvar(cn[1])) {O.cluster_raw=st_data(cl,cn[1]);labs=strofreal(O.cluster_raw:+0,"%21x");}
        else {O.cluster_raw=J(0,1,.);labs=st_sdata(cl,cn[1]);}
        O.recipe.cmap=_cdid_op_cmap(labs,discard);O.cluster_labels=discard
    }
    st_store(rr,ordervar,(1::rows(rr)));st_store(rr,groupvar,g)
    p=crexternal(name);*p=O
    st_local("op_rows",strofreal(rows(rr)));st_local("op_varying_w",strofreal(O.varying_w))
}
/* This invokes the ONE common native FV encoding primitive. The ordinary consumer
   supplies its own fixed global versus first-calendar policy, not a chained union. */
void cdid_op_encode(string scalar name,real scalar grid,string scalar generated,
    string scalar expanded,real rowvector omitted)
{
    pointer scalar p
    struct cdid_op_input scalar O
    struct cdid_reference_inputs scalar S
    struct cdid_ref_design scalar D
    p=findexternal(name);if(p==NULL)_cdid_op_fail("Intake object is absent.",498)
    O=*p;D.keys=(1::st_nobs());D.intercept=O.intercept;D.encoded=0;D.formula=O.formula
    S.designs=D;S.n_slots=1;S.max_entries=O.maxentries;S.retained_entries=0
    cdid_ref_encoded(S,1,generated,expanded,omitted)
    if(grid) {O.recipe.Z=S.designs[1].value;O.zcols=S.designs[1].columns;}
    else {O.recipe.X=S.designs[1].value;O.xcols=S.designs[1].columns;}
    *p=O
}
void cdid_op_gridmark(string scalar name,string scalar variable)
{
    pointer scalar p
    struct cdid_op_input scalar O
    p=findexternal(name);O=*p
    st_store(O.gridrows,variable,J(rows(O.gridrows),1,1))
}
/* No provider tensor is retained. Scalar/draw/family metadata is copied to a
   conditional-specific durable owner only when store() was explicitly requested. */
void cdid_op_execute(string scalar name,real scalar B,real scalar alpha,real scalar chunk,
    real scalar bytes,string scalar rng,string scalar signsmat,string scalar posting,
    string scalar store,real scalar maxmatrix,string scalar stats,string scalar draws,string scalar family,
    string scalar mult)
{
    pointer scalar p,q,p_fn,p_pl
    struct cdid_op_input scalar O
    struct cdid_op_state scalar R
    struct cdid_oc_problem scalar P
    struct cdid_oc_model colvector models
    struct cdid_co_process scalar co
    struct cdid_cpt_budget scalar budget
    struct cdid_oc_mult scalar acc
    struct cdid_oc_result scalar result
    struct cdid_public_sign_stream scalar stream
    struct cdid_public_sign_batch scalar batch
    struct cdid_infer_budget scalar ib
    real scalar cap,c,N,M,K,C,unavail
    real colvector bad_m
    real matrix fixed
    p=findexternal(name);if(p==NULL)_cdid_op_fail("Intake object is absent.",498)
    O=*p;cap=floor(bytes/8)
    if(store!="") {
        if(findexternal(store)!=NULL)_cdid_op_fail("Requested conditional state object already exists.",498)
    }
    budget=cdid_cpt_budget_create(2^31-1,2^31-1,2^31-1,B,2^31-1,32767,32767,
        chunk,2^31-1,chunk,max((chunk,cap)),chunk,chunk,cap,floor(bytes),.)
    fixed=(signsmat=="" ? J(0,0,.) : st_matrix(signsmat))
    if(O.recipe.method=="custom") {
        _cdid_op_validate_symbols(O.callback_name,O.payload_name)
        p_fn=findexternal(O.callback_name+"()")
        p_pl=(O.payload_name!="" ? findexternal(O.payload_name) : NULL)
        if(cols(O.recipe.X)<1 | cols(O.recipe.Z)<1)
            _cdid_op_fail("Ordinary conditional test requires at least one covariate or intercept column.",498)
        N=rows(O.recipe.Z);M=N;K=rows(O.recipe.cells);C=max(O.recipe.cmap)
        if(posting=="matrix" & max((B,K,M))>maxmatrix)
            _cdid_op_fail("Requested matrix results exceed matsize; request posting(state) with store() before drawing.",908)
        ib=cdid_infer_budget_create(N,C,K,B,chunk,max((chunk,cap)),chunk,max((chunk,cap)),cap,floor(bytes))
        /* Validate signs/budget before any callback or caller RNG consumption. */
        stream=cdid_public_rng_start(B,C,ib,rng,fixed)
        co=cdid_co_prepare(O.recipe,budget,O.xcols,p_fn,p_pl)
        R.family=J(K,5,.);R.family_reasons=J(K,1,"")
        for(c=1;c<=K;c++) {
            unavail=any(co.cell_rc[,c]:!=0)
            R.family[c,]=(c,O.recipe.cells[c,1],O.recipe.cells[c,2],O.recipe.cells[c,3],unavail)
            if(unavail) {
                bad_m=_cdid_op_ix(co.cell_rc[,c]:!=0)
                R.family_reasons[c]=co.cell_reason[bad_m[1],c]
            }
        }
        R.stats=(co.N,co.M,co.K,C,B,alpha,.,.,.,co.M,rows(O.rows),O.requestedN,O.inputN)
        R.status="unavailable";R.reason="";R.draws=J(0,1,.)
        R.method="custom";R.formula=O.formula;R.membership=co.membership_law
        R.channel="supplied by the method(custom) callback";R.rng=rng
        R.policy="data preparation as in did 2.5.1; weights as supplied; unweighted grid of all units; user-supplied estimator"
        R.multiplicity=J(co.M,1,1);R.grid_first=(1::co.M);R.retained_rows=O.rows
        R.rowmap=O.rowmap;R.cluster_labels=O.cluster_labels;R.xcols=O.xcols;R.zcols=O.zcols
        R.cluster_positions=O.recipe.cmap;R.cluster_values=O.cluster_raw
        R.preprocess=(O.rows,O.ids,O.recipe.row_period,O.recipe.row_g,O.recipe.row_w,
            (O.panel ? O.recipe.row_unit : (1::rows(O.rows))))
        R.callback_name=O.callback_name;R.payload_name=O.payload_name;R.n_callbacks=co.n_callbacks
        if(!co.available) {
            R.reason="At least one declared pre cell is unavailable; the complete family is retained."
        }
        else {
            acc=cdid_oc_xp_mult_start(co.process,B,budget)
            while(stream.cursor<stream.B) {
                batch=cdid_public_rng_next(stream)
                cdid_oc_xp_mult_add(co.process,acc,batch)
            }
            cdid_public_rng_complete(stream)
            result=cdid_oc_xp_finish(co.process,acc,alpha)
            R.status="available";R.draws=result.CvMb
            R.stats[7..9]=(result.CvM,result.CvMcval,result.CvMpval)
        }
    }
    else {
        P=cdid_oc_create(O.recipe,budget)
        if(posting=="matrix" & max((B,P.K,P.grid.n_distinct))>maxmatrix)
            _cdid_op_fail("Requested matrix results exceed matsize; request posting(state) with store() before drawing.",908)
        ib=cdid_infer_budget_create(P.N,P.C,P.K,B,chunk,max((chunk,cap)),chunk,max((chunk,cap)),cap,floor(bytes))
        /* Validate signs/budget before any nuisance work or caller RNG consumption. */
        stream=cdid_public_rng_start(B,P.C,ib,rng,fixed)
        models=cdid_oc_prepare_all(P)
        R.family=cdid_oc_family(models);R.family_reasons=cdid_oc_family_reasons(models)
        R.stats=(P.N,P.M,P.K,P.C,B,alpha,.,.,.,P.grid.n_distinct,rows(O.rows),O.requestedN,O.inputN)
        R.status="unavailable";R.reason="";R.draws=J(0,1,.)
        R.method=O.recipe.method;R.formula=O.formula;R.membership=P.membership_law
        R.channel=(P.d1_policy=="corrected" ? "corrected doubly robust influence function for repeated cross-sections" : "influence function of the ordinary estimator");R.rng=rng
        R.policy="data preparation as in did 2.5.1; weights as supplied; unweighted grid of all units; intercept-only formula and observation identifiers by default"
        R.multiplicity=P.grid.mult;R.grid_first=P.grid.first;R.retained_rows=O.rows
        R.rowmap=O.rowmap;R.cluster_labels=O.cluster_labels;R.xcols=O.xcols;R.zcols=O.zcols
        R.cluster_positions=O.recipe.cmap;R.cluster_values=O.cluster_raw
        R.preprocess=(O.rows,O.ids,O.recipe.row_period,O.recipe.row_g,O.recipe.row_w,
            (O.panel ? O.recipe.row_unit : (1::rows(O.rows))))
        R.callback_name="";R.payload_name="";R.n_callbacks=0
        if(any(R.family[,5]:!=0)) {
            R.reason="At least one declared pre cell is unavailable; the complete family is retained."
        }
        else {
            /* Exact high-level core composition, retaining the typed model array. */
            for(c=1;c<=P.K;c++) {
                cdid_oc_model_grid(P,models[c])
                cdid_oc_provider(P,models[c])
            }
            acc=cdid_oc_mult_start(P,models,B)
            while(stream.cursor<stream.B) {
                batch=cdid_public_rng_next(stream)
                cdid_oc_mult_add(P,models,acc,batch)
            }
            cdid_public_rng_complete(stream)
            result=cdid_oc_finish(P,models,acc,alpha)
            R.status="available";R.draws=result.CvMb
            R.stats[7..9]=(result.CvM,result.CvMcval,result.CvMpval)
            R.channel=result.channel
        }
    }
    st_matrix(stats,R.stats)
    if(posting=="matrix") {
        if(rows(R.draws))st_matrix(draws,R.draws)
        st_matrix(family,R.family);st_matrix(mult,(R.grid_first,R.multiplicity))
    }
    st_local("op_status",R.status);st_local("op_reason",R.reason)
    st_local("op_xcols",invtokens(R.xcols));st_local("op_zcols",invtokens(R.zcols))
    st_local("op_channel",R.channel);st_local("op_membership",R.membership)
    st_local("op_family_reasons",invtokens(R.family_reasons'," | "))
    st_local("op_n_callbacks",strofreal(R.n_callbacks))
    if(store!="") {q=crexternal(store);*q=R;}
}
void cdid_op_drop(string scalar name)
{
    if(findexternal(name)!=NULL)rmexternal(name)
}
end
