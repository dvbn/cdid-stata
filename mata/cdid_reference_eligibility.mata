*! Cohort-specific reference eligibility for time-varying covariates; StataNow19.5.
*! When a time-varying covariate required at a target cohort's reference date b_g
*! is not supplied for some units, missingreference(exclude) (the default)
*! restricts that cohort's comparisons to the units E_g whose required raw values
*! at b_g are known. Units keep their identities, rows, clusters and full
*! influence-function rows, other cohorts still use them, and the common
*! normalizer N is unchanged. Values are never read from another date.
version 19.5
mata:
struct cdid_mref_state {
    string scalar policy
    real scalar pending, active, anticipation, N
    real colvector target_cohorts, calendar
    string rowvector variables
    real colvector index
    real matrix excluded
    real matrix targets, contexts
    string colvector ctx_by_var, ctx_ids_head, tgt_ids_head
    real matrix occ_sub
    real matrix aux
    string colvector aux_K, aux_reason
    real colvector aux_missing
    string colvector aux_missing_reason
}
real colvector _cdid_mref_ix(real matrix x)
{
    return(vec(selectindex(x)))
}
pointer(struct cdid_mref_state scalar) scalar _cdid_mref_p()
{
    return(findexternal("__cdid_mref_state"))
}
void _cdid_mref_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
/* policy "" removes the request; "exclude" or "error" starts a fresh one. */
void cdid_mref_request(string scalar policy)
{
    pointer(struct cdid_mref_state scalar) scalar p
    struct cdid_mref_state scalar S
    string scalar pol
    pol=strlower(strtrim(policy))
    if(findexternal("__cdid_mref_state")!=NULL) rmexternal("__cdid_mref_state")
    if(pol=="") return
    if(pol!="exclude" & pol!="error") _cdid_mref_fail("missingreference() must be exclude or error.",198)
    S.policy=pol;S.pending=S.active=0;S.anticipation=0;S.N=0
    S.target_cohorts=S.calendar=S.index=J(0,1,.)
    S.variables=J(1,0,"");S.excluded=J(0,0,.)
    S.targets=J(0,9,.);S.contexts=J(0,3,.)
    S.ctx_by_var=S.ctx_ids_head=S.tgt_ids_head=J(0,1,"")
    S.occ_sub=J(0,5,.);S.aux=J(0,4,.);S.aux_K=S.aux_reason=J(0,1,"")
    S.aux_missing=J(0,1,.);S.aux_missing_reason=J(0,1,"")
    p=crexternal("__cdid_mref_state")
    *p=S
}
string scalar cdid_mref_policy()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    return(p==NULL ? "" : p->policy)
}
real scalar cdid_mref_active()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return(0)
    return(p->active)
}
real scalar cdid_mref_pending()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return(0)
    return(p->pending)
}
void _cdid_mref_set_pending()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p!=NULL) p->pending=1
}
void cdid_mref_clear_pending()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p!=NULL) p->pending=0
}
/* Record the request's target cohorts, calendar and anticipation. */
void cdid_mref_set_request(real colvector targets, real colvector calendar, real scalar ant)
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return
    p->target_cohorts=uniqrows(sort(select(targets,targets:<.),1))
    p->calendar=calendar;p->anticipation=ant
}
/* Restricted routes check this before reading any restricted context. */
real scalar cdid_mref_any()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return(0)
    return(p->active | rows(p->aux_missing)>0)
}
void cdid_mref_require_route(string scalar what)
{
    if(!cdid_mref_any()) return
    _cdid_mref_fail("Restricting cohort comparisons to units with known reference covariates (missingreference(exclude)) is not available "+
        what+". Some units lack a time-varying covariate at a date that a comparison requires. Supply those covariates on their "+
        "unit-date rows (the outcome may be missing on such a row); missingreference(error) lists the missing values.",459)
}
/* Requirement scans and catalog: context number (1..n) of each target cohort value; 0 when unrestricted. */
real colvector cdid_mref_ctx_of(real colvector cohort_values)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector out,hit
    real scalar i,b,bi
    out=J(rows(cohort_values),1,0)
    p=_cdid_mref_p()
    if(p==NULL) return(out)
    if(!p->active) return(out)
    for(i=1;i<=rows(cohort_values);i++) {
        if(cohort_values[i]>=.) continue
        b=_cdid_pr_reference(p->calendar,cohort_values[i],p->anticipation)
        if(b>=.) continue
        bi=_cdid_population_lookup(b,p->calendar)
        hit=_cdid_mref_ix(p->index:==bi)
        if(rows(hit)) out[i]=hit[1]
    }
    return(out)
}
/* Calendar index of the eligibility context of each target cohort value (0 when unrestricted). */
real colvector cdid_mref_ctx_date_of(real colvector cohort_values)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector c,out
    real scalar i
    c=cdid_mref_ctx_of(cohort_values);out=J(rows(c),1,0)
    p=_cdid_mref_p()
    for(i=1;i<=rows(c);i++) if(c[i]>0) out[i]=p->index[c[i]]
    return(out)
}
real colvector cdid_mref_index()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return(J(0,1,.))
    if(!p->active) return(J(0,1,.))
    return(p->index)
}
/* N x n_ctx: 1 = excluded from that context. */
real matrix cdid_mref_excluded(real scalar N)
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return(J(N,0,.))
    if(!p->active) return(J(N,0,.))
    if(rows(p->excluded)!=N) _cdid_mref_fail("Reference eligibility contexts do not match the unit roster.",498)
    return(p->excluded)
}
/* Catalog occupancy of the restricted contexts in index space. tg holds the
   catalog's sorted target cohort indices; tg_ctx receives the context of each
   target (0 = owner cells) and sub the excluded units' (context, cohort index,
   s index, t index, count). Returns 1 when restricted contexts exist. */
real scalar cdid_mref_catalog(real colvector calendar, real colvector cohorts, real scalar nt,
    real colvector tg, real colvector tg_ctx, real matrix sub)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector k,s,t,keep
    tg_ctx=J(rows(tg),1,0);sub=J(0,5,.)
    p=_cdid_mref_p()
    if(p==NULL) return(0)
    if(!p->active) return(0)
    tg_ctx=cdid_mref_ctx_of(cohorts[tg])
    if(rows(p->occ_sub)) {
        k=_cdid_cohort_index_vec(cohorts,nt,p->occ_sub[,2],J(rows(p->occ_sub),1,1))
        s=_cdid_population_lookup(p->occ_sub[,3],calendar)
        t=_cdid_population_lookup(p->occ_sub[,4],calendar)
        if(any(k:==0) | any(s:==0) | any(t:==0))
            _cdid_mref_fail("Eligibility occupancy lies outside the catalog cohort and calendar dictionaries.",498)
        sub=(p->occ_sub[,1],k,s,t,p->occ_sub[,5])
        keep=_cdid_mref_ix(sub[,5]:>0)
        sub=(rows(keep) ? sub[keep,] : J(0,5,.))
    }
    return(1)
}
/* Auxiliary comparisons made unavailable because their own earlier-date
   covariate is missing for an eligible unit (q indices into the request). */
real colvector cdid_mref_aux_missing(real scalar Q)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector out
    real scalar i
    out=J(Q,1,0)
    p=_cdid_mref_p()
    if(p==NULL) return(out)
    for(i=1;i<=rows(p->aux_missing);i++) if(p->aux_missing[i]>=1 & p->aux_missing[i]<=Q) out[p->aux_missing[i]]=1
    return(out)
}
string scalar cdid_mref_aux_missing_reason(real scalar q)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector hit
    p=_cdid_mref_p()
    if(p==NULL) return("")
    hit=_cdid_mref_ix(p->aux_missing:==q)
    return(rows(hit) ? p->aux_missing_reason[hit[1]] : "")
}
void _cdid_mref_aux_missing_add(real scalar q, string scalar reason)
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return
    if(rows(p->aux_missing)) if(any(p->aux_missing:==q)) return
    p->aux_missing=p->aux_missing\q;p->aux_missing_reason=p->aux_missing_reason\reason
}
/* Record one unavailable auxiliary comparison (deduplicated by q). */
void cdid_mref_aux_note(real scalar q, real scalar g, real scalar s, real scalar t, string scalar K, string scalar reason)
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return
    if(rows(p->aux)) if(any(p->aux[,1]:==q)) return
    p->aux=p->aux\(q,g,s,t);p->aux_K=p->aux_K\K;p->aux_reason=p->aux_reason\reason
}
/* Cross-pair fits have no same-pair pretrend comparisons: their notes are dropped. */
void cdid_mref_aux_clear()
{
    pointer(struct cdid_mref_state scalar) scalar p
    p=_cdid_mref_p()
    if(p==NULL) return
    p->aux=J(0,4,.);p->aux_K=p->aux_reason=J(0,1,"")
}
/* Units of a restricted target cohort outside its own context (aggregation masses). */
real colvector cdid_mref_share_excluded(struct cdid_population scalar P)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real colvector out,c
    real scalar i
    out=J(P.N,1,0)
    p=_cdid_mref_p()
    if(p==NULL) return(out)
    if(!p->active) return(out)
    if(rows(p->excluded)!=P.N) _cdid_mref_fail("Reference eligibility contexts do not match the unit roster.",498)
    c=cdid_mref_ctx_of(P.cohort)
    for(i=1;i<=P.N;i++) if(c[i]>0 & P.cohort_known[i]==1) out[i]=p->excluded[i,c[i]]
    return(out)
}
/* Stored results: the aggregation column of a restricted fit (empty otherwise). */
real colvector cdid_mref_state_excluded(struct cdid_result_state scalar P)
{
    real colvector ix
    ix=_cdid_mref_ix(P.unit_columns:=="reference_share_excluded")
    if(length(ix)!=1) return(J(0,1,.))
    if(rows(P.unit_metadata)==0) return(J(0,1,.))
    return(P.unit_metadata[,ix[1]])
}
string scalar _cdid_mref_ids(struct cdid_population scalar P, real colvector pos)
{
    string scalar out
    real scalar j,n
    n=min((5,rows(pos)));out=""
    for(j=1;j<=n;j++) out=out+(j>1 ? ", " : "")+
        (P.ids.string_id ? P.ids.string_id_value[pos[j]] : _cdid_pr_num(P.ids.numeric_id[pos[j]]))
    if(rows(pos)>5) out=out+", ..."
    return(out)
}
string scalar _cdid_mref_id1(struct cdid_population scalar P, real scalar u)
{
    return(P.ids.string_id ? P.ids.string_id_value[u] : _cdid_pr_num(P.ids.numeric_id[u]))
}
/* Build the restricted contexts of the affected reference dates (calendar
   indices) from the raw reference values of the dated formulas' variables. */
void _cdid_mref_build(struct cdid_reference_inputs scalar S, real colvector affected, real rowvector vars)
{
    pointer(struct cdid_mref_state scalar) scalar p
    struct cdid_population scalar P
    real colvector ctx,miss,ex,own,tg,tc,obs,rr,ok,ei
    real matrix raw,excl,TG,CX,OS,dd
    real scalar c,j,i,n,g,u,d1,d2,ng,nT,b,bi,ne
    string scalar byvar
    string colvector byv,idh,tgh
    transmorphic scalar cells
    real rowvector key
    p=_cdid_mref_p()
    if(p==NULL) _cdid_mref_fail("Reference eligibility requires an active request.",498)
    P=*S.population
    if(rows(p->calendar)!=rows(P.calendar)) p->calendar=P.calendar
    p->N=P.N;p->variables=S.raw_names[vars]
    ctx=uniqrows(sort(select(affected,affected:>0),1))
    excl=J(P.N,rows(ctx),0);byv=idh=J(rows(ctx),1,"");CX=J(rows(ctx),3,.)
    for(c=1;c<=rows(ctx);c++) {
        raw=cdid_ref_raw_query(S,ctx[c]*(S.N+1):+(1::S.N),vars)
        miss=(rowmissing(raw):>0);excl[,c]=miss
        byvar=""
        for(j=1;j<=cols(vars);j++) byvar=byvar+(j>1 ? ", " : "")+S.raw_names[vars[j]]+"="+strofreal(colmissing(raw[,j]))
        byv[c]=byvar
        ex=_cdid_mref_ix(miss)
        if(rows(ex)) idh[c]=_cdid_mref_ids(P,ex)
        CX[c,]=(c,P.calendar[ctx[c]],sum(miss))
    }
    p->ctx_by_var=byv;p->ctx_ids_head=idh;p->contexts=CX
    /* Only contexts that exclude at least one unit restrict comparisons. */
    ok=(rows(ctx) ? _cdid_mref_ix(CX[,3]:>0) : J(0,1,.))
    if(rows(ok)) {
        p->index=ctx[ok]
        p->excluded=excl[,ok']
    }
    else {
        p->index=J(0,1,.);p->excluded=J(P.N,0,.)
    }
    p->active=(rows(ok)>0)
    tg=p->target_cohorts
    if(rows(tg)==0) tg=uniqrows(sort(select(P.cohort,(P.cohort:<.):&(P.cohort_known:==1)),1))
    ng=rows(tg);TG=J(ng,9,.);tgh=J(ng,1,"")
    for(i=1;i<=ng;i++) {
        g=tg[i];b=_cdid_pr_reference(P.calendar,g,p->anticipation)
        bi=0
        if(b<.) bi=_cdid_population_lookup(b,P.calendar)
        tc=J(0,1,.)
        if(bi>0 & rows(ctx)) tc=_cdid_mref_ix(ctx:==bi)
        own=_cdid_mref_ix((P.cohort:==g):&(P.cohort_known:==1))
        if(rows(tc)) {
            ex=_cdid_mref_ix(excl[,tc[1]])
            ne=0
            if(rows(own)) ne=sum(excl[own,tc[1]])
            TG[i,]=(g,b,tc[1],rows(own),rows(own)-ne,ne,0,0,(rows(ex)>0))
            if(rows(ex)) {
                TG[i,7]=sum((P.cohort[ex]:>=.):&(P.cohort_known[ex]:==1))
                TG[i,8]=sum((P.cohort[ex]:<.):&(P.cohort_known[ex]:==1):&(P.cohort[ex]:!=g))
            }
            if(ne) {
                ei=_cdid_mref_ix(excl[own,tc[1]])
                tgh[i]=_cdid_mref_ids(P,own[ei])
            }
        }
        else TG[i,]=(g,b,.,rows(own),rows(own),0,0,0,0)
    }
    p->targets=TG;p->tgt_ids_head=tgh
    /* Pair occupancy of the excluded units with observed outcomes at both dates. */
    OS=J(0,5,.);nT=rows(P.calendar)
    if(rows(ok)) {
        cells=asarray_create("real",4)
        for(c=1;c<=rows(ok);c++) {
            ex=_cdid_mref_ix(excl[,ok[c]]:&(P.A:==1))
            for(j=1;j<=rows(ex);j++) {
                u=ex[j]
                rr=cdid_population_rows(P,J(nT,1,u),P.calendar)
                obs=_cdid_mref_ix(rr:>0)
                if(rows(obs)==0) continue
                obs=select(obs,P.y[rr[obs]]:<.)
                for(d1=1;d1<rows(obs);d1++) for(d2=d1+1;d2<=rows(obs);d2++) {
                    key=(c,(P.cohort[u]<. ? P.cohort[u] : maxdouble()),obs[d1],obs[d2])
                    if(asarray_contains(cells,key)) asarray(cells,key,asarray(cells,key)+1)
                    else asarray(cells,key,1)
                }
            }
        }
        dd=asarray_keys(cells)
        n=rows(dd)
        if(n) {
            OS=J(n,5,.)
            for(i=1;i<=n;i++) OS[i,]=(dd[i,1],(dd[i,2]==maxdouble() ? . : dd[i,2]),P.calendar[dd[i,3]],P.calendar[dd[i,4]],asarray(cells,dd[i,]))
            OS=sort(OS,(1,3,4,2))
        }
    }
    p->occ_sub=OS
}
/* Requirement scan: the dated formulas' raw variables (dated or unit attributes). */
real rowvector _cdid_mref_vars(struct cdid_reference_inputs scalar S)
{
    real rowvector v
    real scalar r,j
    v=J(1,0,.)
    for(r=1;r<=rows(S.roles);r++) {
        if(!S.roles[r].active | !S.roles[r].dated) continue
        for(j=1;j<=cols(S.roles[r].raw_index);j++) {
            if(S.kind[S.roles[r].raw_index[j]]==3) continue
            if(cols(v)) if(any(v:==S.roles[r].raw_index[j])) continue
            v=v,S.roles[r].raw_index[j]
        }
    }
    return(sort(v',1)')
}
/* The note shown with restricted fits and the exclusion warning. */
void cdid_mref_display()
{
    pointer(struct cdid_mref_state scalar) scalar p
    real scalar i,n,shown
    real colvector r
    string scalar txt
    p=_cdid_mref_p()
    if(p==NULL) return
    if(p->policy!="exclude") return
    if(p->active) {
        r=_cdid_mref_ix(p->targets[,9]:==1);n=rows(r);shown=min((5,n));txt=""
        for(i=1;i<=shown;i++) txt=txt+(i>1 ? ", " : "")+"cohort "+_cdid_pr_num(p->targets[r[i],1])+" (date "+
            _cdid_pr_num(p->targets[r[i],2])+"): "+strofreal(p->targets[r[i],6])+" treated, "+
            strofreal(p->targets[r[i],7])+" never-treated, "+strofreal(p->targets[r[i],8])+" from other cohorts"
        if(n>shown) txt=txt+" and additional cohorts"
        printf("{txt}{p 0 9 2}Warning: Excluded units with missing reference covariates from the affected cohort comparisons. "+
            "Counts across the supplied units at each reference date, including units outside the comparison groups: %s. "+
            "See e(reference_eligibility) for details. If you know a missing covariate, supply it on its unit-date row, "+
            "with a missing value for a missing outcome.{p_end}\n",txt)
    }
    else if(rows(p->aux)) {
        printf("{txt}{p 0 9 2}Warning: Some pretrend comparisons are unavailable because required reference covariates are missing. "+
            "See e(reference_aux_unavailable) for details.{p_end}\n")
    }
}
/* Output note of a restricted fit (also on replay), read from e(). */
void cdid_mref_note()
{
    real matrix T,A
    real colvector r
    string scalar txt
    real scalar i
    if(st_global("e(missing_reference)")!="exclude") return
    if(st_global("e(reference_status)")=="restricted") {
        T=st_matrix("e(reference_eligibility)")
        if(rows(T) & cols(T)==8) {
            r=_cdid_mref_ix(T[,8]:==1);txt=""
            for(i=1;i<=rows(r);i++) txt=txt+(i>1 ? ", " : "")+_cdid_pr_num(T[r[i],1])
            printf("{txt}Restricted comparisons for cohort(s) %s. See e(reference_eligibility).\n",txt)
            r=_cdid_mref_ix(T[,5]:>0);txt=""
            for(i=1;i<=rows(r);i++) txt=txt+(i>1 ? ", " : "")+"cohort "+_cdid_pr_num(T[r[i],1])+" cover "+
                strofreal(T[r[i],4])+" of "+strofreal(T[r[i],3])+" treated units"
            if(rows(r)) printf("{txt}Effects for %s.\n",txt)
        }
    }
    A=st_matrix("e(reference_aux_unavailable)")
    if(rows(A)) printf("{txt}  %f auxiliary pretrend comparison(s) unavailable under missingreference(exclude).\n",rows(A))
}
/* The eligibility record travels with the stored results (saving() and cdid_load). */
void cdid_mref_state_record(struct cdid_result_state scalar S)
{
    pointer(struct cdid_mref_state scalar) scalar p
    real matrix T,Cx,A
    real scalar i,j
    string scalar status
    p=_cdid_mref_p()
    if(p==NULL) return
    if(p->policy!="exclude" | (!p->active & rows(p->aux)==0)) return
    status=(p->active ? "restricted" : "auxiliary_unavailable")
    S.text_settings=S.text_settings\("missing_reference",p->policy)\("reference_status",status)
    if(p->active) {
        S.text_settings=S.text_settings\("reference_variables",invtokens(p->variables))
        T=p->targets[,(1,2,4,5,6,7,8,9)]
        Cx=p->contexts[_cdid_mref_ix(p->contexts[,3]:>0),(2,3)]
    }
    else {T=J(0,8,.);Cx=J(0,2,.);}
    A=(rows(p->aux) ? sort(p->aux[,(2,3,4)],(1,2,3)) : J(0,3,.))
    S.numeric_settings=S.numeric_settings\rows(T)\rows(Cx)\rows(A)
    S.numeric_setting_names=S.numeric_setting_names\"reference_eligibility.rows"\"reference_contexts.rows"\"reference_aux_unavailable.rows"
    for(i=1;i<=rows(T);i++) for(j=1;j<=8;j++) {
        S.numeric_settings=S.numeric_settings\T[i,j]
        S.numeric_setting_names=S.numeric_setting_names\("reference_eligibility."+strofreal(i)+"."+strofreal(j))
    }
    for(i=1;i<=rows(Cx);i++) for(j=1;j<=2;j++) {
        S.numeric_settings=S.numeric_settings\Cx[i,j]
        S.numeric_setting_names=S.numeric_setting_names\("reference_contexts."+strofreal(i)+"."+strofreal(j))
    }
    for(i=1;i<=rows(A);i++) for(j=1;j<=3;j++) {
        S.numeric_settings=S.numeric_settings\A[i,j]
        S.numeric_setting_names=S.numeric_setting_names\("reference_aux_unavailable."+strofreal(i)+"."+strofreal(j))
    }
}
real matrix _cdid_mref_state_matrix(struct cdid_result_state scalar P, string scalar name, real scalar nc)
{
    real colvector h
    real matrix M
    real scalar n,i,j
    h=_cdid_mref_ix(P.numeric_setting_names:==name+".rows")
    if(rows(h)!=1) return(J(0,nc,.))
    n=P.numeric_settings[h[1]];M=J(n,nc,.)
    for(i=1;i<=n;i++) for(j=1;j<=nc;j++) {
        h=_cdid_mref_ix(P.numeric_setting_names:==name+"."+strofreal(i)+"."+strofreal(j))
        if(rows(h)==1) M[i,j]=P.numeric_settings[h[1]]
    }
    return(M)
}
string scalar _cdid_mref_state_text(struct cdid_result_state scalar P, string scalar name)
{
    real colvector h
    if(rows(P.text_settings)==0) return("")
    h=_cdid_mref_ix(P.text_settings[,1]:==name)
    return(rows(h) ? P.text_settings[h[1],2] : "")
}
void _cdid_mref_post_matrix(string scalar name, real matrix M, string colvector cn, string scalar rowprefix, real scalar rowcol)
{
    string matrix rs
    real scalar i
    if(rows(M)==0) return
    rs=J(rows(M),2,"")
    for(i=1;i<=rows(M);i++) rs[i,2]=(rowcol ? rowprefix+_cdid_pr_num(M[i,rowcol]) : rowprefix+strofreal(i))
    st_matrix(name,M);st_matrixrowstripe(name,rs);st_matrixcolstripe(name,(J(rows(cn),1,""),cn))
}
/* cdid_load: repost the eligibility record of stored chained results. */
void cdid_mref_load_post(string scalar token)
{
    struct cdid_result_state scalar P
    real matrix T
    string scalar pol
    P=cdid_state_get(token)
    pol=_cdid_mref_state_text(P,"missing_reference")
    if(pol=="") return
    st_global("e(missing_reference)",pol)
    if(pol!="exclude") return
    st_global("e(reference_status)",_cdid_mref_state_text(P,"reference_status"))
    T=_cdid_mref_state_matrix(P,"reference_eligibility",8)
    if(rows(T)) {
        _cdid_mref_post_matrix("e(reference_eligibility)",T,("cohort"\"reference_date"\"n_treated"\"n_treated_eligible"\
            "n_treated_excluded"\"n_never_treated_excluded"\"n_other_cohort_excluded"\"restricted"),"g",1)
        _cdid_mref_post_matrix("e(reference_contexts)",_cdid_mref_state_matrix(P,"reference_contexts",2),
            ("reference_date"\"n_excluded"),"date",1)
        st_global("e(reference_variables)",_cdid_mref_state_text(P,"reference_variables"))
        st_global("e(reference_restricted)",invtokens(strofreal(select(T[,1],T[,8]:==1)',"%24.17g")))
    }
    _cdid_mref_post_matrix("e(reference_aux_unavailable)",_cdid_mref_state_matrix(P,"reference_aux_unavailable",3),
        ("cohort"\"s"\"t"),"r",0)
}
/* Stored results of a fit that ran with missingreference(exclude). */
void cdid_mref_post()
{
    pointer(struct cdid_mref_state scalar) scalar p
    real matrix T,Cx,A
    string matrix rs,cs
    real scalar i
    string scalar status
    p=_cdid_mref_p()
    if(p==NULL) return
    st_global("e(missing_reference)",p->policy)
    if(p->policy!="exclude") return
    status=(p->active ? "restricted" : (rows(p->aux) ? "auxiliary_unavailable" : "unrestricted"))
    st_global("e(reference_status)",status)
    if(p->active) {
        T=p->targets[,(1,2,4,5,6,7,8,9)]
        rs=J(rows(T),2,"");for(i=1;i<=rows(T);i++) rs[i,2]="g"+_cdid_pr_num(T[i,1])
        cs=(J(8,1,""),("cohort"\"reference_date"\"n_treated"\"n_treated_eligible"\"n_treated_excluded"\"n_never_treated_excluded"\"n_other_cohort_excluded"\"restricted"))
        st_matrix("e(reference_eligibility)",T);st_matrixrowstripe("e(reference_eligibility)",rs);st_matrixcolstripe("e(reference_eligibility)",cs)
        Cx=p->contexts[_cdid_mref_ix(p->contexts[,3]:>0),(2,3)]
        rs=J(rows(Cx),2,"");for(i=1;i<=rows(Cx);i++) rs[i,2]="date"+_cdid_pr_num(Cx[i,1])
        cs=(J(2,1,""),("reference_date"\"n_excluded"))
        st_matrix("e(reference_contexts)",Cx);st_matrixrowstripe("e(reference_contexts)",rs);st_matrixcolstripe("e(reference_contexts)",cs)
        st_global("e(reference_variables)",invtokens(p->variables))
        st_global("e(reference_restricted)",invtokens(strofreal(select(p->targets[,1],p->targets[,9]:==1)',"%24.17g")))
    }
    if(rows(p->aux)) {
        A=p->aux[,(2,3,4)];A=sort(A,(1,2,3))
        rs=J(rows(A),2,"");for(i=1;i<=rows(A);i++) rs[i,2]="r"+strofreal(i)
        cs=(J(3,1,""),("cohort"\"s"\"t"))
        st_matrix("e(reference_aux_unavailable)",A);st_matrixrowstripe("e(reference_aux_unavailable)",rs);st_matrixcolstripe("e(reference_aux_unavailable)",cs)
    }
}
end
