*! Shared reference-date/raw-history/FV preparation; StataNow 19.5.
*! GPL-3 adaptation of David Benatia, Christophe Bellego, Joel Cuerrier and Vincent Dortet-Bernadet.
*! R cdid /, chained_reference_inputs.R.
*! No nuisance fitting, alias reduction, carry-forward, or per-recipe full-F copies.
version 19.5
mata:
struct cdid_ref_role {
    real scalar active, intercept, dated, slot
    string scalar formula
    real rowvector raw_index
    real colvector keys, first_view, fit_keys
}
struct cdid_ref_design {
    string scalar formula
    real scalar intercept, dated, encoded
    string rowvector raw_names, columns, expanded
    real rowvector raw_index, omitted, column_map
    real colvector keys
    real matrix value
}
struct cdid_reference_inputs {
    pointer(struct cdid_population scalar) scalar population
    struct cdid_population_ids scalar original_ids
    real colvector calendar
    real scalar N, R, max_entries, retained_entries, n_slots, n_views, n_recipes
    string scalar route, method, fix_weights, weight_kind
    string rowvector raw_names
    real matrix data, attributes
    real rowvector has_data, has_frame, kind, first_change
    real colvector source_rows, unit_universe, weight_keys, weights, V, view_recipe
    real matrix recipes, view_meta, used_cohorts
    real colvector cohorts
    struct cdid_ref_role colvector roles
    struct cdid_ref_design colvector designs
    real colvector elig_index
    real matrix elig_excluded
}
void _cdid_ref_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_ref_ix(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
real scalar _cdid_ref_cost(real scalar a,real scalar b,real scalar cap,real scalar used)
{
    if(a<0 | b<0 | a>=. | b>=. | a!=floor(a) | b!=floor(b) |
        (a>0 & b>(2^53-1)/a) | a*b>cap-used)
        _cdid_ref_fail("Required reference storage exceeds the explicit numeric budget.",908)
    return(used+a*b)
}
real colvector cdid_ref_lookup(real colvector keys,real colvector table)
{
    return(_cdid_population_lookup(keys,table))
}
real rowvector _cdid_ref_names(string rowvector names,string scalar raw)
{
    string rowvector words
    real rowvector out
    real colvector ix
    real scalar j
    words=tokens(raw);out=J(1,cols(words),0)
    for(j=1;j<=cols(words);j++) {
        ix=_cdid_ref_ix(names:==words[j])
        if(rows(ix)!=1)_cdid_ref_fail("An active formula raw variable is absent or repeated.",498)
        out[j]=ix[1]
    }
    return(out)
}
/* kind1=dateless first supplied attribute,2=actual dated data,3=calendar date.
   All missing values are missing, never implicit zero or a carry-forward. */
/* O(1) lookup guard. The roster check is explicit at construction/activation;
   raw history/designs remain immutable borrowed inputs during a provider pass. */
void _cdid_ref_check_population(struct cdid_reference_inputs scalar S)
{
    struct cdid_population scalar P
    if(S.population==NULL)_cdid_ref_fail("The borrowed original population no longer exists.",498)
    P=*S.population
    if(P.N!=S.N | P.R!=S.R | P.ids.string_id!=S.original_ids.string_id |
        rows(P.calendar)!=rows(S.calendar))_cdid_ref_fail("The original population roster object changed.",498)
}
void cdid_ref_activate(struct cdid_reference_inputs scalar S)
{
    struct cdid_population scalar P
    _cdid_ref_check_population(S);P=*S.population
    if(P.ids.string_id) {
        if(rows(P.ids.string_id_value)!=S.N)_cdid_ref_fail("Original string roster length changed.",498)
        if(!all(P.ids.string_id_value:==S.original_ids.string_id_value))
            _cdid_ref_fail("The original population string identities changed.",498)
    }
    else {
        if(rows(P.ids.numeric_id)!=S.N)_cdid_ref_fail("Original numeric roster length changed.",498)
        if(!all(P.ids.numeric_id:==S.original_ids.numeric_id))
            _cdid_ref_fail("The original population numeric identities changed.",498)
    }
    if(any(P.calendar:!=S.calendar))_cdid_ref_fail("The original calendar changed.",498)
}
real matrix cdid_ref_raw_query(struct cdid_reference_inputs scalar S,
    real colvector keys,real rowvector variables)
{
    struct cdid_population scalar P
    real colvector di,units,rr,sel
    real matrix out
    real scalar j,v
    P=*S.population
    _cdid_ref_check_population(S)
    if(cols(keys)!=1 | missing(keys) | any(keys:!=floor(keys)) | any(keys:<1) |
        any(variables:<1) | any(variables:>cols(S.raw_names)))
        _cdid_ref_fail("Raw reference queries require exact date/unit keys and raw column indices.",498)
    di=floor(keys:/(S.N+1));units=keys-di:*(S.N+1)
    if(any(di:>rows(P.calendar)) | any(units:<1) | any(units:>S.N))
        _cdid_ref_fail("A reference key is outside the original roster/calendar.",498)
    (void)_cdid_ref_cost(rows(keys),cols(variables),S.max_entries,S.retained_entries)
    out=J(rows(keys),cols(variables),.)
    for(j=1;j<=cols(variables);j++) {
        v=variables[j]
        if(S.kind[v]==3) {
            sel=_cdid_ref_ix(di:>0)
            if(rows(sel))out[sel,j]=P.calendar[di[sel]]
        }
        else if(S.kind[v]==1)out[,j]=S.attributes[units,v]
        else if(S.kind[v]==2) {
            sel=_cdid_ref_ix(di:>0)
            if(rows(sel)) {
                rr=cdid_population_rows(P,units[sel],P.calendar[di[sel]])
                sel=sel[_cdid_ref_ix(rr:>0)];rr=select(rr,rr:>0)
                if(rows(sel))out[sel,j]=S.data[rr,v]
            }
        }
    }
    return(out)
}
struct cdid_reference_inputs scalar cdid_ref_source(
    pointer(struct cdid_population scalar) scalar population,
    real colvector source_rows,string rowvector raw_names,
    real matrix data,real matrix frame,real rowvector has_data,real rowvector has_frame,
    string scalar time_name,real scalar max_entries)
{
    struct cdid_reference_inputs scalar S
    struct cdid_population scalar P
    real scalar q,j,i,u,change,seen,used
    real colvector present,obs_u,diff_ix
    if(population==NULL)_cdid_ref_fail("An original population pointer is required.",498)
    P=*population;q=cols(raw_names)
    if(max_entries<1 | max_entries>=. | max_entries!=floor(max_entries) | max_entries>2^53-1)
        _cdid_ref_fail("max_entries must be an explicit positive exact integer.",498)
    if(rows(data)!=P.R | cols(data)!=q | rows(frame)!=P.N | cols(frame)!=q |
        cols(has_data)!=q | cols(has_frame)!=q | rows(source_rows)!=P.R |
        cols(source_rows)!=1 | missing(source_rows) | any(source_rows:<1) |
        any(source_rows:!=floor(source_rows)) | rows(uniqrows(source_rows))!=P.R)
        _cdid_ref_fail("Raw input dimensions must match original data rows, IDs and explicit presence maps.",498)
    if(cols(uniqrows(raw_names'))!=1 | rows(uniqrows(raw_names'))!=q |
        any((has_data:!=0):&(has_data:!=1)) | any((has_frame:!=0):&(has_frame:!=1)))
        _cdid_ref_fail("Raw names must be unique and presence maps binary.",498)
    used=_cdid_ref_cost(P.R,q+3,max_entries,0)
    used=_cdid_ref_cost(P.N,2*q+2,max_entries,used)
    used=_cdid_ref_cost(rows(P.calendar),1,max_entries,used)
    S.population=population
    S.original_ids.string_id=P.ids.string_id;S.original_ids.frame_supplied=P.ids.frame_supplied
    S.original_ids.numeric_id=P.ids.numeric_id;S.original_ids.string_id_value=P.ids.string_id_value
    S.original_ids.row_unit=J(0,1,0);S.calendar=P.calendar
    S.N=P.N;S.R=P.R;S.max_entries=max_entries
    S.raw_names=raw_names;S.data=data;S.source_rows=source_rows
    S.attributes=J(P.N,q,.);S.has_data=has_data;S.has_frame=has_frame
    S.kind=J(1,q,0);S.first_change=J(1,q,0)
    obs_u=_cdid_ref_ix(P.first_row:>0)
    for(j=1;j<=q;j++) {
        if(raw_names[j]==time_name) {S.kind[j]=3;continue;}
        if(!has_data[j] & !has_frame[j])_cdid_ref_fail("An active raw variable is absent from data and frame: "+raw_names[j],459)
        change=0
        if(has_data[j]) {
            if(P.R>0 & !missing(data[,j])) {
                if(rows(obs_u))S.attributes[obs_u,j]=data[P.first_row[obs_u],j]
                diff_ix=_cdid_ref_ix(data[,j]:!=S.attributes[P.ids.row_unit,j])
                change=(rows(diff_ix) ? diff_ix[1] : 0)
            }
            else {
                for(i=1;i<=P.R;i++) {
                    if(data[i,j]>=.)continue
                    u=P.ids.row_unit[i]
                    if(S.attributes[u,j]>=.)S.attributes[u,j]=data[i,j]
                    else if(data[i,j]!=S.attributes[u,j] & change==0)change=i
                }
            }
        }
        S.first_change[j]=change
        if(has_frame[j]) {
            if(has_data[j]) {
                if(P.R>0 & !missing(data[,j])) {
                    diff_ix=_cdid_ref_ix((frame[P.ids.row_unit,j]:>=.):|(data[,j]:!=frame[P.ids.row_unit,j]))
                    if(rows(diff_ix))
                        _cdid_ref_fail("Frame and supplied data values disagree for active variable "+raw_names[j]+" at original data row "+strofreal(source_rows[diff_ix[1]])+".",459)
                }
                else {
                    for(i=1;i<=P.R;i++)if(data[i,j]<.) {
                        u=P.ids.row_unit[i]
                        if(frame[u,j]>=. | data[i,j]!=frame[u,j])
                            _cdid_ref_fail("Frame and supplied data values disagree for active variable "+raw_names[j]+" at original data row "+strofreal(source_rows[i])+".",459)
                    }
                }
            }
            S.attributes[,j]=frame[,j];S.kind[j]=1
        }
        else S.kind[j]=(change>0 ? 2 : 1)
    }
    S.retained_entries=used
    S.roles=cdid_ref_role(3,1);S.designs=cdid_ref_design(0,1)
    S.n_slots=S.n_views=S.n_recipes=0
    return(S)
}
void _cdid_ref_add(transmorphic scalar set,real colvector keys,
    real scalar view,real scalar cap,real scalar used)
{
    real scalar i
    for(i=1;i<=rows(keys);i++)if(!asarray_contains(set,keys[i])) {
        used=_cdid_ref_cost(1,2,cap,used)
        asarray(set,keys[i],view)
    }
}
void _cdid_ref_role_set(struct cdid_ref_role scalar role,string scalar formula,
    real rowvector raw,real scalar active,real scalar intercept,
    struct cdid_reference_inputs scalar S)
{
    role.formula=invtokens(tokens(formula));role.raw_index=raw
    role.active=active;role.intercept=intercept;role.slot=0
    role.keys=role.first_view=role.fit_keys=J(0,1,0)
    role.dated=(active & cols(raw)>0 ? any(S.kind[raw]:==2) | any(S.kind[raw]:==3) : 0)
}
/* Rows of one view and cohort: x rows (observed pair units), weight rows and
   selection rows, optionally restricted to an eligibility context. */
void _cdid_ref_view_rows(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar plan,
    real scalar j,real scalar ck,string scalar route,real colvector U,real colvector elig,real scalar has_grid,
    real colvector xrows,real colvector wrows,real colvector qrows)
{
    real colvector units,rr,rt,vals
    units=_cdid_ref_ix((P.cohort_known:==1):&(P.cohort:==plan.cohorts[ck]))
    if(rows(elig) & rows(units)) units=select(units,elig[units])
    wrows=select(units,(route=="baseline" ? P.A[units] : U[units]))
    qrows=select(units,U[units])
    if(route=="imputation" | route=="augmented") {xrows=units;return;}
    xrows=J(0,1,0)
    if(rows(units)==0) return
    if(has_grid) {
        rr=P.obs_grid[units,plan.raw[j,2]]
        rt=P.obs_grid[units,plan.raw[j,3]]
    }
    else {
        rr=cdid_population_rows(P,units,J(rows(units),1,P.calendar[plan.raw[j,2]]))
        rt=cdid_population_rows(P,units,J(rows(units),1,P.calendar[plan.raw[j,3]]))
    }
    vals=_cdid_ref_ix((rr:>0):&(rt:>0))
    if(rows(vals))xrows=units[vals[_cdid_ref_ix((P.y[rr[vals]]:<.):&(P.y[rt[vals]]:<.))]]
}
/* Uses every scientific post comparison and BOTH auxiliary families of the
   plan; plan.raw is (cohort-index,s-index,t-index), controls are its exact CSR. */
void cdid_ref_require(struct cdid_reference_inputs scalar S,
    struct cdid_baseline_plan_state scalar plan,real scalar anticipation,
    string scalar route,string scalar method,real colvector U,string scalar fix_weights,
    string scalar outcome_formula,string scalar treatment_formula,string scalar selection_formula,
    string scalar outcome_raw,string scalar treatment_raw,string scalar selection_raw,
    real scalar outcome_intercept,real scalar treatment_intercept,real scalar selection_intercept,
    real scalar weight_index)
{
    struct cdid_population scalar P
    transmorphic scalar a1,a2,a3,aw,recseen,usedco
    real scalar j,k,r,d,a,base,ref,anydated,wdated,used,role,c,nv,idx,ck,nc,has_grid,use_dense,use_seenco,maxk
    real scalar elig_on,e,ectx,has_skip,xp_ex
    real colvector ks,units,rr,rt,xrows,wrows,qrows,keys,vals,v1,v2,v3,vw,view_ctx,view_e,elig,elig_index
    pointer(real colvector) colvector coh_u,coh_w,coh_q
    real matrix da,seen_co,elig_excluded
    string scalar key
    P=*S.population
    if(route!="baseline" & route!="selected" & route!="imputation" & route!="augmented")
        _cdid_ref_fail("route must be baseline, selected, imputation or augmented.",498)
    if(method!="ipw" & method!="reg" & method!="dr")_cdid_ref_fail("method must be ipw, reg or dr.",498)
    if(rows(U)!=S.N | cols(U)!=1 | missing(U) | any((U:!=0):&(U:!=1)))
        _cdid_ref_fail("The explicit original-roster selection universe must be binary.",498)
    if(fix_weights!="" & fix_weights!="first_period" & fix_weights!="base_period")
        _cdid_ref_fail("Sampling weights vary within unit; specify fixweights(first_period) or fixweights(base_period).",459)
    if(weight_index<0 | weight_index>cols(S.raw_names) | weight_index!=floor(weight_index))
        _cdid_ref_fail("Invalid raw weight-column index.",498)
    S.route=route;S.method=method;S.fix_weights=fix_weights;S.unit_universe=U
    _cdid_ref_role_set(S.roles[1],outcome_formula,_cdid_ref_names(S.raw_names,outcome_raw),1,outcome_intercept,S)
    _cdid_ref_role_set(S.roles[2],treatment_formula,(method!="reg" ? _cdid_ref_names(S.raw_names,treatment_raw) : J(1,0,0)),method!="reg",treatment_intercept,S)
    _cdid_ref_role_set(S.roles[3],selection_formula,((route=="selected" | route=="augmented") ? _cdid_ref_names(S.raw_names,selection_raw) : J(1,0,0)),(route=="selected" | route=="augmented"),selection_intercept,S)
    anydated=S.roles[1].dated | S.roles[2].dated | S.roles[3].dated
    wdated=(weight_index>0 ? S.kind[weight_index]!=1 : 0)
    if(wdated & fix_weights=="")_cdid_ref_fail("Sampling weights vary within unit; specify fixweights(first_period) or fixweights(base_period).",459)
    S.weight_kind=(wdated ? "dated" : "static")
    S.cohorts=plan.cohorts;S.used_cohorts=J(0,3,0)
    nv=rows(plan.raw);S.n_views=nv;S.view_recipe=J(nv,1,0);S.view_meta=J(nv,6,.)
    /* Restricted eligibility contexts (missingreference(exclude)): the context
       of every view of target g is part of its recipe identity. */
    elig_on=0
    if(anydated & findexternal("__cdid_mref_state")!=NULL) elig_on=cdid_mrri_active()
    S.elig_index=J(0,1,.);S.elig_excluded=J(S.N,0,.)
    view_ctx=view_e=J(nv,1,0)
    if(elig_on) {
        elig_index=cdid_mrri_index();elig_excluded=cdid_mrri_excluded(S.N)
        S.elig_index=elig_index;S.elig_excluded=elig_excluded
        if(nv) {
            view_ctx=cdid_mrri_ctx_of(plan.cohorts[plan.raw[,1]])
            for(j=1;j<=nv;j++) if(view_ctx[j]>0) view_e[j]=elig_index[view_ctx[j]]
        }
    }
    has_skip=(rows(plan.raw_mref_skip)==nv & nv>0)
    S.recipes=J(0,(elig_on ? 3 : 2),.);used=_cdid_ref_cost(nv,7,S.max_entries,S.retained_entries)
    nc=rows(plan.cohorts)
    coh_u=coh_w=coh_q=J(nc,1,NULL)
    for(ck=1;ck<=nc;ck++) {
        units=_cdid_ref_ix((P.cohort_known:==1):&(P.cohort:==plan.cohorts[ck]))
        wrows=select(units,(route=="baseline" ? P.A[units] : U[units]))
        qrows=select(units,U[units])
        coh_u[ck]=&(units:+0);coh_w[ck]=&(wrows:+0);coh_q[ck]=&(qrows:+0)
    }
    has_grid=(rows(P.obs_grid)==P.N & cols(P.obs_grid)==rows(P.calendar))
    maxk=(rows(P.calendar)+1)*(S.N+1)
    use_dense=(maxk<=4000000 & S.N>=0)
    use_seenco=(nv*nc<=1000000)
    if(use_dense) {
        v1=J((S.roles[1].dated ? maxk : S.N),1,0)
        v2=J((S.roles[2].active ? (S.roles[2].dated ? maxk : S.N) : 0),1,0)
        v3=J((S.roles[3].active ? (S.roles[3].dated ? maxk : S.N) : 0),1,0)
        vw=J((wdated ? maxk : S.N),1,0)
    }
    else {
        a1=asarray_create("real",1);a2=asarray_create("real",1);a3=asarray_create("real",1);aw=asarray_create("real",1)
    }
    recseen=asarray_create("string",1)
    if(use_seenco)seen_co=J(nv,nc,0)
    else usedco=asarray_create("real",2)
    for(j=1;j<=nv;j++) {
        ref=_cdid_pr_reference(P.calendar,plan.cohorts[plan.raw[j,1]],anticipation)
        base=cdid_ref_lookup(ref,P.calendar)
        d=(j<=plan.M ? base : plan.raw[j,2])
        if(anydated & d==0)_cdid_ref_fail("A used view has no valid design reference date.",459)
        a=(wdated ? (fix_weights=="first_period" ? 1 : base) : 0)
        if(wdated & a==0)_cdid_ref_fail("A used view has no valid weight anchor.",459)
        e=view_e[j];ectx=view_ctx[j]
        if(elig_on) {
            da=((anydated ? d : 0),a,e)
            key=strofreal(da[1],"%21.0f")+":"+strofreal(a,"%21.0f")+":"+strofreal(e,"%21.0f")
        }
        else {
            da=((anydated ? d : 0),a)
            key=strofreal(da[1],"%21.0f")+":"+strofreal(a,"%21.0f")
        }
        if(!asarray_contains(recseen,key)) {
            used=_cdid_ref_cost(1,2,S.max_entries,used)
            S.recipes=S.recipes\da;asarray(recseen,key,rows(S.recipes))
        }
        r=asarray(recseen,key);S.view_recipe[j]=r
        S.view_meta[j,]=(plan.cohorts[plan.raw[j,1]],P.calendar[plan.raw[j,2]],P.calendar[plan.raw[j,3]],d,a,j<=plan.M)
        /* An auxiliary comparison made unavailable in its eligibility context demands nothing. */
        if(has_skip) if(plan.raw_mref_skip[j]) continue
        elig=(ectx>0 ? !elig_excluded[,ectx] : J(0,1,.))
        ks=plan.raw[j,1]\plan.raw_K_index[plan.raw_K_ptr[j]..plan.raw_K_ptr[j+1]-1]
        for(k=1;k<=rows(ks);k++) {
            ck=ks[k]
            if(use_seenco) {
                if(!seen_co[r,ck]) {
                    used=_cdid_ref_cost(1,3,S.max_entries,used)
                    seen_co[r,ck]=1
                    S.used_cohorts=S.used_cohorts\(r,ck,1+2*S.roles[2].active+4*S.roles[3].active+8)
                }
            }
            else if(!asarray_contains(usedco,(r,ck))) {
                used=_cdid_ref_cost(1,3,S.max_entries,used)
                asarray(usedco,(r,ck),1)
                S.used_cohorts=S.used_cohorts\(r,ck,1+2*S.roles[2].active+4*S.roles[3].active+8)
            }
            if(ectx>0) _cdid_ref_view_rows(P,plan,j,ck,route,U,elig,has_grid,xrows=.,wrows=.,qrows=.)
            else {
            units=*(coh_u[ck]);wrows=*(coh_w[ck]);qrows=*(coh_q[ck])
            if(route=="imputation" | route=="augmented") {xrows=units;}
            else if(has_grid) {
                if(rows(units)) {
                    rr=P.obs_grid[units,plan.raw[j,2]]
                    rt=P.obs_grid[units,plan.raw[j,3]]
                    vals=_cdid_ref_ix((rr:>0):&(rt:>0))
                    if(rows(vals))xrows=units[vals[_cdid_ref_ix((P.y[rr[vals]]:<.):&(P.y[rt[vals]]:<.))]]
                    else xrows=J(0,1,0)
                }
                else xrows=J(0,1,0)
            }
            else {
                rr=cdid_population_rows(P,units,J(rows(units),1,P.calendar[plan.raw[j,2]]))
                rt=cdid_population_rows(P,units,J(rows(units),1,P.calendar[plan.raw[j,3]]))
                xrows=J(0,1,0)
                vals=_cdid_ref_ix((rr:>0):&(rt:>0))
                if(rows(vals))xrows=units[vals[_cdid_ref_ix((P.y[rr[vals]]:<.):&(P.y[rt[vals]]:<.))]]
            }
            }
            if(use_dense) {
                if(rows(xrows)) {
                    keys=(S.roles[1].dated ? d : 0)*(S.N+1):+xrows
                    vals=_cdid_ref_ix(v1[keys]:==0)
                    if(rows(vals)) {
                        used=_cdid_ref_cost(rows(vals),2,S.max_entries,used)
                        v1[keys[vals]]=J(rows(vals),1,j)
                    }
                    if(S.roles[2].active) {
                        keys=(S.roles[2].dated ? d : 0)*(S.N+1):+xrows
                        vals=_cdid_ref_ix(v2[keys]:==0)
                        if(rows(vals)) {
                            used=_cdid_ref_cost(rows(vals),2,S.max_entries,used)
                            v2[keys[vals]]=J(rows(vals),1,j)
                        }
                    }
                }
                if(S.roles[3].active & rows(qrows)) {
                    keys=(S.roles[3].dated ? d : 0)*(S.N+1):+qrows
                    vals=_cdid_ref_ix(v3[keys]:==0)
                    if(rows(vals)) {
                        used=_cdid_ref_cost(rows(vals),2,S.max_entries,used)
                        v3[keys[vals]]=J(rows(vals),1,j)
                    }
                }
                if(rows(wrows)) {
                    keys=a*(S.N+1):+wrows
                    vals=_cdid_ref_ix(vw[keys]:==0)
                    if(rows(vals)) {
                        used=_cdid_ref_cost(rows(vals),2,S.max_entries,used)
                        vw[keys[vals]]=J(rows(vals),1,j)
                    }
                }
            }
            else {
                _cdid_ref_add(a1,(S.roles[1].dated ? d : 0):*(S.N+1):+xrows,j,S.max_entries,used)
                if(S.roles[2].active)_cdid_ref_add(a2,(S.roles[2].dated ? d : 0):*(S.N+1):+xrows,j,S.max_entries,used)
                if(S.roles[3].active)_cdid_ref_add(a3,(S.roles[3].dated ? d : 0):*(S.N+1):+qrows,j,S.max_entries,used)
                _cdid_ref_add(aw,a*(S.N+1):+wrows,j,S.max_entries,used)
            }
        }
    }
    if(use_dense) {
        S.roles[1].keys=_cdid_ref_ix(v1:>0)
        S.roles[1].first_view=(rows(S.roles[1].keys) ? v1[S.roles[1].keys] : J(0,1,0))
        S.roles[2].keys=_cdid_ref_ix(v2:>0)
        S.roles[2].first_view=(rows(S.roles[2].keys) ? v2[S.roles[2].keys] : J(0,1,0))
        S.roles[3].keys=_cdid_ref_ix(v3:>0)
        S.roles[3].first_view=(rows(S.roles[3].keys) ? v3[S.roles[3].keys] : J(0,1,0))
        S.weight_keys=_cdid_ref_ix(vw:>0)
    }
    else {
        S.roles[1].keys=sort(vec(asarray_keys(a1)),1)
        S.roles[2].keys=sort(vec(asarray_keys(a2)),1)
        S.roles[3].keys=sort(vec(asarray_keys(a3)),1)
        for(role=1;role<=3;role++) {
            keys=S.roles[role].keys;S.roles[role].first_view=J(rows(keys),1,0)
            for(j=1;j<=rows(keys);j++) {
                if(role==1)S.roles[role].first_view[j]=asarray(a1,keys[j])
                else if(role==2)S.roles[role].first_view[j]=asarray(a2,keys[j])
                else S.roles[role].first_view[j]=asarray(a3,keys[j])
            }
        }
        S.weight_keys=sort(vec(asarray_keys(aw)),1)
    }
    for(role=1;role<=3;role++)if(S.roles[role].active) {
        if(S.roles[role].dated)S.roles[role].fit_keys=S.roles[role].keys
        else {
            if(role==3)S.roles[role].fit_keys=_cdid_ref_ix(U)
            else if(route=="imputation" | route=="augmented" | (route=="baseline" & method=="ipw"))S.roles[role].fit_keys=(1::S.N)
            else S.roles[role].fit_keys=P.A_to_F
        }
    }
    /* Every required raw reference covariate value must be supplied. */
    if(anydated) {
        /* Cross-pair comparisons read covariates on their own rows. */
        xp_ex=0
        if(findexternal("__cdid_mref_state")!=NULL) if(_cdid_ref_xp_on()) if(cdid_mrri_policy()=="exclude") xp_ex=1
        if(xp_ex) _cdid_ref_xp_check(S,P,plan,anticipation,route,U,has_grid,elig_on,view_ctx)
        else _cdid_ref_check_missing(S,P,plan,anticipation,route,U,has_grid,view_ctx,view_e)
    }
    if(weight_index>0) {
        S.weights=cdid_ref_raw_query(S,S.weight_keys,weight_index)
        S.V=cdid_ref_raw_query(S,(S.kind[weight_index]!=1 ? S.N+1 : 0):+(1::S.N),weight_index)
    }
    else {S.weights=J(rows(S.weight_keys),1,1);S.V=J(S.N,1,1);}
    if(missing(S.weights) | any(S.weights:<=0))_cdid_ref_weight_fail(S,P,weight_index)
    used=_cdid_ref_cost(rows(S.weight_keys)+S.N,1,S.max_entries,used)
    S.retained_entries=used;S.n_recipes=rows(S.recipes)
}
/* A required anchor weight is missing or nonpositive. Weights are never carried
   forward, substituted or excluded; missingreference() concerns covariates only. */
void _cdid_ref_weight_fail(struct cdid_reference_inputs scalar S,struct cdid_population scalar P,real scalar weight_index)
{
    real colvector bad
    real scalar d,unit,key
    string scalar item
    bad=_cdid_ref_ix((S.weights:>=.) :| (S.weights:<=0));item=""
    if(rows(bad) & weight_index>0) {
        key=S.weight_keys[bad[1]];d=floor(key/(S.N+1));unit=key-d*(S.N+1)
        if(d>=1 & d<=length(P.calendar) & unit>=1 & unit<=S.N) {
            item=" (id "+(P.ids.string_id ? P.ids.string_id_value[unit] : _cdid_pr_num(P.ids.numeric_id[unit]))+
                ", date "+_cdid_pr_num(P.calendar[d])+", variable "+S.raw_names[weight_index]+
                (rows(bad)>1 ? "; and "+strofreal(rows(bad)-1)+" more" : "")+")"
        }
    }
    _cdid_ref_fail("fixweights() needs a positive sampling weight at its anchor date for every unit it uses; a weight is missing or nonpositive there"+item+
        ". Supply the weight on that unit-date row (the outcome may be missing on that row); weights are never carried forward or substituted, and missingreference() applies to covariates only.",459)
}
/* Cross-pair comparisons (crosspair) read the covariates of target cohort g at its
   reference date b_g on every repeat observer of g (every unit of g in the selection
   universe with selection()) and on every unit of a control cohort observed at two
   dates before that cohort is treated; the selection design is read on the selection
   universe of g and of its controls. Under missingreference(exclude) a missing value on
   these rows starts the eligibility context of b_g, and the preparation is repeated. */
void _cdid_ref_xp_check(struct cdid_reference_inputs scalar S,struct cdid_population scalar P,
    struct cdid_baseline_plan_state scalar plan,real scalar anticipation,string scalar route,
    real colvector U,real scalar has_grid,real scalar elig_on,real colvector view_ctx)
{
    real scalar i,j,k,t,g,b,gi,ck,T,tmax,cx
    real colvector tg,Ks,units,cu,cf,dem,qd,affected,vx,vs,rr,cnt,el,aux,vv
    string scalar cg
    real matrix raw
    if(plan.M<1) return
    /* Values demanded only by same-pair pretrend comparisons at their own earlier dates are
       not read by the cross-pair comparisons: when such a value is missing, those
       requirement rows are dropped. */
    for(k=1;k<=3;k++) if(S.roles[k].active & S.roles[k].dated & rows(S.roles[k].keys)) {
        aux=(S.roles[k].first_view:>plan.M)
        if(!any(aux)) continue
        vv=select(S.roles[k].raw_index',S.kind[S.roles[k].raw_index]':!=3)
        if(rows(vv)==0) continue
        raw=cdid_ref_raw_query(S,select(S.roles[k].keys,aux),vv')
        if(!missing(raw)) continue
        S.roles[k].keys=select(S.roles[k].keys,!aux)
        S.roles[k].first_view=select(S.roles[k].first_view,!aux)
        S.roles[k].fit_keys=S.roles[k].keys
    }
    T=rows(P.calendar);affected=J(0,1,.)
    cg=_cdid_ref_xp_cg()
    vx=J(0,1,.);vs=J(0,1,.)
    for(k=1;k<=2;k++) if(S.roles[k].active & S.roles[k].dated)
        vx=vx\select(S.roles[k].raw_index',S.kind[S.roles[k].raw_index]':!=3)
    if(rows(vx)) vx=uniqrows(vx)
    if(S.roles[3].active & S.roles[3].dated) vs=select(S.roles[3].raw_index',S.kind[S.roles[3].raw_index]':!=3)
    tg=uniqrows(plan.raw[|1,1\plan.M,1|])
    for(i=1;i<=rows(tg);i++) {
        gi=tg[i];g=plan.cohorts[gi]
        b=cdid_ref_lookup(_cdid_pr_reference(P.calendar,g,anticipation),P.calendar)
        if(b==0) continue
        cx=0
        for(j=1;j<=plan.M;j++) if(plan.raw[j,1]==gi) cx=view_ctx[j]
        /* Clean control cohorts of the cross-pair comparisons. */
        if(cg=="notyettreated") Ks=_cdid_ref_ix(plan.cohorts:>g)
        else Ks=_cdid_ref_ix(plan.cohorts:>=.)
        units=_cdid_ref_ix((P.cohort_known:==1):&(P.cohort:==g))
        /* Imputed comparisons read every unit of the target and of its controls. */
        if(route=="imputation" | route=="augmented") dem=units
        else dem=select(units,(route=="selected" ? U[units] : P.A[units]))
        qd=select(units,U[units])
        for(k=1;k<=rows(Ks);k++) {
            ck=Ks[k]
            cf=_cdid_ref_ix((P.cohort_known:==1):&(P.cohort:==plan.cohorts[ck]))
            cu=select(cf,P.A[cf])
            if(rows(cu)==0) continue
            tmax=T
            if(plan.cohorts[ck]<.) tmax=sum((P.calendar:+anticipation):<plan.cohorts[ck])
            if(tmax<2) continue
            cnt=J(rows(cu),1,0)
            for(t=1;t<=tmax;t++) {
                if(has_grid) rr=P.obs_grid[cu,t]
                else rr=cdid_population_rows(P,cu,J(rows(cu),1,P.calendar[t]))
                rr=editmissing(rr,0)
                cnt=cnt+((rr:>0):&(P.y[rr:+(rr:<=0)]:<.))
            }
            if(!any(cnt:>=2)) continue
            qd=qd\select(cf,U[cf])
            if(route=="imputation" | route=="augmented") dem=dem\cf
            else dem=dem\select(cu,cnt:>=2)
        }
        if(elig_on & cx>0) {
            el=!S.elig_excluded[,cx]
            if(rows(dem)) dem=select(dem,el[dem])
            if(rows(qd)) qd=select(qd,el[qd])
        }
        if(rows(vx) & rows(dem)) {
            raw=cdid_ref_raw_query(S,b*(S.N+1):+uniqrows(dem),vx')
            if(missing(raw)) {affected=affected\b;continue;}
        }
        if(rows(vs) & rows(qd)) {
            raw=cdid_ref_raw_query(S,b*(S.N+1):+uniqrows(qd),vs')
            if(missing(raw)) affected=affected\b
        }
    }
    if(rows(affected)==0) return
    if(elig_on)
        _cdid_ref_fail("A cross-pair comparison demands an unknown reference value inside its eligibility context; the context and the requirement scan disagree.",498)
    _cdid_mrri_build(S,uniqrows(affected),_cdid_mrri_vars(S))
    _cdid_mrri_set_pending()
    exit(459)
}
string scalar _cdid_ref_route_text(string scalar route)
{
    if(route=="selected") return("with selection() and selmethod(ipw)")
    if(route=="imputation") return("with selmethod(imputation)")
    if(route=="augmented") return("with selmethod(dr)")
    return("for this specification")
}
/* Label of a requirement key (as in the requirement scan of the first demanding view). */
string scalar _cdid_ref_key_label(struct cdid_population scalar P,struct cdid_baseline_plan_state scalar plan,
    real scalar j,real scalar unit)
{
    string scalar kc
    kc=(P.cohort_known[unit]==1 ? _cdid_pr_num(P.cohort[unit]) : "unknown")
    return("support ("+_cdid_pr_num(P.calendar[plan.raw[j,2]])+"->"+_cdid_pr_num(P.calendar[plan.raw[j,3]])+
        ", cohort "+kc+") of "+(j<=plan.M ? "post link "+strofreal(j) : "auxiliary report "+strofreal(j-plan.M))+
        " (target "+_cdid_pr_num(plan.cohorts[plan.raw[j,1]])+")")
}
/* Raw required values of the dated formulas. Under missingreference(exclude)
   a missing value demanded by a main comparison starts cohort-specific
   eligibility contexts at the affected reference dates; one demanded only by an
   auxiliary pretrend comparison makes that comparison unavailable. Both are
   signalled to the caller, which repeats the preparation. */
void _cdid_ref_check_missing(struct cdid_reference_inputs scalar S,struct cdid_population scalar P,
    struct cdid_baseline_plan_state scalar plan,real scalar anticipation,string scalar route,
    real colvector U,real scalar has_grid,real colvector view_ctx,real colvector view_e)
{
    real scalar role,i,n,j,k,ck,d,nshow,unit,base,anymain,has_skip
    real colvector vars,keys,miss,mk1,mk3,bad,ks,xrows,wrows,qrows,elig,hit,affected
    real matrix raw,allhits
    string scalar msg,policy,item
    mk1=mk3=J(0,1,.);msg="";n=0;nshow=0;allhits=J(0,3,.)
    for(role=1;role<=3;role++) {
        if(!S.roles[role].active | !S.roles[role].dated) continue
        vars=select(S.roles[role].raw_index',S.kind[S.roles[role].raw_index]':!=3)
        if(rows(vars)==0 | rows(S.roles[role].keys)==0) continue
        raw=cdid_ref_raw_query(S,S.roles[role].keys,vars')
        miss=_cdid_ref_ix(rowmissing(raw):>0)
        if(rows(miss)==0) continue
        keys=S.roles[role].keys[miss]
        if(role==3) mk3=mk3\keys
        else mk1=mk1\keys
        for(i=1;i<=rows(miss);i++) for(k=1;k<=rows(vars);k++) if(raw[miss[i],k]>=.)
            allhits=allhits\(S.roles[role].keys[miss[i]],vars[k],S.roles[role].first_view[miss[i]])
    }
    if(rows(mk1)+rows(mk3)==0) return
    /* One entry per (unit-date row, variable), in row order. */
    allhits=sort(allhits,(1,2,3))
    for(i=1;i<=rows(allhits);i++) {
        if(i>1) if(allhits[i,1]==allhits[i-1,1] & allhits[i,2]==allhits[i-1,2]) continue
        n++
        if(nshow>=5) continue
        d=floor(allhits[i,1]/(S.N+1));unit=allhits[i,1]-d*(S.N+1)
        item="id "+(P.ids.string_id ? P.ids.string_id_value[unit] : _cdid_pr_num(P.ids.numeric_id[unit]))+", date "+(d>0 ? _cdid_pr_num(P.calendar[d]) : "unit attribute")+
            ", variable "+S.raw_names[allhits[i,2]]+", "+_cdid_ref_key_label(P,plan,allhits[i,3],unit)
        msg=msg+(nshow ? "; " : "")+item;nshow++
    }
    if(n>nshow) msg=msg+"; and "+strofreal(n-nshow)+" more"
    policy=(findexternal("__cdid_mref_state")!=NULL ? cdid_mrri_policy() : "")
    if(policy!="exclude")
        _cdid_ref_fail("Missing required reference covariate values (time-varying covariates are never carried forward, substituted or deleted): "+msg+
            ". Supply each missing covariate on its unit-date row; the outcome may be missing on that row.",459)
    if(_cdid_ref_xp_on()) cdid_mrri_require_route("with crosspair")
    /* Views demanding a missing key. */
    mk1=uniqrows(sort(mk1,1));mk3=(rows(mk3) ? uniqrows(sort(mk3,1)) : J(0,1,.))
    has_skip=(rows(plan.raw_mref_skip)==rows(plan.raw) & rows(plan.raw)>0)
    bad=J(rows(plan.raw),1,0)
    for(j=1;j<=rows(plan.raw);j++) {
        if(has_skip) if(plan.raw_mref_skip[j]) continue
        base=cdid_ref_lookup(_cdid_pr_reference(P.calendar,plan.cohorts[plan.raw[j,1]],anticipation),P.calendar)
        d=(j<=plan.M ? base : plan.raw[j,2])
        elig=(view_ctx[j]>0 ? !S.elig_excluded[,view_ctx[j]] : J(0,1,.))
        ks=plan.raw[j,1]\plan.raw_K_index[plan.raw_K_ptr[j]..plan.raw_K_ptr[j+1]-1]
        for(k=1;k<=rows(ks);k++) {
            ck=ks[k]
            _cdid_ref_view_rows(P,plan,j,ck,route,U,elig,has_grid,xrows=.,wrows=.,qrows=.)
            if(rows(mk1) & rows(xrows) & (S.roles[1].dated | S.roles[2].dated)) {
                hit=cdid_ref_lookup(d*(S.N+1):+xrows,mk1)
                if(any(hit:>0)) {bad[j]=1;break;}
            }
            if(rows(mk3) & rows(qrows) & S.roles[3].dated) {
                hit=cdid_ref_lookup(d*(S.N+1):+qrows,mk3)
                if(any(hit:>0)) {bad[j]=1;break;}
            }
        }
    }
    anymain=(plan.M>0 ? any(bad[|1\plan.M|]) : 0)
    if(anymain) {
        if(cdid_mrri_active())
            _cdid_ref_fail("A main comparison demands an unknown reference value inside its eligibility context; the context and the requirement scan disagree.",498)
        affected=J(0,1,.)
        for(j=1;j<=plan.M;j++) if(bad[j])
            affected=affected\cdid_ref_lookup(_cdid_pr_reference(P.calendar,plan.cohorts[plan.raw[j,1]],anticipation),P.calendar)
        _cdid_mrri_build(S,uniqrows(affected),_cdid_mrri_vars(S))
        _cdid_mrri_set_pending()
        exit(459)
    }
    if(!any(bad))
        _cdid_ref_fail("Missing required reference covariate values (time-varying covariates are never carried forward, substituted or deleted): "+msg+".",459)
    for(j=plan.M+1;j<=rows(plan.raw);j++) if(bad[j])
        _cdid_mrri_aux_missing_add(j-plan.M,"Unavailable: this auxiliary comparison reads covariates at its own earlier date "+
            _cdid_pr_num(P.calendar[plan.raw[j,2]])+", where a required raw value of an eligible unit is missing "+
            "(missingreference(exclude); eligibility contexts and main comparisons are unchanged)")
    _cdid_mrri_set_pending()
    exit(459)
}

real scalar _cdid_ref_is_plain(string scalar formula,string rowvector raw_names)
{
    string rowvector tok
    tok=tokens(formula)
    if(cols(tok)!=cols(raw_names))return(0)
    if(cols(tok)==0)return(1)
    return(all(tok:==raw_names))
}
void _cdid_ref_fast_encode(struct cdid_reference_inputs scalar S,real scalar slot)
{
    real matrix raw,X
    string rowvector terms,cn
    real scalar n,w,ic
    if(slot<1 | slot>S.n_slots | slot!=floor(slot) | S.designs[slot].encoded)
        _cdid_ref_fail("A design slot must be encoded exactly once.",498)
    terms=S.designs[slot].raw_names
    n=rows(S.designs[slot].keys);w=cols(terms);ic=S.designs[slot].intercept
    (void)_cdid_ref_cost(n,2*w+4,S.max_entries,S.retained_entries)
    raw=cdid_ref_raw_query(S,S.designs[slot].keys,S.designs[slot].raw_index)
    if(missing(raw))_cdid_ref_fail("A time-varying covariate is missing at a reference date the design requires (a target cohort's last pre-treatment date or a pretrend comparison endpoint); supply the covariate at those dates, including on rows whose outcome is missing.",459)
    cdid_ref_encoding_budget(S,slot,w)
    if(ic) {X=(J(n,1,1),raw);cn=("_cons",terms);}
    else {X=raw;cn=terms;}
    S.retained_entries=_cdid_ref_cost(rows(X),cols(X),S.max_entries,S.retained_entries)
    S.designs[slot].value=X;S.designs[slot].columns=cn
    S.designs[slot].expanded=terms;S.designs[slot].omitted=J(1,w,0)
    S.designs[slot].column_map=(w ? (ic+1..ic+w) : J(1,0,0))
    S.designs[slot].encoded=1
}

/* Lexical raw-variable collection only. Actual FV parsing, expansion and substitution
   are performed by native fvexpand/fvrevar on the canonical rows, once per slot. */
string rowvector cdid_ref_formula_names(string scalar formula)
{
    string rowvector terms,out
    string scalar t,nm,op
    real scalar j,pos
    terms=tokens(subinstr(formula,"#"," "));out=J(1,0,"")
    for(j=1;j<=cols(terms);j++) {
        t=terms[j];pos=strpos(t,".");nm=t
        if(pos>0) {
            op=substr(t,1,pos-1);nm=substr(t,pos+1,.)
            if(op!="c" & op!="i" & !regexm(op,"^ib[n0-9]+$") & !regexm(op,"^[0-9]+[bno]*$"))
                _cdid_ref_fail("Use canonical c./i./ib#/ibn FV tokens; supply generated covariates for time-series/transforms.",198)
        }
        if(!regexm(nm,"^[A-Za-z_][A-Za-z_0-9]*$"))
            _cdid_ref_fail("Reference formulas require full canonical raw variable names.",198)
        if(!any(out:==nm))out=out,nm
    }
    return(out)
}
void cdid_ref_designs(struct cdid_reference_inputs scalar S)
{
    struct cdid_ref_design scalar D
    real scalar role,k,slot,j
    real colvector keys
    string scalar key,other
    S.designs=cdid_ref_design(0,1);S.n_slots=0
    for(role=1;role<=3;role++) {
        if(!S.roles[role].active)continue
        slot=0
        for(k=1;k<=S.n_slots;k++) {
            if(S.designs[k].formula==S.roles[role].formula &
                S.designs[k].intercept==S.roles[role].intercept) {
                if(S.roles[role].dated)slot=k
                else if(rows(S.designs[k].keys)==rows(S.roles[role].fit_keys)) {
                    if(all(S.designs[k].keys:==S.roles[role].fit_keys))slot=k
                }
            }
        }
        if(slot==0) {
            D.formula=S.roles[role].formula;D.intercept=S.roles[role].intercept
            D.dated=S.roles[role].dated;D.encoded=0
            D.raw_index=S.roles[role].raw_index;D.raw_names=S.raw_names[D.raw_index]
            D.keys=J(0,1,0);D.value=J(0,0,.)
            D.columns=D.expanded=J(1,0,"");D.omitted=D.column_map=J(1,0,0)
            S.designs=S.designs\D;S.n_slots=S.n_slots+1;slot=S.n_slots
        }
        S.roles[role].slot=slot
        S.designs[slot].keys=uniqrows(sort(S.designs[slot].keys\S.roles[role].fit_keys,1))
    }
    for(k=1;k<=S.n_slots;k++) {
        if(rows(S.designs[k].keys)==0)
            _cdid_ref_fail("An active formula has no canonical fitting rows in the scientific union.",459)
    }
}
/* Populate a caller-owned EMPTY temporary frame. Original names make the fixed
   expanded tokens usable without a second independent parser. */
void cdid_ref_stage(struct cdid_reference_inputs scalar S,real scalar slot)
{
    real matrix raw
    real scalar j
    if(slot<1 | slot>S.n_slots | slot!=floor(slot))_cdid_ref_fail("Unknown design slot.",498)
    (void)_cdid_ref_cost(rows(S.designs[slot].keys),2*cols(S.designs[slot].raw_index)+4,S.max_entries,S.retained_entries)
    raw=cdid_ref_raw_query(S,S.designs[slot].keys,S.designs[slot].raw_index)
    if(missing(raw))_cdid_ref_fail("A time-varying covariate is missing at a reference date the design requires (a target cohort's last pre-treatment date or a pretrend comparison endpoint); supply the covariate at those dates, including on rows whose outcome is missing.",459)
    if(st_nobs()!=0 | st_nvar()!=0)_cdid_ref_fail("FV staging requires an empty owned frame.",498)
    st_addobs(rows(raw))
    for(j=1;j<=cols(raw);j++) {
        (void)st_addvar("double",S.designs[slot].raw_names[j])
        st_store(.,S.designs[slot].raw_names[j],raw[,j])
    }
}
void cdid_ref_encoded(struct cdid_reference_inputs scalar S,real scalar slot,
    string scalar generated,string scalar expanded,real rowvector omitted)
{
    real matrix X
    real rowvector map
    real scalar j,k
    string rowvector terms,cn
    terms=tokens(expanded)
    if(slot<1 | slot>S.n_slots | S.designs[slot].encoded)
        _cdid_ref_fail("A design slot must be encoded exactly once.",498)
    if(cols(terms)!=cols(omitted))_cdid_ref_fail("FV omission metadata is not aligned.",498)
    if(strtrim(generated)=="")X=J(rows(S.designs[slot].keys),0,.)
    else X=st_data(.,tokens(generated))
    if(rows(X)!=rows(S.designs[slot].keys) | cols(X)!=sum(omitted:==0) | missing(X))
        _cdid_ref_fail("Native FV encoding is nonfinite or lost canonical row/column alignment.",459)
    map=J(1,cols(terms),0);k=S.designs[slot].intercept
    for(j=1;j<=cols(terms);j++)if(!omitted[j]) {k++;map[j]=k;}
    cn=select(terms,omitted:==0)
    if(S.designs[slot].intercept) {X=(J(rows(X),1,1),X);cn=("_cons",cn);}
    S.retained_entries=_cdid_ref_cost(rows(X),cols(X),S.max_entries,S.retained_entries)
    S.designs[slot].value=X;S.designs[slot].columns=cn
    S.designs[slot].expanded=terms;S.designs[slot].omitted=omitted
    S.designs[slot].column_map=map;S.designs[slot].encoded=1
}
real matrix cdid_ref_design_rows(struct cdid_reference_inputs scalar S,
    real scalar role,real colvector keys)
{
    real colvector hit
    real scalar slot
    _cdid_ref_check_population(S)
    if(role<1 | role>3 | !S.roles[role].active)_cdid_ref_fail("The requested design role is inactive.",498)
    slot=S.roles[role].slot
    if(!S.designs[slot].encoded)_cdid_ref_fail("The canonical design is not yet encoded.",498)
    hit=cdid_ref_lookup(keys,S.designs[slot].keys)
    if(any(hit:==0))_cdid_ref_fail("A design query is outside its fixed required union; grid prediction is a separate consumer.",459)
    return(S.designs[slot].value[hit,.])
}
real colvector cdid_ref_weight_rows(struct cdid_reference_inputs scalar S,
    real scalar recipe,real colvector units)
{
    real colvector hit
    _cdid_ref_check_population(S)
    if(recipe<1 | recipe>S.n_recipes | recipe!=floor(recipe) | cols(units)!=1 |
        missing(units) | any(units:<1) | any(units:>S.N) | any(units:!=floor(units)))
        _cdid_ref_fail("Weight queries need an actual design block and original roster positions.",498)
    hit=cdid_ref_lookup(S.recipes[recipe,2]*(S.N+1):+units,S.weight_keys)
    if(any(hit:==0))_cdid_ref_fail("A weight query is outside its actual required support.",459)
    return(vec(S.weights[hit]))
}
/* Ado boundary: take raw history in existing population row order, keep explicit
   field presence, and align an optional frame without altering caller rows. */
void _cdid_ref_capture(string scalar pop_name,string scalar rows_name,string scalar store,
    string scalar ivar,string scalar time,string scalar raw_names,real scalar cap)
{
    pointer() scalar pp,pr,ps
    struct cdid_population scalar P
    struct cdid_reference_inputs scalar S
    real colvector rr
    real rowvector v_idx,j_idx
    string rowvector names
    real scalar j,v,all_obs
    pp=findexternal(pop_name)
    if(pp==NULL | findexternal(store)!=NULL)_cdid_ref_fail("Population is absent or output name already exists.",498)
    P=*pp;names=tokens(raw_names)
    if(rows_name=="")rr=(1::P.R)
    else {pr=findexternal(rows_name);if(pr==NULL)_cdid_ref_fail("Explicit source-row map is absent.",498);rr=*pr;}
    if(rows(rr)!=P.R | cols(rr)!=1 | missing(rr) | any(rr:<1) | any(rr:>st_nobs()) | any(rr:!=floor(rr)))
        _cdid_ref_fail("Source row map must retain exactly the population's actual supplied history rows.",498)
    all_obs=(P.R==st_nobs() & P.R>0 ? all(rr:==(1::P.R)) : 0)
    if(any((all_obs ? st_data(.,time) : st_data(rr,time)):!=P.date))_cdid_ref_fail("Source dates no longer match the original population history.",498)
    if(P.ids.string_id) {
        if(any((all_obs ? st_sdata(.,ivar) : st_sdata(rr,ivar)):!=vec(P.ids.string_id_value[P.ids.row_unit])))_cdid_ref_fail("Source string IDs changed.",498)
    }
    else if(any((all_obs ? st_data(.,ivar) : st_data(rr,ivar)):!=vec(P.ids.numeric_id[P.ids.row_unit])))_cdid_ref_fail("Source numeric IDs changed.",498)
    (void)_cdid_ref_cost(P.R,cols(names)+2,cap,0)
    (void)_cdid_ref_cost(P.N,2*cols(names),cap,P.R*(cols(names)+2))
    S.population=pp;S.N=P.N;S.R=P.R;S.source_rows=rr;S.raw_names=names
    S.data=J(P.R,cols(names),.);S.attributes=J(P.N,cols(names),.)
    S.has_data=S.has_frame=J(1,cols(names),0)
    v_idx=J(1,0,0);j_idx=J(1,0,0)
    for(j=1;j<=cols(names);j++) {
        v=_st_varindex(names[j])
        if(v<. & v>0) {v_idx=(v_idx,v);j_idx=(j_idx,j);}
    }
    if(cols(v_idx)>0) {
        S.data[,j_idx]=(all_obs ? st_data(.,v_idx) : st_data(rr,v_idx))
        S.has_data[j_idx]=J(1,cols(j_idx),1)
    }
    ps=crexternal(store);if(ps==NULL)_cdid_ref_fail("Cannot create the owned raw staging object.",498)
    *ps=S
}
void _cdid_ref_frame(string scalar store,string scalar ivar,string scalar map_name)
{
    pointer() scalar ps,pm
    struct cdid_reference_inputs scalar S
    struct cdid_population scalar P
    struct cdid_population_ids scalar ids
    real colvector align,fid
    real rowvector v_idx,j_idx
    string matrix map
    string scalar nm
    real scalar j,v,i,id_order
    ps=findexternal(store);if(ps==NULL)_cdid_ref_fail("Raw staging object is absent.",498)
    S=*ps;P=*S.population
    if(P.ids.string_id) {
        ids=cdid_population_ids_string(P.ids.string_id_value,st_sdata(.,ivar),1)
        align=ids.row_unit
        if(rows(align)!=S.N | rows(ids.numeric_id)+rows(ids.string_id_value)!=S.N)
            _cdid_ref_fail("The supplied raw frame must have exactly the original roster IDs.",498)
    }
    else {
        fid=st_data(.,ivar)
        if(rows(fid)==S.N & rows(P.ids.numeric_id)==S.N & !missing(fid) & all(fid:==P.ids.numeric_id))
            align=(1::S.N)
        else {
            ids=cdid_population_ids_numeric(P.ids.numeric_id,fid,1)
            align=ids.row_unit
            if(rows(align)!=S.N | rows(ids.numeric_id)+rows(ids.string_id_value)!=S.N)
                _cdid_ref_fail("The supplied raw frame must have exactly the original roster IDs.",498)
        }
    }
    map=J(0,2,"")
    if(map_name!="") {pm=findexternal(map_name);if(pm==NULL)_cdid_ref_fail("Frame raw mapping is absent.",498);map=*pm;}
    if(cols(map)!=2)_cdid_ref_fail("Frame raw mapping must have two string columns.",498)
    v_idx=J(1,0,0);j_idx=J(1,0,0)
    for(j=1;j<=cols(S.raw_names);j++) {
        nm=S.raw_names[j]
        for(i=1;i<=rows(map);i++)if(map[i,1]==S.raw_names[j])nm=map[i,2]
        v=_st_varindex(nm)
        if(v<. & v>0) {v_idx=(v_idx,v);j_idx=(j_idx,j);}
    }
    if(cols(v_idx)>0) {
        id_order=(S.N==st_nobs() & S.N>0 ? all(align:==(1::S.N)) : 0)
        S.attributes[,j_idx]=(id_order ? st_data(.,v_idx) : st_data(align,v_idx))
        S.has_frame[j_idx]=J(1,cols(j_idx),1)
    }
    *ps=S
}
void _cdid_ref_finalize(string scalar stage,string scalar plan_name,string scalar time,
    real scalar anticipation,string scalar route,string scalar method,string scalar universe,
    string scalar fixweights,string scalar xf,string scalar tf,string scalar sf,
    string scalar xr,string scalar tr,string scalar sr,real scalar xi,real scalar ti,real scalar si,
    string scalar weight,real scalar cap)
{
    pointer() scalar ps,pp
    struct cdid_reference_inputs scalar S,out
    struct cdid_population scalar P
    struct cdid_baseline_plan_state scalar plan
    real colvector U
    real rowvector wi
    ps=findexternal(stage);pp=findexternal(plan_name)
    if(ps==NULL | pp==NULL)_cdid_ref_fail("Raw stage or complete scientific plan is absent.",498)
    S=*ps;P=*S.population;plan=*pp
    out=cdid_ref_source(S.population,S.source_rows,S.raw_names,S.data,S.attributes,S.has_data,S.has_frame,time,cap)
    if(universe=="F")U=J(P.N,1,1)
    else if(universe=="B")U=P.B
    else if(universe=="A")U=P.A
    else _cdid_ref_fail("Selection universe must be F, B or A.",498)
    wi=(weight=="" ? 0 : _cdid_ref_names(out.raw_names,weight))
    cdid_ref_require(out,plan,anticipation,route,method,U,fixweights,xf,tf,sf,xr,tr,sr,xi,ti,si,wi)
    cdid_ref_designs(out);*ps=out
}

void cdid_ref_encoding_budget(struct cdid_reference_inputs scalar S,real scalar slot,real scalar width)
{
    if(slot<1 | slot>S.n_slots | width<0 | width>=. | width!=floor(width))
        _cdid_ref_fail("Invalid FV encoding width.",498)
    (void)_cdid_ref_cost(rows(S.designs[slot].keys),3*width+2*cols(S.designs[slot].raw_index)+8,
        S.max_entries,S.retained_entries)
}
void cdid_ref_publish(string scalar stage,string scalar store)
{
    pointer() scalar p,out
    struct cdid_reference_inputs scalar S
    real scalar j
    p=findexternal(stage)
    if(p==NULL | findexternal(store)!=NULL)_cdid_ref_fail("Staging is absent or destination already exists.",498)
    S=*p;cdid_ref_activate(S)
    for(j=1;j<=S.n_slots;j++)if(!S.designs[j].encoded)_cdid_ref_fail("Every active formula must be encoded before publication.",498)
    out=crexternal(store);if(out==NULL)_cdid_ref_fail("Cannot allocate the explicit reference output.",498)
    *out=S
    st_local("ref_N",strofreal(S.N));st_local("ref_recipes",strofreal(S.n_recipes))
    st_local("ref_views",strofreal(S.n_views));st_local("ref_slots",strofreal(S.n_slots))
}

void _cdid_ref_meta(string scalar name,real scalar slot)
{
    pointer() scalar p
    struct cdid_reference_inputs scalar S
    p=findexternal(name);if(p==NULL)_cdid_ref_fail("Reference staging object is absent.",498)
    S=*p;st_local("nslots",strofreal(S.n_slots))
    if(slot>0) {
        if(slot>S.n_slots)_cdid_ref_fail("Unknown reference design slot.",498)
        if(!S.designs[slot].encoded & _cdid_ref_is_plain(S.designs[slot].formula,S.designs[slot].raw_names)) {
            _cdid_ref_fast_encode(S,slot)
            *p=S
        }
        st_local("slot_encoded",strofreal(S.designs[slot].encoded))
        st_local("formula",S.designs[slot].formula)
    }
}
void _cdid_ref_stage_store(string scalar name,real scalar slot)
{
    pointer() scalar p
    struct cdid_reference_inputs scalar S
    p=findexternal(name);if(p==NULL)_cdid_ref_fail("Reference staging object is absent.",498)
    S=*p;cdid_ref_stage(S,slot)
}
void _cdid_ref_encode_store(string scalar name,real scalar slot,string scalar generated,
    string scalar expanded,real rowvector omitted)
{
    pointer() scalar p
    struct cdid_reference_inputs scalar S
    p=findexternal(name);if(p==NULL)_cdid_ref_fail("Reference staging object is absent.",498)
    S=*p;cdid_ref_encoded(S,slot,generated,expanded,omitted);*p=S
}
void _cdid_ref_budget_store(string scalar name,real scalar slot,real scalar width)
{
    pointer() scalar p
    struct cdid_reference_inputs scalar S
    p=findexternal(name);if(p==NULL)_cdid_ref_fail("Reference staging object is absent.",498)
    S=*p;cdid_ref_encoding_budget(S,slot,width)
}

/* Calls into the reference-eligibility module, resolved only when a restricted fit runs. */
void _cdid_mrri_aux_missing_add(transmorphic a1, transmorphic a2)
{
    _cdid_mref_aux_missing_add(a1, a2)
}

void _cdid_mrri_build(transmorphic a1, transmorphic a2, transmorphic a3)
{
    _cdid_mref_build(a1, a2, a3)
}

void _cdid_mrri_set_pending()
{
    _cdid_mref_set_pending()
}

transmorphic _cdid_mrri_vars(transmorphic a1)
{
    return(_cdid_mref_vars(a1))
}

transmorphic cdid_mrri_active()
{
    return(cdid_mref_active())
}

transmorphic cdid_mrri_ctx_of(transmorphic a1)
{
    return(cdid_mref_ctx_of(a1))
}

transmorphic cdid_mrri_excluded(transmorphic a1)
{
    return(cdid_mref_excluded(a1))
}

transmorphic cdid_mrri_index()
{
    return(cdid_mref_index())
}

transmorphic cdid_mrri_policy()
{
    return(cdid_mref_policy())
}

void cdid_mrri_require_route(transmorphic a1)
{
    cdid_mref_require_route(a1)
}

string scalar _cdid_ref_xp_cg()
{
    return(cdid_xp_control_group())
}
real scalar _cdid_ref_xp_on()
{
    return(cdid_xp_requested())
}

end
