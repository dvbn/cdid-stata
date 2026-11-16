*! Private Stata long-data intake for chained IPW, selected routes, and native state.
*! Follows R cdid,
*! chained_did_request.R prepare_inputs; no duplicated nuisance/estimator formula.
version 19.5
mata:
struct cdid_ci_input {
    real scalar string_id, frame_supplied, frame_has_g, input_rows
    real colvector id, date, y, g, w, source_row, frame_id, frame_g, frame_w
    string colvector string_ids, frame_string_ids
    real matrix X, frame_X
    string rowvector x_names, frame_x_names
    real matrix Q, frame_Q
    string rowvector q_names, frame_q_names
    real colvector frame_rho, data_rho
    real scalar frame_has_rho, data_has_rho
    real scalar cl_string, frame_cl_string, in_data_cl, in_frame_cl
    real colvector cl_num, frame_cl_num
    string colvector cl_str, frame_cl_str
}
void _cdid_ci_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
void _cdid_ci_check_static_tokens(string scalar raw, string scalar role)
{
    string rowvector tok
    real scalar j
    if (strtrim(raw) == "") return
    tok = tokens(raw)
    for (j = 1; j <= cols(tok); j++) {
        if (strpos(tok[j], "#") > 0 | strpos(tok[j], ".") > 0) {
            _cdid_ci_fail("Factor-variable, interaction, and time-series operators in " + role + " are not available in this static-numeric estimator.", 459)
        }
    }
}
void _cdid_ci_data(string scalar name, string scalar mask, string scalar idname,
    string scalar tname, string scalar yname, string scalar gname,
    string scalar covars, string scalar wname,
    | string scalar qcovars, string scalar clname, string scalar rhoname)
{
    struct cdid_ci_input scalar I
    pointer() scalar p
    real colvector ix, mvec
    real matrix M_core
    real scalar all_obs
    mvec=st_data(.,mask)
    ix=vec(selectindex(mvec:!=0))
    all_obs=(rows(ix)==st_nobs() & rows(ix)>0)
    I.string_id=(substr(st_vartype(idname),1,3)=="str")
    I.id=I.frame_id=I.frame_g=I.frame_w=J(0,1,.)
    I.string_ids=I.frame_string_ids=J(0,1,"")
    if (I.string_id) {
        I.string_ids=(all_obs ? st_sdata(.,idname) : st_sdata(ix,idname))
        if (rows(ix)>0) {
            M_core=(all_obs ? st_data(.,(tname,yname,gname)) : st_data(ix,(tname,yname,gname)))
            I.date=M_core[,1]; I.y=M_core[,2]; I.g=M_core[,3]
        }
        else {
            I.date=st_data(ix,tname); I.y=st_data(ix,yname); I.g=st_data(ix,gname)
        }
    }
    else {
        if (rows(ix)>0) {
            M_core=(all_obs ? st_data(.,(idname,tname,yname,gname)) : st_data(ix,(idname,tname,yname,gname)))
            I.id=M_core[,1]; I.date=M_core[,2]; I.y=M_core[,3]; I.g=M_core[,4]
        }
        else {
            I.id=st_data(ix,idname); I.date=st_data(ix,tname); I.y=st_data(ix,yname); I.g=st_data(ix,gname)
        }
    }
    I.source_row=ix; I.input_rows=st_nobs()
    _cdid_ci_check_static_tokens(covars, "covars()")
    I.x_names=tokens(covars); I.frame_x_names=J(1,0,"")
    I.X=cols(I.x_names) ? (all_obs ? st_data(.,I.x_names) : st_data(ix,I.x_names)) : J(rows(ix),0,.)
    I.w=(wname=="" ? J(0,1,.) : (all_obs ? st_data(.,wname) : st_data(ix,wname)))
    I.frame_supplied=I.frame_has_g=0; I.frame_X=J(0,0,.)
    I.Q=I.frame_Q=J(0,0,.); I.q_names=I.frame_q_names=J(1,0,"")
    I.frame_rho=I.data_rho=J(0,1,.); I.frame_has_rho=I.data_has_rho=0
    I.cl_string=I.frame_cl_string=I.in_data_cl=I.in_frame_cl=0
    I.cl_num=I.frame_cl_num=J(0,1,.)
    I.cl_str=I.frame_cl_str=J(0,1,"")
    if (args() >= 9 & qcovars != "") {
        _cdid_ci_check_static_tokens(qcovars, "selcovars()")
        I.q_names = tokens(qcovars)
        I.Q = cols(I.q_names) ? (all_obs ? st_data(., I.q_names) : st_data(ix, I.q_names)) : J(rows(ix), 0, .)
    }
    if (args() >= 10 & clname != "") {
        I.in_data_cl = 1
        I.cl_string = (substr(st_vartype(clname), 1, 3) == "str")
        if (I.cl_string) I.cl_str = (all_obs ? st_sdata(., clname) : st_sdata(ix, clname))
        else I.cl_num = (all_obs ? st_data(., clname) : st_data(ix, clname))
    }
    if (args() >= 11 & rhoname != "") {
        I.data_has_rho = 1
        I.data_rho = (all_obs ? st_data(., rhoname) : st_data(ix, rhoname))
    }
    p=crexternal(name)
    if (p==NULL) _cdid_ci_fail("Transient input name is already in use.",498)
    *p=I
}
// Called through Stata's frame prefix: current frame is restored even on error.
void _cdid_ci_frame(string scalar name, string scalar idname, string scalar gname,
    string scalar covars, string scalar wname,
    | string scalar qcovars, string scalar rhoname, string scalar clname)
{
    struct cdid_ci_input scalar I
    pointer() scalar p
    p=findexternal(name); I=*p
    if ((substr(st_vartype(idname),1,3)=="str")!=I.string_id)
        _cdid_ci_fail("Panel/frame ID types must match; no lossy string/numeric coercion is performed.",459)
    I.frame_supplied=1; I.frame_has_g=(gname!="")
    if (I.string_id) I.frame_string_ids=st_sdata(.,idname)
    else I.frame_id=st_data(.,idname)
    if (I.frame_has_g) I.frame_g=st_data(.,gname)
    _cdid_ci_check_static_tokens(covars, "framecovars()")
    I.frame_x_names=tokens(covars)
    if (cols(I.frame_x_names)) I.frame_X=st_data(.,I.frame_x_names)
    if (wname!="") I.frame_w=st_data(.,wname)
    if (args() >= 6 & qcovars != "") {
        _cdid_ci_check_static_tokens(qcovars, "frameselcovars()")
        I.frame_q_names = tokens(qcovars)
        if (cols(I.frame_q_names)) I.frame_Q = st_data(., I.frame_q_names)
    }
    if (args() >= 7 & rhoname != "") {
        I.frame_has_rho = 1
        I.frame_rho = st_data(., rhoname)
    }
    if (args() >= 8 & clname != "") {
        I.in_frame_cl = 1
        I.frame_cl_string = (substr(st_vartype(clname), 1, 3) == "str")
        if (I.frame_cl_string) I.frame_cl_str = st_sdata(., clname)
        else I.frame_cl_num = st_data(., clname)
    }
    *p=I
}
/* A fixed static transport ratio can come from an explicit frame or from
   the first supplied row of each unit when the public frame is data-derived.
   Outside U is not read. This is the literal R data-derived frame getter:
   later long-row rho values are not evaluated as a dated transport model. */
real colvector cdid_ci_transport(struct cdid_population scalar P,
    struct cdid_ci_input scalar I, real colvector U)
{
    real colvector rho
    real scalar u, first
    rho=J(P.N,1,1)
    if(I.frame_has_rho) {
        if(rows(I.frame_rho)!=P.N) _cdid_ci_fail("Frame selection_rho length mismatch.",459)
        for(u=1;u<=P.N;u++) if(U[u]) rho[u]=I.frame_rho[u]
    }
    else if(I.data_has_rho) {
        if(I.frame_supplied | rows(I.data_rho)!=P.R)
            _cdid_ci_fail("Data-only selection_rho must align with the derived unit roster.",459)
        for(u=1;u<=P.N;u++) if(U[u]) {
            first=P.first_row[u]
            if(first<1) _cdid_ci_fail("Derived transport roster has no supplied unit row.",498)
            rho[u]=I.data_rho[first]
        }
    }
    for(u=1;u<=P.N;u++) if(U[u]) {
        if(rho[u]>=. | rho[u]<=0)
            _cdid_ci_fail("Known selection_rho must be finite and positive on selection_universe.",459)
    }
    return(rho)
}

real scalar cdid_ci_clusters(
    struct cdid_population scalar P,
    struct cdid_ci_input scalar I,
    string scalar cluster_kind,
    real colvector cluster_numeric_id,
    string colvector cluster_string_id,
    real colvector unit_cluster,
    real scalar C)
{
    real scalar is_str, n_d, n_f, n_all, i, k, u, new_label
    real colvector all_num, ord, code, d_code, f_code, num_dict, unit_code, seen
    string colvector all_str, str_dict

    if (!I.in_data_cl & !I.in_frame_cl) {
        cluster_kind = "none"
        cluster_numeric_id = J(P.N, 1, .)
        cluster_string_id = J(P.N, 1, "")
        unit_cluster = (1::P.N)
        C = P.N
        return(0)
    }
    if (I.in_data_cl & I.in_frame_cl & (I.cl_string != I.frame_cl_string)) {
        _cdid_ci_fail("Data and frame cluster variables must have the same storage kind.", 459)
    }
    is_str = (I.in_data_cl ? I.cl_string : I.frame_cl_string)
    n_d = (I.in_data_cl ? (is_str ? rows(I.cl_str) : rows(I.cl_num)) : 0)
    n_f = (I.in_frame_cl ? (is_str ? rows(I.frame_cl_str) : rows(I.frame_cl_num)) : 0)
    n_all = n_d + n_f
    d_code = J(n_d, 1, .)
    f_code = J(n_f, 1, .)
    num_dict = J(0, 1, .)
    str_dict = J(0, 1, "")

    if (!is_str) {
        if ((I.in_data_cl & (n_d != P.R | missing(I.cl_num))) |
            (I.in_frame_cl & (n_f != P.N | missing(I.frame_cl_num)))) {
            _cdid_ci_fail("Cluster labels must be plain nonmissing values.", 459)
        }
        all_num = (I.in_data_cl ? I.cl_num : J(0, 1, .)) \
                  (I.in_frame_cl ? I.frame_cl_num : J(0, 1, .))
        all_num = all_num :+ 0
        ord = order(all_num, 1)
        code = J(n_all, 1, .)
        k = 0
        for (i = 1; i <= n_all; i++) {
            if (i == 1) new_label = 1
            else new_label = (all_num[ord[i]] != all_num[ord[i - 1]])
            if (new_label) {
                k = k + 1
                num_dict = num_dict \ all_num[ord[i]]
            }
            code[ord[i]] = k
        }
        if (I.in_data_cl) d_code = code[1..n_d]
        if (I.in_frame_cl) f_code = code[(n_d + 1)..n_all]
    }
    else {
        if ((I.in_data_cl & (n_d != P.R | any(I.cl_str :== ""))) |
            (I.in_frame_cl & (n_f != P.N | any(I.frame_cl_str :== "")))) {
            _cdid_ci_fail("String cluster labels must be non-empty.", 459)
        }
        all_str = (I.in_data_cl ? I.cl_str : J(0, 1, "")) \
                  (I.in_frame_cl ? I.frame_cl_str : J(0, 1, ""))
        ord = order(all_str, 1)
        code = J(n_all, 1, .)
        k = 0
        for (i = 1; i <= n_all; i++) {
            if (i == 1) new_label = 1
            else new_label = (all_str[ord[i]] != all_str[ord[i - 1]])
            if (new_label) {
                k = k + 1
                str_dict = str_dict \ all_str[ord[i]]
            }
            code[ord[i]] = k
        }
        if (I.in_data_cl) d_code = code[1..n_d]
        if (I.in_frame_cl) f_code = code[(n_d + 1)..n_all]
    }

    unit_code = cdid_population_clusters(P, d_code, f_code, I.in_data_cl, I.in_frame_cl)
    seen = J(k, 1, 0)
    C = 0
    unit_cluster = J(P.N, 1, .)
    for (u = 1; u <= P.N; u++) {
        i = unit_code[u]
        if (seen[i] == 0) {
            C = C + 1
            seen[i] = C
        }
        unit_cluster[u] = seen[i]
    }
    if (!is_str) {
        cluster_kind = "numeric"
        cluster_numeric_id = num_dict[unit_code]
        cluster_string_id = J(P.N, 1, "")
    }
    else {
        cluster_kind = "string"
        cluster_numeric_id = J(P.N, 1, .)
        cluster_string_id = str_dict[unit_code]
    }
    return(1)
}
void _cdid_ci_fit(string scalar name, string scalar target_matrix, string scalar control,
    real scalar anticipation, string scalar pt, real scalar ptstart, string scalar fixweights,
    real scalar noconstant, real scalar chunk, real scalar maxentries, real rowvector limits,
    string scalar token_local)
{
    struct cdid_ci_input scalar I
    struct cdid_population_ids scalar ids
    struct cdid_population scalar P
    struct cdid_population_column scalar column
    struct cdid_population_gate scalar gate
    struct cdid_population_weight scalar weight
    struct cdid_chained_ipw_opts scalar options
    struct cdid_catalog_opts scalar co
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_design scalar D
    struct cdid_chained_ipw_result scalar fit
    struct cdid_result_state scalar S
    struct cdid_state_limits scalar L
    struct cdid_agg_frame scalar F
    struct cdid_agg_shares scalar shares
    struct cdid_agg_inference scalar inf
    struct cdid_infer_flags scalar flags
    struct cdid_infer_budget scalar budget
    pointer() scalar p
    real matrix fullX, dataX, request
    real colvector cohorts, targetcoh, first_source, refs, hit
    real scalar data_hasX, frame_hasX, data_hasW, frame_hasW, q, j, u, a, K
    string scalar token, fix
    p=findexternal(name); I=*p
    if (I.string_id) ids=cdid_population_ids_string(I.string_ids,I.frame_string_ids,I.frame_supplied)
    else ids=cdid_population_ids_numeric(I.id,I.frame_id,I.frame_supplied)
    P=cdid_population_prepare(ids,I.date,I.y,I.g,I.frame_g,I.frame_has_g)
    data_hasX=(cols(I.X)>0); frame_hasX=(cols(I.frame_X)>0)
    data_hasW=(rows(I.w)>0); frame_hasW=(rows(I.frame_w)>0)
    fix=(fixweights=="static" ? "" : fixweights)
    weight=cdid_population_weights(P,I.w,I.frame_w,data_hasW,frame_hasW,fix,0)
    if (weight.reduction!="static" & weight.reduction!="first_period_static")
        _cdid_ci_fail("Private static intake does not support context-dated weight design blocks.",459)
    if (data_hasX & frame_hasX & cols(I.X)!=cols(I.frame_X))
        _cdid_ci_fail("Expanded data/frame covariate columns must match in order and dimension.",459)
    q=(frame_hasX ? cols(I.frame_X) : cols(I.X))
    fullX=J(P.N,q,.)
    for (j=1;j<=q;j++) {
        column=cdid_population_column_numeric(P,"expanded column","treatment",
            (data_hasX ? I.X[,j] : J(0,1,.)),(frame_hasX ? I.frame_X[,j] : J(0,1,.)),
            data_hasX,frame_hasX,0,1)
        fullX[,j]=column.unit_value
    }
    dataX=(data_hasX ? I.X : J(0,q,.))
    if (!noconstant) {
        fullX=(J(P.N,1,1),fullX)
        if (data_hasX) dataX=(J(P.R,1,1),dataX)
    }
    gate=cdid_population_static_gate(P,fullX,dataX,J(P.N,1,1),data_hasX,0,0,1)
    if (!gate.accepted) _cdid_ci_fail(gate.message,459)
    options=cdid_chained_ipw_default_opts()
    options.est_method="ipw"; options.gmm_weighting="identity"; options.selection="none"
    options.dated_design=0; options.pt_assumption=pt; options.pt_start=ptstart
    options.fix_weights=fixweights; options.chunk_size=chunk; options.max_matrix_entries=maxentries
    cohorts=cdid_cipw_cohorts(P); targetcoh=cdid_cipw_targets(P)
    request=st_matrix(target_matrix)
    co=cdid_catalog_default_opts(); co.conventions=options.conventions
    co.mem_limit_bytes=min((2^30,8*maxentries)); co.work_cap=options.work_cap
    pairs=cdid_cipw_pairs_opts(P,cohorts,co)
    C=cdid_catalog_build_opts(P.calendar,cohorts,targetcoh,anticipation,pairs.occupancy,control,request,co)
    D=cdid_chained_ipw_build_design(C,request,pt,ptstart)
    fit=cdid_chained_ipw_fit_opts(P,C,D,fullX,weight.owner_weight,request,options)
    K=rows(fit.att)
    F.n=F.C=P.N; F.clustered=0; F.unit_cluster=(1::P.N)
    F.cohort=P.cohort; F.cohort_known=P.cohort_known; F.w=weight.V; F.share_law=""
    F.frame_id=(I.string_id ? ids.string_id_value : strofreal(ids.numeric_id,"%24.17g"))
    budget=cdid_infer_budget_create(P.N,P.N,max((1,K)),1,min((P.N,chunk)),maxentries,
        1,P.N,maxentries,8*maxentries)
    flags=cdid_infer_resolve_flags(0,0,"cells")
    inf=_cdid_agg_supported_inference(fit.att,fit.target_if,fit.target_available,F,.05,flags,
        0,J(0,0,.),budget,"no requested target is supported")
    if (rows(fit.estimable_target_indices) & inf.inference_status!="computed")
        _cdid_ci_fail("cdid_ipw_pipeline_scale_error: "+inf.reason,498)
    shares=_cdid_agg_shares(F,targetcoh)
    S=cdid_state_empty()
    S.id_kind=(I.string_id ? "string" : "numeric"); S.cluster_kind="none"
    S.point=fit.att; S.available=fit.target_available; S.role=J(K,1,3)
    S.reason=fit.target_unavailable_reason; S.IF=fit.target_if
    S.V_estimator=fit.V_estimator; S.V_asymptotic=fit.V_analytical
    refs=J(K,1,.)
    for (j=1;j<=K;j++) {
        hit=vec(selectindex(C.tg:==fit.target_g_index[j]))
        if (rows(hit)) if (C.tg_ref[hit[1]]>0) refs[j]=C.calendar[C.tg_ref[hit[1]]]
    }
    S.ledger=(fit.target_ledger,refs)
    S.ledger_columns=("g","time","g_index","t_index","catalogue_supported","touching_rows",
        "catalogue_reason","in_design","vertex","parameter","design_status","available","point","se","reference")
    S.unit_numeric_id=(I.string_id ? (1::P.N) : ids.numeric_id)
    S.unit_string_id=(I.string_id ? ids.string_id_value : J(P.N,1,""))
    S.cohort=P.cohort; S.cohort_known=P.cohort_known; S.weight=weight.V
    first_source=J(P.N,1,.)
    for (u=1;u<=P.N;u++) if (P.first_row[u]>0) first_source[u]=I.source_row[P.first_row[u]]
    S.unit_metadata=(P.row_count,P.observed_count,P.A,P.B,P.frame_only,P.F_to_A,P.F_to_B,first_source,weight.owner_weight)
    S.unit_columns=("supplied_rows","observed_rows","A","B","frame_only","F_to_A","F_to_B","first_source_row","owner_weight")
    S.cluster_numeric_id=J(P.N,1,.); S.cluster_string_id=J(P.N,1,"")
    S.shares=(shares.cohorts,shares.p,shares.mass_complete,J(rows(targetcoh),1,shares.shares_known))
    S.share_columns=("cohort","share","mass_complete","shares_known")
    S.share_IF=(shares.shares_known ? shares.P : J(0,rows(targetcoh),.))
    S.share_IF_columns=strofreal(targetcoh',"%24.17g")
    S.calendar=P.calendar
    S.numeric_settings=(anticipation\ptstart\.05\options.glm_epsilon\options.glm_maxit\options.condition_trigger\
        options.ridge_relative\chunk\maxentries\P.N\P.R\P.n_A\P.n_B\I.input_rows\fit.p_covariates\noconstant)
    S.numeric_setting_names=("anticipation"\"pt_start"\"alpha"\"glm_epsilon"\"glm_maxit"\"condition_trigger"\
        "ridge_relative"\"chunk"\"max_matrix_entries"\"N_F"\"supplied_rows"\"n_A"\"n_B"\"input_rows"\"p_covariates"\"noconstant")
    S.text_settings=("method","ipw"\"gmm_weighting","identity"\"design","static expanded numeric columns"\
        "data_covars",invtokens(I.x_names)\"frame_covars",invtokens(I.frame_x_names)\
        "inference","analytical pointwise"\"selection","none"\"control_group",control\"pt_assumption",pt\
        "fix_weights_requested",fixweights\"weight_source",weight.source\"weight_reduction",weight.reduction\
        "share_law","static weighted cohort shares"\"share_status",(shares.shares_known ? "known" : shares.share_reason)\
        "role_convention","3=explicit requested post target; no pretrend/reference rows fabricated"\
        "status",fit.status\"inference_status",inf.inference_status)
    S.source_versions=("R_cdid","1.3.0.0"\
        "engine","inverse probability weighting"\"state_format","cdid-native-result-state/v1"\
        "population","all units of the unit table"\"aggregation","cohort shares and inference")
    L.payload_slots=limits[1]; L.payload_string_bytes=limits[2]
    L.registry_slots=limits[3]; L.registry_string_bytes=limits[4]
    st_local("ci_status",fit.status); st_local("ci_N",strofreal(P.N,"%24.17g"))
    st_local("ci_R",strofreal(P.R,"%24.17g")); st_local("ci_input_R",strofreal(I.input_rows,"%24.17g"))
    st_local("ci_missing_y",strofreal(sum(I.y:>=.),"%24.17g"))
    st_local("ci_A",strofreal(P.n_A,"%24.17g")); st_local("ci_B",strofreal(P.n_B,"%24.17g"))
    st_local("ci_frameonly",strofreal(P.n_frame_only,"%24.17g")); st_local("ci_K",strofreal(K,"%24.17g"))
    token=cdid_state_put(S,L)
    st_local(token_local,token)
}
end
