*! Private static numeric data intake; population/reader semantics + selection/cluster intake.
*! Loads after cdid_input_adapter.mata and cdid_public_request.mata.
version 19.5
mata:
struct cdid_public_prepared {
    struct cdid_population scalar population
    struct cdid_population_weight scalar weight
    struct cdid_population_column colvector columns
    real matrix X
    string rowvector x_names
    real colvector required_X, source_row
    string scalar method, design
    string scalar selection, selection_canonical, selection_model, selection_pooling
    string scalar selection_universe, selection_est_method, selection_justification
    string scalar inactive_arguments, gmm_covariance, imputation_pooling, outcome_pooling
    struct cdid_population_column colvector sel_columns
    real matrix Z_sel
    string rowvector sel_names
    real colvector U, sel_rho
    real scalar x_has_intercept, sel_has_intercept, selection_formula_inherited, sel_has_rho
    string scalar cluster_kind
    real colvector cluster_numeric_id, unit_cluster
    string colvector cluster_string_id
    real scalar C, clustered
}
struct cdid_public_bundle {
    struct cdid_public_prepared scalar prepared
    struct cdid_public_request scalar request
    struct cdid_public_settings scalar settings
    struct cdid_reference_inputs scalar reference
    string scalar baseline_static_status, baseline_static_reason
    string scalar selected_static_status, selected_static_reason
}
void _cdid_prep_role_fail(string scalar cls, string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
struct cdid_public_prepared scalar cdid_public_static_prepare(struct cdid_ci_input scalar I,
    string scalar method,real scalar noconstant,string scalar fixweights,
    | string scalar selection_in,string scalar selection_model_in,
      string scalar selection_pooling_in,string scalar selection_universe_in,
      string scalar selection_justification_in,string scalar selection_est_method_in,
      string scalar imputation_pooling_in,string scalar gmm_covariance_in,
      real scalar selnoconstant_in,real scalar selconstant_in,real scalar selzero_in,
      string scalar outcome_pooling_in,real scalar pre_ref_in)
{
    struct cdid_public_prepared scalar B
    struct cdid_population_ids scalar ids
    struct cdid_population_gate scalar gate, sel_gate
    struct cdid_sel_decl scalar decl
    real scalar dx,fx,q,dq,fq,qq,j,u,sel_nc,sel_c,sel_z,has_sel_design_opt,out_active,pre_ref
    real matrix dataX, dataQ
    string scalar sel,smodel,spool,suniv,sjust,smeth,imppool,outpool,gcov,inact
    if(method!="ipw" & method!="reg" & method!="dr") _cdid_pr_fail("cdid_chained_request_input_error","method() must be dr, reg, ipw, or custom.")
    if(noconstant!=0 & noconstant!=1) _cdid_pr_fail("cdid_chained_request_input_error","noconstant must be binary.")
    pre_ref = (args() >= 17 & pre_ref_in != 0 & pre_ref_in < . ? 1 : 0)
    out_active = (method == "reg" | method == "dr")
    if (out_active) {
        outpool = (args() >= 16 ? strlower(strtrim(outcome_pooling_in)) : "")
        if (outpool == "") outpool = "separate"
        if (outpool != "separate" & outpool != "pooled")
            _cdid_prep_role_fail("cdid_chained_request_outcome_pooling", "outcomepooling() must be separate or pooled.", 198)
        inact = ""
    }
    else {
        outpool = ""
        inact = "outcome_pooling"
    }
    sel = (args() >= 5 ? strlower(strtrim(selection_in)) : "")
    if (sel == "") sel = "none"
    if (sel != "none" & sel != "conditional_mtar" & sel != "corollary1")
        _cdid_prep_role_fail("cdid_chained_request_input_error", "selection() must be none or conditional_mtar.", 198)
    smodel = (args() >= 6 ? strlower(strtrim(selection_model_in)) : "")
    if (smodel == "") smodel = "pair_logit"
    spool = (args() >= 7 ? strlower(strtrim(selection_pooling_in)) : "")
    if (spool == "") spool = "separate"
    suniv = (args() >= 8 ? strupper(strtrim(selection_universe_in)) : "")
    if (suniv == "") suniv = "F"
    sjust = (args() >= 9 ? selection_justification_in : "")
    smeth = (args() >= 10 ? strlower(strtrim(selection_est_method_in)) : "")
    if (smeth == "") smeth = "ipw"
    if (smeth != "ipw" & smeth != "imputation" & smeth != "dr")
        _cdid_prep_role_fail("cdid_chained_request_input_error", "selmethod() must be ipw, imputation, or dr.", 198)
    if (smeth == "dr") decl=cdid_aug_route_gate(method,0,sel,suniv,smodel,spool,
        (I.frame_has_rho | I.data_has_rho),sjust)
    if (smeth == "ipw") {
        /* Per R chained_did_request.R: when selection_est_method == "ipw",
           imputation_pooling is an inactive argument and is never evaluated/validated. */
        B.inactive_arguments = (inact == "" ? "imputation_pooling" : inact + " imputation_pooling")
    }
    else {
        B.inactive_arguments = inact
        imppool = (args() >= 11 ? strlower(strtrim(imputation_pooling_in)) : "")
        if (imppool == "") imppool = "separate"
        if (imppool != "separate" & imppool != "within_cohort" & imppool != "across_cohorts")
            _cdid_prep_role_fail("cdid_chained_request_input_error", "imppooling() must be separate, within_cohort, or across_cohorts.", 198)
        if (smeth == "imputation") {
            if (sel == "none" | suniv != "F" | I.frame_has_rho | I.data_has_rho)
                _cdid_pr_fail("cdid_chained_request_unsupported_selection", "Pure imputation requires active conditional_mtar or corollary1, representative F, and null rho=1.")
            B.inactive_arguments = (inact == "" ? "selection_model selection_pooling selection_xformla" :
                inact + " selection_model selection_pooling selection_xformla")
        }
    }
    gcov = (args() >= 12 ? strlower(strtrim(gmm_covariance_in)) : "")
    if (gcov == "") gcov = "tikhonov"
    if (gcov != "tikhonov")
        _cdid_prep_role_fail("cdid_chained_request_input_error", "gmmcovariance() must be tikhonov.", 198)
    sel_nc = (args() >= 13 & selnoconstant_in != 0 & selnoconstant_in < . ? 1 : 0)
    sel_c = (args() >= 14 & selconstant_in != 0 & selconstant_in < . ? 1 : 0)
    sel_z = (args() >= 15 & selzero_in != 0 & selzero_in < . ? 1 : 0)
    dq = (cols(I.Q) > 0); fq = (cols(I.frame_Q) > 0)
    has_sel_design_opt = (dq | fq | sel_nc | sel_c | sel_z)

    if (smeth == "imputation") {
        /* Deliberately never normalize observation model/pooling/design. */
        decl.canonical_assumption = "conditional_mtar"
        decl.selection_model = "none"
        decl.selection_pooling = "inactive"
        decl.selection_universe = "F"
        decl.justification = sjust
    }
    else {
        decl = cdid_sel_norm_decl(
            sel,
            (sel == "none" ? has_sel_design_opt : 1),
            suniv,
            (I.frame_has_rho | I.data_has_rho),
            (strtrim(sjust) != ""),
            smodel,
            spool,
            sjust)
    }

    B.selection = sel
    B.selection_canonical = decl.canonical_assumption
    B.selection_model = decl.selection_model
    B.selection_pooling = decl.selection_pooling
    B.selection_universe = decl.selection_universe
    B.selection_est_method = smeth
    B.imputation_pooling = (smeth != "ipw" ? imppool : "")
    B.outcome_pooling = outpool
    B.selection_justification = decl.justification
    B.gmm_covariance = gcov
    B.sel_has_rho = (I.frame_has_rho | I.data_has_rho)

    if(I.string_id) ids=cdid_population_ids_string(I.string_ids,I.frame_string_ids,I.frame_supplied)
    else ids=cdid_population_ids_numeric(I.id,I.frame_id,I.frame_supplied)
    B.population=cdid_population_prepare(ids,I.date,I.y,I.g,I.frame_g,I.frame_has_g)
    if (B.selection_canonical == "conditional_mtar" & any(B.population.cohort_known :== 0))
        _cdid_pr_fail("cdid_chained_request_cohort_error", "Corollary1 requires known cohorts on every F row, including frame-only rows.")
    if (B.selection_canonical == "none" & method == "ipw" & any(B.population.cohort_known :== 0))
        _cdid_pr_fail("cdid_chained_request_cohort_error", "method(ipw) requires a known gvar() value for every unit, including population-frame units.")

    B.method=method;B.design="static_numeric";B.source_row=I.source_row
    B.x_has_intercept = (noconstant == 0)
    B.weight=cdid_population_weights(B.population,I.w,I.frame_w,rows(I.w)>0,rows(I.frame_w)>0,fixweights,0)
    if(B.weight.reduction!="static" & B.weight.reduction!="first_period_static") _cdid_pr_fail("cdid_reference_support_error","This private static path requires a context-dated weight design block reader.")
    dx=(cols(I.X)>0);fx=(cols(I.frame_X)>0);q=(fx ? cols(I.frame_X) : cols(I.X))
    if(dx & fx & cols(I.X)!=q) _cdid_pr_fail("cdid_chained_request_covariate_error","Data/frame expanded columns disagree in dimension.")
    if (dx & fx) if (any(I.x_names :!= I.frame_x_names)) _cdid_pr_fail("cdid_chained_request_covariate_error","Data/frame covariate names disagree in order.")
    B.columns=J(q,1,cdid_population_column());B.X=J(B.population.N,q,.)
    for(j=1;j<=q;j++) {
        B.columns[j]=cdid_population_column_numeric(B.population,"expanded column","treatment",(dx ? I.X[,j] : J(0,1,.)),(fx ? I.frame_X[,j] : J(0,1,.)),dx,fx,0,1)
        B.X[,j]=B.columns[j].unit_value
    }
    B.x_names=(fx ? I.frame_x_names : I.x_names);dataX=(dx ? I.X : J(0,q,.))
    if(!noconstant) {
        B.X=(J(B.population.N,1,1),B.X)
        B.x_names=("_cons",B.x_names)
        if(dx) dataX=(J(B.population.R,1,1),dataX)
    }
    /* Baseline pure IPW requires X finite on F; baseline reg/dr AND all selected
       routes (ipw/reg/dr under conditional_mtar) require X finite on A only. */
    B.required_X=((smeth!="ipw" | (method=="ipw" & B.selection_canonical=="none")) ? J(B.population.N,1,1) : B.population.A)
    /* Cross-pair comparisons under selection weights reweight the selection universe of
       each target cohort, so the covariates are required there as well. */
    if(B.selection_canonical=="conditional_mtar" & smeth=="ipw" & cdid_xp_requested())
        B.required_X=((strtrim(B.selection_universe)=="F" | strtrim(B.selection_universe)=="") ? J(B.population.N,1,1) :
            (strtrim(B.selection_universe)=="B" ? (B.population.B:|B.population.A) : B.population.A))
    /* preset(cohort_time_fe) averages the covariates over every unit of each cohort. */
    if(cdid_etwfe_requested()) B.required_X=J(B.population.N,1,1)
    gate=cdid_population_static_gate(B.population,B.X,dataX,B.required_X,dx,0,pre_ref,1)
    if(!gate.accepted) _cdid_pr_fail(gate.error_class,gate.message)
    // cp04_public_design retains NA model-matrix rows outside the fitting roster.
    if(any(B.required_X:==0) & cols(B.X)>0) B.X[vec(selectindex(B.required_X:==0)),]=J(sum(B.required_X:==0),cols(B.X),.)

    /* Universe U and selection design Z_sel (selection_xformla) */
    if (B.selection_canonical == "none" | B.selection_universe == "F") B.U = J(B.population.N, 1, 1)
    else if (B.selection_universe == "B") B.U = B.population.B
    else B.U = B.population.A
    if (sum(B.U) < 1) _cdid_pr_fail("cdid_chained_request_input_error", "Selected universe U has zero units.")

    B.sel_rho = J(B.population.N, 1, 1)
    if (B.selection_canonical == "none" | smeth == "imputation") {
        B.sel_columns = J(0, 1, cdid_population_column())
        B.Z_sel = J(B.population.N, 0, .)
        B.sel_names = J(1, 0, "")
        B.sel_has_intercept = (smeth == "imputation" ? 0 : B.x_has_intercept)
        B.selection_formula_inherited = (smeth == "imputation" ? 0 : 1)
    }
    else {
        if (sel_z & (dq | fq | sel_c))
            _cdid_prep_role_fail("cdid_chained_request_input_error", "selzero (~0) cannot be combined with selcovars(), frameselcovars(), or selconstant.", 198)
        if (sel_c & (dq | fq | sel_nc | sel_z))
            _cdid_prep_role_fail("cdid_chained_request_input_error", "selconstant (~1) cannot be combined with selcovars(), frameselcovars(), selnoconstant, or selzero.", 198)

        if (sel_z) {
            /* Explicit zero-column selection design (~0) */
            B.selection_formula_inherited = 0
            B.sel_has_intercept = 0
            B.sel_columns = J(0, 1, cdid_population_column())
            B.Z_sel = J(B.population.N, 0, .)
            B.sel_names = J(1, 0, "")
        }
        else if (sel_c) {
            /* Explicit intercept-only selection design (~1) */
            B.selection_formula_inherited = 0
            B.sel_has_intercept = 1
            B.sel_columns = J(0, 1, cdid_population_column())
            B.Z_sel = J(B.population.N, 1, 1)
            B.sel_names = ("_cons")
            if (any(B.U :== 0)) B.Z_sel[vec(selectindex(B.U :== 0)), ] = J(sum(B.U :== 0), 1, .)
        }
        else if (dq | fq) {
            /* Explicit separate selection covariates Q */
            B.selection_formula_inherited = 0
            B.sel_has_intercept = (sel_nc == 0)
            qq = (fq ? cols(I.frame_Q) : cols(I.Q))
            if (dq & fq & cols(I.Q) != qq)
                _cdid_pr_fail("cdid_chained_request_covariate_error", "Data/frame selection covariate columns disagree in dimension.")
            if (dq & fq) if (any(I.q_names :!= I.frame_q_names))
                _cdid_pr_fail("cdid_chained_request_covariate_error", "Data/frame selection covariate names disagree in order.")
            B.sel_columns = J(qq, 1, cdid_population_column())
            B.Z_sel = J(B.population.N, qq, .)
            for (j = 1; j <= qq; j++) {
                B.sel_columns[j] = cdid_population_column_numeric(B.population, "selection column", "selection",
                    (dq ? I.Q[, j] : J(0, 1, .)), (fq ? I.frame_Q[, j] : J(0, 1, .)), dq, fq, 0, 1)
                B.Z_sel[, j] = B.sel_columns[j].unit_value
            }
            B.sel_names = (fq ? I.frame_q_names : I.q_names)
            dataQ = (dq ? I.Q : J(0, qq, .))
            if (B.sel_has_intercept) {
                B.Z_sel = (J(B.population.N, 1, 1), B.Z_sel)
                B.sel_names = ("_cons", B.sel_names)
                if (dq) dataQ = (J(B.population.R, 1, 1), dataQ)
            }
            sel_gate = cdid_population_static_gate(B.population, B.Z_sel, dataQ, B.U, dq, 0, 1, 1)
            if (!sel_gate.accepted) _cdid_pr_fail(sel_gate.error_class, sel_gate.message)
            if (any(B.U :== 0) & cols(B.Z_sel) > 0)
                B.Z_sel[vec(selectindex(B.U :== 0)), ] = J(sum(B.U :== 0), cols(B.Z_sel), .)
        }
        else {
            /* Inherited from xformla (or selnoconstant override on xformla covariates) */
            B.selection_formula_inherited = (sel_nc == 0)
            B.sel_has_intercept = (sel_nc ? 0 : B.x_has_intercept)
            B.sel_columns = B.columns
            B.Z_sel = J(B.population.N, q, .)
            for (j = 1; j <= q; j++) B.Z_sel[, j] = B.columns[j].unit_value
            B.sel_names = (fx ? I.frame_x_names : I.x_names)
            dataQ = (dx ? I.X : J(0, q, .))
            if (B.sel_has_intercept) {
                B.Z_sel = (J(B.population.N, 1, 1), B.Z_sel)
                B.sel_names = ("_cons", B.sel_names)
                if (dx) dataQ = (J(B.population.R, 1, 1), dataQ)
            }
            sel_gate = cdid_population_static_gate(B.population, B.Z_sel, dataQ, B.U, dx, 0, 1, 1)
            if (!sel_gate.accepted) _cdid_pr_fail(sel_gate.error_class, sel_gate.message)
            if (any(B.U :== 0) & cols(B.Z_sel) > 0)
                B.Z_sel[vec(selectindex(B.U :== 0)), ] = J(sum(B.U :== 0), cols(B.Z_sel), .)
        }

        if (B.sel_has_rho) B.sel_rho=cdid_ci_transport(B.population,I,B.U)

    }

    B.clustered = cdid_ci_clusters(B.population, I, B.cluster_kind,
        B.cluster_numeric_id, B.cluster_string_id, B.unit_cluster, B.C)
    return(B)
}
void _cdid_public_prepare_store(string scalar input_name,string scalar output_name,
    string scalar method,real scalar noconstant,string scalar fixweights,
    real colvector targets,string scalar control,real scalar ant,string scalar base,
    string scalar pt,real scalar ptstart,real scalar max_entries,
    | string scalar selection_in,string scalar selection_model_in,
      string scalar selection_pooling_in,string scalar selection_universe_in,
      string scalar selection_justification_in,string scalar selection_est_method_in,
      string scalar imputation_pooling_in,string scalar gmm_covariance_in,
      real scalar selnoconstant_in,real scalar selconstant_in,real scalar selzero_in,
      string scalar outcome_pooling_in)
{
    struct cdid_ci_input scalar I
    struct cdid_public_bundle scalar B
    pointer() scalar p
    string scalar smeth_lc
    real scalar pre_ref
    p=findexternal(input_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Transient input is absent.")
    if(findexternal(output_name)!=NULL) _cdid_pr_fail("cdid_chained_request_input_error","Requested output name already exists; no user object is overwritten.")
    I=*p
    smeth_lc=(args()>=18 ? strlower(strtrim(selection_est_method_in)) : "ipw")
    pre_ref=((smeth_lc=="imputation" | smeth_lc=="dr") & st_local("referencechecked")=="")
    B.prepared=cdid_public_static_prepare(I,method,noconstant,fixweights,
        (args()>=13 ? selection_in : "none"),
        (args()>=14 ? selection_model_in : "pair_logit"),
        (args()>=15 ? selection_pooling_in : "separate"),
        (args()>=16 ? selection_universe_in : "F"),
        (args()>=17 ? selection_justification_in : ""),
        (args()>=18 ? selection_est_method_in : "ipw"),
        (args()>=19 ? imputation_pooling_in : ""),
        (args()>=20 ? gmm_covariance_in : "tikhonov"),
        (args()>=21 ? selnoconstant_in : 0),
        (args()>=22 ? selconstant_in : 0),
        (args()>=23 ? selzero_in : 0),
        (args()>=24 ? outcome_pooling_in : ""),
        pre_ref)
    B.request=cdid_public_request_build(B.prepared.population,targets,control,ant,base,pt,ptstart,max_entries)
    B.settings=cdid_public_defaults();B.settings.method=method;B.settings.control_group=control
    B.settings.anticipation=ant;B.settings.base_period=base;B.settings.pt_assumption=pt;B.settings.pt_start=ptstart
    B.settings.selection=B.prepared.selection
    B.settings.selection_canonical=B.prepared.selection_canonical
    B.settings.selection_model=B.prepared.selection_model
    B.settings.selection_pooling=B.prepared.selection_pooling
    B.settings.selection_universe=B.prepared.selection_universe
    B.settings.selection_est_method=B.prepared.selection_est_method
    B.settings.selection_justification=B.prepared.selection_justification
    B.settings.imputation_pooling=B.prepared.imputation_pooling
    B.settings.outcome_pooling=B.prepared.outcome_pooling
    B.settings.inactive_arguments=B.prepared.inactive_arguments
    B.settings.gmm_covariance=B.prepared.gmm_covariance
    p=crexternal(output_name)
    if(p==NULL) _cdid_pr_fail("cdid_chained_request_input_error","Cannot create requested output object.")
    *p=B
    st_local("pr_N",strofreal(B.prepared.population.N));st_local("pr_A",strofreal(B.prepared.population.n_A))
    st_local("pr_post",strofreal(rows(B.request.post)));st_local("pr_aux",strofreal(rows(B.request.aux)))
    st_local("pr_report",strofreal(rows(B.request.report)))
}
void _cdid_imp_reference_begin(string scalar bundle,string scalar popname,string scalar planname,
    string scalar rowname,real scalar cap)
{
    struct cdid_public_bundle scalar B
    struct cdid_catalog_opts scalar co
    struct cdid_pairs scalar pairs
    struct cdid_cat scalar C
    struct cdid_baseline_plan_state scalar plan
    pointer() scalar p
    B=*findexternal(bundle)
    co=cdid_catalog_default_opts();co.mem_limit_bytes=8*cap
    co.conventions="D10_D13_admissible_frame|IMP-"+B.prepared.method+"-"+B.prepared.imputation_pooling+"-F-actual-change-fixed-w"
    pairs=cdid_cipw_pairs_opts(B.prepared.population,B.request.cohorts,co)
    C=cdid_catalog_build_opts(B.request.calendar,B.request.cohorts,B.request.target_cohorts,
        B.request.anticipation,pairs.occupancy,B.request.control_group,B.request.post,co)
    plan=cdid_baseline_plan(B.prepared.population,C,B.request.aux,
        B.request.aux_K_ptr,B.request.aux_K_index,cap,co.work_cap)
    p=crexternal(popname);if(p==NULL)_cdid_ref_fail("Transient population name is in use.",498)
    *p=B.prepared.population
    p=crexternal(planname);if(p==NULL)_cdid_ref_fail("Transient plan name is in use.",498)
    *p=plan
    p=crexternal(rowname);if(p==NULL)_cdid_ref_fail("Transient source-row name is in use.",498)
    *p=B.prepared.source_row
    *findexternal(bundle)=B
}
void _cdid_imp_reference_map(string scalar name,string scalar xf,string scalar ff,
    string scalar w,string scalar fw)
{
    string rowvector x,f
    string matrix map
    real scalar j
    pointer() scalar p
    x=cdid_ref_formula_names(xf);f=cdid_ref_formula_names(ff);map=J(0,2,"")
    if(cols(f)) {
        if(cols(f)!=cols(x))_cdid_ref_fail("Data/frame raw formula mappings disagree in dimension.",459)
        for(j=1;j<=cols(x);j++)map=map\(x[j],f[j])
    }
    if(fw!="")map=map\(w,fw)
    p=crexternal(name);if(p==NULL)_cdid_ref_fail("Transient frame map name is in use.",498)
    *p=map
}
void _cdid_imp_reference_attach(string scalar bundle,string scalar refname,string scalar wname)
{
    struct cdid_public_bundle scalar B
    struct cdid_reference_inputs scalar S
    real scalar active,j,slot
    B=*findexternal(bundle);S=*findexternal(refname)
    active=any(S.kind:!=1)
    for(j=1;j<=2;j++)if(S.roles[j].active)
        active=active | strpos(S.roles[j].formula,".")>0 | strpos(S.roles[j].formula,"#")>0
    st_local("reference_active",strofreal(active))
    if(!active)return
    slot=S.roles[1].slot
    if(cols(S.designs[slot].value)==0)
        _cdid_pr_fail("cdid_chained_request_covariate_error","Covariates must be supplied on the fitting domain.")
    B.prepared.design="reference_recipes"
    B.prepared.X=J(0,cols(S.designs[slot].value),.)
    B.prepared.x_names=S.designs[slot].columns
    B.prepared.x_has_intercept=S.roles[1].intercept
    B.prepared.weight.kind=(S.weight_kind=="dated" ? "varying" : "static");B.prepared.weight.fix_weights=S.fix_weights
    B.prepared.weight.reduction="reference_recipes"
    B.prepared.weight.V=B.prepared.weight.owner_weight=S.V
    B.prepared.weight.scientific_static=J(0,1,.)
    B.prepared.weight.source="none"
    if(wname!="") {
        j=selectindex(S.raw_names:==wname)[1]
        B.prepared.weight.source=(S.has_frame[j] ? "frame" : "data")
    }
    S.population=NULL
    B.reference=S
    *findexternal(bundle)=B
}
void _cdid_aug_reference_attach(string scalar bundle,string scalar refname,string scalar wname,
    real scalar inherited)
{
    struct cdid_public_bundle scalar B
    struct cdid_reference_inputs scalar S
    real scalar slot,j,active
    B=*findexternal(bundle);S=*findexternal(refname)
    active=any(S.kind:!=1)
    for(j=1;j<=3;j++)if(S.roles[j].active)
        active=active | strpos(S.roles[j].formula,".")>0 | strpos(S.roles[j].formula,"#")>0
    st_local("reference_active",strofreal(active))
    if(!active)return
    if(S.route!="augmented" | (S.method!="ipw" & S.method!="reg" & S.method!="dr") | S.method!=B.prepared.method)
        _cdid_ref_fail("Augmented public intake has the wrong scientific route.",498)
    slot=S.roles[1].slot
    if(cols(S.designs[slot].value)==0)
        _cdid_pr_fail("cdid_chained_request_covariate_error","Covariates must be supplied on the fitting domain.")
    B.prepared.design="reference_recipes"
    B.prepared.X=J(0,cols(S.designs[slot].value),.)
    B.prepared.x_names=S.designs[slot].columns
    B.prepared.x_has_intercept=S.roles[1].intercept
    B.prepared.Z_sel=J(0,cols(S.designs[S.roles[3].slot].value),.)
    B.prepared.sel_names=S.designs[S.roles[3].slot].columns
    B.prepared.sel_has_intercept=S.roles[3].intercept
    B.prepared.selection_formula_inherited=inherited
    B.prepared.weight.kind=(S.weight_kind=="dated" ? "varying" : "static")
    B.prepared.weight.fix_weights=S.fix_weights;B.prepared.weight.reduction="reference_recipes"
    B.prepared.weight.V=B.prepared.weight.owner_weight=S.V
    B.prepared.weight.scientific_static=J(0,1,.)
    B.prepared.weight.source="none"
    if(wname!="") {
        j=selectindex(S.raw_names:==wname)[1]
        B.prepared.weight.source=(S.has_frame[j] ? "frame" : "data")
    }
    S.population=NULL;B.reference=S
    *findexternal(bundle)=B
}

end
