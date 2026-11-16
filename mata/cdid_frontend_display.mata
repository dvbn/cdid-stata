*! Replay the stored complete scientific report without refitting or RNG draws.
version 19.5
mata:

/* Plain-language text for the status and reason codes stored with results. */
string scalar cdid_text_reason(string scalar code)
{
    string scalar c
    c = strtrim(code)
    if (c == "") return("")
    if (c == "no_declared_nonreference_pretrend_comparison") return("no pre-treatment cells to test")
    if (c == "declared_pretrend_unavailable_no_available_only_subset")
        return("some pre-treatment cells are unavailable, and the test is not computed on a subset of them")
    if (c == "declared_restrictions_exceed_centered_covariance_rank_bound")
        return("more pre-treatment cells than the rank of their estimated covariance")
    if (c == "declared_pretrend_estimates_or_joint_covariance_not_finite")
        return("the pre-treatment estimates or their covariance are not finite")
    if (c == "joint_pretrend_covariance_not_symmetric") return("the covariance of the pre-treatment cells is not symmetric")
    if (c == "full_joint_covariance_has_no_usable_Cholesky_no_pseudoinverse" |
        c == "full_joint_covariance_fails_existing_Cholesky_conditioning_policy" |
        c == "full_joint_covariance_solve_unusable" | c == "source pre-treatment covariance is singular")
        return("the covariance of the pre-treatment cells is singular")
    if (c == "full_set_Wald_not_finite_and_nonnegative") return("the Wald statistic is not finite")
    if (c == "influence functions not requested" | c == "compute_inffunc=FALSE" |
        c == "compute_inffunc_false_point_estimates_only" | c == "point estimates only; no IF or multiplier work")
        return("point estimates only (noinffunc)")
    if (c == "source state has no stored full-family pretest") return("the stored results contain no pretrend test")
    if (c == "source omits Wald: requested extra clustering absent from source analytical covariance")
        return("not computed with cluster(), as in R's did package")
    if (c == "no pre-treatment cells remain after source analytic variance filter")
        return("no pre-treatment cells with positive variance")
    if (c == "source pre-treatment covariance contains NA") return("the covariance of the pre-treatment cells has missing values")
    if (c == "source ordinary filtered Wald") return("pre-treatment cells with positive variance, as in R's did package")
    if (c == "nonconstant_targets_require_positive_marginal_inference_scales")
        return("standard errors are zero for targets that vary")
    if (c == "complete_multiplier_family_refused_or_missing" | c == "complete_multiplier_family_order_or_membership_mismatch")
        return("the multiplier draws for the complete set of targets are not available")
    if (c == "nonfinite_inference_intervals") return("the confidence intervals are not finite")
    if (c == "no_supported_reporting_targets") return("no target can be estimated")
    if (c == "empty_selection_no_draws" | c == "empty_supported_selection") return("no target selected")
    return(subinstr(c, "_", " "))
}

/* Control cohorts as listed in a reason: "Inf" is the never-treated group. */
string scalar _cdid_text_pool(string scalar pool)
{
    string rowvector p
    string scalar out
    real scalar i
    p = tokens(subinstr(pool, ",", " "))
    if (cols(p) == 0) return("none")
    out = ""
    for (i = 1; i <= cols(p); i++) out = out + (i > 1 ? ", " : "") + (p[i] == "Inf" ? "never treated" : p[i])
    return(out)
}

/* Plain statement of the control cohorts available at a date. */
string scalar _cdid_text_pooltxt(string scalar r)
{
    if (ustrregexm(r, "clean control pool at t=([^ ]+) under the declared control group: \{([^}]*)\}")) {
        if (strtrim(ustrregexs(2)) == "") return("; no control cohort is available at date " + ustrregexs(1))
        return("; control cohorts at date " + ustrregexs(1) + ": " + _cdid_text_pool(ustrregexs(2)))
    }
    return("")
}

/* Plain statement for a short status code. */
string scalar _cdid_text_code(string scalar c)
{
    if (c == "unsupported_pair") return("no units of the cohort or of its controls are observed at both dates of this comparison")
    if (c == "nonfinite_point") return("the comparison estimate is not finite")
    if (c == "no_admissible_reference") return("the cohort has no study date before its treatment date minus the anticipation periods, so it has no reference date")
    if (c == "no_target_pair_occupancy") return("no unit of the cohort is observed at both dates of any comparison")
    if (c == "no_nonempty_control_pool") return("no control cohort is available")
    if (c == "no_target_pair_occupancy_in_eligibility_context") return("no unit of the cohort with the required reference covariates is observed at both dates of any comparison")
    if (c == "nonpositive_share_mass") return("the cohort shares sum to zero")
    if (c == "intended_target") return("")
    if (c == "empty_intended_set") return("no cell belongs to this aggregate")
    if (c == "callback_error") return("the callback stopped with an error")
    if (c == "source_empty_group_period") return("the cohort or its controls have no observations at these dates")
    return(cdid_text_reason(c))
}

/* Role of a cell in words. */
string scalar cdid_text_role(string scalar role)
{
    if (role == "pre") return("pre-treatment")
    if (role == "post") return("post-treatment")
    return(role)
}

/* Label of an aggregate in words: event_-1, group_2005, calendar_2010, ATT_simple. */
string scalar cdid_text_bin(string scalar label)
{
    string scalar s
    s = strtrim(label)
    if (s == "ATT_simple") return("simple average")
    if (substr(s, 1, 6) == "event_") return("event time " + substr(s, 7, .))
    if (substr(s, 1, 6) == "group_") return("cohort " + substr(s, 7, .))
    if (substr(s, 1, 9) == "calendar_") return("date " + substr(s, 10, .))
    return(s)
}

/* Lists of aggregate labels in words. */
string scalar _cdid_text_bins(string scalar list)
{
    string rowvector a
    string scalar out
    real scalar i
    a = tokens(subinstr(list, ",", " "))
    out = ""
    for (i = 1; i <= cols(a); i++) out = out + (i > 1 ? ", " : "") + cdid_text_bin(a[i])
    return(out)
}

/* Cells named in a reason, such as pre(g=2010,t=2001) [unavailable_pretrend], in words. */
string scalar _cdid_text_cells(string scalar s)
{
    string scalar r, d
    d = char(36)
    r = ustrregexra(s, " \[[a-z0-9_]+\]", "")
    r = ustrregexra(r, "\bATT\(g=([^,()]+),t=([^()]+)\)", "cohort " + d + "1 at date " + d + "2 (post-treatment)")
    r = ustrregexra(r, "\bpre\(g=([^,()]+),t=([^()]+)\)", "cohort " + d + "1 at date " + d + "2 (pre-treatment)")
    r = ustrregexra(r, "\breference\(g=([^,()]+),t=([^()]+)\)", "cohort " + d + "1 at date " + d + "2 (reference)")
    return(r)
}

/* Plain-language reason for the availability of a cell or aggregate. Stored reasons
   follow R's support_reason texts; this rewrites them for tables and exports with the
   same information. Reasons that are already plain are returned unchanged. */
string scalar cdid_text_unavailable(string scalar reason)
{
    string scalar r, code, rest
    real scalar k, guard
    r = strtrim(reason)
    if (r == "") return("")
    r = subinstr(r, `"missing_reference = "exclude""', "missingreference(exclude)")
    code = ""
    for (guard = 1; guard <= 5; guard++) {
        if (strpos(r, "Unavailable: ") == 1) {
            r = strtrim(substr(r, 14, .))
            continue
        }
        if (ustrregexm(r, "^([a-z0-9]+(_[a-z0-9]+)+)( \([^)]*\))?:(.*)$")) {
            code = ustrregexs(1)
            rest = strtrim(ustrregexs(4))
            if (code == "empty_catalog") {
                r = (rest == "" ? "no comparison is available for this cohort" : "no comparison is available for this cohort: " + _cdid_text_code(rest))
                return(r)
            }
            if (code == "restricted_target") return("with narm, the weights are renormalized over the available cells")
            if (code == "missing_inner_aggregates") return("some aggregates are unavailable (specify narm to average the available ones): " + _cdid_text_bins(rest))
            if (code == "empty_restricted_outer_set") return("every aggregate is unavailable: " + _cdid_text_bins(subinstr(rest, "every intended bin is unavailable: ", "")))
            if (code == "missing_targets") return("some cells of this aggregate are unavailable (specify narm to average the available ones): " + _cdid_text_cells(rest))
            if (code == "empty_restricted_set") return("every cell of this aggregate is unavailable: " + _cdid_text_cells(subinstr(rest, "every intended cell is unavailable: ", "")))
            if (code == "missing_first_period_weight") return("the first-period sampling weights of some cohorts are missing: " + rest)
            if (code == "no_outer_bins") return("there is no aggregate at a nonnegative event time")
            if (code == "no_declared_post_cells") return("there are no post-treatment cells")
            if (rest == "") return(_cdid_text_code(code))
            r = rest
            continue
        }
        if (ustrregexm(r, "^[a-z0-9]+(_[a-z0-9]+)+$")) return(_cdid_text_code(r))
        break
    }
    if (ustrregexm(r, "^no stored calendar date t < g - anticipation for cohort g=(.+)$"))
        return("cohort " + ustrregexs(1) + " has no study date before its treatment date minus the anticipation periods, so it has no reference date")
    if (ustrregexm(r, "^no cohort of the universe \{([^}]*)\} other than g=([^ ]+) is untreated at date ([^ ]+) with anticipation ([^ ]+) \(comparison dates ([^ ]+)->([^ )]+)\)"))
        return("no control cohort is untreated at date " + ustrregexs(3) + (ustrregexs(4) != "0" ? " (allowing for " + ustrregexs(4) + " anticipation periods)" : "") +
            ", so the comparison of dates " + ustrregexs(5) + " and " + ustrregexs(6) + " has no controls; cohorts in the data: " + _cdid_text_pool(ustrregexs(1)))
    if (ustrregexm(r, "^no retained comparison catalog row exists for cohort g=([^ ;]+)"))
        return("no comparison of cohort " + ustrregexs(1) + " with control cohorts is observed" + _cdid_text_pooltxt(r))
    if (strpos(r, "cell is not represented in the structural design") == 1)
        return("no comparison includes this cell" + _cdid_text_pooltxt(r))
    if (ustrregexm(r, "^no retained comparison row has \(g=([^,]+),t=([^)]+)\) as a qualified post endpoint"))
        return("no comparison of cohort " + ustrregexs(1) + " with valid controls ends at date " + ustrregexs(2) + _cdid_text_pooltxt(r))
    if (strpos(r, "cell lies in design component") == 1)
        return("the comparisons that include this cell are not connected to the cohort's reference date, so its effect is not identified" + _cdid_text_pooltxt(r))
    if (r == "level anchored to the cohort reference through retained comparison rows (prop:graph)")
        return("identified by comparisons connected to the cohort's reference date")
    if (r == "level anchored to the reference or a represented PT-all pre-date through retained comparison rows (prop:graph)")
        return("identified by comparisons connected to the cohort's reference date or to a pre-treatment date under ptassumption(all)")
    if (ustrregexm(r, "^direct auxiliary comparison ([^ ]+)->([^ ]+) with controls \{([^}]*)\}"))
        return("pre-treatment comparison of dates " + ustrregexs(1) + " and " + ustrregexs(2) + " with control cohorts: " + _cdid_text_pool(ustrregexs(3)))
    if (r == "universal base-period reference cell: exact zero with zero influence function")
        return("reference cell: zero by construction (baseperiod(universal))")
    if (strpos(r, "unsupported auxiliary comparison:") == 1)
        return("pre-treatment comparison unavailable: " + cdid_text_unavailable(substr(r, 34, .)))
    if (r == "source ATT unavailable") return("the effect could not be estimated for this cell")
    if (r == "weights renormalized over the supported subset") return("with narm, the weights are renormalized over the available cells")
    return(_cdid_text_cells(r))
}

string colvector cdid_text_unavailable_v(string colvector reason)
{
    string colvector out
    real scalar i
    out = reason
    for (i = 1; i <= rows(out); i++) out[i] = cdid_text_unavailable(out[i])
    return(out)
}

/* Plain-language description of how bstrap() and cband() were applied. */
string scalar cdid_text_inference_rule(string scalar law)
{
    if (law == "explicit_bstrap_and_cband") return("bstrap() and cband() applied as specified")
    if (law == "private_bootstrap_iff_cband") return("multiplier bootstrap and simultaneous band together, as requested by cband()")
    if (law == "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band")
        return("bstrap(0) cband(1): analytic standard errors with a multiplier-bootstrap simultaneous band, as R's aggte()")
    if (law == "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws")
        return("bstrap(0) cband(1): analytic standard errors with pointwise intervals, without multiplier draws")
    if (law == "compute_inffunc_false_point_estimates_only") return("point estimates only (noinffunc)")
    if (law == "empty_selection_no_draws") return("no target selected; no inference")
    if (law == "no_supported_reporting_targets") return("no target can be estimated; no inference")
    return(law == "" ? "" : subinstr(law, "_", " "))
}

/* Stores variable names with a result state so that saved results report them
   when loaded. The new token replaces the given one. */
string scalar cdid_state_bind_text(string scalar token, string scalar request,
    string scalar limits_name, string matrix pairs)
{
    struct cdid_result_state scalar P
    struct cdid_bridge_request scalar R
    struct cdid_state_limits scalar L
    pointer() scalar ptr
    string scalar next_token
    real scalar j
    if (request != "") {
        ptr = findexternal(request)
        if (ptr == NULL) _cdid_rb_fail("cdid could not find its estimation request.", 111)
        R = *ptr
        L = R.limits
    }
    else {
        ptr = findexternal(limits_name)
        if (ptr == NULL) _cdid_rb_fail("cdid could not find its storage limits.", 111)
        L = *ptr
    }
    P = cdid_state_get(token)
    for (j = 1; j <= rows(pairs); j++) {
        if (pairs[j, 2] != "") P.text_settings = P.text_settings \ pairs[j, .]
    }
    next_token = cdid_state_put(P, L)
    cdid_state_drop(token)
    return(next_token)
}

/* Variable names and the command line stored with a result state. */
void cdid_post_state_names(string scalar token)
{
    struct cdid_result_state scalar P
    string scalar cl
    P = cdid_state_get(token)
    st_local("state_ivar", _cdid_ca_opt_text(P, "id_name", _cdid_ca_opt_text(P, "ivar_name", "")))
    cl = _cdid_ca_opt_text(P, "cluster_name", _cdid_ca_opt_text(P, "cdid_bridge.cluster_name", ""))
    if (cl == "") cl = _cdid_ca_opt_text(P, "ord.cluster_name", "")
    st_local("state_clustvar", cl)
    st_local("state_depvar", _cdid_ca_opt_text(P, "depvar", ""))
    st_local("state_time", _cdid_ca_opt_text(P, "time_name", ""))
    st_local("state_gvar", _cdid_ca_opt_text(P, "gvar_name", ""))
    st_local("state_covariates", _cdid_ca_opt_text(P, "covariate_names", ""))
    st_local("state_method", _cdid_ca_opt_text(P, "method", _cdid_ca_opt_text(P, "ord.method", "")))
    st_local("state_cmdline", _cdid_ca_opt_text(P, "cmdline", ""))
}

/* Inference labels posted as e(vce), e(vcetype), e(band) and e(inference_rule). */
void cdid_post_inference_labels(string scalar law, real scalar bstrap, real scalar cband,
    real scalar point_only, string scalar endpoint, string scalar mode,
    real scalar critical, real scalar normal)
{
    struct cdid_infer_flags scalar f
    real scalar boot, simul
    string scalar ep
    if (point_only == 1 | law == "compute_inffunc_false_point_estimates_only" |
        law == "empty_selection_no_draws" | law == "no_supported_reporting_targets") {
        st_local("vce", "none"); st_local("vcetype", ""); st_local("band", "none")
        st_local("rule", cdid_text_inference_rule(law == "" ? "compute_inffunc_false_point_estimates_only" : law))
        return
    }
    if (mode == "ordinary" | law == "") {
        boot = (bstrap == 1)
        simul = (critical < . & normal < . & critical > normal * (1 + 1e-12))
        st_local("rule", "bstrap() and cband() applied as specified")
    }
    else {
        ep = (endpoint == "cells" ? "cells" : (endpoint == "simple" ? "simple" : "aggte_bins"))
        f = cdid_infer_resolve_flags((bstrap == 0 | bstrap == 1 ? bstrap : .), (cband == 1), ep)
        boot = f.bootstrap; simul = f.simultaneous
        st_local("rule", cdid_text_inference_rule(law))
    }
    st_local("vce", (boot ? "bootstrap" : "analytic"))
    st_local("vcetype", (boot ? "Bootstrap" : ""))
    st_local("band", (simul ? "simultaneous" : "pointwise"))
}
void _cdid_disp_rule(real scalar width)
{
    printf("{txt}%s\n", width * "-")
}

/* Header line: description on the left, a labeled count on the right. */
void _cdid_disp_line(string scalar left, string scalar label, real scalar value)
{
    if (label == "") printf("{txt}%s\n", left)
    else if (udstrlen(left) > 44) {
        printf("{txt}%s\n", left)
        printf("{txt}%-44s %-22s= {res}%9s\n", "", label, (value >= . ? "." : strofreal(value, "%9.0fc")))
    }
    else printf("{txt}%-44s %-22s= {res}%9s\n", left, label, (value >= . ? "." : strofreal(value, "%9.0fc")))
}

/* Note below a table, wrapped at the line width. */
void _cdid_disp_note(string scalar text)
{
    if (strtrim(text) != "") printf("{p 0 6 2}{txt}%s{p_end}\n", text)
}

/* Reason for an unavailable estimate, without its leading status code. */
string scalar _cdid_disp_reason(string scalar reason)
{
    return(cdid_text_unavailable(reason))
}

void _cdid_disp_unavailable(string scalar reason)
{
    string scalar r
    r = _cdid_disp_reason(reason)
    printf("{p 4 6 2}{txt}Unavailable: %s{p_end}\n", (r == "" ? "not estimable" : r))
}

/* Plain label for a setting value. */
string scalar _cdid_disp_label(string scalar v)
{
    if (v == "nevertreated") return("never treated")
    if (v == "notyettreated") return("not yet treated")
    if (v == "dr") return("doubly robust (dr)")
    if (v == "ipw") return("inverse probability weighting (ipw)")
    if (v == "reg") return("outcome regression (reg)")
    if (v == "custom") return("user-supplied (custom)")
    return(v)
}

string scalar _cdid_disp_interval_head(real scalar level, real scalar simultaneous)
{
    return(sprintf("[%g%% %s]", level, (simultaneous ? "simult. band" : "conf. interval")))
}

/* e() scalar of the current results, missing when not posted. */
real scalar _cdid_disp_e(string scalar name)
{
    real scalar v
    if (st_global("e(cmd)") != "cdid") return(.)
    v = st_numscalar("e(" + name + ")")
    return(rows(v) == 1 ? v : .)
}

/* Note on standard errors and intervals. */
void _cdid_disp_inference(real scalar pointonly, real scalar bootstrap, real scalar simultaneous,
    real scalar level, string scalar what)
{
    real scalar B
    string scalar se, iv
    if (pointonly) {
        _cdid_disp_note("Point estimates only (noinffunc): no standard errors, intervals or pretrend test.")
        return
    }
    B = _cdid_disp_e("B_requested")
    if (bootstrap) se = "multiplier bootstrap (interquartile range" + (B < . ? " of " + strofreal(B) + " draws" : "") + ")"
    else se = "analytic"
    if (simultaneous) iv = sprintf("simultaneous %g%% band over the available %s", level, what)
    else iv = sprintf("pointwise %g%% intervals", level)
    _cdid_disp_note("Std. err.: " + se + "; " + iv + "." + (bootstrap ? " e(V) holds the analytic covariance." : ""))
}

real scalar _cdid_disp_opt_number(struct cdid_result_state scalar P, string scalar key)
{
    real colvector ix
    if (rows(P.numeric_setting_names) == 0) return(.)
    ix = _cdid_rb_index(P.numeric_setting_names :== "cdid_bridge." + key)
    return(rows(ix) == 1 ? P.numeric_settings[ix[1]] : .)
}

void _cdid_disp_wald(struct cdid_result_state scalar P, string scalar status)
{
    real scalar st, df, p
    string scalar reason
    st = _cdid_disp_opt_number(P, "wald_statistic")
    df = _cdid_disp_opt_number(P, "wald_df")
    p = _cdid_disp_opt_number(P, "wald_p")
    _cdid_disp_wald_line(status, st, df, p, _cdid_ca_opt_text(P, "cdid_bridge.wald_reason", ""))
}

void _cdid_disp_wald_line(string scalar status, real scalar st, real scalar df, real scalar p, string scalar reason)
{
    string scalar r
    if (status == "computed" & st < . & df < . & p < .) {
        printf("{txt}Pretrend Wald test (pre-treatment cells): chi2(%g) = {res}%s{txt}, Prob > chi2 = {res}%6.4f\n",
            df, strtrim(strofreal(st, "%10.2f")), p)
        return
    }
    r = cdid_text_reason(reason)
    if (status == "not_applicable") _cdid_disp_note("Pretrend Wald test: not applicable" + (r != "" ? "; " + r : "") + ".")
    else _cdid_disp_note("Pretrend Wald test: not computed" + (r != "" ? "; " + r : "") + ".")
}

/* Header shared by chained cell and aggregation results. */
void _cdid_disp_chained_header(struct cdid_result_state scalar P, string scalar title,
    string scalar count_label, real scalar count, real scalar N)
{
    string scalar method, weighting, sel, cl, line, base
    real scalar C, preset, crosspair
    method = _cdid_ca_opt_text(P, "method", st_global("e(method)"))
    weighting = _cdid_ca_opt_text(P, "gmm_weighting", st_global("e(gmm_weighting)"))
    sel = _cdid_ca_opt_text(P, "selection", "none")
    preset = 0
    crosspair = 0
    if (rows(P.text_settings)) {
        preset = any((P.text_settings[, 1] :== "preset") :& (P.text_settings[, 2] :== "cohort_time_fe"))
        crosspair = any((P.text_settings[, 1] :== "cross_pair") :& (P.text_settings[, 2] :== "1"))
    }
    cl = st_global("e(clustvar)")
    C = _cdid_disp_e("N_clust")
    printf("\n")
    _cdid_disp_line(title, "Number of units", N)
    if (preset & method == "reg") line = "Method: linear regression with fixed effects"
    else line = (method != "" ? "Method: " + _cdid_disp_label(method) : "")
    _cdid_disp_line(line, count_label, count)
    if (preset & method == "reg") line = "Cohort-by-date fixed effects (within)"
    else if (weighting != "") line = (crosspair ? "Cross-pair comparisons; " : "") + "GMM weighting: " + weighting
    else line = (crosspair ? "Cross-pair comparisons" : "")
    if (cl != "" & C < .) _cdid_disp_line(line, "Number of clusters", C)
    else if (line != "") _cdid_disp_line(line, "", .)
    line = _cdid_disp_label(_cdid_ca_opt_text(P, "control_group", st_global("e(control_group)")))
    base = _cdid_ca_opt_text(P, "base_period", st_global("e(base_period)"))
    if (line != "" & base != "") _cdid_disp_line("Control group: " + line + "; base period: " + base, "", .)
    else if (line != "") _cdid_disp_line("Control group: " + line, "", .)
    if (sel != "" & sel != "none") {
        if (_cdid_ca_opt_text(P, "selection_est_method", "ipw") == "imputation")
            line = "Missing outcomes: " + sel + ", change imputation (" + _cdid_ca_opt_text(P, "imputation_pooling", "separate") + " pooling)"
        else line = "Missing outcomes: " + sel + ", " + _cdid_ca_opt_text(P, "selection_est_method", "ipw") + " (" +
            _cdid_ca_opt_text(P, "selection_model", "pair_logit") + ", " + _cdid_ca_opt_text(P, "selection_pooling", "separate") +
            " pooling, universe " + _cdid_ca_opt_text(P, "selection_universe", "F") + ")"
        _cdid_disp_note(line)
    }
    if (cl != "" & C < .) printf("{txt}(Std. err. adjusted for %s clusters in %s)\n", strofreal(C, "%12.0fc"), cl)
}

void cdid_frontend_display(string scalar token)
{
    struct cdid_result_state scalar P
    struct cdid_infer_flags scalar flags
    real scalar K,j,ig,it,ise,ilo,ihi,alpha,pointonly,level
    string scalar role,title
    P=cdid_state_get(token);K=rows(P.point)
    ig=_cdid_rb_column(P.ledger_columns,"g");it=_cdid_rb_column(P.ledger_columns,"time")
    ise=_cdid_rb_column(P.ledger_columns,"cdid_bridge.cells.se")
    ilo=_cdid_rb_column(P.ledger_columns,"cdid_bridge.cells.band_lower")
    ihi=_cdid_rb_column(P.ledger_columns,"cdid_bridge.cells.band_upper")
    alpha=_cdid_rb_stored_number(P,"cells.alpha");pointonly=_cdid_rb_stored_number(P,"point_only")
    level=100*(1-alpha)
    flags=cdid_infer_resolve_flags(_cdid_rb_stored_number(P,"cells.bstrap"),_cdid_rb_stored_number(P,"cells.cband"),"cells")
    title="Chained difference-in-differences"
    _cdid_disp_chained_header(P,title,"Number of cells",K,rows(P.cohort))
    _cdid_disp_rule(78)
    printf("{txt}%8s %8s  %-9s %10s %10s  %23s\n","Group","Time","Role","ATT","Std. err.",
        (pointonly ? "" : _cdid_disp_interval_head(level,flags.simultaneous)))
    _cdid_disp_rule(78)
    for(j=1;j<=K;j++) {
        role=(P.role[j]==3 ? "post" : (P.role[j]==2 ? "pre" : "reference"))
        printf("{txt}%8.0g %8.0g  %-9s {res}%10.7g %10.7g  %11.7g %11.7g\n",P.ledger[j,ig],P.ledger[j,it],role,
            P.point[j],P.ledger[j,ise],P.ledger[j,ilo],P.ledger[j,ihi])
        if(!P.available[j]) _cdid_disp_unavailable(P.reason[j])
    }
    _cdid_disp_rule(78)
    _cdid_disp_inference(pointonly,flags.bootstrap,flags.simultaneous,level,"cells")
    if (any(P.role:==1)) _cdid_disp_note("Reference cells are zero by construction.")
    if(!pointonly) _cdid_disp_wald(P,_cdid_ca_opt_text(P,"cdid_bridge.wald_status","not_computed"))
}
end
