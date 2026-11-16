*! Build one native library from the explicit composition snapshot; current19.5.
version 19.5
set more off
args root candidate library_out
if "`candidate'"=="" | "`library_out'"=="" {
    display as error "supply the candidate source directory and a fresh library output directory"
    exit 198
}
foreach m in cdid_population cdid_catalog cdid_nuisance cdid_ipw_link cdid_or_dr_link cdid_gmm cdid_within ///
    cdid_imputation_model cdid_imputation_dated cdid_imputation_backend ///
    cdid_chained_ipw cdid_callback cdid_selection cdid_selected_link cdid_augmented_backend cdid_inference cdid_aggregation cdid_result_state ///
    cdid_baseline_provider cdid_baseline_system cdid_input_adapter cdid_result_bridge cdid_public_request ///
    cdid_reference_inputs cdid_reference_eligibility cdid_reference_catalog cdid_dated_baseline_provider cdid_selected_reference cdid_public_prepare cdid_public_rng cdid_public_report cdid_public_provider_report cdid_public_control ///
    cdid_public_support cdid_frontend_types cdid_ipw_auxiliary cdid_frontend_ipw cdid_nuisance_report cdid_nuisance_augmented cdid_selected_public ///
    cdid_imputation_public cdid_augmented_public cdid_baseline_reference_public cdid_selected_reference_public ///
    cdid_custom_baseline cdid_custom_selected cdid_custom_selected_reference ///
    cdid_cross_pair cdid_cross_pair_selected cdid_within_general cdid_cross_pair_reconstruction cdid_cross_pair_reconstruction_dated cdid_within_reconstruction cdid_cross_pair_nuisance cdid_chain cdid_etwfe cdid_etwfe_observation cdid_cross_pair_public cdid_public_fit cdid_public_inference cdid_static cdid_frontend_display cdid_ordinary cdid_ordinary_inference ///
    cdid_ordinary_aggregation cdid_ordinary_stream cdid_ordinary_state cdid_ordinary_post ///
    cdid_ordinary_intake cdid_ordinary_frontend cdid_chained_aggregate_public {
    capture confirm file "`candidate'/mata/`m'.mata"
    if !_rc quietly do "`candidate'/mata/`m'.mata"
    else quietly do "`root'/mata/`m'.mata"
}
capture confirm file "`candidate'/dependencies/csdid/csdid.mata"
if !_rc quietly do "`candidate'/dependencies/csdid/csdid.mata"
else quietly do "`root'/dependencies/csdid/csdid.mata"
foreach m in cdid_conditional cdid_ordinary_conditional cdid_custom_ordinary cdid_ordinary_pretest cdid_cpt_baseline ///
    cdid_cpt_selected cdid_cpt_imputation cdid_pretest_result cdid_pretest_grid ///
    cdid_pretest_imputation_owner cdid_pretest_chained cdid_custom_public cdid_custom_ordinary_att ///
    cdid_user_results cdid_chain_report cdid_public_load {
    capture confirm file "`candidate'/mata/`m'.mata"
    if !_rc quietly do "`candidate'/mata/`m'.mata"
    else quietly do "`root'/mata/`m'.mata"
}
mata: mata mlib create lcdid_native, dir("`library_out'") size(2048)
mata: mata mlib add lcdid_native cdid*() _cdid*() csdid*(), dir("`library_out'") complete
display "CDID_COMPOSED_NATIVE_BUILD_PASS"
