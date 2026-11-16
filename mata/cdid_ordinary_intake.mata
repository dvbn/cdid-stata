*! Additional ordinary metadata translation. Estimator/extraction remain unchanged.
*! did2.5.1 faster-mode standardization normalizes BEFORE support/balance settling.
version 19.5
mata:
void cdid_ord_extraction_metadata(string scalar units_name,string scalar sample_name,
    string scalar time_name,real scalar normalization_mean,string scalar out_units,
    string scalar out_calendar)
{
    real matrix units
    real colvector sample,calendar
    units=st_matrix(units_name);sample=st_matrix(sample_name)
    if(cols(units)<3 | normalization_mean<=0 | normalization_mean>=.)
        _cdid_ord_fail("Raw ordinary unit weights and positive source normalization mean are required.",498)
    calendar=uniqrows(sort(st_data(sample,time_name),1))
    st_matrix(out_units,(units[,1..3],units[,3]/normalization_mean))
    st_matrix(out_calendar,calendar)
}
end
