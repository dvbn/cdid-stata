*! Ordinary extraction-to-state composition. No runtime R/Python dependency.
version 19.5
mata:
string scalar cdid_ord_frontend_build(string scalar clustername,string scalar idname,
    real scalar rcs,real scalar compute_if,real scalar bstrap,real scalar cband,
    real scalar B,real scalar alpha,real scalar row_chunk,real scalar workspace,
    string scalar limits_name)
{
    struct cdid_result_state scalar P
    struct cdid_state_limits scalar L
    struct cdid_ord_summary scalar S
    struct cdid_ord_boot scalar M
    struct cdid_ord_inference scalar O
    struct cdid_infer_budget scalar budget
    real matrix cells,F,units,map
    real colvector sample,calendar,cluster,labels,distinct,o,ids,rawcluster
    string colvector scluster,slabels,alls
    real scalar n,k,C,j,cap,factor,d1,lo,hi,mid
    string matrix versions
    pointer() scalar ptr
    cells=st_matrix("r(cells)");F=st_matrix("r(inffunc)");units=st_matrix("r(source_units)")
    map=st_matrix("r(unit_index)");calendar=st_matrix("r(source_calendar)");sample=st_matrix("r(sample_rows)")
    factor=st_numscalar("r(source_weight_factor)");d1=st_numscalar("r(D1_exception)")
    n=rows(F);k=cols(F)
    if(n<1 | k<1 | B<1 | B!=floor(B) | B>=. | workspace<1 | workspace>=. | row_chunk<1 | row_chunk!=floor(row_chunk))
        _cdid_ord_fail("Invalid ordinary inference dimensions or explicit execution controls.",198)
    L.payload_slots=floor(workspace/8);L.payload_string_bytes=floor(workspace)
    L.registry_slots=floor(workspace/8);L.registry_string_bytes=floor(workspace)
    cap=floor(workspace/8)
    budget=cdid_infer_budget_create(n,n,k,B,min((n,row_chunk)),max((k,min((cap,row_chunk*k)))),
        min((B,row_chunk)),max((n,min((cap,row_chunk*n)))),cap,floor(workspace))
    cluster=labels=(1::n);C=n;scluster=J(n,1,"")
    if(clustername!="") {
        if(st_isstrvar(clustername)) {
            scluster=st_sdata(map[,2],clustername);slabels=uniqrows(sort(scluster,1));C=rows(slabels)
            for(j=1;j<=n;j++) cluster[j]=_cdid_ord_ix(slabels:==scluster[j])[1]
        }
        else {
            labels=st_data(map[,2],clustername)
            if(missing(labels)) _cdid_ord_fail("Retained ordinary cluster labels are missing.",459)
            distinct=uniqrows(sort(labels,1));C=rows(distinct)
            for(j=1;j<=n;j++) {
                lo=1;hi=C
                while(lo<=hi) {mid=floor((lo+hi)/2);if(distinct[mid]==labels[j]) break;if(distinct[mid]<labels[j]) lo=mid+1;else hi=mid-1;}
                cluster[j]=mid
            }
        }
        if(!rcs) {
            ids=st_data(sample,idname);o=order(ids,1);ids=ids[o]
            if(st_isstrvar(clustername)) {
                alls=st_sdata(sample[o],clustername)
                for(j=2;j<=rows(ids);j++) if(ids[j]==ids[j-1] & alls[j]!=alls[j-1]) _cdid_ord_fail("Time-varying ordinary cluster variables are not supported.",459)
            }
            else {
                rawcluster=st_data(sample[o],clustername)
                for(j=2;j<=rows(ids);j++) if(ids[j]==ids[j-1] & rawcluster[j]!=rawcluster[j-1]) _cdid_ord_fail("Time-varying ordinary cluster variables are not supported.",459)
            }
        }
    }
    versions=("R_cdid","1.3.0.0"\"did","2.5.1"\"csdid",st_global("r(backend_version)")\"rcs_dr_influence_function",(d1 ? "corrected" : "unaffected"))
    M=cdid_ord_boot_none(k)
    if(compute_if) {
        S=cdid_ord_if_summary(F,cluster,C,clustername!="",bstrap,row_chunk)
        if(bstrap) M=cdid_ord_boot_stream(S,B,alpha,0,budget,"caller",J(0,0,.))
    }
    O=cdid_ord_report_cells(cells[,6],cells[,1],cells[,2],compute_if ? F : J(0,k,.),cluster,C,
        clustername!="",clustername!="",compute_if,bstrap,cband,alpha,row_chunk,M)
    P=cdid_ord_state_build(cells,compute_if ? F : J(0,k,.),map,units[,2],units[,3],units[,4],factor,calendar,
        labels,cluster,C,clustername,d1,versions,O,alpha,B,L)
    if(clustername!="") if(st_isstrvar(clustername)) {P.cluster_kind="string";P.cluster_string_id=scluster;P.cluster_numeric_id=J(n,1,.);}
    _cdid_ord_meta_matrix(P,"ord.source_sample_rows",sample)
    _cdid_ord_meta_text(P,"ord.panel_mode",st_global("r(panel_mode)"));_cdid_ord_meta_text(P,"ord.compute_path",st_global("r(compute_path)"))
    _cdid_ord_meta_text(P,"ord.method",st_global("r(method)"));_cdid_ord_meta_text(P,"ord.control_group",st_global("r(control_group)"))
    _cdid_ord_meta_text(P,"ord.base_period",st_global("r(base_period)"));_cdid_ord_meta_text(P,"ord.fix_weights",st_global("r(fix_weights)"))
    _cdid_ord_meta_number(P,"ord.anticipation",st_numscalar("r(anticipation)"))
    _cdid_ord_meta_text(P,"ord.native_multiplier","accepted native caller stream; increasing draw batches, column-major cluster order inside batch; distinct from R RNG")
    ptr=crexternal(limits_name);*ptr=L
    return(cdid_state_put(P,L))
}
string scalar cdid_ord_consumer_aggregate(string scalar token,string scalar type,
    real scalar min_e,real scalar max_e,real scalar balance_e,real scalar na_rm,
    real scalar bstrap,real scalar cband,real scalar B,real scalar level,
    real scalar workspace,real scalar row_chunk,string scalar requested_cluster,string scalar limits_name)
{
    struct cdid_result_state scalar P
    struct cdid_ord_agg_request scalar R
    struct cdid_ord_aggregate_result scalar A
    struct cdid_infer_budget scalar budget
    struct cdid_state_limits scalar L
    real colvector cluster
    real scalar C,clustered,n,k,cap
    string scalar stored
    transmorphic scalar calls,fixed
    pointer() scalar ptr
    P=cdid_state_get(token);n=rows(P.cohort)
    R.type=type;R.min_e=min_e;R.max_e=max_e;R.balance_e=balance_e;R.na_rm=na_rm
    R.bstrap=(bstrap==-1 ? _cdid_ord_number(P,"ord.cells.bstrap") : bstrap)
    R.cband=(cband==-1 ? _cdid_ord_number(P,"ord.cells.cband") : cband)
    R.B=(B==. ? _cdid_ord_number(P,"ord.cells.B") : B)
    R.alpha=(level==. ? _cdid_ord_number(P,"ord.cells.alpha") : (level==95 ? .05 : 1-level/100))
    A=cdid_ord_aggregate_points(P.IF,P.point,P.ledger[,1],P.ledger[,2],P.cohort,P.weight,P.calendar,R)
    if(A.status!="computed") _cdid_ord_fail(A.reason,459)
    stored=_cdid_ord_text(P,"ord.cluster_name")
    if(requested_cluster=="inherit") requested_cluster=stored
    clustered=requested_cluster!="" & requested_cluster==stored
    cluster=(clustered ? P.unit_metadata[,_cdid_ord_column(P.unit_columns,"ordinary_cluster_position")] : (1::n))
    C=(clustered ? _cdid_ord_number(P,"ord.C") : n);k=max((1,cols(A.IF)))
    cap=floor(workspace/8)
    L.payload_slots=L.registry_slots=cap;L.payload_string_bytes=L.registry_string_bytes=floor(workspace)
    budget=cdid_infer_budget_create(n,C,k,R.B,min((n,row_chunk)),max((k,min((cap,row_chunk*k)))),
        min((R.B,row_chunk)),max((C,min((cap,row_chunk*C)))),cap,floor(workspace))
    fixed=asarray_create("string",1)
    calls=cdid_ord_aggregate_draws(A,R,cluster,C,clustered,row_chunk,budget,"caller",fixed)
    P=cdid_ord_state_aggregate(P,R,requested_cluster,row_chunk,calls,L)
    ptr=crexternal(limits_name);*ptr=L
    return(cdid_state_put(P,L))
}

/* Set varname to 1 on the listed data rows (used for e(sample)). */
void cdid_mark_rows(string scalar varname, real matrix rows)
{
    real colvector r
    r = vec(rows)
    if (rows(r) == 0) return
    if (missing(r) | any(r :< 1) | any(r :> st_nobs()) | any(r :!= floor(r)))
        _error(3300, "cdid_mark_rows(): row indices out of range")
    st_store(r, varname, J(rows(r), 1, 1))
}
end
