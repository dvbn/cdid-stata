*! Private native result payload/registry and binary persistence, current StataNow19.5.
*! Native I/O: installed [M-5] fopen(), fputmatrix()/fgetmatrix(strict=1).
*! Stores caller data only. No estimator, inference, mass or role derivation.
version 19.5
mata:
struct cdid_state_limits {
    real scalar payload_slots, payload_string_bytes
    real scalar registry_slots, registry_string_bytes
}
struct cdid_result_state {
    string scalar schema, id_kind, cluster_kind
    real matrix ledger, IF, V_estimator, V_asymptotic
    string rowvector ledger_columns
    real colvector point, available, role
    string colvector reason
    real colvector unit_numeric_id, cohort, cohort_known, weight
    string colvector unit_string_id
    real matrix unit_metadata
    string rowvector unit_columns
    real colvector cluster_numeric_id
    string colvector cluster_string_id
    real matrix shares, share_IF
    string rowvector share_columns, share_IF_columns
    real colvector calendar, numeric_settings
    string colvector numeric_setting_names
    string matrix text_settings, source_versions
}
struct cdid_state_registry {
    real scalar next
    string scalar stamp
    string colvector tokens
    pointer() colvector payloads
    real matrix costs
}
void _cdid_rs_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
// Empty builder: caller must set ID/cluster kinds and source provenance.
struct cdid_result_state scalar cdid_state_empty()
{
    struct cdid_result_state scalar P
    P.schema="cdid-native-result-state/v1"; P.id_kind=P.cluster_kind=""
    P.ledger=P.IF=P.V_estimator=P.V_asymptotic=J(0,0,.)
    P.point=P.available=P.role=P.unit_numeric_id=P.cohort=P.cohort_known=P.weight=J(0,1,.)
    P.cluster_numeric_id=P.calendar=P.numeric_settings=J(0,1,.)
    P.reason=P.unit_string_id=P.cluster_string_id=P.numeric_setting_names=J(0,1,"")
    P.unit_metadata=P.shares=P.share_IF=J(0,0,.)
    P.ledger_columns=P.unit_columns=P.share_columns=P.share_IF_columns=J(1,0,"")
    P.text_settings=P.source_versions=J(0,2,"")
    return(P)
}
void _cdid_rs_limits(struct cdid_state_limits scalar L)
{
    real rowvector v
    v=(L.payload_slots,L.payload_string_bytes,L.registry_slots,L.registry_string_bytes)
    if (missing(v) | any(v:<0) | any(v:!=floor(v)) | any(v:>2^53-1))
        _cdid_rs_fail("Limits must be explicit finite nonnegative exact integers.",198)
}
real rowvector cdid_state_cost(struct cdid_result_state scalar P)
{
    real scalar n, b
    n=length(P.ledger)+length(P.IF)+length(P.V_estimator)+length(P.V_asymptotic)+
      length(P.point)+length(P.available)+length(P.role)+length(P.unit_numeric_id)+
      length(P.cohort)+length(P.cohort_known)+length(P.weight)+length(P.unit_metadata)+
      length(P.cluster_numeric_id)+length(P.shares)+length(P.share_IF)+length(P.calendar)+length(P.numeric_settings)
    b=strlen(P.schema)+strlen(P.id_kind)+strlen(P.cluster_kind)+sum(strlen(P.ledger_columns))+
      sum(strlen(P.reason))+sum(strlen(P.unit_string_id))+sum(strlen(P.unit_columns))+
      sum(strlen(P.cluster_string_id))+sum(strlen(P.share_columns))+sum(strlen(P.share_IF_columns))+
      sum(strlen(P.numeric_setting_names))+sum(strlen(P.text_settings))+sum(strlen(P.source_versions))
    return((n,b))
}
void _cdid_rs_budget(real rowvector cost, struct cdid_state_limits scalar L)
{
    _cdid_rs_limits(L)
    if (cols(cost)!=2 | rows(cost)!=1 | missing(cost) | any(cost:<0) | any(cost:!=floor(cost)))
        _cdid_rs_fail("Invalid stored payload size header.",498)
    if (cost[1]>L.payload_slots | cost[2]>L.payload_string_bytes)
        _cdid_rs_fail("Payload exceeds caller-supplied storage limits.",908)
}
void cdid_state_validate(struct cdid_result_state scalar P, struct cdid_state_limits scalar L)
{
    real scalar N, K
    _cdid_rs_budget(cdid_state_cost(P),L)
    N=rows(P.cohort); K=rows(P.point)
    if (P.schema!="cdid-native-result-state/v1") _cdid_rs_fail("Unsupported result-state schema.",498)
    if (!any(P.id_kind:==("numeric","string")) | !any(P.cluster_kind:==("none","numeric","string")))
        _cdid_rs_fail("Caller must specify numeric/string IDs and none/numeric/string clusters.",198)
    if (rows(P.ledger)!=K | cols(P.ledger)!=cols(P.ledger_columns) |
        rows(P.available)!=K | rows(P.role)!=K | rows(P.reason)!=K)
        _cdid_rs_fail("Table of cells, point, role, availability and reason rows must align.",498)
    if (rows(P.IF)!=0 & rows(P.IF)!=N | cols(P.IF)!=K)
        _cdid_rs_fail("IF must be full-unit by cell, or explicitly absent 0-by-cell.",498)
    if (!((rows(P.V_estimator)==K & cols(P.V_estimator)==K) | (rows(P.V_estimator)==0 & cols(P.V_estimator)==0)) |
        !((rows(P.V_asymptotic)==K & cols(P.V_asymptotic)==K) | (rows(P.V_asymptotic)==0 & cols(P.V_asymptotic)==0)))
        _cdid_rs_fail("Passed covariances must be cell-square or absent 0-by-0.",498)
    if (rows(P.unit_numeric_id)!=N | rows(P.unit_string_id)!=N | rows(P.cohort_known)!=N |
        rows(P.weight)!=N | rows(P.unit_metadata)!=N | rows(P.cluster_numeric_id)!=N |
        rows(P.cluster_string_id)!=N | cols(P.unit_metadata)!=cols(P.unit_columns))
        _cdid_rs_fail("Full-unit identities and metadata must align.",498)
    if (missing(P.cohort_known) | any((P.cohort_known:!=0):&(P.cohort_known:!=1)))
        _cdid_rs_fail("cohort_known must be supplied explicitly as 0/1.",498)
    if (cols(P.shares)!=cols(P.share_columns) | cols(P.share_IF)!=cols(P.share_IF_columns) |
        (rows(P.share_IF)!=0 & rows(P.share_IF)!=N))
        _cdid_rs_fail("Passed share metadata/IF labels and full-unit rows must align.",498)
    if (rows(P.numeric_settings)!=rows(P.numeric_setting_names) | cols(P.text_settings)!=2 |
        cols(P.source_versions)!=2 | rows(P.source_versions)==0 | any(P.source_versions:==""))
        _cdid_rs_fail("Settings and explicit source-version key/value pairs are required.",498)
}
struct cdid_state_registry scalar _cdid_rs_registry()
{
    pointer() scalar p
    struct cdid_state_registry scalar R
    p=findexternal("CDID_PRIVATE_RESULT_REGISTRY")
    if (p==NULL) {
        R.next=1; R.stamp=st_tempfilename(); R.tokens=J(0,1,""); R.payloads=J(0,1,NULL); R.costs=J(0,2,.)
        p=crexternal("CDID_PRIVATE_RESULT_REGISTRY"); *p=R
    } else {
        if (eltype(*p)!="struct") _cdid_rs_fail("Result registry name is already in use.",498)
        R=*p
        if (structname(R)!="cdid_state_registry") _cdid_rs_fail("Result registry name has another type.",498)
    }
    return(R)
}
real scalar _cdid_rs_index(struct cdid_state_registry scalar R, string scalar token)
{
    real colvector ix
    ix=vec(selectindex(R.tokens:==token))
    if (rows(ix)!=1) _cdid_rs_fail("Result token is not present in this session; load its saved file.",111)
    return(ix[1])
}
string scalar cdid_state_put(struct cdid_result_state scalar P, struct cdid_state_limits scalar L)
{
    struct cdid_state_registry scalar R
    struct cdid_result_state scalar snapshot
    real rowvector cost, retained
    string scalar token
    pointer() scalar p
    cdid_state_validate(P,L); cost=cdid_state_cost(P); R=_cdid_rs_registry()
    retained=rows(R.costs) ? colsum(R.costs) : (0,0)
    if (retained[1]+cost[1]>L.registry_slots | retained[2]+cost[2]>L.registry_string_bytes)
        _cdid_rs_fail("Retained results exceed caller-supplied registry limits.",908)
    if (R.next>=2^53) _cdid_rs_fail("Result token counter is exhausted.",498)
    token="cdidrs|"+R.stamp+"|"+strtrim(strofreal(R.next,"%21.0f")); R.next=R.next+1
    snapshot=P
    R.tokens=R.tokens\token; R.payloads=R.payloads\(&snapshot); R.costs=R.costs\cost
    p=findexternal("CDID_PRIVATE_RESULT_REGISTRY"); *p=R
    return(token)
}
struct cdid_result_state scalar cdid_state_get(string scalar token)
{
    struct cdid_state_registry scalar R
    struct cdid_result_state scalar P
    real scalar ix
    R=_cdid_rs_registry(); ix=_cdid_rs_index(R,token); P=*R.payloads[ix]
    return(P)
}
// Reading a row/cell chunk does not construct a complete result copy or e(matrix).
real matrix cdid_state_if(string scalar token, real colvector units, real colvector cells)
{
    struct cdid_state_registry scalar R
    pointer(struct cdid_result_state scalar) scalar p
    real scalar ix
    R=_cdid_rs_registry(); ix=_cdid_rs_index(R,token); p=R.payloads[ix]
    if (missing(units) | missing(cells) | any(units:<1) | any(units:>rows(p->IF)) |
        any(cells:<1) | any(cells:>cols(p->IF)) | any(units:!=floor(units)) | any(cells:!=floor(cells)))
        _cdid_rs_fail("IF chunk indices are out of range.",3301)
    if (!rows(units) | !rows(cells)) return(J(rows(units),rows(cells),.))
    return(p->IF[units,cells])
}
void cdid_state_drop(string scalar token)
{
    struct cdid_state_registry scalar R
    real scalar ix
    real colvector keep
    pointer() scalar p
    R=_cdid_rs_registry(); ix=_cdid_rs_index(R,token)
    keep=vec(selectindex((1::rows(R.tokens)):!=ix))
    if (rows(keep)) {
        R.tokens=R.tokens[keep,1]
        R.payloads=R.payloads[keep,1]
        R.costs=R.costs[keep,.]
    }
    else {
        R.tokens=J(0,1,"")
        R.payloads=J(0,1,NULL)
        R.costs=J(0,2,.)
    }
    p=findexternal("CDID_PRIVATE_RESULT_REGISTRY"); *p=R
}
void _cdid_rs_write(real scalar fh, transmorphic matrix X)
{
    real scalar rc
    rc=_fputmatrix(fh,X)
    if (rc<0) {
        (void) _fclose(fh)
        _cdid_rs_fail("Native result-state write failed; partial file retained.",-rc)
    }
}
transmorphic matrix _cdid_rs_read(real scalar fh)
{
    transmorphic matrix X
    X=_fgetmatrix(fh,1)
    if (fstatus(fh)!=0) {
        (void) _fclose(fh)
        _cdid_rs_fail("Incomplete or unreadable native result-state file.",498)
    }
    return(X)
}
void cdid_state_save(string scalar token, string scalar path, struct cdid_state_limits scalar L)
{
    struct cdid_state_registry scalar R
    pointer(struct cdid_result_state scalar) scalar p
    real scalar fh, ix, rc
    R=_cdid_rs_registry(); ix=_cdid_rs_index(R,token); p=R.payloads[ix]
    cdid_state_validate(*p,L)
    fh=_fopen(path,"w") // documented no-overwrite mode, including existing empty files
    if (fh<0) _cdid_rs_fail("cannot create the file; saving() never replaces an existing file",-fh)
    _cdid_rs_write(fh,("cdid-native-result-state","1"))
    _cdid_rs_write(fh,cdid_state_cost(*p))
    _cdid_rs_write(fh,*p)
    rc=_fclose(fh)
    if (rc<0) _cdid_rs_fail("Native result-state close failed.",-rc)
}
string scalar cdid_state_load(string scalar path, struct cdid_state_limits scalar L)
{
    real scalar fh, status, has_nuis, has_chain
    transmorphic matrix H, cost, raw, tail, nuis, chain
    struct cdid_result_state scalar P
    string scalar tok
    _cdid_rs_limits(L)
    fh=_fopen(path,"r")
    if (fh<0) _cdid_rs_fail("cannot open the file; check its name and location",-fh)
    H=_cdid_rs_read(fh)
    if (eltype(H)!="string" | rows(H)!=1 | cols(H)!=2) {
        (void) _fclose(fh); _cdid_rs_fail("Not a native result-state header.",498)
    }
    if (H[1]!="cdid-native-result-state" | H[2]!="1") {
        (void) _fclose(fh); _cdid_rs_fail("The file was not saved by cdid or estat, or by a version of cdid that cdid_load cannot read.",498)
    }
    cost=_cdid_rs_read(fh)
    if (eltype(cost)!="real" | rows(cost)!=1 | cols(cost)!=2) {
        (void) _fclose(fh); _cdid_rs_fail("Invalid result-state size header.",498)
    }
    // Check declared storage before reading the payload; local files are trusted.
    if (missing(cost) | any(cost:<0) | any(cost:!=floor(cost)) |
        cost[1]>L.payload_slots | cost[2]>L.payload_string_bytes) {
        (void) _fclose(fh); _cdid_rs_fail("Saved payload exceeds caller limits or has an invalid size header.",908)
    }
    raw=_cdid_rs_read(fh); tail=_fgetmatrix(fh,1); status=fstatus(fh)
    has_nuis=0; has_chain=0
    while (status==0 & eltype(tail)=="string" & rows(tail)==1 & cols(tail)==2) {
        if (tail[1]=="cdid-nuisance-reports" & tail[2]=="1" & !has_nuis) {
            nuis=_cdid_rs_read(fh); has_nuis=1
        }
        else if (tail[1]=="cdid-chain-info" & tail[2]=="1" & !has_chain) {
            chain=_cdid_rs_read(fh); has_chain=1
        }
        else break
        tail=_fgetmatrix(fh,1); status=fstatus(fh)
    }
    fclose(fh)
    if (status!=-1 | eltype(raw)!="struct") _cdid_rs_fail("Invalid native result-state payload/end marker.",498)
    P=raw
    if (structname(P)!="cdid_result_state") _cdid_rs_fail("Native result-state payload type differs.",498)
    cdid_state_validate(P,L)
    if (any(cdid_state_cost(P):!=cost)) _cdid_rs_fail("Stored result-state size does not match payload.",498)
    tok=cdid_state_put(P,L)
    if (has_nuis) cdid_nuis_loaded(tok,nuis)
    if (has_chain) cdid_chain_loaded(tok,chain)
    return(tok)
}
end
