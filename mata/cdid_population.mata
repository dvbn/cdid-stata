*! Population and exact-date input primitives for cdid; Stata 14.2 or later.
*! Native adaptation of David Benatia et al., cdid 1.1.0.0 (GPL-3),
*! :
*! R/chained_did_request.R and R/chained_reference_inputs.R.
version 14.2
mata:
struct cdid_population_ids {
    real scalar string_id, frame_supplied
    real colvector numeric_id, row_unit
    string colvector string_id_value
}
struct cdid_population {
    struct cdid_population_ids scalar ids
    real scalar N, R, n_panel, n_A, n_B, n_frame_only
    real colvector date, y, calendar, cohort, cohort_known
    real colvector first_row, row_count, observed_count, A, B, frame_only
    real colvector row_order, row_start, observed_order, A_to_F, F_to_A, B_to_F, F_to_B
    real matrix obs_grid
}
struct cdid_population_column {
    string scalar name, role, source, kind
    real scalar active, string_value, first_change
    real colvector data_value, unit_value, data_present, unit_present
    string colvector data_label, unit_label
}
struct cdid_population_values {
    real colvector value, present, row_index, missing_unit, missing_date
    string colvector label, support
}
struct cdid_population_gate {
    real scalar accepted, active
    string scalar error_class, message
}
struct cdid_population_weight {
    string scalar kind, fix_weights, reduction, source
    real colvector data_value, scientific_static, V, owner_weight
}
void _cdid_population_fail(string scalar cls, string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_population_index(real matrix x)
{
    return(vec(selectindex(x)))
}
real colvector _cdid_population_lookup(real colvector x, real colvector table)
{
    real colvector ans
    real scalar j, n, nt, lo, hi, mid, v, hit
    n=rows(x); nt=rows(table)
    ans=J(n,1,0)
    if (n==0 | nt==0) return(ans)
    if (n>1) {
        v=x[1]
        if (all(x:==v)) {
            lo=1; hi=nt
            while (lo<=hi) {
                mid=floor((lo+hi)/2)
                if (table[mid]==v) return(J(n,1,mid))
                if (table[mid]<v) lo=mid+1
                else hi=mid-1
            }
            return(ans)
        }
    }
    if (nt<=32 & n>=16) {
        for (mid=1;mid<=nt;mid++) ans=ans:+mid:*(x:==table[mid])
        return(ans)
    }
    hit=0; v=.
    for (j=1;j<=n;j++) {
        if (j>1) {
            if (x[j]==v) { ans[j]=hit; continue; }
            if (hit>0 & hit<nt) {
                if (x[j]==table[hit+1]) { v=x[j]; hit=hit+1; ans[j]=hit; continue; }
            }
        }
        v=x[j]; lo=1; hi=nt; hit=0
        while (lo<=hi) {
            mid=floor((lo+hi)/2)
            if (table[mid]==v) { hit=mid; break; }
            if (table[mid]<v) lo=mid+1
            else hi=mid-1
        }
        ans[j]=hit
    }
    return(ans)
}
void _cdid_mlib_ready() {}
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_panelsetup(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
/* Sorted distinct dates; integer-valued compact calendars avoid two full sorts. */
real colvector _cdid_population_calendar(real colvector date)
{
    real scalar n, lo, hi
    real colvector mark
    n=rows(date)
    if (n>0) {
        lo=min(date); hi=max(date)
        if (lo>=1 & hi<. & hi>lo & hi-lo<=4*n+64 & hi-lo<=1e7) {
            if (all(date:==floor(date))) {
                mark=J(hi-lo+1,1,0)
                mark[date:-(lo-1)]=J(n,1,1)
                return(select((lo::hi),mark))
            }
        }
    }
    return(uniqrows(sort(date,1)))
}
real colvector _cdid_cohort_index_vec(real colvector cohorts, real scalar nt_index,
    real colvector c_vals, real colvector known_mask)
{
    real colvector ans, fin_ix
    real scalar n, n_fin
    n=rows(c_vals)
    if (n==0) return(J(0,1,0))
    ans=J(n,1,0)
    n_fin=(nt_index>0 ? nt_index-1 : rows(cohorts))
    if (nt_index>0) {
        fin_ix=_cdid_population_index((known_mask:!=0):&(c_vals:>=.))
        if (rows(fin_ix)) ans[fin_ix]=J(rows(fin_ix),1,nt_index)
    }
    if (n_fin>0) {
        fin_ix=_cdid_population_index((known_mask:!=0):&(c_vals:<.))
        if (rows(fin_ix)) ans[fin_ix]=_cdid_population_lookup(c_vals[fin_ix],cohorts[|1\n_fin|])
    }
    return(ans)
}
/* Numeric identity is compared as doubles, never formatted into strings. Codes
   are original full-frame positions; the original values remain reversible. */
struct cdid_population_ids scalar cdid_population_ids_numeric(real colvector id,
    real colvector frame_id, real scalar frame_supplied)
{
    struct cdid_population_ids scalar z
    real colvector o, head, first, f, ix
    real scalar n, nf
    n=rows(id)
    if (n==0 | missing(id)) _cdid_population_fail("cdid_chained_request_input_error","IDs must be finite and nonmissing.",459)
    if (n==1) o=1
    else if (all(id[|2\n|]:>=id[|1\n-1|])) o=(1::n)
    else o=order((id,(1::n)),(1,2))
    head=1
    if (n>1) head=head\(1:+_cdid_population_index(id[o[|2\n|]]:!=id[o[|1\n-1|]]))
    first=sort(o[head],1)
    f=frame_supplied ? frame_id : id[first]
    nf=rows(f)
    if (nf==0 | missing(f)) _cdid_population_fail("cdid_chained_request_frame_error","Frame IDs must be finite, unique and nonmissing.",459)
    if (nf==1) o=1
    else if (all(f[|2\nf|]:>f[|1\nf-1|])) o=(1::nf)
    else o=order(f,1)
    if (nf>1) if (any(f[o[|2\nf|]]:==f[o[|1\nf-1|]]))
        _cdid_population_fail("cdid_chained_request_frame_error","Frame IDs must be unique.",459)
    ix=_cdid_population_lookup(id,f[o])
    if (any(ix:==0)) _cdid_population_fail("cdid_chained_request_frame_error","The frame must cover every data unit.",459)
    z.string_id=0; z.frame_supplied=frame_supplied
    z.numeric_id=f; z.string_id_value=J(0,1,""); z.row_unit=o[ix]
    return(z)
}
struct cdid_population_ids scalar cdid_population_ids_string(string colvector id,
    string colvector frame_id, real scalar frame_supplied)
{
    struct cdid_population_ids scalar z
    real colvector o, first, head, ix
    string colvector f
    real scalar n, j, lo, hi, mid
    n=rows(id)
    if (n==0 | any(id:=="")) _cdid_population_fail("cdid_chained_request_input_error","String IDs must be nonempty and nonmissing.",459)
    /* Sort by string, then take the minimum original row in each run. */
    o=order(id,1); first=J(0,1,0); j=1
    while (j<=n) {
        lo=j
        while (j<n) { if (id[o[j+1]]!=id[o[lo]]) break; j++; }
        first=first\min(o[|lo\j|]); j++
    }
    first=sort(first,1); f=frame_supplied ? frame_id : id[first]
    if (rows(f)==0 | any(f:=="")) _cdid_population_fail("cdid_chained_request_frame_error","Frame string IDs must be nonempty and unique.",459)
    o=order(f,1)
    if (rows(f)>1) if (any(f[o[|2\rows(f)|]]:==f[o[|1\rows(f)-1|]]))
        _cdid_population_fail("cdid_chained_request_frame_error","Frame IDs must be unique.",459)
    ix=J(n,1,0)
    for (j=1;j<=n;j++) {
        lo=1; hi=rows(f)
        while (lo<=hi) {
            mid=floor((lo+hi)/2)
            if (f[o[mid]]==id[j]) { ix[j]=o[mid]; break; }
            if (f[o[mid]]<id[j]) lo=mid+1
            else hi=mid-1
        }
    }
    if (any(ix:==0)) _cdid_population_fail("cdid_chained_request_frame_error","The frame must cover every data unit.",459)
    z.string_id=1; z.frame_supplied=frame_supplied
    z.numeric_id=J(0,1,0); z.string_id_value=f; z.row_unit=ix
    return(z)
}
/* All supplied rows are retained, including missing-Y rows. Unknown frame-only
   cohorts have known=0; known never-treated has cohort=. AND known=1. */
struct cdid_population scalar cdid_population_prepare(struct cdid_population_ids scalar ids,
    real colvector date, real colvector y, real colvector g,
    real colvector frame_g, real scalar frame_has_g)
{
    struct cdid_population scalar p
    real scalar i, u, n, nt
    real colvector pg, o, dcol, grid_vec, u_act, u_zero
    real matrix spans
    p.ids=ids; p.R=rows(ids.row_unit)
    p.N=ids.string_id ? rows(ids.string_id_value) : rows(ids.numeric_id)
    n=p.R
    if (rows(date)!=n | rows(y)!=n | rows(g)!=n | missing(date))
        _cdid_population_fail("cdid_chained_request_input_error","time() and gvar() must be nonmissing and finite, and the outcome numeric, in every observation used.",459)
    if (missing(g) | any(g:<0)) _cdid_population_fail("cdid_chained_request_cohort_error","gvar() must be nonmissing and nonnegative for every unit; code never-treated units as 0.",459)
    p.date=date; p.y=y; p.calendar=_cdid_population_calendar(date)
    p.first_row=J(p.N,1,0); p.row_count=J(p.N,1,0); p.observed_count=J(p.N,1,0)
    pg=J(p.N,1,.)
    if ((n==0 ? 0 : (n==1 ? 1 : all(ids.row_unit[|2\n|]:>=ids.row_unit[|1\n-1|])))) {
        spans=_cdid_panelsetup(ids.row_unit,1)
        u_act=ids.row_unit[spans[,1]]
        p.first_row[u_act]=spans[,1]
        pg[u_act]=g[spans[,1]]
        p.row_count[u_act]=spans[,2]:-spans[,1]:+1
        p.observed_count[u_act]=panelsum(y:<.,spans)
    }
    else {
        for (i=1;i<=n;i++) {
            u=ids.row_unit[i]
            if (p.first_row[u]==0) { p.first_row[u]=i; pg[u]=g[i]; }
            p.row_count[u]=p.row_count[u]+1
            if (y[i]<.) p.observed_count[u]=p.observed_count[u]+1
        }
    }
    if (n==0) o=J(0,1,0)
    else if (n==1) o=1
    else if (all((ids.row_unit[|2\n|]:>ids.row_unit[|1\n-1|]):|((ids.row_unit[|2\n|]:==ids.row_unit[|1\n-1|]):&(date[|2\n|]:>=date[|1\n-1|])))) o=(1::n)
    else o=order((ids.row_unit,date),(1,2))
    p.row_order=o
    if (n>1) if (any((ids.row_unit[o[|2\n|]]:==ids.row_unit[o[|1\n-1|]]):&(date[o[|2\n|]]:==date[o[|1\n-1|]])))
        _cdid_population_fail("cdid_chained_request_input_error","ivar() and time() do not uniquely identify the observations; each unit may appear at most once per period.",459)
    if (n>0) if (any(g:!=pg[ids.row_unit]))
        _cdid_population_fail("cdid_chained_request_cohort_error","gvar() varies within unit; each unit must have a single first-treatment date.",459)
    p.frame_only=p.row_count:==0; p.cohort_known=1:-p.frame_only
    if (ids.frame_supplied & frame_has_g) {
        if (rows(frame_g)!=p.N | missing(frame_g) | any(frame_g:<0))
            _cdid_population_fail("cdid_chained_request_cohort_error","Frame cohort column must be complete finite nonnegative numeric.",459)
        if (any((!p.frame_only):&(pg:!=frame_g)))
            _cdid_population_fail("cdid_chained_request_cohort_error","Frame cohorts disagree with data cohorts by id.",459)
        pg=frame_g; p.cohort_known=J(p.N,1,1)
    }
    p.cohort=pg
    u_zero=_cdid_population_index(p.cohort_known:&(pg:==0))
    if (rows(u_zero)) p.cohort[u_zero]=J(rows(u_zero),1,.)
    p.A=p.observed_count:>=2; p.B=p.observed_count:>=1
    p.A_to_F=_cdid_population_index(p.A); p.B_to_F=_cdid_population_index(p.B)
    p.F_to_A=J(p.N,1,0); p.F_to_B=J(p.N,1,0)
    if (rows(p.A_to_F)) p.F_to_A[p.A_to_F]=(1::rows(p.A_to_F))
    if (rows(p.B_to_F)) p.F_to_B[p.B_to_F]=(1::rows(p.B_to_F))
    p.n_panel=sum(p.row_count:>0); p.n_A=sum(p.A); p.n_B=sum(p.B); p.n_frame_only=sum(p.frame_only)
    p.row_start=runningsum(p.row_count):-p.row_count:+1
    p.observed_order=select(o,y[o]:<.)
    p.obs_grid=J(0,0,0)
    nt=rows(p.calendar)
    if (p.R>0 & nt>0 & p.N*nt<=8000000) {
        dcol=_cdid_population_lookup(date,p.calendar)
        grid_vec=J(p.N*nt,1,0)
        grid_vec[(ids.row_unit:-1):*nt:+dcol]=(1::p.R)
        p.obs_grid=colshape(grid_vec,nt)
    }
    return(p)
}
/* Query exact dates inside each unit's sparse supplied-row postings. Zero means
   no actual row, independently of whether a supplied row has missing Y. */
real colvector cdid_population_rows(struct cdid_population scalar p,
    real colvector units, real colvector dates)
{
    real colvector ans, dcols, ok, grid_vec
    real scalar j, k, nt, d0, c0, lo, hi, mid, row, u
    if (rows(units)!=rows(dates) | missing(units) | missing(dates) |
        any(units:<1) | any(units:>p.N) | any(units:!=floor(units)))
        _cdid_population_fail("cdid_reference_support_error","Invalid original-frame positions or dates.",459)
    k=rows(units)
    if (k==0) return(J(0,1,0))
    nt=rows(p.calendar)
    if (rows(p.obs_grid)==p.N & cols(p.obs_grid)==nt) {
        d0=dates[1]
        if (k==1 | all(dates:==d0)) {
            lo=1; hi=nt; c0=0
            while (lo<=hi) {
                mid=floor((lo+hi)/2)
                if (p.calendar[mid]==d0) { c0=mid; break; }
                if (p.calendar[mid]<d0) lo=mid+1
                else hi=mid-1
            }
            return(c0>0 ? p.obs_grid[units,c0] : J(k,1,0))
        }
        dcols=_cdid_population_lookup(dates,p.calendar)
        ans=J(k,1,0)
        ok=_cdid_population_index(dcols:>0)
        if (rows(ok)>0) {
            if (rows(ok)<=256) {
                for (j=1;j<=rows(ok);j++) {
                    u=ok[j]
                    ans[u]=p.obs_grid[units[u],dcols[u]]
                }
            }
            else {
                grid_vec=colshape(p.obs_grid,1)
                ans[ok]=grid_vec[(units[ok]:-1):*nt:+dcols[ok]]
            }
        }
        return(ans)
    }
    ans=J(k,1,0)
    for (j=1;j<=k;j++) {
        u=units[j]; lo=p.row_start[u]; hi=lo+p.row_count[u]-1
        while (lo<=hi) {
            mid=floor((lo+hi)/2); row=p.row_order[mid]
            if (p.date[row]==dates[j]) { ans[j]=row; break; }
            if (p.date[row]<dates[j]) lo=mid+1
            else hi=mid-1
        }
    }
    return(ans)
}
struct cdid_population_column scalar _cdid_population_column_empty(string scalar name,
    string scalar role, real scalar active)
{
    struct cdid_population_column scalar s
    s.name=name; s.role=role; s.active=active; s.string_value=0
    s.source="none"; s.kind=active ? "absent" : "inactive"; s.first_change=.
    s.data_value=s.unit_value=s.data_present=s.unit_present=J(0,1,0)
    s.data_label=s.unit_label=J(0,1,"")
    return(s)
}
/* Classify raw variables ONLY after the role's evaluated static gate fails.
   Data missingness is ignored for detecting change; frame supply is
   dateless and cannot override any disagreeing supplied data value. */
struct cdid_population_column scalar cdid_population_column_numeric(struct cdid_population scalar p,
    string scalar name, string scalar role, real colvector data_value, real colvector frame_value,
    real scalar in_data, real scalar in_frame, real scalar date_variable, real scalar active)
{
    struct cdid_population_column scalar s
    real scalar i, u
    real colvector u_act, bad, pres
    s=_cdid_population_column_empty(name,role,active)
    if (!active) return(s)
    if (date_variable) { s.source="date"; s.kind="date"; return(s); }
    if (in_data & rows(data_value)!=p.R) _cdid_population_fail("cdid_chained_request_covariate_error","Raw data column length mismatch.",459)
    if (p.ids.frame_supplied & in_frame & rows(frame_value)!=p.N) _cdid_population_fail("cdid_chained_request_covariate_error","Raw frame column length mismatch.",459)
    s.data_value=in_data ? data_value : J(p.R,1,.)
    s.data_present=s.data_value:<.
    s.unit_value=J(p.N,1,.); s.unit_present=J(p.N,1,0)
    if (p.ids.frame_supplied & in_frame) {
        s.source="frame"; s.kind="unit_attribute"; s.unit_value=frame_value; s.unit_present=frame_value:<.
        if (in_data & p.R>0) {
            if (all(s.data_present)) {
                bad=_cdid_population_index((!s.unit_present[p.ids.row_unit]):|(data_value:!=frame_value[p.ids.row_unit]))
                if (rows(bad)) { s.kind="disagree"; s.first_change=bad[1]; }
            }
            else {
                pres=_cdid_population_index(s.data_present)
                if (rows(pres)) {
                    bad=_cdid_population_index((!s.unit_present[p.ids.row_unit[pres]]):|(data_value[pres]:!=frame_value[p.ids.row_unit[pres]]))
                    if (rows(bad)) { s.kind="disagree"; s.first_change=pres[bad[1]]; }
                }
            }
        }
        return(s)
    }
    if (!in_data) return(s)
    s.source="data"; s.kind="unit_attribute"
    if (p.R>0 & all(s.data_present)) {
        u_act=_cdid_population_index(!p.frame_only)
        if (rows(u_act)) {
            s.unit_value[u_act]=data_value[p.first_row[u_act]]
            s.unit_present[u_act]=J(rows(u_act),1,1)
        }
        bad=_cdid_population_index(data_value:!=s.unit_value[p.ids.row_unit])
        if (rows(bad)) { s.kind="dated"; s.first_change=bad[1]; }
        return(s)
    }
    for (i=1;i<=p.R;i++) if (s.data_present[i]) {
        u=p.ids.row_unit[i]
        if (!s.unit_present[u]) { s.unit_value[u]=data_value[i]; s.unit_present[u]=1; }
        else if (data_value[i]!=s.unit_value[u] & s.kind!="dated") { s.kind="dated"; s.first_change=i; }
    }
    return(s)
}
/* Presence is explicit for string raw values, so a supplied empty factor label
   can be distinguished from an absent value by the formula/encoding adapter. */
struct cdid_population_column scalar cdid_population_column_string(struct cdid_population scalar p,
    string scalar name, string scalar role, string colvector data_value, string colvector frame_value,
    real colvector data_present, real colvector frame_present,
    real scalar in_data, real scalar in_frame, real scalar active)
{
    struct cdid_population_column scalar s
    real scalar i, u
    s=_cdid_population_column_empty(name,role,active); s.string_value=1
    if (!active) return(s)
    if (in_data & (rows(data_value)!=p.R | rows(data_present)!=p.R)) _cdid_population_fail("cdid_chained_request_covariate_error","Raw string data column length mismatch.",459)
    if (p.ids.frame_supplied & in_frame & (rows(frame_value)!=p.N | rows(frame_present)!=p.N)) _cdid_population_fail("cdid_chained_request_covariate_error","Raw string frame column length mismatch.",459)
    s.data_label=in_data ? data_value : J(p.R,1,""); s.data_present=in_data ? data_present : J(p.R,1,0)
    s.unit_label=J(p.N,1,""); s.unit_present=J(p.N,1,0)
    if (missing(s.data_present) | any((s.data_present:!=0):&(s.data_present:!=1))) _cdid_population_fail("cdid_chained_request_covariate_error","Presence masks must be binary.",459)
    if (p.ids.frame_supplied & in_frame) {
        if (missing(frame_present) | any((frame_present:!=0):&(frame_present:!=1))) _cdid_population_fail("cdid_chained_request_covariate_error","Presence masks must be binary.",459)
        s.source="frame"; s.kind="unit_attribute"; s.unit_label=frame_value; s.unit_present=frame_present
        if (in_data) for (i=1;i<=p.R;i++) if (s.data_present[i]) {
            u=p.ids.row_unit[i]
            if (!s.unit_present[u] | data_value[i]!=frame_value[u]) { s.kind="disagree"; s.first_change=i; break; }
        }
        return(s)
    }
    if (!in_data) return(s)
    s.source="data"; s.kind="unit_attribute"
    for (i=1;i<=p.R;i++) if (s.data_present[i]) {
        u=p.ids.row_unit[i]
        if (!s.unit_present[u]) { s.unit_label[u]=data_value[i]; s.unit_present[u]=1; }
        else if (data_value[i]!=s.unit_label[u] & s.kind!="dated") { s.kind="dated"; s.first_change=i; }
    }
    return(s)
}
/* Evaluated static model matrices are supplied by the formula adapter. Required
   units are A, U or F according to the caller's route; no raw-variable gate is
   added to a static transformed design that the reference accepts. */
struct cdid_population_gate scalar cdid_population_static_gate(struct cdid_population scalar p,
    real matrix frame_design, real matrix data_design, real colvector required_units,
    real scalar compare_data, real scalar split_variables, real scalar allow_zero_col, real scalar active)
{
    struct cdid_population_gate scalar z
    real colvector rr, dr
    z.accepted=1; z.active=active; z.error_class="cdid_chained_request_covariate_error"; z.message=""
    if (!active) return(z)
    if (rows(required_units)!=p.N | missing(required_units) | any((required_units:!=0):&(required_units:!=1)))
        _cdid_population_fail("cdid_reference_support_error","Static fitting roster must be an explicit full-frame binary mask.",459)
    rr=_cdid_population_index(required_units)
    if (rows(frame_design)!=p.N | (!allow_zero_col & cols(frame_design)==0)) {
        z.accepted=0; z.message="Covariates must be supplied on the fitting roster."; return(z)
    }
    if (missing(frame_design[rr,])) { z.accepted=0; z.message="Covariates contain missing or nonfinite values for units used in estimation; drop those units or impute the covariates before calling cdid."; return(z); }
    if (split_variables) { z.accepted=0; z.message="Covariates are split between data and frame."; return(z); }
    if (compare_data) {
        dr=_cdid_population_index(required_units[p.ids.row_unit])
        if (rows(data_design)!=p.R | cols(data_design)!=cols(frame_design)) {
            z.accepted=0; z.message="Panel and frame model-matrix columns disagree."; return(z)
        }
        if (missing(data_design[dr,]) | any(data_design[dr,]:!=frame_design[p.ids.row_unit[dr],])) {
            z.accepted=0; z.message="Covariates vary within unit, have missing values, or disagree with frame by id."; return(z)
        }
    }
    return(z)
}
real scalar cdid_population_reference_gate(struct cdid_population_gate scalar gate,
    string rowvector kinds)
{
    if (!gate.active | gate.accepted) return(0)
    if (!any((kinds:=="dated") :| (kinds:=="date"))) _cdid_population_fail(gate.error_class,gate.message,459)
    if (any(kinds:=="disagree")) _cdid_population_fail(gate.error_class,"Raw data disagree with frame; frame supply cannot override failed agreement.",459)
    if (any(kinds:=="absent")) _cdid_population_fail(gate.error_class,"A required formula variable is absent from data and frame.",459)
    return(1)
}
/* Caller supplies the catalogue's exact required (original-frame row, calendar
   date, support label) rows. No support universe or extraction date is guessed.
   strict=0 returns complete missing diagnostics; strict=1 raises the typed
   missing-reference error. There is no carry-forward or observation deletion. */
struct cdid_population_values scalar cdid_population_read(struct cdid_population scalar p,
    struct cdid_population_column scalar s, real colvector units, real colvector dates,
    string colvector support, real scalar strict)
{
    struct cdid_population_values scalar z
    real colvector ix, bad, ok, r_ok
    real scalar j, k, u
    k=s.active ? rows(units) : 0
    z.value=J(k,1,.); z.label=J(k,1,""); z.present=J(k,1,0); z.row_index=J(k,1,0)
    z.missing_unit=z.missing_date=J(0,1,0); z.support=J(0,1,"")
    if (!s.active) return(z)
    if (rows(support)!=k) _cdid_population_fail("cdid_reference_support_error","One support label is required per requested row.",459)
    if (s.kind=="unit_attribute") {
        if (k>0 & rows(p.calendar)==0) _cdid_population_fail("cdid_reference_support_error","Invalid original-frame positions or dates.",459)
        if (missing(units) | any(units:<1) | any(units:>p.N) | any(units:!=floor(units)))
            _cdid_population_fail("cdid_reference_support_error","Invalid original-frame positions or dates.",459)
        ix=J(k,1,0)
        if (rows(dates)!=k) _cdid_population_fail("cdid_reference_support_error","One date or static anchor marker is required per requested row.",459)
    } else ix=cdid_population_rows(p,units,dates)
    z.row_index=ix
    if (s.kind=="absent" | s.kind=="disagree") _cdid_population_fail("cdid_chained_request_covariate_error","Absent or conflicting raw variable cannot be read.",459)
    if (s.kind=="date") {
        z.value=dates; z.present=_cdid_population_lookup(dates,p.calendar):>0
    } else if (k>0) {
        if (s.kind=="unit_attribute") {
            z.present=s.unit_present[units]
            if (s.string_value) z.label=s.unit_label[units]
            else z.value=s.unit_value[units]
        } else if (s.kind=="dated") {
            ok=_cdid_population_index(ix:>0)
            if (rows(ok)) {
                r_ok=ix[ok]
                z.present[ok]=s.data_present[r_ok]
                if (s.string_value) z.label[ok]=s.data_label[r_ok]
                else z.value[ok]=s.data_value[r_ok]
            }
        }
    }
    bad=_cdid_population_index(!z.present)
    z.missing_unit=units[bad]; z.missing_date=dates[bad]; z.support=support[bad]
    if (strict & rows(bad)) {
        j=bad[1]; u=units[j]
        errprintf("Missing required %s variable `%s`, id ",s.role,s.name)
        if (p.ids.string_id) errprintf("%s",p.ids.string_id_value[u])
        else errprintf("%21.17g",p.ids.numeric_id[u])
        errprintf(", date %21.17g, %s (no carry-forward, substitution or deletion).\n",dates[j],support[j])
        _cdid_population_fail("cdid_reference_missing_error","Required raw reference value is missing; full diagnostics are available with strict=0.",459)
    }
    return(z)
}
/* Static cluster labels may be encoded reversibly by the same numeric/string
   label dictionary as IDs. Validate codes on all F, including frame-only units.
   The formula/ado adapter owns the external labels; no cluster affects points. */
real colvector cdid_population_clusters(struct cdid_population scalar p,
    real colvector data_code, real colvector frame_code, real scalar in_data, real scalar in_frame)
{
    real colvector out, u_act
    out=J(p.N,1,.)
    if (!in_data & !in_frame) _cdid_population_fail("cdid_chained_request_cluster_error","Cluster variable not found in data or frame.",459)
    if (in_data) {
        if (rows(data_code)!=p.R | missing(data_code)) _cdid_population_fail("cdid_chained_request_cluster_error","Cluster labels must be plain nonmissing values.",459)
        if (p.R>0) {
            u_act=_cdid_population_index(!p.frame_only)
            if (rows(u_act)) out[u_act]=data_code[p.first_row[u_act]]
            if (any(data_code:!=out[p.ids.row_unit]))
                _cdid_population_fail("cdid_chained_request_cluster_error","cluster() varies within unit; each unit must belong to a single cluster.",459)
        }
    }
    if (in_frame) {
        if (rows(frame_code)!=p.N | missing(frame_code)) _cdid_population_fail("cdid_chained_request_cluster_error","Frame-only units need nonmissing cluster labels.",459)
        if (in_data & any((!p.frame_only):&(out:!=frame_code)))
            _cdid_population_fail("cdid_chained_request_cluster_error","Frame cluster labels disagree with data by id.",459)
        out=frame_code
    }
    if (missing(out)) _cdid_population_fail("cdid_chained_request_cluster_error","Frame-only units need cluster labels; zero-eligible clusters are retained.",459)
    return(out)
}
struct cdid_population_weight scalar cdid_population_weights(struct cdid_population scalar p,
    real colvector wd, real colvector wf, real scalar in_data, real scalar in_frame,
    string scalar fix_weights, real scalar dated_design)
{
    struct cdid_population_weight scalar w
    struct cdid_population_column scalar s
    real colvector supplied, q, u_act
    real scalar i, u, varying
    if (fix_weights!="" & fix_weights!="first_period" & fix_weights!="base_period")
        _cdid_population_fail("cdid_chained_request_unsupported_weights","fixweights() must be first_period or base_period in chained mode.",459)
    w.kind="static"; w.fix_weights=fix_weights; w.source="none"; w.reduction="static"
    w.scientific_static=w.V=w.owner_weight=J(p.N,1,1); w.data_value=J(0,1,0)
    if (!in_data & !in_frame) {
        if (dated_design) w.reduction="reference_recipes"
        return(w)
    }
    s=cdid_population_column_numeric(p,"weight","sampling weights at the context anchor",wd,wf,in_data,in_frame,0,1)
    /* Frame provenance cannot mask data changes. */
    varying=0
    if (in_data) {
        if (rows(wd)!=p.R) _cdid_population_fail("cdid_chained_request_weight_error","Data weight length mismatch.",459)
        q=J(p.N,1,.)
        if (p.R>0 & !missing(wd)) {
            u_act=_cdid_population_index(!p.frame_only)
            if (rows(u_act)) q[u_act]=wd[p.first_row[u_act]]
            varying=any(wd:!=q[p.ids.row_unit])
        }
        else {
            for (i=1;i<=p.R;i++) if (wd[i]<.) {
                u=p.ids.row_unit[i]
                if (q[u]>=.) q[u]=wd[i]
                else if (wd[i]!=q[u]) varying=1
            }
        }
    }
    if (varying) {
        supplied=select(wd,wd:<.)
        if (any(supplied:<=0)) _cdid_population_fail("cdid_chained_request_weight_error","Sampling weights must be finite and positive.",459)
        if (p.ids.frame_supplied & in_frame) _cdid_population_fail("cdid_chained_request_weight_error","Frame weights cannot override changing data weights.",459)
        if (fix_weights=="") _cdid_population_fail("cdid_chained_request_unsupported_weights","Sampling weights vary within unit; specify fixweights(first_period) or fixweights(base_period).",459)
        w.kind="varying"; w.source="data"; w.data_value=wd; w.V=J(p.N,1,.)
        for (i=1;i<=p.R;i++) if (p.date[i]==p.calendar[1] & wd[i]<.) w.V[p.ids.row_unit[i]]=wd[i]
        w.scientific_static=J(0,1,0)
    } else {
        q=J(p.N,1,1)
        if (in_data) {
            if (missing(wd) | any(wd:<=0)) _cdid_population_fail("cdid_chained_request_weight_error","Sampling weights must be finite and positive, including supplied missing-Y rows.",459)
            u_act=_cdid_population_index(!p.frame_only)
            if (rows(u_act)) q[u_act]=wd[p.first_row[u_act]]
        }
        if (in_frame) {
            if (rows(wf)!=p.N | missing(wf) | any(wf:<=0)) _cdid_population_fail("cdid_chained_request_weight_error","Frame weights must be finite and positive.",459)
            if (in_data & any((!p.frame_only):&(wf:!=q)))
                _cdid_population_fail("cdid_chained_request_weight_error","Frame weights disagree with data weights by id.",459)
            q=wf; w.source="frame"
        } else {
            if (p.n_frame_only>0) _cdid_population_fail("cdid_chained_request_weight_error","Frame-only units need weights; fixed weights are never imputed.",459)
            w.source="data"
        }
        w.scientific_static=w.V=q
    }
    if (!dated_design & w.kind=="static") w.owner_weight=w.scientific_static
    else if (!dated_design & fix_weights=="first_period" & !missing(w.V)) {
        w.owner_weight=w.V; w.reduction="first_period_static"
    } else { w.owner_weight=J(p.N,1,1); w.reduction="reference_recipes"; }
    return(w)
}
real scalar cdid_population_weight_anchor(struct cdid_population scalar p,
    struct cdid_population_weight scalar w, real scalar reference_date)
{
    if (w.kind=="static") return(.)
    if (w.fix_weights=="first_period") return(p.calendar[1])
    if (reference_date>=. | !any(p.calendar:==reference_date))
        _cdid_population_fail("cdid_reference_support_error","Context reference anchor is not a stored calendar date.",459)
    return(reference_date)
}
struct cdid_population_values scalar cdid_population_weight_read(struct cdid_population scalar p,
    struct cdid_population_weight scalar w, real colvector units, real colvector anchors,
    string colvector support, real scalar strict)
{
    struct cdid_population_column scalar s
    s=_cdid_population_column_empty("weight","sampling weights at the context anchor",1)
    s.kind=w.kind=="varying" ? "dated" : "unit_attribute"; s.source=w.source
    s.data_value=w.data_value; s.data_present=w.data_value:<.
    s.unit_value=w.scientific_static; s.unit_present=w.scientific_static:<.
    return(cdid_population_read(p,s,units,anchors,support,strict))
}
end
