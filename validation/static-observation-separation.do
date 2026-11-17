version 19.5
clear all
set more off
set seed 9822437
set obs 48
generate long id = _n
generate byte role = ceil(id/12)
generate byte z = mod(floor((id-1)/6), 2)
generate byte position = mod(id-1, 6) + 1
generate double x = -1.1 + (position-1)*2.3/5
generate double intercept = rnormal()
expand 2
bysort id: generate byte t = _n
generate byte d = cond(t==1, inlist(role,3,4), inlist(role,2,3))
generate double y = intercept + t*(.17+.13*x+.1*z) + d*(1.3+.2*x+.17*z) + rnormal(0,.5)
clonevar complete = y
local opts ivar(id) time(t) dvar(d) method(dr) selection(conditional_mtar) selcovars(z) bstrap(0) cband(0)

foreach model in pair_logit independent_period_logit {
    quietly replace y = complete
    quietly replace y = . if t==1 & position==6 & !(role==1 & z==1)
    quietly replace y = . if t==2 & position==5 & ///
        !(role==1 & z==1 & "`model'"=="pair_logit")
    generate byte shifted_z = z + 2
    foreach correction in ipw dr {
        capture noisily cdid y x z, `opts' selmethod(`correction') selmodel(`model')
        assert _rc==459
    }
    capture noisily cdid y x z, ivar(id) time(t) dvar(d) method(dr) ///
        selection(conditional_mtar) selcovars(shifted_z) selmethod(dr) ///
        selmodel(`model') bstrap(0) cband(0)
    assert _rc==459
    drop shifted_z

    quietly cdid y x z, `opts' selmethod(dr) selmodel(`model') selpooling(pooled)
    assert e(att)<.

    quietly replace y = . if t==1 & position==6 & role==1 & z==1
    quietly replace y = . if t==2 & position==5 & role==1 & z==1
    quietly cdid y x z, `opts' selmethod(dr) selmodel(`model')
    assert e(att)<.
}

quietly replace y = complete if role==1
capture noisily cdid y x z, ivar(id) time(t) dvar(d) method(dr) ///
    selection(conditional_mtar) selconstant selmethod(dr) ///
    selmodel(independent_period_logit) selpooling(pooled) bstrap(0) cband(0)
assert _rc==459
display "CDID_STATIC_OBSERVATION_SEPARATION_PASS"
