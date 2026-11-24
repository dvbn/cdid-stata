# Performance

This page reports computation times of `cdid` 1.0.0.0 on simulated incomplete panels, measured with the do-file shown below. The times describe one computer and one simulated design; they are not a benchmark against other software and say nothing about the statistical properties of the estimators. Run the same commands on your computer to obtain your own times.

## What was measured

The do-file simulates a staggered panel with three treated cohorts, never-treated units, five covariates and outcomes missing with a probability that depends on one covariate. The commands specify `gmmweighting(optimal)` to match the weighting used in version 1.0.0.0. The do-file times complete `cdid` calls, including data checks, estimation, influence functions, the multiplier bootstrap where requested, and the display and posting of results, for four corrections for missing outcomes and two panel sizes.

```stata
capture program drop perf_panel
program define perf_panel
    args n
    clear
    set seed 20261116
    quietly set obs `n'
    generate id = _n
    generate G = cond(runiform() < .3, 0, 2002 + floor(3 * runiform()))
    forvalues k = 1/5 {
        generate x`k' = rnormal()
    }
    quietly expand 6
    bysort id: generate year = 2000 + _n
    generate y = x1 + .5 * x2 + .2 * (year - 2000) + (G > 0 & year >= G) + rnormal()
    quietly replace y = . if runiform() < invlogit(-1.5 + .5 * x1)
end

* One untimed call first, so that one-time loading costs are not timed
perf_panel 500
quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal) bstrap(0) cband(0)

foreach n in 2000 8000 {
    perf_panel `n'
    timer clear
    timer on 1
    quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal)
    timer off 1
    timer on 2
    quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal) selection(conditional_mtar)
    timer off 2
    timer on 3
    quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal) selection(conditional_mtar) selmethod(imputation)
    timer off 3
    timer on 4
    quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal) selection(conditional_mtar) selmethod(dr)
    timer off 4
    timer on 5
    quietly cdid y x1-x5, ivar(id) time(year) gvar(G) gmmweighting(optimal) selection(conditional_mtar) selmethod(dr) bstrap(0) cband(0)
    timer off 5
    display as text _n "Units: `n'"
    timer list
}
display "processors: " c(processors) ", Stata " c(stata_version) ", " c(os) " " c(machine_type)
```

## Results

```text
Units: 2000
   1:      0.08 /        1 =       0.0840
   2:      0.27 /        1 =       0.2660
   3:      0.23 /        1 =       0.2300
   4:      0.18 /        1 =       0.1820
   5:      0.09 /        1 =       0.0910

Units: 8000
   1:      0.25 /        1 =       0.2510
   2:      0.30 /        1 =       0.3050
   3:      0.19 /        1 =       0.1910
   4:      0.31 /        1 =       0.3070
   5:      0.24 /        1 =       0.2380
```

```text
processors: 4, Stata 19.5, Unix Mac (Apple Silicon)
```

Timer 1 runs the version 1.0.0.0 default fit (doubly robust treatment adjustment, no correction for missing outcomes, 1000 bootstrap draws); timers 2 to 4 add observation weighting, change imputation and the augmented correction; timer 5 repeats timer 4 with analytical standard errors. Each line gives the elapsed seconds of one call. Times this short vary from run to run and indicate orders of magnitude only.

## Controlling memory and time

`workspace()` limits the memory used for the main numerical arrays (1 GiB by default), `chunk()` the number of rows processed at once, and `workcap()` the counted work; a fit that would exceed a limit is refused rather than truncated. `posting(state)` keeps results that are too large for Stata matrices. Analytical standard errors (`bstrap(0)`) avoid the bootstrap draws, and `noinffunc` skips the influence functions when only point estimates are needed.

---

[Guide home](README.md) · [Package overview](../README.md)
