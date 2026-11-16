# Investigate missingness and data structure

Before choosing an assumption about missing outcomes or a model to adjust for them, inspect how many outcomes each unit has, the coverage of each cohort and date, the overlap of pairs of dates, and the composition of the covariates. No diagnostic proves that mean changes are representative or that the covariates are sufficient.

## 1. Count observed outcomes, not rows

The examples use the data of the [hands-on do-file](../hands-on.md); run its Sections 0 and 1 first. These data have a row for every unit and date, including missing outcomes. Counting `_N` per unit would therefore count every unit as fully observed. Section 2 of the do-file counts nonmissing `y`:

```stata
* 2a. Count observed outcome dates per unit across all units in the long dataset
tempvar n_obs_unit tag_unit
bysort id: egen long `n_obs_unit' = count(y)
bysort id: generate byte `tag_unit' = (_n == 1)

display as text _n "=== Unit-level count of observed outcome dates (0, 1, >=2) ==="
tabulate `n_obs_unit' first_treat if `tag_unit', missing

* 2b. Inspect outcome observation rates by calendar period and cohort
display as text _n "=== Outcome observation indicator (obs_y) by period and cohort ==="
tabulate period first_treat, summarize(obs_y) means

* 2c. Compare covariate composition across all units vs repeat-observer units (>=2 dates)
display as text _n "=== Covariate means: all population units vs units with >=2 observed dates ==="
summarize x_bin x_cont z_sel if `tag_unit'
summarize x_bin x_cont z_sel if `tag_unit' & `n_obs_unit' >= 2
```

These counts distinguish the study population $F$, units observed at least once ($B$), units observed at least twice ($A$), singletons $B\setminus A$, and never-observed units $F\setminus B$. If units with no observed outcome have no rows in the long data, merge the counts with the unit table and set their count to zero. Do not count omitted or missing outcomes as observed dates.

## 2. Observation models and support

An observation logit with a free intercept has no finite estimate in a cell where every unit, or no unit, is observed. Pooling the slopes keeps the cell intercepts and does not solve this. Pure imputation has no observation logit, but has its own requirements: enough observed changes, full rank, and support for prediction. A complete cell alone therefore does not guarantee a successful imputation analysis.

Independent-period logits also require that observation at the two dates be independent given the cohort and the covariates. Rotating windows and persistent dropout generally violate this product-probability assumption. Pair logits model the joint indicator directly, but still need a correct covariate specification.

## 3. Inspect the results

After a fit, run `cdid_diagnostics`:

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) bstrap(0) cband(0)
cdid_diagnostics, detail
display "units: " r(N_units) ", at least one outcome: " r(N_with_outcome) ", at least two: " r(N_with_two_outcomes)
```

The diagnostics report the population counts, the reasons for unavailable cells, the models and GMM settings, and the inference. A quantity that was not stored with the results is reported as not retained, not as zero. Available unit covariates and observation counts can motivate substantive investigation, but neither logit significance nor a good response-model fit establishes the outcome mean-change assumption or valid transport. See [Missing outcomes](../methods/missing-outcomes.md) and [Sampling information](../methods/sampling-information.md).

---

[Guide home](../README.md) · [Package overview](../../README.md)
