# Working with cdid results

This page explains how results move between `cdid`, `estat` and the commands that use them. [Results](reference/results.md) lists every stored result.

## Replay, store and restore

`cdid` posts its results in `e()`; typing `cdid` alone replays them. `estat` replaces them with the aggregation, while the group-time results remain available to later `estat` commands. After `estat`, `cdid_chain` still reports the comparisons of the fit in the same session; `cdid_nuisance` needs the fit to be the current result, so store the fit before aggregating and restore it afterwards:

```stata
clear
set seed 20261116
set obs 600
generate id = _n
generate G = cond(runiform() < .4, 0, cond(runiform() < .5, 2004, 2006))
generate x = rnormal()
expand 8
bysort id: generate year = 2000 + _n
generate y = x + .2*(year - 2000) + (G > 0 & year >= G)*(1 + .25*(year - G)) + rnormal()
replace y = . if runiform() < .2
cdid y x, ivar(id) time(year) gvar(G) bstrap(0) cband(0)
estimates store fit
estat dynamic
cdid_chain, cohort(2004)
estimates restore fit
cdid
```

## Commands that use the stored results

None of the following commands uses the data in memory, refits a model or draws random numbers, and each leaves the data, `e()` and the random-number state as they were.

| Command | Purpose |
| --- | --- |
| `cdid_export using file, type(...)` | Save the results as a Stata dataset. |
| `cdid_diagnostics [, detail]` | Display the population, the cells, the settings, the inference and the tests; return them in `r()`. |
| `cdid_plot [, ...]` | Graph the estimates and their intervals. |
| `cdid_chain [, ...]` | Show the comparisons and GMM weights behind each effect. |
| `cdid_nuisance [component]` | Report the fitted models kept with `keepnuisance` in chained mode. After ordinary results, returns `r(status) = "unavailable"` and `r(reason)`. |

## Datasets written by cdid_export

- `type(cells)`, the default: every requested cell, including pretreatment, reference and unavailable cells, with the reason a cell is unavailable.
- `type(aggregate)`: the current aggregation of `estat` (bins, overall estimate, and the cells of each bin with their weights).
- `type(units)`: every unit of the study population, including units supplied only through `frame()`, with its identifier, cohort, weight and cluster.
- `type(influence)`: the influence functions in long form, one row per unit and estimate.
- `type(summary)`: every setting and result as key-value rows.

The variable `se` holds the reported standard error, from the multiplier bootstrap when it was used, and `analytic_se` the analytical one. `lower` and `upper` hold the reported pointwise interval or simultaneous band, and `interval_kind` says which. After an aggregation, `type(cells)` still describes the group-time cells. With `noinffunc`, standard errors and intervals are missing. After `dvar()` estimation, only `type(cells)` and `type(aggregate)` are available. The file name must have the extension `.dta` or none, and an existing file is overwritten only with `replace`.

```stata
cdid_export using cdid_cells, type(cells) replace
estat dynamic
cdid_export using cdid_dynamic, type(aggregate) replace
cdid_plot, export(cdid_dynamic.png, replace)
```

## Saved results

`saving(filename)` with `cdid`, `estat` or `cdid_pretest` writes the complete results, including the influence functions, the reports kept with `keepnuisance`, the record of units excluded under `missingreference(exclude)` and, for files saved by `cdid`, the information used by `cdid_chain`. File names created by `tempfile` are accepted. `cdid_load` restores the results of `cdid` and `estat` into `e()` without the data, and `cdid_pretest_load` restores the results of `cdid_pretest` into `r()`. These files differ from the datasets written by `cdid_export`, and `saving()` never overwrites an existing file. Results of `dvar()` estimation cannot be saved.

```stata
tempfile fit
cdid y x, ivar(id) time(year) gvar(G) bstrap(0) cband(0) saving(`fit')
clear
cdid_load using `fit'
estat group
```

After `cdid_load`, `e(cmdline)` holds the original `cdid` command, including when the file contains an aggregation. For older files without the original command line, it holds the `cdid_load` command. `e(state_file)` names the file. `e(sample)` is unset because the results no longer refer to data in memory.

## Large problems

`posting(matrix)`, the default, posts `e(b)` and `e(V)`. When the number of estimates exceeds Stata's matrix size, `posting(state)` keeps the complete results without posting them; `estat`, `cdid_export` and `cdid_plot` still work. `workspace()`, `chunk()` and `workcap()` limit memory and computation; see [Performance](performance.md).

[Guide home](README.md)
