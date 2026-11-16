# Estimate treatment effects from incomplete panels

Longitudinal datasets are often incomplete: units enter late, leave early, rotate in and out across waves, or miss individual interviews. In multi-period difference-in-differences, how you handle missing outcomes determines both which units contribute to estimation and which causal parameter is identified.

This article explains the four common panel observation structures, shows how chained difference-in-differences (Bellégo, Benatia and Dortet-Bernadet 2025) spans missing periods without dropping incomplete units or giving up within-unit differencing, and applies `cdid` to an incomplete panel.

## Four panel observation structures

In the diagrams below, each row is a unit and each column is a calendar period $t \in \{1, 2, 3, 4\}$. An `x` marks an observed outcome $Y_{i,t}$ and a `.` marks a missing outcome.

### 1. Balanced panel

Every unit is observed at every calendar period:

```text
Unit   t=1  t=2  t=3  t=4
  1     x    x    x    x
  2     x    x    x    x
  3     x    x    x    x
```

For a cohort first treated at $g = 3$ (reference period $t = 2$), every unit has both $Y_{i,2}$ and $Y_{i,4}$, so the direct long difference $Y_{i,4} - Y_{i,2}$ is observed for the entire sample. In a balanced panel, chained and ordinary estimates agree with IPW, no covariates and the same settings (`ptassumption(post)`, `controlgroup(nevertreated)`, `selection(none)`); balance alone does not make other specifications agree ([Get started](../get-started.md#compare-with-ordinary-did)).

### 2. Rotating panel

Each wave interviews a fresh cohort of units for a fixed number of consecutive periods (for example, two consecutive dates) before replacing them:

```text
Unit   t=1  t=2  t=3  t=4
  1     x    x    .    .
  2     .    x    x    .
  3     .    .    x    x
```

Here, **no unit** is observed at all four periods, and **no unit** is observed at both $t = 2$ and $t = 4$:

- `mode(ordinary)` with default balanced-panel handling (`allowunbalanced(false)`) drops every unit and cannot run.
- Even a two-period direct comparison between $t = 2$ and $t = 4$ finds zero units with both $Y_{i,2}$ and $Y_{i,4}$.
- However, adjacent changes $Y_{i,2} - Y_{i,1}$, $Y_{i,3} - Y_{i,2}$, and $Y_{i,4} - Y_{i,3}$ are all observed within units. By chaining the $2 \to 3$ and $3 \to 4$ within-unit changes, `mode(chained)` recovers $\text{ATT}(g=3, t=4)$ relative to reference period $t = 2$.

### 3. Attrition and late entry

Some units leave the sample permanently after a date, while others enter after period 1:

```text
Unit   t=1  t=2  t=3  t=4
  1     x    x    x    x
  2     x    x    x    .
  3     .    x    x    x
```

Every ordered pair $(s, t)$ has some units observed at both endpoints, alongside additional units observed on intermediate steps. Chained GMM uses **both** the direct long differences $(2, 4)$ from units 1 and 3 **and** the chained two-step path $(2, 3) + (3, 4)$ that incorporates unit 2 on $(2, 3)$.

### 4. General incomplete panel

Units miss sporadic dates, some units are observed only once (singletons, $B \setminus A$), and some units in the study population unit table $F$ have no observed outcomes at all ($F \setminus B$):

```text
Unit   t=1  t=2  t=3  t=4   Membership
  1     x    .    x    x    A (>= 2 dates)
  2     .    x    x    .    A (>= 2 dates)
  3     x    x    .    x    A (>= 2 dates)
  4     .    x    .    .    B \ A (singleton: 1 date)
  5     .    .    .    .    F \ B (unit table only: 0 dates)
```

- Units 1, 2, and 3 belong to $A$, the units observed at two or more dates, and supply within-unit changes $\Delta_{s,t} Y_i = Y_{i,t} - Y_{i,s}$ on the specific date pairs where both endpoints are observed.
- Unit 4 belongs to $B \setminus A$: it cannot supply any within-unit change $\Delta_{s,t} Y_i$, because differencing requires two observed dates.
- Unit 5 belongs to $F \setminus B$: it appears in the unit table supplied with `frame()` with baseline covariates and cohort membership, but has no observed outcomes.
- In `mode(chained)`, units 4 and 5 are **never discarded from the target population $F$**: they remain in the cohort shares, covariate distributions on $F$, observation-probability models on $F$, cluster sums, and influence functions, which have one row for each of the $N_F$ units ([Handling missing outcomes](../methods/missing-outcomes.md#the-study-population)).

## How chaining works: telescoping and earlier dates

Consider a cohort first treated at $g = 3$ with `anticipation(0)`, so its reference period is $t = 2$, and suppose we want to identify $\text{ATT}(3, 4) = E[Y_4(3) - Y_4(0) \mid G = 3]$.

### Forward telescoping across intermediate dates

In the fixed cohort population $G = k$, expectations are linear:

$$E[Y_4 - Y_2 \mid G = k] = E[Y_3 - Y_2 \mid G = k] + E[Y_4 - Y_3 \mid G = k]$$

Let $\Delta_{s,t}^{\text{DiD}}$ denote the adjusted difference-in-differences comparison between cohort $g = 3$ and the control group over pair $(s, t)$. Under conditional parallel trends from period 2 onward (`ptassumption(post)`) and representative pair changes within each cohort (`selection(none)`, or after missing-outcome adjustment under `selection(conditional_mtar)`):

- $\Delta_{2,3}^{\text{DiD}}$ identifies $\text{ATT}(3, 3)$,
- $\Delta_{3,4}^{\text{DiD}}$ identifies $\text{ATT}(3, 4) - \text{ATT}(3, 3)$,
- $\Delta_{2,4}^{\text{DiD}}$ directly identifies $\text{ATT}(3, 4)$.

Summing the two adjacent steps $\Delta_{2,3}^{\text{DiD}} + \Delta_{3,4}^{\text{DiD}}$ cancels $\text{ATT}(3, 3)$ and identifies $\text{ATT}(3, 4)$ even if no unit is observed at both $t = 2$ and $t = 4$. When both the direct pair $(2, 4)$ and the adjacent pairs $(2, 3)$ and $(3, 4)$ are observed, `cdid` includes all three in an overidentified GMM system and weights them optimally using their joint influence-function covariance ([GMM weighting and covariance regularization](../methods/gmm.md)).

### Cancellation through dates before the reference period

Under `ptassumption(post)`, conditional parallel trends are only assumed from the reference period $g - 1 - \text{anticipation}$ onward ($t = 2$ for $g = 3$). Even so, comparisons starting at an earlier date $t = 1$ remain informative whenever a cohort has observed pairs $(1, 2)$ and $(1, 3)$:

$$E[Y_3 - Y_1 \mid G = k] - E[Y_2 - Y_1 \mid G = k] = E[Y_3 - Y_2 \mid G = k]$$

In the GMM design matrix $A$, the comparison $\Delta_{1,2}^{\text{DiD}}$ loads on the unrestricted pre-reference differential $\delta_{3, 1 \to 2}$, while $\Delta_{1,3}^{\text{DiD}}$ loads on $\delta_{3, 1 \to 2} + \text{ATT}(3, 3)$. Differencing the two moments eliminates $\delta_{3, 1 \to 2}$ without assuming parallel trends between $t = 1$ and $t = 2$ ([Parallel trends](../methods/parallel-trends.md#comparisons-through-earlier-dates-under-ptassumptionpost)).

> [!IMPORTANT]
> **Population identity versus realized sample means:** Telescoping ($E[Y_4 - Y_2] = E[Y_3 - Y_2] + E[Y_4 - Y_3]$) and the cancellation through earlier dates hold for **population expectations** over the same cohort $G = k$ in $F$. When different subsets of units are observed on $(2, 3)$, $(3, 4)$, and $(2, 4)$, the realized sample means on those subsets will generally differ in a finite draw even when each pair is representative in mean change.

## Comparing `mode(chained)` and `mode(ordinary)` in Stata

In incomplete panels, ordinary mode by default keeps only the units observed in every period. With `allowunbalanced`, it uses the repeated cross-section estimators, with standard errors clustered by unit. The two approaches rest on different assumptions about the population and the sampling. Chained mode combines the changes of all units over the study population, under its assumption that the observed changes are representative. Using more units or obtaining a smaller standard error does not by itself show that an approach is valid.

With the data of the hands-on do-file in memory, the following specifications without covariates compare the three approaches. They illustrate the syntax; the data do not follow the observation patterns shown above:

```stata
cdid y, ivar(id) time(period) gvar(first_treat) ///
    method(ipw) frame(cdid_pop) bstrap(0) cband(0)

cdid y, ivar(id) time(period) gvar(first_treat) ///
    mode(ordinary) method(ipw) bstrap(0) cband(0)

cdid y, ivar(id) time(period) gvar(first_treat) ///
    mode(ordinary) allowunbalanced method(ipw) bstrap(0) cband(0)
```

---

**Reading path, step 3 of 22** · Previous: [Moving from did and csdid](../from-did.md) · Next: [Identify effects when the observed units change composition](identification.md)

[Guide home](../README.md) · [Package overview](../../README.md)
