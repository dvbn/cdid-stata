# Identify effects when the observed units change composition

This example builds a small rotating panel in which the treated units observed at each date have a different mix of a fixed, unobserved unit component. Outcome changes within the same unit still identify the treatment effects, while comparisons of mean outcome levels across dates do not. The reasoning is the argument about mean pair changes in [Estimate treatment effects from incomplete panels](incomplete-panels.md#forward-telescoping-across-intermediate-dates), applied to one concrete design.



## Question

A study follows units over five dates. Units in cohort 3 are first treated at date 3, and the others are never treated. Each unit is observed at exactly two consecutive dates. Among treated units, those with a high fixed component tend to be observed early and those with a low one late. Write ATT(3,t) for the average treatment effect of cohort 3 at date t. Can the changes observed within units identify ATT(3,3), ATT(3,4) and ATT(3,5) when the observed treated units differ in composition from date to date?

## The design

- **Cohort.** $G$ is 0 or 3, each with probability 1/2.
- **Fixed component.** $\alpha$ is $-1$ or $+1$, each with probability 1/2, independently of $G$. The estimators never see $\alpha$, and no covariate is supplied.
- **Untreated outcomes.** $Y_t(0) = \alpha + 0.1\,t + \varepsilon_t$. The $\varepsilon_t$ have mean zero, are independent across dates and are independent of $G$, $\alpha$ and observation.
- **Effect.** Let $\tau_t = 1$ when $G = 3$ and $t \ge 3$, and $\tau_t = 0$ otherwise. The actual outcome is $Y_t = Y_t(0) + \tau_t$, so every treated unit gains exactly 1 from date 3 onward and nothing earlier. ATT(3,3), ATT(3,4) and ATT(3,5) all equal 1, for the simulated units and for the population the design describes.
- **Observation.** Each unit is observed in one window $W$ of two consecutive dates, (1,2), (2,3), (3,4) or (4,5), drawn given $G$ and $\alpha$ with the probabilities in Table 1.

**Table 1. Window probabilities given cohort and fixed component (population algebra)**

| Units | Share of units | (1,2) | (2,3) | (3,4) | (4,5) |
| --- | --- | --- | --- | --- | --- |
| $G = 0$, either value of $\alpha$ | 1/2 | 0.25 | 0.25 | 0.25 | 0.25 |
| $G = 3$, $\alpha = +1$ | 1/4 | 0.40 | 0.30 | 0.20 | 0.10 |
| $G = 3$, $\alpha = -1$ | 1/4 | 0.10 | 0.20 | 0.30 | 0.40 |
| $G = 3$, both values together | 1/2 | 0.25 | 0.25 | 0.25 | 0.25 |

Each window contains treated and never treated units with both values of $\alpha$ with positive probability, and each window holds a quarter of each cohort. No unit is observed at all five dates: every unit has exactly two observed outcomes.

Let $R_t = 1$ when date $t$ lies in the unit's window. Write $p_t^{+}$ and $p_t^{-}$ for the probability that $R_t = 1$ in cohort 3 when $\alpha = +1$ and when $\alpha = -1$; each is the sum of the probabilities of the windows that contain $t$. Because the two values of $\alpha$ have equal mass in cohort 3,

$$m_t = E[\alpha \mid G = 3, R_t = 1] = \frac{p_t^{+} - p_t^{-}}{p_t^{+} + p_t^{-}}.$$

**Table 2. Observation and the mean fixed component by date (population algebra)**

| Quantity | $t = 1$ | $t = 2$ | $t = 3$ | $t = 4$ | $t = 5$ |
| --- | --- | --- | --- | --- | --- |
| $p_t^{+}$ | 0.4 | 0.7 | 0.5 | 0.3 | 0.1 |
| $p_t^{-}$ | 0.1 | 0.3 | 0.5 | 0.7 | 0.4 |
| $P(R_t = 1 \mid G = 3) = P(R_t = 1 \mid G = 0)$ | 0.25 | 0.5 | 0.5 | 0.5 | 0.25 |
| $m_t = E[\alpha \mid G = 3, R_t = 1]$ | 0.6 | 0.4 | 0 | $-0.4$ | $-0.6$ |
| $E[\alpha \mid G = 0, R_t = 1]$ | 0 | 0 | 0 | 0 | 0 |

The two cohorts are observed with the same probability at every date, so cohort 3 makes up half of the observed units at each date. What changes is the mix within cohort 3: its observed units have a mean $\alpha$ of 0.4 at date 2 and $-0.4$ at date 4.

## Why changes within units identify the effects

For a unit observed in window $(s, s+1)$, the fixed component cancels from its change:

$$Y_{s+1} - Y_s = 0.1 + (\tau_{s+1} - \tau_s) + \varepsilon_{s+1} - \varepsilon_s.$$

The change depends on the cohort, the dates and the noise, while the window depends only on the cohort and $\alpha$. Within each cohort, the units observed at both dates of a pair therefore have the cohort's mean change,

$$E[Y_{s+1} - Y_s \mid G, W = (s, s+1)] = E[Y_{s+1} - Y_s \mid G],$$

even though they do not represent the cohort in levels. This mean condition is what `selection(none)` requires, as described in [Handling missing outcomes](../methods/missing-outcomes.md#without-an-observation-model-selectionnone-default). Observation depends on $\alpha$, but $\alpha$ never enters a change, so no model of observation is needed.

Write $D_{s,s+1}$ for the mean change of cohort 3 minus the mean change of the never treated units. Every mean change equals 0.1, except the change of cohort 3 from date 2 to date 3, which is 1.1. Hence $D_{1,2} = 0$, $D_{2,3} = 1$, $D_{3,4} = 0$ and $D_{4,5} = 0$.

Untreated outcomes satisfy $Y_t(0) - Y_2(0) = 0.1\,(t - 2) + \varepsilon_t - \varepsilon_2$ in both cohorts, so PT-post from reference date 2 holds, and there is no anticipation. The first link then identifies the effect, and the later links identify its increments:

$$D_{2,3} = \text{ATT}(3,3) = 1, \qquad D_{3,4} = \text{ATT}(3,4) - \text{ATT}(3,3) = 0, \qquad D_{4,5} = \text{ATT}(3,5) - \text{ATT}(3,4) = 0.$$

Summing the links gives $\text{ATT}(3,4) = D_{2,3} + D_{3,4} = 1$ and $\text{ATT}(3,5) = D_{2,3} + D_{3,4} + D_{4,5} = 1$. The four windows are the only pairs anyone is observed at, so chained difference in differences (chained DiD) can only combine these consecutive comparisons. The cell at date 2, reported before treatment, uses $D_{1,2}$.

## Why comparisons of levels do not

An intercept-only repeated cross section (RCS) IPW estimator compares normalized means of outcome levels in each cohort and date arm. Its constant propensity odds cancel within each normalized control mean. In this design the population propensity is 1/2, away from the trimming threshold, so the population comparison for cell (3,t) with reference date 2 is

$$\bigl(E[Y_t \mid G = 3, R_t = 1] - E[Y_2 \mid G = 3, R_2 = 1]\bigr) - \bigl(E[Y_t \mid G = 0, R_t = 1] - E[Y_2 \mid G = 0, R_2 = 1]\bigr).$$

The common trend cancels, and the observed never treated units have a mean $\alpha$ of zero at every date, so for $t \ge 3$ this contrast equals $1 + m_t - m_2$, with $m_2 = 0.4$.

**Table 3. What each comparison targets (population algebra)**

| Cell | Known value | Chained links | Comparison of observed levels |
| --- | --- | --- | --- |
| Date 2, before treatment (dates 1 and 2) | 0 | $D_{1,2} = 0$ | $m_2 - m_1 = 0.4 - 0.6 = -0.2$ |
| ATT(3,3) | 1 | $D_{2,3} = 1$ | $1 + (0 - 0.4) = 0.6$ |
| ATT(3,4) | 1 | $1 + 0 = 1$ | $1 + (-0.4 - 0.4) = 0.2$ |
| ATT(3,5) | 1 | $1 + 0 + 0 = 1$ | $1 + (-0.6 - 0.4) = 0$ |

Callaway and Sant'Anna's repeated cross section setting assumes that the sample at each date represents the population at that date and that the distribution of cohorts and covariates is stable across dates. Here population parallel trends hold and the cohort shares are stable, but the first condition fails for cohort 3. Its observed units come mostly from $\alpha = +1$ at early dates and mostly from $\alpha = -1$ at late dates, so each level contrast adds a change in composition to the effect.

Fixed differences between units do not cause this on their own. The never treated units take the same two values of $\alpha$, but their windows do not depend on it, so their observed levels represent their cohort at every date. If the windows of cohort 3 were also drawn independently of $\alpha$, as in the rotating example of [Estimate treatment effects from incomplete panels](incomplete-panels.md), the observed levels would represent both cohorts, and the RCS approach would target the same effects under its assumptions. Ordinary unbalanced mode uses RCS routines; whether it answers the question depends on its sampling conditions.

## Assumptions in this example

| Item | In this example |
| --- | --- |
| Whose effect | ATT(3,t) for t = 3, 4, 5 in the original study population |
| Comparison group | Units never treated, G=0 |
| Anticipation | None, so the reference period of cohort 3 is date 2 |
| Parallel trends | PT-post from date 2, true in the design |
| Observed changes | `selection(none)`: within each cohort, units observed at both dates of a pair have the cohort's mean change; true because observation depends only on the cohort and $\alpha$ |
| Treatment adjustment | Treatment weighting by inverse probability weighting (IPW) with an intercept only, no covariates; adequate because parallel trends hold without covariates |
| Overlap and support | Every window has positive probability in both cohorts; in a sample, a window with no observed unit makes the effects that depend on it unavailable |
| Unit table | One row per unit with `id` and `g`; no covariates or weights |


## Stata usage

The tables on this page give population values, not estimates. [Incomplete panels](incomplete-panels.md#comparing-modechained-and-modeordinary-in-stata) applies the three estimation approaches to the data of the hands-on do-file, which follow a different design; their estimates are not the values in these tables. Representative mean changes are conditions on population expectations, not equal sample means, equal levels, equal fixed components or equal distributions.

---

**Reading path, step 4 of 22** · Previous: [Estimate treatment effects from incomplete panels](incomplete-panels.md) · Next: [How to compare precision](precision.md)

[Guide home](../README.md) · [Package overview](../../README.md)
