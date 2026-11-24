# GMM weighting and covariance regularization

This page describes the moment conditions of chained estimation, how the comparisons are weighted, and how their covariance is regularized.

## Moment conditions and the quadratic objective

Let $\hat{d}$ denote the $M \times 1$ vector of two-period comparison estimates used by `cdid`. Each entry of $\hat{d}$ compares outcome changes between an ordered pair of dates for a target cohort with its control group (never-treated or not-yet-treated units, including comparisons that use the same pair of dates again with another control group), constructed using the selected treatment adjustment (`method()`) and missing-outcome correction (`selection()` and `selmethod()`).

Let $\theta$ denote the GMM model parameters: post-treatment effects and any pre-period trend differentials left free by the chosen parallel-trends model. PT-all restricts the differentials before the reference period; these model parameters are distinct from the direct, unrestricted pre-period cells reported for the pretrend test. Let $d$ denote the vector of population comparison means and $A$ the known $M \times p$ design mapping $\theta$ to them. Under sufficient parallel-trends, observation, overlap and model conditions, with the required regularity,

$$
d = A\theta, \qquad \hat{d} \xrightarrow{p} d
$$

This is a population moment restriction and a large-sample statement, not a claim that nuisance-adjusted sample comparisons are unbiased. Given a symmetric positive-semidefinite $M \times M$ weighting matrix $W$, the chained GMM estimator minimizes the quadratic form

$$
Q_W(\theta) = (\hat{d} - A\theta)' W (\hat{d} - A\theta),
$$

and recovers each supported post-treatment target $\tau_j = c_j'\theta$ as a linear contrast of the fitted parameter vector $\hat{\theta}$.

Several properties follow directly from this formulation:

- **Sufficient conditions for recovery.** Validity of every comparison used under the assumed parallel-trends, observation-selection, overlap and nuisance-model conditions, identification of $c_j$ from the row space of $A$, and standard GMM regularity are sufficient for consistency. They are not necessary in every special case: misspecification can cancel in a target contrast. Such cancellation is not a general robustness guarantee.
- **What $W$ changes when moments are valid.** When the full moment conditions hold and $W$ converges in probability to a positive-definite limit $W_0$ that preserves the rank of $A$, changing $W$ can change the linear combination of overidentifying comparisons used to reconstruct $\tau_j$ and its asymptotic and finite-sample precision. It does not change the population target.
- **No finite-sample unbiasedness or arbitrary-$W$ guarantee.** Feasible two-step GMM estimates $W$ from the same sample used to form $\hat{d}$. Inserting a data-dependent weighting matrix does not establish exact finite-sample unbiasedness, and asymptotic validity does not hold for arbitrary unstable or data-adaptive choices of $W$.
- **Sensitivity to incompatible moments.** If the population limits of the comparisons used cannot all be represented by one $A\theta$ (for example after invalid trend or observation restrictions), the moments are incompatible. Misspecification alone need not imply incompatibility or bias in every contrast. Under the relevant convergence conditions, the limiting criterion defines a weighted pseudo-target that generally depends on the limiting weights.
- **Moment-covariance ridge is not a parameter penalty.** When the estimated moment covariance is regularized before inversion to form $W$, the positive ridge is added to the *moment covariance matrix*, not to $A'WA$ and not as an $L_2$ penalty on $\theta$. It imposes no shrinkage target on ATT and supplies no identifying restrictions on $\theta$. A particular estimate can nevertheless move closer to or farther from zero when its comparison weights change.

## Weighting and covariance options

Built-in chained cohort IPW, outcome regression and DR default to `gmmweighting(within)`. It weights comparisons using the inverse of the covariance of within-unit changes implied by a working model. A date shared by two comparisons contributes to their covariance. For example, a shock in period 2 enters the change from 1 to 2 positively and the change from 2 to 3 negatively. The working covariance uses this structure and the unit weights. It determines the weights directly and gives the fixed-effects weighting in [matched regression designs](cross-pair.md).

`gmmweighting(optimal)` uses feasible two-step weighting based on the estimated covariance of comparison scores. Use it to reproduce the earlier default for built-in chained cohort fits. `two_step` is a synonym. Switching fits, custom estimators, and the IPW and DR branches of `preset(cohort_time_fe)` default to `optimal`. The regression branch of the preset uses within-unit least squares. `gmmweighting(identity)` minimizes the unweighted criterion over all comparisons where available. Inspect `e(gmm_weighting)` or run `cdid_diagnostics` after a chained cohort fit to see the applied weighting.

`gmmcovariance(tikhonov)` is the only covariance regularization available, with `condition(1e-10)` and `ridge(1e-4)` as fixed defaults. These are numerical settings, not identification thresholds or universally optimal tuning choices.

### Within-unit weighting

Holding the nuisance fits fixed, each comparison is a weighted sum of unit outcome changes between its two dates, $\hat{d}_m = \sum_i \omega_{mi}(Y_{i t_m} - Y_{i s_m})$, where the loadings $\omega_{mi}$ are positive weights on the target cohort's units and the IPW, OR or DR loadings on the control units. Suppose outcomes are a unit effect plus level shocks that are uncorrelated across units and dates, with variance proportional to $1/w_i$ for a unit with sampling weight $w_i$. `gmmweighting(within)` uses the working weighting matrix built from

$$
V_{ml} = 1\{g_m = g_l\}\, \sigma_{ml} \sum_i \frac{\omega_{mi}\,\omega_{li}}{w_i}, \qquad
\sigma_{ml} = 1\{t_m = t_l\} + 1\{s_m = s_l\} - 1\{t_m = s_l\} - 1\{s_m = t_l\},
$$

where $g_m$ is the target cohort of comparison $m$ and $(s_m, t_m)$ its dates. For two comparisons of the same target cohort, $V_{ml}$ is proportional to their covariance under these assumptions. Comparisons of different target cohorts that share control units are also correlated, but $V$ sets these entries to zero: it is a working matrix, block diagonal by target cohort. `gmmweighting(within)` uses $V$ in place of the estimated covariance, after the exact reductions described next; the choice of weighting affects efficiency, not consistency. Comparisons of the same target cohort, control group and set of observed units that form cycles over dates are reduced to a spanning forest of each group, with $G = E'E$ for the removed comparisons. The weighting matrix does not depend on the outcome data, and standard errors still use the estimated unit influence functions, so inference does not rely on the level-shock model. It is available with or without `selection(conditional_mtar)`, with every `selmethod()`, with `crosspair` and with covariates or weights that vary over time; with observation correction, `outcomepooling(pooled)`, or covariates or weights that vary over time, it requires a constant in the covariate design. It is not available with `method(custom)` or `dvar()`.

## Exact reductions and the metric of the full set of comparisons

Before examining the sample covariance of the comparisons, `cdid` removes comparisons that are exact linear combinations of others, when exact identities allow the full set to be reconstructed from the remaining comparisons (for example, with built-in IPW and no correction for missing outcomes, when distinct comparisons share the same cohort, control group and observed units).

Suppose the full set contains $M$ comparisons and exact identities express the $M \times 1$ residual vector in terms of $L \le M$ remaining comparisons through a sparse $M \times L$ reconstruction matrix $E$:

$$
\hat{d} - A\theta = E(\hat{d}_r - A_r\theta), \qquad G = E'E = C'C,
$$

where $G$ is the symmetric positive-definite $L \times L$ Gram matrix of the reduction and $C$ is its upper-triangular Cholesky factor. When no reduction applies, $L = M$ and $G = C = I_M$.

Preserving $G = E'E$ is essential for three reasons:

1. **Exact preservation of the original $M$-comparison objective.** Any weighting matrix $W$ of the full set becomes $E'WE$ on the $L$ remaining comparisons. In particular, under `gmmweighting(identity)`, the unweighted objective over the full set is
   $$
   (\hat{d} - A\theta)' I_M (\hat{d} - A\theta) = (\hat{d}_r - A_r\theta)' E'E (\hat{d}_r - A_r\theta) = \|C(\hat{d}_r - A_r\theta)\|_2^2.
   $$
   The reduced system is therefore whitened by $C$, so that `gmmweighting(identity)` still minimizes the exact sum of squares over all $M$ original comparisons, rather than replacing $G$ by a naive $I_L$ in reduced coordinates.
2. **No arbitrary row dropping or spanning-tree subsetting.** Correlation alone does not justify removing comparison rows; an exact reduction requires an exact identity, which can involve linear combinations rather than identical rows. Discarding correlated rows, deduplicating rows without preserving $G = E'E$, or restricting the set to a minimal set of comparisons that identifies the effects would change the estimator, and `cdid` does none of these.
3. **Algebraic identity is not extra sample information.** When two comparisons evaluate the exact same sample difference on the same units, or when consecutive changes of the same balanced units telescope to a long change, the extra row repeats existing sample variation rather than adding independent information. Exact reduction removes this duplication without changing the criterion.

## Covariance conditioning and Tikhonov regularization

Under `gmmweighting(optimal)` (or `two_step`), let $\widehat{V}_r$ denote the $L \times L$ covariance of the comparisons on the scale used for weighting, formed from the unit-level influence functions of the comparisons. This normalization need not equal the estimated covariance of $\hat d_r$ itself. A common positive rescaling of $\widehat V_r$ rescales its trace-proportional ridge equally and multiplies the criterion by a constant, leaving its minimizer unchanged. To respect the original $M$-comparison geometry when $L < M$, conditioning and regularization are defined in the $C$-transformed coordinates where the identity metric of the full set, $G = C'C$, becomes $I_L$. Define the transformed $L \times L$ symmetric covariance matrix

$$
B = C \widehat{V}_r C'.
$$

Adding an isotropic ridge $\lambda I_M$ in the original $M$-comparison coordinates corresponds in the reduced transformed coordinates to replacing $B$ by $B + \lambda I_L = C(\widehat{V}_r + \lambda G^{-1})C'$. The effective reduced weighting matrix is therefore

$$
W_r(\lambda) = C'(B + \lambda I_L)^{-1} C = \left(\widehat{V}_r + \lambda G^{-1}\right)^{-1},
$$

which reduces to $(\widehat{V}_r + \lambda I_M)^{-1}$ when $G = I_M$, and otherwise scales the ridge inversely with $G$ so that identical original comparisons receive the same metric weight regardless of how exact reduction grouped them.

### How the weighting is chosen

After checking that $\widehat{V}_r$ is finite and exactly symmetric and forming a symmetric finite $B$, the estimator computes a symmetric eigendecomposition of $B$ with minimum and maximum eigenvalues $d_{\min}(B)$ and $d_{\max}(B)$, checks that the eigendecomposition residual lies within floating-point roundoff bounds, and evaluates the reciprocal spectral condition ratio

$$
\rho = \frac{\max(d_{\min}(B), 0)}{d_{\max}(B)}
$$

(with $\rho = 0$ if $B = 0$). Materially negative eigenvalues or diagonal entries abort the fit; small negative eigenvalues consistent with floating-point roundoff are accepted only when $\widehat{V}_r$ was computed as a cross-product of unit influence functions, which makes it positive semidefinite in exact arithmetic.

Estimation stops with an error if the eigendecomposition fails its numerical accuracy checks.

The weighting then follows one of four cases; `cdid_diagnostics` reports which one applied:

| Case | Condition | Effective reduced weight and solution |
| --- | --- | --- |
| **Well-conditioned unregularized solve** | $B$ is positive definite, $\rho \ge \tau_{\mathrm{cond}}$, the reduced dimension $L$ does not exceed the rank bound given by the number of contributing units, and ordinary Cholesky factorization of $\widehat V_r$ and the subsequent solve pass their checks | $W_r(0) = C' B^{-1} C = \widehat{V}_r^{-1}$, implemented by Cholesky whitening of $\widehat V_r$ and pivoted QR factorization on the whitened design without regularization |
| **Tikhonov covariance ridge** | $B \ne 0$ and $\rho < \tau_{\mathrm{cond}}$, $d_{\min}(B) = 0$, $d_{\min}(B) < 0$ within roundoff, $L$ exceeds the sample rank upper bound, or unregularized Cholesky of $\widehat V_r$ fails after valid covariance diagnosis | $W_r(\lambda) = C'(B + \lambda I_L)^{-1}C$ with $s = \mathrm{tr}(B)/M$ and $\lambda = \alpha_{\mathrm{ridge}} s$, solved by Cholesky factorization of $B + \lambda I_L$ and pivoted QR |
| **Exact-zero covariance limit** | $\widehat{V}_r = 0$ and $B = 0$ identically | Scale-free least squares over all comparisons, with $\lambda = 0$ and reduced weight $G = C'C$ (whitened by $C$), paired with the actual zero sandwich covariance and zero-scale inference diagnostics |
| **Error: invalid covariance or failed solution** | $\widehat{V}_r$ or $B$ is nonfinite or materially indefinite, or the Cholesky or whitened-QR factorization fails its conditioning or residual bounds | Estimation stops with an error; no estimate is returned |

Note three implementation details of this rule:

- **The trace scale uses the number $M$ of comparisons in the full set.** In the Tikhonov case, the scale is $s = \mathrm{tr}(B) / M$, where $M$ is the number of comparisons in the full set and $\mathrm{tr}(B) = \mathrm{tr}(\widehat{V}_r G)$ equals the trace of the full $M \times M$ moment covariance.
- **Stable factorizations rather than explicit matrix inversion.** Neither the unregularized case nor the Tikhonov case forms an explicit dense inverse $B^{-1}$ or $(B + \lambda I_L)^{-1}$. The unregularized case whitens with the Cholesky factor of $\widehat V_r$; the ridge case uses that of $B+\lambda I_L$. Both then apply pivoted QR to the whitened design, and verify factor conditioning and relative residuals before accepting the solution.
- **No automatic repairs.** The default values `condition(1e-10)` and `ridge(1e-4)` are numerical implementation settings, not universal optimal tuning parameters and not thresholds for parameter identification. If $B + \lambda I_L$ or the whitened QR factor still fails its conditioning checks, estimation stops with an error rather than increasing $\lambda$ iteratively, clipping eigenvalues, adding random jitter, dropping comparisons, or switching automatically to `gmmweighting(identity)`.

## Covariance singularity versus parameter identification

Singular or near-singular moment covariance $\Sigma = \mathrm{Cov}(\hat{d})$ and non-identification of a target parameter $\tau = c'\theta$ are distinct mathematical phenomena. Moment covariance is a property of the stochastic errors across comparison rows. Parameter identification is a property of the design matrix $A$ and the maintained population restrictions $d = A\theta$.

The two analytical examples below illustrate why neither implies the other. Here $\hat d$ denotes noisy comparisons and $d = E[\hat d]$ their population means. These are algebraic examples.

### Example 1: Singular moment covariance with an identified parameter

Suppose a scalar treatment effect $\theta$ is measured by two comparison moments that carry the exact same zero-mean shock $u$ with $\mathrm{Var}(u) = \sigma^2 > 0$:

$$
\hat d_1 = \theta + u, \qquad \hat d_2 = \theta + u, \qquad A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \qquad \Sigma = \sigma^2 \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}.
$$

Here the design matrix $A$ has column rank $1$, so $\theta = E[\hat d_1] = E[\hat d_2]$ is identified. Yet the $2 \times 2$ covariance matrix $\Sigma$ has rank $1$ and determinant $0$ because $\hat d_1 - \hat d_2 = 0$ with probability one. In chained estimation, a singular or near-singular covariance of this kind can arise when duplicate or telescoping comparison paths span the same dates. If the duplication is an exact structural equivalence, reduction collapses the pair to $L = 1$ with $G = E'E = 2$ and positive scalar variance. If this duplicated pair is kept without reduction, Tikhonov regularization replaces $\Sigma$ with $\Sigma + \lambda I_2$, which has eigenvalues $2\sigma^2 + \lambda$ and $\lambda > 0$, allowing the GMM step to form a well-defined convex combination of $\hat d_1$ and $\hat d_2$. Because both moments share the same shock $u$, the second moment adds no independent information, and the estimator variance remains $\sigma^2$.

### Example 2: Nonsingular moment covariance with an unidentified parameter

Suppose instead that two parameters $\theta = (\theta_1, \theta_2)'$ (for instance, an unrestricted pre-period trend differential $\theta_1$ and a post-treatment effect $\theta_2$) enter a single available comparison moment with $\mathrm{Var}(u) = \sigma^2 > 0$:

$$
\hat d_1 = \theta_1 + \theta_2 + u, \qquad A = \begin{pmatrix} 1 & 1 \end{pmatrix}, \qquad \Sigma = \begin{pmatrix} \sigma^2 \end{pmatrix}.
$$

Here the $1 \times 1$ moment covariance $\Sigma = (\sigma^2)$ is strictly positive definite and well conditioned, but $A$ has rank $1 < 2$. The moment identifies the sum $\theta_1 + \theta_2$, while the individual contrast $\tau = \theta_2 = (0, 1)\theta$ is not in the row space of $A$. Adding a covariance ridge $\lambda > 0$ changes the scalar weight from $1/\sigma^2$ to $1/(\sigma^2 + \lambda)$, so the normal matrix

$$
A'(\Sigma + \lambda)^{-1}A = \frac{1}{\sigma^2 + \lambda} \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
$$

still has rank $1$ and cannot be inverted for $\theta_1$ and $\theta_2$ separately. Because covariance regularization acts on $\Sigma$ rather than adding rows to $A$ or penalizing $\theta$, it never turns an unidentified parameter into an estimable one.

### Support, errors, and precision

In `cdid`, the availability of post-treatment effects is determined from the design $A$, before and independently of covariance regularization. Pretreatment cells, which are direct comparisons, have their own support rules:

- **Unsupported cells.** When a target contrast $c_j'\theta$ does not lie in the row space of $A$ given the comparisons available and the chosen `ptassumption()`, the cell is unavailable: its estimate and standard error are missing, and the results record the reason. This is a statement about the comparisons available in the sample, not a claim that the effect is unidentifiable in the population under every possible assumption.
- **Errors in the models or computations.** Perfect separation, or a rank or support failure that remains in a required model after its permitted preprocessing, or a covariance or QR factorization that fails the applicable checks, stops the call with an error rather than returning an unavailable result for affected cells or dropping comparisons.
- **No uniform precision dominance.** With independent units, valid moments, consistent covariance estimation, standard regularity and a nonsingular limiting covariance after exact reductions, the unregularized optimal weighting has the usual asymptotic efficiency property within the set of comparisons used. That statement does not automatically extend to the ridge case or to dependent units. Covariance ridge changes the weighting to obtain numerical conditioning when $B$ is singular or ill-conditioned, and `gmmweighting(identity)` uses fixed weights over all comparisons. In finite samples, estimation noise in $\widehat{V}_r$ can offset theoretical efficiency gains, so neither `gmmweighting(optimal)`, nor a particular `ridge()`, nor `gmmweighting(identity)` is guaranteed to produce smaller finite-sample variance or mean squared error in every setting.


## Sandwich inference and clustering

The weighting matrix determines the point-estimate loadings. It is not inverted to report parameter uncertainty. The influence functions have one row per unit of the study population, including units without an observed pair, and account for the estimation of the adjustment models and, in aggregations, of the cohort shares.

`e(V)` is analytical estimator covariance. `e(V_asymptotic)` is `e(N_full)` times that covariance; `e(analytic_se)` is the square root of the diagonal of `e(V)`. The reported `e(se)` can instead come from the multiplier bootstrap. With `noinffunc`, the GMM weighting is computed as usual, so the point estimates are unchanged, but no inference is returned and `estat` cannot aggregate the results.

Clustering sums unit influence functions for inference. It does not change unit-based GMM weights or point estimates, and optimal weighting is not claimed to be cluster-optimal. First-order inference requires valid moments, regular estimators of the comparisons, consistent covariance estimates and a stable limiting case of the weighting. No uniform boundary or few-cluster coverage guarantee follows from a successful fit.

## Examples and computational limits

With the data of the hands-on do-file in memory, compare the default within weights, the earlier optimal weights and the unweighted criterion under analytical inference:

```stata
cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) bstrap(0) cband(0)
cdid_diagnostics, detail

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) gmmweighting(optimal) bstrap(0) cband(0)
cdid_diagnostics, detail

cdid y x_bin x_cont, ivar(id) time(period) gvar(first_treat) ///
    method(dr) gmmweighting(identity) bstrap(0) cband(0)
cdid_diagnostics, detail
```

`workspace()` defaults to 1073741824 bytes (1 GiB). Without `chunk()`, the number of rows processed at a time is set automatically. Without `workcap()`, the number of operations is not limited; `workcap()` counts operations conservatively, not seconds. These limits do not bound the total memory or time used by Stata and do not relax any statistical check.

---

**Reading path, step 15 of 22** · Previous: [Add PT-all restrictions](../examples/pt-all-comparison.md) · Next: [Compare model robustness across treatment and missing outcome adjustments](../examples/robustness.md)

[Guide home](../README.md) · [Package overview](../../README.md)
