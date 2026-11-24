# cdid for Stata 1.4.0.0

- Use `gmmweighting(within)` by default for built-in chained cohort IPW, outcome regression, and doubly robust fits, including supported cross-pair and missing-outcome specifications. The weights account for shared dates using the working covariance of within-unit changes under uncorrelated level shocks. This gives the weighting of matched fixed-effects regressions. Use `gmmweighting(optimal)` or its alias `two_step` to reproduce the earlier default based on the empirical covariance of comparison scores. Estimates and standard errors can differ when multiple comparisons identify an effect. Inference uses the full study-unit influence functions.
- Switching fits, custom estimators, and the IPW and DR branches of `preset(cohort_time_fe)` default to `optimal`. The regression branch of the preset uses within-unit least squares. See [GMM weighting](docs/methods/gmm.md) for the supported specifications and how to inspect the applied weighting.

- Correct the ridge decision for within weighting when the number of comparisons exceeds the empirical influence-function covariance rank bound. Regularization of the working covariance uses its own conditioning and factorization checks.

# cdid for Stata 1.0.0.1

- Reject binary switching fits when a selection covariate separates an observation logit. The check applies to direct-pair and independent-date models used for observation weighting and augmented correction.

# cdid for Stata 1.0.0.0

The first Stata release follows [R cdid 1.3.0.0](https://github.com/dvbn/cdid/releases/tag/v1.3.0.0).

- Estimate chained effects with staggered adoption or with treatment entry and exit.
- Use inverse probability weighting, outcome regression or doubly robust adjustment, with optional observation weighting, change imputation or augmented correction for missing outcomes.
- Aggregate effects, plot event studies, test conditional pretrends and inspect comparison and nuisance-model diagnostics.
- Use the separate csdid 2.0.2 installation for ordinary difference-in-differences.
