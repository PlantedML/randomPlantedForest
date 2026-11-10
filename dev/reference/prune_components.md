# Prune the components of a Random Planted Forest

Drops the small components of a regression forest's functional
decomposition and, optionally, refits one coefficient per kept
component. The rule and the threshold are chosen out of bag by default,
from the forest itself without refits, or by K-fold cross-validation on
the training rows.

## Usage

``` r
prune_components(
  object,
  data,
  y = NULL,
  ...,
  rule = c("D", "R"),
  tau = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1),
  select = c("1se", "min"),
  select_by = c("oob", "cv"),
  folds = 5L,
  fold_id = NULL,
  seed = NULL,
  nthreads = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  `[rpf]`: A fitted regression
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- data:

  `[data.frame | matrix]`: The training data of `object`, in the order
  of the fit: the predictors as passed to
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  as `x`, or the `data` of the formula or recipe interface. For a
  formula or recipe fit, `data` may contain the outcome.

- y:

  `[numeric | NULL: NULL]`: The outcome of the training rows. Required
  unless `data` contains the outcome of a formula or recipe fit.

- ...:

  Reserved for future expansion. Unknown arguments are an error.

- rule:

  `[character: c("D", "R")]`: The rules to consider: `"D"` drops, `"R"`
  drops and refits one coefficient per kept component.

- tau:

  `[numeric: c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1)]`: The
  thresholds to consider, as shares of the variance of the prediction
  (`0.01` = 1%), each at least 0 and below 1.

- select:

  `[character(1): "1se"]`: How the pair is chosen: `"1se"` (the largest
  `tau` within one standard error of the lowest RMSE) or `"min"` (the
  lowest RMSE).

- select_by:

  `[character(1): "oob"]`: How the candidates are scored, and R's
  coefficients fitted: `"oob"` (out of bag, from one purification of
  `object` and no refits) or `"cv"` (K-fold CV with refits of the
  forest; R's coefficients in sample); see Details.

- folds:

  `[integer(1): 5]`: The number of folds, used if `fold_id` is `NULL`:
  of the CV, or out of bag of R's refit and of the SE.

- fold_id:

  `[vector | NULL: NULL]`: The fold of every training row: a vector of
  length `nrow(data)` with at least two distinct values. If `NULL`, the
  rows are assigned to `folds` folds at random (after `set.seed(seed)`).

- seed:

  `[integer(1) | NULL: NULL]`: The seed set before the random fold
  assignment and, for the CV, before every fold's refit. If `NULL`, one
  is drawn from the current random number stream and stored in the
  result, so that the choice can be reproduced. The caller's random
  number state is restored afterwards (see Details).

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads for purifying the
  forests (`object`, and the CV's fold fits) and for the fold fits.
  `NULL` uses the `nthreads` `object` was fitted with. The result does
  not depend on it.

- verbose:

  `[logical(1): FALSE]`: Report the progress of the choice.

## Value

An object of class `rpf_pruned`, a list with elements:

- `fit`: `object` (purified).

- `rule`, `tau`: the chosen rule and threshold.

- `intercept`: `alpha` (the forest's intercept `f_0` for `"D"`).

- `beta`: named vector of the coefficients of the kept components (all 1
  for `"D"`). For `"R"` with `select_by` `"oob"`, `intercept` and `beta`
  are fitted on the out-of-bag components (see Details).

- `sets`: the predictor sets of the kept components, as column indices.

- `components`: a `data.frame` of every nonzero component of `object` on
  the training rows: `component` (its name, as in
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)),
  `order`, `share`, `kept`, `beta` (`NA` if dropped).

- `oob`: `NULL` unless the pair was chosen out of bag; otherwise a
  `data.frame` with one row per candidate: `rule`, `tau`, `oob_rmse`
  (corrected; used for the choice), `se_mse` (from its per-fold MSEs),
  `oob_rmse_raw` and `se_mse_raw` (uncorrected), `n_kept`, `min` and
  `one_se` (the two choices). Its attribute `n_rows` is the number of
  rows with an OOB family.

- `cv`: `NULL` unless the pair was chosen by CV; otherwise a
  `data.frame` with one row per candidate: `rule`, `tau`, `cv_rmse`,
  `se_mse` (SE of the CV MSE), `n_kept` (mean over the folds), `min` and
  `one_se`.

- `select`, `select_by`, `fold_id`, `seed`: as used (`fold_id` and
  `seed` are `NULL` if no choice was made). For a given rule and `tau`,
  `select_by` is `"oob"` if R's coefficients were fitted out of bag,
  otherwise `NULL`.

Use
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
and
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
on it, and
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
to save it.

## Details

The forest is decomposed into its purified components `f_S`, one per set
`S` of predictors (see
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
and
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)),
at all orders up to `max_interaction`, and its intercept `f_0`. The size
of a component is its share of the variance of the prediction `f` on the
training rows, `share_S = Var(f_S) / Var(f)`. The shares are descriptive
sizes, not a partition: the components are not orthogonal on the data,
so they need not add up to 1. A component is kept if its share is at
least `tau` (and it is not identically zero on the training rows). Two
rules are available:

- `"D"` (drop): the prediction is `f_0` plus the kept components.
  `rule = "D"` with `tau = 0` is the unpruned forest.

- `"R"` (drop and refit): the prediction is `alpha + sum(beta_S * f_S)`
  over the kept components, with the intercept `alpha` and the
  coefficients `beta_S` by least squares of `y` on the kept components
  (the minimum-norm solution if the components are collinear). The
  components themselves are not changed, only their scale. With
  `select_by = "oob"` (the default, see below), the least squares use
  each training row's out-of-bag components instead of the forest's own,
  for a chosen and for a given rule and `tau`. The forest's own
  components have fitted the training rows' noise, so least squares on
  them scale the components up (on California housing the three-way
  interactions by about 1.3, where new data wants about 0.9). With
  `select_by = "cv"`, and where the default falls back to it (see
  below), the least squares use the forest's own components on the
  training rows.

**The choice.** If more than one pair of `rule` and `tau` is given (the
default: both rules and the grid 0, 0.01%, 0.03%, 0.1%, 0.3%, 1%, 3%,
10%), the pair is chosen by its out-of-bag (`select_by = "oob"`, the
default) or cross-validated (`"cv"`) RMSE. With `select = "1se"` (the
default), the largest `tau` whose MSE is within one standard error of
the minimum is chosen (the SE is the standard deviation of the per-fold
MSEs of the minimum, divided by `sqrt(K)`), and at that `tau` the rule
with the lower RMSE. With `select = "min"`, the pair with the lowest
RMSE is chosen (ties, up to a relative `1e-12`: the larger `tau`, then
`"D"`). `"1se"` gives sparse models: with `"R"` fitted out of bag the
out-of-bag error changes little over a range of `tau`, and `"min"` keeps
small components that add nothing on new data.

**Out of bag** (`select_by = "oob"`). Every tree family of `object` is
grown on a bootstrap sample of the training rows, and the purification
works family by family. So `object` is purified once more (as
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
does), and on the way every family's components are evaluated at the
training rows outside its bootstrap sample. A row's out-of-bag (OOB)
component is the mean over its OOB families (about 37% of them), and its
OOB intercept likewise. The purification uses the predictors of all
training rows, not the outcome, so no outcome information of a row
enters its OOB components.

- The kept components of every `tau` are those of the final model (the
  shares of `object` on the training rows).

- `"D"` is scored by the OOB intercept plus the kept OOB components.

- For `"R"`, `alpha` and `beta_S` are fitted as for the final model, by
  least squares on the OOB components, but only on the rows outside each
  of K folds (`folds`, `fold_id`, `seed`). They are scored with the OOB
  components of the fold's rows, so that no row scores coefficients
  fitted on it. The final model's are fitted on all rows with at least
  one OOB family.

- An OOB mean averages `m_i` families (about 37% of `B = ntrees`), the
  final model all `B`, so the OOB prediction scatters more around its
  expectation. The OOB MSE is therefore corrected: per row it subtracts
  `(1 / m_i - 1 / B) * v_i`, with `v_i` the across-family variance of
  the candidate's prediction, estimated from the OOB means over the even
  and over the odd families (for a row whose OOB families are all even
  or all odd: the mean estimate of the other rows). Without the
  correction, the noise of the OOB components makes the candidates look
  worse than they are. The uncorrected values are reported as
  `oob_rmse_raw`.

- The OOB RMSE is taken over the rows with at least one OOB family, the
  SE for `"1se"` from the per-fold MSEs.

- With few families the OOB means average only one to a few families,
  and many rows have none (16% at 4 families). In a simulation (p = 6, n
  = 800, 60 replicates) the corrected OOB error of the unpruned `"D"`
  was within 1.5% of the test error at 4, 10, 20 and 50 families, that
  of `"R"` within 4.4% except once (+11.5% at 10 families with a 3-way
  interaction), and the OOB choice was as good as the CV's or better
  from 10 families on (within 0.4% in test error, or better by 7-10%
  where the CV's fold fits choose poorly). Below 10 families (`ntrees`),
  the default `select_by` scores by CV instead (with a message); an
  explicit `select_by = "oob"` warns.

- A deterministic forest has no out-of-bag rows, and a forest restored
  by
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  from a blob made by an earlier version has no in-bag counts (with a
  message): with the default `select_by`, they are scored by CV instead.

- A given rule `"R"` and `tau` (no choice) is fitted on the OOB
  components too, with the same fallbacks: below 10 families or without
  in-bag counts its coefficients are fitted in sample (with a message),
  and so they are for a deterministic forest. A given `"D"` needs no OOB
  components.

- Factor predictors are coded by the outcome of all training rows: their
  levels are ordered by the mean outcome (regression), by the share of
  the second class (binary classification) or by the first principal
  component of the class proportions (three or more classes;
  Coppersmith, Hong and Hosking 1999). So a row's own outcome enters the
  coding its OOB families split on: with factors the OOB error is
  somewhat optimistic.

**Cross-validation** (`select_by = "cv"`). For each fold, the forest is
refitted with the parameters of `object` on the other folds (after
`set.seed(seed)`, the same for every fold), its components, shares and
refit are computed on those rows, and every candidate is scored on the
held-out fold. The CV RMSE is `sqrt(sum(SSE_k) / n)`, summing the
squared errors over the K folds.

The chosen (or given) pair is then applied to `object` itself (not
refitted): its components and shares on `data`, and for `"R"` the least
squares on the OOB components of `data`'s rows (out of bag) or on its
components on `data` (by CV). `data` and `y` must therefore be the rows
`object` was trained on, in the order of the fit (checked).

**Purification.** `object` is purified in place at all orders, as
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
does, if it is not purified yet. This does not change its predictions. A
forest purified with a `maxp_interaction` below its order is purified
again at all orders, with a message: the pruning needs every component.
The out-of-bag path purifies `object` (again) in mode 2, which gives the
same components (a mode-1 purification's up to rounding).

**Random numbers.** The folds are drawn after `set.seed(seed)`, and the
CV's fold fits after it too. Afterwards the random number stream is
restored to its state on entry (after drawing `seed` if it is `NULL`),
so the call does not reset the caller's stream.

**Cost.** The out-of-bag choice purifies `object` once and evaluates
each family's components at its OOB rows; it needs no fit. The CV fits
and purifies the forest K more times. Both evaluate the components on
the training rows: only those the purified forest stores (at most
`sum(choose(p, seq_len(max_interaction)))`, usually far fewer: the
subsets of the trees' predictor sets), in one pass over the rows. The
memory is about `n` doubles per stored component, a few times that for
the out-of-bag choice. `nthreads` speeds up the fits and the
purification.

Only regression forests are supported.

## See also

[`predict.rpf_pruned()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).

## Examples

``` r
# A main effect and a pure interaction
set.seed(1)
x <- matrix(runif(600 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(600, sd = 0.3)
fit <- rpf(x, y, max_interaction = 2, ntrees = 20, splits = 50)

# Rule and threshold chosen out of bag (no refits): the largest tau within one SE
# of the lowest error
pruned <- prune_components(fit, x, y, seed = 1)
pruned
#> ── Pruned Random Planted Forest ────────────────────────────────────────────────
#> Rule R (drop + refit one coefficient per component, fitted out of bag), tau = 10% of Var(prediction)
#> Chosen by out-of-bag error (1-SE rule; R refitted in 5 folds, seed 1) among 16 candidates.
#> 
#> Nonzero components: 9
#> Kept:    2 (by order 1: 1, 2: 1)
#> Dropped: 7 (by order 1: 3, 2: 4)
#> Largest dropped: x1:x3 (share 1.8%); dropped shares sum to 4.4%
#> 
#> Intercept: -0.0162521
#> Kept components (largest first), share on the training rows and coefficient:
#>  component share beta 
#>  x1        61%   0.977
#>  x2:x3     19%   1.863
pruned$oob
#>    rule   tau  oob_rmse     se_mse oob_rmse_raw se_mse_raw n_kept   min one_se
#> 1     D 0e+00 0.5361833 0.01395762    0.5638097 0.01187499      9 FALSE  FALSE
#> 2     D 1e-04 0.5361833 0.01395762    0.5638097 0.01187499      9 FALSE  FALSE
#> 3     D 3e-04 0.5361833 0.01395762    0.5638097 0.01187499      9 FALSE  FALSE
#> 4     D 1e-03 0.5333645 0.01401006    0.5607863 0.01162696      7 FALSE  FALSE
#> 5     D 3e-03 0.5204500 0.01389558    0.5467446 0.01264863      5 FALSE  FALSE
#> 6     D 1e-02 0.5201940 0.01503833    0.5452333 0.01395539      4 FALSE  FALSE
#> 7     D 3e-02 0.5044697 0.02196168    0.5187196 0.02147003      2 FALSE  FALSE
#> 8     D 1e-01 0.5044697 0.02196168    0.5187196 0.02147003      2 FALSE  FALSE
#> 9     R 0e+00 0.3836019 0.02250883    0.4397841 0.02277146      9 FALSE  FALSE
#> 10    R 1e-04 0.3836019 0.02250883    0.4397841 0.02277146      9 FALSE  FALSE
#> 11    R 3e-04 0.3836019 0.02250883    0.4397841 0.02277146      9 FALSE  FALSE
#> 12    R 1e-03 0.3857449 0.02245586    0.4404874 0.02275764      7 FALSE  FALSE
#> 13    R 3e-03 0.3858999 0.02053559    0.4374335 0.01928150      5 FALSE  FALSE
#> 14    R 1e-02 0.3841273 0.01885531    0.4359182 0.01767504      4 FALSE  FALSE
#> 15    R 3e-02 0.3800078 0.01757644    0.4315554 0.01610271      2 FALSE  FALSE
#> 16    R 1e-01 0.3800078 0.01757644    0.4315554 0.01610271      2  TRUE   TRUE

# ... or the lowest error
prune_components(fit, x, y, select = "min", seed = 1)
#> ── Pruned Random Planted Forest ────────────────────────────────────────────────
#> Rule R (drop + refit one coefficient per component, fitted out of bag), tau = 10% of Var(prediction)
#> Chosen by out-of-bag error (lowest OOB RMSE; R refitted in 5 folds, seed 1) among 16 candidates.
#> 
#> Nonzero components: 9
#> Kept:    2 (by order 1: 1, 2: 1)
#> Dropped: 7 (by order 1: 3, 2: 4)
#> Largest dropped: x1:x3 (share 1.8%); dropped shares sum to 4.4%
#> 
#> Intercept: -0.0162521
#> Kept components (largest first), share on the training rows and coefficient:
#>  component share beta 
#>  x1        61%   0.977
#>  x2:x3     19%   1.863

# ... or by 5-fold CV
pruned_cv <- prune_components(fit, x, y, select_by = "cv", seed = 1)
pruned_cv$cv
#>    rule   tau   cv_rmse      se_mse n_kept   min one_se
#> 1     D 0e+00 0.4673140 0.026285812   10.0 FALSE  FALSE
#> 2     D 1e-04 0.4673225 0.026283390    9.8 FALSE  FALSE
#> 3     D 3e-04 0.4674804 0.026401513    9.0 FALSE  FALSE
#> 4     D 1e-03 0.4674112 0.025471341    7.6 FALSE  FALSE
#> 5     D 3e-03 0.4638035 0.024529257    5.8 FALSE  FALSE
#> 6     D 1e-02 0.4630849 0.024904421    3.0 FALSE  FALSE
#> 7     D 3e-02 0.4462390 0.018379566    2.0 FALSE  FALSE
#> 8     D 1e-01 0.4462390 0.018379566    2.0 FALSE  FALSE
#> 9     R 0e+00 0.3851965 0.012296486   10.0 FALSE  FALSE
#> 10    R 1e-04 0.3851968 0.012296395    9.8 FALSE  FALSE
#> 11    R 3e-04 0.3849120 0.012053138    9.0 FALSE  FALSE
#> 12    R 1e-03 0.3824749 0.010638693    7.6 FALSE  FALSE
#> 13    R 3e-03 0.3686001 0.009522577    5.8 FALSE  FALSE
#> 14    R 1e-02 0.3506606 0.007501968    3.0 FALSE  FALSE
#> 15    R 3e-02 0.3375723 0.005728376    2.0 FALSE  FALSE
#> 16    R 1e-01 0.3375723 0.005728376    2.0  TRUE   TRUE

# A given rule and threshold: no choice; R's coefficients are still fitted out of bag
prune_components(fit, x, y, rule = "R", tau = 0.01)
#> ── Pruned Random Planted Forest ────────────────────────────────────────────────
#> Rule R (drop + refit one coefficient per component, fitted out of bag), tau = 1% of Var(prediction)
#> Given, not chosen.
#> 
#> Nonzero components: 9
#> Kept:    4 (by order 1: 1, 2: 3)
#> Dropped: 5 (by order 1: 3, 2: 2)
#> Largest dropped: x2 (share 0.72%); dropped shares sum to 1.3%
#> 
#> Intercept: -0.0169604
#> Kept components (largest first), share on the training rows and coefficient:
#>  component share beta    
#>  x1        61%    0.97670
#>  x2:x3     19%    1.86236
#>  x1:x3     1.8%  -0.04332
#>  x1:x2     1.3%   0.06363

predict(pruned, x[1:5, ])
#> # A tibble: 5 × 1
#>    .pred
#>    <dbl>
#> 1 -0.747
#> 2  0.544
#> 3 -1.08 
#> 4  1.04 
#> 5  0.620
```
