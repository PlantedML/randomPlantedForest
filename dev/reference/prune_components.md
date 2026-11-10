# Prune the components of a Random Planted Forest

Drops the small components of a forest's functional decomposition and,
optionally, refits one coefficient per kept component. The rule and the
threshold are chosen out of bag by default, from the forest itself
without refits, or by K-fold cross-validation on the training rows. By
default every kept interaction of three or more predictors is then
coarsened so that it can be shown in at most `max_panels` heat maps, and
the model is pruned again.

## Usage

``` r
prune_components(
  object,
  data,
  y = NULL,
  ...,
  max_panels = 10,
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

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model (regression or classification).

- data:

  `[data.frame | matrix]`: The training data of `object`, in the order
  of the fit: the predictors as passed to
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  as `x`, or the `data` of the formula or recipe interface. For a
  formula or recipe fit, `data` may contain the outcome.

- y:

  `[numeric | factor | NULL: NULL]`: The outcome of the training rows
  (the factor for classification). Required unless `data` contains the
  outcome of a formula or recipe fit.

- ...:

  Reserved for future expansion. Unknown arguments are an error.

- max_panels:

  `[numeric(1): 10]`: The display: the most heat maps a kept component
  of three or more predictors may need; such a component is coarsened to
  at most `max_panels` heat maps, and the model is pruned again (see
  Details). A whole number of at least 1, or `Inf` for no display (the
  first pruning's model).

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

- `fit`: `object`, as it was (not purified by the pruning).

- `rule`, `tau`: the chosen rule and threshold.

- `intercept`: `alpha` (the forest's intercept `f_0` for `"D"`); for
  classification one per score column.

- `beta`: named vector of the coefficients of the kept components (all 1
  for `"D"`). For `"R"` with `select_by` `"oob"`, `intercept` and `beta`
  are fitted on the out-of-bag components (see Details).

- `sets`: the predictor sets of the kept components, as column indices.

- `components`: a `data.frame` of every nonzero component of `object` on
  the training rows (after a second pruning of the display, only of
  those the first pruning kept: see below): `component` (its name, as in
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)),
  `order`, `share`, `kept`, `beta` (`NA` if dropped), and `model_share`:
  the share of the final model's component,
  `Var(beta_S f_S) / Var(prediction of the model)` on the training rows
  (with the display, of `f_S` in its heat-map form; for classification,
  as `share`, over the centred score columns; `NA` if dropped). `share`
  is what `tau` thresholds; `model_share` is the component's importance
  in the model, which
  [print()](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_pruned.md)
  shows. As the components need not be uncorrelated, the shares need not
  sum to 1.

- `oob`: `NULL` unless the pair was chosen out of bag; otherwise a
  `data.frame` with one row per candidate: `rule`, `tau`, `oob_rmse`
  (corrected; used for the choice; 0 where the correction exceeds the
  uncorrected MSE, see Details), `se_mse` (from its per-fold MSEs),
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

- `criterion` (classification only): `"brier"` or `"logloss"`, what the
  candidates were scored by.

- With the display (finite `max_panels`) and a component coarsened, a
  second pruning runs, and these describe it: its candidates are the
  components the first pruning kept, with the coarsened ones in place,
  and `oob` is its table. `components` then lists only these components,
  and their `share` (which the second pruning's `tau` thresholds) is of
  the variance of their model (the intercept plus these components) on
  the training rows, not of the forest's prediction; every nonzero
  component of `object` with its share of the forest's prediction is in
  `display$first$components`. Without a coarsened component, they are
  the first pruning's. With the display, `components` also has `form`
  (`"full"` for one or two predictors, `"asis"` if its own cut points
  need at most `max_panels` heat maps, `"coarsened"`; `NA` for a
  component of three or more predictors the first pruning dropped, so
  only without a second pruning) and `heat_maps`. The element `display`
  is a list: `max_panels`; `first`, the first pruning's model (without
  `fit`); `coarse`, per coarsened component (named as the components):
  `S` (all the component's predictors), `C` (the coarse predictors,
  `|S| - 2` of `S`: the panel predictors), `q` (the numbers of intervals
  aimed at; merged ties can leave fewer, `lengths(cuts) + 1`), `cuts`
  (per coarse predictor), `cells` (the number of heat maps,
  `prod(lengths(cuts) + 1)`), `var_kept` (the share of the component's
  variance kept), `forest_panels`, `kept` (in the final model),
  `candidates` (every candidate's coarse set, `q`, cells and variance
  kept), `stored` (what prediction needs; no training data) and
  `seconds` (the time to choose and build its form); `asis`, per
  component kept as it is (form `"asis"`), its panels with the same
  names and meaning: `S` (all the component's predictors), `C` (the
  panel predictors: the `|S| - 2` of `S` with the fewest own cut points;
  of ties, the first in `S`), `cuts` (per predictor of `C`, its own cut
  points: the leaf bounds of the trees containing `S` inside the range
  of the training values), `forest_panels` (`prod(lengths(cuts) + 1)`)
  and `kept`; `forms`, per component of three or more predictors the
  first pruning kept, its form, heat maps and variance kept; `second`
  (whether a second pruning ran); `seconds` (the time of each step);
  `cache`, an environment holding the coarsened components' search
  indices, built at their first
  [`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
  (a cache, not part of the model: rebuilt when missing or invalid, as
  after saving and loading). In `coarse` and `asis` alike, a panel is a
  cell of `cuts` (a value `x` of predictor `k` of `C` lies in cell
  `findInterval(x, cuts[[k]])`), shown as a heat map over the two
  predictors `setdiff(S, C)` at full resolution.

- `store`: the pieces of the kept components that are not coarsened,
  from which the model predicts (see **Purification**), and `cache`, an
  environment holding them built for prediction.

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
works family by family. So every family's components (see
**Purification** below) are also evaluated at the training rows outside
its bootstrap sample. A row's out-of-bag (OOB) component is the mean
over its OOB families (about 37% of them), and its OOB intercept
likewise. The purification uses the predictors of all training rows, not
the outcome, so no outcome information of a row enters its OOB
components.

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

- If the correction exceeds a candidate's uncorrected MSE, its corrected
  MSE is reported as 0. That is not a perfect candidate but an
  overshooting, unreliable estimate (it happens with very few rows, e.g.
  `"R"` with about as many coefficients as rows, or with almost no
  noise). Such candidates are not chosen, unless all candidates are at 0
  (then they tie).

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

**The display** (`max_panels`, default 10; `Inf` switches it off). A
component of two predictors is one heat map. A kept component `f_S` of
three or more predictors is shown as heat maps over two of its
predictors, one per cell of the others. If its own cut points need more
than `max_panels` such heat maps, it is coarsened:

- Two predictors of `S` stay at full resolution. Each other predictor
  `k` is cut into `q_k` intervals of equal numbers of training rows
  (ties merged), with `prod(q) <= max_panels`; `q_k = 1` averages `k`
  out. Only the maximal `q` vectors (no `q_k` can grow) are considered;
  they are listed directly, not searched among all
  `max_panels^(|S| - 2)` vectors.

- Within a cell, `f_S` is averaged over the coarse predictors with the
  product of their training marginals. The coarsened component is not
  purified again:
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  identifies the components under the joint distribution of the training
  rows, not under the product of the marginals, and the averaging
  changes the component's means along its two full-resolution predictors
  and its mean over the training rows. It is used as it is (no part of
  it is moved into the intercept or into components of lower order).

- Of all choices of the two full predictors and the `q` vector, the one
  keeping the largest share of the variance of `f_S` on the training
  rows is taken (ties: the first in lexicographic order).

Then the same rule prunes again, on the components the first pruning
kept, with the coarsened ones in place of the originals. Its shares,
which `tau` thresholds, are of the variance of that model (the intercept
plus these components) on the training rows, not of the forest's
prediction. A pair chosen out of bag is chosen again out of bag, with
the same candidates and folds, scored with each family's coarsened
components at its out-of-bag rows (R's coefficients out of bag); a given
rule and `tau`, or a pair chosen by CV, is applied again as it is (R's
coefficients out of bag if the first pruning fitted them so, else in
sample; a CV choice is not repeated, as it would need the display in
every fold's refit). The prediction then uses the coarsened components,
so the pruning judges the model that is shown. In simulations (48 data
sets) this second pruning had 0.49% lower test error than pruning and
then coarsening, and 0.57% lower than not coarsening; on California
housing (50 families, order 3) the default keeps 30 components, shown in
75 heat maps. If no kept component needs coarsening, the model is the
first pruning's. The display costs about as much again as the out-of-bag
pruning: a purification restricted to the components being coarsened
(two at a time) for their candidates, and one for their out-of-bag
values; the stored coarsened components take a few MB each on 16,000
rows at order 3, more at higher orders.

**Purification.** The pruning needs every component's values at the
training rows (and out of bag, its OOB means), not the purified forest:
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
keeps every component of every tree family as pieces, which at high
orders takes far more memory than the forest (see **Cost**). So each
family is purified in turn, as
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
does (mode 2, all orders), its components are evaluated at the training
rows and added up, and the family's pieces are freed. The pruned model
then stores the pieces of its kept components only, merged over the
families exactly as
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
stores them, and
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
and
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
evaluate them from there; a coarsened component (the display) is stored
in its coarsened form. The components equal those of the purified forest
up to rounding (about `1e-15` relative), and the stored ones are the
purified forest's bit for bit. `object` itself is not purified: an
earlier
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
at any `maxp_interaction` or `mode`, is left as it is and does not
change the result.

**Random numbers.** The folds are drawn after `set.seed(seed)`, and the
CV's fold fits after it too. Afterwards the random number stream is
restored to its state on entry (after drawing `seed` if it is `NULL`),
so the call does not reset the caller's stream.

**Cost.** The out-of-bag choice purifies the families of `object` once
(see **Purification**); it needs no fit. The CV fits and purifies the
forest K more times (each fold's forest as
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
does). Only the components the forest has are evaluated (at most
`sum(choose(p, seq_len(max_interaction)))`, usually far fewer: the
subsets of the trees' predictor sets). The memory is the forest's plus
about `n` doubles per component, three times that for the out-of-bag
choice, plus one family's components at a time per thread; storing the
kept components passes over the families once more. On California
housing (16,512 rows, 200 families, `split_try = 3`, `t_try = 0.1`, the
fast draw, 8 threads; measured when this column pass was added) the
whole call of `rpf(splits = NULL, prune = TRUE)` peaked at 2.3 GB of
physical footprint (4.1 GB of resident memory, RSS) at order 3, and at
6.4 GB (9.9 GB RSS) at order 8 (all predictors), most of it the 31
coarsened components the display keeps; there
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
alone needs more than 12 GB of RSS. `nthreads` speeds up the fits and
the purification.

**Classification.** A classification forest is pruned the same way, with
its score columns: one for a binary outcome, one per class otherwise
(`K - 1` for `loss = "logit"`, whose first level is the reference).
Every component then has a column per score column, and:

- A component's share sums its variance over the score columns, centred
  across the classes (`"logit"`: with the reference level's score 0
  joined; `"exponential"`: the part all classes share does not count),
  over the same sum for the prediction. One kept list serves all
  classes.

- `"R"` fits one coefficient per kept component, shared by the classes,
  and one intercept per score column (`"exponential"` with more than two
  classes: the first is 0), out of bag as for regression.

- The criterion is the forest's loss's own: for `loss = "L1"` or `"L2"`
  the squared error of the raw scores summed over the columns (the Brier
  score; least squares for `"R"`), for `"logit"` and `"exponential"` the
  log-loss of the probabilities that `predict(type = "prob")` makes
  (maximum likelihood for `"R"`, by Newton's method from rule `"D"`'s
  model). The OOB correction is its second-order form (see
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  "Classification"). The tables then have `oob_loss`, `se_loss`,
  `oob_loss_raw` and `se_loss_raw` (the mean criterion) instead of the
  RMSE columns, and `cv_loss` and `se_loss` for CV; `select = "1se"`
  uses the SE of the mean criterion. A binary forest with `"L1"` or
  `"L2"` has one 0/1 column and is pruned as a regression forest of that
  response, with the RMSE columns.

- The display coarsens a component with one form (`C`, `q`, cuts) for
  all score columns, chosen by the variance kept of the centred columns;
  its stored form is one per score column.

- [`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
  takes `type` as
  [`predict.rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  does, and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  has one column per component and score column, named as for the
  forest.

## See also

[`predict.rpf_pruned()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md);
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
with `prune = TRUE` fits and prunes in one call.

## Examples

``` r
# A main effect and a pure interaction
set.seed(1)
x <- matrix(runif(600 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(600, sd = 0.3)
fit <- rpf(x, y, max_interaction = 2, ntrees = 20, splits = 50, prune = FALSE)

# Rule and threshold chosen out of bag (no refits): the largest tau within one SE
# of the lowest error
pruned <- prune_components(fit, x, y, seed = 1)
pruned
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 2 components.
#> Intercept: -0.00482281
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  x2:x3     58%       
#>  x1        44%       
#> 
#> How it was chosen:
#>   The forest: 20 tree families of 50 splits.
#>   From the forest's 10 nonzero components (purified): those with a share of at least 10% of the forest's Var(prediction), each with one coefficient, fitted out of bag (rule R).
#>   Rule and threshold chosen by out-of-bag error (1-SE rule; 5 folds, seed 1) among 16 candidates.
#>   Left out: 8 components, shares together 2.1% (largest x2, 0.65%).
pruned$oob
#>    rule   tau  oob_rmse      se_mse oob_rmse_raw  se_mse_raw n_kept   min
#> 1     D 0e+00 0.4326435 0.017132961    0.4577946 0.016471299     10 FALSE
#> 2     D 1e-04 0.4326435 0.017132961    0.4577946 0.016471299     10 FALSE
#> 3     D 3e-04 0.4326435 0.017132961    0.4577946 0.016471299     10 FALSE
#> 4     D 1e-03 0.4307297 0.016162161    0.4535609 0.015358376      7 FALSE
#> 5     D 3e-03 0.4274866 0.014897697    0.4476957 0.014815801      5 FALSE
#> 6     D 1e-02 0.4207083 0.015924505    0.4363396 0.016059946      2 FALSE
#> 7     D 3e-02 0.4207083 0.015924505    0.4363396 0.016059946      2 FALSE
#> 8     D 1e-01 0.4207083 0.015924505    0.4363396 0.016059946      2 FALSE
#> 9     R 0e+00 0.3346053 0.008475738    0.3756329 0.009224183     10 FALSE
#> 10    R 1e-04 0.3346053 0.008475738    0.3756329 0.009224183     10 FALSE
#> 11    R 3e-04 0.3346053 0.008475738    0.3756329 0.009224183     10 FALSE
#> 12    R 1e-03 0.3329692 0.009513342    0.3740199 0.010275340      7 FALSE
#> 13    R 3e-03 0.3347209 0.009565213    0.3737766 0.009465612      5 FALSE
#> 14    R 1e-02 0.3328582 0.009525521    0.3726909 0.009540695      2 FALSE
#> 15    R 3e-02 0.3328582 0.009525521    0.3726909 0.009540695      2 FALSE
#> 16    R 1e-01 0.3328582 0.009525521    0.3726909 0.009540695      2  TRUE
#>    one_se
#> 1   FALSE
#> 2   FALSE
#> 3   FALSE
#> 4   FALSE
#> 5   FALSE
#> 6   FALSE
#> 7   FALSE
#> 8   FALSE
#> 9   FALSE
#> 10  FALSE
#> 11  FALSE
#> 12  FALSE
#> 13  FALSE
#> 14  FALSE
#> 15  FALSE
#> 16   TRUE

# ... or the lowest error
prune_components(fit, x, y, select = "min", seed = 1)
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 2 components.
#> Intercept: -0.00482281
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  x2:x3     58%       
#>  x1        44%       
#> 
#> How it was chosen:
#>   The forest: 20 tree families of 50 splits.
#>   From the forest's 10 nonzero components (purified): those with a share of at least 10% of the forest's Var(prediction), each with one coefficient, fitted out of bag (rule R).
#>   Rule and threshold chosen by out-of-bag error (lowest OOB RMSE; 5 folds, seed 1) among 16 candidates.
#>   Left out: 8 components, shares together 2.1% (largest x2, 0.65%).

# ... or by 5-fold CV
pruned_cv <- prune_components(fit, x, y, select_by = "cv", seed = 1)
pruned_cv$cv
#>    rule   tau   cv_rmse      se_mse n_kept   min one_se
#> 1     D 0e+00 0.3675334 0.012220573   10.0 FALSE  FALSE
#> 2     D 1e-04 0.3675334 0.012220573   10.0 FALSE  FALSE
#> 3     D 3e-04 0.3675998 0.012257470    9.8 FALSE  FALSE
#> 4     D 1e-03 0.3666870 0.012325882    7.0 FALSE  FALSE
#> 5     D 3e-03 0.3620997 0.012530279    3.8 FALSE  FALSE
#> 6     D 1e-02 0.3578317 0.010789540    2.0 FALSE  FALSE
#> 7     D 3e-02 0.3578317 0.010789540    2.0 FALSE  FALSE
#> 8     D 1e-01 0.3578317 0.010789540    2.0 FALSE  FALSE
#> 9     R 0e+00 0.3563900 0.009680989   10.0 FALSE  FALSE
#> 10    R 1e-04 0.3563900 0.009680989   10.0 FALSE  FALSE
#> 11    R 3e-04 0.3562247 0.009595002    9.8 FALSE  FALSE
#> 12    R 1e-03 0.3516531 0.009322421    7.0 FALSE  FALSE
#> 13    R 3e-03 0.3354631 0.008332358    3.8 FALSE  FALSE
#> 14    R 1e-02 0.3286997 0.006150586    2.0 FALSE  FALSE
#> 15    R 3e-02 0.3286997 0.006150586    2.0 FALSE  FALSE
#> 16    R 1e-01 0.3286997 0.006150586    2.0  TRUE   TRUE

# A given rule and threshold: no choice; R's coefficients are still fitted out of bag
prune_components(fit, x, y, rule = "R", tau = 0.01)
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 2 components.
#> Intercept: -0.00482281
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  x2:x3     58%       
#>  x1        44%       
#> 
#> How it was chosen:
#>   The forest: 20 tree families of 50 splits.
#>   From the forest's 10 nonzero components (purified): those with a share of at least 1% of the forest's Var(prediction), each with one coefficient, fitted out of bag (rule R).
#>   Rule and threshold given.
#>   Left out: 8 components, shares together 2.1% (largest x2, 0.65%).

predict(pruned, x[1:5, ])
#> # A tibble: 5 × 1
#>    .pred
#>    <dbl>
#> 1 -0.768
#> 2  0.563
#> 3 -1.18 
#> 4  1.07 
#> 5  0.101

# A three-way interaction, shown in at most 4 heat maps
y3 <- x[, 1] + 4 * x[, 2] * x[, 3] * x[, 4] + rnorm(600, sd = 0.3)
fit3 <- rpf(x, y3, max_interaction = 3, ntrees = 20, splits = 100, prune = FALSE)
shown <- prune_components(fit3, x, y3, max_panels = 4, seed = 1)
shown
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 2 components.
#> Intercept: 0.0509343
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  x2:x3:x4  64%       
#>  x1        38%       
#> 
#> How it was chosen:
#>   The forest: 20 tree families of 100 splits.
#>   From the forest's 14 nonzero components (purified): those with a share of at least 10% of the forest's Var(prediction), each with one coefficient, fitted out of bag (rule R).
#>   Rule and threshold chosen by out-of-bag error (1-SE rule; 5 folds, seed 1) among 16 candidates.
#>   Left out: 12 components, shares together 11% (largest x2, 2%).
#>   Components of three or more predictors in at most 4 heat maps each, then rule and threshold chosen again on these (R, tau = 10%; none left out).
#>   Heat-map forms: x2:x3:x4: x4 in 4 intervals (95% of its variance)
shown$display$forms
#>   component order      form forest_panels heat_maps var_kept
#> 1  x2:x3:x4     3 coarsened           145         4 0.948961
```
