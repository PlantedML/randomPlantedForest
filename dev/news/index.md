# Changelog

## randomPlantedForest (development version)

### Breaking changes

- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  has new defaults: a plain `rpf(y ~ ., data)` now returns the pruned
  model (class `rpf_pruned`), for regression and classification. The new
  defaults, each switched by its argument:

  - `n_splits = NULL`: the number of splits is chosen out of bag while
    the forest grows (was 30).
  - `prune = NULL`: the forest is pruned by
    [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
    with its defaults, and every kept component of three or more
    predictors is shown in at most `max_panels = 10` heat maps (was no
    pruning). With `deterministic = TRUE` it means `prune = FALSE` and
    `n_splits = NULL` means 30, as there are no out-of-bag rows.
  - `n_try_split_points = NULL, share_candidates = NULL`: 3 thresholds
    per candidate and 10% of the possible splits per round for
    `split_structure = "leaf"` (was 10 / 0.4, still the default of the
    other split structures). This less greedy search is usually more
    accurate.
  - The split candidates of `"leaf"` are drawn by the fast draw
    (`options(randomPlantedForest.exact_draw = FALSE)`, now the
    default): the same distribution as the exact draw, with other random
    numbers, so the same seed gives another forest. It is faster on
    larger data.
  - `max_interaction = 3` (was 2). Order 3 is often more accurate than
    order 2; no limit (order p) is often much slower, and usually not
    more accurate.
  - `n_families = 200` (was 50): more families usually lower the test
    error and give sparser pruned models. The loss is `"L2"` for
    regression and, as below, now `"logit"` for classification. So every
    default fit changes, regression and classification, and
    [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
    returns an `rpf_pruned` object instead of an `rpf` unless
    `prune = FALSE`. To fit a regression forest as before, give the
    previous defaults and the exact draw:

  ``` r

  options(randomPlantedForest.exact_draw = TRUE)
  rpf(..., max_interaction = 2, n_families = 50, n_splits = 30, n_try_split_points = 10, share_candidates = 0.4, prune = FALSE)
  ```

  Classification fits of earlier versions cannot be reproduced this way:
  they change for every split structure and loss (see Bug fixes:
  classification uses regression’s split structures, the largest value
  of a predictor can be split off, and three or more classes with the
  exponential loss use SAMME). The default call does more work than
  before, so it takes more time and memory.

- The classification losses `"logit"` and `"exponential"` (binary, and
  SAMME for three or more classes) no longer clamp class shares, and
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)’s
  arguments `delta` and `epsilon` are removed: a call that gives either
  is an error (an unused argument). Each new leaf’s step now minimises
  the leaf’s penalised loss: the loss of its rows, plus that of one
  pseudo-row, split evenly over the classes, at the leaf’s current mean
  score
  ([`?rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  “Classification losses”). A split is scored by its gain over the
  parent leaf’s own step, and each new leaf pays for its pseudo-row. So
  a leaf of one class converges to a finite level (`log(2 n + 1)` for
  `n` rows and two classes) and is not moved again, instead of being
  split and pushed further in every round; splitting such a leaf, or
  repeating a split whose leaves are at their optimum, gains nothing and
  is not made; and a round in which no candidate gains makes no split.
  Previously the split criterion credited a split with the parent leaf’s
  own step at a class share clamped to `[delta, 1 - delta]`, while the
  update clamped it to `[epsilon, 1 - epsilon]`, so pure and nearly pure
  leaves won over splits that separate classes, and logit’s leaves reset
  the mean score to the log-odds of the raw share, which is not a step
  of the loss and left fits underconfident. With the defaults
  (`n_splits` out of bag, pruning) the test log-loss usually changes
  little. At fixed split counts the exponential loss no longer overfits
  as it did; on noisy data with many splits both losses can fit the
  noise more than the clamped logit loss did, which the out-of-bag
  choice of `n_splits` avoids.
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  gains `exact_step` (default `TRUE`): logit’s two new leaves get the
  exact minimisers of their penalised loss (Newton’s method to
  convergence on the leaf’s rows); `FALSE` takes one Newton step. The
  exponential loss’s step is exact by itself, so `exact_step` does not
  change its fits. The split search scores a candidate from sums over
  its two sides built as regression builds them (along the leaf’s sorted
  rows, or per bin), so classification fits are faster. Classification
  fits with `"logit"` or `"exponential"` change; `"L1"` and `"L2"` fits
  and regression are unchanged, and a model saved before still loads,
  prints and predicts.

- The default classification loss is now `"logit"` (`loss = NULL`), with
  two classes and with more; it was `"exponential"` (with more than two
  classes SAMME). So every default classification fit changes;
  `loss = "exponential"` gives the previous default, which stays
  available and unchanged. The test log-loss of `"logit"` is usually as
  low as that of `"exponential"` or lower, more so with three or more
  classes. With many classes `"exponential"` can need many more splits.

- A pruned model no longer holds the forest it was pruned from:
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  with pruning (the default) and
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  return the model without it, and their new argument
  `keep_forest = TRUE` keeps it. The model predicts, prints and gives
  its components from what it stores, bit for bit as with the forest, so
  [`predict()`](https://rdrr.io/r/stats/predict.html) (every `type`),
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
  [`print()`](https://rdrr.io/r/base/print.html), the display,
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  and `bundle()` work as before; the fit itself is unchanged. What is
  lost by default is the forest in `$fit`: pruning it again
  (`prune_components(model$fit, data, ...)`),
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
  and the forest’s own
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md)
  are an error that names `keep_forest = TRUE`. `$fit` still describes
  the forest (`$params`, `$oob_path`, `$splits_oob`, the blueprint), and
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)’
  caller still has the forest it passed. To get the forest back, fit
  with `keep_forest = TRUE`, or with `prune = FALSE` and prune with
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md).
  With `prune = FALSE` the argument has no effect. The forest’s memory
  is free when
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  returns, and the model is a fraction of the forest’s size; the model
  without its forest is plain R data, so
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) give a model that
  predicts, and
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  is quick and gives a blob of the model’s size. Models saved with their
  forest by an earlier version
  ([`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md))
  still load, print and predict, and keep their forest.

- With three or more classes, `loss = "logit"` now reports one score per
  class, summing to 0 across the classes, as `"exponential"` does:
  `predict(type = "link")` (and `"numeric"`) has a column
  `.pred_<level>` for every level,
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  a column `<component>__class:<level>` for every component and level,
  and one intercept and one remainder column per level (`target_levels`
  lists all levels), and
  [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md)
  a column `value_<level>` for every level; the same for pruned models,
  whose [`print()`](https://rdrr.io/r/base/print.html) shows one
  intercept per level. Before, the first level was the reference class
  and had no column (`K - 1` columns for `K` levels). The new scores are
  the old ones with the first level’s score 0 joined and centred across
  the `K` classes, so a class’s score is now its log-odds against the
  geometric mean of all classes’ probabilities, not against the first
  class. The fit, its probabilities and classes, the out-of-bag choices
  and the pruned model (which components are kept, and their
  coefficients) do not change, and a model saved before loads and
  predicts the same probabilities; only what is reported changes.
  `plot(class = NULL)` now draws the first level, which has values of
  its own. To get the old scores back, subtract the first level’s column
  from every column and drop it. Two classes are unchanged: one score,
  the log-odds of the second level.

- A pruned model of three or more classes with `loss = "exponential"`
  and rule `"R"` now reports its intercepts centred across the classes,
  so that its scores sum to 0 as the forest’s do:
  `predict(type = "link")`, the `intercept` of
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  and [`print()`](https://rdrr.io/r/base/print.html). Rule `"R"` fits
  them with the first class’s at 0, so its scores did not sum to 0 (by
  up to 3.7 on a simulated three-class set); the probabilities, classes
  and the stored `intercept` do not change.

### Renamed arguments

- Four arguments of
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  and the split structures have clearer names; the prefix gives the unit
  (`n_` a count, `share_` a share in (0, 1\], `max_` a cap):
  - `ntrees` → `n_families` (the number of tree families),
  - `splits` → `n_splits` (splits per tree family),
  - `split_try` → `n_try_split_points` (split points tried per
    candidate),
  - `t_try` → `share_candidates` (share of the possible splits drawn as
    candidates each round);
  - `split_structure`: `"leaves"` → `"leaf"`, `"hist"` →
    `"leaf_binned"`, `"cur_trees_1"` → `"tree_then_leaf"`,
    `"cur_trees_2"` → `"tree_then_all_leaves"`, `"res_trees"` →
    `"target_tree_then_leaf"` (the candidate, then how the leaf to split
    is found).
- The old names still work and give the same fit, with a deprecation
  warning once per session for each old name (class
  `lifecycle_warning_deprecated`;
  `options(lifecycle_verbosity = "warning")` warns every time, `"quiet"`
  never); giving an argument under both names is an error. The same
  holds for `oob_path_checkpoints(splits = )` and for a
  [`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
  ladder with the columns `split_try` and `t_try`.
- The new names are used in
  [`print()`](https://rdrr.io/r/base/print.html), in the tables of the
  new features (the `n_splits` column of `$oob_path`, the tuning’s
  `$tune$search` and `$tune$settings`,
  [`greediness_ladder()`](https://plantedml.com/randomPlantedForest/dev/reference/greediness_ladder.md)),
  and in errors and documentation. A fitted model’s `$params` keeps the
  old names (`ntrees`, `splits`, `split_try`, `t_try` and the old
  `split_structure` values), so code that reads them, and models saved
  by 0.5.0, keep working.

### New features

- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  gains `keep_level_order`, to keep the level order of factor predictors
  instead of ordering their levels by the outcome. `TRUE` keeps it for
  every factor, ordered or not, and a character vector for the named
  factors only; splits then separate neighbouring levels in the order of
  [`levels()`](https://rdrr.io/r/base/levels.html). The default,
  `FALSE`, orders every factor by the outcome as before, and fits are
  unchanged. The coding is stored in `$factor_levels`, which
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  and
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  use, and [`print()`](https://rdrr.io/r/base/print.html) names the
  factors that keep their order. The details of
  [`?rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  describe how factor predictors are coded.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  trees of one tree family to illustrate how a Random Planted Forest
  predicts. `type = "tree"` draws a split diagram in the style of Figure
  1b of the paper, reconstructed from the fitted leaves,
  `type = "boxes"` draws the leaves in variable space, before or after
  purification.
  [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md)
  returns the leaves as a data frame. Requires `ggplot2`, and
  `patchwork` for `type = "boxes"`. The “Getting started” vignette uses
  it to show how a tree family is built and what purification does.
- [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  prunes the purified components of a forest (regression or
  classification; see the classification entry below). A component is
  kept if its share of the variance of the prediction on the training
  rows is at least `tau`. Rule `"D"` drops the others; rule `"R"` also
  refits the intercept and one coefficient per kept component by least
  squares, which changes their scale but not their shape. The rule and
  `tau` are chosen out of bag by default (`select_by = "oob"`), without
  refitting: one purification evaluates every tree family’s components
  at the training rows outside its bootstrap sample, and each candidate
  is scored by its out-of-bag RMSE, corrected for the extra variance of
  averaging fewer families. The out-of-bag choice usually picks the
  cross-validated model, or one about as good, in a fraction of the
  time. `select_by = "cv"` chooses by K-fold cross-validation instead,
  refitting the forest on each fold; by default it is used below 10 tree
  families and for deterministic forests. By default the largest `tau`
  whose error is within one standard error of the lowest is chosen
  (`select = "1se"`; `select = "min"` takes the lowest). With
  `select_by = "oob"`, the coefficients of rule `"R"` are fitted on the
  out-of-bag components, to score each candidate, for the final model
  and for a given rule and `tau`, because the forest’s own components
  have also fitted the training rows’ noise; with `select_by = "cv"`
  they are fitted in sample. This sparse default is usually more
  accurate than the lowest out-of-bag error with in-sample coefficients,
  with far fewer components; it can drop a small true component. By
  default every kept component of three or more predictors is then shown
  in at most `max_panels = 10` heat maps, and the model is pruned again
  with the shown components (`max_panels = Inf` switches this off). Two
  predictors of such a component stay at full resolution; the others are
  cut into intervals of equal numbers of training rows whose counts
  multiply to at most `max_panels` (a count of 1 averages a predictor
  out), the component is averaged within the cells (weighted by the
  product of the coarse predictors’ training marginals; the result is
  used as it is, not purified again), and the cuts that keep most of its
  variance are taken. The same rule then prunes again, so that the
  pruning judges the model that is shown; predictions use the coarsened
  components, and `tau` is then a share of that model’s variance. This
  is usually a little more accurate than coarsening after pruning, or
  than not coarsening. Every fit records how often each training row is
  in each family’s bootstrap sample (fits are unchanged), and
  `rpf_marshal(include_data = TRUE)` stores these in-bag counts. The
  pruned model is often both sparser and more accurate than the unpruned
  forest. The result, of class `rpf_pruned`, has methods for
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  (now an S3 generic) and
  [`print()`](https://rdrr.io/r/base/print.html), and works with
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md),
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  and `bundle()`.
- `rpf(oob_path = TRUE)` gives the out-of-bag error along the splits
  from one fit, in `$oob_path` (for classification the loss’s criterion,
  see below). With the same seed, the forest with `k` splits is exactly
  the first `k` splits of a forest with more, so while each tree family
  grows, its prediction at the training rows outside its bootstrap
  sample is updated split by split. The table has the plain and the
  corrected out-of-bag RMSE (the correction of
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md))
  at four checkpoints per doubling of the splits
  ([`oob_path_checkpoints()`](https://plantedml.com/randomPlantedForest/dev/reference/oob_path_checkpoints.md)),
  or at given ones. Choosing `n_splits` by the lowest corrected error
  among the values of a cross-validation grid usually does as well as
  5-fold cross-validation over that grid. The fitted forest is the same
  as without `oob_path`, which is off by default.
- `rpf(n_splits = NULL)` chooses the number of splits out of bag while
  the forest grows (regression and classification). The forest’s
  corrected out-of-bag error is followed at four checkpoints per
  doubling of the splits, growing stops once its lowest value is one
  doubling back (at most 25,600 splits), and the forest at the lowest
  error is returned: the forest that
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  with that many `n_splits` and the same seed gives. This usually
  chooses the same `n_splits` as the lowest error along a long fixed
  budget of splits, while growing fewer. The tree families grow one at a
  time, `nthreads` together, and keep only their out-of-bag predictions
  and a 40-byte record per split, from which the chosen forest is
  rebuilt without a refit;
  `options(randomPlantedForest.adaptive_pause_gb)` (0.25 GB) bounds the
  memory of families waiting for the stop.
  [`print()`](https://rdrr.io/r/base/print.html) shows the chosen
  `n_splits` and where growing stopped. It is the new default (see
  Breaking changes). A corrected out-of-bag error at or below 0 at a
  checkpoint (a correction that overshot; never observed) is not chosen
  unless all are. A family waiting for the stop pauses at most one
  doubling past the first families’ estimate and grows on if the stop
  comes later, so a flat out-of-bag curve no longer makes every family
  grow far past the stop, with the same forest. Growing also stops once
  the error is on a plateau: once its lowest value has improved over the
  last doubling by at most one standard error of that improvement
  (paired over the training rows;
  `options(randomPlantedForest.adaptive_tolerance)`, 1 by default, 0 for
  the stop on a new minimum alone). The choice is still the lowest error
  among the splits grown, so the plateau stop changes the forest only
  where the lowest error would have come later. It usually chooses the
  same `n_splits` as the stop on a new minimum alone, and grows fewer
  splits. The families grow only to about the plateau stop, so the fit
  often takes less time; where the families waiting for the stop do not
  fit in the pause budget, they grow on to the estimate of the stop on a
  new minimum, as a regrowth from their seed would cost more.
- `rpf(prune = TRUE)` returns the pruned model instead of the forest
  (regression and classification): the forest is pruned by
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  with its defaults on the training rows, including the display at
  `max_panels` (a new argument of
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  10 by default, `Inf` for no coarsening). The result equals
  `prune_components(rpf(..., prune = FALSE), data)` bit for bit from the
  same random number state, for the formula, x/y and recipe interfaces.
  With `n_splits = NULL`, one call chooses the number of splits out of
  bag, prunes and returns the model to show. `prune = TRUE` is an error
  with `deterministic = TRUE`, and `purify` is not used with it. It is
  the new default (`prune = NULL`, see Breaking changes);
  `prune = FALSE` returns the forest as before.
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  gains `min_leaf_rows`, the minimum number of rows on each side of
  every split, so in every leaf. Every split structure takes it, in
  regression and classification; the default 1 imposes no limit. Split
  points are tried only among those that keep this many rows on each
  side, and leaves without one are no longer drawn: as split candidates
  with `"leaf"` and `"leaf_binned"`, by the leaf draw of
  `"tree_then_leaf"` and `"target_tree_then_leaf"` (a leaf’s weight is
  its number of such split points), and `"tree_then_all_leaves"` skips
  them. Binary predictors are split whenever both sides keep enough
  rows. The default `min_leaf_rows = 1` fits the same models as before,
  bit for bit.
- Classification with `split_structure = "leaf"` (the default) takes the
  fast draw of split candidates
  (`options(randomPlantedForest.exact_draw = FALSE)`) as regression
  does, and `min_leaf_rows` with every split structure. Both use
  regression’s draw: the same candidates are left out, and a
  classification fit with `loss = "L2"` on a 0/1 response is the
  regression fit, bit for bit, with either option. With
  `min_leaf_rows = 1` and the exact draw (the default before this
  version) classification fits are unchanged, bit for bit.
- `rpf(tune = TRUE)` tunes out of bag in one call (regression and
  classification): the greediness of the split search
  (`n_try_split_points` and `share_candidates`, on the ladder of the new
  [`greediness_ladder()`](https://plantedml.com/randomPlantedForest/dev/reference/greediness_ladder.md))
  and `min_leaf_rows`. Every forest of the search has `n_splits` chosen
  out of bag and is pruned by
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  with its defaults, and the forest whose pruned model has the lowest
  corrected out-of-bag error wins; the display (`max_panels`) is applied
  once, to the chosen model. The search is coordinate-wise, about 10
  forests: the ladder from 10 / 0.4 to 1 / 0.0125, then `min_leaf_rows`
  1, 5, 20 and 100, then the two ladder neighbours of the chosen point.
  `tune = "quick"` grows two forests, 3 / 0.1 and 1 / 0.025, and the new
  [`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
  sets your own search. The result is the chosen forest’s pruned model,
  with the search table in `$tune`. The full search often lowers the
  test error clearly against 10 / 0.4, and does as well as 5-fold
  cross-validation over the same settings, or better, at a small
  fraction of its cost; the quick preset gets most of the gain. The
  default stays `tune = FALSE`.
- Classification forests go through the same pipeline as regression
  forests: `rpf(oob_path = )`, `rpf(n_splits = NULL)`,
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  with its display (`max_panels`), `rpf(prune = TRUE)` and
  `rpf(tune = )`. Every out-of-bag choice is scored by the loss’s own
  criterion, corrected for the out-of-bag family count: for
  `loss = "L1"` and `"L2"` the squared error of the raw scores summed
  over the classes (the Brier score), for `"logit"` and `"exponential"`
  the log-loss of the probabilities that `predict(type = "prob")`
  returns, with the correction in its second-order form. The corrected
  out-of-bag log-loss is usually close to the test log-loss, while the
  uncorrected one is above it. A component is kept or dropped for all
  classes, sized by its variance summed over the score columns centred
  across the classes; rule `"R"` fits one coefficient per kept
  component, shared by the classes, and one intercept per score column
  (least squares for the Brier score, maximum likelihood for the
  log-loss), out of bag; the display chooses one panel form per
  component for all classes. The pruned model’s
  [`predict()`](https://rdrr.io/r/stats/predict.html) takes `type` as
  [`predict.rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  does, and its
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  has one column per component and class. With two classes and
  `loss = "L2"`, every choice is that of the regression pipeline on the
  0/1 response, bit for bit. With the new defaults (see Breaking
  changes) a classification fit takes all these steps.
- `options(randomPlantedForest.exact_draw = FALSE)` draws the split
  candidates of regression with `split_structure = "leaf"` from a tree
  of weights kept up to date between rounds, instead of one random
  number per possible split in every round
  ([`?rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  “The draw of candidates for `leaves`”). This fits faster on larger
  data. It draws from the same distribution and leaves out the same
  candidates under `min_leaf_rows`, but with other random numbers, so
  the same seed gives another forest. `TRUE` (the default before this
  version) fits the same models as before, bit for bit.
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  gains `interaction_constraints`, a list of feature sets that restricts
  which interactions the forest may form, as in xgboost: a tree on the
  predictors `T` is grown only if `T` is a subset of one of the sets.
  The sets are given by predictor names or column indices (with a
  formula, indices follow the formula’s predictors; with a recipe, names
  only). They apply together with `max_interaction`, for every
  `split_structure`, in regression and classification, and
  [`print()`](https://rdrr.io/r/base/print.html) shows them. The default
  `NULL` leaves fits unchanged.
- [`print()`](https://rdrr.io/r/base/print.html) of a pruned model shows
  the final model: its intercept and its components, largest first, each
  with its importance, its share of the variance of the model’s
  prediction on the training rows (`$components$model_share`, a new
  column: the coefficient times the component, in the form the model
  uses, over the model’s prediction; for classification over the centred
  score columns). `share`, what `tau` thresholds, is the forest’s
  component before the coefficient, so a component with a large
  coefficient is a larger part of the model than its `share`. How the
  components were chosen (the forest’s families and splits, given or
  chosen out of bag, the rule and threshold, what was left out, the
  heat-map forms) follows under “How it was chosen”. No choice or fit
  changes.

### Bug fixes

- `split_structure = "tree_then_leaf"` samples leaves in proportion to
  their available thresholds, as documented. After a split, the sampling
  weights of new leaves counted rows instead of distinct values, and the
  cached prefix sums were built incorrectly. Leaves with tied values
  were then drawn too often, and draws were wasted on leaves that cannot
  be split. With `deterministic = TRUE`, `"tree_then_leaf"` evaluated
  only the first split candidate in every round; it now evaluates the
  first `max_candidates`, as the other split structures do. Fitted
  models change for regression with
  `split_structure = "tree_then_leaf"`, random and deterministic.
- `split_structure = "target_tree_then_leaf"` with
  `deterministic = TRUE` evaluates thresholds evenly spread over all
  (predictor, leaf, threshold) triples of a candidate. Previously it
  tried only the lowest thresholds of each leaf and left most of the
  variance unexplained. Fitted models change for regression with
  `split_structure = "target_tree_then_leaf"` and `deterministic = TRUE`
  only.
- `split_structure = "tree_then_all_leaves"` with `deterministic = TRUE`
  evaluates `n_try_split_points` thresholds evenly spread over each
  leaf, as `"tree_then_leaf"` does. Previously it tried the 1st to 9th
  smallest values of each leaf. Fitted models change for regression with
  `split_structure = "tree_then_all_leaves"` and `deterministic = TRUE`
  only.
- Classification’s split structures were never completed to match
  regression’s. They now use the same candidate pool, the same sampling
  of candidates, leaves and thresholds, and the same aging; only the
  loss differs. With `loss = "L2"` a classification fit on a 0/1
  response now equals the regression fit (same seed). Previously
  `"tree_then_leaf"` was `"tree_then_all_leaves"`;
  `"target_tree_then_leaf"` used another candidate pool and ignored
  `split_decay_rate`; `"leaf_binned"` only ever split the root, so it
  fitted an additive model of stumps whatever `max_interaction` was; and
  `"tree_then_all_leaves"` never split the candidate tree’s own leaves.
  The fixes above for `"tree_then_leaf"` sampling and for
  `deterministic = TRUE` apply to classification too. Classification
  fits change for every `split_structure` and every loss.
- The internal loss `"exponential_2"` no longer crashes.
- Splits can separate the largest value of a predictor within a leaf.
  Previously every split structure except `"leaf"` for regression
  excluded that split, so binary predictors were never split on, and
  classification could not learn from them at all. Fitted models change
  for classification and for regression with `split_structure` other
  than `"leaf"`.
- A tree family that makes no split predicts the mean outcome of its
  bootstrap sample, and for classification its rows’ fit (with
  `loss = "logit"` or `"exponential"` the class shares smoothed by the
  losses’ pseudo-row, see Breaking changes). Previously it predicted 0,
  a probability of 0.5: constant predictors, a single row or
  `min_leaf_rows` above half the rows gave an all-zero model, and
  families that could not split pulled the forest towards 0.
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  puts that value into the intercept. Fitted models change only where a
  family made no split.
- The same holds when `n_splits = NULL` (the default) chooses the number
  of splits out of bag: a family without a split up to the chosen number
  predicts its rows’ fit instead of 0, as a fit with that `n_splits`
  does. The pruned model already got the level from its intercept; the
  forest in it (`$fit`) did not. Without `n_splits = NULL`, and for
  families with a split, fits are unchanged. The pruned model also no
  longer fails for classification when no component has any variance
  (for example, constant predictors).
- Errors while fitting or purifying with `nthreads > 1` are reported as
  R errors. Previously an error in a worker thread terminated the R
  session.
- `predict(type = "link")` works for classification forests with
  `loss = "L1"` or `"L2"`. It returns the raw predictions, as
  `type = "numeric"` does (the link is the identity). Previously it
  stopped with an error about the size of `pred`.
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  stops with a clear error on inputs that gave wrong fits without a
  warning or failed later: an infinite outcome, `t_try = 0`, predictor
  names that are missing, duplicated or contain `:`, and interaction
  terms in the formula.
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  finds interactions itself, so write `y ~ x1 + x2` instead of
  `y ~ x1 * x2`; to fit a product as a predictor of its own, write
  `I(x1 * x2)` or add it to the data as a column. A character or logical
  outcome and a factor outcome with a single class get the package’s own
  message instead of one from `hardhat` or
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). Unused
  levels of a factor outcome are dropped with a warning; previously a
  binary outcome with an unused level was fitted, but
  [`predict()`](https://rdrr.io/r/stats/predict.html) failed.
  [`predict()`](https://rdrr.io/r/stats/predict.html) on zero rows
  returns zero rows instead of an error.
- The C++ code compiles without warnings under `-Wall -pedantic`.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  returns one intercept per target level for multiclass classification,
  named by `target_levels`
  ([\#70](https://github.com/PlantedML/randomPlantedForest/issues/70)).
  Previously only the first level’s intercept was returned, so
  components plus intercept did not sum to the prediction for other
  levels.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  accepts a matrix as `new_data`, like
  [`predict()`](https://rdrr.io/r/stats/predict.html) does. It no longer
  turns the caller’s data.frame into a data.table by reference.
- [`predict()`](https://rdrr.io/r/stats/predict.html) after
  `purify(maxp_interaction = k)` with `k` below the fit’s
  `max_interaction` returns the forest’s prediction. Previously it
  returned the intercept plus the components up to order `k`, leaving
  out the higher orders.
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  now returns their sum as `remainder`, so that its result sums to the
  prediction, and
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  keeps the purification’s cap.
- `loss = "exponential"` is with more than two classes the multiclass
  exponential loss SAMME (Zhu, Zou, Rosset and Hastie, 2009): the split
  criterion, the leaf step and the weights all come from the one loss
  `exp(-y.f / K)`. Each row has one weight, `exp(-g)` with `g` its
  class’s score, and a leaf moves the scores by
  `log(V_k) - mean_j log(V_j)`, with `V` its class weight sums plus
  those of its pseudo-row (see the classification losses under Breaking
  changes); [`predict()`](https://rdrr.io/r/stats/predict.html) takes
  the softmax of the scores as before. Previously each class had its own
  one-vs-rest weights, and its step was `log(p_k / p_0)` with
  `p_0 = 1 - sum_k p_k` mixing shares under the different weights; it is
  0 at the start and was clamped to `epsilon`, so every class’s score
  rose by a common `log(1 / epsilon)` that the softmax ignores but the
  weights did not. The weights then no longer matched the probabilities
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives. With
  `n_splits` chosen out of bag (`n_splits = NULL`) and the forest pruned
  (`prune = TRUE`), the pruned model’s test log-loss is now about that
  of logit, often with far fewer splits chosen and a faster fit. Fitted
  models change for `loss = "exponential"` with more than two classes
  only; binary fits are unchanged.
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  with `mode = 2` no longer drops training rows that lie on a leaf’s
  lower bound. Purified components and the intercept now equal those of
  `mode = 1`.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  returns a `remainder` whenever `m` leaves out components of the
  forest, also when `predictors` leaves out some of the fit’s predictors
  at the fit’s `max_interaction`. Previously the intercept and `m` then
  did not sum to the prediction.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  takes the components’ predictors from the forest’s own columns, in
  their order. A formula with terms such as `log(x1)`, `I(x1^2)` or
  `. - x3` no longer makes it fail, and with factor predictors the
  components come in the order of the fit’s columns, which has the
  factors last.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  on a `new_data` with zero rows returns zero rows with the right
  columns, instead of failing with “Feature vector is empty.”.
- [`print()`](https://rdrr.io/r/base/print.html) counts the forest’s own
  predictor columns: with formula terms such as `log(x1)` or `. - x3` it
  showed the number of raw columns the formula reads.
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  on a forest purified with `maxp_interaction = k` purifies it again
  when asked for more orders (all orders by default). Previously it
  returned the forest unchanged, so the orders above `k` stayed zero in
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md).
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  also checks its arguments when the forest is already purified. A
  capped forest restored by
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  without its training data cannot be purified again; asking for more
  orders is an error that says so.

### Documentation

- [`?purify`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  states how the components are identified (the marginal identification
  over the training rows, Hiabu, Meyer and Wright, 2023), and that the
  intercept is the mean prediction on the training rows.
- [`?rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  says that `deterministic = TRUE` also ignores `split_decay_rate`, how
  `split_structure = "leaves"` draws its thresholds, and that a seed
  reproduces a fit on the same platform only. `splits` counts rounds:
  near saturation a round can make no split. Ties between equally good
  candidates may go to different candidates under the exact and the fast
  draw, so their fits differ in distribution near saturation, though not
  in accuracy.
- [`?print.rpf_forest`](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_forest.md)’s
  example prints a forest (`export_forest = TRUE`), and
  [`?plot.rpf`](https://plantedml.com/randomPlantedForest/dev/reference/plot.rpf.md)
  cites the RPF paper as 2020.

### Performance

- Classification with `loss = "logit"` fits much faster, with bitwise
  the same fits. For every candidate split, the logit loss copied the
  weights of all observations twice, although it reads only those of the
  leaf being split, so a candidate cost time in proportion to the sample
  size instead of the leaf size.
- Fitting regression forests with `split_structure = "leaf"` (the
  default) is faster on larger data. Fitted models are bitwise
  unchanged.
- [`predict()`](https://rdrr.io/r/stats/predict.html) is faster for
  forests that are not purified, with bitwise the same predictions. From
  512 rows on, each tree of more than 64 leaves gets an index of its
  leaves (a k-d tree, built per call), so a row is tested only against
  the leaves that may hold it instead of against every leaf.
- Fitting regression forests with `split_structure = "tree_then_leaf"`
  is much faster and needs less memory. Fitted models are bitwise
  unchanged.
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  keeps the purified components as weighted rectangles instead of dense
  grids over all split points (with `mode = 2`, the default, which
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  and `rpf(purify = TRUE)` use). Purification at `max_interaction = 3`
  and above needs much less memory and no longer runs out of it on
  larger data. The components equal those of `mode = 1` up to rounding,
  and
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  stores the rectangles.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  evaluates all components in one pass over the rows of `new_data`, and
  skips those the forest does not hold, instead of one pass per
  component, which is faster.
- Classification uses the split search of regression: one implementation
  of the five split structures, which differ between the two only in how
  a candidate split is scored. So classification gets the faster
  `"tree_then_leaf"` and `"leaf"` and the smaller finished forests of
  regression, with bitwise the same fits.
- [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md),
  and so `rpf(prune = TRUE)` and `rpf(tune = )`, no longer build the
  purified store of all components. Each tree family is purified in
  turn, its components are evaluated at the training rows (and out of
  bag) and then freed, and the pruned model stores the pieces of its
  kept components only, from which
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  evaluate them.
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  keeps every piece of every component, which grows steeply with the
  order, so the whole call `rpf(n_splits = NULL, prune = TRUE)` needs
  much less memory at high orders, where it could run out of memory
  while purifying, and less at order 3, in the same time. The pruned
  models are the same up to rounding: the same choices, kept components
  and display. The forest in the result is left unpurified (call
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  for all its components), and an earlier purification of it is no
  longer redone.
- A fitted forest takes less memory: the split search’s per-leaf caches
  are kept apart from the leaves and freed when a tree family is
  finished, so a leaf takes 80 instead of 232 bytes; finished
  classification forests now also drop their training rows and caches.
  Fitted models are bitwise unchanged.
- Classification with `loss = "logit"` fits faster again, with bitwise
  the same fits. The terms of the loss under the current weights do not
  depend on the candidate split, so each observation’s terms are now
  computed once and refreshed only when its weights change, and the
  terms of classes the observation is not in (which are zero) are
  skipped: a candidate takes one exponential per logit column and one
  logarithm per observation of the leaf, instead of two exponentials and
  two logarithms per column plus two logarithms.
- Growing with `split_structure = "leaf"` (the default) is much faster,
  with bitwise the same fits. A leaf’s rows are put in order by each
  row’s integer rank in the predictor (ranked once per tree family) with
  a radix sort, instead of a comparison sort of the values, and the leaf
  keeps that order with the first position of each distinct value, so a
  candidate split no longer passes over the leaf’s values; the
  thresholds are drawn without allocating memory. Classification gains
  less.
- Pruning
  ([`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md),
  and so the default call) is faster on larger data, with bitwise the
  same models. The numbers of training rows in the leaves’ boxes,
  counted in the pass over the tree families, are kept for the later
  passes (the display’s coarsening, its out-of-bag sums, the store of
  the kept components) instead of being counted again, for trees of up
  to 4 predictors; and the index of the components’ rectangles finds its
  cuts without sorting.
- Choosing the number of splits out of bag (`n_splits = NULL`, the
  default) needs much less memory while growing, in the same time and
  with the same fits. Each tree family’s out-of-bag predictions are
  added to the checkpoints’ sums as the family delivers them, instead of
  being kept until the stop is decided; the pool of split candidates
  keeps no record of its entries; and
  `options(randomPlantedForest.adaptive_pause_gb)`, the memory of the
  families waiting for the stop, defaults to 0.25 GB instead of 1 GB.
- [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  of a large forest is several times as fast: the forest’s lists of
  trees and leaves are allocated at their final lengths instead of grown
  by one element at a time, which copied each list for every element.
  The blob is the same.

## randomPlantedForest 0.5.0

### Breaking changes

- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  now defaults to `max_interaction = 2`, main effects and pairwise
  interactions, instead of `1`, an additive model
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
  Values above the number of predictors are capped silently instead of
  with a message.
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  now defaults to `loss = "exponential"` for classification instead of
  `"L2"`, so classification fits give proper probability estimates by
  default; regression still uses `"L2"`
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  arguments are reordered by purpose (forest size, split search, loss,
  other). Pass arguments after the data arguments by name
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  error on unknown arguments, such as misspelled ones, instead of
  silently ignoring them
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  is a regular function instead of an S3 generic
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `preprocess_predictors_predict()` is no longer exported; it is an
  internal step of [`predict()`](https://rdrr.io/r/stats/predict.html)
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- The package now requires R \>= 4.1.0
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- For a given seed, classification models differ from 0.4.0 due to the
  fixes and speedup below, but regression models are unchanged
  ([\#57](https://github.com/PlantedML/randomPlantedForest/issues/57),
  [\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- For direct users of the C++ object (`$fit`): `set_data()` no longer
  fits, call `fit()` afterwards, and an unknown classification loss is
  an error instead of a silent fallback to L2
  ([@jyliuu](https://github.com/jyliuu),
  [\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  no longer has a `cv` argument, which never had an effect: the C++
  cross-validation was a no-op
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Bug fixes

- [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  error on missing values in `new_data`, naming the affected columns.
  Previously such rows silently received arbitrary predictions, as
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  has no handling for missing values.
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  likewise reports missing values in the training data by column
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Classification fits with a fixed seed could differ between runs. When
  a split replaced a leaf, its new split candidates could be lost or
  attached to the wrong leaf, depending on memory layout
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `predict(type = "class")` breaks probability ties by level order
  instead of at random, so class predictions are deterministic and no
  longer advance the random number generator
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  returns its input invisibly, as documented
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Logical predictors in the formula interface are used as a single 0/1
  column, as in the x/y interface. Previously they were expanded into
  two redundant indicator columns, which also counted against
  `max_interaction`
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Performance

- Prediction is much faster: leaves are scanned from flat,
  cache-friendly arrays (~15x single-threaded), and rows are split
  across threads
  ([\#65](https://github.com/PlantedML/randomPlantedForest/issues/65)).
  [`predict()`](https://rdrr.io/r/stats/predict.html) gains `nthreads`,
  defaulting to the `nthreads` used for fitting.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  is much faster: purified forests look up only the trees of the
  requested component, and each distinct input row is predicted once,
  which pays off for low-cardinality features
  ([\#59](https://github.com/PlantedML/randomPlantedForest/issues/59),
  [\#65](https://github.com/PlantedML/randomPlantedForest/issues/65)).
- Classification fits are slightly faster, as redundant work during
  model construction was removed
  ([\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).

### Other improvements

- New “Getting started” vignette covering regression, classification,
  component decomposition, recipes, and saving models
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- The [`print()`](https://rdrr.io/r/base/print.html) method is
  overhauled and gains a
  [`format()`](https://rdrr.io/r/base/format.html) counterpart: it
  reports the forest size, interaction degree, purification state, split
  search settings and number of threads
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `deterministic = TRUE` with `ntrees > 1` warns, as all tree families
  are then identical
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Rewritten documentation for
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
  including guidance on choosing parameters and on how purification
  interacts with
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Internals

- The C++ core no longer depends on `Rcpp`: it uses standard C++ types
  and exceptions, and an `Rcpp` layer (`src/rcpp_interface.*`) converts
  at the boundary. This is groundwork for bindings in other languages
  such as Python ([@jyliuu](https://github.com/jyliuu),
  [\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).
- `cli` and `rlang` are now imported, `mvtnorm` is only needed to build
  the `pkgdown` site
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

## randomPlantedForest 0.4.0

It is now possible to serialize and de-serialize an `rpf` object via
marshalling, meaning a fitted model can be saved to disk or dispatched
to a worker process for parallelization or encapsulation, as is done in
`mlr3`. This also re-opens the door for the `mlr3extralearners` wrapper,
which was removed due to the lack of serialization support.

- New
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  /
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  serialize a fitted forest to a plain R list and back, making
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)-based storage of
  rpf models possible
  ([\#52](https://github.com/PlantedML/randomPlantedForest/issues/52)).
  Purified forests restore their purified state directly; training data
  is only embedded with `include_data = TRUE` (required to
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  after restoring). Blobs record the blob format version and the package
  version they were created with; restoring under an older package than
  the one that saved the model warns. Malformed or corrupt blobs error
  instead of reading out of bounds.
- New
  [`rpf_is_valid()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_is_valid.md)
  checks whether an rpf object’s internal model is usable,
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  now give an actionable error for `rpf` objects restored via
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) without marshaling.
- New
  [`bundle::bundle()`](https://rstudio.github.io/bundle/reference/bundle.html)
  method for rpf models wrapping the marshaling API.

### Behavior changes

- The default for `delta` changed from 0 to 0.001. With `delta = 0`,
  splits producing single-class nodes have infinite logit loss and are
  always rejected, which also degraded binary logit fits.

### Fixes and improvements

- Fixed multiclass classification with `loss = "logit"`, which
  effectively never split and predicted near-uniform class probabilities
  ([\#40](https://github.com/PlantedML/randomPlantedForest/issues/40)).
  The C++ logit loss is a reference-class multinomial formulation
  expecting `K-1` indicator columns, but the R wrapper passed a full
  `K`-column one-hot matrix, pinning the implicit reference-class
  probability to zero. Multiclass logit outcomes are now encoded with
  the first factor level as reference class:
  - `predict(type = "prob")` reconstructs all `K` class probabilities
    from the `K-1` logits; rows sum to 1 exactly.
  - `predict(type = "numeric"/"link")` now returns `K-1` columns named
    after the non-reference levels (previously `K` columns).
  - [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
    on multiclass logit fits returns per-class components for the `K-1`
    non-reference levels; `target_levels` reflects this.

## randomPlantedForest 0.3.0

### Major changes ([\#61](https://github.com/PlantedML/randomPlantedForest/issues/61))

- New
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  arguments controlling split-candidate sampling:
  - `split_structure = "leaves"`: Defines what a split candidate is and
    how candidates are drawn. One of `"leaves"` (default), `"hist"`,
    `"cur_trees_1"`, `"cur_trees_2"`, or `"res_trees"`; see
    [`?rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
    for details.
  - `max_candidates = 50`: Maximum number of split candidates sampled
    per iteration.
  - `split_decay_rate = 0.1`: Exponential aging of repeatedly drawn but
    unchosen split candidates. `split_decay_rate = 0` corresponds to no
    aging and uniform sampling.
  - `delete_leaves = TRUE`: Whether a parent leaf is deleted when
    splitting along an existing dimension.
- **Fitting results change**: The new candidate-sampling defaults and a
  reworked internal RNG mean that fits are not reproducible against
  previous versions, even with the same seed. Install an older commit if
  exact reproduction of previous results is required.
- Seeded fits are now reproducible regardless of `nthreads`: per-tree
  seeds are drawn from R’s RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) gives identical
  forests for serial and multithreaded fits.
- Substantial speedups in fitting (cached per-leaf orderings, prefix
  sums) and reduced memory use (training-only buffers are released after
  each tree family is built).
- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  gains arguments:
  - `mode = 2`: Purification algorithm; `2` is a new fast exact method,
    `1` is the legacy grid-based path.
  - `nthreads = NULL`: Purification is now multithreaded, defaulting to
    the fit’s `nthreads` setting.
  - `maxp_interaction = NULL`: Optionally only compute purified
    components up to this interaction order.
- New
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  argument `export_forest = FALSE`: The flattened forest is no longer
  stored in the fitted object by default, so `rpf_object$forest` is
  `NULL` unless `export_forest = TRUE`. This reduces object size;
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  are unaffected.
- `preprocess_predictors_predict()` is now exported.
- Fixed a memory bug in the legacy purification path where the grid was
  sized one element too large, causing out-of-bounds reads (crashes on
  Windows, silently wrong purification results elsewhere).
- Fixed a crash on Windows when fitting with `nthreads > 1`, caused by a
  `thread_local` buffer with a non-trivial destructor being destroyed at
  thread exit.

### Other changes

- Internals in `src/` have been refactored into modular sub-files
  ([\#53](https://github.com/PlantedML/randomPlantedForest/issues/53))
- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  now errors if a regression target is combined with a `loss` other than
  `"L2"`.
- Allow features of type `logical`, which are now converted via
  `as.integer`.
- The `parallel = TRUE|FALSE` argument in
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  has been substituted by an `nthreads = 1L` argument, allowing for more
  flexible parallelization. The previous behavior only allowed for
  either no parallelization or using n-1 of n available cores. The new
  implementation should be reasonably robust and the default behavior
  remains serial execution.
- Remove `SystemRequirements` field from `DESCRIPTION`: Now the default
  C++ version is C++17 and with a minor change to internal use of random
  numbers, `randomPlantedForest` is now compatible with C++11 through
  C++23.
- Add `remainder` term to `predict_components` output for case where
  `max_interaction` supplied is smaller than `max_interaction` in `rpf`
  fit. In that case, the `m` values don’t sum up to the global
  predictions, so we add a remainder to allow reconstruction of that
  property.

## randomPlantedForest 0.2.1

- Add `glex` class to output of
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
  for extended functionality available with
  [`glex`](https://github.com/PlantedML/glex).
- Add `target_levels` vector to output of
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  to aid multiclass handling. Keeping track of levels is somewhat
  awkward since column names of `$m` need to be identifiable regarding
  the target level.

## randomPlantedForest 0.2.0

- Added a `NEWS.md` file to track changes to the package.
