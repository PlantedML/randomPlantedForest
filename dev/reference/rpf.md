# Random Planted Forest

Fits a random planted forest (Hiabu, Mammen and Meyer, 2020) for
regression or classification. The model is a sum of trees that each
split on at most `max_interaction` predictors, so it decomposes into
main effects and interaction components; see
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md).

## Usage

``` r
rpf(x, ...)

# S3 method for class 'formula'
rpf(
  formula,
  data,
  max_interaction = 3,
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
  prune = NULL,
  max_panels = 10,
  keep_forest = FALSE,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
)

# S3 method for class 'data.frame'
rpf(
  x,
  y,
  max_interaction = 3,
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
  prune = NULL,
  max_panels = 10,
  keep_forest = FALSE,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
)

# S3 method for class 'matrix'
rpf(
  x,
  y,
  max_interaction = 3,
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
  prune = NULL,
  max_panels = 10,
  keep_forest = FALSE,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
)

# S3 method for class 'recipe'
rpf(
  x,
  data,
  max_interaction = 3,
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
  prune = NULL,
  max_panels = 10,
  keep_forest = FALSE,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
)
```

## Arguments

- x, y:

  `[data.frame | matrix | recipe]`, `[numeric | factor]`: Predictors and
  outcome for the x/y interface. `x` holds the predictors, or is a
  [`recipe`](https://recipes.tidymodels.org/reference/recipe.html) (then
  supply `data` instead of `y`). The predictors need unique names that
  do not contain `:`, which joins predictor names in component names
  such as `x1:x2`. `y` is the outcome vector: **numeric** for
  regression, **factor** for classification. A numeric outcome must be
  finite. A factor outcome needs at least two classes; levels that do
  not occur are dropped with a warning.

- ...:

  Only the old names of renamed arguments (`ntrees`, `splits`,
  `split_try`, `t_try`), deprecated: they work, with a warning. Any
  other (e.g. misspelled) argument is an error. See "Renamed arguments"
  in Details.

- formula:

  `[formula]`: A formula with the outcome on the left-hand side and the
  predictors on the right-hand side, e.g. `y ~ x1 + x2`. Interaction
  terms such as `x1:x2` or `x1 * x2` are an error: they would add a
  predictor named `x1:x2`, and `:` is how components are named. The
  forest finds interactions itself, up to `max_interaction`, so list the
  predictors as `y ~ x1 + x2`. To fit a product as a predictor of its
  own, write `I(x1 * x2)`, or add it to `data` as a column with its own
  name. Transformed terms of one predictor, such as `log(x1)`, are
  allowed.

- data:

  `[data.frame]`: Data containing the predictors and the outcome, used
  with `formula` or a recipe `x`.

- max_interaction:

  `[integer(1): 3]`: Maximum number of predictors a single tree may
  split on. `1` fits main effects only (an additive model), `2` adds
  pairwise interactions, `3` (the default) three-way interactions, and
  so on. `0` uses all predictors and values above the number of
  predictors are reduced to it.

- n_families:

  `[integer(1): 200]`: Number of tree families, i.e. the size of the
  forest. Each family is grown on a bootstrap sample and predictions are
  averaged. Was `ntrees`.

- n_splits:

  `[integer(1) | NULL: NULL]`: Number of split rounds per tree family.
  The main tuning parameter: more splits fit more complex functions. A
  round splits the best of its candidates; near saturation, when every
  candidate's leaf is too small to split, a round can make no split.
  `NULL` (the default) chooses it out of bag while the forest grows
  (with `deterministic = TRUE`, which has no out-of-bag rows, `NULL`
  means 30). See "Choosing `n_splits` out of bag" in Details. Was
  `splits`.

- split_structure:

  `[character(1): "leaf"]`: How split candidates are formed and sampled;
  one of `"leaf"`, `"leaf_binned"`, `"tree_then_leaf"`,
  `"tree_then_all_leaves"` or `"target_tree_then_leaf"`. See Details.
  The names of 0.5.0 (`"leaves"`, `"hist"`, `"cur_trees_1"`,
  `"cur_trees_2"`, `"res_trees"`) still work, with a deprecation
  warning.

- n_try_split_points:

  `[integer(1) | NULL: NULL]`: Number of split points (thresholds) tried
  per split candidate. `NULL` means 3 for `split_structure = "leaf"` and
  10 for the other structures. Was `split_try`.

- share_candidates:

  `[numeric(1) | NULL: NULL]`: Share in `(0, 1]` of all possible splits
  drawn as candidates in each round, capped at `max_candidates`. `NULL`
  means 0.1 for `split_structure = "leaf"` and 0.4 for the other
  structures. See "Choosing parameters" in Details. Was `t_try`.

- max_candidates:

  `[integer(1): 50]`: Maximum number of split candidates per round.

- split_decay_rate:

  `[numeric(1): 0.1]`: Down-weights possible splits that were sampled as
  candidates but not chosen. Each such round ages a split by one,
  choosing it resets its age to zero, and splits are sampled with weight
  `exp(-split_decay_rate * age)`. `0` samples uniformly.

- delete_leaves:

  `[logical(1): TRUE]`: Whether to delete a parent leaf when splitting
  along an existing dimension.

- min_leaf_rows:

  `[integer(1): 1]`: Minimum number of rows of the bootstrap sample on
  each side of every split, so in every leaf, for every
  `split_structure` (regression and classification). Split points are
  tried only among those that keep this many rows on each side, and
  leaves without one are not drawn (see "`split_structure`" in Details).
  `1` means no limit.

- interaction_constraints:

  `[list | NULL: NULL]`: Feature sets that restrict which interactions
  the forest may form, e.g. `list(c("x1", "x2"), c("x3", "x4"))`. Each
  set is given by predictor names or by column indices: with `x`, the
  columns of `x`; with a formula, the predictors in the order the
  formula names them (`y ~ .`: the order of the columns of `data`); with
  a recipe, names only (its steps decide the columns). A tree on the
  predictors `T` is grown only if `T` is a subset of one of the sets, as
  with xgboost's interaction constraints. The subsets of a set stay
  allowed, so the forest keeps the main effects it grows a pair from. A
  predictor in no set is not used. Applies together with
  `max_interaction`: set `max_interaction = 0` to let the sets alone
  decide. `NULL` means no constraint. One use is to refit the forest on
  the effects a first fit found, e.g. after dropping its small
  components.

- loss:

  `[character(1) | NULL: NULL]`: Loss function. `NULL` uses `"L2"` for
  regression, the only supported loss, and `"logit"` (the multinomial
  log-likelihood) for classification, with two classes and with more.
  Classification also supports `"exponential"` (with more than two
  classes the multiclass exponential loss SAMME), and `"L1"` and `"L2"`,
  which fit class indicators directly and do not yield proper
  probability estimates. See "Classification losses" in Details.

- exact_step:

  `[logical(1): TRUE]`: Only used with `loss = "logit"` (the
  classification default): whether the two leaves of each split get the
  exact minimisers of the penalised loss (Newton's method to convergence
  on the leaf's rows), or one Newton step from the sums the split search
  scored them by (`FALSE`, faster). The split search is the same either
  way. The exponential loss's step is exact by itself, so `exact_step`
  does not change `"exponential"` fits. See "Classification losses" in
  Details.

- prune:

  `[logical(1) | NULL: NULL]`: Whether to return the pruned model
  instead of the forest: the forest is pruned by
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  with its defaults on the training rows, and the result is its
  `rpf_pruned` model. `NULL` (the default) means `TRUE`, and `FALSE`
  with `deterministic = TRUE`, where `TRUE` is an error (pruning needs
  out-of-bag rows). `FALSE` returns the forest. See "Pruning" in
  Details.

- max_panels:

  `[numeric(1): 10]`: Used only when pruning: the most heat maps a kept
  component of three or more predictors may need; such a component is
  coarsened to at most `max_panels` heat maps and the model is pruned
  again (the `max_panels` of
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)).
  A whole number of at least 1, or `Inf` for no coarsening.

- keep_forest:

  `[logical(1): FALSE]`: Used only when pruning: whether the pruned
  model keeps the forest it was pruned from, in `$fit`. By default it
  does not: the model predicts from the components it keeps, so the
  forest's memory is free when `rpf()` returns, and the model is plain R
  data that [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) stores as
  it is. `TRUE` keeps the forest, for pruning it again with other
  settings
  ([`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  on `$fit`),
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
  or its own predictions; such a model is saved with
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md).
  With `prune = FALSE` it has no effect: the result is the forest. See
  "Pruning" in Details.

- tune:

  `[logical(1) | character(1) | rpf_tune_control: FALSE]`: Whether to
  tune the greediness (`n_try_split_points`, `share_candidates`) and
  `min_leaf_rows` out of bag (not with `deterministic = TRUE`, at least
  10 tree families). `TRUE` runs the full search
  ([`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)),
  `"quick"` the two-forest preset, and a
  [`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
  object your own. Every forest of the search has `n_splits` chosen out
  of bag and is pruned, so `tune` implies `n_splits = NULL` and
  `prune = TRUE`, and giving `n_try_split_points`, `share_candidates` or
  `min_leaf_rows` with it is an error. See "Tuning out of bag" in
  Details.

- purify:

  `[logical(1): FALSE]`: Whether to purify the forest after fitting,
  which
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  requires, with the defaults of
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).
  For other settings, call
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  after fitting instead. Not used with `prune = TRUE`: the pruned model
  stores the components it keeps.

- nthreads:

  `[integer(1): 1]`: Number of threads for fitting. Also the default for
  [predict()](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  and
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).

- export_forest:

  `[logical(1): FALSE]`: Whether to store the flattened forest in the
  returned object as `$forest` (with `prune = TRUE` and
  `keep_forest = TRUE`, in the forest: `$fit$forest`; a pruned model
  without its forest has neither). If `FALSE`, it is `NULL`, which saves
  memory.

- deterministic:

  `[logical(1): FALSE]`: Whether to fit without randomness: no
  bootstrap, the first `max_candidates` possible splits as candidates
  (`share_candidates` and `split_decay_rate` are ignored) and evenly
  spaced thresholds. The result does not depend on the seed, and all
  tree families are identical, so use `n_families = 1`. There are no
  out-of-bag rows, so `n_splits = NULL` means 30 and `prune = NULL`
  means `FALSE`. Mainly useful for testing and debugging.

- keep_level_order:

  `[logical(1) | character: FALSE]`: Which factor predictors keep the
  order of their levels instead of being ordered by the outcome. `FALSE`
  orders every factor by the outcome, `TRUE` keeps the level order of
  every factor, and a character vector keeps it for the named factor
  predictors only. See "Factor predictors" below.

- oob_path:

  `[logical(1) | integer: FALSE]`: Whether to track the out-of-bag (OOB)
  error along the splits (not with `deterministic = TRUE`). `TRUE`
  evaluates it at the checkpoints of
  [`oob_path_checkpoints()`](https://plantedml.com/randomPlantedForest/dev/reference/oob_path_checkpoints.md)
  for `n_splits`, four per doubling; an increasing integer vector in
  `1..n_splits` gives the checkpoints. The fitted forest is the same as
  without. See Details. With `n_splits = NULL` the path of the
  checkpoints grown is stored either way, and only `TRUE` or `FALSE` is
  allowed.

## Value

With the defaults, an object of class `rpf_pruned` (see below). With
`prune = FALSE`, an object of class `"rpf"` with model object contained
in `$fit`. With `oob_path`, `$oob_path` holds the OOB path: a data frame
with one row per checkpoint, `n_splits`, the corrected and the plain OOB
root mean squared error (`oob_rmse`, `oob_rmse_raw`), and the standard
errors of the two mean squared errors (`se_mse`, `se_mse_raw`). Its
attribute `n_rows` is the number of training rows with an OOB tree
family. For classification with `loss = "logit"` or `"exponential"`, or
more than two classes, the columns are the mean corrected and plain
criterion (`oob_loss`, `oob_loss_raw`) and their standard errors
(`se_loss`, `se_loss_raw`), and the attribute `criterion` names it (see
"Classification" in Details). With `n_splits = NULL`, `$params$splits`
is the chosen number of splits (`$params` keeps the names of 0.5.0, see
"Renamed arguments" in Details), `$oob_path` the path of the checkpoints
grown, and `$splits_oob` a list with the chosen `n_splits`, `stopped_at`
(the number of splits where growing stopped) and `capped` (whether it
stopped at the cap while the error was still falling, neither stop rule
reached). When pruning (the default), an object of class `rpf_pruned`
(see
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)),
whose element `fit` describes the forest with the elements above: the
pruned model's OOB path is `$fit$oob_path`, and with `n_splits = NULL`
the chosen number of splits `$fit$params$splits` and the stop
`$fit$splits_oob`. Without `keep_forest`, `fit` holds no forest
(`$fit$fit` is `NULL`); with `keep_forest = TRUE` it is the forest (not
purified). With `tune`, the `rpf_pruned` model of the chosen forest,
with an element `tune`: `search`, the search table (one row per forest
grown, in the order grown: `step`, `n_try_split_points`,
`share_candidates`, `min_leaf_rows`, the chosen `n_splits`,
`stopped_at`, `capped`, the forest's corrected OOB RMSE
`oob_rmse_forest` (as the pruning scores it: its rule `"D"` at
`tau = 0`, which keeps every nonzero component; the OOB path's
`oob_rmse` at the chosen `n_splits` leaves the rows without OOB families
of both halves uncorrected and is higher), its pruned model's `rule`,
`tau`, `n_kept` and corrected OOB RMSE `oob_rmse` (the score; for
classification both are the square root of the mean corrected
criterion), `seconds`, whether each step chose it (`chosen_step1`,
`chosen_step2`, `chosen_step3`; `NA` if not a candidate of the step) and
`final`); `steps`, each step's candidates and choice; `settings`, the
chosen `n_try_split_points`, `share_candidates`, `min_leaf_rows` and
`n_splits`; `control`, the
[`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
used; and `seconds`.

## Details

By default, `rpf()` returns the pruned model that can be shown, in one
call:

1.  it grows 200 tree families of up to three-way interactions, with the
    number of splits chosen out of bag while they grow
    (`n_splits = NULL`);

2.  it prunes the forest's components by
    [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
    with its defaults (`prune = NULL`): the sparse rule, chosen out of
    bag;

3.  it coarsens every kept component of three or more predictors to at
    most `max_panels = 10` heat maps and prunes again.

The result, of class `rpf_pruned`, has methods for
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
and [`print()`](https://rdrr.io/r/base/print.html). It holds the kept
components, not the forest (`keep_forest = TRUE` keeps the forest too).
`prune = FALSE` returns the forest (class `rpf`) instead, and
`n_splits = <number>` fits a fixed number of splits.

### Choosing parameters

The defaults were chosen for regression and carried over to
classification.

`max_interaction` sets which effects the model can represent: `1` for an
additive model, `2` to add pairwise interactions, `3` (the default)
three-way interactions. Order 3 is often more accurate than order 2.
Higher orders gain only where the truth has such interactions, cost much
more time (no limit, `0`, is often much slower), and the pruned model
rarely keeps a component of more than three predictors, so it is shown
in at most `max_panels` heat maps anyway. With a higher
`max_interaction`, splits are spread over more possible components.

`n_splits` sets how closely the model follows the data. By default it is
chosen out of bag while the forest grows (see "Choosing `n_splits` out
of bag"), so it need not be tuned.

`n_families` averages over bootstrap samples to reduce variance, and
fitting and prediction time grow linearly with it. Each tree family is
itself a complete model, a sum of trees covering all components, rather
than a single tree, but averaging more of them still pays off: more
families usually lower the test error a little, and the pruned models
are usually sparser. Hence the default of 200, at the price of time; use
fewer for a quick look.

The split search parameters (`split_structure`, `n_try_split_points`,
`share_candidates`, `max_candidates`, `split_decay_rate`) trade speed
for a more thorough, more greedy search. For `"leaf"`, a less greedy
search is usually better, so its default is 3 split points per candidate
and 10% of the possible splits per round, against 10 and 40% for the
other structures. `tune` searches this greediness out of bag.

For classification, `"logit"` (the default) and `"exponential"` fit on
the link scale and are transformed to probabilities; see "Classification
losses". `"L1"` and `"L2"` fit class indicators directly and their
probabilities are truncated to `[0, 1]`.

### Classification losses

`"exponential"` and `"logit"` fit scores on the link scale: with two
classes one score, the log-odds of the second class; with more, one
score per class, summing to 0 across the classes. The probabilities are
the softmax of the scores. `"exponential"` with more than two classes is
the multiclass exponential loss SAMME of Zhu, Zou, Rosset and Hastie
(2009), whose row weight is `exp(-g)`, `g` the score of its class.
`"logit"` is the multinomial logit; it fits `K - 1` scores against the
first level (whose score is 0), which
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
(`type = "link"`),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
[`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md)
and the pruned model report with the first level's 0 joined and centred
across the `K` classes, so both losses report the same layout. In exact
arithmetic the fit and its pruned model do not depend on which level is
first; in floating point, candidate splits that tie exactly (as
duplicated bootstrap rows make) can be broken differently.

Every split moves the scores of the rows of each new leaf by one step,
which minimises the leaf's penalised loss: the loss of its rows, plus
that of one pseudo-row (split evenly over the classes) at the leaf's
current mean score. The pseudo-row pulls the step towards no change, so
a leaf whose rows all have one class moves towards a finite level (about
`log(2 n + 1)` for `n` rows and two classes) and is not moved further
once it is there, and no clamp of the class shares is needed. A split is
scored by its gain over the parent leaf's own step, and each new leaf
pays for its pseudo-row, so splitting a leaf of one class gains nothing,
nor does repeating a split whose leaves are at their optimum, and such
splits are not made. A round in which no candidate gains makes no split.
The candidates' gains come from sums over their two sides (the
second-order Newton gain for `"logit"`, the exact gain for
`"exponential"`), built as regression builds its sums (along the leaf's
sorted rows, so that a split point costs the same whatever the leaf's
size; for `"tree_then_leaf"` one pass over the leaf per split point).
`"exponential"`'s step is exact; with `exact_step = TRUE` (the default)
`"logit"`'s steps are made exact too, otherwise they are one Newton
step.

`"logit"` is the default, with two classes and with more: its test
log-loss is usually as low as that of `"exponential"` or lower, more so
with three or more classes. `"exponential"` can need many more splits
with many classes.

### Factor predictors

Factor predictors, and character predictors, which are converted to
factors, are coded as integers `1, 2, ...` in an order of their levels,
and splits are thresholds on these codes. A split therefore separates
the levels below a cut from those above it in that order.

By default, the levels are ordered by the outcome, so that the levels
most alike in the outcome are neighbours (Wright and König, 2019):

- regression: by the mean outcome of each level;

- binary classification: by the share of the second outcome level;

- more than two classes: by the first principal component of the class
  proportions of each level (Coppersmith et al., 1999).

This applies to ordered factors too. With `keep_level_order`, the chosen
factors are instead coded in the order of
[`levels()`](https://rdrr.io/r/base/levels.html), as given, ordered or
not, so that every split separates neighbouring levels in that order.

The coding is stored in the fit as `$factor_levels`, the levels of each
factor in the order of their codes.
[predict()](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
and
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
code new data by these levels, by name, so the order of the levels in
`new_data` does not matter.

### `split_structure`

In each round, a `share_candidates` share of all possible splits (capped
at `max_candidates`) is drawn as candidates with weights
`exp(-split_decay_rate * age)`. `split_structure` defines what a
candidate is and how its split points are tried. The names say the
candidate, then how the leaf to split is found.

- `"leaf"` (was `"leaves"`): Split candidates are (leaf,
  split-dimension) pairs. For each sampled candidate,
  `n_try_split_points` split points are drawn uniformly without
  replacement from the distinct values of that predictor among the
  leaf's rows, leaving out the smallest and, with `min_leaf_rows`, those
  that leave fewer rows on either side (all of them if there are fewer),
  and evaluated to choose the best split.

- `"leaf_binned"` (was `"hist"`): As `"leaf"`, but split points are
  drawn from boundaries of quantile bins computed once per predictor,
  which makes evaluating them faster on large data. With
  `min_leaf_rows`, only from the boundaries that leave that many rows on
  each side, and candidates without one are not drawn.

- `"tree_then_leaf"` (was `"cur_trees_1"`): Split candidates are
  (current-tree, split-dimension) pairs. For each sampled candidate,
  perform `n_try_split_points` evaluations. Each evaluation samples a
  leaf from the set of valid current trees (with probability
  proportional to its number of available split points) and then
  uniformly samples a single split point within that leaf. With
  `min_leaf_rows`, a leaf's available split points are those that leave
  that many rows on each side, so a leaf without one is not drawn.

- `"tree_then_all_leaves"` (was `"cur_trees_2"`): Split candidates are
  (current-tree, split-dimension) pairs. For each sampled candidate,
  iterate through every valid leaf. Within each leaf, sample
  `n_try_split_points` split points uniformly and evaluate them. With
  `min_leaf_rows`, only split points that leave that many rows on each
  side, and leaves without one are skipped.

- `"target_tree_then_leaf"` (was `"res_trees"`): Split candidates are
  resulting trees. For each sampled candidate, run `n_try_split_points`
  evaluations by sampling a (split-dimension, leaf) pair from all valid
  pairs (with probability proportional to its number of available split
  points), then uniformly sampling one split point within that pair.
  With `min_leaf_rows`, the available split points are those that leave
  that many rows on each side, as for `"tree_then_leaf"`.

### Renamed arguments

Four arguments and the split structures were renamed after version
0.5.0: `ntrees` is now `n_families`, `splits` is `n_splits`, `split_try`
is `n_try_split_points` and `t_try` is `share_candidates`; the split
structures are listed above. The prefix gives the unit: `n_` a count,
`share_` a share in `(0, 1]`, `max_` a cap. The old names still work and
give the same fit, with a deprecation warning once per session for each
old name (as lifecycle's soft deprecation;
`options(lifecycle_verbosity = "warning")` warns every time, `"quiet"`
never); giving an argument under both names is an error. The fitted
model's `$params` keeps the old names (`ntrees`, `splits`, `split_try`,
`t_try`, and the old `split_structure` values), so code that reads them,
and models saved by 0.5.0, keep working.

### The draw of candidates for `"leaf"`

With `split_structure = "leaf"` (regression and classification), the
global option `randomPlantedForest.exact_draw` chooses how the
candidates are drawn. With `FALSE` (the default), the candidates are
drawn from a Fenwick tree over the weights, kept up to date between
rounds, which costs a few random numbers per candidate. With `TRUE`, the
exact draw of earlier versions, every possible split gets a random
number in every round. Both draw from the same distribution, but with
other random numbers, so the same seed gives another (equally good)
forest. Ties between equally good candidates may go to different
candidates under the two draws, so with many splits, near saturation,
their fits differ in distribution, though not in accuracy. With
`min_leaf_rows`, both leave out the same candidates: those whose leaf
has no threshold with that many rows on each side. The fast draw is
faster on large data with many splits. The option is read when the
forest is fitted and stored in `params$exact_draw`; it does not affect
the other split structures. To reproduce a regression fit of version
0.5.0, set it and give 0.5.0's defaults, as the example "A regression
forest of version 0.5.0" does:

    options(randomPlantedForest.exact_draw = TRUE)

Classification fits of 0.5.0 cannot be reproduced: they changed for
every split structure and loss (see NEWS: classification uses
regression's split structures, the largest value of a predictor can be
split off, and three or more classes with the exponential loss use
SAMME).

### Reproducibility

With the same seed, `rpf()` gives the same forest on the same platform,
also with more threads. Another platform can give another forest for the
same seed, because the integer sampling of the C++ standard library
(`std::uniform_int_distribution`) differs between implementations, such
as libc++ and libstdc++.

### The out-of-bag path

With the same seed, the forest grown with `k` splits is exactly the
first `k` splits of every tree family of a forest grown with more:
nothing in the fit depends on the total number of splits. Each family is
grown on a bootstrap sample, so every training row is out of bag for the
families that did not draw it (about 37% of them). With `oob_path`, each
family's prediction at its OOB rows is updated split by split while it
grows, and copied at each checkpoint. So one fit gives the OOB error of
every forest with fewer splits, for example to choose `n_splits` without
cross-validation. This costs `2 * nrow(x) * length(checkpoints)` doubles
and a little time; the fit itself is unchanged.

The OOB prediction of a row is the mean over its OOB families, and the
plain error (`oob_rmse_raw`) is the root mean squared error of these
predictions over the rows with an OOB family. That mean averages fewer
families than the forest, so the plain error overstates the forest's
error by the families' extra variance, more so with more splits. The
corrected error (`oob_rmse`) subtracts this excess,
`(1 / m_i - 1 / n_families) * v_i` for row `i` with `m_i` OOB families,
where `v_i` is the variance across families estimated from the OOB means
over the families of even and of odd index,
`v_i = (g_even - g_odd)^2 / (1 / m_even + 1 / m_odd)`; this is also the
correction of
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)'
OOB choice. A row whose OOB families are all even or all odd is not
corrected here
([`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
gives it the mean estimate of the other rows): about 2% of the rows at
20 families, 18% at 10 and over half at 4 (none at 50), so with few
families the corrected error is biased upwards, and above the corrected
error
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
gives the same forest (rule `"D"` at `tau = 0`). Factor levels are
ordered by the outcome of all training rows: by the mean outcome
(regression), by the share of the second class (binary classification)
or by the first principal component of the class proportions (three or
more classes). So a row's own outcome enters the ordering its OOB
families split on, and the OOB error is somewhat optimistic with factor
predictors, most so with many levels.

### Choosing `n_splits` out of bag

With `n_splits = NULL` (not deterministic), the number of splits is
chosen while the forest grows, by the corrected OOB error of the
out-of-bag path. The families grow to checkpoints, four per doubling of
the splits
([`oob_path_checkpoints()`](https://plantedml.com/randomPlantedForest/dev/reference/oob_path_checkpoints.md)
up to 25,600), and the forest's corrected error is followed checkpoint
by checkpoint. Growing stops once the lowest error is one doubling of
the splits back (for example, no new minimum from 100 to 200 splits), or
once it is on a plateau: the lowest error improved over the last
doubling by at most
`getOption("randomPlantedForest.adaptive_tolerance", 1)` standard errors
of that improvement (paired over the training rows: the standard
deviation of the rows' differences in corrected error between the two
checkpoints, divided by the square root of their number). A tolerance of
0 stops only on the first rule. Growing also stops at 25,600 splits,
with a warning if the error was still falling there. The result is the
forest at the lowest error (ties go to fewer splits; a corrected error
at or below 0, a correction that overshot, is not chosen unless all
are): the forest that `rpf()` with that many `n_splits` and the same
seed returns. The first rule alone usually chooses the same `n_splits`
as the lowest error along a long fixed budget of splits, while growing
fewer. The plateau stop usually chooses the same `n_splits` as the first
rule alone, while growing fewer splits; it matters most on long, flat
out-of-bag curves.

The families grow one at a time, `nthreads` together, each to the stop
estimated from the families grown before it. A grown family keeps only
its OOB predictions at the checkpoints and a record of its splits (40
bytes a split), from which the chosen forest is rebuilt, without a
refit. A family that reaches the estimate waits there with its growing
state while the waiting families fit in
`getOption("randomPlantedForest.adaptive_pause_gb", 0.25)` GB, so that
it can grow on cheaply if the stop comes later; otherwise it is freed,
and grown again from its seed if needed. The option trades memory for
time; the fit depends neither on it nor on `nthreads`. A smaller budget
lowers the peak memory, while a budget of 0 can cost much more time
(families grown again). Growing all families together would hold every
family's growing state until the stop, which on large data needs several
times the memory of growing them one at a time.

### Pruning

By default (`prune = NULL`, or `TRUE`), `rpf()` fits the forest, then
prunes its components with
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
at that function's defaults, on the training rows: the rule and the
threshold are chosen out of bag (by cross-validation below 10 tree
families, with a message), the sparse choice `select = "1se"`, and every
kept component of three or more predictors is shown in at most
`max_panels` heat maps before the model is pruned again. The result is
the pruned model, of class `rpf_pruned`, with methods for
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
and [`print()`](https://rdrr.io/r/base/print.html). It is the same, bit
for bit, as `prune_components(rpf(..., prune = FALSE), data)` from the
same random number state, for every interface. For other settings of the
pruning (`rule`, `tau`, `select`, `select_by`, folds or seed), fit with
`prune = FALSE` and call
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md).
Together with `n_splits = NULL`, one call chooses the number of splits
out of bag, prunes and returns the model to show. Pruning needs
out-of-bag rows, so with `deterministic = TRUE` `prune = NULL` means
`FALSE` and `prune = TRUE` is an error.

The pruned model does not keep the forest (`keep_forest = FALSE`): it
predicts, prints and shows its components from what it stores, the same,
bit for bit, as with the forest, and is a fraction of the forest's size.
What needs the forest is then an error that says so: pruning again,
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
the forest's own
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
and
[plot()](https://plantedml.com/randomPlantedForest/dev/reference/plot.rpf.md).
With `keep_forest = TRUE` the forest is in `$fit`, so that
`prune_components(model$fit, data, ...)` prunes it again.

### Classification

A classification forest goes through the same steps (the OOB path,
`n_splits = NULL`, `prune`, `max_panels` and `tune`) with its score
columns: one for a binary outcome, one per class otherwise. Every
out-of-bag choice is scored by the loss's own criterion, corrected for
the OOB family count as above:

- `loss = "L1"` or `"L2"`: the squared error of the raw scores summed
  over the columns (the Brier score, before the scores are truncated to
  probabilities). With two classes this is the regression criterion of
  the 0/1 outcome, and the choices are the regression ones, bit for bit.

- `loss = "logit"` or `"exponential"`: the log-loss of the probabilities
  that `predict(type = "prob")` returns. The correction is its
  second-order form,
  `(1 / m_i - 1 / n_families) * d' H d / (2 * (1 / m_even + 1 / m_odd))`
  with `d = g_even - g_odd` and `H = diag(p) - p p'` over the score
  columns (for the squared error `H = 2 I`, the correction above).

The pruning keeps or drops a component for all classes, sizes it by its
variance summed over the score columns centred across the classes, and
fits one coefficient per kept component shared by the classes (see
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)).
The defaults were chosen for regression and carried over.

### Tuning out of bag

With `tune`, `rpf()` searches the greediness of the split search
(`n_try_split_points` and `share_candidates`, on the ladder of
[`greediness_ladder()`](https://plantedml.com/randomPlantedForest/dev/reference/greediness_ladder.md))
and the minimum leaf size `min_leaf_rows`, without cross-validation.
Every forest of the search is grown with `n_splits` chosen out of bag
and pruned by
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
with its defaults (without the display), and forests are compared by the
corrected out-of-bag RMSE of their pruned model (for classification by
the root of its mean corrected criterion); the lowest wins. A pruned
model at 0 (every candidate's correction overshot, see
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md))
is not chosen unless all are. The search is coordinate-wise
([`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)):
with `tune = TRUE` the ladder from 10 / 0.4 to 1 / 0.0125 at
`min_leaf_rows = 1`, then `min_leaf_rows` 1, 5, 20 and 100 at the chosen
point, then the point's two ladder neighbours, about 10 forests in all.
`tune = "quick"` grows two forests, 3 / 0.1 and 1 / 0.025 at
`min_leaf_rows = 1`. The chosen forest's pruned model is returned, with
the display at `max_panels`, and the search in its element `tune`. Every
forest is grown from the same random number state, so the result is the
model that `rpf(n_splits = NULL, prune = TRUE)` returns at the chosen
settings from that state.

The full search often lowers the test error clearly against the
greediest point, 10 / 0.4, and does as well as 5-fold cross-validation
over the same settings, or better, at a small fraction of its cost; the
quick preset gets most of the gain. A search costs about as many fits as
it grows forests, and their prunings.

## References

Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
directly interpretable tree ensemble.
[doi:10.48550/arXiv.2012.14563](https://doi.org/10.48550/arXiv.2012.14563)

Wright, M. N. & König, I. R. (2019). Splitting on categorical predictors
in random forests. PeerJ 7:e6339.
[doi:10.7717/peerj.6339](https://doi.org/10.7717/peerj.6339)

Coppersmith, D., Hong, S. J. & Hosking, J. R. (1999). Partitioning
nominal attributes in decision trees. Data Mining and Knowledge
Discovery 3:197–217.
[doi:10.1023/A:1009869804967](https://doi.org/10.1023/A%3A1009869804967)

## See also

[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
and
[print()](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_pruned.md)
of the pruned model,
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
for its components,
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
for other pruning rules,
[`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
for `tune`, and
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
to save a model.

## Examples

``` r
# The default: splits chosen out of bag, pruned, every component shown
set.seed(1)
model <- rpf(mpg ~ ., data = mtcars)
model
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 9 components.
#> Intercept: 20.2843
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  disp      3.5%      
#>  cyl       3.4%      
#>  wt        3.1%      
#>  hp        2.4%      
#>  carb      2.1%      
#>  am        1.1%      
#>  gear      1.1%      
#>  drat      0.67%     
#>  vs        0.47%     
#> 
#> How it was chosen:
#>   The forest: 200 tree families of 15 splits (chosen out of bag; growing stopped at 25).
#>   From the forest's 175 nonzero components (purified): those with a share of at least 0.3% of the forest's Var(prediction), as they are (rule D).
#>   Rule and threshold chosen by out-of-bag error (1-SE rule; 5 folds, seed 682372961) among 16 candidates.
#>   Left out: 166 components, shares together 0.4% (largest qsec, 0.16%).
predict(model, mtcars[1:3, ])
#> # A tibble: 3 × 1
#>   .pred
#>   <dbl>
#> 1  20.5
#> 2  20.2
#> 3  25.6

# Classification, the same steps
set.seed(1)
model <- rpf(Species ~ ., data = iris, n_families = 50)
predict(model, iris[c(1, 51, 101), ], type = "prob")
#> # A tibble: 3 × 3
#>   .pred_setosa .pred_versicolor .pred_virginica
#>          <dbl>            <dbl>           <dbl>
#> 1     0.992            0.00826         0.000191
#> 2     0.00207          0.984           0.0140  
#> 3     0.000134         0.000840        0.999   

# The forest instead of the pruned model, with x and y, pairwise interactions
# and a fixed number of splits
set.seed(1)
rpfit <- rpf(
  x = mtcars[, c("cyl", "wt", "hp")],
  y = mtcars$mpg,
  max_interaction = 2,
  n_splits = 30,
  prune = FALSE
)
rpfit
#> ── Regression Random Planted Forest ────────────────────────────────────────────
#> Predictors: `cyl`, `wt`, and `hp`
#> 200 tree families with 30 splits each on 3 predictors, interactions to degree
#> 2.
#> ℹ Forest is not purified.
#> 
#> ── Tree growing 
#>      split_structure: leaf
#>   n_try_split_points: 3
#>     share_candidates: 0.1
#>       max_candidates: 50
#>     split_decay_rate: 0.1
#>        delete_leaves: TRUE
#> 
#> ℹ Fit using 1 thread, also the default for `predict()` and `purify()`.

# A factor that keeps its level order: splits separate neighbouring levels
cars <- mtcars
cars$gear <- factor(cars$gear, levels = c(3, 4, 5))
rpfit <- rpf(mpg ~ gear + wt, data = cars, keep_level_order = "gear")
rpfit$factor_levels
#> NULL

# The out-of-bag error along the splits, from one fit
set.seed(1)
rpfit <- rpf(
  mpg ~ .,
  data = mtcars,
  n_families = 20,
  n_splits = 100,
  oob_path = TRUE,
  prune = FALSE
)
rpfit$oob_path[which.min(rpfit$oob_path$oob_rmse), ]
#>    n_splits oob_rmse oob_rmse_raw  se_mse se_mse_raw
#> 12       42 2.557867     2.709333 1.91871   1.895507

# \donttest{
# A regression forest of version 0.5.0: its defaults and the exact draw
old <- options(randomPlantedForest.exact_draw = TRUE)
set.seed(1)
rpfit <- rpf(
  mpg ~ cyl + wt,
  data = mtcars,
  max_interaction = 2,
  n_families = 50,
  n_splits = 30,
  n_try_split_points = 10,
  share_candidates = 0.4,
  prune = FALSE
)
options(old)
# }

# \donttest{
# Tuning out of bag: the quick preset, two forests
set.seed(1)
x <- matrix(runif(400 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(400, sd = 0.3)
tuned <- rpf(x, y, n_families = 10, tune = "quick")
tuned$tune$search
#>   step n_try_split_points share_candidates min_leaf_rows n_splits stopped_at
#> 1    1                  3            0.100             1       84        168
#> 2    1                  1            0.025             1      100        168
#>   capped oob_rmse_forest rule tau n_kept  oob_rmse seconds chosen_step1
#> 1  FALSE       0.4025576    R 0.1      2 0.3393394   0.080         TRUE
#> 2  FALSE       0.4066480    R 0.1      2 0.3657293   0.059        FALSE
#>   chosen_step2 chosen_step3 final
#> 1           NA           NA  TRUE
#> 2           NA           NA FALSE
# }
```
