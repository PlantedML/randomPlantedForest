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
  max_interaction = 2,
  ntrees = 50,
  splits = 30,
  split_structure = "leaves",
  split_try = 10,
  t_try = 0.4,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  ...
)

# S3 method for class 'data.frame'
rpf(
  x,
  y,
  max_interaction = 2,
  ntrees = 50,
  splits = 30,
  split_structure = "leaves",
  split_try = 10,
  t_try = 0.4,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  ...
)

# S3 method for class 'matrix'
rpf(
  x,
  y,
  max_interaction = 2,
  ntrees = 50,
  splits = 30,
  split_structure = "leaves",
  split_try = 10,
  t_try = 0.4,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  ...
)

# S3 method for class 'recipe'
rpf(
  x,
  data,
  max_interaction = 2,
  ntrees = 50,
  splits = 30,
  split_structure = "leaves",
  split_try = 10,
  t_try = 0.4,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
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

  Not currently used, but required for extensibility. Unknown (e.g.
  misspelled) arguments are an error.

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

  `[integer(1): 2]`: Maximum number of predictors a single tree may
  split on. `1` fits main effects only (an additive model), `2` adds
  pairwise interactions, and so on. `0` uses all predictors and values
  above the number of predictors are reduced to it.

- ntrees:

  `[integer(1): 50]`: Number of tree families, i.e. the size of the
  forest. Each family is grown on a bootstrap sample and predictions are
  averaged.

- splits:

  `[integer(1): 30]`: Number of split rounds per tree family. The main
  tuning parameter: more splits fit more complex functions. A round
  splits the best of its candidates; near saturation, when every
  candidate's leaf is too small to split, a round can make no split.

- split_structure:

  `[character(1): "leaves"]`: How split candidates are formed and
  sampled; one of `"leaves"`, `"hist"`, `"cur_trees_1"`,
  `"cur_trees_2"`, or `"res_trees"`. See Details.

- split_try:

  `[integer(1): 10]`: Number of thresholds evaluated per split
  candidate.

- t_try:

  `[numeric(1): 0.4]`: Proportion in `(0, 1]` of all possible splits
  sampled as candidates in each round, capped at `max_candidates`.

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
  each side of every split, so in every leaf. Leaves without a threshold
  that keeps this many rows on each side are not drawn as split
  candidates. Only for regression with `split_structure = "leaves"`; `1`
  (no limit) otherwise.

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
  regression, the only supported loss, and `"exponential"` for
  classification. Classification also supports `"logit"`, which gives
  similar results but is slower, and `"L1"` and `"L2"`, which fit class
  indicators directly and do not yield proper probability estimates.

- delta:

  `[numeric(1): 0.001]`: Only used if `loss` is `"logit"` or
  `"exponential"`. Class proportions are truncated to
  `[delta, 1 - delta]` when computing the loss of a split. Must be less
  than 0.5, from where every proportion is truncated to the same value,
  and positive for `"logit"`: with `delta = 0`, nodes containing a
  single class have infinite loss and their splits are always rejected.

- epsilon:

  `[numeric(1): 0.1]`: Only used if `loss` is `"logit"` or
  `"exponential"`. Class proportions are truncated to
  `[epsilon, 1 - epsilon]` when computing the fit in a leaf. Unlike
  `delta`, this caps the size of individual leaf updates and acts as
  regularization: smaller values permit larger jumps on the link scale.
  Must be greater than 0 and less than 0.5: at 0 the leaf values of pure
  leaves are infinite, and from 0.5 on every leaf value is 0.

- purify:

  `[logical(1): FALSE]`: Whether to purify the forest after fitting,
  which
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  requires, with the defaults of
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).
  For other settings, call
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  after fitting instead.

- nthreads:

  `[integer(1): 1]`: Number of threads for fitting. Also the default for
  [predict()](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  and
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).

- export_forest:

  `[logical(1): FALSE]`: Whether to store the flattened forest in the
  returned object as `$forest`. If `FALSE`, `$forest` is `NULL`, which
  saves memory.

- deterministic:

  `[logical(1): FALSE]`: Whether to fit without randomness: no
  bootstrap, the first `max_candidates` possible splits as candidates
  (`t_try` and `split_decay_rate` are ignored) and evenly spaced
  thresholds. The result does not depend on the seed, and all tree
  families are identical, so use `ntrees = 1`. Mainly useful for testing
  and debugging.

- keep_level_order:

  `[logical(1) | character: FALSE]`: Which factor predictors keep the
  order of their levels instead of being ordered by the outcome. `FALSE`
  orders every factor by the outcome, `TRUE` keeps the level order of
  every factor, and a character vector keeps it for the named factor
  predictors only. See "Factor predictors" below.

## Value

Object of class `"rpf"` with model object contained in `$fit`.

## Details

### Choosing parameters

Start with `max_interaction` and `splits`. `max_interaction` sets which
effects the model can represent: `1` for an additive model, `2` (the
default) to add pairwise interactions, which can still be plotted as
heatmaps. Higher values are more flexible, but slower and the components
are harder to interpret. `splits` sets how closely the model follows the
data and is the parameter to tune. With a higher `max_interaction`,
splits are spread over more possible components, so tune both together.

`ntrees` averages over bootstrap samples to reduce variance and fitting
and prediction time grow linearly with it. The default is lower than the
usual 500 trees of a random forest because each tree family is itself a
complete model, a sum of trees covering all components, rather than a
single tree. Averaging therefore stops paying off sooner, usually after
a few dozen families. The split search parameters (`split_structure`,
`split_try`, `t_try`, `max_candidates`, `split_decay_rate`) trade speed
for a more thorough search.

For classification, `"exponential"` (the default) and `"logit"` fit on
the link scale and are transformed to probabilities, `"exponential"`
being faster. `"L1"` and `"L2"` fit class indicators directly and their
probabilities are truncated to `[0, 1]`.

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

In each round, a `t_try` fraction of all possible splits (capped at
`max_candidates`) is drawn as candidates with weights
`exp(-split_decay_rate * age)`. `split_structure` defines what a
candidate is and how its thresholds are evaluated.

- `leaves`: Split candidates are (leaf, split-dimension) pairs. For each
  sampled candidate, `split_try` thresholds are drawn uniformly without
  replacement from the distinct values of that predictor among the
  leaf's rows, leaving out the smallest and, with `min_leaf_rows`, those
  that leave fewer rows on either side (all of them if there are fewer),
  and evaluated to choose the best split.

- `hist`: As `leaves`, but thresholds are drawn from boundaries of
  quantile bins computed once per predictor, which makes evaluating them
  faster on large data.

- `cur_trees_1`: Split candidates are (current-tree, split-dimension)
  pairs. For each sampled candidate, perform `split_try` evaluations.
  Each evaluation samples a leaf from the set of valid current trees
  (with probability proportional to its number of available thresholds)
  and then uniformly samples a single threshold within that leaf.

- `cur_trees_2`: Split candidates are (current-tree, split-dimension)
  pairs. For each sampled candidate, iterate through every valid leaf.
  Within each leaf, sample `split_try` thresholds uniformly and evaluate
  them.

- `res_trees`: Split candidates are resulting trees. For each sampled
  candidate, run `split_try` evaluations by sampling a (split-dimension,
  leaf) pair from all valid pairs (with probability proportional to its
  number of available thresholds), then uniformly sampling one threshold
  within that pair.

### The draw of candidates for `leaves`

For regression with `split_structure = "leaves"`, the global option
`randomPlantedForest.exact_draw` chooses how the candidates are drawn.
With `TRUE` (the default), every possible split gets a random number in
every round. With `FALSE`, the candidates are drawn from a Fenwick tree
over the weights, kept up to date between rounds, which costs a few
random numbers per candidate instead. Both draw from the same
distribution, but with other random numbers, so the same seed gives
another (equally good) forest. Ties between equally good candidates may
go to different candidates under the two draws, so with many splits,
near saturation, their fits differ in distribution, though not in
accuracy. With `min_leaf_rows`, both leave out the same candidates:
those whose leaf has no threshold with that many rows on each side. The
fast draw pays off on large data with many splits. The option is read
when the forest is fitted and stored in `params$exact_draw`; it does not
affect classification or the other split structures.

    options(randomPlantedForest.exact_draw = FALSE)

### Reproducibility

With the same seed, `rpf()` gives the same forest on the same platform,
also with more threads. Another platform can give another forest for the
same seed, because the integer sampling of the C++ standard library
(`std::uniform_int_distribution`) differs between implementations, such
as libc++ and libstdc++.

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

## Examples

``` r
# Regression with formula
rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)

# Regression with x and y, allowing pairwise interactions
rpfit <- rpf(x = mtcars[, c("cyl", "wt", "hp")], y = mtcars$mpg, max_interaction = 2)

# Classification
rpfit <- rpf(Species ~ ., data = iris, loss = "exponential")

# A factor that keeps its level order: splits separate neighbouring levels
cars <- mtcars
cars$gear <- factor(cars$gear, levels = c(3, 4, 5))
rpfit <- rpf(mpg ~ gear + wt, data = cars, keep_level_order = "gear")
rpfit$factor_levels
#> $gear
#> [1] "3" "4" "5"
#> 
```
