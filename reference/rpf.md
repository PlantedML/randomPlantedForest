# Random Planted Forest

Fits a random planted forest (Hiabu, Mammen and Meyer, 2020) for
regression or classification. The model is a sum of trees that each
split on at most `max_interaction` predictors, so it decomposes into
main effects and interaction components; see
[`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md).

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
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
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
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
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
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
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
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  ...
)
```

## Arguments

- x, y:

  `[data.frame | matrix | recipe]`, `[numeric | factor]`: Predictors and
  outcome for the x/y interface. `x` holds the predictors, or is a
  [`recipe`](https://recipes.tidymodels.org/reference/recipe.html) (then
  supply `data` instead of `y`). `y` is the outcome vector: **numeric**
  for regression, **factor** for classification.

- ...:

  Not currently used, but required for extensibility. Unknown (e.g.
  misspelled) arguments are an error.

- formula:

  `[formula]`: A formula with the outcome on the left-hand side and the
  predictors on the right-hand side, e.g. `y ~ x1 + x2`.

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

  `[integer(1): 30]`: Number of splits per tree family. The main tuning
  parameter: more splits fit more complex functions.

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

- loss:

  `[character(1) | NULL: NULL]`: Loss function. `NULL` uses `"L2"` for
  regression, the only supported loss, and `"exponential"` for
  classification. Classification also supports `"logit"`, which gives
  similar results but is slower, and `"L1"` and `"L2"`, which fit class
  indicators directly and do not yield proper probability estimates.

- delta:

  `[numeric(1): 0.001]`: Only used if `loss` is `"logit"` or
  `"exponential"`. Class proportions are truncated to
  `[delta, 1 - delta]` when computing the loss of a split. Should be
  positive for `"logit"`: with `delta = 0`, nodes containing a single
  class have infinite loss and their splits are always rejected.

- epsilon:

  `[numeric(1): 0.1]`: Only used if `loss` is `"logit"` or
  `"exponential"`. Class proportions are truncated to
  `[epsilon, 1 - epsilon]` when computing the fit in a leaf. Unlike
  `delta`, this caps the size of individual leaf updates and acts as
  regularization: smaller values permit larger jumps on the link scale.

- purify:

  `[logical(1): FALSE]`: Whether to purify the forest after fitting,
  which
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  requires, with the defaults of
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md).
  For other settings, call
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  after fitting instead.

- nthreads:

  `[integer(1): 1]`: Number of threads for fitting. Also the default for
  [predict()](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md)
  and
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md).

- export_forest:

  `[logical(1): FALSE]`: Whether to store the flattened forest in the
  returned object as `$forest`. If `FALSE`, `$forest` is `NULL`, which
  saves memory.

- deterministic:

  `[logical(1): FALSE]`: Whether to fit without randomness: no
  bootstrap, the first `max_candidates` possible splits as candidates
  (`t_try` is ignored) and evenly spaced thresholds. The result does not
  depend on the seed, and all tree families are identical, so use
  `ntrees = 1`. Mainly useful for testing and debugging.

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

### `split_structure`

In each round, a `t_try` fraction of all possible splits (capped at
`max_candidates`) is drawn as candidates with weights
`exp(-split_decay_rate * age)`. `split_structure` defines what a
candidate is and how its thresholds are evaluated.

- `leaves`: Split candidates are (leaf, split-dimension) pairs. For each
  sampled candidate, `split_try` thresholds are drawn uniformly from the
  valid range within that leaf and evaluated to choose the best split.

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

## References

Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
directly interpretable tree ensemble.
[doi:10.48550/arXiv.2012.14563](https://doi.org/10.48550/arXiv.2012.14563)

## Examples

``` r
# Regression with formula
rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)

# Regression with x and y, allowing pairwise interactions
rpfit <- rpf(x = mtcars[, c("cyl", "wt", "hp")], y = mtcars$mpg, max_interaction = 2)

# Classification
rpfit <- rpf(Species ~ ., data = iris, loss = "exponential")
```
