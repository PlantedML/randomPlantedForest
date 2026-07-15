# Random Planted Forest

Random Planted Forest

## Usage

``` r
rpf(x, ...)

# S3 method for class 'data.frame'
rpf(
  x,
  y,
  max_interaction = 1,
  ntrees = 50,
  splits = 30,
  split_try = 10,
  t_try = 0.4,
  split_decay_rate = 0.1,
  max_candidates = 50,
  delete_leaves = TRUE,
  deterministic = FALSE,
  nthreads = 1,
  purify = FALSE,
  cv = FALSE,
  loss = "L2",
  delta = 0.001,
  epsilon = 0.1,
  split_structure = "leaves",
  export_forest = FALSE,
  ...
)

# S3 method for class 'matrix'
rpf(
  x,
  y,
  max_interaction = 1,
  ntrees = 50,
  splits = 30,
  split_try = 10,
  t_try = 0.4,
  split_decay_rate = 0.1,
  max_candidates = 50,
  delete_leaves = TRUE,
  deterministic = FALSE,
  nthreads = 1,
  purify = FALSE,
  cv = FALSE,
  loss = "L2",
  delta = 0.001,
  epsilon = 0.1,
  split_structure = "leaves",
  export_forest = FALSE,
  ...
)

# S3 method for class 'formula'
rpf(
  formula,
  data,
  max_interaction = 1,
  ntrees = 50,
  splits = 30,
  split_try = 10,
  t_try = 0.4,
  split_decay_rate = 0.1,
  max_candidates = 50,
  delete_leaves = TRUE,
  deterministic = FALSE,
  nthreads = 1,
  purify = FALSE,
  cv = FALSE,
  loss = "L2",
  delta = 0.001,
  epsilon = 0.1,
  split_structure = "leaves",
  export_forest = FALSE,
  ...
)

# S3 method for class 'recipe'
rpf(
  x,
  data,
  max_interaction = 1,
  ntrees = 50,
  splits = 30,
  split_try = 10,
  t_try = 0.4,
  split_decay_rate = 0.1,
  max_candidates = 50,
  delete_leaves = TRUE,
  deterministic = FALSE,
  nthreads = 1,
  purify = FALSE,
  cv = FALSE,
  loss = "L2",
  delta = 0.001,
  epsilon = 0.1,
  split_structure = "leaves",
  export_forest = FALSE,
  ...
)
```

## Arguments

- x, data:

  Feature `matrix`, or `data.frame`, or
  [`recipe`](https://recipes.tidymodels.org/reference/recipe.html).

- ...:

  (Unused).

- y:

  Target vector for use with `x`. The class of `y` (either `numeric` or
  [`factor`](https://rdrr.io/r/base/factor.html)) determines if
  regression or classification will be performed.

- max_interaction:

  `[1]`: Maximum level of interaction determining maximum number of
  split dimensions for a tree. The default `1` corresponds to main
  effects only. If `0`, the number fo columns in `x` is used, i.e. for
  10 predictors, this is equivalent to setting `max_interaction = 10`.

- ntrees:

  `[50]`: Number of trees generated per family.

- splits:

  `[30]`: Number of splits performed for each tree family.

- split_try:

  `[10]`: Number of split points to be considered when choosing a split
  candidate.

- t_try:

  `[0.4]`: A value in (0,1\] specifying the proportion of viable
  split-candidates in each round.

- split_decay_rate:

  `[0.1]`: Exponential decay factor for aging split-candidates. Possible
  splits are initiated with age=0. Whenever a possible split becomes a
  split_candidate (i.e. it has been drawn when drawing
  max(max_candidates , t_try \* possible options ) times) it ages by +1.
  The age of the single split-candidate with minimal loss is reset to
  zero. Split_candidates are sampled from Possible_splits with weight
  exp(-split_decay_rate\_ \* age). A high split_decay_rate means faster
  aging. split_decay_rate=0 results in no aging and uniform sampling.

- max_candidates:

  `[50]`: Maximum number of split-candidates sampled per iteration.
  Number of split_candidates in each round is given by
  max(max_candidates , t_try \* possible options).

- delete_leaves:

  `[TRUE]`: Whether to delete a parent leaf when splitting along an
  existing dimension.

- deterministic:

  `[FALSE]`: Choose whether approach deterministic or random.

- nthreads:

  `[1L]`: Number of threads used for computation, defaulting to serial
  execution.

- purify:

  `[FALSE]`: Whether the forest should be purified. Set to `TRUE` to
  enable components extract with
  [`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md)
  are valid. Can be achieved after fitting with
  [`purify()`](http://plantedml.com/randomPlantedForest/reference/purify.md).

- cv:

  `[FALSE]`: Determines if cross validation is performed.

- loss:

  `["L2"]`: For regression, only `"L2"` is supported. For
  classification, `"L1"`, `"logit"` and `"exponential"` are also
  available. `"exponential"` yields similar results as `"logit"` while
  being significantly faster.

- delta:

  `[0.001]`: Only used if `loss` is `"logit"` or `"exponential"`.
  Proportion of class membership is truncated to be within
  `[delta, 1-delta]` when calculating the loss to determine the optimal
  split. Should be positive for `"logit"`: with `delta = 0`, nodes
  containing only a single class produce an infinite loss and the
  corresponding splits are always rejected.

- epsilon:

  `[0.1]`: Only used if loss = `"logit"` or `"exponential"`. Proportion
  of class membership is truncated to be within `[epsilon, 1-epsilon]`
  when calculating the fit in a leaf. Unlike `delta` (a numerical guard
  for the split criterion), this caps the magnitude of individual leaf
  updates and acts as regularization: smaller values permit larger
  per-leaf jumps on the link scale.

- split_structure:

  `["leaves"]`: Defines the structure of a possible split and how to
  choose split_candidates. Can be one of "leaves", "hist",
  "cur_trees_1", "cur_trees_2", or "res_trees". Further details are
  given below.

- export_forest:

  `[FALSE]`: Whether to store the flattened forest in the returned
  object as `$forest`. If `FALSE`, `$forest` is `NULL`, reducing memory
  use of the returned object.

- formula:

  Formula specification, e.g. y ~ x1 + x2.

## Value

Object of class `"rpf"` with model object contained in `$fit`.

## Details

### splits

The number of `splits` is the main tuning parameter affecting the
accuracy of predictions.

### split_structure

The `split_structure` argument controls how split candidates are
constructed and sampled. In each round, a `t_try` fraction (capped by
`max_candidates`) is drawn from the pool of all possible splits with
weights `exp(-split_decay_rate * age)`.

- leaves:

  Split candidates are (leaf, split-dimension) pairs. For each sampled
  candidate, `split_try` thresholds are drawn uniformly from the valid
  range within that leaf and evaluated to choose the best split.

- cur_trees_1:

  Split candidates are (current-tree, split-dimension) pairs. For each
  sampled candidate, perform `split_try` evaluations. Each evaluation
  samples a leaf from the set of valid current trees (with probability
  proportional to its number of available thresholds) and then uniformly
  samples a single threshold within that leaf.

- cur_trees_2:

  Split candidates are (current-tree, split-dimension) pairs. For each
  sampled candidate, iterate through every valid leaf. Within each leaf,
  sample `split_try` thresholds uniformly and evaluate them.

- res_trees:

  Split candidates are resulting trees. For each sampled candidate, run
  `split_try` evaluations by sampling a (split-dimension, leaf) pair
  from all valid pairs (with probability proportional to its number of
  available thresholds), then uniformly sampling one threshold within
  that pair.

## Examples

``` r
# Regression with x and y
rpfit <- rpf(x = mtcars[, c("cyl", "wt")], y = mtcars$mpg)

# Regression with formula
rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
```
