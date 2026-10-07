#' Random Planted Forest
#'
#' Fits a random planted forest (Hiabu, Mammen and Meyer, 2020) for regression
#' or classification. The model is a sum of trees that each split on at most
#' `max_interaction` predictors, so it decomposes into main effects and
#' interaction components; see [predict_components()].
#'
#' @param x,y `[data.frame | matrix | recipe]`, `[numeric | factor]`:
#'   Predictors and outcome for the x/y interface.
#'   `x` holds the predictors, or is a [`recipe`][recipes::recipe] (then supply
#'   `data` instead of `y`).
#'   `y` is the outcome vector: __numeric__ for regression, __factor__ for
#'   classification.
#' @param ... Not currently used, but required for extensibility. Unknown
#'   (e.g. misspelled) arguments are an error.
#' @param formula `[formula]`: A formula with the outcome on the left-hand side and the
#'   predictors on the right-hand side, e.g. `y ~ x1 + x2`.
#' @param data `[data.frame]`: Data containing the predictors and the outcome,
#'   used with `formula` or a recipe `x`.
#' @param max_interaction `[integer(1): 2]`: Maximum number of predictors a single tree may
#'   split on. `1` fits main effects only (an additive model), `2` adds
#'   pairwise interactions, and so on. `0` uses all predictors and values above
#'   the number of predictors are reduced to it.
#' @param ntrees `[integer(1): 50]`: Number of tree families, i.e. the size of the forest.
#'   Each family is grown on a bootstrap sample and predictions are averaged.
#' @param splits `[integer(1): 30]`: Number of splits per tree family. The main tuning
#'   parameter: more splits fit more complex functions.
#' @param split_structure `[character(1): "leaves"]`: How split candidates are formed and
#'   sampled; one of `"leaves"`, `"hist"`, `"cur_trees_1"`, `"cur_trees_2"`,
#'   or `"res_trees"`. See Details.
#' @param split_try `[integer(1): 10]`: Number of thresholds evaluated per split candidate.
#' @param t_try `[numeric(1): 0.4]`: Proportion in `(0, 1]` of all possible splits sampled
#'   as candidates in each round, capped at `max_candidates`.
#' @param max_candidates `[integer(1): 50]`: Maximum number of split candidates per round.
#' @param split_decay_rate `[numeric(1): 0.1]`: Down-weights possible splits that were
#'   sampled as candidates but not chosen. Each such round ages a split by one,
#'   choosing it resets its age to zero, and splits are sampled with weight
#'   `exp(-split_decay_rate * age)`. `0` samples uniformly.
#' @param delete_leaves `[logical(1): TRUE]`: Whether to delete a parent leaf when
#'   splitting along an existing dimension.
#' @param loss `[character(1) | NULL: NULL]`: Loss function. `NULL` uses `"L2"` for
#'   regression, the only supported loss, and `"exponential"` for classification.
#'   Classification also supports `"logit"`, which gives similar results but is
#'   slower, and `"L1"` and `"L2"`, which fit class indicators directly and do
#'   not yield proper probability estimates.
#' @param delta `[numeric(1): 0.001]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[delta, 1 - delta]` when computing the
#'   loss of a split. Should be positive for `"logit"`: with `delta = 0`, nodes
#'   containing a single class have infinite loss and their splits are always
#'   rejected.
#' @param epsilon `[numeric(1): 0.1]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[epsilon, 1 - epsilon]` when computing
#'   the fit in a leaf. Unlike `delta`, this caps the size of individual leaf
#'   updates and acts as regularization: smaller values permit larger jumps on
#'   the link scale.
#' @param purify `[logical(1): FALSE]`: Whether to purify the forest after fitting, which
#'   [predict_components()] requires, with the defaults of [purify()]. For other
#'   settings, call [purify()] after fitting instead.
#' @param nthreads `[integer(1): 1]`: Number of threads for fitting. Also the default for
#'   [predict()][predict.rpf()] and [purify()].
#' @param export_forest `[logical(1): FALSE]`: Whether to store the flattened forest in
#'   the returned object as `$forest`. If `FALSE`, `$forest` is `NULL`, which
#'   saves memory.
#' @param deterministic `[logical(1): FALSE]`: Whether to fit without randomness: no
#'   bootstrap, the first `max_candidates` possible splits as candidates
#'   (`t_try` is ignored) and evenly spaced thresholds. The result does not
#'   depend on the seed, and all tree families are identical, so use
#'   `ntrees = 1`. Mainly useful for testing and debugging.
#'
#' @return Object of class `"rpf"` with model object contained in `$fit`.
#' @export
#' @importFrom methods new
#' @importFrom hardhat mold
#' @importFrom hardhat default_xy_blueprint
#' @importFrom hardhat default_formula_blueprint
#' @importFrom hardhat default_recipe_blueprint
#'
#' @details
#' ## Choosing parameters
#'
#' Start with `max_interaction` and `splits`. `max_interaction` sets which
#' effects the model can represent: `1` for an additive model, `2` (the
#' default) to add pairwise interactions, which can still be plotted as
#' heatmaps. Higher values are more flexible, but slower and the components are
#' harder to interpret. `splits` sets how closely the model
#' follows the data and is the parameter to tune. With a higher
#' `max_interaction`, splits are spread over more possible components, so tune
#' both together.
#'
#' `ntrees` averages over bootstrap samples to reduce variance and fitting and
#' prediction time grow linearly with it. The default is lower than the
#' usual 500 trees of a random forest because each tree family is itself a
#' complete model, a sum of trees covering all components, rather than a
#' single tree. Averaging therefore stops paying off sooner, usually after a
#' few dozen families. The split search parameters
#' (`split_structure`, `split_try`, `t_try`, `max_candidates`,
#' `split_decay_rate`) trade speed for a more thorough search.
#'
#' For classification, `"exponential"` (the default) and `"logit"` fit on the
#' link scale and are transformed to probabilities, `"exponential"` being faster.
#' `"L1"` and `"L2"` fit class indicators directly and their probabilities are
#' truncated to `[0, 1]`.
#'
#' ## `split_structure`
#'
#' In each round, a `t_try` fraction of all possible splits (capped at
#' `max_candidates`) is drawn as candidates with weights
#' `exp(-split_decay_rate * age)`. `split_structure` defines what a candidate is
#' and how its thresholds are evaluated.
#'
#' - `leaves`: Split candidates are (leaf, split-dimension) pairs. For each sampled
#'   candidate, `split_try` thresholds are drawn uniformly from the valid range within
#'   that leaf and evaluated to choose the best split.
#'
#' - `hist`: As `leaves`, but thresholds are drawn from boundaries of
#'   quantile bins computed once per predictor, which makes evaluating them
#'   faster on large data.
#'
#' - `cur_trees_1`: Split candidates are (current-tree, split-dimension) pairs. For each
#'   sampled candidate, perform `split_try` evaluations. Each evaluation samples a leaf
#'   from the set of valid current trees (with probability proportional to its number of
#'   available thresholds) and then uniformly samples a single threshold within that leaf.
#'
#' - `cur_trees_2`: Split candidates are (current-tree, split-dimension) pairs. For each
#'   sampled candidate, iterate through every
#'   valid leaf. Within each leaf, sample `split_try` thresholds uniformly and
#'   evaluate them.
#'
#' - `res_trees`: Split candidates are resulting trees. For each sampled candidate, run
#'   `split_try` evaluations by sampling a (split-dimension, leaf) pair from all valid
#'   pairs (with probability proportional to its number of available thresholds), then
#'   uniformly sampling one threshold within that pair.
#'
#' @references
#' Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
#' directly interpretable tree ensemble. \doi{10.48550/arXiv.2012.14563}
#'
#' @examples
#' # Regression with formula
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
#'
#' # Regression with x and y, allowing pairwise interactions
#' rpfit <- rpf(x = mtcars[, c("cyl", "wt", "hp")], y = mtcars$mpg, max_interaction = 2)
#'
#' # Classification
#' rpfit <- rpf(Species ~ ., data = iris, loss = "exponential")
rpf <- function(x, ...) {
  UseMethod("rpf")
}

#' @export
#' @noRd
rpf.default <- function(x, ...) {
  cli::cli_abort("{.fn rpf} is not defined for a {.cls {class(x)[1]}}.")
}

# Formula method
#' @export
#' @rdname rpf
rpf.formula <- function(
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
) {
  rlang::check_dots_empty()
  # model.matrix() would expand logical predictors into two dummy columns; use them
  # as 0/1 like the x/y interface. The outcome stays as is so it is still validated.
  outcome <- all.vars(rlang::f_lhs(formula))
  logical_predictors <- setdiff(names(data)[vapply(data, is.logical, logical(1))], outcome)
  for (col in logical_predictors) {
    data[[col]] <- as.integer(data[[col]])
  }
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic
  )
}

# XY method - data frame
#' @export
#' @rdname rpf
rpf.data.frame <- function(
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
) {
  rlang::check_dots_empty()
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic
  )
}

# XY method - matrix
#' @export
#' @rdname rpf
rpf.matrix <- function(
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
) {
  rlang::check_dots_empty()
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic
  )
}

# Recipe method
#' @export
#' @rdname rpf
rpf.recipe <- function(
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
) {
  rlang::check_dots_empty()
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic
  )
}

# Bridge: validates arguments and calls rpf_impl() with processed input
#' @noRd
#' @param processed `[list]`: Output of `hardhat::mold` from respective rpf methods
#' @importFrom hardhat validate_outcomes_are_univariate
rpf_bridge <- function(
  processed,
  max_interaction,
  ntrees,
  splits,
  split_structure,
  split_try,
  t_try,
  max_candidates,
  split_decay_rate,
  delete_leaves,
  loss,
  delta,
  epsilon,
  purify,
  nthreads,
  export_forest,
  deterministic
) {
  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  check_no_missing(processed$predictors, "The predictors", call = call)
  check_no_missing(processed$outcomes, "The outcome", call = call)
  predictors <- preprocess_predictors_fit(processed)
  checkmate::assert_string(loss, null.ok = TRUE)
  outcomes <- preprocess_outcome(processed, loss)
  loss <- outcomes$loss
  p <- ncol(predictors$predictors_matrix)

  # Check arguments
  checkmate::assert_int(max_interaction, lower = 0)

  # rewrite max_interaction so 0 -> "maximum", e.g. ncol(X)
  if (max_interaction == 0) {
    max_interaction <- p
  }
  # same applies to values > p, silently, as the default of 2 exceeds p for a single predictor
  if (max_interaction > p) {
    max_interaction <- p
  }

  checkmate::assert_int(ntrees, lower = 1)
  checkmate::assert_int(splits, lower = 1)
  checkmate::assert_choice(split_structure, choices = c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees"))
  checkmate::assert_int(split_try, lower = 1)
  checkmate::assert_number(t_try, lower = 0, upper = 1)
  checkmate::assert_int(max_candidates, lower = 1)
  checkmate::assert_number(split_decay_rate, lower = 0)
  checkmate::assert_flag(delete_leaves)

  # "median" loss is implemented but discarded
  loss_functions <- switch(outcomes$mode, "regression" = "L2", "classification" = c("L1", "L2", "logit", "exponential"))
  checkmate::assert_choice(loss, choices = loss_functions)
  checkmate::assert_number(delta, lower = 0, upper = 1)
  checkmate::assert_number(epsilon, lower = 0, upper = 1)

  checkmate::assert_flag(purify)
  checkmate::assert_int(nthreads, lower = 1L)
  checkmate::assert_flag(export_forest)
  checkmate::assert_flag(deterministic)
  if (deterministic && ntrees > 1) {
    cli::cli_warn(c(
      "With {.code deterministic = TRUE}, all {ntrees} tree families are identical.",
      "i" = "Use {.code ntrees = 1} to fit the same model faster."
    ))
  }

  params <- list(
    loss = loss,
    ntrees = ntrees,
    max_interaction = max_interaction,
    splits = splits,
    split_try = split_try,
    t_try = t_try,
    split_decay_rate = split_decay_rate,
    max_candidates = max_candidates,
    delete_leaves = delete_leaves,
    split_structure = split_structure,
    delta = delta,
    epsilon = epsilon,
    deterministic = deterministic,
    nthreads = nthreads,
    purify = purify
  )

  fit <- rpf_impl(
    Y = outcomes$outcomes,
    X = predictors$predictors_matrix,
    mode = outcomes$mode,
    params = params
  )

  # Optionally export a compact R list representation of the forest.
  forest <- NULL
  if (export_forest) {
    forest <- fit$get_model()
    class(forest) <- "rpf_forest"
  }

  new_rpf(
    fit = fit,
    blueprint = processed$blueprint,
    mode = outcomes$mode,
    factor_levels = predictors$factor_levels,
    params = params,
    forest = forest
  )
}

# Intermediate to hold model object with blueprint used for prediction
new_rpf <- function(fit, blueprint, ...) {
  hardhat::new_model(
    fit = fit,
    blueprint = blueprint,
    class = "rpf",
    ...
  )
}

# Assemble the positional parameter vector for the C++ constructors.
# 13 elements for regression, 15 (+ delta, epsilon) for classification.
rpf_param_vector <- function(params, mode) {
  split_mode <- switch(
    params$split_structure,
    res_trees = 0L,
    cur_trees_2 = 1L,
    cur_trees_1 = 2L,
    leaves = 3L,
    hist = 4L
  )
  base <- c(
    params$max_interaction,
    params$ntrees,
    params$splits,
    params$split_try,
    params$t_try,
    params$purify,
    params$deterministic,
    params$nthreads,
    0, # cross-validation slot, unused by the C++ core
    params$split_decay_rate,
    params$max_candidates,
    params$delete_leaves,
    split_mode
  )
  if (mode == "classification") {
    base <- c(base, params$delta, params$epsilon)
  }
  base
}

# Main fitting function and interface to C++ implementation
rpf_impl <- function(Y, X, mode = c("regression", "classification"), params) {
  # Final input validation, should be superfluous
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
  mode <- match.arg(mode)
  pars <- rpf_param_vector(params, mode)

  if (mode == "classification") {
    fit <- new(ClassificationRPF, Y, X, params$loss, pars)
  } else {
    fit <- new(RandomPlantedForest, Y, X, pars)
  }

  fit
}
