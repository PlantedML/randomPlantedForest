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
#'   `data` instead of `y`). The predictors need unique names that do not
#'   contain `:`, which joins predictor names in component names such as `x1:x2`.
#'   `y` is the outcome vector: __numeric__ for regression, __factor__ for
#'   classification. A numeric outcome must be finite. A factor outcome needs at
#'   least two classes; levels that do not occur are dropped with a warning.
#' @param ... Not currently used, but required for extensibility. Unknown
#'   (e.g. misspelled) arguments are an error.
#' @param formula `[formula]`: A formula with the outcome on the left-hand side and the
#'   predictors on the right-hand side, e.g. `y ~ x1 + x2`. Interaction terms
#'   such as `x1:x2` or `x1 * x2` are an error: they would add a predictor named
#'   `x1:x2`, and `:` is how components are named. The forest finds interactions
#'   itself, up to `max_interaction`, so list the predictors as `y ~ x1 + x2`.
#'   To fit a product as a predictor of its own, write `I(x1 * x2)`, or add it to
#'   `data` as a column with its own name. Transformed terms of one predictor,
#'   such as `log(x1)`, are allowed.
#' @param data `[data.frame]`: Data containing the predictors and the outcome,
#'   used with `formula` or a recipe `x`.
#' @param max_interaction `[integer(1): 2]`: Maximum number of predictors a single tree may
#'   split on. `1` fits main effects only (an additive model), `2` adds
#'   pairwise interactions, and so on. `0` uses all predictors and values above
#'   the number of predictors are reduced to it.
#' @param ntrees `[integer(1): 50]`: Number of tree families, i.e. the size of the forest.
#'   Each family is grown on a bootstrap sample and predictions are averaged.
#' @param splits `[integer(1): 30]`: Number of split rounds per tree family. The main tuning
#'   parameter: more splits fit more complex functions. A round splits the best of its
#'   candidates; near saturation, when every candidate's leaf is too small to split, a
#'   round can make no split.
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
#' @param min_leaf_rows `[integer(1): 1]`: Minimum number of rows of the
#'   bootstrap sample on each side of every split, so in every leaf. Leaves
#'   without a threshold that keeps this many rows on each side are not drawn
#'   as split candidates. Only for regression with `split_structure = "leaves"`; `1`
#'   (no limit) otherwise.
#' @param interaction_constraints `[list | NULL: NULL]`: Feature sets that restrict which
#'   interactions the forest may form, e.g. `list(c("x1", "x2"), c("x3", "x4"))`.
#'   Each set is given by predictor names or by column indices: with `x`, the
#'   columns of `x`; with a formula, the predictors in the order the formula
#'   names them (`y ~ .`: the order of the columns of `data`); with a recipe,
#'   names only (its steps decide the columns). A tree on the predictors `T` is
#'   grown only if `T` is a subset of one of the sets, as with xgboost's
#'   interaction constraints. The subsets of a set stay allowed, so the forest
#'   keeps the main effects it grows a pair from. A predictor in no set is not
#'   used. Applies together with `max_interaction`: set `max_interaction = 0` to
#'   let the sets alone decide. `NULL` means no constraint. One use is to refit
#'   the forest on the effects a first fit found, e.g. after dropping its small
#'   components.
#' @param loss `[character(1) | NULL: NULL]`: Loss function. `NULL` uses `"L2"` for
#'   regression, the only supported loss, and `"exponential"` for classification.
#'   Classification also supports `"logit"`, which gives similar results but is
#'   slower, and `"L1"` and `"L2"`, which fit class indicators directly and do
#'   not yield proper probability estimates.
#' @param delta `[numeric(1): 0.001]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[delta, 1 - delta]` when computing the
#'   loss of a split. Must be less than 0.5, from where every proportion is
#'   truncated to the same value, and positive for `"logit"`: with `delta = 0`,
#'   nodes containing a single class have infinite loss and their splits are
#'   always rejected.
#' @param epsilon `[numeric(1): 0.1]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[epsilon, 1 - epsilon]` when computing
#'   the fit in a leaf. Unlike `delta`, this caps the size of individual leaf
#'   updates and acts as regularization: smaller values permit larger jumps on
#'   the link scale. Must be greater than 0 and less than 0.5: at 0 the leaf
#'   values of pure leaves are infinite, and from 0.5 on every leaf value is 0.
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
#'   (`t_try` and `split_decay_rate` are ignored) and evenly spaced thresholds. The result does not
#'   depend on the seed, and all tree families are identical, so use
#'   `ntrees = 1`. Mainly useful for testing and debugging.
#' @param keep_level_order `[logical(1) | character: FALSE]`: Which factor
#'   predictors keep the order of their levels instead of being ordered by the
#'   outcome. `FALSE` orders every factor by the outcome, `TRUE` keeps the level
#'   order of every factor, and a character vector keeps it for the named
#'   factor predictors only. See "Factor predictors" below.
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
#' ## Factor predictors
#'
#' Factor predictors, and character predictors, which are converted to
#' factors, are coded as integers `1, 2, ...` in an order of their levels, and
#' splits are thresholds on these codes. A split therefore separates the levels
#' below a cut from those above it in that order.
#'
#' By default, the levels are ordered by the outcome, so that the levels
#' most alike in the outcome are neighbours (Wright and König, 2019):
#' - regression: by the mean outcome of each level;
#' - binary classification: by the share of the second outcome level;
#' - more than two classes: by the first principal component of the class
#'   proportions of each level (Coppersmith et al., 1999).
#'
#' This applies to ordered factors too. With `keep_level_order`, the chosen
#' factors are instead coded in the order of `levels()`, as given, ordered or
#' not, so that every split separates neighbouring levels in that order.
#'
#' The coding is stored in the fit as `$factor_levels`, the levels of each
#' factor in the order of their codes. [predict()][predict.rpf()],
#' [predict_components()] and [purify()] code new data by these levels, by
#' name, so the order of the levels in `new_data` does not matter.
#'
#' ## `split_structure`
#'
#' In each round, a `t_try` fraction of all possible splits (capped at
#' `max_candidates`) is drawn as candidates with weights
#' `exp(-split_decay_rate * age)`. `split_structure` defines what a candidate is
#' and how its thresholds are evaluated.
#'
#' - `leaves`: Split candidates are (leaf, split-dimension) pairs. For each sampled
#'   candidate, `split_try` thresholds are drawn uniformly without replacement from the
#'   distinct values of that predictor among the leaf's rows, leaving out the smallest
#'   and, with `min_leaf_rows`, those that leave fewer rows on either side (all of them
#'   if there are fewer), and evaluated to choose the best split.
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
#' ## The draw of candidates for `leaves`
#'
#' For regression with `split_structure = "leaves"`, the global option
#' `randomPlantedForest.exact_draw` chooses how the candidates are drawn.
#' With `TRUE` (the default), every possible split gets a random number in
#' every round. With `FALSE`, the candidates are drawn from a Fenwick tree
#' over the weights, kept up to date between rounds, which costs a few random
#' numbers per candidate instead. Both draw from the same distribution, but
#' with other random numbers, so the same seed gives another (equally good)
#' forest. Ties between equally good candidates may go to different candidates
#' under the two draws, so with many splits, near saturation, their fits differ
#' in distribution, though not in accuracy. With `min_leaf_rows`, both leave
#' out the same candidates: those
#' whose leaf has no threshold with that many rows on each side. The fast draw
#' pays off on large data with many splits. The option is read when the forest
#' is fitted and stored in `params$exact_draw`; it does not affect
#' classification or the other split structures.
#'
#' ```r
#' options(randomPlantedForest.exact_draw = FALSE)
#' ```
#'
#' ## Reproducibility
#'
#' With the same seed, `rpf()` gives the same forest on the same platform, also
#' with more threads. Another platform can give another forest for the same
#' seed, because the integer sampling of the C++ standard library
#' (`std::uniform_int_distribution`) differs between implementations, such as
#' libc++ and libstdc++.
#'
#' @references
#' Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
#' directly interpretable tree ensemble. \doi{10.48550/arXiv.2012.14563}
#'
#' Wright, M. N. & König, I. R. (2019). Splitting on categorical predictors in
#' random forests. PeerJ 7:e6339. \doi{10.7717/peerj.6339}
#'
#' Coppersmith, D., Hong, S. J. & Hosking, J. R. (1999). Partitioning nominal
#' attributes in decision trees. Data Mining and Knowledge Discovery 3:197--217.
#' \doi{10.1023/A:1009869804967}
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
#'
#' # A factor that keeps its level order: splits separate neighbouring levels
#' cars <- mtcars
#' cars$gear <- factor(cars$gear, levels = c(3, 4, 5))
#' rpfit <- rpf(mpg ~ gear + wt, data = cars, keep_level_order = "gear")
#' rpfit$factor_levels
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
) {
  rlang::check_dots_empty()
  check_no_interaction_terms(formula, data)
  # model.matrix() would expand logical predictors into two dummy columns; use them
  # as 0/1 like the x/y interface. The outcome stays as is so it is still validated.
  outcome <- all.vars(rlang::f_lhs(formula))
  logical_predictors <- setdiff(names(data)[vapply(data, is.logical, logical(1))], outcome)
  for (col in logical_predictors) {
    data[[col]] <- as.integer(data[[col]])
  }
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  # column indices refer to the predictors in the formula's order (the model matrix puts factors last)
  interaction_constraints <- constraint_indices_to_names(
    interaction_constraints,
    labels = attr(stats::terms(formula, data = data), "term.labels"),
    columns = colnames(processed$predictors)
  )
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
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order
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
) {
  rlang::check_dots_empty()
  check_outcome_type(y)
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
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order
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
) {
  rlang::check_dots_empty()
  check_outcome_type(y)
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
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order
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
) {
  rlang::check_dots_empty()
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  # the recipe's steps decide the predictor columns: indices have no user-given order to refer to
  if (is.list(interaction_constraints) && any(vapply(interaction_constraints, is.numeric, logical(1)))) {
    cli::cli_abort("With a recipe, give {.arg interaction_constraints} by predictor names, not column indices.")
  }
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
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order
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
  min_leaf_rows,
  interaction_constraints,
  loss,
  delta,
  epsilon,
  purify,
  nthreads,
  export_forest,
  deterministic,
  keep_level_order
) {
  # report errors from helpers against the user's rpf() call
  call <- rlang::caller_env()
  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  check_no_missing(processed$predictors, "The predictors", call = call)
  check_no_missing(processed$outcomes, "The outcome", call = call)
  check_predictor_names(names(processed$predictors), call = call)
  processed <- check_outcome(processed, call = call)
  predictors <- preprocess_predictors_fit(processed, keep_level_order = keep_level_order, call = call)
  checkmate::assert_string(loss, null.ok = TRUE)
  outcomes <- preprocess_outcome(processed, loss, call = call)
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
  if (t_try <= 0) {
    cli::cli_abort("{.arg t_try} must be greater than 0, not {t_try}.", call = call)
  }
  checkmate::assert_int(max_candidates, lower = 1)
  checkmate::assert_number(split_decay_rate, lower = 0)
  checkmate::assert_flag(delete_leaves)
  checkmate::assert_int(min_leaf_rows, lower = 1)
  if (min_leaf_rows > 1 && (outcomes$mode != "regression" || split_structure != "leaves")) {
    cli::cli_abort(
      "{.arg min_leaf_rows} > 1 needs regression with {.code split_structure = \"leaves\"}.",
      call = call
    )
  }

  # "median" loss is implemented but discarded
  loss_functions <- switch(outcomes$mode, "regression" = "L2", "classification" = c("L1", "L2", "logit", "exponential"))
  checkmate::assert_choice(loss, choices = loss_functions)
  checkmate::assert_number(delta, lower = 0, upper = 1)
  checkmate::assert_number(epsilon, lower = 0, upper = 1)
  if (loss %in% c("logit", "exponential")) {
    # Both truncate class proportions to [value, 1 - value], which is empty from 0.5 on
    if (delta >= 0.5 || (loss == "logit" && delta <= 0)) {
      range <- if (loss == "logit") "greater than 0" else "at least 0"
      cli::cli_abort(
        "With {.code loss = \"{loss}\"}, {.arg delta} must be {range} and less than 0.5, not {delta}.",
        call = call
      )
    }
    if (epsilon <= 0 || epsilon >= 0.5) {
      cli::cli_abort(
        "With {.code loss = \"{loss}\"}, {.arg epsilon} must be greater than 0 and less than 0.5, not {epsilon}.",
        call = call
      )
    }
  }

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

  interaction_constraints <- check_interaction_constraints(
    interaction_constraints,
    colnames(predictors$predictors_matrix),
    max_interaction,
    call = call
  )

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
    min_leaf_rows = as.integer(min_leaf_rows),
    interaction_constraints = interaction_constraints,
    split_structure = split_structure,
    delta = delta,
    epsilon = epsilon,
    deterministic = deterministic,
    nthreads = nthreads,
    purify = purify,
    exact_draw = rpf_exact_draw(call = call),
    keep_level_order = predictors$keep_level_order
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

# The draw of candidates for split_structure = "leaves": the option
# randomPlantedForest.exact_draw (TRUE, the default: the exact draw; FALSE: the
# fast draw, the same distribution with other random numbers).
rpf_exact_draw <- function(call = rlang::caller_env()) {
  exact_draw <- getOption("randomPlantedForest.exact_draw", TRUE)
  if (!rlang::is_bool(exact_draw)) {
    cli::cli_abort(
      "The option {.code randomPlantedForest.exact_draw} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }
  exact_draw
}

# Assemble the positional parameter vector for the C++ constructors.
# 15 elements for regression (+ min_leaf_rows, exact_draw), 15 (+ delta,
# epsilon) for classification, which has neither. Models saved before these
# options existed have no such entries in their params and were fitted with
# min_leaf_rows = 1 and the exact draw.
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
  } else {
    base <- c(base, params$min_leaf_rows %||% 1L, params$exact_draw %||% TRUE)
  }
  base
}

# Main fitting function and interface to C++ implementation
rpf_impl <- function(Y, X, mode = c("regression", "classification"), params) {
  # Final input validation, should be superfluous
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
  mode <- match.arg(mode)
  pars <- rpf_param_vector(params, mode)

  if (is.null(params$interaction_constraints)) {
    if (mode == "classification") {
      fit <- new(ClassificationRPF, Y, X, params$loss, pars)
    } else {
      fit <- new(RandomPlantedForest, Y, X, pars)
    }
  } else {
    # allowed feature sets as 1-based column indices
    constraints <- lapply(params$interaction_constraints, match, table = colnames(X))
    if (mode == "classification") {
      fit <- new(ClassificationRPF, Y, X, params$loss, pars, constraints)
    } else {
      fit <- new(RandomPlantedForest, Y, X, pars, constraints)
    }
  }

  fit
}

# Formula interface: translate the sets given by column indices into predictor names, index i referring to the
# i-th predictor of the formula (`labels`, its term labels; `y ~ .` gives the data's column order). The model matrix
# (`columns`) orders the predictors differently: hardhat puts factors last. Sets of names are left as they are;
# check_interaction_constraints() validates everything else.
constraint_indices_to_names <- function(interaction_constraints, labels, columns, call = rlang::caller_env()) {
  if (!is.list(interaction_constraints)) {
    return(interaction_constraints)
  }
  lapply(seq_along(interaction_constraints), function(i) {
    s <- interaction_constraints[[i]]
    if (!is.numeric(s) || length(s) == 0 || anyNA(s)) {
      return(s)
    }
    if (any(s != round(s)) || any(s < 1) || any(s > length(labels))) {
      cli::cli_abort(
        paste(
          "{.code interaction_constraints[[{i}]]} must contain indices between 1 and {length(labels)}",
          "of the formula's predictors (or predictor names)."
        ),
        call = call
      )
    }
    nm <- labels[s]
    if (!all(nm %in% columns)) {
      bad <- nm[!nm %in% columns]
      cli::cli_abort(
        c(
          "{.code interaction_constraints[[{i}]]}: the formula's predictor{?s} {.val {bad}} {?is/are} not a single model column.",
          "i" = "Give the sets by column names: {.val {columns}}."
        ),
        call = call
      )
    }
    nm
  })
}

# Validate `interaction_constraints` and convert it to a list of sorted, unique predictor
# names (NULL stays NULL: no constraint).
check_interaction_constraints <- function(interaction_constraints, feature_names, max_interaction, call) {
  if (is.null(interaction_constraints)) {
    return(NULL)
  }
  p <- length(feature_names)
  if (!is.list(interaction_constraints) || length(interaction_constraints) == 0) {
    cli::cli_abort(
      c(
        "{.arg interaction_constraints} must be {.code NULL} or a non-empty list of feature sets.",
        "i" = "For example {.code list(c(\"x1\", \"x2\"), \"x3\")}."
      ),
      call = call
    )
  }
  sets <- lapply(seq_along(interaction_constraints), function(i) {
    s <- interaction_constraints[[i]]
    if (length(s) == 0 || anyNA(s)) {
      cli::cli_abort(
        "{.code interaction_constraints[[{i}]]} must be a non-empty vector without {.code NA}.",
        call = call
      )
    }
    if (is.character(s)) {
      unknown <- setdiff(s, feature_names)
      if (length(unknown) > 0) {
        cli::cli_abort(
          c(
            "{.code interaction_constraints[[{i}]]} names unknown predictor{?s} {.val {unknown}}.",
            "i" = "The predictors are {.val {feature_names}}."
          ),
          call = call
        )
      }
      idx <- match(s, feature_names)
    } else if (is.numeric(s)) {
      if (any(s != round(s)) || any(s < 1) || any(s > p)) {
        cli::cli_abort(
          "{.code interaction_constraints[[{i}]]} must contain column indices between 1 and {p} (or predictor names).",
          call = call
        )
      }
      idx <- as.integer(s)
    } else {
      cli::cli_abort(
        "{.code interaction_constraints[[{i}]]} must be predictor names or column indices, not {.obj_type_friendly {s}}.",
        call = call
      )
    }
    feature_names[sort(unique(idx))]
  })
  if (max(lengths(sets)) > max_interaction) {
    cli::cli_inform(c(
      "i" = paste(
        "{.arg interaction_constraints} has a set of {max(lengths(sets))} predictors,",
        "but {.arg max_interaction} is {max_interaction}."
      ),
      " " = "Trees use at most {max_interaction} predictor{?s}."
    ))
  }
  sets
}
