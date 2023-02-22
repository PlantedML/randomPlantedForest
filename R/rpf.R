#' Random Planted Forest
#'
#' @param x,data Feature `matrix`, or `data.frame`, or [`recipe`][recipes::recipe].
#' @param y Target vector for use with `x`.
#'   The class of `y` (either `numeric` or [`factor`]) determines if regression
#'   or classification will be performed.
#' @param formula Formula specification, e.g. y ~ x1 + x2.
#' @param max_interaction `[1]`: Maximum level of interaction determining maximum
#'   number of split dimensions for a tree.
#'   The default `1` corresponds to main effects only.
#'   If `0`, the number fo columns in `x` is used, i.e. for 10 predictors,
#'   this is equivalent to setting `max_interaction = 10`.
#' @param ntrees `[50]`: Number of trees generated per family.
#' @param splits `[30]`: Number of splits performed for each tree family.
#' @param split_try `[10]`: Number of split points to be considered when choosing a split candidate.
#' @param t_try `[0.4]`: A value in (0,1] specifying the proportion of viable split-candidates in each round.
#' @param deterministic `[FALSE]`: Choose whether approach deterministic or random.
#' @param parallel `[FALSE]`: Perform algorithm in parallel or serialized.
#' @param purify `[FALSE]`: Whether the forest should be purified.
#'   Set to `TRUE` to enable components extract with [`extract_components()`] are valid.
#'   Can be achieved after fitting with [`purify()`].
#' @param cv `[FALSE]`: Determines if cross validation is performed.
#' @param loss `["L2"]`: For regression, only `"L2"` is supported. For
#'   classification, `"L1"`, `"logit"` and "`exponential`" are also available.
#'   "`exponential`" yield similar results as "`logit`" while being significantly faster.
#' @param delta `[0]`: Only used if loss = `"logit"` or `"exponential"`.
#'   Proportion of class membership is truncated to be smaller 1-delta when calculating
#'   the loss to determine the optimal split.
#' @param epsilon `[0.1]`: Only used if loss = `"logit"` or `"exponential"`.
#'   Proportion of class membership is truncated to be smaller 1-epsilon when calculating
#'   the fit in a leaf.
#' @param ... (Unused).
#'
#' @return Object of class `"rpf"` with model object contained in `$fit`.
#' @export
#' @importFrom methods new
#' @importFrom hardhat mold
#' @importFrom hardhat default_xy_blueprint
#' @importFrom hardhat default_formula_blueprint
#' @importFrom hardhat default_recipe_blueprint
#'
#' @examples
#' # Regression with x and y
#' rpfit <- rpf(x = mtcars[, c("cyl", "wt")], y = mtcars$mpg)
#'
#' # Regression with formula
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
rpf <- function(x, ...) {
  UseMethod("rpf")
}

#' @export
#' @noRd
# @rdname rpf
rpf.default <- function(x, ...) {
  stop(
    "`rpf()` is not defined for a '", class(x)[1], "'.",
    call. = FALSE
  )
}

# XY method - data frame
#' @export
#' @rdname rpf
rpf.data.frame <- function(x, y, max_interaction = 1, ntrees = 50, splits = 30,
                           split_try = 10, t_try = 0.4, deterministic = FALSE,
                           parallel = FALSE, purify = FALSE, cv = FALSE,
                           loss = "L2", delta = 0, epsilon = 0.1, ...) {
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed, max_interaction, ntrees, splits,
    split_try, t_try, deterministic,
    parallel, purify, cv,
    loss, delta, epsilon
  )
}

# XY method - matrix
#' @export
#' @rdname rpf
rpf.matrix <- function(x, y, max_interaction = 1, ntrees = 50, splits = 30,
                       split_try = 10, t_try = 0.4, deterministic = FALSE,
                       parallel = FALSE, purify = FALSE, cv = FALSE,
                       loss = "L2", delta = 0, epsilon = 0.1, ...) {
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed, max_interaction, ntrees, splits,
    split_try, t_try, deterministic,
    parallel, purify, cv,
    loss, delta, epsilon
  )}

# Formula method
#' @export
#' @rdname rpf
rpf.formula <- function(formula, data, max_interaction = 1, ntrees = 50, splits = 30,
                        split_try = 10, t_try = 0.4, deterministic = FALSE,
                        parallel = FALSE, purify = FALSE, cv = FALSE,
                        loss = "L2", delta = 0, epsilon = 0.1, ...) {
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  rpf_bridge(
    processed, max_interaction, ntrees, splits,
    split_try, t_try, deterministic,
    parallel, purify, cv,
    loss, delta, epsilon
  )
}

# Recipe method
#' @export
#' @rdname rpf
rpf.recipe <- function(x, data, max_interaction = 1, ntrees = 50, splits = 30,
                       split_try = 10, t_try = 0.4, deterministic = FALSE,
                       parallel = FALSE, purify = FALSE, cv = FALSE,
                       loss = "L2", delta = 0, epsilon = 0.1, ...) {
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  rpf_bridge(
    processed, max_interaction, ntrees, splits,
    split_try, t_try, deterministic,
    parallel, purify, cv,
    loss, delta, epsilon
  )
}

# Bridge: Calls rpf_impl() with processed input
#
#' @param processed Output of `hardhat::mold` from respective rpf methods
#' @importFrom hardhat validate_outcomes_are_univariate
rpf_bridge <- function(processed, max_interaction = 1, ntrees = 50, splits = 30,
                       split_try = 10, t_try = 0.4, deterministic = FALSE,
                       parallel = FALSE, purify = FALSE, cv = FALSE,
                       loss = "L2", delta = 0, epsilon = 0.1) {
  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  predictors <- preprocess_predictors_fit(processed)
  outcomes <- preprocess_outcome(processed, loss)
  p <- ncol(predictors$predictors_matrix)

  # Check arguments
  checkmate::assert_int(max_interaction, lower = 0)

  # rewrite max_interaction so 0 -> "maximum", e.g. ncol(X)
  if (max_interaction == 0) {
    max_interaction <- p
  }
  # same applies to values > p
  if (max_interaction > p) {
    message("`max_interaction` set to ", max_interaction, " but ncol(X) is ", p, ".\n",
            "Setting `max_interaction` to ", p)
    max_interaction <- p
  }

  checkmate::assert_int(ntrees, lower = 1)
  checkmate::assert_int(splits, lower = 1)
  checkmate::assert_int(split_try, lower = 1)

  checkmate::assert_number(t_try, lower = 0, upper = 1)
  checkmate::assert_number(delta, lower = 0, upper = 1)
  checkmate::assert_number(epsilon, lower = 0, upper = 1)

  # "median" is implemented but discarded
  checkmate::assert_choice(
    loss,
    choices = c("L1", "L2", "logit", "exponential")
  )

  checkmate::assert_flag(deterministic)
  checkmate::assert_flag(parallel)
  checkmate::assert_flag(purify)
  checkmate::assert_flag(cv)

  fit <- rpf_impl(
    Y = outcomes$outcomes, X = predictors$predictors_matrix,
    mode = outcomes$mode,
    max_interaction = max_interaction, ntrees = ntrees, splits = splits,
    split_try = split_try, t_try = t_try, deterministic = deterministic,
    parallel = parallel, purify = purify, cv = cv,
    loss = loss, delta = delta, epsilon = epsilon
  )

  forest <- fit$get_model()
  class(forest) <- "rpf_forest"

  new_rpf(
    fit = fit,
    blueprint = processed$blueprint,
    mode = outcomes$mode,
    factor_levels = predictors$factor_levels,
    params = list(
      loss = loss, # FIXME: Dedup, requires changes in tests and predict
      ntrees = ntrees,
      max_interaction = max_interaction,
      splits = splits,
      split_try = split_try, t_try = t_try,
      delta = delta, epsilon = epsilon,
      deterministic = deterministic,
      parallel = parallel, purify = purify, cv = cv
    ),
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

# Main fitting function and interface to C++ implementation
rpf_impl <- function(Y, X, mode = c("regression", "classification"),
                     max_interaction = 1, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4,
                     deterministic = FALSE, parallel = FALSE, purify = FALSE, cv = FALSE,
                     loss = "L2", delta = 0, epsilon = 0.1) {
  # Final input validation, should be superfluous
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
  mode <- match.arg(mode)

  if (mode == "classification") {
    fit <- new(ClassificationRPF, Y, X, loss, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv, delta, epsilon
    ))
  } else if (mode == "regression") {
    fit <- new(RandomPlantedForest, Y, X, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv
    ))
  }

  fit
}
