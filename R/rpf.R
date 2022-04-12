#' Random Planted Forest
#'
# FIXME: Parameters need describing
#' @param x Feature matrix or `data.frame`.
#' @param y Target vector for use with `x`.
#'   The class of `y` (either `numeric` or [`factor`]) determines if regression
#'   or classification will be performed.
#' @param formula Formula specification, e.g. y ~ x1 + x2.
#' @param data A `data.frame` for use with `formula`.
#' @param max_interaction `[1]`: Maximum level of interaction determining maximum
#'   number of split dimensions for a tree.
#' @param ntrees `[50]`: Number of trees generated per family.
#' @param splits `[30]`: Number of performed splits for each tree family.
#' @param split_try `[10]`: Number of split points to be considered when considering a split candidate.
#' @param t_try `[0.4]`: A value in (0,1] specifying the proportion of split-candidetes viable in each round. 
#' @param deterministic `[FALSE]`: Choose whether approach deterministic or random.
#' @param parallel `[FALSE]`: Perform algorithm in parallel or serialized.
#' @param purify `[FALSE]`: Whether the forest should be purified.
#' @param cv `[FALSE]`: Determines if cross validation is performed.
#' @param loss `["L2"]`: For regression, only `"L2"` is supported. For
#'   classification, `"L1"`, `"logit"` and "`exponential`" are also available.
#' @param delta `[0]`: Only used if loss = `"logit"` or `"exponential"`. Proportion of class membership is truncated to be smaller 1-delta when calculating the loss to determin the optimal split. 
#' @param epsilon `[0.1]`: Only used if loss = `"logit"` or `"exponential"`. Proportion of class membership is truncated to be smaller 1-epsilon when calculating the fit in a leave. 
#' @param ... (Ignored).
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
rpf <- function(x, ..., max_interaction = 1, ntrees = 50, splits = 30,
                split_try = 10, t_try = 0.4, deterministic = FALSE,
                parallel = FALSE, purify = FALSE, cv = FALSE,
                loss = "L2", delta = 0, epsilon = 0.1) {
  UseMethod("rpf")
}

#' @export
#' @rdname rpf
rpf.default <- function(x, ...) {
  stop(
    "`rpf()` is not defined for a '", class(x)[1], "'.",
    call. = FALSE
  )
}

# XY method - data frame
#' @export
#' @rdname rpf
rpf.data.frame <- function(x, y, ...) {
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# XY method - matrix
#' @export
#' @rdname rpf
rpf.matrix <- function(x, y, ...) {
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# Formula method
#' @export
#' @rdname rpf
rpf.formula <- function(formula, data, ...) {
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# Recipe method
#' @export
#' @rdname rpf
rpf.recipe <- function(x, data, ...) {
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# Bridge: Calls rpf_impl() with processed input
#' @importFrom hardhat validate_outcomes_are_univariate
rpf_bridge <- function(processed, max_interaction = 1, ntrees = 50, splits = 30,
                       split_try = 10, t_try = 0.4, deterministic = FALSE,
                       parallel = FALSE, purify = FALSE, cv = FALSE,
                       loss = "L2", delta = 0, epsilon = 0.1) {
  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  predictors <- preprocess_predictors_fit(processed)
  outcomes <- preprocess_outcome(processed, loss)

  # FIXME: loss function handling for multiclass is a clunky hack to ensure
  # the user only sees e.g. "logit" as an option
  if (ncol(outcomes$outcomes) > 1) {
    loss <- switch (loss,
      "logit" = "logit_2",
      "exponential" = "exponential_2",
      loss
    )
  }

  # Check arguments
  checkmate::assert_integerish(max_interaction, lower = 1, len = 1)
  checkmate::assert_integerish(ntrees, lower = 1, len = 1)
  checkmate::assert_integerish(splits, lower = 1, len = 1)
  checkmate::assert_integerish(split_try, lower = 1, len = 1)

  checkmate::assert_numeric(t_try, lower = 0, upper = 1, len = 1)
  # FIXME: What is delta/epsilon and what can it look like?
  checkmate::assert_numeric(delta, lower = 0, upper = 1, len = 1)
  checkmate::assert_numeric(epsilon, lower = 0, upper = 1, len = 1)

  # "median" is implemented but discarded
  checkmate::assert_choice(
    loss,
    choices = c("L1", "L2", "logit", "logit_2", "exponential", "exponential_2"), 
    null.ok = FALSE
  )

  checkmate::assert_logical(deterministic, len = 1)
  checkmate::assert_logical(parallel, len = 1)
  checkmate::assert_logical(purify, len = 1)
  checkmate::assert_logical(cv, len = 1)

  fit <- rpf_impl(
    Y = outcomes$outcomes, X = predictors$predictors_matrix,
    mode = outcomes$mode,
    max_interaction = max_interaction, ntrees = ntrees, splits = splits,
    split_try = split_try, t_try = t_try, deterministic = deterministic,
    parallel = parallel, purify = purify, cv = cv,
    loss = loss, delta = delta, epsilon = epsilon
  )

  new_rpf(
    fit = fit,
    blueprint = processed$blueprint,
    mode = outcomes$mode,
    factor_levels = predictors$factor_levels,
    loss = loss
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
  # checkmate::assert_matrix(Y, mode = "numeric", any.missing = FALSE)

  if (mode == "classification") {
    # FIXME: Handling for classification modes, must allow 1/0 or 1/-1 and 
    # be a matrix
    #checkmate::assert_integer(Y, lower = 0)

    fit <- new(ClassificationRPF, Y, X, loss, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv, delta, epsilon
    ))
  } else if (mode == "regression") {
    # FIXME: Loss missing here?
    # Passing loss as arg gives error
    # "no valid constructor available for the argument list"
    # N.B. Neither delta nor epsilon are passed as well
    fit <- new(RandomPlantedForest, Y, X, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv
    ))
  } else {
    stop("y should be either numeric (regression) or factor (classification)")
  }

  fit
}
