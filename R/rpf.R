#' Random Planted Forest
#'
# FIXME: Parameters need describing
#' @param x Feature matrix or `data.frame`.
#' @param y Target vector for use with `x`.
#' @param formula Formula specification, e.g. y ~ x1 + x2.
#' @param data A `data.frame` for use with `formula`.
#' @param max_interaction Maximum level of interaction determining maximum
#'   number of split dimensions for a tree
#' @param ntrees Number of trees generated per family
#' @param splits Number of performed splits for each tree family
#' @param split_try TODO
#' @param t_try TODO
#' @param deterministic Choose whether approach deterministic or random
#' @param parallel Perform algorithm in parallel or serialized
#' @param purify Whether the forest should be purified
#' @param cv Determines if cross validation is performed
#' @param loss "L2" or "logit" / TODO
#' @param delta TODO
#' @param epsilon TODO
#' @param ... (Ignored)
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
#' \dontrun{
#' # Regression with x and y
#' rpfit <- rpf(x = mtcars[, c("cyl", "wt")], y = mtcars$mpg)
#'
#' # Regression with formula
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
#' }
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
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# XY method - matrix
#' @export
#' @rdname rpf
rpf.matrix <- function(x, y, ...) {
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE, indicators = "none")
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
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  rpf_bridge(processed, ...)
}

# Bridge: Calls rpf_impl() with processed input
#' @importFrom hardhat validate_outcomes_are_univariate
#' @importFrom data.table .SD ':=' as.data.table
rpf_bridge <- function(processed, ...) {

  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  outcomes <- processed$outcomes[[1]]
  predictors <- as.data.table(processed$predictors)
  
  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }
  
  # Factor predictors: Order by response (see https://doi.org/10.7717/peerj.6339)
  factor_cols <- names(which(sapply(predictors, is.factor)))
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, order_factor_by_response, y = outcomes), .SDcols = factor_cols]
  }
  
  # Save re-ordered factor levels
  factor_levels <- hardhat::get_levels(predictors)
  
  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }
  predictors_matrix <- as.matrix(predictors)
  
  fit <- rpf_impl(Y = outcomes, X = predictors_matrix, ...)

  new_rpf(
    fit = fit,
    blueprint = processed$blueprint, 
    factor_levels = factor_levels
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
rpf_impl <- function(Y, X, max_interaction = 1, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4,
                deterministic = FALSE, parallel = FALSE, purify = FALSE, cv = FALSE,
                loss = "L2", delta = 0, epsilon = 0.1) {

  # FIXME: Proper classif/regression detection and switching

  # Assume binary Y for classif
  if (length(unique(Y)) == 2) {
    fit <- new(ClassificationRPF, Y, X, loss, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv, delta, epsilon
    ))
  } else {
    fit <- new(RandomPlantedForest, Y, X, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv
    ))
  }

  fit
}
