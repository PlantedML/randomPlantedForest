#' Random Planted Forest
#'
# FIXME: Parameters need describing
# @param Y Vector of target features
#' @param x Feature matrix or `formula`.
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
#' @param (Ignored)
#'
#' @return C++ object of class `"Rcpp_RandomPlantedForest"`.
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
rpf <- function(x, max_interaction = 1, ntrees = 50, splits = 30,
                split_try = 10, t_try = 0.4, deterministic = FALSE,
                parallel = FALSE, purify = FALSE, cv = FALSE,
                loss = "L2", delta = 0, epsilon = 0.1, ...) {
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
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE)
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
rpf_bridge <- function(processed, ...) {

  hardhat::validate_outcomes_are_univariate(processed$outcomes)

  # FIXME: Categorical features need consideration
  # Can't coerce to numeric matrix if categoricals aren't encoded
  predictors <- as.matrix(processed$predictors)
  outcomes <- processed$outcomes[[1]]

  fit <- rpf_impl(Y = outcomes, X = predictors, ...)

  new_rpf(
    fit = fit,
    blueprint = processed$blueprint
  )
}

# Intermediate to hold model object with blueprint used for prediction
new_rpf <- function(fit, blueprint) {

  hardhat::new_model(
    fit = fit,
    blueprint = blueprint,
    class = "rpf"
  )
}

# Main fitting function and interface to C++ implementation
rpf_impl <- function(Y, X, max_interaction = 1, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4,
                deterministic = FALSE, parallel = FALSE, purify = FALSE, cv = FALSE,
                loss = "L2", delta = 0, epsilon = 0.1) {
  if (!missing(loss) | !missing(delta) | !missing(epsilon)) {
    return(new(ClassificationRPF, Y, X, loss, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv, delta, epsilon
    )))
  }
  return(new(RandomPlantedForest, Y, X, c(
    max_interaction, ntrees, splits, split_try, t_try,
    purify, deterministic, parallel, cv
  )))
}
