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
rpf_bridge <- function(
    processed, max_interaction = 1, ntrees = 50, splits = 30,
    split_try = 10, t_try = 0.4, deterministic = FALSE,
    parallel = FALSE, purify = FALSE, cv = FALSE,
    loss = "L2", delta = 0, epsilon = 0.1
  ) {

  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  predictors <- preprocess_predictors_fit(processed)
  outcomes <- processed$outcomes[[1]]
  
  # Check arguments
  checkmate::assert_integerish(max_interaction, lower = 1, len = 1)
  checkmate::assert_integerish(ntrees, lower = 1, len = 1)
  checkmate::assert_integerish(splits, lower = 1, len = 1)
  checkmate::assert_integerish(split_try, lower = 1, len = 1)
  
  checkmate::assert_numeric(t_try, lower = 0, upper = 1, len = 1)
  # FIXME: What is delta/epsilon and what can it look like?
  checkmate::assert_numeric(delta, lower = 0, upper = 1, len = 1)
  checkmate::assert_numeric(epsilon, lower = 0, upper = 1, len = 1)
  
  checkmate::assert_choice(
    loss, choices = c(
      "L1", "L2",
      # "median", # Discarded but present in C++ impl
      "logit", "exponential"
    ), null.ok = FALSE
  )
  
  checkmate::assert_logical(deterministic, len = 1)
  checkmate::assert_logical(parallel, len = 1)
  checkmate::assert_logical(purify, len = 1)
  checkmate::assert_logical(cv, len = 1)
  
  fit <- rpf_impl(
    Y = outcomes, X = predictors$predictors_matrix, 
    max_interaction = max_interaction, ntrees = ntrees, splits = splits,
    split_try = split_try, t_try = t_try, deterministic = deterministic,
    parallel = parallel, purify = purify, cv = cv,
    loss = loss, delta = delta, epsilon = epsilon
  )
  
  new_rpf(
    fit = fit,
    blueprint = processed$blueprint, 
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
rpf_impl <- function(
    Y, X,
    max_interaction = 1, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4,
    deterministic = FALSE, parallel = FALSE, purify = FALSE, cv = FALSE,
    loss = "L2", delta = 0, epsilon = 0.1
  ) {

  # Input validation
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)

  # Task type detection: Could be more concise
  is_binary <- length(unique(Y)) == 2
  is_integerish <- checkmate::test_integerish(Y, any.missing = FALSE)
  is_factor <- checkmate::test_factor(Y, any.missing = FALSE)
  is_numeric <- checkmate::test_numeric(Y, any.missing = FALSE)
  
  if (is_binary & is_integerish) {
    warning("y is binary integer, assuming classification task")
    Y <- as.integer(Y)
  }

  # Assume binary Y for classif
  if (is_factor | (is_binary & is_integerish)) {
    
    # Coerce to integer sequence 1 to nlevels(x)
    if (is_factor) {
      Y <- as.integer(Y)
      
      # Binary factor will result in 1,2, recode to 0,1
      if (is_binary) {
        Y <- Y - 1L
      }
    }
    
    # Recode to 0,1 in case e.g -1,1 is supplied
    if (is_binary & !identical(range(Y), c(0L, 1L))) {
      warning("y is binary integer but not 0,1 - transforming to 0,1")
      Y[Y == min(Y)] <- 0L
      Y[Y == max(Y)] <- 1L
    }
    
    # FIXME: Final assertion just in case should be superfluous
    if (is_binary) {
      checkmate::assert_integer(Y, lower = 0, upper = 1)
    } else {
      checkmate::assert_integer(Y, lower = 1)
    }
    
    fit <- new(ClassificationRPF, Y, X, loss, c(
      max_interaction, ntrees, splits, split_try, t_try,
      purify, deterministic, parallel, cv, delta, epsilon
    ))
  } else if (is_numeric) {
    #FIXME: Loss missing here?
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
