#' Random Planted Forest Predictions
#'
#' @param object A fit object as returned by [`rpf()`].
#' @param new_data Input for new observations.
#' @param type `["numeric"]`: One of `"numeric"` for regression outcomes,
#' `"class"` for class predictions or `"prob"` for probability predictions.
#' @param components `[0]` TODO.
#' @param ... Unused.
#'
#' @return For regression: A data.frame with column `.pred` with the same
#' number of rows as `new_data`.
#' For Classification: TBI
#'
#' @export
#' @importFrom hardhat forge validate_prediction_size
#'
#' @examples
#' \dontrun{
#' rpfit <- rpf(y = mtcars$mpg, x = as.matrix(mtcars[, c("cyl", "wt")]))
#' predict(rpfit, mtcars[, c("cyl", "wt")], components = 0)
#' }
predict.rpf <- function(object, new_data, type = "numeric",
                        components = 0, ...) {

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)

  out <- predict_rpf_bridge(type, object, processed$predictors, components, ...)

  hardhat::validate_prediction_size(out, new_data)

  out
}


# Bridge: Passes new data to corresponding predict function
predict_rpf_bridge <- function(type, object, predictors, ...) {

  type <- match.arg(type, choices = c("numeric", "class", "prob"))
  predictors <- as.data.table(predictors)
  
  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }
  
  # Re-order factor levels according to saved order 
  factor_cols <- names(object$factor_levels)
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := Map(factor, .SD, object$factor_levels, ordered = TRUE), .SDcols = factor_cols]
  }
  
  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }
  predictors_matrix <- as.matrix(predictors)

  switch(
    type,
    numeric = predict_rpf_numeric(object, predictors_matrix, ...),
    class = predict_rpf_class(object, predictors_matrix, ...),
    prob = predict_rpf_prob(object, predictors_matrix, ...)
  )
}

# Predict function for numeric outcome / regression
predict_rpf_numeric <- function(object, new_data, components, ...){
  pred <- object$fit$predict_matrix(new_data, components)
  out <- hardhat::spruce_numeric(pred)

  out
}


# Classification ----------------------------------------------------------

# Predict function for classification: Class prediction
predict_rpf_class <- function(object, new_data, components, ...){
  if (is.integer(object$blueprint$ptypes$outcomes[[1]])) {
    outcome_levels <- c("0", "1")
  } else {
    outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  }
  
  # Predict probability
  pred_prob <- predict_rpf_prob(object, new_data, components, ...)
  
  # For each instance, class with higher probability
  pred_class <- factor(outcome_levels[max.col(pred_prob)], levels = outcome_levels)
  out <- hardhat::spruce_class(pred_class)

  out
}

# Predict function for classification: Probability prediction
predict_rpf_prob <- function(object, new_data, components, ...){
  if (is.integer(object$blueprint$ptypes$outcomes[[1]])) {
    outcome_levels <- c("0", "1")
  } else {
    outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  }
  
  pred_raw <- object$fit$predict_matrix(new_data, components)
  
  if (object$loss %in% c("logit", "exponential")) {
    # logit^-1 transformation for logit/exp loss
    pred_prob <- 1/(1 + exp(-pred_raw))
  } else if (object$loss %in% c("L1", "L2")) {
    # Truncate probabilities at [0,1] for L1/L2 loss
    pred_prob <- pmax(0, pmin(1, pred_raw))
  }
  
  #FIXME: Multiclass
  pred <- cbind(1 - pred_prob, pred_prob)
  out <- hardhat::spruce_prob(outcome_levels, pred)

  out
}
