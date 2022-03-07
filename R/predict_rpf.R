#' Random Planted Forest Predictions
#'
#' @param object A fit object as returned by [`rpf()`].
#' @param new_data Input for new observations.
#' @param components [`0`]: TODO
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

  predictors <- as.matrix(predictors)

  switch(
    type,
    numeric = predict_rpf_numeric(object, predictors, ...),
    class = predict_rpf_class(object, predictors, ...),
    prob = predict_rpf_prob(object, predictors, ...)
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
  pred <- object$fit$predict_matrix(new_data, components)
  out <- hardhat::spruce_class(pred)

  out
}

# Predict function for classification: Probability prediction
predict_rpf_prob <- function(object, new_data, components, ...){
  pred <- object$fit$predict_matrix(new_data, components)
  out <- hardhat::spruce_prob(pred)

  out
}
