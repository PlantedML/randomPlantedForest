#' Extract predicted components from a Random Planted Forest
#'
#' @inheritParams predict.rpf
#' @param predictors [`character`] Vector of one or more column names of predictor variables in `new_data`
#' to extract components for.
#'
#' @return A numeric `matrix`
#' @export
#'
#' @examples
#'
#' # Regression task, only numeric predictors
#' train <-  mtcars[1:20, ]
#' test <-  mtcars[21:32, ]
#'
#' set.seed(23)
#' object <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
#'
#' # One component, first predictor column in test data
#' extract_components(object, test, predictors = c("cyl"))
#'
#' # Two components, 1st and 3rd predictor columns in test data
#' extract_components(object, test, predictors = c("cyl", "hp"))
#'
#' \dontrun{
#' # trying to pass all columns in test data errors since target "mpg" is included
#' extract_components(object, test, predictors = names(test))
#'
#' # passing all predictor names by dropping `"mpg"` first
#' extract_components(object, test, predictors = names(test[, -1]))
#'
#' # If `predictors` is not specified, we don't know what to do
#' extract_components(object, test)
#' }

extract_components <- function(object, new_data, predictors = NULL) {

  # Ensure selected predictors are subset of original predictors
  checkmate::assert_subset(predictors, choices = names(object$blueprint$ptypes$predictors), empty.ok = FALSE)

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)

  # Encode factors to (re-)ordered integers according to information saved during model fit
  new_data <- preprocess_predictors_predict(object, processed$predictors)

  # Indices of selected predictors
  components <- match(predictors, colnames(new_data))

  # subset new_data to contain only selected components
  new_data <- new_data[, components, drop = FALSE]

  object$fit$predict_matrix(new_data, components)
}
