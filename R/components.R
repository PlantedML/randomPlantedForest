#' Extract predicted components from a Random Planted Forest
#'
#' Prediction components are a functional decomposition of the model prediction.
#' The sum of all components equals the overall predicted value for an observation.
#'
#' * `extract_components` extracts all possible components up to `max_interaction` degrees,
#' which is defined when calling [`rpf()`].
#' * `extract_component` allows extracting only a single component for a given predictor
#' or combination of predictors.
#'
#' @inheritParams predict.rpf
#' @param predictors [`character`] Vector of one or more column names of predictor variables
#' in `new_data` to extract components for.
#'
#' @return A [`tibble`][tibble::tibble] with the same number of rows as `new_data` and one
#' column for each main or interaction term.
#' @export
#'
#' @examples
#'
#' # Regression task, only some predictors
#' train <-  mtcars[1:20, 1:4]
#' test <-  mtcars[21:32, 1:4]
#'
#' set.seed(23)
#' rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
#'
#' # Extract all components, including main effects and interaction terms up to `max_interaction`
#' (components <- extract_components(rpfit, test))
#'
#' # sums to prediction
#' cbind(rowSums(components), predict(rpfit, test))
#'
#' # Component for interaction term of two predictors
#' extract_component(rpfit, test, predictors = c("cyl", "hp"))
extract_components <- function(object, new_data) {
  # Get predictor names to keep track of them
  pred_names <- names(object$blueprint$ptypes$predictors)

  # iterate over 1 through max_interaction, get all subsets of predictors,
  # extract the component for each combination and append them column wise
  all_components <- lapply(seq_len(object$params$max_interaction), function(i) {
    combn(pred_names, i, simplify = FALSE) |>
      lapply(function(x) extract_component(object, test, x)) |>
      do.call(cbind, args = _)
  }) |>
    do.call(cbind, args = _) |>
    tibble::as_tibble()

  # Remove components with constant 0s
  all_components[, vapply(all_components, function(x) !all(x == 0), FUN.VALUE = logical(1))]
}

#' @rdname extract_components
#' @export
extract_component <- function(object, new_data, predictors = NULL) {

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

  ret <- object$fit$predict_matrix(new_data, components)
  colnames(ret) <- paste0(predictors, collapse = ":")
  # colnames(ret) <- c(".pred_m")

  tibble::as_tibble(ret)
}


