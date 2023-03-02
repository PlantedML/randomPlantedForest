#' Extract predicted components from a Random Planted Forest
#'
#' Prediction components are a functional decomposition of the model prediction.
#' The sum of all components equals the overall predicted value for an observation.
#'
#' Extracts all possible components up to `max_interaction` degrees,
#' up to the value set when calling [`rpf()`]. The intercept is always included.
#' Optionally `predictors` can be specified to only include components including the given variables.
#' If `max_interaction` is greater than `length(predictors)`, the `max_interaction` will be lowered accordingly.
#'
#' @note
#' Depending on the number of predictors and `max_interaction`, the number of components will
#' increase drastically to `sum(choose(ncol(new_data), seq_len(max_interaction)))`.
#'
#' @inheritParams predict.rpf
#' @param predictors [`character`] or `NULL`: Vector of one or more column names of predictor variables
#'   in `new_data` to extract components for.
#'   If `NULL`, all variables and their interactions are returned.
#' @param max_interaction [`integer`] or `NULL`: Maximum degree of interactions to consider.
#'   Default will use the `max_interaction` parameter from the [`rpf`] object.
#'   Must be between `1` (main effects only) and the `max_interaction` of the [`rpf`] object.
#'
#' @return A `list` with elements:
#' - `m` ([`data.table`][data.table::data.table]): Components for each main effect and
#' interaction term, representing the functional decomposition of the prediction.
#' All components together with the intercept sum up
#' to the prediction.
#' For multiclass classification, the number of output columns is multiplied by
#' the number of levels in the outcome.
#' - `intercept` (`numeric(1)`): Expected value of the prediction.
#' - `x` ([`data.table`][data.table::data.table]): Copy of `new_data` containing predictors selected
#' by `predictors`.
#' - `target_levels` (`character`): For multiclass classification only: Vector of target levels
#' which can be used to disassemble `m`, as names include both term and target level.
#'
#' @export
#' @importFrom hardhat forge
#' @importFrom data.table as.data.table
#' @importFrom utils combn
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
#' (components <- predict_components(rpfit, test))
#'
#' # sums to prediction
#' cbind(
#'   m_sum = rowSums(components$m) + components$intercept,
#'   prediction = predict(rpfit, test)
#' )
#'
#' # Only get components with interactions of a lower degree, ignoring 3-way interactions
#' predict_components(rpfit, test, max_interaction = 2)
#'
#' # Only retrieve main effects
#' (main_effects <- predict_components(rpfit, test, max_interaction = 1))
#'
#' # The difference is the combined contribution of interaction effects
#' cbind(
#'   m_sum = rowSums(main_effects$m) + main_effects$intercept,
#'   prediction = predict(rpfit, test)
#' )
#'
predict_components <- function(object, new_data, max_interaction = NULL, predictors = NULL) {
  checkmate::assert_class(object, classes = "rpf")

  # if max_interaction is not provided, we use the max from the rpf model fit
  if (is.null(max_interaction)) {
    max_interaction <- object$params$max_interaction
  } else {
    # Otherwise it has to be between 1 and the max from the model fit
    checkmate::assert_int(max_interaction, lower = 1, upper = object$params$max_interaction)
  }

  # Either use all predictors in the model or a subset of interest (for speed / compactness)
  if (is.null(predictors)) {
    predictors <- names(object$blueprint$ptypes$predictors)
  } else {
    checkmate::assert_subset(predictors, choices = names(object$blueprint$ptypes$predictors))
  }

  # Check if forest is purified, if not we do that now
  if (!is_purified(object)) purify(object)

  # If max_interaction is greater than number of predictors requested we need to adjust that
  max_interaction <- min(max_interaction, length(predictors))

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)
  # Encode factors to (re-)ordered integers according to information saved during model fit
  # Need matrix version for model predictions, and keep data.table version for return (for plotting)
  new_data_matrix <- preprocess_predictors_predict(object, processed$predictors)
  data.table::setDT(new_data)

  # iterate over 1 through max_interaction, get all subsets of predictors,
  # extract the component for each combination and append them column wise
  all_components <- lapply(seq_len(max_interaction), function(i) {
    combinations <- utils::combn(predictors, i, simplify = FALSE)
    components <- lapply(combinations, function(x) {
      .predict_single_component(object, new_data_matrix, x)
      })
    do.call(cbind, args = components)
  })

  all_components <- do.call(cbind, args = all_components)
  # Get intercept (scalar), using only one row of x as input as we don't need it repeated
  intercept <- .predict_single_component(
    object, new_data_matrix[, 1, drop = FALSE], predictors = NULL
  )

  ret <- list(
    m = data.table::as.data.table(all_components),
    intercept = intercept[[1]],
    x = new_data[, predictors, with = FALSE]
  )

  # Get outcome levels for multiclass handling
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  if (length(outcome_levels) > 2) {
    ret$target_levels <- outcome_levels
  }

  class(ret) <- c("glex", "rpf_components", class(ret))
  ret

}

#' Internal function to extract a single component
#'
#' Extracts one component at a time, consisting of supplied `predictors`.
#' Not fit for general use since the preprocessing of `new_data` is only done in `predict_components`
#' to save on computing time for cases where many interactions/components are extracted.
#'
#' @noRd
#' @keywords internal
#' @return A n x 1 `matrix()`
#' @examples
#' # Regression task, only some predictors
#' train <-  mtcars[1:20, 1:4]
#' test <-  mtcars[21:32, 1:4]
#'
#' set.seed(23)
#' rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
#'
#' # Internal data preprocessing only done in predict_components to save time
#' processed <- hardhat::forge(test, rpfit$blueprint)
#' test <- randomPlantedForest:::preprocess_predictors_predict(rpfit, processed$predictors)
#'
#' # Component for interaction term of two predictors
#' randomPlantedForest:::.predict_single_component(rpfit, test, predictors = c("cyl", "hp"))
#'
#' # Retrieving the intercept
#' randomPlantedForest:::.predict_single_component(rpfit, test, predictors = NULL)
.predict_single_component <- function(object, new_data, predictors = NULL) {

  # Indices of selected predictors
  if (is.null(predictors)) {
    # Intercept is retrieved by passing -1 to predict_matrix() as a special case
    components <- -1
    out_names <- "intercept"
  } else {
    # Indices of predictors need to be in ascending order,
    # such that new_data has columns in the correct order.
    # Otherwise, predict_matrix will return wrong results
    components <- sort.int(match(predictors, colnames(new_data)))
    # subset new_data to contain only selected components
    new_data <- new_data[, components, drop = FALSE]
    # Make names, e.g. intercept, x1, x1:x2 etc. and sort predictors for consistent alphabetic order
    out_names <- paste0(sort(predictors), collapse = ":")
  }

  # Get predicted components from C++
  # new_data must be a matrix and may only contain columns required for components, in the original order
  # components must be integer indices and refer to the column position in the original training data
  ret <- object$fit$predict_matrix(new_data, components)

  # Get outcome levels for multiclass handling
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  if (length(outcome_levels) > 2) {
    # Multiclass needs disambiguation with one column for each predicted class
    colnames(ret) <- paste0(out_names, "__class:", outcome_levels)
  } else {
    # Regression and binary classif: Single-column component
    colnames(ret) <- out_names
  }

  ret
}


