#' Extract predicted components from a Random Planted Forest
#'
#' Prediction components are a functional decomposition of the model prediction.
#' The sum of all components equals the overall predicted value for an observation.
#'
#' * `extract_components` extracts all possible components up to `max_interaction` degrees,
#'  up to the value set when calling [`rpf()`]. The intercept is always included.
#'  Optionally `predictors` can be specified to only include components including the given variables.
#' * `extract_component` allows extracting only a single component for a given predictor
#' or interaction of predictors, including the intercept if `predictors` is `NULL`.
#'
#' @note
#' Depending on the number of predictors and `max_interaction`, the number of components will
#' increase drastically to `sum(choose(ncol(new_data), seq_len(max_interaction)))`.
#'
#' @inheritParams predict.rpf
#' @param predictors [`character`] or `NULL`: Vector of one or more column names of predictor variables
#'   in `new_data` to extract components for.
#'   In `extract_components`: If `NULL`, all variables and their interactions are returned.
#'   in `extract_component`: If `NULL`, the intercept component is returned.
#' @param max_interaction [`integer`] or `NULL`: Maximum degree of interactions to consider.
#'   Default will use the `max_interaction` parameter from the [`rpf`] object.
#'   Must be between `1` (main effects only) and the `max_interaction` of the [`rpf`] object.
#'
#' @return A [`tibble`][tibble::tibble] with the same number of rows as `new_data` and one
#' column for each main and interaction term requested. For multiclass classification, the
#' number of output columns is multiplied by the number of levels in the outcome.
#' @export
#' @importFrom tibble as_tibble
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
#' (components <- extract_components(rpfit, test))
#'
#' # sums to prediction
#' cbind(rowSums(components), predict(rpfit, test))
#'
#' # Only get components with interactions of a lower degree, ignoring 3-way interactions
#' extract_components(rpfit, test, max_interaction = 2)
#'
#' # Only retrieve main effects
#' (main_effects <- extract_components(rpfit, test, max_interaction = 1))
#'
#' # The difference is the combined contribution of interaction effects
#' cbind(rowSums(main_effects), predict(rpfit, test))
#'

extract_components <- function(object, new_data, max_interaction = NULL, predictors = NULL) {

  if (is.null(max_interaction)) {
    max_interaction <- object$params$max_interaction
  } else {
    checkmate::assert_int(max_interaction, lower = 1, upper = object$params$max_interaction)
  }

  checkmate::assert_subset(
    predictors,
    choices = names(object$blueprint$ptypes$predictors),
    empty.ok = TRUE
  )

  if (is.null(predictors)) {
    predictors <- names(object$blueprint$ptypes$predictors)
  }

  # iterate over 1 through max_interaction, get all subsets of predictors,
  # extract the component for each combination and append them column wise
  all_components <- lapply(seq_len(max_interaction), function(i) {
    combinations <- utils::combn(predictors, i, simplify = FALSE)
    components <- lapply(combinations, function(x) extract_component(object, new_data, x))
    do.call(cbind, args = components)
  })

  all_components <- tibble::as_tibble(do.call(cbind, args = all_components))
  intercept <- extract_component(object, new_data, predictors = NULL)
  all_components <- cbind(intercept, all_components)

  data.table::as.data.table(all_components)
}

#' @rdname extract_components
#' @export
#' @examples
#' # Regression task, only some predictors
#' train <-  mtcars[1:20, 1:4]
#' test <-  mtcars[21:32, 1:4]
#'
#' set.seed(23)
#' rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
#'
#' # Component for interaction term of two predictors
#' extract_component(rpfit, test, predictors = c("cyl", "hp"))
#'
#' # Retrieving the intercept
#' extract_component(rpfit, test, predictors = NULL)
extract_component <- function(object, new_data, predictors = NULL) {
  # Ensure selected predictors are subset of original predictors
  # but allow NULL if intercept is requested
  checkmate::assert_subset(
    predictors,
    choices = names(object$blueprint$ptypes$predictors),
    empty.ok = TRUE
  )

  # Check if forest is purified, if not we do that now
  # First, check that purification is reported correctly.
  checkmate::assert_logical(object$fit$is_purified(), any.missing = FALSE, len = 1)
  if (!object$fit$is_purified()) {
    purify(object)
  }

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)
  # Get outcome levels for later multiclass handling
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  # Encode factors to (re-)ordered integers according to information saved during model fit
  new_data <- preprocess_predictors_predict(object, processed$predictors)

  # Indices of selected predictors
  if (is.null(predictors)) {
    # Intercept is retrieved by passing -1 to predict_matrix() as a special case
    components <- -1
  } else {
    # Indices of predictors need to be in ascending order such that new_data has columns in the correct order
    # Otherwise, predict_matrix will return wrong results
    components <- sort(match(predictors, colnames(new_data)))
    # subset new_data to contain only selected components
    new_data <- new_data[, components, drop = FALSE]
  }

  # Get predicted components from C++
  # new_data may only contain columns required for components
  ret <- object$fit$predict_matrix(new_data, components)

  # Make names, e.g. intercept, x1, x1:x2 etc. and sort predictors for consistent alphabetic order
  out_names <- ifelse(is.null(predictors), "intercept", paste0(sort(predictors), collapse = ":"))

  if (length(outcome_levels) > 2) {
    # Multiclass needs disambiguation with one column for each predicted class
    colnames(ret) <- paste0(out_names, "_", outcome_levels)
  } else {
    # Regression and binary classif: Single-column component
    colnames(ret) <- out_names
  }

  ret
}


