#' The components of a Random Planted Forest's prediction
#'
#' Evaluates the components of a model at the rows of `new_data`: its intercept and one column per main effect and
#' interaction. They are a functional decomposition of the prediction: for every row, the intercept plus the
#' components equals the prediction. For a pruned model, the default of [rpf()], these are its kept components.
#'
#' Extracts all possible components up to `max_interaction` degrees,
#' up to the value set when calling [`rpf()`]. The intercept is always included.
#' Optionally `predictors` can be specified to only include components made of the given predictors;
#' the sum of the other components is returned as `remainder`.
#' If `max_interaction` is greater than `length(predictors)`, the `max_interaction` will be lowered accordingly.
#'
#' All components are evaluated in one pass over the rows of `new_data`.
#' Components that the purified forest does not hold, such as sets of predictors
#' that never appear together in a tree, are zero columns.
#'
#' The components are only unique for a purified forest; [purify()] says how they are identified.
#' An unpurified forest is purified in place on first use with the defaults of [purify()].
#' For other settings, such as more threads, call [purify()] before `predict_components()`.
#'
#' @note
#' Depending on the number of predictors and `max_interaction`, the number of components will
#' increase drastically to `sum(choose(ncol(new_data), seq_len(max_interaction)))`.
#'
#' @section Pruned forests:
#' For an `rpf_pruned` object from [prune_components()] or [rpf()], `m` holds only the kept components, each scaled
#' by its coefficient `beta_S`, and `intercept` is the pruned model's intercept `alpha`, so that they sum to
#' [`predict()`][predict.rpf_pruned] of the pruned model. Components that are not kept are absent, not zero columns.
#' If `predictors` or `max_interaction` leave out kept components, their sum is returned as `remainder`.
#'
#' @inheritParams predict.rpf
#' @param object `[rpf | rpf_pruned]`: A fitted [`rpf`] model, or a pruned one from [prune_components()] or [rpf()].
#' @param predictors `[character | NULL: NULL]`: Vector of one or more names of the forest's predictors
#'   to extract components for. These are the columns the forest was fitted on: with a formula, its terms,
#'   such as `log(x1)`.
#'   If `NULL`, all predictors and their interactions are returned, in the order of the fit's columns.
#' @param max_interaction `[integer(1) | NULL: NULL]`: Maximum degree of
#'   interactions to consider. `NULL` uses the `max_interaction` parameter from the [`rpf`] object.
#'   Must be between `1` (main effects only) and the `max_interaction` of the [`rpf`] object.
#'
#' @return A `list` with elements:
#' - `m` ([`data.table`][data.table::data.table]): Components for each main effect and
#' interaction term, representing the functional decomposition of the prediction.
#' All components together with the intercept sum up
#' to the prediction.
#' For multiclass classification, the number of output columns is multiplied by
#' the number of `target_levels`, the levels of the outcome: one column per component and level, named
#' `<component>__class:<level>`. With `loss = "logit"` or `"exponential"` a component's columns sum to 0 across the
#' levels, as the scores of [`predict()`][predict.rpf] with `type = "link"` do.
#' - `intercept` (`numeric`): Expected value of the prediction.
#' For multiclass classification, one value per target level, named by `target_levels`.
#' - `x` ([`data.table`][data.table::data.table]): The predictors selected by `predictors`, as the forest
#' sees them in `new_data` (with a formula, its terms, such as `log(x1)`).
#' - `target_levels` (`character`): For multiclass classification only: Vector of target levels
#' which can be used to disassemble `m`, as names include both term and target level.
#' - `remainder` (`numeric`): Only if `m` leaves out components of the forest: if `predictors` leaves out
#' some of the forest's predictors, if `max_interaction` (lowered to the number of `predictors`) is below that
#' of the [`rpf`] object, or if the forest was purified with a `maxp_interaction` below it (see [purify()]). The prediction minus `intercept` and `m`'s row sums,
#' so that the three sum to the prediction.
#' For multiclass classification, one column per target level.
#'
#' @export
#' @importFrom hardhat forge
#' @importFrom data.table as.data.table
#' @importFrom utils combn
#'
#' @seealso [purify()] for how the components are identified, [prune_components()] for the pruned model's
#'   components. The glex package plots the result with `autoplot()`.
#'
#' @examples
#' # The default model: its kept components add up to its prediction
#' set.seed(1)
#' model <- rpf(mpg ~ wt + hp + cyl, data = mtcars)
#' components <- predict_components(model, mtcars[1:3, ])
#' components$m
#' components$intercept + rowSums(components$m)
#' predict(model, mtcars[1:3, ])
#'
#' # A forest: all its components up to its max_interaction
#' train <-  mtcars[1:20, 1:4]
#' test <-  mtcars[21:32, 1:4]
#'
#' set.seed(23)
#' rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, n_families = 30, prune = FALSE)
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
predict_components <- function(object, new_data, ...) {
  UseMethod("predict_components")
}

#' @rdname predict_components
#' @export
predict_components.rpf <- function(object, new_data, max_interaction = NULL, predictors = NULL, ...) {
  rlang::check_dots_empty()
  check_rpf_alive(object)

  # if max_interaction is not provided, we use the max from the rpf model fit
  if (is.null(max_interaction)) {
    max_interaction <- object$params$max_interaction
  } else {
    # Otherwise it has to be between 1 and the max from the model fit
    checkmate::assert_int(max_interaction, lower = 1, upper = object$params$max_interaction)
  }

  # Either use all predictors in the model or a subset of interest (for speed / compactness).
  # The fit's own columns, in the fit's order: a formula's terms such as log(x1) are columns of their own,
  # and hardhat::mold() puts factors last.
  fit_predictors <- fit_predictor_names(object)
  if (is.null(predictors)) {
    predictors <- fit_predictors
  } else {
    checkmate::assert_subset(predictors, choices = fit_predictors)
  }

  if (!is_purified(object)) {
    purify(object)
  }

  # If max_interaction is greater than number of predictors requested we need to adjust that
  max_interaction <- min(max_interaction, length(predictors))

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)
  check_no_missing(processed$predictors, "{.arg new_data}")
  # Encode factors to (re-)ordered integers according to information saved during model fit
  # Need matrix version for model predictions, and keep data.table version for return (for plotting)
  new_data_matrix <- preprocess_predictors_predict(object, processed$predictors)
  # a copy: setDT() would fail on a matrix and turn the caller's data.frame into a data.table.
  # The forged predictors, so that a formula's terms such as log(x1) are columns too
  new_data_dt <- data.table::as.data.table(processed$predictors)

  # all subsets of predictors up to max_interaction, extracted in one pass (the components the purified forest does
  # not store are zero columns)
  combinations <- unlist(
    lapply(seq_len(max_interaction), function(i) utils::combn(predictors, i, simplify = FALSE)),
    recursive = FALSE
  )
  all_components <- .predict_components(object, new_data_matrix, combinations)
  # Get intercept (scalar), using only one row of x as input as we don't need it repeated
  # (the intercept does not depend on the row; with zero rows, one row of zeros, as C++ needs a row)
  intercept_data <- new_data_matrix[, 1, drop = FALSE]
  if (nrow(intercept_data) == 0) {
    intercept_data <- matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, colnames(new_data_matrix)[1]))
  }
  intercept <- .predict_single_component(
    object,
    intercept_data,
    predictors = NULL
  )

  ret <- list(
    m = data.table::as.data.table(all_components),
    intercept = intercept[[1]],
    x = new_data_dt[, predictors, with = FALSE]
  )

  # Get outcome levels for multiclass handling (one column per level: multiclass logit's scores centred as
  # predict(type = "link") reports them)
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  if (length(outcome_levels) > 2) {
    ret$intercept <- stats::setNames(intercept[1, ], outcome_levels)
    ret$target_levels <- outcome_levels
  }

  # If m leaves out components of the forest (max_interaction here is smaller than that of the model fit,
  # or predictors leaves out some of the fit's predictors), we calculate a remainder term
  # to make m's sum up to SHAPs etc, see https://github.com/PlantedML/glex/issues/11
  # The same after purify(maxp_interaction = k) below the fit's order: the components above k
  # are zero columns, and the remainder holds their sum.

  left_out <- max_interaction < object$params$max_interaction || !all(fit_predictors %in% predictors)
  if ((left_out || object$fit$is_purify_capped()) && nrow(new_data_matrix) == 0) {
    # zero rows: a remainder with zero rows, as predict() cannot run on them
    ret$remainder <- if (length(outcome_levels) > 2) {
      data.table::as.data.table(stats::setNames(rep(list(numeric(0)), length(outcome_levels)), outcome_levels))
    } else {
      numeric(0)
    }
  } else if (left_out || object$fit$is_purify_capped()) {
    pred <- predict.rpf(object, new_data = new_data, type = "numeric")

    # handling differs if multiclass
    if (length(outcome_levels) > 2) {
      ret$remainder <- calc_remainders_multiclass(all_components, outcome_levels, pred, intercept)
    } else {
      # regression and binary classif: straight forward.
      ret$remainder <- pred[[1]] - ret$intercept - rowSums(all_components)
    }
  }

  class(ret) <- c("glex", "rpf_components", class(ret))
  ret
}

#' Internal function to extract several components in one pass
#'
#' The components consisting of each of the predictor sets `sets` (a list of character vectors), at the rows of the
#' preprocessed matrix `new_data` (all predictors, as from `preprocess_predictors_predict()`). Only the components the
#' purified forest stores are evaluated; the others are zero. Each column equals `.predict_single_component()`.
#'
#' @noRd
#' @keywords internal
#' @return A n x `length(sets)` `matrix()` (times the number of outcome levels for multiclass), named as
#'   `.predict_single_component()` names its columns.
.predict_components <- function(object, new_data, sets) {
  idx <- lapply(sets, function(s) sort.int(match(s, colnames(new_data))))
  out_names <- vapply(sets, function(s) paste0(sort(s), collapse = ":"), character(1))
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  if (length(outcome_levels) > 2) {
    out_names <- paste0(rep(out_names, each = length(outcome_levels)), "__class:", outcome_levels)
  }
  if (nrow(new_data) == 0) {
    # the C++ code needs at least one row
    ret <- matrix(numeric(0), nrow = 0, ncol = length(out_names))
  } else {
    # one column per level (multiclass logit: its K - 1 score columns per set centred, reported_scores())
    ret <- reported_scores(object, object$fit$predict_components_matrix(new_data, idx, 0L), length(sets))
  }
  colnames(ret) <- out_names
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
#' rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, n_families = 30, prune = FALSE)
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

  # Component inputs are often low-cardinality, so predict each distinct row once.
  # match(x, x) codes values exactly as the index of their first occurrence;
  # folding codes column by column keeps keys <= n^2 < 2^53.
  first_index <- function(x) match(x, x)
  n <- nrow(new_data)
  key <- first_index(new_data[, 1L])
  for (j in seq_len(ncol(new_data))[-1L]) {
    key <- first_index((key - 1) * n + first_index(new_data[, j]))
  }
  first <- which(key == seq_len(n))

  # Get predicted components from C++
  # new_data must be a matrix and may only contain columns required for components, in the original order
  # components must be integer indices and refer to the column position in the original training data
  ret <- object$fit$predict_matrix(new_data[first, , drop = FALSE], components, 0L)
  # one column per level (multiclass logit: its K - 1 score columns centred, reported_scores())
  ret <- reported_scores(object, ret[match(key, first), , drop = FALSE])

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
