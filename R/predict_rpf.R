#' Random Planted Forest Predictions
#'
#' @param object A fit object of class [`rpf`].
#' @param new_data Data for new observations to predict.
#' @param type `"numeric"` for regression outcomes,
#' `"class"` for class predictions or `"prob"` for probability predictions.
#'
#' For classification and `loss = "L1"` or `"L2"`, `"numeric"` yields raw
#' predictions which are not guaranteed to be valid probabilities in `[0, 1]`.
#' For `type = "prob"`, these are truncated to ensure this property.
#'
#' If `loss` is `"logit"` or `"exponential"`, `type = "link"` is an alias
#' for `type = "numeric"`, as in this case the raw predictions have the
#' additional interpretation similar to the linear predictor in a [`glm`].
#' @param components `[0]` TODO.
#' @param ... Unused.
#'
#' @return For regression: A [`tbl`][tibble::tibble] with column `.pred` with
#' the same number of rows as `new_data`.
#'
#' For classification: A [`tbl`][tibble::tibble] with one column for each
#' level in `y` containing class probabilities if `type = "prob"`.
#' For `type = "class"`, one column `.pred` with class predictions is returned.
#' For `type = "numeric"` or `"link"`, one column `.pred` with raw predictions.
#'
#' @export
#' @importFrom hardhat forge validate_prediction_size
#'
#' @examples
#' # Regression with L2 loss
#' rpfit <- rpf(y = mtcars$mpg, x = mtcars[, c("cyl", "wt")])
#' predict(rpfit, mtcars[, c("cyl", "wt")], components = 0)
predict.rpf <- function(object, new_data,
                        type = ifelse(object$mode == "regression", "numeric", "prob"),
                        components = 0, ...) {

  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)

  out <- predict_rpf_bridge(type, object, processed$predictors, components, ...)

  hardhat::validate_prediction_size(out, new_data)

  out
}


# Bridge: Passes new data to corresponding predict function
predict_rpf_bridge <- function(type, object, predictors, ...) {
  type <- match.arg(type, choices = c("numeric", "class", "prob", "link"))
  predictors <- preprocess_predictors_predict(object, predictors)

  if (object$mode == "regression") {
    if (type != "numeric") {
      warning(
        paste0(
          "Only predict type 'numeric' supported for regression, ",
          "but type is set to '", type, "'. Setting type to 'numeric'."
        )
      )
    }
    type <- "numeric"
  } else if (object$mode == "classification") {
    if (object$loss %in% c("logit", "logit_2", "exponential", "exponential_2")) {
      # numeric yields raw predictions, which is the link in the binary case
      if (type == "link") type <- "numeric"
    }
  }

  switch(type,
    numeric = predict_rpf_numeric(object, predictors, ...),
    class = predict_rpf_class(object, predictors, ...),
    prob = predict_rpf_prob(object, predictors, ...)
  )
}

# Predict function for numeric outcome / regression
predict_rpf_numeric <- function(object, new_data, components, ...) {
  pred <- object$fit$predict_matrix(new_data, components)

  if (ncol(pred) == 1) {
    # Regression or binary case: just return predictions as-s, single column
    # Convert n x 1 matrix to numeric vector
    pred <- as.numeric(pred)
    out <- hardhat::spruce_numeric(pred)
  } else {
    # multiclass case for 'link' type: get prediction matrix, clean up
    # via spruce_prob, but no transformation
    outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
    out <- hardhat::spruce_prob(outcome_levels, pred)
  }

  out
}


# Classification ----------------------------------------------------------

# Predict function for classification: Class prediction
predict_rpf_class <- function(object, new_data, components, ...) {
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  # Predict probability
  pred_prob <- predict_rpf_prob(object, new_data, components, ...)

  # For each instance, class with higher probability
  pred_class <- factor(outcome_levels[max.col(pred_prob)], levels = outcome_levels)
  out <- hardhat::spruce_class(pred_class)

  out
}

# Predict function for classification: Probability prediction
predict_rpf_prob <- function(object, new_data, components, ...) {
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  pred_raw <- object$fit$predict_matrix(new_data, components)
  if (object$loss %in% c("logit", "logit_2", "exponential", "exponential_2")) {
    if (ncol(pred_raw) == 1) {
      # logit^-1 transformation for logit/exp loss
      pred_prob <- 1 / (1 + exp(-pred_raw))
    } else {
      # FIXME:
      # softmax() defined in utils.R, should be identical to logit^-1 for 
      # binary case but not properly tested yet
      pred_prob <- softmax(pred_raw)
    }
  } else if (object$loss %in% c("L1", "L2")) {
    # Truncate probabilities at [0,1] for L1/L2 loss
    # pred_prob <- pmax(0, pmin(1, pred_raw)) # did not work for matrix predictions
    pred_prob <- apply(pred_raw, 2, function(col) pmax(0, pmin(1, col)))
  }

  # FIXME: Hacky solution
  if (ncol(pred_prob) == 1) { # Binary classif yields n x 1 prediction matrix
    pred <- cbind(1 - pred_prob, pred_prob)
    out <- hardhat::spruce_prob(outcome_levels, pred)
  } else { # otherwise one column per level
    out <- hardhat::spruce_prob(outcome_levels, pred_prob)
  }

  out
}
