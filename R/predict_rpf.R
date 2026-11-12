#' Predict with a Random Planted Forest
#'
#' Predictions of a forest (class `rpf`, from [rpf()] with `prune = FALSE`) for new data: the mean of its tree
#' families' predictions. For classification, `type` chooses probabilities, classes or the raw scores. A pruned model,
#' the default of [rpf()], has its own method, [predict.rpf_pruned()], with the same arguments.
#'
#' `new_data` is prepared as the training data were: by the formula or recipe of the fit, and with factor levels coded
#' as in the fit (see "Factor predictors" in [rpf()]). Missing values are an error.
#'
#' @param object `[rpf]`: A fitted [`rpf`] model.
#' @param new_data `[data.frame | matrix]`: New observations to predict.
#' @param type `[character(1)]`: `"numeric"` for regression outcomes,
#' `"class"` for class predictions or `"prob"` for probability predictions.
#' Defaults to `"numeric"` for regression and `"prob"` for classification.
#'
#' For classification and `loss = "L1"` or `"L2"`, `"numeric"` yields raw
#' predictions which are not guaranteed to be valid probabilities in `[0, 1]`.
#' For `type = "prob"`, these are truncated to ensure this property.
#'
#' For classification, `type = "link"` is an alias for `type = "numeric"`.
#' If `loss` is `"logit"` or `"exponential"`, the raw predictions are the
#' scores on the link scale, like the linear predictor of a [`glm`]: with two
#' classes one score, the log-odds of the second class; with more, one score
#' per class (`"exponential"`, summing to 0) or per class but the first
#' (`"logit"`, whose first class has score 0). The probabilities of
#' `type = "prob"` are their softmax (see "Classification losses" in [rpf()]).
#' If `loss` is `"L1"` or `"L2"`, the link is the identity.
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the `nthreads`
#'   the forest was fitted with, capped at the available cores.
#' @param ... Not currently used, but required for extensibility. Unknown
#'   arguments are an error.
#'
#' @return For regression: A [`tbl`][tibble::tibble] with column `.pred` with
#' the same number of rows as `new_data`.
#'
#' For classification: A [`tbl`][tibble::tibble] with one column `.pred_<level>`
#' for each level in `y` containing class probabilities if `type = "prob"`.
#' For `type = "class"`, one column `.pred_class` with class predictions is returned.
#' For `type = "numeric"` or `"link"`, raw predictions are returned: one column
#' `.pred` for binary outcomes (the log-odds of the second level with
#' `loss = "logit"` or `"exponential"`), and one column `.pred_<level>` per level
#' for multiclass outcomes. With `loss = "logit"` or `"exponential"`, these are
#' one score per class that sum to 0 across the classes, and the probabilities
#' are their softmax. (The multiclass logit loss fits `K - 1` scores against the
#' first level; they are reported with the first level's score 0 joined and
#' centred across the `K` classes, which leaves the probabilities as they are.)
#'
#' @seealso [predict.rpf_pruned()] for the pruned model, [predict_components()] for the prediction as a sum of
#'   components.
#' @export
#' @importFrom hardhat forge validate_prediction_size
#'
#' @examples
#' # Regression: a tibble with the column .pred
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, n_families = 20, prune = FALSE)
#' predict(rpfit, mtcars[1:3, ])
#'
#' # Classification: probabilities by default, or classes, or the raw scores
#' rpfit <- rpf(Species ~ ., data = iris, n_families = 20, prune = FALSE)
#' predict(rpfit, iris[c(1, 51, 101), ])
#' predict(rpfit, iris[c(1, 51, 101), ], type = "class")
#' predict(rpfit, iris[c(1, 51, 101), ], type = "link")
predict.rpf <- function(
  object,
  new_data,
  type = ifelse(object$mode == "regression", "numeric", "prob"),
  nthreads = NULL,
  ...
) {
  rlang::check_dots_empty()
  check_rpf_alive(object)
  if (!is.null(nthreads)) {
    checkmate::assert_int(nthreads, lower = 1)
  }
  # 0 tells C++ to use the forest's own nthreads
  nthreads <- as.integer(if (is.null(nthreads)) 0L else nthreads)
  # Enforces column order, type, column names, etc
  processed <- hardhat::forge(new_data, object$blueprint)
  check_no_missing(processed$predictors, "{.arg new_data}")

  out <- predict_rpf_bridge(type, object, processed$predictors, nthreads)

  hardhat::validate_prediction_size(out, new_data)

  out
}


# Bridge: Passes new data to corresponding predict function
predict_rpf_bridge <- function(type, object, predictors, nthreads = 0L) {
  type <- match.arg(type, choices = c("numeric", "class", "prob", "link"))
  predictors <- preprocess_predictors_predict(object, predictors)

  if (object$mode == "regression") {
    if (type != "numeric") {
      cli::cli_warn(c(
        "Only {.code type = \"numeric\"} is supported for regression.",
        "i" = "Using it instead of {.val {type}}."
      ))
    }
    type <- "numeric"
  } else if (object$mode == "classification") {
    # numeric yields the raw predictions, which are the link for every loss
    # (the identity for L1 and L2)
    if (type == "link") type <- "numeric"
  }

  switch(
    type,
    numeric = predict_rpf_numeric(object, predictors, nthreads),
    class = predict_rpf_class(object, predictors, nthreads),
    prob = predict_rpf_prob(object, predictors, nthreads)
  )
}

# The raw predictions, one row per row of new_data. The C++ code needs at least
# one row, so zero rows give a zero-row matrix with the fit's number of columns.
predict_rpf_raw <- function(object, new_data, nthreads = 0L) {
  if (nrow(new_data) == 0) {
    n_levels <- length(levels(object$blueprint$ptypes$outcomes[[1]]))
    n_cols <- if (object$mode == "regression" || n_levels == 2) 1L else length(model_outcome_levels(object))
    return(matrix(numeric(0), nrow = 0, ncol = n_cols))
  }
  object$fit$predict_matrix(new_data, 0, nthreads)
}

# Predict function for numeric outcome / regression
predict_rpf_numeric <- function(object, new_data, nthreads = 0L) {
  predict_rpf_numeric_raw(object, predict_rpf_raw(object, new_data, nthreads))
}

# The output of type "numeric" from the raw predictions `pred` (n x value_size) of the forest `object` or of a model
# built on it (its pruned models)
predict_rpf_numeric_raw <- function(object, pred) {
  if (ncol(pred) == 1) {
    # Regression or binary case: just return predictions as-is, single column
    # Convert n x 1 matrix to numeric vector
    pred <- as.numeric(pred)
    out <- hardhat::spruce_numeric(pred)
  } else {
    # multiclass case for 'link' type: get prediction matrix, clean up
    # via spruce_prob, but no transformation
    # One score per level (multiclass logit: its K - 1 scores centred with the reference's 0, reported_scores())
    out <- hardhat::spruce_prob(levels(object$blueprint$ptypes$outcomes[[1]]), reported_scores(object, pred))
  }

  out
}


# Classification ----------------------------------------------------------

# Predict function for classification: Probability prediction
predict_rpf_prob <- function(object, new_data, nthreads = 0L) {
  predict_rpf_prob_raw(object, predict_rpf_raw(object, new_data, nthreads))
}

# The probabilities of the raw predictions `pred_raw` (n x value_size) of the forest `object` or of a model built on it
# (its pruned models): the link of the forest's loss
predict_rpf_prob_raw <- function(object, pred_raw) {
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  if (length(outcome_levels) == 2) {
    # Binary outcome

    if (object$params$loss %in% c("logit", "exponential")) {
      # logit^-1 transformation for logit/exp loss
      pred_prob <- 1 / (1 + exp(-pred_raw))
    } else if (object$params$loss %in% c("L1", "L2")) {
      # Truncate probabilities at [0,1] for L1/L2 loss
      pred_prob <- pmax(0, pmin(1, pred_raw))
    }

    # Ensure a plain numeric vector in binary case
    pred_prob <- as.numeric(pred_prob)
    # Binary classif yields two columns ordered by outcome levels
    pred_prob <- cbind(1 - pred_prob, pred_prob)
  } else {
    # Multiclass

    if (object$params$loss == "logit") {
      # Raw predictions are K-1 logits relative to the reference class
      # (first outcome level), so p_ref = 1 / (1 + sum(exp(W))) and
      # p_k = exp(W_k) * p_ref. Rows sum to 1 by construction.
      exp_raw <- exp(pred_raw)
      denom <- 1 + rowSums(exp_raw)
      pred_prob <- cbind(1 / denom, exp_raw / denom)
      # A row whose scores overflow exp() (above about 709) would give Inf / Inf: the same probabilities with the
      # row's largest score m taken out, p_ref = exp(-m) / d and p_k = exp(W_k - m) / d, d = exp(-m) + sum(exp(W - m))
      big <- which(is.infinite(denom))
      if (length(big)) {
        W <- pred_raw[big, , drop = FALSE]
        m <- apply(W, 1, max)
        E <- exp(W - m)
        d <- exp(-m) + rowSums(E)
        pred_prob[big, ] <- cbind(exp(-m) / d, E / d)
      }
    } else if (object$params$loss == "exponential") {
      # softmax for multi-class
      pred_prob <- softmax(pred_raw)
    } else if (object$params$loss %in% c("L1", "L2")) {
      # Clamp to [0,1] and renormalize rows
      pred_prob <- pmin(1, pmax(0, pred_raw))
      # pmin/pmax drop dimensions; restore matrix shape explicitly
      dim(pred_prob) <- dim(pred_raw)
      rs <- rowSums(pred_prob)
      rs[!is.finite(rs) | rs <= 0] <- 1
      pred_prob <- pred_prob / rs
    }
  }

  hardhat::spruce_prob(outcome_levels, pred_prob)
}

# Class prediction: Take probability prediction and convert to class labels
predict_rpf_class <- function(object, new_data, nthreads = 0L) {
  predict_rpf_class_of_prob(object, predict_rpf_prob(object, new_data, nthreads))
}

# The classes of the raw predictions `pred_raw` (the most probable), for the pruned models
predict_rpf_class_raw <- function(object, pred_raw) {
  predict_rpf_class_of_prob(object, predict_rpf_prob_raw(object, pred_raw))
}

# The classes of the probabilities `pred_prob` (as predict_rpf_prob() returns them)
predict_rpf_class_of_prob <- function(object, pred_prob) {
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])

  # For each instance, class with higher probability; ties go to the first
  # level so predictions are deterministic and leave the RNG untouched
  pred_class <- factor(
    outcome_levels[max.col(as.matrix(pred_prob), ties.method = "first")],
    levels = outcome_levels
  )
  out <- hardhat::spruce_class(pred_class)

  out
}

# The output of type `type` from the scores `raw` (n x value_size) of a classification model built on the forest
# `object` (an rpf): its pruned models. "link" is "numeric" for every loss: the scores.
predict_rpf_scores <- function(object, raw, type) {
  raw <- as.matrix(raw)
  switch(
    type,
    numeric = ,
    link = predict_rpf_numeric_raw(object, raw),
    class = predict_rpf_class_raw(object, raw),
    prob = predict_rpf_prob_raw(object, raw)
  )
}
