#' Order factor levels by response
#'
#' Regression: Order by mean(y)
#' Binary classification: Order by mean(y==1)
#' Multiclass: Order by first principal component of the weighted covariance
#' matrix of the contingency table
#'
#' See https://doi.org/10.7717/peerj.6339 for details.
#'
#' @param x Factor variable to order
#' @param y Response
#'
#' @return Re-ordered (ordered) factor
#' @noRd
order_factor_by_response <- function(x, y) {
  if (is.numeric(y)) {
    # Regression: Order by mean(y)
    means <- sapply(levels(x), function(lev) {
      mean(y[x == lev])
    })
  } else if (is.factor(y)) {
    if (nlevels(y) <= 2) {
      # Binary classification: Order by mean(y==1)
      means <- sapply(levels(x), function(lev) {
        mean((as.numeric(y) - 1)[x == lev])
      })
    } else {
      # Multiclass: Order by first principal component of the weighted
      # covariance matrix of the contingency table
      means <- pca_order(x, y)
    }
  } else {
    stop(paste(
      "Ordering of factor columns only implemented for regression",
      "and classification outcomes."
    ))
  }

  levels_ordered <- as.character(levels(x)[order(means)])
  factor(x, levels = levels_ordered, ordered = TRUE, exclude = NULL)
}

#' PCA-ordering of factors
#'
#' Order factor levels by first principal component of the weighted covariance
#' matrix of the contingency table
#'
#' @param x Factor variable to order
#' @param y Response
#'
#' @return Order of factor levels
#'
#' @references
#' Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999)
#' Partitioning Nominal Attributes in Decision Trees.
#' Data Min Knowl Discov 3:197. \doi{10.1023/A:1009869804967}.
#' @noRd
pca_order <- function(x, y) {
  if (nlevels(x) < 2) {
    return(seq(1, nlevels(x)))
  }

  # Create contingency table of the nominal outcome with the nominal covariate
  N <- table(x, droplevels(y))

  # PCA of weighted covariance matrix of class probabilites
  P <- N / rowSums(N)
  S <- stats::cov.wt(P, wt = rowSums(N))$cov
  pc1 <- stats::prcomp(S, rank. = 1)$rotation
  score <- P %*% pc1

  # Return ordered factor levels
  order(score)
}

# Sort factor predictors by outcome and re-encode as integer
# save original factor levels for prediction step
# Used in rpf_bridge()
#' @importFrom data.table .SD ':=' as.data.table
preprocess_predictors_fit <- function(processed) {
  predictors <- as.data.table(processed$predictors)

  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }

  # Factor predictors: Order by response (see https://doi.org/10.7717/peerj.6339)
  factor_cols <- names(which(sapply(predictors, is.factor)))
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(
      .SD, order_factor_by_response,
      y = processed$outcomes[[1]]
    ), .SDcols = factor_cols]
  }

  # Save re-ordered factor levels
  factor_levels <- hardhat::get_levels(predictors)

  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }

  list(
    factor_levels = factor_levels,
    predictors_matrix = as.matrix(predictors)
  )
}

# Sort factor predictors using stored level information
# Used in predict_rpf_bridge()
preprocess_predictors_predict <- function(object, predictors) {
  predictors <- as.data.table(predictors)

  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }

  # Re-order factor levels according to saved order
  factor_cols <- names(object$factor_levels)
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := Map(
      factor, .SD, object$factor_levels,
      ordered = TRUE
    ), .SDcols = factor_cols]
  }

  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }

  as.matrix(predictors)
}

# Check outcome to detect mode (classif/regr), return suitable outcome
# Used in rpf_impl()
# Loss is need to transform 1/0 to 1/-1 for exponential
#' @importFrom stats model.matrix
preprocess_outcome <- function(processed, loss) {
  outcomes <- processed$outcomes[[1]]

  # Task type detection: Could be more concise
  is_binary <- length(unique(outcomes)) == 2
  is_integerish <- checkmate::test_integerish(outcomes, any.missing = FALSE)
  is_factor <- checkmate::test_factor(outcomes, any.missing = FALSE)
  is_numeric <- checkmate::test_numeric(outcomes, any.missing = FALSE)

  if (is_binary & is_integerish) {
    warning(paste(
      "y is binary integer, assuming regression task.",
      "Recode y to a factor for classification."
    ))
  }

  if (is_factor) {
    mode <- "classification"

    if (is_binary) {
      # Binary case: Convert to 0, 1 integer
      outcomes <- as.integer(outcomes) - 1L
      # rpf_impl expects Y to be a matrix
      outcomes <- as.matrix(outcomes, ncol = 1)
    } else { # Multiclass
      # One-hot encoding
      outcomes_mat <- stats::model.matrix(~ -1 + outcomes)
      colnames(outcomes_mat) <- levels(outcomes)

      outcomes <- outcomes_mat
    }

    # Exponential expects outcome in 1/-1
    if (loss == "exponential") {
      outcomes[outcomes == 0] <- -1
    }

  } else if (is_numeric) {
    mode <- "regression"
    # rpf_impl expects Y to be a matrix
    outcomes <- as.matrix(outcomes, ncol = 1)
  } else {
    # mode <- "unsupported"
    stop("y should be either numeric (regression) or factor (classification)")
  }

  list(
    outcomes = outcomes,
    mode = mode
  )
}

# Softmax for multiclass prediction using centered logsumexp for numerical stability
# should be identical to sigmoid for scalar input
# Works on matrix input where one row is assumed to be one vector of predictions
# for a single observations
softmax <- function(x) {
  rowmax <- apply(x, 1, max)
  lse <- rowmax + log(rowSums(exp(x - rowmax)))
  exp(x - lse)
}

#' Get remainders where max_interaction requested in `predict_components` is smaller than
#' `max_interaction` set in `rpf`.
#' This is somewhat cumbersome unfortunately, and presumably will have to be partially
#' repeated for other methods in `glex`
#' @noRd
#' @keywords internal
#' @param m Components as calculated in `predict_components`
#' @param levels Outcome levels as stored in `rpf$blueprint$ptypes$outcomes`.
#' @param pred Regular model predictions as returned by `predict.rpf`.
#' @param intercept Intercept as stored in output of `predict_components`.
calc_remainders_multiclass <- function(m, levels, pred, intercept) {
  m <- data.table::as.data.table(m)
  pred <- data.table::as.data.table(pred)
  intercept <- data.table::as.data.table(intercept)

  m[, ".id" := .I]
  pred[, ".id" := .I]
  intercept[, ".id" := .I]

  split_names <- function(mn, split_string = "__class:", target_index = 2) {
    vapply(mn, function(x) {
      unlist(strsplit(x, split = split_string, fixed = TRUE))[[target_index]]
    }, character(1), USE.NAMES = FALSE)
  }

  m_long <- data.table::melt(m, id.vars = ".id", value.name = "m",
                             variable.name = "term", variable.factor = FALSE)

  m_long[, class := split_names(term, split_string = "__class:", target_index = 2)]
  m_sums <- m_long[, list(m_sum = sum(m)), by = c(".id", "class")]

  # long format predictions for calculations
  pred_long <- data.table::melt(pred, id.vars = ".id", value.name = "pred",
                                variable.name = "class", variable.factor = FALSE)
  pred_long[, class := split_names(class, split_string = ".pred_", target_index = 2)]


  # Same game for intercept
  intercept_long <- data.table::melt(intercept, id.vars = ".id", value.name = "intercept",
                                variable.name = "class", variable.factor = FALSE)
  intercept_long[, class := split_names(class, split_string = "__class:", target_index = 2)]

  # merge pred with intercept
  pred_long <- pred_long[intercept_long, on = c(".id", "class")]
  # merge predictions with sum(m)
  merged <- pred_long[m_sums, on = c(".id", "class")]
  # remainder is prediction - sum(m) - intercept, can be positive or negative
  merged[, remainder := pred - m_sum - intercept]
  # Return one column per outcome class
  merged_wide <- data.table::dcast(merged, .id  ~ class, value.var = "remainder")
  # Remove intermediate id column
  merged_wide[, ".id" := NULL]
  # sort columns to original outcome level order, consistent with predict() order
  data.table::setcolorder(merged_wide, levels)

  merged_wide
}
