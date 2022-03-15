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
    stop(paste("Ordering of factor columns only implemented for regression",
               "and classification outcomes."))
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
pca_order <- function(x, y) {
  if (nlevels(x) < 2) {
    return(seq(1, nlevels(x)))
  }
  
  # Create contingency table of the nominal outcome with the nominal covariate
  N <- table(x, droplevels(y))
  
  # PCA of weighted covariance matrix of class probabilites
  P <- N/rowSums(N)
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
      .SD, order_factor_by_response, y = processed$outcomes[[1]]
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
      factor, .SD, object$factor_levels, ordered = TRUE
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
preprocess_outcome <- function(processed) {
  outcomes <- processed$outcomes[[1]]
  
  # Task type detection: Could be more concise
  is_binary <- length(unique(outcomes)) == 2
  is_integerish <- checkmate::test_integerish(outcomes, any.missing = FALSE)
  
  if (is_binary & is_integerish) {
    warning("y is binary integer, assuming classification task")
    outcomes <- factor(outcomes)
  }
  
  is_factor <- checkmate::test_factor(outcomes, any.missing = FALSE)
  is_numeric <- checkmate::test_numeric(outcomes, any.missing = FALSE)
  
  if (is_factor) {
    mode <- "classification"
    outcomes <- as.integer(outcomes) - 1L
  } else if (is_numeric) {
    mode <- "regression"
  } else {
    mode <- "unsupported"
  }
  
  list(
    outcomes = outcomes,
    mode = mode
  )
}
