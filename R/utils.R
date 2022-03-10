

#' Order factor column by response
#' 
#' Regression: Order by mean(y)
#' Binary classification: Order by mean(y==1)
#' Multiclass: Order by first principal component of the weighted covariance matrix of the contingency table
#' Survival: Order by median survival if available or largest quantile available in all strata if median not available
#'
#' See https://doi.org/10.7717/peerj.6339 for details
#' 
#' @param x Factor variable to order
#' @param y Response
#'
#' @return Re-ordered (ordered) factor
order_factor_by_response <- function(x, y) {
  means <- sapply(levels(x), function(lev) {
    mean(y[x == lev])
  })
  levels_ordered <- as.character(levels(x)[order(means)])
  factor(x, levels = levels_ordered, ordered = TRUE, exclude = NULL)
}
