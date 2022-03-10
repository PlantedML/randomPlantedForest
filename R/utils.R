

#' Order factor levels by response
#' 
#' Regression: Order by mean(y)
#' Binary classification: Order by mean(y==1)
#' Multiclass: Order by first principal component of the weighted covariance matrix of the contingency table
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
      # Multiclass: Order by first principal component of the weighted covariance matrix of the contingency table
      means <- pca_order(x, y)
    }
  } else {
    stop("Ordering of factor columns only implemented for regression and classification outcomes.")
  }
  
  levels_ordered <- as.character(levels(x)[order(means)])
  factor(x, levels = levels_ordered, ordered = TRUE, exclude = NULL)
}

# 
# 
#' Order factor levels by first principal component of the weighted covariance matrix of the contingency table
#' 
#' Reference: Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197. \doi{10.1023/A:1009869804967}.
#'
#' @param x Factor variable to order
#' @param y Response
#'
#' @return Order of factor levels
pca_order <- function(x, y) {
  if (nlevels(x) < 2) {
    return(seq(1, nlevels(x)))
  }
  
  # Create contingency table of the nominal outcome with the nominal covariate
  N <- table(x, droplevels(y))
  
  # PCA of weighted covariance matrix of class probabilites
  P <- N/rowSums(N)
  S <- cov.wt(P, wt = rowSums(N))$cov
  pc1 <- prcomp(S, rank. = 1)$rotation
  score <- P %*% pc1
  
  # Return ordered factor levels
  order(score)
}
