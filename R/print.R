#' Print an rpf fit
#'
#' @param x And object of class `rpf`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly: `x`.
#' @seealso [`rpf`].
#' @export
#'
#' @examples
#' rpf(mpg ~., data = mtcars, max_interaction = 2, ntrees = 10)
print.rpf <- function(x, ...) {
  model_formula <- sub("\\s\\+ 0$", "" , deparse(x$blueprint$formula))
  mode <- switch(x$mode, regression = "Regression", classification = "Classification")

  cat("--", mode, "Random Planted Forest --\n\n")

  # Print interaction level (all, main only, n-th degree)
  p <- length(x$blueprint$ptypes$predictors)
  # if max_interaction is p, it's equivalent to 0 -> all interactions
  full <- p == x$params$max_interaction
  maxint <- ifelse(full, 0, x$params$max_interaction)
  degree <- switch(
    # numeric expression doesn't allow default value, so characterify
    as.character(maxint),
    "0" = "all possible interactions.\n",
    "1" = "main effects only.\n",
    paste0(maxint, "-degree interactions.\n")
   )

  cat("Formula:", model_formula, "\n")
  cat("Fit using", p, "predictors and", degree)

  purification <- ifelse(x$fit$is_purified(), "is", "is _not_")
  cat("Forest ", purification, "purified!\n\n")

  param_names <- names(x$params)
  nm_lengths <- nchar(param_names)

  cat("Called with parameters:\n\n")
  for (i in seq_along(x$params)) {
    cat(sprintf(paste0(" %", max(nm_lengths), "s: %s\n"), param_names[[i]], x$params[[i]]))
  }

  invisible(x)
}
