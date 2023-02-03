#' Purify a Random Planted Forest
#'
#' TODO: Explain what this does
#'
#' Unless [`rpf()`] is called with `purify = TRUE`, the forest has to be purified after fit
#' to ensure the components extracted by [`extract_components()`] are valid.
#' [`extract_components()`] will automatically purify a forest if [`is_purified()`] reports `FALSE`.
#'
#' @param x And object of class `rpf`.
#' @param ... (Unused)
#'
#' @return Invisibly: The [`rpf`] object.
#' @export
#'
#' @examples
#' rpfit <- rpf(mpg ~., data = mtcars, max_interaction = 2, ntrees = 10)
#' purify(rpfit)
purify <- function(x, ...) {
  UseMethod("purify")
}

#' @export
#' @rdname purify
purify.default <- function(x, ...) {
  stop(
    "`purify()` is not defined for a '", class(x)[1], "'.",
    call. = FALSE
  )
}

#' @export
#' @rdname purify
#' @importFrom utils capture.output
purify.rpf <- function(x, ...) {
  x$fit$purify()
  x
}

#' Check if a forest is purified
#' @export
#' @rdname purify
is_purified <- function(x) {
  checkmate::assert_class(x, "rpf")
  x$fit$is_purified()
}
