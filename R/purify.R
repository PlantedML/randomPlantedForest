#' Purify a Random Planted Forest
#'
#' Purifies an rpf object.
#'
#' Unless [`rpf()`] is called with `purify = TRUE`, the forest has to be purified after fit
#' to ensure the components extracted by [`predict_components()`] are valid.
#' [`predict_components()`] will automatically purify a forest if [`is_purified()`] reports `FALSE`.
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
    "`purify()` is not defined for a '",
    class(x)[1],
    "'.",
    call. = FALSE
  )
}

#' @param maxp_interaction integer or NULL: Only compute/store purified components
#'   up to this interaction order. Higher-order purified trees are zeroed (not
#'   computed) but still implicitly influence lower orders during purification.
#'   If NULL, purify all orders (default behavior).
#' @param mode integer(1): Purification algorithm mode. 1 = legacy grid path
#'   used by `fit$fit$purify()`; 2 = fast exact KD-tree based path. Defaults to 2.
#' @param nthreads integer or NULL: number of threads to use. If NULL, defaults
#'   to min of the object's configured `nthreads` and available threads.
#' @export
#' @rdname purify
#' @importFrom utils capture.output
purify.rpf <- function(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL) {
  checkmate::assert_class(x, "rpf")
  check_rpf_alive(x)
  checkmate::assert_int(mode, lower = 1, upper = 2)
  if (!is.null(nthreads)) {
    checkmate::assert_int(nthreads, lower = 1)
  }
  if (is.null(maxp_interaction)) {
    # Default: exact cut points, full interaction order
    x$fit$purify_threads(0L, as.integer(if (is.null(nthreads)) 0L else nthreads), as.integer(mode))
  } else {
    checkmate::assert_int(maxp_interaction, lower = 1)
    x$fit$purify_threads(
      as.integer(maxp_interaction),
      as.integer(if (is.null(nthreads)) 0L else nthreads),
      as.integer(mode)
    )
  }
  x
}

#' Check if a forest is purified
#' @export
#' @rdname purify
is_purified <- function(x) {
  checkmate::assert_class(x, "rpf")
  check_rpf_alive(x)
  x$fit$is_purified()
}
