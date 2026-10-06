#' Purify a Random Planted Forest
#'
#' Purification makes the components of the forest's functional decomposition
#' unique, which [predict_components()] relies on. Unless [rpf()] was called
#' with `purify = TRUE`, [predict_components()] purifies the forest on first
#' use.
#'
#' The forest is modified in place: `x` and every copy of it are purified,
#' whether or not the result is assigned.
#'
#' @param x An object of class `rpf`.
#' @param ... Not currently used, but required for extensibility. Unknown
#'   arguments are an error.
#'
#' @return `purify()` returns `x` invisibly. `is_purified()` returns `TRUE` or
#'   `FALSE`.
#' @export
#'
#' @examples
#' rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, ntrees = 10)
#' is_purified(rpfit)
#' purify(rpfit)
#' is_purified(rpfit)
purify <- function(x, ...) {
  UseMethod("purify")
}

#' @export
#' @rdname purify
purify.default <- function(x, ...) {
  cli::cli_abort("{.fn purify} is not defined for a {.cls {class(x)[1]}}.")
}

#' @param maxp_interaction `[NULL]`: Highest interaction order to purify.
#'   Higher-order components are set to zero, but still influence lower orders
#'   during purification. `NULL` purifies all orders.
#' @param mode `[2]`: Purification algorithm: `2` is the fast exact KD-tree
#'   based algorithm, `1` the original grid-based one.
#' @param nthreads `[NULL]`: Number of threads. `NULL` uses the `nthreads` the
#'   forest was fitted with, capped at the available cores.
#' @export
#' @rdname purify
#' @importFrom rlang %||%
purify.rpf <- function(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL) {
  rlang::check_dots_empty()
  checkmate::assert_class(x, "rpf")
  check_rpf_alive(x)
  checkmate::assert_int(maxp_interaction, lower = 1, null.ok = TRUE)
  checkmate::assert_int(mode, lower = 1, upper = 2)
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  # 0 tells C++ to use all orders / the forest's own nthreads
  x$fit$purify_threads(
    as.integer(maxp_interaction %||% 0L),
    as.integer(nthreads %||% 0L),
    as.integer(mode)
  )
  invisible(x)
}

#' Check if a forest is purified
#' @export
#' @rdname purify
is_purified <- function(x) {
  checkmate::assert_class(x, "rpf")
  check_rpf_alive(x)
  x$fit$is_purified()
}
