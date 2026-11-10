#' Purify a Random Planted Forest
#'
#' Purification makes the components of the forest's functional decomposition
#' unique, which [predict_components()] relies on. Unless [rpf()] was called
#' with `purify = TRUE`, [predict_components()] purifies the forest on first
#' use.
#'
#' The forest is modified in place: `x` and every copy of it are purified,
#' whether or not the result is assigned.
#' If the forest is already purified with all the orders asked for, `purify()` returns it unmodified, whatever
#' `mode` and `nthreads` are. A forest purified with a `maxp_interaction` below the order asked for (all orders
#' with the default `NULL`) is purified again, which needs its training data: a forest restored by
#' [rpf_unmarshal()] without it can only be purified up to its earlier cap.
#'
#' @param x `[rpf]`: A fitted [`rpf`] model.
#' @param maxp_interaction `[integer(1) | NULL: NULL]`: Highest interaction order to purify.
#'   Higher-order components are not computed, but still influence lower orders
#'   during purification, which are the same as when purifying all orders.
#'   [predict_components()] then returns zero for the higher orders and their sum as `remainder`.
#'   [predict()][predict.rpf()] still returns the whole forest's prediction.
#'   `NULL` purifies all orders.
#' @param mode `[integer(1): 2]`: Purification algorithm: `2` is the fast exact KD-tree
#'   based algorithm, which keeps the components as weighted rectangles, `1` the
#'   original grid-based one, which stores every component as a dense grid over
#'   the split points and can run out of memory from `max_interaction = 3` on.
#'   Both give the same components up to rounding.
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the `nthreads` the
#'   forest was fitted with, capped at the available cores.
#' @param ... Reserved for future expansion.
#'
#' @return `purify()` returns `x` invisibly. `is_purified()` returns `TRUE` or
#'   `FALSE`.
#' @export
#' @importFrom rlang %||%
#' @examples
#' rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, ntrees = 10)
#' is_purified(rpfit)
#' purify(rpfit)
#' is_purified(rpfit)
purify <- function(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL) {
  checkmate::assert_class(x, "rpf")
  rlang::check_dots_empty()
  check_rpf_alive(x)
  checkmate::assert_int(maxp_interaction, lower = 1, null.ok = TRUE)
  checkmate::assert_int(mode, lower = 1, upper = 2)
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  if (is_purified(x)) {
    # 0: the stored purification has all orders
    cap <- x$fit$get_purify_cap()
    asked <- maxp_interaction %||% Inf
    if (cap == 0 || asked <= cap) {
      return(invisible(x))
    }
    # purified again below, from the leaves, which purification leaves unchanged
  }
  # 0 tells C++ to use all orders / the forest's own nthreads
  tryCatch(
    x$fit$purify_threads(
      as.integer(maxp_interaction %||% 0L),
      as.integer(nthreads %||% 0L),
      as.integer(mode)
    ),
    error = function(e) {
      if (!is_purified(x) || x$fit$get_purify_cap() == 0) {
        stop(e)
      }
      cli::cli_abort(
        "The forest is purified up to order {x$fit$get_purify_cap()} only and cannot be purified to higher orders.",
        parent = e
      )
    }
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
