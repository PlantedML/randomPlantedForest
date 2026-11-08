#' Checkpoints of the out-of-bag path
#'
#' The numbers of splits at which [rpf()] with `oob_path = TRUE` evaluates the
#' out-of-bag (OOB) error: a geometric grid with four points per doubling,
#' `round(25 * 2^(j / 4))` from 6 up to `splits` (it contains 12, 25, 50, ...,
#' 1600, 3200), and `splits` itself; every number up to 12 if `splits <= 12`.
#' The error changes slowly in `log(splits)`, so a geometric grid resolves it
#' evenly, and the path's memory, `2 * nrow(x) * length(checkpoints)` doubles,
#' stays small; every split would cost `2 * nrow(x) * splits`.
#'
#' @param splits `[integer(1)]`: The largest number of splits (the `splits` of
#'   the fit).
#'
#' @return An increasing integer vector.
#' @seealso [rpf()] (`oob_path`)
#' @export
#'
#' @examples
#' oob_path_checkpoints(400)
oob_path_checkpoints <- function(splits) {
  checkmate::assert_int(splits, lower = 1)
  if (splits <= 12) {
    return(seq_len(splits))
  }
  j <- seq(-8, ceiling(4 * log2(splits / 25)))
  cp <- round(25 * 2^(j / 4))
  as.integer(sort(unique(c(cp[cp <= splits], splits))))
}

# The table of a fit's OOB path (rpf()'s `$oob_path`), from the sums the fit tracked (`fit$oob_path()`: checkpoints,
# sum0 / sum1 = each row's OOB predictions summed over the families of even / odd index, one column per checkpoint,
# count = the number of OOB families of each half). One row per checkpoint: splits, oob_rmse (corrected), oob_rmse_raw,
# and the standard errors of the two mean squared errors, se_mse and se_mse_raw. Attribute n_rows: the number of rows
# with an OOB family. The correction is the one of prune_components()' OOB choice (see rpf()'s details).
oob_path_table <- function(path, y, ntrees, call = rlang::caller_env()) {
  y <- as.numeric(y)
  m0 <- path$count[, 1]
  m1 <- path$count[, 2]
  m <- m0 + m1
  ok <- m > 0
  both <- m0 > 0 & m1 > 0
  n_ok <- sum(ok)
  if (n_ok == 0) {
    cli::cli_abort("No training row is out of bag for any tree family.", call = call)
  }
  cp <- path$checkpoints
  out <- data.frame(
    splits = as.integer(cp),
    oob_rmse = NA_real_,
    oob_rmse_raw = NA_real_,
    se_mse = NA_real_,
    se_mse_raw = NA_real_
  )
  for (c in seq_along(cp)) {
    s0 <- path$sum0[, c]
    s1 <- path$sum1[, c]
    e2 <- (y - (s0 + s1) / m)^2
    v <- numeric(length(y))
    v[both] <- (s0[both] / m0[both] - s1[both] / m1[both])^2 / (1 / m0[both] + 1 / m1[both])
    corr <- e2 - (1 / m - 1 / ntrees) * v
    out$oob_rmse[c] <- sqrt(max(mean(corr[ok]), 0))
    out$oob_rmse_raw[c] <- sqrt(mean(e2[ok]))
    out$se_mse[c] <- stats::sd(corr[ok]) / sqrt(n_ok)
    out$se_mse_raw[c] <- stats::sd(e2[ok]) / sqrt(n_ok)
  }
  attr(out, "n_rows") <- n_ok
  out
}

# Validate rpf()'s `oob_path`: NULL (off) or the checkpoints
check_oob_path <- function(oob_path, splits, mode, deterministic, call = rlang::caller_env()) {
  if (isFALSE(oob_path)) {
    return(NULL)
  }
  if (isTRUE(oob_path)) {
    checkpoints <- oob_path_checkpoints(splits)
  } else if (
    checkmate::test_integerish(oob_path, lower = 1, upper = splits, any.missing = FALSE, min.len = 1) &&
      !is.unsorted(oob_path, strictly = TRUE)
  ) {
    checkpoints <- as.integer(oob_path)
  } else {
    cli::cli_abort(
      "{.arg oob_path} must be {.code TRUE}, {.code FALSE} or an increasing integer vector of split counts in
      {.code 1..splits}.",
      call = call
    )
  }
  if (mode != "regression") {
    cli::cli_abort("{.arg oob_path} is available for regression only.", call = call)
  }
  if (deterministic) {
    cli::cli_abort(
      "A deterministic forest has no out-of-bag rows: {.arg oob_path} needs {.code deterministic = FALSE}.",
      call = call
    )
  }
  checkpoints
}
