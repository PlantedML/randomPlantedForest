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
# sum0 / sum1 = each row's OOB predictions summed over the families of even / odd index, value_size columns per
# checkpoint, count = the number of OOB families of each half). One row per checkpoint: splits, oob_rmse (corrected),
# oob_rmse_raw, and the standard errors of the two mean squared errors, se_mse and se_mse_raw. Attribute n_rows: the
# number of rows with an OOB family. The correction is the one of prune_components()' OOB choice (see rpf()'s details),
# except for the rows whose OOB families are all even or all odd: they are left uncorrected here, where the pruning
# gives them the mean variance estimate of the rows with both halves. So the path's corrected error of a forest is
# higher than the one the pruning gives its unpruned model (rule D at tau = 0) by those rows' correction, by more
# with fewer families.
# y: the numeric outcome, or a classification criterion (oob_score_y()), whose table has the mean corrected criterion
# instead (oob_path_frame()).
oob_path_table <- function(path, y, ntrees, call = rlang::caller_env()) {
  n_ok <- sum(path$count[, 1] + path$count[, 2] > 0)
  if (n_ok == 0) {
    cli::cli_abort("No training row is out of bag for any tree family.", call = call)
  }
  cp <- path$checkpoints
  points <- matrix(NA_real_, length(cp), 4)
  if (inherits(y, "rpf_crit")) {
    V <- y$V
    for (c in seq_along(cp)) {
      cols <- (c - 1L) * V + seq_len(V)
      points[c, ] <- oob_path_point_crit(
        path$sum0[, cols, drop = FALSE],
        path$sum1[, cols, drop = FALSE],
        path$count,
        y,
        ntrees
      )
    }
    return(oob_path_frame(cp, points, n_ok, y))
  }
  y <- as.numeric(y)
  for (c in seq_along(cp)) {
    points[c, ] <- oob_path_point(path$sum0[, c], path$sum1[, c], path$count, y, ntrees)
  }
  oob_path_frame(cp, points, n_ok)
}

# One checkpoint of the OOB path, over the rows with an OOB family: the corrected and the plain OOB mean squared
# error and their standard errors, from s0 / s1, the rows' OOB predictions summed over the families of even / odd
# index, and count, the number of OOB families of each half (see rpf()'s details). rows = TRUE: a list of that (point)
# and the rows' corrected errors (rows; over the rows with an OOB family), for the adaptive fit's plateau.
oob_path_point <- function(s0, s1, count, y, ntrees, rows = FALSE) {
  m0 <- count[, 1]
  m1 <- count[, 2]
  m <- m0 + m1
  ok <- m > 0
  both <- m0 > 0 & m1 > 0
  n_ok <- sum(ok)
  e2 <- (y - (s0 + s1) / m)^2
  v <- numeric(length(y))
  v[both] <- (s0[both] / m0[both] - s1[both] / m1[both])^2 / (1 / m0[both] + 1 / m1[both])
  corr <- e2 - (1 / m - 1 / ntrees) * v
  point <- c(mean(corr[ok]), mean(e2[ok]), stats::sd(corr[ok]) / sqrt(n_ok), stats::sd(e2[ok]) / sqrt(n_ok))
  if (rows) list(point = point, rows = corr[ok]) else point
}

# oob_path_point() for a classification criterion `crit` (oob_score_y()): the mean corrected and plain criterion over
# the rows with an OOB family and their standard errors, from s0 / s1 (n x V: the sums of the families of even / odd
# index) and count. Rows whose OOB families are all even or all odd are not corrected (as for regression). rows: as
# for oob_path_point().
oob_path_point_crit <- function(s0, s1, count, crit, ntrees, rows = FALSE) {
  m0 <- count[, 1]
  m1 <- count[, 2]
  m <- m0 + m1
  ok <- m > 0
  n_ok <- sum(ok)
  s0 <- s0[ok, , drop = FALSE]
  s1 <- s1[ok, , drop = FALSE]
  r <- crit_rows(
    crit,
    which(ok),
    (s0 + s1) / m[ok],
    s0 / pmax(m0[ok], 1),
    s1 / pmax(m1[ok], 1),
    m0[ok],
    m1[ok],
    ntrees,
    impute = FALSE
  )
  point <- c(mean(r$corr), mean(r$e), stats::sd(r$corr) / sqrt(n_ok), stats::sd(r$e) / sqrt(n_ok))
  if (rows) list(point = point, rows = r$corr) else point
}

# The table from the checkpoints' oob_path_point() rows. With a classification criterion `crit`, the columns are the
# mean corrected and plain criterion (oob_loss, oob_loss_raw) and their standard errors (se_loss, se_loss_raw), and the
# attribute `criterion` names it ("brier" or "logloss").
oob_path_frame <- function(checkpoints, points, n_rows, crit = NULL) {
  if (!is.null(crit)) {
    out <- data.frame(
      splits = as.integer(checkpoints),
      oob_loss = points[, 1],
      oob_loss_raw = points[, 2],
      se_loss = points[, 3],
      se_loss_raw = points[, 4]
    )
    attr(out, "n_rows") <- n_rows
    attr(out, "criterion") <- crit$type
    return(out)
  }
  out <- data.frame(
    splits = as.integer(checkpoints),
    oob_rmse = sqrt(pmax(points[, 1], 0)),
    oob_rmse_raw = sqrt(points[, 2]),
    se_mse = points[, 3],
    se_mse_raw = points[, 4]
  )
  attr(out, "n_rows") <- n_rows
  out
}

# Validate rpf()'s `oob_path`: NULL (off) or the checkpoints
check_oob_path <- function(oob_path, splits, deterministic, call = rlang::caller_env()) {
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
  if (deterministic) {
    cli::cli_abort(
      "A deterministic forest has no out-of-bag rows: {.arg oob_path} needs {.code deterministic = FALSE}.",
      call = call
    )
  }
  checkpoints
}
