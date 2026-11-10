# The adaptive fit of rpf(splits = NULL): the number of splits chosen out of bag while growing.
#
# The families grow to checkpoints (four per doubling of the splits, up to a cap). At each checkpoint the C++ fit
# (fit_alone() in oob_path.cpp) passes the OOB path there to the stopper made here: the families' predictions at their
# OOB rows, summed over the families of even and of odd index, as rpf(oob_path = ) sums them. The stopper computes the
# forest's corrected OOB error (for classification its criterion: criterion.R) and decides whether growing goes on;
# the C++ fit then returns the forest at the checkpoint it chose.

# The engine of rpf(splits = NULL) (not deterministic): grows to `checkpoints` (the last one is the cap) and stops once
# the lowest corrected OOB error is `patience` doublings of the splits back (round(4 * patience) checkpoints: four per
# doubling by default), or once it has improved by at most `tolerance` paired standard errors over the last `patience`
# doublings (a plateau; tolerance 0: no plateau stop), and at least `floor` splits are grown; or at the cap. The forest
# at the lowest error is rebuilt from the records of its splits, so it is the forest of rpf() with that many splits
# and the same seed. Returns that forest (the C++ object), the chosen number of splits, the OOB path's table of the
# checkpoints grown, where growing stopped, and whether it stopped at the cap with the error still falling (with a
# warning). Y: the C++ outcome matrix (rpf_impl()); score_y: what the error is (oob_score_y(): the numeric outcome, the
# default Y[, 1], or a classification criterion, then the forest's own criterion: the Brier score or the log-loss).
tune_splits_oob <- function(
  Y,
  X,
  params,
  checkpoints = oob_path_checkpoints(25600L),
  patience = 1,
  tolerance = getOption("randomPlantedForest.adaptive_tolerance", 1),
  floor = 0,
  mode = "regression",
  score_y = NULL,
  call = rlang::caller_env()
) {
  pause_gb <- getOption("randomPlantedForest.adaptive_pause_gb", 1)
  if (!checkmate::test_number(pause_gb, lower = 0, finite = TRUE)) {
    cli::cli_abort(
      "{.code options(randomPlantedForest.adaptive_pause_gb)} must be a number of GB, at least 0.",
      call = call
    )
  }
  if (!checkmate::test_number(tolerance, lower = 0, finite = TRUE)) {
    cli::cli_abort(
      "{.code options(randomPlantedForest.adaptive_tolerance)} must be a number of standard errors, at least 0.",
      call = call
    )
  }
  y <- score_y %||% as.numeric(Y[, 1])
  patience_cp <- max(1L, as.integer(round(4 * patience)))
  st <- oob_stopper(y, params$ntrees, checkpoints, patience = patience_cp, tolerance = tolerance, floor = floor)
  params$splits <- max(checkpoints)
  fit <- rpf_impl(
    Y,
    X,
    mode = mode,
    params = params,
    adaptive = list(
      checkpoints = checkpoints,
      stopper = st$step,
      patience = patience_cp,
      tolerance = tolerance,
      floor = floor,
      pause_gb = pause_gb
    )
  )
  if (isTRUE(st$no_oob)) {
    cli::cli_abort("No training row is out of bag for any tree family: {.code splits = NULL} needs some.", call = call)
  }
  res <- st$result()
  if (res$capped) {
    cli::cli_warn(
      c(
        "Growing stopped at the cap ({res$stopped_at} splits) while the out-of-bag error was still falling.",
        "i" = "A fit with more {.arg splits} and {.code oob_path = TRUE} shows whether more splits help."
      ),
      call = call
    )
  }
  c(list(fit = fit), res)
}

# The stopper. y: the training response, or a classification criterion (oob_score_y()); ntrees: the number of
# families; checkpoints: the split counts (the last one is the cap); patience: stop once the lowest error is `patience`
# checkpoints back; tolerance: or once the error is on a plateau (below); floor: but not before `floor` splits.
# step(c, path) is called at checkpoints c = 1, 2, ... with the OOB path there (list(checkpoints, sum0, sum1, count)
# with value_size columns, as `$fit$oob_path()` gives it) and returns 0 to grow on, or the chosen checkpoint (its index:
# the lowest error) to stop: once the patience or the plateau and the floor allow it, or at the last checkpoint. The
# error is the corrected OOB mean squared error, or the mean corrected criterion. result() gives the chosen number of
# splits, the OOB path's table of the checkpoints grown (as oob_path_table() makes it), where growing stopped, and
# whether it stopped at the cap with the error still falling (neither the patience nor the plateau reached).
#
# The plateau (tolerance > 0, from checkpoint patience + 1 on): with b the lowest error up to c and a the lowest up to
# c - patience (the first of ties, as oob_select()), the gain e[a] - e[b] over the last `patience` checkpoints is at
# most tolerance * sd(r_a - r_b) / sqrt(n), the standard error of the paired difference of the rows' corrected errors
# at a and b (n: the rows with an OOB family). Where the lowest error did not move (b <= c - patience) the patience
# stops anyway, so tolerance = 0 is the patience alone. The rows' errors are kept only for the checkpoints that can
# still be a or b: the last `patience` ones and the ties of the lowest error up to c - patience.
oob_stopper <- function(y, ntrees, checkpoints, patience = 4L, tolerance = 0, floor = 0) {
  crit <- if (inherits(y, "rpf_crit")) y
  if (is.null(crit)) {
    y <- as.numeric(y)
  }
  C <- length(checkpoints)
  plateau <- tolerance > 0
  st <- new.env(parent = emptyenv())
  st$c <- 0L
  st$sel <- NA_integer_
  st$capped <- FALSE
  st$no_oob <- FALSE
  st$points <- matrix(NA_real_, C, 4)
  st$rows <- vector("list", C)
  st$step <- function(c, path) {
    if (c != st$c + 1L) {
      cli::cli_abort("The adaptive fit's checkpoints came out of order.", .internal = TRUE)
    }
    st$c <- c
    if (c == 1L) {
      st$n_rows <- sum(path$count[, 1] + path$count[, 2] > 0)
      if (st$n_rows == 0) {
        # (reported by tune_splits_oob() after the fit: an error here would reach R through the C++ fit)
        st$no_oob <- TRUE
        st$sel <- 1L
        return(1L)
      }
    }
    pt <- if (is.null(crit)) {
      oob_path_point(path$sum0[, 1], path$sum1[, 1], path$count, y, ntrees, rows = plateau)
    } else {
      oob_path_point_crit(path$sum0, path$sum1, path$count, crit, ntrees, rows = plateau)
    }
    if (plateau) {
      st$rows[[c]] <- pt$rows
      pt <- pt$point
    }
    st$points[c, ] <- pt
    e <- st$points[seq_len(c), 1]
    best <- oob_select(e)
    done <- c - best >= patience
    if (plateau && c > patience) {
      old <- seq_len(c - patience)
      if (!done) {
        a <- oob_select(e[old])
        se <- stats::sd(st$rows[[a]] - st$rows[[best]]) / sqrt(st$n_rows)
        done <- e[a] - e[best] <= tolerance * se
      }
      # the rows of the checkpoints that can no longer be a or b are freed: all but the ties of the lowest among
      # e[old] (oob_lowest(), as oob_select() chooses a)
      st$rows[old[-oob_lowest(e[old])]] <- list(NULL)
    }
    last <- c == C
    if (last && !done) {
      st$capped <- TRUE
    }
    if (last || (done && checkpoints[c] >= floor)) {
      st$sel <- best
      return(best)
    }
    0L
  }
  st$result <- function() {
    g <- seq_len(st$c)
    list(
      splits = as.integer(checkpoints[st$sel]),
      path = oob_path_frame(checkpoints[g], st$points[g, , drop = FALSE], st$n_rows, crit),
      stopped_at = as.integer(checkpoints[st$c]),
      capped = st$capped
    )
  }
  st
}

# The chosen checkpoint among 1..length(mse) (the corrected OOB mean squared errors grown, or mean criteria): the
# lowest by oob_lowest(); ties go to fewer splits
oob_select <- function(mse) {
  oob_lowest(mse)[1]
}

# The lowest of the corrected out-of-bag errors `x` (mean squared errors, RMSEs or mean criteria), as the adaptive stop
# (oob_select()), prune_components()' out-of-bag choice (prune_choose()) and rpf(tune = ) (tune_oob_choose()) choose:
# the indices of its ties (relative 1e-12), in the order of `x`. A corrected error at or below 0 is a correction that
# overshot (the tables show it clamped at 0), not a perfect model: such values are left out unless all are, and then
# they all tie.
oob_lowest <- function(x) {
  if (all(x <= 0)) {
    return(seq_along(x))
  }
  ok <- which(x > 0)
  ok[x[ok] <= min(x[ok]) * (1 + 1e-12)]
}
