# The adaptive fit of rpf(splits = NULL): the number of splits chosen out of bag while growing.
#
# The families grow to checkpoints (four per doubling of the splits, up to a cap). At each checkpoint the C++ fit
# (fit_alone() in oob_path.cpp) passes the OOB path there to the stopper made here: the families' predictions at their
# OOB rows, summed over the families of even and of odd index, as rpf(oob_path = ) sums them. The stopper computes the
# forest's corrected OOB error and decides whether growing goes on; the C++ fit then returns the forest at the
# checkpoint it chose.

# The engine of rpf(splits = NULL) (regression, not deterministic): grows to `checkpoints` (the last one is the cap)
# and stops once the lowest corrected OOB error is `patience` doublings of the splits back (round(4 * patience)
# checkpoints: four per doubling by default) and at least `floor` splits are grown, or at the cap. The forest at the
# lowest error is rebuilt from the records of its splits, so it is the forest of rpf() with that many splits and the
# same seed. Returns that forest (the C++ object), the chosen number of splits, the OOB path's table of the checkpoints
# grown, where growing stopped, and whether it stopped at the cap with the error still falling (with a warning).
tune_splits_oob <- function(
  Y,
  X,
  params,
  checkpoints = oob_path_checkpoints(25600L),
  patience = 1,
  floor = 0,
  call = rlang::caller_env()
) {
  pause_gb <- getOption("randomPlantedForest.adaptive_pause_gb", 1)
  if (!checkmate::test_number(pause_gb, lower = 0, finite = TRUE)) {
    cli::cli_abort(
      "{.code options(randomPlantedForest.adaptive_pause_gb)} must be a number of GB, at least 0.",
      call = call
    )
  }
  y <- as.numeric(Y[, 1])
  patience_cp <- max(1L, as.integer(round(4 * patience)))
  st <- oob_stopper(y, params$ntrees, checkpoints, patience = patience_cp, floor = floor)
  params$splits <- max(checkpoints)
  fit <- rpf_impl(
    Y,
    X,
    mode = "regression",
    params = params,
    adaptive = list(
      checkpoints = checkpoints,
      stopper = st$step,
      patience = patience_cp,
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

# The stopper. y: the training response; ntrees: the number of families; checkpoints: the split counts (the last one
# is the cap); patience: stop once the lowest error is `patience` checkpoints back; floor: but not before `floor`
# splits. step(c, path) is called at checkpoints c = 1, 2, ... with the OOB path there (list(checkpoints, sum0, sum1,
# count) with one column, as `$fit$oob_path()` gives it) and returns 0 to grow on, or the chosen checkpoint (its
# index) to stop: once the patience and the floor allow it, or at the last checkpoint. result() gives the chosen number
# of splits, the OOB path's table of the checkpoints grown (as oob_path_table() makes it), where growing stopped, and
# whether it stopped at the cap with the error still falling (its lowest value less than `patience` checkpoints back).
oob_stopper <- function(y, ntrees, checkpoints, patience = 4L, floor = 0) {
  y <- as.numeric(y)
  C <- length(checkpoints)
  st <- new.env(parent = emptyenv())
  st$c <- 0L
  st$sel <- NA_integer_
  st$capped <- FALSE
  st$no_oob <- FALSE
  st$points <- matrix(NA_real_, C, 4)
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
    st$points[c, ] <- oob_path_point(path$sum0[, 1], path$sum1[, 1], path$count, y, ntrees)
    best <- oob_select(st$points[seq_len(c), 1])
    last <- c == C
    if (last && c - best < patience) {
      st$capped <- TRUE
    }
    if (last || (c - best >= patience && checkpoints[c] >= floor)) {
      st$sel <- best
      return(best)
    }
    0L
  }
  st$result <- function() {
    g <- seq_len(st$c)
    list(
      splits = as.integer(checkpoints[st$sel]),
      path = oob_path_frame(checkpoints[g], st$points[g, , drop = FALSE], st$n_rows),
      stopped_at = as.integer(checkpoints[st$c]),
      capped = st$capped
    )
  }
  st
}

# The chosen checkpoint among 1..length(mse) (the corrected OOB mean squared errors grown): the lowest; ties (relative
# 1e-12) go to fewer splits
oob_select <- function(mse) {
  which(mse <= min(mse) * (1 + 1e-12) | mse == min(mse))[1]
}
