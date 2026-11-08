#' The greediness ladder of the out-of-bag tuning
#'
#' The points of greediness that `rpf(tune = )` searches: pairs of `split_try` (thresholds evaluated per split
#' candidate) and `t_try` (share of the possible splits drawn as candidates in each round), greediest first. Each step
#' down the ladder halves `t_try` and `split_try` (rounded up, at least 1), from 40 / 1 to 1 / 0.003125. Row 3 is
#' [rpf()]'s default, 10 / 0.4.
#'
#' @return A `data.frame` with columns `split_try` and `t_try`, one row per point, greediest first.
#' @seealso [tune_control()], [rpf()]
#' @export
#' @examples
#' greediness_ladder()
greediness_ladder <- function() {
  data.frame(
    split_try = c(40L, 20L, 10L, 5L, 3L, 2L, 1L, 1L, 1L, 1L),
    t_try = c(1, 0.8, 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125)
  )
}

#' Settings of the out-of-bag tuning
#'
#' What `rpf(tune = )` searches. The search is coordinate-wise, so that it needs few forests (9 to 11 with the
#' defaults):
#' 1. the points `search` of `ladder` at the first `min_leaf_rows`;
#' 2. every `min_leaf_rows` at the point chosen in step 1 (skipped if only one is given);
#' 3. the two ladder neighbours of that point at the `min_leaf_rows` chosen in step 2 (if `neighbours`).
#'
#' Every forest has its number of splits chosen out of bag while it grows (as `rpf(splits = NULL)`) and is pruned
#' by [prune_components()] with its defaults, without the display. A step chooses, among all forests grown at the
#' coordinate it varies, the one whose pruned model has the lowest corrected out-of-bag RMSE (ties: the greedier
#' point, the larger `min_leaf_rows`). Step 3 chooses among every forest grown at the chosen `min_leaf_rows`, so with
#' `min_leaf_rows` 1 it also compares step 1's.
#'
#' `rpf(tune = TRUE)` is `tune_control()`, and `rpf(tune = "quick")` is
#' `tune_control(search = c(5, 7), min_leaf_rows = 1, neighbours = FALSE)`: two forests, 3 / 0.1 and 1 / 0.025.
#'
#' @param ladder `[data.frame: greediness_ladder()]`: The points of greediness, greediest first: columns `split_try`
#'   (whole numbers, at least 1) and `t_try` (in `(0, 1]`).
#' @param search `[integer: 3:8]`: The rows of `ladder` searched in step 1 (on the default ladder from 10 / 0.4 to
#'   1 / 0.0125).
#' @param min_leaf_rows `[integer: c(1, 5, 20, 100)]`: The values of [rpf()]'s `min_leaf_rows` searched in step 2;
#'   the first one is that of steps 1. Values above 1 need `split_structure = "leaves"`.
#' @param neighbours `[logical(1): TRUE]`: Whether to run step 3.
#'
#' @return An object of class `rpf_tune_control`, a list with the four settings, for `rpf(tune = )`.
#' @seealso [greediness_ladder()], [rpf()]
#' @export
#' @examples
#' # The full search and the quick preset
#' tune_control()
#' tune_control(search = c(5, 7), min_leaf_rows = 1, neighbours = FALSE)
#'
#' \donttest{
#' # A small search: two points of greediness, then two leaf sizes
#' set.seed(1)
#' tuned <- rpf(
#'   mpg ~ .,
#'   data = mtcars,
#'   ntrees = 10,
#'   tune = tune_control(search = c(3, 5), min_leaf_rows = c(1, 3), neighbours = FALSE)
#' )
#' tuned$tune$search
#' }
tune_control <- function(
  ladder = greediness_ladder(),
  search = 3:8,
  min_leaf_rows = c(1, 5, 20, 100),
  neighbours = TRUE
) {
  checkmate::assert_data_frame(ladder, min.rows = 1)
  if (!all(c("split_try", "t_try") %in% names(ladder))) {
    cli::cli_abort("{.arg ladder} needs the columns {.code split_try} and {.code t_try}.")
  }
  checkmate::assert_integerish(ladder$split_try, lower = 1, any.missing = FALSE, .var.name = "ladder$split_try")
  checkmate::assert_numeric(ladder$t_try, lower = 0, upper = 1, any.missing = FALSE, .var.name = "ladder$t_try")
  if (any(ladder$t_try <= 0)) {
    cli::cli_abort("{.code ladder$t_try} must be positive.")
  }
  checkmate::assert_integerish(search, lower = 1, any.missing = FALSE, min.len = 1, unique = TRUE)
  if (any(search > nrow(ladder))) {
    cli::cli_abort("{.arg search} must be rows of {.arg ladder}, which has {nrow(ladder)} row{?s}.")
  }
  checkmate::assert_integerish(min_leaf_rows, lower = 1, any.missing = FALSE, min.len = 1, unique = TRUE)
  checkmate::assert_flag(neighbours)
  structure(
    list(
      ladder = data.frame(split_try = as.integer(ladder$split_try), t_try = as.numeric(ladder$t_try)),
      search = sort(as.integer(search)),
      min_leaf_rows = as.integer(min_leaf_rows),
      neighbours = neighbours
    ),
    class = "rpf_tune_control"
  )
}

# rpf()'s `tune`: NULL (no tuning) or an rpf_tune_control
tune_resolve <- function(tune, call) {
  if (inherits(tune, "rpf_tune_control")) {
    return(tune)
  }
  if (isFALSE(tune)) {
    return(NULL)
  }
  if (isTRUE(tune)) {
    return(tune_control())
  }
  if (identical(tune, "quick")) {
    return(tune_control(search = c(5L, 7L), min_leaf_rows = 1L, neighbours = FALSE))
  }
  cli::cli_abort(
    "{.arg tune} must be {.code FALSE}, {.code TRUE}, {.code \"quick\"} or a {.fn tune_control} object.",
    call = call
  )
}

# Which of the settings that `tune` takes over the user gave in the call of an rpf() method (its environment)
tune_given <- function(env) {
  args <- c("splits", "split_try", "t_try", "min_leaf_rows", "prune")
  vapply(args, function(a) !eval(call("missing", as.name(a)), env), logical(1))
}

# The checks of rpf(tune = ): random regression with at least 10 families (the pruning's out-of-bag choice), and none
# of the settings the search chooses given
tune_check_call <- function(control, given, splits, prune, mode, deterministic, ntrees, split_structure, call) {
  if (mode != "regression") {
    cli::cli_abort(
      c("{.arg tune} supports regression only.", "x" = "The outcome is a factor (classification)."),
      call = call
    )
  }
  if (deterministic) {
    cli::cli_abort(
      "{.arg tune} needs out-of-bag rows: it is not available with {.code deterministic = TRUE}.",
      call = call
    )
  }
  set <- intersect(c("split_try", "t_try", "min_leaf_rows"), names(given)[given])
  if (length(set) > 0) {
    cli::cli_abort(
      c(
        "{.arg tune} chooses {.arg split_try}, {.arg t_try} and {.arg min_leaf_rows}: don't give {.arg {set}}.",
        "i" = "Set the values to search with {.fn tune_control}."
      ),
      call = call
    )
  }
  if (given[["splits"]] && !is.null(splits)) {
    cli::cli_abort(
      "{.arg tune} chooses {.arg splits} out of bag ({.code splits = NULL}): don't give a number.",
      call = call
    )
  }
  if (given[["prune"]] && !prune) {
    cli::cli_abort(
      "{.arg tune} scores every forest by its pruned model and returns the pruned model: {.code prune = FALSE} is not possible.",
      call = call
    )
  }
  if (ntrees < prune_oob_min_ntrees) {
    cli::cli_abort(
      c(
        "{.arg tune} needs at least {prune_oob_min_ntrees} tree families, not {ntrees}.",
        "i" = "Every forest is pruned out of bag, which needs out-of-bag rows of several families."
      ),
      call = call
    )
  }
  if (split_structure != "leaves" && any(control$min_leaf_rows > 1)) {
    cli::cli_abort(
      c(
        "{.arg tune} searches {.arg min_leaf_rows} above 1, which needs {.code split_structure = \"leaves\"}.",
        "i" = "Use {.code tune = tune_control(min_leaf_rows = 1)} with {.code split_structure = \"{split_structure}\"}."
      ),
      call = call
    )
  }
  invisible(TRUE)
}

# The search of rpf(tune = ) on rpf_bridge()'s checked arguments (`params` with splits = NA: chosen out of bag).
#
# Every forest is grown from the random number state on entry, by the fit of rpf(splits = NULL) at its point of the
# ladder and min_leaf_rows, and pruned by prune_components()' defaults without the display (rpf_prune(), max_panels =
# Inf), which draws its seed from the stream after the fit as rpf(prune = TRUE) does. Its score is the corrected OOB
# RMSE of that pruned model, the score prune_components() chose it by. The chosen forest is then pruned again from the
# same state, now with the display at `max_panels`: so a one-point search returns what set.seed(s); rpf(splits = NULL,
# prune = TRUE, ...) returns, and leaves the stream where that call leaves it.
#
# Returns the rpf_pruned model of the chosen forest, with `tune`: `search` (tune_oob_table()), `steps` (each step's
# candidates and choice, as keys "<ladder row>_<min_leaf_rows>"), `settings` (the chosen split_try, t_try,
# min_leaf_rows, splits), `control` and `seconds`.
tune_run <- function(processed, outcomes, predictors, params, control, max_panels, export_forest, call) {
  ladder <- control$ladder
  t_all <- proc.time()[["elapsed"]]
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  rng0 <- get(".Random.seed", envir = globalenv())

  rows <- list() # one record per forest grown
  forests <- list() # the forests (rpf) that can still be chosen
  rng_fit <- list() # per forest: the random number state after its fit, from which it is pruned
  key <- function(i, m) paste0(i, "_", m)
  # one forest: the adaptive fit at ladder row i and min_leaf_rows m from the state on entry; the warning that growing
  # stopped at the cap is recorded per forest (and given for the chosen one), not repeated
  fit_at <- function(i, m) {
    assign(".Random.seed", rng0, envir = globalenv())
    p <- params
    p$split_try <- ladder$split_try[i]
    p$t_try <- ladder$t_try[i]
    p$min_leaf_rows <- as.integer(m)
    withCallingHandlers(
      rpf_fit(
        processed,
        outcomes,
        predictors,
        p,
        adaptive = TRUE,
        checkpoints = NULL,
        export_forest = FALSE,
        call = call
      ),
      warning = function(w) {
        if (grepl("Growing stopped at the cap", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  grow <- function(i, m, step) {
    k <- key(i, m)
    if (!is.null(rows[[k]])) {
      return(k)
    }
    t0 <- proc.time()[["elapsed"]]
    res <- fit_at(i, m)
    rng_fit[[k]] <<- get(".Random.seed", envir = globalenv())
    pr <- rpf_prune(res, processed, Inf, call = call)
    # the chosen candidate's row of the out-of-bag table: its corrected OOB RMSE is the forest's score
    j <- which(pr$oob$rule == pr$rule & pr$oob$tau == pr$tau)[1]
    path <- res$oob_path
    rows[[k]] <<- list(
      ladder = i,
      step = step,
      split_try = ladder$split_try[i],
      t_try = ladder$t_try[i],
      min_leaf_rows = as.integer(m),
      splits = as.integer(res$params$splits),
      stopped_at = as.integer(res$splits_oob$stopped_at),
      capped = isTRUE(res$splits_oob$capped),
      oob_rmse_forest = path$oob_rmse[match(res$params$splits, path$splits)],
      rule = pr$rule,
      tau = pr$tau,
      n_kept = length(pr$beta),
      oob_rmse = pr$oob$oob_rmse[j],
      seconds = proc.time()[["elapsed"]] - t0
    )
    forests[[k]] <<- res
    k
  }
  # the choice among the forests `ks`; simpler: their ranks (1 = simplest)
  choose <- function(ks, simpler) {
    ks[tune_oob_choose(vapply(ks, function(k) rows[[k]]$oob_rmse^2, numeric(1)), simpler)]
  }
  # one step: grow the forests (i, m) of the candidates not grown yet, in order, and choose among all of them. Only
  # the forests that can still be chosen are kept: the step's best so far and the previous step's choice
  steps <- list()
  incumbent <- NULL
  run_step <- function(name, step, cand_i, cand_m, simpler) {
    ks <- character(0)
    for (j in seq_along(cand_i)) {
      ks <- c(ks, grow(cand_i[j], cand_m[j], step))
      best <- choose(ks, rank(simpler[seq_along(ks)], ties.method = "first"))
      for (k in setdiff(names(forests), c(best, incumbent))) {
        forests[[k]] <<- NULL
      }
    }
    chosen <- choose(ks, simpler)
    steps[[name]] <<- list(candidates = ks, chosen = chosen)
    incumbent <<- chosen
    for (k in setdiff(names(forests), chosen)) {
      forests[[k]] <<- NULL
    }
    invisible(gc())
    rows[[chosen]]
  }

  # step 1: greediness at the first min_leaf_rows
  search <- control$search
  m1 <- control$min_leaf_rows[1]
  i1 <- run_step("step1", 1L, search, rep(m1, length(search)), rank_greedy(search))$ladder

  # step 2: the minimum leaf size at the chosen point (larger is simpler)
  m2 <- m1
  m_all <- control$min_leaf_rows
  if (length(m_all) > 1) {
    m2 <- run_step("step2", 2L, rep(i1, length(m_all)), m_all, rank(-m_all, ties.method = "first"))$min_leaf_rows
  }

  # step 3: the ladder neighbours of the chosen point at the chosen m; the choice among every forest grown at that m
  i3 <- i1
  if (control$neighbours) {
    nb <- intersect(c(i1 - 1L, i1 + 1L), seq_len(nrow(ladder)))
    grown_at_m <- vapply(rows, function(r) if (r$min_leaf_rows == m2) r$ladder else NA_integer_, integer(1))
    at_m <- sort(unique(c(i1, nb, grown_at_m[!is.na(grown_at_m)])))
    if (length(at_m) > 1) {
      i3 <- run_step("step3", 3L, at_m, rep(m2, length(at_m)), rank_greedy(at_m))$ladder
    }
  }
  final <- key(i3, m2)

  # the final model: the chosen forest, pruned from the state after its fit with the display (a forest dropped by
  # mistake would be grown again, the same forest)
  res <- forests[[final]] %||% fit_at(i3, m2)
  rm(forests)
  if (export_forest) {
    forest <- res$fit$get_model()
    class(forest) <- "rpf_forest"
    res$forest <- forest
  }
  assign(".Random.seed", rng_fit[[final]], envir = globalenv())
  pruned <- rpf_prune(res, processed, max_panels, call = call)

  fr <- rows[[final]]
  if (fr$capped) {
    cli::cli_warn(
      c(
        "The chosen forest stopped growing at the cap ({fr$stopped_at} splits) while its out-of-bag error was still falling.",
        "i" = "A fit with more {.arg splits} and {.code oob_path = TRUE} shows whether more splits help."
      ),
      call = call
    )
  }
  pruned$tune <- list(
    search = tune_oob_table(rows, steps, final),
    steps = steps,
    settings = fr[c("split_try", "t_try", "min_leaf_rows", "splits")],
    control = control,
    seconds = proc.time()[["elapsed"]] - t_all
  )
  pruned
}

# ranks of ladder rows by simplicity: the greediest (smallest row) is the simplest
rank_greedy <- function(idx) rank(idx, ties.method = "first")

# The choice among k forests by their corrected OOB MSEs `mse`: the lowest (ties, relative 1e-12: the simplest by the
# ranks `simpler`, 1 = simplest). Returns its index.
tune_oob_choose <- function(mse, simpler) {
  mse <- unname(mse)
  ties <- which(mse <= min(mse) * (1 + 1e-12) | mse == min(mse))
  ties[order(simpler[ties])][1]
}

# The search table: one row per forest grown, in the order grown, with each step's choice (NA: not a candidate of
# the step) and the final one
tune_oob_table <- function(rows, steps, final) {
  cols <- c(
    "step",
    "split_try",
    "t_try",
    "min_leaf_rows",
    "splits",
    "stopped_at",
    "capped",
    "oob_rmse_forest",
    "rule",
    "tau",
    "n_kept",
    "oob_rmse",
    "seconds"
  )
  tab <- do.call(rbind, lapply(rows, function(r) as.data.frame(r[cols], stringsAsFactors = FALSE)))
  keys <- names(rows)
  for (s in c("step1", "step2", "step3")) {
    st <- steps[[s]]
    tab[[paste0("chosen_", s)]] <- if (is.null(st)) NA else ifelse(keys %in% st$candidates, keys == st$chosen, NA)
  }
  tab$final <- keys == final
  rownames(tab) <- NULL
  tab
}
