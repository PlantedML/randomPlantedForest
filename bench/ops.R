# What gets measured: simulated data and the operations on it. Op functions run
# in the measurement child, so they call the package by namespace and use only
# their arguments and the helpers sourced there (refs.R, measure.R, ops.R).

bench_make_data <- function(task = c("regr", "binary", "multiclass"), n, p, seed = 1L) {
  task <- match.arg(task)
  stopifnot(p >= 3)
  set.seed(seed)
  m <- n + 1000L
  X <- as.data.frame(matrix(stats::runif(m * p), m, p))
  X$f <- factor(sample(letters[1:5], m, replace = TRUE))
  # main effects, a two-way and a three-way interaction, and the factor
  signal <- 2 * X$V1 + sin(3 * X$V2) + 2 * X$V1 * X$V3 + 4 * X$V1 * X$V2 * X$V3 + as.integer(X$f) / 5
  noisy <- signal + stats::rnorm(m, sd = 0.5)
  X$y <- switch(
    task,
    regr = noisy,
    binary = factor(ifelse(noisy > stats::median(noisy), "b", "a")),
    multiclass = cut(noisy, stats::quantile(noisy, 0:3 / 3), labels = c("a", "b", "c"), include.lowest = TRUE)
  )
  list(train = X[seq_len(n), ], test = X[n + seq_len(1000L), ])
}

# Names and values of rpf's arguments before the renaming in 0.5.0.9000.
BENCH_OLD_ARGS <- c(
  n_families = "ntrees",
  n_splits = "splits",
  n_try_split_points = "split_try",
  share_candidates = "t_try"
)
BENCH_OLD_STRUCTURES <- c(
  leaf = "leaves",
  leaf_binned = "hist",
  tree_then_leaf = "cur_trees_1",
  tree_then_all_leaves = "cur_trees_2",
  target_tree_then_leaf = "res_trees"
)
# `available`: the names of rpf()'s arguments in the installed ref. Arguments
# it does not have are dropped (prune before it existed).
bench_rpf_args <- function(args, available) {
  if (!"n_families" %in% available) {
    hit <- names(args) %in% names(BENCH_OLD_ARGS)
    names(args)[hit] <- BENCH_OLD_ARGS[names(args)[hit]]
    if (!is.null(args$split_structure)) {
      args$split_structure <- unname(BENCH_OLD_STRUCTURES[args$split_structure])
    }
  }
  args[names(args) %in% available]
}

# The cell's fit settings, translated for the installed version.
.bench_fit_args <- function(data, cell, ...) {
  args <- c(
    list(
      x = data$train[setdiff(names(data$train), "y")],
      y = data$train$y,
      max_interaction = cell$max_interaction,
      n_families = cell$n_families,
      nthreads = cell$nthreads
    ),
    list(...)
  )
  if (!is.na(cell$loss)) {
    args$loss <- cell$loss
  }
  if (!is.na(cell$split_structure)) {
    args$split_structure <- cell$split_structure
  }
  bench_rpf_args(args, .bench_rpf_formals())
}

.bench_fit <- function(data, cell, ...) {
  do.call(randomPlantedForest::rpf, .bench_fit_args(data, cell, ...))
}

# The ref's own estimate of an op's additional peak memory in MB, if its
# namespace exports `fn`: the largest of the steps it returns. The names are
# placeholders until rpf has estimators (one per function, as in the futureverse
# resources proposal). A failing estimator gives no bound, not a lost cell.
.bench_estimate <- function(fn, ...) {
  ns <- asNamespace("randomPlantedForest")
  if (!fn %in% getNamespaceExports(ns)) {
    return(NA_real_)
  }
  tryCatch(max(unlist(get(fn, envir = ns)(...))) / 1024^2, error = function(e) NA_real_)
}

# Rows whose predictions stand for a fitted model in the output hash.
.bench_probe <- function(data) utils::head(data$test, 100L)

bench_ops <- function() {
  fit_op <- function(requires, ...) {
    extra <- list(...)
    list(
      requires = requires,
      setup = function(data, cell) data,
      bound = function(state, cell) {
        do.call(.bench_estimate, c(list("rpf_resources"), do.call(.bench_fit_args, c(list(state, cell), extra))))
      },
      run = function(state, cell) {
        fit <- do.call(.bench_fit, c(list(state, cell), extra))
        stats::predict(fit, .bench_probe(state))
      }
    )
  }
  fitted <- function(data, cell) {
    c(data, list(fit = .bench_fit(data, cell, n_splits = cell$n_splits, prune = FALSE)))
  }
  list(
    fit_fixed = list(
      setup = function(data, cell) data,
      bound = function(state, cell) {
        do.call(
          .bench_estimate,
          c(list("rpf_resources"), .bench_fit_args(state, cell, n_splits = cell$n_splits, prune = FALSE))
        )
      },
      run = function(state, cell) {
        fit <- .bench_fit(state, cell, n_splits = cell$n_splits, prune = FALSE)
        stats::predict(fit, .bench_probe(state))
      }
    ),
    # oob_path came with n_splits chosen out of bag, prune with the pruned default
    fit_oob = fit_op("oob_path", n_splits = NULL, prune = FALSE),
    fit_default = fit_op("prune"),
    # purify() changes the forest in place, so a second call would time a no-op
    purify = list(
      single = TRUE,
      bound = function(state, cell) .bench_estimate("purify_resources", state$fit),
      setup = fitted,
      run = function(state, cell) {
        p <- randomPlantedForest::purify(state$fit)
        randomPlantedForest::predict_components(p, .bench_probe(state))
      }
    ),
    predict = list(
      bound = function(state, cell) .bench_estimate("predict_resources", state$fit, new_data = state$test),
      setup = fitted,
      run = function(state, cell) stats::predict(state$fit, state$test)
    ),
    predict_components = list(
      bound = function(state, cell) {
        .bench_estimate("predict_components_resources", state$fit, new_data = state$test)
      },
      setup = function(data, cell) {
        s <- fitted(data, cell)
        s$fit <- randomPlantedForest::purify(s$fit)
        s
      },
      run = function(state, cell) randomPlantedForest::predict_components(state$fit, state$test)
    )
  )
}

bench_ref_supports <- function(op, args) {
  all(bench_ops()[[op]]$requires %in% args)
}
