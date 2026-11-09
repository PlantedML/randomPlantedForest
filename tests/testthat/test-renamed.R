# The arguments renamed after 0.5.0 keep their old names as deprecated aliases: the same fit, bit for bit, a
# deprecation warning once per session for each old name, and an error when both names are given.

# Forget the old names warned about, now and when the test ends, with lifecycle's verbosity at its default
local_fresh_deprecations <- function(frame = parent.frame()) {
  rlang::local_options(lifecycle_verbosity = NULL, .frame = frame)
  randomPlantedForest:::rpf_deprecation_reset()
  do.call(on.exit, list(quote(randomPlantedForest:::rpf_deprecation_reset()), add = TRUE), envir = frame)
}

renamed_fit <- function(..., seed = 1) {
  set.seed(seed)
  rpf(mpg ~ ., data = mtcars, max_interaction = 2, prune = FALSE, ...)
}

# the whole fit: every family's trees, leaf values and intervals, and the predictions
renamed_numbers <- function(fit) {
  list(fit$fit$get_model(), predict(fit, mtcars), fit$params)
}

test_that("each old argument name gives the fit of its new name, with a deprecation warning", {
  local_fresh_deprecations()
  new_args <- list(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3)
  old_names <- c(
    n_families = "ntrees",
    n_splits = "splits",
    n_try_split_points = "split_try",
    share_candidates = "t_try"
  )
  ref <- do.call(renamed_fit, new_args)
  for (new in names(new_args)) {
    args <- new_args
    names(args)[names(args) == new] <- old_names[[new]]
    expect_warning(
      fit <- do.call(renamed_fit, args),
      sprintf("Use `%s` instead of `%s`", new, old_names[[new]]),
      class = "rpf_deprecated"
    )
    expect_identical(renamed_numbers(fit), renamed_numbers(ref))
  }
})

test_that("all old names at once give one warning and the same fit", {
  local_fresh_deprecations()
  ref <- renamed_fit(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3)
  warnings <- list()
  fit <- withCallingHandlers(
    renamed_fit(ntrees = 3, splits = 20, split_try = 4, t_try = 0.3, split_structure = "leaves"),
    warning = function(w) {
      warnings[[length(warnings) + 1]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1)
  expect_s3_class(warnings[[1]], "lifecycle_warning_deprecated")
  msg <- conditionMessage(warnings[[1]])
  for (old in c("ntrees", "splits", "split_try", "t_try", "\"leaves\"")) {
    expect_match(msg, old, fixed = TRUE)
  }
  expect_identical(renamed_numbers(fit), renamed_numbers(ref))
})

test_that("the new names give no warning", {
  expect_no_warning(renamed_fit(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3))
})

test_that("each old split structure name gives the fit of its new name", {
  local_fresh_deprecations()
  old_new <- c(
    leaves = "leaf",
    hist = "leaf_binned",
    cur_trees_1 = "tree_then_leaf",
    cur_trees_2 = "tree_then_all_leaves",
    res_trees = "target_tree_then_leaf"
  )
  for (old in names(old_new)) {
    ref <- renamed_fit(n_families = 2, n_splits = 15, split_structure = old_new[[old]])
    expect_warning(
      fit <- renamed_fit(n_families = 2, n_splits = 15, split_structure = old),
      sprintf("split_structure = \"%s\"", old_new[[old]]),
      fixed = TRUE,
      class = "rpf_deprecated"
    )
    expect_identical(renamed_numbers(fit), renamed_numbers(ref))
    # the stored params keep the old code
    expect_identical(fit$params$split_structure, old)
  }
})

test_that("giving an argument under both names is an error", {
  both <- list(
    list(ntrees = 3, n_families = 3),
    list(splits = 20, n_splits = 20),
    list(split_try = 4, n_try_split_points = 4),
    list(t_try = 0.3, share_candidates = 0.3),
    # NULL counts as given, under either name
    list(splits = NULL, n_splits = 20),
    list(splits = 20, n_splits = NULL)
  )
  for (args in both) {
    expect_error(do.call(renamed_fit, args), "were given")
  }
})

test_that("an old name given as NULL is the new name's NULL", {
  local_fresh_deprecations()
  set.seed(1)
  ref <- rpf(mpg ~ ., data = mtcars, n_families = 10, n_splits = NULL, prune = FALSE)
  set.seed(1)
  expect_warning(
    fit <- rpf(mpg ~ ., data = mtcars, n_families = 10, splits = NULL, prune = FALSE),
    class = "rpf_deprecated"
  )
  expect_identical(fit$params$splits, ref$params$splits)
  expect_identical(fit$oob_path, ref$oob_path)
  expect_identical(renamed_numbers(fit), renamed_numbers(ref))
})

test_that("the old names work in every rpf() method, and other unknown arguments are still an error", {
  x <- mtcars[, c("cyl", "wt", "hp")]
  set.seed(1)
  ref <- rpf(x, mtcars$mpg, n_families = 2, n_splits = 10, prune = FALSE)
  local_fresh_deprecations()
  set.seed(1)
  expect_warning(fit <- rpf(x, mtcars$mpg, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, x), predict(ref, x))
  randomPlantedForest:::rpf_deprecation_reset()
  set.seed(1)
  expect_warning(fit <- rpf(as.matrix(x), mtcars$mpg, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, as.matrix(x)), predict(ref, x))
  skip_if_not_installed("recipes")
  rec <- recipes::recipe(mpg ~ cyl + wt + hp, data = mtcars)
  randomPlantedForest:::rpf_deprecation_reset()
  set.seed(1)
  expect_warning(fit <- rpf(rec, mtcars, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, mtcars), predict(ref, x))
  expect_error(rpf(x, mtcars$mpg, n_tree = 2), class = "rlib_error_dots_nonempty")
  expect_error(rpf(x, mtcars$mpg, ntrees = 2, foo = 1), class = "rlib_error_dots_nonempty")
})

test_that("tune: an old name counts as given", {
  expect_error(
    suppressWarnings(rpf(mpg ~ ., data = mtcars, n_families = 10, tune = TRUE, split_try = 3)),
    "n_try_split_points"
  )
})

test_that("oob_path_checkpoints() takes the old name", {
  local_fresh_deprecations()
  expect_warning(
    cp <- oob_path_checkpoints(splits = 100),
    "Use `n_splits` instead of `splits`",
    class = "rpf_deprecated"
  )
  expect_identical(cp, oob_path_checkpoints(n_splits = 100))
  expect_error(suppressWarnings(oob_path_checkpoints(100, splits = 100)), "were given")
  expect_error(oob_path_checkpoints(100, ntrees = 5), "Unused argument")
})

test_that("tune_control() takes a ladder with the old column names", {
  local_fresh_deprecations()
  old <- data.frame(split_try = c(3, 1), t_try = c(0.1, 0.025))
  new <- data.frame(n_try_split_points = c(3, 1), share_candidates = c(0.1, 0.025))
  expect_warning(tc <- tune_control(ladder = old, search = 1:2), class = "rpf_deprecated")
  expect_identical(tc, tune_control(ladder = new, search = 1:2))
  expect_error(tune_control(ladder = cbind(old, new), search = 1:2), "give one")
  # once per session
  expect_no_warning(tune_control(ladder = old, search = 1:2))
})

test_that("each old name warns once per session; other old names still warn", {
  local_fresh_deprecations()
  ref <- renamed_fit(n_families = 3, n_splits = 20)
  expect_warning(renamed_fit(ntrees = 3, n_splits = 20), "`ntrees`", class = "rpf_deprecated")
  # the same old name again: no warning, the same fit (a learner fitting in a loop warns once)
  for (i in 1:3) {
    expect_no_warning(fit <- renamed_fit(ntrees = 3, n_splits = 20))
    expect_identical(renamed_numbers(fit), renamed_numbers(ref))
  }
  # a call with a warned and a new old name warns about the new one only
  w <- expect_warning(renamed_fit(ntrees = 3, splits = 20), class = "rpf_deprecated")
  expect_match(conditionMessage(w), "`splits`", fixed = TRUE)
  expect_no_match(conditionMessage(w), "`ntrees`", fixed = TRUE)
  expect_match(conditionMessage(w), "once per session", fixed = TRUE)
  # an old split structure is its own old name; another function shares the record (splits is warned already)
  expect_warning(renamed_fit(n_families = 2, n_splits = 10, split_structure = "hist"), "leaf_binned")
  expect_no_warning(renamed_fit(n_families = 2, n_splits = 10, split_structure = "hist"))
  expect_no_warning(oob_path_checkpoints(splits = 100))
  # both names stay an error after the warning
  expect_error(renamed_fit(ntrees = 3, n_families = 3), "were given")
})

test_that("lifecycle_verbosity: \"warning\" warns every time, \"quiet\" never", {
  local_fresh_deprecations()
  rlang::local_options(lifecycle_verbosity = "warning")
  for (i in 1:2) {
    expect_warning(renamed_fit(ntrees = 2, n_splits = 10), class = "rpf_deprecated")
  }
  rlang::local_options(lifecycle_verbosity = "quiet")
  randomPlantedForest:::rpf_deprecation_reset()
  expect_no_warning(renamed_fit(ntrees = 2, n_splits = 10))
})

test_that("print() and the result tables use the new names", {
  set.seed(1)
  fit <- rpf(mpg ~ ., data = mtcars, n_families = 10, prune = FALSE)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "n_try_split_points: 3", fixed = TRUE)
  expect_match(out, "share_candidates: 0.1", fixed = TRUE)
  expect_match(out, "split_structure: leaf\n", fixed = TRUE)
  expect_match(out, "`n_splits` chosen out of bag", fixed = TRUE)
  expect_identical(names(fit$oob_path)[1], "n_splits")
  expect_identical(names(fit$splits_oob), c("n_splits", "stopped_at", "capped"))
  expect_named(greediness_ladder(), c("n_try_split_points", "share_candidates"))
})
