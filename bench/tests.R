# Tests for the benchmark suite itself, run from the repository root through the
# runner, so they use the same rig environment as real runs:
#
#   rig run bench/rpf-bench test
#
# Not part of the package's testthat suite: bench/ is in .Rbuildignore.

library(testthat)

# test_file() runs from bench/, the suite expects the repository root.
if (basename(getwd()) == "bench") {
  setwd("..")
}

for (f in setdiff(Sys.glob("bench/*.R"), c("bench/cli.R", "bench/tests.R"))) {
  source(f)
}

test_that("bench runs in its own library", {
  expect_match(.libPaths()[1], "rig")
})

## ---- fixtures ------------------------------------------------------------------

git <- function(repo, ...) {
  out <- system2("git", c("-C", repo, ...), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"))
  }
  out
}

# A throwaway git repo holding an R-only package named randomPlantedForest, so
# ref handling is tested without compiling the real package.
fake_rpf_source <- function(version = "0.0.1", body = "list(...)", extra = character(), args = character()) {
  c(
    paste0("rpf <- function(", paste(c(if (length(args)) paste0(args, " = NULL"), "..."), collapse = ", "), ") ", body),
    extra
  )
}

write_fake_rpf <- function(
  repo,
  version = "0.0.1",
  body = "list(...)",
  extra = character(),
  exports = "rpf",
  args = character()
) {
  dir.create(file.path(repo, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "Package: randomPlantedForest",
      "Title: Fake",
      paste0("Version: ", version),
      "Description: Fake package for bench tests.",
      "License: MIT",
      "Authors@R: person('A', 'B', role = c('aut', 'cre'), email = 'a@b.c')"
    ),
    file.path(repo, "DESCRIPTION")
  )
  writeLines(paste0("export(", exports, ")"), file.path(repo, "NAMESPACE"))
  writeLines(fake_rpf_source(version, body, extra, args), file.path(repo, "R", "rpf.R"))
}

fake_repo <- function() {
  repo <- tempfile("fakerpf")
  dir.create(repo)
  git(repo, "init", "-q", "-b", "main")
  git(repo, "config", "user.email", "bench@example.com")
  git(repo, "config", "user.name", "bench")
  write_fake_rpf(repo)
  writeLines("readme", file.path(repo, "README.md"))
  git(repo, "add", "-A")
  git(repo, "commit", "-q", "-m", "init")
  repo
}

commit_all <- function(repo, msg) {
  git(repo, "add", "-A")
  git(repo, "commit", "-q", "--allow-empty", "-m", msg)
}

## ---- refs ------------------------------------------------------------------------

test_that("refs parse into worktree and git refs", {
  expect_equal(bench_parse_ref("WORKTREE")$type, "worktree")
  expect_equal(bench_parse_ref("dev")$type, "git")
  expect_equal(bench_parse_ref("dev")$rev, "dev")
})

test_that("a ref only on origin resolves to the remote branch", {
  upstream <- fake_repo()
  git(upstream, "branch", "feature")
  clone <- tempfile("clone")
  system2("git", c("clone", "-q", upstream, clone))
  expect_equal(bench_ref_id("feature", clone), bench_ref_id("origin/feature", clone))
  expect_equal(
    bench_install_ref("feature", tempfile("lib"), clone)$commit,
    git(clone, "rev-parse", "origin/feature")
  )
  expect_error(bench_ref_id("no-such-ref", clone), "no-such-ref")
})

test_that("the ref id covers package content only", {
  repo <- fake_repo()
  id1 <- bench_ref_id("HEAD", repo)
  writeLines("changed readme", file.path(repo, "README.md"))
  commit_all(repo, "readme")
  expect_equal(bench_ref_id("HEAD", repo), id1)
  write_fake_rpf(repo, body = "list(..., changed = TRUE)")
  commit_all(repo, "code")
  expect_false(identical(bench_ref_id("HEAD", repo), id1))
})

test_that("the working tree has the same id as HEAD when nothing changed", {
  repo <- fake_repo()
  expect_equal(bench_ref_id("WORKTREE", repo), bench_ref_id("HEAD", repo))
  write_fake_rpf(repo, body = "list(..., uncommitted = TRUE)")
  expect_false(identical(bench_ref_id("WORKTREE", repo), bench_ref_id("HEAD", repo)))
})

test_that("refs with identical package content are measured once", {
  repo <- fake_repo()
  writeLines("changed readme", file.path(repo, "README.md"))
  commit_all(repo, "readme")
  expect_message(refs <- bench_resolve_refs(c("HEAD", "HEAD~1"), repo), "HEAD~1")
  expect_equal(refs, "HEAD")
})

test_that("a ref is installed into its own library and reused", {
  repo <- fake_repo()
  root <- tempfile("lib")
  ref <- bench_install_ref("HEAD", root, repo)
  expect_true(file.exists(file.path(ref$lib, "randomPlantedForest", "DESCRIPTION")))
  expect_equal(ref$version, "0.0.1")
  expect_equal(ref$commit, git(repo, "rev-parse", "HEAD"))
  expect_equal(ref$args, "...")
  log <- paste0(ref$lib, ".log")
  stamp <- file.mtime(log)
  Sys.sleep(1.1)
  again <- bench_install_ref("HEAD", root, repo)
  expect_equal(again$lib, ref$lib)
  expect_equal(file.mtime(log), stamp)
})

test_that("the working tree installs with its uncommitted changes", {
  repo <- fake_repo()
  write_fake_rpf(repo, version = "0.0.2")
  ref <- bench_install_ref("WORKTREE", tempfile("lib"), repo)
  expect_equal(ref$version, "0.0.2")
})

test_that("a failed install names the ref", {
  repo <- fake_repo()
  writeLines("rpf <- function( {", file.path(repo, "R", "rpf.R"))
  commit_all(repo, "broken")
  expect_error(bench_install_ref("HEAD", tempfile("lib"), repo), "HEAD")
})

test_that("a process without the ref library cannot pass the library check", {
  repo <- fake_repo()
  ref <- bench_install_ref("HEAD", tempfile("lib"), repo)
  ok <- callr::r(
    function(lib, helpers) {
      source(helpers)
      .libPaths(c(lib, .libPaths()))
      bench_assert_lib(lib)
      TRUE
    },
    args = list(lib = ref$lib, helpers = normalizePath("bench/refs.R"))
  )
  expect_true(ok)
  expect_error(
    callr::r(
      function(lib, helpers) {
        source(helpers)
        bench_assert_lib(lib)
      },
      args = list(lib = ref$lib, helpers = normalizePath("bench/refs.R"))
    ),
    "randomPlantedForest"
  )
})

test_that("a dirty tree is refused unless the working tree is measured", {
  repo <- fake_repo()
  write_fake_rpf(repo, body = "list(..., dirty = TRUE)")
  expect_error(bench_assert_clean_tree(c("HEAD", "main"), repo), "uncommitted")
  expect_silent(bench_assert_clean_tree(c("WORKTREE", "main"), repo))
})

## ---- measurement -------------------------------------------------------------------

fake_ref <- function(body = "list(...)") {
  repo <- fake_repo()
  write_fake_rpf(repo, body = body)
  commit_all(repo, "fake")
  bench_install_ref("HEAD", tempfile("lib"), repo)
}

alloc_op <- list(
  min_version = "0.0.0",
  setup = function(data, cell) NULL,
  run = function(state, cell) {
    x <- numeric(25e6) # 200 MB
    x[1] <- 1
    sum(x)
  }
)

measure_cell <- list(iters = 2L)

measure_data <- function() {
  path <- tempfile(fileext = ".rds")
  saveRDS(list(), path)
  path
}

test_that("RSS measurement sees a 200 MB allocation", {
  skip_if_not(ps::ps_is_supported())
  ref <- fake_ref()
  m <- bench_measure(ref$lib, alloc_op, measure_cell, measure_data(), metric = "rss")
  expect_null(m$error)
  expect_length(m$time, 2L)
  expect_gt(m$peak_mb - m$floor_mb, 150)
  expect_lt(m$peak_mb - m$floor_mb, 400)
  expect_equal(m$digest, bench_digest(1))
})

test_that("cgroup measurement sees a 200 MB allocation", {
  skip_if(is.null(.bench_cgroup_dir()))
  ref <- fake_ref()
  m <- bench_measure(ref$lib, alloc_op, measure_cell, measure_data(), metric = "cgroup")
  expect_null(m$error)
  expect_gt(m$peak_mb - m$floor_mb, 150)
  expect_lt(m$peak_mb - m$floor_mb, 400)
})

test_that("without memory accounting only time is measured", {
  ref <- fake_ref()
  m <- bench_measure(ref$lib, alloc_op, measure_cell, measure_data(), metric = "none")
  expect_true(is.na(m$peak_mb))
  expect_true(is.na(m$floor_mb))
  expect_length(m$time, 2L)
})

test_that("an error in the measured op comes back instead of crashing", {
  ref <- fake_ref()
  boom <- list(min_version = "0.0.0", setup = function(data, cell) NULL, run = function(state, cell) stop("boom"))
  m <- bench_measure(ref$lib, boom, measure_cell, measure_data(), metric = "none")
  expect_match(m$error, "boom")
})

test_that("the measurement child loads the ref's package with pinned threads", {
  ref <- fake_ref(body = "Sys.getenv('OMP_NUM_THREADS')")
  op <- list(
    min_version = "0.0.0",
    setup = function(data, cell) NULL,
    run = function(state, cell) randomPlantedForest::rpf()
  )
  m <- bench_measure(ref$lib, op, measure_cell, measure_data(), metric = "none")
  expect_null(m$error)
  expect_equal(m$digest, bench_digest("1"))
})

test_that("digests compare values", {
  expect_equal(bench_digest(c(1, 2, 3)), bench_digest(c(1, 2, 3)))
  expect_false(identical(bench_digest(c(1, 2, 3)), bench_digest(c(1, 2, 4))))
  expect_equal(bench_digest(data.frame(a = 1:3)), bench_digest(data.frame(a = 1:3)))
})

test_that("the metric follows the platform and the scope", {
  expect_equal(.bench_choose_metric("Darwin", cgroup = NULL, scoped = FALSE, slurm = FALSE, ps_ok = TRUE), "rss")
  expect_equal(.bench_choose_metric("Linux", cgroup = "/x", scoped = TRUE, slurm = FALSE, ps_ok = TRUE), "cgroup")
  expect_equal(.bench_choose_metric("Linux", cgroup = "/x", scoped = FALSE, slurm = TRUE, ps_ok = TRUE), "cgroup")
  expect_equal(.bench_choose_metric("Linux", cgroup = "/x", scoped = FALSE, slurm = FALSE, ps_ok = TRUE), "none")
  expect_equal(.bench_choose_metric("Darwin", cgroup = NULL, scoped = FALSE, slurm = FALSE, ps_ok = FALSE), "none")
})

test_that("a local Linux run outside a scope asks for one", {
  expect_true(.bench_needs_scope("Linux", scoped = FALSE, slurm = FALSE, has_systemd_run = TRUE))
  expect_false(.bench_needs_scope("Linux", scoped = TRUE, slurm = FALSE, has_systemd_run = TRUE))
  expect_false(.bench_needs_scope("Linux", scoped = FALSE, slurm = TRUE, has_systemd_run = TRUE))
  expect_false(.bench_needs_scope("Darwin", scoped = FALSE, slurm = FALSE, has_systemd_run = FALSE))
})

## ---- data, ops and cells -------------------------------------------------------------

test_that("simulated data has the requested shape", {
  d <- bench_make_data("multiclass", 200, 4)
  expect_s3_class(d$train$y, "factor")
  expect_equal(nlevels(d$train$y), 3L)
  expect_equal(nrow(d$train), 200L)
  expect_equal(nrow(d$test), 1000L)
  expect_equal(ncol(d$train), 4L + 2L) # p numeric, one factor, y
  expect_equal(nlevels(bench_make_data("binary", 100, 3)$train$y), 2L)
  expect_type(bench_make_data("regr", 100, 3)$train$y, "double")
  expect_equal(bench_make_data("regr", 50, 3), bench_make_data("regr", 50, 3))
})

old_api <- c(
  "x",
  "y",
  "max_interaction",
  "ntrees",
  "splits",
  "split_structure",
  "split_try",
  "t_try",
  "nthreads",
  "loss"
)
mid_api <- c(old_api, "oob_path", "prune") # before the rename, version 0.5.0.9000 already
new_api <- c(
  "x",
  "y",
  "max_interaction",
  "n_families",
  "n_splits",
  "split_structure",
  "nthreads",
  "loss",
  "oob_path",
  "prune"
)

test_that("arguments are translated by what rpf() accepts, not by version", {
  args <- list(n_families = 5, n_splits = 10, split_structure = "leaf", prune = FALSE)
  expect_equal(bench_rpf_args(args, old_api), list(ntrees = 5, splits = 10, split_structure = "leaves"))
  expect_equal(bench_rpf_args(args, mid_api), list(ntrees = 5, splits = 10, split_structure = "leaves", prune = FALSE))
  expect_equal(bench_rpf_args(args, new_api), args)
  expect_equal(bench_rpf_args(list(split_structure = "target_tree_then_leaf"), old_api)$split_structure, "res_trees")
})

test_that("ops needing newer features are skipped where rpf() lacks them", {
  expect_false(bench_ref_supports("fit_oob", old_api))
  expect_false(bench_ref_supports("fit_default", old_api))
  expect_true(bench_ref_supports("fit_fixed", old_api))
  expect_true(bench_ref_supports("fit_oob", mid_api))
  expect_true(bench_ref_supports("fit_default", new_api))
})

test_that("every cell names a known op", {
  for (tier in c("quick", "full")) {
    cells <- bench_cells(tier)
    expect_true(all(cells$op %in% names(bench_ops())), label = tier)
    expect_true(all(c("op", "task", "loss", "n", "nthreads", "iters", "mem_reps", "tier") %in% names(cells)))
    expect_true(all(is.na(cells$loss[cells$task == "regr"])))
  }
  quick <- bench_cells("quick")
  expect_equal(nrow(bench_select_cells(quick, ops = "predict")), sum(quick$op == "predict"))
  expect_true(all(bench_select_cells(quick, ops = "predict")$op == "predict"))
  expect_equal(nrow(bench_select_cells(quick, max_cells = 2)), 2L)
  expect_equal(unique(bench_select_cells(quick, sizes = 1000)$n), 1000)
})

test_that("prediction ops scale over the number of predicted rows", {
  expect_true("n_new" %in% BENCH_CELL_COLS)
  for (tier in c("quick", "full")) {
    cells <- bench_cells(tier)
    expect_true(all(cells$n_new[!cells$op %in% c("predict", "predict_components")] == 1000L), label = tier)
    expect_gt(length(unique(cells$n_new[cells$op == "predict"])), 1L)
  }
  full <- bench_cells("full")
  expect_setequal(unique(full$n_new[full$op == "predict_components"]), c(1000L, 10000L, 100000L))
})

test_that("predicted rows are drawn apart from the training rows", {
  a <- bench_make_data("regr", 50, 3, n_new = 10L)
  b <- bench_make_data("regr", 50, 3, n_new = 20L)
  expect_equal(a$train, b$train)
  expect_equal(nrow(a$test), 10L)
  expect_equal(nrow(b$test), 20L)
})

test_that("logit cells are sized to fit the default walltime", {
  full <- bench_cells("full")
  logit <- !is.na(full$loss) & full$loss == "logit"
  expect_true(any(logit))
  expect_true(all(full$mem_reps[logit] == 3L))
  expect_true(all(full$iters[logit & full$op %in% c("fit_fixed", "fit_oob")] == 1L))
  expect_true(all(full$iters[!logit & full$op == "fit_fixed"] == 5L))
})

test_that("choices are matched exactly, never partially", {
  expect_error(bench_cells("ful"), "Unknown tier 'ful'. Use one of quick, full.", fixed = TRUE)
  expect_error(bench_make_registry(tempfile("reg"), "slrum"), "Unknown backend 'slrum'", fixed = TRUE)
})

test_that("a filter value that selects nothing is an error naming the allowed ones", {
  quick <- bench_cells("quick")
  expect_error(
    bench_select_cells(quick, ops = c("predict", "fit_fixd")),
    "Unknown op 'fit_fixd'. Use one or more of",
    fixed = TRUE
  )
  expect_error(bench_select_cells(quick, tasks = "regression"), "Unknown task 'regression'", fixed = TRUE)
  expect_error(bench_select_cells(quick, sizes = "5k"), "Unknown size '5k'", fixed = TRUE)
  expect_error(bench_select_cells(quick, sizes = "123"), "Unknown size '123'", fixed = TRUE)
  expect_error(bench_select_cells(quick, max_cells = 0L), "--max-cells", fixed = TRUE)
})

test_that("the baseline must be one of the refs and survives deduplication", {
  expect_error(bench_order_refs(c("a", "b"), "dev"), "Baseline 'dev' is not among the refs (a, b)", fixed = TRUE)
  expect_equal(bench_order_refs(c("WORKTREE", "dev"), "dev"), c("dev", "WORKTREE"))
})

## ---- one cell across refs ---------------------------------------------------------------

two_refs <- function() {
  repo <- fake_repo()
  write_fake_rpf(repo, version = "0.0.1", body = "1")
  commit_all(repo, "a")
  root <- tempfile("lib")
  a <- bench_install_ref("HEAD", root, repo)
  write_fake_rpf(
    repo,
    version = "0.0.2",
    body = "2",
    extra = "rpf_resources <- function(...) {\n  list(...)\n  c(grow = 100 * 1024^2, prune = 50 * 1024^2)\n}",
    exports = "rpf, rpf_resources",
    args = "new_arg"
  )
  commit_all(repo, "b")
  b <- bench_install_ref("HEAD", root, repo)
  list(a = a, b = b)
}

fake_ops <- list(
  any = list(
    setup = function(data, cell) NULL,
    run = function(state, cell) randomPlantedForest::rpf(),
    bound = function(state, cell) .bench_estimate("rpf_resources")
  ),
  new = list(
    requires = "new_arg",
    setup = function(data, cell) NULL,
    run = function(state, cell) randomPlantedForest::rpf()
  )
)

fake_cell <- function(op) {
  cell <- bench_cells("quick")[1, ]
  cell$op <- op
  cell$n <- 100L
  cell$iters <- 2L
  cell$mem_reps <- 2L
  cell
}

test_that("a cell gives one row per ref with the schema's columns", {
  refs <- two_refs()
  rows <- bench_run_cell(fake_cell("any"), refs, tempfile("data"), ops = fake_ops, metric = "none")
  expect_named(rows, bench_schema())
  expect_equal(nrow(rows), 2L)
  expect_equal(rows$ref, c(refs$a$spec, refs$b$spec))
  expect_equal(rows$rpf_version, c("0.0.1", "0.0.2"))
  expect_false(identical(rows$digest[1], rows$digest[2]))
  expect_true(all(rows$time_median >= 0))
  expect_true(all(is.na(rows$skipped)))
})

test_that("a cell keeps every timed call and every round's memory", {
  rows <- bench_run_cell(fake_cell("any"), two_refs(), tempfile("data"), ops = fake_ops, metric = "none")
  # 2 rounds x 2 timed calls, 2 rounds
  expect_equal(lengths(strsplit(rows$times, ";")), c(4L, 4L))
  expect_equal(lengths(strsplit(rows$peak_deltas_mb, ";")), c(2L, 2L))
})

test_that("a ref without the op's arguments is skipped with a reason", {
  refs <- two_refs()
  rows <- bench_run_cell(fake_cell("new"), refs, tempfile("data"), ops = fake_ops, metric = "none")
  expect_equal(rows$skipped[1], "needs rpf() argument new_arg")
  expect_true(is.na(rows$time_median[1]))
  expect_true(is.na(rows$skipped[2]))
  expect_false(is.na(rows$time_median[2]))
})

test_that("the memory bound comes from refs that provide an estimator", {
  refs <- two_refs()
  rows <- bench_run_cell(fake_cell("any"), refs, tempfile("data"), ops = fake_ops, metric = "none")
  expect_true(is.na(rows$bound_mb[1]))
  expect_equal(rows$bound_mb[2], 100)
})


## ---- collation and report ------------------------------------------------------------------

fake_rows <- function(
  refs = c("dev", "pr"),
  time = c(1, 1.2),
  mem = c(100, 100),
  digest = c("x", "x"),
  bound = c(NA, NA),
  metric = "cgroup",
  skipped = c(NA, NA)
) {
  cell <- fake_cell("fit_fixed")
  rows <- lapply(seq_along(refs), function(i) {
    row <- as.list(setNames(rep(NA, length(bench_schema())), bench_schema()))
    row[BENCH_CELL_COLS] <- as.list(cell[BENCH_CELL_COLS])
    row$ref <- refs[i]
    row$rpf_version <- "0.5.0.9000"
    row$metric <- metric
    row$time_median <- time[i]
    row$peak_delta_mb <- mem[i]
    row$digest <- digest[i]
    row$bound_mb <- bound[i]
    row$skipped <- skipped[i]
    row$tier <- "quick"
    row$host <- "test"
    row$timestamp <- "2026-11-13T12:00:00+0100"
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

test_that("deltas get verdicts from the thresholds", {
  d <- bench_deltas(fake_rows(time = c(1, 1.2)), baseline = "dev")
  expect_equal(nrow(d), 1L)
  expect_equal(d$ref, "pr")
  expect_equal(d$time_pct, 20)
  expect_equal(d$time_verdict, "slower")
  expect_equal(bench_deltas(fake_rows(time = c(1, 1.05)))$time_verdict, "inconclusive")
  expect_equal(bench_deltas(fake_rows(time = c(1, 0.5)))$time_verdict, "faster")
  expect_equal(bench_deltas(fake_rows(mem = c(100, 110)))$mem_verdict, "more")
  expect_equal(bench_deltas(fake_rows(mem = c(100, 101)))$mem_verdict, "same")
  expect_false(bench_deltas(fake_rows())$output_differs)
  expect_true(bench_deltas(fake_rows(digest = c("x", "y")))$output_differs)
})

test_that("the report lists differing outputs and broken bounds", {
  rep <- bench_report(fake_rows(digest = c("x", "y"), mem = c(100, 120), bound = c(NA, 100)), baseline = "dev")
  expect_match(rep, "OUTPUT DIFFERS", all = FALSE)
  expect_match(rep, "exceeds", all = FALSE)
  expect_match(rep, "pr", all = FALSE)
})

test_that("the report lists skipped cells with their reason", {
  rep <- bench_report(fake_rows(skipped = c(NA, "needs rpf >= 0.5.0.9000"), time = c(1, NA), mem = c(100, NA)))
  expect_match(rep, "needs rpf >= 0.5.0.9000", all = FALSE, fixed = TRUE)
})

test_that("a single ref gives a report without deltas", {
  rows <- fake_rows(refs = "dev", time = 1, mem = 100, digest = "x", bound = NA, skipped = NA)
  expect_equal(nrow(bench_deltas(rows, baseline = "dev")), 0L)
  rep <- bench_report(rows, baseline = "dev")
  expect_match(rep, "Only one distinct ref was measured, so there is nothing to compare.", all = FALSE, fixed = TRUE)
})

test_that("the report says when memory was not measured", {
  rep <- bench_report(fake_rows(metric = "none", mem = c(NA, NA)))
  expect_match(rep, "Memory was not measured (no cgroup v2 or RSS).", all = FALSE, fixed = TRUE)
})

test_that("a baseline that was not measured is an error, not replaced", {
  expect_error(
    bench_deltas(fake_rows(refs = c("a", "b")), baseline = "dev"),
    "Baseline 'dev' was not measured. Use one of a, b.",
    fixed = TRUE
  )
})

test_that("results and history are written as plain tables", {
  dir <- tempfile("results")
  paths <- bench_write_results(fake_rows(), dir)
  expect_true(all(file.exists(paths)))
  expect_true(any(grepl("report\\.md$", paths)))
  hist <- tempfile(fileext = ".csv")
  bench_history_append(fake_rows(refs = c("v0.5.0", "dev")), anchors = "v0.5.0", path = hist)
  newer <- fake_rows(refs = c("v0.5.0", "dev"), time = c(2, 2))
  newer$timestamp <- "2026-11-14T12:00:00+0100"
  bench_history_append(newer, anchors = "v0.5.0", path = hist)
  h <- utils::read.csv(hist)
  expect_equal(nrow(h), 1L)
  expect_equal(h$ref, "v0.5.0")
  expect_equal(h$time_median, 2)
})

## ---- registry ---------------------------------------------------------------------------------

test_that("a local registry runs every cell against every ref", {
  refs <- two_refs()
  cells <- rbind(fake_cell("any"), fake_cell("new"))
  reg <- bench_make_registry(tempfile("reg"), "local")
  bench_submit(reg, cells, refs, ops = fake_ops, metric = "none")
  rows <- bench_collect(reg)
  expect_named(rows, bench_schema())
  expect_equal(nrow(rows), 4L)
  expect_equal(sum(!is.na(rows$skipped)), 1L)
})

test_that("cluster resources grow with threads and data size", {
  cells <- bench_cells("full")
  res <- bench_cell_resources(cells)
  expect_equal(nrow(res), nrow(cells))
  expect_true(all(res$ncpus == pmax(cells$nthreads, 1L) + 1L))
  expect_true(all(res$memory > 0))
  big <- which.max(cells$n * cells$n_families)
  small <- which.min(cells$n * cells$n_families)
  expect_gt(res$memory[big], res$memory[small])
  expect_true(all(res$walltime == 28800L))
})

test_that("skipped refs are marked in the delta table, not reported as same output", {
  rows <- fake_rows(skipped = c(NA, "needs rpf >= 0.5.0.9000"), time = c(1, NA), mem = c(100, NA), digest = c("x", NA))
  d <- bench_deltas(rows)
  expect_equal(d$skipped, "needs rpf >= 0.5.0.9000")
  rep <- bench_report(rows)
  # not a row of n/a in the table, but one line per ref, op and reason
  expect_false(any(grepl("| pr |", rep, fixed = TRUE)))
  expect_match(rep, "- pr cannot run fit_fixed (1 cell): needs rpf >= 0.5.0.9000", all = FALSE, fixed = TRUE)
})

test_that("the deltas table shows the values behind each percentage", {
  rep <- bench_report(fake_rows(time = c(1, 1.2), mem = c(100, 120)))
  expect_match(rep, "1 s → 1.2 s (+20.0%)", all = FALSE, fixed = TRUE)
  expect_match(rep, "100 MB → 120 MB (+20.0%)", all = FALSE, fixed = TRUE)
})

test_that("memory below resolution on both sides gets no verdict", {
  d <- bench_deltas(fake_rows(mem = c(0.03, 0)))
  expect_true(d$mem_below)
  expect_true(is.na(d$mem_verdict))
})

test_that("output that differs everywhere is one sentence, not a list", {
  rep <- bench_report(fake_rows(digest = c("x", "y")))
  expect_match(rep, "Every compared cell (1) gives output that differs from dev's.", all = FALSE, fixed = TRUE)
  expect_false(any(grepl(": pr vs dev", rep, fixed = TRUE)))
})

## ---- review fixes: measurement ------------------------------------------------------------

test_that("C1: random draws are seeded, so identical refs give identical output", {
  ref <- fake_ref()
  op <- list(min_version = "0.0.0", setup = function(data, cell) NULL, run = function(state, cell) stats::runif(3))
  a <- bench_measure(ref$lib, op, measure_cell, measure_data(), metric = "none")
  b <- bench_measure(ref$lib, op, measure_cell, measure_data(), metric = "none")
  expect_equal(a$digest, b$digest)
})

test_that("C3: the package is loaded before the floor is taken", {
  ref <- fake_ref()
  op <- list(
    min_version = "0.0.0",
    setup = function(data, cell) isNamespaceLoaded("randomPlantedForest"),
    run = function(state, cell) state
  )
  m <- bench_measure(ref$lib, op, measure_cell, measure_data(), metric = "none")
  expect_equal(m$digest, bench_digest(TRUE))
})

test_that("C3: peak memory comes from the first timed call only", {
  skip_if(is.null(.bench_cgroup_dir()))
  ref <- fake_ref()
  # each call keeps 200 MB alive, like a forest in C++ memory the GC cannot see
  hoard <- list(
    min_version = "0.0.0",
    setup = function(data, cell) new.env(),
    run = function(state, cell) {
      x <- numeric(25e6)
      x[1] <- 1
      state[[as.character(length(ls(state)))]] <- x
      NULL
    }
  )
  m <- bench_measure(ref$lib, hoard, list(iters = 3L), measure_data(), metric = "cgroup")
  expect_null(m$error)
  expect_gt(m$peak_mb - m$floor_mb, 150)
  expect_lt(m$peak_mb - m$floor_mb, 350)
})

test_that("I2: a child that dies is recorded for that ref, not fatal", {
  ref <- fake_ref()
  crash <- list(
    min_version = "0.0.0",
    setup = function(data, cell) NULL,
    run = function(state, cell) tools::pskill(Sys.getpid(), tools::SIGKILL)
  )
  m <- bench_measure(ref$lib, crash, measure_cell, measure_data(), metric = "none")
  expect_match(m$error, "child died")
})

test_that("I2: a measurement without the ref's package still stops the run", {
  empty <- tempfile("emptylib")
  dir.create(empty)
  expect_error(bench_measure(empty, alloc_op, measure_cell, measure_data(), metric = "none"), "randomPlantedForest")
})

## ---- review fixes: cells, registry, scope ---------------------------------------------------

test_that("C1: output that changes between rounds is flagged", {
  refs <- two_refs()
  ops <- list(
    pid = list(
      min_version = "0.0.0",
      setup = function(data, cell) NULL,
      run = function(state, cell) Sys.getpid()
    )
  )
  rows <- bench_run_cell(fake_cell("pid"), refs["a"], tempfile("data"), ops = ops, metric = "none")
  expect_match(rows$error, "differs between rounds")
})

test_that("C2: Slurm jobs choose their memory metric on the compute node", {
  expect_null(bench_job_metric("slurm"))
  expect_true(bench_job_metric("local") %in% c("cgroup", "rss", "none"))
})

test_that("I4: datasets are written before any job starts, atomically", {
  cells <- rbind(fake_cell("any"), fake_cell("new"))
  cells$n[2] <- 200L
  dir <- tempfile("data")
  paths <- bench_prepare_data(cells, dir)
  expect_length(unique(paths), 2L)
  expect_true(all(file.exists(paths)))
  expect_length(list.files(dir), 2L) # no temporary files left behind
})

test_that("I3: collect reports jobs without results and keeps them out of the history", {
  refs <- two_refs()
  cells <- rbind(fake_cell("any"), fake_cell("missing_op"))
  reg <- bench_make_registry(tempfile("reg"), "local")
  suppressMessages(bench_submit(reg, cells, refs, ops = fake_ops, metric = "none"))
  rows <- suppressWarnings(bench_collect(reg))
  missing <- attr(rows, "missing")
  expect_equal(nrow(missing), 1L)
  expect_equal(missing$op, "missing_op")
  dir <- tempfile("results")
  paths <- bench_write_results(rows, dir, baseline = "HEAD")
  expect_match(paste(readLines(paths[2]), collapse = "\n"), "Missing", fixed = TRUE)
  hist <- tempfile(fileext = ".csv")
  expect_message(bench_history_append(rows, anchors = refs$a$spec, path = hist), "incomplete")
  expect_false(file.exists(hist))
})

test_that("I5: resubmit doubles each job's own memory", {
  refs <- two_refs()
  cells <- rbind(fake_cell("missing_op"), fake_cell("missing_op"))
  cells$n <- c(50000L, 1000L)
  cells$n_families <- c(200L, 20L)
  reg <- bench_make_registry(tempfile("reg"), "local")
  suppressMessages(bench_submit(reg, cells, refs, ops = fake_ops, metric = "none"))
  before <- vapply(batchtools::getJobResources(1:2, reg = reg)$resources, `[[`, numeric(1), "memory")
  expect_false(before[1] == before[2])
  # expired and failed ids come unsorted in practice
  expect_message(bench_resubmit(reg, ids = c(2L, 1L), mem_factor = 2), "--yes")
  after <- vapply(batchtools::getJobResources(1:2, reg = reg)$resources, `[[`, numeric(1), "memory")
  expect_equal(after, before)
  suppressMessages(bench_resubmit(reg, ids = c(2L, 1L), mem_factor = 2, yes = TRUE))
  after <- vapply(batchtools::getJobResources(1:2, reg = reg)$resources, `[[`, numeric(1), "memory")
  expect_equal(after, 2 * before)
})

test_that("resubmit scales walltime on its own, and refuses to change nothing", {
  refs <- two_refs()
  reg <- bench_make_registry(tempfile("reg"), "local")
  suppressMessages(bench_submit(reg, fake_cell("missing_op"), refs, ops = fake_ops, metric = "none"))
  res <- function(what) batchtools::getJobResources(1L, reg = reg)$resources[[1]][[what]]
  mem <- res("memory")
  wall <- res("walltime")
  expect_message(bench_resubmit(reg, ids = 1L, walltime_factor = 3), "walltime 8 -> 24 h", fixed = TRUE)
  suppressMessages(bench_resubmit(reg, ids = 1L, walltime_factor = 3, yes = TRUE))
  expect_equal(res("walltime"), 3 * wall)
  expect_equal(res("memory"), mem)
  expect_error(bench_resubmit(reg, ids = 1L), "--mem-factor or --walltime-factor", fixed = TRUE)
})

test_that("the QoS follows the walltime, not the site's default QoS", {
  refs <- two_refs()
  reg <- bench_make_registry(tempfile("reg"), "local")
  # as in the BIPS site configuration
  reg$default.resources$qos <- "medium"
  suppressMessages(bench_submit(reg, fake_cell("missing_op"), refs, ops = fake_ops, metric = "none"))
  qos <- function() batchtools::getJobResources(1L, reg = reg)$resources[[1]]$qos
  expect_null(qos())
  suppressMessages(bench_resubmit(reg, ids = 1L, walltime_factor = 4, yes = TRUE))
  expect_null(qos())
})

test_that("I6: an unusable systemd-run is detected, not trusted", {
  bin <- tempfile("bin")
  dir.create(bin)
  fake <- file.path(bin, "systemd-run")
  writeLines(c("#!/bin/sh", "echo 'Failed to connect to bus: No medium found' >&2", "exit 1"), fake)
  Sys.chmod(fake, "755")
  old <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, old, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old), add = TRUE)
  expect_false(.bench_systemd_usable())
})

test_that("I1: memory too small to measure says so instead of n/a", {
  rep <- bench_report(fake_rows(mem = c(0, 0)))
  expect_match(rep, "below resolution", all = FALSE, fixed = TRUE)
})

## ---- second review round ------------------------------------------------------------------

test_that("bounds are computed where the op's state exists, for any op", {
  refs <- two_refs()
  ops <- list(
    uses_state = list(
      setup = function(data, cell) list(rows = 5),
      run = function(state, cell) NULL,
      bound = function(state, cell) .bench_estimate("rpf_resources", state$rows)
    )
  )
  rows <- bench_run_cell(fake_cell("uses_state"), refs, tempfile("data"), ops = ops, metric = "none")
  expect_equal(rows$bound_mb, c(NA, 100))
})

test_that("an estimator that fails gives no bound instead of losing the cell", {
  refs <- two_refs()
  ops <- list(
    bad = list(
      setup = function(data, cell) NULL,
      run = function(state, cell) NULL,
      bound = function(state, cell) .bench_estimate("rpf_resources", stop_here = stop("estimator broke"))
    )
  )
  rows <- bench_run_cell(fake_cell("bad"), refs["b"], tempfile("data"), ops = ops, metric = "none")
  expect_true(is.na(rows$bound_mb))
  expect_true(is.na(rows$error))
})

test_that("every op has an estimator hook", {
  for (op in names(bench_ops())) {
    expect_true(is.function(bench_ops()[[op]]$bound), label = op)
  }
})

test_that("a missing hash is n/a, never the same output", {
  rows <- fake_rows(digest = c("x", NA))
  expect_true(is.na(bench_deltas(rows)$output_differs))
  rep <- bench_report(rows)
  expect_match(rep, "| n/a |", all = FALSE, fixed = TRUE)
  # the output column is the last one of a row
  expect_false(any(grepl("\\| same \\|$", strsplit(rep, "\n")[[1]])))
})

test_that("rows from different R minor versions are flagged", {
  rows <- fake_rows()
  rows$r_version <- c("4.6.1", "4.6.2")
  expect_false(grepl("R versions differ", bench_report(rows), fixed = TRUE))
  rows$r_version <- c("4.5.2", "4.6.1")
  expect_match(bench_report(rows), "R versions differ", fixed = TRUE)
})

test_that("ref libraries are kept per R minor version", {
  ref <- fake_ref()
  expect_equal(basename(dirname(ref$lib)), paste0("R-", getRversion()$major, ".", getRversion()$minor))
})

## ---- plots and the HTML report -----------------------------------------------------------

test_that("raw values become one row each, summaries stand in for old results", {
  rows <- fake_rows()
  rows$times <- c("1;1.1;0.9", NA)
  rows$peaks_mb <- c("300;301", NA)
  rows$peak_mb <- c(300.5, 250)
  rows$floor_mb <- c(200, 150)
  rows$time_min <- c(0.9, 1.1)
  rows$time_max <- c(1.1, 1.3)
  long <- bench_long(rows, "time")
  expect_equal(long$value[long$ref == "dev"], c(1, 1.1, 0.9))
  expect_true(all(long$raw[long$ref == "dev"]))
  # no raw times: median, min and max, marked as such
  expect_equal(sort(long$value[long$ref == "pr"]), c(1.1, 1.2, 1.3))
  expect_false(any(long$raw[long$ref == "pr"]))
  # memory is the absolute peak, with the level after setup alongside
  mem <- bench_long(rows, "memory")
  expect_equal(mem$value[mem$ref == "dev"], c(300, 301))
  expect_equal(mem$value[mem$ref == "pr"], 250)
  expect_equal(unique(mem$floor[mem$ref == "pr"]), 150)
  # results written before the raw columns existed
  old <- fake_rows()
  expect_equal(nrow(bench_long(old, "time")), 2L)
})

test_that("plots build for every op and for scaling", {
  rows <- rbind(fake_rows(), fake_rows())
  rows$n[3:4] <- 5000L
  rows$times <- c("1;1.1", "1.2;1.3", "4;4.2", "5;5.1")
  for (what in c("time", "memory")) {
    expect_s3_class(ggplot2::ggplot_build(bench_plot_cells(rows, what)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(bench_plot_scaling(rows, what)), "ggplot_built")
  }
  # a setting measured at one size has no line to draw
  single <- fake_rows(refs = c("dev", "pr"))
  single$op <- "predict"
  scaling <- bench_plot_scaling(rbind(rows, single), "time")
  expect_false("predict" %in% scaling$data$op)
  # boxes only from raw values, summaries stay points
  rows$times[1:2] <- NA
  boxes <- ggplot2::layer_data(bench_plot_cells(rows, "time"), 1L)
  expect_equal(nrow(boxes), 2L)
})

test_that("scaling can run over the predicted rows", {
  rows <- rbind(fake_rows(), fake_rows())
  rows$op <- "predict"
  rows$n_new[3:4] <- 10000L
  rows$times <- c("1;1.1", "1.2;1.3", "4;4.2", "5;5.1")
  p <- bench_plot_scaling(rows, "time", over = "n_new")
  expect_setequal(p$data$n_new, c(1000L, 10000L))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("memory scaling shows the absolute peak, which a reused setup cannot hide", {
  rows <- rbind(fake_rows(mem = c(0, 0)), fake_rows(mem = c(0, 0)))
  rows$n[3:4] <- 5000L
  rows$peaks_mb <- c("300;302", "310;312", "400;402", "410;412")
  p <- bench_plot_scaling(rows, "memory")
  expect_equal(nrow(p$data), 4L)
  expect_true(all(p$data$value >= 300))
})

test_that("a cell records the absolute peak and the level after setup of every round", {
  rows <- bench_run_cell(fake_cell("any"), two_refs(), tempfile("data"), ops = fake_ops, metric = "none")
  expect_equal(lengths(strsplit(rows$peaks_mb, ";")), c(2L, 2L))
  expect_equal(lengths(strsplit(rows$floors_mb, ";")), c(2L, 2L))
})

test_that("the setup line sits at the lowest round, so no peak falls below it", {
  rows <- fake_rows()
  rows$peaks_mb <- c("250;310", "300;301")
  rows$floors_mb <- c("240;300", "290;295")
  rows$floor_mb <- c(270, 292.5)
  mem <- bench_long(rows, "memory")
  expect_equal(unique(mem$floor[mem$ref == "dev"]), 240)
  expect_true(all(mem$value >= mem$floor))
})

test_that("the report says how to read its numbers", {
  rep <- bench_report(fake_rows(), baseline = "dev")
  expect_match(rep, "## How to read this", all = FALSE, fixed = TRUE)
  expect_match(rep, "level after setup", all = FALSE, fixed = TRUE)
  expect_match(rep, "training rows", all = FALSE, fixed = TRUE)
})

test_that("the report format is checked before anything is rendered", {
  expect_error(bench_render_report("x.csv", format = "pdf"), "Unknown format 'pdf'. Use one of md, html.", fixed = TRUE)
})

test_that("the HTML report renders from a rows file", {
  skip_if(!nzchar(Sys.which("quarto")), "quarto not on PATH")
  rows <- rbind(fake_rows(), fake_rows())
  rows$n[3:4] <- 5000L
  rows$times <- c("1;1.1", "1.2;1.3", "4;4.2", "5;5.1")
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(rows, csv, row.names = FALSE)
  out <- suppressMessages(bench_render_report(csv, baseline = "dev", format = "html", dir = tempfile("html")))
  expect_true(file.exists(out))
  expect_match(paste(readLines(out, warn = FALSE), collapse = "\n"), "fit_fixed", fixed = TRUE)
})

test_that("rows from before n_new read as 1000 predicted rows", {
  old <- fake_rows()
  old$n_new <- NULL
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(old, csv, row.names = FALSE)
  rows <- bench_read_rows(csv)
  expect_equal(rows$n_new, c(1000L, 1000L))
  expect_named(rows, bench_schema())
})

test_that("cell labels name the predicted rows for prediction ops only", {
  rows <- rbind(fake_rows(), fake_rows())
  rows$op[3:4] <- "predict"
  label <- .bench_cell_label(rows)
  expect_false(any(grepl("predicted=", label[1:2])))
  expect_true(all(grepl("predicted=1000", label[3:4], fixed = TRUE)))
})
