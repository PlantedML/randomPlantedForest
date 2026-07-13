# TEMPORARY DEBUG HARNESS — do not merge.
# Localizes a Windows-only native crash: runs every test file in its own R
# process so a crash kills only the child, and BEGIN/END markers survive in
# testthat.Rout even when a child dies mid-file.
library(testthat)
library(randomPlantedForest)

# ---- crash condition matrix (each config in its own child process) ----
matrix_script <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(TRUE)",
  "cfg <- as.integer(args[[1]])",
  "Sys.setenv(RPF_DEBUG = '1')",
  "library(randomPlantedForest)",
  "step <- function(s) { cat('STEP', s, '\\n'); flush(stdout()) }",
  "fit_once <- function(nt, ...) { set.seed(13); rpf(mpg ~ wt + cyl, data = mtcars, nthreads = nt, ...) }",
  "if (cfg == 1) { step('fit ntrees=1 nthreads=2'); f <- fit_once(2, ntrees = 1); step('done') }",
  "if (cfg == 2) { step('fit serial'); f <- fit_once(1); step('purify mode2 nthreads=2'); purify(f, mode = 2L, nthreads = 2L); step('done') }",
  "if (cfg == 3) { step('fit serial'); f <- fit_once(1); step('purify mode1 nthreads=2'); purify(f, mode = 1L, nthreads = 2L); step('done') }",
  "quit(status = 0, save = 'no')"
), matrix_script)

rscript0 <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
cat("\n===== CRASH MATRIX =====\n")
for (cfg in 1:3) {
  cat("----- config", cfg, "-----\n"); flush(stdout())
  st <- suppressWarnings(system2(rscript0, c(matrix_script, cfg)))
  cat("----- config", cfg, "exit status:", st, "-----\n"); flush(stdout())
}
# ---- end crash matrix ----

files <- list.files("testthat", pattern = "^test-.*\\.R$", full.names = TRUE)

child_script <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(TRUE)",
  "library(testthat)",
  "library(randomPlantedForest)",
  "res <- as.data.frame(testthat::test_file(args[[1]], package = 'randomPlantedForest'))",
  "bad <- sum(res$failed) + sum(res$error)",
  "quit(status = as.integer(min(bad, 100L)), save = 'no')"
), child_script)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

statuses <- integer(length(files))
for (i in seq_along(files)) {
  f <- files[[i]]
  cat("\n===== BEGIN", f, "=====\n")
  flush(stdout())
  statuses[i] <- suppressWarnings(system2(rscript, c(child_script, f)))
  cat("===== END", f, "status:", statuses[i], "=====\n")
  flush(stdout())
}

cat("\n===== SUMMARY =====\n")
for (i in seq_along(files)) cat(sprintf("%-45s status %d\n", files[[i]], statuses[i]))
flush(stdout())

if (any(statuses != 0)) {
  stop("At least one test file failed or crashed - see markers above.")
}
