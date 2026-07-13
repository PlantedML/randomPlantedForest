# TEMPORARY DEBUG HARNESS — do not merge.
# Localizes a Windows-only native crash: runs every test file in its own R
# process so a crash kills only the child, and BEGIN/END markers survive in
# testthat.Rout even when a child dies mid-file.
library(testthat)
library(randomPlantedForest)

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
