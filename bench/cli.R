#| name: rpf-bench
#| title: randomPlantedForest benchmark suite
#| description: |
#|   Compare randomPlantedForest refs for runtime, peak memory and output.
#|   Run from the repository root via `rig run bench/rpf-bench`.
#|   See bench/README.md.

for (f in setdiff(Sys.glob("bench/*.R"), c("bench/cli.R", "bench/tests.R"))) {
  source(f)
}

split_opt <- function(x) if (is.na(x) || !nzchar(x)) NA else trimws(strsplit(x, ",")[[1]])

switch(
  command <- "",

  #| title: Measure refs on the selected cells
  #| description: |
  #|   Installs each ref into bench/lib and measures every cell against every
  #|   ref. Locally it waits and writes bench/results/report.md. On slurm it
  #|   submits and exits: run `collect` once the jobs are done.
  run = {
    #| description: "Which grid: quick (small, for a laptop) or full (large, for the cluster)."
    #| short: t
    tier <- "quick"

    #| description: "Where cells run: local (one after another, in this session) or slurm (one job per cell)."
    #| short: b
    backend <- "local"

    #| description: "Refs to compare, comma separated. Default: WORKTREE,dev, plus the anchors in bench/anchors.csv on the full tier."
    #| short: r
    refs <- ""

    #| description: "Ref the others are compared against."
    baseline <- "dev"

    #| description: "Only these ops, comma separated (fit_fixed, fit_oob, fit_default, purify, predict, predict_components)."
    #| short: o
    ops <- ""

    #| description: "Only these tasks, comma separated (regr, binary, multiclass)."
    tasks <- ""

    #| description: "Only these numbers of rows, comma separated."
    sizes <- ""

    #| description: "Only the first N cells, to check a setup cheaply."
    #| short: n
    max_cells <- NA_integer_

    #| description: "Timed calls per measurement, instead of the tier's default."
    iters <- NA_integer_

    #| description: "Rounds per ref, instead of the tier's default."
    mem_reps <- NA_integer_

    #| description: "Where ref libraries are installed. Must be readable from the compute nodes."
    lib_root <- "bench/lib"

    #| description: "Print what would run and exit, without installing or measuring."
    dry_run <- FALSE

    # everything that can be wrong with the inputs fails here, before any install
    bench_check_choice(backend, c("local", "slurm"), "backend")
    cells <- bench_select_cells(
      bench_cells(tier),
      ops = split_opt(ops),
      tasks = split_opt(tasks),
      sizes = split_opt(sizes),
      max_cells = max_cells
    )
    if (!is.na(iters)) {
      cells$iters <- iters
    }
    if (!is.na(mem_reps)) {
      cells$mem_reps <- mem_reps
    }
    specs <- split_opt(refs)
    if (all(is.na(specs))) {
      specs <- c("WORKTREE", "dev")
      if (tier == "full" && file.exists("bench/anchors.csv")) {
        specs <- unique(c(specs, utils::read.csv("bench/anchors.csv")$ref))
      }
    }
    specs <- bench_order_refs(specs, baseline)
    message(sprintf(
      "rpf bench | tier %s | backend %s | %d cells | refs %s | baseline %s",
      tier,
      backend,
      nrow(cells),
      paste(specs, collapse = ", "),
      baseline
    ))
    if (dry_run) {
      print(cells[, c("op", "task", "loss", "n", "nthreads", "iters", "mem_reps")], row.names = FALSE)
      quit(save = "no", status = 0L)
    }

    if (backend == "local") {
      bench_ensure_scope(commandArgs(TRUE))
    }
    bench_assert_clean_tree(specs)
    specs <- bench_resolve_refs(specs)
    if (length(specs) < 2L) {
      message("Only one distinct ref (", specs, "), so there is nothing to compare. Name refs with --refs.")
    }
    installed <- lapply(specs, bench_install_ref, root = lib_root)
    metric <- bench_job_metric(backend)
    message("Memory metric: ", metric %||% "chosen per job")

    dir <- file.path("bench/registry", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", tier))
    reg <- bench_make_registry(dir, backend)
    bench_submit(reg, cells, installed, metric = metric)
    if (backend == "local") {
      paths <- bench_write_results(bench_collect(reg), baseline = baseline)
      message("Report: ", paths[2])
    } else {
      message("Submitted to ", dir, ". Run `collect` when the jobs are done.")
    }
  },

  #| title: Collect finished jobs into results and a report
  collect = {
    #| description: "Registry directory. Default: the newest under bench/registry."
    registry <- ""

    #| description: "Ref the others are compared against."
    baseline <- "dev"

    reg <- batchtools::loadRegistry(bench_registry_dir(if (nzchar(registry)) registry else NA), writeable = FALSE)
    rows <- bench_collect(reg)
    paths <- bench_write_results(rows, baseline = baseline)
    if (file.exists("bench/anchors.csv")) {
      bench_history_append(rows, utils::read.csv("bench/anchors.csv")$ref)
    }
    message("Report: ", paths[2])
  },

  #| title: Resubmit expired or failed jobs with more memory
  #| description: |
  #|   Lists the jobs and their new memory request. Submits only with --yes.
  resubmit = {
    #| description: "Registry directory. Default: the newest under bench/registry."
    registry <- ""

    #| description: "Multiply the memory request by this factor."
    mem_factor <- 2

    #| description: "Submit. Without it, only list the jobs and their new memory."
    yes <- FALSE

    reg <- batchtools::loadRegistry(bench_registry_dir(if (nzchar(registry)) registry else NA), writeable = TRUE)
    bench_resubmit(reg, mem_factor = mem_factor, yes = yes)
  },

  #| title: Run the suite's own tests
  test = {
    res <- as.data.frame(testthat::test_file("bench/tests.R", reporter = "summary"))
    quit(save = "no", status = as.integer(any(res$failed > 0 | res$error)))
  },

  #| title: Render the report from the newest results
  report = {
    #| description: "Ref the others are compared against."
    baseline <- "dev"

    csvs <- sort(Sys.glob("bench/results/rows-*.csv"), decreasing = TRUE)
    if (!length(csvs)) {
      stop("No results under bench/results.", call. = FALSE)
    }
    rows <- utils::read.csv(csvs[1], stringsAsFactors = FALSE)
    cat(bench_report(rows, baseline), "\n")
  }
)
