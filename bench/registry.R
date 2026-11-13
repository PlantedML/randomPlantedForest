# batchtools runs the cells, locally one after another or one Slurm job each.

BENCH_SOURCES <- c("refs.R", "measure.R", "ops.R", "cells.R", "run-cell.R", "collate.R")
BENCH_WALLTIME <- 28800L # seconds

bench_make_registry <- function(dir, cluster = c("local", "slurm")) {
  cluster <- match.arg(cluster)
  dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
  srcs <- file.path(BENCH_DIR, BENCH_SOURCES)
  if (cluster == "slurm") {
    # The site config (/etc/xdg/batchtools/config.R) supplies the template.
    reg <- batchtools::makeRegistry(file.dir = dir, source = srcs, packages = character(), make.default = FALSE)
    if (is.null(reg$cluster.functions) || identical(reg$cluster.functions$name, "Interactive")) {
      stop("No Slurm cluster functions: batchtools read no site configuration.", call. = FALSE)
    }
    return(reg)
  }
  reg <- batchtools::makeRegistry(
    file.dir = dir,
    source = srcs,
    packages = character(),
    make.default = FALSE,
    conf.file = NA_character_
  )
  # In-process and one at a time: concurrent cells would disturb each other.
  reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive()
  reg
}

# Slurm jobs start a plain Rscript, so the submitting session's library paths
# (rig's environment) go first again before anything is measured.
.bench_job <- function(i, cells, refs, libs, data_dir, ops, metric) {
  .libPaths(unique(c(libs, .libPaths())))
  bench_run_cell(cells[i, ], refs, data_dir, ops = ops %||% bench_ops(), metric = metric %||% bench_metric())
}

# A Slurm job decides on its node: the submit host has neither a job cgroup nor
# a scope, so its choice would be "none" for every job.
bench_job_metric <- function(cluster) {
  if (cluster == "local") bench_metric() else NULL
}

# Every dataset once, before any job starts, so concurrent jobs never read a
# half-written file.
bench_prepare_data <- function(cells, data_dir) {
  vapply(seq_len(nrow(cells)), function(i) .bench_cell_data(cells[i, ], data_dir), character(1))
}

# Conservative defaults, refined once measured peaks exist: memory scales with
# rows times families, one CPU more than rpf's threads for the sampler.
bench_cell_resources <- function(cells) {
  work <- cells$n * cells$n_families / (10000 * 50)
  data.frame(
    ncpus = pmax(cells$nthreads, 1L) + 1L,
    memory = as.integer(pmax(4096, ceiling(4096 * work / 1024) * 1024)),
    walltime = rep(BENCH_WALLTIME, nrow(cells))
  )
}

bench_submit <- function(reg, cells, refs, resources = list(), ops = NULL, metric = NULL) {
  data_dir <- file.path(reg$file.dir, "data")
  bench_prepare_data(cells, data_dir)
  batchtools::batchMap(
    .bench_job,
    i = seq_len(nrow(cells)),
    more.args = list(cells = cells, refs = refs, libs = .libPaths(), data_dir = data_dir, ops = ops, metric = metric),
    reg = reg
  )
  req <- bench_cell_resources(cells)
  for (nm in intersect(names(resources), names(req))) {
    req[[nm]] <- rep(as.integer(resources[[nm]]), nrow(cells))
  }
  # batchtools takes one resource set per submitJobs() call.
  grp <- paste(req$ncpus, req$memory, req$walltime)
  for (g in unique(grp)) {
    i <- which(grp == g)
    batchtools::submitJobs(
      ids = data.frame(job.id = i),
      resources = list(ncpus = req$ncpus[i[1]], memory = req$memory[i[1]], walltime = req$walltime[i[1]]),
      reg = reg
    )
  }
  invisible(reg)
}

# Jobs without results, with their cell, so a partial grid never reads as complete.
.bench_missing_jobs <- function(reg) {
  one <- function(ids, status) data.frame(job.id = ids, status = rep(status, length(ids)))
  status <- rbind(
    one(batchtools::findErrors(reg = reg)$job.id, "failed"),
    one(batchtools::findExpired(reg = reg)$job.id, "expired (walltime or memory)"),
    one(batchtools::findQueued(reg = reg)$job.id, "queued"),
    one(batchtools::findRunning(reg = reg)$job.id, "running"),
    one(batchtools::findNotSubmitted(reg = reg)$job.id, "not submitted")
  )
  status <- status[!duplicated(status$job.id), ]
  if (!nrow(status)) {
    return(cbind(status, data.frame(op = character(), task = character(), n = integer())))
  }
  cells <- lapply(status$job.id, function(id) {
    pars <- batchtools::makeJob(id, reg = reg)$pars
    pars$cells[pars$i, c("op", "task", "n")]
  })
  out <- cbind(status, do.call(rbind, cells))
  rownames(out) <- NULL
  out
}

bench_collect <- function(reg) {
  done <- batchtools::findDone(reg = reg)
  if (!nrow(done)) {
    stop("No job finished. See batchtools::getErrorMessages().", call. = FALSE)
  }
  missing <- .bench_missing_jobs(reg)
  if (nrow(missing)) {
    message(
      nrow(missing),
      " job(s) have no results: ",
      paste(sprintf("%d (%s)", missing$job.id, missing$status), collapse = ", "),
      ". See batchtools::getErrorMessages()."
    )
  }
  out <- do.call(rbind, batchtools::reduceResultsList(ids = done, reg = reg))
  rownames(out) <- NULL
  attr(out, "missing") <- missing
  out
}

# Expired jobs (walltime or memory kill) and failed ones, with more memory.
bench_resubmit <- function(reg, mem_factor = 2, ids = NULL) {
  ids <- ids %||% union(batchtools::findExpired(reg = reg)$job.id, batchtools::findErrors(reg = reg)$job.id)
  if (!length(ids)) {
    message("Nothing to resubmit.")
    return(invisible(integer()))
  }
  # sorted by job id, whatever the order of ids
  res <- batchtools::getJobResources(ids = ids, reg = reg)
  for (k in seq_len(nrow(res))) {
    r <- res$resources[[k]]
    r$memory <- as.integer(r$memory * mem_factor)
    batchtools::submitJobs(ids = data.frame(job.id = res$job.id[k]), resources = r, reg = reg)
  }
  message(length(ids), " job(s) resubmitted with ", mem_factor, " times the memory.")
  invisible(ids)
}

# The newest registry under bench/registry unless one is named.
bench_registry_dir <- function(registry = NA) {
  if (!is.na(registry)) {
    return(registry)
  }
  dirs <- sort(list.dirs("bench/registry", recursive = FALSE), decreasing = TRUE)
  if (!length(dirs)) {
    stop("No registry under bench/registry.", call. = FALSE)
  }
  dirs[1]
}
