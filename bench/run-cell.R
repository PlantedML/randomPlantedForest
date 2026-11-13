# One cell measured against every ref, in one process tree on one machine. The
# ref loop lives here, not in the batchtools grid: refs on different nodes
# would bring back the drift that same-machine comparison removes.

bench_schema <- function() {
  c(
    "ref",
    "rpf_version",
    "commit",
    BENCH_CELL_COLS,
    "metric",
    "peak_mb",
    "floor_mb",
    "peak_delta_mb",
    "mem_limit_mb",
    "time_median",
    "time_min",
    "time_max",
    "digest",
    "skipped",
    "error",
    "bound_mb",
    "tier",
    "host",
    "kernel",
    "r_version",
    "job_id",
    "timestamp"
  )
}

# The cgroup's memory limit in MB, walking up because Slurm sets it on the job
# cgroup while the process sits in a child. NA when unlimited.
.bench_cgroup_limit_mb <- function(dir = .bench_cgroup_dir()) {
  if (is.null(dir)) {
    return(NA_real_)
  }
  lim <- Inf
  d <- dir
  repeat {
    for (nm in c("memory.max", "memory.high")) {
      v <- suppressWarnings(as.numeric(.bench_read(file.path(d, nm))[1]))
      if (length(v) && !is.na(v) && v > 0) {
        lim <- min(lim, v / 1024^2)
      }
    }
    if (nchar(d) <= nchar("/sys/fs/cgroup")) {
      break
    }
    d <- dirname(d)
  }
  if (is.finite(lim)) lim else NA_real_
}

.bench_cell_data <- function(cell, data_dir) {
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(data_dir, sprintf("%s-n%d-p%d.rds", cell$task, cell$n, cell$p))
  if (!file.exists(path)) {
    # rename is atomic: a reader sees no file or a whole one
    tmp <- tempfile("data-", tmpdir = data_dir, fileext = ".tmp")
    saveRDS(bench_make_data(cell$task, cell$n, cell$p), tmp)
    file.rename(tmp, path)
  }
  path
}

bench_run_cell <- function(cell, refs, data_dir, ops = bench_ops(), metric = bench_metric()) {
  op <- ops[[cell$op]]
  if (is.null(op)) {
    stop("Unknown op '", cell$op, "'.", call. = FALSE)
  }
  data_path <- .bench_cell_data(cell, data_dir)
  lacking <- lapply(refs, function(r) setdiff(op$requires, r$args))
  supported <- lengths(lacking) == 0L
  runs <- lapply(refs, function(r) list())
  idx <- which(supported)
  for (round in seq_len(cell$mem_reps)) {
    # not sample(idx): with one index i, sample(i) permutes 1:i
    for (i in idx[sample.int(length(idx))]) {
      runs[[i]][[round]] <- bench_measure(refs[[i]]$lib, op, cell, data_path, metric)
    }
  }
  stamp <- list(
    metric = metric,
    mem_limit_mb = .bench_cgroup_limit_mb(),
    tier = cell$tier %||% NA_character_,
    host = Sys.info()[["nodename"]],
    kernel = Sys.info()[["release"]],
    r_version = format(getRversion()),
    job_id = if (nzchar(Sys.getenv("SLURM_JOB_ID"))) Sys.getenv("SLURM_JOB_ID") else NA_character_,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  rows <- lapply(seq_along(refs), function(i) {
    r <- refs[[i]]
    ms <- runs[[i]]
    times <- unlist(lapply(ms, `[[`, "time"))
    errors <- unlist(lapply(ms, `[[`, "error"))
    digests <- unique(stats::na.omit(unlist(lapply(ms, `[[`, "digest"))))
    if (length(digests) > 1L) {
      errors <- c(errors, "output differs between rounds (not deterministic)")
    }
    num <- function(x) if (length(x)) x else NA_real_
    med <- function(f) num(stats::median(vapply(ms, f, numeric(1)), na.rm = TRUE))
    row <- c(
      list(ref = r$spec, rpf_version = r$version, commit = r$commit),
      as.list(cell[BENCH_CELL_COLS]),
      stamp["metric"],
      list(
        peak_mb = if (length(ms)) med(function(m) m$peak_mb) else NA_real_,
        floor_mb = if (length(ms)) med(function(m) m$floor_mb) else NA_real_,
        peak_delta_mb = if (length(ms)) med(function(m) m$peak_mb - m$floor_mb) else NA_real_
      ),
      stamp["mem_limit_mb"],
      list(
        time_median = if (length(times)) stats::median(times) else NA_real_,
        time_min = if (length(times)) min(times) else NA_real_,
        time_max = if (length(times)) max(times) else NA_real_,
        digest = if (length(ms)) ms[[length(ms)]]$digest else NA_character_,
        skipped = if (supported[i]) {
          NA_character_
        } else {
          paste0("needs rpf() argument ", paste(lacking[[i]], collapse = ", "))
        },
        error = if (length(errors)) errors[1] else NA_character_,
        bound_mb = if (length(ms)) med(function(m) m$bound_mb) else NA_real_
      ),
      stamp[c("tier", "host", "kernel", "r_version", "job_id", "timestamp")]
    )
    as.data.frame(row, stringsAsFactors = FALSE)[bench_schema()]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
