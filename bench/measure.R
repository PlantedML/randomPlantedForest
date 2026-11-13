# Measuring one op for one ref: a fresh R process per measurement, timed calls
# inside it, and memory sampled from outside while it works.

BENCH_DIR <- normalizePath("bench", mustWork = TRUE)

# Thread pools other than rpf's own are pinned, so nthreads is the only knob.
BENCH_THREAD_ENV <- c(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  R_DATATABLE_NUM_THREADS = "1"
)

BENCH_SAMPLE_INTERVAL <- 0.005
# rpf seeds itself from R's RNG, so a fixed seed makes every ref's output comparable.
BENCH_SEED <- 20261113L

## ---- memory metric -------------------------------------------------------------

.bench_read <- function(path) {
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)), error = function(e) character(0))
}

# cgroup v2 directory of this process, or NULL without memory accounting.
.bench_cgroup_dir <- function() {
  cg <- .bench_read("/proc/self/cgroup")
  line <- cg[startsWith(cg, "0::")]
  if (!length(line)) {
    return(NULL)
  }
  dir <- file.path("/sys/fs/cgroup", sub("^0::/?", "", line[1]))
  if (file.exists(file.path(dir, "memory.stat"))) dir else NULL
}

# anon + shmem of the cgroup: every process of the run is charged once.
.bench_cgroup_mb <- function(dir) {
  st <- .bench_read(file.path(dir, "memory.stat"))
  get <- function(key) {
    v <- st[startsWith(st, paste0(key, " "))]
    if (length(v)) as.numeric(sub("^\\S+ ", "", v[1])) else 0
  }
  if (!length(st)) {
    return(NA_real_)
  }
  (get("anon") + get("shmem")) / 1024^2
}

# Only a cgroup that holds this run alone gives clean numbers: a Slurm job, or
# a scope started by bench_ensure_scope(). A desktop session's cgroup does not.
.bench_choose_metric <- function(sysname, cgroup, scoped, slurm, ps_ok) {
  if (sysname == "Linux" && !is.null(cgroup) && (scoped || slurm)) {
    "cgroup"
  } else if (sysname == "Darwin" && ps_ok) {
    "rss"
  } else {
    "none"
  }
}

bench_metric <- function() {
  m <- .bench_choose_metric(
    Sys.info()[["sysname"]],
    cgroup = .bench_cgroup_dir(),
    scoped = nzchar(Sys.getenv("RPF_BENCH_SCOPED")),
    slurm = nzchar(Sys.getenv("SLURM_JOB_ID")),
    ps_ok = ps::ps_is_supported()
  )
  if (m == "none" && !nzchar(Sys.getenv("RPF_BENCH_WARNED"))) {
    Sys.setenv(RPF_BENCH_WARNED = "1")
    warning(
      "Memory is not measured: accurate measurement needs cgroup v2 on Linux ",
      "(a Slurm job, or a local run started through systemd-run).",
      call. = FALSE
    )
  }
  m
}

.bench_needs_scope <- function(sysname, scoped, slurm, has_systemd_run) {
  sysname == "Linux" && !scoped && !slurm && has_systemd_run
}

# The binary alone is not enough: without a user bus (ssh without lingering,
# sudo -u, containers) systemd-run exits with an error.
.bench_systemd_usable <- function() {
  nzchar(Sys.which("systemd-run")) &&
    identical(
      suppressWarnings(system2(
        "systemd-run",
        c("--user", "--scope", "--quiet", "--", "true"),
        stdout = FALSE,
        stderr = FALSE
      )),
      0L
    )
}

# Re-runs the CLI inside its own systemd scope, so the cgroup holds this run
# alone. Returns only when no re-run is needed.
bench_ensure_scope <- function(args) {
  needed <- .bench_needs_scope(
    Sys.info()[["sysname"]],
    scoped = nzchar(Sys.getenv("RPF_BENCH_SCOPED")),
    slurm = nzchar(Sys.getenv("SLURM_JOB_ID")),
    has_systemd_run = .bench_systemd_usable()
  )
  if (!needed) {
    return(invisible())
  }
  message("Re-running inside a systemd scope for clean memory numbers.")
  status <- system2(
    "systemd-run",
    c("--user", "--scope", "--quiet", "--", "env", "RPF_BENCH_SCOPED=1", "rig", "run", "bench/rpf-bench", shQuote(args))
  )
  quit(save = "no", status = status)
}

## ---- output hashes ---------------------------------------------------------------

# Values only: classes and attributes differ between versions (tibble vs data
# frame) without the numbers changing.
.bench_strip <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (is.list(x)) {
    nm <- names(x)
    x <- lapply(unclass(x), .bench_strip)
    names(x) <- nm
    return(x)
  }
  keep <- attributes(x)[intersect(names(attributes(x)), c("names", "dim", "dimnames"))]
  attributes(x) <- if (length(keep)) keep else NULL
  x
}

bench_digest <- function(x) {
  rlang::hash(.bench_strip(x))
}

## ---- one measurement ---------------------------------------------------------------

# Runs in the child. Markers tell the sampling parent which phase it is in.
.bench_child <- function(lib, setup, run, bound, single, cell, data_path, markers, helpers) {
  for (h in helpers) {
    source(h)
  }
  bench_assert_lib(lib)
  tryCatch(
    {
      # loading is not part of any op: about 40 MB and most of a second
      loadNamespace(BENCH_PKG)
      data <- readRDS(data_path)
      set.seed(BENCH_SEED)
      state <- setup(data, cell)
      bound_mb <- if (is.function(bound)) bound(state, cell) else NA_real_
      invisible(gc())
      file.create(file.path(markers, "floor"))
      Sys.sleep(0.25)
      iters <- if (single) 1L else cell$iters
      time <- numeric(iters)
      out <- NULL
      for (i in seq_len(iters)) {
        # memory from the first call only: forests in C++ memory are invisible
        # to R's GC, so later calls would add up earlier ones
        if (i == 1L) {
          file.create(file.path(markers, "timed"))
        }
        set.seed(BENCH_SEED + 1L)
        t0 <- proc.time()[["elapsed"]]
        out <- run(state, cell)
        time[i] <- proc.time()[["elapsed"]] - t0
        if (i == 1L) file.create(file.path(markers, "rest"))
      }
      file.create(file.path(markers, "done"))
      list(time = time, digest = bench_digest(out), bound_mb = bound_mb, error = NULL)
    },
    error = function(e) list(time = numeric(), digest = NA_character_, error = conditionMessage(e))
  )
}

.bench_phase <- function(markers) {
  if (file.exists(file.path(markers, "done"))) {
    "done"
  } else if (file.exists(file.path(markers, "rest"))) {
    "rest"
  } else if (file.exists(file.path(markers, "timed"))) {
    "timed"
  } else if (file.exists(file.path(markers, "floor"))) {
    "floor"
  } else {
    "setup"
  }
}

bench_measure <- function(lib, op, cell, data_path, metric = bench_metric()) {
  markers <- tempfile("bench-markers-")
  dir.create(markers)
  on.exit(unlink(markers, recursive = TRUE), add = TRUE)
  px <- callr::r_bg(
    .bench_child,
    args = list(
      lib = lib,
      setup = op$setup,
      run = op$run,
      bound = op$bound,
      single = isTRUE(op$single),
      cell = cell,
      data_path = data_path,
      markers = markers,
      helpers = file.path(BENCH_DIR, c("refs.R", "measure.R", "ops.R"))
    ),
    libpath = c(lib, .libPaths()),
    env = c(callr::rcmd_safe_env(), BENCH_THREAD_ENV),
    supervise = TRUE
  )
  cgroup <- .bench_cgroup_dir()
  sample_mb <- switch(
    metric,
    cgroup = function() .bench_cgroup_mb(cgroup),
    rss = function() {
      tryCatch(px$get_memory_info()[["rss"]] / 1024^2, error = function(e) NA_real_)
    },
    none = NULL
  )
  values <- numeric()
  phases <- character()
  while (px$is_alive()) {
    if (!is.null(sample_mb)) {
      # phase first: a sample taken after the marker still belongs to that phase
      phases <- c(phases, .bench_phase(markers))
      values <- c(values, sample_mb())
    }
    Sys.sleep(BENCH_SAMPLE_INTERVAL)
  }
  res <- tryCatch(px$get_result(), error = function(e) {
    msg <- paste(conditionMessage(e), conditionMessage(e$parent %||% e))
    # a ref library that does not hold the package is a setup error, not a result
    if (grepl("not to the ref library", msg, fixed = TRUE)) {
      stop(e)
    }
    list(
      time = numeric(),
      digest = NA_character_,
      error = sprintf("child died (exit status %s): %s", px$get_exit_status(), conditionMessage(e))
    )
  })
  floor_mb <- peak_mb <- NA_real_
  if (!is.null(sample_mb) && is.null(res$error)) {
    fl <- values[phases == "floor"]
    tm <- values[phases == "timed"]
    floor_mb <- if (length(fl)) stats::median(fl, na.rm = TRUE) else NA_real_
    peak <- c(tm, utils::tail(fl, 1))
    peak_mb <- if (any(!is.na(peak))) max(peak, na.rm = TRUE) else NA_real_
  }
  list(
    time = res$time,
    peak_mb = peak_mb,
    floor_mb = floor_mb,
    digest = res$digest,
    bound_mb = res$bound_mb %||% NA_real_,
    error = res$error
  )
}
