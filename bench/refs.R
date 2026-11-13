# Refs: which versions of randomPlantedForest are measured, and how each gets a
# library of its own. A ref is a git rev or "WORKTREE" (uncommitted changes).

BENCH_PKG <- "randomPlantedForest"
# What makes two refs the same package for measuring purposes.
BENCH_PKG_PATHS <- c("DESCRIPTION", "NAMESPACE", "R", "src", "inst")

.bench_git <- function(repo, ..., ok_status = 0L) {
  out <- suppressWarnings(system2("git", c("-C", repo, ...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  if (!status %in% ok_status) {
    stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bench_parse_ref <- function(spec) {
  if (identical(spec, "WORKTREE")) {
    list(type = "worktree", rev = "HEAD")
  } else {
    list(type = "git", rev = spec)
  }
}

# A ref that is not a local commit falls back to origin/<ref>, so remote
# branches like dev work without a local checkout.
.bench_rev <- function(spec, repo) {
  rev <- bench_parse_ref(spec)$rev
  for (cand in c(rev, paste0("origin/", rev))) {
    ok <- .bench_git(repo, "rev-parse", "--verify", "--quiet", paste0(cand, "^{commit}"), ok_status = c(0L, 1L))
    if (length(ok) && nzchar(ok[1])) {
      return(cand)
    }
  }
  stop("Unknown ref '", spec, "': neither a commit nor a branch on origin.", call. = FALSE)
}

# Package files of the working tree: tracked and untracked but not ignored, so
# build products like src/*.o never count.
.bench_worktree_files <- function(repo) {
  files <- .bench_git(repo, "ls-files", "-co", "--exclude-standard", "--", BENCH_PKG_PATHS)
  files[file.exists(file.path(repo, files))]
}

# "path blob" lines for the package files of a ref, in git's blob ids either way,
# so WORKTREE and an identical commit get the same id.
.bench_ref_blobs <- function(spec, repo) {
  ref <- bench_parse_ref(spec)
  if (ref$type == "git") {
    out <- .bench_git(repo, "ls-tree", "-r", .bench_rev(spec, repo), "--", BENCH_PKG_PATHS)
    # "<mode> blob <sha>\t<path>"
    sub("^\\S+ \\S+ (\\S+)\t(.*)$", "\\2 \\1", out)
  } else {
    files <- .bench_worktree_files(repo)
    blobs <- .bench_git(repo, "hash-object", "--", files)
    paste(files, blobs)
  }
}

bench_ref_id <- function(spec, repo = ".") {
  substr(rlang::hash(sort(.bench_ref_blobs(spec, repo))), 1, 12)
}

# Baseline first, so that deduplication keeps its name.
bench_order_refs <- function(specs, baseline) {
  if (!baseline %in% specs) {
    stop(
      sprintf(
        "Baseline '%s' is not among the refs (%s). Add it to --refs or choose one of them with --baseline.",
        baseline,
        paste(specs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  c(baseline, setdiff(specs, baseline))
}

bench_resolve_refs <- function(specs, repo = ".") {
  ids <- vapply(specs, bench_ref_id, character(1), repo = repo)
  dup <- duplicated(ids)
  if (any(dup)) {
    for (i in which(dup)) {
      message(
        "Ref ",
        specs[i],
        " has the same package content as ",
        specs[match(ids[i], ids)],
        " and is measured once."
      )
    }
  }
  specs[!dup]
}

# Source tree of a ref in a fresh directory: a detached worktree for git refs,
# a copy of the package files for WORKTREE.
.bench_ref_source <- function(spec, repo, dir) {
  ref <- bench_parse_ref(spec)
  if (ref$type == "git") {
    .bench_git(repo, "worktree", "add", "--detach", "-f", dir, .bench_rev(spec, repo))
  } else {
    files <- .bench_worktree_files(repo)
    for (f in files) {
      dir.create(file.path(dir, dirname(f)), recursive = TRUE, showWarnings = FALSE)
      file.copy(file.path(repo, f), file.path(dir, f))
    }
  }
  dir
}

.bench_ref_commit <- function(spec, repo) {
  ref <- bench_parse_ref(spec)
  sha <- .bench_git(repo, "rev-parse", paste0(.bench_rev(spec, repo), "^{commit}"))
  if (ref$type == "worktree") paste0("WORKTREE@", sha) else sha
}

bench_install_ref <- function(spec, root = "bench/lib", repo = ".") {
  repo <- normalizePath(repo, mustWork = TRUE)
  id <- bench_ref_id(spec, repo)
  # a library built by one R minor version does not belong to another
  root <- file.path(root, paste0("R-", R.version$major, ".", strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1]))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root)
  lib <- file.path(root, id)
  log <- paste0(lib, ".log")
  desc <- file.path(lib, BENCH_PKG, "DESCRIPTION")
  if (!file.exists(desc)) {
    unlink(lib, recursive = TRUE)
    dir.create(lib)
    src <- file.path(tempfile("rpf-src-"), BENCH_PKG)
    dir.create(dirname(src))
    .bench_ref_source(spec, repo, src)
    on.exit(
      {
        if (bench_parse_ref(spec)$type == "git") {
          .bench_git(repo, "worktree", "remove", "--force", src, ok_status = c(0L, 128L))
        }
        unlink(dirname(src), recursive = TRUE)
      },
      add = TRUE
    )
    status <- system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "INSTALL", "--no-test-load", "--no-multiarch", "-l", shQuote(lib), shQuote(src)),
      stdout = log,
      stderr = log,
      env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
    )
    if (status != 0L || !file.exists(desc)) {
      unlink(lib, recursive = TRUE)
      stop(
        "Installing ref '",
        spec,
        "' failed. Last lines of ",
        log,
        ":\n",
        paste(utils::tail(readLines(log), 20), collapse = "\n"),
        call. = FALSE
      )
    }
  }
  list(
    spec = spec,
    id = id,
    lib = lib,
    version = unname(read.dcf(desc, fields = "Version")[1, 1]),
    commit = .bench_ref_commit(spec, repo),
    args = callr::r(.bench_rpf_formals, libpath = c(lib, .libPaths()))
  )
}

# What rpf() accepts decides which API a ref has. The version number does not:
# 0.5.0.9000 predates the argument renaming.
.bench_rpf_formals <- function() {
  ns <- asNamespace("randomPlantedForest")
  f <- if (exists("rpf.data.frame", ns, inherits = FALSE)) get("rpf.data.frame", ns) else get("rpf", ns)
  names(formals(f))
}

# .libPaths() silently drops a missing directory, so a process could load the
# session's package instead of the ref's. Every process that measures checks.
bench_assert_lib <- function(lib) {
  from <- tryCatch(normalizePath(find.package(BENCH_PKG)), error = function(e) "<not found>")
  if (!startsWith(from, normalizePath(lib, mustWork = FALSE))) {
    stop(BENCH_PKG, " resolves to ", from, ", not to the ref library ", lib, call. = FALSE)
  }
  invisible(TRUE)
}

bench_assert_clean_tree <- function(refs, repo = ".") {
  dirty <- .bench_git(repo, "status", "--porcelain", "--untracked-files=no", "--", BENCH_PKG_PATHS)
  if (length(dirty) && !"WORKTREE" %in% refs) {
    stop(
      "The working tree has uncommitted package changes:\n",
      paste(dirty, collapse = "\n"),
      "\nCommit them, or measure them on purpose with the ref WORKTREE.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
