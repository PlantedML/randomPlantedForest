# The argument names and split structures renamed after randomPlantedForest 0.5.0: old name = new name. The old ones
# still work, with a deprecation warning. A fitted model's `$params` keeps the old names and structure codes: it is
# the record of what the C++ forest was built with, and models saved by 0.5.0 have it so.
renamed_args <- c(
  ntrees = "n_families",
  splits = "n_splits",
  split_try = "n_try_split_points",
  t_try = "share_candidates"
)
renamed_split_structures <- c(
  leaves = "leaf",
  hist = "leaf_binned",
  cur_trees_1 = "tree_then_leaf",
  cur_trees_2 = "tree_then_all_leaves",
  res_trees = "target_tree_then_leaf"
)

# The stored code of a `split_structure` (new name), and the name of a stored code
split_structure_code <- function(split_structure) {
  names(renamed_split_structures)[match(split_structure, renamed_split_structures)]
}
split_structure_name <- function(code) {
  unname(renamed_split_structures[code])
}

# The old argument names of a function, given in its `...`: called first in the function as
# `rpf_renamed(..., .env = environment(), .fn = "<its name>")`. Each old argument given is moved to its new name in
# `.env` (the new name must be an argument of the function and must not be given too), and an old `split_structure`
# value is replaced by its new name. The old names not warned about yet in this session get one deprecation warning
# (rpf_deprecation_warn()). Any other argument in `...` is an error, as rlang::check_dots_empty() gives it.
rpf_renamed <- function(..., ntrees, splits, split_try, t_try, .env, .fn = "rpf") {
  rlang::check_dots_empty(call = .env)
  here <- environment()
  old_given <- names(renamed_args)[
    !vapply(names(renamed_args), function(a) eval(call("missing", as.name(a)), here), logical(1))
  ]
  used <- list()
  for (old in old_given) {
    new <- renamed_args[[old]]
    if (!exists(new, envir = .env, inherits = FALSE)) {
      cli::cli_abort("Unused argument {.arg {old}}.", call = .env)
    }
    if (!eval(call("missing", as.name(new)), .env)) {
      cli::cli_abort(
        c(
          "Give {.arg {new}} only: {.arg {old}} is its old name.",
          "x" = "Both {.arg {new}} and {.arg {old}} were given."
        ),
        call = .env
      )
    }
    assign(new, get(old, envir = here), envir = .env)
    used[[old]] <- sprintf("Use {.arg %s} instead of {.arg %s}.", new, old)
  }
  if (exists("split_structure", envir = .env, inherits = FALSE)) {
    s <- get("split_structure", envir = .env)
    if (is.character(s) && length(s) == 1 && s %in% names(renamed_split_structures)) {
      assign("split_structure", renamed_split_structures[[s]], envir = .env)
      used[[paste0("split_structure = ", s)]] <- sprintf(
        "Use {.code split_structure = \"%s\"} instead of {.code \"%s\"}.",
        renamed_split_structures[[s]],
        s
      )
    }
  }
  rpf_deprecation_warn(used, .fn, call = .env)
  invisible(old_given)
}

# The old names warned about in this session. As lifecycle's deprecate_soft() default, each old name warns once per
# session, so that a caller that fits in a loop (a learner in a benchmark) warns once, not on every fit.
rpf_deprecations <- new.env(parent = emptyenv())

# Forget the old names warned about (for tests)
rpf_deprecation_reset <- function() {
  rm(list = ls(rpf_deprecations, all.names = TRUE), envir = rpf_deprecations)
  invisible(NULL)
}

# One deprecation warning for the old names `used` (a list named by a key per old name, each element a bullet) not
# warned about yet in this session. As for lifecycle, options(lifecycle_verbosity = "warning") warns every time and
# "quiet" never. Returns whether it warned.
rpf_deprecation_warn <- function(used, fn, call) {
  verbosity <- getOption("lifecycle_verbosity", "default")
  if (length(used) == 0 || identical(verbosity, "quiet")) {
    return(invisible(FALSE))
  }
  if (!identical(verbosity, "warning")) {
    used <- used[!names(used) %in% ls(rpf_deprecations, all.names = TRUE)]
    if (length(used) == 0) {
      return(invisible(FALSE))
    }
    for (key in names(used)) {
      assign(key, TRUE, envir = rpf_deprecations)
    }
  }
  bullets <- unlist(used, use.names = FALSE)
  names(bullets) <- rep("i", length(bullets))
  cli::cli_warn(
    c(
      "{.fn {fn}}: {cli::qty(length(bullets))}{?a name/names} of randomPlantedForest 0.5.0 {?is/are} deprecated, see {.code ?rpf}.",
      bullets,
      " " = "This warning is shown once per session for each old name."
    ),
    class = c("rpf_deprecated", "lifecycle_warning_deprecated"),
    call = call
  )
  invisible(TRUE)
}
