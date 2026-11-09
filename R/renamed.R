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
# `.env` (the new name must be an argument of the function and must not be given too), an old `split_structure` value
# is replaced by its new name, and the call gets one deprecation warning for all of them. Any other argument in `...`
# is an error, as rlang::check_dots_empty() gives it.
rpf_renamed <- function(..., ntrees, splits, split_try, t_try, .env, .fn = "rpf") {
  rlang::check_dots_empty(call = .env)
  here <- environment()
  old_given <- names(renamed_args)[
    !vapply(names(renamed_args), function(a) eval(call("missing", as.name(a)), here), logical(1))
  ]
  used <- character(0)
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
    used <- c(used, sprintf("Use {.arg %s} instead of {.arg %s}.", new, old))
  }
  if (exists("split_structure", envir = .env, inherits = FALSE)) {
    s <- get("split_structure", envir = .env)
    if (is.character(s) && length(s) == 1 && s %in% names(renamed_split_structures)) {
      assign("split_structure", renamed_split_structures[[s]], envir = .env)
      used <- c(
        used,
        sprintf(
          "Use {.code split_structure = \"%s\"} instead of {.code \"%s\"}.",
          renamed_split_structures[[s]],
          s
        )
      )
    }
  }
  if (length(used) > 0) {
    fn <- .fn
    names(used) <- rep("i", length(used))
    cli::cli_warn(
      c(
        "{.fn {fn}}: {cli::qty(length(used))}{?a name/names} of randomPlantedForest 0.5.0 {?is/are} deprecated, see {.code ?rpf}.",
        used
      ),
      class = c("rpf_deprecated", "lifecycle_warning_deprecated"),
      call = .env
    )
  }
  invisible(old_given)
}
