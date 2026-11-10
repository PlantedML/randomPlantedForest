#' Print an rpf
#'
#' @param x `[rpf]`: A fitted [`rpf`] model.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly: `x`.
#' @seealso [`rpf`].
#' @export
#'
#' @examples
#' rpf(mpg ~ cyl + wt + drat, data = mtcars, max_interaction = 2, ntrees = 10)
print.rpf <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @export
format.rpf <- function(x, ...) {
  mode <- switch(x$mode, regression = "Regression", classification = "Classification")
  # the fit's own columns: a formula's terms such as log(x1) are columns of their own, y ~ . - x3 drops one
  fit_predictors <- fit_predictor_names(x)
  p <- length(fit_predictors)
  # max_interaction == p is equivalent to 0: all interactions
  maxint <- if (p == x$params$max_interaction) 0 else x$params$max_interaction
  degree <- switch(
    as.character(maxint),
    "0" = "{.emph all possible interactions}",
    "1" = "{.emph main effects only}",
    "{.emph interactions to degree {.val {maxint}}}"
  )

  params <- x$params

  # format into lines so print() writes to stdout, not cli's message stream
  cli::cli_format_method({
    cli::cli_rule(left = "{mode} Random Planted Forest")
    # only formula blueprints carry a formula; xy and recipe fits don't
    if (is.null(x$blueprint$formula)) {
      predictors <- cli::cli_vec(
        fit_predictors,
        list("vec-trunc" = 5)
      )
      cli::cli_text("{.field Predictors}: {.var {predictors}}")
    } else {
      model_formula <- sub("\\s\\+ 0$", "", deparse1(x$blueprint$formula))
      cli::cli_text("{.field Formula}: {.code {model_formula}}")
    }
    cli::cli_text(
      "{.val {params$ntrees}} tree famil{?y/ies} with ",
      "{.val {params$splits}} split{?s} each on {.val {p}} predictor{?s}, ",
      degree,
      "."
    )
    if (is_purified(x)) {
      cli::cli_alert_success("Forest is purified.")
    } else {
      cli::cli_alert_info("Forest is not purified.")
    }
    # NULL for fits from versions before keep_level_order
    kept <- params$keep_level_order
    if (length(kept) > 0) {
      sorted <- setdiff(names(x$factor_levels), kept)
      # one quantity per cli string, so the two parts are formatted separately
      msg <- cli::format_inline("Factor{?s} {.var {kept}} {?keeps its/keep their} level order")
      if (length(sorted) > 0) {
        msg <- paste0(msg, "; ", cli::format_inline("{.var {sorted}} {?is/are} ordered by the outcome"))
      }
      cli::cli_alert_info("{msg}.")
    }
    if (params$deterministic) {
      cli::cli_alert_warning("Fit deterministically.")
    }
    if (
      isFALSE(params$exact_draw) &&
        !params$deterministic &&
        x$mode == "regression" &&
        params$split_structure == "leaves"
    ) {
      cli::cli_alert_info("Candidates drawn by the fast draw ({.code randomPlantedForest.exact_draw = FALSE}).")
    }
    if (!is.null(params$interaction_constraints)) {
      sets <- vapply(params$interaction_constraints, function(s) paste0("{", paste(s, collapse = ", "), "}"), "")
      cli::cli_text("{.field Interaction constraints}: {.var {sets}}")
    }

    cli::cli_h3("Tree growing")
    print_params(params[c(
      "split_structure",
      "split_try",
      "t_try",
      "max_candidates",
      "split_decay_rate",
      "delete_leaves",
      # only when set, so that fits without it print as before
      if (!is.null(params$min_leaf_rows) && params$min_leaf_rows > 1) "min_leaf_rows"
    )])
    if (x$mode == "classification") {
      cli::cli_h3("Loss")
      loss_params <- if (params$loss %in% c("logit", "exponential")) {
        c("loss", "delta", "epsilon")
      } else {
        "loss"
      }
      print_params(params[loss_params])
    }

    cli::cli_text("")
    cli::cli_alert_info(
      "Fit using {.val {params$nthreads}} thread{?s}, also the default for {.fn predict} and {.fn purify}."
    )
  })
}

print_params <- function(params) {
  values <- vapply(params, format, character(1))
  cli::cli_verbatim(paste0(
    "  ",
    format(names(values), justify = "right"),
    ": ",
    values
  ))
}

#' Compact printing of forest structures
#'
#' These methods are provided to avoid flooding the console with long nested lists containing tree structures.
#'
#' @param x `[rpf_forest]`: Flattened forest, as in `$forest` of an [`rpf`].
#' @param ... Further arguments passed to or from other methods.
#' @seealso [`rpf`]
#' @export
#' @examples
#'
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 10)
#' print(rpfit$forest)
#' str(rpfit$forest)
print.rpf_forest <- function(x, ...) {
  cli::cat_line(cli::format_inline("<rpf_forest> of {length(x)} tree{?s}"))
  invisible(x)
}

#' @rdname print.rpf_forest
#' @param object `[rpf_forest]`: Flattened forest.
#' @export
str.rpf_forest <- function(object, ...) print(object, ...)
