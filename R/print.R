#' Print a Random Planted Forest
#'
#' Shows a forest (class `rpf`, from [rpf()] with `prune = FALSE`): the kind of model and its formula, the number of
#' tree families and of splits (with how the splits were chosen), the predictors and the interaction order, whether
#' the forest is purified, the settings of the split search and the number of threads. A pruned model, the default of
#' [rpf()], is printed by [print.rpf_pruned()].
#'
#' @param x `[rpf]`: A fitted [`rpf`] model.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly: `x`.
#' @seealso [rpf()], [print.rpf_pruned()].
#' @export
#'
#' @examples
#' rpf(mpg ~ cyl + wt + drat, data = mtcars, max_interaction = 2, n_families = 10, prune = FALSE)
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
    # n_splits = NULL: chosen out of bag while growing
    oob <- x$splits_oob
    if (!is.null(oob)) {
      at <- match(oob$n_splits, x$oob_path$n_splits)
      # (classification: the loss's criterion, the Brier score or the log-loss)
      crit <- attr(x$oob_path, "criterion")
      err <- if (is.null(crit)) {
        paste("RMSE", format(signif(x$oob_path$oob_rmse[at], 4)))
      } else {
        paste(
          switch(crit, brier = "Brier score", logloss = "log-loss"),
          format(signif(x$oob_path$oob_loss[at], 4))
        )
      }
      cli::cli_alert_info(
        "{.arg n_splits} chosen out of bag (lowest corrected OOB {err}); growing stopped at
        {.val {oob$stopped_at}} split{?s}."
      )
      if (oob$capped) {
        cli::cli_alert_warning("Growing stopped at the cap while the out-of-bag error was still falling.")
      }
    }
    if (!rpf_has_forest(x)) {
      # the fit of a pruned model without its forest: the forest's description only
      cli::cli_alert_info("The forest itself was not kept (pruned with {.code keep_forest = FALSE}).")
    } else if (is_purified(x)) {
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
      isTRUE(params$exact_draw %||% TRUE) &&
        !params$deterministic &&
        params$split_structure == "leaves"
    ) {
      cli::cli_alert_info("Candidates drawn by the exact draw ({.code randomPlantedForest.exact_draw = TRUE}).")
    }
    if (!is.null(params$interaction_constraints)) {
      sets <- vapply(params$interaction_constraints, function(s) paste0("{", paste(s, collapse = ", "), "}"), "")
      cli::cli_text("{.field Interaction constraints}: {.var {sets}}")
    }

    cli::cli_h3("Tree growing")
    # under the argument names (the stored params keep the old ones, renamed_args)
    growing <- params[c(
      "split_structure",
      "split_try",
      "t_try",
      "max_candidates",
      "split_decay_rate",
      "delete_leaves",
      # only when set, so that fits without it print as before
      if (!is.null(params$min_leaf_rows) && params$min_leaf_rows > 1) "min_leaf_rows"
    )]
    growing$split_structure <- split_structure_name(growing$split_structure)
    names(growing)[match(c("split_try", "t_try"), names(growing))] <- renamed_args[c("split_try", "t_try")]
    print_params(growing)
    if (x$mode == "classification") {
      cli::cli_h3("Loss")
      # logit: its exact_step (a model saved by an earlier version has none, and delta / epsilon, which are not shown:
      # its losses were replaced); the pseudo-rows only when set by the option
      loss_params <- c(
        "loss",
        if (identical(params$loss, "logit") && !is.null(params$exact_step)) "exact_step",
        if (params$loss %in% c("logit", "exponential") && !is.null(params$prior_rows)) "prior_rows"
      )
      # the level (classif_level()), on the score scale; multiclass per class, centred as the scores are reported
      level <- NULL
      if (!is.null(x$level)) {
        lev <- x$level
        if (length(lev) > 1) {
          if (identical(params$loss, "logit")) {
            lev <- c(0, lev)
          }
          lev <- stats::setNames(lev - mean(lev), levels(x$blueprint$ptypes$outcomes[[1]]))
          level <- list(level = paste0(names(lev), " ", sprintf("%+.3f", lev), collapse = ", "))
        } else {
          level <- list(level = sprintf("%+.3f", lev))
        }
      }
      print_params(c(params[loss_params], level))
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
#' The flattened forest that [rpf()] stores with `export_forest = TRUE` (as `$forest`) is a long nested list of
#' trees. These methods print a one-line summary instead of the whole list.
#'
#' @param x `[rpf_forest]`: Flattened forest, as in `$forest` of an [`rpf`].
#' @param ... Further arguments passed to or from other methods.
#' @return Invisibly: `x` (`object` for `str()`).
#' @seealso [`rpf`]
#' @export
#' @examples
#'
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, n_families = 10, export_forest = TRUE, prune = FALSE)
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
