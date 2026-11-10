# rpf has no missing value handling: splits compare feature values directly,
# so an NA would silently send a row through arbitrary leaves.
check_no_missing <- function(data, what, call = rlang::caller_env()) {
  missing_cols <- names(data)[vapply(data, anyNA, logical(1))]
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      c(
        paste(what, "must not contain missing values."),
        "x" = "Missing values in {.var {missing_cols}}.",
        "i" = "Impute or remove them first, e.g. with a recipe passed to {.fn rpf}."
      ),
      call = call
    )
  }
  invisible(data)
}

abort_outcome_type <- function(y, call = rlang::caller_env()) {
  cli::cli_abort(
    c(
      "The outcome must be numeric (regression) or a factor (classification).",
      "x" = "It is of class {.cls {class(y)}}.",
      "i" = if (is.character(y) || is.logical(y)) "For classification, convert it with {.fn factor}."
    ),
    call = call
  )
}

# The x/y interface: hardhat::mold() cannot standardize an outcome vector that is
# neither numeric nor a factor and would fail with its own message first.
check_outcome_type <- function(y, call = rlang::caller_env()) {
  if (is.atomic(y) && is.null(dim(y)) && !is.numeric(y) && !is.factor(y)) {
    abort_outcome_type(y, call = call)
  }
  invisible(y)
}

# Checks the molded outcome before it is encoded: numeric outcomes must be finite
# (NA is reported by check_no_missing()), factor outcomes need two classes. Levels
# that do not occur are dropped with a warning, in the outcome and in the
# blueprint's ptype, so that predict() returns one column per class of the fit.
check_outcome <- function(processed, call = rlang::caller_env()) {
  y <- processed$outcomes[[1]]
  if (is.factor(y)) {
    used <- levels(y)[levels(y) %in% y]
    if (length(used) < 2) {
      cli::cli_abort(
        c(
          "A factor outcome needs at least two classes.",
          "x" = "Only {.val {used}} occurs."
        ),
        call = call
      )
    }
    unused <- setdiff(levels(y), used)
    if (length(unused) > 0) {
      cli::cli_warn(
        c(
          "Dropped {length(unused)} unused level{?s} of the outcome: {.val {unused}}.",
          "i" = "The forest predicts the {length(used)} classes that occur."
        ),
        call = call
      )
      processed$outcomes[[1]] <- factor(y, levels = used)
      ptypes <- processed$blueprint$ptypes
      ptypes$outcomes[[1]] <- factor(ptypes$outcomes[[1]], levels = used)
      processed$blueprint <- hardhat::update_blueprint(processed$blueprint, ptypes = ptypes)
    }
  } else if (checkmate::test_numeric(y)) {
    if (any(is.infinite(y))) {
      cli::cli_abort(
        c(
          "The outcome must be finite.",
          "x" = "It has {sum(is.infinite(y))} infinite value{?s}."
        ),
        call = call
      )
    }
  } else {
    abort_outcome_type(y, call = call)
  }
  processed
}

# Predictors are identified by name: forge() matches new data by name, and
# component names join predictor names with ":".
check_predictor_names <- function(names, call = rlang::caller_env()) {
  missing <- which(is.na(names) | names == "")
  if (length(missing) > 0) {
    cli::cli_abort(
      c(
        "Every predictor must have a name.",
        "x" = "Column{?s} {as.character(missing)} {?has/have} no name.",
        "i" = "Set column names, e.g. with {.fn colnames}."
      ),
      call = call
    )
  }
  duplicated_names <- unique(names[duplicated(names)])
  if (length(duplicated_names) > 0) {
    cli::cli_abort(
      c(
        "Predictor names must be unique.",
        "x" = "Duplicated: {.val {duplicated_names}}."
      ),
      call = call
    )
  }
  with_colon <- names[grepl(":", names, fixed = TRUE)]
  if (length(with_colon) > 0) {
    cli::cli_abort(
      c(
        "Predictor names must not contain {.code :}.",
        "x" = "Found {.val {with_colon}}.",
        "i" = "Component names join predictor names with {.code :}, as in {.code x1:x2}. Rename {cli::qty(length(with_colon))}the column{?s}."
      ),
      call = call
    )
  }
  invisible(names)
}

# Interaction terms in a formula (x1:x2, x1 * x2, (...)^2) would add product
# columns as predictors, which the components would then mix with the forest's own
# interactions.
check_no_interaction_terms <- function(formula, data, call = rlang::caller_env()) {
  terms <- stats::terms(formula, data = data)
  interactions <- attr(terms, "term.labels")[attr(terms, "order") > 1]
  if (length(interactions) > 0) {
    cli::cli_abort(
      c(
        "The formula must not contain interaction terms.",
        "x" = "Found {.code {interactions}}.",
        "i" = "Such a term becomes a predictor named {.code x1:x2}, and {.code :} is how {.fn rpf} names components.",
        "i" = "{.fn rpf} finds interactions itself, up to {.arg max_interaction}: write {.code y ~ x1 + x2}.",
        "i" = "To fit a product as a predictor of its own, write {.code I(x1 * x2)}, or add it to the data as a column with its own name."
      ),
      call = call
    )
  }
  invisible(formula)
}

#' Order factor levels by response
#'
#' Regression: Order by mean(y)
#' Binary classification: Order by mean(y==1)
#' Multiclass: Order by first principal component of the weighted covariance
#' matrix of the contingency table
#'
#' See https://doi.org/10.7717/peerj.6339 for details.
#'
#' @param x `[factor]`: Variable to order.
#' @param y `[numeric | factor]`: Response.
#'
#' @return Re-ordered (ordered) factor
#' @noRd
order_factor_by_response <- function(x, y) {
  if (is.numeric(y)) {
    # Regression: Order by mean(y)
    means <- sapply(levels(x), function(lev) {
      mean(y[x == lev])
    })
  } else if (is.factor(y)) {
    if (nlevels(y) <= 2) {
      # Binary classification: Order by mean(y==1)
      means <- sapply(levels(x), function(lev) {
        mean((as.numeric(y) - 1)[x == lev])
      })
    } else {
      # Multiclass: Order by first principal component of the weighted
      # covariance matrix of the contingency table
      means <- pca_order(x, y)
    }
  } else {
    cli::cli_abort("Ordering of factor columns only implemented for regression and classification outcomes.")
  }

  levels_ordered <- as.character(levels(x)[order(means)])
  factor(x, levels = levels_ordered, ordered = TRUE, exclude = NULL)
}

#' PCA-ordering of factors
#'
#' Order factor levels by first principal component of the weighted covariance
#' matrix of the contingency table
#'
#' @param x `[factor]`: Variable to order.
#' @param y `[factor]`: Response.
#'
#' @return Order of factor levels
#'
#' @references
#' Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999)
#' Partitioning Nominal Attributes in Decision Trees.
#' Data Min Knowl Discov 3:197. \doi{10.1023/A:1009869804967}.
#' @noRd
pca_order <- function(x, y) {
  if (nlevels(x) < 2) {
    return(seq(1, nlevels(x)))
  }

  # Create contingency table of the nominal outcome with the nominal covariate
  N <- table(x, droplevels(y))

  # PCA of weighted covariance matrix of class probabilities
  P <- N / rowSums(N)
  S <- stats::cov.wt(P, wt = rowSums(N))$cov
  pc1 <- stats::prcomp(S, rank. = 1)$rotation
  score <- P %*% pc1

  # Return ordered factor levels
  order(score)
}

# Which factor predictors keep their level order: `keep_level_order` is a flag
# (TRUE: all factors) or names of factor predictors. Returns the names, in the
# order of the predictors.
resolve_keep_level_order <- function(keep_level_order, predictors, factor_cols, call = rlang::caller_env()) {
  if (is.logical(keep_level_order)) {
    if (!checkmate::test_flag(keep_level_order)) {
      cli::cli_abort(
        "{.arg keep_level_order} must be {.code TRUE}, {.code FALSE} or names of factor predictors.",
        call = call
      )
    }
    return(if (keep_level_order) factor_cols else character(0))
  }
  if (!checkmate::test_character(keep_level_order, any.missing = FALSE)) {
    cli::cli_abort(
      "{.arg keep_level_order} must be {.code TRUE}, {.code FALSE} or names of factor predictors.",
      call = call
    )
  }
  unknown <- setdiff(keep_level_order, names(predictors))
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        "{.arg keep_level_order} must name factor predictors.",
        "x" = "Not a predictor: {.val {unknown}}."
      ),
      call = call
    )
  }
  not_factor <- setdiff(keep_level_order, factor_cols)
  if (length(not_factor) > 0) {
    cli::cli_abort(
      c(
        "{.arg keep_level_order} must name factor predictors.",
        "x" = "Not a factor: {.val {not_factor}}."
      ),
      call = call
    )
  }
  intersect(factor_cols, keep_level_order)
}

# Sort factor predictors by outcome, except those in `keep_level_order`, which
# keep the order of their levels, and re-encode as integer.
# Save the factor levels in that order for the prediction step.
# Used in rpf_bridge()
preprocess_predictors_fit <- function(processed, keep_level_order = FALSE, call = rlang::caller_env()) {
  predictors <- as.data.table(processed$predictors)

  # Convert logicals to binary integers
  logical_cols <- names(which(sapply(predictors, is.logical)))
  if (length(logical_cols) > 0) {
    predictors[, (logical_cols) := lapply(.SD, as.integer), .SDcols = logical_cols]
  }

  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }

  # Factor predictors: Order by response (see https://doi.org/10.7717/peerj.6339),
  # unless they keep the order of their levels
  factor_cols <- names(which(sapply(predictors, is.factor)))
  kept_cols <- resolve_keep_level_order(keep_level_order, predictors, factor_cols, call = call)
  sorted_cols <- setdiff(factor_cols, kept_cols)
  if (length(sorted_cols) > 0) {
    predictors[,
      (sorted_cols) := lapply(
        .SD,
        order_factor_by_response,
        y = processed$outcomes[[1]]
      ),
      .SDcols = sorted_cols
    ]
  }

  # Save re-ordered factor levels
  factor_levels <- hardhat::get_levels(predictors)

  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }

  list(
    factor_levels = factor_levels,
    keep_level_order = kept_cols,
    predictors_matrix = as.matrix(predictors)
  )
}

#' Preprocess predictors for prediction
#'
#' Convert logical and character columns to appropriate types, re-order factor
#' levels to match the ordering learned during fitting (stored in
#' `object$factor_levels`), re-encode factor columns as integers, and return a
#' numeric matrix suitable for the underlying C++ prediction routines.
#'
#' @param object `[rpf]`: A fitted [`rpf`] model.
#' @param predictors `[data.frame | matrix]`: Predictor values to preprocess.
#'
#' @return A numeric matrix with the same number of rows as `predictors`.
#' @noRd
preprocess_predictors_predict <- function(object, predictors) {
  predictors <- as.data.table(predictors)

  # Convert logicals to binary integers
  logical_cols <- names(which(sapply(predictors, is.logical)))
  if (length(logical_cols) > 0) {
    predictors[, (logical_cols) := lapply(.SD, as.integer), .SDcols = logical_cols]
  }

  # Convert characters to factors
  char_cols <- names(which(sapply(predictors, is.character)))
  if (length(char_cols) > 0) {
    predictors[, (char_cols) := lapply(.SD, factor), .SDcols = char_cols]
  }

  # Re-order factor levels according to saved order
  factor_cols <- names(object$factor_levels)
  if (length(factor_cols) > 0) {
    predictors[,
      (factor_cols) := Map(
        factor,
        .SD,
        object$factor_levels,
        ordered = TRUE
      ),
      .SDcols = factor_cols
    ]
  }

  # Convert factors to integer and data to matrix
  if (length(factor_cols) > 0) {
    predictors[, (factor_cols) := lapply(.SD, as.integer), .SDcols = factor_cols]
  }

  as.matrix(predictors)
}

# Check outcome to detect mode (classif/regr), return suitable outcome
# Used in rpf_impl()
# Loss is need to transform 1/0 to 1/-1 for exponential
#' @importFrom stats model.matrix
preprocess_outcome <- function(processed, loss, call = rlang::caller_env()) {
  outcomes <- processed$outcomes[[1]]

  # Task type detection: Could be more concise
  is_binary <- length(unique(outcomes)) == 2
  is_integerish <- checkmate::test_integerish(outcomes, any.missing = FALSE)
  is_factor <- checkmate::test_factor(outcomes, any.missing = FALSE)
  is_numeric <- checkmate::test_numeric(outcomes, any.missing = FALSE)

  if (is_binary & is_integerish) {
    cli::cli_warn(c(
      "The outcome is a binary integer, assuming a regression task.",
      "i" = "Recode it to a factor for classification."
    ))
  }

  if (is_factor) {
    mode <- "classification"
    loss <- loss %||% "exponential"

    if (is_binary) {
      # Binary case: Convert to 0, 1 integer
      outcomes <- as.integer(outcomes) - 1L
      # rpf_impl expects Y to be a matrix
      outcomes <- as.matrix(outcomes, ncol = 1)
    } else {
      # Multiclass
      # One-hot encoding
      outcomes_mat <- stats::model.matrix(~ -1 + outcomes)
      colnames(outcomes_mat) <- levels(outcomes)

      if (loss == "logit") {
        # The C++ logit loss is a reference-class multinomial formulation:
        # it models K-1 logits and treats 1 - rowSums(Y) as the implicit
        # reference class (first factor level), so it needs K-1 indicator
        # columns. Full one-hot pins the reference-class probability to
        # zero and no split ever improves the loss (#40).
        outcomes_mat <- outcomes_mat[, -1, drop = FALSE]
      }

      outcomes <- outcomes_mat
    }

    # Exponential expects outcome in 1/-1
    if (loss == "exponential") {
      outcomes[outcomes == 0] <- -1
    }
  } else if (is_numeric) {
    mode <- "regression"
    loss <- loss %||% "L2"
    # rpf_impl expects Y to be a matrix
    outcomes <- as.matrix(outcomes, ncol = 1)
  } else {
    # mode <- "unsupported"; check_outcome() has rejected these already
    abort_outcome_type(outcomes, call = call)
  }

  list(
    outcomes = outcomes,
    mode = mode,
    loss = loss
  )
}

# Softmax for multiclass prediction using centered logsumexp for numerical stability
# should be identical to sigmoid for scalar input
# Works on matrix input where one row is assumed to be one vector of predictions
# for a single observations
softmax <- function(x) {
  rowmax <- apply(x, 1, max)
  lse <- rowmax + log(rowSums(exp(x - rowmax)))
  exp(x - lse)
}

# The fit's predictor columns, in the fit's order (the columns of the C++ forest).
# Not names(blueprint$ptypes$predictors): those are the raw columns a formula reads, while the fit has a column
# per term (log(x1), I(x1^2); y ~ . - x3 drops x3), and hardhat::mold() puts factors last.
fit_predictor_names <- function(object) {
  ptype <- hardhat::forge(object$blueprint$ptypes$predictors, object$blueprint)$predictors
  colnames(preprocess_predictors_predict(object, ptype))
}

# Outcome levels represented by columns of the C++ prediction matrix.
# Multiclass logit uses reference-class encoding: the first level has no column.
model_outcome_levels <- function(object) {
  outcome_levels <- levels(object$blueprint$ptypes$outcomes[[1]])
  if (object$params$loss == "logit" && length(outcome_levels) > 2) {
    outcome_levels <- outcome_levels[-1]
  }
  outcome_levels
}

#' Get remainders where max_interaction requested in `predict_components` is smaller than
#' `max_interaction` set in `rpf`.
#' This is somewhat cumbersome unfortunately, and presumably will have to be partially
#' repeated for other methods in `glex`
#' @noRd
#' @keywords internal
#' @param m `[data.table]`: Components as calculated in `predict_components`
#' @param levels `[character]`: Outcome levels as stored in `rpf$blueprint$ptypes$outcomes`.
#' @param pred `[data.frame]`: Regular model predictions as returned by `predict.rpf`.
#' @param intercept `[numeric]`: Intercept as stored in output of `predict_components`.
calc_remainders_multiclass <- function(m, levels, pred, intercept) {
  # data.table NSE warnings
  term <- remainder <- m_sum <- NULL

  m <- data.table::as.data.table(m)
  pred <- data.table::as.data.table(pred)
  intercept <- data.table::as.data.table(intercept)

  m[, ".id" := .I]
  pred[, ".id" := .I]
  intercept[, ".id" := .I]

  split_names <- function(mn, split_string = "__class:", target_index = 2) {
    vapply(
      mn,
      function(x) {
        unlist(strsplit(x, split = split_string, fixed = TRUE))[[target_index]]
      },
      character(1),
      USE.NAMES = FALSE
    )
  }

  m_long <- data.table::melt(m, id.vars = ".id", value.name = "m", variable.name = "term", variable.factor = FALSE)

  m_long[, class := split_names(term, split_string = "__class:", target_index = 2)]
  m_sums <- m_long[, list(m_sum = sum(m)), by = c(".id", "class")]

  # long format predictions for calculations
  pred_long <- data.table::melt(
    pred,
    id.vars = ".id",
    value.name = "pred",
    variable.name = "class",
    variable.factor = FALSE
  )
  pred_long[, class := split_names(class, split_string = ".pred_", target_index = 2)]

  # Same game for intercept
  intercept_long <- data.table::melt(
    intercept,
    id.vars = ".id",
    value.name = "intercept",
    variable.name = "class",
    variable.factor = FALSE
  )
  intercept_long[, class := split_names(class, split_string = "__class:", target_index = 2)]

  # merge pred with intercept
  pred_long <- pred_long[intercept_long, on = c(".id", "class")]
  # merge predictions with sum(m)
  merged <- pred_long[m_sums, on = c(".id", "class")]
  # remainder is prediction - sum(m) - intercept, can be positive or negative
  merged[, remainder := pred - m_sum - intercept]
  # Return one column per outcome class
  merged_wide <- data.table::dcast(merged, .id ~ class, value.var = "remainder")
  # Remove intermediate id column
  merged_wide[, ".id" := NULL]
  # sort columns to original outcome level order, consistent with predict() order
  data.table::setcolorder(merged_wide, levels)

  merged_wide
}

# The predictors' names as the forest's columns hold them (the order of the sets' indices: hardhat::mold() puts the
# factors last, so it can differ from the ptype's)
prune_fit_names <- function(fit) {
  ptype <- hardhat::forge(fit$blueprint$ptypes$predictors, fit$blueprint)$predictors
  colnames(preprocess_predictors_predict(fit, ptype))
}
