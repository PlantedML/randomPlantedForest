#' Prune the components of a Random Planted Forest
#'
#' Drops the small components of a regression forest's functional decomposition and, optionally, refits one
#' coefficient per kept component. The rule and the threshold are chosen by K-fold cross-validation on the training
#' rows.
#'
#' The forest is decomposed into its purified components `f_S`, one per set `S` of predictors (see [purify()] and
#' [predict_components()]), at all orders up to `max_interaction`, and its intercept `f_0`. The size of a component is
#' its share of the variance of the prediction `f` on the training rows, `share_S = Var(f_S) / Var(f)`. The shares are
#' descriptive sizes, not a partition: the components are not orthogonal on the data, so they need not add up to 1. A
#' component is kept if its share is at least `tau` (and it is not identically zero on the training rows). Two rules
#' are available:
#'
#' - `"D"` (drop): the prediction is `f_0` plus the kept components. `rule = "D"` with `tau = 0` is the unpruned
#'   forest.
#' - `"R"` (drop and refit): the prediction is `alpha + sum(beta_S * f_S)` over the kept components, with the
#'   intercept `alpha` and the coefficients `beta_S` by least squares of `y` on the kept components over the training
#'   rows (the minimum-norm solution if the components are collinear). The components themselves are not changed,
#'   only their scale.
#'
#' **Cross-validation.** If more than one pair of `rule` and `tau` is given (the default: both rules and the grid
#' 0, 0.01%, 0.03%, 0.1%, 0.3%, 1%, 3%, 10%), the pair is chosen by K-fold cross-validation on the training rows. For
#' each fold, the forest is refitted with the parameters of `object` on the other folds (after `set.seed(seed)`, the
#' same for every fold), its components, shares and refit are computed on those rows, and every candidate is scored
#' on the held-out fold. The CV RMSE is `sqrt(sum(SSE_k) / n)`, summing the squared errors over the K folds. With
#' `select = "min"`, the pair with the lowest CV RMSE is chosen (ties, up to a relative `1e-12`: the larger `tau`,
#' then `"D"`). With `select = "1se"`, the largest `tau` whose CV MSE is within one standard error of the minimum is
#' chosen (the SE is the standard deviation of the per-fold MSEs of the minimum, divided by `sqrt(K)`), and at that
#' `tau` the rule with the lower CV RMSE.
#'
#' The chosen pair is then applied to `object` itself (not refitted): its components and shares on `data`, and for
#' `"R"` the least squares on `data`. `data` and `y` must therefore be the rows `object` was trained on, in the order
#' of the fit (checked).
#'
#' **Purification.** `object` is purified in place at all orders, as [purify()] does, if it is not purified yet. This
#' does not change its predictions. A forest purified with a `maxp_interaction` below its order is purified again at
#' all orders, with a message: the pruning needs every component.
#'
#' **Random numbers.** The folds are drawn after `set.seed(seed)`, and the fold fits after it too. Afterwards the
#' random number stream is restored to its state on entry (after drawing `seed` if it is `NULL`), so the call does not
#' reset the caller's stream.
#'
#' **Cost.** The CV fits and purifies the forest K more times. The components are evaluated on the training rows:
#' only those the purified forest stores (at most `sum(choose(p, seq_len(max_interaction)))`, usually far fewer: the
#' subsets of the trees' predictor sets), in one pass over the rows. The memory is about `n` doubles per stored
#' component. `nthreads` speeds up the fits and the purification.
#'
#' Only regression forests are supported.
#'
#' @param object `[rpf]`: A fitted regression [`rpf`] model.
#' @param data `[data.frame | matrix]`: The training data of `object`, in the order of the fit: the predictors as
#'   passed to [rpf()] as `x`, or the `data` of the formula or recipe interface. For a formula or recipe fit, `data`
#'   may contain the outcome.
#' @param y `[numeric | NULL: NULL]`: The outcome of the training rows. Required unless `data` contains the outcome
#'   of a formula or recipe fit.
#' @param ... Reserved for future expansion. Unknown arguments are an error.
#' @param rule `[character: c("D", "R")]`: The rules to consider: `"D"` drops, `"R"` drops and refits one coefficient
#'   per kept component.
#' @param tau `[numeric: c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1)]`: The thresholds to consider, as shares of
#'   the variance of the prediction (`0.01` = 1%), each at least 0 and below 1.
#' @param select `[character(1): "min"]`: How the CV chooses: `"min"` (the lowest CV RMSE) or `"1se"` (the largest
#'   `tau` within one standard error of it).
#' @param folds `[integer(1): 5]`: The number of CV folds, used if `fold_id` is `NULL`.
#' @param fold_id `[vector | NULL: NULL]`: The fold of every training row: a vector of length `nrow(data)` with at
#'   least two distinct values. If `NULL`, the rows are assigned to `folds` folds at random (after `set.seed(seed)`).
#' @param seed `[integer(1) | NULL: NULL]`: The seed set before the random fold assignment and before every fold's
#'   refit. If `NULL`, one is drawn from the current random number stream and stored in the result, so that the CV
#'   can be reproduced. The caller's random number state is restored afterwards (see Details).
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads for the fold fits and for purifying the forests
#'   (`object`, if needed, and the fold fits). `NULL` uses the `nthreads` `object` was fitted with. The result does
#'   not depend on it.
#' @param verbose `[logical(1): FALSE]`: Report the progress of the CV.
#'
#' @return An object of class `rpf_pruned`, a list with elements:
#' - `fit`: `object` (purified).
#' - `rule`, `tau`: the chosen rule and threshold.
#' - `intercept`: `alpha` (the forest's intercept `f_0` for `"D"`).
#' - `beta`: named vector of the coefficients of the kept components (all 1 for `"D"`).
#' - `sets`: the predictor sets of the kept components, as column indices.
#' - `components`: a `data.frame` of every nonzero component of `object` on the training rows: `component` (its
#'   name, as in [predict_components()]), `order`, `share`, `kept`, `beta` (`NA` if dropped).
#' - `cv`: `NULL` if no CV was run; otherwise a `data.frame` with one row per candidate: `rule`, `tau`, `cv_rmse`,
#'   `se_mse` (SE of the CV MSE), `n_kept` (mean over the folds), `min` and `one_se` (the two choices).
#' - `select`, `fold_id`, `seed`: as used (`fold_id` and `seed` are `NULL` if no CV was run).
#'
#' Use [`predict()`][predict.rpf_pruned] and [predict_components()] on it, and [rpf_marshal()] to save it.
#'
#' @seealso [predict.rpf_pruned()], [predict_components()], [purify()].
#' @export
#'
#' @examples
#' # A main effect and a pure interaction
#' set.seed(1)
#' x <- matrix(runif(600 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(600, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, ntrees = 20, splits = 50)
#'
#' # Rule and threshold chosen by 5-fold CV
#' pruned <- prune_components(fit, x, y, seed = 1)
#' pruned
#' pruned$cv
#'
#' # A given rule and threshold: no CV
#' prune_components(fit, x, y, rule = "R", tau = 0.01)
#'
#' predict(pruned, x[1:5, ])
prune_components <- function(
  object,
  data,
  y = NULL,
  ...,
  rule = c("D", "R"),
  tau = c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1),
  select = c("min", "1se"),
  folds = 5L,
  fold_id = NULL,
  seed = NULL,
  nthreads = NULL,
  verbose = FALSE
) {
  rlang::check_dots_empty()
  checkmate::assert_class(object, "rpf")
  check_rpf_alive(object)
  if (object$mode != "regression") {
    cli::cli_abort(c(
      "{.fn prune_components} supports regression forests only.",
      "x" = "{.arg object} is a {object$mode} forest."
    ))
  }
  if (missing(data) || is.null(data)) {
    cli::cli_abort("{.arg data} is required: the training rows of {.arg object}, as passed to {.fn rpf}.")
  }
  checkmate::assert_character(rule, min.len = 1, any.missing = FALSE)
  checkmate::assert_subset(rule, c("D", "R"))
  rule <- unique(rule)
  checkmate::assert_numeric(tau, lower = 0, any.missing = FALSE, min.len = 1, unique = TRUE)
  if (any(tau >= 1)) {
    cli::cli_abort("{.arg tau} is a share of the variance of the prediction and must be below 1.")
  }
  select <- rlang::arg_match(select)
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  nthreads <- as.integer(nthreads %||% object$params$nthreads)
  checkmate::assert_flag(verbose)

  tr <- prune_training_data(object, data, y)
  n <- length(tr$y)
  # the shares and R's coefficients are those of the training rows
  prune_check_training_rows(object, tr)
  prune_check_purified(object, nthreads)
  cand <- expand.grid(tau = sort(tau), rule = rule, stringsAsFactors = FALSE)[, c("rule", "tau")]

  cv <- NULL
  if (nrow(cand) > 1) {
    if (is.null(seed)) {
      seed <- sample.int(.Machine$integer.max, 1L)
    }
    checkmate::assert_int(seed)
    seed <- as.integer(seed)
    # the folds and the fold fits set the seed: afterwards the caller's random number stream is restored
    rng <- prune_save_rng()
    on.exit(prune_restore_rng(rng), add = TRUE)
    if (is.null(fold_id)) {
      checkmate::assert_int(folds, lower = 2, upper = n)
      set.seed(seed)
      fold_id <- sample(rep_len(seq_len(folds), n))
    } else {
      checkmate::assert_atomic_vector(fold_id, any.missing = FALSE, len = n)
      if (length(unique(fold_id)) < 2) {
        cli::cli_abort("{.arg fold_id} must have at least two distinct folds.")
      }
    }
    cv <- prune_cv(object, tr, cand, fold_id, seed, nthreads, verbose)
    choice <- cv[cv[[c(min = "min", "1se" = "one_se")[[select]]]], ]
    rule_sel <- choice$rule
    tau_sel <- choice$tau
  } else {
    rule_sel <- cand$rule
    tau_sel <- cand$tau
    fold_id <- NULL
    seed <- NULL
  }

  # The chosen rule and tau on the forest itself, on all training rows
  dec <- prune_decompose(object, tr$predictors, rep(TRUE, n), nthreads = nthreads)
  sets <- dec$sets[dec$nonzero]
  rf <- prune_fit_rule(rule_sel, tau_sel, dec, tr$y)

  components <- data.frame(
    component = names(dec$share),
    order = lengths(sets),
    share = unname(dec$share),
    kept = names(dec$share) %in% names(rf$beta),
    beta = unname(rf$beta[names(dec$share)]),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      fit = object,
      rule = rule_sel,
      tau = tau_sel,
      intercept = rf$alpha,
      beta = rf$beta,
      sets = sets[components$kept],
      components = components,
      cv = cv,
      select = select,
      fold_id = fold_id,
      seed = seed
    ),
    class = "rpf_pruned"
  )
}

# The training rows: forged predictors (a tibble, as rpf() gets them from hardhat::mold()) and the outcome
prune_training_data <- function(object, data, y, call = rlang::caller_env()) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data.frame or matrix, not {.obj_type_friendly {data}}.", call = call)
  }
  bp <- object$blueprint
  if (is.null(y)) {
    processed <- tryCatch(
      hardhat::forge(data, bp, outcomes = TRUE),
      error = function(e) {
        cli::cli_abort(
          c(
            "{.arg y} is missing and the outcome could not be taken from {.arg data}.",
            "i" = "Pass the outcome as {.arg y}."
          ),
          call = call
        )
      }
    )
    y <- processed$outcomes[[1]]
  } else {
    processed <- hardhat::forge(data, bp)
  }
  if (!is.numeric(y) || anyNA(y)) {
    cli::cli_abort("{.arg y} must be a numeric outcome without missing values.", call = call)
  }
  y <- as.numeric(y)
  n <- nrow(processed$predictors)
  if (length(y) != n) {
    cli::cli_abort("{.arg y} has {length(y)} value{?s} but {.arg data} has {n} row{?s}.", call = call)
  }
  if (n < 2) {
    cli::cli_abort("{.arg data} must have at least two rows.", call = call)
  }
  check_no_missing(processed$predictors, "{.arg data}", call = call)
  list(predictors = processed$predictors, y = y)
}

# `data` and `y` must be the training rows of `object`, in the order of the fit (bitwise, as the forest stores them)
prune_check_training_rows <- function(object, tr, call = rlang::caller_env()) {
  train <- object$fit$get_data()
  if (length(train$X) == 0) {
    cli::cli_abort(
      c(
        "{.arg object} has no training data, so {.arg data} cannot be checked to be its training rows.",
        "i" = "It was restored by {.fn rpf_unmarshal} from a blob made without {.code include_data = TRUE}."
      ),
      call = call
    )
  }
  X <- preprocess_predictors_predict(object, tr$predictors)
  same <- identical(dim(X), dim(train$X)) && all(X == train$X) && all(tr$y == train$Y)
  if (!same) {
    cli::cli_abort(
      "{.arg data} and {.arg y} must be the training rows of {.arg object}, in the order of the fit.",
      call = call
    )
  }
  invisible()
}

# The random number state of the session (NULL if there is none yet), and its restoration
prune_save_rng <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}
prune_restore_rng <- function(rng) {
  if (!is.null(rng)) {
    assign(".Random.seed", rng, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
}

# Purify `object` (mode 2, all orders) unless it is purified at all its orders. A purification below the forest's
# order (purify(maxp_interaction = k)) leaves the higher components zero, which the pruning would silently drop (and
# the fold fits would not): it is purified again at all orders. predict() does not change.
prune_check_purified <- function(object, nthreads) {
  if (!is_purified(object)) {
    object$fit$purify_threads(0L, nthreads, 2L)
  } else if (object$fit$is_purify_capped()) {
    cli::cli_inform(c(
      "i" = "{.arg object} was purified with a {.arg maxp_interaction} below its order; purifying it at all orders.",
      " " = "{.fn prune_components} prunes the whole forest. The forest's predictions do not change."
    ))
    object$fit$purify_threads(0L, nthreads, 2L)
  }
  invisible()
}

# The predictor sets of the components the purified forest stores (as column indices), up to its order, ordered by
# order and then lexicographically (the order of utils::combn() over all sets)
prune_stored_sets <- function(fit) {
  sets <- lapply(fit$fit$purified_components(), as.integer)
  sets <- sets[lengths(sets) >= 1 & lengths(sets) <= fit$params$max_interaction]
  if (length(sets) == 0) {
    return(list())
  }
  width <- max(lengths(sets))
  key <- vapply(sets, function(s) c(length(s), s, rep(0L, width - length(s))), integer(width + 1))
  sets[do.call(order, as.data.frame(t(key)))]
}

# Names of predictor sets, as predict_components() names its columns
prune_set_names <- function(sets, feature_names) {
  vapply(sets, function(s) paste0(sort(feature_names[s]), collapse = ":"), character(1))
}

# The components of `fit` (purified first if needed) on the preprocessed predictor matrix X, one column per set, in
# one pass
prune_component_matrix <- function(fit, X, sets, nthreads = 0L) {
  if (!is_purified(fit)) {
    fit$fit$purify_threads(0L, nthreads, 2L)
  }
  out <- fit$fit$predict_components_matrix(X, sets, nthreads)
  dimnames(out) <- NULL
  out
}

# Decomposition of `fit` on its own training rows (logical `own`) and, optionally, other rows (logical `other`): the
# components, the intercept, the shares of the nonzero components (named) and which sets they are. Only the
# components the purified forest stores are evaluated (the others are zero), not all 2^p - 1 sets.
prune_decompose <- function(fit, predictors, own, other = NULL, nthreads = 0L) {
  X <- preprocess_predictors_predict(fit, predictors)
  if (!is_purified(fit)) {
    fit$fit$purify_threads(0L, nthreads, 2L)
  }
  sets <- prune_stored_sets(fit)
  rows <- c(which(own), if (!is.null(other)) which(other))
  M <- prune_component_matrix(fit, X[rows, , drop = FALSE], sets, nthreads)
  n_own <- sum(own)
  M_own <- M[seq_len(n_own), , drop = FALSE]
  intercept <- as.numeric(fit$fit$predict_matrix(X[1, 1, drop = FALSE], -1, 0L))[1]
  # the purified forest's prediction is its intercept plus all its components, which are all in M
  pred_own <- intercept + rowSums(M_own)
  var_pop <- function(v) mean((v - mean(v))^2)
  v_comp <- vapply(seq_len(ncol(M_own)), function(k) var_pop(M_own[, k]), numeric(1))
  nonzero <- v_comp > 0
  share <- v_comp[nonzero] / var_pop(pred_own)
  names(share) <- prune_set_names(sets[nonzero], colnames(X))
  M <- M[, nonzero, drop = FALSE]
  colnames(M) <- names(share)
  list(
    sets = sets,
    nonzero = nonzero,
    share = share,
    intercept = intercept,
    own = M[seq_len(n_own), , drop = FALSE],
    other = if (!is.null(other)) M[-seq_len(n_own), , drop = FALSE]
  )
}

# Minimum-norm least squares (as LAPACK's gelsd with the default cutoff eps * max(dim) * largest singular value)
prune_lstsq <- function(X, y) {
  s <- svd(X)
  tol <- max(dim(X)) * .Machine$double.eps * s$d[1]
  d <- ifelse(s$d > tol, 1 / s$d, 0)
  as.numeric(s$v %*% (d * crossprod(s$u, y)))
}

# The kept components of rule and tau, and (alpha, beta) on the decomposition's own rows
prune_fit_rule <- function(rule, tau, dec, y_own) {
  kept <- names(dec$share)[dec$share >= tau]
  if (rule == "D") {
    return(list(alpha = dec$intercept, beta = stats::setNames(rep(1, length(kept)), kept)))
  }
  coef <- prune_lstsq(cbind(1, dec$own[, kept, drop = FALSE]), y_own)
  list(alpha = coef[1], beta = stats::setNames(coef[-1], kept))
}

prune_predict_dec <- function(rf, comp) {
  rf$alpha + as.numeric(comp[, names(rf$beta), drop = FALSE] %*% rf$beta)
}

# Refit `object`'s forest on the rows `rows` with all its parameters (as rpf() would fit them, after set.seed(seed)),
# on `nthreads` threads and without purification. The factor predictors are coded on those rows, as rpf() codes them.
prune_refit <- function(object, predictors, y, rows, seed, nthreads) {
  processed <- list(
    predictors = predictors[rows, , drop = FALSE],
    outcomes = data.frame(.outcome = y[rows])
  )
  pre <- preprocess_predictors_fit(processed)
  params <- object$params
  params$nthreads <- nthreads
  params$purify <- FALSE
  set.seed(seed)
  fit <- rpf_impl(
    Y = matrix(y[rows], ncol = 1),
    X = pre$predictors_matrix,
    mode = "regression",
    params = params
  )
  new_rpf(
    fit = fit,
    blueprint = object$blueprint,
    mode = "regression",
    factor_levels = pre$factor_levels,
    params = params
  )
}

# K-fold CV of every candidate (rule, tau); returns the CV table with the "min" and "1se" choices marked
prune_cv <- function(object, tr, cand, fold_id, seed, nthreads, verbose) {
  n <- length(tr$y)
  fids <- sort(unique(fold_id))
  sse <- matrix(NA_real_, nrow(cand), length(fids))
  n_kept <- matrix(NA_integer_, nrow(cand), length(fids))
  n_k <- integer(length(fids))
  for (k in seq_along(fids)) {
    ho <- fold_id == fids[k]
    t0 <- proc.time()[["elapsed"]]
    fit_k <- prune_refit(object, tr$predictors, tr$y, !ho, seed, nthreads)
    dec <- prune_decompose(fit_k, tr$predictors, !ho, other = ho, nthreads = nthreads)
    y_own <- tr$y[!ho]
    y_ho <- tr$y[ho]
    n_k[k] <- sum(ho)
    cache <- list()
    for (i in seq_len(nrow(cand))) {
      kept <- names(dec$share)[dec$share >= cand$tau[i]]
      key <- paste(cand$rule[i], paste(kept, collapse = ","))
      if (is.null(cache[[key]])) {
        rf <- prune_fit_rule(cand$rule[i], cand$tau[i], dec, y_own)
        cache[[key]] <- sum((prune_predict_dec(rf, dec$other) - y_ho)^2)
      }
      sse[i, k] <- cache[[key]]
      n_kept[i, k] <- length(kept)
    }
    if (verbose) {
      secs <- sprintf("%.1f", proc.time()[["elapsed"]] - t0)
      n_nonzero <- length(dec$share)
      cli::cli_inform(
        "Fold {fids[k]}: {n_k[k]} row{?s} held out, {n_nonzero} nonzero component{?s}, {secs} s."
      )
    }
    rm(fit_k, dec)
  }

  cv_rmse <- sqrt(rowSums(sse) / n)
  mse_k <- sweep(sse, 2, n_k, "/")
  se_mse <- apply(mse_k, 1, stats::sd) / sqrt(length(fids))
  ch <- prune_choose(cv_rmse, se_mse, cand)

  data.frame(
    rule = cand$rule,
    tau = cand$tau,
    cv_rmse = cv_rmse,
    se_mse = se_mse,
    n_kept = rowMeans(n_kept),
    min = ch$min,
    one_se = ch$one_se,
    stringsAsFactors = FALSE
  )
}

# The two choices among the candidates from their RMSE and the SE of their MSE (logical vectors)
prune_choose <- function(rmse, se_mse, cand) {
  rule_rank <- match(cand$rule, c("D", "R"))
  # min: the lowest RMSE; ties (relative 1e-12): the larger tau, then D
  ties <- which(abs(rmse / min(rmse) - 1) <= 1e-12)
  best <- ties[order(-cand$tau[ties], rule_rank[ties])][1]
  # 1se: the largest tau with MSE within one SE of the minimum; at that tau the lower RMSE, then D
  ok <- which(rmse^2 <= (rmse[best]^2 + se_mse[best]) * (1 + 1e-12))
  ok <- ok[cand$tau[ok] == max(cand$tau[ok])]
  one <- ok[order(rmse[ok], rule_rank[ok])][1]
  list(min = seq_along(rmse) == best, one_se = seq_along(rmse) == one)
}

# The kept components (scaled by beta) on new data, as a matrix with one named column per component
prune_kept_matrix <- function(object, new_data, nthreads = 0L, call = rlang::caller_env()) {
  fit <- object$fit
  check_rpf_alive(fit, call = call)
  processed <- hardhat::forge(new_data, fit$blueprint)
  check_no_missing(processed$predictors, "{.arg new_data}", call = call)
  X <- preprocess_predictors_predict(fit, processed$predictors)
  M <- prune_component_matrix(fit, X, object$sets, nthreads)
  M <- sweep(M, 2, object$beta, "*")
  colnames(M) <- names(object$beta)
  M
}

#' Predictions of a pruned Random Planted Forest
#'
#' The prediction of a pruned forest is `alpha + sum(beta_S * f_S(x))` over the kept components `f_S` of the forest,
#' with the intercept `alpha` and the coefficients `beta_S` (see [prune_components()]).
#'
#' @param object `[rpf_pruned]`: A pruned forest from [prune_components()].
#' @param new_data `[data.frame | matrix]`: New observations to predict.
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the `nthreads` the forest was fitted
#'   with.
#' @param ... Not currently used, but required for extensibility. Unknown arguments are an error.
#'
#' @return A [`tbl`][tibble::tibble] with column `.pred` with the same number of rows as `new_data`.
#' @seealso [prune_components()], [predict_components()].
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, ntrees = 10, splits = 30)
#' pruned <- prune_components(fit, x, y, rule = "R", tau = 0.01)
#' predict(pruned, x[1:5, ])
predict.rpf_pruned <- function(object, new_data, nthreads = NULL, ...) {
  rlang::check_dots_empty()
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  # 0 tells C++ to use the forest's own nthreads
  M <- prune_kept_matrix(object, new_data, as.integer(nthreads %||% 0L))
  out <- hardhat::spruce_numeric(object$intercept + rowSums(M))
  hardhat::validate_prediction_size(out, new_data)
  out
}

#' @rdname predict_components
#' @export
predict_components.rpf_pruned <- function(object, new_data, max_interaction = NULL, predictors = NULL, ...) {
  rlang::check_dots_empty()
  all_predictors <- names(object$fit$blueprint$ptypes$predictors)
  if (is.null(predictors)) {
    predictors <- all_predictors
  } else {
    checkmate::assert_subset(predictors, choices = all_predictors)
  }
  max_order <- object$fit$params$max_interaction
  if (is.null(max_interaction)) {
    max_interaction <- max_order
  } else {
    checkmate::assert_int(max_interaction, lower = 1, upper = max_order)
  }

  M <- prune_kept_matrix(object, new_data)
  # the kept components that use only `predictors` and have at most `max_interaction` predictors (the sets index the
  # fit's columns, in which hardhat::mold() puts the factors last)
  fit_names <- prune_fit_names(object$fit)
  feats <- lapply(object$sets, function(s) fit_names[s])
  use <- vapply(feats, function(f) length(f) <= max_interaction && all(f %in% predictors), logical(1))

  ret <- list(
    m = data.table::as.data.table(M[, use, drop = FALSE]),
    intercept = object$intercept,
    # a copy: setDT() would fail on a matrix and turn the caller's data.frame into a data.table
    x = data.table::as.data.table(new_data)[, predictors, with = FALSE]
  )
  if (!all(use)) {
    ret$remainder <- rowSums(M[, !use, drop = FALSE])
  }
  class(ret) <- c("glex", "rpf_components", class(ret))
  ret
}

#' Print a pruned Random Planted Forest
#'
#' Shows the rule and threshold, how they were chosen, what was kept and dropped, and the coefficients.
#'
#' @param x `[rpf_pruned]`: A pruned forest from [prune_components()].
#' @param n_show `[integer(1): 20]`: The largest number of kept components to list.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly: `x`.
#' @seealso [prune_components()].
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, ntrees = 10, splits = 30)
#' prune_components(fit, x, y, rule = "R", tau = 0.01)
print.rpf_pruned <- function(x, n_show = 20, ...) {
  cat(format(x, n_show = n_show, ...), sep = "\n")
  invisible(x)
}

#' @export
format.rpf_pruned <- function(x, n_show = 20, ...) {
  checkmate::assert_int(n_show, lower = 0)
  cmp <- x$components
  rule_name <- c(D = "D (drop)", R = "R (drop + refit one coefficient per component)")[[x$rule]]
  how <- if (is.null(x$cv)) {
    "Given, not chosen by CV."
  } else {
    sprintf(
      "Chosen by %d-fold CV (%s, seed %d) among %d candidates.",
      length(unique(x$fold_id)),
      if (x$select == "min") "lowest CV RMSE" else "1-SE rule",
      x$seed,
      nrow(x$cv)
    )
  }
  kept <- cmp[cmp$kept, ]
  dropped <- cmp[!cmp$kept, ]
  by_order <- function(d) {
    if (nrow(d) == 0) {
      return("0")
    }
    tab <- table(factor(d$order, levels = sort(unique(cmp$order))))
    paste0(nrow(d), " (by order ", paste(names(tab), tab, sep = ": ", collapse = ", "), ")")
  }
  lines <- c(
    sprintf("Rule %s, tau = %s of Var(prediction)", rule_name, format_share(x$tau)),
    how,
    "",
    sprintf("Nonzero components: %d", nrow(cmp)),
    paste("Kept:   ", by_order(kept)),
    paste("Dropped:", by_order(dropped))
  )
  if (nrow(dropped) > 0) {
    big <- dropped[which.max(dropped$share), ]
    lines <- c(
      lines,
      sprintf(
        "Largest dropped: %s (share %s); dropped shares sum to %s",
        big$component,
        format_share(big$share),
        format_share(sum(dropped$share))
      )
    )
  }
  lines <- c(lines, "", sprintf("Intercept: %s", format(x$intercept, digits = 6)))
  if (nrow(kept) > 0 && n_show > 0) {
    kept <- kept[order(-kept$share), ]
    show <- utils::head(kept, n_show)
    tab <- data.frame(
      component = show$component,
      share = vapply(show$share, format_share, ""),
      beta = format(show$beta, digits = 4)
    )
    lines <- c(
      lines,
      "Kept components (largest first), share on the training rows and coefficient:",
      utils::capture.output(print(tab, row.names = FALSE, right = FALSE))
    )
    if (nrow(kept) > n_show) {
      lines <- c(lines, sprintf("... and %d more", nrow(kept) - n_show))
    }
  }

  # format into lines so print() writes to stdout, not cli's message stream
  c(cli::cli_format_method(cli::cli_rule(left = "Pruned Random Planted Forest")), lines)
}

format_share <- function(s) {
  if (s == 0) {
    return("0")
  }
  paste0(trimws(formatC(signif(100 * s, 2), digits = 2, format = "fg")), "%")
}
