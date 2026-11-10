#' Prune the components of a Random Planted Forest
#'
#' Drops the small components of a regression forest's functional decomposition and, optionally, refits one
#' coefficient per kept component. The rule and the threshold are chosen out of bag by default, from the forest itself
#' without refits, or by K-fold cross-validation on the training rows.
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
#' **The choice.** If more than one pair of `rule` and `tau` is given (the default: both rules and the grid 0, 0.01%,
#' 0.03%, 0.1%, 0.3%, 1%, 3%, 10%), the pair is chosen by its out-of-bag (`select_by = "oob"`, the default) or
#' cross-validated (`"cv"`) RMSE. With `select = "min"`, the pair with the lowest RMSE is chosen (ties, up to a
#' relative `1e-12`: the larger `tau`, then `"D"`). With `select = "1se"`, the largest `tau` whose MSE is within one
#' standard error of the minimum is chosen (the SE is the standard deviation of the per-fold MSEs of the minimum,
#' divided by `sqrt(K)`), and at that `tau` the rule with the lower RMSE.
#'
#' **Out of bag** (`select_by = "oob"`). Every tree family of `object` is grown on a bootstrap sample of the training
#' rows, and the purification works family by family. So `object` is purified once more (as [purify()] does), and on
#' the way every family's components are evaluated at the training rows outside its bootstrap sample. A row's
#' out-of-bag (OOB) component is the mean over its OOB families (about 37% of them), and its OOB intercept likewise.
#' The purification uses the predictors of all training rows, not the outcome, so no outcome information of a row
#' enters its OOB components.
#' - The kept components of every `tau` are those of the final model (the shares of `object` on the training rows).
#' - `"D"` is scored by the OOB intercept plus the kept OOB components.
#' - For `"R"`, `alpha` and `beta_S` are fitted as for the final model, by least squares on the components of
#'   `object`, but only on the rows outside each of K folds (`folds`, `fold_id`, `seed`). They are scored with the OOB
#'   components of the fold's rows, so that no row scores coefficients fitted on it.
#' - An OOB mean averages `m_i` families (about 37% of `B = ntrees`), the final model all `B`, so the OOB prediction
#'   scatters more around its expectation. The OOB MSE is therefore corrected: per row it subtracts
#'   `(1 / m_i - 1 / B) * v_i`, with `v_i` the across-family variance of the candidate's prediction, estimated from
#'   the OOB means over the even and over the odd families (for a row whose OOB families are all even or all odd: the
#'   mean estimate of the other rows). Without the correction, the noise of the OOB components (amplified by
#'   `beta_S > 1`) makes `"R"` look worse than it is. The uncorrected values are reported as `oob_rmse_raw`.
#' - The OOB RMSE is taken over the rows with at least one OOB family, the SE for `"1se"` from the per-fold MSEs.
#' - With few families the OOB means average only one to a few families, and many rows have none (16% at 4
#'   families). In a simulation (p = 6, n = 800, 60 replicates) the corrected OOB error of the unpruned `"D"` was
#'   within 1.5% of the test error at 4, 10, 20 and 50 families, that of `"R"` within 4.4% except once (+11.5% at 10
#'   families with a 3-way interaction), and the OOB choice was as good as the CV's or better from 10 families on
#'   (within 0.4% in test error, or better by 7-10% where the CV's fold fits choose poorly). Below 10 families
#'   (`ntrees`), the default `select_by` scores by CV instead (with a message); an explicit `select_by = "oob"` warns.
#' - A deterministic forest has no out-of-bag rows, and a forest restored by [rpf_unmarshal()] from a blob made by an
#'   earlier version has no in-bag counts: with the default `select_by`, they are scored by CV instead.
#' - Factor predictors are coded by the outcome of all training rows: their levels are ordered by the mean outcome
#'   (regression), by the share of the second class (binary classification) or by the first principal component of
#'   the class proportions (three or more classes; Coppersmith, Hong and Hosking 1999). So a row's own outcome enters
#'   the coding its OOB families split on: with factors the OOB error is somewhat optimistic.
#'
#' **Cross-validation** (`select_by = "cv"`). For each fold, the forest is refitted with the parameters of `object` on
#' the other folds (after `set.seed(seed)`, the same for every fold), its components, shares and refit are computed
#' on those rows, and every candidate is scored on the held-out fold. The CV RMSE is `sqrt(sum(SSE_k) / n)`, summing
#' the squared errors over the K folds.
#'
#' The chosen pair is then applied to `object` itself (not refitted): its components and shares on `data`, and for
#' `"R"` the least squares on `data`. `data` and `y` must therefore be the rows `object` was trained on, in the order
#' of the fit (checked).
#'
#' **Purification.** `object` is purified in place at all orders, as [purify()] does, if it is not purified yet. This
#' does not change its predictions. A forest purified with a `maxp_interaction` below its order is purified again at
#' all orders, with a message: the pruning needs every component. The out-of-bag choice purifies `object` (again) in
#' mode 2, which gives the same components (a mode-1 purification's up to rounding).
#'
#' **Random numbers.** The folds are drawn after `set.seed(seed)`, and the CV's fold fits after it too. Afterwards the
#' random number stream is restored to its state on entry (after drawing `seed` if it is `NULL`), so the call does not
#' reset the caller's stream.
#'
#' **Cost.** The out-of-bag choice purifies `object` once and evaluates each family's components at its OOB rows; it
#' needs no fit. The CV fits and purifies the forest K more times. Both evaluate the components on the training rows:
#' only those the purified forest stores (at most `sum(choose(p, seq_len(max_interaction)))`, usually far fewer: the
#' subsets of the trees' predictor sets), in one pass over the rows. The memory is about `n` doubles per stored
#' component, a few times that for the out-of-bag choice. `nthreads` speeds up the fits and the purification.
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
#' @param select `[character(1): "min"]`: How the pair is chosen: `"min"` (the lowest RMSE) or `"1se"` (the largest
#'   `tau` within one standard error of it).
#' @param select_by `[character(1): "oob"]`: How the candidates are scored: `"oob"` (out of bag, from one
#'   purification of `object` and no refits) or `"cv"` (K-fold CV with refits of the forest); see Details.
#' @param folds `[integer(1): 5]`: The number of folds, used if `fold_id` is `NULL`: of the CV, or out of bag of R's
#'   refit and of the SE.
#' @param fold_id `[vector | NULL: NULL]`: The fold of every training row: a vector of length `nrow(data)` with at
#'   least two distinct values. If `NULL`, the rows are assigned to `folds` folds at random (after `set.seed(seed)`).
#' @param seed `[integer(1) | NULL: NULL]`: The seed set before the random fold assignment and, for the CV, before
#'   every fold's refit. If `NULL`, one is drawn from the current random number stream and stored in the result, so
#'   that the choice can be reproduced. The caller's random number state is restored afterwards (see Details).
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads for purifying the forests (`object`, and the CV's
#'   fold fits) and for the fold fits. `NULL` uses the `nthreads` `object` was fitted with. The result does not depend
#'   on it.
#' @param verbose `[logical(1): FALSE]`: Report the progress of the choice.
#'
#' @return An object of class `rpf_pruned`, a list with elements:
#' - `fit`: `object` (purified).
#' - `rule`, `tau`: the chosen rule and threshold.
#' - `intercept`: `alpha` (the forest's intercept `f_0` for `"D"`).
#' - `beta`: named vector of the coefficients of the kept components (all 1 for `"D"`).
#' - `sets`: the predictor sets of the kept components, as column indices.
#' - `components`: a `data.frame` of every nonzero component of `object` on the training rows: `component` (its
#'   name, as in [predict_components()]), `order`, `share`, `kept`, `beta` (`NA` if dropped).
#' - `oob`: `NULL` unless the pair was chosen out of bag; otherwise a `data.frame` with one row per candidate: `rule`,
#'   `tau`, `oob_rmse` (corrected; used for the choice), `se_mse` (from its per-fold MSEs), `oob_rmse_raw` and
#'   `se_mse_raw` (uncorrected), `n_kept`, `min` and `one_se` (the two choices). Its attribute `n_rows` is the number
#'   of rows with an OOB family.
#' - `cv`: `NULL` unless the pair was chosen by CV; otherwise a `data.frame` with one row per candidate: `rule`, `tau`,
#'   `cv_rmse`, `se_mse` (SE of the CV MSE), `n_kept` (mean over the folds), `min` and `one_se`.
#' - `select`, `select_by`, `fold_id`, `seed`: as used (`select_by`, `fold_id` and `seed` are `NULL` if no choice was
#'   made).
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
#' # Rule and threshold chosen out of bag (no refits)
#' pruned <- prune_components(fit, x, y, seed = 1)
#' pruned
#' pruned$oob
#'
#' # ... or by 5-fold CV
#' pruned_cv <- prune_components(fit, x, y, select_by = "cv", seed = 1)
#' pruned_cv$cv
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
  select_by = c("oob", "cv"),
  folds = 5L,
  fold_id = NULL,
  seed = NULL,
  nthreads = NULL,
  verbose = FALSE
) {
  rlang::check_dots_empty()
  default_select_by <- missing(select_by)
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
  select_by <- rlang::arg_match(select_by)
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  nthreads <- as.integer(nthreads %||% object$params$nthreads)
  checkmate::assert_flag(verbose)

  tr <- prune_training_data(object, data, y)
  n <- length(tr$y)
  # the shares, R's coefficients and the out-of-bag rows are those of the training rows
  prune_check_training_rows(object, tr)
  cand <- expand.grid(tau = sort(tau), rule = rule, stringsAsFactors = FALSE)[, c("rule", "tau")]
  # how the candidates are scored; NULL for a given rule and tau
  select_by <- if (nrow(cand) > 1) prune_select_by(object, select_by, default_select_by)
  # the out-of-bag choice purifies the forest itself, on the way to its out-of-bag components
  prune_check_purified(object, nthreads, purify = !identical(select_by, "oob"))

  cv <- NULL
  oob <- NULL
  dec <- NULL
  if (!is.null(select_by)) {
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
    if (select_by == "cv") {
      cv <- prune_cv(object, tr, cand, fold_id, seed, nthreads, verbose)
      table <- cv
    } else {
      sel <- prune_oob(object, tr, cand, fold_id, nthreads, verbose)
      oob <- sel$table
      dec <- sel$dec
      table <- oob
    }
    choice <- table[table[[c(min = "min", "1se" = "one_se")[[select]]]], ]
    rule_sel <- choice$rule
    tau_sel <- choice$tau
  } else {
    rule_sel <- cand$rule
    tau_sel <- cand$tau
    fold_id <- NULL
    seed <- NULL
  }

  # The chosen rule and tau on the forest itself, on all training rows (the out-of-bag choice has decomposed it)
  dec <- dec %||% prune_decompose(object, tr$predictors, rep(TRUE, n), nthreads = nthreads)
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
      oob = oob,
      select = select,
      select_by = select_by,
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
# the fold fits would not): it is purified again at all orders. predict() does not change. With purify = FALSE (the
# out-of-bag choice, whose purify_oob() purifies at all orders anyway) only the message.
prune_check_purified <- function(object, nthreads, purify = TRUE) {
  capped <- is_purified(object) && object$fit$is_purify_capped()
  if (capped) {
    cli::cli_inform(c(
      "i" = "{.arg object} was purified with a {.arg maxp_interaction} below its order; purifying it at all orders.",
      " " = "{.fn prune_components} prunes the whole forest. The forest's predictions do not change."
    ))
  }
  if (purify && (capped || !is_purified(object))) {
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

# The same least squares, faster for wide X: a pivoted QR (LAPACK) where X has full numerical rank (|R_kk| above
# 1e-10 |R_11|), there equal to the minimum-norm solution up to rounding; prune_lstsq() otherwise
prune_lstsq_qr <- function(X, y) {
  q <- qr(X, LAPACK = TRUE)
  d <- abs(diag(q$qr))
  if (length(d) < ncol(X) || !all(is.finite(d)) || min(d) <= 1e-10 * max(d)) {
    return(prune_lstsq(X, y))
  }
  as.numeric(qr.coef(q, y))
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
  # the factors that keep their level order in the forest keep it in the fold fits too
  pre <- preprocess_predictors_fit(processed, keep_level_order = object$params$keep_level_order %||% FALSE)
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

# Below this many tree families the default scores by CV instead of out of bag (T072's simulation: from 10 families
# on the out-of-bag choice is as good as the CV's or better)
prune_oob_min_ntrees <- 10L

# How the candidates are scored. An explicit select_by is used as given (out of bag with a warning below
# prune_oob_min_ntrees families). The default is out of bag, but CV for a deterministic forest (no out-of-bag rows),
# a forest without in-bag counts (restored from a blob of an earlier version) and below prune_oob_min_ntrees families.
prune_select_by <- function(object, select_by, default) {
  if (select_by == "cv") {
    return("cv")
  }
  deterministic <- isTRUE(object$params$deterministic)
  if (default && deterministic) {
    return("cv")
  }
  instead <- "the candidates are scored by CV ({.code select_by = \"cv\"}) instead of out of bag."
  if (default && inherits(try(object$fit$inbag(), silent = TRUE), "try-error")) {
    cli::cli_inform(paste("{.arg object} has no in-bag counts (restored from a blob of an earlier version):", instead))
    return("cv")
  }
  ntrees <- object$params$ntrees
  min_ntrees <- prune_oob_min_ntrees
  if (!deterministic && ntrees < min_ntrees) {
    few <- "{.arg object} has {ntrees} tree famil{?y/ies}, fewer than {min_ntrees}:"
    if (default) {
      cli::cli_inform(paste(few, instead))
      return("cv")
    }
    cli::cli_warn(c(
      paste(few, "few training rows have more than one out-of-bag family, and many have none."),
      "i" = "Consider {.code select_by = \"cv\"}."
    ))
  }
  "oob"
}

# Out-of-bag (OOB) choice among the candidates (rule, tau), from one purification of `object` and no refits.
# `object` is purified (mode 2, all orders, as purify() does) and, on the way, every family's components are evaluated
# at the training rows outside its bootstrap sample (purify_oob()). Returns prune_oob_score()'s table and the forest's
# decomposition on its training rows.
prune_oob <- function(object, tr, cand, fold_id, nthreads, verbose, call = rlang::caller_env()) {
  t0 <- proc.time()[["elapsed"]]
  # `tr` are the training rows (checked by prune_components())
  X <- preprocess_predictors_predict(object, tr$predictors)
  o <- object$fit$purify_oob(0L, nthreads)
  dec <- prune_decompose(object, tr$predictors, rep(TRUE, length(tr$y)), nthreads = nthreads)
  out <- prune_oob_score(
    o,
    dec,
    colnames(X),
    tr$y,
    cand,
    fold_id,
    object$params$ntrees,
    deterministic = isTRUE(object$params$deterministic),
    call = call
  )
  if (verbose) {
    n_use <- out$info$n_use
    n_all <- length(tr$y)
    per_row <- sprintf("%.1f", out$info$mean_m)
    n_comp <- length(dec$share)
    secs <- sprintf("%.1f", proc.time()[["elapsed"]] - t0)
    cli::cli_inform(paste(
      "Out of bag: {n_use} of {n_all} rows with out-of-bag families ({per_row} per row),",
      "{n_comp} nonzero component{?s}, {secs} s."
    ))
  }
  out$info <- NULL
  out
}

# The OOB table of every candidate from the OOB sums `o` (as purify_oob() returns them) and the forest's decomposition
# on its training rows `dec` (prune_decompose()). feature_names: the predictors' names; y_all: the training outcome;
# n_trees: the number of families. The kept components of each tau come from the forest's own shares (those of the
# final model). D is scored with the OOB intercept and components. R's (alpha, beta) are fitted as the final model's,
# by least squares on the forest's own components, but on the rows outside each fold (prune_lstsq_qr()), and scored
# with the OOB components of the fold's rows. Rows without an OOB family are left out. The MSE is corrected for the
# OOB means' larger across-family variance (see below). Returns the table, `dec` and `info` (rows used, mean number of
# OOB families per row).
prune_oob_score <- function(
  o,
  dec,
  feature_names,
  y_all,
  cand,
  fold_id,
  n_trees,
  deterministic = FALSE,
  call = rlang::caller_env()
) {
  m0 <- o$count[, 1]
  m1 <- o$count[, 2]
  use <- m0 + m1 > 0
  n_use <- sum(use)
  if (n_use < 2) {
    why <- if (deterministic) " (it is deterministic)" else ""
    cli::cli_abort(
      c(
        "{.code select_by = \"oob\"} needs training rows outside the bootstrap samples; this forest has {n_use}{why}.",
        "i" = "Use {.code select_by = \"cv\"}."
      ),
      call = call
    )
  }
  m0 <- m0[use]
  m1 <- m1[use]
  m <- m0 + m1
  # the OOB sums of the forest's nonzero components (columns of dec$own), per half of the families
  comp <- names(dec$share)
  idx <- match(comp, prune_set_names(o$sets, feature_names))
  half <- function(s) {
    out <- matrix(0, n_use, length(comp))
    out[, !is.na(idx)] <- s[use, idx[!is.na(idx)], drop = FALSE]
    out
  }
  S0 <- half(o$sum0)
  S1 <- half(o$sum1)
  i0 <- o$intercept0[use, 1]
  i1 <- o$intercept1[use, 1]
  rm(o)
  F_oob <- (S0 + S1) / m
  I_oob <- (i0 + i1) / m
  # the means over each half (0 where a half has no OOB family)
  F0 <- S0 / pmax(m0, 1)
  F1 <- S1 / pmax(m1, 1)
  I0 <- i0 / pmax(m0, 1)
  I1 <- i1 / pmax(m1, 1)
  rm(S0, S1)
  both <- m0 > 0 & m1 > 0
  y <- y_all[use]
  fid <- fold_id[use]
  fids <- sort(unique(fid))
  if (length(fids) < 2) {
    cli::cli_abort(
      "The training rows with out-of-bag families must lie in at least two folds of {.arg fold_id}.",
      call = call
    )
  }
  n_k <- vapply(fids, function(k) sum(fid == k), numeric(1))

  # The table of all candidates, with R's coefficients fitted on F_beta of the rows outside each fold
  score <- function(F_beta) {
    sse <- matrix(NA_real_, nrow(cand), length(fids))
    sse_corr <- sse
    n_kept <- integer(nrow(cand))
    cache <- list()
    for (i in seq_len(nrow(cand))) {
      ki <- which(dec$share >= cand$tau[i])
      key <- paste(cand$rule[i], paste(ki, collapse = ","))
      if (is.null(cache[[key]])) {
        if (cand$rule[i] == "D") {
          g <- I_oob + rowSums(F_oob[, ki, drop = FALSE])
          g0 <- I0 + rowSums(F0[, ki, drop = FALSE])
          g1 <- I1 + rowSums(F1[, ki, drop = FALSE])
        } else {
          g <- g0 <- g1 <- numeric(n_use)
          for (k in fids) {
            te <- fid == k
            coef <- prune_lstsq_qr(cbind(1, F_beta[!te, ki, drop = FALSE]), y[!te])
            pred <- function(Fm) coef[1] + as.numeric(Fm[te, ki, drop = FALSE] %*% coef[-1])
            g[te] <- pred(F_oob)
            g0[te] <- pred(F0)
            g1[te] <- pred(F1)
          }
        }
        e2 <- (y - g)^2
        # The OOB mean has m (not n_trees) families: its variance across families is v / m instead of the forest's
        # v / n_trees. v is estimated from the two halves, E (g0 - g1)^2 = v (1 / m0 + 1 / m1). A row whose OOB
        # families are all even or all odd gets the mean estimate of the rows with both halves: which families a row
        # is out of bag for does not depend on the row, so this is unbiased for the rows' mean v.
        v <- (g0 - g1)^2 / (1 / pmax(m0, 1) + 1 / pmax(m1, 1))
        v[!both] <- if (any(both)) mean(v[both]) else 0
        e2_corr <- e2 - (1 / m - 1 / n_trees) * v
        cache[[key]] <- list(
          sse = vapply(fids, function(k) sum(e2[fid == k]), numeric(1)),
          sse_corr = vapply(fids, function(k) sum(e2_corr[fid == k]), numeric(1))
        )
      }
      sse[i, ] <- cache[[key]]$sse
      sse_corr[i, ] <- cache[[key]]$sse_corr
      n_kept[i] <- length(ki)
    }
    # the choice uses the corrected MSE; the uncorrected one is reported as *_raw
    se_of <- function(x) apply(sweep(x, 2, n_k, "/"), 1, stats::sd) / sqrt(length(fids))
    oob_rmse <- sqrt(pmax(rowSums(sse_corr), 0) / n_use)
    se_mse <- se_of(sse_corr)
    ch <- prune_choose(oob_rmse, se_mse, cand)
    table <- data.frame(
      rule = cand$rule,
      tau = cand$tau,
      oob_rmse = oob_rmse,
      se_mse = se_mse,
      oob_rmse_raw = sqrt(rowSums(sse) / n_use),
      se_mse_raw = se_of(sse),
      n_kept = n_kept,
      min = ch$min,
      one_se = ch$one_se,
      stringsAsFactors = FALSE
    )
    attr(table, "n_rows") <- n_use
    table
  }

  # R's coefficients are fitted as the final model's are, on the forest's own components
  list(table = score(dec$own[use, , drop = FALSE]), dec = dec, info = list(n_use = n_use, mean_m = mean(m)))
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
  how <- if (!is.null(x$oob)) {
    sprintf(
      "Chosen by out-of-bag error (%s; R refitted in %d folds, seed %d) among %d candidates.",
      if (x$select == "min") "lowest OOB RMSE" else "1-SE rule",
      length(unique(x$fold_id)),
      x$seed,
      nrow(x$oob)
    )
  } else if (is.null(x$cv)) {
    "Given, not chosen."
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
