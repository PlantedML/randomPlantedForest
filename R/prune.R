#' Prune the components of a Random Planted Forest
#'
#' Drops the small components of a regression forest's functional decomposition and, optionally, refits one
#' coefficient per kept component. The rule and the threshold are chosen out of bag by default, from the forest itself
#' without refits, or by K-fold cross-validation on the training rows. By default every kept interaction of three or
#' more predictors is then coarsened so that it can be shown in at most `max_panels` heat maps, and the model is pruned
#' again.
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
#'   intercept `alpha` and the coefficients `beta_S` by least squares of `y` on the kept components (the minimum-norm
#'   solution if the components are collinear). The components themselves are not changed, only their scale. With
#'   `select_by = "oob"` (the default, see below), the least squares use each training row's out-of-bag components
#'   instead of the forest's own, for a chosen and for a given rule and `tau`. The forest's own components have fitted
#'   the training rows' noise, so least squares on them scale the components up (on California housing the three-way
#'   interactions by about 1.3, where new data wants about 0.9). With `select_by = "cv"`, and where the default falls
#'   back to it (see below), the least squares use the forest's own components on the training rows.
#'
#' **The choice.** If more than one pair of `rule` and `tau` is given (the default: both rules and the grid 0, 0.01%,
#' 0.03%, 0.1%, 0.3%, 1%, 3%, 10%), the pair is chosen by its out-of-bag (`select_by = "oob"`, the default) or
#' cross-validated (`"cv"`) RMSE. With `select = "1se"` (the default), the largest `tau` whose MSE is within one
#' standard error of the minimum is chosen (the SE is the standard deviation of the per-fold MSEs of the minimum,
#' divided by `sqrt(K)`), and at that `tau` the rule with the lower RMSE. With `select = "min"`, the pair with the
#' lowest RMSE is chosen (ties, up to a relative `1e-12`: the larger `tau`, then `"D"`). `"1se"` gives sparse models:
#' with `"R"` fitted out of bag the out-of-bag error changes little over a range of `tau`, and `"min"` keeps small
#' components that add nothing on new data.
#'
#' **Out of bag** (`select_by = "oob"`). Every tree family of `object` is grown on a bootstrap sample of the training
#' rows, and the purification works family by family. So `object` is purified once more (as [purify()] does), and on
#' the way every family's components are evaluated at the training rows outside its bootstrap sample. A row's
#' out-of-bag (OOB) component is the mean over its OOB families (about 37% of them), and its OOB intercept likewise.
#' The purification uses the predictors of all training rows, not the outcome, so no outcome information of a row
#' enters its OOB components.
#' - The kept components of every `tau` are those of the final model (the shares of `object` on the training rows).
#' - `"D"` is scored by the OOB intercept plus the kept OOB components.
#' - For `"R"`, `alpha` and `beta_S` are fitted as for the final model, by least squares on the OOB components, but
#'   only on the rows outside each of K folds (`folds`, `fold_id`, `seed`). They are scored with the OOB components of
#'   the fold's rows, so that no row scores coefficients fitted on it. The final model's are fitted on all rows with
#'   at least one OOB family.
#' - An OOB mean averages `m_i` families (about 37% of `B = ntrees`), the final model all `B`, so the OOB prediction
#'   scatters more around its expectation. The OOB MSE is therefore corrected: per row it subtracts
#'   `(1 / m_i - 1 / B) * v_i`, with `v_i` the across-family variance of the candidate's prediction, estimated from
#'   the OOB means over the even and over the odd families (for a row whose OOB families are all even or all odd: the
#'   mean estimate of the other rows). Without the correction, the noise of the OOB components makes the candidates
#'   look worse than they are. The uncorrected values are reported as `oob_rmse_raw`.
#' - If the correction exceeds a candidate's uncorrected MSE, its corrected MSE is reported as 0. That is not a perfect
#'   candidate but an overshooting, unreliable estimate (it happens with very few rows, e.g. `"R"` with about as many
#'   coefficients as rows, or with almost no noise). Such candidates are not chosen, unless all candidates are at 0
#'   (then they tie).
#' - The OOB RMSE is taken over the rows with at least one OOB family, the SE for `"1se"` from the per-fold MSEs.
#' - With few families the OOB means average only one to a few families, and many rows have none (16% at 4
#'   families). In a simulation (p = 6, n = 800, 60 replicates) the corrected OOB error of the unpruned `"D"` was
#'   within 1.5% of the test error at 4, 10, 20 and 50 families, that of `"R"` within 4.4% except once (+11.5% at 10
#'   families with a 3-way interaction), and the OOB choice was as good as the CV's or better from 10 families on
#'   (within 0.4% in test error, or better by 7-10% where the CV's fold fits choose poorly). Below 10 families
#'   (`ntrees`), the default `select_by` scores by CV instead (with a message); an explicit `select_by = "oob"` warns.
#' - A deterministic forest has no out-of-bag rows, and a forest restored by [rpf_unmarshal()] from a blob made by an
#'   earlier version has no in-bag counts (with a message): with the default `select_by`, they are scored by CV
#'   instead.
#' - A given rule `"R"` and `tau` (no choice) is fitted on the OOB components too, with the same fallbacks: below 10
#'   families or without in-bag counts its coefficients are fitted in sample (with a message), and so they are for a
#'   deterministic forest. A given `"D"` needs no OOB components.
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
#' The chosen (or given) pair is then applied to `object` itself (not refitted): its components and shares on `data`,
#' and for `"R"` the least squares on the OOB components of `data`'s rows (out of bag) or on its components on `data`
#' (by CV). `data` and `y` must therefore be the rows `object` was trained on, in the order of the fit (checked).
#'
#' **The display** (`max_panels`, default 10; `Inf` switches it off). A component of two predictors is one heat map.
#' A kept component `f_S` of three or more predictors is shown as heat maps over two of its predictors, one per cell of
#' the others. If its own cut points need more than `max_panels` such heat maps, it is coarsened:
#' - Two predictors of `S` stay at full resolution. Each other predictor `k` is cut into `q_k` intervals of equal
#'   numbers of training rows (ties merged), with `prod(q) <= max_panels`; `q_k = 1` averages `k` out. Only the
#'   maximal `q` vectors (no `q_k` can grow) are considered; they are listed directly, not searched among all
#'   `max_panels^(|S| - 2)` vectors.
#' - Within a cell, `f_S` is averaged over the coarse predictors with the product of their training marginals. The
#'   coarsened component is not purified again: [purify()] identifies the components under the joint distribution of
#'   the training rows, not under the product of the marginals, and the averaging changes the component's means along
#'   its two full-resolution predictors and its mean over the training rows. It is used as it is (no part of it is
#'   moved into the intercept or into components of lower order).
#' - Of all choices of the two full predictors and the `q` vector, the one keeping the largest share of the variance
#'   of `f_S` on the training rows is taken (ties: the first in lexicographic order).
#'
#' Then the same rule prunes again, on the components the first pruning kept, with the coarsened ones in place of the
#' originals. Its shares, which `tau` thresholds, are of the variance of that model (the intercept plus these
#' components) on the training rows, not of the forest's prediction. A pair chosen out of bag is chosen again out of
#' bag, with the same candidates and folds, scored with each family's coarsened components at its out-of-bag rows (R's
#' coefficients out of bag); a given rule and `tau`, or a pair chosen by CV, is applied again as it is (R's
#' coefficients out of bag if the first pruning fitted them so, else in sample; a CV choice is not repeated, as it would
#' need the display in every fold's refit). The prediction then uses the coarsened components, so the
#' pruning judges the model that is shown. In simulations (48 data sets) this second pruning had 0.49% lower test
#' error than pruning and then coarsening, and 0.57% lower than not coarsening; on California housing (50 families,
#' order 3) the default keeps 30 components, shown in 75 heat maps. If no kept component needs coarsening, the model is
#' the first pruning's.
#' The display costs about as much again as the out-of-bag pruning: one restricted purification for the coarsened
#' components' candidates and one for their out-of-bag values; the stored coarsened components take a few MB each on
#' 16,000 rows.
#'
#' **Purification.** `object` is purified in place at all orders, as [purify()] does, if it is not purified yet. This
#' does not change its predictions. A forest purified with a `maxp_interaction` below its order is purified again at
#' all orders, with a message: the pruning needs every component. The out-of-bag path purifies `object` (again) in
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
#' @param max_panels `[numeric(1): 10]`: The display: the most heat maps a kept component of three or more predictors
#'   may need; such a component is coarsened to at most `max_panels` heat maps, and the model is pruned again (see
#'   Details). A whole number of at least 1, or `Inf` for no display (the first pruning's model).
#' @param rule `[character: c("D", "R")]`: The rules to consider: `"D"` drops, `"R"` drops and refits one coefficient
#'   per kept component.
#' @param tau `[numeric: c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1)]`: The thresholds to consider, as shares of
#'   the variance of the prediction (`0.01` = 1%), each at least 0 and below 1.
#' @param select `[character(1): "1se"]`: How the pair is chosen: `"1se"` (the largest `tau` within one standard
#'   error of the lowest RMSE) or `"min"` (the lowest RMSE).
#' @param select_by `[character(1): "oob"]`: How the candidates are scored, and R's coefficients fitted: `"oob"` (out
#'   of bag, from one purification of `object` and no refits) or `"cv"` (K-fold CV with refits of the forest; R's
#'   coefficients in sample); see Details.
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
#' - `beta`: named vector of the coefficients of the kept components (all 1 for `"D"`). For `"R"` with `select_by`
#'   `"oob"`, `intercept` and `beta` are fitted on the out-of-bag components (see Details).
#' - `sets`: the predictor sets of the kept components, as column indices.
#' - `components`: a `data.frame` of every nonzero component of `object` on the training rows (after a second pruning
#'   of the display, only of those the first pruning kept: see below): `component` (its name, as in
#'   [predict_components()]), `order`, `share`, `kept`, `beta` (`NA` if dropped).
#' - `oob`: `NULL` unless the pair was chosen out of bag; otherwise a `data.frame` with one row per candidate: `rule`,
#'   `tau`, `oob_rmse` (corrected; used for the choice; 0 where the correction exceeds the uncorrected MSE, see
#'   Details), `se_mse` (from its per-fold MSEs), `oob_rmse_raw` and `se_mse_raw` (uncorrected), `n_kept`, `min` and
#'   `one_se` (the two choices). Its attribute `n_rows` is the number of rows with an OOB family.
#' - `cv`: `NULL` unless the pair was chosen by CV; otherwise a `data.frame` with one row per candidate: `rule`, `tau`,
#'   `cv_rmse`, `se_mse` (SE of the CV MSE), `n_kept` (mean over the folds), `min` and `one_se`.
#' - `select`, `select_by`, `fold_id`, `seed`: as used (`fold_id` and `seed` are `NULL` if no choice was made). For a
#'   given rule and `tau`, `select_by` is `"oob"` if R's coefficients were fitted out of bag, otherwise `NULL`.
#' - With the display (finite `max_panels`) and a component coarsened, a second pruning runs, and these describe it:
#'   its candidates are the components the first pruning kept, with the coarsened ones in place, and `oob` is its
#'   table. `components` then lists only these components, and their `share` (which the second pruning's `tau`
#'   thresholds) is of the variance of their model (the intercept plus these components) on the training rows, not of
#'   the forest's prediction; every nonzero component of `object` with its share of the forest's prediction is in
#'   `display$first$components`. Without a coarsened component, they are the first pruning's. With the display,
#'   `components` also has `form` (`"full"` for one or two predictors, `"asis"` if its own cut points need at most
#'   `max_panels` heat maps, `"coarsened"`; `NA` for a component of three or more predictors the first pruning dropped,
#'   so only without a second pruning) and `heat_maps`. The element `display` is a list:
#'   `max_panels`; `first`, the first pruning's model (without `fit`); `coarse`, per coarsened component (named as the
#'   components): `S` (all the component's predictors), `C` (the coarse predictors, `|S| - 2` of `S`: the panel
#'   predictors), `q` (the numbers of intervals aimed at; merged ties can leave fewer, `lengths(cuts) + 1`), `cuts`
#'   (per coarse predictor), `cells` (the number of heat maps, `prod(lengths(cuts) + 1)`), `var_kept` (the share of
#'   the component's variance kept), `forest_panels`, `kept` (in the final model), `candidates` (every candidate's
#'   coarse set, `q`, cells and variance kept), `stored` (what prediction needs; no training data) and `seconds` (the
#'   time to choose and build its form); `asis`, per component kept as it is (form `"asis"`), its panels with the
#'   same names and meaning: `S` (all the component's predictors), `C` (the panel predictors: the `|S| - 2` of `S`
#'   with the fewest own cut points; of ties, the first in `S`), `cuts` (per predictor of `C`, its own cut points: the
#'   leaf bounds of the trees containing `S` inside the range of the training values), `forest_panels`
#'   (`prod(lengths(cuts) + 1)`) and `kept`; `forms`, per component of three or more predictors the first pruning
#'   kept, its form, heat maps and variance kept; `second` (whether a second pruning ran); `seconds` (the time of each
#'   step); `cache`, an environment holding the coarsened components' search indices, built at their first
#'   [`predict()`][predict.rpf_pruned] (a cache, not part of the model: rebuilt when missing or invalid, as after
#'   saving and loading). In `coarse` and `asis` alike, a panel is a cell of `cuts` (a value `x` of predictor `k` of
#'   `C` lies in cell `findInterval(x, cuts[[k]])`), shown as a heat map over the two predictors `setdiff(S, C)` at
#'   full resolution.
#'
#' Use [`predict()`][predict.rpf_pruned] and [predict_components()] on it, and [rpf_marshal()] to save it.
#'
#' @seealso [predict.rpf_pruned()], [predict_components()], [purify()]; [rpf()] with `prune = TRUE` fits and prunes in one
#'   call.
#' @export
#'
#' @examples
#' # A main effect and a pure interaction
#' set.seed(1)
#' x <- matrix(runif(600 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(600, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, ntrees = 20, splits = 50)
#'
#' # Rule and threshold chosen out of bag (no refits): the largest tau within one SE
#' # of the lowest error
#' pruned <- prune_components(fit, x, y, seed = 1)
#' pruned
#' pruned$oob
#'
#' # ... or the lowest error
#' prune_components(fit, x, y, select = "min", seed = 1)
#'
#' # ... or by 5-fold CV
#' pruned_cv <- prune_components(fit, x, y, select_by = "cv", seed = 1)
#' pruned_cv$cv
#'
#' # A given rule and threshold: no choice; R's coefficients are still fitted out of bag
#' prune_components(fit, x, y, rule = "R", tau = 0.01)
#'
#' predict(pruned, x[1:5, ])
#'
#' # A three-way interaction, shown in at most 4 heat maps
#' y3 <- x[, 1] + 4 * x[, 2] * x[, 3] * x[, 4] + rnorm(600, sd = 0.3)
#' fit3 <- rpf(x, y3, max_interaction = 3, ntrees = 20, splits = 100)
#' shown <- prune_components(fit3, x, y3, max_panels = 4, seed = 1)
#' shown
#' shown$display$forms
prune_components <- function(
  object,
  data,
  y = NULL,
  ...,
  max_panels = 10,
  rule = c("D", "R"),
  tau = c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1),
  select = c("1se", "min"),
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
  checkmate::assert_number(max_panels, lower = 1)
  if (is.finite(max_panels)) {
    checkmate::assert_int(max_panels, lower = 1, upper = .Machine$integer.max)
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

  prune_run(
    object,
    prune_training_data(object, data, y),
    max_panels = max_panels,
    rule = rule,
    tau = tau,
    select = select,
    select_by = select_by,
    default_select_by = default_select_by,
    folds = folds,
    fold_id = fold_id,
    seed = seed,
    nthreads = nthreads,
    verbose = verbose
  )
}

# The pruning of prune_components() on the checked arguments: `tr` holds the training rows (prune_training_data()),
# `default_select_by` whether `select_by` was left at its default. rpf(prune = TRUE) calls it with the rows it was
# fitted on. Errors are reported against `call`.
prune_run <- function(
  object,
  tr,
  max_panels,
  rule,
  tau,
  select,
  select_by,
  default_select_by,
  folds,
  fold_id,
  seed,
  nthreads,
  verbose,
  call = rlang::caller_env()
) {
  n <- length(tr$y)
  # the shares, R's coefficients and the out-of-bag rows are those of the training rows
  prune_check_training_rows(object, tr, call = call)
  cand <- expand.grid(tau = sort(tau), rule = rule, stringsAsFactors = FALSE)[, c("rule", "tau")]
  # a choice among the candidates, or a given rule and tau
  choose <- nrow(cand) > 1
  # How the candidates are scored, "oob" or "cv", and so whether R's coefficients are fitted on the out-of-bag
  # components or in sample. A given R is fitted out of bag too, with the choice's fallbacks; a given D, or a given R
  # fitted in sample, has none (NULL)
  select_by <- if (choose || cand$rule == "R") prune_select_by(object, select_by, default_select_by, given = !choose)
  if (!choose && identical(select_by, "cv")) {
    select_by <- NULL
  }
  # the out-of-bag path purifies the forest itself, on the way to its out-of-bag components
  prune_check_purified(object, nthreads, purify = !identical(select_by, "oob"))

  cv <- NULL
  oob <- NULL
  dec <- NULL
  # the out-of-bag components that R's final coefficients are fitted on, if out of bag
  oob_fit <- NULL
  # the display prunes again out of bag: it needs the forest's OOB sums
  display <- is.finite(max_panels)
  sums <- NULL
  if (choose) {
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
        cli::cli_abort("{.arg fold_id} must have at least two distinct folds.", call = call)
      }
    }
    if (select_by == "cv") {
      cv <- prune_cv(object, tr, cand, fold_id, seed, nthreads, verbose)
      table <- cv
    } else {
      sel <- prune_oob(object, tr, cand, fold_id, nthreads, verbose, keep_sums = display)
      oob <- sel$table
      dec <- sel$dec
      oob_fit <- sel$oob_fit
      sums <- sel$sums
      rm(sel)
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
    if (identical(select_by, "oob")) {
      # a given R: the out-of-bag components for its coefficients, as the choice fits them (no scores, no folds)
      sel <- prune_oob(object, tr, NULL, NULL, nthreads, verbose, keep_sums = display)
      dec <- sel$dec
      oob_fit <- sel$oob_fit
      sums <- sel$sums
      rm(sel)
    }
  }

  # The chosen rule and tau on the forest itself, on all training rows (the out-of-bag path has decomposed it). R's
  # coefficients: on the out-of-bag components with select_by = "oob", as each candidate was scored; else in sample
  dec <- dec %||% prune_decompose(object, tr$predictors, rep(TRUE, n), nthreads = nthreads)
  rf <- prune_fit_rule(rule_sel, tau_sel, dec, tr$y, oob = oob_fit)
  res <- prune_result(object, rule_sel, tau_sel, dec, rf, cv, oob, select, select_by, fold_id, seed)
  if (display) {
    rm(oob_fit)
    res <- prune_display(res, object, tr, dec, sums, cand, as.integer(max_panels), nthreads)
  }
  res
}

# The rpf_pruned model of rule and tau on the decomposition `dec`, with R's (alpha, beta) `rf` (prune_fit_rule())
prune_result <- function(object, rule, tau, dec, rf, cv, oob, select, select_by, fold_id, seed) {
  sets <- dec$sets[dec$nonzero]
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
      rule = rule,
      tau = tau,
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
  intercept <- as.numeric(fit$fit$predict_matrix(X[1, 1, drop = FALSE], -1, 0L))[1]
  prune_dec_matrix(M, sets, intercept, sum(own), colnames(X), !is.null(other))
}

# The decomposition from the components M (rows: the own rows, then the others; columns: `sets`) and the intercept.
# The prediction is the intercept plus all columns of M.
prune_dec_matrix <- function(M, sets, intercept, n_own, feature_names, other = FALSE) {
  M_own <- M[seq_len(n_own), , drop = FALSE]
  # the purified forest's prediction is its intercept plus all its components, which are all in M
  pred_own <- intercept + rowSums(M_own)
  var_pop <- function(v) mean((v - mean(v))^2)
  v_comp <- vapply(seq_len(ncol(M_own)), function(k) var_pop(M_own[, k]), numeric(1))
  nonzero <- v_comp > 0
  share <- v_comp[nonzero] / var_pop(pred_own)
  names(share) <- prune_set_names(sets[nonzero], feature_names)
  M <- M[, nonzero, drop = FALSE]
  colnames(M) <- names(share)
  list(
    sets = sets,
    nonzero = nonzero,
    share = share,
    intercept = intercept,
    own = M[seq_len(n_own), , drop = FALSE],
    other = if (other) M[-seq_len(n_own), , drop = FALSE]
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

# The kept components of rule and tau, and R's (alpha, beta) by least squares on the decomposition's own rows; with
# `oob` (chosen out of bag: prune_oob_score()'s `oob_fit`, list(F = the out-of-bag components of the rows with OOB
# families, columns as names(dec$share), y = their outcome)) on the out-of-bag components instead
prune_fit_rule <- function(rule, tau, dec, y_own, oob = NULL) {
  kept <- names(dec$share)[dec$share >= tau]
  if (rule == "D") {
    return(list(alpha = dec$intercept, beta = stats::setNames(rep(1, length(kept)), kept)))
  }
  coef <- if (is.null(oob)) {
    prune_lstsq(cbind(1, dec$own[, kept, drop = FALSE]), y_own)
  } else {
    prune_lstsq(cbind(1, oob$F[, match(kept, names(dec$share)), drop = FALSE]), oob$y)
  }
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

# The two choices among the candidates from their RMSE and the SE of their MSE (logical vectors). left_out: candidates
# not chosen unless all are left out (prune_oob_score(): those whose corrected MSE was set to 0, an overshooting
# correction); they then tie at 0. This is the rule of oob_lowest(), which the adaptive stop and rpf(tune = ) use
prune_choose <- function(rmse, se_mse, cand, left_out = NULL) {
  if (any(left_out) && !all(left_out)) {
    use <- !left_out
    ch <- prune_choose(rmse[use], se_mse[use], cand[use, , drop = FALSE])
    spread <- function(x) replace(logical(length(rmse)), use, x)
    return(list(min = spread(ch$min), one_se = spread(ch$one_se)))
  }
  rule_rank <- match(cand$rule, c("D", "R"))
  # min: the lowest RMSE by oob_lowest() (ties within relative 1e-12): the larger tau, then D. At a lowest RMSE of 0
  # the candidates at 0 tie (all left out, or a CV RMSE of 0, which is a perfect fit, not an overshooting correction)
  ties <- if (min(rmse) > 0) oob_lowest(rmse) else which(rmse == 0)
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

# How the candidates are scored, and so how R's coefficients are fitted ("oob": on the out-of-bag components; "cv": in
# sample). An explicit select_by is used as given (out of bag with a warning below prune_oob_min_ntrees families). The
# default is out of bag, but CV for a deterministic forest (no out-of-bag rows), a forest without in-bag counts
# (restored from a blob of an earlier version) and below prune_oob_min_ntrees families. given: for a given rule R
# (no candidates to score), whose coefficients are then fitted in sample (the messages say so).
prune_select_by <- function(object, select_by, default, given = FALSE) {
  if (select_by == "cv") {
    return("cv")
  }
  deterministic <- isTRUE(object$params$deterministic)
  if (default && deterministic) {
    return("cv")
  }
  instead <- if (given) {
    "R's coefficients are fitted in sample instead of out of bag."
  } else {
    "the candidates are scored by CV ({.code select_by = \"cv\"}) instead of out of bag."
  }
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
# at the training rows outside its bootstrap sample (purify_oob()). Returns prune_oob_score()'s table, the forest's
# decomposition on its training rows and the out-of-bag components for R's final fit (`oob_fit`). With cand = NULL
# (a given rule R), no table: only `dec` and `oob_fit`. keep_sums: also the OOB sums of purify_oob() (`sums`), which
# the display's second pruning needs.
prune_oob <- function(object, tr, cand, fold_id, nthreads, verbose, keep_sums = FALSE, call = rlang::caller_env()) {
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
  if (keep_sums) {
    out$sums <- o
  }
  out
}

# The OOB table of every candidate from the OOB sums `o` (as purify_oob() returns them) and the forest's decomposition
# on its training rows `dec` (prune_decompose()). feature_names: the predictors' names; y_all: the training outcome;
# n_trees: the number of families. The kept components of each tau come from the forest's own shares (those of the
# final model). D is scored with the OOB intercept and components. R's (alpha, beta) are fitted as the final model's,
# by least squares on the OOB components (not the forest's own, which have fitted the training rows' noise), but on
# the rows outside each fold (prune_lstsq_qr()), and scored with the OOB components of the fold's rows. Rows without
# an OOB family are left out. The MSE is corrected for the OOB means' larger across-family variance (see below).
# Returns the table, `dec`, `oob_fit` (list(F = the OOB components, y = the outcome) of the rows with OOB families,
# for R's final fit by prune_fit_rule()) and `info` (rows used, mean number of OOB families per row). With
# cand = NULL (a given rule R) nothing is scored and the table is NULL.
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
  # a given rule R (cand = NULL, no candidates to score): only the OOB components for its coefficients
  if (is.null(cand)) {
    return(list(
      table = NULL,
      dec = dec,
      oob_fit = list(F = F_oob, y = y_all[use]),
      info = list(n_use = n_use, mean_m = mean(m))
    ))
  }
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
    ch <- prune_choose(oob_rmse, se_mse, cand, left_out = oob_rmse == 0)
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

  # R's coefficients are fitted as the final model's are, on the OOB components
  list(
    table = score(F_oob),
    dec = dec,
    oob_fit = list(F = F_oob, y = y),
    info = list(n_use = n_use, mean_m = mean(m))
  )
}

# The kept components (scaled by beta) on new data, as a matrix with one named column per component
prune_kept_matrix <- function(object, new_data, nthreads = 0L, call = rlang::caller_env()) {
  fit <- object$fit
  check_rpf_alive(fit, call = call)
  processed <- hardhat::forge(new_data, fit$blueprint)
  check_no_missing(processed$predictors, "{.arg new_data}", call = call)
  X <- preprocess_predictors_predict(fit, processed$predictors)
  coarse <- names(object$beta) %in% names(object$display$coarse)
  if (any(coarse)) {
    # the display's coarsened components from their stored form, the others from the forest
    M <- matrix(0, nrow(X), length(object$beta))
    if (!all(coarse)) {
      M[, !coarse] <- prune_component_matrix(fit, X, object$sets[!coarse], nthreads)
    }
    for (k in which(coarse)) {
      M[, k] <- rpf_coarse_eval(prune_display_built(object, names(object$beta)[k]), X)
    }
  } else {
    M <- prune_component_matrix(fit, X, object$sets, nthreads)
  }
  M <- sweep(M, 2, object$beta, "*")
  colnames(M) <- names(object$beta)
  M
}

#' Predictions of a pruned Random Planted Forest
#'
#' The prediction of a pruned forest is `alpha + sum(beta_S * f_S(x))` over the kept components `f_S` of the forest,
#' with the intercept `alpha` and the coefficients `beta_S` (see [prune_components()]).
#'
#' @param object `[rpf_pruned]`: A pruned forest from [prune_components()] or [rpf()].
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
#' @param x `[rpf_pruned]`: A pruned forest from [prune_components()] or [rpf()].
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
  # chosen out of bag, R's coefficients were fitted on the out-of-bag components
  if (x$rule == "R" && identical(x$select_by, "oob")) {
    rule_name <- sub(")$", ", fitted out of bag)", rule_name)
  }
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
  if (!is.null(x$tune)) {
    st <- x$tune$settings
    lines <- c(
      sprintf(
        "Tuned out of bag (rpf(tune = )): %d forests; chosen split_try = %d, t_try = %s, min_leaf_rows = %d, splits = %d.",
        nrow(x$tune$search),
        st$split_try,
        format(st$t_try),
        st$min_leaf_rows,
        st$splits
      ),
      "",
      lines
    )
  }
  if (!is.null(x$display)) {
    lines <- c(lines, "", format_display(x, n_show))
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
    what <- "share on the training rows and coefficient"
    if (!is.null(x$display)) {
      tab$heat_maps <- show$heat_maps
      what <- "share on the training rows, coefficient and heat maps"
    }
    lines <- c(
      lines,
      sprintf("Kept components (largest first), %s:", what),
      utils::capture.output(print(tab, row.names = FALSE, right = FALSE))
    )
    if (nrow(kept) > n_show) {
      lines <- c(lines, sprintf("... and %d more", nrow(kept) - n_show))
    }
  }

  # format into lines so print() writes to stdout, not cli's message stream
  c(cli::cli_format_method(cli::cli_rule(left = "Pruned Random Planted Forest")), lines)
}

# The display's lines of format.rpf_pruned(): the first pruning, the coarsened components, the heat maps
format_display <- function(x, n_show) {
  d <- x$display
  f <- d$first
  how <- if (d$second) "coarsened, then pruned again" else "nothing to coarsen, so not pruned again"
  lines <- c(
    sprintf("Display: every kept component in at most %d heat maps (%s).", d$max_panels, how),
    sprintf("First pruning: rule %s, tau = %s, %d components kept.", f$rule, format_share(f$tau), length(f$beta))
  )
  co <- d$coarse
  if (length(co) > 0) {
    kept_co <- vapply(co, function(z) isTRUE(z$kept), logical(1))
    vk <- round(range(vapply(co, function(z) z$var_kept, numeric(1))), 2)
    lines <- c(
      lines,
      sprintf(
        "Coarsened: %d component%s (%d kept by the second pruning); variance kept %s.",
        length(co),
        if (length(co) == 1) "" else "s",
        sum(kept_co),
        if (vk[1] == vk[2]) format(vk[1]) else paste(vk, collapse = "-")
      )
    )
    if (any(kept_co) && n_show > 0) {
      z <- utils::head(co[kept_co], n_show)
      tab <- data.frame(
        component = names(z),
        coarse = vapply(z, function(v) paste0(v$C, "(", lengths(v$cuts) + 1, ")", collapse = " x "), ""),
        heat_maps = vapply(z, function(v) v$cells, numeric(1)),
        var_kept = vapply(z, function(v) format(v$var_kept, digits = 3), "")
      )
      lines <- c(
        lines,
        "Kept coarsened components: coarse features (intervals), heat maps, variance kept:",
        utils::capture.output(print(tab, row.names = FALSE, right = FALSE))
      )
      if (sum(kept_co) > n_show) {
        lines <- c(lines, sprintf("... and %d more", sum(kept_co) - n_show))
      }
    }
  }
  c(lines, sprintf("Heat maps of the model: %g", sum(x$components$heat_maps[x$components$kept])))
}

format_share <- function(s) {
  if (s == 0) {
    return("0")
  }
  paste0(trimws(formatC(signif(100 * s, 2), digits = 2, format = "fg")), "%")
}
