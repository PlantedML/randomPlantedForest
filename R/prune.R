#' Prune the components of a Random Planted Forest
#'
#' Drops the small components of a forest's functional decomposition and, optionally, refits one coefficient per
#' kept component. The rule and the threshold are chosen out of bag by default, from the forest itself
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
#' rows, and the purification works family by family. So every family's components (see **Purification** below) are
#' also evaluated at the training rows outside its bootstrap sample. A row's out-of-bag (OOB) component is the mean
#' over its OOB families (about 37% of them), and its OOB intercept likewise.
#' The purification uses the predictors of all training rows, not the outcome, so no outcome information of a row
#' enters its OOB components.
#' - The kept components of every `tau` are those of the final model (the shares of `object` on the training rows).
#' - `"D"` is scored by the OOB intercept plus the kept OOB components.
#' - For `"R"`, `alpha` and `beta_S` are fitted as for the final model, by least squares on the OOB components, but
#'   only on the rows outside each of K folds (`folds`, `fold_id`, `seed`). They are scored with the OOB components of
#'   the fold's rows, so that no row scores coefficients fitted on it. The final model's are fitted on all rows with
#'   at least one OOB family.
#' - An OOB mean averages `m_i` families (about 37% of `B = n_families`), the final model all `B`, so the OOB prediction
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
#'   (`n_families`), the default `select_by` scores by CV instead (with a message); an explicit `select_by = "oob"` warns.
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
#' The display costs about as much again as the out-of-bag pruning: a purification restricted to the components being
#' coarsened (two at a time) for their candidates, and one for their out-of-bag values; the stored coarsened components
#' take a few MB each on 16,000 rows at order 3, more at higher orders.
#'
#' **Purification.** The pruning needs every component's values at the training rows (and out of bag, its OOB
#' means), not the purified forest: [purify()] keeps every component of every tree family as pieces, which at high
#' orders takes far more memory than the forest (see **Cost**). So each family is
#' purified in turn, as [purify()] does (mode 2, all orders), its components are evaluated at the training rows and
#' added up, and the family's pieces are freed. The pruned model then stores the pieces of its kept components only,
#' merged over the families exactly as [purify()] stores them, and [`predict()`][predict.rpf_pruned] and
#' [predict_components()] evaluate them from there; a coarsened component (the display) is stored in its coarsened
#' form. The components equal those of the purified forest up to rounding (about `1e-15` relative), and the stored
#' ones are the purified forest's bit for bit. `object` itself is not purified: an earlier [purify()], at any
#' `maxp_interaction` or `mode`, is left as it is and does not change the result.
#'
#' **Random numbers.** The folds are drawn after `set.seed(seed)`, and the CV's fold fits after it too. Afterwards the
#' random number stream is restored to its state on entry (after drawing `seed` if it is `NULL`), so the call does not
#' reset the caller's stream.
#'
#' **Cost.** The out-of-bag choice purifies the families of `object` once (see **Purification**); it needs no fit. The
#' CV fits and purifies the forest K more times (each fold's forest as [purify()] does). Only the components the
#' forest has are evaluated (at most `sum(choose(p, seq_len(max_interaction)))`, usually far fewer: the subsets of the
#' trees' predictor sets). The memory is the forest's plus about `n` doubles per component, three times that for the
#' out-of-bag choice, plus one family's components at a time per thread; storing the kept components passes over the
#' families once more. On California housing (16,512 rows, 200 families, `n_try_split_points = 3`, `share_candidates = 0.1`, the fast
#' draw, 8 threads; measured when this column pass was added) the whole call of `rpf(n_splits = NULL, prune = TRUE)`
#' peaked at 2.3 GB of physical footprint (4.1 GB of resident memory, RSS) at order 3, and at 6.4 GB (9.9 GB RSS) at
#' order 8 (all predictors), most of it the 31 coarsened components the display keeps; there [purify()] alone needs
#' more than 12 GB of RSS. `nthreads` speeds up the fits and the purification.
#'
#' **Classification.** A classification forest is pruned the same way, with its score columns: one for a binary
#' outcome, one per class otherwise (`K - 1` for `loss = "logit"`, whose first level is the reference). Every component
#' then has a column per score column, and:
#' - A component's share sums its variance over the score columns, centred across the classes (`"logit"`: with the
#'   reference level's score 0 joined; `"exponential"`: the part all classes share does not count), over the same sum
#'   for the prediction. One kept list serves all classes.
#' - `"R"` fits one coefficient per kept component, shared by the classes, and one intercept per score column
#'   (`"exponential"` with more than two classes: the first is 0), out of bag as for regression.
#' - The criterion is the forest's loss's own: for `loss = "L1"` or `"L2"` the squared error of the raw scores summed
#'   over the columns (the Brier score; least squares for `"R"`), for `"logit"` and `"exponential"` the log-loss of the
#'   probabilities that `predict(type = "prob")` makes (maximum likelihood for `"R"`, by Newton's method from rule
#'   `"D"`'s model). The OOB correction is its second-order form (see [rpf()], "Classification"). The tables then have
#'   `oob_loss`, `se_loss`, `oob_loss_raw` and `se_loss_raw` (the mean criterion) instead of the RMSE columns, and
#'   `cv_loss` and `se_loss` for CV; `select = "1se"` uses the SE of the mean criterion. A binary forest with
#'   `"L1"` or `"L2"` has one 0/1 column and is pruned as a regression forest of that response, with the RMSE columns.
#' - The display coarsens a component with one form (`C`, `q`, cuts) for all score columns, chosen by the variance
#'   kept of the centred columns; its stored form is one per score column.
#' - [`predict()`][predict.rpf_pruned] takes `type` as [predict.rpf()] does, and [predict_components()] has one column
#'   per component and score column, named as for the forest.
#'
#' @param object `[rpf]`: A fitted [`rpf`] model (regression or classification).
#' @param data `[data.frame | matrix]`: The training data of `object`, in the order of the fit: the predictors as
#'   passed to [rpf()] as `x`, or the `data` of the formula or recipe interface. For a formula or recipe fit, `data`
#'   may contain the outcome.
#' @param y `[numeric | factor | NULL: NULL]`: The outcome of the training rows (the factor for classification).
#'   Required unless `data` contains the outcome of a formula or recipe fit.
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
#' - `fit`: `object`, as it was (not purified by the pruning).
#' - `rule`, `tau`: the chosen rule and threshold.
#' - `intercept`: `alpha` (the forest's intercept `f_0` for `"D"`); for classification one per score column.
#' - `beta`: named vector of the coefficients of the kept components (all 1 for `"D"`). For `"R"` with `select_by`
#'   `"oob"`, `intercept` and `beta` are fitted on the out-of-bag components (see Details).
#' - `sets`: the predictor sets of the kept components, as column indices.
#' - `components`: a `data.frame` of every nonzero component of `object` on the training rows (after a second pruning
#'   of the display, only of those the first pruning kept: see below): `component` (its name, as in
#'   [predict_components()]), `order`, `share`, `kept`, `beta` (`NA` if dropped), and `model_share`:
#'   the share of the final model's component, `Var(beta_S f_S) / Var(prediction of the model)` on the training rows
#'   (with the display, of `f_S` in its heat-map form; for classification, as `share`, over the centred score columns;
#'   `NA` if dropped). `share` is what `tau` thresholds; `model_share` is the component's importance in the model,
#'   which [print()][print.rpf_pruned()] shows. As the components need not be uncorrelated, the shares need not sum to
#'   1.
#' - `oob`: `NULL` unless the pair was chosen out of bag; otherwise a `data.frame` with one row per candidate: `rule`,
#'   `tau`, `oob_rmse` (corrected; used for the choice; 0 where the correction exceeds the uncorrected MSE, see
#'   Details), `se_mse` (from its per-fold MSEs), `oob_rmse_raw` and `se_mse_raw` (uncorrected), `n_kept`, `min` and
#'   `one_se` (the two choices). Its attribute `n_rows` is the number of rows with an OOB family.
#' - `cv`: `NULL` unless the pair was chosen by CV; otherwise a `data.frame` with one row per candidate: `rule`, `tau`,
#'   `cv_rmse`, `se_mse` (SE of the CV MSE), `n_kept` (mean over the folds), `min` and `one_se`.
#' - `select`, `select_by`, `fold_id`, `seed`: as used (`fold_id` and `seed` are `NULL` if no choice was made). For a
#'   given rule and `tau`, `select_by` is `"oob"` if R's coefficients were fitted out of bag, otherwise `NULL`.
#' - `criterion` (classification only): `"brier"` or `"logloss"`, what the candidates were scored by.
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
#' - `store`: the pieces of the kept components that are not coarsened, from which the model predicts (see
#'   **Purification**), and `cache`, an environment holding them built for prediction.
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
#' fit <- rpf(x, y, max_interaction = 2, n_families = 20, n_splits = 50, prune = FALSE)
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
#' fit3 <- rpf(x, y3, max_interaction = 3, n_families = 20, n_splits = 100, prune = FALSE)
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
# fitted on. store = FALSE leaves out the stored pieces of the kept components (rpf(tune = )'s search scores the
# pruned models but does not predict with them). Errors are reported against `call`.
#
# The forest's components come from one pass over its families (prune_columns()), not from a purified store, and the
# result stores the pieces of its kept components only (prune_store()): the forest is left as it is.
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
  store = TRUE,
  call = rlang::caller_env()
) {
  n <- length(tr$y)
  # the leaves' range counts that the pass over the families keeps for the later passes are dropped at the end
  on.exit(object$fit$drop_purify_counts(), add = TRUE)
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
  # the predictor names as the components are named
  fnames <- colnames(preprocess_predictors_predict(object, tr$predictors[1, , drop = FALSE]))

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
      sel <- prune_oob(object, tr, fnames, cand, fold_id, nthreads, verbose, keep_sums = display, call = call)
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
      sel <- prune_oob(object, tr, fnames, NULL, NULL, nthreads, verbose, keep_sums = display, call = call)
      dec <- sel$dec
      oob_fit <- sel$oob_fit
      sums <- sel$sums
      rm(sel)
    }
  }

  # The chosen rule and tau on the forest itself, on all training rows (the out-of-bag path has its components). R's
  # coefficients: on the out-of-bag components with select_by = "oob", as each candidate was scored; else in sample
  dec <- dec %||% prune_columns(object, fnames, nthreads, oob = FALSE, crit = tr$crit)$dec
  rf <- prune_fit_rule(rule_sel, tau_sel, dec, prune_score_y(tr), oob = oob_fit)
  res <- prune_result(object, rule_sel, tau_sel, dec, rf, cv, oob, select, select_by, fold_id, seed)
  if (display) {
    rm(oob_fit)
    res <- prune_display(res, object, tr, dec, sums, cand, as.integer(max_panels), nthreads)
  }
  # classification: what the candidates were scored by
  res$criterion <- tr$crit$type %||% tr$score_type
  if (store) {
    # the pieces of the kept components that are not coarsened (the display stores those in its own form)
    coarse <- names(res$beta) %in% names(res$display$coarse)
    res$store <- prune_store(object, res$sets[!coarse], nthreads)
    res$cache <- new.env(parent = emptyenv())
    # each component's share of the final model, on the training rows
    res$components$model_share <- rep(NA_real_, nrow(res$components))
    X <- preprocess_predictors_predict(object, tr$predictors)
    res$components$model_share[match(names(res$beta), res$components$component)] <-
      prune_model_share(res, X, tr$crit, nthreads)
  }
  res
}

# The share of each component of the model `object` (an rpf_pruned with its store) on the rows X (preprocessed):
# Var(beta_S f_S) / Var(prediction), f_S in its final form (stored or heat-map form). With V > 1 score columns the
# variances are those of the score columns centred as `crit` centres them (crit_centred_var(), as the shares of the
# decomposition). In the order of object$beta.
prune_model_share <- function(object, X, crit, nthreads) {
  V <- pruned_V(object)
  if (length(object$beta) == 0) {
    return(numeric(0))
  }
  M <- prune_kept_matrix_x(object, X, nthreads)
  size <- if (V == 1L) {
    function(A) mean((A - mean(A))^2)
  } else {
    crit <- crit %||% list(centre = "none")
    function(A) crit_centred_var(crit, A)
  }
  total <- size(pruned_scores(object, M))
  vapply(seq_along(object$beta), function(k) size(M[, crit_blk(k, V), drop = FALSE]), numeric(1)) / total
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
# (prune_rows())
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
  check_no_missing(processed$predictors, "{.arg data}", call = call)
  prune_rows(object, processed$predictors, y, call = call)
}

# The training rows `tr` of the pruning from the forged predictors and the outcome y: list(predictors, y), with y the
# numeric outcome. For classification y is the factor outcome (with the forest's levels), and `tr` also holds
# `outcome` (the factor, for the CV's refits), `Y` (the outcome as the C++ forest codes it, for the check that these are
# the training rows) and the criterion (oob_score_y()): a binary L1 / L2 forest's 0/1 column is scored as a regression
# response (y the 0/1 column, score_type "brier"), otherwise `crit` is the classification criterion.
prune_rows <- function(object, predictors, y, call = rlang::caller_env()) {
  n <- nrow(predictors)
  if (identical(object$mode, "classification")) {
    lev <- levels(object$blueprint$ptypes$outcomes[[1]])
    if (!is.factor(y)) {
      y <- factor(y, levels = lev)
    }
    if (!identical(levels(y), lev) || anyNA(y)) {
      cli::cli_abort(
        "{.arg y} must be the factor outcome of {.arg object} (levels {.val {lev}}) without missing values.",
        call = call
      )
    }
  } else if (!is.numeric(y) || anyNA(y)) {
    cli::cli_abort("{.arg y} must be a numeric outcome without missing values.", call = call)
  }
  if (length(y) != n) {
    cli::cli_abort("{.arg y} has {length(y)} value{?s} but {.arg data} has {n} row{?s}.", call = call)
  }
  if (n < 2) {
    cli::cli_abort("{.arg data} must have at least two rows.", call = call)
  }
  if (!identical(object$mode, "classification")) {
    return(list(predictors = predictors, y = as.numeric(y)))
  }
  loss <- object$params$loss
  Y <- preprocess_outcome(list(outcomes = list(y)), loss, call = call)$outcomes
  crit <- oob_criterion("classification", loss, y)
  if (crit$V != ncol(Y)) {
    cli::cli_abort("{.arg y} must have every level of the outcome of {.arg object}.", call = call)
  }
  if (crit_is_scalar_mse(crit)) {
    return(list(predictors = predictors, y = crit$y, outcome = y, Y = Y, score_type = crit$type))
  }
  list(predictors = predictors, y = y, outcome = y, Y = Y, crit = crit)
}

# What the candidates are scored by: the numeric outcome, or the classification criterion of the training rows
prune_score_y <- function(tr) {
  tr$crit %||% tr$y
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
  Y <- tr$Y %||% tr$y
  same <- identical(dim(X), dim(train$X)) && all(X == train$X) && length(Y) == length(train$Y) && all(Y == train$Y)
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

# The forest's components at its training rows from one pass over its families (purify_columns(), C++): each family
# is purified as purify() purifies it (mode 2, all orders), its components are evaluated at all training rows, and
# then it is freed, so the purified store of all components is never built and the forest is left as it is (an
# earlier purification, at any order or mode, stays as it is and is not read). Returns `dec`, the decomposition as
# prune_decompose() gives it (the components in the order of prune_order_sets(), equal to the purified store's up to
# rounding, its intercept bitwise), and with `oob` also `o`, the out-of-bag sums as purify_oob() gives them
# (prune_oob_score()). feature_names: the predictors' names. crit: the classification criterion (prune_rows()), with
# V score columns per component (component-major, as the C++ columns), or NULL.
prune_columns <- function(object, feature_names, nthreads, oob = TRUE, crit = NULL) {
  cc <- object$fit$purify_columns(0L, nthreads, oob)
  V <- length(cc$intercept)
  sets <- lapply(cc$sets, as.integer)
  ord <- prune_order_sets(sets)
  # a forest of order k has no component above k; as prune_stored_sets()
  ord <- ord[lengths(sets[ord]) <= object$params$max_interaction]
  M <- cc$own[, crit_blk(ord, V), drop = FALSE]
  cc$own <- NULL
  intercept <- if (V == 1L) cc$intercept[1] else cc$intercept
  dec <- prune_dec_matrix(M, sets[ord], intercept, nrow(M), feature_names, crit = crit)
  rm(M)
  list(dec = dec, o = if (oob) cc[c("sets", "sum0", "sum1", "intercept0", "intercept1", "count")])
}

# The pieces of the components `sets` (column indices) of `object`, merged exactly as purify() stores them (C++
# store_components(): one pass over the families restricted to these components), as plain R data: what predict()
# needs for them, without the forest's other components
prune_store <- function(object, sets, nthreads) {
  object$fit$store_components(lapply(sets, as.integer), nthreads)
}

# The built store of a pruned model (an external pointer: the store's tables and indices, as the purified forest
# builds them), kept in the model's cache environment; rebuilt if missing or gone (a deserialized model)
prune_store_built <- function(object) {
  cache <- object$cache
  xp <- if (is.environment(cache)) cache[[".store"]] else NULL
  if (is.null(xp) || !rpf_store_valid(xp)) {
    xp <- rpf_store_build(object$store, as.integer(object$fit$params$nthreads))
    if (is.environment(cache)) {
      assign(".store", xp, envir = cache)
    }
  }
  xp
}

# The order of predictor sets (column indices): by order, then lexicographically (the order of utils::combn() over
# all sets)
prune_order_sets <- function(sets) {
  if (length(sets) == 0) {
    return(integer(0))
  }
  width <- max(lengths(sets))
  key <- vapply(sets, function(s) c(length(s), s, rep(0L, width - length(s))), integer(width + 1))
  do.call(order, as.data.frame(t(key)))
}

# The predictor sets of the components the purified forest stores (as column indices), up to its order, ordered by
# order and then lexicographically (the order of utils::combn() over all sets)
prune_stored_sets <- function(fit) {
  sets <- lapply(fit$fit$purified_components(), as.integer)
  sets <- sets[lengths(sets) >= 1 & lengths(sets) <= fit$params$max_interaction]
  sets[prune_order_sets(sets)]
}

# Names of predictor sets, as predict_components() names its columns
prune_set_names <- function(sets, feature_names) {
  vapply(sets, function(s) paste0(sort(feature_names[s]), collapse = ":"), character(1))
}

# Decomposition of `fit` on its own training rows (logical `own`) and, optionally, other rows (logical `other`): the
# components, the intercept, the shares of the nonzero components (named) and which sets they are. Only the
# components the purified forest stores are evaluated (the others are zero), not all 2^p - 1 sets. `fit` is purified
# first if it is not: the CV's fold fits, which need their components at the held-out rows and are discarded
# afterwards (the forest being pruned is decomposed by prune_columns() instead).
prune_decompose <- function(fit, predictors, own, other = NULL, nthreads = 0L, crit = NULL) {
  X <- preprocess_predictors_predict(fit, predictors)
  if (!is_purified(fit)) {
    fit$fit$purify_threads(0L, nthreads, 2L)
  }
  sets <- prune_stored_sets(fit)
  rows <- c(which(own), if (!is.null(other)) which(other))
  M <- fit$fit$predict_components_matrix(X[rows, , drop = FALSE], sets, nthreads)
  dimnames(M) <- NULL
  intercept <- as.numeric(fit$fit$predict_matrix(X[1, 1, drop = FALSE], -1, 0L))
  if (length(intercept) == 1L) {
    intercept <- intercept[1]
  }
  prune_dec_matrix(M, sets, intercept, sum(own), colnames(X), !is.null(other), crit = crit)
}

# The decomposition from the components M (rows: the own rows, then the others; columns: `sets`) and the intercept.
# The prediction is the intercept plus all columns of M. With a criterion of V > 1 score columns (crit), M has V
# columns per set: prune_dec_matrix_v().
prune_dec_matrix <- function(M, sets, intercept, n_own, feature_names, other = FALSE, crit = NULL) {
  if (!is.null(crit) && crit$V > 1L) {
    return(prune_dec_matrix_v(M, sets, intercept, n_own, feature_names, other, crit))
  }
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

# prune_dec_matrix() for V > 1 score columns (classification): M has V columns per set, component-major (column
# (k - 1) * V + v). A component's size is the variance of its score columns centred across the classes, summed
# (crit_centred_var()), and its share that over the same of the prediction. The decomposition also holds V and the
# criterion; its matrices' columns are named <component>#<v>.
prune_dec_matrix_v <- function(M, sets, intercept, n_own, feature_names, other, crit) {
  V <- crit$V
  K <- length(sets)
  M_own <- M[seq_len(n_own), , drop = FALSE]
  pred_own <- matrix(intercept, n_own, V, byrow = TRUE) + crit_sum_comp(M_own, seq_len(K), V)
  v_comp <- vapply(seq_len(K), function(k) crit_centred_var(crit, M_own[, crit_blk(k, V), drop = FALSE]), numeric(1))
  nonzero <- v_comp > 0
  share <- v_comp[nonzero] / crit_centred_var(crit, pred_own)
  names(share) <- prune_set_names(sets[nonzero], feature_names)
  M <- M[, crit_blk(which(nonzero), V), drop = FALSE]
  # (no component with variance, e.g. every family without a split: no columns)
  colnames(M) <- if (length(share)) paste0(rep(names(share), each = V), "#", seq_len(V)) else character(0)
  list(
    sets = sets,
    nonzero = nonzero,
    share = share,
    intercept = intercept,
    V = V,
    crit = crit,
    own = M[seq_len(n_own), , drop = FALSE],
    other = if (other) M[-seq_len(n_own), , drop = FALSE]
  )
}

# The columns of the components named `comps` in a decomposition's matrices (V per component when V > 1)
dec_cols <- function(dec, comps) {
  i <- match(comps, names(dec$share))
  if (is.null(dec$V)) i else crit_blk(i, dec$V)
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
# families, columns as names(dec$share), y = their outcome)) on the out-of-bag components instead. y_own: the outcome,
# or a classification criterion (then crit_fit_R(): one beta per component shared by the classes, one alpha per score
# column; D keeps the forest's intercepts).
prune_fit_rule <- function(rule, tau, dec, y_own, oob = NULL) {
  kept <- names(dec$share)[dec$share >= tau]
  if (inherits(y_own, "rpf_crit") && rule == "R") {
    ki <- match(kept, names(dec$share))
    coef <- if (is.null(oob)) {
      crit_fit_R(y_own, dec$own, ki, seq_len(y_own$n), start_alpha = dec$intercept)
    } else {
      crit_fit_R(oob$y, oob$F, ki, seq_len(oob$y$n), start_alpha = dec$intercept)
    }
    return(list(alpha = coef$alpha, beta = stats::setNames(coef$beta, kept)))
  }
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
# y: the numeric outcome, or the factor outcome of a classification forest (coded as rpf() codes it).
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
  Y <- if (object$mode == "classification") {
    preprocess_outcome(processed, params$loss)$outcomes
  } else {
    matrix(y[rows], ncol = 1)
  }
  set.seed(seed)
  fit <- rpf_impl(
    Y = Y,
    X = pre$predictors_matrix,
    mode = object$mode,
    params = params
  )
  new_rpf(
    fit = fit,
    blueprint = object$blueprint,
    mode = object$mode,
    factor_levels = pre$factor_levels,
    params = params
  )
}

# K-fold CV of every candidate (rule, tau); returns the CV table with the "min" and "1se" choices marked. With a
# classification criterion (tr$crit) the folds' forests are refitted on the factor outcome and the held-out rows are
# scored by the raw criterion: the table has cv_loss (the mean criterion) and se_loss instead of cv_rmse and se_mse.
prune_cv <- function(object, tr, cand, fold_id, seed, nthreads, verbose) {
  n <- length(tr$y)
  crit <- tr$crit
  fids <- sort(unique(fold_id))
  sse <- matrix(NA_real_, nrow(cand), length(fids))
  n_kept <- matrix(NA_integer_, nrow(cand), length(fids))
  n_k <- integer(length(fids))
  for (k in seq_along(fids)) {
    ho <- fold_id == fids[k]
    t0 <- proc.time()[["elapsed"]]
    # (a classification forest is refitted on its factor outcome, also where it is scored as a 0/1 response)
    fit_k <- prune_refit(object, tr$predictors, tr$outcome %||% tr$y, !ho, seed, nthreads)
    dec <- prune_decompose(fit_k, tr$predictors, !ho, other = ho, nthreads = nthreads, crit = crit)
    y_own <- if (is.null(crit)) tr$y[!ho] else crit_subset(crit, which(!ho))
    y_ho <- tr$y[ho]
    n_k[k] <- sum(ho)
    cache <- list()
    for (i in seq_len(nrow(cand))) {
      kept <- names(dec$share)[dec$share >= cand$tau[i]]
      key <- paste(cand$rule[i], paste(kept, collapse = ","))
      if (is.null(cache[[key]])) {
        rf <- prune_fit_rule(cand$rule[i], cand$tau[i], dec, y_own)
        cache[[key]] <- if (is.null(crit)) {
          sum((prune_predict_dec(rf, dec$other) - y_ho)^2)
        } else {
          # the held-out rows' criterion of the scores alpha + sum_S beta_S f_S (no correction: the whole forest)
          G <- crit_predict_R(rf, dec$other, match(names(rf$beta), names(dec$share)), crit$V)
          sum(crit_loss_raw(crit, which(ho), G))
        }
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
  if (!is.null(crit)) {
    out <- data.frame(
      rule = cand$rule,
      tau = cand$tau,
      cv_loss = rowSums(sse) / n,
      se_loss = se_mse,
      n_kept = rowMeans(n_kept),
      min = ch$min,
      one_se = ch$one_se,
      stringsAsFactors = FALSE
    )
    attr(out, "criterion") <- crit$type
    return(out)
  }

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

# Out-of-bag (OOB) choice among the candidates (rule, tau), from one pass over the families of `object` and no
# refits: every family is purified (mode 2, all orders, as purify() does) and its components are evaluated at the
# training rows; those outside its bootstrap sample give the OOB sums (prune_columns()). Returns prune_oob_score()'s
# table, the forest's decomposition on its training rows and the out-of-bag components for R's final fit (`oob_fit`).
# With cand = NULL (a given rule R), no table: only `dec` and `oob_fit`. keep_sums: also the OOB sums (`sums`), which
# the display's second pruning needs. fnames: the predictors' names.
prune_oob <- function(
  object,
  tr,
  fnames,
  cand,
  fold_id,
  nthreads,
  verbose,
  keep_sums = FALSE,
  call = rlang::caller_env()
) {
  t0 <- proc.time()[["elapsed"]]
  # `tr` are the training rows (checked by prune_components())
  pc <- prune_columns(object, fnames, nthreads, crit = tr$crit)
  o <- pc$o
  dec <- pc$dec
  rm(pc)
  out <- prune_oob_score(
    o,
    dec,
    fnames,
    prune_score_y(tr),
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

# The OOB table of every candidate from the OOB sums `o` (as prune_columns() returns them: those of purify_oob()) and
# the forest's decomposition on its training rows `dec` (prune_columns()). feature_names: the predictors' names; y_all: the training outcome;
# n_trees: the number of families. The kept components of each tau come from the forest's own shares (those of the
# final model). D is scored with the OOB intercept and components. R's (alpha, beta) are fitted as the final model's,
# by least squares on the OOB components (not the forest's own, which have fitted the training rows' noise), but on
# the rows outside each fold (prune_lstsq_qr()), and scored with the OOB components of the fold's rows. Rows without
# an OOB family are left out. The MSE is corrected for the OOB means' larger across-family variance (see below).
# Returns the table, `dec`, `oob_fit` (list(F = the OOB components, y = the outcome) of the rows with OOB families,
# for R's final fit by prune_fit_rule()) and `info` (rows used, mean number of OOB families per row). With
# cand = NULL (a given rule R) nothing is scored and the table is NULL. y_all may be a classification criterion:
# prune_oob_score_crit().
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
  if (inherits(y_all, "rpf_crit")) {
    return(prune_oob_score_crit(o, dec, feature_names, y_all, cand, fold_id, n_trees, deterministic, call))
  }
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

# prune_oob_score() for a classification criterion `crit` (criterion.R; V score columns per component,
# component-major): the same candidates, folds, rows and choice. D is scored with the OOB intercepts and components; R
# with crit_fit_R()'s coefficients fitted on the OOB components of the rows outside each fold, applied to the fold's
# OOB means and half means. The per-row criterion is corrected by crit_rows() (the rows without both halves imputed).
# The table has oob_loss (the mean corrected criterion; 0 where the correction overshoots, left out of the choice as for
# regression), se_loss, oob_loss_raw and se_loss_raw; the choice is prune_choose() on it, as on the MSE. oob_fit's y is
# the criterion of the rows with OOB families.
prune_oob_score_crit <- function(o, dec, feature_names, crit, cand, fold_id, n_trees, deterministic, call) {
  V <- crit$V
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
  comp <- names(dec$share)
  idx <- match(comp, prune_set_names(o$sets, feature_names))
  half <- function(s) {
    out <- matrix(0, n_use, length(comp) * V)
    ok <- which(!is.na(idx))
    if (length(ok)) {
      out[, crit_blk(ok, V)] <- s[use, crit_blk(idx[ok], V), drop = FALSE]
    }
    out
  }
  S0 <- half(o$sum0)
  S1 <- half(o$sum1)
  i0 <- o$intercept0[use, , drop = FALSE]
  i1 <- o$intercept1[use, , drop = FALSE]
  rm(o)
  F_oob <- (S0 + S1) / m
  rows <- which(use)
  cu <- crit_subset(crit, rows)
  if (is.null(cand)) {
    return(list(
      table = NULL,
      dec = dec,
      oob_fit = list(F = F_oob, y = cu),
      info = list(n_use = n_use, mean_m = mean(m))
    ))
  }
  I_oob <- (i0 + i1) / m
  F0 <- S0 / pmax(m0, 1)
  F1 <- S1 / pmax(m1, 1)
  I0 <- i0 / pmax(m0, 1)
  I1 <- i1 / pmax(m1, 1)
  rm(S0, S1)
  fid <- fold_id[use]
  fids <- sort(unique(fid))
  if (length(fids) < 2) {
    cli::cli_abort(
      "The training rows with out-of-bag families must lie in at least two folds of {.arg fold_id}.",
      call = call
    )
  }
  n_k <- vapply(fids, function(k) sum(fid == k), numeric(1))
  sse <- matrix(NA_real_, nrow(cand), length(fids))
  sse_corr <- sse
  n_kept <- integer(nrow(cand))
  cache <- list()
  for (i in seq_len(nrow(cand))) {
    ki <- which(dec$share >= cand$tau[i])
    key <- paste(cand$rule[i], paste(ki, collapse = ","))
    if (is.null(cache[[key]])) {
      if (cand$rule[i] == "D") {
        G <- I_oob + crit_sum_comp(F_oob, ki, V)
        G0 <- I0 + crit_sum_comp(F0, ki, V)
        G1 <- I1 + crit_sum_comp(F1, ki, V)
      } else {
        G <- G0 <- G1 <- matrix(0, n_use, V)
        for (k in fids) {
          te <- fid == k
          coef <- crit_fit_R(cu, F_oob, ki, which(!te), start_alpha = dec$intercept)
          sub <- function(Fm) Fm[te, , drop = FALSE]
          G[te, ] <- crit_predict_R(coef, sub(F_oob), ki, V)
          G0[te, ] <- crit_predict_R(coef, sub(F0), ki, V)
          G1[te, ] <- crit_predict_R(coef, sub(F1), ki, V)
        }
      }
      r <- crit_rows(crit, rows, G, G0, G1, m0, m1, n_trees, impute = TRUE)
      cache[[key]] <- list(
        sse = vapply(fids, function(k) sum(r$e[fid == k]), numeric(1)),
        sse_corr = vapply(fids, function(k) sum(r$corr[fid == k]), numeric(1))
      )
    }
    sse[i, ] <- cache[[key]]$sse
    sse_corr[i, ] <- cache[[key]]$sse_corr
    n_kept[i] <- length(ki)
  }
  se_of <- function(x) apply(sweep(x, 2, n_k, "/"), 1, stats::sd) / sqrt(length(fids))
  loss <- pmax(rowSums(sse_corr), 0) / n_use
  se <- se_of(sse_corr)
  ch <- prune_choose(sqrt(loss), se, cand, left_out = loss == 0)
  table <- data.frame(
    rule = cand$rule,
    tau = cand$tau,
    oob_loss = loss,
    se_loss = se,
    oob_loss_raw = rowSums(sse) / n_use,
    se_loss_raw = se_of(sse),
    n_kept = n_kept,
    min = ch$min,
    one_se = ch$one_se,
    stringsAsFactors = FALSE
  )
  attr(table, "n_rows") <- n_use
  attr(table, "criterion") <- crit$type
  list(
    table = table,
    dec = dec,
    oob_fit = list(F = F_oob, y = cu),
    info = list(n_use = n_use, mean_m = mean(m))
  )
}

# The root of the corrected mean error of an OOB table (a pruning's, or rpf()'s OOB path): its OOB RMSE, or for a
# classification criterion the root of the mean corrected criterion (as the tuner compares forests by it)
oob_table_rmse <- function(tab) {
  loss <- tab[["oob_loss", exact = TRUE]]
  if (is.null(loss)) tab[["oob_rmse", exact = TRUE]] else sqrt(pmax(loss, 0))
}

# The number of score columns of a pruned model (1 for regression and binary classification)
pruned_V <- function(object) {
  length(object$intercept)
}

# The kept components (scaled by beta) on new data, as a matrix with one named column per component (and score
# column: V per component, named as predict_components() names them, <component>__class:<level>)
prune_kept_matrix <- function(object, new_data, nthreads = 0L, call = rlang::caller_env()) {
  fit <- object$fit
  check_rpf_alive(fit, call = call)
  processed <- hardhat::forge(new_data, fit$blueprint)
  check_no_missing(processed$predictors, "{.arg new_data}", call = call)
  X <- preprocess_predictors_predict(fit, processed$predictors)
  prune_kept_matrix_x(object, X, nthreads)
}

# prune_kept_matrix() on the preprocessed predictors X (preprocess_predictors_predict())
prune_kept_matrix_x <- function(object, X, nthreads = 0L) {
  fit <- object$fit
  if (nthreads < 1) {
    nthreads <- as.integer(fit$params$nthreads)
  }
  # the display's coarsened components from their stored form, the others from the model's store
  V <- pruned_V(object)
  kept <- names(object$beta)
  coarse <- kept %in% names(object$display$coarse)
  M <- matrix(0, nrow(X), length(kept) * V)
  if (!all(coarse)) {
    M[, crit_blk(which(!coarse), V)] <- rpf_store_eval(
      prune_store_built(object),
      X,
      lapply(object$sets[!coarse], as.integer),
      nthreads
    )
  }
  for (k in which(coarse)) {
    b <- prune_display_built(object, kept[k])
    M[, crit_blk(k, V)] <- if (V == 1L) rpf_coarse_eval(b, X) else vapply(b, rpf_coarse_eval, numeric(nrow(X)), X)
  }
  M <- sweep(M, 2, rep(object$beta, each = V), "*")
  colnames(M) <- if (V == 1L || !length(kept)) {
    kept
  } else {
    paste0(rep(kept, each = V), "__class:", model_outcome_levels(object$fit))
  }
  M
}

# The scores (n x V) of a pruned model from its scaled kept columns M (prune_kept_matrix())
pruned_scores <- function(object, M) {
  V <- pruned_V(object)
  raw <- matrix(object$intercept, nrow(M), V, byrow = TRUE)
  if (ncol(M) > 0) {
    for (v in seq_len(V)) {
      raw[, v] <- raw[, v] + rowSums(M[, seq(v, ncol(M), by = V), drop = FALSE])
    }
  }
  raw
}

#' Predictions of a pruned Random Planted Forest
#'
#' The prediction of a pruned forest is `alpha + sum(beta_S * f_S(x))` over the kept components `f_S` of the forest,
#' with the intercept `alpha` and the coefficients `beta_S` (see [prune_components()]). For classification these are
#' the scores, one per score column, which `type` turns into probabilities or classes as [predict.rpf()] does.
#'
#' @param object `[rpf_pruned]`: A pruned forest from [prune_components()] or [rpf()].
#' @param new_data `[data.frame | matrix]`: New observations to predict.
#' @param type `[character(1)]`: As in [predict.rpf()]: `"numeric"` for regression (the default); for classification
#'   `"prob"` (the default), `"class"`, or `"numeric"` / `"link"` (the scores `alpha + sum(beta_S * f_S(x))`).
#' @param nthreads `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the `nthreads` the forest was fitted
#'   with.
#' @param ... Not currently used, but required for extensibility. Unknown arguments are an error.
#'
#' @return A [`tbl`][tibble::tibble] with the same number of rows as `new_data`: for regression column `.pred`; for
#'   classification the columns of [predict.rpf()] for `type`.
#' @seealso [prune_components()], [predict_components()].
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, n_families = 10, n_splits = 30, prune = FALSE)
#' pruned <- prune_components(fit, x, y, rule = "R", tau = 0.01)
#' predict(pruned, x[1:5, ])
predict.rpf_pruned <- function(
  object,
  new_data,
  type = ifelse(object$fit$mode == "regression", "numeric", "prob"),
  nthreads = NULL,
  ...
) {
  rlang::check_dots_empty()
  checkmate::assert_int(nthreads, lower = 1, null.ok = TRUE)
  type <- rlang::arg_match(type, c("numeric", "class", "prob", "link"))
  # 0 tells C++ to use the forest's own nthreads
  M <- prune_kept_matrix(object, new_data, as.integer(nthreads %||% 0L))
  out <- if (object$fit$mode == "regression") {
    if (type != "numeric") {
      cli::cli_warn(c(
        "Only {.code type = \"numeric\"} is supported for regression.",
        "i" = "Using it instead of {.val {type}}."
      ))
    }
    hardhat::spruce_numeric(object$intercept + rowSums(M))
  } else {
    predict_rpf_scores(object$fit, pruned_scores(object, M), type)
  }
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
  V <- pruned_V(object)
  # the kept components that use only `predictors` and have at most `max_interaction` predictors (and their V columns;
  # the sets index the fit's columns, in which hardhat::mold() puts the factors last)
  fit_names <- prune_fit_names(object$fit)
  feats <- lapply(object$sets, function(s) fit_names[s])
  use <- vapply(feats, function(f) length(f) <= max_interaction && all(f %in% predictors), logical(1))
  usec <- rep(use, each = V)

  ret <- list(
    m = data.table::as.data.table(M[, usec, drop = FALSE]),
    intercept = object$intercept,
    # a copy: setDT() would fail on a matrix and turn the caller's data.frame into a data.table
    x = data.table::as.data.table(new_data)[, predictors, with = FALSE]
  )
  # multiclass: per score column, as predict_components() of the forest names them
  if (V > 1L) {
    lev <- model_outcome_levels(object$fit)
    ret$intercept <- stats::setNames(object$intercept, lev)
    ret$target_levels <- lev
  }
  if (!all(use)) {
    ret$remainder <- if (V == 1L) {
      rowSums(M[, !usec, drop = FALSE])
    } else {
      R <- pruned_scores(list(intercept = rep(0, V)), M[, !usec, drop = FALSE])
      colnames(R) <- model_outcome_levels(object$fit)
      data.table::as.data.table(R)
    }
  }
  class(ret) <- c("glex", "rpf_components", class(ret))
  ret
}

#' Print a pruned Random Planted Forest
#'
#' Shows the model: its intercept and its components, each with its importance, the component's share of the variance
#' of the model's prediction on the training rows (`components$model_share`). How the components were chosen (the
#' forest's families and splits, given or chosen out of bag, the rule and threshold, what was left out, the heat-map
#' forms) follows briefly at the end.
#'
#' @param x `[rpf_pruned]`: A pruned forest from [prune_components()] or [rpf()].
#' @param n_show `[integer(1): 20]`: The largest number of components to list.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly: `x`.
#' @seealso [prune_components()].
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
#' fit <- rpf(x, y, max_interaction = 2, n_families = 10, n_splits = 30, prune = FALSE)
#' prune_components(fit, x, y, rule = "R", tau = 0.01)
print.rpf_pruned <- function(x, n_show = 20, ...) {
  cat(format(x, n_show = n_show, ...), sep = "\n")
  invisible(x)
}

#' @export
format.rpf_pruned <- function(x, n_show = 20, ...) {
  checkmate::assert_int(n_show, lower = 0)
  cmp <- x$components
  kept <- cmp[cmp$kept, ]
  V <- pruned_V(x)
  K <- nrow(kept)
  classif <- identical(x$fit$mode, "classification")

  # the model: what it predicts, its intercept and its components with their importance
  of <- if (K == 1) "1 component" else paste(K, "components")
  lines <- if (classif) {
    lev <- levels(x$fit$blueprint$ptypes$outcomes[[1]])
    sprintf(
      "Classification (%s loss), %d classes: %s. Scores = intercept + the sum of %s.",
      x$fit$params$loss,
      length(lev),
      paste(lev, collapse = ", "),
      of
    )
  } else {
    sprintf("Regression. Prediction = intercept + the sum of %s.", of)
  }
  lines <- c(
    lines,
    if (V > 1L) {
      sprintf(
        "Intercept per class: %s",
        paste(model_outcome_levels(x$fit), format(x$intercept, digits = 6), collapse = ", ")
      )
    } else {
      sprintf("Intercept: %s", format(x$intercept, digits = 6))
    }
  )
  if (K > 0 && n_show > 0) {
    # importance: the component's share of the model's prediction (models from before model_share: of the forest's)
    imp <- kept$model_share
    of_what <- "share of the variance of the prediction on the training rows"
    if (is.null(imp) || anyNA(imp)) {
      imp <- kept$share
      of_what <- "share of the variance of the forest's prediction on the training rows"
    } else if (V > 1L) {
      of_what <- "share of the variance of the scores (centred across the classes) on the training rows"
    }
    o <- order(-imp)
    show <- utils::head(o, n_show)
    tab <- data.frame(component = kept$component[show], importance = vapply(imp[show], format_share, ""))
    lines <- c(
      lines,
      "",
      sprintf("Components (largest first), importance = %s:", of_what),
      utils::capture.output(print(tab, row.names = FALSE, right = FALSE))
    )
    if (K > n_show) {
      lines <- c(lines, sprintf("... and %d more", K - n_show))
    }
  }

  lines <- c(lines, "", "How it was chosen:", paste0("  ", format_pruned_choice(x)))
  # format into lines so print() writes to stdout, not cli's message stream
  c(cli::cli_format_method(cli::cli_rule(left = "Random Planted Forest model")), lines)
}

# The lines of format.rpf_pruned() on how the components were chosen: the rule and threshold on the forest's
# components and how they were chosen, what was left out, the heat-map forms (the display) and the tuning
format_pruned_choice <- function(x) {
  d <- x$display
  # without the display, the one pruning; with it, the first (on the forest's components), then the second
  m <- if (is.null(d)) x else d$first
  cmp <- m$components
  # the forest's size: its families and splits, given or chosen out of bag (rpf(splits = NULL))
  n_fmt <- function(v) format(v, big.mark = ",", scientific = FALSE)
  so <- x$fit$splits_oob
  lines <- sprintf(
    "The forest: %s tree famil%s of %s splits%s.",
    n_fmt(x$fit$params$ntrees),
    if (identical(as.integer(x$fit$params$ntrees), 1L)) "y" else "ies",
    n_fmt(x$fit$params$splits),
    if (is.null(so)) {
      ""
    } else {
      sprintf(
        " (chosen out of bag; growing stopped at %s%s)",
        n_fmt(so$stopped_at),
        if (isTRUE(so$capped)) ", the cap, with the error still falling" else ""
      )
    }
  )
  rule <- c(
    D = "as they are",
    R = if (identical(m$select_by, "oob")) {
      "each with one coefficient, fitted out of bag"
    } else {
      "each with one coefficient"
    }
  )[[m$rule]]
  lines <- c(
    lines,
    sprintf(
      paste(
        "From the forest's %d nonzero components (purified): those with a share of at least %s of the forest's",
        "Var(prediction), %s (rule %s)."
      ),
      nrow(cmp),
      format_share(m$tau),
      rule,
      m$rule
    )
  )
  crit <- pruned_crit_name(x)
  lines <- c(
    lines,
    if (!is.null(m$oob)) {
      sprintf(
        "Rule and threshold chosen by out-of-bag error (%s; %d folds, seed %d) among %d candidates.",
        if (m$select == "min") paste("lowest OOB", crit) else "1-SE rule",
        length(unique(m$fold_id)),
        m$seed,
        nrow(m$oob)
      )
    } else if (!is.null(m$cv)) {
      sprintf(
        "Rule and threshold chosen by %d-fold CV (%s, seed %d) among %d candidates.",
        length(unique(m$fold_id)),
        if (m$select == "min") paste("lowest CV", crit) else "1-SE rule",
        m$seed,
        nrow(m$cv)
      )
    } else {
      "Rule and threshold given."
    }
  )
  out <- cmp[!cmp$kept, ]
  if (nrow(out) > 0) {
    big <- out[which.max(out$share), ]
    lines <- c(
      lines,
      sprintf(
        "Left out: %d component%s, shares together %s (largest %s, %s).",
        nrow(out),
        if (nrow(out) == 1) "" else "s",
        format_share(sum(out$share)),
        big$component,
        format_share(big$share)
      )
    )
  }
  if (!is.null(d)) {
    co <- d$coarse[vapply(d$coarse, function(z) isTRUE(z$kept), logical(1))]
    if (d$second) {
      n_out <- sum(!x$components$kept)
      lines <- c(
        lines,
        sprintf(
          paste(
            "Components of three or more predictors in at most %d heat maps each, then rule and threshold chosen again on",
            "these (%s, tau = %s; %s left out)."
          ),
          d$max_panels,
          x$rule,
          format_share(x$tau),
          if (n_out == 0) "none" else n_out
        )
      )
    }
    if (length(co) > 0) {
      forms <- vapply(
        names(co),
        function(nm) {
          z <- co[[nm]]
          sprintf(
            "%s: %s (%s of its variance)",
            nm,
            paste0(z$C, " in ", lengths(z$cuts) + 1, " intervals", collapse = ", "),
            format_share(z$var_kept)
          )
        },
        ""
      )
      lines <- c(lines, paste0("Heat-map forms: ", forms))
    }
  }
  if (!is.null(x$tune)) {
    st <- x$tune$settings
    lines <- c(
      lines,
      sprintf(
        paste(
          "The forest was tuned out of bag (rpf(tune = )): %d forests; n_try_split_points = %d, share_candidates = %s,",
          "min_leaf_rows = %d, n_splits = %d."
        ),
        nrow(x$tune$search),
        st$n_try_split_points,
        format(st$share_candidates),
        st$min_leaf_rows,
        st$n_splits
      )
    )
  }
  lines
}

# The name of what a pruned model was chosen by: the RMSE, or the classification criterion
pruned_crit_name <- function(x) {
  switch(x$criterion %||% "mse", mse = "RMSE", brier = "Brier score", logloss = "log-loss")
}

format_share <- function(s) {
  if (s == 0) {
    return("0")
  }
  paste0(trimws(formatC(signif(100 * s, 2), digits = 2, format = "fg")), "%")
}
