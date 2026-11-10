#' Random Planted Forest
#'
#' Fits a random planted forest (Hiabu, Mammen and Meyer, 2020) for regression
#' or classification. The model is a sum of trees that each split on at most
#' `max_interaction` predictors, so it decomposes into main effects and
#' interaction components; see [predict_components()].
#'
#' By default, `rpf()` returns the pruned model that can be shown, in one call:
#' 1. it grows 200 tree families of up to three-way interactions, with the
#'    number of splits chosen out of bag while they grow (`splits = NULL`);
#' 2. it prunes the forest's components by [prune_components()] with its
#'    defaults (`prune = NULL`): the sparse rule, chosen out of bag;
#' 3. it coarsens every kept component of three or more predictors to at most
#'    `max_panels = 10` heat maps and prunes again.
#'
#' The result, of class `rpf_pruned`, has methods for
#' [`predict()`][predict.rpf_pruned], [predict_components()] and `print()`.
#' `prune = FALSE` returns the forest (class `rpf`) instead, and
#' `splits = <number>` fits a fixed number of splits.
#'
#' @param x,y `[data.frame | matrix | recipe]`, `[numeric | factor]`:
#'   Predictors and outcome for the x/y interface.
#'   `x` holds the predictors, or is a [`recipe`][recipes::recipe] (then supply
#'   `data` instead of `y`). The predictors need unique names that do not
#'   contain `:`, which joins predictor names in component names such as `x1:x2`.
#'   `y` is the outcome vector: __numeric__ for regression, __factor__ for
#'   classification. A numeric outcome must be finite. A factor outcome needs at
#'   least two classes; levels that do not occur are dropped with a warning.
#' @param ... Not currently used, but required for extensibility. Unknown
#'   (e.g. misspelled) arguments are an error.
#' @param formula `[formula]`: A formula with the outcome on the left-hand side and the
#'   predictors on the right-hand side, e.g. `y ~ x1 + x2`. Interaction terms
#'   such as `x1:x2` or `x1 * x2` are an error: they would add a predictor named
#'   `x1:x2`, and `:` is how components are named. The forest finds interactions
#'   itself, up to `max_interaction`, so list the predictors as `y ~ x1 + x2`.
#'   To fit a product as a predictor of its own, write `I(x1 * x2)`, or add it to
#'   `data` as a column with its own name. Transformed terms of one predictor,
#'   such as `log(x1)`, are allowed.
#' @param data `[data.frame]`: Data containing the predictors and the outcome,
#'   used with `formula` or a recipe `x`.
#' @param max_interaction `[integer(1): 3]`: Maximum number of predictors a single tree may
#'   split on. `1` fits main effects only (an additive model), `2` adds
#'   pairwise interactions, `3` (the default) three-way interactions, and so on. `0` uses all predictors and values above
#'   the number of predictors are reduced to it.
#' @param ntrees `[integer(1): 200]`: Number of tree families, i.e. the size of the forest.
#'   Each family is grown on a bootstrap sample and predictions are averaged.
#' @param splits `[integer(1) | NULL: NULL]`: Number of split rounds per tree family. The main tuning
#'   parameter: more splits fit more complex functions. A round splits the best of its
#'   candidates; near saturation, when every candidate's leaf is too small to split, a
#'   round can make no split. `NULL` (the default) chooses it out of bag while the forest
#'   grows (with `deterministic = TRUE`, which has no out-of-bag rows, `NULL` means 30). See
#'   "Choosing `splits` out of bag" in Details.
#' @param split_structure `[character(1): "leaves"]`: How split candidates are formed and
#'   sampled; one of `"leaves"`, `"hist"`, `"cur_trees_1"`, `"cur_trees_2"`,
#'   or `"res_trees"`. See Details.
#' @param split_try `[integer(1) | NULL: NULL]`: Number of thresholds evaluated per split
#'   candidate. `NULL` means 3 for `split_structure = "leaves"` and 10 for the
#'   other structures.
#' @param t_try `[numeric(1) | NULL: NULL]`: Proportion in `(0, 1]` of all possible splits sampled
#'   as candidates in each round, capped at `max_candidates`. `NULL` means 0.1
#'   for `split_structure = "leaves"` and 0.4 for the other structures. See
#'   "Choosing parameters" in Details.
#' @param max_candidates `[integer(1): 50]`: Maximum number of split candidates per round.
#' @param split_decay_rate `[numeric(1): 0.1]`: Down-weights possible splits that were
#'   sampled as candidates but not chosen. Each such round ages a split by one,
#'   choosing it resets its age to zero, and splits are sampled with weight
#'   `exp(-split_decay_rate * age)`. `0` samples uniformly.
#' @param delete_leaves `[logical(1): TRUE]`: Whether to delete a parent leaf when
#'   splitting along an existing dimension.
#' @param min_leaf_rows `[integer(1): 1]`: Minimum number of rows of the
#'   bootstrap sample on each side of every split, so in every leaf. Leaves
#'   without a threshold that keeps this many rows on each side are not drawn
#'   as split candidates. Only for `split_structure = "leaves"` (regression and
#'   classification); `1` (no limit) otherwise.
#' @param interaction_constraints `[list | NULL: NULL]`: Feature sets that restrict which
#'   interactions the forest may form, e.g. `list(c("x1", "x2"), c("x3", "x4"))`.
#'   Each set is given by predictor names or by column indices: with `x`, the
#'   columns of `x`; with a formula, the predictors in the order the formula
#'   names them (`y ~ .`: the order of the columns of `data`); with a recipe,
#'   names only (its steps decide the columns). A tree on the predictors `T` is
#'   grown only if `T` is a subset of one of the sets, as with xgboost's
#'   interaction constraints. The subsets of a set stay allowed, so the forest
#'   keeps the main effects it grows a pair from. A predictor in no set is not
#'   used. Applies together with `max_interaction`: set `max_interaction = 0` to
#'   let the sets alone decide. `NULL` means no constraint. One use is to refit
#'   the forest on the effects a first fit found, e.g. after dropping its small
#'   components.
#' @param loss `[character(1) | NULL: NULL]`: Loss function. `NULL` uses `"L2"` for
#'   regression, the only supported loss, and `"exponential"` for classification
#'   (with more than two classes the multiclass exponential loss SAMME).
#'   Classification also supports `"logit"`, which gives similar results but is
#'   slower, and `"L1"` and `"L2"`, which fit class indicators directly and do
#'   not yield proper probability estimates.
#' @param delta `[numeric(1): 0.001]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[delta, 1 - delta]` when computing the
#'   loss of a split. Must be less than 0.5, from where every proportion is
#'   truncated to the same value, and positive for `"logit"`: with `delta = 0`,
#'   nodes containing a single class have infinite loss and their splits are
#'   always rejected.
#' @param epsilon `[numeric(1): 0.1]`: Only used if `loss` is `"logit"` or `"exponential"`.
#'   Class proportions are truncated to `[epsilon, 1 - epsilon]` when computing
#'   the fit in a leaf. Unlike `delta`, this caps the size of individual leaf
#'   updates and acts as regularization: smaller values permit larger jumps on
#'   the link scale. Must be greater than 0 and less than 0.5: at 0 the leaf
#'   values of pure leaves are infinite, and from 0.5 on every leaf value is 0.
#' @param prune `[logical(1) | NULL: NULL]`: Whether to return the pruned model
#'   instead of the forest: the forest is pruned by [prune_components()] with
#'   its defaults on the training rows, and the result is its `rpf_pruned`
#'   model. `NULL` (the default) means `TRUE`, and `FALSE` with
#'   `deterministic = TRUE`, where `TRUE` is an error (pruning needs out-of-bag
#'   rows). `FALSE` returns the forest. See "Pruning" in Details.
#' @param max_panels `[numeric(1): 10]`: Used only when pruning: the
#'   most heat maps a kept component of three or more predictors may need; such
#'   a component is coarsened to at most `max_panels` heat maps and the model is
#'   pruned again (the `max_panels` of [prune_components()]). A whole number of
#'   at least 1, or `Inf` for no coarsening.
#' @param tune `[logical(1) | character(1) | rpf_tune_control: FALSE]`: Whether
#'   to tune the greediness (`split_try`, `t_try`) and `min_leaf_rows` out of
#'   bag (not with `deterministic = TRUE`, at least 10 tree families). `TRUE` runs the full search ([tune_control()]), `"quick"` the
#'   two-forest preset, and a [tune_control()] object your own. Every forest of
#'   the search has `splits` chosen out of bag and is pruned, so `tune` implies
#'   `splits = NULL` and `prune = TRUE`, and giving `split_try`, `t_try` or
#'   `min_leaf_rows` with it is an error. See "Tuning out of bag" in Details.
#' @param purify `[logical(1): FALSE]`: Whether to purify the forest after fitting, which
#'   [predict_components()] requires, with the defaults of [purify()]. For other
#'   settings, call [purify()] after fitting instead. Not used with
#'   `prune = TRUE`: the pruned model stores the components it keeps.
#' @param nthreads `[integer(1): 1]`: Number of threads for fitting. Also the default for
#'   [predict()][predict.rpf()] and [purify()].
#' @param export_forest `[logical(1): FALSE]`: Whether to store the flattened forest in
#'   the returned object as `$forest` (with `prune = TRUE`, in the forest:
#'   `$fit$forest`). If `FALSE`, it is `NULL`, which saves memory.
#' @param deterministic `[logical(1): FALSE]`: Whether to fit without randomness: no
#'   bootstrap, the first `max_candidates` possible splits as candidates
#'   (`t_try` and `split_decay_rate` are ignored) and evenly spaced thresholds. The result does not
#'   depend on the seed, and all tree families are identical, so use
#'   `ntrees = 1`. There are no out-of-bag rows, so `splits = NULL` means 30
#'   and `prune = NULL` means `FALSE`. Mainly useful for testing and debugging.
#' @param keep_level_order `[logical(1) | character: FALSE]`: Which factor
#'   predictors keep the order of their levels instead of being ordered by the
#'   outcome. `FALSE` orders every factor by the outcome, `TRUE` keeps the level
#'   order of every factor, and a character vector keeps it for the named
#'   factor predictors only. See "Factor predictors" below.
#' @param oob_path `[logical(1) | integer: FALSE]`: Whether to track the
#'   out-of-bag (OOB) error along the splits (not with
#'   `deterministic = TRUE`). `TRUE` evaluates it at the checkpoints of
#'   [oob_path_checkpoints()] for `splits`, four per doubling; an increasing
#'   integer vector in `1..splits` gives the checkpoints. The fitted forest is
#'   the same as without. See Details. With `splits = NULL` the path of the
#'   checkpoints grown is stored either way, and only `TRUE` or `FALSE` is
#'   allowed.
#'
#' @return With the defaults, an object of class `rpf_pruned` (see below). With
#'   `prune = FALSE`, an object of class `"rpf"` with model object contained in `$fit`. With
#'   `oob_path`, `$oob_path` holds the OOB path: a data frame with one row per
#'   checkpoint, `splits`, the corrected and the plain OOB root mean squared
#'   error (`oob_rmse`, `oob_rmse_raw`), and the standard errors of the two mean
#'   squared errors (`se_mse`, `se_mse_raw`). Its attribute `n_rows` is the
#'   number of training rows with an OOB tree family. For classification with
#'   `loss = "logit"` or `"exponential"`, or more than two classes, the columns
#'   are the mean corrected and plain criterion (`oob_loss`, `oob_loss_raw`) and
#'   their standard errors (`se_loss`, `se_loss_raw`), and the attribute
#'   `criterion` names it (see "Classification" in Details). With `splits = NULL`,
#'   `$params$splits` is the chosen number of splits, `$oob_path` the path of
#'   the checkpoints grown, and `$splits_oob` a list with the chosen `splits`,
#'   `stopped_at` (the number of splits where growing stopped) and `capped`
#'   (whether it stopped at the cap while the error was still falling).
#'   When pruning (the default), an object of class `rpf_pruned` (see
#'   [prune_components()]), whose element `fit` is the forest (not purified),
#'   with the elements above: the pruned model's OOB path is `$fit$oob_path`, and
#'   with `splits = NULL` the chosen number of splits `$fit$params$splits` and
#'   the stop `$fit$splits_oob`.
#'   With `tune`, the `rpf_pruned` model of the chosen forest, with an element
#'   `tune`: `search`, the search table (one row per forest grown, in the order
#'   grown: `step`, `split_try`, `t_try`, `min_leaf_rows`, the chosen `splits`,
#'   `stopped_at`, `capped`, the forest's corrected OOB RMSE `oob_rmse_forest`
#'   (as the pruning scores it: its rule `"D"` at `tau = 0`, which keeps every
#'   nonzero component; the OOB path's `oob_rmse` at the chosen `splits` leaves
#'   the rows without OOB families of both halves uncorrected and is higher),
#'   its pruned model's `rule`, `tau`, `n_kept` and corrected OOB RMSE
#'   `oob_rmse` (the score; for classification both are the square root of the
#'   mean corrected criterion), `seconds`, whether each step chose it
#'   (`chosen_step1`, `chosen_step2`, `chosen_step3`; `NA` if not a candidate of
#'   the step) and `final`); `steps`, each step's candidates and choice;
#'   `settings`, the chosen `split_try`, `t_try`, `min_leaf_rows` and `splits`;
#'   `control`, the [tune_control()] used; and `seconds`.
#' @export
#' @importFrom methods new
#' @importFrom hardhat mold
#' @importFrom hardhat default_xy_blueprint
#' @importFrom hardhat default_formula_blueprint
#' @importFrom hardhat default_recipe_blueprint
#'
#' @details
#' ## Choosing parameters
#'
#' Start with `max_interaction` and `splits`. `max_interaction` sets which
#' effects the model can represent: `1` for an additive model, `2` (the
#' default) to add pairwise interactions, which can still be plotted as
#' heatmaps. Higher values are more flexible, but slower and the components are
#' harder to interpret. `splits` sets how closely the model
#' follows the data and is the parameter to tune. With a higher
#' `max_interaction`, splits are spread over more possible components, so tune
#' both together.
#'
#' `ntrees` averages over bootstrap samples to reduce variance and fitting and
#' prediction time grow linearly with it. The default is lower than the
#' usual 500 trees of a random forest because each tree family is itself a
#' complete model, a sum of trees covering all components, rather than a
#' single tree. Averaging therefore stops paying off sooner, usually after a
#' few dozen families. The split search parameters
#' (`split_structure`, `split_try`, `t_try`, `max_candidates`,
#' `split_decay_rate`) trade speed for a more thorough search.
#'
#' For classification, `"exponential"` (the default) and `"logit"` fit on the
#' link scale and are transformed to probabilities, `"exponential"` being faster.
#' With more than two classes `"exponential"` is the multiclass exponential
#' loss SAMME (Zhu, Zou, Rosset and Hastie, 2009): with one score per class,
#' summing to 0, a row's weight is `exp(-g)`, `g` the score of its class, and
#' the probabilities are the softmax of the scores.
#' `"L1"` and `"L2"` fit class indicators directly and their probabilities are
#' truncated to `[0, 1]`.
#'
#' ## Factor predictors
#'
#' Factor predictors, and character predictors, which are converted to
#' factors, are coded as integers `1, 2, ...` in an order of their levels, and
#' splits are thresholds on these codes. A split therefore separates the levels
#' below a cut from those above it in that order.
#'
#' By default, the levels are ordered by the outcome, so that the levels
#' most alike in the outcome are neighbours (Wright and König, 2019):
#' - regression: by the mean outcome of each level;
#' - binary classification: by the share of the second outcome level;
#' - more than two classes: by the first principal component of the class
#'   proportions of each level (Coppersmith et al., 1999).
#'
#' This applies to ordered factors too. With `keep_level_order`, the chosen
#' factors are instead coded in the order of `levels()`, as given, ordered or
#' not, so that every split separates neighbouring levels in that order.
#'
#' The coding is stored in the fit as `$factor_levels`, the levels of each
#' factor in the order of their codes. [predict()][predict.rpf()],
#' [predict_components()] and [purify()] code new data by these levels, by
#' name, so the order of the levels in `new_data` does not matter.
#'
#' ## `split_structure`
#'
#' In each round, a `t_try` fraction of all possible splits (capped at
#' `max_candidates`) is drawn as candidates with weights
#' `exp(-split_decay_rate * age)`. `split_structure` defines what a candidate is
#' and how its thresholds are evaluated.
#'
#' - `leaves`: Split candidates are (leaf, split-dimension) pairs. For each sampled
#'   candidate, `split_try` thresholds are drawn uniformly without replacement from the
#'   distinct values of that predictor among the leaf's rows, leaving out the smallest
#'   and, with `min_leaf_rows`, those that leave fewer rows on either side (all of them
#'   if there are fewer), and evaluated to choose the best split.
#'
#' - `hist`: As `leaves`, but thresholds are drawn from boundaries of
#'   quantile bins computed once per predictor, which makes evaluating them
#'   faster on large data.
#'
#' - `cur_trees_1`: Split candidates are (current-tree, split-dimension) pairs. For each
#'   sampled candidate, perform `split_try` evaluations. Each evaluation samples a leaf
#'   from the set of valid current trees (with probability proportional to its number of
#'   available thresholds) and then uniformly samples a single threshold within that leaf.
#'
#' - `cur_trees_2`: Split candidates are (current-tree, split-dimension) pairs. For each
#'   sampled candidate, iterate through every
#'   valid leaf. Within each leaf, sample `split_try` thresholds uniformly and
#'   evaluate them.
#'
#' - `res_trees`: Split candidates are resulting trees. For each sampled candidate, run
#'   `split_try` evaluations by sampling a (split-dimension, leaf) pair from all valid
#'   pairs (with probability proportional to its number of available thresholds), then
#'   uniformly sampling one threshold within that pair.
#'
#' ## The draw of candidates for `leaves`
#'
#' With `split_structure = "leaves"` (regression and classification), the
#' global option `randomPlantedForest.exact_draw` chooses how the candidates
#' are drawn.
#' With `TRUE` (the default), every possible split gets a random number in
#' every round. With `FALSE`, the candidates are drawn from a Fenwick tree
#' over the weights, kept up to date between rounds, which costs a few random
#' numbers per candidate instead. Both draw from the same distribution, but
#' with other random numbers, so the same seed gives another (equally good)
#' forest. Ties between equally good candidates may go to different candidates
#' under the two draws, so with many splits, near saturation, their fits differ
#' in distribution, though not in accuracy. With `min_leaf_rows`, both leave
#' out the same candidates: those
#' whose leaf has no threshold with that many rows on each side. The fast draw
#' pays off on large data with many splits. The option is read when the forest
#' is fitted and stored in `params$exact_draw`; it does not affect the other
#' split structures.
#'
#' ```r
#' options(randomPlantedForest.exact_draw = FALSE)
#' ```
#'
#' ## Reproducibility
#'
#' With the same seed, `rpf()` gives the same forest on the same platform, also
#' with more threads. Another platform can give another forest for the same
#' seed, because the integer sampling of the C++ standard library
#' (`std::uniform_int_distribution`) differs between implementations, such as
#' libc++ and libstdc++.
#'
#' ## The out-of-bag path
#'
#' With the same seed, the forest grown with `k` splits is exactly the first
#' `k` splits of every tree family of a forest grown with more: nothing in the
#' fit depends on the total number of splits. Each family is grown on a
#' bootstrap sample, so every training row is out of bag for the families that
#' did not draw it (about 37% of them). With `oob_path`, each family's
#' prediction at its OOB rows is updated split by split while it grows, and
#' copied at each checkpoint. So one fit gives the OOB error of every forest
#' with fewer splits, for example to choose `splits` without cross-validation.
#' This costs `2 * nrow(x) * length(checkpoints)` doubles and a little time;
#' the fit itself is unchanged.
#'
#' The OOB prediction of a row is the mean over its OOB families, and the plain
#' error (`oob_rmse_raw`) is the root mean squared error of these predictions
#' over the rows with an OOB family. That mean averages fewer families than the
#' forest, so the plain error overstates the forest's error by the families'
#' extra variance, more so with more splits. The corrected error (`oob_rmse`)
#' subtracts this excess, `(1 / m_i - 1 / ntrees) * v_i` for row `i` with `m_i`
#' OOB families, where `v_i` is the variance across families estimated from the
#' OOB means over the families of even and of odd index,
#' `v_i = (g_even - g_odd)^2 / (1 / m_even + 1 / m_odd)`; this is also the
#' correction of [prune_components()]' OOB choice. A row whose OOB families
#' are all even or all odd is not corrected here (`prune_components()` gives
#' it the mean estimate of the other rows): about 2% of the rows at 20
#' families, 18% at 10 and over half at 4 (none at 50), so with few families
#' the corrected error is biased upwards, and above the corrected error
#' `prune_components()` gives the same forest (rule `"D"` at `tau = 0`).
#' Factor levels are ordered by the outcome of all training rows: by the mean
#' outcome (regression), by the share of the second class (binary
#' classification) or by the first principal component of the class
#' proportions (three or more classes). So a row's own outcome enters the
#' ordering its OOB families split on, and the OOB error is somewhat optimistic
#' with factor predictors, most so with many levels.
#'
#' ## Choosing `splits` out of bag
#'
#' With `splits = NULL` (not deterministic), the number of splits
#' is chosen while the forest grows, by the corrected OOB error of the
#' out-of-bag path. The families grow to checkpoints, four per doubling of the
#' splits ([oob_path_checkpoints()] up to 25,600), and the forest's corrected
#' error is followed checkpoint by checkpoint. Growing stops once the lowest
#' error is one doubling of the splits back (for example, no new minimum from
#' 100 to 200 splits), or at 25,600 splits, with a warning if the error was
#' still falling there. The result is the forest at the lowest error (ties go
#' to fewer splits; a corrected error at or below 0, a correction that
#' overshot, is not chosen unless all are): the forest that `rpf()` with that
#' many `splits` and the same seed returns. On 52 data sets (44 synthetic, 8 splits of California
#' housing) this chose the same `splits` as the lowest error along a fixed
#' budget of 3,200 splits on 50 of them, and a later, lower minimum on the
#' other two, while growing 0.72 of the splits.
#'
#' The families grow one at a time, `nthreads` together, each to the stop
#' estimated from the families grown before it. A grown family keeps only its
#' OOB predictions at the checkpoints and a record of its splits (40 bytes a
#' split), from which the chosen forest is rebuilt, without a refit. A family
#' that reaches the estimate waits there with its growing state while the
#' waiting families fit in `getOption("randomPlantedForest.adaptive_pause_gb",
#' 1)` GB, so that it can grow on cheaply if the stop comes later; otherwise
#' it is freed, and grown again from its seed if needed. The option trades
#' memory for time; the fit depends neither on it nor on `nthreads`. Growing
#' all families together would hold every family's growing state until the
#' stop: about 18 GB of resident memory (RSS) on California housing (16,512
#' rows) at 200 families and `max_interaction = 3`, where the fit one family at
#' a time peaks at 3.4 GB of RSS (2.3 GB of physical footprint; the adaptive
#' fit alone, measured when it was added: `split_try = 3`, `t_try = 0.1`, the
#' exact draw, 8 threads).
#'
#' ## Pruning
#'
#' With `prune = TRUE`, `rpf()` fits the forest, then prunes its components
#' with [prune_components()] at that function's defaults, on the
#' training rows: the rule and the threshold are chosen out of bag (by
#' cross-validation below 10 tree families, with a message), the sparse choice
#' `select = "1se"`, and every kept component of three or more predictors is
#' shown in at most `max_panels` heat maps before the model is pruned again. The
#' result is the pruned model, of class `rpf_pruned`, with methods for
#' [`predict()`][predict.rpf_pruned], [predict_components()] and `print()`. It
#' is the same, bit for bit, as
#' `prune_components(rpf(..., prune = FALSE), data)` from the same random
#' number state, for every interface. For other settings of the pruning (`rule`,
#' `tau`, `select`, `select_by`, folds or seed), fit with `prune = FALSE` and call
#' [prune_components()]. Together with `splits = NULL`, one call chooses the
#' number of splits out of bag, prunes and returns the model to show. Pruning
#' needs out-of-bag rows, so `prune = TRUE` is an error with
#' `deterministic = TRUE`.
#'
#' ## Classification
#'
#' A classification forest goes through the same steps (the OOB path,
#' `splits = NULL`, `prune`, `max_panels` and `tune`) with its score columns:
#' one for a binary outcome, one per class otherwise (`K - 1` for
#' `loss = "logit"`: the first level is the reference). Every out-of-bag choice
#' is scored by the loss's own criterion, corrected for the OOB family count as
#' above:
#' - `loss = "L1"` or `"L2"`: the squared error of the raw scores summed over
#'   the columns (the Brier score, before the scores are truncated to
#'   probabilities). With two classes this is the regression criterion of the
#'   0/1 outcome, and the choices are the regression ones, bit for bit.
#' - `loss = "logit"` or `"exponential"`: the log-loss of the probabilities
#'   that `predict(type = "prob")` returns. The correction is its second-order
#'   form, `(1 / m_i - 1 / ntrees) * d' H d / (2 * (1 / m_even + 1 / m_odd))`
#'   with `d = g_even - g_odd` and `H = diag(p) - p p'` over the free score
#'   columns (for the squared error `H = 2 I`, the correction above).
#'
#' The pruning keeps or drops a component for all classes, sizes it by its
#' variance summed over the score columns centred across the classes, and fits
#' one coefficient per kept component shared by the classes (see
#' [prune_components()]). The criterion was checked against the test log-loss
#' on simulated data; the defaults were chosen on regression and carried over.
#'
#' ## Tuning out of bag
#'
#' With `tune`, `rpf()` searches the greediness of the split search
#' (`split_try` and `t_try`, on the ladder of [greediness_ladder()]) and the
#' minimum leaf size `min_leaf_rows`, without cross-validation. Every forest of
#' the search is grown with `splits` chosen out of bag and pruned by
#' [prune_components()] with its defaults (without the display), and forests are
#' compared by the corrected out-of-bag RMSE of their pruned model (for
#' classification by the root of its mean corrected criterion); the lowest
#' wins. A pruned model at 0 (every candidate's correction overshot, see
#' [prune_components()]) is not chosen unless all are. The search is
#' coordinate-wise ([tune_control()]): with `tune = TRUE`
#' the ladder from 10 / 0.4 to 1 / 0.0125 at `min_leaf_rows = 1`, then
#' `min_leaf_rows` 1, 5, 20 and 100 at the chosen point, then the point's two
#' ladder neighbours, about 10 forests in all. `tune = "quick"` grows two
#' forests, 3 / 0.1 and 1 / 0.025 at `min_leaf_rows = 1`. The chosen forest's
#' pruned model is returned, with the display at `max_panels`, and the search in
#' its element `tune`. Every forest is grown from the same random number state,
#' so the result is the model that `rpf(splits = NULL, prune = TRUE)` returns
#' at the chosen settings from that state.
#'
#' On 48 data sets (fresh seeds; 44 synthetic, 4 of California housing) the full
#' search had 10.8% lower test error than `rpf()`'s 10 / 0.4 (no family of data
#' sets worse than +0.9%), and 2.8% lower than 5-fold cross-validation over the
#' same settings at 1/24 of its cost; the quick preset had 9.4% lower test
#' error (worst +0.4%). In those studies a forest was scored with the pruning's
#' choice of the time (the lowest out-of-bag error, R's coefficients fitted in
#' sample); the sparse default of [prune_components()] used here was not
#' compared there. A search costs about as many fits as it grows forests, and
#' their prunings.
#'
#' @references
#' Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
#' directly interpretable tree ensemble. \doi{10.48550/arXiv.2012.14563}
#'
#' Wright, M. N. & König, I. R. (2019). Splitting on categorical predictors in
#' random forests. PeerJ 7:e6339. \doi{10.7717/peerj.6339}
#'
#' Coppersmith, D., Hong, S. J. & Hosking, J. R. (1999). Partitioning nominal
#' attributes in decision trees. Data Mining and Knowledge Discovery 3:197--217.
#' \doi{10.1023/A:1009869804967}
#'
#' @examples
#' # Regression with formula
#' rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
#'
#' # Regression with x and y, allowing pairwise interactions
#' rpfit <- rpf(x = mtcars[, c("cyl", "wt", "hp")], y = mtcars$mpg, max_interaction = 2)
#'
#' # Classification
#' rpfit <- rpf(Species ~ ., data = iris, loss = "exponential")
#'
#' # A factor that keeps its level order: splits separate neighbouring levels
#' cars <- mtcars
#' cars$gear <- factor(cars$gear, levels = c(3, 4, 5))
#' rpfit <- rpf(mpg ~ gear + wt, data = cars, keep_level_order = "gear")
#' rpfit$factor_levels
#'
#' # The out-of-bag error along the splits, from one fit
#' set.seed(1)
#' rpfit <- rpf(mpg ~ ., data = mtcars, ntrees = 20, splits = 100, oob_path = TRUE)
#' rpfit$oob_path[which.min(rpfit$oob_path$oob_rmse), ]
#'
#' # The number of splits chosen out of bag while growing
#' set.seed(1)
#' rpfit <- rpf(mpg ~ ., data = mtcars, ntrees = 20, splits = NULL)
#' rpfit
#' rpfit$params$splits
#'
#' # The whole pipeline: splits chosen out of bag, pruned
#' set.seed(1)
#' pruned <- rpf(mpg ~ ., data = mtcars, ntrees = 20, splits = NULL, prune = TRUE)
#' pruned
#'
#' \donttest{
#' # Tuning out of bag: the quick preset, two forests
#' set.seed(1)
#' x <- matrix(runif(400 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(400, sd = 0.3)
#' tuned <- rpf(x, y, ntrees = 10, tune = "quick")
#' tuned$tune$search
#' }
rpf <- function(x, ...) {
  UseMethod("rpf")
}

#' @export
#' @noRd
rpf.default <- function(x, ...) {
  cli::cli_abort("{.fn rpf} is not defined for a {.cls {class(x)[1]}}.")
}

# Formula method
#' @export
#' @rdname rpf
rpf.formula <- function(
  formula,
  data,
  max_interaction = 3,
  ntrees = 200,
  splits = NULL,
  split_structure = "leaves",
  split_try = NULL,
  t_try = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  prune = NULL,
  max_panels = 10,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
) {
  rlang::check_dots_empty()
  check_no_interaction_terms(formula, data)
  # model.matrix() would expand logical predictors into two dummy columns; use them
  # as 0/1 like the x/y interface. The outcome stays as is so it is still validated.
  outcome <- all.vars(rlang::f_lhs(formula))
  logical_predictors <- setdiff(names(data)[vapply(data, is.logical, logical(1))], outcome)
  for (col in logical_predictors) {
    data[[col]] <- as.integer(data[[col]])
  }
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  # column indices refer to the predictors in the formula's order (the model matrix puts factors last)
  interaction_constraints <- constraint_indices_to_names(
    interaction_constraints,
    labels = attr(stats::terms(formula, data = data), "term.labels"),
    columns = colnames(processed$predictors)
  )
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment()),
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order,
    oob_path = oob_path
  )
}

# XY method - data frame
#' @export
#' @rdname rpf
rpf.data.frame <- function(
  x,
  y,
  max_interaction = 3,
  ntrees = 200,
  splits = NULL,
  split_structure = "leaves",
  split_try = NULL,
  t_try = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  prune = NULL,
  max_panels = 10,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
) {
  rlang::check_dots_empty()
  check_outcome_type(y)
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment()),
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order,
    oob_path = oob_path
  )
}

# XY method - matrix
#' @export
#' @rdname rpf
rpf.matrix <- function(
  x,
  y,
  max_interaction = 3,
  ntrees = 200,
  splits = NULL,
  split_structure = "leaves",
  split_try = NULL,
  t_try = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  prune = NULL,
  max_panels = 10,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
) {
  rlang::check_dots_empty()
  check_outcome_type(y)
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment()),
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order,
    oob_path = oob_path
  )
}

# Recipe method
#' @export
#' @rdname rpf
rpf.recipe <- function(
  x,
  data,
  max_interaction = 3,
  ntrees = 200,
  splits = NULL,
  split_structure = "leaves",
  split_try = NULL,
  t_try = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  delta = 0.001,
  epsilon = 0.1,
  prune = NULL,
  max_panels = 10,
  tune = FALSE,
  purify = FALSE,
  nthreads = 1,
  export_forest = FALSE,
  deterministic = FALSE,
  keep_level_order = FALSE,
  oob_path = FALSE,
  ...
) {
  rlang::check_dots_empty()
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  # the recipe's steps decide the predictor columns: indices have no user-given order to refer to
  if (is.list(interaction_constraints) && any(vapply(interaction_constraints, is.numeric, logical(1)))) {
    cli::cli_abort("With a recipe, give {.arg interaction_constraints} by predictor names, not column indices.")
  }
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    ntrees = ntrees,
    splits = splits,
    split_structure = split_structure,
    split_try = split_try,
    t_try = t_try,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    delta = delta,
    epsilon = epsilon,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment()),
    purify = purify,
    nthreads = nthreads,
    export_forest = export_forest,
    deterministic = deterministic,
    keep_level_order = keep_level_order,
    oob_path = oob_path
  )
}

# Bridge: validates arguments and calls rpf_impl() with processed input
#' @noRd
#' @param processed `[list]`: Output of `hardhat::mold` from respective rpf methods
#' @importFrom hardhat validate_outcomes_are_univariate
rpf_bridge <- function(
  processed,
  max_interaction,
  ntrees,
  splits,
  split_structure,
  split_try,
  t_try,
  max_candidates,
  split_decay_rate,
  delete_leaves,
  min_leaf_rows,
  interaction_constraints,
  loss,
  delta,
  epsilon,
  prune,
  max_panels,
  tune,
  given,
  purify,
  nthreads,
  export_forest,
  deterministic,
  keep_level_order,
  oob_path
) {
  # report errors from helpers against the user's rpf() call
  call <- rlang::caller_env()
  hardhat::validate_outcomes_are_univariate(processed$outcomes)
  check_no_missing(processed$predictors, "The predictors", call = call)
  check_no_missing(processed$outcomes, "The outcome", call = call)
  check_predictor_names(names(processed$predictors), call = call)
  processed <- check_outcome(processed, call = call)
  predictors <- preprocess_predictors_fit(processed, keep_level_order = keep_level_order, call = call)
  checkmate::assert_string(loss, null.ok = TRUE)
  outcomes <- preprocess_outcome(processed, loss, call = call)
  loss <- outcomes$loss
  p <- ncol(predictors$predictors_matrix)

  # Check arguments
  checkmate::assert_int(max_interaction, lower = 0)

  # rewrite max_interaction so 0 -> "maximum", e.g. ncol(X)
  if (max_interaction == 0) {
    max_interaction <- p
  }
  # same applies to values > p, silently, as the default of 2 exceeds p for a single predictor
  if (max_interaction > p) {
    max_interaction <- p
  }

  checkmate::assert_int(ntrees, lower = 1)
  checkmate::assert_int(splits, lower = 1, null.ok = TRUE)
  checkmate::assert_choice(split_structure, choices = c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees"))
  checkmate::assert_int(split_try, lower = 1, null.ok = TRUE)
  checkmate::assert_number(t_try, lower = 0, upper = 1, null.ok = TRUE)
  # NULL greediness: 3 / 0.1 for "leaves", 10 / 0.4 (master's) for the other structures
  split_try <- split_try %||% if (split_structure == "leaves") 3 else 10
  t_try <- t_try %||% if (split_structure == "leaves") 0.1 else 0.4
  if (t_try <= 0) {
    cli::cli_abort("{.arg t_try} must be greater than 0, not {t_try}.", call = call)
  }
  checkmate::assert_int(max_candidates, lower = 1)
  checkmate::assert_number(split_decay_rate, lower = 0)
  checkmate::assert_flag(delete_leaves)
  checkmate::assert_int(min_leaf_rows, lower = 1)
  if (min_leaf_rows > 1 && split_structure != "leaves") {
    cli::cli_abort(
      "{.arg min_leaf_rows} > 1 needs {.code split_structure = \"leaves\"}.",
      call = call
    )
  }

  # "median" loss is implemented but discarded
  loss_functions <- switch(outcomes$mode, "regression" = "L2", "classification" = c("L1", "L2", "logit", "exponential"))
  checkmate::assert_choice(loss, choices = loss_functions)
  checkmate::assert_number(delta, lower = 0, upper = 1)
  checkmate::assert_number(epsilon, lower = 0, upper = 1)
  if (loss %in% c("logit", "exponential")) {
    # Both truncate class proportions to [value, 1 - value], which is empty from 0.5 on
    if (delta >= 0.5 || (loss == "logit" && delta <= 0)) {
      range <- if (loss == "logit") "greater than 0" else "at least 0"
      cli::cli_abort(
        "With {.code loss = \"{loss}\"}, {.arg delta} must be {range} and less than 0.5, not {delta}.",
        call = call
      )
    }
    if (epsilon <= 0 || epsilon >= 0.5) {
      cli::cli_abort(
        "With {.code loss = \"{loss}\"}, {.arg epsilon} must be greater than 0 and less than 0.5, not {epsilon}.",
        call = call
      )
    }
  }

  checkmate::assert_flag(prune, null.ok = TRUE)
  checkmate::assert_flag(deterministic)
  # prune = NULL: the pruned model for a random forest, the forest with deterministic = TRUE (no out-of-bag rows)
  prune <- prune %||% !deterministic
  # tune: the search chooses split_try, t_try, min_leaf_rows and splits, and prunes every forest
  control <- tune_resolve(tune, call = call)
  if (!is.null(control)) {
    tune_check_call(control, given, splits, prune, deterministic, ntrees, split_structure, call = call)
    splits <- NULL
    prune <- TRUE
  }
  checkmate::assert_number(max_panels, lower = 1)
  if (is.finite(max_panels)) {
    checkmate::assert_int(max_panels, lower = 1, upper = .Machine$integer.max)
  }
  checkmate::assert_flag(purify)
  checkmate::assert_int(nthreads, lower = 1L)
  checkmate::assert_flag(export_forest)
  if (prune && deterministic) {
    cli::cli_abort(
      c(
        "{.code prune = TRUE} is not available with {.code deterministic = TRUE}.",
        "i" = "Fit with {.code prune = FALSE} and call {.fn prune_components} on the forest."
      ),
      call = call
    )
  }
  if (deterministic && ntrees > 1) {
    cli::cli_warn(c(
      "With {.code deterministic = TRUE}, all {ntrees} tree families are identical.",
      "i" = "Use {.code ntrees = 1} to fit the same model faster."
    ))
  }
  # splits = NULL: chosen out of bag while growing (a random forest), otherwise 30
  adaptive <- is.null(splits) && !deterministic
  if (is.null(splits) && !adaptive) {
    splits <- 30L
  }
  checkpoints <- NULL
  if (!adaptive) {
    checkpoints <- check_oob_path(oob_path, splits, deterministic, call = call)
  } else if (!isTRUE(oob_path) && !isFALSE(oob_path)) {
    cli::cli_abort(
      c(
        "With {.code splits = NULL}, {.arg oob_path} must be {.code TRUE} or {.code FALSE}.",
        "i" = "The out-of-bag path of the checkpoints grown is stored either way."
      ),
      call = call
    )
  }

  interaction_constraints <- check_interaction_constraints(
    interaction_constraints,
    colnames(predictors$predictors_matrix),
    max_interaction,
    call = call
  )

  params <- list(
    loss = loss,
    ntrees = ntrees,
    max_interaction = max_interaction,
    splits = if (adaptive) NA_integer_ else splits,
    split_try = split_try,
    t_try = t_try,
    split_decay_rate = split_decay_rate,
    max_candidates = max_candidates,
    delete_leaves = delete_leaves,
    min_leaf_rows = as.integer(min_leaf_rows),
    interaction_constraints = interaction_constraints,
    split_structure = split_structure,
    delta = delta,
    epsilon = epsilon,
    deterministic = deterministic,
    nthreads = nthreads,
    # the pruning leaves the forest unpurified (prune_components()), so that rpf(prune = TRUE) is its pruning
    purify = purify && !prune,
    exact_draw = rpf_exact_draw(call = call),
    keep_level_order = predictors$keep_level_order
  )

  if (!is.null(control)) {
    return(tune_run(processed, outcomes, predictors, params, control, max_panels, export_forest, call = call))
  }
  res <- rpf_fit(processed, outcomes, predictors, params, adaptive, checkpoints, export_forest, call = call)
  if (prune) {
    res <- rpf_prune(res, processed, max_panels, call = call)
  }
  res
}

# The forest of rpf_bridge()'s checked arguments: `adaptive` whether `splits` is chosen out of bag while growing,
# `checkpoints` those of rpf(oob_path = ) (fixed splits) or NULL. The OOB error is the regression one or the
# classification criterion (oob_score_y()).
rpf_fit <- function(processed, outcomes, predictors, params, adaptive, checkpoints, export_forest, call) {
  score_y <- if (adaptive || !is.null(checkpoints)) {
    oob_score_y(outcomes$mode, params$loss, processed$outcomes[[1]])
  }
  if (adaptive) {
    chosen <- tune_splits_oob(
      outcomes$outcomes,
      predictors$predictors_matrix,
      params,
      mode = outcomes$mode,
      score_y = score_y,
      call = call
    )
    fit <- chosen$fit
    params$splits <- chosen$splits
  } else {
    fit <- rpf_impl(
      Y = outcomes$outcomes,
      X = predictors$predictors_matrix,
      mode = outcomes$mode,
      params = params,
      oob_path = checkpoints
    )
  }

  # Optionally export a compact R list representation of the forest.
  forest <- NULL
  if (export_forest) {
    forest <- fit$get_model()
    class(forest) <- "rpf_forest"
  }

  res <- new_rpf(
    fit = fit,
    blueprint = processed$blueprint,
    mode = outcomes$mode,
    factor_levels = predictors$factor_levels,
    params = params,
    forest = forest
  )
  if (adaptive) {
    res$oob_path <- chosen$path
    res$splits_oob <- chosen[c("splits", "stopped_at", "capped")]
  } else if (!is.null(checkpoints)) {
    res$oob_path <- oob_path_table(fit$oob_path(), score_y, params$ntrees, call = call)
  }
  res
}

# rpf(prune = TRUE): prune_components() with its defaults on the training rows, the same as prune_components(res, data)
# after rpf(prune = FALSE), random numbers included. max_panels: the display's (Inf: none); store: prune_run()'s
rpf_prune <- function(res, processed, max_panels, store = TRUE, call) {
  defaults <- lapply(formals(prune_components)[c("rule", "tau", "select", "select_by", "folds")], eval)
  prune_run(
    res,
    prune_rows(res, processed$predictors, processed$outcomes[[1]]),
    max_panels = max_panels,
    rule = defaults$rule,
    tau = defaults$tau,
    select = defaults$select[[1]],
    select_by = defaults$select_by[[1]],
    default_select_by = TRUE,
    folds = defaults$folds,
    fold_id = NULL,
    seed = NULL,
    nthreads = as.integer(res$params$nthreads),
    verbose = FALSE,
    store = store,
    call = call
  )
}

# Intermediate to hold model object with blueprint used for prediction
new_rpf <- function(fit, blueprint, ...) {
  hardhat::new_model(
    fit = fit,
    blueprint = blueprint,
    class = "rpf",
    ...
  )
}

# The draw of candidates for split_structure = "leaves": the option
# randomPlantedForest.exact_draw (FALSE, the default: the fast draw; TRUE: the
# exact draw, the same distribution with other random numbers).
rpf_exact_draw <- function(call = rlang::caller_env()) {
  exact_draw <- getOption("randomPlantedForest.exact_draw", FALSE)
  if (!rlang::is_bool(exact_draw)) {
    cli::cli_abort(
      "The option {.code randomPlantedForest.exact_draw} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }
  exact_draw
}

# Assemble the positional parameter vector for the C++ constructors.
# 15 elements for regression (+ min_leaf_rows, exact_draw), 17 (+ delta,
# epsilon, min_leaf_rows, exact_draw) for classification. Models saved before
# these options existed have no such entries in their params and were fitted
# with min_leaf_rows = 1 and the exact draw.
rpf_param_vector <- function(params, mode) {
  split_mode <- switch(
    params$split_structure,
    res_trees = 0L,
    cur_trees_2 = 1L,
    cur_trees_1 = 2L,
    leaves = 3L,
    hist = 4L
  )
  base <- c(
    params$max_interaction,
    params$ntrees,
    params$splits,
    params$split_try,
    params$t_try,
    params$purify,
    params$deterministic,
    params$nthreads,
    0, # cross-validation slot, unused by the C++ core
    params$split_decay_rate,
    params$max_candidates,
    params$delete_leaves,
    split_mode
  )
  if (mode == "classification") {
    base <- c(base, params$delta, params$epsilon)
  }
  c(base, params$min_leaf_rows %||% 1L, params$exact_draw %||% TRUE)
}

# Main fitting function and interface to C++ implementation
# `oob_path`: NULL, or the checkpoints of the OOB path (see check_oob_path())
# `adaptive`: NULL, or the adaptive fit's checkpoints, stopper, patience, floor and pause_gb (see tune_splits_oob())
rpf_impl <- function(Y, X, mode = c("regression", "classification"), params, oob_path = NULL, adaptive = NULL) {
  # Final input validation, should be superfluous
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
  mode <- match.arg(mode)
  pars <- rpf_param_vector(params, mode)
  # allowed feature sets as 1-based column indices (NULL: no constraint, the constructors without them)
  constraints <- if (!is.null(params$interaction_constraints)) {
    lapply(params$interaction_constraints, match, table = colnames(X))
  }
  # configure only (no data, no fit): the constructors that fit give the same forest
  configure <- function() {
    if (is.null(constraints)) {
      if (mode == "classification") new(ClassificationRPF, params$loss, pars) else new(RandomPlantedForest, pars)
    } else {
      if (mode == "classification") {
        new(ClassificationRPF, params$loss, pars, constraints)
      } else {
        new(RandomPlantedForest, pars, constraints)
      }
    }
  }

  if (!is.null(adaptive)) {
    # configure, then fit one family at a time: the forest at the chosen number of splits
    fit <- configure()
    fit$set_adaptive(
      as.integer(adaptive$checkpoints),
      adaptive$stopper,
      adaptive$patience,
      adaptive$floor,
      adaptive$pause_gb
    )
    fit$set_data(Y, X)
    fit$fit()
  } else if (!is.null(oob_path)) {
    # configure, then fit: the same forest as the constructor that fits
    fit <- configure()
    fit$set_oob_path(as.integer(oob_path))
    fit$set_data(Y, X)
    fit$fit()
  } else if (is.null(constraints)) {
    if (mode == "classification") {
      fit <- new(ClassificationRPF, Y, X, params$loss, pars)
    } else {
      fit <- new(RandomPlantedForest, Y, X, pars)
    }
  } else {
    if (mode == "classification") {
      fit <- new(ClassificationRPF, Y, X, params$loss, pars, constraints)
    } else {
      fit <- new(RandomPlantedForest, Y, X, pars, constraints)
    }
  }

  fit
}

# Formula interface: translate the sets given by column indices into predictor names, index i referring to the
# i-th predictor of the formula (`labels`, its term labels; `y ~ .` gives the data's column order). The model matrix
# (`columns`) orders the predictors differently: hardhat puts factors last. Sets of names are left as they are;
# check_interaction_constraints() validates everything else.
constraint_indices_to_names <- function(interaction_constraints, labels, columns, call = rlang::caller_env()) {
  if (!is.list(interaction_constraints)) {
    return(interaction_constraints)
  }
  lapply(seq_along(interaction_constraints), function(i) {
    s <- interaction_constraints[[i]]
    if (!is.numeric(s) || length(s) == 0 || anyNA(s)) {
      return(s)
    }
    if (any(s != round(s)) || any(s < 1) || any(s > length(labels))) {
      cli::cli_abort(
        paste(
          "{.code interaction_constraints[[{i}]]} must contain indices between 1 and {length(labels)}",
          "of the formula's predictors (or predictor names)."
        ),
        call = call
      )
    }
    nm <- labels[s]
    if (!all(nm %in% columns)) {
      bad <- nm[!nm %in% columns]
      cli::cli_abort(
        c(
          "{.code interaction_constraints[[{i}]]}: the formula's predictor{?s} {.val {bad}} {?is/are} not a single model column.",
          "i" = "Give the sets by column names: {.val {columns}}."
        ),
        call = call
      )
    }
    nm
  })
}

# Validate `interaction_constraints` and convert it to a list of sorted, unique predictor
# names (NULL stays NULL: no constraint).
check_interaction_constraints <- function(interaction_constraints, feature_names, max_interaction, call) {
  if (is.null(interaction_constraints)) {
    return(NULL)
  }
  p <- length(feature_names)
  if (!is.list(interaction_constraints) || length(interaction_constraints) == 0) {
    cli::cli_abort(
      c(
        "{.arg interaction_constraints} must be {.code NULL} or a non-empty list of feature sets.",
        "i" = "For example {.code list(c(\"x1\", \"x2\"), \"x3\")}."
      ),
      call = call
    )
  }
  sets <- lapply(seq_along(interaction_constraints), function(i) {
    s <- interaction_constraints[[i]]
    if (length(s) == 0 || anyNA(s)) {
      cli::cli_abort(
        "{.code interaction_constraints[[{i}]]} must be a non-empty vector without {.code NA}.",
        call = call
      )
    }
    if (is.character(s)) {
      unknown <- setdiff(s, feature_names)
      if (length(unknown) > 0) {
        cli::cli_abort(
          c(
            "{.code interaction_constraints[[{i}]]} names unknown predictor{?s} {.val {unknown}}.",
            "i" = "The predictors are {.val {feature_names}}."
          ),
          call = call
        )
      }
      idx <- match(s, feature_names)
    } else if (is.numeric(s)) {
      if (any(s != round(s)) || any(s < 1) || any(s > p)) {
        cli::cli_abort(
          "{.code interaction_constraints[[{i}]]} must contain column indices between 1 and {p} (or predictor names).",
          call = call
        )
      }
      idx <- as.integer(s)
    } else {
      cli::cli_abort(
        "{.code interaction_constraints[[{i}]]} must be predictor names or column indices, not {.obj_type_friendly {s}}.",
        call = call
      )
    }
    feature_names[sort(unique(idx))]
  })
  if (max(lengths(sets)) > max_interaction) {
    cli::cli_inform(c(
      "i" = paste(
        "{.arg interaction_constraints} has a set of {max(lengths(sets))} predictors,",
        "but {.arg max_interaction} is {max_interaction}."
      ),
      " " = "Trees use at most {max_interaction} predictor{?s}."
    ))
  }
  sets
}
