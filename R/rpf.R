#' Random Planted Forest
#'
#' Fits a random planted forest (Hiabu, Mammen and Meyer, 2020) for regression
#' or classification. The model is a sum of trees that each split on at most
#' `max_interaction` predictors, so it decomposes into main effects and
#' interaction components; see [predict_components()].
#'
#' By default, `rpf()` returns the pruned model that can be shown, in one call:
#' 1. it grows 200 tree families of up to three-way interactions, with the
#'    number of splits chosen out of bag while they grow (`n_splits = NULL`);
#' 2. it prunes the forest's components by [prune_components()] with its
#'    defaults (`prune = NULL`): the sparse rule, chosen out of bag;
#' 3. it coarsens every kept component of three or more predictors to at most
#'    `max_panels = 10` heat maps and prunes again.
#'
#' The result, of class `rpf_pruned`, has methods for
#' [`predict()`][predict.rpf_pruned], [predict_components()] and `print()`.
#' `prune = FALSE` returns the forest (class `rpf`) instead, and
#' `n_splits = <number>` fits a fixed number of splits.
#'
#' @param x,y `[data.frame | matrix | recipe]`, `[numeric | factor]`:
#'   Predictors and outcome for the x/y interface.
#'   `x` holds the predictors, or is a [`recipe`][recipes::recipe] (then supply
#'   `data` instead of `y`). The predictors need unique names that do not
#'   contain `:`, which joins predictor names in component names such as `x1:x2`.
#'   `y` is the outcome vector: __numeric__ for regression, __factor__ for
#'   classification. A numeric outcome must be finite. A factor outcome needs at
#'   least two classes; levels that do not occur are dropped with a warning.
#' @param ... Only the old names of renamed arguments (`ntrees`, `splits`,
#'   `split_try`, `t_try`), deprecated: they work, with a warning. Any other
#'   (e.g. misspelled) argument is an error. See "Renamed arguments" in Details.
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
#' @param n_families `[integer(1): 200]`: Number of tree families, i.e. the size of the forest.
#'   Each family is grown on a bootstrap sample and predictions are averaged. Was `ntrees`.
#' @param n_splits `[integer(1) | NULL: NULL]`: Number of split rounds per tree family. The main tuning
#'   parameter: more splits fit more complex functions. A round splits the best of its
#'   candidates; near saturation, when every candidate's leaf is too small to split, a
#'   round can make no split. `NULL` (the default) chooses it out of bag while the forest
#'   grows (with `deterministic = TRUE`, which has no out-of-bag rows, `NULL` means 30). See
#'   "Choosing `n_splits` out of bag" in Details. Was `splits`.
#' @param split_structure `[character(1): "leaf"]`: How split candidates are formed and
#'   sampled; one of `"leaf"`, `"leaf_binned"`, `"tree_then_leaf"`,
#'   `"tree_then_all_leaves"` or `"target_tree_then_leaf"`. See Details. The
#'   names of 0.5.0 (`"leaves"`, `"hist"`, `"cur_trees_1"`, `"cur_trees_2"`,
#'   `"res_trees"`) still work, with a deprecation warning.
#' @param n_try_split_points `[integer(1) | NULL: NULL]`: Number of split points (thresholds)
#'   tried per split candidate. `NULL` means 3 for `split_structure = "leaf"` and 10 for
#'   the other structures. Was `split_try`.
#' @param share_candidates `[numeric(1) | NULL: NULL]`: Share in `(0, 1]` of all possible
#'   splits drawn as candidates in each round, capped at `max_candidates`. `NULL`
#'   means 0.1 for `split_structure = "leaf"` and 0.4 for the other structures.
#'   See "Choosing parameters" in Details. Was `t_try`.
#' @param max_candidates `[integer(1): 50]`: Maximum number of split candidates per round.
#' @param split_decay_rate `[numeric(1): 0.1]`: Down-weights possible splits that were
#'   sampled as candidates but not chosen. Each such round ages a split by one,
#'   choosing it resets its age to zero, and splits are sampled with weight
#'   `exp(-split_decay_rate * age)`. `0` samples uniformly.
#' @param delete_leaves `[logical(1): TRUE]`: Whether to delete a parent leaf when
#'   splitting along an existing dimension.
#' @param min_leaf_rows `[integer(1): 1]`: Minimum number of rows of the
#'   bootstrap sample on each side of every split, so in every leaf, for every
#'   `split_structure` (regression and classification). Split points are tried
#'   only among those that keep this many rows on each side, and leaves without
#'   one are not drawn (see "`split_structure`" in Details). `1` means no limit.
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
#'   regression, the only supported loss, and `"logit"` (the multinomial
#'   log-likelihood) for classification, with two classes and with more.
#'   Classification also supports `"exponential"` (with more than two classes
#'   the multiclass exponential loss SAMME), and `"L1"` and `"L2"`, which fit
#'   class indicators directly and do not yield proper probability estimates.
#'   See "Classification losses" in Details.
#' @param exact_step `[logical(1): TRUE]`: Only used with `loss = "logit"` (the
#'   classification default): whether the two leaves of each split get the
#'   exact minimisers of the penalised loss (Newton's method to convergence on
#'   the leaf's rows), or one Newton step from the sums the split search scored
#'   them by (`FALSE`, faster). The split
#'   search is the same either way. The exponential loss's step is exact by
#'   itself, so `exact_step` does not change `"exponential"` fits. See
#'   "Classification losses" in Details.
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
#'   to tune the greediness (`n_try_split_points`, `share_candidates`) and
#'   `min_leaf_rows` out of bag (not with `deterministic = TRUE`, at least 10
#'   tree families). `TRUE` runs the full search ([tune_control()]), `"quick"` the
#'   two-forest preset, and a [tune_control()] object your own. Every forest of
#'   the search has `n_splits` chosen out of bag and is pruned, so `tune` implies
#'   `n_splits = NULL` and `prune = TRUE`, and giving `n_try_split_points`,
#'   `share_candidates` or `min_leaf_rows` with it is an error. See "Tuning out
#'   of bag" in Details.
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
#'   (`share_candidates` and `split_decay_rate` are ignored) and evenly spaced thresholds. The
#'   result does not depend on the seed, and all tree families are identical, so use
#'   `n_families = 1`. There are no out-of-bag rows, so `n_splits = NULL` means 30
#'   and `prune = NULL` means `FALSE`. Mainly useful for testing and debugging.
#' @param keep_level_order `[logical(1) | character: FALSE]`: Which factor
#'   predictors keep the order of their levels instead of being ordered by the
#'   outcome. `FALSE` orders every factor by the outcome, `TRUE` keeps the level
#'   order of every factor, and a character vector keeps it for the named
#'   factor predictors only. See "Factor predictors" below.
#' @param oob_path `[logical(1) | integer: FALSE]`: Whether to track the
#'   out-of-bag (OOB) error along the splits (not with
#'   `deterministic = TRUE`). `TRUE` evaluates it at the checkpoints of
#'   [oob_path_checkpoints()] for `n_splits`, four per doubling; an increasing
#'   integer vector in `1..n_splits` gives the checkpoints. The fitted forest is
#'   the same as without. See Details. With `n_splits = NULL` the path of the
#'   checkpoints grown is stored either way, and only `TRUE` or `FALSE` is
#'   allowed.
#'
#' @return With the defaults, an object of class `rpf_pruned` (see below). With
#'   `prune = FALSE`, an object of class `"rpf"` with model object contained in `$fit`. With
#'   `oob_path`, `$oob_path` holds the OOB path: a data frame with one row per
#'   checkpoint, `n_splits`, the corrected and the plain OOB root mean squared
#'   error (`oob_rmse`, `oob_rmse_raw`), and the standard errors of the two mean
#'   squared errors (`se_mse`, `se_mse_raw`). Its attribute `n_rows` is the
#'   number of training rows with an OOB tree family. For classification with
#'   `loss = "logit"` or `"exponential"`, or more than two classes, the columns
#'   are the mean corrected and plain criterion (`oob_loss`, `oob_loss_raw`) and
#'   their standard errors (`se_loss`, `se_loss_raw`), and the attribute
#'   `criterion` names it (see "Classification" in Details). With `n_splits = NULL`,
#'   `$params$splits` is the chosen number of splits (`$params` keeps the names
#'   of 0.5.0, see "Renamed arguments" in Details), `$oob_path` the path of
#'   the checkpoints grown, and `$splits_oob` a list with the chosen `n_splits`,
#'   `stopped_at` (the number of splits where growing stopped) and `capped`
#'   (whether it stopped at the cap while the error was still falling, neither
#'   stop rule reached).
#'   When pruning (the default), an object of class `rpf_pruned` (see
#'   [prune_components()]), whose element `fit` is the forest (not purified),
#'   with the elements above: the pruned model's OOB path is `$fit$oob_path`, and
#'   with `n_splits = NULL` the chosen number of splits `$fit$params$splits` and
#'   the stop `$fit$splits_oob`.
#'   With `tune`, the `rpf_pruned` model of the chosen forest, with an element
#'   `tune`: `search`, the search table (one row per forest grown, in the order
#'   grown: `step`, `n_try_split_points`, `share_candidates`, `min_leaf_rows`,
#'   the chosen `n_splits`, `stopped_at`, `capped`, the forest's corrected OOB
#'   RMSE `oob_rmse_forest` (as the pruning scores it: its rule `"D"` at
#'   `tau = 0`, which keeps every nonzero component; the OOB path's `oob_rmse` at
#'   the chosen `n_splits` leaves the rows without OOB families of both halves
#'   uncorrected and is higher),
#'   its pruned model's `rule`, `tau`, `n_kept` and corrected OOB RMSE
#'   `oob_rmse` (the score; for classification both are the square root of the
#'   mean corrected criterion), `seconds`, whether each step chose it
#'   (`chosen_step1`, `chosen_step2`, `chosen_step3`; `NA` if not a candidate of
#'   the step) and `final`); `steps`, each step's candidates and choice;
#'   `settings`, the chosen `n_try_split_points`, `share_candidates`,
#'   `min_leaf_rows` and `n_splits`; `control`, the [tune_control()] used; and `seconds`.
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
#' The defaults were chosen on 44 simulated data sets and California housing
#' (regression) and carried over to classification.
#'
#' `max_interaction` sets which effects the model can represent: `1` for an
#' additive model, `2` to add pairwise interactions, `3` (the default)
#' three-way interactions. Order 3 had 15% lower test error than order 2 on
#' the simulated sets and 1.5% lower on California housing. Higher orders gain
#' only where the truth has such interactions, cost much more time (no limit
#' was 11 to 41 times slower on California housing, and worse there), and the
#' pruned model rarely keeps a component of more than three predictors, so it
#' is shown in at most `max_panels` heat maps anyway. With a higher
#' `max_interaction`, splits are spread over more possible components.
#'
#' `n_splits` sets how closely the model follows the data. By default it is
#' chosen out of bag while the forest grows (see "Choosing `n_splits` out of
#' bag"), so it need not be tuned.
#'
#' `n_families` averages over bootstrap samples to reduce variance, and fitting and
#' prediction time grow linearly with it. Each tree family is itself a complete
#' model, a sum of trees covering all components, rather than a single tree,
#' but averaging more of them still pays off: the test error fell by about 2%
#' per doubling of the families up to 200 on the simulated sets, and the
#' pruned models were sparser (California housing: 34 instead of 47 kept
#' components at 200 instead of 50 families). Hence the default of 200, at the
#' price of time; use fewer for a quick look.
#'
#' The split search parameters (`split_structure`, `n_try_split_points`,
#' `share_candidates`, `max_candidates`, `split_decay_rate`) trade speed for a
#' more thorough, more greedy search. For `"leaf"`, a less greedy search is
#' better: 3 split points per candidate and 10% of the possible splits per
#' round (the default) had 5.7% to 6.8% lower test error than 10 and 40% (the
#' default of the other structures, and of `"leaf"` before). `tune` searches
#' this greediness out of bag.
#'
#' For classification, `"logit"` (the default) and `"exponential"` fit on the
#' link scale and are transformed to probabilities; see "Classification
#' losses". `"L1"` and `"L2"` fit class indicators directly and their
#' probabilities are truncated to `[0, 1]`.
#'
#' ## Classification losses
#'
#' `"exponential"` and `"logit"` fit scores on the link scale: with two classes
#' one score, the log-odds of the second class; with more, `"logit"` one score
#' per class but the first (whose score is 0) and `"exponential"` one score per
#' class, summing to 0 (the multiclass exponential loss SAMME of Zhu, Zou,
#' Rosset and Hastie, 2009, whose row weight is `exp(-g)`, `g` the score of its
#' class). The probabilities are the softmax of the scores.
#'
#' Every split moves the scores of the rows of each new leaf by one step, which
#' minimises the leaf's penalised loss: the loss of its rows, plus that of one
#' pseudo-row (split evenly over the classes) at the leaf's current mean score.
#' The pseudo-row pulls the step towards no change, so a leaf whose rows
#' all have one class moves towards a finite level (about `log(2 n + 1)` for
#' `n` rows and two classes) and is not moved further once it is there, and no
#' clamp of the class shares is needed. A split is scored by its gain over the
#' parent leaf's own step, and each new leaf pays for its pseudo-row, so
#' splitting a leaf of one class gains nothing, nor does repeating a split whose
#' leaves are at their optimum, and such splits are not made. A round in which
#' no candidate gains makes no split. The candidates' gains come from sums over
#' their two sides (the second-order Newton gain for `"logit"`, the exact gain
#' for `"exponential"`), built as regression builds its sums (along the leaf's
#' sorted rows, so that a split point costs the same whatever the leaf's size;
#' for `"tree_then_leaf"` one pass over the leaf per split point).
#' `"exponential"`'s step is exact; with `exact_step = TRUE` (the default)
#' `"logit"`'s steps are made exact too, otherwise they are one Newton step.
#'
#' `"logit"` is the default, with two classes and with more. On 8 public data
#' sets (4 binary, 4 with 3 to 7 classes; the default pipeline, 5 repeated
#' splits each) its test log-loss was lower than or tied with that of
#' `"exponential"` on every set, and lower by 0.053 on average with 3 or more
#' classes. `"exponential"` can need many more splits with many classes: on a
#' set of 7 classes it chose 1,345 to 4,525 splits where `"logit"` chose about
#' 22.
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
#' In each round, a `share_candidates` share of all possible splits (capped at
#' `max_candidates`) is drawn as candidates with weights
#' `exp(-split_decay_rate * age)`. `split_structure` defines what a candidate is
#' and how its split points are tried. The names say the candidate, then how
#' the leaf to split is found.
#'
#' - `"leaf"` (was `"leaves"`): Split candidates are (leaf, split-dimension) pairs. For each
#'   sampled candidate, `n_try_split_points` split points are drawn uniformly without replacement
#'   from the distinct values of that predictor among the leaf's rows, leaving out the smallest
#'   and, with `min_leaf_rows`, those that leave fewer rows on either side (all of them if
#'   there are fewer), and evaluated to choose the best split.
#'
#' - `"leaf_binned"` (was `"hist"`): As `"leaf"`, but split points are drawn
#'   from boundaries of quantile bins computed once per predictor, which makes
#'   evaluating them faster on large data. With `min_leaf_rows`, only from the
#'   boundaries that leave that many rows on each side, and candidates without
#'   one are not drawn.
#'
#' - `"tree_then_leaf"` (was `"cur_trees_1"`): Split candidates are
#'   (current-tree, split-dimension) pairs. For each sampled candidate, perform
#'   `n_try_split_points` evaluations. Each evaluation samples a leaf from the
#'   set of valid current trees (with probability proportional to its number of
#'   available split points) and then uniformly samples a single split point
#'   within that leaf. With `min_leaf_rows`, a leaf's available split points
#'   are those that leave that many rows on each side, so a leaf without one is
#'   not drawn.
#'
#' - `"tree_then_all_leaves"` (was `"cur_trees_2"`): Split candidates are
#'   (current-tree, split-dimension) pairs. For each sampled candidate, iterate
#'   through every valid leaf. Within each leaf, sample `n_try_split_points`
#'   split points uniformly and evaluate them. With `min_leaf_rows`, only split
#'   points that leave that many rows on each side, and leaves without one are
#'   skipped.
#'
#' - `"target_tree_then_leaf"` (was `"res_trees"`): Split candidates are
#'   resulting trees. For each sampled candidate, run `n_try_split_points`
#'   evaluations by sampling a (split-dimension, leaf) pair from all valid pairs
#'   (with probability proportional to its number of available split points),
#'   then uniformly sampling one split point within that pair. With
#'   `min_leaf_rows`, the available split points are those that leave that many
#'   rows on each side, as for `"tree_then_leaf"`.
#'
#' ## Renamed arguments
#'
#' Four arguments and the split structures were renamed after version 0.5.0:
#' `ntrees` is now `n_families`, `splits` is `n_splits`, `split_try` is
#' `n_try_split_points` and `t_try` is `share_candidates`; the split structures
#' are listed above. The prefix gives the unit: `n_` a count, `share_` a
#' share in `(0, 1]`, `max_` a cap. The old names still work and give the same
#' fit, with a deprecation warning once per session for each old name (as
#' lifecycle's soft deprecation; `options(lifecycle_verbosity = "warning")`
#' warns every time, `"quiet"` never); giving an argument under both names is
#' an error. The fitted model's `$params` keeps the old names
#' (`ntrees`, `splits`, `split_try`, `t_try`, and the old `split_structure`
#' values), so code that reads them, and models saved by 0.5.0, keep working.
#'
#' ## The draw of candidates for `"leaf"`
#'
#' With `split_structure = "leaf"` (regression and classification), the
#' global option `randomPlantedForest.exact_draw` chooses how the candidates
#' are drawn.
#' With `FALSE` (the default), the candidates are drawn from a Fenwick tree
#' over the weights, kept up to date between rounds, which costs a few random
#' numbers per candidate. With `TRUE`, the exact draw of earlier versions,
#' every possible split gets a random number in every round. Both draw from
#' the same distribution, but with other random numbers, so the same seed gives
#' another (equally good) forest. Ties between equally good candidates may go
#' to different candidates under the two draws, so with many splits, near
#' saturation, their fits differ in distribution, though not in accuracy. With
#' `min_leaf_rows`, both leave out the same candidates: those whose leaf has no
#' threshold with that many rows on each side. The fast draw is about 2.5 times
#' as fast on large data with many splits. The option is read when the forest
#' is fitted and stored in `params$exact_draw`; it does not affect the other
#' split structures. To reproduce a regression fit of version 0.5.0, set it and
#' give 0.5.0's defaults, as the example "A regression forest of version 0.5.0"
#' does:
#'
#' ```r
#' options(randomPlantedForest.exact_draw = TRUE)
#' ```
#'
#' Classification fits of 0.5.0 cannot be reproduced: they changed for every
#' split structure and loss (see NEWS: classification uses regression's split
#' structures, the largest value of a predictor can be split off, and three or
#' more classes with the exponential loss use SAMME).
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
#' with fewer splits, for example to choose `n_splits` without cross-validation.
#' This costs `2 * nrow(x) * length(checkpoints)` doubles and a little time;
#' the fit itself is unchanged.
#'
#' The OOB prediction of a row is the mean over its OOB families, and the plain
#' error (`oob_rmse_raw`) is the root mean squared error of these predictions
#' over the rows with an OOB family. That mean averages fewer families than the
#' forest, so the plain error overstates the forest's error by the families'
#' extra variance, more so with more splits. The corrected error (`oob_rmse`)
#' subtracts this excess, `(1 / m_i - 1 / n_families) * v_i` for row `i` with `m_i`
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
#' ## Choosing `n_splits` out of bag
#'
#' With `n_splits = NULL` (not deterministic), the number of splits
#' is chosen while the forest grows, by the corrected OOB error of the
#' out-of-bag path. The families grow to checkpoints, four per doubling of the
#' splits ([oob_path_checkpoints()] up to 25,600), and the forest's corrected
#' error is followed checkpoint by checkpoint. Growing stops once the lowest
#' error is one doubling of the splits back (for example, no new minimum from
#' 100 to 200 splits), or once it is on a plateau: the lowest error improved
#' over the last doubling by at most
#' `getOption("randomPlantedForest.adaptive_tolerance", 1)` standard errors of
#' that improvement (paired over the training rows: the standard deviation of
#' the rows' differences in corrected error between the two checkpoints,
#' divided by the square root of their number). A tolerance of 0 stops only on
#' the first rule. Growing also stops at 25,600 splits, with a warning if the
#' error was still falling there. The result is the forest at the lowest error
#' (ties go to fewer splits; a corrected error at or below 0, a correction that
#' overshot, is not chosen unless all are): the forest that `rpf()` with that
#' many `n_splits` and the same seed returns. On 52 data sets (44 synthetic, 8 splits of
#' California housing) the first rule alone chose the same `n_splits` as the
#' lowest error along a fixed budget of 3,200 splits on 50 of them, and a
#' later, lower minimum on the other two, while growing 0.72 of the splits.
#' The plateau stop chose the same `n_splits` as the first rule alone on all 44
#' synthetic sets, the 4 California sets tried and 24 classification fits
#' (logit and exponential loss), while growing 10% (regression) to 29%
#' (classification) fewer splits; it matters most on long, flat out-of-bag
#' curves.
#'
#' The families grow one at a time, `nthreads` together, each to the stop
#' estimated from the families grown before it. A grown family keeps only its
#' OOB predictions at the checkpoints and a record of its splits (40 bytes a
#' split), from which the chosen forest is rebuilt, without a refit. A family
#' that reaches the estimate waits there with its growing state while the
#' waiting families fit in `getOption("randomPlantedForest.adaptive_pause_gb",
#' 0.25)` GB, so that it can grow on cheaply if the stop comes later; otherwise
#' it is freed, and grown again from its seed if needed. The option trades
#' memory for time; the fit depends neither on it nor on `nthreads`. A budget
#' of 0.25 GB instead of 1 GB lowered the peak memory of the default call by
#' about 0.9 GB on 80,000 simulated rows in the same time, while a budget of 0
#' can cost much more time (families grown again). Growing
#' all families together would hold every family's growing state until the
#' stop: about 18 GB of resident memory (RSS) on California housing (16,512
#' rows) at 200 families and `max_interaction = 3`, where the fit one family at
#' a time peaks at 3.4 GB of RSS (2.3 GB of physical footprint; the adaptive
#' fit alone, measured when it was added: `n_try_split_points = 3`,
#' `share_candidates = 0.1`, the exact draw, 8 threads).
#'
#' ## Pruning
#'
#' By default (`prune = NULL`, or `TRUE`), `rpf()` fits the forest, then prunes its components
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
#' [prune_components()]. Together with `n_splits = NULL`, one call chooses the
#' number of splits out of bag, prunes and returns the model to show. Pruning
#' needs out-of-bag rows, so with `deterministic = TRUE` `prune = NULL` means
#' `FALSE` and `prune = TRUE` is an error.
#'
#' ## Classification
#'
#' A classification forest goes through the same steps (the OOB path,
#' `n_splits = NULL`, `prune`, `max_panels` and `tune`) with its score columns:
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
#'   form, `(1 / m_i - 1 / n_families) * d' H d / (2 * (1 / m_even + 1 / m_odd))`
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
#' (`n_try_split_points` and `share_candidates`, on the ladder of [greediness_ladder()]) and the
#' minimum leaf size `min_leaf_rows`, without cross-validation. Every forest of
#' the search is grown with `n_splits` chosen out of bag and pruned by
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
#' so the result is the model that `rpf(n_splits = NULL, prune = TRUE)` returns
#' at the chosen settings from that state.
#'
#' On 48 data sets (fresh seeds; 44 synthetic, 4 of California housing) the full
#' search had 10.8% lower test error than 10 / 0.4 (the default of earlier versions; no family of data
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
#' # The default: splits chosen out of bag, pruned, every component shown
#' set.seed(1)
#' model <- rpf(mpg ~ ., data = mtcars)
#' model
#' predict(model, mtcars[1:3, ])
#'
#' # Classification, the same steps
#' set.seed(1)
#' model <- rpf(Species ~ ., data = iris, n_families = 50)
#' predict(model, iris[c(1, 51, 101), ], type = "prob")
#'
#' # The forest instead of the pruned model, with x and y, pairwise interactions
#' # and a fixed number of splits
#' set.seed(1)
#' rpfit <- rpf(
#'   x = mtcars[, c("cyl", "wt", "hp")],
#'   y = mtcars$mpg,
#'   max_interaction = 2,
#'   n_splits = 30,
#'   prune = FALSE
#' )
#' rpfit
#'
#' # A factor that keeps its level order: splits separate neighbouring levels
#' cars <- mtcars
#' cars$gear <- factor(cars$gear, levels = c(3, 4, 5))
#' rpfit <- rpf(mpg ~ gear + wt, data = cars, keep_level_order = "gear")
#' rpfit$factor_levels
#'
#' # The out-of-bag error along the splits, from one fit
#' set.seed(1)
#' rpfit <- rpf(
#'   mpg ~ .,
#'   data = mtcars,
#'   n_families = 20,
#'   n_splits = 100,
#'   oob_path = TRUE,
#'   prune = FALSE
#' )
#' rpfit$oob_path[which.min(rpfit$oob_path$oob_rmse), ]
#'
#' \donttest{
#' # A regression forest of version 0.5.0: its defaults and the exact draw
#' old <- options(randomPlantedForest.exact_draw = TRUE)
#' set.seed(1)
#' rpfit <- rpf(
#'   mpg ~ cyl + wt,
#'   data = mtcars,
#'   max_interaction = 2,
#'   n_families = 50,
#'   n_splits = 30,
#'   n_try_split_points = 10,
#'   share_candidates = 0.4,
#'   prune = FALSE
#' )
#' options(old)
#' }
#'
#' \donttest{
#' # Tuning out of bag: the quick preset, two forests
#' set.seed(1)
#' x <- matrix(runif(400 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
#' y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(400, sd = 0.3)
#' tuned <- rpf(x, y, n_families = 10, tune = "quick")
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
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
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
  old <- rpf_renamed(..., .env = environment())
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
    n_families = n_families,
    n_splits = n_splits,
    split_structure = split_structure,
    n_try_split_points = n_try_split_points,
    share_candidates = share_candidates,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    exact_step = exact_step,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment(), old),
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
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
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
  old <- rpf_renamed(..., .env = environment())
  check_outcome_type(y)
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    n_families = n_families,
    n_splits = n_splits,
    split_structure = split_structure,
    n_try_split_points = n_try_split_points,
    share_candidates = share_candidates,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    exact_step = exact_step,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment(), old),
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
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
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
  old <- rpf_renamed(..., .env = environment())
  check_outcome_type(y)
  blueprint <- hardhat::default_xy_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, y, blueprint = blueprint)
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    n_families = n_families,
    n_splits = n_splits,
    split_structure = split_structure,
    n_try_split_points = n_try_split_points,
    share_candidates = share_candidates,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    exact_step = exact_step,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment(), old),
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
  n_families = 200,
  n_splits = NULL,
  split_structure = "leaf",
  n_try_split_points = NULL,
  share_candidates = NULL,
  max_candidates = 50,
  split_decay_rate = 0.1,
  delete_leaves = TRUE,
  min_leaf_rows = 1,
  interaction_constraints = NULL,
  loss = NULL,
  exact_step = TRUE,
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
  old <- rpf_renamed(..., .env = environment())
  blueprint <- hardhat::default_recipe_blueprint(intercept = FALSE)
  processed <- hardhat::mold(x, data, blueprint = blueprint)
  # the recipe's steps decide the predictor columns: indices have no user-given order to refer to
  if (is.list(interaction_constraints) && any(vapply(interaction_constraints, is.numeric, logical(1)))) {
    cli::cli_abort("With a recipe, give {.arg interaction_constraints} by predictor names, not column indices.")
  }
  rpf_bridge(
    processed,
    max_interaction = max_interaction,
    n_families = n_families,
    n_splits = n_splits,
    split_structure = split_structure,
    n_try_split_points = n_try_split_points,
    share_candidates = share_candidates,
    max_candidates = max_candidates,
    split_decay_rate = split_decay_rate,
    delete_leaves = delete_leaves,
    min_leaf_rows = min_leaf_rows,
    interaction_constraints = interaction_constraints,
    loss = loss,
    exact_step = exact_step,
    prune = prune,
    max_panels = max_panels,
    tune = tune,
    given = tune_given(environment(), old),
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
  n_families,
  n_splits,
  split_structure,
  n_try_split_points,
  share_candidates,
  max_candidates,
  split_decay_rate,
  delete_leaves,
  min_leaf_rows,
  interaction_constraints,
  loss,
  exact_step,
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
  # same applies to values > p, silently, as the default of 3 exceeds p for one or two predictors
  if (max_interaction > p) {
    max_interaction <- p
  }

  checkmate::assert_int(n_families, lower = 1)
  checkmate::assert_int(n_splits, lower = 1, null.ok = TRUE)
  checkmate::assert_choice(split_structure, choices = unname(renamed_split_structures))
  checkmate::assert_int(n_try_split_points, lower = 1, null.ok = TRUE)
  checkmate::assert_number(share_candidates, lower = 0, upper = 1, null.ok = TRUE)
  # NULL greediness: 3 / 0.1 for "leaf", 10 / 0.4 (master's) for the other structures
  n_try_split_points <- n_try_split_points %||% if (split_structure == "leaf") 3 else 10
  share_candidates <- share_candidates %||% if (split_structure == "leaf") 0.1 else 0.4
  if (share_candidates <= 0) {
    cli::cli_abort("{.arg share_candidates} must be greater than 0, not {share_candidates}.", call = call)
  }
  checkmate::assert_int(max_candidates, lower = 1)
  checkmate::assert_number(split_decay_rate, lower = 0)
  checkmate::assert_flag(delete_leaves)
  checkmate::assert_int(min_leaf_rows, lower = 1)

  # "median" loss is implemented but discarded
  loss_functions <- switch(outcomes$mode, "regression" = "L2", "classification" = c("L1", "L2", "logit", "exponential"))
  checkmate::assert_choice(loss, choices = loss_functions)
  checkmate::assert_flag(exact_step)

  checkmate::assert_flag(prune, null.ok = TRUE)
  checkmate::assert_flag(deterministic)
  # prune = NULL: the pruned model for a random forest, the forest with deterministic = TRUE (no out-of-bag rows)
  prune <- prune %||% !deterministic
  # tune: the search chooses n_try_split_points, share_candidates, min_leaf_rows and n_splits, and prunes every forest
  control <- tune_resolve(tune, call = call)
  if (!is.null(control)) {
    tune_check_call(control, given, n_splits, prune, deterministic, n_families, call = call)
    n_splits <- NULL
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
  if (deterministic && n_families > 1) {
    cli::cli_warn(c(
      "With {.code deterministic = TRUE}, all {n_families} tree families are identical.",
      "i" = "Use {.code n_families = 1} to fit the same model faster."
    ))
  }
  # n_splits = NULL: chosen out of bag while growing (a random forest), otherwise 30
  adaptive <- is.null(n_splits) && !deterministic
  if (is.null(n_splits) && !adaptive) {
    n_splits <- 30L
  }
  checkpoints <- NULL
  if (!adaptive) {
    checkpoints <- check_oob_path(oob_path, n_splits, deterministic, call = call)
  } else if (!isTRUE(oob_path) && !isFALSE(oob_path)) {
    cli::cli_abort(
      c(
        "With {.code n_splits = NULL}, {.arg oob_path} must be {.code TRUE} or {.code FALSE}.",
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

  # the stored params keep randomPlantedForest 0.5.0's names and split structure codes (renamed_args)
  params <- list(
    loss = loss,
    ntrees = n_families,
    max_interaction = max_interaction,
    splits = if (adaptive) NA_integer_ else n_splits,
    split_try = n_try_split_points,
    t_try = share_candidates,
    split_decay_rate = split_decay_rate,
    max_candidates = max_candidates,
    delete_leaves = delete_leaves,
    min_leaf_rows = as.integer(min_leaf_rows),
    interaction_constraints = interaction_constraints,
    split_structure = split_structure_code(split_structure),
    deterministic = deterministic,
    nthreads = nthreads,
    # the pruning leaves the forest unpurified (prune_components()), so that rpf(prune = TRUE) is its pruning
    purify = purify && !prune,
    exact_draw = rpf_exact_draw(call = call),
    keep_level_order = predictors$keep_level_order
  )
  if (outcomes$mode == "classification") {
    # logit's steps (exponential's are exact anyway), and the losses' pseudo-rows, by an option for studies only
    params$exact_step <- exact_step
    params$prior_rows <- rpf_prior_rows(call = call)
  }

  if (!is.null(control)) {
    return(tune_run(processed, outcomes, predictors, params, control, max_panels, export_forest, call = call))
  }
  res <- rpf_fit(processed, outcomes, predictors, params, adaptive, checkpoints, export_forest, call = call)
  if (prune) {
    res <- rpf_prune(res, processed, max_panels, call = call)
  }
  res
}

# The forest of rpf_bridge()'s checked arguments: `adaptive` whether `n_splits` is chosen out of bag while growing,
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
    res$splits_oob <- list(n_splits = chosen$splits, stopped_at = chosen$stopped_at, capped = chosen$capped)
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

# The draw of candidates for split_structure = "leaf": the option
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

# The pseudo-rows a of the logit and exponential losses (1): the undocumented
# option randomPlantedForest.prior_rows, for studies; NULL if it is not set.
rpf_prior_rows <- function(call = rlang::caller_env()) {
  a <- getOption("randomPlantedForest.prior_rows")
  if (is.null(a)) {
    return(NULL)
  }
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a) || a <= 0) {
    cli::cli_abort("The option {.code randomPlantedForest.prior_rows} must be a positive number.", call = call)
  }
  as.double(a)
}

# Assemble the positional parameter vector for the C++ constructors.
# 15 elements for regression (the 13, min_leaf_rows, exact_draw), 17 for
# classification (+ exact_step, prior_rows). Models saved before these options
# existed have no such entries in their params and were fitted with
# min_leaf_rows = 1 and the exact draw; the delta and epsilon in the params of
# a classification model saved by an earlier version are not used.
rpf_param_vector <- function(params, mode) {
  # the stored code (renamed_split_structures); a new name is taken too
  structure <- params$split_structure
  if (structure %in% renamed_split_structures) {
    structure <- split_structure_code(structure)
  }
  split_mode <- switch(
    structure,
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
  base <- c(base, params$min_leaf_rows %||% 1L, params$exact_draw %||% TRUE)
  if (mode == "classification") {
    base <- c(base, params$exact_step %||% TRUE, params$prior_rows %||% 1)
  }
  base
}

# Main fitting function and interface to C++ implementation
# `oob_path`: NULL, or the checkpoints of the OOB path (see check_oob_path())
# `adaptive`: NULL, or the adaptive fit's checkpoints, stopper, patience, tolerance, floor and pause_gb (see
# tune_splits_oob())
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
      adaptive$tolerance,
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
