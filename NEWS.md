# randomPlantedForest (development version)

## New features

* `rpf()` gains `keep_level_order`, to keep the level order of factor predictors instead of ordering their levels by the outcome.
  `TRUE` keeps it for every factor, ordered or not, and a character vector for the named factors only; splits then separate neighbouring levels in the order of `levels()`.
  The default, `FALSE`, orders every factor by the outcome as before, and fits are unchanged.
  The coding is stored in `$factor_levels`, which `predict()`, `predict_components()` and `purify()` use, and `print()` names the factors that keep their order.
  The details of `?rpf` describe how factor predictors are coded.
* `plot()` draws the trees of one tree family to illustrate how a Random Planted Forest predicts.
  `type = "tree"` draws a split diagram in the style of Figure 1b of the paper, reconstructed from the fitted leaves, `type = "boxes"` draws the leaves in variable space, before or after purification.
  `rpf_boxes()` returns the leaves as a data frame.
  Requires `ggplot2`, and `patchwork` for `type = "boxes"`.
  The "Getting started" vignette uses it to show how a tree family is built and what purification does.
* `rpf()` gains `min_leaf_rows`, the minimum number of rows on each side of every split, so in every leaf.
  Values above 1 need regression with `split_structure = "leaves"`; the default 1 imposes no limit and is accepted everywhere.
  Leaves that cannot be split with this many rows on each side are no longer drawn as split candidates, and binary predictors are split whenever both sides keep enough rows.
  The default `min_leaf_rows = 1` fits the same models as before, bit for bit.
* `options(randomPlantedForest.exact_draw = FALSE)` draws the split candidates of regression with `split_structure = "leaves"` from a tree of weights kept up to date between rounds, instead of one random number per possible split in every round (`?rpf`, "The draw of candidates for `leaves`").
  This fits about 2.5 times as fast on larger data, for example 218 s → 86 s at 1 thread on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits).
  It draws from the same distribution and leaves out the same candidates under `min_leaf_rows`, but with other random numbers, so the same seed gives another forest.
  The default `TRUE` fits the same models as before, bit for bit.
* `rpf()` gains `interaction_constraints`, a list of feature sets that restricts which interactions the forest may form, as in xgboost: a tree on the predictors `T` is grown only if `T` is a subset of one of the sets.
  The sets are given by predictor names or column indices (with a formula, indices follow the formula's predictors; with a recipe, names only).
  They apply together with `max_interaction`, for every `split_structure`, in regression and classification, and `print()` shows them.
  The default `NULL` leaves fits unchanged.

## Bug fixes

* `split_structure = "cur_trees_1"` samples leaves in proportion to their available thresholds, as documented.
  After a split, the sampling weights of new leaves counted rows instead of distinct values, and the cached prefix sums were built incorrectly.
  Leaves with tied values were then drawn too often, and draws were wasted on leaves that cannot be split.
  With `deterministic = TRUE`, `"cur_trees_1"` evaluated only the first split candidate in every round; it now evaluates the first `max_candidates`, as the other split structures do.
  Fitted models change for regression with `split_structure = "cur_trees_1"`, random and deterministic.
* `split_structure = "res_trees"` with `deterministic = TRUE` evaluates thresholds evenly spread over all (predictor, leaf, threshold) triples of a candidate.
  Previously it tried only the lowest thresholds of each leaf and left most of the variance unexplained.
  Fitted models change for regression with `split_structure = "res_trees"` and `deterministic = TRUE` only.
* `split_structure = "cur_trees_2"` with `deterministic = TRUE` evaluates `split_try` thresholds evenly spread over each leaf, as `"cur_trees_1"` does.
  Previously it tried the 1st to 9th smallest values of each leaf.
  Fitted models change for regression with `split_structure = "cur_trees_2"` and `deterministic = TRUE` only.
* Classification's split structures were never completed to match regression's.
  They now use the same candidate pool, the same sampling of candidates, leaves and thresholds, and the same aging; only the loss differs.
  With `loss = "L2"` a classification fit on a 0/1 response now equals the regression fit (same seed).
  Previously `"cur_trees_1"` was `"cur_trees_2"`; `"res_trees"` used another candidate pool and ignored `split_decay_rate`; `"hist"` only ever split the root, so it fitted an additive model of stumps whatever `max_interaction` was; and `"cur_trees_2"` never split the candidate tree's own leaves.
  The fixes above for `"cur_trees_1"` sampling and for `deterministic = TRUE` apply to classification too.
  Classification fits change for every `split_structure` and every loss.
* The internal loss `"exponential_2"` no longer crashes.
* Splits can separate the largest value of a predictor within a leaf.
  Previously every split structure except `"leaves"` for regression excluded that split, so binary predictors were never split on, and classification could not learn from them at all.
  Fitted models change for classification and for regression with `split_structure` other than `"leaves"`.
* A tree family that makes no split predicts the mean outcome of its bootstrap sample, and for classification its class shares (truncated to `[epsilon, 1 - epsilon]` with `loss = "logit"` or `"exponential"`).
  Previously it predicted 0, a probability of 0.5: constant predictors, a single row or `min_leaf_rows` above half the rows gave an all-zero model, and families that could not split pulled the forest towards 0.
  `purify()` puts that value into the intercept.
  Fitted models change only where a family made no split.
* Errors while fitting or purifying with `nthreads > 1` are reported as R errors.
  Previously an error in a worker thread terminated the R session.
* `predict(type = "link")` works for classification forests with `loss = "L1"` or `"L2"`.
  It returns the raw predictions, as `type = "numeric"` does (the link is the identity).
  Previously it stopped with an error about the size of `pred`.
* `rpf()` stops with a clear error on inputs that gave wrong fits without a warning or failed later: an infinite outcome, `epsilon` outside `(0, 0.5)`, `delta` of 0.5 or more (or 0 with `loss = "logit"`), `t_try = 0`, predictor names that are missing, duplicated or contain `:`, and interaction terms in the formula.
  `rpf()` finds interactions itself, so write `y ~ x1 + x2` instead of `y ~ x1 * x2`; to fit a product as a predictor of its own, write `I(x1 * x2)` or add it to the data as a column.
  A character or logical outcome and a factor outcome with a single class get the package's own message instead of one from `hardhat` or `model.matrix()`.
  Unused levels of a factor outcome are dropped with a warning; previously a binary outcome with an unused level was fitted, but `predict()` failed.
  `predict()` on zero rows returns zero rows instead of an error.
* The C++ code compiles without warnings under `-Wall -pedantic`.
* `predict_components()` returns one intercept per target level for multiclass classification, named by `target_levels` (#70).
  Previously only the first level's intercept was returned, so components plus intercept did not sum to the prediction for other levels.
* `predict_components()` accepts a matrix as `new_data`, like `predict()` does.
  It no longer turns the caller's data.frame into a data.table by reference.
* `predict()` after `purify(maxp_interaction = k)` with `k` below the fit's `max_interaction` returns the forest's prediction.
  Previously it returned the intercept plus the components up to order `k`, leaving out the higher orders.
  `predict_components()` now returns their sum as `remainder`, so that its result sums to the prediction, and `rpf_marshal()` keeps the purification's cap.
* `purify()` with `mode = 2` no longer drops training rows that lie on a leaf's lower bound. Purified components and the intercept now equal those of `mode = 1`.
* `predict_components()` returns a `remainder` whenever `m` leaves out components of the forest, also when `predictors` leaves out some of the fit's predictors at the fit's `max_interaction`.
  Previously the intercept and `m` then did not sum to the prediction.
* `predict_components()` takes the components' predictors from the forest's own columns, in their order.
  A formula with terms such as `log(x1)`, `I(x1^2)` or `. - x3` no longer makes it fail, and with factor predictors the components come in the order of the fit's columns, which has the factors last.
* `predict_components()` on a `new_data` with zero rows returns zero rows with the right columns, instead of failing with "Feature vector is empty.".
* `print()` counts the forest's own predictor columns: with formula terms such as `log(x1)` or `. - x3` it showed the number of raw columns the formula reads.
* `purify()` on a forest purified with `maxp_interaction = k` purifies it again when asked for more orders (all orders by default).
  Previously it returned the forest unchanged, so the orders above `k` stayed zero in `predict_components()`.
  `purify()` also checks its arguments when the forest is already purified.
  A capped forest restored by `rpf_unmarshal()` without its training data cannot be purified again; asking for more orders is an error that says so.

## Documentation

* `?purify` states how the components are identified (the marginal identification over the training rows, Hiabu, Meyer and Wright, 2023), and that the intercept is the mean prediction on the training rows.
* `?rpf` says that `deterministic = TRUE` also ignores `split_decay_rate`, how `split_structure = "leaves"` draws its thresholds, and that a seed reproduces a fit on the same platform only.
* `?print.rpf_forest`'s example prints a forest (`export_forest = TRUE`), and `?plot.rpf` cites the RPF paper as 2020.

## Performance

* Classification with `loss = "logit"` fits 3x to 11x faster, with bitwise the same fits.
  For every candidate split, the logit loss copied the weights of all observations twice, although it reads only those of the leaf being split, so a candidate cost time in proportion to the sample size instead of the leaf size.
  With 300 rows and 150 splits (4 tree families, 1 thread), a binary fit takes 0.37 s instead of 1.9 s and a three-class fit 0.88 s instead of 2.9 s; with 1,000 rows and 500 splits, a binary fit takes 2.3 s instead of 26 s.
* Fitting regression forests with `split_structure = "leaves"` (the default) is about 1.4 times as fast on larger data, for example 310 s → 216 s at 1 thread and 50 s → 34 s at 8 threads on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits).
  Fitted models are bitwise unchanged.
* `predict()` is about 3x faster for forests that are not purified, with bitwise the same predictions.
  From 512 rows on, each tree of more than 64 leaves gets an index of its leaves (a k-d tree, built per call), so a row is tested only against the leaves that may hold it instead of against every leaf.
  On California housing (100 tree families, `max_interaction = 3`, 4,128 rows) prediction takes 1.9 s instead of 5.7 s at 1 thread and 0.45 s instead of 1.26 s at 8 threads; with `delete_leaves = FALSE`, 2.9 s instead of 10.7 s and 0.71 s instead of 2.07 s.
* Fitting regression forests with `split_structure = "cur_trees_1"` is about 6 times as fast and needs less memory, for example 1,299 s → 216 s at 1 thread and 238 s → 38 s at 8 threads on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits), with the peak memory at 8 threads down from 2.5 GB to 0.9 GB.
  Fitted models are bitwise unchanged.
* `purify()` keeps the purified components as weighted rectangles instead of dense grids over all split points (with `mode = 2`, the default, which `predict_components()` and `rpf(purify = TRUE)` use).
  Purification at `max_interaction = 3` and above no longer runs out of memory: on California housing (16,512 rows, 50 tree families, 94,138 leaves) order 3 takes 6 s and 0.25 GB with 4 threads, where it was stopped at 6 GB before, and the full order (8) takes 40 s and 2.3 GB.
  The components equal those of `mode = 1` up to rounding, and `rpf_marshal()` stores the rectangles.
* `predict_components()` evaluates all components in one pass over the rows of `new_data`, and skips those the forest does not hold, instead of one pass per component: on a 10-feature full-order fit (1,023 components, 636 of them not held), 5,000 rows take 0.5 s instead of 2.0 s.
* A fitted forest takes about 40% less memory: the split search's per-leaf caches are kept apart from the leaves and freed when a tree family is finished, so a leaf takes 80 instead of 232 bytes; finished classification forests now also drop their training rows and caches.
  On California housing (`max_interaction = 3`, 10,000 splits, 200 tree families) the fitted forest falls from 821 MB to 490 MB and the fit's peak resident memory from about 1.9 GB to 1.6 GB, in the same time.
  Fitted models are bitwise unchanged.

# randomPlantedForest 0.5.0

## Breaking changes

* `rpf()` now defaults to `max_interaction = 2`, main effects and pairwise interactions, instead of `1`, an additive model (#66).
  Values above the number of predictors are capped silently instead of with a message.
* `rpf()` now defaults to `loss = "exponential"` for classification instead of `"L2"`, so classification fits give proper probability estimates by default; regression still uses `"L2"` (#66).
* `rpf()` arguments are reordered by purpose (forest size, split search, loss, other). Pass arguments after the data arguments by name (#66).
* `rpf()`, `predict()` and `purify()` error on unknown arguments, such as misspelled ones, instead of silently ignoring them (#66).
* `purify()` is a regular function instead of an S3 generic (#66).
* `preprocess_predictors_predict()` is no longer exported; it is an internal step of `predict()` (#66).
* The package now requires R >= 4.1.0 (#66).
* For a given seed, classification models differ from 0.4.0 due to the fixes and speedup below, but regression models are unchanged (#57, #66).
* For direct users of the C++ object (`$fit`): `set_data()` no longer fits, call `fit()` afterwards, and an unknown classification loss is an error instead of a silent fallback to L2 (@jyliuu, #57).
* `rpf()` no longer has a `cv` argument, which never had an effect: the C++ cross-validation was a no-op (#66).

## Bug fixes

* `predict()` and `predict_components()` error on missing values in `new_data`, naming the affected columns.
  Previously such rows silently received arbitrary predictions, as `rpf()` has no handling for missing values.
  `rpf()` likewise reports missing values in the training data by column (#66).
* Classification fits with a fixed seed could differ between runs.
  When a split replaced a leaf, its new split candidates could be lost or attached to the wrong leaf, depending on memory layout (#66).
* `predict(type = "class")` breaks probability ties by level order instead of at random, so class predictions are deterministic and no longer advance the random number generator (#66).
* `purify()` returns its input invisibly, as documented (#66).
* Logical predictors in the formula interface are used as a single 0/1 column, as in the x/y interface.
  Previously they were expanded into two redundant indicator columns, which also counted against `max_interaction` (#66).

## Performance

* Prediction is much faster: leaves are scanned from flat, cache-friendly arrays (~15x single-threaded), and rows are split across threads (#65).
  `predict()` gains `nthreads`, defaulting to the `nthreads` used for fitting.
* `predict_components()` is much faster: purified forests look up only the trees of the requested component, and each distinct input row is predicted once, which pays off for low-cardinality features (#59, #65).
* Classification fits are slightly faster, as redundant work during model construction was removed (#57).

## Other improvements

* New "Getting started" vignette covering regression, classification, component decomposition, recipes, and saving models (#66).
* The `print()` method is overhauled and gains a `format()` counterpart: it reports the forest size, interaction degree, purification state, split search settings and number of threads (#66).
* `deterministic = TRUE` with `ntrees > 1` warns, as all tree families are then identical (#66).
* Rewritten documentation for `rpf()`, `purify()` and `predict_components()`, including guidance on choosing parameters and on how purification interacts with `predict_components()` (#66).

## Internals

* The C++ core no longer depends on `Rcpp`: it uses standard C++ types and exceptions, and an `Rcpp` layer (`src/rcpp_interface.*`) converts at the boundary.
  This is groundwork for bindings in other languages such as Python (@jyliuu, #57).
* `cli` and `rlang` are now imported, `mvtnorm` is only needed to build the `pkgdown` site (#66).

# randomPlantedForest 0.4.0

It is now possible to serialize and de-serialize an `rpf` object via marshalling,
meaning a fitted model can be saved to disk or dispatched to a worker process for 
parallelization or encapsulation, as is done in `mlr3`. This also re-opens the door 
for the `mlr3extralearners` wrapper, which was removed due to the lack of 
serialization support.

* New `rpf_marshal()` / `rpf_unmarshal()` serialize a fitted forest to a plain
  R list and back, making `saveRDS()`-based storage of rpf models possible (#52).
  Purified forests restore their purified state directly; training data is only
  embedded with `include_data = TRUE` (required to `purify()` after restoring).
  Blobs record the blob format version and the package version they were
  created with; restoring under an older package than the one that saved the
  model warns. Malformed or corrupt blobs error instead of reading out of bounds.
* New `rpf_is_valid()` checks whether an rpf object's internal model is usable,
  `predict()`, `purify()` and `predict_components()` now give an actionable
  error for `rpf` objects restored via `readRDS()` without marshaling.
* New `bundle::bundle()` method for rpf models wrapping the marshaling API.

## Behavior changes

* The default for `delta` changed from 0 to 0.001. With `delta = 0`, splits
  producing single-class nodes have infinite logit loss and are always rejected,
  which also degraded binary logit fits.

## Fixes and improvements

* Fixed multiclass classification with `loss = "logit"`, which effectively never
  split and predicted near-uniform class probabilities (#40). The C++ logit loss
  is a reference-class multinomial formulation expecting `K-1` indicator columns,
  but the R wrapper passed a full `K`-column one-hot matrix, pinning the implicit
  reference-class probability to zero. Multiclass logit outcomes are now encoded
  with the first factor level as reference class:
  * `predict(type = "prob")` reconstructs all `K` class probabilities from the
    `K-1` logits; rows sum to 1 exactly.
  * `predict(type = "numeric"/"link")` now returns `K-1` columns named after the
    non-reference levels (previously `K` columns).
  * `predict_components()` on multiclass logit fits returns per-class components
    for the `K-1` non-reference levels; `target_levels` reflects this.

# randomPlantedForest 0.3.0

## Major changes (#61)

* New `rpf()` arguments controlling split-candidate sampling:
  * `split_structure = "leaves"`: Defines what a split candidate is and how candidates are drawn.
    One of `"leaves"` (default), `"hist"`, `"cur_trees_1"`, `"cur_trees_2"`, or `"res_trees"`; see `?rpf` for details.
  * `max_candidates = 50`: Maximum number of split candidates sampled per iteration.
  * `split_decay_rate = 0.1`: Exponential aging of repeatedly drawn but unchosen split candidates.
    `split_decay_rate = 0` corresponds to no aging and uniform sampling.
  * `delete_leaves = TRUE`: Whether a parent leaf is deleted when splitting along an existing dimension.
* **Fitting results change**: The new candidate-sampling defaults and a reworked internal RNG
  mean that fits are not reproducible against previous versions, even with the same seed.
  Install an older commit if exact reproduction of previous results is required.
* Seeded fits are now reproducible regardless of `nthreads`: per-tree seeds are drawn from R's RNG,
  so `set.seed()` gives identical forests for serial and multithreaded fits.
* Substantial speedups in fitting (cached per-leaf orderings, prefix sums) and reduced memory use
  (training-only buffers are released after each tree family is built).
* `purify()` gains arguments:
  * `mode = 2`: Purification algorithm; `2` is a new fast exact method, `1` is the legacy grid-based path.
  * `nthreads = NULL`: Purification is now multithreaded, defaulting to the fit's `nthreads` setting.
  * `maxp_interaction = NULL`: Optionally only compute purified components up to this interaction order.
* New `rpf()` argument `export_forest = FALSE`: The flattened forest is no longer stored in the
  fitted object by default, so `rpf_object$forest` is `NULL` unless `export_forest = TRUE`.
  This reduces object size; `predict()`, `purify()`, and `predict_components()` are unaffected.
* `preprocess_predictors_predict()` is now exported.
* Fixed a memory bug in the legacy purification path where the grid was sized one element too large,
  causing out-of-bounds reads (crashes on Windows, silently wrong purification results elsewhere).
* Fixed a crash on Windows when fitting with `nthreads > 1`, caused by a `thread_local` buffer
  with a non-trivial destructor being destroyed at thread exit.

## Other changes

* Internals in `src/` have been refactored into modular sub-files (#53)
* `rpf()` now errors if a regression target is combined with a `loss` other than `"L2"`.
* Allow features of type `logical`, which are now converted via `as.integer`.
* The `parallel = TRUE|FALSE` argument in `rpf()` has been substituted by an `nthreads = 1L` argument, allowing for more flexible parallelization.
  The previous behavior only allowed for either no parallelization or using n-1 of n available cores. 
  The new implementation should be reasonably robust and the default behavior remains serial execution.
* Remove `SystemRequirements` field from `DESCRIPTION`: Now the default C++ version is C++17 and 
  with a minor change to internal use of random numbers, `randomPlantedForest` is now compatible with C++11 through C++23.
* Add `remainder` term to `predict_components` output for case where `max_interaction` supplied is smaller than `max_interaction` in `rpf` fit.
  In that case, the `m` values don't sum up to the global predictions, so we add a remainder to allow reconstruction of that property.

# randomPlantedForest 0.2.1

* Add `glex` class to output of `predict_components()`, for extended functionality available with [`glex`](https://github.com/PlantedML/glex).
* Add `target_levels` vector to output of `predict_components()` to aid multiclass handling.
Keeping track of levels is somewhat awkward since column names of `$m` need to be identifiable
regarding the target level.

# randomPlantedForest 0.2.0

* Added a `NEWS.md` file to track changes to the package.
