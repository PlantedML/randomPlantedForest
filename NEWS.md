# randomPlantedForest (development version)

## Breaking changes

* `rpf()` has new defaults: a plain `rpf(y ~ ., data)` now returns the pruned model (class `rpf_pruned`), for regression and classification.
  The new defaults, each switched by its argument:
  - `n_splits = NULL`: the number of splits is chosen out of bag while the forest grows (was 30).
  - `prune = NULL`: the forest is pruned by `prune_components()` with its defaults, and every kept component of three or more predictors is shown in at most `max_panels = 10` heat maps (was no pruning). With `deterministic = TRUE` it means `prune = FALSE` and `n_splits = NULL` means 30, as there are no out-of-bag rows.
  - `n_try_split_points = NULL, share_candidates = NULL`: 3 thresholds per candidate and 10% of the possible splits per round for `split_structure = "leaf"` (was 10 / 0.4, still the default of the other split structures). On 48 simulated and benchmark data sets this greediness had 5.7% to 6.8% lower test error than 10 / 0.4.
  - The split candidates of `"leaf"` are drawn by the fast draw (`options(randomPlantedForest.exact_draw = FALSE)`, now the default): the same distribution as the exact draw, with other random numbers, so the same seed gives another forest. It is about 2.5 times as fast on larger data.
  - `max_interaction = 3` (was 2). Order 3 had 15% lower test error than order 2 on simulated data and 1.5% on California housing; no limit (order p) was worse on California housing and 11 to 41 times slower.
  - `n_families = 200` (was 50): more families had lower test error at every order (2% per doubling) and sparser pruned models.
  The loss is `"L2"` for regression and, as below, now `"logit"` for classification.
  So every default fit changes, regression and classification, and `rpf()` returns an `rpf_pruned` object instead of an `rpf` unless `prune = FALSE`.
  To fit a regression forest as before, give master's settings and the exact draw:
  ```r
  options(randomPlantedForest.exact_draw = TRUE)
  rpf(..., max_interaction = 2, n_families = 50, n_splits = 30, n_try_split_points = 10, share_candidates = 0.4, prune = FALSE)
  ```
  Classification fits of earlier versions cannot be reproduced this way: they change for every split structure and loss (see Bug fixes: classification uses regression's split structures, the largest value of a predictor can be split off, and three or more classes with the exponential loss use SAMME).
  The default call does more work than before: 0.2 s on `mtcars`, 0.4 s on `iris`, and about 2.3 minutes with a peak of 4.6 GB of resident memory (RSS; 2.7 GB of physical footprint, which varies by a few tenths of a GB from run to run) on California housing (16,512 rows) at 8 threads, measured with all changes of this version.

* The classification losses `"logit"` and `"exponential"` (binary, and SAMME for three or more classes) no longer clamp class shares, and `rpf()`'s arguments `delta` and `epsilon` are removed: a call that gives either is an error (an unused argument).
  Each new leaf's step now minimises the leaf's penalised loss: the loss of its rows, plus that of one pseudo-row, split evenly over the classes, at the leaf's current mean score (`?rpf`, "Classification losses").
  A split is scored by its gain over the parent leaf's own step, and each new leaf pays for its pseudo-row.
  So a leaf of one class converges to a finite level (`log(2 n + 1)` for `n` rows and two classes) and is not moved again, instead of being split and pushed further in every round (scores of several hundred after 300 rounds); splitting such a leaf, or repeating a split whose leaves are at their optimum, gains nothing and is not made; and a round in which no candidate gains makes no split.
  Previously the split criterion credited a split with the parent leaf's own step at a class share clamped to `[delta, 1 - delta]`, while the update clamped it to `[epsilon, 1 - epsilon]`, so pure and nearly pure leaves won over splits that separate classes, and logit's leaves reset the mean score to the log-odds of the raw share, which is not a step of the loss and left fits underconfident.
  With the defaults (`n_splits` out of bag, pruning), on six data sets (penguins, a 300-row and a 4,000-row simulated set, each with two and three classes; 3 seeds), the test log-loss of the exponential loss changes by -0.018 to +0.010 on five sets and by +0.048 on the 300-row three-class set (worse in all three seeds), and that of the logit loss by -0.017 to +0.016.
  At fixed split counts the exponential loss no longer overfits (300 splits: 0.58 against 3.42 on the 300-row three-class set, 0.08 against 1.00 on two-class penguins); on the noisy 4,000-row sets at 300 splits both fit the noise more than the clamped logit loss did (logit +0.03 and +0.05, exponential +0.05 and +0.13; all still below the old exponential loss), which the out-of-bag choice of `n_splits` avoids.
  `rpf()` gains `exact_step` (default `TRUE`): logit's two new leaves get the exact minimisers of their penalised loss (Newton's method to convergence on the leaf's rows); `FALSE` takes one Newton step.
  The exponential loss's step is exact by itself, so `exact_step` does not change its fits.
  The split search scores a candidate from sums over its two sides built as regression builds them (along the leaf's sorted rows, or per bin), so classification fits are faster: at 1 thread on 10,000 rows (100 splits, 20 tree families) a binary fit takes 3.8 s instead of 7.9 s with the exponential loss and 3.6 s instead of 9.7 s with logit (regression: 2.7 s), a three-class fit 5.2 s instead of 7.1 s and 5.1 s instead of 11.5 s.
  Classification fits with `"logit"` or `"exponential"` change; `"L1"` and `"L2"` fits and regression are unchanged, and a model saved before still loads, prints and predicts.

* The default classification loss is now `"logit"` (`loss = NULL`), with two classes and with more; it was `"exponential"` (with more than two classes SAMME).
  So every default classification fit changes; `loss = "exponential"` gives the previous default, which stays available and unchanged.
  On 8 public data sets (4 binary, 4 with 3 to 7 classes; the default pipeline, 5 repeated 75 / 25 splits each), the test log-loss of `"logit"` was lower than or tied with that of `"exponential"` on every set: a tie on the binary sets (+0.0004 on average) and 0.053 lower with 3 or more classes (lower on 2 of the 4 sets, tied on the others).
  With many classes `"exponential"` can need many more splits: on Dry Bean (7 classes, 2,000 rows) the out-of-bag choice was 1,345 to 4,525 splits (310 to 990 s per fit), where `"logit"` chose 22 on average.

* A pruned model no longer holds the forest it was pruned from: `rpf()` with pruning (the default) and `prune_components()` return the model without it, and their new argument `keep_forest = TRUE` keeps it.
  The model predicts, prints and gives its components from what it stores, bit for bit as with the forest, so `predict()` (every `type`), `predict_components()`, `print()`, the display, `rpf_marshal()` and `bundle()` work as before; the fit itself is unchanged.
  What is lost by default is the forest in `$fit`: pruning it again (`prune_components(model$fit, data, ...)`), `purify()`, and the forest's own `predict()`, `predict_components()`, `plot()` and `rpf_boxes()` are an error that names `keep_forest = TRUE`. `$fit` still describes the forest (`$params`, `$oob_path`, `$splits_oob`, the blueprint), and `prune_components()`' caller still has the forest it passed.
  To get the forest back, fit with `keep_forest = TRUE`, or with `prune = FALSE` and prune with `prune_components()`. With `prune = FALSE` the argument has no effect.
  On California housing (16,512 rows, 8 threads) the forest is 0.6 GB of memory, which is free when `rpf()` returns, and the model 0.19 GB; the model without its forest is plain R data, so `saveRDS()` and `readRDS()` give a model that predicts (6 s, 113 MB), and `rpf_marshal()` takes no time and gives a blob of the model's size, against 13 s and 1.15 GB with the forest.
  Models saved with their forest by an earlier version (`rpf_marshal()`) still load, print and predict, and keep their forest.

## Renamed arguments

* Four arguments of `rpf()` and the split structures have clearer names; the prefix gives the unit (`n_` a count, `share_` a share in (0, 1], `max_` a cap):
  - `ntrees` → `n_families` (the number of tree families),
  - `splits` → `n_splits` (splits per tree family),
  - `split_try` → `n_try_split_points` (split points tried per candidate),
  - `t_try` → `share_candidates` (share of the possible splits drawn as candidates each round);
  - `split_structure`: `"leaves"` → `"leaf"`, `"hist"` → `"leaf_binned"`, `"cur_trees_1"` → `"tree_then_leaf"`, `"cur_trees_2"` → `"tree_then_all_leaves"`, `"res_trees"` → `"target_tree_then_leaf"` (the candidate, then how the leaf to split is found).
* The old names still work and give the same fit, with a deprecation warning once per session for each old name (class `lifecycle_warning_deprecated`; `options(lifecycle_verbosity = "warning")` warns every time, `"quiet"` never); giving an argument under both names is an error.
  The same holds for `oob_path_checkpoints(splits = )` and for a `tune_control()` ladder with the columns `split_try` and `t_try`.
* The new names are used in `print()`, in the tables of the new features (the `n_splits` column of `$oob_path`, the tuning's `$tune$search` and `$tune$settings`, `greediness_ladder()`), and in errors and documentation.
  A fitted model's `$params` keeps the old names (`ntrees`, `splits`, `split_try`, `t_try` and the old `split_structure` values), so code that reads them, and models saved by 0.5.0, keep working.

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
* `prune_components()` prunes the purified components of a forest (regression or classification; see the classification entry below).
  A component is kept if its share of the variance of the prediction on the training rows is at least `tau`.
  Rule `"D"` drops the others; rule `"R"` also refits the intercept and one coefficient per kept component by least squares, which changes their scale but not their shape.
  The rule and `tau` are chosen out of bag by default (`select_by = "oob"`), without refitting: one purification evaluates every tree family's components at the training rows outside its bootstrap sample, and each candidate is scored by its out-of-bag RMSE, corrected for the extra variance of averaging fewer families.
  On 52 data sets (44 synthetic, 8 splits of California housing) the out-of-bag choice picked the cross-validated model on 43, its test error was 0.09% higher on average, and it took a twelfth of the time.
  `select_by = "cv"` chooses by K-fold cross-validation instead, refitting the forest on each fold; by default it is used below 10 tree families and for deterministic forests.
  By default the largest `tau` whose error is within one standard error of the lowest is chosen (`select = "1se"`; `select = "min"` takes the lowest).
  With `select_by = "oob"`, the coefficients of rule `"R"` are fitted on the out-of-bag components, to score each candidate, for the final model and for a given rule and `tau`, because the forest's own components have also fitted the training rows' noise; with `select_by = "cv"` they are fitted in sample.
  On the same 52 data sets this sparse default had 1.3% lower test error than the lowest out-of-bag error with in-sample coefficients, with half the components (29 instead of 83 on California housing) and no false component on the synthetic sets; it can drop a small true component (a four-way interaction in one synthetic design, +2.5% test error there).
  By default every kept component of three or more predictors is then shown in at most `max_panels = 10` heat maps, and the model is pruned again with the shown components (`max_panels = Inf` switches this off).
  Two predictors of such a component stay at full resolution; the others are cut into intervals of equal numbers of training rows whose counts multiply to at most `max_panels` (a count of 1 averages a predictor out), the component is averaged within the cells (weighted by the product of the coarse predictors' training marginals; the result is used as it is, not purified again), and the cuts that keep most of its variance are taken.
  The same rule then prunes again, so that the pruning judges the model that is shown; predictions use the coarsened components, and `tau` is then a share of that model's variance.
  In simulations (48 data sets) this had 0.49% lower test error than coarsening after pruning, and 0.57% lower than not coarsening; on California housing (50 families, order 3) the default keeps 30 components, shown in 75 heat maps.
  Every fit records how often each training row is in each family's bootstrap sample (fits are unchanged), and `rpf_marshal(include_data = TRUE)` stores these in-bag counts.
  On 44 synthetic data sets with known sparse structure, the pruned models kept no false component (except five tiny main effects on one additive set) and 111 of the 113 true interactions, and their test error was 15% lower than that of the unpruned forest tuned by CV, lower on 42 of the 44 sets.
  The result, of class `rpf_pruned`, has methods for `predict()`, `predict_components()` (now an S3 generic) and `print()`, and works with `rpf_marshal()`, `rpf_unmarshal()` and `bundle()`.
* `rpf(oob_path = TRUE)` gives the out-of-bag error along the splits from one fit, in `$oob_path` (for classification the loss's criterion, see below).
  With the same seed, the forest with `k` splits is exactly the first `k` splits of a forest with more, so while each tree family grows, its prediction at the training rows outside its bootstrap sample is updated split by split.
  The table has the plain and the corrected out-of-bag RMSE (the correction of `prune_components()`) at four checkpoints per doubling of the splits (`oob_path_checkpoints()`), or at given ones.
  On 52 data sets, choosing `n_splits` by the lowest corrected error among the values of a cross-validation grid did as well as 5-fold cross-validation over that grid (test error 0.06% lower on average), at 0.9 of its cost.
  The fitted forest is the same as without `oob_path`, which is off by default.
* `rpf(n_splits = NULL)` chooses the number of splits out of bag while the forest grows (regression and classification).
  The forest's corrected out-of-bag error is followed at four checkpoints per doubling of the splits, growing stops once its lowest value is one doubling back (at most 25,600 splits), and the forest at the lowest error is returned: the forest that `rpf()` with that many `n_splits` and the same seed gives.
  On 52 data sets this chose the same `n_splits` as the lowest error along a fixed budget of 3,200 splits on 50, and a later, lower minimum on the other two, while growing 0.72 of the splits.
  The tree families grow one at a time, `nthreads` together, and keep only their out-of-bag predictions and a 40-byte record per split, from which the chosen forest is rebuilt without a refit; `options(randomPlantedForest.adaptive_pause_gb)` (0.25 GB) bounds the memory of families waiting for the stop.
  `print()` shows the chosen `n_splits` and where growing stopped.
  It is the new default (see Breaking changes).
  A corrected out-of-bag error at or below 0 at a checkpoint (a correction that overshot; never observed) is not chosen unless all are.
  A family waiting for the stop pauses at most one doubling past the first families' estimate and grows on if the stop comes later, so a flat out-of-bag curve no longer makes every family grow far past the stop: a three-class fit with `loss = "logit"` and otherwise the defaults on 100 rows took 12 s instead of 35 s, with the same forest (measured when the pause was added, before the plateau stop below).
  Growing also stops once the error is on a plateau: once its lowest value has improved over the last doubling by at most one standard error of that improvement (paired over the training rows; `options(randomPlantedForest.adaptive_tolerance)`, 1 by default, 0 for the stop on a new minimum alone).
  The choice is still the lowest error among the splits grown, so the plateau stop changes the forest only where the lowest error would have come later.
  On recorded paths it chose the same `n_splits` as the stop on a new minimum alone on all 44 synthetic regression sets, 4 splits of California housing and 24 classification fits (logit and exponential loss, binary and three classes), and grew 10% (regression) to 29% (classification) fewer splits.
  The families grow only to about the plateau stop (on 33 of these fits, 4 threads): the adaptive fit took 19% less time with logit, 4% less with the exponential loss, 11% less on 8 synthetic regression sets, and the same time on California housing, where the families waiting for the stop do not fit in the pause budget and grow on to the estimate of the stop on a new minimum, as a regrowth from their seed would cost more.
  A three-class fit on 100 rows with `loss = "logit"` and otherwise the defaults chooses 119 splits instead of 238 and takes 2.5 s instead of 8.7 s at 1 thread; with `loss = "exponential"` (then the default) it chooses 11 splits either way and stops growing at 15 instead of 21 (both measured with all changes of this version).
* `rpf(prune = TRUE)` returns the pruned model instead of the forest (regression and classification): the forest is pruned by `prune_components()` with its defaults on the training rows, including the display at `max_panels` (a new argument of `rpf()`, 10 by default, `Inf` for no coarsening).
  The result equals `prune_components(rpf(..., prune = FALSE), data)` bit for bit from the same random number state, for the formula, x/y and recipe interfaces.
  With `n_splits = NULL`, one call chooses the number of splits out of bag, prunes and returns the model to show.
  `prune = TRUE` is an error with `deterministic = TRUE`, and `purify` is not used with it.
  It is the new default (`prune = NULL`, see Breaking changes); `prune = FALSE` returns the forest as before.
* `rpf()` gains `min_leaf_rows`, the minimum number of rows on each side of every split, so in every leaf.
  Every split structure takes it, in regression and classification; the default 1 imposes no limit.
  Split points are tried only among those that keep this many rows on each side, and leaves without one are no longer drawn: as split candidates with `"leaf"` and `"leaf_binned"`, by the leaf draw of `"tree_then_leaf"` and `"target_tree_then_leaf"` (a leaf's weight is its number of such split points), and `"tree_then_all_leaves"` skips them.
  Binary predictors are split whenever both sides keep enough rows.
  The default `min_leaf_rows = 1` fits the same models as before, bit for bit.
* Classification with `split_structure = "leaf"` (the default) takes the fast draw of split candidates (`options(randomPlantedForest.exact_draw = FALSE)`) as regression does, and `min_leaf_rows` with every split structure.
  Both use regression's draw: the same candidates are left out, and a classification fit with `loss = "L2"` on a 0/1 response is the regression fit, bit for bit, with either option.
  With `min_leaf_rows = 1` and the exact draw (the default before this version) classification fits are unchanged, bit for bit.
* `rpf(tune = TRUE)` tunes out of bag in one call (regression and classification): the greediness of the split search (`n_try_split_points` and `share_candidates`, on the ladder of the new `greediness_ladder()`) and `min_leaf_rows`.
  Every forest of the search has `n_splits` chosen out of bag and is pruned by `prune_components()` with its defaults, and the forest whose pruned model has the lowest corrected out-of-bag error wins; the display (`max_panels`) is applied once, to the chosen model.
  The search is coordinate-wise, about 10 forests: the ladder from 10 / 0.4 to 1 / 0.0125, then `min_leaf_rows` 1, 5, 20 and 100, then the two ladder neighbours of the chosen point.
  `tune = "quick"` grows two forests, 3 / 0.1 and 1 / 0.025, and the new `tune_control()` sets your own search.
  The result is the chosen forest's pruned model, with the search table in `$tune`.
  On 48 data sets (44 synthetic, 4 splits of California housing) the full search had 10.8% lower test error than 10 / 0.4 and 2.8% lower than 5-fold cross-validation over the same settings, at 1/24 of its cost; the quick preset had 9.4% lower test error. Those studies scored each forest with the pruning's choice of the time (the lowest out-of-bag error, R's coefficients in sample), not with the sparse default used here.
  The default stays `tune = FALSE`.
* Classification forests go through the same pipeline as regression forests: `rpf(oob_path = )`, `rpf(n_splits = NULL)`, `prune_components()` with its display (`max_panels`), `rpf(prune = TRUE)` and `rpf(tune = )`.
  Every out-of-bag choice is scored by the loss's own criterion, corrected for the out-of-bag family count: for `loss = "L1"` and `"L2"` the squared error of the raw scores summed over the classes (the Brier score), for `"logit"` and `"exponential"` the log-loss of the probabilities that `predict(type = "prob")` returns, with the correction in its second-order form.
  On simulated data at 100 tree families the corrected out-of-bag log-loss was within -0.8% to +1.4% of the test log-loss (mean per setting), the uncorrected one +0.6% to +4.1% above it.
  A component is kept or dropped for all classes, sized by its variance summed over the score columns centred across the classes; rule `"R"` fits one coefficient per kept component, shared by the classes, and one intercept per score column (least squares for the Brier score, maximum likelihood for the log-loss), out of bag; the display chooses one panel form per component for all classes.
  The pruned model's `predict()` takes `type` as `predict.rpf()` does, and its `predict_components()` has one column per component and class.
  With two classes and `loss = "L2"`, every choice is that of the regression pipeline on the 0/1 response, bit for bit.
  With the new defaults (see Breaking changes) a classification fit takes all these steps.
* `options(randomPlantedForest.exact_draw = FALSE)` draws the split candidates of regression with `split_structure = "leaf"` from a tree of weights kept up to date between rounds, instead of one random number per possible split in every round (`?rpf`, "The draw of candidates for `leaves`").
  This fits about 2.5 times as fast on larger data, for example 218 s → 86 s at 1 thread on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits).
  It draws from the same distribution and leaves out the same candidates under `min_leaf_rows`, but with other random numbers, so the same seed gives another forest.
  `TRUE` (the default before this version) fits the same models as before, bit for bit.
* `rpf()` gains `interaction_constraints`, a list of feature sets that restricts which interactions the forest may form, as in xgboost: a tree on the predictors `T` is grown only if `T` is a subset of one of the sets.
  The sets are given by predictor names or column indices (with a formula, indices follow the formula's predictors; with a recipe, names only).
  They apply together with `max_interaction`, for every `split_structure`, in regression and classification, and `print()` shows them.
  The default `NULL` leaves fits unchanged.
* `print()` of a pruned model shows the final model: its intercept and its components, largest first, each with its importance, its share of the variance of the model's prediction on the training rows (`$components$model_share`, a new column: the coefficient times the component, in the form the model uses, over the model's prediction; for classification over the centred score columns).
  `share`, what `tau` thresholds, is the forest's component before the coefficient: a three-predictor component of share 3.5% with coefficient 3.0 is 19% of the model.
  How the components were chosen (the forest's families and splits, given or chosen out of bag, the rule and threshold, what was left out, the heat-map forms) follows under "How it was chosen". No choice or fit changes.

## Bug fixes

* `split_structure = "tree_then_leaf"` samples leaves in proportion to their available thresholds, as documented.
  After a split, the sampling weights of new leaves counted rows instead of distinct values, and the cached prefix sums were built incorrectly.
  Leaves with tied values were then drawn too often, and draws were wasted on leaves that cannot be split.
  With `deterministic = TRUE`, `"tree_then_leaf"` evaluated only the first split candidate in every round; it now evaluates the first `max_candidates`, as the other split structures do.
  Fitted models change for regression with `split_structure = "tree_then_leaf"`, random and deterministic.
* `split_structure = "target_tree_then_leaf"` with `deterministic = TRUE` evaluates thresholds evenly spread over all (predictor, leaf, threshold) triples of a candidate.
  Previously it tried only the lowest thresholds of each leaf and left most of the variance unexplained.
  Fitted models change for regression with `split_structure = "target_tree_then_leaf"` and `deterministic = TRUE` only.
* `split_structure = "tree_then_all_leaves"` with `deterministic = TRUE` evaluates `n_try_split_points` thresholds evenly spread over each leaf, as `"tree_then_leaf"` does.
  Previously it tried the 1st to 9th smallest values of each leaf.
  Fitted models change for regression with `split_structure = "tree_then_all_leaves"` and `deterministic = TRUE` only.
* Classification's split structures were never completed to match regression's.
  They now use the same candidate pool, the same sampling of candidates, leaves and thresholds, and the same aging; only the loss differs.
  With `loss = "L2"` a classification fit on a 0/1 response now equals the regression fit (same seed).
  Previously `"tree_then_leaf"` was `"tree_then_all_leaves"`; `"target_tree_then_leaf"` used another candidate pool and ignored `split_decay_rate`; `"leaf_binned"` only ever split the root, so it fitted an additive model of stumps whatever `max_interaction` was; and `"tree_then_all_leaves"` never split the candidate tree's own leaves.
  The fixes above for `"tree_then_leaf"` sampling and for `deterministic = TRUE` apply to classification too.
  Classification fits change for every `split_structure` and every loss.
* The internal loss `"exponential_2"` no longer crashes.
* Splits can separate the largest value of a predictor within a leaf.
  Previously every split structure except `"leaf"` for regression excluded that split, so binary predictors were never split on, and classification could not learn from them at all.
  Fitted models change for classification and for regression with `split_structure` other than `"leaf"`.
* A tree family that makes no split predicts the mean outcome of its bootstrap sample, and for classification its rows' fit (with `loss = "logit"` or `"exponential"` the class shares smoothed by the losses' pseudo-row, see Breaking changes).
  Previously it predicted 0, a probability of 0.5: constant predictors, a single row or `min_leaf_rows` above half the rows gave an all-zero model, and families that could not split pulled the forest towards 0.
  `purify()` puts that value into the intercept.
  Fitted models change only where a family made no split.
* The same holds when `n_splits = NULL` (the default) chooses the number of splits out of bag: a family without a split up to the chosen number predicts its rows' fit instead of 0, as a fit with that `n_splits` does.
  The pruned model already got the level from its intercept; the forest in it (`$fit`) did not.
  Without `n_splits = NULL`, and for families with a split, fits are unchanged.
  The pruned model also no longer fails for classification when no component has any variance (for example, constant predictors).
* Errors while fitting or purifying with `nthreads > 1` are reported as R errors.
  Previously an error in a worker thread terminated the R session.
* `predict(type = "link")` works for classification forests with `loss = "L1"` or `"L2"`.
  It returns the raw predictions, as `type = "numeric"` does (the link is the identity).
  Previously it stopped with an error about the size of `pred`.
* `rpf()` stops with a clear error on inputs that gave wrong fits without a warning or failed later: an infinite outcome, `t_try = 0`, predictor names that are missing, duplicated or contain `:`, and interaction terms in the formula.
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
* `loss = "exponential"` is with more than two classes the multiclass exponential loss SAMME (Zhu, Zou, Rosset and Hastie, 2009): the split criterion, the leaf step and the weights all come from the one loss `exp(-y.f / K)`.
  Each row has one weight, `exp(-g)` with `g` its class's score, and a leaf moves the scores by `log(V_k) - mean_j log(V_j)`, with `V` its class weight sums plus those of its pseudo-row (see the classification losses under Breaking changes); `predict()` takes the softmax of the scores as before.
  Previously each class had its own one-vs-rest weights, and its step was `log(p_k / p_0)` with `p_0 = 1 - sum_k p_k` mixing shares under the different weights; it is 0 at the start and was clamped to `epsilon`, so every class's score rose by a common `log(1 / epsilon)` that the softmax ignores but the weights did not.
  The weights then no longer matched the probabilities `predict()` gives (by factors up to 10^10 in small tests).
  With `n_splits` chosen out of bag (`n_splits = NULL`) and the forest pruned (`prune = TRUE`; order 3, 200 families, `n_try_split_points = 3`, `share_candidates = 0.1`, the fast draw) on two simulated three-class data sets, the pruned model's test log-loss was 0.022 and 0.036 above logit's; it is now 0.061 below and 0.003 above, with 15 and 25 splits chosen instead of 21–25 and 238–400, and on the larger set (4,000 rows) the fit is four times as fast.
  Fitted models change for `loss = "exponential"` with more than two classes only; binary fits are unchanged.
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
  `splits` counts rounds: near saturation a round can make no split. Ties between equally good candidates may go to different candidates under the exact and the fast draw, so their fits differ in distribution near saturation, though not in accuracy.
* `?print.rpf_forest`'s example prints a forest (`export_forest = TRUE`), and `?plot.rpf` cites the RPF paper as 2020.

## Performance

* Classification with `loss = "logit"` fits 3x to 11x faster, with bitwise the same fits.
  For every candidate split, the logit loss copied the weights of all observations twice, although it reads only those of the leaf being split, so a candidate cost time in proportion to the sample size instead of the leaf size.
  With 300 rows and 150 splits (4 tree families, 1 thread), a binary fit takes 0.37 s instead of 1.9 s and a three-class fit 0.88 s instead of 2.9 s; with 1,000 rows and 500 splits, a binary fit takes 2.3 s instead of 26 s.
* Fitting regression forests with `split_structure = "leaf"` (the default) is about 1.4 times as fast on larger data, for example 310 s → 216 s at 1 thread and 50 s → 34 s at 8 threads on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits).
  Fitted models are bitwise unchanged.
* `predict()` is about 3x faster for forests that are not purified, with bitwise the same predictions.
  From 512 rows on, each tree of more than 64 leaves gets an index of its leaves (a k-d tree, built per call), so a row is tested only against the leaves that may hold it instead of against every leaf.
  On California housing (100 tree families, `max_interaction = 3`, 4,128 rows) prediction takes 1.9 s instead of 5.7 s at 1 thread and 0.45 s instead of 1.26 s at 8 threads; with `delete_leaves = FALSE`, 2.9 s instead of 10.7 s and 0.71 s instead of 2.07 s.
* Fitting regression forests with `split_structure = "tree_then_leaf"` is about 6 times as fast and needs less memory, for example 1,299 s → 216 s at 1 thread and 238 s → 38 s at 8 threads on California housing (16,512 rows, `max_interaction = 3`, 10,000 splits), with the peak memory at 8 threads down from 2.5 GB to 0.9 GB.
  Fitted models are bitwise unchanged.
* `purify()` keeps the purified components as weighted rectangles instead of dense grids over all split points (with `mode = 2`, the default, which `predict_components()` and `rpf(purify = TRUE)` use).
  Purification at `max_interaction = 3` and above no longer runs out of memory: on California housing (16,512 rows, 50 tree families, 94,138 leaves) order 3 takes 6 s and 0.25 GB with 4 threads, where it was stopped at 6 GB before, and the full order (8) takes 40 s and 2.3 GB.
  The components equal those of `mode = 1` up to rounding, and `rpf_marshal()` stores the rectangles.
* `predict_components()` evaluates all components in one pass over the rows of `new_data`, and skips those the forest does not hold, instead of one pass per component: on a 10-feature full-order fit (1,023 components, 636 of them not held), 5,000 rows take 0.5 s instead of 2.0 s.
* Classification uses the split search of regression: one implementation of the five split structures, which differ between the two only in how a candidate split is scored.
  So classification gets the faster `"tree_then_leaf"` and `"leaf"` and the smaller finished forests of regression, with bitwise the same fits.
  On simulated data with 16,512 rows and 8 predictors (`max_interaction = 3`, 10,000 splits, 50 tree families, 8 threads, binary, `loss = "exponential"`), `"tree_then_leaf"` takes 759 s instead of 1,044 s with a peak memory of 1.2 GB instead of 2.9 GB, and `"leaf"` 97 s instead of 113 s.
  With three classes the peak memory of `"tree_then_leaf"` falls about threefold (2.4 GB to 0.8 GB at 1,000 splits).
* `prune_components()`, and so `rpf(prune = TRUE)` and `rpf(tune = )`, no longer build the purified store of all components.
  Each tree family is purified in turn, its components are evaluated at the training rows (and out of bag) and then freed, and the pruned model stores the pieces of its kept components only, from which `predict()` and `predict_components()` evaluate them.
  `purify()` keeps every piece of every component, which grows steeply with the order: on California housing (16,512 rows, 200 tree families, `n_try_split_points = 3`, `share_candidates = 0.1`, the fast draw, 8 threads; both builds measured when this change was made) the whole call `rpf(n_splits = NULL, prune = TRUE)` at `max_interaction = 8` (all predictors) was stopped at 12 GB of resident memory (RSS) while purifying, and now finishes with a peak physical footprint of 6.4 GB (9.9 GB RSS), most of it the coarsened components of the display; at order 3 the footprint peaks at 2.3 GB instead of 3.2 GB (RSS 4.1 GB instead of 5.4 GB), in the same time.
  The pruned models are the same up to rounding: on 136 pruned models of 36 forests (simulated and benchmark data, orders 2 to 10) the choices, the kept components and the display are the same, and test predictions agree to 3e-15 of their scale.
  The forest in the result is left unpurified (call `purify()` for all its components), and an earlier purification of it is no longer redone.
* A fitted forest takes about 40% less memory: the split search's per-leaf caches are kept apart from the leaves and freed when a tree family is finished, so a leaf takes 80 instead of 232 bytes; finished classification forests now also drop their training rows and caches.
  On California housing (`max_interaction = 3`, 10,000 splits, 200 tree families) the fitted forest falls from 821 MB to 490 MB and the fit's peak resident memory from about 1.9 GB to 1.6 GB, in the same time.
  Fitted models are bitwise unchanged.
* Classification with `loss = "logit"` fits another 1.2x to 1.8x faster, with bitwise the same fits.
  The terms of the loss under the current weights do not depend on the candidate split, so each observation's terms are now computed once and refreshed only when its weights change, and the terms of classes the observation is not in (which are zero) are skipped: a candidate takes one exponential per logit column and one logarithm per observation of the leaf, instead of two exponentials and two logarithms per column plus two logarithms.
  With 4,000 rows, 300 splits, 200 tree families and `max_interaction = 3` (1 thread), a binary fit takes 74 s instead of 95 s and a three-class fit 91 s instead of 129 s; with 300 rows, 8.8 s → 6.2 s (binary) and 22.5 s → 13.6 s (three classes).
* Growing with `split_structure = "leaf"` (the default) is 2 to 5.5 times as fast, with bitwise the same fits.
  A leaf's rows are put in order by each row's integer rank in the predictor (ranked once per tree family) with a radix sort, instead of a comparison sort of the values, and the leaf keeps that order with the first position of each distinct value, so a candidate split no longer passes over the leaf's values; the thresholds are drawn without allocating memory.
  With the default call at 8 threads, growing takes 38 s instead of 78 s on California housing (16,512 rows) and 6.2 s instead of 34.6 s on 20,000 simulated rows with 10 predictors.
  Classification gains less: on 10,000 simulated rows with 10 predictors growing takes 5.0 s instead of 11.2 s with two classes and 9.2 s instead of 15.4 s with three (`loss = "logit"`; the exponential loss alike).
* Pruning (`prune_components()`, and so the default call) is about 1.7 times as fast on larger data, with bitwise the same models.
  The numbers of training rows in the leaves' boxes, counted in the pass over the tree families, are kept for the later passes (the display's coarsening, its out-of-bag sums, the store of the kept components) instead of being counted again, for trees of up to 4 predictors; and the index of the components' rectangles finds its cuts without sorting.
  On California housing (the default call, 8 threads) pruning takes 23 s instead of 39 s; at `max_interaction = 8` (50 tree families) it takes 240 s instead of 270 s, with no more memory.
* Choosing the number of splits out of bag (`n_splits = NULL`, the default) needs about half the memory while growing, in the same time and with the same fits.
  Each tree family's out-of-bag predictions are added to the checkpoints' sums as the family delivers them, instead of being kept until the stop is decided; the pool of split candidates keeps no record of its entries; and `options(randomPlantedForest.adaptive_pause_gb)`, the memory of the families waiting for the stop, defaults to 0.25 GB instead of 1 GB.
  With the default call at 8 threads, growing's peak physical footprint falls from 2.2 GB to 1.1 GB on California housing (16,512 rows; pruning, at about 2 GB, is now the call's peak there), and from 2.0 GB to 0.8 GB on 20,000 simulated rows with 10 predictors, where the call's peak resident memory (RSS) falls from 2.6 GB to 1.5 GB.
  On 80,000 simulated rows the smaller pause budget alone lowers the call's peak footprint from 2.9 GB to 2.0 GB (RSS 4.9 GB to 4.0 GB), in 41.7 s instead of 42.0 s.
* `rpf_marshal()` of a large forest is several times as fast: the forest's lists of trees and leaves are allocated at their final lengths instead of grown by one element at a time, which copied each list for every element.
  On a California housing forest of 50 tree families and 274,759 leaves, `rpf_marshal()` takes 0.25 s instead of 1.16 s; the blob is the same.

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
