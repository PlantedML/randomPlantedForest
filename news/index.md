# Changelog

## randomPlantedForest 0.3.0

### Major changes ([\#61](https://github.com/PlantedML/randomPlantedForest/issues/61))

- New
  [`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  arguments controlling split-candidate sampling:
  - `split_structure = "leaves"`: Defines what a split candidate is and
    how candidates are drawn. One of `"leaves"` (default), `"hist"`,
    `"cur_trees_1"`, `"cur_trees_2"`, or `"res_trees"`; see
    [`?rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
    for details.
  - `max_candidates = 50`: Maximum number of split candidates sampled
    per iteration.
  - `split_decay_rate = 0.1`: Exponential aging of repeatedly drawn but
    unchosen split candidates. `split_decay_rate = 0` corresponds to no
    aging and uniform sampling.
  - `delete_leaves = TRUE`: Whether a parent leaf is deleted when
    splitting along an existing dimension.
- **Fitting results change**: The new candidate-sampling defaults and a
  reworked internal RNG mean that fits are not reproducible against
  previous versions, even with the same seed. Install an older commit if
  exact reproduction of previous results is required.
- Seeded fits are now reproducible regardless of `nthreads`: per-tree
  seeds are drawn from R’s RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) gives identical
  forests for serial and multithreaded fits.
- Substantial speedups in fitting (cached per-leaf orderings, prefix
  sums) and reduced memory use (training-only buffers are released after
  each tree family is built).
- [`purify()`](http://plantedml.com/randomPlantedForest/reference/purify.md)
  gains arguments:
  - `mode = 2`: Purification algorithm; `2` is a new fast exact method,
    `1` is the legacy grid-based path.
  - `nthreads = NULL`: Purification is now multithreaded, defaulting to
    the fit’s `nthreads` setting.
  - `maxp_interaction = NULL`: Optionally only compute purified
    components up to this interaction order.
- New
  [`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  argument `export_forest = FALSE`: The flattened forest is no longer
  stored in the fitted object by default, so `rpf_object$forest` is
  `NULL` unless `export_forest = TRUE`. This reduces object size;
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`purify()`](http://plantedml.com/randomPlantedForest/reference/purify.md),
  and
  [`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md)
  are unaffected.
- [`preprocess_predictors_predict()`](http://plantedml.com/randomPlantedForest/reference/preprocess_predictors_predict.md)
  is now exported.
- Fixed a memory bug in the legacy purification path where the grid was
  sized one element too large, causing out-of-bounds reads (crashes on
  Windows, silently wrong purification results elsewhere).
- Fixed a crash on Windows when fitting with `nthreads > 1`, caused by a
  `thread_local` buffer with a non-trivial destructor being destroyed at
  thread exit.

### Other changes

- Internals in `src/` have been refactored into modular sub-files
  ([\#53](https://github.com/PlantedML/randomPlantedForest/issues/53))
- [`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  now errors if a regression target is combined with a `loss` other than
  `"L2"`.
- Allow features of type `logical`, which are now converted via
  `as.integer`.
- The `parallel = TRUE|FALSE` argument in
  [`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  has been substituted by an `nthreads = 1L` argument, allowing for more
  flexible parallelization. The previous behavior only allowed for
  either no parallelization or using n-1 of n available cores. The new
  implementation should be reasonably robust and the default behavior
  remains serial execution.
- Remove `SystemRequirements` field from `DESCRIPTION`: Now the default
  C++ version is C++17 and with a minor change to internal use of random
  numbers, `randomPlantedForest` is now compatible with C++11 through
  C++23.
- Add `remainder` term to `predict_components` output for case where
  `max_interaction` supplied is smaller than `max_interaction` in `rpf`
  fit. In that case, the `m` values don’t sum up to the global
  predictions, so we add a remainder to allow reconstruction of that
  property.

## randomPlantedForest 0.2.1

- Add `glex` class to output of
  [`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md),
  for extended functionality available with
  [`glex`](https://github.com/PlantedML/glex).
- Add `target_levels` vector to output of
  [`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md)
  to aid multiclass handling. Keeping track of levels is somewhat
  awkward since column names of `$m` need to be identifiable regarding
  the target level.

## randomPlantedForest 0.2.0

- Added a `NEWS.md` file to track changes to the package.
