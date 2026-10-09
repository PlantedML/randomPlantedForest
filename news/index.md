# Changelog

## randomPlantedForest 0.5.0

### Breaking changes

- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  now defaults to `max_interaction = 2`, main effects and pairwise
  interactions, instead of `1`, an additive model
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
  Values above the number of predictors are capped silently instead of
  with a message.
- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  now defaults to `loss = "exponential"` for classification instead of
  `"L2"`, so classification fits give proper probability estimates by
  default; regression still uses `"L2"`
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  arguments are reordered by purpose (forest size, split search, loss,
  other). Pass arguments after the data arguments by name
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  error on unknown arguments, such as misspelled ones, instead of
  silently ignoring them
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  is a regular function instead of an S3 generic
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `preprocess_predictors_predict()` is no longer exported; it is an
  internal step of [`predict()`](https://rdrr.io/r/stats/predict.html)
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- The package now requires R \>= 4.1.0
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- For a given seed, classification models differ from 0.4.0 due to the
  fixes and speedup below, but regression models are unchanged
  ([\#57](https://github.com/PlantedML/randomPlantedForest/issues/57),
  [\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- For direct users of the C++ object (`$fit`): `set_data()` no longer
  fits, call `fit()` afterwards, and an unknown classification loss is
  an error instead of a silent fallback to L2
  ([@jyliuu](https://github.com/jyliuu),
  [\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).
- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  no longer has a `cv` argument, which never had an effect: the C++
  cross-validation was a no-op
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Bug fixes

- [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  returns one intercept per target level for multiclass classification,
  named by `target_levels`. Previously only the first level’s intercept
  was returned, so components plus intercept did not sum to the
  prediction for other levels.
- [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  error on missing values in `new_data`, naming the affected columns.
  Previously such rows silently received arbitrary predictions, as
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  has no handling for missing values.
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  likewise reports missing values in the training data by column
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Classification fits with a fixed seed could differ between runs. When
  a split replaced a leaf, its new split candidates could be lost or
  attached to the wrong leaf, depending on memory layout
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `predict(type = "class")` breaks probability ties by level order
  instead of at random, so class predictions are deterministic and no
  longer advance the random number generator
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  returns its input invisibly, as documented
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Logical predictors in the formula interface are used as a single 0/1
  column, as in the x/y interface. Previously they were expanded into
  two redundant indicator columns, which also counted against
  `max_interaction`
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Performance

- Prediction is much faster: leaves are scanned from flat,
  cache-friendly arrays (~15x single-threaded), and rows are split
  across threads
  ([\#65](https://github.com/PlantedML/randomPlantedForest/issues/65)).
  [`predict()`](https://rdrr.io/r/stats/predict.html) gains `nthreads`,
  defaulting to the `nthreads` used for fitting.
- [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  is much faster: purified forests look up only the trees of the
  requested component, and each distinct input row is predicted once,
  which pays off for low-cardinality features
  ([\#59](https://github.com/PlantedML/randomPlantedForest/issues/59),
  [\#65](https://github.com/PlantedML/randomPlantedForest/issues/65)).
- Classification fits are slightly faster, as redundant work during
  model construction was removed
  ([\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).

### Other improvements

- New “Getting started” vignette covering regression, classification,
  component decomposition, recipes, and saving models
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- The [`print()`](https://rdrr.io/r/base/print.html) method is
  overhauled and gains a
  [`format()`](https://rdrr.io/r/base/format.html) counterpart: it
  reports the forest size, interaction degree, purification state, split
  search settings and number of threads
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- `deterministic = TRUE` with `ntrees > 1` warns, as all tree families
  are then identical
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).
- Rewritten documentation for
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md),
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md),
  including guidance on choosing parameters and on how purification
  interacts with
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

### Internals

- The C++ core no longer depends on `Rcpp`: it uses standard C++ types
  and exceptions, and an `Rcpp` layer (`src/rcpp_interface.*`) converts
  at the boundary. This is groundwork for bindings in other languages
  such as Python ([@jyliuu](https://github.com/jyliuu),
  [\#57](https://github.com/PlantedML/randomPlantedForest/issues/57)).
- `cli` and `rlang` are now imported, `mvtnorm` is only needed to build
  the `pkgdown` site
  ([\#66](https://github.com/PlantedML/randomPlantedForest/issues/66)).

## randomPlantedForest 0.4.0

It is now possible to serialize and de-serialize an `rpf` object via
marshalling, meaning a fitted model can be saved to disk or dispatched
to a worker process for parallelization or encapsulation, as is done in
`mlr3`. This also re-opens the door for the `mlr3extralearners` wrapper,
which was removed due to the lack of serialization support.

- New
  [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
  /
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
  serialize a fitted forest to a plain R list and back, making
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)-based storage of
  rpf models possible
  ([\#52](https://github.com/PlantedML/randomPlantedForest/issues/52)).
  Purified forests restore their purified state directly; training data
  is only embedded with `include_data = TRUE` (required to
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  after restoring). Blobs record the blob format version and the package
  version they were created with; restoring under an older package than
  the one that saved the model warns. Malformed or corrupt blobs error
  instead of reading out of bounds.
- New
  [`rpf_is_valid()`](https://plantedml.com/randomPlantedForest/reference/rpf_is_valid.md)
  checks whether an rpf object’s internal model is usable,
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  now give an actionable error for `rpf` objects restored via
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) without marshaling.
- New
  [`bundle::bundle()`](https://rstudio.github.io/bundle/reference/bundle.html)
  method for rpf models wrapping the marshaling API.

### Behavior changes

- The default for `delta` changed from 0 to 0.001. With `delta = 0`,
  splits producing single-class nodes have infinite logit loss and are
  always rejected, which also degraded binary logit fits.

### Fixes and improvements

- Fixed multiclass classification with `loss = "logit"`, which
  effectively never split and predicted near-uniform class probabilities
  ([\#40](https://github.com/PlantedML/randomPlantedForest/issues/40)).
  The C++ logit loss is a reference-class multinomial formulation
  expecting `K-1` indicator columns, but the R wrapper passed a full
  `K`-column one-hot matrix, pinning the implicit reference-class
  probability to zero. Multiclass logit outcomes are now encoded with
  the first factor level as reference class:
  - `predict(type = "prob")` reconstructs all `K` class probabilities
    from the `K-1` logits; rows sum to 1 exactly.
  - `predict(type = "numeric"/"link")` now returns `K-1` columns named
    after the non-reference levels (previously `K` columns).
  - [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
    on multiclass logit fits returns per-class components for the `K-1`
    non-reference levels; `target_levels` reflects this.

## randomPlantedForest 0.3.0

### Major changes ([\#61](https://github.com/PlantedML/randomPlantedForest/issues/61))

- New
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  arguments controlling split-candidate sampling:
  - `split_structure = "leaves"`: Defines what a split candidate is and
    how candidates are drawn. One of `"leaves"` (default), `"hist"`,
    `"cur_trees_1"`, `"cur_trees_2"`, or `"res_trees"`; see
    [`?rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
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
- [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  gains arguments:
  - `mode = 2`: Purification algorithm; `2` is a new fast exact method,
    `1` is the legacy grid-based path.
  - `nthreads = NULL`: Purification is now multithreaded, defaulting to
    the fit’s `nthreads` setting.
  - `maxp_interaction = NULL`: Optionally only compute purified
    components up to this interaction order.
- New
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  argument `export_forest = FALSE`: The flattened forest is no longer
  stored in the fitted object by default, so `rpf_object$forest` is
  `NULL` unless `export_forest = TRUE`. This reduces object size;
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md),
  and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  are unaffected.
- `preprocess_predictors_predict()` is now exported.
- Fixed a memory bug in the legacy purification path where the grid was
  sized one element too large, causing out-of-bounds reads (crashes on
  Windows, silently wrong purification results elsewhere).
- Fixed a crash on Windows when fitting with `nthreads > 1`, caused by a
  `thread_local` buffer with a non-trivial destructor being destroyed at
  thread exit.

### Other changes

- Internals in `src/` have been refactored into modular sub-files
  ([\#53](https://github.com/PlantedML/randomPlantedForest/issues/53))
- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  now errors if a regression target is combined with a `loss` other than
  `"L2"`.
- Allow features of type `logical`, which are now converted via
  `as.integer`.
- The `parallel = TRUE|FALSE` argument in
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
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
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md),
  for extended functionality available with
  [`glex`](https://github.com/PlantedML/glex).
- Add `target_levels` vector to output of
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  to aid multiclass handling. Keeping track of levels is somewhat
  awkward since column names of `$m` need to be identifiable regarding
  the target level.

## randomPlantedForest 0.2.0

- Added a `NEWS.md` file to track changes to the package.
