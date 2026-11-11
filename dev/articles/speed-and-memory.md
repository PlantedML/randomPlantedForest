# Speed, memory and saving models

``` r

library(randomPlantedForest)

make_data <- function(n) {
  d <- data.frame(matrix(runif(n * 5, -1, 1), ncol = 5, dimnames = list(NULL, paste0("x", 1:5))))
  d$y <- sin(pi * d$x1) + d$x2 + 2 * d$x1 * d$x3 + 3 * d$x2 * d$x3 * d$x4 + rnorm(n, sd = 0.3)
  d
}
```

## Threads

`nthreads` sets the number of threads for fitting, and is the default
for [`predict()`](https://rdrr.io/r/stats/predict.html) and
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
of the model. The default is 1. Threads never change the result: with
the same seed, the forest is the same at any number of threads.

``` r

set.seed(1)
train <- make_data(1000)

time_fit <- function(threads) {
  set.seed(2)
  t <- system.time(m <- rpf(y ~ ., data = train, nthreads = threads))[["elapsed"]]
  list(model = m, seconds = t)
}
one <- time_fit(1)
two <- time_fit(2)
c(one_thread = one$seconds, two_threads = two$seconds)
#>  one_thread two_threads 
#>       9.028       5.095
identical(predict(one$model, train), predict(two$model, train))
#> [1] TRUE
```

The times are those of the machine that built this page.

## What the time grows with

The default call does three things: it grows the forest until the
out-of-bag error stops falling, purifies and prunes it, and coarsens its
three-way components. Its time grows with:

- **the number of rows**, directly and through the number of splits:
  more data support more splits before the out-of-bag error stops
  falling;
- **the number of tree families** (`n_families`, 200): linearly, in
  fitting and in prediction. Fewer families are fine for a first look,
  and more are usually more accurate;
- **the order** (`max_interaction`, 3): a family can grow into many more
  trees at a higher order, and purification and pruning work on many
  more components. No limit (`0`) is often much slower than order 3, and
  usually not more accurate;
- **the number of predictors**: the possible components number
  `sum(choose(p, 1:3))` at order 3, 175 for 10 predictors and 1,350 for
  20.

The number of splits chosen, and the time, at a few sizes of the data
above:

``` r

sizes <- sapply(c(250, 500, 1000), function(n) {
  set.seed(1)
  d <- make_data(n)
  set.seed(2)
  t <- system.time(m <- rpf(y ~ ., data = d, nthreads = 2))[["elapsed"]]
  c(rows = n, splits = m$fit$splits_oob$n_splits, seconds = round(t, 1))
})
t(sizes)
#>      rows splits seconds
#> [1,]  250    200     1.9
#> [2,]  500    238     2.8
#> [3,] 1000    336     5.1
```

Ways to make a fit faster, at some cost in accuracy or in what the model
can show:

- more threads (`nthreads`);
- fewer families (`n_families`), or a lower order (`max_interaction`);
- `split_structure = "leaf_binned"`, which tries split points from
  quantile bins computed once per predictor;
- `n_splits = <number>` and `prune = FALSE`, which skip the choice and
  the pruning (the model is then the whole forest).

## Memory

While the number of splits is chosen, families grow one at a time,
`nthreads` together, and keep only their out-of-bag predictions and a
short record of their splits; the chosen forest is rebuilt from these
records. A family that reaches the estimated stop can wait there with
its growing state, so that it can grow on cheaply if the stop comes
later; `options(randomPlantedForest.adaptive_pause_gb = 1)` bounds the
memory of the waiting families (in GB). More memory saves time; less is
freed and regrown when needed. The result does not depend on it.

Pruning purifies one family at a time and keeps only the pieces of the
components it keeps, so it needs far less memory than
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
of the whole forest. A pruned model stores its kept components, not the
forest, whose memory is freed when
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
returns. `keep_forest = TRUE` keeps the forest (unpurified) in
`model$fit`, to prune it again or to use it.

## Saving and restoring models

A pruned model (the default) is plain R data:
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) keep it as it is.

``` r

model <- two$model
tmp <- tempfile(fileext = ".rds")
saveRDS(model, tmp)
all.equal(predict(readRDS(tmp), train), predict(model, train))
#> [1] TRUE
```

A forest (`prune = FALSE`), or a pruned model fitted with
`keep_forest = TRUE`, holds a C++ object, which
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) does not preserve:

``` r

set.seed(3)
forest <- rpf(y ~ ., data = train, prune = FALSE, n_families = 20, nthreads = 2)
saveRDS(forest, tmp)
broken <- readRDS(tmp)
rpf_is_valid(broken)
#> [1] FALSE
predict(broken, head(train))
#> Error in `predict()`:
#> ! The C++ forest behind this <rpf> object is gone.
#> ℹ Most likely it was saved with `saveRDS()` and restored with `readRDS()`.
#> ℹ Use `blob <- rpf_marshal(x)` before saving and `rpf_unmarshal(blob)` after
#>   loading, see `?rpf_marshal()`.
```

[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
turns it into a plain R list before saving, and
[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
restores it:

``` r

saveRDS(rpf_marshal(forest), tmp)
restored <- rpf_unmarshal(readRDS(tmp))
all.equal(predict(restored, train), predict(forest, train))
#> [1] TRUE
```

[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
works for pruned models too, with or without their forest, so code that
saves models of both kinds can always use it. A restored model predicts
without its training data. `rpf_marshal(include_data = TRUE)` also
stores the training data and the bootstrap samples, which a restored
forest needs to be purified or pruned.

The [bundle](https://rstudio.github.io/bundle/) package does the same
through its generic interface, which tidymodels uses:

``` r

b <- bundle::bundle(forest)
saveRDS(b, tmp)
restored <- bundle::unbundle(readRDS(tmp))
all.equal(predict(restored, train), predict(forest, train))
#> [1] TRUE
```

## tidymodels

[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
follows the conventions of tidymodels’
[hardhat](https://hardhat.tidymodels.org/): it takes a formula, x and y,
or a [recipe](https://recipes.tidymodels.org/), and
[`predict()`](https://rdrr.io/r/stats/predict.html) returns a tibble
with `.pred` (regression), `.pred_<level>` (probabilities) or
`.pred_class` columns, one row per row of `new_data`.

A recipe bundles preprocessing with the model, for example imputing
missing values, which
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
itself does not accept:

``` r

library(recipes)

# rows with a missing outcome are dropped by hand
peng <- subset(palmerpenguins::penguins, !is.na(body_mass_g))

rec <- recipe(body_mass_g ~ ., data = peng) |>
  step_impute_median(all_numeric_predictors()) |>
  step_impute_mode(all_nominal_predictors())

set.seed(1)
by_recipe <- rpf(rec, data = peng, nthreads = 2)
predict(by_recipe, head(peng))
#> # A tibble: 6 × 1
#>   .pred
#>   <dbl>
#> 1 3746.
#> 2 3444.
#> 3 3568.
#> 4 3507.
#> 5 4017.
#> 6 3400.
```

The recipe’s steps are stored with the model and applied in
[`predict()`](https://rdrr.io/r/stats/predict.html).

## mlr3

[mlr3extralearners](https://mlr3extralearners.mlr-org.com/) has the
learners `regr.rpf` and `classif.rpf`.
