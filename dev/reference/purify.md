# Purify a Random Planted Forest

Purification makes the components of the forest's functional
decomposition unique, which
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
relies on. Unless
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
was called with `purify = TRUE`,
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
purifies the forest on first use.

## Usage

``` r
purify(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL)

is_purified(x)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- ...:

  Reserved for future expansion.

- maxp_interaction:

  `[integer(1) | NULL: NULL]`: Highest interaction order to purify.
  Higher-order components are set to zero, but still influence lower
  orders during purification, and
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  then returns zero for them. `NULL` purifies all orders.

- mode:

  `[integer(1): 2]`: Purification algorithm: `2` is the fast exact
  KD-tree based algorithm, `1` the original grid-based one.

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with, capped at the available cores.

## Value

`purify()` returns `x` invisibly. `is_purified()` returns `TRUE` or
`FALSE`.

## Details

The forest is modified in place: `x` and every copy of it are purified,
whether or not the result is assigned. `purify()` is idempotent, meaning
if the forest is already purified it just returns it unmodified.

## Examples

``` r
rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, ntrees = 10)
is_purified(rpfit)
#> [1] FALSE
purify(rpfit)
is_purified(rpfit)
#> [1] TRUE
```
