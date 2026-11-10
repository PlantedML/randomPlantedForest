# Predictions of a pruned Random Planted Forest

The prediction of a pruned forest is `alpha + sum(beta_S * f_S(x))` over
the kept components `f_S` of the forest, with the intercept `alpha` and
the coefficients `beta_S` (see
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)).

## Usage

``` r
# S3 method for class 'rpf_pruned'
predict(object, new_data, nthreads = NULL, ...)
```

## Arguments

- object:

  `[rpf_pruned]`: A pruned forest from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  or
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).

- new_data:

  `[data.frame | matrix]`: New observations to predict.

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with.

- ...:

  Not currently used, but required for extensibility. Unknown arguments
  are an error.

## Value

A [`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with
column `.pred` with the same number of rows as `new_data`.

## See also

[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md).

## Examples

``` r
set.seed(1)
x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
fit <- rpf(x, y, max_interaction = 2, ntrees = 10, splits = 30)
pruned <- prune_components(fit, x, y, rule = "R", tau = 0.01)
predict(pruned, x[1:5, ])
#> # A tibble: 5 × 1
#>    .pred
#>    <dbl>
#> 1 -0.121
#> 2 -1.46 
#> 3  0.156
#> 4  0.528
#> 5 -1.18 
```
