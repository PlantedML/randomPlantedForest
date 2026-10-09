# Random Planted Forest Predictions

Random Planted Forest Predictions

## Usage

``` r
# S3 method for class 'rpf'
predict(
  object,
  new_data,
  type = ifelse(object$mode == "regression", "numeric", "prob"),
  nthreads = NULL,
  ...
)
```

## Arguments

- object:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- new_data:

  `[data.frame | matrix]`: New observations to predict.

- type:

  `[character(1)]`: `"numeric"` for regression outcomes, `"class"` for
  class predictions or `"prob"` for probability predictions. Defaults to
  `"numeric"` for regression and `"prob"` for classification.

  For classification and `loss = "L1"` or `"L2"`, `"numeric"` yields raw
  predictions which are not guaranteed to be valid probabilities in
  `[0, 1]`. For `type = "prob"`, these are truncated to ensure this
  property.

  If `loss` is `"logit"` or `"exponential"`, `type = "link"` is an alias
  for `type = "numeric"`, as in this case the raw predictions have the
  additional interpretation similar to the linear predictor in a
  [`glm`](https://rdrr.io/r/stats/glm.html).

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with, capped at the available cores.

- ...:

  Not currently used, but required for extensibility. Unknown arguments
  are an error.

## Value

For regression: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with column
`.pred` with the same number of rows as `new_data`.

For classification: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with one
column `.pred_<level>` for each level in `y` containing class
probabilities if `type = "prob"`. For `type = "class"`, one column
`.pred_class` with class predictions is returned. For `type = "numeric"`
or `"link"`, raw predictions are returned: one column `.pred` for binary
outcomes, and one column `.pred_<level>` per level for multiclass
outcomes. With `loss = "logit"`, the first level is the reference class
and has no column, so `K - 1` columns are returned for `K` levels.

## Examples

``` r
# Regression with L2 loss
rpfit <- rpf(y = mtcars$mpg, x = mtcars[, c("cyl", "wt")])
predict(rpfit, mtcars[, c("cyl", "wt")])
#> # A tibble: 32 × 1
#>    .pred
#>    <dbl>
#>  1  20.7
#>  2  20.5
#>  3  25.0
#>  4  21.0
#>  5  17.8
#>  6  18.3
#>  7  14.9
#>  8  23.5
#>  9  22.5
#> 10  18.6
#> # ℹ 22 more rows
```
