# Random Planted Forest Predictions

Random Planted Forest Predictions

## Usage

``` r
# S3 method for class 'rpf'
predict(
  object,
  new_data,
  type = ifelse(object$mode == "regression", "numeric", "prob"),
  ...
)
```

## Arguments

- object:

  A fit object of class
  [`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md).

- new_data:

  Data for new observations to predict.

- type:

  `"numeric"` for regression outcomes, `"class"` for class predictions
  or `"prob"` for probability predictions.

  For classification and `loss = "L1"` or `"L2"`, `"numeric"` yields raw
  predictions which are not guaranteed to be valid probabilities in
  `[0, 1]`. For `type = "prob"`, these are truncated to ensure this
  property.

  If `loss` is `"logit"` or `"exponential"`, `type = "link"` is an alias
  for `type = "numeric"`, as in this case the raw predictions have the
  additional interpretation similar to the linear predictor in a
  [`glm`](https://rdrr.io/r/stats/glm.html).

- ...:

  Unused.

## Value

For regression: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with column
`.pred` with the same number of rows as `new_data`.

For classification: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with one
column for each level in `y` containing class probabilities if
`type = "prob"`. For `type = "class"`, one column `.pred` with class
predictions is returned. For `type = "numeric"` or `"link"`, one column
`.pred` with raw predictions.

## Examples

``` r
# Regression with L2 loss
rpfit <- rpf(y = mtcars$mpg, x = mtcars[, c("cyl", "wt")])
predict(rpfit, mtcars[, c("cyl", "wt")])
#> # A tibble: 32 × 1
#>    .pred
#>    <dbl>
#>  1  20.7
#>  2  20.3
#>  3  25.2
#>  4  21.2
#>  5  17.1
#>  6  18.6
#>  7  15.1
#>  8  24.1
#>  9  23.0
#> 10  19.4
#> # ℹ 22 more rows
```
