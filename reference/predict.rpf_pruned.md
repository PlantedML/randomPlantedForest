# Predictions of a pruned Random Planted Forest

The prediction of a pruned forest is `alpha + sum(beta_S * f_S(x))` over
the kept components `f_S` of the forest, with the intercept `alpha` and
the coefficients `beta_S` (see
[`prune_components()`](https://plantedml.com/randomPlantedForest/reference/prune_components.md)).
For classification these are the scores, one per score column, which
`type` turns into probabilities or classes as
[`predict.rpf()`](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md)
does.

## Usage

``` r
# S3 method for class 'rpf_pruned'
predict(
  object,
  new_data,
  type = ifelse(object$fit$mode == "regression", "numeric", "prob"),
  nthreads = NULL,
  ...
)
```

## Arguments

- object:

  `[rpf_pruned]`: A pruned forest from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/reference/prune_components.md)
  or
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md).

- new_data:

  `[data.frame | matrix]`: New observations to predict.

- type:

  `[character(1)]`: As in
  [`predict.rpf()`](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md):
  `"numeric"` for regression (the default); for classification `"prob"`
  (the default), `"class"`, or `"numeric"` / `"link"` (the scores
  `alpha + sum(beta_S * f_S(x))`; with `loss = "logit"` or
  `"exponential"` and more than two classes one per class, summing to 0
  across the classes).

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with.

- ...:

  Not currently used, but required for extensibility. Unknown arguments
  are an error.

## Value

A [`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with the
same number of rows as `new_data`: for regression column `.pred`; for
classification the columns of
[`predict.rpf()`](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md)
for `type`.

## See also

[`prune_components()`](https://plantedml.com/randomPlantedForest/reference/prune_components.md),
[`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md).

## Examples

``` r
set.seed(1)
x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
fit <- rpf(x, y, max_interaction = 2, n_families = 10, n_splits = 30, prune = FALSE)
pruned <- prune_components(fit, x, y, rule = "R", tau = 0.01)
predict(pruned, x[1:5, ])
#> # A tibble: 5 × 1
#>     .pred
#>     <dbl>
#> 1 -0.192 
#> 2 -1.73  
#> 3 -0.0109
#> 4  0.610 
#> 5 -1.45  
```
