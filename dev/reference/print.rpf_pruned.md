# Print a pruned Random Planted Forest

Shows the model: its intercept and its components, each with its
importance, the component's share of the variance of the model's
prediction on the training rows (`components$model_share`). How the
components were chosen (the forest's families and splits, given or
chosen out of bag, the rule and threshold, what was left out, the
heat-map forms) follows briefly at the end.

## Usage

``` r
# S3 method for class 'rpf_pruned'
print(x, n_show = 20, ...)
```

## Arguments

- x:

  `[rpf_pruned]`: A pruned forest from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  or
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).

- n_show:

  `[integer(1): 20]`: The largest number of components to list.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly: `x`.

## See also

[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md).

## Examples

``` r
set.seed(1)
x <- matrix(runif(300 * 3, -1, 1), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(300, sd = 0.3)
fit <- rpf(x, y, max_interaction = 2, n_families = 10, n_splits = 30, prune = FALSE)
prune_components(fit, x, y, rule = "R", tau = 0.01)
#> ── Random Planted Forest model ─────────────────────────────────────────────────
#> Regression. Prediction = intercept + the sum of 3 components.
#> Intercept: -0.0494822
#> 
#> Components (largest first), importance = share of the variance of the prediction on the training rows:
#>  component importance
#>  x2:x3     51%       
#>  x1        37%       
#>  x3        0.093%    
#> 
#> How it was chosen:
#>   The forest: 10 tree families of 30 splits.
#>   From the forest's 6 nonzero components (purified): those with a share of at least 1% of the forest's Var(prediction), each with one coefficient, fitted out of bag (rule R).
#>   Rule and threshold given.
#>   Left out: 3 components, shares together 1.5% (largest x2, 0.75%).
```
