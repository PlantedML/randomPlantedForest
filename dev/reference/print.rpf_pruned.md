# Print a pruned Random Planted Forest

Shows the rule and threshold, how they were chosen, what was kept and
dropped, and the coefficients.

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

  `[integer(1): 20]`: The largest number of kept components to list.

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
fit <- rpf(x, y, max_interaction = 2, ntrees = 10, splits = 30)
prune_components(fit, x, y, rule = "R", tau = 0.01)
#> ── Pruned Random Planted Forest ────────────────────────────────────────────────
#> Rule R (drop + refit one coefficient per component, fitted out of bag), tau = 1% of Var(prediction)
#> Given, not chosen.
#> 
#> Nonzero components: 6
#> Kept:    6 (by order 1: 3, 2: 3)
#> Dropped: 0
#> 
#> Display: every kept component in at most 10 heat maps (nothing to coarsen, so not pruned again).
#> First pruning: rule R, tau = 1%, 6 components kept.
#> Heat maps of the model: 6
#> 
#> Intercept: -0.0251072
#> Kept components (largest first), share on the training rows, coefficient and heat maps:
#>  component share beta     heat_maps
#>  x1        57%    0.90346 1        
#>  x2:x3     12%    1.83593 1        
#>  x1:x2     2.8%  -0.06836 1        
#>  x3        2.6%   0.41588 1        
#>  x2        2.4%   0.25804 1        
#>  x1:x3     1.4%  -0.35157 1        
```
