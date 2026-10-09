# Print an rpf

Print an rpf

## Usage

``` r
# S3 method for class 'rpf'
print(x, ...)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly: `x`.

## See also

[`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).

## Examples

``` r
rpf(mpg ~ cyl + wt + drat, data = mtcars, max_interaction = 2, ntrees = 10)
#> ── Regression Random Planted Forest ────────────────────────────────────────────
#> Formula: `mpg ~ cyl + wt + drat`
#> 10 tree families with 30 splits each on 3 predictors, interactions to degree 2.
#> ℹ Forest is not purified.
#> 
#> ── Tree growing 
#>    split_structure: leaves
#>          split_try: 10
#>              t_try: 0.4
#>     max_candidates: 50
#>   split_decay_rate: 0.1
#>      delete_leaves: TRUE
#> 
#> ℹ Fit using 1 thread, also the default for `predict()` and `purify()`.
```
