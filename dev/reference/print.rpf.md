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
rpf(mpg ~ cyl + wt + drat, data = mtcars, max_interaction = 2, n_families = 10, prune = FALSE)
#> ── Regression Random Planted Forest ────────────────────────────────────────────
#> Formula: `mpg ~ cyl + wt + drat`
#> 10 tree families with 6 splits each on 3 predictors, interactions to degree 2.
#> ℹ `n_splits` chosen out of bag (lowest corrected OOB RMSE 2.723); growing stopped at
#> 12 splits.
#> ℹ Forest is not purified.
#> 
#> ── Tree growing 
#>      split_structure: leaf
#>   n_try_split_points: 3
#>     share_candidates: 0.1
#>       max_candidates: 50
#>     split_decay_rate: 0.1
#>        delete_leaves: TRUE
#> 
#> ℹ Fit using 1 thread, also the default for `predict()` and `purify()`.
```
