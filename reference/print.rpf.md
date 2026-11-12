# Print a Random Planted Forest

Shows a forest (class `rpf`, from
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
with `prune = FALSE`): the kind of model and its formula, the number of
tree families and of splits (with how the splits were chosen), the
predictors and the interaction order, whether the forest is purified,
the settings of the split search and the number of threads. A pruned
model, the default of
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md),
is printed by
[`print.rpf_pruned()`](https://plantedml.com/randomPlantedForest/reference/print.rpf_pruned.md).

## Usage

``` r
# S3 method for class 'rpf'
print(x, ...)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  model.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly: `x`.

## See also

[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md),
[`print.rpf_pruned()`](https://plantedml.com/randomPlantedForest/reference/print.rpf_pruned.md).

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
