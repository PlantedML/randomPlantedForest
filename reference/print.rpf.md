# Print an rpf fit

Print an rpf fit

## Usage

``` r
# S3 method for class 'rpf'
print(x, ...)
```

## Arguments

- x:

  And object of class `rpf`.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly: `x`.

## See also

[`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md).

## Examples

``` r
rpf(mpg ~ cyl + wt + drat, data = mtcars, max_interaction = 2, ntrees = 10)
#> -- Regression Random Planted Forest --
#> 
#> Formula: mpg ~ cyl + wt + drat 
#> Fit using 3 predictors and 2-degree interactions.
#> Forest is _not_ purified!
#> 
#> Called with parameters:
#> 
#>              loss: L2
#>            ntrees: 10
#>   max_interaction: 2
#>            splits: 30
#>         split_try: 10
#>             t_try: 0.4
#>  split_decay_rate: 0.1
#>    max_candidates: 50
#>     delete_leaves: TRUE
#>   split_structure: leaves
#>             delta: 0.001
#>           epsilon: 0.1
#>     deterministic: FALSE
#>          nthreads: 1
#>            purify: FALSE
#>                cv: FALSE
```
