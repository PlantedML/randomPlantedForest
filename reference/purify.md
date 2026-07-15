# Purify a Random Planted Forest

Purifies an rpf object.

## Usage

``` r
purify(x, ...)

# Default S3 method
purify(x, ...)

# S3 method for class 'rpf'
purify(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL)

is_purified(x)
```

## Arguments

- x:

  And object of class `rpf`.

- ...:

  (Unused)

- maxp_interaction:

  integer or NULL: Only compute/store purified components up to this
  interaction order. Higher-order purified trees are zeroed (not
  computed) but still implicitly influence lower orders during
  purification. If NULL, purify all orders (default behavior).

- mode:

  integer(1): Purification algorithm mode. 1 = legacy grid path used by
  `fit$fit$purify()`; 2 = fast exact KD-tree based path. Defaults to 2.

- nthreads:

  integer or NULL: number of threads to use. If NULL, defaults to min of
  the object's configured `nthreads` and available threads.

## Value

Invisibly: The
[`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
object.

## Details

Unless
[`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md) is
called with `purify = TRUE`, the forest has to be purified after fit to
ensure the components extracted by
[`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md)
are valid.
[`predict_components()`](http://plantedml.com/randomPlantedForest/reference/predict_components.md)
will automatically purify a forest if `is_purified()` reports `FALSE`.

## Examples

``` r
rpfit <- rpf(mpg ~., data = mtcars, max_interaction = 2, ntrees = 10)
purify(rpfit)
#> -- Regression Random Planted Forest --
#> 
#> Formula: mpg ~ . 
#> Fit using 10 predictors and 2-degree interactions.
#> Forest is purified!
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
