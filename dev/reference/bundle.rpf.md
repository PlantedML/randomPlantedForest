# Bundle an rpf model

Method for
[`bundle::bundle()`](https://rstudio.github.io/bundle/reference/bundle.html)
wrapping
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)/[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md),
so rpf models work with the standard tidymodels serialization workflow.
Training data is not included; see
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
for the implications. Pruned forests from
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
or
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
are bundled the same way.

## Usage

``` r
# S3 method for class 'rpf'
bundle(x, ...)

# S3 method for class 'rpf_pruned'
bundle(x, ...)
```

## Arguments

- x:

  `[rpf | rpf_pruned]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model, or a pruned one from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md).

- ...:

  Unused.

## Value

An object of class `bundled_rpf` (`bundled_rpf_pruned`) for `bundle()`,
or the restored
[rpf](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
(`rpf_pruned`) object for `unbundle()`.

## Examples

``` r
fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10, prune = FALSE)
b <- bundle::bundle(fit)
tmp <- tempfile(fileext = ".rds")
saveRDS(b, tmp)
restored <- bundle::unbundle(readRDS(tmp))
predict(restored, mtcars)
#> # A tibble: 32 × 1
#>    .pred
#>    <dbl>
#>  1  20.7
#>  2  20.6
#>  3  23.3
#>  4  20.2
#>  5  16.5
#>  6  18.2
#>  7  15.8
#>  8  23.7
#>  9  22.7
#> 10  18.6
#> # ℹ 22 more rows
```
