# Bundle an rpf model

Method for
[`bundle::bundle()`](https://rstudio.github.io/bundle/reference/bundle.html)
wrapping
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)/[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md),
so rpf models work with the standard tidymodels serialization workflow.
Training data is not included; see
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
for the implications.

## Usage

``` r
# S3 method for class 'rpf'
bundle(x, ...)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- ...:

  Unused.

## Value

An object of class `bundled_rpf` for `bundle()`, or the restored
[rpf](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
object for `unbundle()`.

## Examples

``` r
fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10)
b <- bundle::bundle(fit)
tmp <- tempfile(fileext = ".rds")
saveRDS(b, tmp)
restored <- bundle::unbundle(readRDS(tmp))
predict(restored, mtcars)
#> # A tibble: 32 × 1
#>    .pred
#>    <dbl>
#>  1  20.4
#>  2  20.1
#>  3  23.8
#>  4  20.5
#>  5  18.1
#>  6  18.3
#>  7  14.7
#>  8  24.2
#>  9  22.4
#> 10  18.7
#> # ℹ 22 more rows
```
