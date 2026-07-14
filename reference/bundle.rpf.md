# Bundle an rpf model

Method for
[`bundle::bundle()`](https://rstudio.github.io/bundle/reference/bundle.html)
wrapping
[`rpf_marshal()`](http://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)/[`rpf_unmarshal()`](http://plantedml.com/randomPlantedForest/reference/rpf_marshal.md),
so rpf models work with the standard tidymodels serialization workflow.
Training data is not included; see
[`rpf_marshal()`](http://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
for the implications.

## Usage

``` r
# S3 method for class 'rpf'
bundle(x, ...)
```

## Arguments

- x:

  An [rpf](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  object.

- ...:

  Unused.

## Value

An object of class `bundled_rpf` for `bundle()`, or the restored
[rpf](http://plantedml.com/randomPlantedForest/reference/rpf.md) object
for `unbundle()`.

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
#>  1  21.1
#>  2  19.4
#>  3  23.8
#>  4  20.1
#>  5  17.5
#>  6  18.5
#>  7  14.7
#>  8  23.5
#>  9  22.5
#> 10  19.4
#> # ℹ 22 more rows
```
