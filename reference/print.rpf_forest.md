# Compact printing of forest structures

The flattened forest that
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
stores with `export_forest = TRUE` (as `$forest`) is a long nested list
of trees. These methods print a one-line summary instead of the whole
list.

## Usage

``` r
# S3 method for class 'rpf_forest'
print(x, ...)

# S3 method for class 'rpf_forest'
str(object, ...)
```

## Arguments

- x:

  `[rpf_forest]`: Flattened forest, as in `$forest` of an
  [`rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md).

- ...:

  Further arguments passed to or from other methods.

- object:

  `[rpf_forest]`: Flattened forest.

## Value

Invisibly: `x` (`object` for
[`str()`](https://rdrr.io/r/utils/str.html)).

## See also

[`rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md)

## Examples

``` r

rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, n_families = 10, export_forest = TRUE, prune = FALSE)
print(rpfit$forest)
#> <rpf_forest> of 10 trees
str(rpfit$forest)
#> <rpf_forest> of 10 trees
```
