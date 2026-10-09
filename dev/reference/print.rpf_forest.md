# Compact printing of forest structures

These methods are provided to avoid flooding the console with long
nested lists containing tree structures.

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
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).

- ...:

  Further arguments passed to or from other methods.

- object:

  `[rpf_forest]`: Flattened forest.

## See also

[`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)

## Examples

``` r

rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 10)
print(rpfit$forest)
#> NULL
str(rpfit$forest)
#>  NULL
```
