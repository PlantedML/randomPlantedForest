# Compact printing of forest structures

These methods are provided to avoid flooding the console with long
nested lists containing tree structures. Note

## Usage

``` r
# S3 method for class 'rpf_forest'
print(x, ...)

# S3 method for class 'rpf_forest'
str(object, ...)
```

## Arguments

- x:

  Object of class `rpf_forest`

- ...:

  Further arguments passed to or from other methods.

- object:

  Object of class `rpf_forest`

## See also

[`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md)

## Examples

``` r

rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 10)
print(rpfit$forest)
#> NULL
str(rpfit$forest)
#>  NULL
```
