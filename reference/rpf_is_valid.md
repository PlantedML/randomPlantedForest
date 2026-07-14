# Check whether an rpf object's C++ forest is still alive

The C++ forest does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); an `rpf` object
restored via [`readRDS()`](https://rdrr.io/r/base/readRDS.html) without
[`rpf_marshal()`](http://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)/[`rpf_unmarshal()`](http://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
is unusable.

## Usage

``` r
rpf_is_valid(x)
```

## Arguments

- x:

  An object of class `rpf`.

## Value

`TRUE` if the underlying model can be used, `FALSE` otherwise.
