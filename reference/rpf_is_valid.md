# Check whether an rpf object's C++ forest is still alive

The C++ forest does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); an `rpf` object
restored via [`readRDS()`](https://rdrr.io/r/base/readRDS.html) without
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)/[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
is unusable.

## Usage

``` r
rpf_is_valid(x)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  model, possibly restored from disk.

## Value

`TRUE` if the underlying model can be used, `FALSE` otherwise.
