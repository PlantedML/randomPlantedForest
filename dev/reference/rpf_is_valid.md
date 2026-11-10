# Check whether an rpf object's C++ forest is still alive

The C++ forest does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); an `rpf` object
restored via [`readRDS()`](https://rdrr.io/r/base/readRDS.html) without
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)/[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
is unusable.

## Usage

``` r
rpf_is_valid(x)
```

## Arguments

- x:

  `[rpf | rpf_pruned]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model, or a pruned one from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md),
  possibly restored from disk.

## Value

`TRUE` if the underlying model can be used, `FALSE` otherwise.
