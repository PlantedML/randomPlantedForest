# Serialize and restore Random Planted Forests

`rpf_marshal()` converts an
[rpf](http://plantedml.com/randomPlantedForest/reference/rpf.md) object
into a plain R list safe for
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); `rpf_unmarshal()`
reverses it. The C++ forest behind an `rpf` object does not survive R
serialization, so this explicit round-trip is required to store or
transfer fitted models.

## Usage

``` r
rpf_marshal(x, include_data = FALSE)

rpf_unmarshal(blob)
```

## Arguments

- x:

  An object of class `rpf`.

- include_data:

  `[FALSE]`: Store training data in the blob, enabling
  [`purify()`](http://plantedml.com/randomPlantedForest/reference/purify.md)
  after restoring.

- blob:

  An object of class `rpf_marshaled` created by `rpf_marshal()`.

## Value

`rpf_marshal()` returns a list of class `rpf_marshaled`;
`rpf_unmarshal()` returns an object of class
[rpf](http://plantedml.com/randomPlantedForest/reference/rpf.md).

## Details

Training data is only included with `include_data = TRUE`; without it, a
restored forest can predict (including purified prediction if the forest
was purified before marshaling) but can never be purified afterwards.

## Examples

``` r
fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10)
blob <- rpf_marshal(fit)
tmp <- tempfile(fileext = ".rds")
saveRDS(blob, tmp)
restored <- rpf_unmarshal(readRDS(tmp))
all.equal(predict(fit, mtcars), predict(restored, mtcars))
#> [1] TRUE
```
