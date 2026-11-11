# Check whether a model's C++ forest is still alive

The C++ forest does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); an `rpf` object, or
an `rpf_pruned` one fitted with `keep_forest = TRUE`, restored via
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) without
[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)/[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
is unusable, and
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
on it is an error that says so. A pruned model without its forest
(`keep_forest = FALSE`, the default) has no C++ forest and is always
valid; the `fit` of such a model is not (it has no forest to use).

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

`TRUE` if the underlying model can be used, `FALSE` otherwise (also for
the `fit` of a pruned model without its forest).

## See also

[`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
to save a model so that it can be restored.

## Examples

``` r
fit <- rpf(mpg ~ wt + cyl, data = mtcars, n_families = 10, prune = FALSE)
tmp <- tempfile(fileext = ".rds")

# saveRDS() alone loses the C++ forest
saveRDS(fit, tmp)
rpf_is_valid(readRDS(tmp))
#> [1] FALSE

# rpf_marshal() keeps it
saveRDS(rpf_marshal(fit), tmp)
rpf_is_valid(rpf_unmarshal(readRDS(tmp)))
#> [1] TRUE
unlink(tmp)
```
