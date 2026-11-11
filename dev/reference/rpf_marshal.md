# Serialize and restore Random Planted Forests

`rpf_marshal()` converts an
[rpf](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
object into a plain R list safe for
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

  `[rpf | rpf_pruned]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model, or a pruned one from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md).

- include_data:

  `[logical(1): FALSE]`: Store training data in the blob, enabling
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  after restoring.

- blob:

  `[rpf_marshaled | rpf_pruned_marshaled]`: Object created by
  `rpf_marshal()`.

## Value

`rpf_marshal()` returns a list of class `rpf_marshaled`
(`rpf_pruned_marshaled` for a pruned forest); `rpf_unmarshal()` returns
an object of class
[rpf](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
(`rpf_pruned`).

## Details

Training data is only included with `include_data = TRUE`; without it, a
restored forest can predict (including purified prediction if the forest
was purified before marshaling) but can never be purified afterwards.
With it, the in-bag counts of the bootstrap samples are stored too, so
that
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
can choose out of bag after restoring.

A pruned forest from
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
or
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
(class `rpf_pruned`) is marshaled the same way: its forest as above, and
the pruning (rule, coefficients, components, the table of the choice,
the stored kept components, the display's coarsened components) as it
is. The restored model predicts as before, from its stored components,
so it needs no training data for that. A pruned model without its forest
(the default of
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
and
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md),
`keep_forest = FALSE`) is plain R data:
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) stores it as it is,
and its blob is the model without the caches of its built indices
(`include_data` has nothing to include).

## See also

[`rpf_is_valid()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_is_valid.md);
[`bundle()`](https://plantedml.com/randomPlantedForest/dev/reference/bundle.rpf.md)
for the same through the bundle package.

## Examples

``` r
fit <- rpf(mpg ~ wt + cyl, data = mtcars, n_families = 10, prune = FALSE)
blob <- rpf_marshal(fit)
tmp <- tempfile(fileext = ".rds")
saveRDS(blob, tmp)
restored <- rpf_unmarshal(readRDS(tmp))
all.equal(predict(fit, mtcars), predict(restored, mtcars))
#> [1] TRUE
unlink(tmp)
```
