# Checkpoints of the out-of-bag path

The numbers of splits at which
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
with `oob_path = TRUE` evaluates the out-of-bag (OOB) error: a geometric
grid with four points per doubling, `round(25 * 2^(j / 4))` from 6 up to
`n_splits` (it contains 12, 25, 50, ..., 1600, 3200), and `n_splits`
itself; every number up to 12 if `n_splits <= 12`. The error changes
slowly in `log(n_splits)`, so a geometric grid resolves it evenly, and
the path's memory, `2 * nrow(x) * length(checkpoints)` doubles, stays
small; every split would cost `2 * nrow(x) * n_splits`.

## Usage

``` r
oob_path_checkpoints(n_splits, ...)
```

## Arguments

- n_splits:

  `[integer(1)]`: The largest number of splits (the `n_splits` of the
  fit). Was `splits`.

- ...:

  Only the old name `splits`, deprecated: it works, with a warning.

## Value

An increasing integer vector.

## See also

[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
(`oob_path`)

## Examples

``` r
oob_path_checkpoints(400)
#>  [1]   6   7   9  11  12  15  18  21  25  30  35  42  50  59  71  84 100 119 141
#> [20] 168 200 238 283 336 400
```
