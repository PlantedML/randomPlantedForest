# The greediness ladder of the out-of-bag tuning

The points of greediness that `rpf(tune = )` searches: pairs of
`split_try` (thresholds evaluated per split candidate) and `t_try`
(share of the possible splits drawn as candidates in each round),
greediest first. Each step down the ladder halves `t_try` and
`split_try` (rounded up, at least 1), from 40 / 1 to 1 / 0.003125,
except the first, which takes `t_try` from 1 to 0.8. Row 3 is
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)'s
default, 10 / 0.4.

## Usage

``` r
greediness_ladder()
```

## Value

A `data.frame` with columns `split_try` and `t_try`, one row per point,
greediest first.

## See also

[`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md),
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)

## Examples

``` r
greediness_ladder()
#>    split_try    t_try
#> 1         40 1.000000
#> 2         20 0.800000
#> 3         10 0.400000
#> 4          5 0.200000
#> 5          3 0.100000
#> 6          2 0.050000
#> 7          1 0.025000
#> 8          1 0.012500
#> 9          1 0.006250
#> 10         1 0.003125
```
