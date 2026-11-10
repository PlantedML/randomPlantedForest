# The greediness ladder of the out-of-bag tuning

The points of greediness that `rpf(tune = )` searches: pairs of
`n_try_split_points` (split points tried per split candidate) and
`share_candidates` (share of the possible splits drawn as candidates in
each round), greediest first. Each step down the ladder halves
`share_candidates` and `n_try_split_points` (rounded up, at least 1),
from 40 / 1 to 1 / 0.003125, except the first, which takes
`share_candidates` from 1 to 0.8. Row 3 is the other split structures'
default, 10 / 0.4; row 5 is
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)'s
default for `split_structure = "leaf"`, 3 / 0.1.

## Usage

``` r
greediness_ladder()
```

## Value

A `data.frame` with columns `n_try_split_points` and `share_candidates`,
one row per point, greediest first.

## See also

[`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md),
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)

## Examples

``` r
greediness_ladder()
#>    n_try_split_points share_candidates
#> 1                  40         1.000000
#> 2                  20         0.800000
#> 3                  10         0.400000
#> 4                   5         0.200000
#> 5                   3         0.100000
#> 6                   2         0.050000
#> 7                   1         0.025000
#> 8                   1         0.012500
#> 9                   1         0.006250
#> 10                  1         0.003125
```
