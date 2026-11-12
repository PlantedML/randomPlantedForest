# Settings of the out-of-bag tuning

What `rpf(tune = )` searches. The search is coordinate-wise, so that it
needs few forests (9 to 11 with the defaults):

1.  the points `search` of `ladder` at the first `min_leaf_rows`;

2.  every `min_leaf_rows` at the point chosen in step 1 (skipped if only
    one is given);

3.  the two ladder neighbours of that point at the `min_leaf_rows`
    chosen in step 2 (if `neighbours`).

## Usage

``` r
tune_control(
  ladder = greediness_ladder(),
  search = 3:8,
  min_leaf_rows = c(1, 5, 20, 100),
  neighbours = TRUE
)
```

## Arguments

- ladder:

  `[data.frame: greediness_ladder()]`: The points of greediness,
  greediest first: columns `n_try_split_points` (whole numbers, at
  least 1) and `share_candidates` (in `(0, 1]`). The column names of
  0.5.0, `split_try` and `t_try`, still work, with a deprecation
  warning.

- search:

  `[integer: 3:8]`: The rows of `ladder` searched in step 1 (on the
  default ladder from 10 / 0.4 to 1 / 0.0125).

- min_leaf_rows:

  `[integer: c(1, 5, 20, 100)]`: The values of
  [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)'s
  `min_leaf_rows` searched in step 2; the first one is that of step 1.

- neighbours:

  `[logical(1): TRUE]`: Whether to run step 3.

## Value

An object of class `rpf_tune_control`, a list with the four settings,
for `rpf(tune = )`.

## Details

Every forest has its number of splits chosen out of bag while it grows
(as `rpf(n_splits = NULL)`) and is pruned by
[`prune_components()`](https://plantedml.com/randomPlantedForest/reference/prune_components.md)
with its defaults, without the display. A step chooses, among all
forests grown at the coordinate it varies, the one whose pruned model
has the lowest corrected out-of-bag RMSE (ties: the greedier point, the
larger `min_leaf_rows`). Step 3 chooses among every forest grown at the
chosen `min_leaf_rows`, so with `min_leaf_rows` 1 it also compares step
1's.

`rpf(tune = TRUE)` is `tune_control()`, and `rpf(tune = "quick")` is
`tune_control(search = c(5, 7), min_leaf_rows = 1, neighbours = FALSE)`:
two forests, 3 / 0.1 and 1 / 0.025.

## See also

[`greediness_ladder()`](https://plantedml.com/randomPlantedForest/reference/greediness_ladder.md),
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)

## Examples

``` r
# The full search and the quick preset
tune_control()
#> $ladder
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
#> 
#> $search
#> [1] 3 4 5 6 7 8
#> 
#> $min_leaf_rows
#> [1]   1   5  20 100
#> 
#> $neighbours
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "rpf_tune_control"
tune_control(search = c(5, 7), min_leaf_rows = 1, neighbours = FALSE)
#> $ladder
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
#> 
#> $search
#> [1] 5 7
#> 
#> $min_leaf_rows
#> [1] 1
#> 
#> $neighbours
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "rpf_tune_control"

# \donttest{
# A small search: two points of greediness, then two leaf sizes
set.seed(1)
x <- matrix(runif(400 * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(400, sd = 0.3)
tuned <- rpf(
  x,
  y,
  max_interaction = 2,
  n_families = 10,
  tune = tune_control(search = c(3, 5), min_leaf_rows = c(1, 20), neighbours = FALSE)
)
tuned$tune$search
#>   step n_try_split_points share_candidates min_leaf_rows n_splits stopped_at
#> 1    1                 10              0.4             1       84        100
#> 2    1                  3              0.1             1       59         71
#> 3    2                 10              0.4            20      168        200
#>   capped oob_rmse_forest rule tau n_kept  oob_rmse seconds chosen_step1
#> 1  FALSE       0.4844060    R 0.1      2 0.3666121   0.059         TRUE
#> 2  FALSE       0.4459042    R 0.1      2 0.3752002   0.034        FALSE
#> 3  FALSE       0.3832057    R 0.1      2 0.3539268   0.130           NA
#>   chosen_step2 chosen_step3 final
#> 1        FALSE           NA FALSE
#> 2           NA           NA FALSE
#> 3         TRUE           NA  TRUE
# }
```
