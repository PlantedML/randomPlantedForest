# Extract predicted components from a Random Planted Forest

Prediction components are a functional decomposition of the model
prediction. The sum of all components equals the overall predicted value
for an observation.

## Usage

``` r
predict_components(object, new_data, max_interaction = NULL, predictors = NULL)
```

## Arguments

- object:

  A fit object of class
  [`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md).

- new_data:

  Data for new observations to predict.

- max_interaction:

  [`integer`](https://rdrr.io/r/base/integer.html) or `NULL`: Maximum
  degree of interactions to consider. Default will use the
  `max_interaction` parameter from the
  [`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  object. Must be between `1` (main effects only) and the
  `max_interaction` of the
  [`rpf`](http://plantedml.com/randomPlantedForest/reference/rpf.md)
  object.

- predictors:

  [`character`](https://rdrr.io/r/base/character.html) or `NULL`: Vector
  of one or more column names of predictor variables in `new_data` to
  extract components for. If `NULL`, all variables and their
  interactions are returned.

## Value

A `list` with elements:

- `m`
  ([`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)):
  Components for each main effect and interaction term, representing the
  functional decomposition of the prediction. All components together
  with the intercept sum up to the prediction. For multiclass
  classification, the number of output columns is multiplied by the
  number of levels in the outcome.

- `intercept` (`numeric(1)`): Expected value of the prediction.

- `x`
  ([`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)):
  Copy of `new_data` containing predictors selected by `predictors`.

- `target_levels` (`character`): For multiclass classification only:
  Vector of target levels which can be used to disassemble `m`, as names
  include both term and target level.

## Details

Extracts all possible components up to `max_interaction` degrees, up to
the value set when calling
[`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md).
The intercept is always included. Optionally `predictors` can be
specified to only include components including the given variables. If
`max_interaction` is greater than `length(predictors)`, the
`max_interaction` will be lowered accordingly.

## Note

Depending on the number of predictors and `max_interaction`, the number
of components will increase drastically to
`sum(choose(ncol(new_data), seq_len(max_interaction)))`.

## Examples

``` r

# Regression task, only some predictors
train <-  mtcars[1:20, 1:4]
test <-  mtcars[21:32, 1:4]

set.seed(23)
rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)

# Extract all components, including main effects and interaction terms up to `max_interaction`
(components <- predict_components(rpfit, test))
#> $m
#>            cyl       disp         hp    cyl:disp       cyl:hp     disp:hp
#>          <num>      <num>      <num>       <num>        <num>       <num>
#>  1:  0.9013797  3.2511514  1.3776657  0.08511109 -0.002332250 -0.11682939
#>  2: -0.6706529 -0.9186116 -0.4097990 -0.14816938 -0.079604722  0.10864962
#>  3: -0.6706529 -0.9186116 -0.4097990 -0.14816938 -0.079604722  0.10864962
#>  4: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#>  5: -0.6706529 -1.1532004 -0.3559023 -0.07494483 -0.084604722  0.12038493
#>  6:  0.9013797  4.6372466  4.7893287  0.17994554 -0.091148415  0.36694210
#>  7:  0.9013797  3.2511514  4.7893287  0.08511109 -0.091148415  0.07237951
#>  8:  0.9013797  4.6372466  0.6132023  0.17994554  0.006048960 -0.46138786
#>  9: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#> 10: -0.1689229  0.9540856 -0.3559023 -0.08618539 -0.045141470 -0.09680110
#> 11: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#> 12:  0.9013797  3.2511514  0.9815403  0.08511109 -0.005884567  0.03742733
#>     cyl:disp:hp
#>           <num>
#>  1: -0.06166953
#>  2:  0.04480135
#>  3:  0.04480135
#>  4:  0.07122433
#>  5:  0.06090002
#>  6: -0.13176489
#>  7: -0.10430258
#>  8: -0.02536408
#>  9:  0.07122433
#> 10: -0.01180950
#> 11:  0.07122433
#> 12: -0.06367401
#> 
#> $intercept
#> [1] 20.49833
#> 
#> $x
#>       cyl  disp    hp
#>     <num> <num> <num>
#>  1:     4 120.1    97
#>  2:     8 318.0   150
#>  3:     8 304.0   150
#>  4:     8 350.0   245
#>  5:     8 400.0   175
#>  6:     4  79.0    66
#>  7:     4 120.3    91
#>  8:     4  95.1   113
#>  9:     8 351.0   264
#> 10:     6 145.0   175
#> 11:     8 301.0   335
#> 12:     4 121.0   109
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# sums to prediction
cbind(
  m_sum = rowSums(components$m) + components$intercept,
  prediction = predict(rpfit, test)
)
#>       m_sum    .pred
#> 1  25.93281 25.93281
#> 2  18.42495 18.42495
#> 3  18.42495 18.42495
#> 4  15.06680 15.06680
#> 5  18.34031 18.34031
#> 6  31.15026 31.15026
#> 7  29.40223 29.40223
#> 8  26.34941 26.34941
#> 9  15.06680 15.06680
#> 10 20.68766 20.68766
#> 11 15.06680 15.06680
#> 12 25.68539 25.68539

# Only get components with interactions of a lower degree, ignoring 3-way interactions
predict_components(rpfit, test, max_interaction = 2)
#> $m
#>            cyl       disp         hp    cyl:disp       cyl:hp     disp:hp
#>          <num>      <num>      <num>       <num>        <num>       <num>
#>  1:  0.9013797  3.2511514  1.3776657  0.08511109 -0.002332250 -0.11682939
#>  2: -0.6706529 -0.9186116 -0.4097990 -0.14816938 -0.079604722  0.10864962
#>  3: -0.6706529 -0.9186116 -0.4097990 -0.14816938 -0.079604722  0.10864962
#>  4: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#>  5: -0.6706529 -1.1532004 -0.3559023 -0.07494483 -0.084604722  0.12038493
#>  6:  0.9013797  4.6372466  4.7893287  0.17994554 -0.091148415  0.36694210
#>  7:  0.9013797  3.2511514  4.7893287  0.08511109 -0.091148415  0.07237951
#>  8:  0.9013797  4.6372466  0.6132023  0.17994554  0.006048960 -0.46138786
#>  9: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#> 10: -0.1689229  0.9540856 -0.3559023 -0.08618539 -0.045141470 -0.09680110
#> 11: -0.6706529 -0.9186116 -3.6554527 -0.14816938 -0.017505790 -0.09236340
#> 12:  0.9013797  3.2511514  0.9815403  0.08511109 -0.005884567  0.03742733
#> 
#> $intercept
#> [1] 20.49833
#> 
#> $x
#>       cyl  disp    hp
#>     <num> <num> <num>
#>  1:     4 120.1    97
#>  2:     8 318.0   150
#>  3:     8 304.0   150
#>  4:     8 350.0   245
#>  5:     8 400.0   175
#>  6:     4  79.0    66
#>  7:     4 120.3    91
#>  8:     4  95.1   113
#>  9:     8 351.0   264
#> 10:     6 145.0   175
#> 11:     8 301.0   335
#> 12:     4 121.0   109
#> 
#> $remainder
#>  [1] -0.06166953  0.04480135  0.04480135  0.07122433  0.06090002 -0.13176489
#>  [7] -0.10430258 -0.02536408  0.07122433 -0.01180950  0.07122433 -0.06367401
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# Only retrieve main effects
(main_effects <- predict_components(rpfit, test, max_interaction = 1))
#> $m
#>            cyl       disp         hp
#>          <num>      <num>      <num>
#>  1:  0.9013797  3.2511514  1.3776657
#>  2: -0.6706529 -0.9186116 -0.4097990
#>  3: -0.6706529 -0.9186116 -0.4097990
#>  4: -0.6706529 -0.9186116 -3.6554527
#>  5: -0.6706529 -1.1532004 -0.3559023
#>  6:  0.9013797  4.6372466  4.7893287
#>  7:  0.9013797  3.2511514  4.7893287
#>  8:  0.9013797  4.6372466  0.6132023
#>  9: -0.6706529 -0.9186116 -3.6554527
#> 10: -0.1689229  0.9540856 -0.3559023
#> 11: -0.6706529 -0.9186116 -3.6554527
#> 12:  0.9013797  3.2511514  0.9815403
#> 
#> $intercept
#> [1] 20.49833
#> 
#> $x
#>       cyl  disp    hp
#>     <num> <num> <num>
#>  1:     4 120.1    97
#>  2:     8 318.0   150
#>  3:     8 304.0   150
#>  4:     8 350.0   245
#>  5:     8 400.0   175
#>  6:     4  79.0    66
#>  7:     4 120.3    91
#>  8:     4  95.1   113
#>  9:     8 351.0   264
#> 10:     6 145.0   175
#> 11:     8 301.0   335
#> 12:     4 121.0   109
#> 
#> $remainder
#>  [1] -0.09572008 -0.07432313 -0.07432313 -0.18681424  0.02173540  0.32397433
#>  [7] -0.03796040 -0.30075744 -0.18681424 -0.23993745 -0.18681424  0.05297984
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# The difference is the combined contribution of interaction effects
cbind(
  m_sum = rowSums(main_effects$m) + main_effects$intercept,
  prediction = predict(rpfit, test)
)
#>       m_sum    .pred
#> 1  26.02853 25.93281
#> 2  18.49927 18.42495
#> 3  18.49927 18.42495
#> 4  15.25362 15.06680
#> 5  18.31858 18.34031
#> 6  30.82629 31.15026
#> 7  29.44019 29.40223
#> 8  26.65016 26.34941
#> 9  15.25362 15.06680
#> 10 20.92759 20.68766
#> 11 15.25362 15.06680
#> 12 25.63241 25.68539
```
