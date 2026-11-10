# Extract predicted components from a Random Planted Forest

Prediction components are a functional decomposition of the model
prediction. The sum of all components equals the overall predicted value
for an observation.

## Usage

``` r
predict_components(object, new_data, ...)

# S3 method for class 'rpf'
predict_components(
  object,
  new_data,
  max_interaction = NULL,
  predictors = NULL,
  ...
)

# S3 method for class 'rpf_pruned'
predict_components(
  object,
  new_data,
  max_interaction = NULL,
  predictors = NULL,
  ...
)
```

## Arguments

- object:

  `[rpf | rpf_pruned]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model, or a pruned one from
  [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  or
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).

- new_data:

  `[data.frame | matrix]`: New observations to predict.

- ...:

  Not currently used, but required for extensibility. Unknown arguments
  are an error.

- max_interaction:

  `[integer(1) | NULL: NULL]`: Maximum degree of interactions to
  consider. `NULL` uses the `max_interaction` parameter from the
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  object. Must be between `1` (main effects only) and the
  `max_interaction` of the
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  object.

- predictors:

  `[character | NULL: NULL]`: Vector of one or more names of the
  forest's predictors to extract components for. These are the columns
  the forest was fitted on: with a formula, its terms, such as
  `log(x1)`. If `NULL`, all predictors and their interactions are
  returned, in the order of the fit's columns.

## Value

A `list` with elements:

- `m`
  ([`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)):
  Components for each main effect and interaction term, representing the
  functional decomposition of the prediction. All components together
  with the intercept sum up to the prediction. For multiclass
  classification, the number of output columns is multiplied by the
  number of `target_levels`: the levels of the outcome, or for
  `loss = "logit"` all levels but the first (the reference level).

- `intercept` (`numeric`): Expected value of the prediction. For
  multiclass classification, one value per target level, named by
  `target_levels`.

- `x`
  ([`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)):
  The predictors selected by `predictors`, as the forest sees them in
  `new_data` (with a formula, its terms, such as `log(x1)`).

- `target_levels` (`character`): For multiclass classification only:
  Vector of target levels which can be used to disassemble `m`, as names
  include both term and target level.

- `remainder` (`numeric`): Only if `m` leaves out components of the
  forest: if `predictors` leaves out some of the forest's predictors, if
  `max_interaction` (lowered to the number of `predictors`) is below
  that of the
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  object, or if the forest was purified with a `maxp_interaction` below
  it (see
  [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)).
  The prediction minus `intercept` and `m`'s row sums, so that the three
  sum to the prediction. For multiclass classification, one column per
  target level.

## Details

Extracts all possible components up to `max_interaction` degrees, up to
the value set when calling
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md).
The intercept is always included. Optionally `predictors` can be
specified to only include components made of the given predictors; the
sum of the other components is returned as `remainder`. If
`max_interaction` is greater than `length(predictors)`, the
`max_interaction` will be lowered accordingly.

All components are evaluated in one pass over the rows of `new_data`.
Components that the purified forest does not hold, such as sets of
predictors that never appear together in a tree, are zero columns.

The components are only unique for a purified forest;
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
says how they are identified. An unpurified forest is purified in place
on first use with the defaults of
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md).
For other settings, such as more threads, call
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
before `predict_components()`.

## Note

Depending on the number of predictors and `max_interaction`, the number
of components will increase drastically to
`sum(choose(ncol(new_data), seq_len(max_interaction)))`.

## Pruned forests

For an `rpf_pruned` object from
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
or
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
`m` holds only the kept components, each scaled by its coefficient
`beta_S`, and `intercept` is the pruned model's intercept `alpha`, so
that they sum to
[`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
of the pruned model. Components that are not kept are absent, not zero
columns. If `predictors` or `max_interaction` leave out kept components,
their sum is returned as `remainder`.

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
#>             cyl       disp         hp    cyl:disp      cyl:hp      disp:hp
#>           <num>      <num>      <num>       <num>       <num>        <num>
#>  1:  1.00484882  3.3276221  1.3422678  0.03379633 -0.02689605 -0.151403172
#>  2: -0.56783894 -0.9030339 -0.4562557 -0.14324922 -0.10287188  0.090570784
#>  3: -0.56783894 -0.9030339 -0.4562557 -0.14324922 -0.10287188  0.090570784
#>  4: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#>  5: -0.56783894 -1.1176398 -0.4023590 -0.09728299 -0.10787188  0.113283341
#>  6:  1.00484882  4.7249899  4.7966008  0.11735816 -0.14405930  0.292556145
#>  7:  1.00484882  3.3276221  4.7966008  0.03379633 -0.14405930  0.005578045
#>  8:  1.00484882  4.7249899  0.5737534  0.11735816 -0.01448719 -0.503522742
#>  9: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#> 10: -0.06601026  0.9930589 -0.4023590 -0.10044105 -0.06850736 -0.124462799
#> 11: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#> 12:  1.00484882  3.3276221  0.9461424  0.03379633 -0.03044837  0.002853543
#>     cyl:disp:hp
#>           <num>
#>  1: -0.06416397
#>  2:  0.04088736
#>  3:  0.04088736
#>  4:  0.06802170
#>  5:  0.05328418
#>  6: -0.10877009
#>  7: -0.08889227
#>  8: -0.02027404
#>  9:  0.06802170
#> 10: -0.01036044
#> 11:  0.06802170
#> 12: -0.06616846
#> 
#> $intercept
#> [1] 20.46674
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
#>             cyl       disp         hp    cyl:disp      cyl:hp      disp:hp
#>           <num>      <num>      <num>       <num>       <num>        <num>
#>  1:  1.00484882  3.3276221  1.3422678  0.03379633 -0.02689605 -0.151403172
#>  2: -0.56783894 -0.9030339 -0.4562557 -0.14324922 -0.10287188  0.090570784
#>  3: -0.56783894 -0.9030339 -0.4562557 -0.14324922 -0.10287188  0.090570784
#>  4: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#>  5: -0.56783894 -1.1176398 -0.4023590 -0.09728299 -0.10787188  0.113283341
#>  6:  1.00484882  4.7249899  4.7966008  0.11735816 -0.14405930  0.292556145
#>  7:  1.00484882  3.3276221  4.7966008  0.03379633 -0.14405930  0.005578045
#>  8:  1.00484882  4.7249899  0.5737534  0.11735816 -0.01448719 -0.503522742
#>  9: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#> 10: -0.06601026  0.9930589 -0.4023590 -0.10044105 -0.06850736 -0.124462799
#> 11: -0.56783894 -0.9030339 -3.7116755 -0.14324922 -0.04570978 -0.096450642
#> 12:  1.00484882  3.3276221  0.9461424  0.03379633 -0.03044837  0.002853543
#> 
#> $intercept
#> [1] 20.46674
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
#>  [1] -0.06416397  0.04088736  0.04088736  0.06802170  0.05328418 -0.10877009
#>  [7] -0.08889227 -0.02027404  0.06802170 -0.01036044  0.06802170 -0.06616846
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# Only retrieve main effects
(main_effects <- predict_components(rpfit, test, max_interaction = 1))
#> $m
#>             cyl       disp         hp
#>           <num>      <num>      <num>
#>  1:  1.00484882  3.3276221  1.3422678
#>  2: -0.56783894 -0.9030339 -0.4562557
#>  3: -0.56783894 -0.9030339 -0.4562557
#>  4: -0.56783894 -0.9030339 -3.7116755
#>  5: -0.56783894 -1.1176398 -0.4023590
#>  6:  1.00484882  4.7249899  4.7966008
#>  7:  1.00484882  3.3276221  4.7966008
#>  8:  1.00484882  4.7249899  0.5737534
#>  9: -0.56783894 -0.9030339 -3.7116755
#> 10: -0.06601026  0.9930589 -0.4023590
#> 11: -0.56783894 -0.9030339 -3.7116755
#> 12:  1.00484882  3.3276221  0.9461424
#> 
#> $intercept
#> [1] 20.46674
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
#>  [1] -0.20866686 -0.11466295 -0.11466295 -0.21738794 -0.03858734  0.15708491
#>  [7] -0.19357719 -0.42092581 -0.21738794 -0.30377165 -0.21738794 -0.05996695
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# The difference is the combined contribution of interaction effects
cbind(
  m_sum = rowSums(main_effects$m) + main_effects$intercept,
  prediction = predict(rpfit, test)
)
#>       m_sum    .pred
#> 1  26.14148 25.93281
#> 2  18.53961 18.42495
#> 3  18.53961 18.42495
#> 4  15.28419 15.06680
#> 5  18.37890 18.34031
#> 6  30.99318 31.15026
#> 7  29.59581 29.40223
#> 8  26.77033 26.34941
#> 9  15.28419 15.06680
#> 10 20.99143 20.68766
#> 11 15.28419 15.06680
#> 12 25.74535 25.68539
```
