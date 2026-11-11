# The components of a Random Planted Forest's prediction

Evaluates the components of a model at the rows of `new_data`: its
intercept and one column per main effect and interaction. They are a
functional decomposition of the prediction: for every row, the intercept
plus the components equals the prediction. For a pruned model, the
default of
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
these are its kept components.

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
  number of `target_levels`, the levels of the outcome: one column per
  component and level, named `<component>__class:<level>`. With
  `loss = "logit"` or `"exponential"` a component's columns sum to 0
  across the levels, as the scores of
  [`predict()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  with `type = "link"` do.

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

## See also

[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
for how the components are identified,
[`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
for the pruned model's components. The glex package plots the result
with [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Examples

``` r
# The default model: its kept components add up to its prediction
set.seed(1)
model <- rpf(mpg ~ wt + hp + cyl, data = mtcars)
components <- predict_components(model, mtcars[1:3, ])
components$m
#>            wt        hp        cyl
#>         <num>     <num>      <num>
#> 1:  0.3989752 0.2785334 -0.3703708
#> 2: -0.1242620 0.2785334 -0.3703708
#> 3:  1.6816513 1.1853078  2.4229759
components$intercept + rowSums(components$m)
#> [1] 20.50687 19.98363 25.48966
predict(model, mtcars[1:3, ])
#> # A tibble: 3 × 1
#>   .pred
#>   <dbl>
#> 1  20.5
#> 2  20.0
#> 3  25.5

# A forest: all its components up to its max_interaction
train <-  mtcars[1:20, 1:4]
test <-  mtcars[21:32, 1:4]

set.seed(23)
rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, n_families = 30, prune = FALSE)

# Extract all components, including main effects and interaction terms up to `max_interaction`
(components <- predict_components(rpfit, test))
#> $m
#>            cyl       disp         hp  cyl:disp      cyl:hp     disp:hp
#>          <num>      <num>      <num>     <num>       <num>       <num>
#>  1:  2.0588524  1.6373900  1.5625285 0.4989844 -0.24174939 -0.07853976
#>  2: -1.5316997 -1.2941107 -0.3316503 0.2765535  0.03790460  0.26998016
#>  3: -1.5316997 -1.2941107 -0.3316503 0.2765535  0.03790460  0.26998016
#>  4: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#>  5: -1.5316997 -1.0895381 -0.5228127 0.2383035  0.07278092  0.41129734
#>  6:  2.0588524  3.1270869  3.0127215 0.8284196  0.29270117  0.84404186
#>  7:  2.0588524  1.6373900  3.0127215 0.4989844  0.29270117  0.48801491
#>  8:  2.0588524  3.1270869  0.4660429 0.8284196 -0.15597044 -0.16323737
#>  9: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#> 10: -0.4691853  0.7975219 -0.5228127 0.1072327  0.15896496 -0.22308406
#> 11: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#> 12:  2.0588524  1.6373900  1.1297472 0.4989844 -0.21032282 -0.01398476
#>     cyl:disp:hp
#>           <num>
#>  1: -0.03937451
#>  2:  0.02770790
#>  3:  0.02770790
#>  4:  0.01243666
#>  5:  0.02770790
#>  6:  0.08060198
#>  7: -0.04889802
#>  8:  0.07121579
#>  9:  0.01243666
#> 10:  0.01560123
#> 11:  0.01243666
#> 12: -0.04470784
#> 
#> $intercept
#> [1] 20.35293
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
#> 1  25.75102 25.75102
#> 2  17.80761 17.80761
#> 3  17.80761 17.80761
#> 4  14.49275 14.49275
#> 5  17.95897 17.95897
#> 6  30.59735 30.59735
#> 7  28.29269 28.29269
#> 8  26.58534 26.58534
#> 9  14.49275 14.49275
#> 10 20.21717 20.21717
#> 11 14.49275 14.49275
#> 12 25.40889 25.40889

# Only get components with interactions of a lower degree, ignoring 3-way interactions
predict_components(rpfit, test, max_interaction = 2)
#> $m
#>            cyl       disp         hp  cyl:disp      cyl:hp     disp:hp
#>          <num>      <num>      <num>     <num>       <num>       <num>
#>  1:  2.0588524  1.6373900  1.5625285 0.4989844 -0.24174939 -0.07853976
#>  2: -1.5316997 -1.2941107 -0.3316503 0.2765535  0.03790460  0.26998016
#>  3: -1.5316997 -1.2941107 -0.3316503 0.2765535  0.03790460  0.26998016
#>  4: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#>  5: -1.5316997 -1.0895381 -0.5228127 0.2383035  0.07278092  0.41129734
#>  6:  2.0588524  3.1270869  3.0127215 0.8284196  0.29270117  0.84404186
#>  7:  2.0588524  1.6373900  3.0127215 0.4989844  0.29270117  0.48801491
#>  8:  2.0588524  3.1270869  0.4660429 0.8284196 -0.15597044 -0.16323737
#>  9: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#> 10: -0.4691853  0.7975219 -0.5228127 0.1072327  0.15896496 -0.22308406
#> 11: -1.5316997 -1.2941107 -3.1737806 0.2765535 -0.28273580  0.13315752
#> 12:  2.0588524  1.6373900  1.1297472 0.4989844 -0.21032282 -0.01398476
#> 
#> $intercept
#> [1] 20.35293
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
#>  [1] -0.03937451  0.02770790  0.02770790  0.01243666  0.02770790  0.08060198
#>  [7] -0.04889802  0.07121579  0.01243666  0.01560123  0.01243666 -0.04470784
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# Only retrieve main effects
(main_effects <- predict_components(rpfit, test, max_interaction = 1))
#> $m
#>            cyl       disp         hp
#>          <num>      <num>      <num>
#>  1:  2.0588524  1.6373900  1.5625285
#>  2: -1.5316997 -1.2941107 -0.3316503
#>  3: -1.5316997 -1.2941107 -0.3316503
#>  4: -1.5316997 -1.2941107 -3.1737806
#>  5: -1.5316997 -1.0895381 -0.5228127
#>  6:  2.0588524  3.1270869  3.0127215
#>  7:  2.0588524  1.6373900  3.0127215
#>  8:  2.0588524  3.1270869  0.4660429
#>  9: -1.5316997 -1.2941107 -3.1737806
#> 10: -0.4691853  0.7975219 -0.5228127
#> 11: -1.5316997 -1.2941107 -3.1737806
#> 12:  2.0588524  1.6373900  1.1297472
#> 
#> $intercept
#> [1] 20.35293
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
#>  [1] 0.13932070 0.61214613 0.61214613 0.13941184 0.75008963 2.04576466
#>  [7] 1.23080242 0.58042762 0.13941184 0.05871484 0.13941184 0.22996894
#> 
#> attr(,"class")
#> [1] "glex"           "rpf_components" "list"          

# The difference is the combined contribution of interaction effects
cbind(
  m_sum = rowSums(main_effects$m) + main_effects$intercept,
  prediction = predict(rpfit, test)
)
#>       m_sum    .pred
#> 1  25.61170 25.75102
#> 2  17.19547 17.80761
#> 3  17.19547 17.80761
#> 4  14.35334 14.49275
#> 5  17.20888 17.95897
#> 6  28.55159 30.59735
#> 7  27.06189 28.29269
#> 8  26.00491 26.58534
#> 9  14.35334 14.49275
#> 10 20.15845 20.21717
#> 11 14.35334 14.49275
#> 12 25.17892 25.40889
```
