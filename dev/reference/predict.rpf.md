# Predict with a Random Planted Forest

Predictions of a forest (class `rpf`, from
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
with `prune = FALSE`) for new data: the mean of its tree families'
predictions. For classification, `type` chooses probabilities, classes
or the raw scores. A pruned model, the default of
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md),
has its own method,
[`predict.rpf_pruned()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md),
with the same arguments.

## Usage

``` r
# S3 method for class 'rpf'
predict(
  object,
  new_data,
  type = ifelse(object$mode == "regression", "numeric", "prob"),
  nthreads = NULL,
  ...
)
```

## Arguments

- object:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- new_data:

  `[data.frame | matrix]`: New observations to predict.

- type:

  `[character(1)]`: `"numeric"` for regression outcomes, `"class"` for
  class predictions or `"prob"` for probability predictions. Defaults to
  `"numeric"` for regression and `"prob"` for classification.

  For classification and `loss = "L1"` or `"L2"`, `"numeric"` yields raw
  predictions which are not guaranteed to be valid probabilities in
  `[0, 1]`. For `type = "prob"`, these are truncated to ensure this
  property.

  For classification, `type = "link"` is an alias for
  `type = "numeric"`. If `loss` is `"logit"` or `"exponential"`, the raw
  predictions are the scores on the link scale, like the linear
  predictor of a [`glm`](https://rdrr.io/r/stats/glm.html): with two
  classes one score, the log-odds of the second class; with more, one
  score per class (`"exponential"`, summing to 0) or per class but the
  first (`"logit"`, whose first class has score 0). The probabilities of
  `type = "prob"` are their softmax (see "Classification losses" in
  [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)).
  If `loss` is `"L1"` or `"L2"`, the link is the identity.

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with, capped at the available cores.

- ...:

  Not currently used, but required for extensibility. Unknown arguments
  are an error.

## Value

For regression: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with column
`.pred` with the same number of rows as `new_data`.

For classification: A
[`tbl`](https://tibble.tidyverse.org/reference/tibble.html) with one
column `.pred_<level>` for each level in `y` containing class
probabilities if `type = "prob"`. For `type = "class"`, one column
`.pred_class` with class predictions is returned. For `type = "numeric"`
or `"link"`, raw predictions are returned: one column `.pred` for binary
outcomes (the log-odds of the second level with `loss = "logit"` or
`"exponential"`), and one column `.pred_<level>` per level for
multiclass outcomes. With `loss = "logit"` or `"exponential"`, these are
one score per class that sum to 0 across the classes, and the
probabilities are their softmax. (The multiclass logit loss fits `K - 1`
scores against the first level; they are reported with the first level's
score 0 joined and centred across the `K` classes, which leaves the
probabilities as they are.)

## Details

`new_data` is prepared as the training data were: by the formula or
recipe of the fit, and with factor levels coded as in the fit (see
"Factor predictors" in
[`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)).
Missing values are an error.

## See also

[`predict.rpf_pruned()`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
for the pruned model,
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
for the prediction as a sum of components.

## Examples

``` r
# Regression: a tibble with the column .pred
rpfit <- rpf(mpg ~ cyl + wt, data = mtcars, n_families = 20, prune = FALSE)
predict(rpfit, mtcars[1:3, ])
#> # A tibble: 3 × 1
#>   .pred
#>   <dbl>
#> 1  21.1
#> 2  20.2
#> 3  24.2

# Classification: probabilities by default, or classes, or the raw scores
rpfit <- rpf(Species ~ ., data = iris, n_families = 20, prune = FALSE)
predict(rpfit, iris[c(1, 51, 101), ])
#> # A tibble: 3 × 3
#>   .pred_setosa .pred_versicolor .pred_virginica
#>          <dbl>            <dbl>           <dbl>
#> 1      0.996            0.00275        0.000944
#> 2      0.0125           0.965          0.0226  
#> 3      0.00316          0.00630        0.991   
predict(rpfit, iris[c(1, 51, 101), ], type = "class")
#> # A tibble: 3 × 1
#>   .pred_class
#>   <fct>      
#> 1 setosa     
#> 2 versicolor 
#> 3 virginica  
predict(rpfit, iris[c(1, 51, 101), ], type = "link")
#> # A tibble: 3 × 3
#>   .pred_setosa .pred_versicolor .pred_virginica
#>          <dbl>            <dbl>           <dbl>
#> 1         4.28            -1.61           -2.68
#> 2        -1.65             2.70           -1.05
#> 3        -2.15            -1.46            3.60
```
