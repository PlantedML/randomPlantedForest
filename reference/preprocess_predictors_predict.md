# Preprocess predictors for prediction

Convert logical and character columns to appropriate types, re-order
factor levels to match the ordering learned during fitting (stored in
`object$factor_levels`), re-encode factor columns as integers, and
return a numeric matrix suitable for the underlying C++ prediction
routines.

## Usage

``` r
preprocess_predictors_predict(object, predictors)
```

## Arguments

- object:

  An object of class `rpf` returned by
  [`rpf()`](http://plantedml.com/randomPlantedForest/reference/rpf.md).

- predictors:

  A data frame or matrix of predictor values to preprocess.

## Value

A numeric matrix with the same number of rows as `predictors`.

## Details

This is primarily an internal utility used by
[`predict()`](https://rdrr.io/r/stats/predict.html) methods but is
exported to support advanced users and tests.

## Examples

``` r
rpfit <- rpf(x = mtcars[, c("cyl", "wt")], y = mtcars$mpg)
processed <- hardhat::forge(mtcars[, c("cyl", "wt")], rpfit$blueprint)
X <- preprocess_predictors_predict(rpfit, processed$predictors)
dim(X)
#> [1] 32  2
```
