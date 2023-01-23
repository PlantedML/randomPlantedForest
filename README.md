
<!-- README.md is generated from README.Rmd. Please edit that file -->

# randomPlantedForest

<!-- badges: start -->

[![R-CMD-check](https://github.com/PlantedML/randomPlantedForest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PlantedML/randomPlantedForest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/PlantedML/randomPlantedForest/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PlantedML/randomPlantedForest?branch=master)
<!-- badges: end -->

`randomPlantedForest` …

## Installation

You can install the development version of `randomPlantedForest` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("PlantedML/randomPlantedForest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(randomPlantedForest)

rpfit <- rpf(mpg ~ cyl + wt, data = mtcars)
rpfit
#> -- Regression Random Planted Forest --
#> 
#> Formula: mpg ~ cyl + wt 
#> Fit using 2 predictors and main effects only.
#> Called with parameters:
#> 
#>             loss: L2
#>           ntrees: 50
#>  max_interaction: 1
#>           splits: 30
#>        split_try: 10
#>            t_try: 0.4
#>            delta: 0
#>          epsilon: 0.1
#>    deterministic: FALSE
#>         parallel: FALSE
#>           purify: FALSE
#>               cv: FALSE

predict(rpfit, new_data = mtcars[1:10, ])
#> # A tibble: 10 × 1
#>    .pred
#>    <dbl>
#>  1  20.4
#>  2  20.4
#>  3  25.6
#>  4  20.5
#>  5  17.0
#>  6  18.4
#>  7  15.2
#>  8  22.9
#>  9  23.4
#> 10  19.0
```
