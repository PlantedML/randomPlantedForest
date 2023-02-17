
<!-- README.md is generated from README.Rmd. Please edit that file -->

# randomPlantedForest

<!-- badges: start -->

[![R-CMD-check](https://github.com/PlantedML/randomPlantedForest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PlantedML/randomPlantedForest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/PlantedML/randomPlantedForest/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PlantedML/randomPlantedForest?branch=master)
<!-- badges: end -->

`randomPlantedForest` implements “Random Planted Forest”, a directly
interpretable tree ensemble [(arxiv)](https://arxiv.org/abs/2012.14563).

## Installation

You can install the development version of `randomPlantedForest` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("PlantedML/randomPlantedForest")
```

## Example

Model fitting uses a familiar interface:

``` r
library(randomPlantedForest)

rpfit <- rpf(mpg ~ cyl + wt + hp, data = mtcars, ntrees = 25, max_interaction = 2)
rpfit
#> -- Regression Random Planted Forest --
#> 
#> Formula: mpg ~ cyl + wt + hp 
#> Fit using 3 predictors and 2-degree interactions.
#> Forest is _not_ purified!
#> 
#> Called with parameters:
#> 
#>             loss: L2
#>           ntrees: 25
#>  max_interaction: 2
#>           splits: 30
#>        split_try: 10
#>            t_try: 0.4
#>            delta: 0
#>          epsilon: 0.1
#>    deterministic: FALSE
#>         parallel: FALSE
#>           purify: FALSE
#>               cv: FALSE

predict(rpfit, new_data = mtcars) |>
  cbind(mpg = mtcars$mpg) |>
  head()
#>      .pred  mpg
#> 1 20.96153 21.0
#> 2 20.73825 21.0
#> 3 26.43672 22.8
#> 4 20.89290 21.4
#> 5 18.09168 18.7
#> 6 18.19773 18.1
```

Prediction components can be accessed via `extract_components`,
including the intercept, main effects, and interactions up to a
specified degree:

``` r
m <- extract_components(rpfit, new_data = mtcars) 

str(m)
#> Classes 'data.table' and 'data.frame':   32 obs. of  7 variables:
#>  $ intercept: num  20.2 20.2 20.2 20.2 20.2 ...
#>  $ cyl      : num  -0.72 -0.72 2.35 -0.72 -1.93 ...
#>  $ wt       : num  0.16 -0.28 4.111 -0.301 -1.502 ...
#>  $ hp       : num  0.51 0.51 0.509 0.51 0.478 ...
#>  $ cyl:wt   : num  0.432 0.398 0.401 0.581 0.282 ...
#>  $ cyl:hp   : num  0.417 0.417 -0.479 0.417 0.5 ...
#>  $ hp:wt    : num  -0.0355 0.2149 -0.6549 0.2074 0.0648 ...
#>  - attr(*, ".internal.selfref")=<externalptr>
```
