
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
  cbind(mpg = mtcars$mpg)
#>       .pred  mpg
#> 1  20.95646 21.0
#> 2  20.90617 21.0
#> 3  26.14025 22.8
#> 4  20.47584 21.4
#> 5  18.16133 18.7
#> 6  18.39644 18.1
#> 7  14.30703 14.3
#> 8  25.22488 24.4
#> 9  22.98006 22.8
#> 10 18.60618 19.2
#> 11 18.60618 17.8
#> 12 16.10386 16.4
#> 13 16.82663 17.3
#> 14 16.03545 15.2
#> 15 12.42553 10.4
#> 16 11.67859 10.4
#> 17 13.50277 14.7
#> 18 30.64132 32.4
#> 19 31.32540 30.4
#> 20 32.43137 33.9
#> 21 22.54169 21.5
#> 22 16.19924 15.5
#> 23 16.89015 15.2
#> 24 13.52060 13.3
#> 25 17.57316 19.2
#> 26 29.72913 27.3
#> 27 27.40128 26.0
#> 28 29.09396 30.4
#> 29 15.83796 15.8
#> 30 19.89757 19.7
#> 31 14.69952 15.0
#> 32 23.01672 21.4
```

Prediction components can be accessed via `extract_components`,
including the intercept, main effects, and interactions up to a
specified degree:

``` r
m <- extract_components(rpfit, new_data = mtcars) 
m
#>     intercept        cyl         wt          hp     cyl:wt       cyl:hp
#>         <num>      <num>      <num>       <num>      <num>        <num>
#>  1:  20.37915 -0.7568261  0.4725960  0.64748659  0.3808161  0.084108813
#>  2:  20.37915 -0.7568261  0.4572379  0.64748659  0.3049704  0.084108813
#>  3:  20.37915  2.1206344  1.1162665  1.98008781  0.4373314  0.089280739
#>  4:  20.37915 -0.7568261 -0.4704765  0.64748659  0.6325637  0.084108813
#>  5:  20.37915 -1.9178420 -1.0476970 -0.10151336  0.6240551  0.112309152
#>  6:  20.37915 -0.7568261 -1.2139024  0.03164736  0.1773348 -0.008819774
#>  7:  20.37915 -1.9178420 -1.7271224 -3.06977724  0.5609756 -0.092673146
#>  8:  20.37915  2.1206344 -0.5136003  3.61740719 -0.1993157  0.029001176
#>  9:  20.37915  2.1206344  0.1435016  0.83031706 -0.1340190 -0.102349280
#> 10:  20.37915 -0.7568261 -1.0476970 -0.35763152  0.2072313  0.022276237
#> 11:  20.37915 -0.7568261 -1.0476970 -0.35763152  0.2072313  0.022276237
#> 12:  20.37915 -1.9178420 -1.6356193 -1.31959065  0.5303389  0.138716810
#> 13:  20.37915 -1.9178420 -1.4738336 -1.31959065  0.5775431  0.138716810
#>           hp:wt
#>           <num>
#>  1: -0.25087030
#>  2: -0.20995715
#>  3:  0.01750237
#>  4: -0.04016400
#>  5:  0.11287083
#>  6: -0.21214157
#>  7:  0.17431676
#>  8: -0.20839316
#>  9: -0.25716823
#> 10:  0.15968229
#> 11:  0.15968229
#> 12: -0.07129199
#> 13:  0.44249048
#>  [ reached getOption("max.print") -- omitted 20 rows ]
```
