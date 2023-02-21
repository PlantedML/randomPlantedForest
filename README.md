
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
library(ggplot2)
library(patchwork)

mtcars$cyl <- factor(mtcars$cyl)
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
#> 1 21.27391 21.0
#> 2 21.13418 21.0
#> 3 24.05330 22.8
#> 4 20.59945 21.4
#> 5 17.89879 18.7
#> 6 18.73767 18.1
```

Prediction components can be accessed via `extract_components`,
including the intercept, main effects, and interactions up to a
specified degree. The returned object also contains the original data as
`x`, which is required for visualization.

``` r
components <- extract_components(rpfit, new_data = mtcars) 

str(components)
#> List of 2
#>  $ m:Classes 'data.table' and 'data.frame':  32 obs. of  7 variables:
#>   ..$ intercept: num [1:32] 20.2 20.2 20.2 20.2 20.2 ...
#>   ..$ cyl      : num [1:32] 0.383 0.383 1.505 0.383 -1.245 ...
#>   ..$ wt       : num [1:32] 0.515 0.152 1.567 -0.382 -1.137 ...
#>   ..$ hp       : num [1:32] 0.318 0.318 0.665 0.318 -0.122 ...
#>   ..$ cyl:wt   : num [1:32] -0.0773 -0.0123 0.2817 0.1982 -0.0894 ...
#>   ..$ cyl:hp   : num [1:32] -0.2054 -0.2054 -0.1045 -0.2054 -0.0994 ...
#>   ..$ hp:wt    : num [1:32] 0.0901 0.2489 -0.11 0.0375 0.3418 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ x:Classes 'data.table' and 'data.frame':  32 obs. of  3 variables:
#>   ..$ cyl: Factor w/ 3 levels "4","6","8": 2 2 1 2 3 2 3 1 1 2 ...
#>   ..$ wt : num [1:32] 2.62 2.88 2.32 3.21 3.44 ...
#>   ..$ hp : num [1:32] 110 110 93 110 175 105 245 62 95 123 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  - attr(*, "class")= chr [1:2] "rpf_components" "list"
```

Various visualization options are available, e.g. for main and
second-order interaction effects:

``` r
p1 <- autoplot(components, "wt")
p2 <- autoplot(components, "hp")
p3 <- autoplot(components, "cyl")
p4 <- autoplot(components, c("wt", "hp"))

(p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Selected effects for mtcars",
    caption = "(It's a tiny dataset but it has to fit in a README, okay?)"
  )
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

See the [Bikesharing: Decomposition
Example](http://plantedml.com/randomPlantedForest/articles/Bikesharing-Decomposition-Example.html)
article for more examples.
