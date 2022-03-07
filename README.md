# randomPlantedForest

<!-- badges: start -->
[![R-CMD-check](https://github.com/PlantedML/randomPlantedForest/workflows/R-CMD-check/badge.svg)](https://github.com/PlantedML/randomPlantedForest/actions)
<!-- badges: end -->

The goal of randomPlantedForest is to ...

## Installation

You can install the development version of randomPlantedForest from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PlantedML/randomPlantedForest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(randomPlantedForest)

rp_fit <- rpf(mpg ~ cyl + wt, data = mtcars)
```

