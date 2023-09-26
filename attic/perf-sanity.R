# Local install
# devtools::install(".")
# Remote install
# remotes::install_github("PlantedML/randomPlantedForest")

library(randomPlantedForest)
library(microbenchmark)
library(mlbench)
data("PimaIndiansDiabetes")

bm <- microbenchmark(logit = rpf(diabetes ~ ., PimaIndiansDiabetes, ntrees = 10, splits = 10, loss = "logit"),
                     L1 = rpf(diabetes ~ ., PimaIndiansDiabetes, ntrees = 10, splits = 10, loss = "L1"),
                     times = 1)

bm



# other -----------------------------------------------------------------------------------------------------------


library(randomPlantedForest)
rp <- rpf(mpg ~ ., data = mtcars)

rp$fit$get_parameters()
