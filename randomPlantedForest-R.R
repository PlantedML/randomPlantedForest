# load external libraries ------------------------
library(inline)
library(Rcpp)
library(RcppParallel)
library(rbenchmark)
sourceCpp(paste(rpf_path, '/randomPlantedForest-R.cpp', sep=''))


# wrap class constructor in R function
new_rpf <- function(Y, X, max_interaction=1, ntrees=50, splits=30, split_try=10, t_try=0.4, deterministic=FALSE, parallel=FALSE, purify=FALSE, cv=FALSE){
  rpf_cpp <- new(RandomPlantedForest, Y, X, max_interaction, ntrees, splits, c(split_try, t_try, purify, deterministic, parallel, cv))
  return(rpf_cpp)
}
