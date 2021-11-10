# load external libraries ------------------------
library(inline)
library(Rcpp)
library(RcppParallel)
sourceCpp(paste(rpf_path, '/randomPlantedForest.cpp', sep=''))


# wrap class constructor in R function
new_rpf <- function(Y, X, max_interaction=1, ntrees=50, splits=30, split_try=10, t_try=0.4, 
                    deterministic=FALSE, parallel=FALSE, purify=FALSE, cv=FALSE,
                    loss='L2', delta=0, epsilon=0.1){
  if(!missing(loss) | !missing(delta) | !missing(epsilon)){
    return(new(ClassificationRPF, Y, X, loss, c(max_interaction, ntrees, splits, split_try, t_try, 
                                                purify, deterministic, parallel, cv, delta, epsilon)))
  }
  return(new(RandomPlantedForest, Y, X, c(max_interaction, ntrees, splits, split_try, t_try, 
                                          purify, deterministic, parallel, cv)))
}
