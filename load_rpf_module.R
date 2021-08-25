library(inline)
library(Rcpp)


# enable compiler optimization
myPlugin <- getPlugin("Rcpp")
myPlugin$PKG_CXXFLAGS <- c(myPlugin$PKG_CXXFLAGS, "-O3 -funroll-loops -fvectorize")


# set path to code
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)


# source code
source(paste(rpf_path, '/ARCHIVE_Codes_used_in_paper/generate_data.R', sep=''))
source(paste(rpf_path, '/predict_rpf.R', sep=''))
source(paste(rpf_path, '/purify_rpf.R', sep=''))
source(paste(rpf_path, '/rpf.R', sep=''))
sourceCpp(paste(rpf_path, '/C-Code.cpp', sep=''))
sourceCpp(paste(rpf_path, '/randomPlantedForest-R.cpp', sep=''))


# generate test data
data <- generate_data(Model=1)
samples_x <- data$X # extract samples
samples_y <- data$Y_true


# create new class object
n_splits <- 15
max_inter <- 2
n_trees <- 1
t_try <- 0.5
d <- TRUE
rpf_cpp <- new(RandomPlantedForest, data$Y_start, data$X, max_inter, n_trees, n_splits, t_try)
rpf_cpp$set_deterministic(d)
rpf_R <- rpf(data$Y_start, data$X, max_interaction = max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, deterministic=d)


# cross-validation
rpf_cpp$cross_validation(2, c(5,50), c(0.2,0.5,0.7,0.9), c(1,2,5,10))


# purify
rpf_R_purified <- rpf_purify(rpf_R)
rpf_cpp$purify()


# view result
rpf_cpp$print()
rpf_R 


# predict
predictions_cpp <- rpf_cpp$predict_matrix(samples_x, c(0))
predictions_R <- predict_rpf(samples_x, rpf_R, c(0))

predictions_cpp <- rpf_cpp$predict_matrix(samples_x[, 1:2], c(1,2))
predictions_R <- predict_rpf(samples_x[, 1:2], rpf_R, c(1,2))

predictions_cpp <- rpf_cpp$predict_vector(samples_x[, 1], c(1))
predictions_R <- predict_rpf(samples_x[, 1], rpf_R, c(1))

predictions_cpp <- rpf_cpp$predict_vector(samples_x[1, ], c(0))
predictions_R <- predict_rpf(samples_x[1, ], rpf_R, c(0))

predictions_cpp <- rpf_cpp$predict_matrix(samples_x, c(0))
predictions_R <- predict_rpf(samples_x, rpf_R, c(0))


# accuracy
variation <- mean(data$Y_true^2)
MSE_cpp <- sum((samples_y - predictions_cpp)^2) / length(samples_y)
MSE_R <- sum((samples_y - predictions_R)^2) / length(samples_y)
MSE_cpp
MSE_R


# benchmark
library(rbenchmark)
benchmark( "rpf_cpp" = {
                rpf_cpp <- new(RandomPlantedForest, data$Y_start, data$X, max_inter, n_trees, n_splits, t_try)
                predictions_cpp <- rpf_cpp$predict_matrix(samples_x, c(0))
           },
           "rpf_R" = {
                rpf_R <- rpf(data$Y_start, data$X, max_interaction = max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, deterministic=FALSE)
                predictions_R <- predict_rpf(samples_x, rpf_R, c(0))
           },
           replications=1
)


