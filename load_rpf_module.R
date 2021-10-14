library(inline)
library(Rcpp)
library(RcppParallel)


# enable compiler optimization ------------------------
myPlugin <- getPlugin("Rcpp")
myPlugin$PKG_CXXFLAGS <- c(myPlugin$PKG_CXXFLAGS, "-O3 -funroll-loops")


# set path to code ------------------------
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)


# source code ------------------------
source(paste(rpf_path, '/ARCHIVE_Codes_used_in_paper/generate_data.R', sep=''))
source(paste(rpf_path, '/predict_rpf.R', sep=''))
source(paste(rpf_path, '/purify_rpf.R', sep=''))
source(paste(rpf_path, '/rpf.R', sep=''))
sourceCpp(paste(rpf_path, '/C-Code.cpp', sep=''))
sourceCpp(paste(rpf_path, '/randomPlantedForest-R.cpp', sep=''))


# generate test data ------------------------
sample_size = 500
data <- generate_data(Model=1, p=4, n=sample_size)
test_size = floor(length(data$Y_true) / 5)
x_test <- data$X[1:test_size, ] # extract samples
y_test <- data$Y_true[1:test_size]
x_train <- data$X[(test_size+1):sample_size, ] # extract samples
y_train <- data$Y_start[(test_size+1):sample_size]


# set parameters ------------------------
n_splits <- 15
max_inter <- 1
n_trees <- 50
split_try <- 10
t_try <- 0.5
deterministic_forest <- FALSE
parallel <- FALSE
purify_forest <- FALSE


# train models ------------------------
rpf_cpp <- new(RandomPlantedForest, y_train, x_train, max_inter, n_trees, n_splits, c(split_try, t_try, purify_forest, deterministic_forest, parallel))
#rpf_cpp$set_deterministic(d)
rpf_R <- rpf(y_train, x_train, max_interaction=max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, split_try = split_try, deterministic=deterministic_forest)


# cross-validation ------------------------
rpf_cpp$cross_validation(5, c(5, 20), c(0.2,0.5,0.7,0.9), c(1,2,5,10))


# purify ------------------------
rpf_R_purified <- rpf_purify(rpf_R)
rpf_cpp$purify()


# view result ------------------------
rpf_cpp$print()
rpf_cpp$get_parameters()
rpf_R 


# predict ------------------------

predictions_cpp <- rpf_cpp$predict_matrix(x_test[, 1:2], c(1,2))
predictions_R <- predict_rpf(x_test[, 1:2], rpf_R, c(1,2))

predictions_cpp <- rpf_cpp$predict_vector(x_test[, 1], c(1))
predictions_R <- predict_rpf(x_test[, 1], rpf_R, c(1))

predictions_cpp <- rpf_cpp$predict_vector(x_test[1, ], c(0))
predictions_R <- predict_rpf(x_test[1, ], rpf_R, c(0))

predictions_cpp <- rpf_cpp$predict_matrix(x_test, c(0))
predictions_R <- predict_rpf(x_test, rpf_R, c(0))

predictions_cpp <- rpf_cpp$predict_matrix(x_test, c(0))
predictions_R <- predict_rpf(x_test, rpf_R, c(0))


# accuracy ------------------------
variation <- mean(y_test^2)
MSE_cpp <- rpf_cpp$MSE(predictions_cpp, y_test) 
MSE_R <- sum((y_test - predictions_R)^2) / length(y_test)
MSE_cpp
MSE_R
variation


# benchmark ------------------------
library(rbenchmark)
benchmark( "rpf_cpp" = {
                rpf_cpp <- new(RandomPlantedForest, y_train, x_train, max_inter, n_trees, n_splits, t_try)
           },
           "rpf_cpp_predict" = {
             predictions_cpp <- rpf_cpp$predict_matrix(x_test, c(0))
           },
           "rpf_R" = {
                rpf_R <- rpf(y_train, x_train, max_interaction = max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, deterministic=FALSE)
           },
           "rpf_R_predict" = {
                predictions_R <- predict_rpf(x_test, rpf_R, c(0))
           },
           replications=1
)

benchmark( "rpf_cpp_cv" = {
              rpf_cpp$cross_validation(3, c(5,50), c(0.2,0.5,0.7,0.9), c(1,2,5,10))
            },
            replications=1
)

benchmark(  "rpf_cpp_sequential" = {
                rpf_cpp$set_parallel(FALSE);
                rpf_cpp$set_data(y_train, x_train)
            },
            "rpf_cpp_parallel" = {
                rpf_cpp$set_parallel(TRUE);
                rpf_cpp$set_data(y_train, x_train)
            },
            "rpf_R" = {
                rpf_R <- rpf(y_train, x_train, max_interaction = max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, deterministic=FALSE)
            },
            replications=1
)


