
# source code ------------------------
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)
source(paste(rpf_path, '/randomPlantedForest.R', sep=''))

# generate test data ------------------------

sample_size <- 500
data <- generate_data(Model=1, p=4, n=sample_size)
test_size <- floor(length(data$Y_true) / 5)
x_test <- data$X[1:test_size, ] # extract samples
y_test <- data$Y_true[1:test_size]
x_train <- data$X[(test_size+1):sample_size, ] # extract samples
y_train <- data$Y_start[(test_size+1):sample_size]


# set parameters ------------------------
n_splits <- 15
max_inter <- 2
n_trees <- 50
split_try <- 10
t_try <- 0.5
deterministic_forest <- TRUE
parallel <- TRUE
purify_forest <- FALSE
loss <- 'logit' 
delta <- 0.1
epsilon <- 0

# train models ------------------------

start_time <- Sys.time()
rpf_cpp <- new_rpf(y_train, x_train,  max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try = split_try, deterministic=deterministic_forest, parallel=parallel)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
rpf_R <- rpf(y_train, x_train, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try=split_try, deterministic=deterministic_forest)
end_time <- Sys.time()
end_time - start_time

# change parameters ------------------------
rpf_cpp$set_parameters("deterministic", 1)
rpf_cpp$set_parameters(c("t_try", "split_try"), c(0.4, 5))

# cross-validation ------------------------
rpf_cpp$cross_validation(5, c(5, 20), c(0.2,0.5,0.7,0.9), c(1,2,5,10))

# purify ------------------------
rpf_cpp$purify()
rpf_cpp$new_purify()
rpf_R_purified <- rpf_purify(rpf_R)

# view result ------------------------
rpf_cpp$print()
rpf_cpp$get_parameters()
rpf_cpp_model <- rpf_cpp$get_model()
rpf_cpp_model
rpf_R

# predict ------------------------

predictions_cpp <- predict(rpf_cpp, x_test[, 1:2], c(1,2))
predictions_R <- predict_rpf(x_test[, 1:2], rpf_R, c(1,2))

predictions_cpp <- predict(rpf_cpp, x_test[, 1], c(1))
predictions_R <- predict_rpf(x_test[, 1], rpf_R, c(1))

predictions_cpp <- predict(rpf_cpp, x_test[1, ], c(0))
predictions_R <- predict_rpf(x_test[1, ], rpf_R, c(0))

predictions_cpp <- predict(rpf_cpp, x_test, c(0))
predictions_R <- predict_rpf(x_test, rpf_R, c(0))

predictions_cpp <- predict(rpf_cpp, x_test, c(0))
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

timings <- benchmark( "rpf_cpp_sequential" = {
                rpf_cpp <- new_rpf(y_train, x_train,  max_interaction=max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, split_try = split_try, deterministic=deterministic_forest, parallel = FALSE)
            },
            "rpf_cpp_parallel" = {
                rpf_cpp <- new_rpf(y_train, x_train,  max_interaction=max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, split_try = split_try, deterministic=deterministic_forest, parallel = TRUE)
            },
           "rpf_cpp_predict" = {
                predictions_cpp <- rpf_cpp$predict_matrix(x_test, c(0))
           },
           "rpf_R" = {
                rpf_R <- rpf(y_train, x_train, max_interaction=max_inter, t_try=t_try, ntrees = n_trees, splits = n_splits, split_try = split_try, deterministic=deterministic_forest)
           },
           "rpf_R_predict" = {
                predictions_R <- predict_rpf(x_test, rpf_R, c(0))
           },
           replications=10
)

benchmark_rpf()

