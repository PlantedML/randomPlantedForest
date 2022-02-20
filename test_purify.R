
# source code ------------------------
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)
source('randomPlantedForest.R')

# generate test data ------------------------

sample_size <- 500
data <- generate_data(Model=2, p=4, n=sample_size)
test_size <- floor(length(data$Y_true) / 5)
x_test <- data$X[1:test_size, ] # extract samples
y_test <- data$Y_true[1:test_size]
x_train <- data$X[(test_size+1):sample_size, ] # extract samples
y_train <- data$Y_start[(test_size+1):sample_size]


# set parameters ------------------------
n_splits <- 15
max_inter <- 4
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

#testRPF()

start_time <- Sys.time()
rpf_cpp <- new_rpf(y_train, x_train,  max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try = split_try, deterministic=deterministic_forest, parallel=parallel)
end_time <- Sys.time()
end_time - start_time
rpf_cpp_model_old <- rpf_cpp$get_model()

rpf_cpp$new_purify()
rpf_cpp_model_new <- rpf_cpp$get_model()

start_time <- Sys.time()
rpf_R <- rpf(y_train, x_train, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try=split_try, deterministic=deterministic_forest)
end_time <- Sys.time()
end_time - start_time


# purify ------------------------
rpf_cpp$purify()
rpf_cpp$new_purify()
rpf_R_purified <- rpf_purify(rpf_R)


# predict ------------------------
predictions_cpp <- predict(rpf_cpp, x_test, c(0))
predictions_R <- predict_rpf(x_test, rpf_R, c(0))


# accuracy ------------------------
variation <- mean(y_test^2)
MSE_cpp <- rpf_cpp$MSE(predictions_cpp, y_test) 
MSE_R <- sum((y_test - predictions_R)^2) / length(y_test)
MSE_cpp
MSE_R
variation



