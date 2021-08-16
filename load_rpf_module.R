library(inline)
library(Rcpp)

rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'

source(paste(rpf_path, '_/generate_data.R', sep=''))

# generate test data
data <- generate_data(Model=1)

# extract samples
samples_x <- data$X[1:10, ]
samples_y <- data$Y_true[1:10]

source(paste(rpf_path, '_/predict_rpf.R', sep=''))
source(paste(rpf_path, '/purify_rpf.R', sep=''))
source(paste(rpf_path, '_/rpf.R', sep=''))
sourceCpp(paste(rpf_path, '_/C-Code.cpp', sep=''))
sourceCpp(paste(rpf_path, '/randomPlantedForest-R.cpp', sep=''))

# create new class object
cat("\014")
n_splits <- 20;
max_inter <- 2;
n_trees <- 1;
rpf_cpp <- new(RandomPlantedForest, data$Y_start, data$X, max_inter, n_trees, n_splits, 0.4)
rpf_R <- rpf(data$Y_start, data$X, max_interaction = max_inter, ntrees = n_trees, splits = n_splits, deterministic=TRUE)
#rpf_purify(rpf_R)
#rpf_cpp$purify()
#rpf_cpp$print()

# predict
predictions_cpp <- rpf_cpp$predict(samples_x, c(0))
predictions_R <- predict_rpf(samples_x, rpf_R)

# accuracy
MSE <- mean(data$Y_true^2)
MSE_cpp <- sum((samples_y - predictions_cpp)^2) / length(samples_y)
MSE_R <- sum((samples_y - predictions_R)^2) / length(samples_y)
MSE_cpp
MSE_R

