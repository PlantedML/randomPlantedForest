library(inline)
library(Rcpp)

rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'

source(paste(rpf_path, '_/generate_data.R', sep=''))
source(paste(rpf_path, '_/rpf.R', sep=''))
source(paste(rpf_path, '_/predict_rpf.R', sep=''))
source(paste(rpf_path, '/purify_rpf.R', sep=''))
sourceCpp(paste(rpf_path, '_/C-Code.cpp', sep=''))
sourceCpp(paste(rpf_path, '/randomPlantedForest-R.cpp', sep=''))

# generate test data
data <- generate_data(Model=7)

# extract samples
samples_x <- data$X[1:10, ]
samples_y <- data$Y_true[1:10]

# create new class object
cat("\014")
rpf_cpp <- new(RandomPlantedForest, data$Y_true, data$X, 1, 1, 50, 0.4)
rpf_R <- rpf(data$Y_true, data$X, max_interaction = 1, ntrees = 1)
#rpf_purify(rpf_R)
#rpf_cpp$purify()
rpf_cpp$print()

# predict
predictions_cpp <- rpf_cpp$predict(samples_x, c(0))
predictions_R <- predict_rpf(samples_x, rpf_R)

# accuracy
MSE <- mean(data$Y_true^2)
MSE_cpp <- sum((samples_y - predictions_cpp)^2) / length(samples_y)
MSE_R <- sum((samples_y - predictions_R)^2) / length(samples_y)
MSE_cpp
MSE_R

