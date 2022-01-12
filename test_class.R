# source code ------------------------
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)
source(paste(rpf_path, '/randomPlantedForest.R', sep=''))

# generate test data ------------------------

# generate data for regression
sample_size <- 500
data <- generate_data(Model=1, p=4, n=sample_size)
test_size <- floor(length(data$Y_true) / 5)
x_test <- data$X[1:test_size, ] # extract samples
y_test <- data$Y_true[1:test_size]
x_train <- data$X[(test_size+1):sample_size, ] # extract samples
y_train <- data$Y_start[(test_size+1):sample_size]

# load iris as classification dataset
library(datasets)

# prepare dataset
data_class <- as.list(subset(iris, Species != "virginica")) # exclude third class 
data_class$Y[data_class$Species=='setosa'] <- 1 # express class with value
data_class$Y[data_class$Species=='versicolor'] <- 0
data_class$X <- cbind(data_class$Sepal.Length, data_class$Sepal.Width, data_class$Petal.Length, data_class$Petal.Width) # save X as matrix
data_class[c('Species', 'Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width')] <- NULL
sample_size_2 <- length(data_class$Y)
test_size_2 <- floor(sample_size_2 / 5)

# shuffle data
order <- sample(1:sample_size_2)
data_class$X <- data_class$X[order, ]
data_class$Y <- data_class$Y[order]

# extract samples
x_test_class <- data_class$X[1:test_size_2, ] 
y_test_class <- data_class$Y[1:test_size_2]
x_train_class <- data_class$X[(test_size_2+1):sample_size_2, ]
y_train_class <- data_class$Y[(test_size_2+1):sample_size_2]

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

start_time <- Sys.time()
rpf_cpp_class <- new_rpf(y_train_class, x_train_class, loss=loss, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits = n_splits, split_try=split_try, deterministic=deterministic_forest, parallel=F)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
rpf_R_class <- rpf(y_train_class, x_train_class, loss=loss, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try=split_try, deterministic=deterministic_forest)
end_time <- Sys.time()
end_time - start_time

predictions_cpp_class <- rpf_cpp_class$predict_matrix(x_test_class, c(0))
predictions_R_class <- predict_rpf(x_test_class, rpf_R_class, c(0))

variation_class <- mean(y_test_class^2)
MSE_cpp_class <- rpf_cpp_class$MSE(sigma(predictions_cpp_class), y_test_class) 
MSE_R_class <- sum((y_test_class - sigma(predictions_R_class))^2) / length(y_test_class)
MSE_cpp_class
MSE_R_class
variation_class