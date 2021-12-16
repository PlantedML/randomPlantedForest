library(data.table)
library(batchtools)
library(ggplot2)
setwd("~/Sync/Projects/Random Planted Forest/111-neue_Version/Marvins_Version")

set.seed(42)

repls <- 1

# Data
n <- 500
sparse_p <- c(4)
p <- c(4)
Model <- 1

# RF
ntree <- 500
mtry <- c(1/4) # will be multiplied by p
maxnodes <- c(40,50)

# xgboost
eta <- c(0.005)
nrounds <- c(100)
max.depth <- c(1,2)

# RPF
ntrees <- 50
splits <-  c(10)
split_try <- c(2)
t_try <- c(0.25)
max_interaction <- c(1,2)

# ranger
num.trees <- 500
mtry_rg <- c(1/4) # will be multiplied by p
max.depth_rg <- c(1)
replace <- c(TRUE, FALSE)

# Backfitting

bandwidth <- c(0.1, 0.15)


#BART

bart_pow <- 1
bart_a <- c(0.6)
bart_ntree <- c(50)
bart_ndpost <- 1600

#MARS

mars_degree <- 1:1
mars_penalty <- 1:2

# Registry ----------------------------------------------------------------
reg_name <- "rpf_tune_all"
reg_dir <- file.path("test2")
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = c("randomForest",
                                    "xgboost",
                                    "BART",
                                    "interpret",
                                    "survival",
                                    "randomForest",
                                    "nlme",
                                    "mgcv",
                                    "mda",
                                    "ranger"),
                       source = c("rpf.R", 
                                  "predict_rpf.R", 
                                  "generate_data.R",
                                  "SBF_reg_Additive.R",
                                  "mhdata_reg_Additive.R"
                       ))

# Problems -----------------------------------------------------------
myprob <- function(job, data, ...) {
  dat_train <- generate_data(rho=0.3,sparsity=2,sigma=1, covariates='normal', ...)
  dat_test <- generate_data(rho=0.3,sparsity=2,sigma=1, covariates='normal', ...)
  
  list(train = dat_train, 
       test = dat_test)
}
addProblem(name = "myprob", fun = myprob, seed = 43)

# Algorithms -----------------------------------------------------------
run_rf <- function(data, job, instance, mtry, ...) {
  fit <- randomForest(x = instance$train$X,
                      y = instance$train$Y_start,
                      mtry = mtry * ncol(instance$train$X), 
                      ...)
  pred <- predict(fit, instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "rf", fun = run_rf)

run_xgboost <- function(data, job, instance, ...) {
  fit <- xgboost(data = instance$train$X,
                 label = instance$train$Y_start,
                 nthread = 1,
                 early_stopping_rounds = NULL,
                 verbose = F,
                 objective = "reg:squarederror", 
                 ...)
  pred <- predict(fit, instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "xgboost", fun = run_xgboost)

run_rpf <- function(data, job, instance, ...) {
  fit <- rpf(X=instance$train$X,
             Y=instance$train$Y_start,
             variables=NULL,
             min_leaf_size=1, 
             ...)
  pred <- predict_rpf(forest_res = fit, X = instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "rpf", fun = run_rpf)

run_ranger <- function(data, job, instance, mtry, ...) {
  colnames(instance$train$X) <- paste0("X", 1:ncol(instance$train$X))
  colnames(instance$test$X) <- paste0("X", 1:ncol(instance$test$X))
  fit <- ranger(x = instance$train$X,
                y = instance$train$Y_start,
                mtry = mtry * ncol(instance$train$X), 
                ...)
  pred <- predict(fit, instance$test$X)$predictions
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "ranger", fun = run_ranger)


run_gam <- function(data, job, instance, ...){
  
  p <- dim(instance$train$X)[2]
  
  train.data <- data.frame(cbind(instance$train$Y_start,instance$train$X))
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  
  test.data <- data.frame(instance$test$X)
  names(test.data) <- paste0("V",1:p)
  
  pred <- paste0("s(","V", 1:p,")")
  frmla <- reformulate(pred, "Y")
  
  fit <- gam(frmla, data=train.data,select=TRUE,method="REML")
  pred <- predict(fit, test.data)
  mse=mean((pred-instance$test$Y_true)^2)
}
addAlgorithm(name = "gam", fun = run_gam)

run_sbf <- function(data, job, instance, ...){
  
  p <- dim(instance$train$X)[2]
  pred <- paste0("V", 1:p)
  frmla <- reformulate(pred,"Y")
  
  train.data <- data.frame(cbind(instance$train$Y_start,instance$train$X))
  
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  test.data <- data.frame(instance$test$X)
  
  x.grid <- lapply(1:p, function(k) sort(unique(c(instance$test$X[,k],instance$train$X[,k]))    ))
  
  model_backfit.LL<-SBF.reg.LL(frmla,
                               train.data,
                               weight='sw',
                               it=15,
                               x.grid=x.grid,
                               integral.approx='midd',
                               kcorr=TRUE,
                               LC=FALSE,
                               wrong=FALSE,
                               classic.backfit=FALSE, 
                               ...)
  
  fits_backfit.LL=0
  
  for(i in 1:p){
    
    eval_points <- sapply(  1:length(instance$test$X[,i]), function(l) which(instance$test$X[l,i]==model_backfit.LL$x.grid[[i]])    )
    
    fits_backfit.LL=fits_backfit.LL+model_backfit.LL$f_backfit[[i]][eval_points]

  }

  mse=mean((fits_backfit.LL-instance$test$Y_true)^2)
}
addAlgorithm(name = "sbf", fun = run_sbf)


run_bart <- function(data, job, instance, ...) {
  
  fit <- wbart(x.train=instance$train$X, y.train=instance$train$Y_start, power=bart_pow, a=bart_a, ntree=bart_ntree, ndpost = bart_ndpost, sparse = T)
  
  pred <- predict(fit, instance$test$X)
  
  Y_rep <- matrix(rep(instance$test$Y_true,dim(pred)[1]), ncol=500, byrow= T)
  
  mse <- mean(colMeans(pred-Y_rep)^2)
  mse
}
addAlgorithm(name = "bart", fun = run_bart)

run_mars <- function(data, job, instance, ...) {
  
  fit <- mars(x=instance$train$X, y=instance$train$Y_start, degree=mars_degree, penalty = mars_penalty)
  
  pred <- predict(fit, instance$test$X)
  
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "mars", fun = run_mars)

run_average <- function(data, job, instance, ...){
  
  mse <- mean((mean(instance$train$Y_start)-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "average", fun = run_average)

run_nearestneighbor <- function(data, job, instance, ...){
  
  fit <- knn(instance$train$X, instance$test$X, instance$train$Y_start)
  pred <- as.numeric(paste(fit))
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "nearestneighbor", fun = run_nearestneighbor)

# Experiments -----------------------------------------------------------
prob_design <- list(myprob = rbind(expand.grid(n = n, 
                                               p = sparse_p, 
                                               Model = Model,
                                               sparse = T,
                                               stringsAsFactors = FALSE),
                                   expand.grid(n = n, 
                                               p = p, 
                                               Model = Model,
                                               sparse = F,
                                               stringsAsFactors = FALSE)))
{
algo_design <- list(rf = expand.grid(ntree = ntree, 
                                     mtry = mtry, 
                                     maxnodes = maxnodes,
                                     stringsAsFactors = FALSE), 
                    xgboost = expand.grid(eta = eta,
                                          nrounds = nrounds, 
                                          max.depth = max.depth, 
                                          stringsAsFactors = FALSE), 
                     rpf = expand.grid(ntrees = ntrees, 
                                       splits = splits,
                                       split_try = split_try, 
                                       t_try = t_try,
                                       max_interaction = max_interaction,
                                       stringsAsFactors = FALSE), 
                    ranger = expand.grid(num.trees = num.trees, 
                                         mtry = mtry_rg, 
                                         max.depth = max.depth_rg,
                                         replace = replace,
                                         stringsAsFactors = FALSE),
                     gam = data.frame(),
                     sbf = data.frame(bandwidth=bandwidth),
                     # bart = data.frame(bart_pow=bart_pow,
                     #                   bart_a=bart_a,
                     #                   bart_ntree=bart_ntree,
                     #                   bart_ndpost=bart_ndpost),
                     mars = data.frame(mars_degree=mars_degree,
                                       mars_penalty=mars_penalty),
                     average = data.frame(),
                     nearestneighbor = data.frame()
                    )
}
addExperiments(prob_design, algo_design, repls = repls)
summarizeExperiments()

# Test jobs -----------------------------------------------------------
#testJob(id = 1)
#testJob(id = 500)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = reg_name, chunks.as.arrayjobs = TRUE, 
                              ncpus = 1, memory = 6000, walltime = 10*24*3600, 
                              max.concurrent.jobs = 400))
} else {
  submitJobs()
}
waitForJobs()

# Get results -------------------------------------------------------------
res <-  flatten(ijoin(reduceResultsDataTable(), getJobPars()))
res[, mse := result.1]

# Save result
saveRDS(res, "tune_result_all.Rds")

# Average over repls
res_mean <- res[, mean(mse), by = .(problem, algorithm, n, p, Model, ntree, mtry, maxnodes, eta, nrounds, max.depth, ntrees, splits, split_try, t_try, max_interaction, num.trees, replace)]
res_mean[, mse := V1]

# Get best parameters per method
best_params <- res_mean[ , .SD[which.min(mse)], by = .(algorithm, n, p, Model, max_interaction, max.depth)]
saveRDS(best_params, "best_params.Rds")

# Results
# RF: 
# xgboost additive (max.depth=2):
# xgboost interaction: 
# RPF additive: 
# RPF interaction: 