# Es fehlt noch: 
# Alles nach Rename
# Alle CV-Programme

library(data.table)
library(batchtools)
library(ggplot2)

set.seed(1042)

repls <- 10#100 # Bug in batchtools adds repls^2 jobs 10 -> 100 (https://github.com/mllg/batchtools/issues/226)

# Registry ----------------------------------------------------------------
reg_name <- "rpf_sim_all"
reg_dir <- file.path("registries", reg_name)
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
                                  "mhdata_reg_Additive.R",
                                  "planted_forest_CV.R",
                                  "BART_CV.R",
                                  "xgboost_CV.R"
                       ))

# Problems -----------------------------------------------------------
myprob <- function(job, data, ...) {
  
  dat_train <- generate_data(rho=0.3,sparsity=2,sigma=1,  covariates='normal', ...)
  dat_test <- generate_data(rho=0.3,sparsity=2,sigma=1, covariates='normal', ...)
  
  list(train = dat_train, 
       test = dat_test)
}
addProblem(name = "myprob", fun = myprob, seed = 1043)

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
  mse
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
  
  pred<-SBF.reg.LL(frmla,
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
  fit=0
  
  for(i in 1:p){
    
    eval_points <- sapply(  1:length(instance$test$X[,i]), function(l) which(instance$test$X[l,i]==pred$x.grid[[i]])    )
    fit=fit+pred$f_backfit[[i]][eval_points]
  }
  
  mse=mean((fit-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "sbf", fun = run_sbf)

run_bart <- function(data, job, instance, ...) {
  
  fit <- wbart(x.train=instance$train$X, y.train=instance$train$Y_start, sparse = T)
  pred <- predict(fit, instance$test$X)
  
  Y_rep <- matrix(rep(instance$test$Y_true,dim(pred)[1]), ncol=500, byrow= T)
  
  mse <- mean(colMeans(pred-Y_rep)^2)
  mse
}
addAlgorithm(name = "bart", fun = run_bart)

run_mars <- function(data, job, instance, ...) {
  
  fit <- mars(x=instance$train$X, y=instance$train$Y_start)
  pred <- predict(fit, instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "mars", fun = run_mars)

run_rpf_CV <- function(data, job, instance, ...){
  
  parameters <- rpf_CV(Y=instance$train$Y_start, instance$train$X)

  fit <- rpf(X=instance$train$X,
             Y=instance$train$Y_start,
             variables=NULL,
             min_leaf_size=1, 
             splits=parameters$splits,
             split_try=parameters$split_try, 
             t_try=parameters$t_try, 
             max_interaction=parameters$max_interaction,
             ...)
  pred <- predict_rpf(forest_res = fit, X = instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "rpf_CV", fun = run_rpf_CV)

run_xgboost_CV <- function(data, job, instance, ...){
  
  parameters <- xgboost_CV(Y=instance$train$Y_start, X=instance$train$X)
  
  p <- dim(instance$train$X)[2]
  
  train.data <- data.frame(cbind(instance$train$Y_start,instance$train$X))
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  
  test.data <- data.frame(instance$test$X)
  names(test.data) <- paste0("V",1:p)
  
  pred<- xgboost(data = as.matrix(train.data[,2:(p+1)]),
                 label = as.vector(train.data[,"Y"]),
                 max.depth = parameters$interaction_depth,
                 eta = parameters$shrinkage,
                 nthread = 1,
                 nrounds = parameters$n.trees,
                 early_stopping_rounds = NULL,
                 verbose = F,
                 objective = "reg:squarederror")
  fits=predict(pred, as.matrix(test.data))
  mse=mean((fits-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "xgboost_CV", fun = run_xgboost_CV)

run_bart_CV <- function(data, job, instance, ...){
  
  parameters <- BART_CV(Y=instance$train$Y_start, X=instance$train$X)
  
  fit <- wbart(x.train=instance$train$X, y.train=instance$train$Y_start, power=parameters$n_pow, a=parameters$n_spars, ntree=parameters$n_tree, ndpost = 1600, sparse = T)
  pred <- predict(fit, instance$test$X)
  
  Y_rep <- matrix(rep(instance$test$Y_true,dim(pred)[1]), ncol=500, byrow= T)
  
  mse <- mean(colMeans(pred-Y_rep)^2)
  mse
}
addAlgorithm(name = "bart_CV", fun = run_bart_CV)

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
# Generate experiment design from best tuning parameters
best_params <- readRDS("best_params.Rds")
best_params <- best_params[order(algorithm, max_interaction, max.depth, Model, n, p, sparse), ]

# Problems
prob_design <- list(myprob = best_params[algorithm == "rf", .(n, p, Model,sparse)])

# Algorithms
algo_design <- list(
  rf = best_params[algorithm == "rf", .(ntree, mtry, maxnodes)], 
  xgboost = best_params[algorithm == "xgboost", .(eta, nrounds, max.depth)], 
  rpf = best_params[algorithm == "rpf", .(ntrees, splits, split_try, t_try, max_interaction)],
  ranger = best_params[algorithm == "ranger", .(num.trees, mtry, max.depth, replace)],
  gam = data.frame(),
  sbf = best_params[algorithm == "sbf", .(bandwidth)],
  bart = best_params[algorithm == "bart", .(power, a, ntree, ndpost)],
  mars = best_params[algorithm == "mars", .(degree, penalty)],
  rpf_CV = data.frame(),
  xgboost_CV = data.frame(),
  bart_CV = data.frame(),
  average = data.frame(),
  nearestneighbor = data.frame())

addExperiments(prob_design, algo_design, repls = repls, combine = "bind")
summarizeExperiments()

# Test jobs -----------------------------------------------------------
#testJob(id = 1)
#testJob(id = 500)
#testJob(id=2281)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotDone()
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
saveRDS(res, "sim_result_all.Rds")

# Rename
res[, Method := factor(paste(algorithm, max.depth, max_interaction), 
                       levels = c("rf NA NA",  "xgboost 1 NA", "xgboost 2 NA", "xgboost 3 NA", "xgboost 4 NA", "rpf NA 1", "rpf NA 2", "rpf NA 3", "rpf NA 4", "rpf NA 30", "ranger 1 NA", "ranger 2 NA", "ranger 3 NA", "ranger 4 NA", "ranger 0 NA", "gam NA NA", "bart NA NA", "sbf NA NA", "mars NA NA", "rpf_CV NA NA", "xgboost_CV NA NA", "bart_CV NA NA","average NA NA", "nearestneighbor NA NA"), 
                       labels = c("RF", "xgboost additive", "xgboost interaction 2", "xgboost interaction 3",  "xgboost interaction 4", "RPF additive", "RPF interaction 2", "RPF interaction 3", "RPF interaction 4", "RPF interaction 30", "ranger additive", "ranger interaction 2", "ranger interaction 3", "ranger interaction 4", "ranger", "gam", "bart", "sbf", "mars", "rpf_CV", "xgboost_CV", "bart_CV", "average", "nearestneighbor"))]

# Average over repls
res_mean <- res[, mean(mse), by = .(problem, algorithm, Model, sparse, Method, n, p, ntree, mtry, maxnodes, eta, nrounds, max.depth, ntrees, splits, split_try, t_try, max_interaction, num.trees, replace, bandwidth, power, a, ntree, ndpost, degree, penalty)]
res_mean[, mse := V1]
res_mean
res_mean[Model == 2 & p == 4 & sparse, .(Method, round(V1,3))]

# Plot results -------------------------------------------------------------
ggplot(res[sparse == TRUE, ], aes(x = Method, y = mse)) +
  geom_boxplot() + 
  coord_flip() + 
  facet_grid(p ~ Model, scales = "free") + 
  ylab("MSE")
ggsave("sim_sparse.pdf", width = 15, height = 10)

ggplot(res[sparse == FALSE, ], aes(x = Method, y = mse)) +
  geom_boxplot() + 
  coord_flip() + 
  facet_grid(p ~ Model, scales = "free") + 
  ylab("MSE")
ggsave("sim_nosparse.pdf", width = 15, height = 10)

# Export results -------------------------------------------------------------
res_export <- copy(res)
res_export[, problem := NULL]
res_export[, result.1 := NULL]
saveRDS(res_export, "sim_export.Rds")
