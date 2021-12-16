library(data.table)
library(batchtools)
library(ggplot2)

set.seed(1042)

repls <- 2

n <- 500
p <- 4
sparse <- F
Model <- 1

# RF
ntree <- 500
mtry <- 22
maxnodes <- 500

# xgboost
eta <- c(0.04, 0.02)
nrounds <- c(7000, 7000)
max.depth <- c(2, 3)

# RPF
ntrees <- 50
splits <- 5
split_try <- 2
t_try <- 0.75
max_interaction <- 1

# Backfitting

b.grid <- 0.2

#Pspline

df <- 2

#BART

bart_pow <- 1
bart_a <- 0.6
bart_ntree <- 50
bart_ndpost <- 400

#MARS

mars_degree <- 1
mars_penalty <- 1

# Registry ----------------------------------------------------------------
reg_dir <- file.path("testtest")
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
  dat_train <- generate_data(rho=0.3,sparsity=2,sigma=1,covariates='normal', ...)
  dat_test <- generate_data(rho=0.3,sparsity=2,sigma=1,covariates='normal', ...)
  
  list(train = dat_train, 
       test = dat_test)
}
addProblem(name = "myprob", fun = myprob, seed = 1043)

# Algorithms -----------------------------------------------------------
run_rf <- function(data, job, instance, ...) {
  fit <- randomForest(x = instance$train$X,
                      y = instance$train$Y_true,
                      ...)
  pred <- predict(fit, instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "rf", fun = run_rf)

run_xgboost <- function(data, job, instance, ...) {
  fit <- xgboost(data = instance$train$X,
                 label = instance$train$Y_true,
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
  fit <- rpf(X = instance$train$X,
             Y = instance$train$Y_true,
             variables = NULL,
             min_leaf_size = 1,
             ...)
  pred <- predict_rpf(forest_res = fit, X = instance$test$X)
  mse <- mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "rpf", fun = run_rpf)

# run_backfit.LL <- function(data, job, instance, ...) {
# 
#   p <- dim(instance$train$X)[2]
# 
#   train.data <- data.frame(cbind(instance$train$Y_start,instance$train$X))
#   names(train.data)[1] <- "Y"
#   names(train.data)[2:(p+1)] <- paste0("V", 1:p)
# 
#   pred <- paste0("V", 1:p)
#   frmla <- reformulate(pred,"Y")
# 
#   grid <- numeric(p)
#   grid[1:(p+1)] <- rep(c(b.grid),p+1)
# 
#   model<-SBF.reg.LL(frmla,
#                     train.data,
#                     grid,
#                     weight='sw',
#                     it=15,
#                     x.grid=NULL,
#                     integral.approx='midd',
#                     kcorr=TRUE,
#                     LC=FALSE,
#                     wrong=FALSE,
#                     classic.backfit=FALSE)
# 
#   pred <- 0
# 
#   for(i in 1:p){
# 
#     names(model$f_backfit[[i]])=names(model$x.grid[[i]])
# 
#     pred <- pred+model$f_backfit[[i]][as.character(sort(as.numeric(names(model$x.grid[[i]]))))]
#   }
# 
#   mse=mean((pred-instance$train$Y_true)^2)
#   mse
# }
# addAlgorithm(name = "backfit.LL", fun = run_backfit.LL)

run_pspline <- function(data, job, instance, ...) {
  
  p <- dim(instance$train$X)[2]
  
  train.data <- data.frame(cbind(instance$train$Y_start,instance$train$X))
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  
  test.data <- data.frame(instance$test$X)
  names(test.data) <- paste0("V",1:p)
  
  pred <- paste0("pspline(","V", 1:p,",df=",c(df),")")
  frmla <- reformulate(pred,"Y")
  
  fit <- lm(frmla, train.data)
  pred <- predict(fit, test.data)
  mse=mean((pred-instance$test$Y_true)^2)
  mse
}
addAlgorithm(name = "pspline", fun = run_pspline)

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
prob_design <- list(myprob = expand.grid(n = n, 
                                         p = p, 
                                         sparse = sparse,
                                         Model = Model,
                                         stringsAsFactors = FALSE))
algo_design <- list(rf = data.frame(ntree = ntree, 
                                    mtry = mtry, 
                                    maxnodes = maxnodes,
                                    stringsAsFactors = FALSE), 
                    xgboost = data.frame(eta = eta,
                                         nrounds = nrounds, 
                                         max.depth = max.depth, 
                                         stringsAsFactors = FALSE), 
                    rpf = data.frame(ntrees = ntrees,
                                     splits = splits,
                                     split_try = split_try,
                                     t_try = t_try,
                                     max_interaction = max_interaction,
                                     stringsAsFactors = FALSE),
                    pspline = data.frame(df = df),
                    gam = data.frame(),
                    bart = data.frame(bart_pow=bart_pow,
                                      bart_a=bart_a,
                                      bart_ntree=bart_ntree,
                                      bart_ndpost=bart_ndpost),
                    mars = data.frame(mars_degree=mars_degree,
                                      mars_penalty=mars_penalty),
                    average = data.frame(),
                    nearestneighbor = data.frame())
addExperiments(prob_design, algo_design, repls = repls)
summarizeExperiments()

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
saveRDS(res, "sim_result.Rds")

# Average over repls
res_mean <- res[, mean(mse), by = .(problem, algorithm, n, p, ntree, mtry, maxnodes, eta, nrounds, max.depth, ntrees, splits, split_try, t_try, max_interaction, df, bart_pow, bart_a, bart_ntree, bart_ndpost)]
res_mean[, mse := V1]
res_mean