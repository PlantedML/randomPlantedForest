source("EBM-Regression.R")
source("SBF_reg_Additive.R")
source("mhdata_reg_Additive.R")
source("generate_data.R")
source("planted_forest_main_simu.R")
source("planted_forest_CV.R")

library(parallel)

Monte_Carlo=100
Kerne=20

simulate_function <- function(Monte_Carlo){
  
  library(Rcpp)
  sourceCpp("C-Code.cpp")
  library(interpret)
  library(survival)
  library(randomForest)
  library(nlme)
  library(mgcv)
  library(xgboost)
  
  # Generate Data
  
  Data<-generate_data(n=500,p=10,rho=0.3,sparsity=3,sigma=1, Model=4, covariates='normal')
  
  p <- dim(Data$X)[2]
  
  train.data <- data.frame(cbind(Data$Y_start,Data$X))
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  
  pred <- paste0("V", 1:p)
  frmla <- reformulate(pred,"Y")
  
  #Estimator including interaction terms
  
  model_xgboost_interaction<- xgboost(data = as.matrix(train.data[,2:(p+1)]),
                                      label = as.vector(train.data[,"Y"]),
                                      max.depth = 2,
                                      eta = 0.08,
                                      nthread = 1,
                                      nrounds = 600,
                                      early_stopping_rounds = NULL,
                                      verbose = F,
                                      objective = "reg:squarederror")
  
  fits=predict(model_xgboost_interaction, as.matrix(train.data[,2:(p+1)]))
  
  MSE_xgboost_interaction=mean((fits-Data$Y_true)^2)
  
  #EBM
  
  # Additive Estimator
  
  model_EBM_additive <- ebm_regression(data.frame(Data$X),
                                       Data$Y_start,
                                       num_outer_bags = 16,
                                       validation_size = 0.15,
                                       max_epochs = 3000,
                                       num_early_stopping_run_length = 50,
                                       learning_rate = 0.04,
                                       max_tree_splits = 2,
                                       min_instances_for_split = 2,
                                       random_state = sample(1:999999,1))
  
  fits=ebm_prediction(model_EBM_additive, data.frame(Data$X))
  
  MSE_EBM_additive=mean((fits-Data$Y_true)^2)
  
  model_randomForests <- randomForest(x=Data$X,
                                      y=Data$Y_start,
                                      ntree = 500,
                                      m_try = 7,
                                      nodesize = 3,
                                      maxnodes = 500)
  
  fits=predict(model_randomForests)
  
  MSE_randomForests=mean((fits-Data$Y_true)^2)
  
  #PF
  
  fits=planted_forest_simu(Y=Data$Y_start,
                        X=Data$X,
                        max_interaction=1,
                        ntrees=50,
                        splits=125,
                        split_try=20,
                        t_try=0.25,
                        variables=NULL,
                        leaf_size=rep(1,p),
                        alternative=T)
  
  MSE_PF_alternative=mean((fits-Data$Y_true)^2)
  
  #Additive estimator
  
  fits=planted_forest_simu(Y=Data$Y_start,
                        X=Data$X,
                        max_interaction=1,
                        ntrees=50,
                        splits=125,
                        split_try=20,
                        t_try=0.25,
                        variables=NULL,
                        leaf_size=rep(1,p))
  
  MSE_PF_additive=mean((fits-Data$Y_true)^2)
  
  #Estimator including interaction terms
  
  fits=planted_forest_simu(Y=Data$Y_start,
                        X=Data$X,
                        max_interaction=2,
                        ntrees=50,
                        splits=150,
                        split_try=20,
                        t_try=0.5,
                        variables=NULL,
                        leaf_size=rep(1,p))
  
  MSE_PF_interaction=mean((fits-Data$Y_true)^2)
  
  #Estimator including arbitrary interaction terms
  
  fits=planted_forest_simu(Y=Data$Y_start,
                        X=Data$X,
                        max_interaction=4,
                        ntrees=50,
                        splits=125,
                        split_try=20,
                        t_try=0.75,
                        variables=NULL,
                        leaf_size=rep(1,p))
  
  MSE_PF_fullinteraction=mean((fits-Data$Y_true)^2)
  
  #BF
  
  b.LL.grid <- b.LC.grid <- b.BF.grid <- numeric(p)
  
  b.LL.grid[1:(p+1)] <- rep(c(0.1),p+1)
  
  model_backfit.LL<-SBF.reg.LL(frmla,
                               train.data,
                               b.LL.grid,
                               weight='sw',
                               it=15,
                               x.grid=NULL,
                               integral.approx='midd',
                               kcorr=TRUE,
                               LC=FALSE,
                               wrong=FALSE,
                               classic.backfit=FALSE)
  
  b.LC.grid[1:(p+1)] <- rep(c(0.1),p+1)  
  
  model_backfit.LC<-SBF.reg.LL(frmla,
                               train.data,
                               b.LC.grid,
                               weight='sw',
                               it=15,
                               x.grid=NULL,
                               integral.approx='midd',
                               kcorr=TRUE,
                               LC=TRUE,
                               wrong=FALSE,
                               classic.backfit=FALSE)
  
  b.BF.grid[1:(p+1)] <- rep(c(0.1),p+1)
  
  model_backfit.BF<-SBF.reg.LL(frmla,
                               train.data,
                               b.BF.grid,
                               weight='sw',
                               it=15,
                               x.grid=NULL,
                               integral.approx='midd',
                               kcorr=TRUE,
                               LC=TRUE,
                               wrong=FALSE,
                               classic.backfit=TRUE)
  
  fits_backfit.LL=0
  fits_backfit.LC=0
  fits_backfit.BF=0
  
  for(i in 1:p){
    
    names(model_backfit.LL$f_backfit[[i]])=names(model_backfit.LL$x.grid[[i]])
    names(model_backfit.LC$f_backfit[[i]])=names(model_backfit.LC$x.grid[[i]])
    names(model_backfit.BF$f_backfit[[i]])=names(model_backfit.BF$x.grid[[i]])
    
    fits_backfit.LL=fits_backfit.LL+model_backfit.LL$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.LL$x.grid[[i]]))))]
    fits_backfit.LC=fits_backfit.LC+model_backfit.LC$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.LC$x.grid[[i]]))))]
    fits_backfit.BF=fits_backfit.BF+model_backfit.BF$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.BF$x.grid[[i]]))))]
  }
  
  MSE_backfit.LL=mean((fits_backfit.LL-Data$Y_true)^2)
  MSE_backfit.LC=mean((fits_backfit.LC-Data$Y_true)^2)
  MSE_backfit.BF=mean((fits_backfit.BF-Data$Y_true)^2)
  
  #Pspline
  
  pred2 <- paste0("pspline(","V", 1:p,",df=",c(4),")")
  frmla2 <- reformulate(pred2,"Y")
  model_pspline=lm(frmla2, train.data)
  fits=predict(model_pspline)
  
  MSE_pspline=mean((fits-Data$Y_true)^2)
  
  pred3 <- paste0("s(","V", 1:p,")")
  frmla3 <- reformulate(pred3, "Y")
  model_gam <- gam(frmla3, data=train.data)
  fits <- predict(model_gam, train.data)
  
  if(!is.null(fits)){MSE_gam=mean((fits-Data$Y_true)^2)} else {MSE_gam=1000}
  
  #Gradient boosting
  
  model_xgboost_additive <- xgboost(data = as.matrix(train.data[,2:(p+1)]),
                                    label = as.vector(train.data[,"Y"]),
                                    max.depth = 1,
                                    eta = 0.08,
                                    nthread = 1,
                                    verbose = F,
                                    early_stopping_rounds = NULL,
                                    nrounds = 7000,
                                    objective = "reg:squarederror")
  
  fits=predict(model_xgboost_additive, as.matrix(train.data[,2:(p+1)]))
  
  MSE_xgboost_additive=mean((fits-Data$Y_true)^2)
  
  return(list(X=Data$X,
              Y=Data$Y_start,
              MSE_xgboost_additive=MSE_xgboost_additive,
              MSE_xgboost_interaction=MSE_xgboost_interaction,
              MSE_EBM_additive=MSE_EBM_additive,
              MSE_PF_alternative=MSE_PF_alternative,
              MSE_PF_additive=MSE_PF_additive,
              MSE_PF_interaction=MSE_PF_interaction,
              MSE_PF_fullinteraction=MSE_PF_fullinteraction,
              MSE_backfit.LL=MSE_backfit.LL,
              MSE_backfit.LC=MSE_backfit.LC,
              MSE_backfit.BF=MSE_backfit.BF,
              MSE_pspline=MSE_pspline,
              MSE_gam=MSE_gam,
              MSE_randomForests=MSE_randomForests))
}

cl <- makeCluster(Kerne)

clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
set.seed(23)

Results <- parSapply(cl,1:Monte_Carlo, simulate_function)

stopCluster(cl)

# save(Results,file="Model4-K=10_nosparse_FinalResults.Rdata")