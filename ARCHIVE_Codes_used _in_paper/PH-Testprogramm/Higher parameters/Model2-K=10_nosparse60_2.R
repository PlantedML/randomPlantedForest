source("PH-Testprogramm/EBM-Regression.R")
source("PH-Testprogramm/SBF_reg_Additive.R")
source("PH-Testprogramm/mhdata_reg_Additive.R")
source("PH-Testprogramm/generate_data.R")
# source("EBM-Regression.R")
# source("SBF_reg_Additive.R")
# source("mhdata_reg_Additive.R")
# source("generate_data.R")
library(parallel)

Monte_Carlo=30
Kerne=Monte_Carlo

simulate_function <- function(Monte_Carlo){
  
  TimeAll=Sys.time()
  
  library(Rcpp)
  sourceCpp("PH-Testprogramm/C-Code_Number_2.cpp")
  # sourceCpp("C-Code_Number_2.cpp")
  library(interpret)
  library(survival)
  library(randomForest)
  library(nlme)
  library(mgcv)
  library(parallel)
  library(xgboost)
  
  # Generate Data
  
  Data<-generate_data(n=500,p=10,rho=0.3,sparcity=9,sigma=1, Model=2, covariates='normal')
  
  p <- dim(Data$X)[2]
  
  train.data <- data.frame(cbind(Data$Y_start,Data$X))
  names(train.data)[1] <- "Y"
  names(train.data)[2:(p+1)] <- paste0("V", 1:p)
  
  pred <- paste0("V", 1:p)
  frmla <- reformulate(pred,"Y")
  
  # All required Parameters
  
  # Random Forests
  
  m_try=c(floor(p/4), floor(p/2), floor(3*p/4))
  maxnodes=c(40,60,80,100,120,dim(Data$X)[1])
  
  # Random Mixed Forests
  
  splits=(p/2)*c(60)
  m=c(5,20)
  Itersplit_try=c(0.1,0.25,0.75)
  
  # Backfitting
  
  b.grid=c(0.1,0.2,0.3)
  
  # Pspline
  
  df=c(2,3,4)
  
  # xgBoost/Gradient Boosting
  
  shrinkage=c(0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
  n.trees=c(100, 300, 600, 1000 , 3000, 5000, 7000)
  interaction_depth=c(2,3,4)
  
  #Random Forests
  
  #Estimator including interaction terms
  
  # TimeGradientBoosting_interaction=Sys.time()
  # 
  # MSE_xgboost_interactionES<-MSE_xgboost_interaction<-MSE_EBM_additive<-MSE_EBM_additive2<-rep(0,length(shrinkage)*length(n.trees)*length(interaction_depth))
  # dim(MSE_xgboost_interactionES)<-dim(MSE_xgboost_interaction)<-dim(MSE_EBM_additive)<-dim(MSE_EBM_additive2)<-c(length(shrinkage),length(n.trees),length(interaction_depth))
  # dimnames(MSE_xgboost_interactionES)<-dimnames(MSE_xgboost_interaction)<-dimnames(MSE_EBM_additive)<-dimnames(MSE_EBM_additive2)<-list(shrinkage,n.trees,interaction_depth)
  # 
  # for(gbm_shrinkage in shrinkage){
  #   for(gbm_n.trees in n.trees){
  #     for(gbm_interaction_depth in interaction_depth){
  #       
  #       model<- xgboost(data = as.matrix(train.data[,2:(p+1)]),
  #                       label = as.vector(train.data[,"Y"]),
  #                       max.depth = gbm_interaction_depth,
  #                       eta = gbm_shrinkage,
  #                       nthread = 1,
  #                       nrounds = gbm_n.trees,
  #                       early_stopping_rounds = 50,
  #                       verbose = F,
  #                       objective = "reg:squarederror")
  #       
  #       fits=predict(model, as.matrix(train.data[,2:(p+1)]))
  #       
  #       MSE_xgboost_interactionES[as.character(gbm_shrinkage),as.character(gbm_n.trees),as.character(gbm_interaction_depth)]=mean((fits-Data$Y_true)^2)
  #       
  #       model<- xgboost(data = as.matrix(train.data[,2:(p+1)]),
  #                       label = as.vector(train.data[,"Y"]),
  #                       max.depth = gbm_interaction_depth,
  #                       eta = gbm_shrinkage,
  #                       nthread = 1,
  #                       nrounds = gbm_n.trees,
  #                       early_stopping_rounds = NULL,
  #                       verbose = F,
  #                       objective = "reg:squarederror")
  #       
  #       fits=predict(model, as.matrix(train.data[,2:(p+1)]))
  #       
  #       MSE_xgboost_interaction[as.character(gbm_shrinkage),as.character(gbm_n.trees),as.character(gbm_interaction_depth)]=mean((fits-Data$Y_true)^2)
  #       
  #       #EBM
  #       
  #       # Additive Estimator
  #       
  #       model <- ebm_regression(data.frame(Data$X),
  #                               Data$Y_start,
  #                               num_outer_bags = 16,
  #                               validation_size = 0.15,
  #                               max_epochs = gbm_n.trees,
  #                               num_early_stopping_run_length = 50,
  #                               learning_rate = gbm_shrinkage,
  #                               max_tree_splits = gbm_interaction_depth,
  #                               min_instances_for_split = 2,
  #                               random_state = sample(1:999999,1))
  #       
  #       fits=ebm_prediction(model, data.frame(Data$X))
  #       
  #        MSE_EBM_additive[as.character(gbm_shrinkage),as.character(gbm_n.trees),as.character(gbm_interaction_depth)]=mean((fits-Data$Y_true)^2)
  #       # 
  #       # model <- ebm_regression(data.frame(Data$X),
  #       #                         Data$Y_start,
  #       #                         num_outer_bags = 16,
  #       #                         validation_size = 0,
  #       #                         max_epochs = gbm_n.trees,
  #       #                         num_early_stopping_run_length = gbm_n.trees,
  #       #                         learning_rate = gbm_shrinkage,
  #       #                         max_tree_splits = gbm_interaction_depth,
  #       #                         min_instances_for_split = 2,
  #       #                         random_state = sample(1:999999,1))
  #       # 
  #       # fits=ebm_prediction(model, data.frame(Data$X))
  #       # 
  #       # MSE_EBM_additive2[as.character(gbm_shrinkage),as.character(gbm_n.trees),as.character(gbm_interaction_depth)]=mean((fits-Data$Y_true)^2)
  #       # 
  #       
  #     }
  #   }
  # }
  # 
  # TimeGradientBoosting_interaction=Sys.time()-TimeGradientBoosting_interaction
  # 
  # TimeRandomForests<-Sys.time()
  # 
  # MSE_randomForests=rep(0,length(m_try)*length(maxnodes))
  # dim(MSE_randomForests)=c(length(m_try),length(maxnodes))
  # dimnames(MSE_randomForests)=list(m_try,maxnodes)
  # 
  # for(randomForests_m_try in m_try){
  #   for(randomForests_maxnodes in maxnodes){
  #     
  #     model <- randomForest(x=Data$X,
  #                           y=Data$Y_start,
  #                           ntree = 500,
  #                           m_try = randomForests_m_try,
  #                           nodesize = 3,
  #                           maxnodes = randomForests_maxnodes)
  #     
  #     fits=predict(model)
  #     
  #     MSE_randomForests[as.character(randomForests_m_try),as.character(randomForests_maxnodes)]=mean((fits-Data$Y_true)^2)
  #   }
  # }
  # 
  # TimeRandomForests=Sys.time()-TimeRandomForests
  # 
  # #PF
  
  TimePF=Sys.time()
  
  MSE_PF_alternative<-MSE_PF_additive<-MSE_PF_interaction<-MSE_PF_fullinteraction<-rep(0,length(splits)*length(m)*length(Itersplit_try))
  dim(MSE_PF_alternative)<-dim(MSE_PF_additive)<-dim(MSE_PF_interaction)<-dim(MSE_PF_fullinteraction)<-c(length(splits),length(m),length(Itersplit_try))
  dimnames(MSE_PF_alternative)<-dimnames(MSE_PF_additive)<-dimnames(MSE_PF_interaction)<-dimnames(MSE_PF_fullinteraction)<-list(splits,m,Itersplit_try)
  
  for(PF_splits in splits){
    for(PF_m in m){
      for(PF_Itersplit_try in Itersplit_try){
        
        # TimePF_alternative=Sys.time()
        # 
        # fits=planted_forest_2(Y=Data$Y_start,
        #                       X=Data$X,
        #                       max_interaction=1,
        #                       Baum=50,
        #                       splits=PF_splits,
        #                       m=PF_m,
        #                       Itersplit_try=PF_Itersplit_try,
        #                       variables=NULL,
        #                       Blattgroesse=rep(1,p),
        #                       alternative=T)
        # 
        # MSE_PF_alternative[as.character(PF_splits),as.character(PF_m),as.character(PF_Itersplit_try)]=mean((fits-Data$Y_true)^2)
        # 
        # TimePF_alternative=Sys.time()-TimePF_alternative
        # 
        # #Additive estimator
        # 
        # TimePF_additive=Sys.time()
        # 
        # fits=planted_forest_2(Y=Data$Y_start,
        #                       X=Data$X,
        #                       max_interaction=1,
        #                       Baum=50,
        #                       splits=PF_splits,
        #                       m=PF_m,
        #                       Itersplit_try=PF_Itersplit_try,
        #                       variables=NULL,
        #                       Blattgroesse=rep(1,p))
        # 
        # MSE_PF_additive[as.character(PF_splits),as.character(PF_m),as.character(PF_Itersplit_try)]=mean((fits-Data$Y_true)^2)
        # 
        # TimePF_additive=Sys.time()-TimePF_additive
        # 
        # #Estimator including interaction terms
        # 
        TimePF_interaction=Sys.time()

        fits=planted_forest_2(Y=Data$Y_start,
                              X=Data$X,
                              max_interaction=2,
                              Baum=50,
                              splits=PF_splits,
                              m=PF_m,
                              Itersplit_try=PF_Itersplit_try,
                              variables=NULL,
                              Blattgroesse=rep(1,p))

        MSE_PF_interaction[as.character(PF_splits),as.character(PF_m),as.character(PF_Itersplit_try)]=mean((fits-Data$Y_true)^2)

        TimePF_interaction=Sys.time()-TimePF_interaction
        # 
        #Estimator including arbitrary interaction terms
        
        TimePF_fullinteraction=Sys.time()
        
        fits=planted_forest_2(Y=Data$Y_start,
                              X=Data$X,
                              max_interaction=PF_splits,
                              Baum=50,
                              splits=PF_splits,
                              m=PF_m,
                              Itersplit_try=PF_Itersplit_try,
                              variables=NULL,
                              Blattgroesse=rep(1,p))
        
        MSE_PF_fullinteraction[as.character(PF_splits),as.character(PF_m),as.character(PF_Itersplit_try)]=mean((fits-Data$Y_true)^2)
        
        TimePF_fullinteraction=Sys.time()-TimePF_fullinteraction
        
      }
    }
  }
  
  TimePF=Sys.time()-TimePF
  
  #BF
  # 
  # TimeBackfit=Sys.time()
  # 
  # MSE_backfit.LL=rep(0,length(b.grid))
  # MSE_backfit.LC=rep(0,length(b.grid))
  # MSE_backfit.BF=rep(0,length(b.grid))
  # dim(MSE_backfit.LL)<-dim(MSE_backfit.LC)<-dim(MSE_backfit.BF)<-c(1,length(b.grid))
  # dimnames(MSE_backfit.LL)<-dimnames(MSE_backfit.LC)<-dimnames(MSE_backfit.BF)<-list(1,b.grid)
  # 
  # for(backfit_b.grid in b.grid){
  #   
  #   b.LL.grid <- b.LC.grid <- b.BF.grid <- numeric(p)
  #   b.LL.grid[1:(p+1)] <- b.LC.grid[1:(p+1)] <- b.BF.grid[1:(p+1)] <- rep(c(backfit_b.grid),p+1)
  #   
  #   model_backfit.LL<-SBF.reg.LL(frmla,
  #                                train.data,
  #                                b.LL.grid,
  #                                weight='sw',
  #                                it=15,
  #                                x.grid=NULL,
  #                                integral.approx='midd',
  #                                kcorr=TRUE,
  #                                LC=FALSE,
  #                                wrong=FALSE,
  #                                classic.backfit=FALSE)
  #   
  #   model_backfit.LC<-SBF.reg.LL(frmla,
  #                                train.data,
  #                                b.LC.grid,
  #                                weight='sw',
  #                                it=15,
  #                                x.grid=NULL,
  #                                integral.approx='midd',
  #                                kcorr=TRUE,
  #                                LC=TRUE,
  #                                wrong=FALSE,
  #                                classic.backfit=FALSE)
  #   
  #   model_backfit.BF<-SBF.reg.LL(frmla,
  #                                train.data,
  #                                b.BF.grid,
  #                                weight='sw',
  #                                it=15,
  #                                x.grid=NULL,
  #                                integral.approx='midd',
  #                                kcorr=TRUE,
  #                                LC=TRUE,
  #                                wrong=FALSE,
  #                                classic.backfit=TRUE)
  #   
  #   fits_backfit.LL=0
  #   fits_backfit.LC=0
  #   fits_backfit.BF=0
  #   
  #   for(i in 1:p){
  #     
  #     names(model_backfit.LL$f_backfit[[i]])=names(model_backfit.LL$x.grid[[i]])
  #     names(model_backfit.LC$f_backfit[[i]])=names(model_backfit.LC$x.grid[[i]])
  #     names(model_backfit.BF$f_backfit[[i]])=names(model_backfit.BF$x.grid[[i]])
  #     
  #     fits_backfit.LL=fits_backfit.LL+model_backfit.LL$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.LL$x.grid[[i]]))))]
  #     fits_backfit.LC=fits_backfit.LC+model_backfit.LC$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.LC$x.grid[[i]]))))]
  #     fits_backfit.BF=fits_backfit.BF+model_backfit.BF$f_backfit[[i]][as.character(sort(as.numeric(names(model_backfit.BF$x.grid[[i]]))))]
  #   }
  #   
  #   MSE_backfit.LL[1,as.character(backfit_b.grid)]=mean((fits_backfit.LL-Data$Y_true)^2)
  #   MSE_backfit.LC[1,as.character(backfit_b.grid)]=mean((fits_backfit.LC-Data$Y_true)^2)
  #   MSE_backfit.BF[1,as.character(backfit_b.grid)]=mean((fits_backfit.BF-Data$Y_true)^2)
  # }
  # 
  # TimeBackfit=Sys.time()-TimeBackfit
  # 
  # #Pspline
  # 
  # TimePspline=Sys.time()
  # 
  # MSE_pspline=rep(0,length(df))
  # dim(MSE_pspline)<-c(1,length(df))
  # dimnames(MSE_pspline)<-list(1,df)
  # 
  # for(pspline_df in df){
  #   
  #   pred2 <- paste0("pspline(","V", 1:p,",df=",c(pspline_df),")")
  #   frmla2 <- reformulate(pred2,"Y")
  #   model=lm(frmla2, train.data)
  #   fits=predict(model)
  #   
  #   MSE_pspline[1,as.character(pspline_df)]=mean((fits-Data$Y_true)^2)
  # }
  # 
  # pred3 <- paste0("s(","V", 1:p,")")
  # frmla3 <- reformulate(pred3, "Y")
  # model <- gam(frmla3, data=train.data)
  # fits <- predict(model, train.data)
  # 
  # if(!is.null(fits)){MSE_gam=mean((fits-Data$Y_true)^2)} else {MSE_gam=1000}
  # 
  # TimePspline=Sys.time()-TimePspline
  # 
  # #Gradient boosting
  # 
  # #Additive estimator
  # 
  # TimeGradientBoosting_additive=Sys.time()
  # 
  # MSE_xgboost_additiveES<-MSE_xgboost_additive<-rep(0,length(shrinkage)*length(n.trees))
  # dim(MSE_xgboost_additiveES)<-dim(MSE_xgboost_additive)<-c(length(shrinkage),length(n.trees))
  # dimnames(MSE_xgboost_additiveES)<-dimnames(MSE_xgboost_additive)<-list(shrinkage,n.trees)
  # 
  # for(gbm_shrinkage in shrinkage){
  #   for(gbm_n.trees in n.trees){
  #     
  #     model <- xgboost(data = as.matrix(train.data[,2:(p+1)]),
  #                      label = as.vector(train.data[,"Y"]),
  #                      max.depth = 1,
  #                      eta = gbm_shrinkage,
  #                      nthread = 1,
  #                      verbose = F,
  #                      early_stopping_rounds = NULL,
  #                      nrounds = gbm_n.trees,
  #                      objective = "reg:squarederror")
  #     
  #     fits=predict(model, as.matrix(train.data[,2:(p+1)]))
  #     
  #     MSE_xgboost_additive[as.character(gbm_shrinkage),as.character(gbm_n.trees)]=mean((fits-Data$Y_true)^2)
  #     
  #     model <- xgboost(data = as.matrix(train.data[,2:(p+1)]),
  #                      label = as.vector(train.data[,"Y"]),
  #                      max.depth = 1,
  #                      eta = gbm_shrinkage,
  #                      nthread = 1,
  #                      verbose = F,
  #                      early_stopping_rounds = 50,
  #                      nrounds = gbm_n.trees,
  #                      objective = "reg:squarederror")
  #     
  #     fits=predict(model, as.matrix(train.data[,2:(p+1)]))
  #     
  #     MSE_xgboost_additiveES[as.character(gbm_shrinkage),as.character(gbm_n.trees)]=mean((fits-Data$Y_true)^2)
  #   }
  # }
  # 
  # TimeGradientBoosting_additive=Sys.time()-TimeGradientBoosting_additive
  # 
  # TimeAll=Sys.time()-TimeAll
#   
   return(list(
     #MSE_xgboost_additive=MSE_xgboost_additive,
#               MSE_xgboost_additiveES=MSE_xgboost_additiveES,
#               MSE_xgboost_interactionES=MSE_xgboost_interactionES,
#               MSE_xgboost_interaction=MSE_xgboost_interaction,
#               MSE_EBM_additive2=MSE_EBM_additive2,
#               MSE_EBM_additive=MSE_EBM_additive,
#               MSE_PF_alternative=MSE_PF_alternative,
#               MSE_PF_additive=MSE_PF_additive,
               MSE_PF_interaction=MSE_PF_interaction,
               MSE_PF_fullinteraction=MSE_PF_fullinteraction
#               MSE_backfit.LL=MSE_backfit.LL,
#               MSE_backfit.LC=MSE_backfit.LC,
#               MSE_backfit.BF=MSE_backfit.BF,
#               MSE_pspline=MSE_pspline,
#               MSE_gam=MSE_gam,
#               MSE_randomForests=MSE_randomForests,
#              TimeAll=TimeAll
  ))
}

planted_forest_2<- function(Y, X, max_interaction=2, Baum=50, splits=30, m=10, Itersplit_try=0.4, variables=NULL, Blattgroesse=rep(1,p), alternative=F)
{
  force(Itersplit_try)
  force(X)
  force(Y)
  
  # Baum= Anzahl der B?um-famililien im "Forest"
  # splits= Anzahl der Iterationsschritte.
  # m= number of possible split points
  # Blattgroesse= minimal leaf size
  
  p <- ncol(X)
  n <- nrow(X)
  
  Schleife <- function(run){
    
    X_orig=X
    
    #???K_active=0
    # Berechnen des bootstap sample
    
    subsample <- sample(n,n,replace=TRUE)
    
    X <- X[subsample,]
    Y <- Y[subsample]
    
    # Definition der verschiedenen Bl?cke
    
    # Koordinatengruppenindizees des i-ten Eintrags
    # variables[[i]] is a vector of variables  in tree i
    
    if (is.null(variables)){
      variables=list()
      for(i in 1:p){
        variables[[i]]=i
      }
    }
    
    # Definition der Partitionen
    # individuals[[i]][[k]] is a vector of individuals in leaf k of tree i
    
    individuals=list()
    
    for(i in 1:length(variables)){
      individuals[[i]]=list(1:n)
    }
    
    individuals_orig=list()
    
    for(i in 1:length(variables)){
      individuals_orig[[i]]=list(1:n)
    }
    
    Y_fitted=rep(0,n)
    
    Possible_Splits=list()
    for(i_1 in 1:length(variables)){
      for(i_2 in variables[[i_1]]){
        Possible_Splits[[i_1]]=list(i_2,variables[[i_1]])
      }
    }
    
    while(splits>0){
      
      splits=splits-1
      
      split_try=ceiling(Itersplit_try*length(Possible_Splits))
      
      split_candidates <- sample(Possible_Splits, split_try)
      
      # Hilfsvariablen
      
      # R=Matrix der Residuen bzgl. des besten splits f?r verschiedenen Koordinatengruppen K[[i]],
      # Ikx=Indizes der besten zu splitenden Menge, der zu splitenden Koordinate und des zugeh?rigen Splitpunkts f?r verschiedene Koordinatengruppen K.
      
      R=Calc_Optimal_split2(Y, X, m, variables, individuals, Blattgroesse, split_candidates)
      
      R_opt <- R[1]
      
      Ikx_opt <- R[2:5]
      
      if(R_opt<Inf){
        
        if(max_interaction>1){
          
          if(length(individuals[[Ikx_opt[1]]])==1){
            
            for(i_1 in ((1:p)[-variables[[Ikx_opt[1]]]]) ){
              
              Possible_exists=0
              
              for(i_2 in 1:length(Possible_Splits)){
                
                if(Possible_Splits[[i_2]][[1]]==i_1 & length(Possible_Splits[[i_2]][[2]])==length(variables[[Ikx_opt[1]]])+1){
                  
                  if(all(sort(unique(i_1,variables[[Ikx_opt[1]]])) == Possible_Splits[[i_2]][[2]])){
                    
                    Possible_exists=1
                  }
                } 
              }
              if(Possible_exists==0){  
                
                Possible_Splits[[length(Possible_Splits)+1]]=list(i_1,sort(c(i_1,variables[[Ikx_opt[1]]])))
              }
            }
          }
        }
        
        I_2<-individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]],Ikx_opt[3]]<Ikx_opt[4]]
        I_1<-individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]],Ikx_opt[3]]>=Ikx_opt[4]]
        
        I_22<-individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]][X_orig[individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]],Ikx_opt[3]]<Ikx_opt[4]]
        I_12<-individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]][X_orig[individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]],Ikx_opt[3]]>=Ikx_opt[4]]
        
        
        y_2 <- mean(Y[I_2])
        y_1 <- mean(Y[I_1])
        
        # Updaten der Y
        Y[I_2] <- Y[I_2]-y_2
        Y[I_1] <- Y[I_1]-y_1
        
        Y_fitted[I_22] <- Y_fitted[I_22]+y_2
        Y_fitted[I_12] <- Y_fitted[I_12]+y_1
        
        if(Ikx_opt[3] %in% variables[[Ikx_opt[1]]]){  ### if split variable is already in tree to be split
          
          individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])+1]] <- I_2      #### individuals for  new leaf 
          individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_1                       #### new split-leaf = remaining individuals
          
          individuals_orig[[Ikx_opt[1]]][[length(individuals_orig[[Ikx_opt[1]]])+1]] <- I_22      #### individuals for  new leaf 
          individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_12
          
          if(alternative){
            
            for(j in 1:length(individuals[[Ikx_opt[1]]])){
              
              # Updaten der Y
              Y[individuals[[Ikx_opt[1]]][[j]]] <- Y[individuals[[Ikx_opt[1]]][[j]]]-mean(Y[individuals[[Ikx_opt[1]]][[j]]])
              
              Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]] <- Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]]-mean(Y[individuals_orig[[Ikx_opt[1]]][[j]]])
            }
            
          }
          
        } else {
          
          TreeExists=0
          
          for(i in 1:length(variables)){
            
            if(length(variables[[i]])==length(variables[[Ikx_opt[1]]])+1){ ###if number of variables in tree i is same as num of var in new split tree
              if(all(variables[[i]]==sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3])))){ ### if variables in tree i are same as  var in new split tree
                
                TreeExists=1  ###  tree i is the same as new split tree
                
                individuals[[i]][[length(individuals[[i]])+1]] <- I_2
                individuals[[i]][[length(individuals[[i]])+1]] <- I_1
                
                individuals_orig[[i]][[length(individuals_orig[[i]])+1]] <- I_22
                individuals_orig[[i]][[length(individuals_orig[[i]])+1]] <- I_12
              }
            }
          }
          
          if(TreeExists==0){
            
            variables[[length(variables)+1]] <- sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3]))
            individuals[[length(individuals)+1]] <- list()
            individuals[[length(individuals)]][[1]] <- I_2
            individuals[[length(individuals)]][[2]] <- I_1
            
            individuals_orig[[length(individuals_orig)+1]] <- list()
            individuals_orig[[length(individuals_orig)]][[1]] <- I_22
            individuals_orig[[length(individuals_orig)]][[2]] <- I_12
            
            if(length(variables[[length(variables)]])<max_interaction){
              
              for(i_1 in ((1:p)[-variables[[length(variables)]]]) ){
                Possible_exists=0
                for(i_2 in 1:length(Possible_Splits)){
                  if(Possible_Splits[[i_2]][[1]]==i_1 & length(Possible_Splits[[i_2]][[2]])==length(variables[[length(variables)]])+1){
                    if(all(sort(unique(i_1,variables[[length(variables)]]))==Possible_Splits[[i_2]][[2]])){
                      Possible_exists=1
                    }
                  } 
                }
                if(Possible_exists==0){  
                  Possible_Splits[[length(Possible_Splits)+1]]=list(i_1,sort(c(i_1,variables[[length(variables)]])))
                }
              }
              
              for(i_1 in variables[[length(variables)]]){
                Possible_exists=0
                for(i_2 in 1:length(Possible_Splits)){
                  if(Possible_Splits[[i_2]][[1]]==i_1 & length(Possible_Splits[[i_2]][[2]])==length(variables[[length(variables)]])){
                    if(all(variables[[length(variables)]] == Possible_Splits[[i_2]][[2]])){
                      Possible_exists=1
                    }
                  } 
                }
                if(Possible_exists==0){  
                  Possible_Splits[[length(Possible_Splits)+1]]=list(i_1,variables[[length(variables)]])
                }
              }
            }
          }
        }
      }
    }
    
    return(Y_fitted)
  }
  
  forest_res <- sapply(1:Baum, Schleife) 
  
  Y_fitted=0
  
  for(i in 1:Baum){
    Y_fitted=Y_fitted+forest_res[,i]
  }
  
  Y_fitted=Y_fitted/Baum
  
  return(Y_fitted)
}

cl <- makeCluster(Kerne)

clusterExport(cl, varlist = ls(), envir = .GlobalEnv)

Results <- parSapply(cl,1:Monte_Carlo, simulate_function)

stopCluster(cl)

save(Results,file="Model2-K=10_nosparse_ResultsHP60.Rdata")