xgboost_CV<- function(Y, X, CV=2)
{
  force(X)
  force(Y)
  
  p <- ncol(X)
  
  shrinkage=c(0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
  n.trees=c(100, 300, 600, 1000 , 3000, 5000, 7000)
  interaction_depth=c(1,2,3,4)
  
  # shrinkage=c(0.005, 0.01)
  # n.trees=c(100, 300)
  # interaction_depth=c(1,2)
  
  MSE_matrix<-rep(0,length(shrinkage)*length(n.trees)*length(interaction_depth))
  dim(MSE_matrix)<-c(length(shrinkage),length(n.trees),length(interaction_depth))
  dimnames(MSE_matrix)<-list(shrinkage,n.trees,interaction_depth)
  
  # splits_test=c(10,15)
  # m_test=c(2,5)
  # Itersplit_try_test=c(0.25,0.5)
  # max_interaction_test=c(1,2)
  
  # Baum= Anzahl der B?um-famililien im "Forest"
  # splits= Anzahl der Iterationsschritte.
  # m= number of possible split points
  # Blattgroesse= minimal leaf size
  
  X_orig=X
  
  Y_orig=Y
  
  MSE_CV=Inf
  
  for(gbm_shrinkage in shrinkage){
    
    for(gbm_n.trees in n.trees){
      
      for(gbm_interaction_depth in interaction_depth){
          
          MSE=0
          
          for(fold in 1:CV){
            
            indizees_CV=(1:nrow(X_orig))[(((fold-1)*nrow(X_orig)/CV+1):(fold*nrow(X_orig)/CV))]
            
            X=X_orig[-indizees_CV,]
            
            Y=Y_orig[-indizees_CV]
            
            p <- dim(X)[2]
            
            train.data <- data.frame(cbind(Y,X))
            names(train.data)[1] <- "Y"
            names(train.data)[2:(p+1)] <- paste0("V", 1:p)
            
            test.data <- data.frame(cbind(Y_orig[indizees_CV],X_orig[indizees_CV,]))
            names(test.data)[1] <- "Y"
            names(test.data)[2:(p+1)] <- paste0("V", 1:p)
            
            model<- xgboost(data = as.matrix(train.data[,2:(p+1)]),
                            label = as.vector(train.data[,"Y"]),
                            max.depth = gbm_interaction_depth,
                            eta = gbm_shrinkage,
                            nthread = 1,
                            nrounds = gbm_n.trees,
                            early_stopping_rounds = NULL,
                            verbose = F,
                            objective = "reg:squarederror")
            
            fits=predict(model, as.matrix(test.data[,2:(p+1)]))
            
            MSE=MSE+sum((fits-Y_orig[indizees_CV])^2)
          }
          
          MSE_matrix[as.character(gbm_shrinkage),as.character(gbm_n.trees),as.character(gbm_interaction_depth)]=MSE
          
          if(MSE<MSE_CV){
            MSE_CV=MSE
            shrinkage_solution=gbm_shrinkage
            n.trees_solution=gbm_n.trees
            interaction_depth_solution=gbm_interaction_depth
          }
      }
    }
  }
  
  return(list(shrinkage=shrinkage_solution, n.trees=n.trees_solution, interaction_depth=interaction_depth_solution, matrix=MSE_matrix))
}
