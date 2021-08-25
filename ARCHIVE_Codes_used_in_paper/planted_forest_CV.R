# Calculate optimal parameters for the random planted forests algorithm using Cross-Validation
# Input: (X,Y) = Data, ntrees = Number of families of trees used, 
#        variables = list of trees in each family in the beginning (default: one tree for each coordinate), 
#        leaf_size = minimum number of nodes in each leaf, alternative = alternative updating, CV = Number of folds          
# Output: list: [1] = Number of splits, [2] = t_try, [3] = maximum oder of interaction

planted_forest_CV<- function(Y, X, ntrees=50, variables=NULL, leaf_size=rep(1,p), alternative=F, CV=2)
{
  force(X)
  force(Y)
  
  p <- ncol(X)
  
  splits_test=c(10,15,20,25,30)
  split_try_test=c(2,5,10,20)
  t_try_test=c(0.25,0.5,0.75)
  max_interaction_test=c(1,2,p)
  
  X_orig=X
  
  Y_orig=Y
  
  # function to create a family of trees
  
  tree_fam <- function(run){
    
    # Calculate the bootstap sample
    
    n <- nrow(X)
    
    subsample <- sample(n,n,replace=TRUE)
    
    X <- X[subsample,]
    Y <- Y[subsample]
    
    # variables[[i]] is a vector of variables  in tree i
    
    if (is.null(variables)){
      variables=list()
      for(i in 1:p){
        variables[[i]]=i
      }
    }
    
    # individuals[[i]][[k]] is a vector of indices of the bootstrap sample in leaf k of tree i
    
    individuals=list()
    
    for(i in 1:length(variables)){
      individuals[[i]]=list(1:n)
    }
    
    # individuals_orig[[i]][[k]] is a vector of indices of the original sample in leaf k of tree i 
    
    individuals_orig=list()
    
    for(i in 1:length(variables)){
      individuals_orig[[i]]=list(1:nrow(X_orig))
    }
    
    Y_fitted=rep(0,nrow(X_orig))
    
    # Possible_Splits corresponds to the set of viable splits. 
    # Possible_Splits[[i_1]][[1]] is the coordinate used for splitting, Possible_Splits[[i_1]][[2]] is the tree which results from the splitting 
    
    Possible_Splits=list()
    for(i_1 in 1:length(variables)){
      for(i_2 in variables[[i_1]]){
        Possible_Splits[[i_1]]=list(i_2,variables[[i_1]])
      }
    }
    
    while(splits>0){
      
      splits=splits-1
      
      m_try=ceiling(t_try*length(Possible_Splits))
      
      split_candidates <- sample(Possible_Splits, m_try)
      
      R=Calc_Optimal_split2(Y, X, split_try, variables, individuals, leaf_size, split_candidates)
      
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
        
        Y[I_2] <- Y[I_2]-y_2
        Y[I_1] <- Y[I_1]-y_1
        
        Y_fitted[I_22] <- Y_fitted[I_22]+y_2
        Y_fitted[I_12] <- Y_fitted[I_12]+y_1
        
        if(Ikx_opt[3] %in% variables[[Ikx_opt[1]]]){  ### if split variable is already in tree to be split
          
          individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])+1]] <- I_2      #### individuals for  new leaf 
          individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_1                       #### new split-leaf = remaining individuals
          
          individuals_orig[[Ikx_opt[1]]][[length(individuals_orig[[Ikx_opt[1]]])+1]] <- I_22      #### individuals for new leaf 
          individuals_orig[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_12
          
          if(alternative){
            
            for(j in 1:length(individuals[[Ikx_opt[1]]])){
              
              Y[individuals[[Ikx_opt[1]]][[j]]] <- Y[individuals[[Ikx_opt[1]]][[j]]]-mean(Y[individuals[[Ikx_opt[1]]][[j]]])
              
              Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]] <- Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]]-mean(Y[individuals_orig[[Ikx_opt[1]]][[j]]])
            }
            
          }
          
        } else {
          
          TreeExists=0
          
          for(i in 1:length(variables)){
            
            if(length(variables[[i]])==length(variables[[Ikx_opt[1]]])+1){ ### if number of variables in tree i is same as num of var in new split tree
              if(all(variables[[i]]==sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3])))){ ### if variables in tree i are same as var in new split tree
                
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
  
  MSE_CV=Inf
  
  for(split_try in split_try_test){
    
    for(splits in splits_test){
    
      for(t_try in t_try_test){
          
        for(max_interaction in max_interaction_test){
        
          MSE=0
            
          for(fold in 1:CV){
              
            indizees_CV=(1:nrow(X_orig))[-(((fold-1)*nrow(X_orig)/CV):(fold*nrow(X_orig)/CV))]
              
            X=X_orig[-indizees_CV,]
              
            Y=Y_orig[-indizees_CV]
              
            forest_res <- sapply(1:ntrees, tree_fam) 
              
            Y_fitted=0
              
            for(i in 1:ntrees){ Y_fitted=Y_fitted+forest_res[,i] }
              
            Y_fitted=Y_fitted/ntrees
              
            MSE=MSE+sum((Y_fitted[indizees_CV]-Y_orig[indizees_CV])^2)
          }
            
          if(MSE<MSE_CV){

            MSE_CV=MSE
            splits_solution=splits
            split_try_solution=split_try
            t_try_solution=t_try
            max_interaction_solution=max_interaction
          }
        }
      }
    }
  }
  
  return(list(splits=splits_solution, split_try=split_try_solution, t_try=t_try_solution, max_interaction=max_interaction_solution))
}