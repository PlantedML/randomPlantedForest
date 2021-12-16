# Implementation of the random planted forest algorithm used for the simulations where NO plots where made 
# Input: (X,Y) = Data, max_interaction = maximum order of interaction, ntrees = Number of families of trees used, splits = number of splits,
#        split_try = number of considered splitpoints for each interval in each iteration step, t_try = percentage of viable trees considered in each iteration step,
#        variables = list of trees in each family in the beginning (default: one tree for each coordinate), 
#        leaf_size = minimum number of nodes in each leaf, alternative = alternative updating       
# Output: Fitted Y_values

planted_forest_simu<- function(Y, X, max_interaction=2, ntrees=50, splits=30, split_try=10, t_try=0.4, variables=NULL, leaf_size=rep(1,p), alternative=F)
{
  force(t_try)
  force(X)
  force(Y)
  
  p <- ncol(X)
  n <- nrow(X)
  
  tree_fam <- function(run){
    
    X_orig=X
    
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
      individuals_orig[[i]]=list(1:n)
    }
    
    Y_fitted=rep(0,n)
    
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
          
          individuals_orig[[Ikx_opt[1]]][[length(individuals_orig[[Ikx_opt[1]]])+1]] <- I_22      #### individuals for  new leaf 
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
              if(all(variables[[i]]==sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3])))){ ### if variables in tree i are same as  var in new split tree
                
                TreeExists=1  ### tree i is the same as new split tree
                
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
  
  forest_res <- sapply(1:ntrees, tree_fam) 
  
  Y_fitted=0
  
  for(i in 1:ntrees){
    Y_fitted=Y_fitted+forest_res[,i]
  }
  
  Y_fitted=Y_fitted/ntrees
  
  return(Y_fitted)
}

# Implementation of the random planted forest algorithm used for the simulations where plots were made
# Input: (X,Y) = Data, max_interaction = maximum order of interaction, ntrees = Number of families of trees used, splits = number of splits,
#        split_try = number of considered splitpoints for each interval in each iteration step, t_try = percentage of viable trees considered in each iteration step,
#        variables = list of trees in each family in the beginning (default: one tree for each coordinate), 
#        leaf_size = minimum number of nodes in each leaf, alternative = alternative updating       
# Output: list: [1] = Fitted Y_values, 
#               [2] = planted forest model: [2][1] Final leaves of the trees, [2][2] = estimated values corresponding to the leaves,
#               [2][3] = coordinates of the trees 

planted_forest_simu_plots<- function(Y, X, max_interaction=2, ntrees=50, splits=30, split_try=10, t_try=0.4, variables=NULL, leaf_size=rep(1,p), alternative=F){
  
  force(t_try)
  
  p <- ncol(X)
  
  a <- apply(X,2,min)     ## lower bounds
  b <- apply(X,2,max)     ### upper bounds
  
  tree_fam <- function(run){
    
    subsample <- sample(nrow(X),nrow(X),replace=TRUE)
    
    X <- X[subsample,]
    Y <- Y[subsample]
    
    # variables[[i]] is a vector of variables  in tree i
    if (is.null(variables)){
      variables=list()
      for(i in 1:p){
        variables[[i]]=i
      }
    }
    # intervals[[i]][[k]]= is list of matrices describing intervals for tree i, leaf k  
    
    intervals=list()
    
    for(i in 1:length(variables)){
      intervals[[i]] <- list()
      intervals[[i]][[1]] <- matrix(nrow=2, ncol=p)
      for(j in 1:p){
        intervals[[i]][[1]][,j]=c(a[j],b[j])
      }
    }
    
    # values[[i]] is a vector of predictions for each leaf in tree i
    
    values=list()
    
    for(i in 1:length(variables)){
      values[[i]]=0
    }
    
    # individuals[[i]][[k]] is a vector of individuals in leaf k of tree i
    
    individuals=list()
    
    for(i in 1:length(variables)){
      individuals[[i]]=list(1:nrow(X))
    }
    
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
        
        y_2 <- mean(Y[I_2])
        y_1 <- mean(Y[I_1])
        
        Y[I_2] <- Y[I_2]-y_2
        Y[I_1] <- Y[I_1]-y_1
        
        if(Ikx_opt[3] %in% variables[[Ikx_opt[1]]]){  ### if split variable is already in tree to be split
          
          individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])+1]] <- I_2      #### individuals for  new leaf 
          individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_1                       #### new split-leaf = remaining individuals
          
          intervals[[Ikx_opt[1]]][[length(intervals[[Ikx_opt[1]]])+1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]] #add one leaf to intervals
          
          intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][2,Ikx_opt[3]] <- Ikx_opt[4]  ### new leaf has new interval at split variable: (split value= upper bopund, lower bound remains)
          intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][1,Ikx_opt[3]] <- Ikx_opt[4] ### split leaf has new interval at split variable: (split value= lower bound, upper bound remains)
          
          values[[Ikx_opt[1]]][length(individuals[[Ikx_opt[1]]])] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_2
          values[[Ikx_opt[1]]][Ikx_opt[2]] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_1
          
          if(alternative){
            
            for(j in 1:length(intervals[[Ikx_opt[1]]])){
              
              y <- mean(Y[individuals[[Ikx_opt[1]]][[j]]])
              
              Y[individuals[[Ikx_opt[1]]][[j]]] <- Y[individuals[[Ikx_opt[1]]][[j]]]-y
              
              values[[Ikx_opt[1]]][j] <- values[[Ikx_opt[1]]][j]+y
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
                
                intervals[[i]][[length(individuals[[i]])]] <- intervals[[i]][[length(individuals[[i]])-1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]]
                
                intervals[[i]][[length(individuals[[i]])-1]][2,Ikx_opt[3]] <- Ikx_opt[4]
                intervals[[i]][[length(individuals[[i]])]][1,Ikx_opt[3]] <- Ikx_opt[4]
                
                values[[i]][length(individuals[[i]])-1]=y_2
                values[[i]][length(individuals[[i]])]=y_1
              }
            }
          }
          
          if(TreeExists==0){
            
            variables[[length(variables)+1]] <- sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3]))
            individuals[[length(individuals)+1]] <- list()
            individuals[[length(individuals)]][[1]] <- I_2
            individuals[[length(individuals)]][[2]] <- I_1
            
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
            
            intervals[[length(intervals)+1]] <- list()
            intervals[[length(intervals)]][[1]] <- intervals[[length(intervals)]][[2]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]]
            
            intervals[[length(intervals)]][[1]][2,Ikx_opt[3]] <- Ikx_opt[4]
            intervals[[length(intervals)]][[2]][1,Ikx_opt[3]] <- Ikx_opt[4]
            
            values[[length(variables)]] <- y_2
            values[[length(variables)]][2] <- y_1
            
          }
        }
      }
    }
    
    # Purity Algorithm
    
    constant=0
    
    for(l in max_interaction:1){
      
      for(i_0 in 1:length(variables)){
        
        if(length(variables[[i_0]])==l){
          
          if(l==1){
            
            i_1=variables[[i_0]]
            
            intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[1]]
            
            intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
            
            values[[i_0]][[length(intervals[[i_0]])]]=0
            
            for(i_3 in 1:(length(intervals[[i_0]])-1)){
              
              values[[i_0]][[length(intervals[[i_0]])]]=values[[i_0]][[length(intervals[[i_0]])]]-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
              constant=constant+values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
            }
            
          } else {
            
            for(i_1 in 1:dim(X)[1]){
              
              if(i_1 %in% variables[[i_0]]){
                
                exists=0
                
                for(i_2 in 1:length(variables)){
                  
                  if(length(variables[[i_2]])==l-1){
                    
                    if(all(variables[[i_2]]==variables[[i_0]][variables[[i_0]]!=i_1])){
                      
                      exists=1
                      
                      for(i_3 in 1:length(intervals[[i_0]])){
                        
                        exists_2=0
                        
                        intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[i_3]]
                        
                        intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                        
                        values[[i_0]][[length(intervals[[i_0]])]]=-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        
                        for(i_4 in 1:length(intervals[[i_2]])){
                          
                          Cube=intervals[[i_0]][[i_3]]
                          
                          Cube[,i_1]=c(a[i_1],b[i_1])
                          
                          if(all(intervals[[i_2]][[i_4]]==Cube) & exists_2==0){
                            
                            exists_2=1
                            
                            values[[i_2]][[i_4]]=values[[i_2]][[i_4]]+values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                          } 
                        }
                        
                        if(exists_2==0){
                          
                          intervals[[i_2]][[length(intervals[[i_2]])+1]]=intervals[[i_0]][[length(intervals[[i_0]])]]
                          
                          values[[i_2]][[length(intervals[[i_2]])]]=values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        }
                      }
                    }
                  }
                }
                
                if(exists==0){
                  
                  variables[[length(variables)+1]]=variables[[i_0]][variables[[i_0]]!=i_1]
                  
                  intervals[[length(variables)]]=list()
                  
                  values[[length(variables)]]=vector()
                  
                  for(i_3 in 1:length(intervals[[i_0]])){
                    
                    intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[i_3]]
                    
                    intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                    
                    values[[i_0]][[length(intervals[[i_0]])]]=-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                    
                    intervals[[length(variables)]][[length(intervals[[length(variables)]])+1]]=intervals[[i_0]][[length(intervals[[i_0]])]]
                    
                    values[[length(variables)]][[length(intervals[[length(variables)]])]]=-values[[i_0]][[length(intervals[[i_0]])]]
                    
                  }
                }
              }
            }
          }
        }
      }
    }
    
    return(list(intervals=intervals, values=values, variables=variables, constant=constant))
  }
  
  forest_res <- sapply(1:ntrees, tree_fam) 
  
  # Calculate the estimator for a value from a single family of trees
  # Input: x = Input vector, s = Index of the family of trees, forest_res = planted forest model  
  # Output: Estimated function f(x)
  
  F_single_fam=function(x,s=1,forest_res){
    
    f=0
    
    for(i in 1:length(forest_res[[3,s]])){
      
      for(j in 1:length(forest_res[[1,s]][[i]])){
        
        if(prod(forest_res[[1,s]][[i]][[j]][1,forest_res[[3,s]][[i]]]<= x[forest_res[[3,s]][[i]]] & (forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]>x[forest_res[[3,s]][[i]]]|forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]==b[forest_res[[3,s]][[i]]]))){
          
          f=f+forest_res[[2,s]][[i]][j]
        }   
      }
    }
    
    return(f+forest_res[["constant",s]])  
    
  }
  
  # Calculate the average of the estimators of the individual families of trees
  # Input:  x = Input vector
  # Output: Estimated function f(x)
  
  F_average=function(x){
    
    y=0
    
    for(s in 1:ntrees){ y=y+F_single_fam(x,s,forest_res) }
    
    return(y/ntrees)
  }
  
  Y_fitted=0
  
  for(j in 1:nrow(X)){ Y_fitted[j]=F_average(X[j,]) }
  
  return(list(Y_fitted=Y_fitted,forest_res=forest_res))
}