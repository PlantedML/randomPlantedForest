# Calculate the estimator for a value from a single family of trees
# Input: x = Input vector, s = Index of the family of trees, i = Coordinates of the component to be estimated, forest_res = planted forest model
# Output: Estimated function f(x) if i=0, otherwise estimated component f_i(x)

F_single_fam=function(x,s=1,i=0,forest_res){
  
  f=0
  
  a <- forest_res[,s][[8]]
  b <- forest_res[,s][[4]]
  
  x<- sapply(1:length(a), function(j) max(a[j],x[j]))
  x<- sapply(1:length(a), function(j) min(b[j],x[j]))
  
  
  if(length(i)==1){
    
    if(i==0){
      
      for(i in 1:length(forest_res[[3,s]])){  ###  tree 3=variables  --> i= trees
        
        variables <- forest_res[,s][[3]][[i]]
        
        for(j in 1:length(forest_res[[1,s]][[i]])){ ## 1= intervals --> j is leaf of tree i
          
          
          in_leaf <- sapply(variables, function (l) {
            
            categorical = is.element(l,forest_res[[6,s]])
            if(categorical){
              
              return ( is.element(x[l] , forest_res[[1,s]][[i]][[j]][,l]))
              
            } else {
              
              return((forest_res[[1,s]][[i]][[j]][1,l]<= x[l])&  (forest_res[[1,s]][[i]][[j]][2,l] >x[l] ))
            }
            
            
          })
          if (all(in_leaf))   f=f+forest_res[[2,s]][[i]][j]
        }   
        
      }
      # f<-mean(f[-1])
      if(is.null(forest_res[,s]$constant)){ return(f) }
      
      return(f+forest_res[["constant",s]])
    }
  }
  
  for(i_0 in 1:length(forest_res[[3,s]])){
    
    if(length(forest_res[[3,s]][[i_0]])==length(i)){
      
      if(prod(forest_res[[3,s]][[i_0]]==sort(i))){  ## variables in tree i_0  match i
        
        for(j in 1:length(forest_res[[1,s]][[i_0]])){  ### j is leaf in tree i_0
          
          categorical = is.element(i,forest_res[[6,s]])
          
          if(all(forest_res[[1,s]][[i_0]][[j]][1,forest_res[[3,s]][[i_0]][!categorical]]<= x[!categorical] & (forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]][!categorical]]>x[!categorical]|forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]][!categorical]]==forest_res[[4,s]][forest_res[[3,s]][[i_0]][!categorical]]))& min(1,as.numeric(sapply( i[categorical] , function (k) is.element(x[which(i==k)],  forest_res[[1,s]][[i]][[j]][,k]   ))))){
            
            f=f+forest_res[[2,s]][[i_0]][j]
          }   
        }
      }
    }
  }
  
  return(f) 
}

# Calculate the average of the estimators of the individual families of trees
# Input:  x = Input vector, forest_res = random planted forest model, i = Coordinates of the component to be estimated
# Output: Estimated function f(x) if i=0, otherwise estimated component f_i(x)

F_average=function(x, forest_res, i=0){
  
  ntrees=dim(forest_res)[2]
  
  y=0
  
  for(s in 1:ntrees){ y=y + F_single_fam(x, s, i, forest_res) }
  
  return(y/ntrees)
}

# Predict a value using the random planted forest model
# Input:  x = Input vector, forest_res = random planted forest model, i = Coordinates of the component to be estimated
# Output: Estimated function f(x) if i=0, otherwise estimated component f_i(x)

predict_rpf=function(X, forest_res, i=0){
  
  if(length(i)==1){
    
    if(i==0){
      
      p = dim(forest_res[,1]$intervals[[1]][[1]])[2] 
      
      if(length(X)==p){ return(F_average(X,forest_res)) }
      
      if(!is.matrix(X)){ return("The input X has the wrong dimension in order to calculate f(x)") }
      
      if(dim(X)[2] != p){ return("The input X has the wrong dimension in order to calculate f(x)") }
      
      Y=apply(X,1,F_average, forest_res=forest_res)
      
      return(Y)
    }
    
    if(length(X)>1){
      
      dim(X)=c(length(X),1)
    }
  }
  
  if(is.null(forest_res[,1]$constant)){ print("Note: The estimator has not been purified.") }
  
  if(length(X)==length(i)){ return(F_average(X,forest_res, i=i)) }
  
  if(!is.matrix(X)){ return("The input X has the wrong dimension in order to calculate f(x)") }
  
  if(dim(X)[2]!=length(i)){ return("The input X has the wrong dimension in order to calculate f_i(x)") }
  
  Y=apply(X,1,F_average, forest_res=forest_res, i=i)
  
  return(Y)
}