# Calculate the estimator for a value from a single family of trees
# Input: x = Input vector, s = Index of the family of trees, i = Coordinates of the component to be estimated, forest_res = planted forest model
# Output: Estimated function f(x) if i=0, otherwise estimated component f_i(x)

f_single_fam=function(x,s=1,i=0,forest_res){
  
  f=0
  
  if(length(i)==1){
    
    if(i==0){
      
      for(i in 1:length(forest_res[[3,s]])){
        
        for(j in 1:length(forest_res[[1,s]][[i]])){
          
          if(prod(forest_res[[1,s]][[i]][[j]][1,forest_res[[3,s]][[i]]]<= x[forest_res[[3,s]][[i]]] & (forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]>x[forest_res[[3,s]][[i]]]|forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]==forest_res[[4,s]][forest_res[[3,s]][[i]]]))){
            
            f=f+forest_res[[2,s]][[i]][j]
          }   
        }
      }
      
      if(is.null(forest_res[,s]$constant)){ return(f) }
      
      return(f+forest_res[["constant",s]])
    }
  }
  
  for(i_0 in 1:length(forest_res[[3,s]])){
    
    if(length(forest_res[[3,s]][[i_0]])==length(i)){
      
      if(prod(forest_res[[3,s]][[i_0]]==sort(i))){
        
        for(j in 1:length(forest_res[[1,s]][[i_0]])){
          
          if(prod(forest_res[[1,s]][[i_0]][[j]][1,forest_res[[3,s]][[i_0]]]<= x & (forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]]]>x|forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]]]==forest_res[[4,s]][forest_res[[3,s]][[i_0]]]))){
            
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

planted_forest_f=function(x, forest_res,f_single_fam, i=0){
  
  ntrees=dim(forest_res)[2]
  
  y=0
  
  for(s in 1:ntrees){ y=y+f_single_fam(x, s, i, forest_res) }
  
  return(y/ntrees)
}
