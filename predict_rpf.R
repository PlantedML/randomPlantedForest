
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
