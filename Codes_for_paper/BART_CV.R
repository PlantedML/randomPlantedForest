BART_CV<- function(Y, X, CV=2)
{
  force(X)
  force(Y)
  
  p <- ncol(X)
  
  power=1:5
  sparsity_parameter=c(0.6,0.75,0.9)
  ntree=c(50,100,150,200,250,300)
  
  MSE_matrix<-1:(length(power)*length(sparsity_parameter)*length(ntree))
  dim(MSE_matrix)<-c(length(power),length(sparsity_parameter),length(ntree))
  dimnames(MSE_matrix)<-list(power,sparsity_parameter,ntree)
  
  X_orig=X
  
  Y_orig=Y
  
  MSE_CV=Inf
  
  for(n_pow in power){
    
    for(n_spars in sparsity_parameter){
      
      for(n_tree in ntree){
        
        MSE=0
        
        for(fold in 1:CV){
          
          indizees_CV=(1:nrow(X_orig))[(((fold-1)*nrow(X_orig)/CV+1):(fold*nrow(X_orig)/CV))]
          
          X=X_orig[-indizees_CV,]
          
          Y=Y_orig[-indizees_CV]
          
          p <- dim(X)[2]
          
          invisible(capture.output(post <- wbart(x.train=X, y.train=Y, power=n_pow, a=n_spars, ntree=n_tree, ndpost = 1600)))
          
          invisible(capture.output(fits <- predict(post, X_orig[indizees_CV,])))
          
          Y_rep=matrix(rep(Y_orig[indizees_CV],dim(fits)[1]), ncol=length(indizees_CV), byrow= T)
          
          MSE=MSE+mean(colMeans(fits-Y_rep)^2)
        }
        
        MSE_matrix[as.character(n_pow),as.character(n_spars),as.character(n_tree)]=MSE
        
        if(MSE<MSE_CV){
          MSE_CV=MSE
          n_pow_solution=n_pow
          n_spars_solution=n_spars
          n_tree_solution=n_tree
        }
      }
    }
  }
  
  return(list(n_pow=n_pow_solution, n_spars=n_spars_solution, n_tree=n_tree_solution, matrix=MSE_matrix))
}
