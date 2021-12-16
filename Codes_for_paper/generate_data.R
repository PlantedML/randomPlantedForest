# Data generation for simulation
# Input: p = Dimension of X_i, n = Number of Data points, rho = Correlation coefficient, sparsity = number of components active
# Output: list: [1] = Y-values, [2] = Denoised Y-values, [3] = X-values

generate_data <-function(n=500,p=4,rho=0.3,sparsity=2,sigma=1, sparse = T, Model=1, covariates='normal'){
  
  if(covariates=='normal'){
    
    z_1 <- rnorm(p*n,0,(1-rho)^(1/2))
    z_0 <- rnorm(n,0,rho^(1/2))
    X <- 2.5*atan(z_1+z_0)/pi
    dim(X) <- c(n,p)
  }
  
  if(covariates=='uniform'){
    
    X <- matrix( runif(p*n,-1,1), nrow=n, ncol=p)
  }
  
  epsilon <- rnorm(nrow(X),0,sigma)
  
  ### Model1 = additive 
  
  if (Model==1){
    
    F_1=function(x){
      
      return(2*sin(pi*x))
    }
    
    if(sparse == F && !(p==30)){ sparsity = p }
    
    F_ursprung=function(x){
      
      m <- rep(0,sparsity)
      
      for(j in 1:sparsity){ m[j]=(-1)^j*F_1(x[j]) }
      
      return(sum(m))
    }
    
  }
  
  ### Model2 = Interaction 
  
  if (Model==2){
    
    F_1=function(x){
      
      return(2*sin(pi*x))
    }
    
    if(sparse == F && !(p==30)){ sparsity = p-1}
    
    F_ursprung=function(x){
      
      m <- rep(0,sparsity)
      
      for(j in 1:sparsity){ m[j]=(-1)^j*F_1(x[j]*x[j+1]) }
      
      return(sum(m))
    }
    
  }
  
  ### Model3 = hierarchical
  
  if (Model==3){
    
    if(sparse == F && !(p==30)){ sparsity = p-1 }
    
    F_1=function(x){
      
      return(2*sin(pi*x))
    }
    
    F_ursprung=function(x){
      
      m <- m_2<- rep(0,(sparsity+1))
      
      for(j in 1:(sparsity+1)){ m[j]=(-1)^j*F_1(x[j]) }
      
      for(j in 1:sparsity){
        
        m_2[j]=(-1)^j*F_1(x[j]*x[j+1])
      }
      
      return(sum(m)+sum(m_2))
    }
    
  }
  
  ### Model4 = additive + Jump 
  
  if (Model==4){
    
    F_1=function(x){
      
      if(x>0){ return(2*sin(pi*x)-2) }
      if(x<0){ return(2*sin(pi*x)+2) }
      
      0
    }
    
    if(sparse == F && !(p==30)){ sparsity = p }
    
    F_ursprung=function(x){
      
      m <- rep(0,sparsity)
      
      for(j in 1:sparsity){ m[j]=(-1)^j*F_1(x[j]) }
      
      return(sum(m))
    }
    
  }
  
  ### Model5 = interaction + Jump 
  
  if (Model==5){
    
    F_1=function(x){
      
      if(x>0){ return(2*sin(pi*x)-2) }
      if(x<0){ return(2*sin(pi*x)+2) }
      
      0
    }
    
    if(sparse == F && !(p==30)){ sparsity = p-1 }
    
    F_ursprung=function(x){
      
      m <- rep(0,sparsity)
      
      for(j in 1:sparsity){ m[j]=(-1)^j*F_1(x[j]*x[j+1]) }
      
      return(sum(m))
    }
    
  }
  
  ### Model6 = hierarchical + Jump
  
  if (Model==6){
    
    F_1=function(x){
      
      if(x>0){ return(2*sin(pi*x)-2) }
      if(x<0){ return(2*sin(pi*x)+2) }
      
      0
    }
    
    if(sparse == F && !(p==30)){ sparsity = p-1 }
    
    F_ursprung=function(x){
      
      m <- m_2<-rep(0,sparsity)
      
      for(j in 1:sparsity){ m[j]=(-1)^j*F_1(x[j]*x[j+1]) }
      for(j in 1:(sparsity+1)){ m_2[j]=(-1)^j*2*sin(pi*x[j]) }
      
      return(sum(m)+sum(m_2))
    }
    
  }
  
  Y_true=0
  
  for(i in 1:nrow(X)){ Y_true[i]=F_ursprung(X[i,]) }
  
  Y_start=Y_true+epsilon
  
  return(list(Y_start=Y_start,Y_true=Y_true,X=X))
  
}