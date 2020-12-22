# Implementation of the random planted forest algorithm
# Input: (Y,X) = Data, max_interaction = maximum order of interaction, ntrees = Number of families of trees used, splits = number of splits,
#        split_try = number of considered splitpoints for each interval in each iteration step, t_try = percentage of viable trees considered in each iteration step,
#        variables = list of trees in each family in the beginning (default: one tree for each coordinate), 
#        leaf_size = minimum number of nodes in each leaf, alternative = alternative updating       
# Output: list of families of trees: [i][1] Final leaves of the trees, [i][2] = estimated values corresponding to the leaves, [i][3] = coordinates of the trees 
#                                     (for the i-th family, i=1,...,ntrees)

planted_forest<- function(Y, X, max_interaction=2, ntrees=50, splits=30, split_try=10, t_try=0.4, variables=NULL, leaf_size=rep(1,p), alternative=F){
  
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
    
    return(list(intervals=intervals, values=values, variables=variables, b=b))
  }
  
  forest_res <- sapply(1:ntrees, tree_fam) 
  
  return(forest_res)
}

# Purify a model created by the random planted forest algorithm
# Input: forest_res = random planted forest model 
# Output: purified random planted forest model

purify=function(forest_res){  
  
  p=length(forest_res[,1]$intervals[[1]][[1]][1,])
  
  a=c(forest_res[,1]$intervals[[2]][[1]][1,1],forest_res[,1]$intervals[[1]][[1]][1,2:p])
  
  b=c(forest_res[,1]$intervals[[2]][[1]][2,1],forest_res[,1]$intervals[[1]][[1]][2,2:p])
  
  print("Note: The model is purified with respect to the hypercube")
  
  print(matrix(c(a,b), nrow=2, byrow=T))
  
  forest_res=rbind(forest_res, rep(0,50))
  
  dimnames(forest_res)[[1]][5]="constant"
  
  for(s in 1:dim(forest_res)[2]){
    
    max_interaction=0
    
    for(i in 1:length(forest_res[,s]$variables)){ max_interaction=max(c(max_interaction,length(forest_res[,s]$variables[[i]]))) }
    
    for(l in max_interaction:1){
      
      for(i_0 in 1:length(forest_res[,s]$variables)){
        
        if(length(forest_res[,s]$variables[[i_0]])==l){
          
          if(l==1){
            
            i_1=forest_res[,s]$variables[[i_0]]
            
            forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])+1]]=forest_res[,s]$intervals[[i_0]][[1]]
            
            forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
            
            forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]=0
            
            for(i_3 in 1:(length(forest_res[,s]$intervals[[i_0]])-1)){
              
              # Here it becomes interesting
              
              forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]=forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]-forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
              forest_res[,s]$constant=forest_res[,s]$constant+forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
            }
            
          } else {
            
            for(i_1 in 1:p){
              
              if(i_1 %in% forest_res[,s]$variables[[i_0]]){
                
                exists=0
                
                for(i_2 in 1:length(forest_res[,s]$variables)){
                  
                  if(length(forest_res[,s]$variables[[i_2]])==l-1){
                    
                    if(all(forest_res[,s]$variables[[i_2]]==forest_res[,s]$variables[[i_0]][forest_res[,s]$variables[[i_0]]!=i_1])){
                      
                      exists=1
                      
                      for(i_3 in 1:length(forest_res[,s]$intervals[[i_0]])){
                        
                        exists_2=0
                        
                        forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])+1]]=forest_res[,s]$intervals[[i_0]][[i_3]]
                        
                        forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                        
                        #here 2
                        
                        forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]=-forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        
                        for(i_4 in 1:length(forest_res[,s]$intervals[[i_2]])){
                          
                          Cube=forest_res[,s]$intervals[[i_0]][[i_3]]
                          
                          Cube[,i_1]=c(a[i_1],b[i_1])
                          
                          if(all(forest_res[,s]$intervals[[i_2]][[i_4]]==Cube) & exists_2==0){
                            
                            exists_2=1
                            
                            forest_res[,s]$values[[i_2]][[i_4]]=forest_res[,s]$values[[i_2]][[i_4]]+forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                          } 
                        }
                        
                        if(exists_2==0){
                          
                          forest_res[,s]$intervals[[i_2]][[length(forest_res[,s]$intervals[[i_2]])+1]]=forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]
                          
                          forest_res[,s]$values[[i_2]][[length(forest_res[,s]$intervals[[i_2]])]]=forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        }
                      }
                    }
                  }
                }
                
                if(exists==0){
                  
                  forest_res[,s]$variables[[length(forest_res[,s]$variables)+1]]=forest_res[,s]$variables[[i_0]][forest_res[,s]$variables[[i_0]]!=i_1]
                  
                  forest_res[,s]$intervals[[length(forest_res[,s]$variables)]]=list()
                  
                  forest_res[,s]$values[[length(forest_res[,s]$variables)]]=vector()
                  
                  for(i_3 in 1:length(forest_res[,s]$intervals[[i_0]])){
                    
                    forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])+1]]=forest_res[,s]$intervals[[i_0]][[i_3]]
                    
                    forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                    
                    forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]=-forest_res[,s]$values[[i_0]][[i_3]]*(forest_res[,s]$intervals[[i_0]][[i_3]][2,i_1]-forest_res[,s]$intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                    
                    forest_res[,s]$intervals[[length(forest_res[,s]$variables)]][[length(forest_res[,s]$intervals[[length(forest_res[,s]$variables)]])+1]]=forest_res[,s]$intervals[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]
                    
                    forest_res[,s]$values[[length(forest_res[,s]$variables)]][[length(forest_res[,s]$intervals[[length(forest_res[,s]$variables)]])]]=-forest_res[,s]$values[[i_0]][[length(forest_res[,s]$intervals[[i_0]])]]
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(forest_res)
}

# Calculate the estimator for a value from a single family of trees
# Input: x = Input vector, s = Index of the family of trees, i = Coordinates of the component to be estimated, forest_res = planted forest model
# Output: Estimated function f(x) if i=0, otherwise estimated component f_i(x)

F_single_fam=function(x,s=1,i=0,forest_res){
  
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

F_average=function(x, forest_res, i=0){
  
  ntrees=dim(forest_res)[2]
  
  y=0
  
  for(s in 1:ntrees){ y=y+F_single_fam(x, s, i, forest_res) }
  
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

# plot the partial dependence plot of a random planted forest model (up to 2 dimensional) 
# Input: forest_res = random planted forest model, x = Vector of values creating the support of the plot, i = Coordinate to be plotted 
# Output: plot of a random planted forest model on the support x. If x=NULL the whole support is considered.

plot_rpf=function(forest_res, x=NULL, i=1, Gridsize = 30, type = "p",  xlim = NULL, ylim = NULL,
                  log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                  ann = par("ann"), axes = TRUE, frame.plot = axes,
                  panel.first = NULL, panel.last = NULL, asp = NA,
                  xgap.axis = NA, ygap.axis = NA, breaks = -15:15/7,
                  cex.axis = 1.5,cex.lab=1.5,cex.main=1.5){
  
  p=length(forest_res[,1]$intervals[[1]][[1]][1,])
  
  if(length(i)==2){
    
    if(max(i)>p){ return(paste("The coordinate i can only take on integers smaller then",p))  }
    
    if(!is.null(x)){ print("The X value is ignored when creating heatmaps!") }
      
    a=c(forest_res[,1]$intervals[[2]][[1]][1,1],forest_res[,1]$intervals[[1]][[1]][1,2:p])
      
    b=c(forest_res[,1]$intervals[[2]][[1]][2,1],forest_res[,1]$intervals[[1]][[1]][2,2:p])
    
    x=seq(a[i[1]],b[i[2]], length.out = Gridsize)
    y=seq(a[i[1]],b[i[2]], length.out = Gridsize)
    
    f_1=function(x,y){ return(predict_rpf(c(x,y), forest_res, i)) }
    f_2=function(x,y){ return(mapply(f_1,x,y)) }
    
    z=outer(x,y,f_2)
      
    image(x,
          y,
          z, 
          col=grey.colors(30),
          xlab=xlab,
          ylab=xlab,
          main=main,
          breaks = breaks,cex.axis = cex.axis,cex.lab=cex.lab,cex.main=cex.main)
  } else {
  
    if(length(i)!=1|i[1]>p){ return(paste("The coordinate i can only take on integers smaller then",p)) }
    
    if(is.null(x)){
      
      if(i==1){
        
        a=forest_res[,1]$intervals[[2]][[1]][1,1]
        
        b=forest_res[,1]$intervals[[2]][[1]][2,1]
      } else {
        
        a=forest_res[,1]$intervals[[1]][[1]][1,i]
        
        b=forest_res[,1]$intervals[[1]][[1]][2,i]
      }
      
      x=seq(a,b,length.out = 100)
    } else {
      
      x=sort(x)
    }
    
    dim(x)=c(length(x),1)
    
    plot(x,predict_rpf(x,forest_res,i=i), type = type,  xlim = xlim, ylim = ylim,
         log = log, main = main, sub = sub, xlab = xlab, ylab = ylab,
         ann = ann, axes = axes, frame.plot = frame.plot,
         panel.first = panel.first, panel.last = panel.last, asp = asp,
         xgap.axis = xgap.axis, ygap.axis = ygap.axis)
  }
}
