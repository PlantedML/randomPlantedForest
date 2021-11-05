#### input: values[[i]][[j]]        intervals[[i]][[j]]              variables[[i]][[j]] 
####      res[,s][[2]][[i]][[j]]    res[,s][[1]][[i]][[j]]    res[,s][[3]][[i]]


purify <- function(res,X,cores=1){

  X<-as.matrix(X)
values <- individuals <- list()

if (cores==1) purified <- lapply(1:ncol(res), function(s)
{
  
  ### lim_list is a list giving for each variable all interval end-points
  
  lim_list <- lapply(1:ncol(X),  ### go through all variables of the component
                     function(x){
                       my.var.tree.index2 <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) is.element(x,as.numeric(res[,s][[3]][[i]])))]  ## all trees that have x as variable
                       
                       if (sum(x==res[,s][[6]])==0){ ## not categorical variable
                         
                         bounds <- sort(unique( na.omit( as.numeric(unlist(sapply(1:length(my.var.tree.index2), ### go through all relevant trees
                                                                                  
                                                                                  function(i){ 
                                                                                    as.numeric(sapply(1:length( res[,s][[1]][[my.var.tree.index2[i]]]),
                                                                                                      function(j){ ### go through all leaves and get interval ends of variable
                                                                                                        res[,s][[1]][[my.var.tree.index2[i]]][[j]][,x]
                                                                                                      }))
                                                                                  }))))))
                         # bounds[length(bounds)] <- bounds[length(bounds)] +0.001 ### may need fix
                         return(bounds)
                       }else{ ## categorical variable
                         
                         temp <- lapply(my.var.tree.index2, ### go through all relevant trees
                                        
                                        function(i){ 
                                          sapply(1:length( res[,s][[1]][[i]]),
                                                 function(j){ ### go through all leaves and get their sets 
                                                   res[,s][[1]][[i]][[j]][,x]
                                                 })
                                        })
                         
                         
                         for(i in 1:length(temp )) { 
                           if (i==1) temp2 <- temp[[1]] else temp2 <- cbind(temp2,temp[[i]])
                         }
                         
                         temp <- temp2
                         
                         if (!is.null(ncol(temp2))){
                           rm(temp)
                           for (j in 1:ncol(temp2)){
                             if (j==1) temp <- matrix(temp2[,1]) else{
                               intercepts <-  sapply(1:ncol(temp), function(jj){
                                 sapply(1:nrow(temp2), function(k){
                                   any(temp2[k,j]==temp[,jj], na.rm=TRUE)
                                 })
                               })
                               
                               for(jj in 1:ncol(as.matrix(intercepts))){
                                 if(any(as.matrix(intercepts)[,jj])){
                                   #                             if(sum(as.matrix(intercepts)[,jj])!=length(unique(temp2[as.matrix(intercepts)[,jj],j]))){
                                   temp3 <- temp2[as.matrix(intercepts)[,jj],j]
                                   temp <- cbind(temp,c(temp3, rep(NA, nrow(temp2)-sum(as.matrix(intercepts)[,jj])) ))
                                   temp[t(sapply(1:length(temp3), function(jjj) which(temp3[jjj]==temp[,jj]))),jj]<-NA
                                   #}
                                 }
                               }
                               intercepts <- apply(intercepts, 1 , function(j) !any(j))
                               temp<-cbind(temp,c(temp2[intercepts,j],rep(NA,sum(!intercepts))))
                               temp<-as.matrix(temp[, colSums(is.na(temp)) != nrow(temp)])
                             }
                           }
                         }
                         
                         return(temp)
                         
                       }})
  
  my.functions <- unique(res[,s][[3]])  ### all function components
  
  
  
  ####
  ######### Get values and indiviudals
  
  for (l in 1:length( my.functions)){  ### go through all  function components
    
    my.var.tree.index <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) setequal(as.numeric(res[,s][[3]][[i]]),my.functions[[l]]))]  ### all trees equal function component
    
    
    
    
    
    ### Leaves are defined through lim_list
    #### valaues is a list giving for each leaf a value. 
    #### individuals is a list giving for each leaf the number of individuals. 
    
    if (l==1){values <- individuals<-list()}
    values[[l]] <- individuals[[l]] <- array(0, dim=sapply(1:length(my.functions[[l]]), function(m) {
      
      if (!any(my.functions[[l]][m]==res[,s][[6]])) {length(lim_list[[my.functions[[l]][m]]] )-1} else{
        ncol(as.matrix(lim_list[[my.functions[[l]][m]]]))
      }
      
    }
    
    ))
    
    
    
    
    ## x is index through all leaves
    x<-expand.grid(   lapply(1:length(my.functions[[l]]), function(j) 
      
      if (!any(my.functions[[l]][j]==res[,s][[6]])) {
        1:(length(lim_list[[my.functions[[l]][j]]])-1)} else{
          1:ncol(as.matrix(lim_list[[my.functions[[l]][j]]]))
        }
      
    )     ) 
    
    # if (is.null(dim(x))){
    # 
    #   for (a in 1:(length(x)-1))
    #   {
    #     
    #     for (i in my.var.tree.index){
    #       for(j in 1:length(res[,s][[2]][[i]])){ ### go through all leaves
    #       
    #         in_leave <- all(sapply(1:length(my.functions[[l]]),  function(k) {  
    #           res[,s][[1]][[i]][j][1,k]>=lim_list[[1]][a] & res[,s][[1]][[i]][j][2,k]<=lim_list[[1]][a+1] 
    #         }
    #         ))
    #         
    #     values[[s]][[l]][a] <- values[[s]][[l]][a] + res[,s][[2]][[i]][j] * in_leave
    #     individuals[[s]][[l]][a] <- individuals[[s]][[l]][a] + length(res[,s][[7]][[i]][j]) * in_leave
    #       
    #        
    #       } 
    #      } 
    #   }
    #   
    #   
    # }else
    {
      
      for (a in 1:(nrow(x)))
      {
        for (i in my.var.tree.index){
          for(j in 1:length(res[,s][[2]][[i]])){
            
            
            in_leave <- all(sapply(1:length(my.functions[[l]]),  function(k) {
              kk<-my.functions[[l]][k]
              if (!any(my.functions[[l]][k]==res[,s][[6]])){
                res[,s][[1]][[i]][[j]][1,my.functions[[l]][k]]<=lim_list[[kk]][x[a,k]] & res[,s][[1]][[i]][[j]][2,my.functions[[l]][k]]>=lim_list[[kk]][[x[a,k]+1]]
              }else{
                all(is.element(sort(lim_list[[kk]][,x[a,k]]),res[,s][[1]][[i]][[j]][,my.functions[[l]][k]])) ### any may be enough
              }
            }
            ))
            values[[l]][as.matrix((x[a,]))] <- values[[l]][as.matrix((x[a,]))] + res[,s][[2]][[i]][j] * in_leave 
            #individuals[[s]][[l]][as.matrix((x[a,]))] <- individuals[[s]][[l]][as.matrix((x[a,]))] + length(res[,s][[7]][[i]][j]) * in_leave
          } 
        }
        
        individuals[[l]][as.matrix((x[a,]))] <- sum(apply(
          sapply(1:length(my.functions[[l]]), function(k){
            kk <- my.functions[[l]][k]
            if (!any(my.functions[[l]][k]==res[,s][[6]])){
              #res[,s][[1]][[i]][j][1,k]>=lim_list[[k]][x[a,k]] & res[,s][[1]][[i]][j][2,k]<=lim_list[[k]][[x[a+1,k]]]
              up <- lim_list[[kk]][[x[a,k]+1]]
              #  if ((x[a,k]+1)==length(lim_list[[kk]])) up <- lim_list[[kk]][[x[a,k]+1]] + 0.01
              X[,my.functions[[l]][k]]>=lim_list[[kk]][x[a,k]] & X[,my.functions[[l]][k]]<up
            }else{
              is.element(X[,my.functions[[l]][k]],lim_list[[kk]][,x[a,k]]) 
            }
            
            
          }), 1, all))
      }
      
    }
    
    
  }
  
  
  
  
  #go through my.functions from top to bottom. 
  
  l_max = length(my.functions)  ### length of original components
  
  my.functions[[l_max +1]] <-  0   ### add constant function to the list of functions
  values[[l_max +1]] <-  array(0,dim=1)  
  
  
  #### getnew components 
  
  all_components_created <- 0
  while(  all_components_created == 0){
    
    for (l in length( my.functions):1){
      
      
      
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l_max +1
        
        if (length(l2)==0)  {
          
          my.functions[[length(my.functions)+1]] <- my.functions[[l]][-j] 
          values[[length(my.functions)]] <-  array(0, dim=dim(values[[l]])[-j])
          
          
        }
      }
    }
    
    all_components_created <- 1
    for (l in length( my.functions):1){
      
      
      
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l_max +1
        
        if (length(l2)==0)  {
          
          all_components_created <- 0
          break
          
        }
      }
    }
    
    
  }
  
  #### get individuals for leaves in new components 
  if ((l_max +2)<=length( my.functions)){
    for (l in  (l_max +2):length( my.functions)){
      
      
      x<-expand.grid( lapply(my.functions[[l]], function(j){
        if (!any(j==res[,s][[6]])) {
          1:(length(lim_list[[j]])-1)} else{
            1:ncol(as.matrix(lim_list[[j]]))
          }}
      )     
      )
      
      individuals[[l]] <- array(0, dim=sapply(1:length(my.functions[[l]]), function(m) {
        
        if (!any(my.functions[[l]][m]==res[,s][[6]])) {length(lim_list[[my.functions[[l]][m]]] )-1} else{
          ncol(as.matrix(lim_list[[my.functions[[l]][m]]]))
        }
        
      }
      ))
      
      for (a in 1:(nrow(x)))
      {
       # print(    dim(individuals[[l]]))
        individuals[[l]][as.matrix((x[a,]))] <- sum(apply(
          sapply(1:length(my.functions[[l]]), function(k){
            kk <- my.functions[[l]][k]
            if (!any(my.functions[[l]][k]==res[,s][[6]])){
              #res[,s][[1]][[i]][j][1,k]>=lim_list[[k]][x[a,k]] & res[,s][[1]][[i]][j][2,k]<=lim_list[[k]][[x[a+1,k]]]
              up <- lim_list[[kk]][[x[a,k]+1]]
              #  if ((x[a,k]+1)==length(lim_list[[kk]])) up <- lim_list[[kk]][[x[a,k]+1]] + 0.01
              X[,my.functions[[l]][k]]>=lim_list[[kk]][x[a,k]] & X[,my.functions[[l]][k]]<up
            }else{
              is.element(X[,my.functions[[l]][k]],lim_list[[kk]][,x[a,k]]) 
            }
            
            
          }), 1, all))
      }
      
    }
  }
  
  
  
  
  #################
  #######################         Purify
  
  
  ##### sort components
  values <- values[order(sapply(my.functions,function(x) length(x) ))]
  individuals <- individuals[order(sapply(my.functions,function(x) length(x) ))]
  
  my.functions <- my.functions[order(sapply(my.functions,function(x) length(x) ))]
  
  ### where is the constant function
  l0 <- which(sapply(my.functions,"[[",1)==0)
  
  
  
  #### start purify
  tol=rep(1,length( my.functions ))
  for (l in (length( my.functions):1)[(length( my.functions):1)!=l0]   ){ ## go through all functions
    
    #print(c("start",s,l))
    x<-expand.grid( lapply(my.functions[[l]], function(j){ ## x= all leaves in function
      if (!any(j==res[,s][[6]])) {
        1:(length(lim_list[[j]])-1)} else{
          1:ncol(as.matrix(lim_list[[j]]))
        }}
    )     
    )
    
    it=0
    while (tol[l]>0.00000000001 & it<100){
      if(it==99) print(tol[l])
      it<- it+1
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l0
        
        
        xj       <- as.matrix(unique( x[ ,-j ]))
        
        n.values <-  nrow(xj)
        n.values.j <- dim(values[[l]])[j]
        
        for(k in 1:n.values){
          
          n.values.j
          if (length(dim(values[[l]]))>1){
            if (j==1) slice <- t(sapply(1:n.values.j, function(z) c(z,xj[k,]) )) else if
            (j==length(dim(values[[l]]))){ slice <-t(sapply(1:n.values.j, function(z) c(xj[k,],z) ))} else{
              slice <- t(sapply(1:n.values.j, function(z) c(xj[k,1:(j-1)],z,xj[k,j:(ncol(xj))]) ))
            }
            
            if (is.null(dim(slice))) slice<-t(slice)
            
            #slice <-   x[ x[,j]==k,]
            if (sum(individuals[[l]][as.matrix(slice)])!=0){
              weight  <- individuals[[l]][as.matrix(slice)]/ sum(individuals[[l]][as.matrix(slice)])
              avg     <- t(weight) %*% values[[l]][as.matrix(slice)]
              values[[l]][as.matrix(slice)] <- values[[l]][as.matrix(slice)] - as.numeric(avg) 
              values[[l2]][t(as.matrix(xj[k,]))] <- values[[l2]][t(as.matrix(xj[k,]))] + as.numeric(avg) 
            }
          }
          
          if(length(dim(values[[l]]))==1) {
            weight <- individuals[[l]] / sum(individuals[[l]])
            avg   <- t(weight) %*% values[[l]]
            values[[l]]<- values[[l]] - as.numeric(avg) 
            values[[l0]] <- values[[l0]] + as.numeric(avg) 
          } 
        }
      } 
      
      if(length(my.functions[[l]])>1) tol[l] <- max(c(rowSums(values[[l]]*individuals[[l]]), colSums(values[[l]]*individuals[[l]])))
      if(length(my.functions[[l]])==1) tol[l] <- max(sum(values[[l]]*individuals[[l]]))
    }
  #  print(c("end",s,l))   
  }
  
  return(list(values=values, lim_list =lim_list ,individuals=individuals, my.functions= my.functions,categorical_variables=categorical_variables  ))
}) else {
purified <- mclapply(1:ncol(res), function(s)
{
  
  ### lim_list is a list giving for each variable all interval end-points
  
  lim_list <- lapply(1:ncol(X),  ### go through all variables of the component
                     function(x){
                       my.var.tree.index2 <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) is.element(x,as.numeric(res[,s][[3]][[i]])))]  ## all trees that have x as variable
                       
                       if (sum(x==res[,s][[6]])==0){ ## not categorical variable
                         
                         bounds <- sort(unique( na.omit( as.numeric(unlist(sapply(1:length(my.var.tree.index2), ### go through all relevant trees
                                                                                  
                                                                                  function(i){ 
                                                                                    as.numeric(sapply(1:length( res[,s][[1]][[my.var.tree.index2[i]]]),
                                                                                                      function(j){ ### go through all leaves and get interval ends of variable
                                                                                                        res[,s][[1]][[my.var.tree.index2[i]]][[j]][,x]
                                                                                                      }))
                                                                                  }))))))
                         # bounds[length(bounds)] <- bounds[length(bounds)] +0.001 ### may need fix
                         return(bounds)
                       }else{ ## categorical variable
                         
                         temp <- lapply(my.var.tree.index2, ### go through all relevant trees
                                        
                                        function(i){ 
                                          sapply(1:length( res[,s][[1]][[i]]),
                                                 function(j){ ### go through all leaves and get their sets 
                                                   res[,s][[1]][[i]][[j]][,x]
                                                 })
                                        })
                         
                         
                         for(i in 1:length(temp )) { 
                           if (i==1) temp2 <- temp[[1]] else temp2 <- cbind(temp2,temp[[i]])
                         }
                         
                         temp <- temp2
                         
                         if (!is.null(ncol(temp2))){
                           rm(temp)
                           for (j in 1:ncol(temp2)){
                             if (j==1) temp <- matrix(temp2[,1]) else{
                               intercepts <-  sapply(1:ncol(temp), function(jj){
                                 sapply(1:nrow(temp2), function(k){
                                   any(temp2[k,j]==temp[,jj], na.rm=TRUE)
                                 })
                               })
                               
                               for(jj in 1:ncol(as.matrix(intercepts))){
                                 if(any(as.matrix(intercepts)[,jj])){
                                   #                             if(sum(as.matrix(intercepts)[,jj])!=length(unique(temp2[as.matrix(intercepts)[,jj],j]))){
                                   temp3 <- temp2[as.matrix(intercepts)[,jj],j]
                                   temp <- cbind(temp,c(temp3, rep(NA, nrow(temp2)-sum(as.matrix(intercepts)[,jj])) ))
                                   temp[t(sapply(1:length(temp3), function(jjj) which(temp3[jjj]==temp[,jj]))),jj]<-NA
                                   #}
                                 }
                               }
                               intercepts <- apply(intercepts, 1 , function(j) !any(j))
                               temp<-cbind(temp,c(temp2[intercepts,j],rep(NA,sum(!intercepts))))
                               temp<-as.matrix(temp[, colSums(is.na(temp)) != nrow(temp)])
                             }
                           }
                         }
                         
                         return(temp)
                         
                       }})
  
  my.functions <- unique(res[,s][[3]])  ### all function components
  
  
  
  ####
  ######### Get values and indiviudals
  
  for (l in 1:length( my.functions)){  ### go through all  function components
    
    my.var.tree.index <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) setequal(as.numeric(res[,s][[3]][[i]]),my.functions[[l]]))]  ### all trees equal function component
    
    
    
    
    
    ### Leaves are defined through lim_list
    #### valaues is a list giving for each leaf a value. 
    #### individuals is a list giving for each leaf the number of individuals. 
    
    if (l==1){values <- individuals<-list()}
    values[[l]] <- individuals[[l]] <- array(0, dim=sapply(1:length(my.functions[[l]]), function(m) {
      
      if (!any(my.functions[[l]][m]==res[,s][[6]])) {length(lim_list[[my.functions[[l]][m]]] )-1} else{
        ncol(as.matrix(lim_list[[my.functions[[l]][m]]]))
      }
      
    }
    
    ))
    
    
    
    
    ## x is index through all leaves
    x<-expand.grid(   lapply(1:length(my.functions[[l]]), function(j) 
      
      if (!any(my.functions[[l]][j]==res[,s][[6]])) {
        1:(length(lim_list[[my.functions[[l]][j]]])-1)} else{
          1:ncol(as.matrix(lim_list[[my.functions[[l]][j]]]))
        }
      
    )     ) 
    
    # if (is.null(dim(x))){
    # 
    #   for (a in 1:(length(x)-1))
    #   {
    #     
    #     for (i in my.var.tree.index){
    #       for(j in 1:length(res[,s][[2]][[i]])){ ### go through all leaves
    #       
    #         in_leave <- all(sapply(1:length(my.functions[[l]]),  function(k) {  
    #           res[,s][[1]][[i]][j][1,k]>=lim_list[[1]][a] & res[,s][[1]][[i]][j][2,k]<=lim_list[[1]][a+1] 
    #         }
    #         ))
    #         
    #     values[[s]][[l]][a] <- values[[s]][[l]][a] + res[,s][[2]][[i]][j] * in_leave
    #     individuals[[s]][[l]][a] <- individuals[[s]][[l]][a] + length(res[,s][[7]][[i]][j]) * in_leave
    #       
    #        
    #       } 
    #      } 
    #   }
    #   
    #   
    # }else
    {
      
      for (a in 1:(nrow(x)))
      {
        for (i in my.var.tree.index){
          for(j in 1:length(res[,s][[2]][[i]])){
            
            
            in_leave <- all(sapply(1:length(my.functions[[l]]),  function(k) {
              kk<-my.functions[[l]][k]
              if (!any(my.functions[[l]][k]==res[,s][[6]])){
                res[,s][[1]][[i]][[j]][1,my.functions[[l]][k]]<=lim_list[[kk]][x[a,k]] & res[,s][[1]][[i]][[j]][2,my.functions[[l]][k]]>=lim_list[[kk]][[x[a,k]+1]]
              }else{
                all(is.element(sort(lim_list[[kk]][,x[a,k]]),res[,s][[1]][[i]][[j]][,my.functions[[l]][k]])) ### any may be enough
              }
            }
            ))
            values[[l]][as.matrix((x[a,]))] <- values[[l]][as.matrix((x[a,]))] + res[,s][[2]][[i]][j] * in_leave 
            #individuals[[s]][[l]][as.matrix((x[a,]))] <- individuals[[s]][[l]][as.matrix((x[a,]))] + length(res[,s][[7]][[i]][j]) * in_leave
          } 
        }
        
        individuals[[l]][as.matrix((x[a,]))] <- sum(apply(
          sapply(1:length(my.functions[[l]]), function(k){
            kk <- my.functions[[l]][k]
            if (!any(my.functions[[l]][k]==res[,s][[6]])){
              #res[,s][[1]][[i]][j][1,k]>=lim_list[[k]][x[a,k]] & res[,s][[1]][[i]][j][2,k]<=lim_list[[k]][[x[a+1,k]]]
              up <- lim_list[[kk]][[x[a,k]+1]]
            #  if ((x[a,k]+1)==length(lim_list[[kk]])) up <- lim_list[[kk]][[x[a,k]+1]] + 0.01
              X[,my.functions[[l]][k]]>=lim_list[[kk]][x[a,k]] & X[,my.functions[[l]][k]]<up
            }else{
              is.element(X[,my.functions[[l]][k]],lim_list[[kk]][,x[a,k]]) 
            }
            
            
          }), 1, all))
      }
      
    }
    
    
  }
  
  
  
  
  #go through my.functions from top to bottom. 
  
  l_max = length(my.functions)  ### length of original components
  
  my.functions[[l_max +1]] <-  0   ### add constant function to the list of functions
  values[[l_max +1]] <-  array(0,dim=1)  
  
  
  #### getnew components 
  
  all_components_created <- 0
  while(  all_components_created == 0){
    
    for (l in length( my.functions):1){
      
      
      
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l_max +1
        
        if (length(l2)==0)  {
          
          my.functions[[length(my.functions)+1]] <- my.functions[[l]][-j] 
          values[[length(my.functions)]] <-  array(0, dim=dim(values[[l]])[-j])
          
          
        }
      }
    }
    
    all_components_created <- 1
    for (l in length( my.functions):1){
      
      
      
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l_max +1
        
        if (length(l2)==0)  {
          
          all_components_created <- 0
          break
          
        }
      }
    }
    
    
  }
  
  #### get individuals for leaves in new components 
  if ((l_max +2)<=length( my.functions)){
    for (l in  (l_max +2):length( my.functions)){
      
      
      x<-expand.grid( lapply(my.functions[[l]], function(j){
        if (!any(j==res[,s][[6]])) {
          1:(length(lim_list[[j]])-1)} else{
            1:ncol(as.matrix(lim_list[[j]]))
          }}
      )     
      )
      
      individuals[[l]] <- array(0, dim=sapply(1:length(my.functions[[l]]), function(m) {
        
        if (!any(my.functions[[l]][m]==res[,s][[6]])) {length(lim_list[[my.functions[[l]][m]]] )-1} else{
          ncol(as.matrix(lim_list[[my.functions[[l]][m]]]))
        }
        
      }
      ))
      
      for (a in 1:(nrow(x)))
      {
        #print(    dim(individuals[[l]]))
        individuals[[l]][as.matrix((x[a,]))] <- sum(apply(
          sapply(1:length(my.functions[[l]]), function(k){
            kk <- my.functions[[l]][k]
            if (!any(my.functions[[l]][k]==res[,s][[6]])){
              #res[,s][[1]][[i]][j][1,k]>=lim_list[[k]][x[a,k]] & res[,s][[1]][[i]][j][2,k]<=lim_list[[k]][[x[a+1,k]]]
              up <- lim_list[[kk]][[x[a,k]+1]]
            #  if ((x[a,k]+1)==length(lim_list[[kk]])) up <- lim_list[[kk]][[x[a,k]+1]] + 0.01
              X[,my.functions[[l]][k]]>=lim_list[[kk]][x[a,k]] & X[,my.functions[[l]][k]]<up
            }else{
              is.element(X[,my.functions[[l]][k]],lim_list[[kk]][,x[a,k]]) 
            }
            
            
          }), 1, all))
      }
      
    }
  }
  
  
  
  
  #################
  #######################         Purify
  
  
  ##### sort components
  values <- values[order(sapply(my.functions,function(x) length(x) ))]
  individuals <- individuals[order(sapply(my.functions,function(x) length(x) ))]
  
  my.functions <- my.functions[order(sapply(my.functions,function(x) length(x) ))]
  
  ### where is the constant function
  l0 <- which(sapply(my.functions,"[[",1)==0)
  
  
  
  #### start purify
  tol=rep(1,length( my.functions ))
  for (l in (length( my.functions):1)[(length( my.functions):1)!=l0]   ){ ## go through all functions
    
    #print(c("start",s,l))
    x<-expand.grid( lapply(my.functions[[l]], function(j){ ## x= all leaves in function
      if (!any(j==res[,s][[6]])) {
        1:(length(lim_list[[j]])-1)} else{
          1:ncol(as.matrix(lim_list[[j]]))
        }}
    )     
    )
    
    it=0
    while (tol[l]>0.00000000001 & it<100){
      if(it==99) print(tol[l])
      it<- it+1
      for (j in (1:length(dim(values[[l]])))){ ## j=which variable to integrate over
        
        
        
        #l2 is tree that has same variables as l minus j-variable
        
        if (length(my.functions[[l]])!= 1){ 
          l2 <- which(sapply(1:length( my.functions ), function (z) 
          {length(my.functions[[z]])==length(my.functions[[l]][-j])&
              all(is.element(my.functions[[z]],my.functions[[l]][-j]))}
          )
          )
        } else l2 <- l0
        
        
        xj       <- as.matrix(unique( x[ ,-j ]))
        
        n.values <-  nrow(xj)
        n.values.j <- dim(values[[l]])[j]
        
        for(k in 1:n.values){
          
          n.values.j
          if (length(dim(values[[l]]))>1){
            if (j==1) slice <- t(sapply(1:n.values.j, function(z) c(z,xj[k,]) )) else if
            (j==length(dim(values[[l]]))){ slice <-t(sapply(1:n.values.j, function(z) c(xj[k,],z) ))} else{
              slice <- t(sapply(1:n.values.j, function(z) c(xj[k,1:(j-1)],z,xj[k,j:(ncol(xj))]) ))
            }
            
            if (is.null(dim(slice))) slice<-t(slice)
            
            #slice <-   x[ x[,j]==k,]
            if (sum(individuals[[l]][as.matrix(slice)])!=0){
              weight  <- individuals[[l]][as.matrix(slice)]/ sum(individuals[[l]][as.matrix(slice)])
              avg     <- t(weight) %*% values[[l]][as.matrix(slice)]
              values[[l]][as.matrix(slice)] <- values[[l]][as.matrix(slice)] - as.numeric(avg) 
              values[[l2]][t(as.matrix(xj[k,]))] <- values[[l2]][t(as.matrix(xj[k,]))] + as.numeric(avg) 
            }
          }
          
          if(length(dim(values[[l]]))==1) {
            weight <- individuals[[l]] / sum(individuals[[l]])
            avg   <- t(weight) %*% values[[l]]
            values[[l]]<- values[[l]] - as.numeric(avg) 
            values[[l0]] <- values[[l0]] + as.numeric(avg) 
          } 
        }
      } 
      
      if(length(my.functions[[l]])>1) tol[l] <- max(c(rowSums(values[[l]]*individuals[[l]]), colSums(values[[l]]*individuals[[l]])))
      if(length(my.functions[[l]])==1) tol[l] <- max(sum(values[[l]]*individuals[[l]]))
    }
    #print(c("end",s,l))   
  }
  
  return(list(values=values, lim_list =lim_list ,individuals=individuals, my.functions= my.functions,categorical_variables=categorical_variables  ))
}, mc.cores=cores)
}

return(purified)

}