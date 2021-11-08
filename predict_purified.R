##### make prediction

pred_pur<-function(x,pur_res, interpret=TRUE){
#sapply(1:nrow(x), function(i){
  
  categorical_variables=pur_res$categorical_variables
  
  if (is.null(dim(x)))  x<-t(as.matrix(x))
  x <-as.matrix(x)
 
  if(!interpret){ 
  pred <- apply(sapply(1:length(pur_res),
              function(s){
                l0 <- which(sapply(pur_res[[s]]$my.functions,"[[",1)==0)
                return(rowSums(as.matrix(sapply((1:length(pur_res[[s]]$values))[-l0], ### go through all non-zero compoennts
                                  function(k){
                                    my.var <-pur_res[[s]]$my.functions[[k]]  
                                    
                                    pos <- numeric(length(my.var))
                                    y <- sapply(1:nrow(x),function(i){
                                    for(kk in (1:length(pos))){
                                      
                                      if (!is.element(my.var[kk], categorical_variables)){
                                        bounds <-  unlist((pur_res[[s]]$lim_list)[[my.var[kk]]])
                                        bounds[length(bounds)]<- bounds[length(bounds)]
                                        pos[kk] <- which(x[i,my.var[kk]]<bounds)[1]-1
                                        
                                        if (is.na(pos[kk])) {print("warning: observation is out of training range; extrapolating") 
                                          pos[kk]<-length(bounds)-1} else if (pos[kk]==0) {print("warning: observation is out of training range; extrapolating") 
                                            pos[kk]<-1}
                                        
                                        
                                      } else {    leaves <-  pur_res[[s]]$lim_list[[my.var[kk]]]
                                      pos[kk] <- which(sapply(1:ncol(leaves), function(j) is.element(x[i,my.var[kk]],  leaves[,j])))
                                      if (length(pos[kk])==0) {print("warning: observation is out of training range; extrapolating") }
                                      }
                                      
                                    }
                                      return(pur_res[[s]]$values[[k]][t(as.matrix(pos))]) 
                                    })
                                    as.numeric(y)}
                )),na.rm=TRUE)+as.numeric(pur_res[[s]]$values[[l0]]))
                
              }
              
              
  ),1,mean)
  return(pred)
  
 # }
#)
  } 
  
if (interpret){
  
  my_components <- unique(unlist(sapply(1:length(pur_res), function(s)
    pur_res[[s]]$my.functions[sapply(1:length(pur_res[[s]]$my.functions), function(z) sum(abs(pur_res[[s]]$values[[z]]))>0)]
  ),recursive = FALSE)) 
  
  my_components <- my_components[order(sapply(my_components, function(x){ if (length(x)>1) {length(x)+10} else x   }      ))]
  
  var_names <- sapply(1:length(my_components), function(k){
    paste(
      paste(colnames(x)[as.numeric(my_components[[k]])], collapse = ","),"(",
      paste(as.numeric(my_components[[k]]), collapse = ","),")",sep = "")
  }
  )
  
  pred <- sapply(1:nrow(x),function(i){
    
  
    
    
    res1 <-  sapply(1:length(pur_res),
               function(s){
                 l0 <- which(sapply(pur_res[[s]]$my.functions,"[[",1)==0)
                 return(sapply(1:length(my_components ), ### go through all compoennts
                                                 function(k){
                                                   
                                                   
                                                   component <- which(sapply(1:length(pur_res[[s]]$my.functions), function(l) setequal(pur_res[[s]]$my.functions[[l]],my_components[[k]])))
                                                   
                                                   if (length(component)==0) return(0)
                                                   
                                                   my.var <-pur_res[[s]]$my.functions[[component]]  
                                                   
                                                   if(sum(my.var==0)) return(pur_res[[s]]$values[[component]])
                                                   
                                                   pos <- numeric(length(my.var))
                                                   
                                                   
                                                  {
                                                     for(kk in (1:length(pos))){
                                                       
                                                       if (!is.element(my.var[kk], categorical_variables)){
                                                         bounds <-  unlist((pur_res[[s]]$lim_list)[[my.var[kk]]])
                                                         bounds[length(bounds)]<- bounds[length(bounds)]
                                                         pos[kk] <- which(x[i,my.var[kk]]<bounds)[1]-1
                                                         
                                                         if (is.na(pos[kk])) {print("warning: observation is out of training range; extrapolating") 
                                                           pos[kk]<-length(bounds)-1} else if (pos[kk]==0) {print("warning: observation is out of training range; extrapolating") 
                                                             pos[kk]<-1}
                                                         
                                                         
                                                       } else {    leaves <-  pur_res[[s]]$lim_list[[my.var[kk]]]
                                                       pos[kk] <- which(sapply(1:ncol(leaves), function(j) is.element(x[i,my.var[kk]],  leaves[,j])))
                                                       if (length(pos[kk])==0) {print("warning: observation is out of training range; extrapolating") }
                                                       }
                                                       
                                                     }
                                                    y <- pur_res[[s]]$values[[component]][t(as.matrix(pos))]
                                                    return(y)
                                                   }
                                                   }
                 ))
               }
              )
    
    res1 <- apply(  res1  ,1,mean)
    
    res1 <- c(sum(res1,na.rm=TRUE), res1)
    
    names(res1)=c("Y.hat", var_names)
  

    return(res1)
  })
  
  return(t(pred))
  # }
  #)
} 
  
  
  
}
  
  



