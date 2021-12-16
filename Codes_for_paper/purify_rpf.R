# Purify a model created by the random planted forest algorithm
# Input: forest_res = random planted forest model 
# Output: purified random planted forest model

rpf_purify=function(forest_res){  
  
  p=length(forest_res[,1]$intervals[[1]][[1]][1,])
  
  a=c(forest_res[,1]$intervals[[2]][[1]][1,1],forest_res[,1]$intervals[[1]][[1]][1,2:p])
  
  b=c(forest_res[,1]$intervals[[2]][[1]][2,1],forest_res[,1]$intervals[[1]][[1]][2,2:p])
  
  print("Note: The model is purified with respect to the hypercube")
  
  print(matrix(c(a,b), nrow=2, byrow=T))
  
  forest_res=rbind(forest_res, rep(0, dim(forest_res)[2]))
  
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