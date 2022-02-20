planted_forest_2<- function(Y, X, max_interaction=2, Baum=50, splits=30, m=10, Itersplit_try=0.4, variables=NULL, Blattgroesse=rep(1,p), alternative=F)
{
  force(Itersplit_try)
  force(X)
  force(Y)
  
  # Baum= Anzahl der B?um-famililien im "Forest"
  # splits= Anzahl der Iterationsschritte.
  # m= number of possible split points
  # Blattgroesse= minimal leaf size
  
  p <- ncol(X)
  n <- nrow(X)
  
  Schleife <- function(run){
    
    X_orig=X
    
    #???K_active=0
    # Berechnen des bootstap sample
    
    subsample <- sample(n,n,replace=TRUE)
    
    X <- X[subsample,]
    Y <- Y[subsample]
    
    # Definition der verschiedenen Bl?cke
    
    # Koordinatengruppenindizees des i-ten Eintrags
    # variables[[i]] is a vector of variables  in tree i
    
    if (is.null(variables)){
      variables=list()
      for(i in 1:p){
        variables[[i]]=i
      }
    }
    
    # Definition der Partitionen
    # individuals[[i]][[k]] is a vector of individuals in leaf k of tree i
    
    individuals=list()
    
    for(i in 1:length(variables)){
      individuals[[i]]=list(1:n)
    }
    
    individuals_orig=list()
    
    for(i in 1:length(variables)){
      individuals_orig[[i]]=list(1:n)
    }
    
    Y_fitted=rep(0,n)
    
    Possible_Splits=list()
    for(i_1 in 1:length(variables)){
      for(i_2 in variables[[i_1]]){
        Possible_Splits[[i_1]]=list(i_2,variables[[i_1]])
      }
    }
    
    while(splits>0){
      
      splits=splits-1
      
      split_try=ceiling(Itersplit_try*length(Possible_Splits))
      
      split_candidates <- sample(Possible_Splits, split_try)
      
      # Hilfsvariablen
      
      # R=Matrix der Residuen bzgl. des besten splits f?r verschiedenen Koordinatengruppen K[[i]],
      # Ikx=Indizes der besten zu splitenden Menge, der zu splitenden Koordinate und des zugeh?rigen Splitpunkts f?r verschiedene Koordinatengruppen K.
      
      R=Calc_Optimal_split2(Y, X, m, variables, individuals, Blattgroesse, split_candidates)
      
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
        
        # Updaten der Y
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
              
              # Updaten der Y
              Y[individuals[[Ikx_opt[1]]][[j]]] <- Y[individuals[[Ikx_opt[1]]][[j]]]-mean(Y[individuals[[Ikx_opt[1]]][[j]]])
              
              Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]] <- Y_fitted[individuals_orig[[Ikx_opt[1]]][[j]]]-mean(Y[individuals_orig[[Ikx_opt[1]]][[j]]])
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
  
  forest_res <- sapply(1:Baum, Schleife) 
  
  Y_fitted=0
  
  for(i in 1:Baum){
    Y_fitted=Y_fitted+forest_res[,i]
  }
  
  Y_fitted=Y_fitted/Baum
  
  return(Y_fitted)
}