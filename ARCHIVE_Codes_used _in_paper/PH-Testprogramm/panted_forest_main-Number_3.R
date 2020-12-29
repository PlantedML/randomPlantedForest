planted_forest_3<- function(Y, X, max_interaction=2, m_try=3, t_try=3, Baum=50, splits=30, m=10, Itersplit=1, RandomIterSplit=0, Itert_try=0, n.cores=NULL, single_tree_ignore_m_try=FALSE, start_ignore_m_try=TRUE, m_try_after_t_try=FALSE, only_splitable_tree=FALSE, split_try=0, Itersplit_try=0,  tree_from_m_try=FALSE, variables=NULL, new_trees=TRUE, Blattgroesse=rep(1,p))
{
  
  force(RandomIterSplit)
  force(Itersplit_try)
  
  # Baum= Anzahl der B?um-famililien im "Forest"
  # splits= Anzahl der Iterationsschritte.
  # m= number of possible split points
  # Blattgroesse= minimal leaf size
  
  p <- ncol(X)
  
  a <- apply(X,2,min)     ## lower bounds
  b <- apply(X,2,max)     ### upper bounds
  
  Schleife <- function(run){
    
    #???K_active=0
    # Berechnen des bootstap sample
    
    subsample <- sample(nrow(X),nrow(X),replace=TRUE)
    
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
    # intervals[[i]][[k]]= is list of matrices describing intervals  for tree i, leaf k  
    
    intervals=list()
    
    for(i in 1:length(variables)){
      intervals[[i]] <- list()
      intervals[[i]][[1]] <- matrix(nrow=2, ncol=p)
      for(j in 1:p){
        intervals[[i]][[1]][,j]=c(a[j],b[j])
      }
    }
    
    # Funktionswerte ben?tigt f?r das arithmetische Mittel
    
    # values[[i]] is a vector of predictions for each leaf in tree i
    
    values=list()
    
    for(i in 1:length(variables)){
      values[[i]]=0
    }
    
    # Definition der Partitionen
    # individuals[[i]][[k]] is a vector of individuals in leaf k of tree i
    
    individuals=list()
    
    for(i in 1:length(variables)){
      individuals[[i]]=list(1:nrow(X))
    }
    
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
      
      R_opt <- Inf
      
      # R=Matrix der Residuen bzgl. des besten splits f?r verschiedenen Koordinatengruppen K[[i]],
      # Ikx=Indizes der besten zu splitenden Menge, der zu splitenden Koordinate und des zugeh?rigen Splitpunkts f?r verschiedene Koordinatengruppen K.
      
      for(i_1 in 1:length(variables)){
        
        for(i_3 in 1:length(split_candidates)){
          
          k=split_candidates[[i_3]][[1]]
          
          if(length(unique(c(variables[[i_1]],k)))==length(split_candidates[[i_3]][[2]])){
            
            if(all(sort(unique(c(variables[[i_1]],k)))==split_candidates[[i_3]][[2]])){
              
              for(i_2 in 1:length(individuals[[i_1]])){
                
                I=individuals[[i_1]][[i_2]]
                
                for(i in 1:m){
                  
                  x <- sample(X[individuals[[i_1]][[i_2]], k], 1)
                  
                  if(min(c(sum(X[I,k]>=x),sum(X[I,k]<x)))>=Blattgroesse[k]){
                    
                    #R=ResiAverage(X,Y,x,I,k)  ### residual sum of squares for split at x[k,i]
                    
                    R=sum((Y[I][X[I,k]>=x]-mean(Y[I][X[I,k]>=x]))^2)+sum((Y[I][X[I,k]<x]-mean(Y[I][X[I,k]<x]))^2)+sum(Y[-I]^2)
                    
                    if(R<R_opt){
                      
                      R_opt <- R
                      Ikx_opt <- c(i_1,i_2,k,x)
                    }
                  }
                }
              }
            }
          }
        }
      }
      
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
      
      # Updaten der Y
      Y[I_2] <- Y[I_2]-y_2
      Y[I_1] <- Y[I_1]-y_1
      
      if(Ikx_opt[3] %in% variables[[Ikx_opt[1]]]){  ### if split variable is already in tree to be split
        
        individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])+1]] <- I_2      #### individuals for  new leaf 
        individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_1                       #### new split-leaf = remaining individuals
        
        # Updaten der mehrdimensionalen Intervalle
        intervals[[Ikx_opt[1]]][[length(intervals[[Ikx_opt[1]]])+1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]] #add one leaf to intervals
        
        intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][2,Ikx_opt[3]] <- Ikx_opt[4]  ### new leaf has new interval at split variable: (split value= upper bopund, lower bound remains)
        intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][1,Ikx_opt[3]] <- Ikx_opt[4] ### split leaf has new interval at split variable: (split value= lower bound, upper bound remains)
        
        values[[Ikx_opt[1]]][length(individuals[[Ikx_opt[1]]])] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_2
        values[[Ikx_opt[1]]][Ikx_opt[2]] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_1
        
      } else {
        
        TreeExists=0
        
        for(i in 1:length(variables)){
          
          if(length(variables[[i]])==length(variables[[Ikx_opt[1]]])+1){ ###if number of variables in tree i is same as num of var in new split tree
            if(all(variables[[i]]==sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3])))){ ### if variables in tree i are same as  var in new split tree
              
              TreeExists=1  ###  tree i is the same as new split tree
              
              individuals[[i]][[length(individuals[[i]])+1]] <- I_2
              individuals[[i]][[length(individuals[[i]])+1]] <- I_1
              
              # Updaten der mehrdimensionalen Intervalle
              
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
    
    return(list(intervals=intervals, values=values, variables=variables,  individuals=individuals))
  }
  
  # Parallelisieren
  # Erstelle ein Cluster
  if(is.null(n.cores)) Kerne <- detectCores() else Kerne <- n.cores
  
  if (Kerne!=1){
    cl <- makeCluster(Kerne)
    
    
    clusterExport(cl,c("Y","X","max_interaction", "m_try", "t_try", "Baum", "splits","m", "Itersplit", "Itert_try", "a","b","p","x","Blattgroesse"), envir=environment())
    
    forest_res <- parSapply(cl, 1:Baum, Schleife) # Berechne die n?tigen Informationen f?r den Sch?tzer
    
    
    stopCluster(cl) # Cluster wieder entfernen
  } else forest_res <- sapply(1:Baum, Schleife) 
  
  # Definition der resultierenden Funktionen 
  # Definition der Funktion eines Baumes s
  F_tree_Mittel=function(x,s=1,i=0,forest_res){
    
    # x=Eingabevektor, s=Index des Baumes, i=Index des Summaden, den man berechnen m?chte. Bei i=0 wird die gesammte   Funktion berechnet 
    
    f=0
    
    #Fall 1: Man berechnet f(x).
    
    if(length(i)==1){
      
      if(i==0){
        
        if(length(x) != p){
          
          print("The vector x has the wrong dimension in order to calculate f(x)")  
        }
        
        for(i in 1:length(forest_res[[3,s]])){
          
          for(j in 1:length(forest_res[[1,s]][[i]])){
            
            if(prod(forest_res[[1,s]][[i]][[j]][1,forest_res[[3,s]][[i]]]<= x[forest_res[[3,s]][[i]]] & (forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]>x[forest_res[[3,s]][[i]]]|forest_res[[1,s]][[i]][[j]][2,forest_res[[3,s]][[i]]]==b[forest_res[[3,s]][[i]]]))){
              
              f=f+forest_res[[2,s]][[i]][j]
            }   
          }
        }
        
        return(f)   
      }
    }
    
    # Fall 1: Man berechnet den Wert des i-ten Summanden f_i(x).
    
    for(i_0 in 1:length(forest_res[[3,s]])){
      
      if(length(forest_res[[3,s]][[i_0]])==length(i)){
        
        if(prod(forest_res[[3,s]][[i_0]]==sort(i))){
          
          if(length(x) != length(forest_res[[3,s]][[i_0]])){
            
            print("The vector x has the wrong dimension in order to calculate f_i(x)")  
          }
          
          for(j in 1:length(forest_res[[1,s]][[i_0]])){
            
            if(prod(forest_res[[1,s]][[i_0]][[j]][1,forest_res[[3,s]][[i_0]]]<= x & (forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]]]>x|forest_res[[1,s]][[i_0]][[j]][2,forest_res[[3,s]][[i_0]]]==b[forest_res[[3,s]][[i_0]]]))){
              
              f=f+forest_res[[2,s]][[i_0]][j]
            }   
          }
        }
      }
    }
    
    return(f)
  }
  
  # Definition der Funktion des Waldes 
  
  F_Mittel1 <- forest_res
  
  Baum_mittel <- Baum
  
  F_average=function(x,i=0){
    
    y=0
    for(s in 1:Baum_mittel){
      y=y+F_tree_Mittel(x,s,i,F_Mittel1)
    }
    return(y/Baum_mittel)
  }
  # Berechnung des MSEs als Indikator f?r die Performance
  
  Y_fitted=0
  for(j in 1:nrow(X)){
    Y_fitted[j]=F_average(X[j,])
  }
  
  #MSE_average=1/n*sum((Y_true-S)^2)
  
  return(list(Y_fitted=Y_fitted,forest_res=forest_res))
}