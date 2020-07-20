planted_forest<- function(Y, X, max_interaction=2, m_try=3, t_try, Baum=50, splits=30, m=10, Itersplit=1,n.cores=NULL, single_tree_ignore_m_try=TRUE, tree_from_m_try=FALSE, variables=NULL, new_trees=TRUE, Blattgroesse=rep(1,p))
{
  
  library(parallel)
  # Baum= Anzahl der B?um-famililien im "Forest"
  # splits= Anzahl der Iterationsschritte.
  # m= number of possible split points
  # Blattgroesse= minimal leaf size
  TimeMittel<-Sys.time()
  
  p <- ncol(X)
    
  a <- apply(X,2,min)     ## lower bounds
  b <- apply(X,2,max)     ### upper bounds
  
  x=t(X[1:m,])
  for(i in 1:p){
    h=1:m
    x[i,]=a[i]+(b[i]-a[i])*h/m
  }
  
  Schleife <- function(run){
    
    #???K_active=0
    # Berechnen des bootstap sample
    
    subsample <- sample(nrow(X),nrow(X),replace=TRUE)
    
    X_alt <- X
    Y_alt <- Y
    
    for(i in nrow(X)){
      X[i,]=X_alt[subsample[i],]
      Y[i]=Y_alt[subsample[i]]
    }
    
    
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
    
    #  F5=? F6=? 
    
    F<-list()
    F[[5]]=0
    
    F[[6]]=0
    
    # Residuen f?r das arithmetische Mittel
    
    ResiAverage <- function(X,Y,x,I,k){
      
      # X=Daten, Y=Aktuelles Y, x=Splitpunkt, I=Aktuelle Teilmenge der Indizes, k=Aktuelle Koordinate
      
      y_1 <- mean(Y[I][X[I,k]>=x])  ### mean Y of those in I with  X[I,k]>=x
      y_2 <- mean(Y[I][X[I,k]<x])   ### mean Y of those in I with  X[I,k]<x
      
      # y_1,y_2=Sch?tzer nach dem Split.
      
      Y[I][X[I,k]>=x] <- Y[I][X[I,k]>=x]-y_1
      Y[I][X[I,k]<x] <- Y[I][X[I,k]<x]-y_2
      
      return(sum(Y^2))
    }
    
    # Gegeben der aktuellen Koordinatenmenge K, wo ist der beste Split. 
    # Es wird das minimale Residuum sowie die Stelle des zugeh?rigen splits zur?ckgegeben.
    
    # Bemerkung: Dies wird f?r endlich viele m?gliche Splitpunkte durchgef?hrt.
    
    SplitTest <- function(X,Y,x,I,K){
      
      # X=Daten, Y=aktuelles Y, x=Matrix der m?glichen Splitpunkte, I=Aktuelle Teilmenge der Indizes, 
      # K=Aktuelle Koordinaten.
      
      R <- matrix(max(Y^2)*length(Y)+1, nrow=nrow(x), ncol=ncol(x))
      
      # R=Matrix der Residuen bzgl. der verschiedenen splits.
      # Bemerkung: Die Residuen werden nur f?r die relevanten Stellen berechnet. Auf die anderen wird in return auch nicht zugegriffen.
      
      for(k in K){
        for(i in 1:ncol(x)){
          if(min(c(sum(X[I,k]>=x[k,i]),sum(X[I,k]<x[k,i])))>=Blattgroesse[k]){
            R[k,i]=ResiAverage(X,Y,x[k,i],I,k)  ### residual sum of squares for split at x[k,i]
          }
        }
      }
      
      
      # Kurzer Test, ob ein Fehler vorherrscht
      
      if(sum(R>max(Y^2)*length(Y)+1)>0){
        print("Residuen drehen durch in SplitTest")
      }
      
      return(c(min(R),which(R == min(R), arr.ind = TRUE)[1,]))
    }
    
    # Iteration
    
    while(splits>0){
      
      splits=splits-1
      
     
      k_use=sample(1:p,m_try)
      if (m_try==0) k_use <- 1:p 
      
      if(tree_from_m_try==FALSE){
      tree_use <- sample(1:length(values), min(t_try,length(values)))
      if (t_try ==0) tree_use <- 1:length(values)
      } else{
         tree_candidates <- sapply(1:length(variables), function(i) {sum(k_use %in% variables[[i]])>=1}) 
         #if (m_try==4) write(tree_candidates, file="bla.txt",append=TRUE)
        tree_use <- sample(which(tree_candidates==TRUE), min(t_try,sum(tree_candidates)))
        if (t_try ==0) tree_use <- 1:length(values)
      }
      
      
      
      
      
      
      # Hilfsvariablen
      
      R <- rep(Inf,length(variables))
      Ikx <- matrix(nrow=3,ncol=length(variables))
      
      # R=Matrix der Residuen bzgl. des besten splits f?r verschiedenen Koordinatengruppen K[[i]],
      # Ikx=Indizes der besten zu splitenden Menge, der zu splitenden Koordinate und des zugeh?rigen Splitpunkts f?r verschiedene Koordinatengruppen K.
      
      
      for(i_1 in tree_use){
        
        #i_1 is index of current tree
        
        # Nur Teilsample verwenden
        
        if(single_tree_ignore_m_try==TRUE & length(variables[[i_1]])==1) k_neu <- union(k_use,variables[[i_1]]) else  k_neu <- k_use
        
        if(length(individuals[[i_1]])==1 | length(variables[[i_1]])==max_interaction | new_trees==FALSE ){
          ### if tree has only one leaf or already max number of interactions or new trees are not allowed
          k_neu <- intersect(k_neu,variables[[i_1]])
        }
        
        
        # Residuumsberechnung
        R_1 <- rep(sum(Y^2),length(individuals[[i_1]]))  # R_1=Matrix der Residuen des besten splits f?r versch. Blaetter
        kx_1<-matrix(nrow=2, ncol=length(individuals[[i_1]])) # kx_1=Index der zu splitenden Koordinate und der Splitpunkt.
        
        
        
        if(length(k_neu)>0){
          
          for(i_2 in 1:length(individuals[[i_1]])){
            
            #i_2 is index of current leaf
            # Berechnung der Splitpunkte.
            #i_3 is variable (in leaf i_2 in tree i_1)
            
            if(Itersplit == 1){
              for(i_3 in  variables[[i_1]]){
                h <- 1:m
                x[i_3,] <- quantile(X[individuals[[i_1]][[i_2]], i_3], h/m)
              }
            }
            
            h <- SplitTest(X,Y,x,individuals[[i_1]][[i_2]],k_neu)
            if(h[1]<sum(Y^2)){
              R_1[i_2] <- h[1]
              kx_1[,i_2] <- c(h[2],x[h[2],h[3]])
            }
          }
          R[i_1] <- min(R_1)
          Ikx[,i_1] <- c(which.min(R_1),kx_1[,which.min(R_1)])  # Ikx: each column for different tree. first row=which leaf to split; second row=which which variable to split, third row= where(=value) to split the variable
          
        } 
        
        # Speicherung des Wertes des optimalen Residuums sowie dessen Index, die zu splitende Koordinate und den Splitpunkt
        if (i_1==tail(tree_use,1)){
          if ( sum(R==Inf)==length(R) ) splits <- splits+1
        }
        
      }
      
      
      # Definition der optimalen Werte
      
      R_opt   <- min(R)
      Ikx_opt <- c(which.min(R), Ikx[,which.min(R)]) # best tree, leaf, variable, split-valaue
      
      
      # Mindestabbruchbedingung
      
      if(min(R) < sum(Y^2)){
        
        F[[6]] <- c(F[[6]], c(Ikx_opt[1], Ikx_opt[3]) ) 
        
        # Updaten der Terme
        I_12 <-  which(X[,Ikx_opt[3]]<Ikx_opt[4])   #### individuals who fulfil split criteria 2
        I_11 <-  which(X[,Ikx_opt[3]]>=Ikx_opt[4])   #### individuals who fulfil split criteria 1
        I_1 <-   individuals[[Ikx_opt[1]]][[Ikx_opt[2]]]  ### individuals in split-leaf 
        
        if(Ikx_opt[3] %in% variables[[Ikx_opt[1]]]){  ### if split variable is already in tree to be split
          
          
          individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])+1]] <- intersect(I_12, I_1 )      #### individuals for  new leaf 
          individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- intersect( I_11 , I_1)                         #### new split-leaf = remaining individuals
          
          # Updaten der mehrdimensionalen Intervalle
          intervals[[Ikx_opt[1]]][[length(intervals[[Ikx_opt[1]]])+1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]] #add one leaf to intervals
          
          intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][2,Ikx_opt[3]] <- Ikx_opt[4]  ### new leaf has new interval at split variable: (split value= upper bopund, lower bound remains)
          intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][1,Ikx_opt[3]] <- Ikx_opt[4] ### split leaf has new interval at split variable: (split value= lower bound, upper bound remains)
          
          # Berechnung des Differenzterms
          y_2 <- mean(Y[individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]]])
          y_1 <- mean(Y[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]]])
          
          
          # Updaten der Y
          Y[individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]]] <- Y[individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]]]-y_2
          Y[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]]] <- Y[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]]]-y_1
          
          # Updaten des zugeh?rigen Funktionswerts
          values[[Ikx_opt[1]]][length(individuals[[Ikx_opt[1]]])] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_2
          values[[Ikx_opt[1]]][Ikx_opt[2]] <- values[[Ikx_opt[1]]][Ikx_opt[2]]+y_1
          
        } else {
          
          TreeExists=0
          
          for(i in 1:length(variables)){
            
            if(length(variables[[i]])==length(variables[[Ikx_opt[1]]])+1){ ###if number of variables in tree i is same as num of var in new split tree
              if(prod(variables[[i]]==sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3])))){ ### if variables in tree i are same as  var in new split tree
                
                TreeExists=1  ###  tree i is the same as new split tree
                
                individuals[[i]][[length(individuals[[i]])+1]] <- intersect(I_12, I_1)
                individuals[[i]][[length(individuals[[i]])+1]] <- intersect(I_11, I_1)
                
                # Updaten der mehrdimensionalen Intervalle
                
                intervals[[i]][[length(individuals[[i]])]] <- intervals[[i]][[length(individuals[[i]])-1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]]
                
                intervals[[i]][[length(individuals[[i]])-1]][2,Ikx_opt[3]] <- Ikx_opt[4]
                intervals[[i]][[length(individuals[[i]])]][1,Ikx_opt[3]] <- Ikx_opt[4]
                
                # Berechnung des Differenzterms
                
                y_2 <- mean(Y[individuals[[i]][[length(individuals[[i]])-1]]])
                y_1 <- mean(Y[individuals[[i]][[length(individuals[[i]])]]])
                
                # Updaten der Y
                
                Y[individuals[[i]][[length(individuals[[i]])-1]]] <- Y[individuals[[i]][[length(individuals[[i]])-1]]]-y_2
                Y[individuals[[i]][[length(individuals[[i]])]]] <- Y[individuals[[i]][[length(individuals[[i]])]]]-y_1
                
                # Updaten des zugeh?rigen Funktionswerts
                
                values[[i]][length(individuals[[i]])-1]=y_2
                values[[i]][length(individuals[[i]])]=y_1
                
                
              }
            }
          }
          
          if(TreeExists==0){
            
            variables[[length(variables)+1]] <- sort(c(variables[[Ikx_opt[1]]],Ikx_opt[3]))
            individuals[[length(individuals)+1]] <- list()
            individuals[[length(individuals)]][[1]] <- intersect(which(X[,Ikx_opt[3]]<Ikx_opt[4]), individuals[[Ikx_opt[1]]][[Ikx_opt[2]]])
            individuals[[length(individuals)]][[2]] <- intersect(which(X[,Ikx_opt[3]]>=Ikx_opt[4]), individuals[[Ikx_opt[1]]][[Ikx_opt[2]]])
            
            # Updaten der mehrdimensionalen Intervalle
            
            intervals[[length(intervals)+1]] <- list()
            intervals[[length(intervals)]][[1]] <- intervals[[length(intervals)]][[2]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]]
            
            intervals[[length(intervals)]][[1]][2,Ikx_opt[3]] <- Ikx_opt[4]
            intervals[[length(intervals)]][[2]][1,Ikx_opt[3]] <- Ikx_opt[4]
            
            # Berechnung des Differenzterms
            
            y_2 <- mean(Y[individuals[[length(variables)]][[1]]])
            y_1 <- mean(Y[individuals[[length(variables)]][[2]]])
            
            # Updaten der Y
            
            Y[individuals[[length(variables)]][[1]]] <- Y[individuals[[length(variables)]][[1]]]-y_2
            Y[individuals[[length(variables)]][[2]]] <- Y[individuals[[length(variables)]][[2]]]-y_1
            
            values[[length(variables)]] <- y_2
            values[[length(variables)]][2] <- y_1
            
          }
          
        }
        
      }
    }
    
    # Normieren der Komponenten
    
    F[[5]]=0
    
    # Purity Algorithmus
    
    for(l in max_interaction:1){
      
      for(i_0 in 1:length(variables)){
        
        if(length(variables[[i_0]])==l){
          
          if(l==1){
            
            i_1=variables[[i_0]]
            
            intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[1]]
            
            intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
            
            values[[i_0]][[length(intervals[[i_0]])]]=0
            
            for(i_3 in 1:(length(intervals[[i_0]])-1)){
              
              values[[i_0]][[length(intervals[[i_0]])]]=values[[i_0]][[length(intervals[[i_0]])]]-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
              F[[5]]=F[[5]]+values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
              
            }
            
          } else {
            
            for(i_1 in 1:dim(X)[1]){
              
              if(i_1 %in% variables[[i_0]]){
                
                exists=0
                
                for(i_2 in 1:length(variables)){
                  
                  if(length(variables[[i_2]])==l-1){
                    
                    if(all(variables[[i_2]]==variables[[i_0]][variables[[i_0]]!=i_1])){
                      
                      exists=1
                      
                      for(i_3 in 1:length(intervals[[i_0]])){
                        
                        exists_2=0
                        
                        intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[i_3]]
                        
                        intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                        
                        values[[i_0]][[length(intervals[[i_0]])]]=-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        
                        for(i_4 in 1:length(intervals[[i_2]])){
                          
                          Cube=intervals[[i_0]][[i_3]]
                          
                          Cube[,i_1]=c(a[i_1],b[i_1])
                          
                          if(all(intervals[[i_2]][[i_4]]==Cube) & exists_2==0){
                            
                            exists_2=1
                            
                            values[[i_2]][[i_4]]=values[[i_2]][[i_4]]+values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                          } 
                        }
                        
                        if(exists_2==0){
                          
                          intervals[[i_2]][[length(intervals[[i_2]])+1]]=intervals[[i_0]][[length(intervals[[i_0]])]]
                          
                          values[[i_2]][[length(intervals[[i_2]])]]=values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        }
                      }
                    }
                  }
                }
                
                if(exists==0){
                  
                  variables[[length(variables)+1]]=variables[[i_0]][variables[[i_0]]!=i_1]
                  
                  intervals[[length(variables)]]=list()
                  
                  values[[length(variables)]]=vector()
                  
                  for(i_3 in 1:length(intervals[[i_0]])){
                    
                    intervals[[i_0]][[length(intervals[[i_0]])+1]]=intervals[[i_0]][[i_3]]
                    
                    intervals[[i_0]][[length(intervals[[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                    
                    values[[i_0]][[length(intervals[[i_0]])]]=-values[[i_0]][[i_3]]*(intervals[[i_0]][[i_3]][2,i_1]-intervals[[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                    
                    intervals[[length(variables)]][[length(intervals[[length(variables)]])+1]]=intervals[[i_0]][[length(intervals[[i_0]])]]
                    
                    values[[length(variables)]][[length(intervals[[length(variables)]])]]=-values[[i_0]][[length(intervals[[i_0]])]]
                    
                  }
                }
              }
            }
          }
        }
      }
    }
    
    return(list(intervals=intervals, values=values, variables=variables,  individuals=individuals,  F5=F[[5]], F6=F[[6]]))
  }
  
  # Parallelisieren
  # Erstelle ein Cluster
  if(is.null(n.cores)) Kerne <- detectCores() else Kerne <- n.cores
  cl <- makeCluster(Kerne)
  
  
  clusterExport(cl,c("Y","X","max_interaction", "m_try", "t_try", "Baum", "splits", "m", "Itersplit","a","b","p","x","Blattgroesse"), envir=environment())
  
  forest_res <- parSapply(cl, 1:Baum, Schleife) # Berechne die n?tigen Informationen f?r den Sch?tzer
  
  
  stopCluster(cl) # Cluster wieder entfernen
  
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
        
        return(f+forest_res[[5,s]])   
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
  
  TimeMittel=Sys.time()-TimeMittel
  return(list(Y_fitted=Y_fitted,forest_res=forest_res))
}

