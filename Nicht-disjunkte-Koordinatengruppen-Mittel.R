TimeMittel<-Sys.time()

Y=Y_start

# Baum=Anzahl der B?ume f?r random Forests, K=Anzahl der potentiell maximierbaren Koordinaten

Schleife=function(durchlauf){
  
  K_active=0
  
  # Berechnen des bootstap sample
  
  subsample=sample(dim(X)[2],dim(X)[2],replace=TRUE)
  X_alt=X
  Y_alt=Y
  for(i in dim(X)[2]){
    X[,i]=X_alt[,subsample[i]]
    Y[i]=Y_alt[subsample[i]]
  }
  
  # x=X
  
  # Definition der verschiedenen Bl?cke
  
  F=list()
  
  # F[[1]][[i]][[k]]= is list of matrices describing intervals  for tree i, leaf k  (which tree tree i given in F3)
  
  F[[1]]=list()
  
  for(i in 1:length(K)){
    F[[1]][[i]]=list()
    F[[1]][[i]][[1]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
    for(j in 1:(dim(X)[1])){
      F[[1]][[i]][[1]][,j]=c(a[j],b[j])
    }
  }
  
  # Funktionswerte ben?tigt f?r das arithmetische Mittel
  
  # F[[2]][[i]] is a vector of predictions for each leaf in tree i
  
  F[[2]]=list()
  
  for(i in 1:length(K)){
    F[[2]][[i]]=0
  }
  
  # Koordinatengruppenindizees des i-ten Eintrags
  # F[[3]][[i]] is a vector of variables  in tree i
  
  F[[3]]=K
  
  
  # Definition der Partitionen
  # F[[4]][[i]][[k]] is a vector of individuals in leaf k of tree i
  
  F[[4]]=list()
  
  for(i in 1:length(K)){
    F[[4]][[i]]=list(1:(dim(X)[2]))
  }
 
 #  F5=? F6=? 
  
  
  F[[5]]=0
  
  F[[6]]=0
  
  # Residuen f?r das arithmetische Mittel
  
  ResiAverage <- function(X,Y,x,I,k){
    
    # X=Daten, Y=Aktuelles Y, x=Splitpunkt, I=Aktuelle Teilmenge der Indizes, k=Aktuelle Koordinate
    
    y_1=mean(Y[I][X[k,I]>=x])  ### mean Y of those in I with  X[k,I]>=x
    y_2=mean(Y[I][X[k,I]<x])   ### mean Y of those in I with  X[k,I]<x
    
    # y_1,y_2=Sch?tzer nach dem Split.
    
    Y[I][X[k,I]>=x]=Y[I][X[k,I]>=x]-y_1
    Y[I][X[k,I]<x]=Y[I][X[k,I]<x]-y_2
    
    return(sum(Y^2))
  }
  
  # Gegeben der aktuellen Koordinatenmenge K, wo ist der beste Split. 
  # Es wird das minimale Residuum sowie die Stelle des zugeh?rigen Splits zur?ckgegeben.
  
  # Bemerkung: Dies wird f?r endlich viele m?gliche Splitpunkte durchgef?hrt.
  
  SplitTest <- function(X,Y,x,I,K){
    
    # X=Daten, Y=aktuelles Y, x=Matrix der m?glichen Splitpunkte, I=Aktuelle Teilmenge der Indizes, 
    # K=Aktuelle Koordinaten.
    
    R=rep(max(Y^2)*length(Y)+1,dim(x)[2]*dim(x)[1])
    dim(R)=c(dim(x)[1],dim(x)[2])
    
    # R=Matrix der Residuen bzgl. der verschiedenen Splits.
    # Bemerkung: Die Residuen werden nur f?r die relevanten Stellen berechnet. Auf die anderen wird in return auch nicht zugegriffen.
    
    for(k in 1:dim(x)[1]){
      if(k %in% K){
        for(i in 1:(dim(x)[2])){
          if(min(c(sum(X[k,I]>=x[k,i]),sum(X[k,I]<x[k,i])))>=Blattgroesse[k]){
            R[k,i]=ResiAverage(X,Y,x[k,i],I,k)  ### residual sum of squares for split at x[k,i]
          }
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
  
  while(Splits>0){
    
    Splits=Splits-1
    
    K_use=sample(dim(X)[1],K_anzahl)
    
    # Hilfsvariablen
    
    R=1:length(F[[3]])
    Ikx=1:(3*length(F[[3]]))
    dim(Ikx)=c(3,length(F[[3]]))
    
    # R=Matrix der Residuen bzgl. des besten Splits f?r verschiedenen Koordinatengruppen K[[i]],
    # Ikx=Indizes der besten zu splitenden Menge, der zu splitenden Koordinate und des zugeh?rigen Splitpunkts f?r verschiedene Koordinatengruppen K.
    
    for(i_1 in 1:length(F[[3]])){
     
       #i_1 is index of current tree
      
      # Nur Teilsample verwenden
      
      K_neu=K_use
      
      if(length(F[[4]][[i_1]])==1 | length(F[[3]][[i_1]])==Komplex ){
        
        K_neu=intersect(K_use,F[[3]][[i_1]])
        
      }
      
      
      # Residuumsberechnung
      
      R_1=rep(sum(Y^2),length(F[[4]][[i_1]]))  # R_1=Matrix der Residuen des besten Splits f?r versch. Teilmengen der Datenmenge, 
      
     ## kx_1=1:(2*length(F[[4]][[i_1]]))   # kx_1=Index der zu splitenden Koordinate und der Splitpunkt. 
     
     #??????  dim(kx_1)=c(2,length(F[[4]][[i_1]]))
      
      kx_1<-matrix(nrow=2, ncol=length(F[[4]][[i_1]]))
      
      if(length(K_neu)>0){
        
       
        for(i_2 in 1:length(F[[4]][[i_1]])){
          
          #i_2 is index of current leaf
          
          # Berechnung der Splitpunkte.
          
          if(IterSplit == 1){
            for(i_3 in 1:(dim(X)[1])){
              if(i_3 %in% F[[3]][[i_1]]){
                h=1:m
                x[i_3,]=quantile(X[i_3,F[[4]][[i_1]][[i_2]]], h/m)
              }
            }
          }
          
          h=SplitTest(X,Y,x,F[[4]][[i_1]][[i_2]],K_neu)
          if(h[1]<sum(Y^2)){
            R_1[i_2]=h[1]
            kx_1[,i_2]=c(h[2],x[h[2],h[3]])
          }
        }
        
      }
      
      # Speicherung des Wertes des optimalen Residuums sowie dessen Index, die zu splitende Koordinate und den Splitpunkt
      
      R[i_1]=min(R_1)
      Ikx[,i_1]=c(which.min(R_1),kx_1[,which.min(R_1)])
    }
    
    # Mindestabbruchbedingung
    
    if(min(R) < sum(Y^2)){
      
      # Definition der optimalen Werte
      
      i_opt=which.min(R)
      Ikx_opt=Ikx[,i_opt]
      
      F[[6]]=c(F[[6]],c(i_opt,Ikx_opt[2]))
      
      # Updaten der Terme
      
      if(Ikx_opt[2] %in% F[[3]][[i_opt]]){
        
        F[[4]][[i_opt]][[length(F[[4]][[i_opt]])+1]]=intersect(which(X[Ikx_opt[2],]<Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
        
        # Updaten der mehrdimensionalen Intervalle
        
        F[[1]][[i_opt]][[length(F[[4]][[i_opt]])]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
        for(j in 1:(dim(X)[1])){
          F[[1]][[i_opt]][[length(F[[4]][[i_opt]])]][,j]=F[[1]][[i_opt]][[Ikx_opt[1]]][,j]
        }
        F[[1]][[i_opt]][[length(F[[4]][[i_opt]])]][,Ikx_opt[2]]=c(F[[1]][[i_opt]][[Ikx_opt[1]]][1,Ikx_opt[2]],Ikx_opt[3])
        
        # Berechnung des Differenzterms
        
        y_2=mean(Y[F[[4]][[i_opt]][[length(F[[4]][[i_opt]])]]])
        
        # Updaten der Y
        
        Y[F[[4]][[i_opt]][[length(F[[4]][[i_opt]])]]]=Y[F[[4]][[i_opt]][[length(F[[4]][[i_opt]])]]]-y_2
        
        # Updaten des zugeh?rigen Funktionswerts
        
        F[[2]][[i_opt]][length(F[[4]][[i_opt]])]=F[[2]][[i_opt]][Ikx_opt[1]]+y_2
        
        # Updaten der mehrdimensionalen Intervalle
        # Bemerkung: Das updaten der mehrdimensionalen Intervalle geschieht nur, falls auch wirklich das Datenset geteilt wird.
        
        F[[1]][[i_opt]][[Ikx_opt[1]]][1,Ikx_opt[2]]=Ikx_opt[3]
        
        # Updaten der I    
        
        F[[4]][[i_opt]][[Ikx_opt[1]]]=intersect(which(X[Ikx_opt[2],]>=Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
        
        # Berechnung des Differenzterms
        
        y_1=mean(Y[F[[4]][[i_opt]][[Ikx_opt[1]]]])
        
        # Updaten der Y
        
        Y[F[[4]][[i_opt]][[Ikx_opt[1]]]]=Y[F[[4]][[i_opt]][[Ikx_opt[1]]]]-y_1
        
        # Updaten des zugeh?rigen Funktionswerts
        
        F[[2]][[i_opt]][Ikx_opt[1]]=F[[2]][[i_opt]][Ikx_opt[1]]+y_1
        
      } else {
        
        TreeExists=0
        
        for(i in 1:length(F[[3]])){
          
          if(length(F[[3]][[i]])==length(F[[3]][[i_opt]])+1){
            
            if(prod(F[[3]][[i]]==sort(c(F[[3]][[i_opt]],Ikx_opt[2])))){
              
              TreeExists=1
              
              F[[4]][[i]][[length(F[[4]][[i]])+1]]=intersect(which(X[Ikx_opt[2],]<Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
              
              # Updaten der mehrdimensionalen Intervalle
              
              F[[1]][[i]][[length(F[[4]][[i]])]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
              for(j in 1:(dim(X)[1])){
                F[[1]][[i]][[length(F[[4]][[i]])]][,j]=F[[1]][[i_opt]][[Ikx_opt[1]]][,j]
              }
              F[[1]][[i]][[length(F[[4]][[i]])]][,Ikx_opt[2]]=c(F[[1]][[i_opt]][[Ikx_opt[1]]][1,Ikx_opt[2]],Ikx_opt[3])
              
              # Berechnung des Differenzterms
              
              y_2=mean(Y[F[[4]][[i]][[length(F[[4]][[i]])]]])
              
              # Updaten der Y
              
              Y[F[[4]][[i]][[length(F[[4]][[i]])]]]=Y[F[[4]][[i]][[length(F[[4]][[i]])]]]-y_2
              
              # Updaten des zugeh?rigen Funktionswerts
              
              F[[2]][[i]][length(F[[4]][[i]])]=y_2
              
              # Updaten der mehrdimensionalen Intervalle
              # Bemerkung: Das updaten der mehrdimensionalen Intervalle geschieht nur, falls auch wirklich das Datenset geteilt wird.
              
              F[[4]][[i]][[length(F[[4]][[i]])+1]]=intersect(which(X[Ikx_opt[2],]>=Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
              
              # Updaten der mehrdimensionalen Intervalle
              
              F[[1]][[i]][[length(F[[4]][[i]])]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
              for(j in 1:(dim(X)[1])){
                F[[1]][[i]][[length(F[[4]][[i]])]][,j]=F[[1]][[i_opt]][[Ikx_opt[1]]][,j]
              }
              F[[1]][[i]][[length(F[[4]][[i]])]][,Ikx_opt[2]]=c(Ikx_opt[3],F[[1]][[i_opt]][[Ikx_opt[1]]][2,Ikx_opt[2]])
              
              # Berechnung des Differenzterms
              
              y_1=mean(Y[F[[4]][[i]][[length(F[[4]][[i]])]]])
              
              # Updaten der Y
              
              Y[F[[4]][[i]][[length(F[[4]][[i]])]]]=Y[F[[4]][[i]][[length(F[[4]][[i]])]]]-y_1
              
              # Updaten des zugeh?rigen Funktionswerts
              
              F[[2]][[i]][length(F[[4]][[i]])]=y_1
              
            }
          }
        }
        
        if(TreeExists==0){
          
          F[[3]][[length(F[[3]])+1]]=sort(c(F[[3]][[i_opt]],Ikx_opt[2]))
          F[[4]][[length(F[[3]])]]=list()
          F[[4]][[length(F[[3]])]][[1]]=intersect(which(X[Ikx_opt[2],]<Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
          
          # Updaten der mehrdimensionalen Intervalle
          
          F[[1]][[length(F[[3]])]]=list()
          F[[1]][[length(F[[3]])]][[1]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
          for(j in 1:(dim(X)[1])){
            F[[1]][[length(F[[3]])]][[1]][,j]=F[[1]][[i_opt]][[Ikx_opt[1]]][,j]
          }
          F[[1]][[length(F[[3]])]][[1]][,Ikx_opt[2]]=c(F[[1]][[i_opt]][[Ikx_opt[1]]][1,Ikx_opt[2]],Ikx_opt[3])
          
          # Berechnung des Differenzterms
          
          y_2=mean(Y[F[[4]][[length(F[[3]])]][[1]]])
          
          # Updaten der Y
          
          Y[F[[4]][[length(F[[3]])]][[1]]]=Y[F[[4]][[length(F[[3]])]][[1]]]-y_2
          
          # Updaten des zugeh?rigen Funktionswerts
          
          F[[2]][[length(F[[3]])]]=y_2
          
          # Updaten der mehrdimensionalen Intervalle
          # Bemerkung: Das updaten der mehrdimensionalen Intervalle geschieht nur, falls auch wirklich das Datenset geteilt wird.
          
          F[[4]][[length(F[[3]])]][[2]]=intersect(which(X[Ikx_opt[2],]>=Ikx_opt[3]), F[[4]][[i_opt]][[Ikx_opt[1]]])
          
          # Updaten der mehrdimensionalen Intervalle
          
          F[[1]][[length(F[[3]])]][[2]]=matrix(1:(dim(X)[1]*2), nrow=2, ncol=dim(X)[1])
          for(j in 1:(dim(X)[1])){
            F[[1]][[length(F[[3]])]][[2]][,j]=F[[1]][[i_opt]][[Ikx_opt[1]]][,j]
          }
          F[[1]][[length(F[[3]])]][[2]][,Ikx_opt[2]]=c(Ikx_opt[3],F[[1]][[i_opt]][[Ikx_opt[1]]][2,Ikx_opt[2]])
          
          # Berechnung des Differenzterms
          
          y_1=mean(Y[F[[4]][[length(F[[3]])]][[2]]])
          
          # Updaten der Y
          
          Y[F[[4]][[length(F[[3]])]][[2]]]=Y[F[[4]][[length(F[[3]])]][[2]]]-y_1
          
          # Updaten des zugeh?rigen Funktionswerts
          
          F[[2]][[length(F[[3]])]][2]=y_1
          
        }
        
      }
      
    }
  }
  
  # Normieren der Komponenten
  
  F[[5]]=0
  
  # Purity Algorithmus
  
  for(l in Komplex:1){
    
    for(i_0 in 1:length(F[[3]])){
      
      if(length(F[[3]][[i_0]])==l){
        
        if(l==1){
          
          i_1=F[[3]][[i_0]]
          
          F[[1]][[i_0]][[length(F[[1]][[i_0]])+1]]=F[[1]][[i_0]][[1]]
          
          F[[1]][[i_0]][[length(F[[1]][[i_0]])]][,i_1]=c(a[i_1],b[i_1])
          
          F[[2]][[i_0]][[length(F[[1]][[i_0]])]]=0
        
          for(i_3 in 1:(length(F[[1]][[i_0]])-1)){
            
            F[[2]][[i_0]][[length(F[[1]][[i_0]])]]=F[[2]][[i_0]][[length(F[[1]][[i_0]])]]-F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
            
            F[[5]]=F[[5]]+F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
            
          }
          
        } else {
        
          for(i_1 in 1:dim(X)[1]){
            
            if(i_1 %in% F[[3]][[i_0]]){
              
              exists=0
              
              for(i_2 in 1:length(F[[3]])){
                
                if(length(F[[3]][[i_2]])==l-1){
                  
                  if(all(F[[3]][[i_2]]==F[[3]][[i_0]][F[[3]][[i_0]]!=i_1])){
                    
                    exists=1
                    
                    for(i_3 in 1:length(F[[1]][[i_0]])){
                      
                      exists_2=0
                      
                      F[[1]][[i_0]][[length(F[[1]][[i_0]])+1]]=F[[1]][[i_0]][[i_3]]
                      
                      F[[1]][[i_0]][[length(F[[1]][[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                      
                      F[[2]][[i_0]][[length(F[[1]][[i_0]])]]=-F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                      
                      for(i_4 in 1:length(F[[1]][[i_2]])){
                        
                        Cube=F[[1]][[i_0]][[i_3]]
                        
                        Cube[,i_1]=c(a[i_1],b[i_1])
                        
                        if(all(F[[1]][[i_2]][[i_4]]==Cube) & exists_2==0){
                          
                          exists_2=1
                          
                          F[[2]][[i_2]][[i_4]]=F[[2]][[i_2]][[i_4]]+F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                        } 
                      }
                      
                      if(exists_2==0){
                        
                        F[[1]][[i_2]][[length(F[[1]][[i_2]])+1]]=F[[1]][[i_0]][[length(F[[1]][[i_0]])]]
                        
                        F[[2]][[i_2]][[length(F[[1]][[i_2]])]]=F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                      }
                    }
                  }
                }
              }
              
              if(exists==0){
                
                F[[3]][[length(F[[3]])+1]]=F[[3]][[i_0]][F[[3]][[i_0]]!=i_1]
                
                F[[1]][[length(F[[3]])]]=list()
                
                F[[2]][[length(F[[3]])]]=vector()
                
                for(i_3 in 1:length(F[[1]][[i_0]])){
                  
                  F[[1]][[i_0]][[length(F[[1]][[i_0]])+1]]=F[[1]][[i_0]][[i_3]]
                  
                  F[[1]][[i_0]][[length(F[[1]][[i_0]])]][,i_1]=c(a[i_1],b[i_1])
                  
                  F[[2]][[i_0]][[length(F[[1]][[i_0]])]]=-F[[2]][[i_0]][[i_3]]*(F[[1]][[i_0]][[i_3]][2,i_1]-F[[1]][[i_0]][[i_3]][1,i_1])/(b[i_1]-a[i_1])
                  
                  F[[1]][[length(F[[3]])]][[length(F[[1]][[length(F[[3]])]])+1]]=F[[1]][[i_0]][[length(F[[1]][[i_0]])]]
                  
                  F[[2]][[length(F[[3]])]][[length(F[[1]][[length(F[[3]])]])]]=-F[[2]][[i_0]][[length(F[[1]][[i_0]])]]
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(F)
}

# Parallelisieren

# Erstelle ein Cluster

cl=makeCluster(Kerne)

clusterExport(cl,c("X","Y","K","a","b","IterSplit","x","m","K_anzahl","Blattgroesse","Komplex","Splits"))

# Berechne die n?tigen Informationen f?r den Sch?tzer

F=parSapply(cl, 1:Baum, Schleife)

# Cluster wieder entfernen

stopCluster(cl)

# Definition der resultierenden Funktionen 

# Definition der Funktion eines Baumes s

F_tree_Mittel=function(x,s=1,i=0,F){
  
  # x=Eingabevektor, s=Index des Baumes, i=Index des Summaden, den man berechnen m?chte. Bei i=0 wird die gesammte Funktion berechnet 
  
  f=0
    
  #Fall 1: Man berechnet f(x).
  
  if(length(i)==1){
    
    if(i==0){
      
      if(length(x) != dim(X)[1]){
        
        print("The vector x has the wrong dimension in order to calculate f(x)")  
      }
      
      for(i in 1:length(F[[3,s]])){
        
        for(j in 1:length(F[[1,s]][[i]])){
          
          if(prod(F[[1,s]][[i]][[j]][1,F[[3,s]][[i]]]<= x[F[[3,s]][[i]]] & (F[[1,s]][[i]][[j]][2,F[[3,s]][[i]]]>x[F[[3,s]][[i]]]|F[[1,s]][[i]][[j]][2,F[[3,s]][[i]]]==b[F[[3,s]][[i]]]))){
            
            f=f+F[[2,s]][[i]][j]
          }   
        }
      }
      
      return(f+F[[5,s]])   
    }
  }

  # Fall 1: Man berechnet den Wert des i-ten Summanden f_i(x).
  
  for(i_0 in 1:length(F[[3,s]])){
    
    if(length(F[[3,s]][[i_0]])==length(i)){
      
      if(prod(F[[3,s]][[i_0]]==sort(i))){
        
        if(length(x) != length(F[[3,s]][[i_0]])){
    
          print("The vector x has the wrong dimension in order to calculate f_i(x)")  
        }
        
        for(j in 1:length(F[[1,s]][[i_0]])){
          
          if(prod(F[[1,s]][[i_0]][[j]][1,F[[3,s]][[i_0]]]<= x & (F[[1,s]][[i_0]][[j]][2,F[[3,s]][[i_0]]]>x|F[[1,s]][[i_0]][[j]][2,F[[3,s]][[i_0]]]==b[F[[3,s]][[i_0]]]))){
            
            f=f+F[[2,s]][[i_0]][j]
          }   
        }
      }
    }
  }
  
  return(f)
}

# Definition der Funktion des Waldes 

F_Mittel1=F

Baum_mittel=Baum

F_average=function(x,i=0){
  
  y=0
  for(s in 1:Baum_mittel){
    y=y+F_tree_Mittel(x,s,i,F_Mittel1)
  }
  return(y/Baum_mittel)
}

# Berechnung des MSEs als Indikator f?r die Performance

S=0
for(j in 1:(dim(X)[2])){
  S[j]=F_average(X[,j])
}

MSE_average=1/n*sum((Y_true-S)^2)

TimeMittel=Sys.time()-TimeMittel