library(randomForest)
library(mgcv)
library(pdp)

F_1=function(x){
  return(-2*sin(pi*x))
}

F_2=function(x){
  if(x>0){ return(-(2*sin(pi*x)-2)) }
  if(x<0){ return(-(2*sin(pi*x)+2)) }
  0
}

F_3=function(x){
  if(x>0){ return(-(2*sin(pi*x)+1)) }
  if(x<0){ return(-(2*sin(pi*x)-1)) }
  0
}

F_inter <- function(x,y){
-2*sin(pi*x)-2*sin(pi*x*y)
}

Plots_Simulation <- function(feature=1, res=Results, num_Monte_Carlo=Monte_Carlo, TrueFunc=F_1, model1="PF_additive", model2=NULL, model3=NULL, plot_name="planted forest"){
  
  num_Monte_Carlo = min(Monte_Carlo, num_Monte_Carlo)
  
  feature<<-feature
  
  p<-dim(res["X",1][[1]])[2]
  
  col_function <- function(model){
    
    if(model==model1){
      return("grey")
    }
    
    if(model==model2){
      return("blue")
    }
    
    if(model==model3){
      return("yellow")
    }
  }
  
  #Generate empty plot
  
  YLim=c(-3,3)
  
  XLim=c(-1,1)
  
  plot(1, type="n", xlim=XLim, ylim=YLim , 
       main=plot_name, ylab="m_1", xlab=paste0("x_",feature),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  
  for(model in c(model1,model2,model3)){
    
    #Plots for PF
    
    if(model=="PF_alternative" | model=="PF_additive" | model=="PF_interaction" | model=="PF_fullinteraction"){
      
      hello=model
      
      Baum_mittel <- dim(res[paste0("model_",model),1][[1]])[2]
      
      for(j in 1:num_Monte_Carlo){
        
        a <- apply(res["X",j][[1]],2,min)
        b <- apply(res["X",j][[1]],2,max)
        
        F_Mittel1 <- res[paste0("model_",model),j][[1]]
        
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
        
        F_average=function(x,i=0){
          
          y=0
          for(s in 1:Baum_mittel){
            y=y+F_tree_Mittel(x,s,i,F_Mittel1)
          }
          return(y/Baum_mittel)
        }
        
        x=sort(res["X",j][[1]][,feature])
        
        F_average0=function(x){
          return(F_average(x,feature))
        }
        
        lines(x, sapply(x,F_average0), col=col_function(model))
      }
    }
    
    #Plots for backfit
    
    if(model=="backfit.LL" | model=="backfit.LC" | model=="backfit.BF"){
      
      for(j in 1:num_Monte_Carlo){
        
        model_backfit <- paste0("model_",model)
        lines(res[model_backfit,j][[1]]$x.grid[[feature]],res[model_backfit,j][[1]]$f_backfit[[feature]], col=col_function(model))
      }
    }
    
    if(model=="pspline"){
      
      for(j in 1:num_Monte_Carlo){
        
        model_pspline <- termplot(res["model_pspline",j][[1]], term=feature, plot=FALSE)
        lines(model_pspline[[1]][[1]], model_pspline[[1]][[2]], col=col_function(model))
      }
    }
    
    if(model=="xgboost_additive" | model=="xgboost_interaction"){
      
      hello="xgboost"
      
      for(j in 1:num_Monte_Carlo){
        
        X_boost=data.frame(res["X",j][[1]])
        model_boost=res[paste0("model_",model),j][[1]]
        
        names(X_boost)[1:p] <- paste0("V", 1:p)
        
        feature_xgboost=paste0("V",feature)
        
        fits=partial(model_boost, pred.var = feature_xgboost, ice = FALSE, center = FALSE, type="regression",
                    plot = FALSE, rug = FALSE, train=X_boost, grid.resolution = 200)
        x_1=fits[[1]]
        x_2=fits[[2]]
        lines(x_1,x_2, col=col_function(model))
      }
    }
    
    if(model=="randomForests"){
      
      for(j in 1:num_Monte_Carlo){
        
        x=res["model_randomForests",j][[1]]
        
        partialPlot(x,res[,j]$X, x.var= as.double(feature), add=T, col=col_function(model))
      }
    }
    
    if(model=="EBM_additive"){
      
      for(j in 1:num_Monte_Carlo){
        
        x=res["model_EBM_additive",j][[1]]
        lines(x$bin_edges[[feature]],x$model[[feature]][2:256],col=col_function(model))
      }
    }
    
    if(model=="gam"){
      
      hello="gam"
      
      for(j in 1:num_Monte_Carlo){
        
        model_gam=res["model_gam",j][[1]]
        par(new = T)
        plot(model_gam, se=F, select=feature, xlim=XLim, ylim=YLim ,  xlab=NA, ylab=NA, col=col_function(model),xaxt="n", yaxt="n")
      }
    }
  }
  
  #Plot TrueFunc
  

  if(is.function(TrueFunc)){
    
    x <- -100:100/100
    lines(x, sapply(x,TrueFunc), lwd=2)
    
  }
  
#  legend("topright", c("True Func", "hello", model2, model3), fill=c("black", "grey", "blue", "yellow"), cex=1.2)
  
}