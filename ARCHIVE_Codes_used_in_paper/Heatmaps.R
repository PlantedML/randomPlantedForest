Features=c(1,2)
Gridsize=200
j=1

library(iml)
library(pdp)
library(xgboost)
library(Rcpp)

# load("Model4-K=30_FinalResults_plots.Rdata",envir=.GlobalEnv)


Data<-list(X=Results["X",j][[1]],Y_start=Results["Y",j][[1]])

p <- dim(Data$X)[2]

train.data <- data.frame(cbind(Data$Y_start,Data$X))
names(train.data)[1] <- "Y"
names(train.data)[2:(p+1)] <- paste0("V", 1:p)

a <- apply(Data$X,2,min)
b <- apply(Data$X,2,max)

model_PF_interaction=Results["model_PF_interaction",j][[1]]

Baum_mittel=50
F_Mittel1=model_PF_interaction

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

i=Features

x_PF1=seq(a[1],b[1], by=(b[Features[1]]-a[Features[1]])/(Gridsize-1))
x_PF2=seq(a[2],b[2], by=(b[Features[2]]-a[Features[2]])/(Gridsize-1))

# f_1=function(x,y){return(F_regression(c(x,y),i))}
f_1=function(x,y){return(F_average(c(x,y),i))}
f_2=function(x,y){return(mapply(f_1,x,y))}

z_PF=outer(x_PF1,x_PF2,f_2)

model_xgboost_interaction<- Results["model_xgboost_interaction",j][[1]]

X_boost=train.data[2:(p+1)]
model_boost=model_xgboost_interaction

names(X_boost)[1:p] <- paste0("V", 1:p)

pred <- function(model, newdata)  {
  results <- predict(model, as.matrix(newdata))
  return(results)
}

prededed=Predictor$new(model=model_boost, data=train.data, predict.function = pred, y="Y")

twodim=FeatureEffect$new(
  prededed,
  Features,
  method = "pdp",
  center.at = 0,
  grid.size = Gridsize
)

firstdim=FeatureEffect$new(
  prededed,
  Features[1],
  method = "pdp",
  center.at = 0,
  grid.size = Gridsize
)

seconddim=FeatureEffect$new(
  prededed,
  Features[2],
  method = "pdp",
  center.at = 0,
  grid.size = Gridsize
)

singledim=outer(firstdim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),2],seconddim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),2],function(x,y){return(x+y)})

twodimen=twodim$results[,3]

dim(twodimen)=c(Gridsize,Gridsize)
values=twodimen-singledim
pdf(file = "Model1-heatmap.pdf",  width = 9, height =3) 
par(mfrow=c(1,3),oma = c(0, 0, 3, 0))

TrueModel=outer(((-Gridsize):Gridsize)/Gridsize,
                ((-Gridsize):Gridsize)/Gridsize,
                function(x,y){return(-2*sin(pi*x*y))})

# image(firstdim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),1],
#       seconddim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),1],
#       TrueModel, 
#       col=grey.colors(10),
#       xlab=paste0("x_",as.character(Features[1])),
#       ylab=paste0("x_",as.character(Features[2])),
#       main=paste0("m_12(x_1,x_2)=-2 sin(","x_",as.character(Features[1]),"x_",as.character(Features[2]),")"))

# image(((-100):100)/100,
#       ((-100):100)/100,
#       TrueModel, 
#       col=grey.colors(10),
#       xlab=paste0("x_",as.character(Features[1])),
#       ylab=paste0("x_",as.character(Features[2])),
#       main=paste0("m_12(x_1,x_2)= -2sin(","x_",as.character(Features[1]),"x_",as.character(Features[2]),")"))

image(((-Gridsize):Gridsize)/Gridsize,
      ((-Gridsize):Gridsize)/Gridsize,
      TrueModel, 
      col=grey.colors(30),
      xlab=paste0("x_",as.character(Features[1])),
      ylab=paste0("x_",as.character(Features[2])),
      main="true function",
      breaks = -15:15/7,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5)

image(x_PF1,
      x_PF2,
      z_PF, 
      col=grey.colors(30),
      xlab=paste0("x_",as.character(Features[1])),
      ylab=paste0("x_",as.character(Features[2])),
      main="planted forest (interaction(2))",
      breaks = -15:15/7,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5)

image(firstdim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),1],
      seconddim$results[c(1:(Gridsize/2),(Gridsize/2+2):(Gridsize+1)),1],
      values, 
      col=grey.colors(30),
      xlab=paste0("x_",as.character(Features[1])),
      ylab=paste0("x_",as.character(Features[2])),
      main="xgboost (interaction(2))",
      breaks = -15:15/7,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5)

mtext("Model 3: hierarchical-interaction+sparse+jump, d=30", outer = TRUE, cex = 1.5)
dev.off()

