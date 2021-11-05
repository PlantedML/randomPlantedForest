#auc(Y[test],as.numeric(res2[[1]][[1]]))

##### plot main effects

plot_pur <- function(pur_res,X, categorical_variables, k)
{
  library(RColorBrewer)
  cols = brewer.pal(4, "Blues")
  # Define colour pallete
  pal = colorRampPalette(c("blue", "red"))
  # Use the following line with RColorBrewer
  pal = colorRampPalette(cols)
  # Rank variable for colour assignment
  
  if (length(k)==1){
#### which component to plot
#k <- 5
################



if (!is.element(k, categorical_variables)) { 

x <- seq(min(X[,k]),max(X[,k]), length.out=200)
y<-sapply(1:length(x), function(i)
mean(sapply(1:length(pur_res), function(s){
bounds <-  unlist((pur_res[[s]]$lim_list)[[k]])
bounds[length(bounds)]<- bounds[length(bounds)]+0.000001
pos <- which(x[i]<bounds)[1]-1
(pur_res[[s]]$values)[[k]][pos]
}
)
)
)
plot(x,y, main=colnames(X)[k])
}


if (is.element(k, categorical_variables)) { 
  
  x <- unique(X[,k])
  y<-sapply(1:length(x), function(i)
    mean(sapply(1:length(pur_res), function(s){
      leaves <-  ((pur_res[[s]]$lim_list)[[k]])
      pos <- sapply(1:ncol(leaves), function(j) is.element(x[i],  leaves[,j]))
      if(sum(pos)==0) return(NA) else (pur_res[[s]]$values)[[k]][pos]
    }
    )
    ,na.rm=TRUE)
  )
  barplot(y, names=x, main=colnames(X)[k])
}


}

  if (length(k)==2){

######################
##### plot interaction effects
#############################
#my.interactions<- unique(unlist(sapply(1:length(pur_res), function(s)   unique(res[,s][[3]])),recursive = FALSE))[-(1:10)]
#k <- my.interactions[[1]]

#############
#k=c(7,10)  ## which component to plot
################


if (!is.element(k[1], categorical_variables))  x1 <- seq(min(X[,k[1]]),max(X[,k[1]]), length.out=50)
if (is.element(k[1], categorical_variables))    x1 <- unique(X[,k[1]])  
  
if (!is.element(k[2], categorical_variables))  x2 <- seq(min(X[,k[2]]),max(X[,k[2]]), length.out=50)
if (is.element(k[2], categorical_variables))    x2 <- unique(X[,k[2]])  

x <- expand.grid(x1,x2)
y <- sapply(1:nrow(x), function(i)
  mean(sapply(1:length(pur_res), function(s){
    
my.var.tree.index <- (1:length(pur_res[[s]]$my.functions))[ sapply(1:length(pur_res[[s]]$my.functions), function(i) setequal(as.numeric(pur_res[[s]]$my.functions[[i]]),k))]  ### all trees equal function component

if(length(my.var.tree.index )==0) return(NA)

if (!is.element(k[1], categorical_variables)){
bounds <-  unlist((pur_res[[s]]$lim_list)[[k[1]]])
bounds[length(bounds)]<- bounds[length(bounds)]+0.000001
pos1 <- which(x[i,1]<bounds)[1]-1
} else {
  leaves <-  ((pur_res[[s]]$lim_list)[[k[1]]])
pos1 <- sapply(1:ncol(leaves), function(j) is.element(x[i,1],  leaves[,j]))}


if (!is.element(k[2], categorical_variables)){
  bounds <-  unlist((pur_res[[s]]$lim_list)[[k[2]]])
  bounds[length(bounds)]<- bounds[length(bounds)]+0.000001
  pos2 <- which(x[i,2]<bounds)[1]-1
} else {
  leaves <-  ((pur_res[[s]]$lim_list)[[k[2]]])
  pos2 <-  sapply(1:ncol(leaves), function(j) is.element(x[i,2],  leaves[,j]))
}

if(sum(pos1)==0|sum(pos2)==0) return(NA) else (pur_res[[s]]$values)[[my.var.tree.index]][pos1,pos2]


  }
  ),na.rm=TRUE
  )
)


if (!all(is.element(k, categorical_variables))){
yy= findInterval(y, sort(y))
plot(x[,1], x[,2], col=pal(length(yy))[yy] , main=colnames(X)[k])
legend("topright", col=pal(2), pch=19,
       legend=c(round(range(y), 5)))
}else{
#order(x[1:11,1])

xy <- cbind(x,y)

#xy<- xy[order(xy[,1]),]


xyT<-xy
xyT[,1]<-findInterval(xy[,1], unique(sort(xy[,1])))
xyT[,2]<-findInterval(xy[,2], unique(sort(xy[,2])))

xy_k <- matrix(ncol=length(unique(xy[,2])), nrow=length(unique(xy[,1])))
for(j in 1:length(unique(xy[,2]))){
  for(i in 1:length(unique(xy[,1]))){
    xy_k[i,j] <- xy[apply(c(i,j)==t(xyT[,1:2]),2,all),3]
}}

rownames(xy_k)<- unique(sort(xy[,1]))
colnames(xy_k)<- unique(sort(xy[,2]))

barplot(xy_k, beside = TRUE ,main=colnames(X)[k], xlab=colnames(X)[k[2]])
legend("topright", fill=gray.colors(length(unique(sort(xy[,1])))),
       legend=unique(sort(xy[,1])), 5, title=colnames(X)[k[1]])

}
}

}
