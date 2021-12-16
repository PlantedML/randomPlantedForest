# plot the partial dependence plot of a random planted forest model (up to 2 dimensional) 
# Input: forest_res = random planted forest model, x = Vector of values creating the support of the plot, i = Coordinate to be plotted 
# Output: plot of a random planted forest model on the support x. If x=NULL the whole support is considered.

plot_rpf=function(forest_res, x=NULL, i=1, Gridsize = 30, type = "p",  xlim = NULL, ylim = NULL,
                  log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                  ann = par("ann"), axes = TRUE, frame.plot = axes,
                  panel.first = NULL, panel.last = NULL, asp = NA,
                  xgap.axis = NA, ygap.axis = NA, breaks = -15:15/7,
                  cex.axis = 1.5,cex.lab=1.5,cex.main=1.5){
  
  p=length(forest_res[,1]$intervals[[1]][[1]][1,])
  
  if(length(i)==2){
    
    if(max(i)>p){ return(paste("The coordinate i can only take on integers smaller then",p))  }
    
    if(!is.null(x)){ print("The X value is ignored when creating heatmaps!") }
    
    a=c(forest_res[,1]$intervals[[2]][[1]][1,1],forest_res[,1]$intervals[[1]][[1]][1,2:p])
    
    b=c(forest_res[,1]$intervals[[2]][[1]][2,1],forest_res[,1]$intervals[[1]][[1]][2,2:p])
    
    x=seq(a[i[1]],b[i[2]], length.out = Gridsize)
    y=seq(a[i[1]],b[i[2]], length.out = Gridsize)
    
    f_1=function(x,y){ return(predict_rpf(c(x,y), forest_res, i)) }
    f_2=function(x,y){ return(mapply(f_1,x,y)) }
    
    z=outer(x,y,f_2)
    
    image(x,
          y,
          z, 
          col=grey.colors(30),
          xlab=xlab,
          ylab=xlab,
          main=main,
          breaks = breaks,cex.axis = cex.axis,cex.lab=cex.lab,cex.main=cex.main)
  } else {
    
    if(length(i)!=1|i[1]>p){ return(paste("The coordinate i can only take on integers smaller then",p)) }
    
    if(is.null(x)){
      
      if(i==1){
        
        a=forest_res[,1]$intervals[[2]][[1]][1,1]
        
        b=forest_res[,1]$intervals[[2]][[1]][2,1]
      } else {
        
        a=forest_res[,1]$intervals[[1]][[1]][1,i]
        
        b=forest_res[,1]$intervals[[1]][[1]][2,i]
      }
      
      x=seq(a,b,length.out = 100)
    } else {
      
      x=sort(x)
    }
    
    dim(x)=c(length(x),1)
    
    plot(x,predict_rpf(x,forest_res,i=i), type = type,  xlim = xlim, ylim = ylim,
         log = log, main = main, sub = sub, xlab = xlab, ylab = ylab,
         ann = ann, axes = axes, frame.plot = frame.plot,
         panel.first = panel.first, panel.last = panel.last, asp = asp,
         xgap.axis = xgap.axis, ygap.axis = ygap.axis)
  }
}