
SBF.reg.LL<-function(formula,data,bandwidth,weight='sw',x.grid=NULL,n.grid.additional=0, x.min=NULL, x.max=NULL, integral.approx='right',it=100,kern=function(u){return(0.75*(1-u^2                         )*(abs(u)<1))},initial=NULL,kcorr=kcorr,LC,wrong,classic.backfit)
{  
  
p = ncol(data) - 1  
bandwidth=rep(bandwidth,p)
"div" <- function(x,y) ifelse(y==0&x==0,0,base:::"/"(x,y))

Terms <- terms(x=formula,data=data)
mm <- na.omit(get_all_vars(formula(Terms),data=data))
if (NROW(mm) == 0) stop("No (non-missing) observations")

response <- model.response(model.frame(update(formula,".~1"),data=mm))
X       <- prodlim::model.design(Terms,
                                 data=mm,
                                 maxOrder=1,
                                 dropIntercept=TRUE)[[1]]


Y <- response
#status <- as.vector(response[,"status"])


# }}}




taylor.f<-function(f,f.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
{
  
  taylor.f.i <- array(dim=c(n,d))

  for (k in 1:d){
    
    taylor.f.i[,k] <-    as.numeric( t(as.matrix(dx[[k]])) %*% ((f[[k]]+f.1[[k]]*dX.b[[k]])* K.X.b[[k]])  )
    
  }
  

  
  return(list( taylor.f.i= taylor.f.i))
}

get.f.new<-function(O1,O3,D,f,f.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
{

  if(classic.backfit==FALSE){
  
    taylor.f <- taylor.f(f,f.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if(dim(taylor.f$taylor.f.i)[2]==2) {taylor.f.minusk <-  taylor.f$taylor.f.i[,-k]} else
    taylor.f.minusk <-  rowSums(taylor.f$taylor.f.i[,-k]) 
    
  }
  
  
  if(classic.backfit==TRUE){
    temp <- array(dim=c(n,d)) 
    for (j in 1:d){
      A <- (dX.b[[j]])+1
      A <- A==1
      A  <- A*1
     temp[,j]<- f[[j]] %*% A
    }
    
    taylor.f.minusk<-array(dim=c(n,n.grid[1])) 
    taylor.f.minusk <-  rowSums(temp[,-k]) 
    
  }
    
    
  
      O2 <-  K.X.b[[k]]  %*%  taylor.f.minusk    
    
    
    

    O <- O1[[k]] - f.1[[k]]*O3[[k]]-O2


  
  
  return(as.numeric(div(O,D[[k]])))
} 

get.f.1.new<-function(O1.1,O3.1,D.1,f,f.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)
{

    taylor.f <- taylor.f(f,f.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    
    if(dim(taylor.f$taylor.f.i)[2]==2) {taylor.f.minusk <-  taylor.f$taylor.f.i[,-k]} else
      taylor.f.minusk <-  rowSums(taylor.f$taylor.f.i[,-k]) 
    
    
    

    
    O2.1 <- ( dX.b[[k]]*K.X.b[[k]] ) %*%  (taylor.f.minusk)   
    


    
    O.1 <- O1.1[[k]]-f[[k]]*O3.1[[k]]-O2.1
    
  
  
    return(as.numeric(div(O.1,D.1[[k]])))
} 

get.D<-function(k,dx,K.X.b,K.b,Y){
  
 
    D <-   rowSums(K.X.b[[k]]) 
  
  return(D)
}
  
get.O1<-function(k,Y,K.X.b){
  
  O1 <- as.numeric(  K.X.b[[k]] %*% Y )
  
  return(O1)
}


get.O3<-function(k,dx,K.X.b,K.b,Y){
  
  O3 <-   rowSums( dX.b[[k]]*K.X.b[[k]] )  
  
  return(O3) 
}

get.D.1<-function(k,dx,K.X.b,K.b,dX.b,dX0.b,Y){
  
 
    D <-   rowSums(dX.b[[k]]^2*K.X.b[[k]])  
    
  return(D)
}
    

get.O1.1<-function(k,Y,K.X.b,dX.b){
  O1<- as.numeric((dX.b[[k]]*K.X.b[[k]]) %*% Y)
  return(O1)
}
  
  
get.O3.1<-function(k,dx,K.X.b,K.b,dX.b,dX0.b,Y){
  
O3 <- rowSums( dX.b[[k]]*K.X.b[[k]] )
  
  return(O3)
}
    
    
d <- ncol(X) 
n <- nrow(X)

if(is.null(x.grid)) x.grid<-lapply(1:d,function(k) X[order(X[,k]),k])
if(is.null(x.min)) x.min<-sapply(x.grid,head,1)
if(is.null(x.max)) x.max<-sapply(x.grid,tail,1)



x.grid.additional <- lapply(1:d, function(k) seq(x.min[k],x.max[k], length=n.grid.additional))
x.grid <- lapply(1:d, function(k) sort(c(x.grid[[k]], x.grid.additional[[k]])))
n.grid <- sapply(x.grid, length)

if (integral.approx=='midd'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  (ddx[1]/2)+(x.grid[[k]][1]-x.min[k])  ,(ddx[-(n.grid[k]-1)]+ddx[-1])/2,  
      (ddx[n.grid[k]-1]/2)+(x.max[k]-x.grid[[k]][n.grid[k]])  
  )
  }
  )
}

if (integral.approx=='left'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  ddx,  x.max[k]-x.grid[[k]][n.grid[k]])  
  }
  )
}

if (integral.approx=='right'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  (x.grid[[k]][1]-x.min[k])  ,ddx)
  }
  )
}



dX.b<-K.X.b<-k.X.b<-list()
for( k in 1:d){
  K.X.b[[k]]<-array(0,dim=c(n,n.grid[k]))
  k.X.b[[k]]<-numeric(n)
  
  x.grid.array<-matrix(-X[,k],nrow =n.grid[k], ncol =n ,byrow=TRUE)   # 
  u<- x.grid.array+x.grid[[k]]      ####  u[,i]=x_k - X_{ik}
  K.X.b[[k]] <- apply(u/bandwidth[k],1:2,kern)/(bandwidth[k])
  k.X.b[[k]]<-colSums(dx[[k]]*K.X.b[[k]])   
  if (kcorr==FALSE) k.X.b[[k]] <- rep(1,n)
  dX.b[[k]]<- u               #### dX.b[[k]] [,i] = x_k-X_{ik} 
 # K.X.bC[[k]] <- K.X.b[[k]]/ k.X.b[[k]]
  K.X.b[[k]]<- K.X.b[[k]] %*%  diag(1/k.X.b[[k]])  ### row-wise division ---> col integration=1
  #(dx[[k]]*K.X.b[[k]]) =1
}


D<-O1<-O3<-D1.1<-O1.1<-O3.1<-list()

D <- lapply(1:d, get.D, dx,K.X.b,K.b,Y)
O1 <- lapply(1:d, get.O1, Y,K.X.b)
O3 <- lapply(1:d, get.O3, dx,K.X.b,K.b,Y)

D.1 <- lapply(1:d, get.D.1, dx,K.X.b,K.b,dX.b,dX0.b,Y)
O1.1 <- lapply(1:d, get.O1.1,Y,K.X.b,dX.b)
O3.1 <- lapply(1:d, get.O3.1,dx,K.X.b,K.b,dX.b,dX0.b,Y)

f_backfit<-list()
f.1_backfit<-list()

if (is.null(initial)){
  for(k in 1:d){
    f_backfit[[k]]<-rep(0, n.grid[k])
  }
} else  f_backfit<-initial

for(k in 1:d){
  f.1_backfit[[k]]<-rep(0, n.grid[k])
  
}


count<-rep(1,d)
for (l in 2:it)
{  
  
  
  
  f_backfit_old<-f_backfit
  f.1_backfit_old<-f.1_backfit
  for(k in 1:d)
  {
  
  
    f_backfit[[k]]<- get.f.new(O1,O3,D,f_backfit,f.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
    
  
     #f_backfit[[k]][is.nan(f_backfit[[k]])]<-0
    #f_backfit[[k]][is.na(f_backfit[[k]])]<-0
     # f_backfit[[k]][f_backfit[[k]]==Inf]<-0
     # f_backfit[[k]][f_backfit[[k]]==-Inf]<-0
  # f_backfit[[k]][f_backfit[[k]]<0]<-0.0001
    #  }
    #  for(k in 1:d)
    #   {
    #    if (k!=1) {
    
    #     }
    if (LC==TRUE) f.1_backfit[[k]]<-rep(0,n.grid[k]) else{
      f.1_backfit[[k]]<- get.f.1.new(O1.1,O3.1,D.1,f_backfit,f.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)}
  #  f.1_backfit[[k]][is.nan(f.1_backfit[[k]])]<-0
  #  f.1_backfit[[k]][is.na(f.1_backfit[[k]])]<-0
    # f.1_backfit[[k]][f.1_backfit[[k]]>30]<-30
    #  f.1_backfit[[k]][f.1_backfit[[k]]<=(-30)]<--30
    #   print(!(prod(f.1_backfit_old[[k]] == 0) ))
    
    # 
    # if ((prod(f.1_backfit_old[[k]] == 0) ))  {### if old vaues all zero
    #   if (max(abs(f.1_backfit_old[[k]]- f.1_backfit[[k]])) >(100))    {f.1_backfit[[k]] <-rep(0, n.grid[k])
    #   #  f.1_backfit[[k]]<-rep(0, n.grid[k])
    #   count[k]<-1
    #   }}
    # 
    # if (!(prod(f.1_backfit_old[[k]] == 0) )){ ### if old vaues NOT all zero
    #   if (max(abs(f.1_backfit_old[[k]]- f.1_backfit[[k]])) >(100/log(count[k])))    f.1_backfit[[k]] <-rep(0, n.grid[k])
    #   #  f.1_backfit[[k]]<-rep(0, n.grid[k])
    # }
    
    count[k]<-count[k]+1
    # plot(x.grid[[k]],f.1_backfit[[k]],lty=3,col=1,lwd=2) 
    
  }
  #}
  
  
  
  #plot(x.grid[[2]],phi[[1]](x.grid[[2]]),lty=3,col=1,lwd=2) 
  #  if (l==2) plot(x.grid[[2]],log(f_backfit[[2]]),ylim=c(-2,2) ) else lines(x.grid[[2]],log(f_backfit[[2]] ) )
  #plot(x.grid[[1]], log(f_backfit[[1]]),col='blue',lwd=2) 
  
  
    
    
  if (max(max(abs(unlist(f_backfit_old)-unlist(f_backfit)),na.rm=TRUE),max(abs(unlist(f.1_backfit_old)-unlist(f.1_backfit)),na.rm=TRUE))<= 0.001) break
  for(k in 1:d){
    print(c(l,max(abs(unlist((f_backfit_old[[k]]))-unlist((f_backfit[[k]]))),na.rm=TRUE),max(abs(unlist((f.1_backfit_old[[k]]))-unlist((f.1_backfit[[k]]
    ))),na.rm=TRUE)))
  }
}


for(k in 1:d){
  #f_backfit[[1]]<-f_backfit[[1]]+ (t(f_backfit[[k]])%*%dx[[k]])
 # print(f_backfit[[k]]%*%dx[[k]])
  f_backfit[[k]]<-f_backfit[[k]]- ((dx[[k]]%*%f_backfit[[k]])/sum(dx[[k]]))
  #print(f_backfit[[k]]%*%dx[[k]])
}

O1[[k]]<-((O1[[k]]/D[[k]])- ((dx[[k]]%*%(O1[[k]]/D[[k]]))/sum(dx[[k]])))*D[[k]]
return(list(f_backfit=f_backfit,f.1_backfit=f.1_backfit,l=l,x.grid=x.grid,dx=dx, O1=O1, D=D))
}
