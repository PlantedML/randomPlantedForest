data(mtcars)
source("rpf.R")
source("purify_2.R")
source("plot_purified.R")
source("predict_purified.R")
source("variable_importance.R")
source("predict_rpf.R")


categorical_variables <- c(7,8)



res <- rpf(as.numeric(mtcars[1:20,1]), as.matrix(mtcars[1:20,-1]), max_interaction=3, ntrees=50, splits=20, split_try=5, t_try=0.5, variables=NULL, min_leaf_size=1, alternative=F, loss="L2", epsilon=0.1, categorical_variables=categorical_variables, delta=0, cores=4)


y.hat <- predict_rpf(as.matrix(mtcars[21:30,-1]),res) ## get predictions without purifying



pur_res <- purify(res,as.matrix(mtcars[1:20,-1]),cores=4)


## a list of all non-zero trees. 
unique(unlist(sapply(1:length(pur_res), function(s)
   pur_res[[s]]$my.functions[sapply(1:length(pur_res[[s]]$my.functions), function(z) sum(abs(pur_res[[s]]$values[[z]]))>0)]
  ),recursive = FALSE)) 





y.hat2 <- pred_pur(as.matrix(mtcars[21:30,-1]),pur_res, interpret=FALSE)  #### get predictions after puryfing, note y.hat=y.hat2 


y.hat3 <- pred_pur(as.matrix(mtcars[21:30,-1]),pur_res, interpret=TRUE) #### get predictions after puryfing, winth local explanation



## get variable importance
X <- as.matrix(mtcars[1:20,-1])
var_importance(pur_res,X, plot=TRUE, left_margin=7)

  
k <- c(1,5)  

plot_pur(pur_res,X, categorical_variables, k)



#### compare mse on test set (last 10 rows) eith lm and xgboost

data<-data.frame(mtcars[1:20,-1])
lm.res <- lm(as.numeric(mtcars[1:20,1])~., data)
sum((as.numeric(mtcars[21:30,1])-predict(lm.res,data.frame(mtcars[21:30,-1])))^2)/30

sum((as.numeric(mtcars[21:30,1])-pred_pur(as.matrix(mtcars[21:30,-1]),pur_res, interpret=FALSE))^2)/30
sum((as.numeric(mtcars[21:30,1])-predict_rpf(as.matrix(mtcars[21:30,-1]),res))^2)/30

library(xgboost)
param <- list(max_depth =8, eta = 0.1, verbose = 2, nthread = 1)
dtrain <- xgb.DMatrix(as.matrix(mtcars[1:20,-1]), label = as.numeric(mtcars[1:20,1]))
bst <- xgb.train(param, dtrain, nrounds = 1000)

#bst<-xgboost(as.matrix(mtcars[1:20,-1]),label = as.numeric(mtcars[1:20,1]),max_depth = 5, eta = 1, verbose = 2, nthread = 1, nrounds = 10)
               
sum((as.numeric(mtcars[21:30,1])-predict(bst, as.matrix(mtcars[21:30,-1])))^2)/30



