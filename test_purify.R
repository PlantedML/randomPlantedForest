
# source code ------------------------
rpf_path <- '/home/maike/Dokumente/HiWi/Planted_Forest'
setwd(rpf_path)
source(paste(rpf_path, '/randomPlantedForest.R', sep=''))

# generate test data ------------------------

sample_size <- 500
data <- generate_data(Model=1, p=4, n=sample_size)
test_size <- floor(length(data$Y_true) / 5)
x_test <- data$X[1:test_size, ] # extract samples
y_test <- data$Y_true[1:test_size]
x_train <- data$X[(test_size+1):sample_size, ] # extract samples
y_train <- data$Y_start[(test_size+1):sample_size]


# set parameters ------------------------
n_splits <- 15
max_inter <- 2
n_trees <- 50
split_try <- 10
t_try <- 0.5
deterministic_forest <- TRUE
parallel <- TRUE
purify_forest <- FALSE
loss <- 'logit' 
delta <- 0.1
epsilon <- 0


# train models ------------------------

start_time <- Sys.time()
rpf_cpp <- new_rpf(y_train, x_train,  max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try = split_try, deterministic=deterministic_forest, parallel=parallel)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
rpf_R <- rpf(y_train, x_train, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try=split_try, deterministic=deterministic_forest)
end_time <- Sys.time()
end_time - start_time


# purify ------------------------
rpf_cpp$purify()
rpf_cpp$new_purify()
rpf_R_purified <- rpf_purify(rpf_R)


# predict ------------------------
predictions_cpp <- predict(rpf_cpp, x_test, c(0))
predictions_R <- predict_rpf(x_test, rpf_R, c(0))


# accuracy ------------------------
variation <- mean(y_test^2)
MSE_cpp <- rpf_cpp$MSE(predictions_cpp, y_test) 
MSE_R <- sum((y_test - predictions_R)^2) / length(y_test)
MSE_cpp
MSE_R
variation


# --------------------------

res <- rpf_R
X <- x_train
s <- 1
X<-as.matrix(X)
values <- individuals <- list()

lim_list <- lapply(1:ncol(X),  ### go through all variables of the component
                   function(x){
                     my.var.tree.index2 <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) is.element(x,as.numeric(res[,s][[3]][[i]])))]  ## all trees that have x as variable

                     bounds <- sort(unique( na.omit( as.numeric(unlist(sapply(1:length(my.var.tree.index2), ### go through all relevant trees
                                                                              
                                                                              function(i){ 
                                                                                as.numeric(sapply(1:length( res[,s][[1]][[my.var.tree.index2[i]]]),
                                                                                                  function(j){ ### go through all leaves and get interval ends of variable
                                                                                                    res[,s][[1]][[my.var.tree.index2[i]]][[j]][,x]
                                                                                                  }))
                                                                              }))))))
                     # bounds[length(bounds)] <- bounds[length(bounds)] +0.001 ### may need fix
                     return(bounds)
                   })

my.functions <- unique(res[,s][[3]]) ### all function components


for (l in 1:length( my.functions)){  ### go through all  function components
  
  my.var.tree.index <- (1:length(res[,s][[3]]))[ sapply(1:length(res[,s][[3]]), function(i) setequal(as.numeric(res[,s][[3]][[i]]),my.functions[[l]]))]  ### all trees equal function component
  # print(my.var.tree.index)
  
  if (l==1){values <- individuals<-list()}
  values[[l]] <- individuals[[l]] <- array(0, dim=sapply(1:length(my.functions[[l]]), function(m) {
      if (!any(my.functions[[l]][m]==res[,s][[6]])) {length(lim_list[[my.functions[[l]][m]]] )-1} else{
        ncol(as.matrix(lim_list[[my.functions[[l]][m]]]))
      }
    }
  ))
  
  ## x is index through all leaves
  x<-expand.grid(   lapply(1:length(my.functions[[l]]), function(j)
    {
      print(l)
      1:(length(lim_list[[my.functions[[l]][j]]])-1)
    }
    
  )) 
  
}


