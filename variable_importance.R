
var_importance <- function(pur_res,X, plot=TRUE, left_margin=7){

my_components <- unique(unlist(sapply(1:length(pur_res), function(s)
  pur_res[[s]]$my.functions[sapply(1:length(pur_res[[s]]$my.functions), function(z) sum(abs(pur_res[[s]]$values[[z]]))>0)]
),recursive = FALSE)) 



my_components <- my_components[order(sapply(my_components, function(x){ if (length(x)>1) {length(x)+10} else x   }      ))]

var_names <- sapply(1:length(my_components), function(k){
  paste(
  paste(colnames(X)[as.numeric(my_components[[k]])], collapse = ","),"(",
  paste(as.numeric(my_components[[k]]), collapse = ","),")",sep = "")
  }
  )

#var <- sapply(1:length(my_components), function(k){paste(as.numeric(my_components[[k]]), collapse = ",")})

  
importance <- sapply(1:length(my_components), function(k){
 mean( sapply(1:length(pur_res), function(s){
  
    my.var.tree.index <- (1:length(pur_res[[s]]$my.functions))[ sapply(1:length(pur_res[[s]]$my.functions), function(i) setequal(as.numeric(pur_res[[s]]$my.functions[[i]]),my_components[[k]]))]  ### all trees equal function component
    
    if(length(my.var.tree.index )==0) return(0) 
    if (sum(my_components[[k]])==0) return(abs(pur_res[[s]]$values[[my.var.tree.index]]))
    return(sum(pur_res[[s]]$individuals[[my.var.tree.index]]*abs(pur_res[[s]]$values[[my.var.tree.index]]))/sum(pur_res[[1]]$individuals[[1]]))
    
}))
})

res <- cbind(var_names,importance )
res2 <- res[order(as.numeric(res[,2])),]
par(mar = c(3, left_margin, 3, 3)) # Set the margin on all sides to 6
if (plot==TRUE) barplot(height=as.numeric(res2[,2]),main="Variable importance", horiz=TRUE,names.arg=res2[,1],las=1)

res[order(-as.numeric(res[,2])),]

}

