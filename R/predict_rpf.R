predict <- function(rpf, X, components){
  if(is.matrix(X)) return(rpf$predict_matrix(X, components))
  if(is.vector(X)) return(rpf$predict_vector(X, components))
}
