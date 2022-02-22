#' Random Planted Forest Predictions
#'
#' @param object A fit object as returned by [`rpf()`].
#' @param new_data Matrix or vector input for new observatoins.
#' @param components TODO
#' @param ... Unused.
#'
#' @return Vector of predictions (FIXME: Type, dimension, response vs. prob)
#' @export
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
predict.Rcpp_RandomPlantedForest <- function(object, new_data, components, ...){
  if(is.matrix(new_data)) return(object$predict_matrix(new_data, components))
  if(is.vector(new_data)) return(object$predict_vector(new_data, components))
}
