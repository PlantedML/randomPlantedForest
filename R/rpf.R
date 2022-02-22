#' Random Planted Forest
#'
# FIXME: Parameters need describing
#' @param Y Vector of target features
#' @param X Feature matrix
#' @param max_interaction Maximum level of interaction determining maximum
#'   number of split dimensions for a tree
#' @param ntrees Number of trees generated per family
#' @param splits Number of performed splits for each tree family
#' @param split_try TODO
#' @param t_try TODO
#' @param deterministic Choose whether approach deterministic or random
#' @param parallel Perform algorithm in parallel or serialized
#' @param purify Whether the forest should be purified
#' @param cv Determines if cross validation is performed
#' @param loss "L2" or "logit" / TODO
#' @param delta TODO
#' @param epsilon TODO
#'
#' @return Object
#' @export
#' @importFrom methods new
#' @examples
#' \dontrun{
#' # TODO (see tests)
#' }
rpf <- function(Y, X, max_interaction=1, ntrees=50, splits=30, split_try=10, t_try=0.4,
                    deterministic=FALSE, parallel=FALSE, purify=FALSE, cv=FALSE,
                    loss='L2', delta=0, epsilon=0.1){
  if(!missing(loss) | !missing(delta) | !missing(epsilon)){
    return(new(ClassificationRPF, Y, X, loss, c(max_interaction, ntrees, splits, split_try, t_try,
                                                purify, deterministic, parallel, cv, delta, epsilon)))
  }
  return(new(RandomPlantedForest, Y, X, c(max_interaction, ntrees, splits, split_try, t_try,
                                          purify, deterministic, parallel, cv)))
}
