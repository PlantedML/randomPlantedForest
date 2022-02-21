#' Random Planted Forest
#'
# FIXME: Parameters need describing - possibly also doable in C++? See also https://gallery.rcpp.org/articles/documenting-rcpp-packages/
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
#'
#' @examples
#' sample_size <- 500
#' data <- generate_data(Model=1, p=4, n=sample_size)
#' test_size <- floor(length(data$Y_true) / 5)
#' x_test <- data$X[1:test_size, ] # extract samples
#' y_test <- data$Y_true[1:test_size]
#' x_train <- data$X[(test_size+1):sample_size, ] # extract samples
#' y_train <- data$Y_start[(test_size+1):sample_size]
#'
#'
#' # set parameters ------------------------
#' n_splits <- 15
#' max_inter <- 2
#' n_trees <- 50
#' split_try <- 10
#' t_try <- 0.5
#' deterministic_forest <- TRUE
#' parallel <- TRUE
#' purify_forest <- FALSE
#' loss <- 'logit'
#' delta <- 0.1
#' epsilon <- 0
#'
#' # train models ------------------------
#' rpf_fit <- rpf(
#'   y_train, x_train,  max_interaction=max_inter, t_try=t_try,
#'   ntrees=n_trees, splits=n_splits, split_try = split_try,
#'   deterministic=deterministic_forest, parallel=parallel
#' )
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
