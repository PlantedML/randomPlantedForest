# ebm_feature_combination <- function(count_features_in_combination = 1) {
#   count_features_in_combination <- as.double(count_features_in_combination)
#   ret <- structure(list(count_features_in_combination = count_features_in_combination), class = "ebm_feature_combination")
#   return(ret)
# }

ebm_regression <- function(
  X, 
  y, 
  num_outer_bags = 16, 
  validation_size = 0.15, 
  max_epochs = 2000, 
  num_early_stopping_run_length = 50, 
  learning_rate = 0.01, 
  max_tree_splits = 2, 
  min_instances_for_split = 2, 
  random_state = 42
) {
  col_names = colnames(X)
  
  bin_edges <- vector(mode = "list") #, ncol(X))
  # TODO: I know this binning is buggy.  Review
  for(col_name in col_names) bin_edges[[col_name]] <- unique(quantile(X[[col_name]], seq(0,1, 1.0 / 256)))
  for(col_name in col_names) bin_edges[[col_name]] <- bin_edges[[col_name]][2:(length(bin_edges[[col_name]])-1)]
  for(col_name in col_names) X[[col_name]] <- as.integer(findInterval(X[[col_name]], bin_edges[[col_name]]))
  features <- lapply(col_names, function(col_name) { interpret:::ebm_feature(n_bins = length(bin_edges[[col_name]]) + 1) })
  feature_combinations <- create_main_feature_combinations(features)
  
  set.seed(random_state)
  val_indexes = sample(1:length(y), ceiling(length(y) * validation_size))
  
  X_train = X[-val_indexes,]
  y_train = y[-val_indexes]
  X_val = X[val_indexes,]
  y_val = y[val_indexes] 
  
  X_train_vec <- vector(mode = "numeric") # , ncol(X_train) * nrow(X_train)
  for(col_name in col_names) X_train_vec[(length(X_train_vec) + 1):(length(X_train_vec) + length(X_train[[col_name]]))] <- X_train[[col_name]]
  
  X_val_vec <- vector(mode = "numeric") # , ncol(X_val) * nrow(X_val)
  for(col_name in col_names) X_val_vec[(length(X_val_vec) + 1):(length(X_val_vec) + length(X_val[[col_name]]))] <- X_val[[col_name]]
  
  result_list = interpret:::cyclic_gradient_boost(
    "regression",
    2,
    features,
    feature_combinations,
    X_train_vec,
    y_train,
    NULL,
    X_val_vec,
    y_val,
    NULL,
    0,
    random_state,
    learning_rate,
    max_tree_splits, 
    min_instances_for_split, 
    max_epochs,
    num_early_stopping_run_length
  )
  model <- vector(mode = "list")
  for(i in seq_along(col_names)) {
    model[[col_names[[i]]]] <- result_list$model_update[[i]]
  }
  
  return(list(bin_edges = bin_edges, model = model))
}

create_main_feature_combinations <- function(features) {
  feature_combinations <- lapply(seq_along(features), function(i) { ebm_feature_combination(i) })
  return(feature_combinations)
}

ebm_prediction <- function (model, X) {
  col_names = colnames(X)
  X_binned <- vector(mode = "list") #, ncol(X))
  for(col_name in col_names) X_binned[[col_name]] <- as.integer(findInterval(X[[col_name]], model$bin_edges[[col_name]]) + 1)
  
  scores <- vector(mode = "numeric", nrow(X))
  for(col_name in col_names) {
    bin_vals <- model$model[[col_name]]
    bin_indexes <- X_binned[[col_name]]
    update_scores <- bin_vals[bin_indexes]
    scores <- scores + update_scores
  }
  
  return(scores)
}