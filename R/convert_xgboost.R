
# Function to get leaf bounds from xgboost
#FIXME: This seems to be the bottleneck, could do it in C++
get_leaf_bounds <- function(trees, tree,x) {
  max_node <- trees[Tree == tree, max(Node)]
  num_nodes <- max_node + 1
  lb <- matrix(-Inf, nrow = num_nodes, ncol = ncol(x))
  ub <- matrix(Inf, nrow = num_nodes, ncol = ncol(x))
  for (nn in 0:max_node) {
    if (trees[Tree == tree & Node == nn, !is.na(Yes)]) {
      left_child <- trees[Tree == tree & Node == nn, Yes]
      right_child <- trees[Tree == tree & Node == nn, No]
      splitvar <- trees[Tree == tree & Node == nn, Feature_num]

      # Children inherit bounds
      ub[left_child + 1, ] <- ub[nn + 1, ]
      ub[right_child + 1, ] <- ub[nn + 1, ]
      lb[left_child + 1, ] <- lb[nn + 1, ]
      lb[right_child + 1, ] <- lb[nn + 1, ]

      # Restrict by new split
      ub[left_child + 1, splitvar] <- trees[Tree == tree & Node == nn, Split]
      lb[right_child + 1, splitvar] <- trees[Tree == tree & Node == nn, Split]
    }
  }

  # Return bounds of leaves only
  leaves <- trees[Tree == tree & Feature == "Leaf", Node+1]
  list(lower = lb[leaves, ],
       upper = ub[leaves, ])
}

#' Convert xgboost to rpf object
#'
#' @param xg xgboost object
#' @param x data used to train the xgboost model
#' @param y target used to train the xgboost model
#'
#' @return rpf object
#' @importFrom xgboost xgb.model.dt.tree
#' @export
convert_xgboost_rpf <- function(xg, x, y) {
  trees <- xgboost::xgb.model.dt.tree(model = xg, use_int_id = TRUE)
  trees[, Feature_num := as.integer(factor(Feature, levels = c("Leaf", colnames(x)))) - 1L]
  num_trees <- trees[, max(Tree)+1]

  # create a dummy rpf
  rpfit <- rpf(x = x, y = y, max_interaction = 0, ntrees = num_trees, splits = 1,
               purify = FALSE)

  # Overwrite rpf trees
  for (t in seq_len(num_trees)) {
    # xgboost adds 0.5 to prediction
    rpfit$forest[[t]]$values[[1]][[1]] <- 0.5
    rpfit$forest[[t]]$variables[[2]] <- trees[Tree == t-1 & Feature_num > 0, sort(unique(Feature_num))]
    rpfit$forest[[t]]$values[[2]] <- as.list(as.numeric(num_trees)*trees[Tree == t-1 & Feature == "Leaf", Quality])

    rpfit$forest[[t]]$intervals[[2]] <- rep(rpfit$forest[[t]]$intervals[[1]], length(rpfit$forest[[t]]$values[[2]]))

    # Get leaf bounds
    leaf_bounds <- get_leaf_bounds(trees, t-1,x)
    leaves <- trees[Tree == t-1 & Feature == "Leaf", Node+1]
    for (i in seq_along(leaves)) {
      rpfit$forest[[t]]$intervals[[2]][[i]][1, ] <- pmax(rpfit$forest[[t]]$intervals[[2]][[i]][1, ],
                                                         leaf_bounds$lower[i, ])
      rpfit$forest[[t]]$intervals[[2]][[i]][2, ] <- pmin(rpfit$forest[[t]]$intervals[[2]][[i]][2, ],
                                                         leaf_bounds$upper[i, ])
    }
  }

  # Also overwrite C++ forest
  rpfit$fit$set_model(rpfit$forest)

  # Return manipulated rpf object
  rpfit
}
