# Implementation of the random planted forest algorithm
# Input: (Y,X) = Data, max_interaction = maximum order of interaction, ntrees = Number of families of trees used, splits = number of splits,
#        split_try = number of considered splitpoints for each interval in each iteration step, t_try = percentage of viable trees considered in each iteration step,
#        variables = list of trees in each family in the beginning (default: one tree for each coordinate),
#        min_leaf_size = minimum number of nodes in each leaf, alternative = alternative updating
# Output: list of families of trees: [i][1] Final leaves of the trees, [i][2] = estimated values corresponding to the leaves, [i][3] = coordinates of the trees
#                                     (for the i-th family, i=1,...,ntrees)


oldrpf <- function(Y, X, max_interaction = 2, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4, variables = NULL, min_leaf_size = 1, alternative = F, loss = "L2", epsilon = 0.1, deterministic = FALSE, categorical_variables = NULL, delta = 0, cores = 1) {
  force(t_try)

  p <- ncol(X)
  n <- length(Y)
  eps <- 0.001 * apply(X, 2, function(s) min(diff(sort(unique(s)))))

  a <- apply(X, 2, min) ## lower bounds
  b <- apply(X, 2, max) + eps ### upper bounds

  if (loss == "exponential") {
    Y[Y == 0] <- -1
  }

  min_leaf_size <- rep(min_leaf_size, p)

  if (!is.null(categorical_variables)) {
    max_categorical <- max(sapply(1:length(categorical_variables), function(j) (length(unique(as.matrix(X[, categorical_variables[j]]))))))
  } else {
    max_categorical <- 2
  }


  tree_fam <- function(run) {
    #library(Rcpp)
    Rcpp::sourceCpp(here::here("oldrpf_compare", "C-Code.cpp"))

    subsample <- 1:n

    if (!deterministic) {
      subsample <- sample(n, n, replace = TRUE)


      X <- X[subsample, ]
      Y <- Y[subsample]
    }

    if (deterministic) t_try <- 1

    W <- rep(1, n)
    if (loss == "logit") W <- rep(0, n)

    # variables[[i]] is a vector of variables  in tree i
    variables <- list()
    variables[[1]] <- 0

    # Hier anderer start

    # intervals[[i]][[k]]= is list of matrices describing intervals for tree i, leaf k

    intervals <- list()

    intervals[[1]] <- list()
    intervals[[1]][[1]] <- matrix(nrow = max_categorical, ncol = p)
    for (j in 1:p) {
      intervals[[1]][[1]][1:2, j] <- c(a[j], b[j])
    }

    # values[[i]] is a vector of predictions for each leaf in tree i

    values <- list()
    values[[1]] <- 0

    # individuals[[i]][[k]] is a vector of individuals in leaf k of tree i

    individuals <- list()
    individuals[[1]] <- list(1:nrow(X))

    # Possible_Splits corresponds to the set of viable splits.
    # Possible_Splits[[i_1]][[1]] is the coordinate used for splitting, Possible_Splits[[i_1]][[2]] is the tree which results from the splitting

    Possible_Splits <- list()
    for (i in 1:p) {
      Possible_Splits[[length(Possible_Splits) + 1]] <- list(i, i)
    }

    while (splits > 0) {
      splits <- splits - 1

      m_try <- ceiling(t_try * length(Possible_Splits))

      if (deterministic) {
        split_candidates <- Possible_Splits
      } else {
        split_candidates <- sample(Possible_Splits, m_try)
      }

      if (!is.null(categorical_variables)) {
        R <- Calc_Optimal_split2(Y, W, as.matrix(X), split_try, variables, individuals, min_leaf_size, split_candidates, loss, categorical_variables, max_categorical, delta)
      } else {
        R <- Calc_Optimal_split2(Y, W, as.matrix(X), split_try, variables, individuals, min_leaf_size, split_candidates, loss, p + 1, 0, delta, deterministic)
      }

      R_opt <- R[1]

      Ikx_opt <- R[2:length(R)]
      Ikx_opt[6:length(Ikx_opt)][Ikx_opt[6:length(Ikx_opt)] == Inf] <- NA

      # print(c("Min Split Sum", R[1]))
      # print(c("Min Split Point", R[5]))

      if (R_opt < Inf) {
        if (Ikx_opt[5] == 1) {
          I_2 <- individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]], Ikx_opt[3]] < Ikx_opt[4]]
          I_1 <- individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]], Ikx_opt[3]] >= Ikx_opt[4]]
        } else {
          I_2 <- individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][is.element(X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]], Ikx_opt[3]], Ikx_opt[6:length(Ikx_opt)])]
          I_1 <- individuals[[Ikx_opt[1]]][[Ikx_opt[2]]][!is.element(X[individuals[[Ikx_opt[1]]][[Ikx_opt[2]]], Ikx_opt[3]], Ikx_opt[6:length(Ikx_opt)])]
        }


        if (min(length(I_2), length(I_1)) > 0) {
          if (loss == "median") {
            y_2 <- median(Y[I_2])
            y_1 <- median(Y[I_1])

            Y[I_2] <- Y[I_2] - y_2
            Y[I_1] <- Y[I_1] - y_1
          } else if (loss == "L1" | loss == "L2") {
            y_2 <- mean(Y[I_2])
            y_1 <- mean(Y[I_1])

            Y[I_2] <- Y[I_2] - y_2
            Y[I_1] <- Y[I_1] - y_1
          } else if (loss == "exponential") {
            #          R21 = mean(((Y[I_1]+1)/2) *(W[I_1])/sum(W[I_1]))
            #          R31 = mean(((Y[I_2]+1)/2) *(W[I_2])/sum(W[I_2]))
            R21 <- sum(((Y[I_1] + 1) / 2) * W[I_1]) / sum(W[I_1])
            R31 <- sum(((Y[I_2] + 1) / 2) * W[I_2]) / sum(W[I_2])

            # print(c("M_s = ", R31))
            # print(c("M_b = ", R21))

            R21 <- min(1 - epsilon, max(epsilon, R21))
            R31 <- min(1 - epsilon, max(epsilon, R31))

            if (sum(W[I_1]) == 0) {
              y_1 <- 0
            } else {
              y_1 <- log(R21 / (1 - R21))
            }
            W[I_1] <- W[I_1] * exp(-0.5 * Y[I_1] * y_1)

            if (sum(W[I_2]) == 0) {
              y_2 <- 0
            } else {
              y_2 <- log(R31 / (1 - R31))
            }
            W[I_2] <- W[I_2] * exp(-0.5 * Y[I_2] * y_2)


            W[I_1][is.infinite(y_1)] <- 0
            W[I_2][is.infinite(y_2)] <- 0
          } else if (loss == "logit") {
            y_22 <- min(1 - epsilon, max(epsilon, mean(Y[I_2])))
            y_11 <- min(1 - epsilon, max(epsilon, mean(Y[I_1])))

            y_2 <- log(y_22 / (1 - y_22)) - mean(W[I_2])
            y_1 <- log(y_11 / (1 - y_11)) - mean(W[I_1])

            W[I_1] <- W[I_1] + y_1
            W[I_2] <- W[I_2] + y_2
          }

          # print(y_1)
          # print(y_2)

          if (Ikx_opt[3] %in% variables[[Ikx_opt[1]]]) {
            ### if split variable is already in tree to be split

            individuals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]]) + 1]] <- I_2 #### individuals for  new leaf
            individuals[[Ikx_opt[1]]][[Ikx_opt[2]]] <- I_1 #### new split-leaf = remaining individuals

            intervals[[Ikx_opt[1]]][[length(intervals[[Ikx_opt[1]]]) + 1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]] # add one leaf to intervals


            if (Ikx_opt[5] == 1) {
              intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][2, Ikx_opt[3]] <- Ikx_opt[4] ### new leaf has new interval at split variable: (split value= upper bopund, lower bound remains)
              intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][1, Ikx_opt[3]] <- Ikx_opt[4] ### split leaf has new interval at split variable: (split value= lower bound, upper bound remains)
            } else {
              ## old leaf
              intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][, Ikx_opt[3]][is.element(intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][, Ikx_opt[3]], Ikx_opt[6:length(Ikx_opt)])] <- NA
              ## new leaf
              intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][, Ikx_opt[3]] <- NA
              intervals[[Ikx_opt[1]]][[length(individuals[[Ikx_opt[1]]])]][1:length(Ikx_opt[6:length(Ikx_opt)]), Ikx_opt[3]] <- Ikx_opt[6:length(Ikx_opt)]
            }
            values[[Ikx_opt[1]]][length(individuals[[Ikx_opt[1]]])] <- values[[Ikx_opt[1]]][Ikx_opt[2]] + y_2
            values[[Ikx_opt[1]]][Ikx_opt[2]] <- values[[Ikx_opt[1]]][Ikx_opt[2]] + y_1

            if (alternative) {
              for (j in 1:length(intervals[[Ikx_opt[1]]])) {
                y <- mean(Y[individuals[[Ikx_opt[1]]][[j]]])

                Y[individuals[[Ikx_opt[1]]][[j]]] <- Y[individuals[[Ikx_opt[1]]][[j]]] - y

                values[[Ikx_opt[1]]][j] <- values[[Ikx_opt[1]]][j] + y
              }
            }
          } else {
            TreeExists <- 0

            for (i in 1:length(variables)) {
              if ((Ikx_opt[1] != 1 && length(variables[[i]]) == length(variables[[Ikx_opt[1]]]) + 1) || (Ikx_opt[1] == 1 && length(variables[[i]]) == 1)) {
                ### if number of variables in tree i is same as num of var in new split tree
                if ((Ikx_opt[1] != 1 && all(variables[[i]] == sort(c(variables[[Ikx_opt[1]]], Ikx_opt[3])))) || (Ikx_opt[1] == 1 && variables[[i]] == Ikx_opt[3])) {
                  ### if variables in tree i are same as  var in new split tree

                  TreeExists <- 1 ###  tree i is the same as new split tree

                  individuals[[i]][[length(individuals[[i]]) + 1]] <- I_2 ## add one leaf
                  individuals[[i]][[length(individuals[[i]]) + 1]] <- I_1 ## add second leaf

                  intervals[[i]][[length(individuals[[i]])]] <- intervals[[i]][[length(individuals[[i]]) - 1]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]] ## add two leaves

                  if (Ikx_opt[5] == 1) {
                    intervals[[i]][[length(individuals[[i]]) - 1]][2, Ikx_opt[3]] <- Ikx_opt[4]
                    intervals[[i]][[length(individuals[[i]])]][1, Ikx_opt[3]] <- Ikx_opt[4]
                  } else {
                    ## second leaf
                    intervals[[i]][[length(individuals[[i]])]][, Ikx_opt[3]][is.element(intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][, Ikx_opt[3]], Ikx_opt[6:length(Ikx_opt)])] <- NA
                    ## first leaf
                    intervals[[i]][[length(individuals[[i]]) - 1]][, Ikx_opt[3]] <- NA
                    intervals[[i]][[length(individuals[[i]]) - 1]][1:length(Ikx_opt[6:length(Ikx_opt)]), Ikx_opt[3]] <- Ikx_opt[6:length(Ikx_opt)]
                  }
                  values[[i]][length(individuals[[i]]) - 1] <- y_2
                  values[[i]][length(individuals[[i]])] <- y_1
                }
              }
            }

            if (TreeExists == 0) {
              if (Ikx_opt[1] == 1) {
                variables[[length(variables) + 1]] <- Ikx_opt[3]
              } else {
                variables[[length(variables) + 1]] <- sort(c(variables[[Ikx_opt[1]]], Ikx_opt[3]))
              }
              individuals[[length(individuals) + 1]] <- list()
              individuals[[length(individuals)]][[1]] <- I_2
              individuals[[length(individuals)]][[2]] <- I_1

              if (length(variables[[length(variables)]]) < max_interaction) {
                for (i_1 in ((1:p)[-variables[[length(variables)]]])) {
                  Possible_exists <- 0
                  for (i_2 in 1:length(Possible_Splits)) {
                    if (Possible_Splits[[i_2]][[1]] == i_1 & length(Possible_Splits[[i_2]][[2]]) == length(variables[[length(variables)]]) + 1) {
                      if (all(sort(unique(c(i_1, variables[[length(variables)]]))) == Possible_Splits[[i_2]][[2]])) {
                        Possible_exists <- 1
                      }
                    }
                  }
                  if (Possible_exists == 0) {
                    Possible_Splits[[length(Possible_Splits) + 1]] <- list(i_1, sort(c(i_1, variables[[length(variables)]])))
                  }
                }
              }

              for (i_1 in variables[[length(variables)]]) {
                Possible_exists <- 0
                for (i_2 in 1:length(Possible_Splits)) {
                  if (Possible_Splits[[i_2]][[1]] == i_1 & length(Possible_Splits[[i_2]][[2]]) == length(variables[[length(variables)]])) {
                    if (all(variables[[length(variables)]] == Possible_Splits[[i_2]][[2]])) {
                      Possible_exists <- 1
                    }
                  }
                }
                if (Possible_exists == 0) {
                  Possible_Splits[[length(Possible_Splits) + 1]] <- list(i_1, variables[[length(variables)]])
                }
              }

              intervals[[length(intervals) + 1]] <- list()

              intervals[[length(intervals)]][[1]] <- intervals[[length(intervals)]][[2]] <- intervals[[Ikx_opt[1]]][[Ikx_opt[2]]]

              if (Ikx_opt[5] == 1) {
                intervals[[length(intervals)]][[1]][2, Ikx_opt[3]] <- Ikx_opt[4]
                intervals[[length(intervals)]][[2]][1, Ikx_opt[3]] <- Ikx_opt[4]
              } else {
                ## second leaf
                intervals[[length(intervals)]][[2]][, Ikx_opt[3]][is.element(intervals[[Ikx_opt[1]]][[Ikx_opt[2]]][, Ikx_opt[3]], Ikx_opt[6:length(Ikx_opt)])] <- NA
                ## first leaf
                intervals[[length(intervals)]][[1]][, Ikx_opt[3]] <- NA
                intervals[[length(intervals)]][[1]][1:length(Ikx_opt[6:length(Ikx_opt)]), Ikx_opt[3]] <- Ikx_opt[6:length(Ikx_opt)]
              }
              values[[length(variables)]] <- y_2
              values[[length(variables)]][2] <- y_1
            }
          }
        }
      }
    }

    return(list(intervals = intervals, values = values, variables = variables, b = b - eps, subsample = subsample, categorical_variables = categorical_variables, individuals = individuals, a = a))

    # return(list(intervals=intervals[[2:length(intervals)]], values=values[[2:length(values)]], variables=variables[[2:length(variables)]], b=b-eps, subsample=subsample, categorical_variables=categorical_variables, individuals=individuals[[2:length(individuals)]], a=a))
  }

  if (cores == 1) {
    forest_res <- sapply(1:ntrees, tree_fam)
  } else {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = ls(), envir = environment())
    forest_res <- parSapply(cl, 1:ntrees, tree_fam)
    stopCluster(cl)
  }
  # clusterExport
  # Y_hat=rep(0,n)
  # for(s in 1:ntrees){ Y_hat <- Y_hat + forest_res[[s]]$Y_hat }
  # Y_hat <- Y_hat/ntrees

  return(forest_res)
}






#
# ##### iris data
# data("iris")
# iris$Species <- as.numeric(iris$Species)
#  iris <- iris[iris$Species<=2, ]
#  iris$Species [iris$Species==2] <- 0
#
# forest_res <-   rpf(iris$Species, as.matrix(iris[,1:4]), max_interaction=2, ntrees=50, splits=20, split_try=10, t_try=0.4, variables=NULL, min_leaf_size=rep(1,4), alternative=FALSE, loss="logit",epsilon=0.001, categorical_variables=NULL)
#
# f<-predict_rpf(as.matrix(iris[,1:4]),forest_res )
#
# library(pROC)
# auc(as.numeric(iris[,5]), f)
#
# y_hat<-1/(1+exp(-f))
#
# View(cbind(y_hat,iris[,5]))
#
