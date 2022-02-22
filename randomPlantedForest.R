sigma <- function(x){
  1 / (1 + exp(-x))
}

benchmark_rpf <- function(start=1, end=100){
  problem_size <- seq(start, end, floor((end-start)/10)) # todo: exponential increase
  timings_c_seq <- c()
  timings_c_par <- c()
  timings_r <- c()
  for(N in problem_size){
    # generate data for regression
    data <- generate_data(Model=1, p=4, n=N)

    # specify parameters
    n_splits <- 15
    max_inter <- 2
    n_trees <- 50
    split_try <- 10
    t_try <- 0.5
    deterministic_forest <- TRUE
    purify_forest <- FALSE
    loss <- 'L2'
    delta <- 0.1
    epsilon <- 0

    # train models and measure time
    start_time <- Sys.time()
    rpf_cpp <- new_rpf(data$Y_start, data$X,  max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try = split_try, deterministic=deterministic_forest, parallel=F)
    end_time <- Sys.time()
    timings_c_seq[length(timings_c_seq)+1] <- end_time - start_time

    start_time <- Sys.time()
    rpf_cpp <- new_rpf(data$Y_start, data$X,  max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try = split_try, deterministic=deterministic_forest, parallel=T)
    end_time <- Sys.time()
    timings_c_par[length(timings_c_par)+1] <- end_time - start_time

    start_time <- Sys.time()
    rpf_R <- rpf(data$Y_start, data$X, max_interaction=max_inter, t_try=t_try, ntrees=n_trees, splits=n_splits, split_try=split_try, deterministic=deterministic_forest)
    end_time <- Sys.time()
    timings_r[length(timings_r)+1] <- end_time - start_time
  }

  plot(problem_size, timings_r, type="b", col="red", xlab='N', ylab='time')
  lines(problem_size, timings_c_seq, type="b", col="blue")
  lines(problem_size, timings_c_par, type="b", col="green")
  legend(x="topleft", legend=c("R", "Cpp_sequential", "Cpp_parallel"), col=c("red", "blue", "green"), lty=1:2, cex=0.8)

  return(list(problem_size, timings_c, timings_r_seq, timings_c_par))
}
