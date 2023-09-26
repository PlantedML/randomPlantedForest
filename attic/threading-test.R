library(randomPlantedForest)
library(mlbench)

xdat <- as.data.frame(mlbench.peak(n = 100, d = 10))

bm <- bench::mark(
  threads1 = rpf(y ~ ., data = xdat, nthreads = 1, ntrees = 100),
  threads2 = rpf(y ~ ., data = xdat, nthreads = 2, ntrees = 100),
  threads3 = rpf(y ~ ., data = xdat, nthreads = 3, ntrees = 100),
  threads4 = rpf(y ~ ., data = xdat, nthreads = 4, ntrees = 100),
  threads5 = rpf(y ~ ., data = xdat, nthreads = 5, ntrees = 100),
  min_iterations = 5,
  check = FALSE
)
bm

rpfit <- rpf(bikers ~ ., data = ISLR2::Bikeshare)


rpf(mpg ~ ., data = mtcars, nthreads = 1, ntrees = 10)
rpf(mpg ~ ., data = mtcars, nthreads = 2, ntrees = 3)



threads1 = rpf(y ~ ., data = xdat, nthreads = 1, ntrees = 100, deterministic = TRUE)
threads2 = rpf(y ~ ., data = xdat, nthreads = 2, ntrees = 100, deterministic = TRUE)
threads3 = rpf(y ~ ., data = xdat, nthreads = 3, ntrees = 100, deterministic = TRUE)
threads4 = rpf(y ~ ., data = xdat, nthreads = 4, ntrees = 100, deterministic = TRUE)
threads5 = rpf(y ~ ., data = xdat, nthreads = 5, ntrees = 100, deterministic = TRUE)

pred_threads1 = predict(threads1, xdat)
pred_threads2 = predict(threads2, xdat)
pred_threads3 = predict(threads3, xdat)
pred_threads4 = predict(threads4, xdat)
pred_threads5 = predict(threads5, xdat)

all.equal(pred_threads1, pred_threads2)
all.equal(pred_threads1, pred_threads3)
all.equal(pred_threads1, pred_threads4)
all.equal(pred_threads1, pred_threads5)

threads1$fit$get_parameters()
threads5$fit$get_parameters()


# iterate over families of trees and modify
if(nthreads > 1){
  # int nthreads = std::thread::hardware_concurrency() - 1;
  for(int n = 0; n < n_trees; n += nthreads){
    # Rcout << "n = " << n << std::endl;
    if(n >= (n_trees - nthreads)) {
      # This can lead to 0 threads, which is not ideal as the loop never finishes then
      nthreads = n_trees % nthreads;
      # Rcout << "Division gets " << nthreads << " threads" << std::endl;
    }
    # Rcout << "Setting " << nthreads << " threads" << std::endl;
    # Shoddy fix to ensure we have at leats one thread
    if (nthreads == 0) {
      nthreads += 1;
    }

    std::vector<std::thread> threads(nthreads);

    for(int t=0; t<nthreads; ++t){
      # Rcout << (n + t) << "/" << n_trees << std::endl;
      threads[t] = std::thread(&RandomPlantedForest::create_tree_family, this, std::ref(initial_leaves), n + t);
    }

    for(auto& t: threads){
      if(t.joinable()) t.join();
    }
  }
}
