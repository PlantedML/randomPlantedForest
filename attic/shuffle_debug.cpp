#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Previous cpp11 approach
inline int randWrapper(const int n) { return floor(R::runif(0, 1) * n); }

template <typename Iter>
void shuffle_vector_cpp11(Iter first, Iter last) {
  std::random_shuffle(first, last, randWrapper);
}

// [[Rcpp::export]]
NumericVector shuffle_cpp11(NumericVector x) {
  shuffle_vector_cpp11(x.begin(), x.end());

  return x;
}

int random_index(const int n) { return static_cast<int>(R::runif(0, 1) * n); }

// cpp17 compatible
template <typename Iter>
void shuffle_vector_cpp17(Iter first, Iter last) {
  int n = std::distance(first, last);
  while (n > 1) {
    int k = random_index(n--);
    std::swap(*(first + n), *(first + k));
  }
}

// [[Rcpp::export]]
NumericVector shuffle_cpp17(NumericVector x) {
  shuffle_vector_cpp17(x.begin(), x.end());

  return x;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
shuffle_cpp11(1:10)
shuffle_cpp17(1:10)
*/
