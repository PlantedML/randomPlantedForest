#include "helper.hpp"
#include "internal_utils.hpp"

using namespace Rcpp;

namespace utils {

// Helper function to generate random number using R's RNG
// this replaces the previous randWrapper and later use of std::random_shuffle,
// as the latter is removed in C++17 and I couldn't figure out an easy replacement.
int random_index(const int n) { return static_cast<int>(rpf_utils::rng_runif01() * n); }

//  ----------------- functions for converting R and Cpp types -----------------


/**
 * \brief Convert the std container set of type int into an IntegerVector
 * from rcpp.
 *
 * \param v the vector that is converted.
 */
Rcpp::IntegerVector from_std_set(std::set<int> v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

/**
 * \brief Convert the std container vector of type int into an IntegerVector
 * from rcpp.
 *
 * \param v the vector that is converted.
 */
Rcpp::IntegerVector from_std_vec(std::vector<int> v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

/**
 * \brief Convert the std container vector of type double into a NumericVector
 * from rcpp.
 *
 * \param v the vector that is converted.
 */
Rcpp::NumericVector from_std_vec(std::vector<double> v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

/**
 * \brief Convert the nested std container vector containing a vector itself
 * of type double into a NumericMatrix from rcpp.
 *
 * Predefines a NumericMatrix of respective size. Afterwards iterates over the
 * outer vector, then transforms the inner vector using 'from_std_vec'-function
 * and inserts row into NumericMatrix at the correct position.
 *
 * \param v the vector of vectors that is converted.
 */
Rcpp::NumericMatrix from_std_vec(std::vector<std::vector<double>> v) {
  if(v.empty()) return Rcpp::NumericMatrix();
  Rcpp::NumericMatrix m(v.size(), v[0].size());
  for(unsigned int row=0; row<v.size(); ++row){
    m(row, _) = Rcpp::NumericMatrix(1, v[row].size(), from_std_vec(v[row]).begin());
  }
  return m;
}

std::vector<int> to_std_vec(Rcpp::IntegerVector rv) {
  return std::vector<int>(rv.begin(), rv.end());
}

std::vector<double> to_std_vec(Rcpp::NumericVector rv) {
  return std::vector<double>(rv.begin(), rv.end());
}

std::vector<std::vector<double>> to_std_vec(Rcpp::NumericMatrix rv) {
  // Efficient, allocation-minimal copy:
  // - Avoids per-row temporary NumericVector creation
  // - Copies column-major R memory into row-major nested vectors
  const int rows = rv.nrow();
  const int cols = rv.ncol();
  std::vector<std::vector<double>> X((size_t)rows, std::vector<double>((size_t)cols));
  if (rows == 0 || cols == 0) return X;
  const double *data = rv.begin(); // column-major, column stride = rows
  for (int j = 0; j < cols; ++j) {
    const double *colptr = data + (size_t)j * (size_t)rows;
    for (int i = 0; i < rows; ++i) {
      X[(size_t)i][(size_t)j] = colptr[(size_t)i];
    }
  }
  return X;
}

std::set<int> to_std_set(Rcpp::NumericVector rv) {
  return std::set<int>(rv.begin(), rv.end());
}

std::set<int> to_std_set(Rcpp::IntegerVector rv) {
  return std::set<int>(rv.begin(), rv.end());
}


}
