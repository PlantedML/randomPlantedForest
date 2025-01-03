#include "random_utils.hpp"

namespace utils
{
  thread_local std::mt19937 StdRNGBackend::generator;
  bool StdRNGBackend::seeded = false;
  RNGBackend *RandomGenerator::backend = &RandomGenerator::std_backend;
  StdRNGBackend RandomGenerator::std_backend;
  RcppRNGBackend RandomGenerator::rcpp_backend;
}
