#include "helper.hpp"
#include "internal_utils.hpp"

namespace utils {

// Uniform index in [0, n) from the calling thread's per-tree RNG stream.
int random_index(const int n) { return static_cast<int>(rpf_utils::rng_runif01() * n); }

}
