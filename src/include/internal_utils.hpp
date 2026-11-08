// Internal utility helpers extracted from rpf.cpp to declutter large files.
// Kept minimal and header-only where templating is required.

#ifndef INTERNAL_UTILS_HPP
#define INTERNAL_UTILS_HPP

#include <vector>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <limits>
#include <random>

#include "trees.hpp"

namespace rpf_utils {

// RNG helpers
double rng_runif01();
double rng_runif(double a, double b);
int rng_randint(int left_inclusive, int right_exclusive);
// Swap the thread-local RNG pointer; returns previous pointer
std::mt19937_64 *swap_tls_rng(std::mt19937_64 *new_ptr);

// Leaf/order/prefix helpers
void ensure_order_and_sorted_vals_for_leaf(
    const std::vector<std::vector<double>> &X,
    Leaf &leaf,
    int k,
    std::vector<size_t> &order_out,
    std::vector<double> &sorted_vals_out
);

std::vector<double> compute_unique_sorted_values(const std::vector<double> &sorted_vals);

void build_prefix_and_total_given_order(
    const std::vector<std::vector<double>> &Y,
    const Leaf &leaf,
    const std::vector<size_t> &order,
    size_t value_size,
    std::vector<std::vector<double>> &prefix_out,
    std::vector<double> &total_out
);

void finalize_split_from_sums(Split &winner, const std::vector<std::vector<double>> &X, size_t value_size);

// Sampling helpers
std::vector<size_t> sample_weighted_indices_filtered(const std::vector<double> &weights, size_t n_candidates);

std::vector<int> compute_even_spread_indices(int left_inclusive, int right_exclusive, size_t max_draws);
std::vector<int> sample_unique_ints_uniform_R(int left_inclusive, int right_exclusive, size_t k);

// Fenwick helpers used by cur_trees_1 sampling cache
void fenwick_add(std::vector<double> &bit, size_t idx1, double delta);
size_t fenwick_find_by_prefix(const std::vector<double> &bit, double target);
// Append an element with value w (the new node also sums the elements of its range before it)
void fenwick_push_back(std::vector<double> &bit, double w);

// Sampling weight of a leaf in dimension k for res_trees: its number of available thresholds,
// max(0, unique values of X[, k] in the leaf - 2 * leaf_size + 1), the width of the range [leaf_size,
// unique - leaf_size] the samplers draw from. Caches the unique count in leaf.unique_count_cache. cur_trees_1
// computes the same weight from the family's ranks (rpf_ct1::leaf_weight in cur_trees_1.hpp), in both its lazy
// Fenwick rebuild and the incremental update after a split, so that the two agree.
double leaf_threshold_weight(const std::vector<std::vector<double>> &X, Leaf &leaf, int k, int leaf_size);

// Aging helper must be header (templated)
template <typename CandidateT>
inline void age_pool_by_sample(const std::vector<size_t> &sample_idxs, int best_idx, std::vector<CandidateT> &pool)
{
  for (size_t idx : sample_idxs) {
    if (static_cast<int>(idx) != best_idx)
      pool[idx].age += 1.0;
    else
      pool[idx].age = 0.0;
  }
}

} // namespace rpf_utils

// Thread-local pointer to the working-set bin cache used by histogram split
// mode (mode 4). The buffer itself is owned by create_tree_family's stack;
// only a trivially-destructible pointer lives in TLS. A thread_local object
// with a non-trivial destructor (e.g. std::vector) must NOT be used here:
// its emutls destructor crashes at std::thread exit on Windows (MinGW) when
// the package DLL is loaded via LoadLibrary.
extern thread_local std::vector<std::vector<int>> *tls_working_bin_id_ptr;

#endif // INTERNAL_UTILS_HPP
