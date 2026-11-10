// Internal utility helpers extracted from rpf.cpp to declutter large files.
// Kept minimal and header-only where templating is required.

#ifndef INTERNAL_UTILS_HPP
#define INTERNAL_UTILS_HPP

#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <limits>
#include <random>
#include <utility>

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

// The fast draw of `leaves` (exact_draw = false): the same age-weighted draw
// without replacement as the exact draw, by a Fenwick tree over the pool's
// weights, kept up to date between the splits of a tree family. Drawing k
// entries costs O(k log pool) instead of a random number per pool entry.
// Entries are drawn one at a time, each with probability proportional to its
// weight among those not yet drawn: the same distribution as the k smallest
// exponential keys -log(u) / w of the exact draw, but other random numbers.
// Methods in splits_leaves.cpp.
struct LeavesSampler {
  const void *owner = nullptr; // the pool the tree mirrors
  size_t synced = 0;           // pool entries [0, synced) are in the tree
  std::vector<double> w;       // the entries' weights
  std::vector<double> bit;     // their Fenwick tree (1-based nodes in bit[0..])
  size_t updates = 0;          // set() calls since the last rebuild()
  std::vector<double> taken;   // draw(): the drawn entries' weights, restored after it
  // min_leaf_rows: entries below `synced` whose leaf changed since the last
  // draw (their rows_ok reset), so their weights are set again before the next
  std::vector<size_t> dirty;
  // Forget the weights and mirror the pool `o` from its first entry.
  void reset(const void *o);
  // Recompute the tree from the weights (bounds the running sums' rounding).
  void rebuild();
  // Change entry i's weight.
  void set(size_t i, double nw);
  // Append an entry.
  void push(double nw);
  // The sum of all weights.
  double total() const;
  // Append to `out` k distinct entries with positive weight, drawn one at a
  // time in proportion to their weights (fewer if fewer have one).
  void draw(size_t k, std::vector<size_t> &out);
};

// Scratch buffers of the `leaves` split search (splits_leaves.cpp), reused
// across the candidates and splits of one tree family instead of allocated per
// candidate. Owned by create_tree_family's stack; TLS holds only a pointer to
// it, for the same reason as above. The fast draw's sampler lives here too, so
// it lasts exactly one tree family.
struct LeavesScratch {
  std::vector<double> unique, prefix, total;
  std::vector<size_t> first_pos;
  std::vector<std::pair<double, size_t>> keys;
  LeavesSampler sampler;
};

namespace rpf_utils {
// For tests of the two draws of `leaves`: how often each entry of a pool with
// the given ages is among the n_candidates drawn, over `reps` draws (exact or
// fast), with the random numbers of a generator seeded with `seed`. The fast
// draw's tree is kept across the draws, as in a fit. rows_ok: empty, or 0 / 1
// per entry, as min_leaf_rows sets it (0: never drawn).
std::vector<double> leaves_draw_counts(
    const std::vector<double> &ages,
    double decay,
    size_t n_candidates,
    size_t reps,
    bool exact,
    const std::vector<int> &rows_ok,
    uint64_t seed
);
} // namespace rpf_utils
extern thread_local LeavesScratch *tls_leaves_scratch_ptr;

#endif // INTERNAL_UTILS_HPP
