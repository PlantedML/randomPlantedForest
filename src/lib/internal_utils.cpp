// Internal utilities for split sampling, RNG, caching, and prefix sums shared
// across split modes and both regression/classification flows.
//
// These helpers centralize frequently reused logic and are intentionally kept
// low-level and stateless, using thread-local state only for RNG where needed.
#include "internal_utils.hpp"
#include <random>
#include <cmath>

namespace {
  // Thread-local RNG pointer used in worker threads for reproducible randomness
  thread_local std::mt19937_64* tls_rng_ptr = nullptr;
}

namespace rpf_utils {
void fenwick_add(std::vector<double> &bit, size_t idx1, double delta)
{
  // bit is 1-based; idx1 in [1, bit.size()]
  size_t n = bit.size();
  while (idx1 <= n) { bit[idx1 - 1] += delta; idx1 += idx1 & (~idx1 + 1); }
}

size_t fenwick_find_by_prefix(const std::vector<double> &bit, double target)
{
  // Return smallest i such that sum(i) >= target; 1-based index
  size_t n = bit.size();
  size_t idx = 0; double sum = 0.0;
  // Largest power of two <= n
  size_t step = 1ULL << (63 - __builtin_clzll((unsigned long long)std::max<size_t>(1, n)));
  while (step) {
    size_t next = idx + step; if (next <= n) {
      double val = bit[next - 1];
      if (sum + val < target) { sum += val; idx = next; }
    }
    step >>= 1;
  }
  return std::min(n, idx + 1);
}

std::mt19937_64* swap_tls_rng(std::mt19937_64* new_ptr)
{
  std::mt19937_64* old = tls_rng_ptr;
  tls_rng_ptr = new_ptr;
  return old;
}

double rng_runif01()
{
  if (tls_rng_ptr) {
    return std::generate_canonical<double, 53>(*tls_rng_ptr);
  }
  static thread_local std::mt19937_64 fallback_rng(0x9E3779B97F4A7C15ULL);
  return std::generate_canonical<double, 53>(fallback_rng);
}

double rng_runif(double a, double b)
{
  double u = rng_runif01();
  return a + u * (b - a);
}

int rng_randint(int left_inclusive, int right_exclusive)
{
  if (right_exclusive <= left_inclusive) return left_inclusive;
  if (tls_rng_ptr) {
    std::uniform_int_distribution<int> dist(left_inclusive, right_exclusive - 1);
    return dist(*tls_rng_ptr);
  }
  static thread_local std::mt19937_64 fallback_rng(0xD1B54A32D192ED03ULL);
  std::uniform_int_distribution<int> dist(left_inclusive, right_exclusive - 1);
  return dist(fallback_rng);
}

void ensure_order_and_sorted_vals_for_leaf(
    const std::vector<std::vector<double>> &X,
    Leaf &leaf,
    int k,
    std::vector<size_t> &order_out,
    std::vector<double> &sorted_vals_out)
{
  const size_t m = leaf.individuals.size();
  if (leaf.order_cache.count(k) && leaf.order_cache[k].size() == m) {
    order_out = leaf.order_cache[k];
  } else {
    order_out.resize(m);
    std::iota(order_out.begin(), order_out.end(), 0);
    std::stable_sort(order_out.begin(), order_out.end(), [&](size_t a, size_t b){
      return X[leaf.individuals[a]][k] < X[leaf.individuals[b]][k];
    });
    leaf.order_cache[k] = order_out;
  }
  if (leaf.sorted_vals_cache.count(k) && leaf.sorted_vals_cache[k].size() == m) {
    sorted_vals_out = leaf.sorted_vals_cache[k];
  } else {
    sorted_vals_out.resize(m);
    for (size_t i = 0; i < m; ++i)
      sorted_vals_out[i] = X[leaf.individuals[order_out[i]]][k];
    leaf.sorted_vals_cache[k] = sorted_vals_out;
  }
}

std::vector<double> compute_unique_sorted_values(const std::vector<double> &sorted_vals)
{
  std::vector<double> unique;
  unique.reserve(sorted_vals.size());
  if (!sorted_vals.empty()) {
    unique.push_back(sorted_vals[0]);
    for (size_t i = 1; i < sorted_vals.size(); ++i)
      if (sorted_vals[i] != unique.back()) unique.push_back(sorted_vals[i]);
  }
  return unique;
}

void build_prefix_and_total_given_order(
    const std::vector<std::vector<double>> &Y,
    const Leaf &leaf,
    const std::vector<size_t> &order,
    size_t value_size,
    std::vector<std::vector<double>> &prefix_out,
    std::vector<double> &total_out)
{
  const size_t m = leaf.individuals.size();
  prefix_out.assign(value_size, std::vector<double>(m, 0.0));
  for (size_t p = 0; p < value_size; ++p) {
    double acc = 0.0;
    for (size_t i = 0; i < m; ++i) {
      acc += Y[leaf.individuals[order[i]]][p];
      prefix_out[p][i] = acc;
    }
  }
  total_out.assign(value_size, 0.0);
  for (size_t p = 0; p < value_size; ++p)
    total_out[p] = prefix_out[p][m - 1];
}

void finalize_split_from_sums(
    Split &winner,
    const std::vector<std::vector<double>> &X,
    size_t value_size)
{
  if (std::isinf(winner.min_sum) || winner.leaf_index == nullptr) return;
  const int kfin = winner.split_coordinate - 1;
  Leaf &leaf_fin = *winner.leaf_index;
  const double sp_fin = winner.split_point;
  winner.I_s.clear(); winner.I_b.clear();
  for (int ind : leaf_fin.individuals) {
    if (X[ind][kfin] < sp_fin) winner.I_s.push_back(ind); else winner.I_b.push_back(ind);
  }
  winner.M_s.assign(value_size, 0.0);
  winner.M_b.assign(value_size, 0.0);
  if (!winner.I_s.empty()) for (size_t p = 0; p < value_size; ++p)
    winner.M_s[p] = winner.sum_s[p] / static_cast<double>(winner.I_s.size());
  if (!winner.I_b.empty()) for (size_t p = 0; p < value_size; ++p)
    winner.M_b[p] = winner.sum_b[p] / static_cast<double>(winner.I_b.size());
}

std::vector<size_t> sample_weighted_indices_filtered(
    const std::vector<double> &weights,
    size_t n_candidates)
{
  std::vector<size_t> pos_idx; pos_idx.reserve(weights.size());
  std::vector<double> pos_w;   pos_w.reserve(weights.size());
  for (size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0.0) { pos_idx.push_back(i); pos_w.push_back(weights[i]); }
  const size_t P = pos_idx.size();
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (P == 0) {
    std::vector<size_t> all(weights.size()); std::iota(all.begin(), all.end(), 0);
    size_t k = std::min(n_candidates, all.size());
    for (size_t i = 0; i < k; ++i) {
      size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i));
      if (j >= all.size()) j = all.size() - 1; std::swap(all[i], all[j]);
    }
    for (size_t i = 0; i < k; ++i) sample_idxs.push_back(all[i]);
  } else {
    size_t k = std::min(n_candidates, P);
    std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
    for (size_t i = 0; i < P; ++i) {
      double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min();
      double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]);
    }
    if (k < keys.size()) {
      std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
      keys.resize(k);
    }
    for (auto &kv : keys) sample_idxs.push_back(kv.second);
  }
  return sample_idxs;
}

std::vector<int> compute_even_spread_indices(int left_inclusive, int right_exclusive, size_t max_draws)
{
  std::vector<int> result;
  int range = right_exclusive - left_inclusive; if (range <= 0) return result;
  size_t draws = std::min<size_t>(max_draws, static_cast<size_t>(range));
  if (draws == 0) return result;
  result.reserve(draws);
  for (size_t j = 1; j <= draws; ++j) {
    int pos = left_inclusive + static_cast<int>(std::floor((static_cast<double>(j) * range) / static_cast<double>(draws + 1)));
    if (pos < left_inclusive) pos = left_inclusive;
    if (pos >= right_exclusive) pos = right_exclusive - 1;
    if (!result.empty() && pos <= result.back()) pos = std::min(right_exclusive - 1, result.back() + 1);
    result.push_back(pos);
  }
  return result;
}

std::vector<int> sample_unique_ints_uniform_R(int left_inclusive, int right_exclusive, size_t k)
{
  std::vector<int> result; int range = right_exclusive - left_inclusive; if (range <= 0) return result;
  k = std::min<size_t>(k, static_cast<size_t>(range));
  if (k == 0) return result;
  if (k * 4 >= static_cast<size_t>(range)) {
    std::vector<int> all(range); std::iota(all.begin(), all.end(), left_inclusive);
    for (size_t i = 0; i < k; ++i) {
      size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i)); if (j >= all.size()) j = all.size() - 1;
      std::swap(all[i], all[j]);
    }
    result.assign(all.begin(), all.begin() + static_cast<long>(k));
    std::sort(result.begin(), result.end());
    return result;
  }
  std::unordered_set<int> used; result.reserve(k);
  while (result.size() < k) {
    int s = rng_randint(left_inclusive, right_exclusive); if (s >= right_exclusive) s = right_exclusive - 1;
    if (used.insert(s).second) result.push_back(s);
  }
  std::sort(result.begin(), result.end());
  return result;
}

} // namespace rpf_utils


