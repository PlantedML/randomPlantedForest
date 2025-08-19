#include "rpf.hpp"
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <random>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <thread>

// Internal helpers to reduce duplication across modes (no behavior changes)
namespace {
  // Thread-local RNG pointer used in worker threads for reproducible randomness
  thread_local std::mt19937_64* tls_rng_ptr = nullptr;
  inline double rng_runif01()
  {
    if (tls_rng_ptr) {
      // 53-bit precision uniform in [0,1)
      return std::generate_canonical<double, 53>(*tls_rng_ptr);
    }
    // Fallback: deterministic static RNG (used only if caller forgot to set tls_rng_ptr)
    static thread_local std::mt19937_64 fallback_rng(0x9E3779B97F4A7C15ULL);
    return std::generate_canonical<double, 53>(fallback_rng);
  }
  inline double rng_runif(double a, double b)
  {
    double u = rng_runif01();
    return a + u * (b - a);
  }
  inline int rng_randint(int left_inclusive, int right_exclusive)
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
  // Ensure per-leaf stable order and its sorted projection along feature k.
  // Reuses per-leaf caches when available and returns the materialized buffers.
  inline void ensure_order_and_sorted_vals_for_leaf(
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

  // Compute unique values from a sorted buffer (stable, no side effects)
  inline std::vector<double> compute_unique_sorted_values(const std::vector<double> &sorted_vals)
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

  // Build prefix sums and totals along a precomputed order
  inline void build_prefix_and_total_given_order(
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

  // Finalize the chosen split by building index sets and means from stored sums
  inline void finalize_split_from_sums(
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

  // Weighted sampling without replacement using Efraimidis–Spirakis with R RNG
  inline std::vector<size_t> sample_weighted_indices_filtered(
      const std::vector<double> &weights,
      size_t n_candidates)
  {
    std::vector<size_t> pos_idx; pos_idx.reserve(weights.size());
    std::vector<double> pos_w;   pos_w.reserve(weights.size());
    for (size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0.0) { pos_idx.push_back(i); pos_w.push_back(weights[i]); }
    const size_t P = pos_idx.size();
    std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
    if (P == 0) {
      // uniform fallback over all indices using partial Fisher–Yates via R RNG
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

  inline std::vector<int> compute_even_spread_indices(int left_inclusive, int right_exclusive, size_t max_draws)
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

  inline std::vector<int> sample_unique_ints_uniform_R(int left_inclusive, int right_exclusive, size_t k)
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

  // Thread-local cache for histogram mode per working set (per tree-family build)
  // Avoids races on the class member when building families in parallel
  thread_local std::vector<std::vector<int>> tls_working_bin_id;

  // Age helper for both SplitCandidate and ResultingTreeCandidate pools
  template <typename CandidateT>
  inline void age_pool_by_sample(const std::vector<size_t> &sample_idxs, int best_idx, std::vector<CandidateT> &pool)
  {
    for (size_t idx : sample_idxs) {
      if (static_cast<int>(idx) != best_idx) pool[idx].age += 1.0; else pool[idx].age = 0.0;
    }
  }
}

// Utilities shared across modes
bool RandomPlantedForest::possibleExists(
    int dim,
    const std::vector<SplitCandidate>& possible_splits,
    const std::set<int>& resulting_dims)
{
  for (const auto& c : possible_splits) {
    if (c.dim == dim && c.tree && c.tree->split_dims == resulting_dims)
      return true;
  }
  return false;
}

bool RandomPlantedForest::leafCandidateExists(
    const std::vector<SplitCandidate>& possible_splits,
    const std::shared_ptr<DecisionTree>& tree,
    size_t leaf_idx,
    int dim)
{
  for (const auto& c : possible_splits) {
    if (c.dim == dim && c.tree.get() == tree.get() && c.leaf_idx == leaf_idx)
      return true;
  }
  return false;
}

bool RandomPlantedForest::is_purified()
{
  return purified;
}

void RandomPlantedForest::L2_loss(Split &split)
{
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.min_sum = 0;
  for (size_t p = 0; p < value_size; ++p)
  {
    const double Ms = split.M_s[p];
    const double Mb = split.M_b[p];
    split.min_sum += -2 * Ms * split.sum_s[p] + split.I_s.size() * (Ms * Ms);
    split.min_sum += -2 * Mb * split.sum_b[p] + split.I_b.size() * (Mb * Mb);
  }
}

// constructor (parsing includes split_structure)
RandomPlantedForest::RandomPlantedForest(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                                         const NumericVector parameters)
{
  Rcpp::RNGScope scope;
  std::vector<double> pars = to_std_vec(parameters);
  if (pars.size() != 12 && pars.size() != 13)
  {
    Rcpp::stop("RandomPlantedForest requires 12 or 13 parameters, got %d", pars.size());
  }
  else
  {
    this->max_interaction = pars[0];
    this->n_trees = pars[1];
    this->n_splits = pars[2];
    this->split_try = pars[3];
    this->t_try = pars[4];
    this->purify_forest = pars[5];
    this->deterministic = pars[6];
    this->nthreads = pars[7];
    this->cross_validate = pars[8];
    this->split_decay_rate_ = pars[9];
    this->max_candidates_   = static_cast<size_t>(pars[10]);
    this->delete_leaves   = (pars[11] != 0);
    // map: 0=res_trees, 1=cur_trees_2, 2=cur_trees_1, 3=leaves, 4=hist
    this->split_structure_mode_ = (pars.size() >= 13) ? static_cast<int>(pars[12]) : 3;
  }
  this->set_data(samples_Y, samples_X);
}

// --------------- calcOptimalSplit per mode ---------------

// Mode 3: leaves (per-leaf candidates)

Split RandomPlantedForest::calcOptimalSplit_leaves(const std::vector<std::vector<double>> &Y,
                                                   const std::vector<std::vector<double>> &X,
                                                   std::vector<SplitCandidate> &possible_splits,
                                                   TreeFamily &curr_family)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  if (possible_splits.empty()) return min_split;

  // Sample candidate indices by age-weight, capped by max_candidates_
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs = this->deterministic ? std::vector<size_t>() : sample_weighted_indices_filtered(weights, n_candidates);
  if (this->deterministic) { for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i) sample_idxs.push_back(i); }

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx);
    int k = it->dim - 1;
    if (!it->tree || it->leaf_idx >= it->tree->leaves.size()) continue;
    Leaf* leafPtr = &it->tree->leaves[it->leaf_idx];

    const int leaf_size = this->n_leaves[k];
    const size_t m = leafPtr->individuals.size();
    if (m == 0) continue;

    // Build or reuse per-leaf caches and derive unique thresholds
    std::vector<size_t> order; std::vector<double> sorted_vals;
    ensure_order_and_sorted_vals_for_leaf(X, *leafPtr, k, order, sorted_vals);
    std::vector<double> unique = compute_unique_sorted_values(sorted_vals);

    if (unique.size() < 2 * static_cast<size_t>(leaf_size)) continue; // not enough room for both sides

    // Choose sample threshold indices within [leaf_size, unique.size() - leaf_size)
    std::vector<int> samples;
    int left = leaf_size; int right = (int)unique.size() - leaf_size;
    samples = this->deterministic ? compute_even_spread_indices(left, right, (size_t)this->split_try)
                                  : sample_unique_ints_uniform_R(left, right, (size_t)this->split_try);

    // Always use prefix sums
    std::vector<std::vector<double>> prefix; // [p][i]
    std::vector<double> total;               // [p]
    build_prefix_and_total_given_order(Y, *leafPtr, order, this->value_size, prefix, total);

    for (size_t si = 0; si < samples.size(); ++si) {
      const double sp = unique[samples[si]];
      const size_t pos = static_cast<size_t>(std::lower_bound(sorted_vals.begin(), sorted_vals.end(), sp) - sorted_vals.begin());
      if (pos == 0 || pos >= m) continue;
      if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size)) continue;

      double loss = 0.0;
      for (size_t p = 0; p < this->value_size; ++p) {
        const double sum_s_base = prefix[p][pos - 1];
        const double sum_b_base = total[p] - sum_s_base;
        loss -= (sum_s_base * sum_s_base) / static_cast<double>(pos);
        loss -= (sum_b_base * sum_b_base) / static_cast<double>(m - pos);
      }

      if (loss < min_split.min_sum) {
        min_split.min_sum = loss;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        min_split.split_point = sp;
        best_idx = (int)idx;
        // Prepare sums for Ms/Mb later via prefix
        min_split.sum_s.assign(this->value_size, 0.0);
        min_split.sum_b.assign(this->value_size, 0.0);
        for (size_t p = 0; p < this->value_size; ++p) {
          const double sum_s_base = prefix[p][pos - 1];
          const double sum_b_base = total[p] - sum_s_base;
          min_split.sum_s[p] = sum_s_base;
          min_split.sum_b[p] = sum_b_base;
        }
      }
    }
  }

  // Age candidates and finalize winner
  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);

  return min_split;
}

// Mode 1: cur_trees_2 (try every leaf within current/pred trees, sampling split points incrementally)
Split RandomPlantedForest::calcOptimalSplit_curTrees2(const std::vector<std::vector<double>> &Y,
                                                      const std::vector<std::vector<double>> &X,
                                                      std::vector<SplitCandidate> &possible_splits,
                                                      TreeFamily &curr_family)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));

  std::vector<double> weights(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!this->deterministic) {
    
    std::vector<size_t> pos_idx; pos_idx.reserve(possible_splits.size());
    std::vector<double> pos_w;   pos_w.reserve(possible_splits.size());
    for (size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0.0) { pos_idx.push_back(i); pos_w.push_back(weights[i]); }
    const size_t P = pos_idx.size();
    if (P == 0) {
      std::vector<size_t> all(possible_splits.size()); std::iota(all.begin(), all.end(), 0);
      size_t k = std::min<size_t>(n_candidates, all.size());
      for (size_t i = 0; i < k; ++i) { size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i)); if (j >= all.size()) j = all.size() - 1; std::swap(all[i], all[j]); }
      for (size_t i = 0; i < k; ++i) sample_idxs.push_back(all[i]);
    } else if (n_candidates * 8 < P) {
      // Use filtered positives + E-S keys path below for consistency and speed
      size_t k2 = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k2 < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k2, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k2); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    } else {
      size_t k = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    }
  } else { for (size_t i=0;i<n_candidates && i<possible_splits.size();++i) sample_idxs.push_back(i); }

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx);
    int k = it->dim - 1;
    int leaf_size = this->n_leaves[k];

    std::set<int> tree_dims = it->tree->split_dims;
    tree_dims.erase(k + 1); tree_dims.erase(0);

    std::vector<std::shared_ptr<DecisionTree>> curr_trees;
    if (tree_dims.empty()) {
      auto itZero = curr_family.find(std::set<int>{0});
      if (itZero != curr_family.end() && itZero->second) curr_trees.push_back(itZero->second);
    }
    if (auto itS = curr_family.find(tree_dims); itS != curr_family.end() && itS->second) curr_trees.push_back(itS->second);
    if (auto itD = curr_family.find(it->tree->split_dims); itD != curr_family.end() && itD->second) {
      if (curr_trees.empty() || curr_trees.back().get() != itD->second.get()) curr_trees.push_back(itD->second);
    }

    for (auto &curr_tree : curr_trees) {
      if (curr_tree->leaves.empty()) continue;
      for (auto &leaf : curr_tree->leaves) {
        // Reuse per-leaf-feature cache for order/sorted values; derive unique thresholds
        std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf;
        ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
        std::vector<double> unique_samples = compute_unique_sorted_values(sorted_vals_cf);
        if (unique_samples.size() < 2 * static_cast<size_t>(leaf_size)) continue;

        // Local prefix sums only if we need to evaluate more than one threshold
        const size_t m = leaf.individuals.size();
        std::vector<int> samples;
        if (this->deterministic) {
          int maxp = std::min<int>((int)unique_samples.size() - 1, 9);
          samples.resize(maxp); std::iota(samples.begin(), samples.end(), 1);
        } else {
          samples.resize(this->split_try);
          for (size_t i = 0; i < samples.size(); ++i) samples[i] = rng_randint(leaf_size, (int)unique_samples.size() - leaf_size);
          std::sort(samples.begin(), samples.end());
        }
        const bool single_eval = (samples.size() == 1);
        std::vector<std::vector<double>> prefix_cf; // [value_size][m]
        std::vector<double> total_cf;               // [value_size]
        if (!single_eval) build_prefix_and_total_given_order(Y, leaf, order_cf, this->value_size, prefix_cf, total_cf);

        for (size_t si = 0; si < samples.size(); ++si) {
          const double sp = unique_samples[samples[si]];
          // locate position in sorted values
          size_t pos = static_cast<size_t>(std::lower_bound(sorted_vals_cf.begin(), sorted_vals_cf.end(), sp) - sorted_vals_cf.begin());
          if (pos == 0 || pos >= m) continue;
          if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size)) continue;
          // loss computation: use prefix if available, else scan
          double loss = 0.0;
          if (!single_eval) {
            for (size_t p = 0; p < this->value_size; ++p) {
              const double sum_s_base = prefix_cf[p][pos - 1];
              const double sum_b_base = total_cf[p] - sum_s_base;
              loss -= (sum_s_base * sum_s_base) / static_cast<double>(pos);
              loss -= (sum_b_base * sum_b_base) / static_cast<double>(m - pos);
            }
          } else {
            // Fallback for split_try == 1: scan the leaf once for this threshold
            size_t ns = 0, nb = 0;
            std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
            for (int ind : leaf.individuals) {
              const bool left_side = (X[ind][k] < sp);
              if (left_side) {
                ++ns;
                for (size_t p = 0; p < this->value_size; ++p) { double v = Y[ind][p]; sum_s_adj[p] += v; }
              } else {
                ++nb;
                for (size_t p = 0; p < this->value_size; ++p) { double v = Y[ind][p]; sum_b_adj[p] += v; }
              }
            }
            if (ns == 0 || nb == 0) { continue; }
            for (size_t p = 0; p < this->value_size; ++p) {
              loss -= (sum_s_adj[p] * sum_s_adj[p]) / static_cast<double>(ns);
              loss -= (sum_b_adj[p] * sum_b_adj[p]) / static_cast<double>(nb);
            }
          }
          if (loss < min_split.min_sum) {
            min_split.min_sum = loss;
            min_split.tree_index = curr_tree;
            min_split.leaf_index = &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point = sp;
            best_idx = (int)idx;
            // store adjusted sums for Ms/Mb later
            min_split.sum_s.assign(this->value_size, 0.0);
            min_split.sum_b.assign(this->value_size, 0.0);
            if (!single_eval) {
              for (size_t p = 0; p < this->value_size; ++p) {
                const double sum_s_base = prefix_cf[p][pos - 1];
                const double sum_b_base = total_cf[p] - sum_s_base;
                min_split.sum_s[p] = sum_s_base;
                min_split.sum_b[p] = sum_b_base;
              }
            } else {
              // Recompute adjusted sums for this threshold via scan
              for (int ind : leaf.individuals) {
                if (X[ind][k] < sp) {
                  for (size_t p = 0; p < this->value_size; ++p) { double v = Y[ind][p]; min_split.sum_s[p] += v; }
                } else {
                  for (size_t p = 0; p < this->value_size; ++p) { double v = Y[ind][p]; min_split.sum_b[p] += v; }
                }
              }
            }
          }
        }
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}

// Mode 2: cur_trees_1 (pair-sampling within predecessor/current trees)
Split RandomPlantedForest::calcOptimalSplit_curTrees1(const std::vector<std::vector<double>> &Y,
                                                      const std::vector<std::vector<double>> &X,
                                                      std::vector<SplitCandidate> &possible_splits,
                                                      TreeFamily &curr_family)
{
  Split curr_split, min_split; min_split.min_sum = std::numeric_limits<double>::infinity(); curr_split.Y = &Y;

  unsigned int raw = (unsigned int)std::ceil(this->t_try * possible_splits.size());
  unsigned int upper = std::min<unsigned int>((unsigned int)this->max_candidates_, (unsigned int)possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw, upper));
  std::vector<double> weights(possible_splits.size());

  for (size_t i = 0; i < possible_splits.size(); ++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);

  size_t positive_count_ = 0; for (double w : weights) if (w > 0.0) ++positive_count_;

  if (positive_count_ == 0) { n_candidates = 1; }
  else { if (n_candidates > positive_count_) n_candidates = static_cast<unsigned int>(positive_count_); }
  
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!this->deterministic) {
    std::vector<size_t> pos_idx; pos_idx.reserve(possible_splits.size());
    std::vector<double> pos_w;   pos_w.reserve(possible_splits.size());
    for (size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0.0) { pos_idx.push_back(i); pos_w.push_back(weights[i]); }
    const size_t P = pos_idx.size();
    if (P == 0) {
      std::vector<size_t> all(possible_splits.size()); std::iota(all.begin(), all.end(), 0);
      size_t k = std::min<size_t>(n_candidates, all.size());
      for (size_t i = 0; i < k; ++i) { size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i)); if (j >= all.size()) j = all.size() - 1; std::swap(all[i], all[j]); }
      for (size_t i = 0; i < k; ++i) sample_idxs.push_back(all[i]);
    } else if (n_candidates * 8 < P) {
      size_t k2 = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k2 < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k2, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k2); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    } else {
      size_t k = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    }
  } else { for (size_t i=0;i<n_candidates && i<possible_splits.size();++i) sample_idxs.push_back(i); }
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx); int k = it->dim - 1; int leaf_size = this->n_leaves[k];
    std::set<int> Dprime_minus_k = it->tree->split_dims; Dprime_minus_k.erase(k + 1); Dprime_minus_k.erase(0);
    std::vector<std::shared_ptr<DecisionTree>> sources; sources.reserve(2);
    if (Dprime_minus_k.empty()) { if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end()) sources.push_back(itZero->second); }
    else { if (auto itS = curr_family.find(Dprime_minus_k); itS != curr_family.end()) sources.push_back(itS->second); }
    if (auto itD = curr_family.find(it->tree->split_dims); itD != curr_family.end()) if (sources.empty() || sources.back().get() != itD->second.get()) sources.push_back(itD->second);

    // Fast path when split_try is small relative to number of feasible leaves
    // Heuristic: use fast path if split_try <= max(10, floor(1 * L)), where L is feasible leaf count
    if (!this->deterministic) {
      struct LeafEnv { std::shared_ptr<DecisionTree> tree; Leaf* leaf; int left; int right_env; };
      std::vector<LeafEnv> env; env.reserve(sources.size()*8);
      std::vector<double> w_env; w_env.reserve(sources.size()*8);
      // Build envelope weights W = max(0, m - 2*leaf_size)
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty()) continue;
        for (auto &leaf : src_tree->leaves) {
          const int m = static_cast<int>(leaf.individuals.size());
          const int left = leaf_size;
          const int right_env = m - leaf_size;
          const int width = right_env - left;
          if (width > 0) { env.push_back({src_tree, &leaf, left, right_env}); w_env.push_back(static_cast<double>(width)); }
        }
      }
      if (!env.empty()) {
        // Weighted sampling WITH replacement using a Fenwick tree over envelope weights,
        // updating weights after every draw (picked leaf loses one feasible threshold).
        struct Fenwick {
          std::vector<double> bit; size_t n;
          explicit Fenwick(size_t n_) : bit(n_ + 1, 0.0), n(n_) {}
          void add(size_t idx, double delta) { size_t i = idx + 1; while (i <= n) { bit[i] += delta; i += i & (~i + 1); } }
          void build(const std::vector<double>& w) { for (size_t i = 0; i < n; ++i) if (w[i] != 0.0) add(i, w[i]); }
          size_t lower_bound(double target) const {
            size_t idx = 0; double acc = 0.0; size_t bitMask = 1; while ((bitMask << 1) <= n) bitMask <<= 1;
            for (size_t step = bitMask; step != 0; step >>= 1) {
              size_t next = idx + step; if (next <= n && acc + bit[next] <= target) { acc += bit[next]; idx = next; }
            }
            return (idx >= n ? n - 1 : idx);
          }
        };
        const size_t L = env.size();
        std::vector<double> w_curr(L, 0.0); double total_w = 0.0;
        for (size_t i = 0; i < L; ++i) { double w = w_env[i]; if (w > 0.0) { w_curr[i] = w; total_w += w; } }
        if (total_w <= 0.0) { continue; }
        Fenwick ft(L); ft.build(w_curr);
        struct DynLeaf {
          std::vector<double> unique;
          std::vector<char> used_flags;
          int left=0; int right_true=0; int remaining=-1; bool initialized=false;
          // caches to enable O(1) threshold evaluations
          std::vector<size_t> order_cf;
          std::vector<double> sorted_vals;
          std::vector<std::vector<double>> prefix_cf; // [value_size][m]
          std::vector<double> total_cf;               // [value_size]
        };
        std::vector<DynLeaf> dyn(L);
        for (size_t i = 0; i < L; ++i) dyn[i].left = env[i].left;
        size_t draws = std::min((size_t)this->split_try, (size_t)std::floor(total_w));
        for (size_t t = 0; t < draws; ++t) {
          if (total_w <= 0.0) break;
          // sample leaf index proportional to current weights using R RNG
          double u = rng_runif01(); if (u >= 1.0) u = 1.0 - std::numeric_limits<double>::min();
          double r = u * total_w; size_t pos = ft.lower_bound(r);
          // lazily initialize selected leaf's unique thresholds and correct its weight
          if (!dyn[pos].initialized) {
            auto &E = env[pos];
            if (E.leaf->order_cache.count(k) && E.leaf->order_cache[k].size() == E.leaf->individuals.size()) { dyn[pos].order_cf = E.leaf->order_cache[k]; }
            else { dyn[pos].order_cf.resize(E.leaf->individuals.size()); std::iota(dyn[pos].order_cf.begin(), dyn[pos].order_cf.end(), 0); std::stable_sort(dyn[pos].order_cf.begin(), dyn[pos].order_cf.end(), [&](size_t a, size_t b){ return X[E.leaf->individuals[a]][k] < X[E.leaf->individuals[b]][k]; }); E.leaf->order_cache[k] = dyn[pos].order_cf; }
            if (E.leaf->sorted_vals_cache.count(k) && E.leaf->sorted_vals_cache[k].size() == E.leaf->individuals.size()) { dyn[pos].sorted_vals = E.leaf->sorted_vals_cache[k]; }
            else { dyn[pos].sorted_vals.resize(E.leaf->individuals.size()); for (size_t ii=0; ii<E.leaf->individuals.size(); ++ii) dyn[pos].sorted_vals[ii] = X[E.leaf->individuals[dyn[pos].order_cf[ii]]][k]; E.leaf->sorted_vals_cache[k] = dyn[pos].sorted_vals; }
            dyn[pos].unique.clear(); dyn[pos].unique.reserve(dyn[pos].sorted_vals.size()); if (!dyn[pos].sorted_vals.empty()) { dyn[pos].unique.push_back(dyn[pos].sorted_vals[0]); for (size_t ii=1; ii<dyn[pos].sorted_vals.size(); ++ii) if (dyn[pos].sorted_vals[ii] != dyn[pos].unique.back()) dyn[pos].unique.push_back(dyn[pos].sorted_vals[ii]); }
            dyn[pos].right_true = static_cast<int>(dyn[pos].unique.size()) - dyn[pos].left;
            int corrected = std::max(0, dyn[pos].right_true - dyn[pos].left);
            dyn[pos].remaining = corrected;
            dyn[pos].used_flags.assign((size_t)std::max(0, corrected), 0);
            dyn[pos].initialized = true;
            // build prefix sums once for this leaf along the cached order
            const size_t m_build = E.leaf->individuals.size();
            dyn[pos].prefix_cf.assign(this->value_size, std::vector<double>(m_build, 0.0));
            for (size_t p = 0; p < this->value_size; ++p) {
              double accv = 0.0; for (size_t ii = 0; ii < m_build; ++ii) { accv += Y[E.leaf->individuals[dyn[pos].order_cf[ii]]][p]; dyn[pos].prefix_cf[p][ii] = accv; }
            }
            dyn[pos].total_cf.assign(this->value_size, 0.0);
            for (size_t p = 0; p < this->value_size; ++p) dyn[pos].total_cf[p] = dyn[pos].prefix_cf[p][m_build - 1];
            // adjust Fenwick weight from approximate to true remaining
            double delta = (double)corrected - w_curr[pos]; if (delta != 0.0) { ft.add(pos, delta); total_w += delta; w_curr[pos] += delta; if (total_w <= 0.0) break; }
            if (dyn[pos].remaining <= 0) { --t; continue; } // resample this draw if leaf has no feasible threshold
          }
          if (dyn[pos].remaining <= 0) { --t; continue; }
          auto &E = env[pos];
          // sample an unused threshold uniformly within [left, right_true)
          int left = dyn[pos].left; int right_true = dyn[pos].right_true;
          int s_idx; do { s_idx = rng_randint(left, right_true); if (s_idx >= right_true) s_idx = right_true - 1; } while (dyn[pos].used_flags[(size_t)(s_idx - left)]);
          dyn[pos].used_flags[(size_t)(s_idx - left)] = 1; dyn[pos].remaining -= 1; ft.add(pos, -1.0); w_curr[pos] -= 1.0; total_w -= 1.0;
          double sp = dyn[pos].unique[(size_t)s_idx];
          // Evaluate split via prefix sums (preserve sampling behavior, avoid rescans)
          const size_t m_eval = E.leaf->individuals.size();
          size_t pos_in_sorted = static_cast<size_t>(std::lower_bound(dyn[pos].sorted_vals.begin(), dyn[pos].sorted_vals.end(), sp) - dyn[pos].sorted_vals.begin());
          if (pos_in_sorted == 0 || pos_in_sorted >= m_eval) { if (total_w <= 0.0) break; else continue; }
          double loss = 0.0;
          std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
          for (size_t p = 0; p < this->value_size; ++p) {
            const double sum_s_base = dyn[pos].prefix_cf[p][pos_in_sorted - 1];
            const double sum_b_base = dyn[pos].total_cf[p] - sum_s_base;
            sum_s_adj[p] = sum_s_base; sum_b_adj[p] = sum_b_base;
            loss -= (sum_s_adj[p] * sum_s_adj[p]) / static_cast<double>(pos_in_sorted);
            loss -= (sum_b_adj[p] * sum_b_adj[p]) / static_cast<double>(m_eval - pos_in_sorted);
          }
          if (loss < min_split.min_sum) {
            min_split.min_sum = loss;
            min_split.tree_index = E.tree; min_split.leaf_index = E.leaf; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx;
            min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj;
          }
          if (total_w <= 0.0) break;
        }
        continue; // proceed to next candidate idx using fast path only
      }
    }

    struct LeafUnit {
      std::shared_ptr<DecisionTree> tree; Leaf* leaf; std::vector<double> unique;
      int left=0,right=0,remaining=0; int used_count=0; std::vector<char> used_flags;
      // caches to enable O(1) threshold evaluations
      std::vector<size_t> order_cf;
      std::vector<double> sorted_vals;
      std::vector<std::vector<double>> prefix_cf; // [value_size][m]
      std::vector<double> total_cf;               // [value_size]
    };
    std::vector<LeafUnit> units; units.reserve(sources.size()*8);
    for (const auto &src_tree : sources) {
      if (src_tree->leaves.empty()) continue; for (auto &leaf : src_tree->leaves) {
        LeafUnit u; u.tree = src_tree; u.leaf = &leaf; u.left = (int)leaf_size;
        // Use caches to estimate number of unique thresholds without materializing
        std::vector<size_t> order_cf;
        if (leaf.order_cache.count(k) && leaf.order_cache[k].size() == leaf.individuals.size()) { order_cf = leaf.order_cache[k]; }
        else { order_cf.resize(leaf.individuals.size()); std::iota(order_cf.begin(), order_cf.end(), 0); std::stable_sort(order_cf.begin(), order_cf.end(), [&](size_t a, size_t b){ return X[leaf.individuals[a]][k] < X[leaf.individuals[b]][k]; }); leaf.order_cache[k] = order_cf; }
        std::vector<double> sorted_vals_cf;
        if (leaf.sorted_vals_cache.count(k) && leaf.sorted_vals_cache[k].size() == leaf.individuals.size()) { sorted_vals_cf = leaf.sorted_vals_cache[k]; }
        else { sorted_vals_cf.resize(leaf.individuals.size()); for (size_t i=0;i<leaf.individuals.size();++i) sorted_vals_cf[i] = X[leaf.individuals[order_cf[i]]][k]; leaf.sorted_vals_cache[k] = sorted_vals_cf; }
        size_t unique_count = 0; if (!sorted_vals_cf.empty()) { unique_count = 1; for (size_t i=1;i<sorted_vals_cf.size(); ++i) if (sorted_vals_cf[i] != sorted_vals_cf[i-1]) ++unique_count; }
        u.right = (int)unique_count - (int)leaf_size; if (u.right > u.left) { u.remaining = u.right - u.left; u.used_count = 0; /* used_flags and unique filled lazily */ units.push_back(std::move(u)); }
      }
    }
    size_t total_remaining = 0; for (auto &u : units) total_remaining += (size_t)u.remaining; if (total_remaining == 0) continue;
    size_t draws = std::min((size_t)this->split_try, total_remaining);
    auto draw_leaf_index = [&](double r)->size_t { size_t acc=0; for (size_t i=0;i<units.size();++i){ if (units[i].remaining<=0) continue; acc += (size_t)units[i].remaining; if (r < (double)acc) return i; } for (size_t i=0;i<units.size();++i) if (units[i].remaining>0) return i; return units.size()-1; };
    for (size_t t=0; t<draws; ++t) {
      size_t leaf_i; if (this->deterministic) { double step = (double)total_remaining / (double)draws; double target = step * (t + 0.5); if (target >= (double)total_remaining) target = (double)(total_remaining - 1); leaf_i = draw_leaf_index(target); }
      else { double r = rng_runif(0.0, (double)total_remaining); leaf_i = draw_leaf_index(r); }
      auto &sel = units[leaf_i];
      // Lazily materialize unique thresholds and used flags for selected leaf
      if (sel.unique.empty()) {
        if (sel.leaf->order_cache.count(k) && sel.leaf->order_cache[k].size() == sel.leaf->individuals.size()) { sel.order_cf = sel.leaf->order_cache[k]; }
        else { sel.order_cf.resize(sel.leaf->individuals.size()); std::iota(sel.order_cf.begin(), sel.order_cf.end(), 0); std::stable_sort(sel.order_cf.begin(), sel.order_cf.end(), [&](size_t a, size_t b){ return X[sel.leaf->individuals[a]][k] < X[sel.leaf->individuals[b]][k]; }); sel.leaf->order_cache[k] = sel.order_cf; }
        if (sel.leaf->sorted_vals_cache.count(k) && sel.leaf->sorted_vals_cache[k].size() == sel.leaf->individuals.size()) { sel.sorted_vals = sel.leaf->sorted_vals_cache[k]; }
        else { sel.sorted_vals.resize(sel.leaf->individuals.size()); for (size_t i=0;i<sel.leaf->individuals.size();++i) sel.sorted_vals[i] = X[sel.leaf->individuals[sel.order_cf[i]]][k]; sel.leaf->sorted_vals_cache[k] = sel.sorted_vals; }
        sel.unique.clear(); sel.unique.reserve(sel.sorted_vals.size()); if (!sel.sorted_vals.empty()) { sel.unique.push_back(sel.sorted_vals[0]); for (size_t i=1;i<sel.sorted_vals.size(); ++i) if (sel.sorted_vals[i] != sel.unique.back()) sel.unique.push_back(sel.sorted_vals[i]); }
        if (sel.used_flags.empty()) sel.used_flags.assign((size_t)(sel.right - sel.left), 0);
        // build prefix sums once for this leaf along the cached order
        const size_t mloc_build = sel.leaf->individuals.size();
        sel.prefix_cf.assign(this->value_size, std::vector<double>(mloc_build, 0.0));
        for (size_t p = 0; p < this->value_size; ++p) {
          double accv = 0.0; for (size_t ii = 0; ii < mloc_build; ++ii) { accv += Y[sel.leaf->individuals[sel.order_cf[ii]]][p]; sel.prefix_cf[p][ii] = accv; }
        }
        sel.total_cf.assign(this->value_size, 0.0);
        for (size_t p = 0; p < this->value_size; ++p) sel.total_cf[p] = sel.prefix_cf[p][mloc_build - 1];
      } else if (sel.used_flags.empty()) { sel.used_flags.assign((size_t)(sel.right - sel.left), 0); }
      int s_idx; if (this->deterministic) {
        int range = sel.right - sel.left; int guess = sel.left + (int)(((double)sel.used_count+0.5)/((double)sel.remaining+0.5) * range); if (guess >= sel.right) guess = sel.right - 1; int lo=guess, hi=guess; bool found=false; while (lo>=sel.left || hi<sel.right) { if (lo>=sel.left && !sel.used_flags[lo - sel.left]) { s_idx = lo; found=true; break; } if (hi<sel.right && !sel.used_flags[hi - sel.left]) { s_idx = hi; found=true; break; } --lo; ++hi; } if (!found) for (int p=sel.left;p<sel.right;++p) if (!sel.used_flags[p - sel.left]) { s_idx=p; break; }
      } else { do { s_idx = rng_randint(sel.left, sel.right); if (s_idx >= sel.right) s_idx = sel.right - 1; } while (sel.used_flags[s_idx - sel.left]); }
      sel.used_flags[s_idx - sel.left] = 1; sel.used_count += 1; sel.remaining -= 1; total_remaining -= 1;
      double sp = sel.unique[s_idx];
      // Evaluate via prefix sums (preserve sampling behavior, avoid rescans)
      const size_t mloc2 = sel.leaf->individuals.size();
      size_t pos_in_sorted2 = static_cast<size_t>(std::lower_bound(sel.sorted_vals.begin(), sel.sorted_vals.end(), sp) - sel.sorted_vals.begin());
      if (pos_in_sorted2 < 1 || pos_in_sorted2 >= mloc2) { if (total_remaining == 0) break; else continue; }
      double loss2 = 0.0;
      std::vector<double> sum_s_adj2(this->value_size, 0.0), sum_b_adj2(this->value_size, 0.0);
      for (size_t p = 0; p < this->value_size; ++p) {
        const double sum_s_base2 = sel.prefix_cf[p][pos_in_sorted2 - 1];
        const double sum_b_base2 = sel.total_cf[p] - sum_s_base2;
        sum_s_adj2[p] = sum_s_base2; sum_b_adj2[p] = sum_b_base2;
        loss2 -= (sum_s_adj2[p] * sum_s_adj2[p]) / static_cast<double>(pos_in_sorted2);
        loss2 -= (sum_b_adj2[p] * sum_b_adj2[p]) / static_cast<double>(mloc2 - pos_in_sorted2);
      }
      if (loss2 < min_split.min_sum) { min_split.min_sum = loss2; min_split.tree_index = sel.tree; min_split.leaf_index = sel.leaf; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; min_split.sum_s = sum_s_adj2; min_split.sum_b = sum_b_adj2; }
      if (total_remaining == 0) break;
    }
  }
  for (size_t idx : sample_idxs) { if ((int)idx != best_idx) possible_splits[idx].age += 1.0; else possible_splits[idx].age = 0.0; }
  // finalize winner: build indices and means from stored sums
  if (!std::isinf(min_split.min_sum) && min_split.leaf_index != nullptr) {
    int kfin = min_split.split_coordinate - 1;
    auto &leaf_fin = *min_split.leaf_index;
    double sp_fin = min_split.split_point;
    min_split.I_s.clear(); min_split.I_b.clear();
    for (int ind : leaf_fin.individuals) {
      if (X[ind][kfin] < sp_fin) min_split.I_s.push_back(ind); else min_split.I_b.push_back(ind);
    }
    min_split.M_s.assign(this->value_size, 0.0);
    min_split.M_b.assign(this->value_size, 0.0);
    if (!min_split.I_s.empty()) for (size_t p = 0; p < this->value_size; ++p) min_split.M_s[p] = min_split.sum_s[p] / static_cast<double>(min_split.I_s.size());
    if (!min_split.I_b.empty()) for (size_t p = 0; p < this->value_size; ++p) min_split.M_b[p] = min_split.sum_b[p] / static_cast<double>(min_split.I_b.size());
  }
  return min_split;
}

// Mode 0: res_trees (operate on resulting trees pool)
bool RandomPlantedForest::resultingTreeExists(const std::vector<RandomPlantedForest::ResultingTreeCandidate>& pool, const std::set<int>& dims) {
  for (const auto &c : pool) if (c.tree->get_split_dims() == dims) return true; return false;
}

Split RandomPlantedForest::calcOptimalSplit_resTrees(const std::vector<std::vector<double>> &Y,
                                                     const std::vector<std::vector<double>> &X,
                                                     std::vector<ResultingTreeCandidate> &possible_trees,
                                                     TreeFamily &curr_family)
{
  Split curr_split, min_split; min_split.min_sum = std::numeric_limits<double>::infinity(); curr_split.Y = &Y;
  unsigned int raw = (unsigned int)std::ceil(this->t_try * possible_trees.size());
  unsigned int upper = std::min<unsigned int>((unsigned int)this->max_candidates_, (unsigned int)possible_trees.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw, upper));
  std::vector<double> weights(possible_trees.size()); for (size_t i=0;i<possible_trees.size();++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_trees[i].age);
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!this->deterministic) {
    std::vector<size_t> pos_idx; pos_idx.reserve(possible_trees.size());
    std::vector<double> pos_w;   pos_w.reserve(possible_trees.size());
    for (size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0.0) { pos_idx.push_back(i); pos_w.push_back(weights[i]); }
    const size_t P = pos_idx.size();
    if (P == 0) {
      std::vector<size_t> all(possible_trees.size()); std::iota(all.begin(), all.end(), 0);
      size_t k = std::min<size_t>(n_candidates, all.size());
      for (size_t i = 0; i < k; ++i) { size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i)); if (j >= all.size()) j = all.size() - 1; std::swap(all[i], all[j]); }
      for (size_t i = 0; i < k; ++i) sample_idxs.push_back(all[i]);
    } else if (n_candidates * 8 < P) {
      size_t k2 = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k2 < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k2, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k2); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    } else {
      size_t k = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double,size_t>> keys; keys.reserve(P);
      for (size_t i = 0; i < P; ++i) { double u = rng_runif01(); if (u <= 0.0) u = std::numeric_limits<double>::min(); double key = -std::log(u) / pos_w[i]; keys.emplace_back(key, pos_idx[i]); }
      if (k < keys.size()) { std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto& a, const auto& b){ return a.first < b.first; }); keys.resize(k); }
      for (auto &kv : keys) sample_idxs.push_back(kv.second);
    }
  } else { for (size_t i=0;i<n_candidates && i<possible_trees.size();++i) sample_idxs.push_back(i); }
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto itCand = possible_trees.begin(); std::advance(itCand, idx);
    const std::set<int>& Dprime = itCand->tree->split_dims;
    std::shared_ptr<DecisionTree> tree_Dprime; if (auto itD = curr_family.find(Dprime); itD != curr_family.end()) tree_Dprime = itD->second;
    for (int k_dim : Dprime) {
      const int k = k_dim - 1; int leaf_size = this->n_leaves[k]; std::vector<std::shared_ptr<DecisionTree>> sources; sources.reserve(2);
      std::set<int> S = Dprime; S.erase(k_dim);
      if (S.empty()) { if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end()) sources.push_back(itZero->second); }
      else { if (auto itS = curr_family.find(S); itS != curr_family.end()) sources.push_back(itS->second); }
      if (tree_Dprime) { if (sources.empty() || sources.back().get() != tree_Dprime.get()) sources.push_back(tree_Dprime); }
      struct LeafUnit {
        std::shared_ptr<DecisionTree> tree; Leaf* leaf;
        int left=0,right=0,remaining=0; int used_count=0; std::vector<char> used_flags;
        // Lazy caches for fast per-threshold evaluation (built on first use)
        bool initialized=false;
        std::vector<size_t> order_cf;
        std::vector<double> sorted_vals;
        std::vector<double> unique; // unique sorted feature values within this leaf for k
        std::vector<std::vector<double>> prefix_cf; // [value_size][m]
        std::vector<double> total_cf;               // [value_size]
      };
      std::vector<LeafUnit> units; units.reserve(sources.size()*8);
      for (const auto &src_tree : sources) { if (src_tree->leaves.empty()) continue; for (auto &leaf : src_tree->leaves) {
        LeafUnit u; u.tree = src_tree; u.leaf = &leaf;
        // Ensure per-leaf-feature order and sorted projection are cached
        std::vector<size_t> tmp_order; std::vector<double> tmp_sorted;
        ensure_order_and_sorted_vals_for_leaf(X, leaf, k, tmp_order, tmp_sorted);
        const auto &sorted_vals_cf_ref = leaf.sorted_vals_cache[k];
        size_t unique_count = 0; if (!sorted_vals_cf_ref.empty()) { unique_count = 1; for (size_t i=1;i<sorted_vals_cf_ref.size(); ++i) if (sorted_vals_cf_ref[i] != sorted_vals_cf_ref[i-1]) ++unique_count; }
        u.left = (int)leaf_size; u.right = (int)unique_count - (int)leaf_size; if (u.right > u.left) { u.remaining = u.right - u.left; u.used_count=0; /* used_flags filled lazily */ units.push_back(std::move(u)); }
      } }
      size_t total_remaining = 0; for (auto &u : units) total_remaining += (size_t)u.remaining; if (total_remaining == 0) continue; size_t draws = std::min((size_t)this->split_try, total_remaining);
      auto draw_leaf_index = [&](double r)->size_t { size_t acc=0; for (size_t i=0;i<units.size();++i){ if (units[i].remaining<=0) continue; acc += (size_t)units[i].remaining; if (r < (double)acc) return i; } for (size_t i=0;i<units.size();++i) if (units[i].remaining>0) return i; return units.size()-1; };
      for (size_t t=0; t<draws; ++t) {
        size_t leaf_i; if (this->deterministic) { double step = (double)total_remaining / (double)draws; double target = step * (t + 0.5); if (target >= (double)total_remaining) target = (double)(total_remaining - 1); leaf_i = draw_leaf_index(target); } else { double r = rng_runif(0.0, (double)total_remaining); leaf_i = draw_leaf_index(r); }
        auto &sel = units[leaf_i]; int s_idx; if (this->deterministic) {
          int range = sel.right - sel.left; int guess = sel.left + (int)(((double)sel.used_count+0.5)/((double)sel.remaining+0.5)*range); if (guess >= sel.right) guess = sel.right - 1; int lo=guess, hi=guess; bool found=false; while (lo>=sel.left || hi<sel.right) { if (lo>=sel.left && (!sel.used_flags.empty() ? !sel.used_flags[lo - sel.left] : true)) { s_idx=lo; found=true; break;} if (hi<sel.right && (!sel.used_flags.empty() ? !sel.used_flags[hi - sel.left] : true)) { s_idx=hi; found=true; break;} --lo; ++hi; } if (!found) for (int p=sel.left;p<sel.right;++p) if (sel.used_flags.empty() || !sel.used_flags[p - sel.left]) { s_idx=p; break; }
        } else { do { s_idx = rng_randint(sel.left, sel.right); } while (!sel.used_flags.empty() && sel.used_flags[s_idx - sel.left]); }
        if (sel.used_flags.empty()) sel.used_flags.assign((size_t)(sel.right - sel.left), 0);
        sel.used_flags[s_idx - sel.left] = 1; sel.used_count += 1; sel.remaining -= 1; total_remaining -= 1;

        // Lazily initialize caches for this leaf (order, sorted values, unique list, and prefix sums)
        if (!sel.initialized) {
          // order and sorted values
          ensure_order_and_sorted_vals_for_leaf(X, *sel.leaf, k, sel.order_cf, sel.sorted_vals);
          sel.unique.clear(); sel.unique.reserve(sel.sorted_vals.size()); if (!sel.sorted_vals.empty()) { sel.unique.push_back(sel.sorted_vals[0]); for (size_t i=1;i<sel.sorted_vals.size(); ++i) if (sel.sorted_vals[i] != sel.unique.back()) sel.unique.push_back(sel.sorted_vals[i]); }
          const size_t mloc = sel.leaf->individuals.size();
          build_prefix_and_total_given_order(Y, *sel.leaf, sel.order_cf, this->value_size, sel.prefix_cf, sel.total_cf);
          sel.initialized = true;
        }

        double sp = sel.unique[(size_t)s_idx];
        // Evaluate split via prefix sums
        const size_t m_eval = sel.leaf->individuals.size();
        size_t pos_in_sorted = static_cast<size_t>(std::lower_bound(sel.sorted_vals.begin(), sel.sorted_vals.end(), sp) - sel.sorted_vals.begin());
        if (pos_in_sorted == 0 || pos_in_sorted >= m_eval) { if (total_remaining == 0) break; else continue; }
        double loss = 0.0;
        std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
        for (size_t p = 0; p < this->value_size; ++p) {
          const double sum_s_base = sel.prefix_cf[p][pos_in_sorted - 1];
          const double sum_b_base = sel.total_cf[p] - sum_s_base;
          // Offsets disabled here as in original path
          sum_s_adj[p] = sum_s_base; sum_b_adj[p] = sum_b_base;
          loss -= (sum_s_adj[p] * sum_s_adj[p]) / static_cast<double>(pos_in_sorted);
          loss -= (sum_b_adj[p] * sum_b_adj[p]) / static_cast<double>(m_eval - pos_in_sorted);
        }
        if (loss < min_split.min_sum) {
          min_split.min_sum = loss;
          min_split.tree_index = sel.tree; min_split.leaf_index = sel.leaf; min_split.split_coordinate = k_dim; min_split.split_point = sp; best_idx = (int)idx;
          min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj;
        }
        if (total_remaining == 0) break;
      }
    }
  }
  age_pool_by_sample(sample_idxs, best_idx, possible_trees);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}

// Mode 4: histogram-binned evaluation (per-leaf candidates like mode 3)
Split RandomPlantedForest::calcOptimalSplit_hist(const std::vector<std::vector<double>> &Y,
                                                 const std::vector<std::vector<double>> &X,
                                                 std::vector<SplitCandidate> &possible_splits,
                                                 TreeFamily &curr_family)
{
  Split min_split; min_split.min_sum = std::numeric_limits<double>::infinity();
  if (possible_splits.empty()) return min_split;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs = this->deterministic ? std::vector<size_t>() : sample_weighted_indices_filtered(weights, n_candidates);
  if (this->deterministic) { for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i) sample_idxs.push_back(i); }

  // Use per-feature effective bin count based on actual cut count for stability
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx);
    if (!it->tree || it->leaf_idx >= it->tree->leaves.size()) continue;
    const int k_dim = it->dim; // 1-based
    const int k = k_dim - 1;
    Leaf* leafPtr = &it->tree->leaves[it->leaf_idx];
    const int leaf_min = this->n_leaves[k];
    const size_t m = leafPtr->individuals.size();
    if (m == 0) continue;

    // Build histogram for this leaf and feature k using cached working bin ids
    const auto &cuts_k = (k >= 0 && k < (int)feature_cut_points_.size()) ? feature_cut_points_[k] : std::vector<double>{};
    size_t Kf = cuts_k.size() + 1; if (Kf < 2) continue; // cannot split without at least 2 bins
    std::vector<int> cnt(Kf, 0);
    std::vector<std::vector<double>> sum(Kf, std::vector<double>(this->value_size, 0.0));
    const bool have_cached = (split_structure_mode_ == 4) && ((size_t)k < tls_working_bin_id.size());
    if (have_cached) {
      const std::vector<int> &bin_k = tls_working_bin_id[(size_t)k];
      for (int ind : leafPtr->individuals) {
        int b = bin_k[(size_t)ind];
        cnt[(size_t)b] += 1;
        for (size_t p = 0; p < this->value_size; ++p) sum[(size_t)b][p] += Y[ind][p];
      }
    } else {
      for (int ind : leafPtr->individuals) {
        double v = X[ind][k];
        int b = 0;
        if (!cuts_k.empty()) {
          auto itb = std::upper_bound(cuts_k.begin(), cuts_k.end(), v);
          b = (int)std::distance(cuts_k.begin(), itb);
          if (b < 0) b = 0; if ((size_t)b >= Kf) b = (int)Kf - 1;
        }
        cnt[(size_t)b] += 1;
        for (size_t p = 0; p < this->value_size; ++p) sum[(size_t)b][p] += Y[ind][p];
      }
    }

    // Single sweep over bin boundaries, no extra prefix storage
    const int total_n = (int)m;
    std::vector<double> total_sum(this->value_size, 0.0);
    for (size_t b = 0; b < Kf; ++b) {
      for (size_t p = 0; p < this->value_size; ++p) total_sum[p] += sum[b][p];
    }
    int left_n = 0;
    std::vector<double> left_sum(this->value_size, 0.0);
    for (size_t b_left = 0; b_left + 1 <= Kf - 1; ++b_left) {
      // Move boundary after bin b_left: left includes bins [0..b_left]
      left_n += cnt[b_left];
      for (size_t p = 0; p < this->value_size; ++p) left_sum[p] += sum[b_left][p];
      int right_n = total_n - left_n;
      if (left_n < leaf_min || right_n < leaf_min) continue;
      double loss = 0.0;
      for (size_t p = 0; p < this->value_size; ++p) {
        double ls = left_sum[p];
        double rs = total_sum[p] - ls;
        loss -= (ls * ls) / (double)left_n;
        loss -= (rs * rs) / (double)right_n;
      }
      if (loss < min_split.min_sum) {
        min_split.min_sum = loss;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        // Boundary index is b_left+1 in terms of bins on the right side
        double sp = 0.0;
        if (k >= 0 && k < (int)feature_cut_points_.size() && !feature_cut_points_[k].empty()) {
          const auto &cuts = feature_cut_points_[k];
          size_t cp_idx = (size_t)std::min<size_t>(b_left, cuts.size() - 1);
          sp = cuts[cp_idx];
        } else {
          sp = 0.5 * (leafPtr->intervals[k].first + leafPtr->intervals[k].second);
        }
        min_split.split_point = sp;
        best_idx = (int)idx;
        // Store sums for this boundary
        min_split.sum_s.assign(this->value_size, 0.0);
        min_split.sum_b.assign(this->value_size, 0.0);
        for (size_t p = 0; p < this->value_size; ++p) { min_split.sum_s[p] = left_sum[p]; min_split.sum_b[p] = total_sum[p] - left_sum[p]; }
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}

// Dispatcher used by create_tree_family
Split RandomPlantedForest::calcOptimalSplit(const std::vector<std::vector<double>> &Y,
                                            const std::vector<std::vector<double>> &X,
                                            std::vector<SplitCandidate> &possible_splits,
                                            TreeFamily &curr_family)
{
  if (split_structure_mode_ == 3) {
    return this->calcOptimalSplit_leaves(Y, X, possible_splits, curr_family);
  } else if (split_structure_mode_ == 2) {
    return this->calcOptimalSplit_curTrees1(Y, X, possible_splits, curr_family);
  } else if (split_structure_mode_ == 1) {
    return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family);
  } else if (split_structure_mode_ == 4) {
    return this->calcOptimalSplit_hist(Y, X, possible_splits, curr_family);
  } else {
    // Not used for res_trees; a separate path below uses its own pool type
    return Split{};
  }
}

void RandomPlantedForest::set_data(const NumericMatrix &samples_Y, const NumericMatrix &samples_X)
{
  this->Y = to_std_vec(samples_Y);
  this->X = to_std_vec(samples_X);
  if (Y.empty()) throw std::invalid_argument("Y empty - no data provided.");
  if (X.empty()) throw std::invalid_argument("X empty - no data provided.");
  this->feature_size = X[0].size();
  this->value_size = Y[0].size();
  for (const auto &vec : X) if (vec.size() != (size_t)feature_size) throw std::invalid_argument("Feature dimensions of X not uniform.");
  if (Y.size() != X.size()) throw std::invalid_argument("X and Y are not of the same length!");
  this->n_leaves = std::vector<int>(feature_size, 1);
  this->sample_size = X.size();
  this->upper_bounds = std::vector<double>(feature_size);
  this->lower_bounds = std::vector<double>(feature_size);
  for (int i = 0; i < feature_size; ++i) {
    double minVal = X[0][i], maxVal = X[0][i];
    for (size_t j = 0; j < sample_size; ++j) { double currVal = X[j][i]; if (currVal < minVal) minVal = currVal; if (currVal > maxVal) maxVal = currVal; }
    this->upper_bounds[i] = maxVal + 2 * eps; this->lower_bounds[i] = minVal;
  }
  // Prepare histogram bins if histogram mode is requested
  if (this->split_structure_mode_ == 4) {
    const size_t K = std::max<size_t>(2, std::min<size_t>(num_bins_, static_cast<size_t>(std::max(2, (int)std::sqrt((double)sample_size)))));
    this->num_bins_ = K;
    feature_cut_points_.assign((size_t)feature_size, std::vector<double>());
    sample_bin_id_.assign((size_t)feature_size, std::vector<int>(sample_size, 0));
    // For each feature, compute quantile cuts using sorted sample values
    for (int k = 0; k < feature_size; ++k) {
      std::vector<double> vals(sample_size);
      for (size_t i = 0; i < sample_size; ++i) vals[i] = X[i][k];
      std::sort(vals.begin(), vals.end());
      vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
      size_t unique_n = vals.size();
      size_t cuts = (K >= 2) ? (K - 1) : 1;
      if (unique_n <= 1 || cuts == 0) { feature_cut_points_[k].clear(); feature_cut_points_[k].shrink_to_fit(); continue; }
      feature_cut_points_[k].resize(cuts);
      for (size_t c = 1; c <= cuts; ++c) {
        double q = (double)c / (double)K; size_t idx = static_cast<size_t>(std::floor(q * (double)(unique_n - 1)));
        if (idx >= unique_n) idx = unique_n - 1; feature_cut_points_[k][c - 1] = vals[idx];
      }
      // Assign bin ids for all samples in original X for this feature
      for (size_t i = 0; i < sample_size; ++i) {
        double v = X[i][k];
        auto &cuts_k = feature_cut_points_[k];
        int bin = 0;
        if (!cuts_k.empty()) {
          auto itb = std::upper_bound(cuts_k.begin(), cuts_k.end(), v);
          bin = (int)std::distance(cuts_k.begin(), itb);
        }
        sample_bin_id_[k][i] = bin;
      }
    }
  }
  this->fit();
  if (cross_validate) { this->cross_validation(); }
}

void RandomPlantedForest::create_tree_family(std::vector<Leaf> initial_leaves, size_t n)
{
  TreeFamily curr_family;
  curr_family.insert({std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves))});

  // res_trees uses a separate pool
  if (split_structure_mode_ == 0) {
    std::vector<ResultingTreeCandidate> possible_trees;
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_trees.emplace_back(treePtr);
    }

    // Bootstrap samples
    int sample_index; std::vector<std::vector<double>> samples_X, samples_Y;
    if (deterministic) { samples_X = X; samples_Y = Y; this->t_try = 1; }
    else {
      samples_X = std::vector<std::vector<double>>(sample_size); samples_Y = std::vector<std::vector<double>>(sample_size);
      for (size_t i = 0; i < sample_size; ++i) { sample_index = rng_randint(0, (int)sample_size); samples_Y[i] = Y[sample_index]; samples_X[i] = X[sample_index]; }
    }

    Split curr_split;
    for (int split_count = 0; split_count < n_splits; ++split_count) {
      curr_split = this->calcOptimalSplit_resTrees(samples_Y, samples_X, possible_trees, curr_family);
      if (!std::isinf(curr_split.min_sum)) {
        // ensure D' and its one-step supersets are in pool
        std::set<int> Dprime = curr_split.tree_index->split_dims; Dprime.insert(curr_split.split_coordinate); Dprime.erase(0);
        if (!resultingTreeExists(possible_trees, Dprime)) { if (auto found = treeExists(Dprime, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({Dprime, std::make_shared<DecisionTree>(DecisionTree(Dprime))}); possible_trees.emplace_back(curr_family[Dprime]); } }
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> U = Dprime; U.insert(feature_dim); if (U.size() == Dprime.size()) continue; if (max_interaction >= 0 && U.size() > (size_t)max_interaction) continue; if (resultingTreeExists(possible_trees, U)) continue; if (auto found = treeExists(U, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({U, std::make_shared<DecisionTree>(DecisionTree(U))}); possible_trees.emplace_back(curr_family[U]); }
        }

        // Mutate residuals (restore old behavior)
        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            samples_Y[individual] -= curr_split.M_s;
          else
            samples_Y[individual] -= curr_split.M_b;
        }
        Leaf leaf_s, leaf_b; leaf_s.individuals = curr_split.I_s; leaf_b.individuals = curr_split.I_b; leaf_s.value = curr_split.M_s; leaf_b.value = curr_split.M_b; leaf_s.intervals = curr_split.leaf_index->intervals; leaf_b.intervals = curr_split.leaf_index->intervals; leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point; leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
        std::set<int> resulting_dims = curr_split.tree_index->split_dims; resulting_dims.insert(curr_split.split_coordinate); resulting_dims.erase(0);
        std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);
        if (!found_tree) {
          curr_family.insert({resulting_dims, std::make_shared<DecisionTree>(DecisionTree(resulting_dims))});
          found_tree = curr_family[resulting_dims];
        }
        if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate)) && delete_leaves) { leaf_s.value += curr_split.leaf_index->value; leaf_b.value += curr_split.leaf_index->value; *curr_split.leaf_index = leaf_b; curr_split.tree_index->leaves.push_back(leaf_s); }
        else { found_tree->leaves.push_back(leaf_s); found_tree->leaves.push_back(leaf_b); }
      }
    }

    auto keys = getKeys(curr_family);
    for (auto &key : keys) {
      if (curr_family[key]->leaves.size() == 0) { curr_family.erase(key); continue; }
      for (auto &leaf : curr_family[key]->leaves) leaf.individuals.clear();
    }
    tree_families[n] = curr_family; return;
  }

  // Non-res_trees modes use SplitCandidate pool
  std::vector<SplitCandidate> possible_splits;
  if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
    // leaves: seed with leaf-level candidates from null tree (single leaf at index 0)
    auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
        std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
        if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
        if (!leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
      }
    };
    auto null_tree = curr_family[{0}];
    if (!null_tree->leaves.empty()) add_leaf_candidates(null_tree, 0);

    // bootstrap
    int sample_index; std::vector<std::vector<double>> samples_X, samples_Y; std::vector<int> boot_idx(sample_size);
    if (deterministic) { samples_X = X; samples_Y = Y; this->t_try = 1; for (size_t i=0;i<sample_size;++i) boot_idx[i] = static_cast<int>(i); }
    else {
      samples_X = std::vector<std::vector<double>>(sample_size); samples_Y = std::vector<std::vector<double>>(sample_size);
      for (size_t i=0;i<sample_size;++i){ sample_index = rng_randint(0, (int)sample_size); boot_idx[i] = sample_index; samples_Y[i] = Y[sample_index]; samples_X[i] = X[sample_index]; }
    }

    // In histogram mode, cache per-feature bin ids for the working (bootstrapped) dataset
    if (split_structure_mode_ == 4) {
      tls_working_bin_id.assign((size_t)feature_size, std::vector<int>(sample_size, 0));
      for (int k = 0; k < feature_size; ++k) {
        // Reuse global precomputed bin ids via bootstrap index mapping
        if (!feature_cut_points_.empty() && (size_t)k < sample_bin_id_.size()) {
          for (size_t i = 0; i < sample_size; ++i) tls_working_bin_id[k][i] = sample_bin_id_[k][(size_t)boot_idx[i]];
        } else {
          // Fallback: compute on-the-fly (should be rare if cuts are available)
          const auto &cuts_k = (k >= 0 && k < (int)feature_cut_points_.size()) ? feature_cut_points_[k] : std::vector<double>{};
          for (size_t i = 0; i < sample_size; ++i) {
            int bin = 0; if (!cuts_k.empty()) { auto itb = std::upper_bound(cuts_k.begin(), cuts_k.end(), samples_X[i][k]); bin = (int)std::distance(cuts_k.begin(), itb); }
            tls_working_bin_id[k][i] = bin;
          }
        }
      }
    }

    Split curr_split;
    for (int split_count = 0; split_count < n_splits; ++split_count) {
      if (split_structure_mode_ == 4) curr_split = this->calcOptimalSplit_hist(samples_Y, samples_X, possible_splits, curr_family);
      else curr_split = this->calcOptimalSplit_leaves(samples_Y, samples_X, possible_splits, curr_family);
      if (!std::isinf(curr_split.min_sum)) {
        // Mutate residuals (restore old behavior)
        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) samples_Y[individual] -= curr_split.M_s; else samples_Y[individual] -= curr_split.M_b;
        }
        Leaf leaf_s, leaf_b; leaf_s.individuals = curr_split.I_s; leaf_b.individuals = curr_split.I_b; leaf_s.value = curr_split.M_s; leaf_b.value = curr_split.M_b; leaf_s.intervals = curr_split.leaf_index->intervals; leaf_b.intervals = curr_split.leaf_index->intervals; leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point; leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
        std::set<int> resulting_dims = curr_split.tree_index->split_dims; resulting_dims.insert(curr_split.split_coordinate); resulting_dims.erase(0);
        std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);
        if (!found_tree) {
          curr_family.insert({resulting_dims, std::make_shared<DecisionTree>(DecisionTree(resulting_dims))});
          found_tree = curr_family[resulting_dims];
        }
        auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
            if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
            if (!leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
          }
        };
        if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate)) && delete_leaves) {
          leaf_s.value += curr_split.leaf_index->value; leaf_b.value += curr_split.leaf_index->value;
          size_t idx_b = static_cast<size_t>(curr_split.leaf_index - &curr_split.tree_index->leaves[0]);
          *curr_split.leaf_index = leaf_b; curr_split.tree_index->leaves.push_back(leaf_s); size_t idx_s = curr_split.tree_index->leaves.size() - 1;
          add_leaf_candidates(curr_split.tree_index, idx_b);
          add_leaf_candidates(curr_split.tree_index, idx_s);
        } else {
          found_tree->leaves.push_back(leaf_s); found_tree->leaves.push_back(leaf_b); size_t idx_s = found_tree->leaves.size() - 2; size_t idx_b = found_tree->leaves.size() - 1;
          add_leaf_candidates(found_tree, idx_s); add_leaf_candidates(found_tree, idx_b);
        }
      }
    }
    auto keys = getKeys(curr_family);
    for (auto &key : keys) { if (curr_family[key]->leaves.size() == 0) { curr_family.erase(key); continue; } for (auto &leaf : curr_family[key]->leaves) leaf.individuals.clear(); }
    tree_families[n] = curr_family; return;
  }

  // cur_trees_1 and cur_trees_2: initialize with {j} trees
  for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
    auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
    curr_family.insert({{feature_dim}, treePtr});
    // leaf_idx unused for these modes
    possible_splits.emplace_back(feature_dim, treePtr, static_cast<size_t>(0));
  }

  // bootstrap
  int sample_index; std::vector<std::vector<double>> samples_X, samples_Y;
  if (deterministic) { samples_X = X; samples_Y = Y; this->t_try = 1; }
  else { samples_X = std::vector<std::vector<double>>(sample_size); samples_Y = std::vector<std::vector<double>>(sample_size); for (size_t i=0;i<sample_size;++i){ sample_index = rng_randint(0, (int)sample_size); samples_Y[i] = Y[sample_index]; samples_X[i] = X[sample_index]; } }

  Split curr_split;
  for (int split_count = 0; split_count < n_splits; ++split_count) {
    if (split_structure_mode_ == 2) curr_split = this->calcOptimalSplit_curTrees1(samples_Y, samples_X, possible_splits, curr_family);
    else curr_split = this->calcOptimalSplit_curTrees2(samples_Y, samples_X, possible_splits, curr_family);
    if (!std::isinf(curr_split.min_sum)) {
      // Update possible_splits like tryeveryleaf/splittrynew
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
        std::set<int> curr_dims = curr_split.tree_index->split_dims; curr_dims.insert(curr_split.split_coordinate); curr_dims.insert(feature_dim); curr_dims.erase(0);
        if (possibleExists(feature_dim, possible_splits, curr_dims)) continue;
        if (max_interaction >= 0 && curr_dims.size() > (size_t)max_interaction) continue;
        if (auto found = treeExists(curr_dims, curr_family)) possible_splits.emplace_back(feature_dim, found, static_cast<size_t>(0));
        else { curr_family.insert({curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))}); possible_splits.emplace_back(feature_dim, curr_family[curr_dims], static_cast<size_t>(0)); }
      }

      for (int individual : curr_split.leaf_index->individuals) {
        if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) samples_Y[individual] -= curr_split.M_s; else samples_Y[individual] -= curr_split.M_b;
      }
      Leaf leaf_s, leaf_b; leaf_s.individuals = curr_split.I_s; leaf_b.individuals = curr_split.I_b; leaf_s.value = curr_split.M_s; leaf_b.value = curr_split.M_b; leaf_s.intervals = curr_split.leaf_index->intervals; leaf_b.intervals = curr_split.leaf_index->intervals; leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point; leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
      std::set<int> resulting_dims = curr_split.tree_index->split_dims; resulting_dims.insert(curr_split.split_coordinate); resulting_dims.erase(0);
      std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);
      if (!found_tree) { curr_family.insert({resulting_dims, std::make_shared<DecisionTree>(DecisionTree(resulting_dims))}); found_tree = curr_family[resulting_dims]; }
      if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate)) && delete_leaves) { leaf_s.value += curr_split.leaf_index->value; leaf_b.value += curr_split.leaf_index->value; *curr_split.leaf_index = leaf_b; curr_split.tree_index->leaves.push_back(leaf_s); }
      else { found_tree->leaves.push_back(leaf_s); found_tree->leaves.push_back(leaf_b); }
    }
  }

  auto keys = getKeys(curr_family);
  for (auto &key : keys) { if (curr_family[key]->leaves.size() == 0) { curr_family.erase(key); continue; } for (auto &leaf : curr_family[key]->leaves) leaf.individuals.clear(); }
  tree_families[n] = curr_family;
}

// fit forest to new data
void RandomPlantedForest::fit()
{
  // setup initial set of individuals
  std::vector<int> initial_individuals(sample_size);
  std::iota(initial_individuals.begin(), initial_individuals.end(), 0);

  // initialize intervals with lower and upper bounds
  std::vector<Interval> initial_intervals(feature_size);
  for (int i = 0; i < feature_size; ++i)
    initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};

  // set properties of first leaf
  Leaf initial_leaf;
  {
    initial_leaf.value = std::vector<double>(value_size, 0);
    initial_leaf.individuals = initial_individuals;
    initial_leaf.intervals = initial_intervals;
    
  }
  std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf

  // initialize tree families
  this->tree_families = std::vector<TreeFamily>(n_trees);

  unsigned int threads_to_use = static_cast<unsigned int>(nthreads);
  // Prepare per-tree seeds derived from R's RNG to ensure reproducibility under threading (main thread only)
  {
    Rcpp::RNGScope scope; // ensure safe access to R RNG on main thread
    tree_seeds_.assign((size_t)std::max(1, n_trees), 0ULL);
    for (int i = 0; i < std::max(1, n_trees); ++i) {
      uint32_t a = (uint32_t)std::floor(R::runif(0.0, 4294967296.0));
      uint32_t b = (uint32_t)std::floor(R::runif(0.0, 4294967296.0));
      unsigned long long seed = ((unsigned long long)a << 32) ^ (unsigned long long)b;
      if (seed == 0ULL) seed = 88172645463393265ULL;
      tree_seeds_[(size_t)i] = seed;
    }
  }
  // We support multithreading in non-deterministic mode by using per-tree thread-local RNG seeded from R on the main thread.

  if (threads_to_use > 1)
  {
    if (threads_to_use > std::thread::hardware_concurrency())
    {
      Rcout << "Requested " << threads_to_use << " threads but only " << std::thread::hardware_concurrency() << " available" << std::endl;
    }
    for (int start = 0; start < n_trees; start += (int)threads_to_use)
    {
      int batch = std::min<int>((int)threads_to_use, n_trees - start);
      if (batch <= 0) break;
      std::vector<std::thread> threads((size_t)batch);
      for (int i = 0; i < batch; ++i)
      {
        int tree_index = start + i;
        threads[(size_t)i] = std::thread([this, &initial_leaves](int tree_index_inner){
          // Initialize thread-local RNG from pre-generated per-tree seed
          std::mt19937_64 rng_local;
          std::mt19937_64* prev_ptr = tls_rng_ptr;
          if (!tree_seeds_.empty() && (size_t)tree_index_inner < tree_seeds_.size()) {
            rng_local.seed(tree_seeds_[(size_t)tree_index_inner]);
          } else {
            rng_local.seed(88172645463393265ULL ^ (unsigned long long)tree_index_inner);
          }
          tls_rng_ptr = &rng_local;
          this->create_tree_family(initial_leaves, (size_t)tree_index_inner);
          tls_rng_ptr = prev_ptr; // restore
        }, tree_index);
      }
      for (auto &th : threads)
      {
        if (th.joinable()) th.join();
      }
    }
  }
  else
  {
    for (int n = 0; n < n_trees; ++n)
    {
      create_tree_family(initial_leaves, n);
    }
  }

  if (purify_forest)
  {
    this->purify_3();
  }
  else
  {
    purified = false;
  }
}

// predict single feature vector (from leaves variant)
std::vector<double> RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index)
{
  std::vector<double> total_res = std::vector<double>(value_size, 0);
  if (!purified)
  {
    if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
        for (auto &tree : tree_family)
          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (auto &dim : tree.first)
            {
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              {
                valid = false;
                break;
              }
            }
            if (valid) total_res += leaf.value;
          }
    }
    else
    {
      for (auto &tree_family : this->tree_families)
        for (auto &tree : tree_family)
        {
          if (tree.first != component_index) continue;
          std::vector<int> dims; for (auto dim : tree.first) dims.push_back(dim);
          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (unsigned int i = 0; i < dims.size(); ++i)
            {
              int dim = dims[i];
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[i] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[i] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              { valid = false; break; }
            }
            if (valid) total_res += leaf.value;
          }
        }
    }
  }
  else
  {
    if (component_index == std::set<int>{-1})
    {
      for (auto &tree_family : this->tree_families)
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0})
          { leaf_index = std::vector<int>(tree.first.size(), 0); total_res += tree.second->GridLeaves.values[leaf_index]; }
        }
    }
    else if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0}) { leaf_index = std::vector<int>(tree.first.size(), 0); }
          else {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              int dim = 0; { auto dim_pnt = tree.first.begin(); std::advance(dim_pnt, dim_index); dim = *dim_pnt; --dim; }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2) { leaf_index[dim_index] = 0; continue; }
              auto it = std::upper_bound(bounds.begin(), bounds.end(), X[dim]);
              int c = static_cast<int>(std::distance(bounds.begin(), it));
              leaf_index[dim_index] = std::min(std::max(0, c - 1), (int)bounds.size() - 2);
            }
          }
          for (int &index : leaf_index) index = std::max(0, index);
          total_res += tree.second->GridLeaves.values[leaf_index];
        }
    }
    else
    {
      for (auto &tree_family : this->tree_families)
        for (auto &tree : tree_family)
        {
          if (tree.first != component_index) continue;
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0}) { leaf_index = std::vector<int>(tree.first.size(), 0); }
          else {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              int dim = 0; { auto dim_pnt = tree.first.begin(); std::advance(dim_pnt, dim_index); dim = *dim_pnt; --dim; }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2) { leaf_index[dim_index] = 0; continue; }
              auto it = std::upper_bound(bounds.begin(), bounds.end(), X[dim_index]);
              int c = static_cast<int>(std::distance(bounds.begin(), it));
              leaf_index[dim_index] = std::min(std::max(0, c - 1), (int)bounds.size() - 2);
            }
          }
          for (int &index : leaf_index) index = std::max(0, index);
          total_res += tree.second->GridLeaves.values[leaf_index];
        }
    }
  }
  return total_res / n_trees;
}

Rcpp::NumericMatrix RandomPlantedForest::predict_matrix(const NumericMatrix &X, const NumericVector components)
{
  std::vector<std::vector<double>> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions;
  if (feature_vec.empty()) throw std::invalid_argument("Feature vector is empty.");
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec[0].size() != (size_t)this->feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  if (component_index != std::set<int>{0} && component_index != std::set<int>{-1} && component_index.size() != feature_vec[0].size())
    throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");
  for (auto &vec : feature_vec) predictions.push_back(predict_single(vec, component_index));
  return from_std_vec(predictions);
}

Rcpp::NumericMatrix RandomPlantedForest::predict_vector(const NumericVector &X, const NumericVector components)
{
  std::vector<double> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions; Rcpp::NumericMatrix res;
  if (feature_vec.empty()) { Rcout << "Feature vector is empty." << std::endl; return res; }
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec.size() != (size_t)this->feature_size) { Rcout << "Feature vector has wrong dimension." << std::endl; return res; }
  if (component_index == std::set<int>{0}) { predictions.push_back(predict_single(feature_vec, component_index)); }
  else { for (auto vec : feature_vec) predictions.push_back(predict_single(std::vector<double>{vec}, component_index)); }
  res = from_std_vec(predictions); return res;
}

double RandomPlantedForest::MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true)
{ return sum(Rcpp::pow(Y_true - Y_predicted, 2)) / Y_true.size(); }

double RandomPlantedForest::MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true)
{
  double sumv = 0; int Y_size = Y_predicted.size();
  for (int i = 0; i < Y_size; ++i) sumv += MSE_vec(Y_predicted(i, _), Y_true(i, _));
  return sumv / Y_size;
}

void RandomPlantedForest::print()
{
  for (int n = 0; n < n_trees; ++n)
  {
    TreeFamily family = tree_families[n]; auto keys = getKeys(family);
    for (size_t m = 0; m < keys.size(); ++m)
    {
      DecisionTree tree = *(family[keys[m]]);
      Rcout << m + 1 << " Tree: "; Rcout << "Dims="; for (const auto &dim : tree.split_dims) Rcout << dim << ",";
      Rcout << std::endl << "Leaves: (" << tree.leaves.size() << ")" << std::endl;
      for (const auto &leaf : tree.leaves)
      {
        Rcout << "Intervals="; for (const auto &interval : leaf.intervals) { Rcout << interval.first << "," << interval.second << "/"; }
        Rcout << " Value="; for (const auto &val : leaf.value) Rcout << val << ", "; Rcout << std::endl;
      }
      Rcout << std::endl;
    }
    Rcout << std::endl << std::endl;
  }
}

void RandomPlantedForest::get_parameters()
{
  Rcout << "Parameters: n_trees=" << n_trees << ", n_splits=" << n_splits << ", max_interaction=" << max_interaction << ", t_try=" << t_try
        << ", split_decay_rate=" << split_decay_rate_<< ", max_candidates="  << max_candidates_
        << ", split_try=" << split_try << ", purified=" << purified << ", deterministic=" << deterministic << ", nthreads=" << nthreads
        << ", feature_size=" << feature_size << ", sample_size=" << sample_size
        << ", split_structure_mode=" << split_structure_mode_ << std::endl;
}

void RandomPlantedForest::set_parameters(StringVector keys, NumericVector values)
{
  if (keys.size() != values.size()) { Rcout << "Size of input vectors is not the same. " << std::endl; return; }
  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    if (keys[i] == "deterministic") this->deterministic = values[i];
    else if (keys[i] == "nthreads") this->nthreads = values[i];
    else if (keys[i] == "purify") this->purify_forest = values[i];
    else if (keys[i] == "n_trees") this->n_trees = values[i];
    else if (keys[i] == "n_splits") this->n_splits = values[i];
    else if (keys[i] == "t_try") this->t_try = values[i];
    else if (keys[i] == "split_try") this->split_try = values[i];
    else if (keys[i] == "max_interaction") this->max_interaction = values[i];
    else if (keys[i] == "cv") this->cross_validate = values[i];
    else if (keys[i] == "split_decay_rate") this->split_decay_rate_ = values[i];
    else if (keys[i] == "max_candidates") this->max_candidates_ = static_cast<size_t>(values[i]);
    else if (keys[i] == "delete_leaves") this->delete_leaves = static_cast<bool>(values[i]);
    else if (keys[i] == "leaf_feature_cache_cap") this->leaf_feature_cache_cap_ = static_cast<size_t>(values[i]);
    
    else if (keys[i] == "split_structure_mode") this->split_structure_mode_ = static_cast<int>(values[i]);
    else Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
  }
  this->fit();
}

List RandomPlantedForest::get_model()
{
  List model;
  for (const auto &family : tree_families)
  {
    List variables, family_values, family_intervals;
    for (const auto &tree : family)
    {
      List tree_values; List tree_intervals; variables.push_back(from_std_set(tree.first));
      for (const auto &leaf : tree.second->leaves)
      {
        NumericMatrix leaf_values; for (const auto &val : leaf.value) leaf_values.push_back(val);
        tree_values.push_back(leaf_values);
        NumericVector intervals; for (const auto &interval : leaf.intervals) { intervals.push_back(interval.first); intervals.push_back(interval.second); }
        NumericMatrix leaf_intervals(2, feature_size, intervals.begin()); tree_intervals.push_back(leaf_intervals);
      }
      family_intervals.push_back(tree_intervals); family_values.push_back(tree_values);
    }
    model.push_back(List::create(Named("variables") = variables, _["values"] = family_values, _["intervals"] = family_intervals));
  }
  return (model);
}



void RandomPlantedForest::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries)
{

  /*
   bool cv_tmp = this->cross_validate;
   this->cross_validate = false;
   if(deterministic) {
   Rcout << "Note: Set model to non-deterministic. " << std::endl;
   deterministic = false;
   }
   std::set<int> splits_vec = to_std_set(splits);
   std::vector<int> split_tries_vec = to_std_vec(split_tries);
   std::vector<double> t_tries_vec = to_std_vec(t_tries);
   if(splits_vec.size()!=2) {Rcout << "Min and max needed for number of splits." << std::endl; return;}
   // remember optimal parameter set and MSE
   double  MSE_sum = 0, curr_MSE = 0, MSE_min = INF, optimal_split = INF, optimal_t_try = INF, optimal_split_try = INF;
   int optimal_inter = 1;
   std::vector<int> order(sample_size);
   std::iota(order.begin(), order.end(), 0);
   std::random_shuffle(order.begin(), order.end(), randWrapper);
   double tmp = double(sample_size)/double(n_sets);
   int set_size = round(tmp);
   // remember original data samples
   NumericMatrix X_original = from_std_vec(X);
   NumericVector Y_original = from_std_vec(Y);
   // set level of interactions
   std::set<int> interactions{1};
   if(feature_size >= 2){
   interactions.insert(2);
   interactions.insert(feature_size);
   }
   // go through all parameter combinations
   for(int inter: interactions){
   this->max_interaction = inter;
   for(int splits=*splits_vec.begin(); splits<=*--splits_vec.end(); splits=ceil(splits*1.2)){
   this->n_splits = splits;
   for(auto t: t_tries){
   this->t_try = t;
   for(auto s: split_tries){
   this->split_try = s;
   // k-fold cross-validation: go over all possible combinations as test set
   MSE_sum = 0;
   for(int n_set=0; n_set<n_sets; ++n_set){
   // split data into training and test sets
   int test_size = set_size;
   if(n_set == n_sets-1) test_size = order.size() - (n_sets-1) * set_size;
   int train_size = order.size() - test_size, i = 0, j = 0;
   NumericVector Y_train(train_size), Y_test_true(test_size), Y_test_predicted;
   NumericMatrix X_train(train_size, feature_size), X_test(test_size, feature_size);
   for(int index=0; index<order.size(); ++index){
   if( (index >= (n_set * set_size)) && (index < ((n_set + 1) * set_size))){
   Y_test_true[i] = Y_original[order[index]];
   X_test(i, _ ) = X_original(order[index], _ );
   ++i;
   }else{
   Y_train[j] = Y_original[order[index]];
   X_train(j, _ ) = X_original(order[index], _ );
   ++j;
   }
   }
   // fit to training data
   this->set_data(Y_train, X_train);
   // predict with test set and determine mse
   Y_test_predicted = this->predict_matrix(X_test);
   MSE_sum += this->MSE(Y_test_predicted, Y_test_true);
   }
   // average
   curr_MSE = MSE_sum / n_sets;
   Rcout << inter << ", " << splits << ", " << t << ", " << s << ": MSE=" << curr_MSE << std::endl;
   // update optimal
   if(curr_MSE < MSE_min){
   MSE_min = curr_MSE;
   optimal_split = splits;
   optimal_t_try = t;
   optimal_split_try = s;
   optimal_inter = inter;
   }
   }
   }
   }
   }
   // reset X&Y to original and fit with optimal pars
   this->n_splits = optimal_split;
   this->t_try = optimal_t_try;
   this->split_try = optimal_split_try;
   this->max_interaction = optimal_inter;
   this->set_data(Y_original, X_original);
   this->cross_validate = cv_tmp;
   Rcout << "Optimal parameters: " << optimal_inter << ", " << optimal_split << ", " << optimal_t_try << ", " << optimal_split_try << ": MSE=" << MSE_min << std::endl;
   */
}



void RandomPlantedForest::purify_1()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (auto tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    while (curr_max >= 1)
    {

      // go through split dimensions of all trees
      auto keys = getKeys(curr_family);
      std::vector<std::set<int>>::reverse_iterator key = keys.rbegin();
      while (key != keys.rend())
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (curr_max == 1)
              {
                tree = curr_family[std::set<int>{0}];
              }
              else
              {
                if (!tree)
                {
                  curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                  tree = curr_family[tree_dims];
                }
              }

              // go through leaves of current tree
              int n_leaves = curr_tree->leaves.size();
              for (int l = 0; l < n_leaves; ++l)
              {
                auto &curr_leaf = curr_tree->leaves[l];

                double multiplier = (curr_leaf.intervals[feature_dim - 1].second - curr_leaf.intervals[feature_dim - 1].first) / (upper_bounds[feature_dim - 1] - lower_bounds[feature_dim - 1]);

                // new leaf including intervals and value
                Leaf new_leaf = curr_leaf; // initialize intervals with first leaf
                new_leaf.intervals[feature_dim - 1].first = lower_bounds[feature_dim - 1];
                new_leaf.intervals[feature_dim - 1].second = upper_bounds[feature_dim - 1];
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = -curr_leaf.value[i] * multiplier; // update value of new leaf

                // append new leaf
                if (!leafExists(new_leaf.intervals, curr_tree))
                  curr_tree->leaves.push_back(new_leaf);
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = curr_leaf.value[i] * multiplier; // update value of new leaf
                if (!leafExists(new_leaf.intervals, tree))
                  tree->leaves.push_back(new_leaf);
              }
            }
          }
        }
        key++;
      }

      // update currently considered dimension size
      --curr_max;
    }
  }

  purified = true;
}

void RandomPlantedForest::purify_2()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);

    // go through all variables of the component
    for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
    {
      std::vector<double> bounds;

      // go through trees of family
      for (const auto &curr_tree : curr_family)
      {

        // consider only relevant trees that have current dimension as variable
        if (!curr_tree.first.count(curr_dim))
          continue;

        // go through leaves of tree
        for (const auto &curr_leaf : curr_tree.second->leaves)
        {
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
        }
      }
      std::sort(bounds.begin(), bounds.end());
      bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
      lim_list[curr_dim - 1] = bounds;
    }

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids(curr_family.size() - 1);
    std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
    std::vector<std::set<int>> variables(curr_family.size() - 1);

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {

      if (curr_tree.first == std::set<int>{0})
        continue; // ignore null tree

      // fill space with dimensions
      std::vector<int> dimensions;
      for (const auto &dim : curr_tree.first)
      {
        dimensions.push_back(lim_list[dim - 1].size() - 1); // size - 1 ?
      }

      // setup grid for leaf indices
      auto grid = grid::NDGrid(dimensions);

      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
      values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
      variables[tree_index] = curr_tree.first;

      // fill grid points with individuals and values
      while (!grid.nextPoint())
      {

        std::vector<int> gridPoint = grid.getPoint();

        bool in_leaf = true;

        // go through sample points to sum up individuals
        for (const auto &feature_vec : X)
        {
          int dim_index = 0;
          in_leaf = true;
          for (const auto &dim : curr_tree.first)
          {
            double val = feature_vec[dim - 1];
            if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // consider individuals only if all in
          if (in_leaf)
            individuals[tree_index][gridPoint] += 1;
        }

        // go through leaves of tree to sum up values
        for (const auto &leaf : curr_tree.second->get_leaves())
        {

          in_leaf = true;
          int dim_index = 0;
          for (const auto &dim : curr_tree.first)
          {
            // consider values only if all in
            if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // sum up values
          if (in_leaf)
            values[tree_index][gridPoint] += leaf.value; // todo: multiclass
        }
      }

      ++tree_index;
    }

    // ------------- create new trees -------------

    // insert null tree
    grids.insert(grids.begin(), grid::NDGrid());
    values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
    variables.insert(variables.begin(), std::set<int>{0});

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (const auto &tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    auto keys = getKeys(curr_family);
    while (curr_max > 1)
    {

      // go through split dimensions of all trees
      for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          int dim_index = 0;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (!tree)
              {

                // get index of old and new tree
                auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
                curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));

                // remove matrix dimension of respective variable
                std::vector<int> matrix_dimensions = values[old_tree_index].dims;
                matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);

                // initialize data for new tree
                auto grid = grid::NDGrid(matrix_dimensions);
                grids.insert(grids.begin() + tree_index, grid);
                // initialize with zero vector of length value_size
                values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
                individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
                variables.insert(variables.begin() + tree_index, tree_dims);

                // fill individuals of new trees
                while (!grid.nextPoint())
                {

                  std::vector<int> gridPoint = grid.getPoint();
                  bool in_leaf = true;

                  // go through sample points to sum up individuals
                  for (const auto &feature_vec : X)
                  {
                    int dim_index = 0;
                    in_leaf = true;
                    for (const auto &dim : tree_dims)
                    {
                      double val = feature_vec[dim - 1];
                      if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
                        in_leaf = false;
                      ++dim_index;
                    }

                    // consider individuals only if all in
                    if (in_leaf)
                      individuals[tree_index][gridPoint] += 1;
                  }
                }
              }

              dim_index++;
            }
          }
        }
      }

      // update currently considered dimension size
      --curr_max;
    }

    // ------------- purify -------------

    // measure tolerance and number of iterations
    std::vector<double> tol(curr_family.size(), 1);
    int iter;

    // iterate backwards through tree family
    int curr_tree_index = curr_family.size() - 1;
    for (TreeFamily::reverse_iterator curr_tree = curr_family.rbegin(); curr_tree != curr_family.rend(); ++curr_tree)
    {
      iter = 0;
      std::set<int> curr_dims = curr_tree->second->get_split_dims();

      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;

      // repeat until tolerance small enough and (?) maximum number of iterations reached
      while ((tol[curr_tree_index] > 0.00000000001) && (iter < 100))
      {

        // go through feature dims
        int curr_dim_index = 0;
        for (const auto &feature_dim : curr_dims)
        {

          // get tree that has same variables as curr_tree minus j-variable
          std::set<int> tree_dims = curr_dims;
          tree_dims.erase(tree_dims.find(feature_dim));
          int tree_index = 0; // if tree not exist, set to null tree
          if (curr_family.find(tree_dims) != curr_family.end())
            tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims)) - 1;

          // update values
          if (grids[curr_tree_index].dimensions.size() == 1)
          { // one dimensional case

            int sum_ind = 0;
            std::vector<double> avg(value_size, 0);

            // get sum of individuals
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              sum_ind += individuals[curr_tree_index][tmp];
            }
            if (sum_ind == 0)
              continue;

            // calc avg
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              avg += (individuals[curr_tree_index][tmp] * values[curr_tree_index][tmp]) / sum_ind;
            }

            // update values of one dimensional and null tree
            for (int i = 0; i < values[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              values[curr_tree_index][tmp] -= avg;
            }
            std::vector<int> tmp{0};
            values[tree_index][tmp] += avg;
          }
          else
          { // higher dimensional case

            // setup new grid without dimension j
            std::vector<int> new_dimensions = grids[curr_tree_index].dimensions;
            int j_dim = new_dimensions[curr_dim_index];
            new_dimensions.erase(new_dimensions.begin() + curr_dim_index);
            grid::NDGrid grid = grid::NDGrid(new_dimensions);

            // go through values without dimension j
            while (!grid.nextPoint())
            {
              auto gridPoint = grid.getPoint();
              gridPoint.push_back(0);

              int sum_ind = 0;
              std::vector<double> avg(value_size, 0);

              // go through slice to sum up individuals
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // get sum of individuals
                sum_ind += individuals[curr_tree_index][gridPoint];
              }

              // go through slice to calc avg
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // calc avg
                avg += (individuals[curr_tree_index][gridPoint] * values[curr_tree_index][gridPoint]) / sum_ind;
              }

              // go through slice to update values
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // update values of current slice
                values[curr_tree_index][gridPoint] -= avg;
              }

              // update lower dimensional tree
              gridPoint.pop_back();
              values[tree_index][gridPoint] += avg;
            }
          }

          ++curr_dim_index;
        }

        // update tolerance
        if (variables[curr_tree_index].size() == 1)
        {
          tol[curr_tree_index] = 1; // todo
        }
        else
        {
          tol[curr_tree_index] = 1;
        }

        ++iter;
      }

      --curr_tree_index;
    }

    // ------------- attach to rpf class -------------

    // fill with new trees
    for (size_t tree_index = 0; tree_index < variables.size(); ++tree_index)
    {
      LeafGrid curr_gridLeaf;
      curr_gridLeaf.grid = grids[tree_index];
      curr_gridLeaf.individuals = individuals[tree_index];
      curr_gridLeaf.lim_list = lim_list;
      curr_gridLeaf.values = values[tree_index];
      curr_family[variables[tree_index]]->GridLeaves = curr_gridLeaf;
    }
  }

  purified = true;
}

void RandomPlantedForest::purify_3()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);

    // go through all variables of the component
    for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
    {
      std::vector<double> bounds;

      // go through trees of family
      for (const auto &curr_tree : curr_family)
      {

        // consider only relevant trees that have current dimension as variable
        if (!curr_tree.first.count(curr_dim))
          continue;

        // go through leaves of tree
        for (const auto &curr_leaf : curr_tree.second->leaves)
        {
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].first);
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
        }
      }
      std::sort(bounds.begin(), bounds.end());
      bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
      // int i_last = bounds.size()-1;
      // double bibi = bounds[i_last] + 0.0001;
      // bounds[i_last] = bounds[i_last] + 0.0001;
      lim_list[curr_dim - 1] = bounds;
    }

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids(curr_family.size() - 1);
    std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values_old(curr_family.size() - 1);
    std::vector<std::set<int>> variables(curr_family.size() - 1);

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {

      if (curr_tree.first == std::set<int>{0})
      {

        // values[tree_index] = rpf::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
        continue; // ignore null tree
      }

      // fill space with dimensions
      std::vector<int> dimensions;
      for (const auto &dim : curr_tree.first)
      {
        // lim_list holds interval endpoints; number of cells is (#endpoints - 1)
        int num_bounds = static_cast<int>(lim_list[dim - 1].size());
        dimensions.push_back(std::max(1, num_bounds - 1));
      }

      // setup grid for leaf indices
      auto grid = grid::NDGrid(dimensions);

      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
      values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));     // changed
      values_old[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
      variables[tree_index] = curr_tree.first;

      // fill grid points with individuals and values
      while (!grid.nextPoint())
      {

        std::vector<int> gridPoint = grid.getPoint();

        bool in_leaf = true;

        // go through sample points to sum up individuals
        for (const auto &feature_vec : X)
        {
          int dim_index = 0;
          in_leaf = true;
          for (const auto &dim : curr_tree.first)
          {
            double val = feature_vec[dim - 1];
            if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // consider individuals only if all in
          if (in_leaf)
            individuals[tree_index][gridPoint] += 1;
        }

        // go through leaves of tree to sum up values
        for (const auto &leaf : curr_tree.second->get_leaves())
        {

          in_leaf = true;
          int dim_index = 0;
          for (const auto &dim : curr_tree.first)
          {
            // consider values only if all in
            if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // sum up values
          if (in_leaf)
          {

            values[tree_index][gridPoint] += leaf.value;     // todo: multiclass
            values_old[tree_index][gridPoint] += leaf.value; // todo: multiclass
          }
        }
      }

      ++tree_index;
    }

    // Rcout << variables.size();
    // for(int i = 0; i<variables.size(); ++i){
    //
    //   // Rcout << variables[i].size();
    //
    //   for(auto dim: variables[i]) Rcout << dim << ",";
    //
    //   //  Rcout << variables[i][j] << ",";
    //   //}
    //
    //   Rcout << std::endl;
    // }

    // ------------- create new trees -------------

    // insert null tree
    grids.insert(grids.begin(), grid::NDGrid());
    values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    values_old.insert(values_old.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
    variables.insert(variables.begin(), std::set<int>{0});

    // recap maximum number of dimensions of current family
    unsigned int curr_max = curr_family.rbegin()->first.size();

    while (curr_max > 1)
    {

      auto keys = getKeys(curr_family);
      // go through split dimensions of all trees
      for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
      {
        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;
        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {
          // go through feature dims
          int dim_index = 0;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {
            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {
              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree
              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (!tree)
              {
                // get index of old and new tree
                auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
                curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
                // remove matrix dimension of respective variable
                std::vector<int> matrix_dimensions = values[old_tree_index].dims;
                // std::vector<int> matrix_dimensions = values_old[old_tree_index].dims;

                // Rcout << typeof(matrix_dimensions.begin()) << std::endl;

                matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);
                // initialize data for new tree
                auto grid = grid::NDGrid(matrix_dimensions);
                grids.insert(grids.begin() + tree_index, grid);
                values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
                values_old.insert(values_old.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
                individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
                variables.insert(variables.begin() + tree_index, tree_dims);
                // fill individuals of new trees
                while (!grid.nextPoint())
                {
                  std::vector<int> gridPoint = grid.getPoint();
                  bool in_leaf = true;
                  // go through sample points to sum up individuals
                  for (const auto &feature_vec : X)
                  {
                    int dim_index2 = 0;
                    in_leaf = true;
                    for (const auto &dim : tree_dims)
                    {
                      double val = feature_vec[dim - 1];
                      if (!((val >= lim_list[dim - 1][gridPoint[dim_index2]]) && (val < lim_list[dim - 1][gridPoint[dim_index2] + 1])))
                        in_leaf = false;
                      ++dim_index2;
                    }
                    // consider individuals only if all in
                    if (in_leaf)
                      individuals[tree_index][gridPoint] += 1;
                  }
                }
              }
              dim_index++;
            }
          }
        }
      }
      // update currently considered dimension size
      --curr_max;
    }

    // Rcout << std::endl;
    // Rcout << std::endl;
    // Rcout << std::endl;
    //
    // for(int i = 0; i<variables.size(); ++i){
    //
    //   // Rcout << variables[i].size();
    //
    //   for(auto dim: variables[i]) Rcout << dim << ",";
    //
    //   //  Rcout << variables[i][j] << ",";
    //   //}
    //
    //   Rcout << std::endl;
    // }

    // ------------- purify -------------
    // iterate backwards through tree family
    int tree_index_t = curr_family.size() - 1;
    for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t)
    {
      std::set<int> curr_dims = *tree_t;
      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;
      // Rcout << std::endl << tree_index_t << " - T: ";
      //  Rcout << "tree_t:";
      //  for(auto dim: curr_dims) Rcout << dim << ", ";
      //  Rcout << std::endl;

      auto grid = grids[tree_index_t];
      //     Rcout << "Grid dimensions of T: ";
      //     for(auto dim: grid.dimensions) Rcout << dim << ", ";
      //     Rcout << std::endl;
      // go through subtrees of t
      int tree_index_u = variables.size();
      for (auto tree_u = variables.rbegin(); tree_u != variables.rend(); ++tree_u)
      {
        --tree_index_u;
        // j_dims = dims of t without u
        std::set<int> j_dims = curr_dims;
        if (tree_u->size() > curr_dims.size())
          continue;
        // check if subset
        bool subset = true;
        for (const auto dim : *tree_u)
        {
          if (tree_t->count(dim) == 0)
          {
            subset = false;
            break;
          }
          j_dims.erase(dim);
        }
        if (!subset)
          continue;

        // Rcout << "Hello";
        // Rcout << "   " << tree_index_u << " - U: ";
        // for(auto dim: *tree_u) Rcout << dim << ", ";
        // Rcout << std::endl;
        // Rcout << "   Individuals: ";

        double tot_sum = 0;
        grid = grids[tree_index_u];
        while (!grid.nextPoint())
        {
          auto gridPoint = grid.getPoint();
          //     Rcout << individuals[tree_index_u][gridPoint] << ", ";
          tot_sum += individuals[tree_index_u][gridPoint];
        }
        // Rcout << "Total sum: " << tot_sum << std::endl;
        // Rcout << std::endl;

        grid = grids[tree_index_u];
        //     Rcout << "      Grid dimensions of U: ";
        //     for(auto dim: grid.dimensions) Rcout << dim << ", ";
        //     Rcout << std::endl;

        // Rcout<< "j_dims: "<<j_dims.size() << std::endl;;

        std::vector<double> update(value_size, 0);

        if (j_dims.size() == 0)
        {

          // grid = grids[tree_index_u];
          while (!grid.nextPoint())
          {
            auto gridPoint_i = grid.getPoint();
            //     Rcout << "         " << "i: ";
            //     for(auto p: gridPoint_i) Rcout << p << ", ";
            //     Rcout << std::endl << "         ";
            double curr_sum = individuals[tree_index_u][gridPoint_i];
            //     Rcout << ", Current Sum: " << curr_sum << std::endl;
            //     Rcout << std::endl << "         " << "i, j: ";
            update += (curr_sum / tot_sum) * values_old[tree_index_t][gridPoint_i];
            //     Rcout << std::endl;
          }

          int tree_index_s = variables.size();
          for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
          {

            // Rcout << "tree_s:";
            // for(auto dim: *tree_s) Rcout << dim << ", ";
            // Rcout << std::endl;

            --tree_index_s;
            if (*tree_s == std::set<int>{0})
            {

              auto gridPoint_0 = std::vector<int>{0};
              values[tree_index_s][gridPoint_0] += update;
              //     Rcout << std::endl;
              //}

              /*
               for(auto tree_0: curr_family){

               if(tree_0.first == std::set<int>{0}){

               Rcout << tree_0.first.size();
               std::vector<int> leaf_index(tree_0.first.size(), 0);
               std::vector<int> leaf_index(tree_0.second->GridLeaves.values.size(), 0);

               int Test = tree_0.second->GridLeaves.values.size();
               Rcout << Test;
               tree_0.second->GridLeaves.values[leaf_index] += update;
               }
               }
               */
            }
            else
            {

              // check if S subset of T

              bool subset = true;
              for (const auto dim : *tree_s)
              {
                if (tree_t->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              if (!subset)
                continue;

              // Rcout << pow(-1, (*tree_s).size()) << std::endl;

              auto grid_k = grids[tree_index_s];
              while (!grid_k.nextPoint())
              {
                auto gridPoint_k = grid_k.getPoint();
                //
                //      if((*tree_s).size()>2){
                //      Rcout << std::endl << "            " << "j, k: ";
                //      for(auto p: gridPoint_k) Rcout << p << ", ";
                //      Rcout << std::endl;
                //      }
                //
                //      Rcout << pow(-1, (*tree_s).size()) * update << std::endl;
                values[tree_index_s][gridPoint_k] += pow(-1, (*tree_s).size()) * update;
              }
            }
          }
          // Rcout << std::endl;
        }
        else
        {

          std::vector<int> j_sizes(j_dims.size(), 0);
          for (size_t j = 0; j < j_dims.size(); ++j)
          {
            auto tmp = j_dims.begin();
            std::advance(tmp, j);
            int j_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*tmp));
            j_sizes[j] = grids[tree_index_t].dimensions[j_index];
          }

          // Rcout<<"Hello 1";

          grid::NDGrid grid_j = grid::NDGrid(j_sizes);
          while (!grid_j.nextPoint())
          {

            std::vector<double> update(value_size, 0);
            auto gridPoint_j = grid_j.getPoint();
            //     Rcout << "         " << "j: ";
            //     for(auto p: gridPoint_j) Rcout << p << ", ";
            //     Rcout << std::endl;
            // calc update
            grid = grids[tree_index_u];
            while (!grid.nextPoint())
            {
              auto gridPoint_i = grid.getPoint();
              //     Rcout << "         " << "i: ";
              //     for(auto p: gridPoint_i) Rcout << p << ", ";
              //     Rcout << std::endl << "         ";
              double curr_sum = individuals[tree_index_u][gridPoint_i];
              //     Rcout << ", Current Sum: " << curr_sum << std::endl;
              std::vector<int> gridPoint_ij(tree_t->size(), 0);
              for (size_t j = 0; j < gridPoint_j.size(); ++j)
              {
                auto j_dim = j_dims.begin();
                std::advance(j_dim, j);
                int j_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*j_dim));
                //     Rcout << "         j_dim=" << *j_dim << ", j_index=" << j_index;
                gridPoint_ij[j_index] = gridPoint_j[j];
              }
              for (size_t i = 0; i < gridPoint_i.size(); ++i)
              {
                auto i_dim = tree_u->begin();
                std::advance(i_dim, i);
                int i_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*i_dim));
                //     Rcout << "         i_dim=" << *i_dim << ", i_index=" << i_index;
                gridPoint_ij[i_index] = gridPoint_i[i];
              }
              //     Rcout << std::endl << "         " << "i, j: ";
              //     for(auto p: gridPoint_ij) Rcout << p << ", ";
              //     Rcout << std::endl;
              update += (curr_sum / tot_sum) * values_old[tree_index_t][gridPoint_ij];
              //     Rcout << std::endl;
            }

            // Rcout << "Hello_2";
            // update trees
            int tree_index_s = variables.size();
            for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
            {
              --tree_index_s;
              // check if T\U=j_dims subset of S and S subset of T
              bool subset = true;
              for (const auto dim : j_dims)
              {
                if (tree_s->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              for (const auto dim : *tree_s)
              {
                if (tree_t->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              if (!subset)
                continue;
              //     Rcout << "         " << "S: ";
              //     for(auto dim: *tree_s) Rcout << dim << ", ";
              //     Rcout << std::endl;
              // S cap U
              std::set<int> k_dims = *tree_s;
              std::set<int> k_dims_h1 = *tree_s;
              std::set<int> k_dims_h2 = *tree_u;
              for (const auto dim : *tree_u)
                k_dims.insert(dim);
              for (const auto dim : *tree_s)
                k_dims_h2.erase(dim);
              for (const auto dim : *tree_u)
                k_dims_h1.erase(dim);
              for (const auto dim : k_dims_h1)
                k_dims.erase(dim);
              for (const auto dim : k_dims_h2)
                k_dims.erase(dim);

              // std::set<int> k_dims = *tree_s;
              // for(const auto dim: *tree_t) k_dims.erase(dim);
              // for(const auto dim: *tree_u) k_dims.insert(dim);

              //     Rcout << "         " << "k_dims: ";
              //     for(auto dim: k_dims) Rcout << dim << ", ";
              //     Rcout << std::endl;

              if (k_dims.size() == 0)
              {

                values[tree_index_s][gridPoint_j] += pow(-1, (*tree_s).size() - j_dims.size()) * update;
              }
              else
              {

                // Rcout <<"k_dims :";
                // for(auto dim: k_dims) Rcout << dim << ", ";
                // Rcout << std::endl;

                std::vector<int> k_sizes(k_dims.size(), 0);
                for (size_t k = 0; k < k_dims.size(); ++k)
                {
                  auto tmp = k_dims.begin();
                  std::advance(tmp, k);
                  int k_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*tmp));
                  k_sizes[k] = grids[tree_index_t].dimensions[k_index];
                }
                // Rcout << "         " << "k_sizes: ";
                // for(auto dim: k_sizes) Rcout << dim << ", ";
                // Rcout << std::endl;
                grid::NDGrid grid_k = grid::NDGrid(k_sizes);
                while (!grid_k.nextPoint())
                {
                  auto gridPoint_k = grid_k.getPoint();
                  // Rcout << "            " << "k: ";
                  // for(auto p: gridPoint_k) Rcout << p << ", ";
                  // Rcout << std::endl << "         ";
                  std::vector<int> gridPoint_jk(tree_s->size(), 0);
                  for (size_t j = 0; j < gridPoint_j.size(); ++j)
                  {
                    auto j_dim = j_dims.begin();
                    std::advance(j_dim, j);
                    int j_index = std::distance(variables[tree_index_s].begin(), variables[tree_index_s].find(*j_dim));
                    // Rcout << "         j_dim=" << *j_dim << ", j_index=" << j_index;
                    gridPoint_jk[j_index] = gridPoint_j[j];
                  }
                  for (size_t k = 0; k < gridPoint_k.size(); ++k)
                  {
                    auto k_dim = k_dims.begin();
                    std::advance(k_dim, k);
                    int k_index = std::distance(variables[tree_index_s].begin(), variables[tree_index_s].find(*k_dim));
                    // Rcout << "         k_dim=" << *k_dim << ", k_index=" << k_index;
                    gridPoint_jk[k_index] = gridPoint_k[k];
                  }
                  // Rcout << std::endl << "            " << "j, k: ";
                  // for(auto p: gridPoint_jk) Rcout << p << ", ";
                  // Rcout << std::endl;

                  // Rcout << pow(-1, (*tree_s).size() - j_dims.size()) * update[0];
                  values[tree_index_s][gridPoint_jk] += pow(-1, (*tree_s).size() - j_dims.size()) * update;
                }
              }
            }
          }
        }
      }
      --tree_index_t;
    }

    // ------------- attach to rpf class -------------

    // fill with new trees
    for (size_t tree_index = 0; tree_index < variables.size(); ++tree_index)
    {
      LeafGrid curr_gridLeaf;
      curr_gridLeaf.grid = grids[tree_index];
      curr_gridLeaf.individuals = individuals[tree_index];
      curr_gridLeaf.lim_list = lim_list;
      curr_gridLeaf.values = values[tree_index];
      curr_family[variables[tree_index]]->GridLeaves = curr_gridLeaf;
    }
  }

  purified = true;
}
