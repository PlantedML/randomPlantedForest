// Split-mode: cur_trees_1. Samples feasible leaves proportionally to their
// number of candidate thresholds, then evaluates a single threshold.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

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
  } else { for (size_t i=0;i+n_candidates<=possible_splits.size() && i< n_candidates;++i) sample_idxs.push_back(i); }
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx); int k = it->dim - 1; int leaf_size = this->n_leaves[k];
    std::set<int> Dprime_minus_k = it->tree->split_dims; Dprime_minus_k.erase(k + 1); Dprime_minus_k.erase(0);
    std::vector<std::shared_ptr<DecisionTree>> sources; sources.reserve(2);
    if (Dprime_minus_k.empty()) { if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end()) sources.push_back(itZero->second); }
    else { if (auto itS = curr_family.find(Dprime_minus_k); itS != curr_family.end()) sources.push_back(itS->second); }
    if (auto itD = curr_family.find(it->tree->split_dims); itD != curr_family.end()) if (sources.empty() || sources.back().get() != itD->second.get()) sources.push_back(itD->second);

    if (!this->deterministic) {
      struct LeafEnv { std::shared_ptr<DecisionTree> tree; Leaf* leaf; int left; int right_env; };
      std::vector<LeafEnv> env; env.reserve(sources.size()*8);
      std::vector<double> w_env; w_env.reserve(sources.size()*8);
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty()) continue;
        for (auto &leaf : src_tree->leaves) {
          const int m = static_cast<int>(leaf.individuals.size());
          const int left = leaf_size; const int right_env = m - leaf_size; const int width = right_env - left;
          if (width > 0) { env.push_back({src_tree, &leaf, left, right_env}); w_env.push_back(static_cast<double>(width)); }
        }
      }
      if (env.empty()) continue;
      double total = 0.0; for (double w : w_env) total += w;
      for (size_t t = 0; t < (size_t)this->split_try; ++t) {
        double r = rng_runif(0.0, total); size_t i = 0; while (i + 1 < w_env.size() && r > w_env[i]) { r -= w_env[i]; ++i; }
        auto &sel = env[i]; if (sel.right_env - sel.left <= 1) continue; int s_idx = rng_randint(sel.left, sel.right_env);
        double sp = [&]{
          std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf; ensure_order_and_sorted_vals_for_leaf(X, *sel.leaf, k, order_cf, sorted_vals_cf);
          size_t pos = (size_t)std::max(sel.left, std::min(sel.right_env - 1, s_idx)); return sorted_vals_cf[pos];
        }();

        size_t ns = 0, nb = 0; std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
        for (int ind : sel.leaf->individuals) {
          if (X[ind][k] < sp) { ++ns; for (size_t p = 0; p < this->value_size; ++p) sum_s_adj[p] += Y[ind][p]; }
          else { ++nb; for (size_t p = 0; p < this->value_size; ++p) sum_b_adj[p] += Y[ind][p]; }
        }
        if (ns == 0 || nb == 0) continue; double loss = 0.0;
        for (size_t p = 0; p < this->value_size; ++p) { loss -= (sum_s_adj[p] * sum_s_adj[p]) / (double)ns; loss -= (sum_b_adj[p] * sum_b_adj[p]) / (double)nb; }
        if (loss < min_split.min_sum) { min_split.min_sum = loss; min_split.tree_index = sel.tree; min_split.leaf_index = sel.leaf; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj; }
      }
    } else {
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty()) continue;
        for (auto &leaf : src_tree->leaves) {
          std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf; ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
          if ((int)sorted_vals_cf.size() <= 2 * leaf_size) continue; int left = leaf_size; int right = (int)sorted_vals_cf.size() - leaf_size;
          std::vector<int> samples = compute_even_spread_indices(left, right, (size_t)this->split_try);
          for (int s_idx : samples) {
            const double sp = sorted_vals_cf[(size_t)s_idx]; size_t ns = 0, nb = 0; std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
            for (int ind : leaf.individuals) { if (X[ind][k] < sp) { ++ns; for (size_t p = 0; p < this->value_size; ++p) sum_s_adj[p] += Y[ind][p]; } else { ++nb; for (size_t p = 0; p < this->value_size; ++p) sum_b_adj[p] += Y[ind][p]; } }
            if (ns == 0 || nb == 0) continue; double loss = 0.0; for (size_t p = 0; p < this->value_size; ++p) { loss -= (sum_s_adj[p] * sum_s_adj[p]) / (double)ns; loss -= (sum_b_adj[p] * sum_b_adj[p]) / (double)nb; }
            if (loss < min_split.min_sum) { min_split.min_sum = loss; min_split.tree_index = src_tree; min_split.leaf_index = &leaf; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj; }
          }
        }
      }
    }
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}


