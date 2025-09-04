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
      auto ensure_weights_cache = [&](const std::shared_ptr<DecisionTree>& tree, int kdim){
        // Lazy-size vectors to feature_size once
        if ((int)tree->weights_epoch_by_dim_v.size() < this->feature_size) {
          tree->weights_epoch_by_dim_v.assign((size_t)this->feature_size, -1);
          tree->fenwick_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
          tree->leaf_weights_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
          tree->weights_total_by_dim_v.assign((size_t)this->feature_size, 0.0);
        }
        // Recompute BIT if epoch mismatches or size changed
        bool need = true;
        if (tree->weights_epoch_by_dim_v[(size_t)kdim] == tree->weights_epoch) {
          if (tree->fenwick_by_dim_v[(size_t)kdim].size() == tree->leaves.size()) need = false;
        }
        if (!need) return; // cache fresh
        const size_t L = tree->leaves.size();
        std::vector<double> bit(L, 0.0), wts(L, 0.0);
        double total = 0.0;
        for (size_t li = 0; li < L; ++li) {
          auto &leaf = tree->leaves[li];
          // Determine number of unique thresholds available in this leaf for kdim
          size_t unique_count = 0;
          auto it_uc = leaf.unique_count_cache.find(kdim);
          if (it_uc != leaf.unique_count_cache.end()) {
            unique_count = it_uc->second;
          } else {
            // Build or reuse sorted values, then count uniques
            std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf;
            ensure_order_and_sorted_vals_for_leaf(X, leaf, kdim, order_cf, sorted_vals_cf);
            if (!sorted_vals_cf.empty()) {
              unique_count = 1;
              for (size_t i = 1; i < sorted_vals_cf.size(); ++i)
                if (sorted_vals_cf[i] != sorted_vals_cf[i - 1]) ++unique_count;
            }
            leaf.unique_count_cache[kdim] = unique_count;
          }
          // Weight = number of unique thresholds that respect min leaf size
          const long width_unique = (long)unique_count - 2L * (long)leaf_size;
          const double w = (width_unique > 0L) ? static_cast<double>(width_unique) : 0.0;
          wts[li] = w; total += w; if (w != 0.0) rpf_utils::fenwick_add(bit, li + 1, w);
        }
        tree->fenwick_by_dim_v[(size_t)kdim] = std::move(bit);
        tree->leaf_weights_by_dim_v[(size_t)kdim] = std::move(wts);
        tree->weights_total_by_dim_v[(size_t)kdim] = total;
        tree->weights_epoch_by_dim_v[(size_t)kdim] = tree->weights_epoch;
      };

      struct SourceInfo { std::shared_ptr<DecisionTree> tree; double total; };
      std::vector<SourceInfo> src_info; src_info.reserve(sources.size());
      double total_all = 0.0;
      for (const auto &src_tree : sources) {
        if (!src_tree || src_tree->leaves.empty()) continue;
        ensure_weights_cache(src_tree, k);
        double tot = src_tree->weights_total_by_dim_v[(size_t)k];
        if (tot <= 0.0) continue;
        src_info.push_back({src_tree, tot});
        total_all += tot;
      }
      if (src_info.empty() || total_all <= 0.0) continue;

      for (size_t t = 0; t < (size_t)this->split_try; ++t) {
        // Sample a source tree proportionally to its total weight
        double r_src = rng_runif(0.0, total_all);
        size_t si = 0;
        while (si + 1 < src_info.size() && r_src > src_info[si].total) { r_src -= src_info[si].total; ++si; }
        auto &sel_tree = src_info[si].tree;

        // Sample a leaf within the selected tree using prefix sums
        const auto &bit = sel_tree->fenwick_by_dim_v[(size_t)k];
        double tot_leaf = src_info[si].total;
        if (bit.empty() || tot_leaf <= 0.0) continue;
        double r_leaf = rng_runif(0.0, tot_leaf);
        size_t leaf_idx_sel = rpf_utils::fenwick_find_by_prefix(bit, r_leaf);
        if (leaf_idx_sel == 0) continue; // safety
        leaf_idx_sel -= 1; // to 0-based
        if (leaf_idx_sel >= sel_tree->leaves.size()) continue;
        Leaf *leaf_ptr = &sel_tree->leaves[leaf_idx_sel];
        // Sample by unique thresholds: build/reuse unique values for this leaf and dim
        std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf;
        ensure_order_and_sorted_vals_for_leaf(X, *leaf_ptr, k, order_cf, sorted_vals_cf);
        size_t unique_count = 0; std::vector<double>* unique_ptr = nullptr;
        if (leaf_ptr->unique_vals_cache.count(k)) {
          unique_ptr = &leaf_ptr->unique_vals_cache[k];
          unique_count = unique_ptr->size();
          leaf_ptr->unique_count_cache[k] = unique_count;
        } else {
          auto uniques = compute_unique_sorted_values(sorted_vals_cf);
          unique_count = uniques.size();
          leaf_ptr->unique_count_cache[k] = unique_count;
          leaf_ptr->unique_vals_cache[k] = std::move(uniques);
          unique_ptr = &leaf_ptr->unique_vals_cache[k];
        }
        const int left = leaf_size; const int right_exclusive = (int)unique_count - leaf_size;
        if (right_exclusive - left <= 1) continue;
        int s_idx = rng_randint(left, right_exclusive);
        double sp = (*unique_ptr)[(size_t)s_idx];

        size_t ns = 0, nb = 0; std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
        for (int ind : leaf_ptr->individuals) {
          if (X[ind][k] < sp) { ++ns; for (size_t p = 0; p < this->value_size; ++p) sum_s_adj[p] += Y[ind][p]; }
          else { ++nb; for (size_t p = 0; p < this->value_size; ++p) sum_b_adj[p] += Y[ind][p]; }
        }
        if (ns == 0 || nb == 0) continue; double loss = 0.0;
        for (size_t p = 0; p < this->value_size; ++p) { loss -= (sum_s_adj[p] * sum_s_adj[p]) / (double)ns; loss -= (sum_b_adj[p] * sum_b_adj[p]) / (double)nb; }
        if (loss < min_split.min_sum) { min_split.min_sum = loss; min_split.tree_index = sel_tree; min_split.leaf_index = leaf_ptr; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj; }
      }
    } else {
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty()) continue;
        for (auto &leaf : src_tree->leaves) {
          std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf; ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
          // Build/reuse unique values for deterministic sampling across unique thresholds
          size_t unique_count = 0; std::vector<double>* unique_ptr = nullptr;
          if (leaf.unique_vals_cache.count(k)) {
            unique_ptr = &leaf.unique_vals_cache[k];
            unique_count = unique_ptr->size();
            leaf.unique_count_cache[k] = unique_count;
          } else {
            auto uniques = compute_unique_sorted_values(sorted_vals_cf);
            unique_count = uniques.size();
            leaf.unique_count_cache[k] = unique_count;
            leaf.unique_vals_cache[k] = std::move(uniques);
            unique_ptr = &leaf.unique_vals_cache[k];
          }
          if ((int)unique_count <= 2 * leaf_size) continue; int left = leaf_size; int right = (int)unique_count - leaf_size;
          std::vector<int> samples = compute_even_spread_indices(left, right, (size_t)this->split_try);
          for (int s_idx : samples) {
            const double sp = (*unique_ptr)[(size_t)s_idx]; size_t ns = 0, nb = 0; std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
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


