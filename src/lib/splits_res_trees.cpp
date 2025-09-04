// Split-mode: res_trees. Operates on the pool of resulting trees constructed
// by expanding dimension sets, evaluating one threshold per leaf via prefix sums.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

Split RandomPlantedForest::calcOptimalSplit_resTrees(const std::vector<std::vector<double>> &Y,
                                                     const std::vector<std::vector<double>> &X,
                                                     std::vector<ResultingTreeCandidate> &possible_trees,
                                                     TreeFamily &curr_family)
{
  Split curr_split, min_split; min_split.min_sum = std::numeric_limits<double>::infinity(); curr_split.Y = &Y;

  if (possible_trees.empty()) return min_split;
  unsigned int raw_candidates = (unsigned int)std::ceil(this->t_try * possible_trees.size());
  unsigned int upper = std::min<unsigned int>((unsigned int)this->max_candidates_, (unsigned int)possible_trees.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));

  std::vector<double> weights(possible_trees.size());
  for (size_t i = 0; i < possible_trees.size(); ++i) weights[i] = std::exp(-this->split_decay_rate_ * possible_trees[i].age);
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
    if (idx >= possible_trees.size()) continue;
    auto &cand = possible_trees[idx];
    auto treePtr = cand.tree; if (!treePtr) continue;

    // Ensure per-tree, per-dimension weight caches (Fenwick + totals) like cur_trees_1
    auto ensure_weights_cache = [&](const std::shared_ptr<DecisionTree>& tree, int kdim){
      if (!tree) return;
      // Lazy-size vectors to feature_size once
      if ((int)tree->weights_epoch_by_dim_v.size() < this->feature_size) {
        tree->weights_epoch_by_dim_v.assign((size_t)this->feature_size, -1);
        tree->fenwick_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
        tree->leaf_weights_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
        tree->weights_total_by_dim_v.assign((size_t)this->feature_size, 0.0);
      }
      bool need = true;
      if (tree->weights_epoch_by_dim_v[(size_t)kdim] == tree->weights_epoch) {
        if (tree->fenwick_by_dim_v[(size_t)kdim].size() == tree->leaves.size()) need = false;
      }
      if (!need) return;
      const int k = kdim; const size_t L = tree->leaves.size();
      std::vector<double> bit(L, 0.0), wts(L, 0.0);
      double total = 0.0;
      const int leaf_size_local = this->n_leaves[k];
      for (size_t li = 0; li < L; ++li) {
        auto &leaf = tree->leaves[li];
        size_t unique_count = 0;
        auto it_uc = leaf.unique_count_cache.find(kdim);
        if (it_uc != leaf.unique_count_cache.end()) {
          unique_count = it_uc->second;
        } else {
          std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf;
          ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
          if (!sorted_vals_cf.empty()) {
            unique_count = 1;
            for (size_t i = 1; i < sorted_vals_cf.size(); ++i)
              if (sorted_vals_cf[i] != sorted_vals_cf[i - 1]) ++unique_count;
          }
          leaf.unique_count_cache[kdim] = unique_count;
        }
        const long width_unique = (long)unique_count - 2L * (long)leaf_size_local;
        const double w = (width_unique > 0L) ? static_cast<double>(width_unique) : 0.0;
        wts[li] = w; total += w; if (w != 0.0) rpf_utils::fenwick_add(bit, li + 1, w);
      }
      tree->fenwick_by_dim_v[(size_t)kdim] = std::move(bit);
      tree->leaf_weights_by_dim_v[(size_t)kdim] = std::move(wts);
      tree->weights_total_by_dim_v[(size_t)kdim] = total;
      tree->weights_epoch_by_dim_v[(size_t)kdim] = tree->weights_epoch;
    };

    // Per-candidate local state for leaves and dimensions used
    struct LeafDimState {
      bool initialized=false; int left=0; int right=0; size_t used_count=0; std::vector<char> used_flags;
      std::vector<size_t> order_cf; std::vector<double> sorted_vals; std::vector<double>* unique_ptr=nullptr; std::vector<std::vector<double>> prefix_cf; std::vector<double> total_cf;
    };
    // Keyed by Leaf* then by k (dimension index)
    std::unordered_map<Leaf*, std::unordered_map<int, LeafDimState>> local_states;

    // Buckets over (kdim, src_tree) with lazy local mutable copies of weights
    struct Bucket { int kdim; std::shared_ptr<DecisionTree> tree; const std::vector<double>* bit_src=nullptr; const std::vector<double>* wts_src=nullptr; std::vector<double> bit; std::vector<double> wts; bool has_local=false; double total=0.0; };
    std::vector<Bucket> buckets;
    size_t grand_total_remaining = 0;

    for (int kdim : treePtr->split_dims) {
      if (kdim == 0) continue; const int k = kdim - 1; const int leaf_size = this->n_leaves[k];

      std::vector<std::shared_ptr<DecisionTree>> sources; sources.reserve(2);
      std::set<int> S = treePtr->split_dims; S.erase(kdim);
      if (S.empty()) { if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end()) sources.push_back(itZero->second); }
      else { if (auto itS = curr_family.find(S); itS != curr_family.end()) sources.push_back(itS->second); }
      if (auto itD = curr_family.find(treePtr->split_dims); itD != curr_family.end()) if (sources.empty() || sources.back().get() != itD->second.get()) sources.push_back(itD->second);

      for (const auto &src_tree : sources) {
        if (!src_tree || src_tree->leaves.empty()) continue;
        ensure_weights_cache(src_tree, k);
        double tot = src_tree->weights_total_by_dim_v[(size_t)k];
        if (tot <= 0.0) continue;
        Bucket b; b.kdim = kdim; b.tree = src_tree; b.total = tot;
        b.bit_src = &src_tree->fenwick_by_dim_v[(size_t)k];
        b.wts_src = &src_tree->leaf_weights_by_dim_v[(size_t)k];
        buckets.push_back(std::move(b));
        grand_total_remaining += (size_t)std::llround(std::max(0.0, tot));
        (void)leaf_size; // silence unused if compiled with warnings
      }
    }
    if (buckets.empty() || grand_total_remaining == 0) continue;

    // Fenwick over bucket totals for O(log B) selection and updates
    std::vector<double> bucket_bit(buckets.size(), 0.0);
    for (size_t i=0;i<buckets.size(); ++i) if (buckets[i].total > 0.0) rpf_utils::fenwick_add(bucket_bit, i+1, buckets[i].total);
    auto fenwick_prefix_sum = [&](const std::vector<double>& bit, size_t idx1)->double { double s=0.0; while (idx1>0) { s += bit[idx1-1]; idx1 -= idx1 & (~idx1 + 1); } return s; };

    const double total_all0 = std::accumulate(buckets.begin(), buckets.end(), 0.0, [](double s, const Bucket& b){ return s + std::max(0.0, b.total); });
    double bucket_total_all = total_all0;
    size_t draws = std::min((size_t)this->split_try, (size_t)std::llround(total_all0));

    for (size_t t=0; t<draws; ++t) {
      // Select a bucket according to current totals
      size_t b_idx;
      if (this->deterministic) {
        if (total_all0 <= 0.0) break;
        double step = total_all0 / (double)draws;
        double target = step * (t + 0.5); if (target >= total_all0) target = std::max(0.0, total_all0 - 1.0);
        b_idx = rpf_utils::fenwick_find_by_prefix(bucket_bit, target);
        if (b_idx == 0) continue; --b_idx;
      } else {
        if (bucket_total_all <= 0.0) break;
        double r = rng_runif(0.0, bucket_total_all);
        b_idx = rpf_utils::fenwick_find_by_prefix(bucket_bit, r);
        if (b_idx == 0) continue; --b_idx;
      }

      auto &bucket = buckets[b_idx];
      const int kdim = bucket.kdim; const int k = kdim - 1;
      if (bucket.total <= 0.0) { continue; }

      // Sample a leaf in the bucket via local Fenwick
      double r_leaf;
      if (this->deterministic) {
        double step = (total_all0 <= 0.0) ? 0.0 : (total_all0 / (double)draws);
        double target_global = step * (t + 0.5); if (target_global >= total_all0) target_global = std::max(0.0, total_all0 - 1.0);
        double before = fenwick_prefix_sum(bucket_bit, b_idx);
        double inside = target_global - before; if (inside < 0.0) inside = 0.0; if (inside >= bucket.total) inside = std::max(0.0, bucket.total - 1.0);
        r_leaf = inside;
      } else {
        r_leaf = rng_runif(0.0, std::max(0.0, bucket.total));
      }
      const std::vector<double>& bit_view = bucket.has_local ? bucket.bit : *(bucket.bit_src);
      size_t leaf_idx_sel = rpf_utils::fenwick_find_by_prefix(bit_view, r_leaf);
      if (leaf_idx_sel == 0) continue; leaf_idx_sel -= 1; size_t wts_size = bucket.has_local ? bucket.wts.size() : (bucket.wts_src ? bucket.wts_src->size() : 0); if (leaf_idx_sel >= wts_size || leaf_idx_sel >= bucket.tree->leaves.size()) continue;
      Leaf *leaf_ptr = &bucket.tree->leaves[leaf_idx_sel];

      // Prepare local per-leaf-per-dim state lazily
      auto &state = local_states[leaf_ptr][k];
      if (!state.initialized) {
        // Build order and sorted values
        ensure_order_and_sorted_vals_for_leaf(X, *leaf_ptr, k, state.order_cf, state.sorted_vals);
        // Compute or reuse unique values and left/right bounds
        size_t unique_count = 0;
        if (leaf_ptr->unique_vals_cache.count(k)) {
          state.unique_ptr = &leaf_ptr->unique_vals_cache[k];
          unique_count = state.unique_ptr->size();
          leaf_ptr->unique_count_cache[k] = unique_count;
        } else {
          auto uniques = compute_unique_sorted_values(state.sorted_vals);
          unique_count = uniques.size();
          leaf_ptr->unique_count_cache[k] = unique_count;
          leaf_ptr->unique_vals_cache[k] = std::move(uniques);
          state.unique_ptr = &leaf_ptr->unique_vals_cache[k];
        }
        const int leaf_size_here = this->n_leaves[k];
        state.left = leaf_size_here; state.right = (int)unique_count - leaf_size_here;
        if (state.right < state.left) { state.left = 0; state.right = 0; }
        if (state.right > state.left) state.used_flags.assign((size_t)(state.right - state.left), 0);
        // Build prefix sums for fast evaluation
        build_prefix_and_total_given_order(Y, *leaf_ptr, state.order_cf, this->value_size, state.prefix_cf, state.total_cf);
        state.initialized = true;
      }
      const std::vector<double>& wts_view = bucket.has_local ? bucket.wts : *(bucket.wts_src);
      if (state.right <= state.left || wts_view[leaf_idx_sel] <= 0.0) { continue; }

      // Select threshold index within [left, right) avoiding repeats
      int s_idx;
      if (this->deterministic) {
        int range = state.right - state.left;
        int remaining_here = (int)wts_view[leaf_idx_sel];
        int guess = state.left + (int)(((double)state.used_count + 0.5) / ((double)remaining_here + 0.5) * range);
        if (guess >= state.right) guess = state.right - 1;
        int lo = guess, hi = guess; bool found = false;
        while (lo >= state.left || hi < state.right) {
          if (lo >= state.left && (state.used_flags.empty() || !state.used_flags[lo - state.left])) { s_idx = lo; found = true; break; }
          if (hi < state.right && (state.used_flags.empty() || !state.used_flags[hi - state.left])) { s_idx = hi; found = true; break; }
          --lo; ++hi;
        }
        if (!found) {
          for (int p = state.left; p < state.right; ++p) {
            if (state.used_flags.empty() || !state.used_flags[p - state.left]) { s_idx = p; break; }
          }
        }
      } else {
        do { s_idx = rng_randint(state.left, state.right); } while (!state.used_flags.empty() && state.used_flags[s_idx - state.left]);
      }
      if (!state.used_flags.empty()) state.used_flags[(size_t)(s_idx - state.left)] = 1;
      state.used_count += 1;

      if (!bucket.has_local) { bucket.bit = *(bucket.bit_src); bucket.wts = *(bucket.wts_src); bucket.has_local = true; }
      // Evaluate loss at chosen threshold
      double sp = (*state.unique_ptr)[(size_t)s_idx];
      const size_t m_eval = leaf_ptr->individuals.size();
      size_t pos_in_sorted = static_cast<size_t>(std::lower_bound(state.sorted_vals.begin(), state.sorted_vals.end(), sp) - state.sorted_vals.begin());
      if (pos_in_sorted == 0 || pos_in_sorted >= m_eval) { continue; }

      double loss = 0.0; std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
      for (size_t p = 0; p < this->value_size; ++p) {
        const double sum_s_base = state.prefix_cf[p][pos_in_sorted - 1]; const double sum_b_base = state.total_cf[p] - sum_s_base;
        sum_s_adj[p] = sum_s_base; sum_b_adj[p] = sum_b_base;
        loss -= (sum_s_adj[p] * sum_s_adj[p]) / static_cast<double>(pos_in_sorted);
        loss -= (sum_b_adj[p] * sum_b_adj[p]) / static_cast<double>(m_eval - pos_in_sorted);
      }
      if (loss < min_split.min_sum) {
        min_split.min_sum = loss; min_split.tree_index = bucket.tree; min_split.leaf_index = leaf_ptr; min_split.split_coordinate = kdim; min_split.split_point = sp; best_idx = (int)idx; min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj;
      }

      // Consume one threshold from this leaf locally: update BIT, wts, totals
      bucket.wts[leaf_idx_sel] -= 1.0; if (bucket.wts[leaf_idx_sel] < 0.0) bucket.wts[leaf_idx_sel] = 0.0;
      rpf_utils::fenwick_add(bucket.bit, leaf_idx_sel + 1, -1.0);
      bucket.total -= 1.0; if (bucket.total < 0.0) bucket.total = 0.0;
      rpf_utils::fenwick_add(bucket_bit, b_idx + 1, -1.0);
      bucket_total_all -= 1.0; if (bucket_total_all < 0.0) bucket_total_all = 0.0;
    }
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_trees);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}
