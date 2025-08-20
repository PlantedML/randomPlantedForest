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
    // Evaluate splits only along the dimensions present in this candidate tree
    for (int kdim : treePtr->split_dims) {
      if (kdim == 0) continue; // skip null-dim marker if present
      const int k = kdim - 1;
      const int leaf_size = this->n_leaves[k];

      // Determine source trees: predecessor S = D'\{k}, null (if S empty), and current D'
      std::vector<std::shared_ptr<DecisionTree>> sources; sources.reserve(2);
      std::set<int> S = treePtr->split_dims; S.erase(kdim);
      if (S.empty()) { if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end()) sources.push_back(itZero->second); }
      else { if (auto itS = curr_family.find(S); itS != curr_family.end()) sources.push_back(itS->second); }
      if (auto itD = curr_family.find(treePtr->split_dims); itD != curr_family.end()) if (sources.empty() || sources.back().get() != itD->second.get()) sources.push_back(itD->second);

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
          ensure_order_and_sorted_vals_for_leaf(X, *sel.leaf, k, sel.order_cf, sel.sorted_vals);
          sel.unique.clear(); sel.unique.reserve(sel.sorted_vals.size()); if (!sel.sorted_vals.empty()) { sel.unique.push_back(sel.sorted_vals[0]); for (size_t i=1;i<sel.sorted_vals.size(); ++i) if (sel.sorted_vals[i] != sel.unique.back()) sel.unique.push_back(sel.sorted_vals[i]); }
          const size_t mloc = sel.leaf->individuals.size();
          build_prefix_and_total_given_order(Y, *sel.leaf, sel.order_cf, this->value_size, sel.prefix_cf, sel.total_cf);
          sel.initialized = true;
        }

        double sp = sel.unique[(size_t)s_idx];
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
          min_split.tree_index = sel.tree; min_split.leaf_index = sel.leaf; min_split.split_coordinate = kdim; min_split.split_point = sp; best_idx = (int)idx;
          min_split.sum_s = sum_s_adj; min_split.sum_b = sum_b_adj;
        }
        if (total_remaining == 0) break;
      }
    }
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_trees);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}


