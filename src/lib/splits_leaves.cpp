// Split-mode: leaves. Evaluates per-leaf candidate splits using cached
// per-leaf orders and prefix sums, with age-weighted candidate sampling.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

Split RandomPlantedForest::calcOptimalSplit_leaves(const std::vector<std::vector<double>> &Y,
                                                    const std::vector<std::vector<double>> &X,
                                                    std::vector<SplitCandidate> &possible_splits,
                                                    TreeFamily &curr_family)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  if (possible_splits.empty()) return min_split;

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

    std::vector<size_t> order; std::vector<double> sorted_vals;
    ensure_order_and_sorted_vals_for_leaf(X, *leafPtr, k, order, sorted_vals);
    std::vector<double> unique = compute_unique_sorted_values(sorted_vals);

    if (unique.size() < 2 * static_cast<size_t>(leaf_size)) continue;

    std::vector<int> samples;
    int left = leaf_size; int right = (int)unique.size() - leaf_size;
    samples = this->deterministic ? compute_even_spread_indices(left, right, (size_t)this->split_try)
                                  : sample_unique_ints_uniform_R(left, right, (size_t)this->split_try);

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

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}


