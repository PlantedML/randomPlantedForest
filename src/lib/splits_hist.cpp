// Split-mode: histogram-binned evaluation (mode 4). Mirrors leaves mode but
// evaluates candidate thresholds at per-feature global bin boundaries.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

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
    const bool have_cached = (tls_working_bin_id_ptr != nullptr && (size_t)k < tls_working_bin_id_ptr->size());
    if (have_cached) {
      const std::vector<int> &bin_k = (*tls_working_bin_id_ptr)[(size_t)k];
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

    // Build prefix across bins then sample only split_try boundaries
    const int total_n = (int)m;
    std::vector<double> total_sum(this->value_size, 0.0);
    for (size_t b = 0; b < Kf; ++b) {
      for (size_t p = 0; p < this->value_size; ++p) total_sum[p] += sum[b][p];
    }
    std::vector<int> prefix_cnt(Kf, 0);
    std::vector<std::vector<double>> prefix_sum(Kf, std::vector<double>(this->value_size, 0.0));
    for (size_t b = 0; b < Kf; ++b) {
      prefix_cnt[b] = cnt[b] + (b > 0 ? prefix_cnt[b - 1] : 0);
      for (size_t p = 0; p < this->value_size; ++p)
        prefix_sum[b][p] = sum[b][p] + (b > 0 ? prefix_sum[b - 1][p] : 0.0);
    }

    // Valid boundaries are b_left in [0, Kf-2] such that both sides satisfy leaf_min
    int first_valid = -1, last_valid = -1;
    for (size_t b_left = 0; b_left + 1 <= Kf - 1; ++b_left) {
      int ln = prefix_cnt[b_left];
      int rn = total_n - ln;
      if (ln >= leaf_min && rn >= leaf_min) {
        if (first_valid < 0) first_valid = (int)b_left;
        last_valid = (int)b_left;
      }
    }
    if (first_valid < 0 || last_valid < first_valid) continue;

    // Sample boundary indices within [first_valid, last_valid]
    std::vector<int> samples = this->deterministic
      ? compute_even_spread_indices(first_valid, last_valid + 1, (size_t)this->split_try)
      : sample_unique_ints_uniform_R(first_valid, last_valid + 1, (size_t)this->split_try);

    for (size_t si = 0; si < samples.size(); ++si) {
      int b_left = samples[si];
      if (b_left < first_valid || b_left > last_valid) continue;
      int left_n = prefix_cnt[(size_t)b_left];
      int right_n = total_n - left_n;
      if (left_n < leaf_min || right_n < leaf_min) continue;

      double loss = 0.0;
      for (size_t p = 0; p < this->value_size; ++p) {
        double ls = prefix_sum[(size_t)b_left][p];
        double rs = total_sum[p] - ls;
        loss -= (ls * ls) / (double)left_n;
        loss -= (rs * rs) / (double)right_n;
      }
      if (loss < min_split.min_sum) {
        min_split.min_sum = loss;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        // Map boundary index to actual split point using precomputed cuts
        double sp = 0.0;
        if (k >= 0 && k < (int)feature_cut_points_.size() && !feature_cut_points_[k].empty()) {
          const auto &cuts = feature_cut_points_[k];
          size_t cp_idx = (size_t)std::min<size_t>((size_t)b_left, cuts.size() - 1);
          sp = cuts[cp_idx];
        } else {
          sp = 0.5 * (leafPtr->intervals[k].first + leafPtr->intervals[k].second);
        }
        min_split.split_point = sp;
        best_idx = (int)idx;
        // Store sums for this boundary
        min_split.sum_s.assign(this->value_size, 0.0);
        min_split.sum_b.assign(this->value_size, 0.0);
        for (size_t p = 0; p < this->value_size; ++p) { min_split.sum_s[p] = prefix_sum[(size_t)b_left][p]; min_split.sum_b[p] = total_sum[p] - prefix_sum[(size_t)b_left][p]; }
      }
    }
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}


