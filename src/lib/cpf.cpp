
#include "cpf.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include "internal_utils.hpp"
using namespace rpf_utils;
#include <limits>
#include <unordered_set>

// ----------------- rpf subclass for classification -----------------

/**
 * \brief Create a prediction model based on Random Forests for classification data sets.
 */

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// Maps a loss-name string to `loss` and `calcLoss` (verbatim move of the
// former inline if/else chain from the data-fitting constructor).
void ClassificationRPF::set_loss_function(const std::string &loss)
{
  if (loss == "L1") {
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  } else if (loss == "L2") {
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  } else if (loss == "median") {
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  } else if (loss == "logit") {
    this->loss = LossType::logit;
    this->calcLoss = &ClassificationRPF::logit_loss;
  } else if (loss == "logit_2") {
    this->loss = LossType::logit_2;
    this->calcLoss = &ClassificationRPF::logit_loss_2;
  } else if (loss == "logit_3") {
    this->loss = LossType::logit_3;
    this->calcLoss = &ClassificationRPF::logit_loss_3;
  } else if (loss == "logit_4") {
    this->loss = LossType::logit_4;
    this->calcLoss = &ClassificationRPF::logit_loss_4;
  } else if (loss == "exponential") {
    this->loss = LossType::exponential;
    this->calcLoss = &ClassificationRPF::exponential_loss;
  } else if (loss == "exponential_2") {
    this->loss = LossType::exponential_2;
    this->calcLoss = &ClassificationRPF::exponential_loss_2;
  } else if (loss == "exponential_3") {
    this->loss = LossType::exponential_3;
    this->calcLoss = &ClassificationRPF::exponential_loss_3;
  } else {
    throw std::invalid_argument("Unknown loss function '" + loss + "'.");
  }
}

ClassificationRPF::ClassificationRPF(const std::string &loss, const RPFParams &params, double delta, double epsilon)
    : RandomPlantedForest(params), delta(delta), epsilon(epsilon)
{
  set_loss_function(loss);
}

// The fit of one leaf holding the rows of a family that made no split (fit_null_leaf()): the leaf update of
// create_tree_family() for that leaf, with the initial weights (0 for logit, 1 for exponential). For L1 and L2
// that is the leaf's mean of the class indicators, for logit and exponential the link of the class shares,
// truncated to [epsilon, 1 - epsilon] as the update truncates them. Each loss's prediction is then the leaf's
// class shares (truncated for logit and exponential).
std::vector<double>
ClassificationRPF::null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows)
{
  if (loss == LossType::median)
    return calcMedian(Y, rows);
  std::vector<double> M = calcMean(Y, rows);
  if (loss == LossType::L1 || loss == LossType::L2)
    return M;
  // the class shares: logit codes the classes as 0 / 1, exponential as -1 / 1
  if (loss == LossType::exponential || loss == LossType::exponential_2 || loss == LossType::exponential_3)
    for (auto &m : M)
      m = (m + 1) / 2;
  const double M_p = 1 - std::accumulate(M.begin(), M.end(), 0.0); // the reference class
  auto truncated = [this](double m) { return std::min(std::max(epsilon, m), 1 - epsilon); };
  auto above_eps = [this](double m) { return std::max(epsilon, m); };
  std::vector<double> value(value_size, 0);
  switch (this->loss) {
  case LossType::logit:
  case LossType::exponential:
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(truncated(M[p]) / truncated(M_p));
    break;
  case LossType::logit_2:
  case LossType::logit_4:
  case LossType::exponential_2:
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(above_eps(M[p]) / above_eps(1 - M[p]));
    break;
  case LossType::logit_3:
  case LossType::exponential_3: {
    double mean_log = log(above_eps(M_p));
    for (size_t p = 0; p < value_size; ++p)
      mean_log += log(above_eps(M[p]));
    mean_log /= (double)(value_size + 1);
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(above_eps(M[p])) - mean_log;
    break;
  }
  default:
    break;
  }
  return value;
}

// ----------------- split modes (classification variants) -----------------
// Each mode mirrors its regression version in splits_*.cpp: the same candidate pool, the same
// sampling of candidates, leaves and thresholds (and therefore the same RNG draws), and the same
// aging. Only the loss differs: every evaluated split is scored with calcLoss, which needs the
// individuals on both sides (I_s / I_b) and the weights W. The sums sum_s / sum_b are computed in
// the same order as in regression, so that with loss "L2" a classification fit reproduces the
// regression fit on the same 0/1 response.

// Draw this round's split candidates from a pool with the given ages, with weights
// exp(-split_decay_rate * age), as the regression split modes do.
std::vector<size_t> ClassificationRPF::sampleCandidates(const std::vector<double> &ages, bool clamp_to_positive)
{
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(t_try * ages.size()));
  unsigned int upper = std::min<size_t>(max_candidates_, ages.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights_vec(ages.size());
  for (size_t i = 0; i < ages.size(); ++i)
    weights_vec[i] = std::exp(-split_decay_rate_ * ages[i]);
  if (clamp_to_positive) { // as in splits_cur_trees_1.cpp
    size_t positive = 0;
    for (double w : weights_vec)
      if (w > 0.0)
        ++positive;
    if (positive == 0)
      n_candidates = 1;
    else if (n_candidates > positive)
      n_candidates = (unsigned int)positive;
  }
  std::vector<size_t> sample_idxs;
  if (deterministic) {
    for (size_t i = 0; i < n_candidates && i < ages.size(); ++i)
      sample_idxs.push_back(i);
  } else
    sample_idxs = sample_weighted_indices_filtered(weights_vec, n_candidates);
  return sample_idxs;
}

// Score splitting `leaf` in dimension k at split_point. The caller sets curr_split.sum_s / sum_b;
// I_s / I_b are filled in the order of leaf.individuals (as finalize_split_from_sums does).
void ClassificationRPF::evaluateSplit(
    Split &curr_split, const Leaf &leaf, int k, double split_point, const std::vector<std::vector<double>> &X
)
{
  curr_split.I_s.clear();
  curr_split.I_b.clear();
  for (int ind : leaf.individuals) {
    if (X[ind][k] < split_point)
      curr_split.I_s.push_back(ind);
    else
      curr_split.I_b.push_back(ind);
  }
  (this->*ClassificationRPF::calcLoss)(curr_split);
}

// Per-tree, per-dimension Fenwick cache of the number of available thresholds of each leaf, used to
// sample leaves in cur_trees_1 and res_trees (same as the lazy cache in splits_cur_trees_1.cpp).
void ClassificationRPF::ensureThresholdWeights(
    const std::shared_ptr<DecisionTree> &tree, int kdim, const std::vector<std::vector<double>> &X
)
{
  if (!tree)
    return;
  if ((int)tree->weights_epoch_by_dim_v.size() < this->feature_size) {
    tree->weights_epoch_by_dim_v.assign((size_t)this->feature_size, -1);
    tree->fenwick_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
    tree->leaf_weights_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
    tree->weights_total_by_dim_v.assign((size_t)this->feature_size, 0.0);
  }
  if (tree->weights_epoch_by_dim_v[(size_t)kdim] == tree->weights_epoch &&
      tree->fenwick_by_dim_v[(size_t)kdim].size() == tree->leaves.size())
    return; // cache fresh
  const size_t L = tree->leaves.size();
  const int leaf_size = this->n_leaves[kdim];
  std::vector<double> bit(L, 0.0), wts(L, 0.0);
  double total = 0.0;
  for (size_t li = 0; li < L; ++li) {
    const double w = leaf_threshold_weight(X, tree->leaves[li], kdim, leaf_size);
    wts[li] = w;
    total += w;
    if (w != 0.0)
      fenwick_add(bit, li + 1, w);
  }
  tree->fenwick_by_dim_v[(size_t)kdim] = std::move(bit);
  tree->leaf_weights_by_dim_v[(size_t)kdim] = std::move(wts);
  tree->weights_total_by_dim_v[(size_t)kdim] = total;
  tree->weights_epoch_by_dim_v[(size_t)kdim] = tree->weights_epoch;
}

// After a split (cur_trees_1 / cur_trees_2): update the cached weight of a replaced leaf, or append the
// weight of a new leaf, in every initialised dimension (as create_tree_family does in rpf.cpp).
void ClassificationRPF::updateThresholdWeight(
    const std::shared_ptr<DecisionTree> &tree, size_t leaf_idx, bool append, const std::vector<std::vector<double>> &X
)
{
  if ((int)tree->fenwick_by_dim_v.size() < this->feature_size)
    return;
  for (int kdim = 0; kdim < this->feature_size; ++kdim) {
    auto &bit = tree->fenwick_by_dim_v[(size_t)kdim];
    auto &wts = tree->leaf_weights_by_dim_v[(size_t)kdim];
    if (bit.empty() || wts.size() != bit.size())
      continue;
    const double w_new = leaf_threshold_weight(X, tree->leaves[leaf_idx], kdim, this->n_leaves[kdim]);
    if (append) {
      fenwick_push_back(bit, w_new);
      wts.push_back(w_new);
      tree->weights_total_by_dim_v[(size_t)kdim] += w_new;
    } else {
      const double delta_w = w_new - (leaf_idx < wts.size() ? wts[leaf_idx] : 0.0);
      if (leaf_idx < wts.size())
        wts[leaf_idx] = w_new;
      if (delta_w != 0.0) {
        fenwick_add(bit, leaf_idx + 1, delta_w);
        tree->weights_total_by_dim_v[(size_t)kdim] += delta_w;
      }
    }
  }
}

// Mode 1: cur_trees_2 (classification variant of splits_cur_trees_2.cpp)
Split ClassificationRPF::calcOptimalSplit_curTrees2(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;

  std::vector<double> ages(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i)
    ages[i] = possible_splits[i].age;
  std::vector<size_t> sample_idxs = sampleCandidates(ages);

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    int leaf_size = n_leaves[k];

    // Leaves of the tree D' \ {k} (or of the null tree) and of the candidate tree D' itself
    std::set<int> tree_dims = it->tree->split_dims;
    tree_dims.erase(k + 1);
    tree_dims.erase(0);

    std::vector<std::shared_ptr<DecisionTree>> curr_trees;
    if (tree_dims.empty()) {
      auto itZero = curr_family.find(std::set<int>{0});
      if (itZero != curr_family.end() && itZero->second)
        curr_trees.push_back(itZero->second);
    }
    if (auto itS = curr_family.find(tree_dims); itS != curr_family.end() && itS->second)
      curr_trees.push_back(itS->second);
    if (auto itD = curr_family.find(it->tree->split_dims); itD != curr_family.end() && itD->second) {
      if (curr_trees.empty() || curr_trees.back().get() != itD->second.get())
        curr_trees.push_back(itD->second);
    }

    for (auto &curr_tree : curr_trees) {
      if (curr_tree->leaves.empty())
        continue;
      for (auto &leaf : curr_tree->leaves) {
        std::vector<size_t> order_cf;
        std::vector<double> sorted_vals_cf;
        ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
        std::vector<double> *unique_ptr = nullptr;
        size_t unique_count = 0;
        if (leaf.cache().unique_count_cache.count(k)) {
          unique_count = leaf.cache().unique_count_cache[k];
          if (unique_count != 0 && leaf.cache().unique_vals_cache.count(k))
            unique_ptr = &leaf.cache().unique_vals_cache[k];
        }
        if (!unique_ptr) {
          auto uniques = compute_unique_sorted_values(sorted_vals_cf);
          unique_count = uniques.size();
          leaf.cache().unique_count_cache[k] = unique_count;
          leaf.cache().unique_vals_cache[k] = std::move(uniques);
          unique_ptr = &leaf.cache().unique_vals_cache[k];
        }
        if (unique_count < 2 * static_cast<size_t>(leaf_size))
          continue;

        const size_t m = leaf.individuals.size();
        std::vector<int> samples;
        if (deterministic) {
          // split_try thresholds evenly spread over the same range the random draws use
          samples = compute_even_spread_indices(leaf_size, (int)unique_count - leaf_size + 1, (size_t)split_try);
          if (samples.empty())
            continue;
        } else {
          samples.resize(split_try);
          for (size_t i = 0; i < samples.size(); ++i)
            samples[i] = rng_randint(leaf_size, (int)unique_count - leaf_size + 1);
          std::sort(samples.begin(), samples.end());
        }
        const bool single_eval = (samples.size() == 1);
        std::vector<std::vector<double>> prefix_cf;
        std::vector<double> total_cf;
        if (!single_eval)
          build_prefix_and_total_given_order(Y, leaf, order_cf, value_size, prefix_cf, total_cf);

        for (size_t si = 0; si < samples.size(); ++si) {
          const double sp = (*unique_ptr)[samples[si]];
          size_t pos = static_cast<size_t>(
              std::lower_bound(sorted_vals_cf.begin(), sorted_vals_cf.end(), sp) - sorted_vals_cf.begin()
          );
          if (pos == 0 || pos >= m)
            continue;
          if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
            continue;
          curr_split.sum_s.assign(value_size, 0.0);
          curr_split.sum_b.assign(value_size, 0.0);
          if (!single_eval) {
            for (size_t p = 0; p < value_size; ++p) {
              curr_split.sum_s[p] = prefix_cf[p][pos - 1];
              curr_split.sum_b[p] = total_cf[p] - curr_split.sum_s[p];
            }
          } else {
            for (int ind : leaf.individuals) {
              if (X[ind][k] < sp) {
                for (size_t p = 0; p < value_size; ++p)
                  curr_split.sum_s[p] += Y[ind][p];
              } else {
                for (size_t p = 0; p < value_size; ++p)
                  curr_split.sum_b[p] += Y[ind][p];
              }
            }
          }
          evaluateSplit(curr_split, leaf, k, sp, X);
          if (curr_split.min_sum < min_split.min_sum) {
            min_split = curr_split;
            min_split.tree_index = curr_tree;
            min_split.leaf_index = &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point = sp;
            best_idx = (int)idx;
          }
        }
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  return min_split;
}

// Mode 3: leaves (classification variant of splits_leaves.cpp)
Split ClassificationRPF::calcOptimalSplit_leaves(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;
  if (possible_splits.empty())
    return min_split;

  std::vector<double> ages(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i)
    ages[i] = possible_splits[i].age;
  std::vector<size_t> sample_idxs = sampleCandidates(ages);

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    if (!it->tree || it->leaf_idx >= it->tree->leaves.size())
      continue;
    Leaf *leafPtr = &it->tree->leaves[it->leaf_idx];

    const int leaf_size = n_leaves[k];
    const size_t m = leafPtr->individuals.size();
    if (m == 0)
      continue;
    if (m < static_cast<size_t>(2 * leaf_size))
      continue;

    std::vector<size_t> order;
    std::vector<double> sorted_vals;
    ensure_order_and_sorted_vals_for_leaf(X, *leafPtr, k, order, sorted_vals);
    std::vector<double> unique = compute_unique_sorted_values(sorted_vals);
    std::vector<size_t> first_pos; // first position of each unique value in sorted_vals
    first_pos.reserve(unique.size());
    if (!sorted_vals.empty()) {
      first_pos.push_back(0);
      for (size_t i = 1; i < sorted_vals.size(); ++i)
        if (sorted_vals[i] != sorted_vals[i - 1])
          first_pos.push_back(i);
    }
    if (unique.size() < 2 * static_cast<size_t>(leaf_size))
      continue;

    int left = leaf_size;
    int right_exclusive = (int)unique.size() - leaf_size + 1;
    std::vector<int> samples = deterministic ? compute_even_spread_indices(left, right_exclusive, (size_t)split_try)
                                             : sample_unique_ints_uniform_R(left, right_exclusive, (size_t)split_try);

    std::vector<std::vector<double>> prefix;
    std::vector<double> total;
    build_prefix_and_total_given_order(Y, *leafPtr, order, value_size, prefix, total);

    for (size_t si = 0; si < samples.size(); ++si) {
      const size_t uidx = static_cast<size_t>(samples[si]);
      if (uidx >= unique.size() || uidx >= first_pos.size())
        continue;
      const double sp = unique[uidx];
      const size_t pos = first_pos[uidx];
      if (pos == 0 || pos >= m)
        continue;
      if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
        continue;
      curr_split.sum_s.assign(value_size, 0.0);
      curr_split.sum_b.assign(value_size, 0.0);
      for (size_t p = 0; p < value_size; ++p) {
        curr_split.sum_s[p] = prefix[p][pos - 1];
        curr_split.sum_b[p] = total[p] - curr_split.sum_s[p];
      }
      evaluateSplit(curr_split, *leafPtr, k, sp, X);
      if (curr_split.min_sum < min_split.min_sum) {
        min_split = curr_split;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        min_split.split_point = sp;
        best_idx = (int)idx;
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  return min_split;
}

// Mode 4: histogram-binned (classification variant of splits_hist.cpp)
Split ClassificationRPF::calcOptimalSplit_hist(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;
  if (possible_splits.empty())
    return min_split;

  std::vector<double> ages(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i)
    ages[i] = possible_splits[i].age;
  std::vector<size_t> sample_idxs = sampleCandidates(ages);

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    if (!it->tree || it->leaf_idx >= it->tree->leaves.size())
      continue;
    const int k = it->dim - 1;
    Leaf *leafPtr = &it->tree->leaves[it->leaf_idx];
    const int leaf_min = n_leaves[k];
    const size_t m = leafPtr->individuals.size();
    if (m == 0)
      continue;

    // Histogram of this leaf over the global bins of feature k (bin = number of cut points <= x)
    const auto &cuts_k =
        (k >= 0 && k < (int)feature_cut_points_.size()) ? feature_cut_points_[k] : std::vector<double>{};
    size_t Kf = cuts_k.size() + 1;
    if (Kf < 2)
      continue;
    std::vector<int> cnt(Kf, 0);
    std::vector<std::vector<double>> sum(Kf, std::vector<double>(value_size, 0.0));
    for (int ind : leafPtr->individuals) {
      double v = X[ind][k];
      int b = 0;
      if (!cuts_k.empty()) {
        auto itb = std::upper_bound(cuts_k.begin(), cuts_k.end(), v);
        b = (int)std::distance(cuts_k.begin(), itb);
        if (b < 0)
          b = 0;
        if ((size_t)b >= Kf)
          b = (int)Kf - 1;
      }
      cnt[(size_t)b] += 1;
      for (size_t p = 0; p < value_size; ++p)
        sum[(size_t)b][p] += Y[ind][p];
    }

    const int total_n = (int)m;
    std::vector<double> total_sum(value_size, 0.0);
    for (size_t b = 0; b < Kf; ++b)
      for (size_t p = 0; p < value_size; ++p)
        total_sum[p] += sum[b][p];
    std::vector<int> prefix_cnt(Kf, 0);
    std::vector<std::vector<double>> prefix_sum(Kf, std::vector<double>(value_size, 0.0));
    for (size_t b = 0; b < Kf; ++b) {
      prefix_cnt[b] = cnt[b] + (b > 0 ? prefix_cnt[b - 1] : 0);
      for (size_t p = 0; p < value_size; ++p)
        prefix_sum[b][p] = sum[b][p] + (b > 0 ? prefix_sum[b - 1][p] : 0.0);
    }

    // Valid boundaries b_left in [0, Kf-2] leave at least leaf_min individuals on both sides
    int first_valid = -1, last_valid = -1;
    for (size_t b_left = 0; b_left + 1 <= Kf - 1; ++b_left) {
      int ln = prefix_cnt[b_left];
      int rn = total_n - ln;
      if (ln >= leaf_min && rn >= leaf_min) {
        if (first_valid < 0)
          first_valid = (int)b_left;
        last_valid = (int)b_left;
      }
    }
    if (first_valid < 0 || last_valid < first_valid)
      continue;

    // Sample split_try boundaries within [first_valid, last_valid]
    std::vector<int> samples = deterministic
                                   ? compute_even_spread_indices(first_valid, last_valid + 1, (size_t)split_try)
                                   : sample_unique_ints_uniform_R(first_valid, last_valid + 1, (size_t)split_try);

    for (size_t si = 0; si < samples.size(); ++si) {
      int b_left = samples[si];
      if (b_left < first_valid || b_left > last_valid)
        continue;
      int left_n = prefix_cnt[(size_t)b_left];
      int right_n = total_n - left_n;
      if (left_n < leaf_min || right_n < leaf_min)
        continue;

      // Boundary b_left splits at its cut point: x < cut <=> bin(x) <= b_left
      double sp;
      if (!cuts_k.empty())
        sp = cuts_k[std::min<size_t>((size_t)b_left, cuts_k.size() - 1)];
      else
        sp = 0.5 * (leafPtr->intervals[k].first + leafPtr->intervals[k].second);

      curr_split.sum_s.assign(value_size, 0.0);
      curr_split.sum_b.assign(value_size, 0.0);
      for (size_t p = 0; p < value_size; ++p) {
        curr_split.sum_s[p] = prefix_sum[(size_t)b_left][p];
        curr_split.sum_b[p] = total_sum[p] - prefix_sum[(size_t)b_left][p];
      }
      evaluateSplit(curr_split, *leafPtr, k, sp, X);
      if (curr_split.min_sum < min_split.min_sum) {
        min_split = curr_split;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        min_split.split_point = sp;
        best_idx = (int)idx;
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  return min_split;
}

// Mode 2: cur_trees_1 (classification variant of splits_cur_trees_1.cpp)
Split ClassificationRPF::calcOptimalSplit_curTrees1(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;

  std::vector<double> ages(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i)
    ages[i] = possible_splits[i].age;
  // (deterministic: the first n_candidates candidates, as in the other modes)
  std::vector<size_t> sample_idxs = sampleCandidates(ages, true);

  // Score the split of `leaf` at sp; sums in the order of leaf.individuals as in regression
  auto try_split = [&](Leaf &leaf, const std::shared_ptr<DecisionTree> &tree, int k, double sp, size_t idx, int &best) {
    curr_split.sum_s.assign(value_size, 0.0);
    curr_split.sum_b.assign(value_size, 0.0);
    size_t ns = 0, nb = 0;
    for (int ind : leaf.individuals) {
      if (X[ind][k] < sp) {
        ++ns;
        for (size_t p = 0; p < value_size; ++p)
          curr_split.sum_s[p] += Y[ind][p];
      } else {
        ++nb;
        for (size_t p = 0; p < value_size; ++p)
          curr_split.sum_b[p] += Y[ind][p];
      }
    }
    if (ns == 0 || nb == 0)
      return;
    evaluateSplit(curr_split, leaf, k, sp, X);
    if (curr_split.min_sum < min_split.min_sum) {
      min_split = curr_split;
      min_split.tree_index = tree;
      min_split.leaf_index = &leaf;
      min_split.split_coordinate = k + 1;
      min_split.split_point = sp;
      best = (int)idx;
    }
  };

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    int leaf_size = n_leaves[k];

    // Source trees: D' \ {k} (or the null tree) and D' itself
    std::set<int> Dprime_minus_k = it->tree->split_dims;
    Dprime_minus_k.erase(k + 1);
    Dprime_minus_k.erase(0);
    std::vector<std::shared_ptr<DecisionTree>> sources;
    sources.reserve(2);
    if (Dprime_minus_k.empty()) {
      if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end())
        sources.push_back(itZero->second);
    } else {
      if (auto itS = curr_family.find(Dprime_minus_k); itS != curr_family.end())
        sources.push_back(itS->second);
    }
    if (auto itD = curr_family.find(it->tree->split_dims); itD != curr_family.end())
      if (sources.empty() || sources.back().get() != itD->second.get())
        sources.push_back(itD->second);

    if (!deterministic) {
      struct SourceInfo {
        std::shared_ptr<DecisionTree> tree;
        double total;
      };
      std::vector<SourceInfo> src_info;
      src_info.reserve(sources.size());
      double total_all = 0.0;
      for (const auto &src_tree : sources) {
        if (!src_tree || src_tree->leaves.empty())
          continue;
        ensureThresholdWeights(src_tree, k, X);
        double tot = src_tree->weights_total_by_dim_v[(size_t)k];
        if (tot <= 0.0)
          continue;
        src_info.push_back({src_tree, tot});
        total_all += tot;
      }
      if (src_info.empty() || total_all <= 0.0)
        continue;

      // split_try draws: a source tree and a leaf proportionally to the leaf's available thresholds,
      // then one threshold uniformly within the leaf
      for (size_t t = 0; t < (size_t)split_try; ++t) {
        double r_src = rng_runif(0.0, total_all);
        size_t si = 0;
        while (si + 1 < src_info.size() && r_src > src_info[si].total) {
          r_src -= src_info[si].total;
          ++si;
        }
        auto &sel_tree = src_info[si].tree;

        const auto &bit = sel_tree->fenwick_by_dim_v[(size_t)k];
        double tot_leaf = src_info[si].total;
        if (bit.empty() || tot_leaf <= 0.0)
          continue;
        double r_leaf = rng_runif(0.0, tot_leaf);
        size_t leaf_idx_sel = fenwick_find_by_prefix(bit, r_leaf);
        if (leaf_idx_sel == 0)
          continue;
        leaf_idx_sel -= 1;
        if (leaf_idx_sel >= sel_tree->leaves.size())
          continue;
        Leaf *leaf_ptr = &sel_tree->leaves[leaf_idx_sel];
        std::vector<size_t> order_cf;
        std::vector<double> sorted_vals_cf;
        ensure_order_and_sorted_vals_for_leaf(X, *leaf_ptr, k, order_cf, sorted_vals_cf);
        size_t unique_count = 0;
        std::vector<double> *unique_ptr = nullptr;
        if (leaf_ptr->cache().unique_vals_cache.count(k)) {
          unique_ptr = &leaf_ptr->cache().unique_vals_cache[k];
          unique_count = unique_ptr->size();
          leaf_ptr->cache().unique_count_cache[k] = unique_count;
        } else {
          auto uniques = compute_unique_sorted_values(sorted_vals_cf);
          unique_count = uniques.size();
          leaf_ptr->cache().unique_count_cache[k] = unique_count;
          leaf_ptr->cache().unique_vals_cache[k] = std::move(uniques);
          unique_ptr = &leaf_ptr->cache().unique_vals_cache[k];
        }
        const int left = leaf_size;
        const int right_exclusive = (int)unique_count - leaf_size + 1;
        if (right_exclusive <= left)
          continue;
        int s_idx = rng_randint(left, right_exclusive);
        try_split(*leaf_ptr, sel_tree, k, (*unique_ptr)[(size_t)s_idx], idx, best_idx);
      }
    } else {
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty())
          continue;
        for (auto &leaf : src_tree->leaves) {
          std::vector<size_t> order_cf;
          std::vector<double> sorted_vals_cf;
          ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
          size_t unique_count = 0;
          std::vector<double> *unique_ptr = nullptr;
          if (leaf.cache().unique_vals_cache.count(k)) {
            unique_ptr = &leaf.cache().unique_vals_cache[k];
            unique_count = unique_ptr->size();
            leaf.cache().unique_count_cache[k] = unique_count;
          } else {
            auto uniques = compute_unique_sorted_values(sorted_vals_cf);
            unique_count = uniques.size();
            leaf.cache().unique_count_cache[k] = unique_count;
            leaf.cache().unique_vals_cache[k] = std::move(uniques);
            unique_ptr = &leaf.cache().unique_vals_cache[k];
          }
          if ((int)unique_count < 2 * leaf_size)
            continue;
          int left = leaf_size;
          int right = (int)unique_count - leaf_size + 1;
          std::vector<int> samples = compute_even_spread_indices(left, right, (size_t)split_try);
          for (int s_idx : samples)
            try_split(leaf, src_tree, k, (*unique_ptr)[(size_t)s_idx], idx, best_idx);
        }
      }
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  return min_split;
}

// Mode 0: res_trees (classification variant of splits_res_trees.cpp)
Split ClassificationRPF::calcOptimalSplit_resTrees(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<RandomPlantedForest::ResultingTreeCandidate> &possible_trees,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;
  if (possible_trees.empty())
    return min_split;

  std::vector<double> ages(possible_trees.size());
  for (size_t i = 0; i < possible_trees.size(); ++i)
    ages[i] = possible_trees[i].age;
  std::vector<size_t> sample_idxs = sampleCandidates(ages);

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    if (idx >= possible_trees.size())
      continue;
    auto treePtr = possible_trees[idx].tree;
    if (!treePtr)
      continue;

    // Per-candidate state of each (leaf, dimension) pair drawn so far
    struct LeafDimState {
      bool initialized = false;
      int left = 0;
      int right = 0;
      size_t used_count = 0;
      std::vector<char> used_flags;
      std::vector<size_t> order_cf;
      std::vector<double> sorted_vals;
      std::vector<double> *unique_ptr = nullptr;
      std::vector<std::vector<double>> prefix_cf;
      std::vector<double> total_cf;
    };
    std::unordered_map<Leaf *, std::unordered_map<int, LeafDimState>> local_states;

    // Buckets over (split dimension k in D', source tree) with lazy local copies of the leaf weights
    struct Bucket {
      int kdim;
      std::shared_ptr<DecisionTree> tree;
      const std::vector<double> *bit_src = nullptr;
      const std::vector<double> *wts_src = nullptr;
      std::vector<double> bit;
      std::vector<double> wts;
      bool has_local = false;
      double total = 0.0;
    };
    std::vector<Bucket> buckets;
    size_t grand_total_remaining = 0;

    for (int kdim : treePtr->split_dims) {
      if (kdim == 0)
        continue;
      const int k = kdim - 1;
      std::vector<std::shared_ptr<DecisionTree>> sources;
      sources.reserve(2);
      std::set<int> S = treePtr->split_dims;
      S.erase(kdim);
      if (S.empty()) {
        if (auto itZero = curr_family.find(std::set<int>{0}); itZero != curr_family.end())
          sources.push_back(itZero->second);
      } else {
        if (auto itS = curr_family.find(S); itS != curr_family.end())
          sources.push_back(itS->second);
      }
      if (auto itD = curr_family.find(treePtr->split_dims); itD != curr_family.end())
        if (sources.empty() || sources.back().get() != itD->second.get())
          sources.push_back(itD->second);

      for (const auto &src_tree : sources) {
        if (!src_tree || src_tree->leaves.empty())
          continue;
        ensureThresholdWeights(src_tree, k, X);
        double tot = src_tree->weights_total_by_dim_v[(size_t)k];
        if (tot <= 0.0)
          continue;
        Bucket b;
        b.kdim = kdim;
        b.tree = src_tree;
        b.total = tot;
        b.bit_src = &src_tree->fenwick_by_dim_v[(size_t)k];
        b.wts_src = &src_tree->leaf_weights_by_dim_v[(size_t)k];
        buckets.push_back(std::move(b));
        grand_total_remaining += (size_t)std::llround(std::max(0.0, tot));
      }
    }
    if (buckets.empty() || grand_total_remaining == 0)
      continue;

    std::vector<double> bucket_bit(buckets.size(), 0.0);
    for (size_t i = 0; i < buckets.size(); ++i)
      if (buckets[i].total > 0.0)
        fenwick_add(bucket_bit, i + 1, buckets[i].total);

    const double total_all0 = std::accumulate(buckets.begin(), buckets.end(), 0.0, [](double s, const Bucket &b) {
      return s + std::max(0.0, b.total);
    });
    double bucket_total_all = total_all0;
    size_t draws = std::min((size_t)split_try, (size_t)std::llround(total_all0));
    // Deterministic: the split_try thresholds evenly spread over all (dimension, leaf, threshold) triples of this
    // candidate, in bucket, leaf and threshold order (the analogue of drawing them uniformly without replacement)
    std::vector<int> det_targets;
    std::vector<double> det_bucket_start;
    if (deterministic) {
      det_targets = compute_even_spread_indices(0, (int)std::llround(total_all0), draws);
      double start = 0.0;
      for (const auto &b : buckets) {
        det_bucket_start.push_back(start);
        start += std::max(0.0, b.total);
      }
    }
    // deterministic: the selected leaf, and its threshold counted from the leaf's left bound
    int det_offset = 0;
    size_t det_leaf = 0;

    // split_try draws of a (dimension, leaf) pair proportionally to its remaining thresholds, then
    // one unused threshold uniformly within it
    for (size_t t = 0; t < draws; ++t) {
      size_t b_idx;
      if (deterministic) {
        if (t >= det_targets.size())
          break;
        const double g = (double)det_targets[t];
        b_idx =
            (size_t)(std::upper_bound(det_bucket_start.begin(), det_bucket_start.end(), g) - det_bucket_start.begin());
        if (b_idx == 0)
          continue;
        --b_idx;
      } else {
        if (bucket_total_all <= 0.0)
          break;
        double r = rng_runif(0.0, bucket_total_all);
        b_idx = fenwick_find_by_prefix(bucket_bit, r);
        if (b_idx == 0)
          continue;
        --b_idx;
      }

      auto &bucket = buckets[b_idx];
      const int kdim = bucket.kdim;
      const int k = kdim - 1;
      if (bucket.total <= 0.0)
        continue;

      double r_leaf;
      if (deterministic) {
        // The leaf holding global threshold g in the bucket's original leaf order (weights before this round)
        double inside = (double)det_targets[t] - det_bucket_start[b_idx];
        const std::vector<double> &wts0 = *(bucket.wts_src);
        size_t li = 0;
        while (li + 1 < wts0.size() && inside >= wts0[li]) {
          inside -= wts0[li];
          ++li;
        }
        det_offset = (int)inside;
        det_leaf = li;
        r_leaf = 0.0;
      } else {
        r_leaf = rng_runif(0.0, std::max(0.0, bucket.total));
      }
      const std::vector<double> &bit_view = bucket.has_local ? bucket.bit : *(bucket.bit_src);
      size_t leaf_idx_sel = deterministic ? det_leaf + 1 : fenwick_find_by_prefix(bit_view, r_leaf);
      if (leaf_idx_sel == 0)
        continue;
      leaf_idx_sel -= 1;
      size_t wts_size = bucket.has_local ? bucket.wts.size() : (bucket.wts_src ? bucket.wts_src->size() : 0);
      if (leaf_idx_sel >= wts_size || leaf_idx_sel >= bucket.tree->leaves.size())
        continue;
      Leaf *leaf_ptr = &bucket.tree->leaves[leaf_idx_sel];

      auto &state = local_states[leaf_ptr][k];
      if (!state.initialized) {
        ensure_order_and_sorted_vals_for_leaf(X, *leaf_ptr, k, state.order_cf, state.sorted_vals);
        size_t unique_count = 0;
        if (leaf_ptr->cache().unique_vals_cache.count(k)) {
          state.unique_ptr = &leaf_ptr->cache().unique_vals_cache[k];
          unique_count = state.unique_ptr->size();
          leaf_ptr->cache().unique_count_cache[k] = unique_count;
        } else {
          auto uniques = compute_unique_sorted_values(state.sorted_vals);
          unique_count = uniques.size();
          leaf_ptr->cache().unique_count_cache[k] = unique_count;
          leaf_ptr->cache().unique_vals_cache[k] = std::move(uniques);
          state.unique_ptr = &leaf_ptr->cache().unique_vals_cache[k];
        }
        const int leaf_size_here = n_leaves[k];
        state.left = leaf_size_here;
        state.right = (int)unique_count - leaf_size_here + 1;
        if (state.right < state.left) {
          state.left = 0;
          state.right = 0;
        }
        if (state.right > state.left)
          state.used_flags.assign((size_t)(state.right - state.left), 0);
        build_prefix_and_total_given_order(Y, *leaf_ptr, state.order_cf, value_size, state.prefix_cf, state.total_cf);
        state.initialized = true;
      }
      const std::vector<double> &wts_view = bucket.has_local ? bucket.wts : *(bucket.wts_src);
      if (state.right <= state.left || wts_view[leaf_idx_sel] <= 0.0)
        continue;

      int s_idx = -1;
      if (deterministic) {
        s_idx = state.left + det_offset;
        if (s_idx >= state.right)
          continue;
      } else {
        do {
          s_idx = rng_randint(state.left, state.right);
        } while (!state.used_flags.empty() && state.used_flags[s_idx - state.left]);
      }
      // every threshold in range already used
      if (s_idx < 0)
        continue;
      if (!state.used_flags.empty())
        state.used_flags[(size_t)(s_idx - state.left)] = 1;
      state.used_count += 1;

      if (!bucket.has_local) {
        bucket.bit = *(bucket.bit_src);
        bucket.wts = *(bucket.wts_src);
        bucket.has_local = true;
      }
      double sp = (*state.unique_ptr)[(size_t)s_idx];
      const size_t m_eval = leaf_ptr->individuals.size();
      size_t pos_in_sorted = static_cast<size_t>(
          std::lower_bound(state.sorted_vals.begin(), state.sorted_vals.end(), sp) - state.sorted_vals.begin()
      );
      if (pos_in_sorted == 0 || pos_in_sorted >= m_eval)
        continue;

      curr_split.sum_s.assign(value_size, 0.0);
      curr_split.sum_b.assign(value_size, 0.0);
      for (size_t p = 0; p < value_size; ++p) {
        curr_split.sum_s[p] = state.prefix_cf[p][pos_in_sorted - 1];
        curr_split.sum_b[p] = state.total_cf[p] - curr_split.sum_s[p];
      }
      evaluateSplit(curr_split, *leaf_ptr, k, sp, X);
      if (curr_split.min_sum < min_split.min_sum) {
        min_split = curr_split;
        min_split.tree_index = bucket.tree;
        min_split.leaf_index = leaf_ptr;
        min_split.split_coordinate = kdim;
        min_split.split_point = sp;
        best_idx = (int)idx;
      }

      // Consume this threshold: update the local leaf weights and the bucket totals
      bucket.wts[leaf_idx_sel] -= 1.0;
      if (bucket.wts[leaf_idx_sel] < 0.0)
        bucket.wts[leaf_idx_sel] = 0.0;
      fenwick_add(bucket.bit, leaf_idx_sel + 1, -1.0);
      bucket.total -= 1.0;
      if (bucket.total < 0.0)
        bucket.total = 0.0;
      fenwick_add(bucket_bit, b_idx + 1, -1.0);
      bucket_total_all -= 1.0;
      if (bucket_total_all < 0.0)
        bucket_total_all = 0.0;
    }
  }

  age_pool_by_sample(sample_idxs, best_idx, possible_trees);
  return min_split;
}

// Dispatcher selecting by split_structure_mode_
Split ClassificationRPF::calcOptimalSplit(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    std::vector<std::vector<double>> &weights
)
{
  if (split_structure_mode_ == 4)
    return this->calcOptimalSplit_hist(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 3)
    return this->calcOptimalSplit_leaves(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 2)
    return this->calcOptimalSplit_curTrees1(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 1)
    return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family, weights);
  return Split{};
}

void ClassificationRPF::create_tree_family(std::vector<Leaf> initial_leaves, size_t n)
{

  TreeFamily curr_family;
  curr_family.insert(
      std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))
  ); // save tree with one leaf in the beginning

  // Seed per mode
  std::vector<SplitCandidate> possible_splits;
  std::vector<ResultingTreeCandidate> possible_trees;
  if (split_structure_mode_ == 0) {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      if (!this->dims_allowed({feature_dim}))
        continue;
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_trees.emplace_back(treePtr);
    }
  } else if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
    auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree> &T, size_t li) {
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
        std::set<int> res_dims = T->split_dims;
        res_dims.insert(feature_dim);
        res_dims.erase(0);
        if (!this->dims_allowed(res_dims))
          continue;
        if (!this->leafCandidateExists(possible_splits, T, li, feature_dim))
          possible_splits.emplace_back(feature_dim, T, li);
      }
    };
    auto null_tree = curr_family[{0}];
    if (!null_tree->leaves.empty())
      add_leaf_candidates(null_tree, 0);
  } else {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      if (!this->dims_allowed({feature_dim}))
        continue;
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_splits.emplace_back(feature_dim, treePtr, (size_t)0);
    }
  }

  // sample data points with replacement
  int sample_index;
  init_inbag(n);
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;

  // deterministic
  if (deterministic) {
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  } else {
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);

    // bagging/subsampling
    for (size_t i = 0; i < (size_t)sample_size; ++i) {
      sample_index = rpf_utils::rng_randint(0, (int)sample_size);
      add_inbag(n, sample_index);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }

  // initialize weights
  std::vector<std::vector<double>> weights;
  switch (this->loss) {
  case LossType::logit:
  case LossType::logit_2:
  case LossType::logit_3:
  case LossType::logit_4:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 0);
    break;
  case LossType::exponential:
  case LossType::exponential_2:
  case LossType::exponential_3:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 1);
    break;
  default:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 0);
  }

  // modify existing or add new trees through splitting
  Split curr_split;
  for (int split_count = 0; split_count < n_splits; ++split_count) {

    // find optimal split
    if (split_structure_mode_ == 0)
      curr_split = this->calcOptimalSplit_resTrees(samples_Y, samples_X, possible_trees, curr_family, weights);
    else
      curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family, weights);

    // continue only if we get a significant result
    if (!std::isinf(curr_split.min_sum)) {

      // update pools by mode
      if (split_structure_mode_ == 0) {
        std::set<int> Dprime = curr_split.tree_index->split_dims;
        Dprime.insert(curr_split.split_coordinate);
        Dprime.erase(0);
        if (!this->resultingTreeExists(possible_trees, Dprime)) {
          if (auto found = treeExists(Dprime, curr_family))
            possible_trees.emplace_back(found);
          else {
            curr_family.insert({Dprime, std::make_shared<DecisionTree>(DecisionTree(Dprime))});
            possible_trees.emplace_back(curr_family[Dprime]);
          }
        }
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> U = Dprime;
          U.insert(feature_dim);
          if (U.size() == Dprime.size())
            continue;
          if (!this->dims_allowed(U))
            continue;
          if (this->resultingTreeExists(possible_trees, U))
            continue;
          if (auto found = treeExists(U, curr_family))
            possible_trees.emplace_back(found);
          else {
            curr_family.insert({U, std::make_shared<DecisionTree>(DecisionTree(U))});
            possible_trees.emplace_back(curr_family[U]);
          }
        }
      } else if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
        // Leaf-level candidates are added after leaf construction below (we need indices)
      } else {
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> curr_dims = curr_split.tree_index->split_dims;
          curr_dims.insert(curr_split.split_coordinate);
          curr_dims.insert(feature_dim);
          curr_dims.erase(0);
          if (possibleExists(feature_dim, possible_splits, curr_dims))
            continue;
          if (!this->dims_allowed(curr_dims))
            continue;
          if (auto found = treeExists(curr_dims, curr_family))
            possible_splits.emplace_back(feature_dim, found, (size_t)0);
          else {
            curr_family.insert({curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))});
            possible_splits.emplace_back(feature_dim, curr_family[curr_dims], (size_t)0);
          }
        }
      }

      // update values of individuals of split interval
      std::vector<double> update_s = curr_split.M_s, update_b = curr_split.M_b;
      switch (this->loss) {
      case LossType::L1:
      case LossType::L2:
      case LossType::median: {
        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
            samples_Y[individual] -= update_s;
          } else {
            samples_Y[individual] -= update_b;
          }
        }
        break;
      }
      case LossType::logit: {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::min(std::max(epsilon, M), 1 - epsilon); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::min(std::max(epsilon, M), 1 - epsilon); });

        double M_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double M_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p) {
          update_s[p] = log(M_s[p] / M_sp) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_bp) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
            weights[individual] += update_s;
          } else {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_2: {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p) {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
            weights[individual] += update_s;
          } else {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_3: {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });

        std::for_each(M_s.begin(), M_s.end(), [&](double &M) { M = log(M); });
        std::for_each(M_b.begin(), M_b.end(), [&](double &M) { M = log(M); });

        double M_sp = std::max(epsilon, curr_split.M_sp);
        double M_bp = std::max(epsilon, curr_split.M_bp);

        M_sp = log(M_sp);
        M_bp = log(M_bp);

        double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
        double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (unsigned int p = 0; p < M_s.size(); ++p) {
          update_s[p] = M_s[p] - sum_s - W_s_mean[p];
          update_b[p] = M_b[p] - sum_b - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
            weights[individual] += update_s;
          } else {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_4: {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p) {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals) {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
            weights[individual] += update_s;
          } else {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::exponential: {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) {
          S = std::min(std::max(epsilon, S), 1 - epsilon);
        });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) {
          S = std::min(std::max(epsilon, S), 1 - epsilon);
        });

        double sum_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double sum_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        for (unsigned int p = 0; p < sum_s.size(); ++p) {
          update_s[p] = log(sum_s[p] / sum_sp);
          update_b[p] = log(sum_b[p] / sum_bp);
        }

        for (int individual : curr_split.leaf_index->individuals) {
          for (unsigned int p = 0; p < update_s.size(); ++p) {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            } else {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_2: {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;
        std::vector<double> sum_s2 = curr_split.M_s;
        std::vector<double> sum_b2 = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = std::max(epsilon, S); });

        std::for_each(sum_s2.begin(), sum_s2.end(), [this](double &S) { S = std::max(epsilon, 1 - S); });
        std::for_each(sum_b2.begin(), sum_b2.end(), [this](double &S) { S = std::max(epsilon, 1 - S); });

        for (size_t p = 0; p < value_size; ++p) {
          update_s[p] = log(sum_s[p] / sum_s2[p]);
          update_b[p] = log(sum_b[p] / sum_b2[p]);
        }

        for (int individual : curr_split.leaf_index->individuals) {
          for (size_t p = 0; p < value_size; ++p) {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            } else {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_3: {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = std::max(epsilon, S); });

        std::for_each(sum_s.begin(), sum_s.end(), [&](double &S) { S = log(S); });
        std::for_each(sum_b.begin(), sum_b.end(), [&](double &S) { S = log(S); });

        double sum_sp = std::max(epsilon, curr_split.M_sp);
        double sum_bp = std::max(epsilon, curr_split.M_bp);

        sum_sp = log(sum_sp);
        sum_bp = log(sum_bp);

        sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
        sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);

        sum_sp = sum_sp / (sum_s.size() + 1);
        sum_bp = sum_bp / (sum_b.size() + 1);

        for (size_t p = 0; p < sum_s.size(); ++p) {
          update_s[p] = sum_s[p] - sum_sp;
          update_b[p] = sum_b[p] - sum_bp;
        }

        for (int individual : curr_split.leaf_index->individuals) {
          for (size_t p = 0; p < update_s.size(); ++p) {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            } else {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      }

      // construct new leaves
      Leaf leaf_s, leaf_b;
      {
        leaf_s.individuals = curr_split.I_s;
        leaf_b.individuals = curr_split.I_b;

        leaf_s.value = update_s;
        leaf_b.value = update_b;

        // initialize interval with split interval
        leaf_s.intervals = curr_split.leaf_index->intervals;
        leaf_b.intervals = curr_split.leaf_index->intervals;

        // interval of leaf with smaller individuals has new upper bound in splitting dimension
        leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point;
        // interval of leaf with bigger individuals has new lower bound in splitting dimension
        leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
      }

      // construct split_dims of resulting tree when splitting in split_coordinate
      std::set<int> resulting_dims = curr_split.tree_index->split_dims;
      resulting_dims.insert(curr_split.split_coordinate);
      resulting_dims.erase(0);

      // check if resulting tree already exists in family
      std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);

      // leaves / hist: add the (leaf, dimension) candidates of a new leaf, as in rpf.cpp
      auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree> &T, size_t li) {
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> res_dims = T->split_dims;
          res_dims.insert(feature_dim);
          res_dims.erase(0);
          if (!this->dims_allowed(res_dims))
            continue;
          if (!this->leafCandidateExists(possible_splits, T, li, feature_dim))
            possible_splits.emplace_back(feature_dim, T, li);
        }
      };

      // determine which tree is modified
      if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate)) &&
          delete_leaves) { // if split variable is already in tree to be split
        // change values
        {
          leaf_s.value += curr_split.leaf_index->value;
          leaf_b.value += curr_split.leaf_index->value;
        }
        // index before push_back, which may reallocate leaves and invalidate leaf_index
        size_t idx_b = (size_t)(curr_split.leaf_index - &curr_split.tree_index->leaves[0]);
        *curr_split.leaf_index = leaf_b;                 // replace old interval
        curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
        size_t idx_s = curr_split.tree_index->leaves.size() - 1;
        if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
          add_leaf_candidates(curr_split.tree_index, idx_b);
          add_leaf_candidates(curr_split.tree_index, idx_s);
        } else if (split_structure_mode_ == 1 || split_structure_mode_ == 2) {
          updateThresholdWeight(curr_split.tree_index, idx_b, false, samples_X);
          updateThresholdWeight(curr_split.tree_index, idx_s, true, samples_X);
        }
      } else { // otherwise
        if (!found_tree) {
          curr_family.insert({resulting_dims, std::make_shared<DecisionTree>(DecisionTree(resulting_dims))});
          found_tree = curr_family[resulting_dims];
        }
        found_tree->leaves.push_back(leaf_s); // append new leaves
        found_tree->leaves.push_back(leaf_b);
        size_t idx_s = found_tree->leaves.size() - 2;
        size_t idx_b = found_tree->leaves.size() - 1;
        if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
          add_leaf_candidates(found_tree, idx_s);
          add_leaf_candidates(found_tree, idx_b);
        } else if (split_structure_mode_ == 1 || split_structure_mode_ == 2) {
          updateThresholdWeight(found_tree, idx_s, true, samples_X);
          updateThresholdWeight(found_tree, idx_b, true, samples_X);
        }
      }
    }
  }

  fit_null_leaf(curr_family, samples_Y);

  // remove empty trees, and drop the rows and split-search caches of each leaf (needed only while growing)
  auto keys = getKeys(curr_family);
  for (auto &key : keys) {
    if (curr_family[key]->leaves.size() == 0) {
      curr_family.erase(key);
      continue;
    }
    for (auto &leaf : curr_family[key]->leaves) {
      leaf.individuals.clear();
      leaf.individuals.shrink_to_fit();
      leaf.free_cache();
    }
  }

  tree_families[n] = curr_family;
}

// fit forest to new data
void ClassificationRPF::fit()
{
  // Use the base class multithreaded trainer with RNG seeding identical to regression
  RandomPlantedForest::fit();
}

/*  retrospectively change parameters of existing class object,
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data',
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values)
{
  if (keys.size() != values.size()) {
    warn("Size of input vectors is not the same. ");
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i) {
    if (keys[i] == "deterministic") {
      this->deterministic = values[i];
    } else if (keys[i] == "nthreads") {
      this->nthreads = values[i];
    } else if (keys[i] == "purify") {
      this->purify_forest = values[i];
    } else if (keys[i] == "n_trees") {
      this->n_trees = values[i];
    } else if (keys[i] == "n_splits") {
      this->n_splits = values[i];
    } else if (keys[i] == "t_try") {
      this->t_try = values[i];
    } else if (keys[i] == "split_try") {
      this->split_try = values[i];
    } else if (keys[i] == "max_interaction") {
      this->max_interaction = values[i];
    } else if (keys[i] == "cv") {
      this->cross_validate = values[i];
    } else if (keys[i] == "loss") {
      if (keys[i] == "L1") {
        this->loss = LossType::L1;
        this->calcLoss = &ClassificationRPF::L1_loss;
      } else if (keys[i] == "L2") {
        this->loss = LossType::L2;
        this->calcLoss = &ClassificationRPF::L2_loss;
      } else if (keys[i] == "median") {
        this->loss = LossType::median;
        this->calcLoss = &ClassificationRPF::median_loss;
      } else if (keys[i] == "logit") {
        this->loss = LossType::logit;
        this->calcLoss = &ClassificationRPF::logit_loss;
      } else if (keys[i] == "exponential") {
        this->loss = LossType::exponential;
        this->calcLoss = &ClassificationRPF::exponential_loss;
      } else {
        warn("Unkown loss function.");
      }
    } else if (keys[i] == "delta") {
      this->delta = values[i];
    } else if (keys[i] == "epsilon") {
      this->epsilon = values[i];
    } else if (keys[i] == "split_decay_rate") {
      this->split_decay_rate_ = values[i];
    } else if (keys[i] == "max_candidates") {
      this->max_candidates_ = static_cast<size_t>(values[i]);
    } else {
      warn("Unkown parameter key  '" + keys[i] + "' .");
    }
  }
  this->fit();
}
