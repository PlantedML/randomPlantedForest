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
#include "internal_utils.hpp"

// Use utilities via namespace alias
using namespace rpf_utils;

namespace {
  // Thread-local cache for histogram mode per working set (per tree-family build)
  // Avoids races on the class member when building families in parallel
  thread_local std::vector<std::vector<int>> tls_working_bin_id;
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

// Mode 3: leaves implementation moved to lib/splits_leaves.cpp

// Mode 1: cur_trees_2 moved to lib/splits_cur_trees_2.cpp

// Mode 2: cur_trees_1 (pair-sampling within predecessor/current trees)
// Mode 2: cur_trees_1 moved to lib/splits_cur_trees_1.cpp

// Mode 0: res_trees (operate on resulting trees pool)
bool RandomPlantedForest::resultingTreeExists(const std::vector<RandomPlantedForest::ResultingTreeCandidate>& pool, const std::set<int>& dims) {
  for (const auto &c : pool) if (c.tree->get_split_dims() == dims) return true; return false;
}

// Mode 0: res_trees moved to lib/splits_res_trees.cpp

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

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
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
// fit() moved to lib/training.cpp

// predict single feature vector (from leaves variant)
// predict_single moved to lib/predict.cpp

// predict_matrix moved to lib/predict.cpp
// predict_vector moved to lib/predict.cpp

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



// purify_1 moved to lib/purify.cpp

// purify_2 moved to lib/purify.cpp

// purify_3 moved to lib/purify.cpp
