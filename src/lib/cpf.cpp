
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
void ClassificationRPF::set_loss_function(const String &loss)
{
  if (loss == "L1")
  {
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  }
  else if (loss == "L2")
  {
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  }
  else if (loss == "median")
  {
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  }
  else if (loss == "logit")
  {
    this->loss = LossType::logit;
    this->calcLoss = &ClassificationRPF::logit_loss;
  }
  else if (loss == "logit_2")
  {
    this->loss = LossType::logit_2;
    this->calcLoss = &ClassificationRPF::logit_loss_2;
  }
  else if (loss == "logit_3")
  {
    this->loss = LossType::logit_3;
    this->calcLoss = &ClassificationRPF::logit_loss_3;
  }
  else if (loss == "logit_4")
  {
    this->loss = LossType::logit_4;
    this->calcLoss = &ClassificationRPF::logit_loss_4;
  }
  else if (loss == "exponential")
  {
    this->loss = LossType::exponential;
    this->calcLoss = &ClassificationRPF::exponential_loss;
  }
  else if (loss == "exponential_2")
  {
    this->loss = LossType::exponential_2;
    this->calcLoss = &ClassificationRPF::exponential_loss_2;
  }
  else if (loss == "exponential_3")
  {
    this->loss = LossType::exponential_3;
    this->calcLoss = &ClassificationRPF::exponential_loss_3;
  }
  else
  {
    Rcout << "Unkown loss function, set to default (L2)." << std::endl;
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  }
}

// constructor with parameters split_try, t_try, purify_forest, deterministic, nthreads
ClassificationRPF::ClassificationRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                                     const String loss, const NumericVector parameters)
    : RandomPlantedForest(
        samples_Y,
        samples_X,
        // pass first 13 parameters to base (includes split_structure)
        parameters.size() >= 13 ? parameters[Rcpp::Range(0, 12)] : parameters[Rcpp::Range(0, 11)]
      )
{

  // Ensure correct Rcpp RNG state
  Rcpp::RNGScope scope;

  // initialize class members
  std::vector<double> pars = to_std_vec(parameters);
  set_loss_function(loss);
  if (pars.size() != 15)
  {
    Rcout << "Wrong number of parameters - set to default." << std::endl;
    this->max_interaction = 1;
    this->n_trees = 50;
    this->n_splits = 30;
    this->split_try = 10;
    this->t_try = 0.4;
    this->purify_forest = 0;
    this->deterministic = 0;
    this->nthreads = 1;
    this->cross_validate = 0;
    this->split_decay_rate_ = 0.1;
    this->max_candidates_   = 50;
    this->delete_leaves   = 1;
    this->delta = 0.1;
    this->epsilon = 0;
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
    this->delete_leaves   =  pars[11];
    // pars[12] is split_structure for base; already consumed by base
    this->delta = pars[13];
    this->epsilon = pars[14];   
  }

  // set data and data related members
  this->set_data(samples_Y, samples_X);
}

// Params-only constructor: parses configuration and loss but loads no data
// and does not fit. Used by rpf_unmarshal() to rebuild a serialized forest.
ClassificationRPF::ClassificationRPF(const String loss, const NumericVector parameters)
    : RandomPlantedForest(
          parameters.size() >= 13 ? NumericVector(parameters[Rcpp::Range(0, 12)])
                                  : NumericVector(parameters[Rcpp::Range(0, 11)]))
{
  std::vector<double> pars = to_std_vec(parameters);
  set_loss_function(loss);
  if (pars.size() == 15) {
    this->delta = pars[13];
    this->epsilon = pars[14];
  }
}

// Mode 1: cur_trees_2 (classification variant)
Split ClassificationRPF::calcOptimalSplit_curTrees2(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{

  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;
  std::set<int> tree_dims;
  std::vector<double> unique_samples;
  int k;
  unsigned int n = 0;
  double leaf_size;

  // sample possible splits
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));  
  std::vector<size_t> split_candidates;

      // 1) Build weights = exp(-decay_rate * age)
    std::vector<double> weights_vec(possible_splits.size());
    for (size_t i = 0; i < possible_splits.size(); ++i) {
      weights_vec[i] = std::exp(-split_decay_rate_ * possible_splits[i].age);
    }

    // 2) Sample n_candidates indices *without* replacement
    std::vector<size_t> sample_idxs;
    sample_idxs.reserve(n_candidates);

    if (!deterministic) {
      // Use weighted reservoir sampling driven by thread-local RNG
      std::vector<bool> used(possible_splits.size(), false);
      std::vector<double> w = weights_vec;
      while (sample_idxs.size() < n_candidates) {
        double tot = 0.0; for (double v : w) tot += (v > 0.0 ? v : 0.0);
        if (tot <= 0.0) break;
        double u = rpf_utils::rng_runif(0.0, tot);
        double acc = 0.0; size_t pick = 0;
        for (size_t i = 0; i < w.size(); ++i) { acc += (w[i] > 0.0 ? w[i] : 0.0); if (u <= acc) { pick = i; break; } }
        if (!used[pick]) { used[pick] = true; sample_idxs.push_back(pick); w[pick] = 0.0; }
      }
    } else {
      // deterministic fallback: first n_candidates
      for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i)
        sample_idxs.push_back(i);
    }

    split_candidates = sample_idxs;

  // track which one gave us the best split
  size_t chosen_idx = std::numeric_limits<size_t>::max();

  // consider a fraction of possible splits
  while (n < n_candidates)
  {

    // since size of possible splits changes, check if candidate in range
    if (possible_splits.empty())
      break;
    if (split_candidates[n] >= 0 && (size_t)split_candidates[n] >= possible_splits.size())
    { ++n; continue; }

    auto candidate = possible_splits.begin();
    std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
    k = candidate->dim - 1;                     // split dim of  candidate, converted to index starting at 0
    leaf_size = n_leaves[k];

    // Test if splitting in the  tree w.r.t. the coordinate "k" is an element of candidate tree
    tree_dims = candidate->tree->split_dims;
    tree_dims.erase(k + 1);
    tree_dims.erase(0);

    std::vector<std::shared_ptr<DecisionTree>> curr_trees;
    if (tree_dims.size() == 0)
      curr_trees.push_back(curr_family[std::set<int>{0}]);
    if (curr_family.find(tree_dims) != curr_family.end())
      curr_trees.push_back(curr_family[tree_dims]);
    if (curr_family.find(candidate->tree->split_dims) != curr_family.end())

      // go through all trees in current family
      for (auto &curr_tree : curr_trees)
      {

        // skip if tree has no leaves
        if (curr_tree->leaves.size() == 0)
          continue;

        // go through all leaves of current tree
     /*    for (auto &leaf : curr_tree->leaves)
        {

          std::vector<double> tot_sum(value_size, 0);

          // extract sample points according to individuals from X and Y
          unique_samples = std::vector<double>(leaf.individuals.size());
          for (unsigned int i = 0; i < leaf.individuals.size(); ++i)
          {
            unique_samples[i] = X[leaf.individuals[i]][k];
          }
          std::sort(unique_samples.begin(), unique_samples.end());
          unique_samples.erase(std::unique(unique_samples.begin(), unique_samples.end()), unique_samples.end());

          // check if number of sample points is within limit
          if (unique_samples.size() < 2 * leaf_size)
            continue;

          // consider split_try-number of samples
          std::vector<int> samples;
          if (deterministic)
          { // sequential samples if deterministic
            samples = std::vector<int>(std::min((int)unique_samples.size() - 1, 9));
            std::iota(samples.begin(), samples.end(), 1);
          }
          else
          { // randomly picked samples otherwise
            samples = std::vector<int>(split_try);
            for (size_t i = 0; i < samples.size(); ++i)
              samples[i] = rpf_utils::rng_randint((int)leaf_size, (int)unique_samples.size() - (int)leaf_size);
            std::sort(samples.begin(), samples.end());
          }

          // go through samples
          for (size_t sample_pos = 0; sample_pos < samples.size(); ++sample_pos)
          {

            // get samplepoint
            sample_point = unique_samples[samples[sample_pos]];

            // clear current split
            {
              curr_split.I_s.clear();
              curr_split.I_b.clear();
              curr_split.I_s.reserve(leaf.individuals.size());
              curr_split.I_b.reserve(leaf.individuals.size());
              curr_split.M_s = std::vector<double>(value_size, 0);
              curr_split.M_b = std::vector<double>(value_size, 0);
            }

            // get samples greater/smaller than samplepoint
            if (sample_pos == 0)
            {
              curr_split.sum_s = std::vector<double>(value_size, 0);
              curr_split.sum_b = std::vector<double>(value_size, 0);

              for (int individual : leaf.individuals)
              {
                if (X[individual][k] < sample_point)
                {
                  curr_split.I_s.push_back(individual);
                  curr_split.sum_s += Y[individual];
                }
                else
                {
                  curr_split.I_b.push_back(individual);
                  curr_split.sum_b += Y[individual];
                }
              }

              tot_sum = curr_split.sum_s + curr_split.sum_b;
            }
            else
            {

              for (int individual : leaf.individuals)
              {
                if (X[individual][k] < sample_point)
                {
                  if (X[individual][k] >= unique_samples[samples[sample_pos - 1]])
                  {
                    curr_split.sum_s += Y[individual];
                  }
                  curr_split.I_s.push_back(individual);
                }
                else
                {
                  curr_split.I_b.push_back(individual);
                }
              }

              curr_split.sum_b = tot_sum - curr_split.sum_s;
            }

            // accumulate squared mean and get mean
            (this->*ClassificationRPF::calcLoss)(curr_split);

            // update split if squared sum is smaller
            if (curr_split.min_sum < min_split.min_sum)
            {
              min_split = curr_split;
              min_split.tree_index = curr_tree;
              min_split.leaf_index = &leaf;
              min_split.split_coordinate = k + 1;
              min_split.split_point = sample_point;
              chosen_idx = split_candidates[n];
            }
          }
        } */

        // Mirror regression: traverse all leaves, sample split_try positions per leaf
        for (auto &leaf : curr_tree->leaves) {
          std::vector<size_t> order_cf; std::vector<double> sorted_vals_cf;
          ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
          std::vector<double> unique_vals = compute_unique_sorted_values(sorted_vals_cf);
          if (unique_vals.size() < 2 * static_cast<size_t>(leaf_size)) continue;

          const size_t m = leaf.individuals.size();
          std::vector<int> samples;
          if (this->deterministic) {
            int maxp = std::min<int>((int)unique_vals.size() - 1, 9);
            samples.resize(maxp); std::iota(samples.begin(), samples.end(), 1);
          } else {
            samples.resize(this->split_try);
            for (size_t i = 0; i < samples.size(); ++i)
              samples[i] = rpf_utils::rng_randint(leaf_size, (int)unique_vals.size() - (int)leaf_size);
            std::sort(samples.begin(), samples.end());
          }

          for (size_t si = 0; si < samples.size(); ++si) {
            const double sp = unique_vals[(size_t)samples[si]];
            size_t pos = (size_t)(std::lower_bound(sorted_vals_cf.begin(), sorted_vals_cf.end(), sp) - sorted_vals_cf.begin());
            if (pos == 0 || pos >= m) continue;
            if (pos < (size_t)leaf_size || (m - pos) < (size_t)leaf_size) continue;

            // Build I_s/I_b and sums for classification loss
            curr_split.I_s.clear(); curr_split.I_b.clear();
            curr_split.I_s.reserve(m); curr_split.I_b.reserve(m);
            curr_split.sum_s.assign(value_size, 0.0); curr_split.sum_b.assign(value_size, 0.0);
            for (int ind : leaf.individuals) {
              if (X[ind][k] < sp) { curr_split.I_s.push_back(ind); curr_split.sum_s += Y[ind]; }
              else { curr_split.I_b.push_back(ind); curr_split.sum_b += Y[ind]; }
            }

            (this->*ClassificationRPF::calcLoss)(curr_split);
            if (curr_split.min_sum < min_split.min_sum) {
              min_split = curr_split;
              min_split.tree_index       = curr_tree;
              min_split.leaf_index       = &leaf;
              min_split.split_coordinate = k + 1;
              min_split.split_point      = sp;
              chosen_idx                 = split_candidates[n];
            }
          }
        }

      }

    ++n;
  }

  for (size_t idx : split_candidates) {
    if (idx == chosen_idx) {
     possible_splits[idx].age = 0.0;        // reset for the winner
    } else {
     possible_splits[idx].age += 1.0;      // age everyone else
    }
  }

  return min_split;
}

// Mode 3: leaves (classification variant)
Split ClassificationRPF::calcOptimalSplit_leaves(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  Split curr_split, min_split; min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y; curr_split.W = &weights;
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights_vec(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights_vec[i] = std::exp(-split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!deterministic) {
    std::vector<bool> used(possible_splits.size(), false);
    std::vector<double> w = weights_vec;
    while (sample_idxs.size() < n_candidates) { double tot = 0.0; for (double v:w) tot += (v>0.0? v:0.0); if (tot<=0.0) break; double u=rpf_utils::rng_runif(0.0, tot); double acc=0.0; size_t pick=0; for (size_t i=0;i<w.size();++i){ acc += (w[i]>0.0? w[i]:0.0); if (u<=acc){ pick=i; break; } } if (!used[pick]) { used[pick]=true; sample_idxs.push_back(pick); w[pick]=0.0; } }
  } else {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i) sample_idxs.push_back(i);
  }
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx);
    int k = it->dim - 1; int leaf_size = n_leaves[k];
    auto treePtr = it->tree; if (treePtr->leaves.empty() || it->leaf_idx >= treePtr->leaves.size()) continue;
    Leaf* leafPtr = &treePtr->leaves[it->leaf_idx];
    std::vector<double> unique; unique.reserve(leafPtr->individuals.size());
    for (int ind : leafPtr->individuals) unique.push_back(X[ind][k]);
    std::sort(unique.begin(), unique.end()); unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    int left = (int)leaf_size; int right = (int)unique.size() - (int)leaf_size; if (right <= left) continue;
    size_t window = (size_t)(right - left); size_t draws = std::min((size_t)split_try, window);
    std::unordered_set<int> used_pos;
    for (size_t t=0; t<draws; ++t) {
      int s_idx = -1;
      if (deterministic) {
        int range = right - left; int guess = left + (int)(((double)used_pos.size()+0.5)/((double)range+0.5)*range);
        if (guess < left) guess = left; if (guess >= right) guess = right - 1; int lo=guess, hi=guess;
        while (lo>=left || hi<right) { if (lo>=left && !used_pos.count(lo)) { s_idx=lo; break; } if (hi<right && !used_pos.count(hi)) { s_idx=hi; break; } --lo; ++hi; }
        if (s_idx == -1) for (int p=left;p<right;++p) if (!used_pos.count(p)) { s_idx=p; break; }
      } else { do { s_idx = rpf_utils::rng_randint(left, right); } while (used_pos.count(s_idx)); }
      used_pos.insert(s_idx);
      double sp = unique[(size_t)s_idx];
      curr_split.I_s.clear(); curr_split.I_b.clear();
      curr_split.M_s.assign(value_size, 0); curr_split.M_b.assign(value_size, 0);
      curr_split.sum_s.assign(value_size, 0); curr_split.sum_b.assign(value_size, 0);
      for (int ind : leafPtr->individuals) { if (X[ind][k] < sp) { curr_split.I_s.push_back(ind); curr_split.sum_s += Y[ind]; } else { curr_split.I_b.push_back(ind); curr_split.sum_b += Y[ind]; } }
      (this->*ClassificationRPF::calcLoss)(curr_split);
      if (curr_split.min_sum < min_split.min_sum) { min_split = curr_split; min_split.tree_index = treePtr; min_split.leaf_index = leafPtr; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; }
    }
  }
  for (size_t idx : sample_idxs) { if ((int)idx != best_idx) possible_splits[idx].age += 1.0; else possible_splits[idx].age = 0.0; }
  return min_split;
}

// Mode 4: histogram-binned (classification variant)
Split ClassificationRPF::calcOptimalSplit_hist(const std::vector<std::vector<double>> &Y,
                                           const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits,
                                           TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  Split min_split; min_split.min_sum = std::numeric_limits<double>::infinity();
  if (possible_splits.empty()) return min_split;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights_vec(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights_vec[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!deterministic) {
    std::vector<bool> used(possible_splits.size(), false);
    std::vector<double> w = weights_vec;
    while (sample_idxs.size() < n_candidates) { double tot = 0.0; for (double v:w) tot += (v>0.0? v:0.0); if (tot<=0.0) break; double u=rpf_utils::rng_runif(0.0, tot); double acc=0.0; size_t pick=0; for (size_t i=0;i<w.size();++i){ acc += (w[i]>0.0? w[i]:0.0); if (u<=acc){ pick=i; break; } } if (!used[pick]) { used[pick]=true; sample_idxs.push_back(pick); w[pick]=0.0; } }
  } else {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i) sample_idxs.push_back(i);
  }

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

    // Build histogram for this leaf and feature k using global cut points from base
    const auto &cuts_k = (k >= 0 && k < (int)feature_cut_points_.size()) ? feature_cut_points_[k] : std::vector<double>{};
    size_t Kf = cuts_k.size() + 1; if (Kf < 2) continue;
    std::vector<int> cnt(Kf, 0);
    std::vector<std::vector<double>> sum(Kf, std::vector<double>(this->value_size, 0.0));
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

    // Single sweep over bin boundaries
    const int total_n = (int)m;
    std::vector<double> total_sum(this->value_size, 0.0);
    for (size_t b = 0; b < Kf; ++b) for (size_t p = 0; p < this->value_size; ++p) total_sum[p] += sum[b][p];
    int left_n = 0; std::vector<double> left_sum(this->value_size, 0.0);
    for (size_t b_left = 0; b_left + 1 <= Kf - 1; ++b_left) {
      left_n += cnt[b_left];
      for (size_t p = 0; p < this->value_size; ++p) left_sum[p] += sum[b_left][p];
      int right_n = total_n - left_n;
      if (left_n < leaf_min || right_n < leaf_min) continue;

      // Fill curr split buffers for loss calculation
      Split curr_split; curr_split.Y = &Y; curr_split.W = &weights;
      curr_split.I_s.clear(); curr_split.I_b.clear();
      curr_split.sum_s.assign(this->value_size, 0.0);
      curr_split.sum_b.assign(this->value_size, 0.0);
      for (size_t p = 0; p < this->value_size; ++p) {
        curr_split.sum_s[p] = left_sum[p];
        curr_split.sum_b[p] = total_sum[p] - left_sum[p];
      }

      // For classification losses we still need I_s/I_b indices; build once per boundary
      curr_split.I_s.reserve(m); curr_split.I_b.reserve(m);
      double sp_val;
      if (k >= 0 && k < (int)feature_cut_points_.size() && !feature_cut_points_[k].empty()) {
        const auto &cuts = feature_cut_points_[k];
        size_t cp_idx = (size_t)std::min<size_t>(b_left, cuts.size() - 1);
        sp_val = cuts[cp_idx];
      } else {
        sp_val = 0.5 * (leafPtr->intervals[k].first + leafPtr->intervals[k].second);
      }
      for (int ind : leafPtr->individuals) {
        if (X[ind][k] < sp_val) curr_split.I_s.push_back(ind); else curr_split.I_b.push_back(ind);
      }

      // Compute classification loss
      (this->*ClassificationRPF::calcLoss)(curr_split);

      if (curr_split.min_sum < min_split.min_sum) {
        min_split = curr_split;
        min_split.tree_index = it->tree;
        min_split.leaf_index = leafPtr;
        min_split.split_coordinate = k + 1;
        min_split.split_point = sp_val;
        best_idx = (int)idx;
      }
    }
  }

  for (size_t idx : sample_idxs) { if ((int)idx != best_idx) possible_splits[idx].age += 1.0; else possible_splits[idx].age = 0.0; }
  return min_split;
}

// Mode 2: cur_trees_1 (classification variant)
Split ClassificationRPF::calcOptimalSplit_curTrees1(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  // reuse current implementation by sampling per-leaf candidates across predecessor/current trees
  // We delegate to the old flow by temporarily constructing the same sampling but using loss with W
  // For brevity, call the curTrees2 variant which already samples leaves within available trees
  return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family, weights);
}

// Mode 0: res_trees (classification variant)
Split ClassificationRPF::calcOptimalSplit_resTrees(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<RandomPlantedForest::ResultingTreeCandidate> &possible_trees, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  // Classification loss evaluation on res_trees follows the base structure; to keep changes minimal here,
  // we adopt the cur_trees_1 sampling over the trees in possible_trees' dims using our calcLoss and W.
  // Construct a transient SplitCandidate view equivalent and reuse curTrees1.
  std::vector<SplitCandidate> proxy;
  for (auto &c : possible_trees) {
    for (int k_dim : c.tree->split_dims) {
      proxy.emplace_back(k_dim, c.tree, (size_t)0);
    }
  }
  return this->calcOptimalSplit_curTrees1(Y, X, proxy, curr_family, weights);
}

// Dispatcher selecting by split_structure_mode_
Split ClassificationRPF::calcOptimalSplit(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  if (split_structure_mode_ == 4) return this->calcOptimalSplit_hist(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 3) return this->calcOptimalSplit_leaves(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 2) return this->calcOptimalSplit_curTrees1(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 1) return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family, weights);
  return Split{};
}

void ClassificationRPF::create_tree_family(std::vector<Leaf> initial_leaves, size_t n)
{

  TreeFamily curr_family;
  curr_family.insert(std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))); // save tree with one leaf in the beginning

  // Seed per mode
  std::vector<SplitCandidate> possible_splits;
  std::vector<ResultingTreeCandidate> possible_trees;
  if (split_structure_mode_ == 0) {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_trees.emplace_back(treePtr);
    }
  } else if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
    auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
        std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
        if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
        if (!this->leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
      }
    };
    auto null_tree = curr_family[{0}];
    if (!null_tree->leaves.empty()) add_leaf_candidates(null_tree, 0);
  } else {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_splits.emplace_back(feature_dim, treePtr, (size_t)0);
    }
  }

  // sample data points with replacement
  int sample_index;
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;

  // deterministic
  if (deterministic)
  {
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  }
  else
  {
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);

    // bagging/subsampling
    for (size_t i = 0; i < sample_size; ++i)
    {
      sample_index = rpf_utils::rng_randint(0, (int)sample_size);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }

  // initialize weights
  std::vector<std::vector<double>> weights;
  switch (this->loss)
  {
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
  for (int split_count = 0; split_count < n_splits; ++split_count)
  {

    // find optimal split
    if (split_structure_mode_ == 0) curr_split = this->calcOptimalSplit_resTrees(samples_Y, samples_X, possible_trees, curr_family, weights);
    else curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family, weights);

    // continue only if we get a significant result
    if (!std::isinf(curr_split.min_sum))
    {

      // update pools by mode
      if (split_structure_mode_ == 0) {
        std::set<int> Dprime = curr_split.tree_index->split_dims; Dprime.insert(curr_split.split_coordinate); Dprime.erase(0);
        if (!this->resultingTreeExists(possible_trees, Dprime)) { if (auto found = treeExists(Dprime, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({Dprime, std::make_shared<DecisionTree>(DecisionTree(Dprime))}); possible_trees.emplace_back(curr_family[Dprime]); } }
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> U = Dprime; U.insert(feature_dim); if (U.size() == Dprime.size()) continue; if (max_interaction >= 0 && U.size() > (size_t)max_interaction) continue; if (this->resultingTreeExists(possible_trees, U)) continue; if (auto found = treeExists(U, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({U, std::make_shared<DecisionTree>(DecisionTree(U))}); possible_trees.emplace_back(curr_family[U]); }
        }
      } else if (split_structure_mode_ == 3 || split_structure_mode_ == 4) {
        // Leaf-level candidates are added after leaf construction below (we need indices)
      } else {
        if (curr_split.tree_index->split_dims.count(curr_split.split_coordinate) == 0) {
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> curr_dims = curr_split.tree_index->split_dims; curr_dims.insert(curr_split.split_coordinate); curr_dims.insert(feature_dim); curr_dims.erase(0);
            if (possibleExists(feature_dim, possible_splits, curr_dims)) continue;
            if (max_interaction >= 0 && curr_dims.size() > (size_t)max_interaction) continue;
            if (auto found = treeExists(curr_dims, curr_family)) possible_splits.emplace_back(feature_dim, found, (size_t)0);
            else { curr_family.insert({curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))}); possible_splits.emplace_back(feature_dim, curr_family[curr_dims], (size_t)0); }
          }
        }
      }

      // update values of individuals of split interval
      std::vector<double> update_s = curr_split.M_s, update_b = curr_split.M_b;
      switch (this->loss)
      {
      case LossType::L1:
      case LossType::L2:
      case LossType::median:
      {
        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            samples_Y[individual] -= update_s;
          }
          else
          {
            samples_Y[individual] -= update_b;
          }
        }
        break;
      }
      case LossType::logit:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::min(std::max(epsilon, M), 1 - epsilon); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::min(std::max(epsilon, M), 1 - epsilon); });

        double M_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double M_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_sp) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_bp) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_2:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_3:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s.begin(), M_s.end(), [&](double &M)
                      { M = log(M); });
        std::for_each(M_b.begin(), M_b.end(), [&](double &M)
                      { M = log(M); });

        double M_sp = std::max(epsilon, curr_split.M_sp);
        double M_bp = std::max(epsilon, curr_split.M_bp);

        M_sp = log(M_sp);
        M_bp = log(M_bp);

        double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
        double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (unsigned int p = 0; p < M_s.size(); ++p)
        {
          update_s[p] = M_s[p] - sum_s - W_s_mean[p];
          update_b[p] = M_b[p] - sum_b - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_4:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::exponential:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::min(std::max(epsilon, S), 1 - epsilon); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::min(std::max(epsilon, S), 1 - epsilon); });

        double sum_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double sum_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        for (unsigned int p = 0; p < sum_s.size(); ++p)
        {
          update_s[p] = log(sum_s[p] / sum_sp);
          update_b[p] = log(sum_b[p] / sum_bp);
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (unsigned int p = 0; p < update_s.size(); ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_2:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;
        std::vector<double> sum_s2 = curr_split.M_s;
        std::vector<double> sum_b2 = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::max(epsilon, S); });

        std::for_each(sum_s2.begin(), sum_s2.end(), [this](double &S)
                      { S = std::max(epsilon, 1 - S); });
        std::for_each(sum_b2.begin(), sum_b2.end(), [this](double &S)
                      { S = std::max(epsilon, 1 - S); });

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(sum_s[p] / sum_s2[p]);
          update_b[p] = log(sum_b[p] / sum_b2[p]);
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (size_t p = 0; p < value_size; ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_3:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::max(epsilon, S); });

        std::for_each(sum_s.begin(), sum_s.end(), [&](double &S)
                      { S = log(S); });
        std::for_each(sum_b.begin(), sum_b.end(), [&](double &S)
                      { S = log(S); });

        double sum_sp = std::max(epsilon, curr_split.M_sp);
        double sum_bp = std::max(epsilon, curr_split.M_bp);

        sum_sp = log(sum_sp);
        sum_bp = log(sum_bp);

        sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
        sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);

        sum_sp = sum_sp / (sum_s.size() + 1);
        sum_bp = sum_bp / (sum_b.size() + 1);

        for (size_t p = 0; p < sum_s.size(); ++p)
        {
          update_s[p] = sum_s[p] - sum_sp;
          update_b[p] = sum_b[p] - sum_bp;
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (size_t p = 0; p < update_s.size(); ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
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

      // determine which tree is modified
      if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate))&& delete_leaves)
      { // if split variable is already in tree to be split
        // change values
        {
          leaf_s.value += curr_split.leaf_index->value;
          leaf_b.value += curr_split.leaf_index->value;
        }
        *curr_split.leaf_index = leaf_b;                 // replace old interval
        curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
        if (split_structure_mode_ == 3) {
          size_t idx_b = (size_t)(curr_split.leaf_index - &curr_split.tree_index->leaves[0]);
          size_t idx_s = curr_split.tree_index->leaves.size() - 1;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> res_dims_b = curr_split.tree_index->split_dims; res_dims_b.insert(feature_dim); res_dims_b.erase(0);
            if (max_interaction < 0 || res_dims_b.size() <= (size_t)max_interaction) {
              if (!this->leafCandidateExists(possible_splits, curr_split.tree_index, idx_b, feature_dim)) {
                possible_splits.emplace_back(feature_dim, curr_split.tree_index, idx_b);
              }
            }
            std::set<int> res_dims_s = curr_split.tree_index->split_dims; res_dims_s.insert(feature_dim); res_dims_s.erase(0);
            if (max_interaction < 0 || res_dims_s.size() <= (size_t)max_interaction) {
              if (!this->leafCandidateExists(possible_splits, curr_split.tree_index, idx_s, feature_dim)) {
                possible_splits.emplace_back(feature_dim, curr_split.tree_index, idx_s);
              }
            }
          }
        }
      }
      else
      {                                       // otherwise
        if (!found_tree) {
          curr_family.insert({resulting_dims, std::make_shared<DecisionTree>(DecisionTree(resulting_dims))});
          found_tree = curr_family[resulting_dims];
        }
        found_tree->leaves.push_back(leaf_s); // append new leaves
        found_tree->leaves.push_back(leaf_b);
        if (split_structure_mode_ == 3) {
          size_t idx_s = found_tree->leaves.size() - 2; size_t idx_b = found_tree->leaves.size() - 1;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> res_dims_s = found_tree->split_dims; res_dims_s.insert(feature_dim); res_dims_s.erase(0);
            if (max_interaction < 0 || res_dims_s.size() <= (size_t)max_interaction) {
              if (!this->leafCandidateExists(possible_splits, found_tree, idx_s, feature_dim)) {
                possible_splits.emplace_back(feature_dim, found_tree, idx_s);
              }
            }
            std::set<int> res_dims_b = found_tree->split_dims; res_dims_b.insert(feature_dim); res_dims_b.erase(0);
            if (max_interaction < 0 || res_dims_b.size() <= (size_t)max_interaction) {
              if (!this->leafCandidateExists(possible_splits, found_tree, idx_b, feature_dim)) {
                possible_splits.emplace_back(feature_dim, found_tree, idx_b);
              }
            }
          }
        }
      }
    }
  }

  // remove empty trees & clear individuals of each tree
  auto keys = getKeys(curr_family);
  for (auto &key : keys)
  {
    if (curr_family[key]->leaves.size() == 0)
    {
      curr_family.erase(key);
      continue;
    }
    for (auto &leaf : curr_family[key]->leaves)
    {
      leaf.individuals.clear();
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
void ClassificationRPF::set_parameters(StringVector keys, NumericVector values)
{
  if (keys.size() != values.size())
  {
    Rcout << "Size of input vectors is not the same. " << std::endl;
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    if (keys[i] == "deterministic")
    {
      this->deterministic = values[i];
    }
    else if (keys[i] == "nthreads")
    {
      this->nthreads = values[i];
    }
    else if (keys[i] == "purify")
    {
      this->purify_forest = values[i];
    }
    else if (keys[i] == "n_trees")
    {
      this->n_trees = values[i];
    }
    else if (keys[i] == "n_splits")
    {
      this->n_splits = values[i];
    }
    else if (keys[i] == "t_try")
    {
      this->t_try = values[i];
    }
    else if (keys[i] == "split_try")
    {
      this->split_try = values[i];
    }
    else if (keys[i] == "max_interaction")
    {
      this->max_interaction = values[i];
    }
    else if (keys[i] == "cv")
    {
      this->cross_validate = values[i];
    }
    else if (keys[i] == "loss")
    {
      if (keys[i] == "L1")
      {
        this->loss = LossType::L1;
        this->calcLoss = &ClassificationRPF::L1_loss;
      }
      else if (keys[i] == "L2")
      {
        this->loss = LossType::L2;
        this->calcLoss = &ClassificationRPF::L2_loss;
      }
      else if (keys[i] == "median")
      {
        this->loss = LossType::median;
        this->calcLoss = &ClassificationRPF::median_loss;
      }
      else if (keys[i] == "logit")
      {
        this->loss = LossType::logit;
        this->calcLoss = &ClassificationRPF::logit_loss;
      }
      else if (keys[i] == "exponential")
      {
        this->loss = LossType::exponential;
        this->calcLoss = &ClassificationRPF::exponential_loss;
      }
      else
      {
        Rcout << "Unkown loss function." << std::endl;
      }
    }
    else if (keys[i] == "delta")
    {
      this->delta = values[i];
    }
    else if (keys[i] == "epsilon")
    {
      this->epsilon = values[i];
    }
    else if (keys[i] == "split_decay_rate") 
    {
      this->split_decay_rate_ = values[i];
    }
    else if (keys[i] == "max_candidates") 
    {
      this->max_candidates_ = static_cast<size_t>(values[i]);
    }
    else
    {
      Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}
