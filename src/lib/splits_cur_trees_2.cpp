// Split-mode: cur_trees_2. Tries random thresholds across all leaves of
// predecessor/current trees, using age decay for candidate sampling.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

Split RandomPlantedForest::calcOptimalSplit_curTrees2(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family,
    FamilyLoss *loss_fn
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));

  std::vector<double> weights(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i)
    weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs;
  sample_idxs.reserve(n_candidates);
  if (!this->deterministic) {
    std::vector<size_t> pos_idx;
    pos_idx.reserve(possible_splits.size());
    std::vector<double> pos_w;
    pos_w.reserve(possible_splits.size());
    for (size_t i = 0; i < weights.size(); ++i)
      if (weights[i] > 0.0) {
        pos_idx.push_back(i);
        pos_w.push_back(weights[i]);
      }
    const size_t P = pos_idx.size();
    if (P == 0) {
      std::vector<size_t> all(possible_splits.size());
      std::iota(all.begin(), all.end(), 0);
      size_t k = std::min<size_t>(n_candidates, all.size());
      for (size_t i = 0; i < k; ++i) {
        size_t j = i + static_cast<size_t>(rng_runif01() * (double)(all.size() - i));
        if (j >= all.size())
          j = all.size() - 1;
        std::swap(all[i], all[j]);
      }
      for (size_t i = 0; i < k; ++i)
        sample_idxs.push_back(all[i]);
    } else if (n_candidates * 8 < P) {
      size_t k2 = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double, size_t>> keys;
      keys.reserve(P);
      for (size_t i = 0; i < P; ++i) {
        double u = rng_runif01();
        if (u <= 0.0)
          u = std::numeric_limits<double>::min();
        double key = -std::log(u) / pos_w[i];
        keys.emplace_back(key, pos_idx[i]);
      }
      if (k2 < keys.size()) {
        std::nth_element(keys.begin(), keys.begin() + k2, keys.end(), [](const auto &a, const auto &b) {
          return a.first < b.first;
        });
        keys.resize(k2);
      }
      for (auto &kv : keys)
        sample_idxs.push_back(kv.second);
    } else {
      size_t k = std::min<size_t>(n_candidates, P);
      std::vector<std::pair<double, size_t>> keys;
      keys.reserve(P);
      for (size_t i = 0; i < P; ++i) {
        double u = rng_runif01();
        if (u <= 0.0)
          u = std::numeric_limits<double>::min();
        double key = -std::log(u) / pos_w[i];
        keys.emplace_back(key, pos_idx[i]);
      }
      if (k < keys.size()) {
        std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto &a, const auto &b) {
          return a.first < b.first;
        });
        keys.resize(k);
      }
      for (auto &kv : keys)
        sample_idxs.push_back(kv.second);
    }
  } else {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i)
      sample_idxs.push_back(i);
  }

  const int min_rows = this->min_leaf_rows_;
  const bool rows_check = min_rows > 1;
  int best_idx = -1;
  const size_t scan_T = loss_fn ? loss_fn->scan_terms() : 0;
  std::vector<double> scan_prefix, scan_s, scan_b;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    int leaf_size = this->n_leaves[k];

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
        // Reuse cached order and sorted values
        std::vector<size_t> order_cf;
        std::vector<double> sorted_vals_cf;
        ensure_order_and_sorted_vals_for_leaf(X, leaf, k, order_cf, sorted_vals_cf);
        // Unique count & values caching
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
        // the thresholds drawn from: [leaf_size, distinct - leaf_size], with min_leaf_rows only those with min_rows
        // rows on each side (a leaf without one is skipped and draws nothing)
        int left = leaf_size;
        int right_exclusive = (int)unique_count - leaf_size + 1;
        if (rows_check) {
          min_rows_range_sorted(sorted_vals_cf, *unique_ptr, leaf_size, min_rows, left, right_exclusive);
          if (right_exclusive <= left)
            continue;
        }
        std::vector<int> samples;
        if (this->deterministic) {
          // split_try thresholds evenly spread over the same range the random draws use
          samples = compute_even_spread_indices(left, right_exclusive, (size_t)this->split_try);
          if (samples.empty())
            continue;
        } else {
          samples.resize(this->split_try);
          for (size_t i = 0; i < samples.size(); ++i)
            samples[i] = rng_randint(left, right_exclusive);
          std::sort(samples.begin(), samples.end());
        }
        const bool single_eval = (samples.size() == 1);
        // The loss's running sums (FamilyLoss::scan_terms): the prefix sums of its row terms along the leaf's sorted
        // rows [i * T + q] (one threshold: the two sides' sums by a pass, as for Y)
        if (scan_T > 0) {
          scan_prefix.resize(single_eval ? scan_T : m * scan_T);
          scan_s.assign(scan_T, 0.0);
          scan_b.assign(scan_T, 0.0);
          size_t ns_single = 0;
          if (single_eval) {
            const double sp = (*unique_ptr)[samples[0]];
            for (int ind : leaf.individuals) {
              if (X[ind][k] < sp) {
                loss_fn->add_terms(ind, scan_s.data());
                ++ns_single;
              } else {
                loss_fn->add_terms(ind, scan_b.data());
              }
            }
            for (size_t q = 0; q < scan_T; ++q)
              scan_prefix[q] = scan_s[q] + scan_b[q];
          } else {
            loss_fn->prefix_terms(leaf, order_cf.data(), m, scan_prefix.data());
          }
          const double *tot = single_eval ? scan_prefix.data() : &scan_prefix[(m - 1) * scan_T];
          const double vp = loss_fn->region_value(tot, m);
          for (size_t si = 0; si < samples.size(); ++si) {
            const double sp = (*unique_ptr)[samples[si]];
            size_t pos = static_cast<size_t>(
                std::lower_bound(sorted_vals_cf.begin(), sorted_vals_cf.end(), sp) - sorted_vals_cf.begin()
            );
            if (pos == 0 || pos >= m)
              continue;
            if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
              continue;
            if (!single_eval) {
              const double *s = &scan_prefix[(pos - 1) * scan_T];
              std::copy(s, s + scan_T, scan_s.begin());
              for (size_t q = 0; q < scan_T; ++q)
                scan_b[q] = tot[q] - s[q];
            } else if (ns_single != pos) {
              continue;
            }
            const double loss = loss_fn->scan_score(scan_s.data(), scan_b.data(), pos, m - pos, vp);
            if (loss < min_split.min_sum) {
              min_split.min_sum = loss;
              min_split.tree_index = curr_tree;
              min_split.leaf_index = &leaf;
              min_split.split_coordinate = k + 1;
              min_split.split_point = sp;
              best_idx = (int)idx;
            }
          }
          continue;
        }
        std::vector<std::vector<double>> prefix_cf; // [value_size][m]
        std::vector<double> total_cf;               // [value_size]
        if (!single_eval)
          build_prefix_and_total_given_order(Y, leaf, order_cf, this->value_size, prefix_cf, total_cf);

        for (size_t si = 0; si < samples.size(); ++si) {
          const double sp = (*unique_ptr)[samples[si]];
          size_t pos = static_cast<size_t>(
              std::lower_bound(sorted_vals_cf.begin(), sorted_vals_cf.end(), sp) - sorted_vals_cf.begin()
          );
          if (pos == 0 || pos >= m)
            continue;
          if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
            continue;
          double loss = 0.0;
          if (loss_fn) { // the loss's criterion (classification), from the same sums
            curr_split.sum_s.assign(this->value_size, 0.0);
            curr_split.sum_b.assign(this->value_size, 0.0);
            if (!single_eval) {
              for (size_t p = 0; p < this->value_size; ++p) {
                curr_split.sum_s[p] = prefix_cf[p][pos - 1];
                curr_split.sum_b[p] = total_cf[p] - curr_split.sum_s[p];
              }
            } else {
              for (int ind : leaf.individuals) {
                if (X[ind][k] < sp) {
                  for (size_t p = 0; p < this->value_size; ++p)
                    curr_split.sum_s[p] += Y[ind][p];
                } else {
                  for (size_t p = 0; p < this->value_size; ++p)
                    curr_split.sum_b[p] += Y[ind][p];
                }
              }
            }
            loss = loss_fn->score(curr_split, leaf, k, sp, X);
          } else if (!single_eval) {
            for (size_t p = 0; p < this->value_size; ++p) {
              const double sum_s_base = prefix_cf[p][pos - 1];
              const double sum_b_base = total_cf[p] - sum_s_base;
              loss -= (sum_s_base * sum_s_base) / static_cast<double>(pos);
              loss -= (sum_b_base * sum_b_base) / static_cast<double>(m - pos);
            }
          } else {
            size_t ns = 0, nb = 0;
            std::vector<double> sum_s_adj(this->value_size, 0.0), sum_b_adj(this->value_size, 0.0);
            for (int ind : leaf.individuals) {
              const bool left_side = (X[ind][k] < sp);
              if (left_side) {
                ++ns;
                for (size_t p = 0; p < this->value_size; ++p) {
                  double v = Y[ind][p];
                  sum_s_adj[p] += v;
                }
              } else {
                ++nb;
                for (size_t p = 0; p < this->value_size; ++p) {
                  double v = Y[ind][p];
                  sum_b_adj[p] += v;
                }
              }
            }
            if (ns == 0 || nb == 0) {
              continue;
            }
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
              for (int ind : leaf.individuals) {
                if (X[ind][k] < sp) {
                  for (size_t p = 0; p < this->value_size; ++p) {
                    double v = Y[ind][p];
                    min_split.sum_s[p] += v;
                  }
                } else {
                  for (size_t p = 0; p < this->value_size; ++p) {
                    double v = Y[ind][p];
                    min_split.sum_b[p] += v;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  if (loss_fn)
    loss_fn->finalize(min_split, X);
  else
    finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}
