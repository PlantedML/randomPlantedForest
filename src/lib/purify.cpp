#include "rpf.hpp"
#include "kdtree.hpp"
#include <functional>
#include <numeric>
#include <thread>
#include <exception>
#include <limits>
#include <cmath>
#include <unordered_map>

// Mode 2 for one family: computes the purified components exactly (KD-tree range counts, one per leaf and subset)
// and adds them to `acc` as weighted rectangles grouped by (S, j); see rects.hpp. Nothing is rasterised and the
// family is not modified. Components above maxp_interaction are not produced (they are zero).
void RandomPlantedForest::purify_fast_exact_family(
    const TreeFamily &curr_family,
    int maxp_interaction,
    const rpf_kd::KDTree &kdt,
    const std::vector<std::vector<double>> &cuts,
    rpf_rect::Accumulator &acc
) const
{
  // Normalize cap: treat 0 (or out-of-range) as full order p = feature_size
  if (maxp_interaction <= 0 || maxp_interaction > feature_size)
    maxp_interaction = feature_size;

  auto popcount32 = [](unsigned int x) -> int {
    x = x - ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    return (int)((((x + (x >> 4)) & 0x0F0F0F0Fu) * 0x01010101u) >> 24);
  };

  // 1) Build lim_list (unique cut endpoints per feature)
  std::vector<std::vector<double>> lim_list(feature_size);
  for (int d = 1; d <= feature_size; ++d) {
    std::vector<double> bounds;
    for (const auto &kv : curr_family) {
      if (!kv.first.count(d))
        continue;
      for (const auto &leaf : kv.second->leaves) {
        bounds.push_back(leaf.intervals[d - 1].first);
        bounds.push_back(leaf.intervals[d - 1].second);
      }
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    lim_list[d - 1] = bounds;
  }

  acc.intercept.assign(value_size, 0.0);

  // Precompute tot(U) with half-open domain [front, back)
  std::map<std::set<int>, double, utils::setComp> tot_cache;
  auto get_tot_for_U = [&](const std::set<int> &U) -> double {
    auto it = tot_cache.find(U);
    if (it != tot_cache.end())
      return it->second;
    for (int u : U)
      if ((int)lim_list[u - 1].size() < 2) {
        tot_cache.insert({U, 0.0});
        return 0.0;
      }
    std::vector<rpf_kd::RangeConstraint> consU;
    consU.reserve(U.size());
    for (int u : U) {
      const auto &lims = lim_list[u - 1];
      consU.push_back({u - 1, lims.front(), lims.back()});
    }
    size_t cnt = consU.empty() ? (size_t)sample_size : kdt.range_count(consU);
    double tot = (double)cnt;
    tot_cache.insert({U, tot});
    return tot;
  };

  // Exact cache for KD range_count queries, keyed by the ordered (dim, lo_idx, hi_idx) triples of the constrained dims.
  // The key is the triples themselves, so two different queries never share a cache entry.
  struct KeyHash {
    size_t operator()(const std::vector<int> &v) const
    {
      unsigned long long acc = 1469598103934665603ULL;
      for (int x : v) {
        acc ^= (unsigned long long)(unsigned int)x + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
        acc = (acc ^ (acc >> 33)) * 0xff51afd7ed558ccdULL;
      }
      return (size_t)acc;
    }
  };
  std::unordered_map<std::vector<int>, size_t, KeyHash> kd_cache;
  kd_cache.reserve(1u << 15);

  // 2) Accumulate leaf contributions from ALL trees T (any order), enumerating only S up to maxp
  for (const auto &kv : curr_family) {
    const std::set<int> &T = kv.first;
    if (T == std::set<int>{0})
      continue;
    const auto &leaves = kv.second->leaves;
    if (leaves.empty())
      continue;
    std::vector<int> Tvec(T.begin(), T.end());
    const int tdim = (int)Tvec.size();

    // The terms a leaf of T writes: for every S subset T (1 <= |S| <= maxp) and j subset S, the weight of j times
    // (-1)^(|S| - |j|) on the leaf's box restricted to j, into group (S, j); j = {} goes to the constant of S.
    struct Term {
      int jmask;
      double sign;
      std::vector<double> *constant;
      rpf_rect::RectList *rects;
      std::vector<int> jpos;
    };
    std::vector<Term> terms;
    for (int smask = 1; smask < (1 << tdim); ++smask) {
      int s_size = popcount32((unsigned)smask);
      if (s_size > maxp_interaction)
        continue;
      std::vector<int> S, spos;
      for (int b = 0; b < tdim; ++b)
        if (smask & (1 << b)) {
          S.push_back(Tvec[b]);
          spos.push_back(b);
        }
      auto &entry = acc.comps[S];
      if (entry.first.empty())
        entry.first.assign(value_size, 0.0);
      for (int jS = 0; jS < (1 << s_size); ++jS) {
        Term t;
        t.jmask = 0;
        t.constant = nullptr;
        t.rects = nullptr;
        for (int b = 0; b < s_size; ++b)
          if (jS & (1 << b)) {
            t.jmask |= (1 << spos[b]);
            t.jpos.push_back(spos[b]);
          }
        t.sign = ((s_size - (int)t.jpos.size()) % 2) ? -1.0 : 1.0;
        if (jS == 0)
          t.constant = &entry.first;
        else {
          auto g = entry.second.find(jS);
          if (g == entry.second.end())
            g = entry.second.emplace(jS, rpf_rect::RectList((int)t.jpos.size(), (int)value_size)).first;
          t.rects = &g->second;
        }
        terms.push_back(std::move(t));
      }
    }

    std::vector<std::vector<double>> contrib((size_t)1 << tdim);
    std::vector<char> has_contrib((size_t)1 << tdim, 0);
    std::vector<int32_t> lo_cell(tdim), hi_cell(tdim), jlo(tdim), jhi(tdim);
    std::vector<rpf_kd::RangeConstraint> rc_by_pos(tdim);
    std::vector<char> rc_ok(tdim, 0);
    std::vector<int> rc_lo_idx(tdim, -1), rc_hi_idx(tdim, -1);
    std::vector<rpf_kd::RangeConstraint> cons;
    cons.reserve((size_t)tdim);
    std::vector<int> key;
    key.reserve(3 * (size_t)tdim);

    for (const auto &leaf : leaves) {
      for (int b = 0; b < tdim; ++b) {
        int d = Tvec[b];
        const auto &lims = lim_list[d - 1];
        const auto &cd = cuts[d - 1];
        double a = leaf.intervals[d - 1].first, e = leaf.intervals[d - 1].second;
        // the box in cells of `cuts`; the family's outermost cut points are open (as the dense grids clamp)
        lo_cell[b] = (a <= lims.front()) ? 0 : (int32_t)(std::lower_bound(cd.begin(), cd.end(), a) - cd.begin()) + 1;
        hi_cell[b] = (e >= lims.back()) ? (int32_t)cd.size() + 1
                                        : (int32_t)(std::lower_bound(cd.begin(), cd.end(), e) - cd.begin()) + 1;
        // KD constraints on [l, r) and their exact positions in lim_list for caching
        rc_ok[b] = 0;
        if ((int)lims.size() < 2)
          continue;
        double l = std::max(a, lims.front());
        double r = std::min(e, lims.back());
        if (r > l) {
          rc_by_pos[b] = {d - 1, l, r};
          rc_ok[b] = 1;
          rc_lo_idx[b] = (int)(std::lower_bound(lims.begin(), lims.end(), l) - lims.begin());
          rc_hi_idx[b] = (int)(std::lower_bound(lims.begin(), lims.end(), r) - lims.begin());
        }
      }

      // E[f_leaf | X_j] weights for j with |j| <= maxp: cnt(box on U = T \ j) / tot(U) * value
      std::fill(has_contrib.begin(), has_contrib.end(), 0);
      for (int jmask = 0; jmask < (1 << tdim); ++jmask) {
        if (popcount32((unsigned)jmask) > maxp_interaction)
          continue;
        std::set<int> U;
        cons.clear();
        key.clear();
        bool empty_range = false;
        for (int b = 0; b < tdim; ++b) {
          if (jmask & (1 << b))
            continue;
          U.insert(Tvec[b]);
          if (!rc_ok[b]) {
            empty_range = true;
            continue;
          }
          cons.push_back(rc_by_pos[b]);
          key.push_back(Tvec[b] - 1);
          key.push_back(rc_lo_idx[b]);
          key.push_back(rc_hi_idx[b]);
        }
        size_t cnt = 0;
        if (empty_range)
          cnt = 0;
        else if (cons.empty())
          cnt = (size_t)sample_size;
        else {
          auto kIt = kd_cache.find(key);
          if (kIt != kd_cache.end())
            cnt = kIt->second;
          else {
            cnt = kdt.range_count(cons);
            kd_cache.emplace(key, cnt);
          }
        }
        double totU = get_tot_for_U(U);
        if (totU <= 0.0)
          continue;
        contrib[jmask] = ((double)cnt / totU) * leaf.value;
        has_contrib[jmask] = 1;
      }

      if (has_contrib[0])
        acc.intercept += contrib[0];
      for (const auto &t : terms) {
        if (!has_contrib[t.jmask])
          continue;
        const std::vector<double> &w = contrib[t.jmask];
        if (t.constant) {
          for (size_t k = 0; k < value_size; ++k)
            (*t.constant)[k] += t.sign * w[k];
          continue;
        }
        for (size_t i = 0; i < t.jpos.size(); ++i) {
          jlo[i] = lo_cell[t.jpos[i]];
          jhi[i] = hi_cell[t.jpos[i]];
        }
        t.rects->push(jlo.data(), jhi.data(), w.data(), t.sign);
      }
    }
  }

  // 3) Merge equal boxes within each group
  acc.sort_merge();
}

// Unified purifier entry: mode 1 = grid path, mode 2 = fast exact path (weighted rectangles, see rects.hpp)
void RandomPlantedForest::purify(int maxp_interaction, int nthreads_param, int mode)
{
  if (sample_size == 0 || X.empty())
    throw std::runtime_error("Cannot purify: no training data available.");

  // Determine threads: if user provided >0, use it; otherwise default to
  // min(object-configured nthreads, hardware concurrency)
  unsigned int threads_to_use = 0;
  if (nthreads_param > 0) {
    threads_to_use = static_cast<unsigned int>(nthreads_param);
  } else {
    unsigned int avail = std::thread::hardware_concurrency();
    unsigned int obj = static_cast<unsigned int>(std::max(1, nthreads));
    unsigned int eff_avail = (avail > 0 ? avail : 1u);
    threads_to_use = std::min(obj, eff_avail);
  }

  // Below the fit's order, the higher-order components stay zero, so they no longer sum to the forest's prediction
  const bool capped = maxp_interaction > 0 && maxp_interaction < std::min(max_interaction, feature_size);

  if (threads_to_use > 1) {
    unsigned int avail = std::thread::hardware_concurrency();
    if (avail > 0 && threads_to_use > avail) {
      warn("Requested " + std::to_string(threads_to_use) + " threads but only " + std::to_string(avail) + " available");
    }
  }
  threads_to_use = std::max(1u, threads_to_use);

  // Runs fn(family index) for all families, at most threads_to_use at a time, in batches of increasing index.
  // after_batch(start, batch) runs on the calling thread after each batch.
  auto run_families = [this, threads_to_use](
                          const std::function<void(size_t)> &fn, const std::function<void(size_t, size_t)> &after_batch
                      ) {
    for (size_t start = 0; start < this->tree_families.size(); start += (size_t)threads_to_use) {
      size_t batch = std::min<size_t>((size_t)threads_to_use, this->tree_families.size() - start);
      if (batch == 1)
        fn(start);
      else {
        std::vector<std::thread> threads(batch);
        std::vector<std::exception_ptr> errors(batch);
        for (size_t i = 0; i < batch; ++i) {
          // an exception escaping a std::thread would terminate R
          threads[i] = std::thread(
              [&fn, &errors, i](size_t f) {
                try {
                  fn(f);
                } catch (...) {
                  errors[i] = std::current_exception();
                }
              },
              start + i
          );
        }
        for (auto &th : threads)
          if (th.joinable())
            th.join();
        for (auto &e : errors)
          if (e)
            std::rethrow_exception(e);
      }
      after_batch(start, batch);
    }
  };

  // Drop the store of an earlier mode-2 purification first: it is not read here, and holding both would double the
  // memory. The forest counts as unpurified until the new components are in place.
  rect_forest_.reset();
  purified = false;
  purify_capped = false;

  if (mode != 2) {
    run_families(
        [this, maxp_interaction](size_t f) { this->purify_3_family(this->tree_families[f], maxp_interaction); },
        [](size_t, size_t) {}
    );
    purified = true;
    purify_capped = capped;
    return;
  }

  // Mode 2. Drop the dense grids of an earlier mode-1 purification.
  for (auto &fam : this->tree_families)
    for (auto &tree : fam)
      tree.second->GridLeaves = LeafGrid();

  auto forest = std::make_shared<rpf_rect::PurifiedForest>();
  // Cell coordinates shared by all families: the distinct interval end points of every feature
  forest->cuts.assign(feature_size, std::vector<double>());
  for (const auto &fam : this->tree_families)
    for (const auto &tree : fam)
      for (int d : tree.first) {
        if (d < 1)
          continue;
        for (const auto &leaf : tree.second->leaves) {
          forest->cuts[d - 1].push_back(leaf.intervals[d - 1].first);
          forest->cuts[d - 1].push_back(leaf.intervals[d - 1].second);
        }
      }
  for (auto &c : forest->cuts) {
    std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());
    c.shrink_to_fit();
  }

  std::vector<int> all_idx(sample_size);
  std::iota(all_idx.begin(), all_idx.end(), 0);
  rpf_kd::KDTree kdt(&X, all_idx, feature_size);

  // Families are merged in index order through a binary counter (pairwise merges of equal size), so the result
  // does not depend on the number of threads and the working set stays near the size of the merged result.
  std::vector<std::pair<int, rpf_rect::Accumulator>> stack;
  std::vector<rpf_rect::Accumulator> batch_acc(threads_to_use);
  run_families(
      [&](size_t f) {
        this->purify_fast_exact_family(
            this->tree_families[f], maxp_interaction, kdt, forest->cuts, batch_acc[f % threads_to_use]
        );
      },
      [&](size_t, size_t batch) {
        for (size_t i = 0; i < batch; ++i) {
          stack.emplace_back(0, std::move(batch_acc[i]));
          batch_acc[i] = rpf_rect::Accumulator();
          while (stack.size() >= 2 && stack[stack.size() - 2].first == stack.back().first) {
            stack[stack.size() - 2].second.absorb(stack.back().second);
            stack[stack.size() - 2].first += 1;
            stack.pop_back();
          }
        }
      }
  );
  rpf_rect::Accumulator total;
  for (auto &s : stack)
    total.absorb(s.second);
  stack.clear();
  forest->maxp = (maxp_interaction <= 0 || maxp_interaction > feature_size) ? feature_size : maxp_interaction;
  forest->build_from(total, (int)value_size, threads_to_use);
  rect_forest_ = forest;
  purified = true;
  purify_capped = capped;
}

// The components of purify(mode = 2) of one family on the dense grid over the family's cut points, in the format of
// mode 1's grids (for get_grid_family(), e.g. rpf_boxes()). Recomputes the family's rectangles, so it needs the
// training data, and the grid is as large as mode 1's.
GridFamilyExport RandomPlantedForest::family_grid_from_rects(size_t f) const
{
  if (sample_size == 0 || X.empty())
    throw std::runtime_error("The per-family grids of purify(mode = 2) need the training data.");
  const TreeFamily &fam = tree_families[f];
  const int maxp = rect_forest_->maxp;
  std::vector<int> all_idx(sample_size);
  std::iota(all_idx.begin(), all_idx.end(), 0);
  rpf_kd::KDTree kdt(&X, all_idx, feature_size);
  rpf_rect::Accumulator acc;
  purify_fast_exact_family(fam, maxp, kdt, rect_forest_->cuts, acc);
  rpf_rect::PurifiedForest pf;
  pf.cuts = rect_forest_->cuts;
  pf.build_from(acc, (int)value_size);

  GridFamilyExport out;
  out.lim_list.assign(feature_size, std::vector<double>());
  for (int d = 1; d <= feature_size; ++d) {
    auto &b = out.lim_list[d - 1];
    for (const auto &kv : fam) {
      if (!kv.first.count(d))
        continue;
      for (const auto &leaf : kv.second->leaves) {
        b.push_back(leaf.intervals[d - 1].first);
        b.push_back(leaf.intervals[d - 1].second);
      }
    }
    std::sort(b.begin(), b.end());
    b.erase(std::unique(b.begin(), b.end()), b.end());
  }
  // the family's trees and their subsets up to maxp (as mode 1 creates them); trees above maxp stay zero
  std::set<std::set<int>, utils::setComp> keys;
  for (const auto &kv : fam) {
    keys.insert(kv.first);
    if (kv.first == std::set<int>{0})
      continue;
    std::vector<int> T(kv.first.begin(), kv.first.end());
    for (int mask = 1; mask < (1 << (int)T.size()); ++mask) {
      std::set<int> S;
      for (size_t b = 0; b < T.size(); ++b)
        if (mask & (1 << b))
          S.insert(T[b]);
      if ((int)S.size() <= maxp)
        keys.insert(S);
    }
  }
  keys.insert(std::set<int>{0});
  for (const auto &S : keys) {
    GridTreeExport t;
    t.variables = S;
    if (S == std::set<int>{0}) {
      t.dims = {1};
      t.values.assign(1, pf.intercept);
      out.trees.push_back(std::move(t));
      continue;
    }
    std::vector<int> Sv(S.begin(), S.end());
    size_t n_cells = 1;
    for (int d : Sv) {
      t.dims.push_back(std::max(1, (int)out.lim_list[d - 1].size()));
      n_cells *= (size_t)t.dims.back();
    }
    t.values.assign(n_cells, std::vector<double>(value_size, 0.0));
    auto it = pf.comps.find(Sv);
    if (it != pf.comps.end()) {
      std::vector<int32_t> cells(feature_size, 0);
      std::vector<int> idx(Sv.size(), 0);
      for (size_t e = 0; e < n_cells; ++e) {
        // column-major position e -> grid index; the last index of each dimension is not a cell and stays zero
        size_t rem = e;
        bool inside = true;
        for (size_t q = 0; q < Sv.size(); ++q) {
          idx[q] = (int)(rem % (size_t)t.dims[q]);
          rem /= (size_t)t.dims[q];
          inside = inside && idx[q] < t.dims[q] - 1;
        }
        if (!inside)
          continue;
        for (size_t q = 0; q < Sv.size(); ++q)
          cells[Sv[q] - 1] = pf.cell(Sv[q] - 1, out.lim_list[Sv[q] - 1][idx[q]]);
        pf.add_component(it->second, cells.data(), t.values[e].data(), (int)value_size);
      }
    }
    out.trees.push_back(std::move(t));
  }
  return out;
}

// Purify a single family, but only materialize outputs up to maxp_interaction.
// Higher-order trees (|dims| > maxp_interaction) are left with zero-valued grids,
// but are still used as sources during purification so that lower-order components
// are computed correctly.
void RandomPlantedForest::purify_3_family(TreeFamily &curr_family, int maxp_interaction)
{
  // Normalize cap: treat 0 (or out-of-range) as full order p = feature_size
  if (maxp_interaction <= 0 || maxp_interaction > feature_size)
    maxp_interaction = feature_size;

  // lim_list is a list giving for each variable all interval end-points
  std::vector<std::vector<double>> lim_list(feature_size);

  // go through all variables of the component
  for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim) {
    std::vector<double> bounds;

    // go through trees of family
    for (const auto &curr_tree : curr_family) {
      // consider only relevant trees that have current dimension as variable
      if (!curr_tree.first.count(curr_dim))
        continue;
      // go through leaves of tree
      for (const auto &curr_leaf : curr_tree.second->leaves) {
        // get interval ends of variable
        bounds.push_back(curr_leaf.intervals[curr_dim - 1].first);
        bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
      }
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    lim_list[curr_dim - 1] = bounds;
  }

  // Precompute per-sample bin indices for each feature based on lim_list
  // -1 means the sample falls outside the covered bounds for that feature
  std::vector<std::vector<int>> sample_bins;
  if (sample_size > 0 && feature_size > 0) {
    sample_bins.assign(sample_size, std::vector<int>(feature_size, -1));
    for (int s = 0; s < sample_size; ++s) {
      const auto &xrow = X[s];
      for (int d = 1; d <= feature_size; ++d) {
        const auto &lims = lim_list[d - 1];
        if (lims.empty())
          continue;
        const double val = xrow[d - 1];
        auto it = std::upper_bound(lims.begin(), lims.end(), val);
        int pos = static_cast<int>(it - lims.begin());
        if (pos == 0 || pos >= static_cast<int>(lims.size())) {
          sample_bins[s][d - 1] = -1; // outside
        } else {
          sample_bins[s][d - 1] = pos - 1; // interval index in [0, lims.size()-2]
        }
      }
    }
  }

  // initialize values and individuals for each tree in family
  std::vector<grid::NDGrid> grids(curr_family.size() - 1);
  std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
  std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
  std::vector<utils::Matrix<std::vector<double>>> values_old(curr_family.size() - 1);
  std::vector<std::set<int>> variables(curr_family.size() - 1);

  //  ------------- setup finer grid  -------------
  int tree_index = 0;
  for (const auto &curr_tree : curr_family) {
    if (curr_tree.first == std::set<int>{0}) {
      continue; // ignore null tree
    }

    // fill space with dimensions
    std::vector<int> dimensions;
    dimensions.reserve(curr_tree.first.size());
    for (const auto &dim : curr_tree.first) {
      dimensions.push_back(lim_list[dim - 1].size());
    }

    // setup grid for leaf indices
    auto grid = grid::NDGrid(dimensions);

    // initialize data for current tree
    grids[tree_index] = grid;
    individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
    values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));
    values_old[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));
    variables[tree_index] = curr_tree.first;

    // 1) Fill individuals using precomputed sample bins
    if (!curr_tree.first.empty()) {
      std::vector<int> point;
      point.reserve(curr_tree.first.size());
      for (int s = 0; s < sample_size; ++s) {
        point.clear();
        bool outside = false;
        for (const auto &dim : curr_tree.first) {
          int b = sample_bins.empty() ? -1 : sample_bins[s][dim - 1];
          if (b < 0) {
            outside = true;
            break;
          }
          point.push_back(b);
        }
        if (!outside) {
          individuals[tree_index][point] += 1;
        }
      }
    }

    // 2) Values accumulation: leaf-centric rectangular updates over the grid
    if (!curr_tree.first.empty()) {
      const size_t nd = curr_tree.first.size();
      // For each leaf, determine covered index ranges along each dim, then add leaf.value to all covered grid cells
      for (const auto &leaf : curr_tree.second->get_leaves()) {
        std::vector<int> start(nd, 0), stop(nd, -1);
        size_t idx_dim = 0;
        bool empty = false;
        for (const auto &dim : curr_tree.first) {
          const auto &lims = lim_list[dim - 1];
          const int dim_len = static_cast<int>(grids[tree_index].dimensions[idx_dim]);
          const int cell_max = (dim_len >= 2) ? (dim_len - 2) : -1;
          const double left = leaf.intervals[dim - 1].first;
          const double right = leaf.intervals[dim - 1].second;
          int k_low = static_cast<int>(std::lower_bound(lims.begin(), lims.end(), left) - lims.begin());
          int ub = static_cast<int>(std::upper_bound(lims.begin(), lims.end(), right) - lims.begin());
          int k_high = ub - 2; // we need lims[k+1] <= right
          if (k_low < 0)
            k_low = 0;
          if (k_high > cell_max)
            k_high = cell_max;
          if (k_low > k_high) {
            empty = true;
            break;
          }
          start[idx_dim] = k_low;
          stop[idx_dim] = k_high;
          ++idx_dim;
        }
        if (empty)
          continue;

        // Iterate over cartesian product of [start[d], stop[d]] for all dims d
        std::vector<int> gridPoint = start;
        while (true) {
          values[tree_index][gridPoint] += leaf.value;
          values_old[tree_index][gridPoint] += leaf.value;
          // increment like odometer
          if (nd == 0)
            break;
          size_t pos = nd;
          while (pos > 0) {
            --pos;
            if (gridPoint[pos] < stop[pos]) {
              ++gridPoint[pos];
              break;
            }
            gridPoint[pos] = start[pos];
          }
          if (pos == 0 && gridPoint[pos] == start[pos])
            break; // finished full cycle
        }
      }
    }

    ++tree_index;
  }

  // ------------- create new trees -------------
  grids.insert(grids.begin(), grid::NDGrid());
  values.insert(
      values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0))
  );
  values_old.insert(
      values_old.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0))
  );
  individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
  variables.insert(variables.begin(), std::set<int>{0});

  unsigned int curr_max = curr_family.rbegin()->first.size();
  while (curr_max > 1) {
    auto keys = getKeys(curr_family);
    for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key) {
      auto &curr_tree = curr_family[(*key)];
      std::set<int> curr_dims = curr_tree->split_dims;
      if (curr_dims.size() == curr_max) {
        int dim_index2 = 0;
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          if (curr_tree->split_dims.count(feature_dim) != 0) {
            std::set<int> tree_dims = curr_tree->split_dims;
            tree_dims.erase(tree_dims.find(feature_dim));
            std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
            if (!tree) {
              auto old_tree_index =
                  std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
              curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
              auto tree_index2 = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
              std::vector<int> matrix_dimensions = values[old_tree_index].dims;
              matrix_dimensions.erase(matrix_dimensions.begin() + dim_index2);
              auto grid = grid::NDGrid(matrix_dimensions);
              grids.insert(grids.begin() + tree_index2, grid);
              values.insert(
                  values.begin() + tree_index2,
                  utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0))
              );
              values_old.insert(
                  values_old.begin() + tree_index2,
                  utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0))
              );
              individuals.insert(individuals.begin() + tree_index2, utils::Matrix<int>(matrix_dimensions));
              variables.insert(variables.begin() + tree_index2, tree_dims);
              // fill individuals of new trees using precomputed sample bins
              if (!tree_dims.empty()) {
                std::vector<int> point2;
                point2.reserve(tree_dims.size());
                for (int s = 0; s < sample_size; ++s) {
                  point2.clear();
                  bool outside2 = false;
                  for (const auto &dim2 : tree_dims) {
                    int b2 = sample_bins.empty() ? -1 : sample_bins[s][dim2 - 1];
                    if (b2 < 0) {
                      outside2 = true;
                      break;
                    }
                    point2.push_back(b2);
                  }
                  if (!outside2) {
                    individuals[tree_index2][point2] += 1;
                  }
                }
              }
            }
            dim_index2++;
          }
        }
      }
    }
    --curr_max;
  }

  // ------------- purify -------------
  std::vector<std::vector<int>> dim_to_pos(variables.size(), std::vector<int>(feature_size + 1, -1));
  for (size_t idx = 0; idx < variables.size(); ++idx) {
    int pos = 0;
    for (const auto dim : variables[idx]) {
      if (dim >= 0 && dim <= feature_size)
        dim_to_pos[idx][dim] = pos++;
    }
  }

  std::vector<double> total_individuals(variables.size(), 0.0);
  for (size_t idx = 0; idx < variables.size(); ++idx) {
    double tot = 0.0;
    if (variables[idx] == std::set<int>{0}) {
      std::vector<int> only{0};
      tot += individuals[idx][only];
    } else {
      auto grid_sum = grids[idx];
      while (!grid_sum.nextPoint()) {
        auto gp = grid_sum.getPoint();
        tot += individuals[idx][gp];
      }
    }
    total_individuals[idx] = tot;
  }

  int tree_index_t = curr_family.size() - 1;
  for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t) {
    std::set<int> curr_dims = *tree_t;
    if (curr_dims == std::set<int>{0})
      continue;

    auto grid = grids[tree_index_t];
    int tree_index_u = variables.size();
    for (auto tree_u = variables.rbegin(); tree_u != variables.rend(); ++tree_u) {
      --tree_index_u;
      std::set<int> j_dims = curr_dims;
      if (tree_u->size() > curr_dims.size())
        continue;
      bool subset = true;
      for (const auto dim : *tree_u) {
        if (tree_t->count(dim) == 0) {
          subset = false;
          break;
        }
        j_dims.erase(dim);
      }
      if (!subset)
        continue;

      double tot_sum = total_individuals[tree_index_u];
      if (tot_sum == 0.0)
        continue;
      const double inv_tot_sum = 1.0 / tot_sum;

      grid = grids[tree_index_u];
      std::vector<double> update(value_size, 0);

      if (j_dims.size() == 0) {
        while (!grid.nextPoint()) {
          auto gridPoint_i = grid.getPoint();
          double curr_sum = individuals[tree_index_u][gridPoint_i];
          update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_i];
        }

        int tree_index_s = variables.size();
        for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s) {
          --tree_index_s;
          if (*tree_s == std::set<int>{0}) {
            auto gridPoint_0 = std::vector<int>{0};
            values[tree_index_s][gridPoint_0] += update;
          } else {
            bool subset2 = true;
            for (const auto dim : *tree_s) {
              if (tree_t->count(dim) == 0) {
                subset2 = false;
                break;
              }
            }
            if (!subset2)
              continue;
            if (maxp_interaction > 0 && tree_s->size() > (size_t)maxp_interaction)
              continue; // skip materializing > cap
            auto grid_k = grids[tree_index_s];
            while (!grid_k.nextPoint()) {
              auto gridPoint_k = grid_k.getPoint();
              int sign0 = ((*tree_s).size() % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_k] += sign0 * update;
            }
          }
        }
      } else {
        std::vector<int> j_sizes(j_dims.size(), 0);
        for (size_t j = 0; j < j_dims.size(); ++j) {
          auto tmp = j_dims.begin();
          std::advance(tmp, j);
          int j_index = dim_to_pos[tree_index_t][*tmp];
          j_sizes[j] = grids[tree_index_t].dimensions[j_index];
        }
        grid::NDGrid grid_j = grid::NDGrid(j_sizes);
        while (!grid_j.nextPoint()) {
          std::vector<double> update(value_size, 0);
          auto gridPoint_j = grid_j.getPoint();
          grid = grids[tree_index_u];
          while (!grid.nextPoint()) {
            auto gridPoint_i = grid.getPoint();
            double curr_sum = individuals[tree_index_u][gridPoint_i];
            std::vector<int> gridPoint_ij(tree_t->size(), 0);
            for (size_t j = 0; j < gridPoint_j.size(); ++j) {
              auto j_dim = j_dims.begin();
              std::advance(j_dim, j);
              int j_index = dim_to_pos[tree_index_t][*j_dim];
              gridPoint_ij[j_index] = gridPoint_j[j];
            }
            for (size_t i = 0; i < gridPoint_i.size(); ++i) {
              auto i_dim = tree_u->begin();
              std::advance(i_dim, i);
              int i_index = dim_to_pos[tree_index_t][*i_dim];
              gridPoint_ij[i_index] = gridPoint_i[i];
            }
            update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_ij];
          }

          int tree_index_s = variables.size();
          for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s) {
            --tree_index_s;
            bool subset2 = true;
            for (const auto dim : j_dims) {
              if (tree_s->count(dim) == 0) {
                subset2 = false;
                break;
              }
            }
            for (const auto dim : *tree_s) {
              if (tree_t->count(dim) == 0) {
                subset2 = false;
                break;
              }
            }
            if (!subset2)
              continue;
            // Skip writing for components above the cap
            if (maxp_interaction > 0 && tree_s->size() > (size_t)maxp_interaction)
              continue;

            std::set<int> k_dims = *tree_s;
            std::set<int> k_dims_h1 = *tree_s;
            std::set<int> k_dims_h2 = *tree_u;
            for (const auto dim : *tree_u)
              k_dims.insert(dim);
            for (const auto dim : *tree_s)
              k_dims_h2.erase(dim);
            for (const auto dim : *tree_u)
              k_dims_h1.erase(dim);
            for (const auto dim : k_dims_h1)
              k_dims.erase(dim);
            for (const auto dim : k_dims_h2)
              k_dims.erase(dim);

            if (k_dims.size() == 0) {
              size_t diff = (*tree_s).size() - j_dims.size();
              int sign = (diff % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_j] += sign * update;
            } else {
              std::vector<int> k_sizes(k_dims.size(), 0);
              for (size_t k = 0; k < k_dims.size(); ++k) {
                auto tmp = k_dims.begin();
                std::advance(tmp, k);
                int k_index = dim_to_pos[tree_index_t][*tmp];
                k_sizes[k] = grids[tree_index_t].dimensions[k_index];
              }
              grid::NDGrid grid_k = grid::NDGrid(k_sizes);
              while (!grid_k.nextPoint()) {
                auto gridPoint_k = grid_k.getPoint();
                std::vector<int> gridPoint_jk(tree_s->size(), 0);
                for (size_t j = 0; j < gridPoint_j.size(); ++j) {
                  auto j_dim = j_dims.begin();
                  std::advance(j_dim, j);
                  int j_index = dim_to_pos[tree_index_s][*j_dim];
                  gridPoint_jk[j_index] = gridPoint_j[j];
                }
                for (size_t k = 0; k < gridPoint_k.size(); ++k) {
                  auto k_dim = k_dims.begin();
                  std::advance(k_dim, k);
                  int k_index = dim_to_pos[tree_index_s][*k_dim];
                  gridPoint_jk[k_index] = gridPoint_k[k];
                }
                size_t diff = (*tree_s).size() - j_dims.size();
                int sign2 = (diff % 2 == 0) ? 1 : -1;
                values[tree_index_s][gridPoint_jk] += sign2 * update;
              }
            }
          }
        }
      }
    }
    --tree_index_t;
  }

  // ------------- attach to rpf class -------------
  for (size_t tree_index3 = 0; tree_index3 < variables.size(); ++tree_index3) {
    LeafGrid curr_gridLeaf;
    curr_gridLeaf.grid = grids[tree_index3];
    curr_gridLeaf.individuals = individuals[tree_index3];
    curr_gridLeaf.lim_list = lim_list;
    // If this tree exceeds the cap, attach a zero-valued matrix of the correct shape
    if (maxp_interaction > 0 && variables[tree_index3] != std::set<int>{0} &&
        variables[tree_index3].size() > (size_t)maxp_interaction) {
      curr_gridLeaf.values =
          utils::Matrix<std::vector<double>>(grids[tree_index3].dimensions, std::vector<double>(value_size, 0));
    } else {
      curr_gridLeaf.values = values[tree_index3];
    }
    curr_family[variables[tree_index3]]->GridLeaves = curr_gridLeaf;
  }
}
