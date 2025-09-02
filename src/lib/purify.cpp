#include "rpf.hpp"
#include "kdtree.hpp"
#include "diffbuf.hpp"
#include <limits>
#include <cmath>
#include <unordered_map>

// Generates the next combination of k indices from a set of n elements.
static inline bool next_combination(std::vector<int> &p, int n)
{
  int k = (int)p.size();
  for (int i = k - 1; i >= 0; --i)
  {
    if (p[i] < n - k + i)
    {
      p[i]++;
      for (int j = i + 1; j < k; ++j)
      {
        p[j] = p[j - 1] + 1;
      }
      return true;
    }
  }
  return false;
}

void RandomPlantedForest::purify_1()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (auto tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    while (curr_max >= 1)
    {

      // go through split dimensions of all trees
      auto keys = getKeys(curr_family);
      std::vector<std::set<int>>::reverse_iterator key = keys.rbegin();
      while (key != keys.rend())
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (curr_max == 1)
              {
                tree = curr_family[std::set<int>{0}];
              }
              else
              {
                if (!tree)
                {
                  curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                  tree = curr_family[tree_dims];
                }
              }

              // go through leaves of current tree
              int n_leaves = curr_tree->leaves.size();
              for (int l = 0; l < n_leaves; ++l)
              {
                auto &curr_leaf = curr_tree->leaves[l];

                double multiplier = (curr_leaf.intervals[feature_dim - 1].second - curr_leaf.intervals[feature_dim - 1].first) / (upper_bounds[feature_dim - 1] - lower_bounds[feature_dim - 1]);

                // new leaf including intervals and value
                Leaf new_leaf = curr_leaf; // initialize intervals with first leaf
                new_leaf.intervals[feature_dim - 1].first = lower_bounds[feature_dim - 1];
                new_leaf.intervals[feature_dim - 1].second = upper_bounds[feature_dim - 1];
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = -curr_leaf.value[i] * multiplier; // update value of new leaf

                // append new leaf
                if (!leafExists(new_leaf.intervals, curr_tree))
                  curr_tree->leaves.push_back(new_leaf);
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = curr_leaf.value[i] * multiplier; // update value of new leaf
                if (!leafExists(new_leaf.intervals, tree))
                  tree->leaves.push_back(new_leaf);
              }
            }
          }
        }
        key++;
      }

      // update currently considered dimension size
      --curr_max;
    }
  }

  purified = true;
}


void RandomPlantedForest::purify_fast_exact_family(TreeFamily &curr_family, int maxp_interaction)
{
  // Normalize cap: treat 0 (or out-of-range) as full order p = feature_size
  if (maxp_interaction <= 0 || maxp_interaction > feature_size) maxp_interaction = feature_size;
  
  // Portable 32-bit popcount to avoid compiler-specific builtins
  auto popcount32 = [](unsigned int x) -> int {
    x = x - ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    return (int)((((x + (x >> 4)) & 0x0F0F0F0Fu) * 0x01010101u) >> 24);
  };
  auto nextDown = [](double x) { return std::nextafter(x, -std::numeric_limits<double>::infinity()); };

  // 0) Ensure all subset components exist in the family (sources and targets)
  {
    auto base_keys = getKeys(curr_family);
    for (const auto &T : base_keys) {
      if (T == std::set<int>{0}) continue;
      std::vector<int> dims; dims.reserve(T.size());
      for (int d : T) dims.push_back(d);
      int k = (int)dims.size();
      for (int mask = 1; mask < (1 << k); ++mask) {
        if (maxp_interaction > 0) {
          int bits = popcount32((unsigned)mask);
          if (bits > maxp_interaction) continue;
        }
        std::set<int> S;
        for (int b = 0; b < k; ++b) if (mask & (1 << b)) S.insert(dims[b]);
        if (curr_family.find(S) == curr_family.end()) {
          curr_family.insert({S, std::make_shared<DecisionTree>(DecisionTree(S))});
        }
      }
    }
    if (curr_family.find(std::set<int>{0}) == curr_family.end()) {
      curr_family.insert({std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}))});
    }
  }

  // 1) Build lim_list (unique cut endpoints per feature)
  std::vector<std::vector<double>> lim_list(feature_size);
  for (int d = 1; d <= feature_size; ++d) {
    std::vector<double> bounds;
    for (const auto &kv : curr_family) {
      if (!kv.first.count(d)) continue;
      for (const auto &leaf : kv.second->leaves) {
        bounds.push_back(leaf.intervals[d - 1].first);
        bounds.push_back(leaf.intervals[d - 1].second);
      }
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    lim_list[d - 1] = bounds;
  }

  // Precompute number of cells per feature (endpoints - 1), clamped at 0
  std::vector<int> cells_by_dim(feature_size + 1, 0);
  for (int d = 1; d <= feature_size; ++d) cells_by_dim[d] = std::max(0, (int)lim_list[d - 1].size() - 1);

  // 2) Prepare per-S diff buffers (emit only S with |S|<=maxp; keep intercept)
  auto keys = getKeys(curr_family);
  std::vector<std::set<int>> S_vars; S_vars.reserve(keys.size());
  std::vector<rpf_diff::NDArray<std::vector<double>>> diff_S; diff_S.reserve(keys.size());
  std::vector<double> intercept(value_size, 0.0);
  for (const auto &S : keys) {
    if (S != std::set<int>{0} && maxp_interaction > 0 && (int)S.size() > maxp_interaction) continue;
    S_vars.push_back(S);
    if (S == std::set<int>{0}) diff_S.emplace_back(rpf_diff::NDArray<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    else {
      std::vector<int> diff_dims; diff_dims.reserve(S.size());
      for (int d : S) { int K = (int)lim_list[d - 1].size(); int cells = std::max(0, K - 1); diff_dims.push_back(cells + 1); }
      diff_S.emplace_back(rpf_diff::NDArray<std::vector<double>>(diff_dims, std::vector<double>(value_size, 0)));
    }
  }

  std::map<std::set<int>, int, utils::setComp> s_index_map; for (size_t i = 0; i < S_vars.size(); ++i) s_index_map[S_vars[i]] = (int)i;
  auto set_to_vec = [](const std::set<int> &S){ std::vector<int> v; v.reserve(S.size()); for (int x : S) v.push_back(x); return v; };

  // 3) KD-tree over all samples
  std::vector<int> all_idx(sample_size); for (int i = 0; i < sample_size; ++i) all_idx[i] = i;
  rpf_kd::KDTree kdt(&X, all_idx, feature_size);

  // Precompute tot(U) with half-open domain [front, back)
  std::map<std::set<int>, double, utils::setComp> tot_cache;
  auto get_tot_for_U = [&](const std::set<int>& U)->double {
    auto it = tot_cache.find(U); if (it != tot_cache.end()) return it->second;
    for (int u : U) if ((int)lim_list[u - 1].size() < 2) { tot_cache.insert({U, 0.0}); return 0.0; }
    std::vector<rpf_kd::RangeConstraint> consU; consU.reserve(U.size());
    for (int u : U) { const auto &lims = lim_list[u - 1]; double lo = lims.front(); double hi = nextDown(lims.back()); consU.push_back({u - 1, lo, hi}); }
    size_t cnt = consU.empty() ? (size_t)sample_size : kdt.range_count(consU);
    double tot = (double)cnt; tot_cache.insert({U, tot}); return tot;
  };

  // Exact cache for KD range_count queries keyed by (dim, lo_idx, hi_idx) triples per constrained dim
  // Key construction: 64-bit hash mixed from ordered triples to avoid building strings/sets
  auto mix64 = [](unsigned long long x){
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x;
  };
  auto pack3 = [&](unsigned long long acc, int d, int lo, int hi){
    unsigned long long k = ((unsigned long long)(unsigned int)d << 32) ^ ((unsigned long long)(unsigned int)lo << 16) ^ (unsigned long long)(unsigned int)hi;
    acc ^= mix64(k + 0x9e3779b97f4a7c15ULL + (acc<<6) + (acc>>2));
    return acc;
  };
  std::unordered_map<unsigned long long, size_t> kd_cache; kd_cache.reserve(1u << 15);

  // 4) Accumulate leaf contributions from ALL trees T (any order), enumerating only S up to maxp
  for (const auto &kv : curr_family) {
    const std::set<int> &T = kv.first; if (T == std::set<int>{0}) continue;
    const auto &leaves = kv.second->get_leaves();
    std::vector<int> Tvec = set_to_vec(T); const int tdim = (int)Tvec.size();
    // map from dimension id -> position index in Tvec
    std::vector<int> pos_in_T(feature_size + 1, -1);
    for (int i = 0; i < tdim; ++i) pos_in_T[Tvec[i]] = i;

    for (const auto &leaf : leaves) {
      // Pre-cache per-dim grid cell ranges for this leaf
      std::vector<int> lo_cached(feature_size + 1, 0), hi_cached(feature_size + 1, 0);
      for (int d : T) {
        const auto &lims = lim_list[d - 1]; int cells = std::max(0, (int)lims.size() - 1);
        int k_low = (int)(std::lower_bound(lims.begin(), lims.end(), leaf.intervals[d - 1].first) - lims.begin());
        int ub = (int)(std::upper_bound(lims.begin(), lims.end(), leaf.intervals[d - 1].second) - lims.begin());
        int k_high_cell = std::max(0, ub - 2);
        lo_cached[d] = std::max(0, k_low);
        hi_cached[d] = std::min(cells, k_high_cell + 1);
      }

      // Precompute per-dimension KD constraints for this leaf
      std::vector<rpf_kd::RangeConstraint> rc_by_dim(feature_size + 1);
      std::vector<char> rc_ok(feature_size + 1, 0);
      // Also store lim_list boundary indices for exact caching
      std::vector<int> rc_lo_idx(feature_size + 1, -1);
      std::vector<int> rc_hi_idx(feature_size + 1, -1);
      for (int d : T) {
        const auto &lims = lim_list[d - 1];
        if ((int)lims.size() < 2) { rc_ok[d] = 0; continue; }
        double l = std::max(leaf.intervals[d - 1].first, lims.front());
        double r = std::min(leaf.intervals[d - 1].second, lims.back());
        double hi = nextDown(r);
        if (!(hi >= l)) { rc_ok[d] = 0; }
        else {
          rc_by_dim[d] = {d - 1, l, hi}; rc_ok[d] = 1;
          int lidx = (int)(std::lower_bound(lims.begin(), lims.end(), l) - lims.begin());
          int ridx = (int)(std::lower_bound(lims.begin(), lims.end(), r) - lims.begin());
          rc_lo_idx[d] = lidx; rc_hi_idx[d] = ridx;
        }
      }

      // Precompute E[f_T | X_j] only for j with |j| <= maxp_interaction by enumerating combinations
      // Store by mask over positions in T (0..tdim-1) to avoid building a full 2^tdim table
      std::unordered_map<int, std::vector<double>> contrib_by_mask;
      contrib_by_mask.reserve(32u);
      std::vector<rpf_kd::RangeConstraint> cons; cons.reserve((size_t)tdim);

      auto compute_for_j = [&](const std::vector<int> &j_pos){
        // Build complement U = T \ j and corresponding KD constraints
        std::set<int> U;
        U.clear();
        int jmask = 0;
        std::vector<char> is_in_j((size_t)tdim, 0);
        for (int pos : j_pos) { if (pos >= 0 && pos < tdim) { is_in_j[(size_t)pos] = 1; jmask |= (1 << pos); } }
        for (int b = 0; b < tdim; ++b) if (!is_in_j[(size_t)b]) U.insert(Tvec[b]);
        cons.clear(); cons.reserve(U.size());
        bool empty_range = false;
        // Build exact cache key from ordered (dim, lo_idx, hi_idx)
        unsigned long long key = 1469598103934665603ULL; // FNV offset basis-ish seed
        for (int u : U) {
          if (!rc_ok[u]) { empty_range = true; break; }
          cons.push_back(rc_by_dim[u]);
          key = pack3(key, u - 1, rc_lo_idx[u], rc_hi_idx[u]);
        }
        size_t cnt = 0;
        if (empty_range) cnt = 0;
        else if (cons.empty()) cnt = (size_t)sample_size;
        else {
          auto kIt = kd_cache.find(key);
          if (kIt != kd_cache.end()) cnt = kIt->second;
          else { cnt = kdt.range_count(cons); kd_cache.emplace(key, cnt); }
        }
        double totU = get_tot_for_U(U); if (totU <= 0.0) return;
        contrib_by_mask[jmask] = ((double)cnt / totU) * leaf.value;
      };

      // j size = 0
      compute_for_j(std::vector<int>{});
      // j sizes 1..min(tdim, maxp_interaction)
      int maxk = std::min(tdim, maxp_interaction);
      for (int k = 1; k <= maxk; ++k) {
        std::vector<int> p(k); for (int i = 0; i < k; ++i) p[i] = i;
        do { compute_for_j(p); } while (next_combination(p, tdim));
      }

      // Efficiently iterate directly over target subsets S up to size maxp_interaction
      for (int k = 0; k <= std::min(tdim, maxp_interaction); ++k) {
        std::vector<int> p(k); for (int i = 0; i < k; ++i) p[i] = i;
        if (k == 0) {
          // j = {} corresponds to mask 0
          auto it0 = contrib_by_mask.find(0);
          if (it0 != contrib_by_mask.end()) intercept += it0->second;
        } else {
          do {
            std::set<int> S; for (int idx : p) S.insert(Tvec[idx]);
            // Inclusion-exclusion over all j subset S, writing per-term rectangles
            auto itS = s_index_map.find(S); if (itS == s_index_map.end()) { /* nothing to write */ }
            else {
              int s_idx = itS->second;
              std::vector<int> Svec = set_to_vec(S); const int s_dim = (int)Svec.size();
              for (int sm = 0; sm < (1 << s_dim); ++sm) {
                int jmask_on_T = 0; int jcount = 0;
                for (int b = 0; b < s_dim; ++b) {
                  if (sm & (1 << b)) { ++jcount; int d = Svec[b]; int pos = pos_in_T[d]; if (pos >= 0) jmask_on_T |= (1 << pos); }
                }
                auto jit = contrib_by_mask.find(jmask_on_T);
                if (jit == contrib_by_mask.end()) continue;
                const std::vector<double> &contrib_j = jit->second;
                int sign_flip = ((int)S.size() - jcount) % 2;
                std::vector<double> signed_contrib = sign_flip ? (contrib_j * (-1)) : contrib_j;
                // Build rectangle: restrict dims in j to leaf's range; others span entire domain
                std::vector<int> lo; lo.reserve(S.size()); std::vector<int> hi; hi.reserve(S.size());
                for (int di = 0; di < s_dim; ++di) {
                  int d = Svec[di];
                  if (sm & (1 << di)) { lo.push_back(lo_cached[d]); hi.push_back(hi_cached[d]); }
                  else { lo.push_back(0); hi.push_back(cells_by_dim[d]); }
                }
                rpf_diff::add_rect(diff_S[s_idx], lo, hi, signed_contrib);
              }
            }
          } while (next_combination(p, tdim));
        }
      }
    }
  }

  // 5) Finalize per S
  for (size_t i = 0; i < S_vars.size(); ++i) {
    const auto &S = S_vars[i]; LeafGrid gl; gl.lim_list = lim_list;
    if (S == std::set<int>{0}) {
      gl.grid = grid::NDGrid(); gl.values = utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)); gl.individuals = utils::Matrix<int>(std::vector<int>{1}, 0);
      std::vector<int> idx0{0}; gl.values[idx0] = intercept;
    } else {
      std::vector<int> dims_end; std::vector<int> cells_dims; for (int d : S) { int K = (int)lim_list[d - 1].size(); dims_end.push_back(std::max(1, K)); cells_dims.push_back(std::max(0, K - 1)); }
      rpf_diff::inclusive_scan_inplace(diff_S[i]); gl.grid = grid::NDGrid(dims_end);
      gl.values = utils::Matrix<std::vector<double>>(dims_end, std::vector<double>(value_size, 0)); gl.individuals = utils::Matrix<int>(dims_end, 0);
      auto g = grid::NDGrid(cells_dims); while (!g.nextPoint()) { auto point = g.getPoint(); gl.values[point] = diff_S[i].at(point); }
    }
    curr_family[S]->GridLeaves = gl;
  }

  // 6) Overwrite high orders with zeros if capped
  if (maxp_interaction > 0) {
    for (const auto &S : keys) {
      if (S == std::set<int>{0} || (int)S.size() <= maxp_interaction) continue;
      LeafGrid gl; gl.lim_list = lim_list; std::vector<int> dims_end; for (int d : S) { int K = (int)lim_list[d - 1].size(); dims_end.push_back(std::max(1, K)); }
      gl.grid = grid::NDGrid(dims_end); gl.values = utils::Matrix<std::vector<double>>(dims_end, std::vector<double>(value_size, 0)); gl.individuals = utils::Matrix<int>(dims_end, 0);
      curr_family[S]->GridLeaves = gl;
    }
  }
}





// Unified purifier entry: mode 1 = grid path, mode 2 = fast exact path
void RandomPlantedForest::purify(int maxp_interaction, int nthreads_param, int mode)
{
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

  auto worker = [this, maxp_interaction, mode](TreeFamily &fam){
    if (mode == 2) this->purify_fast_exact_family(fam, maxp_interaction);
    else this->purify_3_family(fam, maxp_interaction);
  };

  if (threads_to_use > 1)
  {
    unsigned int avail = std::thread::hardware_concurrency();
    if (avail > 0 && threads_to_use > avail)
    {
      Rcout << "Requested " << threads_to_use << " threads but only " << avail << " available" << std::endl;
    }
    for (size_t start = 0; start < this->tree_families.size(); start += (size_t)threads_to_use)
    {
      size_t batch = std::min<size_t>((size_t)threads_to_use, this->tree_families.size() - start);
      if (batch == 0) break;
      std::vector<std::thread> threads(batch);
      for (size_t i = 0; i < batch; ++i)
      {
        size_t fam_index = start + i;
        threads[i] = std::thread([&worker](TreeFamily *fam_ptr){ worker(*fam_ptr); }, &this->tree_families[fam_index]);
      }
      for (auto &th : threads)
      {
        if (th.joinable()) th.join();
      }
    }
    purified = true;
    return;
  }

  for (auto &fam : this->tree_families) worker(fam);
  purified = true;
}



// Purify a single family, but only materialize outputs up to maxp_interaction.
// Higher-order trees (|dims| > maxp_interaction) are left with zero-valued grids,
// but are still used as sources during purification so that lower-order components
// are computed correctly.
void RandomPlantedForest::purify_3_family(TreeFamily &curr_family, int maxp_interaction)
{
  // Normalize cap: treat 0 (or out-of-range) as full order p = feature_size
  if (maxp_interaction <= 0 || maxp_interaction > feature_size) maxp_interaction = feature_size;
  
  // lim_list is a list giving for each variable all interval end-points
  std::vector<std::vector<double>> lim_list(feature_size);

  // go through all variables of the component
  for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
  {
    std::vector<double> bounds;

    // go through trees of family
    for (const auto &curr_tree : curr_family)
    {
      // consider only relevant trees that have current dimension as variable
      if (!curr_tree.first.count(curr_dim))
        continue;
      // go through leaves of tree
      for (const auto &curr_leaf : curr_tree.second->leaves)
      {
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
  if (sample_size > 0 && feature_size > 0)
  {
    sample_bins.assign(sample_size, std::vector<int>(feature_size, -1));
    for (int s = 0; s < sample_size; ++s)
    {
      const auto &xrow = X[s];
      for (int d = 1; d <= feature_size; ++d)
      {
        const auto &lims = lim_list[d - 1];
        if (lims.empty()) continue;
        const double val = xrow[d - 1];
        auto it = std::upper_bound(lims.begin(), lims.end(), val);
        int pos = static_cast<int>(it - lims.begin());
        if (pos == 0 || pos >= static_cast<int>(lims.size()))
        {
          sample_bins[s][d - 1] = -1; // outside
        }
        else
        {
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
  for (const auto &curr_tree : curr_family)
  {
    if (curr_tree.first == std::set<int>{0})
    {
      continue; // ignore null tree
    }

    // fill space with dimensions
    std::vector<int> dimensions;
    dimensions.reserve(curr_tree.first.size());
    for (const auto &dim : curr_tree.first)
    {
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
    if (!curr_tree.first.empty())
    {
      std::vector<int> point; point.reserve(curr_tree.first.size());
      for (int s = 0; s < sample_size; ++s)
      {
        point.clear(); bool outside = false;
        for (const auto &dim : curr_tree.first)
        {
          int b = sample_bins.empty() ? -1 : sample_bins[s][dim - 1];
          if (b < 0) { outside = true; break; }
          point.push_back(b);
        }
        if (!outside) { individuals[tree_index][point] += 1; }
      }
    }

    // 2) Values accumulation: leaf-centric rectangular updates over the grid
    if (!curr_tree.first.empty())
    {
      const size_t nd = curr_tree.first.size();
      // For each leaf, determine covered index ranges along each dim, then add leaf.value to all covered grid cells
      for (const auto &leaf : curr_tree.second->get_leaves())
      {
        std::vector<int> start(nd, 0), stop(nd, -1);
        size_t idx_dim = 0;
        bool empty = false;
        for (const auto &dim : curr_tree.first)
        {
          const auto &lims = lim_list[dim - 1];
          const int dim_len = static_cast<int>(grids[tree_index].dimensions[idx_dim]);
          const int cell_max = (dim_len >= 2) ? (dim_len - 2) : -1;
          const double left = leaf.intervals[dim - 1].first;
          const double right = leaf.intervals[dim - 1].second;
          int k_low = static_cast<int>(std::lower_bound(lims.begin(), lims.end(), left) - lims.begin());
          int ub = static_cast<int>(std::upper_bound(lims.begin(), lims.end(), right) - lims.begin());
          int k_high = ub - 2; // we need lims[k+1] <= right
          if (k_low < 0) k_low = 0;
          if (k_high > cell_max) k_high = cell_max;
          if (k_low > k_high) { empty = true; break; }
          start[idx_dim] = k_low;
          stop[idx_dim] = k_high;
          ++idx_dim;
        }
        if (empty) continue;

        // Iterate over cartesian product of [start[d], stop[d]] for all dims d
        std::vector<int> gridPoint = start;
        while (true)
        {
          values[tree_index][gridPoint] += leaf.value;
          values_old[tree_index][gridPoint] += leaf.value;
          // increment like odometer
          if (nd == 0) break;
          size_t pos = nd;
          while (pos > 0)
          {
            --pos;
            if (gridPoint[pos] < stop[pos]) { ++gridPoint[pos]; break; }
            gridPoint[pos] = start[pos];
          }
          if (pos == 0 && gridPoint[pos] == start[pos]) break; // finished full cycle
        }
      }
    }

    ++tree_index;
  }

  // ------------- create new trees -------------
  grids.insert(grids.begin(), grid::NDGrid());
  values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
  values_old.insert(values_old.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
  individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
  variables.insert(variables.begin(), std::set<int>{0});

  unsigned int curr_max = curr_family.rbegin()->first.size();
  while (curr_max > 1)
  {
    auto keys = getKeys(curr_family);
    for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
    {
      auto &curr_tree = curr_family[(*key)];
      std::set<int> curr_dims = curr_tree->split_dims;
      if (curr_dims.size() == curr_max)
      {
        int dim_index2 = 0;
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
        {
          if (curr_tree->split_dims.count(feature_dim) != 0)
          {
            std::set<int> tree_dims = curr_tree->split_dims;
            tree_dims.erase(tree_dims.find(feature_dim));
            std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
            if (!tree)
            {
              auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
              curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
              auto tree_index2 = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
              std::vector<int> matrix_dimensions = values[old_tree_index].dims;
              matrix_dimensions.erase(matrix_dimensions.begin() + dim_index2);
              auto grid = grid::NDGrid(matrix_dimensions);
              grids.insert(grids.begin() + tree_index2, grid);
              values.insert(values.begin() + tree_index2, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
              values_old.insert(values_old.begin() + tree_index2, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
              individuals.insert(individuals.begin() + tree_index2, utils::Matrix<int>(matrix_dimensions));
              variables.insert(variables.begin() + tree_index2, tree_dims);
              // fill individuals of new trees using precomputed sample bins
              if (!tree_dims.empty())
              {
                std::vector<int> point2; point2.reserve(tree_dims.size());
                for (int s = 0; s < sample_size; ++s)
                {
                  point2.clear(); bool outside2 = false;
                  for (const auto &dim2 : tree_dims)
                  {
                    int b2 = sample_bins.empty() ? -1 : sample_bins[s][dim2 - 1];
                    if (b2 < 0) { outside2 = true; break; }
                    point2.push_back(b2);
                  }
                  if (!outside2) { individuals[tree_index2][point2] += 1; }
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
  for (size_t idx = 0; idx < variables.size(); ++idx)
  {
    int pos = 0;
    for (const auto dim : variables[idx])
    {
      if (dim >= 0 && dim <= feature_size) dim_to_pos[idx][dim] = pos++;
    }
  }

  std::vector<double> total_individuals(variables.size(), 0.0);
  for (size_t idx = 0; idx < variables.size(); ++idx)
  {
    double tot = 0.0;
    if (variables[idx] == std::set<int>{0})
    {
      std::vector<int> only{0};
      tot += individuals[idx][only];
    }
    else
    {
      auto grid_sum = grids[idx];
      while (!grid_sum.nextPoint())
      {
        auto gp = grid_sum.getPoint();
        tot += individuals[idx][gp];
      }
    }
    total_individuals[idx] = tot;
  }

  int tree_index_t = curr_family.size() - 1;
  for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t)
  {
    std::set<int> curr_dims = *tree_t;
    if (curr_dims == std::set<int>{0})
      continue;

    auto grid = grids[tree_index_t];
    int tree_index_u = variables.size();
    for (auto tree_u = variables.rbegin(); tree_u != variables.rend(); ++tree_u)
    {
      --tree_index_u;
      std::set<int> j_dims = curr_dims;
      if (tree_u->size() > curr_dims.size())
        continue;
      bool subset = true;
      for (const auto dim : *tree_u)
      {
        if (tree_t->count(dim) == 0)
        {
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

      if (j_dims.size() == 0)
      {
        while (!grid.nextPoint())
        {
          auto gridPoint_i = grid.getPoint();
          double curr_sum = individuals[tree_index_u][gridPoint_i];
          update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_i];
        }

        int tree_index_s = variables.size();
        for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
        {
          --tree_index_s;
          if (*tree_s == std::set<int>{0})
          {
            auto gridPoint_0 = std::vector<int>{0};
            values[tree_index_s][gridPoint_0] += update;
          }
          else
          {
            bool subset2 = true;
            for (const auto dim : *tree_s)
            {
              if (tree_t->count(dim) == 0)
              {
                subset2 = false;
                break;
              }
            }
            if (!subset2)
              continue;
            if (maxp_interaction > 0 && tree_s->size() > (size_t)maxp_interaction) continue; // skip materializing > cap
            auto grid_k = grids[tree_index_s];
            while (!grid_k.nextPoint())
            {
              auto gridPoint_k = grid_k.getPoint();
              int sign0 = ((*tree_s).size() % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_k] += sign0 * update;
            }
          }
        }
      }
      else
      {
        std::vector<int> j_sizes(j_dims.size(), 0);
        for (size_t j = 0; j < j_dims.size(); ++j)
        {
          auto tmp = j_dims.begin();
          std::advance(tmp, j);
          int j_index = dim_to_pos[tree_index_t][*tmp];
          j_sizes[j] = grids[tree_index_t].dimensions[j_index];
        }
        grid::NDGrid grid_j = grid::NDGrid(j_sizes);
        while (!grid_j.nextPoint())
        {
          std::vector<double> update(value_size, 0);
          auto gridPoint_j = grid_j.getPoint();
          grid = grids[tree_index_u];
          while (!grid.nextPoint())
          {
            auto gridPoint_i = grid.getPoint();
            double curr_sum = individuals[tree_index_u][gridPoint_i];
            std::vector<int> gridPoint_ij(tree_t->size(), 0);
            for (size_t j = 0; j < gridPoint_j.size(); ++j)
            {
              auto j_dim = j_dims.begin();
              std::advance(j_dim, j);
              int j_index = dim_to_pos[tree_index_t][*j_dim];
              gridPoint_ij[j_index] = gridPoint_j[j];
            }
            for (size_t i = 0; i < gridPoint_i.size(); ++i)
            {
              auto i_dim = tree_u->begin();
              std::advance(i_dim, i);
              int i_index = dim_to_pos[tree_index_t][*i_dim];
              gridPoint_ij[i_index] = gridPoint_i[i];
            }
            update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_ij];
          }

          int tree_index_s = variables.size();
          for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
          {
            --tree_index_s;
            bool subset2 = true;
            for (const auto dim : j_dims)
            {
              if (tree_s->count(dim) == 0)
              {
                subset2 = false;
                break;
              }
            }
            for (const auto dim : *tree_s)
            {
              if (tree_t->count(dim) == 0)
              {
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

            if (k_dims.size() == 0)
            {
              size_t diff = (*tree_s).size() - j_dims.size();
              int sign = (diff % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_j] += sign * update;
            }
            else
            {
              std::vector<int> k_sizes(k_dims.size(), 0);
              for (size_t k = 0; k < k_dims.size(); ++k)
              {
                auto tmp = k_dims.begin();
                std::advance(tmp, k);
                int k_index = dim_to_pos[tree_index_t][*tmp];
                k_sizes[k] = grids[tree_index_t].dimensions[k_index];
              }
              grid::NDGrid grid_k = grid::NDGrid(k_sizes);
              while (!grid_k.nextPoint())
              {
                auto gridPoint_k = grid_k.getPoint();
                std::vector<int> gridPoint_jk(tree_s->size(), 0);
                for (size_t j = 0; j < gridPoint_j.size(); ++j)
                {
                  auto j_dim = j_dims.begin();
                  std::advance(j_dim, j);
                  int j_index = dim_to_pos[tree_index_s][*j_dim];
                  gridPoint_jk[j_index] = gridPoint_j[j];
                }
                for (size_t k = 0; k < gridPoint_k.size(); ++k)
                {
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
  for (size_t tree_index3 = 0; tree_index3 < variables.size(); ++tree_index3)
  {
    LeafGrid curr_gridLeaf;
    curr_gridLeaf.grid = grids[tree_index3];
    curr_gridLeaf.individuals = individuals[tree_index3];
    curr_gridLeaf.lim_list = lim_list;
    // If this tree exceeds the cap, attach a zero-valued matrix of the correct shape
    if (maxp_interaction > 0 && variables[tree_index3] != std::set<int>{0} && variables[tree_index3].size() > (size_t)maxp_interaction)
    {
      curr_gridLeaf.values = utils::Matrix<std::vector<double>>(grids[tree_index3].dimensions, std::vector<double>(value_size, 0));
    }
    else
    {
      curr_gridLeaf.values = values[tree_index3];
    }
    curr_family[variables[tree_index3]]->GridLeaves = curr_gridLeaf;
  }
}




void RandomPlantedForest::purify_2()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);

    // go through all variables of the component
    for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
    {
      std::vector<double> bounds;

      // go through trees of family
      for (const auto &curr_tree : curr_family)
      {

        // consider only relevant trees that have current dimension as variable
        if (!curr_tree.first.count(curr_dim))
          continue;

        // go through leaves of tree
        for (const auto &curr_leaf : curr_tree.second->leaves)
        {
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
        }
      }
      std::sort(bounds.begin(), bounds.end());
      bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
      lim_list[curr_dim - 1] = bounds;
    }

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids(curr_family.size() - 1);
    std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
    std::vector<std::set<int>> variables(curr_family.size() - 1);

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {

      if (curr_tree.first == std::set<int>{0})
        continue; // ignore null tree

      // fill space with dimensions
      std::vector<int> dimensions;
      for (const auto &dim : curr_tree.first)
      {
        dimensions.push_back(lim_list[dim - 1].size() - 1); // size - 1 ?
      }

      // setup grid for leaf indices
      auto grid = grid::NDGrid(dimensions);

      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
      values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
      variables[tree_index] = curr_tree.first;

      // fill grid points with individuals and values
      while (!grid.nextPoint())
      {

        std::vector<int> gridPoint = grid.getPoint();

        bool in_leaf = true;

        // go through sample points to sum up individuals
        for (const auto &feature_vec : X)
        {
          int dim_index = 0;
          in_leaf = true;
          for (const auto &dim : curr_tree.first)
          {
            double val = feature_vec[dim - 1];
            if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
            {
              in_leaf = false;
              break;
            }
            ++dim_index;
          }

          // consider individuals only if all in
          if (in_leaf)
            individuals[tree_index][gridPoint] += 1;
        }

        // go through leaves of tree to sum up values
        for (const auto &leaf : curr_tree.second->get_leaves())
        {

          in_leaf = true;
          int dim_index = 0;
          for (const auto &dim : curr_tree.first)
          {
            // consider values only if all in
            if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
            {
              in_leaf = false;
              break;
            }
            ++dim_index;
          }

          // sum up values
          if (in_leaf)
            values[tree_index][gridPoint] += leaf.value; // todo: multiclass
        }
      }

      ++tree_index;
    }

    // ------------- create new trees -------------

    // insert null tree
    grids.insert(grids.begin(), grid::NDGrid());
    values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
    variables.insert(variables.begin(), std::set<int>{0});

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (const auto &tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    auto keys = getKeys(curr_family);
    while (curr_max > 1)
    {

      // go through split dimensions of all trees
      for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          int dim_index = 0;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (!tree)
              {

                // get index of old and new tree
                auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
                curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));

                // remove matrix dimension of respective variable
                std::vector<int> matrix_dimensions = values[old_tree_index].dims;
                matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);

                // initialize data for new tree
                auto grid = grid::NDGrid(matrix_dimensions);
                grids.insert(grids.begin() + tree_index, grid);
                values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(0, value_size)));
                individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
                variables.insert(variables.begin() + tree_index, tree_dims);

                // fill individuals of new trees
                while (!grid.nextPoint())
                {

                  std::vector<int> gridPoint = grid.getPoint();
                  bool in_leaf = true;

                  // go through sample points to sum up individuals
                  for (const auto &feature_vec : X)
                  {
                    int dim_index = 0;
                    in_leaf = true;
                    for (const auto &dim : tree_dims)
                    {
                      double val = feature_vec[dim - 1];
                      if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
                        in_leaf = false;
                      ++dim_index;
                    }

                    // consider individuals only if all in
                    if (in_leaf)
                      individuals[tree_index][gridPoint] += 1;
                  }
                }
              }

              dim_index++;
            }
          }
        }
      }

      // update currently considered dimension size
      --curr_max;
    }

    // ------------- purify -------------

    // measure tolerance and number of iterations
    std::vector<double> tol(curr_family.size(), 1);
    int iter;

    // iterate backwards through tree family
    int curr_tree_index = curr_family.size() - 1;
    for (TreeFamily::reverse_iterator curr_tree = curr_family.rbegin(); curr_tree != curr_family.rend(); ++curr_tree)
    {
      iter = 0;
      std::set<int> curr_dims = curr_tree->second->get_split_dims();

      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;

      // repeat until tolerance small enough and (?) maximum number of iterations reached
      while ((tol[curr_tree_index] > 0.00000000001) && (iter < 100))
      {

        // go through feature dims
        int curr_dim_index = 0;
        for (const auto &feature_dim : curr_dims)
        {

          // get tree that has same variables as curr_tree minus j-variable
          std::set<int> tree_dims = curr_dims;
          tree_dims.erase(tree_dims.find(feature_dim));
          int tree_index = 0; // if tree not exist, set to null tree
          if (curr_family.find(tree_dims) != curr_family.end())
            tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims)) - 1;

          // update values
          if (grids[curr_tree_index].dimensions.size() == 1)
          { // one dimensional case

            int sum_ind = 0;
            std::vector<double> avg(value_size, 0);

            // get sum of individuals
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              sum_ind += individuals[curr_tree_index][tmp];
            }
            if (sum_ind == 0)
              continue;

            // calc avg
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              avg += (individuals[curr_tree_index][tmp] * values[curr_tree_index][tmp]) / sum_ind;
            }

            // update values of one dimensional and null tree
            for (int i = 0; i < values[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              values[curr_tree_index][tmp] -= avg;
            }
            std::vector<int> tmp{0};
            values[tree_index][tmp] += avg;
          }
          else
          { // higher dimensional case

            // setup new grid without dimension j
            std::vector<int> new_dimensions = grids[curr_tree_index].dimensions;
            int j_dim = new_dimensions[curr_dim_index];
            new_dimensions.erase(new_dimensions.begin() + curr_dim_index);
            grid::NDGrid grid = grid::NDGrid(new_dimensions);

            // go through values without dimension j
            while (!grid.nextPoint())
            {
              auto gridPoint = grid.getPoint();
              gridPoint.push_back(0);

              int sum_ind = 0;
              std::vector<double> avg(value_size, 0);

              // go through slice to sum up individuals
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // get sum of individuals
                sum_ind += individuals[curr_tree_index][gridPoint];
              }

              // go through slice to calc avg
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // calc avg
                avg += (individuals[curr_tree_index][gridPoint] * values[curr_tree_index][gridPoint]) / sum_ind;
              }

              // go through slice to update values
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // update values of current slice
                values[curr_tree_index][gridPoint] -= avg;
              }

              // update lower dimensional tree
              gridPoint.pop_back();
              values[tree_index][gridPoint] += avg;
            }
          }

          ++curr_dim_index;
        }

        // update tolerance
        if (variables[curr_tree_index].size() == 1)
        {
          tol[curr_tree_index] = 1; // todo
        }
        else
        {
          tol[curr_tree_index] = 1;
        }

        ++iter;
      }

      --curr_tree_index;
    }

    // ------------- attach to rpf class -------------

    // fill with new trees
    for (size_t tree_index = 0; tree_index < variables.size(); ++tree_index)
    {
      LeafGrid curr_gridLeaf;
      curr_gridLeaf.grid = grids[tree_index];
      curr_gridLeaf.individuals = individuals[tree_index];
      curr_gridLeaf.lim_list = lim_list;
      curr_gridLeaf.values = values[tree_index];
      curr_family[variables[tree_index]]->GridLeaves = curr_gridLeaf;
    }
  }

  purified = true;
}

