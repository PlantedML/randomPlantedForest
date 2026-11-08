// Split-mode: leaves. Evaluates per-leaf candidate splits using cached
// per-leaf orders and prefix sums, with age-weighted candidate sampling.
#include "rpf.hpp"
#include "internal_utils.hpp"
#include <cmath>

using namespace rpf_utils;

namespace {
// The leaf's cached stable order and sorted values of feature k, by reference.
// Computed and cached as ensure_order_and_sorted_vals_for_leaf() does when
// missing (the same caches, the same conditions), without copying them out.
// The references stay valid while no entry of these two maps is erased
// (unordered_map nodes do not move).
void leaf_order_refs(
    const std::vector<std::vector<double>> &X,
    Leaf &leaf,
    int k,
    const std::vector<size_t> *&order,
    const std::vector<double> *&sorted_vals
)
{
  const size_t m = leaf.individuals.size();
  auto it = leaf.order_cache.find(k);
  if (it == leaf.order_cache.end() || it->second.size() != m) {
    std::vector<size_t> &o = leaf.order_cache[k];
    o.resize(m);
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) {
      return X[leaf.individuals[a]][k] < X[leaf.individuals[b]][k];
    });
    it = leaf.order_cache.find(k);
  }
  order = &it->second;
  auto is = leaf.sorted_vals_cache.find(k);
  if (is == leaf.sorted_vals_cache.end() || is->second.size() != m) {
    std::vector<double> &sv = leaf.sorted_vals_cache[k];
    sv.resize(m);
    for (size_t i = 0; i < m; ++i)
      sv[i] = X[leaf.individuals[(*order)[i]]][k];
    is = leaf.sorted_vals_cache.find(k);
  }
  sorted_vals = &is->second;
}

// The age-weighted draw of n_candidates pool entries without replacement, as
// sample_weighted_indices_filtered() on the weights exp(-decay * age), fused
// into one pass: the same weights (each entry's exp cached until its age
// changes), the same random draws in the same order (one per entry with a
// positive weight), the same keys -log(u) / w in the same order and the same
// nth_element, so the same sample in the same order.
//
// With rows_check, entries whose rows_ok is 0 are left out, as an entry with
// weight 0 is (no random draw).
template <typename Cand>
std::vector<size_t> draw_leaf_candidates(
    std::vector<Cand> &pool,
    size_t n_candidates,
    double decay,
    std::vector<std::pair<double, size_t>> &keys,
    bool rows_check
)
{
  keys.clear();
  for (size_t i = 0; i < pool.size(); ++i) {
    Cand &c = pool[i];
    if (rows_check && c.rows_ok == 0)
      continue;
    if (!(c.w_age == c.age)) {
      c.w = std::exp(-decay * c.age);
      c.w_age = c.age;
    }
    if (c.w > 0.0) {
      double u = rng_runif01();
      if (u <= 0.0)
        u = std::numeric_limits<double>::min();
      keys.emplace_back(-std::log(u) / c.w, i);
    }
  }
  if (keys.empty()) // no positive weight: uniformly, as sample_weighted_indices_filtered()
    return sample_weighted_indices_filtered(std::vector<double>(pool.size(), 0.0), n_candidates);
  const size_t k = std::min(n_candidates, keys.size());
  if (k < keys.size())
    std::nth_element(keys.begin(), keys.begin() + k, keys.end(), [](const auto &a, const auto &b) {
      return a.first < b.first;
    });
  std::vector<size_t> out;
  out.reserve(k);
  for (size_t j = 0; j < k; ++j)
    out.push_back(keys[j].second);
  return out;
}

// The distinct-value indices u of a leaf's sorted values (first positions
// first_pos of the m rows) whose threshold leaves at least leaf_size distinct
// values and at least min_rows rows on each side: [lo, hi_excl). first_pos is
// increasing, so these u are contiguous. With min_rows = 1 this is
// [leaf_size, distinct - leaf_size + 1), the range without a row limit.
void min_rows_range(const std::vector<size_t> &first_pos, size_t m, int leaf_size, int min_rows, int &lo, int &hi_excl)
{
  const int n_distinct = (int)first_pos.size();
  const size_t mr = (size_t)min_rows;
  const int a = (int)(std::lower_bound(first_pos.begin(), first_pos.end(), mr) - first_pos.begin());
  const int b = (m >= mr) ? (int)(std::upper_bound(first_pos.begin(), first_pos.end(), m - mr) - first_pos.begin()) : 0;
  lo = std::max(leaf_size, a);
  hi_excl = std::min(n_distinct - leaf_size + 1, b);
}

// Whether the leaf has a threshold on feature k with at least min_rows rows and
// at least leaf_size distinct values on each side.
bool min_rows_feasible(const std::vector<std::vector<double>> &X, const Leaf &leaf, int k, int leaf_size, int min_rows)
{
  const size_t m = leaf.individuals.size();
  if (m < 2 * (size_t)min_rows || m < 2 * (size_t)leaf_size)
    return false;
  std::vector<double> v(m);
  for (size_t i = 0; i < m; ++i)
    v[i] = X[leaf.individuals[i]][k];
  if (leaf_size == 1) { // a threshold exists iff the min_rows-th smallest value is below the min_rows-th largest
    std::nth_element(v.begin(), v.begin() + (min_rows - 1), v.end());
    const double a = v[min_rows - 1];
    std::nth_element(v.begin(), v.begin() + (m - min_rows), v.end());
    return a < v[m - min_rows];
  }
  std::sort(v.begin(), v.end());
  std::vector<size_t> first_pos(1, 0);
  for (size_t i = 1; i < m; ++i)
    if (v[i] != v[i - 1])
      first_pos.push_back(i);
  int lo, hi_excl;
  min_rows_range(first_pos, m, leaf_size, min_rows, lo, hi_excl);
  return lo < hi_excl;
}
} // namespace

Split RandomPlantedForest::calcOptimalSplit_leaves(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  if (possible_splits.empty())
    return min_split;

  // Scratch buffers: the tree family's (create_tree_family), else local ones
  LeavesScratch local_scratch;
  LeavesScratch &sc = tls_leaves_scratch_ptr ? *tls_leaves_scratch_ptr : local_scratch;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  // min_leaf_rows: pool entries whose leaf has no threshold with min_rows rows
  // on each side are not drawn (checked once per leaf and feature)
  const int min_rows = this->min_leaf_rows_;
  const bool rows_check = min_rows > 1;
  std::vector<size_t> sample_idxs;
  if (this->deterministic) {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i)
      sample_idxs.push_back(i);
  } else {
    if (rows_check)
      for (auto &c : possible_splits)
        if (c.rows_ok < 0) {
          const bool ok =
              c.tree && c.leaf_idx < c.tree->leaves.size() &&
              min_rows_feasible(X, c.tree->leaves[c.leaf_idx], c.dim - 1, this->n_leaves[c.dim - 1], min_rows);
          c.rows_ok = ok ? 1 : 0;
        }
    sample_idxs = draw_leaf_candidates(possible_splits, n_candidates, this->split_decay_rate_, sc.keys, rows_check);
  }

  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    if (!it->tree || it->leaf_idx >= it->tree->leaves.size())
      continue;
    Leaf *leafPtr = &it->tree->leaves[it->leaf_idx];

    const int leaf_size = this->n_leaves[k];
    const size_t m = leafPtr->individuals.size();
    if (m == 0)
      continue;
    // Quick infeasibility check: cannot split if fewer than 2*leaf_size individuals
    if (m < static_cast<size_t>(2 * leaf_size) || m < static_cast<size_t>(2 * min_rows))
      continue;

    // The leaf's cached order and sorted values by reference; its distinct
    // values and their first positions in one pass, into reused buffers
    const std::vector<size_t> *order_p;
    const std::vector<double> *sorted_p;
    leaf_order_refs(X, *leafPtr, k, order_p, sorted_p);
    const std::vector<size_t> &order = *order_p;
    const std::vector<double> &sorted_vals = *sorted_p;
    std::vector<double> &unique = sc.unique;
    std::vector<size_t> &first_pos = sc.first_pos;
    unique.clear();
    first_pos.clear();
    if (!sorted_vals.empty()) {
      unique.push_back(sorted_vals[0]);
      first_pos.push_back(0);
      for (size_t i = 1; i < sorted_vals.size(); ++i) {
        if (sorted_vals[i] != unique.back())
          unique.push_back(sorted_vals[i]);
        if (sorted_vals[i] != sorted_vals[i - 1])
          first_pos.push_back(i);
      }
    }

    if (unique.size() < 2 * static_cast<size_t>(leaf_size))
      continue;

    std::vector<int> samples;
    int left = leaf_size;
    int right_exclusive = (int)unique.size() - leaf_size + 1;
    if (rows_check) { // only thresholds with at least min_rows rows on each side
      min_rows_range(first_pos, m, leaf_size, min_rows, left, right_exclusive);
      if (left >= right_exclusive)
        continue;
    }
    samples = this->deterministic ? compute_even_spread_indices(left, right_exclusive, (size_t)this->split_try)
                                  : sample_unique_ints_uniform_R(left, right_exclusive, (size_t)this->split_try);

    // Prefix sums once per candidate evaluation, into a reused buffer [p * m + i]
    std::vector<double> &prefix = sc.prefix;
    std::vector<double> &total = sc.total;
    prefix.resize(this->value_size * m);
    total.assign(this->value_size, 0.0);
    for (size_t p = 0; p < this->value_size; ++p) {
      double acc = 0.0;
      double *pp = prefix.data() + p * m;
      for (size_t i = 0; i < m; ++i) {
        acc += Y[leafPtr->individuals[order[i]]][p];
        pp[i] = acc;
      }
      total[p] = pp[m - 1];
    }

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
      if (pos < static_cast<size_t>(min_rows) || (m - pos) < static_cast<size_t>(min_rows))
        continue;

      double loss = 0.0;
      for (size_t p = 0; p < this->value_size; ++p) {
        const double sum_s_base = prefix[p * m + pos - 1];
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
          const double sum_s_base = prefix[p * m + pos - 1];
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
