// Split-mode: leaves. Evaluates per-leaf candidate splits using cached
// per-leaf orders and prefix sums, with age-weighted candidate sampling.
#include "rpf.hpp"
#include "internal_utils.hpp"
#include <cmath>

using namespace rpf_utils;

namespace rpf_utils {
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
} // namespace rpf_utils

namespace {
// Each working row's rank in feature k among the family's rows (sc.rank[k]), computed when the feature is first
// searched: rows with equal values share a rank, a smaller value has a smaller rank. "Equal" is the comparator's
// (neither value below the other), so -0 and +0 share a rank, as they tie in a sort by value.
const std::vector<uint32_t> &family_rank(const std::vector<std::vector<double>> &X, int k, LeavesScratch &sc)
{
  if (sc.rank_owner != static_cast<const void *>(&X) || sc.rank.size() != X[0].size()) {
    sc.rank_owner = &X;
    sc.rank.assign(X[0].size(), std::vector<uint32_t>());
    sc.rank_max.assign(X[0].size(), 0);
  }
  std::vector<uint32_t> &rank = sc.rank[(size_t)k];
  if (rank.size() == X.size())
    return rank;
  const size_t n = X.size();
  std::vector<std::pair<double, uint32_t>> v(n);
  for (size_t i = 0; i < n; ++i)
    v[i] = {X[i][(size_t)k], (uint32_t)i};
  std::sort(v.begin(), v.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
  rank.resize(n);
  uint32_t r = 0;
  for (size_t j = 0; j < n; ++j) {
    if (j > 0 && v[j - 1].first < v[j].first)
      ++r;
    rank[v[j].second] = r;
  }
  sc.rank_max[(size_t)k] = r;
  return rank;
}

// The leaf's rows stably sorted by feature k and the first position of each distinct value, cached in the leaf
// until its rows change. The order is the one std::stable_sort gives the leaf's rows by X[row][k]: the rows are
// sorted by (the row's rank in the family, its position in the leaf), as integers. Leaves of up to 256 rows are
// sorted by comparison; larger ones by a radix sort on the rank (least significant digit first, which keeps the
// positions' order among equal ranks).
const LeafOrder &leaf_order(const std::vector<std::vector<double>> &X, Leaf &leaf, int k, LeavesScratch &sc)
{
  const size_t m = leaf.individuals.size();
  auto &cache = leaf.cache().leaf_order_cache;
  auto it = cache.find(k);
  if (it != cache.end() && it->second.rows.size() == m)
    return it->second;
  const std::vector<uint32_t> &rank = family_rank(X, k, sc);
  std::vector<uint64_t> &key = sc.sort_keys;
  key.resize(m);
  for (size_t i = 0; i < m; ++i)
    key[i] = ((uint64_t)rank[(size_t)leaf.individuals[i]] << 32) | (uint64_t)i;
  if (m <= 256) {
    std::sort(key.begin(), key.end());
  } else {
    const int bits = 11;
    const uint64_t mask = (1u << bits) - 1;
    std::vector<uint64_t> &tmp = sc.sort_tmp;
    tmp.resize(m);
    uint64_t *src = key.data(), *dst = tmp.data();
    uint32_t count[1u << bits];
    for (int shift = 32; shift < 64 && (sc.rank_max[(size_t)k] >> (shift - 32)) != 0; shift += bits) {
      std::fill(count, count + (1u << bits), 0u);
      for (size_t i = 0; i < m; ++i)
        ++count[(src[i] >> shift) & mask];
      uint32_t at = 0;
      for (uint32_t d = 0; d < (1u << bits); ++d) {
        const uint32_t c = count[d];
        count[d] = at;
        at += c;
      }
      for (size_t i = 0; i < m; ++i)
        dst[count[(src[i] >> shift) & mask]++] = src[i];
      std::swap(src, dst);
    }
    if (src != key.data())
      std::copy(src, src + m, key.data());
  }
  LeafOrder &lo = cache[k];
  lo.rows.resize(m);
  lo.first_pos.clear();
  uint64_t prev = 0;
  for (size_t i = 0; i < m; ++i) {
    lo.rows[i] = (int32_t)leaf.individuals[(size_t)(key[i] & 0xffffffffu)];
    const uint64_t r = key[i] >> 32;
    if (i == 0 || r != prev)
      lo.first_pos.push_back((int32_t)i);
    prev = r;
  }
  lo.first_pos.shrink_to_fit();
  return lo;
}

// min_rows_range() on a leaf's cached first positions
void min_rows_range_of(
    const std::vector<int32_t> &first_pos, size_t m, int leaf_size, int min_rows, int &lo, int &hi_excl
)
{
  const int n_distinct = (int)first_pos.size();
  // (positions are below m <= INT32_MAX: compared as 64-bit integers)
  const int64_t mr = (int64_t)min_rows, mm = (int64_t)m;
  const auto beg = first_pos.begin(), end = first_pos.end();
  const int a = (int)(std::lower_bound(beg, end, mr) - beg);
  const int b = mm >= mr ? (int)(std::upper_bound(beg, end, mm - mr) - beg) : 0;
  lo = std::max(leaf_size, a);
  hi_excl = std::min(n_distinct - leaf_size + 1, b);
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

// A pool entry's weight for the fast draw: exp(-decay * age), cached in the
// entry until its age changes, as draw_leaf_candidates() does. With
// rows_check, 0 for an entry whose rows_ok is 0, which the exact draw leaves
// out.
template <typename Cand> double leaf_weight(Cand &c, double decay, bool rows_check)
{
  if (rows_check && c.rows_ok == 0)
    return 0.0;
  if (!(c.w_age == c.age)) {
    c.w = std::exp(-decay * c.age);
    c.w_age = c.age;
  }
  return c.w;
}

// The fast draw: brings the sampler's tree up to date with the pool, then
// draws from it. Brought up to date here: the entries appended since the last
// draw and, with rows_check, the entries whose leaf changed since (fs.dirty,
// filled where create_tree_family() resets their rows_ok). An entry whose
// rows_ok is unknown (-1) gets feasible(entry) first, so that every entry's
// rows_ok is known at the draw, as before the exact draw. The drawn entries'
// new ages are set after the split, in calcOptimalSplit_leaves().
template <typename Cand, typename Feasible>
std::vector<size_t> draw_leaf_candidates_fast(
    std::vector<Cand> &pool, size_t n_candidates, double decay, LeavesSampler &fs, bool rows_check, Feasible feasible
)
{
  if (fs.owner != static_cast<const void *>(&pool) || fs.synced > pool.size())
    fs.reset(&pool);
  auto weight = [&](Cand &c) {
    if (rows_check && c.rows_ok < 0)
      c.rows_ok = feasible(c) ? 1 : 0;
    return leaf_weight(c, decay, rows_check);
  };
  for (size_t i : fs.dirty)
    if (i < fs.synced)
      fs.set(i, weight(pool[i]));
  fs.dirty.clear();
  for (size_t i = fs.synced; i < pool.size(); ++i)
    fs.push(weight(pool[i]));
  fs.synced = pool.size();
  std::vector<size_t> out;
  fs.draw(n_candidates, out);
  if (out.empty()) // no positive weight: uniformly, as the exact draw
    out = sample_weighted_indices_filtered(std::vector<double>(pool.size(), 0.0), n_candidates);
  return out;
}

// The lowest set bit of i: the number of entries a Fenwick node i holds.
inline size_t lowbit(size_t i)
{
  return i & (~i + 1);
}
} // namespace

void LeavesSampler::reset(const void *o)
{
  owner = o;
  synced = 0;
  w.clear();
  bit.clear();
  updates = 0;
  dirty.clear();
}

void LeavesSampler::rebuild()
{
  const size_t n = w.size();
  bit = w;
  for (size_t i = 1; i <= n; ++i) {
    const size_t j = i + lowbit(i);
    if (j <= n)
      bit[j - 1] += bit[i - 1];
  }
  updates = 0;
}

void LeavesSampler::set(size_t i, double nw)
{
  const double d = nw - w[i];
  w[i] = nw;
  if (d != 0.0)
    fenwick_add(bit, i + 1, d);
  if (++updates > 4 * w.size() + 1024) // bounds the rounding drift of the running sums
    rebuild();
}

void LeavesSampler::push(double nw)
{
  // node i (1-based) holds the entries (i - lowbit(i), i]: nw plus the nodes
  // i - 1, i - 2, i - 4, ... below lowbit(i)
  w.push_back(nw);
  const size_t i = bit.size() + 1;
  const size_t low = lowbit(i);
  double node = nw;
  for (size_t step = 1; step < low; step <<= 1)
    node += bit[i - step - 1];
  bit.push_back(node);
}

double LeavesSampler::total() const
{
  double s = 0.0;
  for (size_t i = bit.size(); i > 0; i -= lowbit(i))
    s += bit[i - 1];
  return s;
}

void LeavesSampler::draw(size_t k, std::vector<size_t> &out)
{
  taken.clear();
  const size_t start = out.size();
  for (size_t j = 0; j < k; ++j) {
    double tot = total();
    if (!(tot > 0.0)) {
      rebuild();
      tot = total();
      if (!(tot > 0.0))
        break;
    }
    size_t id = SIZE_MAX;
    for (int attempt = 0; attempt < 2 && id == SIZE_MAX; ++attempt) {
      const double r = (1.0 - rng_runif01()) * tot; // in (0, tot]
      const size_t i1 = fenwick_find_by_prefix(bit, r);
      if (i1 >= 1 && i1 <= w.size() && w[i1 - 1] > 0.0) {
        id = i1 - 1;
      } else { // rounding of the running sums: rebuild and draw again
        rebuild();
        tot = total();
        if (!(tot > 0.0))
          break;
      }
    }
    if (id == SIZE_MAX) { // rounding again: the last entry with a positive weight
      for (size_t i = w.size(); i > 0; --i)
        if (w[i - 1] > 0.0) {
          id = i - 1;
          break;
        }
      if (id == SIZE_MAX)
        break;
    }
    out.push_back(id);
    taken.push_back(w[id]);
    set(id, 0.0); // without replacement
  }
  for (size_t j = 0; j < taken.size(); ++j)
    set(out[start + j], taken[j]);
}

namespace rpf_utils {
std::vector<double> leaves_draw_counts(
    const std::vector<double> &ages,
    double decay,
    size_t n_candidates,
    size_t reps,
    bool exact,
    const std::vector<int> &rows_ok,
    uint64_t seed
)
{
  struct Entry {
    double age, w = 0.0, w_age = std::numeric_limits<double>::quiet_NaN();
    signed char rows_ok = -1;
  };
  const bool rows_check = !rows_ok.empty();
  if (rows_check && rows_ok.size() != ages.size())
    throw std::invalid_argument("leaves_draw_counts: rows_ok must be empty or as long as ages.");
  std::vector<Entry> pool;
  for (size_t i = 0; i < ages.size(); ++i) {
    pool.push_back(Entry{ages[i]});
    if (rows_check)
      pool.back().rows_ok = rows_ok[i] != 0 ? 1 : 0;
  }
  auto feasible = [](const Entry &) { return true; }; // every rows_ok is known
  std::mt19937_64 rng(seed);
  struct RngGuard {
    std::mt19937_64 *old;
    ~RngGuard() { swap_tls_rng(old); }
  } guard{swap_tls_rng(&rng)};
  LeavesSampler fs;
  std::vector<std::pair<double, size_t>> keys;
  std::vector<double> counts(pool.size(), 0.0);
  for (size_t r = 0; r < reps; ++r) {
    const std::vector<size_t> drawn =
        exact ? draw_leaf_candidates(pool, n_candidates, decay, keys, rows_check)
              : draw_leaf_candidates_fast(pool, n_candidates, decay, fs, rows_check, feasible);
    for (size_t i : drawn)
      counts[i] += 1.0;
  }
  return counts;
}
} // namespace rpf_utils

// This round's candidates of `leaves`, for regression and classification
// alike (the same entries from the same pool and random numbers). With
// min_leaf_rows, pool entries whose leaf has no threshold with min_rows rows
// on each side are not drawn (checked once per leaf and feature), by either
// draw.
std::vector<size_t> RandomPlantedForest::draw_leaves_candidates(
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    size_t n_candidates,
    LeavesScratch &sc
)
{
  const int min_rows = this->min_leaf_rows_;
  const bool rows_check = min_rows > 1;
  auto feasible = [&](const SplitCandidate &c) {
    return c.tree && c.leaf_idx < c.tree->leaves.size() &&
           min_rows_feasible(X, c.tree->leaves[c.leaf_idx], c.dim - 1, this->n_leaves[c.dim - 1], min_rows);
  };
  std::vector<size_t> sample_idxs;
  if (this->deterministic) {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i)
      sample_idxs.push_back(i);
  } else if (this->exact_draw_) {
    if (rows_check)
      for (auto &c : possible_splits)
        if (c.rows_ok < 0)
          c.rows_ok = feasible(c) ? 1 : 0;
    sample_idxs = draw_leaf_candidates(possible_splits, n_candidates, this->split_decay_rate_, sc.keys, rows_check);
  } else {
    sample_idxs = draw_leaf_candidates_fast(
        possible_splits, n_candidates, this->split_decay_rate_, sc.sampler, rows_check, feasible
    );
  }
  return sample_idxs;
}

// The fast draw's tree: the drawn entries' weights at their new ages.
void RandomPlantedForest::leaves_sampler_after_aging(
    std::vector<SplitCandidate> &possible_splits, const std::vector<size_t> &sample_idxs, LeavesScratch &sc
)
{
  if (this->deterministic || this->exact_draw_)
    return;
  const bool rows_check = this->min_leaf_rows_ > 1;
  for (size_t idx : sample_idxs)
    sc.sampler.set(idx, leaf_weight(possible_splits[idx], this->split_decay_rate_, rows_check));
}

Split RandomPlantedForest::calcOptimalSplit_leaves(
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

  if (possible_splits.empty())
    return min_split;

  // Scratch buffers: the tree family's (create_tree_family), else local ones
  LeavesScratch local_scratch;
  LeavesScratch &sc = tls_leaves_scratch_ptr ? *tls_leaves_scratch_ptr : local_scratch;

  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(this->t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(this->max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  const int min_rows = this->min_leaf_rows_;
  const bool rows_check = min_rows > 1;
  const size_t scan_T = loss_fn ? loss_fn->scan_terms() : 0;
  std::vector<size_t> sample_idxs = draw_leaves_candidates(X, possible_splits, n_candidates, sc);

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

    // The leaf's rows sorted by the feature, and the first position of each distinct value (cached in the leaf)
    const LeafOrder &lord = leaf_order(X, *leafPtr, k, sc);
    const int32_t *rows = lord.rows.data();
    const std::vector<int32_t> &first_pos = lord.first_pos;
    const size_t n_unique = first_pos.size();

    if (n_unique < 2 * static_cast<size_t>(leaf_size))
      continue;

    std::vector<int> &samples = sc.samples;
    int left = leaf_size;
    int right_exclusive = (int)n_unique - leaf_size + 1;
    if (rows_check) { // only thresholds with at least min_rows rows on each side
      min_rows_range_of(first_pos, m, leaf_size, min_rows, left, right_exclusive);
      if (left >= right_exclusive)
        continue;
    }
    if (this->deterministic)
      samples = compute_even_spread_indices(left, right_exclusive, (size_t)this->split_try);
    else
      sample_unique_ints_uniform_into(left, right_exclusive, (size_t)this->split_try, samples, sc.samples_all);

    // The loss's running sums (FamilyLoss::scan_terms): the prefix sums of its row terms along the sorted rows,
    // [i * T + q], and the leaf's own value; a threshold then costs O(T), whatever the leaf's rows
    if (scan_T > 0) {
      std::vector<double> &tp = sc.prefix;
      tp.resize(scan_T * m);
      loss_fn->prefix_terms_rows(rows, m, tp.data());
      const double *tot = tp.data() + (m - 1) * scan_T;
      const double vp = loss_fn->region_value(tot, m);
      sc.total.resize(scan_T);
      for (size_t si = 0; si < samples.size(); ++si) {
        const size_t uidx = static_cast<size_t>(samples[si]);
        if (uidx >= n_unique)
          continue;
        const size_t pos = (size_t)first_pos[uidx];
        const double sp = X[(size_t)rows[pos]][(size_t)k]; // the distinct value: the first of its run
        if (pos == 0 || pos >= m)
          continue;
        if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
          continue;
        if (pos < static_cast<size_t>(min_rows) || (m - pos) < static_cast<size_t>(min_rows))
          continue;
        const double *s = tp.data() + (pos - 1) * scan_T;
        for (size_t q = 0; q < scan_T; ++q)
          sc.total[q] = tot[q] - s[q];
        const double loss = loss_fn->scan_score(s, sc.total.data(), pos, m - pos, vp);
        if (loss < min_split.min_sum) {
          min_split.min_sum = loss;
          min_split.tree_index = it->tree;
          min_split.leaf_index = leafPtr;
          min_split.split_coordinate = k + 1;
          min_split.split_point = sp;
          best_idx = (int)idx;
        }
      }
      continue;
    }

    // Prefix sums once per candidate evaluation, into a reused buffer [p * m + i]
    std::vector<double> &prefix = sc.prefix;
    std::vector<double> &total = sc.total;
    prefix.resize(this->value_size * m);
    total.assign(this->value_size, 0.0);
    for (size_t p = 0; p < this->value_size; ++p) {
      double acc = 0.0;
      double *pp = prefix.data() + p * m;
      for (size_t i = 0; i < m; ++i) {
        acc += Y[(size_t)rows[i]][p];
        pp[i] = acc;
      }
      total[p] = pp[m - 1];
    }

    for (size_t si = 0; si < samples.size(); ++si) {
      const size_t uidx = static_cast<size_t>(samples[si]);
      if (uidx >= n_unique)
        continue;
      const size_t pos = (size_t)first_pos[uidx];
      const double sp = X[(size_t)rows[pos]][(size_t)k]; // the distinct value: the first of its run
      if (pos == 0 || pos >= m)
        continue;
      if (pos < static_cast<size_t>(leaf_size) || (m - pos) < static_cast<size_t>(leaf_size))
        continue;
      if (pos < static_cast<size_t>(min_rows) || (m - pos) < static_cast<size_t>(min_rows))
        continue;

      double loss = 0.0;
      if (loss_fn) { // the loss's criterion (classification), from the same sums
        curr_split.sum_s.resize(this->value_size);
        curr_split.sum_b.resize(this->value_size);
        for (size_t p = 0; p < this->value_size; ++p) {
          curr_split.sum_s[p] = prefix[p * m + pos - 1];
          curr_split.sum_b[p] = total[p] - curr_split.sum_s[p];
        }
        loss = loss_fn->score(curr_split, *leafPtr, k, sp, X);
      } else {
        for (size_t p = 0; p < this->value_size; ++p) {
          const double sum_s_base = prefix[p * m + pos - 1];
          const double sum_b_base = total[p] - sum_s_base;
          loss -= (sum_s_base * sum_s_base) / static_cast<double>(pos);
          loss -= (sum_b_base * sum_b_base) / static_cast<double>(m - pos);
        }
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
  leaves_sampler_after_aging(possible_splits, sample_idxs, sc);
  if (loss_fn)
    loss_fn->finalize(min_split, X);
  else
    finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}
