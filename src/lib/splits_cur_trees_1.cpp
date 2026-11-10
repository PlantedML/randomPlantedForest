// Split-mode: cur_trees_1. Samples feasible leaves proportionally to their
// number of candidate thresholds, then evaluates a single threshold.
//
// How the search is organised (the fit is that of a row-by-row search, bit for bit):
// - All tries of a split search (split_try per candidate) are drawn first, with the same random draws in the same
//   order, and then evaluated; the first strict minimum in candidate and try order wins.
// - A try's pass over its leaf's rows is branch-free and compares the rows' ranks (cur_trees_1.hpp) instead of their
//   values. Four tries run in lockstep (rpf_ct1::lane_sums), so that the additions of their sums overlap. Each sum is
//   still accumulated row by row in the leaf's row order, so the losses, the winner and the leaf values are the same.
// - Leaf weights and distinct values come from the family's ranks: no sorted order is kept per leaf and feature.
#include "rpf.hpp"
#include "internal_utils.hpp"
#include "cur_trees_1.hpp"

using namespace rpf_utils;

namespace rpf_ct1 {

uint32_t FamilyRanks::next_epoch()
{
  if (++epoch == 0) {
    std::fill(stamp.begin(), stamp.end(), 0u);
    epoch = 1;
  }
  return epoch;
}

FamilyRanks &family_ranks(std::shared_ptr<FamilyRanks> &slot, const std::vector<std::vector<double>> &X)
{
  if (slot && slot->X == &X && slot->n == X.size())
    return *slot;
  auto fr = std::make_shared<FamilyRanks>();
  fr->X = &X;
  fr->n = X.size();
  fr->p = X.empty() ? 0 : X[0].size();
  fr->rank.assign(fr->p, std::vector<int32_t>(fr->n, 0));
  fr->value.assign(fr->p, std::vector<double>());
  std::vector<int32_t> idx(fr->n);
  size_t max_distinct = 1;
  for (size_t k = 0; k < fr->p; ++k) {
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int32_t a, int32_t b) { return X[(size_t)a][k] < X[(size_t)b][k]; });
    auto &rk = fr->rank[k];
    auto &val = fr->value[k];
    for (size_t i = 0; i < fr->n; ++i) {
      const double v = X[(size_t)idx[i]][k];
      if (val.empty() || val.back() < v)
        val.push_back(v);
      rk[(size_t)idx[i]] = (int32_t)val.size() - 1;
    }
    val.shrink_to_fit();
    max_distinct = std::max(max_distinct, val.size());
  }
  fr->stamp.assign(max_distinct, 0u);
  fr->epoch = 0;
  slot = std::move(fr);
  return *slot;
}

size_t unique_count(FamilyRanks &fr, Leaf &leaf, int k)
{
  auto it = leaf.cache().unique_count_cache.find(k);
  if (it != leaf.cache().unique_count_cache.end())
    return it->second;
  const int32_t *rk = fr.rank[(size_t)k].data();
  uint32_t *st = fr.stamp.data();
  const uint32_t e = fr.next_epoch();
  size_t c = 0;
  for (int ind : leaf.individuals) {
    const int32_t r = rk[ind];
    c += (st[r] != e) ? 1 : 0;
    st[r] = e;
  }
  leaf.cache().unique_count_cache[k] = c;
  return c;
}

double leaf_weight(FamilyRanks &fr, Leaf &leaf, int k, int leaf_size)
{
  const size_t unique = unique_count(fr, leaf, k);
  const long width_unique = (long)unique - 2L * (long)leaf_size + 1L;
  return (width_unique > 0L) ? static_cast<double>(width_unique) : 0.0;
}

const std::vector<double> &unique_values(FamilyRanks &fr, const std::vector<std::vector<double>> &X, Leaf &leaf, int k)
{
  auto it = leaf.cache().unique_vals_cache.find(k);
  if (it != leaf.cache().unique_vals_cache.end())
    return it->second;
  const int32_t *rk = fr.rank[(size_t)k].data();
  uint32_t *st = fr.stamp.data();
  const uint32_t e = fr.next_epoch();
  auto &buf = fr.buf;
  buf.resize(leaf.individuals.size());
  size_t c = 0;
  for (int ind : leaf.individuals) {
    const int32_t r = rk[ind];
    buf[c] = r;
    c += (st[r] != e) ? 1 : 0;
    st[r] = e;
  }
  buf.resize(c);
  std::sort(buf.begin(), buf.end());
  const auto &val = fr.value[(size_t)k];
  std::vector<double> u(buf.size());
  for (size_t i = 0; i < buf.size(); ++i) {
    double v = val[(size_t)buf[i]];
    // -0 and +0 share a rank: the run is represented by its first row in the leaf's order, as a stable sort gives
    if (v == 0.0)
      for (int ind : leaf.individuals)
        if (X[(size_t)ind][(size_t)k] == 0.0) {
          v = X[(size_t)ind][(size_t)k];
          break;
        }
    u[i] = v;
  }
  leaf.cache().unique_count_cache[k] = u.size();
  auto &slot = leaf.cache().unique_vals_cache[k];
  slot = std::move(u);
  return slot;
}

int32_t rank_of(const FamilyRanks &fr, int k, double v)
{
  // std::lower_bound without branches: the first i with val[i] >= v
  const auto &val = fr.value[(size_t)k];
  if (val.empty())
    return 0;
  const double *base = val.data();
  size_t n = val.size();
  while (n > 1) {
    const size_t half = n / 2;
    base = (base[half] < v) ? base + half : base;
    n -= half;
  }
  return (int32_t)((base - val.data()) + ((*base < v) ? 1 : 0));
}

namespace {
// L lanes in lockstep over rows [i0, i1) of each. A row adds its y to its lane's "below" sum if its rank is below the
// lane's threshold, and to the "not below" sum otherwise; the other sum gets + 0.0, which leaves it unchanged (a sum
// starts at +0 and so is never -0). So each sum is the row-by-row sum over its rows in the leaf's order, as a pass
// with a branch per row gives; the lanes only interleave independent sums.
template <size_t L> void lanes_step(Lane *const *ln, size_t i0, size_t i1, const std::vector<std::vector<double>> &Y)
{
  const int *rows[L];
  const int32_t *rk[L];
  int32_t thr[L];
  double S[L], B[L];
  size_t N[L];
  for (size_t j = 0; j < L; ++j) {
    rows[j] = ln[j]->rows->data();
    rk[j] = ln[j]->rk;
    thr[j] = ln[j]->thr;
    S[j] = ln[j]->s;
    B[j] = ln[j]->b;
    N[j] = ln[j]->ns;
  }
  for (size_t i = i0; i < i1; ++i) {
    for (size_t j = 0; j < L; ++j) {
      const int ind = rows[j][i];
      const double y = Y[(size_t)ind][0];
      const bool l = rk[j][ind] < thr[j];
      S[j] += l ? y : 0.0;
      B[j] += l ? 0.0 : y;
      N[j] += l ? 1 : 0;
    }
  }
  for (size_t j = 0; j < L; ++j) {
    ln[j]->s = S[j];
    ln[j]->b = B[j];
    ln[j]->ns = N[j];
  }
}
} // namespace

void lane_sums(std::vector<Lane> &lanes, const std::vector<std::vector<double>> &Y)
{
  // Longest first, in groups of 4 of about the same length; within a group, when the shortest lane ends, the others
  // go on with one lane fewer. The grouping changes no lane's sums, so the order of equal lengths does not matter.
  std::vector<Lane *> ord(lanes.size());
  for (size_t i = 0; i < lanes.size(); ++i) {
    lanes[i].ns = 0;
    lanes[i].s = 0.0;
    lanes[i].b = 0.0;
    ord[i] = &lanes[i];
  }
  std::sort(ord.begin(), ord.end(), [](const Lane *a, const Lane *b) { return a->rows->size() > b->rows->size(); });
  for (size_t g = 0; g < ord.size(); g += 4) {
    const size_t L = std::min<size_t>(4, ord.size() - g);
    Lane *grp[4];
    for (size_t j = 0; j < L; ++j)
      grp[j] = ord[g + L - 1 - j]; // ascending length
    size_t from = 0;
    for (size_t j = 0; j < L; ++j) {
      const size_t to = grp[j]->rows->size();
      if (to > from) {
        switch (L - j) {
        case 4:
          lanes_step<4>(grp + j, from, to, Y);
          break;
        case 3:
          lanes_step<3>(grp + j, from, to, Y);
          break;
        case 2:
          lanes_step<2>(grp + j, from, to, Y);
          break;
        default:
          lanes_step<1>(grp + j, from, to, Y);
          break;
        }
        from = to;
      }
    }
  }
}

} // namespace rpf_ct1

namespace {
// One threshold of one leaf: the row-by-row pass, for responses with several columns.
void sums_rowwise(
    const std::vector<std::vector<double>> &X,
    const std::vector<std::vector<double>> &Y,
    const Leaf &leaf,
    int k,
    double sp,
    size_t value_size,
    size_t &ns,
    size_t &nb,
    std::vector<double> &sum_s,
    std::vector<double> &sum_b
)
{
  ns = 0;
  nb = 0;
  sum_s.assign(value_size, 0.0);
  sum_b.assign(value_size, 0.0);
  for (int ind : leaf.individuals) {
    if (X[ind][k] < sp) {
      ++ns;
      for (size_t p = 0; p < value_size; ++p)
        sum_s[p] += Y[ind][p];
    } else {
      ++nb;
      for (size_t p = 0; p < value_size; ++p)
        sum_b[p] += Y[ind][p];
    }
  }
}
} // namespace

Split RandomPlantedForest::calcOptimalSplit_curTrees1(
    const std::vector<std::vector<double>> &Y,
    const std::vector<std::vector<double>> &X,
    std::vector<SplitCandidate> &possible_splits,
    TreeFamily &curr_family
)
{
  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;

  unsigned int raw = (unsigned int)std::ceil(this->t_try * possible_splits.size());
  unsigned int upper =
      std::min<unsigned int>((unsigned int)this->max_candidates_, (unsigned int)possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw, upper));
  std::vector<double> weights(possible_splits.size());

  for (size_t i = 0; i < possible_splits.size(); ++i)
    weights[i] = std::exp(-this->split_decay_rate_ * possible_splits[i].age);

  size_t positive_count_ = 0;
  for (double w : weights)
    if (w > 0.0)
      ++positive_count_;

  if (positive_count_ == 0) {
    n_candidates = 1;
  } else {
    if (n_candidates > positive_count_)
      n_candidates = static_cast<unsigned int>(positive_count_);
  }

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

  // The family's ranks of its working rows, held by its null tree
  auto itNull = curr_family.find(std::set<int>{0});
  if (itNull == curr_family.end())
    throw std::logic_error("cur_trees_1: the tree family has no null tree.");
  rpf_ct1::FamilyRanks &fr = rpf_ct1::family_ranks(itNull->second->ct1_ranks, X);
  const bool one_column = (this->value_size == 1);

  // A drawn (or, if deterministic, listed) threshold sp of a leaf in feature k, for candidate idx, and its sums
  struct Try {
    std::shared_ptr<DecisionTree> tree;
    Leaf *leaf;
    int k;
    size_t idx;
    double sp;
    int32_t r_sp;
    size_t ns, nb;
    double s, b;
    std::vector<double> sum_s, sum_b;
  };
  std::vector<Try> tries;
  std::vector<rpf_ct1::Lane> lanes;
  // The sums of tries[from..): with one column in lanes (rpf_ct1::lane_sums), with several row by row
  auto evaluate = [&](size_t from) {
    if (one_column) {
      lanes.resize(tries.size() - from);
      for (size_t i = from; i < tries.size(); ++i) {
        const Try &t = tries[i];
        lanes[i - from] = rpf_ct1::Lane{&t.leaf->individuals, fr.rank[(size_t)t.k].data(), t.r_sp, 0, 0.0, 0.0};
      }
      rpf_ct1::lane_sums(lanes, Y);
      for (size_t i = from; i < tries.size(); ++i) {
        Try &t = tries[i];
        const rpf_ct1::Lane &l = lanes[i - from];
        t.ns = l.ns;
        t.nb = t.leaf->individuals.size() - l.ns;
        t.s = l.s;
        t.b = l.b;
      }
    } else {
      for (size_t i = from; i < tries.size(); ++i) {
        Try &t = tries[i];
        sums_rowwise(X, Y, *t.leaf, t.k, t.sp, this->value_size, t.ns, t.nb, t.sum_s, t.sum_b);
      }
    }
  };
  // Compare tries[from..) in order (candidate order, then try order) with the best so far: the first strict minimum
  // wins
  int best_idx = -1;
  auto compare = [&](size_t from) {
    for (size_t i = from; i < tries.size(); ++i) {
      const Try &t = tries[i];
      if (t.ns == 0 || t.nb == 0)
        continue;
      double loss = 0.0;
      if (one_column) {
        loss -= (t.s * t.s) / (double)t.ns;
        loss -= (t.b * t.b) / (double)t.nb;
      } else {
        for (size_t p = 0; p < this->value_size; ++p) {
          loss -= (t.sum_s[p] * t.sum_s[p]) / (double)t.ns;
          loss -= (t.sum_b[p] * t.sum_b[p]) / (double)t.nb;
        }
      }
      if (loss < min_split.min_sum) {
        min_split.min_sum = loss;
        min_split.tree_index = t.tree;
        min_split.leaf_index = t.leaf;
        min_split.split_coordinate = t.k + 1;
        min_split.split_point = t.sp;
        best_idx = (int)t.idx;
        if (one_column) {
          min_split.sum_s.assign(1, t.s);
          min_split.sum_b.assign(1, t.b);
        } else {
          min_split.sum_s = t.sum_s;
          min_split.sum_b = t.sum_b;
        }
      }
    }
  };

  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin();
    std::advance(it, idx);
    int k = it->dim - 1;
    int leaf_size = this->n_leaves[k];
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

    if (!this->deterministic) {
      auto ensure_weights_cache = [&](const std::shared_ptr<DecisionTree> &tree, int kdim) {
        // Lazy-size vectors to feature_size once
        if ((int)tree->weights_epoch_by_dim_v.size() < this->feature_size) {
          tree->weights_epoch_by_dim_v.assign((size_t)this->feature_size, -1);
          tree->fenwick_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
          tree->leaf_weights_by_dim_v.assign((size_t)this->feature_size, std::vector<double>());
          tree->weights_total_by_dim_v.assign((size_t)this->feature_size, 0.0);
        }
        // Recompute BIT if epoch mismatches or size changed
        bool need = true;
        if (tree->weights_epoch_by_dim_v[(size_t)kdim] == tree->weights_epoch) {
          if (tree->fenwick_by_dim_v[(size_t)kdim].size() == tree->leaves.size())
            need = false;
        }
        if (!need)
          return; // cache fresh
        const size_t L = tree->leaves.size();
        std::vector<double> bit(L, 0.0), wts(L, 0.0);
        double total = 0.0;
        for (size_t li = 0; li < L; ++li) {
          // Weight = number of unique thresholds that respect min leaf size
          const double w = rpf_ct1::leaf_weight(fr, tree->leaves[li], kdim, leaf_size);
          wts[li] = w;
          total += w;
          if (w != 0.0)
            rpf_utils::fenwick_add(bit, li + 1, w);
        }
        tree->fenwick_by_dim_v[(size_t)kdim] = std::move(bit);
        tree->leaf_weights_by_dim_v[(size_t)kdim] = std::move(wts);
        tree->weights_total_by_dim_v[(size_t)kdim] = total;
        tree->weights_epoch_by_dim_v[(size_t)kdim] = tree->weights_epoch;
      };

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
        ensure_weights_cache(src_tree, k);
        double tot = src_tree->weights_total_by_dim_v[(size_t)k];
        if (tot <= 0.0)
          continue;
        src_info.push_back({src_tree, tot});
        total_all += tot;
      }
      if (src_info.empty() || total_all <= 0.0)
        continue;

      // Draw the tries; they are evaluated with the other candidates' tries below (evaluating draws nothing)
      for (size_t t = 0; t < (size_t)this->split_try; ++t) {
        // Sample a source tree proportionally to its total weight
        double r_src = rng_runif(0.0, total_all);
        size_t si = 0;
        while (si + 1 < src_info.size() && r_src > src_info[si].total) {
          r_src -= src_info[si].total;
          ++si;
        }
        auto &sel_tree = src_info[si].tree;

        // Sample a leaf within the selected tree using prefix sums
        const auto &bit = sel_tree->fenwick_by_dim_v[(size_t)k];
        double tot_leaf = src_info[si].total;
        if (bit.empty() || tot_leaf <= 0.0)
          continue;
        double r_leaf = rng_runif(0.0, tot_leaf);
        size_t leaf_idx_sel = rpf_utils::fenwick_find_by_prefix(bit, r_leaf);
        if (leaf_idx_sel == 0)
          continue;        // safety
        leaf_idx_sel -= 1; // to 0-based
        if (leaf_idx_sel >= sel_tree->leaves.size())
          continue;
        Leaf *leaf_ptr = &sel_tree->leaves[leaf_idx_sel];
        // Sample by unique thresholds: build/reuse unique values for this leaf and dim
        const std::vector<double> &uniques = rpf_ct1::unique_values(fr, X, *leaf_ptr, k);
        const size_t unique_count = uniques.size();
        const int left = leaf_size;
        const int right_exclusive = (int)unique_count - leaf_size + 1;
        if (right_exclusive <= left)
          continue;
        int s_idx = rng_randint(left, right_exclusive);
        Try tr;
        tr.tree = sel_tree;
        tr.leaf = leaf_ptr;
        tr.k = k;
        tr.idx = idx;
        tr.sp = uniques[(size_t)s_idx];
        tr.r_sp = rpf_ct1::rank_of(fr, k, tr.sp);
        tries.push_back(std::move(tr));
      }
    } else {
      const size_t from = tries.size();
      for (const auto &src_tree : sources) {
        if (src_tree->leaves.empty())
          continue;
        for (auto &leaf : src_tree->leaves) {
          // Build/reuse unique values for deterministic sampling across unique thresholds
          const std::vector<double> &uniques = rpf_ct1::unique_values(fr, X, leaf, k);
          const size_t unique_count = uniques.size();
          if ((int)unique_count < 2 * leaf_size)
            continue;
          int left = leaf_size;
          int right = (int)unique_count - leaf_size + 1;
          std::vector<int> samples = compute_even_spread_indices(left, right, (size_t)this->split_try);
          for (int s_idx : samples) {
            Try tr;
            tr.tree = src_tree;
            tr.leaf = &leaf;
            tr.k = k;
            tr.idx = idx;
            tr.sp = uniques[(size_t)s_idx];
            tr.r_sp = rpf_ct1::rank_of(fr, k, tr.sp);
            tries.push_back(std::move(tr));
          }
        }
      }
      // Deterministic: the thresholds of all of a candidate's leaves at once
      evaluate(from);
      compare(from);
      tries.clear();
    }
  }
  // Random: every candidate's tries at once (all draws are done)
  if (!this->deterministic) {
    evaluate(0);
    compare(0);
  }

  rpf_utils::age_pool_by_sample(sample_idxs, best_idx, possible_splits);
  finalize_split_from_sums(min_split, X, this->value_size);
  return min_split;
}
