// The pruning without the purified store of all components (prune_components()). purify(mode = 2) splits every leaf
// into one piece per subset of its tree's features and keeps all of them, so at high orders the store does not fit in
// memory (California housing at order 8: more than 24 GB), while the pruning needs only each component's values at the
// training rows and its out-of-bag sums, and the pruned model only the pieces of its kept components.
// - purify_columns(): one pass over the families; each family is purified and its components are evaluated at the
//   training rows (family_columns()), then freed.
// - store_components(): the pieces of chosen components, merged as purify() merges and stores them; BuiltStore
//   evaluates them as predict_components() evaluates the purified forest.
#include "rpf.hpp"
#include "subset_count.hpp"
#include <array>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace {
int popcount32(unsigned int x)
{
  x = x - ((x >> 1) & 0x55555555u);
  x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
  return (int)((((x + (x >> 4)) & 0x0F0F0F0Fu) * 0x01010101u) >> 24);
}

// A set of 1-based features, sorted; throws unless distinct and in 1..p
std::vector<int> checked_set(std::vector<int> S, int p)
{
  std::sort(S.begin(), S.end());
  if (S.empty() || S.front() < 1 || S.back() > p || std::adjacent_find(S.begin(), S.end()) != S.end())
    throw std::invalid_argument("A component must be a set of distinct features in 1..p.");
  return S;
}
} // namespace

// One family's purified components at n rows, the same purification as purify_fast_exact_family(): a leaf l of tree T
// writes, for every S subset of T and j subset of S, the term (-1)^(|S| - |j|) c_j(l) 1{x_j in box_j(l)} into f_S,
// with c_j(l) = cnt(box on T \ j) / tot(T \ j) * value(l) (the same range counts and totals). The term does not depend
// on S, so instead of one rectangle per (S, j) (3^|T| per leaf) the tree's g_{T,j}(x) = sum_l c_j(l) 1{x_j in
// box_j(l)} are evaluated at the rows, one per j (2^|T|), and f_S^T(x) = sum_{j subset S} (-1)^(|S| - |j|) g_{T,j}(x)
// follows by a subset Moebius transform per row. g_{T,j} is a step function on the tree's own grid over j (its leaves'
// boundaries): a dense table built by a difference array if the grid is small, otherwise a k-d index over the tree's
// leaves. The values equal the rectangle store's up to rounding.
void RandomPlantedForest::family_columns(
    const TreeFamily &curr_family,
    int maxp_interaction,
    const std::vector<std::vector<double>> &cuts,
    const int32_t *cells,
    size_t stride,
    size_t n,
    const rpf_cnt::SubsetCounter &counter,
    std::map<std::vector<int>, std::vector<double>> &out,
    std::vector<double> &intercept
) const
{
  if (maxp_interaction <= 0 || maxp_interaction > feature_size)
    maxp_interaction = feature_size;
  const int VS = (int)value_size;
  const size_t dense_max = std::max<size_t>(4096, 4 * n);

  // as purify_fast_exact_family(): the family's cut end points, the totals and the cached range counts
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
  intercept.assign(VS, 0.0);
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
    for (int u : U) {
      const auto &lims = lim_list[u - 1];
      consU.push_back({u - 1, lims.front(), lims.back()});
    }
    size_t cnt = consU.empty() ? (size_t)sample_size : counter.range_count(consU);
    double tot = (double)cnt;
    tot_cache.insert({U, tot});
    return tot;
  };
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

  std::vector<double> G; // per j (2^|T|): n x VS
  for (const auto &kv : curr_family) {
    const std::set<int> &T = kv.first;
    if (T == std::set<int>{0})
      continue;
    const auto &leaves = kv.second->leaves;
    if (leaves.empty())
      continue;
    std::vector<int> Tvec(T.begin(), T.end());
    const int tdim = (int)Tvec.size();
    const size_t nm = (size_t)1 << tdim, L = leaves.size();

    // per leaf: its box in forest cells and its weights c_j (as purify_fast_exact_family())
    std::vector<int32_t> LO(L * tdim), HI(L * tdim);
    std::vector<double> C(L * nm * VS, 0.0);
    std::vector<rpf_kd::RangeConstraint> rc_by_pos(tdim), cons;
    std::vector<char> rc_ok(tdim, 0);
    std::vector<int> rc_lo_idx(tdim, -1), rc_hi_idx(tdim, -1), key;
    // tot(T \ j) per j, once per tree
    std::vector<double> totm(nm, 0.0);
    for (size_t jmask = 0; jmask < nm; ++jmask) {
      if (popcount32((unsigned)jmask) > maxp_interaction)
        continue;
      std::set<int> U;
      for (int b = 0; b < tdim; ++b)
        if (!(jmask & ((size_t)1 << b)))
          U.insert(Tvec[b]);
      totm[jmask] = get_tot_for_U(U);
    }
    for (size_t l = 0; l < L; ++l) {
      const auto &leaf = leaves[l];
      for (int b = 0; b < tdim; ++b) {
        int d = Tvec[b];
        const auto &lims = lim_list[d - 1];
        const auto &cd = cuts[d - 1];
        double a = leaf.intervals[d - 1].first, e = leaf.intervals[d - 1].second;
        LO[l * tdim + b] =
            (a <= lims.front()) ? 0 : (int32_t)(std::lower_bound(cd.begin(), cd.end(), a) - cd.begin()) + 1;
        HI[l * tdim + b] = (e >= lims.back()) ? (int32_t)cd.size() + 1
                                              : (int32_t)(std::lower_bound(cd.begin(), cd.end(), e) - cd.begin()) + 1;
        rc_ok[b] = 0;
        if ((int)lims.size() < 2)
          continue;
        double lo = std::max(a, lims.front());
        double hi = std::min(e, lims.back());
        if (hi > lo) {
          rc_by_pos[b] = {d - 1, lo, hi};
          rc_ok[b] = 1;
          rc_lo_idx[b] = (int)(std::lower_bound(lims.begin(), lims.end(), lo) - lims.begin());
          rc_hi_idx[b] = (int)(std::lower_bound(lims.begin(), lims.end(), hi) - lims.begin());
        }
      }
      for (size_t jmask = 0; jmask < nm; ++jmask) {
        if (popcount32((unsigned)jmask) > maxp_interaction)
          continue;
        const double totU = totm[jmask];
        if (totU <= 0.0)
          continue;
        cons.clear();
        key.clear();
        bool empty_range = false;
        for (int b = 0; b < tdim; ++b) {
          if (jmask & ((size_t)1 << b))
            continue;
          if (!rc_ok[b]) {
            empty_range = true;
            break;
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
            cnt = counter.range_count(cons);
            kd_cache.emplace(key, cnt);
          }
        }
        for (int v = 0; v < VS; ++v)
          C[(l * nm + jmask) * VS + v] = ((double)cnt / totU) * leaf.value[v];
      }
      for (int v = 0; v < VS; ++v)
        intercept[v] += C[(l * nm) * VS + v];
    }

    // the tree's grid: per position b the distinct box boundaries; rows and boxes in its coordinates (cell c lies in
    // [lo, hi) iff pos(lo) <= pos(c) < pos(hi), pos(c) = #boundaries <= c)
    std::vector<std::vector<int32_t>> Bd(tdim), P(tdim), PLO(tdim), PHI(tdim);
    for (int b = 0; b < tdim; ++b) {
      auto &B = Bd[b];
      for (size_t l = 0; l < L; ++l) {
        B.push_back(LO[l * tdim + b]);
        B.push_back(HI[l * tdim + b]);
      }
      std::sort(B.begin(), B.end());
      B.erase(std::unique(B.begin(), B.end()), B.end());
      auto pos = [&B](int32_t c) { return (int32_t)(std::upper_bound(B.begin(), B.end(), c) - B.begin()); };
      P[b].resize(n);
      const int d0 = Tvec[b] - 1;
      for (size_t i = 0; i < n; ++i)
        P[b][i] = pos(cells[i * stride + d0]);
      PLO[b].resize(L);
      PHI[b].resize(L);
      for (size_t l = 0; l < L; ++l) {
        PLO[b][l] = pos(LO[l * tdim + b]);
        PHI[b][l] = pos(HI[l * tdim + b]);
      }
    }

    G.assign(nm * n * VS, 0.0);
    // j = {}: the constant sum_l c_{}(l)
    {
      std::vector<double> c0(VS, 0.0);
      for (size_t l = 0; l < L; ++l)
        for (int v = 0; v < VS; ++v)
          c0[v] += C[(l * nm) * VS + v];
      for (size_t i = 0; i < n; ++i)
        for (int v = 0; v < VS; ++v)
          G[i * VS + v] = c0[v];
    }
    std::vector<int> jb;
    std::vector<size_t> sz, st;
    std::vector<double> tab;
    std::vector<int32_t> c(tdim), lo(tdim), hi(tdim);
    for (size_t jmask = 1; jmask < nm; ++jmask) {
      const int m = popcount32((unsigned)jmask);
      if (m > maxp_interaction)
        continue;
      jb.clear();
      for (int b = 0; b < tdim; ++b)
        if (jmask & ((size_t)1 << b))
          jb.push_back(b);
      double cells_prod = 1;
      sz.assign(m, 0);
      st.assign(m, 0);
      for (int q = m - 1; q >= 0; --q) {
        sz[q] = Bd[jb[q]].size() + 1;
        st[q] = (q == m - 1) ? 1 : st[q + 1] * sz[q + 1];
        cells_prod *= (double)sz[q];
      }
      double *Gj = G.data() + jmask * n * VS;
      if (cells_prod <= (double)dense_max) {
        // difference array over the grid: each box adds its weight at its 2^m corners with signs
        const size_t ncell = (size_t)cells_prod;
        tab.assign(ncell * VS, 0.0);
        for (size_t l = 0; l < L; ++l) {
          const double *w = C.data() + (l * nm + jmask) * VS;
          bool nz = false;
          for (int v = 0; v < VS; ++v)
            nz = nz || w[v] != 0.0;
          if (!nz)
            continue;
          bool empty = false;
          for (int q = 0; q < m; ++q)
            empty = empty || PLO[jb[q]][l] >= PHI[jb[q]][l];
          if (empty)
            continue;
          for (size_t corner = 0; corner < ((size_t)1 << m); ++corner) {
            size_t idx = 0;
            bool out_of = false;
            int sgn = 1;
            for (int q = 0; q < m; ++q) {
              const bool up = (corner >> q) & 1;
              const size_t x = up ? (size_t)PHI[jb[q]][l] : (size_t)PLO[jb[q]][l];
              if (x >= sz[q])
                out_of = true; // the upper end of the grid: no correction needed
              idx += x * st[q];
              if (up)
                sgn = -sgn;
            }
            if (out_of)
              continue;
            for (int v = 0; v < VS; ++v)
              tab[idx * VS + v] += sgn * w[v];
          }
        }
        for (int q = 0; q < m; ++q)
          for (size_t idx = 0; idx < ncell; ++idx)
            if ((idx / st[q]) % sz[q] != 0)
              for (int v = 0; v < VS; ++v)
                tab[idx * VS + v] += tab[(idx - st[q]) * VS + v];
        for (size_t i = 0; i < n; ++i) {
          size_t idx = 0;
          for (int q = 0; q < m; ++q)
            idx += (size_t)P[jb[q]][i] * st[q];
          for (int v = 0; v < VS; ++v)
            Gj[i * VS + v] = tab[idx * VS + v];
        }
      } else {
        rpf_rect::RectList rl(m, VS);
        for (size_t l = 0; l < L; ++l) {
          const double *w = C.data() + (l * nm + jmask) * VS;
          bool nz = false;
          for (int v = 0; v < VS; ++v)
            nz = nz || w[v] != 0.0;
          bool empty = false;
          for (int q = 0; q < m; ++q) {
            lo[q] = PLO[jb[q]][l];
            hi[q] = PHI[jb[q]][l];
            empty = empty || lo[q] >= hi[q];
          }
          if (nz && !empty)
            rl.push(lo.data(), hi.data(), w, 1.0);
        }
        std::vector<int32_t> nc(m);
        for (int q = 0; q < m; ++q)
          nc[q] = (int32_t)sz[q];
        rpf_rect::KdIndex index;
        index.build(rl, nc);
        for (size_t i = 0; i < n; ++i) {
          for (int q = 0; q < m; ++q)
            c[q] = P[jb[q]][i];
          index.add_at(rl, c.data(), Gj + i * VS);
        }
      }
    }
    // f_S^T = sum_{j subset S} (-1)^(|S| - |j|) g_{T,j}: the subset Moebius transform, for |S| <= maxp
    for (int b = 0; b < tdim; ++b)
      for (size_t mask = 0; mask < nm; ++mask) {
        if (!(mask & ((size_t)1 << b)) || popcount32((unsigned)mask) > maxp_interaction)
          continue;
        double *dst = G.data() + mask * n * VS;
        const double *src = G.data() + (mask ^ ((size_t)1 << b)) * n * VS;
        for (size_t e = 0; e < n * VS; ++e)
          dst[e] -= src[e];
      }
    std::vector<int> S;
    for (size_t mask = 1; mask < nm; ++mask) {
      if (popcount32((unsigned)mask) > maxp_interaction)
        continue;
      S.clear();
      for (int b = 0; b < tdim; ++b)
        if (mask & ((size_t)1 << b))
          S.push_back(Tvec[b]);
      std::vector<double> &o = out[S];
      if (o.empty())
        o.assign(n * VS, 0.0);
      const double *src = G.data() + mask * n * VS;
      for (size_t e = 0; e < n * VS; ++e)
        o[e] += src[e];
    }
  }
}

ComponentColumns RandomPlantedForest::purify_columns(int maxp_interaction, int nthreads_param, bool oob)
{
  if (sample_size == 0 || X.empty())
    throw std::runtime_error("Cannot purify: no training data available.");
  if (oob) {
    bool ok = inbag_.size() == tree_families.size() && !tree_families.empty();
    for (const auto &inb : inbag_)
      ok = ok && inb.size() == (size_t)sample_size;
    if (!ok)
      throw std::invalid_argument("The forest has no in-bag counts.");
  }
  const unsigned threads_to_use = purify_threads(nthreads_param);
  rpf_rect::PurifiedForest pf; // only its cells: those of purify()'s store
  pf.cuts = forest_cuts();
  const size_t n = (size_t)sample_size, B = tree_families.size(), vs = value_size;
  const size_t stride = (size_t)std::max(1, feature_size);
  std::vector<int32_t> cells(n * stride, 0);
  for (size_t i = 0; i < n; ++i)
    for (int d = 0; d < feature_size; ++d)
      cells[i * stride + d] = pf.cell(d, X[i][d]);
  rpf_cnt::SubsetCounter counter(&X, feature_size);

  // per component: the sum over all families at every training row, and with oob over the OOB families of each half
  // (n x vs, row-major)
  struct Col {
    std::vector<double> own;
    std::array<std::vector<double>, 2> half;
  };
  std::map<std::vector<int>, Col> cols;
  std::array<std::vector<double>, 2> icpt;
  std::vector<int> count;
  if (oob) {
    for (int h = 0; h < 2; ++h)
      icpt[h].assign(n * vs, 0.0);
    count.assign(n * 2, 0);
  }
  // one family's components at all training rows, between their evaluation and their addition
  struct Slot {
    std::map<std::vector<int>, std::vector<double>> vals;
    std::vector<double> intercept;
  };
  std::vector<Slot> slots(threads_to_use);
  // the families' intercepts summed through purify_mode2()'s binary counter, so that the intercept is the store's
  std::vector<std::pair<int, std::vector<double>>> istack;

  run_families(
      B,
      threads_to_use,
      [&](size_t f) {
        Slot &sl = slots[f % threads_to_use];
        sl = Slot();
        this->family_columns(
            this->tree_families[f], maxp_interaction, pf.cuts, cells.data(), stride, n, counter, sl.vals, sl.intercept
        );
      },
      [&](size_t start, size_t batch) {
        // in family order, on the calling thread: the sums do not depend on the number of threads
        for (size_t i = 0; i < batch; ++i) {
          const size_t f = start + i;
          const int h = (int)(f % 2);
          Slot &sl = slots[i];
          for (auto it = sl.vals.begin(); it != sl.vals.end(); it = sl.vals.erase(it)) {
            Col &c = cols[it->first];
            if (c.own.empty()) {
              c.own.assign(n * vs, 0.0);
              if (oob) {
                c.half[0].assign(n * vs, 0.0);
                c.half[1].assign(n * vs, 0.0);
              }
            }
            const double *v = it->second.data();
            double *own = c.own.data();
            for (size_t e = 0; e < n * vs; ++e)
              own[e] += v[e];
            if (oob) {
              const std::vector<uint8_t> &inb = this->inbag_[f];
              double *hf = c.half[h].data();
              for (size_t r = 0; r < n; ++r)
                if (inb[r] == 0)
                  for (size_t q = 0; q < vs; ++q)
                    hf[r * vs + q] += v[r * vs + q];
            }
          }
          if (oob) {
            const std::vector<uint8_t> &inb = this->inbag_[f];
            for (size_t r = 0; r < n; ++r) {
              if (inb[r] != 0)
                continue;
              for (size_t q = 0; q < vs; ++q)
                icpt[h][r * vs + q] += sl.intercept[q];
              count[r * 2 + h] += 1;
            }
          }
          istack.emplace_back(0, std::move(sl.intercept));
          while (istack.size() >= 2 && istack[istack.size() - 2].first == istack.back().first) {
            auto &a = istack[istack.size() - 2].second;
            const auto &b = istack.back().second;
            for (size_t q = 0; q < a.size(); ++q)
              a[q] += b[q];
            istack[istack.size() - 2].first += 1;
            istack.pop_back();
          }
          sl = Slot();
        }
      }
  );
  std::vector<double> total;
  for (const auto &s : istack) {
    if (total.empty())
      total = s.second;
    else
      for (size_t q = 0; q < total.size(); ++q)
        total[q] += s.second[q];
  }
  if (total.empty())
    total.assign(vs, 0.0);

  // to column-major matrices; each component's sums are freed once copied
  ComponentColumns out;
  const size_t K = cols.size();
  out.own.assign(n * K * vs, 0.0);
  if (oob) {
    out.oob.sum0.assign(n * K * vs, 0.0);
    out.oob.sum1.assign(n * K * vs, 0.0);
  }
  size_t k = 0;
  for (auto it = cols.begin(); it != cols.end(); it = cols.erase(it), ++k) {
    out.oob.sets.push_back(it->first);
    for (size_t q = 0; q < vs; ++q)
      for (size_t i = 0; i < n; ++i) {
        const size_t e = i * vs + q, col = (k * vs + q) * n + i;
        out.own[col] = it->second.own[e] / n_trees;
        if (oob) {
          out.oob.sum0[col] = it->second.half[0][e];
          out.oob.sum1[col] = it->second.half[1][e];
        }
      }
  }
  if (oob) {
    out.oob.intercept0.assign(n * vs, 0.0);
    out.oob.intercept1.assign(n * vs, 0.0);
    out.oob.count.assign(n * 2, 0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t q = 0; q < vs; ++q) {
        out.oob.intercept0[q * n + i] = icpt[0][i * vs + q];
        out.oob.intercept1[q * n + i] = icpt[1][i * vs + q];
      }
      out.oob.count[i] = count[i * 2];
      out.oob.count[n + i] = count[i * 2 + 1];
    }
  }
  out.intercept.assign(vs, 0.0);
  for (size_t q = 0; q < vs; ++q)
    out.intercept[q] = total[q] / n_trees;
  return out;
}

ComponentStore RandomPlantedForest::store_components(const std::vector<std::vector<int>> &sets, int nthreads) const
{
  if (sample_size == 0 || X.empty())
    throw std::runtime_error("Cannot purify: no training data available.");
  std::set<std::vector<int>> only;
  for (const auto &s : sets)
    only.insert(checked_set(s, feature_size));
  ComponentStore st;
  st.p = feature_size;
  st.value_size = (int)value_size;
  st.n_trees = n_trees;
  st.cuts = forest_cuts();
  st.pieces = merged_components(only, std::max(1, nthreads), st.cuts);
  st.pieces.intercept.clear();
  // as PurifiedForest::build_from(): the exact zeros dropped, empty groups left out
  for (auto &c : st.pieces.comps) {
    if (c.second.first.empty())
      c.second.first.assign(value_size, 0.0);
    auto &groups = c.second.second;
    for (auto g = groups.begin(); g != groups.end();) {
      g->second.drop_zeros();
      if (g->second.size() == 0)
        g = groups.erase(g);
      else
        ++g;
    }
  }
  return st;
}

BuiltStore::BuiltStore(ComponentStore &store, unsigned build_threads)
{
  auto fail = []() { throw std::invalid_argument("Corrupt component store."); };
  const int p = store.p;
  vs = store.value_size;
  n_trees = store.n_trees;
  if (p < 1 || vs < 1 || n_trees < 1 || (int)store.cuts.size() != p)
    fail();
  for (const auto &c : store.cuts)
    for (size_t i = 1; i < c.size(); ++i)
      if (!(c[i - 1] < c[i]))
        fail();
  for (const auto &kv : store.pieces.comps) {
    const std::vector<int> &S = kv.first;
    if (S.empty() || S.size() > 30 || checked_set(S, p) != S || (int)kv.second.first.size() != vs)
      fail();
    for (const auto &g : kv.second.second) {
      const int jmask = g.first, m = popcount32((unsigned)jmask);
      if (jmask < 1 || jmask >= (1 << S.size()))
        fail();
      const rpf_rect::RectList &rl = g.second;
      const size_t R = rl.size();
      if (rl.m != m || rl.vs != vs || rl.box.size() != R * 2 * (size_t)m || rl.w.size() != R * (size_t)vs)
        fail();
      std::vector<int32_t> nc;
      for (size_t b = 0; b < S.size(); ++b)
        if ((jmask >> b) & 1)
          nc.push_back((int32_t)store.cuts[S[b] - 1].size() + 1);
      for (size_t r = 0; r < R; ++r)
        for (int q = 0; q < m; ++q)
          if (rl.lo(r)[q] < 0 || rl.lo(r)[q] >= rl.hi(r)[q] || rl.hi(r)[q] > nc[q])
            fail();
    }
  }
  pf.cuts = std::move(store.cuts);
  pf.maxp = p;
  pf.build_from(store.pieces, vs, std::max(1u, build_threads));
}

std::vector<double>
BuiltStore::eval(const double *X, int n, int p, const std::vector<std::vector<int>> &sets, int nthreads) const
{
  if (p != (int)pf.cuts.size())
    throw std::invalid_argument("Feature vector has wrong dimension.");
  std::vector<const rpf_rect::Component *> comps;
  for (const auto &s : sets) {
    auto it = pf.comps.find(checked_set(s, p));
    if (it == pf.comps.end())
      throw std::invalid_argument("A component is not in the store.");
    comps.push_back(&it->second);
  }
  std::vector<double> out((size_t)n * comps.size() * vs, 0.0);
  if (n > 0 && !comps.empty())
    rects_components(pf, comps, X, n, vs, n_trees, (unsigned)std::max(1, nthreads), out.data());
  return out;
}
