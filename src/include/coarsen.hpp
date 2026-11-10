// The panel form of a purified component (prune_components(max_panels = P)): component f_S of order |S| >= 3 shown in
// at most P heat maps. |S| - 2 coarse features C are cut into intervals (cuts t_k, an interval I = [t_{i-1},
// t_i), the outer ones open); the two others (D) stay at full resolution. Within a cell the component is averaged with
// product-marginal weights: each rectangle's indicator on a coarse axis k is replaced by mu_k(B & I) / mu_k(I) for the
// cell's interval I (mu_k: the training rows' marginal of feature k). This is linear in the rectangles, so it applies
// to the merged store (the forest's component) and to one family's rectangles alike. Averaging over a coarse feature
// under mu_k commutes with it; the means along the two full-resolution features and over the training rows change, and
// as purify() identifies the components under the training rows' joint distribution (not the product of marginals),
// the coarsened component is not purified (it is not re-centred either).
//
// Coordinates: forest cells as in rects.hpp (feature d: cell c(x) = #{cuts[d] <= x}, 0..K_d; a box [lo, hi) of cells,
// lo = 0 / hi = K_d + 1 open). A rectangle of the group (S, j) restricts only the features of j; here every rectangle
// is written as an |S|-box with the whole range [0, K_d + 1) on the features outside j (indicator 1, and on a coarse
// axis mu_k(I) / mu_k(I) = 1, so nothing changes).
#ifndef RPF_COARSEN_HPP
#define RPF_COARSEN_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "rects.hpp"

namespace rpf_coarse {

// All rectangles of one component as |S|-boxes (positions in S), with weights (value size 1)
struct Boxes {
  int k = 0;
  std::vector<int32_t> lo, hi; // R x k, row-major
  std::vector<double> w;
  size_t size() const { return w.size(); }
  void push(const int32_t *l, const int32_t *h, double v)
  {
    lo.insert(lo.end(), l, l + k);
    hi.insert(hi.end(), h, h + k);
    w.push_back(v);
  }
};

// One family's (or the merged) accumulator entry of S -> Boxes; ncell[q] = K_d + 1 of the feature at position q
inline Boxes boxes_from_groups(
    int k, double constant, const std::map<int, rpf_rect::RectList> &groups, const std::vector<int32_t> &ncell
)
{
  Boxes B;
  B.k = k;
  std::vector<int32_t> l(k), h(k);
  if (constant != 0.0) {
    for (int q = 0; q < k; ++q) {
      l[q] = 0;
      h[q] = ncell[q];
    }
    B.push(l.data(), h.data(), constant);
  }
  for (const auto &g : groups) {
    std::vector<int> pos;
    for (int b = 0; b < k; ++b)
      if (g.first & (1 << b))
        pos.push_back(b);
    const rpf_rect::RectList &rl = g.second;
    for (size_t r = 0; r < rl.size(); ++r) {
      if (rl.w[r] == 0.0)
        continue;
      for (int q = 0; q < k; ++q) {
        l[q] = 0;
        h[q] = ncell[q];
      }
      for (size_t i = 0; i < pos.size(); ++i) {
        l[pos[i]] = rl.lo(r)[i];
        h[pos[i]] = rl.hi(r)[i];
      }
      B.push(l.data(), h.data(), rl.w[r]);
    }
  }
  return B;
}

// The merged store's component (its groups built: tables and k-d indices) -> Boxes, the same function: a table is
// written as runs of equal value, an index as its node constants (on the node's region) and its leaves' rectangles
// (each clipped to the leaf's region, where add_at() tests it)
inline Boxes boxes_from_component(const rpf_rect::Component &c, const std::vector<int32_t> &ncell)
{
  const int k = (int)c.S.size();
  Boxes B;
  B.k = k;
  std::vector<int32_t> l(k), h(k);
  auto full = [&]() {
    for (int q = 0; q < k; ++q) {
      l[q] = 0;
      h[q] = ncell[q];
    }
  };
  if (!c.constant.empty() && c.constant[0] != 0.0) {
    full();
    B.push(l.data(), h.data(), c.constant[0]);
  }
  for (const auto &g : c.groups) {
    const int m = (int)g.pos.size();
    if (!g.table.empty()) {
      const int p0 = g.pos[0];
      const size_t nc = g.table.size();
      size_t a = 0;
      while (a < nc) {
        size_t b = a + 1;
        while (b < nc && g.table[b] == g.table[a])
          ++b;
        if (g.table[a] != 0.0) {
          full();
          l[p0] = (int32_t)a;
          h[p0] = (int32_t)b;
          B.push(l.data(), h.data(), g.table[a]);
        }
        a = b;
      }
      continue;
    }
    const rpf_rect::KdIndex &ix = g.index;
    if (ix.nodes.empty())
      continue;
    struct Item {
      int node;
      std::vector<int32_t> rlo, rhi;
    };
    std::vector<Item> stack;
    {
      Item it;
      it.node = 0;
      for (int q = 0; q < m; ++q) {
        it.rlo.push_back(0);
        it.rhi.push_back(ncell[g.pos[q]]);
      }
      stack.push_back(std::move(it));
    }
    while (!stack.empty()) {
      Item it = std::move(stack.back());
      stack.pop_back();
      const rpf_rect::KdIndex::Node &nd = ix.nodes[it.node];
      if (nd.cst >= 0 && ix.csts[nd.cst] != 0.0) {
        full();
        for (int q = 0; q < m; ++q) {
          l[g.pos[q]] = it.rlo[q];
          h[g.pos[q]] = it.rhi[q];
        }
        B.push(l.data(), h.data(), ix.csts[nd.cst]);
      }
      if (nd.left < 0) {
        for (int i = nd.r0; i < nd.r1; ++i) {
          const size_t r = (size_t)ix.ids[i];
          if (g.rects.w[r] == 0.0)
            continue;
          full();
          bool empty = false;
          for (int q = 0; q < m; ++q) {
            int32_t a = std::max(g.rects.lo(r)[q], it.rlo[q]), b = std::min(g.rects.hi(r)[q], it.rhi[q]);
            if (a >= b)
              empty = true;
            l[g.pos[q]] = a;
            h[g.pos[q]] = b;
          }
          if (!empty)
            B.push(l.data(), h.data(), g.rects.w[r]);
        }
        continue;
      }
      Item L{nd.left, it.rlo, it.rhi}, Rr{nd.right, it.rlo, it.rhi};
      L.rhi[nd.dim] = std::min(L.rhi[nd.dim], (int32_t)nd.split);
      Rr.rlo[nd.dim] = std::max(Rr.rlo[nd.dim], (int32_t)nd.split);
      stack.push_back(std::move(Rr));
      stack.push_back(std::move(L));
    }
  }
  return B;
}

// A coarse axis: feature `feat` (0-based) at position `pos` of S, cut at t; PJ[c * q + i] = the training rows with
// forest cell < c in interval i (c = 0..K + 1), cntI[i] = the rows in interval i
struct Axis {
  int pos = 0, feat = 0, q = 1;
  std::vector<double> t, PJ, cntI;
  // mu(B & I_i) / mu(I_i) for the box [lo, hi) of forest cells
  double a(int32_t lo, int32_t hi, int i) const { return (PJ[(size_t)hi * q + i] - PJ[(size_t)lo * q + i]) / cntI[i]; }
  int interval(double x) const { return (int)(std::upper_bound(t.begin(), t.end(), x) - t.begin()); }
};

// cells: the training rows' forest cells of the feature (n), x: their values; ncell = K + 1
inline Axis make_axis(
    int pos,
    int feat,
    const std::vector<double> &t,
    const std::vector<int32_t> &cells,
    const std::vector<double> &x,
    int32_t ncell
)
{
  Axis A;
  A.pos = pos;
  A.feat = feat;
  A.t = t;
  A.q = (int)t.size() + 1;
  const int q = A.q;
  std::vector<double> cnt((size_t)ncell * q, 0.0);
  A.cntI.assign(q, 0.0);
  for (size_t i = 0; i < x.size(); ++i) {
    int iv = A.interval(x[i]);
    cnt[(size_t)cells[i] * q + iv] += 1.0;
    A.cntI[iv] += 1.0;
  }
  A.PJ.assign((size_t)(ncell + 1) * q, 0.0);
  for (int32_t c = 0; c < ncell; ++c)
    for (int i = 0; i < q; ++i)
      A.PJ[(size_t)(c + 1) * q + i] = A.PJ[(size_t)c * q + i] + cnt[(size_t)c * q + i];
  return A;
}

// The coarsened component: rectangles over the full-resolution features D (positions dpos in S) with one weight per
// coarse cell (row-major over the coarse axes in the order of `axes`), sorted and merged; `index` built on a copy
inline const Axis &ax(const Axis &a)
{
  return a;
}
inline const Axis &ax(const Axis *a)
{
  return *a;
}

template <class V> inline int n_cells_of(const V &axes)
{
  int nc = 1;
  for (const auto &a : axes)
    nc *= ax(a).q;
  return nc;
}

// The coarse weights of rectangle r of B: out[0..ncells) (axes: Axis or const Axis *)
template <class V> inline void coarse_weights(const Boxes &B, size_t r, const V &axes, double *out, int nc)
{
  out[0] = B.w[r];
  int len = 1;
  for (const auto &a : axes) {
    const Axis &A = ax(a);
    const int32_t lo = B.lo[r * B.k + A.pos], hi = B.hi[r * B.k + A.pos];
    // new index = old * q + i (row-major, earlier axes more significant)
    for (int c = len - 1; c >= 0; --c) {
      const double v = out[c];
      for (int i = A.q - 1; i >= 0; --i)
        out[c * A.q + i] = v * A.a(lo, hi, i);
    }
    len *= A.q;
  }
  (void)nc;
}

// B coarsened on `axes`, as a RectList over the positions dpos (m = |dpos|, vs = ncells), sorted and merged, exact
// zeros dropped
template <class V> inline rpf_rect::RectList coarse_rects(const Boxes &B, const V &axes, const std::vector<int> &dpos)
{
  const int nc = n_cells_of(axes), m = (int)dpos.size();
  rpf_rect::RectList rl(m, nc);
  std::vector<double> w(nc);
  std::vector<int32_t> l(std::max(1, m)), h(std::max(1, m));
  if (m >= 1 && m <= 2) {
    // Equal boxes merged as they come (their weights added in the order of B, as sort_merge() adds them), then
    // written in sorted order: the same list as below, bitwise, without one rectangle per box of B held at once
    struct KeyHash {
      size_t operator()(const std::pair<uint64_t, uint64_t> &k) const
      {
        uint64_t x = k.first * 0x9e3779b97f4a7c15ULL ^ (k.second + 0x7f4a7c159e3779b9ULL + (k.first << 6));
        x = (x ^ (x >> 31)) * 0xbf58476d1ce4e5b9ULL;
        return (size_t)(x ^ (x >> 29));
      }
    };
    auto u = [](int32_t v) { return (uint64_t)(uint32_t)v; };
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t, KeyHash> at;
    std::vector<std::array<int32_t, 4>> keys; // lo[0..m), hi[0..m), zero-padded
    std::vector<double> acc;
    for (size_t r = 0; r < B.size(); ++r) {
      coarse_weights(B, r, axes, w.data(), nc);
      std::array<int32_t, 4> key{{0, 0, 0, 0}};
      for (int q = 0; q < m; ++q) {
        key[q] = B.lo[r * B.k + dpos[q]];
        key[m + q] = B.hi[r * B.k + dpos[q]];
      }
      const std::pair<uint64_t, uint64_t> kk{(u(key[0]) << 32) | u(key[1]), (u(key[2]) << 32) | u(key[3])};
      auto it = at.find(kk);
      if (it == at.end()) {
        at.emplace(kk, keys.size());
        keys.push_back(key);
        acc.insert(acc.end(), w.begin(), w.end());
      } else {
        double *dst = acc.data() + it->second * nc;
        for (int c = 0; c < nc; ++c)
          dst[c] += w[c];
      }
    }
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t, KeyHash>().swap(at);
    std::vector<size_t> perm(keys.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b) {
      return std::lexicographical_compare(
          keys[a].begin(), keys[a].begin() + 2 * m, keys[b].begin(), keys[b].begin() + 2 * m
      );
    });
    rl.box.reserve(keys.size() * 2 * m);
    rl.w.reserve(keys.size() * nc);
    for (size_t i : perm)
      rl.push(keys[i].data(), keys[i].data() + m, acc.data() + i * nc, 1.0);
    rl.drop_zeros();
    return rl;
  }
  for (size_t r = 0; r < B.size(); ++r) {
    coarse_weights(B, r, axes, w.data(), nc);
    for (int q = 0; q < m; ++q) {
      l[q] = B.lo[r * B.k + dpos[q]];
      h[q] = B.hi[r * B.k + dpos[q]];
    }
    rl.push(l.data(), h.data(), w.data(), 1.0);
  }
  rl.sort_merge();
  rl.drop_zeros();
  return rl;
}

// A built coarsened component for point queries
struct Built {
  int m = 0, nc = 1;
  rpf_rect::RectList rl; // after the index build (only the rectangles its leaves test)
  rpf_rect::KdIndex index;
  std::vector<double> cst; // m = 0: the sum of all weights per cell
  void build(const rpf_rect::RectList &stored, const std::vector<int32_t> &ncell_d)
  {
    m = stored.m;
    nc = stored.vs;
    if (m == 0) {
      cst.assign(nc, 0.0);
      for (size_t r = 0; r < stored.size(); ++r)
        for (int c = 0; c < nc; ++c)
          cst[c] += stored.w[r * nc + c];
      return;
    }
    rl = stored;
    index.build(rl, ncell_d);
  }
  // the value at D-cells dc[0..m) and coarse cell cc; buf: nc doubles of scratch
  double at(const int32_t *dc, int cc, double *buf) const
  {
    if (m == 0)
      return cst[cc];
    std::fill(buf, buf + nc, 0.0);
    index.add_at(rl, dc, buf);
    return buf[cc];
  }
};

// The maximal q vectors (q_i >= 1, product <= P, no q_i can grow) of n axes, in lexicographic order
inline void q_vectors(int n, long long P, std::vector<std::vector<int>> &out)
{
  out.clear();
  if (n == 0) {
    out.push_back({});
    return;
  }
  std::vector<int> q(n, 1);
  std::function<void(int, long long)> rec = [&](int i, long long prod) {
    if (i == n) {
      for (int a = 0; a < n; ++a)
        if (prod / q[a] * (q[a] + 1) <= P)
          return;
      out.push_back(q);
      return;
    }
    for (int v = 1; prod * v <= P; ++v) {
      q[i] = v;
      rec(i + 1, prod * v);
    }
    q[i] = 1;
  };
  rec(0, 1);
}

// Equal-mass cuts of the sorted values xs into q intervals (ties merged; only cuts above the smallest value)
inline std::vector<double> eq_cuts(const std::vector<double> &xs, int q)
{
  std::vector<double> t;
  const long long n = (long long)xs.size();
  for (long long i = 1; i < q; ++i) {
    long long idx = (long long)std::ceil((double)(i * n) / (double)q);
    if (idx > n - 1)
      idx = n - 1;
    t.push_back(xs[(size_t)idx]);
  }
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());
  std::vector<double> out;
  for (double v : t)
    if (v > xs[0])
      out.push_back(v);
  return out;
}

inline double var_pop(const std::vector<double> &v)
{
  if (v.empty())
    return 0.0;
  double s = 0.0;
  for (double x : v)
    s += x;
  const double mu = s / (double)v.size();
  double ss = 0.0;
  for (double x : v)
    ss += (x - mu) * (x - mu);
  return ss / (double)v.size();
}

// ---- The results of RandomPlantedForest::coarsen_select() / coarsen_oob() and the predictor of a stored form (std
// types; the R binding converts them)

// The stored form of a coarsened component: what predict needs, no training data. The rectangles are over the two
// full-resolution features `dfeat` (in forest cells of `dcuts`), with one weight per coarse cell (row-major over the
// coarse features `cfeat`, cut at `ccuts`). Features are 1-based. A value is the sum over the families: divide by
// n_trees.
struct Form {
  std::vector<int> S, cfeat, dfeat;
  std::vector<std::vector<double>> ccuts, dcuts;
  rpf_rect::RectList rects;
  int n_trees = 1;
};

// One candidate (C, q) of a component: the coarse set (index into the (|S| - 2)-subsets of S in lexicographic order,
// 0-based), the q vector, its cells and the share of the component's training-row variance it keeps
struct Candidate {
  int combo = 0;
  std::vector<int> q;
  int cells = 1;
  double var_kept = 0;
};

// One component's panel form
struct Selected {
  std::vector<int> S; // 1-based, ascending
  int form = 0;       // 0 absent (zero), 1 as it is (its own cut points need at most P panels), 2 coarsened
  double forest_panels = 0, var_forest = 0, var_kept = 1, seconds = 0;
  // coarsened: positions in S of the coarse and the full features, q, cuts per coarse feature; as it is: cpos and cuts
  // are its panels (the |S| - 2 features with the fewest own cut points, and these cut points)
  std::vector<int> cpos, dpos, q;
  std::vector<std::vector<double>> cuts;
  std::vector<Candidate> cands;
  Form stored;
  std::vector<double> g; // the coarsened component at the training rows (divided by n_trees)
  int n_diff = 0, n_kd = 0;
};

// A coarsening to apply: component S coarsened on C (1-based features of S, in this order) with these cuts
struct Spec {
  std::vector<int> S, C;
  std::vector<std::vector<double>> cuts;
};

// Each family's coarsened components at its out-of-bag rows, summed over the families of even (sum0) and odd (sum1)
// index, in family order: column-major n x specs.size()
struct OobSums {
  std::vector<double> sum0, sum1;
};

// A stored form built for point queries (k-d index over its rectangles)
class Model {
public:
  // Throws std::invalid_argument if the form is inconsistent
  explicit Model(const Form &f);
  // The values at the n rows of column-major X (n x p, all features)
  std::vector<double> eval(const double *X, int n, int p) const;

private:
  Built bt;
  std::vector<int> dfeat, cfeat; // 0-based
  std::vector<std::vector<double>> dcuts, ccuts;
  double n_trees = 1;
};

} // namespace rpf_coarse

#endif
