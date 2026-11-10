// Purified components stored as sums of weighted axis-parallel rectangles (purify mode 2).
//
// Mode 2 writes the purified component f_S of a family as a sum of terms
//   sign * w * 1{x_j in box_j(leaf)},   j subset of S,
// one per leaf of a tree T containing S and per j (the indicator only restricts the features in j). Instead of
// rasterising these terms onto a dense grid over all of a family's cut points, they are kept as rectangles, grouped
// by (S, j), merged across leaves and families (equal boxes add their weights) and indexed for point queries.
//
// Coordinates are integer cells. For feature d, `cuts[d]` holds the sorted distinct interval end points of all
// leaves; a value x lies in cell c(x) = #{k : cuts[d][k] <= x} in [0, K_d], and a box [a, b) covers the cells
// [index(a) + 1, index(b) + 1). An end point at the family's outermost cut is open (cell 0 or K_d + 1), which
// reproduces the clamping of the dense grids for points outside the training range.
#ifndef RPF_RECTS_HPP
#define RPF_RECTS_HPP

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace rpf_rect {
// A list of boxes over m features with weights of length vs. box = lo[0..m), hi[0..m) per rectangle.
struct RectList {
  int m = 0;
  int vs = 1;
  std::vector<int32_t> box;
  std::vector<double> w;

  RectList() {}
  RectList(int m_, int vs_) : m(m_), vs(vs_) {}
  size_t size() const { return vs > 0 ? w.size() / (size_t)vs : 0; }
  const int32_t *lo(size_t r) const { return box.data() + r * 2 * (size_t)m; }
  const int32_t *hi(size_t r) const { return box.data() + r * 2 * (size_t)m + m; }
  void push(const int32_t *lo_, const int32_t *hi_, const double *w_, double sign)
  {
    box.insert(box.end(), lo_, lo_ + m);
    box.insert(box.end(), hi_, hi_ + m);
    for (int k = 0; k < vs; ++k)
      w.push_back(sign * w_[k]);
  }
  bool key_less(size_t a, size_t b) const
  {
    const int32_t *ka = box.data() + a * 2 * (size_t)m, *kb = box.data() + b * 2 * (size_t)m;
    return std::lexicographical_compare(ka, ka + 2 * m, kb, kb + 2 * m);
  }
  bool key_equal(size_t a, const RectList &o, size_t b) const
  {
    return std::equal(
        box.data() + a * 2 * (size_t)m, box.data() + (a + 1) * 2 * (size_t)m, o.box.data() + b * 2 * (size_t)m
    );
  }
  // Sort by box and add the weights of equal boxes (in their original order, so the result is deterministic).
  void sort_merge()
  {
    size_t R = size();
    std::vector<size_t> perm(R);
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(perm.begin(), perm.end(), [&](size_t a, size_t b) { return key_less(a, b); });
    RectList out(m, vs);
    out.box.reserve(box.size());
    out.w.reserve(w.size());
    for (size_t i = 0; i < R; ++i) {
      size_t r = perm[i];
      if (out.size() > 0 && key_equal(r, out, out.size() - 1)) {
        double *dst = out.w.data() + (out.size() - 1) * vs;
        for (int k = 0; k < vs; ++k)
          dst[k] += w[r * vs + k];
      } else
        out.push(lo(r), hi(r), w.data() + r * vs, 1.0);
    }
    *this = std::move(out);
  }
  // Merge two sorted, duplicate-free lists; the weights of a come first in each sum.
  static RectList merge(const RectList &a, const RectList &b)
  {
    RectList out(a.m, a.vs);
    out.box.reserve(a.box.size() + b.box.size());
    out.w.reserve(a.w.size() + b.w.size());
    size_t i = 0, j = 0, na = a.size(), nb = b.size();
    const int m2 = 2 * a.m;
    while (i < na || j < nb) {
      int c;
      if (i == na)
        c = 1;
      else if (j == nb)
        c = -1;
      else {
        const int32_t *ka = a.box.data() + i * m2, *kb = b.box.data() + j * m2;
        c = std::lexicographical_compare(ka, ka + m2, kb, kb + m2) ? -1 : (std::equal(ka, ka + m2, kb) ? 0 : 1);
      }
      if (c < 0) {
        out.push(a.lo(i), a.hi(i), a.w.data() + i * a.vs, 1.0);
        ++i;
      } else if (c > 0) {
        out.push(b.lo(j), b.hi(j), b.w.data() + j * b.vs, 1.0);
        ++j;
      } else {
        out.push(a.lo(i), a.hi(i), a.w.data() + i * a.vs, 1.0);
        double *dst = out.w.data() + (out.size() - 1) * out.vs;
        for (int k = 0; k < out.vs; ++k)
          dst[k] += b.w[j * b.vs + k];
        ++i;
        ++j;
      }
    }
    return out;
  }
  void drop_zeros()
  {
    RectList out(m, vs);
    for (size_t r = 0; r < size(); ++r) {
      bool nz = false;
      for (int k = 0; k < vs; ++k)
        nz = nz || w[r * vs + k] != 0.0;
      if (nz)
        out.push(lo(r), hi(r), w.data() + r * vs, 1.0);
    }
    *this = std::move(out);
  }
};

// Point-query index over a RectList: a k-d decomposition of the cell space. A node adds the weights of the
// rectangles that cover its whole region; the others go to the child regions they overlap (rectangles across the
// split go to both). A leaf tests its rectangles one by one.
// The split of a node is chosen among the features and up to 15 quantiles of the rectangle boundaries inside its
// region, for the least expected number of rectangles a query meets below it (share of the region's cells on each
// side times the rectangles overlapping that side); a node is split only if that is below `stop` times its own
// rectangles. Nodes are split largest first, while the rectangle references held by the nodes stay within
// `budget` times the number of rectangles: this bounds the memory when rectangles overlap a lot (high orders).
struct KdIndex {
  struct Node {
    int dim = -1, split = 0, left = -1, right = -1, cst = -1, r0 = 0, r1 = 0;
  };
  std::vector<Node> nodes;
  std::vector<double> csts;
  std::vector<int> ids;

  // Afterwards rl holds only the rectangles the leaves test, in the order of the leaves (for locality); the others
  // are in the node constants.
  void build(
      RectList &rl, const std::vector<int32_t> &n_cells, size_t leaf_size = 32, double stop = 0.9, double budget = 16.0
  )
  {
    nodes.clear();
    csts.clear();
    ids.clear();
    const int m = rl.m, vs = rl.vs;
    struct Pending {
      int node;
      std::vector<int32_t> rlo, rhi;
      std::vector<int> idx;
    };
    std::vector<Pending> pend;
    auto cmp = [&](size_t a, size_t b) { return pend[a].idx.size() < pend[b].idx.size(); };
    std::vector<size_t> heap, free_slots;
    auto add_pending = [&](std::vector<int32_t> rlo, std::vector<int32_t> rhi, std::vector<int> idx) {
      int node = (int)nodes.size();
      nodes.emplace_back();
      size_t slot;
      if (!free_slots.empty()) {
        slot = free_slots.back();
        free_slots.pop_back();
      } else {
        slot = pend.size();
        pend.emplace_back();
      }
      pend[slot].node = node;
      pend[slot].rlo = std::move(rlo);
      pend[slot].rhi = std::move(rhi);
      pend[slot].idx = std::move(idx);
      heap.push_back(slot);
      std::push_heap(heap.begin(), heap.end(), cmp);
      return node;
    };
    {
      std::vector<int> idx(rl.size());
      std::iota(idx.begin(), idx.end(), 0);
      add_pending(std::vector<int32_t>(m, 0), std::vector<int32_t>(n_cells.begin(), n_cells.end()), std::move(idx));
    }
    const double max_refs = budget * (double)std::max<size_t>(1, rl.size());
    double refs = (double)rl.size(); // references held by pending nodes and leaves
    std::vector<double> cst(vs);
    std::vector<int> rest, left, right;
    std::vector<int32_t> bnd, los, his;
    while (!heap.empty()) {
      std::pop_heap(heap.begin(), heap.end(), cmp);
      size_t slot = heap.back();
      heap.pop_back();
      Pending &P = pend[slot];
      const int me = P.node;
      std::fill(cst.begin(), cst.end(), 0.0);
      bool has_cst = false;
      rest.clear();
      for (int r : P.idx) {
        const int32_t *lo = rl.lo(r), *hi = rl.hi(r);
        bool covers = true;
        for (int d = 0; d < m && covers; ++d)
          covers = lo[d] <= P.rlo[d] && hi[d] >= P.rhi[d];
        if (covers) {
          for (int k = 0; k < vs; ++k)
            cst[k] += rl.w[(size_t)r * vs + k];
          has_cst = true;
        } else
          rest.push_back(r);
      }
      refs -= (double)(P.idx.size() - rest.size());
      std::vector<int>().swap(P.idx);
      if (has_cst) {
        nodes[me].cst = (int)csts.size();
        csts.insert(csts.end(), cst.begin(), cst.end());
      }
      int best_d = -1;
      int32_t s = 0;
      double best_L = 0, best_R = 0;
      if (rest.size() > leaf_size) {
        double best_cost = stop * (double)rest.size();
        for (int d = 0; d < m; ++d) {
          bnd.clear();
          los.clear();
          his.clear();
          for (int r : rest) {
            int32_t a = rl.lo(r)[d], b = rl.hi(r)[d];
            los.push_back(a);
            his.push_back(b);
            if (a > P.rlo[d] && a < P.rhi[d])
              bnd.push_back(a);
            if (b > P.rlo[d] && b < P.rhi[d])
              bnd.push_back(b);
          }
          if (bnd.empty())
            continue;
          std::sort(bnd.begin(), bnd.end());
          bnd.erase(std::unique(bnd.begin(), bnd.end()), bnd.end());
          std::sort(los.begin(), los.end());
          std::sort(his.begin(), his.end());
          const size_t nq = std::min<size_t>(15, bnd.size());
          const double width = (double)(P.rhi[d] - P.rlo[d]);
          for (size_t q = 0; q < nq; ++q) {
            int32_t c = bnd[(size_t)(((double)q + 0.5) / (double)nq * (double)bnd.size())];
            double L = (double)(std::lower_bound(los.begin(), los.end(), c) - los.begin());
            double R = (double)(his.end() - std::upper_bound(his.begin(), his.end(), c));
            double pl = (double)(c - P.rlo[d]) / width;
            double cost = pl * L + (1.0 - pl) * R;
            if (cost < best_cost) {
              best_cost = cost;
              best_d = d;
              s = c;
              best_L = L;
              best_R = R;
            }
          }
        }
        if (best_d >= 0 && refs + best_L + best_R - (double)rest.size() > max_refs)
          best_d = -1;
      }
      if (best_d < 0) {
        nodes[me].r0 = (int)ids.size();
        ids.insert(ids.end(), rest.begin(), rest.end());
        nodes[me].r1 = (int)ids.size();
        std::vector<int32_t>().swap(P.rlo);
        std::vector<int32_t>().swap(P.rhi);
        free_slots.push_back(slot);
        continue;
      }
      left.clear();
      right.clear();
      for (int r : rest) {
        if (rl.lo(r)[best_d] < s)
          left.push_back(r);
        if (rl.hi(r)[best_d] > s)
          right.push_back(r);
      }
      refs += (double)(left.size() + right.size()) - (double)rest.size();
      nodes[me].dim = best_d;
      nodes[me].split = s;
      std::vector<int32_t> lrlo = P.rlo, lrhi = P.rhi, rrlo = std::move(P.rlo), rrhi = std::move(P.rhi);
      free_slots.push_back(slot); // P is not used below: add_pending may reuse or move the slot
      lrhi[best_d] = s;
      rrlo[best_d] = s;
      int l = add_pending(std::move(lrlo), std::move(lrhi), left);
      int rr = add_pending(std::move(rrlo), std::move(rrhi), right);
      nodes[me].left = l;
      nodes[me].right = rr;
    }
    std::vector<int> new_id(rl.size(), -1);
    RectList kept(m, vs);
    for (int &id : ids) {
      if (new_id[id] < 0) {
        new_id[id] = (int)kept.size();
        kept.push(rl.lo(id), rl.hi(id), rl.w.data() + (size_t)id * vs, 1.0);
      }
      id = new_id[id];
    }
    rl = std::move(kept);
    rl.box.shrink_to_fit(); // kept grew by push(): no spare capacity in the stored components
    rl.w.shrink_to_fit();
    ids.shrink_to_fit();
    nodes.shrink_to_fit();
    csts.shrink_to_fit();
  }

  // Adds the sum at the point with cells c[0..m) to out[0..vs).
  void add_at(const RectList &rl, const int32_t *c, double *out) const
  {
    if (nodes.empty())
      return;
    const int m = rl.m, vs = rl.vs;
    int n = 0;
    while (true) {
      const Node &nd = nodes[n];
      if (nd.cst >= 0)
        for (int k = 0; k < vs; ++k)
          out[k] += csts[nd.cst + k];
      if (nd.left < 0) {
        for (int i = nd.r0; i < nd.r1; ++i) {
          size_t r = (size_t)ids[i];
          const int32_t *lo = rl.lo(r), *hi = rl.hi(r);
          bool in = true;
          for (int d = 0; d < m && in; ++d)
            in = lo[d] <= c[d] && c[d] < hi[d];
          if (in)
            for (int k = 0; k < vs; ++k)
              out[k] += rl.w[r * vs + k];
        }
        return;
      }
      n = c[nd.dim] < nd.split ? nd.left : nd.right;
    }
  }
};

// The part of a component f_S that depends on the features j (a subset of S, given as positions in S).
struct Group {
  std::vector<int> pos;  // positions in S of the features in j
  std::vector<int> feat; // the features in j, 0-based
  RectList rects;
  std::vector<double> table; // |j| = 1: the sum per cell (n_cells * vs); the rectangles are not needed then
  KdIndex index;             // |j| >= 2
};

struct Component {
  std::vector<int> S;           // 1-based features, ascending
  std::vector<double> constant; // the j = {} terms
  std::vector<Group> groups;
};

// Accumulates (S, j) rectangle groups; used per family, then merged across families.
struct Accumulator {
  std::vector<double> intercept;
  // S (1-based, ascending) -> (constant, mask of j over the positions of S -> rectangles)
  std::map<std::vector<int>, std::pair<std::vector<double>, std::map<int, RectList>>> comps;

  void sort_merge()
  {
    for (auto &c : comps)
      for (auto &g : c.second.second)
        g.second.sort_merge();
  }
  // this += o (both sorted and merged); the terms of this come first in every sum
  void absorb(Accumulator &o)
  {
    if (intercept.empty())
      intercept = o.intercept;
    else
      for (size_t k = 0; k < intercept.size(); ++k)
        intercept[k] += o.intercept[k];
    for (auto &c : o.comps) {
      auto it = comps.find(c.first);
      if (it == comps.end()) {
        comps.emplace(c.first, std::move(c.second));
        continue;
      }
      auto &mine = it->second;
      for (size_t k = 0; k < mine.first.size(); ++k)
        mine.first[k] += c.second.first[k];
      for (auto &g : c.second.second) {
        auto gt = mine.second.find(g.first);
        if (gt == mine.second.end())
          mine.second.emplace(g.first, std::move(g.second));
        else {
          gt->second = RectList::merge(gt->second, g.second);
          g.second = RectList();
        }
      }
    }
    o.comps.clear();
  }
};

// All purified components of a forest (summed over the families, not yet divided by the number of trees).
struct PurifiedForest {
  std::vector<std::vector<double>> cuts; // per feature (0-based)
  std::vector<double> intercept;
  std::map<std::vector<int>, Component> comps;
  size_t n_rects = 0;
  int maxp = 0; // maxp_interaction of the purify() call (normalised: 1..p)

  int32_t cell(int feature, double x) const
  {
    const auto &c = cuts[feature];
    return (int32_t)(std::upper_bound(c.begin(), c.end(), x) - c.begin());
  }

  // Checks a restored store (e.g. from a marshaled forest) for consistency; throws std::invalid_argument.
  void validate(int vs) const
  {
    auto fail = [](const std::string &what) { throw std::invalid_argument("Corrupt purified components: " + what); };
    const int p = (int)cuts.size();
    if ((int)intercept.size() != vs)
      fail("intercept length");
    for (const auto &kv : comps) {
      const Component &c = kv.second;
      if (c.S != kv.first || c.S.empty() || (int)c.constant.size() != vs)
        fail("component header");
      for (size_t i = 0; i < c.S.size(); ++i)
        if (c.S[i] < 1 || c.S[i] > p || (i > 0 && c.S[i] <= c.S[i - 1]))
          fail("component features");
      for (const Group &g : c.groups) {
        const int m = (int)g.feat.size();
        if (m < 1 || g.pos.size() != g.feat.size())
          fail("group features");
        for (int q = 0; q < m; ++q)
          if (g.pos[q] < 0 || g.pos[q] >= (int)c.S.size() || c.S[g.pos[q]] - 1 != g.feat[q])
            fail("group features");
        if (m == 1 && !g.table.empty()) {
          if (g.table.size() != (cuts[g.feat[0]].size() + 1) * (size_t)vs)
            fail("table size");
          continue;
        }
        const RectList &rl = g.rects;
        const size_t R = rl.size();
        if (rl.m != m || rl.vs != vs || rl.box.size() != R * 2 * (size_t)m || rl.w.size() != R * (size_t)vs)
          fail("rectangle list size");
        for (size_t r = 0; r < R; ++r)
          for (int q = 0; q < m; ++q) {
            const int32_t a = rl.lo(r)[q], b = rl.hi(r)[q], nc = (int32_t)cuts[g.feat[q]].size() + 1;
            if (a < 0 || a >= b || b > nc)
              fail("rectangle bounds");
          }
        const KdIndex &ix = g.index;
        const int nn = (int)ix.nodes.size();
        for (int i = 0; i < nn; ++i) {
          const KdIndex::Node &nd = ix.nodes[i];
          if ((nd.left < 0) != (nd.right < 0))
            fail("index node");
          if (nd.left >= 0 &&
              (nd.left <= i || nd.right <= i || nd.left >= nn || nd.right >= nn || nd.dim < 0 || nd.dim >= m))
            fail("index node");
          if (nd.cst < -1 || (nd.cst >= 0 && (size_t)nd.cst + vs > ix.csts.size()))
            fail("index constant");
          if (nd.r0 < 0 || nd.r0 > nd.r1 || (size_t)nd.r1 > ix.ids.size())
            fail("index leaf");
        }
        for (int id : ix.ids)
          if (id < 0 || (size_t)id >= R)
            fail("index leaf");
      }
    }
  }

  // Takes the merged groups out of acc (emptying it) and builds the query tables and indices, the latter on up to
  // build_threads threads.
  void build_from(Accumulator &acc, int vs, unsigned build_threads = 1)
  {
    intercept = acc.intercept;
    if (intercept.empty())
      intercept.assign(vs, 0.0);
    comps.clear();
    n_rects = 0;
    std::vector<Group *> todo;
    for (auto it = acc.comps.begin(); it != acc.comps.end(); it = acc.comps.erase(it)) {
      Component comp;
      comp.S = it->first;
      comp.constant = it->second.first;
      for (auto &g : it->second.second) {
        Group grp;
        for (size_t b = 0; b < comp.S.size(); ++b)
          if (g.first & (1 << b)) {
            grp.pos.push_back((int)b);
            grp.feat.push_back(comp.S[b] - 1);
          }
        grp.rects = std::move(g.second);
        grp.rects.drop_zeros();
        if (grp.rects.size() == 0)
          continue;
        n_rects += grp.rects.size();
        comp.groups.push_back(std::move(grp));
      }
      auto &stored = comps.emplace(comp.S, std::move(comp)).first->second;
      for (auto &grp : stored.groups)
        todo.push_back(&grp);
    }
    std::atomic<size_t> next(0);
    auto work = [&]() {
      for (size_t i = next++; i < todo.size(); i = next++)
        build_group(*todo[i], vs);
    };
    build_threads = std::max(1u, std::min<unsigned>(build_threads, (unsigned)todo.size()));
    if (build_threads <= 1)
      work();
    else {
      std::vector<std::thread> threads;
      for (unsigned t = 0; t < build_threads; ++t)
        threads.emplace_back(work);
      for (auto &th : threads)
        th.join();
    }
  }

  // |j| = 1: the sum per cell (the rectangles are dropped); |j| >= 2: the k-d index.
  void build_group(Group &grp, int vs) const
  {
    if (grp.feat.size() == 1) {
      size_t nc = cuts[grp.feat[0]].size() + 1;
      std::vector<double> diff((nc + 1) * vs, 0.0);
      for (size_t r = 0; r < grp.rects.size(); ++r)
        for (int k = 0; k < vs; ++k) {
          diff[(size_t)grp.rects.lo(r)[0] * vs + k] += grp.rects.w[r * vs + k];
          diff[(size_t)grp.rects.hi(r)[0] * vs + k] -= grp.rects.w[r * vs + k];
        }
      grp.table.assign(nc * vs, 0.0);
      std::vector<double> run(vs, 0.0);
      for (size_t c = 0; c < nc; ++c)
        for (int k = 0; k < vs; ++k) {
          run[k] += diff[c * vs + k];
          grp.table[c * vs + k] = run[k];
        }
      grp.rects = RectList(1, vs);
    } else {
      std::vector<int32_t> n_cells;
      for (int f : grp.feat)
        n_cells.push_back((int32_t)cuts[f].size() + 1);
      grp.index.build(grp.rects, n_cells);
    }
  }

  // Adds f_S at n points to out (n x vs, row-major). cells is n x stride, row-major, indexed by 0-based feature.
  // Loops over the groups outside and the points inside, so that a group's index stays in cache.
  void
  add_component_rows(const Component &comp, const int32_t *cells, size_t stride, size_t n, double *out, int vs) const
  {
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < vs; ++k)
        out[i * vs + k] += comp.constant[k];
    std::vector<int32_t> c;
    for (const auto &g : comp.groups) {
      if (!g.table.empty()) {
        const int f = g.feat[0];
        for (size_t i = 0; i < n; ++i) {
          const double *t = g.table.data() + (size_t)cells[i * stride + f] * vs;
          for (int k = 0; k < vs; ++k)
            out[i * vs + k] += t[k];
        }
        continue;
      }
      const size_t m = g.feat.size();
      c.resize(m);
      for (size_t i = 0; i < n; ++i) {
        for (size_t q = 0; q < m; ++q)
          c[q] = cells[i * stride + g.feat[q]];
        g.index.add_at(g.rects, c.data(), out + i * vs);
      }
    }
  }

  // Adds f_S at the point whose cells (indexed by 0-based feature) are `cells` to out[0..vs).
  void add_component(const Component &comp, const int32_t *cells, double *out, int vs) const
  {
    for (int k = 0; k < vs; ++k)
      out[k] += comp.constant[k];
    int32_t c[64];
    std::vector<int32_t> cbig;
    for (const auto &g : comp.groups) {
      if (!g.table.empty()) {
        const double *t = g.table.data() + (size_t)cells[g.feat[0]] * vs;
        for (int k = 0; k < vs; ++k)
          out[k] += t[k];
        continue;
      }
      int32_t *cp = c;
      if (g.feat.size() > 64) {
        cbig.resize(g.feat.size());
        cp = cbig.data();
      }
      for (size_t i = 0; i < g.feat.size(); ++i)
        cp[i] = cells[g.feat[i]];
      g.index.add_at(g.rects, cp, out);
    }
  }
};
} // namespace rpf_rect

#endif
