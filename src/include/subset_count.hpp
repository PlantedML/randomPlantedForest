// Counts of training rows in axis-parallel boxes, for purify mode 2's leaf weights (cnt(box on U) / tot(U)).
//
// rpf_kd::KDTree answers a box over a subset U of the p features from one tree over all p features, which prunes only
// at the levels that split a feature of U. SubsetCounter instead answers a box over U from a k-d tree over the rows'
// coordinates on U alone (built when U is first queried and kept), a box over one feature by binary search, and drops
// the constraints that every row meets. The counts are the same integers as rpf_kd::KDTree::range_count's (a row is
// in [left, right) on every constrained feature), so purifications that use it are identical.
#ifndef RPF_SUBSET_COUNT_HPP
#define RPF_SUBSET_COUNT_HPP

#include "kdtree.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <vector>

namespace rpf_cnt {
// A k-d tree over n points in m dimensions (stored in tree order), for counts of the points in boxes [l, r)
struct BoxKD {
  int m = 0;
  std::vector<double> pts;          // n x m
  std::vector<double> nlo, nhi;     // per node: the bounding box of its points (m each)
  std::vector<int32_t> first, last; // per node: its points [first, last)
  std::vector<int32_t> left, right; // per node: children (-1: a leaf)

  void build(const std::vector<std::vector<double>> &X, const std::vector<int> &U, size_t leaf_size = 16)
  {
    m = (int)U.size();
    const size_t n = X.size();
    std::vector<int32_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> c(n * (size_t)m);
    for (size_t i = 0; i < n; ++i)
      for (int q = 0; q < m; ++q)
        c[i * m + q] = X[i][(size_t)U[q]];
    struct Job {
      int node, a, b;
    };
    std::vector<Job> stack;
    auto new_node = [&](int a, int b) {
      int id = (int)first.size();
      first.push_back(a);
      last.push_back(b);
      left.push_back(-1);
      right.push_back(-1);
      for (int q = 0; q < m; ++q) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (int k = a; k < b; ++k) {
          double v = c[(size_t)idx[k] * m + q];
          lo = std::min(lo, v);
          hi = std::max(hi, v);
        }
        nlo.push_back(lo);
        nhi.push_back(hi);
      }
      return id;
    };
    stack.push_back({new_node(0, (int)n), 0, (int)n});
    while (!stack.empty()) {
      Job j = stack.back();
      stack.pop_back();
      if ((size_t)(j.b - j.a) <= leaf_size)
        continue;
      int axis = 0;
      double span = -1;
      for (int q = 0; q < m; ++q) {
        double s = nhi[(size_t)j.node * m + q] - nlo[(size_t)j.node * m + q];
        if (s > span) {
          span = s;
          axis = q;
        }
      }
      if (span <= 0)
        continue;
      const int mid = (j.a + j.b) / 2;
      std::nth_element(idx.begin() + j.a, idx.begin() + mid, idx.begin() + j.b, [&](int32_t x, int32_t y) {
        return c[(size_t)x * m + axis] < c[(size_t)y * m + axis];
      });
      int l = new_node(j.a, mid), r = new_node(mid, j.b);
      left[(size_t)j.node] = l;
      right[(size_t)j.node] = r;
      stack.push_back({l, j.a, mid});
      stack.push_back({r, mid, j.b});
    }
    pts.resize(n * (size_t)m);
    for (size_t k = 0; k < n; ++k)
      for (int q = 0; q < m; ++q)
        pts[k * m + q] = c[(size_t)idx[k] * m + q];
  }

  // The number of points with l[q] <= x[q] < r[q] for every q
  size_t count(const double *l, const double *r) const
  {
    size_t cnt = 0;
    int stack[128];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
      const int nd = stack[--top];
      const double *lo = &nlo[(size_t)nd * m], *hi = &nhi[(size_t)nd * m];
      bool out = false, in = true;
      for (int q = 0; q < m; ++q) {
        if (hi[q] < l[q] || lo[q] >= r[q]) {
          out = true;
          break;
        }
        if (lo[q] < l[q] || hi[q] >= r[q])
          in = false;
      }
      if (out)
        continue;
      if (in) {
        cnt += (size_t)(last[nd] - first[nd]);
        continue;
      }
      if (left[nd] < 0) {
        cnt += count_leaf(nd, l, r);
        continue;
      }
      if (top + 2 > 128) { // deeper than any balanced tree here; finish this subtree by recursion
        cnt += count_from(left[nd], l, r) + count_from(right[nd], l, r);
        continue;
      }
      stack[top++] = left[nd];
      stack[top++] = right[nd];
    }
    return cnt;
  }

private:
  size_t count_leaf(int nd, const double *l, const double *r) const
  {
    size_t cnt = 0;
    for (int k = first[nd]; k < last[nd]; ++k) {
      const double *x = &pts[(size_t)k * m];
      bool ok = true;
      for (int q = 0; q < m && ok; ++q)
        ok = x[q] >= l[q] && x[q] < r[q];
      cnt += ok ? 1 : 0;
    }
    return cnt;
  }

  size_t count_from(int nd, const double *l, const double *r) const
  {
    const double *lo = &nlo[(size_t)nd * m], *hi = &nhi[(size_t)nd * m];
    bool in = true;
    for (int q = 0; q < m; ++q) {
      if (hi[q] < l[q] || lo[q] >= r[q])
        return 0;
      if (lo[q] < l[q] || hi[q] >= r[q])
        in = false;
    }
    if (in)
      return (size_t)(last[nd] - first[nd]);
    if (left[nd] < 0)
      return count_leaf(nd, l, r);
    return count_from(left[nd], l, r) + count_from(right[nd], l, r);
  }
};

class SubsetCounter {
public:
  SubsetCounter(const std::vector<std::vector<double>> *X, int p) : X_(X)
  {
    const size_t n = X->size();
    sorted_.assign((size_t)p, std::vector<double>(n));
    for (int d = 0; d < p; ++d) {
      for (size_t i = 0; i < n; ++i)
        sorted_[(size_t)d][i] = (*X)[i][(size_t)d];
      std::sort(sorted_[(size_t)d].begin(), sorted_[(size_t)d].end());
    }
  }

  // The same count as rpf_kd::KDTree::range_count(cons) over the rows of X (thread-safe)
  size_t range_count(const std::vector<rpf_kd::RangeConstraint> &cons) const
  {
    const size_t n = X_->size();
    if (n == 0)
      return 0;
    // the constraints every row meets are dropped; the others sorted by feature
    std::vector<const rpf_kd::RangeConstraint *> keep;
    for (const auto &c : cons) {
      const auto &s = sorted_[(size_t)c.dim];
      if (!(c.left < c.right))
        return 0;
      if (c.left <= s.front() && c.right > s.back())
        continue;
      keep.push_back(&c);
    }
    if (keep.empty())
      return n;
    if (keep.size() == 1) {
      const auto &s = sorted_[(size_t)keep[0]->dim];
      return (size_t)(std::lower_bound(s.begin(), s.end(), keep[0]->right) -
                      std::lower_bound(s.begin(), s.end(), keep[0]->left));
    }
    std::sort(keep.begin(), keep.end(), [](const rpf_kd::RangeConstraint *a, const rpf_kd::RangeConstraint *b) {
      return a->dim < b->dim;
    });
    std::vector<int> U;
    std::vector<double> l, r;
    for (const auto *c : keep) {
      if (!U.empty() && U.back() == c->dim) { // the same feature twice: intersect
        l.back() = std::max(l.back(), c->left);
        r.back() = std::min(r.back(), c->right);
        if (!(l.back() < r.back()))
          return 0;
        continue;
      }
      U.push_back(c->dim);
      l.push_back(c->left);
      r.push_back(c->right);
    }
    return tree(U).count(l.data(), r.data());
  }

private:
  const std::vector<std::vector<double>> *X_;
  std::vector<std::vector<double>> sorted_; // per feature: the training rows' values, ascending
  mutable std::mutex mu_;
  mutable std::map<std::vector<int>, std::shared_ptr<BoxKD>> trees_;

  // The k-d tree of the features U (0-based, ascending), built on first use
  const BoxKD &tree(const std::vector<int> &U) const
  {
    std::lock_guard<std::mutex> lk(mu_);
    auto it = trees_.find(U);
    if (it != trees_.end())
      return *it->second;
    auto t = std::make_shared<BoxKD>();
    t->build(*X_, U);
    return *trees_.emplace(U, t).first->second;
  }
};
} // namespace rpf_cnt

#endif
