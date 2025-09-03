#ifndef RPF_KDTREE_HPP
#define RPF_KDTREE_HPP

#include <vector>
#include <limits>
#include <algorithm>
#include <memory>

// Lightweight KD-tree for orthogonal range counts.
// - Header-only to avoid build system changes
// - Supports arbitrary dimensionality
// - Query provides constraints only for a subset of dimensions; others are unconstrained

namespace rpf_kd
{
  struct RangeConstraint
  {
    int dim;            // 0-based feature index
    double left;        // inclusive lower bound
    double right;       // exclusive upper bound
  };

  namespace detail
  {
    struct Node
    {
      // Bounding box for quick acceptance/rejection
      std::vector<double> minv;
      std::vector<double> maxv;
      int axis = -1;          // split axis; -1 means leaf
      double split_value = 0; // split threshold
      size_t size = 0;        // number of points in subtree
      std::unique_ptr<Node> left;
      std::unique_ptr<Node> right;
      std::vector<int> idxs;  // indices when leaf
    };
  }

  class KDTree
  {
  public:
    KDTree() = default;

    KDTree(const std::vector<std::vector<double>> *X_ptr,
           const std::vector<int> &all_indices,
           int dims,
           size_t leaf_size = 32)
    {
      build(X_ptr, all_indices, dims, leaf_size);
    }

    void build(const std::vector<std::vector<double>> *X_ptr,
               const std::vector<int> &all_indices,
               int dims,
               size_t leaf_size = 32)
    {
      X_ = X_ptr;
      dims_ = dims;
      leaf_size_ = leaf_size;
      root_ = build_recursive(all_indices);
    }

    // Count number of points with constraints on a subset of dims
    size_t range_count(const std::vector<RangeConstraint> &constraints) const
    {
      return range_count_recursive(root_.get(), constraints);
    }

  private:
    const std::vector<std::vector<double>> *X_ = nullptr;
    int dims_ = 0;
    size_t leaf_size_ = 32;
    std::unique_ptr<detail::Node> root_;

    std::unique_ptr<detail::Node> build_recursive(const std::vector<int> &idxs)
    {
      auto node = std::make_unique<detail::Node>();
      node->size = idxs.size();
      node->minv.assign(dims_, std::numeric_limits<double>::infinity());
      node->maxv.assign(dims_, -std::numeric_limits<double>::infinity());
      for (int i : idxs)
      {
        for (int d = 0; d < dims_; ++d)
        {
          double v = (*X_)[i][d];
          if (v < node->minv[d]) node->minv[d] = v;
          if (v > node->maxv[d]) node->maxv[d] = v;
        }
      }

      if (idxs.size() <= leaf_size_)
      {
        node->axis = -1; node->idxs = idxs; return node;
      }

      // Choose split axis by widest spread
      int axis = 0; double best_span = -1.0;
      for (int d = 0; d < dims_; ++d)
      {
        double span = node->maxv[d] - node->minv[d];
        if (span > best_span) { best_span = span; axis = d; }
      }
      node->axis = axis;

      // Median split on chosen axis
      std::vector<int> left_idxs, right_idxs; left_idxs.reserve(idxs.size()); right_idxs.reserve(idxs.size());
      std::vector<int> tmp = idxs;
      size_t mid = tmp.size() / 2;
      std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end(), [&](int a, int b){ return (*X_)[a][axis] < (*X_)[b][axis]; });
      double split = (*X_)[tmp[mid]][axis];
      node->split_value = split;
      for (int i : idxs)
      {
        if ((*X_)[i][axis] < split) left_idxs.push_back(i); else right_idxs.push_back(i);
      }
      if (left_idxs.empty() || right_idxs.empty())
      {
        // Fallback: make leaf if degenerate split
        node->axis = -1; node->idxs = idxs; return node;
      }
      node->left = build_recursive(left_idxs);
      node->right = build_recursive(right_idxs);
      return node;
    }

    static inline bool box_outside(const std::vector<double> &minv, const std::vector<double> &maxv,
                                   const std::vector<RangeConstraint> &C)
    {
      for (const auto &rc : C)
      {
        if (maxv[rc.dim] <= rc.left) return true;
        if (minv[rc.dim] >= rc.right) return true;
      }
      return false;
    }

    static inline bool box_inside(const std::vector<double> &minv, const std::vector<double> &maxv,
                                  const std::vector<RangeConstraint> &C)
    {
      for (const auto &rc : C)
      {
        if (minv[rc.dim] < rc.left) return false;
        if (maxv[rc.dim] > rc.right) return false;
      }
      return true;
    }

    size_t range_count_recursive(const detail::Node *node, const std::vector<RangeConstraint> &C) const
    {
      if (!node) return 0;
      if (!C.empty())
      {
        if (box_outside(node->minv, node->maxv, C)) return 0;
        if (box_inside(node->minv, node->maxv, C)) return node->size;
      }
      if (node->axis == -1)
      {
        size_t cnt = 0;
        for (int i : node->idxs)
        {
          bool inside = true;
          for (const auto &rc : C)
          {
            double v = (*X_)[i][rc.dim];
            if (!(v >= rc.left && v < rc.right)) { inside = false; break; }
          }
          if (inside) ++cnt;
        }
        return cnt;
      }
      return range_count_recursive(node->left.get(), C) + range_count_recursive(node->right.get(), C);
    }
  };
}

#endif // RPF_KDTREE_HPP


