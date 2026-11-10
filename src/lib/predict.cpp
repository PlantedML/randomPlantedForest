// Prediction entry points split out from rpf.cpp for readability and reuse.
#include "rpf.hpp"
#include <algorithm>
#include <cstdint>
#include <exception>
#include <iterator>
#include <limits>
#include <thread>

// Prediction from the rectangle components of purify(mode = 2): {-1} the intercept, {0} the whole forest, otherwise
// the component S = component_index, with X holding the features of S in ascending order.
std::vector<double>
RandomPlantedForest::predict_single_rect(const std::vector<double> &X, const std::set<int> &component_index) const
{
  const rpf_rect::PurifiedForest &pf = *rect_forest_;
  const int vs = (int)value_size;
  std::vector<double> total_res(value_size, 0);
  if (component_index == std::set<int>{-1}) {
    total_res = pf.intercept;
  } else if (component_index == std::set<int>{0}) {
    std::vector<int32_t> cells(feature_size);
    for (int d = 0; d < feature_size; ++d)
      cells[d] = pf.cell(d, X[d]);
    total_res = pf.intercept;
    for (const auto &c : pf.comps)
      pf.add_component(c.second, cells.data(), total_res.data(), vs);
  } else {
    std::vector<int> S(component_index.begin(), component_index.end());
    auto it = pf.comps.find(S);
    if (it != pf.comps.end()) {
      std::vector<int32_t> cells(feature_size, 0);
      for (size_t i = 0; i < S.size(); ++i)
        cells[S[i] - 1] = pf.cell(S[i] - 1, X[i]);
      pf.add_component(it->second, cells.data(), total_res.data(), vs);
    }
  }
  return total_res / n_trees;
}

// predict single feature vector
std::vector<double> RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index)
{
  std::vector<double> total_res = std::vector<double>(value_size, 0);

  // After a capped purify() the purified components hold only the low orders: the whole forest is predicted from the
  // leaves, which purification leaves unchanged
  const bool from_leaves = !purified || (purify_capped && component_index == std::set<int>{0});

  if (!from_leaves && rect_forest_)
    return predict_single_rect(X, component_index);

  if (from_leaves) {
    // consider all components
    if (component_index == std::set<int>{0}) {
      for (auto &tree_family : this->tree_families) {
        for (auto &tree : tree_family) {
          for (auto &leaf : tree.second->leaves) {
            bool valid = true;
            for (auto &dim : tree.first) {
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[std::max(0, dim - 1)] ||
                     leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) &&
                    (leaf.intervals[std::max(0, dim - 1)].second > X[std::max(0, dim - 1)] ||
                     leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)]))) {
                valid = false;
                break;
              }
            }
            if (valid) {
              for (size_t p = 0; p < value_size && p < leaf.value.size(); ++p) {
                total_res[p] += leaf.value[p];
              }
            }
          }
        }
      }
    } else { // choose components for prediction
      for (auto &tree_family : this->tree_families) {
        for (auto &tree : tree_family) {
          // only consider trees with same dimensions as component_index
          if (tree.first != component_index)
            continue;

          std::vector<int> dims;
          for (auto dim : tree.first) {
            dims.push_back(dim);
          }

          for (auto &leaf : tree.second->leaves) {
            bool valid = true;
            for (unsigned int i = 0; i < dims.size(); ++i) {
              int dim = dims[i];
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[i] ||
                     leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) &&
                    (leaf.intervals[std::max(0, dim - 1)].second > X[i] ||
                     leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)]))) {
                valid = false;
                break;
              }
            }
            if (valid) {
              for (size_t p = 0; p < value_size && p < leaf.value.size(); ++p) {
                total_res[p] += leaf.value[p];
              }
            }
          }
        }
      }
    }
  } else {
    if (component_index == std::set<int>{-1}) {
      for (auto &tree_family : this->tree_families) {
        for (auto &tree : tree_family) {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0}) {
            leaf_index = std::vector<int>(tree.first.size(), 0);

            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p) {
              total_res[p] += vals[p];
            }
          }
        }
      }
    } else if (component_index == std::set<int>{0}) {
      for (auto &tree_family : this->tree_families) {
        for (auto &tree : tree_family) {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0}) {
            leaf_index = std::vector<int>(tree.first.size(), 0);
          } else {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index) {
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // convert to 0-based original feature index
              }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2) {
                leaf_index[dim_index] = 0;
                continue;
              }
              // Use the original feature index into X, not the position within the tree's dim set
              auto it = std::upper_bound(bounds.begin(), bounds.end(), X[dim]);
              int c = static_cast<int>(std::distance(bounds.begin(), it));
              leaf_index[dim_index] = std::min(std::max(0, c - 1), (int)bounds.size() - 2);
            }
          }
          for (int &index : leaf_index)
            index = std::max(0, index);
          {
            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p) {
              total_res[p] += vals[p];
            }
          }
        }
      }
    } else {
      for (auto &tree_family : this->tree_families) {
        for (auto &tree : tree_family) {
          if (tree.first != component_index)
            continue;
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0}) {
            leaf_index = std::vector<int>(tree.first.size(), 0);
          } else {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index) {
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // 0-based original feature index for bounds lookup only
              }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2) {
                leaf_index[dim_index] = 0;
                continue;
              }
              // For component-specific prediction, X contains only the selected dims in ascending order.
              // Use the position within the selected dims (dim_index) to read the value.
              auto it = std::upper_bound(bounds.begin(), bounds.end(), X[dim_index]);
              int c = static_cast<int>(std::distance(bounds.begin(), it));
              leaf_index[dim_index] = std::min(std::max(0, c - 1), (int)bounds.size() - 2);
            }
          }
          for (int &index : leaf_index)
            index = std::max(0, index);
          {
            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p) {
              total_res[p] += vals[p];
            }
          }
        }
      }
    }
  }

  return total_res / n_trees;
}

namespace {
// One tree's leaves in contiguous arrays, so the per-row scan streams through memory.
struct FlatTree {
  std::vector<int> cols;      // column of X compared against each tree dimension
  std::vector<double> lo, hi; // n_leaves x cols.size(), row-major; open bounds stored as -/+inf
  std::vector<double> values; // n_leaves x value_size
  size_t n_leaves = 0;
  // The leaf index, if built (LeafIndexer): a k-d tree. An inner node sends x to node a if x[dim] < split, else
  // (also for NaN) to node b. A bucket (dim < 0) holds the leaves items[a, b), in ascending order.
  struct Node {
    double split;
    uint32_t a, b;
    int dim;
  };
  std::vector<Node> nodes;
  std::vector<uint32_t> items;
};

// The leaf rule of predict_single(): an interval edge at the training bound is open (stored as -/+inf), so new data
// outside the training range still lands in a leaf. The explicit inf checks keep NaN/inf inputs matching it.
inline bool in_leaf(const double *lo, const double *hi, const double *x, size_t nd)
{
  const double inf = std::numeric_limits<double>::infinity();
  for (size_t j = 0; j < nd; ++j)
    if (!((lo[j] <= x[j] || lo[j] == -inf) && (x[j] < hi[j] || hi[j] == inf)))
      return false;
  return true;
}

// From kIndexRows rows on, the trees of more than kIndexLeaves leaves are indexed before the rows are predicted. With
// fewer rows building the index costs more than it saves, and in smaller trees the scan is as fast (CA0, 100
// families). A bucket holds kIndexBucket leaves or fewer unless no cut separates them.
constexpr int kIndexRows = 512;
constexpr size_t kIndexLeaves = 64;
constexpr size_t kIndexBucket = 8;
constexpr size_t kIndexCopies = 4;  // at most this many index entries per leaf
constexpr size_t kIndexSample = 16; // a node's cut is taken from the bounds of at most this many of its leaves
constexpr int kIndexDepth = 48;

// The leaf index of a tree: a k-d tree that cuts the space, not the set of leaves. A node's cut on dimension j at s
// sends the leaves with lo < s to the left node and those with hi > s to the right one (lo == -inf and hi == inf
// count as such), so a leaf across the cut goes to both. A row x goes left if x[j] < s, else right. Every leaf
// that in_leaf() accepts for x is then in the bucket that x reaches: lo <= x < s gives lo < s; x >= s with x < hi
// gives hi > s, and x = NaN or +inf needs hi == inf. The leaves stay in ascending order in every node, so a bucket's
// matching leaves are added in the scan's order.
// The cut on each dimension is the median of the leaf bounds inside the node's cell (of a sample of its leaves); the
// dimension whose cut leaves the fewest leaves on the larger side is taken. A node becomes a bucket if no cut makes
// both sides smaller, or if the leaves across the cut would exceed the node's share of the copies (kIndexCopies per
// leaf in all, shared between two nodes in proportion to their leaves). Built per call, freed with the FlatTree.
class LeafIndexer {
public:
  // Index the tree (if it has more than kIndexLeaves leaves). One indexer per thread: its scratch is reused.
  void index(FlatTree &tree)
  {
    const size_t L = tree.n_leaves;
    nd = tree.cols.size();
    if (L <= kIndexLeaves || nd == 0 || L > std::numeric_limits<uint32_t>::max() / kIndexCopies)
      return;
    ft = &tree;
    const double inf = std::numeric_limits<double>::infinity();
    cell_lo.assign(nd, -inf);
    cell_hi.assign(nd, inf);
    scratch.resize(2 * kIndexDepth + 1);
    std::vector<uint32_t> &all = scratch[2 * kIndexDepth];
    all.resize(L);
    for (size_t l = 0; l < L; ++l)
      all[l] = (uint32_t)l;
    build(all, (kIndexCopies - 1) * L, 0);
    tree.nodes.shrink_to_fit();
    tree.items.shrink_to_fit();
  }

private:
  FlatTree *ft = nullptr;
  size_t nd = 0;
  std::vector<double> cell_lo, cell_hi, bounds;
  std::vector<std::vector<uint32_t>> scratch; // the leaves of the two children of the node at each depth

  uint32_t build(const std::vector<uint32_t> &leaves, size_t spare, int depth)
  {
    const uint32_t id = (uint32_t)ft->nodes.size();
    ft->nodes.push_back(FlatTree::Node{0.0, 0, 0, -1});
    const size_t m = leaves.size();
    const double *lo = ft->lo.data(), *hi = ft->hi.data();
    size_t best_j = nd, best_large = m, best_sum = 0;
    double best_s = 0.0;
    for (size_t j = 0; m > kIndexBucket && depth < kIndexDepth && j < nd; ++j) {
      bounds.clear();
      const size_t n_sample = std::min(m, kIndexSample);
      for (size_t i = 0; i < n_sample; ++i) {
        const size_t l = leaves[i * m / n_sample];
        const double a = lo[l * nd + j], b = hi[l * nd + j];
        if (a > cell_lo[j] && a < cell_hi[j])
          bounds.push_back(a);
        if (b > cell_lo[j] && b < cell_hi[j])
          bounds.push_back(b);
      }
      if (bounds.empty())
        continue;
      const auto mid = bounds.begin() + (std::ptrdiff_t)(bounds.size() / 2);
      std::nth_element(bounds.begin(), mid, bounds.end());
      const double s = *mid;
      size_t n_left = 0, n_right = 0;
      for (uint32_t l : leaves) {
        n_left += lo[l * nd + j] < s;
        n_right += hi[l * nd + j] > s;
      }
      const size_t large = std::max(n_left, n_right), sum = n_left + n_right;
      if (large < best_large || (large == best_large && best_j < nd && sum < best_sum)) {
        best_j = j;
        best_s = s;
        best_large = large;
        best_sum = sum;
      }
    }
    const size_t copies = best_sum > m ? best_sum - m : 0; // the leaves across the cut
    if (best_j == nd || copies > spare) {
      ft->nodes[id].a = (uint32_t)ft->items.size();
      ft->items.insert(ft->items.end(), leaves.begin(), leaves.end());
      ft->nodes[id].b = (uint32_t)ft->items.size();
      return id;
    }
    spare -= copies;

    std::vector<uint32_t> &left = scratch[2 * (size_t)depth], &right = scratch[2 * (size_t)depth + 1];
    left.clear();
    right.clear();
    for (uint32_t l : leaves) {
      if (lo[l * nd + best_j] < best_s)
        left.push_back(l);
      if (hi[l * nd + best_j] > best_s)
        right.push_back(l);
    }
    ft->nodes[id].dim = (int)best_j;
    ft->nodes[id].split = best_s;
    const size_t spare_left = (size_t)((double)spare * (double)left.size() / (double)best_sum);
    const double upper = cell_hi[best_j], lower = cell_lo[best_j];
    cell_hi[best_j] = best_s;
    const uint32_t a = build(left, spare_left, depth + 1);
    cell_hi[best_j] = upper;
    cell_lo[best_j] = best_s;
    const uint32_t b = build(right, spare - spare_left, depth + 1);
    cell_lo[best_j] = lower;
    ft->nodes[id].a = a;
    ft->nodes[id].b = b;
    return id;
  }
};

// Runs f(begin, end) on contiguous row ranges, one per thread.
template <class F> void parallel_rows(int n, unsigned int threads, F f)
{
  threads = std::max(1u, std::min<unsigned int>(threads, (unsigned int)n));
  if (threads == 1) {
    f(0, n);
    return;
  }
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(threads);
  int chunk = (n + (int)threads - 1) / (int)threads;
  for (unsigned int t = 0; t < threads; ++t) {
    int begin = (int)t * chunk, end = std::min(n, begin + chunk);
    pool.emplace_back([&f, &errors, t, begin, end] {
      // an exception escaping a std::thread would terminate R
      try {
        f(begin, end);
      } catch (...) {
        errors[t] = std::current_exception();
      }
    });
  }
  for (auto &th : pool)
    th.join();
  for (auto &e : errors)
    if (e)
      std::rethrow_exception(e);
}
} // namespace

// predict multiple feature vectors
std::vector<double>
RandomPlantedForest::predict_matrix(const double *X, int n, int p, const std::set<int> &component_index, int nthreads)
{
  if (n == 0 || p == 0)
    throw std::invalid_argument("Feature vector is empty.");
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && p != this->feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  if (component_index != std::set<int>{0} && component_index != std::set<int>{-1} &&
      component_index.size() != (size_t)p)
    throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");

  unsigned int threads = nthreads > 0 ? (unsigned int)nthreads : (unsigned int)std::max(1, this->nthreads);
  threads = std::min(threads, std::max(1u, std::thread::hardware_concurrency()));

  const size_t vs = value_size;
  std::vector<double> out((size_t)n * vs, 0.0);
  double *res = out.data(); // column-major: res[k * n + row]

  const bool all_components = component_index == std::set<int>{0};
  const bool intercept_only = component_index == std::set<int>{-1};

  // After a capped purify() the purified components hold only the low orders: the whole forest is predicted from the
  // leaves below, which purification leaves unchanged
  const bool from_leaves = !purified || (purify_capped && all_components);

  if (!from_leaves && rect_forest_) {
    // Components of purify(mode = 2): rectangles (rects.hpp). Same sums in the same order as predict_single_rect.
    const rpf_rect::PurifiedForest &pf = *rect_forest_;
    std::vector<const rpf_rect::Component *> comps;
    std::vector<int> col_feature; // 0-based feature of each column of X
    if (all_components) {
      for (const auto &c : pf.comps)
        comps.push_back(&c.second);
      for (int d = 0; d < feature_size; ++d)
        col_feature.push_back(d);
    } else if (!intercept_only) {
      std::vector<int> S(component_index.begin(), component_index.end());
      auto it = pf.comps.find(S);
      if (it != pf.comps.end())
        comps.push_back(&it->second);
      for (int d : S)
        col_feature.push_back(d - 1);
    }
    const bool with_intercept = all_components || intercept_only;
    const double *x = X;
    parallel_rows(n, threads, [&](int begin, int end) {
      const size_t chunk = 2048, stride = (size_t)std::max(1, feature_size);
      std::vector<int32_t> cells(chunk * stride, 0);
      std::vector<double> acc(chunk * vs);
      for (int r0 = begin; r0 < end; r0 += (int)chunk) {
        const size_t m = (size_t)std::min(end - r0, (int)chunk);
        for (size_t q = 0; q < col_feature.size(); ++q) {
          const int f = col_feature[q];
          for (size_t i = 0; i < m; ++i)
            cells[i * stride + f] = pf.cell(f, x[q * (size_t)n + r0 + i]);
        }
        std::fill(acc.begin(), acc.end(), 0.0);
        if (with_intercept)
          for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < vs; ++k)
              acc[i * vs + k] = pf.intercept[k];
        for (const auto *comp : comps)
          pf.add_component_rows(*comp, cells.data(), stride, m, acc.data(), (int)vs);
        for (size_t i = 0; i < m; ++i)
          for (size_t k = 0; k < vs; ++k)
            res[k * (size_t)n + r0 + i] = acc[i * vs + k] / n_trees;
      }
    });
    return out;
  }

  if (!from_leaves) {
    // Purified trees are grids: per tree dimension, binary-search the cell
    // boundaries. Matching trees are resolved once here, not per row.
    struct GridTree {
      utils::Matrix<std::vector<double>> *values;
      std::vector<const std::vector<double> *> bounds; // cell limits per tree dimension
      std::vector<int> cols;                           // column of X per tree dimension
    };
    std::vector<GridTree> grids;
    for (auto &tree_family : this->tree_families) {
      for (auto &tree : tree_family) {
        const bool is_intercept = tree.first == std::set<int>{0};
        if (intercept_only ? !is_intercept : (!all_components && tree.first != component_index))
          continue;
        GridTree g{&tree.second->GridLeaves.values, {}, {}};
        if (!is_intercept) {
          int pos = 0;
          for (int dim : tree.first) {
            g.bounds.push_back(&tree.second->GridLeaves.lim_list[(size_t)(dim - 1)]);
            g.cols.push_back(all_components ? dim - 1 : pos);
            ++pos;
          }
        }
        grids.push_back(std::move(g));
      }
    }

    const double *x = X;
    parallel_rows(n, threads, [&](int begin, int end) {
      std::vector<int> idx;
      for (const GridTree &g : grids) {
        const size_t nd = g.cols.size();
        idx.assign(std::max<size_t>(1, nd), 0);
        for (int r = begin; r < end; ++r) {
          for (size_t j = 0; j < nd; ++j) {
            const std::vector<double> &b = *g.bounds[j];
            if (b.size() < 2) {
              idx[j] = 0;
              continue;
            }
            int c = (int)(std::upper_bound(b.begin(), b.end(), x[(size_t)g.cols[j] * n + r]) - b.begin());
            idx[j] = std::min(std::max(0, c - 1), (int)b.size() - 2);
          }
          const std::vector<double> &vals = (*g.values)[idx];
          for (size_t k = 0; k < vs && k < vals.size(); ++k)
            res[k * n + r] += vals[k];
        }
      }
      for (size_t k = 0; k < vs; ++k)
        for (int r = begin; r < end; ++r)
          res[k * n + r] /= n_trees;
    });
    return out;
  }

  // Same leaf membership rule as predict_single: an interval edge at the
  // training bound is open, so new data outside the training range still lands in a leaf.
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<FlatTree> trees;
  for (auto &tree_family : this->tree_families) {
    for (auto &tree : tree_family) {
      if (!all_components && tree.first != component_index)
        continue;
      FlatTree ft;
      std::vector<int> dims;
      int pos = 0;
      for (int dim : tree.first) {
        dims.push_back(std::max(0, dim - 1));
        ft.cols.push_back(all_components ? std::max(0, dim - 1) : pos);
        ++pos;
      }
      const auto &leaves = tree.second->leaves;
      ft.n_leaves = leaves.size();
      ft.lo.reserve(ft.n_leaves * dims.size());
      ft.hi.reserve(ft.n_leaves * dims.size());
      ft.values.assign(ft.n_leaves * vs, 0.0);
      for (size_t l = 0; l < ft.n_leaves; ++l) {
        for (int d : dims) {
          const Interval &iv = leaves[l].intervals[(size_t)d];
          ft.lo.push_back(iv.first == lower_bounds[(size_t)d] ? -inf : iv.first);
          ft.hi.push_back(iv.second == upper_bounds[(size_t)d] ? inf : iv.second);
        }
        for (size_t k = 0; k < vs && k < leaves[l].value.size(); ++k)
          ft.values[l * vs + k] = leaves[l].value[k];
      }
      trees.push_back(std::move(ft));
    }
  }

  // From kIndexRows rows on, the leaves of each tree of more than kIndexLeaves leaves are indexed (trees spread over
  // the threads). A row's matching leaves are then added in leaf order, as the scan adds them, so the predictions are
  // bitwise the same.
  if (n >= kIndexRows)
    parallel_rows((int)trees.size(), threads, [&](int begin, int end) {
      LeafIndexer indexer;
      for (int t = begin; t < end; ++t)
        indexer.index(trees[(size_t)t]);
    });

  const double *x = X; // column-major: x[c * n + row]
  parallel_rows(n, threads, [&](int begin, int end) {
    // Tree-outer loop keeps one tree's leaves hot in cache across all rows of the chunk.
    std::vector<double> xr;
    for (const FlatTree &ft : trees) {
      const size_t nd = ft.cols.size();
      xr.resize(nd);
      for (int r = begin; r < end; ++r) {
        for (size_t j = 0; j < nd; ++j)
          xr[j] = x[(size_t)ft.cols[j] * n + r];
        if (!ft.nodes.empty()) {
          // the bucket that holds every leaf of this row, then the scan's test on its leaves, in leaf order
          const FlatTree::Node *node = &ft.nodes[0];
          while (node->dim >= 0)
            node = &ft.nodes[xr[(size_t)node->dim] < node->split ? node->a : node->b];
          for (uint32_t q = node->a; q < node->b; ++q) {
            const size_t l = ft.items[q];
            if (in_leaf(&ft.lo[l * nd], &ft.hi[l * nd], xr.data(), nd))
              for (size_t k = 0; k < vs; ++k)
                res[k * n + r] += ft.values[l * vs + k];
          }
          continue;
        }
        for (size_t l = 0; l < ft.n_leaves; ++l)
          if (in_leaf(&ft.lo[l * nd], &ft.hi[l * nd], xr.data(), nd))
            for (size_t k = 0; k < vs; ++k)
              res[k * n + r] += ft.values[l * vs + k];
      }
    }
    for (size_t k = 0; k < vs; ++k)
      for (int r = begin; r < end; ++r)
        res[k * n + r] /= n_trees;
  });
  return out;
}

Matrix2D
RandomPlantedForest::predict_vector(const std::vector<double> &feature_vec, const std::set<int> &component_index)
{
  Matrix2D predictions;
  if (feature_vec.empty()) {
    warn("Feature vector is empty.");
    return predictions;
  }
  if (component_index == std::set<int>{0} && this->feature_size >= 0 &&
      feature_vec.size() != (size_t)this->feature_size) {
    warn("Feature vector has wrong dimension.");
    return predictions;
  }
  if (component_index == std::set<int>{0}) {
    predictions.push_back(predict_single(feature_vec, component_index));
  } else {
    for (auto vec : feature_vec)
      predictions.push_back(predict_single(std::vector<double>{vec}, component_index));
  }
  return predictions;
}

std::vector<std::vector<int>> RandomPlantedForest::purified_components() const
{
  if (!purified)
    throw std::invalid_argument("The forest is not purified.");
  std::set<std::vector<int>> sets;
  if (rect_forest_) {
    for (const auto &c : rect_forest_->comps)
      sets.insert(c.first);
  } else {
    // mode 1: the trees of the families hold the purified grids (including the subsets of the grown trees)
    for (const auto &family : tree_families)
      for (const auto &tree : family) {
        std::vector<int> S;
        for (int d : tree.first)
          if (d > 0)
            S.push_back(d);
        if (!S.empty())
          sets.insert(S);
      }
  }
  return std::vector<std::vector<int>>(sets.begin(), sets.end());
}

std::vector<double> RandomPlantedForest::predict_components(
    const double *X, int n, int p, const std::vector<std::vector<int>> &sets, int nthreads
)
{
  if (p != feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  const size_t vs = value_size, K = sets.size();
  std::vector<std::vector<int>> S(sets);
  for (auto &s : S) {
    std::sort(s.begin(), s.end());
    if (s.empty() || s.front() < 1 || s.back() > feature_size || std::adjacent_find(s.begin(), s.end()) != s.end())
      throw std::invalid_argument("A component must be a set of distinct features in 1..p.");
  }
  std::vector<double> out((size_t)n * K * vs, 0.0); // column-major: out[(k * vs + c) * n + row]
  if (n == 0 || K == 0)
    return out;

  if (!(purified && rect_forest_)) {
    // leaf scan (unpurified) or mode-1 grids; a purified forest's components without a tree are zero
    std::set<std::set<int>> trees;
    for (const auto &family : tree_families)
      for (const auto &tree : family)
        trees.insert(tree.first);
    std::vector<double> sub;
    for (size_t k = 0; k < K; ++k) {
      const std::set<int> comp(S[k].begin(), S[k].end());
      if (purified && !trees.count(comp))
        continue;
      sub.resize((size_t)n * S[k].size());
      for (size_t q = 0; q < S[k].size(); ++q)
        std::copy(X + (size_t)(S[k][q] - 1) * n, X + (size_t)S[k][q] * n, sub.begin() + q * (size_t)n);
      std::vector<double> r = predict_matrix(sub.data(), n, (int)S[k].size(), comp, nthreads);
      std::copy(r.begin(), r.end(), out.begin() + k * vs * (size_t)n);
    }
    return out;
  }

  // Rectangle store: the same sums in the same order as predict_matrix
  const rpf_rect::PurifiedForest &pf = *rect_forest_;
  std::vector<const rpf_rect::Component *> comps(K, nullptr);
  for (size_t k = 0; k < K; ++k) {
    auto it = pf.comps.find(S[k]);
    if (it != pf.comps.end())
      comps[k] = &it->second;
  }
  unsigned int threads = nthreads > 0 ? (unsigned int)nthreads : (unsigned int)std::max(1, this->nthreads);
  rects_components(pf, comps, X, n, (int)vs, n_trees, threads, out.data());
  return out;
}

void rects_components(
    const rpf_rect::PurifiedForest &pf,
    const std::vector<const rpf_rect::Component *> &comps,
    const double *X,
    int n,
    int vs,
    int n_trees,
    unsigned threads,
    double *out
)
{
  const size_t K = comps.size();
  const int p = (int)pf.cuts.size();
  threads = std::min(threads, std::max(1u, std::thread::hardware_concurrency()));
  // per chunk of rows the cells of all features are computed once, then the components are the outer loop
  parallel_rows(n, threads, [&](int begin, int end) {
    const size_t chunk = 2048, stride = (size_t)std::max(1, p);
    std::vector<int32_t> cells(chunk * stride, 0);
    std::vector<double> acc(chunk * vs);
    for (int r0 = begin; r0 < end; r0 += (int)chunk) {
      const size_t m = (size_t)std::min(end - r0, (int)chunk);
      for (int f = 0; f < p; ++f)
        for (size_t i = 0; i < m; ++i)
          cells[i * stride + f] = pf.cell(f, X[(size_t)f * n + r0 + i]);
      for (size_t k = 0; k < K; ++k) {
        if (!comps[k])
          continue;
        std::fill(acc.begin(), acc.end(), 0.0);
        pf.add_component_rows(*comps[k], cells.data(), stride, m, acc.data(), vs);
        for (int c = 0; c < vs; ++c)
          for (size_t i = 0; i < m; ++i)
            out[(k * vs + c) * (size_t)n + r0 + i] = acc[i * vs + c] / n_trees;
      }
    }
  });
}
