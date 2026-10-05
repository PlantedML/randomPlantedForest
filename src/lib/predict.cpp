// Prediction entry points split out from rpf.cpp for readability and reuse.
#include "rpf.hpp"
#include <algorithm>
#include <exception>
#include <iterator>
#include <limits>
#include <thread>

// predict single feature vector
std::vector<double> RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index)
{
  std::vector<double> total_res = std::vector<double>(value_size, 0);

  if (!purified)
  {
    // consider all components
    if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (auto &dim : tree.first)
            {
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              {
                valid = false;
                break;
              }
            }
            if (valid)
            {
              for (size_t p = 0; p < value_size && p < leaf.value.size(); ++p)
              {
                total_res[p] += leaf.value[p];
              }
            }
          }
        }
      }
    }
    else
    { // choose components for prediction
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          // only consider trees with same dimensions as component_index
          if (tree.first != component_index)
            continue;

          std::vector<int> dims;
          for (auto dim : tree.first)
          {
            dims.push_back(dim);
          }

          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (unsigned int i = 0; i < dims.size(); ++i)
            {
              int dim = dims[i];
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[i] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[i] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              {
                valid = false;
                break;
              }
            }
            if (valid)
            {
              for (size_t p = 0; p < value_size && p < leaf.value.size(); ++p)
              {
                total_res[p] += leaf.value[p];
              }
            }
          }
        }
      }
    }
  }
  else
  {
    if (component_index == std::set<int>{-1})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0})
          {
            leaf_index = std::vector<int>(tree.first.size(), 0);
            
            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p)
            {
              total_res[p] += vals[p];
            }
          }
        }
      }
    }
    else if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0})
          {
            leaf_index = std::vector<int>(tree.first.size(), 0);
          }
          else
          {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // convert to 0-based original feature index
              }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2)
              {
                leaf_index[dim_index] = 0;
                continue;
              }
              // Use the original feature index into X, not the position within the tree's dim set
              auto it = std::upper_bound(bounds.begin(), bounds.end(), X[dim]);
              int c = static_cast<int>(std::distance(bounds.begin(), it));
              leaf_index[dim_index] = std::min(std::max(0, c - 1), (int)bounds.size() - 2);
            }
          }
          for (int &index : leaf_index) index = std::max(0, index);
          {
            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p)
            {
              total_res[p] += vals[p];
            }
          }
        }
      }
    }
    else
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          if (tree.first != component_index)
            continue;
          std::vector<int> leaf_index(tree.first.size(), -1);
          if (tree.first == std::set<int>{0})
          {
            leaf_index = std::vector<int>(tree.first.size(), 0);
          }
          else
          {
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // 0-based original feature index for bounds lookup only
              }
              auto &bounds = tree.second->GridLeaves.lim_list[dim];
              if (bounds.size() < 2)
              {
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
          for (int &index : leaf_index) index = std::max(0, index);
          {
            const auto &vals = tree.second->GridLeaves.values[leaf_index];
            for (size_t p = 0; p < value_size && p < vals.size(); ++p)
            {
              total_res[p] += vals[p];
            }
          }
        }
      }
    }
  }

  return total_res / n_trees;
}

namespace
{
// One tree's leaves in contiguous arrays, so the per-row scan streams through memory.
struct FlatTree
{
  std::vector<int> cols;      // column of X compared against each tree dimension
  std::vector<double> lo, hi; // n_leaves x cols.size(), row-major; open bounds stored as -/+inf
  std::vector<double> values; // n_leaves x value_size
  size_t n_leaves = 0;
};

// Runs f(begin, end) on contiguous row ranges, one per thread.
template <class F>
void parallel_rows(int n, unsigned int threads, F f)
{
  threads = std::max(1u, std::min<unsigned int>(threads, (unsigned int)n));
  if (threads == 1)
  {
    f(0, n);
    return;
  }
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(threads);
  int chunk = (n + (int)threads - 1) / (int)threads;
  for (unsigned int t = 0; t < threads; ++t)
  {
    int begin = (int)t * chunk, end = std::min(n, begin + chunk);
    pool.emplace_back([&f, &errors, t, begin, end]
                      {
      // an exception escaping a std::thread would terminate R
      try { f(begin, end); }
      catch (...) { errors[t] = std::current_exception(); } });
  }
  for (auto &th : pool)
    th.join();
  for (auto &e : errors)
    if (e)
      std::rethrow_exception(e);
}
} // namespace

// predict multiple feature vectors
std::vector<double> RandomPlantedForest::predict_matrix(const double *X, int n, int p, const std::set<int> &component_index, int nthreads)
{
  if (n == 0 || p == 0)
    throw std::invalid_argument("Feature vector is empty.");
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && p != this->feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  if (component_index != std::set<int>{0} && component_index != std::set<int>{-1} && component_index.size() != (size_t)p)
    throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");

  unsigned int threads = nthreads > 0 ? (unsigned int)nthreads : (unsigned int)std::max(1, this->nthreads);
  threads = std::min(threads, std::max(1u, std::thread::hardware_concurrency()));

  const size_t vs = value_size;
  std::vector<double> out((size_t)n * vs, 0.0);
  double *res = out.data(); // column-major: res[k * n + row]

  const bool all_components = component_index == std::set<int>{0};
  const bool intercept_only = component_index == std::set<int>{-1};

  if (purified)
  {
    // Purified trees are grids: per tree dimension, binary-search the cell
    // boundaries. Matching trees are resolved once here, not per row.
    struct GridTree
    {
      utils::Matrix<std::vector<double>> *values;
      std::vector<const std::vector<double> *> bounds; // cell limits per tree dimension
      std::vector<int> cols;                           // column of X per tree dimension
    };
    std::vector<GridTree> grids;
    for (auto &tree_family : this->tree_families)
    {
      for (auto &tree : tree_family)
      {
        const bool is_intercept = tree.first == std::set<int>{0};
        if (intercept_only ? !is_intercept : (!all_components && tree.first != component_index))
          continue;
        GridTree g{&tree.second->GridLeaves.values, {}, {}};
        if (!is_intercept)
        {
          int pos = 0;
          for (int dim : tree.first)
          {
            g.bounds.push_back(&tree.second->GridLeaves.lim_list[(size_t)(dim - 1)]);
            g.cols.push_back(all_components ? dim - 1 : pos);
            ++pos;
          }
        }
        grids.push_back(std::move(g));
      }
    }

    const double *x = X;
    parallel_rows(n, threads, [&](int begin, int end)
                  {
      std::vector<int> idx;
      for (const GridTree &g : grids) {
        const size_t nd = g.cols.size();
        idx.assign(std::max<size_t>(1, nd), 0);
        for (int r = begin; r < end; ++r) {
          for (size_t j = 0; j < nd; ++j) {
            const std::vector<double> &b = *g.bounds[j];
            if (b.size() < 2) { idx[j] = 0; continue; }
            int c = (int)(std::upper_bound(b.begin(), b.end(), x[(size_t)g.cols[j] * n + r]) - b.begin());
            idx[j] = std::min(std::max(0, c - 1), (int)b.size() - 2);
          }
          const std::vector<double> &vals = (*g.values)[idx];
          for (size_t k = 0; k < vs && k < vals.size(); ++k) res[k * n + r] += vals[k];
        }
      }
      for (size_t k = 0; k < vs; ++k)
        for (int r = begin; r < end; ++r) res[k * n + r] /= n_trees; });
    return out;
  }

  // Same leaf membership rule as predict_single: an interval edge at the
  // training bound is open, so new data outside the training range still lands in a leaf.
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<FlatTree> trees;
  for (auto &tree_family : this->tree_families)
  {
    for (auto &tree : tree_family)
    {
      if (!all_components && tree.first != component_index)
        continue;
      FlatTree ft;
      std::vector<int> dims;
      int pos = 0;
      for (int dim : tree.first)
      {
        dims.push_back(std::max(0, dim - 1));
        ft.cols.push_back(all_components ? std::max(0, dim - 1) : pos);
        ++pos;
      }
      const auto &leaves = tree.second->leaves;
      ft.n_leaves = leaves.size();
      ft.lo.reserve(ft.n_leaves * dims.size());
      ft.hi.reserve(ft.n_leaves * dims.size());
      ft.values.assign(ft.n_leaves * vs, 0.0);
      for (size_t l = 0; l < ft.n_leaves; ++l)
      {
        for (int d : dims)
        {
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

  const double *x = X; // column-major: x[c * n + row]
  parallel_rows(n, threads, [&](int begin, int end)
                {
    // Tree-outer loop keeps one tree's leaves hot in cache across all rows of the chunk.
    std::vector<double> xr;
    for (const FlatTree &ft : trees) {
      const size_t nd = ft.cols.size();
      xr.resize(nd);
      for (int r = begin; r < end; ++r) {
        for (size_t j = 0; j < nd; ++j) xr[j] = x[(size_t)ft.cols[j] * n + r];
        for (size_t l = 0; l < ft.n_leaves; ++l) {
          const double *lo = &ft.lo[l * nd], *hi = &ft.hi[l * nd];
          bool valid = true;
          for (size_t j = 0; j < nd; ++j) {
            // explicit inf checks keep NaN/inf inputs matching predict_single
            if (!((lo[j] <= xr[j] || lo[j] == -inf) && (xr[j] < hi[j] || hi[j] == inf))) {
              valid = false;
              break;
            }
          }
          if (valid)
            for (size_t k = 0; k < vs; ++k) res[k * n + r] += ft.values[l * vs + k];
        }
      }
    }
    for (size_t k = 0; k < vs; ++k)
      for (int r = begin; r < end; ++r) res[k * n + r] /= n_trees; });
  return out;
}

Matrix2D RandomPlantedForest::predict_vector(const std::vector<double> &feature_vec, const std::set<int> &component_index)
{
  Matrix2D predictions;
  if (feature_vec.empty()) { warn("Feature vector is empty."); return predictions; }
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec.size() != (size_t)this->feature_size) { warn("Feature vector has wrong dimension."); return predictions; }
  if (component_index == std::set<int>{0}) { predictions.push_back(predict_single(feature_vec, component_index)); }
  else { for (auto vec : feature_vec) predictions.push_back(predict_single(std::vector<double>{vec}, component_index)); }
  return predictions;
}

