// Prediction entry points split out from rpf.cpp for readability and reuse.
#include "rpf.hpp"
#include <algorithm>
#include <iterator>

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

// predict multiple feature vectors
Rcpp::NumericMatrix RandomPlantedForest::predict_matrix(const NumericMatrix &X, const NumericVector components)
{
  std::vector<std::vector<double>> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions;
  if (feature_vec.empty())
    throw std::invalid_argument("Feature vector is empty.");
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec[0].size() != (size_t)this->feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  if (component_index != std::set<int>{0} && component_index != std::set<int>{-1} && component_index.size() != feature_vec[0].size())
    throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");
  for (auto &vec : feature_vec)
  {
    predictions.push_back(predict_single(vec, component_index));
  }
  return from_std_vec(predictions);
}

Rcpp::NumericMatrix RandomPlantedForest::predict_vector(const NumericVector &X, const NumericVector components)
{
  std::vector<double> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions; Rcpp::NumericMatrix res;
  if (feature_vec.empty()) { Rcout << "Feature vector is empty." << std::endl; return res; }
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec.size() != (size_t)this->feature_size) { Rcout << "Feature vector has wrong dimension." << std::endl; return res; }
  if (component_index == std::set<int>{0}) { predictions.push_back(predict_single(feature_vec, component_index)); }
  else { for (auto vec : feature_vec) predictions.push_back(predict_single(std::vector<double>{vec}, component_index)); }
  res = from_std_vec(predictions); return res;
}

