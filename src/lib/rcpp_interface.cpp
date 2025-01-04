#include "rcpp_interface.hpp"
#include "random_utils.hpp"

// Helper to convert NumericMatrix to std::vector<std::vector<double>>
static std::vector<std::vector<double>> toStd2D(const Rcpp::NumericMatrix &mat)
{
  std::vector<std::vector<double>> vec(mat.nrow(), std::vector<double>(mat.ncol()));
  for (int i = 0; i < mat.nrow(); ++i)
  {
    for (int j = 0; j < mat.ncol(); ++j)
    {
      vec[i][j] = mat(i, j);
    }
  }
  return vec;
}

// Helper to convert NumericVector to std::vector<double>
static std::vector<double> toStd1D(const Rcpp::NumericVector &vec)
{
  return std::vector<double>(vec.begin(), vec.end());
}

// New helper function to convert tree families to R List
static List treeFamiliesToRList(const std::vector<TreeFamily> &families, int feature_size)
{
  List model;
  for (const auto &family : families)
  {
    List variables, family_values, family_intervals;
    for (const auto &tree : family)
    {
      List tree_values;
      List tree_intervals;
      variables.push_back(from_std_set(tree.first));
      for (const auto &leaf : tree.second->get_leaves())
      {
        NumericMatrix leaf_values;
        for (const auto &val : leaf.value)
        {
          leaf_values.push_back(val);
        }
        tree_values.push_back(leaf_values);

        NumericVector intervals;
        for (const auto &interval : leaf.intervals)
        {
          intervals.push_back(interval.first);
          intervals.push_back(interval.second);
        }
        NumericMatrix leaf_intervals(2, feature_size, intervals.begin());
        tree_intervals.push_back(leaf_intervals);
      }
      family_intervals.push_back(tree_intervals);
      family_values.push_back(tree_values);
    }
    model.push_back(List::create(Named("variables") = variables, _["values"] = family_values, _["intervals"] = family_intervals));
  }
  return model;
}

RcppRPF::RcppRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                 const NumericVector parameters)
    : RandomPlantedForest(toStd2D(samples_Y), toStd2D(samples_X),
                          toStd1D(parameters))
{
  utils::RandomGenerator::use_r_random();

  this->fit();

  if (cross_validate)
  {
    RandomPlantedForest::cross_validation();
  }
}

NumericMatrix RcppRPF::predict_matrix(const NumericMatrix &X, const NumericVector components)
{
  auto result = RandomPlantedForest::predict_matrix(toStd2D(X), toStd1D(components));
  NumericMatrix rResult(result.size(), result[0].size());
  for (size_t i = 0; i < result.size(); ++i)
  {
    for (size_t j = 0; j < result[i].size(); ++j)
    {
      rResult(i, j) = result[i][j];
    }
  }
  return rResult;
}

NumericMatrix RcppRPF::predict_vector(const NumericVector &X, const NumericVector components)
{
  auto result = RandomPlantedForest::predict_vector(toStd1D(X), toStd1D(components));
  NumericMatrix rResult(result.size(), result[0].size());
  for (size_t i = 0; i < result.size(); ++i)
  {
    for (size_t j = 0; j < result[i].size(); ++j)
    {
      rResult(i, j) = result[i][j];
    }
  }
  return rResult;
}

void RcppRPF::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries)
{
  RandomPlantedForest::cross_validation(n_sets, std::vector<int>(splits.begin(), splits.end()),
                                        std::vector<double>(t_tries.begin(), t_tries.end()),
                                        std::vector<int>(split_tries.begin(), split_tries.end()));
}

double RcppRPF::MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true)
{
  return RandomPlantedForest::MSE(toStd2D(Y_predicted), toStd2D(Y_true));
}

void RcppRPF::set_parameters(StringVector keys, NumericVector values)
{
  RandomPlantedForest::set_parameters(std::vector<std::string>(keys.begin(), keys.end()),
                                      std::vector<double>(values.begin(), values.end()));
}

double RcppRPF::MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true)
{
  return RandomPlantedForest::MSE_vec(toStd1D(Y_predicted), toStd1D(Y_true));
}

void RcppRPF::purify_3()
{
  RandomPlantedForest::purify_3();
}

void RcppRPF::print()
{
  RandomPlantedForest::print();
}

void RcppRPF::get_parameters()
{
  RandomPlantedForest::get_parameters();
}

List RcppRPF::get_model()
{
  return treeFamiliesToRList(tree_families, feature_size);
}

bool RcppRPF::is_purified()
{
  return RandomPlantedForest::is_purified();
}

RcppCPF::RcppCPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                 const std::string loss, const NumericVector parameters)
    : ClassificationRPF(toStd2D(samples_Y), toStd2D(samples_X), loss, toStd1D(parameters))
{
  utils::RandomGenerator::use_r_random();
  RandomPlantedForest::fit();

  if (cross_validate)
  {
    ClassificationRPF::cross_validation();
  }
}

void RcppCPF::set_parameters(StringVector keys, NumericVector values)
{
  ClassificationRPF::set_parameters(std::vector<std::string>(keys.begin(), keys.end()),
                                    std::vector<double>(values.begin(), values.end()));
}

NumericMatrix RcppCPF::predict_matrix(const NumericMatrix &X, const NumericVector components)
{
  auto result = ClassificationRPF::predict_matrix(toStd2D(X), toStd1D(components));
  NumericMatrix rResult(result.size(), result[0].size());
  for (size_t i = 0; i < result.size(); ++i)
  {
    for (size_t j = 0; j < result[i].size(); ++j)
    {
      rResult(i, j) = result[i][j];
    }
  }
  return rResult;
}

NumericMatrix RcppCPF::predict_vector(const NumericVector &X, const NumericVector components)
{
  auto result = ClassificationRPF::predict_vector(toStd1D(X), toStd1D(components));
  NumericMatrix rResult(result.size(), result[0].size());
  for (size_t i = 0; i < result.size(); ++i)
  {
    for (size_t j = 0; j < result[i].size(); ++j)
    {
      rResult(i, j) = result[i][j];
    }
  }
  return rResult;
}

void RcppCPF::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries)
{
  ClassificationRPF::cross_validation(n_sets, std::vector<int>(splits.begin(), splits.end()),
                                      std::vector<double>(t_tries.begin(), t_tries.end()),
                                      std::vector<int>(split_tries.begin(), split_tries.end()));
}

double RcppCPF::MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true)
{
  return ClassificationRPF::MSE(toStd2D(Y_predicted), toStd2D(Y_true));
}

void RcppCPF::purify_3()
{
  RandomPlantedForest::purify_3();
}

void RcppCPF::print()
{
  RandomPlantedForest::print();
}

void RcppCPF::get_parameters()
{
  RandomPlantedForest::get_parameters();
}

List RcppCPF::get_model()
{
  return treeFamiliesToRList(tree_families, feature_size);
}

bool RcppCPF::is_purified()
{
  return RandomPlantedForest::is_purified();
}
