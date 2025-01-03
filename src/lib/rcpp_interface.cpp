#include "rcpp_interface.hpp"

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

RcppRPF::RcppRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                 const NumericVector parameters)
    : RandomPlantedForest(toStd2D(samples_Y), toStd2D(samples_X),
                          toStd1D(parameters))
{
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
  return RandomPlantedForest::get_model();
}

bool RcppRPF::is_purified()
{
  return RandomPlantedForest::is_purified();
}

RcppCPF::RcppCPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                 const std::string loss, const NumericVector parameters)
    : ClassificationRPF(toStd2D(samples_Y), toStd2D(samples_X), loss,
                        toStd1D(parameters))
{
}

void RcppCPF::set_parameters(StringVector keys, NumericVector values)
{
  ClassificationRPF::set_parameters(std::vector<std::string>(keys.begin(), keys.end()),
                                    std::vector<double>(values.begin(), values.end()));
}
