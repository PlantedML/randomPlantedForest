#ifndef RCPP_INTERFACE_H
#define RCPP_INTERFACE_H

#include "cpf.hpp"
#include "rpf.hpp"
#include <Rcpp.h>

using namespace Rcpp;

class RcppRPF : public RandomPlantedForest {

public:
  RcppRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
          const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0});
  RcppRPF() {};
  void set_data(const NumericMatrix &samples_Y, const NumericMatrix &samples_X);
  NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
  NumericMatrix predict_vector(const NumericVector &X, const NumericVector components = {0});
  void cross_validation(int n_sets = 4, IntegerVector splits = {5, 50}, NumericVector t_tries = {0.2, 0.5, 0.7, 0.9}, IntegerVector split_tries = {1, 2, 5, 10});
  double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true);
  void set_parameters(StringVector keys, NumericVector values);
  void purify_3();
  void print();
  void get_parameters();
  List get_model();
  bool is_purified();

protected:
  double MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true);
};

#endif // RCPP_INTERFACE_H
