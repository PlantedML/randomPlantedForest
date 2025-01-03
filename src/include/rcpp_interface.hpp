#ifndef RCPP_INTERFACE_H
#define RCPP_INTERFACE_H

#include "cpf.hpp"
#include "rpf.hpp"
#include <Rcpp.h>

using namespace Rcpp;

class RcppInterface
{
public:
  virtual void set_parameters(StringVector keys, NumericVector values) = 0;
  virtual NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components) = 0;
  virtual NumericMatrix predict_vector(const NumericVector &X, const NumericVector components) = 0;
  virtual void cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries) = 0;
  virtual double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true) = 0;
};

class RcppRPF : public RandomPlantedForest, public RcppInterface
{

public:
  RcppRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
          const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0});
  RcppRPF() {};
  NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components = {0}) override;
  NumericMatrix predict_vector(const NumericVector &X, const NumericVector components = {0}) override;
  void cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries) override;
  double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true) override;
  void set_parameters(StringVector keys, NumericVector values) override;
  void purify_3();
  void print();
  void get_parameters();
  List get_model();
  bool is_purified();

protected:
  double MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true);
};

class RcppCPF : public ClassificationRPF, public RcppInterface
{
public:
  RcppCPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
          const std::string loss = "L2", const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0});
  void set_parameters(StringVector keys, NumericVector values) override;
  NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components = {0}) override;
  NumericMatrix predict_vector(const NumericVector &X, const NumericVector components = {0}) override;
  void cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries) override;
  double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true) override;

  void purify_3();
  void print();
  void get_parameters();
  List get_model();
  bool is_purified();
};

#endif // RCPP_INTERFACE_H
