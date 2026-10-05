#include "rcpp_interface.h"

// Methods shared by the regression and classification classes
template <class T>
void expose_methods(Rcpp::class_<T> &cls)
{
  cls.method("set_data", &T::set_data)
      .method("fit", &T::fit)
      .method("set_shape", &T::set_shape)
      .method("set_training_data", &T::set_training_data)
      .method("get_data", &T::get_data)
      .method("get_bounds", &T::get_bounds)
      .method("get_shape", &T::get_shape)
      .method("set_model", &T::set_model)
      .method("get_parameters", &T::get_parameters)
      .method("cross_validation", &T::cross_validation)
      .method("predict_matrix", &T::predict_matrix)
      .method("predict_vector", &T::predict_vector)
      .method("MSE", &T::MSE)
      .method("purify_threads", &T::purify)
      .method("print", &T::print)
      .method("set_parameters", &T::set_parameters)
      .method("get_model", &T::get_model)
      .method("get_grid_leaves", &T::get_grid_leaves)
      .method("set_grid_leaves", &T::set_grid_leaves)
      .method("is_purified", &T::is_purified);
}

RCPP_MODULE(mod_rpf)
{
  Rcpp::class_<RcppRPF> rpf("RandomPlantedForest");
  rpf.constructor<const Rcpp::NumericMatrix, const Rcpp::NumericMatrix, const Rcpp::NumericVector>()
      .constructor<const Rcpp::NumericVector>();
  expose_methods(rpf);

  Rcpp::class_<RcppCPF> cpf("ClassificationRPF");
  cpf.constructor<const Rcpp::NumericMatrix, const Rcpp::NumericMatrix, std::string, const Rcpp::NumericVector>()
      .constructor<std::string, const Rcpp::NumericVector>();
  expose_methods(cpf);
}
