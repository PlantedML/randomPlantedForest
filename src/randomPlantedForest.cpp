#include "rcpp_interface.h"

// Methods shared by the regression and classification classes
template <class T> void expose_methods(Rcpp::class_<T> &cls)
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
      .method("predict_components_matrix", &T::predict_components_matrix)
      .method("purified_components", &T::purified_components)
      .method("MSE", &T::MSE)
      .method("purify_threads", &T::purify)
      .method("print", &T::print)
      .method("set_parameters", &T::set_parameters)
      .method("get_model", &T::get_model)
      .method("get_grid_leaves", &T::get_grid_leaves)
      .method("set_grid_leaves", &T::set_grid_leaves)
      .method("get_grid_family", &T::get_grid_family)
      .method("get_purified_rects", &T::get_purified_rects)
      .method("set_purified_rects", &T::set_purified_rects)
      .method("purify_stats", &T::purify_stats)
      .method("is_purified", &T::is_purified)
      .method("is_purify_capped", &T::is_purify_capped)
      .method("set_purify_capped", &T::set_purify_capped)
      .method("get_purify_cap", &T::get_purify_cap)
      .method("set_purify_cap", &T::set_purify_cap)
      .method("inbag", &T::inbag)
      .method("set_inbag", &T::set_inbag)
      .method("purify_oob", &T::purify_oob)
      .method("purify_columns", &T::purify_columns)
      .method("store_components", &T::store_components)
      // the out-of-bag path and the adaptive fit
      .method("set_oob_path", &T::set_oob_path)
      .method("oob_path", &T::oob_path)
      .method("set_adaptive", &T::set_adaptive)
      .method("adaptive_stats", &T::adaptive_stats)
      .method("oob_row_errors", &T::oob_row_errors)
      // the panel form of purified components
      .method("coarsen_select", &T::coarsen_select)
      .method("coarsen_oob", &T::coarsen_oob)
      .method("coarsen_eval", &T::coarsen_eval);
}

RCPP_MODULE(mod_rpf)
{
  // internal, for tests: inclusion counts of the two draws of split_structure = "leaves"
  Rcpp::function("rpf_leaves_draw_counts", &leaves_draw_counts);
  // internal, for tests: the logit and exponential losses' criterion and steps on one split (cpf.hpp)
  Rcpp::function("rpf_classif_check", &classif_check);
  // a coarsened component from its stored form (prune_components(max_panels = ))
  Rcpp::function("rpf_coarse_build", &rpf_coarse_build);
  Rcpp::function("rpf_coarse_eval", &rpf_coarse_eval);
  Rcpp::function("rpf_coarse_valid", &rpf_coarse_valid);
  Rcpp::function("rpf_coarse_predict", &rpf_coarse_predict);
  // the stored pieces of a pruned model's kept components (prune_components())
  Rcpp::function("rpf_store_build", &rpf_store_build);
  Rcpp::function("rpf_store_eval", &rpf_store_eval);
  Rcpp::function("rpf_store_valid", &rpf_store_valid);

  Rcpp::class_<RcppRPF> rpf("RandomPlantedForest");
  rpf.constructor<const Rcpp::NumericMatrix, const Rcpp::NumericMatrix, const Rcpp::NumericVector>()
      .constructor<const Rcpp::NumericMatrix, const Rcpp::NumericMatrix, const Rcpp::NumericVector, const Rcpp::List>()
      .constructor<const Rcpp::NumericVector>()
      .constructor<const Rcpp::NumericVector, const Rcpp::List>();
  expose_methods(rpf);

  Rcpp::class_<RcppCPF> cpf("ClassificationRPF");
  cpf.constructor<const Rcpp::NumericMatrix, const Rcpp::NumericMatrix, std::string, const Rcpp::NumericVector>()
      .constructor<
          const Rcpp::NumericMatrix,
          const Rcpp::NumericMatrix,
          std::string,
          const Rcpp::NumericVector,
          const Rcpp::List>()
      .constructor<std::string, const Rcpp::NumericVector>()
      .constructor<std::string, const Rcpp::NumericVector, const Rcpp::List>();
  expose_methods(cpf);
  // internal, for tests: the scores and loss weights of each family's sample rows (set_parameters("keep_weights", 1))
  cpf.method("kept_weights", &RcppCPF::kept_weights);
}
