#include "rpf.hpp"
#include "cpf.hpp"

// ----------------- Rcpp include  -----------------


RCPP_MODULE(mod_rpf)
{

  class_<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericMatrix, const NumericMatrix, const NumericVector>()
      .constructor<const NumericVector>()
      .method("set_data", &RandomPlantedForest::set_data)
      .method("set_shape", &RandomPlantedForest::set_shape)
      .method("set_training_data", &RandomPlantedForest::set_training_data)
      .method("get_data", &RandomPlantedForest::get_data)
      .method("get_bounds", &RandomPlantedForest::get_bounds)
      .method("get_shape", &RandomPlantedForest::get_shape)
      .method("set_model", &RandomPlantedForest::set_model)
      .method("get_parameters", &RandomPlantedForest::get_parameters)
      .method("cross_validation", &RandomPlantedForest::cross_validation)
      .method("predict_matrix", &RandomPlantedForest::predict_matrix)
      .method("predict_vector", &RandomPlantedForest::predict_vector)
      .method("MSE", &RandomPlantedForest::MSE)
      .method("purify_threads", static_cast<void (RandomPlantedForest::*)(int,int,int)>(&RandomPlantedForest::purify))
      .method("print", &RandomPlantedForest::print)
      .method("set_parameters", &RandomPlantedForest::set_parameters)
      .method("get_model", &RandomPlantedForest::get_model)
      .method("get_grid_leaves", &RandomPlantedForest::get_grid_leaves)
      .method("set_grid_leaves", &RandomPlantedForest::set_grid_leaves)
      .method("is_purified", &RandomPlantedForest::is_purified);

  class_<ClassificationRPF>("ClassificationRPF")
      .derives<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericMatrix, const NumericMatrix, const String, const NumericVector>()
      .constructor<const String, const NumericVector>()
      .method("set_parameters", &ClassificationRPF::set_parameters);
}
