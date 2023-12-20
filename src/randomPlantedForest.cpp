#include "rpf.hpp"
#include "cpf.hpp"

// ----------------- Rcpp include  -----------------


RCPP_MODULE(mod_rpf)
{

  class_<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericMatrix, const NumericMatrix, const NumericVector>()
      .method("set_data", &RandomPlantedForest::set_data)
      .method("cross_validation", &RandomPlantedForest::cross_validation)
      .method("predict_matrix", &RandomPlantedForest::predict_matrix)
      .method("predict_vector", &RandomPlantedForest::predict_vector)
      .method("MSE", &RandomPlantedForest::MSE)
      .method("purify", &RandomPlantedForest::purify_3)
      .method("print", &RandomPlantedForest::print)
      .method("get_parameters", &RandomPlantedForest::get_parameters)
      .method("set_parameters", &RandomPlantedForest::set_parameters)
      .method("get_model", &RandomPlantedForest::get_model)
      .method("is_purified", &RandomPlantedForest::is_purified);

  class_<ClassificationRPF>("ClassificationRPF")
      .derives<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericMatrix, const NumericMatrix, const String, const NumericVector>()
      .method("set_parameters", &ClassificationRPF::set_parameters);
}
