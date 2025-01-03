#include "rpf.hpp"
#include "cpf.hpp"
#include "rcpp_interface.hpp"

// ----------------- Rcpp include  -----------------

RCPP_MODULE(mod_rpf)
{

    class_<RcppRPF>("RandomPlantedForest")
        .constructor<const NumericMatrix, const NumericMatrix, const NumericVector>()
        .method("cross_validation", &RcppRPF::cross_validation)
        .method("predict_matrix", &RcppRPF::predict_matrix)
        .method("predict_vector", &RcppRPF::predict_vector)
        .method("MSE", &RcppRPF::MSE)
        .method("purify", &RcppRPF::purify_3)
        .method("print", &RcppRPF::print)
        .method("get_parameters", &RcppRPF::get_parameters)
        .method("set_parameters", &RcppRPF::set_parameters)
        .method("get_model", &RcppRPF::get_model)
        .method("is_purified", &RcppRPF::is_purified);

    class_<RcppCPF>("ClassificationRPF")
        .constructor<const NumericMatrix, const NumericMatrix, const String, const NumericVector>()
        .method("set_parameters", &RcppCPF::set_parameters)
        .method("predict_matrix", &RcppCPF::predict_matrix)
        .method("predict_vector", &RcppCPF::predict_vector)
        .method("cross_validation", &RcppCPF::cross_validation)
        .method("MSE", &RcppCPF::MSE)
        .method("purify", &RcppCPF::purify_3)
        .method("print", &RcppCPF::print)
        .method("get_parameters", &RcppCPF::get_parameters)
        .method("set_parameters", &RcppCPF::set_parameters)
        .method("get_model", &RcppCPF::get_model)
        .method("is_purified", &RcppCPF::is_purified);
}
