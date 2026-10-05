// R binding for the core forest classes: converts Rcpp types at the boundary,
// routes core warnings to the R console and draws tree seeds from R's RNG.
// Everything under include/ and lib/ is free of Rcpp; a binding for another
// language would mirror this file.
#ifndef RCPP_INTERFACE_H
#define RCPP_INTERFACE_H

#include <Rcpp.h>
#include "cpf.hpp"

Matrix2D to_matrix2d(const Rcpp::NumericMatrix &m);
Rcpp::NumericMatrix from_matrix2d(const Matrix2D &m);
std::set<int> to_int_set(const Rcpp::NumericVector &v);

// Positional parameter vectors as built by rpf_param_vector() in R/rpf.R
RPFParams parse_rpf_params(const std::vector<double> &pars);
RandomPlantedForest make_rpf(const Rcpp::NumericVector &parameters);
ClassificationRPF make_cpf(const std::string &loss, const Rcpp::NumericVector &parameters, bool fitting);

Rcpp::List model_to_list(const std::vector<FamilyExport> &model, int feature_size);
std::vector<FamilyExport> list_to_model(const Rcpp::List &model, int feature_size);
Rcpp::List grid_to_list(const std::vector<GridFamilyExport> &grid, int value_size);
std::vector<GridFamilyExport> list_to_grid(const Rcpp::List &grid);

// Route core output to the R console and seeds to R's RNG.
void connect_to_r(RandomPlantedForest &core);

template <class Core>
class RcppForest
{
public:
  // Construct and fit (regression)
  RcppForest(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X,
             const Rcpp::NumericVector parameters)
      : core(make_rpf(parameters))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  // Construct and fit (classification)
  RcppForest(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X,
             const std::string loss, const Rcpp::NumericVector parameters)
      : core(make_cpf(loss, parameters, true))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  // Params-only constructors: no data, no fit. Used by rpf_unmarshal().
  explicit RcppForest(const Rcpp::NumericVector parameters) : core(make_rpf(parameters)) { connect_to_r(core); }
  RcppForest(const std::string loss, const Rcpp::NumericVector parameters)
      : core(make_cpf(loss, parameters, false)) { connect_to_r(core); }

  void set_data(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X)
  {
    core.set_data(to_matrix2d(samples_Y), to_matrix2d(samples_X));
  }
  void fit()
  {
    Rcpp::RNGScope scope;
    core.fit();
  }
  void set_shape(int feature_size, int value_size, int sample_size,
                 const Rcpp::NumericVector lower, const Rcpp::NumericVector upper)
  {
    core.set_shape(feature_size, value_size, sample_size,
                   Rcpp::as<std::vector<double>>(lower), Rcpp::as<std::vector<double>>(upper));
  }
  void set_training_data(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X)
  {
    core.set_training_data(to_matrix2d(samples_Y), to_matrix2d(samples_X));
  }
  Rcpp::List get_data()
  {
    return Rcpp::List::create(Rcpp::Named("X") = from_matrix2d(core.get_X()),
                              Rcpp::Named("Y") = from_matrix2d(core.get_Y()));
  }
  Rcpp::List get_bounds()
  {
    return Rcpp::List::create(Rcpp::Named("lower") = Rcpp::wrap(core.get_lower_bounds()),
                              Rcpp::Named("upper") = Rcpp::wrap(core.get_upper_bounds()));
  }
  Rcpp::List get_shape()
  {
    return Rcpp::List::create(Rcpp::Named("feature_size") = core.get_feature_size(),
                              Rcpp::Named("value_size") = core.get_value_size(),
                              Rcpp::Named("sample_size") = core.get_sample_size());
  }
  Rcpp::List get_model() { return model_to_list(core.get_model(), core.get_feature_size()); }
  void set_model(Rcpp::List &model) { core.set_model(list_to_model(model, core.get_feature_size())); }
  Rcpp::List get_grid_leaves() { return grid_to_list(core.get_grid_leaves(), core.get_value_size()); }
  void set_grid_leaves(Rcpp::List &grid) { core.set_grid_leaves(list_to_grid(grid)); }

  Rcpp::NumericMatrix predict_matrix(const Rcpp::NumericMatrix &X, const Rcpp::NumericVector components, int nthreads)
  {
    std::vector<double> pred = core.predict_matrix(REAL(X), X.nrow(), X.ncol(), to_int_set(components), nthreads);
    Rcpp::NumericMatrix out(X.nrow(), core.get_value_size());
    std::copy(pred.begin(), pred.end(), out.begin());
    return out;
  }
  Rcpp::NumericMatrix predict_vector(const Rcpp::NumericVector &X, const Rcpp::NumericVector components)
  {
    return from_matrix2d(core.predict_vector(Rcpp::as<std::vector<double>>(X), to_int_set(components)));
  }
  double MSE(const Rcpp::NumericMatrix &Y_predicted, const Rcpp::NumericMatrix &Y_true)
  {
    return core.MSE(to_matrix2d(Y_predicted), to_matrix2d(Y_true));
  }
  void purify(int maxp_interaction, int nthreads, int mode)
  {
    if (core.get_X().empty())
      Rcpp::stop("Cannot purify: no training data available. If this forest was "
                 "restored with rpf_unmarshal(), marshal it with include_data = TRUE.");
    core.purify(maxp_interaction, nthreads, mode);
  }
  void cross_validation(int n_sets, Rcpp::IntegerVector splits, Rcpp::NumericVector t_tries, Rcpp::IntegerVector split_tries)
  {
    core.cross_validation(n_sets, Rcpp::as<std::vector<int>>(splits), Rcpp::as<std::vector<double>>(t_tries),
                          Rcpp::as<std::vector<int>>(split_tries));
  }
  void print() { core.print(Rcpp::Rcout); }
  void get_parameters() { core.get_parameters(Rcpp::Rcout); }
  // Refits with the updated parameters.
  void set_parameters(Rcpp::StringVector keys, Rcpp::NumericVector values)
  {
    Rcpp::RNGScope scope;
    core.set_parameters(Rcpp::as<std::vector<std::string>>(keys), Rcpp::as<std::vector<double>>(values));
  }
  bool is_purified() { return core.is_purified(); }

private:
  Core core;

  void fit_data(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X,
                const Rcpp::NumericVector &parameters)
  {
    set_data(samples_Y, samples_X);
    fit();
    if (parameters.size() > 8 && parameters[8] != 0)
      core.cross_validation();
  }
};

typedef RcppForest<RandomPlantedForest> RcppRPF;
typedef RcppForest<ClassificationRPF> RcppCPF;

#endif // RCPP_INTERFACE_H
