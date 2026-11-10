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
RandomPlantedForest make_rpf(const Rcpp::NumericVector &parameters, const Rcpp::List &constraints = Rcpp::List());
ClassificationRPF make_cpf(
    const std::string &loss,
    const Rcpp::NumericVector &parameters,
    bool fitting,
    const Rcpp::List &constraints = Rcpp::List()
);

Rcpp::List model_to_list(const std::vector<FamilyExport> &model, int feature_size);
std::vector<FamilyExport> list_to_model(const Rcpp::List &model, int feature_size);
Rcpp::List grid_to_list(const std::vector<GridFamilyExport> &grid, int value_size);
std::vector<GridFamilyExport> list_to_grid(const Rcpp::List &grid);
Rcpp::List rects_to_list(const rpf_rect::PurifiedForest *pf);
rpf_rect::PurifiedForest list_to_rects(const Rcpp::List &x, int value_size);
Rcpp::List purify_stats_to_list(const PurifyStats &st);
Rcpp::List oob_path_to_list(const OobPath &o);
Rcpp::List adaptive_stats_to_list(const AdaptiveStats &st);

// Route core output to the R console and seeds to R's RNG.
void connect_to_r(RandomPlantedForest &core);

// For tests: rpf_utils::leaves_draw_counts() (the exact or the fast draw of
// `leaves`) with a seed from R's RNG. rows_ok: empty, or 0 / 1 per entry (0: no
// valid threshold under min_leaf_rows, never drawn).
Rcpp::NumericVector leaves_draw_counts(
    const Rcpp::NumericVector &ages,
    double decay,
    int n_candidates,
    int reps,
    bool exact,
    const Rcpp::IntegerVector &rows_ok
);

template <class Core> class RcppForest {
public:
  // Construct and fit (regression)
  RcppForest(
      const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X, const Rcpp::NumericVector parameters
  )
      : core(make_rpf(parameters))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  // Construct and fit (classification)
  RcppForest(
      const Rcpp::NumericMatrix &samples_Y,
      const Rcpp::NumericMatrix &samples_X,
      const std::string loss,
      const Rcpp::NumericVector parameters
  )
      : core(make_cpf(loss, parameters, true))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  // Construct and fit with interaction constraints (a list of allowed feature sets, 1-based column indices)
  RcppForest(
      const Rcpp::NumericMatrix &samples_Y,
      const Rcpp::NumericMatrix &samples_X,
      const Rcpp::NumericVector parameters,
      const Rcpp::List constraints
  )
      : core(make_rpf(parameters, constraints))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  RcppForest(
      const Rcpp::NumericMatrix &samples_Y,
      const Rcpp::NumericMatrix &samples_X,
      const std::string loss,
      const Rcpp::NumericVector parameters,
      const Rcpp::List constraints
  )
      : core(make_cpf(loss, parameters, true, constraints))
  {
    connect_to_r(core);
    fit_data(samples_Y, samples_X, parameters);
  }
  // Params-only constructors: no data, no fit. Used by rpf_unmarshal().
  explicit RcppForest(const Rcpp::NumericVector parameters) : core(make_rpf(parameters)) { connect_to_r(core); }
  RcppForest(const std::string loss, const Rcpp::NumericVector parameters) : core(make_cpf(loss, parameters, false))
  {
    connect_to_r(core);
  }
  // Params-only with interaction constraints: configured for set_data() + fit() (the OOB path, the adaptive fit).
  RcppForest(const Rcpp::NumericVector parameters, const Rcpp::List constraints)
      : core(make_rpf(parameters, constraints))
  {
    connect_to_r(core);
  }
  RcppForest(const std::string loss, const Rcpp::NumericVector parameters, const Rcpp::List constraints)
      : core(make_cpf(loss, parameters, false, constraints))
  {
    connect_to_r(core);
  }

  void set_data(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X)
  {
    core.set_data(to_matrix2d(samples_Y), to_matrix2d(samples_X));
  }
  void fit()
  {
    Rcpp::RNGScope scope;
    core.fit();
  }
  void set_shape(
      int feature_size,
      int value_size,
      int sample_size,
      const Rcpp::NumericVector lower,
      const Rcpp::NumericVector upper
  )
  {
    core.set_shape(
        feature_size,
        value_size,
        sample_size,
        Rcpp::as<std::vector<double>>(lower),
        Rcpp::as<std::vector<double>>(upper)
    );
  }
  void set_training_data(const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X)
  {
    core.set_training_data(to_matrix2d(samples_Y), to_matrix2d(samples_X));
  }
  Rcpp::List get_data()
  {
    return Rcpp::List::create(
        Rcpp::Named("X") = from_matrix2d(core.get_X()), Rcpp::Named("Y") = from_matrix2d(core.get_Y())
    );
  }
  Rcpp::List get_bounds()
  {
    return Rcpp::List::create(
        Rcpp::Named("lower") = Rcpp::wrap(core.get_lower_bounds()),
        Rcpp::Named("upper") = Rcpp::wrap(core.get_upper_bounds())
    );
  }
  Rcpp::List get_shape()
  {
    return Rcpp::List::create(
        Rcpp::Named("feature_size") = core.get_feature_size(),
        Rcpp::Named("value_size") = core.get_value_size(),
        Rcpp::Named("sample_size") = core.get_sample_size()
    );
  }
  Rcpp::List get_model() { return model_to_list(core.get_model(), core.get_feature_size()); }
  void set_model(Rcpp::List &model) { core.set_model(list_to_model(model, core.get_feature_size())); }
  Rcpp::List get_grid_leaves() { return grid_to_list(core.get_grid_leaves(), core.get_value_size()); }
  void set_grid_leaves(Rcpp::List &grid) { core.set_grid_leaves(list_to_grid(grid)); }
  // One element of get_grid_leaves(); family is 1-based, as in R
  Rcpp::List get_grid_family(int family)
  {
    if (core.get_purified_rects() != nullptr && core.get_X().empty())
      Rcpp::stop(
          "The boxes of a forest purified with mode = 2 are computed from its training data. If this forest was "
          "restored with rpf_unmarshal(), marshal it with include_data = TRUE."
      );
    Rcpp::List one = grid_to_list({core.get_grid_family(family - 1)}, core.get_value_size());
    return one[0];
  }
  Rcpp::List get_purified_rects() { return rects_to_list(core.get_purified_rects()); }
  void set_purified_rects(Rcpp::List &x) { core.set_purified_rects(list_to_rects(x, core.get_value_size())); }
  Rcpp::List purify_stats() { return purify_stats_to_list(core.purify_stats()); }

  Rcpp::NumericMatrix predict_matrix(const Rcpp::NumericMatrix &X, const Rcpp::NumericVector components, int nthreads)
  {
    std::vector<double> pred = core.predict_matrix(REAL(X), X.nrow(), X.ncol(), to_int_set(components), nthreads);
    Rcpp::NumericMatrix out(X.nrow(), core.get_value_size());
    std::copy(pred.begin(), pred.end(), out.begin());
    return out;
  }
  // sets: a list of 1-based feature indices per component
  Rcpp::NumericMatrix predict_components_matrix(const Rcpp::NumericMatrix &X, const Rcpp::List sets, int nthreads)
  {
    std::vector<std::vector<int>> S;
    for (R_xlen_t k = 0; k < sets.size(); ++k)
      S.push_back(Rcpp::as<std::vector<int>>(sets[k]));
    std::vector<double> pred = core.predict_components(REAL(X), X.nrow(), X.ncol(), S, nthreads);
    Rcpp::NumericMatrix out(X.nrow(), (int)(S.size() * core.get_value_size()));
    std::copy(pred.begin(), pred.end(), out.begin());
    return out;
  }
  Rcpp::List purified_components()
  {
    std::vector<std::vector<int>> sets = core.purified_components();
    Rcpp::List out(sets.size());
    for (size_t i = 0; i < sets.size(); ++i)
      out[(R_xlen_t)i] = Rcpp::IntegerVector(sets[i].begin(), sets[i].end());
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
      Rcpp::stop(
          "Cannot purify: no training data available. If this forest was "
          "restored with rpf_unmarshal(), marshal it with include_data = TRUE."
      );
    core.purify(maxp_interaction, nthreads, mode);
  }
  void
  cross_validation(int n_sets, Rcpp::IntegerVector splits, Rcpp::NumericVector t_tries, Rcpp::IntegerVector split_tries)
  {
    core.cross_validation(
        n_sets,
        Rcpp::as<std::vector<int>>(splits),
        Rcpp::as<std::vector<double>>(t_tries),
        Rcpp::as<std::vector<int>>(split_tries)
    );
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
  bool is_purify_capped() { return core.is_purify_capped(); }
  void set_purify_capped(bool capped) { core.set_purify_capped(capped); }
  int get_purify_cap() { return core.get_purify_cap(); }
  void set_purify_cap(int cap) { core.set_purify_cap(cap); }
  // In-bag counts: an n x ntrees integer matrix (how often each training row is in each family's bootstrap sample)
  Rcpp::IntegerMatrix inbag()
  {
    const auto &inb = core.get_inbag();
    const int n = core.get_sample_size();
    if (inb.empty())
      Rcpp::stop(
          "The forest has no in-bag counts. If it was restored with rpf_unmarshal(), marshal it with "
          "include_data = TRUE."
      );
    Rcpp::IntegerMatrix out(n, (int)inb.size());
    for (size_t f = 0; f < inb.size(); ++f)
      std::copy(inb[f].begin(), inb[f].end(), out.begin() + f * (size_t)n);
    return out;
  }
  void set_inbag(const Rcpp::IntegerMatrix &x)
  {
    std::vector<std::vector<uint8_t>> inb((size_t)x.ncol());
    for (int f = 0; f < x.ncol(); ++f)
      for (int i = 0; i < x.nrow(); ++i) {
        const int c = x(i, f);
        if (c < 0 || c > 255)
          Rcpp::stop("In-bag counts must lie in 0..255.");
        inb[(size_t)f].push_back((uint8_t)c);
      }
    core.set_inbag(std::move(inb));
  }
  // purify(maxp_interaction, nthreads, 2), and the out-of-bag components of the training rows (see OobComponents)
  Rcpp::List purify_oob(int maxp_interaction, int nthreads)
  {
    if (core.get_X().empty())
      Rcpp::stop(
          "Cannot purify: no training data available. If this forest was "
          "restored with rpf_unmarshal(), marshal it with include_data = TRUE."
      );
    if (core.get_inbag().empty())
      Rcpp::stop(
          "The forest has no in-bag counts. If it was restored with rpf_unmarshal(), marshal it with "
          "include_data = TRUE."
      );
    OobComponents o = core.purify_oob(maxp_interaction, nthreads);
    const int n = core.get_sample_size(), vs = core.get_value_size(), K = (int)o.sets.size();
    Rcpp::List sets(K);
    for (int k = 0; k < K; ++k)
      sets[k] = Rcpp::IntegerVector(o.sets[k].begin(), o.sets[k].end());
    auto num = [n](std::vector<double> &v, int ncol) {
      Rcpp::NumericMatrix m(n, ncol);
      std::copy(v.begin(), v.end(), m.begin());
      std::vector<double>().swap(v);
      return m;
    };
    Rcpp::IntegerMatrix count(n, 2);
    std::copy(o.count.begin(), o.count.end(), count.begin());
    return Rcpp::List::create(
        Rcpp::Named("sets") = sets,
        Rcpp::Named("sum0") = num(o.sum0, K * vs),
        Rcpp::Named("sum1") = num(o.sum1, K * vs),
        Rcpp::Named("intercept0") = num(o.intercept0, vs),
        Rcpp::Named("intercept1") = num(o.intercept1, vs),
        Rcpp::Named("count") = count
    );
  }
  // The out-of-bag path (regression): set the checkpoints before fit(), read the sums after it (see OobPath)
  void set_oob_path(const Rcpp::IntegerVector &checkpoints)
  {
    core.set_oob_path(Rcpp::as<std::vector<int>>(checkpoints));
  }
  Rcpp::List oob_path()
  {
    if (!core.has_oob_path())
      Rcpp::stop("The forest has no OOB path: fit it with `oob_path`.");
    return oob_path_to_list(core.oob_path());
  }
  // The adaptive fit (regression): the next fit() grows to `checkpoints` (the last one the forest's splits) and calls
  // stopper(c, path) at each checkpoint c (1-based) with the OOB path there (as oob_path() gives it, one column); it
  // returns 0 to grow on or the index of the chosen checkpoint (see AdaptiveStopper)
  void set_adaptive(
      const Rcpp::IntegerVector &checkpoints, Rcpp::Function stopper, int patience, double floor, double pause_gb
  )
  {
    core.set_adaptive(
        Rcpp::as<std::vector<int>>(checkpoints),
        [stopper](size_t c, const OobPath &path) -> int {
          return Rcpp::as<int>(stopper((int)c + 1, oob_path_to_list(path)));
        },
        patience,
        floor,
        pause_gb
    );
  }
  Rcpp::List adaptive_stats() { return adaptive_stats_to_list(core.adaptive_stats()); }

private:
  Core core;

  void fit_data(
      const Rcpp::NumericMatrix &samples_Y, const Rcpp::NumericMatrix &samples_X, const Rcpp::NumericVector &parameters
  )
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
