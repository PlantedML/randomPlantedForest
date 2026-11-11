#include "rcpp_interface.h"
#include "internal_utils.hpp"

using namespace Rcpp;

Matrix2D to_matrix2d(const NumericMatrix &m)
{
  // Copies column-major R memory into row-major nested vectors
  const int rows = m.nrow(), cols = m.ncol();
  Matrix2D out((size_t)rows, std::vector<double>((size_t)cols));
  const double *data = REAL(m);
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
      out[(size_t)i][(size_t)j] = data[(size_t)j * rows + i];
  return out;
}

NumericMatrix from_matrix2d(const Matrix2D &m)
{
  if (m.empty())
    return NumericMatrix();
  NumericMatrix out(m.size(), m[0].size());
  for (size_t i = 0; i < m.size(); ++i)
    for (size_t j = 0; j < m[i].size(); ++j)
      out(i, j) = m[i][j];
  return out;
}

std::set<int> to_int_set(const NumericVector &v)
{
  return std::set<int>(v.begin(), v.end());
}

RPFParams parse_rpf_params(const std::vector<double> &pars)
{
  if (pars.size() < 12 || pars.size() > 15)
    Rcpp::stop("RandomPlantedForest requires 12 to 15 parameters, got %d", (int)pars.size());
  RPFParams p;
  p.max_interaction = pars[0];
  p.n_trees = pars[1];
  p.n_splits = pars[2];
  p.split_try = pars[3];
  p.t_try = pars[4];
  p.purify_forest = pars[5];
  p.deterministic = pars[6];
  p.nthreads = pars[7];
  p.cross_validate = pars[8];
  p.split_decay_rate = pars[9];
  p.max_candidates = static_cast<size_t>(pars[10]);
  p.delete_leaves = (pars[11] != 0);
  if (pars.size() >= 13)
    p.split_structure_mode = static_cast<int>(pars[12]);
  if (pars.size() >= 14) // older vectors (marshaled models) lack these two: 1 and the exact draw
    p.min_leaf_rows = static_cast<int>(pars[13]);
  if (pars.size() >= 15)
    p.exact_draw = (pars[14] != 0);
  return p;
}

std::vector<std::set<int>> parse_constraints(const List &constraints)
{
  std::vector<std::set<int>> out;
  for (R_xlen_t i = 0; i < constraints.size(); ++i) {
    std::vector<int> dims = as<std::vector<int>>(constraints[i]);
    for (int d : dims)
      if (d < 1)
        Rcpp::stop("interaction constraints: feature indices must be >= 1, got %d", d);
    out.emplace_back(dims.begin(), dims.end());
  }
  return out;
}

RandomPlantedForest make_rpf(const NumericVector &parameters, const List &constraints)
{
  RPFParams params = parse_rpf_params(as<std::vector<double>>(parameters));
  params.interaction_constraints = parse_constraints(constraints);
  return RandomPlantedForest(params);
}

// Classification vectors are regression's 15 values (the 13, min_leaf_rows,
// exact_draw), then exact_step and prior_rows (the pseudo-rows a of logit and
// exponential, cpf.hpp); a vector without them gets TRUE and 1.
ClassificationRPF
make_cpf(const std::string &loss, const NumericVector &parameters, bool fitting, const List &constraints)
{
  std::vector<double> pars = as<std::vector<double>>(parameters);
  const bool sized = pars.size() >= 15 && pars.size() <= 17;
  RPFParams params;
  bool exact_step = true;
  double prior_rows = 1;
  if (sized) {
    params = parse_rpf_params(std::vector<double>(pars.begin(), pars.begin() + 15));
    if (pars.size() >= 16)
      exact_step = pars[15] != 0;
    if (pars.size() >= 17)
      prior_rows = pars[16];
  } else {
    if (!fitting)
      Rcpp::stop("ClassificationRPF requires 15 to 17 parameters, got %d", (int)pars.size());
    Rcout << "Wrong number of parameters - set to default." << std::endl;
    if (pars.size() >= 13)
      params.split_structure_mode = static_cast<int>(pars[12]);
  }
  if (!(prior_rows > 0) || !std::isfinite(prior_rows))
    Rcpp::stop("ClassificationRPF: the prior rows must be positive and finite, got %g", prior_rows);
  params.interaction_constraints = parse_constraints(constraints);
  return ClassificationRPF(loss, params, exact_step, prior_rows);
}

List model_to_list(const std::vector<FamilyExport> &model, int feature_size)
{
  List out;
  for (const FamilyExport &family : model) {
    List variables, family_values, family_intervals;
    for (const TreeExport &tree : family) {
      variables.push_back(IntegerVector(tree.variables.begin(), tree.variables.end()));
      List tree_values, tree_intervals;
      for (size_t k = 0; k < tree.values.size(); ++k) {
        tree_values.push_back(NumericVector(tree.values[k].begin(), tree.values[k].end()));
        NumericMatrix leaf_intervals(2, feature_size);
        for (int l = 0; l < feature_size; ++l) {
          leaf_intervals(0, l) = tree.intervals[k][l].first;
          leaf_intervals(1, l) = tree.intervals[k][l].second;
        }
        tree_intervals.push_back(leaf_intervals);
      }
      family_values.push_back(tree_values);
      family_intervals.push_back(tree_intervals);
    }
    out.push_back(
        List::create(
            Named("variables") = variables, Named("values") = family_values, Named("intervals") = family_intervals
        )
    );
  }
  return out;
}

std::vector<FamilyExport> list_to_model(const List &model, int feature_size)
{
  std::vector<FamilyExport> out;
  for (int i = 0; i < model.size(); ++i) {
    List family = model[i];
    List variables = family["variables"];
    List values = family["values"];
    List intervals = family["intervals"];
    FamilyExport fam;
    for (int j = 0; j < variables.size(); ++j) {
      TreeExport tree;
      IntegerVector tree_variables = variables[j];
      tree.variables = std::set<int>(tree_variables.begin(), tree_variables.end());
      List tree_values = values[j];
      List tree_intervals = intervals[j];
      for (int k = 0; k < tree_values.size(); ++k) {
        tree.values.push_back(as<std::vector<double>>(tree_values[k]));
        NumericMatrix leaf_intervals = tree_intervals[k];
        if (leaf_intervals.nrow() < 2 || leaf_intervals.ncol() < feature_size)
          Rcpp::stop(
              "Corrupt model data: leaf interval matrix has dimensions %dx%d, expected 2x%d.",
              leaf_intervals.nrow(),
              leaf_intervals.ncol(),
              feature_size
          );
        std::vector<Interval> ivs(feature_size);
        for (int l = 0; l < feature_size; ++l)
          ivs[l] = Interval{leaf_intervals(0, l), leaf_intervals(1, l)};
        tree.intervals.push_back(ivs);
      }
      fam.push_back(std::move(tree));
    }
    out.push_back(std::move(fam));
  }
  return out;
}

List grid_to_list(const std::vector<GridFamilyExport> &grid, int value_size)
{
  List families;
  for (const GridFamilyExport &family : grid) {
    List lim_list;
    for (const auto &v : family.lim_list)
      lim_list.push_back(wrap(v));
    List trees;
    for (const GridTreeExport &tree : family.trees) {
      NumericMatrix values(tree.values.size(), value_size);
      for (size_t e = 0; e < tree.values.size(); ++e)
        for (int p = 0; p < value_size; ++p)
          values(e, p) = tree.values[e][p];
      trees.push_back(
          List::create(
              Named("variables") = IntegerVector(tree.variables.begin(), tree.variables.end()),
              Named("dims") = wrap(tree.dims),
              Named("values") = values
          )
      );
    }
    families.push_back(List::create(Named("lim_list") = lim_list, Named("trees") = trees));
  }
  return families;
}

std::vector<GridFamilyExport> list_to_grid(const List &grid)
{
  std::vector<GridFamilyExport> out;
  for (int i = 0; i < grid.size(); ++i) {
    List family = grid[i];
    GridFamilyExport fam;
    List lim_list = family["lim_list"];
    for (int l = 0; l < lim_list.size(); ++l)
      fam.lim_list.push_back(as<std::vector<double>>(lim_list[l]));
    List trees = family["trees"];
    for (int j = 0; j < trees.size(); ++j) {
      List tr = trees[j];
      GridTreeExport tree;
      IntegerVector variables = tr["variables"];
      tree.variables = std::set<int>(variables.begin(), variables.end());
      tree.dims = as<std::vector<int>>(tr["dims"]);
      tree.values = to_matrix2d(tr["values"]);
      fam.trees.push_back(std::move(tree));
    }
    out.push_back(std::move(fam));
  }
  return out;
}

namespace {
// Two 32-bit draws from R's RNG composed into one 64-bit seed
std::uint64_t seed_from_r()
{
  std::uint64_t hi = static_cast<std::uint64_t>(R::runif(0.0, 4294967296.0));
  std::uint64_t lo = static_cast<std::uint64_t>(R::runif(0.0, 4294967296.0));
  return (hi << 32) ^ lo;
}
} // namespace

List rects_to_list(const rpf_rect::PurifiedForest *pf)
{
  if (pf == nullptr)
    return List();
  List cuts;
  for (const auto &c : pf->cuts)
    cuts.push_back(wrap(c));
  List comps;
  for (const auto &kv : pf->comps) {
    const rpf_rect::Component &c = kv.second;
    List groups;
    for (const rpf_rect::Group &g : c.groups) {
      std::vector<int> nodes;
      nodes.reserve(g.index.nodes.size() * 7);
      for (const auto &nd : g.index.nodes) {
        int f[7] = {nd.dim, nd.split, nd.left, nd.right, nd.cst, nd.r0, nd.r1};
        nodes.insert(nodes.end(), f, f + 7);
      }
      groups.push_back(
          List::create(
              Named("feat") = wrap(g.feat),
              Named("table") = wrap(g.table),
              Named("box") = IntegerVector(g.rects.box.begin(), g.rects.box.end()),
              Named("w") = wrap(g.rects.w),
              Named("nodes") = wrap(nodes),
              Named("csts") = wrap(g.index.csts),
              Named("ids") = wrap(g.index.ids)
          )
      );
    }
    comps.push_back(
        List::create(Named("S") = wrap(c.S), Named("constant") = wrap(c.constant), Named("groups") = groups)
    );
  }
  return List::create(
      Named("cuts") = cuts,
      Named("intercept") = wrap(pf->intercept),
      Named("n_rects") = (double)pf->n_rects,
      Named("maxp") = pf->maxp,
      Named("components") = comps
  );
}

rpf_rect::PurifiedForest list_to_rects(const List &x, int value_size)
{
  rpf_rect::PurifiedForest pf;
  List cuts = x["cuts"];
  for (int i = 0; i < cuts.size(); ++i)
    pf.cuts.push_back(as<std::vector<double>>(cuts[i]));
  pf.intercept = as<std::vector<double>>(x["intercept"]);
  pf.n_rects = (size_t)as<double>(x["n_rects"]);
  pf.maxp = as<int>(x["maxp"]);
  List comps = x["components"];
  for (int i = 0; i < comps.size(); ++i) {
    List cl = comps[i];
    rpf_rect::Component c;
    c.S = as<std::vector<int>>(cl["S"]);
    c.constant = as<std::vector<double>>(cl["constant"]);
    List groups = cl["groups"];
    for (int j = 0; j < groups.size(); ++j) {
      List gl = groups[j];
      rpf_rect::Group g;
      g.feat = as<std::vector<int>>(gl["feat"]);
      for (int f : g.feat) {
        auto it = std::find(c.S.begin(), c.S.end(), f + 1);
        g.pos.push_back(it == c.S.end() ? -1 : (int)(it - c.S.begin()));
      }
      g.table = as<std::vector<double>>(gl["table"]);
      g.rects = rpf_rect::RectList((int)g.feat.size(), value_size);
      IntegerVector box = gl["box"];
      g.rects.box.assign(box.begin(), box.end());
      g.rects.w = as<std::vector<double>>(gl["w"]);
      std::vector<int> nodes = as<std::vector<int>>(gl["nodes"]);
      if (nodes.size() % 7 != 0)
        throw std::invalid_argument("Corrupt purified components: index nodes.");
      g.index.nodes.reserve(nodes.size() / 7);
      for (size_t k = 0; k < nodes.size(); k += 7) {
        rpf_rect::KdIndex::Node nd;
        nd.dim = nodes[k];
        nd.split = nodes[k + 1];
        nd.left = nodes[k + 2];
        nd.right = nodes[k + 3];
        nd.cst = nodes[k + 4];
        nd.r0 = nodes[k + 5];
        nd.r1 = nodes[k + 6];
        g.index.nodes.push_back(nd);
      }
      g.index.csts = as<std::vector<double>>(gl["csts"]);
      g.index.ids = as<std::vector<int>>(gl["ids"]);
      c.groups.push_back(std::move(g));
    }
    std::vector<int> key = c.S;
    pf.comps.emplace(key, std::move(c));
  }
  return pf;
}

// An empty list without a rectangle store
List purify_stats_to_list(const PurifyStats &st)
{
  if (!st.rectangles)
    return List();
  return List::create(
      Named("n_components") = st.n_components,
      Named("n_rects") = st.n_rects,
      Named("n_groups") = st.n_groups,
      Named("n_index_nodes") = st.n_index_nodes,
      Named("n_table_cells") = st.n_table_cells,
      Named("bytes") = st.bytes
  );
}

// checkpoints, the n x (checkpoints * value_size) sums sum0 / sum1 (column-major, as R's matrices) and the n x 2 count
List oob_path_to_list(const OobPath &o)
{
  const int n = (int)(o.count.size() / 2);
  const int ncol = n > 0 ? (int)(o.sum0.size() / (size_t)n) : 0;
  NumericMatrix s0(n, ncol), s1(n, ncol);
  std::copy(o.sum0.begin(), o.sum0.end(), s0.begin());
  std::copy(o.sum1.begin(), o.sum1.end(), s1.begin());
  IntegerMatrix count(n, 2);
  std::copy(o.count.begin(), o.count.end(), count.begin());
  return List::create(
      Named("checkpoints") = IntegerVector(o.checkpoints.begin(), o.checkpoints.end()),
      Named("sum0") = s0,
      Named("sum1") = s1,
      Named("count") = count
  );
}

List adaptive_stats_to_list(const AdaptiveStats &st)
{
  return List::create(
      Named("pilot_bound") = st.pilot_bound,
      Named("growths") = st.growths,
      Named("regrowths") = st.regrowths,
      Named("splits_grown") = st.splits_grown,
      Named("pauses") = st.pauses,
      Named("resumes") = st.resumes,
      Named("max_paused_gb") = st.max_paused_gb,
      Named("seconds_grow") = st.seconds_grow,
      Named("seconds_stopper") = st.seconds_stopper,
      Named("seconds_rebuild") = st.seconds_rebuild,
      Named("bound_trace") = IntegerVector(st.bound_trace.begin(), st.bound_trace.end())
  );
}

// ---- The pruning without the purified store (purify_columns(), store_components())

List columns_to_list(ComponentColumns &cc, int n, int value_size)
{
  const int K = (int)cc.oob.sets.size();
  List sets(K);
  for (int k = 0; k < K; ++k)
    sets[k] = IntegerVector(cc.oob.sets[k].begin(), cc.oob.sets[k].end());
  // each vector is freed once copied
  auto num = [n](std::vector<double> &v, int ncol) {
    NumericMatrix m(n, ncol);
    std::copy(v.begin(), v.end(), m.begin());
    std::vector<double>().swap(v);
    return m;
  };
  List out = List::create(
      Named("sets") = sets,
      Named("own") = num(cc.own, K * value_size),
      Named("intercept") = NumericVector(cc.intercept.begin(), cc.intercept.end())
  );
  if (!cc.oob.count.empty()) {
    IntegerMatrix count(n, 2);
    std::copy(cc.oob.count.begin(), cc.oob.count.end(), count.begin());
    out["sum0"] = num(cc.oob.sum0, K * value_size);
    out["sum1"] = num(cc.oob.sum1, K * value_size);
    out["intercept0"] = num(cc.oob.intercept0, value_size);
    out["intercept1"] = num(cc.oob.intercept1, value_size);
    out["count"] = count;
  }
  return out;
}

List store_to_list(ComponentStore &st)
{
  const int vs = st.value_size;
  List comps;
  for (auto it = st.pieces.comps.begin(); it != st.pieces.comps.end(); it = st.pieces.comps.erase(it)) {
    List groups;
    for (const auto &g : it->second.second) {
      const rpf_rect::RectList &rl = g.second;
      const int R = (int)rl.size(), m = rl.m;
      IntegerMatrix box(R, 2 * m);
      NumericMatrix w(R, vs);
      for (int r = 0; r < R; ++r) {
        for (int q = 0; q < 2 * m; ++q)
          box(r, q) = rl.box[(size_t)r * 2 * m + q];
        for (int c = 0; c < vs; ++c)
          w(r, c) = rl.w[(size_t)r * vs + c];
      }
      groups.push_back(List::create(Named("jmask") = g.first, Named("box") = box, Named("w") = w));
    }
    comps.push_back(
        List::create(
            Named("S") = IntegerVector(it->first.begin(), it->first.end()),
            Named("constant") = wrap(it->second.first),
            Named("groups") = groups
        )
    );
  }
  List cuts;
  for (const auto &c : st.cuts)
    cuts.push_back(wrap(c));
  return List::create(
      Named("p") = st.p,
      Named("value_size") = vs,
      Named("n_trees") = st.n_trees,
      Named("cuts") = cuts,
      Named("comps") = comps
  );
}

ComponentStore list_to_store(const List &x)
{
  ComponentStore st;
  st.p = as<int>(x["p"]);
  st.value_size = as<int>(x["value_size"]);
  st.n_trees = as<int>(x["n_trees"]);
  const int vs = st.value_size;
  List cuts = x["cuts"];
  for (R_xlen_t d = 0; d < cuts.size(); ++d)
    st.cuts.push_back(as<std::vector<double>>(cuts[d]));
  List comps = x["comps"];
  for (R_xlen_t k = 0; k < comps.size(); ++k) {
    List c = comps[k];
    auto &entry = st.pieces.comps[as<std::vector<int>>(c["S"])];
    entry.first = as<std::vector<double>>(c["constant"]);
    List groups = c["groups"];
    for (R_xlen_t g = 0; g < groups.size(); ++g) {
      List gr = groups[g];
      IntegerMatrix box = gr["box"];
      NumericMatrix w = gr["w"];
      if (box.ncol() % 2 != 0 || box.nrow() != w.nrow() || w.ncol() != vs)
        stop("Corrupt component store.");
      const int R = box.nrow(), m = box.ncol() / 2;
      rpf_rect::RectList rl(m, vs);
      rl.box.resize((size_t)R * 2 * m);
      rl.w.resize((size_t)R * vs);
      for (int r = 0; r < R; ++r) {
        for (int q = 0; q < 2 * m; ++q)
          rl.box[(size_t)r * 2 * m + q] = box(r, q);
        for (int v = 0; v < vs; ++v)
          rl.w[(size_t)r * vs + v] = w(r, v);
      }
      if (!entry.second.emplace(as<int>(gr["jmask"]), std::move(rl)).second)
        stop("Corrupt component store.");
    }
  }
  return st;
}

SEXP rpf_store_build(const List store, int nthreads)
{
  ComponentStore st = list_to_store(store);
  return XPtr<BuiltStore>(new BuiltStore(st, (unsigned)std::max(1, nthreads)), true);
}

NumericMatrix rpf_store_eval(SEXP store, const NumericMatrix &X, const List sets, int nthreads)
{
  if (!rpf_store_valid(store))
    stop("The built component store is gone (e.g. after saveRDS()).");
  XPtr<BuiltStore> B(store);
  std::vector<std::vector<int>> S;
  for (R_xlen_t k = 0; k < sets.size(); ++k)
    S.push_back(as<std::vector<int>>(sets[k]));
  std::vector<double> v = B->eval(REAL(X), X.nrow(), X.ncol(), S, nthreads);
  NumericMatrix out(X.nrow(), (int)S.size() * B->value_size());
  std::copy(v.begin(), v.end(), out.begin());
  return out;
}

bool rpf_store_valid(SEXP store)
{
  return TYPEOF(store) == EXTPTRSXP && R_ExternalPtrAddr(store) != nullptr;
}

void connect_to_r(RandomPlantedForest &core)
{
  core.verbose_out = &Rcout;
  core.seed_source = seed_from_r;
  core.check_interrupt = [] { Rcpp::checkUserInterrupt(); };
}

NumericVector leaves_draw_counts(
    const NumericVector &ages, double decay, int n_candidates, int reps, bool exact, const IntegerVector &rows_ok
)
{
  if (n_candidates < 1 || reps < 0)
    Rcpp::stop("n_candidates must be positive and reps non-negative");
  if (rows_ok.size() != 0 && rows_ok.size() != ages.size())
    Rcpp::stop("rows_ok must be empty or as long as ages");
  Rcpp::RNGScope scope;
  const std::vector<double> counts = rpf_utils::leaves_draw_counts(
      as<std::vector<double>>(ages),
      decay,
      (size_t)n_candidates,
      (size_t)reps,
      exact,
      as<std::vector<int>>(rows_ok),
      seed_from_r()
  );
  return NumericVector(counts.begin(), counts.end());
}

List classif_check(
    const std::string &loss,
    const NumericMatrix &Y,
    const NumericMatrix &F,
    const IntegerVector &side,
    bool exact_step,
    double prior_rows
)
{
  const size_t n = (size_t)Y.nrow(), K = (size_t)Y.ncol();
  if ((size_t)F.nrow() != n || (size_t)F.ncol() != K || (size_t)side.size() != n || K == 0)
    Rcpp::stop("classif_check: Y, F and side must have the same rows, Y and F the same columns");
  if (loss != "logit" && loss != "exponential")
    Rcpp::stop("classif_check: loss must be \"logit\" or \"exponential\"");
  const DriftFree df(
      loss == "logit" ? DriftFree::logit : (K == 1 ? DriftFree::exp_binary : DriftFree::samme), K, prior_rows
  );
  const Matrix2D y = to_matrix2d(Y), f = to_matrix2d(F);
  std::vector<double> ts(df.T, 0.0), tb(df.T, 0.0), tp(df.T, 0.0), p(K);
  std::vector<int> rows_s, rows_b;
  for (size_t i = 0; i < n; ++i) {
    double *out = side[i] == 0 ? ts.data() : tb.data();
    (side[i] == 0 ? rows_s : rows_b).push_back((int)i);
    if (df.kind == DriftFree::logit) {
      double mx = 0;
      for (size_t q = 0; q < K; ++q)
        mx = std::max(mx, f[i][q]);
      double den = std::exp(-mx);
      for (size_t q = 0; q < K; ++q) {
        p[q] = std::exp(f[i][q] - mx);
        den += p[q];
      }
      for (size_t q = 0; q < K; ++q)
        p[q] /= den;
      df.add_logit(y[i].data(), f[i].data(), p.data(), out);
    } else {
      // binary exp(-y F / 2); SAMME exp(-F_c), c the class coded 1 (the last if none is)
      double w;
      if (K == 1) {
        w = std::exp(-0.5 * y[i][0] * f[i][0]);
      } else {
        size_t c = 0;
        while (c + 1 < K && !(y[i][c] > 0))
          ++c;
        w = std::exp(-f[i][c]);
      }
      df.add_exp(y[i].data(), f[i].data(), w, out);
    }
  }
  for (size_t q = 0; q < df.T; ++q)
    tp[q] = ts[q] + tb[q];
  const double ns = (double)rows_s.size(), nb = (double)rows_b.size();
  NumericVector newton_s(K), newton_b(K);
  double scs = 0, scb = 0;
  DriftFree::Scratch sc;
  const double vs = df.value(ts.data(), ns, newton_s.begin(), &scs, sc);
  const double vb = df.value(tb.data(), nb, newton_b.begin(), &scb, sc);
  const double vp = df.value(tp.data(), ns + nb, nullptr, nullptr, sc);
  NumericVector step_s = Rcpp::clone(newton_s), step_b = Rcpp::clone(newton_b);
  if (df.kind == DriftFree::logit && exact_step) {
    std::vector<double> c;
    df.logit_exact_step(rows_s, f, y, c);
    step_s = NumericVector(c.begin(), c.end());
    df.logit_exact_step(rows_b, f, y, c);
    step_b = NumericVector(c.begin(), c.end());
  }
  return List::create(
      Named("score") = df.score(ts.data(), tb.data(), ns, nb, vp, sc),
      Named("raw") = vs + vb - vp,
      Named("step_s") = step_s,
      Named("step_b") = step_b,
      Named("newton_s") = newton_s,
      Named("newton_b") = newton_b
  );
}

// ---- The panel form of purified components (coarsen.hpp)

rpf_coarse::Spec list_to_spec(const List &sp)
{
  rpf_coarse::Spec s;
  s.S = as<std::vector<int>>(sp["S"]);
  s.C = as<std::vector<int>>(sp["C"]);
  List cuts = sp["cuts"];
  for (R_xlen_t i = 0; i < cuts.size(); ++i)
    s.cuts.push_back(as<std::vector<double>>(cuts[i]));
  return s;
}

std::vector<rpf_coarse::Spec> list_to_specs(const List &specs)
{
  std::vector<rpf_coarse::Spec> out;
  for (R_xlen_t k = 0; k < specs.size(); ++k)
    out.push_back(list_to_spec(specs[k]));
  return out;
}

static List cuts_to_list(const std::vector<std::vector<double>> &cuts)
{
  List out(cuts.size());
  for (size_t i = 0; i < cuts.size(); ++i)
    out[(R_xlen_t)i] = NumericVector(cuts[i].begin(), cuts[i].end());
  return out;
}

static std::vector<std::vector<double>> list_to_cuts(const List &x)
{
  std::vector<std::vector<double>> out;
  for (R_xlen_t i = 0; i < x.size(); ++i)
    out.push_back(as<std::vector<double>>(x[i]));
  return out;
}

// The stored form as plain R data: S, cfeat, ccuts, dfeat, dcuts (features 1-based), box (R x 2m: lo, then hi, in
// forest cells of dcuts), w (R x cells), n_trees
List form_to_list(const rpf_coarse::Form &f)
{
  const rpf_rect::RectList &rl = f.rects;
  const int m = rl.m, nc = rl.vs;
  const size_t R = rl.size();
  IntegerMatrix box((int)R, 2 * m);
  NumericMatrix w((int)R, nc);
  for (size_t r = 0; r < R; ++r) {
    for (int q = 0; q < m; ++q) {
      box((int)r, q) = rl.lo(r)[q];
      box((int)r, m + q) = rl.hi(r)[q];
    }
    for (int c = 0; c < nc; ++c)
      w((int)r, c) = rl.w[r * nc + c];
  }
  return List::create(
      Named("S") = IntegerVector(f.S.begin(), f.S.end()),
      Named("cfeat") = IntegerVector(f.cfeat.begin(), f.cfeat.end()),
      Named("ccuts") = cuts_to_list(f.ccuts),
      Named("dfeat") = IntegerVector(f.dfeat.begin(), f.dfeat.end()),
      Named("dcuts") = cuts_to_list(f.dcuts),
      Named("box") = box,
      Named("w") = w,
      Named("n_trees") = f.n_trees
  );
}

rpf_coarse::Form list_to_form(const List &x)
{
  rpf_coarse::Form f;
  f.S = as<std::vector<int>>(x["S"]);
  f.cfeat = as<std::vector<int>>(x["cfeat"]);
  f.dfeat = as<std::vector<int>>(x["dfeat"]);
  f.ccuts = list_to_cuts(x["ccuts"]);
  f.dcuts = list_to_cuts(x["dcuts"]);
  f.n_trees = as<int>(x["n_trees"]);
  IntegerMatrix box = x["box"];
  NumericMatrix w = x["w"];
  const int m = (int)f.dfeat.size(), nc = w.ncol();
  const size_t R = (size_t)w.nrow();
  if (box.ncol() != 2 * m || (size_t)box.nrow() != R)
    stop("Corrupt coarsened component.");
  f.rects = rpf_rect::RectList(m, nc);
  std::vector<int32_t> l(std::max(1, m)), h(std::max(1, m));
  std::vector<double> wr(nc);
  for (size_t r = 0; r < R; ++r) {
    for (int q = 0; q < m; ++q) {
      l[q] = box((int)r, q);
      h[q] = box((int)r, m + q);
    }
    for (int c = 0; c < nc; ++c)
      wr[c] = w((int)r, c);
    f.rects.push(l.data(), h.data(), wr.data(), 1.0);
  }
  return f;
}

// One component's panel form: S, form ("absent", "asis", "coarsened"), forest_panels, var_forest, seconds; if
// coarsened also C, q, cuts, cells, var_kept, g (the training-row values), candidates (a data.frame: combo (1-based),
// q, cells, var_kept), stored (form_to_list()) and n_eval (candidates evaluated by difference array, by k-d index)
List selected_to_list(const rpf_coarse::Selected &R)
{
  List o = List::create(
      Named("S") = IntegerVector(R.S.begin(), R.S.end()),
      Named("form") = R.form == 0 ? "absent" : (R.form == 1 ? "asis" : "coarsened"),
      Named("forest_panels") = R.forest_panels,
      Named("var_forest") = R.var_forest,
      Named("seconds") = R.seconds
  );
  if (R.form == 1) {
    // kept as it is: its panels, the features of S with the fewest own cut points and these cut points
    IntegerVector pfeat;
    for (int a : R.cpos)
      pfeat.push_back(R.S[a]);
    o["C"] = pfeat;
    o["cuts"] = cuts_to_list(R.cuts);
  }
  if (R.form != 2)
    return o;
  IntegerVector cfeat, q(R.q.begin(), R.q.end());
  for (size_t a = 0; a < R.cpos.size(); ++a)
    cfeat.push_back(R.S[R.cpos[a]]);
  const size_t nc = R.cands.size();
  IntegerVector cc(nc), ccells(nc);
  NumericVector cvk(nc);
  CharacterVector cq(nc);
  for (size_t a = 0; a < nc; ++a) {
    cc[a] = R.cands[a].combo + 1;
    ccells[a] = R.cands[a].cells;
    cvk[a] = R.cands[a].var_kept;
    std::string s;
    for (size_t b = 0; b < R.cands[a].q.size(); ++b)
      s += (b ? "," : "") + std::to_string(R.cands[a].q[b]);
    cq[a] = s;
  }
  int cells = 1;
  for (const auto &t : R.cuts)
    cells *= (int)t.size() + 1;
  o["C"] = cfeat;
  o["q"] = q;
  o["cuts"] = cuts_to_list(R.cuts);
  o["cells"] = cells;
  o["var_kept"] = R.var_kept;
  // one score column: a vector and a form; V > 1: an n x V matrix and a list of V forms
  const int V = (int)R.stored.size();
  if (V == 1)
    o["g"] = NumericVector(R.g.begin(), R.g.end());
  else {
    NumericMatrix g((int)(R.g.size() / (size_t)V), V);
    std::copy(R.g.begin(), R.g.end(), g.begin());
    o["g"] = g;
  }
  o["candidates"] = DataFrame::create(
      Named("combo") = cc,
      Named("q") = cq,
      Named("cells") = ccells,
      Named("var_kept") = cvk,
      Named("stringsAsFactors") = false
  );
  if (V == 1)
    o["stored"] = form_to_list(R.stored[0]);
  else {
    List forms(V);
    for (int v = 0; v < V; ++v)
      forms[v] = form_to_list(R.stored[(size_t)v]);
    o["stored"] = forms;
  }
  o["n_eval"] = IntegerVector::create(R.n_diff, R.n_kd);
  return o;
}

// A stored form built once for point queries: an external pointer (NULL after serialization)
SEXP rpf_coarse_build(const List stored)
{
  return XPtr<rpf_coarse::Model>(new rpf_coarse::Model(list_to_form(stored)), true);
}

NumericVector rpf_coarse_eval(SEXP model, const NumericMatrix &X)
{
  XPtr<rpf_coarse::Model> M(model);
  if (M.get() == nullptr)
    stop("The built coarsened component is gone (e.g. after saveRDS()).");
  std::vector<double> v = M->eval(REAL(X), X.nrow(), X.ncol());
  return NumericVector(v.begin(), v.end());
}

bool rpf_coarse_valid(SEXP model)
{
  return TYPEOF(model) == EXTPTRSXP && R_ExternalPtrAddr(model) != nullptr;
}

NumericVector rpf_coarse_predict(const List stored, const NumericMatrix &X)
{
  rpf_coarse::Model M(list_to_form(stored));
  std::vector<double> v = M.eval(REAL(X), X.nrow(), X.ncol());
  return NumericVector(v.begin(), v.end());
}
