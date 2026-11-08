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

RandomPlantedForest make_rpf(const NumericVector &parameters)
{
  return RandomPlantedForest(parse_rpf_params(as<std::vector<double>>(parameters)));
}

// Classification vectors append delta and epsilon to the 13 regression values.
ClassificationRPF make_cpf(const std::string &loss, const NumericVector &parameters, bool fitting)
{
  std::vector<double> pars = as<std::vector<double>>(parameters);
  std::vector<double> base(pars.begin(), pars.begin() + std::min<size_t>(pars.size(), pars.size() >= 13 ? 13 : 12));
  RPFParams params = parse_rpf_params(base);
  double delta = pars.size() == 15 ? pars[13] : 0.1;
  double epsilon = pars.size() == 15 ? pars[14] : 0;
  if (pars.size() != 15) {
    if (!fitting)
      Rcpp::stop("ClassificationRPF requires 15 parameters, got %d", (int)pars.size());
    Rcout << "Wrong number of parameters - set to default." << std::endl;
    const int mode = params.split_structure_mode;
    params = RPFParams();
    params.split_structure_mode = mode;
  }
  return ClassificationRPF(loss, params, delta, epsilon);
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

void connect_to_r(RandomPlantedForest &core)
{
  core.verbose_out = &Rcout;
  core.seed_source = seed_from_r;
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
