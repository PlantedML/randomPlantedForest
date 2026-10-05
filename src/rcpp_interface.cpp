#include "rcpp_interface.h"

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
  if (pars.size() != 12 && pars.size() != 13)
    Rcpp::stop("RandomPlantedForest requires 12 or 13 parameters, got %d", (int)pars.size());
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
  if (pars.size() == 13)
    p.split_structure_mode = static_cast<int>(pars[12]);
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
  if (pars.size() != 15)
  {
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
  for (const FamilyExport &family : model)
  {
    List variables, family_values, family_intervals;
    for (const TreeExport &tree : family)
    {
      variables.push_back(IntegerVector(tree.variables.begin(), tree.variables.end()));
      List tree_values, tree_intervals;
      for (size_t k = 0; k < tree.values.size(); ++k)
      {
        tree_values.push_back(NumericVector(tree.values[k].begin(), tree.values[k].end()));
        NumericMatrix leaf_intervals(2, feature_size);
        for (int l = 0; l < feature_size; ++l)
        {
          leaf_intervals(0, l) = tree.intervals[k][l].first;
          leaf_intervals(1, l) = tree.intervals[k][l].second;
        }
        tree_intervals.push_back(leaf_intervals);
      }
      family_values.push_back(tree_values);
      family_intervals.push_back(tree_intervals);
    }
    out.push_back(List::create(Named("variables") = variables, Named("values") = family_values,
                               Named("intervals") = family_intervals));
  }
  return out;
}

std::vector<FamilyExport> list_to_model(const List &model, int feature_size)
{
  std::vector<FamilyExport> out;
  for (int i = 0; i < model.size(); ++i)
  {
    List family = model[i];
    List variables = family["variables"];
    List values = family["values"];
    List intervals = family["intervals"];
    FamilyExport fam;
    for (int j = 0; j < variables.size(); ++j)
    {
      TreeExport tree;
      IntegerVector tree_variables = variables[j];
      tree.variables = std::set<int>(tree_variables.begin(), tree_variables.end());
      List tree_values = values[j];
      List tree_intervals = intervals[j];
      for (int k = 0; k < tree_values.size(); ++k)
      {
        tree.values.push_back(as<std::vector<double>>(tree_values[k]));
        NumericMatrix leaf_intervals = tree_intervals[k];
        if (leaf_intervals.nrow() < 2 || leaf_intervals.ncol() < feature_size)
          Rcpp::stop("Corrupt model data: leaf interval matrix has dimensions %dx%d, expected 2x%d.",
                     leaf_intervals.nrow(), leaf_intervals.ncol(), feature_size);
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
  for (const GridFamilyExport &family : grid)
  {
    List lim_list;
    for (const auto &v : family.lim_list)
      lim_list.push_back(wrap(v));
    List trees;
    for (const GridTreeExport &tree : family.trees)
    {
      NumericMatrix values(tree.values.size(), value_size);
      for (size_t e = 0; e < tree.values.size(); ++e)
        for (int p = 0; p < value_size; ++p)
          values(e, p) = tree.values[e][p];
      trees.push_back(List::create(Named("variables") = IntegerVector(tree.variables.begin(), tree.variables.end()),
                                   Named("dims") = wrap(tree.dims),
                                   Named("values") = values));
    }
    families.push_back(List::create(Named("lim_list") = lim_list, Named("trees") = trees));
  }
  return families;
}

std::vector<GridFamilyExport> list_to_grid(const List &grid)
{
  std::vector<GridFamilyExport> out;
  for (int i = 0; i < grid.size(); ++i)
  {
    List family = grid[i];
    GridFamilyExport fam;
    List lim_list = family["lim_list"];
    for (int l = 0; l < lim_list.size(); ++l)
      fam.lim_list.push_back(as<std::vector<double>>(lim_list[l]));
    List trees = family["trees"];
    for (int j = 0; j < trees.size(); ++j)
    {
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

void connect_to_r(RandomPlantedForest &core)
{
  core.verbose_out = &Rcout;
  core.seed_source = []
  {
    // Two 32-bit draws from R's RNG composed into one 64-bit seed
    std::uint64_t hi = static_cast<std::uint64_t>(R::runif(0.0, 4294967296.0));
    std::uint64_t lo = static_cast<std::uint64_t>(R::runif(0.0, 4294967296.0));
    return (hi << 32) ^ lo;
  };
}
