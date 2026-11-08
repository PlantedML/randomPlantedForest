// Out-of-bag (OOB) path along the splits (regression). While a family grows, its prediction at each of its OOB rows
// is updated split by split, as the residuals of its bootstrap rows are: every tracked row sits in a list per leaf, a
// split moves the rows of the split leaf to the two new leaves and adds M_s or M_b to their predictions. At each
// checkpoint the predictions are copied; after the family is grown the copies are added to per-row sums on the calling
// thread, in family order. The forest grown with k splits and the same seed is exactly the first k splits of a longer
// fit (nothing in the fit depends on the total number of splits), so the sums at checkpoint k are the OOB sums of
// that forest.
#include "path_tracker.hpp"
#include <numeric>

void RandomPlantedForest::set_oob_path(const std::vector<int> &checkpoints)
{
  for (size_t c = 0; c < checkpoints.size(); ++c) {
    if (checkpoints[c] < 1 || checkpoints[c] > n_splits || (c > 0 && checkpoints[c] <= checkpoints[c - 1]))
      throw std::invalid_argument("OOB path checkpoints must be increasing split counts in 1..splits.");
  }
  path_checkpoints_ = checkpoints;
  path_sum_[0].clear();
  path_sum_[1].clear();
}

void RandomPlantedForest::path_reset()
{
  if (!path_checkpoints_.empty() && path_checkpoints_.back() > n_splits)
    throw std::invalid_argument("OOB path checkpoints must be increasing split counts in 1..splits.");
  path_family_.assign(path_checkpoints_.empty() ? 0 : (size_t)n_trees, nullptr);
  const size_t len = path_checkpoints_.size() * (size_t)sample_size * value_size;
  for (int h = 0; h < 2; ++h)
    path_sum_[h].assign(len, 0.0);
}

std::shared_ptr<RandomPlantedForest::PathTracker> RandomPlantedForest::path_start(size_t n)
{
  if (path_checkpoints_.empty())
    return nullptr;
  auto t = std::make_shared<PathTracker>();
  for (int i = 0; i < sample_size; ++i)
    if (inbag_[n][(size_t)i] == 0)
      t->rows.push_back(i);
  t->pred.assign(t->rows.size() * value_size, 0.0);
  std::vector<int> all(t->rows.size());
  std::iota(all.begin(), all.end(), 0);
  t->in_leaf[std::set<int>{0}] = {all};
  t->snaps.reserve(path_checkpoints_.size() * t->pred.size());
  path_family_[n] = t;
  return t;
}

void RandomPlantedForest::path_split(PathTracker &trk, const Split &split) const
{
  const DecisionTree &T = *split.tree_index;
  auto &lists = trk.in_leaf[T.split_dims];
  const size_t li = (size_t)(split.leaf_index - &T.leaves[0]);
  if (lists.size() < T.leaves.size())
    lists.resize(T.leaves.size());
  const size_t k = (size_t)split.split_coordinate - 1;
  std::vector<int> S, B;
  for (int p : lists[li]) {
    // the rule of create_tree_family: below the split point to the smaller half
    const bool small = X[(size_t)trk.rows[(size_t)p]][k] < split.split_point;
    (small ? S : B).push_back(p);
    const std::vector<double> &M = small ? split.M_s : split.M_b;
    for (size_t v = 0; v < value_size; ++v)
      trk.pred[(size_t)p * value_size + v] += M[v];
  }
  // the leaf bookkeeping of create_tree_family: the split leaf becomes the larger half and the smaller half is
  // appended to its tree, or both halves are appended (smaller first) to the tree of the resulting dims
  if (T.split_dims.count(split.split_coordinate) && delete_leaves) {
    lists[li] = std::move(B);
    lists.push_back(std::move(S));
  } else {
    std::set<int> dims = T.split_dims;
    dims.insert(split.split_coordinate);
    dims.erase(0);
    auto &dst = trk.in_leaf[dims];
    dst.push_back(std::move(S));
    dst.push_back(std::move(B));
  }
}

void RandomPlantedForest::path_step(PathTracker &trk, int k) const
{
  if (trk.next < path_checkpoints_.size() && path_checkpoints_[trk.next] == k) {
    trk.snaps.insert(trk.snaps.end(), trk.pred.begin(), trk.pred.end());
    ++trk.next;
  }
}

void RandomPlantedForest::path_merge(size_t n)
{
  if (n >= path_family_.size() || !path_family_[n])
    return;
  PathTracker &t = *path_family_[n];
  path_step(t, n_splits);
  const size_t C = path_checkpoints_.size(), no = t.rows.size(), vs = value_size, N = (size_t)sample_size;
  if (t.snaps.size() != C * no * vs)
    throw std::logic_error("OOB path: a family missed a checkpoint.");
  std::vector<double> &sum = path_sum_[n % 2];
  for (size_t c = 0; c < C; ++c)
    for (size_t v = 0; v < vs; ++v) {
      double *col = &sum[(c * vs + v) * N];
      for (size_t p = 0; p < no; ++p)
        col[(size_t)t.rows[p]] += t.snaps[(c * no + p) * vs + v];
    }
  path_family_[n].reset();
}

bool RandomPlantedForest::has_oob_path() const
{
  const size_t len = path_checkpoints_.size() * (size_t)sample_size * value_size;
  return !path_checkpoints_.empty() && path_sum_[0].size() == len && inbag_.size() == (size_t)n_trees;
}

OobPath RandomPlantedForest::oob_path() const
{
  if (!has_oob_path())
    throw std::invalid_argument("The forest has no OOB path (fit it with checkpoints).");
  OobPath out;
  out.checkpoints = path_checkpoints_;
  out.sum0 = path_sum_[0];
  out.sum1 = path_sum_[1];
  const size_t N = (size_t)sample_size;
  out.count.assign(2 * N, 0);
  for (size_t f = 0; f < inbag_.size(); ++f)
    for (size_t i = 0; i < N; ++i)
      if (inbag_[f][i] == 0)
        ++out.count[(f % 2) * N + i];
  return out;
}
