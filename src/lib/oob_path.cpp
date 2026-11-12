// Out-of-bag (OOB) path along the splits. While a family grows, its prediction at each of its OOB rows is updated split
// by split, as the residuals of its bootstrap rows are: every tracked row sits in a list per leaf, a split moves the
// rows of the split leaf to the two new leaves and adds M_s or M_b to their predictions (a classification forest: its
// leaves' increments, value_size columns; for logit and exponential on the link scale). At each
// checkpoint the predictions are copied; after the family is grown the copies are added to per-row sums on the calling
// thread, in family order. The forest grown with k splits and the same seed is exactly the first k splits of a longer
// fit (nothing in the fit depends on the total number of splits), so the sums at checkpoint k are the OOB sums of
// that forest.
//
// The adaptive fit (set_adaptive(), fit_alone() below) uses the path to choose the number of splits while growing.
#include "path_tracker.hpp"
#include "internal_utils.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <numeric>
#include <thread>
#ifdef __APPLE__
#include <malloc/malloc.h>
#endif

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
  // the adaptive fit passes the predictions to its stopper and keeps no sums
  const size_t C = alone_ ? 0 : path_checkpoints_.size();
  const size_t len = C * (size_t)sample_size * value_size;
  for (int h = 0; h < 2; ++h)
    path_sum_[h].assign(len, 0.0);
}

std::shared_ptr<RandomPlantedForest::PathTracker>
RandomPlantedForest::path_start(size_t n, TreeFamily *fam, const std::vector<std::vector<double>> *Y)
{
  if (path_checkpoints_.empty())
    return nullptr;
  auto t = std::make_shared<PathTracker>();
  // the fit of the null leaf, as fit_null_leaf() gives it to a family that made no split (from the same rows and Y):
  // the path tracks a family without a split at it (path_snapshot()), and the adaptive fit's rebuild gives it to a
  // family without a split at the chosen checkpoint
  if (Y && fam) {
    auto null_tree = fam->find(std::set<int>{0});
    if (null_tree != fam->end() && !null_tree->second->leaves.empty() &&
        !null_tree->second->leaves[0].individuals.empty())
      t->null_value = null_leaf_value(*Y, null_tree->second->leaves[0].individuals);
  }
  t->family = n;
  t->fam = fam;
  for (int i = 0; i < sample_size; ++i)
    if (inbag_[n][(size_t)i] == 0)
      t->rows.push_back(i);
  t->pred.assign(t->rows.size() * value_size, 0.0);
  std::vector<int> all(t->rows.size());
  std::iota(all.begin(), all.end(), 0);
  t->in_leaf[std::set<int>{0}] = {all};
  if (!alone_)
    t->snaps.reserve(path_checkpoints_.size() * t->pred.size());
  path_family_[n] = t;
  return t;
}

void RandomPlantedForest::path_split(
    PathTracker &trk, const Split &split, const std::vector<double> *inc_s, const std::vector<double> *inc_b
) const
{
  trk.has_split = true;
  const std::vector<double> &Ms = inc_s ? *inc_s : split.M_s, &Mb = inc_b ? *inc_b : split.M_b;
  const DecisionTree &T = *split.tree_index;
  auto &lists = trk.in_leaf[T.split_dims];
  const size_t li = (size_t)(split.leaf_index - &T.leaves[0]);
  if (lists.size() < T.leaves.size())
    lists.resize(T.leaves.size());
  const size_t k = (size_t)split.split_coordinate - 1;
  if (alone_) {
    // the adaptive fit: the record to rebuild the family's leaves at a checkpoint (alone_rebuild())
    auto it = trk.dims_id.find(T.split_dims);
    if (it == trk.dims_id.end()) {
      it = trk.dims_id.emplace(T.split_dims, (int32_t)trk.dims_list.size()).first;
      trk.dims_list.push_back(T.split_dims);
    }
    trk.rec.push_back(
        RpfSplitRecord{it->second, (int32_t)li, (int32_t)split.split_coordinate, split.split_point, Ms[0], Mb[0]}
    );
    if (value_size > 1) {
      trk.rec_vals.insert(trk.rec_vals.end(), Ms.begin(), Ms.begin() + (std::ptrdiff_t)value_size);
      trk.rec_vals.insert(trk.rec_vals.end(), Mb.begin(), Mb.begin() + (std::ptrdiff_t)value_size);
    }
  }
  std::vector<int> S, B;
  for (int p : lists[li]) {
    // the rule of create_tree_family: below the split point to the smaller half
    const bool small = X[(size_t)trk.rows[(size_t)p]][k] < split.split_point;
    (small ? S : B).push_back(p);
    const std::vector<double> &M = small ? Ms : Mb;
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

// The family's predictions at its OOB rows now, for a checkpoint. A family that has made no split predicts its null
// leaf's fit (fit_null_leaf(): the fit of its rows, e.g. their mean, or logit's log-odds of its sample), not the 0
// its rows start from: the path holds what a forest with so many splits predicts. A family with a split: the tracked
// predictions (its null leaf stays at 0).
std::vector<double> RandomPlantedForest::path_snapshot(const PathTracker &trk) const
{
  if (trk.has_split || trk.null_value.size() != value_size)
    return trk.pred;
  std::vector<double> out(trk.pred.size());
  for (size_t p = 0; p < trk.rows.size(); ++p)
    for (size_t v = 0; v < value_size; ++v)
      out[p * value_size + v] = trk.null_value[v];
  return out;
}

bool RandomPlantedForest::path_step(PathTracker &trk, int k) const
{
  if (alone_)
    return alone_step(trk, (size_t)k);
  if (trk.next < path_checkpoints_.size() && path_checkpoints_[trk.next] == k) {
    const std::vector<double> now = path_snapshot(trk);
    trk.snaps.insert(trk.snaps.end(), now.begin(), now.end());
    ++trk.next;
  }
  return true;
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

// -------------------------------------------------------------------------------------------------------------------
// The adaptive fit: choose the number of splits out of bag while growing, one family at a time.
//
// The stopper follows the forest's corrected OOB error checkpoint by checkpoint, as if all families grew together, and
// stops once its lowest value is `patience` checkpoints back, or once that value has improved by at most `tolerance`
// paired standard errors over the last `patience` checkpoints (a plateau; and at least `floor` splits are grown), or at
// the last checkpoint. Growing all families together would hold every family's growing state (bootstrap copy,
// residuals, leaves with their rows and caches, candidate pool) until the stop. Here at most `nthreads` families grow
// at a time, each alone and to a bound (a number of checkpoints). A family keeps its OOB predictions at the checkpoints
// it reached and a record of its splits (RpfSplitRecord, 40 bytes a split, and 16 bytes per score column more with
// value_size > 1). A family's splits depend only on its seed (its own random stream; the budget, the last checkpoint,
// is the same), so its predictions at a checkpoint are those of growing all families together, and the record rebuilds
// its leaves there.
//
// 1. Pilot: the first `nthreads` families grow together, checkpoint by checkpoint, until the stopper's rule on their
//    own corrected OOB error (AloneFit::estimate(), the same patience, tolerance and floor) stops. That is the first
//    bound.
// 2. The other families grow alone to the bound, `nthreads` at a time; after each one the estimate is made again from
//    all families so far, and the bound moves to it. A family that reaches the bound pauses there (its thread waits,
//    holding its growing state with the per-leaf caches freed) while the paused families fit in the pause budget
//    (`pause_gb`; their size is estimated), and otherwise finishes: its growing state is freed. When the bound rises,
//    the paused families below it grow on (cheaply); the finished ones are grown again from their seed (the same
//    splits, further) once every family has been grown. A family pauses at most `patience` checkpoints (one doubling)
//    past the pilot's bound: on a flat OOB curve the estimate from part of the families keeps finding new minima and
//    would grow every family far past the stop (a three-class logit forest of 200 families on 100 rows grew 2.7 times
//    the splits it needed), while a paused family that falls short of the stop grows on cheaply in step 3. A family
//    that cannot pause grows to the uncapped estimate, as a regrowth from its seed is not cheap.
// 3. The stopper is then called checkpoint by checkpoint with the OOB path there (the families' predictions summed in
//    family order, as fit() sums them). Where it goes on past the bound, the families grow further first (to the
//    earliest stop the curve so far allows).
// 4. The forest at the checkpoint the stopper chose is rebuilt from the records (alone_rebuild()).
// So the stopper's calls and decisions do not depend on the bound, the threads or the pause budget, and the forest is
// the one a fit with the chosen number of splits and the same seed gives. The bound decides only how far families are
// grown and how often. The memory is that of `nthreads` growing families, the paused ones within the budget and the
// records, instead of every family's growing state.

namespace {

// Give the memory freed by the fit back to the system (macOS keeps freed small blocks in the process, so the growing
// state freed after the fit would otherwise stay in its resident size)
void release_free_memory()
{
#ifdef __APPLE__
  malloc_zone_pressure_relief(nullptr, 0);
#endif
}

} // namespace

struct RandomPlantedForest::AloneFit {
  // A family of the adaptive fit. Once finished it keeps its OOB predictions at the checkpoints it reached and the
  // record of its splits, not its leaves.
  struct Family {
    std::vector<int> rows;                  /**< its OOB rows */
    std::vector<std::vector<double>> snaps; /**< per checkpoint reached: the predictions at its OOB rows */
    std::vector<RpfSplitRecord> rec;        /**< its splits in order */
    std::vector<double> rec_vals;           /**< value_size > 1: per record its 2 * value_size values */
    std::vector<uint32_t> rec_at;           /**< per checkpoint reached: the number of records up to it */
    std::vector<std::set<int>> dims;        /**< the split trees' dims (RpfSplitRecord::dims) */
    std::vector<double> null_value;         /**< the null leaf's fit, for a family without a split (path_start()) */
    size_t reached = 0;                     /**< the number of checkpoints reached */
    int state = 0;                          /**< 0 waiting, 1 growing, 2 grown, 3 queued to grow again, 4 paused */
    bool pilot = false;                     /**< one of the first `nthreads` families */
    bool permit = false;                    /**< it holds one of the `nthreads` permits (it is growing) */
    bool regrow = false;                    /**< it is grown again from its seed */
    bool release = false;                   /**< (paused) finish instead of growing on: the end of the fit */
    double pause_bytes = 0;                 /**< (paused) the estimated size of its growing state */
  };
  // the settings (set_adaptive())
  AdaptiveStopper stopper;
  int patience = 4; // the stopper's patience (checkpoints), tolerance (paired SEs) and floor (splits), for the bound
  double tolerance = 0;
  double floor = 0;
  double pause_budget = 0; // bytes
  // the fit
  std::mutex m;
  std::condition_variable cv; // the pool and the pilot wait on it
  // a paused family waits on its own (wake()): with one for all, every change woke every paused family's thread, and
  // with many families paused the system time of those wake-ups was many times the fit's, also at nthreads = 1.
  // A vector, not new[]: its length check bounds the allocation, where GCC's LTO warned -Walloc-size-larger-than
  std::vector<std::condition_variable> paused_cv;
  std::atomic<size_t> bound{0}; // a family grows to this many checkpoints (read at each checkpoint)
  std::atomic<size_t> bound_free{
      0};                   // the estimate without the cap on pausing (phase 2), for families that cannot pause
  size_t need_patience = 0; // estimate()'s bound by the patience alone where the curve ends before a stop
  std::atomic<bool> cancel{false};
  std::exception_ptr err;
  std::vector<Family> fam;
  // the pilot: the first families grow together, checkpoint by checkpoint, until the bound's estimate stops
  bool pilot_on = false;
  size_t pilot_n = 0, pilot_arrived = 0, pilot_gen = 0;
  // the pool: the next family not yet grown, the finished families queued to grow again, the families growing and
  // the free permits; one thread per family (a paused family waits on its own thread)
  size_t next_new = 0, running = 0;
  int permits = 1;
  std::deque<size_t> queue;
  std::vector<std::thread> threads;
  std::vector<size_t> start; // the families started first (the pilot)
  double paused_bytes = 0;
  bool estimating = true; // re-estimate the bound after each family (off once the stopper decides)
  size_t consumed = 0;    // checkpoints already passed to the stopper (their predictions are freed)
  // the bound's estimate: per checkpoint and row (and score column: [i * V + v]), the sums of the counted families' OOB
  // predictions over the families of even and of odd index and their numbers (empty: no family counted there yet), and
  // the sums of squares over all
  std::vector<std::vector<double>> s0, s1, q;
  std::vector<std::vector<int>> m0, m1;
  std::vector<size_t> counted;                 // per family: the checkpoints added to the estimate
  std::vector<size_t> n_counted;               // per checkpoint: the number of families counted there
  const RandomPlantedForest *forest = nullptr; // its oob_row_error(): a row's error
  size_t V = 1;                                // the score columns
  std::vector<int> cp;
  size_t B = 0, N = 0;
  AdaptiveStats stats;

  // The stopper's OOB sums, kept up to date instead of every family's predictions at every checkpoint until the
  // stop is decided. Per parity h of the family index and checkpoint c: the sum of the families' predictions there
  // (N * V, score column v at [v * N + i]; empty until the first family is added), added in family order as fit()
  // adds the OOB path: next_[h][c] is the next family to add. A family's predictions at a checkpoint are added, and
  // freed, as soon as it has delivered them (it paused or finished: Family::reached) and every earlier family of its
  // parity has been added there. flushed[b]: family b's checkpoints added so far (a prefix: an earlier checkpoint is
  // delivered no later than a later one).
  std::vector<std::vector<double>> sum_[2];
  std::vector<size_t> next_[2];
  std::vector<size_t> flushed;

  // Add what can be added (under the lock, after count(): a family's predictions are counted before they are freed)
  void flush()
  {
    for (size_t c = consumed; c < cp.size(); ++c)
      for (size_t h = 0; h < 2; ++h) {
        size_t &b = next_[h][c];
        while (b < B) {
          Family &F = fam[b];
          if (c >= F.reached || c >= F.snaps.size())
            break;
          std::vector<double> &pr = F.snaps[c];
          if (pr.size() != F.rows.size() * V)
            break; // (not there: nothing is added out of order)
          if (!F.rows.empty()) {
            std::vector<double> &sum = sum_[h][c];
            if (sum.empty())
              sum.assign(N * V, 0.0);
            for (size_t i = 0; i < F.rows.size(); ++i)
              for (size_t v = 0; v < V; ++v)
                sum[v * N + (size_t)F.rows[i]] += pr[i * V + v];
          }
          std::vector<double>().swap(pr);
          flushed[b] = c + 1;
          b += 2;
        }
      }
  }

  // Add family n's predictions at checkpoints [counted[n], to) to the estimate (under the lock)
  void count(size_t n, const std::vector<int> &rows, const std::vector<std::vector<double>> &snaps, size_t to)
  {
    for (size_t c = counted[n]; c < to && c < snaps.size(); ++c) {
      if (s0[c].empty()) {
        s0[c].assign(N * V, 0.0);
        s1[c].assign(N * V, 0.0);
        m0[c].assign(N, 0);
        m1[c].assign(N, 0);
        q[c].assign(N * V, 0.0);
      }
      std::vector<double> &s = (n % 2 ? s1 : s0)[c];
      std::vector<int> &mm = (n % 2 ? m1 : m0)[c];
      const std::vector<double> &p = snaps[c];
      if (p.size() < rows.size() * V)
        break; // (freed: already passed to the stopper)
      for (size_t i = 0; i < rows.size(); ++i) {
        const size_t r = (size_t)rows[i] * V;
        for (size_t v = 0; v < V; ++v) {
          const double x = p[i * V + v];
          s[r + v] += x;
          q[c][r + v] += x * x;
        }
        ++mm[(size_t)rows[i]];
      }
      ++n_counted[c];
      counted[n] = c + 1;
    }
  }

  // The bound: the number of checkpoints the forest would grow to by the stopper's rule (the lowest corrected OOB
  // error `patience` checkpoints back, or its gain over the last `patience` checkpoints within `tolerance` paired
  // standard errors; at least `floor` splits) on the families counted so far (at each checkpoint, those that reached
  // it; the curve ends where fewer than min(4, pilot size) are counted). Where every family is counted, the error is
  // the stopper's (so the bound is its stop, up to the rounding of the sums); elsewhere v_i is the sample variance over
  // the row's OOB families, on the rows with at least two (with few families the even / odd estimate, and the rows it
  // leaves uncorrected, stop too early). Where the curve ends before a stop, the earliest stop that could follow by
  // the patience alone. A row's error is the forest's oob_row_error() (regression: (y - g)^2 - (1 / m - 1 / B) v; with
  // V > 1 columns and the sample variances, the diagonal of their covariance).
  size_t estimate()
  {
    const size_t C = cp.size();
    std::vector<double> e;
    const size_t min_fam = std::min<size_t>(4, std::max<size_t>(1, std::min(B, pilot_n)));
    std::vector<double> g(V), d(V), var(V);
    // row i's error at checkpoint c in `err`; false if the row is not counted there
    auto row_error = [&](size_t c, size_t i, double &err) -> bool {
      const bool all = n_counted[c] == B || n_counted[c] < 2; // (one family: no variance across families)
      const int a = m0[c][i], b = m1[c][i], m = a + b;
      if (m == 0 || (!all && m < 2))
        return false;
      for (size_t v = 0; v < V; ++v) {
        const double g0 = s0[c][i * V + v], g1 = s1[c][i * V + v];
        g[v] = (g0 + g1) / m;
        if (!all)
          var[v] = (q[c][i * V + v] - m * g[v] * g[v]) / (m - 1);
        else if (a > 0 && b > 0)
          d[v] = g0 / a - g1 / b;
      }
      const double coef = 1.0 / m - 1.0 / (double)B;
      if (!all)
        err = forest->oob_row_error(i, g.data(), nullptr, 0.0, var.data(), coef);
      else if (a > 0 && b > 0)
        err = forest->oob_row_error(i, g.data(), d.data(), 1.0 / a + 1.0 / b, nullptr, coef);
      else
        err = forest->oob_row_error(i, g.data(), nullptr, 0.0, nullptr, coef);
      return true;
    };
    for (size_t c = 0; c < C && !s0[c].empty() && n_counted[c] >= min_fam; ++c) {
      double sum = 0, err = 0;
      size_t n_ok = 0;
      for (size_t i = 0; i < N; ++i) {
        if (!row_error(c, i, err))
          continue;
        sum += err;
        ++n_ok;
      }
      if (n_ok == 0)
        break;
      e.push_back(sum / (double)n_ok);
    }
    // the plateau: the gain of the lowest error from checkpoint a to checkpoint b is at most `tolerance` standard
    // errors of the rows' paired differences (over the rows counted at both)
    auto plateau = [&](size_t a, size_t b) -> bool {
      const double gain = e[a] - e[b];
      if (gain <= 0)
        return true;
      double mean = 0, m2 = 0, ea = 0, eb = 0;
      size_t n = 0;
      for (size_t i = 0; i < N; ++i) {
        if (!row_error(a, i, ea) || !row_error(b, i, eb))
          continue;
        const double x = ea - eb, dx = x - mean;
        ++n;
        mean += dx / (double)n;
        m2 += dx * (x - mean);
      }
      if (n < 2)
        return false;
      return gain <= tolerance * std::sqrt(m2 / (double)(n - 1) / (double)n);
    };
    size_t need = C;
    const size_t L = e.size();
    if (L == 0)
      return need_patience = std::min<size_t>(C, 1);
    double mn = e[0], mn_lag = e[0];
    size_t best = 0, best_lag = 0; // the first lowest error up to c, and up to c - patience
    // the stop by the patience alone (need_patience) and by the patience or the plateau (need). Before every family is
    // counted the rows' errors are noisier, so the plateau's SE is larger and its estimate can stop early: the bound
    // for families that cannot pause (a regrowth from the seed is not cheap) is the patience's.
    bool stopped = false, stopped_patience = false;
    need_patience = C;
    for (size_t c = 0; c < L && !stopped_patience; ++c) {
      mn = std::min(mn, e[c]);
      while (!(e[best] <= mn * (1 + 1e-12) || e[best] == mn))
        ++best;
      if ((c - best >= (size_t)patience && cp[c] >= floor) || c + 1 == C) {
        need_patience = c + 1;
        stopped_patience = true;
        if (!stopped)
          need = c + 1;
        stopped = true;
      } else if (!stopped && tolerance > 0 && c >= (size_t)patience) {
        mn_lag = std::min(mn_lag, e[c - patience]);
        while (!(e[best_lag] <= mn_lag * (1 + 1e-12) || e[best_lag] == mn_lag))
          ++best_lag;
        if (cp[c] >= floor && plateau(best_lag, best)) {
          need = c + 1;
          stopped = true;
        }
      }
    }
    // where the curve ends before a stop: the earliest stop by the patience alone; with the plateau on, the
    // checkpoint after the next one (the next one, where the plateau could stop first, made families pause and resume
    // at every checkpoint on long paths: three times the pauses and up to 14% more time on a three-class exponential
    // fit of 566 splits; two further gave back most of the time saved on short logit paths)
    auto earliest = [&](size_t c) {
      while (c < C && cp[c] < floor)
        ++c;
      return std::min(C, c + 1);
    };
    if (!stopped_patience)
      need_patience = earliest(std::max(L, best + (size_t)patience));
    if (!stopped)
      need = tolerance > 0 ? std::min(earliest(L + 1), need_patience) : need_patience;
    stats.bound_trace.push_back((int)need);
    return need;
  }

  // After a change (under the lock): wake the pool and the pilot, and of the paused families those whose wait ends: all
  // of them at the end of the fit, otherwise those below the bound, in family order, as many as there are free permits
  void wake()
  {
    cv.notify_all();
    if (paused_cv.empty())
      return;
    int free = permits;
    for (size_t b = 0; b < fam.size(); ++b) {
      const Family &F = fam[b];
      if (F.state != 4)
        continue;
      if (cancel.load() || F.release)
        paused_cv[b].notify_one();
      else if (free > 0 && bound.load() > F.reached) {
        paused_cv[b].notify_one();
        --free;
      }
    }
  }

  // After a family is counted (under the lock): the bound moves to the new estimate, at most `patience` checkpoints
  // past the pilot's bound (bound_free keeps the estimate by the patience alone, uncapped), and the finished families
  // short of it are queued to grow again (the paused ones resume on their own)
  void after_count()
  {
    if (!estimating || pilot_on)
      return;
    size_t need = estimate();
    bound_free = need_patience;
    if (stats.pilot_bound > 0)
      need = std::min(need, (size_t)stats.pilot_bound + (size_t)patience);
    bound = need;
    for (size_t b = 0; b < fam.size(); ++b)
      if (fam[b].state == 2 && fam[b].reached < need) {
        fam[b].state = 3;
        queue.push_back(b);
      }
  }
};

void RandomPlantedForest::set_adaptive(
    const std::vector<int> &checkpoints,
    AdaptiveStopper stopper,
    int patience,
    double tolerance,
    double floor,
    double pause_gb
)
{
  if (deterministic)
    throw std::invalid_argument("The adaptive fit needs out-of-bag rows (deterministic = FALSE).");
  if (checkpoints.empty() || checkpoints.back() != n_splits)
    throw std::invalid_argument("The adaptive fit's last checkpoint must be its number of splits.");
  if (!stopper || patience < 1 || !(tolerance >= 0) || pause_gb < 0)
    throw std::invalid_argument("The adaptive fit needs a stopper, patience >= 1, tolerance >= 0 and pause_gb >= 0.");
  set_oob_path(checkpoints);
  alone_ = std::make_shared<AloneFit>();
  alone_->stopper = std::move(stopper);
  alone_->patience = patience;
  alone_->tolerance = tolerance;
  alone_->floor = floor;
  alone_->pause_budget = pause_gb * 1e9;
}

AdaptiveStats RandomPlantedForest::adaptive_stats() const
{
  return alone_stats_;
}

void RandomPlantedForest::free_leaf_caches(TreeFamily &family)
{
  for (auto &kv : family)
    for (auto &leaf : kv.second->leaves)
      leaf.free_cache();
}

// `k` splits are done. At a checkpoint c: keep the predictions and the number of records; a pilot family waits there
// for the other pilot families. Below the bound (or at the last checkpoint) the family grows on; at the bound it
// pauses if the paused families fit in the budget, otherwise it grows one checkpoint further and finishes.
bool RandomPlantedForest::alone_step(PathTracker &trk, size_t k) const
{
  AloneFit &A = *alone_;
  const size_t C = path_checkpoints_.size();
  auto step = [&]() -> bool {
    if (trk.next >= C || (size_t)path_checkpoints_[trk.next] != k)
      return !A.cancel;
    const size_t c = trk.next++;
    trk.csnaps.push_back(path_snapshot(trk));
    trk.rec_at.push_back((uint32_t)trk.rec.size());
    if (A.cancel)
      return false;
    const size_t n = trk.family;
    if (A.fam[n].pilot) {
      std::unique_lock<std::mutex> lk(A.m);
      if (A.pilot_on) {
        const size_t gen = A.pilot_gen;
        if (++A.pilot_arrived >= A.pilot_n) {
          // the last pilot family at checkpoint c: count them all there, and end the pilot if the estimate stops
          for (size_t b = 0; b < A.pilot_n; ++b) {
            const PathTracker &t = *path_family_[b];
            A.count(b, t.rows, t.csnaps, c + 1);
          }
          if (A.estimate() <= c + 1) {
            A.bound = c + 1;
            A.stats.pilot_bound = (int)(c + 1);
            A.pilot_on = false;
          }
          A.pilot_arrived = 0;
          ++A.pilot_gen;
          A.wake();
        } else
          A.cv.wait(lk, [&] { return A.pilot_gen != gen || A.cancel.load(); });
      }
    }
    if (A.cancel)
      return false;
    if (c + 1 < A.bound.load() || c + 1 == C)
      return c + 1 < C;
    // at the bound: the growing state without the per-leaf caches (the bootstrap copy and residuals, and per leaf its
    // struct, intervals, rows and roughly its entries in the candidate pool)
    double bytes = 0;
    if (trk.fam) {
      bytes = (double)sample_size * (double)(64 + 8 * feature_size + 48);
      for (const auto &kv : *trk.fam)
        for (const auto &leaf : kv.second->leaves)
          bytes += (double)sizeof(Leaf) + 16.0 * (double)leaf.intervals.size() + 4.0 * (double)leaf.individuals.size() +
                   64.0 + 136.0 * (double)std::min(feature_size, 5);
    }
    std::unique_lock<std::mutex> lk(A.m);
    AloneFit::Family &F = A.fam[n];
    // no room to pause: finish one checkpoint past the bound, uncapped while the bound is estimated (growing a
    // checkpoint further costs a fifth more splits; falling short of the stop, a regrowth from the seed)
    if (A.paused_bytes + bytes > A.pause_budget || !trk.fam) {
      const size_t to = A.estimating ? std::max(A.bound.load(), A.bound_free.load()) : A.bound.load();
      return c + 1 < to + 1 && !A.cancel;
    }
    // pause: the predictions so far go to the family's record (the stopper reads them there), the bound is estimated
    // again, and the thread gives up its permit until the bound rises past it (or the fit ends)
    F.snaps.resize(std::max(F.snaps.size(), trk.csnaps.size()));
    for (size_t j = std::max(A.consumed, A.flushed[n]); j < trk.csnaps.size(); ++j)
      if (!trk.csnaps[j].empty())
        F.snaps[j] = std::move(trk.csnaps[j]);
    if (F.rows.empty())
      F.rows = trk.rows;
    F.reached = c + 1;
    F.state = 4;
    F.pause_bytes = bytes;
    A.paused_bytes += bytes;
    A.stats.max_paused_gb = std::max(A.stats.max_paused_gb, A.paused_bytes / 1e9);
    ++A.stats.pauses;
    F.permit = false;
    ++A.permits;
    --A.running;
    A.count(n, F.rows, F.snaps, F.reached);
    A.flush();
    A.after_count();
    A.wake();
    lk.unlock();
    free_leaf_caches(*trk.fam);
    lk.lock();
    A.paused_cv[n].wait(lk, [&] { return A.cancel.load() || F.release || (A.bound.load() > c + 1 && A.permits > 0); });
    A.paused_bytes -= F.pause_bytes;
    F.pause_bytes = 0;
    if (A.cancel || F.release)
      return false;
    // grow on
    --A.permits;
    ++A.running;
    F.permit = true;
    F.state = 1;
    ++A.stats.resumes;
    return true;
  };
  const bool go = step();
  if (!go)
    trk.stopped = true;
  return go;
}

void RandomPlantedForest::alone_finish(size_t n)
{
  AloneFit &A = *alone_;
  std::shared_ptr<PathTracker> tp = path_family_[n];
  // the last checkpoint (the split loop ends before it), unless the family stopped earlier
  if (tp && !tp->stopped && !A.cancel)
    alone_step(*tp, (size_t)n_splits);
  path_family_[n].reset();
  TreeFamily().swap(tree_families[n]); // the family is kept as its record
  std::lock_guard<std::mutex> lk(A.m);
  AloneFit::Family &f = A.fam[n];
  if (f.permit) {
    f.permit = false;
    ++A.permits;
    --A.running;
  }
  if (A.cancel || !tp) {
    A.wake();
    return;
  }
  PathTracker &t = *tp;
  if (f.regrow)
    ++A.stats.regrowths;
  f.regrow = false;
  ++A.stats.growths;
  f.state = 2;
  f.rows = std::move(t.rows);
  f.snaps.resize(std::max(f.snaps.size(), t.csnaps.size()));
  for (size_t j = std::max(A.consumed, A.flushed[n]); j < t.csnaps.size(); ++j)
    if (!t.csnaps[j].empty())
      f.snaps[j] = std::move(t.csnaps[j]);
  f.rec = std::move(t.rec);
  f.rec.shrink_to_fit();
  f.rec_vals = std::move(t.rec_vals);
  f.rec_vals.shrink_to_fit();
  f.rec_at = std::move(t.rec_at);
  f.dims = std::move(t.dims_list);
  f.null_value = std::move(t.null_value);
  f.reached = t.next;
  A.stats.splits_grown += f.reached ? (double)path_checkpoints_[f.reached - 1] : 0.0;
  A.count(n, f.rows, f.snaps, f.reached);
  A.flush();
  A.after_count();
  A.wake();
}

// Grow the families until every one has reached the bound (paused or finished): the pilot first (if on), then the
// families not yet grown, then the finished ones short of the bound, grown again; at most `nthreads` at a time (the
// permits), each on a thread of its own (a paused family keeps its thread). Resumed families go before new ones.
void RandomPlantedForest::alone_pool(const std::vector<Leaf> &initial_leaves)
{
  AloneFit &A = *alone_;
  const size_t B = (size_t)n_trees;
  auto run = [this, &initial_leaves](size_t n) {
    AloneFit &A = *alone_;
    std::mt19937_64 rng_local(tree_seeds_[n]);
    std::mt19937_64 *prev = rpf_utils::swap_tls_rng(&rng_local);
    // an exception escaping a std::thread would terminate R: it stops the fit and is rethrown on the calling thread
    try {
      create_tree_family(initial_leaves, n);
      rpf_utils::swap_tls_rng(prev);
      alone_finish(n);
    } catch (...) {
      rpf_utils::swap_tls_rng(prev);
      std::lock_guard<std::mutex> lk(A.m);
      if (!A.err)
        A.err = std::current_exception();
      A.cancel = true;
      if (A.fam[n].permit) {
        A.fam[n].permit = false;
        ++A.permits;
        --A.running;
      }
      A.wake();
    }
  };
  auto launch = [&](size_t n) { // under the lock, with a permit taken
    AloneFit::Family &F = A.fam[n];
    F.permit = true;
    F.state = 1;
    ++A.running;
    if (A.threads[n].joinable())
      A.threads[n].join(); // (its earlier growth has finished)
    A.threads[n] = std::thread(run, n);
  };
  // a paused family short of the bound (it resumes on its own as soon as a permit is free)
  auto resume_due = [&]() {
    for (size_t b = 0; b < B; ++b)
      if (A.fam[b].state == 4 && A.fam[b].reached < A.bound.load())
        return true;
    return false;
  };
  const auto t0 = std::chrono::steady_clock::now();
  try {
    for (;;) {
      {
        std::unique_lock<std::mutex> lk(A.m);
        for (;;) {
          if (A.cancel)
            break;
          if (A.permits > 0 && !resume_due()) {
            if (!A.start.empty()) {
              const size_t n = A.start.back();
              A.start.pop_back();
              --A.permits;
              launch(n);
              continue;
            }
            if (!A.pilot_on && A.next_new < B) {
              --A.permits;
              launch(A.next_new++);
              continue;
            }
            if (A.next_new >= B && !A.queue.empty()) {
              const size_t n = A.queue.front();
              A.queue.pop_front();
              if (A.fam[n].reached >= A.bound.load()) {
                A.fam[n].state = 2;
                continue;
              }
              A.fam[n].regrow = true;
              --A.permits;
              launch(n);
              continue;
            }
          }
          break;
        }
        if (A.cancel)
          break;
        if (A.start.empty() && A.next_new >= B && A.queue.empty() && A.running == 0 && !resume_due())
          break;
        A.cv.wait_for(lk, std::chrono::milliseconds(100));
      }
      if (check_interrupt)
        check_interrupt();
    }
  } catch (...) {
    {
      std::lock_guard<std::mutex> lk(A.m);
      A.cancel = true;
      A.wake();
    }
    for (auto &th : A.threads)
      if (th.joinable())
        th.join();
    throw;
  }
  A.stats.seconds_grow += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (A.err) {
    for (auto &th : A.threads)
      if (th.joinable())
        th.join();
    std::rethrow_exception(A.err);
  }
}

void RandomPlantedForest::alone_end()
{
  AloneFit &A = *alone_;
  {
    std::lock_guard<std::mutex> lk(A.m);
    A.estimating = false;
    for (auto &F : A.fam)
      F.release = true;
    A.wake();
  }
  for (auto &th : A.threads)
    if (th.joinable())
      th.join();
}

// Every family's leaves at checkpoint c (values and intervals), rebuilt from its record with create_tree_family()'s
// bookkeeping of a split: the forest a fit with that many splits and the same seed gives (the leaves' rows and caches,
// which a finished family frees, are not rebuilt)
void RandomPlantedForest::alone_rebuild(size_t c)
{
  AloneFit &A = *alone_;
  const auto t0 = std::chrono::steady_clock::now();
  const size_t B = (size_t)n_trees;
  Leaf root;
  root.value.assign(value_size, 0.0);
  root.intervals.resize((size_t)feature_size);
  for (int i = 0; i < feature_size; ++i)
    root.intervals[(size_t)i] = Interval{lower_bounds[(size_t)i], upper_bounds[(size_t)i]};
  tree_families.assign(B, TreeFamily());
  drop_purify_counts();
  std::atomic<size_t> next{0};
  std::exception_ptr err;
  std::mutex err_m;
  auto work = [&]() {
    try {
      for (size_t b; (b = next++) < B;) {
        const AloneFit::Family &F = A.fam[b];
        if (F.reached <= c)
          throw std::logic_error("Adaptive fit: a family has not reached the chosen checkpoint.");
        TreeFamily fam;
        auto r0 = std::make_shared<DecisionTree>(std::set<int>{0});
        r0->leaves.push_back(root);
        fam.insert({std::set<int>{0}, r0});
        const size_t R = F.rec_at[c];
        // no split up to the checkpoint: the null leaf's fit, as create_tree_family() gives it (fit_null_leaf())
        if (R == 0 && F.null_value.size() == value_size)
          r0->leaves[0].value = F.null_value;
        for (size_t i = 0; i < R; ++i) {
          const RpfSplitRecord &r = F.rec[i];
          const std::set<int> &dims = F.dims[(size_t)r.dims];
          const std::shared_ptr<DecisionTree> T = fam.at(dims);
          const Leaf &parent = T->leaves.at((size_t)r.leaf);
          Leaf ls, lb;
          if (value_size == 1) {
            ls.value.assign(1, r.ms);
            lb.value.assign(1, r.mb);
          } else {
            const double *rv = F.rec_vals.data() + i * 2 * value_size;
            ls.value.assign(rv, rv + value_size);
            lb.value.assign(rv + value_size, rv + 2 * value_size);
          }
          ls.intervals = parent.intervals;
          lb.intervals = parent.intervals;
          ls.intervals[(size_t)r.coord - 1].second = r.point;
          lb.intervals[(size_t)r.coord - 1].first = r.point;
          if (dims.count(r.coord) && delete_leaves) {
            for (size_t v = 0; v < value_size; ++v) {
              ls.value[v] += parent.value[v];
              lb.value[v] += parent.value[v];
            }
            T->leaves[(size_t)r.leaf] = std::move(lb);
            T->leaves.push_back(std::move(ls));
          } else {
            std::set<int> res = dims;
            res.insert(r.coord);
            res.erase(0);
            std::shared_ptr<DecisionTree> &D = fam[res];
            if (!D)
              D = std::make_shared<DecisionTree>(res);
            D->leaves.push_back(std::move(ls));
            D->leaves.push_back(std::move(lb));
          }
        }
        for (auto &kv : fam)
          kv.second->leaves.shrink_to_fit(); // no spare capacity in the finished forest
        tree_families[b] = std::move(fam);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lk(err_m);
      if (!err)
        err = std::current_exception();
      next = B;
    }
  };
  const size_t W = std::max<size_t>(1, std::min<size_t>((size_t)std::max(1, nthreads), B));
  std::vector<std::thread> threads;
  for (size_t w = 1; w < W; ++w)
    threads.emplace_back(work);
  work();
  for (auto &th : threads)
    th.join();
  if (err)
    std::rethrow_exception(err);
  n_splits = path_checkpoints_[c];
  A.stats.seconds_rebuild += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void RandomPlantedForest::fit_alone(const std::vector<Leaf> &initial_leaves)
{
  AloneFit &A = *alone_;
  const size_t B = (size_t)n_trees, C = path_checkpoints_.size(), N = (size_t)sample_size, V = value_size;
  A.B = B;
  A.N = N;
  A.cp = path_checkpoints_;
  A.forest = this;
  A.V = V;
  A.fam.assign(B, AloneFit::Family());
  A.paused_cv = std::vector<std::condition_variable>(B);
  A.counted.assign(B, 0);
  for (auto *v : {&A.s0, &A.s1, &A.q})
    v->assign(C, {});
  for (auto *v : {&A.m0, &A.m1})
    v->assign(C, {});
  A.n_counted.assign(C, 0);
  for (size_t h = 0; h < 2; ++h) {
    A.sum_[h].assign(C, {});
    A.next_[h].assign(C, h);
  }
  A.flushed.assign(B, 0);
  A.bound = C;
  A.threads.clear();
  A.threads.resize(B);
  // the pilot: the first min(nthreads, ntrees) families
  const size_t W = std::max<size_t>(1, std::min<size_t>((size_t)std::max(1, nthreads), B));
  A.pilot_n = A.next_new = W;
  A.permits = (int)W;
  A.running = 0;
  A.pilot_on = true;
  A.start.clear();
  for (size_t b = W; b-- > 0;) {
    A.fam[b].pilot = true;
    A.start.push_back(b);
  }
  // the end, kept or not: the paused families finish, the threads are joined, the growing state and the predictions
  // are freed (with `keep`, the records stay for the rebuild)
  auto end = [&](bool keep) {
    if (!keep) {
      std::lock_guard<std::mutex> lk(A.m);
      A.cancel = true;
    }
    alone_end();
    std::vector<std::thread>().swap(A.threads);
    if (keep) {
      for (auto &F : A.fam) {
        std::vector<std::vector<double>>().swap(F.snaps);
        std::vector<int>().swap(F.rows);
      }
    } else
      std::vector<AloneFit::Family>().swap(A.fam);
    for (auto *v : {&A.s0, &A.s1, &A.q})
      std::vector<std::vector<double>>().swap(*v);
    for (auto *v : {&A.m0, &A.m1})
      std::vector<std::vector<int>>().swap(*v);
    for (size_t h = 0; h < 2; ++h)
      std::vector<std::vector<double>>().swap(A.sum_[h]);
    A.forest = nullptr;
    A.queue.clear();
    for (auto &t : path_family_)
      t.reset();
  };
  // after the fit, kept or not: the forest is a plain one (no path, no adaptive state), the stopper is released
  auto off = [&]() {
    alone_stats_ = A.stats;
    path_checkpoints_.clear();
    path_family_.clear();
    alone_.reset();
    release_free_memory();
  };
  size_t chosen = 0;
  try {
    alone_pool(initial_leaves);
    if (A.pilot_on)
      throw std::logic_error("Adaptive fit: the pilot did not stop.");
    // the OOB path at one checkpoint at a time, as the stopper gets it
    OobPath path;
    path.count.assign(2 * N, 0);
    for (size_t f = 0; f < B; ++f)
      for (size_t i = 0; i < N; ++i)
        if (inbag_[f][i] == 0)
          ++path.count[(f % 2) * N + i];
    bool stopped = false;
    for (size_t c = 0; c < C && !stopped; ++c) {
      bool lack = false;
      for (size_t b = 0; b < B; ++b)
        if (A.fam[b].reached <= c)
          lack = true;
      A.consumed = c;
      if (lack) {
        // the stopper goes on past the bound: every family grows further, to the earliest stop the curve so far
        // allows (the paused ones resume, the finished ones are grown again; then by the patience alone, as a
        // regrowth from the seed is not cheap)
        {
          std::lock_guard<std::mutex> lk(A.m);
          A.estimating = false;
          size_t need = A.estimate();
          for (size_t b = 0; b < B; ++b)
            if (A.fam[b].state == 2 && A.fam[b].reached <= c)
              need = A.need_patience;
          A.bound = std::min(C, std::max(c + 1, need));
          for (size_t b = 0; b < B; ++b)
            if (A.fam[b].state == 2 && A.fam[b].reached < A.bound) {
              A.fam[b].state = 3;
              A.queue.push_back(b);
            }
          A.wake();
        }
        alone_pool(initial_leaves);
      }
      const auto t0 = std::chrono::steady_clock::now();
      // every family's predictions at checkpoint c, summed over the families of even and of odd index in family order
      // (as fit() sums the OOB path; score column v at [v * N + i])
      path.checkpoints.assign(1, path_checkpoints_[c]);
      {
        // every family has delivered checkpoint c: its sums are complete (AloneFit::flush())
        std::lock_guard<std::mutex> lk(A.m);
        A.flush();
        if (A.next_[0][c] < B || A.next_[1][c] < B)
          throw std::logic_error("Adaptive fit: a family's predictions are missing at a checkpoint.");
        path.sum0 = std::move(A.sum_[0][c]);
        path.sum1 = std::move(A.sum_[1][c]);
        std::vector<double>().swap(A.sum_[0][c]);
        std::vector<double>().swap(A.sum_[1][c]);
      }
      if (path.sum0.empty())
        path.sum0.assign(N * V, 0.0);
      if (path.sum1.empty())
        path.sum1.assign(N * V, 0.0);
      const int sel = A.stopper(c, path);
      A.stats.seconds_stopper += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      for (auto &F : A.fam)
        if (c < F.snaps.size())
          std::vector<double>().swap(F.snaps[c]);
      A.consumed = c + 1;
      if (check_interrupt)
        check_interrupt();
      if (sel < 0 || sel > (int)c + 1)
        throw std::invalid_argument("The adaptive fit's stopper must return 0 or a checkpoint grown.");
      if (sel > 0) {
        chosen = (size_t)sel - 1;
        stopped = true;
      } else if (c + 1 == C)
        throw std::invalid_argument("The adaptive fit's stopper must choose a checkpoint at the last one.");
    }
    end(true);
    alone_rebuild(chosen);
  } catch (...) {
    if (!A.threads.empty())
      end(false);
    tree_families.assign((size_t)n_trees, TreeFamily());
    off();
    throw;
  }
  off();
}

// The regression criterion of a row: the squared error, corrected as the OOB path's table corrects it (see rpf.hpp)
double RandomPlantedForest::oob_row_error(
    size_t i, const double *g, const double *d, double inv, const double *var, double coef
) const
{
  const double r = Y[i][0] - g[0];
  double v = 0;
  if (d)
    v = d[0] * d[0] / inv;
  else if (var)
    v = var[0];
  return r * r - coef * v;
}

std::vector<double> RandomPlantedForest::oob_row_errors(
    const std::vector<double> &G,
    const std::vector<double> &D,
    const std::vector<double> &inv,
    const std::vector<double> &coef
) const
{
  const size_t n = (size_t)sample_size, V = value_size;
  if (Y.size() != n || G.size() != n * V || D.size() != n * V || inv.size() != n || coef.size() != n)
    throw std::invalid_argument("oob_row_errors(): G and D must be n x value_size, inv and coef of length n.");
  std::vector<double> out(n), g(V), d(V);
  for (size_t i = 0; i < n; ++i) {
    bool has_d = true;
    for (size_t v = 0; v < V; ++v) {
      g[v] = G[v * n + i];
      d[v] = D[v * n + i];
      if (std::isnan(d[v]))
        has_d = false;
    }
    out[i] = oob_row_error(i, g.data(), has_d ? d.data() : nullptr, inv[i], nullptr, coef[i]);
  }
  return out;
}
