// The panel form of purified components (coarsen.hpp): for prune_components(max_panels = P), every kept component
// of order >= 3 is shown in at most P heat maps. Candidates (C, q) in lexicographic order, equal-mass cuts of the
// training rows, the choice by the variance kept. The components' rectangles come from passes over the families
// restricted to them (merged_components()), their values at the training rows from the caller (purify_columns()), so
// no purified store is needed.
#include "rpf.hpp"
#include "subset_count.hpp"
#include <atomic>
#include <chrono>
#include <exception>
#include <memory>
#include <set>
#include <stdexcept>
#include <thread>

using rpf_coarse::Axis;
using rpf_coarse::Boxes;
using rpf_coarse::Built;
using rpf_coarse::Candidate;
using rpf_coarse::Selected;

namespace {
// The training rows per feature: values, sorted values, forest cells
struct Train {
  size_t n = 0;
  int p = 0;
  std::vector<std::vector<double>> x, xs;
  std::vector<std::vector<int32_t>> cells;
  std::vector<int32_t> ncell;
};

Train make_train(const std::vector<std::vector<double>> &X, const rpf_rect::PurifiedForest &pf, int p)
{
  Train T;
  T.n = X.size();
  T.p = p;
  T.x.assign(p, std::vector<double>(T.n));
  T.cells.assign(p, std::vector<int32_t>(T.n));
  T.ncell.assign(p, 0);
  for (int d = 0; d < p; ++d) {
    for (size_t i = 0; i < T.n; ++i) {
      T.x[d][i] = X[i][d];
      T.cells[d][i] = pf.cell(d, X[i][d]);
    }
    T.ncell[d] = (int32_t)pf.cuts[d].size() + 1;
  }
  T.xs = T.x;
  for (auto &v : T.xs)
    std::sort(v.begin(), v.end());
  return T;
}

std::vector<int32_t> ncell_of(const Train &T, const std::vector<int> &S)
{
  std::vector<int32_t> nc;
  for (int d : S)
    nc.push_back(T.ncell[d - 1]);
  return nc;
}

std::vector<Axis> axes_of(
    const Train &T,
    const std::vector<int> &S,
    const std::vector<int> &cpos,
    const std::vector<std::vector<double>> &cuts
)
{
  std::vector<Axis> A;
  for (size_t i = 0; i < cpos.size(); ++i) {
    const int f = S[cpos[i]] - 1;
    A.push_back(rpf_coarse::make_axis(cpos[i], f, cuts[i], T.cells[f], T.x[f], T.ncell[f]));
  }
  return A;
}

// The coarse cell of training row i (row-major over the axes)
template <class V> inline int coarse_cell(const V &A, const Train &T, size_t i)
{
  int cc = 0;
  for (const auto &a0 : A) {
    const Axis &a = rpf_coarse::ax(a0);
    cc = cc * a.q + a.interval(T.x[a.feat][i]);
  }
  return cc;
}

// Values of a built coarsened component at the training rows (sums over the families, not divided)
template <class V>
std::vector<double> eval_train(const Built &b, const V &A, const std::vector<int> &dfeat, const Train &T)
{
  std::vector<double> out(T.n), buf(b.nc);
  std::vector<int32_t> dc(std::max<size_t>(1, dfeat.size()));
  for (size_t i = 0; i < T.n; ++i) {
    for (size_t q = 0; q < dfeat.size(); ++q)
      dc[q] = T.cells[dfeat[q]][i];
    out[i] = b.at(dc.data(), coarse_cell(A, T, i), buf.data());
  }
  return out;
}

// Grid limit (doubles) of the difference-array evaluation of a candidate; beyond it the k-d index
const double DIFF_MAX = 8e6;

// fval: the forest's f_S at the training rows (T.n values); pf: the forest's cells
void select_one(
    Selected &R,
    const rpf_rect::PurifiedForest &pf,
    const Train &T,
    double fp,
    const rpf_rect::Accumulator &raw,
    int n_trees,
    int P,
    const double *fval
)
{
  auto t0 = std::chrono::steady_clock::now();
  const std::vector<int> &S = R.S;
  const int k = (int)S.size();
  R.forest_panels = fp;

  bool nonzero = false;
  for (size_t i = 0; i < T.n && !nonzero; ++i)
    nonzero = fval[i] != 0.0;
  if (!nonzero) {
    R.form = 0;
    return;
  }
  std::vector<int32_t> ncell = ncell_of(T, S);
  R.var_forest = rpf_coarse::var_pop(std::vector<double>(fval, fval + T.n));
  if (fp <= (double)P) {
    R.form = 1;
    R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return;
  }
  R.form = 2;
  // the component's rectangles as the families wrote them, merged
  auto rw = raw.comps.find(S);
  if (rw == raw.comps.end())
    throw std::logic_error("coarsen_select(): a nonzero component without rectangles.");
  const Boxes B =
      rpf_coarse::boxes_from_groups(k, rw->second.first.empty() ? 0.0 : rw->second.first[0], rw->second.second, ncell);
  const size_t NR = B.size();
  std::vector<std::vector<int>> qv;
  rpf_coarse::q_vectors(k - 2, P, qv);
  // the coarse sets C: (k - 2)-subsets of the positions, lexicographic (itertools.combinations)
  std::vector<std::vector<int>> combos;
  {
    std::vector<int> c(k - 2);
    std::function<void(int, int)> rec = [&](int i, int from) {
      if (i == k - 2) {
        combos.push_back(c);
        return;
      }
      for (int v = from; v < k; ++v) {
        c[i] = v;
        rec(i + 1, v + 1);
      }
    };
    rec(0, 0);
  }
  // equal-mass cuts per (position, q)
  std::map<std::pair<int, int>, std::vector<double>> cut_cache;
  auto cuts_of = [&](int pos, int q) -> const std::vector<double> & {
    auto key = std::make_pair(pos, q);
    auto c = cut_cache.find(key);
    if (c != cut_cache.end())
      return c->second;
    std::vector<double> t = q > 1 ? rpf_coarse::eq_cuts(T.xs[S[pos] - 1], q) : std::vector<double>();
    return cut_cache.emplace(key, t).first->second;
  };

  std::map<std::pair<int, int>, Axis> axis_cache; // the coarse axis of (position, q)
  double best = -1;
  std::vector<double> gbuf(T.n);
  std::vector<double> W;
  for (size_t ci = 0; ci < combos.size(); ++ci) {
    const std::vector<int> &cpos = combos[ci];
    std::vector<int> dpos;
    for (int q = 0; q < k; ++q)
      if (std::find(cpos.begin(), cpos.end(), q) == cpos.end())
        dpos.push_back(q);
    const int a = dpos[0], b = dpos[1];
    const int fa = S[a] - 1, fb = S[b] - 1;
    // compressed grid of D's boundaries (0 always in it)
    std::vector<int32_t> Ba{0}, Bb{0};
    for (size_t r = 0; r < NR; ++r) {
      Ba.push_back(B.lo[r * k + a]);
      Ba.push_back(B.hi[r * k + a]);
      Bb.push_back(B.lo[r * k + b]);
      Bb.push_back(B.hi[r * k + b]);
    }
    std::sort(Ba.begin(), Ba.end());
    Ba.erase(std::unique(Ba.begin(), Ba.end()), Ba.end());
    std::sort(Bb.begin(), Bb.end());
    Bb.erase(std::unique(Bb.begin(), Bb.end()), Bb.end());
    const size_t Na = Ba.size(), Nb = Bb.size();
    std::vector<int32_t> ia0(NR), ia1(NR), ib0(NR), ib1(NR);
    auto ix = [](const std::vector<int32_t> &v, int32_t x) {
      return (int32_t)(std::lower_bound(v.begin(), v.end(), x) - v.begin());
    };
    for (size_t r = 0; r < NR; ++r) {
      ia0[r] = ix(Ba, B.lo[r * k + a]);
      ia1[r] = ix(Ba, B.hi[r * k + a]);
      ib0[r] = ix(Bb, B.lo[r * k + b]);
      ib1[r] = ix(Bb, B.hi[r * k + b]);
    }
    std::vector<size_t> rowidx(T.n);
    for (size_t i = 0; i < T.n; ++i) {
      size_t ra = (size_t)(std::upper_bound(Ba.begin(), Ba.end(), T.cells[fa][i]) - Ba.begin()) - 1;
      size_t rb = (size_t)(std::upper_bound(Bb.begin(), Bb.end(), T.cells[fb][i]) - Bb.begin()) - 1;
      rowidx[i] = ra * Nb + rb;
    }
    for (const auto &qs : qv) {
      std::vector<std::vector<double>> cuts;
      for (size_t i = 0; i < cpos.size(); ++i)
        cuts.push_back(cuts_of(cpos[i], qs[i]));
      std::vector<const Axis *> A;
      for (size_t i = 0; i < cpos.size(); ++i) {
        auto key = std::make_pair(cpos[i], qs[i]);
        auto c = axis_cache.find(key);
        if (c == axis_cache.end()) {
          const int f = S[cpos[i]] - 1;
          c = axis_cache.emplace(key, rpf_coarse::make_axis(cpos[i], f, cuts[i], T.cells[f], T.x[f], T.ncell[f])).first;
        }
        A.push_back(&c->second);
      }
      const int nc = rpf_coarse::n_cells_of(A);
      W.assign(nc, 0.0);
      if ((double)Na * (double)Nb * (double)nc <= DIFF_MAX) {
        R.n_diff++;
        std::vector<double> G(Na * Nb * (size_t)nc, 0.0);
        for (size_t r = 0; r < NR; ++r) {
          rpf_coarse::coarse_weights(B, r, A, W.data(), nc);
          double *g00 = G.data() + ((size_t)ia0[r] * Nb + ib0[r]) * nc;
          double *g10 = G.data() + ((size_t)ia1[r] * Nb + ib0[r]) * nc;
          double *g01 = G.data() + ((size_t)ia0[r] * Nb + ib1[r]) * nc;
          double *g11 = G.data() + ((size_t)ia1[r] * Nb + ib1[r]) * nc;
          const bool ha = (size_t)ia1[r] < Na, hb = (size_t)ib1[r] < Nb;
          for (int c = 0; c < nc; ++c) {
            g00[c] += W[c];
            if (ha)
              g10[c] -= W[c];
            if (hb)
              g01[c] -= W[c];
            if (ha && hb)
              g11[c] += W[c];
          }
        }
        const size_t rowlen = Nb * (size_t)nc;
        for (size_t i = 1; i < Na; ++i) {
          double *cur = G.data() + i * rowlen, *prev = cur - rowlen;
          for (size_t j = 0; j < rowlen; ++j)
            cur[j] += prev[j];
        }
        for (size_t i = 0; i < Na; ++i)
          for (size_t j = 1; j < Nb; ++j) {
            double *cur = G.data() + (i * Nb + j) * nc, *prev = cur - nc;
            for (int c = 0; c < nc; ++c)
              cur[c] += prev[c];
          }
        for (size_t i = 0; i < T.n; ++i)
          gbuf[i] = G[rowidx[i] * nc + coarse_cell(A, T, i)] / n_trees;
      } else {
        R.n_kd++;
        rpf_rect::RectList rl = rpf_coarse::coarse_rects(B, A, dpos);
        Built bt;
        bt.build(rl, {T.ncell[fa], T.ncell[fb]});
        gbuf = eval_train(bt, A, {fa, fb}, T);
        for (auto &v : gbuf)
          v /= n_trees;
      }
      const double vk = R.var_forest > 0 ? rpf_coarse::var_pop(gbuf) / R.var_forest : 0.0;
      R.cands.push_back(Candidate{(int)ci, qs, nc, vk});
      if (best < 0 || vk > best + 1e-15) {
        best = vk;
        R.cpos = cpos;
        R.dpos = dpos;
        R.q = qs;
        R.cuts = cuts;
        R.var_kept = vk;
      }
    }
  }
  // the chosen candidate's stored form and its training-row values (as predict computes them)
  std::vector<Axis> A = axes_of(T, S, R.cpos, R.cuts);
  rpf_coarse::Form &F = R.stored;
  F.rects = rpf_coarse::coarse_rects(B, A, R.dpos);
  F.S = S;
  F.n_trees = n_trees;
  Built bt;
  std::vector<int> dfeat;
  std::vector<int32_t> nd;
  for (int q : R.dpos) {
    dfeat.push_back(S[q] - 1);
    nd.push_back(T.ncell[S[q] - 1]);
    F.dfeat.push_back(S[q]);
    F.dcuts.push_back(pf.cuts[S[q] - 1]);
  }
  for (size_t i = 0; i < R.cpos.size(); ++i) {
    F.cfeat.push_back(S[R.cpos[i]]);
    F.ccuts.push_back(R.cuts[i]);
  }
  bt.build(F.rects, nd);
  R.g = eval_train(bt, A, dfeat, T);
  for (auto &v : R.g)
    v /= n_trees;
  R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Run fn(i) for i in 0..n on up to nthreads threads (dynamic); an exception in fn is rethrown on the calling thread
void parallel_for(size_t n, int nthreads, const std::function<void(size_t)> &fn)
{
  const size_t nt = std::max<size_t>(1, std::min<size_t>((size_t)std::max(1, nthreads), n));
  if (nt <= 1) {
    for (size_t i = 0; i < n; ++i)
      fn(i);
    return;
  }
  std::atomic<size_t> next(0);
  std::vector<std::exception_ptr> errors(nt);
  std::vector<std::thread> th;
  for (size_t t = 0; t < nt; ++t)
    // an exception escaping a std::thread would terminate R
    th.emplace_back([&, t]() {
      try {
        for (size_t i = next++; i < n; i = next++)
          fn(i);
      } catch (...) {
        errors[t] = std::current_exception();
        next = n;
      }
    });
  for (auto &x : th)
    x.join();
  for (auto &e : errors)
    if (e)
      std::rethrow_exception(e);
}

// A set as given (1-based features), sorted
std::vector<int> sorted_set(std::vector<int> S)
{
  std::sort(S.begin(), S.end());
  return S;
}

void check_coarsen(size_t value_size, bool has_data)
{
  if (value_size != 1)
    throw std::invalid_argument("Coarsening supports regression forests only.");
  if (!has_data)
    throw std::invalid_argument("Coarsening needs the forest's training data.");
}
} // namespace

rpf_rect::Accumulator RandomPlantedForest::merged_components(
    const std::set<std::vector<int>> &only, int nthreads, const std::vector<std::vector<double>> &cuts
) const
{
  rpf_rect::Accumulator total;
  if (only.empty())
    return total;
  rpf_kd::KDTree kdt;
  rpf_cnt::SubsetCounter counter(&X, feature_size);
  const unsigned nt = (unsigned)std::max(1, nthreads);
  std::vector<std::pair<int, rpf_rect::Accumulator>> stack;
  std::vector<rpf_rect::Accumulator> batch(nt);
  run_families(
      tree_families.size(),
      nt,
      [&](size_t f) {
        rpf_rect::Accumulator &acc = batch[f % nt];
        acc = rpf_rect::Accumulator();
        purify_fast_exact_family(tree_families[f], feature_size, kdt, cuts, acc, &counter, &only);
      },
      [&](size_t, size_t nb) {
        for (size_t i = 0; i < nb; ++i) {
          stack.emplace_back(0, std::move(batch[i]));
          batch[i] = rpf_rect::Accumulator();
          while (stack.size() >= 2 && stack[stack.size() - 2].first == stack.back().first) {
            stack[stack.size() - 2].second.absorb(stack.back().second);
            stack[stack.size() - 2].first += 1;
            stack.pop_back();
          }
        }
      }
  );
  for (auto &st : stack)
    total.absorb(st.second);
  return total;
}

std::vector<Selected> RandomPlantedForest::coarsen_select(
    const std::vector<std::vector<int>> &sets, int P, int nthreads, const double *fvals, int n_rows, int n_cols
)
{
  check_coarsen(value_size, !X.empty());
  if (P < 1)
    throw std::invalid_argument("P must be at least 1.");
  const size_t K = sets.size();
  if (n_rows != sample_size || n_cols != (int)K || (K > 0 && fvals == nullptr))
    throw std::invalid_argument("`fvals` must have one column per set and one row per training row.");
  rpf_rect::PurifiedForest pf; // the cells only
  pf.cuts = forest_cuts();
  const Train T = make_train(X, pf, feature_size);
  std::vector<Selected> res(K);
  for (size_t i = 0; i < K; ++i) {
    res[i].S = sorted_set(sets[i]);
    if (res[i].S.size() < 3 || res[i].S.front() < 1 || res[i].S.back() > feature_size)
      throw std::invalid_argument("Each set must have at least 3 features in 1..p.");
  }
  // the panels the forest's f_S needs: the product of the |S| - 2 smallest interval counts of its interior cut points
  // over the leaves of all trees >= S, all families. If at most P (f_S kept as it is), R.cpos and R.cuts are its
  // panels: the |S| - 2 features with the fewest interior cut points (the first of ties), in the order of S, and
  // those cut points.
  auto panels_of = [&](Selected &R) -> double {
    const std::vector<int> &S = R.S;
    const int k = (int)S.size();
    std::vector<std::vector<double>> own(k);
    for (const auto &fam : tree_families)
      for (const auto &tree : fam) {
        bool sup = true;
        for (int d : S)
          sup = sup && tree.first.count(d);
        if (!sup)
          continue;
        for (const auto &leaf : tree.second->leaves)
          for (int q = 0; q < k; ++q) {
            own[q].push_back(leaf.intervals[S[q] - 1].first);
            own[q].push_back(leaf.intervals[S[q] - 1].second);
          }
      }
    std::vector<double> qn;
    for (int q = 0; q < k; ++q) {
      auto &v = own[q];
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      const auto &xs = T.xs[S[q] - 1];
      size_t c = 0;
      for (double u : v)
        c += (u > xs.front() && u <= xs.back());
      qn.push_back((double)c + 1);
    }
    const std::vector<double> counts = qn;
    std::sort(qn.begin(), qn.end());
    double fp = 1;
    for (int q = 0; q < k - 2; ++q)
      fp *= qn[q];
    if (fp <= (double)P) {
      std::vector<int> pos;
      for (int q = 0; q < k; ++q)
        pos.push_back(q);
      std::stable_sort(pos.begin(), pos.end(), [&](int a, int b) { return counts[a] < counts[b]; });
      pos.resize(k - 2);
      std::sort(pos.begin(), pos.end());
      R.cpos = pos;
      for (int q : pos) {
        const auto &xs = T.xs[S[q] - 1];
        std::vector<double> t;
        for (double u : own[q])
          if (u > xs.front() && u <= xs.back())
            t.push_back(u);
        R.cuts.push_back(t);
      }
    }
    return fp;
  };
  std::vector<double> fps(K);
  parallel_for(K, nthreads, [&](size_t i) { fps[i] = panels_of(res[i]); });
  std::set<std::vector<int>> need;
  for (size_t i = 0; i < K; ++i)
    if (fps[i] > (double)P)
      need.insert(res[i].S);
  const rpf_rect::Accumulator raw = merged_components(need, nthreads, pf.cuts);
  parallel_for(K, nthreads, [&](size_t i) {
    select_one(res[i], pf, T, fps[i], raw, n_trees, P, fvals + i * (size_t)sample_size);
  });
  return res;
}

namespace {
// A coarsening on the training rows: the coarse axes and the positions in S
struct TrainSpec {
  std::vector<int> S, cpos, dpos, dfeat;
  std::vector<Axis> A;
  std::vector<int32_t> ncell, nd;
};

TrainSpec make_spec(const rpf_coarse::Spec &sp, const Train &T, int p)
{
  TrainSpec s;
  s.S = sorted_set(sp.S);
  if (s.S.empty() || s.S.front() < 1 || s.S.back() > p)
    throw std::invalid_argument("S must be a set of features in 1..p.");
  const std::vector<int> &C = sp.C;
  if (sp.cuts.size() != C.size())
    throw std::invalid_argument("`cuts` must have one vector per feature of C.");
  for (int q = 0; q < (int)s.S.size(); ++q) {
    auto f = std::find(C.begin(), C.end(), s.S[q]);
    if (f == C.end()) {
      s.dpos.push_back(q);
      s.dfeat.push_back(s.S[q] - 1);
      s.nd.push_back(T.ncell[s.S[q] - 1]);
    }
  }
  // the coarse axes in the order of C as given
  for (size_t i = 0; i < C.size(); ++i) {
    auto f = std::find(s.S.begin(), s.S.end(), C[i]);
    if (f == s.S.end())
      throw std::invalid_argument("C must be a subset of S.");
    s.cpos.push_back((int)(f - s.S.begin()));
  }
  s.A = axes_of(T, s.S, s.cpos, sp.cuts);
  for (const auto &a : s.A)
    for (double c : a.cntI)
      if (!(c > 0))
        throw std::invalid_argument("A coarse interval holds no training row.");
  s.ncell = ncell_of(T, s.S);
  return s;
}
} // namespace

rpf_coarse::OobSums RandomPlantedForest::coarsen_oob(const std::vector<rpf_coarse::Spec> &specs, int nthreads)
{
  check_coarsen(value_size, !X.empty());
  bool ok = inbag_.size() == tree_families.size() && !tree_families.empty();
  for (const auto &inb : inbag_)
    ok = ok && inb.size() == (size_t)sample_size;
  if (!ok)
    throw std::invalid_argument("The forest has no in-bag counts.");
  rpf_rect::PurifiedForest pf; // the cells only
  pf.cuts = forest_cuts();
  const Train T = make_train(X, pf, feature_size);
  const size_t K = specs.size(), n = T.n;
  std::vector<TrainSpec> sp;
  std::set<std::vector<int>> only;
  for (size_t k = 0; k < K; ++k) {
    sp.push_back(make_spec(specs[k], T, feature_size));
    only.insert(sp.back().S);
  }
  rpf_coarse::OobSums out;
  out.sum0.assign(n * K, 0.0);
  out.sum1.assign(n * K, 0.0);
  if (K == 0)
    return out;
  rpf_kd::KDTree kdt;
  rpf_cnt::SubsetCounter counter(&X, feature_size);
  const unsigned nt = (unsigned)std::max(1, nthreads);
  // per slot: the family's OOB rows and its values there, per spec
  std::vector<std::vector<int>> rows(nt);
  std::vector<std::vector<std::vector<double>>> vals(nt);
  run_families(
      tree_families.size(),
      nt,
      [&](size_t f) {
        const size_t slot = f % nt;
        rpf_rect::Accumulator acc;
        purify_fast_exact_family(tree_families[f], feature_size, kdt, pf.cuts, acc, &counter, &only);
        std::vector<int> &rw = rows[slot];
        rw.clear();
        for (size_t i = 0; i < n; ++i)
          if (inbag_[f][i] == 0)
            rw.push_back((int)i);
        vals[slot].assign(K, std::vector<double>());
        for (size_t k = 0; k < K; ++k) {
          const TrainSpec &s = sp[k];
          std::vector<double> &v = vals[slot][k];
          v.assign(rw.size(), 0.0);
          auto it = acc.comps.find(s.S);
          if (it == acc.comps.end())
            continue;
          Boxes Bx = rpf_coarse::boxes_from_groups(
              (int)s.S.size(), it->second.first.empty() ? 0.0 : it->second.first[0], it->second.second, s.ncell
          );
          rpf_rect::RectList rl = rpf_coarse::coarse_rects(Bx, s.A, s.dpos);
          Built bt;
          bt.build(rl, s.nd);
          std::vector<double> buf(bt.nc);
          std::vector<int32_t> dc(std::max<size_t>(1, s.dfeat.size()));
          for (size_t r = 0; r < rw.size(); ++r) {
            const size_t i = (size_t)rw[r];
            for (size_t q = 0; q < s.dfeat.size(); ++q)
              dc[q] = T.cells[s.dfeat[q]][i];
            v[r] = bt.at(dc.data(), coarse_cell(s.A, T, i), buf.data());
          }
        }
      },
      [&](size_t start, size_t batch) {
        // in family order, on the calling thread: the sums do not depend on the number of threads
        for (size_t i = 0; i < batch; ++i) {
          std::vector<double> &sum = ((start + i) % 2 == 0) ? out.sum0 : out.sum1;
          for (size_t k = 0; k < K; ++k)
            for (size_t r = 0; r < rows[i].size(); ++r)
              sum[k * n + (size_t)rows[i][r]] += vals[i][k][r];
        }
      }
  );
  return out;
}

std::vector<double> RandomPlantedForest::coarsen_eval(const rpf_coarse::Spec &spec, const double *Xn, int n, int p)
{
  if (!(purified && rect_forest_))
    throw std::invalid_argument("The forest is not purified (mode 2).");
  check_coarsen(value_size, !X.empty());
  if (p != feature_size)
    throw std::invalid_argument("X has the wrong number of columns.");
  const rpf_rect::PurifiedForest &pf = *rect_forest_;
  const Train T = make_train(X, pf, feature_size);
  TrainSpec s = make_spec(spec, T, feature_size);
  std::vector<double> out((size_t)n, 0.0);
  auto it = pf.comps.find(s.S);
  if (it == pf.comps.end())
    return out;
  Boxes Bx = rpf_coarse::boxes_from_component(it->second, s.ncell);
  rpf_rect::RectList rl = rpf_coarse::coarse_rects(Bx, s.A, s.dpos);
  Built bt;
  bt.build(rl, s.nd);
  std::vector<double> buf(bt.nc);
  std::vector<int32_t> dc(std::max<size_t>(1, s.dfeat.size()));
  for (size_t i = 0; i < (size_t)n; ++i) {
    for (size_t q = 0; q < s.dfeat.size(); ++q)
      dc[q] = pf.cell(s.dfeat[q], Xn[(size_t)s.dfeat[q] * n + i]);
    int cc = 0;
    for (const auto &a : s.A)
      cc = cc * a.q + a.interval(Xn[(size_t)a.feat * n + i]);
    out[i] = bt.at(dc.data(), cc, buf.data()) / n_trees;
  }
  return out;
}

rpf_coarse::Model::Model(const Form &f)
{
  const int m = (int)f.dfeat.size(), nc = f.rects.vs;
  const size_t R = f.rects.size();
  if (f.rects.m != m || (int)f.dcuts.size() != m || f.ccuts.size() != f.cfeat.size() || f.n_trees < 1)
    throw std::invalid_argument("Corrupt coarsened component.");
  n_trees = f.n_trees;
  std::vector<int32_t> nd(m);
  for (int q = 0; q < m; ++q) {
    if (f.dfeat[q] < 1)
      throw std::invalid_argument("Corrupt coarsened component.");
    dfeat.push_back(f.dfeat[q] - 1);
    dcuts.push_back(f.dcuts[q]);
    nd[q] = (int32_t)dcuts[q].size() + 1;
    for (size_t r = 0; r < R; ++r)
      if (f.rects.lo(r)[q] < 0 || f.rects.lo(r)[q] >= f.rects.hi(r)[q] || f.rects.hi(r)[q] > nd[q])
        throw std::invalid_argument("Corrupt coarsened component.");
  }
  int ncc = 1;
  for (size_t a = 0; a < f.cfeat.size(); ++a) {
    if (f.cfeat[a] < 1)
      throw std::invalid_argument("Corrupt coarsened component.");
    cfeat.push_back(f.cfeat[a] - 1);
    ccuts.push_back(f.ccuts[a]);
    ncc *= (int)ccuts[a].size() + 1;
  }
  if (ncc != nc)
    throw std::invalid_argument("Corrupt coarsened component.");
  bt.build(f.rects, nd);
}

std::vector<double> rpf_coarse::Model::eval(const double *X, int n, int p) const
{
  const int m = (int)dfeat.size();
  for (int f : dfeat)
    if (f >= p)
      throw std::invalid_argument("X has too few columns.");
  for (int f : cfeat)
    if (f >= p)
      throw std::invalid_argument("X has too few columns.");
  std::vector<double> out((size_t)n);
  std::vector<double> buf(bt.nc);
  std::vector<int32_t> dc(std::max(1, m));
  for (size_t i = 0; i < (size_t)n; ++i) {
    for (int q = 0; q < m; ++q) {
      const auto &c = dcuts[q];
      dc[q] = (int32_t)(std::upper_bound(c.begin(), c.end(), X[(size_t)dfeat[q] * n + i]) - c.begin());
    }
    int cc = 0;
    for (size_t a = 0; a < cfeat.size(); ++a) {
      const auto &t = ccuts[a];
      cc = cc * ((int)t.size() + 1) +
           (int)(std::upper_bound(t.begin(), t.end(), X[(size_t)cfeat[a] * n + i]) - t.begin());
    }
    out[i] = bt.at(dc.data(), cc, buf.data()) / n_trees;
  }
  return out;
}
