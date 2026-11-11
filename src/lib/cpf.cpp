
#include "cpf.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include "internal_utils.hpp"
using namespace rpf_utils;
#include <limits>
#include <unordered_set>

// ----------------- rpf subclass for classification -----------------

/**
 * \brief Create a prediction model based on Random Forests for classification data sets.
 */

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// loss moved to lib/losses_*.cpp

// Maps a loss-name string to `loss` and `calcLoss` (verbatim move of the
// former inline if/else chain from the data-fitting constructor).
void ClassificationRPF::set_loss_function(const std::string &loss)
{
  if (loss == "L1") {
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  } else if (loss == "L2") {
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  } else if (loss == "median") {
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  } else if (loss == "logit") {
    this->loss = LossType::logit;
    this->calcLoss = &ClassificationRPF::logit_loss;
  } else if (loss == "logit_2") {
    this->loss = LossType::logit_2;
    this->calcLoss = &ClassificationRPF::logit_loss_2;
  } else if (loss == "logit_3") {
    this->loss = LossType::logit_3;
    this->calcLoss = &ClassificationRPF::logit_loss_3;
  } else if (loss == "logit_4") {
    this->loss = LossType::logit_4;
    this->calcLoss = &ClassificationRPF::logit_loss_4;
  } else if (loss == "exponential") {
    this->loss = LossType::exponential;
    this->calcLoss = &ClassificationRPF::exponential_loss;
  } else if (loss == "exponential_2") {
    this->loss = LossType::exponential_2;
    this->calcLoss = &ClassificationRPF::exponential_loss_2;
  } else if (loss == "exponential_3") {
    this->loss = LossType::exponential_3;
    this->calcLoss = &ClassificationRPF::exponential_loss_3;
  } else {
    throw std::invalid_argument("Unknown loss function '" + loss + "'.");
  }
}

ClassificationRPF::ClassificationRPF(const std::string &loss, const RPFParams &params, double delta, double epsilon)
    : RandomPlantedForest(params), delta(delta), epsilon(epsilon)
{
  set_loss_function(loss);
}

// The fit of one leaf holding the rows of a family that made no split (fit_null_leaf()): the leaf update of
// create_tree_family() for that leaf, with the initial weights (0 for logit, 1 for exponential). For L1 and L2
// that is the leaf's mean of the class indicators, for logit and exponential the link of the class shares,
// truncated to [epsilon, 1 - epsilon] as the update truncates them. Each loss's prediction is then the leaf's
// class shares (truncated for logit and exponential).
std::vector<double>
ClassificationRPF::null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows)
{
  if (loss == LossType::median)
    return calcMedian(Y, rows);
  std::vector<double> M = calcMean(Y, rows);
  if (loss == LossType::L1 || loss == LossType::L2)
    return M;
  // the class shares: logit codes the classes as 0 / 1, exponential as -1 / 1
  if (loss == LossType::exponential || loss == LossType::exponential_2 || loss == LossType::exponential_3)
    for (auto &m : M)
      m = (m + 1) / 2;
  const double M_p = 1 - std::accumulate(M.begin(), M.end(), 0.0); // the reference class
  auto truncated = [this](double m) { return std::min(std::max(epsilon, m), 1 - epsilon); };
  auto above_eps = [this](double m) { return std::max(epsilon, m); };
  std::vector<double> value(value_size, 0);
  switch (this->loss) {
  case LossType::logit:
  case LossType::exponential:
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(truncated(M[p]) / truncated(M_p));
    break;
  case LossType::logit_2:
  case LossType::logit_4:
  case LossType::exponential_2:
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(above_eps(M[p]) / above_eps(1 - M[p]));
    break;
  case LossType::logit_3:
  case LossType::exponential_3: {
    double mean_log = log(above_eps(M_p));
    for (size_t p = 0; p < value_size; ++p)
      mean_log += log(above_eps(M[p]));
    mean_log /= (double)(value_size + 1);
    for (size_t p = 0; p < value_size; ++p)
      value[p] = log(above_eps(M[p])) - mean_log;
    break;
  }
  default:
    break;
  }
  return value;
}

// ----------------- the loss of a tree family -----------------
// The split searches (splits_*.cpp) and create_tree_family() (rpf.cpp) are regression's: the same candidate pool, the
// same sampling of candidates, leaves and thresholds (and therefore the same random draws), the same aging and the
// same sums of Y over each try's two sides. A classification family scores each try with calcLoss, which needs the
// individuals on both sides (I_s / I_b, in the order of the leaf's individuals) and the weights W, and after each
// split it updates the weights (or, for L1 / L2 / median, the working rows' Y) and sets the leaves' values from the
// loss's step. For L2, calcLoss is L2_loss (rpf.cpp): regression's expression in regression's order. With loss "L2"
// on a 0/1 response a classification fit is therefore the regression fit, bit for bit, ties in loss included.

class ClassificationRPF::Loss : public FamilyLoss {
public:
  Loss(ClassificationRPF &forest, std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx)
      : f(forest), Y(samples_Y)
  {
    const size_t vs = f.value_size;
    switch (f.loss) {
    case LossType::exponential:
    case LossType::exponential_2:
    case LossType::exponential_3:
      weights.assign((size_t)f.sample_size, std::vector<double>(vs, 1));
      break;
    default: // logit (the log-odds start at 0), L1 / L2 / median (unused)
      weights.assign((size_t)f.sample_size, std::vector<double>(vs, 0));
    }
    if (f.keep_weights_)
      kept_rows = boot_idx;
    // logit: the candidate-independent terms of the loss (losses_logit.cpp); TLS holds a pointer to them while the
    // family grows
    if (f.loss == LossType::logit) {
      logit_terms.reset(weights, Y, vs);
      tls_logit_terms_ptr = &logit_terms;
    }
  }
  ~Loss() override { tls_logit_terms_ptr = nullptr; }

  double score(Split &trial, const Leaf &leaf, int k, double sp, const std::vector<std::vector<double>> &X) override
  {
    trial.Y = &Y;
    trial.W = &weights;
    trial.I_s.clear();
    trial.I_b.clear();
    for (int ind : leaf.individuals) {
      if (X[ind][k] < sp)
        trial.I_s.push_back(ind);
      else
        trial.I_b.push_back(ind);
    }
    (f.*f.calcLoss)(trial);
    return trial.min_sum;
  }

  void finalize(Split &winner, const std::vector<std::vector<double>> &X) override
  {
    if (std::isinf(winner.min_sum) || winner.leaf_index == nullptr)
      return;
    score(winner, *winner.leaf_index, winner.split_coordinate - 1, winner.split_point, X);
  }

  void apply(
      const Split &split,
      const std::vector<std::vector<double>> &X,
      std::vector<std::vector<double>> &samples_Y,
      std::vector<double> &update_s,
      std::vector<double> &update_b
  ) override;

  void finish(size_t n) override
  {
    // (tests) the weights as the family's growth leaves them (cpf.hpp)
    if (f.keep_weights_ && n < f.kept_weights_.size())
      f.kept_weights_[n] = KeptWeights{std::move(kept_rows), std::move(weights)};
  }

private:
  ClassificationRPF &f;
  std::vector<std::vector<double>> &Y; // the family's working rows' Y
  std::vector<std::vector<double>> weights;
  LogitTerms logit_terms;
  std::vector<int> kept_rows; // (keep_weights_, tests) each working row's row in the data
};

std::unique_ptr<FamilyLoss>
ClassificationRPF::family_loss(std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx)
{
  return std::unique_ptr<FamilyLoss>(new Loss(*this, samples_Y, boot_idx));
}

// The leaves' increments from the loss's step at the split's statistics (M_s, M_b, ...), and the update of the
// weights of the split leaf's individuals (L1 / L2 / median: of their Y)
void ClassificationRPF::Loss::apply(
    const Split &curr_split,
    const std::vector<std::vector<double>> &samples_X,
    std::vector<std::vector<double>> &samples_Y,
    std::vector<double> &update_s,
    std::vector<double> &update_b
)
{
  const size_t value_size = f.value_size;
  const double epsilon = f.epsilon;
  update_s = curr_split.M_s;
  update_b = curr_split.M_b;
  switch (f.loss) {
  case LossType::L1:
  case LossType::L2:
  case LossType::median: {
    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
        samples_Y[individual] -= update_s;
      } else {
        samples_Y[individual] -= update_b;
      }
    }
    break;
  }
  case LossType::logit: {

    std::vector<double> M_s = curr_split.M_s;
    std::vector<double> M_b = curr_split.M_b;

    std::for_each(M_s.begin(), M_s.end(), [epsilon](double &M) { M = std::min(std::max(epsilon, M), 1 - epsilon); });
    std::for_each(M_b.begin(), M_b.end(), [epsilon](double &M) { M = std::min(std::max(epsilon, M), 1 - epsilon); });

    double M_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
    double M_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

    std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
    std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

    for (size_t p = 0; p < value_size; ++p) {
      update_s[p] = log(M_s[p] / M_sp) - W_s_mean[p];
      update_b[p] = log(M_b[p] / M_bp) - W_b_mean[p];
    }

    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
        weights[individual] += update_s;
      } else {
        weights[individual] += update_b;
      }
      if (tls_logit_terms_ptr == &logit_terms)
        logit_terms.update_row((size_t)individual);
    }

    break;
  }
  case LossType::logit_2: {

    std::vector<double> M_s = curr_split.M_s;
    std::vector<double> M_b = curr_split.M_b;

    std::vector<double> M_s2 = curr_split.M_s;
    std::vector<double> M_b2 = curr_split.M_b;

    std::for_each(M_s.begin(), M_s.end(), [epsilon](double &M) { M = std::max(epsilon, M); });
    std::for_each(M_b.begin(), M_b.end(), [epsilon](double &M) { M = std::max(epsilon, M); });

    std::for_each(M_s2.begin(), M_s2.end(), [epsilon](double &M) { M = std::max(epsilon, 1 - M); });
    std::for_each(M_b2.begin(), M_b2.end(), [epsilon](double &M) { M = std::max(epsilon, 1 - M); });

    std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
    std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

    for (size_t p = 0; p < value_size; ++p) {
      update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
      update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
    }

    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
        weights[individual] += update_s;
      } else {
        weights[individual] += update_b;
      }
    }

    break;
  }
  case LossType::logit_3: {

    std::vector<double> M_s = curr_split.M_s;
    std::vector<double> M_b = curr_split.M_b;

    std::for_each(M_s.begin(), M_s.end(), [epsilon](double &M) { M = std::max(epsilon, M); });
    std::for_each(M_b.begin(), M_b.end(), [epsilon](double &M) { M = std::max(epsilon, M); });

    std::for_each(M_s.begin(), M_s.end(), [&](double &M) { M = log(M); });
    std::for_each(M_b.begin(), M_b.end(), [&](double &M) { M = log(M); });

    double M_sp = std::max(epsilon, curr_split.M_sp);
    double M_bp = std::max(epsilon, curr_split.M_bp);

    M_sp = log(M_sp);
    M_bp = log(M_bp);

    double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
    double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);

    std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
    std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

    for (unsigned int p = 0; p < M_s.size(); ++p) {
      update_s[p] = M_s[p] - sum_s - W_s_mean[p];
      update_b[p] = M_b[p] - sum_b - W_b_mean[p];
    }

    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
        weights[individual] += update_s;
      } else {
        weights[individual] += update_b;
      }
    }

    break;
  }
  case LossType::logit_4: {

    std::vector<double> M_s = curr_split.M_s;
    std::vector<double> M_b = curr_split.M_b;

    std::vector<double> M_s2 = curr_split.M_s;
    std::vector<double> M_b2 = curr_split.M_b;

    std::for_each(M_s.begin(), M_s.end(), [epsilon](double &M) { M = std::max(epsilon, M); });
    std::for_each(M_b.begin(), M_b.end(), [epsilon](double &M) { M = std::max(epsilon, M); });

    std::for_each(M_s2.begin(), M_s2.end(), [epsilon](double &M) { M = std::max(epsilon, 1 - M); });
    std::for_each(M_b2.begin(), M_b2.end(), [epsilon](double &M) { M = std::max(epsilon, 1 - M); });

    std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
    std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

    for (size_t p = 0; p < value_size; ++p) {
      update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
      update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
    }

    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
        weights[individual] += update_s;
      } else {
        weights[individual] += update_b;
      }
    }

    break;
  }
  case LossType::exponential: {

    // K >= 3 classes, SAMME (cpf.hpp): the step from the class shares the criterion left in M_s / M_b, clamped
    // to [epsilon, 1 - epsilon]; each row's one weight times exp(-a_c), c its class
    if (value_size > 1) {
      f.samme_step(curr_split.M_s, epsilon, update_s);
      f.samme_step(curr_split.M_b, epsilon, update_b);
      std::vector<double> fac_s(value_size), fac_b(value_size);
      for (size_t k = 0; k < value_size; ++k) {
        fac_s[k] = exp(-update_s[k]);
        fac_b[k] = exp(-update_b[k]);
      }
      for (int individual : curr_split.leaf_index->individuals) {
        const std::vector<double> &y = samples_Y[individual];
        size_t k = 0;
        while (k + 1 < value_size && !(y[k] > 0))
          ++k;
        const double f =
            samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point ? fac_s[k] : fac_b[k];
        for (double &w : weights[individual])
          w *= f;
      }
      break;
    }

    // binary
    std::vector<double> sum_s = curr_split.M_s;
    std::vector<double> sum_b = curr_split.M_b;

    std::for_each(sum_s.begin(), sum_s.end(), [epsilon](double &S) {
      S = std::min(std::max(epsilon, S), 1 - epsilon);
    });
    std::for_each(sum_b.begin(), sum_b.end(), [epsilon](double &S) {
      S = std::min(std::max(epsilon, S), 1 - epsilon);
    });

    double sum_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
    double sum_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

    for (unsigned int p = 0; p < sum_s.size(); ++p) {
      update_s[p] = log(sum_s[p] / sum_sp);
      update_b[p] = log(sum_b[p] / sum_bp);
    }

    for (int individual : curr_split.leaf_index->individuals) {
      for (unsigned int p = 0; p < update_s.size(); ++p) {
        if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
        } else {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
        }
      }
    }

    break;
  }
  case LossType::exponential_2: {

    std::vector<double> sum_s = curr_split.M_s;
    std::vector<double> sum_b = curr_split.M_b;
    std::vector<double> sum_s2 = curr_split.M_s;
    std::vector<double> sum_b2 = curr_split.M_b;

    std::for_each(sum_s.begin(), sum_s.end(), [epsilon](double &S) { S = std::max(epsilon, S); });
    std::for_each(sum_b.begin(), sum_b.end(), [epsilon](double &S) { S = std::max(epsilon, S); });

    std::for_each(sum_s2.begin(), sum_s2.end(), [epsilon](double &S) { S = std::max(epsilon, 1 - S); });
    std::for_each(sum_b2.begin(), sum_b2.end(), [epsilon](double &S) { S = std::max(epsilon, 1 - S); });

    for (size_t p = 0; p < value_size; ++p) {
      update_s[p] = log(sum_s[p] / sum_s2[p]);
      update_b[p] = log(sum_b[p] / sum_b2[p]);
    }

    for (int individual : curr_split.leaf_index->individuals) {
      for (size_t p = 0; p < value_size; ++p) {
        if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
        } else {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
        }
      }
    }

    break;
  }
  case LossType::exponential_3: {

    std::vector<double> sum_s = curr_split.M_s;
    std::vector<double> sum_b = curr_split.M_b;

    std::for_each(sum_s.begin(), sum_s.end(), [epsilon](double &S) { S = std::max(epsilon, S); });
    std::for_each(sum_b.begin(), sum_b.end(), [epsilon](double &S) { S = std::max(epsilon, S); });

    std::for_each(sum_s.begin(), sum_s.end(), [&](double &S) { S = log(S); });
    std::for_each(sum_b.begin(), sum_b.end(), [&](double &S) { S = log(S); });

    double sum_sp = std::max(epsilon, curr_split.M_sp);
    double sum_bp = std::max(epsilon, curr_split.M_bp);

    sum_sp = log(sum_sp);
    sum_bp = log(sum_bp);

    sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
    sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);

    sum_sp = sum_sp / (sum_s.size() + 1);
    sum_bp = sum_bp / (sum_b.size() + 1);

    for (size_t p = 0; p < sum_s.size(); ++p) {
      update_s[p] = sum_s[p] - sum_sp;
      update_b[p] = sum_b[p] - sum_bp;
    }

    for (int individual : curr_split.leaf_index->individuals) {
      for (size_t p = 0; p < update_s.size(); ++p) {
        if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point) {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
        } else {
          weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
        }
      }
    }

    break;
  }
  }
}

// fit forest to new data
void ClassificationRPF::fit()
{
  kept_weights_.assign(keep_weights_ ? (size_t)n_trees : 0, KeptWeights());
  // Use the base class multithreaded trainer with RNG seeding identical to regression
  RandomPlantedForest::fit();
}

/*  retrospectively change parameters of existing class object,
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data',
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values)
{
  if (keys.size() != values.size()) {
    warn("Size of input vectors is not the same. ");
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i) {
    if (keys[i] == "deterministic") {
      this->deterministic = values[i];
    } else if (keys[i] == "nthreads") {
      this->nthreads = values[i];
    } else if (keys[i] == "purify") {
      this->purify_forest = values[i];
    } else if (keys[i] == "n_trees") {
      this->n_trees = values[i];
    } else if (keys[i] == "n_splits") {
      this->n_splits = values[i];
    } else if (keys[i] == "t_try") {
      this->t_try = values[i];
    } else if (keys[i] == "split_try") {
      this->split_try = values[i];
    } else if (keys[i] == "max_interaction") {
      this->max_interaction = values[i];
    } else if (keys[i] == "cv") {
      this->cross_validate = values[i];
    } else if (keys[i] == "loss") {
      if (keys[i] == "L1") {
        this->loss = LossType::L1;
        this->calcLoss = &ClassificationRPF::L1_loss;
      } else if (keys[i] == "L2") {
        this->loss = LossType::L2;
        this->calcLoss = &ClassificationRPF::L2_loss;
      } else if (keys[i] == "median") {
        this->loss = LossType::median;
        this->calcLoss = &ClassificationRPF::median_loss;
      } else if (keys[i] == "logit") {
        this->loss = LossType::logit;
        this->calcLoss = &ClassificationRPF::logit_loss;
      } else if (keys[i] == "exponential") {
        this->loss = LossType::exponential;
        this->calcLoss = &ClassificationRPF::exponential_loss;
      } else {
        warn("Unkown loss function.");
      }
    } else if (keys[i] == "delta") {
      this->delta = values[i];
    } else if (keys[i] == "epsilon") {
      this->epsilon = values[i];
    } else if (keys[i] == "split_decay_rate") {
      this->split_decay_rate_ = values[i];
    } else if (keys[i] == "max_candidates") {
      this->max_candidates_ = static_cast<size_t>(values[i]);
    } else if (keys[i] == "keep_weights") { // internal, for tests (cpf.hpp)
      this->keep_weights_ = values[i] != 0;
    } else {
      warn("Unkown parameter key  '" + keys[i] + "' .");
    }
  }
  this->fit();
}

// The per-row criterion of a classification forest (rpf.hpp; the same as R's crit_rows()). L1 / L2 (and median): the
// squared error summed over the score columns (Y holds their 0/1 coding), H = 2I. logit and exponential: the log-loss
// of the probabilities predict() makes: binary p = 1 / (1 + exp(-g)); multiclass logit the softmax of the K - 1 scores
// and the reference level's score 0; multiclass exponential the softmax of the K scores. H = diag(p) - pp' over the
// free score columns.
double ClassificationRPF::oob_row_error(
    size_t i, const double *g, const double *d, double inv, const double *var, double coef
) const
{
  const size_t V = value_size;
  const std::vector<double> &y = Y[i];
  const bool logit_like =
      loss == LossType::logit || loss == LossType::logit_2 || loss == LossType::logit_3 || loss == LossType::logit_4;
  const bool expo = loss == LossType::exponential || loss == LossType::exponential_2 || loss == LossType::exponential_3;
  if (!logit_like && !expo) {
    double e = 0, q = 0;
    for (size_t v = 0; v < V; ++v) {
      const double r = y[v] - g[v];
      e += r * r;
      if (d)
        q += d[v] * d[v] / inv;
      else if (var)
        q += var[v];
    }
    return e - coef * q;
  }
  // the probabilities of the free score columns, and the log-loss of the row's class
  std::vector<double> p(V);
  double e = 0;
  if (V == 1) {
    // exponential codes the classes as -1 / 1, logit as 0 / 1
    const bool second = expo ? y[0] > 0 : y[0] > 0.5;
    p[0] = 1 / (1 + std::exp(-g[0]));
    // -log p(y): log(1 + exp(-g)) for the second class, log(1 + exp(g)) for the first, stably
    const double z = second ? -g[0] : g[0];
    e = z > 0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
  } else {
    // logit: the reference score 0 joins the softmax; exponential: the K = V scores alone
    double mx = logit_like ? 0.0 : -std::numeric_limits<double>::infinity();
    for (size_t v = 0; v < V; ++v)
      mx = std::max(mx, g[v]);
    double den = logit_like ? std::exp(-mx) : 0.0;
    for (size_t v = 0; v < V; ++v)
      den += std::exp(g[v] - mx);
    const double lse = mx + std::log(den);
    int cls = -1; // the column of the row's class (-1: the reference level of logit)
    for (size_t v = 0; v < V; ++v)
      if (y[v] > (expo ? 0.0 : 0.5))
        cls = (int)v;
    for (size_t v = 0; v < V; ++v)
      p[v] = std::exp(g[v] - lse);
    e = cls < 0 ? lse : lse - g[cls];
  }
  double q = 0;
  if (d) {
    double a = 0, b = 0;
    for (size_t v = 0; v < V; ++v) {
      a += p[v] * d[v] * d[v];
      b += p[v] * d[v];
    }
    q = 0.5 * (a - b * b) / inv;
  } else if (var)
    for (size_t v = 0; v < V; ++v)
      q += 0.5 * p[v] * (1 - p[v]) * var[v];
  return e - coef * q;
}

int ClassificationRPF::score_centring() const
{
  if (value_size < 2)
    return 0;
  if (loss == LossType::logit || loss == LossType::logit_2 || loss == LossType::logit_3 || loss == LossType::logit_4)
    return 1;
  if (loss == LossType::exponential || loss == LossType::exponential_2 || loss == LossType::exponential_3)
    return 2;
  return 0;
}
