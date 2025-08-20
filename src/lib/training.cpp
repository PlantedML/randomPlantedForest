// Training orchestration split out from rpf.cpp. Builds tree families,
// manages bootstrapping and threading, and handles optional purification.
#include "rpf.hpp"
#include "internal_utils.hpp"

using namespace rpf_utils;

void RandomPlantedForest::fit()
{
  std::vector<int> initial_individuals(sample_size);
  std::iota(initial_individuals.begin(), initial_individuals.end(), 0);

  std::vector<Interval> initial_intervals(feature_size);
  for (int i = 0; i < feature_size; ++i)
    initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};

  Leaf initial_leaf;
  {
    initial_leaf.value = std::vector<double>(value_size, 0);
    initial_leaf.individuals = initial_individuals;
    initial_leaf.intervals = initial_intervals;
  }
  std::vector<Leaf> initial_leaves{initial_leaf};

  this->tree_families = std::vector<TreeFamily>(n_trees);

  // Generate per-tree seeds from R's RNG to ensure reproducibility across runs
  // when the user sets the R seed. These seeds will be used regardless of
  // threading mode.
  tree_seeds_.assign((size_t)std::max(0, n_trees), 0ULL);
  for (int i = 0; i < n_trees; ++i) {
    // Two 32-bit chunks composed into a 64-bit seed using R's RNG
    unsigned long long hi = static_cast<unsigned long long>(R::runif(0.0, 4294967296.0));
    unsigned long long lo = static_cast<unsigned long long>(R::runif(0.0, 4294967296.0));
    tree_seeds_[(size_t)i] = (hi << 32) ^ lo ^ static_cast<unsigned long long>(i);
  }

  unsigned int threads_to_use = static_cast<unsigned int>(nthreads);
  if (threads_to_use == 0) threads_to_use = 1;
  if (threads_to_use > 1)
  {
    if (threads_to_use > std::thread::hardware_concurrency())
    {
      Rcout << "Requested " << threads_to_use << " threads but only " << std::thread::hardware_concurrency() << " available" << std::endl;
    }
    for (int start = 0; start < n_trees; start += (int)threads_to_use)
    {
      int batch = std::min<int>((int)threads_to_use, n_trees - start);
      if (batch <= 0) break;
      std::vector<std::thread> threads((size_t)batch);
      for (int i = 0; i < batch; ++i)
      {
        int tree_index = start + i;
        threads[(size_t)i] = std::thread([this, &initial_leaves](int tree_index_inner){
          std::mt19937_64 rng_local;
          std::mt19937_64* prev_ptr = rpf_utils::swap_tls_rng(nullptr);
          if (!tree_seeds_.empty() && (size_t)tree_index_inner < tree_seeds_.size()) {
            rng_local.seed(tree_seeds_[(size_t)tree_index_inner]);
          } else {
            rng_local.seed(88172645463393265ULL ^ (unsigned long long)tree_index_inner);
          }
          rpf_utils::swap_tls_rng(&rng_local);
          this->create_tree_family(initial_leaves, (size_t)tree_index_inner);
          rpf_utils::swap_tls_rng(prev_ptr);
        }, tree_index);
      }
      for (auto &th : threads)
      {
        if (th.joinable()) th.join();
      }
    }
  }
  else
  {
    // Single-threaded: still drive randomness from per-tree seeds
    std::mt19937_64 rng_local;
    std::mt19937_64* prev_ptr = rpf_utils::swap_tls_rng(nullptr);
    for (int n = 0; n < n_trees; ++n)
    {
      if (!tree_seeds_.empty() && (size_t)n < tree_seeds_.size()) {
        rng_local.seed(tree_seeds_[(size_t)n]);
      } else {
        rng_local.seed(88172645463393265ULL ^ (unsigned long long)n);
      }
      rpf_utils::swap_tls_rng(&rng_local);
      create_tree_family(initial_leaves, n);
    }
    rpf_utils::swap_tls_rng(prev_ptr);
  }

  if (purify_forest)
  {
    this->purify_3();
  }
  else
  {
    purified = false;
  }
}


