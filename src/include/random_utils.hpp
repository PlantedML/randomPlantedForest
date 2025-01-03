#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <chrono>
#include <Rcpp.h>

namespace utils
{

  // Abstract base class for random number generators
  class RNGBackend
  {
  public:
    virtual ~RNGBackend() = default;
    virtual void begin_rng() {} // Called before a sequence of random numbers
    virtual void end_rng() {}   // Called after a sequence of random numbers
    virtual int random_int(int n) = 0;
    virtual double random_double() = 0;
  };

  // Standard C++ random number generator backend
  class StdRNGBackend : public RNGBackend
  {
  private:
    static thread_local std::mt19937 generator;
    static bool seeded;

  public:
    void seed(uint32_t seed)
    {
      generator.seed(seed);
      seeded = true;
    }

    void initialize()
    {
      if (!seeded)
      {
        auto now = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         now.time_since_epoch())
                         .count();
        generator.seed(static_cast<uint32_t>(nanos));
        seeded = true;
      }
    }

    int random_int(int n) override
    {
      initialize();
      std::uniform_int_distribution<int> dist(0, n - 1);
      return dist(generator);
    }

    double random_double() override
    {
      initialize();
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      return dist(generator);
    }
  };

  // R random number generator backend
  class RcppRNGBackend : public RNGBackend
  {
  private:
    std::unique_ptr<Rcpp::RNGScope> rng_scope;

  public:
    void begin_rng() override
    {
      if (!rng_scope)
      {
        rng_scope = std::make_unique<Rcpp::RNGScope>();
      }
    }

    void end_rng() override
    {
      rng_scope.reset();
    }

    int random_int(int n) override
    {
      return static_cast<int>(unif_rand() * n);
    }

    double random_double() override
    {
      return unif_rand();
    }
  };

  class RandomGenerator
  {
  private:
    static RNGBackend *backend;
    static StdRNGBackend std_backend;
    static RcppRNGBackend rcpp_backend;

  public:
    // RAII class to handle RNG state
    class RNGScope
    {
    private:
      RNGBackend *backend;

    public:
      RNGScope() : backend(RandomGenerator::backend)
      {
        backend->begin_rng();
      }
      ~RNGScope()
      {
        backend->end_rng();
      }
    };

    // Switch to using R's RNG
    static void use_r_random()
    {
      backend = &rcpp_backend;
    }

    // Switch to using C++ standard RNG
    static void use_std_random()
    {
      backend = &std_backend;
    }

    // Initialize the standard generator with a seed
    static void seed(uint32_t seed)
    {
      std_backend.seed(seed);
    }

    // Generate random integer in range [0, n)
    static int random_index(int n)
    {
      RNGScope scope;
      return backend->random_int(n);
    }

    // Generate random double in range [0, 1)
    static double random_double()
    {
      RNGScope scope;
      return backend->random_double();
    }

    // Shuffle a range of elements
    template <typename Iter>
    static void shuffle(Iter first, Iter last)
    {
      RNGScope scope;
      auto n = std::distance(first, last);
      for (auto i = n - 1; i > 0; --i)
      {
        std::swap(*(first + i), *(first + backend->random_int(i + 1)));
      }
    }

    // Sample n elements with replacement
    template <typename T>
    static std::vector<T> sample_with_replacement(const std::vector<T> &population, size_t n)
    {
      RNGScope scope;
      std::vector<T> result;
      result.reserve(n);
      for (size_t i = 0; i < n; ++i)
      {
        result.push_back(population[backend->random_int(population.size())]);
      }
      return result;
    }
  };

} // namespace utils

#endif // RANDOM_UTILS_H
