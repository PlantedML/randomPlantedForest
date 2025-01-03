#include "helper.hpp"
#include "random_utils.hpp"
#include <vector>
#include <iostream>
#include <set>

namespace utils
{

  // Helper function to generate random number using our RandomGenerator
  int random_index(const int n) { return RandomGenerator::random_index(n); }

  //  ----------------- functions for converting Cpp types -----------------

  /**
   * \brief Convert the std container set of type int into a std::vector<int>.
   *
   * \param v the set that is converted.
   */
  std::vector<int> from_std_set(std::set<int> v)
  {
    return std::vector<int>(v.begin(), v.end());
  }

  /**
   * \brief Convert the std container vector of type int into a std::vector<int>.
   *
   * \param v the vector that is converted.
   */
  std::vector<int> from_std_vec(std::vector<int> v)
  {
    return v;
  }

  /**
   * \brief Convert the std container vector of type double into a std::vector<double>.
   *
   * \param v the vector that is converted.
   */
  std::vector<double> from_std_vec(std::vector<double> v)
  {
    return v;
  }

  /**
   * \brief Convert the nested std container vector containing a vector itself
   * of type double into a std::vector<std::vector<double>>.
   *
   * \param v the vector of vectors that is converted.
   */
  std::vector<std::vector<double>> from_std_vec(std::vector<std::vector<double>> v)
  {
    return v;
  }

  std::vector<int> to_std_vec(std::vector<int> rv)
  {
    return rv;
  }

  std::vector<double> to_std_vec(std::vector<double> rv)
  {
    return rv;
  }

  std::vector<std::vector<double>> to_std_vec(std::vector<std::vector<double>> rv)
  {
    return rv;
  }

  std::set<int> to_std_set(std::vector<int> rv)
  {
    return std::set<int>(rv.begin(), rv.end());
  }

  std::set<int> to_std_set(std::vector<double> rv)
  {
    return std::set<int>(rv.begin(), rv.end());
  }

}
