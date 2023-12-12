#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <Rcpp.h>
#include <thread>
#include <assert.h>

#ifndef UTILS_H
#define UTILS_H

namespace utils
{
  // ----------------- custom data types -----------------

  struct setComp
  {
    bool operator()(const std::set<int> &a, const std::set<int> &b) const
    {
      if (a == b)
        return false; // what if same?
      if (a.size() == b.size())
      {
        std::set<int>::iterator it2 = b.begin();
        for (std::set<int>::iterator it1 = a.begin(); it1 != a.end(); it1++)
        {
          if (*it1 < *it2)
          {
            return true;
          }
          else if (*it1 > *it2)
          {
            return false;
          }
          it2++;
        }
      }
      return (a.size() < b.size());
    }
  };

  template <typename T>
  struct Matrix
  {
    Matrix() {}
    Matrix(std::vector<int> dimensions, T initial_value = 0)
        : dims(dimensions)
    {
      for (auto d : dimensions)
        n_entries *= d;
      entries = std::vector<T>(n_entries, initial_value);
      if (n_entries != entries.size())
        throw std::invalid_argument("Invalid matrix size.");
    }
    T &operator[](std::vector<int> &indices)
    {
      if (indices.size() != dims.size())
        throw std::invalid_argument("Invalid number of indices.");
      int index = 0;
      for (size_t i = indices.size() - 1; i > 0; --i)
      {
        if (indices[i] >= dims[i])
          throw std::invalid_argument("Index out of range.");
        int a = indices[i];
        for (size_t d = 0; d < dims.size(); ++d)
        {
          if (d < i)
          {
            a *= dims[d];
          }
        }
        index += a;
      }
      index += indices[0];
      if (index > n_entries || index < 0)
        throw std::invalid_argument("Index out of range.");
      return entries[index];
    }
    std::vector<int> dims;
    size_t n_entries = 1;

  private:
    std::vector<T> entries;
  };

  typedef std::pair<double, double> Interval;

  const double eps = 0.01;

  const double INF = std::numeric_limits<double>::infinity();

  int random_index(const int n);

  template <typename Iter>
  void shuffle_vector(Iter first, Iter last)
  {
    int n = std::distance(first, last);
    while (n > 1)
    {
      int k = random_index(n--);
      std::swap(*(first + n), *(first + k));
    }
  };

  std::vector<int> to_std_vec(Rcpp::IntegerVector rv);
  std::vector<double> to_std_vec(Rcpp::NumericVector rv);
  std::vector<std::vector<double>> to_std_vec(Rcpp::NumericMatrix rv);
  std::set<int> to_std_set(Rcpp::NumericVector rv);
  std::set<int> to_std_set(Rcpp::IntegerVector rv);

  // functions for converting R and Cpp types

  Rcpp::IntegerVector from_std_set(std::set<int> v);
  Rcpp::IntegerVector from_std_vec(std::vector<int> v);
  Rcpp::NumericVector from_std_vec(std::vector<double> v);
  Rcpp::NumericMatrix from_std_vec(std::vector<std::vector<double>> v);

  //  ----------------- overload of vector operators -----------------

  // enable to subtract scalar from vector entries
  template <typename T, typename S>
  std::vector<T> operator-(const std::vector<T> &vec_a, S val)
  {
    std::vector<T> res = vec_a;
    for (auto &entry : res)
      entry -= val;
    return res;
  }

  // enable to subtract scalar from vector entries inplace
  template <typename T, typename S>
  void operator-=(std::vector<T> &vec, S val)
  {
    for (auto &entry : vec)
      entry -= val;
  }

  // enable to add scalar to vector entries
  template <typename T, typename S>
  std::vector<T> operator+(const std::vector<T> &vec_a, S val)
  {
    std::vector<T> res = vec_a;
    for (auto &entry : res)
      entry += val;
    return res;
  }

  // enable to add scalar to vector entries inplace
  template <typename T, typename S>
  void operator+=(std::vector<T> &vec, S val)
  {
    for (auto &entry : vec)
      entry += val;
  }

  // enable to multiply vector entries by scalar
  template <typename T, typename S>
  std::vector<T> operator*(const std::vector<T> &vec_a, S val)
  {
    std::vector<T> res = vec_a;
    for (auto &entry : res)
      entry *= val;
    return res;
  }

  // enable to multiply vector entries by scalar
  template <typename T, typename S>
  std::vector<T> operator*(S val, const std::vector<T> &vec_a)
  {
    std::vector<T> res = vec_a;
    for (auto &entry : res)
      entry *= val;
    return res;
  }

  // enable to multiply vector entries by scalar inplace
  template <typename T, typename S>
  void operator*=(std::vector<T> &vec, S val)
  {
    for (auto &entry : vec)
      entry *= val;
  }

  // enable to divide vector entries by scalar
  template <typename T, typename S>
  std::vector<T> operator/(const std::vector<T> &vec_a, S val)
  {
    std::vector<T> res = vec_a;
    for (auto &entry : res)
      entry /= val;
    return res;
  }

  // enable to divide vector entries by scalar inplace
  template <typename T, typename S>
  void operator/=(std::vector<T> &vec, S val)
  {
    for (auto &entry : vec)
      entry /= val;
  }

  // element-wise exp() of vector
  template <typename T>
  std::vector<T> exp(const std::vector<T> &vec)
  {
    std::vector<T> res = vec;
    for (auto &entry : res)
      entry = exp(entry);
    return res;
  }

  /**
   * \brief Extract keys from a std::map as vector of arbitrary type.
   *
   * \param map with arbitrary key and value type.
   */
  template <typename KT, typename VT>
  std::vector<KT> getKeys(std::map<KT, VT, setComp> &m)
  {
    std::vector<KT> keys;
    for (const auto &entry : m)
    {
      keys.push_back(entry.first);
    }
    return keys;
  }

  template <typename VT>
  std::vector<std::vector<VT>> transpose(std::vector<std::vector<VT>> &mat)
  {

    if (mat.size() <= 0)
      throw std::invalid_argument("Matrix is empty.");
    size_t value_size = mat[0].size();

    std::vector<std::vector<VT>> columns(value_size);
    for (auto vec : mat)
    {
      for (size_t n_col = 0; n_col < value_size; ++n_col)
      {
        columns[n_col].push_back(vec[n_col]);
      }
    }
    return columns;
  }

  /**
   * \brief Calculates median of the vector.
   *
   * \param vec a vector of arbitrary type.
   */
  template <typename VT>
  VT calcMedian(std::vector<VT> &vec)
  {

    // sort vector
    std::sort(vec.begin(), vec.end());
    size_t s = vec.size();

    // differ between even and odd case
    if (s % 2 == 0)
      return (vec[s / 2 - 1] + vec[s / 2]) / 2;
    return vec[s / 2];
  }

  /**
   * \brief Calculates row- or columnwise median of a matrix.
   *
   * \param mat a matrix of arbitrary type.
   */
  template <typename VT>
  std::vector<VT> calcMedian(std::vector<std::vector<VT>> &mat, bool colwise = true)
  {

    std::vector<VT> res;

    if (colwise)
    {
      std::vector<std::vector<VT>> columns = transpose(mat);

      for (auto vec : columns)
      {
        res.push_back(calcMedian(vec));
      }
    }
    else
    {
      for (auto vec : mat)
      {
        res.push_back(calcMedian(vec));
      }
    }

    return res;
  }

  /**
   * \brief Calculates row- or columnwise median of a matrix.
   *
   * \param mat a matrix of arbitrary type.
   * \param indices a vector giving the rows/columns to calculate for.
   * \param colwise a bool indicating if calculate row or columnwise.
   */
  template <typename VT>
  std::vector<VT> calcMedian(const std::vector<std::vector<VT>> &mat, const std::vector<int> &indices, bool colwise = true)
  {

    if (mat.size() == 0)
      throw std::invalid_argument("calcMedian: Matrix empty - no data provided.");

    int colSize = mat[0].size(), rowSize = mat.size(), index;

    std::vector<VT> res;
    std::vector<VT> tmp = std::vector<VT>(indices.size(), 0);

    if (colwise)
    {
      res = std::vector<VT>(colSize, 0);

      for (int col = 0; col < colSize; ++col)
      {

        index = 0;
        for (int row : indices)
        {
          if (row >= rowSize)
            throw std::invalid_argument("Row index out of range.");
          tmp[index] = mat[row][col];
          ++index;
        }
        res[col] = calcMedian(tmp);
      }
    }
    else
    {
      // todo
    }

    return res;
  }

  /**
   * \brief Calculate mean of a vector.
   *
   * \param vec a vector of arbitrary type.
   */
  template <typename VT>
  VT calcMean(std::vector<VT> &vec)
  {
    if (vec.empty())
      return 0;
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  }

  /**
   * \brief Calculate  row- or columnwise mean of a matrix.
   *
   * \param mat a matrix of arbitrary type.
   */
  template <typename VT>
  std::vector<VT> calcMean(std::vector<std::vector<VT>> &mat, bool colwise = true)
  {

    if (mat.size() == 0)
      throw std::invalid_argument("calcMean: Matrix empty - no data provided.");

    int colSize = mat[0].size(), rowSize = mat.size();
    std::vector<VT> res(std::max(colSize, rowSize), 0);

    if (colwise)
    {
      res = std::vector<VT>(colSize, 0);
      for (int col = 0; col < colSize; ++col)
      {
        for (int row = 0; row < rowSize; ++row)
        {
          res[col] += mat[row][col];
        }
        res[col] /= rowSize;
      }
    }
    else
    {
      res = std::vector<VT>(rowSize, 0);
      for (int row = 0; row < rowSize; ++row)
      {
        for (int col = 0; col < colSize; ++col)
        {
          res[row] += mat[row][col];
        }
        res[row] /= colSize;
      }
    }

    return res;
  }

  template <typename VT>
  std::vector<VT> calcMean(const std::vector<std::vector<VT>> &mat, const std::vector<int> &indices, bool colwise = true)
  {

    if (mat.size() == 0)
      throw std::invalid_argument("calcMean: Matrix empty - no data provided.");

    int colSize = mat[0].size(), rowSize = mat.size();
    std::vector<VT> res;

    if (colwise)
    {
      res = std::vector<VT>(colSize, 0);
      rowSize = indices.size();
      for (int col = 0; col < colSize; ++col)
      {
        for (int row : indices)
        {
          if (row >= (int)mat.size() || row < 0)
            throw std::invalid_argument("Row index out of range.");
          res[col] += mat[row][col];
        }
        res[col] /= rowSize;
      }
    }
    else
    {
      res = std::vector<VT>(rowSize, 0);
      colSize = indices.size();
      for (int row = 0; row < rowSize; ++row)
      {
        for (int col : indices)
        {
          if (col >= (int)mat[0].size() || col < 0)
            throw std::invalid_argument("Column index out of range.");
          res[row] += mat[row][col];
        }
        res[row] /= colSize;
      }
    }

    return res;
  }

  // enable to add two vectors
  template <typename T>
  std::vector<T> operator+(const std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    std::vector<T> res = vec_a;
    for (unsigned int i = 0; i < vec_b.size(); ++i)
      res[i] += vec_b[i];
    return res;
  }

  // enable to add two vectors inplace
  template <typename T>
  void operator+=(std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    for (unsigned int i = 0; i < vec_b.size(); ++i)
      vec_a[i] += vec_b[i];
  }

  // enable to subtract two vectors
  template <typename T>
  std::vector<T> operator-(const std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    std::vector<T> res = vec_a;
    for (unsigned int i = 0; i < vec_b.size(); ++i)
      res[i] -= vec_b[i];
    return res;
  }

  // enable to subtract two vectors inplace
  template <typename T>
  void operator-=(std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    for (unsigned int i = 0; i < vec_b.size(); ++i)
      vec_a[i] -= vec_b[i];
  }

  // enable to multiply two vectors
  template <typename T>
  std::vector<T> operator*(const std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    std::vector<T> res = vec_a;
    for (unsigned int i = 0; i < vec_b.size(); ++i)
      res[i] *= vec_b[i];
    return res;
  }

  // enable to multiply two vectors inplace
  template <typename T>
  void operator*=(std::vector<T> &vec_a, const std::vector<T> &vec_b)
  {
    if (vec_a.size() != vec_b.size())
      throw std::invalid_argument("The two vectors are not of same size.");

    for (unsigned int i = 0; i < vec_b.size(); ++i)
      vec_a[i] *= vec_b[i];
  }

} // namespace utils

#endif // UTILS_H
