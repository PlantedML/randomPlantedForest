#ifndef GRID_H
#define GRID_H

#include "helper.hpp"

using namespace utils;

namespace grid
{
  typedef std::vector<int> Dimension;
  typedef std::vector<Dimension> Space;
  typedef std::vector<typename Dimension::iterator> Point;

  class NDGrid
  {

  private:
    Space space;
    bool first = true;
    Point current;

  public:
    Dimension dimensions;

  public:
    NDGrid();
    NDGrid(std::vector<int> dims);
    std::vector<int> getPoint();
    bool nextPoint();
  };

  struct LeafGrid
  {
    NDGrid grid;
    std::vector<std::vector<double>> lim_list;
    utils::Matrix<int> individuals;
    utils::Matrix<std::vector<double>> values;
    // rpf::Matrix<std::vector<double>> values_old;
  };
};

#endif