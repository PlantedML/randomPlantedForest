#include "grid.hpp"
#include <vector>
#include <iostream>

using namespace grid;

NDGrid::NDGrid() {}

NDGrid::NDGrid(std::vector<int> dims) : dimensions(dims)
{
  // fill space with dimensions
  for (const auto &dim : dims)
  {
    Dimension d;
    for (int i = 0; i < dim; ++i)
    {
      d.push_back(i);
    }
    space.push_back(d);
  }
}

std::vector<int> NDGrid::getPoint()
{
  std::vector<int> gridPoint(current.size());
  for (unsigned int i = 0; i < current.size(); ++i)
  {
    gridPoint[i] = *current[i];
  }
  return gridPoint;
}

bool NDGrid::nextPoint()
{

  // get first point in N-dimensional space
  if (first)
  {
    first = false;
    for (Space::iterator dims_it = space.begin(); dims_it != space.end(); ++dims_it)
    {
      current.push_back((*dims_it).begin());
    }
    return false;
  }

  // go to next point in space
  Space::iterator dims_it = space.begin();
  Point::iterator cur_it = current.begin();
  for (; dims_it != space.end(); ++dims_it, ++cur_it)
  {

    // check if next in dimension is at the end
    if (++(*cur_it) == (*dims_it).end())
    {

      // check if we have finished whole space
      if (dims_it == space.end() - 1)
      {
        // reset setup for next time of iteration
        first = true;
        current.clear();
        // stop running now
        return true;
      }
      // reset this dimension to begin
      // and go to next dimension
      *cur_it = (*dims_it).begin();
    }
    else
    {
      // next point is valid
      break;
    }
  }
  return false;
}
