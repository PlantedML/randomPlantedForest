#ifndef RPF_DIFFBUF_HPP
#define RPF_DIFFBUF_HPP

#include <vector>
#include <cstddef>
#include <algorithm>

// N-D difference buffer for axis-aligned rectangular range updates
// and reconstruction via inclusive prefix scans along each dimension.

namespace rpf_diff
{
  template <typename T>
  struct NDArray
  {
    std::vector<int> dims;      // logical dimensions
    std::vector<T> data;        // flat row-major data

    NDArray() {}
    explicit NDArray(const std::vector<int> &d, const T &init = T()) : dims(d)
    {
      size_t n = 1; for (int v : d) n *= (size_t)v; data.assign(n, init);
    }

    inline size_t offset(const std::vector<int> &idx) const
    {
      size_t off = 0; size_t stride = 1;
      for (size_t k = 0; k < dims.size(); ++k)
      {
        off += (size_t)idx[k] * stride; stride *= (size_t)dims[k];
      }
      return off;
    }

    inline T &at(const std::vector<int> &idx) { return data[offset(idx)]; }
    inline const T &at(const std::vector<int> &idx) const { return data[offset(idx)]; }
  };

  // Apply a constant add v onto a closed-open hyper-rectangle [lo, hi) via difference corners
  template <typename T>
  void add_rect(NDArray<T> &diff, const std::vector<int> &lo, const std::vector<int> &hi, const T &v)
  {
    const size_t d = diff.dims.size();
    // iterate over 2^d corners
    std::vector<int> corner(d, 0);
    for (;;) {
      int flips = 0; for (size_t k = 0; k < d; ++k) if (corner[k]) ++flips;
      T sign = (flips % 2 == 0) ? v : (v * (-1));
      std::vector<int> idx(d);
      for (size_t k = 0; k < d; ++k) idx[k] = corner[k] ? hi[k] : lo[k];
      diff.at(idx) += sign;
      size_t pos = 0;
      while (pos < d) { if (corner[pos] == 0) { corner[pos] = 1; break; } corner[pos] = 0; ++pos; }
      if (pos == d) break;
    }
  }

  // Inclusive prefix scan along each dimension in-place (converts diff -> values)
  template <typename T>
  void inclusive_scan_inplace(NDArray<T> &arr)
  {
    const size_t d = arr.dims.size();
    if (d == 0) return;
    for (size_t axis = 0; axis < d; ++axis)
    {
      // number of slabs orthogonal to axis
      size_t nslab = 1; for (size_t k = 0; k < d; ++k) if (k != axis) nslab *= (size_t)arr.dims[k];
      std::vector<int> slab_idx(d, 0);
      for (size_t s = 0; s < nslab; ++s)
      {
        // decode slab index into coordinates for all dims except axis
        size_t tmp = s;
        for (size_t k = 0; k < d; ++k)
        {
          if (k == axis) continue;
          slab_idx[k] = (int)(tmp % (size_t)arr.dims[k]);
          tmp /= (size_t)arr.dims[k];
        }
        // inclusive scan along axis
        std::vector<int> run = slab_idx; run[axis] = 0;
        T acc = arr.at(run); acc -= acc; // zero of correct shape
        for (int t = 0; t < arr.dims[axis]; ++t)
        {
          run[axis] = t; acc += arr.at(run); arr.at(run) = acc;
        }
      }
    }
  }
}

#endif // RPF_DIFFBUF_HPP


