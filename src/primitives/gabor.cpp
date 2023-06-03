#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>

#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

#define SVD_RANK 2

namespace hmap
{

Array gabor(std::vector<int> shape,
            float            kw,
            float            angle,
            float            footprint_threshold)
{
  Array array = Array(shape);

  std::vector<float> x = linspace(-1.f, 1.f, array.shape[0]);
  std::vector<float> y = linspace(-1.f, 1.f, array.shape[1]);

  float width = std::sqrt(-0.5f * M_PI / std::log(footprint_threshold));
  float iw2 = 1.f / (width * width);
  float ca = std::cos(angle / 180.f * M_PI);
  float sa = std::sin(angle / 180.f * M_PI);

  for (int i = 0; i < array.shape[0]; i++)
    for (int j = 0; j < array.shape[1]; j++)
      array(i, j) = std::exp(-M_PI * (x[i] * x[i] + y[j] * y[j]) * 0.5f * iw2) *
                    std::cos(kw * (x[i] * ca + y[j] * sa));

  return array;
}

Array gabor_noise(std::vector<int> shape,
                  float            kw,
                  float            angle,
                  int              width,
                  float            density,
                  uint             seed)
{
  Array weight = white_sparse(shape, 0.f, 1.f, density, seed);
  Array kernel = gabor({width, width}, kw, angle);

  return convolve2d_svd(weight, kernel, SVD_RANK);
}

} // namespace hmap
