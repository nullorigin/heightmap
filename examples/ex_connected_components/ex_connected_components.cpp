#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  uint                     seed = 5;

  hmap::Array z = hmap::perlin(shape, res, seed);
  hmap::clamp_min(z, 0.f);

  hmap::Array labels = hmap::connected_components(z);

  z.to_png("ex_connected_components0.png", hmap::cmap::inferno);
  labels.to_png("ex_connected_components1.png", hmap::cmap::nipy_spectral);
}