#include <iostream>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/noise.hpp"
#include "highmap/op.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  int                      seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  hmap::remap(z);
  hmap::gamma_correction(z, 0.5f);
  z.to_png("ex_gamma_correction.png", hmap::cmap::gray);
}