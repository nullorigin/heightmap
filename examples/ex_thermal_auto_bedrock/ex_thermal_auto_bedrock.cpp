#include <iostream>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/erosion.hpp"
#include "highmap/noise.hpp"
#include "highmap/op.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {2.f, 2.f};
  int                      seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  auto        z0 = z;

  hmap::thermal_auto_bedrock(z, 0.3f / shape[0]);

  z0.to_png("ex_thermal_auto_bedrock0.png", hmap::cmap::gray);
  z.to_png("ex_thermal_auto_bedrock1.png", hmap::cmap::gray);
}