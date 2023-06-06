#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  int                      seed = 1;

  auto z1 = hmap::perlin(shape, res, seed);
  auto z2 = hmap::perlin_billow(shape, res, seed);
  auto z3 = hmap::perlin_mix(shape, res, seed);

  hmap::export_banner_png("ex_perlin.png", {z1, z2, z3}, hmap::cmap::viridis);
}
