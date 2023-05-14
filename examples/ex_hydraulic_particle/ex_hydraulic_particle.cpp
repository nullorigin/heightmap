#include <iostream>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/erosion.hpp"
#include "highmap/noise.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  int                      nparticles = 10000;
  int                      seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  auto        z0 = z;

  hmap::hydraulic_particle(z, nparticles, seed);
  z0.to_png("ex_hydraulic_particle0.png", hmap::cmap::gray);
  z.to_png("ex_hydraulic_particle1.png", hmap::cmap::gray);
}