#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/io.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  auto        z0 = z;

  // talus limit defined locally
  auto talus = hmap::constant(shape, 0.2f / (float)shape.x);
  hmap::sediment_deposition(z, talus);

  hmap::export_banner_png("ex_sediment_deposition.png",
                          {z0, z},
                          hmap::cmap::terrain,
                          true);
}
