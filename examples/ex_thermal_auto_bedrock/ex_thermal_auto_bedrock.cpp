#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/io.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {2.f, 2.f};
  int               seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  auto        z0 = z;

  hmap::thermal_auto_bedrock(z, 0.1f / shape.x);

  hmap::export_banner_png("ex_thermal_auto_bedrock.png",
                          {z0, z},
                          hmap::cmap::terrain,
                          true);
}
