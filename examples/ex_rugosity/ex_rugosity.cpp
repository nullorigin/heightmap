#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;
  int               ir = 5;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  hmap::Array c = hmap::rugosity(z, ir);

  z.to_png("ex_rugosity0.png", hmap::cmap::viridis);
  c.to_png("ex_rugosity1.png", hmap::cmap::viridis);
}
