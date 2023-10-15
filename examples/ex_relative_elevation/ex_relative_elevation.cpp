#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  hmap::remap(z);

  int  ir = 32;
  auto zr = hmap::relative_elevation(z, ir);

  hmap::export_banner_png("ex_relative_elevation.png",
                          {z, zr},
                          hmap::cmap::terrain);
}
