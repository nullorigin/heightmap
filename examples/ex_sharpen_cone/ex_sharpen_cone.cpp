#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z0 = hmap::fbm_perlin(shape, res, seed);
  auto        z1 = z0;

  int ir = 16;
  hmap::sharpen_cone(z1, ir);

  hmap::export_wavefront_obj("out0.obj", z0);
  hmap::export_wavefront_obj("out1.obj", z1);

  hmap::export_banner_png("ex_sharpen_clone.png", {z0, z1}, hmap::cmap::jet);
}