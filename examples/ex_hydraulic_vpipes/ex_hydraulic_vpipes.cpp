#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {512, 512};
  hmap::Vec2<float> res = {2.f, 2.f};
  int               seed = 2;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  hmap::remap(z);
  auto z0 = z;

  hmap::hydraulic_vpipes(z, 300);

  z.infos();

  hmap::export_banner_png("ex_hydraulic_vpipes.png",
                          {z0, z},
                          hmap::cmap::terrain);
}
