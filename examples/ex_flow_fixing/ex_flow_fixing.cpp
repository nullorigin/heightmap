#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int> shape = {256, 256};
  shape = {512, 512};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 2;

  hmap::Array z = hmap::fbm_simplex(shape, res, seed);
  hmap::remap(z);
  hmap::Array zf = hmap::flow_fixing(z);

  hmap::export_banner_png("ex_flow_fixing.png",
                          {z, zf},
                          hmap::cmap::terrain,
                          true);
}
