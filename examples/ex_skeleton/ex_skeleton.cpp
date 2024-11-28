#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {2.f, 8.f};
  int               seed = 2;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z, -1.f, 0.4f);

  hmap::make_binary(z);

  auto sk = hmap::skeleton(z);

  hmap::export_banner_png("ex_skeleton.png", {z, sk}, hmap::Cmap::GRAY);
}
