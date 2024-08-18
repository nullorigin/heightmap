#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);
  auto z0 = z;

  float vmax = 0.8f;
  float radius = 0.1f;
  float k_smoothing = 0.5f;

  hmap::hydraulic_blur(z, vmax, radius, k_smoothing);

  hmap::export_banner_png("ex_hydraulic_blur.png",
                          {z0, z},
                          hmap::Cmap::TERRAIN,
                          true);
}
