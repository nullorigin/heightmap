#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int> shape = {256, 256};
  shape = {2048, 2048};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);

  hmap::clamp_min(z, 0.f);

  hmap::Timer::Start("exact");
  auto d0 = hmap::distance_transform(z);
  hmap::Timer::Stop("exact");

  hmap::Timer::Start("approx.");
  auto d1 = hmap::distance_transform_approx(z);
  hmap::Timer::Stop("approx.");

  hmap::Timer::Start("manhattan");
  auto d2 = hmap::distance_transform_manhattan(z);
  hmap::Timer::Stop("manhattan");

  z.to_png("ex_distance_transform0.png", hmap::Cmap::VIRIDIS);
  d0.to_png("ex_distance_transform1.png", hmap::Cmap::VIRIDIS);
  d1.to_png("ex_distance_transform2.png", hmap::Cmap::VIRIDIS);
  d2.to_png("ex_distance_transform3.png", hmap::Cmap::VIRIDIS);
}
