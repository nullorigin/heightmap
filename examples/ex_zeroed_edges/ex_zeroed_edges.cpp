#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  float sigma = 1.f;

  auto z1 = z;
  auto z2 = z;
  auto z3 = z;
  hmap::zeroed_edges(z1, sigma, hmap::DistanceFunction::EUCLIDIAN);
  hmap::zeroed_edges(z2, sigma, hmap::DistanceFunction::CHEBYSHEV);
  hmap::zeroed_edges(z3, sigma, hmap::DistanceFunction::MANHATTAN);

  hmap::export_banner_png("ex_zeroed_edges.png",
                          {z, z1, z2, z3},
                          hmap::cmap::inferno);
}
