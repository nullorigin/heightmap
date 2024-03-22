#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::n_perlin, shape, res, seed);

  int         ir = 8;
  hmap::Array w = hmap::valley_width(z, ir);

  z.to_png("ex_valley_width0.png", hmap::cmap::terrain);
  w.to_png("ex_valley_width1.png", hmap::cmap::inferno);
}
