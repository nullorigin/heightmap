#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {2.f, 2.f};
  int                      seed = 1;

  hmap::Array z = hmap::ridged_perlin(shape, res, seed);

  z.to_png("ex_ridged_perlin.png", hmap::cmap::terrain, true);
}