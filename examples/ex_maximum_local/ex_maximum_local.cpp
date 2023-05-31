#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  int                      seed = 1;
  int                      radius = 5;

  hmap::Array z = hmap::fbm_perlin(shape, res, seed);
  hmap::Array zmin = hmap::minimum_local(z, radius);
  hmap::Array zmax = hmap::maximum_local(z, radius);
  hmap::Array zdisk = hmap::maximum_local_disk(z, radius);

  hmap::export_banner_png("ex_maximum_local.png",
                          {z, zmin, zmax, zdisk},
                          hmap::cmap::viridis);
}
