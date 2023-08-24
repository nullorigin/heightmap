#include "highmap.hpp"
#include "highmap/dbg.hpp"

int main(void)
{
  hmap::Timer timer = hmap::Timer();

  const hmap::Vec2<int>   shape = {512, 512};
  const hmap::Vec2<float> res = {4.f, 4.f};
  int                     seed = 2;

  timer.start("fbm_perlin");
  hmap::Array z = hmap::fbm_perlin(shape, res, seed, 8);
  timer.stop("fbm_perlin");

  z.to_png("out.png", hmap::cmap::gray);
}
