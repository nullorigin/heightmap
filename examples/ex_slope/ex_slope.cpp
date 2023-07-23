#include "highmap.hpp"

int main(void)
{
  const std::vector<int>   shape = {256, 256};
  const std::vector<float> res = {4.f, 4.f};
  int                      seed = 1;

  hmap::Array noise = 0.2f * hmap::fbm_perlin(shape, res, seed);

  float       talus = 1.f / shape[0];
  hmap::Array sx = hmap::slope_x(shape, talus);
  hmap::Array sy = hmap::slope_y(shape, talus);

  hmap::Array oblique = hmap::slope_x(shape, talus, &noise) +
                        hmap::slope_y(shape, talus, &noise);

  hmap::Array valley = maximum_smooth(hmap::slope_x(shape, talus, &noise),
                                      hmap::slope_y(shape, talus, &noise),
                                      0.1f);

  hmap::export_banner_png("ex_slope.png",
                          {sx, sy, oblique, valley},
                          hmap::cmap::terrain);
}