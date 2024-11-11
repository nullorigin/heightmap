#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {2.f, 2.f};
  int               seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);
  auto z0 = z;

  int   iterations = 10;
  float intensity = 0.001f;
  float talus = 2.f / shape.x;

  hmap::thermal_schott(z, talus, iterations, intensity);

  // align talus constraint with the elevation
  hmap::Array talus_map = z;
  hmap::remap(talus_map, talus / 100.f, talus);

  hmap::thermal_schott(z, talus_map, iterations, intensity);

  hmap::export_banner_png("ex_thermal_schott.png",
                          {z0, z},
                          hmap::Cmap::TERRAIN,
                          true);  
}
