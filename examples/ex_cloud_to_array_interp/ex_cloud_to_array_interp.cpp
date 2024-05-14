#include "highmap.hpp"

int main(void)
{
  hmap::Vec2<int> shape = {256, 256};
  int             seed = 1;

  hmap::Vec4<float> bbox = {-1.f, 0.f, 0.5f, 1.5f};

  hmap::Cloud cloud = hmap::Cloud(10, seed, bbox);

  hmap::Array z0 = hmap::Array(shape);
  cloud.to_array(z0, bbox);

  hmap::Array z1 = hmap::Array(shape);
  cloud.to_array_interp(z1, bbox);

  hmap::Array nx = 0.5f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                          shape,
                                          {2.f, 2.f},
                                          seed++);
  hmap::Array ny = 0.5f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                          shape,
                                          {2.f, 2.f},
                                          seed++);

  int         interpolation_method = 0; // 0 linear, 1 thinplate
  hmap::Array z2 = hmap::Array(shape);
  cloud.to_array_interp(z2, bbox, interpolation_method, &nx, &ny);

  hmap::Array z3 = hmap::Array(shape);
  interpolation_method = 1;
  cloud.to_array_interp(z3, bbox, interpolation_method, &nx, &ny);

  hmap::export_banner_png("ex_cloud_to_array_interp.png",
                          {z0, z1, z2, z3},
                          hmap::cmap::inferno);
}
