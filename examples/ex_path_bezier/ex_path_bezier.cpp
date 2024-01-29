#include "highmap/array.hpp"
#include "highmap/geometry.hpp"
#include "highmap/io.hpp"

int main(void)
{
  hmap::Vec2<int> shape = {256, 256};
  int             seed = 3;

  hmap::Vec4<float> bbox = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Path        path = hmap::Path(10, seed, {1.2f, 1.8f, -0.3, 0.3f});
  path.reorder_nns();
  path.closed = true;

  auto z1 = hmap::Array(shape);
  path.to_array(z1, bbox);

  auto z2 = hmap::Array(shape);
  path.bezier(0.5f, 10); // curvature, point density
  path.to_array(z2, bbox);

  auto z3 = hmap::Array(shape);
  path.resample(0.05f);
  path.fractalize(2, seed);
  path.to_array(z3, bbox);

  hmap::export_banner_png("ex_path_bezier.png",
                          {z1, z2, z3},
                          hmap::cmap::inferno);
}
