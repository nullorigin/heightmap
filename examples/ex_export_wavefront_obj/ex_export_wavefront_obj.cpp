#include "highmap.hpp"

int main(void)
{

  // heightmaps

  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {2.f, 2.f};
  int               seed = 1;

  hmap::Array z = hmap::fbm_simplex(shape, res, seed);
  hmap::clamp_min_smooth(z, 0.f, 0.2f);
  hmap::remap(z);

  hmap::export_wavefront_obj("hmap.obj", z, hmap::mesh_type::tri_optimized);
  hmap::export_wavefront_obj("hmap_quad.obj", z, hmap::mesh_type::quad);
  hmap::export_wavefront_obj("hmap_tri.obj", z, hmap::mesh_type::tri);

  z.to_png("out.png", hmap::cmap::gray);

  // lines

  hmap::Vec4<float> bbox = {0.f, 1.f, 0.f, 1.f};

  int        npoints = 50;
  bool       closed = false;
  hmap::Path path = hmap::Path(npoints, seed, bbox, closed);
  path.set_values(0.f);
  path.set_values_from_array(z, bbox);
  path.reorder_nns();

  hmap::export_wavefront_obj("path.obj", path);
}