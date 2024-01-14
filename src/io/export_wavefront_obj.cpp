/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>

#include "WavefrontObjWriter/obj.h"
#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"
#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/io.hpp"
#include "highmap/op.hpp"

namespace hmap
{

void export_wavefront_obj(std::string  fname,
                          const Array &array,
                          int          mesh_type,
                          float        elevation_scaling,
                          float        max_error)
{
  // remove extension if already provided
  size_t      lastindex = fname.find_last_of(".");
  std::string fname_raw = fname.substr(0, lastindex);

  wow::Obj obj;

  std::vector<float> x = linspace(0.f, 1.f, array.shape.x);
  std::vector<float> y = linspace(0.f, 1.f, array.shape.y);

  switch (mesh_type)
  {

  case hmap::mesh_type::quad:
  {
    LOG_DEBUG("exporting OBJ (quad)");

    for (int i = 0; i < array.shape.x - 1; i++)
      for (int j = 0; j < array.shape.y - 1; j++)
      {
        obj.appendVertex(x[i], elevation_scaling * array(i, j), y[j]);
        obj.appendVertex(x[i + 1], elevation_scaling * array(i + 1, j), y[j]);
        obj.appendVertex(x[i + 1],
                         elevation_scaling * array(i + 1, j + 1),
                         y[j + 1]);
        obj.appendVertex(x[i], elevation_scaling * array(i, j + 1), y[j + 1]);
        obj.closeFace();
      }
  }
  break;

  case hmap::mesh_type::tri:
  {
    LOG_DEBUG("exporting OBJ (tri)");

    for (int i = 0; i < array.shape.x; i++)
      for (int j = 0; j < array.shape.y; j++)
        obj.appendVertex(x[i], elevation_scaling * array(i, j), y[j]);

    for (int i = 0; i < array.shape.x - 1; i++)
      for (int j = 0; j < array.shape.y - 1; j++)
      {
        obj.appendFace(wow::Face(array.linear_index(i, j),
                                 array.linear_index(i + 1, j),
                                 array.linear_index(i, j + 1)));

        obj.appendFace(wow::Face(array.linear_index(i + 1, j),
                                 array.linear_index(i + 1, j + 1),
                                 array.linear_index(i, j + 1)));
      }
  }
  break;

  case hmap::mesh_type::tri_optimized:
  {
    LOG_DEBUG("exporting OBJ (tri optimized)");

    const auto   p_hmap = std::make_shared<Heightmap>(array.shape.x,
                                                    array.shape.y,
                                                    array.get_vector());
    Triangulator tri(p_hmap);
    const int    max_triangles = 0;
    const int    max_points = 0;

    LOG_DEBUG("remeshing (Delaunay)");
    tri.Run(max_error, max_triangles, max_points);

    auto points = tri.Points(elevation_scaling);
    auto triangles = tri.Triangles();

    LOG_DEBUG("points: %ld (max: %d)", points.size(), array.size());
    LOG_DEBUG("triangles: %ld", triangles.size());

    float ax = 1.f / (float)array.shape.x;
    float ay = 1.f / (float)array.shape.y;

    for (size_t k = 0; k < points.size(); k++)
      obj.appendVertex(ax * points[k].x, points[k].z, ay * points[k].y);

    for (size_t k = 0; k < triangles.size(); k++)
      obj.appendFace(wow::Face(triangles[k].x, triangles[k].y, triangles[k].z));
  }
  break;
  }

  obj.output(fname_raw);
}

} // namespace hmap
