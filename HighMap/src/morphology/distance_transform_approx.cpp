/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/math.hpp"

namespace hmap
{

Array distance_transform_approx(const Array &array,
                                bool         return_squared_distance)
{
  Vec2<int> shape = array.shape;

  Array edt(shape);

  // initialize the output with INF for background and 0 for foreground
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      edt(i, j) = array(i, j) > 0.f ? 0.f : std::numeric_limits<float>::max();

  // first pass: top-left to nottom-right
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (edt(i, j) > 0)
      {
        if (i > 0) edt(i, j) = std::min(edt(i, j), edt(i - 1, j) + 1.f);
        if (j > 0) edt(i, j) = std::min(edt(i, j), edt(i, j - 1) + 1.f);
        if (i > 0 && j > 0)
          edt(i, j) = std::min(edt(i, j), edt(i - 1, j - 1) + M_SQRT2f);
      }

  // second pass: bottom-right to top-left
  for (int j = shape.y - 1; j >= 0; --j)
    for (int i = shape.x - 1; i >= 0; --i)
      if (edt(i, j) > 0)
      {
        if (i < shape.x - 1)
          edt(i, j) = std::min(edt(i, j), edt(i + 1, j) + 1.f);
        if (j < shape.y - 1)
          edt(i, j) = std::min(edt(i, j), edt(i, j + 1) + 1.f);
        if (i < shape.x - 1 && j < shape.y - 1)
          edt(i, j) = std::min(edt(i, j), edt(i + 1, j + 1) + M_SQRT2f);
      }

  if (return_squared_distance)
    return edt * edt;
  else
    return edt;
}

Array distance_transform_manhattan(const Array &array,
                                   bool         return_squared_distance)
{
  Vec2<int> shape = array.shape;

  Array edt(shape);

  // initialize the output with INF for background and 0 for foreground
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      edt(i, j) = array(i, j) > 0.f ? 0.f : std::numeric_limits<float>::max();

  // first pass: top-left to nottom-right
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (edt(i, j) > 0)
      {
        if (i > 0) edt(i, j) = std::min(edt(i, j), edt(i - 1, j) + 1.f);
        if (j > 0) edt(i, j) = std::min(edt(i, j), edt(i, j - 1) + 1.f);
      }

  // second pass: bottom-right to top-left
  for (int j = shape.y - 1; j >= 0; --j)
    for (int i = shape.x - 1; i >= 0; --i)
      if (edt(i, j) > 0)
      {
        if (i < shape.x - 1)
          edt(i, j) = std::min(edt(i, j), edt(i + 1, j) + 1.f);
        if (j < shape.y - 1)
          edt(i, j) = std::min(edt(i, j), edt(i, j + 1) + 1.f);
      }

  if (return_squared_distance)
    return edt * edt;
  else
    return edt;
}

} // namespace hmap
