/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/geometry.hpp"
#include "highmap/kernels.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

namespace hmap
{

void alter_elevation(Array      &array,
                     Cloud      &cloud,
                     int         ir,
                     float       footprint_ratio,
                     Vec2<float> shift,
                     Vec2<float> scale)
{
  Array amp = Array(array.shape);

  for (auto &p : cloud.points)
  {
    int ic = (int)((p.x - shift.x) / scale.x * array.shape.x);
    int jc = (int)((p.y - shift.y) / scale.y * array.shape.y);

    // kernel size
    int   nk = (int)((2 * ir + 1) * (1.f + footprint_ratio * std::abs(p.v)));
    Array kernel = cubic_pulse(Vec2(nk, nk));

    // --
    // truncate kernel to make it fit into the heightmap array
    int nk_i0 = (int)(std::floor(0.5f * kernel.shape.x)); // left
    int nk_i1 = kernel.shape.x - nk_i0;                   // right
    int nk_j0 = (int)(std::floor(0.5f * kernel.shape.y));
    int nk_j1 = kernel.shape.y - nk_j0;

    int ik0 = std::max(0, nk_i0 - ic);
    int jk0 = std::max(0, nk_j0 - jc);
    int ik1 = std::min(kernel.shape.x,
                       kernel.shape.x - (ic + nk_i1 - amp.shape.x));
    int jk1 = std::min(kernel.shape.y,
                       kernel.shape.y - (jc + nk_j1 - amp.shape.y));

    // where it goes in the array
    int i0 = std::max(ic - nk_i0, 0);
    int j0 = std::max(jc - nk_j0, 0);
    // i1 = std::min(i + nk1, amp.shape.x);
    // j1 = std::min(j + nk1, amp.shape.y);

    float array_min = std::numeric_limits<float>::max();
    float array_max = std::numeric_limits<float>::min();

    for (int i = ik0; i < ik1; i++)
      for (int j = jk0; j < jk1; j++)
      {
        array_min = std::min(array(i, j), array_min);
        array_max = std::max(array(i, j), array_max);
      }
    float array_ptp = array_max - array_min;

    // define scaling by roughly determining if the modification is
    // applied to a peak or a valley
    float scaling;
    if (array(ic, jc) > 0.5f * array_ptp + array_min)
      scaling = (array(ic, jc) - array_min) / array_ptp;
    else
      scaling = (array_max - array(ic, jc)) / array_ptp;

    if (array_ptp > 0.f)
    {
      for (int i = ik0; i < ik1; i++)
        for (int j = jk0; j < jk1; j++)
          amp(i - ik0 + i0, j - jk0 + j0) += p.v * kernel(i, j) * scaling;
    }
  }

  array += amp;
}

} // namespace hmap
