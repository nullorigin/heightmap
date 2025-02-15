/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include "macrologger.h"

#define LAPLACE_SIGMA 0.2f
#define LAPLACE_ITERATIONS 1

namespace hmap
{

//----------------------------------------------------------------------
// Main operator
//----------------------------------------------------------------------

void thermal_flatten(Array       &z,
                     const Array &talus,
                     const Array &bedrock,
                     int          iterations)
{
  std::vector<int>   di = DI;
  std::vector<int>   dj = DJ;
  std::vector<float> c = CD;
  const uint         nb = di.size();

  // main loop
  for (int it = 0; it < iterations; it++)
  {
    // modify neighbor search at each iterations to limit numerical
    // artifacts
    std::rotate(di.begin(), di.begin() + 1, di.end());
    std::rotate(dj.begin(), dj.begin() + 1, dj.end());
    std::rotate(c.begin(), c.begin() + 1, c.end());

    for (int j = 1; j < z.shape.y - 1; j++)
      for (int i = 1; i < z.shape.x - 1; i++)
      {
        if (z(i, j) > bedrock(i, j))
        {
          float dmax = 0.f;
          int   ka = -1;

          for (uint k = 0; k < nb; k++)
          {
            float dz = z(i, j) - z(i + di[k], j + dj[k]);
            if (dz > dmax)
            {
              dmax = dz;
              ka = k;
            }
          }

          if (dmax > 0.f and dmax < talus(i, j))
          {
            float amount = 0.5f * dmax * c[ka];
            z(i, j) -= amount;
            z(i + di[ka], j + dj[ka]) += amount;
          }
        }
      }
  }

  // clean-up: fix boundaries, remove spurious oscillations and make
  // sure final elevation is not lower than the bedrock
  extrapolate_borders(z);
  laplace(z, LAPLACE_SIGMA, LAPLACE_ITERATIONS);
  clamp_min(z, bedrock);
}

//----------------------------------------------------------------------
// Overloading
//----------------------------------------------------------------------

void thermal_flatten(Array &z, float talus, int iterations)
{
  Array talus_map = constant(z.shape, talus);
  Array bedrock = constant(z.shape, z.min() - z.ptp());
  thermal_flatten(z, talus_map, bedrock, iterations);
}

} // namespace hmap
