#include <cmath>

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/hydrology.hpp"
#include "highmap/op.hpp"

#include "macrologger.h"

namespace hmap
{

void hydraulic_ridge(Array &z,
                     float  talus,
                     float  intensity,
                     float  smoothing_factor,
                     float  noise_ratio,
                     uint   seed)
{
  // erosion depth
  float ze_max = 2.5f;
  Array ze = flow_accumulation_dinf(z, talus);
  ze = log10(ze);
  minimum_smooth(ze, ze_max, ze_max);
  remap(ze);

  float landing_width_ratio = 0.1f;

  thermal_scree(ze,
                talus,
                seed,
                noise_ratio,
                0.f, // ze_min
                2.f * ze_max,
                smoothing_factor,
                landing_width_ratio,
                false);

  z -= intensity * ze;
}

} // namespace hmap
