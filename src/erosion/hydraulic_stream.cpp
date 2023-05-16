#include <cmath>

#include "highmap/array.hpp"
#include "highmap/hydrology.hpp"
#include "highmap/op.hpp"

#include <iostream>

namespace hmap
{

void hydraulic_stream(Array &z,
                      Array &z_bedrock,
                      float  c_erosion,
                      float  talus_ref,
                      float  clipping_ratio)
{
  // use flow accumulation to determine erosion intensity
  Array facc = flow_accumulation_dinf(z, talus_ref);

  // clip large flow accumulation values using a value loosely based
  // on the standard deviation (of an equivalent symetric
  // distribution)
  float vmax = clipping_ratio * std::pow(facc.sum() / (float)facc.size(), 0.5f);
  clamp(facc, 0.f, vmax);
  remap(facc);

  Array ze = z - c_erosion * facc;
  z = maximum(z_bedrock, ze);
}

} // namespace hmap