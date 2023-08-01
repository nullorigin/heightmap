/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "FastNoiseLite.h"

#include "highmap/array.hpp"
#include "highmap/primitives.hpp"

namespace hmap
{

hmap::Array worley(Vec2<int>   shape,
                   Vec2<float> kw,
                   uint        seed,
                   Array      *p_noise_x,
                   Array      *p_noise_y,
                   Vec2<float> shift,
                   Vec2<float> scale)
{
  hmap::Array   array = hmap::Array(shape);
  FastNoiseLite noise(seed);

  noise.SetFrequency(1.0f);
  noise.SetNoiseType(FastNoiseLite::NoiseType_Cellular);
  noise.SetCellularReturnType(FastNoiseLite::CellularReturnType_Distance);
  noise.SetCellularJitter(1.f);

  array = helper_get_noise(array,
                           kw,
                           p_noise_x,
                           p_noise_y,
                           shift,
                           scale,
                           [&noise](float x, float y)
                           { return noise.GetNoise(x, y); });
  return array;
}

} // namespace hmap
