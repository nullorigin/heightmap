/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "FastNoiseLite.h"

#include "highmap/array.hpp"
#include "highmap/primitives.hpp"

namespace hmap
{

Array perlin(std::vector<int>   shape,
             std::vector<float> kw,
             uint               seed,
             Array             *p_noise_x,
             Array             *p_noise_y,
             std::vector<float> shift,
             std::vector<float> scale)
{
  Array         array = Array(shape);
  FastNoiseLite noise(seed);

  noise.SetFrequency(1.0f);
  noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);

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

Array perlin_billow(std::vector<int>   shape,
                    std::vector<float> kw,
                    uint               seed,
                    Array             *p_noise_x,
                    Array             *p_noise_y,
                    std::vector<float> shift,
                    std::vector<float> scale)
{
  Array         array = Array(shape);
  FastNoiseLite noise(seed);

  noise.SetFrequency(1.0f);
  noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);

  array = helper_get_noise(array,
                           kw,
                           p_noise_x,
                           p_noise_y,
                           shift,
                           scale,
                           [&noise](float x, float y) {
                             return 2.f * std::abs(noise.GetNoise(x, y)) - 1.f;
                           });
  return array;
}

Array perlin_mix(std::vector<int>   shape,
                 std::vector<float> kw,
                 uint               seed,
                 Array             *p_noise_x,
                 Array             *p_noise_y,
                 std::vector<float> shift,
                 std::vector<float> scale)
{
  Array         array = Array(shape);
  FastNoiseLite noise(seed);

  noise.SetFrequency(1.0f);
  noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);

  array = helper_get_noise(array,
                           kw,
                           p_noise_x,
                           p_noise_y,
                           shift,
                           scale,
                           [&noise](float x, float y)
                           {
                             return 0.5f * noise.GetNoise(x, y) +
                                    std::abs(noise.GetNoise(x, y)) - 0.5f;
                           });
  return array;
}

} // namespace hmap
