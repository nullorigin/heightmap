#include <cmath>

#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/primitives.hpp"

namespace hmap
{

Array caldera(std::vector<int>   shape,
              float              radius,
              float              sigma_inner,
              float              sigma_outer,
              float              z_bottom,
              Array             &noise,
              float              noise_r_amp,
              float              noise_z_ratio,
              std::vector<float> shift)
{
  Array z = Array(shape);
  int   ic = (int)((0.5f - shift[0]) * z.shape[0]);
  int   jc = (int)((0.5f - shift[1]) * z.shape[1]);

  float si2 = sigma_inner * sigma_inner;
  float so2 = sigma_outer * sigma_outer;

  for (int i = 0; i < z.shape[0]; i++)
    for (int j = 0; j < z.shape[1]; j++)
    {
      float r = std::hypot((float)(i - ic), (float)(j - jc)) - radius;

      r += noise_r_amp * (2 * noise(i, j) - 1);

      if (r < 0.f)
        z(i, j) = z_bottom + std::exp(-0.5f * r * r / si2) * (1 - z_bottom);
      else
        z(i, j) = 1 / (1 + r * r / so2);

      z(i, j) *= 1.f + noise_z_ratio * (2.f * noise(i, j) - 1.f);
    }

  return z;
}

Array caldera(std::vector<int>   shape,
              float              radius,
              float              sigma_inner,
              float              sigma_outer,
              float              z_bottom,
              std::vector<float> shift)
{
  Array noise = constant(shape, 0.f);
  Array z = caldera(shape,
                    radius,
                    sigma_inner,
                    sigma_outer,
                    z_bottom,
                    noise,
                    0.f,
                    0.f,
                    shift);
  return z;
}

Array crater(std::vector<int>   shape,
             float              radius,
             float              depth,
             float              lip_decay,
             float              lip_height_ratio,
             std::vector<float> shift)
{
  Array z = Array(shape);
  int   ic = (int)((0.5f - shift[0]) * z.shape[0]);
  int   jc = (int)((0.5f - shift[1]) * z.shape[1]);

  for (int i = 0; i < z.shape[0]; i++)
    for (int j = 0; j < z.shape[1]; j++)
    {
      float r = std::hypot((float)(i - ic), (float)(j - jc));

      z(i, j) = std::min(r * r / (radius * radius),
                         1.f + lip_height_ratio *
                                   std::exp(-(r - radius) / lip_decay));
      z(i, j) -= 1.f;
      z(i, j) *= depth;
    }

  return z;
}

Array peak(std::vector<int>   shape,
           float              radius,
           Array             &noise,
           float              noise_r_amp,
           float              noise_z_ratio,
           std::vector<float> shift)
{
  Array z = Array(shape);
  int   ic = (int)((0.5f - shift[0]) * z.shape[0]);
  int   jc = (int)((0.5f - shift[1]) * z.shape[1]);

  for (int i = 0; i < z.shape[0]; i++)
    for (int j = 0; j < z.shape[1]; j++)
    {
      float r = std::hypot((float)(i - ic), (float)(j - jc)) / radius;
      r += r * noise_r_amp / radius * (2 * noise(i, j) - 1);

      if (r < 1.f)
        z(i, j) = 1.f - r * r * (3.f - 2.f * r);

      z(i, j) *= 1.f + noise_z_ratio * (2.f * noise(i, j) - 1.f);
    }

  return z;
}

} // namespace hmap
