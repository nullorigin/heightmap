#include <cmath>
#include <random>

#include "Interpolate.hpp"
#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/op.hpp"
#include "highmap/primitives.hpp"

#define NSIGMA 2

namespace hmap
{

void expand(Array &array, int ir)
{
  Array array_new = array;
  int   ni = array.shape[0];
  int   nj = array.shape[1];
  Array k = cubic_pulse({2 * ir + 1, 2 * ir + 1});

  for (int i = 0; i < ni; i++)
  {
    int p1 = std::max(0, i - ir) - i;
    int p2 = std::min(ni, i + ir + 1) - i;
    for (int j = 0; j < nj; j++)
    {
      int q1 = std::max(0, j - ir) - j;
      int q2 = std::min(nj, j + ir + 1) - j;
      for (int p = p1; p < p2; p++)
        for (int q = q1; q < q2; q++)
        {
          float v = array(i + p, j + q) * k(p + ir, q + ir);
          if (v > array_new(i, j))
            array_new(i, j) = v;
        }
    }
  }

  array = array_new;
}

void gain(Array &array, float gain)
{
  auto lambda = [&gain](float x)
  {
    return x < 0.5 ? 0.5f * std::pow(2.f * x, gain)
                   : 1.f - 0.5f * std::pow(2.f * (1.f - x), gain);
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void gamma_correction(Array &array, float gamma)
{
  auto lambda = [&gamma](float x) { return std::pow(x, gamma); };
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void gamma_correction_local(Array &array, float gamma, int ir, float k)
{
  Array amin = minimum_local(array, ir);
  Array amax = maximum_local(array, ir);

  if (k != 0) // with smoothing
  {
    for (int i = 0; i < array.shape[0]; i++)
      for (int j = 0; j < array.shape[1]; j++)
      {
        float v = (array(i, j) - amin(i, j)) / (amax(i, j) - amin(i, j));
        v = std::sqrt(v * v + k);
        array(i, j) =
            std::pow(v, gamma) * (amax(i, j) - amin(i, j)) + amin(i, j);
      }
  }
  else // without smoothing
  {
    for (int i = 0; i < array.shape[0]; i++)
      for (int j = 0; j < array.shape[1]; j++)
      {
        float v = (array(i, j) - amin(i, j)) / (amax(i, j) - amin(i, j));
        array(i, j) =
            std::pow(v, gamma) * (amax(i, j) - amin(i, j)) + amin(i, j);
      }
  }
}

void laplace(Array &array, float sigma, int iterations)
{
  for (int it = 0; it < iterations; it++)
  {
    Array delta = laplacian(array);
    array += sigma * delta;
  }
}

void laplace_edge_preserving(Array &array,
                             float  talus,
                             float  sigma,
                             int    iterations)
{
  for (int it = 0; it < iterations; it++)
  {
    Array c = gradient_norm(array);
    c = 1.f / (1.f + c * c / (talus * talus));

    Array dcx = gradient_x(c);
    Array dcy = gradient_y(c);
    Array dzx = gradient_x(array);
    Array dzy = gradient_y(array);
    Array delta = laplacian(array);

    array += sigma * (dcx * dzx + dcy * dzy + c * delta);
  }
}

void low_pass_high_order(Array &array, int order, float sigma)
{
  Array df = array;

  // filtering coefficients
  std::vector<float> kernel;

  switch (order)
  {
  case (5):
    kernel = {0.0625f, -0.25f, 0.375f, -0.25f, 0.0625f};
    break;

  case (7):
    kernel = {-0.015625f,
              0.09375f,
              -0.234375f,
              0.3125f,
              -0.234375f,
              0.09375f,
              -0.015625f};
    break;

  case (9):
    kernel = {0.00390625f,
              -0.03125f,
              0.109375f,
              -0.21875f,
              0.2734375f,
              -0.21875f,
              0.109375f,
              -0.03125f,
              0.00390625f};
    break;
  }

  df = convolve1d_i(df, kernel);
  df = convolve1d_j(df, kernel);

  array = array - sigma * df;
}

void make_binary(Array &array, float threshold)
{
  auto lambda = [&threshold](float a)
  { return std::abs(a) > threshold ? 1.f : 0.f; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recast_canyon(Array &array, const Array &vcut, float gamma)
{
  auto lambda = [&gamma](float a, float b)
  { return a > b ? a : b * std::pow(a / b, gamma); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vcut.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

void recast_canyon(Array &array, float vcut, float gamma)
{
  auto lambda = [&vcut, &gamma](float a)
  { return a > vcut ? a : vcut * std::pow(a / vcut, gamma); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recast_peak(Array &array, int ir, float gamma, float k)
{
  Array ac = array;
  smooth_cpulse(ac, ir);
  array = maximum_smooth(array, ac, k);
  clamp_min(array, 0.f);
  array = ac * pow(array, gamma);
}

void recast_rocky_slopes(Array &array,
                         float  talus,
                         int    ir,
                         float  amplitude,
                         uint   seed,
                         float  kw,
                         float  gamma,
                         Array *p_noise)
{
  // slope-based criteria
  Array c = select_gradient_binary(array, talus);
  smooth_cpulse(c, ir);

  if (!p_noise)
  {
    Array noise = fbm_perlin(array.shape, {kw, kw}, seed, 4, 0.f);
    gamma_correction_local(noise, gamma, ir);
    {
      int ir2 = (int)(ir / 4.f);
      if (ir2 > 1)
        gamma_correction_local(noise, gamma, ir2);
    }
    array += amplitude * noise * c;
  }
  else
    array += amplitude * (*p_noise) * c;
}

void recurve(Array                    &array,
             const std::vector<float> &t,
             const std::vector<float> &v)
{
  _1D::MonotonicInterpolator<float> interp;
  interp.setData(t, v);

  auto lambda = [&interp](float a) { return interp(a); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recurve_bexp(Array &array, float tau)
{
  float c = -1.f / tau;
  auto  lambda = [&c](float a) { return 1.f - std::exp(c * a); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recurve_exp(Array &array, float tau)
{
  float c = -1.f / tau;
  auto  lambda = [&c](float a) { return std::exp(c * (1.f - a)); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recurve_s(Array &array)
{
  auto lambda = [](float a) { return a * a * (3.f - 2.f * a); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void recurve_smoothstep_rational(Array &array, float n)
{
  auto lambda = [&n](float a)
  {
    float an = std::pow(a, n);
    return an / (an + std::pow(1.f - a, n));
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void sharpen(Array &array, float ratio)
{
  Array lp = Array(array.shape);

  for (int i = 1; i < array.shape[0] - 1; i++)
    for (int j = 1; j < array.shape[1] - 1; j++)
    {
      lp(i, j) = 5.f * array(i, j) - array(i + 1, j) - array(i - 1, j) -
                 array(i, j - 1) - array(i, j + 1);
    }
  extrapolate_borders(lp);
  array = (1.f - ratio) * array + ratio * lp;
}

void shrink(Array &array, int ir)
{
  float amax = array.max();
  array = amax - array;
  expand(array, ir);
  array = amax - array;
}

void smooth_cpulse(Array &array, int ir)
{
  // define kernel
  const int          nk = 2 * ir + 1;
  std::vector<float> k(nk);

  float sum = 0.f;
  float x0 = (float)nk / 2.f;
  for (int i = 0; i < nk; i++)
  {
    float x = std::abs((float)i - x0) / (float)ir;
    k[i] = 1.f - x * x * (3.f - 2.f * x);
    sum += k[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k[i] /= sum;
  }

  // eventually convolve
  array = convolve1d_i(array, k);
  array = convolve1d_j(array, k);
}

void smooth_gaussian(Array &array, int ir)
{
  // define Gaussian kernel (we keep NSIGMA standard deviations of the
  // kernel support)
  const int          nk = NSIGMA * (2 * ir + 1);
  std::vector<float> k(nk);

  float sum = 0.f;
  float sig2 = (float)(ir * ir);
  float x0 = (float)nk / 2.f;
  for (int i = 0; i < nk; i++)
  {
    float x = (float)i - x0;
    k[i] = std::exp(-0.5f * std::pow(x, 2.f) / sig2);
    sum += k[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k[i] /= sum;
  }

  // eventually convolve
  array = convolve1d_i(array, k);
  array = convolve1d_j(array, k);
}

void smooth_fill(Array &array, int ir, float k)
{
  Array array_smooth = array;
  smooth_cpulse(array_smooth, ir);
  array = maximum_smooth(array, array_smooth, k);
}

void smooth_fill_holes(Array &array, int ir)
{
  Array array_smooth = mean_local(array, ir);

  // mask based on concave regions
  Array mask = curvature_mean(array_smooth);
  clamp_min(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 1)
    smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void smooth_fill_smear_peaks(Array &array, int ir)
{
  Array array_smooth = mean_local(array, ir);

  // mask based on concave regions
  Array mask = curvature_mean(array_smooth);
  clamp_max(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 0)
    smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void steepen(Array &array, float scale, int ir)
{
  Array dx = gradient_x(array) * ((float)array.shape[0] * -scale);
  Array dy = gradient_y(array) * ((float)array.shape[1] * -scale);

  smooth_cpulse(dx, ir);
  smooth_cpulse(dy, ir);

  warp(array, dx, dy);
}

void steepen_convective(Array &array,
                        float  angle,
                        int    iterations,
                        int    ir,
                        float  dt)
{
  Array dx = Array(array.shape);
  Array dy = Array(array.shape);
  float alpha = angle / 180.f * M_PI;
  float ca = std::cos(alpha);
  float sa = std::sin(alpha);

  for (int it = 0; it < iterations; it++)
  {
    if (ir > 0)
    {
      Array array_filtered = array;
      smooth_cpulse(array_filtered, ir);
      gradient_x(array_filtered, dx);
      gradient_y(array_filtered, dy);
    }
    else
    {
      gradient_x(array, dx);
      gradient_y(array, dy);
    }
    array *= 1.f - dt * (ca * dx + sa * dy); // == du / dt = - u * du / dx
  }
}

} // namespace hmap
