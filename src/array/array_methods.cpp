/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <vector>

#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/op.hpp"
#include "highmap/vector.hpp"

namespace hmap

{

std::vector<float> Array::col_to_vector(int j)
{
  std::vector<float> vec(this->shape.x);
  for (int i = 0; i < this->shape.x; i++)
    vec[i] = (*this)(i, j);
  return vec;
}

void Array::depose_amount_bilinear_at(int   i,
                                      int   j,
                                      float u,
                                      float v,
                                      float amount)
{
  (*this)(i, j) += amount * (1 - u) * (1 - v);
  (*this)(i + 1, j) += amount * u * (1 - v);
  (*this)(i, j + 1) += amount * (1 - u) * v;
  (*this)(i + 1, j + 1) += amount * u * v;
}

void Array::depose_amount_kernel_bilinear_at(int   i,
                                             int   j,
                                             float u,
                                             float v,
                                             int   ir,
                                             float amount)
{
  Array kernel = Array({2 * ir + 1, 2 * ir + 1});

  // compute kernel first
  for (int p = -ir; p < ir + 1; p++)
  {
    for (int q = -ir; q < ir + 1; q++)
    {
      float x = (float)p - u;
      float y = (float)q - v;
      float r = std::max(0.f, 1.f - std::hypot(x, y));
      kernel(p + ir, q + ir) = r;
    }
  }
  kernel.normalize();

  // perform deposition
  this->depose_amount_kernel_at(i, j, kernel, amount);
}

void Array::depose_amount_kernel_at(int i, int j, Array &kernel, float amount)
{
  const int ir = (kernel.shape.x - 1) / 2;
  const int jr = (kernel.shape.y - 1) / 2;

  for (int p = 0; p < kernel.shape.x; p++)
  {
    for (int q = 0; q < kernel.shape.y; q++)
    {
      (*this)(i + p - ir, j + q - jr) += amount * kernel(p, q);
    }
  }
}

Array Array::extract_slice(Vec4<int> idx)
{
  Array array_out = Array({idx.b - idx.a, idx.d - idx.c});

  for (int i = idx.a; i < idx.b; i++)
    for (int j = idx.c; j < idx.d; j++)
      array_out(i - idx.a, j - idx.c) = (*this)(i, j);

  return array_out;
}

float Array::get_gradient_x_at(int i, int j) const
{
  return 0.5f * ((*this)(i + 1, j) - (*this)(i - 1, j));
}

float Array::get_gradient_y_at(int i, int j) const
{
  return 0.5f * ((*this)(i, j + 1) - (*this)(i, j - 1));
}

float Array::get_gradient_x_bilinear_at(int i, int j, float u, float v) const
{
  float f00 = (*this)(i, j) - (*this)(i - 1, j);
  float f10 = (*this)(i + 1, j) - (*this)(i, j);
  float f01 = (*this)(i, j + 1) - (*this)(i - 1, j + 1);
  float f11 = (*this)(i + 1, j + 1) - (*this)(i, j + 1);

  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

float Array::get_gradient_y_bilinear_at(int i, int j, float u, float v) const
{
  float f00 = (*this)(i, j) - (*this)(i, j - 1);
  float f10 = (*this)(i + 1, j) - (*this)(i + 1, j - 1);
  float f01 = (*this)(i, j + 1) - (*this)(i, j);
  float f11 = (*this)(i + 1, j + 1) - (*this)(i + 1, j);

  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

Vec3<float> Array::get_normal_at(int i, int j) const
{
  Vec3<float> normal;

  normal.x = -this->get_gradient_x_at(i, j);
  normal.y = -this->get_gradient_y_at(i, j);
  normal.z = 1.f;

  float norm = std::hypot(normal.x, normal.y, normal.z);
  normal /= norm;

  return normal;
}

float Array::get_value_bilinear_at(int i, int j, float u, float v) const
{
  float a10 = (*this)(i + 1, j) - (*this)(i, j);
  float a01 = (*this)(i, j + 1) - (*this)(i, j);
  float a11 = (*this)(i + 1, j + 1) - (*this)(i + 1, j) - (*this)(i, j + 1) +
              (*this)(i, j);

  return (*this)(i, j) + a10 * u + a01 * v + a11 * u * v;
}

Array hstack(const Array &array1, const Array &array2) // friend function
{
  Array array_out = Array({array1.shape.x + array2.shape.x, array1.shape.y});

  for (int i = 0; i < array1.shape.x; i++)
    for (int j = 0; j < array1.shape.y; j++)
      array_out(i, j) = array1(i, j);

  for (int i = 0; i < array2.shape.x; i++)
    for (int j = 0; j < array1.shape.y; j++)
      array_out(i + array1.shape.x, j) = array2(i, j);

  return array_out;
}

int Array::linear_index(int i, int j)
{
  return i * this->shape.y + j;
}

float Array::max() const
{
  return *std::max_element(this->vector.begin(), this->vector.end());
}

float Array::min() const
{
  return *std::min_element(this->vector.begin(), this->vector.end());
};

void Array::normalize()
{
  float sum = this->sum();

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 this->vector.begin(),
                 [&sum](float v) { return v / sum; });
}

float Array::ptp() const
{
  return this->max() - this->min();
}

Array Array::resample_to_shape(Vec2<int> new_shape)
{
  Array array_out = Array(new_shape);

  // interpolation grid scaled to the starting grid to ease seeking of
  // the reference (i, j) indices during bilinear interpolation
  std::vector<float> x = linspace(0.f, (float)this->shape.x - 1, new_shape.x);
  std::vector<float> y = linspace(0.f, (float)this->shape.y - 1, new_shape.y);

  for (int i = 0; i < new_shape.x; i++)
  {
    int iref = (int)x[i];
    for (int j = 0; j < new_shape.y; j++)
    {
      int   jref = (int)y[j];
      float u = x[i] - (float)iref;
      float v = y[j] - (float)jref;
      array_out(i, j) = this->get_value_bilinear_at(iref, jref, u, v);
    }
  }

  return array_out;
}

std::vector<float> Array::row_to_vector(int i)
{
  std::vector<float> vec(this->shape.y);
  for (int j = 0; j < this->shape.y; j++)
    vec[j] = (*this)(i, j);
  return vec;
}

void Array::set_slice(Vec4<int> idx, float value)
{
  for (int i = idx.a; i < idx.b; i++)
    for (int j = idx.c; j < idx.d; j++)
      (*this)(i, j) = value;
}

int Array::size() const
{
  return this->shape.x * this->shape.y;
}

float Array::sum()
{
  return std::accumulate(this->vector.begin(), this->vector.end(), 0.f);
}

Array vstack(const Array &array1, const Array &array2) // friend function
{
  Array array_out = Array({array1.shape.x, array1.shape.y + array2.shape.y});

  for (int i = 0; i < array1.shape.x; i++)
  {
    for (int j = 0; j < array1.shape.y; j++)
      array_out(i, j) = array1(i, j);

    for (int j = 0; j < array2.shape.y; j++)
      array_out(i, j + array1.shape.y) = array2(i, j);
  }

  return array_out;
}

} // namespace hmap
