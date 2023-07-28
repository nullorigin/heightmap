/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cmath>

#include "highmap/array.hpp"

namespace hmap
{

void extrapolate_borders(Array &array, int nbuffer)
{
  const int ni = array.shape[0];
  const int nj = array.shape[1];

  for (int j = 0; j < nj; j++)
    for (int k = nbuffer - 1; k > -1; k--)
    {
      array(k, j) = 2.f * array(k + 1, j) - array(k + 2, j);
      array(ni - 1 - k, j) = 2.f * array(ni - 2 - k, j) - array(ni - 3 - k, j);
    }

  for (int i = 0; i < ni; i++)
    for (int k = nbuffer - 1; k > -1; k--)
    {
      array(i, k) = 2.f * array(i, k + 1) - array(i, k + 2);
      array(i, nj - 1 - k) = 2.f * array(i, nj - 2 - k) - array(i, nj - 3 - k);
    }
}

void fill_borders(Array &array)
{
  const int ni = array.shape[0];
  const int nj = array.shape[1];

  for (int j = 0; j < nj; j++)
  {
    array(0, j) = array(1, j);
    array(ni - 1, j) = array(ni - 2, j);
  }

  for (int i = 0; i < ni; i++)
  {
    array(i, 0) = array(i, 1);
    array(i, nj - 1) = array(i, nj - 2);
  }
}

Array generate_buffered_array(const Array &array, std::vector<int> buffers)
{
  Array array_out = Array({array.shape[0] + buffers[0] + buffers[1],
                           array.shape[1] + buffers[2] + buffers[3]});

  for (int i = 0; i < array.shape[0]; i++)
    for (int j = 0; j < array.shape[1]; j++)
      array_out(i + buffers[0], j + buffers[2]) = array(i, j);

  int i1 = buffers[0];
  int i2 = buffers[1];
  int j1 = buffers[2];
  int j2 = buffers[3];

  for (int i = 0; i < i1; i++)
    for (int j = j1; j < array_out.shape[1] - j2; j++)
      array_out(i, j) = array_out(2 * i1 - i, j);

  for (int i = array_out.shape[0] - i2; i < array_out.shape[0]; i++)
    for (int j = j1; j < array_out.shape[1] - j2; j++)
      array_out(i, j) = array_out(2 * (array_out.shape[0] - i2) - i - 1, j);

  for (int i = 0; i < array_out.shape[0]; i++)
    for (int j = 0; j < j1; j++)
      array_out(i, j) = array_out(i, 2 * j1 - j);

  for (int i = 0; i < array_out.shape[0]; i++)
    for (int j = array_out.shape[1] - j2; j < array_out.shape[1]; j++)
      array_out(i, j) = array_out(i, 2 * (array_out.shape[1] - j2) - j - 1);

  return array_out;
}

void make_periodic(Array &array, int nbuffer)
{
  int ni = array.shape[0];
  int nj = array.shape[1];

  Array a1 = array;
  for (int i = 0; i < nbuffer; i++)
  {
    float r = 0.5f * (float)i / ((float)nbuffer - 1.f);
    for (int j = 0; j < nj; j++)
    {
      a1(i, j) = (0.5f + r) * array(i, j) + (0.5f - r) * array(ni - 1 - i, j);
      a1(ni - 1 - i, j) = (0.5f + r) * array(ni - 1 - i, j) +
                          (0.5f - r) * array(i, j);
    }
  }

  Array a2 = a1;
  for (int j = 0; j < nbuffer; j++)
  {
    float r = 0.5f * (float)j / ((float)nbuffer - 1);
    for (int i = 0; i < ni; i++)
    {
      a2(i, j) = (0.5 + r) * a1(i, j) + (0.5 - r) * a1(i, nj - 1 - j);
      a2(i, nj - 1 - j) = (0.5 + r) * a1(i, nj - 1 - j) + (0.5 - r) * a1(i, j);
    }
  }

  array = a2;
}

void set_borders(Array             &array,
                 std::vector<float> border_values,
                 std::vector<int>   buffer_sizes)
{
  // west
  for (int i = 0; i < buffer_sizes[0]; i++)
    for (int j = 0; j < array.shape[1]; j++)
    {
      float r = (float)i / (float)buffer_sizes[0];
      r = r * r * (3.f - 2.f * r);
      array(i, j) = (1.f - r) * border_values[0] + r * array(i, j);
    }

  // east
  for (int i = array.shape[0] - buffer_sizes[1]; i < array.shape[0]; i++)
    for (int j = 0; j < array.shape[1]; j++)
    {
      float r = 1.f - (float)(i - array.shape[0] + buffer_sizes[1]) /
                          (float)buffer_sizes[1];
      r = r * r * (3.f - 2.f * r);
      array(i, j) = (1.f - r) * border_values[1] + r * array(i, j);
    }

  // south
  for (int i = 0; i < array.shape[0]; i++)
    for (int j = 0; j < buffer_sizes[2]; j++)
    {
      float r = (float)j / (float)buffer_sizes[2];
      r = r * r * (3.f - 2.f * r);
      array(i, j) = (1.f - r) * border_values[2] + r * array(i, j);
    }

  // north
  for (int i = 0; i < array.shape[0]; i++)
    for (int j = array.shape[1] - buffer_sizes[3]; j < array.shape[1]; j++)
    {
      float r = 1.f - (float)(j - array.shape[1] + buffer_sizes[3]) /
                          (float)buffer_sizes[3];
      r = r * r * (3.f - 2.f * r);
      array(i, j) = (1.f - r) * border_values[3] + r * array(i, j);
    }
}

void set_borders(Array &array, float border_values, int buffer_sizes)
{
  std::vector<float> bv = {border_values,
                           border_values,
                           border_values,
                           border_values};
  std::vector<int>   bs = {buffer_sizes,
                           buffer_sizes,
                           buffer_sizes,
                           buffer_sizes};
  set_borders(array, bv, bs);
}

void sym_borders(Array &array, std::vector<int> buffer_sizes)
{
  const int i1 = buffer_sizes[0];
  const int i2 = buffer_sizes[1];
  const int j1 = buffer_sizes[2];
  const int j2 = buffer_sizes[3];

  // fill-in the blanks...
  for (int i = 0; i < i1; i++)
  {
    for (int j = j1; j < array.shape[1] - j2; j++)
    {
      array(i, j) = array(2 * i1 - i, j);
    }
  }

  for (int i = array.shape[0] - i2; i < array.shape[0]; i++)
  {
    for (int j = j1; j < array.shape[1] - j2; j++)
    {
      array(i, j) = array(2 * (array.shape[0] - i2) - i - 1, j);
    }
  }

  for (int i = 0; i < array.shape[0]; i++)
  {
    for (int j = 0; j < j1; j++)
    {
      array(i, j) = array(i, 2 * j1 - j);
    }
  }

  for (int i = 0; i < array.shape[0]; i++)
  {
    for (int j = array.shape[1] - j2; j < array.shape[1]; j++)
    {
      array(i, j) = array(i, 2 * (array.shape[1] - j2) - j - 1);
    }
  }
}

void zeroed_borders(Array &array)
{
  const int ni = array.shape[0];
  const int nj = array.shape[1];

  for (int j = 0; j < nj; j++)
  {
    array(0, j) = 0.f;
    array(ni - 1, j) = 0.f;
  }

  for (int i = 0; i < ni; i++)
  {
    array(i, 0) = 0.f;
    array(i, nj - 1) = 0.f;
  }
}

} // namespace hmap
