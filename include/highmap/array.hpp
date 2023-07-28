/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/**
 * @file array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 * @version 0.1
 * @date 2023-05-07
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace hmap
{

/**
 * @brief Array class, helper to manipulate 2D float array with "(i, j)"
 * indexing.
 *
 */
class Array
{
public:
  /**
   * @brief Array shape {ni, nj}.
   *
   */
  std::vector<int> shape;

  /**
   * @brief Vector for data storage, size shape[0] * shape[1].
   *
   */
  std::vector<float> vector;

  /**
   * @brief Construct a new Array object.
   *
   * @param shape Array shape {ni, nj}.
   * @param value Array filling value at creation.
   *
   * **Example**
   * @include ex_array.cpp
   *
   *
   */
  Array(std::vector<int> shape);

  Array(std::vector<int> shape, float value); ///< @overload

  //----------------------------------------
  // accessors
  //----------------------------------------

  /**
   * @brief Get the shape object.
   *
   * @return std::vector<int> Shape {ni, nj}.
   */
  std::vector<int> get_shape();

  /**
   * @brief Get the vector object.
   *
   * @return std::vector<float> Vector of size shape[0] * shape[1].
   */
  std::vector<float> get_vector();

  /**
   * @brief Set the array shape.
   *
   * @param new_shape New shape.
   */
  void set_shape(std::vector<int> new_shape);

  //----------------------------------------
  // overload
  //----------------------------------------

  /**
   * @brief Assignment overloading (scalar assignement).
   *
   * @param value Scalar value.
   * @return Array Reference to the current object.
   */
  Array operator=(const float value);

  /**
   * @brief Division-assignment overloading (scalar).
   *
   * @param value Scalar value.
   */
  Array operator*=(const float value);

  Array operator*=(const Array &array); ///< @overload

  /**
   * @brief Division-assignment overloading (scalar).
   *
   * @param value Scalar value.
   */
  Array operator/=(const float value);

  Array operator/=(const Array &array); ///< @overload

  /**
   * @brief Division-assignment overloading (scalar).
   *
   * @param value Scalar value.
   */
  Array operator+=(const float value);

  Array operator+=(const Array &array); ///< @overload

  /**
   * @brief Division-assignment overloading (scalar).
   *
   * @param value Scalar value.
   */
  Array operator-=(const float value);

  Array operator-=(const Array &array); ///< @overload

  /**
   * @brief Multiplication overloading (right multiply by a scalar).
   *
   * @param value Scalar value.
   * @return Array Reference to the current object.
   */
  Array operator*(const float value) const;

  /**
   * @brief Multiplication overloading (element-wise product by another array).
   *
   * @param array Another Array.
   * @return Array Reference to the resulting object.
   */
  Array operator*(const Array &array) const;

  /**
   * @brief Multiplication overloading (left multiply by a scalar).
   *
   * @param value Scalar value.
   * @param array Another Array.
   * @return Array Reference to the resulting object.
   */
  friend Array operator*(const float value, const Array &array);

  /**
   * @brief Division overloading (right divide by a scalar).
   *
   * @param value Scalar value.
   * @return Array Reference to the current object.
   */
  Array operator/(const float value) const;

  /**
   * @brief Division overloading (element-wise division by another array).
   *
   * @param array Another Array.
   * @return Array
   */
  Array operator/(const Array &array) const;

  /**
   * @brief Division overloading (left divide by a scalar).
   *
   * @param value Scalar value.
   * @param array Another Array.
   * @return Array Reference to the resulting object.
   */
  friend Array operator/(const float value, const Array &array);

  /**
   * @brief Addition overloading (right add by a scalar).
   *
   * @param value Scalar value.
   * @return Array Reference to the current object.
   */
  Array operator+(const float value) const;

  /**
   * @brief Addition overloading (element-wise addition by another array).
   *
   * @param array Another Array.
   * @return Array
   */
  Array operator+(const Array &array) const;

  /**
   * @brief Addition overloading (left add by a scalar).
   *
   * @param value Scalar value.
   * @param array Another Array.
   * @return Array Reference to the resulting object.
   */
  friend Array operator+(const float value, const Array &array);

  /**
   * @brief Unary minus overloading.
   *
   * @return Array Reference to the current object.
   */
  Array operator-() const;

  /**
   * @brief Subtraction overloading (right substract by a scalar).
   *
   * @param value Scalar value.
   * @return Array Reference to the current object.
   */
  Array operator-(const float value) const;

  /**
   * @brief Subtraction overloading (element-wise substract by another array).
   *
   * @param array Another Array.
   * @return Array
   */
  Array operator-(const Array &array) const;

  /**
   * @brief Subtraction overloading (left substract by a scalar).
   *
   * @param value Scalar value.
   * @param array Another Array.
   * @return Array Reference to the resulting object.
   */
  friend const Array operator-(const float value, const Array &array);

  /**
   * @brief Call overloading, return array value at index (i, j).
   *
   * @param i 'i' index.
   * @param j 'j' index.
   * @return float& Array value at index (i, j).
   */

  float &operator()(int i, int j)
  {
    return this->vector[i * this->shape[1] + j];
  }

  const float &operator()(int i, int j) const ///< @overload
  {
    return this->vector[i * this->shape[1] + j];
  }

  //----------------------------------------
  // methods
  //----------------------------------------

  /**
   * @brief Return a column 'j' as a std::vector.
   *
   * @param j Colunm index.
   * @return std::vector<float>
   */
  std::vector<float> col_to_vector(int j);

  /**
   * @brief Distribute a value 'amount' around the four cells (i, j), (i + 1,
   * j), (i, j + 1), (i + 1, j + 1) by "reversing" the bilinear interpolation.
   *
   * @param i Index.
   * @param j Index.
   * @param u 'u' interpolation parameter, expected to be in [0, 1[.
   * @param v 'v' interpolation parameter, expected to be in [0, 1[.
   * @param amount Amount to be deposited.
   */
  void depose_amount_bilinear_at(int i, int j, float u, float v, float amount);

  void depose_amount_kernel_bilinear_at(int   i,
                                        int   j,
                                        float u,
                                        float v,
                                        int   ir,
                                        float amount);

  /**
   * @brief Distribute a value 'amount' around the cell (i, j) using a
   * a 1D deposition kernel (applied to both direction).
   *
   * @param i Index.
   * @param j Index.
   * @param kernel Deposition kernel (1D), must have an odd number of elements.
   * @param amount Amount to be deposited.
   */
  void depose_amount_kernel_at(int i, int j, Array &kernel, float amount);

  /**
   * @brief Extract the value of a slice {i1, i2, j1, j2} to create a new array.
   *
   * @param idx Slice extent indices: {i1, i2, j1, j2}.
   * @return Array Resulting array.
   */
  Array extract_slice(std::vector<int> idx);

  /**
   * @brief Return the gradient in the 'x' (or 'i' index) of at the index (i,
   * j).
   *
   * @warning Based on a 2nd order central difference scheme, cannot be used
   * at the borders, i.e. for i = 0, j = 0, i = shape[0] - 1 or j =
   * shape[1] - 1.
   *
   * @param i Index, expected to be in [1, shape[0] - 2].
   * @param j Index, expected to be in [1, shape[1] - 2].
   * @return float
   */
  float get_gradient_x_at(int i, int j) const;

  /**
   * @brief Return the gradient in the 'y' (or 'j' index) of at the index (i,
   * j).
   *
   * @warning Based on a 2nd order central difference scheme, cannot be used
   * at the borders, i.e. for i = 0, j = 0, i = shape[0] - 1 or j =
   * shape[1] - 1.
   *
   * @param i Index, expected to be in [1, shape[0] - 2].
   * @param j Index, expected to be in [1, shape[1] - 2].
   * @return float
   */
  float get_gradient_y_at(int i, int j) const;

  /**
   * @brief Return the gradient in the 'x' (or 'i' index) of at the location (x,
   * y) near the index (i, j) using bilinear interpolation.
   *
   * @warning Based on a 2nd order central difference scheme, cannot be used
   * at the borders, i.e. for i = 0, j = 0, i = shape[0] - 1 or j =
   * shape[1] - 1.
   *
   * @param i Index, expected to be in [1, shape[0] - 2].
   * @param j Index, expected to be in [1, shape[1] - 2].
   * @param u 'u' interpolation parameter, expected to be in [0, 1[.
   * @param v 'v' interpolation parameter, expected to be in [0, 1[.
   * @return float
   */
  float get_gradient_x_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Return the gradient in the 'y' (or 'j' index) of at the location (x,
   * y) near the index (i, j) using bilinear interpolation.
   *
   * @warning Based on a 2nd order central difference scheme, cannot be used at
   * the borders, i.e. for i = 0, j = 0, i = shape[0] - 1 or j = shape[1] - 1.
   *
   * @param i Index, expected to be in [1, shape[0] - 2].
   * @param j Index, expected to be in [1, shape[1] - 2].
   * @param u 'u' interpolation parameter, expected to be in [0, 1[.
   * @param v 'v' interpolation parameter, expected to be in [0, 1[.
   * @return float
   */
  float get_gradient_y_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Return the surface normal at the index (i, j).
   *
   * @param i Index.
   * @param j Index.
   * @return std::vector<float> Normal vector (3 components).
   */
  std::vector<float> get_normal_at(int i, int j) const;

  /**
   * @brief Return the array value at the location (x, y) near the index (i, j)
   * using bilinear interpolation.
   *
   * @warning Based on bilinear interpolation, cannot be used at the upper
   * borders, i.e. for i = shape[0] - 1 or j = shape[1] - 1.
   *
   * @param i Index, expected to be in [0, shape[0] - 2].
   * @param j Index, expected to be in [0, shape[1] - 2].
   * @param u 'u' interpolation parameter, expected to be in [0, 1[.
   * @param v 'v' interpolation parameter, expected to be in [0, 1[.
   * @return float
   */
  float get_value_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Return stacked arrays in sequence horizontally (column wise).
   *
   * @param array1 1st array.
   * @param array2 2st array.
   * @return Array Reference to the resulting object.
   */
  friend Array hstack(const Array &array1, const Array &array2);

  /**
   * @brief Display a bunch of infos on the array.
   *
   */
  void infos(std::string msg = "") const;

  /**
   * @brief Return the linear index corresponding to the (i, j) cell.
   *
   * @param i 'i' index.
   * @param j 'j' index.
   * @return int Linear index.
   */
  int linear_index(int i, int j);

  /**
   * @brief Return the value of the greastest element in the array.
   *
   * @return float
   */
  float max() const;

  /**
   * @brief Return the value of the smallest element in the array.
   *
   * @return float
   */
  float min() const;

  /**
   * @brief Normalize array values so that the array sum is equal to 1.
   *
   */
  void normalize();

  /**
   * @brief Print vector values to stdout.
   *
   */
  void print();

  /**
   * @brief Return the peak-to-peak amplitude (i.e. max - min) of the array
   * values.
   *
   * @return float
   */
  float ptp() const;

  /**
   * @brief Return a resampled array of shape 'new_shape' using bilinear
   * interpolation.
   *
   * @param new_shape Target shape.
   * @return Array Resampled array.
   *
   * **Example**
   * @include ex_resample_to_shape.cpp
   *
   * **Result**
   * @image html ex_resample_to_shape.png
   */
  Array resample_to_shape(std::vector<int> new_shape);

  /**
   * @brief Return a column 'i' as a std::vector.
   *
   * @param i Row index.
   * @return std::vector<float>
   */
  std::vector<float> row_to_vector(int i);

  /**
   * @brief Set the value of a slice {i1, i2, j1, j2} of data.
   *
   * @param idx Slice extent indices: {i1, i2, j1, j2}.
   * @param value New value.
   */
  void set_slice(std::vector<int> idx, float value);

  /**
   * @brief Return the array size (number of elements).
   *
   * @return int
   */
  int size() const;

  /**
   * @brief Return of the array values.
   *
   * @return float
   */
  float sum();

  /**
   * @brief Export array a raw binary file.
   *
   * @param fname File name.
   */
  void to_file(std::string fname);

  /**
   * @brief Export array as png image file.
   *
   * @param fname File name.
   * @param cmap Colormap (@see cmap).
   * @param hillshading Activate hillshading.
   *
   * **Example**
   * @include ex_perlin.cpp
   */
  void to_png(std::string fname, int cmap, bool hillshading = false);

  /**
   * @brief Return stacked arrays in sequence vertically (row wise).
   *
   * @param array1 1st array.
   * @param array2 2st array.
   * @return Array Reference to the resulting object.
   */
  friend Array vstack(const Array &array1, const Array &array2);
};

} // namespace hmap
