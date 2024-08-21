/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/**
 * @file colorize.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 * @version 0.1
 * @date 2023-05-08
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <png.h>

#include "highmap/array.hpp"
#include "highmap/array3.hpp"

namespace hmap
{

enum Cmap : int; // highmap/colormap.hpp

/**
 * @brief Apply hillshading to an image.
 *
 * @param img Input image.
 * @param array Elevation array.
 * @param exponent Power exponent applied to the hillshade values.
 */
void apply_hillshade(Array3            &img,
                     const hmap::Array &array,
                     float              vmin = 0.f,
                     float              vmax = 1.f,
                     float              exponent = 1.f);

/**
 * @brief Apply hillshading to an image.
 *
 * @param img Input image.
 * @param array Elevation array.
 * @param exponent Power exponent applied to the hillshade values.
 * @param is_img_rgba Whether the input image as an alpha channel or not.
 */
void apply_hillshade(std::vector<uint8_t> &img,
                     const hmap::Array    &array,
                     float                 vmin = 0.f,
                     float                 vmax = 1.f,
                     float                 exponent = 1.f,
                     bool                  is_img_rgba = false);

/**
 * @brief
 * @param array
 * @param vmin
 * @param vmax
 * @param cmap
 * @param hillshading
 * @param reverse
 * @param p_noise
 * @return
 */
Array3 colorize(Array &array,
                float  vmin,
                float  vmax,
                int    cmap,
                bool   hillshading,
                bool   reverse = false,
                Array *p_noise = nullptr);

/**
 * @brief Export array values to a 8 bit grayscale image.
 *
 * @param array Input array.
 * @return std::vector<uint8_t> Output image.
 */
Array3 colorize_grayscale(const Array &array);

/**
 * @brief Export array values to a 8 bit grayscale histogram image.
 *
 * @param array Input array.
 */
Array3 colorize_histogram(const Array &array);

/**
 * @brief Export a pair of arrays to a 8 bit colored image.
 * @param array1 Input array.
 * @param array2 Input array.
 * @return Output image.
 *
 * **Example**
 * @include ex_colorize_vec2.cpp
 *
 * **Result**
 * @image html ex_colorize_vec2.png
 */
Array3 colorize_vec2(const Array &array1, const Array &array2);

} // namespace hmap
