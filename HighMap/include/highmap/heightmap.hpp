/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/**
 * @file heightmap.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <functional>

#include "highmap/array.hpp"
#include "highmap/export.hpp"

namespace hmap
{

enum NormalMapBlendingMethod : int
{
  NMAP_LINEAR,
  NMAP_DERIVATIVE,
  NMAP_UDN,
  NMAP_UNITY,
  NMAP_WHITEOUT
};

static std::map<std::string, int> normal_map_blending_method_as_string = {
    {"Linear", NMAP_LINEAR},
    {"Partial derivative", NMAP_DERIVATIVE},
    {"Unreal Developer Network", NMAP_UDN},
    {"Unity", NMAP_UNITY},
    {"Whiteout", NMAP_WHITEOUT},
};

// --- forward declarations
class HeightMapRGBA;
HeightMapRGBA mix_heightmap_rgba(HeightMapRGBA &rgba1,
                                 HeightMapRGBA &rgba2,
                                 bool           use_sqrt_avg = true);
HeightMapRGBA mix_heightmap_rgba(std::vector<HeightMapRGBA *> p_rgba_list,
                                 bool use_sqrt_avg = true);

// --- classes

/**
 * @brief Tile class, to manipulate a restricted region of an heightmap (with
 * contextual informations).
 */
class Tile : public Array
{
public:
  /**
   * @brief Tile shift in each direction, assuming the global domain is a unit
   * square.
   *
   * For example, if the tiling is {4, 2}, the shift of tile {3, 2} is {0.75,
   * 0.5}.
   */
  Vec2<float> shift;

  /**
   * @brief Scale of the tile in each direction, assuming the global domain is
   * a unit square.
   *
   * For example, if the tiling is {4, 2} without overlap, the scale is {0.25,
   * 0.5}.
   */
  Vec2<float> scale;

  /**
   * @brief Tile bounding box {xmin, xmax, ymin, ymax}.
   */
  Vec4<float> bbox;

  /**
   * @brief Construct a new Tile object.
   *
   * @param shape Shape.
   * @param shift Shift.
   */
  Tile(Vec2<int> shape, Vec2<float> shift, Vec2<float> scale, Vec4<float> bbox);

  Tile(); ///< @overload

  /**
   * @brief Assignment overloading (array).
   *
   * @param array
   */
  void operator=(const Array &array);

  /**
   * @brief Fill tile values by interpolating (bilinear) values from
   * another array.
   *
   * @param array Input array.
   */
  void from_array_interp(Array &array);

  /**
   * @brief Fill tile values by interpolating (nearest) values from
   * another array.
   *
   * @param array Input array.
   */
  void from_array_interp_nearest(Array &array);

  /**
   * @brief Print some informations about the object.
   */
  void infos() const;
};

/**
 * @brief HeightMap class, to manipulate heightmap (with contextual
 * informations).
 */
class HeightMap
{
public:
  /**
   * @brief Heightmap global shape.
   */
  Vec2<int> shape;

  /**
   * @brief Heightmap bounding box {xmin, xmax, ymin, ymax}.
   */
  Vec4<float> bbox = {0.f, 1.f, 0.f, 1.f};

  /**
   * @brief Tiling setup (number of tiles in each direction).
   */
  Vec2<int> tiling = {1, 1};

  /**
   * @brief Tile overlapping, in [0, 1[.
   */
  float overlap = 0.f;

  /**
   * @brief Tile storage.
   */
  std::vector<Tile> tiles = {};

  /**
   * @brief Construct a new HeightMap object.
   *
   * @param shape Shape.
   * @param bbox Bounding box.
   * @param tiling Tiling setup.
   * @param overlap Tile overlapping.
   */
  HeightMap(Vec2<int> shape, Vec4<float> bbox, Vec2<int> tiling, float overlap);

  HeightMap(Vec2<int> shape, Vec2<int> tiling,
            float overlap); ///< @overload

  HeightMap(Vec2<int> shape,
            Vec2<int> tiling,
            float     overlap,
            float     fill_value); ///< @overload

  HeightMap(Vec2<int> shape, Vec4<float> bbox,
            Vec2<int> tiling); ///< @overload

  HeightMap(Vec2<int> shape, Vec4<float> bbox); ///< @overload

  HeightMap(Vec2<int> shape, Vec2<int> tiling); ///< @overload

  HeightMap(Vec2<int> shape); ///< @overload

  HeightMap(); ///< @overload

  //----------------------------------------
  // accessors
  //----------------------------------------

  //----------------------------------------
  // methods
  //----------------------------------------

  /**
   * @brief Get the number of tiles
   *
   * @return size_t Number of tiles.
   */
  size_t get_ntiles();

  /**
   * @brief Get the tile linear index.
   *
   * @param i Tile i index.
   * @param j Tile j index
   * @return int Linear index.
   */
  int get_tile_index(int i, int j);

  float get_value_bilinear(float x, float y);

  float get_value_nearest(float x, float y);

  void set_bbox(Vec4<float> new_bbox);

  /**
   * @brief Set the tile overlapping.
   *
   * @param new_overlap New overlap.
   */
  void set_overlap(float new_overlap);

  /**
   * @brief Set the heightmap shape.
   *
   * @param new_shape New shape.
   */
  void set_shape(Vec2<int> new_shape);

  /**
   * @brief Set the shape / tiling / overlap in one pass.
   *
   * @param new_shape New shape.
   * @param new_tiling New tiling.
   * @param new_overlap New overlap.
   */
  void set_sto(Vec2<int> new_shape, Vec2<int> new_tiling, float new_overlap);

  /**
   * @brief Set the tiling setup.
   *
   * @param new_tiling New tiling.
   */
  void set_tiling(Vec2<int> new_tiling);

  /**
   * @brief Fill tile values by interpolating (bilinear) values from
   * another array.
   *
   * @param array Input array.
   *
   * **Example**
   * @include ex_heightmap_from_array.cpp
   *
   * **Result**
   * @image html ex_heightmap_from_array0.png
   * @image html ex_heightmap_from_array1.png
   */
  void from_array_interp(Array &array);

  /**
   * @brief Fill tile values by interpolating (nearest neighbor) values from
   * another array.
   *
   * @param array Input array.
   */
  void from_array_interp_nearest(Array &array);

  /**
   * @brief Print some informations about the object.
   */
  void infos();

  /**
   * @brief Inverse the heightmap values (max - values).
   */
  void inverse();

  /**
   * @brief Return the value of the greatest element in the heightmap data.
   *
   * @return float
   */
  float max();

  /**
   * @brief Return the mean of the heightmap data.
   *
   * @return float
   */
  float mean();

  /**
   * @brief Return the value of the smallest element in the heightmap data.
   *
   * @return float
   */
  float min();

  /**
   * @brief Remap heightmap elements from a starting range to a target range.
   *
   * By default the starting range is taken to be [min(), max()] of the input
   * array.
   *
   * @param array Input array.
   * @param vmin The lower bound of the range to remap to.
   * @param vmax The lower bound of the range to remap to.
   * @param from_min The lower bound of the range to remap from.
   * @param from_max The upper bound of the range to remap from.
   */
  void remap(float vmin = 0.f, float vmax = 1.f);

  void remap(float vmin, float vmax, float from_min, float from_max);

  /**
   * @brief Smooth the transitions between each tiles (when overlap > 0).
   */
  void smooth_overlap_buffers();

  /**
   * @brief Return the sum of the heightmap data.
   *
   * @return float
   */
  float sum();

  /**
   * @brief Return the heightmap as an array.
   *
   * @param shape_export Array shape.
   * @return Array Resulting array.
   */
  Array to_array(Vec2<int> shape_export);

  Array to_array(); ///< @overload

  /**
   * @brief Converts the heightmap to a 16-bit grayscale representation.
   *
   * @return A `std::vector<uint16_t>` containing the 16-bit grayscale image
   * data.
   */
  std::vector<uint16_t> to_grayscale_image_16bit();

  std::vector<uint16_t> to_grayscale_image_16bit_multithread();

  std::vector<uint8_t> to_grayscale_image_8bit();

  /**
   * @brief Returns the unique elements of the heightmap.
   *
   * @return std::vector<float> Unique values.
   */
  std::vector<float> unique_values();

  /**
   * @brief Update tile parameters.
   */
  void update_tile_parameters();
};

/**
 * @brief HeightMap class, to manipulate a set of RGB heightmap for heightmap
 * texturing.
 */
struct HeightMapRGB
{
  /**
   * @brief RGB component heightmap storage.
   *
   * **Example**
   * @include ex_ex_heightmap_rgb.cpp
   *
   * **Result**
   * @image html ex_ex_heightmap_rgb0.png
   * @image html ex_ex_heightmap_rgb1.png
   * @image html ex_ex_heightmap_rgb2.png
   * @image html ex_ex_heightmap_rgb3.png
   */
  std::vector<HeightMap> rgb;

  /**
   * @brief Shape.
   */
  Vec2<int> shape = {0, 0};

  /**
   * @brief Constructor.
   * @param r Heightmap for R (red) component.
   * @param g Heightmap for G (green) component.
   * @param b Heightmap for B (blue) component.
   */
  HeightMapRGB(HeightMap r, HeightMap g, HeightMap b);

  HeightMapRGB(); ///< @overload

  /**
   * @brief Set the shape / tiling / overlap in one pass.
   *
   * @param new_shape New shape.
   * @param new_tiling New tiling.
   * @param new_overlap New overlap.
   */
  void set_sto(Vec2<int> new_shape, Vec2<int> new_tiling, float new_overlap);

  /**
   * @brief Fill RGB heightmap components based on a colormap and an input
   * reference heightmap.
   * @param h Input heightmap.
   * @param vmin Lower bound for scaling to array [0, 1].
   * @param vmax Upper bound for scaling to array [0, 1]
   * @param cmap Colormap (see {@link cmap}).
   * @param reverse Reverse colormap.
   */
  void colorize(HeightMap &h,
                float      vmin,
                float      vmax,
                int        cmap,
                bool       reverse = false);

  /**
   * @brief Fill RGB heightmap components based on a colormap and an input
   * reference heightmap.
   * @param h Input heightmap.
   * @param vmin Lower bound for scaling to array [0, 1].
   * @param vmax Upper bound for scaling to array [0, 1]
   * @param colormap_colors Colormap RGB colors as a vector of RGB colors.
   * @param reverse Reverse colormap.
   */
  void colorize(HeightMap                      &h,
                float                           vmin,
                float                           vmax,
                std::vector<std::vector<float>> colormap_colors,
                bool                            reverse = false);

  /**
   * @brief Normalize RGB heightmaps amplitude.
   */
  void normalize();

  /**
   * @brief Convert the RGB heightmap to a 8bit RGB image.
   * @param shape_img Resulting image shape.
   * @return Image data.
   */
  std::vector<uint8_t> to_img_8bit(Vec2<int> shape_img = {0, 0});

  /**
   * @brief Export the RGB heightmap to a 16bit png file.
   * @param fname File name.
   */
  void to_png(const std::string &fname, int depth = CV_8U);

  /**
   * @brief Mix two RGB heightmap using linear interpolation.
   * @param rgb1 1st RGB heightmap.
   * @param rgb2 2st RGB heightmap.
   * @param t Mixing parameter, in [0, 1].
   * @return RGB heightmap.
   */
  friend HeightMapRGB mix_heightmap_rgb(HeightMapRGB &rgb1,
                                        HeightMapRGB &rgb2,
                                        HeightMap    &t);

  friend HeightMapRGB mix_heightmap_rgb(HeightMapRGB &rgb1,
                                        HeightMapRGB &rgb2,
                                        float         t); ///< @overload

  /**
   * @brief Mix two RGB heightmap using weighted quadratic averaging.
   * @param rgb1 1st RGB heightmap.
   * @param rgb2 2st RGB heightmap.
   * @param t Mixing parameter, in [0, 1].
   * @return RGB heightmap.
   */
  friend HeightMapRGB mix_heightmap_rgb_sqrt(HeightMapRGB &rgb1,
                                             HeightMapRGB &rgb2,
                                             HeightMap    &t);

  friend HeightMapRGB mix_heightmap_rgb_sqrt(HeightMapRGB &rgb1,
                                             HeightMapRGB &rgb2,
                                             float         t);
};

/**
 * @brief HeightMap class, to manipulate a set of RGBA heightmap for heightmap
 * texturing.
 *
 * **Example**
 * @include ex_ex_heightmap_rgba.cpp
 *
 * **Result**
 * @image html ex_ex_heightmap_rgba0.png
 * @image html ex_ex_heightmap_rgba1.png
 * @image html ex_ex_heightmap_rgba2.png
 * @image html ex_ex_heightmap_rgba3.png
 * @image html ex_ex_heightmap_rgba4.png
 */
struct HeightMapRGBA
{
  /**
   * @brief RGBA component heightmap storage.
   */
  std::vector<HeightMap> rgba;

  /**
   * @brief Shape.
   */
  Vec2<int> shape = {0, 0};

  /**
   * @brief Constructor.
   * @param r Heightmap for R (red) component.
   * @param g Heightmap for G (green) component.
   * @param b Heightmap for B (blue) component.
   * @param a Heightmap for A (alpha) component.
   */
  HeightMapRGBA(HeightMap r, HeightMap g, HeightMap b, HeightMap a);

  HeightMapRGBA(Vec2<int> shape,
                Vec2<int> tiling,
                float     overlap,
                Array     array_r,
                Array     array_g,
                Array     array_b,
                Array     array_a);

  HeightMapRGBA(Vec2<int> shape, Vec2<int> tiling, float overlap);

  HeightMapRGBA(); ///< @overload

  /**
   * @brief Set the alpha channel.
   * @param new_alpha Alpha value(s).
   */
  void set_alpha(HeightMap new_alpha);

  void set_alpha(float new_alpha);

  /**
   * @brief Set the shape / tiling / overlap in one pass.
   *
   * @param new_shape New shape.
   * @param new_tiling New tiling.
   * @param new_overlap New overlap.
   */
  void set_sto(Vec2<int> new_shape, Vec2<int> new_tiling, float new_overlap);

  /**
   * @brief Export the RGB heightmap to a 8bit png file.
   * @param fname File name.
   */
  void to_png(const std::string &fname, int depth = CV_8U);

  /**
   * @brief Fill RGB heightmap components based on a colormap and an
   input
   * reference heightmap.
   * @param color_level Input heightmap for color level.
   * @param vmin Lower bound for scaling to array [0, 1].
   * @param vmax Upper bound for scaling to array [0, 1]
   * @param cmap Colormap (see {@link cmap}).
   * @param p_alpha Reference to input heightmap for alpha channel, expected in
   [0, 1].
   * @param reverse Reverse colormap.
   */
  void colorize(HeightMap &color_level,
                float      vmin,
                float      vmax,
                int        cmap,
                HeightMap *p_alpha = nullptr,
                bool       reverse = false,
                HeightMap *p_noise = nullptr);

  /**
   * @brief Fill RGBA heightmap components based on a colormap and
   input reference heightmaps for the color level and the transparency.
   * @param color_level Input heightmap for color level.
   * @param vmin Lower bound for scaling to array [0, 1].
   * @param vmax Upper bound for scaling to array [0, 1]
   * @param colormap_colors Colormap RGB colors as a vector of RGB
   colors.
   * @param p_alpha Reference to input heightmap for alpha channel, expected in
   [0, 1].
   * @param reverse Reverse colormap.
   */
  void colorize(HeightMap                      &color_level,
                float                           vmin,
                float                           vmax,
                std::vector<std::vector<float>> colormap_colors,
                HeightMap                      *p_alpha = nullptr,
                bool                            reverse = false,
                HeightMap                      *p_noise = nullptr);

  /**
   * @brief Computes the luminance of an RGBA height map.
   *
   * This method creates a grayscale `HeightMap` based on the luminance
   * values calculated from the red, green, and blue channels of the RGBA height
   * map. The luminance is computed using the standard formula: \f$ L = 0.299
   * \times R + 0.587 \times G + 0.114 \times B \f$.
   *
   * @return A `HeightMap` representing the grayscale luminance of the current
   * RGBA height map.
   *
   * @see https://stackoverflow.com/questions/596216 for details on the
   * luminance calculation.
   */
  HeightMap luminance();

  /**
   * @brief Mix two RGBA heightmap using alpha compositing ("over").
   * @param rgba1 1st RGBA heightmap.
   * @param rgba2 2st RGBA heightmap.
   * @param use_sqrt_avg Whether to use or not square averaging.
   * @return RGBA heightmap.
   */
  friend HeightMapRGBA mix_heightmap_rgba(HeightMapRGBA &rgba1,
                                          HeightMapRGBA &rgba2,
                                          bool           use_sqrt_avg);

  /**
   * @brief Mix two RGBA heightmap using alpha compositing ("over").
   * @param rgba_plist Heightmap reference list.
   * @param use_sqrt_avg Whether to use or not square averaging.
   * @return RGBA heightmap.
   */
  friend HeightMapRGBA mix_heightmap_rgba(
      std::vector<HeightMapRGBA *> rgba_plist,
      bool                         use_sqrt_avg);

  /**
   * @brief Normalize RGBA heightmaps amplitude.
   */
  void normalize();

  /**
   * @brief Convert the RGB heightmap to a 8bit RGB image.
   * @param shape_img Resulting image shape.
   * @return Image data.
   */
  std::vector<uint8_t> to_img_8bit(Vec2<int> shape_img = {0, 0});
};

/**
 * @brief Mixes two normal maps in RGBA format to create a blended normal map.
 *
 * This function blends a base normal map and a detail normal map using a
 * specified blending method. The detail map can be scaled to control the
 * intensity of its effect.
 *
 * @param nmap_base       Reference to the base normal map in RGBA format.
 * @param nmap_detail     Reference to the detail normal map in RGBA format.
 * @param detail_scaling  Scaling factor for the detail normal map intensity in
 * [-1.f, 1.f]. Default is 1.0f.
 * @param blending_method Method to blend the two normal maps. Options are
 * specified by the NormalMapBlendingMethod enum (e.g., NMAP_DERIVATIVE).
 *
 * @return A HeightMapRGBA object that contains the result of blending the base
 * and detail normal maps.
 *
 * **Example**
 * @include ex_mix_normal_map_rgba.cpp
 *
 * **Result**
 * @image html ex_mix_normal_map_rgba0.png
 * @image html ex_mix_normal_map_rgba1.png
 * @image html ex_mix_normal_map_rgba2.png
 * @image html ex_mix_normal_map_rgba3.png
 * @image html ex_mix_normal_map_rgba4.png
 * @image html ex_mix_normal_map_rgba5.png
 * @image html ex_mix_normal_map_rgba6.png
 */
HeightMapRGBA mix_normal_map_rgba(HeightMapRGBA          &nmap_base,
                                  HeightMapRGBA          &nmap_detail,
                                  float                   detail_scaling = 1.f,
                                  NormalMapBlendingMethod blending_method =
                                      NormalMapBlendingMethod::NMAP_DERIVATIVE);

// shape, shift, scale, noise_x, noise_y

/**
 * @brief
 *
 * @param h
 * @param p_noise_x
 * @param p_noise_y
 * @param nullary_op
 *
 * **Example**
 * @include ex_heightmap_fill.cpp
 *
 * **Result**
 * @image html ex_heightmap_fill0.png
 * @image html ex_heightmap_fill1.png
 */
void fill(HeightMap &h,
          HeightMap *p_noise_x,
          HeightMap *p_noise_y,
          std::function<
              Array(Vec2<int>, Vec4<float>, Array *p_noise_x, Array *p_noise_y)>
              nullary_op);

void fill(HeightMap                          &h,
          HeightMap                          &hin,
          HeightMap                          *p_noise_x,
          HeightMap                          *p_noise_y,
          std::function<Array(hmap::Array &,
                              Vec2<int>,
                              Vec4<float>,
                              hmap::Array *,
                              hmap::Array *)> nullary_op);

void fill(HeightMap                          &h,
          HeightMap                          *p_noise_x,
          HeightMap                          *p_noise_y,
          HeightMap                          *p_stretching,
          std::function<Array(Vec2<int>,
                              Vec4<float>,
                              hmap::Array *,
                              hmap::Array *,
                              hmap::Array *)> nullary_op);

// shape, shift, scale and noise
void fill(
    HeightMap                                                   &h,
    HeightMap                                                   *p_noise,
    std::function<Array(Vec2<int>, Vec4<float>, Array *p_noise)> nullary_op);

// shape, shift and scale
void fill(HeightMap                                   &h,
          std::function<Array(Vec2<int>, Vec4<float>)> nullary_op);

// shape only
void fill(HeightMap &h, std::function<Array(Vec2<int>)> nullary_op);

void transform(HeightMap &h, std::function<void(Array &)> unary_op);

void transform(HeightMap                                &h,
               std::function<void(Array &, Vec4<float>)> unary_op);

void transform(HeightMap                                         &h,
               HeightMap                                         *p_noise_x,
               std::function<void(Array &, Vec4<float>, Array *)> unary_op);

void transform(
    HeightMap                                                  &h,
    HeightMap                                                  *p_noise_x,
    HeightMap                                                  *p_noise_y,
    std::function<void(Array &, Vec4<float>, Array *, Array *)> unary_op);

void transform(HeightMap                                             &h,
               std::function<void(Array &, Vec2<float>, Vec2<float>)> unary_op);

// input array and mask
void transform(HeightMap                            &h,
               HeightMap                            *p_mask,
               std::function<void(Array &, Array *)> unary_op);

// for erosion
void transform(
    HeightMap       &h,
    hmap::HeightMap *p_1,
    hmap::HeightMap *p_2,
    hmap::HeightMap *p_3,
    hmap::HeightMap *p_4,
    hmap::HeightMap *p_5,
    std::function<void(Array &, Array *, Array *, Array *, Array *, Array *)>
        unary_op);

void transform(
    HeightMap                                              &h,
    hmap::HeightMap                                        *p_1,
    hmap::HeightMap                                        *p_2,
    hmap::HeightMap                                        *p_3,
    std::function<void(Array &, Array *, Array *, Array *)> unary_op);

void transform(HeightMap                                     &h,
               hmap::HeightMap                               *p_1,
               hmap::HeightMap                               *p_2,
               std::function<void(Array &, Array *, Array *)> unary_op);

void transform(HeightMap                            &h1,
               HeightMap                            &h2,
               std::function<void(Array &, Array &)> binary_op);

void transform(HeightMap                                         &h1,
               HeightMap                                         &h2,
               std::function<void(Array &, Array &, Vec4<float>)> binary_op);

void transform(HeightMap                                     &h1,
               HeightMap                                     &h2,
               HeightMap                                     &h3,
               std::function<void(Array &, Array &, Array &)> ternary_op);

void transform(
    HeightMap                                                  &h1,
    HeightMap                                                  &h2,
    HeightMap                                                  &h3,
    std::function<void(Array &, Array &, Array &, Vec4<float>)> ternary_op);

// with returned array
void transform(HeightMap                    &h_out, // output
               HeightMap                    &h1,    // in 1
               std::function<Array(Array &)> unary_op);

void transform(HeightMap                             &h_out, // output
               HeightMap                             &h1,    // in 1
               HeightMap                             &h2,    // in 2
               std::function<Array(Array &, Array &)> binary_op);

void transform(
    HeightMap                                              &h1,
    HeightMap                                              &h2,
    HeightMap                                              &h3,
    HeightMap                                              &h4,
    std::function<void(Array &, Array &, Array &, Array &)> ternary_op);

void transform(
    HeightMap &h1,
    HeightMap &h2,
    HeightMap &h3,
    HeightMap &h4,
    HeightMap &h5,
    HeightMap &h6,
    std::function<void(Array &, Array &, Array &, Array &, Array &, Array &)>
        op);

} // namespace hmap
