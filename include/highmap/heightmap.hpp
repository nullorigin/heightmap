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

namespace hmap
{

class Tile : public Array
{
public:
  std::vector<float> shift;
  std::vector<float> bbox;

  Tile();

  Tile(std::vector<int> shape, std::vector<float> shift);

  void operator=(const Array &array);

  void infos() const;
};

/**
 * @brief HeightMap class, to manipulate heightmap (with contextual
 * informations).
 *
 */
class HeightMap
{
public:
  std::vector<int>   shape;
  std::vector<float> bbox = {0.f, 1.f, 0.f, 1.f};
  std::vector<int>   tiling = {1, 1};
  std::vector<float> tile_scale = {1.f, 1.f};
  int                ntiles;
  std::vector<Tile>  tiles = {};
  std::vector<int>   shape_tile = {};

  HeightMap(std::vector<int>   shape,
            std::vector<float> bbox,
            std::vector<int>   tiling);

  HeightMap(std::vector<int> shape, std::vector<float> bbox); /// @overload

  HeightMap(std::vector<int> shape, std::vector<int> tiling); /// @overload

  HeightMap(std::vector<int> shape); /// @overload

  HeightMap(); /// @overload

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

  /**
   * @brief Set the array shape.
   *
   * @param new_shape New shape.
   */
  void set_shape(std::vector<int> new_shape);

  /**
   * @brief Print some informations about the object.
   *
   */
  void infos();

  /**
   * @brief Return the value of the smallest element in the heightmap data.
   *
   * @return float
   */
  float min();

  /**
   * @brief Return the value of the greatest element in the heightmap data.
   *
   * @return float
   */
  float max();

  /**
   * @brief Remap heightmap values.
   *
   * @param vmin The lower bound of the range to remap to.
   * @param vmax The lower bound of the range to remap to.
   */
  void remap(float vmin = 0.f, float vmax = 1.f);

  /**
   * @brief Return the heightmap as an array.
   *
   * @param shape_export Array shape.
   * @return Array Resulting array.
   */
  Array to_array(std::vector<int> shape_export);

  Array to_array(); // @overload

  /**
   * @brief Update tile parameters.
   *
   */
  void update_tile_parameters();

private:
};

void fill(
    HeightMap &h,
    HeightMap &noise,
    std::function<Array(std::vector<int>, std::vector<float>, Array *p_noise)>
        nullary_op); // shape and shift and noise

void fill(HeightMap &h,
          std::function<Array(std::vector<int>, std::vector<float>)>
              nullary_op); // shape and shift

void fill(HeightMap                             &h,
          std::function<Array(std::vector<int>)> nullary_op); // shape only

void transform(HeightMap &h, std::function<void(Array &)> unary_op);

void transform(HeightMap                            &h1,
               HeightMap                            &h2,
               std::function<void(Array &, Array &)> binary_op);

void transform(HeightMap                            &h1,
               HeightMap                            &h2,
               HeightMap                            &h3,
               std::function<void(Array &, Array &, Array &)> ternary_op);

} // namespace hmap