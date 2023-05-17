#include <cmath>

#include "highmap/array.hpp"

namespace hmap
{

Array atan(const Array &array)
{
  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::atan(v); });
  return array_out;
}

Array cos(const Array &array)
{
  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::cos(v); });
  return array_out;
}

Array hypot(Array &array1, Array &array2)
{
  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return std::hypot(a, b); });
  return array_out;
}

Array lerp(Array &array1, Array &array2, Array &t)
{
  return array1 * (1.f - t) + array2 * t;
}

Array pow(const Array &array, float exp)
{
  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&exp](float v) { return std::pow(v, exp); });
  return array_out;
}

Array sin(const Array &array)
{
  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::sin(v); });
  return array_out;
}

} // namespace hmap