/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>

#include "macrologger.h"

#include "highmap/array.hpp"

namespace hmap
{

Array read_to_array(std::string fname)
{
  cv::Mat mat = cv::imread(fname, cv::IMREAD_GRAYSCALE);

  if (mat.data == nullptr)
  {
    LOG_ERROR("error while reading the image file: %s", fname.c_str());
    return Array();
  }
  else
  {
    cv::rotate(mat, mat, cv::ROTATE_90_CLOCKWISE);
    bool remap = true;
    return cv_mat_to_array(mat, remap);
  }
}

} // namespace hmap
