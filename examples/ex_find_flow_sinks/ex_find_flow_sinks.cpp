#include "highmap.hpp"
#include "highmap/dbg.hpp"

int main(void)
{
  hmap::Vec2<int>   shape = {256, 256};
  hmap::Vec2<float> res = {4.f, 4.f};
  int               seed = 1;

  hmap::Array z = hmap::noise(hmap::NoiseType::n_simplex2, shape, res, seed);

  std::vector<int> is, js;
  hmap::find_flow_sinks(z, is, js);

  for (size_t k = 0; k < is.size(); k++)
    LOG_INFO("%d %d", is[k], js[k]);
}
