#include "highmap/array.hpp"
#include "highmap/geometry.hpp"
#include "highmap/io.hpp"

int main(void)
{
  int                seed = 1;
  std::vector<float> bbox = {-1.f, 2.f, 0.f, 5.f};

  // create a cloud of points and convert it to a graph using Delaunay
  // triangulation
  int         npoints = 15;
  hmap::Cloud cloud = hmap::Cloud(npoints, seed, bbox);
  hmap::Graph graph = cloud.to_graph_delaunay();

  graph.print();

  graph.update_adjacency_matrix();
  graph.update_connectivity();

  int              i_point_start = 0;
  int              i_point_end = 4;
  std::vector<int> route = graph.dijkstra(i_point_start, i_point_end);

  std::cout << "Route (point indices):\n";
  for (auto &p : route)
    std::cout << p << "\n";
}
