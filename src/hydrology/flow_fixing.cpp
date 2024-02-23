/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/dbg.hpp"
#include "highmap/hydrology.hpp"
#include "highmap/op.hpp"

namespace hmap
{

Array flow_fixing(const Array &z)
{
  int   ir = 8;
  int   n_dummy_nodes = 1000;
  uint  seed = 0;
  float alpha = 0.5f;

  Array z_out = Array(z.shape);

  // --- find flow sinks based on a smoothed version of the heightmap
  // --- to avoid small localized sinks

  std::vector<int> is, js; // sinks
  std::vector<int> ie, je; // exits
  {
    Array zf = z;
    smooth_cpulse(zf, ir);

    // sinks
    find_flow_sinks(zf, is, js);

    // exits: determine border cells that can be used as flow exits
    for (int j = 1; j < z.shape.y - 1; j++)
      for (int i : {0, z.shape.x - 1})
        if (zf(i, j + 1) > zf(i, j) && zf(i, j - 1) > zf(i, j))
        {
          ie.push_back(i);
          je.push_back(j);
        }

    for (int i = 1; i < z.shape.x - 1; i++)
      for (int j : {0, z.shape.y - 1})
        if (zf(i + 1, j) > zf(i, j) && zf(i - 1, j) > zf(i, j))
        {
          ie.push_back(i);
          je.push_back(j);
        }

    LOG_DEBUG("sinks");
    for (size_t k = 0; k < is.size(); k++)
      LOG_INFO("%d %d", is[k], js[k]);

    LOG_DEBUG("exits");
    for (size_t k = 0; k < ie.size(); k++)
      LOG_INFO("%d %d", ie[k], je[k]);
  }

  // create a graph using the sinks, the exits and some dummy nodes to
  // generate an arbitrary mesh connecting these nodes
  Graph graph = Graph();
  {
    Cloud cloud;

    for (size_t k = 0; k < is.size(); k++)
      cloud.add_point(Point((float)is[k], (float)js[k]));

    for (size_t k = 0; k < ie.size(); k++)
      cloud.add_point(Point((float)ie[k], (float)je[k]));

    std::vector<float> x(n_dummy_nodes);
    std::vector<float> y(n_dummy_nodes);
    Vec4<float>        bbox = {0.f, z.shape.x - 1.f, 0.f, z.shape.y - 1.f};
    random_grid_jittered(x, y, 0.4f, seed, bbox);

    // TODO random that fits the borders

    for (size_t k = 0; k < x.size(); k++)
      cloud.add_point(Point(x[k], y[k]));

    // Delanay triangulation
    graph = cloud.to_graph_delaunay();
    graph.set_values_from_array(z, bbox);
    graph.update_connectivity();

    graph.to_png("graph0.png");
  }

  // --- define flow network

  Mat<int> is_river = Mat<int>(
      Vec2<int>((int)graph.get_npoints(), (int)graph.get_npoints()));

  // compute adjacency matrix based on the Euclidian distance between
  // points...
  graph.update_adjacency_matrix();

  // ...and add elevation difference
  for (size_t i = 0; i < graph.get_npoints(); i++)
    for (size_t r = 0; r < graph.connectivity[i].size(); r++)
    {
      int j = graph.connectivity[i][r];
      if (j > (int)i)
      {
        float dz = graph.points[i].v - graph.points[j].v;
        dz = std::abs(dz) * z.shape.x;
        graph.adjacency_matrix[{i, j}] += dz + 0.5f *
                                                   std::max(graph.points[i].v,
                                                            graph.points[j].v) *
                                                   z.shape.x;
        graph.adjacency_matrix[{j, i}] = graph.adjacency_matrix[{i, j}];
      }
    }

  // for each sink, find a way out for the water
  size_t ns = is.size();
  size_t ne = ie.size();

  for (size_t k = 0; k < graph.get_npoints(); k++)
    LOG_DEBUG("%d %f %f", k, graph.points[k].x, graph.points[k].y);

  for (size_t k = 0; k < ns; k++)
  {
    // for each sinks, find the cheapest path to connect that sink to
    // one of the exits
    float            best_cost = std::numeric_limits<float>::max();
    std::vector<int> best_path = {};

    for (size_t p = 0; p < ne; p++)
    {
      std::vector<int> path = graph.dijkstra(k, ns + p);
      float            cost = 0.f;

      if (path.size() > 2)
      {
        for (size_t i = 0; i < path.size() - 1; i++)
          cost += graph.adjacency_matrix[{path[i], path[i + 1]}];

        LOG_DEBUG("in loop %ld %ld %ld %f", k, ns + p, path.size(), cost);
        LOG_DEBUG("%f %f %f %f",
                  graph.points[k].x,
                  graph.points[k].y,
                  graph.points[ns + p].x,
                  graph.points[ns + p].y);

        if (cost < best_cost)
        {
          best_path = path;
          best_cost = cost;
        }
      }
    }

    LOG_DEBUG("%ld %ld", k, best_path.size());

    // update river/non-river status
    for (size_t i = 0; i < best_path.size() - 1; i++)
    {
      int i1 = std::min(best_path[i], best_path[i + 1]);
      int i2 = std::max(best_path[i], best_path[i + 1]);
      is_river(i1, i2) += 1;
    }

    // weight adjacency matrix using elevation difference and
    // road/non-road type of the edge
    for (size_t i = 0; i < graph.get_npoints(); i++)
      for (size_t r = 0; r < graph.connectivity[i].size(); r++) // neighbors
      {
        int j = graph.connectivity[i][r];
        if (j > (int)i && is_river(i, j) == 1)
        {
          graph.adjacency_matrix[{i, j}] *= alpha;
          graph.adjacency_matrix[{j, i}] = graph.adjacency_matrix[{i, j}];
        }
      }
  }

  //--- remove orphan edges and rebuild road network graph
  Graph network = Graph(graph.get_x(), graph.get_y());

  for (size_t i = 0; i < graph.get_npoints(); i++)
    for (size_t r = 0; r < graph.connectivity[i].size(); r++)
    {
      int j = graph.connectivity[i][r];
      // LOG_DEBUG("%ld %d %d", i, j, is_river(i, j));
      if ((j > (int)i) and (is_river(i, j) > 0))
        network.add_edge({(int)i, j}, (float)is_river((int)i, j));
    }

  // // store city size in node value (equals to 0 if the node is not a
  // // city)
  // for (size_t i = 0; i < network.get_npoints(); i++)
  //   if (i < network.get_npoints() - nc)
  //     network.points[i].v = 0.f;
  //   else
  //     network.points[i].v = size[i - network.get_npoints() + nc];

  // final clean-up
  network = network.remove_orphan_points();
  network.update_adjacency_matrix();

  network.set_values(1.f);

  Array       zn = Array(z.shape);
  Vec4<float> bbox = {0.f, z.shape.x - 1.f, 0.f, z.shape.y - 1.f};
  network.to_array(zn, bbox, false);

  network.to_png("graph1.png");

  // network.print();

  return z - 0.1f * zn;
}

} // namespace hmap
