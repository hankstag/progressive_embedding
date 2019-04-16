// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Zhongshi Jiang <jiangzs@nyu.edu> and 2016 Alec Jacboson
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "dijkstra.h"
#include <igl/adjacency_list.h>

template <typename IndexType, typename DerivedD, typename DerivedP>
IGL_INLINE int igl::dijkstra(
  const IndexType &source,
  const std::set<IndexType> &targets,
  const std::vector<std::vector<IndexType> >& VV,
  const std::vector<double>& weights,
  Eigen::PlainObjectBase<DerivedD> &min_distance,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
  int numV = VV.size();
  min_distance.setConstant(numV, 1, std::numeric_limits<typename DerivedD::Scalar>::infinity());
  min_distance[source] = 0;
  previous.setConstant(numV, 1, -1);
  std::set<std::pair<typename DerivedD::Scalar, IndexType> > vertex_queue;
  vertex_queue.insert(std::make_pair(min_distance[source], source));

  while (!vertex_queue.empty())
  {
    typename DerivedD::Scalar dist = vertex_queue.begin()->first;
    IndexType u = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());

    if (targets.find(u)!= targets.end())
      return u;

    // Visit each edge exiting u
    const std::vector<int> &neighbors = VV[u];
    for (std::vector<int>::const_iterator neighbor_iter = neighbors.begin();
         neighbor_iter != neighbors.end();
         neighbor_iter++)
    {
      IndexType v = *neighbor_iter;
      typename DerivedD::Scalar distance_through_u = dist + weights[u];
      if (distance_through_u < min_distance[v]) {
        vertex_queue.erase(std::make_pair(min_distance[v], v));

        min_distance[v] = distance_through_u;
        previous[v] = u;
        vertex_queue.insert(std::make_pair(min_distance[v], v));

      }

    }
  }
  //we should never get here
  return -1;
}

template <typename IndexType, typename DerivedD, typename DerivedP>
IGL_INLINE int igl::dijkstra(
  const IndexType &source,
  const std::set<IndexType> &targets,
  const std::vector<std::vector<IndexType> >& VV,
  Eigen::PlainObjectBase<DerivedD> &min_distance,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
  std::vector<double> weights(VV.size(), 1.0);
  return dijkstra(source, targets, VV, weights, min_distance, previous);
}

template <typename IndexType, typename DerivedP>
IGL_INLINE void igl::dijkstra(
  const IndexType &vertex,
  const Eigen::MatrixBase<DerivedP> &previous,
  std::vector<IndexType> &path)
{
  IndexType source = vertex;
  path.clear();
  for ( ; source != -1; source = previous[source])
    path.push_back(source);
}


template <typename IndexType, typename DerivedV, typename DerivedF,
typename DerivedD, typename DerivedP>
IGL_INLINE void igl::dijkstra(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const std::vector<std::vector<IndexType> >& VV,
  const IndexType &source,
  const std::set<IndexType> &targets,
  Eigen::PlainObjectBase<DerivedD> &min_distance,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
  int numV = VV.size();
  if (numV ==0) igl::adjacency_list(F,VV);
  numV = VV.size();

  min_distance.setConstant(numV, 1, std::numeric_limits<typename DerivedD::Scalar>::infinity());
  min_distance[source] = 0;
  previous.setConstant(numV, 1, -1);
  std::set<std::pair<typename DerivedD::Scalar, IndexType> > vertex_queue;
  vertex_queue.insert(std::make_pair(min_distance[source], source));

  while (!vertex_queue.empty())
  {
    typename DerivedD::Scalar dist = vertex_queue.begin()->first;
    IndexType u = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());

    if (targets.find(u)!= targets.end())
      return u;

    // Visit each edge exiting u
    const std::vector<int> &neighbors = VV[u];
    for (std::vector<int>::const_iterator neighbor_iter = neighbors.begin();
         neighbor_iter != neighbors.end();
         neighbor_iter++)
    {
      IndexType v = *neighbor_iter;
      typename DerivedD::Scalar distance_through_u = dist + (V.row(u) - V.row(v)).norm();
      if (distance_through_u < min_distance[v]) {
        vertex_queue.erase(std::make_pair(min_distance[v], v));

        min_distance[v] = distance_through_u;
        previous[v] = u;
        vertex_queue.insert(std::make_pair(min_distance[v], v));

      }

    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::dijkstra<int, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, std::set<int, std::less<int>, std::allocator<int> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::dijkstra<int, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<int, std::allocator<int> >&);
#endif
