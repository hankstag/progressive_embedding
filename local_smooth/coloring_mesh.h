#ifndef COLORING_MESH_H
#define COLORING_MESH_H

#include <Eigen/Core>

// greedy graph coloring algorithm
// [https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/]
void coloring_mesh(const Eigen::MatrixXi& F, Eigen::VectorXi& C);

#endif