#ifndef DECOMPOSE_POLYGON
#define DECOMPOSE_POLYGON

#include "shor.h"

void decompose_polygon(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F, 
  std::vector<std::vector<int>>& L
);

#endif