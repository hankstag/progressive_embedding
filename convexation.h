#ifndef CONVEXATION_H
#define CONVEXATION_H

#include <Eigen/Core>

void convexify(
  const Eigen::Matrix<double,Eigen::Dynamic,2>& P,
  Eigen::Matrix<double,Eigen::Dynamic,2>& Pn
);

void reverse(
  Eigen::MatrixXd& uv
);

#endif