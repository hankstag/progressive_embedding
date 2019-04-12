#ifndef ENERGY_H
#define ENERGY_H

#include <Eigen/Core>

void grad_operator(
    const Eigen::Matrix<double,3,3>& T, // reference shape
    Eigen::Matrix<double,3,3>& g
);

void gradient_equilateral(
  double target_area,
  Eigen::Matrix<double,2,3>& G_t
);

#endif