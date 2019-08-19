#ifndef VALIDITY_CHECK_H
#define VALIDITY_CHECK_H

#include <Eigen/Core>

bool is_face_flipped(const Eigen::Matrix<double,3,2>& T);

void grad_to_eqtri(
  double target_area,
  Eigen::Matrix<double,2,3>& G_t
);

// given a face in 3d, return the gradient operator for any double function defined on it
void grad_operator(
    const Eigen::Matrix3d& T, // reference shape
    Eigen::Matrix3d& g
);

void flipped_elements(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
);

#endif
