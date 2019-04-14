#ifndef LOCAL_SMOOTH_H
#define LOCAL_SMOOTH_H

// perform local smoothing on an existing map
// specifically: use a equalateral of size of 4A as reference shape for 
//               calculation of the Jacobian for triangles with area smaller 
//               than A

#include <Eigen/Core>
#include "../flipped_elements.h"

template <typename Scalar>
void grad_operator(
    const Eigen::Matrix<Scalar,3,3>& T, // reference shape
    Eigen::Matrix<Scalar,3,3>& g
);

void local_smoothing(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& M, // constraints marked as 1
    Eigen::MatrixXd& uv,
    int loop,
    double eps,
    double target_area0
);

void compute_dirichlet_energy(
  const double target_area,
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F,
  const std::vector<int>& L,
  const double eps,
  Eigen::VectorXd& E
);

#endif