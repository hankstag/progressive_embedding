#ifndef LOCAL_SMOOTH_H
#define LOCAL_SMOOTH_H

// perform local smoothing on an existing map
// specifically: use a equalateral of size of 4A as reference shape for 
//               calculation of the Jacobian for triangles with area smaller 
//               than A

#include <Eigen/Core>

void local_smoothing(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& U, // constrained vertices marked as 1
    Eigen::MatrixXd& uv,
    int loop,
    double eps,
    double avg
);


#endif