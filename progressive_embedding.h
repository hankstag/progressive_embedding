#ifndef PROGRESSIVE_EMBEDDING
#define PROGRESSIVE_EMBEDDING

#include <Eigen/Core>

bool progressive_embedding(
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  const Eigen::VectorXi& bi,
  const Eigen::MatrixXd& b,
  double eps
);

bool progressive_fix(
    const Eigen::VectorXi& cs,
    Eigen::VectorXi& bi,
    Eigen::MatrixXd& b,
    const Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& uv
);

int expand_to_boundary(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& B,
  const Eigen::MatrixXi& dEF_s,
  const int fid,
  Eigen::VectorXi& local_F
);

void check_result(
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F,
  const Eigen::Matrix<double,2,3>& G
);

#endif