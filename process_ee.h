#ifndef PROCESS_EE_H
#define PROCESS_EE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>


void process_ee(
  const Eigen::MatrixXi& ee,
  const Eigen::MatrixXi& Fuv,
  Eigen::MatrixXi& cut
);

void buildAeq(
  const Eigen::MatrixXi& cut,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& Fuv,
  const Eigen::MatrixXd& uv,
  Eigen::SparseMatrix<double>& Aeq
);

bool read_cone_ee(
  std::string &fname, 
  Eigen::MatrixXi& RI,
  Eigen::MatrixXi& EE
);


#endif