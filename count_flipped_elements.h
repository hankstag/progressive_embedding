#ifndef COUNT_FLIPPED_ELEMENTS
#define COUNT_FLIPPED_ELEMENTS

#include <Eigen/Core>

void test_flip(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F, 
  std::vector<int>& fp
);

int test_flip(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F
);

void count_flipped_element(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
);

#endif