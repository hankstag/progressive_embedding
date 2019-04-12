#ifndef FLIPPED_ELEMENTS
#define FLIPPED_ELEMENTS

#include <Eigen/Core>
#include <vector>

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

void flipped_elements(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
);

#endif