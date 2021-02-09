#ifndef EDGE_SPLIT
#define EDGE_SPLIT

#include <Eigen/Core>
#include <vector>
#include "fData.h"

bool edge_split(
    std::vector<std::vector<double>>& V0,
    std::vector<std::vector<int>>& F0,
    std::vector<std::vector<int>>& FF_vec,
    std::vector<std::vector<int>>& FFi_vec,
    int f0,
    int e0
);

bool edge_split(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXi& FF,
    Eigen::MatrixXi& FFi,
    int f0,
    int e0
);

bool edge_split(
  std::vector<std::vector<double>>& V0,
  std::vector<std::vector<int>>& F0,
  std::vector<std::vector<int>>& FF_vec,
  std::vector<std::vector<int>>& FFi_vec,
  std::vector<fData>& props,
  int f0,
  int e0
);

#endif