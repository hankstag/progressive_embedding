#ifndef EDGE_SPLIT
#define EDGE_SPLIT

#include <Eigen/Core>
#include <vector>
#include "field/field.h"


bool edge_split(
  std::vector<std::vector<double>>& V,
  std::vector<std::vector<int>>& F,
  std::vector<std::vector<int>>& FF_vec,
  std::vector<std::vector<int>>& FFi_vec,
  std::vector<Property>& props,
  int f0,
  int e0
);

bool edge_split_m2(
  std::vector<std::vector<double>>& V0,
  std::vector<std::vector<int>>& F0,
  std::vector<std::vector<int>>& FF_vec,
  std::vector<std::vector<int>>& FFi_vec,
  std::vector<Property>& props,
  int f0,
  int e0
);


#endif
