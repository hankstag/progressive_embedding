#ifndef PATH_TRACING_H
#define PATH_TRACING_H

#include <igl/adjacency_list.h>
#include <set>

// tracing path in mesh
bool path_tracing(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::pair<int,int>& T,
    const std::set<int>& A, // vertices needs to avoid
    const std::vector<std::pair<int,int>>& Mask,
    const std::set<int>& no_enter_f,
    const Eigen::MatrixXi& TT_x,
    Eigen::MatrixXd& Vn,
    Eigen::MatrixXi& Fn,
    std::vector<int>& E
);

void direct_geodesic(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::vector<std::pair<int,int>>& Mask, // [E*2] matrix, [faceid, index of edge that's a cut]
    const std::set<int>& no_enter_f,
    const Eigen::MatrixXi& TT_x,
    int s,
    int t,
    Eigen::MatrixXd& Vn,
    Eigen::MatrixXi& Fn,
    std::vector<int>& L
);

#endif
