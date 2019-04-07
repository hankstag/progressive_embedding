#ifndef MST_H
#define MST_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <set>

void mst(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    int s,
    std::vector<int>& parent,
    std::set<int>& mst_set,
    Eigen::SparseMatrix<int>& graph
);

#endif