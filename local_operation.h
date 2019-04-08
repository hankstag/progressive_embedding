#ifndef LOCAL_OPERATIONS_H
#define LOCAL_OPERATIONS_H

#include <Eigen/Core>
#include <igl/edge_lengths.h>
#include "cut_mesh/edge_flaps.h"
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/unique.h>
#include <igl/intersect.h>
#include <igl/setunion.h>
#include <igl/list_to_matrix.h>
#include <igl/setdiff.h>
#include <queue>
#include <set>
#include <map>
#include <iostream>

#define IGL_COLLAPSE_EDGE_NULL 0

void neighbor(const int e, bool op, 
    const Eigen::MatrixXi& F, 
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE, 
    std::vector<int>& NF
);

std::vector<int> circulation(
    const int e,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE
);

void neighbor(const int e, bool op, 
    const Eigen::MatrixXi& F, 
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE, 
    std::vector<int>& NF
);

void remove_empty_faces(Eigen::MatrixXi& F);

bool edge_collapse_is_valid(
    const int ae,
    const Eigen::RowVectorXd & p,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE,
    std::vector<int>& NF
);

bool edge_flip(
    const Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    int e,
    Eigen::MatrixXi& dEF,
    Eigen::MatrixXi& dEI,
    Eigen::VectorXi& EE,
    Eigen::MatrixXi& allE
);

void collapsing(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& uv,
    double eps
);

bool collapse_edge(
    const int ae,
    const Eigen::RowVectorXd & p,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXi& dEF,
    Eigen::MatrixXi& dEI,
    Eigen::VectorXi& EE,
    Eigen::MatrixXi& allE,
    const std::set<int>& prior={}
);

#endif