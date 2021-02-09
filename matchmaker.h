#ifndef MATCHMAKER
#define MATCHMAKER

#include "shor.h"
#include "path_tracing.h"
#include "Property.h"
#include <igl/harmonic.h>

// given a polygon P and a source mesh M (their boundaries matches in #v)
// using tutte to flattening the M into P
void match_maker(
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    std::vector<Property>& props,
    Eigen::MatrixXd &uv,
    const Eigen::MatrixXd &c,
    const Eigen::VectorXi &ci,
    const Eigen::VectorXi &R,
    const Eigen::VectorXi &T,
    const Eigen::MatrixXd &P,
    const Eigen::VectorXi &mark
);

void remove_ears(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Property>& props);

#endif