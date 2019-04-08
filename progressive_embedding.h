#ifndef PROGRESSIVE_TUTTE
#define PROGRESSIVE_TUTTE

#include <igl/harmonic.h>
bool progressive_fix(
    const Eigen::VectorXi& c,
    Eigen::VectorXi& bi,
    Eigen::MatrixXd& b,
    const Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& uv
);

template <typename DerivedT>
double doublearea(const Eigen::MatrixBase<DerivedT>& T);

#endif