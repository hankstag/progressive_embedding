#ifndef TARGET_POLYGON
#define TARGET_POLYGON

#include <Eigen/Core>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>
#include <igl/slice.h>
#include <unordered_set>
#include <map>
#include <vector>

const int REC_EDGE_NUM = 100;

int match_vertex_in_range(Eigen::MatrixXd& V, Eigen::MatrixXi& F, const int v, std::pair<int,int>& r);

void sample_polyline(
    const Eigen::MatrixXd& L,
    const Eigen::VectorXd& rt,
    Eigen::MatrixXd& SL
);

void line_ratio(
    const Eigen::MatrixXd& py_line,
    Eigen::VectorXd& rt
);


// connect all outer-bd constraints to 
// bd of boundary and extract polygon
void outer_boundary(
    const Eigen::VectorXi& ci,
    std::map<int,int>& m,
    const Eigen::MatrixXd& bnd,
    const Eigen::VectorXi& bndi,
    Eigen::MatrixXd& V2,
    Eigen::MatrixXi& F2,
    Eigen::MatrixXd& P0,
    Eigen::VectorXi& T
);

void generate_hilbert(
    int order,
    const Eigen::MatrixXd& c,
    const Eigen::VectorXi& ci,
    std::map<int,int>& m,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
);

void background_mesh(
    int N,                     
    const Eigen::MatrixXd& c,
    const Eigen::VectorXi& ci,
    std::map<int,int>& m,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
);

void target_polygon(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& p, 
    const Eigen::VectorXi& pi,
    std::vector<Eigen::MatrixXd>& poly
);

#endif
