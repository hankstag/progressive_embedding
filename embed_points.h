#ifndef EMBED_POINTS
#define EMBED_POINTS

#include <Eigen/Core>

void embed_points(
    const Eigen::MatrixXd& C,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
);
#endif


// #include <unordered_map>
// #include <igl/boundary_loop.h>
// #include <igl/unproject_onto_mesh.h>

// void sequence(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& F,
//     const Eigen::VectorXi& I,
//     const Eigen::MatrixXi& E,    
//     const std::vector<Eigen::MatrixXd>& P_set,
//     Eigen::VectorXi& Seq,
//     Eigen::VectorXi& mask // indicating whether tracing a path between i and i+1
// );

// void split(
//     Eigen::MatrixXd& V,
//     Eigen::MatrixXi& F,
//     const int N,
//     std::vector<int>& path,
//     std::unordered_map<int,int>& m
// );

// void cut_mesh_to_disk(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& F,
//     const Eigen::VectorXi& Seq,
//     const Eigen::VectorXi& mask,    
//     Eigen::MatrixXd& Vn,
//     Eigen::MatrixXi& Fn,
//     std::unordered_map<int,int>& m,
//     std::vector<std::vector<int>>& cuts
// );

// void embed_constraints(
//     Eigen::MatrixXd& V,
//     Eigen::MatrixXi& F,
//     std::vector<Eigen::MatrixXd>& P_set,
//     const Eigen::MatrixXd& c,
//     const Eigen::VectorXi& ci,
//     const Eigen::MatrixXi& E,
//     std::vector<Eigen::VectorXi>& T,
//     Eigen::VectorXi& T0,
//     Eigen::VectorXi& R,
//     Eigen::MatrixXd& Pn,
//     Eigen::MatrixXi& En
// );

// // void embed_constraints(
// //     Eigen::MatrixXd& V,
// //     Eigen::MatrixXi& F,
// //     const Eigen::MatrixXd& P,
// //     const Eigen::MatrixXd& c,
// //     const Eigen::VectorXi& cp,
// //     const Eigen::VectorXi& ci,
// //     const Eigen::MatrixXi& E,
// //     Eigen::VectorXi& R,
// //     Eigen::VectorXi& T,
// //     Eigen::MatrixXd& Pn
// // );


// #endif