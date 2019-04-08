#include "coloring_mesh.h"
#include <igl/adjacency_list.h>
#include <iostream>
void coloring_mesh(const Eigen::MatrixXi& F, Eigen::VectorXi& C){
    
    int n = F.maxCoeff()+1;
    // meanless to color an empty mesh
    assert(n>0);

    std::vector<std::vector<int>> N; // Adjacency list
    igl::adjacency_list(F,N);
    // initialize color for each vertex
    C.setConstant(n,-1);

    size_t vm = 0;// max valence
    for(auto l: N)
        vm = std::max(l.size(),vm);
    
    // assign color 0 for vertex 0
    C[0] = 0;
    // for every 1 ring, mark the color that has been used
    // in theory, it need (vm+1) colors at most
    Eigen::VectorXi G;
    G.setZero(vm+1);
    int nc = 1; // # of color used
    for(int i=1;i<n;i++){
        for(int k=0;k<N[i].size();k++){
            if(C[N[i][k]]!=-1)
                G(C[N[i][k]]) = 1;
        }
        for(int j=0;j<G.rows();j++){
            if(G(j)==0){
                C[i] = j;
                break;
            }
        }
        G.setZero();
    }
}