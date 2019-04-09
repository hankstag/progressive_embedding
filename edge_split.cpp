#include "edge_split.h"
#include <vector>

bool edge_split(
    std::vector<std::vector<double>>& V0,
    std::vector<std::vector<int>>& F0,
    std::vector<std::vector<int>>& FF_vec,
    std::vector<std::vector<int>>& FFi_vec,
    int f0,
    int e0
){
    int f1 = FF_vec[f0][e0];
    if(f1 == -1) return false;
    int e1 = int(FFi_vec[f0][e0]);
    int e01 = (e0 + 1) % 3;
    int e02 = (e0 + 2) % 3;
    int e11 = (e1 + 1) % 3;
    int e12 = (e1 + 2) % 3;
    int f01 = int(FF_vec[f0][e01]);
    int f02 = int(FF_vec[f0][e02]);
    int f11 = int(FF_vec[f1][e11]);
    int f12 = int(FF_vec[f1][e12]);

    int u1 = F0[f0][e01];
    int u0 = F0[f1][e11];
    int v0 = F0[f0][e02];
    int v1 = F0[f1][e12];
    // pushback middle
    V0.push_back({(V0[u0][0]+V0[u1][0])/2,
                 (V0[u0][1]+V0[u1][1])/2,
                 (V0[u0][2]+V0[u1][2])/2
                 });
    F0.push_back(F0[f0]);
    F0.push_back(F0[f1]);
    int ux = V0.size()-1;

    F0[f1][e1] = ux;
    F0[f0][e0] = ux;
    int fx1 = F0.size()-1;
    int fx0 = F0.size()-2;
    F0[fx1][e11] = ux;
    F0[fx0][e01] = ux;

    if(f12 != -1){
        FF_vec[f12][FFi_vec[f1][e12]] = fx1;
    }
    if(f02 != -1){
        FF_vec[f02][FFi_vec[f0][e02]] = fx0;
    }

    FF_vec.push_back(FF_vec[f0]);
    FF_vec.push_back(FF_vec[f1]);
    FF_vec[f0][e02] = fx0;
    FF_vec[f0][e0] = fx1;
    FF_vec[f1][e12] = fx1;
    FF_vec[f1][e1] = fx0;
    FF_vec[fx1][e11] = f1;
    FF_vec[fx0][e01] = f0;

    FFi_vec.push_back(FFi_vec[f0]);
    FFi_vec.push_back(FFi_vec[f1]);
    FFi_vec[f0][e02] = e01;
    FFi_vec[f1][e12] = e11;
    FFi_vec[fx0][e01] = e02;
    FFi_vec[fx1][e11] = e12;
    return true;
}

bool edge_split(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXi& FF,
    Eigen::MatrixXi& FFi,
    int f0,
    int e0
){
    int f1 = FF(f0, e0);
    if(f1 == -1) return false;
    int e1 = int(FFi(f0, e0));
    int e01 = (e0 + 1) % 3;
    int e02 = (e0 + 2) % 3;
    int e11 = (e1 + 1) % 3;
    int e12 = (e1 + 2) % 3;
    int f01 = int(FF(f0, e01));
    int f02 = int(FF(f0, e02));
    int f11 = int(FF(f1, e11));
    int f12 = int(FF(f1, e12));

    int u1 = F(f0, e01);
    int u0 = F(f1, e11);
    int v0 = F(f0, e02);
    int v1 = F(f1, e12);

    // pushback middle
    V.conservativeResize(V.rows()+1,V.cols());
    V.bottomRows(1)<<(V.row(u0)+V.row(u1))/2;
    F.conservativeResize(F.rows()+2,3);
    int ux = V.rows()-1;
    F.bottomRows(2) << F.row(f0),F.row(f1);
    F(f1,e1) = ux;
    F(f0,e0) = ux;
    int fx1 = F.rows()-1;
    int fx0 = F.rows()-2;
    F(fx1,e11) = ux;
    F(fx0,e01) = ux;

    if(f12 != -1){
        FF(f12,FFi(f1,e12)) = fx1;
    }
    if(f02 != -1){
        FF(f02,FFi(f0,e02)) = fx0;
    }

    FF.conservativeResize(FF.rows()+2,3);
    FF.bottomRows(2)<<FF.row(f0),FF.row(f1);
    FF(f0,e02) = fx0;
    FF(f0,e0) = fx1;
    FF(f1,e12) = fx1;
    FF(f1,e1) = fx0;
    FF(fx1,e11) = f1;
    FF(fx0,e01) = f0;

    FFi.conservativeResize(FFi.rows()+2,3);
    FFi.bottomRows(2)<<FFi.row(f0),FFi.row(f1);
    FFi(f0,e02) = e01;
    FFi(f1,e12) = e11;
    FFi(fx0,e01) = e02;
    FFi(fx1,e11) = e12;
    return true;
}