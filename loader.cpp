#include "loader.h"
#include <igl/boundary_loop.h>
#include <igl/readOBJ.h>
#include "shor.h"

void load_model(
    std::string model, 
    Eigen::MatrixXd& V,
    Eigen::MatrixXd& uv,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& P,
    Eigen::VectorXi& R,
    Eigen::VectorXi& T
){
    Eigen::MatrixXd CN;
    Eigen::MatrixXi Fuv,FN;
    igl::readOBJ(model,V,uv,CN,F,Fuv,FN);

    // copy face from uv
    std::vector<std::vector<int>> D1,D2;
    igl::boundary_loop(F  ,D1);    
    igl::boundary_loop(Fuv,D2);
    std::cout<<"model size: #F("<<F.rows()<<"), #V("<<V.rows()<<")"<<std::endl;
    std::cout<<"source #BD("<<D1.size()<<")"<<std::endl;
    std::cout<<"target #BD("<<D2.size()<<")"<<std::endl;
    #define USEUV
    #ifdef USEUV
    Eigen::MatrixXd nV(uv.rows(),3);
    for(int i=0;i<F.rows();i++){
        nV.row(Fuv(i,0)) << V.row(F(i,0));
        nV.row(Fuv(i,1)) << V.row(F(i,1));
        nV.row(Fuv(i,2)) << V.row(F(i,2));
    }
    F = Fuv;
    V = nV;
    #endif
    //fill_in_holes(V,F);
    auto D = D2[0];
    P.resize(D.size(),2);
    T.resize(P.rows());
    set_rotation_index(uv,Fuv,R);
    for(int i=0;i<P.rows();i++){
        P.row(i)<<uv.row(D[i]);
        T(i) = D[i];
    }
}

void load_model_with_seam(
    const std::string model,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& polygon,
    Eigen::VectorXi& bd
){
    Eigen::MatrixXd uv,CN;
    Eigen::MatrixXi Fuv,FN;
    igl::readOBJ(model,V,uv,CN,F,Fuv,FN);
    if(uv.rows()!=0){
        Eigen::MatrixXd nV(uv.rows(),3);
        for(int i=0;i<F.rows();i++){
            nV.row(Fuv(i,0)) << V.row(F(i,0));
            nV.row(Fuv(i,1)) << V.row(F(i,1));
            nV.row(Fuv(i,2)) << V.row(F(i,2));
        }
        F = Fuv;
        V = nV;
    }
    igl::boundary_loop(F,bd);
    if(uv.rows()!=0)
        igl::slice(uv,bd,1,polygon);
}